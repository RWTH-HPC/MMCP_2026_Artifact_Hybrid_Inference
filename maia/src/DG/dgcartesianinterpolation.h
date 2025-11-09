// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef DGINTERPOLATION_H_
#define DGINTERPOLATION_H_
//#define CALC_POLY_DERIV_SORTING

#include <algorithm>
#include <iomanip>

#include "INCLUDE/maiaconstants.h"
#include "INCLUDE/maiatypes.h"
#include "UTIL/maiamath.h"
#include "UTIL/tensor.h"
#include "enums.h"

/**
 * \brief Class stores precalculated values for interpolation & integration on
 *        the reference interval [-1,1].
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
 * \date 2012-11-05
 *
 * \details The following two sources may be cited: \n\n
 *          Kopriva09: David A. Kopriva,
 *                     Implementing Spectral Methods for Partial Differential
 *                     Equations, Springer, 2009\n
 *          HesthavenWarburton08: J. S. Hesthaven and T. Warburton,
 *                                Nodal Discontinuous Galerkin Methods,
 *                                Springer, 2008
 */
class DgInterpolation {
  // Member methods
 public:
  DgInterpolation();
  DgInterpolation(const MInt polyDeg, const DgPolynomialType polyType, const MInt noNodes,
                  const DgIntegrationMethod intMethod, const MBool sbpMode, const MString sbpOperator);
  ~DgInterpolation();

  void init(const MInt polyDeg, const DgPolynomialType polyType, const MInt noNodes,
            const DgIntegrationMethod intMethod, const MBool sbpMode, const MString sbpOperator);

  void initInterpolation(const MInt polyDeg, const DgPolynomialType polyType, const MInt noNodes,
                         const DgIntegrationMethod intMethod, const MBool sbpMode, const MString sbpOperator);

  // Member variables
 public:
  // SBP Mode
  MBool m_sbpMode;
  // SBP Operator
  MString m_sbpOperator;
  // Polynomial degree.
  MInt m_polyDeg;
  // No nodes 1D
  MInt m_noNodes;
  // Polynomial type (i.e. Legendre).
  MInt m_polyType;
  // Quadrature method (i.e. Gauss or Gauss-Lobatto).
  MInt m_intMethod;
  // Nodes for integration & interpolation.
  MFloatVector m_nodes;
  // Quadrature weights for integration.
  MFloatVector m_wInt;
  // Barycentric weights for Lagrange interpolation.
  MFloatVector m_wBary;
  // Derivative matrix normalized by integration weights.
  MFloatMatrix m_Dhat;
  // Values of Lagrange polynomials at x = -1.0 (index 0) and +1.0 (index 1)
  MFloatVector m_LFace[2];
  // Values of normalized Lagrange polynomials at x = -1.0 (index 0)
  // and +1.0 (index 1)
  MFloatVector m_LhatFace[2];
  // mutable spline m_spline;
};


namespace maia {
namespace dg {

/**
 * \brief Holds helper functions for the interpolation.
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
 * \date 2012-11-05
 *
 * \details The following two sources may be cited: \n\n
 *          Kopriva09: David A. Kopriva,
 *                     Implementing Spectral Methods for Partial Differential
 *                     Equations, Springer, 2009\n
 *          HesthavenWarburton08: J. S. Hesthaven and T. Warburton,
 *                                Nodal Discontinuous Galerkin Methods,
 *                                Springer, 2008
 */
namespace interpolation {


/**
 * \brief Evaluates the Legendre polynomial and its derivative of degree Nmax at
 *        point x.
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
 * \date 2012-11-06
 *
 * \param[in] Nmax The polynomial degree.
 * \param[in] x The evaluation point.
 * \param[out] polynomial The resulting value of the Legendre polynomial.
 * \param[out] derivative The resulting value of the derivative.
 *
 * \details Taken from Kopriva09, p. 63, algorithm 22.
 */
inline void calcLegendrePolyAndDeriv(MInt Nmax, MFloat x, MFloat* polynomial, MFloat* derivative) {
  maia::math::calculateLegendrePolyAndDeriv(Nmax, x, polynomial, derivative);
}

/**
 * \brief Calculate the Gauss integration nodes and weight for the Legendre
 *        polynomials on the interval [-1,1].
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
 * \date 2012-11-07
 *
 * \param[in] Nmax Maximum degree of the polynomials.
 * \param[out] nodes The resulting integration nodes.
 * \param[out] wInt The resulting integration weights.
 *
 * \details Taken from Kopriva09, p. 64, algorithm 23.
 */
inline void calcLegendreGaussNodesAndWeights(MInt Nmax, MFloat* nodes, MFloat* wInt) {
  maia::math::calculateLegendreGaussNodesAndWeights(Nmax, nodes, wInt);
}


/**
 * \brief Auxiliary function (only used by
 *        calcLegendreGaussLobattoNodesAndWeights())
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
 * \date 2012-11-07
 *
 * \param[in] Nmax The polynomial degree.
 * \param[in] x The evaluation point.
 * \param[out] q The resulting q value.
 * \param[out] qDeriv The resulting q' value.
 * \param[out] poly The resulting value of the Legendre polynomial.
 *
 * \details Taken from Kopriva09, p. 65, algorithm 24.
 */
inline void calcQandL(MInt Nmax, MFloat x, MFloat& q, MFloat& qDeriv, MFloat& poly) {
  TRACE();
  // This may only be called for Nmax >= 2 (see Kopriva09)
  if(Nmax < 2) {
    mTerm(1, AT_, "Gauss-Lobatto needs at least two nodes!");
  }

  MFloat polyLast1, polyLast2, derivLast1, derivLast2, deriv;

  polyLast2 = F1;
  polyLast1 = x;
  derivLast2 = F0;
  derivLast1 = F1;

  for(MInt k = 2; k <= Nmax; k++) {
    poly = (F2 * k - F1) / k * x * polyLast1 - (k - F1) / k * polyLast2;
    deriv = derivLast2 + (F2 * k - F1) * polyLast1;
    polyLast2 = polyLast1;
    polyLast1 = poly;
    derivLast2 = derivLast1;
    derivLast1 = deriv;
  }

  MFloat polyNp1, derivNp1;
  MInt k = Nmax + 1;

  polyNp1 = (F2 * k - F1) / k * x * poly - (k - F1) / k * polyLast2;
  derivNp1 = derivLast2 + (F2 * k - F1) * polyLast1;

  q = polyNp1 - polyLast2;
  qDeriv = derivNp1 - derivLast2;
}


/**
 * \brief Calculate the Gauss-Lobatto integration nodes and weight for the
 *        Legendre polynomials on the interval [-1,1].
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
 * \date 2012-11-07
 *
 * \param[in] Nmax Maximum degree of the polynomials.
 * \param[out] nodes The resulting integration nodes.
 * \param[out] wInt The resulting integration weights.
 *
 * \details Taken from Kopriva09, p. 66, algorithm 25.
 */
inline void calcLegendreGaussLobattoNodesAndWeights(MInt Nmax, MFloat* nodes, MFloat* wInt) {
  TRACE();
  // Gauss-Lobatto only valid for degree >= 1
  ASSERT(Nmax > 0, "Number of nodes must be greater than zero!");

  // Reset nodes and weights
  const MInt noNodes = Nmax + 1;
  std::fill_n(&nodes[0], noNodes, F0);
  std::fill_n(&wInt[0], noNodes, F0);

  // Return values for polynomial and derivative
  MFloat q, qDeriv, poly;

  // Set tolerance and number of iterations. According to Kopriva09,
  // 1.0E-15 and 10 should be more than sufficient.
  const MFloat tol = F4 * MFloatEps;
  const MInt noIterations = 10;

  // Catch simple cases before going into the full loop
  if(Nmax == 1) {
    nodes[0] = -F1;
    wInt[0] = F1;
    nodes[1] = F1;
    wInt[1] = wInt[0];
  } else {
    nodes[0] = -F1;
    wInt[0] = F2 / (Nmax * (Nmax + F1));
    nodes[Nmax] = F1;
    wInt[Nmax] = wInt[0];

    for(MInt j = 1; j < (Nmax + 1) / 2; j++) {
      // Calculate starting guess for Newton method
      nodes[j] = -cos((j + F1B4) * PI / Nmax - F3B8 / Nmax / PI * F1 / (j + F1B4));

      // Use Newton method to find root of Legendre polynomial
      // -> this is also the integration node
      for(MInt k = 0; k < noIterations; k++) {
        calcQandL(Nmax, nodes[j], q, qDeriv, poly);
        const MFloat delta = -q / qDeriv;
        nodes[j] += delta;

        // Stop iterations if error is small enough
        if(fabs(delta) <= tol * fabs(nodes[j])) {
          break;
        }
      }

      // Calculate weight
      calcQandL(Nmax, nodes[j], q, qDeriv, poly);
      wInt[j] = F2 / (Nmax * (Nmax + F1) * poly * poly);

      // Set nodes and weights according to symmetry properties
      nodes[Nmax - j] = -nodes[j];
      wInt[Nmax - j] = wInt[j];
    }
  }

  // If odd number of nodes (noNodes = Nmax + 1), set center node toi
  // origin (0.0) and calculate weight
  if(Nmax % 2 == 0) {
    calcQandL(Nmax, F0, q, qDeriv, poly);
    nodes[Nmax / 2] = F0;
    wInt[Nmax / 2] = F2 / (Nmax * (Nmax + F1) * poly * poly);
  }
}

/**
 * \brief Calculates the barycentric weights for Lagrange interpolation at thei
 *        specified nodes.
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
 * \date 2012-11-05
 *
 * \param[in] Nmax The polynomial degree.
 * \param[in] nodes Nodes of the interpolating polynomial.
 * \param[out] weights The barycentric weights.
 *
 * \details Taken from Kopriva09, p. 75, algorithm 30
 */
inline void calcBarycentricWeights(MInt Nmax, const MFloat* nodes, MFloat* weights) {
  TRACE();

  const MInt noNodes = Nmax + 1;
  std::fill_n(&weights[0], noNodes, F1);

  // Main iteration loop to calculate inverse weights
  for(MInt j = 1; j < noNodes; j++) {
    for(MInt k = 0; k < j; k++) {
      weights[k] = weights[k] * (nodes[k] - nodes[j]);
      weights[j] = weights[j] * (nodes[j] - nodes[k]);
    }
  }

  // Invert previous results to get final weights
  for(MInt j = 0; j < noNodes; j++) {
    weights[j] = F1 / weights[j];
  }
}


/**
 * \brief Calculates the interpolated value at point x given a set of nodes,
 *        values, and weights.
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
 * \date 2012-11-07
 *
 * \param[in] x The point which should be evaluated.
 * \param[in] Nmax The maximum polynomial degree.
 * \param[in] nodes The set of Nmax+1 node locations.
 * \param[in] values The set of values at the nodes.
 * \param[in] wBary The barycentric weights for the given set of nodes.
 *
 * \return The interpolated value.
 *
 * \details Taken from Kopriva09, p. 75, algorithm 31.
 */
inline MFloat getLagrangeInterpolation(MFloat x, MInt Nmax, const MFloat* nodes, const MFloat* values,
                                       const MFloat* wBary) {
  // TRACE();

  const MInt noNodes = Nmax + 1;
  MFloat numerator = F0;
  MFloat denominator = F0;

  for(MInt j = 0; j < noNodes; j++) {
    if(approx(x, nodes[j], MFloatEps)) {
      return values[j];
    }
    MFloat t = wBary[j] / (x - nodes[j]);
    numerator += t * values[j];
    denominator += t;
  }

  return numerator / denominator;
}


/**
 * \brief Calculates the values of the Lagrangian polynomials l_j for a given
 *        point x in [-1,1].
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
 * \date 2012-11-05
 *
 * \param[in] x The point at which the polynomials should be evaluated.
 * \param[in] noNodes Number of nodes..
 * \param[in] nodes The nodes of the Lagrange polynomials.
 * \param[in] wBary The barycentric weights of the Lagrange polynomials.
 * \param[out] polynomials The Nmax+1 Lagrange polynomials l_j, evaluated at x.
 *
 * \details Taken from Kopriva09, p. 77, algorithm 34.
 */
inline void calcLagrangeInterpolatingPolynomials(const MFloat x, const MInt polyDeg, const MFloat* nodes,
                                                 const MFloat* wBary, MFloat* polynomials) {
  // TRACE();
  const MInt noNodes = polyDeg + 1;
  std::fill_n(&polynomials[0], noNodes, F0);

  // Catch situation where evaluation point is on node, which saves considerable
  // computational time
  for(MInt j = 0; j < noNodes; j++) {
    if(approx(x, nodes[j], MFloatEps)) {
      polynomials[j] = F1;
      return;
    }
  }

  // Otherwise calculate Lagrange interpolating polynomials according to
  // Eq. (3.39) in Kopriva09
  MFloat s = F0;
  MFloat t = F0;
  for(MInt j = 0; j < noNodes; j++) {
    t = wBary[j] / (x - nodes[j]);
    polynomials[j] = t;
    s += t;
  }
  for(MInt j = 0; j < noNodes; j++) {
    polynomials[j] = polynomials[j] / s;
  }
}

/**
 * \brief Calculates the Lagrange polynomials evaluated at point x in [-1,1]
 *        and pre-divides them by the integration weights.
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
 * \date 2012-11-08
 *
 * \param[in] x Evaluation point.
 * \param[in] polyDeg Maximum polynomial degree.
 * \param[in] nodes The integration nodes.
 * \param[in] wInt The integration weights.
 * \param[in] wBary The barycentric weights.
 * \param[out] Lhat The Lagrange polynomials divided by the integration weights,
 *                  i.e. Lhat_j(x) = L_j(x) / w(j).
 */
inline void calcLhat(const MFloat x, const MInt polyDeg, const MFloat* nodes, const MFloat* wInt, const MFloat* wBary,
                     MFloat* Lhat) {
  TRACE();
  const MInt noNodes = polyDeg + 1;

  // Calculate the Lagrange polynomials and evaluate them at x
  calcLagrangeInterpolatingPolynomials(x, polyDeg, &nodes[0], &wBary[0], Lhat);
  // Normalize by integration weights
  for(MInt j = 0; j < noNodes; j++) {
    Lhat[j] = Lhat[j] / wInt[j];
  }
}

#ifdef CALC_POLY_DERIV_SORTING
inline MBool ascendingAbsVal(MFloat i, MFloat j) { return (std::fabs(i) < std::fabs(j)); }
#endif


/**
 * \brief Calculates the first derivative approximation matrix.
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
 * \date 2012-11-05
 *
 * \param[in] Nmax Maximum polynomial degree.
 * \param[in] nodes The nodes of the Lagrange polynomial.
 * \param[in] wBary The barycentric weights of the Lagrange polynomial.
 * \param[out] derivMatrix The derivative matrix.
 *
 * \details Taken from Kopriva09, p. 82, algorithm 37.
 */
inline void calcPolynomialDerivativeMatrix(MInt Nmax, const MFloat* nodes, const MFloat* wBary, MFloat* derivMatrix) {
  TRACE();

#ifdef CALC_POLY_DERIV_SORTING
  const MInt noNodes = Nmax + 1;
  MFloatMatrix D(&derivMatrix[0], noNodes, noNodes);
  std::vector<MFloat> column(noNodes);
  D.set(F0);

  // Iterate over all indices and calculate the D_ij according to Eq. (3.48).
  // For i == j, use the negative sum trick.
  for(MInt i = 0; i < noNodes; i++) {
    for(MInt j = 0; j < noNodes; j++) {
      if(j != i) {
        D(i, j) = wBary[j] / wBary[i] * F1 / (nodes[i] - nodes[j]);
        column[j] = D(i, j);
      } else {
        column[j] = F0;
      }
    }
    std::sort(column.begin(), column.end(), ascendingAbsVal);
    for(MInt j = 0; j < noNodes; j++) {
      D(i, i) -= column[j];
    }
  }


#else
  const MInt noNodes = Nmax + 1;
  MFloatMatrix D(&derivMatrix[0], noNodes, noNodes);
  D.set(F0);

  // Iterate over all indices and calculate the D_ij according to Eq. (3.48).
  // For i == j, use the negative sum trick.
  for(MInt i = 0; i < noNodes; i++) {
    for(MInt j = 0; j < noNodes; j++) {
      if(j != i) {
        D(i, j) = wBary[j] / wBary[i] * F1 / (nodes[i] - nodes[j]);
        D(i, i) -= D(i, j);
      }
    }
  }
#endif
}

/**
 * \brief Calculates the polynomial derivative matrix and normalizes it using
 *        the integration weights.
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
 * \date 2012-11-08
 *
 * \param[in] Nmax Maximum polynomial degree.
 * \param[in] nodes The integration nodes.
 * \param[in] wInt The integration weights.
 * \param[in] wBary The barycentric weights.
 *
 * \return The normalized derivative matrix.
 *
 * \details This method returns Dhat(j,n) = -D(n,j) * w(n) / w(j) as defined in
 *          Kopriva 09, p. 137, Eq. (4.139).
 */
inline void calcDhat(MInt Nmax, const MFloat* nodes, const MFloat* wInt, const MFloat* wBary, MFloat* dhatMatrix) {
  TRACE();

  // Get number of nodes and standard derivative matrix
  const MInt noNodes = Nmax + 1;
  MFloatMatrix Dhat(&dhatMatrix[0], noNodes, noNodes);
  calcPolynomialDerivativeMatrix(Nmax, &nodes[0], &wBary[0], &Dhat[0]);

  Dhat.transpose();

  // Add normalization using integration weights
  for(MInt j = 0; j < noNodes; j++) {
    for(MInt n = 0; n < noNodes; n++) {
      Dhat(j, n) = -Dhat(j, n) * wInt[n] / wInt[j];
    }
  }
}

/**
 * \brief Extrapolates ("prolongs") the Gauss node values of an element to a
 *        given face.
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
 * \date 2012-12-07
 *
 * \tparam nDim Number of spatial dimensions of the element.
 * \tparam noVariables Number of variables per node.
 * \tparam T,U,V Any container type that overloads operator[] for element access
 *               (including pointers).
 * \param[in] source Pointer to memory where the values at the nodes are
 *                   located.
 * \param[in] faceId Face to which the values are prolonged
 *                   ([-x, x, -y, y] = [0, 1, 2, 3])
 * \param[in] polyDeg Polynomial degree inside the element.
 * \param[in] LFace Array with Lagrange polynomials at x = -1.0 (index 0) and
 *                  +1.0 (index 1)
 * \param[out] destination Pointer to memory where the extrapolated values
 *                         should be stored.
 */
template <MInt nDim, MInt noVariables>
void prolongToFaceGauss(const MFloat* const source,
                        const MInt faceId,
                        const MInt noNodes1D,
                        const MFloat* const LFaceN,
                        const MFloat* const LFaceP,
                        MFloat* const destination) {
  // TRACE();

  const MInt noNodes1D3 = (nDim == 3) ? noNodes1D : 1;
  const MFloatTensor src(const_cast<MFloat*>(source), noNodes1D, noNodes1D, noNodes1D3, noVariables);
  MFloatTensor dest(destination, noNodes1D, noNodes1D3, noVariables);

  // Reset destination to zero
  std::fill_n(destination, noNodes1D * noNodes1D3 * noVariables, F0);

  // Index for x-direction: i
  // Index for y-direction: j
  // Index for z-direction: k
  // Index for variables:   l
  switch(faceId) {
    case 0: // prolong to negative x-direction
      for(MInt i = 0; i < noNodes1D; i++) {
        for(MInt j = 0; j < noNodes1D; j++) {
          for(MInt k = 0; k < noNodes1D3; k++) {
            for(MInt l = 0; l < noVariables; l++) {
              dest(j, k, l) += src(i, j, k, l) * LFaceN[i];
            }
          }
        }
      }
      break;
    case 1: // prolong to positive x-direction
      for(MInt i = 0; i < noNodes1D; i++) {
        for(MInt j = 0; j < noNodes1D; j++) {
          for(MInt k = 0; k < noNodes1D3; k++) {
            for(MInt l = 0; l < noVariables; l++) {
              dest(j, k, l) += src(i, j, k, l) * LFaceP[i];
            }
          }
        }
      }
      break;
    case 2: // prolong to negative y-direction
      for(MInt i = 0; i < noNodes1D; i++) {
        for(MInt j = 0; j < noNodes1D; j++) {
          for(MInt k = 0; k < noNodes1D3; k++) {
            for(MInt l = 0; l < noVariables; l++) {
              dest(i, k, l) += src(i, j, k, l) * LFaceN[j];
            }
          }
        }
      }
      break;
    case 3: // prolong to positive y-direction
      for(MInt i = 0; i < noNodes1D; i++) {
        for(MInt j = 0; j < noNodes1D; j++) {
          for(MInt k = 0; k < noNodes1D3; k++) {
            for(MInt l = 0; l < noVariables; l++) {
              dest(i, k, l) += src(i, j, k, l) * LFaceP[j];
            }
          }
        }
      }
      break;
    case 4: // prolong to negative z-direction
      for(MInt i = 0; i < noNodes1D; i++) {
        for(MInt j = 0; j < noNodes1D; j++) {
          for(MInt k = 0; k < noNodes1D; k++) {
            for(MInt l = 0; l < noVariables; l++) {
              dest(i, j, l) += src(i, j, k, l) * LFaceN[k];
            }
          }
        }
      }
      break;
    case 5: // prolong to positive z-direction
      for(MInt i = 0; i < noNodes1D; i++) {
        for(MInt j = 0; j < noNodes1D; j++) {
          for(MInt k = 0; k < noNodes1D; k++) {
            for(MInt l = 0; l < noVariables; l++) {
              dest(i, j, l) += src(i, j, k, l) * LFaceP[k];
            }
          }
        }
      }
      break;
    default:
      mTerm(1, AT_, "Bad face id.");
  }
}


/**
 * \brief Extrapolates ("prolongs") the Gauss-Lobatto node values of an element
 *        to a given face.
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
 * \date 2012-12-07
 *
 * \tparam nDim Number of spatial dimensions of the element.
 * \tparam noVariables Number of variables per node.
 * \tparam T,U,V,W Any container type that overloads operator[] for element
 *                 access (including pointers).
 * \param[in] source Pointer to memory where the values at the nodes are
 *                   located.
 * \param[in] faceId Face to which the values are prolonged
 *                   ([-x, x, -y, y] = [0, 1, 2, 3])
 * \param[in] polyDeg Polynomial degree inside the element.
 * \param[out] destination Pointer to memory where the extrapolated values
 *                         should be stored.
 */
template <MInt nDim, MInt noVariables, class T, class U>
void prolongToFaceGaussLobatto(const T source, const MInt faceId, const MInt noNodes1D, U destination) {
  // TRACE();

  const MInt noNodes1D3 = (nDim == 3) ? noNodes1D : 1;
  const MFloatTensor src(&source[0], noNodes1D, noNodes1D, noNodes1D3, noVariables);
  MFloatTensor dest(&destination[0], noNodes1D, noNodes1D3, noVariables);

  // Index for x-direction: i
  // Index for y-direction: j
  // Index for z-direction: k
  // Index for variables:   l
  switch(faceId) {
    case 0: // prolong to negative x-direction
      for(MInt j = 0; j < noNodes1D; j++) {
        for(MInt k = 0; k < noNodes1D3; k++) {
          for(MInt l = 0; l < noVariables; l++) {
            dest(j, k, l) = src(0, j, k, l);
          }
        }
      }
      break;
    case 1: // prolong to positive x-direction
      for(MInt j = 0; j < noNodes1D; j++) {
        for(MInt k = 0; k < noNodes1D3; k++) {
          for(MInt l = 0; l < noVariables; l++) {
            dest(j, k, l) = src(noNodes1D - 1, j, k, l);
          }
        }
      }
      break;
    case 2: // prolong to negative y-direction
      for(MInt i = 0; i < noNodes1D; i++) {
        for(MInt k = 0; k < noNodes1D3; k++) {
          for(MInt l = 0; l < noVariables; l++) {
            dest(i, k, l) = src(i, 0, k, l);
          }
        }
      }
      break;
    case 3: // prolong to positive y-direction
      for(MInt i = 0; i < noNodes1D; i++) {
        for(MInt k = 0; k < noNodes1D3; k++) {
          for(MInt l = 0; l < noVariables; l++) {
            dest(i, k, l) = src(i, noNodes1D - 1, k, l);
          }
        }
      }
      break;
    case 4: // prolong to negative z-direction
      for(MInt i = 0; i < noNodes1D; i++) {
        for(MInt j = 0; j < noNodes1D; j++) {
          for(MInt l = 0; l < noVariables; l++) {
            dest(i, j, l) = src(i, j, 0, l);
          }
        }
      }
      break;
    case 5: // prolong to positive z-direction
      for(MInt i = 0; i < noNodes1D; i++) {
        for(MInt j = 0; j < noNodes1D; j++) {
          for(MInt l = 0; l < noVariables; l++) {
            dest(i, j, l) = src(i, j, noNodes1D - 1, l);
          }
        }
      }
      break;
    default:
      mTerm(1, AT_, "Bad face id.");
  }
}


/**
 * \brief Calculate the polynomial interpolation matrix (Vandermonde) to
 *        interpolate from one set of nodes to another.
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
 * \date 2013-05-14
 *
 * \tparam T,U,V,W Any container type that overloads operator[] for element
 *                 access (including pointers).
 * \param[in] noNodesIn
 * \param[in] nodesIn
 * \param[in] noNodesOut
 * \param[in] nodesOut
 * \param[in] wBary
 * \param[out] vandermonde
 */
template <class T, class U, class V, class W>
void calcPolynomialInterpolationMatrix(MInt noNodesIn, const T nodesIn, MInt noNodesOut, const U nodesOut,
                                       const V wBary, W vandermonde) {
  TRACE();

  MFloatMatrix vdm(&vandermonde[0], noNodesOut, noNodesIn);

  for(MInt k = 0; k < noNodesOut; k++) {
    MBool rowHasMatch = false;
    for(MInt j = 0; j < noNodesIn; j++) {
      vdm(k, j) = F0;
      if(approx(nodesOut[k], nodesIn[j], MFloatEps)) {
        rowHasMatch = true;
        vdm(k, j) = F1;
      }
    }

    if(rowHasMatch == false) {
      MFloat s = 0;
      for(MInt j = 0; j < noNodesIn; j++) {
        MFloat t = wBary[j] / (nodesOut[k] - nodesIn[j]);
        vdm(k, j) = t;
        s += t;
      }
      for(MInt j = 0; j < noNodesIn; j++) {
        vdm(k, j) = vdm(k, j) / s;
      }
    }
  }
}


/**
 * \brief
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
 * \date 2013-05-14
 *
 * \tparam nDim Number of spatial dimensions of the element.
 * \tparam T,U,V Any container type that overloads operator[] for element access
 *               (including pointers).
 * \param[in] source
 * \param[in] vandermonde
 * \param[in] noNodesIn
 * \param[in] noNodesOut
 * \param[in] noVariables
 * \param[out] destination
 */
template <MInt nDim, class T, class U, class V>
void interpolateNodes(const T source, U vandermonde, MInt noNodesIn, MInt noNodesOut, MInt noVariables, V destination) {
  // TRACE();
  const MInt noNodesIn1D = noNodesIn;
  const MInt noNodesIn1D2 = (nDim == 2 || nDim == 3) ? noNodesIn1D : 1;
  const MInt noNodesIn1D3 = (nDim == 3) ? noNodesIn1D : 1;
  const MInt noNodesOut1D = noNodesOut;
  const MInt noNodesOut1D2 = (nDim == 2 || nDim == 3) ? noNodesOut1D : 1;
  const MInt noNodesOut1D3 = (nDim == 3) ? noNodesOut1D : 1;
  const MFloatTensor src(&source[0], noNodesIn, noNodesIn1D2, noNodesIn1D3, noVariables);
  MFloatTensor dest(&destination[0], noNodesOut, noNodesOut1D2, noNodesOut1D3, noVariables);
  const MFloatMatrix vdm(&vandermonde[0], noNodesOut, noNodesIn);
  dest.set(F0);

  // Allocate buffers (if necessary) for intermediate results
  MFloat* buf1Ptr = nullptr;
  MFloat* buf2Ptr = nullptr;
  MFloatTensor buffer1, buffer2;
  switch(nDim) {
    case 1: {
      buf1Ptr = &destination[0];
      buf2Ptr = nullptr;
    } break;

    case 2: {
      buffer1.resize(noNodesOut, noNodesIn1D2, noNodesIn1D3, noVariables);
      buf1Ptr = &buffer1[0];
      buf2Ptr = &destination[0];
    } break;

    case 3: {
      buffer1.resize(noNodesOut, noNodesIn1D2, noNodesIn1D3, noVariables);
      buffer2.resize(noNodesOut, noNodesOut, noNodesIn, noVariables);
      buf1Ptr = &buffer1[0];
      buf2Ptr = &buffer2[0];
    } break;

    default:
      mTerm(1, AT_, "Bad dimension - must be 1, 2 or 3.");
  }
  MFloatTensor buf1(buf1Ptr, noNodesOut, noNodesIn1D2, noNodesIn1D3, noVariables);
  buf1.set(F0);

  // x-direction
  for(MInt i = 0; i < noNodesOut1D; i++) {
    for(MInt j = 0; j < noNodesIn1D2; j++) {
      for(MInt k = 0; k < noNodesIn1D3; k++) {
        for(MInt ii = 0; ii < noNodesIn1D; ii++) {
          for(MInt l = 0; l < noVariables; l++) {
            buf1(i, j, k, l) += vdm(i, ii) * src(ii, j, k, l);
          }
        }
      }
    }
  }

  // Return here for 1D since buf2Ptr is a nullptr
  IF_CONSTEXPR(nDim == 1) { return; }

  MFloatTensor buf2(buf2Ptr, noNodesOut, noNodesOut1D2, noNodesIn1D3, noVariables);
  buf2.set(F0);

  // Interpolate in y-direction only for 2D & 3D
  IF_CONSTEXPR(nDim == 2 || nDim == 3) {
    // y-direction
    for(MInt i = 0; i < noNodesOut1D; i++) {
      for(MInt j = 0; j < noNodesOut1D; j++) {
        for(MInt k = 0; k < noNodesIn1D3; k++) {
          for(MInt jj = 0; jj < noNodesIn1D; jj++) {
            for(MInt l = 0; l < noVariables; l++) {
              buf2(i, j, k, l) += vdm(j, jj) * buf1(i, jj, k, l);
            }
          }
        }
      }
    }
  }

  // Interpolate in z-direction only for 3D
  IF_CONSTEXPR(nDim == 3) {
    // z-direction
    for(MInt i = 0; i < noNodesOut1D; i++) {
      for(MInt j = 0; j < noNodesOut1D; j++) {
        for(MInt k = 0; k < noNodesOut1D; k++) {
          for(MInt kk = 0; kk < noNodesIn1D; kk++) {
            for(MInt l = 0; l < noVariables; l++) {
              dest(i, j, k, l) += vdm(k, kk) * buf2(i, j, kk, l);
            }
          }
        }
      }
    }
  }
}

} // namespace interpolation
} // namespace dg
} // namespace maia

#endif /* DGINTERPOLATION_H_ */
