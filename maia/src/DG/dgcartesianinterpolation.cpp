// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "dgcartesianinterpolation.h"
#include "UTIL/timer.h"
#include "sbpcartesianinterpolation.h"

// Use the namespace that includes all the static helper functions for
// interpolation
using namespace maia::dg::interpolation;

// Standard namespace
using namespace std;

/**
 * \brief Default constructor only sets default values for member variables.
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
 * \date 2012-11-07
 */
DgInterpolation::DgInterpolation()
  : m_sbpMode(-1), m_sbpOperator(), m_polyDeg(-1), m_polyType(-1), m_intMethod(-1), m_nodes(), m_wInt(), m_wBary() {
  TRACE();
  // empty constructor
}


/**
 * \brief Constructor passes arguments to init().
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
 * \date 2012-11-07
 *
 * \param[in] polyDeg Maximum polynomial degree.
 * \param[in] polyType Polynomial type.
 * \param[in] intMethod Integration method.
 */
DgInterpolation::DgInterpolation(const MInt polyDeg,
                                 const DgPolynomialType polyType,
                                 const MInt noNodes1D,
                                 const DgIntegrationMethod intMethod,
                                 const MBool sbpMode,
                                 const MString sbpOperator)
  : m_sbpMode(false),
    m_sbpOperator(""),
    m_polyDeg(-1),
    m_noNodes(-1),
    m_polyType(-1),
    m_intMethod(-1),
    m_nodes(),
    m_wInt(),
    m_wBary() {
  TRACE();

  init(polyDeg, polyType, noNodes1D, intMethod, sbpMode, sbpOperator);
}


/**
 * \brief Destructor clears all member variables.
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
 * \date 2012-11-07
 */
DgInterpolation::~DgInterpolation() {
  TRACE();
  m_sbpMode = false;
  m_sbpOperator = "";
  m_polyDeg = -1;
  m_polyType = -1;
  m_intMethod = -1;
  m_nodes.clear();
  m_wInt.clear();
  m_wBary.clear();
  m_Dhat.clear();
  m_LFace[0].clear();
  m_LFace[1].clear();
  m_LhatFace[0].clear();
  m_LhatFace[1].clear();
}


/**
 * \brief Sets the member variables and calls the appropriate functions to
 *        calculate the nodes and weights etc.
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
 * \date 2012-11-07
 *
 * \param[in] polyDeg Maximum polynomial degree.
 * \param[in] polyType Polynomial type.
 * \param[in] intMethod Integration method.
 * \param[in] noNodes Number of nodes
 * \param[in] sbpMode SBP mode
 * \param[in] sbpOperator Used SBP operator
 */
void DgInterpolation::init(const MInt polyDeg, const DgPolynomialType polyType, const MInt noNodes,
                           const DgIntegrationMethod intMethod, const MBool sbpMode, const MString sbpOperator) {
  TRACE();

  // Set member variables
  m_polyDeg = polyDeg;
  m_noNodes = noNodes;
  m_polyType = polyType;
  m_intMethod = intMethod;
  m_sbpMode = sbpMode;
  m_sbpOperator = sbpOperator;

  if(!m_sbpMode) {
    ASSERT(polyDeg + 1 == noNodes, "polyDeg+1 not equal to noNodes in DG mode! (" << polyDeg << " " << noNodes << ")");
  }

  m_nodes.resize(noNodes);
  m_wInt.resize(noNodes);
  m_wBary.resize(noNodes);
  m_Dhat.resize(noNodes, noNodes);

  if(m_sbpMode) {
    // DDRP case: Operators are loaded as full matrix and not as set of coefficients
    // from which these are constructed from ...
    if(m_sbpOperator == "go4/DDRP307_N32") {
      readDDRP(&m_nodes[0], &m_wInt[0], &m_Dhat[0]);
    } else {
      // Get Coefficients
      MFloatVector sbpA, sbpP, sbpQ;
      // Check if available in sbpoperator.h otherwise try to read from file
      MBool exists = getSBPOperator(noNodes, sbpOperator, sbpA, sbpP, sbpQ);
      if(!exists) {
        readSbpOperator(noNodes, sbpOperator, sbpA, sbpP, sbpQ);
      }

      // TODO labels:DG,toremove Remove, dg_integration method not needed for SBP
      if(m_polyType == DG_POLY_LEGENDRE && m_intMethod == DG_INTEGRATE_GAUSS) {
        mTerm(1, "LG NOT AVAILABLE FOR SBP MODE");
      } else if(m_polyType == DG_POLY_LEGENDRE && m_intMethod == DG_INTEGRATE_GAUSS_LOBATTO) {
        calcSBPNodes(noNodes, &m_nodes[0]);
        calcSBPWeights(noNodes, sbpP, &m_wInt[0]);
      }

      // Calculate the polynomial derivative matrix D, and normalize it with the
      // integration weights to get Dhat.
      calcDhatSBP(noNodes, sbpA, sbpQ, &m_wInt[0], &m_Dhat[0]);
    }
  } else {
    // Go through all possible combinations of the integration method and the
    // polynomial type, and call the corresponding function to calculate the
    // integration nodes and weights.
    if(m_polyType == DG_POLY_LEGENDRE && m_intMethod == DG_INTEGRATE_GAUSS) {
      calcLegendreGaussNodesAndWeights(m_polyDeg, &m_nodes[0], &m_wInt[0]);
    } else if(m_polyType == DG_POLY_LEGENDRE && m_intMethod == DG_INTEGRATE_GAUSS_LOBATTO) {
      calcLegendreGaussLobattoNodesAndWeights(m_polyDeg, &m_nodes[0], &m_wInt[0]);
    }

    // Calculate and set the barycentric weights for Lagrange interpolation
    calcBarycentricWeights(m_polyDeg, &m_nodes[0], &m_wBary[0]);

    // Calculate the polynomial derivative matrix D, and normalize it with the
    // integration weights to get Dhat.
    calcDhat(m_polyDeg, &m_nodes[0], &m_wInt[0], &m_wBary[0], &m_Dhat[0]);
  }

  // Calculate the values of the Lagrange polynomials at the boundaries of the
  // interval [-1,1]
  m_LFace[0].resize(noNodes);
  m_LFace[1].resize(noNodes);

  if(m_sbpMode) {
    calcLinearInterpolationBase(-F1, noNodes, &m_nodes[0], &m_LFace[0][0]);
    calcLinearInterpolationBase(+F1, noNodes, &m_nodes[0], &m_LFace[1][0]);
  } else {
    calcLagrangeInterpolatingPolynomials(-F1, m_polyDeg, &m_nodes[0], &m_wBary[0], &m_LFace[0][0]);
    calcLagrangeInterpolatingPolynomials(+F1, m_polyDeg, &m_nodes[0], &m_wBary[0], &m_LFace[1][0]);
  }

  // Calculate the values of the Lagrange polynomials at the boundaries of the
  // interval [-1,1], and normalize them with the integration weights.
  m_LhatFace[0].resize(noNodes);
  m_LhatFace[1].resize(noNodes);

  if(m_sbpMode) {
    for(MInt i = 0; i < noNodes; i++) {
      m_LhatFace[0][i] = m_LFace[0][i] / m_wInt[i];
      m_LhatFace[1][i] = m_LFace[1][i] / m_wInt[i];
    }
  } else {
    calcLhat(-F1, m_polyDeg, &m_nodes[0], &m_wInt[0], &m_wBary[0], &m_LhatFace[0][0]);
    calcLhat(+F1, m_polyDeg, &m_nodes[0], &m_wInt[0], &m_wBary[0], &m_LhatFace[1][0]);
  }
}

/**
 * \brief Sets the member variables neccessary for the interpolation between sets of nodes.
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
 * \date 2012-11-07
 *
 * \param[in] polyDeg Maximum polynomial degree.
 * \param[in] polyType Polynomial type.
 * \param[in] intMethod Integration method.
 * \param[in] noNodes Number of nodes
 * \param[in] sbpMode SBP mode
 * \param[in] sbpOperator Used SBP operator
 */
void DgInterpolation::initInterpolation(const MInt polyDeg, const DgPolynomialType polyType, const MInt noNodes,
                                        const DgIntegrationMethod intMethod, const MBool sbpMode,
                                        const MString sbpOperator) {
  TRACE();

  if(sbpMode) {
    // Set member variables
    m_polyDeg = polyDeg;
    m_noNodes = noNodes;
    // m_polyType = polyType;
    // m_intMethod = intMethod;
    // m_sbpMode = sbpMode;
    // m_sbpOperator = sbpOperator;

    m_nodes.resize(noNodes);
    // m_wInt.resize(noNodes);
    // m_wBary.resize(noNodes);
    // m_Dhat.resize(noNodes, noNodes);

    calcSBPNodes(noNodes, &m_nodes[0]);
  } else {
    init(polyDeg, polyType, noNodes, intMethod, sbpMode, sbpOperator);
  }
}
