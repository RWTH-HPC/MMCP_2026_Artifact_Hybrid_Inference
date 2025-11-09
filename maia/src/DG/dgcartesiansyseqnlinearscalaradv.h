// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef DGSYSEQNLINEARSCALARADV_H_
#define DGSYSEQNLINEARSCALARADV_H_

#include <algorithm>
#include <limits>
#include <numeric>
#include "INCLUDE/maiaconstants.h"
#include "INCLUDE/maiatypes.h"
#include "UTIL/tensor.h"
#include "dgcartesiansyseqn.h"

template <MInt nDim>
class DgSysEqnLinearScalarAdv : public DgSysEqn<nDim, DgSysEqnLinearScalarAdv<nDim>> {
  // Declare parent a friend so that CRTP can access this class's private
  // methods/members
  friend class DgSysEqn<nDim, DgSysEqnLinearScalarAdv>;

  // Member functions
 public:
  explicit DgSysEqnLinearScalarAdv(MInt solverId);
  void calcInitialCondition(const MFloat t, const MFloat* x, MFloat* nodeVars, MFloat* u) const;
  void calcFlux(const MFloat* nodeVars, const MFloat* u, const MInt noNodes1D, MFloat* flux) const;
  void calcSource(const MFloat* nodeVars, const MFloat* u, const MInt noNodes1D, const MFloat t, const MFloat* x,
                  MFloat* src) const;
  void calcSpongeSource(const MFloat* nodeVars, const MFloat* u, const MInt noNodes1D, const MFloat* eta,
                        MFloat* src) const;
  MFloat getTimeStep(const MFloat* nodeVars, const MFloat* u, const MInt noNodes1D, const MFloat invJacobian,
                     const MInt sbpMode) const;
  void calcRiemann(const MFloat* nodeVarsL, const MFloat* nodeVarsR, const MFloat* stateL, const MFloat* stateR,
                   const MInt noNodes1D, const MInt dirId, MFloat* flux) const;
  void primToCons(const MFloat* prim, MFloat* cons) const;
  void consToPrim(const MFloat* cons, MFloat* prim) const;

  MFloat cflFactor(const MInt polyDeg) const;

  void getDefaultNodeVars(MFloat* const NotUsed(nodeVars)) const {};
  MBool extendNodeVar(const MInt NotUsed(varId)) const { return true; };

  // Member variables
 private:
  static const MString s_sysEqnName;

  static const MInt s_noVariables = 1;
  static const MString s_consVarNames[s_noVariables];
  static const MString s_primVarNames[s_noVariables];

  static const MInt s_noNodeVars = 0;
  static const MBool s_hasTimeDependentNodeVars = false;
  static constexpr const MString* s_nodeVarNames = nullptr;

  MFloat m_advectionVelocity[nDim]{};

  MString m_dgIntegrationMethod;
  MInt m_dgTimeIntegrationScheme;
  static const MInt s_maxPolyDeg = 31;
  std::array<std::array<MFloat, s_maxPolyDeg + 1>, 2> m_cflFactor;
};

// Initialize static member variables
template <MInt nDim>
const MString DgSysEqnLinearScalarAdv<nDim>::s_sysEqnName = "DG_SYSEQN_LINEARSCALARADV";

template <MInt nDim>
const MString DgSysEqnLinearScalarAdv<nDim>::s_consVarNames[] = {"scalar"};

template <MInt nDim>
const MString DgSysEqnLinearScalarAdv<nDim>::s_primVarNames[] = {"scalar"};


/**
 * \brief Constructor calls parent constructor & loads all necessary properties
 *        for this equation
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
 * \date 2013-04-18
 *
 * \param[in] solverId Current solver id.
 */
template <MInt nDim>
inline DgSysEqnLinearScalarAdv<nDim>::DgSysEqnLinearScalarAdv(MInt solverId)
  : DgSysEqn<nDim, DgSysEqnLinearScalarAdv>(solverId) {
  // Read and set advection velocity
  for(MInt i = 0; i < nDim; i++) {
    m_advectionVelocity[i] = F1;
    /*! \property
    \page propertyPageDG DG
    \section advectionVelocity
    <code>MFloat DgSysEqnLinearScalarAdv::m_advectionVelocity</code> \n
    default = <code>1</code> \n \n
    Sets the advection velocity for the Linear scalar advection equation.
    Keywords: <i>DG, LINEAR_SCALAR_ADVECTION</i>
*/
    m_advectionVelocity[i] =
        Context::getSolverProperty<MFloat>("advectionVelocity", this->m_solverId, AT_, &m_advectionVelocity[i], i);
  }

  // Get the quadrature method being used
  MString dgIntegrationMethod = "DG_INTEGRATE_GAUSS";
  m_dgIntegrationMethod =
      Context::getSolverProperty<MString>("dgIntegrationMethod", this->m_solverId, AT_, &dgIntegrationMethod);

  // Get the time integration scheme being used
  MString dgTimeIntegrationScheme = "DG_TIMEINTEGRATION_CARPENTER_4_5";
  m_dgTimeIntegrationScheme = string2enum(
      Context::getSolverProperty<MString>("dgTimeIntegrationScheme", this->m_solverId, AT_, &dgTimeIntegrationScheme));


  // cfl correction factor table - Obtained by running different test cases
  // to get the maximum stable CFL number solution for polynomial degrees
  // varying from 0 to 15. 1st row -> 2D cases & 2nd row -> 3D cases.
  // Note: The cfl factors are not calibrated and thus not suitable to be used
  //       with RiemannRoe's Flux.
  switch(m_dgTimeIntegrationScheme) {
    case DG_TIMEINTEGRATION_CARPENTER_4_5: {
      default:
        if(m_dgIntegrationMethod == "DG_INTEGRATE_GAUSS") {
          // Legendre-Gauss nodes
          m_cflFactor = {{{{2.448730, 2.420654, 2.171631, 1.834717, 1.652832, 1.497803, 1.369629, 1.246338, 1.151123,
                            1.071777, 1.015625, 0.959473, 0.904541, 0.856934, 0.814209, 0.783691}},
                          {{2.891846, 2.716064, 2.454112, 2.081513, 1.877617, 1.702831, 1.523950, 1.382303, 1.281892,
                            1.203783, 1.132922, 1.060420, 1.002581, 0.960691, 0.912857, 0.868814}}}};
        } else {
          // Legendre-Gauss-Lobatto nodes
          m_cflFactor = {{{{5.277100, 5.277100, 3.908691, 3.056641, 2.852783, 2.486572, 2.232660, 2.036133, 1.865234,
                            1.748047, 1.624756, 1.542969, 1.455078, 1.365967, 1.307373, 1.243896}},
                          {{5.848389, 5.848389, 3.948975, 3.621826, 3.244629, 2.897949, 2.618408, 2.298584, 2.114258,
                            1.972656, 1.834717, 1.741943, 1.669922, 1.539307, 1.470947, 1.406250}}}};
        }
        break;
    }
    case DG_TIMEINTEGRATION_TOULORGEC_4_8: {
      // Using CFL factors used for CARPENTER 4/5
      if(m_dgIntegrationMethod == "DG_INTEGRATE_GAUSS") {
        // Legendre-Gauss nodes
        m_cflFactor = {{{{2.448730, 2.420654, 2.171631, 1.834717, 1.652832, 1.497803, 1.369629, 1.246338, 1.151123,
                          1.071777, 1.015625, 0.959473, 0.904541, 0.856934, 0.814209, 0.783691}},
                        {{2.891846, 2.716064, 2.454112, 2.081513, 1.877617, 1.702831, 1.523950, 1.382303, 1.281892,
                          1.203783, 1.132922, 1.060420, 1.002581, 0.960691, 0.912857, 0.868814}}}};
      } else {
        // Legendre-Gauss-Lobatto nodes
        m_cflFactor = {{{{5.277100, 5.277100, 3.908691, 3.056641, 2.852783, 2.486572, 2.232660, 2.036133, 1.865234,
                          1.748047, 1.624756, 1.542969, 1.455078, 1.365967, 1.307373, 1.243896}},
                        {{5.848389, 5.848389, 3.948975, 3.621826, 3.244629, 2.897949, 2.618408, 2.298584, 2.114258,
                          1.972656, 1.834717, 1.741943, 1.669922, 1.539307, 1.470947, 1.406250}}}};
      }
      break;
    }
    case DG_TIMEINTEGRATION_NIEGEMANN_4_14: {
      if(m_dgIntegrationMethod == "DG_INTEGRATE_GAUSS") {
        // Legendre-Gauss nodes
        m_cflFactor = {{{{9.878051, 7.782532, 6.381652, 5.689331, 5.249816, 4.883361, 4.852050, 4.584167, 4.302368,
                          4.062316, 3.855896, 3.674987, 3.481323, 3.322448, 3.161255, 3.032532}},
                        {{11.289367, 10.013732, 7.770935, 7.036865, 6.197266, 5.667297, 5.417968, 5.106018, 4.787108,
                          4.493713, 4.266418, 4.072753, 3.840820, 3.692383, 3.507995, 3.385070}}}};
      } else {
        // Legendre-Gauss-Lobatto nodes
        m_cflFactor = {{{{16.500915, 16.500915, 10.074035, 8.384399, 7.637573, 6.993957, 6.730712, 6.372375, 6.059265,
                          5.831970, 5.558288, 5.319396, 5.126892, 4.957580, 4.782471, 4.585327}},
                        {{18.998840, 18.998840, 11.339233, 9.595093, 8.869141, 8.107239, 7.406799, 6.981201, 6.632140,
                          6.465149, 6.123046, 5.956054, 5.783263, 5.485229, 5.329834, 5.112976}}}};
      }
      break;
    }
    case DG_TIMEINTEGRATION_NIEGEMANN_4_13: {
      // Using CFL factors used for CARPENTER 4/5
      if(m_dgIntegrationMethod == "DG_INTEGRATE_GAUSS") {
        // Legendre-Gauss nodes
        m_cflFactor = {{{{2.448730, 2.420654, 2.171631, 1.834717, 1.652832, 1.497803, 1.369629, 1.246338, 1.151123,
                          1.071777, 1.015625, 0.959473, 0.904541, 0.856934, 0.814209, 0.783691}},
                        {{2.891846, 2.716064, 2.454112, 2.081513, 1.877617, 1.702831, 1.523950, 1.382303, 1.281892,
                          1.203783, 1.132922, 1.060420, 1.002581, 0.960691, 0.912857, 0.868814}}}};
      } else {
        // Legendre-Gauss-Lobatto nodes
        m_cflFactor = {{{{5.277100, 5.277100, 3.908691, 3.056641, 2.852783, 2.486572, 2.232660, 2.036133, 1.865234,
                          1.748047, 1.624756, 1.542969, 1.455078, 1.365967, 1.307373, 1.243896}},
                        {{5.848389, 5.848389, 3.948975, 3.621826, 3.244629, 2.897949, 2.618408, 2.298584, 2.114258,
                          1.972656, 1.834717, 1.741943, 1.669922, 1.539307, 1.470947, 1.406250}}}};
      }
      break;
    }
    case DG_TIMEINTEGRATION_TOULORGEC_3_7: {
      // Using CFL factors used for CARPENTER 4/5
      if(m_dgIntegrationMethod == "DG_INTEGRATE_GAUSS") {
        // Legendre-Gauss nodes
        m_cflFactor = {{{{2.448730, 2.420654, 2.171631, 1.834717, 1.652832, 1.497803, 1.369629, 1.246338, 1.151123,
                          1.071777, 1.015625, 0.959473, 0.904541, 0.856934, 0.814209, 0.783691}},
                        {{2.891846, 2.716064, 2.454112, 2.081513, 1.877617, 1.702831, 1.523950, 1.382303, 1.281892,
                          1.203783, 1.132922, 1.060420, 1.002581, 0.960691, 0.912857, 0.868814}}}};
      } else {
        // Legendre-Gauss-Lobatto nodes
        m_cflFactor = {{{{5.277100, 5.277100, 3.908691, 3.056641, 2.852783, 2.486572, 2.232660, 2.036133, 1.865234,
                          1.748047, 1.624756, 1.542969, 1.455078, 1.365967, 1.307373, 1.243896}},
                        {{5.848389, 5.848389, 3.948975, 3.621826, 3.244629, 2.897949, 2.618408, 2.298584, 2.114258,
                          1.972656, 1.834717, 1.741943, 1.669922, 1.539307, 1.470947, 1.406250}}}};
      }
      break;
    }
    case DG_TIMEINTEGRATION_TOULORGEF_4_8: {
      // Using CFL factors used for CARPENTER 4/5
      if(m_dgIntegrationMethod == "DG_INTEGRATE_GAUSS") {
        // Legendre-Gauss nodes
        m_cflFactor = {{{{2.448730, 2.420654, 2.171631, 1.834717, 1.652832, 1.497803, 1.369629, 1.246338, 1.151123,
                          1.071777, 1.015625, 0.959473, 0.904541, 0.856934, 0.814209, 0.783691}},
                        {{2.891846, 2.716064, 2.454112, 2.081513, 1.877617, 1.702831, 1.523950, 1.382303, 1.281892,
                          1.203783, 1.132922, 1.060420, 1.002581, 0.960691, 0.912857, 0.868814}}}};
      } else {
        // Legendre-Gauss-Lobatto nodes
        m_cflFactor = {{{{5.277100, 5.277100, 3.908691, 3.056641, 2.852783, 2.486572, 2.232660, 2.036133, 1.865234,
                          1.748047, 1.624756, 1.542969, 1.455078, 1.365967, 1.307373, 1.243896}},
                        {{5.848389, 5.848389, 3.948975, 3.621826, 3.244629, 2.897949, 2.618408, 2.298584, 2.114258,
                          1.972656, 1.834717, 1.741943, 1.669922, 1.539307, 1.470947, 1.406250}}}};
      }
      break;
    }
  }
}

/**
 * \brief Calculate the initial condition for a certain point in space.
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
 * \date 2012-11-28
 *
 * \param[in] t Time at which the initial condition should be applied.
 * \param[in] x A pointer to the node coordinates.
 * \param[out] u A pointer to the solution storage.
 */
template <MInt nDim>
void DgSysEqnLinearScalarAdv<nDim>::calcInitialCondition(const MFloat t, const MFloat* x, MFloat* NotUsed(nodeVars),
                                                         MFloat* u) const {
  /// \property
  /// \page dgICLs List of initial conditions for solvertype `MAIA_DISCONTINUOUS_GALERKIN` and
  /// `DG_SYSEQN_LINEARSCALARADV` List of all available `initialCondition` in
  /// DgSysEqnLinearScalarAdv<nDim>::calcInitialCondition():
  switch(this->m_initialCondition) {
    case 1: {
      /// \dgSwitchCase{dgICLs, 1, convergence test}
      /// \note The convergence test requires a square/cube domain. If another
      ///       geometry is set, no guarantee can be made
      ///       for the validity of the results.

      // Set domain-specific values and user-defined properties
      const MFloat frequency = this->m_initialNumberWaves; // = no. full oscillations in the domain
      const MFloat A = F1B5;                               // = max. amplitude of the i.c.
      const MFloat length = F2;                            // = side length of square/cube domain
      const MFloat origin[] = {-F1, -F1, -F1};             // lower left corner

      // Calculate initialization value
      const MFloat omega = F2 * PI * frequency / length;
      MFloat compound = sin(omega * (x[0] - origin[0] - m_advectionVelocity[0] * t));
      for(MInt i = 1; i < nDim; i++) {
        compound *= cos(omega * (x[i] - origin[i] - F1B4 * length / frequency - m_advectionVelocity[i] * t));
      }
      const MFloat init = F2 + A * compound;

      // Set to init
      u[0] = init;
      break;
    }
    case 2: {
      /// \dgSwitchCase{dgICLs, 2, constant}
      u[0] = F2;
      break;
    }
    case 4: {
      /// \dgSwitchCase{dgICLs, 4, linear x-direction}
      /// \f$ u = 2 + x_1 - a_1 t \f$ <br>

      u[0] = F2;
      u[0] += x[0] - m_advectionVelocity[0] * t;
      break;
    }
    case 5: {
      /// \dgSwitchCase{dgICLs, 5, Gauss function}
      /// \f$ c = [0, 0, 0] \f$ <br>
      /// \f$ \text{pow} = -25.0 * \sum_{i=1}^{\text{nDim}} { \left( x_i - a_i t \right) }^{2} \f$ <br>
      /// \f$ u = 2^{\text{pow}} \f$ <br>

      const MFloat origin[] = {F0, F0, F0}; // center of the Gauss pulse
      MFloat xNormalized = F0;
      for(MInt i = 0; i < nDim; i++) {
        xNormalized += pow(x[i] - origin[i] - m_advectionVelocity[i] * t, F2);
      }
      u[0] = pow(F2, -25.0 * xNormalized);
      break;
    }
    default:
      mTerm(1, AT_, "The specified initial condition (" + std::to_string(this->m_initialCondition) + ")is not valid!");
  }
}


/**
 * \brief Calculates the physical fluxes in all dimensions for all integrations
 *        points within a cell.
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
 * \date 2013-04-18
 *
 * \param[in] u Pointer to solution variables (size: noNodes1D^nDim*noVars)).
 * \param[in] noNodes1D Number of nodes 1D in the cell.
 * \param[out] flux Calculated flux (size: noNodes1D^nDim*noVars*nDim)).
 */
template <MInt nDim>
void DgSysEqnLinearScalarAdv<nDim>::calcFlux(const MFloat* NotUsed(nodeVars), const MFloat* u, const MInt noNodes1D,
                                             MFloat* flux) const {
  const MInt noNodes1D3 = (nDim == 3) ? noNodes1D : 1;
  const MFloatTensor U(const_cast<MFloat*>(u), noNodes1D, noNodes1D, noNodes1D3);
  MFloatTensor f(flux, noNodes1D, noNodes1D, noNodes1D3, nDim);

  for(MInt i = 0; i < noNodes1D; i++) {
    for(MInt j = 0; j < noNodes1D; j++) {
      for(MInt k = 0; k < noNodes1D3; k++) {
        for(MInt l = 0; l < nDim; l++) {
          f(i, j, k, l) = U(i, j, k) * m_advectionVelocity[l];
        }
      }
    }
  }
}


/**
 * Yield an empirically-derived maximum stable CFL number.
 * The user should be able to specify a CFL number of "1.0" in virtually all cases.
 *
 * \author Rodrigo Miguez (rodrigo) rodrigo.miguez@rwth-aachen.de
 * \date 2018-02-28
 *
 * \param[in] polyDeg Polynomial degree in the cell.
 * \param[out] cfl factor obtained from m_cflFactor.
 *
 */
template <MInt nDim>
MFloat DgSysEqnLinearScalarAdv<nDim>::cflFactor(const MInt polyDeg) const {
  if(polyDeg > s_maxPolyDeg) {
    TERMM(1, "Polynomial degree exceeds maximum supported (" + std::to_string(s_maxPolyDeg) + ")!");
  }
  // Note: The factor "0.95" is used to avoid numerical instabilities when using the full
  //       cfl factor value.
  const MFloat factor = 0.95 * m_cflFactor[nDim - 2][polyDeg];

  return factor;
}


/**
 * \brief Calculates the source terms for all integration points within a cell.
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
 * \date 2013-04-18
 *
 * \param[in] u Pointer to solution variables (size: noNodes1D^nDim*noVars)).
 * \param[in] noNodes1D Number of nodes 1D in the cell.
 * \param[in] t Current time.
 * \param[in] x Coordinates of the integration points of the current cell.
 * \param[out] src Pointer to where source terms are stored
 *                 (size: noNodes1D^nDim*noVars)).
 */
template <MInt nDim>
void DgSysEqnLinearScalarAdv<nDim>::calcSource(const MFloat* NotUsed(nodeVars), const MFloat* NotUsed(u),
                                               const MInt noNodes1D, const MFloat NotUsed(t), const MFloat* NotUsed(x),
                                               MFloat* src) const {
  const MInt noNodes1D3 = (nDim == 3) ? noNodes1D : 1;
  MFloatTensor sourceTerms(src, noNodes1D, noNodes1D, noNodes1D3);

  /// \property
  /// \page dgSourceLs List of source terms for solvertype `MAIA_DISCONTINUOUS_GALERKIN` and `DG_SYSEQN_LINEARSCALARADV`
  /// List of all available `sourceTerm` in DgSysEqnLinearScalarAdv<nDim>::calcSource():
  switch(this->m_sourceTerm) {
    case 0: {
      /// \dgSwitchCase{dgSourceLs, 0, no source term}
      sourceTerms.set(F0);
      break;
    }
    default:
      mTerm(1, AT_, "The specified source term is not valid!");
  }
}

/**
 * \brief Calculates the sponge source terms for all integration points within a cell.
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
 * \date 2013-04-18
 *
 * \param[in] nodeVars Pointer to node variables (size: noNodes1D^nDim*noNodeVars)).
 * \param[in] u Pointer to solution variables (size: noNodes1D^nDim*noVars)
 * \param[in] noNodes1D Number of nodes 1D in the cell.
 * \param[in] t Current time.
 * \param[in] x Coordinates of the integration points of the current cell.
 * \param[out] src Pointer to where source terms are stored
 *                 (size: noNodes1D^nDim*noVars)).
 */
template <MInt nDim>
void DgSysEqnLinearScalarAdv<nDim>::calcSpongeSource(const MFloat* NotUsed(nodeVars),
                                                     const MFloat* NotUsed(u),
                                                     const MInt NotUsed(noNodes1D),
                                                     const MFloat* NotUsed(eta),
                                                     MFloat* NotUsed(src)) const {
  mTerm(1, AT_, "Sponge: Implementation for linear scalar adv. is missing");
}

/**
 * \brief Calculate the time step for an explicit time stepping scheme for a
 *        given element.
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
 * \date 2013-04-18
 *
 * \param[in] u Pointer to element variables.
 * \param[in] noNodes1D Current polynomial degree.
 * \param[in] invJacobian Inverse Jacobian of element.
 * \return Calculated time step for the element.
 */
template <MInt nDim>
MFloat DgSysEqnLinearScalarAdv<nDim>::getTimeStep(const MFloat* NotUsed(nodeVars), const MFloat* const NotUsed(u),
                                                  const MInt noNodes1D, const MFloat invJacobian,
                                                  const MInt sbpMode) const {
  // Calculate maximum propagation speed (lambda)
  MFloat maxLambda = 0;
  for(MInt n = 0; n < nDim; n++) {
    maxLambda += fabs(m_advectionVelocity[n]);
  }

  // Calculate time step for current element
  MFloat dt;
  if(sbpMode) {
    dt = this->cfl() * F2 / (invJacobian * maxLambda * (noNodes1D - 1));
  } else {
    const MInt polyDeg = noNodes1D - 1;
    dt = this->cflScaled(polyDeg) * cflFactor(polyDeg) * F2 / (invJacobian * maxLambda);
  }
  return dt;
}


/**
 * \brief Calculates the numerical flux at a surface given two states (left and
 *        right).
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
 * \date 2013-04-18
 *
 * \param[in] stateL Surface variables on the left side of the surface
 *                   (i.e. in -x/-y/-z-direction).
 * \param[in] stateR Surface variables on the right side of the surface
 *                   (i.e. in +x/+y/+z-direction).
 * \param[in] noNodes1D Number of nodes 1D.
 * \param[in] dirId Direction (0 = x-direction, 1 = y-direction,
 *                  2 = z-direction).
 * \param[out] flux Calculated Riemann flux (size: noNodes1D^(nDim-1)*noVars).
 */
template <MInt nDim>
inline void DgSysEqnLinearScalarAdv<nDim>::calcRiemann(const MFloat* NotUsed(nodeVarsL),
                                                       const MFloat* NotUsed(nodeVarsR),
                                                       const MFloat* stateL,
                                                       const MFloat* stateR,
                                                       const MInt noNodes1D,
                                                       const MInt dirId,
                                                       MFloat* flux) const {
  const MInt noNodes1D3 = (nDim == 3) ? noNodes1D : 1;
  const MFloatTensor uL(const_cast<MFloat*>(stateL), noNodes1D, noNodes1D3);
  const MFloatTensor uR(const_cast<MFloat*>(stateR), noNodes1D, noNodes1D3);

  // Calculate maximum eigenvalue
  const MFloat maxLambda = m_advectionVelocity[dirId];

  // Classical upwind flux...
  MFloatTensor riemann(flux, noNodes1D, noNodes1D3);
  for(MInt i = 0; i < noNodes1D; i++) {
    for(MInt j = 0; j < noNodes1D3; j++) {
      riemann(i, j) = F1B2 * ((maxLambda + fabs(maxLambda)) * uL(i, j) + (maxLambda - fabs(maxLambda)) * uR(i, j));
    }
  }
}


/**
 * \brief Calculates a set of primitive variables from a set of conservative
 *        variables.
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
 * \date 2013-04-18
 *
 * \param[in] cons Conservative variable storage (size: noVars).
 * \param[out] prim Primitive variable storage (size: noVars).
 */
template <MInt nDim>
inline void DgSysEqnLinearScalarAdv<nDim>::consToPrim(const MFloat* cons, MFloat* prim) const {
  // Very easy...
  std::copy(cons, cons + 1, prim);
}


/**
 * \brief Calculates a set of conservative variables from a set of primitive
 *        variables.
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
 * \date 2013-04-18
 *
 * \param[in] prim Primitive variable storage (size: noVars).
 * \param[out] cons Conservative variable storage (size: noVars).
 */
template <MInt nDim>
inline void DgSysEqnLinearScalarAdv<nDim>::primToCons(const MFloat* prim, MFloat* cons) const {
  // Very easy...
  std::copy(prim, prim + 1, cons);
}

#endif // DGSYSEQNLINEARSCALARADV_H_
