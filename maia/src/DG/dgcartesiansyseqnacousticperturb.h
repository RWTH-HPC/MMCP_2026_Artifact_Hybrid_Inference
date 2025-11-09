// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef DGSYSEQNACOUSTICPERTURB_H_
#define DGSYSEQNACOUSTICPERTURB_H_

#include <algorithm>
#include <array>
#include <iterator>
#include <limits>
#include <numeric>
#include "INCLUDE/maiaconstants.h"
#include "INCLUDE/maiatypes.h"
#include "UTIL/tensor.h"
#include "dgcartesiansyseqn.h"
#include "filter.h"

// Only needed (really: possible to use) in single-threaded applications
#ifndef _OPENMP
#include "MEMORY/scratch.h"
#endif

template <MInt nDim>
class DgSysEqnAcousticPerturb : public DgSysEqn<nDim, DgSysEqnAcousticPerturb<nDim>> {
  // Declare parent a friend so that CRTP can access private methods/members
  friend class DgSysEqn<nDim, DgSysEqnAcousticPerturb>;

  // Member functions
 public:
  explicit DgSysEqnAcousticPerturb(MInt solverId);
  void calcInitialCondition(const MFloat t, const MFloat* x, MFloat* const nodeVars, MFloat* u) const;
  void calcFlux(const MFloat* const nodeVars, const MFloat* const q, const MInt noNodes1D, MFloat* const flux) const;
  void calcSource(const MFloat* const nodeVars, const MFloat* const u, const MInt noNodes1D, const MFloat t,
                  const MFloat* const x, MFloat* const src) const;
  void calcSpongeSource(const MFloat* nodeVars, const MFloat* u, const MInt noNodes1D, const MFloat* eta,
                        MFloat* src) const;
  MFloat getTimeStep(const MFloat* nodeVars, const MFloat* const u, const MInt noNodes1D, const MFloat invJacobian,
                     const MInt sbpMode) const;
  void calcRiemann(const MFloat* nodeVarsL, const MFloat* nodeVarsR, const MFloat* stateL, const MFloat* stateR,
                   const MInt noNodes1D, const MInt dirId, MFloat* flux) const;
  void calcRiemannLaxFriedich(const MFloat* nodeVarsL, const MFloat* nodeVarsR, const MFloat* stateL,
                              const MFloat* stateR, const MInt noNodes1D, const MInt dirId, MFloat* flux) const;
  void calcRiemannRoe(const MFloat* nodeVarsL, const MFloat* nodeVarsR, const MFloat* stateL, const MFloat* stateR,
                      const MInt noNodes1D, const MInt dirId, MFloat* flux) const;

  void primToCons(const MFloat* prim, MFloat* cons) const;
  void consToPrim(const MFloat* cons, MFloat* prim) const;

  void getDefaultNodeVars(MFloat* const nodeVars) const;
  void getDefaultNodeVarsBody(MFloat* const nodeVars) const;
  MBool extendNodeVar(const MInt varId) const;

 private:
  void calcFlux1D(const MFloat* const nodeVars, const MFloat* const q, const MInt noNodes1D, const MInt dirId,
                  MFloat* const flux) const;

  MFloat cflFactor(const MInt polyDeg) const;

  // Member variables
 private:
  static const MString s_sysEqnName;

  static const MInt s_noVariables = nDim + 1;
  static const MString s_consVarNames[s_noVariables];
  static const MString s_primVarNames[s_noVariables];

  // Node vars: U0, V0, W0, RHO0, C0, DC0_DX, DC0_DY, DC0_DZ
  static const MInt s_noNodeVars = 2 * nDim + 2;
  static const MBool s_hasTimeDependentNodeVars = false;
  static const MString s_nodeVarNames[s_noNodeVars];

  std::array<MFloat, nDim> m_meanVelocity;
  MFloat m_meanDensity;
  MFloat m_meanSpeedOfSound;
  MBool m_constantSpeedOfSound = false;
  MFloat m_spongePressureInfy;
  MFloat m_spongeSigma;
  MBool m_compressibleSourceTerm = true;
  MString m_dgIntegrationMethod;
  MInt m_dgTimeIntegrationScheme;
  static const MInt s_maxPolyDeg = 31;
  std::array<std::array<MFloat, s_maxPolyDeg + 1>, 2> m_cflFactor;

 public:
  // Hold indices for primitive and conservative variables
  struct CV;
};

// Initialize static member variables
template <MInt nDim>
const MString DgSysEqnAcousticPerturb<nDim>::s_sysEqnName = "DG_SYSEQN_ACOUSTICPERTURB";

// For variable names, see dgsyseqnacousticperturb.cpp

// Use these classes to access position of member variables
template <>
struct DgSysEqnAcousticPerturb<2>::CV {
  static constexpr const MInt nDim = 2;

  // Conservative variables
  static constexpr const MInt U = 0;
  static constexpr const MInt V = 1;
  static constexpr const MInt UU[nDim] = {0, 1};
  static constexpr const MInt P = nDim;

  // Node variables
  static constexpr const MInt U0 = 0;
  static constexpr const MInt V0 = 1;
  static constexpr const MInt UU0[nDim] = {0, 1};
  static constexpr const MInt RHO0 = 2;
  static constexpr const MInt C0 = 3;
  static constexpr const MInt DC0_DX = 4;
  static constexpr const MInt DC0_DY = 5;
  static constexpr const MInt DC0[nDim] = {4, 5};
};

template <>
struct DgSysEqnAcousticPerturb<3>::CV {
  static constexpr const MInt nDim = 3;

  // Conservative variables
  static constexpr const MInt U = 0;
  static constexpr const MInt V = 1;
  static constexpr const MInt W = 2;
  static constexpr const MInt UU[nDim] = {0, 1, 2};
  static constexpr const MInt P = nDim;

  // Node variables
  static constexpr const MInt U0 = 0;
  static constexpr const MInt V0 = 1;
  static constexpr const MInt W0 = 2;
  static constexpr const MInt UU0[nDim] = {0, 1, 2};
  static constexpr const MInt RHO0 = 3;
  static constexpr const MInt C0 = 4;
  static constexpr const MInt DC0_DX = 5;
  static constexpr const MInt DC0_DY = 6;
  static constexpr const MInt DC0_DZ = 7;
  static constexpr const MInt DC0[nDim] = {5, 6, 7};
};


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
DgSysEqnAcousticPerturb<nDim>::DgSysEqnAcousticPerturb(MInt solverId)
  : DgSysEqn<nDim, DgSysEqnAcousticPerturb>(solverId) {
  // Read and set mean velocity
  for(MInt i = 0; i < nDim; i++) {
    m_meanVelocity[i] = 0.0;
    m_meanVelocity[i] =
        Context::getSolverProperty<MFloat>("meanVelocity", this->m_solverId, AT_, &m_meanVelocity[i], i);
  }

  /*! \property
    \page propertyPageDG DG
    \section meanDensity
    <code>MFloat DgSysEqnAcousticPerturb::m_meanDensity</code>\n
    default = <code>1.0</code>\n \n
    Specify mean density in the acoustic perturbation equations.\n
    Possible values are:
    <ul>
      <li>any float > 0</li>
    </ul>
    Keywords: <i>DISCONTINUOUS_GALERKIN, ACOUSTIC_PERTURBATION_EQUATIONS</i>
  */
  m_meanDensity = 1.0;
  m_meanDensity = Context::getSolverProperty<MFloat>("meanDensity", this->m_solverId, AT_, &m_meanDensity);

  /*! \property
    \page propertyPageDG DG
    \section meanSpeedOfSound
    <code>MFloat DgSysEqnAcousticPerturb::m_meanSpeedOfSound</code>\n
    default = <code>1.0</code>\n \n
    Specify mean speed of sound in the acoustic perturbation equations.\n
    Note: as only a constant value for c0 can be specified, its derivatives
    will always be zero thus setting a nonzero value for the derivatives does
    not make sense.\n
    Possible values are:
    <ul>
      <li>any float > 0</li>
    </ul>
    Keywords: <i>DISCONTINUOUS_GALERKIN, ACOUSTIC_PERTURBATION_EQUATIONS</i>
  */
  m_meanSpeedOfSound = 1.0;
  m_meanSpeedOfSound =
      Context::getSolverProperty<MFloat>("meanSpeedOfSound", this->m_solverId, AT_, &m_meanSpeedOfSound);

  /*! \property
    \page propertyPageDG DG
    \section constantSpeedOfSound
    <code>MBool DgSysEqnAcousticPerturb::m_constantSpeedOfSound</code>\n
    default = <code>false</code>\n \n
    Overwrite the speed of sound loaded from a mean file in the coupling class with the specified
    constant value in the whole domain. Currently only relevant for IC 790.\n
    Keywords: <i>DISCONTINUOUS_GALERKIN, ACOUSTIC_PERTURBATION_EQUATIONS</i>
  */
  m_constantSpeedOfSound = false;
  m_constantSpeedOfSound =
      Context::getSolverProperty<MBool>("constantSpeedOfSound", this->m_solverId, AT_, &m_constantSpeedOfSound);

  m_log << "APE default mean variables: ";
  for(MInt i = 0; i < nDim; i++) {
    m_log << "u_" << i << " = " << m_meanVelocity[i] << ", ";
  }
  m_log << "rho0 = " << m_meanDensity << ", c0 = " << m_meanSpeedOfSound
        << " (c0 = constant: " << m_constantSpeedOfSound << ")" << std::endl;

  /*! \property
    \page propertyPageDG DG
    \section spongePressureInfy
    <code>MFloat DgSysEqnAcousticPerturb::m_spongePressureInfy</code>\n
    default = <code>0.0</code>\n \n
    Desired pressure value on infinity (effectively on the boundary) for sponge boundaries.\n
    Possible values are:
    <ul>
      <li>any float</li>
    </ul>
    Keywords: <i>DISCONTINUOUS_GALERKIN, ACOUSTIC_PERTURBATION_EQUATIONS, SPONGE</i>
  */
  m_spongePressureInfy = 0.0;
  m_spongePressureInfy =
      Context::getSolverProperty<MFloat>("spongePressureInfy", this->m_solverId, AT_, &m_spongePressureInfy);

  /*! \property
    \page propertyPageDG DG
    \section spongeSigma
    <code>MFloat DgSysEqnAcousticPerturb::m_spongeSigma</code>\n
    default = <code>1.0</code>\n \n
    Heuristic value for sponge factor - has to be estimated, usually 0.5 - 1.\n
    Possible values are:
    <ul>
      <li>any float > 0</li>
    </ul>
    Keywords: <i>DISCONTINUOUS_GALERKIN, ACOUSTIC_PERTURBATION_EQUATIONS, SPONGE</i>
  */
  m_spongeSigma = 1.0;
  m_spongeSigma = Context::getSolverProperty<MFloat>("spongeSigma", this->m_solverId, AT_, &m_spongeSigma);

  // Enables compressible source term (due to recasting the APE in conservative form)
  m_compressibleSourceTerm = true;
  m_compressibleSourceTerm =
      Context::getSolverProperty<MBool>("compressibleSourceTerm", this->m_solverId, AT_, &m_compressibleSourceTerm);

  const MString statusMsg = (m_compressibleSourceTerm) ? "enabled" : "disabled";
  m_log << "APE: compressible source term " << statusMsg << std::endl;

  // Get the quadrature method being used
  MString dgIntegrationMethod = "DG_INTEGRATE_GAUSS";
  m_dgIntegrationMethod =
      Context::getSolverProperty<MString>("dgIntegrationMethod", this->m_solverId, AT_, &dgIntegrationMethod);

  // Get the time integration scheme being used (USED WITH SWITCH)
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
          m_cflFactor = {{{{2.497559, 2.390137, 1.978760, 1.652832, 1.413574, 1.230469, 1.097168, 0.984375, 0.890625,
                            0.812988, 0.747070, 0.691406, 0.643066, 0.602051, 0.565430, 0.533203}},
                          {{2.976074, 2.574463, 2.055664, 1.697998, 1.444092, 1.251221, 1.123047, 1.003418, 0.904541,
                            0.827637, 0.758057, 0.700684, 0.650635, 0.607910, 0.570068, 0.537109}}}};
        } else {
          // Legendre-Gauss-Lobatto nodes
          m_cflFactor = {{{{5.916748, 5.916748, 3.868408, 2.753906, 2.116699, 1.722412, 1.458984, 1.259766, 1.107422,
                            0.990234, 0.895020, 0.815918, 0.750000, 0.692871, 0.644531, 0.602051}},
                          {{7.445068, 7.445068, 4.113770, 2.825928, 2.176514, 1.760254, 1.502686, 1.293945, 1.131592,
                            1.008301, 0.906982, 0.823975, 0.756836, 0.698242, 0.649414, 0.606689}}}};
        }
        break;
    }
    case DG_TIMEINTEGRATION_TOULORGEC_4_8: {
      // Using CFL factors used for CARPENTER 4/5
      if(m_dgIntegrationMethod == "DG_INTEGRATE_GAUSS") {
        // Legendre-Gauss nodes
        m_cflFactor = {{{{2.497559, 2.390137, 1.978760, 1.652832, 1.413574, 1.230469, 1.097168, 0.984375, 0.890625,
                          0.812988, 0.747070, 0.691406, 0.643066, 0.602051, 0.565430, 0.533203}},
                        {{2.976074, 2.574463, 2.055664, 1.697998, 1.444092, 1.251221, 1.123047, 1.003418, 0.904541,
                          0.827637, 0.758057, 0.700684, 0.650635, 0.607910, 0.570068, 0.537109}}}};
      } else {
        // Legendre-Gauss-Lobatto nodes
        m_cflFactor = {{{{5.916748, 5.916748, 3.868408, 2.753906, 2.116699, 1.722412, 1.458984, 1.259766, 1.107422,
                          0.990234, 0.895020, 0.815918, 0.750000, 0.692871, 0.644531, 0.602051}},
                        {{7.445068, 7.445068, 4.113770, 2.825928, 2.176514, 1.760254, 1.502686, 1.293945, 1.131592,
                          1.008301, 0.906982, 0.823975, 0.756836, 0.698242, 0.649414, 0.606689}}}};
      }
      break;
    }
    case DG_TIMEINTEGRATION_NIEGEMANN_4_14: {
      if(m_dgIntegrationMethod == "DG_INTEGRATE_GAUSS") {
        // Legendre-Gauss nodes
        m_cflFactor = {{{{9.669311, 9.232117, 7.577270, 6.309753, 5.397094, 4.675781, 4.187561, 3.739929, 3.382751,
                          3.085876, 2.836548, 2.625488, 2.443419, 2.283385, 2.144226, 2.021301}},
                        {{13.695678, 10.427733, 7.956481, 6.546325, 5.566406, 4.791747, 4.311645, 3.826904, 3.447692,
                          3.163574, 2.896850, 2.667236, 2.480529, 2.314696, 2.169739, 2.042175}}}};
      } else {
        // Legendre-Gauss-Lobatto nodes
        m_cflFactor = {{{{18.998840, 18.998840, 12.540649, 9.887329, 7.777893, 6.397888, 5.569885, 4.792907, 4.213073,
                          3.758484, 3.392029, 3.090515, 2.840026, 2.625488, 2.442260, 2.283385}},
                        {{18.998840, 18.998840, 17.312683, 11.183837, 8.501526, 6.825806, 5.831970, 5.017882, 4.361510,
                          3.876770, 3.475524, 3.154296, 2.889892, 2.662597, 2.473572, 2.307739}}}};
      }
      break;
    }
    case DG_TIMEINTEGRATION_NIEGEMANN_4_13: {
      // Using CFL factors used for CARPENTER 4/5
      if(m_dgIntegrationMethod == "DG_INTEGRATE_GAUSS") {
        // Legendre-Gauss nodes
        m_cflFactor = {{{{2.497559, 2.390137, 1.978760, 1.652832, 1.413574, 1.230469, 1.097168, 0.984375, 0.890625,
                          0.812988, 0.747070, 0.691406, 0.643066, 0.602051, 0.565430, 0.533203}},
                        {{2.976074, 2.574463, 2.055664, 1.697998, 1.444092, 1.251221, 1.123047, 1.003418, 0.904541,
                          0.827637, 0.758057, 0.700684, 0.650635, 0.607910, 0.570068, 0.537109}}}};
      } else {
        // Legendre-Gauss-Lobatto nodes
        m_cflFactor = {{{{5.916748, 5.916748, 3.868408, 2.753906, 2.116699, 1.722412, 1.458984, 1.259766, 1.107422,
                          0.990234, 0.895020, 0.815918, 0.750000, 0.692871, 0.644531, 0.602051}},
                        {{7.445068, 7.445068, 4.113770, 2.825928, 2.176514, 1.760254, 1.502686, 1.293945, 1.131592,
                          1.008301, 0.906982, 0.823975, 0.756836, 0.698242, 0.649414, 0.606689}}}};
      }
      break;
    }
    case DG_TIMEINTEGRATION_TOULORGEC_3_7: {
      // Using CFL factors used for CARPENTER 4/5
      if(m_dgIntegrationMethod == "DG_INTEGRATE_GAUSS") {
        // Legendre-Gauss nodes
        m_cflFactor = {{{{2.497559, 2.390137, 1.978760, 1.652832, 1.413574, 1.230469, 1.097168, 0.984375, 0.890625,
                          0.812988, 0.747070, 0.691406, 0.643066, 0.602051, 0.565430, 0.533203}},
                        {{2.976074, 2.574463, 2.055664, 1.697998, 1.444092, 1.251221, 1.123047, 1.003418, 0.904541,
                          0.827637, 0.758057, 0.700684, 0.650635, 0.607910, 0.570068, 0.537109}}}};
      } else {
        // Legendre-Gauss-Lobatto nodes
        m_cflFactor = {{{{5.916748, 5.916748, 3.868408, 2.753906, 2.116699, 1.722412, 1.458984, 1.259766, 1.107422,
                          0.990234, 0.895020, 0.815918, 0.750000, 0.692871, 0.644531, 0.602051}},
                        {{7.445068, 7.445068, 4.113770, 2.825928, 2.176514, 1.760254, 1.502686, 1.293945, 1.131592,
                          1.008301, 0.906982, 0.823975, 0.756836, 0.698242, 0.649414, 0.606689}}}};
      }
      break;
    }
    case DG_TIMEINTEGRATION_TOULORGEF_4_8: {
      // Using CFL factors used for CARPENTER 4/5
      if(m_dgIntegrationMethod == "DG_INTEGRATE_GAUSS") {
        // Legendre-Gauss nodes
        m_cflFactor = {{{{2.497559, 2.390137, 1.978760, 1.652832, 1.413574, 1.230469, 1.097168, 0.984375, 0.890625,
                          0.812988, 0.747070, 0.691406, 0.643066, 0.602051, 0.565430, 0.533203}},
                        {{2.976074, 2.574463, 2.055664, 1.697998, 1.444092, 1.251221, 1.123047, 1.003418, 0.904541,
                          0.827637, 0.758057, 0.700684, 0.650635, 0.607910, 0.570068, 0.537109}}}};
      } else {
        // Legendre-Gauss-Lobatto nodes
        m_cflFactor = {{{{5.916748, 5.916748, 3.868408, 2.753906, 2.116699, 1.722412, 1.458984, 1.259766, 1.107422,
                          0.990234, 0.895020, 0.815918, 0.750000, 0.692871, 0.644531, 0.602051}},
                        {{7.445068, 7.445068, 4.113770, 2.825928, 2.176514, 1.760254, 1.502686, 1.293945, 1.131592,
                          1.008301, 0.906982, 0.823975, 0.756836, 0.698242, 0.649414, 0.606689}}}};
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
 *
 * Note: If you (for some reason) use nodeVariables information to calculate the
 * initial condition for the conservative variables, you also *must* set the
 * nodeVariables here, otherwise the calcErrorNorms() method of the DG solver
 * procudes faulty results.
 */
template <MInt nDim>
void DgSysEqnAcousticPerturb<nDim>::calcInitialCondition(const MFloat t,
                                                         const MFloat* x,
                                                         MFloat* const nodeVars,
                                                         MFloat* u) const {
  /// \property
  /// \page dgICApe List of initial conditions for solvertype `MAIA_DISCONTINUOUS_GALERKIN` and
  /// `DG_SYSEQN_ACOUSTICPERTURB` List of all available `initialCondition` in
  /// DgSysEqnAcousticPerturb<nDim>::calcInitialCondition():
  switch(this->m_initialCondition) {
    case 0: /// \dgSwitchCase{dgICApe, 0, zero}
    case 1: /// \dgSwitchCase{dgICApe, 1, constant}
    {
      const MFloat c = (this->m_initialCondition == 0) ? 0.0 : F1B2;
      std::fill_n(u, s_noVariables, c);

      // Initialize mean velocities, mean density, mean speed of sound and its
      // derivatives
      std::copy_n(std::begin(m_meanVelocity), nDim, &nodeVars[CV::UU0[0]]);
      nodeVars[CV::RHO0] = m_meanDensity;
      nodeVars[CV::C0] = m_meanSpeedOfSound;
      std::fill_n(&nodeVars[CV::DC0[0]], nDim, 0.0);
      break;
    }

    case 2: {
      /// \dgSwitchCase{dgICApe, 2, zero, mean flow in x-direction with boundary layer (2D)}
      std::fill_n(u, s_noVariables, 0.0);

      // freestream mean velocity in x-direction
      const MFloat freestreamMa = 0.3;
      const MFloat yPos = x[1];
      MFloat u0 = 0.0;
      // boundary layer profile from y=0 to y=1
      // u0 = 0.3*(2*y - 2*y^3 + y^4)
      if(yPos < 1.0 && yPos > 0.0) {
        u0 = freestreamMa * (2 * yPos - 2 * pow(yPos, 3) + pow(yPos, 4));
      } else if(yPos >= 1.0) {
        u0 = freestreamMa;
      }
      nodeVars[CV::UU0[0]] = u0;
      nodeVars[CV::UU0[1]] = 0.0;

      // Mean density, mean speed of sound and its derivatives
      nodeVars[CV::RHO0] = 1.0;
      nodeVars[CV::C0] = 1.0;
      std::fill_n(&nodeVars[CV::DC0[0]], nDim, 0.0);

      break;
    }

    case 3: {
      /// \dgSwitchCase{dgICApe, 3, zero, mean flow in x-direction}
      std::fill_n(u, s_noVariables, 0.0);

      nodeVars[CV::UU0[0]] = m_meanVelocity[0];
      nodeVars[CV::UU0[1]] = m_meanVelocity[1];

      // Mean density, mean speed of sound and its derivatives
      nodeVars[CV::RHO0] = 1.0;
      nodeVars[CV::C0] = 1.0;
      std::fill_n(&nodeVars[CV::DC0[0]], nDim, 0.0);

      break;
    }

    case 5: {
      /// \dgSwitchCase{dgICApe, 5, convergence test}
      /// This convergence test uses the following initial condition:<br>
      /// \f$ \omega = 2 \pi \dfrac{\text{frequency}}{\text{length}} \f$ <br>
      /// \f$ \text{init} = c + A \sin[\omega (-a t + x_1 + x_2 + x_3)] \f$ <br>
      /// \f$ u, v, w = \text{init} \f$ <br>
      /// \f$ p = \text{init}*\text{init}\f$<br>

      // Set domain-specific values and user-defined properties
      const MFloat frequency = this->m_initialNumberWaves; // = rate of change / oscillations
      const MFloat A = F1B5;                               // = max. amplitude of the i.c.
      const MFloat length = F2;                            // = charactersitic length
      const MFloat a = F1B2;                               // = advection velocity
      const MFloat c = F2;                                 // = constant offset

      // Calculate initialization value
      const MFloat omega = F2 * PI * frequency / length;
      const MFloat init = c + A * sin(omega * std::accumulate(x, x + nDim, -a * t));

      // Intialize state variables
      for(MInt i = 0; i < nDim; i++) {
        u[CV::UU[i]] = init;
      }
      u[CV::P] = init * init;

      // Initialize mean velocities, mean density, mean speed of sound and its
      // derivatives
      std::copy_n(std::begin(m_meanVelocity), nDim, &nodeVars[CV::UU0[0]]);
      nodeVars[CV::RHO0] = m_meanDensity;
      nodeVars[CV::C0] = m_meanSpeedOfSound;
      std::fill_n(&nodeVars[CV::DC0[0]], nDim, 0.0);

      break;
    }

    case 22: {
      /// \dgSwitchCase{dgICApe, 22, ...}
      /// Set freestream quantities based on given Mach number (meanVelocity = M_inf)
      MFloat mach = 0.0;
      // Note: mean velocity is used here as free stream mach number M=u/c_inf, not u/c_0
      for(MInt i = 0; i < nDim; i++) {
        mach += POW2(m_meanVelocity[i]);
      }
      mach = sqrt(mach);
      const MFloat kappa = 1.4;

      // Freestream values non-dim with stagnation state
      const MFloat T_inf = 1.0 / (1.0 + (kappa - 1.0) * 0.5 * mach * mach);
      const MFloat c_inf = sqrt(T_inf);
      const MFloat rho_inf = pow(T_inf, 1.0 / (kappa - 1.0));

      std::fill_n(u, s_noVariables, 0.0); // zero IC

      nodeVars[CV::RHO0] = rho_inf;
      nodeVars[CV::C0] = c_inf;
      nodeVars[CV::UU0[0]] = m_meanVelocity[0] * c_inf;
      nodeVars[CV::UU0[1]] = m_meanVelocity[1] * c_inf;
      IF_CONSTEXPR(nDim == 3) { nodeVars[CV::UU0[2]] = m_meanVelocity[2] * c_inf; }
      std::fill_n(&nodeVars[CV::DC0[0]], nDim, 0.0);

      break;
    }

    case 100: {
      /// \dgSwitchCase{dgICApe, 100, First CAA Workshop Category 3 Problem 1}
      /// \f$ p = \exp \left[ - \dfrac{\log 2}{9} \left( x_1^2 + x_2^2 \right) \right] \f$ <br>
      /// \f$ u, v = 0 \f$ <br>

      u[CV::U] = 0.0;
      u[CV::V] = 0.0;
      // u[ CV::U ] = 0.04 * x[1] * exp( -0.04 * log(2.0) *
      // (x[0] * x[0] + x[1]* x[1]) );
      // u[ CV::V ] = -0.04 * ( x[0] - 67.0) *
      // exp( -0.04 * log(2.0) * (x[0] * x[0] + x[1] * x[1] ) );

      u[CV::P] = exp(-log(2.0) / 9.0 * (x[0] * x[0] + x[1] * x[1]));

      // Initialize mean velocities, mean density, mean speed of sound and its
      // derivatives
      std::copy_n(std::begin(m_meanVelocity), nDim, &nodeVars[CV::UU0[0]]);
      nodeVars[CV::RHO0] = m_meanDensity;
      nodeVars[CV::C0] = m_meanSpeedOfSound;
      std::fill_n(&nodeVars[CV::DC0[0]], nDim, 0.0);
      break;
    }

    case 101: {
      /// \dgSwitchCase{dgICApe, 101, Gauss pulse for sponge test}
      /// \f$  p = 2.0 - \exp \left[ -4.0 \left( x_1^2 + x_2^2 \right) \right] \f$ <br>
      /// \f$ u, v = 0 \f$ <br>

      u[CV::U] = 0.0;
      u[CV::V] = 0.0;
      u[CV::P] = 2.0 - exp(-1.0 / 0.25 * (x[0] * x[0] + x[1] * x[1]));

      // Initialize mean velocities, mean density, mean speed of sound and its
      // derivatives
      std::copy_n(std::begin(m_meanVelocity), nDim, &nodeVars[CV::UU0[0]]);
      nodeVars[CV::RHO0] = m_meanDensity;
      nodeVars[CV::C0] = m_meanSpeedOfSound;
      std::fill_n(&nodeVars[CV::DC0[0]], nDim, 0.0);
      break;
    }

    case 102: {
      /// \dgSwitchCase{dgICApe, 102, Gauss pulse for sponge test modified}
      /// \f$ p = 2.0 - \exp \left[ -4.0 \left( { \left( x_1 + 7.5 \right) }^{2} + x_2^2 \right) \right] \f$ <br>
      /// \f$ u, v = 0 \f$ <br>

      u[CV::U] = 0.0;
      u[CV::V] = 0.0;
      u[CV::P] = 2.0 - exp(-1.0 / 0.25 * ((x[0] + 7.5) * (x[0] + 7.5) + x[1] * x[1]));

      // Initialize mean velocities, mean density, mean speed of sound and its
      // derivatives
      std::copy_n(std::begin(m_meanVelocity), nDim, &nodeVars[CV::UU0[0]]);
      nodeVars[CV::RHO0] = m_meanDensity;
      nodeVars[CV::C0] = m_meanSpeedOfSound;
      std::fill_n(&nodeVars[CV::DC0[0]], nDim, 0.0);
      break;
    }

    case 105: {
      /// \dgSwitchCase{dgICApe, 105, First CAA Workshop Category 3 Problem 1 in 3D}
      /// \f$  p = 2.0 - \exp \left[ -4.0 \left( x_1^2 + x_2^2 \right) \right] \f$ <br>
      /// \f$ {uu}_i = 0, \quad i=0..\text{ndim} \f$ <br>

      // Reset velocities
      for(MInt i = 0; i < nDim; i++) {
        u[CV::UU[i]] = 0.0;
      }
      // u[ CV::U ] = 0.04 * x[1] * exp( -0.04 * log(2.0) *
      // (x[0] * x[0] + x[1]* x[1]) );
      // u[ CV::V ] = -0.04 * ( x[0] - 67.0) *
      // exp( -0.04 * log(2.0) * (x[0] * x[0] + x[1] * x[1] ) );

      u[CV::P] = 2.0 - exp(-1.0 / 2.0 * (x[0] * x[0] + x[1] * x[1] + x[2] * x[2]));

      // Initialize mean velocities, mean density, mean speed of sound and its
      // derivatives
      std::copy_n(std::begin(m_meanVelocity), nDim, &nodeVars[CV::UU0[0]]);
      nodeVars[CV::RHO0] = m_meanDensity;
      nodeVars[CV::C0] = m_meanSpeedOfSound;
      std::fill_n(&nodeVars[CV::DC0[0]], nDim, 0.0);
      break;
    }

    case 106: {
      /// \dgSwitchCase{dgICApe, 106, Gauss pulse in 3D}
      /// \f$  p = \exp \left[ -0.5 \left( x_1^2 + x_2^2 + x_3^2 \right) \right] \f$ <br>
      /// \f$ {uu}_i = 0, \quad i=0..\text{ndim} \f$ <br>

      // Reset velocities
      for(MInt i = 0; i < nDim; i++) {
        u[CV::UU[i]] = 0.0;
      }

      u[CV::P] = exp(-1.0 / 2.0 * (x[0] * x[0] + x[1] * x[1] + x[2] * x[2]));

      // Initialize mean velocities, mean density, mean speed of sound and its
      // derivatives
      std::copy_n(std::begin(m_meanVelocity), nDim, &nodeVars[CV::UU0[0]]);
      nodeVars[CV::RHO0] = m_meanDensity;
      nodeVars[CV::C0] = m_meanSpeedOfSound;
      std::fill_n(&nodeVars[CV::DC0[0]], nDim, 0.0);
      break;
    }

    case 107: {
      /// \dgSwitchCase{dgICApe, 107, First CAA Workshop Category 3 Problem 1 in 3D modified}
      /// \f$  p = 2.0 - \exp \left[ -0.5 \left( { \left( x_1 + 7.5 \right) }^{2} + x_2^2 + x_3^2 \right) \right] \f$
      /// <br> \f$ {uu}_i = 0, \quad i=0..\text{ndim} \f$ <br>

      // Reset velocities
      for(MInt i = 0; i < nDim; i++) {
        u[CV::UU[i]] = 0.0;
      }

      u[CV::P] = 2.0 - exp(-1.0 / 2.0 * ((x[0] + 7.5) * (x[0] + 7.5) + x[1] * x[1] + x[2] * x[2]));

      // Initialize mean velocities, mean density, mean speed of sound and its
      // derivatives
      std::copy_n(std::begin(m_meanVelocity), nDim, &nodeVars[CV::UU0[0]]);
      nodeVars[CV::RHO0] = m_meanDensity;
      nodeVars[CV::C0] = m_meanSpeedOfSound;
      std::fill_n(&nodeVars[CV::DC0[0]], nDim, 0.0);
      break;
    }

    case 115: {
      /// \dgSwitchCase{dgICApe, 115, Spherical pressure pulse}
      /// \f$ r = { x_1^2 + x_2^2 + \left( x_3 + 45.0 \right) }^{2} \f$ <br>
      /// \f$ p = \exp \left[ 1.0 + 0.28 \exp \left( \dfrac{-r}{2.5} \right) \right] - \exp [1.0]  \f$ <br>
      /// \f$ {uu}_i = 0, \quad i=0..\text{ndim} \f$ <br>

      // Reset velocities
      for(MInt i = 0; i < nDim; i++) {
        u[CV::UU[i]] = 0.0;
      }

      // Set pressure distribution
      const MFloat r = x[0] * x[0] + x[1] * x[1] + (x[2] + 45.0) * (x[2] + 45.0);
      u[CV::P] = exp(1 + 0.28 * exp(-r / (2.5))) - exp(1.0);

      // Initialize mean velocities, mean density, mean speed of sound and its
      // derivatives
      std::copy_n(std::begin(m_meanVelocity), nDim, &nodeVars[CV::UU0[0]]);
      nodeVars[CV::RHO0] = m_meanDensity;
      nodeVars[CV::C0] = m_meanSpeedOfSound;
      std::fill_n(&nodeVars[CV::DC0[0]], nDim, 0.0);
      break;
    }

    case 201: {
      /// \dgSwitchCase{dgICApe, 201, Spherical pressure pulse for wall reflection}
      /// \f$ p= \exp \left[ - \dfrac{\log 2.0}{25} \left( {x_1^2+ \left( x_2-25.0 \right) }^{2} \right) \right] \f$
      /// <br> \f$ {uu}_i = 0, \quad i=0..\text{ndim} \f$ <br>

      // Reset velocities
      for(MInt i = 0; i < nDim; i++) {
        u[CV::UU[i]] = 0.0;
      }

      // Set pressure distribution
      u[CV::P] = exp(-1.0 * log(2.0) / 25 * (x[0] * x[0] + (x[1] - 25.0) * (x[1] - 25.0)));

      // Initialize mean velocities, mean density, mean speed of sound and its
      // derivatives
      std::copy_n(std::begin(m_meanVelocity), nDim, &nodeVars[CV::UU0[0]]);
      nodeVars[CV::RHO0] = m_meanDensity;
      nodeVars[CV::C0] = m_meanSpeedOfSound;
      std::fill_n(&nodeVars[CV::DC0[0]], nDim, 0.0);
      break;
    }

    case 202: {
      /// \dgSwitchCase{dgICApe, 202, 1D sine wave}
      /// \f$ p = \sin (\pi x_0) \f$ <br>
      /// \f$ u = \sin (\pi x_0) \f$ <br>
      /// \f$ v = 0 \f$ <br>

      u[CV::U] = sin(PI * x[0]);
      u[CV::V] = 0.0;
      u[CV::P] = sin(PI * x[0]);

      // Initialize mean velocities, mean density, mean speed of sound and its
      // derivatives
      std::copy_n(std::begin(m_meanVelocity), nDim, &nodeVars[CV::UU0[0]]);
      nodeVars[CV::RHO0] = m_meanDensity;
      nodeVars[CV::C0] = m_meanSpeedOfSound;
      std::fill_n(&nodeVars[CV::DC0[0]], nDim, 0.0);
      break;
    }


    case 790: {
      /// \dgSwitchCase{dgICApe, 790, ...}
      /// Spinning vortices and other coupled simulations with zero
      /// initial condition

      // Initial condition for setting everything to zero
      // Averaged quantities are read from the given mean file
      std::fill_n(u, s_noVariables, 0.0);

      // Set constant speed of sound and reset its derivatives
      if(m_constantSpeedOfSound) {
        nodeVars[CV::C0] = m_meanSpeedOfSound;
        for(MInt i = 0; i < nDim; i++) {
          nodeVars[CV::DC0[i]] = 0.0;
        }
      }
      break;
    }

    case 1661: {
      /// \dgSwitchCase{dgICApe, 1661, ...}
      /// test condition for extension of nodevars feature, setting nodevars inside a box,
      /// tanh profile U0
      std::fill_n(u, s_noVariables, 0.0); // zero initial condition

      // Set default mean velocitites everywhere
      std::copy_n(std::begin(m_meanVelocity), nDim, &nodeVars[CV::UU0[0]]);

      // Initialize mean velocities, mean density, mean speed of sound and its derivatives inside a
      // box, use zero/rho0/c0 values outside
      MBool pointInBox = true;
      for(MInt dim = 0; dim < nDim; dim++) {
        if(fabs(x[dim]) > 0.75) {
          pointInBox = false;
        }
      }
      if(pointInBox) {
        const MFloat uj = 0.9 * m_meanSpeedOfSound;
        const MFloat r0 = 0.1;
        const MFloat delOmega = 0.05 * r0;
        MFloat rTemp = 0.0;
        for(MInt dim = 1; dim < nDim; dim++) {
          rTemp += x[dim] * x[dim];
        }
        const MFloat r = sqrt(rTemp);
        nodeVars[CV::U0] = uj * (0.5 + 0.5 * tanh((r0 - r) / (2.0 * delOmega)));
        nodeVars[CV::V0] = 0.3 * (1.0 + 0.1 * x[1]);
        nodeVars[CV::RHO0] = m_meanDensity * (1.0 + 0.1 * x[0]);
        nodeVars[CV::C0] = m_meanSpeedOfSound * (1.0 + 0.1 * x[1]);

        std::fill_n(&nodeVars[CV::DC0[0]], nDim, 0.1);
        nodeVars[CV::DC0[0]] = 0.01;
      } else {
        std::fill_n(&nodeVars[CV::UU0[0]], nDim, 0.0);
        nodeVars[CV::RHO0] = m_meanDensity;
        nodeVars[CV::C0] = m_meanSpeedOfSound;
        std::fill_n(&nodeVars[CV::DC0[0]], nDim, 0.0);
      }
      break;
    }

    case 1662: {
      /// \dgSwitchCase{dgICApe, 1662, ...}
      /// test condition for extension of nodevars feature, set nodevars for x>0.5, tanh
      /// profile for U0
      std::fill_n(u, s_noVariables, 0.0); // zero initial condition

      const MFloat uj = 0.9 * m_meanSpeedOfSound;
      const MFloat r0 = 0.1;
      const MFloat delOmega = 0.05 * r0;
      MFloat rTemp = 0.0;
      for(MInt dim = 1; dim < nDim; dim++) {
        rTemp += x[dim] * x[dim];
      }
      const MFloat r = sqrt(rTemp);

      // Set default node variables
      std::copy_n(std::begin(m_meanVelocity), nDim, &nodeVars[CV::UU0[0]]);
      nodeVars[CV::RHO0] = m_meanDensity;
      nodeVars[CV::C0] = m_meanSpeedOfSound;
      std::fill_n(&nodeVars[CV::DC0[0]], nDim, 0.0);

      // Initialize mean velocities, mean density, mean speed of sound and its derivatives for
      // x > 0.5
      if(x[0] >= 0.5) {
        nodeVars[CV::U0] = uj * (0.5 + 0.5 * tanh((r0 - r) / (2.0 * delOmega)));
        nodeVars[CV::C0] = m_meanSpeedOfSound * (1.0 + 0.1 * x[1]);
        nodeVars[CV::DC0[0]] = 0.01;
        // Set density only for y > 0.5
        if(x[1] >= 0.5) {
          nodeVars[CV::RHO0] = m_meanDensity * (1.0 + 0.1 * x[0]);
        }
      }

      break;
    }
    default:
      mTerm(1, AT_, "The specified initial condition (" + std::to_string(this->m_initialCondition) + ") is not valid!");
  }
}


/**
 * \brief Calculates the physical fluxes in all dimensions for all integrations
 *        points within a cell.
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
 * \date 2013-04-18
 *
 * \param[in] q Pointer to solution variables (size: noNodes1D^nDim*noVars)).
 * \param[in] noNodes1D Number of nodes 1D in the cell.
 * \param[out] flux Calculated flux (size: noNodes1D^nDim*noVars*nDim)).
 */
template <MInt nDim>
void DgSysEqnAcousticPerturb<nDim>::calcFlux(const MFloat* const nodeVars,
                                             const MFloat* const q,
                                             const MInt noNodes1D,
                                             MFloat* const flux) const {
  // TRACE();

  const MInt noNodes1D3 = (nDim == 3) ? noNodes1D : 1;
  const MFloatTensor U(const_cast<MFloat*>(q), noNodes1D, noNodes1D, noNodes1D3, s_noVariables);
  const MFloatTensor U0(const_cast<MFloat*>(nodeVars), noNodes1D, noNodes1D, noNodes1D3, s_noNodeVars);
  MFloatTensor f(flux, noNodes1D, noNodes1D, noNodes1D3, nDim, s_noVariables);

  // The following flux equations are calculated here (if 2D, consider all
  // z-components to be zero; for an incompressible case: cBar=1, rhoBar=1):

  // Flux in x-direction:
  // f_x(CV[u]) = u*uBar + v*vBar + w*wBar + p/rhoBar
  // f_x(CV[v]) = 0
  // f_x(CV[w]) = 0
  // f_x(CV[p]) = rhoBar*cBar^2*u + uBar*p

  // Flux in y-direction:
  // f_y(CV[u]) = 0
  // f_y(CV[v]) = u*uBar + v*vBar + w*wBar + p/rhoBar
  // f_y(CV[w]) = 0
  // f_y(CV[p]) = rhoBar*cBar^2*v + vBar*p

  // Flux in z-direction:
  // f_z(CV[u]) = 0
  // f_z(CV[v]) = 0
  // f_z(CV[w]) = u*uBar + v*vBar + w*wBar + p/rhoBar
  // f_z(CV[p]) = rhoBar*cBar^2*w + wBar*p

  for(MInt i = 0; i < noNodes1D; i++) {
    for(MInt j = 0; j < noNodes1D; j++) {
      for(MInt k = 0; k < noNodes1D3; k++) {
        for(MInt d = 0; d < nDim; d++) {
          // Momentum fluxes
          for(MInt fd = 0; fd < nDim; fd++) {
            f(i, j, k, d, CV::UU[fd]) = 0.0;
          }
          f(i, j, k, d, CV::UU[d]) = std::inner_product(&U0(i, j, k, CV::UU0[0]),
                                                        &U0(i, j, k, CV::UU0[0]) + nDim,
                                                        &U(i, j, k, CV::UU[0]),
                                                        U(i, j, k, CV::P) / U0(i, j, k, CV::RHO0));

          // Energy flux
          f(i, j, k, d, CV::P) = U0(i, j, k, CV::RHO0) * POW2(U0(i, j, k, CV::C0)) * U(i, j, k, CV::UU[d])
                                 + U0(i, j, k, CV::UU0[d]) * U(i, j, k, CV::P);
        }
      }
    }
  }
}


/**
 * \brief Calculates the physical fluxes in dimension `dirId` for all
 *        integrations points on an element face.
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
 * \date 2013-04-18
 *
 * \param[in] q Pointer to solution variables
 *              (size: noNodes1D^(nDim-1)*noVars)).
 * \param[in] noNodes1D Number of nodes 1D in the cell.
 * \param[in] dirId Orientation of the face (0 - x, 1 - y, 2 - z).
 * \param[out] flux Calculated flux (size: noNodes1D^(nDim-1)*noVars)).
 */
template <MInt nDim>
void DgSysEqnAcousticPerturb<nDim>::calcFlux1D(const MFloat* const nodeVars,
                                               const MFloat* const q,
                                               const MInt noNodes1D,
                                               const MInt dirId,
                                               MFloat* const flux) const {
  // TRACE();

  const MInt noNodes1D3 = (nDim == 3) ? noNodes1D : 1;
  const MFloatTensor U(const_cast<MFloat*>(q), noNodes1D, noNodes1D3, s_noVariables);
  const MFloatTensor U0(const_cast<MFloat*>(nodeVars), noNodes1D, noNodes1D3, s_noNodeVars);
  MFloatTensor f(flux, noNodes1D, noNodes1D3, s_noVariables);

  for(MInt i = 0; i < noNodes1D; i++) {
    for(MInt j = 0; j < noNodes1D3; j++) {
      // Momentum fluxes
      for(MInt d = 0; d < nDim; d++) {
        if(d == dirId) {
          f(i, j, CV::UU[d]) = std::inner_product(&U0(i, j, CV::UU0[0]),
                                                  &U0(i, j, CV::UU0[0]) + nDim,
                                                  &U(i, j, CV::UU[0]),
                                                  U(i, j, CV::P) / U0(i, j, CV::RHO0));
        } else {
          f(i, j, CV::UU[d]) = 0.0;
        }
      }

      // Energy flux
      f(i, j, CV::P) = U0(i, j, CV::RHO0) * POW2(U0(i, j, CV::C0)) * U(i, j, CV::UU[dirId])
                       + U0(i, j, CV::UU0[dirId]) * U(i, j, CV::P);
    }
  }
}


/**
 *  Yield an empirically-derived maximum stable CFL number.
 *  The user should be able to specify a CFL number of "1.0" in virtually all cases.
 *
 * \author Rodrigo Miguez (rodrigo) rodrigo.miguez@rwth-aachen.de
 * \date 2018-02-28
 *
 * \param[in] polyDeg Polynomial degree in the cell.
 * \param[out] cfl factor obtained from m_cflFactor.
 *
 */
template <MInt nDim>
MFloat DgSysEqnAcousticPerturb<nDim>::cflFactor(const MInt polyDeg) const {
  if(polyDeg > s_maxPolyDeg) {
    TERMM(1, "Polynomial degree exceeds maximum supported");
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
 * \param[in] noNodes1D in the cell.
 * \param[in] t Current time.
 * \param[in] x Coordinates of the integration points of the current cell.
 * \param[out] src Pointer to where source terms are stored
 *                 (size: noNodes1D^nDim*noVars)).
 */
template <MInt nDim>
void DgSysEqnAcousticPerturb<nDim>::calcSource(const MFloat* const nodeVars, const MFloat* const u,
                                               const MInt noNodes1D, const MFloat t, const MFloat* const x,
                                               MFloat* const src) const {
  const MInt noNodes1D3 = (nDim == 3) ? noNodes1D : 1;
  MFloatTensor s(src, noNodes1D, noNodes1D, noNodes1D3, s_noVariables);
  const MFloatTensor X(const_cast<MFloat*>(x), noNodes1D, noNodes1D, noNodes1D3, nDim);
  const MFloatTensor U(const_cast<MFloat*>(u), noNodes1D, noNodes1D, noNodes1D3, s_noVariables);
  const MFloatTensor U0(const_cast<MFloat*>(nodeVars), noNodes1D, noNodes1D, noNodes1D3, s_noNodeVars);

  /// \property
  /// \page dgSourceApe List of source terms for solvertype `MAIA_DISCONTINUOUS_GALERKIN` and
  /// `DG_SYSEQN_ACOUSTICPERTURB` List of all available `sourceTerm` in DgSysEqnAcousticPerturb<nDim>::calcSource():
  switch(this->m_sourceTerm) {
    case 0: {
      /// \dgSwitchCase{dgSourceApe, 0, no source term}
      s.set(F0);
      break;
    }

    case 1: {
      /// \dgSwitchCase{dgSourceApe, 1, APE-1}
      mTerm(1, AT_, "The APE-1 sources are not yet implemented!");
      break;
    }

    case 2: {
      /// \dgSwitchCase{dgSourceApe, 2, APE-2}
      mTerm(1, AT_, "The APE-2 sources are not yet implemented!");
      break;
    }

    case 3: {
      /// \dgSwitchCase{dgSourceApe, 3, APE-3}
      mTerm(1, AT_, "The APE-3 sources are not yet implemented!");
      break;
    }

    case 4: {
      /// \dgSwitchCase{dgSourceApe, 4, APE-4}
      mTerm(1, AT_, "The APE-4 sources are not yet implemented!");
      break;
    }

    case 5: {
      /// \dgSwitchCase{dgSourceApe, 5, Convergence test}
      /// For this convergence test, these source terms are implemented:
      /// \f$
      /// u, v, w = A \omega \cos[\omega (-a t + x1 + x2 + x3)] *
      ///           (-a + 2 c + \mta{u1} + \mta{u2} + \mta{u3} +
      ///            2 A \sin[\omega (-a t + x1 + x2 + x3)])
      /// \f$\f$
      /// p = A \omega \cos[\omega (-a t + x1 + x2 + x3)] *
      ///     (3 + 2 c (-a + \mta{u1} + \mta{u2} + \mta{u3}) +
      ///      2 A (-a + \mta{u1} + \mta{u2} + \mta{u3}) *
      ///      \sin[\omega (-a t + x1 + x2 + x3)])
      /// \f$

      // Set domain-specific values and user-defined properties
      const MFloat frequency = this->m_initialNumberWaves; // = rate of change / oscillations
      const MFloat A = F1B5;                               // = max. amplitude of the i.c.
      const MFloat length = F2;                            // = charactersitic length
      const MFloat a = F1B2;                               // = advection velocity
      const MFloat c = F2;                                 // = constant offset

      // Calculate initialization value
      const MFloat omega = F2 * PI * frequency / length;

      for(MInt i = 0; i < noNodes1D; i++) {
        for(MInt j = 0; j < noNodes1D; j++) {
          for(MInt k = 0; k < noNodes1D3; k++) {
            const MFloat sumX = std::accumulate(&X(i, j, k, 0), &X(i, j, k, 0) + nDim, F0);
            const MFloat inner = omega * (sumX - a * t);
            const MFloat sumVelBarMa = std::accumulate(&U0(i, j, k, CV::UU0[0]), &U0(i, j, k, CV::UU0[0]) + nDim, -a);
            const MFloat F2sinA = F2 * sin(inner) * A;
            const MFloat cosA = cos(inner) * A;

            // Set velocity sources
            for(MInt n = 0; n < nDim; n++) {
              s(i, j, k, CV::UU[n]) = omega * cosA * (F2 * c + sumVelBarMa + F2sinA);
            }
            // Set pressure source
            s(i, j, k, CV::P) = omega * cosA * (nDim + F2 * c * sumVelBarMa + sumVelBarMa * F2sinA);
          }
        }
      }
      break;
    }

    case 6:
      /// \dgSwitchCase{dgSourceApe, 6, Artificial source (same as in the APE4Solver)}
      /// periodically oscillating Gaussian-shaped source of the form
      /// \f$p = A \exp(-((x-x0)^2+(y-y0)^2+(z-z0)^2)/(2*r^2)) \sin(2*pi*f*t)\f$
      {
        const MFloat pos[3] = {0.0, 0.0, 0.0};       // Source center
        const MFloat r = 0.1;                        // Radius
        const MFloat A = 1.0;                        // Amplitude
        const MFloat f = this->m_initialNumberWaves; // Frequency
        // Calculate constant factors
        const MFloat factor = A * sin(2 * PI * f * t);
        const MFloat expFactor = -1.0 / (2 * r * r);

        for(MInt i = 0; i < noNodes1D; i++) {
          for(MInt j = 0; j < noNodes1D; j++) {
            for(MInt k = 0; k < noNodes1D3; k++) {
              // Set velocity sources to zero
              std::fill_n(&s(i, j, k, CV::UU[0]), nDim, F0);

              // Compute squared distance between point and source center
              MFloat distSquared = F0;
              for(MInt n = 0; n < nDim; n++) {
                distSquared += POW2(X(i, j, k, n) - pos[n]);
              }

              // Set pressure source
              s(i, j, k, CV::P) = factor * exp(expFactor * distSquared);
            }
          }
        }
        break;
      }

    case 60:
      /// \dgSwitchCase{dgSourceApe, 60, Artificial source (same as in the APE4Solver)}
      /// periodically oscillating Gaussian-shaped source of the form
      /// \f$p = A \exp(-((x-x0)^2+(y-y0)^2+(z-z0)^2)/(2*r^2)) \sin(2*pi*f*t)\f$
      {
        const MFloat pos[3] = {25.0, 0.0, 0.0};
        const MFloat r = 0.25; // Radius
        const MFloat A = 0.1;  // Amplitude
        const MFloat f = 0.5;  // Frequency
        // Calculate constant factors
        const MFloat factor = A * sin(2 * PI * f * t);
        const MFloat expFactor = -1.0 / (2 * r * r);

        for(MInt i = 0; i < noNodes1D; i++) {
          for(MInt j = 0; j < noNodes1D; j++) {
            for(MInt k = 0; k < noNodes1D3; k++) {
              // Set velocity sources to zero
              std::fill_n(&s(i, j, k, CV::UU[0]), nDim, F0);

              // Compute squared distance between point and source center
              MFloat distSquared = F0;
              for(MInt n = 0; n < nDim; n++) {
                distSquared += POW2(X(i, j, k, n) - pos[n]);
              }

              // Set pressure source
              s(i, j, k, CV::P) = factor * exp(expFactor * distSquared);
            }
          }
        }
        break;
      }

    case 7: {
      TERMM(1, "dipole source needs testing");
      /// \dgSwitchCase{dgSourceApe, 7, Dipole}
      /// Dipole source created by two monopole sources of equal strength and opposite phase
      /// and separated by a small distance (monopole sources are similar to case #6)
      /// Ref: Russell, Titlow, Bemmen: "Acoustic monopoles, dipoles and quadrupoles: An experiment
      /// revisited", 1998.
      const MFloat d = 0.01;                       // Separation distance of the two monopole sources
      const MFloat pos1[3] = {-0.5 * d, 0.0, 0.0}; // Source center monopole #1
      const MFloat pos2[3] = {0.5 * d, 0.0, 0.0};  // Source center monopole #2

      const MFloat r = 0.1; // Radius
      const MFloat A = 1.0; // Amplitude
      const MFloat f = 1.0; // Frequency

      // Calculate constant factors
      const MFloat factor1 = A * sin(2 * PI * f * t);
      const MFloat factor2 = A * sin(2 * PI * (f * t + 0.5)); // Opposite phase for monopole #2

      const MFloat expFactor = -1.0 / (2 * r * r);

      for(MInt i = 0; i < noNodes1D; i++) {
        for(MInt j = 0; j < noNodes1D; j++) {
          for(MInt k = 0; k < noNodes1D3; k++) {
            // Set velocity sources to zero
            std::fill_n(&s(i, j, k, CV::UU[0]), nDim, F0);

            // Compute squared distances between point and source centers
            MFloat distSquared1 = F0;
            MFloat distSquared2 = F0;
            for(MInt n = 0; n < nDim; n++) {
              distSquared1 += POW2(X(i, j, k, n) - pos1[n]);
              distSquared2 += POW2(X(i, j, k, n) - pos2[n]);
            }

            // Set pressure source
            s(i, j, k, CV::P) = factor1 * exp(expFactor * distSquared1) + factor2 * exp(expFactor * distSquared2);
          }
        }
      }
      break;
    }

    case 70: {
      /// \dgSwitchCase{dgSourceApe, 70, analytical source term in 2D - S(x,y) for both}
      const MFloat freq = 20.0; // 1.0/30.0;
      const MFloat omega = 2.0 * PI * freq;
      const MFloat eps = 0.5;
      const MFloat alpha = std::log(2.0) / 2.0;
      const MFloat x_s = 125.0;
      const MFloat y_s = 0.0;

      for(MInt i = 0; i < noNodes1D; i++) {
        for(MInt j = 0; j < noNodes1D; j++) {
          for(MInt k = 0; k < noNodes1D3; k++) {
            // Calculate source term
            const MFloat x_norm = X(i, j, k, 0) - x_s;
            const MFloat y_norm = X(i, j, k, 1) - y_s;
            const MFloat f = eps * std::exp(-alpha * (x_norm * x_norm + y_norm * y_norm));

            // Apply source term
            const MFloat source = f * std::sin(omega * t);
            s(i, j, k, CV::U) = source;
            s(i, j, k, CV::V) = source;
            s(i, j, k, CV::P) = 0.0;
          }
        }
      }
      break;
    }

    case 71: {
      /// \dgSwitchCase{dgSourceApe, 71, analytical source term in 2D - S(x,y) only for u}
      const MFloat freq = 1.0 / 5.0; // 20.0;//1.0/30.0;
      const MFloat omega = 2.0 * PI * freq;
      const MFloat eps = 1e-4; // 0.5;
      const MFloat alpha = std::log(2.0) / 2.0;
      const MFloat x_s = 125.0;
      const MFloat y_s = 0.0;

      for(MInt i = 0; i < noNodes1D; i++) {
        for(MInt j = 0; j < noNodes1D; j++) {
          for(MInt k = 0; k < noNodes1D3; k++) {
            // Calculate source term
            const MFloat x_norm = X(i, j, k, 0) - x_s;
            const MFloat y_norm = X(i, j, k, 1) - y_s;
            const MFloat f = eps * std::exp(-alpha * (x_norm * x_norm + y_norm * y_norm));

            // Apply source term
            const MFloat source = f * std::sin(omega * t);
            s(i, j, k, CV::U) = source;
            s(i, j, k, CV::V) = 0.0;
            s(i, j, k, CV::P) = 0.0;
          }
        }
      }
      break;
    }

    case 800: {
      TERMM_IF_COND(nDim == 3, "only useful in 2D");
      /// \dgSwitchCase{dgSourceApe, 800, Spinning vortices - Ewert source terms (2D only)}
      /// This implements the qm equation on page 388 in:
      ///
      /// Ewert, Schroeder: Acoustic perturbation equations based on flow
      ///                   decomposition via source filtering, 2003.

      // Determine position-independent quantities
      const MFloat gamma = 1.0;
      const MFloat sigma = 1.0;
      const MFloat omega = gamma / 4. / PI;
      const MFloat alpha = gamma * gamma / 8 / PI / PI / sigma / sigma;
      const MFloat theta = omega * t;
      const MFloat bx = std::cos(theta);
      const MFloat by = std::sin(theta);

      /*! \property
        \page propertyPageDG DG
        \section filterSlopeWidth
        <code>MFloat filterSlopeWidth</code>\n
        Default: None \n
        Sets the width of the slope filters used in the APE equations.\n
        Possible values are:
        <ul>
          <li>Any floating point value.</li>
        </ul>
        Keywords: <i>DISCONTINUOUS_GALERKIN, APE</i>
      */
      // Read filter properties
      const MFloat filterSlopeWidth = Context::getSolverProperty<MFloat>("filterSlopeWidth", this->m_solverId, AT_);

      std::array<MFloat, nDim> filterRegionMin{};
      std::array<MFloat, nDim> filterRegionMax{};
      for(MInt i = 0; i < 2; i++) {
        filterRegionMin[i] = Context::getSolverProperty<MFloat>("filterRegionMin", this->m_solverId, AT_, i);
        filterRegionMax[i] = Context::getSolverProperty<MFloat>("filterRegionMax", this->m_solverId, AT_, i);
      }

      // Determine position-dependent quantities and apply sources
      for(MInt i = 0; i < noNodes1D; i++) {
        for(MInt j = 0; j < noNodes1D; j++) {
          // Calculate source term
          const MFloat rxPos = X(i, j, 0, 0) - bx;
          const MFloat ryPos = X(i, j, 0, 1) - by;
          const MFloat rPos = std::sqrt(rxPos * rxPos + ryPos * ryPos);
          const MFloat rxNeg = X(i, j, 0, 0) + bx;
          const MFloat ryNeg = X(i, j, 0, 1) + by;
          const MFloat rNeg = std::sqrt(rxNeg * rxNeg + ryNeg * ryNeg);
          const MFloat qmEwert =
              alpha * (std::exp(-rPos * rPos / 2. / sigma / sigma) - std::exp(-rNeg * rNeg / 2. / sigma / sigma));
          const MFloat qmEwertX = qmEwert * std::cos(theta);
          const MFloat qmEwertY = qmEwert * std::sin(theta);

          // Determine filter
          const MFloat filter = maia::filter::slope::cosbox<nDim>(&filterRegionMin[0], &filterRegionMax[0],
                                                                  filterSlopeWidth, &X(i, j, 0, 0));

          // Apply source term
          s(i, j, 0, CV::U) = filter * qmEwertX;
          s(i, j, 0, CV::V) = filter * qmEwertY;
          s(i, j, 0, CV::P) = 0.0;
        }
      }
      break;
    }

    default:
      mTerm(1, AT_, "The specified source term is not valid!");
  }

  // Return if compresible source term should not be considered
  if(!m_compressibleSourceTerm) {
    return;
  }

  // Add compressible source term to p'
  // q_cons = 2 * (rhoBar * cBar * u' + uBar * 1/cBar * p') * grad(cBar)
  for(MInt i = 0; i < noNodes1D; i++) {
    for(MInt j = 0; j < noNodes1D; j++) {
      for(MInt k = 0; k < noNodes1D3; k++) {
        for(MInt dim = 0; dim < nDim; dim++) {
          s(i, j, k, CV::P) += 2
                               * (U0(i, j, k, CV::RHO0) * U0(i, j, k, CV::C0) * U(i, j, k, CV::UU[dim])
                                  + U0(i, j, k, CV::UU0[dim]) / U0(i, j, k, CV::C0) * U(i, j, k, CV::P))
                               * U0(i, j, k, CV::DC0[dim]);
        }
      }
    }
  }
}

/**
 * \brief Calculates the sponge source terms for all integration points within a cell.
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
 * \date 2013-04-18
 *
 * \param[in] nodeVars Pointer to node variables (size: noNodes1D^nDim*noNodeVars)
 * \param[in] u Pointer to solution variables (size: noNodes1D^nDim*noVars)).
 * \param[in] noNodes1D in the cell.
 * \param[in] eta Pointer to eta used for the sponge calculations.
 * \param[out] src Pointer to where source terms are stored
 *                 (size: noNodes1D^nDim*noVars)).
 */
template <MInt nDim>
void DgSysEqnAcousticPerturb<nDim>::calcSpongeSource(const MFloat* NotUsed(nodeVars), const MFloat* u,
                                                     const MInt noNodes1D, const MFloat* eta, MFloat* src) const {
  const MInt noNodes1D3 = (nDim == 3) ? noNodes1D : 1;

  // gamma = heat capacity ratio for monoatomic gas
  MFloat gammaMinusOne = 5.0 / 3.0 - 1.0;
  MFloat FgammaMinusOne = 1.0 / gammaMinusOne;

  const MFloatTensor spongeEta(const_cast<MFloat*>(eta), noNodes1D, noNodes1D, noNodes1D3);
  MFloatTensor uState(const_cast<MFloat*>(u), noNodes1D, noNodes1D, noNodes1D3, s_noVariables);
  MFloatTensor spongeSource(src, noNodes1D, noNodes1D, noNodes1D3, s_noVariables);

  for(MInt i = 0; i < noNodes1D; i++) {
    for(MInt j = 0; j < noNodes1D; j++) {
      for(MInt k = 0; k < noNodes1D3; k++) {
        MFloat deltaP = (m_spongePressureInfy - uState(i, j, k, CV::P)) * FgammaMinusOne;

        spongeSource(i, j, k, CV::P) = m_spongeSigma * spongeEta(i, j, k) * deltaP;
        for(MInt l = 0; l < nDim; l++) {
          spongeSource(i, j, k, CV::UU[l]) = 0.0;
        }
      }
    }
  }
}


/**
 * \brief Calculate the time step for an explicit time stepping scheme for a
 *        given element.
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
 * \date 2013-04-18
 *
 * \param[in] u Pointer to element variables.
 * \param[in] noNodes1D Current number of nodes 1D.
 * \param[in] invJacobian Inverse Jacobian of element.
 * \return Calculated time step for the element.
 */
template <MInt nDim>
MFloat DgSysEqnAcousticPerturb<nDim>::getTimeStep(const MFloat* const nodeVars,
                                                  const MFloat* const NotUsed(u),
                                                  const MInt noNodes1D,
                                                  const MFloat invJacobian,
                                                  const MInt sbpMode) const {
  // TRACE();

  const MInt noNodes1D3 = (nDim == 3) ? noNodes1D : 1;
  const MFloat inf = std::numeric_limits<MFloat>::infinity();
  const MFloatTensor U0(const_cast<MFloat*>(nodeVars), noNodes1D, noNodes1D, noNodes1D3, s_noNodeVars);

  MFloat maxLambda[] = {-inf, -inf, -inf};
  for(MInt i = 0; i < noNodes1D; i++) {
    for(MInt j = 0; j < noNodes1D; j++) {
      for(MInt k = 0; k < noNodes1D3; k++) {
        for(MInt d = 0; d < nDim; d++) {
          maxLambda[d] = std::max(maxLambda[d], fabs(U0(i, j, k, CV::UU0[d])) + U0(i, j, k, CV::C0));
        }
      }
    }
  }

  MFloat dt;
  if(sbpMode) {
    dt = this->cfl() * F2 / (invJacobian * std::accumulate(&maxLambda[0], &maxLambda[0] + nDim, F0) * (noNodes1D - 1));
  } else {
    const MInt polyDeg = noNodes1D - 1;
    dt = this->cflScaled(polyDeg) * cflFactor(polyDeg) * F2
         / (invJacobian * std::accumulate(&maxLambda[0], &maxLambda[0] + nDim, F0));
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
 * \param[in] nodeVarsL Surface node variables on the left side of the surface
 *                      (i.e. in -x/-y/-z-direction).
 * \param[in] nodeVarsR Surface node variables on the right side of the surface
 *                      (i.e. in +x/+y/+z-direction).
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
void DgSysEqnAcousticPerturb<nDim>::calcRiemann(const MFloat* nodeVarsL,
                                                const MFloat* nodeVarsR,
                                                const MFloat* stateL,
                                                const MFloat* stateR,
                                                const MInt noNodes1D,
                                                const MInt dirId,
                                                MFloat* flux) const {
  // TRACE();

  // Solve Riemann problem
  switch(this->m_riemannSolver) {
    case 0: // local Lax-Friedrichs flux
    {
      calcRiemannLaxFriedich(nodeVarsL, nodeVarsR, stateL, stateR, noNodes1D, dirId, flux);
      break;
    }

    case 1: // Roe's flux
    {
      TERMM(1, "Implementation does not work for internal fluxes, see comments");
      calcRiemannRoe(nodeVarsL, nodeVarsR, stateL, stateR, noNodes1D, dirId, flux);
      break;
    }

    default:
      mTerm(1, AT_, "The specified Riemann solver is not valid!");
      break;
  }
}


/// \brief Calculate Riemann flux using the local Lax-Friedrichs flux scheme.
///
/// \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
/// \date 2015-03-13
///
/// Note: for argument description, see calcRieman(...).
template <MInt nDim>
void DgSysEqnAcousticPerturb<nDim>::calcRiemannLaxFriedich(const MFloat* nodeVarsL,
                                                           const MFloat* nodeVarsR,
                                                           const MFloat* stateL,
                                                           const MFloat* stateR,
                                                           const MInt noNodes1D,
                                                           const MInt dirId,
                                                           MFloat* flux) const {
  // TRACE();

  const MInt noNodes1D3 = (nDim == 3) ? noNodes1D : 1;
  const MFloatTensor uL(const_cast<MFloat*>(stateL), noNodes1D, noNodes1D3, s_noVariables);
  const MFloatTensor uR(const_cast<MFloat*>(stateR), noNodes1D, noNodes1D3, s_noVariables);

  const MFloatTensor nL(const_cast<MFloat*>(nodeVarsL), noNodes1D, noNodes1D3, s_noNodeVars);
  const MFloatTensor nR(const_cast<MFloat*>(nodeVarsR), noNodes1D, noNodes1D3, s_noNodeVars);

#ifndef _OPENMP
  MFloatScratchSpace maxLambda(noNodes1D, noNodes1D3, AT_, "maxLambda");
#else
  MFloatTensor maxLambda(noNodes1D, noNodes1D3);
#endif

  // Calculate maximum eigenvalue
  for(MInt i = 0; i < noNodes1D; i++) {
    for(MInt j = 0; j < noNodes1D3; j++) {
      maxLambda(i, j) = std::max(fabs(nL(i, j, CV::UU0[dirId])) + nL(i, j, CV::C0),
                                 fabs(nR(i, j, CV::UU0[dirId])) + nR(i, j, CV::C0));
    }
  }

#ifndef _OPENMP
  MFloatScratchSpace fluxL(noNodes1D, noNodes1D3, s_noVariables, AT_, "fluxL");
  MFloatScratchSpace fluxR(noNodes1D, noNodes1D3, s_noVariables, AT_, "fluxR");
#else
  MFloatTensor fluxL(noNodes1D, noNodes1D3, s_noVariables);
  MFloatTensor fluxR(noNodes1D, noNodes1D3, s_noVariables);
#endif

  // Calculate flux from left and right state
  calcFlux1D(nodeVarsL, stateL, noNodes1D, dirId, &fluxL[0]);
  calcFlux1D(nodeVarsR, stateR, noNodes1D, dirId, &fluxR[0]);

  // Solve Riemann problem
  MFloatTensor riemann(flux, noNodes1D, noNodes1D3, s_noVariables);
  for(MInt i = 0; i < noNodes1D; i++) {
    for(MInt j = 0; j < noNodes1D3; j++) {
      for(MInt n = 0; n < s_noVariables; n++) {
        riemann(i, j, n) = 0.5 * ((fluxL(i, j, n) + fluxR(i, j, n)) - maxLambda(i, j) * (uR(i, j, n) - uL(i, j, n)));
      }
    }
  }
}


/*---2D Fluxes in Primitive form---
 * The fluxes computed at surface, we could derive
 * inward and outward waves compositions via characteristic analysis,
 * and enforce inward outwaves to zero.
 *
 * Ref. Bauer's paper:
 * "-Application of a Discontinuous Galerkin Method to Discretize
 * Acoustic Perturbation Equations, AIAA, May, 2011
 *
 * \author H.-J. Cheng  <h-j.cheng@aia.rwth-aachen.de>
 * \date 2014-12-25
 *
 * labels:DG IMPORTANT: Please note that the mean flow matrix calculation is just a hack
 * right now and only works if this method is used for a boundary condition, not
 * for internal fluxes!
 */
template <>
inline void DgSysEqnAcousticPerturb<2>::calcRiemannRoe(const MFloat* nodeVarsL,
                                                       const MFloat* nodeVarsR,
                                                       const MFloat* stateL,
                                                       const MFloat* stateR,
                                                       const MInt noNodes1D,
                                                       const MInt dirId,
                                                       MFloat* flux) const {
  // TRACE();

  const MInt nDim = 2;
  const MFloatTensor uL(const_cast<MFloat*>(stateL), noNodes1D, noVars());
  const MFloatTensor uR(const_cast<MFloat*>(stateR), noNodes1D, noVars());
  const MFloatTensor nL(const_cast<MFloat*>(nodeVarsL), noNodes1D, noNodeVars());
  const MFloatTensor nR(const_cast<MFloat*>(nodeVarsR), noNodes1D, noNodeVars());
  MFloatTensor f(flux, noNodes1D, s_noVariables);

  // TODO labels:DG use node variables instead
  // Define mean flow constants
  const MFloat c = 1.0;
  const MFloat rho = 1.0;

  // Normal direction
  const std::array<MInt, nDim> n = {{1 - dirId, dirId}};

  // Characteristic boundary conditions
  for(MInt i = 0; i < noNodes1D; i++) {
    // FIXME labels:DG,totest The following is just a hack and probably does not work when the
    // Roe Riemann solver is used for internal fluxes
    // Determine mean flow
    const std::array<MFloat, nDim> meanVelocity = {
        {0.5 * (nL(i, CV::U0) + nR(i, CV::U0)), 0.5 * (nL(i, CV::V0) + nR(i, CV::V0))}};

    // Mean flow matrix with simplified transformation matrix
    const std::array<MFloat, nDim> um = {
        {meanVelocity[0] * n[0] + meanVelocity[1] * n[1], -meanVelocity[0] * n[1] + meanVelocity[1] * n[0]}};

    // stateL(vecL) and stateR(vecR) variables multiplied by transformation
    // matrix (see reference)
    const std::array<MFloat, nDim> vecL = {
        {uL(i, CV::UU[0]) * n[0] + uL(i, CV::UU[1]) * n[1], -uL(i, CV::UU[0]) * n[1] + uL(i, CV::UU[1]) * n[0]}};
    const std::array<MFloat, nDim> vecR = {
        {uR(i, CV::UU[0]) * n[0] + uR(i, CV::UU[1]) * n[1], -uR(i, CV::UU[0]) * n[1] + uR(i, CV::UU[1]) * n[0]}};

    // Calculate characteristic waves
    // L1 and L2 composite from vecL and vecR
    const MFloat L1 = 0.5 * (uR(i, CV::P) + rho * c * (-vecR[0] + um[1] / (um[0] - c) * (-vecR[1])));
    const MFloat L2 = 0.5 * (uL(i, CV::P) + rho * c * (vecL[0] + um[1] / (um[0] + c) * vecL[1]));

    // Calculate Riemann flux
    f(i, CV::P) = (um[0] - c) * (L1) + (um[0] + c) * (L2);
    f(i, CV::UU[0]) = 1 / rho / c * (-1 * ((um[0] - c) * L1) + ((um[0] + c) * L2)) * n[0];
    f(i, CV::UU[1]) = 1 / rho / c * (-1 * ((um[0] - c) * L1) + ((um[0] + c) * L2)) * n[1];
  }
}


/*---3D Fluxes in Primitive form---
 * The fluxes computed at surface, we could derive
 * inward and outward waves compositions via characteristic analysis,
 * and enforce inward outwaves to zero.
 * Ref. Bauer Paper.
 *
 * Ref. Bauer's paper:
 * "-Application of a Discontinuous Galerkin Method to Discretize
 * Acoustic Perturbation Equations, AIAA, May, 2011
 *
 * \author H.-J. Cheng  <h-j.cheng@aia.rwth-aachen.de>
 * \date 2014-12-25
 *
 * labels:DG IMPORTANT: Please note that the mean flow matrix calculation is just a hack
 * right now and only works if this method is used for a boundary condition, not
 * for internal fluxes!
 */
template <>
inline void DgSysEqnAcousticPerturb<3>::calcRiemannRoe(const MFloat* nodeVarsL,
                                                       const MFloat* nodeVarsR,
                                                       const MFloat* stateL,
                                                       const MFloat* stateR,
                                                       const MInt noNodes1D,
                                                       const MInt dirId,
                                                       MFloat* flux) const {
  // TRACE();

  const MInt nDim = 3;
  const MFloatTensor uL(const_cast<MFloat*>(stateL), noNodes1D, noNodes1D, noVars());
  const MFloatTensor uR(const_cast<MFloat*>(stateR), noNodes1D, noNodes1D, noVars());
  const MFloatTensor nL(const_cast<MFloat*>(nodeVarsL), noNodes1D, noNodes1D, noNodeVars());
  const MFloatTensor nR(const_cast<MFloat*>(nodeVarsR), noNodes1D, noNodes1D, noNodeVars());
  MFloatTensor f(flux, noNodes1D, noNodes1D, s_noVariables);

  // TODO labels:DG use node variables instead
  // Define mean flow constants
  const MFloat c = 1.0;
  const MFloat rho = 1.0;

  // Normal direction
  const std::array<MInt, nDim> n = {{dirId == 0 ? 1 : 0, dirId == 1 ? 1 : 0, dirId == 2 ? 1 : 0}};

  // Characteristic boundary conditions
  // stateL and stateR multiplied by transformation matrix (see reference)
  for(MInt i = 0; i < noNodes1D; i++) {
    for(MInt j = 0; j < noNodes1D; j++) {
      // FIXME labels:DG,totest The following is just a hack and probably does not work when the
      // Roe Riemann solver is used for internal fluxes
      // Determine mean flow
      const std::array<MFloat, nDim> meanVelocity = {{0.5 * (nL(i, j, CV::U0) + nR(i, j, CV::U0)),
                                                      0.5 * (nL(i, j, CV::V0) + nR(i, j, CV::V0)),
                                                      0.5 * (nL(i, j, CV::W0) + nR(i, j, CV::W0))}};

      // Mean flow matrix with transformation matrix (see ref.)
      const std::array<MFloat, nDim> um = {{meanVelocity[0] * n[0] + meanVelocity[1] * n[1] + meanVelocity[2] * n[2],
                                            meanVelocity[0] * n[1] + meanVelocity[1] * n[2] + meanVelocity[2] * n[0],
                                            meanVelocity[0] * n[2] + meanVelocity[1] * n[0] + meanVelocity[2] * n[1]}};

      const std::array<MFloat, nDim> vecL = {
          {uL(i, j, CV::UU[0]) * n[0] + uL(i, j, CV::UU[1]) * n[1] + uL(i, j, CV::UU[2]) * n[2],
           uL(i, j, CV::UU[0]) * n[1] + uL(i, j, CV::UU[1]) * n[2] + uL(i, j, CV::UU[2]) * n[0],
           uL(i, j, CV::UU[0]) * n[2] + uL(i, j, CV::UU[1]) * n[0] + uL(i, j, CV::UU[2]) * n[1]}};

      const std::array<MFloat, nDim> vecR = {
          {uR(i, j, CV::UU[0]) * n[0] + uR(i, j, CV::UU[1]) * n[1] + uR(i, j, CV::UU[2]) * n[2],
           uR(i, j, CV::UU[0]) * n[1] + uR(i, j, CV::UU[1]) * n[2] + uR(i, j, CV::UU[2]) * n[0],
           uR(i, j, CV::UU[0]) * n[2] + uR(i, j, CV::UU[1]) * n[0] + uR(i, j, CV::UU[2]) * n[1]}};

      // Calculate characteristic waves
      // L1 and L2 composite from vecL and vecR
      const MFloat L1 =
          0.5
          * (uR(i, j, CV::P)
             + rho * c * (-vecL[0] + um[1] / (um[0] - c) * (-vecL[1]) + um[2] / (um[0] - c) * (-vecL[2])));
      const MFloat L2 =
          0.5
          * (uL(i, j, CV::P) + rho * c * (vecR[0] + um[1] / (um[0] + c) * vecR[1] + um[2] / (um[0] + c) * (vecR[2])));
      // Calculate Riemann flux
      f(i, j, CV::P) = (um[0] - c) * (L1) + (um[0] + c) * (L2);
      f(i, j, CV::UU[0]) = 1 / rho / c * (-1 * ((um[0] - c) * L1) + ((um[0] + c) * L2)) * n[0];
      f(i, j, CV::UU[1]) = 1 / rho / c * (-1 * ((um[0] - c) * L1) + ((um[0] + c) * L2)) * n[1];
      f(i, j, CV::UU[2]) = 1 / rho / c * (-1 * ((um[0] - c) * L1) + ((um[0] + c) * L2)) * n[2];
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
void DgSysEqnAcousticPerturb<nDim>::consToPrim(const MFloat* cons, MFloat* prim) const {
  // Very easy...
  std::copy(cons, cons + s_noVariables, prim);
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
void DgSysEqnAcousticPerturb<nDim>::primToCons(const MFloat* prim, MFloat* cons) const {
  // Very easy...
  std::copy(prim, prim + s_noVariables, cons);
}


/// \brief Return the default node variables.
///
/// \author Ansgar Niemoeller (ansgar) <a.niemoeller@aia.rwth-aachen.de>
/// \date 2016-07-30
///
/// \param[out] nodeVars Contains the default node variables ordered by CV.
template <MInt nDim>
void DgSysEqnAcousticPerturb<nDim>::getDefaultNodeVars(MFloat* const nodeVars) const {
  // Mean velocities
  std::copy_n(std::begin(m_meanVelocity), nDim, &nodeVars[CV::UU0[0]]);
  // Mean density
  nodeVars[CV::RHO0] = m_meanDensity;
  // Mean speed of sound
  nodeVars[CV::C0] = m_meanSpeedOfSound;
  // Derivatives of mean speed of sound, always zero
  std::fill_n(&nodeVars[CV::DC0[0]], nDim, 0.0);
}


/// \brief Return the default node variables inside a body.
template <MInt nDim>
void DgSysEqnAcousticPerturb<nDim>::getDefaultNodeVarsBody(MFloat* const nodeVars) const {
  // Mean velocities
  std::fill_n(&nodeVars[CV::UU0[0]], nDim, 0.0);
  // Mean density
  nodeVars[CV::RHO0] = m_meanDensity;
  // Mean speed of sound
  nodeVars[CV::C0] = m_meanSpeedOfSound;
  // Derivatives of mean speed of sound, always zero
  std::fill_n(&nodeVars[CV::DC0[0]], nDim, 0.0);
}


/// \brief Return if the given node variable should be extended
template <MInt nDim>
MBool DgSysEqnAcousticPerturb<nDim>::extendNodeVar(const MInt varId) const {
  // Do not extend...
  if(varId == CV::C0) return false; // mean speed of sound ...
  for(MInt i = 0; i < nDim; i++) {
    if(varId == CV::DC0[i]) { // and its derivatives
      return false;
    }
  }
  return true;
}


#endif // DGSYSEQNACOUSTICPERTURB_H_
