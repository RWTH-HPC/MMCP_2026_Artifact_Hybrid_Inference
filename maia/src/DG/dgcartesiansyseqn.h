// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef DGSYSEQN_H_
#define DGSYSEQN_H_

#include <vector>
#include "INCLUDE/maiatypes.h"
#include "IO/context.h"
#include "UTIL/debug.h"
#include "UTIL/functions.h"
#include "property.h"


/**
 * \brief Base class for concrete system-of-equations classes (CRTP).
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
 * \date 2013-08-08
 *
 * \tparam nDim Number of spatial dimensions
 * \tparam SysEqn System-of-equations class
 */
template <MInt nDim, class SysEqn>
class DgSysEqn {
  // CRTP
 private:
  SysEqn& derived() { return static_cast<SysEqn&>(*this); }
  const SysEqn& derived() const { return static_cast<const SysEqn&>(*this); }
  MFloat m_cfl;

  // Methods
 public:
  static constexpr MInt noVars() { return SysEqn::s_noVariables; }
  static constexpr MInt noNodeVars() { return SysEqn::s_noNodeVars; }
  static constexpr MBool hasTimeDependentNodeVars() { return SysEqn::s_hasTimeDependentNodeVars; }
  static const MString& sysEqnName() { return SysEqn::s_sysEqnName; }
  static const MString& consVarNames(MInt i) { return SysEqn::s_consVarNames[i]; }
  static const MString& primVarNames(MInt i) { return SysEqn::s_primVarNames[i]; }
  static const MString& nodeVarNames(MInt i) { return SysEqn::s_nodeVarNames[i]; }
  MFloat cfl() const { return m_cfl; }

  MFloat cflScaled(const MInt polyDeg) const {
    // cf. Cockburn and Shu, Runge-Kutta Discontinuous Galerkin Methods for
    // Convection-Dominated Problems, J. Sci. Comp. 16(3), 2001.
    return cfl() / (2.0 * polyDeg + 1.0);
  }

 protected:
  explicit DgSysEqn(MInt solverId);


  // Member variables
 public:
  MInt m_initialCondition;
  MFloat m_initialNumberWaves;
  MInt m_sourceTerm;
  MInt m_riemannSolver;


 protected:
  const MInt m_solverId;
};


template <MInt nDim, class SysEqn>
DgSysEqn<nDim, SysEqn>::DgSysEqn(MInt solverId) : m_solverId(solverId) {
  TRACE();

  // Read and set initial condition
  m_initialCondition = Context::getSolverProperty<MInt>("initialCondition", m_solverId, AT_);

  // Read and set initial number of waves in the domain
  m_initialNumberWaves = 2.0;
  m_initialNumberWaves =
      Context::getSolverProperty<MFloat>("initialNumberWaves", m_solverId, AT_, &m_initialNumberWaves);

  // Read and set source term
  /*! \property
      \page propertyPageDG DG
      \section sourceTerm
      <code>MBool DgSysEqn::m_sourceTerm</code>\n
      default = <code>0</code>\n\n
      This property stores which source term should be applied (if any). This
      property might have different meanings in different systems of equations.
      Refer to DgSysEqnXXX::calcSource() to see its meaning.\n\n
      Possible values are:
      <ul>
      <li><code>0</code> (no source term)</li>
      <li><code>1...n</code> Apply source term n</li>
      </ul>
      Keywords: <i>DISCONTINUOUS_GALERKIN</i>
  */
  m_sourceTerm = 0;
  m_sourceTerm = Context::getSolverProperty<MInt>("sourceTerm", m_solverId, AT_, &m_sourceTerm);

  // Read and set Riemann solver term
  /*! \property
      \page propertyPageDG DG
      \section riemannSolver
      <code>MBool DgSysEqn::m_riemannSolver</code>\n
      default = <code>0</code>\n\n
      This property stores which Riemann solver should be used for the numerical
      fluxes. The property might have different meanings for different systems
      of equations. Refer to DgSysEqnXXX::calcRiemann() to see its meaning.
      n\n
      Possible values are:
      <ul>
      <li><code>0</code> Use local Lax-Friedrichs flux (should be implemented
          for each system of equations)</li>
      <li><code>1...n</code> Use Riemann solver n</li>
      </ul>
      Keywords: <i>DISCONTINUOUS_GALERKIN, RIEMANN</i>
  */
  m_riemannSolver = 0;
  m_riemannSolver = Context::getSolverProperty<MInt>("riemannSolver", m_solverId, AT_, &m_riemannSolver);


  // Read and set CFL number for the time step calculation
  m_cfl = Context::getSolverProperty<MFloat>("cfl", m_solverId, AT_);
}

#endif // DGSYSEQN_H_
