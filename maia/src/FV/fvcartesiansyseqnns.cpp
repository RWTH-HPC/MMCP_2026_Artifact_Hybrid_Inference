// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "fvcartesiansyseqnns.h"

template <MInt nDim>
FvSysEqnNS<nDim>::FvSysEqnNS(const MInt solverId, const MInt noSpecies) : m_solverId(solverId) {
  CV = new ConservativeVariables(noSpecies);
  PV = new PrimitiveVariables(noSpecies);
  FV = new FluxVariables(noSpecies);

  readProperties();
}

template <MInt nDim>
FvSysEqnNS<nDim>::ConservativeVariables::ConservativeVariables(const MInt noSpecies)
  : m_noSpecies(noSpecies), noVariables(nDim + 2 + noSpecies) {
  if(m_noSpecies > 0) {
    mAlloc(RHO_Y, m_noSpecies, "FvSysEqnNS::ConservativeVariables::RHO_Y", AT_);
    for(MUint i = 0; i < m_noSpecies; ++i) {
      RHO_Y[i] = RHO_C + i;
    }
  }
}

template <MInt nDim>
FvSysEqnNS<nDim>::FluxVariables::FluxVariables(const MInt noSpecies) : ConservativeVariables(noSpecies) {}

template <MInt nDim>
FvSysEqnNS<nDim>::ConservativeVariables::~ConservativeVariables() {
  mDeallocate(RHO_Y);
}

template <MInt nDim>
FvSysEqnNS<nDim>::PrimitiveVariables::PrimitiveVariables(const MInt noSpecies)
  : m_noSpecies(noSpecies), noVariables(nDim + 2 + noSpecies) {
  if(m_noSpecies > 0) {
    mAlloc(Y, m_noSpecies, "FvSysEqnNS::PrimitiveVariables::Y", AT_);
    for(MUint i = 0; i < m_noSpecies; ++i) {
      Y[i] = C + i;
    }
  }
}

template <MInt nDim>
FvSysEqnNS<nDim>::PrimitiveVariables::~PrimitiveVariables() {
  mDeallocate(Y);
}

template <MInt nDim>
void FvSysEqnNS<nDim>::PrimitiveVariables::getPrimitiveVariableNames(MString* names) {
  TRACE();

  for(MInt i = 0; i < noVariables; i++) {
    names[i] = varNames[i];
  }
}

template <MInt nDim>
void FvSysEqnNS<nDim>::readProperties() {
  m_gamma = Context::getSolverProperty<MFloat>("gamma", m_solverId, AT_, &m_gamma);
  m_gammaMinusOne = m_gamma - 1.0;
  m_F1BGammaMinusOne = 1 / m_gammaMinusOne;

  m_Pr = Context::getSolverProperty<MFloat>("Pr", m_solverId, AT_, &m_Pr);
  m_F1BPr = 1 / m_Pr;

  m_referenceTemperature = 273.15;
  /*! \page propertyPage1
    \section referenceTemperature
    <code>MFloat FvCartesianSolver::m_referenceTemperature </code>\n
    default = <code>273.15</code>\n \n
  Reference temperature \f$ T_{\mathrm{ref}}\f$
  Used to scale the Sutherland's constant as follows: \f$ S/T_{\mathrm{ref}} \f$
  Also used for the computation of the reference sound speed and combustion (TF) related quantities
    possible values are:
    <ul>
    <li>Non-negative floating point values</li>
    </ul>
    Keywords: <i>FINITE_VOLUME, VARIABLES</i>
  */
  m_referenceTemperature =
      Context::getSolverProperty<MFloat>("referenceTemperature", m_solverId, AT_, &m_referenceTemperature);

  m_sutherlandConstant = 110.4;
  /*! \page propertyPage1
    \section sutherlandConstant
    <code>MFloat FvCartesianSolver::m_sutherlandConstant </code>\n
    default = <code>110.4 K</code>\n \n
    Sutherland's constant. Used by Sutherland's law.
    possible values are:
    <ul>
    <li>Non-negative floating point values</li>
    </ul>
    Keywords: <i>FINITE_VOLUME, VARIABLES</i>
  */
  m_sutherlandConstant =
      Context::getSolverProperty<MFloat>("sutherlandConstant", m_solverId, AT_, &m_sutherlandConstant);

  m_sutherlandConstantThermal = m_sutherlandConstant; // default value assumes a constant Prandtl number
  /*! \page propertyPage1
    \section sutherlandConstantThermal
    <code>MFloat FvCartesianSolver::m_sutherlandConstantThermal </code>\n
    default = <code>110.4 K</code>\n \n
    Sutherland's constant for thermal conductivity. Recommended value: 194.0. See 'Viscous Fluid Flow' by F.M. White.
    possible values are:
    <ul>
    <li>Non-negative floating point values</li>
    </ul>
    Keywords: <i>FINITE_VOLUME, VARIABLES</i>
  */
  m_sutherlandConstantThermal =
      Context::getSolverProperty<MFloat>("sutherlandConstantThermal", m_solverId, AT_, &m_sutherlandConstantThermal);

  m_sutherlandConstant /= m_referenceTemperature;
  m_sutherlandPlusOne = m_sutherlandConstant + F1;
  m_sutherlandConstantThermal /= m_referenceTemperature;
  m_sutherlandPlusOneThermal = m_sutherlandConstantThermal + F1;

  m_gFGMOrPr = m_gamma * m_F1BGammaMinusOne * m_F1BPr;

  // Stencil specific properties
  m_enhanceThreePointViscFluxFactor = Context::getSolverProperty<MFloat>("enhanceThreePointViscFluxFactor", m_solverId,
                                                                         AT_, &m_enhanceThreePointViscFluxFactor);

  /*! \page propertyPage1
    \section viscousFluxScheme
    <code>MString FvCartesianSolver::m_viscousFluxScheme </code>\n
    default = <code>FIVE_POINT</code>\n
    Scheme for the calculation of the viscous flux\n
    Possible values are:
    <ul>
      <li>THREE_POINT </li>
      <li>FIVE_POINT </li>
      <li>FIVE_POINT_STABILIZED </li>
    </ul>
    Keywords: <i>FINITE VOLUME, NUMERICS, FLUX</i>
  */
  MString m_viscousFluxScheme = "FIVE_POINT";
  m_viscousFluxScheme = Context::getSolverProperty<MString>("viscousFluxScheme", m_solverId, AT_, &m_viscousFluxScheme);

  /*! \page propertyPage1
    \section viscousFluxScheme
    <code>MFloat FvCartesianSolver::m_enhanceThreePointViscFluxFactor </code>\n
    default = <code>0.1</code>\n
    FIVE_POINT_STABILIZED combines the THREE_POINT stencil and the FIVE_POINT stencil for the viscous flux computation\n
    This property provides further control and the final stencil of the viscous flux is
    (1-enhanceThreePointViscFluxFactor)*FIVE_POINT + * enhanceThreePointViscFluxFactor*THREE_POINT
    \n
    Possible values are:
    <ul>
      Floats between 0 and 1
    </ul>
    Keywords: <i>FINITE VOLUME, NUMERICS, FLUX</i>
  */
  if(string2enum(m_viscousFluxScheme) == FIVE_POINT_STABILIZED
     && Context::propertyExists("enhanceThreePointViscFluxFactor", m_solverId)) {
    m_enhanceThreePointViscFluxFactor = Context::getSolverProperty<MFloat>(
        "enhanceThreePointViscFluxFactor", m_solverId, AT_, &m_enhanceThreePointViscFluxFactor);
  } else {
    m_enhanceThreePointViscFluxFactor = 0.1;
  }

  /*! \page propertiesFV
  \section centralizeSurfaceVariablesFactor
  <code>MFloat FvCartesianSolver::m_centralizeSurfaceVariablesFactor </code>\n
  default = <code>0</code>\n
  upwind factor computation based on this factor
  Keywords: <i>FINITE_VOLUME, AUSM, NUMERICS</i>
  */
  m_centralizeSurfaceVariablesFactor = 0.0;

  m_centralizeSurfaceVariablesFactor = Context::getSolverProperty<MFloat>(
      "centralizeSurfaceVariablesFactor", m_solverId, AT_, &m_centralizeSurfaceVariablesFactor);
  if(m_centralizeSurfaceVariablesFactor > 0.0) {
    m_log << "centralize surface variable factor: " << m_centralizeSurfaceVariablesFactor << std::endl;
  }
}

template <>
constexpr MUint FvSysEqnNS<2>::index0[2]{1, 0};
template <>
constexpr MUint FvSysEqnNS<3>::index0[3]{1, 2, 0};
template <>
constexpr MUint FvSysEqnNS<2>::index1[2]{};
template <>
constexpr MUint FvSysEqnNS<3>::index1[3]{2, 0, 1};

template class FvSysEqnNS<2>;
template class FvSysEqnNS<3>;
