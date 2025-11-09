// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef FvSysEqnDetChem_H_
#define FvSysEqnDetChem_H_

#include <algorithm>
#include <array>
#include <iterator>
#include <memory>
#include <numeric>
#include "INCLUDE/maiaconstants.h"
#include "INCLUDE/maiatypes.h"
#include "MEMORY/alloc.h"
#include "filter.h"
#include "fvcartesiansyseqnns.h"

#if defined(WITH_CANTERA)

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wfloat-equal"
#include "cantera/base/Solution.h"
#include "cantera/thermo/IdealGasPhase.h"
#include "cantera/thermo/Species.h"
#include "cantera/thermo/SpeciesThermoInterpType.h"
#include "cantera/transport.h"
#include "cantera/zerodim.h"
#pragma GCC diagnostic pop
using namespace Cantera;

/**
 * @brief Class containing the special methods of the detailed chemistry formulation. Inherits from the Navier-Stokes
 * equation class.
 *
 * @tparam nDim Number of dimensions
 */
template <MInt nDim>
class FvSysEqnDetChem : public FvSysEqnNS<nDim> {
  // Declare parent a friend so that CRTP can access private methods/members
  // friend class FvSysEqn<nDim>;

 public:
  FvSysEqnDetChem(const MInt solverId, const MInt noSpecies);

  template <MInt stencil = AUSM>
  inline void Ausm(const MInt orientation, const MFloat upwindCoefficient, const MFloat A, const MFloat* const leftVars,
                   const MFloat* const rightVars, const MFloat* const srfcCoeff, MFloat* const flux);

  inline void AusmBndryCorrection(const MInt orientation, const MFloat A, const MFloat* const leftVars,
                                  const MFloat* const rightVars, MFloat* const flux);

  // Dispatch function for classic five-point-style viscous flux schemes
  template <MInt stencil>
  inline void viscousFlux(const MInt orientation, const MFloat A, const MFloat* const vars0, const MFloat* const vars1,
                          const MFloat* const slope0, const MFloat* const slope1, const MFloat* const srfcCoeff,
                          const MFloat f0, const MFloat f1, MFloat* const flux) {
    if(stencil == FIVE_POINT) {
      viscousFluxFivePoint(orientation, A, vars0, vars1, slope0, slope1, srfcCoeff, f0, f1, flux);
    } else {
      TERMM(1, "Viscous Flux Scheme not implemented.");
    }
  }

  // Dispatch function for compact viscous flux schemes
  template <MInt stencil>
  inline void viscousFlux(const MInt, const MFloat, const MBool, const MFloat* const, const MFloat* const,
                          const MFloat* const, const MFloat* const, const MFloat* const, const MFloat* const,
                          const MFloat* const, const MFloat* const, const MFloat* const, const MFloat, const MFloat,
                          MFloat* const) {
    if(stencil == THREE_POINT) {
      TERMM(1, "Three point viscous Flux Scheme not implemented for SysEqnDetChem.");
    } else {
      TERMM(1, "Viscous Flux Scheme not implemented.");
    }
  }

  inline void viscousFluxFivePoint(const MInt orientation, const MFloat A, const MFloat* const vars0,
                                   const MFloat* const vars1, const MFloat* const slope0, const MFloat* const slope1,
                                   const MFloat* const srfcCoeff, const MFloat f0, const MFloat f1, MFloat* const flux);

  inline void getSpeciesDiffusionMassFluxes(const MInt orientation, const MFloat* const vars0,
                                            const MFloat* const vars1, const MFloat* const slope0,
                                            const MFloat* const slope1, const MFloat* const srfcCoeff, const MFloat f0,
                                            const MFloat f1, std::vector<MFloat>& J, std::vector<MFloat>& dXdn,
                                            MFloat& dTdn, const MBool soretEffect);

  inline void computeSurfaceCoefficients(const MInt RKStep, const MFloat* const vars0, const MFloat* const vars1,
                                         MFloat* const srfcCoeff, std::shared_ptr<Cantera::ThermoPhase> gas,
                                         std::shared_ptr<Cantera::Transport> trans);

  inline void computeSpeciesReactionRates(const MFloat& m_timeStep, const MFloat& cellVolume,
                                          const MFloat* const pvarsCell, const MFloat* const avarsCell,
                                          MFloat* const reactionRatesCell, Cantera::IdealGasReactor* zeroD_reactor,
                                          Cantera::ReactorNet* zeroD_reactorNet, std::shared_ptr<Cantera::Solution> sol,
                                          std::shared_ptr<Cantera::ThermoPhase> gas);

  inline void computePrimitiveVariables(const MFloat* const cvarsCell, MFloat* const pvarsCell,
                                        const MFloat* const avarsCell);

  inline void computeConservativeVariables(const MFloat* const pvarsCell, MFloat* const cvarsCell,
                                           const MFloat* const avarsCell);

  inline void evaluateSensibleEnergy(MFloat& sensibleEnergy, const MFloat* const pvarsCell,
                                     const MFloat& meanMolarMass);

  inline void iterateTemperature(MFloat& T, const MFloat* const cvarsCell, const MFloat* const pvarsCell,
                                 const MFloat* const avarsCell, const MFloat& velPOW2);

  inline void computeMeanMolarWeight_PV(const MFloat* const pvarsCell, MFloat* const avarsCell);

  inline void computeMeanMolarWeight_CV(const MFloat* const pvarsCell, MFloat* const avarsCell);

  inline void computeGamma(const MFloat* const pvarsCell, MFloat* const avarsCell,
                           std::shared_ptr<Cantera::ThermoPhase> gas);

  inline MFloat computePhi(const MFloat* const Y);

  inline std::vector<std::vector<MFloat>> conservativeSlopes(const MFloat* const pvarsCell,
                                                             const MFloat* const cvarsCell,
                                                             const MFloat* const avarsCell,
                                                             const MFloat* const slopesCell);

 private:
  using Base = FvSysEqnNS<nDim>;
  using Base::index0;
  using Base::index1;
  using Base::m_noSpecies;
  using Base::m_solverId;

  // Equation specific values
  using Base::m_F1BGammaMinusOne;
  using Base::m_gamma;
  using Base::m_gammaMinusOne;
  using Base::m_gFGMOrPr;
  using Base::m_referenceTemperature;
  using Base::m_sutherlandConstant;
  using Base::m_sutherlandConstantThermal;
  using Base::m_sutherlandPlusOne;
  using Base::m_sutherlandPlusOneThermal;

  using Base::getArray012;

  MFloat m_gasConstant;
  MFloat m_fGasConstant;


  void readProperties();
  void allocateNoDiffusionCoefficients(const MInt);

  // To be read from properties file
  MString m_reactionMechanism;
  MString m_transportModel;
  MString m_phaseName;

  // Computed in allocateNoDiffusionCoefficients(noSpecies)
  MInt m_noSurfaceCoefficients;
  MInt m_noDiffusionCoefficients;
  MInt m_noThermalDiffusionCoefficients;
  MBool m_multiDiffusion;
  MBool m_computeSrfcCoeffsEveryRKStep;

 public:
  using Base::m_muInfinity;
  using Base::m_Re0;

  // Stores species NASA coefficients and physical properties
  struct NASACoefficients;
  struct SpeciesProperties;

  SpeciesProperties* m_species = nullptr;
  NASACoefficients* m_NASA = nullptr;
  // Hold indices for primitive and conservative variables
  struct ConservativeVariables;
  struct FluxVariables;
  struct PrimitiveVariables;
  struct AdditionalVariables;
  struct SurfaceCoefficients;
  static constexpr MBool hasAV = true;
  static constexpr MBool hasSC = true; // add to other sysEqn
  ConservativeVariables* CV = nullptr;
  FluxVariables* FV = nullptr;
  PrimitiveVariables* PV = nullptr;
  AdditionalVariables* AV = nullptr;
  SurfaceCoefficients* SC = nullptr;
  MBool m_soretEffect;

  MString m_oxidizer;
  MString m_fuel;
  MFloat m_fuelOxidizerStochiometricRatio;
};

/** @brief Static indices for accessing conservative variables in nDim spatial dimensions
 */
template <MInt nDim>
struct FvSysEqnDetChem<nDim>::ConservativeVariables {
  static const MInt Segfault = std::numeric_limits<MInt>::min();

  static constexpr MInt RHO_U = 0;
  static constexpr MInt RHO_V = 1;
  static constexpr MInt RHO_W = (nDim == 3) ? 2 : Segfault;
  static constexpr std::array<MInt, nDim> RHO_VV = getArray012();
  static constexpr MInt RHO_E = nDim;
  static constexpr MInt RHO = nDim + 1;

  static constexpr MInt RHO_C = nDim + 2;

  MUint m_noSpecies;
  const MInt noVariables;

  MInt* RHO_Y = nullptr;

  ConservativeVariables(const MInt noSpecies);
  ~ConservativeVariables();
};

/** @brief Static indices for accessing flux variables. In this SysEqn identical to the conservative variables.
 */
template <MInt nDim>
struct FvSysEqnDetChem<nDim>::FluxVariables : ConservativeVariables {
  FluxVariables(const MInt noSpecies);
};

/** @brief Static indices for accessing primitive variables in nDim spatial dimensions
 */
template <MInt nDim>
struct FvSysEqnDetChem<nDim>::PrimitiveVariables {
  static const MInt Segfault = std::numeric_limits<MInt>::min();

  static constexpr MInt U = 0;
  static constexpr MInt V = 1;
  static constexpr MInt W = (nDim == 3) ? 2 : Segfault;
  static constexpr std::array<MInt, nDim> VV = getArray012();
  static constexpr MInt RHO = nDim;
  static constexpr MInt P = nDim + 1;

  static constexpr MInt C = nDim + 2;

  MUint m_noSpecies;
  const MInt noVariables;

  const std::array<MString, nDim + 3> varNames = [] {
    IF_CONSTEXPR(nDim == 2) {
      std::array<MString, nDim + 3> a = {"u", "v", "rho", "p", "c"};
      return a;
    }
    else {
      std::array<MString, nDim + 3> a = {"u", "v", "w", "rho", "p", "c"};
      return a;
    }
  }();

  MInt* Y = nullptr;

  PrimitiveVariables(const MInt noSpecies);
  void getPrimitiveVariableNames(MString* names) {
    for(MInt i = 0; i < noVariables; i++) {
      names[i] = varNames[i];
    }
  };
  ~PrimitiveVariables();
};


/**
 * @brief Struct for storing all relevant species information, which is read from a given mechanism file.
 *
 * @tparam nDim Number of dimensions
 */
template <MInt nDim>
struct FvSysEqnDetChem<nDim>::SpeciesProperties {
  const MString majorSpecies = "N2";
  MInt majorSpeciesIndex;

  const MFloat referenceTemp = 298.0;

  std::map<std::string, MInt> speciesMap;

  std::vector<MString> speciesName;
  MFloat* molarMass = nullptr;
  MFloat* fMolarMass = nullptr;
  MFloat* specificGasConstant = nullptr;
  MFloat* fSpecificGasConstant = nullptr;
  MFloat* standardHeatFormation = nullptr;

  SpeciesProperties(const MInt noSpecies, const FvSysEqnDetChem<nDim>& sysEqn);
  ~SpeciesProperties();

 protected:
  void getSpeciesProperties(const MInt noSpecies, const FvSysEqnDetChem<nDim>& sysEqn);
};

/**
 * @brief Stores all NASA coefficients. NASA coefficients are used to compute the cp and cv values. Additional, modified
 * NASA coefficients are stored, which are used to compute the sensible energy. For now only valid for NASA-7
 * coefficients, which have low temperature and high temperature regions, each one described by a different polynomial.
 *
 * @tparam nDim Number of dimensions
 */
template <MInt nDim>
struct FvSysEqnDetChem<nDim>::NASACoefficients {
  friend class FvSysEqnDetChem<nDim>;
  static constexpr MFloat referenceTemp = 298.0;
  static constexpr MFloat transitionTemp = 1000.0;

  static constexpr MInt totalNumberCoefficientsPerSpecies = 15;
  static constexpr MInt noNASACoefficients = 7;
  static constexpr MInt noNASACoefficientsCpPolynomial = 5;

  MFloat* lowTemp = nullptr;
  MFloat* highTemp = nullptr;

  MFloat* integralLowTemp = nullptr;
  MFloat* integralHighTemp = nullptr;

  MFloat* lowTempIntegrationConstantsEnergy = nullptr;
  MFloat* highTempIntegrationConstantsEnergy = nullptr;

  MFloat* lowTempIntegrationConstantsEnthalpy = nullptr;
  MFloat* highTempIntegrationConstantsEnthalpy = nullptr;

  NASACoefficients(const MInt noSpecies, const FvSysEqnDetChem<nDim>& sysEqn);
  ~NASACoefficients();

 private:
  void getNASACoefficients(const MInt noSpecies, const FvSysEqnDetChem<nDim>& sysEqn);
  void computeSensibleEnergyIntegrationConstants(const FvSysEqnDetChem<nDim>& sysEqn);
};

/**
 * @brief  Static indices for accessing surface coefficients
 *
 * @tparam nDim Number of dimensions
 */
template <MInt nDim>
struct FvSysEqnDetChem<nDim>::SurfaceCoefficients {
  static const MInt Segfault = std::numeric_limits<MInt>::min();

  const MInt m_noDiffusionCoefficients;
  const MInt m_noThermalDiffusionCoefficients;
  const MInt m_noSurfaceCoefficients;

  static constexpr MInt MU = 0;
  static constexpr MInt LAMBDA = 1;
  static constexpr MInt CP = 2;
  static constexpr MInt W_MEAN = 3;
  static constexpr MInt D0 = 4;

  MInt* D = nullptr;
  MInt* DT = nullptr;

  SurfaceCoefficients(const MInt noSpecies, const MInt noDiffusionCoefficients,
                      const MInt noThermalDiffusionCoefficients);
  ~SurfaceCoefficients();
};

/**
 * @brief Static indices for accessing additional variables
 *
 * @tparam nDim Number of dimensions
 */
template <MInt nDim>
struct FvSysEqnDetChem<nDim>::AdditionalVariables {
  static constexpr MInt noVariables = 2;

  static constexpr MInt GAMMA = 0;
  static constexpr MInt W_MEAN = 1;
};


template <MInt nDim>
template <MInt stencil>
inline void FvSysEqnDetChem<nDim>::Ausm(const MInt orientation, const MFloat upwindCoefficient, const MFloat A,
                                        const MFloat* const leftVars, const MFloat* const rightVars,
                                        const MFloat* const srfcCoeff, MFloat* const flux) {
  // Multispecies specific variables and coefficients
  const MFloat CP = srfcCoeff[SC->CP];

  // compute detailed chemistry additional variables
  MFloat fMeanMolarWeightL = 0.0;
  MFloat fMeanMolarWeightR = 0.0;

  for(MUint s = 0; s < PV->m_noSpecies; ++s) {
    fMeanMolarWeightL += leftVars[PV->Y[s]] * m_species->fMolarMass[s];
    fMeanMolarWeightR += rightVars[PV->Y[s]] * m_species->fMolarMass[s];
  }

  const MFloat meanMolarWeightL = F1 / fMeanMolarWeightL;
  const MFloat meanMolarWeightR = F1 / fMeanMolarWeightR;
  const MFloat specificGasConstantL = m_gasConstant * fMeanMolarWeightL;
  const MFloat specificGasConstantR = m_gasConstant * fMeanMolarWeightR;
  const MFloat CVL = CP - specificGasConstantL;
  const MFloat CVR = CP - specificGasConstantR;
  const MFloat gammaL = CP / CVL;
  const MFloat gammaR = CP / CVR;

  // catch the primitive variables rho and p,
  // compute speed of sound, and interface mach number
  const MFloat RHOL = leftVars[PV->RHO];
  const MFloat PL = leftVars[PV->P];
  const MFloat AL = sqrt(gammaL * mMax(MFloatEps, PL / mMax(MFloatEps, RHOL)));
  const MFloat ML = leftVars[orientation] / AL;

  const MFloat RHOR = rightVars[PV->RHO];
  const MFloat PR = rightVars[PV->P];
  const MFloat AR = sqrt(gammaR * mMax(MFloatEps, PR / mMax(MFloatEps, RHOR)));
  const MFloat MR = rightVars[orientation] / AR;

  // calculation of the resulting pressure and mach number on the surface
  const MFloat MLR = 0.5 * (ML + MR);
  const MFloat PLR = PL * (0.5 + upwindCoefficient * ML) + PR * (0.5 - upwindCoefficient * MR);

  // calculation of the left and right rho*a
  const MFloat RHO_AL = RHOL * AL;
  const MFloat RHO_AR = RHOR * AR;

  // calculation of the resulting mass flux through the surface
  const MFloat RHO_U2 = 0.25 * (MLR * (RHO_AL + RHO_AR) + fabs(MLR) * (RHO_AL - RHO_AR));
  const MFloat AbsRHO_U2 = fabs(RHO_U2);

  // velocity magnitudes
  const MFloat U2L = std::inner_product(&leftVars[PV->U], &leftVars[PV->U] + nDim, &leftVars[PV->U], 0.0);
  const MFloat U2R = std::inner_product(&rightVars[PV->U], &rightVars[PV->U] + nDim, &rightVars[PV->U], 0.0);

  ///\TODO labels:FV,totest check sensible energy vs sensible enthalpy
  // Compute left sensible energy value
  MFloat sensibleEnergyL = F0;
  this->evaluateSensibleEnergy(sensibleEnergyL, leftVars, meanMolarWeightL);

  // Compute right sensible energy value
  MFloat sensibleEnergyR = F0;
  this->evaluateSensibleEnergy(sensibleEnergyR, rightVars, meanMolarWeightR);

  const MFloat PLfRHOL = PL / RHOL;
  const MFloat PRfRHOR = PR / RHOR;

  const MFloat e0 = sensibleEnergyL + 0.5 * U2L + PLfRHOL;
  const MFloat e1 = sensibleEnergyR + 0.5 * U2R + PRfRHOR;

  std::array<MFloat, nDim> pFactor{};
  pFactor[orientation] = 1.0;

  for(MUint n = 0; n < nDim; n++) {
    flux[FV->RHO_VV[n]] = (RHO_U2 * (leftVars[PV->VV[n]] + rightVars[PV->VV[n]])
                           + AbsRHO_U2 * (leftVars[PV->VV[n]] - rightVars[PV->VV[n]]) + PLR * pFactor[n])
                          * A;
  }

  flux[FV->RHO_E] = (RHO_U2 * (e0 + e1) + AbsRHO_U2 * (e0 - e1)) * A;
  flux[FV->RHO] = 2.0 * RHO_U2 * A;

  for(MUint s = 0; s < PV->m_noSpecies; ++s) {
    const MFloat YsL = leftVars[PV->Y[s]];
    const MFloat YsR = rightVars[PV->Y[s]];
    flux[FV->RHO_Y[s]] = (RHO_U2 * (YsL + YsR) + AbsRHO_U2 * (YsL - YsR)) * A;
  }
}

template <MInt nDim>
inline void FvSysEqnDetChem<nDim>::AusmBndryCorrection(const MInt orientation, const MFloat A,
                                                       const MFloat* const leftVars, const MFloat* const rightVars,
                                                       MFloat* const flux) {
  const MFloat PL = leftVars[PV->P];
  const MFloat PR = rightVars[PV->P];
  const MFloat PLR = 0.5 * (PR + PL);

  std::array<MFloat, nDim> pFactor{};
  pFactor[orientation] = 1.0;

  for(MUint n = 0; n < nDim; n++) {
    flux[FV->RHO_VV[n]] = PLR * pFactor[n] * A;
  }

  flux[FV->RHO_E] = 0.0;
  flux[FV->RHO] = 0.0;

  // Flux calculation for species transport
  // TODO labels:FV @Julian, Make noSpecies constexpr
  for(MUint s = 0; s < PV->m_noSpecies; ++s) {
    flux[FV->RHO_Y[s]] = 0.0;
  }
}

template <MInt nDim>
inline void FvSysEqnDetChem<nDim>::viscousFluxFivePoint(const MInt orientation, const MFloat A,
                                                        const MFloat* const vars0, const MFloat* const vars1,
                                                        const MFloat* const slope0, const MFloat* const slope1,
                                                        const MFloat* const srfcCoeff, const MFloat f0, const MFloat f1,
                                                        MFloat* const flux) {
  std::vector<MFloat> dXdn(PV->m_noSpecies, F0);
  std::vector<MFloat> J(PV->m_noSpecies, F0);
  std::vector<MFloat> speciesSensibleEnthalpy(PV->m_noSpecies, F0);


  // Constant coefficients defined at the surface from left and right values
  const MFloat mue = srfcCoeff[SC->MU];
  const MFloat lambda = srfcCoeff[SC->LAMBDA];
  const MFloat meanMolarMass = srfcCoeff[SC->W_MEAN];

  MFloat mue_wm = 0.0;

  const std::array<MFloat, nDim> velocity = [&] {
    std::array<MFloat, nDim> tmp{};
    for(MUint n = 0; n < nDim; n++) {
      tmp[n] = F1B2 * (vars0[PV->VV[n]] + vars1[PV->VV[n]]);
    }
    return tmp;
  }();

  const MFloat rho = F1B2 * (vars0[PV->RHO] + vars1[PV->RHO]);
  const MFloat Frho = F1 / rho;
  const MFloat p = F1B2 * (vars0[PV->P] + vars1[PV->P]);

  // Compute the temperature  through equation of state (ideal gas)
  const MFloat T = p * Frho * meanMolarMass * m_fGasConstant;

  // Indices for the orientations
  const MUint id0 = orientation;
  const MUint id1 = index0[orientation];
  const MUint id2 = index1[orientation];

  const MBool isLowTempRegion = (T < m_NASA->transitionTemp);
  const MFloat* const NASAIntegralCoeffs =
      isLowTempRegion ? ALIGNED_F(m_NASA->integralLowTemp) : ALIGNED_F(m_NASA->integralHighTemp);
  const MFloat* const NASAIntegralConstants = isLowTempRegion ? ALIGNED_F(m_NASA->lowTempIntegrationConstantsEnthalpy)
                                                              : ALIGNED_F(m_NASA->highTempIntegrationConstantsEnthalpy);

  // Compute species sensible enthalpy
  for(MUint s = 0; s < PV->m_noSpecies; ++s) {
    const MFloat speciesSpecificGasConstant = m_species->specificGasConstant[s];

    MFloat sensibleEnthalpy = F0;
    MInt offsetIntegralCoeffs = m_NASA->noNASACoefficientsCpPolynomial * s;

    for(MInt i = 4; (i < 5) && (i >= 0); --i) {
      sensibleEnthalpy = sensibleEnthalpy * T + NASAIntegralCoeffs[offsetIntegralCoeffs + i];
    }

    sensibleEnthalpy *= speciesSpecificGasConstant * T;
    sensibleEnthalpy -= NASAIntegralConstants[s];

    speciesSensibleEnthalpy[s] = sensibleEnthalpy;
  }


  MFloat dTdn = F0;
  MFloat diffEnthalpyFlux = F0;

  std::fill(J.begin(), J.end(), F0);
  std::fill(dXdn.begin(), dXdn.end(), F0);

  getSpeciesDiffusionMassFluxes(orientation, vars0, vars1, slope0, slope1, srfcCoeff, f0, f1, J, dXdn, dTdn,
                                m_soretEffect);


  for(MUint s = 0; s < PV->m_noSpecies; ++s) {
    diffEnthalpyFlux += speciesSensibleEnthalpy[s] * J[s];
  }

  const MFloat q = lambda * dTdn - diffEnthalpyFlux;

  // Compute the stress terms
  const MUint s00 = id0 * nDim + id0;
  const MUint s01 = id0 * nDim + id1;
  const MUint s10 = id1 * nDim + id0;
  const MUint s11 = id1 * nDim + id1;

  std::array<MFloat, nDim> tau{};
  IF_CONSTEXPR(nDim == 2) {
    tau[id0] = (mue + mue_wm)
               * (f0 * (F4B3 * slope0[s00] - F2B3 * slope0[s11]) + f1 * (F4B3 * slope1[s00] - F2B3 * slope1[s11]));
    tau[id1] = (mue + mue_wm) * (f0 * (slope0[s01] + slope0[s10]) + f1 * (slope1[s01] + slope1[s10]));
  }
  else IF_CONSTEXPR(nDim == 3) {
    const MUint s22 = id2 * nDim + id2;
    const MUint s02 = id0 * nDim + id2;
    const MUint s20 = id2 * nDim + id0;
    tau[id0] = (mue + mue_wm)
               * (f0 * (F4B3 * slope0[s00] - F2B3 * (slope0[s11] + slope0[s22]))
                  + f1 * (F4B3 * slope1[s00] - F2B3 * (slope1[s11] + slope1[s22])));
    tau[id1] = (mue + mue_wm) * (f0 * (slope0[s01] + slope0[s10]) + f1 * (slope1[s01] + slope1[s10]));
    tau[id2] = (mue + mue_wm) * (f0 * (slope0[s02] + slope0[s20]) + f1 * (slope1[s02] + slope1[s20]));
  }

  // Compute the flux
  for(MUint n = 0; n < nDim; n++) {
    flux[FV->RHO_VV[n]] -= A * tau[n];
  }
  flux[FV->RHO_E] -= A * (std::inner_product(velocity.begin(), velocity.end(), tau.begin(), 0.0) + q);

  // Species transport equations
  ///\TODO labels:FV,totest check if negativ expression is correct
  for(MUint s = 0; s < PV->m_noSpecies; ++s) {
    flux[FV->RHO_Y[s]] -= -A * J[s]; // remove
  }
}

template <MInt nDim>
inline void FvSysEqnDetChem<nDim>::getSpeciesDiffusionMassFluxes(
    const MInt orientation, const MFloat* const vars0, const MFloat* const vars1, const MFloat* const slope0,
    const MFloat* const slope1, const MFloat* const srfcCoeff, const MFloat f0, const MFloat f1, std::vector<MFloat>& J,
    std::vector<MFloat>& dXdn, MFloat& dTdn, const MBool soretEffect) {
  std::vector<MFloat> speciesDiffusionMassFlux(PV->m_noSpecies, F0); // j_k
  std::vector<MFloat> d(PV->m_noSpecies, F0);
  std::vector<MFloat> VX(PV->m_noSpecies, F0);

  const MFloat rho = F1B2 * (vars0[PV->RHO] + vars1[PV->RHO]);
  const MFloat Frho = F1 / rho;
  const MFloat p = F1B2 * (vars0[PV->P] + vars1[PV->P]);
  const MFloat Fp = F1 / p;

  const MFloat meanMolarMass = srfcCoeff[SC->W_MEAN];
  const MFloat fMeanMolarMass = F1 / meanMolarMass;

  // Compute the temperature through equation of state (ideal gas)
  const MFloat T = p * Frho * meanMolarMass * m_fGasConstant;
  const MFloat FT = F1 / T;

  // Indices for the orientations
  const MUint id0 = orientation;

  const MUint sq0 = PV->P * nDim + id0;
  const MUint sq1 = PV->RHO * nDim + id0;

  const MFloat dPdn = f0 * slope0[sq0] + f1 * slope1[sq0];
  const MFloat dRhodn = f0 * slope0[sq1] + f1 * slope1[sq1];

  MFloat summSpeciesGradients = F0;
  for(MUint s = 0; s < PV->m_noSpecies; s++) {
    const MUint sqs = PV->Y[s] * nDim + id0;
    summSpeciesGradients += m_species->fMolarMass[s] * (f0 * slope0[sqs] + f1 * slope1[sqs]);
  }

  const MFloat dWdn = -POW2(srfcCoeff[SC->W_MEAN]) * summSpeciesGradients;

  for(MUint s = 0; s < PV->m_noSpecies; s++) {
    const MUint sqs = PV->Y[s] * nDim + id0;
    const MFloat Ys = F1B2 * (vars0[PV->Y[s]] + vars1[PV->Y[s]]);
    dXdn[s] = m_species->fMolarMass[s] * (meanMolarMass * (f0 * slope0[sqs] + f1 * slope1[sqs]) + dWdn * Ys);
  }

  dTdn = m_fGasConstant * (p * Frho * dWdn + (Frho * meanMolarMass * (dPdn - p * Frho * dRhodn)));

  for(MUint s = 0; s < PV->m_noSpecies; ++s) {
    const MFloat Ys = F1B2 * (vars0[PV->Y[s]] + vars1[PV->Y[s]]);
    d[s] = dXdn[s] + (meanMolarMass * m_species->fMolarMass[s] * Ys - Ys) * Fp * dPdn;
  }

  if(m_multiDiffusion) {
    for(MUint s = 0; s < PV->m_noSpecies; ++s) {
      for(MUint j = 0; j < PV->m_noSpecies; ++j) {
        if(j == s) continue;
        VX[s] += m_species->molarMass[j] * srfcCoeff[SC->D[PV->m_noSpecies * j + s]] * d[j];
      }

      VX[s] *= fMeanMolarMass;

      if(soretEffect) VX[s] -= srfcCoeff[SC->DT[s]] * meanMolarMass * Frho * m_species->fMolarMass[s] * FT * dTdn;

      J[s] = rho * m_species->molarMass[s] * fMeanMolarMass * VX[s];
    }
  } else {
    for(MUint s = 0; s < PV->m_noSpecies; ++s) {
      // TODO labels:FV check
      VX[s] = -srfcCoeff[SC->D[s]] * dXdn[s];

      J[s] = rho * m_species->molarMass[s] * fMeanMolarMass * VX[s];
    }
  }
}

template <MInt nDim>
inline void FvSysEqnDetChem<nDim>::computePrimitiveVariables(const MFloat* const cvarsCell, MFloat* const pvarsCell,
                                                             const MFloat* const avarsCell) {
  const MFloat fRho = F1 / cvarsCell[CV->RHO];

  MFloat velPOW2 = F0;
  for(MInt n = 0; n < nDim; ++n) { // compute velocity
    pvarsCell[PV->VV[n]] = cvarsCell[CV->RHO_VV[n]] * fRho;
    velPOW2 += POW2(pvarsCell[PV->VV[n]]);
  }

  // Density
  pvarsCell[PV->RHO] = cvarsCell[CV->RHO];

  // Compute the species
  for(MUint s = 0; s < PV->m_noSpecies; ++s) {
    pvarsCell[PV->Y[s]] = mMax(F0, mMin(F1, cvarsCell[CV->RHO_Y[s]] * fRho));
  }

  MFloat T = F0;
  iterateTemperature(T, cvarsCell, pvarsCell, avarsCell, velPOW2);

  pvarsCell[PV->P] = pvarsCell[PV->RHO] * T * m_gasConstant / avarsCell[AV->W_MEAN];

#ifndef NDEBUG
  if(isnan(pvarsCell[PV->P]) || isnan(pvarsCell[PV->RHO]) || isnan(T)) {
    std::cout << "FvSysEqnDetChem<nDim>::computePrimitiveVariables has NaN value. T: " << T
              << ", rho: " << pvarsCell[PV->RHO] << ", p: " << pvarsCell[PV->P] << std::endl;
  }

  if(pvarsCell[PV->P] < 0 || pvarsCell[PV->RHO] < 0 || T < 0) {
    std::cout << "FvSysEqnDetChem<nDim>::computePrimitiveVariables has negative value. T: " << T
              << ", rho: " << pvarsCell[PV->RHO] << ", p: " << pvarsCell[PV->P] << std::endl;
  }
#endif
}

/**
 * @brief Iterates the temperature from the sensible energy. Iteration is necessary, since the heat capacities are a
function of temperature. The iteration is performed with a Newton iterative Method.
 *
 * @tparam nDim Number of dimensions
 * @param T Temperature, is given as
 * @param cvarsCell Pointer to the conservative variables at the cell
 * @param pvarsCell Pointer to the primitive variables at the cell
 * @param avarsCell Pointer to the additional variables at the cell
 * @param velPOW2 Absolute value of the velocity squared
 */
template <MInt nDim>
inline void FvSysEqnDetChem<nDim>::iterateTemperature(MFloat& T, const MFloat* const cvarsCell,
                                                      const MFloat* const NotUsed(pvarsCell),
                                                      const MFloat* const avarsCell, const MFloat& velPOW2) {
  const MFloat energyZeroLevel = m_gasConstant * m_species->referenceTemp / avarsCell[AV->W_MEAN];

  MFloat LHS = cvarsCell[CV->RHO_E] / cvarsCell[CV->RHO] - F1B2 * velPOW2 + energyZeroLevel;

  // TODO labels:FV Change to old variables and convergence check
  // Supply initial temperature value to start iteration
  MFloat T_n = 300.0;

  for(MInt it = 0; it < 10; ++it) {
    const MBool isLowTempRegion = (T_n < m_NASA->transitionTemp);
    // Assigns relevant low or high temperature coefficients and integration constants
    const MFloat* const NASACoeffs = isLowTempRegion ? ALIGNED_F(m_NASA->lowTemp) : ALIGNED_F(m_NASA->highTemp);
    const MFloat* const NASAIntegralCoeffs =
        isLowTempRegion ? ALIGNED_F(m_NASA->integralLowTemp) : ALIGNED_F(m_NASA->integralHighTemp);
    const MFloat* const NASAIntegralConstants = isLowTempRegion ? ALIGNED_F(m_NASA->lowTempIntegrationConstantsEnergy)
                                                                : ALIGNED_F(m_NASA->highTempIntegrationConstantsEnergy);

    MFloat T_nPlus1 = F0, f_xi = F0, fS_xi = F0;
    for(MUint s = 0; s < PV->m_noSpecies; s++) {
      const MInt offsetIntegralCoeffs = m_NASA->noNASACoefficientsCpPolynomial * s;
      const MInt offsetNASAPolynomial = m_NASA->noNASACoefficients * s;
      const MFloat Y_s = cvarsCell[CV->RHO_Y[s]] / cvarsCell[CV->RHO];
      const MFloat specificGasConstant = m_species->specificGasConstant[s];

      // Calculate NASA polynomial for species s with the Horner's rule
      MFloat NASAIntegralPolynomial = F0, NASAPolynomial = F0;
      for(MInt i = 4; (i < 5) && (i >= 0); i--) {
        NASAIntegralPolynomial = NASAIntegralPolynomial * T_n + NASAIntegralCoeffs[offsetIntegralCoeffs + i];
        NASAPolynomial = NASAPolynomial * T_n + NASACoeffs[offsetNASAPolynomial + i];
      }

      NASAIntegralPolynomial *= specificGasConstant;
      NASAPolynomial *= specificGasConstant;

      NASAIntegralPolynomial = (NASAIntegralPolynomial - specificGasConstant) * T_n;
      NASAPolynomial -= specificGasConstant;

      f_xi += NASAIntegralPolynomial * Y_s;
      f_xi -= NASAIntegralConstants[s] * Y_s;
      fS_xi += NASAPolynomial * Y_s;
    }

    f_xi -= LHS;
    T_nPlus1 = T_n - (f_xi / fS_xi);
    T_n = T_nPlus1;
  }

  T = T_n;
}

/**
 * @brief Computes the sensible energy. Species sensible energy is computed as: e_s = Int(cv(T)dT) - R*T_ref/W_mean
 *
 * @tparam nDim Number of dimensions
 * @param sensibleEnergy Sensible energy
 * @param pvarsCell Pointer to the primitive variables at the cell
 * @param meanMolarMass Mean molar mass at the cell
 */
template <MInt nDim>
inline void FvSysEqnDetChem<nDim>::evaluateSensibleEnergy(MFloat& sensibleEnergy, const MFloat* const pvarsCell,
                                                          const MFloat& meanMolarMass) {
  const MFloat T = pvarsCell[PV->P] / pvarsCell[PV->RHO] * meanMolarMass * m_fGasConstant;
  const MBool isLowTempRegion = (T < m_NASA->transitionTemp);

  const MFloat* const NASAIntegralCoeffs =
      isLowTempRegion ? ALIGNED_F(m_NASA->integralLowTemp) : ALIGNED_F(m_NASA->integralHighTemp);
  const MFloat* const NASAIntegralConstants = isLowTempRegion ? ALIGNED_F(m_NASA->lowTempIntegrationConstantsEnergy)
                                                              : ALIGNED_F(m_NASA->highTempIntegrationConstantsEnergy);

  for(MUint s = 0; s < PV->m_noSpecies; s++) {
    MInt offsetIntegralCoeffs = m_NASA->noNASACoefficientsCpPolynomial * s;

    // Calculate NASA polynomial for species s with the Horner's rule
    MFloat sensibleEnergySpecies = F0, NASAIntegralPolynomial = F0;
    for(MInt i = 4; (i < 5) && (i >= 0); --i) {
      NASAIntegralPolynomial = NASAIntegralPolynomial * T + NASAIntegralCoeffs[offsetIntegralCoeffs + i];
    }

    NASAIntegralPolynomial *= m_species->specificGasConstant[s];
    NASAIntegralPolynomial = (NASAIntegralPolynomial - m_species->specificGasConstant[s]) * T;
    sensibleEnergySpecies = NASAIntegralPolynomial - NASAIntegralConstants[s];
    sensibleEnergy += sensibleEnergySpecies * pvarsCell[PV->Y[s]];
  }

  // Substract zero level sensible energy
  sensibleEnergy -= m_gasConstant * m_NASA->referenceTemp / meanMolarMass;
}

template <MInt nDim>
inline void FvSysEqnDetChem<nDim>::computeConservativeVariables(const MFloat* const pvarsCell,
                                                                MFloat* const cvarsCell,
                                                                const MFloat* const avarsCell) {
  // Compute rho v
  MFloat velPOW2 = F0;
  for(MInt vel = 0; vel < nDim; ++vel) {
    const MFloat v = pvarsCell[vel];
    cvarsCell[vel] = v * pvarsCell[PV->RHO];
    velPOW2 += POW2(v);
  }

  // Copy the density
  cvarsCell[CV->RHO] = pvarsCell[PV->RHO];

  MFloat sensibleEnergy = F0;
  evaluateSensibleEnergy(sensibleEnergy, pvarsCell, avarsCell[AV->W_MEAN]);

  cvarsCell[CV->RHO_E] = pvarsCell[PV->RHO] * sensibleEnergy + F1B2 * pvarsCell[PV->RHO] * velPOW2;

  // compute rho Yi
  for(MUint s = 0; s < PV->m_noSpecies; ++s) {
    cvarsCell[CV->RHO_Y[s]] = pvarsCell[PV->Y[s]] * pvarsCell[PV->RHO];
  }
}

/**
 * @brief Computes the species reaction rates. A reactor is constructed at each cell, with equal thermodynamic state.
 * The reactor content is then integrated in time and advanced by the same amount as the flow simulation. From the
 * initial and final concentrations of species, a time averaged species reaction rate is computed, which can then be
 * used in the species equations in the flow solver.
 *
 * @tparam nDim
 * @param m_timeStep Current time step of the flow solver
 * @param cellVolume Volume of the computational cell
 * @param pvarsCell Pointer to the primitive variables at the cell
 * @param avarsCell Pointer to the additional variables at the cell
 * @param reactionRatesCell Pointer to the reaction rates at the cell
 * @param zeroD_reactor Pointer to the Cantera::Reactor object
 * @param zeroD_reactorNet Pointer to the Cantera::ReactorNet object
 * @param sol Pointer to the Cantera::Solution object
 * @param gas Pointer to the Cantera::ThermoPhase object
 */
template <MInt nDim>
inline void FvSysEqnDetChem<nDim>::computeSpeciesReactionRates(
    const MFloat& m_timeStep, const MFloat& NotUsed(cellVolume), const MFloat* const pvarsCell,
    const MFloat* const avarsCell, MFloat* const reactionRatesCell, Cantera::IdealGasReactor* zeroD_reactor,
    Cantera::ReactorNet* zeroD_reactorNet, std::shared_ptr<Cantera::Solution> sol,
    std::shared_ptr<Cantera::ThermoPhase> gas) {
  const MFloat p = pvarsCell[PV->P];
  const MFloat rho = pvarsCell[PV->RHO];
  const MFloat meanMolarMass = avarsCell[AV->W_MEAN];
  const MFloat T = p / rho * meanMolarMass * m_fGasConstant;

  gas->setState_TPY(T, p, &pvarsCell[PV->Y[0]]);

  zeroD_reactor->insert(sol);
  zeroD_reactorNet->setInitialTime(0.0);

  const MFloat timeStep = m_timeStep;
  const MFloat reactorDensity_t0 = zeroD_reactor->density();

  zeroD_reactorNet->advance(timeStep);

  auto reactorMassFractions = zeroD_reactor->massFractions();
  const MFloat reactorDensity = zeroD_reactor->density();

  MFloat C_k_t_0, C_k_deltaT;
  for(MUint s = 0; s < PV->m_noSpecies; s++) {
    C_k_t_0 = pvarsCell[PV->Y[s]] * reactorDensity_t0 * m_species->fMolarMass[s];
    C_k_deltaT = reactorMassFractions[s] * reactorDensity * m_species->fMolarMass[s];
    reactionRatesCell[s] = ((C_k_deltaT - C_k_t_0) / (timeStep)) * m_species->molarMass[s];
  }
}

/**
 * @brief Computes the transport coefficients at the surface. These are then used in the computation of the surface
 * diffusion fluxes
 *
 * @tparam nDim Number of dimensions
 * @param RKStep Current Runge-Kutta step
 * @param vars0 Pointer to the primitive variables at one side of the surface
 * @param vars1 Pointer to the primitive variables at the other side of the surface
 * @param srfcCoeff Pointer to the surface coefficients at the surface
 * @param gas Pointer to the Cantera::ThermoPhase object
 * @param trans Pointer to the Cantera::Transport object
 */
template <MInt nDim>
inline void FvSysEqnDetChem<nDim>::computeSurfaceCoefficients(const MInt RKStep, const MFloat* const vars0,
                                                              const MFloat* const vars1, MFloat* const srfcCoeff,
                                                              std::shared_ptr<Cantera::ThermoPhase> gas,
                                                              std::shared_ptr<Cantera::Transport> trans) {
  // Compute mean molar mass and store species vector
  std::vector<MFloat> Y_k(PV->m_noSpecies, F0);
  MFloat fMeanMolarMass = F0;
  for(MUint s = 0; s < PV->m_noSpecies; ++s) {
    Y_k[s] = F1B2 * (vars0[PV->Y[s]] + vars1[PV->Y[s]]);
    fMeanMolarMass += Y_k[s] * m_species->fMolarMass[s];
  }

  srfcCoeff[SC->W_MEAN] = F1 / fMeanMolarMass;

  // Compute the rest of surface coefficients only at RKStep 0
  if(RKStep == 0 || m_computeSrfcCoeffsEveryRKStep) {
    const MFloat p = F1B2 * (vars0[PV->P] + vars1[PV->P]);
    const MFloat rho = F1B2 * (vars0[PV->RHO] + vars1[PV->RHO]);
    const MFloat T = p / rho * srfcCoeff[SC->W_MEAN] * m_fGasConstant;

    // Set thermodynamic state of the gas
    gas->setState_TPY(T, p, &Y_k[0]);

    srfcCoeff[SC->CP] = gas->cp_mass();

    srfcCoeff[SC->MU] = trans->viscosity();
    srfcCoeff[SC->LAMBDA] = trans->thermalConductivity();

    if(m_multiDiffusion) {
      // fortran ordering for multicomponent coefficients!
      trans->getMultiDiffCoeffs(PV->m_noSpecies, &srfcCoeff[SC->D[0]]);
      trans->getThermalDiffCoeffs(&srfcCoeff[SC->DT[0]]);
    } else {
      trans->getMixDiffCoeffs(&srfcCoeff[SC->D[0]]);
    }
  }
}

/**
 * @brief Computes gamma at the cell. Used for the computation of the time step.
 *
 * @tparam nDim Number of dimensions
 * @param pvarsCell Pointer to the primitive variables at the cell
 * @param avarsCell Pointer to the additional variables at the cell
 * @param gas Pointer to the Cantera::ThermoPhase object
 */
template <MInt nDim>
inline void FvSysEqnDetChem<nDim>::computeGamma(const MFloat* const pvarsCell, MFloat* const avarsCell,
                                                std::shared_ptr<Cantera::ThermoPhase> gas) {
  const MFloat T = pvarsCell[PV->P] / pvarsCell[PV->RHO] * avarsCell[AV->W_MEAN] * m_fGasConstant;

  gas->setState_TPY(T, pvarsCell[PV->P], &pvarsCell[PV->Y[0]]);

  const MFloat Cp = gas->cp_mass();
  const MFloat Cv = Cp - m_gasConstant / avarsCell[AV->W_MEAN];

  avarsCell[AV->GAMMA] = Cp / Cv;
}

/**
 * @brief Reconstructs the conservative slopes from the primitive variables and primitive slopes at the cell. Used in
 * the interpolation of the conservative variables unto the newly created leaf cells during the grid refinement step
 *
 * @tparam nDim Number of dimensions
 * @param pvarsCell Pointer to the primitive variables at the cell
 * @param cvarsCell Pointer to the conservative variables at the cell
 * @param avarsCell Pointer to the additional variables at the cell
 * @param slopesCell Pointer to the primitive slopes at the cell
 * @return std::vector<std::vector<MFloat>> Returns a vector of vectors with dimensions noVariables * nDim
 */
template <MInt nDim>
inline std::vector<std::vector<MFloat>>
FvSysEqnDetChem<nDim>::conservativeSlopes(const MFloat* const pvarsCell, const MFloat* const NotUsed(cvarsCell),
                                          const MFloat* const avarsCell, const MFloat* const slopesCell) {
  MFloat U2 = F0;
  for(MInt i = 0; i < nDim; i++) {
    U2 += POW2(pvarsCell[PV->VV[i]]);
  }

  std::vector<std::vector<MFloat>> dQ(PV->noVariables, std::vector<MFloat>(nDim));
  for(MInt d = 0; d < nDim; d++) {
    MFloat term1 = F0, term2 = F0, term3 = F0;
    const MFloat meanMolarMass = avarsCell[AV->W_MEAN];
    const MFloat T = pvarsCell[PV->P] / pvarsCell[PV->RHO] * meanMolarMass * m_fGasConstant;
    const MBool isLowTempRegion = (T < m_NASA->transitionTemp);

    const MFloat* const NASACoeffs = isLowTempRegion ? ALIGNED_F(m_NASA->lowTemp) : ALIGNED_F(m_NASA->highTemp);
    const MFloat* const NASAIntegralCoeffs =
        isLowTempRegion ? ALIGNED_F(m_NASA->integralLowTemp) : ALIGNED_F(m_NASA->integralHighTemp);
    const MFloat* const NASAIntegralConstants = isLowTempRegion ? ALIGNED_F(m_NASA->lowTempIntegrationConstantsEnergy)
                                                                : ALIGNED_F(m_NASA->highTempIntegrationConstantsEnergy);


    MFloat dT0Term = F0;
    MFloat sensibleEnergy = F0;
    dT0Term += slopesCell[PV->RHO * nDim + d] * m_gasConstant * m_NASA->referenceTemp / meanMolarMass;
    for(MUint s = 0; s < PV->m_noSpecies; s++) {
      dT0Term += pvarsCell[PV->RHO] * m_gasConstant * m_NASA->referenceTemp * slopesCell[PV->Y[s] * nDim + d]
                 * m_species->fMolarMass[s];
    }

    MFloat summ_YIntCvdT = F0;
    MFloat summ_dYIntCvdT = F0;
    for(MUint s = 0; s < PV->m_noSpecies; s++) {
      const MInt offsetIntegralCoeffs = m_NASA->noNASACoefficientsCpPolynomial * s;

      // Calculate NASA polynomial for species s
      MFloat IntCvdT = F0;
      MFloat NASAIntegralPolynomial = F0;

      // Horner's rule for polynomials
      for(MInt i = 4; i >= 0; --i) {
        NASAIntegralPolynomial = NASAIntegralPolynomial * T + NASAIntegralCoeffs[offsetIntegralCoeffs + i];
      }

      NASAIntegralPolynomial *= m_species->specificGasConstant[s];
      NASAIntegralPolynomial = (NASAIntegralPolynomial - m_species->specificGasConstant[s]) * T;
      IntCvdT = NASAIntegralPolynomial - NASAIntegralConstants[s];
      sensibleEnergy += pvarsCell[PV->Y[s]] * IntCvdT;
      summ_YIntCvdT += IntCvdT * pvarsCell[PV->Y[s]];
      summ_dYIntCvdT += IntCvdT * slopesCell[PV->Y[s] * nDim + d];
    }

    term1 = slopesCell[PV->RHO * nDim + d] * summ_YIntCvdT;
    term2 = pvarsCell[PV->RHO] * summ_dYIntCvdT;

    //#####################################
    MFloat summSpeciesGradients = F0;
    for(MUint s = 0; s < PV->m_noSpecies; s++) {
      summSpeciesGradients += m_species->fMolarMass[s] * slopesCell[PV->Y[s] * nDim + d];
    }

    const MFloat dW_dx = -(POW2(avarsCell[AV->W_MEAN])) * summSpeciesGradients;

    const MFloat dT_dx = m_fGasConstant
                         * (pvarsCell[PV->P] / pvarsCell[PV->RHO] * dW_dx
                            + F1 / pvarsCell[PV->RHO] * meanMolarMass
                                  * (slopesCell[PV->P * nDim + d]
                                     - pvarsCell[PV->P] / pvarsCell[PV->RHO] * slopesCell[PV->RHO * nDim + d]));

    // MFloat term3sub1 = F0, term3sub2 = F0;
    MFloat term3sub2 = F0;
    for(MUint s = 0; s < PV->m_noSpecies; s++) {
      // const MInt offsetIntegralCoeffs = m_NASA->noNASACoefficientsCpPolynomial * s;
      const MInt offsetNASAPolynomial = m_NASA->noNASACoefficients * s;
      /*
      MFloat intdCvdT = F0;
      MFloat integrationConstant = F0;
      //Cp = F0;
      MFloat t2 = F0, t1 = F0;
      for(MInt i = 1; i < m_NASA->noNASACoefficientsCpPolynomial; ++i) {
        if (isLowTempRegion) {
        //intdCvdT = intdCvdT * T + NASACoeffs[offsetNASAPolynomial + i];
        //integrationConstant = integrationConstant * m_NASA->referenceTemp + NASACoeffs[offsetNASAPolynomial + i];
        intdCvdT += m_species->specificGasConstant[s] * NASACoeffs[offsetNASAPolynomial + i] * std::pow(T, MFloat(i));
        integrationConstant += m_species->specificGasConstant[s] * NASACoeffs[offsetNASAPolynomial + i] *
      std::pow(m_NASA->referenceTemp, MFloat(i)); } else { intdCvdT += m_species->specificGasConstant[s] *
      NASACoeffs[offsetNASAPolynomial + i] * std::pow(T, MFloat(i)); integrationConstant +=
      m_species->specificGasConstant[s] * NASACoeffs[offsetNASAPolynomial + i] * std::pow(m_NASA->transitionTemp,
      MFloat(i)); t2 += m_species->specificGasConstant[s] * m_NASA->lowTemp[offsetNASAPolynomial + i] *
      std::pow(m_NASA->transitionTemp, MFloat(i)); t1 += m_species->specificGasConstant[s] *
      m_NASA->lowTemp[offsetNASAPolynomial + i] * std::pow(m_NASA->referenceTemp, MFloat(i)); integrationConstant += t1
      - t2;
        }
      }

      term3sub1 += pvarsCell[PV->RHO] * pvarsCell[PV->Y[s]] * dT_dx * (intdCvdT - integrationConstant);*/

      MFloat Cv = F0;
      for(MInt i = 4; i >= 0; --i) {
        Cv = Cv * T + NASACoeffs[offsetNASAPolynomial + i];
      }

      // this computes heat capacity an constant pressure (cp)
      Cv *= m_species->specificGasConstant[s];

      // this corrects the value to heat capacity at constant volume (cv = cp - R_k)
      Cv -= m_species->specificGasConstant[s];

      term3sub2 += pvarsCell[PV->RHO] * pvarsCell[PV->Y[s]] * Cv * dT_dx;
    }

    term3 = term3sub2;

    MFloat sensibleEnergyTerm = term1 + term2 + term3 - dT0Term;

    MFloat velocityTerm = F1B2 * U2 * slopesCell[PV->RHO * nDim + d];

    for(MInt j = 0; j < nDim; j++) {
      velocityTerm += pvarsCell[PV->RHO] * pvarsCell[PV->VV[j]] * slopesCell[PV->VV[j] * nDim + d]; // to change
    }

    dQ[CV->RHO_E][d] = sensibleEnergyTerm + velocityTerm;

    dQ[CV->RHO][d] = slopesCell[PV->RHO * nDim + d];

    for(MInt j = 0; j < nDim; j++) {
      dQ[CV->RHO_VV[j]][d] =
          pvarsCell[PV->VV[j]] * slopesCell[PV->RHO * nDim + d] + pvarsCell[PV->RHO] * slopesCell[PV->VV[j] * nDim + d];
    }
    for(MUint s = 0; s < PV->m_noSpecies; ++s) {
      dQ[CV->RHO_Y[s]][d] =
          pvarsCell[PV->Y[s]] * slopesCell[PV->RHO * nDim + d] + pvarsCell[PV->RHO] * slopesCell[PV->Y[s] * nDim + d];
    }
  }

  return dQ;
}

/**
 * @brief Computes the mean molar weight from the conservative variables at the cell
 *
 * @tparam nDim Number of dimensions
 * @param cvarsCell Pointer to the conservative variables at the cell
 * @param avarsCell Pointer to the additional variables at the cell
 */
template <MInt nDim>
inline void FvSysEqnDetChem<nDim>::computeMeanMolarWeight_CV(const MFloat* const cvarsCell, MFloat* const avarsCell) {
  MFloat fMeanMolarWeight = F0;
  for(MUint s = 0; s < PV->m_noSpecies; ++s) {
    MFloat fRho = F1 / cvarsCell[CV->RHO];
    const MFloat Y_s = fRho * cvarsCell[CV->RHO_Y[s]];
    fMeanMolarWeight += Y_s * m_species->fMolarMass[s];
  }
  avarsCell[AV->W_MEAN] = (F1 / fMeanMolarWeight);
}

/**
 * @brief Computes the mean molar weight from the primitive variables at the cell
 *
 * @tparam nDim Number of dimensions
 * @param pvarsCell Pointer to the primitive variables at the cell
 * @param avarsCell Pointer to the additional variables at the cell
 */
template <MInt nDim>
inline void FvSysEqnDetChem<nDim>::computeMeanMolarWeight_PV(const MFloat* const pvarsCell, MFloat* const avarsCell) {
  MFloat fMeanMolarWeight = F0;
  for(MUint s = 0; s < PV->m_noSpecies; ++s) {
    fMeanMolarWeight += pvarsCell[PV->Y[s]] * m_species->fMolarMass[s];
  }
  avarsCell[AV->W_MEAN] = (F1 / fMeanMolarWeight);
}

template <MInt nDim>
MFloat FvSysEqnDetChem<nDim>::computePhi(const MFloat* const Y) {
  MInt indexOxidizer = m_species->speciesMap[m_oxidizer];
  MInt indexFuel = m_species->speciesMap[m_fuel];

  MFloat stochiometricMolarMassRatio = m_species->molarMass[indexFuel] / m_species->molarMass[indexOxidizer];

  MFloat phi = (Y[indexFuel] / Y[indexOxidizer]) / (m_fuelOxidizerStochiometricRatio * stochiometricMolarMassRatio);

  return phi;
}

#else

// To allow compilation of SysEqnDetChem when CANTERA is not defined
template <MInt nDim>
class FvSysEqnDetChem : public FvSysEqnNS<nDim> {
 public:
  FvSysEqnDetChem(const MInt solverId, const MInt noSpecies);

  void evaluateSensibleEnergy(MFloat&, const MFloat* const, const MFloat&);
  void computeMeanMolarWeight_PV(const MFloat* const, MFloat* const);
  void computeMeanMolarWeight_CV(const MFloat* const, MFloat* const);

  struct AdditionalVariables;
  struct SurfaceCoefficients;
  struct NASACoefficients;
  struct SpeciesProperties;
  SpeciesProperties* m_species;
  NASACoefficients* m_NASA;
  static constexpr MBool hasAV = true;
  static constexpr MBool hasSC = true; // add to other sysEqn
  AdditionalVariables* AV;
  SurfaceCoefficients* SC;
};

template <MInt nDim>
void FvSysEqnDetChem<nDim>::computeMeanMolarWeight_CV(const MFloat* const, MFloat* const) {}

template <MInt nDim>
void FvSysEqnDetChem<nDim>::computeMeanMolarWeight_PV(const MFloat* const, MFloat* const) {}

template <MInt nDim>
void FvSysEqnDetChem<nDim>::evaluateSensibleEnergy(MFloat&, const MFloat* const, const MFloat&) {}

// Stores species properties
template <MInt nDim>
struct FvSysEqnDetChem<nDim>::SpeciesProperties {
  const MString majorSpecies = "N2";
  MInt majorSpeciesIndex;

  const MFloat referenceTemp = 298.0;

  std::map<std::string, MInt> speciesMap;

  std::vector<MString> speciesName;
  MFloat* molarMass = nullptr;
  MFloat* fMolarMass = nullptr;
  MFloat* specificGasConstant = nullptr;
  MFloat* fSpecificGasConstant = nullptr;
  MFloat* standardHeatFormation = nullptr;
};

/// \brief Static indices for accessing surface coefficients
template <MInt nDim>
struct FvSysEqnDetChem<nDim>::SurfaceCoefficients {
  static const MInt Segfault = std::numeric_limits<MInt>::min();

  const MInt m_noDiffusionCoefficients;
  const MInt m_noThermalDiffusionCoefficients;
  const MInt m_noSurfaceCoefficients;

  static constexpr MInt MU = 0;
  static constexpr MInt LAMBDA = 1;
  static constexpr MInt CP = 2;
  static constexpr MInt W_MEAN = 3;
  static constexpr MInt D0 = 4;

  MInt* D = nullptr;
  MInt* DT = nullptr;
};

/// \brief Static indices for accessing additional variables
template <MInt nDim>
struct FvSysEqnDetChem<nDim>::AdditionalVariables {
  static constexpr MInt noVariables = 2;

  static constexpr MInt GAMMA = 0;
  static constexpr MInt W_MEAN = 1;
};

#endif // WITH_CANTERA

#endif // FvSysEqnDetChem_H_
