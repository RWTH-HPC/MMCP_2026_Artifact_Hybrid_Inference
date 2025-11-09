// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef FVSYSEQNNS_H_
#define FVSYSEQNNS_H_

#include <algorithm>
#include <array>
#include <iomanip>
#include <iterator>
#include <limits>
#include <numeric>

#include "INCLUDE/maiaconstants.h"
//#include "fvsyseqn.h"
#include "INCLUDE/maiatypes.h"
#include "MEMORY/alloc.h"
#include "enums.h"
#include "filter.h"

// Only needed (really: possible to use) in single-threaded applications
#ifndef _OPENMP
#include "MEMORY/scratch.h"
#endif


template <MInt nDim>
class FvSysEqnNS {
  // Declare parent a friend so that CRTP can access private methods/members
  // friend class FvSysEqn<nDim>;

 public:
  FvSysEqnNS(const MInt solverId, const MInt noSpecies);

  template <MInt scheme = AUSM>
  inline void Ausm(const MInt orientation, const MFloat upwindCoefficient, const MFloat A, const MFloat* const leftVars,
                   const MFloat* const rightVars, const MFloat* const srfcCoeff, MFloat* const flux) {
    IF_CONSTEXPR(scheme == AUSM) { Ausm_(orientation, upwindCoefficient, A, leftVars, rightVars, srfcCoeff, flux); }
    else IF_CONSTEXPR(scheme == AUSMPLUS) {
      AusmPlus_(orientation, upwindCoefficient, A, leftVars, rightVars, srfcCoeff, flux);
    }
    else IF_CONSTEXPR(scheme == SLAU) {
      Slau_(orientation, upwindCoefficient, A, leftVars, rightVars, srfcCoeff, flux);
    }
  }

  inline void Ausm_(const MInt orientation, const MFloat upwindCoefficient, const MFloat A,
                    const MFloat* const leftVars, const MFloat* const rightVars, const MFloat* const NotUsed(srfcCoeff),
                    MFloat* const flux);
  inline void AusmPlus_(const MInt orientation, const MFloat upwindCoefficient, const MFloat A,
                        const MFloat* const leftVars, const MFloat* const rightVars,
                        const MFloat* const NotUsed(srfcCoeff), MFloat* const flux);
  inline void Slau_(const MInt orientation, const MFloat NotUsed(upwindCoefficient), const MFloat A,
                    const MFloat* const leftVars, const MFloat* const rightVars, const MFloat* const NotUsed(srfcCoeff),
                    MFloat* const flux);

  inline void AusmBndryCorrection(const MInt orientation, const MFloat A, const MFloat* const leftVars,
                                  const MFloat* const rightVars, MFloat* const flux);

  inline void AusmALECorrection(const MInt orientation, const MFloat A, MFloat* const flux, MFloat* const surfVars,
                                const MFloat* const bndrySurfVars);

  template <MInt centralizeScheme>
  inline void centralizeSurfaceVariables(MFloat* const varL, MFloat* const varR,
                                         [[maybe_unused]] const MInt orientation,
                                         [[maybe_unused]] const MFloat levelFac);

  // Dispatch function for classic five-point-style viscous flux schemes
  template <MInt stencil>
  inline void viscousFlux(const MInt orientation, const MFloat A, const MFloat* const vars0, const MFloat* const vars1,
                          const MFloat* const slope0, const MFloat* const slope1, const MFloat* const srfcCoeff,
                          const MFloat f0, const MFloat f1, MFloat* const flux) {
    IF_CONSTEXPR(stencil == FIVE_POINT) {
      viscousFluxFivePoint(orientation, A, vars0, vars1, slope0, slope1, srfcCoeff, f0, f1, flux);
    }
    else {
      TERMM(1, "Viscous Flux Scheme not implemented.");
    }
  }

  // Dispatch function for compact viscous flux schemes
  template <MInt stencil>
  inline void viscousFlux(const MInt orientation, const MFloat A, const MBool isBndry,
                          const MFloat* const surfaceCoords, const MFloat* const coord0, const MFloat* const coord1,
                          const MFloat* const cellVars0, const MFloat* const cellVars1, const MFloat* const vars0,
                          const MFloat* const vars1, const MFloat* const slope0, const MFloat* const slope1,
                          const MFloat f0, const MFloat f1, MFloat* const flux) {
    IF_CONSTEXPR(stencil == THREE_POINT) {
      viscousFluxThreePoint(orientation, A, isBndry, surfaceCoords, coord0, coord1, cellVars0, cellVars1, vars0, vars1,
                            slope0, slope1, f0, f1, flux);
    }
    else IF_CONSTEXPR(stencil == FIVE_POINT_STABILIZED) {
      viscousFluxStabilized(orientation, A, isBndry, surfaceCoords, coord0, coord1, cellVars0, cellVars1, vars0, vars1,
                            slope0, slope1, f0, f1, flux);
    }
    else {
      TERMM(1, "Viscous Flux Scheme not implemented.");
    }
  }

  inline void viscousFluxFivePoint(const MInt orientation, const MFloat A, const MFloat* const vars0,
                                   const MFloat* const vars1, const MFloat* const slope0, const MFloat* const slope1,
                                   const MFloat* const NotUsed(srfcCoeff), const MFloat f0, const MFloat f1,
                                   MFloat* const flux);

  inline void viscousFluxThreePoint(const MInt orientation, const MFloat A, const MBool isBndry,
                                    const MFloat* const surfaceCoords, const MFloat* const coord0,
                                    const MFloat* const coord1, const MFloat* const cellVars0,
                                    const MFloat* const cellVars1, const MFloat* const vars0, const MFloat* const vars1,
                                    const MFloat* const slope0, const MFloat* const slope1, const MFloat f0,
                                    const MFloat f1, MFloat* const flux);

  inline void viscousFluxStabilized(const MInt orientation, const MFloat A, const MBool isBndry,
                                    const MFloat* const surfaceCoords, const MFloat* const coord0,
                                    const MFloat* const coord1, const MFloat* const cellVars0,
                                    const MFloat* const cellVars1, const MFloat* const vars0, const MFloat* const vars1,
                                    const MFloat* const slope0, const MFloat* const slope1, const MFloat f0,
                                    const MFloat f1, MFloat* const flux);

  // Dispatch function for wall-model viscous flux correction with five-point stencil
  template <MInt stencil>
  inline void wmViscousFluxCorrection(const MInt orientation, const MFloat A, const MFloat* const vars0,
                                      const MFloat* const vars1, const MFloat* const slope0, const MFloat* const slope1,
                                      const MFloat f0, const MFloat f1, MFloat* const flux, MFloat const mue_wm) {
    IF_CONSTEXPR(stencil == FIVE_POINT) {
      wmViscousFluxCorrectionFivePoint(orientation, A, vars0, vars1, slope0, slope1, f0, f1, flux, mue_wm);
    }
    else {
      TERMM(1, "Viscous Flux Scheme not implemented.");
    }
  }

  inline void wmViscousFluxCorrectionFivePoint(const MInt orientation, const MFloat A, const MFloat* const vars0,
                                               const MFloat* const vars1, const MFloat* const slope0,
                                               const MFloat* const slope1, const MFloat f0, const MFloat f1,
                                               MFloat* const flux, MFloat const mue_wm);


  inline void computePrimitiveVariables(const MFloat* const cvarsCell, MFloat* const pvarsCell,
                                        const MFloat* const NotUsed(avarsCell));

  inline void computeConservativeVariables(const MFloat* const pvarsCell, MFloat* const cvarsCell,
                                           const MFloat* const NotUsed(avarsCell));

  inline std::vector<std::vector<MFloat>> conservativeSlopes(const MFloat* const pvarsCell,
                                                             const MFloat* const cvarsCell,
                                                             const MFloat* const avarsCell,
                                                             const MFloat* const slopesCell);

  void computeVolumeForces(const MInt, MFloat*, MFloat*, const MFloat, const MFloat* const, const MInt,
                           const MInt* const, const MFloat* const, const MInt, MFloat){};

  /// Speed of sound: a = sqrt(gamma * pressure / density)
  inline constexpr MFloat speedOfSound(const MFloat density, const MFloat pressure) {
    return std::sqrt(m_gamma * pressure / density);
  };

  /// speed of sound squared a^2 = gamma * pressure / density
  inline constexpr MFloat speedOfSoundSquared(const MFloat density, const MFloat pressure) {
    return m_gamma * pressure / density;
  };

  /// Speed of sound: a = sqrt(T)
  inline constexpr MFloat speedOfSound(const MFloat temperature) { return std::sqrt(temperature); };

  /// Temperature: T = gamma * pressure / density (equation of state - ideal gas law)
  inline constexpr MFloat temperature_ES(const MFloat density, const MFloat pressure) {
    return m_gamma * pressure / density;
  };
  /// pressure: p = rho * T / gamma (equation of state - ideal gas law)
  inline constexpr MFloat pressure_ES(const MFloat temperture, const MFloat density) {
    return density * temperture / m_gamma;
  }

  /// density: rho = gamma * p / T (equation of state - ideal gas law)
  inline constexpr MFloat density_ES(const MFloat pressure, const MFloat temperature) {
    return m_gamma * pressure / temperature;
  }

  /// Temperature: T = 1 / (1 + (gamma - 1)/2 * Ma^2) (isentropic relationship)
  inline constexpr MFloat temperature_IR(const MFloat Ma) { return 1 / (1.0 + 0.5 * m_gammaMinusOne * POW2(Ma)); }

  /// pressure: p = ( T ^(gamma/(gamma -1 ))) / gamma (isentropic relationship)
  /// NOTE: this is under consideration of the selected pressure non-dimensionalization!
  inline constexpr MFloat pressure_IR(const MFloat temperature) {
    return pow(temperature, m_FGammaBGammaMinusOne) / m_gamma;
  }

  /// pressure/density iteration based on mass-flux and p_old
  /// NOTE: see thesis of Ingolf Hoerschler
  inline constexpr MFloat pressure_IRit(const MFloat pressure, const MFloat massFlux) {
    return pow(1.0 - m_gammaMinusOne * F1B2 * pow(pressure, -2.0 * m_F1BGamma) * POW2(massFlux),
               m_FGammaBGammaMinusOne);
  }

  /// density: rho = T^(1/(gamma -1 )) (isentropic relationship)
  inline constexpr MFloat density_IR(const MFloat temperature) { return pow(temperature, m_F1BGammaMinusOne); }

  /// density: rho = (p * gamma)^(1/gamma) (isentropic relationship)
  inline constexpr MFloat density_IR_P(const MFloat pressure) { return pow(pressure * m_gamma, m_F1BGamma); }

  /// pressure from conservative variables: cv = (rho, ||rho_u||^2, rhoE)
  /// p = (gamma - 1) * (rhoE - 0.5 ||rho_u||^2 / rho )
  inline constexpr MFloat pressure(const MFloat density, const MFloat momentumDensitySquared,
                                   const MFloat energyDensity) {
    return m_gammaMinusOne * (energyDensity - 0.5 / density * momentumDensitySquared);
  };

  /// energy from primitive variables: pv = (p, rho, u)
  /// iE = p / (gamma -1 ) + 0.5 * rho * ||u||^2
  inline constexpr MFloat internalEnergy(const MFloat pressure, const MFloat density, const MFloat velocitySquared) {
    return pressure * m_F1BGammaMinusOne + 0.5 * density * velocitySquared;
  }

  /// energy from primitive pressure variable (or delta p)
  /// pE = p / (gamma -1 )
  inline constexpr MFloat pressureEnergy(const MFloat pressure) { return pressure * m_F1BGammaMinusOne; }

  /// enthalpy from primitive variables
  inline constexpr MFloat enthalpy(const MFloat pressure, const MFloat density) {
    return pressure / density * m_FGammaBGammaMinusOne;
  }

  /// entropy from primitive variables
  inline constexpr MFloat entropy(const MFloat pressure, const MFloat density) {
    return pressure / pow(density, m_gamma);
  }

  /// Crocco-Busemann relation
  inline constexpr MFloat CroccoBusemann(const MFloat Ma, const MFloat x) {
    return (1 + 0.5 * m_gammaMinusOne * POW2(Ma) * x * (1.0 - x));
  }

  /// van-Driest Transformation (correspods to R*H)
  inline constexpr MFloat vanDriest(const MFloat Ma) {
    return (sqrt((m_gammaMinusOne / 2.0 * POW2(Ma) * pow(m_Pr, F1B3))
                 / (1.0 + m_gammaMinusOne / 2.0 * POW2(Ma) * pow(m_Pr, F1B3))));
  }

  /// Sutherland-law of viscosity
  /// mue = T^3/2 * (1+S/T_0)(T + S/T_0)
  /// with the default values S=110.4 k and T_0 = 273.15 K
  inline constexpr MFloat sutherlandLaw(const MFloat T) {
    return (T * sqrt(T) * m_sutherlandPlusOne) / (T + m_sutherlandConstant);
  }

  /// Computes the time-step from the CFL condition for the Euler-eqts as follows:
  ///
  /// dt =   C * dx / (||u|| + a)
  inline MFloat computeTimeStepEulerMagnitude(const MFloat rho, const std::array<MFloat, nDim> u, const MFloat p,
                                              const MFloat C, const MFloat dx) {
    const MFloat a = speedOfSound(rho, p);
    const MFloat u_mag = std::sqrt(std::inner_product(&u[0], &u[nDim], &u[0], 0.0));
    return C * dx / (u_mag + a);
  };

  /// reference heat capacity ratio used in the non-dimensioanlization
  /// (i.e. at the reference condition with reference mixture/state/material)
  inline constexpr MFloat gamma_Ref() { return m_gamma; }

  /// reference heat capacity at const. pressure as used in the non-dimensioanlization
  /// (i.e. at the reference condition with reference mixture/state/material)
  inline constexpr MFloat cp_Ref() { return m_F1BGammaMinusOne; }

  /// reference heat capacity at const. volume as used in the non-dimensioanlization
  /// (i.e. at the reference condition with reference mixture/state/material)
  inline constexpr MFloat cv_Ref() { return m_F1BGamma * m_F1BGammaMinusOne; }

  /// reference pressure as used in the non-dimensioanlization
  /// (i.e. at the reference condition with reference mixture/state/material)
  inline constexpr MFloat p_Ref() { return m_F1BGamma; }

  /// Computes the time-step from the stability condition for the Diffusion Eqt.
  ///
  /// dt = C * dx^2 / diffusion_coefficient  where C = (0, 1/4]
  inline constexpr MFloat computeTimeStepDiffusion(const MFloat diffusion_coefficient, const MFloat C,
                                                   const MFloat dx) {
    return C * POW2(dx) / diffusion_coefficient;
  };

 public:
  static constexpr MInt m_noRansEquations = 0;
  static constexpr MInt m_ransModel = NORANS;

 protected:
  MUint m_noSpecies;

  MInt m_solverId;

  // Equation specific values
  MFloat m_gamma = 1.4;
  MFloat m_gammaMinusOne = m_gamma - 1;
  MFloat m_F1BGammaMinusOne = 1 / m_gammaMinusOne;
  MFloat m_FGammaBGammaMinusOne = m_gamma / m_gammaMinusOne;
  MFloat m_F1BGamma = 1 / m_gamma;
  MFloat m_Pr = 0.72;
  MFloat m_F1BPr = 1 / m_Pr;
  MFloat m_referenceTemperature = {};
  MFloat m_sutherlandPlusOne = {};
  MFloat m_sutherlandConstant = {};
  MFloat m_sutherlandPlusOneThermal = {};
  MFloat m_sutherlandConstantThermal = {};
  MFloat m_gFGMOrPr = {};

  // Stencil specific properties
  MFloat m_enhanceThreePointViscFluxFactor{};

  MFloat m_centralizeSurfaceVariablesFactor = 0.0;

  void readProperties();
  static inline MFloat sgn(MFloat val) { return (val < F0) ? -F1 : F1; }

  static constexpr std::array<MInt, nDim> getArray012() {
    IF_CONSTEXPR(nDim == 2) {
      std::array<MInt, 2> a = {0, 1};
      return a;
    }
    else {
      std::array<MInt, 3> a = {0, 1, 2};
      return a;
    }
  }

 public:
  // Equation specific values
  MFloat m_Re0 = {};
  MFloat m_muInfinity = {};

  // Hold indices for primitive and conservative variables
  struct ConservativeVariables;
  struct FluxVariables; // these are the variables for which the fluxes are calculated
                        // usually the same as conservative variables
  struct PrimitiveVariables;
  struct AdditionalVariables; // these are additional variables that are saved in the cell collector
                              // they might be used in SysEqn functions etc.
  struct SurfaceCoefficients;
  static constexpr MBool hasAV = false;
  static constexpr MBool hasSC = false;
  ConservativeVariables* CV = nullptr;
  FluxVariables* FV = nullptr;
  PrimitiveVariables* PV = nullptr;
  AdditionalVariables* AV = nullptr;
  SurfaceCoefficients* SC = nullptr;

  // Index helper constants for viscousFlux
  // template <MInt nDim>
  static const MUint index0[nDim];

  // template <MInt nDim>
  static const MUint index1[nDim];
};


/// \brief Static indices for accessing conservative variables
/// in nDim spatial dimensions
template <MInt nDim>
struct FvSysEqnNS<nDim>::ConservativeVariables {
  static const MInt Segfault = std::numeric_limits<MInt>::min();

  static constexpr MInt RHO_U = 0;
  static constexpr MInt RHO_V = 1;
  static constexpr MInt RHO_W = nDim == 3 ? 2 : Segfault;
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

/// \brief Static indices for accessing flux variables
/// in this SysEqn identical to the conservative variables
template <MInt nDim>
struct FvSysEqnNS<nDim>::FluxVariables : ConservativeVariables {
  FluxVariables(const MInt noSpecies);
};

/// \brief Static indices for accessing primitive variables
/// in nDim spatial dimensions
template <MInt nDim>
struct FvSysEqnNS<nDim>::PrimitiveVariables {
  static const MInt Segfault = std::numeric_limits<MInt>::min();

  static constexpr MInt U = 0;
  static constexpr MInt V = 1;
  static constexpr MInt W = nDim == 3 ? 2 : Segfault;
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
  void getPrimitiveVariableNames(MString* names);
  ~PrimitiveVariables();
};

/// \brief No additional variables are used in this SysEqn
template <MInt nDim>
struct FvSysEqnNS<nDim>::AdditionalVariables {
  static constexpr MInt noVariables = 0;
};

template <MInt nDim>
struct FvSysEqnNS<nDim>::SurfaceCoefficients {
  const MInt m_noSurfaceCoefficients = 0;
};

template <MInt nDim>
inline void FvSysEqnNS<nDim>::Ausm_(const MInt orientation, const MFloat upwindCoefficient, const MFloat A,
                                    const MFloat* const leftVars, const MFloat* const rightVars,
                                    const MFloat* const NotUsed(srfcCoeff), MFloat* const flux) {
  // catch the primitive variables rho and p,
  // compute speed of sound, and interface mach number
  const MFloat RHOL = leftVars[PV->RHO];
  const MFloat PL = leftVars[PV->P];
  const MFloat AL = sqrt(m_gamma * mMax(MFloatEps, PL / mMax(MFloatEps, RHOL)));
  const MFloat ML = leftVars[orientation] / AL;

  const MFloat RHOR = rightVars[PV->RHO];
  const MFloat PR = rightVars[PV->P];
  const MFloat AR = sqrt(m_gamma * mMax(MFloatEps, PR / mMax(MFloatEps, RHOR)));
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

  // calculation of the energy:
  const MFloat PLfRHOL = PL / RHOL;
  const MFloat PRfRHOR = PR / RHOR;

  // velocity magnitudes
  const MFloat U2L = std::inner_product(&leftVars[PV->U], &leftVars[PV->U] + nDim, &leftVars[PV->U], 0.0);
  const MFloat U2R = std::inner_product(&rightVars[PV->U], &rightVars[PV->U] + nDim, &rightVars[PV->U], 0.0);

  const MFloat e0 = PLfRHOL * m_F1BGammaMinusOne + 0.5 * U2L + PLfRHOL;
  const MFloat e1 = PRfRHOR * m_F1BGammaMinusOne + 0.5 * U2R + PRfRHOR;

  std::array<MFloat, nDim> pFactor{};
  pFactor[orientation] = 1.0;

  for(MUint n = 0; n < nDim; n++) {
    flux[FV->RHO_VV[n]] = (RHO_U2 * (leftVars[PV->VV[n]] + rightVars[PV->VV[n]])
                           + AbsRHO_U2 * (leftVars[PV->VV[n]] - rightVars[PV->VV[n]]) + PLR * pFactor[n])
                          * A;
  }

  flux[FV->RHO_E] = (RHO_U2 * (e0 + e1) + AbsRHO_U2 * (e0 - e1)) * A;
  flux[FV->RHO] = 2.0 * RHO_U2 * A;

  // Flux calculation for species transport
  // TODO labels:FV @Julian, Make noSpecies constexpr
  for(MUint s = 0; s < PV->m_noSpecies; ++s) {
    const MFloat YsL = leftVars[PV->Y[s]];
    const MFloat YsR = rightVars[PV->Y[s]];
    flux[FV->RHO_Y[s]] = (RHO_U2 * (YsL + YsR) + AbsRHO_U2 * (YsL - YsR)) * A;
  }
}


template <MInt nDim>
inline void FvSysEqnNS<nDim>::AusmPlus_(const MInt orientation, const MFloat upwindCoefficient, const MFloat A,
                                        const MFloat* const leftVars, const MFloat* const rightVars,
                                        const MFloat* const NotUsed(srfcCoeff), MFloat* const flux) {
  if(PV->m_noSpecies > 0) {
    mTerm(1, "Not yet implemented for multiple-species!");
  }

  std::array<MFloat, nDim> pFactor = {};
  pFactor[orientation] = F1;

  MFloat UL = leftVars[PV->U];
  MFloat VL = leftVars[PV->V];
  MFloat WL = 0.0;
  IF_CONSTEXPR(nDim == 3) { WL = leftVars[PV->W]; }
  MFloat RL = leftVars[PV->RHO];
  MFloat PL = leftVars[PV->P];

  MFloat UR = rightVars[PV->U];
  MFloat VR = rightVars[PV->V];
  MFloat WR = 0.0;
  IF_CONSTEXPR(nDim == 3) { WR = rightVars[PV->W]; }
  MFloat RR = rightVars[PV->RHO];
  MFloat PR = rightVars[PV->P];

  MFloat ql = leftVars[orientation];
  MFloat AL = sqrt(m_gamma * PL / RL);
  MFloat KL = 0.0;
  KL = F1B2 * (UL * UL + VL * VL);
  IF_CONSTEXPR(nDim == 3) { KL = F1B2 * (UL * UL + VL * VL + WL * WL); }
  MFloat HL = (m_gamma * PL) / (m_gammaMinusOne * RL) + KL;

  MFloat qr = rightVars[orientation];
  MFloat AR = sqrt(m_gamma * PR / RR);
  MFloat KR = 0.0;
  KR = F1B2 * (UR * UR + VR * VR);
  IF_CONSTEXPR(nDim == 3) { KR = F1B2 * (UR * UR + VR * VR + WR * WR); }
  MFloat HR = (m_gamma * PR) / (m_gammaMinusOne * RR) + KR;

  MFloat ML = ql / AL;
  MFloat MR = qr / AR;

  MFloat MLP = (fabs(ML) > F1) ? F1B2 * (ML + fabs(ML)) : F1B4 * POW2(ML + F1);
  MFloat MRM = (fabs(MR) > F1) ? F1B2 * (MR - fabs(MR)) : -F1B4 * POW2(MR - F1);
  MFloat M = MLP + MRM;
  MFloat P = F1B2 * (PL + PR) + upwindCoefficient * (ML * PL - MR * PR);

  flux[FV->RHO] = F1B2 * (M * (RL * AL + RR * AR) + fabs(M) * (RL * AL - RR * AR)) * A;
  flux[FV->RHO_U] =
      (F1B2 * (M * (RL * AL * UL + RR * AR * UR) + fabs(M) * (RL * AL * UL - RR * AR * UR)) + P * pFactor[0]) * A;
  flux[FV->RHO_V] =
      (F1B2 * (M * (RL * AL * VL + RR * AR * VR) + fabs(M) * (RL * AL * VL - RR * AR * VR)) + P * pFactor[1]) * A;
  IF_CONSTEXPR(nDim == 3) {
    flux[FV->RHO_W] =
        (F1B2 * (M * (RL * AL * WL + RR * AR * WR) + fabs(M) * (RL * AL * WL - RR * AR * WR)) + P * pFactor[2]) * A;
  }
  flux[FV->RHO_E] = F1B2 * (M * (RL * AL * HL + RR * AR * HR) + fabs(M) * (RL * AL * HL - RR * AR * HR)) * A;
}


template <MInt nDim>
inline void FvSysEqnNS<nDim>::Slau_(const MInt orientation, const MFloat NotUsed(upwindCoefficient), const MFloat A,
                                    const MFloat* const leftVars, const MFloat* const rightVars,
                                    const MFloat* const NotUsed(srfcCoeff), MFloat* const flux) {
  if(PV->m_noSpecies > 0) {
    mTerm(1, "Not yet implemented for multiple-species!");
  }

  std::array<MFloat, nDim> pFactor = {};
  pFactor[orientation] = F1;

  const MFloat UL = leftVars[PV->U];
  const MFloat VL = leftVars[PV->V];
  MFloat WL = 0.0;
  IF_CONSTEXPR(nDim == 3) { WL = leftVars[PV->W]; }
  const MFloat RL = leftVars[PV->RHO];
  const MFloat PL = leftVars[PV->P];

  const MFloat UR = rightVars[PV->U];
  const MFloat VR = rightVars[PV->V];
  MFloat WR = 0.0;
  IF_CONSTEXPR(nDim == 3) { WR = rightVars[PV->W]; }
  const MFloat RR = rightVars[PV->RHO];
  const MFloat PR = rightVars[PV->P];

  const MFloat ql = leftVars[orientation];
  const MFloat AL = sqrt(m_gamma * PL / RL);
  MFloat KL = 0.0;
  IF_CONSTEXPR(nDim == 2) { KL = F1B2 * (UL * UL + VL * VL); }
  IF_CONSTEXPR(nDim == 3) { KL = F1B2 * (UL * UL + VL * VL + WL * WL); }

  const MFloat HL = (m_gamma * PL) / (m_gammaMinusOne * RL) + KL;

  const MFloat qr = rightVars[orientation];
  const MFloat AR = sqrt(m_gamma * PR / RR);
  MFloat KR = 0.0;
  IF_CONSTEXPR(nDim == 2) { KR = F1B2 * (UR * UR + VR * VR); }
  IF_CONSTEXPR(nDim == 3) { KR = F1B2 * (UR * UR + VR * VR + WR * WR); }

  const MFloat HR = (m_gamma * PR) / (m_gammaMinusOne * RR) + KR;

  const MFloat ALR = F1B2 * (AL + AR);
  const MFloat ML = ql / ALR;
  const MFloat MR = qr / ALR;

  const MFloat BL = (fabs(ML) > F1) ? F1B2 * (F1 + sgn(ML)) : F1B4 * POW2(ML + F1) * (F2 - ML);
  const MFloat BR = (fabs(MR) > F1) ? F1B2 * (F1 - sgn(MR)) : F1B4 * POW2(MR - F1) * (F2 + MR);

  const MFloat MH = mMin(F1, sqrt(KL + KR) / ALR);
  const MFloat XI = POW2(F1 - MH);
  // const MFloat P = F1B2*( PL+PR + (BL-BR)*(PL-PR) + (F1-XI)*(BL+BR-F1)*(PL+PR) ); //SLAU
  const MFloat P = F1B2 * (PL + PR + (BL - BR) * (PL - PR) + (BL + BR - F1) * sqrt(KL + KR) * (RL + RR) * ALR); // SLAU2
  const MFloat g = -mMax(mMin(ML, F0), -F1) * mMin(mMax(MR, F0), F1);
  const MFloat Q = (RL * fabs(ql) + RR * fabs(qr)) / (RL + RR);
  const MFloat QL = (F1 - g) * Q + g * fabs(ql);
  const MFloat QR = (F1 - g) * Q + g * fabs(qr);
  const MFloat M = F1B2 * (RL * (ql + fabs(QL)) + RR * (qr - fabs(QR)) - XI * (PR - PL) / ALR);

  flux[FV->RHO] = M * A;
  flux[FV->RHO_U] = (F1B2 * ((M + fabs(M)) * UL + (M - fabs(M)) * UR) + P * pFactor[0]) * A;
  flux[FV->RHO_V] = (F1B2 * ((M + fabs(M)) * VL + (M - fabs(M)) * VR) + P * pFactor[1]) * A;
  IF_CONSTEXPR(nDim == 3) { flux[FV->RHO_W] = (F1B2 * ((M + fabs(M)) * WL + (M - fabs(M)) * WR) + P * pFactor[2]) * A; }
  flux[FV->RHO_E] = F1B2 * ((M + fabs(M)) * HL + (M - fabs(M)) * HR) * A;
}


template <MInt nDim>
inline void FvSysEqnNS<nDim>::AusmBndryCorrection(const MInt orientation, const MFloat A, const MFloat* const leftVars,
                                                  const MFloat* const rightVars, MFloat* const flux) {
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
inline void FvSysEqnNS<nDim>::AusmALECorrection(const MInt orientation,
                                                const MFloat A,
                                                MFloat* const flux,
                                                MFloat* const surfVars,
                                                const MFloat* const bndrySurfVars) {
  surfVars[PV->VV[orientation]] = bndrySurfVars[PV->VV[orientation]];
  surfVars[PV->noVariables + PV->VV[orientation]] = bndrySurfVars[PV->VV[orientation]];

  for(MInt v = 0; v < FV->noVariables; ++v) {
    flux[v] = 0.0;
  }

  flux[CV->RHO_VV[orientation]] = bndrySurfVars[PV->P] * A;
  flux[CV->RHO_E] = bndrySurfVars[PV->P] * bndrySurfVars[PV->VV[orientation]] * A;
}

template <MInt nDim>
template <MInt centralizeScheme>
inline void FvSysEqnNS<nDim>::centralizeSurfaceVariables(MFloat* const varL,
                                                         MFloat* const varR,
                                                         [[maybe_unused]] const MInt orientation,
                                                         [[maybe_unused]] const MFloat levelFac) {
  static_assert(centralizeScheme > 0 && centralizeScheme < 6, "Unintended/unimplemented function-call!");

  IF_CONSTEXPR(centralizeScheme == 1) {
    const MFloat z = mMin(F1, mMax(fabs(varL[orientation]) / sqrt(m_gamma * varL[PV->P] / varL[PV->RHO]),
                                   fabs(varR[orientation]) / sqrt(m_gamma * varR[PV->P] / varR[PV->RHO])));
    for(MInt k = 0; k < nDim; k++) {
      const MInt velId = PV->VV[k];
      const MFloat ML = varL[velId];
      const MFloat MR = varR[velId];
      varL[velId] = F1B2 * ((F1 + z) * ML + (F1 - z) * MR);
      varR[velId] = F1B2 * ((F1 + z) * MR + (F1 - z) * ML);
    }
  }
  else IF_CONSTEXPR(centralizeScheme == 2) {
    const MFloat z = mMin(F1, (fabs(varL[orientation]) / sqrt(m_gamma * varL[PV->P] / varL[PV->RHO]))
                                  * (fabs(varR[orientation]) / sqrt(m_gamma * varR[PV->P] / varR[PV->RHO])));
    for(MInt k = 0; k < nDim; k++) {
      const MInt velId = PV->VV[k];
      const MFloat ML = varL[velId];
      const MFloat MR = varR[velId];
      varL[velId] = F1B2 * ((F1 + z) * ML + (F1 - z) * MR);
      varR[velId] = F1B2 * ((F1 + z) * MR + (F1 - z) * ML);
    }
  }
  else IF_CONSTEXPR(centralizeScheme == 3) {
    for(MInt v = 0; v < PV->noVariables; v++) {
      const MFloat Vmean = F1B2 * (varL[v] + varR[v]);
      varL[v] = Vmean;
      varR[v] = Vmean;
    }
  }
  else IF_CONSTEXPR(centralizeScheme == 4) {
    const MFloat z = pow(mMin(F1, (fabs(varL[orientation]) / sqrt(m_gamma * varL[PV->P] / varL[PV->RHO]))
                                      * (fabs(varR[orientation]) / sqrt(m_gamma * varR[PV->P] / varR[PV->RHO]))),
                         levelFac);
    for(MInt k = 0; k < nDim; k++) {
      const MInt velId = PV->VV[k];
      const MFloat ML = varL[velId];
      const MFloat MR = varR[velId];
      varL[velId] = F1B2 * ((F1 + z) * ML + (F1 - z) * MR);
      varR[velId] = F1B2 * ((F1 + z) * MR + (F1 - z) * ML);
    }
  }
  else {
    const MFloat z = mMin(F1, m_centralizeSurfaceVariablesFactor
                                  * mMax((fabs(varL[orientation]) / sqrt(m_gamma * varL[PV->P] / varL[PV->RHO])),
                                         (fabs(varR[orientation]) / sqrt(m_gamma * varR[PV->P] / varR[PV->RHO]))));
    for(MInt k = 0; k < nDim; k++) {
      const MInt velId = PV->VV[k];
      const MFloat ML = varL[velId];
      const MFloat MR = varR[velId];
      varL[velId] = F1B2 * ((F1 + z) * ML + (F1 - z) * MR);
      varR[velId] = F1B2 * ((F1 + z) * MR + (F1 - z) * ML);
    }
  }
}

template <MInt nDim>
inline void FvSysEqnNS<nDim>::viscousFluxFivePoint(const MInt orientation, const MFloat A, const MFloat* const vars0,
                                                   const MFloat* const vars1, const MFloat* const slope0,
                                                   const MFloat* const slope1, const MFloat* const, const MFloat f0,
                                                   const MFloat f1, MFloat* const flux) {
  static const MFloat rhoUDth = m_muInfinity * m_F1BPr;
  static const MFloat F1BRe0 = F1 / m_Re0;

  // calculate the primitve variables on the surface u,v,rho,p
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

  // Temperature on the surface T = gamma * p / rho
  const MFloat T = m_gamma * p * Frho;

  // Indices for the orientations
  const MUint id0 = orientation;
  const MUint id1 = index0[orientation];
  const MUint id2 = index1[orientation];

  // Compute A / Re
  const MFloat dAOverRe = A * F1BRe0;

  // calculate the viscosity with the sutherland law mue = T^3/2 * (1+S/T_0)(T + S/T_0)
  const MFloat mue = (T * sqrt(T) * m_sutherlandPlusOne) / (T + m_sutherlandConstant);

  // calculate the heat flux (T_x = gamma*(p_x/rho - p/rho^2 * rho_x))

  const MFloat lambda = (T * sqrt(T) * m_sutherlandPlusOneThermal) / (T + m_sutherlandConstantThermal);

  const MUint sq0 = PV->P * nDim + id0;
  const MUint sq1 = PV->RHO * nDim + id0;
  const MFloat q = (lambda)*m_gFGMOrPr * Frho
                   * ((f0 * slope0[sq0] + f1 * slope1[sq0]) - p * Frho * (f0 * slope0[sq1] + f1 * slope1[sq1]));

  // Compute the stress terms
  const MUint s00 = id0 * nDim + id0;
  const MUint s01 = id0 * nDim + id1;
  const MUint s10 = id1 * nDim + id0;
  const MUint s11 = id1 * nDim + id1;

  std::array<MFloat, nDim> tau{};
  IF_CONSTEXPR(nDim == 2) {
    tau[id0] = mue * (f0 * (F4B3 * slope0[s00] - F2B3 * slope0[s11]) + f1 * (F4B3 * slope1[s00] - F2B3 * slope1[s11]));
    tau[id1] = mue * (f0 * (slope0[s01] + slope0[s10]) + f1 * (slope1[s01] + slope1[s10]));
  }
  else IF_CONSTEXPR(nDim == 3) {
    const MUint s22 = id2 * nDim + id2;
    const MUint s02 = id0 * nDim + id2;
    const MUint s20 = id2 * nDim + id0;
    tau[id0] = mue
               * (f0 * (F4B3 * slope0[s00] - F2B3 * (slope0[s11] + slope0[s22]))
                  + f1 * (F4B3 * slope1[s00] - F2B3 * (slope1[s11] + slope1[s22])));
    tau[id1] = mue * (f0 * (slope0[s01] + slope0[s10]) + f1 * (slope1[s01] + slope1[s10]));
    tau[id2] = mue * (f0 * (slope0[s02] + slope0[s20]) + f1 * (slope1[s02] + slope1[s20]));
  }

  // Compute the flux
  for(MUint n = 0; n < nDim; n++) {
    flux[FV->RHO_VV[n]] -= dAOverRe * tau[n];
  }
  flux[FV->RHO_E] -= dAOverRe * (std::inner_product(velocity.begin(), velocity.end(), tau.begin(), 0.0) + q);

  // progress variable
  const MFloat c = dAOverRe * rhoUDth;
  for(MUint s = 0; s < PV->m_noSpecies; ++s) {
    const MUint i = PV->Y[s] * nDim + id0;
    flux[FV->RHO_Y[s]] -= c * (f0 * slope0[i] + f1 * slope1[i]);
  }
}

template <MInt nDim>
inline void FvSysEqnNS<nDim>::wmViscousFluxCorrectionFivePoint(const MInt orientation, const MFloat A,
                                                               const MFloat* const vars0, const MFloat* const vars1,
                                                               const MFloat* const slope0, const MFloat* const slope1,
                                                               const MFloat f0, const MFloat f1, MFloat* const flux,
                                                               const MFloat mue_wm) {
  static const MFloat F1BRe0 = F1 / m_Re0;


  // calculate the primitve variables on the surface u,v,rho,p
  const std::array<MFloat, nDim> velocity = [&] {
    std::array<MFloat, nDim> tmp{};
    for(MUint n = 0; n < nDim; n++) {
      tmp[n] = F1B2 * (vars0[PV->VV[n]] + vars1[PV->VV[n]]);
    }
    return tmp;
  }();

  // Indices for the orientations
  const MUint id0 = orientation;
  const MUint id1 = index0[orientation];
  const MUint id2 = index1[orientation];

  // Compute A / Re
  const MFloat dAOverRe = A * F1BRe0;

  // Compute the stress terms
  const MUint s00 = id0 * nDim + id0;
  const MUint s01 = id0 * nDim + id1;
  const MUint s10 = id1 * nDim + id0;
  const MUint s11 = id1 * nDim + id1;

  std::array<MFloat, nDim> tau{};
  IF_CONSTEXPR(nDim == 2) {
    tau[id0] =
        mue_wm * (f0 * (F4B3 * slope0[s00] - F2B3 * slope0[s11]) + f1 * (F4B3 * slope1[s00] - F2B3 * slope1[s11]));
    tau[id1] = mue_wm * (f0 * (slope0[s01] + slope0[s10]) + f1 * (slope1[s01] + slope1[s10]));
  }
  else IF_CONSTEXPR(nDim == 3) {
    const MUint s22 = id2 * nDim + id2;
    const MUint s02 = id0 * nDim + id2;
    const MUint s20 = id2 * nDim + id0;
    tau[id0] = mue_wm
               * (f0 * (F4B3 * slope0[s00] - F2B3 * (slope0[s11] + slope0[s22]))
                  + f1 * (F4B3 * slope1[s00] - F2B3 * (slope1[s11] + slope1[s22])));
    tau[id1] = mue_wm * (f0 * (slope0[s01] + slope0[s10]) + f1 * (slope1[s01] + slope1[s10]));
    tau[id2] = mue_wm * (f0 * (slope0[s02] + slope0[s20]) + f1 * (slope1[s02] + slope1[s20]));
  }

  // Correct the fluxes
  for(MUint n = 0; n < nDim; n++) {
    flux[CV->RHO_VV[n]] -= dAOverRe * tau[n];
  }
  flux[CV->RHO_E] -= dAOverRe * (std::inner_product(velocity.begin(), velocity.end(), tau.begin(), 0.0));
}

template <MInt nDim>
inline void FvSysEqnNS<nDim>::viscousFluxThreePoint(const MInt orientation, const MFloat A, const MBool isBndry,
                                                    const MFloat* const surfaceCoords, const MFloat* const coord0,
                                                    const MFloat* const coord1, const MFloat* const cellVars0,
                                                    const MFloat* const cellVars1, const MFloat* const vars0,
                                                    const MFloat* const vars1, const MFloat* const slope0,
                                                    const MFloat* const slope1, const MFloat f0, const MFloat f1,
                                                    MFloat* const flux) {
  static const MFloat rhoUDth = m_muInfinity * m_F1BPr;
  static const MFloat F1BRe0 = F1 / m_Re0;

  // calculate the primitve variables on the surface u,v,rho,p
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

  // Temperature on the surface T = gamma * p / rho
  const MFloat T = m_gamma * p * Frho;

  // Indices for the orientations
  const MUint id0 = orientation;
  const MUint id1 = index0[orientation];
  const MUint id2 = index1[orientation];

  // Compute A / Re
  const MFloat dAOverRe = A * F1BRe0;

  // calculate the viscosity with the sutherland law mue = T^3/2 * (1+S/T_0)(T + S/T_0)
  const MFloat mue = (T * sqrt(T) * m_sutherlandPlusOne) / (T + m_sutherlandConstant);

  const MFloat lambda = (T * sqrt(T) * m_sutherlandPlusOneThermal) / (T + m_sutherlandConstantThermal);

  // TODO labels:FV @Julian, change to stack alloc
  std::vector<MFloat> slopes(PV->noVariables * nDim);

  /*#ifdef _VISCOUS_SLOPES_COMPACT_VARIANT_
      const MFloat dx = F1B2 * (coord1[id0] - coord0[id0]);
  #endif*/
  for(MInt var = 0; var < PV->noVariables; var++) {
    slopes[nDim * var + id0] =
        cellVars1[var]
        - cellVars0[var]
        /*#ifdef _VISCOUS_SLOPES_COMPACT_VARIANT_
                                          + (surfaceCoords[id0] - coord1[id0] + dx) * slope1[var * nDim + id0]
                                          - (surfaceCoords[id0] - coord0[id0] - dx) * slope0[var * nDim + id0]
        #endif*/
        + (surfaceCoords[id1] - coord1[id1]) * slope1[var * nDim + id1]
        - (surfaceCoords[id1] - coord0[id1]) * slope0[var * nDim + id1];
    IF_CONSTEXPR(nDim == 3) {
      slopes[nDim * var + id0] += (surfaceCoords[id2] - coord1[id2]) * slope1[var * nDim + id2]
                                  - (surfaceCoords[id2] - coord0[id2]) * slope0[var * nDim + id2];
    }
    /*#ifdef _VISCOUS_SLOPES_COMPACT_VARIANT_
                                      )
                                     / (F2 * dx);
    #else*/
    // )
    slopes[nDim * var + id0] /= (coord1[id0] - coord0[id0]);


    //#endif
    slopes[nDim * var + id1] = f0 * slope0[var * nDim + id1] + f1 * slope1[var * nDim + id1];
    IF_CONSTEXPR(nDim == 3) {
      slopes[nDim * var + id2] = f0 * slope0[var * nDim + id2] + f1 * slope1[var * nDim + id2];
    }
    if(isBndry) {
      slopes[nDim * var + id0] = f0 * slope0[var * nDim + id0] + f1 * slope1[var * nDim + id0];
    }
  }

  const MUint sq0 = PV->P * nDim + id0;
  const MUint sq1 = PV->RHO * nDim + id0;

  const MFloat q = lambda * m_gFGMOrPr * Frho * (slopes[sq0] - p * Frho * slopes[sq1]);

  // Compute the stress terms
  const MUint s00 = id0 * nDim + id0;
  const MUint s01 = id0 * nDim + id1;
  const MUint s10 = id1 * nDim + id0;
  const MUint s11 = id1 * nDim + id1;

  std::array<MFloat, nDim> tau{};
  IF_CONSTEXPR(nDim == 2) {
    tau[id0] = mue * (F4B3 * slopes[s00] - F2B3 * slopes[s11]);
    tau[id1] = mue * (slopes[s01] + slopes[s10]);
  }
  else IF_CONSTEXPR(nDim == 3) {
    const MUint s22 = id2 * nDim + id2;
    const MUint s02 = id0 * nDim + id2;
    const MUint s20 = id2 * nDim + id0;
    tau[id0] = mue * (F4B3 * slopes[s00] - F2B3 * (slopes[s11] + slopes[s22]));
    tau[id1] = mue * (slopes[s01] + slopes[s10]);
    tau[id2] = mue * (slopes[s02] + slopes[s20]);
  }

  // Compute the flux
  for(MUint n = 0; n < nDim; n++) {
    flux[FV->RHO_VV[n]] -= dAOverRe * tau[n];
  }
  flux[FV->RHO_E] -= dAOverRe * (std::inner_product(velocity.begin(), velocity.end(), tau.begin(), 0.0) + q);

  // progress variable
  const MFloat c = dAOverRe * rhoUDth;
  for(MUint s = 0; s < PV->m_noSpecies; ++s) {
    const MUint i = PV->Y[s] * nDim + id0;
    flux[FV->RHO_Y[s]] -= c * (f0 * slope0[i] + f1 * slope1[i]);
  }
}

/**
 * \brief Computes the viscous fluxes using a five-point stencil (less dissipative) blended with a compact stencil
 * (increased stability). Short: (1-enhanceThreePointViscFluxFactor)*FIVE_POINT +
 * enhanceThreePointViscFluxFactor*THREE_POINT Default for centralizeViscousFlux is 0.1 if not defined in property file
 *
 *
 * \author Lennart Schneiders, Konstantin Froehlich, refactored and moved here: Julian Vorspohl
 */
template <MInt nDim>
inline void FvSysEqnNS<nDim>::viscousFluxStabilized(const MInt orientation, const MFloat A, const MBool isBndry,
                                                    const MFloat* const surfaceCoords, const MFloat* const coord0,
                                                    const MFloat* const coord1, const MFloat* const cellVars0,
                                                    const MFloat* const cellVars1, const MFloat* const vars0,
                                                    const MFloat* const vars1, const MFloat* const slope0,
                                                    const MFloat* const slope1, const MFloat f0, const MFloat f1,
                                                    MFloat* const flux) {
  static const MFloat rhoUDth = m_muInfinity * m_F1BPr;
  static const MFloat F1BRe0 = F1 / m_Re0;
  const MFloat fac = m_enhanceThreePointViscFluxFactor; // increased stability

  // calculate the primitve variables on the surface u,v,rho,p
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

  // Temperature on the surface T = gamma * p / rho
  const MFloat T = m_gamma * p * Frho;

  // Indices for the orientations
  const MUint id0 = orientation;
  const MUint id1 = index0[orientation];
  const MUint id2 = index1[orientation];

  // Compute A / Re
  const MFloat dAOverRe = A * F1BRe0;

  // calculate the viscosity with the sutherland law mue = T^3/2 * (1+S/T_0)(T + S/T_0)
  const MFloat mue = (T * sqrt(T) * m_sutherlandPlusOne) / (T + m_sutherlandConstant);

  const MFloat lambda = (T * sqrt(T) * m_sutherlandPlusOneThermal) / (T + m_sutherlandConstantThermal);

  // TODO labels:FV @Julian, change to stack alloc
  std::vector<MFloat> slopes(PV->noVariables * nDim);

  /*#ifdef _VISCOUS_SLOPES_COMPACT_VARIANT_
      const MFloat dx = F1B2 * (coord1[id0] - coord0[id0]);
  #endif*/
  for(MInt var = 0; var < PV->noVariables; var++) {
    slopes[nDim * var + id0] =
        cellVars1[var]
        - cellVars0[var]
        /*#ifdef _VISCOUS_SLOPES_COMPACT_VARIANT_
                                          + (surfaceCoords[id0] - coord1[id0] + dx) * slope1[var * nDim + id0]
                                          - (surfaceCoords[id0] - coord0[id0] - dx) * slope0[var * nDim + id0]
        #endif*/
        + (surfaceCoords[id1] - coord1[id1]) * slope1[var * nDim + id1]
        - (surfaceCoords[id1] - coord0[id1]) * slope0[var * nDim + id1];
    IF_CONSTEXPR(nDim == 3) {
      slopes[nDim * var + id0] += (surfaceCoords[id2] - coord1[id2]) * slope1[var * nDim + id2]
                                  - (surfaceCoords[id2] - coord0[id2]) * slope0[var * nDim + id2];
    }
    /*#ifdef _VISCOUS_SLOPES_COMPACT_VARIANT_
                                      )
                                     / (F2 * dx);
    #else*/
    // )
    slopes[nDim * var + id0] /= (coord1[id0] - coord0[id0]);


    //#endif
    slopes[nDim * var + id1] = f0 * slope0[var * nDim + id1] + f1 * slope1[var * nDim + id1];
    IF_CONSTEXPR(nDim == 3) {
      slopes[nDim * var + id2] = f0 * slope0[var * nDim + id2] + f1 * slope1[var * nDim + id2];
    }

    slopes[nDim * var + id0] =
        fac * slopes[nDim * var + id0] + (1.0 - fac) * (f0 * slope0[var * nDim + id0] + f1 * slope1[var * nDim + id0]);

    if(isBndry) {
      slopes[nDim * var + id0] = f0 * slope0[var * nDim + id0] + f1 * slope1[var * nDim + id0];
    }
  }

  const MUint sq0 = PV->P * nDim + id0;
  const MUint sq1 = PV->RHO * nDim + id0;

  const MFloat q = lambda * m_gFGMOrPr * Frho * (slopes[sq0] - p * Frho * slopes[sq1]);

  // Compute the stress terms
  const MUint s00 = id0 * nDim + id0;
  const MUint s01 = id0 * nDim + id1;
  const MUint s10 = id1 * nDim + id0;
  const MUint s11 = id1 * nDim + id1;

  std::array<MFloat, nDim> tau{};
  IF_CONSTEXPR(nDim == 2) {
    tau[id0] = mue * (F4B3 * slopes[s00] - F2B3 * slopes[s11]);
    tau[id1] = mue * (slopes[s01] + slopes[s10]);
  }
  else IF_CONSTEXPR(nDim == 3) {
    const MUint s22 = id2 * nDim + id2;
    const MUint s02 = id0 * nDim + id2;
    const MUint s20 = id2 * nDim + id0;
    tau[id0] = mue * (F4B3 * slopes[s00] - F2B3 * (slopes[s11] + slopes[s22]));
    tau[id1] = mue * (slopes[s01] + slopes[s10]);
    tau[id2] = mue * (slopes[s02] + slopes[s20]);
  }

  // Compute the flux
  for(MUint n = 0; n < nDim; n++) {
    flux[FV->RHO_VV[n]] -= dAOverRe * tau[n];
  }
  flux[FV->RHO_E] -= dAOverRe * (std::inner_product(velocity.begin(), velocity.end(), tau.begin(), 0.0) + q);

  // progress variable
  const MFloat c = dAOverRe * rhoUDth;
  for(MUint s = 0; s < PV->m_noSpecies; ++s) {
    const MUint i = PV->Y[s] * nDim + id0;
    flux[FV->RHO_Y[s]] -= c * (f0 * slope0[i] + f1 * slope1[i]);
  }
}

template <MInt nDim>
inline void FvSysEqnNS<nDim>::computePrimitiveVariables(const MFloat* const cvarsCell, MFloat* const pvarsCell,
                                                        const MFloat* const NotUsed(avarsCell)) {
  const MFloat fRho = F1 / cvarsCell[CV->RHO];
  MFloat velPOW2 = F0;
  for(MInt n = 0; n < nDim; ++n) { // compute velocity
    pvarsCell[PV->VV[n]] = cvarsCell[CV->RHO_VV[n]] * fRho;
    velPOW2 += POW2(pvarsCell[PV->VV[n]]);
  }

  // density and pressure:
  pvarsCell[PV->RHO] = cvarsCell[CV->RHO]; // density
  pvarsCell[PV->P] = m_gammaMinusOne * (cvarsCell[CV->RHO_E] - F1B2 * pvarsCell[PV->RHO] * velPOW2);

  // compute the species
  for(MUint s = 0; s < PV->m_noSpecies; ++s) {
    pvarsCell[PV->Y[s]] = cvarsCell[CV->RHO_Y[s]] * fRho;
  }
}

template <MInt nDim>
inline void FvSysEqnNS<nDim>::computeConservativeVariables(const MFloat* const pvarsCell, MFloat* const cvarsCell,
                                                           const MFloat* const NotUsed(avarsCell)) {
  // compute rho v
  MFloat velPOW2 = F0;
  for(MInt vel = 0; vel < nDim; ++vel) {
    cvarsCell[vel] = pvarsCell[vel] * pvarsCell[PV->RHO];
    velPOW2 += POW2(pvarsCell[vel]);
  }

  // compute rho e
  cvarsCell[CV->RHO_E] = pvarsCell[PV->P] * m_F1BGammaMinusOne + F1B2 * pvarsCell[PV->RHO] * velPOW2;

  // copy the density
  cvarsCell[CV->RHO] = pvarsCell[PV->RHO];

  // compute rho Yi
  for(MUint s = 0; s < PV->m_noSpecies; ++s) {
    cvarsCell[CV->RHO_Y[s]] = pvarsCell[PV->Y[s]] * pvarsCell[PV->RHO];
  }
}

template <MInt nDim>
inline std::vector<std::vector<MFloat>>
FvSysEqnNS<nDim>::conservativeSlopes(const MFloat* const pvarsCell, const MFloat* const NotUsed(cvarsCell),
                                     const MFloat* const NotUsed(avarsCell), const MFloat* const slopesCell) {
  MFloat U2 = F0;
  for(MInt i = 0; i < nDim; i++) {
    U2 += POW2(pvarsCell[PV->VV[i]]);
  }

  std::vector<std::vector<MFloat>> dQ(CV->noVariables, std::vector<MFloat>(nDim));

  for(MInt d = 0; d < nDim; d++) {
    dQ[CV->RHO][d] = slopesCell[PV->RHO * nDim + d];
    dQ[CV->RHO_E][d] = slopesCell[PV->P * nDim + d] / (m_gamma - F1) + F1B2 * U2 * slopesCell[PV->RHO * nDim + d];
    for(MInt j = 0; j < nDim; j++) {
      dQ[CV->RHO_VV[j]][d] =
          pvarsCell[PV->VV[j]] * slopesCell[PV->RHO * nDim + d] + pvarsCell[PV->RHO] * slopesCell[PV->VV[j] * nDim + d];
      dQ[CV->RHO_E][d] += pvarsCell[PV->RHO] * pvarsCell[PV->VV[j]] * slopesCell[PV->VV[j] * nDim + d];
    }
    for(MUint s = 0; s < CV->m_noSpecies; ++s) {
      dQ[CV->RHO_Y[s]][d] =
          pvarsCell[PV->Y[s]] * slopesCell[PV->RHO * nDim + d] + pvarsCell[PV->RHO] * slopesCell[PV->Y[s] * nDim + d];
    }
  }

  return dQ;
}


#endif // FVSYSEQNNS_H_
