// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef FVSYSEQNEEGAS_H_
#define FVSYSEQNEEGAS_H_

#include <algorithm>
#include <array>
#include <iterator>
#include <numeric>
#include "INCLUDE/maiaconstants.h"
#include "INCLUDE/maiatypes.h"
#include "MEMORY/alloc.h"
#include "filter.h"
#include "fvcartesiansyseqnns.h"


template <MInt nDim>
class FvSysEqnEEGas : public FvSysEqnNS<nDim> {
  // Declare parent a friend so that CRTP can access private methods/members
  // friend class FvSysEqn<nDim>;

 public:
  FvSysEqnEEGas(const MInt solverId, const MInt noSpecies);

  template <MInt stencil = AUSM>
  inline void Ausm(const MInt orientation, const MFloat upwindCoefficient, const MFloat A, const MFloat* const leftVars,
                   const MFloat* const rightVars, const MFloat* const NotUsed(srfcCoeff), MFloat* const flux);

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
      TERMM(1, "Three point viscous Flux Scheme not implemented for SysEqnEEGas.");
    } else {
      TERMM(1, "Viscous Flux Scheme not implemented.");
    }
  }

  inline void viscousFluxFivePoint(const MInt orientation, const MFloat A, const MFloat* const vars0,
                                   const MFloat* const vars1, const MFloat* const slope0, const MFloat* const slope1,
                                   const MFloat* const srfcCoeff, const MFloat f0, const MFloat f1, MFloat* const flux);

  inline void computePrimitiveVariables(const MFloat* const cvarsCell, MFloat* const pvarsCell,
                                        const MFloat* const avarsCell);

  inline void computeConservativeVariables(const MFloat* const pvarsCell, MFloat* const cvarsCell,
                                           const MFloat* const NotUsed(avarsCell));

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

  void readProperties();
  MFloat m_EEGasEps;

 public:
  using Base::m_muInfinity;
  using Base::m_Re0;

  // Hold indices for primitive and conservative variables
  struct ConservativeVariables;
  struct FluxVariables;
  struct PrimitiveVariables;
  struct AdditionalVariables;
  static constexpr MBool hasAV = true;
  static constexpr MBool hasSC = false;
  ConservativeVariables* CV = nullptr;
  FluxVariables* FV = nullptr;
  PrimitiveVariables* PV = nullptr;
  AdditionalVariables* AV = nullptr;
};


/// \brief Static indices for accessing conservative variables
/// in nDim spatial dimensions
template <MInt nDim>
struct FvSysEqnEEGas<nDim>::ConservativeVariables {
  static const MInt Segfault = std::numeric_limits<MInt>::min();

  static constexpr MUint m_noSpecies = 1;
  static constexpr MInt noVariables = nDim + 1;

  // these are the actual conservative variables. A is the void fraction alpha.
  static constexpr MInt A_RHO_U = 0;
  static constexpr MInt A_RHO_V = 1;
  static constexpr MInt A_RHO_W = nDim == 3 ? 2 : Segfault;
  static constexpr std::array<MInt, nDim> A_RHO_VV = getArray012();
  static constexpr MInt A_RHO = nDim;

  // the other ones are neccesary to compile the code.
  static constexpr MInt RHO_U = 0;
  static constexpr MInt RHO_V = 1;
  static constexpr MInt RHO_W = nDim == 3 ? 2 : Segfault;
  static constexpr std::array<MInt, nDim> RHO_VV = getArray012();
  static constexpr MInt RHO_E = Segfault;
  static constexpr MInt RHO = nDim;

  static constexpr MInt RHO_Y[1] = {nDim};
};

/// \brief Static indices for accessing flux variables
template <MInt nDim>
struct FvSysEqnEEGas<nDim>::FluxVariables : ConservativeVariables {
  // the Flux variables with P are used for the pressure fluxes in the gas-liquid E-E equations
  static constexpr MInt P_RHO_U = ConservativeVariables::noVariables + 0;
  static constexpr MInt P_RHO_V = ConservativeVariables::noVariables + 1;
  static constexpr MInt P_RHO_W = nDim == 3 ? ConservativeVariables::noVariables + 2 : ConservativeVariables::Segfault;
  static constexpr std::array<MInt, nDim> P_RHO_VV = [] {
    IF_CONSTEXPR(nDim == 2) {
      std::array<MInt, 2> a = {ConservativeVariables::noVariables + 0, ConservativeVariables::noVariables + 1};
      return a;
    }
    else {
      std::array<MInt, 3> a = {ConservativeVariables::noVariables + 0, ConservativeVariables::noVariables + 1,
                               ConservativeVariables::noVariables + 2};
      return a;
    }
  }();
  static constexpr MInt noVariables = ConservativeVariables::noVariables + 3;
};

/// \brief Static indices for accessing primitive variables
/// in nDim spatial dimensions
template <MInt nDim>
struct FvSysEqnEEGas<nDim>::PrimitiveVariables {
  static const MInt Segfault = std::numeric_limits<MInt>::min();
  static constexpr MUint m_noSpecies = 1;
  static constexpr MInt noVariables = nDim + 3;

  static constexpr MInt U = 0;
  static constexpr MInt V = 1;
  static constexpr MInt W = nDim == 3 ? 2 : Segfault;
  static constexpr std::array<MInt, nDim> VV = getArray012();
  static constexpr MInt RHO = nDim;
  static constexpr MInt P = nDim + 1;

  // A is the void fraction alpha
  static constexpr MInt A = nDim + 2;
  static constexpr MInt Y[1] = {A};

  const std::array<MString, nDim + 3> varNames = [] {
    IF_CONSTEXPR(nDim == 2) {
      std::array<MString, nDim + 3> a = {"u", "v", "rho", "p", "alpha"};
      return a;
    }
    else {
      std::array<MString, nDim + 3> a = {"u", "v", "w", "rho", "p", "alpha"};
      return a;
    }
  }();
  void getPrimitiveVariableNames(MString* names) {
    for(MInt i = 0; i < noVariables; i++) {
      names[i] = varNames[i];
    }
  };
};

/// \brief Static indices for accessing additional variables
template <MInt nDim>
struct FvSysEqnEEGas<nDim>::AdditionalVariables {
  static constexpr MInt noVariables = 1;

  static constexpr MInt DC = 0; // DepthCorrectionValues
};

template <MInt nDim>
template <MInt stencil>
inline void FvSysEqnEEGas<nDim>::Ausm(const MInt orientation, const MFloat upwindCoefficient, const MFloat A,
                                      const MFloat* const leftVars, const MFloat* const rightVars,
                                      const MFloat* const NotUsed(srfcCoeff), MFloat* const flux) {
  // catch the primitive variables rho and p,
  // compute speed of sound, and interface mach number
  const MFloat RHOL = leftVars[PV->RHO];
  const MFloat PL = leftVars[PV->P];
  const MFloat AL = sqrt(m_gamma * mMax(MFloatEps, PL / mMax(MFloatEps, RHOL)));
  const MFloat ML = leftVars[orientation] / AL;

  const MFloat ALPHAL = mMax(mMin(leftVars[PV->Y[0]], F1), F0);

  const MFloat RHOR = rightVars[PV->RHO];
  const MFloat PR = rightVars[PV->P];
  const MFloat AR = sqrt(m_gamma * mMax(MFloatEps, PR / mMax(MFloatEps, RHOR)));
  const MFloat MR = rightVars[orientation] / AR;

  const MFloat ALPHAR = mMax(mMin(rightVars[PV->Y[0]], F1), F0);

  // calculation of the resulting pressure and mach number on the surface
  const MFloat MLR = 0.5 * (ML + MR);
  const MFloat PLR = PL * (0.5 + upwindCoefficient * ML) + PR * (0.5 - upwindCoefficient * MR);

  const MFloat ALPHALR = 0.5 * (ALPHAL + ALPHAR);

  // calculation of the left and right rho*a
  const MFloat RHO_AL = RHOL * AL;
  const MFloat RHO_AR = RHOR * AR;

  // calculation of the resulting mass flux through the surface
  const MFloat RHO_U2 = 0.25 * (MLR * (RHO_AL + RHO_AR) + fabs(MLR) * (RHO_AL - RHO_AR));
  const MFloat AbsRHO_U2 = fabs(RHO_U2);

  std::array<MFloat, nDim> pFactor{};
  pFactor[orientation] = 1.0;

  for(MUint n = 0; n < nDim; n++) {
    flux[FV->A_RHO_VV[n]] = ALPHALR
                            * (RHO_U2 * (leftVars[PV->VV[n]] + rightVars[PV->VV[n]])
                               + AbsRHO_U2 * (leftVars[PV->VV[n]] - rightVars[PV->VV[n]]))
                            * A;
    flux[FV->P_RHO_VV[n]] = PLR * pFactor[n] * A;
  }

  flux[FV->A_RHO] = ALPHALR * 2.0 * RHO_U2 * A;
}

template <MInt nDim>
inline void FvSysEqnEEGas<nDim>::AusmBndryCorrection(const MInt /* orientation */, const MFloat /* A */,
                                                     const MFloat* const /* leftVars */,
                                                     const MFloat* const /* rightVars */, MFloat* const flux) {
  for(MInt n = 0; n < nDim; n++) {
    flux[FV->A_RHO_VV[n]] = F0;
  }
  flux[FV->A_RHO] = F0;
}

template <MInt nDim>
inline void FvSysEqnEEGas<nDim>::viscousFluxFivePoint(const MInt orientation, const MFloat A, const MFloat* const vars0,
                                                      const MFloat* const vars1, const MFloat* const slope0,
                                                      const MFloat* const slope1,
                                                      const MFloat* const NotUsed(srfcCoeff), const MFloat f0,
                                                      const MFloat f1, MFloat* const flux) {
  static const MFloat F1BRe0 = F1 / m_Re0;

  const MFloat rho = F1B2 * (vars0[PV->RHO] + vars1[PV->RHO]);
  const MFloat F1Brho = F1 / rho;
  const MFloat p = F1B2 * (vars0[PV->P] + vars1[PV->P]);
  const MFloat alpha = F1B2 * (vars0[PV->A] + vars1[PV->A]);

  // Temperature on the surface T = gamma * p / rho
  const MFloat T = m_gamma * p * F1Brho;

  // Indices for the orientations
  const MUint id0 = orientation;
  const MUint id1 = index0[orientation];
  const MUint id2 = index1[orientation];

  // Compute A / Re
  const MFloat dAOverRe = A * F1BRe0;

  // calculate the viscosity with the sutherland law mue = T^3/2 * (1+S/T_0)(T + S/T_0)
  const MFloat mue = (T * sqrt(T) * m_sutherlandPlusOne) / (T + m_sutherlandConstant);

  // Compute the stress terms
  const MUint s00 = id0 * nDim + id0;
  const MUint s01 = id0 * nDim + id1;
  const MUint s10 = id1 * nDim + id0;
  const MUint s11 = id1 * nDim + id1;

  std::array<MFloat, nDim> tau{};
  if(nDim == 2) {
    tau[id0] =
        mue * alpha * (f0 * (F4B3 * slope0[s00] - F2B3 * slope0[s11]) + f1 * (F4B3 * slope1[s00] - F2B3 * slope1[s11]));
    tau[id1] = mue * alpha * (f0 * (slope0[s01] + slope0[s10]) + f1 * (slope1[s01] + slope1[s10]));
  } else if(nDim == 3) {
    const MUint s22 = id2 * nDim + id2;
    const MUint s02 = id0 * nDim + id2;
    const MUint s20 = id2 * nDim + id0;
    tau[id0] = mue * alpha
               * (f0 * (F4B3 * slope0[s00] - F2B3 * (slope0[s11] + slope0[s22]))
                  + f1 * (F4B3 * slope1[s00] - F2B3 * (slope1[s11] + slope1[s22])));
    tau[id1] = mue * alpha * (f0 * (slope0[s01] + slope0[s10]) + f1 * (slope1[s01] + slope1[s10]));
    tau[id2] = mue * alpha * (f0 * (slope0[s02] + slope0[s20]) + f1 * (slope1[s02] + slope1[s20]));
  }

  // Compute the flux
  for(MUint n = 0; n < nDim; n++) {
    flux[FV->A_RHO_VV[n]] -= dAOverRe * tau[n];
  }
}

template <MInt nDim>
inline void FvSysEqnEEGas<nDim>::computePrimitiveVariables(const MFloat* const cvarsCell, MFloat* const pvarsCell,
                                                           const MFloat* const avarsCell) {
  MFloat fRhoAlpha = F0;
  if(cvarsCell[CV->A_RHO] > -m_EEGasEps) {
    fRhoAlpha = F1 / mMax(cvarsCell[CV->A_RHO], m_EEGasEps);
  } else {
    fRhoAlpha = F1 / cvarsCell[CV->A_RHO];
  }

  for(MInt vel = 0; vel < nDim; ++vel) {
    pvarsCell[PV->VV[vel]] = cvarsCell[CV->A_RHO_VV[vel]] * fRhoAlpha;
  }
  pvarsCell[PV->RHO] = pvarsCell[PV->P] * m_gamma + avarsCell[AV->DC]; // density

  pvarsCell[PV->A] = cvarsCell[CV->A_RHO] / pvarsCell[PV->RHO]; // alpha
}

template <MInt nDim>
inline void FvSysEqnEEGas<nDim>::computeConservativeVariables(const MFloat* const pvarsCell,
                                                              MFloat* const cvarsCell,
                                                              const MFloat* const NotUsed(avarsCell)) {
  cvarsCell[CV->A_RHO] = pvarsCell[PV->RHO] * pvarsCell[PV->A];
  for(MInt vel = 0; vel < nDim; vel++) {
    cvarsCell[CV->A_RHO_VV[vel]] = cvarsCell[CV->A_RHO] * pvarsCell[PV->VV[vel]];
  }
}

template <MInt nDim>
inline std::vector<std::vector<MFloat>>
FvSysEqnEEGas<nDim>::conservativeSlopes(const MFloat* const NotUsed(pvarsCell), const MFloat* const NotUsed(cvarsCell),
                                        const MFloat* const NotUsed(avarsCell),
                                        const MFloat* const NotUsed(slopesCell)) {
  mTerm(1, AT_,
        "SysEqnEEGas has not been tested with adaptation. The computation of the conservative variables should be "
        "adjusted. Terminating.");
  std::vector<std::vector<MFloat>> dQ(CV->noVariables, std::vector<MFloat>(nDim));
  return dQ;
}

#endif // FVSYSEQNEEGAS_H_
