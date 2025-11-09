// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef FVSYSEQNRANS_H_
#define FVSYSEQNRANS_H_

#include <algorithm>
#include <array>
#include <iterator>
#include <limits>
#include <numeric>
#include "INCLUDE/maiaconstants.h"
#include "INCLUDE/maiatypes.h"
#include "MEMORY/alloc.h"
#include "filter.h"
#include "fvcartesiansyseqnns.h"
#include "fvransmodelconstants.h"

// Only needed (really: possible to use) in single-threaded applications
#ifndef _OPENMP
#include "MEMORY/scratch.h"
#endif


template <MInt nDim, class RANSModel>
class FvSysEqnRANS : public FvSysEqnNS<nDim> {
  // Declare parent a friend so that CRTP can access private methods/members
  // friend class FvSysEqn<nDim>;

 public:
  FvSysEqnRANS(const MInt solverId, const MInt noSpecies);

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
  inline void viscousFlux(const MInt orientation, const MFloat A, const MBool isBndry,
                          const MFloat* const surfaceCoords, const MFloat* const coord0, const MFloat* const coord1,
                          const MFloat* const cellVars0, const MFloat* const cellVars1, const MFloat* const vars0,
                          const MFloat* const vars1, const MFloat* const slope0, const MFloat* const slope1,
                          const MFloat f0, const MFloat f1, MFloat* const flux) {
    if(stencil == THREE_POINT) {
      viscousFluxThreePoint(orientation, A, isBndry, surfaceCoords, coord0, coord1, cellVars0, cellVars1, vars0, vars1,
                            slope0, slope1, f0, f1, flux);
    } else {
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

  inline void computePrimitiveVariables(const MFloat* const cvarsCell, MFloat* const pvarsCell,
                                        const MFloat* const NotUsed(avarsCell));

  inline void computeConservativeVariables(const MFloat* const pvarsCell, MFloat* const cvarsCell,
                                           const MFloat* const NotUsed(avarsCell));

  inline std::vector<std::vector<MFloat>> conservativeSlopes(const MFloat* const pvarsCell,
                                                             const MFloat* const cvarsCell,
                                                             const MFloat* const avarsCell,
                                                             const MFloat* const slopesCell);

  void computeVolumeForces(const MInt cellId, MFloat* vars, MFloat* rhs, const MFloat cellVol,
                           const MFloat* const slopes, const MInt recData, const MInt* const recNghbr,
                           const MFloat* const recConst, const MInt noRecNghbr, MFloat dist) {
    IF_CONSTEXPR(m_ransModel == RANS_SA_DV) {
      computeVolumeForcesRANS_SA(cellId, vars, rhs, cellVol, slopes, recData, recNghbr, recConst, noRecNghbr, dist);
    }
    IF_CONSTEXPR(m_ransModel == RANS_FS) {
      computeVolumeForcesRANS_FS(cellId, vars, rhs, cellVol, slopes, recData, recNghbr, recConst, noRecNghbr, dist);
    }
    IF_CONSTEXPR(m_ransModel == RANS_KOMEGA) {
      computeVolumeForcesRANS_KOMEGA(cellId, vars, rhs, cellVol, slopes, recData, recNghbr, recConst, noRecNghbr, dist);
    }
  }

  void computeVolumeForcesRANS_SA(const MInt cellId, MFloat* vars, MFloat* rhs, const MFloat cellVol,
                                  const MFloat* const slopes, const MInt recData, const MInt* const /* recNghbr */,
                                  const MFloat* const /* recConst */, const MInt noRecNghbr, MFloat dist);
  void computeVolumeForcesRANS_FS(const MInt cellId, MFloat* vars, MFloat* rhs, MFloat cellVol,
                                  const MFloat* const slopes, MInt recData, const MInt* const recNghbr,
                                  const MFloat* const recConst, MInt noRecNghbr, MFloat /* dist */);
  void computeVolumeForcesRANS_KOMEGA(const MInt cellId, MFloat* vars, MFloat* rhs, const MFloat cellVol,
                                      const MFloat* const slopes, const MInt /* recData */,
                                      const MInt* const /* recNghbr */, const MFloat* const /* recConst */,
                                      const MInt /* noRecNghbr */, MFloat /* dist */);

 public:
  static constexpr MInt m_noRansEquations = RANSModel::noRansEquations;
  static constexpr MInt m_ransModel = RANSModel::ransModel;

 private:
  using Base = FvSysEqnNS<nDim>;
  using Base::index0;
  using Base::index1;
  using Base::m_noSpecies;
  using Base::m_solverId;

  // Equation specific values
  using Base::m_F1BGammaMinusOne;
  using Base::m_F1BPr;
  using Base::m_gamma;
  using Base::m_gammaMinusOne;
  using Base::m_gFGMOrPr;
  using Base::m_Pr;
  using Base::m_referenceTemperature;
  using Base::m_sutherlandConstant;
  using Base::m_sutherlandConstantThermal;
  using Base::m_sutherlandPlusOne;
  using Base::m_sutherlandPlusOneThermal;

  using Base::getArray012;

  const MFloat m_Pr_t = 0.9;

  MBool m_SACurvature = false;

 public:
  using Base::m_muInfinity;
  using Base::m_Re0;

  // Hold indices for primitive and conservative variables
  struct ConservativeVariables;
  struct FluxVariables;
  struct PrimitiveVariables;
  static constexpr MBool hasSC = false;
  ConservativeVariables* CV = nullptr;
  FluxVariables* FV = nullptr;
  PrimitiveVariables* PV = nullptr;
};


/// \brief Static indices for accessing conservative variables
/// in nDim spatial dimensions
template <MInt nDim, class RANSModel>
struct FvSysEqnRANS<nDim, RANSModel>::ConservativeVariables {
  static constexpr MInt m_noRansEquations = RANSModel::noRansEquations;

  static constexpr MInt Segfault = std::numeric_limits<MInt>::min();

  static constexpr MInt RHO_U = 0;
  static constexpr MInt RHO_V = 1;
  static constexpr MInt RHO_W = nDim == 3 ? 2 : Segfault;
  static constexpr std::array<MInt, nDim> RHO_VV = getArray012();
  static constexpr MInt RHO_E = nDim;
  static constexpr MInt RHO = nDim + 1;
  static constexpr MInt RHO_N = nDim + 2;
  static constexpr MInt RHO_K = nDim + 2;
  static constexpr MInt RHO_OMEGA = nDim + 3;
  static constexpr MInt RHO_C = nDim + 2 + m_noRansEquations;

  MUint m_noSpecies;
  const MInt noVariables;

  MInt* RHO_NN = nullptr;
  MInt* RHO_Y = nullptr;

  ConservativeVariables(const MInt noSpecies);
  ~ConservativeVariables();
};

template <MInt nDim, class RANSModel>
struct FvSysEqnRANS<nDim, RANSModel>::FluxVariables : ConservativeVariables {
  FluxVariables(const MInt noSpecies);
};

/// \brief Static indices for accessing primitive variables
/// in nDim spatial dimensions
template <MInt nDim, class RANSModel>
struct FvSysEqnRANS<nDim, RANSModel>::PrimitiveVariables {
  static constexpr MInt m_noRansEquations = RANSModel::noRansEquations;

  static const MInt Segfault = std::numeric_limits<MInt>::min();

  static constexpr MInt U = 0;
  static constexpr MInt V = 1;
  static constexpr MInt W = nDim == 3 ? 2 : Segfault;
  static constexpr std::array<MInt, nDim> VV = getArray012();
  static constexpr MInt RHO = nDim;
  static constexpr MInt P = nDim + 1;
  static constexpr MInt T = nDim + 1;

  static constexpr MInt N = nDim + 2;
  static constexpr MInt K = nDim + 2;
  static constexpr MInt OMEGA = nDim + 3;
  static constexpr MInt C = nDim + 2 + m_noRansEquations;

  MUint m_noSpecies;
  const MInt noVariables;

  static std::vector<MString> varNames;

  MInt* NN = nullptr;
  MInt* Y = nullptr;

  PrimitiveVariables(const MInt noSpecies);
  void getPrimitiveVariableNames(MString* names);
  ~PrimitiveVariables();
};

template <MInt nDim, class RANSModel>
template <MInt stencil>
inline void FvSysEqnRANS<nDim, RANSModel>::Ausm(const MInt orientation, const MFloat upwindCoefficient, const MFloat A,
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

  // calculation of the velocity magnitude
  const MFloat U2L = std::inner_product(&leftVars[PV->U], &leftVars[PV->U] + nDim, &leftVars[PV->U], 0.0);
  const MFloat U2R = std::inner_product(&rightVars[PV->U], &rightVars[PV->U] + nDim, &rightVars[PV->U], 0.0);

  const MFloat e0 = PLfRHOL * m_F1BGammaMinusOne + 0.5 * U2L + PLfRHOL;
  const MFloat e1 = PRfRHOR * m_F1BGammaMinusOne + 0.5 * U2R + PRfRHOR;

  std::array<MFloat, nDim> pFactor{};
  pFactor[orientation] = 1.0;

  for(MInt n = 0; n < nDim; n++) {
    flux[FV->RHO_VV[n]] = (RHO_U2 * (leftVars[PV->VV[n]] + rightVars[PV->VV[n]])
                           + AbsRHO_U2 * (leftVars[PV->VV[n]] - rightVars[PV->VV[n]]) + PLR * pFactor[n])
                          * A;
  }

  flux[FV->RHO_E] = (RHO_U2 * (e0 + e1) + AbsRHO_U2 * (e0 - e1)) * A;
  flux[FV->RHO] = 2.0 * RHO_U2 * A;

  // RANS specific flux calculations
  for(MInt r = 0; r < PV->m_noRansEquations; ++r) {
    const MFloat NL = leftVars[PV->NN[r]];
    const MFloat NR = rightVars[PV->NN[r]];
    flux[FV->RHO_NN[r]] = (RHO_U2 * (NL + NR) + AbsRHO_U2 * (NL - NR)) * A;
  }

  // Flux calculation for species transport
  // TODO labels:FV Make noSpecies constexpr
  for(MUint s = 0; s < PV->m_noSpecies; ++s) {
    const MFloat YsL = leftVars[PV->Y[s]];
    const MFloat YsR = rightVars[PV->Y[s]];
    flux[FV->RHO_Y[s]] = (RHO_U2 * (YsL + YsR) + AbsRHO_U2 * (YsL - YsR)) * A;
  }
}

template <MInt nDim, class RANSModel>
inline void FvSysEqnRANS<nDim, RANSModel>::AusmBndryCorrection(const MInt orientation, const MFloat A,
                                                               const MFloat* const leftVars,
                                                               const MFloat* const rightVars, MFloat* const flux) {
  const MFloat PL = leftVars[PV->P];
  const MFloat PR = rightVars[PV->P];
  const MFloat PLR = 0.5 * (PR + PL);

  std::array<MFloat, nDim> pFactor{};
  pFactor[orientation] = 1.0;

  for(MInt n = 0; n < nDim; n++) {
    flux[FV->RHO_VV[n]] = PLR * pFactor[n] * A;
  }

  flux[FV->RHO_E] = 0.0;
  flux[FV->RHO] = 0.0;

  // RANS specific flux calculations
  for(MInt r = 0; r < PV->m_noRansEquations; ++r) {
    flux[FV->RHO_NN[r]] = 0.0;
  }

  // Flux calculation for species transport
  // TODO labels:FV Make noSpecies constexpr
  for(MUint s = 0; s < PV->m_noSpecies; ++s) {
    flux[FV->RHO_Y[s]] = 0.0;
  }
}

template <MInt nDim, class RANSModel>
inline void FvSysEqnRANS<nDim, RANSModel>::viscousFluxFivePoint(const MInt orientation, const MFloat A,
                                                                const MFloat* const vars0, const MFloat* const vars1,
                                                                const MFloat* const slope0, const MFloat* const slope1,
                                                                const MFloat* const NotUsed(srfcCoeff), const MFloat f0,
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

  // Turbulence modelling
  MFloat mue_t = 0.0;

  IF_CONSTEXPR(m_ransModel == RANS_SA_DV) {
    const MFloat nuTilde = F1B2 * (vars0[PV->N] + vars1[PV->N]);
    const MFloat nuLaminar = mue / rho;
    const MFloat chi = nuTilde / (nuLaminar);
    const MFloat fv1 = pow(chi, 3) / (pow(chi, 3) + RM_SA_DV::cv1to3);
    const MFloat nuTurb = fv1 * nuTilde;
    if(nuTilde > 0.0) {
      mue_t = rho * nuTurb;
    }
  }
  IF_CONSTEXPR(m_ransModel == RANS_FS) {
    const MFloat nuTilde = F1B2 * (vars0[PV->N] + vars1[PV->N]);
    const MFloat nuLaminar = mue / rho;
    const MFloat chi = nuTilde / (nuLaminar);
    const MFloat fv1 = pow(chi, 3) / (pow(chi, 3) + RM_FS::facv1to3);
    const MFloat nuTurb = fv1 * nuTilde;
    if(nuTilde > 0.0) {
      mue_t = rho * nuTurb;
    }
  }
  IF_CONSTEXPR(m_ransModel == RANS_KOMEGA) {
    const MFloat k = F1B2 * (vars0[PV->K] + vars1[PV->K]);
    const MFloat o = F1B2 * (vars0[PV->OMEGA] + vars1[PV->OMEGA]);
    const MFloat nuTurb = k / o;
    if(nuTurb > 0.0) {
      mue_t = rho * nuTurb;
    }
  }


  const MFloat lambda_t = mue_t / m_Pr_t;

  // Wall-modeling with precomputed, approximated mue_wm
  MFloat mue_wm = 0.0;
  /*if(m_wmLES) {
    if(*srfcs[srfcId].m_bndryCndId == 3399) {
      MInt cellId = nghbrCellIds[2 * srfcId];
      if(!a_isBndryCell(cellId)) cellId = nghbrCellIds[2 * srfcId + 1];
      if(a_isBndryCell(cellId) && !a_isHalo(cellId)) {
        mue_wm = computeWMViscositySpalding(cellId);
      }
    }
  }*/

  const MFloat lambda = (T * sqrt(T) * m_sutherlandPlusOneThermal) / (T + m_sutherlandConstantThermal);

  const MUint sq0 = PV->P * nDim + id0;
  const MUint sq1 = PV->RHO * nDim + id0;
  const MFloat q = (lambda + lambda_t) * m_gFGMOrPr * Frho
                   * ((f0 * slope0[sq0] + f1 * slope1[sq0]) - p * Frho * (f0 * slope0[sq1] + f1 * slope1[sq1]));

  // Compute the stress terms
  const MUint s00 = id0 * nDim + id0;
  const MUint s01 = id0 * nDim + id1;
  const MUint s10 = id1 * nDim + id0;
  const MUint s11 = id1 * nDim + id1;

  std::array<MFloat, nDim> tau{};
  IF_CONSTEXPR(nDim == 2) {
    tau[id0] = (mue + mue_t + mue_wm)
               * (f0 * (F4B3 * slope0[s00] - F2B3 * slope0[s11]) + f1 * (F4B3 * slope1[s00] - F2B3 * slope1[s11]));
    tau[id1] = (mue + mue_t + mue_wm) * (f0 * (slope0[s01] + slope0[s10]) + f1 * (slope1[s01] + slope1[s10]));
  }
  else IF_CONSTEXPR(nDim == 3) {
    const MUint s22 = id2 * nDim + id2;
    const MUint s02 = id0 * nDim + id2;
    const MUint s20 = id2 * nDim + id0;
    tau[id0] = (mue + mue_t + mue_wm)
               * (f0 * (F4B3 * slope0[s00] - F2B3 * (slope0[s11] + slope0[s22]))
                  + f1 * (F4B3 * slope1[s00] - F2B3 * (slope1[s11] + slope1[s22])));
    tau[id1] = (mue + mue_t + mue_wm) * (f0 * (slope0[s01] + slope0[s10]) + f1 * (slope1[s01] + slope1[s10]));
    tau[id2] = (mue + mue_t + mue_wm) * (f0 * (slope0[s02] + slope0[s20]) + f1 * (slope1[s02] + slope1[s20]));
  }

  // Compute the flux
  for(MUint n = 0; n < nDim; n++) {
    flux[FV->RHO_VV[n]] -= dAOverRe * tau[n];
  }
  flux[FV->RHO_E] -= dAOverRe * (std::inner_product(velocity.begin(), velocity.end(), tau.begin(), 0.0) + q);

  IF_CONSTEXPR(m_ransModel == RANS_SA_DV) {
    const MUint n0 = PV->N * nDim + id0;
    const MUint rho0 = PV->RHO * nDim + id0;

    const MFloat nuTilde = F1B2 * (vars0[PV->N] + vars1[PV->N]);
    // const MFloat nuLaminar = mue/rho;
    // source: Density Corrections for Turbulence Models," Aerospace Science and Technology, Vol. 4, 2000, pp. 1--11
    MFloat viscflux_nu = RM_SA_DV::Fsigma
                         * (mue * (f0 * slope0[n0] + f1 * slope1[n0])
                            + sqrt(rho) * nuTilde
                                  * (sqrt(rho) * (f0 * slope0[n0] + f1 * slope1[n0])
                                     + nuTilde * F1 / (F2 * sqrt(rho)) * (f0 * slope0[rho0] + f1 * slope1[rho0])));

    // no correction model
    // MFloat viscflux_nu = RM_SA_DV::Fsigma * (nuLaminar + nuTilde) * (f0 * slope0[n0] + f1 * slope1[n0]);

    flux[FV->RHO_N] -= dAOverRe * (viscflux_nu);
  }

  IF_CONSTEXPR(m_ransModel == RANS_FS) {
    const MUint n0 = PV->N * nDim + id0;
    const MFloat nuTilde = F1B2 * (vars0[PV->N] + vars1[PV->N]);
    const MFloat nuLaminar = mue / rho;
    MFloat viscflux_nu = rho * (nuLaminar + RM_FS::fasigma * nuTilde) * (f0 * slope0[n0] + f1 * slope1[n0]);
    flux[FV->RHO_N] -= dAOverRe * (viscflux_nu);
  }

  IF_CONSTEXPR(m_ransModel == RANS_KOMEGA) {
    const MUint n0 = PV->N * nDim + id0;
    const MUint rho0 = PV->RHO * nDim + id0;

    const MFloat nuTilde = F1B2 * (vars0[PV->N] + vars1[PV->N]);
    MFloat viscflux_nu = RM_SA_DV::Fsigma
                         * (mue * (f0 * slope0[n0] + f1 * slope1[n0])
                            + sqrt(rho) * nuTilde
                                  * (sqrt(rho) * (f0 * slope0[n0] + f1 * slope1[n0])
                                     + nuTilde * F1 / (F2 * sqrt(rho)) * (f0 * slope0[rho0] + f1 * slope1[rho0])));

    flux[FV->RHO_N] -= dAOverRe * (viscflux_nu);
    const MUint nk = PV->K * nDim + id0;
    const MUint no = PV->OMEGA * nDim + id0;

    const MFloat k = F1B2 * (vars0[PV->K] + vars1[PV->K]);
    const MFloat omega = F1B2 * (vars0[PV->OMEGA] + vars1[PV->OMEGA]);
    const MFloat mueTurb = rho * k / omega;

    MFloat viscflux_k = (mue + RM_KOMEGA::sigma_k * mueTurb) * (f0 * slope0[nk] + f1 * slope1[nk]);
    MFloat viscflux_omega = (mue + RM_KOMEGA::sigma_o * mueTurb) * (f0 * slope0[no] + f1 * slope1[no]);

    flux[FV->RHO_K] -= dAOverRe * (viscflux_k);
    flux[FV->RHO_OMEGA] -= dAOverRe * (viscflux_omega);
  }


  // progress variable
  const MFloat c = dAOverRe * rhoUDth;
  for(MUint s = 0; s < PV->m_noSpecies; ++s) {
    const MUint i = PV->Y[s] * nDim + id0;
    flux[FV->RHO_Y[s]] -= c * (f0 * slope0[i] + f1 * slope1[i]);
  }
}

template <MInt nDim, class RANSModel>
inline void FvSysEqnRANS<nDim, RANSModel>::viscousFluxThreePoint(
    const MInt orientation, const MFloat A, const MBool isBndry, const MFloat* const surfaceCoords,
    const MFloat* const coord0, const MFloat* const coord1, const MFloat* const cellVars0,
    const MFloat* const cellVars1, const MFloat* const vars0, const MFloat* const vars1, const MFloat* const slope0,
    const MFloat* const slope1, const MFloat f0, const MFloat f1, MFloat* const flux) {
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

template <MInt nDim, class RANSModel>
inline void FvSysEqnRANS<nDim, RANSModel>::computePrimitiveVariables(const MFloat* const cvarsCell,
                                                                     MFloat* const pvarsCell,
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

  // RANS variables
  for(MInt r = 0; r < PV->m_noRansEquations; r++) {
    pvarsCell[PV->NN[r]] = cvarsCell[CV->RHO_NN[r]] * fRho;
  }

  // compute the species
  for(MUint s = 0; s < PV->m_noSpecies; ++s) {
    pvarsCell[PV->Y[s]] = cvarsCell[CV->RHO_Y[s]] * fRho;
  }
}

template <MInt nDim, class RANSModel>
inline void FvSysEqnRANS<nDim, RANSModel>::computeConservativeVariables(const MFloat* const pvarsCell,
                                                                        MFloat* const cvarsCell,
                                                                        const MFloat* const NotUsed(avarsCell)) {
  // compute rho v
  MFloat velPOW2 = F0;
  for(MInt vel = 0; vel < nDim; ++vel) {
    const MFloat v = pvarsCell[vel];
    cvarsCell[vel] = v * pvarsCell[PV->RHO];
    velPOW2 += POW2(v);
  }

  // compute rho e
  cvarsCell[CV->RHO_E] = pvarsCell[PV->P] * m_F1BGammaMinusOne + F1B2 * pvarsCell[PV->RHO] * velPOW2;

  // copy the density
  cvarsCell[CV->RHO] = pvarsCell[PV->RHO];

  for(MInt r = 0; r < PV->m_noRansEquations; r++) {
    cvarsCell[CV->RHO_NN[r]] = pvarsCell[PV->NN[r]] * pvarsCell[PV->RHO];
  }

  // compute rho Yi
  for(MUint s = 0; s < PV->m_noSpecies; ++s) {
    cvarsCell[CV->RHO_Y[s]] = pvarsCell[PV->Y[s]] * pvarsCell[PV->RHO];
  }
}

template <MInt nDim, class RANSModel>
inline std::vector<std::vector<MFloat>>
FvSysEqnRANS<nDim, RANSModel>::conservativeSlopes(const MFloat* const pvarsCell, const MFloat* const NotUsed(cvarsCell),
                                                  const MFloat* const NotUsed(avarsCell),
                                                  const MFloat* const slopesCell) {
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

    for(MInt r = 0; r < m_noRansEquations; r++) {
      dQ[CV->RHO_NN[r]][d] =
          pvarsCell[PV->NN[r]] * slopesCell[PV->RHO * nDim + d] + pvarsCell[PV->RHO] * slopesCell[PV->NN[r] * nDim + d];
    }

    for(MUint s = 0; s < CV->m_noSpecies; ++s) {
      dQ[CV->RHO_Y[s]][d] =
          pvarsCell[PV->Y[s]] * slopesCell[PV->RHO * nDim + d] + pvarsCell[PV->RHO] * slopesCell[PV->Y[s] * nDim + d];
    }
  }

  return dQ;
}

template <MInt nDim, class RANSModel>
void FvSysEqnRANS<nDim, RANSModel>::computeVolumeForcesRANS_SA(const MInt cellId, MFloat* vars, MFloat* rhs,
                                                               const MFloat cellVol, const MFloat* const slopes,
                                                               const MInt recData, const MInt* const /* recNghbr */,
                                                               const MFloat* const /* recConst */,
                                                               const MInt noRecNghbr, MFloat dist) {
  // return;

  std::ignore = recData;
  std::ignore = noRecNghbr;

  const MFloat rRe = F1 / m_Re0;
  const MFloat rho = vars[cellId * PV->noVariables + PV->RHO];
  const MFloat p = vars[cellId * PV->noVariables + PV->P];
  const MFloat nuTilde = vars[cellId * PV->noVariables + PV->N];
  const MFloat T = m_gamma * p / rho;
  const MFloat mue = (T * sqrt(T) * m_sutherlandPlusOne) / (T + m_sutherlandConstant);
  const MFloat nuLaminar = mue / rho;

  MFloat s = 0.0;
  MFloat diff = 0.0;

  // calculation of d(sqrt(rho) nu_tilde)/dx
  // source: "Density Corrections for Turbulence Models", Aerospace Science and Technology, Vol. 4, 2000, pp. 1--11
  for(MInt d = 0; d < nDim; d++) {
    MInt indexRho = cellId * PV->noVariables * nDim + PV->RHO * nDim + d;
    MInt indexNT = cellId * PV->noVariables * nDim + PV->N * nDim + d;
    diff += rRe * RM_SA_DV::Fsigma * RM_SA_DV::cb2
            * POW2((sqrt(rho) * slopes[indexNT] + nuTilde * slopes[indexRho] / (2 * sqrt(rho))));
    // no correction model
    // diff += rRe * RM_SA_DV::Fsigma * (rho * RM_SA_DV::cb2 * POW2(slopes[indexNT]) - (nuTilde + nuLaminar) *
    // slopes[indexNT] * slopes[indexRho]);
  }


  for(MInt d1 = 0; d1 < nDim; d1++) {
    for(MInt d2 = 0; d2 < nDim; d2++) {
      MInt indexVV12 = cellId * PV->noVariables * nDim + PV->VV[d1] * nDim + d2;
      MInt indexVV21 = cellId * PV->noVariables * nDim + PV->VV[d2] * nDim + d1;
      MFloat wij = 0.5 * (slopes[indexVV12] - slopes[indexVV21]);
      s += wij * wij;
    }
  }
  s = sqrt(2 * s);

  const MFloat chi = nuTilde / (nuLaminar);
  const MFloat fv1 = pow(chi, 3) / (pow(chi, 3) + RM_SA_DV::cv1to3);
  const MFloat Fv2 = F1 - (chi / (F1 + chi * fv1));
  const MFloat Fdist2 = 1.0 / (dist * dist);
  const MFloat term = nuTilde * Fdist2 * RM_SA_DV::Fkap2;
  const MFloat stilde = s + term * Fv2 * rRe;
  const MFloat r = mMin(10.0, rRe * term / stilde);
  const MFloat g = r + RM_SA_DV::cw2 * (pow(r, 6) - r);
  const MFloat Fwterm = (1 + RM_SA_DV::cw3to6) / (pow(g, 6) + RM_SA_DV::cw3to6);
  const MFloat Fw = g * pow(Fwterm, (1.0 / 6.0));
  const MFloat Ft2 = RM_SA_DV::Ft2; // RM_SA_DV::ct3 * exp(-RM_SA_DV::ct4*POW2(chi));

  // RM_SA_DV::Ft2 = 0.0 (because simulations with no tripping term)
  MFloat P = (F1 - Ft2) * RM_SA_DV::cb1 * /* fabs( */ nuTilde /* ) */ * mMax(stilde, 0.3 * s);
  MFloat D = rRe * (RM_SA_DV::cw1 * Fw - RM_SA_DV::cb1 * RM_SA_DV::Fkap2 * Ft2) * pow(nuTilde, 2.0) * Fdist2;
  MFloat prodDest = rho * (P - D);

  rhs[CV->RHO_N] -= cellVol * (prodDest + diff);
}

template <MInt nDim, class RANSModel>
void FvSysEqnRANS<nDim, RANSModel>::computeVolumeForcesRANS_FS(const MInt cellId, MFloat* vars, MFloat* rhs,
                                                               const MFloat cellVol, const MFloat* const slopes,
                                                               const MInt recData, const MInt* const recNghbr,
                                                               const MFloat* const recConst, const MInt noRecNghbr,
                                                               MFloat /* dist */) {
  const MFloat rRe = F1 / m_Re0;

  const MFloat rho = vars[cellId * PV->noVariables + PV->RHO];
  const MFloat p = vars[cellId * PV->noVariables + PV->P];
  const MFloat nuTilde = vars[cellId * PV->noVariables + PV->N];
  const MFloat T = m_gamma * p / rho;
  const MFloat mue = (T * sqrt(T) * m_sutherlandPlusOne) / (T + m_sutherlandConstant);
  const MFloat nuLaminar = mue / rho;

  MFloat SijSij = F0;
  MFloat OijOjkSki = F0;
  MFloat diff2 = F0;

  for(MInt d1 = 0; d1 < nDim; d1++) {
    MInt indexN1 = cellId * PV->noVariables * nDim + PV->N * nDim + d1;
    MInt indexRHO1 = cellId * PV->noVariables * nDim + PV->RHO * nDim + d1;
    diff2 += slopes[indexN1] * slopes[indexRHO1];
    for(MInt d2 = 0; d2 < nDim; d2++) {
      MInt indexVV12 = cellId * PV->noVariables * nDim + PV->VV[d1] * nDim + d2;
      MInt indexVV21 = cellId * PV->noVariables * nDim + PV->VV[d2] * nDim + d1;
      MFloat sij = 0.5 * (slopes[indexVV12] + slopes[indexVV21]);
      MFloat wij = 0.5 * (slopes[indexVV12] - slopes[indexVV21]);
      SijSij += sij * sij;
      for(MInt d3 = 0; d3 < nDim; d3++) {
        MInt indexVV13 = cellId * PV->noVariables * nDim + PV->VV[d1] * nDim + d3;
        MInt indexVV31 = cellId * PV->noVariables * nDim + PV->VV[d3] * nDim + d1;
        MInt indexVV23 = cellId * PV->noVariables * nDim + PV->VV[d2] * nDim + d3;
        MInt indexVV32 = cellId * PV->noVariables * nDim + PV->VV[d3] * nDim + d2;
        MFloat ski = 0.5 * (slopes[indexVV13] + slopes[indexVV31]);
        MFloat wjk = 0.5 * (slopes[indexVV23] - slopes[indexVV32]);
        OijOjkSki += wij * wjk * ski;
      }
    }
  }

  diff2 = rRe * (nuLaminar + RM_FS::fasigma * nuTilde) * diff2;

  const MFloat omega = mMax(1e-16, F1 / sqrt(RM_FS::fabetcs) * sqrt(F2 * SijSij));
  MFloat psi_o = fabs(OijOjkSki / (pow(RM_FS::fabetcs * omega, 3.0)));

  MFloat term1 = 0.0;
  MFloat term2 = 0.0;
  MFloat domega[nDim]{};
  MFloat dnutilde[nDim]{};

  for(MInt nghbr = 0; nghbr < noRecNghbr; nghbr++) {
    const MInt nghbrId = recNghbr[nghbr];

    // calculation of omega for nghbr
    MFloat SijSijNgbhr = 0;
    for(MInt d1 = 0; d1 < nDim; d1++) {
      for(MInt d2 = 0; d2 < nDim; d2++) {
        MInt indexVV12 = nghbrId * PV->noVariables * nDim + PV->VV[d1] * nDim + d2;
        MInt indexVV21 = nghbrId * PV->noVariables * nDim + PV->VV[d2] * nDim + d1;
        MFloat sijNgbhr = 0.5 * (slopes[indexVV12] + slopes[indexVV21]);
        SijSijNgbhr += sijNgbhr * sijNgbhr;
      }
    }
    const MFloat omegaNgbhr = 1 / sqrt(RM_FS::fabetcs) * sqrt(2 * SijSijNgbhr);

    const MFloat nutildeNgbhr = vars[nghbrId * PV->noVariables + PV->N];
    for(MInt d = 0; d < nDim; d++) {
      const MFloat recConst_ = recConst[nDim * (recData + nghbr) + d];
      const MFloat deltaOmega = omegaNgbhr - omega;
      const MFloat deltaNutilde = nutildeNgbhr - nuTilde;
      domega[d] += recConst_ * deltaOmega;
      dnutilde[d] += recConst_ * deltaNutilde;
    }
  }

  for(MInt d = 0; d < nDim; ++d) {
    term1 += domega[d] * domega[d];
    term2 += dnutilde[d] * domega[d];
  }

  // psi_o = 0;
  const MFloat psi_k = mMax(F0, F1 / pow(omega, F3) * (nuTilde * term1 + omega * term2));
  const MFloat f_betc = (1 + RM_FS::fapsio1 * psi_o) / (1 + RM_FS::fapsio2 * psi_o);
  const MFloat f_betcs = (1 + RM_FS::fapsik1 * psi_k) / (1 + RM_FS::fapsik2 * psi_k);
  const MFloat beta = f_betc / f_betc * RM_FS::fabetc;
  const MFloat beta_s = f_betcs / f_betcs * RM_FS::fabetcs;

  MFloat Sijduidxj = F0;
  MFloat duidxi = F0;
  for(MInt d1 = 0; d1 < nDim; d1++) {
    for(MInt d2 = 0; d2 < nDim; d2++) {
      MInt indexVV12 = cellId * PV->noVariables * nDim + PV->VV[d1] * nDim + d2;
      MInt indexVV21 = cellId * PV->noVariables * nDim + PV->VV[d2] * nDim + d1;
      Sijduidxj += 0.5 * slopes[indexVV12] * (slopes[indexVV12] + slopes[indexVV21]);
    }
    MInt indexVV11 = cellId * PV->noVariables * nDim + PV->VV[d1] * nDim + d1;
    duidxi += slopes[indexVV11];
  }

  // MFloat P_FS = mMin(0.0, 2 * (F1 - RM_FS::faalpha) * nuTilde * (Sijduidxj / omega - F1B3 * duidxi));
  MFloat P_FS = 2 * (F1 - RM_FS::faalpha) * nuTilde * (Sijduidxj / omega - F1B3 * duidxi);
  MFloat D_FS = (beta_s - beta) * nuTilde * omega;
  MFloat diff1 = mMin(0.0, rRe * F2 * rho * (nuLaminar + RM_FS::fasigma * nuTilde) / omega * term2);

  const MFloat prodDest_FS = rho * (P_FS - D_FS);
  // const MFloat TU = m_turbulenceDegree;
  // const MFloat T = m_TInfinity;
  /* const MFloat mue_init = */
  /*     (m_TInfinity * sqrt(m_TInfinity) * m_sutherlandPlusOne) / (m_TInfinity + m_sutherlandConstant); */
  /* const MFloat nu_init = mue_init / m_rhoInfinity; */
  /* const MFloat k = Re * F3B2 * rho * pow(TU * m_UInfinity, F2) * RM_FS::faphi * pow(nuTilde / nu_init, F5); */

  const MFloat k = 0.0;

  const MFloat diff_FS = diff1 - diff2;

  rhs[CV->RHO_N] -= cellVol * (prodDest_FS + diff_FS - k);
}

template <MInt nDim, class RANSModel>
void FvSysEqnRANS<nDim, RANSModel>::computeVolumeForcesRANS_KOMEGA(const MInt cellId, MFloat* vars, MFloat* rhs,
                                                                   const MFloat cellVol, const MFloat* const slopes,
                                                                   const MInt /* recData */,
                                                                   const MInt* const /* recNghbr */,
                                                                   const MFloat* const /* recConst */,
                                                                   const MInt /* noRecNghbr */, MFloat /* dist */) {
  const MFloat Re = m_Re0;
  const MFloat rRe = F1 / m_Re0;
  const MFloat rho = vars[cellId * PV->noVariables + PV->RHO];     // vars(cellId,PV->RHO);//
  const MFloat k = vars[cellId * PV->noVariables + PV->K];         // vars(cellId,PV->K);//
  const MFloat omega = vars[cellId * PV->noVariables + PV->OMEGA]; // vars(cellId,PV->OMEGA);//

  MFloat dukdxk = F0;
  for(MInt d = 0; d < nDim; d++) {
    dukdxk += slopes[PV->VV[d] * nDim + d];
  }

  MFloat omega_hat = omega; // max(omega,rRe*RM_KOMEGA::C_lim/sqrt(RM_KOMEGA::betas)*sqrt(2*Sij_Sij_));
  // MFloat omega_hat = RM_KOMEGA::C_lim/sqrt(RM_KOMEGA::betas)*sqrt(2*SijSij));

  const MFloat mue_Turb = /*m_Re0 */ rho * k / omega_hat;

  MFloat tau = F0;
  MFloat P = F0;
  MFloat term = F0;
  for(MInt d1 = 0; d1 < nDim; d1++) {
    for(MInt d2 = 0; d2 < nDim; d2++) {
      tau = rRe * mue_Turb * (slopes[PV->VV[d1] * nDim + d2] + slopes[PV->VV[d2] * nDim + d1]);
      if(d1 == d2) {
        tau -= (/*rRe * mue_Turb * F2B3 * dukdxk +*/ F2B3 * rho * k);
      }
      P += slopes[PV->VV[d1] * nDim + d2] * tau;
    }
    term += slopes[PV->K * nDim + d1] * slopes[PV->OMEGA * nDim + d1];
  }

  MFloat OijOjkSki = F0;
  for(MInt d1 = 0; d1 < nDim; d1++) {
    for(MInt d2 = 0; d2 < nDim; d2++) {
      MFloat wij = 0.5 * (slopes[PV->VV[d1] * nDim + d2] - slopes[PV->VV[d2] * nDim + d1]);
      for(MInt d3 = 0; d3 < nDim; d3++) {
        MFloat ski = 0.5 * (slopes[PV->VV[d3] * nDim + d1] + slopes[PV->VV[d1] * nDim + d3]);
        MFloat wjk = 0.5 * (slopes[PV->VV[d2] * nDim + d3] - slopes[PV->VV[d3] * nDim + d2]);
        OijOjkSki += wij * wjk * ski;
      }
    }
  }

  MFloat X_omega = fabs(OijOjkSki / pow(RM_KOMEGA::betas * omega, F3));
  MFloat f_beta = (F1 + RM_KOMEGA::fbeta1 * X_omega) / (F1 + RM_KOMEGA::fbeta2 * X_omega);
  MFloat beta = RM_KOMEGA::beta * f_beta;

  const MFloat P_k = P;
  const MFloat P_o = RM_KOMEGA::alpha * omega / k * P;

  const MFloat D_k = Re * RM_KOMEGA::betas * rho * k * omega;
  const MFloat D_o = Re * beta * rho * POW2(omega);

  // MFloat diff_o = F0;
  // if(term > 0) {
  //   diff_o = rRe * rho * RM_KOMEGA::sigma_d / omega * term;
  // }

  const MFloat prodDest_K = P_k - D_k;
  const MFloat prodDest_OMEGA = P_o - D_o;

  rhs[CV->RHO_K] -= cellVol * (prodDest_K);
  rhs[CV->RHO_OMEGA] -= cellVol * (prodDest_OMEGA);
}

#endif // FVSYSEQNRANS_H_
