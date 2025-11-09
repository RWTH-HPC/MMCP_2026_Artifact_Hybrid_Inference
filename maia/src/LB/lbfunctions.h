// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef LBFUNCTIONS_H
#define LBFUNCTIONS_H
#include <array>
#include <numeric>
#include "lbconstants.h"
#include "lblatticedescriptor.h"

// Within this namespace inline functions are provided for several LB
// calculations given as templates
// Definitions must be independent from anything besides lbconstants.h.
namespace lbfunc {
using namespace lbconstants;

/** \brief  Calculate equilibrium function (incompressible version). Squared velocity has already been computed.
 *  \author (?, refactored by Miro Gondrum)
 *  \date   14.01.2020
 *  \param[in]  p_rho               macroscopic rho
 *  \param[in]  p_squaredVelocity   macroscopic velocity squared
 *  \param[in]  p_u                 macroscopic velocity
 *  \param[out] eqDist              equilibrium distributions based on given macrosopic variables
 */
template <MInt nDim, MInt nDist>
inline void calcEqDistsIncompressible(const MFloat& p_rho, const MFloat& p_squaredVelocity, MFloat const* const p_u,
                                      MFloat* const eqDist) {
  using Ld = LbLatticeDescriptor<nDim, nDist>;
  const MFloat squaredVelocityB2 = p_squaredVelocity * F1B2;

  const MFloat rhomulCSsq = p_rho * CSsq;

  const MFloat preTerm = rhomulCSsq - squaredVelocityB2;

  MFloat b[2 * nDim];
  for(MInt n = 0; n < nDim; n++) {
    b[2 * n] = -p_u[n];
    b[2 * n + 1] = p_u[n];
  }

  // Calculation of distributions for directions with only one component
  constexpr MFloat lb_tp_coef1 = Ld::tp(1) * F1BCSsq;
  for(MInt j = 0; j < Ld::distFld(0); j++) {
    eqDist[j] = lb_tp_coef1 * (preTerm + b[j] + b[j] * b[j] * F1B2mulF1BCSsq);
  }

  // Calculation of distributions for directions with two components
  constexpr MFloat lb_tp_coef2 = Ld::tp(2) * F1BCSsq;
  for(MInt j = 0; j < Ld::distFld(1); j++) {
    const MFloat tmp = (b[Ld::mFld1(2 * j)] + b[Ld::mFld1(2 * j + 1)]);
    eqDist[Ld::distFld(0) + j] = lb_tp_coef2 * (preTerm + tmp + tmp * tmp * F1B2mulF1BCSsq);
  }

  // Calculation of distributions for directions with three components
  constexpr MFloat lb_tp_coef3 = Ld::tp(3) * F1BCSsq;
  for(MInt j = 0; j < Ld::distFld(2); j++) {
    const MFloat tmp = (b[Ld::mFld2(3 * j)] + b[Ld::mFld2(3 * j + 1)] + b[Ld::mFld2(3 * j + 2)]);
    eqDist[Ld::distFld(0) + Ld::distFld(1) + j] = lb_tp_coef3 * (preTerm + tmp + tmp * tmp * F1B2mulF1BCSsq);
  }

  // Calculation of distribution for rest particle distribution (center)
  eqDist[Ld::lastId()] = Ld::tp(0) * (p_rho - F1B2mulF1BCSsq * p_squaredVelocity);
}

/** \brief  Calculate equilibrium function (compressible version). Squared velocity has already been computed.
 *  \author (?, refactored by Miro Gondrum, compressible changes by Daniel Lauwers)
 *  \date   14.01.2020
 *  \param[in]  p_rho     macroscopic rho
 *  \param[in]  p_squaredVelocity  macroscopic velocity squared
 *  \param[in]  p_u       macroscopic velocity
 *  \param[out] eqDist    equilibrium distributions based on given macrosopic variables
 */
template <MInt nDim, MInt nDist>
inline void calcEqDistsCompressible(const MFloat& p_rho, const MFloat& p_squaredVelocity, MFloat const* const p_u,
                                    MFloat* const eqDist) {
  using Ld = LbLatticeDescriptor<nDim, nDist>;
  const MFloat squaredVelocityB2 = p_squaredVelocity * F1B2;

  MFloat b[2 * nDim];
  for(MInt n = 0; n < nDim; n++) {
    b[2 * n] = -p_u[n];
    b[2 * n + 1] = p_u[n];
  }

  // Calculation of distributions for directions with only one component
  constexpr MFloat lb_tp_coef1 = Ld::tp(1) * F1BCSsq;
  for(MInt j = 0; j < Ld::distFld(0); j++) {
    eqDist[j] = lb_tp_coef1 * (CSsq + b[j] + b[j] * b[j] * F1B2mulF1BCSsq - squaredVelocityB2) * p_rho;
  }

  // Calculation of distributions for directions with two components
  constexpr MFloat lb_tp_coef2 = Ld::tp(2) * F1BCSsq;
  for(MInt j = 0; j < Ld::distFld(1); j++) {
    const MFloat tmp = (b[Ld::mFld1(2 * j)] + b[Ld::mFld1(2 * j + 1)]);
    eqDist[Ld::distFld(0) + j] = lb_tp_coef2 * (CSsq + tmp + tmp * tmp * F1B2mulF1BCSsq - squaredVelocityB2) * p_rho;
  }

  // Calculation of distributions for directions with three components
  constexpr MFloat lb_tp_coef3 = Ld::tp(3) * F1BCSsq;
  for(MInt j = 0; j < Ld::distFld(2); j++) {
    const MFloat tmp = (b[Ld::mFld2(3 * j)] + b[Ld::mFld2(3 * j + 1)] + b[Ld::mFld2(3 * j + 2)]);
    eqDist[Ld::distFld(0) + Ld::distFld(1) + j] =
        lb_tp_coef3 * (CSsq + tmp + tmp * tmp * F1B2mulF1BCSsq - squaredVelocityB2) * p_rho;
  }

  // Calculation of distribution for rest particle distribution (center)
  eqDist[Ld::lastId()] = Ld::tp(0) * (1.0 - F1B2mulF1BCSsq * p_squaredVelocity) * p_rho;
}

/** \brief  Calculate equilibrium function. Squared velocity has already been computed.
 *  \author Daniel Lauwers
 *  \date   01.12.2021
 *  \param[in]  p_rho               macroscopic rho
 *  \param[in]  p_squaredVelocity   macroscopic velocity squared
 *  \param[in]  p_u                 macroscopic velocity
 *  \param[out] eqDist              equilibrium distributions based on given macrosopic variables
 */
template <MInt nDim, MInt nDist, MBool compressible = false>
inline void calcEqDists(const MFloat& p_rho, const MFloat& p_squaredVelocity, MFloat const* const p_u,
                        MFloat* const eqDist) {
  if constexpr(compressible)
    calcEqDistsCompressible<nDim, nDist>(p_rho, p_squaredVelocity, p_u, eqDist);
  else
    calcEqDistsIncompressible<nDim, nDist>(p_rho, p_squaredVelocity, p_u, eqDist);
}

/** \brief  Calculate equilibrium function
 *  \author (?, refactored by Miro Gondrum)
 *  \date   14.01.2020
 *  \param[in]  p_rho   macroscopic rho
 *  \param[in]  p_u     macroscopic velocity
 *  \param[out] eqDist  equilibrium distributions based on given macrosopic variables
 */
template <MInt nDim, MInt nDist, MBool compressible = false>
inline void calcEqDists(const MFloat& p_rho, MFloat const* const p_u, MFloat* const eqDist) {
  const MFloat p_squaredVelocity = std::inner_product(&p_u[0], &p_u[nDim], &p_u[0], .0);

  calcEqDists<nDim, nDist, compressible>(p_rho, p_squaredVelocity, p_u, eqDist);
}

#ifdef WAR_NVHPC_PSTL
/** \brief  Calculate equilibrium function (incompressible version, GPU version). Squared velocity has already been
 * computed. \author (?, refactored by Moritz Waldmann) \date   14.01.2020 \param[in]  p_rho       macroscopic rho
 *  \param[in]  p_squaredVelocity   macroscopic velocity squared
 *  \param[in]  p_u                 macroscopic velocity
 *  \param[out] eqDist              equilibrium distributions based on given macrosopic variables
 */
template <MInt nDim, MInt nDist>
inline void calcEqDistsIncompressible(const MFloat& p_rho, const MFloat& p_squaredVelocity, MFloat const* const p_u,
                                      MFloat* const eqDist, const MInt* myMFld1, const MInt* myMFld2,
                                      const MFloat* myTp, const MInt* myDistFld) {
  using Ld = LbLatticeDescriptor<nDim, nDist>;
  const MFloat squaredVelocityB2 = p_squaredVelocity * F1B2;

  const MFloat rhomulCSsq = p_rho * CSsq;

  MFloat b[2 * nDim];
  for(MInt n = 0; n < nDim; n++) {
    b[2 * n] = -p_u[n];
    b[2 * n + 1] = p_u[n];
  }

  // Calculation of distributions for directions with only one component
  const MFloat lb_tp_coef1 = myTp[1] * F1BCSsq;
  for(MInt j = 0; j < myDistFld[0]; j++) {
    eqDist[j] = lb_tp_coef1 * (rhomulCSsq + b[j] + b[j] * b[j] * F1B2mulF1BCSsq - squaredVelocityB2);
  }

  // Calculation of distributions for directions with two components
  const MFloat lb_tp_coef2 = myTp[2] * F1BCSsq;
  for(MInt j = 0; j < myDistFld[1]; j++) {
    const MFloat tmp = (b[myMFld1[2 * j]] + b[myMFld1[2 * j + 1]]);
    eqDist[myDistFld[0] + j] = lb_tp_coef2 * (rhomulCSsq + tmp + tmp * tmp * F1B2mulF1BCSsq - squaredVelocityB2);
  }

  // Calculation of distributions for directions with three components
  const MFloat lb_tp_coef3 = myTp[3] * F1BCSsq;
  for(MInt j = 0; j < myDistFld[2]; j++) {
    const MFloat tmp = (b[myMFld2[3 * j]] + b[myMFld2[3 * j + 1]] + b[myMFld2[3 * j + 2]]);
    eqDist[myDistFld[0] + myDistFld[1] + j] =
        lb_tp_coef3 * (rhomulCSsq + tmp + tmp * tmp * F1B2mulF1BCSsq - squaredVelocityB2);
  }

  // Calculation of distribution for rest particle distribution (center)
  eqDist[Ld::lastId()] = myTp[0] * (p_rho - F1B2mulF1BCSsq * p_squaredVelocity);
}

/** \brief  Calculate equilibrium function (compressible version, GPU version). Squared velocity has already been
 * computed. \author (?, refactored by Moritz Waldmann) \date   14.01.2020 \param[in]  p_rho       macroscopic rho
 *  \param[in]  p_squaredVelocity   macroscopic velocity squared
 *  \param[in]  p_u                 macroscopic velocity
 *  \param[out] eqDist              equilibrium distributions based on given macrosopic variables
 */
template <MInt nDim, MInt nDist>
inline void calcEqDistsCompressible(const MFloat& p_rho, const MFloat& p_squaredVelocity, MFloat const* const p_u,
                                    MFloat* const eqDist, const MInt* myMFld1, const MInt* myMFld2, const MFloat* myTp,
                                    const MInt* myDistFld) {
  using Ld = LbLatticeDescriptor<nDim, nDist>;
  const MFloat squaredVelocityB2 = p_squaredVelocity * F1B2;

  MFloat b[2 * nDim];
  for(MInt n = 0; n < nDim; n++) {
    b[2 * n] = -p_u[n];
    b[2 * n + 1] = p_u[n];
  }

  // Calculation of distributions for directions with only one component
  const MFloat lb_tp_coef1 = myTp[1] * F1BCSsq;
  for(MInt j = 0; j < myDistFld[0]; j++) {
    eqDist[j] = lb_tp_coef1 * (CSsq + b[j] + b[j] * b[j] * F1B2mulF1BCSsq - squaredVelocityB2) * p_rho;
  }

  // Calculation of distributions for directions with two components
  const MFloat lb_tp_coef2 = myTp[2] * F1BCSsq;
  for(MInt j = 0; j < myDistFld[1]; j++) {
    const MFloat tmp = (b[myMFld1[2 * j]] + b[myMFld1[2 * j + 1]]);
    eqDist[myDistFld[0] + j] = lb_tp_coef2 * (CSsq + tmp + tmp * tmp * F1B2mulF1BCSsq - squaredVelocityB2) * p_rho;
  }

  // Calculation of distributions for directions with three components
  const MFloat lb_tp_coef3 = myTp[3] * F1BCSsq;
  for(MInt j = 0; j < myDistFld[2]; j++) {
    const MFloat tmp = (b[myMFld2[3 * j]] + b[myMFld2[3 * j + 1]] + b[myMFld2[3 * j + 2]]);
    eqDist[myDistFld[0] + myDistFld[1] + j] =
        lb_tp_coef3 * (CSsq + tmp + tmp * tmp * F1B2mulF1BCSsq - squaredVelocityB2) * p_rho;
  }

  // Calculation of distribution for rest particle distribution (center)
  eqDist[Ld::lastId()] = myTp[0] * (1.0 - F1B2mulF1BCSsq * p_squaredVelocity) * p_rho;
}

/** \brief  Calculate equilibrium function. Squared velocity has already been computed.
 *  \author Daniel Lauwers
 *  \date   01.12.2021
 *  \param[in]  p_rho               macroscopic rho
 *  \param[in]  p_squaredVelocity   macroscopic velocity squared
 *  \param[in]  p_u                 macroscopic velocity
 *  \param[out] eqDist              equilibrium distributions based on given macrosopic variables
 */
template <MInt nDim, MInt nDist, MBool compressible = false>
inline void calcEqDists(const MFloat& p_rho, const MFloat& p_squaredVelocity, MFloat const* const p_u,
                        MFloat* const eqDist, const MInt* myMFld1, const MInt* myMFld2, const MFloat* myTp,
                        const MInt* myDistFld) {
  if constexpr(compressible)
    calcEqDistsCompressible<nDim, nDist>(p_rho, p_squaredVelocity, p_u, eqDist, myMFld1, myMFld2, myTp, myDistFld);
  else
    calcEqDistsIncompressible<nDim, nDist>(p_rho, p_squaredVelocity, p_u, eqDist, myMFld1, myMFld2, myTp, myDistFld);
}

/** \brief  Calculate equilibrium function
 *  \author (?, refactored by Miro Gondrum)
 *  \date   14.01.2020
 *  \param[in]  p_rho   macroscopic rho
 *  \param[in]  p_u     macroscopic velocity
 *  \param[out] eqDist  equilibrium distributions based on given macrosopic variables
 */
template <MInt nDim, MInt nDist, MBool compressible = false>
inline void calcEqDists(const MFloat& p_rho, MFloat const* const p_u, MFloat* const eqDist, const MInt* myMFld1,
                        const MInt* myMFld2, const MFloat* myTp, const MInt* myDistFld) {
  const MFloat p_squaredVelocity = std::inner_product(&p_u[0], &p_u[nDim], &p_u[0], .0);

  calcEqDists<nDim, nDist, compressible>(p_rho, p_squaredVelocity, p_u, eqDist, myMFld1, myMFld2, myTp, myDistFld);
}
#endif

/**
 * \brief Calculate the non-equilibrium part for given macroscopic variables and distribution
 *
 * Remember: This calculation is only meaningful for the post-propagation state!
 *
 * \author Julian Vorspohl <j.vorspohl@aia.rwth-aachen.de>
 *
 * \param[in]  rho       Macroscopic density
 * \param[in]  u         Macroscopic velocity
 * \param[in]  dist      Distribution
 * \param[out] nonEqDist Non-equilibrium part of the distribution
 */
template <MInt nDim, MInt nDist, MBool compressible = false>
inline void calcNonEqDists(const MFloat& rho, const MFloat* const u, const MFloat* const dist,
                           MFloat* const nonEqDist) {
  // Calculate equlibirum distribution
  std::array<MFloat, nDist> eqDist{};
  calcEqDists<nDim, nDist, compressible>(rho, u, eqDist.data());

  // Calculate non-equilibium part
  for(MInt j = 0; j < nDist; j++) {
    nonEqDist[j] = dist[j] - eqDist[j];
  }
}

/**
 * \brief Calculate the non-equilibrium part for given distribution
 *
 * Remember: This calculation is only meaningful for the post-propagation state!
 *
 * \author Julian Vorspohl <j.vorspohl@aia.rwth-aachen.de>
 *
 * \param[in]  dist      Distribution
 * \param[out] nonEqDist Non-equilibrium part of the distribution
 */
template <MInt nDim, MInt nDist, MBool compressible = false>
inline void calcNonEqDists(const MFloat* const dist, MFloat* const nonEqDist) {
  // Calculate macroscopic variables
  MFloat rho{};
  std::array<MFloat, nDim> u{};
  calcMacroscopicVars(dist, rho, u.data());

  calcNonEqDists<nDim, nDist, compressible>(rho, u.data(), dist, nonEqDist);
}

/**
 * \brief Calculate the momentum flux for a given non-equilibirum distribution
 *
 * Remember: This calculation is only meaningful for the post-propagation state!
 *
 * \author Julian Vorspohl <j.vorspohl@aia.rwth-aachen.de>
 *
 * \param[in]  nonEqDist    Non-Equilibrium distribution
 * \param[out] momentumFlux Momentum flux tensor
 */
template <MInt nDim, MInt nDist>
inline void calcMomentumFluxFromNonEq(const MFloat* const nonEqDist, MFloat* const momentumFlux) {
  using Ld = LbLatticeDescriptor<nDim, nDist>;

  for(MInt d = 0; d < nDim * nDim; d++) {
    momentumFlux[d] = F0;
  }
  for(MInt j = 0; j < nDist; j++) {
    for(MInt k = 0; k < nDim; k++) {
      for(MInt l = 0; l < nDim; l++) {
        momentumFlux[k * nDim + l] += nonEqDist[j] * Ld::ppdfDir(j, k) * Ld::ppdfDir(j, l);
      }
    }
  }
}

/**
 * \brief Calculate the momentum flux for given macroscopic variables and distribution
 *
 * Remember: This calculation is only meaningful for the post-propagation state!
 *
 * \author Julian Vorspohl <j.vorspohl@aia.rwth-aachen.de>
 *
 * \param[in]  rho          Macroscopic density
 * \param[in]  u            Macroscopic velocity
 * \param[in]  dist         Distribution
 * \param[out] momentumFlux Momentum flux tensor
 */
template <MInt nDim, MInt nDist, MBool compressible = false>
inline void calcMomentumFlux(const MFloat& rho, const MFloat* const u, const MFloat* const dist,
                             MFloat* const momentumFlux) {
  std::array<MFloat, nDist> nonEqDist{};
  calcNonEqDists<nDim, nDist, compressible>(rho, u, dist, nonEqDist.data());

  calcMomentumFluxFromNonEq<nDim, nDist>(nonEqDist.data(), momentumFlux);
}

/**
 * \brief Calculate the momentum flux for given macroscopic variables and distribution
 *
 * Remember: This calculation is only meaningful for the post-propagation state!
 *
 * \author Julian Vorspohl <j.vorspohl@aia.rwth-aachen.de>
 *
 * \param[in]  dist         Distribution
 * \param[out] momentumFlux Momentum flux tensor
 */
template <MInt nDim, MInt nDist, MBool compressible = false>
inline void calcMomentumFlux(const MFloat* const dist, MFloat* const momentumFlux) {
  std::array<MFloat, nDist> nonEqDist{};
  calcNonEqDists<nDim, nDist, compressible>(dist, nonEqDist.data());

  calcMomentumFluxFromNonEq<nDim, nDist>(nonEqDist.data(), momentumFlux);
}

/**
 * \brief Calculate thermal equilibrium distribution
 *
 * Attention: This is the basic thermal LBM
 *
 * Use this version if the squaredVelocity is already precomputed
 *
 * \author Julian Vorspohl <j.vorspohl@aia.rwth-aachen.de>
 * \date 29.09.2021
 *
 * \param[in]  T      Macroscopic temperature
 * \param[out] eqDist Resulting equilibrium distribution
 */
template <MInt nDim, MInt nDist>
inline void calcEqDistsThermal(const MFloat& T, const MFloat squaredVelocity, const MFloat* const u,
                               MFloat* const eqDist) {
  calcEqDists<nDim, nDist, true>(T, squaredVelocity, u, eqDist);
}

/**
 * \brief Calculate thermal equilibrium distribution
 *
 * Attention: This is the basic thermal LBM
 *
 * \author Julian Vorspohl <j.vorspohl@aia.rwth-aachen.de>
 * \date 29.09.2021
 *
 * \param[in]  T      Macroscopic temperature
 * \param[out] eqDist Resulting equilibrium distribution
 */
template <MInt nDim, MInt nDist>
inline void calcEqDistsThermal(const MFloat& T, const MFloat* const u, MFloat* const eqDist) {
  const MFloat squaredVelocity = std::inner_product(&u[0], &u[nDim], &u[0], .0);

  calcEqDists<nDim, nDist, true>(T, squaredVelocity, u, eqDist);
}

#ifdef WAR_NVHPC_PSTL
/**
 * \brief Calculate thermal equilibrium distribution
 *
 * Attention: This is the basic thermal LBM
 *
 * Use this version if the squaredVelocity is already precomputed
 *
 * \author Julian Vorspohl <j.vorspohl@aia.rwth-aachen.de>
 * \date 29.09.2021
 *
 * \param[in]  T      Macroscopic temperature
 * \param[out] eqDist Resulting equilibrium distribution
 */
template <MInt nDim, MInt nDist>
inline void calcEqDistsThermal(const MFloat& T, const MFloat squaredVelocity, const MFloat* const u,
                               MFloat* const eqDist, const MInt* myMFld1, const MInt* myMFld2, const MFloat* myTp,
                               const MInt* myDistFld) {
  calcEqDists<nDim, nDist, true>(T, squaredVelocity, u, eqDist, myMFld1, myMFld2, myTp, myDistFld);
}

/**
 * \brief Calculate thermal equilibrium distribution
 *
 * Attention: This is the basic thermal LBM
 *
 * \author Julian Vorspohl <j.vorspohl@aia.rwth-aachen.de>
 * \date 29.09.2021
 *
 * \param[in]  T      Macroscopic temperature
 * \param[out] eqDist Resulting equilibrium distribution
 */
template <MInt nDim, MInt nDist>
inline void calcEqDistsThermal(const MFloat& T, const MFloat* const u, MFloat* const eqDist, const MInt* myMFld1,
                               const MInt* myMFld2, const MFloat* myTp, const MInt* myDistFld) {
  const MFloat squaredVelocity = std::inner_product(&u[0], &u[nDim], &u[0], .0);

  calcEqDists<nDim, nDist, true>(T, squaredVelocity, u, eqDist, myMFld1, myMFld2, myTp, myDistFld);
}
#endif

template <MInt nDim, MInt nDist>
inline void calcEqDistsInnerEnergy(const MFloat& T, const MFloat& rho, const MFloat squaredVelocity,
                                   const MFloat* const u, MFloat* const eqDist) {
  using Ld = LbLatticeDescriptor<nDim, nDist>;

  const MFloat F1BD = F1 / nDim;
  const MFloat squaredVelocityB2 = 0.5 * squaredVelocity;
  const MFloat innerEnergy = rho * T * F1B2 * (MFloat)nDim;

  std::array<MFloat, 2 * nDim> b{};
  for(MInt d = 0; d < nDim; d++) {
    b[2 * d] = -u[d];
    b[2 * d + 1] = u[d];
  }

  // Calculate all equilibrium distributions
  // Calculation of eq. distributions for directions with only one component
  for(MInt j = 0; j < Ld::distFld(0); j++) {
    const MFloat l_innerterm_ie =
        F1BCSsq * F1BD * CSsq + (F1BCSsq * F1BD - F2 * F1BD) * b[j] + b[j] * b[j] * F1B2 * F1BCSsq - squaredVelocityB2;
    eqDist[j] = Ld::tp(1) * F1BCSsq * innerEnergy * l_innerterm_ie;
  }

  // Calculation of eq. distributions for directions with two components
  for(MInt j = 0; j < Ld::distFld(1); j++) {
    const MFloat tmp = (b[Ld::mFld1(2 * j)] + b[Ld::mFld1(2 * j + 1)]);
    const MFloat l_innerterm_ie = F2 * F1BCSsq * F1BD * CSsq + (F2 * F1BCSsq * F1BD - F2 * F1BD) * tmp
                                  + tmp * tmp * F1B2 * F1BCSsq - squaredVelocityB2;
    eqDist[Ld::distFld(0) + j] = Ld::tp(2) * F1BCSsq * innerEnergy * l_innerterm_ie;
  }

  // Calculation of eq. distributions for directions with three components
  for(MInt j = 0; j < Ld::distFld(2); j++) {
    const MFloat tmp = (b[Ld::mFld2(3 * j)] + b[Ld::mFld2(3 * j + 1)] + b[Ld::mFld2(3 * j + 2)]);
    const MFloat l_innerterm_ie = F3 * F1BCSsq * F1BD * CSsq + (F3 * F1BCSsq * F1BD - F2 * F1BD) * tmp
                                  + tmp * tmp * F1B2 * F1BCSsq - squaredVelocityB2;
    eqDist[Ld::distFld(0) + Ld::distFld(1) + j] = Ld::tp(3) * F1BCSsq * innerEnergy * l_innerterm_ie;
  }

  // Rest distribution
  eqDist[Ld::lastId()] = -Ld::tp(0) * innerEnergy * F1B2 * F1BCSsq * squaredVelocity;
}

template <MInt nDim, MInt nDist>
inline void calcEqDistsInnerEnergy(const MFloat& T, const MFloat& rho, const MFloat* const u, MFloat* const eqDist) {
  const MFloat squaredVelocity = std::inner_product(&u[0], &u[nDim], &u[0], .0);

  calcEqDistsInnerEnergy<nDim, nDist>(T, rho, squaredVelocity, u, eqDist);
}

#ifdef WAR_NVHPC_PSTL
template <MInt nDim, MInt nDist>
inline void calcEqDistsInnerEnergy(const MFloat& T, const MFloat& rho, const MFloat squaredVelocity,
                                   const MFloat* const u, MFloat* const eqDist, const MInt* myMFld1,
                                   const MInt* myMFld2, const MFloat* myTp, const MInt* myDistFld) {
  using Ld = LbLatticeDescriptor<nDim, nDist>;

  const MFloat F1BD = F1 / nDim;
  const MFloat squaredVelocityB2 = 0.5 * squaredVelocity;
  const MFloat innerEnergy = rho * T * F1B2 * nDim;

  std::array<MFloat, 2 * nDim> b{};
  for(MInt d = 0; d < nDim; d++) {
    b[2 * d] = -u[d];
    b[2 * d + 1] = u[d];
  }

  // Calculate all equilibrium distributions
  // Calculation of eq. distributions for directions with only one component
  const MFloat lb_tp_coef1 = myTp[1] * F1BCSsq;
  for(MInt j = 0; j < myDistFld[0]; j++) {
    const MFloat l_innerterm_ie =
        F1BCSsq * F1BD * CSsq + (F1BCSsq * F1BD - F2 * F1BD) * b[j] + b[j] * b[j] * F1B2 * F1BCSsq - squaredVelocityB2;
    eqDist[j] = lb_tp_coef1 * innerEnergy * l_innerterm_ie;
  }

  // Calculation of eq. distributions for directions with two components
  const MFloat lb_tp_coef2 = myTp[2] * F1BCSsq;
  for(MInt j = 0; j < myDistFld[1]; j++) {
    const MFloat tmp = (b[myMFld1[2 * j]] + b[myMFld1[2 * j + 1]]);
    const MFloat l_innerterm_ie = F2 * F1BCSsq * F1BD * CSsq + (F2 * F1BCSsq * F1BD - F2 * F1BD) * tmp
                                  + tmp * tmp * F1B2 * F1BCSsq - squaredVelocityB2;
    eqDist[myDistFld[0] + j] = lb_tp_coef2 * innerEnergy * l_innerterm_ie;
  }

  // Calculation of eq. distributions for directions with three components
  const MFloat lb_tp_coef3 = myTp[3] * F1BCSsq;
  for(MInt j = 0; j < myDistFld[2]; j++) {
    const MFloat tmp = (b[myMFld2[3 * j]] + b[myMFld2[3 * j + 1]] + b[myMFld2[3 * j + 2]]);
    const MFloat l_innerterm_ie = F3 * F1BCSsq * F1BD * CSsq + (F3 * F1BCSsq * F1BD - F2 * F1BD) * tmp
                                  + tmp * tmp * F1B2 * F1BCSsq - squaredVelocityB2;
    eqDist[myDistFld[0] + myDistFld[1] + j] = lb_tp_coef3 * innerEnergy * l_innerterm_ie;
  }

  // Rest distribution
  eqDist[Ld::lastId()] = -myTp[0] * innerEnergy * F1B2 * F1BCSsq * squaredVelocity;
}

template <MInt nDim, MInt nDist>
inline void calcEqDistsInnerEnergy(const MFloat& T, const MFloat& rho, const MFloat* const u, MFloat* const eqDist,
                                   const MInt* myMFld1, const MInt* myMFld2, const MFloat* myTp,
                                   const MInt* myDistFld) {
  const MFloat squaredVelocity = std::inner_product(&u[0], &u[nDim], &u[0], .0);

  calcEqDistsInnerEnergy<nDim, nDist>(T, rho, squaredVelocity, u, eqDist, myMFld1, myMFld2, myTp, myDistFld);
}
#endif

template <MInt nDim, MInt nDist>
inline void calcEqDistsTotalEnergy(const MFloat& T, const MFloat& rho, const MFloat squaredVelocity,
                                   const MFloat* const u, MFloat* const eqDistT) {
  using Ld = LbLatticeDescriptor<nDim, nDist>;

  const MFloat squaredVelocityB2 = 0.5 * squaredVelocity;
  const MFloat p0 = rho * CSsq;
  const MFloat totalEnergy = T * nDim * F1B2 + squaredVelocityB2;

  std::array<MFloat, 2 * nDim> b{};
  for(MInt d = 0; d < nDim; d++) {
    b[d * 2] = -u[d];
    b[d * 2 + 1] = u[d];
  }

  // Calculate all equilibrium distributions
  // Calculation of eq. distributions for directions with only one component
  constexpr MFloat lb_tp_coef1 = Ld::tp(1) * F1BCSsq;
  for(MInt j = 0; j < Ld::distFld(0); j++) {
    const MFloat l_innerterm = CSsq + b[j] + b[j] * b[j] * F1B2mulF1BCSsq - squaredVelocityB2;
    const MFloat eq_dist = lb_tp_coef1 * rho * l_innerterm;
    const MFloat l_innerterm_te = b[j] + b[j] * b[j] * F1BCSsq - squaredVelocity + F1B2 * (F1 - nDim * CSsq);
    eqDistT[j] = lb_tp_coef1 * p0 * l_innerterm_te + totalEnergy * eq_dist;
  }

  // Calculation of eq. distributions for directions with two components
  constexpr MFloat lb_tp_coef2 = Ld::tp(2) * F1BCSsq;
  for(MInt j = 0; j < Ld::distFld(1); j++) {
    const MFloat tmp = (b[Ld::mFld1(2 * j)] + b[Ld::mFld1(2 * j + 1)]);
    const MFloat l_innerterm = CSsq + tmp + tmp * tmp * F1B2mulF1BCSsq - squaredVelocityB2;
    const MFloat eq_dist = lb_tp_coef2 * rho * l_innerterm;
    const MFloat l_innerterm_te = tmp + tmp * tmp * F1BCSsq - squaredVelocity + F1B2 * (F2 - nDim * CSsq);
    eqDistT[Ld::distFld(0) + j] = lb_tp_coef2 * p0 * l_innerterm_te + totalEnergy * eq_dist;
  }

  // Calculation of eq. distributions for directions with three components
  constexpr MFloat lb_tp_coef3 = Ld::tp(3) * F1BCSsq;
  for(MInt j = 0; j < Ld::distFld(2); j++) {
    const MFloat tmp = (b[Ld::mFld2(3 * j)] + b[Ld::mFld2(3 * j + 1)] + b[Ld::mFld2(3 * j + 2)]);
    const MFloat l_innerterm = CSsq + tmp + tmp * tmp * F1B2mulF1BCSsq - squaredVelocityB2;
    const MFloat eq_dist = lb_tp_coef3 * rho * l_innerterm;
    const MFloat l_innerterm_te = tmp + tmp * tmp * F1BCSsq - squaredVelocity + F1B2 * (F3 - nDim * CSsq);
    eqDistT[Ld::distFld(0) + Ld::distFld(1) + j] = lb_tp_coef3 * p0 * l_innerterm_te + totalEnergy * eq_dist;
  }

  // Rest distribution
  constexpr MFloat lb_tp_coef0 = Ld::tp(0);
  const MFloat l_innerterm = F1 - F1B2mulF1BCSsq * squaredVelocity;
  const MFloat eq_dist = lb_tp_coef0 * rho * l_innerterm;
  const MFloat l_innerterm_te = -F1BCSsq * squaredVelocity + F1B2mulF1BCSsq * (F0 - nDim * CSsq);
  eqDistT[Ld::lastId()] = lb_tp_coef0 * p0 * l_innerterm_te + totalEnergy * eq_dist;
}

template <MInt nDim, MInt nDist>
inline void calcEqDistsTotalEnergy(const MFloat& T, const MFloat& rho, const MFloat* const u, MFloat* const eqDist) {
  const MFloat squaredVelocity = std::inner_product(&u[0], &u[nDim], &u[0], .0);

  calcEqDistsTotalEnergy<nDim, nDist>(T, rho, squaredVelocity, u, eqDist);
}

#ifdef WAR_NVHPC_PSTL
template <MInt nDim, MInt nDist>
inline void calcEqDistsTotalEnergy(const MFloat& T, const MFloat& rho, const MFloat squaredVelocity,
                                   const MFloat* const u, MFloat* const eqDistT, const MInt* myMFld1,
                                   const MInt* myMFld2, const MFloat* myTp, const MInt* myDistFld) {
  using Ld = LbLatticeDescriptor<nDim, nDist>;

  const MFloat squaredVelocityB2 = 0.5 * squaredVelocity;
  const MFloat p0 = rho * CSsq;
  const MFloat totalEnergy = T * nDim * F1B2 + squaredVelocityB2;

  std::array<MFloat, 2 * nDim> b{};
  for(MInt d = 0; d < nDim; d++) {
    b[2 * d] = -u[d];
    b[2 * d + 1] = u[d];
  }

  // Calculate all equilibrium distributions
  // Calculation of eq. distributions for directions with only one component
  const MFloat lb_tp_coef1 = myTp[1] * F1BCSsq;
  for(MInt j = 0; j < myDistFld[0]; j++) {
    const MFloat l_innerterm = CSsq + b[j] + b[j] * b[j] * F1B2mulF1BCSsq - squaredVelocityB2;
    const MFloat eq_dist = lb_tp_coef1 * rho * l_innerterm;
    const MFloat l_innerterm_te = b[j] + b[j] * b[j] * F1BCSsq - squaredVelocity + F1B2 * (F1 - nDim * CSsq);
    eqDistT[j] = lb_tp_coef1 * p0 * l_innerterm_te + totalEnergy * eq_dist;
  }

  // Calculation of eq. distributions for directions with two components
  const MFloat lb_tp_coef2 = myTp[2] * F1BCSsq;
  for(MInt j = 0; j < myDistFld[1]; j++) {
    const MFloat tmp = (b[myMFld1[2 * j]] + b[myMFld1[2 * j + 1]]);
    const MFloat l_innerterm = CSsq + tmp + tmp * tmp * F1B2mulF1BCSsq - squaredVelocityB2;
    const MFloat eq_dist = lb_tp_coef2 * rho * l_innerterm;
    const MFloat l_innerterm_te = tmp + tmp * tmp * F1BCSsq - squaredVelocity + F1B2 * (F2 - nDim * CSsq);
    eqDistT[myDistFld[0] + j] = lb_tp_coef2 * p0 * l_innerterm_te + totalEnergy * eq_dist;
  }

  // Calculation of eq. distributions for directions with three components
  const MFloat lb_tp_coef3 = myTp[3] * F1BCSsq;
  for(MInt j = 0; j < myDistFld[2]; j++) {
    const MFloat tmp = (b[myMFld2[3 * j]] + b[myMFld2[3 * j + 1]] + b[myMFld2[3 * j + 2]]);
    const MFloat l_innerterm = CSsq + tmp + tmp * tmp * F1B2mulF1BCSsq - squaredVelocityB2;
    const MFloat eq_dist = lb_tp_coef3 * rho * l_innerterm;
    const MFloat l_innerterm_te = tmp + tmp * tmp * F1BCSsq - squaredVelocity + F1B2 * (F3 - nDim * CSsq);
    eqDistT[myDistFld[0] + myDistFld[1] + j] = lb_tp_coef3 * p0 * l_innerterm_te + totalEnergy * eq_dist;
  }

  // Rest distribution
  const MFloat l_innerterm = F1 - F1B2mulF1BCSsq * squaredVelocity;
  const MFloat eq_dist = myTp[0] * rho * l_innerterm;
  const MFloat l_innerterm_te = -F1BCSsq * squaredVelocity + F1B2mulF1BCSsq * (F0 - nDim * CSsq);
  eqDistT[Ld::lastId()] = myTp[0] * p0 * l_innerterm_te + totalEnergy * eq_dist;
}

template <MInt nDim, MInt nDist>
inline void calcEqDistsTotalEnergy(const MFloat& T, const MFloat& rho, const MFloat* const u, MFloat* const eqDist,
                                   const MInt* myMFld1, const MInt* myMFld2, const MFloat* myTp,
                                   const MInt* myDistFld) {
  const MFloat squaredVelocity = std::inner_product(&u[0], &u[nDim], &u[0], .0);

  calcEqDistsTotalEnergy<nDim, nDist>(T, rho, squaredVelocity, u, eqDist, myMFld1, myMFld2, myTp, myDistFld);
}
#endif

/**
 * \brief Calculate transport equilibrium distribution
 *
 * Attention: This is the basic transport LBM
 *
 * Use this version if the squaredVelocity is already precomputed
 *
 * \author Shota Ito, Julian Vorspohl
 * \date 07.06.2022
 *
 * \param[in]  C      Macroscopic concentration
 * \param[out] eqDist Resulting equilibrium distribution
 */
template <MInt nDim, MInt nDist>
inline void calcEqDistsTransport(const MFloat& C, const MFloat squaredVelocity, const MFloat* const u,
                                 MFloat* const eqDist) {
  calcEqDists<nDim, nDist, true>(C, squaredVelocity, u, eqDist);
}

/**
 * \brief Calculate transport equilibrium distribution
 *
 * Attention: This is the basic transport LBM
 *
 * \author Shota Ito, Julian Vorspohl
 * \date 07.06.2022
 *
 * \param[in]  C      Macroscopic concentration
 * \param[out] eqDist Resulting equilibrium distribution
 */
template <MInt nDim, MInt nDist>
inline void calcEqDistsTransport(const MFloat& C, const MFloat* const u, MFloat* const eqDist) {
  const MFloat squaredVelocity = std::inner_product(&u[0], &u[nDim], &u[0], .0);

  calcEqDists<nDim, nDist, true>(C, squaredVelocity, u, eqDist);
}

#ifdef WAR_NVHPC_PSTL
/**
 * \brief Calculate transport equilibrium distribution
 *
 * Attention: This is the basic transport LBM
 *
 * Use this version if the squaredVelocity is already precomputed
 *
 * \author Julian Vorspohl <j.vorspohl@aia.rwth-aachen.de>
 * \date 29.09.2021
 *
 * \param[in]  T      Macroscopic temperature
 * \param[out] eqDist Resulting equilibrium distribution
 */
template <MInt nDim, MInt nDist>
inline void calcEqDistsTransport(const MFloat& C, const MFloat squaredVelocity, const MFloat* const u,
                                 MFloat* const eqDist, const MInt* myMFld1, const MInt* myMFld2, const MFloat* myTp,
                                 const MInt* myDistFld) {
  calcEqDists<nDim, nDist, true>(C, squaredVelocity, u, eqDist, myMFld1, myMFld2, myTp, myDistFld);
}

/**
 * \brief Calculate transport equilibrium distribution
 *
 * Attention: This is the basic transport LBM
 *
 * \author Julian Vorspohl <j.vorspohl@aia.rwth-aachen.de>
 * \date 29.09.2021
 *
 * \param[in]  T      Macroscopic temperature
 * \param[out] eqDist Resulting equilibrium distribution
 */
template <MInt nDim, MInt nDist>
inline void calcEqDistsTransport(const MFloat& C, const MFloat* const u, MFloat* const eqDist, const MInt* myMFld1,
                                 const MInt* myMFld2, const MFloat* myTp, const MInt* myDistFld) {
  const MFloat squaredVelocity = std::inner_product(&u[0], &u[nDim], &u[0], .0);

  calcEqDists<nDim, nDist, true>(C, squaredVelocity, u, eqDist, myMFld1, myMFld2, myTp, myDistFld);
}
#endif

/** \brief  Calculate macroscopic variables from distribution
 *  \author Miro Gondrum
 *  \date   14.01.2020
 *  \param[in]  p_dist  distributions
 *  \param[out] p_rho   macroscopic density
 *  \param[out] p_u     macroscopic velocity (u*rho)
 */
template <MInt nDim, MInt nDist, MBool compressible = false>
inline void calcMacroVars(MFloat const* const p_dist, MFloat& p_rho, MFloat* const p_u) {
  using Ld = LbLatticeDescriptor<nDim, nDist>;
  // density
  p_rho = 0.0;
  for(MInt j = 0; j < nDist; j++) {
    p_rho += p_dist[j];
  }
  // velocities [rho*u]
  for(MInt i = 0; i < nDim; i++) {
    p_u[i] = 0.0;
    for(MInt j = 0; j < Ld::dxQyFld(); j++) {
      p_u[i] += p_dist[Ld::pFld(i, j)];
      p_u[i] -= p_dist[Ld::nFld(i, j)];
    }
    if constexpr(compressible) p_u[i] /= p_rho;
  }
}

#ifdef WAR_NVHPC_PSTL
/** \brief  Calculate macroscopic variables from distribution
 *  \author Miro Gondrum
 *  \date   14.01.2020
 *  \param[in]  p_dist  distributions
 *  \param[out] p_rho   macroscopic density
 *  \param[out] p_u     macroscopic velocity
 */
template <MInt nDim, MInt nDist, MBool compressible = false>
inline void calcMacroVars(MFloat const* const p_dist, MFloat& p_rho, MFloat* const p_u, MInt* const myPFld,
                          MInt* const myNFld, const MInt fldlen) {
  using Ld = LbLatticeDescriptor<nDim, nDist>;
  // density
  p_rho = 0.0;
  for(MInt j = 0; j < nDist; j++) {
    p_rho += p_dist[j];
  }
  // velocities [rho*u]
  for(MInt i = 0; i < nDim; i++) {
    p_u[i] = 0.0;
    for(MInt j = 0; j < fldlen; j++) {
      p_u[i] += p_dist[myPFld[i * fldlen + j]];
      p_u[i] -= p_dist[myNFld[i * fldlen + j]];
    }
    if constexpr(compressible) p_u[i] /= p_rho;
  }
}
#endif
} // namespace lbfunc

#endif
