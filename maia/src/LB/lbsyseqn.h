// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#pragma once

#include "INCLUDE/maiatypes.h"
#include "lbfunctions.h"

namespace maia::lb {

template <MInt nDim, MInt nDist>
class LbSysEqn {
 public:
  virtual void calcEqDists(const MFloat zerothMoment, const MFloat* const firstMoment,
                           MFloat* const p_eqDist) const = 0;
  virtual void calcEqDists(const MFloat zerothMoment, const MFloat squaredFirstMoment, const MFloat* const firstMoment,
                           MFloat* const p_eqDist) const = 0;
#ifdef WAR_NVHPC_PSTL
  virtual void calcEqDists(const MFloat zerothMoment, const MFloat* const firstMoment, MFloat* const p_eqDist,
                           const MInt* mFld1, const MInt* mFld2, const MFloat* mTp, const MInt* mDistFld) const = 0;
  virtual void calcEqDists(const MFloat zerothMoment, const MFloat squaredFirstMoment, const MFloat* const firstMoment,
                           MFloat* const p_eqDist, const MInt* mFld1, const MInt* mFld2, const MFloat* mTp,
                           const MInt* mDistFld) const = 0;
#endif
  virtual void calcMacroVars(MFloat const* const p_dist, MFloat& p_rho, MFloat* const p_u) const = 0;
  virtual void calcPrimitiveVars(MFloat const* const p_dist, MFloat& p_rho, MFloat* const p_u) const = 0;

  void calcNonEqDists(const MFloat zerothMoment, const MFloat* const firstMoment, const MFloat* const p_dist,
                      MFloat* const p_nonEqDist) const {
    // Calculate equlibirum distribution
    std::array<MFloat, nDist> eqDist{};
    calcEqDists(zerothMoment, firstMoment, eqDist.data());

    // Calculate non-equilibium part
    for(MInt j = 0; j < nDist; j++) {
      p_nonEqDist[j] = p_dist[j] - eqDist[j];
    }
  };

  void calcMomentumFlux(const MFloat& zerothMoment, const MFloat* const firstMoments, const MFloat* const dist,
                        MFloat* const momentumFlux) const {
    using Ld = LbLatticeDescriptor<nDim, nDist>;

    std::array<MFloat, nDist> nonEqDist{};
    calcNonEqDists(zerothMoment, firstMoments, dist, nonEqDist.data());

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
};

template <MInt nDim, MInt nDist>
class LbSysEqnIncompressible : public LbSysEqn<nDim, nDist> {
 public:
  inline void calcEqDists(const MFloat p_rho, const MFloat* const p_u, MFloat* const p_eqDist) const override {
    lbfunc::calcEqDists<nDim, nDist, false>(p_rho, p_u, p_eqDist);
  }
  inline void calcEqDists(const MFloat p_rho, const MFloat squaredVelocity, const MFloat* const p_u,
                          MFloat* const p_eqDist) const override {
    lbfunc::calcEqDists<nDim, nDist, false>(p_rho, squaredVelocity, p_u, p_eqDist);
  }
#ifdef WAR_NVHPC_PSTL
  inline void calcEqDists(const MFloat p_rho, MFloat const* const p_u, MFloat* const eqDist, const MInt* mFld1,
                          const MInt* mFld2, const MFloat* mTp, const MInt* mDistFld) const override {
    lbfunc::calcEqDists<nDim, nDist, false>(p_rho, p_u, eqDist, mFld1, mFld2, mTp, mDistFld);
  }
  inline void calcEqDists(const MFloat p_rho, const MFloat squaredVelocity, MFloat const* const p_u,
                          MFloat* const eqDist, const MInt* mFld1, const MInt* mFld2, const MFloat* mTp,
                          const MInt* mDistFld) const override {
    lbfunc::calcEqDists<nDim, nDist, false>(p_rho, squaredVelocity, p_u, eqDist, mFld1, mFld2, mTp, mDistFld);
  }
#endif
  inline void calcMacroVars(MFloat const* const p_dist, MFloat& p_rho, MFloat* const p_u) const override {
    lbfunc::calcMacroVars<nDim, nDist>(p_dist, p_rho, p_u);
  }
  inline void calcPrimitiveVars(MFloat const* const p_dist, MFloat& p_rho, MFloat* const p_u) const override {
    calcMacroVars(p_dist, p_rho, p_u);
  }
};

template <MInt nDim, MInt nDist>
class LbSysEqnCompressible : public LbSysEqn<nDim, nDist> {
 public:
  inline void calcEqDists(const MFloat p_rho, const MFloat* const p_u, MFloat* const p_eqDist) const override {
    lbfunc::calcEqDists<nDim, nDist, true>(p_rho, p_u, p_eqDist);
  }
  inline void calcEqDists(const MFloat p_rho, const MFloat squaredVelocity, const MFloat* const p_u,
                          MFloat* const p_eqDist) const override {
    lbfunc::calcEqDists<nDim, nDist, true>(p_rho, squaredVelocity, p_u, p_eqDist);
  }
#ifdef WAR_NVHPC_PSTL
  inline void calcEqDists(const MFloat p_rho, MFloat const* const p_u, MFloat* const eqDist, const MInt* mFld1,
                          const MInt* mFld2, const MFloat* mTp, const MInt* mDistFld) const override {
    lbfunc::calcEqDists<nDim, nDist, true>(p_rho, p_u, eqDist, mFld1, mFld2, mTp, mDistFld);
  }
  inline void calcEqDists(const MFloat p_rho, const MFloat squaredVelocity, MFloat const* const p_u,
                          MFloat* const eqDist, const MInt* mFld1, const MInt* mFld2, const MFloat* mTp,
                          const MInt* mDistFld) const override {
    lbfunc::calcEqDists<nDim, nDist, true>(p_rho, squaredVelocity, p_u, eqDist, mFld1, mFld2, mTp, mDistFld);
  }
#endif
  inline void calcMacroVars(MFloat const* const p_dist, MFloat& p_rho, MFloat* const p_u) const override {
    lbfunc::calcMacroVars<nDim, nDist, true>(p_dist, p_rho, p_u);
  }
  inline void calcPrimitiveVars(MFloat const* const p_dist, MFloat& p_rho, MFloat* const p_u) const override {
    calcMacroVars(p_dist, p_rho, p_u);
  }
};

} // namespace maia::lb
