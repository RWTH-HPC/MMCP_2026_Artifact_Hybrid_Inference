// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef RANSMODELCONSTANTS_H
#define RANSMODELCONSTANTS_H

#include "INCLUDE/maiatypes.h"
#include "enums.h"

template <RansMethod _>
struct RANSModelConstants {
  static constexpr MInt ransModel = NORANS;
  static constexpr MInt noRansEquations = 0;
};

template <>
struct RANSModelConstants<RANS_SA_DV> {
  static constexpr MInt ransModel = RANS_SA_DV; // 1;
  static constexpr MInt noRansEquations = 1;
  static constexpr MFloat chi = 0.1;
  static constexpr MFloat Fsigma = 1.5;
  static constexpr MFloat cb2 = 0.622;
  static constexpr MFloat Fkap2 = 5.94884;
  static constexpr MFloat cb1 = 0.208;  // for lower Reynolds numbers; // choose 0.1355 for higher Reynolds numbers
  static constexpr MFloat cw1 = 3.2391; // choose 3.6704 for lower Reynolds numbers cb1/kappa^2+(1+cb2)/sigma
  static constexpr MFloat cw2 = 0.3;
  static constexpr MFloat cw3to6 = 64.0;
  static constexpr MFloat cv1to3 = 357.911; // 220.0;
  static constexpr MFloat Ft2 = 0.0;
};
using RM_SA_DV = RANSModelConstants<RANS_SA_DV>;

template <>
struct RANSModelConstants<RANS_FS> {
  static constexpr MInt ransModel = RANS_FS;
  static constexpr MInt noRansEquations = 1;
  static constexpr MFloat chi = 0.1;
  static constexpr MFloat fabetcs = 0.09;
  static constexpr MFloat fabetc = 0.072;
  static constexpr MFloat fapsik1 = 650.0;
  static constexpr MFloat fapsik2 = 100.0;
  static constexpr MFloat fapsio1 = 70.0; // 85.0;
  static constexpr MFloat fapsio2 = 80.0; // 100.0;
  static constexpr MFloat fasigma = 1.2;  // 1.2;
  static constexpr MFloat faalpha = 0.29; // 0.26;
  static constexpr MFloat facv1 = 9.1;
  static constexpr MFloat facv1to3 = 753.571;
  static constexpr MFloat faphi = 0.018; //(fabetcs-fabetc)
};
using RM_FS = RANSModelConstants<RANS_FS>;

template <>
struct RANSModelConstants<RANS_KOMEGA> {
  static constexpr MInt ransModel = RANS_KOMEGA;
  static constexpr MInt noRansEquations = 3;
  static constexpr MFloat alpha = 5.0 / 9.0; // 13.0/25.0;//;
  static constexpr MFloat beta = 0.075;
  static constexpr MFloat beta0 = 0.0708;
  static constexpr MFloat betas = 0.09;
  static constexpr MFloat sigma_k = 0.5;
  static constexpr MFloat sigma_o = 0.5;
  static constexpr MFloat sigma_d = 1.0 / 8.0;
  static constexpr MFloat fbeta1 = 85.0;
  static constexpr MFloat fbeta2 = 100.0;
  static constexpr MFloat C_lim = 7.0 / 8.0;
};
using RM_KOMEGA = RANSModelConstants<RANS_KOMEGA>;

template <>
struct RANSModelConstants<RANS_SST> {
  static constexpr MInt ransModel = RANS_SST;
  static constexpr MInt noRansEquations = 2;
  static constexpr MFloat a1 = 0.31;
  static constexpr MFloat beta1 = 0.075;
  static constexpr MFloat betas = 0.09;
  static constexpr MFloat kappa = 0.41;
  static constexpr MFloat sigma_k1 = 0.85;
  static constexpr MFloat sigma_k2 = 1.0;
  static constexpr MFloat sigma_o1 = 0.5;
  static constexpr MFloat sigma_o2 = 0.856;
};
using RM_SST = RANSModelConstants<RANS_SST>;

template <>
struct RANSModelConstants<RANS_KEPSILON> {
  static constexpr MInt noRansEquations = RANS_KEPSILON;
  static constexpr MFloat rsigma_k = 1.0;
  static constexpr MFloat rsigma_eps = 0.7692307692; //=1/1.3
  static constexpr MFloat C_mu = 0.09;
  static constexpr MFloat C1 = 1.35;
  static constexpr MFloat C2 = 1.8;
  static constexpr MFloat C3 = 0.0115;
  static constexpr MFloat C4 = 0.5;
};
using RM_KEPS = RANSModelConstants<RANS_KEPSILON>;

inline MInt noRansEquations(RansMethod ransMethod) {
  switch(ransMethod) {
    case NORANS:
      return RANSModelConstants<NORANS>::noRansEquations;
    case RANS_SA_DV:
      return RANSModelConstants<RANS_SA_DV>::noRansEquations;
    case RANS_FS:
      return RANSModelConstants<RANS_FS>::noRansEquations;
    case RANS_KOMEGA:
      return RANSModelConstants<RANS_KOMEGA>::noRansEquations;
    case RANS_KEPSILON:
      return RANSModelConstants<RANS_KEPSILON>::noRansEquations;
    default:
      return 0;
      // mTerm(1, AT_, "");
  }
}
#endif
