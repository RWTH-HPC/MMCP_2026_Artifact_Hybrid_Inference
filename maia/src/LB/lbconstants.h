// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef LBCONSTANTS_H
#define LBCONSTANTS_H

#include "INCLUDE/maiaconstants.h"

namespace lbconstants {

const MInt lb_extrapolationPairs[6][10] = {
    {0, 1, 6, 8, 7, 9, 10, 12, 11, 13}, {1, 0, 8, 6, 9, 7, 12, 10, 13, 11},     {2, 3, 6, 7, 8, 9, 14, 16, 15, 17},
    {3, 2, 7, 6, 9, 8, 16, 14, 17, 15}, {4, 5, 10, 11, 12, 13, 14, 15, 16, 17}, {5, 4, 11, 10, 13, 12, 15, 14, 17, 16}};

template <MUint D>
const MString lb_ReStressMeanNames(MUint i) {
  if(D == 2) {
    const MString tmp[3]{"ReStress_uu_mean", "ReStress_uv_mean", "ReStress_vv_mean"};
    return tmp[i];
  } else {
    const MString tmp[6]{"ReStress_uu_mean", "ReStress_uv_mean", "ReStress_uw_mean",
                         "ReStress_vv_mean", "ReStress_vw_mean", "ReStress_ww_mean"};
    return tmp[i];
  }
}

const MFloat lb_gamma2 = 1.4;

} // namespace lbconstants

#endif
