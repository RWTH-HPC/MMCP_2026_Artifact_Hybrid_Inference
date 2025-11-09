// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef SBPOPERATORS_H_
#define SBPOPERATORS_H_

#include <vector>

#include "INCLUDE/maiatypes.h"
#include "UTIL/tensor.h"

typedef struct {
  MString name;
  MInt globalOrder;

  std::vector<MFloat> a;
  std::vector<MFloat> p;
  std::vector<MFloat> q;

} SBPOperator;

static SBPOperator s306{"go4/s306",
                        4,
                        {7.500000000000000000e-01, -1.499999999999999944e-01, 1.666666666666666644e-02},
                        {3.159490740740625858e-01, 1.390393518518553861e+00, 6.275462962962442548e-01,
                         1.240509259259299890e+00, 9.116898148147989378e-01, 1.013912037037039582e+00},
                        {5.929259109463120847e-01, 1.533619426345958214e-01, -4.394121269149411924e-01,
                         2.310857080667841346e-01, -3.796143473275082059e-02, -9.577947078394533076e-02,
                         1.350220207000018213e+00, -8.388659654351942052e-01, 1.773511401654332964e-01,
                         -6.099231390983241852e-01, 9.759918119376018719e-01, -3.084862009886272793e-01,
                         1.502081835377615249e-01, 1.340100907823247312e-01, 6.517530714402868242e-01}};
#endif
