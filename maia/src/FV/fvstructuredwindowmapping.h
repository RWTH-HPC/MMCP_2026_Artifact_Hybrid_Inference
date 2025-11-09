// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef FVSOLVERSTRUCTWINDOWMAPPING
#define FVSOLVERSTRUCTWINDOWMAPPING

#include "INCLUDE/maiatypes.h"

template <MInt nDim>
class StructuredWindowMap {
 public:
  StructuredWindowMap();
  ~StructuredWindowMap(){};
  MInt Id1 = -1;
  MInt start1[nDim]{};
  MInt end1[nDim]{};
  MInt step1[nDim]{};
  MInt Id2 = -1;
  MInt start2[nDim]{};
  MInt end2[nDim]{};
  MInt step2[nDim]{};
  MInt order[nDim]{};
  MInt BC = -1;
  MInt BCsingular[6]{};
  MInt face = -1;
  MInt dir = -1;
  MInt dc1 = -1;
  MInt dc2 = -1;
  MInt originShape = -1;

  // sponge properties:
  MBool hasSponge = false;
  MFloat spongeThickness = F0;
  MFloat beta = F0;
  MFloat sigma = F0;

  // singularity
  MInt Nstar = -1;
  MInt SingularId = -1;
  MInt SingularBlockId[4];

  // fluid-porous interface
  MBool isFluidPorousInterface = false;
};

#endif
