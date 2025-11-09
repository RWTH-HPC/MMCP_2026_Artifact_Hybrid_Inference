// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef MAIA_PARTICLECOLLDATA_H
#define MAIA_PARTICLECOLLDATA_H

#include "globals.h"

template <MInt nDim>
class LPTSpherical;

template <MInt nDim>
class collStruct {
 public:
  MFloat collTime = -1.0;
  typename std::vector<LPTSpherical<nDim>>::iterator particle0;
  typename std::vector<LPTSpherical<nDim>>::iterator particle1;
  MInt bndryId;

  collStruct(){};

  // constructor for new collision event
  collStruct(const collStruct& copy) {
    collTime = copy.collTime;
    particle0 = copy.particle0;
    particle1 = copy.particle1;
    bndryId = copy.bndryId;
  }

  // destructor
  ~collStruct() = default;

  // assigns value to existing collision event, required by sort

  collStruct& operator=(const collStruct& rhs) {
    this->collTime = rhs.collTime;
    this->particle0 = rhs.particle0;
    this->particle1 = rhs.particle1;
    this->bndryId = rhs.bndryId;
    return *this;
  }

  // compares collision times of two collision events, required by sort

  MInt operator==(const collStruct& rhs) const {
    if(!approx(this->collTime, rhs.collTime, MFloatEps)) {
      return 0;
    }
    return 1;
  }

  // finds the smaller of two collision times, required by sort

  MInt operator<(const collStruct& rhs) const {
    if(this->collTime < rhs.collTime) {
      return 1;
    }
    return 0;
  }

  typename std::vector<LPTSpherical<nDim>>::iterator returnParticle0() const { return particle0; }

  typename std::vector<LPTSpherical<nDim>>::iterator returnParticle1() const { return particle1; }
};

#endif // MAIA_PARTICLECOLLDATA_H
