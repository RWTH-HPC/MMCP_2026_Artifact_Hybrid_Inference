// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef LBSRCTERM_H
#define LBSRCTERM_H

#include "INCLUDE/maiatypes.h"

template <MInt nDim, MInt nDist, class SysEqn>
class LbSolverDxQy;

namespace maia::lb {

/** \brief  Abstract class for lb source terms
 *  \author Miro Gondrum
 *  \date   01.02.2022
 */
template <MInt nDim_, MInt nDist_, class SysEqn_>
class LbSrcTerm {
 protected:
  virtual void readProperties() = 0;

 public:
  static constexpr MInt nDim = nDim_;
  static constexpr MInt nDist = nDist_;
  using SysEqn = SysEqn_;
  virtual void init() = 0;
  virtual void apply_preCollision() = 0;
  virtual void apply_postCollision() = 0;
  virtual void apply_postPropagation() = 0;
  virtual ~LbSrcTerm() = default;
};

} // namespace maia::lb

#endif
