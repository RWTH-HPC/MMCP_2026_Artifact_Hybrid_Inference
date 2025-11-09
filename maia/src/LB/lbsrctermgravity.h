// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "lbsolverdxqy.h"
#include "lbsrcterm.h"

namespace maia::lb {

template <MInt nDim, MInt nDist, class SysEqn>
class LbSrcTermGravity : public LbSrcTerm<nDim, nDist, SysEqn> {
 public:
  LbSrcTermGravity(LbSolverDxQy<nDim, nDist, SysEqn>* p_solver) : m_solver(p_solver) { readProperties(); };

  void readProperties() override;
  void init() override;

  void apply_preCollision() override{};
  void apply_postCollision() override{};
  void apply_postPropagation() override;

 private:
  LbSolverDxQy<nDim, nDist, SysEqn>* m_solver;

  enum struct Mode { DIRECT, GALILEO } m_mode;

  std::array<MFloat, nDim> m_Ga{};

  std::array<MFloat, nDim> m_acceleration{};
  std::array<MFloat, nDist> m_forcing{};
};

} // namespace maia::lb
