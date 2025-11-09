// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "lbsrctermcontroller.h"
#include <sstream>
#include "UTIL/functions.h"
#include "lbsyseqn.h"

namespace maia::lb {

template <MInt nDim, MInt nDist, class SysEqn>
LbSrcTermFactory<nDim, nDist, SysEqn>* LbSrcTermFactory<nDim, nDist, SysEqn>::instance() {
  static LbSrcTermFactory<nDim, nDist, SysEqn> fact;
  return &fact;
}

template <MInt nDim, MInt nDist, class SysEqn>
LbSrcTerm<nDim, nDist, SysEqn>*
LbSrcTermFactory<nDim, nDist, SysEqn>::create_srcTerm(const std::string& p_name,
                                                      LbSolverDxQy<nDim, nDist, SysEqn>* p_solver) {
  LbSrcTerm<nDim, nDist, SysEqn>* instance = nullptr;
  auto it = m_function_reg.find(p_name);
  if(it != m_function_reg.end()) {
    // call constructor with relevant input args
    instance = it->second(p_solver);
  } else {
    std::stringstream ss;
    ss << "lbSrcTerms: " << p_name << " is not defined!" << std::endl;
    ss << "Defined are the following lbSrcTerm for d" << nDim << "q" << nDist << ": " << std::endl;
    for(const auto& i : m_function_reg) {
      ss << "     " << i.first << std::endl;
    }
    TERMM(1, ss.str());
  }
  return instance;
}

template class LbSrcTermFactory<2, 9, LbSysEqnIncompressible<2, 9>>;
template class LbSrcTermFactory<3, 19, LbSysEqnIncompressible<3, 19>>;
template class LbSrcTermFactory<3, 27, LbSysEqnIncompressible<3, 27>>;
template class LbSrcTermFactory<2, 9, LbSysEqnCompressible<2, 9>>;
template class LbSrcTermFactory<3, 19, LbSysEqnCompressible<3, 19>>;
template class LbSrcTermFactory<3, 27, LbSysEqnCompressible<3, 27>>;
} // namespace maia::lb
