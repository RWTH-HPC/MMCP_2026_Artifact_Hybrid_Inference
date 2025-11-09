// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "variables.h"

#include "MEMORY/alloc.h"

namespace maia {
namespace fv {

////////////////////////////////////////////////////////////////////////////////
/// ConservativeVariables:

template <MInt nd>
const std::array<MInt, 3> ConservativeVariables<nd>::RHO_VV;

template <MInt nd>
ConservativeVariables<nd>::ConservativeVariables(const MInt noSpecies, const MInt noRans)
  : RHO_Z(nd + 2 + noRans),
    RHO_C(nd + 2 + noRans),
    m_noSpecies(noSpecies),
    m_noRansEquations(noRans),
    noVariables(nd + 2 + noRans + noSpecies) {
  if(m_noRansEquations > 0) {
    mAlloc(RHO_NN, m_noRansEquations, "maia::fv::ConservativeVariables::RHO_NN", AT_);
    for(MInt i = 0; i < m_noRansEquations; ++i) {
      RHO_NN[i] = RHO_N + i;
    }
  }
  if(m_noSpecies > 0) {
    mAlloc(RHO_Y, m_noSpecies, "maia::fv::ConservativeVariables::RHO_Y", AT_);
    for(MInt i = 0; i < m_noSpecies; ++i) {
      RHO_Y[i] = RHO_Z + i;
    }
  }
}


template <MInt nd>
ConservativeVariables<nd>::~ConservativeVariables() {
  mDeallocate(RHO_Y);
}

template struct ConservativeVariables<2>;
template struct ConservativeVariables<3>;

////////////////////////////////////////////////////////////////////////////////
/// PrimitiveVariables:

template <MInt nd>
const std::array<MInt, 3> PrimitiveVariables<nd>::VV;

template <MInt nd>
PrimitiveVariables<nd>::PrimitiveVariables(const MInt noSpecies, const MInt noRans)
  : Z(nd + 2 + noRans),
    C(nd + 2 + noRans),
    m_noSpecies(noSpecies),
    m_noRansEquations(noRans),
    noVariables(nd + 2 + noRans + noSpecies) {
  if(m_noRansEquations > 0) {
    mAlloc(NN, m_noRansEquations, "maia::fv::PrimitiveVariables::NN", AT_);
    for(MInt i = 0; i < m_noRansEquations; ++i) {
      NN[i] = N + i;
    }
  }
  if(m_noSpecies > 0) {
    mAlloc(Y, m_noSpecies, "maia::fv::PrimitiveVariables::Y", AT_);
    for(MInt i = 0; i < m_noSpecies; ++i) {
      Y[i] = Z + i;
    }
  }
}

template <MInt nd>
PrimitiveVariables<nd>::~PrimitiveVariables() {
  mDeallocate(Y);
}

template struct PrimitiveVariables<2>;
template struct PrimitiveVariables<3>;

} // namespace fv
} // namespace maia
