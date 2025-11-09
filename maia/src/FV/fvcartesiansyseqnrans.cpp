// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "fvcartesiansyseqnrans.h"
#include "MEMORY/alloc.h"


template <>
std::vector<MString> FvSysEqnRANS<2, RANSModelConstants<RANS_SA_DV>>::PrimitiveVariables::varNames = {"u", "v",  "rho",
                                                                                                      "p", "nu", "c"};
template <>
std::vector<MString> FvSysEqnRANS<2, RANSModelConstants<RANS_FS>>::PrimitiveVariables::varNames = {"u", "v",  "rho",
                                                                                                   "p", "nu", "c"};
template <>
std::vector<MString> FvSysEqnRANS<2, RANSModelConstants<RANS_KOMEGA>>::PrimitiveVariables::varNames = {
    "u", "v", "w", "rho", "p", "k", "omega", "c"};
template <>
std::vector<MString> FvSysEqnRANS<3, RANSModelConstants<RANS_SA_DV>>::PrimitiveVariables::varNames = {
    "u", "v", "w", "rho", "p", "nu", "c"};
template <>
std::vector<MString> FvSysEqnRANS<3, RANSModelConstants<RANS_FS>>::PrimitiveVariables::varNames = {
    "u", "v", "w", "rho", "p", "nu", "c"};
template <>
std::vector<MString> FvSysEqnRANS<3, RANSModelConstants<RANS_KOMEGA>>::PrimitiveVariables::varNames = {
    "u", "v", "w", "rho", "p", "k", "omega", "c"};

template <MInt nDim, class RANSModel>
FvSysEqnRANS<nDim, RANSModel>::FvSysEqnRANS(const MInt solverId, const MInt noSpecies)
  : FvSysEqnNS<nDim>(solverId, noSpecies) {
  CV = new ConservativeVariables(noSpecies);
  PV = new PrimitiveVariables(noSpecies);
  FV = new FluxVariables(noSpecies);
}

template <MInt nDim, class RANSModel>
FvSysEqnRANS<nDim, RANSModel>::ConservativeVariables::ConservativeVariables(const MInt noSpecies)
  : m_noSpecies(noSpecies), noVariables(nDim + 2 + m_noRansEquations + noSpecies) {
  if(m_noRansEquations > 0) {
    mAlloc(RHO_NN, m_noRansEquations, "FvSysEqnRANS::ConservativeVariables::RHO_NN", AT_);
    for(MInt i = 0; i < m_noRansEquations; ++i) {
      RHO_NN[i] = RHO_N + i;
    }
  }
  if(m_noSpecies > 0) {
    mAlloc(RHO_Y, m_noSpecies, "FvSysEqnRANS::ConservativeVariables::RHO_Y", AT_);
    for(MUint i = 0; i < m_noSpecies; ++i) {
      RHO_Y[i] = RHO_C + i;
    }
  }
}

template <MInt nDim, class RANSModel>
FvSysEqnRANS<nDim, RANSModel>::ConservativeVariables::~ConservativeVariables() {
  mDeallocate(RHO_NN);
  mDeallocate(RHO_Y);
}

template <MInt nDim, class RANSModel>
FvSysEqnRANS<nDim, RANSModel>::FluxVariables::FluxVariables(const MInt noSpecies) : ConservativeVariables(noSpecies) {}

template <MInt nDim, class RANSModel>
FvSysEqnRANS<nDim, RANSModel>::PrimitiveVariables::PrimitiveVariables(const MInt noSpecies)
  : m_noSpecies(noSpecies), noVariables(nDim + 2 + m_noRansEquations + noSpecies) {
  if(m_noRansEquations > 0) {
    mAlloc(NN, m_noRansEquations, "FvSysEqnRANS::PrimitiveVariables::NN", AT_);
    for(MInt i = 0; i < m_noRansEquations; ++i) {
      NN[i] = N + i;
    }
  }
  if(m_noSpecies > 0) {
    mAlloc(Y, m_noSpecies, "FvSysEqnRANS::PrimitiveVariables::Y", AT_);
    for(MUint i = 0; i < m_noSpecies; ++i) {
      Y[i] = C + i;
    }
  }
}

template <MInt nDim, class RANSModel>
FvSysEqnRANS<nDim, RANSModel>::PrimitiveVariables::~PrimitiveVariables() {
  mDeallocate(NN);
  mDeallocate(Y);
}

template <MInt nDim, class RANSModel>
void FvSysEqnRANS<nDim, RANSModel>::PrimitiveVariables::getPrimitiveVariableNames(MString* names) {
  TRACE();

  for(MInt i = 0; i < noVariables; i++) {
    names[i] = varNames[i];
  }
}

template class FvSysEqnRANS<2, RANSModelConstants<RANS_SA_DV>>;
template class FvSysEqnRANS<3, RANSModelConstants<RANS_SA_DV>>;
template class FvSysEqnRANS<2, RANSModelConstants<RANS_FS>>;
template class FvSysEqnRANS<3, RANSModelConstants<RANS_FS>>;
template class FvSysEqnRANS<2, RANSModelConstants<RANS_KOMEGA>>;
template class FvSysEqnRANS<3, RANSModelConstants<RANS_KOMEGA>>;
