// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "fvcartesiansyseqneegas.h"
#include "MEMORY/alloc.h"

template <MInt nDim>
FvSysEqnEEGas<nDim>::FvSysEqnEEGas(const MInt solverId, const MInt noSpecies) : FvSysEqnNS<nDim>(solverId, noSpecies) {
  CV = new ConservativeVariables();
  PV = new PrimitiveVariables();
  FV = new FluxVariables();

  readProperties();
}

template <MInt nDim>
void FvSysEqnEEGas<nDim>::readProperties() {
  m_EEGasEps = 1.0e-10;
  m_EEGasEps = Context::getSolverProperty<MFloat>("EEGasEps", m_solverId, AT_, &m_EEGasEps);
}

template class FvSysEqnEEGas<3>;
