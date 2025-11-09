// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "fvcartesianbndrycndxd.cpp"
#include "fvcartesiansyseqneegas.h"

// The following is required since the class is implemented in the fvbdnrycnd.cpp file
template class Bc1601Class<3>;

template <>
void FvBndryCndXD<3, FvSysEqnEEGas<3>>::cbc1099_1091d(MInt) {
  TERMM(-1, "requires CV->RHO_E");
}

template <>
void FvBndryCndXD<3, FvSysEqnEEGas<3>>::sbc2801x(MInt) {
  TERMM(-1, "requires CV->RHO_E");
}

// Explicit instantiation for 3D
template class FvBndryCndXD<3, FvSysEqnEEGas<3>>;
