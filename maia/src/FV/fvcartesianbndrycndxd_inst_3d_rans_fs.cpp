// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "fvcartesianbndrycndxd.cpp"
#include "fvcartesiansyseqnrans.h"

template class Bc1601Class<3>;
template class FvBndryCndXD<3, FvSysEqnRANS<3, RANSModelConstants<RANS_FS>>>;
