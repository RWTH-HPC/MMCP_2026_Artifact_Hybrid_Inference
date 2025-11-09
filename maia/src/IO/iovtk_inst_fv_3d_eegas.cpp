// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "FV/fvcartesiansyseqneegas.h"
#include "iovtk.cpp"

// Explicit instantiation for 3D
template class VtkIo<3, FvSysEqnEEGas<3>>;
