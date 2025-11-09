// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "COMM/mpioverride.h"
#include "fvstructuredsolverwindowinfo.h"
#include "globals.h"

template <MInt nDim>
StructuredWindowMap<nDim>::StructuredWindowMap() {}

template class StructuredWindowMap<2>;
template class StructuredWindowMap<3>;
