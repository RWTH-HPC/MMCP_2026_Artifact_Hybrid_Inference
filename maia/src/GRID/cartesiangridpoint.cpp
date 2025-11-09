// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "cartesiangridpoint.h"

template <MInt nDim>
MFloat* CartesianGridPoint<nDim>::m_tmpPointer;

// Explicit instantiations for 2D and 3D
template class CartesianGridPoint<2>;
template class CartesianGridPoint<3>;
