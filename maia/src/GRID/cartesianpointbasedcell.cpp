// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "cartesianpointbasedcell.h"

template <MInt nDim>
PointBasedCell<nDim>::~PointBasedCell() {
  //  TRACE();

  m_pointIds = 0;
}

// Explicit instantiations for 2D and 3D
template class PointBasedCell<2>;
template class PointBasedCell<3>;
