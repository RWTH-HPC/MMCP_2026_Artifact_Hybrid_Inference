// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "cartesiangridgencell.h"
#include <climits>
#include "globals.h"

template <MInt nDim>
MInt GridgenCell<nDim>::s_maxNoCells;
template <MInt nDim>
MInt GridgenCell<nDim>::s_dimensions;
template <MInt nDim>
MInt GridgenCell<nDim>::s_noChilds;
template <MInt nDim>
MInt GridgenCell<nDim>::s_noNghbrs;

using namespace maia::collector_memory;

template <MInt nDim>
void GridgenCell<nDim>::allocateElements(void*, void* basePtr, const MInt cellId) {
  // init
  m_rfnDistance_ = std::numeric_limits<MInt>::max();
  b_properties_.fill(0);
  b_solverAffiliation_.fill(0);
  b_solverBoundary_.fill(0);
  b_solverToRefine_.fill(0);

  // Member pointers are set to the cell memory here:
  moveElements(nullptr, basePtr, cellId);
}

template <MInt nDim>
void GridgenCell<nDim>::moveElements(void*, void* basePtr, const MInt cellId) {
  rowMajor1D(m_level_, basePtr, cellId, 1, s_maxNoCells);
  rowMajor1D(m_parentId_, basePtr, cellId, 1, s_maxNoCells);
  rowMajor1D(m_globalId_, basePtr, cellId, 1, s_maxNoCells);
  rowMajor1D(m_noChildIds_, basePtr, cellId, 1, s_maxNoCells);
  rowMajor1D(m_childIds_, basePtr, cellId, s_noChilds, s_maxNoCells);
  rowMajor1D(m_nghbrIds_, basePtr, cellId, s_noNghbrs, s_maxNoCells);
  rowMajor1D(m_coordinates_, basePtr, cellId, s_dimensions, s_maxNoCells);
  rowMajor1D(m_noSolidLayer_, basePtr, cellId, s_maxNoSolvers, s_maxNoCells);
}

template <MInt nDim>
void GridgenCell<nDim>::memCopyElements(void* basePtr, void*, void*, MLong) {
  const GridgenCell<nDim>* const from = reinterpret_cast<GridgenCell<nDim>*>(basePtr);

  copyElements1D(m_level_, from->m_level_, 1);
  copyElements1D(m_parentId_, from->m_parentId_, 1);
  copyElements1D(m_globalId_, from->m_globalId_, 1);
  copyElements1D(m_noChildIds_, from->m_noChildIds_, 1);
  copyElements1D(m_childIds_, from->m_childIds_, s_noChilds);
  copyElements1D(m_nghbrIds_, from->m_nghbrIds_, s_noNghbrs);
  copyElements1D(m_coordinates_, from->m_coordinates_, s_dimensions);
  copyElements1D(m_noSolidLayer_, from->m_noSolidLayer_, s_maxNoSolvers);
}


// Explicit instantiations for 2D and 3D
template class GridgenCell<2>;
template class GridgenCell<3>;
