// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "lptbndrycell.h"

#include "MEMORY/collector.h"
#include "UTIL/debug.h"

template <MInt nDim>
MInt LPTBndryCell<nDim>::m_maxNoSurfaces;

using namespace maia::collector_memory;

template <MInt nDim>
void LPTBndryCell<nDim>::init(MInt NotUsed(dimensions), MInt NotUsed(noSpecies), MInt NotUsed(noRansEquations),
                              MInt /*maxNoCells*/, MInt maxNoSurfaces) {
  TRACE();

  m_maxNoSurfaces = maxNoSurfaces;
}

template <MInt nDim>
void LPTBndryCell<nDim>::allocateElements(void* cellPtr, void*, const MInt) {
  m_noSrfcs = 1;

  // Member pointers are set to the cell memory here:
  moveElements(cellPtr);
}

template <MInt nDim>
void LPTBndryCell<nDim>::moveElements(void* cellPtr) {
  unaligned_cell_wise::rowMajor1D(m_coordinates, cellPtr, nDim);

  unaligned_cell_wise::rowMajor2D(m_srfcs, cellPtr, m_maxNoSurfaces, 1);
  for(MInt i = 0; i < m_maxNoSurfaces; ++i) {
    unaligned_cell_wise::rowMajor1D(m_srfcs[i]->m_coordinates, cellPtr, nDim);
    unaligned_cell_wise::rowMajor1D(m_srfcs[i]->m_normal, cellPtr, nDim);
    unaligned_cell_wise::rowMajor1D(m_srfcs[i]->m_velocity, cellPtr, nDim);
    unaligned_cell_wise::rowMajor2D(m_srfcs[i]->m_planeCoordinates, cellPtr, nDim, nDim);
  }
}

// Explicit instantiations for 2D and 3D
template class LPTBndryCell<2>;
template class LPTBndryCell<3>;
