// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "fvcartesianbndrycell.h"

#include "MEMORY/collector.h"
#include "UTIL/debug.h"
#include "fvcartesiansyseqndetchem.h"
#include "fvcartesiansyseqneegas.h"
#include "fvcartesiansyseqnns.h"
#include "fvcartesiansyseqnrans.h"
//#include "UTIL/timer.h"

template <MInt nDim, class SysEqn>
MInt FvBndryCell<nDim, SysEqn>::m_noSpecies;
template <MInt nDim, class SysEqn>
MInt FvBndryCell<nDim, SysEqn>::m_noRansEquations;
template <MInt nDim, class SysEqn>
MInt FvBndryCell<nDim, SysEqn>::m_noEdges;
template <MInt nDim, class SysEqn>
MInt FvBndryCell<nDim, SysEqn>::m_maxNoSurfaces;
template <MInt nDim, class SysEqn>
MInt FvBndryCell<nDim, SysEqn>::m_noNghbrs;
template <MInt nDim, class SysEqn>
MInt FvBndryCell<nDim, SysEqn>::m_noVariables;

using namespace maia::collector_memory;

template <MInt nDim, class SysEqn>
void FvBndryCell<nDim, SysEqn>::init(MInt NotUsed(dimensions), MInt noSpecies, MInt noRansEquations,
                                     MInt /*maxNoCells*/, MInt maxNoSurfaces) {
  TRACE();

  m_noSpecies = noSpecies;
  if(noRansEquations == -1) {
    m_noRansEquations = 0;
  } else {
    m_noRansEquations = noRansEquations;
  }
  m_noEdges = (nDim == 2) ? 4 : 24;
  m_maxNoSurfaces = maxNoSurfaces;
  m_noNghbrs = 2 * nDim;
  m_noVariables = nDim + 2 + m_noSpecies + m_noRansEquations;
}


template <MInt nDim, class SysEqn>
void FvBndryCell<nDim, SysEqn>::allocateElements(void* cellPtr, void*, const MInt) {
  //  TRACE();

  // Initialize non-static member variables:
  m_linkedCellId = -1;
  m_noSrfcs = 1;

  // Member pointers are set to the cell memory here:
  moveElements(cellPtr);
}

/// Sets cell pointers to the memory location cellPtr:
/// NOTE: this cell is not aligned! cellPtr points to the beginning of the cell
/// memory solver within the rawMemory solver. Conservative aligning requires
/// approx. O(maxNoVarsPerCell*maxNoCells) extra memory.
template <MInt nDim, class SysEqn>
void FvBndryCell<nDim, SysEqn>::moveElements(void* cellPtr) {
  //  TRACE();

  unaligned_cell_wise::rowMajor1D(m_externalFaces, cellPtr, m_noNghbrs);
  unaligned_cell_wise::rowMajor1D(m_associatedSrfc, cellPtr, m_noNghbrs);
  unaligned_cell_wise::rowMajor1D(m_coordinates, cellPtr, nDim);
  unaligned_cell_wise::rowMajor1D(m_masterCoordinates, cellPtr, nDim);

  unaligned_cell_wise::rowMajor2D(m_srfcs, cellPtr, m_maxNoSurfaces, 1);
  for(MInt i = 0; i < m_maxNoSurfaces; ++i) {
    unaligned_cell_wise::rowMajor1D(m_srfcs[i]->m_cutEdge, cellPtr, 2 * m_noEdges);
    unaligned_cell_wise::rowMajor1D(m_srfcs[i]->m_bodyId, cellPtr, 2 * m_noEdges);
    unaligned_cell_wise::rowMajor1D(m_srfcs[i]->m_coordinates, cellPtr, nDim);
    unaligned_cell_wise::rowMajor1D(m_srfcs[i]->m_normalVector, cellPtr, nDim);
    unaligned_cell_wise::rowMajor1D(m_srfcs[i]->m_normalVectorCentroid, cellPtr, nDim);
    unaligned_cell_wise::rowMajor1D(m_srfcs[i]->m_planeVector0, cellPtr, nDim);
    unaligned_cell_wise::rowMajor1D(m_srfcs[i]->m_planeVector1, cellPtr, nDim);
    unaligned_cell_wise::rowMajor2D(m_srfcs[i]->m_cutCoordinates, cellPtr, 2 * m_noEdges, nDim);
  }

  unaligned_cell_wise::rowMajor2D(m_srfcVariables, cellPtr, m_maxNoSurfaces, 1);
  for(MInt i = 0; i < m_maxNoSurfaces; ++i) {
    unaligned_cell_wise::rowMajor1D(m_srfcVariables[i]->m_srfcId, cellPtr, nDim);
    unaligned_cell_wise::rowMajor1D(m_srfcVariables[i]->m_imageCoordinates, cellPtr, nDim);
    unaligned_cell_wise::rowMajor1D(m_srfcVariables[i]->m_imageVariables, cellPtr, m_noVariables);
    unaligned_cell_wise::rowMajor1D(m_srfcVariables[i]->m_variablesType, cellPtr, m_noVariables);
    unaligned_cell_wise::rowMajor1D(m_srfcVariables[i]->m_primVars, cellPtr, m_noVariables);
    unaligned_cell_wise::rowMajor1D(m_srfcVariables[i]->m_normalDeriv, cellPtr, m_noVariables);
    // Invoke placement new operator to call constructor for std::vector on
    // place in memory that was allocated before. If this is not done,
    // std::vector will not work properly
    new(&m_srfcVariables[i]->m_imagePointRecConst) std::vector<MFloat>();
  }
}

// Explicit instantiations for 2D and 3D
template class FvBndryCell<2, FvSysEqnNS<2>>;
template class FvBndryCell<3, FvSysEqnNS<3>>;
template class FvBndryCell<2, FvSysEqnDetChem<2>>;
template class FvBndryCell<3, FvSysEqnDetChem<3>>;
template class FvBndryCell<2, FvSysEqnRANS<2, RANSModelConstants<RANS_SA_DV>>>;
template class FvBndryCell<3, FvSysEqnRANS<3, RANSModelConstants<RANS_SA_DV>>>;
template class FvBndryCell<2, FvSysEqnRANS<2, RANSModelConstants<RANS_FS>>>;
template class FvBndryCell<3, FvSysEqnRANS<3, RANSModelConstants<RANS_FS>>>;
template class FvBndryCell<2, FvSysEqnRANS<2, RANSModelConstants<RANS_KOMEGA>>>;
template class FvBndryCell<3, FvSysEqnRANS<3, RANSModelConstants<RANS_KOMEGA>>>;
template class FvBndryCell<3, FvSysEqnEEGas<3>>;
