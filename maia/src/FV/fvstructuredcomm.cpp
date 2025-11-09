// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "fvstructuredcomm.h"
#include "INCLUDE/maiatypes.h"
#include "enums.h"
#include "globals.h"

template <MInt nDim>
StructuredComm<nDim>::StructuredComm(const MInt noVars_,
                                     MFloat* const* const variables_,
                                     const MInt noCells_,
                                     const MInt noPoints_,
                                     const StructuredCommType commType_)
  : noVars(noVars_), variables(variables_), noCells(noCells_), noPoints(noPoints_), commType(commType_) {
  cellBufferSize = noVars * noCells;
  pointBufferSize = nDim * noPoints;

  cellBuffer = std::make_unique<MFloat[]>(mMax(cellBufferSize, 0));
  pointBuffer = std::make_unique<MFloat[]>(mMax(pointBufferSize, 0));

  // mAlloc(cellBuffer, cellBufferSize, "cellBuffer", F0, AT_);
  // mAlloc(pointBuffer, pointBufferSize, "cellBuffer", F0, AT_);

  mpi_request = MPI_REQUEST_NULL;
}

// template <MInt nDim>
// StructuredComm<nDim>::~StructuredComm() {
//   mDeallocate(cellBuffer);
//   mDeallocate(pointBuffer);
// }

template class StructuredComm<2>;
template class StructuredComm<3>;
