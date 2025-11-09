// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef FVSTRUCTUREDCOMM
#define FVSTRUCTUREDCOMM

#include <memory>
#include "COMM/mpioverride.h"
#include "INCLUDE/maiatypes.h"
#include "enums.h"

// This class contains the information needed for the solver to communicate via MPI with other processes.
// It also contains the send/receive buffers

template <MInt nDim>
class StructuredComm {
 public:
  StructuredComm(const MInt, MFloat* const* const, const MInt, const MInt, const StructuredCommType);
  ~StructuredComm(){};

  const MInt noVars;
  MFloat* const* const variables;
  const MInt noCells;
  const MInt noPoints;
  const StructuredCommType commType;

  MInt bcId{};
  MInt tagHelper{};
  MInt nghbrId{};
  MInt cellBufferSize{};
  MInt pointBufferSize{};

  std::unique_ptr<MFloat[]> cellBuffer{};
  std::unique_ptr<MFloat[]> pointBuffer{};

  /* MFloat* cellBuffer = nullptr; */
  /* MFloat* pointBuffer= nullptr; */


  MPI_Request mpi_request;
  MPI_Status mpi_status;

  std::array<MInt, nDim> startInfoCells{};
  std::array<MInt, nDim> endInfoCells{};
  std::array<MInt, nDim> startInfoPoints{};
  std::array<MInt, nDim> endInfoPoints{};

  std::array<MInt, nDim> orderInfo{};
  std::array<MInt, nDim> stepInfo{};
};


#endif
