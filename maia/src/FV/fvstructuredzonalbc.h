// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef STRUCTUREDZONALBC_H
#define STRUCTUREDZONALBC_H

#include <memory>
#include <vector>
#include "fvstructuredinterpolation.h"
#include "globals.h"

class StructuredZonalComm {
 public:
  StructuredZonalComm(const MInt, const MInt, const MInt);

  const MInt bufferSize;
  const MInt localId;
  const MInt noVars;
  MFloat* buffer = nullptr;
  MInt* mapCellId = nullptr;
};

class StructuredZonalBC {
 public:
  StructuredZonalBC();
  ~StructuredZonalBC();

  MInt receiverBlockId;
  MInt noCellsGlobalBC;
  MInt noCellsLocalBC;
  MBool hasLocalBCMap;
  MBool hasSTG;

  MInt* globalReceiverIds = nullptr;
  MInt* globalLocalMapCellIds = nullptr;
  MInt* globalRcvDomainIds = nullptr;
  MInt* globalSndDomainIds = nullptr;
  MInt* start = nullptr;
  MInt* end = nullptr;

  std::unique_ptr<StructuredInterpolation<3>> interpolation;

  MInt* hasPartnerLocalBC = nullptr;
  MFloat** interpolatedVars = nullptr;
  MFloat* interpolatedVarsAV = nullptr;

  // communicator

  MInt noGlobalRcvDomains;
  MInt noGlobalSndDomains;
  MInt noSndNghbrDomains;
  MInt noRcvNghbrDomains;
  MInt noBufferSndSize;
  MInt* localCommReceiverIds = nullptr;
  MInt* localBufferMapCellIds = nullptr;
  MInt* localBufferIndexCellIds = nullptr;
  MFloat** coordinatesGlobalBC = nullptr;
  MInt noZonalVariables;

  std::vector<std::unique_ptr<StructuredZonalComm>> sndComm;
  std::vector<std::unique_ptr<StructuredZonalComm>> rcvComm;

  MPI_Status* sndStatus = nullptr;
  MPI_Request* sndRequest = nullptr;
  MPI_Status* rcvStatus = nullptr;
  MPI_Request* rcvRequest = nullptr;
};


#endif
