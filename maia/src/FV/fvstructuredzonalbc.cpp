// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "fvstructuredzonalbc.h"


StructuredZonalComm::StructuredZonalComm(const MInt _noCells, const MInt _localId, const MInt _noVars)
  : bufferSize(_noCells), localId(_localId), noVars(_noVars) {
  mAlloc(mapCellId, bufferSize, "bufferMapCellId", -1, AT_);
  mAlloc(buffer, noVars * bufferSize, "buffer", -1.0, AT_);
}

StructuredZonalBC::StructuredZonalBC()
  : receiverBlockId(-1),
    noCellsGlobalBC(0),
    noCellsLocalBC(0),
    hasLocalBCMap(false),
    hasSTG(false),
    noGlobalRcvDomains(0),
    noGlobalSndDomains(0),
    noSndNghbrDomains(0),
    noRcvNghbrDomains(0),
    noZonalVariables(0) {
  mAlloc(start, 3, "start", 0, AT_);
  mAlloc(end, 3, "end", 0, AT_);
}

StructuredZonalBC::~StructuredZonalBC() {}
