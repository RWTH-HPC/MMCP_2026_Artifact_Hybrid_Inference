// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef FVAPESOLVER3D_H_
#define FVAPESOLVER3D_H_

#include "fvcartesiansolverxd.h"
#include "fvcartesiansyseqnns.h"
//#include "fvcartesiansyseqnrans.h"


class FvApeSolver3D : public FvCartesianSolverXD<3, FvSysEqnNS<3>> {
 public:
  FvApeSolver3D(MInt, MInt, MBool*, maia::grid::Proxy<3>& gridProxy_, Geometry<3>& geometry_, const MPI_Comm comm);
};

#endif // ifndef FVAPESOLVER3D_H_
