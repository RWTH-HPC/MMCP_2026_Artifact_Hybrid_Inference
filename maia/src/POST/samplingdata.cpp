// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "samplingdata.h"
#include "DG/dgcartesiansolver.h"
#include "DG/dgcartesiansyseqnacousticperturb.h"
#include "DG/dgcartesiansyseqnlinearscalaradv.h"
#include "FV/fvcartesiansolverxd.h"

// Instantiations for DG solver
template class PointData<2, DgCartesianSolver<2, DgSysEqnAcousticPerturb<2>>>;
template class PointData<3, DgCartesianSolver<3, DgSysEqnAcousticPerturb<3>>>;

template class PointData<2, DgCartesianSolver<2, DgSysEqnLinearScalarAdv<2>>>;
template class PointData<3, DgCartesianSolver<3, DgSysEqnLinearScalarAdv<3>>>;

// Instantiations for FV solver (NS and DetChem)
template class PointData<2, FvCartesianSolverXD<2, FvSysEqnNS<2>>>;
template class PointData<3, FvCartesianSolverXD<3, FvSysEqnNS<3>>>;

template class PointData<2, FvCartesianSolverXD<2, FvSysEqnDetChem<2>>>;
template class PointData<3, FvCartesianSolverXD<3, FvSysEqnDetChem<3>>>;
