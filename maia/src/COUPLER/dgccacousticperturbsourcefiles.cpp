// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "dgccacousticperturbsourcefiles.h"

const constexpr MInt DgCcAcousticPerturb<2, FvSysEqnNS<2>>::MV::LAMB0[DgCcAcousticPerturb<2, FvSysEqnNS<2>>::MV::nDim];
const constexpr MInt DgCcAcousticPerturb<3, FvSysEqnNS<3>>::MV::LAMB0[DgCcAcousticPerturb<3, FvSysEqnNS<3>>::MV::nDim];

const constexpr MInt DgCcAcousticPerturb<2, FvSysEqnNS<2>>::MV::UU0[DgCcAcousticPerturb<2, FvSysEqnNS<2>>::MV::nDim];
const constexpr MInt DgCcAcousticPerturb<3, FvSysEqnNS<3>>::MV::UU0[DgCcAcousticPerturb<3, FvSysEqnNS<3>>::MV::nDim];

const constexpr MInt
    DgCcAcousticPerturb<2, FvSysEqnNS<2>>::MV::VORT0[DgCcAcousticPerturb<2, FvSysEqnNS<2>>::MV::noVorticities];
const constexpr MInt
    DgCcAcousticPerturb<3, FvSysEqnNS<3>>::MV::VORT0[DgCcAcousticPerturb<3, FvSysEqnNS<3>>::MV::noVorticities];

const constexpr MInt DgCcAcousticPerturb<2, FvSysEqnNS<2>>::MV::DC0[DgCcAcousticPerturb<2, FvSysEqnNS<2>>::MV::nDim];
const constexpr MInt DgCcAcousticPerturb<3, FvSysEqnNS<3>>::MV::DC0[DgCcAcousticPerturb<3, FvSysEqnNS<3>>::MV::nDim];

const constexpr MInt
    DgCcAcousticPerturb<2, FvSysEqnNS<2>>::MV::GRADU[POW2(DgCcAcousticPerturb<2, FvSysEqnNS<2>>::MV::nDim)];
const constexpr MInt
    DgCcAcousticPerturb<3, FvSysEqnNS<3>>::MV::GRADU[POW2(DgCcAcousticPerturb<3, FvSysEqnNS<3>>::MV::nDim)];

const constexpr MInt DgCcAcousticPerturb<2, FvSysEqnNS<2>>::MV::DU[DgCcAcousticPerturb<2, FvSysEqnNS<2>>::MV::nDim];
const constexpr MInt DgCcAcousticPerturb<3, FvSysEqnNS<3>>::MV::DU[DgCcAcousticPerturb<3, FvSysEqnNS<3>>::MV::nDim];

const constexpr MInt DgCcAcousticPerturb<2, FvSysEqnNS<2>>::MV::DRHO[DgCcAcousticPerturb<2, FvSysEqnNS<2>>::MV::nDim];
const constexpr MInt DgCcAcousticPerturb<3, FvSysEqnNS<3>>::MV::DRHO[DgCcAcousticPerturb<3, FvSysEqnNS<3>>::MV::nDim];

const constexpr MInt DgCcAcousticPerturb<2, FvSysEqnNS<2>>::MV::DP[DgCcAcousticPerturb<2, FvSysEqnNS<2>>::MV::nDim];
const constexpr MInt DgCcAcousticPerturb<3, FvSysEqnNS<3>>::MV::DP[DgCcAcousticPerturb<3, FvSysEqnNS<3>>::MV::nDim];

const constexpr MInt DgCcAcousticPerturb<2, FvSysEqnNS<2>>::MV::RHO0;
const constexpr MInt DgCcAcousticPerturb<3, FvSysEqnNS<3>>::MV::RHO0;

const constexpr MInt DgCcAcousticPerturb<2, FvSysEqnNS<2>>::MV::P0;
const constexpr MInt DgCcAcousticPerturb<3, FvSysEqnNS<3>>::MV::P0;

const constexpr MInt DgCcAcousticPerturb<2, FvSysEqnNS<2>>::MV::C0;
const constexpr MInt DgCcAcousticPerturb<3, FvSysEqnNS<3>>::MV::C0;

const constexpr MInt
    DgCcAcousticPerturb<2, FvSysEqnNS<2>>::MV::RHODIVU[DgCcAcousticPerturb<2, FvSysEqnNS<2>>::MV::nDim];
const constexpr MInt
    DgCcAcousticPerturb<3, FvSysEqnNS<3>>::MV::RHODIVU[DgCcAcousticPerturb<3, FvSysEqnNS<3>>::MV::nDim];

const constexpr MInt
    DgCcAcousticPerturb<2, FvSysEqnNS<2>>::MV::UGRADRHO[DgCcAcousticPerturb<2, FvSysEqnNS<2>>::MV::nDim];
const constexpr MInt
    DgCcAcousticPerturb<3, FvSysEqnNS<3>>::MV::UGRADRHO[DgCcAcousticPerturb<3, FvSysEqnNS<3>>::MV::nDim];

const constexpr MInt
    DgCcAcousticPerturb<2, FvSysEqnNS<2>>::MV::GRADPRHO[DgCcAcousticPerturb<2, FvSysEqnNS<2>>::MV::nDim];
const constexpr MInt
    DgCcAcousticPerturb<3, FvSysEqnNS<3>>::MV::GRADPRHO[DgCcAcousticPerturb<3, FvSysEqnNS<3>>::MV::nDim];

const constexpr MInt DgCcAcousticPerturb<2, FvSysEqnNS<2>>::MV::UGRADU[DgCcAcousticPerturb<2, FvSysEqnNS<2>>::MV::nDim];
const constexpr MInt DgCcAcousticPerturb<3, FvSysEqnNS<3>>::MV::UGRADU[DgCcAcousticPerturb<3, FvSysEqnNS<3>>::MV::nDim];

// ####################################################################################

const constexpr MInt
    DgCcAcousticPerturb<2, FvSysEqnDetChem<2>>::MV::LAMB0[DgCcAcousticPerturb<2, FvSysEqnDetChem<2>>::MV::nDim];
const constexpr MInt
    DgCcAcousticPerturb<3, FvSysEqnDetChem<3>>::MV::LAMB0[DgCcAcousticPerturb<3, FvSysEqnDetChem<3>>::MV::nDim];

const constexpr MInt
    DgCcAcousticPerturb<2, FvSysEqnDetChem<2>>::MV::UU0[DgCcAcousticPerturb<2, FvSysEqnDetChem<2>>::MV::nDim];
const constexpr MInt
    DgCcAcousticPerturb<3, FvSysEqnDetChem<3>>::MV::UU0[DgCcAcousticPerturb<3, FvSysEqnDetChem<3>>::MV::nDim];

const constexpr MInt DgCcAcousticPerturb<
    2, FvSysEqnDetChem<2>>::MV::VORT0[DgCcAcousticPerturb<2, FvSysEqnDetChem<2>>::MV::noVorticities];
const constexpr MInt DgCcAcousticPerturb<
    3, FvSysEqnDetChem<3>>::MV::VORT0[DgCcAcousticPerturb<3, FvSysEqnDetChem<3>>::MV::noVorticities];

const constexpr MInt
    DgCcAcousticPerturb<2, FvSysEqnDetChem<2>>::MV::DC0[DgCcAcousticPerturb<2, FvSysEqnDetChem<2>>::MV::nDim];
const constexpr MInt
    DgCcAcousticPerturb<3, FvSysEqnDetChem<3>>::MV::DC0[DgCcAcousticPerturb<3, FvSysEqnDetChem<3>>::MV::nDim];

const constexpr MInt
    DgCcAcousticPerturb<2, FvSysEqnDetChem<2>>::MV::GRADU[POW2(DgCcAcousticPerturb<2, FvSysEqnDetChem<2>>::MV::nDim)];
const constexpr MInt
    DgCcAcousticPerturb<3, FvSysEqnDetChem<3>>::MV::GRADU[POW2(DgCcAcousticPerturb<3, FvSysEqnDetChem<3>>::MV::nDim)];

const constexpr MInt
    DgCcAcousticPerturb<2, FvSysEqnDetChem<2>>::MV::DU[DgCcAcousticPerturb<2, FvSysEqnDetChem<2>>::MV::nDim];
const constexpr MInt
    DgCcAcousticPerturb<3, FvSysEqnDetChem<3>>::MV::DU[DgCcAcousticPerturb<3, FvSysEqnDetChem<3>>::MV::nDim];

const constexpr MInt
    DgCcAcousticPerturb<2, FvSysEqnDetChem<2>>::MV::DRHO[DgCcAcousticPerturb<2, FvSysEqnDetChem<2>>::MV::nDim];
const constexpr MInt
    DgCcAcousticPerturb<3, FvSysEqnDetChem<3>>::MV::DRHO[DgCcAcousticPerturb<3, FvSysEqnDetChem<3>>::MV::nDim];

const constexpr MInt
    DgCcAcousticPerturb<2, FvSysEqnDetChem<2>>::MV::DP[DgCcAcousticPerturb<2, FvSysEqnDetChem<2>>::MV::nDim];
const constexpr MInt
    DgCcAcousticPerturb<3, FvSysEqnDetChem<3>>::MV::DP[DgCcAcousticPerturb<3, FvSysEqnDetChem<3>>::MV::nDim];

const constexpr MInt DgCcAcousticPerturb<2, FvSysEqnDetChem<2>>::MV::RHO0;
const constexpr MInt DgCcAcousticPerturb<3, FvSysEqnDetChem<3>>::MV::RHO0;

const constexpr MInt DgCcAcousticPerturb<2, FvSysEqnDetChem<2>>::MV::P0;
const constexpr MInt DgCcAcousticPerturb<3, FvSysEqnDetChem<3>>::MV::P0;

const constexpr MInt DgCcAcousticPerturb<2, FvSysEqnDetChem<2>>::MV::C0;
const constexpr MInt DgCcAcousticPerturb<3, FvSysEqnDetChem<3>>::MV::C0;

const constexpr MInt
    DgCcAcousticPerturb<2, FvSysEqnDetChem<2>>::MV::RHODIVU[DgCcAcousticPerturb<2, FvSysEqnDetChem<2>>::MV::nDim];
const constexpr MInt
    DgCcAcousticPerturb<3, FvSysEqnDetChem<3>>::MV::RHODIVU[DgCcAcousticPerturb<3, FvSysEqnDetChem<3>>::MV::nDim];

const constexpr MInt
    DgCcAcousticPerturb<2, FvSysEqnDetChem<2>>::MV::UGRADRHO[DgCcAcousticPerturb<2, FvSysEqnDetChem<2>>::MV::nDim];
const constexpr MInt
    DgCcAcousticPerturb<3, FvSysEqnDetChem<3>>::MV::UGRADRHO[DgCcAcousticPerturb<3, FvSysEqnDetChem<3>>::MV::nDim];

const constexpr MInt
    DgCcAcousticPerturb<2, FvSysEqnDetChem<2>>::MV::GRADPRHO[DgCcAcousticPerturb<2, FvSysEqnDetChem<2>>::MV::nDim];
const constexpr MInt
    DgCcAcousticPerturb<3, FvSysEqnDetChem<3>>::MV::GRADPRHO[DgCcAcousticPerturb<3, FvSysEqnDetChem<3>>::MV::nDim];

const constexpr MInt
    DgCcAcousticPerturb<2, FvSysEqnDetChem<2>>::MV::UGRADU[DgCcAcousticPerturb<2, FvSysEqnDetChem<2>>::MV::nDim];
const constexpr MInt
    DgCcAcousticPerturb<3, FvSysEqnDetChem<3>>::MV::UGRADU[DgCcAcousticPerturb<3, FvSysEqnDetChem<3>>::MV::nDim];

template class DgCcAcousticPerturb<2, FvSysEqnNS<2>>;
template class DgCcAcousticPerturb<3, FvSysEqnNS<3>>;
template class DgCcAcousticPerturb<2, FvSysEqnDetChem<2>>;
template class DgCcAcousticPerturb<3, FvSysEqnDetChem<3>>;
