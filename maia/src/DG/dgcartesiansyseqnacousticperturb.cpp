// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "dgcartesiansyseqnacousticperturb.h"

template <>
const MString DgSysEqnAcousticPerturb<2>::s_consVarNames[] = {"u", "v", "p"};
template <>
const MString DgSysEqnAcousticPerturb<2>::s_primVarNames[] = {"u", "v", "p"};
template <>
const MString DgSysEqnAcousticPerturb<2>::s_nodeVarNames[] = {"u0", "v0", "rho0", "c0", "dc0_dx", "dc0_dy"};
template <>
const MString DgSysEqnAcousticPerturb<3>::s_consVarNames[] = {"u", "v", "w", "p"};
template <>
const MString DgSysEqnAcousticPerturb<3>::s_primVarNames[] = {"u", "v", "w", "p"};
template <>
const MString DgSysEqnAcousticPerturb<3>::s_nodeVarNames[] = {"u0", "v0",     "w0",     "rho0",
                                                              "c0", "dc0_dx", "dc0_dy", "dc0_dz"};

const constexpr MInt DgSysEqnAcousticPerturb<2>::CV::UU[DgSysEqnAcousticPerturb<2>::CV::nDim];
const constexpr MInt DgSysEqnAcousticPerturb<3>::CV::UU[DgSysEqnAcousticPerturb<3>::CV::nDim];

const constexpr MInt DgSysEqnAcousticPerturb<2>::CV::UU0[DgSysEqnAcousticPerturb<2>::CV::nDim];
const constexpr MInt DgSysEqnAcousticPerturb<3>::CV::UU0[DgSysEqnAcousticPerturb<3>::CV::nDim];

const constexpr MInt DgSysEqnAcousticPerturb<2>::CV::DC0[DgSysEqnAcousticPerturb<2>::CV::nDim];
const constexpr MInt DgSysEqnAcousticPerturb<3>::CV::DC0[DgSysEqnAcousticPerturb<3>::CV::nDim];
