// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "lbbndcnddxqy.cpp"

// Explicit instantiation
template class LbBndCndDxQy<2, 9, maia::lb::LbSysEqnIncompressible<2, 9>>;
template class LbBndCndDxQy<3, 19, maia::lb::LbSysEqnIncompressible<3, 19>>;
template class LbBndCndDxQy<3, 27, maia::lb::LbSysEqnIncompressible<3, 27>>;
template class LbBndCndDxQy<2, 9, maia::lb::LbSysEqnCompressible<2, 9>>;
template class LbBndCndDxQy<3, 19, maia::lb::LbSysEqnCompressible<3, 19>>;
template class LbBndCndDxQy<3, 27, maia::lb::LbSysEqnCompressible<3, 27>>;
