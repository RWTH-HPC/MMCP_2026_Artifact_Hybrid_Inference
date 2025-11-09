// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "postprocessingfv.cpp"

template class PostProcessingFv<3, FvSysEqnNS<3>>;
template class PostProcessingFv<3, FvSysEqnDetChem<3>>;
template class PostProcessingFv<3, FvSysEqnRANS<3, RANSModelConstants<RANS_SA_DV>>>;
template class PostProcessingFv<3, FvSysEqnRANS<3, RANSModelConstants<RANS_FS>>>;
template class PostProcessingFv<3, FvSysEqnRANS<3, RANSModelConstants<RANS_KOMEGA>>>;
template class PostProcessingFv<3, FvSysEqnEEGas<3>>;