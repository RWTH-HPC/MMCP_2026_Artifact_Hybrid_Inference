// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef DGGRIDMAP_H_
#define DGGRIDMAP_H_

#include "INCLUDE/maiatypes.h"

namespace maia {
namespace dg {

struct GridMapOffset {
  MInt m_firstTargetCellId;  //!< Cell id on current domain (target grid)
  MInt m_noTargetCells;      //!< Number of consecutive mapped target cells
  MInt m_firstDonorGlobalId; //!< First cell id of mapped (donor) grid
  MInt m_noDonorCells;       //!< Number of consecutive mapped donor cells
  // For each target cell, store the number of offspring of the donor cell that
  // are mapped to this cell
  std::vector<MInt> m_noOffspring;
};

} // namespace dg
} // namespace maia

#endif // ifndef DGGRIDMAP_H_
