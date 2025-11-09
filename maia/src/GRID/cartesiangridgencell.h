// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef GRIDGENCELL_H
#define GRIDGENCELL_H

#include <array>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <string.h>

#include "INCLUDE/maiaconstants.h"
#include "INCLUDE/maiamacro.h"
#include "INCLUDE/maiatypes.h"

template <MInt nDim>
class GridgenPar;

template <MInt nDim>
class GridgenCell {
  template <MInt nDim_>
  friend class GridgenPar;

 public:
  void allocateElements(void*, void* basePtr, const MInt cellId);
  void moveElements(void*, void* basePtr, const MInt cellId);
  void memCopyElements(void* basePtr, void*, void*, MLong);

  static MInt s_maxNoCells;
  static MInt s_dimensions;
  static MInt s_noChilds;
  static MInt s_noNghbrs;

 private:
  static const MUint s_maxNoSolvers = 8;

  MFloat* m_coordinates_;
  MInt* m_level_;
  MLong* m_parentId_;
  MLong* m_globalId_;
  MInt* m_noChildIds_;
  MLong* m_nghbrIds_;
  MLong* m_childIds_;
  MInt m_rfnDistance_;
  MInt* m_noSolidLayer_;

  // bit user  property
  // -----------------------
  // 0   grd   is inside cell
  // 1   grd   is boundary cell
  // 2   grd   is periodic cell
  // 3   grd   is window cell
  // 4   grd   is halo cell
  // 5   grd   is a refine cell
  // 6   grd   is a visited cell
  // 7   grd   is a moved cell, future deletion
  std::array<char, s_maxNoSolvers> b_properties_;

  std::array<char, s_maxNoSolvers> b_solverAffiliation_;
  std::array<char, s_maxNoSolvers> b_solverBoundary_;
  std::array<char, s_maxNoSolvers> b_solverToRefine_;

 public:
  static void init(MInt dimensions, const MInt, MInt maxNoCells) {
    s_dimensions = dimensions;
    s_noChilds = IPOW2(s_dimensions);
    s_noNghbrs = 2 * s_dimensions;
    s_maxNoCells = maxNoCells;
  }

  static MInt staticElementSize() {
    return (sizeof(MLong)
                * (2 +                        // m_globalId, m_parentId
                   2 * s_dimensions +         // m_nghbrIds
                   s_noChilds)                // m_childIds
            + sizeof(MInt) * (2)              // m_level, m_noChildIds
            + sizeof(MFloat) * (s_dimensions) // m_coordinates_
            + sizeof(MInt) * (s_maxNoSolvers) // m_noSolidLayer_ for max of 8 solvers!
    );
  }
};

#endif
