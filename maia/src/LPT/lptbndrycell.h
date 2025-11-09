// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef LPTBNDRYCELL_H
#define LPTBNDRYCELL_H

#include <vector>
#include "INCLUDE/maiaconstants.h"
#include "INCLUDE/maiatypes.h"

template <MInt nDim>
class LPTBndryCell {
 public:
  static void init(MInt, MInt, MInt, MInt, MInt);

  static MInt m_maxNoSurfaces;

  MInt m_cellId;
  MFloat m_volume;
  MFloat* m_coordinates = nullptr;
  MInt m_noSrfcs;

  struct BodySurface {
    MFloat* m_coordinates = nullptr;
    MFloat* m_normal = nullptr;
    MFloat* m_velocity = nullptr;
    MFloat** m_planeCoordinates = nullptr;
  };

  BodySurface** m_srfcs = nullptr;

  static MInt staticElementSize() {
    return (sizeof(MInt) * (1 + 1)    // cell id and no-surfaces
            + sizeof(MFloat) * (1)    // volume
            + sizeof(MFloat) * (nDim) // coordinates
            + sizeof(MFloat) * m_maxNoSurfaces
                  * (nDim                                       // normal
                     + nDim                                     // coordinates
                     + nDim                                     // velocity
                     + nDim * nDim)                             // planeCoordinates
            + sizeof(MFloat*) * m_maxNoSurfaces * (nDim * nDim) // planeCoordinates
            + sizeof(BodySurface*) * m_maxNoSurfaces + sizeof(BodySurface) * m_maxNoSurfaces);
  }

  void allocateElements(void*, void*, const MInt);
  void moveElements(void*);
};

#endif // ifndef LPTBNDRYCELL_H
