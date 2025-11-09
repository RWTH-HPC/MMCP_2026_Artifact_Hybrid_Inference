// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef CARTESIANGRIDPOINT_H
#define CARTESIANGRIDPOINT_H

#include "INCLUDE/maiaconstants.h"
#include "INCLUDE/maiamacro.h"
#include "INCLUDE/maiatypes.h"

/** \brief This class defines a grid point, grid cells consist
 *        of.
 */
template <MInt nDim>
class CartesianGridPoint {
 private:
  static MFloat* m_tmpPointer;

 public:
  static void init(MInt NotUsed(dimensions), MInt /*distributions*/, MInt /*maxSize*/){};

  /** \brief This member stores the coordinates
   *        of the point defining the grid point.
   */
  MFloat* m_coordinates;

  /** \brief This member stores the cellIds of
   *         adjacent Cells
   */
  MInt* m_cellIds;

  MInt m_noAdjacentCells;

  static MInt staticElementSize() { return (sizeof(MFloat) * nDim + sizeof(MInt) * IPOW2(nDim)); }


  void allocateElements(void* memPointer, void* /*dummy*/, MInt& /*dummy2*/) {
    m_tmpPointer = (MFloat*)memPointer;
    m_coordinates = m_tmpPointer;
    m_tmpPointer += nDim;
    m_cellIds = (MInt*)m_tmpPointer;
  };
};
#endif // CARTESIANGRIDPOINT_H
