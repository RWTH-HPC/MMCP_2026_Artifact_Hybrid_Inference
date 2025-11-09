// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef POINTBASEDCELL_H
#define POINTBASEDCELL_H

#include "INCLUDE/maiatypes.h"
#include "UTIL/debug.h"

/** \brief This class defines a grid cell which is a internal
 *        data structure for cartesian cells needed for an efficient
 *        meshing procedure.
 *
 *        The grid cells consists of grid points and the grid cells
 *        point to the corresponding cartesian cell they are representing.
 */
template <MInt nDim>
class PointBasedCell {
 public:
  static void init(MInt NotUsed(dimensions), MInt /*distributions*/, MInt /*distributions1*/, MInt /*maxNoCells*/){
      //    TRACE();
  };

  PointBasedCell(){};
  ~PointBasedCell();

  /** \brief This member stores the id of the
   *        cartesian cell in the collector of all
   *        cartesian cells the grid cell points to.
   */
  MInt m_cellId;

  /** \brief This member stores the ids of the grid points
   *        the grid cell conists of.
   */
  MInt* m_pointIds;
  static MInt staticElementSize() { return (sizeof(MInt) * (nDim - 1) * 4); }
  void allocateElements(void* memPointer, void* /*basePointer*/, MInt /*cellId*/) {
    //    TRACE();

    m_pointIds = (MInt*)memPointer;
  };
};
#endif // POINTBASEDCELL_H
