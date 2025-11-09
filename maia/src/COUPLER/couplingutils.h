// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef COUPLINGUTILS_H_
#define COUPLINGUTILS_H_

#include <string>

#include <INCLUDE/maiatypes.h>
#include <UTIL/functions.h>

// Collection of functions which are useful for coupling multiple solvers

/** \brief Conversion from solverA id to the solverB id on the same-level only!
 *    \author Tim Wegmann, Julian Vorspohl
 */
template <class SolverA, class SolverB>
MInt convertId(SolverA& solverA, SolverB& solverB, const MInt solverAId) {
  if(solverAId >= solverA.c_noCells() || solverAId < 0) return -1;
  solverA.assertValidGridCellId(solverAId);

  const MInt gridId = solverA.grid().tree().solver2grid(solverAId);
  ASSERT(solverA.grid().solverFlag(gridId, solverA.solverId()), "");

  if(!solverB.grid().solverFlag(gridId, solverB.solverId())) return -1;

  const MInt solverBId = solverB.grid().tree().grid2solver(gridId);

  if(solverBId > 0) {
    ASSERT(solverA.a_level(solverAId) == solverB.a_level(solverBId),
           std::to_string(solverAId) + " " + std::to_string(solverA.a_level(solverAId)) + " "
               + std::to_string(solverB.a_level(solverBId)));
  }

  return solverBId;
}

/** \brief Conversion from solverA id to the solverB id
 *  If no cell on the same level is found, a connection through the Cell-parents is used!
 *  \author Tim Wegmann, Julian Vorspohl
 */
template <class SolverA, class SolverB>
MInt convertIdParent(SolverA& solverA, SolverB& solverB, const MInt solverAId) {
  if(solverAId >= solverA.c_noCells() || solverAId < 0) return -1; // not even a grid-cell, still return -1
  solverA.assertValidGridCellId(solverAId);

  MInt gridId = solverA.grid().tree().solver2grid(solverAId);
  ASSERT(solverA.grid().solverFlag(gridId, solverA.solverId()), "");

  if(solverB.grid().solverFlag(gridId, solverB.solverId())) { // connection on the same level!
    return solverB.grid().tree().grid2solver(gridId);

  } else { // going through all parents! to find a connection

    MInt solverAParentId = solverA.c_parentId(solverAId);
    MInt solverBParentId = convertId(solverA, solverB, solverAParentId);

    while(solverBParentId < 0 && solverAParentId > -1 && solverAParentId < solverA.c_noCells()) {
      solverAParentId = solverA.c_parentId(solverAParentId);
      solverBParentId = convertId(solverA, solverB, solverAParentId);
    }
    return solverBParentId;
  }
}

#endif // COUPLINGUTILS_H_
