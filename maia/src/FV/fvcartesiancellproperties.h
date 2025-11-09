// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef FVCELLPROPERTIES_H_
#define FVCELLPROPERTIES_H_

#include <bitset>
#include <type_traits>

/// FV cell Property Labels.
enum class FvCell {
  IsDummy,
  IsInterface, ///< interface cells, set by cartesiangrid for the solvers
  IsCutOff,
  IsNotGradient, ///< FV Multisolver: Halo cells requiring no gradient computation
  IsInvalid,
  IsFlux,       ///< FV Multisolver: Window or Halo cells participating in
                ///< the flux computation (see: FvCartesianSolver::tagCellsNeededForSurfaceFlux)
  IsActive,     ///< FV: see FvCartesianSolver::writeListOfActiveFlowCells
  IsGhost,      ///< ghost cells, set by the solvers for cartesiangrid
  IsSplitChild, ///< cell is a split child
  IsOnCurrentMGLevel,
  IsInSpongeLayer,
  IsPeriodicWithRot,  ///< periodic cell-->new periodic bc
  IsSplitCell,        ///< cell is a split cell
  HasSplitFace,       ///< cell has a split face
  IsTempLinked,       ///< cell is temporarily linked
  IsInactive,         ///< cell is outside fluid domain
  IsGapCell,          ///< cell is inside gap
  WasGapCell,         ///< cell was inside gap
  NearWall,           ///< cell is near solid boundary
  WasInactive,        ///< cell was outside fluid domain
  AtStructuredRegion, ///< cell has 2*nDim regular neighbors on the same level
  IsHalo,             ///< cell is halo-cell
  IsWindow,           ///< cell is window cell for current solver
  IsPeriodic,         ///< cell is periodic-cell
  IsSplitClone,       ///< cell is an older-version of a splitchild
  IsBndryActive,      ///< cell which is active but has an inactive neighbor (=MB-Boundary-Cell)
  IsMovingBnd,        ///< cell at moving boundary
  IsWMImgCell,
  IsSandpaperTripCell,
  IsInsideReactionZone,
  HasCoarseNghbr, ///< cell has a coarse nghbr (only set for cells with IsNotGradient==false)
  // <<< add new properties here
  // ---
  NumProperties
};


namespace maia {
namespace fv {
namespace cell {

/// Converts property name to underlying integer value
constexpr std::underlying_type<FvCell>::type p(const FvCell property) {
  return static_cast<std::underlying_type<FvCell>::type>(property);
}

using BitsetType = std::bitset<p(FvCell::NumProperties)>;

} // namespace cell
} // namespace fv
} // namespace maia

#endif // FVCELLPROPERTIES_H_
