// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef GRIDCELLPROPERTIES_H_
#define GRIDCELLPROPERTIES_H_

#include <bitset>
#include <type_traits>

/// Grid cell Property Labels.
enum class GridCell {
  IsPeriodic,
  IsHalo,
  // IsGhost,
  IsWindow,
  IsToDelete,
  IsPartLvlAncestor, ///< cell is ancestor of partition level
  IsPartitionCell,   ///< cell is a partition cell
  WasNewlyCreated,   ///< cell was recently created
  WasCoarsened,      ///< cell was recently coarsened
  WasRefined,        ///< cell was recently refined
  // <<< add new properties here
  // ---
  NumProperties
};


namespace maia {
namespace grid {
namespace cell {

/// Converts property name to underlying integer value
constexpr std::underlying_type<GridCell>::type p(const GridCell property) {
  return static_cast<std::underlying_type<GridCell>::type>(property);
}

using BitsetType = std::bitset<p(GridCell::NumProperties)>;

} // namespace cell
} // namespace grid
} // namespace maia

#endif // GRIDCELLPROPERTIES_H_
