// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef LBCELLPROPERTIES_H_
#define LBCELLPROPERTIES_H_

#include <bitset>
#include <type_traits>

/// LB cell Property Labels.
enum class LbCell {
  IsInterface, ///< interface cells, set by cartesiangrid for the solvers
  IsExchange,  ///< FV Particle: Multisolver window or halo cells in
               ///< which particles need to be exchanged
  IsHalo,
  IsWindow,
  IsBndryCell,       ///< LB: cell is a boundary cell
  IsGhost,           ///< LB: cell is a ghost cell
  IsInterfaceParent, ///< LB: cell is a parent of an interface cell
  IsInterfaceChild,  ///< LB: cell is a child of an interface cell
  OnlyBoundary,      ///< LB: ???
  IsActive,
  WasActive,
  WasNewlyCreated,
  // <<< add new properties here
  // ---
  NumProperties
};


namespace maia {
namespace lb {
namespace cell {

/// Converts property name to underlying integer value
constexpr std::underlying_type<LbCell>::type p(const LbCell property) {
  return static_cast<std::underlying_type<LbCell>::type>(property);
}

using BitsetType = std::bitset<p(LbCell::NumProperties)>;

} // namespace cell
} // namespace lb
} // namespace maia

#endif // LBCELLPROPERTIES_H_
