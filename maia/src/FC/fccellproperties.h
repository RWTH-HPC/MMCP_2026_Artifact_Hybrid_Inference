// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef MAIAFCCELLPROPERTIES_H_
#define MAIAFCCELLPROPERTIES_H_

#include <bitset>
#include <type_traits>

/// FC cell Property Labels.
enum class FcCell {
  IsHalo,
  IsWindow,
  IsBndryCell,
  IsActive,
  WasActive,
  NeedsSubCells,
  // <<< add new properties here
  // ---
  NumProperties
};


namespace maia {
namespace fc {
namespace cell {

/// Converts property name to underlying integer value
constexpr std::underlying_type<FcCell>::type p(const FcCell property) {
  return static_cast<std::underlying_type<FcCell>::type>(property);
}

using BitsetType = std::bitset<p(FcCell::NumProperties)>;

} // namespace cell
} // namespace fc
} // namespace maia

#endif // MAIAFCCELLPROPERTIES_H_
