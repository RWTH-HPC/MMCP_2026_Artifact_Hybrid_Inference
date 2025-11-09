// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef LSCELLPROPERTIES_H_
#define LSCELLPROPERTIES_H_

#include <bitset>
#include <type_traits>

/// LS cell Property Labels.
enum class LsCell { RegridTrigger, NearGap, IsHalo, IsWindow, IsBndryG, NumProperties };

/// LS set Property Labels.
//  (for each set!)
enum class LsSet { InBand, IsGBndryCell, IsGZero, WasGZero, HasPositiveSign, NumSetProperties };


namespace maia {
namespace ls {
namespace cell {

/// Converts property name to underlying integer value
constexpr std::underlying_type<LsCell>::type p(const LsCell property) {
  return static_cast<std::underlying_type<LsCell>::type>(property);
}

using BitsetType = std::bitset<p(LsCell::NumProperties)>;


/// Converts property name to underlying integer value
constexpr std::underlying_type<LsSet>::type setP(const LsSet property) {
  return static_cast<std::underlying_type<LsSet>::type>(property);
}

using BitsetTypeSet = std::bitset<setP(LsSet::NumSetProperties)>;


} // namespace cell
} // namespace ls
} // namespace maia

#endif // LSCELLPROPERTIES_H_
