// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef LPTCELLPROPERTIES_H_
#define LPTCELLPROPERTIES_H_

#include <bitset>
#include <type_traits>

/// Lpt cell Property Labels.
enum class LptCell { IsHalo, IsWindow, IsValid, RegridTrigger, NumProperties };


namespace maia {
namespace lpt {
namespace cell {

/// Converts property name to underlying integer value
constexpr std::underlying_type<LptCell>::type p(const LptCell property) {
  return static_cast<std::underlying_type<LptCell>::type>(property);
}

using BitsetType = std::bitset<p(LptCell::NumProperties)>;


} // namespace cell
} // namespace lpt
} // namespace maia

#endif // LPTCELLPROPERTIES_H_
