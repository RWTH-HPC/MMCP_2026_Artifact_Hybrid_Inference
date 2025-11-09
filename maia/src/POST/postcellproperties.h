// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef POSTCELLPROPERTIES_H_
#define POSTCELLPROPERTIES_H_

#include <bitset>
#include <type_traits>

/// POST cell Property Labels.
enum class PostCell {
  IsHalo,   ///< cell is halo-cell
  IsWindow, ///< cell is window-cell
  NumProperties
};

namespace maia {
namespace post {
namespace cell {

/// Converts property name to underlying integer value
constexpr std::underlying_type<PostCell>::type p(const PostCell property) {
  return static_cast<std::underlying_type<PostCell>::type>(property);
}

using BitsetType = std::bitset<p(PostCell::NumProperties)>;

} // namespace cell
} // namespace post
} // namespace maia

#endif // FVCELLPROPERTIES_H_
