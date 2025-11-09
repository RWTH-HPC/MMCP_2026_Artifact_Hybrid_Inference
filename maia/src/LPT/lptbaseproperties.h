// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef LPTBASEPROPERTIES_H_
#define LPTBASEPROPERTIES_H_

#include <bitset>
#include <type_traits>

/// Lpt particle Property Labels.
enum class LptBaseProperty {
  IsWindow,        // particle is inside a window-cell
  ReqSend,         // particle is inside a halo-cell and needs to be send to the neighborDomain!
  WasSend,         // particle has just been send to the new window-cell
  ReqBroadcast,    // particle is ouside the current-domain boundingBox
                   // and needs to be communicated to a domain which is not a neighborDomain!
  IsInvalid,       // particle is inside an invalid cell
  HasCollided,     // particle had a collision
  FirstStep,       // particle first time step, not fully living yet...
  ToBeDeleted,     // particle can be deleted
  ToBeRespawn,     // particle can be respawn at different location/domain
  FullyEvaporated, // particle just fully evaporated
  HadWallColl,     // particle undergoing wall-collision in this TS
  NumProperties
};


namespace maia {
namespace lpt {
namespace baseProperty {

/// Converts property name to underlying integer value
constexpr std::underlying_type<LptBaseProperty>::type p(const LptBaseProperty property) {
  return static_cast<std::underlying_type<LptBaseProperty>::type>(property);
}

using BitsetType = std::bitset<p(LptBaseProperty::NumProperties)>;


} // namespace baseProperty
} // namespace lpt
} // namespace maia

#endif // LPTBASEPROPERTIES_H_
