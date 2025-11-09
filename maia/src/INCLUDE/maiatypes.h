// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only


// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only


// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only


// Copyright (C) 2019 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier:    LGPL-3.0-only


#ifndef TYPES_H
#define TYPES_H
////////////////////////////////////////////////////////////////////////////////
/// \file \brief This files defines the types used in MAIA
////////////////////////////////////////////////////////////////////////////////
#include <cassert>
#include <cstdint>
#include <string>
////////////////////////////////////////////////////////////////////////////////
// Note: If you change a type, make sure to check if typetraits.h needs to be
//       updated as well.
using MFloat = double;
using MLongFloat = long double;

using MString = std::basic_string<char>;
using MChar = char;
using MUchar = unsigned char;
using MBool = bool;

using MShort = int_least16_t;
using MUshort = uint_least16_t;
using MInt = int32_t;
using MUint = uint32_t;
using MLong = int64_t;
using MUlong = uint64_t;

////////////////////////////////////////////////////////////////////////////////
//! define array structures
using MZone = struct z {
  MString name;
  MInt id;
  MInt* solvers;
  MInt noSolvers;
};
////////////////////////////////////////////////////////////////////////////////
template <class T, MInt bits, MBool toggled = false>
class MTXbit {
  static_assert(bits <= sizeof(T) * 8);
  T storage = toggled ? (T)~0 : 0;

 public:
  using type = T;
  static constexpr const MInt MAX = (1 << bits) - 1;
  constexpr MTXbit(T val) noexcept : storage(val) {}
  constexpr MTXbit() = default;
  void set(const MInt index, const MUint value) {
    assert((int)sizeof(T) * 8 / bits > index);
    assert(value <= MAX);
    storage = (storage & (T) ~(MAX << bits * index)) | (value << bits * index);
  }
  MInt get(const MInt index) const {
    assert((int)sizeof(T) * 8 / bits > index);
    return (storage >> (bits * index)) & MAX;
  }
  MBool any() const { return storage != 0; }
  MBool all() const { return storage == static_cast<T>(~0); }
  MTXbit& operator&=(const MTXbit& rhs) {
    this->storage &= rhs.storage;
    return *this;
  }
  MTXbit& operator|=(const MTXbit& rhs) {
    this->storage |= rhs.storage;
    return *this;
  }
  T data() const { return storage; }
};
template <MBool toggled = false>
using M32X4bit = MTXbit<uint_fast32_t, 4, toggled>; // save 8 values [0...15] in 4bytes
template <MBool toggled = false>
using M16X2bit = MTXbit<uint_fast16_t, 2, toggled>; // save 8 values [0...3] in 2bytes
template <MBool toggled = false>
using M8X1bit = MTXbit<uint_fast8_t, 1, toggled>; // save 8 bools in 1byte
#endif
