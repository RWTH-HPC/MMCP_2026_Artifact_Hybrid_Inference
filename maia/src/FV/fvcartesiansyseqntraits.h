// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef FVCARTESIANSYSEQNTRAITS_H
#define FVCARTESIANSYSEQNTRAITS_H

#include "fvcartesiansyseqndetchem.h"
#include "fvcartesiansyseqneegas.h"
#include "fvcartesiansyseqnns.h"
#include "fvcartesiansyseqnrans.h"

// SFINAE tests. Used for constexpr checks if certain PVs exist.

/// \brief Checks if the SysEqn is SysEqnNS
template <class SysEqn>
constexpr bool isNS = std::is_same_v<SysEqn, FvSysEqnNS<2>> || std::is_same_v<SysEqn, FvSysEqnNS<3>>;

/// \brief Checks if the SysEqn is SysEqnRANS
template <class SysEqn>
constexpr bool
    isRANS =
        std::is_same_v<SysEqn,
                       FvSysEqnRANS<2, RANSModelConstants<RANS_SA_DV>>> || std::is_same_v<SysEqn, FvSysEqnRANS<3, RANSModelConstants<RANS_SA_DV>>> || std::is_same_v<SysEqn, FvSysEqnRANS<2, RANSModelConstants<RANS_FS>>> || std::is_same_v<SysEqn, FvSysEqnRANS<3, RANSModelConstants<RANS_FS>>> || std::is_same_v<SysEqn, FvSysEqnRANS<2, RANSModelConstants<RANS_KOMEGA>>> || std::is_same_v<SysEqn, FvSysEqnRANS<3, RANSModelConstants<RANS_KOMEGA>>>;

/// \brief Checks if the SysEqn is SysEqnEEGas
template <class SysEqn>
constexpr bool isEEGas = std::is_same_v<SysEqn, FvSysEqnEEGas<2>> || std::is_same_v<SysEqn, FvSysEqnEEGas<3>>;

/// \brief Checks if the SysEqn is SysEqnDetChem
template <class SysEqn>
constexpr bool isDetChem = std::is_same_v<SysEqn, FvSysEqnDetChem<2>> || std::is_same_v<SysEqn, FvSysEqnDetChem<3>>;

/// \brief Checks if the primitive variable C exists
template <class SysEqn_>
class hasPV_C {
  typedef char one;
  struct two {
    char x[2];
  };
  // the function test returns a char if PV->C exists, two chars else. These functions are never actually called.
  template <class SysEqn>
  static one test(decltype(&SysEqn::PV->C));
  template <class SysEqn>
  static two test(...);

 public:
  enum { value = sizeof(test<SysEqn_>(0)) == sizeof(one) }; // ::value is constexpr true if PV->C exists, else false
};

/// \brief Checks if the primitive variable N exists
template <class SysEqn_>
class hasPV_N {
  typedef char one;
  struct two {
    char x[2];
  };
  // the function test returns a char if PV->N exists, two chars else. These functions are never actually called.
  template <class SysEqn>
  static one test(decltype(&SysEqn::PV->N));
  template <class SysEqn>
  static two test(...);

 public:
  enum { value = sizeof(test<SysEqn_>(0)) == sizeof(one) }; // ::value is constexpr true if PV->N exists, else false
};

template <class SysEqn_>
class hasPV_A {
  typedef char one;
  struct two {
    char x[2];
  };
  // the function test returns a char if PV->A exists, two chars else. These functions are never actually called.
  template <class SysEqn>
  static one test(decltype(&SysEqn::PV->A));
  template <class SysEqn>
  static two test(...);

 public:
  enum { value = sizeof(test<SysEqn_>(0)) == sizeof(one) }; // ::value is constexpr true if PV->A exists, else false
};

template <class SysEqn>
constexpr bool hasE = SysEqn::ConservativeVariables::RHO_E > -1;

#endif // FVCARTESIANSYSEQNTRAITS_H
