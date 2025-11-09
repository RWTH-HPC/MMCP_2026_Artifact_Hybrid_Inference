// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef TYPETRAITS_H
#define TYPETRAITS_H

#include "COMM/mpioverride.h"
#include "INCLUDE/maiatypes.h"

namespace maia {

#ifndef MAIA_MS_COMPILER
// Disallow usage of non-specified traits
namespace detail_ {
template <class T>
struct traits_error {};
} // namespace detail_
template <class T>
struct type_traits {
  using error = typename detail_::traits_error<T>::ERROR_CANNOT_USE_NON_SPECIFIED_TRAITS;
};
#else
template <class T>
struct type_traits {};

template <>
struct type_traits<long> {
  static MString name() { return "MLong"; };
  static MPI_Datatype mpiType() { return MPI_LONG; };
};
#endif

// Specialized templates for known types

template <>
struct type_traits<MFloat> {
  static MString name() { return "MFloat"; };
  static MPI_Datatype mpiType() { return MPI_DOUBLE; };
};

template <>
struct type_traits<MString> {
  static MString name() { return "MString"; };
  // MString has no corresponding MPI data type
};

template <>
struct type_traits<MChar> {
  static MString name() { return "MChar"; };
  static MPI_Datatype mpiType() { return MPI_CHAR; };
};

template <>
struct type_traits<MInt> {
  static MString name() { return "MInt"; };
  static MPI_Datatype mpiType() { return MPI_INT; };
};

template <>
struct type_traits<MUint> {
  static MString name() { return "MUint"; };
  static MPI_Datatype mpiType() { return MPI_UNSIGNED; };
};

template <>
struct type_traits<MLong> {
  static MString name() { return "MLong"; };
  static MPI_Datatype mpiType() { return MPI_LONG; };
};

template <>
struct type_traits<MUlong> {
  static MString name() { return "MUlong"; };
  static MPI_Datatype mpiType() { return MPI_UNSIGNED_LONG; };
};

template <>
struct type_traits<MBool> {
  static MString name() { return "MBool"; };
  static MPI_Datatype mpiType() { return MPI_C_BOOL; };
};

template <>
struct type_traits<uint_fast8_t> {
  static MString name() { return "uint_fast8_t"; };
  static MPI_Datatype mpiType() { return MPI_UINT8_T; };
};

} /* namespace maia */


#endif /* ifndef TYPETRAITS_H */
