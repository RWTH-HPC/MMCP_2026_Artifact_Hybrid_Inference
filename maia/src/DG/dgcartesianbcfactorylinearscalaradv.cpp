// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "INCLUDE/maiatypes.h"
#include "dgcartesianboundarycondition.h"
#include "dgcartesianboundaryconditionfactory.h"
#include "dgcartesiansyseqnlinearscalaradv.h"

// Includes for boundary conditions start here
#include "dgcartesianbcexact.h"

namespace maia {
namespace dg {
namespace bc {
namespace factory {

template <MInt nDim>
struct Init<nDim, DgSysEqnLinearScalarAdv<nDim>> {
  using SysEqn = DgSysEqnLinearScalarAdv<nDim>;
  using FactoryType = DgBoundaryConditionFactory<nDim, SysEqn>;
  static void init(FactoryType& factory);
};

template <MInt nDim>
void Init<nDim, DgSysEqnLinearScalarAdv<nDim>>::init(FactoryType& factory) {
  // Import namespace for ease-of-use
  using namespace maia::dg::bc;

  // Add all boundary conditions with a unique boundary condition id
  factory.add(0, Type<DgBcExact<nDim, SysEqn>>());

  // boundary condition ids for sponge bcs
}

/// Explicit instantiations for 2D and 3D
template struct Init<2, DgSysEqnLinearScalarAdv<2>>;
template struct Init<3, DgSysEqnLinearScalarAdv<3>>;

} // namespace factory
} // namespace bc
} // namespace dg
} // namespace maia

/// Explicit instantiations for 2D and 3D
template class DgBoundaryConditionFactory<2, DgSysEqnLinearScalarAdv<2>>;
template class DgBoundaryConditionFactory<3, DgSysEqnLinearScalarAdv<3>>;
