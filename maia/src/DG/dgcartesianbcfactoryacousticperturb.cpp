// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "INCLUDE/maiatypes.h"
#include "dgcartesianboundarycondition.h"
#include "dgcartesianboundaryconditionfactory.h"
#include "dgcartesiansyseqnacousticperturb.h"

// Includes for boundary conditions start here
#include "dgcartesianbcacousticperturbcbc.h"
#include "dgcartesianbcacousticperturbrbc.h"
#include "dgcartesianbcacousticperturbsolidwall.h"
#include "dgcartesianbcacousticperturbstraightductexit.h"
#include "dgcartesianbcexact.h"


namespace maia {
namespace dg {
namespace bc {
namespace factory {

template <MInt nDim>
struct Init<nDim, DgSysEqnAcousticPerturb<nDim>> {
  using SysEqn = DgSysEqnAcousticPerturb<nDim>;
  using FactoryType = DgBoundaryConditionFactory<nDim, SysEqn>;
  static void init(FactoryType& factory);
};

template <MInt nDim>
void Init<nDim, DgSysEqnAcousticPerturb<nDim>>::init(FactoryType& factory) {
  // Import namespace for ease-of-use
  using namespace maia::dg::bc;

  /// \property
  /// \page dgBCApe List of all boundary conditions for solvertype `MAIA_DISCONTINUOUS_GALERKIN` and
  /// `DG_SYSEQN_ACOUSTICPERTURB` List of all available BCs in Init<nDim, DgSysEqnAcousticPerturb<nDim>>::init():

  // Add all boundary conditions with a unique boundary condition id
  /// \dgSwitchCase{dgBCApe, 0, ...}
  factory.add(0, Type<DgBcExact<nDim, SysEqn>>());

  /// \dgSwitchCase{dgBCApe, 20\,200-203, No-slip wall}
  factory.add(20, Type<DgBcAcousticPerturbSolidWall<nDim, SysEqn, false>>());
  factory.add(200, Type<DgBcAcousticPerturbSolidWall<nDim, SysEqn, false>>());
  factory.add(201, Type<DgBcAcousticPerturbSolidWall<nDim, SysEqn, false>>());
  factory.add(202, Type<DgBcAcousticPerturbSolidWall<nDim, SysEqn, false>>());
  factory.add(203, Type<DgBcAcousticPerturbSolidWall<nDim, SysEqn, false>>());

  /// \dgSwitchCase{dgBCApe, 21, Slip wall}
  factory.add(21, Type<DgBcAcousticPerturbSolidWall<nDim, SysEqn, true>>());
  /// \dgSwitchCase{dgBCApe, 22, Outlet of straight duct combined in solid wall}
  factory.add(22, Type<DgBcAcousticPerturbStraightDuctExit<nDim, SysEqn>>());

  /// \dgSwitchCase{dgBCApe, 301\,304, CBC}
  factory.add(301, Type<DgBcAcousticPerturbCBC<nDim>>()); // CBC
  factory.add(304, Type<DgBcAcousticPerturbCBC<nDim>>()); // CBC

  /// \dgSwitchCase{dgBCApe, 400-405, Radiation boundary condition (RBC)}
  factory.add(400, Type<DgBcAcousticPerturbRBC<nDim>>());
  factory.add(401, Type<DgBcAcousticPerturbRBC<nDim>>());
  factory.add(402, Type<DgBcAcousticPerturbRBC<nDim>>());
  factory.add(403, Type<DgBcAcousticPerturbRBC<nDim>>());
  factory.add(404, Type<DgBcAcousticPerturbRBC<nDim>>());
  factory.add(405, Type<DgBcAcousticPerturbRBC<nDim>>());

  /// \dgSwitchCase{dgBCApe, 30-49, boundary condition ids for sponge bcs}
  factory.add(30, Type<DgBcExact<nDim, SysEqn>>());
  factory.add(31, Type<DgBcExact<nDim, SysEqn>>());
  factory.add(32, Type<DgBcExact<nDim, SysEqn>>());
  factory.add(33, Type<DgBcExact<nDim, SysEqn>>());
  factory.add(34, Type<DgBcExact<nDim, SysEqn>>());
  factory.add(35, Type<DgBcExact<nDim, SysEqn>>());
  factory.add(36, Type<DgBcExact<nDim, SysEqn>>());
  factory.add(37, Type<DgBcExact<nDim, SysEqn>>());
  factory.add(38, Type<DgBcExact<nDim, SysEqn>>());
  factory.add(39, Type<DgBcExact<nDim, SysEqn>>());
  factory.add(40, Type<DgBcExact<nDim, SysEqn>>());
  factory.add(41, Type<DgBcExact<nDim, SysEqn>>());
  factory.add(42, Type<DgBcExact<nDim, SysEqn>>());
  factory.add(43, Type<DgBcExact<nDim, SysEqn>>());
  factory.add(44, Type<DgBcExact<nDim, SysEqn>>());
  factory.add(45, Type<DgBcExact<nDim, SysEqn>>());
  factory.add(46, Type<DgBcExact<nDim, SysEqn>>());
  factory.add(47, Type<DgBcExact<nDim, SysEqn>>());
  factory.add(48, Type<DgBcExact<nDim, SysEqn>>());
  factory.add(49, Type<DgBcExact<nDim, SysEqn>>());
}

/// Explicit instantiations for 2D and 3D
template struct Init<2, DgSysEqnAcousticPerturb<2>>;
template struct Init<3, DgSysEqnAcousticPerturb<3>>;

} // namespace factory
} // namespace bc
} // namespace dg
} // namespace maia

/// Explicit instantiations for 2D and 3D
template class DgBoundaryConditionFactory<2, DgSysEqnAcousticPerturb<2>>;
template class DgBoundaryConditionFactory<3, DgSysEqnAcousticPerturb<3>>;
