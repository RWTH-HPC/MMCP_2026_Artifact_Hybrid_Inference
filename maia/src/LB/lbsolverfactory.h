// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef LBSOLVERFACTORY_H
#define LBSOLVERFACTORY_H

#include <memory>
#include "solver.h"

template <MInt nDim>
class Geometry;

namespace maia {

namespace grid {
template <MInt nDim>
class Proxy;
}

namespace lb {

template <MInt nDim>
class LbSolverFactory {
 public:
  static std::unique_ptr<Solver> create(const MInt solverId, maia::grid::Proxy<nDim>& gridProxy,
                                        Geometry<nDim>& geometry, const MPI_Comm comm);
};

} // namespace lb
} // namespace maia

#endif // LBSOLVERFACTORY_H
