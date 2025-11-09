// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef DGBOUNDARYCONDITIONFACTORY_H_
#define DGBOUNDARYCONDITIONFACTORY_H_

#include <memory>
#include "INCLUDE/maiatypes.h"
#include "factory.h"

template <MInt nDim, class SysEqn>
class DgCartesianSolver;
template <MInt nDim, class SysEqn>
class DgBoundaryCondition;

namespace maia {
namespace dg {
namespace bc {

/// Simple type-to-type mapper for function overloading
template <typename T>
struct Type {
  using type = T;
};

} // namespace bc
} // namespace dg
} // namespace maia


/// \brief Class to handle creation of boundary condition objects.
///
/// \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
/// \date 2014-04-05
///
/// \tparam nDim See DgCartesianSolver for definition.
/// \tparam SysEqn See DgCartesianSolver for definition.
///
/// This is a class following the factory pattern as described in Alexandrescu
/// (2001).
///
/// References:
///   Andrei Alexandrescu (2001): Modern C++ Design.
template <MInt nDim, class SysEqn>
class DgBoundaryConditionFactory {
 public:
  using IdType = MInt;
  using ReturnType = std::unique_ptr<DgBoundaryCondition<nDim, SysEqn>>;
  using SolverType = DgCartesianSolver<nDim, SysEqn>;
  using ProductCreator = std::function<ReturnType(SolverType&, IdType)>;

 public:
  /// Constructor only saves solver reference to member variable.
  explicit DgBoundaryConditionFactory(SolverType& solver) : m_solver(solver) {}

  /// Add new type to factory.
  template <class T>
  MBool add(const IdType id, maia::dg::bc::Type<T>) {
    return m_factory.add(id, [](SolverType& solver, IdType bcId) { return ReturnType(new T(solver, bcId)); });
  }

  /// Create type from identifier.
  ReturnType create(const IdType id) { return m_factory.create(id, m_solver, id); }

 private:
  SolverType& m_solver;
  MFactory<DgBoundaryCondition<nDim, SysEqn>, IdType, ReturnType, ProductCreator, SolverType&, IdType> m_factory;
};

namespace maia {
namespace dg {
namespace bc {
namespace factory {

/// Declare (but not define) a template class factory initializer.
///
/// This allows it to be seen by the solver without instantiating it. The actual
/// implementation is in the ...factory<SysEqn>.cpp files.
/// We need this quirky trick in order to be able to partially specialize it for
/// the system of equations, which is only possible for class templates, not
/// function templates.
template <MInt nDim, class SysEqn>
struct Init {
  using FactoryType = DgBoundaryConditionFactory<nDim, SysEqn>;
  static void init(FactoryType& factory);
};

} // namespace factory
} // namespace bc
} // namespace dg
} // namespace maia

#endif // DGBOUNDARYCONDITIONFACTORY_H_
