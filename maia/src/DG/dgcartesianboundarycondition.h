// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef DGBOUNDARYCONDITION_H_
#define DGBOUNDARYCONDITION_H_

#include "INCLUDE/maiatypes.h"
#include "dgcartesianelementcollector.h"
#include "dgcartesianhelementcollector.h"
#include "dgcartesianinterpolation.h"
#include "dgcartesiansolver.h"
#include "dgcartesiansurfacecollector.h"

template <MInt nDim, class SysEqn>
class DgCartesianSolver;

template <MInt nDim, class SysEqn_>
class DgBoundaryCondition {
  // Typedefs
 public:
  using SysEqn = SysEqn_;
  using SolverType = DgCartesianSolver<nDim, SysEqn>;
  using ElementCollector = maia::dg::collector::ElementCollector<nDim, SysEqn>;
  using HElementCollector = maia::dg::collector::HElementCollector<nDim, SysEqn>;
  using SurfaceCollector = maia::dg::collector::SurfaceCollector<nDim, SysEqn>;

  // Interface methods
 public:
  /// Destructor must be virtual
  virtual ~DgBoundaryCondition() {}

  /// Init method to initialize boundary condition for range of surfaces
  void init(const MInt begin_, const MInt end_) {
    m_begin = begin_;
    m_end = end_;
    init();
  }

  /// Apply method to apply boundary condition
  virtual void apply(const MFloat time) = 0;

  /// Returns name of boundary condition
  virtual MString name() const = 0;

  /// Constructor saves arguments to member variables
  DgBoundaryCondition(SolverType& solver_, MInt bcId) : m_solver(solver_), m_bcId(bcId) {}

  /// Return boundary condition if of this boundary condition
  MInt id() const { return m_bcId; }

  /// Return index of first surface
  MInt begin() const { return m_begin; }

  /// Return index of one-past-last surface
  MInt end() const { return m_end; }

  /// Return number of boundary surfaces
  MInt count() const { return end() - begin(); }

  /// Return number of restart variables that need to be stored/loaded. If a
  /// derived class does not need additional restart variables, it may omit an
  /// implementation.
  virtual MInt noRestartVars() const { return 0; }

  /// Return local number of nodes
  virtual MInt getLocalNoNodes() const { return 0; }

  /// Return name of restart variable
  virtual MString restartVarName(const MInt NotUsed(id)) const { TERMM(1, "Not implemented in derived class!"); }

  /// Copy restart variable data from pointer to boundary condition class
  virtual void setRestartVariable(const MInt NotUsed(id), const MFloat* const NotUsed(data)) {
    TERMM(1, "Not implemented in derived class!");
  }

  /// Copy restart variable data from boundary condition class to pointer
  virtual void getRestartVariable(const MInt NotUsed(id), MFloat* const NotUsed(data)) const {
    TERMM(1, "Not implemented in derived class!");
  }

  virtual MInt noBcElements() const { return 0; }
  virtual MBool hasBcElement(const MInt NotUsed(elementId)) const { return false; }

  // DLB methods
  virtual MInt noCellDataDlb() const { return 0; }
  virtual MInt cellDataTypeDlb(const MInt NotUsed(dataId)) const { return -1; }
  virtual MInt cellDataSizeDlb(const MInt NotUsed(dataId), const MInt NotUsed(cellId)) const { return 0; }
  virtual void getCellDataDlb(const MInt NotUsed(dataId), MFloat* const NotUsed(data)) const {
    TERMM(1, "Not implemented in derived class!");
  }
  virtual void getCellDataDlb(const MInt NotUsed(dataId), MInt* const NotUsed(data)) const {
    TERMM(1, "Not implemented in derived class!");
  }
  virtual void setCellDataDlb(const MInt NotUsed(dataId), const MFloat* const NotUsed(data)) {
    TERMM(1, "Not implemented in derived class!");
  }
  virtual void setCellDataDlb(const MInt NotUsed(dataId), const MInt* const NotUsed(data)) {
    TERMM(1, "Not implemented in derived class!");
  }

  // Methods for derived classes
 protected:
  /// Return reference to solver
  SolverType& solver() { return m_solver; }

  /// Return reference to SysEqn object
  SysEqn& sysEqn() { return m_solver.m_sysEqn; }

  /// Return pointer to surface flux
  MFloat* flux(const MInt i) { return &surfaces().flux(i); }

  /// Return reference to elements
  ElementCollector& elements() { return m_solver.m_elements; }

  /// Return element id corresponding to given cell id
  MInt getElementByCellId(const MInt cellId) const { return m_solver.m_elements.getElementByCellId(cellId); }

  /// Return reference to h-elements
  HElementCollector& helements() { return m_solver.m_helements; }

  /// Return reference to surfaces
  SurfaceCollector& surfaces() { return m_solver.m_surfaces; }

  /// Return true if surface is a MPI surface
  MBool isMpiSurface(const MInt id_) const { return m_solver.isMpiSurface(id_); };

  /// Return if h-element is needed for given cell
  MBool needHElementForCell(const MInt cellId) { return m_solver.needHElementForCell(cellId); };

  /// Return h-element id for an element
  MInt getHElementId(const MInt elementId) { return m_solver.getHElementId(elementId); };

  /// Return interpolation
  const DgInterpolation& interpolation(const MInt polyDeg, const MInt noNodes1D) const {
    return m_solver.m_interpolation[polyDeg][noNodes1D];
  }

  /// Return integration method
  MInt integrationMethod() const { return m_solver.m_dgIntegrationMethod; }

  /// Return time integration scheme
  MInt timeIntegrationScheme() const { return m_solver.m_dgTimeIntegrationScheme; }

  /// Return maximum polynomial degree
  MInt maxPolyDeg() const { return m_solver.m_maxPolyDeg; }

  /// Return maximum number of nodes
  MInt maxNoNodes1D() const { return m_solver.m_maxNoNodes1D; }

  /// Return current time step size
  MFloat dt() const { return m_solver.m_dt; }

  /// Return if a restart is performed
  MBool isRestart() const { return m_solver.m_restart; }

  /// Access to time integration method.
  void subTimeStepRk(const MFloat dt_, const MInt stage, const MInt totalSize, const MFloat* const rhs,
                     MFloat* const variables, MFloat* const timeIntStorage) {
    m_solver.subTimeStepRk(dt_, stage, totalSize, rhs, variables, timeIntStorage);
  }

  // Access to DG-operator methods

  template <class F>
  void calcVolumeIntegral(const MInt noElements, ElementCollector& elem, F& fluxFct) {
    m_solver.calcVolumeIntegral(noElements, elem, fluxFct);
  }

  void resetBuffer(const MInt totalSize, MFloat* const buffer) { m_solver.resetBuffer(totalSize, buffer); }

  void applyJacobian(const MInt noElements, ElementCollector& elem) { m_solver.applyJacobian(noElements, elem); }

  template <class F>
  void calcSourceTerms(const MFloat t, const MInt noElements, ElementCollector& elem, F& sourceFct) {
    m_solver.calcSourceTerms(t, noElements, elem, sourceFct);
  }

  void calcSurfaceIntegral(const MInt begin_, const MInt end_, ElementCollector& elem, SurfaceCollector& surf,
                           HElementCollector& helem, const MInt noHElements) {
    m_solver.calcSurfaceIntegral(begin_, end_, elem, surf, helem, noHElements);
  }

  template <class F>
  void calcRegularSurfaceFlux(const MInt begin_, const MInt end_, SurfaceCollector& surf, F& riemannFct) {
    m_solver.calcRegularSurfaceFlux(begin_, end_, surf, riemannFct);
  }

  // Methods that should not be called from anywhere else
 private:
  /// Init method that may be implemented by the derived classes. This is called
  /// by the init(MInt, MInt) method after the begin/end indices are saved and
  /// accessible by begin()/end(). If a derived class does not need to
  /// initialize anything, it just does not provide an implementation.
  virtual void init() {}

 private:
  /// Store a reference to the solver
  SolverType& m_solver;

  /// The boundary condition id of this boundary condition
  const MInt m_bcId;

  /// Id of first surface of this boundary condition
  MInt m_begin;

  /// Id of one-past-last surface of this boundary condition
  MInt m_end;
};


namespace maia {
namespace dg {
namespace bc {

template <class Functor, class Class, class IdType, class... Args>
void loop(Functor&& fun, Class* object, const IdType begin, const IdType end, Args... args) {
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(IdType i = begin; i < end; i++) {
    (object->*fun)(i, std::forward<Args>(args)...);
  }
}

} // namespace bc
} // namespace dg
} // namespace maia

#endif // DGBOUNDARYCONDITION_H_
