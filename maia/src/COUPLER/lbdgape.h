// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef LBDGAPE_H
#define LBDGAPE_H

#include "coupling.h"
#include "couplingdgape.h"

template <MInt nDim, MInt nDist, class SysEqn>
class LbDgApe final : public CouplingDgApe<nDim, CouplingLB<nDim, nDist, SysEqn>> {
 public:
  // Typedefs
  using Base = Coupling;
  // Typedefs (LB)
  using BaseLb = CouplingLB<nDim, nDist, SysEqn>;
  using LbSolver = typename BaseLb::solverType;
  // Typedefs (DG)
  using SysEqnDg = DgSysEqnAcousticPerturb<nDim>;
  using BaseDg = CouplingDgApe<nDim, BaseLb>;
  using DgSolver = typename BaseDg::DgCartesianSolverType;

  using BaseLb::lbSolver;

  //--Ctor & Dtor
  LbDgApe(const MInt couplingId, LbSolver* const lb, DgSolver* const dg)
    : Base(couplingId), BaseDg(couplingId, dg, lb) {}
  ~LbDgApe(){};

  //--Methods
 public:
  // Main coupling functions
  void init() override {
    initConversionFactors();
    BaseDg::init();
  };
  // Virtual methods to override
  void finalizeCouplerInit() override;

 private:
  void initConversionFactors();
  // Virtual methods to override
  LbSolver& donorSolver(const MInt xSolverId = 0) const override { return lbSolver(xSolverId); };
  void getDonorVelocityAndVorticity(const std::vector<MInt>& donorCellIds, MFloatScratchSpace& p_velocity,
                                    MFloatScratchSpace& p_vorticity) override;

  virtual void performUnitConversion(const MString& /*name*/, const MInt /*count*/, const MInt /*stride*/,
                                     MFloat* /*data*/) override;
  void calcSourceLambLinearized(const MFloat* const /*velocity*/, const MFloat* const /*vorticity*/,
                                MFloat* /*sourceTerms*/) override {
    mTerm(1, AT_, "Error: calcSourceLambLinearized source term not implemented, yet, for lbdgape!");
  };
  void calcSourceLambNonlinear(const MFloat* const velocity, const MFloat* const vorticity,
                               MFloat* const sourceTerms) override;
  void calcSourceQmII(const MFloat* const /*velocity*/, MFloat* const /*sourceTerms*/) override {
    mTerm(1, AT_, "Error: calcSourceQmII source term not implemented, yet, for lbdgape!");
  };
  void calcSourceQmIII(const MFloat* const /*velocity*/, MFloat* /*sourceTerms*/) override {
    mTerm(1, AT_, "Error: calcSourceQmIII source term not implemented, yet, for lbdgape!");
  };
  void calcSourceQe(const MFloat* const /*velocity*/, const MFloat /*time*/, MFloat* const /*sourceTerms*/) override {
    mTerm(1, AT_, "Error: calcSourceQe source term not implemented, yet, for lbdgape!");
  };
  void calcSourceQc(const MFloat* const /*velocity*/, MFloat* const /*sourceTerms*/, const MFloat /*time*/,
                    const MInt /*timeStep*/) override {
    mTerm(1, AT_, "Error: calcSourceQc source term not implemented, yet, for lbdgape!");
  };

  //---Variables----------------------------------------------------------------
 private:
  // Conversion-Factors
  struct ConversionFactors {
    MFloat time{};
    MFloat length{};
    MFloat density{};
    MFloat velocity{};
    MFloat vorticity{};
    MFloat lamb{};
  } m_conversionLb2Dg;
};
#endif // LBDGAPE_H
