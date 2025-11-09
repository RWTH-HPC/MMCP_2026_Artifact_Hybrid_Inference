// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef COUPLINGLPT_H_
#define COUPLINGLPT_H_

#include "LPT/lpt.h"
#include "coupling.h"
#include "couplingutils.h"
#include "solver.h"

// Forward declarations
template <MInt nDim>
class CouplingParticle;

template <MInt nDim, class CouplingFlowSolver>
class CouplingLpt : public CouplingParticle<nDim>, public CouplingFlowSolver {
 public:
  // Typedefs
  using FlowSolver = typename CouplingFlowSolver::solverType;
  using Lpt = typename CouplingParticle<nDim>::solverType;
  using CouplingParticle<nDim>::lpt;

  // Conversion factor type
  struct ConversionFactor {
    MFloat velocity{};
    MFloat velocitySlope{};
    MFloat density{};
    MFloat pressure{};
    MFloat length{};
    MFloat temperature{};
    MFloat time{};
    MFloat force{};
    MFloat energy{};
    MFloat mass{};
    MFloat viscosity{};
    MFloat momentum{};
  };

  ConversionFactor conversionFlowLpt;
  ConversionFactor conversionLptFlow;

  // Constructor & Destructor
  CouplingLpt(const MInt couplingId, Lpt* particle, FlowSolver* flowSolver)
    : Coupling(couplingId), CouplingParticle<nDim>(couplingId, particle), CouplingFlowSolver(couplingId, flowSolver){};
  virtual ~CouplingLpt() = default;
};

#endif // ifndef COUPLINGLPT_H_
