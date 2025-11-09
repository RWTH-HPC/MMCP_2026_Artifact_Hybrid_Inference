// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "oneDFlame.h"
#include "MEMORY/alloc.h"

#if defined(WITH_CANTERA)
using PtrCanteraSolution = typename std::shared_ptr<Cantera::Solution>;
using PtrCanteraThermo = typename std::shared_ptr<Cantera::ThermoPhase>;
using PtrCanteraKinetics = typename std::shared_ptr<Cantera::Kinetics>;
using PtrCanteraTransport = typename std::shared_ptr<Cantera::Transport>;
#endif

OneDFlame::OneDFlame(MFloat T, MFloat p, MFloat* Y, MInt domainId, MInt solverId, std::vector<MString> speciesName)
  : m_TInlet(T),
    m_pInlet(p),
    m_YInlet(Y),
    m_domainId(domainId),
    m_solverId(solverId),
    m_speciesName(std::move(speciesName)) {
// These lines allow Intel compilation when Cantera is not defined
#if not defined(WITH_CANTERA)
  (void)m_TInlet;
  (void)m_pInlet;
  (void)m_YInlet;
  (void)m_rhoInlet;
#endif

  readProperties();
}

OneDFlame::~OneDFlame() = default;

void OneDFlame::run() {
#if defined(WITH_CANTERA)

  const MUint noSpecies = m_canteraThermo->nSpecies();

  const MFloat UInlet = 0.1;

  // MFloat phi = m_equivalenceRatio;
  std::vector<MFloat> X(noSpecies, 0.0);

  // given as molar fraction
  // gas->setEquivalenceRatio(phi, "H2", "O2:0.21,N2:0.79");

  std::vector<MFloat> YInlet(noSpecies, F0);
  for(MUint s = 0; s < noSpecies; ++s) {
    YInlet[s] = m_YInlet[s];
  }

  // set gas to correct thermodynamic state
  m_canteraThermo->setState_TPY(m_TInlet, m_pInlet, &YInlet[0]);
  m_canteraThermo->getMoleFractions(X.data());

  m_rhoInlet = m_canteraThermo->density();

  m_canteraThermo->equilibrate("HP");
  std::vector<MFloat> YOutlet(noSpecies);
  m_canteraThermo->getMassFractions(&YOutlet[0]);

  MFloat rhoOutlet = m_canteraThermo->density();
  MFloat adiabaticT = m_canteraThermo->temperature();

  //-------- step 1: create the flow -------------
  Cantera::StFlow flow(m_canteraThermo);
  flow.setFreeFlow();

  // create an initial grid
  MInt nz = m_noGridPoints;
  MFloat lz = m_domainLength;
  std::vector<MFloat> z(nz);
  MFloat dz = lz / ((MFloat)(nz - 1));
  for(int iz = 0; iz < nz; iz++) {
    z[iz] = ((MFloat)iz) * dz;
  }

  flow.setupGrid(nz, &z[0]);

  std::unique_ptr<Cantera::Transport> trmix(
      Cantera::newTransportMgr(m_transportModel, m_canteraSolution->thermo().get()));

  flow.setTransport(*trmix);
  flow.setKinetics(*m_canteraSolution->kinetics());
  flow.setPressure(m_pInlet);

  //------- step 2: create the inlet  -----------------------
  Cantera::Inlet1D inlet;

  inlet.setMoleFractions(X.data());
  double massFlow = UInlet * m_rhoInlet;
  inlet.setMdot(massFlow);
  inlet.setTemperature(m_TInlet);

  //------- step 3: create the outlet  ---------------------
  Cantera::Outlet1D outlet;

  //=================== create the container and insert the domains =====
  std::vector<Cantera::Domain1D*> domains{&inlet, &flow, &outlet};
  Cantera::Sim1D flame(domains);

  //----------- Supply initial guess----------------------
  std::vector<MFloat> locs{0.0, 0.3, 0.7, 1.0};
  std::vector<MFloat> value;

  const double UOutlet = inlet.mdot() / rhoOutlet;
  value = {UInlet, UInlet, UOutlet, UOutlet};
  flame.setInitialGuess("velocity", locs, value);
  value = {m_TInlet, m_TInlet, adiabaticT, adiabaticT};
  flame.setInitialGuess("T", locs, value);

  for(size_t i = 0; i < noSpecies; i++) {
    value = {m_YInlet[i], m_YInlet[i], YOutlet[i], YOutlet[i]};
    flame.setInitialGuess(m_canteraThermo->speciesName(i), locs, value);
  }

  inlet.setMoleFractions(X.data());
  inlet.setMdot(massFlow);
  inlet.setTemperature(m_TInlet);

  MInt flowdomain = 1;
  MFloat gridRatio = m_gridRatio; // 10
  MFloat gridSlope = m_gridSlope; // 0.1
  MFloat gridCurve = m_gridCurve; // 0.01
  MFloat gridPrune = m_gridPrune;
  MInt loglevel = (m_domainId == 0) ? 1 : 0;

  flame.setRefineCriteria(flowdomain, gridRatio, gridSlope, gridCurve, gridPrune);

  flame.setMaxTimeStepCount(10000);
  flame.setMaxGridPoints(flowdomain, 10000);

  flame.setFixedTemperature(0.5 * (m_TInlet + adiabaticT));

  flame.solve(loglevel, false);
  flow.solveEnergyEqn();
  flame.solve(loglevel, true);

  MFloat flameSpeedNew = flame.value(flowdomain, flow.componentIndex("velocity"), 0);
  MFloat flameSpeedOld = 0.0;

  gridRatio = 10.0;
  MInt noIterations = 10;

  // iterative 1D grid convergence
  for(MInt iteration = 0; iteration < noIterations; ++iteration) {
    gridSlope *= 0.5;
    gridCurve *= 0.35;

    flame.setRefineCriteria(flowdomain, gridRatio, gridSlope, gridCurve, gridPrune);

    flame.solve(loglevel, true);

    flameSpeedOld = flameSpeedNew;

    flameSpeedNew = flame.value(flowdomain, flow.componentIndex("velocity"), 0);

    MFloat relError = fabs(flameSpeedNew - flameSpeedOld) / fabs(flameSpeedNew);

    if(relError < m_allowedRelError) break;
  }

  MInt noGridPoints = flow.nPoints();

  // probably a better way
  mAlloc(m_profile.velocity, noGridPoints, "velocity", -9999.9, AT_);
  mAlloc(m_profile.temperature, noGridPoints, "temperature", -9999.9, AT_);
  mAlloc(m_profile.massFractions, noSpecies * noGridPoints, "massFractions", -9999.9, AT_);

  for(MInt n = 0; n < noGridPoints; ++n) {
    m_profile.velocity[n] = flame.value(flowdomain, flow.componentIndex("velocity"), n);
    m_profile.temperature[n] = flame.value(flowdomain, flow.componentIndex("T"), n);
    m_grid.push_back(flow.grid(n));

    for(MUint s = 0; s < noSpecies; ++s) {
      m_profile.massFractions[noSpecies * n + s] = flame.value(flowdomain, flow.componentIndex(m_speciesName[s]), n);
    }
  }

  ASSERT((MInt)m_grid.size() == noGridPoints,
         "Vector m_oneGridSize has incorrect number of values, terminating now...");

  std::vector<MFloat> dTdx(noGridPoints, 0.0);
  MFloat maximumSlope = 0.0;
  for(MInt n = 1; n < (noGridPoints - 1); ++n) {
    MFloat deltaX = m_grid[n + 1] - m_grid[n - 1];
    MFloat deltaT = m_profile.temperature[n + 1] - m_profile.temperature[n - 1];
    dTdx[n] = deltaT / deltaX;
    if(dTdx[n] > maximumSlope) maximumSlope = dTdx[n];
  }

  m_laminarFlameThickness = (m_profile.temperature[noGridPoints - 1] - m_profile.temperature[0]) / maximumSlope;

  m_fixedTemperatureLocation = flame.fixedTemperatureLocation();
#endif
}

void OneDFlame::log(MFloat cellLengthAtMaxRfnLevel) {
  const MFloat requiredFlameResolution = 10.0;
  if(m_domainId == 0) {
    std::cerr << "Computed flame thickness: " << m_laminarFlameThickness << std::endl;
    std::cerr << "Resolution at maximum refinement level: " << cellLengthAtMaxRfnLevel << std::endl;
    std::cerr << "Flame front resolved by a total of " << m_laminarFlameThickness / cellLengthAtMaxRfnLevel << " cells"
              << std::endl;
    if(m_laminarFlameThickness / cellLengthAtMaxRfnLevel < requiredFlameResolution) {
      std::cerr
          << "Flame front is underresolved! Increase the refinement level to get an accurate solution! A resolution "
             "of at least "
          << requiredFlameResolution << "cells inside the flame front is required." << std::endl;
    }
  }
}

#if defined(WITH_CANTERA)
void OneDFlame::setCanteraObjects(PtrCanteraSolution canteraSolution,
                                  PtrCanteraThermo canteraThermo,
                                  PtrCanteraKinetics canteraKinetics,
                                  PtrCanteraTransport canteraTransport) {
  m_canteraSolution = canteraSolution;
  m_canteraThermo = canteraThermo;
  m_canteraKinetics = canteraKinetics;
  m_canteraTransport = canteraTransport;
}
#endif

void OneDFlame::readProperties() {
  // Default grid values through trial and error of succesfully converged one dimensional simulations
  m_gridRatio = 5.0;
  m_gridRatio = Context::getSolverProperty<MFloat>("gridRatio", m_solverId, AT_, &m_gridRatio);

  m_gridSlope = 0.5;
  m_gridSlope = Context::getSolverProperty<MFloat>("gridSlope", m_solverId, AT_, &m_gridSlope);

  m_gridCurve = 0.3;
  m_gridCurve = Context::getSolverProperty<MFloat>("gridCurve", m_solverId, AT_, &m_gridCurve);

  m_gridPrune = -0.00005;
  m_gridPrune = Context::getSolverProperty<MFloat>("gridPrune", m_solverId, AT_, &m_gridPrune);

  m_noGridPoints = 10;
  m_noGridPoints = Context::getSolverProperty<MInt>("noGridPoints", m_solverId, AT_, &m_noGridPoints);

  m_domainLength = 0.03;
  m_domainLength = Context::getSolverProperty<MFloat>("domainLength", m_solverId, AT_, &m_domainLength);

  // Controls the number of iterations of the one dimensional simulation for a "converged" result (based on laminar
  // burning speed)
  m_allowedRelError = 0.001;
  m_allowedRelError = Context::getSolverProperty<MFloat>("allowedRelError", m_solverId, AT_, &m_allowedRelError);

  m_transportModel = Context::getSolverProperty<MString>("transportModel", m_solverId, AT_);
}
