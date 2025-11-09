// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef OneDFlame_H_
#define OneDFlame_H_

#include "INCLUDE/maiaconstants.h"
#include "INCLUDE/maiatypes.h"
#include "MEMORY/alloc.h"
#include "filter.h"

#if defined(WITH_CANTERA)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wfloat-equal"
#include "cantera/base/Solution.h"
#include "cantera/base/global.h"
#include "cantera/oneD/Boundary1D.h"
#include "cantera/oneD/Sim1D.h"
#include "cantera/oneD/StFlow.h"
#include "cantera/thermo.h"
#include "cantera/thermo/IdealGasPhase.h"
#include "cantera/thermo/Species.h"
#include "cantera/thermo/SpeciesThermoInterpType.h"
#include "cantera/transport.h"
#include "cantera/zerodim.h"
#pragma GCC diagnostic pop

using namespace Cantera;

using PtrCanteraSolution = typename std::shared_ptr<Cantera::Solution>;
using PtrCanteraThermo = typename std::shared_ptr<Cantera::ThermoPhase>;
using PtrCanteraKinetics = typename std::shared_ptr<Cantera::Kinetics>;
using PtrCanteraTransport = typename std::shared_ptr<Cantera::Transport>;
#endif

class OneDFlame {
 public:
  OneDFlame(MFloat, MFloat, MFloat*, MInt, MInt, std::vector<MString>);
  ~OneDFlame();

  void run();
  void log(MFloat);

#if defined(WITH_CANTERA)
  void setCanteraObjects(PtrCanteraSolution, PtrCanteraThermo, PtrCanteraKinetics, PtrCanteraTransport);
#endif

  void readProperties();

  std::vector<MFloat> m_grid;

  MFloat m_fixedTemperatureLocation = F0;
  MFloat m_laminarFlameThickness = F0;

  struct {
    MFloat* velocity = nullptr;
    MFloat* temperature = nullptr;
    MFloat* massFractions = nullptr;
  } m_profile;

 private:
#if defined(WITH_CANTERA)
  PtrCanteraSolution m_canteraSolution;
  PtrCanteraThermo m_canteraThermo;
  PtrCanteraKinetics m_canteraKinetics;
  PtrCanteraTransport m_canteraTransport;
#endif

  MFloat m_rhoInlet = F0;
  const MFloat m_TInlet;
  const MFloat m_pInlet;
  const MFloat* m_YInlet;
  const MInt m_domainId;
  const MInt m_solverId;

  MFloat m_gridRatio = F0;
  MFloat m_gridSlope = F0;
  MFloat m_gridCurve = F0;
  MFloat m_gridPrune = F0;
  MInt m_noGridPoints = F0;
  MFloat m_domainLength = F0;

  MFloat m_allowedRelError = F0;

  const std::vector<MString> m_speciesName;

  MString m_transportModel = "";
};
#endif