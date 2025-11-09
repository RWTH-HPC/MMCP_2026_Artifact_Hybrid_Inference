// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef STRUCTUREDCELL
#define STRUCTUREDCELL

// contains the structs to store cell information and surface information
class StructuredCell {
 public:
  StructuredCell() = default;
  ~StructuredCell() = default;
  MFloat** coordinates = nullptr;
  MFloat** mgOldCoordinates = nullptr;
  MFloat** variables = nullptr;
  MFloat** pvariables = nullptr;
  MFloat** oldVariables = nullptr;
  MFloat* temperature = nullptr;

  MFloat** rightHandSide = nullptr;

  MFloat** eFlux = nullptr;
  MFloat** fFlux = nullptr;
  MFloat** gFlux = nullptr;
  MFloat** viscousFlux = nullptr;
  MFloat** flux = nullptr; // contains the convective flux over the surface
  MFloat** dxt = nullptr;  // volume fluxes for the three cell surfaces

  MFloat** dT = nullptr;

  MFloat** dss = nullptr;
  MFloat** ql = nullptr;
  MFloat** qr = nullptr;

  MFloat* cellJac = nullptr;
  MFloat* oldCellJac = nullptr;
  MFloat* cornerJac = nullptr;
  MFloat* surfJac = nullptr;

  MFloat** cellMetrics = nullptr;
  MFloat** cornerMetrics = nullptr;
  MFloat** surfaceMetrics = nullptr;
  MFloat** surfaceMetricsSingularity = nullptr;
  MFloat** surfaceDist = nullptr;
  MFloat** cellLength = nullptr;

  MFloat* localTimeStep = nullptr;

  // fq field
  MFloat** fq = nullptr;
  MFloat** stg_fq = nullptr;

  // spongeLayer
  MFloat* spongeFactor = nullptr;

  // auxillary data
  MFloat* cf = nullptr;
  MFloat* cp = nullptr;
  MFloat* powerVisc = nullptr;
  MFloat* powerPres = nullptr;
  MInt* cfOffsets = nullptr;
  MInt* cpOffsets = nullptr;
  MInt* powerOffsets = nullptr;

  // least squares
  MFloat** reconstructionConstants = nullptr;
  MInt* nghbr = nullptr;
  MInt* numOfNghbr = nullptr;

  // RANS
  MFloat** saFlux1 = nullptr;
  MFloat** saFlux2 = nullptr;
  MFloat* prodDest = nullptr;
  MFloat* P_keps = nullptr;
  MFloat* turbTimeScale = nullptr;
  MBool* isAnomCell = nullptr;

 private:
  // no private data yet
};


#endif
