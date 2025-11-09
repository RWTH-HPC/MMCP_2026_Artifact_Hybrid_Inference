// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef FVSTRUCTUREDSOLVER_2DRANS
#define FVSTRUCTUREDSOLVER_2DRANS

#include "INCLUDE/maiatypes.h"
#include "MEMORY/scratch.h"
#include "fvransmodelconstants.h"
#include "fvstructuredsolver.h"
#include "fvstructuredsolver2d.h"

class FvStructuredSolver2D;

class FvStructuredSolver2DRans {
 public:
  FvStructuredSolver2DRans(FvStructuredSolver2D* solver);
  ~FvStructuredSolver2DRans();

  // void setBndryCndObject(class StructuredBndryCnd2D* object, MInt noSpecies);
  void initFluxMethod();

  template <MInt noVars, MInt noRANS, MBool rans2eq_production_mode = false>
  void Muscl_AusmDV();
  template <MInt noVars, MInt noRANS, MBool rans2eq_production_mode = false>
  void Muscl_AusmDV_Limited();
  template <MInt noVars, MInt noRANS, MBool rans2eq_production_mode = false>
  void Muscl_Ausm_Limited();
  void (FvStructuredSolver2DRans::*reconstructSurfaceData)();
  void Muscl(MInt timerId = -1);

  void viscousFluxRANS();
  void viscousFlux_SA(); // SPALART-ALLMARAS
  void viscousFlux_KEPSILON();
  void viscousFlux_KEPSILON2();
  void (FvStructuredSolver2DRans::*viscFluxMethod)();
  void (FvStructuredSolver2DRans::*compTurbVisc)();
  void computeTurbViscosity();
  void computeTurbViscosity_SA();
  void computeTurbViscosity_KEPSILON();
  void setAndAllocate_KEPSILON();
  void diffusiveFluxCorrection();

 protected:
  class StructuredBndryCnd2DRans* m_structuredBndryCndRans;


 private:
  RansMethod m_ransMethod;
  FvStructuredSolver2D* m_solver;

  MPI_Comm m_StructuredComm;
  MInt m_solverId;
  MInt* m_nCells;
  MInt* m_nPoints;
  MInt m_noCells;
  StructuredCell* m_cells;
  std::unique_ptr<MConservativeVariables<2>>& CV;
  std::unique_ptr<MPrimitiveVariables<2>>& PV;
  std::unique_ptr<StructuredFQVariables>& FQ;
  MInt m_noGhostLayers;
  MFloat m_eps;
  MFloat m_chi;
  MBool m_P_keps = false;
  static constexpr const MInt nDim = 2;
  const MFloat m_sutherlandConstant;
  const MFloat m_sutherlandPlusOne;

  const MInt xsd = 0;
  const MInt ysd = 1;

  MBool m_dsIsComputed;

  inline MInt cellIndex(MInt i, MInt j);
  inline MInt getCellIdfromCell(MInt origin, MInt incI, MInt incJ);
  inline MFloat getPSI(MInt, MInt);
};

#endif
