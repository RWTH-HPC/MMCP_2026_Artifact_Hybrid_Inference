// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef FVSTRUCTUREDSOLVER_3DRANS
#define FVSTRUCTUREDSOLVER_3DRANS

#include "INCLUDE/maiatypes.h"
#include "MEMORY/scratch.h"
#include "fvransmodelconstants.h"
#include "fvstructuredsolver.h"
#include "fvstructuredsolver3d.h"

class FvStructuredSolver3D;

class FvStructuredSolver3DRans {
 public:
  FvStructuredSolver3DRans(FvStructuredSolver3D* solver);
  ~FvStructuredSolver3DRans();

  void initFluxMethod();

  template <MInt noVars>
  void Muscl_AusmDV();
  template <MInt noVars>
  void Muscl_AusmDV_Limited();
  void (FvStructuredSolver3DRans::*reconstructSurfaceData)();
  void Muscl(MInt timerId = -1);

  void viscousFluxRANS();
  void viscousFlux_SA(); // SPALART-ALLMARAS
  void viscousFlux_FS(); // FARES-SCHROEDER
  void (FvStructuredSolver3DRans::*viscFluxMethod)();
  void (FvStructuredSolver3DRans::*compTurbVisc)();
  void computeTurbViscosity();
  void computeTurbViscosity_SA();
  void computeTurbViscosity_FS();

 protected:
  class StructuredBndryCnd3DRans* m_structuredBndryCndRans;

 private:
  RansMethod m_ransMethod;
  FvStructuredSolver3D* m_solver;

  MPI_Comm m_StructuredComm;
  MInt m_solverId;
  MInt* m_nCells;
  MInt* m_nPoints;
  MInt m_noCells;
  StructuredCell* m_cells;
  std::unique_ptr<MConservativeVariables<3>>& CV;
  std::unique_ptr<MPrimitiveVariables<3>>& PV;
  std::unique_ptr<StructuredFQVariables>& FQ;
  MInt m_noGhostLayers;
  MFloat m_eps;
  MFloat m_chi;
  static constexpr const MInt nDim = 3;
  const MFloat m_sutherlandConstant;
  const MFloat m_sutherlandPlusOne;
  MBool m_dsIsComputed;

  const MInt xsd = 0;
  const MInt ysd = 1;
  const MInt zsd = 2;

  inline MInt cellIndex(MInt i, MInt j, MInt k);
  inline MInt getCellIdfromCell(MInt origin, MInt incI, MInt incJ, MInt incK);
  inline MFloat getPSI(MInt, MInt);
};

#endif
