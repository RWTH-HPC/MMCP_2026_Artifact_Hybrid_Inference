// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef FVSTRUCTUREDSOLVER2D_H
#define FVSTRUCTUREDSOLVER2D_H

#include "fvstructuredbndrycnd2d.h"
#include "fvstructuredsolver.h"
#include "fvstructuredsolver2drans.h"

extern MInt globalTimeStep;

class FvStructuredSolver2DRans;

/** \brief 2D structured solver class
 *
 */
class FvStructuredSolver2D : public FvStructuredSolver<2> {
  template <MBool isRans>
  friend class StructuredBndryCnd2D;
  friend class FvStructuredSolver2DRans;

 public:
  FvStructuredSolver2D(MInt, StructuredGrid<2>*, MBool*, const MPI_Comm comm);
  ~FvStructuredSolver2D();


  virtual void initialCondition();
  virtual void initSolutionStep(MInt mode);

  void computeCellCentreCoordinates();
  void addGhostPointCoordinateValues();
  void extrapolateGhostPointCoordinatesBC();


  void initMovingGrid();
  void moveGrid(const MBool isRestart, const MBool zeroPos) override;

  void initFluxMethod();
  void assignBndryCells();
  void initBndryCnds();
  void allocateSingularities();
  void applyBoundaryCondition();
  void loadRestartBC2600();
  void computePrimitiveVariables() override;
  template <MFloat (FvStructuredSolver::*)(MInt) const = &FvStructuredSolver::dummy>
  void computePrimitiveVariables_();
  void convertSA2KEPS();

  void gather(const MBool, std::vector<std::unique_ptr<StructuredComm<2>>>&) override;
  void scatter(const MBool, std::vector<std::unique_ptr<StructuredComm<2>>>&) override;

  void computeCellLength();

  virtual void computeTimeStep();
  MBool maxResidual();
  MBool rungeKuttaStep();
  void updateSpongeLayer();
  void viscousFlux();
  void viscousFluxRANS();
  template <MBool twoEqRans = false>
  void viscousFluxLES();
  template <MBool twoEqRans = false>
  void viscousFluxLESCompact();
  void (FvStructuredSolver2D::*viscFluxMethod)();

  // save the skin-friction etc ...
  virtual void computeFrictionPressureCoef(MBool computePower) override {
    m_structuredBndryCnd->computeFrictionPressureCoef(computePower);
  }

  inline MInt getPointIdFromCell(MInt i, MInt j);
  inline MInt cellIndex(MInt i, MInt j);
  inline MInt pointIndex(MInt i, MInt j);
  inline MInt getPointIdFromPoint(MInt origin, MInt incI, MInt incJ);
  inline MInt getCellIdFromCell(MInt origin, MInt incI, MInt incJ);
  inline MFloat crossProduct(MFloat vec1[2], MFloat vec2[2]);
  MFloat pressure(MInt cellId);
  inline MFloat getPSI(MInt, MInt);

  void Ausm();
  inline void AusmLES(MFloat* QLeft, MFloat* QRight, const MInt dim, const MInt cellId);
  inline void AusmLES_PTHRC(MFloat* QLeft, MFloat* QRight, MInt dim, MInt cellId);
  inline void AusmDV(MFloat* QLeft, MFloat* QRight, const MInt dim, const MInt cellId);
  inline void AusmLES_MovingGrid(MFloatScratchSpace& QLeft, MFloatScratchSpace& QRight, MInt dim, MInt cellId);

  typedef void (FvStructuredSolver2D::*fluxmethod)(MFloat*, MFloat*, MInt, MInt);
  template <fluxmethod ausm, MInt noVars>
  void Muscl_();
  template <fluxmethod ausm, MInt noVars>
  void MusclStretched_();

  // function pointer for the Muscl_scheme
  void MusclRANS();
  void (FvStructuredSolver2D::*reconstructSurfaceData)();

  // different Muscl schemes
  void Muscl(MInt timerId = -1) override;
  void MusclAlbada();
  void MusclNoLimiter();

  void computeReconstructionConstantsSVD();

  template <MBool twoEqRans = false>
  void viscousFluxCorrection();
  template <MBool twoEqRans = false>
  void viscousFluxCompactCorrection();

  // Porous
  virtual void computePorousRHS(MBool /*isRans*/) override;
  void computePorousRHSCorrection();
  void exchange6002();

  FvStructuredSolver2DRans* m_ransSolver;

 protected:
  StructuredBndryCnd<2>* m_structuredBndryCnd;
  static constexpr const MInt nDim = 2;
};

#endif
