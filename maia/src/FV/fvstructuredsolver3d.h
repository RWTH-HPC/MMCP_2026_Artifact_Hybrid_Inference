// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef FVSTRUCTUREDSOLVER3D_H
#define FVSTRUCTUREDSOLVER3D_H

#include "GRID/structuredgrid.h"
#include "fvstructuredbndrycnd3d.h"
#include "fvstructuredsolver.h"
#include "fvstructuredsolver3drans.h"
// fftw is needed for the channel flow
#include <complex>
#include <fftw3.h>

extern MInt globalTimeStep;
class FvStructuredSolver3DRans;

enum FLUXMETHOD {
  AUSM_LES,
  AUSM_LES_PTHRC,
  AUSMDV,
};

class ParallelIoHdf5;

/** \brief 3D structured solver class
 *
 */
class FvStructuredSolver3D : public FvStructuredSolver<3> {
  template <MBool isRans>
  friend class StructuredBndryCnd3D;
  friend class FvStructuredSolver3DRans;

 public:
  FvStructuredSolver3D(MInt, StructuredGrid<3>*, MBool*, const MPI_Comm comm);
  ~FvStructuredSolver3D();

  void addGhostPointCoordinateValues();
  void extrapolateGhostPointCoordinatesBC();

  void initSolutionStep(MInt mode);
  void initialCondition() override;
  void assignBndryCells();
  void nonReflectingBC();
  void applyBoundaryCondition() override;
  void calcSurfaceMetrics();
  void initBndryCnds();
  void allocateSingularities();
  void applySandpaperTrip();
  void applySandpaperTripAirfoil();
  void initSandpaperTrip();
  void tripForceCoefficients(MFloat*, MFloat*, MFloat*, MInt, MInt);
  void tripFourierCoefficients(MFloat*, MInt, MFloat, MFloat);
  void computeDomainWidth() override;

  void computeZonalConnections();
  void zonalExchange() override;
  void zonalAllreduce();
  void zonalRealInterpolation();
  void zonalGather();
  void zonalSend();
  void zonalReceive();
  void zonalCheck();
  void zonalBufferAverage();
  void zonalScatter();
  std::vector<std::unique_ptr<StructuredZonalBC>> m_zonalBC;
  void spanwiseAvgZonal(std::vector<MFloat*>&) override;

  void distributeFluxToCells();
  void computeTimeStep() override;
  void computeVorticity() override;
  void computeLambda2Criterion() override;
  MFloat computeTotalKineticEngergy();
  MFloat computeTotalPressure();
  void initializeNeighbourCellIds();
  void initFluxMethod();
  void computePrimitiveVariables() override;
  template <MFloat (FvStructuredSolver::*)(MInt) const = &FvStructuredSolver::dummy>
  void computePrimitiveVariables_();
  MFloat pressure(MInt);
  inline MFloat getPSI(MInt, MInt);

  void Ausm();
  inline void AusmLES(MFloat* QLeft, MFloat* QRight, const MInt dim, const MInt cellId);
  inline void AusmLES_PTHRC(MFloat* QLeft, MFloat* QRight, MInt dim, MInt cellId);
  inline void AusmDV(MFloat* QLeft, MFloat* QRight, MInt dim, MInt cellId);

  typedef void (FvStructuredSolver3D::*fluxmethod)(MFloat*, MFloat*, MInt, MInt);
  template <MInt noVars>
  void Muscl_AusmLES();
  template <MInt noVars>
  void Muscl_AusmLES_PTHRC();
  template <MInt noVars>
  void Muscl_AusmDV();
  template <fluxmethod ausm, MInt noVars>
  void MusclStretched_();

  void Muscl(MInt timerId = -1) override;
  void MusclRANS();
  void (FvStructuredSolver3D::*reconstructSurfaceData)();
  // different Muscl schemes:
  void MusclVenkatakrishnan3D();
  void MusclAlbada();
  void MusclMinModLimiter();
  void computeCellLength();

  void computeReconstructionConstantsSVD();

  // limiter functions for MusclVenkatakrishnan3D()
  void VENKATAKRISHNAN_MOD_FCT(MFloat effNghbrDelta, MFloat srfcDelta, MFloat dxEpsSqr, MInt cellPos, MInt var,
                               MFloatScratchSpace& minPhi);
  void VENKATAKRISHNAN_FCT(MFloat effNghbrDelta, MFloat srfcDelta, MFloat dxEpsSqr, MInt cellPos, MInt var,
                           MFloatScratchSpace& minPhi);
  void BARTH_JESPERSON_FCT(MFloat effNghbrDelta, MFloat srfcDelta, MFloat dxEpsSqr, MInt cellPos, MInt var,
                           MFloatScratchSpace& minPhi);
  void (FvStructuredSolver3D::*Venkatakrishnan_function)(MFloat, MFloat, MFloat, MInt, MInt, MFloatScratchSpace&);

  MBool maxResidual();
  MBool rungeKuttaStep();
  void updateSpongeLayer();
  void addDisturbance();
  void viscousFlux();
  void viscousFluxRANS();
  template <MBool twoEqRans = false>
  void viscousFluxLES();
  void (FvStructuredSolver3D::*viscFluxMethod)();
  void moveGrid(const MBool isRestart, const MBool zeroPos) override;
  void applyBodyForce(const MBool isRestart, const MBool zeroPos) override;
  void initMovingGrid() override;
  void initBodyForce() override;
  // for the treatment of boundary condition
  void applyViscousBoundaryCondition();
  void applyInviscidBoundaryCondition();
  // save the skin-friction etc ...
  virtual void computeFrictionPressureCoef(MBool computePower) override {
    m_structuredBndryCnd->computeFrictionPressureCoef(computePower);
  }

  // fftw is needed for the channel flow:
  void initFFTW(fftw_complex* uPhysField, fftw_complex* vPhysField, fftw_complex* wPhysField, MInt lx, MInt ly, MInt lz,
                MInt noPeakModes);
  void getFourierCoefficients(MFloat* k, MFloat k0, std::complex<MFloat>* fourierCoefficient);
  MFloat randnormal(MFloat mu, MFloat sigma);
  void loadRestartBC2600();
  void loadRestartBC2601();
  void loadRestartSTG(MBool);

  void gather(const MBool, std::vector<std::unique_ptr<StructuredComm<3>>>&) override;
  void scatter(const MBool, std::vector<std::unique_ptr<StructuredComm<3>>>&) override;

  // for traveling wave averaging
  void waveGather();
  void waveScatter();
  void waveSend(std::vector<MPI_Request>&);
  void waveReceive(std::vector<MPI_Request>&);
  void waveExchange() override;
  void spanwiseWaveReorder() override;

  // for exchange of postprocessing field
  void gcFillGhostCells(std::vector<MFloat*>&);
  void gcExtrapolate(std::vector<MFloat*>&);

  void viscousFluxCorrection();

  // For time averaging of flow variables
  void computeCumulativeAverage(MBool forceReset) override;

  // Line and plane interpolation
  void initInterpolatedPoints();
  void saveInterpolatedPoints() override;
  void saveNodalBoxes() override;
  void savePointsToAsciiFile(MBool) override;
  void initPointsToAsciiFile() override;

  // Porous
  virtual void computePorousRHS(MBool /*isRans*/) override;
  void exchange6002();

 protected:
  std::unique_ptr<StructuredBndryCnd<3>> m_structuredBndryCnd;

  // index variables
  static constexpr MInt xsd = 0;
  static constexpr MInt ysd = 1;
  static constexpr MInt zsd = 2;

  MFloat m_kineticEOld;

  inline void crossProduct(MFloat*, const MFloat*, const MFloat*);
  inline MInt cellIndex(const MInt, const MInt, const MInt);
  inline MInt pointIndex(const MInt, const MInt, const MInt);
  inline MFloat dist(MFloat* a, MFloat* b);
  // routine to determine point ID from given cell ( point (0,0,0) for unit cube )
  inline MInt getPointIdFromCell(const MInt, const MInt, const MInt);
  // routine to determine point ID from given point
  inline MInt getPointIdfromPoint(const MInt, const MInt, const MInt, const MInt);
  inline MInt getCellIdfromCell(const MInt, const MInt, const MInt, const MInt);
  inline MInt surfId(MInt point, MInt isd, MInt dim);

  std::unique_ptr<FvStructuredSolver3DRans> m_ransSolver;


  // interpolation
  std::unique_ptr<StructuredInterpolation<3>> m_structuredInterpolation;
  std::unique_ptr<StructuredInterpolation<3>> m_pointInterpolation;
  void manualInterpolationCorrection();
  void interpolateFromDonor();

  // postprocessing helper functions
  void loadSampleFile(MString) override;
  void getSampleVariables(MInt, MFloat*) override;
  MFloat getSampleVorticity(MInt, MInt) override;
  void loadAverageRestartFile(const MChar*, MFloat**, MFloat**, MFloat**, MFloat**) override;
  void loadAveragedVariables(const MChar*) override;
  void shiftAverageCellValuesRestart();
  void shiftAverageCellValues();

  MFloat dvardxyz(MInt, MInt, MFloat*) override;
  MFloat dvardx(MInt, MFloat*) override;

  MFloat getCellLengthY(MInt, MInt, MInt) override;
  MFloat getCellCoordinate(MInt, MInt) override;

  static constexpr const MInt nDim = 3;
};


#endif
