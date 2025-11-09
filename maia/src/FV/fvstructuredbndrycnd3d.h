// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef STRUCTUREDBNDRYCND3D
#define STRUCTUREDBNDRYCND3D

#include "fvstructuredbndrycnd.h"
#include "fvstructuredsolver.h"
#include "fvstructuredsolver3d.h"

class FvStructuredSolver3D;

/** \brief Class for the 3D stuctured boundary conditions
 *
 */
template <MBool isRans>
class StructuredBndryCnd3D : public StructuredBndryCnd<3> {
 public:
  friend class FvStructuredSolver3D;

  using Timers = maia::structured::Timers_;

  StructuredBndryCnd3D(FvStructuredSolver<3>* solver, StructuredGrid<3>* grid);
  ~StructuredBndryCnd3D();

  void createSpongeAtBndryCnd();
  void correctBndryCndIndices();

  // void         applyDirichletNeumannBC();
  void readAndDistributeSpongeCoordinates();
  void updateSpongeLayer();
  void computeWallDistances();
  void setRotationalBCProperties();


  inline MFloat dist(MFloat* a, MFloat* b);

  void bc1000(MInt); // wall no slip
  void bc1003(MInt); // isothermal no slip wall
  void bc1004(MInt); // moving adiabatic wall
  void bc1006(MInt); // moving isothermal wall
  void bc1007(MInt); // oscillating wall

  void bc2001(MInt);   // subsonic inflow
  void bc2002(MInt);   // supersonic inflow
  void bc2003(MInt);   // simple subsonic in/outflow
  void bc2004(MInt);   // subsonic outflow
  void bc2005(MInt);   // supersonic outflow
  void bc2009(MInt);   // supersonic outflow after shock
  void bc2012(MInt){}; // characteristic inflow
  void bc2013(MInt){}; // characteristic outflow
  void bc2014(MInt);   // subsonic rotational inflow
  void bc2020(MInt);   // poiseulle flow inflow
  void bc2097(MInt);   // plenum inflow
  void bc2099(MInt);   // subsonic inflow (u=(y/d)^(1/7))
  void bc2222(MInt);   // subsonic RANS outflow bc //junoh
  void bc2402(MInt);   // channel flow
  void bc2500(MInt);   // Rescaling: recycle station
  void bc2501(MInt);   // Rescaling: inlet station
  void bc2510(MInt);   // Rescaling: recycle station RANS
  void bc2511(MInt);   // Rescaling: inlet station RANS
  void bc2600(MInt);   // Prescribing profile
  void bc2601(MInt);   // Prescribing profile
  void bc2700(MInt);   // mode inflow
  void bc2730(MInt);   // fsc outflow
  void bc2888(MInt);   // fsc inflow
  void bc2999(MInt);   // blasius inflow
  void bc2900(MInt);   // Jet Inlet Freund

  void bc3000(MInt); // symmetry
  void bc3001(MInt); // streamline symmetry

  void bc6000(MInt);                  // communication
  virtual void bc6002(MInt) override; // Fluid-porous interface
  void bc7909(MInt);                  // synthetic turbulence generation

  void initBc1000(MInt);   // wall no slip
  void initBc1003(MInt);   // isothermal no slip wall
  void initBc1004(MInt);   // moving adiabatic wall
  void initBc1006(MInt);   // moving isothermal wall
  void initBc1007(MInt);   // oscillating wall
  void initBc2001(MInt){}; // subsonic inflow
  void initBc2002(MInt){}; // supersonic inflow
  void initBc2003(MInt){}; // simple subsonic in/outflow
  void initBc2004(MInt);   // subsonic outflow
  void initBc2005(MInt){}; // supersonic outflow
  void initBc2009(MInt);   // supersonic outflow after shock
  void initBc2012(MInt){}; // characteristic inflow
  void initBc2013(MInt){}; // characteristic outflow
  void initBc2014(MInt){}; // subsonic rotational bc
  void initBc2020(MInt){}; // poiseulle flow
  void initBc2222(MInt);
  void initBc2097(MInt);                  // plenum inflow
  void initBc2099(MInt){};                // subsonic inflow (u=(y/d)^(1/7))
  void initBc2402(MInt);                  // channel flow
  void initBc2500(MInt);                  // Rescaling: recycle station
  void initBc2501(MInt){};                // Rescaling: inlet station
  void initBc2600(MInt);                  // Prescribing profile
  void initBc2601(MInt);                  // Prescribing profile
  void initBc2700(MInt);                  // mode inflow
  void initBc2900(MInt){};                // Jet Inlet Freund
  void initBc3000(MInt){};                // symmetry
  void initBc3001(MInt){};                // streamline symmetry
  void initBc6000(MInt){};                // communication
  void initBc6002(MInt) override{};       // Fluid-porous interface
  void initBc7909(MInt);                  // synthetic turbulence generation

  // empty BC if BC==-1
  void bc9999(MInt){};
  void initBc9999(MInt){};

  inline void crossProduct(MFloat*, MFloat*, MFloat*);
  MInt cellIndex(MInt i, MInt j, MInt k);
  MInt cellIndexBC(MInt i, MInt j, MInt k); // For STG, only the first three rows are used

  virtual void computeFrictionPressureCoef(MBool computePower) override {
    if(computePower)
      computeFrictionPressureCoef_<true>();
    else
      computeFrictionPressureCoef_<false>();
  }
  template <MBool computePower>
  void computeFrictionPressureCoef_(const MBool auxDataWindows = false);
  void computeMomentCoef();

 protected:
  FvStructuredSolver3D* m_solver;
  MInt m_startCommPeriodic;
  MInt m_endCommPeriodic;
  MInt m_periodicS;
  MInt m_startCommChannel;
  MInt m_endCommChannel;
  MInt m_channelInflowRank;

  // periodic
  MInt m_noPeriodicConnections;


  // 7909 synthetic turbulence generation
  //  MFloat** m_stgEddies;
  MFloat generate_rand();
  MFloat generate_rand_weighted();
  MFloat* m_stgVbStart = nullptr;
  MFloat* m_stgVbEnd = nullptr;
  MFloat* m_stgMaxVel = nullptr;
  MFloat** m_stgGlobalLengthScales = nullptr;

  // 2700 mode bc
  MFloat m_isothermalWallTemperature;
  MFloat* m_modeAmp = nullptr; // for mode inflow
  MFloat m_modeSr;
  MInt* m_modeType = nullptr;
  MFloat* m_modePhi = nullptr;
  MFloat m_modeAngle;
  MInt* m_nmbrOfModes = nullptr;
  MInt m_modes;
  MFloat* m_modeOmega = nullptr;
  MFloat* m_modeEtaMin = nullptr;
  MFloat** m_modeK = nullptr;

  // 2601
  MFloat** m_2601effConst = nullptr;
  MFloat* m_2601streamwisePos = nullptr;
  MBool m_2601wave;
  MInt m_2601noSteps;
  MInt m_2601noCoeff;
  MInt m_2601noPos;

  // 2500 rescaling bc
  MFloat m_rescalingBLT;

  // routine to determine point ID from given cell
  inline MInt getPointIdFromCell(MInt i, MInt j, MInt k);
  inline MInt pointIndex(MInt i, MInt j, MInt k);
  // routine to determine point ID from given point
  inline MInt getPointIdfromPoint(MInt origin, MInt incI, MInt incJ, MInt incK);
  inline MFloat pressure(MInt cellId);
  inline MFloat pressure(MInt i, MInt j, MInt k);
  inline MFloat temperature(MInt cellId);
  inline MInt getReverseCellId(MInt i, MInt j, MInt k, MInt face);
  inline MInt getExtrNghbrId(MInt cellId, MInt face);
  inline std::pair<MInt, MInt> getMirrorCellIdPair(MInt i, MInt j, MInt k, MInt face);

  static constexpr MInt m_reverseCellIdDim[18] = {-1, 1, 1, 1, 1, 1, 1, -1, 1, 1, 1, 1, 1, 1, -1, 1, 1, 1};

  static constexpr MInt m_reverseCellIdGC[18] = {1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0};
  static constexpr const MInt nDim = 3;
};

#endif
