// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef STRUCTUREDBNDRYCND
#define STRUCTUREDBNDRYCND

#include "GRID/structuredgrid.h"
#include "fvstructuredcell.h"
#include "fvstructuredfqvariables.h"
#include "fvstructuredsolver.h"
#include "fvstructuredsolverwindowinfo.h"
#include "fvstructuredwindowmapping.h"
#include "globals.h"


// Forward declarations
template <MInt nDim>
class FvStructuredSolver;

/** \brief Base class of the structured boundary conditions
 *
 */
template <MInt nDim>
class StructuredBndryCnd {
  template <MInt nDim_>
  friend class FvStructuredSolver;
  friend class FvStructuredSolver3D;
  template <MInt nDim_>
  friend class FvStructuredSolverWindowInfo;

 public:
  // functions
  StructuredBndryCnd(FvStructuredSolver<nDim>* solver, StructuredGrid<nDim>* grid);
  virtual ~StructuredBndryCnd();
  void applyNonReflectingBC();
  void assignBndryCnds();
  void applyDirichletNeumannBC();
  virtual void correctBndryCndIndices();
  // virtual void periodicExchange(){};

  // MPI Communicator
  MPI_Comm m_StructuredComm;

  // variables
  MInt* m_nCells = nullptr;
  MInt* m_nPoints = nullptr;
  StructuredCell* m_cells = nullptr;
  FvStructuredSolver<nDim>* m_solver = nullptr;
  StructuredGrid<nDim>* m_grid = nullptr;
  std::unique_ptr<MConservativeVariables<nDim>>& CV;
  std::unique_ptr<MPrimitiveVariables<nDim>>& PV;
  std::unique_ptr<StructuredFQVariables>& FQ;
  MInt m_solverId;
  // stores the number of boundary conditions;
  MInt m_noBndryCndIds;

  MInt m_noGhostLayers;
  // for Boundary conditions
  std::vector<std::unique_ptr<StructuredWindowMap<nDim>>>& m_physicalBCMap;
  std::vector<std::unique_ptr<StructuredWindowMap<nDim>>>& m_auxDataMap;
  std::vector<std::unique_ptr<StructuredWindowMap<nDim>>>& m_globalStructuredBndryMaps;

  MInt m_noSpongeDomainInfos;
  MInt* m_spongeBcWindowInfo = nullptr;
  // MInt       m_useSponge;
  MInt m_spongeLayerType;
  MFloat* m_spongeLayerThickness = nullptr;
  MFloat* m_sigmaSponge = nullptr;
  MFloat* m_betaSponge = nullptr;
  MFloat m_targetDensityFactor;
  MFloat m_sigma;

  void saveAuxData();

  virtual void bc0(MInt){};    // nothing to do bc
  virtual void bc1000(MInt){}; // wall no slip
  virtual void bc1001(MInt){}; // euler wall
  virtual void bc1003(MInt){}; // isothermal no slip wall
  virtual void bc1004(MInt){}; // moving adiabatic wall
  virtual void bc1006(MInt){}; // moving isothermal wall
  virtual void bc1007(MInt){}; // oscillating wall
  virtual void bc2001(MInt){}; // subsonic inflow
  virtual void bc2002(MInt){}; // supersonic inflow
  virtual void bc2003(MInt){}; // simple subsonic in/outflow
  virtual void bc2004(MInt){}; // subsonic outflow
  virtual void bc2024(MInt){}; // subsonic outflow
  virtual void bc2005(MInt){}; // supersonic outflow
  virtual void bc2006(MInt){}; // subsonic in/outflow with zero velocity
  virtual void bc2007(MInt){}; // supersonic outflow
  virtual void bc2009(MInt){}; // supersonic outflow after shock
  virtual void bc2020(MInt){}; // poiseulle flow inflow
  virtual void bc2021(MInt){}; // shear flow inflow
  virtual void bc2097(MInt){}; // plenum inflow
  virtual void bc2099(MInt){}; // subsonic inflow (u=(y/d)^(1/7))

  virtual void bc2221(MInt){}; // zonal with STG
  virtual void bc2222(MInt){}; // zonal without STG

  virtual void bc2199(MInt){};     // subsonic inflow compressible bernoulli (tfs2099)x
  virtual void bc2402(MInt){};     // channel flow
  virtual void bc3000(MInt){};     // symmetry
  virtual void bc3001(MInt){};     // streamline symmetry
  virtual void bc6000(MInt){};     // communication
  virtual void bc6002(MInt){};     // Fluid-porous interface
  virtual void bc2012(MInt){};     // characteristic inflow
  virtual void bc2014(MInt){};     // subsonic rotational inflow
  virtual void bc2015(MInt){};     // non reflecting outflow poinsot lele
  virtual void bc7909(MInt){};     // synthetic turbulence generation
  virtual void bc2013(MInt){};     // characteristic outflow
  virtual void bc2500(MInt){};     // Rescaling: recycle station
  virtual void bc2511(MInt){};     // Rescaling: inlet station RANS
  virtual void bc2510(MInt){};     // Rescaling: recycle station RANS
  virtual void bc2501(MInt){};     // Rescaling: inlet station
  virtual void bc2600(MInt){};     // Prescribing profile
  virtual void bc2601(MInt){};     // Prescribing profile
  virtual void bc2700(MInt){};     // mode inflow
  virtual void bc2730(MInt){};     // fsc outflow
  virtual void bc2888(MInt){};     // fsc inflow
  virtual void bc2999(MInt){};     // blasius inflow
  virtual void bc2900(MInt){};     // Jet Freund inlet
  virtual void initBc0(MInt){};    // wall no slip
  virtual void initBc1000(MInt){}; // wall no slip
  virtual void initBc1001(MInt){}; // euler wall
  virtual void initBc1003(MInt){}; // isothermal wall no slip
  virtual void initBc1004(MInt){}; // moving adiabatic wall
  virtual void initBc1006(MInt){}; // moving isothermal wall
  virtual void initBc1007(MInt){}; // oscillating wall
  virtual void initBc2001(MInt){}; // subsonic inflow
  virtual void initBc2003(MInt){}; // simple in/outflow
  virtual void initBc2004(MInt){}; // subsonic outflow
  virtual void initBc2024(MInt){}; // subsonic outflow
  virtual void initBc2002(MInt){}; // supersonic inflow
  virtual void initBc2005(MInt){}; // supersonic outflow
  virtual void initBc2006(MInt){};
  virtual void initBc2007(MInt){}; // supersonic outflow
  virtual void initBc2009(MInt){}; // supersonic outflow after shock
  virtual void initBc2020(MInt){}; // poiseulle flow inflow
  virtual void initBc2021(MInt){}; // shear flow inflow
  virtual void initBc2097(MInt){}; // plenum inflow
  virtual void initBc2099(MInt){}; // subsonic inflow (u=(y/d)^(1/7))

  virtual void initBc2221(MInt){}; // zonal with STG
  virtual void initBc2222(MInt){}; // zonal without STG

  virtual void initBc2199(MInt){}; // subsonic inflow compressible bernoulli
  virtual void initBc2402(MInt){}; // channel flow
  virtual void initBc3000(MInt){}; // symmetry
  virtual void initBc3001(MInt){}; // streamline symmetry
  virtual void initBc6000(MInt){}; // communication
  virtual void initBc6002(MInt){}; // Fluid-porous interface
  virtual void initBc2012(MInt){}; // characteristic inflow
  virtual void initBc7909(MInt){}; // synthetic turbulence generation
  virtual void initBc2013(MInt){}; // characteristic outflow
  virtual void initBc2014(MInt){}; // subsonic rotational inflow
  virtual void initBc2015(MInt){}; // non reflecting outflow poinsot lele
  virtual void initBc2500(MInt){}; // Rescaling: recycle station
  virtual void initBc2501(MInt){}; // Rescaling: inlet station
  virtual void initBc2510(MInt){}; // Rescaling: inlet station
  virtual void initBc2600(MInt){}; // Prescribing profile
  virtual void initBc2601(MInt){}; // Prescribing profile
  virtual void initBc2700(MInt){}; // mode inflow
  virtual void initBc2730(MInt){}; // fsc outflow
  virtual void initBc2888(MInt){}; // fsc inflow
  virtual void initBc2999(MInt){}; // blasius inflow
  virtual void initBc2900(MInt){}; // Jet Freund inlet
  virtual void bc9999(MInt){};
  virtual void initBc9999(MInt){};

  /* virtual void exchangePointsPeriodic(){}; */
  /* virtual void exchangePointsPeriodicS(){}; */


  virtual void computeWallDistances(){};
  virtual void computeLocalWallDistances(){};
  virtual void computeLocalExtendedDistancesAndSetComm() { TERMM(1, "For your nDim it is not implemented yet!"); };
  virtual void updateSpongeLayer(){};
  virtual void computeFrictionCoef() { mTerm(-1, "Method does not exist"); };
  virtual void computeFrictionPressureCoef(MBool /*computePower*/) = 0;
  virtual void distributeWallAndFPProperties(){};


  // function pointer to the boundary method used;
  typedef void (StructuredBndryCnd::*BndryCndHandler)(MInt);

  // Dirichlet Conditions
  BndryCndHandler* nonReflectingBoundaryCondition = nullptr;
  BndryCndHandler* bndryCndHandler = nullptr;
  BndryCndHandler* initBndryCndHandler = nullptr;
  // compute cf
  BndryCndHandler* skinFrictionHandler = nullptr;

 protected:
  // For communication of nearest wall cell properties
  std::vector<MInt> m_wallSendCells;
  std::vector<MInt> m_wallSendcounts;
  std::map<MInt, std::tuple<MInt, MInt, MFloat>> m_wallCellId2recvCell;
  // For communication of nearest fluid-porous interface properties
  std::vector<MInt> m_FPSendCells;
  std::vector<MInt> m_FPSendcounts;
  std::map<MInt, std::tuple<MInt, MInt, MFloat>> m_FPCellId2recvCell;
  //  std::map<MInt, MInt> m_FPCellId2recvCell;
  std::vector<MInt> m_FPSendCells_porous;
  std::vector<MInt> m_FPSendcounts_porous;
  std::map<MInt, std::tuple<MInt, MInt, MFloat>> m_FPCellId2recvCell_porous;

  // for the channel flow we need the surfaces:
  // sorted for in and out if one partition contains both parts of the channel


  MInt m_noCells;
  MFloat m_channelSurfaceIn;
  MFloat m_channelSurfaceOut;
  MInt* m_channelSurfaceIndexMap = nullptr;

  // plenum boundary condition
  MInt* m_plenumSurfaceIndexMap = nullptr;
  MFloat m_plenumSurface;

  MFloat m_sutherlandPlusOne;
  MFloat m_sutherlandConstant;

  // infintiy Values!!!!
  MFloat m_UInfinity;
  MFloat m_VInfinity;
  MFloat m_WInfinity;
  MFloat m_PInfinity;
  MFloat m_TInfinity;
  MFloat m_DthInfinity;
  MFloat m_muInfinity;
  MFloat m_DInfinity;
  MFloat m_VVInfinity[3];
  MFloat m_rhoUInfinity;
  MFloat m_rhoVInfinity;
  MFloat m_rhoWInfinity;
  MFloat m_rhoEInfinity;
  MFloat m_rhoInfinity;
  MFloat m_rhoVVInfinity[3];
};

#endif
