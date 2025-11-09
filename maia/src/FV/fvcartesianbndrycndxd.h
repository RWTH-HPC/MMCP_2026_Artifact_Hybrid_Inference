// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef FVBNDRYCNDXD_H
#define FVBNDRYCNDXD_H

#include <map>
#include <random>
#include <set>
#include <sys/stat.h>
#include <type_traits>
#include <vector>
#include "COMM/mpioverride.h"
#include "GEOM/geometryintersection.h"
#include "GRID/cartesiangridcellproperties.h"
#include "MEMORY/scratch.h"
#include "UTIL/maiamath.h"
#include "fvcartesiancellcollector.h"
#include "fvcartesiansurfacecollector.h"
#include "fvcartesiansyseqntraits.h"

template <MInt nDim, class SysEqn>
class FvCartesianSolverXD;
template <class T>
class Collector;
template <class T>
class List;
template <MInt nDim>
class FvSurface;
template <MInt nDim>
class FvSurfaceCollector;
template <MInt nDim, class SysEqn>
class FvBndryCell;
template <MInt nDim>
class CartesianGrid;
template <MInt nDim, SolverType SolverTypeR, SolverType SolverTypeL>
class MSTG;
template <MInt nDim, class SysEqn>
class FvZonalSTG;
template <MInt nDim>
class FvSysEqnNS;
template <MInt nDim, class RANSModel>
class FvSysEqnRANS;

template <MInt nDim>
class Bc1601Class {
  // members used by turbulent inflow boundary condition 1601
  MInt m_oldMode;
  MBool m_commValues;
  MBool m_useUnif;
  MInt m_noOfModes;
  MBool m_regenerateSeeding;
  MInt m_regenerationInterval;
  MInt m_regenerationCounter;
  std::mt19937 randNumGen;
  // mpi communicator
  MPI_Comm m_comm_bc1601;
  MInt m_rank_bc1601{};
  // Domain id of solver
  const MInt m_domainId;
  // to store random numbers
  MFloat* m_omega = nullptr;
  MFloat* m_p1 = nullptr;
  MFloat* m_p2 = nullptr;
  MFloat* m_p3 = nullptr;
  MFloat* m_q1 = nullptr;
  MFloat* m_q2 = nullptr;
  MFloat* m_q3 = nullptr;
  MFloat* m_dhat1 = nullptr;
  MFloat* m_dhat2 = nullptr;
  MFloat* m_dhat3 = nullptr;
  // to store Reynolds stress tensor:
  MFloat uuref;
  MFloat uvref;
  MFloat uwref;
  MFloat vvref;
  MFloat vwref;
  MFloat wwref;
  // to store cholesky decomposed Reynolds stress tensor:
  MFloat* aa = nullptr;

  // Shin's addition
  MBool smirnov;
  MFloat eigvec[3][3]{};
  MFloat c1, c2, c3;

  // end members used by bc1601

  // functions
  void generateAndCommRandomNumbers();
  MBool loadRandomNumbers();

 public:
  // precalculated values
  MFloat m_l_b;
  MFloat m_v_b;
  MFloat m_tau_b;
  MFloat m_invSigmaSponge;

  // functions
  Bc1601Class(MPI_Comm& communicator, const MInt m_solverId, const MInt domainId, MFloat& u_total,
              MFloat& invSigmaSponge);
  ~Bc1601Class();
  void calculateFlucts(const MFloat that, const MFloat xhat, const MFloat yhat, const MFloat zhat, MFloat* fluctChol);
  void checkRegeneration(const MFloat time);
};


template <MInt nDim, class SysEqn>
class FvBndryCndXD {
 public:
  // Type for cell properties
  using Cell = GridCell;
  using SolverCell = FvCell;

  MInt m_noDirs = 2 * nDim;
  MInt m_noEdges = nDim == 2 ? 4 : 12;
  MInt m_noCorners = nDim == 2 ? 4 : 8;

  friend class FvCartesianSolverXD<nDim, SysEqn>;
  friend class FvZonalSTG<nDim, SysEqn>;

 private:
  SysEqn* m_sysEqn;

 public:
  // collectors...
  maia::fv::collector::FvCellCollector<nDim>& m_cells;
  maia::fv::surface_collector::FvSurfaceCollector<nDim>& m_surfaces;
  List<MInt>* m_smallBndryCells = nullptr;
  List<MInt>* m_sortedBndryCells = nullptr;
  List<MInt>* m_sortedSpongeBndryCells = nullptr;
  std::vector<List<MInt>*> m_sortedCutOffCells{};
  // pointers
  FvCartesianSolverXD<nDim, SysEqn>* m_solver = nullptr;
  FvBndryCell<nDim, SysEqn>* m_bndryCell = nullptr;
  MInt* m_splitSurfaces = nullptr;
  MInt* m_bndryNghbrs = nullptr;
  MInt* m_splitParents = nullptr;
  MInt* m_splitChildren = nullptr;
  MInt m_solverId = -1;
  MBool m_changeAdiabBCToTemp;

  // claudia addition -> not nice, but necessary:
  MInt m_noLevelSetsUsedForMb;
  MBool m_complexBoundaryMB;
  MBool m_cellCoordinatesCorrected;

  typename SysEqn::ConservativeVariables* CV{};
  typename SysEqn::FluxVariables* FV{};
  typename SysEqn::PrimitiveVariables* PV{};
  typename SysEqn::AdditionalVariables* AV{};

  MInt m_minLevel;
  MInt m_maxLevel;

  FvBndryCndXD(FvCartesianSolverXD<nDim, SysEqn>* solver);

  virtual ~FvBndryCndXD();

  Collector<FvBndryCell<nDim, SysEqn>>* m_bndryCells = nullptr;

  SysEqn sysEqn() const { return *m_sysEqn; };

  virtual void writeStlFileOfCell(MInt, const MChar*){};
  void addBoundarySurfaces();
  void addBoundarySurfacesMGC();
  // virtual void correctBoundarySurfaceVariables() = 0;
  void correctBoundarySurfaceVariablesMGC();
  void correctBoundarySurfaceVariablesMGCSurface();
  //  void applyNonReflectingBC();
  void applyNonReflectingBCCutOff();
  void applyNonReflectingBCAfterTreatmentCutOff();
  void computeGhostCells();
  void computeGhostCellsMGC();
  //  void correctGhostCells();
  void computeMirrorCoordinates(MInt, MFloat*, MInt);
  void computeMirrorCoordinates(MInt, MFloat*);
  void computeReverseMap();

  // WMLES
  void initWMBndryCells();
  void initWMSurfaces();
  template <MBool MGC>
  void bc3399(MInt);
  std::list<std::pair<MInt, MInt>> m_wmSrfcToCellId;

  // virtual void copySlopesToSmallCells() = 0;
  void copyRHSIntoGhostCells();
  void correctCellCoordinates();
  void createSortedBndryCellList();
  void createBndryCndHandler();
  void createSpongeAtSpongeBndryCnds();
  void recorrectCellCoordinates();
  void rerecorrectCellCoordinates();
  void applyNeumannBoundaryCondition();
  void storeBoundaryVariables();
  void updateGhostCellVariables();
  void updateCutOffCellVariables();
  void initCutOffBndryCnds();
  void allocateCutOffMemory();
  void exchangeCutOffBoundaryCells();
  void updateGhostCellSlopesInviscid();
  void updateCutOffSlopesInviscid();
  void updateGhostCellSlopesViscous();
  void correctGhostCellSlopesViscous();
  void updateCutOffSlopesViscous();
  void initBndryCnds();
  void setBCTypes(MInt updateOnlyBndryCndId = -1);
  virtual void cmptGhostCells(){};

  void saveBc7901(){};

  void bc0(MInt){/*do nothing*/};
  void bc0Var(MInt){/*do nothing*/};

  void bcInit0001(MInt);
  template <MBool MGC>
  void bcInit0002(MInt);
  template <MBool MGC>
  void bcInit4000(MInt);
  void bcInit0004(MInt);
  void bcNeumann(MInt);
  void bcNeumannIso(MInt);
  void initBndryCommunications();
  void createBoundaryAtCutoff();
  MInt isCutOffInterface(MInt cellId);
  void markCutOff(MIntScratchSpace& cutOffCells);

  template <MBool MGC>
  void sbc2000(MInt);
  template <MInt dir>
  void sbc2001(MInt);
  template <MInt dir>
  void sbc2901(MInt);

  void bc17110(MInt);
  void bc17516(MInt);
  void bc19516(MInt);
  void bc1753(MInt);
  void bc1755(MInt);

  void bc1251(MInt);

  virtual void bc2770(MInt){};
  virtual void sbc2710co(MInt){};
  virtual void bc1001coflowY(MInt){};
  virtual void bc10910(MInt){};
  virtual void bc10990(MInt){};
  virtual void bc1901(MInt){};
  virtual void bc1801(MInt){};
  virtual void bc19520(MInt){};
  virtual void bc1003(MInt){};
  virtual void bc1004(MInt){};
  virtual void bc1005(MInt){};
  virtual void bc1401(MInt){};
  virtual void bc1791(MInt){};
  virtual void bc2001(MInt){};
  virtual void bc2002(MInt){};
  virtual void bc2003(MInt){};

  void bc3002(MInt);
  void bc30021(MInt);
  template <MBool MGC>
  void bc3003(MInt);
  void bc3011(MInt);
  void bc4000(MInt);
  void bc4001(MInt);

  // pointers to the boundary condition methods
  typedef void (FvBndryCndXD::*BndryCndHandler)(MInt);
  typedef void (FvBndryCndXD::*BndryCndHandlerVar)(MInt);
  BndryCndHandler* bndryCndHandlerInit = nullptr;
  BndryCndHandler* bndryCndHandlerCutOffInit = nullptr;
  BndryCndHandlerVar* bndryCndHandlerNeumann = nullptr;
  // variables
  BndryCndHandler* bndryCndHandlerVariables = nullptr;
  BndryCndHandler* bndryCndHandlerCutOffVariables = nullptr;
  BndryCndHandler* bndryCndHandlerSpongeVariables = nullptr;
  // BndryCndHandler* nonReflectingBoundaryCondition = nullptr;
  BndryCndHandler* nonReflectingBoundaryConditionAfterTreatmentCutOff = nullptr;
  BndryCndHandler* nonReflectingCutOffBoundaryCondition = nullptr;
  // gradient (inviscid flux)
  BndryCndHandler* bndryCndHandlerSlopesInviscid = nullptr;
  BndryCndHandler* bndryCndHandlerCutOffSlopesInviscid = nullptr;
  // gradient (viscous flux)
  BndryCndHandler* bndryViscousSlopes = nullptr;
  BndryCndHandler* bndryCutOffViscousSlopes = nullptr;

  // reconstruction
  MFloat** m_reconstructionConstants = nullptr;
  MInt** m_reconstructionNghbrs = nullptr;
  std::vector<MInt> m_smallCutCells;
  std::vector<MInt>* m_nearBoundaryWindowCells = nullptr;
  std::vector<MInt>* m_nearBoundaryHaloCells = nullptr;
  // Azimuthal priodicity near boundary exchange
  std::vector<std::vector<MInt>> m_azimuthalNearBoundaryWindowCells;
  std::vector<std::vector<MInt>> m_azimuthalNearBoundaryWindowMap;
  std::vector<std::vector<MInt>> m_azimuthalNearBoundaryHaloCells;

  MBool m_cellMerging = false;
  MBool m_secondOrderRec;
  MInt m_noFluxRedistributionLayers;

  // properties
  MFloat m_4000timeStepOffset;
  MFloat m_4000timeInterval;
  MFloat m_wFactor;
  MFloat m_pressureRatioChannel;

  MFloat m_sigmaNonRefl;
  MFloat m_sigmaNonReflInflow;
  MBool m_createBoundaryAtCutoff;
  MBool m_createSpongeBoundary;
  MBool m_outputIGPoints;
  MInt* m_bndryCndIds = nullptr;
  MInt* m_cutOffBndryCndIds = nullptr;
  MInt* m_cellsInsideSpongeLayer = nullptr;
  MInt m_noCellsInsideSpongeLayer;
  MInt* m_spongeBndryCndIds = nullptr; ///< holds the sponge boundary IDs
  MFloat* m_spongeFactor = nullptr;
  MInt* m_spongeDirections = nullptr;
  MFloat* m_sigmaSpongeBndryId = nullptr;
  MFloat* m_sigmaEndSpongeBndryId = nullptr;
  MFloat m_spongeLayerThickness;
  MInt m_spongeLayerLayout;
  MBool m_spongeTimeDep;
  MFloat* m_spongeStartIteration = nullptr;
  MFloat* m_spongeEndIteration = nullptr;
  MInt* m_spongeTimeDependent = nullptr;
  MInt* m_bndryCndCells = nullptr;
  MInt* m_spongeBndryCells = nullptr;
  MInt* m_boundarySurfaces = nullptr;
  MInt m_multipleGhostCells = 0;
  MInt m_ipVariableIterative;
  MInt m_surfaceGhostCell;
  MInt m_noBoundarySurfaces{};
  MInt m_noBndryCndIds;
  MInt m_noSpongeBndryCndIds; ///< number of sponge boundary condition IDs
  MInt m_noCutOffBndryCndIds;
  MBool m_cbcCutOff = false;
  MBool m_cbcSmallCellCorrection = false;
  MInt m_noMaxSpongeBndryCells;
  MFloat m_spongeBeta;
  MFloat* m_spongeCoord = nullptr;
  MInt m_maxNoBndryCells;
  MInt m_maxNoBndryCndIds;
  MInt m_noImagePointIterations;
  MFloat m_volumeLimitWall;
  MBool m_smallCellRHSCorrection = false;
  MFloat m_volumeLimitOther;
  MFloat m_meanCoord[3]{};
  // ---
  MFloat** m_shockBcVars = nullptr;
  MInt m_noShockBcCells{};
  MBool m_shockFromInnerSolution{};
  MInt* m_Bc2770TargetCells = nullptr;
  MFloat m_sigmaShock{};
  MFloat m_ys{};
  // ---
  MFloat* m_modeOmega = nullptr;
  MFloat* m_modeAmp = nullptr;
  MInt* m_modeType = nullptr;
  MFloat** m_modeK = nullptr;
  MBool m_clusterCutOffBcs{};
  MInt* m_nmbrOfModes = nullptr;
  MFloat* m_modePhi = nullptr;
  MFloat* m_modeEtaMin = nullptr;
  MInt m_modes{};
  MInt m_besselModes;
  MFloat* m_besselTrig = nullptr;
  MFloat m_Bc3011WallTemperature = NAN;

  // members used by combustion boundary conditions
  MFloat m_radiusFlameTube;
  MFloat m_radiusVelFlameTube;
  MFloat m_shearLayerThickness;
  MFloat m_jetHeight;
  MFloat m_primaryJetRadius;
  MFloat m_secondaryJetRadius;
  MFloat m_targetVelocityFactor;
  MFloat m_momentumThickness;
  MFloat m_shearLayerStrength;
  MFloat m_Ma;
  MFloat m_deltaP{};
  MFloat m_deltaPL{};
  MFloat m_inflowTemperatureRatio;
  MInt m_noSpecies;
  MInt m_noRansEquations;
  MBool m_combustion;
  MBool m_isEEGas;
  MFloat m_nu[10]{};

  MInt m_bc1601_bcId = -1;
  Bc1601Class<nDim>* m_bc1601 = nullptr;
  MBool m_bc1601MoveGenOutOfSponge = false;
  MBool m_jetInletTurbulence = false;

  // BC1251
  MFloat m_bc1251ForcingAmplitude;
  MFloat m_bc1251ForcingWavelength;
  MFloat m_bc1251ForcingFrequency;
  MInt m_bc1251ForcingDirection;

  // BC7901
  MInt m_7901faceNormalDir;
  MInt m_7901BcActive;
  MInt m_7901StartTimeStep;
  MInt m_7901periodicDir;
  MInt m_7901wallDir;
  MInt m_7901globalNoWallNormalLocations = -1;
  MInt* m_7901globalNoPeriodicLocations = nullptr;
  MFloat** m_7901LESAverage = nullptr;
  MFloat** m_7901LESAverageOld = nullptr;
  MInt* m_7901periodicIndex = nullptr;
  // std::vector<MInt> m_7901BcCells;
  std::vector<MFloat>* m_7901periodicLocations = nullptr;
  std::vector<MFloat> m_7901wallNormalLocations;
  std::vector<MFloat> m_7901globalWallNormalLocations;
  /* MInt m_7901StartTimeStep; */
  /* std::vector<MFloat>* m_7901LESAverage = nullptr; */

  // BC7902
  MInt m_rntRoot;
  MInt m_7902faceNormalDir;
  MInt m_7902BcActive;
  MInt m_7902StartTimeStep;
  MInt m_7902periodicDir;
  MInt m_7902wallDir;
  MInt m_7902globalNoWallNormalLocations = -1;
  MInt* m_7902globalNoPeriodicLocations = nullptr;
  MFloat** m_7902LESAverage = nullptr;
  MFloat** m_7902LESAverageOld = nullptr;
  MInt* m_7902periodicIndex = nullptr;
  // std::vector<MInt> m_7902BcCells;
  std::vector<MFloat>* m_7902periodicLocations = nullptr;
  std::vector<MFloat> m_7902wallNormalLocations;
  std::vector<MFloat> m_7902globalWallNormalLocations;

  // BC7903
  // BC7904
  // BC7909
  std::map<MInt, MSTG<nDim, MAIA_FINITE_VOLUME, MAIA_FINITE_VOLUME>*> m_stgBC;
  std::map<MInt, MSTG<nDim, MAIA_STRUCTURED, MAIA_FINITE_VOLUME>*> m_stgBCStrcd;
  std::map<MInt, MBool> m_stgLocal;
  std::vector<MInt> m_stgBcCells;
  MInt m_startSTGTimeStep;
  /* MBool m_cylinderTransformation = false; */
  // std::vector<std::pair<MInt,MInt>> m_stgBcCellIds;
  MPI_Comm m_commStg;

  // cbc
  std::vector<MFloat> m_cbcInflowArea;
  std::vector<std::vector<MFloat>> m_cbcReferencePoint;
  std::vector<MInt> m_cbcBndryCndIds;
  std::vector<std::vector<MFloat>> m_cbcRelax;
  std::vector<std::vector<MFloat>> m_dirTangent;
  std::vector<std::vector<MFloat>> m_dirNormal;
  std::vector<MFloat> m_cbcLref;
  std::vector<std::vector<MInt>> m_cbcDir;
  std::vector<MInt> m_cbcDomainMin;
  MBool m_cbcTurbulence = false;
  MBool m_cbcViscous = false;
  MFloat** m_oldFluctChol = nullptr;
  MFloat m_oldTime = F0;

 protected:
  // parallelization
  MString m_gridCutTest;

  MPI_Comm* m_comm_bc = nullptr;
  MInt m_comm_bc_init = 0;
  MInt* m_bc_comm_pointer = nullptr;
  MPI_Comm* m_comm_bcCo = nullptr;
  MInt m_comm_bcCo_init = 0;
  MInt* m_bcCo_comm_pointer = nullptr;

  std::vector<CutCandidate<nDim>> m_cutCandidates;
  GeometryIntersection<nDim>* m_geometryIntersection;

  /// Return the MPI communicator used by the corresponding solver
  MPI_Comm mpiComm() const { return m_solver->mpiComm(); }

  /// Return the domain id of this solver on the current MPI communicator
  MInt domainId() const { return m_solver->domainId(); }

  /// Return the total number of domains (total number of ranks in current MPI communicator)
  MInt noDomains() const { return m_solver->noDomains(); }

  // Data members used to replace local static variables
 private:
  MBool m_firstUseSetBCTypes = true;
  // Data members used to replace local static variables for fvbndrycndxd.h

 protected:
  // sbc1000co
  MBool m_static_sbc1000co_first = true;
  static constexpr MInt s_sbc1000co_fixedMaxNoBndryCndIds = 10;
  MInt m_static_sbc1000co_directions[s_sbc1000co_fixedMaxNoBndryCndIds]{};
  // initSmallCellCorrection
  MBool m_static_initSmallCellCorrection_firstRun = true;
  // m_static_computeImagePointRecConst_firstRun
  MBool m_static_computeImagePointRecConst_firstRun = true;

 public:
  // output routines to check the cut cell generation
  void plotAllCutPoints();
  /// Only 3D ///
  template <class X = void, std::enable_if_t<nDim == 3, X*> = nullptr>
  void writeStlFileOfCell(MInt, const char*);
  template <class X = void, std::enable_if_t<nDim == 3, X*> = nullptr>
  void plotSurface();
  template <class X = void, std::enable_if_t<nDim == 3, X*> = nullptr>
  void plotEdges(MInt&, MFloat**&);
  template <class X = void, std::enable_if_t<nDim == 3, X*> = nullptr>
  void plotIntersectionPoints(MInt*, MFloat***&);
  template <class X = void, std::enable_if_t<nDim == 3, X*> = nullptr>
  void writeStlOfNodes(MInt, MInt*&, const char*);
  template <class X = void, std::enable_if_t<nDim == 3, X*> = nullptr>
  void plotTriangle(std::ofstream&, MFloat*, MFloat*, MFloat*, MFloat*);

  // TODO labels:FV cleanup and replace with calls to maiamath.h
  // geometric help routines
  template <class X = void, std::enable_if_t<nDim == 3, X*> = nullptr>
  MFloat* vecSub(MFloat*, MFloat*, MFloat*);
  template <class X = void, std::enable_if_t<nDim == 3, X*> = nullptr>
  MFloat* vecScalarMul(MFloat, MFloat*, MFloat*);
  template <class X = void, std::enable_if_t<nDim == 3, X*> = nullptr>
  void computeTri(MFloat*, MFloat*, MFloat*, MFloat*, MFloat*);
  template <class X = void, std::enable_if_t<nDim == 3, X*> = nullptr>
  void computeTri(MFloat*, MFloat*, MFloat*, MFloat*, MFloat*, MFloat*);
  template <class X = void, std::enable_if_t<nDim == 3, X*> = nullptr>
  void computeTrapez(MFloat*, MFloat*, MFloat*, MFloat*, MFloat*, MFloat*, MFloat*);
  template <class X = void, std::enable_if_t<nDim == 3, X*> = nullptr>
  void computePoly3(MFloat*, MFloat*, MFloat*, MFloat*, MFloat*, MFloat*);
  template <class X = void, std::enable_if_t<nDim == 3, X*> = nullptr>
  void computePoly4(MFloat*, MFloat*, MFloat*, MFloat*, MFloat*, MFloat*, MFloat*);
  template <class X = void, std::enable_if_t<nDim == 3, X*> = nullptr>
  void computePoly5(MFloat*, MFloat*, MFloat*, MFloat*, MFloat*, MFloat*, MFloat*, MFloat*);
  template <class X = void, std::enable_if_t<nDim == 3, X*> = nullptr>
  void computePoly6(MFloat*, MFloat*, MFloat*, MFloat*, MFloat*, MFloat*, MFloat*, MFloat*, MFloat*);
  template <class X = void, std::enable_if_t<nDim == 3, X*> = nullptr>
  void computeTetra(MFloat*, MFloat*, MFloat*, MFloat*, MFloat*, MFloat*);
  template <class X = void, std::enable_if_t<nDim == 3, X*> = nullptr>
  void computePyra(MFloat*, MFloat*, MFloat*, MFloat*, MFloat*, MFloat*);
  template <class X = void, std::enable_if_t<nDim == 3, X*> = nullptr>
  void correctFace(MFloat*, MFloat*, MFloat*, MFloat*);
  template <class X = void, std::enable_if_t<nDim == 3, X*> = nullptr>
  void correctCell(MFloat*, MFloat*, MFloat*, MFloat*);
  template <class X = void, std::enable_if_t<nDim == 3, X*> = nullptr>
  void correctNormal(MFloat*);

  // other help routines for cut cell generation
  template <class X = void, std::enable_if_t<nDim == 3, X*> = nullptr>
  void getFeatureEdges(MInt&, MFloat**&, MInt, MInt*&, MFloat*&, MFloat*&);
  template <class X = void, std::enable_if_t<nDim == 3, X*> = nullptr>
  void getIntersectionPoints(MFloat**&, MFloat**&, MInt, MInt, MFloat**&, MInt&);
  template <class X = void, std::enable_if_t<nDim == 3, X*> = nullptr>
  void getSortedElements(const std::vector<MInt>&, MInt&, MInt*&, MInt&, MInt*&, MInt);
  template <class X = void, std::enable_if_t<nDim == 3, X*> = nullptr>
  void correctInflowBoundary(MInt);
  template <class X = void, std::enable_if_t<nDim == 3, X*> = nullptr>
  void correctInflowBoundary(MInt, MFloat*&, MFloat*&);
  template <class X = void, std::enable_if_t<nDim == 2, X*> = nullptr>
  void correctInflowBoundary(MInt, MBool = false) {}
  template <class X = void, std::enable_if_t<nDim == 2, X*> = nullptr>
  void correctInflowBoundary(MInt, MFloat*&, MFloat*&, MBool = false) {}
  template <class X = void, std::enable_if_t<nDim == 3, X*> = nullptr>
  MInt createSplitCell_MGC(MInt, MInt);
  MBool checkInside(MInt);
  MBool checkOutside(MInt);
  MBool checkOutside(const MFloat*, const MInt);

  // cut cell generation main routines
  virtual void generateBndryCells();
  void createBndryCells();
  void computeCutPoints();
  void exchangeComputedCutPoints();
  void checkCutPointsValidity();
  void checkCutPointsValidityParGeom();
  template <class X = void, std::enable_if_t<nDim == 2, X*> = nullptr>
  void createCutFace();
  template <class X = void, std::enable_if_t<nDim == 3, X*> = nullptr>
  void createCutFace();
  template <class X = void, std::enable_if_t<nDim == 3, X*> = nullptr>
  void createCutFaceMGC();
  template <class X = void, std::enable_if_t<nDim == 2, X*> = nullptr>
  void createCutFaceMGC() {
    TERMM(1, "MGC in 2D not implemented yet!");
  }
  void checkBoundaryCells();
  void computePlaneVectors();
  void deleteBndryCell(MInt);

  // other boundary cell related routines
  ATTRIBUTES2(ATTRIBUTE_HOT, ATTRIBUTE_FLATTEN) void copySlopesToSmallCells();
  // ATTRIBUTES2(ATTRIBUTE_HOT, ATTRIBUTE_ALWAYS_INLINE) inline void copySlopesToSmallCells_(const MUint);
  ATTRIBUTES2(ATTRIBUTE_HOT, ATTRIBUTE_FLATTEN) void copyVarsToSmallCells();
  // ATTRIBUTES2(ATTRIBUTE_HOT, ATTRIBUTE_ALWAYS_INLINE) inline void copyVarsToSmallCells_(const MUint);
  ATTRIBUTES2(ATTRIBUTE_HOT, ATTRIBUTE_FLATTEN) void correctBoundarySurfaceVariables();
  // ATTRIBUTES2(ATTRIBUTE_HOT, ATTRIBUTE_ALWAYS_INLINE) inline void correctBoundarySurfaceVariables_(const MUint);
  void correctCoarseBndryCells();
  void correctMasterSlaveSurfaces();
  void detectSmallBndryCells();
  void detectSmallBndryCellsMGC();
  void setNearBoundaryRecNghbrs(MInt updateOnlyBndryCndId = -1);
  void computeImagePointRecConst(MInt updateOnlyBndryCndId = -1);
  void initSmallCellCorrection(MInt updateOnlyBndryCndId = -1);
  void initSmallCellRHSCorrection(MInt updateOnlyBndryCndId = -1);
  void mergeCells();
  void mergeCellsMGC();
  void updateRHSSmallCells();
  void resetCutOff();
  void resetCutOffFirst();
  void resetBndryCommunication();
  void cutOffBcMissingNeighbor(const MInt cellId, const MString bcName);

  void setGapGhostCellVariables(MInt bcId);

  template <class X = void, std::enable_if_t<nDim == 2, X*> = nullptr>
  void computePolygon(MFloat* x, const MInt N, MFloat* centroid, MFloat* area);
  // boundary conditions
  void computeNeumannLSConstants(MInt);
  void computeReconstructionConstants_interpolation();
  void bcNeumann3600(MInt);
  void bcNeumannIsothermal(MInt);
  void bcNeumannIsothermalBurntProfile(MInt);
  void bcNeumannIsothermalBurntProfileH(MInt);
  void bcNeumannIsothermalUnburnt(MInt);
  void bcNeumannIsothermalUnburntProfile(MInt);
  void bcNeumannIsothermalUnburntProfileH(MInt);
  void bcNeumannMb(MInt);
  MFloat updateImagePointVariables(MInt);

  // slopes
  ATTRIBUTES2(ATTRIBUTE_HOT, ATTRIBUTE_FLATTEN) void sbc1000(const MInt);
  void sbc1000co(const MInt);
  void sbc00co(const MInt);
  // ATTRIBUTES2(ATTRIBUTE_HOT, ATTRIBUTE_ALWAYS_INLINE) inline void sbc1000_(const MInt, const MUint);
  void sbc1002(MInt);
  void sbc2801x(MInt);
  void sbc2801y(MInt);

  void bc1001(MInt);
  void bc100100(MInt);
  void bc1000(MInt);
  void bc1009(MInt);
  void bc01(MInt);
  void bc1002(MInt);
  void bc1091(MInt);
  void bc1101(MInt);
  void bc1102(MInt);
  void bc1601(MInt);
  void bc1602(MInt);
  void bc1603(MInt);
  void bc1604(MInt);
  void bc1606(MInt);
  void bcInit1601(MInt);
  void bc10970(MInt);
  void bc10980(MInt);
  void bc11110(MInt);
  void bc16010(MInt);
  void bc16011(MInt);
  void bc16012(MInt);
  void bc16013(MInt);
  void bc16014(MInt);
  void bc16015(MInt);

  // STG
  // Only 3D TODO labels:FV template std::enable_if_t can be removed once FV/fvstg.h is trully generalized
  //               for the templated argument 'nDim'. Otherwise, leads to compilation errors.
  //               std::enable_if_t<nDim == (3) or (2)... is used here to avoid such errors.
  template <class _ = void, std::enable_if_t<hasPV_N<SysEqn>::value, _*> = nullptr,
            std::enable_if_t<nDim == 3, _*> = nullptr>
  void bcInit7901(MInt);
  template <class _ = void, std::enable_if_t<hasPV_N<SysEqn>::value, _*> = nullptr,
            std::enable_if_t<nDim == 3, _*> = nullptr>
  void bcInit7902(MInt);
  template <class _ = void, std::enable_if_t<hasPV_N<SysEqn>::value, _*> = nullptr,
            std::enable_if_t<nDim == 3, _*> = nullptr>
  void bcInit7905(MInt);
  template <class _ = void, std::enable_if_t<!hasPV_N<SysEqn>::value, _*> = nullptr,
            std::enable_if_t<nDim == 3, _*> = nullptr>
  void bcInit7909(MInt);
  template <class _ = void, std::enable_if_t<!hasPV_N<SysEqn>::value, _*> = nullptr,
            std::enable_if_t<nDim == 3, _*> = nullptr>
  void bcInit7809(MInt);

  template <class _ = void, std::enable_if_t<hasPV_N<SysEqn>::value, _*> = nullptr,
            std::enable_if_t<nDim == 3, _*> = nullptr>
  void bc7901(MInt);
  template <class _ = void, std::enable_if_t<hasPV_N<SysEqn>::value, _*> = nullptr,
            std::enable_if_t<nDim == 3, _*> = nullptr>
  void bc7902(MInt);
  template <class _ = void, std::enable_if_t<hasPV_N<SysEqn>::value, _*> = nullptr,
            std::enable_if_t<nDim == 3, _*> = nullptr>
  void bc7905(MInt);
  template <class _ = void, std::enable_if_t<!hasPV_N<SysEqn>::value, _*> = nullptr,
            std::enable_if_t<nDim == 3, _*> = nullptr>
  void bc7903(MInt);

  template <class _ = void, std::enable_if_t<!hasPV_N<SysEqn>::value, _*> = nullptr,
            std::enable_if_t<nDim == 3, _*> = nullptr>
  void bc7909(MInt);
  template <class _ = void, std::enable_if_t<!hasPV_N<SysEqn>::value, _*> = nullptr,
            std::enable_if_t<nDim == 3, _*> = nullptr>
  void bc7809(MInt);


  // Only 2D
  template <class _ = void, std::enable_if_t<hasPV_N<SysEqn>::value, _*> = nullptr,
            std::enable_if_t<nDim == 2, _*> = nullptr>
  void bcInit7901(MInt) {
    TERMM(-1, "");
  }
  template <class _ = void, std::enable_if_t<hasPV_N<SysEqn>::value, _*> = nullptr,
            std::enable_if_t<nDim == 2, _*> = nullptr>
  void bcInit7902(MInt) {
    TERMM(-1, "");
  }
  template <class _ = void, std::enable_if_t<hasPV_N<SysEqn>::value, _*> = nullptr,
            std::enable_if_t<nDim == 2, _*> = nullptr>
  void bcInit7905(MInt) {
    TERMM(-1, "");
  }
  template <class _ = void, std::enable_if_t<!hasPV_N<SysEqn>::value, _*> = nullptr,
            std::enable_if_t<nDim == 2, _*> = nullptr>
  void bcInit7909(MInt) {
    TERMM(-1, "");
  }
  template <class _ = void, std::enable_if_t<!hasPV_N<SysEqn>::value, _*> = nullptr,
            std::enable_if_t<nDim == 2, _*> = nullptr>
  void bcInit7809(MInt) {
    TERMM(-1, "");
  }

  template <class _ = void, std::enable_if_t<hasPV_N<SysEqn>::value, _*> = nullptr,
            std::enable_if_t<nDim == 2, _*> = nullptr>
  void bc7901(MInt) {
    TERMM(-1, "");
  }
  template <class _ = void, std::enable_if_t<hasPV_N<SysEqn>::value, _*> = nullptr,
            std::enable_if_t<nDim == 2, _*> = nullptr>
  void bc7902(MInt) {
    TERMM(-1, "");
  }
  template <class _ = void, std::enable_if_t<hasPV_N<SysEqn>::value, _*> = nullptr,
            std::enable_if_t<nDim == 2, _*> = nullptr>
  void bc7905(MInt) {
    TERMM(-1, "");
  }
  template <class _ = void, std::enable_if_t<!hasPV_N<SysEqn>::value, _*> = nullptr,
            std::enable_if_t<nDim == 2, _*> = nullptr>
  void bc7903(MInt) {
    TERMM(-1, "");
  }
  template <class _ = void, std::enable_if_t<!hasPV_N<SysEqn>::value, _*> = nullptr,
            std::enable_if_t<nDim == 2, _*> = nullptr>
  void bc7909(MInt) {
    TERMM(-1, "");
  }
  template <class _ = void, std::enable_if_t<!hasPV_N<SysEqn>::value, _*> = nullptr,
            std::enable_if_t<nDim == 2, _*> = nullptr>
  void bc7809(MInt) {
    TERMM(-1, "");
  }

  void bcInit30022(MInt);
  void bc30022(MInt);

  void bc1792(MInt);
  void bc1156(MInt);
  void bc1952(MInt);

  void bcInit1050(MInt);

  void initModes(MInt);
  void addModes(MInt);
  void initBesselModes(MInt);
  void addBesselModes(MInt);
  void calcBesselFractions(const MFloat, const MFloat, const MFloat, MFloat&, MFloat&);
  void precomputeBesselTrigonometry(MInt);
  void bc2700(MInt);
  void bcInit1251(MInt);
  void bcInit2770(MInt);
  void bcInit2700(MInt);
  void bc2710(MInt);
  void bc2720(MInt);
  void sbc2720co(MInt);
  void bc3006(MInt);
  void bc3007(MInt);
  void bc2907(MInt);
  void bc29050(MInt);
  void bc3466(MInt);
  void bc3600(MInt);

  MFloat computeCutoffBoundaryGeometry(const MInt, const MInt, MFloat*);

  // cbc
  void cbcTurbulenceInjection(MInt, MFloat*, MInt);
  void cbcRHS(MInt, MInt, MFloat*, MFloat*, MFloat*);
  void cbcDampingOutflow(MInt, MInt, MFloat, MFloat*);
  void cbcDampingInflow(MInt, MInt, MFloat, MFloat*, MString);
  template <unsigned char vTerms>
  void cbcViscousTerms(MInt, MInt, MFloat*, MFloat*, MFloat*, MFloat*, MInt*, MFloat*);
  template <unsigned char tTerms>
  void cbcTransversalTerms(MInt, MInt, MFloat*, MFloat*, MFloat*, MFloat*, MFloat*);
  template <MInt side>
  void cbcOutgoingAmplitudeVariation(MInt, MInt, MFloat*, MFloat*, MFloat*, MFloat*, MFloat*);
  void cbcGradients(MInt, MInt, MFloat*, MFloat*, MFloat*, MFloat*);
  void cbcGradientsViscous(MInt, MInt, MFloat*, MFloat*, MFloat*, MFloat*, MInt*);
  void cbcTauQ(MInt, MFloat*, MFloat*, MInt*);
  void cbcMachCo(MInt, MFloat*);
  void cbcMeanPressureCo(MInt, MFloat*);
  void bcInitCbc(MInt);


  void cbc1099(MInt);
  void cbc109910(MInt);
  void cbc109911(MInt);
  void cbc109921(MInt);
  void cbc1099a(MInt);
  void cbc1099a_after(MInt);
  void cbc1091(MInt);
  void cbc1091a(MInt);
  void cbc3091a(MInt);
  void cbc1091b(MInt);
  void cbc1091c(MInt);
  void cbc1091d(MInt);
  void cbc1091b_after(MInt);
  void cbc1091c_after(MInt);
  void cbc1091d_after(MInt);
  void cbc1091e(MInt);
  void cbc1091e_after(MInt);
  void cbc1099_1091_local(MInt);
  void cbc1099_1091_local_comb(MInt);
  void cbc1291(MInt);
  void cbc1291a(MInt);
  void cbc1291b(MInt);
  void cbc1291tm(MInt);
  void cbc1291tma(MInt);
  void cbc1291tmb(MInt);
  void cbc1291tmc(MInt);
  void cbc1299(MInt);
  void cbc1299tm(MInt);
  void cbc1299a(MInt);
  void cbc2091a(MInt);
  void cbc2091b(MInt);
  void cbc2091b_after(MInt);
  void cbc2091d(MInt);
  void cbc2091d_after(MInt);
  void cbc2099_1091_local_comb(MInt);
  void cbc1099_1091_engine(MInt);

  template <class X = void, std::enable_if_t<nDim == 3, X*> = nullptr>
  void cbc1099_1091_engineOld(MInt);

  virtual void bc3037MGC(MInt){};
  void bc1091MGC(MInt);
  void bc1099MGC(MInt);
  void cbc1099_1091d(MInt);
  void cbc1099_1091d_after(MInt);

  // Data members used to replace local static variables
 private:
  /// Only 3D ///
  // plotEdges
  MInt m_static_plotEdges_iter = 0;
  // plotIntersectionPoints
  MInt m_static_plotIntersectionPoints_iter = 0;
  // writeStlOfNodes
  MInt m_static_writeStlOfNodes_iter = 0;
  // correctInflowBoundary
  MInt m_static_correctInflowBoundary_iter = 0;
  // bc1091MGC
  MInt m_static_bc1091MGC_minTimeSteps = 0;
  MInt m_static_bc1091MGC_nghbrDir = -1;
  MInt m_static_bc1091MGC_edgeCellCounter = 0;
  MInt m_static_bc1091MGC_first = true;
  MInt m_static_bc1091MGC_first2 = true;
  // bc1099MGC
  MFloat m_static_bc1099MGC_timeOfMaxPdiff = 0;
  bool m_static_bc1099MGC_first = true;
  /// Only 3D ends ///
  // bc10970
  MBool m_firstUseBc10970 = true;
  MInt m_dirNBc10970 = -1;
  MInt m_pModeBc10970 = 0;
  // bc10980
  MBool m_firstUseBc10980 = true;
  MInt m_dirNBc10980 = -1;
  // bc11110
  MFloat** m_targetValuesBC11110 = nullptr;
  MInt* m_dirNBc11110 = nullptr;
  MBool m_firstUseBc11110 = true;
  // cbc1099
  MBool m_static_cbc1099_first = true;
  MInt m_static_cbc1099_dirN = -1;
  MInt m_static_cbc1099_dimN = -1;
  MInt m_static_cbc1099_dimT1 = -1;
  MInt m_static_cbc1099_dimT2 = -1;
  MFloat m_static_cbc1099_outFlowArea;
  MFloat m_static_cbc1099_referencePoint[3];
  MInt m_static_cbc1099_domainMin;
  // cbc1099_1091_local
  MBool m_static_cbc1099_1091_local_first = true;
  MInt m_static_cbc1099_1091_local_dirN = -1;
  MInt m_static_cbc1099_1091_local_dimN = -1;
  MInt m_static_cbc1099_1091_local_dimT1 = -1;
  MInt m_static_cbc1099_1091_local_dimT2 = -1;
  MFloat m_static_cbc1099_1091_local_inflowArea;
  MFloat m_static_cbc1099_1091_local_referencePoint[3];
  MFloat m_static_cbc1099_1091_local_targetPressure;
  MInt m_static_cbc1099_1091_local_domainMin = 0;
  MFloat m_static_cbc1099_1091_local_R = 0; // 3D
  MFloat m_static_cbc1099_1091_local_H;     // 2D
  // cbc1099_1091_local_comb
  MBool m_static_cbc1099_1091_local_comb_first = true;
  MInt m_static_cbc1099_1091_local_comb_dirN = -1;
  MInt m_static_cbc1099_1091_local_comb_dimN = -1;
  MInt m_static_cbc1099_1091_local_comb_dimT1 = -1;
  MInt m_static_cbc1099_1091_local_comb_dimT2 = -1;
  MFloat m_static_cbc1099_1091_local_comb_outFlowArea;
  MFloat m_static_cbc1099_1091_local_comb_referencePoint[3];
  // cbc1091
  MBool m_static_cbc1091_first = true;
  MInt m_static_cbc1091_dirN = -1;
  MInt m_static_cbc1091_dimN = -1;
  MInt m_static_cbc1091_dimT1 = -1;
  MInt m_static_cbc1091_dimT2 = -1;
  MFloat m_static_cbc1091_inflowArea;
  MFloat m_static_cbc1091_referencePoint[3];
  MInt m_static_cbc1091_domainMin;
  // cbc1091a
  MBool m_static_cbc1091a_solverProfile = false;
  // cbc3091a
  MBool m_static_cbc3091a_first = true;
  MInt m_static_cbc3091a_dirN = -1;
  MInt m_static_cbc3091a_dimN = -1;
  MInt m_static_cbc3091a_dimT1 = -1;
  MInt m_static_cbc3091a_dimT2 = -1;
  MBool m_static_cbc3091a_solverProfile = false;
  MFloat m_static_cbc3091a_inflowArea;
  MFloat m_static_cbc3091a_R;
  MFloat m_static_cbc3091a_referencePoint[3];
  // cbc1091b
  MBool m_static_cbc1091b_first = true;
  MInt m_static_cbc1091b_dirN = -1;
  MInt m_static_cbc1091b_dimN = -1;
  MInt m_static_cbc1091b_dimT1 = -1;
  MBool m_static_cbc1091b_solverProfile = false;
  MFloat m_static_cbc1091b_inflowArea;
  MFloat m_static_cbc1091b_referencePoint[3];
  // cbc1091c
  MBool m_static_cbc1091c_first = true;
  MInt m_static_cbc1091c_dirN = -1;
  MInt m_static_cbc1091c_dimN = -1;
  MInt m_static_cbc1091c_dimT1 = -1;
  MFloat m_static_cbc1091c_inflowArea;
  MFloat m_static_cbc1091c_referencePoint[3];
  // cbc1091c_after
  MBool m_static_cbc1091c_after_first = true;
  MInt m_static_cbc1091c_after_dirN = -1;
  MInt m_static_cbc1091c_after_dimN = -1;
  MInt m_static_cbc1091c_after_dimT1 = -1;
  MInt m_static_cbc1091c_after_dimT2 = -1;
  // cbc1091d
  MBool m_static_cbc1091d_first = true;
  MInt m_static_cbc1091d_dirN = -1;
  MInt m_static_cbc1091d_dimN = -1;
  MInt m_static_cbc1091d_dimT1 = -1;
  MFloat m_static_cbc1091d_inflowArea;
  MFloat m_static_cbc1091d_referencePoint[3];
  MBool m_static_cbc1091d_solverProfile = false; // 2D
  MInt m_static_cbc1091d_minDom = 0;
  // cbc1091d_after
  MBool m_static_cbc1091d_after_first = true;
  MInt m_static_cbc1091d_after_dirN = -1;
  MInt m_static_cbc1091d_after_dimN = -1;
  MInt m_static_cbc1091d_after_dimT1 = -1;
  MInt m_static_cbc1091d_after_dimT2 = -1;
  // cbc1091e
  MBool m_static_cbc1091e_first = true;
  MInt m_static_cbc1091e_dirN = -1;
  MInt m_static_cbc1091e_dimN = -1;
  MInt m_static_cbc1091e_dimT1 = -1;
  MBool m_static_cbc1091e_solverProfile = false;
  MFloat m_static_cbc1091e_inflowArea;
  // cbc2091a
  MBool m_static_cbc2091a_first = true;
  MInt m_static_cbc2091a_dirN = -1;
  MInt m_static_cbc2091a_dimN = -1;
  MInt m_static_cbc2091a_dimT1 = -1;
  MBool m_static_cbc2091a_solverProfile = false;
  MFloat m_static_cbc2091a_inflowArea;
  MFloat m_static_cbc2091a_H;
  MFloat m_static_cbc2091a_referencePoint[nDim];
  // cbc2091b
  MBool m_static_cbc2091b_first = true;
  MInt m_static_cbc2091b_dirN = -1;
  MInt m_static_cbc2091b_dimN = -1;
  MInt m_static_cbc2091b_dimT1 = -1;
  MBool m_static_cbc2091b_solverProfile = false;
  MFloat m_static_cbc2091b_inflowArea;
  MFloat m_static_cbc2091b_H;
  MFloat m_static_cbc2091b_referencePoint[nDim];
  // cbc2091d
  MBool m_static_cbc2091d_first = true;
  MInt m_static_cbc2091d_dirN = -1;
  MInt m_static_cbc2091d_dimN = -1;
  MInt m_static_cbc2091d_dimT1 = -1;
  MFloat m_static_cbc2091d_inflowArea;
  MFloat m_static_cbc2091d_H;
  MFloat m_static_cbc2091d_referencePoint[nDim];
  MBool m_static_cbc2091d_solverProfile = false;
  // cbc2091d_after
  MBool m_static_cbc2091d_after_first = true;
  MInt m_static_cbc2091d_after_dirN = -1;
  MInt m_static_cbc2091d_after_dimN = -1;
  MInt m_static_cbc2091d_after_dimT1 = -1;
  // cbc2099_1091_local_comb
  MBool m_static_cbc2099_1091_local_comb_first = true;
  MInt m_static_cbc2099_1091_local_comb_dirN = -1;
  MInt m_static_cbc2099_1091_local_comb_dimN = -1;
  MInt m_static_cbc2099_1091_local_comb_dimT1 = -1;
  MFloat m_static_cbc2099_1091_local_comb_outFlowArea;
  MFloat m_static_cbc2099_1091_local_comb_referencePoint[nDim];
  // cbc1099_1091d
  MBool m_static_cbc1099_1091d_first = true;
  MInt m_static_cbc1099_1091d_dirN = -1;
  MInt m_static_cbc1099_1091d_dimN = -1;
  MInt m_static_cbc1099_1091d_dimT1 = -1;
  MInt m_static_cbc1099_1091d_dimT2 = -1;
  MFloat m_static_cbc1099_1091d_inflowArea;
  MFloat m_static_cbc1099_1091d_referencePoint[nDim];
  MFloat m_static_cbc1099_1091d_targetPressure;
  // cbc1099_1091d_after
  MBool m_static_cbc1099_1091d_after_first = true;
  MInt m_static_cbc1099_1091d_after_dirN = -1;
  MInt m_static_cbc1099_1091d_after_dimN = -1;
  MFloat m_static_cbc1099_1091d_after_interpolationFactor;
  // cbc1099_1091_engine
  MBool m_static_cbc1099_1091_engine_first = true;
  MInt m_static_cbc1099_1091_engine_dirN;
  MInt m_static_cbc1099_1091_engine_dimN;
  MInt m_static_cbc1099_1091_engine_dimT1;
  MFloat m_static_cbc1099_1091_engine_inflowArea;
  MInt m_static_cbc1099_1091_engine_domainMin = 0;
  MInt m_static_cbc1099_1091_engine_dimT2;
  MFloat m_static_cbc1099_1091_engine_normal[3];
  MFloat m_static_cbc1099_1091_engine_tangent[3];

  // JANNIK: 109910, 109911
  std::pair<MFloat, MFloat>* m_unTargetData = nullptr;
  std::pair<MFloat, MFloat>* m_vnTargetData = nullptr;
  MInt m_unTargetDataCount;
  MInt m_vnTargetDataCount;
  MFloat** m_horTargetData = nullptr;
  MInt m_horTargetDataCount;
};

#endif
