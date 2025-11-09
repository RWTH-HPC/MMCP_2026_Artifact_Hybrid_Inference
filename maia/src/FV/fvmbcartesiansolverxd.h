// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef FVMBCARTESIANSOLVERXD_H
#define FVMBCARTESIANSOLVERXD_H

#include <algorithm>
#include <map>
#include <random>
#include <set>
#include <unordered_map>
#include <utility>
#include <vector>
#include "GEOM/geometryintersection.h"
#include "GRID/cartesiangridgenpar.h"
#include "GRID/partition.h"
#include "UTIL/kdtree.h"
#include "UTIL/pointbox.h"
#include "fvcartesiancellcollector.h"
#include "fvcartesiansolverxd.h"

//#include "/pds/opt/fftw-3.3.4-mpi/include/fftw3-mpi.h"

#ifndef MAIA_MS_COMPILER
#include <netdb.h>
#include <sys/socket.h>
#endif

#define FINITE_VOLUME_METHOD
#include "GRID/cartesiannetcdf.h"
#undef FINITE_VOLUME_METHOD

#ifndef MB_EXTRA_STUFF_
#define MB_EXTRA_STUFF_
//#define DEBUG(a,b) {if((globalTimeStep%100==0)||(globalTimeStep%100==99)){if(b==16)m_log << "in:
//"<< a<<endl; else if(b==32)m_log << "out: "<< a<<endl;}} #define DEBUG(a,b) {if(b==16)m_log <<
//"in: "<< a<<endl; else if(b==32)m_log << "out: "<< a<<endl;} #define DEBUG(a,b) {if(b==16){cerr
//<< "in: "<<g_domainId<<" "<<a<<endl;cerr.flush();}else if(b==32){cerr << "out: "<<g_domainId<<"
//"<< a<<endl;cerr.flush();}} #define DEBUG_LOG(a) {m_log << a << endl;}
#define DEBUG_LOG(a)

const MInt m_revDir[6] = {1, 0, 3, 2, 5, 4};

#endif

template <MInt nDim, class SysEqn>
class FvCartesianSolverXD;

template <MInt nDim, class SysEqn>
class FvMbCartesianSolverXD : public FvCartesianSolverXD<nDim, SysEqn> {
 public:
  using Base = FvCartesianSolverXD<nDim, SysEqn>;
  using SolverCell = FvCell;
  using Cell = typename maia::grid::tree::Tree<nDim>::Cell;
  using Base::a_associatedBodyIds;
  using Base::a_bndryId;
  using Base::a_cellVolume;
  using Base::a_coordinate;
  using Base::a_hasNeighbor;
  using Base::a_hasProperty;
  using Base::a_isBndryGhostCell;
  using Base::a_isGapCell;
  using Base::a_isHalo;
  using Base::a_isPeriodic;
  using Base::a_isWindow;
  using Base::a_level;
  using Base::a_levelSetValuesMb;
  using Base::a_noCells;
  using Base::a_noReconstructionNeighbors;
  using Base::a_noSurfaces;
  using Base::a_oldVariable;
  using Base::a_properties;
  using Base::a_pvariable;
  using Base::a_rightHandSide;
  using Base::a_slope;
  using Base::a_spongeFactor;
  using Base::a_surfaceArea;
  using Base::a_surfaceBndryCndId;
  using Base::a_surfaceCoordinate;
  using Base::a_surfaceDeltaX;
  using Base::a_surfaceFactor;
  using Base::a_surfaceFlux;
  using Base::a_surfaceNghbrCellId;
  using Base::a_surfaceOrientation;
  using Base::a_surfaceUpwindCoefficient;
  using Base::a_surfaceVariable;
  using Base::a_variable;
  using Base::allocateCommunicationMemory;
  using Base::assertValidGridCellId;
  using Base::c_cellLengthAtCell;
  using Base::c_cellLengthAtLevel;
  using Base::c_childId;
  using Base::c_globalId;
  using Base::c_isLeafCell;
  using Base::c_isToDelete;
  using Base::c_neighborId;
  using Base::c_noCells;
  using Base::c_noChildren;
  using Base::c_parentId;
  using Base::calcSlopesAfterStep;
  using Base::computeDomainAndSpongeDimensions;
  using Base::computeDomainLength;
  using Base::computePrimitiveVariables;
  using Base::copyGridProperties;
  using Base::CV;
  using Base::domainId;
  using Base::entropy;
  using Base::exchangeData;
  using Base::exchangeDataFV;
  using Base::exchangeGapInfo;
  using Base::findNghbrIds;
  using Base::FV;
  using Base::getAssociatedInternalCell;
  using Base::getIdentifier;
  using Base::grid;
  using Base::haloCellId;
  using Base::initNearBoundaryExchange;
  using Base::initSpongeLayer;
  using Base::isActive;
  using Base::isMultilevel;
  using Base::isMultilevelLowestSecondary;
  using Base::isMultilevelPrimary;
  using Base::m_activeCellIds;
  using Base::m_adaptation;
  using Base::m_adaptationInterval;
  using Base::m_associatedBodyIds;
  using Base::m_bndryCells;
  using Base::m_bndryCellSurfacesOffset;
  using Base::m_bndryGhostCellsOffset;
  using Base::m_bndrySurfacesOffset;
  using Base::m_bodyAcceleration;
  using Base::m_bodyAngularAcceleration;
  using Base::m_bodyAngularVelocity;
  using Base::m_bodyCenter;
  using Base::m_bodyHeatFlux;
  using Base::m_bodyTemperature;
  using Base::m_bodyTemperatureDt1;
  using Base::m_bodyVelocity;
  using Base::m_cells;
  using Base::m_cellsInsideSpongeLayer;
  using Base::m_cfl;
  using Base::m_chi;
  using Base::m_closeGaps;
  using Base::m_constructGField;
  using Base::m_createSpongeBoundary;
  using Base::m_deleteNeighbour;
  using Base::m_DInfinity;
  using Base::m_dragOutputInterval;
  using Base::m_DthInfinity;
  using Base::m_dualTimeStepping;
  using Base::m_eps;
  using Base::m_euler;
  using Base::m_gamma;
  using Base::m_gapCellId;
  using Base::m_gapCells;
  using Base::m_gapInitMethod;
  using Base::m_geometry;
  using Base::m_globalUpwindCoefficient;
  using Base::m_hInfinity;
  using Base::m_identNghbrIds;
  using Base::m_initialCondition;
  using Base::m_internalBodyId;
  using Base::m_kInfinityFactor;
  using Base::m_levelSet;
  using Base::m_levelSetAdaptationScheme;
  using Base::m_levelSetMb;
  using Base::m_levelSetRans;
  using Base::m_levelSetValuesMb;
  using Base::m_loadBalancingReinitStage;
  using Base::m_LsRotate;
  using Base::m_Ma;
  using Base::m_maxLevelBeforeAdaptation;
  using Base::m_mpi_receiveRequest;
  using Base::m_mpi_sendRequest;
  using Base::m_noActiveCells;
  using Base::m_noActiveHaloCellOffset;
  using Base::m_noDirs;
  using Base::m_noEmbeddedBodies;
  using Base::m_noGapRegions;
  using Base::m_noLevelSetsUsedForMb;
  using Base::m_noPeriodicGhostBodies;
  using Base::m_noRansEquations;
  using Base::m_noRKSteps;
  using Base::m_noSets;
  using Base::m_noSpecies;
  using Base::m_nuTildeInfinity;
  using Base::m_omegaInfinity;
  using Base::m_omegaInfinityFactor;
  using Base::m_orderOfReconstruction;
  using Base::m_outerBandWidth;
  using Base::m_outputOffset;
  using Base::m_physicalTime;
  using Base::m_physicalTimeStep;
  using Base::m_PInfinity;
  using Base::m_pipeRadius;
  using Base::m_Pr;
  using Base::m_Re;
  using Base::m_reconstructionDataSize;
  using Base::m_recordBodyData;
  using Base::m_referenceLength;
  using Base::m_restart;
  using Base::m_restartInterval;
  using Base::m_rhoInfinity;
  using Base::m_rhoUInfinity;
  using Base::m_rhoVVInfinity;
  using Base::m_rhoWInfinity;
  using Base::m_rhs0;
  using Base::m_SInfinity;
  using Base::m_slopeMemory;
  using Base::m_solutionOffset;
  using Base::m_solutionOutput;
  using Base::m_solverId;
  using Base::m_splitCells;
  using Base::m_splitChilds;
  using Base::m_splitChildToSplitCell;
  using Base::m_splitMpiCommRecv;
  using Base::m_static_getDistanceSplitSphere_first;
  using Base::m_static_getDistanceSplitSphere_h;
  using Base::m_static_updateBodyProperties_c453_firstRun;
  using Base::m_static_updateBodyProperties_firstTime;
  using Base::m_storeNghbrIds;
  using Base::m_surfaces;
  using Base::m_sutherlandPlusOne;
  using Base::m_sysEqn;
  using Base::m_time;
  using Base::m_timeRef;
  using Base::m_totalnosplitchilds;
  using Base::m_UInfinity;
  using Base::m_useCentralDifferencingSlopes;
  using Base::m_VInfinity;
  using Base::m_volumeAcceleration;
  using Base::m_volumeForcingDir;
  using Base::m_volumeThreshold;
  using Base::m_VVInfinity;
  using Base::m_wasAdapted;
  using Base::m_wasBalancedZonal;
  using Base::m_weightActiveCell;
  using Base::m_weightBaseCell;
  using Base::m_weightBndryCell;
  using Base::m_weightLeafCell;
  using Base::m_weightMulitSolverFactor;
  using Base::m_weightNearBndryCell;
  using Base::m_WInfinity;
  using Base::maxNoGridCells;
  using Base::maxUniformRefinementLevel;
  using Base::minLevel;
  using Base::mpiComm;
  using Base::noDomains;
  using Base::noHaloLayers;
  using Base::noInternalCells;
  using Base::noNeighborDomains;
  using Base::PV;
  using Base::resetCutOffCells;
  using Base::setActiveFlag;
  using Base::setConservativeVariables;
  using Base::setPrimitiveVariables;
  using Base::sysEqn;
  using Base::timeStep;

  using Base::a_copyPropertiesSolver;
  using Base::a_cutCellLevel;
  using Base::a_dt1Variable;
  using Base::a_dt2Variable;
  using Base::a_externalSource;
  using Base::a_FcellVolume;
  using Base::a_isBndryCell;
  using Base::a_reconstructionData;
  using Base::a_reconstructionNeighborId;
  using Base::a_resetPropertiesSolver;
  using Base::a_wasGapCell;
  using Base::applyInitialCondition;
  using Base::azimuthalNearBoundaryExchange;
  using Base::buildLeastSquaresStencilSimple;
  using Base::c_coordinate;
  using Base::c_level;
  using Base::checkCells;
  using Base::checkNeighborActive;
  using Base::collectParameters;
  using Base::collectVariables;
  using Base::computeCellSurfaceDistanceVectors;
  using Base::computeCellVolumes;
  using Base::computeConservativeVariables;
  using Base::computePV;
  using Base::computeRecConstSVD;
  using Base::computeVolumeForces;
  using Base::createBoundaryCells;
  using Base::cutOffBoundaryCondition;
  using Base::determineStructuredCells;
  using Base::domainOffset;
  using Base::exchangeAll;
  using Base::exchangeFloatDataAzimuthal;
  using Base::extractPointIdsFromGrid;
  using Base::forceTimeStep;
  using Base::getAdjacentLeafCells;
  using Base::getVorticity;
  using Base::initAzimuthalCartesianHaloInterpolation;
  using Base::initAzimuthalMaxLevelExchange;
  using Base::initAzimuthalNearBoundaryExchange;
  using Base::initAzimuthalReconstruction;
  using Base::initCutOffBoundaryCondition;
  using Base::initViscousFluxComputation;
  using Base::m_A;
  using Base::m_adaptationDampingDistance;
  using Base::m_adaptationLevel;
  using Base::m_adaptationSinceLastRestart;
  using Base::m_adaptationSinceLastRestartBackup;
  using Base::m_angle;
  using Base::m_associatedInternalCells;
  using Base::m_averageVorticity;
  using Base::m_azimuthalBndrySide;
  using Base::m_azimuthalCutRecCoord;
  using Base::m_azimuthalHaloActive;
  using Base::m_azimuthalMaxLevelHaloCells;
  using Base::m_azimuthalMaxLevelWindowCells;
  using Base::m_azimuthalMaxLevelWindowMap;
  using Base::m_azimuthalNearBndryInit;
  using Base::m_azimuthalNearBoundaryBackupMaxCount;
  using Base::m_azimuthalRecConsts;
  using Base::m_azimuthalRecConstSet;
  using Base::m_azimuthalReconstNghbrIds;
  using Base::m_azimuthalRemappedNeighborDomains;
  using Base::m_azimuthalRemappedWindowCells;
  using Base::m_bandWidth;
  using Base::m_bndryLevelJumps;
  using Base::m_bodyIdOutput;
  using Base::m_cellSurfaceMapping;
  using Base::m_computeViscousFlux;
  using Base::m_currentGridFileName;
  using Base::m_dataBlockSize;
  using Base::m_engineSetup;
  using Base::m_externalSourceDt1;
  using Base::m_extractedCells;
  using Base::m_forceRestartGrid;
  using Base::m_forcing;
  using Base::m_gridPoints;
  using Base::m_hasExternalSource;
  using Base::m_levelSetOutput;
  using Base::m_linerLvlJump;
  using Base::m_lsCutCellBaseLevel;
  using Base::m_masterCellIds;
  using Base::m_maxIterations;
  using Base::m_maxLevelHaloCells;
  using Base::m_maxLevelWindowCells;
  using Base::m_maxNearestBodies;
  using Base::m_maxNoAzimuthalRecConst;
  using Base::m_maxNoSets;
  using Base::m_movingAvgInterval;
  using Base::m_multilevel;
  using Base::m_multipleFvSolver;
  using Base::m_noAzimuthalReconstNghbrs;
  using Base::m_noCellsInsideSpongeLayer;
  using Base::m_noCorners;
  using Base::m_noEdges;
  using Base::m_noMaxLevelHaloCells;
  using Base::m_noMaxLevelWindowCells;
  using Base::m_nonBlockingComm;
  using Base::m_noOuterBndryCells;
  using Base::m_noSamples;
  using Base::m_noTimeStepsBetweenSamples;
  using Base::m_outputFormat;
  using Base::m_recalcIds;
  using Base::m_receiveBuffers;
  using Base::m_receiveBuffersNoBlocking;
  using Base::m_reconstructionCellIds;
  using Base::m_reconstructionConstants;
  using Base::m_reconstructionNghbrIds;
  using Base::m_reconstructSurfaceData;
  using Base::m_residualInterval;
  using Base::m_restartFile;
  using Base::m_restartFileOutputTimeStep;
  using Base::m_restartOffset;
  using Base::m_restartOldVariables;
  using Base::m_restartTimeBc2800;
  using Base::m_restartTimeStep;
  using Base::m_revDir;
  using Base::m_rhoEInfinity;
  using Base::m_rhoVInfinity;
  using Base::m_RKalpha;
  using Base::m_RKStep;
  using Base::m_rotIndVarsCV;
  using Base::m_rotIndVarsPV;
  using Base::m_saveVorticityToRestart;
  using Base::m_sendBuffers;
  using Base::m_sendBuffersNoBlocking;
  using Base::m_setToBodiesTable;
  using Base::m_sigmaSponge;
  using Base::m_singleAdaptation;
  using Base::m_smallCellIds;
  using Base::m_solutionInterval;
  using Base::m_solutionTimeSteps;
  using Base::m_splitSurfaces;
  using Base::m_spongeLayerThickness;
  using Base::m_static_advanceSolution_dragCnt;
  using Base::m_static_advanceSolution_firstRun;
  using Base::m_static_advanceSolution_meanDrag;
  using Base::m_static_advanceSolution_meanDragCoeff;
  using Base::m_static_applyBoundaryCondition_ERhoL1;
  using Base::m_static_applyBoundaryCondition_ERhoL2;
  using Base::m_static_applyBoundaryCondition_ERhoLoo;
  using Base::m_static_applyBoundaryCondition_EVelL1;
  using Base::m_static_applyBoundaryCondition_EVelL2;
  using Base::m_static_applyBoundaryCondition_EVelLoo;
  using Base::m_static_applyBoundaryCondition_firstRun;
  using Base::m_static_applyBoundaryCondition_oldMass;
  using Base::m_static_applyBoundaryCondition_oldVol2;
  using Base::m_static_applyBoundaryCondition_refMass;
  using Base::m_static_constructGFieldPredictor_adaptiveGravity;
  using Base::m_static_constructGFieldPredictor_firstRun;
  using Base::m_static_logCellxd_firstRun;
  using Base::m_static_logData_firstRun4;
  using Base::m_static_logData_ic45299_amplitude;
  using Base::m_static_logData_ic45299_cutOffAngle;
  using Base::m_static_logData_ic45299_freqFactor;
  using Base::m_static_logData_ic45299_maxA;
  using Base::m_static_logData_ic45299_maxF;
  using Base::m_static_logData_ic45299_xCutOff;
  using Base::m_static_logData_ic45301_containingCellIds;
  using Base::m_static_logData_ic45301_first;
  using Base::m_static_logData_ic45301_freqFactor;
  using Base::m_static_logData_ic45301_maxF;
  using Base::m_static_logData_ic45301_noPressurePoints;
  using Base::m_static_logData_ic45301_pressurePoints;
  using Base::m_static_logData_ic45301_Strouhal;
  using Base::m_static_redistributeMass_firstRun;
  using Base::m_static_saveSolverSolutionxd_firstRun;
  using Base::m_static_updateBodyProperties_c455_firstRun;
  using Base::m_static_updateSpongeLayer_first;
  using Base::m_static_updateSpongeLayer_mbSpongeLayer;
  using Base::m_static_writeVtkXmlFiles_firstCall;
  using Base::m_static_writeVtkXmlFiles_firstCall2;
  using Base::m_surfaceVarMemory;
  using Base::m_sutherlandConstant;
  using Base::m_sutherlandConstantThermal;
  using Base::m_sutherlandPlusOneThermal;
  using Base::m_targetDensityFactor;
  using Base::m_timerGroup;
  using Base::m_timers;
  using Base::m_TInfinity;
  using Base::m_totalnoghostcells;
  using Base::m_useNonSpecifiedRestartFile;
  using Base::m_variablesName;
  using Base::m_vorticityName;
  using Base::m_vorticitySize;
  using Base::m_vtuCoordinatesThreshold;
  using Base::m_vtuDomainIdOutput;
  using Base::m_vtuGeometryOutput;
  using Base::m_vtuGeometryOutputExtended;
  using Base::m_vtuGlobalIdOutput;
  using Base::m_vtuLevelSetOutput;
  using Base::m_vtuLevelThreshold;
  using Base::m_vtuVelocityGradientOutput;
  using Base::m_vtuWriteGeometryFile;
  using Base::m_writeOutData;
  using Base::maxLevel;
  using Base::maxRefinementLevel;
  using Base::minCell;
  using Base::neighborDomain;
  using Base::noHaloCells;
  using Base::noMinCells;
  using Base::noVariables;
  using Base::noWindowCells;
  using Base::outputDir;
  using Base::rebuildAzimuthalReconstructionConstants;
  using Base::reduceData;
  using Base::reduceVariables;
  using Base::restartDir;
  using Base::saveGridFlowVarsPar;
  using Base::setInfinityState;
  using Base::setRestartFileOutputTimeStep;
  using Base::setTimeStep;
  using Base::setUpwindCoefficient;
  using Base::solverId;
  using Base::tagCellsNeededForSurfaceFlux;
  using Base::windowCellId;
  using Base::writeListOfActiveFlowCells;

  using Timers = maia::fv::Timers_;

  using Base::computeTimeStepEulerDirectional;
  using Base::m_reComputedBndry;
  using Base::m_resetInitialCondition;

  using Base::a_localTimeStep;
  using Base::m_localTS;
  using Base::m_refineDiagonals;
  using Base::reInitActiveCellIdsMemory;

  //================================== WMLES ======================================
  using Base::m_wmDistance;
  using Base::m_wmDomainId;
  using Base::m_wmGlobalNoSrfcProbeIds;
  using Base::m_wmLES;
  using Base::m_wmNoDomains;
  using Base::m_wmOutput;
  using Base::m_wmSurfaceProbeInterval;
  using Base::m_wmUseInterpolation;

  using Base::m_wmIterator;

  using Base::m_comm_wm;

  using Base::m_mpi_wmRecvReq;
  using Base::m_mpi_wmRequest;
  using Base::m_mpi_wmSendReq;
  using Base::m_noWMImgPointsRecv;
  using Base::m_noWMImgPointsSend;
  using Base::m_wmImgRecvBuffer;
  using Base::m_wmImgRecvIdMap;
  using Base::m_wmImgSendBuffer;
  using Base::m_wmLocalNoSrfcProbeIds;
  using Base::m_wmSrfcProbeRecvBuffer;
  using Base::m_wmSrfcProbeSendBuffer;
  using Base::m_wmSurfaces;

  using Base::m_wmImgCellIds;
  using Base::m_wmImgCoords;
  using Base::m_wmImgWMSrfcIds;

  using Base::m_wmSurfaceProbeIds;
  using Base::m_wmSurfaceProbeSrfcs;

  using Base::m_sweptVolume;
  using Base::m_sweptVolumeBal;

  // zonal method
  using Base::determineLESAverageCells;
  using Base::m_LESAverageCells;
  using Base::m_LESNoVarAverage;
  using Base::m_LESValues;
  using Base::m_LESVarAverage;
  using Base::m_LESVarAverageBal;
  using Base::m_noLESVariables;
  using Base::m_noRANSVariables;
  using Base::m_RANSValues;
  using Base::m_zonal;
  using Base::m_zonalAveragingTimeStep;
  using Base::resetZonalLESAverage;
  using Base::resetZonalSolverData;

  KDtree<nDim>* m_bodyTree = nullptr;
  KDtree<nDim>* m_bodyTreeLocal = nullptr;
  std::vector<Point<nDim>> m_bodyCenters;
  std::vector<Point<nDim>> m_bodyCentersLocal;

  friend FvBndryCndXD<nDim, SysEqn>;
  FvMbCartesianSolverXD(MInt, MInt, const MBool*, maia::grid::Proxy<nDim>& gridProxy_, Geometry<nDim>& geometry_,
                        const MPI_Comm comm);
  ~FvMbCartesianSolverXD() override;


  // destructor / constructor
  void initializeMb();
  void initAnalyticalLevelSet();
  void allocateBodyMemory(MInt mode = 0);

  template <class X = void, std::enable_if_t<nDim == 2, X*> = nullptr>
  static MInt locatenear(const Point<nDim>& /*pt*/, MFloat /*r*/, MInt* /*list*/, MInt /*nmax*/,
                         MBool returnCellId = true) {
    std::ignore = returnCellId;
    return 0;
  }

  template <class X = void, std::enable_if_t<nDim == 3, X*> = nullptr>
  MInt locatenear(const Point<nDim>& pt, MFloat r, MInt* list, MInt nmax, MBool returnCellId = true) {
    return m_bodyTree->locatenear(pt, r, &list[0], nmax, returnCellId);
  }

  std::vector<MFloat> m_oldBodyPosition;

  std::vector<MFloat> m_movingPosDiff;

  MInt m_tCutGroup{};
  MInt m_tCutGroupTotal{};

  // settings
  MBool m_trackMovingBndry{};
  MBool m_logBoundaryData{};
  MInt m_motionEquation{};
  MInt m_movingBndryCndId{};
  MInt m_loggingInterval{};
  MInt m_bodySamplingInterval{};
  MInt m_particleSamplingInterval{};
  MInt m_onlineRestartInterval{};
  MInt m_maxBndryLayerLevel{};
  MInt m_maxBndryLayerWidth{};
  MInt m_trackMbStart{};
  MInt m_FSIStart{};
  MInt m_trackMbEnd{};
  MInt m_bodyTypeMb{};
  MFloat m_FSIToleranceU{};
  MFloat m_FSIToleranceX{};
  MFloat m_FSIToleranceW{};
  MFloat m_FSIToleranceR{};
  MFloat m_outsideFactor{};
  MInt m_killSwitchCheckInterval{};
  MInt m_solutionDiverged{};
  SolverType m_solverType;
  MBool m_maxLevelDecrease = false;
  MFloat m_Fr{};
  MFloat m_g{};
  MBool m_haloCellOutput{};
  MBool m_complexBoundary{};
  MBool m_conservationCheck{};
  MBool m_writeCenterLineData{};
  MBool m_gclIntermediate{};
  MBool m_onlineRestart{};
  MFloat m_physicalTimeDt1{};
  MInt m_timeStepAdaptationStart{};
  MInt m_timeStepAdaptationEnd{};
  MFloat m_cflInitial{};
  MFloat m_cflTarget{};
  MFloat m_previousTimeStep{};
  MInt m_centralizeSurfaceVariables{};
  MBool m_generateOuterBndryCells{};
  MBool m_applyCollisionModel{};
  MBool m_applyRotationalCollisionModel{};
  MBool m_LsMovement{};
  MBool m_alwaysResetCutOff{};

  MBool m_standardRK = false;
  MInt m_reConstSVDWeightMode{};

  MInt m_lsCutCellMinLevel;
  MInt* m_lsCutCellLevel = nullptr;
  MBool m_maxLevelChange = false;
  MBool m_dynamicStencil{};

  std::map<MInt, MFloat> m_oldGeomBndryCells;
  MBool* m_geometryChange = nullptr;

  // Gap-Handling
  std::vector<MInt>* m_gapWindowCells = nullptr;
  std::vector<MInt>* m_gapHaloCells = nullptr;
  MBool m_gapCellExchangeInit = false;
  MBool m_noneGapRegions = true;
  MBool m_initGapCell = true;
  MInt* m_gapState = nullptr;


  // offsets & counters
  MInt m_noFlowCells{};
  MInt m_noGCells{};
  MInt m_noSurfaces{};
  MInt m_noBndryCells{};
  MInt m_maxNoEmbeddedBodiesPeriodic{};
  MInt m_initialSurfacesOffset{};
  MInt m_noOuterBoundarySurfaces{};
  MInt m_noLsMbBndryCells{};
  MInt m_noEmergedCells{};
  MInt m_noEmergedWindowCells{};
  MInt m_noBndryCandidates{};
  MInt m_noAzimuthalBndryCandidates{};
  MInt m_noNearBndryCells{};
  MInt m_reconstructionOffset{};
  MInt m_structureStep{};
  MInt m_maxStructureSteps{};
  MInt m_noSurfacePointSamples{};
  MInt m_maxNoSurfacePointSamples{};
  MInt m_noPointParticles{};
  MInt m_noPointParticlesLocal{};
  MInt* m_sampleNghbrs = nullptr;
  MInt* m_sampleNghbrOffsets = nullptr;
  MBool m_gapOpened{};
  MBool m_trackBodySurfaceData = true;
  MBool m_buildCollectedLevelSetFunction = false;
  MInt m_startSet{};
  // MInt m_noSets{};
  MInt* m_noBodiesInSet = nullptr;
  MInt* m_bodyToSetTable = nullptr;

  // constants and pseudo constants
  const clock_t m_t0;
  static constexpr MFloat FD = nDim;
  static constexpr MInt m_noCellNodes = IPOW2(nDim);
  const MInt m_noCVars;
  const MInt m_noFVars;
  const MInt m_noPVars;
  const MInt m_noSlopes;
  const MFloat m_gammaMinusOne;
  MInt m_maxNoSampleNghbrs{};
  MFloat m_U2{};
  MFloat m_rhoU2{};
  MFloat m_bodyDistThreshold{};
  MFloat m_periodicGhostBodyDist{};
  MFloat* m_gridCellArea = nullptr;
  MFloat* m_gravity = nullptr;

  // body properties
  MFloat* m_bodyRadius = nullptr;
  MFloat m_minBodyRadius{};
  MFloat m_maxBodyRadius{};
  MFloat* m_bodyRadii = nullptr;
  MFloat* m_bodyDiameter = nullptr;
  MFloat* m_bodyVolume = nullptr;
  MFloat* m_bodyMass = nullptr;
  MFloat* m_projectedArea = nullptr;
  MFloat m_addedMassCoefficient{};
  MFloat m_densityRatio{};
  MFloat m_capacityConstantVolumeRatio{};
  MFloat m_couplingRate{};
  MFloat** m_coupling = nullptr;
  MFloat m_couplingRateLinear{};
  MFloat m_couplingRatePressure{};
  MFloat m_couplingRateViscous{};
  MFloat* m_bodyQuaternion = nullptr;
  MFloat* m_bodyQuaternionDt1 = nullptr;
  MFloat* m_bodyDensity = nullptr;
  MFloat* m_bodyMomentOfInertia = nullptr;
  MFloat* m_bodyAccelerationDt1 = nullptr;
  MFloat* m_bodyAccelerationDt2 = nullptr;
  MFloat* m_bodyAccelerationDt3 = nullptr;
  MFloat* m_bodyVelocityDt1 = nullptr;
  MFloat* m_bodyVelocityDt2 = nullptr;
  MFloat* m_bodyTorque = nullptr;
  MFloat* m_bodyTorqueDt1 = nullptr;
  MFloat* m_bodyForce = nullptr;
  MFloat* m_bodyForceDt1 = nullptr;
  MFloat* m_bodyAngularVelocityDt1 = nullptr;
  MFloat* m_bodyAngularAccelerationDt1 = nullptr;
  MFloat* m_hydroForce = nullptr;
  MFloat* m_bodyCenterDt1 = nullptr;
  MFloat* m_bodyCenterDt2 = nullptr;
  MFloat* m_bodyCenterInitial = nullptr;
  MFloat* m_bodyTerminalVelocity = nullptr;
  MInt* m_fixedBodyComponents = nullptr;
  MInt* m_fixedBodyComponentsRotation = nullptr;
  MFloat* m_bodyNeutralCenter = nullptr;
  MInt* m_bodyEquation = nullptr;
  MInt* m_bodyInCollision = nullptr;
  MBool* m_bodyNearDomain = nullptr;
  MFloat* m_bodyDataAverage = nullptr;
  MFloat* m_bodyDataAverage2 = nullptr;
  MFloat* m_bodyDataDevAverage = nullptr;
  MFloat* m_bodySumAverage = nullptr;
  MBool* m_bodyDataCollision = nullptr;

  MBool m_forceAdaptation = false;

  // Lagrange particles
  MInt m_pointParticleTwoWayCoupling{};
  MInt m_pointParticleType{};
  MFloat* m_particleTerminalVelocity = nullptr;
  std::vector<MFloat> m_particleCoords;
  std::vector<MFloat> m_particleCoordsDt1;
  std::vector<MFloat> m_particleCoordsInitial;
  std::vector<MFloat> m_particleVelocity;
  std::vector<MFloat> m_particleTemperature;
  std::vector<MFloat> m_particleTemperatureDt1;
  std::vector<MFloat> m_particleHeatFlux;
  std::vector<MFloat> m_particleVelocityDt1;
  std::vector<MFloat> m_particleVelocityFluid;
  std::vector<MFloat> m_particleFluidTemperature;
  std::vector<MFloat> m_particleAcceleration;
  std::vector<MFloat> m_particleAccelerationDt1;
  std::vector<MFloat> m_particleQuaternions;
  std::vector<MFloat> m_particleQuaternionsDt1;
  std::vector<MFloat> m_particleAngularVelocity;
  std::vector<MFloat> m_particleAngularVelocityDt1;
  std::vector<MFloat> m_particleVelocityGradientFluid;
  std::vector<MFloat> m_particleAngularAcceleration;
  std::vector<MFloat> m_particleAngularAccelerationDt1;
  std::vector<MFloat> m_particleShapeParams;
  std::vector<MFloat> m_particleRadii;
  std::vector<MInt> m_particleCellLink;
  std::vector<MInt> m_particleGlobalId;
  std::vector<MInt>* m_periodicGhostBodies = nullptr;
  std::set<MInt> m_bodyWasInCollision;
  std::set<MInt> m_bodyWasActuallyInCollision;

  // collectors
  std::vector<MInt> m_bndryCandidates;
  std::vector<CutCandidate<nDim>> m_cutCandidates;
  std::vector<MInt> m_bndryCandidateIds;
  MBool m_wasBalanced = false;
  std::vector<MInt> m_azimuthalBndryCandidates;
  std::vector<MInt> m_azimuthalBndryCandidateIds;
  std::vector<MFloat> m_candidateNodeValues;
  std::vector<MInt> m_candidateNodeSet;
  std::vector<MInt> m_bndryLayerCells;
  MFloat* m_bodyReducedVelocity = nullptr;
  MFloat* m_bodyReducedMass = nullptr;
  MFloat* m_bodyReducedFrequency = nullptr;
  MFloat* m_bodyDampingCoefficient = nullptr;
  MFloat** m_maxBndryLayerDistances = nullptr;
  MFloat* m_volumeFraction = nullptr;
  MFloat m_stableVolumeFraction{};
  MFloat* m_stableCellVolume = nullptr;
  // MFloat* m_sweptVolume = nullptr;
  MFloat* m_sweptVolumeDt1 = nullptr;
  MFloat** m_rhsViscous = nullptr;
  std::vector<MInt>* m_linkedWindowCells = nullptr;
  std::vector<MInt>* m_linkedHaloCells = nullptr;

  MFloat* m_cellVolumesDt1 = nullptr;
  MFloat* m_cellVolumesDt2 = nullptr;
  MFloat* m_boundaryPressureDt1 = nullptr;
  MFloat* m_boundaryDensityDt1 = nullptr;
  MBool** m_pointIsInside = nullptr;
  MInt* m_noCutCellFaces = nullptr;
  MInt** m_noCutFacePoints = nullptr;
  MInt** m_cutFacePointIds = nullptr;
  MFloat** m_cutFaceArea = nullptr;
  MFloat* m_sampleCoordinates = nullptr;
  MFloat* m_sampleNormals = nullptr;
  MFloat* m_bbox = nullptr;
  MFloat* m_bboxLocal = nullptr;
  MInt* m_nghbrList = nullptr;
  //  MInt* m_periodicDir;
  MPI_Request* g_mpiRequestMb = nullptr;
  MInt* m_sendBufferSize = nullptr;
  MInt* m_receiveBufferSize = nullptr;
  std::set<MInt> m_freeSurfaceIndices;

  std::map<MInt, MFloat> m_oldBndryCells;
  std::map<MInt, std::vector<MFloat>> m_nearBoundaryBackup;
  std::vector<std::tuple<MInt, MInt, MInt>> m_temporarilyLinkedCells;
  std::vector<MInt> m_massRedistributionIds;
  std::vector<MFloat> m_massRedistributionVariables;
  std::vector<MFloat> m_massRedistributionRhs;
  std::vector<MFloat> m_massRedistributionVolume;
  std::vector<MFloat> m_massRedistributionSweptVol;
  std::vector<MInt> m_refOldBndryCells;
  std::set<MInt> m_coarseOldBndryCells;

  std::vector<std::vector<MInt>> m_azimuthalLinkedHaloCells;
  std::vector<std::vector<MInt>> m_azimuthalLinkedWindowCells;
  std::multimap<MLong, std::pair<std::vector<MFloat>, std::vector<MUlong>>> m_azimuthalNearBoundaryBackup;
  std::vector<MInt> m_azimuthalWasNearBndryIds;
  MFloat* m_azimuthalNearBoundaryBackupBalFloat;
  MLong* m_azimuthalNearBoundaryBackupBalLong;
  MInt m_noFloatDataBalance = 0;
  MInt m_noFloatDataAdaptation = 0;
  MInt m_noLongData = 0;
  MInt m_noLongDataBalance = 0;
  const MFloat m_azimuthalCornerEps = 0.01; // Used to determine if recCoords match

  // array pointers
  MFloat* m_slope = nullptr;
  MFloat* m_area = nullptr;
  MFloat* m_surfaceVariables = nullptr;
  MFloat* m_cellCoordinates = nullptr;
  MFloat* m_surfaceCoordinates = nullptr;
  MFloat* m_cellSlopes = nullptr;
  MFloat* m_vars = nullptr;
  MFloat* m_oldVars = nullptr;

  // post processing
  MBool m_firstStats = true;
  static constexpr MInt m_noMeanStatistics = 6;
  static constexpr MFloat m_distThresholdStat[m_noMeanStatistics] = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0};
  MFloat m_oldMeanState[m_noMeanStatistics][4]{};
  MBool m_printHeaderCorr = true;
  MBool m_printHeaderSlip = true;
  static constexpr MInt m_noAngleSetups = 3;
  static constexpr MInt m_noDistSetups = 3;
  MBool m_saveSlipData;
  MInt m_slipInterval;
  MInt m_saveSlipInterval;
  MInt m_noSlipDataOutputs;
  std::vector<MInt> m_particleOffsets;
  MInt* m_slipDataParticleCollision = nullptr;
  MFloat* m_slipDataParticleQuaternion = nullptr;
  MFloat* m_slipDataParticleVel = nullptr;
  MFloat* m_slipDataParticleAngularVel = nullptr;
  MFloat* m_slipDataParticleForce = nullptr;
  MFloat* m_slipDataParticleTorque = nullptr;
  MFloat* m_slipDataParticleFluidVel = nullptr;
  MFloat* m_slipDataParticleFluidVelGrad = nullptr;
  MFloat* m_slipDataParticleFluidVelRot = nullptr;
  MFloat* m_slipDataParticleFluidVelRnd = nullptr;
  MFloat* m_slipDataParticleFluidVelGradRnd = nullptr;
  MFloat* m_slipDataParticleFluidVelRotRnd = nullptr;
  MFloat* m_slipDataParticlePosition = nullptr;
  std::vector<MFloat> m_slipDataTimeSteps;

  // timw-engine specific:
  MFloat* m_gapAngleClose = nullptr;
  MFloat* m_gapAngleOpen = nullptr;
  MFloat* m_gapSign = nullptr;
  MInt m_forceNoGaps;


  // global functions
  void setLevelSetMbCellProperties();
  void writeListOfActiveFlowCells() override;
  void setTimeStep() override;
  void initSolver() override;
  void finalizeInitSolver() override;
  void initSolutionStep(MInt mode) override;
  void reInitSolutionStep(const MInt mode);
  void updateInfinityVariables();
  void advanceSolution();
  void advanceBodies();
  void advanceTimeStep();
  void prepareNextTimeStep() override;
  void computeBoundarySurfaceForces();
  void applyInitialCondition() override;
  void applyExternalSource() override;
  void applyExternalOldSource() override;
  void advanceExternalSource() override;
  MBool maxResidual(MInt mode = 0) override;
  MLong getNumberOfCells(MInt mode);

  // initialisation
  void initGField();
  void initBndryLayer();

  void updateGeometry();
  void constructGField(MBool updateBody = true);
  MInt constructDistance(const MFloat deltaMax, MIntScratchSpace& nearestBodies, MFloatScratchSpace& nearestDist);
  void transferLevelSetFieldValues(MBool);
  void updateLevelSetOutsideBand();
  MFloat CalculateLSV(MInt, MIntScratchSpace&);
  void setParticleFluidVelocities(MFloat* = nullptr);
  void advancePointParticles();
  MBool constructGFieldCorrector();
  void constructGFieldPredictor();
  void initBodyProperties();
  void initPointParticleProperties();
  void initBodyVelocities();
  void updateBodyProperties();
  void transferBodyState(MInt k, MInt b);
  void createPeriodicGhostBodies();
  template <class X = void, std::enable_if_t<nDim == 3, X*> = nullptr>
  void createBodyTree();
  template <class X = void, std::enable_if_t<nDim == 2, X*> = nullptr>
  void createBodyTree();
  void computeBodyMomentOfInertia();
  static void computeRotationMatrix(MFloatScratchSpace&, MFloat* q);
  static void matrixProduct(MFloatScratchSpace& C, MFloatScratchSpace& A, MFloatScratchSpace& B);
  static void matrixProductTranspose(MFloatScratchSpace& C, MFloatScratchSpace& A, MFloatScratchSpace& B);
  static void matrixProductTranspose2(MFloatScratchSpace& C, MFloatScratchSpace& A, MFloatScratchSpace& B);
  static void matrixVectorProduct(MFloat* c, MFloatScratchSpace& A, MFloat* b);
  static void matrixVectorProductTranspose(MFloat* c, MFloatScratchSpace& A, MFloat* b);
  void integrateBodyRotation();
  void getBodyRotation(const MInt bodyId, MFloat* bodyRotation);
  void getBodyRotationDt1(const MInt bodyId, MFloat* bodyRotation);
  void setBodyQuaternions(const MInt bodyId, MFloat* bodyRotation);

  void setGapOpened();
  void gapHandling();
  void initGapCellExchange();
  void gapCellExchange(MInt);
  MInt checkNeighborActivity(MInt);
  void updateGapBoundaryCells();
  void checkGapCells();
  void interpolateGapBodyVelocity();
  MFloat interpolateLevelSet(MInt*, MFloat*, MInt);
  void setGapCellId();

  // boundary generation
  void logBoundaryData(const MChar* fileName, MBool forceOutput);
  void generateBndryCellsMb(const MInt mode);
  void checkHaloBndryCells(MBool sweptVol);
  void determineBndryCandidates();
  void computeNodalLSValues();
  void exchangeNodalLSValues();

  void initAzimuthalLinkedHaloExc();
  void exchangeLinkedHaloCellsForAzimuthalReconstruction();
  void storeAzimuthalPeriodicData(MInt mode = 0);
  void commAzimuthalPeriodicData(MInt mode = 0);
  void computeAzimuthalReconstructionConstants(MInt mode = 0) override;
  void buildAdditionalAzimuthalReconstructionStencil(MInt mode = 0);
  void updateAzimuthalMaxLevelRecCoords();
  void updateAzimuthalNearBoundaryExchange();
  void computeAzimuthalHaloNodalValues();
  void exchangeAzimuthalOuterNodalValues(std::vector<CutCandidate<nDim>>& candidates, std::vector<MInt>& candidateIds);

  void computeCutPointsMb_MGC();
  MBool createCutFaceMb();
  template <class X = void, std::enable_if_t<nDim == 2, X*> = nullptr>
  void createCutFaceMb_MGC();
  template <class X = void, std::enable_if_t<nDim == 3, X*> = nullptr>
  void createCutFaceMb_MGC();
  template <class X = void, std::enable_if_t<nDim == 3, X*> = nullptr>
  void checkNormalVectors();
  void determineNearBndryCells();
  void computeGhostCellsMb();
  void setOuterBoundaryGhostCells();
  void setOuterBoundarySurfaces();
  void setupBoundaryCandidatesAnalytical();
  void checkBoundaryCells();
  void linkBndryCells();
  void redistributeMass();
  MInt createSplitCell(MInt cellId, MInt noSplitChilds);
  MFloat updateCellVolumeGCL(MInt bndryId);
  void correctRefinedBndryCell();
  void correctRefinedBndryCellVolume();
  void correctCoarsenedBndryCellVolume();

  void rebuildKDTreeLocal();
  void rebuildKDTree();

  // reinitialization routines
  void resetSolverFull() override;
  void resetSolverMb();
  void resetSolver() override;
  void resetSurfaces() override;
  void balance(const MInt* const noCellsToReceiveByDomain, const MInt* const noCellsToSendByDomain,
               const MInt* const targetDomainsByCell, const MInt oldNoCells) override;
  void balancePre() override;
  void balancePost() override;
  void finalizeBalance() override;

  MBool hasSplitBalancing() const override { return true; }
  void localToGlobalIds() override;

  MInt noCellDataDlb() const override {
    if(grid().isActive()) {
      return CellDataDlb::count;
    } else {
      return 0;
    }
  };
  MInt cellDataTypeDlb(const MInt dataId) const override {
    MInt dataType = -1;
    if(dataId < CellDataDlb::count) {
      dataType = s_cellDataTypeDlb[dataId];
    } else {
      TERMM(1, "solverCelldataType: invalid data id " + std::to_string(dataId));
    }
    return dataType;
  };
  MInt cellDataSizeDlb(const MInt dataId, const MInt gridCellId) override;

  /// Return solver data for DLB
  void getCellDataDlb(const MInt dataId, const MInt oldNoCells, const MInt* const bufferIdToCellId,
                      MFloat* const data) override;
  void getCellDataDlb(const MInt dataId, const MInt oldNoCells, const MInt* const bufferIdToCellId, MLong* const data);
  /// Set solver data for DLB
  void setCellDataDlb(const MInt dataId, const MFloat* const data) override;
  void setCellDataDlb(const MInt dataId, const MLong* const data);

  MBool adaptationTrigger() override;
  void resetMbBoundaryCells();
  void deleteNeighbourLinks();
  void restoreNeighbourLinks();
  void initializeEmergedCells();
  void initEmergingGapCells();
  void updateCellSurfaceDistanceVectors();
  void updateCellSurfaceDistanceVector(MInt srfcId);
  void updateViscousFluxComputation();
  void leastSquaresReconstruction();
  void buildAdditionalReconstructionStencil();
  void computeVolumeFraction();
  void computeReconstructionConstants() override;

  void prepareAdaptation() override;
  void setSensors(std::vector<std::vector<MFloat>>& sensors,
                  std::vector<MFloat>& sensorWeight,
                  std::vector<std::bitset<64>>& sensorCellFlag,
                  std::vector<MInt>& sensorSolverId) override;
  void postAdaptation() override;
  void finalizeAdaptation() override;
  void sensorPatch(std::vector<std::vector<MFloat>>&, std::vector<std::bitset<64>>&, std::vector<MFloat>&, MInt,
                   MInt) override;

  void removeChilds(const MInt) override;
  void refineCell(const MInt) override;
  void setCellWeights(MFloat*) override;
  void descendLevelSetValue(const MInt cellId, const MInt* const bodyIds, const MInt bodyCnt);
  void descendDistance(const MInt cellId, const MInt* const bodyIds, const MInt bodyCnt,
                       MIntScratchSpace& nearestBodies, MFloatScratchSpace& nearestDist);
  void getBoundaryDistance(MFloatScratchSpace&) override;
  void swapCells(const MInt, const MInt) override;
  void moveSurface(const MInt toSrfcId, const MInt fromSrfcId);
  void compactSurfaces();
  void correctMasterSlaveSurfaces();
  void restoreSurfaces(const MInt cellId);
  void resetSurfaces(MInt cellId);
  void removeSurfaces(MInt cellId);
  void deleteSurface(MInt srfcId);
  void createSurface(MInt srfcId, MInt nghbrId0, MInt nghbrId1, MInt orientation);
  void createSurfaceSplit(MInt srfcId, MInt cellId, MInt splitChildId, MInt direction);
  MInt getNewSurfaceId();
  void rebuildReconstructionConstants(MInt cellId);
  void setCellProperties() override;
  static void setAdditionalCellProperties();
  void setAdditionalActiveFlag(MIntScratchSpace&) override;
  void computeLocalBoundingBox();

  // surface generation
  void createInitialSrfcs();
  void correctSrfcsMb();
  void addBoundarySurfacesMb();

  // MUSCL
  void Muscl(MInt timerId = -1) override;
  ATTRIBUTES1(ATTRIBUTE_HOT) void LSReconstructCellCenter() override;
  ATTRIBUTES1(ATTRIBUTE_HOT) void computeSurfaceValues(MInt timerId = -1) override;
  ATTRIBUTES1(ATTRIBUTE_HOT) void computeSurfaceValuesLimited(MInt timerId = -1) override;

  // base solver
  ATTRIBUTES1(ATTRIBUTE_HOT) void Ausm() override;
  ATTRIBUTES1(ATTRIBUTE_HOT) void viscousFlux() override;
  ATTRIBUTES1(ATTRIBUTE_HOT) void resetRHS() override;
  ATTRIBUTES1(ATTRIBUTE_HOT) void distributeFluxToCells() override;
  void initializeRungeKutta() override;
  void setRungeKuttaFunctionPointer();
  MBool rungeKuttaStep() override;
  MBool (FvMbCartesianSolverXD::*execRungeKuttaStep)();
  template <MInt timeStepMethod>
  ATTRIBUTES1(ATTRIBUTE_HOT)
  MBool rungeKuttaStepStandard();
  template <MInt timeStepMethod>
  ATTRIBUTES1(ATTRIBUTE_HOT)
  MBool rungeKuttaStepNew();
  ATTRIBUTES1(ATTRIBUTE_HOT) MBool rungeKuttaStepNewLocalTimeStepping();
  void computeTimeStep();
  void adaptTimeStep();
  void updateSplitParentVariables() override;
  void updateSpongeLayer() override;
  virtual void generateBndryCells();
  ATTRIBUTES1(ATTRIBUTE_HOT) void exchange() override;
  void initializeMaxLevelExchange() override;
  void updateMultiSolverInformation(MBool fullReset = false) override;
  void exchangeLevelSetData();
  void applyALECorrection();
  void resetRHSCutOffCells() override;

  // boundary conditions
  void applyBoundaryCondition() override;
  void applyBoundaryConditionMb();
  void applyBoundaryConditionSlope();
  void updateGhostCellSlopesInviscid();
  void updateGhostCellSlopesViscous();
  void applyDirichletCondition();
  void applyNeumannCondition();
  void setBoundaryVelocity();
  void copyVarsToSmallCells() override;
  void correctBoundarySurfaceVariablesMb();
  void computeBodySurfaceData(MFloat* pressureForce = nullptr);
  void resetRHSNonInternalCells() override;
  void updateLinkedCells();

  // input/output functions
  void saveSolverSolution(const MBool forceOutput = false, const MBool finalTimeStep = false) override;
  void recordBodyData(const MBool& firstRun);
  void crankAngleSolutionOutput();
  void writeStencil(const MInt bndryId);
  void computeForceCoefficients(MFloat*&);
  template <class X = void, std::enable_if_t<nDim == 2, X*> = nullptr>
  MInt writeGeometryToVtkXmlFile(const MString& fileName);
  template <class X = void, std::enable_if_t<nDim == 3, X*> = nullptr>
  static MInt writeGeometryToVtkXmlFile(const MString& fileName);
  void saveRestartFile(const MBool) override;
  void loadRestartFile() override;
  void loadBodyRestartFile(MInt);
  void setOldGeomBndryCellVolume();
  void saveBodyRestartFile(const MBool backup);
  void saveBodySamples();
  void saveParticleSamples();
  void writeVtkDebug(const MString suffix);
  template <class X = void, std::enable_if_t<nDim == 2, X*> = nullptr>
  void writeVtkXmlOutput(const MString& fileName, MBool debugOutput = false);
  template <class X = void, std::enable_if_t<nDim == 3, X*> = nullptr>
  void writeVtkXmlOutput(const MString& fileName, MBool debugOutput = false);
  template <class X = void, std::enable_if_t<nDim == 2, X*> = nullptr>
  void writeVtkXmlOutput_MGC(const MString& fileName);
  void writeVtkXmlFiles(const MString fileName, const MString GFileName, MBool regularOutput, MBool diverged) override;
  void writeVtkErrorFile();
  void printDynamicCoefficients(MBool firstRun, MFloat* pressureForce);
  std::vector<MInt> m_forceFVMBStatistics;
  void computeFlowStatistics(MBool force);
  void saveParticleSlipData();
  void computeSlipStatistics(const MIntScratchSpace& nearestBodies, const MFloatScratchSpace& nearestDist,
                             const MFloat maxDistConstructed);
  void computeNearBodyFluidVelocity(const MIntScratchSpace& nearestBodies, const MFloatScratchSpace& nearestDist,
                                    MFloat* vel, MFloat* velGrad, MFloat* rotation,
                                    const std::vector<MFloat>& setup = std::vector<MFloat>(),
                                    MInt* skipBodies = nullptr, MFloat* meanBodyState = nullptr,
                                    MFloat* pressure = nullptr, MFloat* velDt = nullptr, MFloat* rotationDt = nullptr,
                                    MFloat* nearestFac = nullptr);
  void computeBodyVolume(MFloatScratchSpace& partVol);
  void saveInitCorrData(const std::vector<MInt>& saveInitCorrDataTimeSteps, const MFloatScratchSpace& bodyVelocityInit,
                        const MFloatScratchSpace& bodyAngularVelocityInit, const MFloatScratchSpace& bodyQuaternionInit,
                        const MFloatScratchSpace& velMeanInit, const MFloatScratchSpace& velGradMeanInit,
                        const MFloatScratchSpace& particleFluidRotationMeanInit, const MIntScratchSpace& corrBodies,
                        const MFloatScratchSpace& bodyForceInit, const MFloatScratchSpace& bodyTorqueInit,
                        const MIntScratchSpace& bodyOffsets);
  MInt loadInitCorrData(const std::vector<MInt>& saveInitCorrDataTimeSteps, MFloatScratchSpace& bodyVelocityInit,
                        MFloatScratchSpace& bodyAngularVelocityInit, MFloatScratchSpace& bodyQuaternionInit,
                        MFloatScratchSpace& velMeanInit, MFloatScratchSpace& velGradMeanInit,
                        MFloatScratchSpace& particleFluidRotationMeanInit, MIntScratchSpace& corrBodies,
                        MFloatScratchSpace& bodyForceInit, MFloatScratchSpace& bodyTorqueInit,
                        const MIntScratchSpace& bodyOffsets);
  MFloat determineCoupling(MFloatScratchSpace& coupling);
  void writeCenterLineVel(const MChar* fileName);
  void readStlFile(const MChar* fileName, MBool readNormals);
  MBool prepareRestart(MBool, MBool&) override;
  void reIntAfterRestart(MBool) override;

  //  time-step and solution Step functions:
  void preTimeStep() override;
  void preSolutionStep(const MInt = -1) override;
  MBool solutionStep() override;
  MBool postSolutionStep() override;
  void postTimeStep() override;
  // postTimeStep is not overloaded from fv-cartesian-solver

  // help functions
  void logOutput();
  void resetSlopes();
  void logData();
  MInt createBndryCellMb(MInt cellId);
  MInt getAdjacentCells(MInt cellId, MInt* adjacentCells);
  MInt getAdjacentGridCells(MInt cellId, MInt* adjacentCells);
  MInt getAdjacentCellsAllLevels(MInt cellId, MInt* adjacentCells);
  MInt getAdjacentCellsExtended(MInt cellId, MInt* adjacentCells);
  MInt getEqualLevelNeighbors(MInt cellId, MInt (&nghbrs)[27]);
  MInt broadcastSignal(const MInt sender, const MInt signal);
  inline MInt getFacingNghbrs(const MInt cellId, const MBool includeAllChilds);
  inline MFloat temperature(MInt cellId);
  static inline void crossProduct(MFloat* c, MFloat* a, MFloat* b);
  inline MFloat filterFloat(MFloat number);
  static inline MFloat CdLaw(MFloat Re) { return 24.0 * (F1 + (0.15 * pow(Re, 0.687))) / Re; }
  MBool gridPointIsInside(MInt, MInt) override;
  MString getConservativeVarName(MInt i);
  MString getPrimitiveVarName(MInt i);
  MInt inverseGJ();
  template <typename T>
  MBool isNan(T val);

  MFloat getDistance(const MFloat* const, const MInt);
  MFloat getDistance(const MInt, const MInt);
  MFloat getDistanceSphere(const MFloat* const, const MInt);
  MFloat getDistanceSphere(const MInt, const MInt);
  MFloat getDistancePiston(const MFloat* const, const MInt);
  MFloat getDistancePiston(const MInt, const MInt);
  MFloat getDistanceEllipsoid(const MFloat* const, const MInt);
  MFloat getDistanceEllipsoid(const MInt, const MInt);
  MFloat getDistanceNaca00XX(const MInt, const MInt);
  MFloat getDistanceSplitSphere(const MInt, const MInt);
  MFloat getDistanceTetrahedron(const MFloat* const, const MInt);
  void getNormal(const MInt, const MInt, MFloat[]);
  void getNormalSphere(const MInt, const MInt, MFloat[]);
  void getNormalEllipsoid(const MInt, const MInt, MFloat[]);
  MFloat getLevelSetValueNaca00XX(const MFloat* const, const MFloat sign);
  MFloat distancePointEllipsoidSpecial(const MFloat e[3], const MFloat y[3], MFloat x[3]);
  MFloat distancePointEllipsoidSpecial2(const MFloat e[3], const MFloat y[3], MFloat x[3]);
  MFloat distancePointEllipsoid(const MFloat e[3], const MFloat y[3], MFloat x[3]);
  static MFloat distancePointEllipseSpecial(const MFloat e[2], const MFloat y[2], MFloat x[2]);
  MFloat distEllipsoidEllipsoid(const MInt k0, const MInt k1, MFloat* xc0, MFloat* xc1);

  void checkCellState(); // claudia
  inline MBool inside(MFloat x, MFloat a, MFloat b);
  MString printTime(const MFloat t);
  MBool fileExists(const MChar* fileName);
  MInt copyFile(const MChar* fromName, const MChar* toName);
  void checkHaloCells(MInt mode = 0);

  GeometryIntersection<nDim>* m_geometryIntersection = nullptr;
  FvBndryCndXD<nDim, SysEqn>* m_fvBndryCnd = nullptr;
  MInt returnNoActiveCorners(MInt);

 private:
  static inline MFloat sgn(MFloat val) { return (val < F0) ? -F1 : F1; }
  void checkDiv() override;

  class surfBase {
   public:
    MInt bodyId{};
    MFloat area{};
    MFloat normal[nDim]{};
    MFloat center[nDim]{};
    surfBase() = default;
    explicit surfBase(MInt _bodyId) : bodyId(_bodyId) {}
  };

  class polyVertex {
   public:
    MFloat coordinates[nDim]{};
    MInt pointId{};
    MInt pointType{}; // 0: corner point; 1: cut points; 2: additional MC vertex; 3: additional
                      // polyhedron clipping vertex
    MInt cartSrfcId{};
    std::vector<MInt> edges;
    std::set<MInt> surfaceIdentificators;
    polyVertex(MInt pId, MInt pType) : pointId(pId), pointType(pType), cartSrfcId(-1) {}
    polyVertex(MFloat* coords, MInt pId, MInt pType) : polyVertex(pId, pType) {
      std::copy(coords, coords + nDim, coordinates);
    }
  };

  class polyEdgeBase {
   public:
    MInt vertices[2];
    MInt edgeType; // 0: Cartesian uncut, 1: Cartesian cut, 2: cut line on face, 3: cut line in cell
    MInt edgeId;   // for edgeType 0/1: Cartesian edge, for edgeType 2/3: cartesian face or -1 if
                   // pure body line
    polyEdgeBase(MInt v0, MInt v1, MInt eId, MInt eType) : vertices{v0, v1}, edgeType(eType), edgeId(eId) {}
  };

  class polyEdge2D : public surfBase, public polyEdgeBase {
   public:
    MInt cutCell;
    MFloat w{};
    polyEdge2D(MInt v0, MInt v1, MInt eId, MInt eType, MInt _bodyId)
      : surfBase(_bodyId), polyEdgeBase(v0, v1, eId, eType), cutCell(-1) {}
  };

  class polyEdge3D : public polyEdgeBase {
   public:
    MInt face[2];
    polyEdge3D(MInt v0, MInt v1, MInt eId, MInt eType) : polyEdgeBase(v0, v1, eId, eType), face{-1, -1} {}
  };

  class polyFace : public surfBase {
   public:
    std::vector<std::pair<MInt, MInt>>
        edges; // contains 1. the edge and 2. the direction of the edge (1: as is, -1: reversed)
    MInt cutCell;
    MInt faceId{};
    MInt faceType{}; // 0: cartesian face; 1: body face
    MFloat w{};
    MInt tmpSetIndex;
    MInt isLine;
    polyFace(MInt fId, MInt fType, MInt _bodyId)
      : surfBase(_bodyId), cutCell(-1), faceId(fId), faceType(fType), tmpSetIndex(-1), isLine(false) {}

   private:
    polyFace() : surfBase(-1), cutCell(-1), tmpSetIndex(-1), isLine(false) {}
  };

  // 3D
  class polyMultiFace : public surfBase {
   public:
    std::vector<MInt> faces;
    std::list<std::pair<MInt, MInt>>
        edges; // contains 1. the edge and 2. the direction of the edge (1: as is, -1: reversed)
  };

  class polyCutCellBase {
   public:
    MInt cartesianCell;
    MFloat volume{};
    MFloat center[nDim]{};
    polyCutCellBase(MInt cartCell, MFloat* _center) : cartesianCell(cartCell) {
      std::copy(_center, _center + nDim, center);
    }
  };

  class polyCutCell : public polyCutCellBase {
   public:
    // The following naming is sick, but I don't care: in 2D these are edges and 3D faces
    std::vector<MInt> faces_edges;
    polyCutCell(MInt cartCell, MFloat* _center) : polyCutCellBase(cartCell, _center) {}
  };

  // 3D
  class splitCartesianFace {
   public:
    MInt direction;
    std::vector<MInt> srfcIds;
    explicit splitCartesianFace(MInt dir) : direction(dir) {}
  };

  // 3D
  class cellWithSplitFace {
   public:
    MInt cellId;
    std::vector<splitCartesianFace> splitFaces;
    explicit cellWithSplitFace(MInt id) : cellId(id) {}
  };


  MBool test_face(MInt, MInt, MInt);
  template <class X = void, std::enable_if_t<nDim == 2, X*> = nullptr>
  void writeVTKFileOfCell(MInt, const std::vector<polyEdge2D>*, const std::vector<polyVertex>*, MInt);
  template <class X = void, std::enable_if_t<nDim == 2, X*> = nullptr>
  void writeVTKFileOfCutCell(MInt, std::vector<polyCutCell>*, const std::vector<polyEdge2D>*,
                             const std::vector<polyVertex>*, MInt);
  template <class X = void, std::enable_if_t<nDim == 2, X*> = nullptr>
  void outputPolyData(const MInt, const std::vector<polyCutCell>*, const std::vector<polyEdge2D>*,
                      const std::vector<polyVertex>*, MInt);
  template <class X = void, std::enable_if_t<nDim == 3, X*> = nullptr>
  void outputPolyData(const MInt, const std::vector<polyCutCell>*, const std::vector<polyFace>*,
                      const std::vector<polyEdge3D>*, const std::vector<polyVertex>*, MInt);
  template <class X = void, std::enable_if_t<nDim == 2, X*> = nullptr>
  void compVolumeIntegrals(std::vector<polyCutCell>*, std::vector<polyEdge2D>*, const std::vector<polyVertex>*);
  template <class X = void, std::enable_if_t<nDim == 2, X*> = nullptr>
  void compFaceIntegrals(polyCutCell*, std::vector<polyEdge2D>*, const std::vector<polyVertex>*, MInt, MInt);
  template <class X = void, std::enable_if_t<nDim == 2, X*> = nullptr>
  void compProjectionIntegrals(polyCutCell*, std::vector<polyEdge2D>*, const std::vector<polyVertex>*, MInt, MInt,
                               MFloat*);


  /** \brief returns the normal corresponding to the triangle abc and returns the result in res
   *
   * \author Claudia Guenther, November 2013
   */
  template <class X = void, std::enable_if_t<nDim == 2, X*> = nullptr>
  static inline void computeNormal(const MFloat* p0, const MFloat* p1, MFloat* res, MFloat& w) {
    const MFloat dx = p1[0] - p0[0];
    const MFloat dy = p1[1] - p0[1];
    res[0] = -dy;
    res[1] = dx;
    const MFloat abs = sqrt(res[0] * res[0] + res[1] * res[1]);
    res[0] /= abs;
    res[1] /= abs;
    w = -res[0] * p0[0] - res[1] * p0[1];
  }


  /** \brief returns the normal corresponding to the triangle abc and returns the result in res
   *
   * \author Claudia Guenther, September 2013
   */
  template <class X = void, std::enable_if_t<nDim == 3, X*> = nullptr>
  static inline MFloat computeNormal(const MFloat* p0, const MFloat* p1, const MFloat* p2, MFloat* res, MFloat& w) {
    const MFloat a0 = p1[0] - p0[0], a1 = p1[1] - p0[1], a2 = p1[2] - p0[2], b0 = p2[0] - p0[0], b1 = p2[1] - p0[1],
                 b2 = p2[2] - p0[2];
    res[0] = a1 * b2 - a2 * b1;
    res[1] = a2 * b0 - a0 * b2;
    res[2] = a0 * b1 - a1 * b0;
    const MFloat abs = sqrt(res[0] * res[0] + res[1] * res[1] + res[2] * res[2]);
    res[0] /= abs;
    res[1] /= abs;
    res[2] /= abs;
    w = -res[0] * p0[0] - res[1] * p0[1] - res[2] * p0[2];
    return abs;
  }

  inline MFloat checkCentroidDiff(MInt srfcId1, MInt srfcId2);
  inline MFloat checkAreaDiff(MInt srfcId1, MInt srfcId2);


  class CsgVector {
   public:
    MFloat xx[3] = {F0, F0, F0};
    CsgVector() { std::fill(xx, xx + nDim, F0); }
    explicit CsgVector(const MFloat* X) { std::copy(X, X + nDim, xx); }
    CsgVector(const CsgVector& Y) { std::copy(Y.xx, Y.xx + nDim, xx); }
    CsgVector(CsgVector&&) noexcept = default;
    CsgVector& operator=(const CsgVector&) = default;
    CsgVector& operator=(CsgVector&&) = default;

    inline CsgVector clone() const {
      CsgVector tmp;
      std::copy(this->xx, this->xx + nDim, tmp.xx);
      return tmp;
    }
    inline CsgVector plus(const CsgVector& a) const {
      CsgVector tmp;
      for(MInt i = 0; i < nDim; i++)
        tmp.xx[i] = this->xx[i] + a.xx[i];
      return tmp;
    }
    inline CsgVector minus(const CsgVector& a) const {
      CsgVector tmp;
      for(MInt i = 0; i < nDim; i++)
        tmp.xx[i] = this->xx[i] - a.xx[i];
      return tmp;
    }
    inline CsgVector times(const MFloat a) const {
      CsgVector tmp;
      for(MInt i = 0; i < nDim; i++)
        tmp.xx[i] = this->xx[i] * a;
      return tmp;
    }
    inline CsgVector dividedBy(const MFloat a) const {
      CsgVector tmp;
      for(MInt i = 0; i < nDim; i++)
        tmp.xx[i] = this->xx[i] / a;
      return tmp;
    }
    inline MFloat dot(const CsgVector& a) const {
      MFloat dotp = F0;
      for(MInt i = 0; i < nDim; i++)
        dotp += this->xx[i] * a.xx[i];
      return dotp;
    }
    inline CsgVector lerp(const CsgVector& a, const MFloat t) const { return this->plus((a.minus(*this)).times(t)); }
    inline MFloat length() const { return sqrt(this->dot(*this)); }
    inline CsgVector unit() const { return this->dividedBy(this->length()); }
    template <class X = void, std::enable_if_t<nDim == 3, X*> = nullptr>
    inline CsgVector cross(const CsgVector& a) const {
      CsgVector tmp;
      tmp.xx[0] = this->xx[1] * a.xx[2] - this->xx[2] * a.xx[1];
      tmp.xx[1] = this->xx[2] * a.xx[0] - this->xx[0] * a.xx[2];
      tmp.xx[2] = this->xx[0] * a.xx[1] - this->xx[1] * a.xx[0];
      return tmp;
    }
    template <class X = void, std::enable_if_t<nDim == 2, X*> = nullptr>
    inline CsgVector cross(const CsgVector&) const {
      TERMM(1, "Get your ass out of this function. This is 3D code!");
    }
    inline void negate() {
      for(MInt i = 0; i < nDim; i++)
        this->xx[i] = -this->xx[i];
    }
  };


  class CsgVertex {
   public:
    CsgVector pos;
    MInt vertexId;
    MInt setIndex;
    CsgVertex(const CsgVector& _pos, MInt _vertexId, MInt _setIndex) : pos(_pos) {
      vertexId = _vertexId;
      setIndex = _setIndex;
    }
    CsgVertex clone() const {
      CsgVertex tmp(this->pos.clone(), vertexId, setIndex);
      return tmp;
    }
    CsgVertex interpolate(CsgVertex* other, MFloat t) const { return CsgVertex(this->pos.lerp(other->pos, t), -1, -1); }
  };

  class CsgPolygon;

  class CsgPlane {
   public:
    static constexpr MFloat eps = (nDim == 3) ? 10000 * MFloatEps : 1000 * MFloatEps;
    CsgVector normal;
    MFloat w{};

    CsgPlane() = default;
    CsgPlane(CsgVector* _normal, MFloat _w) {
      std::copy(_normal->xx, _normal->xx + nDim, normal.xx);
      w = _w;
    }
    explicit CsgPlane(CsgVector abc[nDim]) {
      IF_CONSTEXPR(nDim == 3) {
        CsgVector tmp = abc[1].minus(abc[0]).cross(abc[nDim - 1 /*2*/].minus(abc[0])).unit();
        normal = tmp;
        w = normal.dot(abc[0]);
      }
      else IF_CONSTEXPR(nDim == 2) {
        CsgVector tmp = abc[1].minus(abc[0]).unit();
        normal.xx[0] = -tmp.xx[1];
        normal.xx[1] = tmp.xx[0];
        w = normal.dot(abc[0]);
      }
    }
    ~CsgPlane() = default;
    CsgPlane clone() {
      CsgPlane tmp(&this->normal, this->w);
      return tmp;
    }
    void flip() {
      this->normal.negate();
      this->w = -this->w;
    }
    void splitPolygon(CsgPolygon* polygon,
                      std::vector<CsgPolygon>* coplanarFront,
                      std::vector<CsgPolygon>* coplanarBack,
                      std::vector<CsgPolygon>* front,
                      std::vector<CsgPolygon>* back) const;
    void insertCoplanarPolygon(CsgPolygon* polygon,
                               std::vector<CsgPolygon>* coplanarFront,
                               std::vector<CsgPolygon>* coplanarBack) const;
  };

  class CsgPolygon {
   public:
    std::vector<CsgVertex> vertices;
    // helpers: edgeType[0] corresponds to edge between vertices[0] and vertices[1],...
    CsgPlane plane;
    MInt setIndex{};
    MInt bodyId{};
    MInt faceId{};
    MInt faceType{}; // 0: cartesian face; 1: body face

    CsgPolygon(std::vector<CsgVertex> _vertices, MInt _setIndex, MInt _faceId, MInt _faceType, MInt _bodyId) {
      for(MInt i = 0; (unsigned)i < _vertices.size(); i++) {
        vertices.push_back(_vertices[i].clone());
      }
      setIndex = _setIndex;
      faceId = _faceId;
      faceType = _faceType;
      bodyId = _bodyId;
      CsgVector abc[nDim];
      for(MInt i = 0; i < nDim; i++)
        abc[i] = vertices[i].pos;
      CsgPlane tmp = CsgPlane(abc);
      plane = tmp;
    }

    CsgPolygon(std::vector<CsgVertex> _vertices, MInt _setIndex, MInt _faceId, MInt _faceType, MInt _bodyId,
               CsgPlane _plane) {
      for(MInt i = 0; (unsigned)i < _vertices.size(); i++) {
        vertices.push_back(_vertices[i].clone());
      }
      setIndex = _setIndex;
      faceId = _faceId;
      faceType = _faceType;
      bodyId = _bodyId;
      plane = _plane.clone();
    }

    ~CsgPolygon() = default;

    CsgPolygon clone() const {
      CsgPolygon tmp(vertices, setIndex, faceId, faceType, bodyId, plane);
      return tmp;
    }

    void flip() {
      this->plane.flip();
      std::vector<CsgVertex> tmp;
      for(auto rit = vertices.rbegin(); rit != vertices.rend(); ++rit) {
        tmp.push_back((*rit));
      }
      vertices.swap(tmp);
    }

   private:
    CsgPolygon() = default;
  };

  class CsgNode {
   public:
    CsgPlane plane;
    MBool planeValid;
    CsgNode* front;
    CsgNode* back;
    std::vector<CsgPolygon> polygons;

    explicit CsgNode(std::vector<CsgPolygon> _polygons) {
      planeValid = false;
      front = nullptr;
      back = nullptr;
      polygons.clear();
      this->build(std::move(_polygons));
    }

    CsgNode() {
      planeValid = false;
      front = nullptr;
      back = nullptr;
      polygons.clear();
    }

    ~CsgNode() {
      this->polygons.clear();
      delete this->front;
      delete this->back;
    }

    void invert() {
      for(MInt i = 0; (unsigned)i < this->polygons.size(); i++) {
        this->polygons[i].flip();
      }
      this->plane.flip();
      if(this->front) this->front->invert();
      if(this->back) this->back->invert();
      CsgNode* temp = this->front;
      this->front = this->back;
      this->back = temp;
    }

    std::vector<CsgPolygon> clipPolygons(std::vector<CsgPolygon> _polygons) {
      if(!planeValid) {
        return _polygons;
      }
      std::vector<CsgPolygon> _front;
      std::vector<CsgPolygon> _back;
      for(MInt i = 0; (unsigned)i < _polygons.size(); i++) {
        this->plane.splitPolygon(&_polygons[i], &_front, &_back, &_front, &_back);
      }
      if(this->front) _front = this->front->clipPolygons(_front);
      if(this->back)
        _back = this->back->clipPolygons(_back);
      else
        _back.clear();
      for(MInt i = 0; (unsigned)i < _back.size(); i++)
        _front.push_back(_back[i]);
      return _front;
    }

    void clipTo(CsgNode& bsp) {
      this->polygons = bsp.clipPolygons(this->polygons);
      if(this->front) this->front->clipTo(bsp);
      if(this->back) this->back->clipTo(bsp);
    }

    std::vector<CsgPolygon> allPolygons() {
      std::vector<CsgPolygon> _polygons(this->polygons);
      if(this->front) {
        std::vector<CsgPolygon> _polygons_front = this->front->allPolygons();
        for(MInt i = 0; (unsigned)i < _polygons_front.size(); i++)
          _polygons.push_back(_polygons_front[i]);
      }
      if(this->back) {
        std::vector<CsgPolygon> _polygons_back = this->back->allPolygons();
        for(MInt i = 0; (unsigned)i < _polygons_back.size(); i++)
          _polygons.push_back(_polygons_back[i]);
      }
      return _polygons;
    }

    void build(std::vector<CsgPolygon> _polygons);

    // This function originally only existed in 3D
    void plot(MInt tabCounter, const MChar* fileName) {
      std::ofstream ofl;
      if(tabCounter == 0) {
        ofl.open(fileName);
        ofl.close();
      }

      if(this->front) this->front->plot(tabCounter + 1, fileName);

      ofl.open(fileName, std::ofstream::app);

      if(ofl) {
        for(MInt count = 0; count < tabCounter; count++)
          ofl << "      ";
        ofl << "node plane normal, w:";
        for(MInt dim = 0; dim < nDim; ++dim)
          ofl << " " << this->plane.normal.xx[dim];
        ofl << ", " << this->plane.w << std::endl;

        for(MInt p = 0; (unsigned)p < this->polygons.size(); p++) {
          for(MInt count = 0; count < tabCounter; count++)
            ofl << "      ";
          ofl << "polygon " << p << " plane normal, w:";
          for(MInt dim = 0; dim < nDim; ++dim)
            ofl << " " << this->polygons[p].plane.normal.xx[dim];
          ofl << ", " << this->polygons[p].plane.w << std::endl;

          for(MInt count = 0; count < tabCounter; count++)
            ofl << "      ";
          ofl << "polygon setIndex: " << this->polygons[p].setIndex << ", bodyId: " << this->polygons[p].bodyId
              << ", faceId: " << this->polygons[p].faceId << ", faceType: " << this->polygons[p].faceType << std::endl;
        }
        ofl << std::endl << std::endl;
        ofl.close();
      }

      if(this->back) this->back->plot(tabCounter + 1, fileName);
    }
  };

  class Csg {
   public:
    std::vector<CsgPolygon> polygons;

    Csg() { polygons.clear(); }

    explicit Csg(std::vector<CsgPolygon> _polygons) {
      for(MInt i = 0; (unsigned)i < _polygons.size(); i++) {
        polygons.push_back(_polygons[i].clone());
      }
    }

    ~Csg() = default;

    std::vector<CsgPolygon> intersect(const Csg& csg) const {
      CsgNode a(this->polygons);
      CsgNode b(csg.polygons);

      b.invert();
      b.clipTo(a);
      b.invert();
      a.clipTo(b);
      b.clipTo(a);
      a.build(b.allPolygons());

      return a.allPolygons();
    }

    void inverse() {
      for(MInt i = 0; (unsigned)i < this->polygons.size(); i++) {
        this->polygons[i].flip();
      }
    }
  };


  // Data members used to replace local static variables
  MBool m_static_initSolutionStep_firstRun = true;
  MBool m_static_initSolutionStep_frstrn = true;
  MInt m_static_createCutFaceMb_MGC_bodyFaceJoinMode = 1;
  MFloat m_static_createCutFaceMb_MGC_maxA = 0.1;
  MBool m_static_createCutFaceMb_MGC_first = true;

 private:
  struct CellDataDlb;
  static const std::array<MInt, CellDataDlb::count> s_cellDataTypeDlb;
};


// Struct to differentiate cell (element) data
template <MInt nDim_, class SysEqn>
struct FvMbCartesianSolverXD<nDim_, SysEqn>::CellDataDlb {
  static constexpr const MInt count = 8;

  static constexpr const MInt ELEM_AVARIABLE = 0;
  static constexpr const MInt ELEM_CELLVOLUMES = 1;
  static constexpr const MInt ELEM_RIGHTHANDSIDE = 2;
  static constexpr const MInt ELEM_CELLVOLUMESDT1 = 3;
  static constexpr const MInt ELEM_SWEPTVOLUMES = 4;
  static constexpr const MInt ELEM_PROPERTIES = 5;
  static constexpr const MInt ELEM_AZIMUTHAL_LONG = 6;
  static constexpr const MInt ELEM_AZIMUTHAL_FLOAT = 7;
  static constexpr const MInt ELEM_LESAVERAGE = 8;
};


// Data types of cell data
template <MInt nDim, class SysEqn>
const std::array<MInt, FvMbCartesianSolverXD<nDim, SysEqn>::CellDataDlb::count>
    FvMbCartesianSolverXD<nDim, SysEqn>::s_cellDataTypeDlb = {
        {MFLOAT, MFLOAT, MFLOAT, MFLOAT, MFLOAT, MLONG, MLONG, MFLOAT}};

/// \brief Return data size to be communicated during DLB for a grid cell and given data id
template <MInt nDim_, class SysEqn>
MInt FvMbCartesianSolverXD<nDim_, SysEqn>::cellDataSizeDlb(const MInt dataId, const MInt gridCellId) {
  // Inactive ranks do not have any data to communicate
  if(!isActive()) {
    return 0;
  }

  // Convert to solver cell id and check
  const MInt cellId = grid().tree().grid2solver(gridCellId);
  if(cellId < 0 || cellId >= noInternalCells()) {
    return 0;
  }

  MInt dataSize = 1;

  switch(dataId) {
    case CellDataDlb::ELEM_AVARIABLE: {
      dataSize = m_noCVars;
      break;
    }
    case CellDataDlb::ELEM_CELLVOLUMES: {
      dataSize = 1;
      break;
    }
    case CellDataDlb::ELEM_RIGHTHANDSIDE: {
      dataSize = m_noFVars;
      break;
    }
    case CellDataDlb::ELEM_CELLVOLUMESDT1: {
      dataSize = 1;
      break;
    }
    case CellDataDlb::ELEM_SWEPTVOLUMES: {
      dataSize = 1;
      break;
    }
    case CellDataDlb::ELEM_PROPERTIES: {
      dataSize = 1;
      break;
    }
    case CellDataDlb::ELEM_AZIMUTHAL_LONG: {
      if(grid().azimuthalPeriodicity()) {
        dataSize = m_azimuthalNearBoundaryBackupMaxCount * m_noLongDataBalance;
      } else {
        dataSize = 0;
      }
      break;
    }
    case CellDataDlb::ELEM_AZIMUTHAL_FLOAT: {
      if(grid().azimuthalPeriodicity()) {
        dataSize = m_azimuthalNearBoundaryBackupMaxCount * (nDim_ + m_noFloatDataBalance);
      } else {
        dataSize = 0;
      }
      break;
    }
    default: {
      TERMM(1, "Unknown data id.");
      break;
    }
  }

  return dataSize;
}


/// \brief Store the solver data for a given data id ordered in the given buffer for DLB
template <MInt nDim_, class SysEqn_>
void FvMbCartesianSolverXD<nDim_, SysEqn_>::getCellDataDlb(const MInt dataId, const MInt oldNoCells,
                                                           const MInt* const bufferIdToCellId, MFloat* const data) {
  TRACE();

  MInt localBufferId = 0;
  for(MInt i = 0; i < oldNoCells; i++) {
    const MInt gridCellId = bufferIdToCellId[i];

    if(gridCellId < 0) continue;

    const MInt cellId = grid().tree().grid2solver(gridCellId);
    if(cellId < 0 || cellId >= noInternalCells()) {
      continue;
    }

    switch(dataId) {
      case CellDataDlb::ELEM_AVARIABLE: {
        std::copy_n(&a_variable(cellId, 0), m_noCVars, &data[localBufferId * m_noCVars]);
        break;
      }
      case CellDataDlb::ELEM_CELLVOLUMES: {
        std::copy_n(&a_cellVolume(cellId), 1, &data[localBufferId]);
        break;
      }
      case CellDataDlb::ELEM_RIGHTHANDSIDE: {
        std::copy_n(&a_rightHandSide(cellId, 0), m_noFVars, &data[localBufferId * m_noFVars]);
        break;
      }
      case CellDataDlb::ELEM_CELLVOLUMESDT1: {
        std::copy_n(&m_cellVolumesDt1[cellId], 1, &data[localBufferId]);
        break;
      }
      case CellDataDlb::ELEM_SWEPTVOLUMES: {
        MInt bndryId_ = a_bndryId(cellId);
        if(bndryId_ > -1) {
          data[localBufferId] = m_sweptVolume[bndryId_];
        } else {
          data[localBufferId] = -1;
        }
        break;
      }
      case CellDataDlb::ELEM_AZIMUTHAL_FLOAT: {
        if(grid().azimuthalPeriodicity()) {
          MInt noFloatData = m_azimuthalNearBoundaryBackupMaxCount * (nDim_ + m_noFloatDataBalance);
          std::copy_n(&m_azimuthalNearBoundaryBackupBalFloat[cellId * noFloatData], noFloatData,
                      &data[localBufferId * noFloatData]);
        }
        break;
      }
      default:
        TERMM(1, "Unknown data id.");
        break;
    }
    localBufferId++;
  }
}


/// \brief Store the solver data for a given data id ordered in the given buffer for DLB
template <MInt nDim_, class SysEqn_>
void FvMbCartesianSolverXD<nDim_, SysEqn_>::getCellDataDlb(const MInt dataId, const MInt oldNoCells,
                                                           const MInt* const bufferIdToCellId, MLong* const data) {
  TRACE();

  MInt localBufferId = 0;
  for(MInt i = 0; i < oldNoCells; i++) {
    const MInt gridCellId = bufferIdToCellId[i];

    if(gridCellId < 0) continue;

    const MInt cellId = grid().tree().grid2solver(gridCellId);
    if(cellId < 0 || cellId >= noInternalCells()) {
      continue;
    }

    switch(dataId) {
      case CellDataDlb::ELEM_PROPERTIES: {
        data[localBufferId] = (MLong)a_properties(cellId).to_ulong();
        break;
      }
      case CellDataDlb::ELEM_AZIMUTHAL_LONG: {
        if(grid().azimuthalPeriodicity()) {
          MInt noLongData = m_azimuthalNearBoundaryBackupMaxCount * m_noLongDataBalance;
          std::copy_n(&m_azimuthalNearBoundaryBackupBalLong[cellId * noLongData], noLongData,
                      &data[localBufferId * noLongData]);
        }
        break;
      }
      default:
        TERMM(1, "Unknown data id.");
        break;
    }
    localBufferId++;
  }
}


/// \brief Set the solver cell data after DLB
template <MInt nDim_, class SysEqn_>
void FvMbCartesianSolverXD<nDim_, SysEqn_>::setCellDataDlb(const MInt dataId, const MFloat* const data) {
  TRACE();

  // Nothing to do if solver is not active
  if(!isActive()) {
    return;
  }

  // Set the variables if this is the correct reinitialization stage
  if(m_loadBalancingReinitStage == 0) {
    switch(dataId) {
      case CellDataDlb::ELEM_AVARIABLE: {
        std::copy_n(data, noInternalCells() * m_noCVars, &a_variable(0, 0));
        break;
      }
      case CellDataDlb::ELEM_CELLVOLUMES: {
        std::copy_n(data, noInternalCells(), &a_cellVolume(0));
        break;
      }
      case CellDataDlb::ELEM_RIGHTHANDSIDE: {
        std::copy_n(data, noInternalCells() * m_noFVars, &a_rightHandSide(0, 0));
        break;
      }
      case CellDataDlb::ELEM_CELLVOLUMESDT1: {
        std::copy_n(data, noInternalCells(), &m_cellVolumesDt1[0]);
        break;
      }
      case CellDataDlb::ELEM_SWEPTVOLUMES: {
        std::copy_n(data, noInternalCells(), &m_sweptVolumeBal[0]);
        break;
      }
      case CellDataDlb::ELEM_AZIMUTHAL_FLOAT: {
        if(grid().azimuthalPeriodicity()) {
          MInt noFloatData = m_azimuthalNearBoundaryBackupMaxCount * (nDim_ + m_noFloatDataBalance);
          std::copy_n(data, noInternalCells() * noFloatData, &m_azimuthalNearBoundaryBackupBalFloat[0]);
        }
        break;
      }
      default:
        TERMM(1, "Unknown data id.");
    }
  }
}


/// \brief Set the solver cell data after DLB
template <MInt nDim_, class SysEqn_>
void FvMbCartesianSolverXD<nDim_, SysEqn_>::setCellDataDlb(const MInt dataId, const MLong* const data) {
  TRACE();

  // Nothing to do if solver is not active
  if(!isActive()) {
    return;
  }

  // Set the variables if this is the correct reinitialization stage
  if(m_loadBalancingReinitStage == 0) {
    switch(dataId) {
      case CellDataDlb::ELEM_PROPERTIES: {
        for(MInt i = 0; i < noInternalCells(); i++) {
          a_properties(i) = maia::fv::cell::BitsetType((MUlong)data[i]);
        }
        break;
      }
      case CellDataDlb::ELEM_AZIMUTHAL_LONG: {
        if(grid().azimuthalPeriodicity()) {
          MInt noLongData = m_azimuthalNearBoundaryBackupMaxCount * m_noLongDataBalance;
          std::copy_n(data, noInternalCells() * noLongData, &m_azimuthalNearBoundaryBackupBalLong[0]);
        }
        break;
      }
      default:
        TERMM(1, "Unknown data id.");
    }
  }
}

#endif
