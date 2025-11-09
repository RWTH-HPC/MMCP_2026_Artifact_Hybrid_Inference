// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef LBBNDCND_H
#define LBBNDCND_H

#include <fstream>
#include <unordered_map>
#include <vector>
#include "GEOM/geometryelement.h"
#include "GRID/cartesiangridcellproperties.h"
#include "IO/parallelio.h"
#include "globals.h"
#include "lbgridboundarycell.h"
#include "lbmbcellcollector.h"

template <MInt nDim>
class LbSolver;
template <class T>
class Collector;
template <MInt nDim>
class CartesianGrid;

struct LbBndCndData {
  MInt noBndCells = 0;                        ///< this is number of non-halo (!) boundary cells -> stored for restart
  MInt noBndCellsWithHalos = 0;               ///< number of boundary cells represented in data
  MInt noVars = 0;                            ///< number of variables that are stored for BC
  ParallelIo::size_type globalOffset = 0;     ///< global offset for writing restart file
  MFloat* data = nullptr;                     ///< pointer the relevant data
};

namespace maia::lb {
struct CalcWallForceContainer {
  MBool isRoot = false;
  MInt noComm = 0;
  MString fileName = "force.log";
  MPI_Comm comm;
};
} // namespace maia::lb

template <MInt nDim>
class LbBndCnd {
 public:
  template <MInt nDim_>
  friend class LbSolver;
  template <MInt nDim_, MInt nDist_, class SysEqn_>
  friend class CouplingLB;

  // Types
  using Cell = GridCell;

 protected:
  using MbCellCollector = maia::lb::collector::LbMbCellCollector<nDim>;
  MbCellCollector m_boundaryCellsMb{};

  //! pointer to a member function data type
  LbSolver<nDim>* m_solver;
  typedef void (LbBndCnd::*BndCndHandler)(MInt set);
  std::vector<BndCndHandler> bndCndHandlerVariables;
  std::vector<BndCndHandler> bndCndHandlerRHS;
  std::vector<BndCndHandler> bndCndInitiator;
  typedef void (LbBndCnd::*BndCndHandler_Mb)(MInt set);
  std::vector<BndCndHandler> bndCndHandlerVariables_MB;
  std::vector<BndCndHandler_Mb> bndCndHandlerRHS_MB;
  MInt m_noInternalCells;
  MInt m_initialNoCells[3]{};

  std::map<MInt, MInt> m_boundaryCellMappingMb;

  template <typename K, typename V>
  V GetWithDef(const std::map<K, V>& m, const K& key, const V& defval) {
    typename std::map<K, V>::const_iterator it = m.find(key);
    if(it == m.end()) {
      return defval;
    } else {
      return it->second;
    }
  }

  MInt a_boundaryCellMb(const MInt cellId) { return GetWithDef(m_boundaryCellMappingMb, cellId, -1); }

  // Members needed for segment velocity and normal vector calculation
  MInt m_noInOutSegments{};
  MInt m_noPeriodicSegments;
  MInt m_noSegments{};
  MFloat** m_initialVelocityVecs = nullptr;
  MFloat** m_bndNormals = nullptr;
  MInt* m_inOutSegmentsIds = nullptr;
  MInt* m_periodicSegmentsIds;
  MInt* m_bndNormalDirection = nullptr;
  MInt m_lbNoMovingWalls;
  MInt* m_segIdMovingWalls{};
  MFloat* m_lbWallVelocity{};
  MBool m_calcWallForces;
  MInt m_calcWallForcesInterval;
  MInt m_lbNoHeatedWalls;
  MInt* m_segIdHeatedWalls{};
  MFloat* m_lbWallTemperature{};

  // external forcing terms
  MFloat* m_Fext;
  MInt m_maxNoG0CellsLb;
  //  MInt m_currentNoG0Cells;
  MInt m_maxNoDistributionsInDim;

  // Read from property file
  MString m_initVelocityMethod;
  MString m_bndNormalMethod;
  MInt m_fastParallelGeomNormals{};
  MString m_interpolationDistMethod;
  MBool m_outputWallDistanceField = false;
  MString m_multiBCTreatment;

  // Data for coordinate transformation
  MFloat* m_phi = nullptr;
  MFloat* m_theta = nullptr;
  // MFloat m_cosTheta, m_sinTheta, m_cosPhi, m_sinPhi;
  // MFloat m_cosThetaSq, m_sec2Theta, m_sec2Phi, m_cos2Theta, m_sinThetaSq, m_tan2Theta, m_cotTheta;

  // Data for extrapolation in arbitrary directions
  // MInt** m_axesDirs;
  MInt** m_exDirs = nullptr;
  MFloat** m_exWeights = nullptr;

  //! index of the array where the boundary conditions are stored in
  MInt m_solverId;
  std::vector<LbGridBoundaryCell<nDim>> m_bndCells;

  MString m_gridCutTest;

  MPrimitiveVariables<nDim>* PV;
  // Defines the discretization Model
  MInt m_noDistributions;
  MInt m_methodId;

  //! physical non-dimensional reference length
  MFloat m_omega{};
  MFloat m_Ma;
  MFloat m_Re{};
  MFloat m_nu{};
  MFloat m_gradient{};
  MFloat m_referenceLength;
  MFloat m_domainLength;
  MInt m_densityFluctuations;

  MFloat m_rho1{}; // value of rho at inflow boundary
  MFloat m_rho2{}; // value of rho at outflow boundary

  /*   MInt m_tanhInit; */
  /*   MFloat m_initRe; */
  /*   MInt m_initTime; */
  /*   MInt m_initStartTime; */

  MInt m_lbControlInflow{};
  std::vector<MInt> m_bndCndIds;            // Holds the different BC ids
  std::vector<MInt> m_bndCndOffsets;        // stores the starting positions
  std::vector<MInt> m_bndCndSegIds;         // Holds the different segment ids
  std::vector<MInt> m_mapSegIdsInOutCnd;    // maps segments to in/out segments
  std::vector<MInt> m_mapBndCndSegId2Index; // maps global segmentId to local index (-1 if not available)
  std::vector<MInt> m_mapIndex2BndCndSegId; // maps local index to global segmentId (-1 if not available)
  std::vector<MInt> m_noBndCellsPerSegment; // number of cells that carry a certain boundary condition

  MFloat** m_blasius{};
  MFloat m_blasius_delta{};

  // pulsatile
  MFloat m_pulsatileFrequency;

  // latentHeat
  MBool m_latentHeat = false;
  MBool m_calcSublayerDist = false;

  // BC data
  std::unordered_map<MInt, LbBndCndData> m_bndCndData{};     ///< Stores BC specific data mapped by boundary index
  std::vector<MInt> m_segIdUseBcData{};                      ///< hold number of domains using BC data for this segment

 public:
  LbBndCnd(LbSolver<nDim>* solver);

  virtual ~LbBndCnd();
  virtual void setBndCndHandler();
  virtual void updateVariables();
  virtual void updateRHS();
  void initializeBcData();
  virtual void createMBComm();
  virtual void postCouple();
  virtual void initializeBndMovingBoundaries();

  // void calculateBodyForces(MInt bodyId, MFloat* bodyForces);

 private:
  virtual void createBoundaryCells();
  // ANDI
  virtual void initMembers();
  virtual void calculateVectors();
  MBool (LbBndCnd::*retrieveNormal)(GeometryElement<nDim> ge, MFloat* normal);
  inline MBool calculateNormalFromTriangle(GeometryElement<nDim> ge, MFloat* normal);
  inline MBool getNormalFromSTL(GeometryElement<nDim> ge, MFloat* normal);
  void calculateBndNormals();
  virtual void calculateAveragedNormalAndCenter(MInt segmentId, MFloat* const normal, MFloat* const center);
  virtual void fastParallelGeomNormals3D(std::vector<std::pair<MInt, MInt>> own_segments);
  virtual void normalParallelGeomNormals3D(std::vector<std::pair<MInt, MInt>> own_segments);
  virtual void updateBndNormals(MInt segId, MBool inside, MFloat* avg_normal);
  virtual void checkBndNormalDirection(){};
  virtual void printBndVelInfo();
  virtual void processAngles();

  void bcDataAllocate(MInt index, MInt noVars);
  void bcDataWriteRestartHeader(ParallelIo& parallelio);
  void bcDataWriteRestartData(ParallelIo& parallelIo);
  void bcDataReadRestartData(ParallelIo& parallelIo);
  void applyDirectionChangeInflow(MInt index);
  void applyDirectionChangeOutflow(MInt index);

  // TIMERS
  MInt m_t_BCAll;

  MFloat* sendBuffersMB = nullptr;
  MFloat* receiveBuffersMB = nullptr;
  MInt* haloInformation = nullptr;
  MInt* windowInformation = nullptr;
  MPI_Request* mpi_request = nullptr;

 protected:
  // ANDI
  virtual MInt findBndCnd(MInt index);

  virtual void solveBlasiusZ(MInt index);

  void sortBoundaryCells();

  virtual void calculateForces(MInt /*BndCndId*/){};
  virtual void calculateWallDistances(){};
  virtual void addWallDistanceFieldToOutputFile(ParallelIo& parallelio, const MBool writeHeader,
                                                const MBool writeData) = 0;
  virtual void createChannelBoundaries(){};

  virtual void bc0(MInt /*index*/){};

  virtual void bc10000(MInt /*index*/){};
  virtual void bc10001(MInt /*index*/){};
  virtual void bc10002(MInt /*index*/){};
  virtual void bc10004(MInt /*index*/){};
  virtual void bc10022(MInt /*index*/){};
  virtual void bc10010(MInt /*index*/){};
  virtual void bc10020(MInt /*index*/){};
  virtual void bc10040(MInt /*index*/){};
  virtual void bc10041(MInt /*index*/){};
  virtual void bc10042(MInt /*index*/){};
  virtual void bc10043(MInt /*index*/){};
  virtual void bc10044(MInt /*index*/){};
  virtual void bc10045(MInt /*index*/){};
  virtual void bc10046(MInt /*index*/){};

  virtual void bc10050(MInt /*index*/){};
  virtual void bc10060(MInt /*index*/){};
  virtual void bc10061(MInt /*index*/){};

  virtual void bc10090(MInt /*index*/){};

  virtual void bc10070(MInt /*index*/){};

  virtual void bc10080(MInt /*index*/){};

  virtual void bc20000(MInt /*index*/){};
  virtual void bc20001(MInt /*index*/){};
  virtual void bc20002(MInt /*index*/){};
  virtual void bc20003(MInt /*index*/){};
  virtual void bc20004(MInt /*index*/){};
  virtual void bc20005(MInt /*index*/){};
  virtual void bc20006(MInt /*index*/){};
  virtual void bc20010(MInt /*index*/){};
  virtual void bc20020(MInt /*index*/){};
  virtual void bc20220(MInt /*index*/){};

  virtual void bc20226(MInt /*index*/){};
  virtual void bc20227(MInt /*index*/){};
  virtual void bc20228(MInt /*index*/){};
  virtual void bc20501(MInt /*index*/){};
  virtual void bc20501_init(MInt /*index*/){};

  virtual void bcIBBNeumannInit(MInt /*index*/){};

  virtual void bc20022(MInt /*index*/){};
  virtual void bc20023(MInt /*index*/){};
  virtual void bc20024(MInt /*index*/){};
  virtual void bc20025(MInt /*index*/){};
  virtual void bc20026(MInt /*index*/){};
  virtual void bc20027(MInt /*index*/){};
  virtual void bc20030(MInt /*index*/){};
  virtual void bc20230(MInt /*index*/){};

  virtual void bc20050(MInt /*index*/){};
  virtual void bc20051(MInt /*index*/){};
  virtual void bc20052(MInt /*index*/){};
  virtual void bc20053(MInt /*index*/){};
  virtual void bc20054(MInt /*index*/){};
  virtual void bc20055(MInt /*index*/){};

  virtual void bc30000(MInt /*index*/){};

  virtual void bc30010(MInt /*index*/){};
  virtual void bc30011(MInt /*index*/){};
  virtual void bc30012(MInt /*index*/){};
  virtual void bc30013(MInt /*index*/){};
  virtual void bc30014(MInt /*index*/){};
  virtual void bc30015(MInt /*index*/){};

  virtual void bc30020(MInt /*index*/){};
  virtual void bc30021(MInt /*index*/){};
  virtual void bc30022(MInt /*index*/){};
  virtual void bc30023(MInt /*index*/){};
  virtual void bc30024(MInt /*index*/){};
  virtual void bc30025(MInt /*index*/){};

  virtual void bc30030(MInt /*index*/){};
  virtual void bc30031(MInt /*index*/){};
  virtual void bc30032(MInt /*index*/){};
  virtual void bc30033(MInt /*index*/){};
  virtual void bc30034(MInt /*index*/){};
  virtual void bc30035(MInt /*index*/){};

  virtual void bc30040(MInt /*index*/){};
  virtual void bc30041(MInt /*index*/){};
  virtual void bc30042(MInt /*index*/){};
  virtual void bc30043(MInt /*index*/){};
  virtual void bc30044(MInt /*index*/){};
  virtual void bc30045(MInt /*index*/){};

  virtual void bc30050(MInt /*index*/){};

  virtual void bc40000(MInt /*index*/){};
  virtual void bc40001(MInt /*index*/){};
  virtual void bc40010(MInt /*index*/){};
  virtual void bc40020(MInt /*index*/){};
  virtual void bc40030(MInt /*index*/){};

  virtual void bc40040(MInt /*index*/){};
  virtual void bc40041(MInt /*index*/){};
  virtual void bc40042(MInt /*index*/){};
  virtual void bc40043(MInt /*index*/){};
  virtual void bc40044(MInt /*index*/){};
  virtual void bc40045(MInt /*index*/){};
  virtual void bc40046(MInt /*index*/){};

  virtual void bc40060(MInt /*index*/){};

  virtual void bc40070(MInt /*index*/){};
  virtual void bc40071(MInt /*index*/){};
  virtual void bc40072(MInt /*index*/){};
  virtual void bc40073(MInt /*index*/){};
  virtual void bc40080(MInt /*index*/){};
  virtual void bc40081(MInt /*index*/){};
  virtual void bc40082(MInt /*index*/){};

  virtual void bc40072_40082_init(MInt /*index*/){};

  virtual void bc40100(MInt /*index*/){};
  virtual void bc40110(MInt /*index*/){};

  virtual void bc40120(MInt /*index*/){};

  virtual void bc40130(MInt /*index*/){};

  virtual void bc10111(MInt /*index*/){};

  virtual void bc66666(MInt /*set*/){};
  virtual void bc66667(MInt /*set*/){};
  virtual void bc66668(MInt /*set*/){};

  // for neighbor communicators
  virtual MInt checkForCommBC();
  virtual MBool checkForCommForce();
  virtual void setBCNeighborCommunicator();
  virtual void prepareBC(MInt index, MInt BCCounter, MInt segId);
  virtual void prepareBC4073(MInt BCCounter, MInt segId);
  virtual void setBCWallNeighborCommunicator();
  virtual void prepareBC2000();

  MFloat m_deltaRho{};
  MFloat m_ReLast{};
  MFloat m_rhoLast{};
  MFloat m_lRho{};
  MFloat m_maxDeltaRho{};
  MInt m_currentTimeStep{};

  MPI_Group* tmp_group{};
  MPI_Group* BCGroup{};
  MPI_Comm* m_BCComm{};
  MInt* m_mapBndCndIdSegId{};
  MInt* m_totalNoBcCells{};

  MString* m_BCOutputFileName{};
  MInt* m_allBoundaryIds{};
  MInt m_noAllBoundaryIds{};
  MInt** m_allDomainsHaveBC{};
  MInt m_numberOfCommBCs = 0;
  MInt* m_noBCNeighbors{};
  std::vector<std::vector<MInt>> m_BCneighbors;
  MBool m_calcBcResidual;
  std::ofstream* m_BCResidualStream{};
  MFloat* m_localReCutPoint{};
  MFloat* m_localReCutNormal{};
  MFloat m_localReCutDistance{};
  MFloat m_localReCutRe{};
  MFloat m_localReCutDiameter{};
  MInt* m_allDomainsCalcForceMB{};
  MBool m_hasCommForce;
  MPI_Comm m_BCWallMBComm;
  MInt m_noBCWallMBNeighbors;

  std::unordered_map<MInt, maia::lb::CalcWallForceContainer> m_mapWallForceContainer;
  MString m_forceFile;
  std::vector<MInt> m_BCWallMBNeighbors;
  std::vector<MInt> m_localReCutCells;
  MBool m_hasLocalReCut{};
  MInt m_totalNoDomainsReCut{};
  MInt m_localReCutInterval{};
  MFloat m_localReCutAdpPerc{};
  MInt m_localReCutReportInterval{};
  MInt m_firstBCinComm{};
  MInt m_noReactivatedCells;
};

#endif
