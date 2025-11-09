// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef FVSOLVERSTRUCTWINDOWSINFO
#define FVSOLVERSTRUCTWINDOWSINFO

#include <cstring>
#include <map>
#include <set>
#include <vector>
#include "GRID/structuredgrid.h"
#include "GRID/structuredpartition.h"
#include "INCLUDE/maiatypes.h"
#include "fvstructuredcomm.h"
#include "fvstructuredwindowmapping.h"

template <MInt nDim>
struct pointType {
  MInt BC[nDim]{};
  MInt blockId;
  MInt pos[nDim];
  MInt found;
};

class connectionNode {
 public:
  MInt nDim;
  MInt blockId1;
  MInt pos1[3];
  MInt blockId2;
  MInt pos2[3];
  MInt BC;
  MInt Nstar;
  connectionNode(MInt _nDim) {
    BC = -1;
    blockId1 = -1;
    blockId2 = -1;
    nDim = _nDim;
    for(MInt dim = 0; dim < nDim; dim++) {
      pos1[dim] = -1;
      pos2[dim] = -1;
    }
    Nstar = -1;
  };

  connectionNode(MInt a, MInt b1, MInt* p1, MInt b2, MInt* p2, MInt _nDim) {
    BC = a;
    blockId1 = b1;
    blockId2 = b2;
    nDim = _nDim;
    for(MInt dim = 0; dim < nDim; dim++) {
      pos1[dim] = p1[dim];
      pos2[dim] = p2[dim];
    }
    Nstar = -1;
  };

  connectionNode(MInt a, MInt b1, MInt* p1, MInt b2, MInt* p2, MBool enableSwap, MInt _nDim) {
    MBool swap = false;
    nDim = _nDim;

    if(enableSwap) {
      // ensure a fixed order, i.e., lower blockId first or lower pos first if the blockIds match
      if(b2 < b1) swap = true;
      if(b1 == b2) {
        for(MInt countDim = 0; countDim < nDim; ++countDim) {
          if(p1[countDim] < p2[countDim]) break;
          if(p1[countDim] > p2[countDim]) {
            swap = true;
            break;
          }
        }
      }
    }

    if(swap) {
      BC = a;
      blockId1 = b2;
      blockId2 = b1;
      for(MInt dim = 0; dim < nDim; dim++) {
        pos1[dim] = p2[dim];
        pos2[dim] = p1[dim];
      }
    } else {
      BC = a;
      blockId1 = b1;
      blockId2 = b2;
      for(MInt dim = 0; dim < nDim; dim++) {
        pos1[dim] = p1[dim];
        pos2[dim] = p2[dim];
      }
    }
    Nstar = -1;
  };

  void print() {
    std::cout << "============================" << std::endl;
    std::cout << BC << std::endl;
    std::cout << blockId1;
    for(MInt dim = 0; dim < nDim; dim++) {
      std::cout << " " << pos1[dim];
    }
    std::cout << std::endl;
    std::cout << blockId2;
    for(MInt dim = 0; dim < nDim; dim++) {
      std::cout << " " << pos2[dim];
    }
    std::cout << std::endl;
    std::cout << "Nstar: " << Nstar << std::endl;
    std::cout << "============================" << std::endl;
  };
  MBool operator<(const connectionNode& entry2) const {
    if(BC < entry2.BC) {
      return true;
    }
    if(BC > entry2.BC) {
      return false;
    }
    if(blockId1 < entry2.blockId1) {
      return true;
    }
    if(blockId1 > entry2.blockId1) {
      return false;
    }
    if(blockId2 < entry2.blockId2) {
      return true;
    }
    if(blockId2 > entry2.blockId2) {
      return false;
    }
    for(MInt i = 0; i < nDim; ++i) {
      if(pos1[i] < entry2.pos1[i]) {
        return true;
      }
      if(pos1[i] > entry2.pos1[i]) {
        return false;
      }
      if(pos2[i] < entry2.pos2[i]) {
        return true;
      }
      if(pos2[i] > entry2.pos2[i]) {
        return false;
      }
    }
    return false;
  };
  MBool operator==(const connectionNode& entry2) const {
    if(BC != entry2.BC) {
      return false;
    }
    if(blockId1 != entry2.blockId1) {
      return false;
    }
    if(blockId2 != entry2.blockId2) {
      return false;
    }
    for(MInt i = 0; i < nDim; ++i) {
      if(pos1[i] != entry2.pos1[i]) {
        return false;
      }
      if(pos2[i] != entry2.pos2[i]) {
        return false;
      }
    }
    return true;
  };
};


template <MInt nDim>
class windowInformation {
 public:
  windowInformation(){};
  ~windowInformation(){};
  MInt fixDir = -1;
  MInt inDir1 = -1;
  MInt inDir2 = -1;
  MInt startindex[nDim]{};
  MInt endindex[nDim]{};
  MInt domainId = -1;
  MInt BC = -1;
  MInt blockId = -1;
  MInt windowId = -1;
};


class ParallelIoHdf5;

template <MInt nDim>
class FvStructuredSolverWindowInfo {
  template <MInt nDim_>
  friend class FvStructuredSolver;
  friend class FvStructuredSolver3D;
  friend class FvStructuredSolver2D;
  template <MInt nDim_>
  friend class StructuredBndryCnd;
  template <MInt nDim_>
  friend class StructuredDecomposition;

 public:
  FvStructuredSolverWindowInfo(StructuredGrid<nDim>*, MPI_Comm, const MInt, const MInt, const MInt);
  ~FvStructuredSolverWindowInfo();
  void readWindowInfo();
  void mapCreate(MInt Id1, MInt* start1, MInt* end1, MInt* step1, MInt Id2, MInt* start2, MInt* end2, MInt* step2,
                 MInt* order, MInt BC, std::unique_ptr<StructuredWindowMap<nDim>>& output);

  void setSpongeInformation(MInt noSpongeInfo, MFloat* beta, MFloat* sigma, MFloat* thickness, MInt* bcInfo,
                            MInt informationType);
  void setWallInformation();
  void setLocalWallInformation();
  void setZonalBCInformation();

  MBool mapCheck(std::unique_ptr<StructuredWindowMap<nDim>>& input);
  MBool mapCheck0d(std::unique_ptr<StructuredWindowMap<nDim>>& input);
  MBool mapCheck1d(std::unique_ptr<StructuredWindowMap<nDim>>& input);
  MBool mapCheck3d(std::unique_ptr<StructuredWindowMap<nDim>>& input);
  MBool mapCheck2d(std::unique_ptr<StructuredWindowMap<nDim>>& input);
  MBool mapCheckWave(std::unique_ptr<StructuredWindowMap<nDim>>& input);
  void mapPrint(const std::unique_ptr<StructuredWindowMap<nDim>>& input);
  void mapPrintSimple(std::unique_ptr<StructuredWindowMap<nDim>>& input);
  void mapZero(std::unique_ptr<StructuredWindowMap<nDim>>& output);
  MInt mapCompare(std::unique_ptr<StructuredWindowMap<nDim>>& map1, std::unique_ptr<StructuredWindowMap<nDim>>& map2);
  MInt mapCompare11(const std::unique_ptr<StructuredWindowMap<nDim>>& map1,
                    const std::unique_ptr<StructuredWindowMap<nDim>>& map2);
  void mapInvert(std::unique_ptr<StructuredWindowMap<nDim>>& input, std::unique_ptr<StructuredWindowMap<nDim>>& output);
  void mapInvert1(std::unique_ptr<StructuredWindowMap<nDim>>& output);
  void mapNormalize1(std::unique_ptr<StructuredWindowMap<nDim>>& input,
                     std::unique_ptr<StructuredWindowMap<nDim>>& output);
  void mapNormalize2(std::unique_ptr<StructuredWindowMap<nDim>>& input,
                     std::unique_ptr<StructuredWindowMap<nDim>>& output);
  void mapNormalize3(std::unique_ptr<StructuredWindowMap<nDim>>& output);
  void mapCombine11(std::unique_ptr<StructuredWindowMap<nDim>>& input1,
                    std::unique_ptr<StructuredWindowMap<nDim>>& input2,
                    std::unique_ptr<StructuredWindowMap<nDim>>& output);
  void mapCombine12(std::unique_ptr<StructuredWindowMap<nDim>>& input1,
                    std::unique_ptr<StructuredWindowMap<nDim>>& input2,
                    std::unique_ptr<StructuredWindowMap<nDim>>& output);
  void mapCombine21(std::unique_ptr<StructuredWindowMap<nDim>>& input1,
                    std::unique_ptr<StructuredWindowMap<nDim>>& input2,
                    std::unique_ptr<StructuredWindowMap<nDim>>& output);
  void mapCombine22(std::unique_ptr<StructuredWindowMap<nDim>>& input1,
                    std::unique_ptr<StructuredWindowMap<nDim>>& input2,
                    std::unique_ptr<StructuredWindowMap<nDim>>& output);
  void mapCombineWave(std::unique_ptr<StructuredWindowMap<nDim>>& input1,
                      std::unique_ptr<StructuredWindowMap<nDim>>& input2,
                      std::unique_ptr<StructuredWindowMap<nDim>>& output);
  void mapCombineCell11(std::unique_ptr<StructuredWindowMap<nDim>>& input1,
                        std::unique_ptr<StructuredWindowMap<nDim>>& input2,
                        std::unique_ptr<StructuredWindowMap<nDim>>& output);
  void mapCombineCell12(std::unique_ptr<StructuredWindowMap<nDim>>& input1,
                        std::unique_ptr<StructuredWindowMap<nDim>>& input2,
                        std::unique_ptr<StructuredWindowMap<nDim>>& output);
  void mapCombineCell21(std::unique_ptr<StructuredWindowMap<nDim>>& input1,
                        std::unique_ptr<StructuredWindowMap<nDim>>& input2,
                        std::unique_ptr<StructuredWindowMap<nDim>>& output);
  void mapCombineCell22(std::unique_ptr<StructuredWindowMap<nDim>>& input1,
                        std::unique_ptr<StructuredWindowMap<nDim>>& input2,
                        std::unique_ptr<StructuredWindowMap<nDim>>& output);

  void mapCpy(std::unique_ptr<StructuredWindowMap<nDim>>& input, std::unique_ptr<StructuredWindowMap<nDim>>& output);
  void readMapFromArray(std::unique_ptr<StructuredWindowMap<nDim>>& map, MInt* array);
  void writeMapToArray(std::unique_ptr<StructuredWindowMap<nDim>>& map, MInt* array);

  void initGlobals();
  void writeConnectionWindowInformation3D(MFloat* periodicDisplacements);
  void writeConnectionWindowInformation2D(MFloat* /*periodicDisplacements*/){};
  void readConnectionWindowInformation3D(MFloat* periodicDisplacments);
  void readConnectionWindowInformation2D(MFloat* /*periodicDisplacments*/){};
  void createWindowMapping(MPI_Comm* channelIn, MPI_Comm* channelOut, MPI_Comm* channelWorld, MInt* channelRoots,
                           MPI_Comm* commStg, MInt* commStgRoot, MInt* commStgRootGlobal, MPI_Comm* commBC2600,
                           MInt* commBC2600Root, MInt* commBC2600RootGlobal, MPI_Comm* rescalingCommGrComm,
                           MInt* rescalingCommGrRoot, MInt* rescalingCommGrRootGlobal, MPI_Comm* commPerRotOne,
                           MPI_Comm* commPerRotTwo, MPI_Comm* commPerRotWorld, MInt* rotationRoots, MInt& perRotGroup,
                           SingularInformation* singularity, MInt* hasSingularity, MPI_Comm* plenumComm,
                           MInt* plenumRoots);
  void createWaveWindowMapping(MInt);
  void createCommunicationExchangeFlags(std::vector<std::unique_ptr<StructuredComm<nDim>>>&,
                                        std::vector<std::unique_ptr<StructuredComm<nDim>>>&, const MInt,
                                        MFloat* const* const);
  void createWaveCommunicationExchangeFlags(std::vector<std::unique_ptr<StructuredComm<nDim>>>&,
                                            std::vector<std::unique_ptr<StructuredComm<nDim>>>&,
                                            const MInt noVariables);
  std::vector<std::unique_ptr<StructuredWindowMap<nDim>>> localStructuredBndryCndMaps;
  std::vector<std::unique_ptr<StructuredWindowMap<nDim>>> channelSurfaceIndices;
  std::vector<std::unique_ptr<StructuredWindowMap<nDim>>> plenumInletSurfaceIndices;
  // for the sponge
  std::vector<std::unique_ptr<StructuredWindowMap<nDim>>> m_spongeInfoMap;
  std::vector<std::unique_ptr<StructuredWindowMap<nDim>>> m_wallDistInfoMap;
  // zonal
  std::vector<std::unique_ptr<StructuredWindowMap<nDim>>> m_zonalBCMaps;
  MBool checkZonalBCMaps(std::unique_ptr<StructuredWindowMap<nDim>>&, std::unique_ptr<StructuredWindowMap<nDim>>&);
  //
  void createAuxDataMap(const MInt, const MString, const std::vector<MInt>&, const std::vector<MInt>&, const MBool);


  //=============================================================================================
  // new functions and variables
  void readWindowCoordinates(MFloat* periodicDisplacements);
  void deleteDuplicateWindows(std::vector<std::unique_ptr<StructuredWindowMap<nDim>>>&,
                              std::vector<std::unique_ptr<StructuredWindowMap<nDim>>>&);
  void deleteDuplicateCommMaps();
  void deleteDuplicateBCMaps(std::vector<std::unique_ptr<StructuredWindowMap<nDim>>>&);
  void setBCsingular(std::unique_ptr<StructuredWindowMap<nDim>>&, std::unique_ptr<StructuredWindowMap<nDim>>&,
                     const MInt);
  void periodicPointsChange(MFloat* pt, MInt type, MFloat* periodicDisplacements); // transform point coords

  // for mutilsolver and periodic connections
  MBool addConnection(MInt connectiontype, MInt b1, MInt* p1, MInt b2, MInt* p2);
  MBool findConnection(connectionNode a);
  void removeConnection(connectionNode a);

  // for special connecitons
  MBool addConnection(MInt connectiontype, MInt b1, MInt* p1, MInt b2, MInt* p2, MInt Nstar);
  MBool findConnection(connectionNode a, MInt Nstar);
  void removeConnection(connectionNode a, MInt Nstar);
  void multiBlockAssembling();
  void singularityAssembling();

  std::multiset<connectionNode> connectionset;         // CHANGE_SET
  std::multiset<connectionNode> singularconnectionset; // CHANGE_SET
  std::vector<std::unique_ptr<StructuredWindowMap<nDim>>> window2d;
  std::vector<std::unique_ptr<StructuredWindowMap<nDim>>> window1d;
  std::vector<std::unique_ptr<StructuredWindowMap<nDim>>> window0d;
  std::vector<std::unique_ptr<StructuredWindowMap<nDim>>> singularwindow;
  std::vector<std::unique_ptr<StructuredWindowMap<nDim>>> localSingularMap;
  //==============================================================================================

  // Position in globalStructuredBndryCndMaps and windowId of all windows for which auxData is requested according to
  // property file
  std::map<MInt, MInt> m_auxDataWindowIds;

 private:
  StructuredGrid<nDim>* m_grid;

  std::vector<std::unique_ptr<StructuredWindowMap<nDim>>> globalStructuredBndryCndMaps;
  std::vector<std::unique_ptr<StructuredWindowMap<nDim>>> localStructuredDomainMaps;
  std::vector<std::unique_ptr<StructuredWindowMap<nDim>>> localStructuredDomainBndryMaps;

  // diagonal communication
  std::vector<std::unique_ptr<StructuredWindowMap<nDim>>> m_partitionMapsWithGC;
  std::vector<std::unique_ptr<StructuredWindowMap<nDim>>> m_partitionMapsWithoutGC;
  std::vector<std::unique_ptr<StructuredWindowMap<nDim>>> rcvMap;
  std::vector<std::unique_ptr<StructuredWindowMap<nDim>>> sndMap;
  std::vector<std::unique_ptr<StructuredWindowMap<nDim>>> rcvMapPeriodic;
  std::vector<std::unique_ptr<StructuredWindowMap<nDim>>> sndMapPeriodic;
  std::vector<std::unique_ptr<StructuredWindowMap<nDim>>> physicalBCMap;
  std::vector<std::unique_ptr<StructuredWindowMap<nDim>>> physicalAuxDataMap;
  std::vector<std::unique_ptr<StructuredWindowMap<nDim>>> waveRcvMap;
  std::vector<std::unique_ptr<StructuredWindowMap<nDim>>> waveSndMap;
  std::vector<std::unique_ptr<windowInformation<nDim>>> inputWindows;
  MInt noInputWindowInformation;
  MInt noInputWindowConnections;
  MInt noInputBndryCnds;
  std::unique_ptr<StructuredWindowMap<nDim>> m_myMapWithGC = nullptr;
  std::unique_ptr<StructuredWindowMap<nDim>> m_myMapWithoutGC = nullptr;

  // Communication related
  MInt noDomains() const { return m_noDomains; }
  MInt domainId() const { return m_domainId; }
  MInt solverId() const { return m_solverId; }

  const MInt m_noBlocks;
  const MInt m_blockId;
  const MInt m_noDomains;
  const MInt m_domainId;
  const MInt m_solverId;
  MInt m_noGhostLayers;
  MPI_Comm m_StructuredComm;
};

#endif
