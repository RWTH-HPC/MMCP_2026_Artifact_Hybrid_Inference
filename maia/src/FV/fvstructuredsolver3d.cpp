// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "fvstructuredsolver3d.h"
#include <cmath>
#include "COMM/mpioverride.h"
#include "INCLUDE/maiaconstants.h"
#include "IO/parallelio_hdf5.h"
#include "UTIL/parallelfor.h"
#include "globals.h"

using namespace std;

FvStructuredSolver3D::FvStructuredSolver3D(MInt solverId, StructuredGrid<3>* grid_, MBool* propertiesGroups,
                                           const MPI_Comm comm)
  : FvStructuredSolver<3>(solverId, grid_, propertiesGroups, comm) {
  TRACE();
  const MLong oldAllocatedBytes = allocatedBytes();

  // count the no of necessary FQ fields and allocate
  initializeFQField();

  // compute the cell center coordinates from point coordinates
  m_grid->computeCellCenterCoordinates();

  if(m_rans) {
    m_structuredBndryCnd = make_unique<StructuredBndryCnd3D<true>>(this, m_grid);
  } else {
    m_structuredBndryCnd = make_unique<StructuredBndryCnd3D<false>>(this, m_grid);
  }

  allocateSingularities();

  // allocate memory for aux data maps (cf,cp)
  allocateAuxDataMaps();

  // assign coordinates to all ghost points
  addGhostPointCoordinateValues();

  // if we are RANS we should allocate a new RANS solver
  if(m_rans == true) {
    m_ransSolver = make_unique<FvStructuredSolver3DRans>(this);
  }

  initFluxMethod();

  m_grid->computeCellCenterCoordinates();
  RECORD_TIMER_START(m_timers[Timers::ComputeMetrics]);
  m_grid->computeMetrics();
  RECORD_TIMER_STOP(m_timers[Timers::ComputeMetrics]);
  RECORD_TIMER_START(m_timers[Timers::ComputeJacobian]);
  m_grid->computeJacobian();
  RECORD_TIMER_STOP(m_timers[Timers::ComputeJacobian]);

  // TODO_SS labels:FV,totest By now I am not sure if Code performs correctly for wrongly oriented meshes
  for(MInt k = 0; k < m_nCells[0]; k++) {
    for(MInt j = 0; j < m_nCells[1]; j++) {
      for(MInt i = 0; i < m_nCells[2]; i++) {
        const MInt cellId = cellIndex(i, j, k);
        if(m_cells->cellJac[cellId] < 0.0) mTerm(1, "Negative Jacobian found! Check first if code can cope this!");
      }
    }
  }

  m_convergence = false;

  if(m_zonal) {
    computeZonalConnections();
  }

  // Assign handlers to the correct boundary conditions
  assignBndryCells();

  // Computation of modified wall distance in porous computation requires to set the porosity,
  // Da-number etc. first; on the other hand we need to wait for initializeFQField to be called
  if(m_porous) {
    initPorous();
    // exchange6002();
  }

  if(m_zonal || m_rans) {
    if(m_ransMethod == RANS_SA_DV) {
      m_structuredBndryCnd->computeWallDistances();
    }
  }

  if(m_intpPointsOutputInterval > 0) {
    initInterpolatedPoints();
  }

  if(m_pointsToAsciiOutputInterval > 0) {
    initPointsToAsciiFile();
  }

  printAllocatedMemory(oldAllocatedBytes, "FvStructuredSolver3D", m_StructuredComm);

  RECORD_TIMER_STOP(m_timers[Timers::Constructor]);
}

FvStructuredSolver3D::~FvStructuredSolver3D() { TRACE(); }

void FvStructuredSolver3D::computeZonalConnections() {
  if(!m_rans) {
    // do a spanwise average of the pressure in
    // all LES zones
    m_zonalSpanwiseAvgVars.push_back(m_cells->fq[FQ->AVG_P]);
  }

  //////////////////////////////////////
  ////Find the number of ZonalBC////////
  //////////////////////////////////////
  const MInt noZonalBCMaps = m_windowInfo->m_zonalBCMaps.size();
  m_zonalBC.resize(noZonalBCMaps);

  ///////////////////////////////////
  ///// Start ZonalBC Loop //////////
  ///////////////////////////////////
  if(domainId() == 0) {
    cout << "///////////////////////////////////////////////////////////////////" << endl
         << "Starting zonal bc creation for " << noZonalBCMaps << " zonal bc maps" << endl;
  }

  for(MInt id = 0; id < noZonalBCMaps; id++) {
    if(domainId() == 0) {
      cout << "///////////////////////////////////////////////////////////////////" << endl;
      cout << "//////////// ZONAL BC " << id << " ///////////////////////////////////////////" << endl;
      cout << "///////////////////////////////////////////////////////////////////" << endl;
    }

    m_zonalBC[id] = make_unique<StructuredZonalBC>();
    m_zonalBC[id]->receiverBlockId = m_windowInfo->m_zonalBCMaps[id]->Id1;

    if(domainId() == 0) {
      cout << "Zonal BC" << id << ": Initialization ..." << endl;
    }

    // set up the data of the zonal domain and noCellsBC
    MIntScratchSpace noZonalCells(noDomains(), AT_, "noZonalCells");
    MInt hasRcvDomain = 0;
    m_zonalBC[id]->noZonalVariables = 6;
    m_zonalBC[id]->hasSTG = false;
    m_zonalBC[id]->noCellsGlobalBC = 0;
    m_zonalBC[id]->noCellsLocalBC = 0;
    m_zonalBC[id]->hasLocalBCMap = false;

    for(MInt bcId = 0; bcId < abs((MInt)m_windowInfo->physicalBCMap.size()); ++bcId) {
      m_zonalBC[id]->hasLocalBCMap =
          m_windowInfo->checkZonalBCMaps(m_windowInfo->m_zonalBCMaps[id], m_windowInfo->physicalBCMap[bcId]);

      MInt stgIP = 0;
      if(m_zonalBC[id]->hasLocalBCMap) {
        if(m_windowInfo->physicalBCMap[bcId]->BC == 2221) {
          stgIP = 1; // in this case use one more cell, required for STG
          m_zonalBC[id]->hasSTG = true;
        }

        MInt* start = m_structuredBndryCnd->m_physicalBCMap[bcId]->start1;
        MInt* end = m_structuredBndryCnd->m_physicalBCMap[bcId]->end1;
        m_zonalBC[id]->noCellsLocalBC += (end[0] + stgIP - start[0]) * (end[1] - start[1]) * (end[2] - start[2]);
        m_zonalBC[id]->start[0] = start[0];
        m_zonalBC[id]->start[1] = start[1];
        m_zonalBC[id]->start[2] = start[2];
        m_zonalBC[id]->end[0] = end[0] + stgIP;
        m_zonalBC[id]->end[1] = end[1];
        m_zonalBC[id]->end[2] = end[2];
        hasRcvDomain = 1;
        break;
      }
    }

    if(domainId() == 0) {
      cout << "Zonal BC" << id << ": Initialization ... Finished!" << endl;
      cout << "Zonal BC" << id << ": Collecting local coordinates ... " << endl;
    }

    MFloatScratchSpace localCoordinatesBC(nDim, std::max(1, m_zonalBC[id]->noCellsLocalBC), AT_, "localCoordinatesBC");
    MIntScratchSpace localReceiverIds(std::max(1, m_zonalBC[id]->noCellsLocalBC), AT_, "localReceiverIds");
    MIntScratchSpace localMapCellIds(std::max(1, m_zonalBC[id]->noCellsLocalBC), AT_, "m_localMapCellIds");

    if(m_zonalBC[id]->hasLocalBCMap) {
      MInt bcCellId = 0;
      for(MInt k = m_zonalBC[id]->start[2]; k < m_zonalBC[id]->end[2]; k++) {
        for(MInt j = m_zonalBC[id]->start[1]; j < m_zonalBC[id]->end[1]; j++) {
          for(MInt i = m_zonalBC[id]->start[0]; i < m_zonalBC[id]->end[0]; i++) {
            const MInt cellId = cellIndex(i, j, k);

            if(m_zonalBC[id]->hasSTG) {
              localMapCellIds[bcCellId] = bcCellId;
            } else {
              localMapCellIds[bcCellId] = cellId;
            }

            for(MInt dim = 0; dim < nDim; ++dim) {
              localCoordinatesBC(dim, bcCellId) = m_cells->coordinates[dim][cellId];
            }

            localReceiverIds[bcCellId] = domainId();
            ++bcCellId;
          }
        }
      }
    }

    if(domainId() == 0) {
      cout << "Zonal BC" << id << ": Collecting local coordinates ... Finished!" << endl;
      cout << "Zonal BC" << id << ": Exchanging information about local coordinates ..." << endl;
    }

    MPI_Allreduce(&m_zonalBC[id]->noCellsLocalBC, &m_zonalBC[id]->noCellsGlobalBC, 1, MPI_INT, MPI_SUM,
                  m_StructuredComm, AT_, "m_zonalBC[id]->noCellsLocalBC", "m_zonalBC[id]->noCellsGlobalBC");
    MPI_Allgather(&m_zonalBC[id]->noCellsLocalBC, 1, MPI_INT, &noZonalCells[0], 1, MPI_INT, m_StructuredComm, AT_,
                  "m_zonalBC[id]->noCellsLocalBC", "noZonalCells[0]");

    if(domainId() == 0) {
      cout << "Zonal BC" << id << ": Exchanging information about local coordinates ... Finished!" << endl;
    }

    //////////////////////////////////////////////////////////////////////////
    ////// Gathering the cellCoordinatesGlobalBC from CoordinatesLocalBC//////
    //////////////////////////////////////////////////////////////////////////
    if(domainId() == 0) {
      MInt noCells = 0;
      for(MInt i = 0; i < noDomains(); i++) {
        noCells += noZonalCells[i];
      }
      cout << "Zonal BC" << id << ": Gathering coordinates for " << noCells << " cells on all partitions..." << endl;
    }

    MIntScratchSpace noCellsBCOffsets(noDomains(), AT_, "noCellsBCOffsets");
    mAlloc(m_zonalBC[id]->coordinatesGlobalBC, 3, std::max(1, m_zonalBC[id]->noCellsGlobalBC), "coordinatesGlobalBC",
           0.0, AT_);
    mAlloc(m_zonalBC[id]->globalLocalMapCellIds, std::max(1, m_zonalBC[id]->noCellsGlobalBC), "globalLocalMapCellIds",
           -1, AT_);
    mAlloc(m_zonalBC[id]->globalReceiverIds, std::max(1, m_zonalBC[id]->noCellsGlobalBC), "globalReceiverIds", -1, AT_);

    noCellsBCOffsets[0] = 0;
    for(MInt i = 1; i < noDomains(); i++) {
      noCellsBCOffsets[i] = noCellsBCOffsets[i - 1] + noZonalCells[i - 1];
    }


    // Distribute all zonal cell information (all coordinates, cellIds, domainIds) to all domains
    for(MInt dim = 0; dim < nDim; ++dim) {
      MPI_Allgatherv(&localCoordinatesBC(dim, 0), noZonalCells[domainId()], MPI_DOUBLE,
                     &m_zonalBC[id]->coordinatesGlobalBC[dim][0], &noZonalCells[0], &noCellsBCOffsets[0], MPI_DOUBLE,
                     m_StructuredComm, AT_, "localCoordinatesBC(0,0)", "coordinatesGlobalBC[0][0]");
    }
    MPI_Allgatherv(&localReceiverIds[0], noZonalCells[domainId()], MPI_INT, &m_zonalBC[id]->globalReceiverIds[0],
                   &noZonalCells[0], &noCellsBCOffsets[0], MPI_INT, m_StructuredComm, AT_, "localReceiverIds[0]",
                   "m_zonalBC[id]->globalReceiverIds[0]");
    MPI_Allgatherv(&localMapCellIds[0], noZonalCells[domainId()], MPI_INT, &m_zonalBC[id]->globalLocalMapCellIds[0],
                   &noZonalCells[0], &noCellsBCOffsets[0], MPI_INT, m_StructuredComm, AT_, "localMapCellIds[0]",
                   "m_zonalBC[id]->globalLocalMapCellIds[0]");

    if(domainId() == 0) {
      MInt noCells = 0;
      for(MInt i = 0; i < noDomains(); i++) {
        noCells += noZonalCells[i];
      }
      cout << "Zonal BC" << id << ": Gathering coordinates for " << noCells << " cells on all partitions... Finished!"
           << endl;
      cout << "Zonal BC" << id << ": Searching for suitable interpolation partners..." << endl;
    }
    /////////////////////////////////////////
    ////Set Up PrepareInterpolation  ////////
    /////////////////////////////////////////
    MIntScratchSpace hasPartnerLocalBC(m_zonalBC[id]->noCellsGlobalBC, AT_, "hasPartnerLocalBC");
    MBool hasInterpolationPartnerDomain = true;

    if(m_blockId == m_zonalBC[id]->receiverBlockId) {
      hasInterpolationPartnerDomain = false;
    }

    m_zonalBC[id]->interpolation =
        make_unique<StructuredInterpolation<nDim>>(m_nCells, m_cells->coordinates, m_StructuredComm);

    m_zonalBC[id]->interpolation->prepareZonalInterpolation(m_zonalBC[id]->noCellsGlobalBC,
                                                            m_zonalBC[id]->coordinatesGlobalBC,
                                                            hasPartnerLocalBC.begin(),
                                                            hasInterpolationPartnerDomain);

    if(domainId() == 0) {
      cout << "Zonal BC" << id << ": Searching for suitable interpolation partners... Finished!" << endl;
      cout << "Zonal BC" << id << ": Creating sending and receiving information..." << endl;
    }

    /////////////////////////////////////////////////////////////////////////////////////////
    //////////////////// global sndDomain,rcvDomain Info in each zonalBC///////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////
    // to find out which domains will be rcvDomains

    m_zonalBC[id]->noGlobalRcvDomains = 0;
    MPI_Allreduce(&hasRcvDomain, &m_zonalBC[id]->noGlobalRcvDomains, 1, MPI_INT, MPI_SUM, m_StructuredComm, AT_,
                  "hasRcvDomain", "m_zonalBC[id]->noGlobalRcvDomains");

    mAlloc(m_zonalBC[id]->globalRcvDomainIds, std::max(1, m_zonalBC[id]->noGlobalRcvDomains), "m_globalRcvDomainIds",
           -1, AT_);

    MInt pos = 0;
    for(MInt j = 0; j < noDomains(); j++) {
      if(noZonalCells[j] > 0) {
        m_zonalBC[id]->globalRcvDomainIds[pos] = j;
        ++pos;
      }
    }

    // to find out which domains will be sndDomains
    MInt hasSndDomain = false;
    for(MInt cellIdBC = 0; cellIdBC < m_zonalBC[id]->noCellsGlobalBC; cellIdBC++) {
      if(hasPartnerLocalBC[cellIdBC]) {
        hasSndDomain = true;
        break;
      }
    }

    m_zonalBC[id]->noGlobalSndDomains = 0;
    MPI_Allreduce(&hasSndDomain, &m_zonalBC[id]->noGlobalSndDomains, 1, MPI_INT, MPI_SUM, m_StructuredComm, AT_,
                  "hasSndDomain", "m_zonalBC[id]->noGlobalSndDomains");


    MIntScratchSpace hasSndDomainInfo(noDomains(), AT_, "hasSndDomainInfo");
    MPI_Allgather(&hasSndDomain, 1, MPI_INT, &hasSndDomainInfo[0], 1, MPI_INT, m_StructuredComm, AT_, "hasSndDomain",
                  "hasSndDomainInfo[0]");

    mAlloc(m_zonalBC[id]->globalSndDomainIds, std::max(1, m_zonalBC[id]->noGlobalSndDomains), "m_globalSndDomainIds",
           -1, AT_);

    pos = 0;
    for(MInt j = 0; j < noDomains(); j++) {
      if(hasSndDomainInfo[j] > 0) {
        m_zonalBC[id]->globalSndDomainIds[pos] = j;
        pos++;
      }
    }
    /////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////local send and receive domains Infomation//////////////////
    /////////////////////////////////////////////////////////////////////////////////////////
    // compute no of cells which our domain can interpolate and has to send
    m_zonalBC[id]->noBufferSndSize = 0;
    for(MInt cellIdBC = 0; cellIdBC < m_zonalBC[id]->noCellsGlobalBC; cellIdBC++) {
      if(hasPartnerLocalBC[cellIdBC]) {
        m_zonalBC[id]->noBufferSndSize++;
      }
    }

    mAlloc(m_zonalBC[id]->localCommReceiverIds, std::max(1, m_zonalBC[id]->noBufferSndSize), "localCommReceiverIds",
           AT_);
    mAlloc(m_zonalBC[id]->localBufferMapCellIds, std::max(1, m_zonalBC[id]->noBufferSndSize), "localCommReceiverIds",
           AT_);
    mAlloc(m_zonalBC[id]->localBufferIndexCellIds, std::max(1, m_zonalBC[id]->noBufferSndSize), "localCommReceiverIds",
           AT_);

    // collect the receiver domainIds for each cell that our domain can interpolate
    pos = 0;
    for(MInt cellIdBC = 0; cellIdBC < m_zonalBC[id]->noCellsGlobalBC; cellIdBC++) {
      if(hasPartnerLocalBC[cellIdBC]) {
        m_zonalBC[id]->localCommReceiverIds[pos] = m_zonalBC[id]->globalReceiverIds[cellIdBC];
        pos++;
      }
    }

    // count to how many different domains our domain needs to send data
    m_zonalBC[id]->noSndNghbrDomains = 0;
    for(MInt i = 0; i < m_zonalBC[id]->noGlobalSndDomains; i++) {
      if(m_zonalBC[id]->globalSndDomainIds[i] == domainId()) {
        for(MInt d = 0; d < noDomains(); d++) {
          for(MInt j = 0; j < m_zonalBC[id]->noBufferSndSize; j++) {
            if(m_zonalBC[id]->localCommReceiverIds[j] == d) {
              m_zonalBC[id]->noSndNghbrDomains++;
              break;
            }
          }
        }
      }
    }

    MIntScratchSpace localRcvId(std::max(1, m_zonalBC[id]->noSndNghbrDomains), AT_, "localRcvId");
    MIntScratchSpace localBufferSndSize(std::max(1, m_zonalBC[id]->noSndNghbrDomains), AT_, "localBufferSndSize");

    // save each receiving domainId in localRcvId
    pos = 0;
    for(MInt i = 0; i < m_zonalBC[id]->noGlobalSndDomains; i++) {
      if(m_zonalBC[id]->globalSndDomainIds[i] == domainId()) {
        for(MInt d = 0; d < noDomains(); d++) {
          for(MInt j = 0; j < m_zonalBC[id]->noBufferSndSize; j++) {
            if(m_zonalBC[id]->localCommReceiverIds[j] == d) {
              localRcvId[pos] = d;
              pos++;
              break;
            }
          }
        }
      }
    }


    // save the cellId and mapCellId (the cellId of the receiving cell on the nghbr Domain)
    pos = 0;
    for(MInt cellIdBC = 0; cellIdBC < m_zonalBC[id]->noCellsGlobalBC; cellIdBC++) {
      if(hasPartnerLocalBC[cellIdBC]) {
        const MInt mapCellId = m_zonalBC[id]->globalLocalMapCellIds[cellIdBC];
        m_zonalBC[id]->localBufferMapCellIds[pos] = mapCellId;
        m_zonalBC[id]->localBufferIndexCellIds[pos] = cellIdBC; // this is needed to get InterpolatedVars in zonalGather
        pos++;
      }
    }

    // compute the buffer size for each domain our domain needs to send data to
    for(MInt i = 0; i < m_zonalBC[id]->noSndNghbrDomains; i++) {
      for(MInt j = 0; j < m_zonalBC[id]->noBufferSndSize; j++) {
        if(localRcvId[i] == m_zonalBC[id]->localCommReceiverIds[j]) {
          localBufferSndSize[i]++;
        }
      }
    }

    // get m_zonalBC[id]->noRcvNghbrDomains and SndId
    MIntScratchSpace sndBufferRcvSize(noDomains(), AT_, "sndBufferRcvSize");

    for(MInt noRcv = 0; noRcv < m_zonalBC[id]->noGlobalRcvDomains; noRcv++) {
      MInt bufferRcvSize = 0;
      for(MInt cellIdBC = 0; cellIdBC < m_zonalBC[id]->noCellsGlobalBC; cellIdBC++) {
        if(m_zonalBC[id]->globalRcvDomainIds[noRcv] == m_zonalBC[id]->globalReceiverIds[cellIdBC]) {
          if(hasPartnerLocalBC[cellIdBC]) {
            bufferRcvSize++;
          }
        }
      }

      // get the bufferSndSize
      sndBufferRcvSize[m_zonalBC[id]->globalRcvDomainIds[noRcv]] = bufferRcvSize;
    }

    MIntScratchSpace rcvBufferRcvSize(noDomains(), AT_, "rcvBufferRcvSize");
    MPI_Alltoall(&sndBufferRcvSize[0], 1, MPI_INT, &rcvBufferRcvSize[0], 1, MPI_INT, m_StructuredComm, AT_,
                 "sndBufferRcvSize[0]", "rcvBufferRcvSize[0]");

    // get m_zonalBC[id]->noRcvNghbrDomains, localSndId, localbufferRcvSize.
    // These are from which domains each rcvDoamin will receive
    // and  how much size it is.
    m_zonalBC[id]->noRcvNghbrDomains = 0;
    if(m_zonalBC[id]->hasLocalBCMap) {
      for(MInt i = 0; i < noDomains(); i++) {
        if(rcvBufferRcvSize[i] > 0) { // get noRcvNghbrDoamins
          m_zonalBC[id]->noRcvNghbrDomains++;
        }
      }
    }

    MIntScratchSpace localSndId(std::max(1, m_zonalBC[id]->noRcvNghbrDomains), AT_, "localSndId");
    MIntScratchSpace localBufferRcvSize(std::max(1, m_zonalBC[id]->noRcvNghbrDomains), AT_, "localBufferRcvSize");

    // get localSndId from which domains each rcvDoamin will receive data
    pos = 0;
    if(m_zonalBC[id]->hasLocalBCMap) {
      for(MInt i = 0; i < noDomains(); i++) {
        if(rcvBufferRcvSize[i] > 0) { // sndId is the info the receiver Domains have from which domains I will receive
          localSndId[pos] = i;
          pos++;
        }
      }
    }
    // get localbufferRcvSize
    pos = 0;
    if(m_zonalBC[id]->hasLocalBCMap) {
      for(MInt i = 0; i < noDomains(); i++) {
        if(rcvBufferRcvSize[i] > 0) {
          localBufferRcvSize[pos] = rcvBufferRcvSize[i];
          pos++;
        }
      }
    }

    //////////////////////////////////////////////////////
    ///////////// SET UP COMMUNCATION ////////////////////
    //////////////////////////////////////////////////////

    mAlloc(m_zonalBC[id]->sndRequest, std::max(1, m_zonalBC[id]->noSndNghbrDomains), "sndRequest", AT_);
    mAlloc(m_zonalBC[id]->sndStatus, std::max(1, m_zonalBC[id]->noSndNghbrDomains), "sndStatus", AT_);
    mAlloc(m_zonalBC[id]->rcvRequest, std::max(1, m_zonalBC[id]->noRcvNghbrDomains), "rcvRequest", AT_);
    mAlloc(m_zonalBC[id]->rcvStatus, std::max(1, m_zonalBC[id]->noRcvNghbrDomains), "rcvStatus", AT_);


    for(MInt i = 0; i < m_zonalBC[id]->noSndNghbrDomains; i++) {
      unique_ptr<StructuredZonalComm> sndCommPtr =
          make_unique<StructuredZonalComm>(localBufferSndSize[i], localRcvId[i], m_zonalBC[id]->noZonalVariables);
      m_zonalBC[id]->sndComm.push_back(std::move(sndCommPtr));
      m_zonalBC[id]->sndRequest[i] = MPI_REQUEST_NULL;
    }

    for(MInt i = 0; i < m_zonalBC[id]->noRcvNghbrDomains; i++) {
      unique_ptr<StructuredZonalComm> rcvCommPtr =
          make_unique<StructuredZonalComm>(localBufferRcvSize[i], localSndId[i], m_zonalBC[id]->noZonalVariables);
      m_zonalBC[id]->rcvComm.push_back(std::move(rcvCommPtr));
      m_zonalBC[id]->rcvRequest[i] = MPI_REQUEST_NULL;
    }


    ///////////////////////////////////////////////////////////////
    /////////////////////MapCellId exchange////////////////////////
    ///////////////////////////////////////////////////////////////

    for(MInt noSnd = 0; noSnd < m_zonalBC[id]->noSndNghbrDomains; noSnd++) {
      pos = 0;

      for(MInt cellId = 0; cellId < m_zonalBC[id]->noBufferSndSize; cellId++) {
        if(m_zonalBC[id]->localCommReceiverIds[cellId] == m_zonalBC[id]->sndComm[noSnd]->localId) {
          m_zonalBC[id]->sndComm[noSnd]->mapCellId[pos] = m_zonalBC[id]->localBufferMapCellIds[cellId];
          pos++;
        }
      }

      MInt err = MPI_Isend((void*)m_zonalBC[id]->sndComm[noSnd]->mapCellId, m_zonalBC[id]->sndComm[noSnd]->bufferSize,
                           MPI_INT, m_zonalBC[id]->sndComm[noSnd]->localId, 0, m_StructuredComm,
                           &m_zonalBC[id]->sndRequest[noSnd], AT_, "m_zonalBC[id]->sndComm[noSnd]->intBbuffer");
      if(err)
        cout << "rank " << domainId() << " zonal sending to " << m_zonalBC[id]->sndComm[noSnd]->localId
             << " throws error " << endl;
    }

    for(MInt noRcv = 0; noRcv < m_zonalBC[id]->noRcvNghbrDomains; noRcv++) {
      MInt err =
          MPI_Irecv((void*)&m_zonalBC[id]->rcvComm[noRcv]->mapCellId[0], m_zonalBC[id]->rcvComm[noRcv]->bufferSize,
                    MPI_INT, m_zonalBC[id]->rcvComm[noRcv]->localId, 0, m_StructuredComm,
                    &m_zonalBC[id]->rcvRequest[noRcv], AT_, "(void*)&m_zonalBC[id]->rcvComm[noRcv]->mapCellId");
      if(err)
        cout << "rank " << domainId() << " zonal sending to " << m_zonalBC[id]->sndComm[noRcv]->localId
             << " throws error " << endl;
    }

    MPI_Waitall(m_zonalBC[id]->noSndNghbrDomains, m_zonalBC[id]->sndRequest, m_zonalBC[id]->sndStatus, AT_);
    MPI_Waitall(m_zonalBC[id]->noRcvNghbrDomains, m_zonalBC[id]->rcvRequest, m_zonalBC[id]->rcvStatus, AT_);

    MPI_Barrier(m_StructuredComm, AT_);

    if(domainId() == 0) {
      cout << "Zonal BC" << id << ": Creating sending and receiving information... Finished!" << endl;
      cout << "Zonal BC" << id << ": Completed!" << endl;
    }
  }
}

void FvStructuredSolver3D::zonalExchange() {
  zonalGather();
  zonalSend();
  zonalReceive();
  zonalScatter();
}


void FvStructuredSolver3D::zonalGather() {
  // gather sending variables in buffer
  for(auto const& zBC : m_zonalBC) {
    for(MInt noSnd = 0; noSnd < zBC->noSndNghbrDomains; noSnd++) {
      const MInt noCells = zBC->sndComm[noSnd]->bufferSize;
      MInt pos = 0;
      for(MInt cellId = 0; cellId < zBC->noBufferSndSize; cellId++) {
        if(zBC->localCommReceiverIds[cellId] == zBC->sndComm[noSnd]->localId) {
          const MInt cellIdBC = zBC->localBufferIndexCellIds[cellId];
          if(m_rans) {
            zBC->sndComm[noSnd]->buffer[pos + 0 * noCells] =
                zBC->interpolation->interpolateVariableZonal(m_cells->pvariables[PV->U], cellIdBC);
            zBC->sndComm[noSnd]->buffer[pos + 1 * noCells] =
                zBC->interpolation->interpolateVariableZonal(m_cells->pvariables[PV->V], cellIdBC);
            zBC->sndComm[noSnd]->buffer[pos + 2 * noCells] =
                zBC->interpolation->interpolateVariableZonal(m_cells->pvariables[PV->W], cellIdBC);
            zBC->sndComm[noSnd]->buffer[pos + 3 * noCells] =
                zBC->interpolation->interpolateVariableZonal(m_cells->pvariables[PV->RHO], cellIdBC);
            zBC->sndComm[noSnd]->buffer[pos + 4 * noCells] =
                zBC->interpolation->interpolateVariableZonal(m_cells->pvariables[PV->P], cellIdBC);
            zBC->sndComm[noSnd]->buffer[pos + 5 * noCells] =
                zBC->interpolation->interpolateVariableZonal(m_cells->fq[FQ->NU_T], cellIdBC);
          } else {
            zBC->sndComm[noSnd]->buffer[pos + 0 * noCells] =
                zBC->interpolation->interpolateVariableZonal(m_cells->fq[FQ->AVG_U], cellIdBC);
            zBC->sndComm[noSnd]->buffer[pos + 1 * noCells] =
                zBC->interpolation->interpolateVariableZonal(m_cells->fq[FQ->AVG_V], cellIdBC);
            zBC->sndComm[noSnd]->buffer[pos + 2 * noCells] =
                zBC->interpolation->interpolateVariableZonal(m_cells->fq[FQ->AVG_W], cellIdBC);
            zBC->sndComm[noSnd]->buffer[pos + 3 * noCells] =
                zBC->interpolation->interpolateVariableZonal(m_cells->fq[FQ->AVG_RHO], cellIdBC);
            zBC->sndComm[noSnd]->buffer[pos + 4 * noCells] =
                zBC->interpolation->interpolateVariableZonal(m_cells->fq[FQ->AVG_P], cellIdBC);
            zBC->sndComm[noSnd]->buffer[pos + 5 * noCells] =
                zBC->interpolation->interpolateVariableZonal(m_cells->fq[FQ->NU_T], cellIdBC);
          }
          pos++;
        }
      }
    }
  }
}

void FvStructuredSolver3D::zonalSend() {
  for(auto const& zBC : m_zonalBC) {
    // access only sndDomains
    for(MInt i = 0; i < zBC->noGlobalSndDomains; i++) {
      if(zBC->globalSndDomainIds[i] == domainId()) {
        for(MInt noSnd = 0; noSnd < zBC->noSndNghbrDomains; noSnd++) {
          MInt err = MPI_Isend(&zBC->sndComm[noSnd]->buffer[0], zBC->sndComm[noSnd]->bufferSize * zBC->noZonalVariables,
                               MPI_DOUBLE, zBC->sndComm[noSnd]->localId, 0, m_StructuredComm, &zBC->sndRequest[noSnd],
                               AT_, "zBC->buffer[noSnd][0]");
          if(err) cout << "rank " << domainId() << " zonal sending throws error " << endl;
        }
      }
    }
  }
}

void FvStructuredSolver3D::zonalReceive() {
  for(auto const& zBC : m_zonalBC) {
    for(MInt i = 0; i < zBC->noGlobalRcvDomains; i++) {
      // access only rcvDomains
      if(zBC->globalRcvDomainIds[i] == domainId()) {
        for(MInt noRcv = 0; noRcv < zBC->noRcvNghbrDomains; noRcv++) {
          MInt err = MPI_Irecv(&zBC->rcvComm[noRcv]->buffer[0], zBC->rcvComm[noRcv]->bufferSize * zBC->noZonalVariables,
                               MPI_DOUBLE, zBC->rcvComm[noRcv]->localId, 0, m_StructuredComm, &zBC->rcvRequest[noRcv],
                               AT_, "&zBC->bufferRcvZonal[noRcv][0]");
          if(err) cout << "rank " << domainId() << " zonal receiving throws error " << endl;
        }
      }
    }
    MPI_Waitall(zBC->noSndNghbrDomains, zBC->sndRequest, zBC->sndStatus, AT_);
    MPI_Waitall(zBC->noRcvNghbrDomains, zBC->rcvRequest, zBC->rcvStatus, AT_);
  }
}


void FvStructuredSolver3D::zonalScatter() {
  for(auto const& zBC : m_zonalBC) {
    for(MInt noRcv = 0; noRcv < zBC->noRcvNghbrDomains; noRcv++) {
      const MInt noCells = zBC->rcvComm[noRcv]->bufferSize;
      for(MInt i = 0; i < noCells; i++) {
        const MInt localMapCellIds = zBC->rcvComm[noRcv]->mapCellId[i];
        if(zBC->hasSTG) {
          m_cells->stg_fq[PV->U][localMapCellIds] = zBC->rcvComm[noRcv]->buffer[0 * noCells + i];
          m_cells->stg_fq[PV->V][localMapCellIds] = zBC->rcvComm[noRcv]->buffer[1 * noCells + i];
          m_cells->stg_fq[PV->W][localMapCellIds] = zBC->rcvComm[noRcv]->buffer[2 * noCells + i];
          m_cells->stg_fq[PV->RHO][localMapCellIds] = zBC->rcvComm[noRcv]->buffer[3 * noCells + i];
          m_cells->stg_fq[PV->P][localMapCellIds] = zBC->rcvComm[noRcv]->buffer[4 * noCells + i];
          m_cells->stg_fq[FQ->NU_T][localMapCellIds] = zBC->rcvComm[noRcv]->buffer[5 * noCells + i];
        } else {
          m_cells->fq[FQ->AVG_U][localMapCellIds] = zBC->rcvComm[noRcv]->buffer[0 * noCells + i];
          m_cells->fq[FQ->AVG_V][localMapCellIds] = zBC->rcvComm[noRcv]->buffer[1 * noCells + i];
          m_cells->fq[FQ->AVG_W][localMapCellIds] = zBC->rcvComm[noRcv]->buffer[2 * noCells + i];
          m_cells->fq[FQ->AVG_RHO][localMapCellIds] = zBC->rcvComm[noRcv]->buffer[3 * noCells + i];
          m_cells->fq[FQ->AVG_P][localMapCellIds] = zBC->rcvComm[noRcv]->buffer[4 * noCells + i];
          m_cells->fq[FQ->NU_T][localMapCellIds] = zBC->rcvComm[noRcv]->buffer[5 * noCells + i];
        }
      }
    }
  }
}


void FvStructuredSolver3D::initFluxMethod() {
  // choose the convective flux method
  // for performance issues, it is templated
  //(we do not need to write the muscl stuff too often
  // we need to put the limiter also into the template
  // in order to reduce this function!!
  if(m_rans == true) {
    reconstructSurfaceData = &FvStructuredSolver3D::MusclRANS;
    viscFluxMethod = &FvStructuredSolver3D::viscousFluxRANS;
  } else {
    viscFluxMethod = &FvStructuredSolver3D::viscousFluxLES<>;
    // check if limiter
    if(m_limiter) {
      switch(string2enum(m_limiterMethod)) {
        case VENKATAKRISHNAN_MOD: {
          /*! \property
            \page propertiesFVSTRCTRD
            \section venkFactor
            <code>MInt FvStructuredSolver::m_venkFactor </code>\n
            default = <code> 0 </code>\n \n
            Factor for the Venkatakrishnan Limiter.\n
            Possible values are:\n
            <ul>
            <li>Float > 0.0</li>
            </ul>
            Keywords: <i>LIMITER, STRUCTURED</i>
          */
          m_venkFactor = Context::getSolverProperty<MFloat>("venkFactor", m_solverId,
                                                            AT_); // reads the customizable parameter from properties

          m_log << "Using VENKATAKRISHNAN MOD limiter with VENK factor of " << m_venkFactor << " !" << endl;
          reconstructSurfaceData = &FvStructuredSolver3D::MusclVenkatakrishnan3D;
          Venkatakrishnan_function = &FvStructuredSolver3D::VENKATAKRISHNAN_MOD_FCT;
          break;
        }

        case VENKATAKRISHNAN: {
          m_log << "Using VENKATAKRISHNAN limiter!" << endl;
          reconstructSurfaceData = &FvStructuredSolver3D::MusclVenkatakrishnan3D;
          Venkatakrishnan_function = &FvStructuredSolver3D::VENKATAKRISHNAN_FCT;
          break;
        }

        case BARTH_JESPERSON: {
          m_log << "Using BARTH JESPERSON limiter!" << endl;
          reconstructSurfaceData = &FvStructuredSolver3D::MusclVenkatakrishnan3D;
          Venkatakrishnan_function = &FvStructuredSolver3D::BARTH_JESPERSON_FCT;
          break;
        }
        case MINMOD: {
          m_log << "Using MINMOD limiter!" << endl;
          reconstructSurfaceData = &FvStructuredSolver3D::MusclMinModLimiter;
          break;
        }
        case ALBADA: {
          m_log << "Using VAN ALBADA limiter!" << endl;
          reconstructSurfaceData = &FvStructuredSolver3D::MusclAlbada;
          break;
        }
        default: {
          stringstream errorMessage;
          errorMessage << "Limiter function " << m_limiterMethod << " not implemented!" << endl;
          mTerm(1, AT_, errorMessage.str());
        }
      }
    } else if(m_musclScheme == "Standard") {
      m_log << "Using unlimited MUSCL! (standard Formulation)" << endl;
      if(m_ausmScheme == "Standard") {
        m_log << "Using standard AUSM central" << endl;
        m_dsIsComputed = false;
        switch(CV->noVariables) {
          case 5: {
            reconstructSurfaceData = &FvStructuredSolver3D::Muscl_AusmLES<5>;
            break;
          }
          case 6: {
            reconstructSurfaceData = &FvStructuredSolver3D::Muscl_AusmLES<6>;
            break;
          }
          case 7: {
            reconstructSurfaceData = &FvStructuredSolver3D::Muscl_AusmLES<7>;
            break;
          }
          default: {
            stringstream errorMessage;
            errorMessage << "Number of Variables " << CV->noVariables << " not implemented in template AUSM!" << endl;
            mTerm(1, AT_);
          }
        }
      } else if(m_ausmScheme == "PTHRC") {
        m_log << "Using AUSM PTHRC" << endl;
        m_dsIsComputed = false;
        switch(CV->noVariables) {
          case 5: {
            reconstructSurfaceData = &FvStructuredSolver3D::Muscl_AusmLES_PTHRC<5>;
            break;
          }
          case 6: {
            reconstructSurfaceData = &FvStructuredSolver3D::Muscl_AusmLES_PTHRC<6>;
            break;
          }
          case 7: {
            reconstructSurfaceData = &FvStructuredSolver3D::Muscl_AusmLES_PTHRC<7>;
            break;
          }
          default: {
            stringstream errorMessage;
            errorMessage << "Number of Variables " << CV->noVariables << " not implemented in template AUSM!" << endl;
            mTerm(1, AT_);
          }
        }
      } else if(m_ausmScheme == "AUSMDV") {
        m_log << "Using AUSMDV" << endl;
        m_dsIsComputed = false;
        switch(CV->noVariables) {
          case 5: {
            reconstructSurfaceData = &FvStructuredSolver3D::Muscl_AusmDV<5>;
            break;
          }
          case 6: {
            reconstructSurfaceData = &FvStructuredSolver3D::Muscl_AusmDV<6>;
            break;
          }
          case 7: {
            reconstructSurfaceData = &FvStructuredSolver3D::Muscl_AusmDV<7>;
            break;
          }
          default: {
            stringstream errorMessage;
            errorMessage << "Number of Variables " << CV->noVariables << " not implemented in template AUSM!" << endl;
            mTerm(1, AT_);
          }
        }
      }
    } else if(m_musclScheme == "Stretched") {
      m_log << "Using unlimited MUSCL (streched Grids)";
      mAlloc(m_cells->cellLength, nDim, m_noCells, "m_cells->cellLength", -F1, AT_);
      computeCellLength();
      m_log << "Using standard AUSM central" << endl;
      if(m_ausmScheme == "Standard") {
        switch(CV->noVariables) {
          case 5: {
            reconstructSurfaceData = &FvStructuredSolver3D::MusclStretched_<&FvStructuredSolver3D::AusmLES, 5>;
            break;
          }
          case 6: {
            reconstructSurfaceData = &FvStructuredSolver3D::MusclStretched_<&FvStructuredSolver3D::AusmLES, 6>;
            break;
          }
          case 7: {
            reconstructSurfaceData = &FvStructuredSolver3D::MusclStretched_<&FvStructuredSolver3D::AusmLES, 7>;
            break;
          }
          default: {
            stringstream errorMessage;
            errorMessage << "Number of Variables " << CV->noVariables << " not implemented in template AUSM!" << endl;
            mTerm(1, AT_);
          }
        }
      } else if(m_ausmScheme == "PTHRC") {
        m_log << "Using AUSM PTHRC" << endl;
        switch(CV->noVariables) {
          case 5: {
            reconstructSurfaceData = &FvStructuredSolver3D::MusclStretched_<&FvStructuredSolver3D::AusmLES_PTHRC, 5>;
            break;
          }
          case 6: {
            reconstructSurfaceData = &FvStructuredSolver3D::MusclStretched_<&FvStructuredSolver3D::AusmLES_PTHRC, 6>;
            break;
          }
          case 7: {
            reconstructSurfaceData = &FvStructuredSolver3D::MusclStretched_<&FvStructuredSolver3D::AusmLES_PTHRC, 7>;
            break;
          }
          default: {
            stringstream errorMessage;
            errorMessage << "Number of Variables " << CV->noVariables << " not implemented in template AUSM!" << endl;
            mTerm(1, AT_);
          }
        }
      } else if(m_ausmScheme == "AUSMDV") {
        m_log << "Using AUSMDV" << endl;
        switch(CV->noVariables) {
          case 5: {
            reconstructSurfaceData = &FvStructuredSolver3D::MusclStretched_<&FvStructuredSolver3D::AusmDV, 5>;
            break;
          }
          case 6: {
            reconstructSurfaceData = &FvStructuredSolver3D::MusclStretched_<&FvStructuredSolver3D::AusmDV, 6>;
            break;
          }
          case 7: {
            reconstructSurfaceData = &FvStructuredSolver3D::MusclStretched_<&FvStructuredSolver3D::AusmDV, 7>;
            break;
          }
          default: {
            stringstream errorMessage;
            errorMessage << "Number of Variables " << CV->noVariables << " not implemented in template AUSM!" << endl;
            mTerm(1, AT_);
          }
        }
      }
    }
  }
}


void FvStructuredSolver3D::computeCellLength() {
  // this function can be moved into the MusclSchemeStreched later but for testing it is easier
  // REMEMBER: FOR MOVINg GRIDS THIS NEEDS TO BE CALLED EACH TIME

  MInt P1 = -1, P2 = -1, P3 = -1, P4 = -1, P5 = -1, P6 = -1, P7 = -1, P8 = -1;
  MFloat f1x = F0, f1y = F0, f1z = F0, f2x = F0, f2y = F0, f2z = F0;
  for(MInt k = 0; k < m_nCells[0]; k++) {
    for(MInt j = 0; j < m_nCells[1]; j++) {
      for(MInt i = 0; i < m_nCells[2]; i++) {
        const MInt cellId = cellIndex(i, j, k);
        P1 = getPointIdFromCell(i, j, k);
        P2 = getPointIdfromPoint(P1, 1, 0, 0);
        P3 = getPointIdfromPoint(P1, 1, 1, 0);
        P4 = getPointIdfromPoint(P1, 1, 0, 1);
        P5 = getPointIdfromPoint(P1, 1, 1, 1);
        P6 = getPointIdfromPoint(P1, 0, 1, 0);
        P7 = getPointIdfromPoint(P1, 0, 1, 1);
        P8 = getPointIdfromPoint(P1, 0, 0, 1);
        //----------Idirection
        // face 1
        f1x = F1B4
              * (m_grid->m_coordinates[0][P1] + m_grid->m_coordinates[0][P6] + m_grid->m_coordinates[0][P7]
                 + m_grid->m_coordinates[0][P8]);
        f1y = F1B4
              * (m_grid->m_coordinates[1][P1] + m_grid->m_coordinates[1][P6] + m_grid->m_coordinates[1][P7]
                 + m_grid->m_coordinates[1][P8]);
        f1z = F1B4
              * (m_grid->m_coordinates[2][P1] + m_grid->m_coordinates[2][P6] + m_grid->m_coordinates[2][P7]
                 + m_grid->m_coordinates[2][P8]);
        // face 2
        f2x = F1B4
              * (m_grid->m_coordinates[0][P2] + m_grid->m_coordinates[0][P3] + m_grid->m_coordinates[0][P4]
                 + m_grid->m_coordinates[0][P5]);
        f2y = F1B4
              * (m_grid->m_coordinates[1][P2] + m_grid->m_coordinates[1][P3] + m_grid->m_coordinates[1][P4]
                 + m_grid->m_coordinates[1][P5]);
        f2z = F1B4
              * (m_grid->m_coordinates[2][P2] + m_grid->m_coordinates[2][P3] + m_grid->m_coordinates[2][P4]
                 + m_grid->m_coordinates[2][P5]);
        m_cells->cellLength[0][cellId] = sqrt(POW2(f2x - f1x) + POW2(f2y - f1y) + POW2(f2z - f1z));
        //----------Jdirection
        // face 1
        f1x = F1B4
              * (m_grid->m_coordinates[0][P1] + m_grid->m_coordinates[0][P2] + m_grid->m_coordinates[0][P4]
                 + m_grid->m_coordinates[0][P8]);
        f1y = F1B4
              * (m_grid->m_coordinates[1][P1] + m_grid->m_coordinates[1][P2] + m_grid->m_coordinates[1][P4]
                 + m_grid->m_coordinates[1][P8]);
        f1z = F1B4
              * (m_grid->m_coordinates[2][P1] + m_grid->m_coordinates[2][P2] + m_grid->m_coordinates[2][P4]
                 + m_grid->m_coordinates[2][P8]);
        // face 2
        f2x = F1B4
              * (m_grid->m_coordinates[0][P3] + m_grid->m_coordinates[0][P5] + m_grid->m_coordinates[0][P6]
                 + m_grid->m_coordinates[0][P7]);
        f2y = F1B4
              * (m_grid->m_coordinates[1][P3] + m_grid->m_coordinates[1][P5] + m_grid->m_coordinates[1][P6]
                 + m_grid->m_coordinates[1][P7]);
        f2z = F1B4
              * (m_grid->m_coordinates[2][P3] + m_grid->m_coordinates[2][P5] + m_grid->m_coordinates[2][P6]
                 + m_grid->m_coordinates[2][P7]);
        m_cells->cellLength[1][cellId] = sqrt(POW2(f2x - f1x) + POW2(f2y - f1y) + POW2(f2z - f1z));
        //----------Kdirection
        // face 1
        f1x = F1B4
              * (m_grid->m_coordinates[0][P1] + m_grid->m_coordinates[0][P2] + m_grid->m_coordinates[0][P3]
                 + m_grid->m_coordinates[0][P6]);
        f1y = F1B4
              * (m_grid->m_coordinates[1][P1] + m_grid->m_coordinates[1][P2] + m_grid->m_coordinates[1][P3]
                 + m_grid->m_coordinates[1][P6]);
        f1z = F1B4
              * (m_grid->m_coordinates[2][P1] + m_grid->m_coordinates[2][P2] + m_grid->m_coordinates[2][P3]
                 + m_grid->m_coordinates[2][P6]);
        // face 2
        f2x = F1B4
              * (m_grid->m_coordinates[0][P4] + m_grid->m_coordinates[0][P5] + m_grid->m_coordinates[0][P7]
                 + m_grid->m_coordinates[0][P8]);
        f2y = F1B4
              * (m_grid->m_coordinates[1][P4] + m_grid->m_coordinates[1][P5] + m_grid->m_coordinates[1][P7]
                 + m_grid->m_coordinates[1][P8]);
        f2z = F1B4
              * (m_grid->m_coordinates[2][P5] + m_grid->m_coordinates[2][P5] + m_grid->m_coordinates[2][P7]
                 + m_grid->m_coordinates[2][P8]);
        m_cells->cellLength[2][cellId] = sqrt(POW2(f2x - f1x) + POW2(f2y - f1y) + POW2(f2z - f1z));
      }
    }
  }
}

MFloat FvStructuredSolver3D::getCellLengthY(MInt i, MInt j, MInt k) {
  const MInt P1 = getPointIdFromCell(i, j, k);
  const MInt P2 = getPointIdfromPoint(P1, 1, 0, 0);
  const MInt P3 = getPointIdfromPoint(P1, 1, 1, 0);
  const MInt P4 = getPointIdfromPoint(P1, 1, 0, 1);
  const MInt P5 = getPointIdfromPoint(P1, 1, 1, 1);
  const MInt P6 = getPointIdfromPoint(P1, 0, 1, 0);
  const MInt P7 = getPointIdfromPoint(P1, 0, 1, 1);
  const MInt P8 = getPointIdfromPoint(P1, 0, 0, 1);

  const MFloat lowerY = F1B4
                        * (m_grid->m_coordinates[1][P1] + m_grid->m_coordinates[1][P2] + m_grid->m_coordinates[1][P4]
                           + m_grid->m_coordinates[1][P8]);
  const MFloat upperY = F1B4
                        * (m_grid->m_coordinates[1][P3] + m_grid->m_coordinates[1][P5] + m_grid->m_coordinates[1][P6]
                           + m_grid->m_coordinates[1][P7]);

  return upperY - lowerY;
}

MFloat FvStructuredSolver3D::getCellCoordinate(MInt cellId, MInt dim) { return m_cells->coordinates[dim][cellId]; }


void FvStructuredSolver3D::initSolutionStep(MInt mode) {
  TRACE();

  std::ignore = mode;

  // Compute infinity values from property file
  // and (if no restart) fill cells according
  // to the initialCondition property
  initialCondition();

  if(m_restartInterpolation) {
    interpolateFromDonor();
  } else if(m_restart) {
    loadRestartFile();
  }

  // timestep will be computed and
  // set if globalTimeStep == 0 or
  // constantTimeStep = 0
  setTimeStep();

  // initialize moving grid
  // functions and move grid
  // to correct position
  if(m_movingGrid) {
    RECORD_TIMER_START(m_timers[Timers::MovingGrid]);
    RECORD_TIMER_START(m_timers[Timers::MGMoveGrid]);
    initMovingGrid();
    RECORD_TIMER_STOP(m_timers[Timers::MGMoveGrid]);
    RECORD_TIMER_STOP(m_timers[Timers::MovingGrid]);
  }

  if(m_bodyForce) {
    initBodyForce();
  }

  if(m_useSandpaperTrip) {
    initSandpaperTrip();
  }

  // Get the correct values
  // in the exchange ghostcells
  exchange();

  if(m_rans) {
    m_ransSolver->computeTurbViscosity();
  }

  if(m_zonal) {
    computeCumulativeAverage(false);
    spanwiseAvgZonal(m_zonalSpanwiseAvgVars);
    zonalExchange();
  }

  // Call the init function of each BC
  initBndryCnds();

  // Apply boundary conditions
  // and fill the non-exchange ghostcells
  applyBoundaryCondition();

  if(m_rans) {
    m_ransSolver->computeTurbViscosity();
  }

  // Check for NaNs
  checkNans();

  computeConservativeVariables();
}

/**
 * \brief Interpolates the flow field from a given donor file
 *
 * Instead of a restart from a restart file
 * this method will interpolate from a donor solution
 * onto the flow field
 *
 * \author Marian Albers
 */
void FvStructuredSolver3D::interpolateFromDonor() {
  m_structuredInterpolation = make_unique<StructuredInterpolation<3>>(m_StructuredComm);
  m_structuredInterpolation->prepareInterpolationField(m_nCells, m_cells->coordinates);

  if(m_zonal) {
    for(MInt var = 0; var < m_maxNoVariables; var++) {
      m_structuredInterpolation->interpolateField(m_pvariableNames[var], m_cells->pvariables[var]);
    }

    m_structuredInterpolation->interpolateField(FQ->fqNames[FQ->NU_T], m_cells->fq[FQ->NU_T]);
  } else {
    MBool donorConservative = false;
    if(Context::propertyExists("donorConservative", m_solverId)) {
      donorConservative = Context::getSolverProperty<MBool>("donorConservative", false, AT_);
    }

    if(donorConservative) {
      for(MInt var = 0; var < CV->noVariables; var++) {
        m_structuredInterpolation->interpolateField(m_variableNames[var], m_cells->variables[var]);
      }
      computePrimitiveVariables();
    } else {
      for(MInt var = 0; var < PV->noVariables; var++) {
        m_structuredInterpolation->interpolateField(m_pvariableNames[var], m_cells->pvariables[var]);
      }
    }

    // interpolate NU_T for STG if this is an initial start (otherwise not needed)
    if(m_stgIsActive && m_stgInitialStartup) {
      m_structuredInterpolation->interpolateField(FQ->fqNames[FQ->NU_T], m_cells->fq[FQ->NU_T]);
    }
  }

  if(m_useSponge && m_spongeLayerType == 2) {
    m_structuredInterpolation->interpolateField("rho", m_cells->fq[FQ->SPONGE_RHO]);
    m_structuredInterpolation->interpolateField("rhoE", m_cells->fq[FQ->SPONGE_RHO_E]);
  }

  if(m_useSponge && m_spongeLayerType == 4) {
    m_structuredInterpolation->interpolateField("rho", m_cells->fq[FQ->SPONGE_RHO]);
  }

  manualInterpolationCorrection();
}


/**
 * \brief Computation of infinity values for the conservative and primitive variables
 *
 * Initialization ot the entire flow field
 *
 * \author Pascal S. Meysonnat
 * \date 16.02.2011
 */
void FvStructuredSolver3D::initialCondition() {
  TRACE();
  const MFloat gammaMinusOne = m_gamma - 1.0;
  MFloat UT;
  MFloat pressureCH = F0;

  PV->TInfinity = 1.0 / (1.0 + F1B2 * gammaMinusOne * POW2(m_Ma));
  UT = m_Ma * sqrt(PV->TInfinity);
  PV->UInfinity = UT * cos(m_angle[0]) * cos(m_angle[1]);
  PV->VInfinity = UT * sin(m_angle[0]) * cos(m_angle[1]);
  PV->WInfinity = UT * sin(m_angle[1]);
  PV->VVInfinity[0] = PV->UInfinity;
  PV->VVInfinity[1] = PV->VInfinity;
  PV->VVInfinity[2] = PV->WInfinity;
  PV->PInfinity = pow(PV->TInfinity, (m_gamma / gammaMinusOne)) / m_gamma;

  // compute conservative variables
  CV->rhoInfinity = pow(PV->TInfinity, (1.0 / gammaMinusOne));
  CV->rhoUInfinity = CV->rhoInfinity * PV->UInfinity;
  CV->rhoVInfinity = CV->rhoInfinity * PV->VInfinity;
  CV->rhoWInfinity = CV->rhoInfinity * PV->WInfinity;
  CV->rhoVVInfinity[0] = CV->rhoUInfinity;
  CV->rhoVVInfinity[1] = CV->rhoVInfinity;
  CV->rhoVVInfinity[2] = CV->rhoWInfinity;
  CV->rhoEInfinity = PV->PInfinity / gammaMinusOne + CV->rhoInfinity * (F1B2 * POW2(UT));

  // internal Reynolds number Re0 = Re / ( rho8*M*sqrt(T8)/T8^F072)
  m_Re0 = m_Re * SUTHERLANDLAW(PV->TInfinity) / (CV->rhoInfinity * m_Ma * sqrt(PV->TInfinity));

  // reference enthalpies (needed for combustion computations)
  m_hInfinity = PV->PInfinity / CV->rhoInfinity * m_gamma / gammaMinusOne;

  // reference time (convection time)
  m_timeRef = UT / m_referenceLength;

  m_deltaP = F0;
  // pressure loss per unit length dp = rho_00 u_tau^2 L / D ) here: D=1.0, L=1;
  // channel: dp = rho_00 u_tau^2 L / D )
  // m_deltaP = POW2( m_Ma * m_ReTau * sqrt(PV->TInfinity) / m_Re  ) * CV->rhoInfinity / m_referenceLength;x
  // result is obtained by making deltap dimensionless with a_0^2 and rho_0
  // pipe: dp = lambda * L/D * rho/2 * u^2, lambda = 0.3164 Re^(-1/4) (Blasius)

  if(m_rans) {
    const MFloat lamVisc = SUTHERLANDLAW(PV->TInfinity);
    const MFloat chi = 0.1;
    CV->ransInfinity[0] = chi * (lamVisc);
    PV->ransInfinity[0] = chi * (lamVisc / CV->rhoInfinity);
  }

  m_log << "=================================================" << endl;
  m_log << "           INITIAL CONDITION SUMMARY" << endl;
  m_log << "=================================================" << endl;
  m_log << "Re = " << m_Re << endl;
  m_log << "Re0 = " << m_Re0 << endl;
  m_log << "Ma = " << m_Ma << endl;
  m_log << "TInfinity = " << PV->TInfinity << endl;
  m_log << "UInfinity = " << PV->UInfinity << endl;
  m_log << "VInfinity = " << PV->VInfinity << endl;
  m_log << "WInfinity = " << PV->WInfinity << endl;
  m_log << "PInfinity = " << PV->PInfinity << endl;
  m_log << "rhoInfinity = " << CV->rhoInfinity << endl;
  m_log << "rhoEInfinity = " << CV->rhoEInfinity << endl;
  m_log << "referenceTime = " << m_timeRef << endl;

  if(domainId() == 0) {
    cout << "////////////////////////////////////////////////" << endl;
    cout << "////////// Initial Condition summary ///////////" << endl;
    cout << "////////////////////////////////////////////////" << endl;
    cout << "Re = " << m_Re << endl;
    cout << "Re0 = " << m_Re0 << endl;
    cout << "Ma = " << m_Ma << endl;
    cout << "TInfinity = " << PV->TInfinity << endl;
    cout << "UInfinity = " << PV->UInfinity << endl;
    cout << "VInfinity = " << PV->VInfinity << endl;
    cout << "WInfinity = " << PV->WInfinity << endl;
    cout << "Angle = " << m_angle[0] << "  " << m_angle[1] << endl;
    cout << "PInfinity = " << PV->PInfinity << endl;
    cout << "rhoInfinity = " << CV->rhoInfinity << endl;
    cout << "rhoEInfinity = " << CV->rhoEInfinity << endl;
    cout << "referenceTime = " << m_timeRef << endl;
    cout << " zonal = " << m_zonal << endl;
  }

  if(!m_restart) {
    // inflow condition
    // ----------------
    switch(m_initialCondition) {
      case 0: {
        // parallel inflow field
        for(MInt cellid = 0; cellid < m_noCells; cellid++) {
          // go through every cell
          m_cells->pvariables[PV->RHO][cellid] = CV->rhoInfinity;
          for(MInt i = 0; i < nDim; i++) {
            m_cells->pvariables[PV->VV[i]][cellid] = PV->VVInfinity[i];
          }

          m_cells->pvariables[PV->P][cellid] = PV->PInfinity;

          if(m_rans) {
            m_cells->pvariables[PV->RANS_VAR[0]][cellid] = PV->ransInfinity[0];
          }
        }

        break;
      }
      case 43: {
        // parallel inflow field with pressure peak in the middle of the domain
        for(MInt cellid = 0; cellid < m_noCells; cellid++) {
          // go through every cell
          m_cells->pvariables[PV->RHO][cellid] = CV->rhoInfinity;
          for(MInt i = 0; i < nDim; i++) {
            m_cells->pvariables[PV->VV[i]][cellid] = F0;
          }

          m_cells->pvariables[PV->P][cellid] = PV->PInfinity;

          MFloat radius =
              sqrt(POW2(m_cells->coordinates[0][cellid] - 0.5) + POW2(m_cells->coordinates[1][cellid] - 0.5));

          // impose pressure peak in the middle of the domain
          if(radius <= 0.05) {
            MFloat pAmp = 0.005;
            MFloat pressureSignal = sin(radius / 0.05 * PI) * pAmp + PV->PInfinity;
            m_cells->pvariables[PV->P][cellid] = pressureSignal;
          }
        }
        break;
      }
      case 333: {
        // parallel inflow field
        for(MInt cellid = 0; cellid < m_noCells; cellid++) {
          // go through every cell
          m_cells->pvariables[PV->RHO][cellid] = CV->rhoInfinity;
          for(MInt i = 0; i < nDim; i++) {
            m_cells->pvariables[PV->VV[i]][cellid] = PV->VVInfinity[i];
          }

          m_cells->pvariables[PV->P][cellid] = PV->PInfinity;

          // impose pressure peak in the middle of the domain
          if(m_cells->coordinates[0][cellid] > 0.4 && m_cells->coordinates[0][cellid] < 0.5) {
            MFloat pAmp = 0.005;
            MFloat xCoordinate = m_cells->coordinates[0][cellid] - 0.4;
            MFloat pressureSignal = sin(xCoordinate / 0.1 * PI) * pAmp + PV->PInfinity;
            m_cells->pvariables[PV->P][cellid] = pressureSignal;
          }
        }
        break;
      }
      case 314: // stagnating flow field
      {
        for(MInt cellid = 0; cellid < m_noCells; cellid++) {
          m_cells->pvariables[PV->RHO][cellid] = CV->rhoInfinity;
          for(MInt i = 0; i < nDim; i++) {
            m_cells->pvariables[PV->VV[i]][cellid] = F0;
          }

          m_cells->pvariables[PV->P][cellid] = PV->PInfinity;
        }
        cout << "I.C. stagnating flow field was applied! " << endl;
        break;
      }
      case 315: // Poiseuille flow
      {
        MFloat x = F0, y = F0; // p=F0, T=F0, T0=F0;
        MFloat y_max = F1;     // channel height

        for(MInt cellid = 0; cellid < m_noCells; cellid++) {
          x = m_cells->coordinates[0][cellid];
          y = m_cells->coordinates[1][cellid];

          for(MInt i = 0; i < nDim; i++) {
            m_cells->pvariables[PV->VV[i]][cellid] = F0;
          }
          // all velocity components are 0 except u, Poiseuille distribution:
          m_cells->pvariables[PV->VV[0]][cellid] =
              -(F3 / F2) * PV->UInfinity * (POW2(y - y_max / F2) - POW2(y_max / F2)) / POW2(y_max / F2);

          // the pressure is defined through the axial pressure gradient:
          m_cells->pvariables[PV->P][cellid] =
              PV->PInfinity - F3 * (x + 15.0) * SUTHERLANDLAW(PV->TInfinity) * PV->UInfinity * POW2(F2 / y_max) / m_Re0;

          // compressible Poiseuille Temperature distribution
          // T = T0 - (m_gamma-F1)*m_Pr * POW2(PV->UInfinity*F3/F2) * (F1B2*(F1+ pow( (y-y_max/F2)/(y_max/F2), F4 )
          // )-POW2((y-y_max/F2)/(y_max/F2)));
          m_cells->pvariables[PV->RHO][cellid] = CV->rhoInfinity; // m_gamma * m_cells->variables[ PV->P ][cellid] /
          // T;//
        }

        break;
      }
      case 101: // TAYLOR_GREEN_VORTEX
      {
        // domain boundaries are all 2*pi
        // rho=1.0;
        // u=A*SIN(x)*COS(y)*COS(z)
        // v=- A*COS(x)*SIN(y)*COS(z)
        // w=0.0
        // p=A*A*rho*( 1./(Ms*Ms*kappa) + 1./16.*(COS(2*x)*COS(2.*z)+ 2.*COS(2.*y) +2.*COS(2.*x) +COS(2*y)*COS(2.*z)))
        // Ms =0.1 maximum Mach number
        // A= speed magnitude set to 1.0
        MInt cellId = 0;
        MFloat A = PV->UInfinity;
        MFloat x = F0;
        MFloat y = F0;
        MFloat z = F0;
        for(MInt k = 0; k < m_nCells[0]; k++) {
          for(MInt j = 0; j < m_nCells[1]; j++) {
            for(MInt i = 0; i < m_nCells[2]; i++) {
              cellId = cellIndex(i, j, k);
              x = m_cells->coordinates[0][cellId];
              y = m_cells->coordinates[1][cellId];
              z = m_cells->coordinates[2][cellId];
              m_cells->pvariables[PV->RHO][cellId] = CV->rhoInfinity;
              m_cells->pvariables[PV->VV[0]][cellId] = A * sin(x) * cos(y) * cos(z);
              m_cells->pvariables[PV->VV[1]][cellId] = -A * cos(x) * sin(y) * cos(z);
              m_cells->pvariables[PV->VV[2]][cellId] = 0.0;
              m_cells->pvariables[PV->P][cellId] =
                  PV->PInfinity
                  + F1B16 * (POW2(A) * CV->rhoInfinity) * (cos(2.0 * x) + cos(2.0 * y)) * (2.0 + cos(2.0 * z));
            }
          }
        }
        break;
      }
      case 1234: {
        // laminar channel flow
        MInt cellId = 0;
        m_channelHeight = 2.0;
        m_channelLength = 6.2831;
        // m_deltaP=32.0*SUTHERLANDLAW(PV->TInfinity)*m_Ma*sqrt(PV->TInfinity)*m_channelLength/(m_Re0*m_channelHeight);
        m_deltaP =
            -12.0 * PV->UInfinity * SUTHERLANDLAW(PV->TInfinity) * m_channelLength / (POW2(m_channelHeight) * m_Re0);

        m_channelPresInlet = PV->PInfinity;
        m_channelPresOutlet = PV->PInfinity + m_deltaP;
        MFloat u =
            PV->UInfinity; // POW2(m_channelHeight)*m_deltaP*m_Re0*m_referenceLength*m_referenceLength/m_channelLength;
        for(MInt k = m_noGhostLayers; k < m_nCells[0] - m_noGhostLayers; k++) {
          for(MInt j = m_noGhostLayers; j < m_nCells[1] - m_noGhostLayers; j++) {
            for(MInt i = m_noGhostLayers; i < m_nCells[2] - m_noGhostLayers; i++) {
              cellId = cellIndex(i, j, k);

              // channel height is in j direction
              // so prescribe a mean profile in u(y)
              MFloat x = m_cells->coordinates[0][cellId];
              pressureCH = m_channelPresInlet + (m_channelPresOutlet - m_channelPresInlet) * x / m_channelLength;
              // density:
              MFloat rho = pressureCH * m_gamma / PV->TInfinity;
              m_cells->pvariables[PV->RHO][cellId] = rho;
              m_cells->pvariables[PV->U][cellId] = u;
              m_cells->pvariables[PV->V][cellId] = F0;
              m_cells->pvariables[PV->W][cellId] = F0;
              m_cells->pvariables[PV->P][cellId] = pressureCH;
            }
          }
        }

        break;
      }
      case 1233: {
        // turbulent channel with perturbations
        // calculate the Pressure loss;

        // for the law of the wall
        const MFloat C1 = m_channelC1;
        const MFloat C2 = m_channelC2;
        const MFloat C3 = m_channelC3;
        const MFloat C4 = m_channelC4;
        MFloat xINIT = m_channelInflowPlaneCoordinate;
        MInt cellId = 0;
        MFloat yplus = F0;
        MFloat uTau = m_ReTau * m_Ma * sqrt(PV->TInfinity) / m_Re;
        MFloat prefactor = m_ReTau; // uTau*m_Re0/SUTHERLANDLAW(PV->TInfinity);

        m_deltaP = -CV->rhoInfinity * POW2(uTau) * F2 * (m_channelLength) / m_channelHeight;

        m_log << "uTau: " << uTau << " channelLength: " << m_channelLength << endl;
        // mean velocity profile
        m_channelPresInlet = PV->PInfinity;
        m_channelPresOutlet =
            PV->PInfinity - CV->rhoInfinity * POW2(uTau) * F2 * (xINIT + m_channelLength) / m_channelHeight;
        MFloat deltaP = m_channelPresInlet - m_channelPresOutlet;
        MFloat cfTheo = 2 * POW2(uTau / PV->UInfinity);
        MFloat cdTheo = cfTheo * 2.0 * m_channelWidth * m_channelLength;
        m_log << "deltaP: " << deltaP << " cfTheoretisch: " << cfTheo << " cdTheo: " << cdTheo << endl;


        for(MInt k = 0; k < m_nCells[0]; k++) {
          for(MInt j = 0; j < m_nCells[1]; j++) {
            for(MInt i = 0; i < m_nCells[2]; i++) {
              // channel height is in j direction
              // so prescribe a mean profile in u(y)
              cellId = cellIndex(i, j, k);
              // channel starts at x=0 or we need to prescribe an offset
              MFloat y = m_cells->coordinates[1][cellId];
              MFloat x = m_cells->coordinates[0][cellId];
              pressureCH =
                  m_channelPresInlet + ((m_channelPresOutlet - m_channelPresInlet) * (x - xINIT) / m_channelLength);

              MFloat rho = pressureCH * m_gamma / PV->TInfinity;
              MFloat velFactor = 2.0;
              if(y < m_channelHeight / 2) {
                yplus = prefactor * y;
              } else {
                yplus = prefactor * (m_channelHeight - y);
              }

              // C1 etc are defined in maiaconstants.h
              if(yplus <= 5.0) {
                m_cells->pvariables[PV->U][cellId] = 0.5 * uTau * yplus * velFactor;
              } else if(yplus <= 30 && yplus > 5.0) {
                m_cells->pvariables[PV->U][cellId] = 0.5 * uTau * (C1 * log(yplus) + C2) * velFactor;
              } else if(yplus > 30) {
                m_cells->pvariables[PV->U][cellId] = 0.5 * uTau * (C3 * log(yplus) + C4) * velFactor;
              }

              m_cells->pvariables[PV->RHO][cellId] = rho;
              m_cells->pvariables[PV->V][cellId] = F0;
              m_cells->pvariables[PV->W][cellId] = F0;
              m_cells->pvariables[PV->P][cellId] = pressureCH;
            }
          }
        }
        // create the fluctuations:


        // this way of creating the
        // fluctuations is only possible for
        // quite small mounts of cells
        // if too large this approach will not work anymore
        // 125000000 is a random integer number which has to be tested
        // else the approach with white noise will be employed
        const MInt totalBlockCells =
            (m_grid->getMyBlockNoCells(0) * m_grid->getMyBlockNoCells(1) * m_grid->getMyBlockNoCells(2));

        if(totalBlockCells <= 12500000) {
          fftw_complex *uPhysField, *vPhysField, *wPhysField;

          // we need the total number of points
          MInt lx = m_grid->getMyBlockNoCells(2), ly = m_grid->getMyBlockNoCells(1), lz = m_grid->getMyBlockNoCells(0);

          // field of velocities from positve frequencies
          uPhysField = (fftw_complex*)fftw_malloc(lx * ly * lz * sizeof(fftw_complex));
          vPhysField = (fftw_complex*)fftw_malloc(lx * ly * lz * sizeof(fftw_complex));
          wPhysField = (fftw_complex*)fftw_malloc(lx * ly * lz * sizeof(fftw_complex));

          // no real parallel computation is done
          // should be implemented later !!!!

          if(noDomains() > 1) {
            MFloatScratchSpace sendRecvBufferU(lx * ly * lz, AT_, "sendRecvBufferU");
            MFloatScratchSpace sendRecvBufferV(lx * ly * lz, AT_, "sendRecvBufferV");
            MFloatScratchSpace sendRecvBufferW(lx * ly * lz, AT_, "sendRecvBufferW");
            if(domainId() == 0) {
              MInt m_noPeakModes = 100;
              initFFTW(uPhysField, vPhysField, wPhysField, lx, ly, lz, m_noPeakModes);
              // copy values into the sendRCVbuffer

              for(MInt id = 0; id < lz * ly * lz; id++) {
                sendRecvBufferU[id] = uPhysField[id][0];
                sendRecvBufferV[id] = vPhysField[id][0];
                sendRecvBufferW[id] = wPhysField[id][0];
              }
            }
            MPI_Bcast(&sendRecvBufferU[0], lx * ly * lz, MPI_DOUBLE, 0, m_StructuredComm, AT_, "sendRecvBufferU[0]");
            MPI_Bcast(&sendRecvBufferV[0], lx * ly * lz, MPI_DOUBLE, 0, m_StructuredComm, AT_, "sendRecvBufferV[0]");
            MPI_Bcast(&sendRecvBufferW[0], lx * ly * lz, MPI_DOUBLE, 0, m_StructuredComm, AT_, "sendRecvBufferW[0]");

            if(domainId() != 0) {
              for(MInt id = 0; id < lz * ly * lz; id++) {
                uPhysField[id][0] = sendRecvBufferU[id];
                vPhysField[id][0] = sendRecvBufferV[id];
                wPhysField[id][0] = sendRecvBufferW[id];
              }
            }
          } else {
            // create velocity field
            MInt m_noPeakModes = 100;
            initFFTW(uPhysField, vPhysField, wPhysField, lx, ly, lz, m_noPeakModes);
          }
          // now we need to distribute the
        } else {
          MFloat amp = 0.15;
          for(MInt k = m_noGhostLayers; k < m_nCells[0] - m_noGhostLayers; k++) {
            for(MInt j = m_noGhostLayers; j < m_nCells[1] - m_noGhostLayers; j++) {
              for(MInt i = m_noGhostLayers; i < m_nCells[2] - m_noGhostLayers; i++) {
                cellId = cellIndex(i, j, k);
                m_cells->pvariables[PV->V][cellId] +=
                    amp * 2.0 * (0.5 - (1.0 * rand() / (RAND_MAX + 1.0))) * m_cells->pvariables[PV->U][cellId];
                m_cells->pvariables[PV->W][cellId] +=
                    amp * 2.0 * (0.5 - (1.0 * rand() / (RAND_MAX + 1.0))) * m_cells->pvariables[PV->U][cellId];
                m_cells->pvariables[PV->U][cellId] +=
                    amp * 2.0 * (0.5 - (1.0 * rand() / (RAND_MAX + 1.0))) * m_cells->pvariables[PV->U][cellId];
              }
            }
          }
        }

        break;
      }
      case 1236: {
        // pipe with perturbations
        // calculate the Pressure loss;
        // for the law of the wall
        // const MFloat C1      = m_channelC1;
        // const MFloat C2      = m_channelC2;
        const MFloat C3 = m_channelC3;
        const MFloat C4 = m_channelC4;
        MInt cellId = 0;
        MFloat uTau = m_ReTau * m_Ma * sqrt(PV->TInfinity) / m_Re;
        // MFloat prefactor=m_ReTau;//uTau*m_Re0/SUTHERLANDLAW(PV->TInfinity);

        m_deltaP = -4.0 * CV->rhoInfinity * POW2(uTau) * (m_channelLength) / m_channelHeight;
        m_channelPresInlet = PV->PInfinity;
        m_channelPresOutlet = PV->PInfinity + m_deltaP;
        const MFloat bulkVel = uTau * (C3 * log(m_ReTau / 2) + C4);

        m_log << "=========== Turb. Pipe Flow Inital Condition Summary =========== " << endl;
        m_log << "-->Turbulent pipe flow deltaP: " << m_deltaP << endl;
        m_log << "-->pipe friciton velocity: " << uTau << endl;
        m_log << "-->pipe pressure inflow: " << m_channelPresInlet << endl;
        m_log << "-->pipe pressure outflow: " << m_channelPresOutlet << endl;
        m_log << "--> bulk velocity (u_max)" << bulkVel << endl;
        m_log << "=========== Turb. Pipe Flow Initial Condition Summary Finished =========== " << endl;

        for(MInt k = 0; k < m_nCells[0]; k++) {
          for(MInt j = 0; j < m_nCells[1]; j++) {
            for(MInt i = 0; i < m_nCells[2]; i++) {
              // IMPORTANT PIPE LENGTH HAS TO GO INTO X-DIRECTION
              // so prescribe a mean profile in u(r)
              // centerline is assumed at y=0.0,z=0.0

              cellId = cellIndex(i, j, k);
              // determine the radius
              MFloat r = sqrt(POW2(m_cells->coordinates[1][cellId]) + POW2(m_cells->coordinates[2][cellId]));

              MFloat x = m_cells->coordinates[0][cellId];
              // determine the pressure drop
              pressureCH = m_deltaP / m_channelLength * (x - m_channelInflowPlaneCoordinate) + PV->PInfinity;
              MFloat rho = pressureCH * m_gamma / PV->TInfinity;
              // important: viscous sublayer is not build explicitly. The flow has to build it by itself
              MFloat vel = bulkVel + C3 * log(F1 - min(((F2 * r) / m_channelHeight), 0.999999999)) * uTau;
              m_cells->pvariables[PV->RHO][cellId] = rho;
              m_cells->pvariables[PV->U][cellId] = vel;
              m_cells->pvariables[PV->V][cellId] = F0;
              m_cells->pvariables[PV->W][cellId] = F0;
              m_cells->pvariables[PV->P][cellId] = pressureCH;
            }
          }
        }

        // create the fluctuations:
        // this way of creating the
        // fluctuations is only possible for
        // quite small mounts of cells
        // if too large this approach will not work anymore
        // 125000000 is a random integer number which has to be tested
        // else the approach with white noise will be employed
        const MInt totalBlockCells =
            (m_grid->getMyBlockNoCells(0) * m_grid->getMyBlockNoCells(1) * m_grid->getMyBlockNoCells(2));

        if(totalBlockCells <= 12500000) {
          fftw_complex *uPhysField, *vPhysField, *wPhysField;

          // we need the total number of points
          MInt lx = m_grid->getMyBlockNoCells(2), ly = m_grid->getMyBlockNoCells(1), lz = m_grid->getMyBlockNoCells(0);

          // field of velocities from positve frequencies
          uPhysField = (fftw_complex*)fftw_malloc(lx * ly * lz * sizeof(fftw_complex));
          vPhysField = (fftw_complex*)fftw_malloc(lx * ly * lz * sizeof(fftw_complex));
          wPhysField = (fftw_complex*)fftw_malloc(lx * ly * lz * sizeof(fftw_complex));

          // no real parallel computation is done
          // should be implemented later !!!!

          if(noDomains() > 1) {
            MFloatScratchSpace sendRecvBufferU(lx * ly * lz, AT_, "sendRecvBufferU");
            MFloatScratchSpace sendRecvBufferV(lx * ly * lz, AT_, "sendRecvBufferV");
            MFloatScratchSpace sendRecvBufferW(lx * ly * lz, AT_, "sendRecvBufferW");
            if(domainId() == 0) {
              MInt m_noPeakModes = 100;
              initFFTW(uPhysField, vPhysField, wPhysField, lx, ly, lz, m_noPeakModes);
              // copy values into the sendRCVbuffer

              for(MInt id = 0; id < lz * ly * lz; id++) {
                sendRecvBufferU[id] = uPhysField[id][0];
                sendRecvBufferV[id] = vPhysField[id][0];
                sendRecvBufferW[id] = wPhysField[id][0];
              }
            }
            MPI_Bcast(&sendRecvBufferU[0], lx * ly * lz, MPI_DOUBLE, 0, m_StructuredComm, AT_, "sendRecvBufferU[0]");
            MPI_Bcast(&sendRecvBufferV[0], lx * ly * lz, MPI_DOUBLE, 0, m_StructuredComm, AT_, "sendRecvBufferV[0]");
            MPI_Bcast(&sendRecvBufferW[0], lx * ly * lz, MPI_DOUBLE, 0, m_StructuredComm, AT_, "sendRecvBufferW[0]");

            if(domainId() != 0) {
              for(MInt id = 0; id < lz * ly * lz; id++) {
                uPhysField[id][0] = sendRecvBufferU[id];
                vPhysField[id][0] = sendRecvBufferV[id];
                wPhysField[id][0] = sendRecvBufferW[id];
              }
            }
          } else {
            // create velocity field
            MInt m_noPeakModes = 100;
            initFFTW(uPhysField, vPhysField, wPhysField, lx, ly, lz, m_noPeakModes);
          }
          // now we need to distribute the
        } else {
          MFloat amp = 0.15;
          for(MInt k = m_noGhostLayers; k < m_nCells[0] - m_noGhostLayers; k++) {
            for(MInt j = m_noGhostLayers; j < m_nCells[1] - m_noGhostLayers; j++) {
              for(MInt i = m_noGhostLayers; i < m_nCells[2] - m_noGhostLayers; i++) {
                cellId = cellIndex(i, j, k);
                m_cells->pvariables[PV->V][cellId] +=
                    amp * 2.0 * (0.5 - (1.0 * rand() / (RAND_MAX + 1.0))) * m_cells->pvariables[PV->U][cellId];
                m_cells->pvariables[PV->W][cellId] +=
                    amp * 2.0 * (0.5 - (1.0 * rand() / (RAND_MAX + 1.0))) * m_cells->pvariables[PV->U][cellId];
                m_cells->pvariables[PV->U][cellId] +=
                    amp * 2.0 * (0.5 - (1.0 * rand() / (RAND_MAX + 1.0))) * m_cells->pvariables[PV->U][cellId];
              }
            }
          }
        }

        break;
      }
      case 79092: {
        // approximate mean turbulent boundary layer profile
        const MFloat epss = 1e-10;
        const MFloat reTheta = 1000.000;
        const MFloat theta = 1.0;
        const MFloat delta0 = 72.0 / 7.0 * theta;
        const MFloat kappa = 0.4;
        const MFloat C1 = 3.573244189003983;
        const MFloat PI1 = 0.55;
        const MFloat cf = 0.024 / pow(reTheta, 0.25);
        const MFloat uTau = sqrt(cf / 2.0) * PV->UInfinity;
        const MFloat nu = SUTHERLANDLAW(PV->TInfinity);
        for(MInt k = 0; k < m_nCells[0]; k++) {
          for(MInt j = 0; j < m_nCells[1]; j++) {
            for(MInt i = 0; i < m_nCells[2]; i++) {
              const MInt cellId = cellIndex(i, j, k);
              m_cells->pvariables[PV->RHO][cellId] = CV->rhoInfinity;
              const MFloat y = m_cells->coordinates[1][cellId];
              MFloat yPlus = mMax(uTau * y * m_Re0 / nu, F0);

              if(y > delta0) {
                m_cells->pvariables[PV->U][cellId] = PV->UInfinity;
              } else if(yPlus <= 10.0) {
                m_cells->pvariables[PV->U][cellId] = yPlus * uTau;
              } else if(yPlus <= 30 && yPlus > 10.0) {
                m_cells->pvariables[PV->U][cellId] = uTau * ((F1 / kappa) * log(max(yPlus, epss)) + C1);
              } else if(yPlus > 30.0 && y <= delta0) {
                m_cells->pvariables[PV->U][cellId] =
                    uTau
                    * ((F1 / kappa) * log(max(yPlus, epss)) + C1 + 2 * (PI1 / kappa) * (3 * y * y - 2 * y * y * y));
              }

              m_cells->pvariables[PV->V][cellId] = PV->VInfinity;
              m_cells->pvariables[PV->W][cellId] = PV->WInfinity;
              m_cells->pvariables[PV->RHO][cellId] = CV->rhoInfinity;
              m_cells->pvariables[PV->P][cellId] = PV->PInfinity;

              if(m_rans) {
                m_cells->pvariables[PV->RANS_VAR[0]][cellId] = PV->ransInfinity[0];
              }
            }
          }
        }

        break;
      }

      case 79091: // Turbulent plate
      {
        const MFloat epss = 1e-10;
        const MFloat reTheta = 1000.0;
        const MFloat theta = 1.0;
        const MFloat delta0 = 72.0 / 7.0 * theta;
        const MFloat K = 0.4;
        const MFloat C1 = 3.573244189003983; // With coles
        const MFloat PI1 = 0.55;
        const MFloat cf = 0.024 / pow(reTheta, 0.25);

        for(MInt k = 0; k < m_nCells[0]; k++) {
          for(MInt j = 0; j < m_nCells[1]; j++) {
            for(MInt i = 0; i < m_nCells[2]; i++) {
              const MInt cellId = cellIndex(i, j, k);
              const MFloat mu = SUTHERLANDLAW(PV->TInfinity);
              const MFloat utau = sqrt(cf / 2.0) * m_Ma * sqrt(PV->TInfinity);
              const MFloat yplus = m_cells->coordinates[1][cellId] * sqrt(cf / 2.) * CV->rhoUInfinity / mu * m_Re0;
              const MFloat eta = m_cells->coordinates[1][cellId] / delta0; // y/delta

              // 1-7th profile
              // log-law + wake
              if(m_cells->coordinates[1][cellId] > delta0) {
                m_cells->pvariables[PV->U][cellId] = PV->UInfinity; // Outside BL
              } else if(yplus < 10) {
                m_cells->pvariables[PV->U][cellId] = utau * yplus;
              } else {
                m_cells->pvariables[PV->U][cellId] = mMin(
                    utau
                        * ((1. / K) * log(max(yplus, epss)) + C1 + 2 * PI1 / K * (3 * eta * eta - 2 * eta * eta * eta)),
                    PV->UInfinity);
              }

              m_cells->pvariables[PV->V][cellId] = PV->VInfinity;
              m_cells->pvariables[PV->W][cellId] = PV->WInfinity;
              m_cells->pvariables[PV->RHO][cellId] = CV->rhoInfinity;
              m_cells->pvariables[PV->P][cellId] = PV->PInfinity;

              if(m_rans) {
                m_cells->pvariables[PV->RANS_VAR[0]][cellId] = PV->ransInfinity[0];
              }
            }
          }
        }

        break;
      }
      case 11: // point source in the middle
      {
        // contain an initial perturbation in the middle
        for(MInt cellid = 0; cellid < m_noCells; cellid++) {
          // go through every cell
          m_cells->pvariables[PV->RHO][cellid] = CV->rhoInfinity;
          for(MInt i = 0; i < nDim; i++) {
            m_cells->pvariables[PV->VV[i]][cellid] = PV->VVInfinity[i];
          }

          MFloat amp = 0.0001;
          MFloat x = m_cells->coordinates[0][cellid];
          MFloat y = m_cells->coordinates[1][cellid];
          MFloat z = m_cells->coordinates[2][cellid];
          MFloat r = 0.025;
          MFloat a1 = POW2(x - 0.5);
          MFloat a2 = POW2(y - 0.5);
          MFloat a3 = POW2(z - 0.5);
          MFloat disturb = amp * exp(-(a1 + a2 + a3) / POW2(r) / 2.0);
          m_cells->pvariables[PV->P][cellid] = PV->PInfinity + disturb;
        }
        break;
      }
      /* TESTCASE from C. Bogey, C. Bailly, Three-dimensional non-reflective boundary conditions
       *  for acoustic simulations: far field formulation and validation test cases
       *  Acta acustica united with acustica Volume 88 (2002), 463-471
       */
      case 111: {
        MFloat amp = 0.01;
        MFloat alpha = log(2) / 9;
        for(MInt cellid = 0; cellid < m_noCells; cellid++) {
          // go through every cell
          m_cells->pvariables[PV->RHO][cellid] = CV->rhoInfinity;
          MFloat x = m_cells->coordinates[0][cellid];
          MFloat y = m_cells->coordinates[1][cellid];
          MFloat z = m_cells->coordinates[2][cellid];
          MFloat fluc = amp * exp(-alpha * (POW2(x) + POW2(y) + POW2(z)));
          m_cells->pvariables[PV->RHO][cellid] = CV->rhoInfinity + fluc;
          m_cells->pvariables[PV->P][cellid] = (PV->PInfinity + fluc);
          for(MInt i = 0; i < nDim; i++) {
            m_cells->pvariables[PV->VV[i]][cellid] = PV->VVInfinity[i];
          }
        }
        break;
      }
        /* TESTCASE from C. Bogey, C. Bailly, Three-dimensional non-reflective boundary conditions
         *  for acoustic simulations: far field formulation and validation test cases
         *  Acta acustica united with acustica Volume 88 (2002), 463-471 the vortex
         */
      case 112: {
        MFloat amp = 0.003;
        MFloat b = 0.025;
        MFloat alpha = log(2) / POW2(b);
        MFloat r0 = 0.1; // changed!!!
        MFloat r = F0;
        MFloat phi = F0;
        MFloat v = F0;
        for(MInt cellId = 0; cellId < m_noCells; cellId++) {
          // go through every cell
          m_cells->pvariables[PV->RHO][cellId] = CV->rhoInfinity;
          m_cells->pvariables[PV->P][cellId] = PV->PInfinity;
          m_cells->pvariables[PV->U][cellId] = PV->UInfinity;
          m_cells->pvariables[PV->V][cellId] = PV->VInfinity;
          m_cells->pvariables[PV->W][cellId] = PV->WInfinity;
          MFloat x = m_cells->coordinates[0][cellId] - 0.5;
          MFloat y = m_cells->coordinates[1][cellId] - 0.5;
          MFloat z = m_cells->coordinates[2][cellId] - 0.25;
          r = sqrt(POW2(y) + POW2(z));
          phi = atan(z / y);
          m_cells->pvariables[PV->U][cellId] +=
              amp * (r0 / r) * (r - r0) * exp(-1.0 * alpha * (POW2(x) + POW2(r - r0)));
          v = -1.0 * amp * (r0 / r) * x * exp(-1.0 * alpha * (POW2(x) + POW2(r - r0)));
          m_cells->pvariables[PV->V][cellId] += v * r * sin(phi);
          m_cells->pvariables[PV->W][cellId] += v * r * cos(phi);
        }
        break;
      }
      case 113: {
        MFloat amp = 0.003;
        MFloat b = 0.25 / 4;
        MFloat alpha = log(2) / POW2(b);
        MFloat r0 = 0.25; // changed!!!
        MFloat r = F0;
        MFloat phi = F0;
        MFloat v = F0;
        for(MInt cellId = 0; cellId < m_noCells; cellId++) {
          // go through every cell
          m_cells->pvariables[PV->RHO][cellId] = CV->rhoInfinity;
          m_cells->pvariables[PV->P][cellId] = PV->PInfinity;
          m_cells->pvariables[PV->U][cellId] = PV->UInfinity;
          m_cells->pvariables[PV->V][cellId] = PV->VInfinity;
          m_cells->pvariables[PV->W][cellId] = PV->WInfinity;
          MFloat x = m_cells->coordinates[0][cellId] + 0.5;
          MFloat y = m_cells->coordinates[1][cellId];
          MFloat z = m_cells->coordinates[2][cellId] - 0.25;
          r = sqrt(POW2(y) + POW2(z));
          phi = atan(z / y);
          m_cells->pvariables[PV->U][cellId] +=
              amp * (r0 / r) * (r - r0) * exp(-1.0 * alpha * (POW2(x) + POW2(r - r0)));
          v = -1.0 * amp * (r0 / r) * x * exp(-1.0 * alpha * (POW2(x) + POW2(r - r0)));
          m_cells->pvariables[PV->V][cellId] += v * r * sin(phi);
          m_cells->pvariables[PV->W][cellId] += v * r * cos(phi);
        }
        break;
      }
      case 2: // shear flow is prescribed
      {
        MInt cellId = 0;
        MFloat x = F0;
        for(MInt k = 0; k < m_nCells[0]; k++) {
          for(MInt j = 0; j < m_nCells[1]; j++) {
            for(MInt i = 0; i < m_nCells[2]; i++) {
              cellId = cellIndex(i, j, k);
              x = m_cells->coordinates[0][cellId];
              m_cells->pvariables[PV->RHO][cellId] = CV->rhoInfinity;
              if(x < 0.5) {
                m_cells->pvariables[PV->U][cellId] = F0;
                m_cells->pvariables[PV->V][cellId] = 0.5 * PV->UInfinity;
                m_cells->pvariables[PV->W][cellId] = F0;
              } else {
                m_cells->pvariables[PV->U][cellId] = F0;
                m_cells->pvariables[PV->V][cellId] = PV->UInfinity;
                m_cells->pvariables[PV->W][cellId] = F0;
              }
              m_cells->pvariables[PV->P][cellId] = PV->PInfinity;
            }
          }
        }
        break;
      }
      case 900: { // jet Freund
        for(MInt k = 0; k < m_nCells[0]; k++) {
          for(MInt j = 0; j < m_nCells[1]; j++) {
            for(MInt i = 0; i < m_nCells[2]; i++) {
              MInt cellId = cellIndex(i, j, k);
              MFloat r = sqrt(POW2(m_cells->coordinates[0][cellId]) + POW2(m_cells->coordinates[1][cellId])
                              + POW2(m_cells->coordinates[2][cellId]));
              MFloat u = F1B2 * (F1 - tanh(12.5 * (fabs(r / 0.5) - fabs(0.5 / r)))) * PV->VVInfinity[0];
              m_cells->pvariables[PV->RHO][cellId] = CV->rhoInfinity;
              m_cells->pvariables[PV->U][cellId] = u;
              m_cells->pvariables[PV->V][cellId] = PV->VInfinity;
              m_cells->pvariables[PV->W][cellId] = PV->WInfinity;
              m_cells->pvariables[PV->P][cellId] = PV->PInfinity;
            }
          }
        }
        break;
      }
      case 4001: { // test periodic rotation boundary conditions
        for(MInt k = m_noGhostLayers; k < m_nCells[0] - m_noGhostLayers; ++k) {
          for(MInt j = m_noGhostLayers; j < m_nCells[1] - m_noGhostLayers; ++j) {
            for(MInt i = m_noGhostLayers; i < m_nCells[2] - m_noGhostLayers; ++i) {
              MInt cellId = cellIndex(i, j, k);
              // MFloat x= m_cells->coordinates[0][cellId];
              MFloat y = m_cells->coordinates[1][cellId];
              MFloat z = m_cells->coordinates[2][cellId];
              MFloat phi = atan2(y, z);
              MFloat r = sqrt(POW2(y) + POW2(z));
              MFloat rmax = 10.0;
              m_cells->pvariables[PV->RHO][cellId] = CV->rhoInfinity;
              m_cells->pvariables[PV->U][cellId] = PV->UInfinity;
              m_cells->pvariables[PV->V][cellId] = -(r / rmax) * cos(phi) * 0.1 * PV->UInfinity;
              m_cells->pvariables[PV->W][cellId] = (r / rmax) * sin(phi) * 0.1 * PV->UInfinity;
              m_cells->pvariables[PV->P][cellId] = PV->PInfinity;
            }
          }
        }
        break;
      }
      case 42: {
        for(MInt cellId = 0; cellId < m_noCells; cellId++) {
          m_cells->pvariables[PV->RHO][cellId] = CV->rhoInfinity;
          m_cells->pvariables[PV->V][cellId] = F0;
          m_cells->pvariables[PV->W][cellId] = F0;

          if(m_cells->coordinates[1][cellId] < 0.5) {
            m_cells->pvariables[PV->U][cellId] = 1.001 * PV->UInfinity;
          } else {
            m_cells->pvariables[PV->U][cellId] = 0.999 * PV->UInfinity;
          }

          m_cells->pvariables[PV->P][cellId] = PV->PInfinity;
        }
        break;
      }
      case 44: {
        MFloat amp = 0.003;
        MFloat b = 0.25 / 4;
        MFloat alpha = log(2) / POW2(b);
        MFloat r0 = 0.25; // changed!!!
        MFloat r = F0;
        MFloat phi = F0;
        MFloat v = F0;
        for(MInt cellId = 0; cellId < m_noCells; cellId++) {
          // go through every cell
          m_cells->pvariables[PV->RHO][cellId] = CV->rhoInfinity;
          m_cells->pvariables[PV->P][cellId] = PV->PInfinity;
          m_cells->pvariables[PV->U][cellId] = PV->UInfinity;
          m_cells->pvariables[PV->V][cellId] = PV->VInfinity;
          m_cells->pvariables[PV->W][cellId] = PV->WInfinity;
          MFloat x = m_cells->coordinates[0][cellId] - 0.61;
          MFloat y = m_cells->coordinates[1][cellId] - 0.55;
          MFloat z = m_cells->coordinates[2][cellId] - 0.43;
          r = sqrt(POW2(y) + POW2(z));
          phi = atan(z / y);
          m_cells->pvariables[PV->U][cellId] +=
              amp * (r0 / mMax(r, 0.00001)) * (r - r0) * exp(-1.0 * alpha * (POW2(x) + POW2(r - r0)));
          v = -1.0 * amp * (r0 / mMax(r, 0.00001)) * x * exp(-1.0 * alpha * (POW2(x) + POW2(r - r0)));
          m_cells->pvariables[PV->V][cellId] += v * r * sin(phi);
          m_cells->pvariables[PV->W][cellId] += v * r * cos(phi);
        }
        break;
      }
      case 777: {
        // fsc
        m_log << "falkner skan cooke initial condition (incompressible)" << endl;
        if(!m_fsc) mTerm(1, "property fsc not set. Refer to the description of the property");
        if(true && !domainId()) {
          ofstream fscf;
          // write pressure to file
          fscf.open("pressure.dat", ios::trunc);
          if(fscf) {
            fscf << "#x_maia x_fsc p" << endl;
            for(MInt i = 0; i < 95; i++) {
              const MFloat x = i * F1;
              fscf << x << " " << m_fsc_x0 + x << " " << getFscPressure(x) / PV->PInfinity << endl;
            }
            fscf.close();
          }
          // write velocity to file
          fscf.open("velocity_x0.dat", ios::trunc);
          if(fscf) {
            fscf << "#y eta u v w" << endl;
            for(MInt i = 0; i < 200; i++) {
              // coord
              const MFloat y = i * 0.05;
              fscf << y << " " << getFscEta(F0, y);
              // velocity
              MFloat vel[nDim];
              getFscVelocity(F0, y, vel);
              for(MInt dim = 0; dim < nDim; dim++)
                fscf << " " << vel[dim] / PV->UInfinity;
              fscf << endl;
            }

            fscf.close();
          }
        }

        for(MInt cellid = 0; cellid < m_noCells; cellid++) {
          MFloat vel[nDim];
          getFscVelocity(cellid, vel);
          for(MInt i = 0; i < nDim; i++) {
            m_cells->pvariables[PV->VV[i]][cellid] = vel[i];
          }
          m_cells->pvariables[PV->P][cellid] = getFscPressure(cellid);
          m_cells->pvariables[PV->RHO][cellid] = CV->rhoInfinity;
        }

        break;
      }
      case 999: {
        // blasius laminar bl
        m_log << "Blasius initial condition (incompressible)" << endl;
        if(!m_useBlasius) mTerm(1, "property Blasius not set. Refer to the description of the property");
        if(domainId() == 0) {
          ofstream blasiusf;
          // write velocity to file
          blasiusf.open("velocity_x0.dat", ios::trunc);
          if(blasiusf) {
            blasiusf << "#y eta u v" << endl;
            MFloat d0 = 0.0;
            MFloat d1 = 0.0;
            MFloat d2 = 0.0;
            MBool d0Set = false;

            for(MInt i = 0; i < m_blasius_noPoints; i++) {
              // coord
              const MFloat y = i * 0.05;
              blasiusf << y << " " << getBlasiusEta(F0, y);
              // velocity
              MFloat vel[nDim];
              getBlasiusVelocity(F0, y, vel);
              for(MInt dim = 0; dim < nDim; dim++)
                blasiusf << " " << vel[dim] / PV->UInfinity;
              blasiusf << endl;

              if(!d0Set && vel[0] >= 0.99 * PV->UInfinity) {
                d0 = y;
                d0Set = true;
              }

              if(y < 10.0) {
                d1 += (1 - vel[0] / PV->UInfinity) * 0.05;
                d2 += vel[0] / PV->UInfinity * (1 - vel[0] / PV->UInfinity) * 0.05;
              }
            }
            blasiusf.close();

            cout << "x0: " << m_blasius_x0 << endl;
            cout << "d0: " << d0 << " d1: " << d1 << " d2: " << d2 << endl;
          }
        }

        for(MInt cellid = 0; cellid < m_noCells; cellid++) {
          MFloat vel[nDim];
          getBlasiusVelocity(cellid, vel);
          for(MInt i = 0; i < nDim; i++) {
            m_cells->pvariables[PV->VV[i]][cellid] = vel[i];
          }
          m_cells->pvariables[PV->P][cellid] = PV->PInfinity;
          m_cells->pvariables[PV->RHO][cellid] = CV->rhoInfinity;
        }

        break;
      }
      default: {
        // put the parallel flow field input in here
        // force output that no specific initial condition was chosen
        mTerm(1, AT_, "No (correct) initial Condition is given!");
        break;
      }
    }
  }
}


/**
 * \brief Compute all the interpolation coefficients necessary for the interpolation output
 */
void FvStructuredSolver3D::initInterpolatedPoints() {
  TRACE();
  if(!m_movingGrid) {
    m_pointInterpolation =
        make_unique<StructuredInterpolation<3>>(m_nCells, m_cells->coordinates, m_cells->pvariables, m_StructuredComm);
    // allocate domains points of the lines
    m_pointInterpolation->prepareInterpolation(m_intpPointsNoPointsTotal, m_intpPointsCoordinates,
                                               m_intpPointsHasPartnerLocal);
    MPI_Allreduce(m_intpPointsHasPartnerLocal, m_intpPointsHasPartnerGlobal, m_intpPointsNoPointsTotal, MPI_INT,
                  MPI_SUM, m_StructuredComm, AT_, "m_intpPointsHasPartnerLocal", "m_intpPointsHasPartnerGlobal");
  }
}


void FvStructuredSolver3D::initPointsToAsciiFile() {
  TRACE();
  if(!m_movingGrid) {
    m_pointsToAsciiInterpolation =
        make_unique<StructuredInterpolation<3>>(m_nCells, m_cells->coordinates, m_cells->pvariables, m_StructuredComm);

    // allocate domains points of the lines
    m_pointsToAsciiInterpolation->prepareInterpolation(m_pointsToAsciiNoPoints, m_pointsToAsciiCoordinates,
                                                       m_pointsToAsciiHasPartnerLocal);
    MPI_Allreduce(m_pointsToAsciiHasPartnerLocal, m_pointsToAsciiHasPartnerGlobal, m_pointsToAsciiNoPoints, MPI_INT,
                  MPI_SUM, m_StructuredComm, AT_, "m_pointsToAsciiHasPartnerLocal", "m_pointsToAsciiHasPartnerGlobal");

    for(MInt pointId = 0; pointId < m_pointsToAsciiNoPoints; pointId++) {
      cout << "domainId: " << domainId() << " pointId: " << pointId
           << " hasPartnerLocal: " << m_pointsToAsciiHasPartnerLocal[pointId]
           << " hasPartnerGlobal: " << m_pointsToAsciiHasPartnerGlobal[pointId]
           << " x: " << m_pointsToAsciiCoordinates[0][pointId] << " y: " << m_pointsToAsciiCoordinates[1][pointId]
           << " z: " << m_pointsToAsciiCoordinates[2][pointId] << endl;
    }
  }
}


/**
 * \brief Saves variables of given cells to ASCII file
 *
 * \author Marian Albers
 * \date  2020
 *
 */
void FvStructuredSolver3D::savePointsToAsciiFile(MBool forceWrite) {
  if(m_pointsToAsciiLastComputationStep != globalTimeStep) {
    if(m_movingGrid) {
      for(MInt pointId = 0; pointId < m_pointsToAsciiNoPoints; pointId++) {
        m_intpPointsHasPartnerLocal[pointId] = 0;
        m_intpPointsHasPartnerGlobal[pointId] = 0;
        for(MInt var = 0; var < PV->noVariables; var++) {
          m_intpPointsVarsLocal[var][pointId] = F0;
          m_intpPointsVarsGlobal[var][pointId] = F0;
        }
      }

      m_pointsToAsciiInterpolation = make_unique<StructuredInterpolation<3>>(m_nCells, m_cells->coordinates,
                                                                             m_cells->pvariables, m_StructuredComm);
      // allocate domains points of the lines
      m_pointsToAsciiInterpolation->prepareInterpolation(m_pointsToAsciiNoPoints, m_pointsToAsciiCoordinates,
                                                         m_pointsToAsciiHasPartnerLocal);
      MPI_Allreduce(m_pointsToAsciiHasPartnerLocal, m_pointsToAsciiHasPartnerGlobal, m_intpPointsNoPointsTotal, MPI_INT,
                    MPI_SUM, m_StructuredComm, AT_, "m_pointsToAsciiHasPartnerLocal",
                    "m_pointsToAsciiHasPartnerGlobal");
    }

    // calculation of interpolated variables only for the points in domain
    for(MInt pointId = 0; pointId < m_pointsToAsciiNoPoints; pointId++) {
      m_pointsToAsciiVars[m_pointsToAsciiCounter][3 + pointId] = F0;
      if(m_pointsToAsciiHasPartnerLocal[pointId]) {
        m_pointsToAsciiVars[m_pointsToAsciiCounter][3 + pointId] =
            m_pointsToAsciiInterpolation->getInterpolatedVariable(pointId, m_pointsToAsciiVarId);
      }
    }

    MPI_Allreduce(MPI_IN_PLACE, &m_pointsToAsciiVars[m_pointsToAsciiCounter][3], m_pointsToAsciiNoPoints, MPI_DOUBLE,
                  MPI_SUM, m_StructuredComm, AT_, "m_pointsToAsciiVars[0][0]", "m_pointsToAsciiVars[0][0]");

    for(MInt pointId = 0; pointId < m_pointsToAsciiNoPoints; pointId++) {
      if(m_pointsToAsciiHasPartnerGlobal[pointId] > 1) {
        m_pointsToAsciiVars[m_pointsToAsciiCounter][3 + pointId] /= (MFloat)m_pointsToAsciiHasPartnerGlobal[pointId];
      }
    }

    m_pointsToAsciiVars[m_pointsToAsciiCounter][0] = globalTimeStep;
    m_pointsToAsciiVars[m_pointsToAsciiCounter][1] = m_time;
    m_pointsToAsciiVars[m_pointsToAsciiCounter][2] = m_physicalTime;

    m_pointsToAsciiCounter++;
    m_pointsToAsciiLastComputationStep = globalTimeStep;
  }


  if((m_pointsToAsciiCounter >= m_pointsToAsciiOutputInterval || forceWrite)
     && m_pointsToAsciiLastOutputStep != globalTimeStep) {
    if(domainId() == 0) {
      cout << "globalTimeStep: " << globalTimeStep << " writing point data out to ascii file" << endl;
      MString filename = "./pointVars.dat";
      FILE* f_forces;
      f_forces = fopen(filename.c_str(), "a+");

      for(MInt j = 0; j < m_pointsToAsciiCounter; j++) {
        fprintf(f_forces, "%d", (MInt)m_pointsToAsciiVars[j][0]);
        fprintf(f_forces, " %.8f", m_pointsToAsciiVars[j][1]);
        fprintf(f_forces, " %.8f", m_pointsToAsciiVars[j][2]);
        for(MInt i = 0; i < m_pointsToAsciiNoPoints; i++) {
          fprintf(f_forces, " %.8f", m_pointsToAsciiVars[j][3 + i]);
        }
        fprintf(f_forces, "\n");
      }
      fclose(f_forces);
    }

    m_pointsToAsciiCounter = 0;
    m_pointsToAsciiLastOutputStep = globalTimeStep;
  }
}


/**
 * \brief Manually correct errors made by the restart interpolation
 *
 *  In case some cells did not get correct values in
 *  the interpolation process (due to no matching donor cells)
 *  Put your routines that will correct these cells here.
 */
void FvStructuredSolver3D::manualInterpolationCorrection() {
  /*! \property
    \page propertiesFVSTRCTRD
    \section interpolationCorrection
    <code>MBool interpolationCorrection </code>\n
    default = <code> 0 </code>\n \n
    Trigger a manual correction of the interpolation\n
    from a donor grid.\n
    Possible values are:\n
    <ul>
    <li>true/false</li>
    </ul>
    Keywords: <i>LIMITER, STRUCTURED</i>
  */
  MBool interpolationCorrection = false;
  if(Context::propertyExists("interpolationCorrection", m_solverId)) {
    interpolationCorrection =
        Context::getSolverProperty<MBool>("interpolationCorrection", m_solverId, AT_, &interpolationCorrection);
  }
  if(m_nOffsetCells[1] == 0 && m_rans && interpolationCorrection) {
    cout << "Correcting interpolation" << endl;
    for(MInt i = 0; i < m_nCells[2]; i++) {
      for(MInt k = 0; k < m_nCells[0]; k++) {
        MInt cellIdA1 = cellIndex(i, 2, k);
        MInt cellIdA2 = cellIndex(i, 3, k);
        MInt cellIdA3 = cellIndex(i, 4, k);

        for(MInt var = 0; var < PV->noVariables; var++) {
          m_cells->pvariables[var][cellIdA1] =
              m_cells->pvariables[var][cellIdA2]
              + (m_cells->coordinates[1][cellIdA1] - m_cells->coordinates[1][cellIdA2])
                    / (m_cells->coordinates[1][cellIdA2] - m_cells->coordinates[1][cellIdA3])
                    * (m_cells->variables[var][cellIdA2] - m_cells->pvariables[var][cellIdA3]);
        }
      }
    }
  }
}


/**
 * \brief Returns a normal distributed random-number with mu=mean and sigma=standard deviation
 *
 */
MFloat FvStructuredSolver3D::randnormal(MFloat mu, MFloat sigma) {
  TRACE();
  static MBool deviateAvailable = false; //        flag
  static float storedDeviate;            //        deviate from previous calculation
  MFloat polar, rsquared, var1, var2;

  // If no deviate has been stored, the polar Box-Muller transformation is
  // performed, producing two independent normally-distributed random
  // deviates.  One is stored for the next round, and one is returned.
  if(!deviateAvailable) {
    // choose pairs of uniformly distributed deviates, discarding those
    // that don't fall within the unit circle
    do {
      var1 = 2.0 * (MFloat(rand()) / MFloat(RAND_MAX)) - 1.0;
      var2 = 2.0 * (MFloat(rand()) / MFloat(RAND_MAX)) - 1.0;
      rsquared = var1 * var1 + var2 * var2;
    } while(rsquared >= 1.0 || approx(rsquared, F0, m_eps));

    // calculate polar tranformation for each deviate
    polar = sqrt(-2.0 * log(rsquared) / rsquared);

    // store first deviate and set flag
    storedDeviate = var1 * polar;
    deviateAvailable = true;

    // return second deviate
    return var2 * polar * sigma + mu;

    // If a deviate is available from a previous call to this function, it is
    // eturned, and the flag is set to false.
  } else {
    deviateAvailable = false;
    return storedDeviate * sigma + mu;
  }
}

void FvStructuredSolver3D::nonReflectingBC() {
  TRACE();
  m_structuredBndryCnd->applyNonReflectingBC();
}

void FvStructuredSolver3D::applyBoundaryCondition() {
  TRACE();
  // treat Dirichlet and Neumann BC in one go!!!
  m_structuredBndryCnd->applyDirichletNeumannBC();
}


void FvStructuredSolver3D::MusclRANS() { m_ransSolver->Muscl(); }


void FvStructuredSolver3D::MusclMinModLimiter() {
  TRACE();
  // stencil identifier
  const MInt IJK[3] = {1, m_nCells[2], m_nCells[1] * m_nCells[2]};
  // switch for order
  const MInt sword = 2;

  // reduce to onedimensional arrays
  MFloat* __restrict x = &m_cells->coordinates[0][0];
  MFloat* __restrict y = &m_cells->coordinates[1][0];
  MFloat* __restrict z = &m_cells->coordinates[2][0];
  MFloat* const* const RESTRICT flux = ALIGNED_F(m_cells->flux);

  const MUint noCells = m_noCells;
  const MFloat* const RESTRICT cellVariables = ALIGNED_F(m_cells->pvariables[0]);

  for(MInt dim = 0; dim < nDim; ++dim) {
    for(MInt k = m_noGhostLayers - 1; k < m_nCells[0] - m_noGhostLayers; ++k) {
      for(MInt j = m_noGhostLayers - 1; j < m_nCells[1] - m_noGhostLayers; ++j) {
        for(MInt i = m_noGhostLayers - 1; i < m_nCells[2] - m_noGhostLayers; ++i) {
          // cell ids
          const MInt I = cellIndex(i, j, k);
          const MInt IP1 = I + IJK[dim];
          const MInt IM1 = I - IJK[dim];
          const MInt IP2 = I + 2 * IJK[dim];

          // distances q_i+1 - q_i
          const MFloat DS = sqrt(POW2(x[IP1] - x[I]) + POW2(y[IP1] - y[I]) + POW2(z[IP1] - z[I]));
          // distances q_i - q_i-1
          const MFloat DSM1 = sqrt(POW2(x[I] - x[IM1]) + POW2(y[I] - y[IM1]) + POW2(z[I] - z[IM1]));
          const MFloat DSP1 = sqrt(POW2(x[IP2] - x[IP1]) + POW2(y[IP2] - y[IP1]) + POW2(z[IP2] - z[IP1]));
          const MFloat DSP = DS / POW2(DSP1 + DS);
          const MFloat DSM = DS / POW2(DSM1 + DS);

          if(sword == 2) {
            for(MInt var = 0; var < PV->noVariables; ++var) {
              const MUint offset = var * noCells;
              const MFloat* const RESTRICT vars = ALIGNED_F(cellVariables + offset);

              const MFloat DQ = vars[IP1] - vars[I];
              const MFloat DQP1 = vars[IP2] - vars[IP1];
              const MFloat DQM1 = vars[I] - vars[IM1];

              const MFloat ri = DQM1 / DQ;
              const MFloat rip = DQ / DQP1;

              const MFloat phii = mMax(F0, mMin(F1, ri));
              const MFloat phiip = mMax(F0, mMin(F1, rip));

              m_QLeft[var] = vars[I] + (DQ * DSM1 + DQM1 * DS) * DSM * phii;
              m_QRight[var] = vars[IP1] - (DQP1 * DS + DQ * DSP1) * DSP * phiip;
            }
          }

          AusmLES(m_QLeft, m_QRight, dim, I);
        }
      }
    }

    // FLUX BALANCE
    for(MInt k = m_noGhostLayers; k < m_nCells[0] - m_noGhostLayers; ++k) {
      for(MInt j = m_noGhostLayers; j < m_nCells[1] - m_noGhostLayers; ++j) {
        for(MInt i = m_noGhostLayers; i < m_nCells[2] - m_noGhostLayers; ++i) {
          const MInt I = cellIndex(i, j, k);
          const MInt IM1 = I - IJK[dim];

          for(MInt v = 0; v < CV->noVariables; ++v) {
            m_cells->rightHandSide[v][I] += flux[v][IM1] - flux[v][I];
          }
        }
      }
    }
  }
}

// Muscl reconstruction with Albada limiter
void FvStructuredSolver3D::MusclAlbada() {
  TRACE();
  // stencil identifier
  const MInt IJK[3] = {1, m_nCells[2], m_nCells[1] * m_nCells[2]};

  // reduce to onedimensional arrays
  MFloat* RESTRICT x = &m_cells->coordinates[0][0];
  MFloat* RESTRICT y = &m_cells->coordinates[1][0];
  MFloat* RESTRICT z = &m_cells->coordinates[2][0];
  MFloat* const* const RESTRICT flux = ALIGNED_F(m_cells->flux);
  MFloat** RESTRICT pvars = m_cells->pvariables;

  /////////IMPORTANT PARAMETER
  // MFloat epsi=F1;
  // MFloat kappa=F1B3;
  /////////END IMPORTANT PARAMETER
  for(MInt dim = 0; dim < nDim; dim++) {
    for(MInt k = m_noGhostLayers - 1; k < m_nCells[0] - m_noGhostLayers; k++) {
      for(MInt j = m_noGhostLayers - 1; j < m_nCells[1] - m_noGhostLayers; j++) {
        for(MInt i = m_noGhostLayers - 1; i < m_nCells[2] - m_noGhostLayers; i++) {
          // cell ids
          const MInt I = cellIndex(i, j, k);
          const MInt IP1 = I + IJK[dim];
          const MInt IM1 = I - IJK[dim];
          const MInt IP2 = I + 2 * IJK[dim];

          // distances q_i+1 - q_i
          const MFloat DS = sqrt(POW2(x[IP1] - x[I]) + POW2(y[IP1] - y[I]) + POW2(z[IP1] - z[I]));
          // distances q_i - q_i-1
          const MFloat DSM1 = sqrt(POW2(x[I] - x[IM1]) + POW2(y[I] - y[IM1]) + POW2(z[I] - z[IM1]));
          const MFloat DSP1 = sqrt(POW2(x[IP2] - x[IP1]) + POW2(y[IP2] - y[IP1]) + POW2(z[IP2] - z[IP1]));
          const MFloat DSP = DS / POW2(DSP1 + DS);
          const MFloat DSM = DS / POW2(DSM1 + DS);

          const MFloat pIM2 = pvars[PV->P][IM1];
          const MFloat pIM1 = pvars[PV->P][I];
          const MFloat pIP2 = pvars[PV->P][IP2];
          const MFloat pIP1 = pvars[PV->P][IP1];

          const MFloat smps = DS * DSP1;
          const MFloat dummy = fabs(pIM2 - F2 * pIM1 + pIP1) / (pIM2 + F2 * pIM1 + pIP1);
          const MFloat dummy1 = fabs(pIM1 - F2 * pIP1 + pIP2) / (pIM1 + F2 * pIP1 + pIP2);
          const MFloat psi = mMin(F1, F6 * mMax(dummy, dummy1));
          const MFloat epsLim = mMax(m_eps, pow(F1B2 * smps, F5));

          for(MInt var = 0; var < PV->noVariables; ++var) {
            const MFloat DQ = pvars[var][IP1] - pvars[var][I];
            const MFloat DQP1 = pvars[var][IP2] - pvars[var][IP1];
            const MFloat DQM1 = pvars[var][I] - pvars[var][IM1];
            const MFloat phi =
                F1B2
                - (F1B2
                   - mMax(F0, (DQP1 * DQM1 * smps + F1B2 * epsLim) / (POW2(DQP1 * DS) + POW2(DQM1 * DSP1) + epsLim)))
                      * psi;

            m_QLeft[var] = pvars[var][I] + DSM * (DSM1 * DQ + DS * DQM1) * phi;
            m_QRight[var] = pvars[var][IP1] - DSP * (DS * DQP1 + DSP1 * DQ) * phi;
          }

          AusmLES(m_QLeft, m_QRight, dim, I); // Flux balance in AUSM
        }
      }
    }

    // FLUX BALANCE
    for(MInt v = 0; v < CV->noVariables; v++) {
      for(MInt k = m_noGhostLayers; k < m_nCells[0] - m_noGhostLayers; k++) {
        for(MInt j = m_noGhostLayers; j < m_nCells[1] - m_noGhostLayers; j++) {
          for(MInt i = m_noGhostLayers; i < m_nCells[2] - m_noGhostLayers; i++) {
            const MInt I = cellIndex(i, j, k);
            const MInt IM1 = I - IJK[dim];
            m_cells->rightHandSide[v][I] += flux[v][IM1] - flux[v][I];
          }
        }
      }
    }
  }
}


/**
 * \brief MUSCL with Venkatakrishan limiter
 *
 * Here, MUSCL and AUSM are run through separately
 * The values of QLeft and QRight (of every cell) from the MUSCL are stored in the scratchspace
 * before passing it to the AUSM scheme.
 * Pros: the limiter can include the values of the neighbouring cells, better results for waves
 * Cons: higher computational effort
 *
 * \author Leo Hoening
 * \date 2015
 */
void FvStructuredSolver3D::MusclVenkatakrishnan3D() {
  TRACE();
  const MUint noCells = m_noCells;
  const MInt IJK[3] = {1, m_nCells[2], m_nCells[1] * m_nCells[2]};
  // reduce to onedimensional arrays
  const MFloat* const RESTRICT x = ALIGNED_F(m_cells->coordinates[0]);
  const MFloat* const RESTRICT y = ALIGNED_F(m_cells->coordinates[1]);
  const MFloat* const RESTRICT z = ALIGNED_F(m_cells->coordinates[2]);
  const MFloat* const RESTRICT cellVariables = ALIGNED_F(m_cells->pvariables[0]);
  MFloat* const* const RESTRICT flux = ALIGNED_F(m_cells->flux);

  MFloatScratchSpace QLeft(m_noCells, PV->noVariables, 3, AT_, "QLeft");
  MFloatScratchSpace QRight(m_noCells, PV->noVariables, 3, AT_, "QRight");
  MFloatScratchSpace minPhi(PV->noVariables, 2, AT_, "minPhi");

  QLeft.fill(F0);
  QRight.fill(F0);

  for(MInt k = m_noGhostLayers - 1; k < m_nCells[0] - m_noGhostLayers; k++) {
    for(MInt j = m_noGhostLayers - 1; j < m_nCells[1] - m_noGhostLayers; j++) {
      for(MInt i = m_noGhostLayers - 1; i < m_nCells[2] - m_noGhostLayers; i++) {
        for(MInt var = 0; var < nDim + 2; var++) {
          const MUint offset = var * noCells;
          const MFloat* const RESTRICT pvars = ALIGNED_F(cellVariables + offset);

          MFloat minNghbrDelta = F0;
          MFloat maxNghbrDelta = F0;
          MFloat effNghbrDelta = F0;

          // 1. get the abs min value of the max and min value from the reconstruction neighbours
          for(MInt dim1 = 0; dim1 < nDim; dim1++) {
            for(MInt rcnstructnNghbr = 0; rcnstructnNghbr < 2; rcnstructnNghbr++) {
              const MInt cellId = cellIndex(i, j, k);
              const MInt IP1 = cellId + IJK[dim1];
              const MInt IM1 = cellId - IJK[dim1];

              MFloat rcnstrctnNghbrValue = F0;
              if(rcnstructnNghbr == 0) {
                rcnstrctnNghbrValue = pvars[IP1]; //(i/j/k)+1
              } else {
                rcnstrctnNghbrValue = pvars[IM1]; //(i/j/k)-1
              }

              const MFloat tmpDelta = rcnstrctnNghbrValue - pvars[cellId]; // i
              maxNghbrDelta = mMax(maxNghbrDelta, tmpDelta);
              minNghbrDelta = mMin(minNghbrDelta, tmpDelta);
            }
          }

          effNghbrDelta = mMin(maxNghbrDelta, abs(minNghbrDelta));

          MFloat srfcDelta = F0;
          MFloat dxEpsSqr = F1;
          for(MInt dim1 = 0; dim1 < nDim; dim1++) {
            const MInt cellId = cellIndex(i, j, k);
            const MInt IP1 = cellId + IJK[dim1];
            const MInt IM1 = cellId - IJK[dim1];

            const MFloat DS = sqrt(POW2(x[IP1] - x[cellId]) + POW2(y[IP1] - y[cellId]) + POW2(z[IP1] - z[cellId]));
            // distances q_i - q_i-1
            const MFloat DSM1 = sqrt(POW2(x[cellId] - x[IM1]) + POW2(y[cellId] - y[IM1]) + POW2(z[cellId] - z[IM1]));
            const MFloat DSM = DS / POW2(DSM1 + DS);

            // 2. get srfcDelta and compute the minimum phi
            const MFloat dx1 =
                DSM
                * (DSM1 * sqrt(POW2(x[IP1] - x[cellId]) + POW2(y[IP1] - y[cellId]) + POW2(z[IP1] - z[cellId]))
                   + DS * sqrt(POW2(x[cellId] - x[IM1]) + POW2(y[cellId] - y[IM1]) + POW2(z[cellId] - z[IM1])));
            const MFloat DQ =
                (pvars[IP1] - pvars[IM1]) / sqrt(POW2(x[IP1] - x[IM1]) + POW2(y[IP1] - y[IM1]) + POW2(z[IP1] - z[IM1]));
            srfcDelta += abs(DQ * dx1);
            dxEpsSqr *= dx1;
          }

          MInt cellPos = 0;

          // calling limiter function
          (this->*Venkatakrishnan_function)(effNghbrDelta, srfcDelta, dxEpsSqr, cellPos, var, minPhi);

          minNghbrDelta = F0;
          maxNghbrDelta = F0;
          effNghbrDelta = F0;

          // 1. get the abs min value of the max and min value from the reconstruction neighbours
          for(MInt dim1 = 0; dim1 < nDim; dim1++) {
            for(MInt rcnstructnNghbr = 0; rcnstructnNghbr < 2; rcnstructnNghbr++) {
              const MInt cellId = cellIndex(i, j, k);
              const MInt IP1 = cellId + IJK[dim1];
              const MInt IP2 = cellId + 2 * IJK[dim1];

              MFloat rcnstrctnNghbrValue = F0;
              if(rcnstructnNghbr == 0) {
                rcnstrctnNghbrValue = pvars[IP2]; //(i/j/k)+2
              } else {
                rcnstrctnNghbrValue = pvars[cellId]; //(i/j/k)
              }

              const MFloat tmpDelta = rcnstrctnNghbrValue - pvars[IP1]; //(i/j/k)+1

              maxNghbrDelta = mMax(maxNghbrDelta, tmpDelta);
              minNghbrDelta = mMin(minNghbrDelta, tmpDelta);
            }
          }

          effNghbrDelta = mMin(maxNghbrDelta, abs(minNghbrDelta));

          srfcDelta = F0;
          dxEpsSqr = F1;
          for(MInt dim1 = 0; dim1 < nDim; dim1++) {
            const MInt cellId = cellIndex(i, j, k);
            const MInt IP1 = cellId + IJK[dim1];
            const MInt IP2 = cellId + 2 * IJK[dim1];

            const MFloat DS = sqrt(POW2(x[IP1] - x[cellId]) + POW2(y[IP1] - y[cellId]) + POW2(z[IP1] - z[cellId]));
            // distances q_i - q_i-1
            const MFloat DSP1 = sqrt(POW2(x[IP2] - x[IP1]) + POW2(y[IP2] - y[IP1]) + POW2(z[IP2] - z[IP1]));
            const MFloat DSP = DS / POW2(DSP1 + DS);

            // 2. get srfcDelta and compute the minimum phi

            const MFloat dx2 =
                DSP
                * (DS * sqrt(POW2(x[IP2] - x[IP1]) + POW2(y[IP2] - y[IP1]) + POW2(z[IP2] - z[IP1]))
                   + DSP1 * sqrt(POW2(x[IP1] - x[cellId]) + POW2(y[IP1] - y[cellId]) + POW2(z[IP1] - z[cellId])));
            const MFloat DQ = (pvars[IP2] - pvars[cellId])
                              / sqrt(POW2(x[IP2] - x[cellId]) + POW2(y[IP2] - y[cellId]) + POW2(z[IP2] - z[cellId]));

            srfcDelta += abs(DQ * dx2);
            dxEpsSqr *= dx2;
          }

          cellPos = 1;
          (this->*Venkatakrishnan_function)(effNghbrDelta, srfcDelta, dxEpsSqr, cellPos, var, minPhi); // calling Venk

          for(MInt dim1 = 0; dim1 < nDim; dim1++) {
            const MInt cellId = cellIndex(i, j, k);
            const MInt IP1 = cellId + IJK[dim1];
            const MInt IM1 = cellId - IJK[dim1];
            const MInt IP2 = cellId + 2 * IJK[dim1];

            const MFloat DS = sqrt(POW2(x[IP1] - x[cellId]) + POW2(y[IP1] - y[cellId]) + POW2(z[IP1] - z[cellId]));
            // distances q_i - q_i-1
            const MFloat DSM1 = sqrt(POW2(x[cellId] - x[IM1]) + POW2(y[cellId] - y[IM1]) + POW2(z[cellId] - z[IM1]));
            const MFloat DSP1 = sqrt(POW2(x[IP2] - x[IP1]) + POW2(y[IP2] - y[IP1]) + POW2(z[IP2] - z[IP1]));
            const MFloat DSP = DS / POW2(DSP1 + DS);
            const MFloat DSM = DS / POW2(DSM1 + DS);

            QLeft(cellId, var, dim1) =
                pvars[cellId]
                + DSM * (DSM1 * (pvars[IP1] - pvars[cellId]) + DS * (pvars[cellId] - pvars[IM1])) * minPhi(var, 0);
            QRight(cellId, var, dim1) =
                pvars[IP1]
                - DSP * (DS * (pvars[IP2] - pvars[IP1]) + DSP1 * (pvars[IP1] - pvars[cellId])) * minPhi(var, 1);
          }
        }
      }
    }
  }

  ////// AUSM /////////
  for(MInt dim = 0; dim < nDim; dim++) {
    for(MInt k = m_noGhostLayers - 1; k < m_nCells[0] - m_noGhostLayers; k++) {
      for(MInt j = m_noGhostLayers - 1; j < m_nCells[1] - m_noGhostLayers; j++) {
        for(MInt i = m_noGhostLayers - 1; i < m_nCells[2] - m_noGhostLayers; i++) {
          // cell ids
          const MInt cellId = cellIndex(i, j, k);

          for(MInt v = 0; v < PV->noVariables; v++) {
            m_QLeft[v] = QLeft(cellId, v, dim);
            m_QRight[v] = QRight(cellId, v, dim);
          }

          AusmLES(m_QLeft, m_QRight, dim, cellId);
        }
      }
    }


    // FLUX BALANCE
    for(MInt v = 0; v < CV->noVariables; v++) {
      for(MInt k = m_noGhostLayers; k < m_nCells[0] - m_noGhostLayers; k++) {
        for(MInt j = m_noGhostLayers; j < m_nCells[1] - m_noGhostLayers; j++) {
          for(MInt i = m_noGhostLayers; i < m_nCells[2] - m_noGhostLayers; i++) {
            const MInt I = cellIndex(i, j, k);
            const MInt IM1 = I - IJK[dim];
            m_cells->rightHandSide[v][I] += flux[v][IM1] - flux[v][I];
          }
        }
      }
    }
  }
}


/**
 * \brief Venkatakrishnan limiter, modified for better results
 *
 * \author Leo Hoening
 * \date 2015
 */
void FvStructuredSolver3D::VENKATAKRISHNAN_MOD_FCT(MFloat effNghbrDelta, MFloat srfcDelta, MFloat dxEpsSqr,
                                                   MInt cellPos, MInt var, MFloatScratchSpace& minPhi) {
  MFloat epsSqr = pow(m_venkFactor, F3) * dxEpsSqr;
  minPhi(var, cellPos) = (pow(effNghbrDelta, F2) + epsSqr + F2 * effNghbrDelta * srfcDelta)
                         / (pow(effNghbrDelta, F2) + F2 * pow(srfcDelta, F2) + effNghbrDelta * srfcDelta + epsSqr);
}

/**
 * \brief Standard Venkatakrishnan limiter
 * \author Leo Hoening
 * \date 2015
 */
void FvStructuredSolver3D::VENKATAKRISHNAN_FCT(MFloat effNghbrDelta, MFloat srfcDelta, MFloat dxEpsSqr, MInt cellPos,
                                               MInt var, MFloatScratchSpace& minPhi) {
  (void)dxEpsSqr;
  const MFloat eps = 1e-12;
  MFloat yps1 = effNghbrDelta / (srfcDelta + eps);
  minPhi(var, cellPos) = mMin((yps1 * yps1 + F2 * yps1) / (yps1 * yps1 + yps1 + F2), F1);
}

/**
 * \brief Barth-Jesperson Limiter
 * \author Leo Hoening
 * \date 2015
 */
void FvStructuredSolver3D::BARTH_JESPERSON_FCT(MFloat effNghbrDelta, MFloat srfcDelta, MFloat dxEpsSqr, MInt cellPos,
                                               MInt var, MFloatScratchSpace& minPhi) {
  (void)dxEpsSqr;
  const MFloat eps = 1e-12;
  MFloat phi_max = effNghbrDelta / (srfcDelta + eps);
  minPhi(var, cellPos) = mMin(phi_max, F1);
}

/**
 * \brief AUSM central
 *
 *  Can be used for moving grids, dxt term is included
 */
// inline void FvStructuredSolver3D::AusmNew(MFloat* QLeft, MFloat* QRight, const MInt dim, const MInt
// cellId)
inline void FvStructuredSolver3D::AusmLES(MFloat* RESTRICT QLeft, MFloat* RESTRICT QRight, const MInt dim,
                                          const MInt I) {
  // MFloat pFactor[3]={F0,F0,F0};
  const MFloat gamma = m_gamma;
  const MFloat gammaMinusOne = gamma - 1.0;
  const MFloat FgammaMinusOne = F1 / gammaMinusOne;

  const MFloat surf0 = m_cells->surfaceMetrics[dim * 3 + 0][I];
  const MFloat surf1 = m_cells->surfaceMetrics[dim * 3 + 1][I];
  const MFloat surf2 = m_cells->surfaceMetrics[dim * 3 + 2][I];

  const MFloat dxdtau = m_cells->dxt[dim][I];

  // calculate pressure
  const MFloat PL = QLeft[PV->P];
  const MFloat UL = QLeft[PV->U];
  const MFloat VL = QLeft[PV->V];
  const MFloat WL = QLeft[PV->W];
  const MFloat RHOL = QLeft[PV->RHO];

  const MFloat PR = QRight[PV->P];
  const MFloat UR = QRight[PV->U];
  const MFloat VR = QRight[PV->V];
  const MFloat WR = QRight[PV->W];
  const MFloat RHOR = QRight[PV->RHO];

  // compute lenght of metric vector for normalization
  const MFloat metricLength = sqrt(POW2(surf0) + POW2(surf1) + POW2(surf2));
  const MFloat fMetricLength = F1 / metricLength;

  // scale by metric length to get velocity in the new basis (get normalized basis vectors)
  const MFloat UUL = ((UL * surf0 + VL * surf1 + WL * surf2) - dxdtau) * fMetricLength;


  const MFloat UUR = ((UR * surf0 + VR * surf1 + WR * surf2) - dxdtau) * fMetricLength;


  // speed of sound
  const MFloat AL = sqrt(gamma * mMax(m_eps, PL / mMax(m_eps, RHOL)));
  const MFloat AR = sqrt(gamma * mMax(m_eps, PR / mMax(m_eps, RHOR)));

  const MFloat MAL = UUL / AL;
  const MFloat MAR = UUR / AR;

  const MFloat MALR = F1B2 * (MAL + MAR);
  const MFloat PLR = PL * (F1B2 + m_chi * MAL) + PR * (F1B2 - m_chi * MAR);

  const MFloat RHO_AL = RHOL * AL;
  const MFloat RHO_AR = RHOR * AR;

  const MFloat PLfRHOL = PL / RHOL;
  const MFloat PRfRHOR = PR / RHOR;

  const MFloat e0 = PLfRHOL * FgammaMinusOne + 0.5 * (POW2(UL) + POW2(VL) + POW2(WL)) + PLfRHOL;
  const MFloat e1 = PRfRHOR * FgammaMinusOne + 0.5 * (POW2(UR) + POW2(VR) + POW2(WR)) + PRfRHOR;

  const MFloat RHOU = F1B2 * (MALR * (RHO_AL + RHO_AR) + fabs(MALR) * (RHO_AL - RHO_AR)) * metricLength;
  const MFloat RHOU2 = F1B2 * RHOU;
  // multiply by metric length to take surface area into account
  const MFloat AbsRHO_U2 = fabs(RHOU2);

  //==>fluxes:
  MFloat* const* const RESTRICT flux = ALIGNED_F(m_cells->flux);
  flux[CV->RHO_U][I] = RHOU2 * (UL + UR) + AbsRHO_U2 * (UL - UR) + PLR * surf0;
  flux[CV->RHO_V][I] = RHOU2 * (VL + VR) + AbsRHO_U2 * (VL - VR) + PLR * surf1;
  flux[CV->RHO_W][I] = RHOU2 * (WL + WR) + AbsRHO_U2 * (WL - WR) + PLR * surf2;
  flux[CV->RHO_E][I] = RHOU2 * (e0 + e1) + AbsRHO_U2 * (e0 - e1) + PLR * dxdtau;
  flux[CV->RHO][I] = RHOU;
}


/**
 * \brief AUSM PTHRC
 *
 * Same AUSM scheme as AusmLES with additional damping controlled
 * by the 4th order pressure derivative. Pressure needs to computed
 * beforehand.
 *
 */
inline void FvStructuredSolver3D::AusmLES_PTHRC(MFloat* QLeft, MFloat* QRight, MInt dim, MInt I) {
  const MFloat gamma = m_gamma;
  const MFloat FgammaMinusOne = m_fgammaMinusOne;

  const MFloat surf0 = m_cells->surfaceMetrics[dim * 3 + 0][I];
  const MFloat surf1 = m_cells->surfaceMetrics[dim * 3 + 1][I];
  const MFloat surf2 = m_cells->surfaceMetrics[dim * 3 + 2][I];

  const MFloat* const RESTRICT p = ALIGNED_F(m_cells->pvariables[PV->P]);

  const MFloat dxdtau = m_cells->dxt[dim][I];

  // calculate pressure
  const MFloat PL = QLeft[PV->P];
  const MFloat UL = QLeft[PV->U];
  const MFloat VL = QLeft[PV->V];
  const MFloat WL = QLeft[PV->W];
  const MFloat RHOL = QLeft[PV->RHO];

  const MFloat PR = QRight[PV->P];
  const MFloat UR = QRight[PV->U];
  const MFloat VR = QRight[PV->V];
  const MFloat WR = QRight[PV->W];
  const MFloat RHOR = QRight[PV->RHO];

  // compute lenght of metric vector for normalization
  const MFloat metricLength = sqrt(POW2(surf0) + POW2(surf1) + POW2(surf2));
  const MFloat fMetricLength = F1 / metricLength;

  // scale by metric length to get velocity in the new basis (get normalized basis vectors)
  const MFloat UUL = ((UL * surf0 + VL * surf1 + WL * surf2) - dxdtau) * fMetricLength;


  const MFloat UUR = ((UR * surf0 + VR * surf1 + WR * surf2) - dxdtau) * fMetricLength;

  // speed of sound
  const MFloat AL = sqrt(gamma * mMax(m_eps, PL / mMax(m_eps, RHOL)));
  const MFloat AR = sqrt(gamma * mMax(m_eps, PR / mMax(m_eps, RHOR)));

  const MFloat MAL = UUL / AL;
  const MFloat MAR = UUR / AR;

  const MFloat MALR = F1B2 * (MAL + MAR);

  // 4th order pressure damping
  const MInt IPJK = getCellIdfromCell(I, 1, 0, 0);
  const MInt IMJK = getCellIdfromCell(I, -1, 0, 0);
  const MInt IP2JK = getCellIdfromCell(I, 2, 0, 0);
  const MInt IM2JK = getCellIdfromCell(I, -2, 0, 0);

  const MInt IJPK = getCellIdfromCell(I, 0, 1, 0);
  const MInt IJMK = getCellIdfromCell(I, 0, -1, 0);
  const MInt IJP2K = getCellIdfromCell(I, 0, 2, 0);
  const MInt IJM2K = getCellIdfromCell(I, 0, -2, 0);

  const MInt IJKP = getCellIdfromCell(I, 0, 0, 1);
  const MInt IJKM = getCellIdfromCell(I, 0, 0, -1);
  const MInt IJKP2 = getCellIdfromCell(I, 0, 0, 2);
  const MInt IJKM2 = getCellIdfromCell(I, 0, 0, -2);

  const MFloat p4I4 = F4 * (p[IPJK] + p[IMJK]) - F6 * (p[I]) - p[IP2JK] - p[IM2JK];
  const MFloat p4J4 = F4 * (p[IJPK] + p[IJMK]) - F6 * (p[I]) - p[IJP2K] - p[IJM2K];
  const MFloat p4K4 = F4 * (p[IJKP] + p[IJKM]) - F6 * (p[I]) - p[IJKP2] - p[IJKM2];

  const MFloat pfac = fabs(p4I4) + fabs(p4J4) + fabs(p4K4);
  const MFloat facl = 1.0 / 1.3 * pfac;
  const MFloat fac = min(1 / 128.0, facl);

  const MFloat PLR = PL * (F1B2 + fac * MAL) + PR * (F1B2 - fac * MAR);

  const MFloat RHO_AL = RHOL * AL;
  const MFloat RHO_AR = RHOR * AR;

  const MFloat PLfRHOL = PL / RHOL;
  const MFloat PRfRHOR = PR / RHOR;

  const MFloat e0 = PLfRHOL * FgammaMinusOne + 0.5 * (POW2(UL) + POW2(VL) + POW2(WL)) + PLfRHOL;
  const MFloat e1 = PRfRHOR * FgammaMinusOne + 0.5 * (POW2(UR) + POW2(VR) + POW2(WR)) + PRfRHOR;

  const MFloat RHOU = F1B2 * (MALR * (RHO_AL + RHO_AR) + fabs(MALR) * (RHO_AL - RHO_AR)) * metricLength;
  const MFloat RHOU2 = F1B2 * RHOU;
  // multiply by metric length to take surface area into account
  const MFloat AbsRHO_U2 = fabs(RHOU2);

  MFloat* const* const RESTRICT flux = ALIGNED_F(m_cells->flux);
  flux[CV->RHO_U][I] = RHOU2 * (UL + UR) + AbsRHO_U2 * (UL - UR) + PLR * surf0;
  flux[CV->RHO_V][I] = RHOU2 * (VL + VR) + AbsRHO_U2 * (VL - VR) + PLR * surf1;
  flux[CV->RHO_W][I] = RHOU2 * (WL + WR) + AbsRHO_U2 * (WL - WR) + PLR * surf2;
  flux[CV->RHO_E][I] = RHOU2 * (e0 + e1) + AbsRHO_U2 * (e0 - e1) + PLR * dxdtau;
  flux[CV->RHO][I] = RHOU;
}

void FvStructuredSolver3D::AusmDV(MFloat* QLeft, MFloat* QRight, const MInt dim, const MInt I) {
  const MFloat gamma = m_gamma;
  const MFloat gammaMinusOne = gamma - 1.0;
  const MFloat FgammaMinusOne = F1 / gammaMinusOne;

  const MFloat surf0 = m_cells->surfaceMetrics[dim * 3 + 0][I];
  const MFloat surf1 = m_cells->surfaceMetrics[dim * 3 + 1][I];
  const MFloat surf2 = m_cells->surfaceMetrics[dim * 3 + 2][I];

  const MFloat dxdtau = m_cells->dxt[dim][I];

  // left side
  const MFloat RHOL = QLeft[PV->RHO];
  const MFloat FRHOL = F1 / RHOL;
  MFloat UL = QLeft[PV->U];
  MFloat VL = QLeft[PV->V];
  MFloat WL = QLeft[PV->W];
  const MFloat PL = QLeft[PV->P];

  // right side
  const MFloat RHOR = QRight[PV->RHO];
  const MFloat FRHOR = F1 / RHOR;
  MFloat UR = QRight[PV->U];
  MFloat VR = QRight[PV->V];
  MFloat WR = QRight[PV->W];
  const MFloat PR = QRight[PV->P];

  const MFloat PLfRHOL = PL / RHOL;
  const MFloat PRfRHOR = PR / RHOR;
  const MFloat e0 = PLfRHOL * FgammaMinusOne + 0.5 * (POW2(UL) + POW2(VL) + POW2(WL)) + PLfRHOL;
  const MFloat e1 = PRfRHOR * FgammaMinusOne + 0.5 * (POW2(UR) + POW2(VR) + POW2(WR)) + PRfRHOR;


  // compute lenght of metric vector for normalization
  const MFloat DGRAD = sqrt(POW2(surf0) + POW2(surf1) + POW2(surf2));
  const MFloat FDGRAD = F1 / DGRAD;

  // scale by metric length to get velocity in the new basis (get normalized basis vectors)
  const MFloat UUL = ((UL * surf0 + VL * surf1 + WL * surf2) - dxdtau) * FDGRAD;


  const MFloat UUR = ((UR * surf0 + VR * surf1 + WR * surf2) - dxdtau) * FDGRAD;

  MFloat AL = FRHOL * PL;
  MFloat AR = FRHOR * PR;

  const MFloat FALR = 2.0 / (AL + AR);
  const MFloat ALPHAL = AL * FALR;
  const MFloat ALPHAR = AR * FALR;

  AL = sqrt(gamma * AL);
  AR = sqrt(gamma * AR);
  AL = mMax(AL, AR);
  AR = AL;

  const MFloat XMAL = UUL / AL;
  const MFloat XMAR = UUR / AR;

  AL = AL * DGRAD;
  AR = AR * DGRAD;

  const MFloat RHOAL = AL * RHOL;
  const MFloat RHOAR = AR * RHOR;

  const MInt IJK[2] = {1, m_nCells[1]};
  const MInt IP1 = I + IJK[dim];

  const MFloat FDV = 0.3;
  const MFloat DXDXEZ = m_cells->coordinates[0][IP1] - m_cells->coordinates[0][I];
  const MFloat DYDXEZ = m_cells->coordinates[1][IP1] - m_cells->coordinates[1][I];
  const MFloat DZDXEZ = m_cells->coordinates[2][IP1] - m_cells->coordinates[2][I];
  MFloat SV = 2.0 * DGRAD / (m_cells->cellJac[I] + m_cells->cellJac[IP1]) * (FDV + (F1 - FDV) * getPSI(I, dim));
  const MFloat SV1 = F0 * SV * DXDXEZ;
  const MFloat SV2 = F0 * SV * DYDXEZ;
  const MFloat SV3 = F0 * SV * DZDXEZ;

  const MFloat XMAL1 = mMin(F1, mMax(-F1, XMAL));
  const MFloat XMAR1 = mMin(F1, mMax(-F1, XMAR));

  MFloat FXMA = F1B2 * (XMAL1 + fabs(XMAL1));
  const MFloat XMALP = ALPHAL * (F1B4 * POW2(XMAL1 + F1) - FXMA) + FXMA + (mMax(F1, XMAL) - F1);
  FXMA = F1B2 * (XMAR1 - fabs(XMAR1));
  const MFloat XMARM = ALPHAR * (-F1B4 * POW2(XMAR1 - F1) - FXMA) + FXMA + (mMin(-F1, XMAR) + F1);

  const MFloat FLP = PL * ((F2 - XMAL1) * POW2(F1 + XMAL1));
  const MFloat FRP = PR * ((F2 + XMAR1) * POW2(F1 - XMAR1));
  const MFloat PLR = F1B4 * (FLP + FRP);

  const MFloat RHOUL = XMALP * RHOAL;
  const MFloat RHOUR = XMARM * RHOAR;
  const MFloat RHOU = RHOUL + RHOUR;
  const MFloat RHOU2 = F1B2 * RHOU;
  const MFloat ARHOU2 = fabs(RHOU2);

  const MFloat UUL2 = SV1 * UUL;
  const MFloat UUR2 = SV1 * UUR;
  UL = UL - UUL2;
  UR = UR - UUR2;
  const MFloat UUL3 = SV2 * UUL;
  const MFloat UUR3 = SV2 * UUR;
  VL = VL - UUL3;
  VR = VR - UUR3;
  const MFloat UUL4 = SV3 * UUL;
  const MFloat UUR4 = SV3 * UUR;
  WL = WL - UUL4;
  WR = WR - UUR4;

  MFloat* const* const RESTRICT flux = ALIGNED_F(m_cells->flux);

  flux[CV->RHO_U][I] = RHOU2 * (UL + UR) + ARHOU2 * (UL - UR) + PLR * surf0 + RHOUL * UUL2 + RHOUR * UUR2;
  flux[CV->RHO_V][I] = RHOU2 * (VL + VR) + ARHOU2 * (VL - VR) + PLR * surf1 + RHOUL * UUL3 + RHOUR * UUR3;
  flux[CV->RHO_W][I] = RHOU2 * (WL + WR) + ARHOU2 * (WL - WR) + PLR * surf2 + RHOUL * UUL4 + RHOUR * UUR4;
  flux[CV->RHO_E][I] = RHOU2 * (e0 + e1) + ARHOU2 * (e0 - e1) + PLR * dxdtau;
  flux[CV->RHO][I] = RHOU;
}

template <MInt noVars>
void FvStructuredSolver3D::Muscl_AusmLES() {
  TRACE();

  const MUint noCells = m_noCells;
  const MInt IJK[3] = {1, m_nCells[2], m_nCells[1] * m_nCells[2]};

  const MFloat* const RESTRICT x = ALIGNED_F(m_cells->coordinates[0]);
  const MFloat* const RESTRICT y = ALIGNED_F(m_cells->coordinates[1]);
  const MFloat* const RESTRICT z = ALIGNED_F(m_cells->coordinates[2]);
  const MFloat* const* const RESTRICT vars = ALIGNED_F(m_cells->pvariables);
  MFloat* const* const RESTRICT dss = ALIGNED_F(m_cells->dss);
  MFloat* const RESTRICT cellRhs = ALIGNED_MF(m_cells->rightHandSide[0]);
  MFloat* const* const RESTRICT flux = ALIGNED_F(m_cells->flux);

  const MFloat gamma = m_gamma;
  const MFloat gammaMinusOne = gamma - 1.0;
  const MFloat FgammaMinusOne = F1 / gammaMinusOne;



  if(m_movingGrid || !m_dsIsComputed) {
    for(MInt dim = 0; dim < nDim; dim++) {
      maia::parallelFor<true, nDim>(beginP0(), endM1(), [=](const MInt& i, const MInt& j, const MInt& k) {
        const MInt I = cellIndex(i, j, k);
        const MInt IP1 = I + IJK[dim];
        dss[dim][I] = sqrt(POW2(x[IP1] - x[I]) + POW2(y[IP1] - y[I]) + POW2(z[IP1] - z[I]));
      });
    }

    m_dsIsComputed = true;
  }


  for(MInt dim = 0; dim < nDim; dim++) {
    maia::parallelFor<true, nDim>(beginP1(), endM2(), [=](const MInt& i, const MInt& j, const MInt& k) {
      const MInt I = cellIndex(i, j, k);
      const MInt IP1 = I + IJK[dim];
      const MInt IM1 = I - IJK[dim];
      const MInt IP2 = I + 2 * IJK[dim];

      const MFloat DS = dss[dim][I];
      const MFloat DSM1 = dss[dim][IM1];
      const MFloat DSP1 = dss[dim][IP1];

      const MFloat DSP = DS / POW2(DSP1 + DS);
      const MFloat DSM = DS / POW2(DSM1 + DS);

      // unrolled the loop so the compiler
      // can optimize better
      const MFloat DQU = vars[PV->U][IP1] - vars[PV->U][I];
      const MFloat DQPU = vars[PV->U][IP2] - vars[PV->U][IP1];
      const MFloat DQMU = vars[PV->U][I] - vars[PV->U][IM1];
      const MFloat UL = vars[PV->U][I] + DSM * (DSM1 * DQU + DS * DQMU);
      const MFloat UR = vars[PV->U][IP1] - DSP * (DS * DQPU + DSP1 * DQU);

      const MFloat DQV = vars[PV->V][IP1] - vars[PV->V][I];
      const MFloat DQPV = vars[PV->V][IP2] - vars[PV->V][IP1];
      const MFloat DQMV = vars[PV->V][I] - vars[PV->V][IM1];
      const MFloat VL = vars[PV->V][I] + DSM * (DSM1 * DQV + DS * DQMV);
      const MFloat VR = vars[PV->V][IP1] - DSP * (DS * DQPV + DSP1 * DQV);

      const MFloat DQW = vars[PV->W][IP1] - vars[PV->W][I];
      const MFloat DQPW = vars[PV->W][IP2] - vars[PV->W][IP1];
      const MFloat DQMW = vars[PV->W][I] - vars[PV->W][IM1];
      const MFloat WL = vars[PV->W][I] + DSM * (DSM1 * DQW + DS * DQMW);
      const MFloat WR = vars[PV->W][IP1] - DSP * (DS * DQPW + DSP1 * DQW);

      const MFloat DQP = vars[PV->P][IP1] - vars[PV->P][I];
      const MFloat DQPP = vars[PV->P][IP2] - vars[PV->P][IP1];
      const MFloat DQMP = vars[PV->P][I] - vars[PV->P][IM1];
      const MFloat PL = vars[PV->P][I] + DSM * (DSM1 * DQP + DS * DQMP);
      const MFloat PR = vars[PV->P][IP1] - DSP * (DS * DQPP + DSP1 * DQP);

      const MFloat DQRHO = vars[PV->RHO][IP1] - vars[PV->RHO][I];
      const MFloat DQPRHO = vars[PV->RHO][IP2] - vars[PV->RHO][IP1];
      const MFloat DQMRHO = vars[PV->RHO][I] - vars[PV->RHO][IM1];
      const MFloat RHOL = vars[PV->RHO][I] + DSM * (DSM1 * DQRHO + DS * DQMRHO);
      const MFloat RHOR = vars[PV->RHO][IP1] - DSP * (DS * DQPRHO + DSP1 * DQRHO);

      const MFloat surf0 = m_cells->surfaceMetrics[dim * 3 + 0][I];
      const MFloat surf1 = m_cells->surfaceMetrics[dim * 3 + 1][I];
      const MFloat surf2 = m_cells->surfaceMetrics[dim * 3 + 2][I];
      const MFloat dxdtau = m_cells->dxt[dim][I];

      // compute length of metric vector for normalization
      const MFloat metricLength = sqrt(POW2(surf0) + POW2(surf1) + POW2(surf2));
      const MFloat fMetricLength = F1 / metricLength;

      // scale by metric length to get velocity in the new basis (get normalized basis vectors)
      const MFloat UUL = ((UL * surf0 + VL * surf1 + WL * surf2) - dxdtau) * fMetricLength;


      const MFloat UUR = ((UR * surf0 + VR * surf1 + WR * surf2) - dxdtau) * fMetricLength;


      // speed of sound
      const MFloat AL = sqrt(gamma * max(m_eps, (PL / max(m_eps, RHOL))));
      const MFloat AR = sqrt(gamma * max(m_eps, (PR / max(m_eps, RHOR))));

      const MFloat MAL = UUL / AL;
      const MFloat MAR = UUR / AR;

      const MFloat MALR = F1B2 * (MAL + MAR);
      const MFloat PLR = PL * (F1B2 + m_chi * MAL) + PR * (F1B2 - m_chi * MAR);

      const MFloat RHO_AL = RHOL * AL;
      const MFloat RHO_AR = RHOR * AR;

      const MFloat PLfRHOL = PL / RHOL;
      const MFloat PRfRHOR = PR / RHOR;

      const MFloat e0 = PLfRHOL * FgammaMinusOne + 0.5 * (POW2(UL) + POW2(VL) + POW2(WL)) + PLfRHOL;
      const MFloat e1 = PRfRHOR * FgammaMinusOne + 0.5 * (POW2(UR) + POW2(VR) + POW2(WR)) + PRfRHOR;

      const MFloat RHOU = F1B2 * (MALR * (RHO_AL + RHO_AR) + fabs(MALR) * (RHO_AL - RHO_AR)) * metricLength;
      const MFloat RHOU2 = F1B2 * RHOU;
      // multiply by metric length to take surface area into account
      const MFloat AbsRHO_U2 = fabs(RHOU2);

      flux[CV->RHO_U][I] = RHOU2 * (UL + UR) + AbsRHO_U2 * (UL - UR) + PLR * surf0;
      flux[CV->RHO_V][I] = RHOU2 * (VL + VR) + AbsRHO_U2 * (VL - VR) + PLR * surf1;
      flux[CV->RHO_W][I] = RHOU2 * (WL + WR) + AbsRHO_U2 * (WL - WR) + PLR * surf2;
      flux[CV->RHO_E][I] = RHOU2 * (e0 + e1) + AbsRHO_U2 * (e0 - e1) + PLR * dxdtau;
      flux[CV->RHO][I] = RHOU;
    });

    // FLUX BALANCE
    for(MUint v = 0; v < noVars; v++) {
      maia::parallelFor<true, nDim>(beginP2(), endM2(), [=](const MInt& i, const MInt& j, const MInt& k) {
        const MInt I = cellIndex(i, j, k);
        const MInt IM1 = I - IJK[dim];
        const MUint offset = v * noCells;
        MFloat* const RESTRICT rhs = ALIGNED_F(cellRhs + offset);
        rhs[I] += flux[v][IM1] - flux[v][I];
      });
    }
  }
}

template void FvStructuredSolver3D::Muscl_AusmLES<5>();
template void FvStructuredSolver3D::Muscl_AusmLES<6>();
template void FvStructuredSolver3D::Muscl_AusmLES<7>();

template <MInt noVars>
void FvStructuredSolver3D::Muscl_AusmLES_PTHRC() {
  TRACE();

  // stencil identifier
  const MInt IJK[3] = {1, m_nCells[2], m_nCells[1] * m_nCells[2]};

  const MFloat* const RESTRICT x = ALIGNED_F(m_cells->coordinates[0]);
  const MFloat* const RESTRICT y = ALIGNED_F(m_cells->coordinates[1]);
  const MFloat* const RESTRICT z = ALIGNED_F(m_cells->coordinates[2]);
  const MFloat* const* const RESTRICT vars = ALIGNED_F(m_cells->pvariables);
  const MFloat* const RESTRICT p = ALIGNED_F(m_cells->pvariables[PV->P]);
  MFloat* const* const RESTRICT dss = ALIGNED_F(m_cells->dss);
  MFloat* const RESTRICT cellRhs = ALIGNED_MF(m_cells->rightHandSide[0]);
  MFloat* const* const RESTRICT flux = ALIGNED_F(m_cells->flux);

  const MUint noCells = m_noCells;
  const MFloat gamma = m_gamma;
  const MFloat gammaMinusOne = gamma - 1.0;
  const MFloat FgammaMinusOne = F1 / gammaMinusOne;

  const MInt noCellsI = m_nCells[2] - 2;
  const MInt noCellsJ = m_nCells[1] - 2;
  const MInt noCellsK = m_nCells[0] - 2;

  const MInt noCellsIP1 = m_nCells[2] - 1;
  const MInt noCellsJP1 = m_nCells[1] - 1;
  const MInt noCellsKP1 = m_nCells[0] - 1;

  if(m_movingGrid || !m_dsIsComputed) {
    for(MInt dim = 0; dim < nDim; dim++) {
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for(MInt k = 0; k < noCellsKP1; k++) {
        for(MInt j = 0; j < noCellsJP1; j++) {
          for(MInt i = 0; i < noCellsIP1; i++) {
            const MInt I = cellIndex(i, j, k);
            const MInt IP1 = I + IJK[dim];
            dss[dim][I] = sqrt(POW2(x[IP1] - x[I]) + POW2(y[IP1] - y[I]) + POW2(z[IP1] - z[I]));
          }
        }
      }
    }

    m_dsIsComputed = true;
  }

  for(MInt dim = 0; dim < nDim; dim++) {
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(MInt k = 1; k < noCellsK; k++) {
      for(MInt j = 1; j < noCellsJ; j++) {
#if defined(MAIA_INTEL_COMPILER)
#pragma ivdep
#pragma vector always
#endif
        for(MInt i = 1; i < noCellsI; i++) {
          const MInt I = cellIndex(i, j, k);
          const MInt IP1 = I + IJK[dim];
          const MInt IM1 = I - IJK[dim];
          const MInt IP2 = I + 2 * IJK[dim];

          const MFloat DS = dss[dim][I];
          const MFloat DSM1 = dss[dim][IM1];
          const MFloat DSP1 = dss[dim][IP1];

          const MFloat DSP = DS / POW2(DSP1 + DS);
          const MFloat DSM = DS / POW2(DSM1 + DS);

          // unrolled the loop so the compiler
          // can optimize better
          const MFloat DQU = vars[PV->U][IP1] - vars[PV->U][I];
          const MFloat DQPU = vars[PV->U][IP2] - vars[PV->U][IP1];
          const MFloat DQMU = vars[PV->U][I] - vars[PV->U][IM1];
          const MFloat UL = vars[PV->U][I] + DSM * (DSM1 * DQU + DS * DQMU);
          const MFloat UR = vars[PV->U][IP1] - DSP * (DS * DQPU + DSP1 * DQU);

          const MFloat DQV = vars[PV->V][IP1] - vars[PV->V][I];
          const MFloat DQPV = vars[PV->V][IP2] - vars[PV->V][IP1];
          const MFloat DQMV = vars[PV->V][I] - vars[PV->V][IM1];
          const MFloat VL = vars[PV->V][I] + DSM * (DSM1 * DQV + DS * DQMV);
          const MFloat VR = vars[PV->V][IP1] - DSP * (DS * DQPV + DSP1 * DQV);

          const MFloat DQW = vars[PV->W][IP1] - vars[PV->W][I];
          const MFloat DQPW = vars[PV->W][IP2] - vars[PV->W][IP1];
          const MFloat DQMW = vars[PV->W][I] - vars[PV->W][IM1];
          const MFloat WL = vars[PV->W][I] + DSM * (DSM1 * DQW + DS * DQMW);
          const MFloat WR = vars[PV->W][IP1] - DSP * (DS * DQPW + DSP1 * DQW);

          const MFloat DQP = vars[PV->P][IP1] - vars[PV->P][I];
          const MFloat DQPP = vars[PV->P][IP2] - vars[PV->P][IP1];
          const MFloat DQMP = vars[PV->P][I] - vars[PV->P][IM1];
          const MFloat PL = vars[PV->P][I] + DSM * (DSM1 * DQP + DS * DQMP);
          const MFloat PR = vars[PV->P][IP1] - DSP * (DS * DQPP + DSP1 * DQP);

          const MFloat DQRHO = vars[PV->RHO][IP1] - vars[PV->RHO][I];
          const MFloat DQPRHO = vars[PV->RHO][IP2] - vars[PV->RHO][IP1];
          const MFloat DQMRHO = vars[PV->RHO][I] - vars[PV->RHO][IM1];
          const MFloat RHOL = vars[PV->RHO][I] + DSM * (DSM1 * DQRHO + DS * DQMRHO);
          const MFloat RHOR = vars[PV->RHO][IP1] - DSP * (DS * DQPRHO + DSP1 * DQRHO);

          const MFloat surf0 = m_cells->surfaceMetrics[dim * 3 + 0][I];
          const MFloat surf1 = m_cells->surfaceMetrics[dim * 3 + 1][I];
          const MFloat surf2 = m_cells->surfaceMetrics[dim * 3 + 2][I];
          const MFloat dxdtau = m_cells->dxt[dim][I];

          // compute lenght of metric vector for normalization
          const MFloat metricLength = sqrt(POW2(surf0) + POW2(surf1) + POW2(surf2));
          const MFloat fMetricLength = F1 / metricLength;

          // scale by metric length to get velocity in the new basis (get normalized basis vectors)
          const MFloat UUL = ((UL * surf0 + VL * surf1 + WL * surf2) - dxdtau) * fMetricLength;


          const MFloat UUR = ((UR * surf0 + VR * surf1 + WR * surf2) - dxdtau) * fMetricLength;


          // speed of sound
          const MFloat AL = sqrt(gamma * max(m_eps, (PL / max(m_eps, RHOL))));
          const MFloat AR = sqrt(gamma * max(m_eps, (PR / max(m_eps, RHOR))));

          const MFloat MAL = UUL / AL;
          const MFloat MAR = UUR / AR;

          const MFloat MALR = F1B2 * (MAL + MAR);

          // 4th order pressure damping
          const MInt IPJK = getCellIdfromCell(I, 1, 0, 0);
          const MInt IMJK = getCellIdfromCell(I, -1, 0, 0);
          const MInt IP2JK = getCellIdfromCell(I, 2, 0, 0);
          const MInt IM2JK = getCellIdfromCell(I, -2, 0, 0);

          const MInt IJPK = getCellIdfromCell(I, 0, 1, 0);
          const MInt IJMK = getCellIdfromCell(I, 0, -1, 0);
          const MInt IJP2K = getCellIdfromCell(I, 0, 2, 0);
          const MInt IJM2K = getCellIdfromCell(I, 0, -2, 0);

          const MInt IJKP = getCellIdfromCell(I, 0, 0, 1);
          const MInt IJKM = getCellIdfromCell(I, 0, 0, -1);
          const MInt IJKP2 = getCellIdfromCell(I, 0, 0, 2);
          const MInt IJKM2 = getCellIdfromCell(I, 0, 0, -2);

          const MFloat p4I4 = F4 * (p[IPJK] + p[IMJK]) - F6 * (p[I]) - p[IP2JK] - p[IM2JK];
          const MFloat p4J4 = F4 * (p[IJPK] + p[IJMK]) - F6 * (p[I]) - p[IJP2K] - p[IJM2K];
          const MFloat p4K4 = F4 * (p[IJKP] + p[IJKM]) - F6 * (p[I]) - p[IJKP2] - p[IJKM2];

          const MFloat pfac = fabs(p4I4) + fabs(p4J4) + fabs(p4K4);
          const MFloat facl = 1.0 / 1.3 * pfac;
          const MFloat fac = min(1.0 / 128.0, facl);

          const MFloat PLR = PL * (F1B2 + fac * MAL) + PR * (F1B2 - fac * MAR);

          const MFloat RHO_AL = RHOL * AL;
          const MFloat RHO_AR = RHOR * AR;

          const MFloat PLfRHOL = PL / RHOL;
          const MFloat PRfRHOR = PR / RHOR;

          const MFloat e0 = PLfRHOL * FgammaMinusOne + 0.5 * (POW2(UL) + POW2(VL) + POW2(WL)) + PLfRHOL;
          const MFloat e1 = PRfRHOR * FgammaMinusOne + 0.5 * (POW2(UR) + POW2(VR) + POW2(WR)) + PRfRHOR;

          const MFloat RHOU = F1B2 * (MALR * (RHO_AL + RHO_AR) + fabs(MALR) * (RHO_AL - RHO_AR)) * metricLength;
          const MFloat RHOU2 = F1B2 * RHOU;
          // multiply by metric length to take surface area into account
          const MFloat AbsRHO_U2 = fabs(RHOU2);

          flux[CV->RHO_U][I] = RHOU2 * (UL + UR) + AbsRHO_U2 * (UL - UR) + PLR * surf0;
          flux[CV->RHO_V][I] = RHOU2 * (VL + VR) + AbsRHO_U2 * (VL - VR) + PLR * surf1;
          flux[CV->RHO_W][I] = RHOU2 * (WL + WR) + AbsRHO_U2 * (WL - WR) + PLR * surf2;
          flux[CV->RHO_E][I] = RHOU2 * (e0 + e1) + AbsRHO_U2 * (e0 - e1) + PLR * dxdtau;
          flux[CV->RHO][I] = RHOU;
        }
      }
    }

    // FLUX BALANCE
    for(MUint v = 0; v < noVars; v++) {
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for(MInt k = m_noGhostLayers; k < m_nCells[0] - m_noGhostLayers; k++) {
        for(MInt j = m_noGhostLayers; j < m_nCells[1] - m_noGhostLayers; j++) {
#if defined(MAIA_INTEL_COMPILER)
#pragma ivdep
#pragma vector always
#endif
          for(MInt i = m_noGhostLayers; i < m_nCells[2] - m_noGhostLayers; i++) {
            const MInt I = cellIndex(i, j, k);
            const MInt IM1 = I - IJK[dim];
            const MUint offset = v * noCells;
            MFloat* const RESTRICT rhs = ALIGNED_F(cellRhs + offset);
            rhs[I] += flux[v][IM1] - flux[v][I];
          }
        }
      }
    }
  }
}

template void FvStructuredSolver3D::Muscl_AusmLES_PTHRC<5>();
template void FvStructuredSolver3D::Muscl_AusmLES_PTHRC<6>();
template void FvStructuredSolver3D::Muscl_AusmLES_PTHRC<7>();


template <MInt noVars>
void FvStructuredSolver3D::Muscl_AusmDV() {
  TRACE();

  // stencil identifier
  const MInt IJK[3] = {1, m_nCells[2], m_nCells[1] * m_nCells[2]};

  const MFloat* const RESTRICT x = ALIGNED_F(m_cells->coordinates[0]);
  const MFloat* const RESTRICT y = ALIGNED_F(m_cells->coordinates[1]);
  const MFloat* const RESTRICT z = ALIGNED_F(m_cells->coordinates[2]);
  const MFloat* const* const RESTRICT vars = ALIGNED_F(m_cells->pvariables);
  MFloat* const* const RESTRICT dss = ALIGNED_F(m_cells->dss);
  MFloat* const* const RESTRICT flux = ALIGNED_F(m_cells->flux);
  MFloat* const RESTRICT cellRhs = ALIGNED_MF(m_cells->rightHandSide[0]);

  const MUint noCells = m_noCells;
  const MFloat gamma = m_gamma;
  const MFloat gammaMinusOne = gamma - 1.0;
  const MFloat FgammaMinusOne = F1 / gammaMinusOne;

  const MInt noCellsI = m_nCells[2] - 2;
  const MInt noCellsJ = m_nCells[1] - 2;
  const MInt noCellsK = m_nCells[0] - 2;

  const MInt noCellsIP1 = m_nCells[2] - 1;
  const MInt noCellsJP1 = m_nCells[1] - 1;
  const MInt noCellsKP1 = m_nCells[0] - 1;

  if(m_movingGrid || !m_dsIsComputed) {
    for(MInt dim = 0; dim < nDim; dim++) {
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for(MInt k = 0; k < noCellsKP1; k++) {
        for(MInt j = 0; j < noCellsJP1; j++) {
          for(MInt i = 0; i < noCellsIP1; i++) {
            const MInt I = cellIndex(i, j, k);
            const MInt IP1 = I + IJK[dim];
            dss[dim][I] = sqrt(POW2(x[IP1] - x[I]) + POW2(y[IP1] - y[I]) + POW2(z[IP1] - z[I]));
          }
        }
      }
    }

    m_dsIsComputed = true;
  }


  for(MInt dim = 0; dim < nDim; dim++) {
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(MInt k = 1; k < noCellsK; k++) {
      for(MInt j = 1; j < noCellsJ; j++) {
#if defined(MAIA_INTEL_COMPILER)
#pragma ivdep
#pragma vector always
#endif
        for(MInt i = 1; i < noCellsI; i++) {
          const MInt I = cellIndex(i, j, k);
          const MInt IP1 = I + IJK[dim];
          const MInt IM1 = I - IJK[dim];
          const MInt IP2 = I + 2 * IJK[dim];

          const MFloat DS = dss[dim][I];
          const MFloat DSM1 = dss[dim][IM1];
          const MFloat DSP1 = dss[dim][IP1];

          const MFloat DSP = DS / POW2(DSP1 + DS);
          const MFloat DSM = DS / POW2(DSM1 + DS);

          // unrolled the loop so the compiler
          // can optimize better
          const MFloat DQU = vars[PV->U][IP1] - vars[PV->U][I];
          const MFloat DQPU = vars[PV->U][IP2] - vars[PV->U][IP1];
          const MFloat DQMU = vars[PV->U][I] - vars[PV->U][IM1];
          MFloat UL = vars[PV->U][I] + DSM * (DSM1 * DQU + DS * DQMU);
          MFloat UR = vars[PV->U][IP1] - DSP * (DS * DQPU + DSP1 * DQU);

          const MFloat DQV = vars[PV->V][IP1] - vars[PV->V][I];
          const MFloat DQPV = vars[PV->V][IP2] - vars[PV->V][IP1];
          const MFloat DQMV = vars[PV->V][I] - vars[PV->V][IM1];
          MFloat VL = vars[PV->V][I] + DSM * (DSM1 * DQV + DS * DQMV);
          MFloat VR = vars[PV->V][IP1] - DSP * (DS * DQPV + DSP1 * DQV);

          const MFloat DQW = vars[PV->W][IP1] - vars[PV->W][I];
          const MFloat DQPW = vars[PV->W][IP2] - vars[PV->W][IP1];
          const MFloat DQMW = vars[PV->W][I] - vars[PV->W][IM1];
          MFloat WL = vars[PV->W][I] + DSM * (DSM1 * DQW + DS * DQMW);
          MFloat WR = vars[PV->W][IP1] - DSP * (DS * DQPW + DSP1 * DQW);

          const MFloat DQP = vars[PV->P][IP1] - vars[PV->P][I];
          const MFloat DQPP = vars[PV->P][IP2] - vars[PV->P][IP1];
          const MFloat DQMP = vars[PV->P][I] - vars[PV->P][IM1];
          const MFloat PL = vars[PV->P][I] + DSM * (DSM1 * DQP + DS * DQMP);
          const MFloat PR = vars[PV->P][IP1] - DSP * (DS * DQPP + DSP1 * DQP);

          const MFloat DQRHO = vars[PV->RHO][IP1] - vars[PV->RHO][I];
          const MFloat DQPRHO = vars[PV->RHO][IP2] - vars[PV->RHO][IP1];
          const MFloat DQMRHO = vars[PV->RHO][I] - vars[PV->RHO][IM1];
          const MFloat RHOL = vars[PV->RHO][I] + DSM * (DSM1 * DQRHO + DS * DQMRHO);
          const MFloat RHOR = vars[PV->RHO][IP1] - DSP * (DS * DQPRHO + DSP1 * DQRHO);

          const MFloat surf0 = m_cells->surfaceMetrics[dim * 3 + 0][I];
          const MFloat surf1 = m_cells->surfaceMetrics[dim * 3 + 1][I];
          const MFloat surf2 = m_cells->surfaceMetrics[dim * 3 + 2][I];
          const MFloat dxdtau = m_cells->dxt[dim][I];

          const MFloat FRHOL = F1 / RHOL;
          const MFloat FRHOR = F1 / RHOR;

          const MFloat PLfRHOL = PL / RHOL;
          const MFloat PRfRHOR = PR / RHOR;
          const MFloat e0 = PLfRHOL * FgammaMinusOne + 0.5 * (POW2(UL) + POW2(VL) + POW2(WL)) + PLfRHOL;
          const MFloat e1 = PRfRHOR * FgammaMinusOne + 0.5 * (POW2(UR) + POW2(VR) + POW2(WR)) + PRfRHOR;


          // compute lenght of metric vector for normalization
          const MFloat DGRAD = sqrt(POW2(surf0) + POW2(surf1) + POW2(surf2));
          const MFloat FDGRAD = F1 / DGRAD;

          // scale by metric length to get velocity in the new basis (get normalized basis vectors)
          const MFloat UUL = ((UL * surf0 + VL * surf1 + WL * surf2) - dxdtau) * FDGRAD;


          const MFloat UUR = ((UR * surf0 + VR * surf1 + WR * surf2) - dxdtau) * FDGRAD;

          MFloat AL = FRHOL * PL;
          MFloat AR = FRHOR * PR;

          const MFloat FALR = 2.0 / (AL + AR);
          const MFloat ALPHAL = AL * FALR;
          const MFloat ALPHAR = AR * FALR;

          AL = sqrt(gamma * AL);
          AR = sqrt(gamma * AR);
          AL = mMax(AL, AR);
          AR = AL;

          const MFloat XMAL = UUL / AL;
          const MFloat XMAR = UUR / AR;

          AL = AL * DGRAD;
          AR = AR * DGRAD;

          const MFloat RHOAL = AL * RHOL;
          const MFloat RHOAR = AR * RHOR;

          const MFloat FDV = 0.3;
          const MFloat DXDXEZ = m_cells->coordinates[0][IP1] - m_cells->coordinates[0][I];
          const MFloat DYDXEZ = m_cells->coordinates[1][IP1] - m_cells->coordinates[1][I];
          const MFloat DZDXEZ = m_cells->coordinates[2][IP1] - m_cells->coordinates[2][I];
          MFloat SV = 2.0 * DGRAD / (m_cells->cellJac[I] + m_cells->cellJac[IP1]) * (FDV + (F1 - FDV) * getPSI(I, dim));
          const MFloat SV1 = F0 * SV * DXDXEZ;
          const MFloat SV2 = F0 * SV * DYDXEZ;
          const MFloat SV3 = F0 * SV * DZDXEZ;

          const MFloat XMAL1 = mMin(F1, mMax(-F1, XMAL));
          const MFloat XMAR1 = mMin(F1, mMax(-F1, XMAR));

          MFloat FXMA = F1B2 * (XMAL1 + fabs(XMAL1));
          const MFloat XMALP = ALPHAL * (F1B4 * POW2(XMAL1 + F1) - FXMA) + FXMA + (mMax(F1, XMAL) - F1);
          FXMA = F1B2 * (XMAR1 - fabs(XMAR1));
          const MFloat XMARM = ALPHAR * (-F1B4 * POW2(XMAR1 - F1) - FXMA) + FXMA + (mMin(-F1, XMAR) + F1);

          const MFloat FLP = PL * ((F2 - XMAL1) * POW2(F1 + XMAL1));
          const MFloat FRP = PR * ((F2 + XMAR1) * POW2(F1 - XMAR1));
          const MFloat PLR = F1B4 * (FLP + FRP);

          const MFloat RHOUL = XMALP * RHOAL;
          const MFloat RHOUR = XMARM * RHOAR;
          const MFloat RHOU = RHOUL + RHOUR;
          const MFloat RHOU2 = F1B2 * RHOU;
          const MFloat ARHOU2 = fabs(RHOU2);

          const MFloat UUL2 = SV1 * UUL;
          const MFloat UUR2 = SV1 * UUR;
          UL = UL - UUL2;
          UR = UR - UUR2;
          const MFloat UUL3 = SV2 * UUL;
          const MFloat UUR3 = SV2 * UUR;
          VL = VL - UUL3;
          VR = VR - UUR3;
          const MFloat UUL4 = SV3 * UUL;
          const MFloat UUR4 = SV3 * UUR;
          WL = WL - UUL4;
          WR = WR - UUR4;

          flux[CV->RHO_U][I] = RHOU2 * (UL + UR) + ARHOU2 * (UL - UR) + PLR * surf0 + RHOUL * UUL2 + RHOUR * UUR2;
          flux[CV->RHO_V][I] = RHOU2 * (VL + VR) + ARHOU2 * (VL - VR) + PLR * surf1 + RHOUL * UUL3 + RHOUR * UUR3;
          flux[CV->RHO_W][I] = RHOU2 * (WL + WR) + ARHOU2 * (WL - WR) + PLR * surf2 + RHOUL * UUL4 + RHOUR * UUR4;
          flux[CV->RHO_E][I] = RHOU2 * (e0 + e1) + ARHOU2 * (e0 - e1) + PLR * dxdtau;
          flux[CV->RHO][I] = RHOU;
        }
      }
    }

    // FLUX BALANCE
    for(MUint v = 0; v < noVars; v++) {
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for(MInt k = m_noGhostLayers; k < m_nCells[0] - m_noGhostLayers; k++) {
        for(MInt j = m_noGhostLayers; j < m_nCells[1] - m_noGhostLayers; j++) {
#if defined(MAIA_INTEL_COMPILER)
#pragma ivdep
#pragma vector always
#endif
          for(MInt i = m_noGhostLayers; i < m_nCells[2] - m_noGhostLayers; i++) {
            const MInt I = cellIndex(i, j, k);
            const MInt IM1 = I - IJK[dim];
            const MUint offset = v * noCells;
            MFloat* const RESTRICT rhs = ALIGNED_F(cellRhs + offset);
            rhs[I] += flux[v][IM1] - flux[v][I];
          }
        }
      }
    }
  }
}

template void FvStructuredSolver3D::Muscl_AusmDV<5>();
template void FvStructuredSolver3D::Muscl_AusmDV<6>();
template void FvStructuredSolver3D::Muscl_AusmDV<7>();


template <FvStructuredSolver3D::fluxmethod ausm, MInt noVars>
void FvStructuredSolver3D::MusclStretched_() {
  TRACE();

  // stencil identifier
  const MUint noCells = m_noCells;
  const MInt IJK[3] = {1, m_nCells[2], m_nCells[1] * m_nCells[2]};
  const MFloat* const RESTRICT cellVariables = ALIGNED_F(m_cells->pvariables[0]);
  const MFloat* const RESTRICT cellLength = ALIGNED_F(m_cells->cellLength[0]);
  MFloat* const RESTRICT cellRhs = ALIGNED_MF(m_cells->rightHandSide[0]);
  MFloat* const RESTRICT qleft = ALIGNED_MF(m_QLeft);
  MFloat* const RESTRICT qright = ALIGNED_MF(m_QRight);
  MFloat* const* const RESTRICT flux = ALIGNED_F(m_cells->flux);
  /////////IMPORTANT PARAMETER
  // MFloat epsi=F1;
  const MFloat phi = F1;
  const MFloat kappa = F0; // F1B3;
  /////////END IMPORTANT PARAMETER
  for(MInt dim = 0; dim < nDim; dim++) {
    const MUint dimOffset = dim * m_noCells;
    const MFloat* const RESTRICT length = ALIGNED_F(cellLength + dimOffset);

    for(MInt k = m_noGhostLayers - 1; k < m_nCells[0] - m_noGhostLayers; k++) {
      for(MInt j = m_noGhostLayers - 1; j < m_nCells[1] - m_noGhostLayers; j++) {
        for(MInt i = m_noGhostLayers - 1; i < m_nCells[2] - m_noGhostLayers; i++) {
          const MInt I = cellIndex(i, j, k);
          const MInt IP1 = I + IJK[dim];
          const MInt IM1 = I - IJK[dim];
          const MInt IP2 = I + 2 * IJK[dim];

          const MFloat rp = (length[I] + length[IP1]) / (F2 * length[I]);
          const MFloat rm = (length[I] + length[IM1]) / (F2 * length[I]);
          const MFloat f = phi / (F2 * (rp + rm));
          const MFloat f1 = (rm + kappa * phi) / rp;
          const MFloat f2 = (rp - kappa * phi) / rm;

          const MFloat rp1 = (length[IP1] + length[IP2]) / (F2 * length[IP1]);
          const MFloat rm1 = (length[IP1] + length[I]) / (F2 * length[IP1]);
          const MFloat fa = phi / (F2 * (rp1 + rm1));
          const MFloat fb = (rm1 - kappa * phi) / rp1;
          const MFloat fc = (rp1 + kappa * phi) / rm1;

          for(MUint v = 0; v < noVars; v++) {
            const MUint offset = v * m_noCells;
            const MFloat* const RESTRICT vars = ALIGNED_F(cellVariables + offset);
            // left variables
            const MFloat DQ = (vars[IP1] - vars[I]);
            const MFloat DQM1 = (vars[I] - vars[IM1]);
            qleft[v] = vars[I] + f * (f1 * DQ + f2 * DQM1);

            // right variables
            const MFloat DQP1 = (vars[IP2] - vars[IP1]);
            const MFloat DQ1 = (vars[IP1] - vars[I]);
            qright[v] = vars[IP1] - fa * (fb * DQP1 + fc * DQ1);
          }

          (this->*ausm)(m_QLeft, m_QRight, dim, I);
        }
      }
    }

    // FLUX BALANCE
    for(MUint v = 0; v < noVars; v++) {
      for(MInt k = m_noGhostLayers; k < m_nCells[0] - m_noGhostLayers; k++) {
        for(MInt j = m_noGhostLayers; j < m_nCells[1] - m_noGhostLayers; j++) {
          for(MInt i = m_noGhostLayers; i < m_nCells[2] - m_noGhostLayers; i++) {
            const MInt I = cellIndex(i, j, k);
            const MInt IM1 = I - IJK[dim];
            const MUint offset = v * noCells;
            MFloat* const RESTRICT rhs = ALIGNED_F(cellRhs + offset);
            rhs[I] += flux[v][IM1] - flux[v][I];
          }
        }
      }
    }
  }
}
// standard Ausm
template void FvStructuredSolver3D::MusclStretched_<&FvStructuredSolver3D::AusmLES, 5>();
template void FvStructuredSolver3D::MusclStretched_<&FvStructuredSolver3D::AusmLES, 6>();
template void FvStructuredSolver3D::MusclStretched_<&FvStructuredSolver3D::AusmLES, 7>();
// pthrc
template void FvStructuredSolver3D::MusclStretched_<&FvStructuredSolver3D::AusmLES_PTHRC, 5>();
template void FvStructuredSolver3D::MusclStretched_<&FvStructuredSolver3D::AusmLES_PTHRC, 6>();
template void FvStructuredSolver3D::MusclStretched_<&FvStructuredSolver3D::AusmLES_PTHRC, 7>();
// AusmDV
template void FvStructuredSolver3D::MusclStretched_<&FvStructuredSolver3D::AusmDV, 5>();
template void FvStructuredSolver3D::MusclStretched_<&FvStructuredSolver3D::AusmDV, 6>();
template void FvStructuredSolver3D::MusclStretched_<&FvStructuredSolver3D::AusmDV, 7>();


void FvStructuredSolver3D::Muscl(MInt NotUsed(timerId)) {
  TRACE();


  if(m_movingGrid) {
    RECORD_TIMER_START(m_timers[Timers::MovingGrid]);
    if(m_RKStep == 0) {
      RECORD_TIMER_START(m_timers[Timers::MGSaveGrid]);
      m_grid->saveGrid();
      m_grid->saveCellJacobian();
      RECORD_TIMER_STOP(m_timers[Timers::MGSaveGrid]);
    }

    RECORD_TIMER_START(m_timers[Timers::MGMoveGrid]);
    moveGrid(false, false);
    RECORD_TIMER_STOP(m_timers[Timers::MGMoveGrid]);

    // compute the volume fluxes
    RECORD_TIMER_START(m_timers[Timers::MGVolumeFlux]);
    m_grid->computeDxt(m_timeStep, m_RKalpha, m_RKStep);
    RECORD_TIMER_STOP(m_timers[Timers::MGVolumeFlux]);
    RECORD_TIMER_STOP(m_timers[Timers::MovingGrid]);
  }

  if(m_bodyForce) {
    applyBodyForce(false, false);
  }

  RECORD_TIMER_START(m_timers[Timers::ConvectiveFlux]);
  (this->*reconstructSurfaceData)();
  RECORD_TIMER_STOP(m_timers[Timers::ConvectiveFlux]);

  if(m_useSandpaperTrip) {
    RECORD_TIMER_START(m_timers[Timers::SandpaperTrip]);
    if(m_tripAirfoil) {
      applySandpaperTripAirfoil();
    } else {
      applySandpaperTrip();
    }
    RECORD_TIMER_STOP(m_timers[Timers::SandpaperTrip]);
  }
}

void FvStructuredSolver3D::Ausm() {
  // Ausm routines have been moved and are called from inside Muscl (better performance)
}

void FvStructuredSolver3D::computeDomainWidth() {
  // if we only compute cd,cl for a part of the domain (m_auxDataCoordinateLimits == true)
  // only compute the average with the width of this section
  if(m_auxDataCoordinateLimits) {
    m_globalDomainWidth = fabs(m_auxDataLimits[3] - m_auxDataLimits[2]);
  } else {
    MFloat minCoordinate = F0, maxCoordinate = F0, minCoordinateGlobal = F0, maxCoordinateGlobal = F0;
    MInt lowPoint = -1, highPoint = -1;

    if(m_forceAveragingDir == 0) {
      highPoint = getPointIdFromCell(m_nCells[2] - m_noGhostLayers, 0, 0);
      lowPoint = getPointIdFromCell(m_noGhostLayers, 0, 0);
    } else if(m_forceAveragingDir == 1) {
      highPoint = getPointIdFromCell(0, m_nCells[1] - m_noGhostLayers, 0);
      lowPoint = getPointIdFromCell(0, m_noGhostLayers, 0);
    } else {
      highPoint = getPointIdFromCell(m_noGhostLayers, m_noGhostLayers, m_nCells[0] - m_noGhostLayers);
      lowPoint = getPointIdFromCell(m_noGhostLayers, m_noGhostLayers, m_noGhostLayers);
    }

    minCoordinate = m_grid->m_coordinates[m_forceAveragingDir][lowPoint];
    maxCoordinate = m_grid->m_coordinates[m_forceAveragingDir][highPoint];
    MPI_Allreduce(&minCoordinate, &minCoordinateGlobal, 1, MPI_DOUBLE, MPI_MIN, m_StructuredComm, AT_, "minCoordinate",
                  "minCoordinateGlobal");
    MPI_Allreduce(&maxCoordinate, &maxCoordinateGlobal, 1, MPI_DOUBLE, MPI_MAX, m_StructuredComm, AT_, "maxCoordinate",
                  "maxCoordinateGlobal");
    m_globalDomainWidth = fabs(maxCoordinateGlobal - minCoordinateGlobal);
    m_log << "Global domain width: " << m_globalDomainWidth << endl;
  }
}

void FvStructuredSolver3D::initSandpaperTrip() {
  TRACE();

  m_tripNoTrips = Context::propertyLength("tripDelta1", m_solverId);

  mAlloc(m_tripDelta1, m_tripNoTrips, "m_tripDelta1", 1.0, AT_);
  mAlloc(m_tripXOrigin, m_tripNoTrips, "m_tripXOrigin", 30.0, AT_);
  mAlloc(m_tripXLength, m_tripNoTrips, "m_tripXLength", 4.0, AT_);
  mAlloc(m_tripYOrigin, m_tripNoTrips, "m_tripYOrigin", 0.0, AT_);
  mAlloc(m_tripYHeight, m_tripNoTrips, "m_tripYHeight", 1.0, AT_);
  mAlloc(m_tripCutoffZ, m_tripNoTrips, "m_tripCutoffZ", 1.7, AT_);
  mAlloc(m_tripMaxAmpSteady, m_tripNoTrips, "m_tripMaxAmpSteady", 0.0, AT_);
  mAlloc(m_tripMaxAmpFluc, m_tripNoTrips, "m_tripMaxAmpFluc", 0.005, AT_);
  mAlloc(m_tripDeltaTime, m_tripNoTrips, "m_tripDeltaTime", 4.0, AT_);
  mAlloc(m_tripTimeStep, m_tripNoTrips, "m_tripTimeStep", 0, AT_);

  for(MInt i = 0; i < m_tripNoTrips; ++i) {
    /*! \property
      \page propertiesFVSTRCTRD
      \section tripDelta1
      <code>MInt FvStructuredSolver::m_tripDelta1 </code>\n
      default = <code> 1.0 </code>\n \n
      Delta1 boundary layer thickness at position of tripping.\n
      Possible values are:\n
      <ul>
      <li>Float > 0.0</li>
      </ul>
      Keywords: <i>TRIP, BOUNDARYLAYER, STRUCTURED</i>
    */
    m_tripDelta1[i] = Context::getSolverProperty<MFloat>("tripDelta1", m_solverId, AT_, &m_tripDelta1[i], i);

    /*! \property
      \page propertiesFVSTRCTRD
      \section tripXOrigin
      <code>MInt FvStructuredSolver::m_tripXOrigin </code>\n
      default = <code> 30.0 </code>\n \n
      Streamwise center position of the trip forcing.\n
      Possible values are:\n
      <ul>
      <li>Float <> 0.0</li>
      </ul>
      Keywords: <i>TRIP, BOUNDARYLAYER, STRUCTURED</i>
    */
    m_tripXOrigin[i] = Context::getSolverProperty<MFloat>("tripXOrigin", m_solverId, AT_, &m_tripXOrigin[i], i);

    /*! \property
      \page propertiesFVSTRCTRD
      \section tripXLength
      <code>MInt FvStructuredSolver::m_tripXLength </code>\n
      default = <code> 1.0 </code>\n \n
      Streamwise extent of the trip forcing as a factor of delta1.\n
      Possible values are:\n
      <ul>
      <li>Float > 0.0</li>
      </ul>
      Keywords: <i>TRIP, BOUNDARYLAYER, STRUCTURED</i>
    */
    MFloat tripX = 4.0;
    if(Context::propertyExists("tripXLength", m_solverId)) {
      tripX = Context::getSolverProperty<MFloat>("tripXLength", m_solverId, AT_, &tripX);
    }
    m_tripXLength[i] = m_tripDelta1[i] * tripX;

    /*! \property
      \page propertiesFVSTRCTRD
      \section tripYOrigin
      <code>MInt FvStructuredSolver::m_tripYOrigin </code>\n
      default = <code> 0.5 </code>\n \n
      Wall-normal center position of the trip forcing.\n
      Possible values are:\n
      <ul>
      <li>Float > 0.0</li>
      </ul>
      Keywords: <i>TRIP, BOUNDARYLAYER, STRUCTURED</i>
    */
    if(Context::propertyExists("tripYOrigin", m_solverId)) {
      m_tripYOrigin[i] = Context::getSolverProperty<MFloat>("tripYOrigin", m_solverId, AT_, &m_tripYOrigin[i], i);
    }

    /*! \property
      \page propertiesFVSTRCTRD
      \section tripYHeight
      <code>MInt FvStructuredSolver::m_tripYHeight </code>\n
      default = <code> 0.5 </code>\n \n
      Wall-normal extent of the trip forcing.\n
      Possible values are:\n
      <ul>
      <li>Float > 0.0</li>
      </ul>
      Keywords: <i>TRIP, BOUNDARYLAYER, STRUCTURED</i>
    */
    MFloat tripY = 1.0;
    if(Context::propertyExists("tripYHeight", m_solverId)) {
      tripY = Context::getSolverProperty<MFloat>("tripYHeight", m_solverId, AT_, &tripY);
    }
    m_tripYHeight[i] = m_tripDelta1[i] * tripY;

    /*! \property
      \page propertiesFVSTRCTRD
      \section tripCutoffZ
      <code>MInt FvStructuredSolver::m_tripCutoffZ </code>\n
      default = <code> 0.5 </code>\n \n
      Cutoff length in z-direction, all wavenumbers higher than 2*PI*cutoffZ will be set to zero.
      Relative value of delta1\n
      Possible values are:\n
      <ul>
      <li>Float > 0.0</li>
      </ul>
      Keywords: <i>TRIP, BOUNDARYLAYER, STRUCTURED</i>
    */
    MFloat tripZ = 1.7;
    if(Context::propertyExists("tripCutoffZ", m_solverId)) {
      tripZ = Context::getSolverProperty<MFloat>("tripCutoffZ", m_solverId, AT_, &tripZ);
    }
    m_tripCutoffZ[i] = m_tripDelta1[i] * tripZ;


    /*! \property
      \page propertiesFVSTRCTRD
      \section tripMaxAmpSteady
      <code>MInt FvStructuredSolver::m_tripMaxAmpSteady </code>\n
      default = <code> 0.0 </code>\n \n
      Strength of the steady forcing amplitude.\n
      Possible values are:\n
      <ul>
      <li>Float >= 0.0</li>
      </ul>
      Keywords: <i>TRIP, BOUNDARYLAYER, STRUCTURED</i>
    */
    m_tripMaxAmpSteady[i] =
        Context::getSolverProperty<MFloat>("tripMaxAmpSteady", m_solverId, AT_, &m_tripMaxAmpSteady[i], i);

    /*! \property
      \page propertiesFVSTRCTRD
      \section tripMaxAmpFluc
      <code>MInt FvStructuredSolver::m_tripMaxAmpFluc </code>\n
      default = <code> 0.005 </code>\n \n
      Strength of the fluctuating forcing amplitude.\n
      Possible values are:\n
      <ul>
      <li>Float > 0.0</li>
      </ul>
      Keywords: <i>TRIP, BOUNDARYLAYER, STRUCTURED</i>
    */
    m_tripMaxAmpFluc[i] =
        Context::getSolverProperty<MFloat>("tripMaxAmpFluc", m_solverId, AT_, &m_tripMaxAmpFluc[i], i);

    m_tripNoModes = 100;

    /*! \property
      \page propertiesFVSTRCTRD
      \section tripDeltaTime
      <code>MInt FvStructuredSolver::m_tripDeltaTime </code>\n
      default = <code> 2.0 </code>\n \n
      Delta t of the tripping, i.e. the\n
      time step to change the Fourier coeffients.\n
      Possible values are:\n
      <ul>
      <li>Float > 0.0</li>
      </ul>
      Keywords: <i>TRIP, BOUNDARYLAYER, STRUCTURED</i>
    */
    MFloat timeCutoff = 4.0;
    if(Context::propertyExists("tripDeltaTime", m_solverId)) {
      timeCutoff = Context::getSolverProperty<MFloat>("tripDeltaTime", m_solverId, AT_, &timeCutoff);
    }
    m_tripDeltaTime[i] = timeCutoff * m_tripDelta1[i] / PV->UInfinity;
    m_tripTimeStep[i] = (MInt)(m_time / m_tripDeltaTime[i]);
  }

  if(Context::propertyExists("RNGSeed") || Context::propertyExists("seedRNGWithTime"))
    mTerm(1, "Properties RNGSeed or seedRNGWithTime not compatible with SandpaperTrip!");
  m_tripSeed = 70;
  srand(m_tripSeed);
  m_tripDomainWidth = 0.0;
  m_tripNoCells = m_nCells[0];


  m_tripUseRestart = true;
  if(Context::propertyExists("tripUseRestart", m_solverId)) {
    m_tripUseRestart = Context::getSolverProperty<MBool>("tripUseRestart", m_solverId, AT_, &m_tripUseRestart);
  }


  MFloat localDomainWidth = m_cells->coordinates[2][cellIndex(0, 0, m_nCells[0] - 2)];
  MPI_Allreduce(&localDomainWidth, &m_tripDomainWidth, 1, MPI_DOUBLE, MPI_MAX, m_StructuredComm, AT_,
                "localDomainWidth", "m_tripDomainWidth");

  ////////////////////////////////////////////////////
  ///////////// AIRFOIL SPECIFIC PART START //////////
  ////////////////////////////////////////////////////

  m_tripAirfoil = false;
  if(Context::propertyExists("tripAirfoil", m_solverId)) {
    m_tripAirfoil = Context::getSolverProperty<MBool>("tripAirfoil", m_solverId, AT_, &m_tripAirfoil);
  }

  if(m_tripAirfoil) {
    if(!m_movingGrid) {
      const MInt blockId = 0; // use information from inputSolver 0
      m_airfoilNoWallPoints = m_grid->getBlockNoPoints(blockId, 2) + 2 * m_noGhostLayers;

      if(domainId() == 0) {
        cout << "NoWallPoints: " << m_airfoilNoWallPoints << endl;
      }

      MFloatScratchSpace localCoords(nDim, m_airfoilNoWallPoints, AT_, "localCoords");
      MFloatScratchSpace localNormalVec(nDim, m_airfoilNoWallPoints, AT_, "localNormalVec");
      localCoords.fill(-99999.0);
      localNormalVec.fill(-99999.0);
      mAlloc(m_airfoilCoords, nDim * m_airfoilNoWallPoints, "m_airfoilCoords", -99999.0, AT_);
      mAlloc(m_airfoilNormalVec, nDim * m_airfoilNoWallPoints, "m_airfoilNormalVec", -99999.0, AT_);

      if(m_blockId == blockId) {
        if(m_nOffsetCells[1] == 0 && m_nOffsetCells[0] == 0) {
          // only collect variables from points at the wall
          for(MInt i = 0; i < m_nPoints[2]; i++) {
            const MInt offset = m_nOffsetCells[2];
            const MInt localPointId = i; // new postion in Array
            const MInt globalPointId = offset + i;

            const MInt pIJK = pointIndex(i, m_noGhostLayers, m_noGhostLayers);
            const MInt pIJPK = pointIndex(i, m_noGhostLayers + 1, m_noGhostLayers);

            // compute the normal vector
            MFloat normalVec[3] = {F0, F0, F0};
            for(MInt dim = 0; dim < nDim; dim++) {
              normalVec[dim] = m_grid->m_coordinates[dim][pIJPK] - m_grid->m_coordinates[dim][pIJK];
            }

            // normalize normalVec
            const MFloat normalLength = sqrt(POW2(normalVec[0]) + POW2(normalVec[1]) + POW2(normalVec[2]));
            for(MInt dim = 0; dim < nDim; dim++) {
              normalVec[dim] /= normalLength;
            }

            localCoords(globalPointId) = m_grid->m_coordinates[0][localPointId];
            localCoords(m_airfoilNoWallPoints + globalPointId) = m_grid->m_coordinates[1][localPointId];
            localCoords(2 * m_airfoilNoWallPoints + globalPointId) = m_grid->m_coordinates[2][localPointId];

            localNormalVec(globalPointId) = normalVec[0];
            localNormalVec(m_airfoilNoWallPoints + globalPointId) = normalVec[1];
            localNormalVec(2 * m_airfoilNoWallPoints + globalPointId) = normalVec[2];
          }
        }
      }

      MPI_Allreduce(&localCoords[0], &m_airfoilCoords[0], nDim * m_airfoilNoWallPoints, MPI_DOUBLE, MPI_MAX,
                    m_StructuredComm, AT_, "localCoords[0]", "m_airfoilCoords[0]");
      MPI_Allreduce(&localNormalVec[0], &m_airfoilNormalVec[0], nDim * m_airfoilNoWallPoints, MPI_DOUBLE, MPI_MAX,
                    m_StructuredComm, AT_, "localNormalVec[0]", "m_airfoilNormalVec[0]");
    }

    // compute distance to wall
    vector<Point<3>> pts;
    for(MInt i = 0; i < m_airfoilNoWallPoints - 1; ++i) {
      const MInt gIJK = i;
      const MInt gIPJK = i + 1;
      const MFloat xWall = F1B2 * (m_airfoilCoords[0 + gIJK] + m_airfoilCoords[0 + gIPJK]);
      const MFloat yWall =
          F1B2 * (m_airfoilCoords[m_airfoilNoWallPoints + gIJK] + m_airfoilCoords[m_airfoilNoWallPoints + gIPJK]);

      Point<3> a(xWall, yWall, 0.0);
      pts.push_back(a);
    }

    KDtree<3> tree(pts);
    mAlloc(m_airfoilWallDist, m_noCells, "m_airfoilWallDist", -9999.0, AT_);

    for(MInt cellId = 0; cellId < m_noCells; cellId++) {
      const MFloat x = m_cells->coordinates[0][cellId];
      const MFloat y = m_cells->coordinates[1][cellId];

      Point<3> pt(x, y, 0);
      MFloat distance = -1.111111111111;
      (void)tree.nearest(pt, distance);
      m_airfoilWallDist[cellId] = distance;
    }
  }

  ////////////////////////////////////////////////////
  ///////////// AIRFOIL SPECIFIC PART END ////////////
  ////////////////////////////////////////////////////

  m_log << "=================================================" << endl
        << "           SANDPAPER TRIP PROPERTIES             " << endl
        << "=================================================" << endl;
  for(MInt i = 0; i < m_tripNoTrips; ++i) {
    m_log << "######### TRIP NUMBER " << i << " ###########" << endl
          << "tripXOrigin: " << m_tripXOrigin[i] << endl
          << "tripXLength: " << m_tripXLength[i] << endl
          << "tripYOrigin: " << m_tripYOrigin[i] << endl
          << "tripYHeight: " << m_tripYHeight[i] << endl
          << "tripMaxAmpFluc: " << m_tripMaxAmpFluc[i] << endl
          << "tripNoModes: " << m_tripNoModes << endl
          << "tripDeltaTime: " << m_tripDeltaTime[i] << endl
          << "tripTimeStep: " << m_tripTimeStep[i] << endl
          << "tripDomainWidth: " << m_tripDomainWidth << endl
          << "###########################################" << endl;
  }
  m_log << "=================================================" << endl;

  mAlloc(m_tripCoords, m_tripNoCells, "m_tripCoords", F0, AT_);
  mAlloc(m_tripG, m_tripNoTrips * m_tripNoCells, "m_tripG", F0, AT_);
  mAlloc(m_tripH1, m_tripNoTrips * m_tripNoCells, "m_tripH1", F0, AT_);
  mAlloc(m_tripH2, m_tripNoTrips * m_tripNoCells, "m_tripH2", F0, AT_);
  mAlloc(m_tripModesG, m_tripNoTrips * 2 * m_tripNoModes, "m_tripModesG", F0, AT_);
  mAlloc(m_tripModesH1, m_tripNoTrips * 2 * m_tripNoModes, "m_tripModesH1", F0, AT_);
  mAlloc(m_tripModesH2, m_tripNoTrips * 2 * m_tripNoModes, "m_tripModesH2", F0, AT_);


  for(MInt k = 0; k < m_tripNoCells; k++) {
    m_tripCoords[k] = m_cells->coordinates[2][cellIndex(0, 0, k)];
  }

  if(m_restart && m_tripUseRestart) {
    stringstream tripPath;
    tripPath << "/trip";
    ParallelIo::size_type dummyOffset = 0;
    ParallelIo::size_type dataSize = m_tripNoTrips * m_tripNoModes * 2;

    if(domainId() == 0) {
      stringstream restartFileName;
      MString restartFile = Context::getSolverProperty<MString>("restartVariablesFileName", m_solverId, AT_);
      restartFileName << outputDir() << restartFile;
      ParallelIoHdf5 pio(restartFileName.str(), maia::parallel_io::PIO_READ, MPI_COMM_SELF);
      pio.readArray(m_tripModesG, tripPath.str(), "tripModesG", 1, &dummyOffset, &dataSize);
      pio.readArray(m_tripModesH1, tripPath.str(), "tripModesH1", 1, &dummyOffset, &dataSize);
      pio.readArray(m_tripModesH2, tripPath.str(), "tripModesH2", 1, &dummyOffset, &dataSize);
    }

    MPI_Bcast(m_tripModesG, dataSize, MPI_DOUBLE, 0, m_StructuredComm, AT_, "m_tripModesG");
    MPI_Bcast(m_tripModesH1, dataSize, MPI_DOUBLE, 0, m_StructuredComm, AT_, "m_tripModesH1");
    MPI_Bcast(m_tripModesH2, dataSize, MPI_DOUBLE, 0, m_StructuredComm, AT_, "m_tripModesH2");
  } else {
    for(MInt i = 0; i < m_tripNoTrips; ++i) {
      const MInt offsetModes = i * 2 * m_tripNoModes;
      tripFourierCoefficients(&m_tripModesG[offsetModes], m_tripNoModes, m_tripDomainWidth, m_tripCutoffZ[i]);
      tripFourierCoefficients(&m_tripModesH1[offsetModes], m_tripNoModes, m_tripDomainWidth, m_tripCutoffZ[i]);
      tripFourierCoefficients(&m_tripModesH2[offsetModes], m_tripNoModes, m_tripDomainWidth, m_tripCutoffZ[i]);
    }
  }


  for(MInt i = 0; i < m_tripNoTrips; ++i) {
    const MInt offset = i * m_tripNoCells;
    const MInt offsetModes = i * 2 * m_tripNoModes;
    tripForceCoefficients(&m_tripModesG[offsetModes], &m_tripG[offset], m_tripCoords, m_tripNoCells, m_tripNoModes);
    tripForceCoefficients(&m_tripModesH1[offsetModes], &m_tripH1[offset], m_tripCoords, m_tripNoCells, m_tripNoModes);
    tripForceCoefficients(&m_tripModesH2[offsetModes], &m_tripH2[offset], m_tripCoords, m_tripNoCells, m_tripNoModes);
  }
}

void FvStructuredSolver3D::applySandpaperTrip() {
  TRACE();

  for(MInt ii = 0; ii < m_tripNoTrips; ++ii) {
    const MFloat t = m_time + m_timeStep * m_RKalpha[m_RKStep];
    const MInt tripTime = (MInt)(t / m_tripDeltaTime[ii]);
    const MFloat p = t / m_tripDeltaTime[ii] - tripTime;
    const MFloat b = 3 * pow(p, 2) - 2 * pow(p, 3);
    const MInt offset = ii * m_tripNoCells;
    const MInt offsetModes = ii * 2 * m_tripNoModes;

    if(tripTime > m_tripTimeStep[ii]) {
      m_tripTimeStep[ii] = tripTime;

      // copy old values from H2 to H1
      for(MInt k = 0; k < m_tripNoCells; k++) {
        m_tripH1[offset + k] = m_tripH2[offset + k];
      }

      // also copy the old mode coefficients
      for(MInt n = 0; n < 2 * m_tripNoModes; n++) {
        m_tripModesH1[offsetModes + n] = m_tripModesH2[offsetModes + n];
      }

      // compute new fourier coefficients
      tripFourierCoefficients(&m_tripModesH2[offsetModes], m_tripNoModes, m_tripDomainWidth, m_tripCutoffZ[ii]);
      tripForceCoefficients(&m_tripModesH2[offsetModes], &m_tripH2[offset], m_tripCoords, m_tripNoCells, m_tripNoModes);
    }

    for(MInt k = m_noGhostLayers; k < m_nCells[0] - m_noGhostLayers; k++) {
      const MFloat forceStrength =
          (m_tripMaxAmpSteady[ii] * m_tripG[offset + k]
           + m_tripMaxAmpFluc[ii] * ((1.0 - b) * m_tripH1[offset + k] + b * m_tripH2[offset + k]));

      for(MInt j = m_noGhostLayers; j < m_nCells[1] - m_noGhostLayers; j++) {
        for(MInt i = m_noGhostLayers; i < m_nCells[2] - m_noGhostLayers; i++) {
          const MInt cellId = cellIndex(i, j, k);
          const MFloat x = m_cells->coordinates[0][cellId];
          const MFloat y = m_cells->coordinates[1][cellId];

          MFloat force = 0.0;
          if(x > m_tripXOrigin[ii] - 2 * m_tripXLength[ii] && x < m_tripXOrigin[ii] + 2 * m_tripXLength[ii]
             && y > m_tripYOrigin[ii] - 2 * m_tripYHeight[ii] && y < m_tripYOrigin[ii] + 2 * m_tripYHeight[ii]) {
            force = exp(-POW2((x - m_tripXOrigin[ii]) / m_tripXLength[ii])
                        - POW2((y - m_tripYOrigin[ii]) / m_tripYHeight[ii]))
                    * forceStrength;
          }

          m_cells->rightHandSide[CV->RHO_V][cellId] +=
              m_cells->variables[CV->RHO][cellId] * force * m_cells->cellJac[cellId];
          m_cells->rightHandSide[CV->RHO_E][cellId] +=
              m_cells->variables[CV->RHO_V][cellId] * force * m_cells->cellJac[cellId];
        }
      }
    }
  }
}

void FvStructuredSolver3D::applySandpaperTripAirfoil() {
  TRACE();

  for(MInt ii = 0; ii < m_tripNoTrips; ++ii) {
    const MFloat t = m_time + m_timeStep * m_RKalpha[m_RKStep];
    const MInt tripTime = (MInt)(t / m_tripDeltaTime[ii]);
    const MFloat p = t / m_tripDeltaTime[ii] - tripTime;
    const MFloat b = 3 * pow(p, 2) - 2 * pow(p, 3);
    const MInt offset = ii * m_tripNoCells;
    const MInt offsetModes = ii * 2 * m_tripNoModes;

    if(tripTime > m_tripTimeStep[ii]) {
      m_tripTimeStep[ii] = tripTime;

      // copy old values from H2 to H1
      for(MInt k = 0; k < m_tripNoCells; k++) {
        m_tripH1[offset + k] = m_tripH2[offset + k];
      }

      // also copy the old mode coefficients
      for(MInt n = 0; n < 2 * m_tripNoModes; n++) {
        m_tripModesH1[offsetModes + n] = m_tripModesH2[offsetModes + n];
      }

      // compute new fourier coefficients
      tripFourierCoefficients(&m_tripModesH2[offsetModes], m_tripNoModes, m_tripDomainWidth, m_tripCutoffZ[ii]);
      tripForceCoefficients(&m_tripModesH2[offsetModes], &m_tripH2[offset], m_tripCoords, m_tripNoCells, m_tripNoModes);
    }

    for(MInt k = m_noGhostLayers; k < m_nCells[0] - m_noGhostLayers; k++) {
      const MFloat forceStrength =
          (m_tripMaxAmpSteady[ii] * m_tripG[offset + k]
           + m_tripMaxAmpFluc[ii] * ((1.0 - b) * m_tripH1[offset + k] + b * m_tripH2[offset + k]));

      for(MInt j = m_noGhostLayers; j < m_nCells[1] - m_noGhostLayers; j++) {
        for(MInt i = m_noGhostLayers; i < m_nCells[2] - m_noGhostLayers; i++) {
          const MInt cellId = cellIndex(i, j, k);
          const MFloat x = m_cells->coordinates[0][cellId];
          const MFloat y = m_cells->coordinates[1][cellId];
          const MFloat wallDist = m_airfoilWallDist[cellId];

          MFloat force = 0.0;
          if(x > m_tripXOrigin[ii] - 2 * m_tripXLength[ii] && x < m_tripXOrigin[ii] + 2 * m_tripXLength[ii]
             && wallDist < 2 * m_tripYHeight[ii]) {
            MFloat factor = 0.0;
            if(ii == 0 && y > 0.0) {
              factor = exp(-POW2((x - m_tripXOrigin[ii]) / m_tripXLength[ii]) - POW2(wallDist / m_tripYHeight[ii]));
            } else if(ii == 1 && y < 0.0) {
              factor = exp(-POW2((x - m_tripXOrigin[ii]) / m_tripXLength[ii]) - POW2(wallDist / m_tripYHeight[ii]));
            }

            force = forceStrength * factor;
          }

          const MInt offsetPoints = m_nOffsetCells[2];
          const MInt globalPointId = offsetPoints + i;
          const MFloat wallNormalVec[3] = {m_airfoilNormalVec[0 + globalPointId],
                                           m_airfoilNormalVec[m_airfoilNoWallPoints + globalPointId],
                                           m_airfoilNormalVec[2 * m_airfoilNoWallPoints + globalPointId]};

          const MFloat density = m_cells->pvariables[PV->RHO][cellId];
          const MFloat u = m_cells->pvariables[PV->U][cellId];
          const MFloat u_force = wallNormalVec[0] * force;

          m_cells->rightHandSide[CV->RHO_U][cellId] += density * u_force * m_cells->cellJac[cellId];
          m_cells->rightHandSide[CV->RHO_E][cellId] += density * u * u_force * m_cells->cellJac[cellId];

          const MFloat v = m_cells->pvariables[PV->V][cellId];
          const MFloat v_force = wallNormalVec[1] * force;

          m_cells->rightHandSide[CV->RHO_V][cellId] += density * v_force * m_cells->cellJac[cellId];
          m_cells->rightHandSide[CV->RHO_E][cellId] += density * v * v_force * m_cells->cellJac[cellId];
        }
      }
    }
  }
}


void FvStructuredSolver3D::tripForceCoefficients(MFloat* modes, MFloat* forceCoef, MFloat* coords, MInt noCells,
                                                 MInt noModes) {
  MFloat maxLocalValue = -999999.0;
  MFloat minLocalValue = 999999.0;
  MFloat* ak = &modes[0];
  MFloat* phik = &modes[noModes];

  for(MInt k = 0; k < noCells; k++) {
    const MFloat z = coords[k];
    for(MInt n = 0; n < noModes; n++) {
      forceCoef[k] += sin(z * ak[n] + phik[n]);
    }

    maxLocalValue = mMax(maxLocalValue, forceCoef[k]);
    minLocalValue = mMin(minLocalValue, forceCoef[k]);
  }

  // find out min and max to normalize coefficients to interval [-1,1]
  MFloat maxGlobalValue = 0.0;
  MFloat minGlobalValue = 0.0;
  MPI_Allreduce(&maxLocalValue, &maxGlobalValue, 1, MPI_DOUBLE, MPI_MAX, m_StructuredComm, AT_, "maxLocalValue",
                "maxGlobalValue");
  MPI_Allreduce(&minLocalValue, &minGlobalValue, 1, MPI_DOUBLE, MPI_MIN, m_StructuredComm, AT_, "minLocalValue",
                "minGlobalValue");

  // normalize the series
  for(MInt k = 0; k < noCells; k++) {
    forceCoef[k] = 2 * (forceCoef[k] - minGlobalValue) / (maxGlobalValue - minGlobalValue) - 1.0;
  }
}


void FvStructuredSolver3D::tripFourierCoefficients(MFloat* modes, MInt noModes, MFloat maxWaveLength,
                                                   MFloat minWaveLength) {
  const MFloat minWavenumber = 2 * PI / maxWaveLength;
  const MFloat maxWavenumber = 2 * PI / minWaveLength;

  MFloat* ak = &modes[0];
  MFloat* phik = &modes[noModes];
  if(domainId() == 0) {
    for(MInt n = 0; n < noModes; n++) {
      ak[n] = (maxWavenumber - minWavenumber) * (rand() / MFloat(RAND_MAX)) + minWavenumber;
      phik[n] = 2 * PI * rand() / MFloat(RAND_MAX);
    }
  }

  MPI_Bcast(&ak[0], noModes, MPI_DOUBLE, 0, m_StructuredComm, AT_, "ak[0]");
  MPI_Bcast(&phik[0], noModes, MPI_DOUBLE, 0, m_StructuredComm, AT_, "phik[0]");
}


void FvStructuredSolver3D::distributeFluxToCells() { TRACE(); }

void FvStructuredSolver3D::computeCumulativeAverage(MBool forceReset) {
  TRACE();
  if(m_RKStep == 0 && m_zoneType == "LES") {
    if(globalTimeStep == 0 || forceReset) {
      for(MInt cellId = 0; cellId < m_noCells; cellId++) {
        m_cells->fq[FQ->AVG_RHO][cellId] = m_cells->pvariables[PV->RHO][cellId];
        m_cells->fq[FQ->AVG_U][cellId] = m_cells->pvariables[PV->U][cellId];
        m_cells->fq[FQ->AVG_V][cellId] = m_cells->pvariables[PV->V][cellId];
        m_cells->fq[FQ->AVG_W][cellId] = m_cells->pvariables[PV->W][cellId];
        m_cells->fq[FQ->AVG_P][cellId] = m_cells->pvariables[PV->P][cellId];
      }
    } else {
      const MFloat timeFacRho = m_Ma * (1.0 / m_zonalAveragingFactor) * sqrt(PV->TInfinity) * m_physicalTimeStep;
      const MFloat timeFacU = m_Ma * (1.0 / m_zonalAveragingFactor) * sqrt(PV->TInfinity) * m_physicalTimeStep;
      const MFloat timeFacV = m_Ma * (1.0 / m_zonalAveragingFactor) * sqrt(PV->TInfinity) * m_physicalTimeStep;
      const MFloat timeFacW = m_Ma * (1.0 / m_zonalAveragingFactor) * sqrt(PV->TInfinity) * m_physicalTimeStep;
      const MFloat timeFacE = m_Ma * (1.0 / m_zonalAveragingFactor) * sqrt(PV->TInfinity) * m_physicalTimeStep;

      //! do the time average of the flow variables for LES and store them.
      for(MInt cellId = 0; cellId < m_noCells; cellId++) {
        // Exponential averaging
        m_cells->fq[FQ->AVG_RHO][cellId] =
            timeFacRho * m_cells->pvariables[PV->RHO][cellId] + (1.0 - timeFacRho) * m_cells->fq[FQ->AVG_RHO][cellId];
        m_cells->fq[FQ->AVG_U][cellId] =
            timeFacU * m_cells->pvariables[PV->U][cellId] + (1.0 - timeFacU) * m_cells->fq[FQ->AVG_U][cellId];
        m_cells->fq[FQ->AVG_V][cellId] =
            timeFacV * m_cells->pvariables[PV->V][cellId] + (1.0 - timeFacV) * m_cells->fq[FQ->AVG_V][cellId];
        m_cells->fq[FQ->AVG_W][cellId] =
            timeFacW * m_cells->pvariables[PV->W][cellId] + (1.0 - timeFacW) * m_cells->fq[FQ->AVG_W][cellId];
        m_cells->fq[FQ->AVG_P][cellId] =
            timeFacE * m_cells->pvariables[PV->P][cellId] + (1.0 - timeFacE) * m_cells->fq[FQ->AVG_P][cellId];
      }
    }
  }
}

void FvStructuredSolver3D::spanwiseAvgZonal(vector<MFloat*>& variables) {
  if(!variables.empty()) {
    MInt totalNoCellsIJ = (m_grid->getMyBlockNoCells(2) * m_grid->getMyBlockNoCells(1));
    const MInt noVars = variables.size();
    MFloatScratchSpace localSpannwiseVars(totalNoCellsIJ, noVars, AT_, "localSpannwiseVars");
    MFloatScratchSpace globalSpannwiseVars(totalNoCellsIJ, noVars, AT_, "globalSpannwiseVars");

    for(MInt varPos = 0; varPos < noVars; varPos++) {
      for(MInt k = 0; k < m_nActiveCells[0]; k++) {
        for(MInt j = 0; j < m_nActiveCells[1]; j++) {
          for(MInt i = 0; i < m_nActiveCells[2]; i++) {
            MInt localCellId3D =
                i + m_noGhostLayers + ((j + m_noGhostLayers) + (k + m_noGhostLayers) * m_nCells[1]) * m_nCells[2];
            MInt globalCellId2D = (i + m_nOffsetCells[2]) + (j + m_nOffsetCells[1]) * m_grid->getMyBlockNoCells(2);
            localSpannwiseVars(globalCellId2D, varPos) += variables[varPos][localCellId3D];
          }
        }
      }
    }

    MPI_Allreduce(&localSpannwiseVars(0, 0), &globalSpannwiseVars(0, 0), totalNoCellsIJ * noVars, MPI_DOUBLE, MPI_SUM,
                  m_commZonal[m_blockId], AT_, "localSpannwiseVars(0", "0)");

    for(MInt varPos = 0; varPos < noVars; varPos++) {
      for(MInt cellId = 0; cellId < totalNoCellsIJ; cellId++) {
        globalSpannwiseVars(cellId, varPos) /= m_grid->getMyBlockNoCells(0);
      }
    }

    for(MInt varPos = 0; varPos < noVars; varPos++) {
      for(MInt k = 0; k < m_nActiveCells[0]; k++) {
        for(MInt j = 0; j < m_nActiveCells[1]; j++) {
          for(MInt i = 0; i < m_nActiveCells[2]; i++) {
            MInt localCellId3D =
                i + m_noGhostLayers + ((j + m_noGhostLayers) + (k + m_noGhostLayers) * m_nCells[1]) * m_nCells[2];
            MInt globalCellId2D = (i + m_nOffsetCells[2]) + (j + m_nOffsetCells[1]) * m_grid->getMyBlockNoCells(2);
            variables[varPos][localCellId3D] = globalSpannwiseVars(globalCellId2D, varPos);
          }
        }
      }
    }
  }
}


void FvStructuredSolver3D::computeTimeStep() {
  TRACE();

  m_timeStep = 1000.0;
  const MFloat* const RESTRICT dxtx = ALIGNED_F(m_cells->dxt[0]);
  const MFloat* const RESTRICT dxty = ALIGNED_F(m_cells->dxt[1]);
  const MFloat* const RESTRICT dxtz = ALIGNED_F(m_cells->dxt[2]);
  const MFloat* const* const RESTRICT metric = m_cells->cellMetrics;

  for(MInt k = m_noGhostLayers; k < m_nCells[0] - m_noGhostLayers; k++) {
    for(MInt j = m_noGhostLayers; j < m_nCells[1] - m_noGhostLayers; j++) {
      for(MInt i = m_noGhostLayers; i < m_nCells[2] - m_noGhostLayers; i++) {
        const MInt cellId = cellIndex(i, j, k);

        const MFloat Frho = F1 / m_cells->pvariables[PV->RHO][cellId];

        // compute the speed of sound
        const MFloat speedOfSound = sqrt(m_gamma * m_cells->pvariables[PV->P][cellId] * Frho);

        // no need for simplified metrics, since information is already contained
        // in cell metrics
        const MFloat lenXi = sqrt(POW2(metric[0][cellId]) + POW2(metric[1][cellId]) + POW2(metric[2][cellId]));
        const MFloat lenEt = sqrt(POW2(metric[3][cellId]) + POW2(metric[4][cellId]) + POW2(metric[5][cellId]));
        const MFloat lenZe = sqrt(POW2(metric[6][cellId]) + POW2(metric[7][cellId]) + POW2(metric[8][cellId]));

        // contravariant velocities
        MFloat U_c = F0;
        MFloat V_c = F0;
        MFloat W_c = F0;

        for(MInt isd = xsd; isd < nDim; isd++) {
          U_c += m_cells->pvariables[PV->VV[isd]][cellId] * metric[xsd * nDim + isd][cellId];
          V_c += m_cells->pvariables[PV->VV[isd]][cellId] * metric[ysd * nDim + isd][cellId];
          W_c += m_cells->pvariables[PV->VV[isd]][cellId] * metric[zsd * nDim + isd][cellId];
        }

        // subtract grid velocity
        U_c -= dxtx[cellId];
        V_c -= dxty[cellId];
        W_c -= dxtz[cellId];

        U_c = fabs(U_c);
        V_c = fabs(V_c);
        W_c = fabs(W_c);

        // has area information in it due to metric terms
        const MFloat eigenvalue = U_c + V_c + W_c + speedOfSound * (lenXi + lenEt + lenZe);

        // divide volume information (jacobian) through area to get characteristic length for CFL
        const MFloat deltaT = m_cfl * m_cells->cellJac[cellId] / eigenvalue;

        if(m_localTimeStep) {
          m_cells->localTimeStep[cellId] = deltaT;
          m_timeStep = F1;
          m_timeRef = F1;
        } else {
          m_timeStep = mMin(m_timeStep, deltaT);
        }
      }
    }
  }
}


void FvStructuredSolver3D::updateSpongeLayer() {
  TRACE();
  if(m_useSponge) m_structuredBndryCnd->updateSpongeLayer();
}


MBool FvStructuredSolver3D::rungeKuttaStep() {
  TRACE();
  const MInt noVars = CV->noVariables;
  const MUint noCells = m_noCells;
  const MFloat rkAlpha = m_RKalpha[m_RKStep];
  const MFloat rkFactor = rkAlpha * m_timeStep;

  MFloat* const RESTRICT oldVars = ALIGNED_F(m_cells->oldVariables[0]);
  MFloat* const RESTRICT vars = ALIGNED_F(m_cells->variables[0]);
  MFloat* const RESTRICT oldCellJac = ALIGNED_MF(m_cells->oldCellJac);
  const MFloat* const RESTRICT cellJac = ALIGNED_MF(m_cells->cellJac);
  const MFloat* const RESTRICT rhs = ALIGNED_MF(m_cells->rightHandSide[0]);

  // set old variables
  if(m_RKStep == 0) {
    for(MInt v = 0; v < noVars; v++) {
      const MUint offset = v * noCells;
      MFloat* const RESTRICT oldCellVars = ALIGNED_F(oldVars + offset);
      const MFloat* const RESTRICT cellVars = ALIGNED_F(vars + offset);
      maia::parallelFor<true>(0, m_noCells, [=](MInt cellId) { oldCellVars[cellId] = cellVars[cellId]; });
    }
  }

  switch(m_rungeKuttaOrder) {
    case 2: {
      // for moving grids we take the old Jacobian into account
      if(m_localTimeStep) {
        for(MInt v = 0; v < noVars; v++) {
          const MUint cellOffset = v * noCells;
          MFloat* const RESTRICT cellVars = vars + cellOffset;
          const MFloat* const RESTRICT oldCellVars = oldVars + cellOffset;
          const MFloat* const RESTRICT cellRhs = rhs + cellOffset;

          maia::parallelFor<true, nDim>(beginP2(), endM2(), [=](const MInt& i, const MInt& j, const MInt& k) {
            const MInt cellId = i + (j + k * m_nCells[1]) * m_nCells[2];
            const MFloat localRkFactor = rkAlpha * m_cells->localTimeStep[cellId];
            const MFloat factor = localRkFactor / m_cells->cellJac[cellId];
            cellVars[cellId] = oldCellVars[cellId] + factor * cellRhs[cellId];
          });
        }
      } else if(m_movingGrid) {
        for(MInt v = 0; v < noVars; v++) {
          const MUint cellOffset = v * noCells;
          MFloat* const RESTRICT cellVars = vars + cellOffset;
          const MFloat* const RESTRICT oldCellVars = oldVars + cellOffset;
          const MFloat* const RESTRICT cellRhs = rhs + cellOffset;

          maia::parallelFor<true, nDim>(beginP2(), endM2(), [=](const MInt& i, const MInt& j, const MInt& k) {
            const MInt cellId = i + (j + k * m_nCells[1]) * m_nCells[2];
            cellVars[cellId] =
                (oldCellVars[cellId] * oldCellJac[cellId] + rkFactor * cellRhs[cellId]) / cellJac[cellId];
          });
        }
      } else {
        for(MInt v = 0; v < noVars; v++) {
          const MUint cellOffset = v * noCells;
          MFloat* const RESTRICT cellVars = vars + cellOffset;
          const MFloat* const RESTRICT oldCellVars = oldVars + cellOffset;
          const MFloat* const RESTRICT cellRhs = rhs + cellOffset;

          maia::parallelFor<true, nDim>(beginP2(), endM2(), [=](const MInt& i, const MInt& j, const MInt& k) {
            const MInt cellId = i + (j + k * m_nCells[1]) * m_nCells[2];
            const MFloat factor = rkFactor / m_cells->cellJac[cellId];
            cellVars[cellId] = oldCellVars[cellId] + factor * cellRhs[cellId];
          });
        }
      }
      break;
    }
    case 3: {
      for(MInt v = 0; v < noVars; v++) {
        const MUint cellOffset = v * noCells;
        MFloat* const RESTRICT cellVars = vars + cellOffset;
        const MFloat* const RESTRICT oldCellVars = oldVars + cellOffset;
        const MFloat* const RESTRICT cellRhs = rhs + cellOffset;

        maia::parallelFor<true, nDim>(beginP2(), endM2(), [=](const MInt& i, const MInt& j, const MInt& k) {
          const MInt cellId = i + (j + k * m_nCells[1]) * m_nCells[2];
          const MFloat factor = rkFactor / m_cells->cellJac[cellId];
          cellVars[cellId] =
              rkAlpha * cellVars[cellId] + (F1 - rkAlpha) * oldCellVars[cellId] - factor * cellRhs[cellId];
        });
      }
      break;
    }
    default: {
      stringstream errorMessage;
      errorMessage << "Given RungeKutta Order " << m_rungeKuttaOrder << " not implemented! " << endl;
      mTerm(1, AT_, errorMessage.str());
    }
  }

  ++m_RKStep;

  if(m_RKStep == m_noRKSteps) {
    m_physicalTime += m_timeStep * m_timeRef;
    m_time += m_timeStep;

    m_RKStep = 0;

    return true;
  } else {
    return false;
  }
}

void FvStructuredSolver3D::addDisturbance() {
  TRACE();
  cout << "enterin addDisturbance " << endl;
  cout << " m_time = " << m_time << endl;
  cout << "m_timeStep = " << m_timeStep << endl;
  cout << " global t= " << globalTimeStep << endl;
  if(m_time >= 0) {
    MFloat pi2 = 8.0 * atan(1);
    MFloat period = 0.1;
    MFloat amp1 = 0.00001 * sin((pi2 / period) * m_time);

    // MFloat T=(m_cells->variables[CV->RHO_E][364])*(m_gamma-F1)/287.1500;
    // MFloat fluc_p = rhsdist*287.15000*T;
    // MFloat fluc_rhoE= (F1/(m_gamma-1))*fluc_p;
    // rhsdist=fluc_rhoE;
    MInt cellId = 0;
    cout << "sin " << amp1 << endl;
    MFloat centerloc = 3.1415;
    // go through all cells and adopt disturbance smoothly!!
    // m_cells->rightHandSide[CV->RHO_E][2191941]+=rhsdist;
    for(MInt k = m_noGhostLayers; k < m_nCells[0] - m_noGhostLayers; k++) {
      for(MInt j = m_noGhostLayers; j < m_nCells[1] - m_noGhostLayers; j++) {
        for(MInt i = m_noGhostLayers; i < m_nCells[2] - m_noGhostLayers; i++) {
          cellId = i + (j + k * m_nCells[1]) * m_nCells[2];
          MFloat x = m_cells->coordinates[0][cellId];
          MFloat y = m_cells->coordinates[1][cellId];
          MFloat z = m_cells->coordinates[2][cellId];
          MFloat r = 0.025;
          MFloat a1 = POW2(x - centerloc);
          MFloat a2 = POW2(y - centerloc);
          MFloat a3 = POW2(z - centerloc);
          MFloat disturb = amp1 * exp(-(a1 + a2 + a3) / POW2(r) / 2.0);
          m_cells->rightHandSide[CV->RHO_E][cellId] += disturb;
          /*    MFloat x=m_cells->coordinates[0][cellId]-centreloc;
                MFloat y=m_cells->coordinates[1][cellId]-centreloc;
                MFloat z=m_cells->coordinates[2][cellId]-centreloc;
                m_cells->rightHandSide[CV->RHO_E][cellId]+=(F1+tanh(x/0.009))*(F1-tanh(x/0.009))*(F1+tanh(y/0.009))*(F1-tanh(y/0.009))*(F1+tanh(z/0.009))*(F1-tanh(z/0.009))*rhsdist;*/
        }
      }
    }
  }

  cout << "leaving addDisturbance " << endl;
}

void FvStructuredSolver3D::assignBndryCells() {
  TRACE();
  m_structuredBndryCnd->assignBndryCnds();
}

void FvStructuredSolver3D::initBndryCnds() {
  TRACE();
  m_structuredBndryCnd->correctBndryCndIndices();
}

/**
 * \brief Extrapolates and exchanges ghost point coordinates
 * \author Pascal Meysonnat
 */
void FvStructuredSolver3D::addGhostPointCoordinateValues() {
  TRACE();
  // > for debugging only: save cellId/pointId and domainId for all cells and points
  if(m_debugOutput) {
    for(MInt k = 0; k < (m_nCells[0]); k++) {
      for(MInt j = 0; j < (m_nCells[1]); j++) {
        for(MInt i = 0; i < (m_nCells[2]); i++) {
          MInt cellId = cellIndex(i, j, k);
          m_cells->fq[FQ->CELLID][cellId] = cellId;
          m_cells->fq[FQ->BLOCKID][cellId] = domainId();
        }
      }
    }
  }
  //< end debugging only

  m_grid->extrapolateGhostPointCoordinates();

  m_grid->exchangePoints(m_sndComm, m_rcvComm, PARTITION_NORMAL);
  m_grid->exchangePoints(m_sndComm, m_rcvComm, PERIODIC_BC);

  extrapolateGhostPointCoordinatesBC();

  m_grid->computeCellCenterCoordinates();

  // MUST be done after cell center computation!!!
  m_grid->exchangePoints(m_sndComm, m_rcvComm, SINGULAR);
  m_grid->exchangePoints(m_sndComm, m_rcvComm, PERIODIC_BC_SINGULAR);

  if(m_hasSingularity > 0) {
    computeReconstructionConstantsSVD();
  }

  if(m_savePartitionOutput) {
    m_grid->writePartitionedGrid();
  }
}

void FvStructuredSolver3D::extrapolateGhostPointCoordinatesBC() {
  for(MInt bcId = 0; bcId < (MInt)m_structuredBndryCnd->m_physicalBCMap.size(); ++bcId) {
    // all the periodic BCs are NOT included.
    // also skip the channel bc
    if(m_structuredBndryCnd->m_physicalBCMap[bcId]->BC == 2401
       || m_structuredBndryCnd->m_physicalBCMap[bcId]->BC == 2402
       || (m_structuredBndryCnd->m_physicalBCMap[bcId]->BC >= 6000
           && m_structuredBndryCnd->m_physicalBCMap[bcId]->BC < 6010)) {
      continue;
    }

    MInt* start = m_structuredBndryCnd->m_physicalBCMap[bcId]->start1;
    MInt* end = m_structuredBndryCnd->m_physicalBCMap[bcId]->end1;
    MInt index = m_structuredBndryCnd->m_physicalBCMap[bcId]->face / 2;
    MInt step = m_structuredBndryCnd->m_physicalBCMap[bcId]->face % 2;
    MInt pos[3], fix[3], mirror[3], ijk[3], extendijk[3];
    MInt pointId, FixPointId, MirrorPointId;

    extendijk[0] = 1;
    extendijk[1] = 1;
    extendijk[2] = 1;
    extendijk[index] = 0;

    for(ijk[2] = start[2]; ijk[2] < end[2] + extendijk[2]; ++ijk[2]) {
      for(ijk[1] = start[1]; ijk[1] < end[1] + extendijk[1]; ++ijk[1]) {
        for(ijk[0] = start[0]; ijk[0] < end[0] + extendijk[0]; ++ijk[0]) {
          for(MInt m = 0; m < 3; ++m) {
            if(index == m) {
              if(step == 1) {
                pos[m] = ijk[m] + 1;
                fix[m] = start[m];
                mirror[m] = 2 * fix[m] - pos[m];
              } else {
                pos[m] = ijk[m];
                fix[m] = end[m];
                mirror[m] = 2 * fix[m] - pos[m];
              }
            } else {
              pos[m] = ijk[m];
              fix[m] = ijk[m];
              mirror[m] = ijk[m];
            }
          } // m

          pointId = pointIndex(pos[0], pos[1], pos[2]);
          FixPointId = pointIndex(fix[0], fix[1], fix[2]);
          MirrorPointId = pointIndex(mirror[0], mirror[1], mirror[2]);

          for(MInt dim = 0; dim < nDim; dim++) {
            m_grid->m_coordinates[dim][pointId] =
                (2 * m_grid->m_coordinates[dim][FixPointId] - m_grid->m_coordinates[dim][MirrorPointId]);
          }
        } // ijk
      }
    }
  } // bcid
}


//============================================================================================================
//====================================COMMUNICATIONS==========================================================


void FvStructuredSolver3D::gather(const MBool periodicExchange,
                                  std::vector<std::unique_ptr<StructuredComm<nDim>>>& sndComm) {
  for(auto& snd : sndComm) {
    if(isPeriodicComm(snd) && !periodicExchange) continue;
    if(periodicExchange && skipPeriodicDirection(snd)) continue;

    std::array<MInt, nDim> begin{snd->startInfoCells[0], snd->startInfoCells[1], snd->startInfoCells[2]};
    std::array<MInt, nDim> end{snd->endInfoCells[0], snd->endInfoCells[1], snd->endInfoCells[2]};
    std::array<MInt, nDim> size{snd->endInfoCells[0] - snd->startInfoCells[0],
                                snd->endInfoCells[1] - snd->startInfoCells[1],
                                snd->endInfoCells[2] - snd->startInfoCells[2]};
    const MInt totalSize = size[0] * size[1] * size[2];

    // Aliasing the unique_pointer in snd to a raw point is needed for PSTL on NVHPC
    auto* cellBuffer = snd->cellBuffer.get();
    auto* variables = &(snd->variables[0]);
    for(MInt var = 0; var < snd->noVars; var++) {
      maia::parallelFor<true, nDim>(begin, end, [=](const MInt& i, const MInt& j, const MInt& k) {
        const MInt cellId = i + (j + k * m_nCells[1]) * m_nCells[2];
        const MInt bufferId = totalSize * var + (i - begin[0]) + ((j - begin[1]) + (k - begin[2]) * size[1]) * size[0];
        cellBuffer[bufferId] = variables[var][cellId];
      });
    }
  }
}


void FvStructuredSolver3D::scatter(const MBool periodicExchange,
                                   std::vector<std::unique_ptr<StructuredComm<nDim>>>& rcvComm) {
  // the ordering of the grid points can be different from
  // sending instance ==> reorder it and copy it to the
  // right place

  for(auto& rcv : rcvComm) {
    if(isPeriodicComm(rcv) && !periodicExchange) continue;
    if(periodicExchange && skipPeriodicDirection(rcv)) continue;

    std::array<MInt, nDim> begin{rcv->startInfoCells[0], rcv->startInfoCells[1], rcv->startInfoCells[2]};
    std::array<MInt, nDim> end{rcv->endInfoCells[0], rcv->endInfoCells[1], rcv->endInfoCells[2]};
    std::array<MInt, nDim> size{rcv->endInfoCells[0] - rcv->startInfoCells[0],
                                rcv->endInfoCells[1] - rcv->startInfoCells[1],
                                rcv->endInfoCells[2] - rcv->startInfoCells[2]};
    const MInt totalSize = size[0] * size[1] * size[2];

    std::array<MInt, nDim> stepBuffer{0};
    std::array<MInt, nDim> startBuffer{0};
    std::array<MInt, nDim> endBuffer{0};
    std::array<MInt, nDim> sizeBuffer{0};

    for(MInt j = 0; j < nDim; j++) {
      stepBuffer[rcv->orderInfo[j]] = rcv->stepInfo[j];
    }

    for(MInt j = 0; j < nDim; j++) {
      endBuffer[j] = size[j] - 1;
      sizeBuffer[rcv->orderInfo[j]] = size[j];
      if(stepBuffer[j] < 0) {
        std::swap(startBuffer[j], endBuffer[j]);
      }
    }

    // Aliasing the unique_pointer in rcv to a raw point is needed for PSTL on NVHPC
    auto* orderInfo = rcv->orderInfo.data();
    auto* startInfoCells = rcv->startInfoCells.data();
    auto* cellBuffer = rcv->cellBuffer.get();
    auto* variables = &(rcv->variables[0]);
    for(MInt var = 0; var < rcv->noVars; var++) {
      maia::parallelFor<true, nDim>(begin, end, [=](const MInt& i, const MInt& j, const MInt& k) {
        std::array<MInt, nDim> start{};
        start[orderInfo[0]] = startBuffer[0] + (i - startInfoCells[0]) * stepBuffer[0];
        start[orderInfo[1]] = startBuffer[1] + (j - startInfoCells[1]) * stepBuffer[1];
        start[orderInfo[2]] = startBuffer[2] + (k - startInfoCells[2]) * stepBuffer[2];

        const MInt bufferId = var * totalSize + start[0] + (start[1] + start[2] * sizeBuffer[1]) * sizeBuffer[0];
        const MInt cellId = i + (j + k * m_nCells[1]) * m_nCells[2];
        variables[var][cellId] = cellBuffer[bufferId];
      });
    }
  }
}

void FvStructuredSolver3D::waveExchange() {
  std::vector<MPI_Request> sndRequests;
  std::vector<MPI_Request> rcvRequests;
  std::vector<MPI_Status> sndStatus;
  std::vector<MPI_Status> rcvStatus;
  sndRequests.reserve(m_waveSndComm.size());
  rcvRequests.reserve(m_waveRcvComm.size());

  waveGather();
  waveSend(sndRequests);
  waveReceive(rcvRequests);

  sndStatus.resize(sndRequests.size());
  MPI_Waitall(sndRequests.size(), &sndRequests[0], &sndStatus[0], AT_);
  rcvStatus.resize(rcvRequests.size());
  MPI_Waitall(rcvRequests.size(), &rcvRequests[0], &rcvStatus[0], AT_);

  waveScatter();
}


void FvStructuredSolver3D::waveGather() {
  for(auto& snd : m_waveSndComm) {
    MInt pos = 0;

    for(MInt var = 0; var < PV->noVariables; var++) {
      for(MInt k = snd->startInfoCells[2]; k < snd->endInfoCells[2]; k++) {
        for(MInt j = snd->startInfoCells[1]; j < snd->endInfoCells[1]; j++) {
          for(MInt i = snd->startInfoCells[0]; i < snd->endInfoCells[0]; i++) {
            const MInt cellId = i + (j + k * m_nCells[1]) * m_nCells[2];
            snd->cellBuffer[pos] = m_cells->pvariables[var][cellId];
            pos++;
          }
        }
      }
    }

    if(m_averageVorticity) {
      for(MInt var = 0; var < nDim; var++) {
        for(MInt k = snd->startInfoCells[2]; k < snd->endInfoCells[2]; k++) {
          for(MInt j = snd->startInfoCells[1]; j < snd->endInfoCells[1]; j++) {
            for(MInt i = snd->startInfoCells[0]; i < snd->endInfoCells[0]; i++) {
              const MInt cellId = i + (j + k * m_nCells[1]) * m_nCells[2];
              snd->cellBuffer[pos] = m_cells->fq[FQ->VORTICITY[var]][cellId];
              pos++;
            }
          }
        }
      }
    }
  }
}

void FvStructuredSolver3D::waveSend(std::vector<MPI_Request>& sndRequests) {
  for(auto& snd : m_waveSndComm) {
    MPI_Request request{};
    const MInt tag = domainId() + (snd->tagHelper) * noDomains();
    const MInt err = MPI_Isend((void*)&snd->cellBuffer[0], snd->cellBufferSize, MPI_DOUBLE, snd->nghbrId, tag,
                               m_StructuredComm, &request, AT_, "snd->cellBuffer");
    sndRequests.push_back(request);
    if(err) cout << "rank " << domainId() << " sending throws error " << endl;
  }
}

void FvStructuredSolver3D::waveReceive(std::vector<MPI_Request>& rcvRequests) {
  for(auto& rcv : m_waveRcvComm) {
    MPI_Request request{};
    const MInt tag = rcv->nghbrId + (rcv->tagHelper) * noDomains();
    const MInt err = MPI_Irecv((void*)&rcv->cellBuffer[0], rcv->cellBufferSize, MPI_DOUBLE, rcv->nghbrId, tag,
                               m_StructuredComm, &request, AT_, "rcv->cellBuffer");
    rcvRequests.push_back(request);
    if(err) cout << "rank " << domainId() << " sending throws error " << endl;
  }
}

void FvStructuredSolver3D::waveScatter() {
  for(auto& rcv : m_waveRcvComm) {
    MInt pos = 0;

    for(MInt var = 0; var < PV->noVariables; var++) {
      for(MInt k = rcv->startInfoCells[2]; k < rcv->endInfoCells[2]; k++) {
        for(MInt j = rcv->startInfoCells[1]; j < rcv->endInfoCells[1]; j++) {
          for(MInt i = rcv->startInfoCells[0]; i < rcv->endInfoCells[0]; i++) {
            const MInt cellId = i + (j + k * m_nCells[1]) * m_nCells[2];
            m_tempWaveSample[var][cellId] = rcv->cellBuffer[pos];
            pos++;
          }
        }
      }
    }

    if(m_averageVorticity) {
      for(MInt var = 0; var < (2 * nDim - 3); var++) {
        for(MInt k = rcv->startInfoCells[2]; k < rcv->endInfoCells[2]; k++) {
          for(MInt j = rcv->startInfoCells[1]; j < rcv->endInfoCells[1]; j++) {
            for(MInt i = rcv->startInfoCells[0]; i < rcv->endInfoCells[0]; i++) {
              const MInt cellId = i + (j + k * m_nCells[1]) * m_nCells[2];
              m_tempWaveSample[PV->noVariables + var][cellId] = rcv->cellBuffer[pos];
              pos++;
            }
          }
        }
      }
    }
  }
}

void FvStructuredSolver3D::spanwiseWaveReorder() {
  RECORD_TIMER_START(m_timers[Timers::WaveSpanwiseReordering]);

  MInt allCellsK = m_grid->getBlockNoCells(0, 0);
  MInt waveZeroPos = ((MInt)round((globalTimeStep - m_movingGridStepOffset) / m_waveNoStepsPerCell)) % allCellsK;

  if(m_waveSpeed < 0.0) {
    waveZeroPos =
        allCellsK - ((MInt)round((globalTimeStep - m_movingGridStepOffset) / m_waveNoStepsPerCell)) % allCellsK;
  }

  m_windowInfo->createWaveWindowMapping(waveZeroPos);

  MInt noVars = PV->noVariables;
  if(m_averageVorticity) {
    noVars += (2 * nDim - 3);
    computeVorticity();
  }

  m_windowInfo->createWaveCommunicationExchangeFlags(m_waveSndComm, m_waveRcvComm, noVars);
  waveExchange();

  RECORD_TIMER_STOP(m_timers[Timers::WaveSpanwiseReordering]);
}


void FvStructuredSolver3D::gcFillGhostCells(vector<MFloat*>& variables) {
  gcExtrapolate(variables);

  std::vector<std::unique_ptr<StructuredComm<nDim>>> gcSndComm;
  std::vector<std::unique_ptr<StructuredComm<nDim>>> gcRcvComm;

  MFloat* const* const varPtr = &variables[0];

  m_windowInfo->createCommunicationExchangeFlags(gcSndComm, gcRcvComm, (MInt)variables.size(), varPtr);
  exchange(gcSndComm, gcRcvComm);
}

void FvStructuredSolver3D::gcExtrapolate(vector<MFloat*>& variables) {
  TRACE();
  // i-direction
  MInt cellId, cellIdAdj1, cellIdAdj2;
  const MInt noVars = variables.size();
  for(MInt k = 0; k < m_nCells[0]; k++) {
    for(MInt j = 0; j < m_nCells[1]; j++) {
      for(MInt i = 0; i < m_noGhostLayers; i++) {
        cellId = cellIndex(m_noGhostLayers - 1 - i, j, k);
        cellIdAdj1 = cellIndex(m_noGhostLayers - i, j, k);
        cellIdAdj2 = cellIndex(m_noGhostLayers + 1 - i, j, k);

        for(MInt varPos = 0; varPos < noVars; varPos++) {
          variables[varPos][cellId] = (F2 * variables[varPos][cellIdAdj1] - variables[varPos][cellIdAdj2]);
        }

        cellId = cellIndex(m_nCells[2] - m_noGhostLayers + i, j, k);
        cellIdAdj1 = cellIndex(m_nCells[2] - m_noGhostLayers - 1 + i, j, k);
        cellIdAdj2 = cellIndex(m_nCells[2] - m_noGhostLayers - 2 + i, j, k);

        for(MInt varPos = 0; varPos < noVars; varPos++) {
          variables[varPos][cellId] = (F2 * variables[varPos][cellIdAdj1] - variables[varPos][cellIdAdj2]);
        }
      }
    }
  }

  // j-direction
  for(MInt k = 0; k < m_nCells[0]; k++) {
    for(MInt i = 0; i < m_nCells[2]; i++) {
      for(MInt j = 0; j < m_noGhostLayers; j++) {
        cellId = cellIndex(i, m_noGhostLayers - 1 - j, k);
        cellIdAdj1 = cellIndex(i, m_noGhostLayers - j, k);
        cellIdAdj2 = cellIndex(i, m_noGhostLayers + 1 - j, k);

        for(MInt varPos = 0; varPos < noVars; varPos++) {
          variables[varPos][cellId] = (F2 * variables[varPos][cellIdAdj1] - variables[varPos][cellIdAdj2]);
        }

        cellId = cellIndex(i, m_nCells[1] - m_noGhostLayers + j, k); // pointId in Array
        cellIdAdj1 = cellIndex(i, m_nCells[1] - m_noGhostLayers - 1 + j, k);
        cellIdAdj2 = cellIndex(i, m_nCells[1] - m_noGhostLayers - 2 + j, k);

        for(MInt varPos = 0; varPos < noVars; varPos++) {
          variables[varPos][cellId] = (F2 * variables[varPos][cellIdAdj1] - variables[varPos][cellIdAdj2]);
        }
      }
    }
  }

  // k-direction
  for(MInt j = 0; j < m_nCells[1]; j++) {
    for(MInt i = 0; i < m_nCells[2]; i++) {
      for(MInt k = 0; k < m_noGhostLayers; k++) {
        cellId = cellIndex(i, j, m_noGhostLayers - 1 - k);
        cellIdAdj1 = cellIndex(i, j, m_noGhostLayers - k);
        cellIdAdj2 = cellIndex(i, j, m_noGhostLayers + 1 - k);

        for(MInt varPos = 0; varPos < noVars; varPos++) {
          variables[varPos][cellId] = (F2 * variables[varPos][cellIdAdj1] - variables[varPos][cellIdAdj2]);
        }

        cellId = cellIndex(i, j, m_nCells[0] - m_noGhostLayers + k);
        cellIdAdj1 = cellIndex(i, j, m_nCells[0] - m_noGhostLayers - 1 + k);
        cellIdAdj2 = cellIndex(i, j, m_nCells[0] - m_noGhostLayers - 2 + k);

        for(MInt varPos = 0; varPos < noVars; varPos++) {
          variables[varPos][cellId] = (F2 * variables[varPos][cellIdAdj1] - variables[varPos][cellIdAdj2]);
        }
      }
    }
  }
}


//====================================COMMUNICATIONS==========================================================
//============================================================================================================
//============================================================================================================


inline MFloat FvStructuredSolver3D::dist(MFloat* a, MFloat* b) {
  MFloat dist1 = F0;
  for(MInt dim = 0; dim < 3; dim++) {
    dist1 += POW2(a[dim * m_noCells] - b[dim * m_noCells]);
  }
  return sqrt(dist1);
}

inline MInt FvStructuredSolver3D::cellIndex(const MInt i, const MInt j, const MInt k) {
  return i + (j + k * m_nCells[1]) * m_nCells[2];
}

inline MInt FvStructuredSolver3D::getCellIdfromCell(const MInt origin, const MInt incI, const MInt incJ,
                                                    const MInt incK) {
  return origin + incI + incJ * m_nCells[2] + incK * m_nCells[2] * m_nCells[1];
}

inline MInt FvStructuredSolver3D::pointIndex(const MInt i, const MInt j, const MInt k) {
  return i + (j + k * m_nPoints[1]) * m_nPoints[2];
}


void FvStructuredSolver3D::viscousFlux() { (this->*viscFluxMethod)(); }

void FvStructuredSolver3D::viscousFluxRANS() { m_ransSolver->viscousFluxRANS(); }

/**
 * \brief Viscous flux computation
 */
template <MBool twoEqRans>
void FvStructuredSolver3D::viscousFluxLES() {
  TRACE();
  const MFloat rPrLam = F1 / m_Pr;
  constexpr MFloat rPrTurb = F1 / 0.9;
  const MFloat rRe = F1 / m_Re0;
  const MFloat gammaMinusOne = m_gamma - 1.0;
  const MFloat FgammaMinusOne = F1 / gammaMinusOne;

  const MFloat* const RESTRICT u = ALIGNED_F(&m_cells->pvariables[PV->U][0]);
  const MFloat* const RESTRICT v = ALIGNED_F(&m_cells->pvariables[PV->V][0]);
  const MFloat* const RESTRICT w = ALIGNED_F(&m_cells->pvariables[PV->W][0]);
  const MFloat* const RESTRICT p = ALIGNED_F(&m_cells->pvariables[PV->P][0]);
  const MFloat* const RESTRICT rho = ALIGNED_F(&m_cells->pvariables[PV->RHO][0]);
  MFloat* const RESTRICT T = ALIGNED_F(&m_cells->temperature[0]);
  MFloat* const RESTRICT muLam = ALIGNED_F(&m_cells->fq[FQ->MU_L][0]);
  MFloat* const RESTRICT muTurb = ALIGNED_F(&m_cells->fq[FQ->MU_T][0]); // this is zero for LES

  MFloat* const* const RESTRICT eflux = ALIGNED_F(m_cells->eFlux);
  MFloat* const* const RESTRICT fflux = ALIGNED_F(m_cells->fFlux);
  MFloat* const* const RESTRICT gflux = ALIGNED_F(m_cells->gFlux);
  MFloat* const* const RESTRICT vflux = ALIGNED_F(m_cells->viscousFlux);

  // only relevant for 2-eq Rans model (Waiting for if constexpr)
  const MFloat* const RESTRICT TKE = (twoEqRans) ? ALIGNED_F(m_cells->pvariables[PV->RANS_VAR[0]]) : nullptr;

  maia::parallelFor<true, nDim>(beginP1(), endM1(), [=](const MInt& i, const MInt& j, const MInt& k) {
    const MInt I = cellIndex(i, j, k);
    T[I] = m_gamma * p[I] / rho[I];
    muLam[I] = SUTHERLANDLAW(T[I]);
  });

  maia::parallelFor<true, nDim>(beginP1(), endM1(), [=](const MInt& i, const MInt& j, const MInt& k) {
    // get the adjacent cells;
    const MInt IJK = cellIndex(i, j, k);
    const MInt IPJK = cellIndex((i + 1), j, k);
    const MInt IPJPK = cellIndex((i + 1), (j + 1), k);
    const MInt IJPK = cellIndex(i, (j + 1), k);
    const MInt IJKP = cellIndex(i, j, (k + 1));
    const MInt IPJKP = cellIndex((i + 1), j, (k + 1));
    const MInt IPJPKP = cellIndex((i + 1), (j + 1), (k + 1));
    const MInt IJPKP = cellIndex(i, (j + 1), (k + 1));

    const MFloat cornerMetrics[9] = {
        m_cells->cornerMetrics[0][IJK], m_cells->cornerMetrics[1][IJK], m_cells->cornerMetrics[2][IJK],
        m_cells->cornerMetrics[3][IJK], m_cells->cornerMetrics[4][IJK], m_cells->cornerMetrics[5][IJK],
        m_cells->cornerMetrics[6][IJK], m_cells->cornerMetrics[7][IJK], m_cells->cornerMetrics[8][IJK]};


    const MFloat dudxi = F1B4 * (u[IPJPKP] + u[IPJPK] + u[IPJKP] + u[IPJK] - u[IJPKP] - u[IJPK] - u[IJKP] - u[IJK]);
    const MFloat dudet = F1B4 * (u[IPJPKP] + u[IJPKP] + u[IPJPK] + u[IJPK] - u[IPJKP] - u[IJKP] - u[IPJK] - u[IJK]);
    const MFloat dudze = F1B4 * (u[IPJPKP] + u[IJPKP] + u[IPJKP] + u[IJKP] - u[IPJPK] - u[IJPK] - u[IPJK] - u[IJK]);

    const MFloat dvdxi = F1B4 * (v[IPJPKP] + v[IPJPK] + v[IPJKP] + v[IPJK] - v[IJPKP] - v[IJPK] - v[IJKP] - v[IJK]);
    const MFloat dvdet = F1B4 * (v[IPJPKP] + v[IJPKP] + v[IPJPK] + v[IJPK] - v[IPJKP] - v[IJKP] - v[IPJK] - v[IJK]);
    const MFloat dvdze = F1B4 * (v[IPJPKP] + v[IJPKP] + v[IPJKP] + v[IJKP] - v[IPJPK] - v[IJPK] - v[IPJK] - v[IJK]);

    const MFloat dwdxi = F1B4 * (w[IPJPKP] + w[IPJPK] + w[IPJKP] + w[IPJK] - w[IJPKP] - w[IJPK] - w[IJKP] - w[IJK]);
    const MFloat dwdet = F1B4 * (w[IPJPKP] + w[IJPKP] + w[IPJPK] + w[IJPK] - w[IPJKP] - w[IJKP] - w[IPJK] - w[IJK]);
    const MFloat dwdze = F1B4 * (w[IPJPKP] + w[IJPKP] + w[IPJKP] + w[IJKP] - w[IPJPK] - w[IJPK] - w[IPJK] - w[IJK]);

    const MFloat dTdxi = F1B4 * (T[IPJPKP] + T[IPJPK] + T[IPJKP] + T[IPJK] - T[IJPKP] - T[IJPK] - T[IJKP] - T[IJK]);
    const MFloat dTdet = F1B4 * (T[IPJPKP] + T[IJPKP] + T[IPJPK] + T[IJPK] - T[IPJKP] - T[IJKP] - T[IPJK] - T[IJK]);
    const MFloat dTdze = F1B4 * (T[IPJPKP] + T[IJPKP] + T[IPJKP] + T[IJKP] - T[IPJPK] - T[IJPK] - T[IPJK] - T[IJK]);

    const MFloat uAvg = F1B8 * (u[IPJPKP] + u[IJPKP] + u[IJPK] + u[IPJPK] + u[IPJKP] + u[IJKP] + u[IJK] + u[IPJK]);
    const MFloat vAvg = F1B8 * (v[IPJPKP] + v[IJPKP] + v[IJPK] + v[IPJPK] + v[IPJKP] + v[IJKP] + v[IJK] + v[IPJK]);
    const MFloat wAvg = F1B8 * (w[IPJPKP] + w[IJPKP] + w[IJPK] + w[IPJPK] + w[IPJKP] + w[IJKP] + w[IJK] + w[IPJK]);

    const MFloat muLamAvg = F1B8
                            * (muLam[IPJPKP] + muLam[IJPKP] + muLam[IJPK] + muLam[IPJPK] + muLam[IPJKP] + muLam[IJKP]
                               + muLam[IJK] + muLam[IPJK]);

    // turbulent viscosity is set to zero for LES
    // and is only non-zero for RANS
    const MFloat muTurbAvg = F1B8
                             * (muTurb[IPJPKP] + muTurb[IJPKP] + muTurb[IJPK] + muTurb[IPJPK] + muTurb[IPJKP]
                                + muTurb[IJKP] + muTurb[IJK] + muTurb[IPJK]);

    MFloat TKEcorner = F0;
    if(twoEqRans)
      TKEcorner =
          -2 / 3 * F1B8
          * (rho[IPJPKP] * TKE[IPJPKP] + rho[IJPKP] * TKE[IJPKP] + rho[IJPK] * TKE[IJPK] + rho[IPJPK] * TKE[IPJPK]
             + rho[IPJKP] * TKE[IPJKP] + rho[IJKP] * TKE[IJKP] + rho[IJK] * TKE[IJK] + rho[IPJK] * TKE[IPJK]);

    // compute tau1 = 2 du/dx - 2/3 ( du/dx + dv/dy + dw/dz )

    // tau_xx = 4/3*( du/dxi * dxi/dx + du/deta * deta/dx + du/dzeta * dzeta/dx )
    //            - 2/3*( dv/dxi * dxi/dy + dv/deta * deta/dy + dv/dzeta * dzeta/dy)
    //            - 2/3*( dw/dxi * dxi/dz + dw/deta * deta/dz + dw/dzeta * dzeta/dz )
    MFloat tau1 = F4B3
                      * (dudxi * cornerMetrics[xsd * 3 + xsd] + dudet * cornerMetrics[ysd * 3 + xsd]
                         + dudze * cornerMetrics[zsd * 3 + xsd])
                  -

                  F2B3
                      * (dvdxi * cornerMetrics[xsd * 3 + ysd] + dvdet * cornerMetrics[ysd * 3 + ysd]
                         + dvdze * cornerMetrics[zsd * 3 + ysd])
                  -

                  F2B3
                      * (dwdxi * cornerMetrics[xsd * 3 + zsd] + dwdet * cornerMetrics[ysd * 3 + zsd]
                         + dwdze * cornerMetrics[zsd * 3 + zsd]);

    // compute tau2 = du/dy + dv/dx

    // tau_xy = du/dxi * dxi/dy + du/deta * deta/dy + du/dzeta * dzeta/dy
    //        + dv/dxi * dxi/dx + dv/deta * deta/dx + dv/dzeta * dzeta/dx
    MFloat tau2 = dudxi * cornerMetrics[xsd * 3 + ysd] + dudet * cornerMetrics[ysd * 3 + ysd]
                  + dudze * cornerMetrics[zsd * 3 + ysd] +

                  dvdxi * cornerMetrics[xsd * 3 + xsd] + dvdet * cornerMetrics[ysd * 3 + xsd]
                  + dvdze * cornerMetrics[zsd * 3 + xsd];

    // compute tau3 = du/dz + dw/dx

    // tau_xz = du/dxi * dxi/dz + du/deta * deta/dz + du/dzeta * dzeta/dz
    //        + dw/dxi * dxi/dx + dw/deta * deta/dx + dw/dzeta * dzeta/dx
    MFloat tau3 = dudxi * cornerMetrics[xsd * 3 + zsd] + dudet * cornerMetrics[ysd * 3 + zsd]
                  + dudze * cornerMetrics[zsd * 3 + zsd] +

                  dwdxi * cornerMetrics[xsd * 3 + xsd] + dwdet * cornerMetrics[ysd * 3 + xsd]
                  + dwdze * cornerMetrics[zsd * 3 + xsd];

    // compute tau4 = 2 dv/dy - 2/3 ( du/dx + dv/dy + dw/dz )

    // tau_yy = 4/3*( dv/dxi * dxi/dy + dv/deta * deta/dy + dv/dzeta * dzeta/dy )
    //            - 2/3*( du/dxi * dxi/dx + du/deta * deta/dx + du/dzeta * dzeta/dx)
    //            - 2/3*( dw/dxi * dxi/dz + dw/deta * deta/dz + dw/dzeta * dzeta/dz )
    MFloat tau4 = F4B3
                      * (dvdxi * cornerMetrics[xsd * 3 + ysd] + dvdet * cornerMetrics[ysd * 3 + ysd]
                         + dvdze * cornerMetrics[zsd * 3 + ysd])
                  -

                  F2B3
                      * (dudxi * cornerMetrics[xsd * 3 + xsd] + dudet * cornerMetrics[ysd * 3 + xsd]
                         + dudze * cornerMetrics[zsd * 3 + xsd])
                  -

                  F2B3
                      * (dwdxi * cornerMetrics[xsd * 3 + zsd] + dwdet * cornerMetrics[ysd * 3 + zsd]
                         + dwdze * cornerMetrics[zsd * 3 + zsd]);

    // compute tau5 = dv/dz + dw/dy

    // tau_yz = dv/dxi * dxi/dz + dv/deta * deta/dz + dv/dzeta * dzeta/dz
    //        + dw/dxi * dxi/dy + dw/deta * deta/dy + dw/dzeta * dzeta/dy
    MFloat tau5 = dvdxi * cornerMetrics[xsd * 3 + zsd] + dvdet * cornerMetrics[ysd * 3 + zsd]
                  + dvdze * cornerMetrics[zsd * 3 + zsd] +

                  dwdxi * cornerMetrics[xsd * 3 + ysd] + dwdet * cornerMetrics[ysd * 3 + ysd]
                  + dwdze * cornerMetrics[zsd * 3 + ysd];

    // compute tau6 = 2 dw/dz - 2/3 ( du/dx + dv/dy + dw/dz )

    // tau_zz = 4/3*( dw/dxi * dxi/dz + dw/deta * deta/dz + dw/dzeta * dzeta/dz )
    //            - 2/3*( du/dxi * dxi/dx + du/deta * deta/dx + du/dzeta * dzeta/dx)
    //            - 2/3*( dv/dxi * dxi/dy + dv/deta * deta/dy + dv/dzeta * dzeta/dy)
    MFloat tau6 = F4B3
                      * (dwdxi * cornerMetrics[xsd * 3 + zsd] + dwdet * cornerMetrics[ysd * 3 + zsd]
                         + dwdze * cornerMetrics[zsd * 3 + zsd])
                  -

                  F2B3
                      * (dudxi * cornerMetrics[xsd * 3 + xsd] + dudet * cornerMetrics[ysd * 3 + xsd]
                         + dudze * cornerMetrics[zsd * 3 + xsd])
                  -

                  F2B3
                      * (dvdxi * cornerMetrics[xsd * 3 + ysd] + dvdet * cornerMetrics[ysd * 3 + ysd]
                         + dvdze * cornerMetrics[zsd * 3 + ysd]);


    const MFloat dTdx = dTdxi * cornerMetrics[xsd * 3 + xsd] + dTdet * cornerMetrics[ysd * 3 + xsd]
                        + dTdze * cornerMetrics[zsd * 3 + xsd];

    const MFloat dTdy = dTdxi * cornerMetrics[xsd * 3 + ysd] + dTdet * cornerMetrics[ysd * 3 + ysd]
                        + dTdze * cornerMetrics[zsd * 3 + ysd];

    const MFloat dTdz = dTdxi * cornerMetrics[xsd * 3 + zsd] + dTdet * cornerMetrics[ysd * 3 + zsd]
                        + dTdze * cornerMetrics[zsd * 3 + zsd];

    const MFloat fJac = 1.0 / m_cells->cornerJac[IJK];
    const MFloat mueOverRe = rRe * fJac * (muLamAvg + muTurbAvg);
    tau1 = mueOverRe * tau1 + TKEcorner;
    tau2 *= mueOverRe;
    tau3 *= mueOverRe;
    tau4 = mueOverRe * tau4 + TKEcorner;
    tau5 *= mueOverRe;
    tau6 = mueOverRe * tau6 + TKEcorner;

    const MFloat muCombined = FgammaMinusOne * rRe * fJac * (rPrLam * muLamAvg + rPrTurb * muTurbAvg);

    const MFloat qx = muCombined * dTdx + uAvg * tau1 + vAvg * tau2 + wAvg * tau3;
    const MFloat qy = muCombined * dTdy + uAvg * tau2 + vAvg * tau4 + wAvg * tau5;
    const MFloat qz = muCombined * dTdz + uAvg * tau3 + vAvg * tau5 + wAvg * tau6;

    // efluxes
    eflux[0][IJK] =
        tau1 * cornerMetrics[xsd * 3 + xsd] + tau2 * cornerMetrics[xsd * 3 + ysd] + tau3 * cornerMetrics[xsd * 3 + zsd];

    eflux[1][IJK] =
        tau2 * cornerMetrics[xsd * 3 + xsd] + tau4 * cornerMetrics[xsd * 3 + ysd] + tau5 * cornerMetrics[xsd * 3 + zsd];

    eflux[2][IJK] =
        tau3 * cornerMetrics[xsd * 3 + xsd] + tau5 * cornerMetrics[xsd * 3 + ysd] + tau6 * cornerMetrics[xsd * 3 + zsd];

    eflux[3][IJK] =
        qx * cornerMetrics[xsd * 3 + xsd] + qy * cornerMetrics[xsd * 3 + ysd] + qz * cornerMetrics[xsd * 3 + zsd];

    // ffluxes
    fflux[0][IJK] =
        tau1 * cornerMetrics[ysd * 3 + xsd] + tau2 * cornerMetrics[ysd * 3 + ysd] + tau3 * cornerMetrics[ysd * 3 + zsd];

    fflux[1][IJK] =
        tau2 * cornerMetrics[ysd * 3 + xsd] + tau4 * cornerMetrics[ysd * 3 + ysd] + tau5 * cornerMetrics[ysd * 3 + zsd];

    fflux[2][IJK] =
        tau3 * cornerMetrics[ysd * 3 + xsd] + tau5 * cornerMetrics[ysd * 3 + ysd] + tau6 * cornerMetrics[ysd * 3 + zsd];

    fflux[3][IJK] =
        qx * cornerMetrics[ysd * 3 + xsd] + qy * cornerMetrics[ysd * 3 + ysd] + qz * cornerMetrics[ysd * 3 + zsd];

    // gfluxes
    gflux[0][IJK] =
        tau1 * cornerMetrics[zsd * 3 + xsd] + tau2 * cornerMetrics[zsd * 3 + ysd] + tau3 * cornerMetrics[zsd * 3 + zsd];

    gflux[1][IJK] =
        tau2 * cornerMetrics[zsd * 3 + xsd] + tau4 * cornerMetrics[zsd * 3 + ysd] + tau5 * cornerMetrics[zsd * 3 + zsd];

    gflux[2][IJK] =
        tau3 * cornerMetrics[zsd * 3 + xsd] + tau5 * cornerMetrics[zsd * 3 + ysd] + tau6 * cornerMetrics[zsd * 3 + zsd];

    gflux[3][IJK] =
        qx * cornerMetrics[zsd * 3 + xsd] + qy * cornerMetrics[zsd * 3 + ysd] + qz * cornerMetrics[zsd * 3 + zsd];
  });


  // viscous flux correction for the singular points
  // m_hasSingularity=0 means no singular points in this block, otherwise do flux correction
  if(m_hasSingularity > 0) {
    viscousFluxCorrection();
  }

  for(MInt var = 0; var < 4; var++) {
    const std::array<MInt, nDim> beginIM1{m_noGhostLayers - 1, m_noGhostLayers, m_noGhostLayers};
    maia::parallelFor<true, nDim>(beginIM1, endM2(), [=](const MInt& i, const MInt& j, const MInt& k) {
      const MInt IJK = cellIndex(i, j, k);
      const MInt IJMK = cellIndex(i, (j - 1), k);
      const MInt IJKM = cellIndex(i, j, (k - 1));
      const MInt IJMKM = cellIndex(i, (j - 1), (k - 1));

      vflux[0][IJK] = F1B4 * (eflux[var][IJK] + eflux[var][IJKM] + eflux[var][IJMK] + eflux[var][IJMKM]);
    });


    const std::array<MInt, nDim> beginJM1{m_noGhostLayers, m_noGhostLayers - 1, m_noGhostLayers};
    maia::parallelFor<true, nDim>(beginJM1, endM2(), [=](const MInt& i, const MInt& j, const MInt& k) {
      const MInt IJK = cellIndex(i, j, k);
      const MInt IMJK = cellIndex((i - 1), j, k);
      const MInt IJKM = cellIndex(i, j, (k - 1));
      const MInt IMJKM = cellIndex((i - 1), j, (k - 1));

      vflux[1][IJK] = F1B4 * (fflux[var][IJK] + fflux[var][IJKM] + fflux[var][IMJK] + fflux[var][IMJKM]);
    });

    const std::array<MInt, nDim> beginKM1{m_noGhostLayers, m_noGhostLayers, m_noGhostLayers - 1};
    maia::parallelFor<true, nDim>(beginKM1, endM2(), [=](const MInt& i, const MInt& j, const MInt& k) {
      const MInt IJK = cellIndex(i, j, k);
      const MInt IMJK = cellIndex((i - 1), j, k);
      const MInt IJMK = cellIndex(i, (j - 1), k);
      const MInt IMJMK = cellIndex((i - 1), (j - 1), k);

      vflux[2][IJK] = F1B4 * (gflux[var][IJK] + gflux[var][IMJK] + gflux[var][IJMK] + gflux[var][IMJMK]);
    });

    maia::parallelFor<true, nDim>(beginP2(), endM2(), [=](const MInt& i, const MInt& j, const MInt& k) {
      const MInt IJK = cellIndex(i, j, k);
      const MInt IMJK = cellIndex(i - 1, j, k);
      const MInt IJMK = cellIndex(i, j - 1, k);
      const MInt IJKM = cellIndex(i, j, k - 1);
      m_cells->rightHandSide[var][IJK] +=
          vflux[0][IJK] - vflux[0][IMJK] + vflux[1][IJK] - vflux[1][IJMK] + vflux[2][IJK] - vflux[2][IJKM];
    });
  }
}


void FvStructuredSolver3D::viscousFluxCorrection() {
  const MFloat rPr = F1 / m_Pr;
  const MFloat rRe = F1 / m_Re0;
  const MFloat gammaMinusOne = m_gamma - 1.0;
  const MFloat FgammaMinusOne = F1 / gammaMinusOne;

  MFloat* RESTRICT rhou = &m_cells->variables[CV->RHO_U][0];
  MFloat* RESTRICT rhov = &m_cells->variables[CV->RHO_V][0];
  MFloat* RESTRICT rhow = &m_cells->variables[CV->RHO_W][0];
  MFloat* RESTRICT rhoE = &m_cells->variables[CV->RHO_E][0];
  MFloat* RESTRICT rho = &m_cells->variables[CV->RHO][0];

  MFloat* const* const RESTRICT eflux = ALIGNED_F(m_cells->eFlux);
  MFloat* const* const RESTRICT fflux = ALIGNED_F(m_cells->fFlux);
  MFloat* const* const RESTRICT gflux = ALIGNED_F(m_cells->gFlux);

  MInt dim = 0;
  MInt start[3], end[3], nghbr[20];
  MInt len1[3];
  MInt totalCells;

  for(MInt i = 0; i < m_hasSingularity; ++i) {
    // only correct for bc 6000 not for bc 4000-5000
    if(m_singularity[i].BC == -6000) {
      totalCells = 1;
      for(MInt j = 0; j < nDim; j++) {
        len1[j] = m_singularity[i].end[j] - m_singularity[i].start[j];
        if(len1[j] != 0) totalCells *= len1[j];
      }

      for(MInt n = 0; n < 3; ++n) {
        if(m_singularity[i].end[n] - m_singularity[i].start[n] > 1) {
          dim = n;
          // start[n]=m_singularity[i].start[n]+1;
          start[n] = m_singularity[i].start[n] + 1;
          end[n] = m_singularity[i].end[n] - 1;
        } else {
          start[n] = m_singularity[i].start[n];
          end[n] = m_singularity[i].end[n];
        }
      }

      MFloat u[20], v[20], w[20], T[20];
      MFloat U, V, W, t, dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz, dTdx, dTdy, dTdz;
      MFloat qx, qy, qz;
      MFloat tau1, tau2, tau3, tau4, tau5, tau6;
      MFloat mueOverRe, mue, mueH;
      for(MInt kk = start[2]; kk < end[2]; ++kk) {
        for(MInt jj = start[1]; jj < end[1]; ++jj) {
          for(MInt ii = start[0]; ii < end[0]; ++ii) {
            MInt count = 0;
            MInt temp[3] = {0, 0, 0};
            MInt IJK = cellIndex(ii + m_singularity[i].Viscous[0], jj + m_singularity[i].Viscous[1],
                                 kk + m_singularity[i].Viscous[2]);

            const MFloat cornerMetrics[9] = {
                m_cells->cornerMetrics[0][IJK], m_cells->cornerMetrics[1][IJK], m_cells->cornerMetrics[2][IJK],
                m_cells->cornerMetrics[3][IJK], m_cells->cornerMetrics[4][IJK], m_cells->cornerMetrics[5][IJK],
                m_cells->cornerMetrics[6][IJK], m_cells->cornerMetrics[7][IJK], m_cells->cornerMetrics[8][IJK]};

            temp[dim] = 1;
            nghbr[count++] = cellIndex(ii, jj, kk);
            nghbr[count++] = cellIndex(ii + temp[0], jj + temp[1], kk + temp[2]);

            for(MInt m = 0; m < m_singularity[i].Nstar - 1; ++m) {
              MInt* change = m_singularity[i].displacement[m];
              nghbr[count++] = cellIndex(ii + change[0], jj + change[1], kk + change[2]);
              nghbr[count++] = cellIndex(ii + temp[0] + change[0], jj + temp[1] + change[1], kk + temp[2] + change[2]);
            }

            if(count != m_singularity[i].Nstar * 2) {
              cout << "what the hell! it is wrong!!!" << endl;
            }

            for(MInt m = 0; m < m_singularity[i].Nstar * 2; ++m) {
              u[m] = rhou[nghbr[m]] / rho[nghbr[m]];
              v[m] = rhov[nghbr[m]] / rho[nghbr[m]];
              w[m] = rhow[nghbr[m]] / rho[nghbr[m]];
              T[m] = (m_gamma * gammaMinusOne
                      * (rhoE[nghbr[m]] - F1B2 * rho[nghbr[m]] * (POW2(u[m]) + POW2(v[m]) + POW2(w[m]))))
                     / rho[nghbr[m]];
            }

            U = F0;
            V = F0;
            W = F0;
            t = F0;
            dudx = F0;
            dudy = F0;
            dudz = F0;
            dvdx = F0;
            dvdy = F0;
            dvdz = F0;
            dwdx = F0;
            dwdy = F0;
            dwdz = F0;
            dTdx = F0;
            dTdy = F0;
            dTdz = F0;

            MInt id2 = ii - start[0] + ((jj - start[1]) + (kk - start[2]) * len1[1]) * len1[0];

            for(MInt n = 0; n < count; n++) {
              MInt ID = id2 * count + n;
              U += m_singularity[i].ReconstructionConstants[0][ID] * u[n];
              dudx += m_singularity[i].ReconstructionConstants[1][ID] * u[n];
              dudy += m_singularity[i].ReconstructionConstants[2][ID] * u[n];
              dudz += m_singularity[i].ReconstructionConstants[3][ID] * u[n];

              V += m_singularity[i].ReconstructionConstants[0][ID] * v[n];
              dvdx += m_singularity[i].ReconstructionConstants[1][ID] * v[n];
              dvdy += m_singularity[i].ReconstructionConstants[2][ID] * v[n];
              dvdz += m_singularity[i].ReconstructionConstants[3][ID] * v[n];

              W += m_singularity[i].ReconstructionConstants[0][ID] * w[n];
              dwdx += m_singularity[i].ReconstructionConstants[1][ID] * w[n];
              dwdy += m_singularity[i].ReconstructionConstants[2][ID] * w[n];
              dwdz += m_singularity[i].ReconstructionConstants[3][ID] * w[n];

              t += m_singularity[i].ReconstructionConstants[0][ID] * T[n];
              dTdx += m_singularity[i].ReconstructionConstants[1][ID] * T[n];
              dTdy += m_singularity[i].ReconstructionConstants[2][ID] * T[n];
              dTdz += m_singularity[i].ReconstructionConstants[3][ID] * T[n];
            }

            tau1 = 2 * dudx - 2 / 3 * (dudx + dvdy + dwdz);
            tau2 = dudy + dvdx;
            tau3 = dudz + dwdx;
            tau4 = 2 * dvdy - 2 / 3 * (dudx + dvdy + dwdz);
            tau5 = dvdz + dwdy;
            tau6 = 2 * dwdz - 2 / 3 * (dudx + dvdy + dwdz);

            mue = SUTHERLANDLAW(t);
            mueOverRe = mue * rRe;
            tau1 *= mueOverRe;
            tau2 *= mueOverRe;
            tau3 *= mueOverRe;
            tau4 *= mueOverRe;
            tau5 *= mueOverRe;
            tau6 *= mueOverRe;
            mueH = FgammaMinusOne * mueOverRe * rPr;

            qx = mueH * dTdx + U * tau1 + V * tau2 + W * tau3;
            qy = mueH * dTdy + U * tau2 + V * tau4 + W * tau5;
            qz = mueH * dTdz + U * tau3 + V * tau5 + W * tau6;


            // efluxes
            eflux[0][IJK] = tau1 * cornerMetrics[xsd * 3 + xsd] + tau2 * cornerMetrics[xsd * 3 + ysd]
                            + tau3 * cornerMetrics[xsd * 3 + zsd];

            eflux[1][IJK] = tau2 * cornerMetrics[xsd * 3 + xsd] + tau4 * cornerMetrics[xsd * 3 + ysd]
                            + tau5 * cornerMetrics[xsd * 3 + zsd];

            eflux[2][IJK] = tau3 * cornerMetrics[xsd * 3 + xsd] + tau5 * cornerMetrics[xsd * 3 + ysd]
                            + tau6 * cornerMetrics[xsd * 3 + zsd];

            eflux[3][IJK] = qx * cornerMetrics[xsd * 3 + xsd] + qy * cornerMetrics[xsd * 3 + ysd]
                            + qz * cornerMetrics[xsd * 3 + zsd];

            // ffluxes
            fflux[0][IJK] = tau1 * cornerMetrics[ysd * 3 + xsd] + tau2 * cornerMetrics[ysd * 3 + ysd]
                            + tau3 * cornerMetrics[ysd * 3 + zsd];

            fflux[1][IJK] = tau2 * cornerMetrics[ysd * 3 + xsd] + tau4 * cornerMetrics[ysd * 3 + ysd]
                            + tau5 * cornerMetrics[ysd * 3 + zsd];

            fflux[2][IJK] = tau3 * cornerMetrics[ysd * 3 + xsd] + tau5 * cornerMetrics[ysd * 3 + ysd]
                            + tau6 * cornerMetrics[ysd * 3 + zsd];

            fflux[3][IJK] = qx * cornerMetrics[ysd * 3 + xsd] + qy * cornerMetrics[ysd * 3 + ysd]
                            + qz * cornerMetrics[ysd * 3 + zsd];

            // gfluxes
            gflux[0][IJK] = tau1 * cornerMetrics[zsd * 3 + xsd] + tau2 * cornerMetrics[zsd * 3 + ysd]
                            + tau3 * cornerMetrics[zsd * 3 + zsd];

            gflux[1][IJK] = tau2 * cornerMetrics[zsd * 3 + xsd] + tau4 * cornerMetrics[zsd * 3 + ysd]
                            + tau5 * cornerMetrics[zsd * 3 + zsd];

            gflux[2][IJK] = tau3 * cornerMetrics[zsd * 3 + xsd] + tau5 * cornerMetrics[zsd * 3 + ysd]
                            + tau6 * cornerMetrics[zsd * 3 + zsd];

            gflux[3][IJK] = qx * cornerMetrics[zsd * 3 + xsd] + qy * cornerMetrics[zsd * 3 + ysd]
                            + qz * cornerMetrics[zsd * 3 + zsd];
          }
        }
      }
    }
  }
}


void FvStructuredSolver3D::computePorousRHS(MBool isRans) {
  TRACE();

  const MFloat rRe0 = 1.0 / m_Re0;

  const MFloat* const* const RESTRICT pvars = m_cells->pvariables;
  const MFloat* const RESTRICT muLam = &m_cells->fq[FQ->MU_L][0];
  const MFloat* const RESTRICT por = &m_cells->fq[FQ->POROSITY][0];
  const MFloat* const RESTRICT Da = &m_cells->fq[FQ->DARCY][0];
  const MFloat* const RESTRICT cf = &m_cells->fq[FQ->FORCH][0];

  if(isRans) {
    mTerm(1, "Porous stuff for 3D RANS not implemented yet!");
  } else { // LES
    for(MInt k = m_noGhostLayers; k < m_nCells[0] - m_noGhostLayers; ++k) {
      for(MInt j = m_noGhostLayers; j < m_nCells[1] - m_noGhostLayers; ++j) {
        for(MInt i = m_noGhostLayers; i < m_nCells[2] - m_noGhostLayers; ++i) {
          const MInt IJK = cellIndex(i, j, k);
          MFloat velAbs = F0;
          for(MInt dim = 0; dim < nDim; ++dim) {
            velAbs += POW2(pvars[PV->VV[dim]][IJK]);
          }
          velAbs = sqrt(velAbs);
          const MFloat rDa = 1.0 / Da[IJK];
          const MFloat porPOW2 = POW2(por[IJK]);
          for(MInt dim = 0; dim < nDim; ++dim) {
            m_cells->rightHandSide[CV->RHO_VV[dim]][IJK] +=
                -(rRe0 * rDa * por[IJK] * muLam[IJK] + porPOW2 * sqrt(rDa) * cf[IJK] * velAbs * pvars[PV->RHO][IJK])
                * pvars[dim /*PV->VV[dim]  <-- why is this giving linking error???*/][IJK] * m_cells->cellJac[IJK];
          }
        }
      }
    }
  }
}


/**
 * \brief Computation of the maximum residual
 *
 * This function computes the maxResidual using
 * Res = deltaT/(CFL*VolOfCell) * |RHS|
 *
 * with deltaT depending on local or global time stepping
 * is used.
 * checks if the computed max density residual
 * is below the convergence criterion and returns
 * boolean variable
 *
 * \author Pascal Meysonnat
 *
 */

MBool FvStructuredSolver3D::maxResidual() {
  TRACE();

  if(globalTimeStep % m_residualInterval != 0) return true;
  MFloat epsilon = pow(10.0, -10.0);
  m_avrgResidual = F0;
  MInt cellId = F0;
  MFloat tmpResidual = F0;
  MFloat maxResidual1 = F0;
  MInt maxResIndex[3];
  // MInt localCounter=F0;
  MFloat maxResidualOrg = F0;
  MFloat localMaxResidual = F0;
  MFloat localAvrgResidual = F0;
  MFloat accumAvrgResidual = F0;
  MFloat globalMaxResidual = F0;
  m_workload = F0;
  // MInt accumCounter=0;
  for(MInt dim = 0; dim < nDim; dim++) {
    maxResIndex[dim] = F0;
  }

  for(MInt k = m_noGhostLayers; k < m_nCells[0] - m_noGhostLayers; k++) {
    for(MInt j = m_noGhostLayers; j < m_nCells[1] - m_noGhostLayers; j++) {
      for(MInt i = m_noGhostLayers; i < m_nCells[2] - m_noGhostLayers; i++) {
        cellId = cellIndex(i, j, k);
        // cerr << cellId << endl;
        tmpResidual = m_timeStep / (m_cfl * m_cells->cellJac[cellId]) * fabs(m_cells->rightHandSide[CV->RHO][cellId]);
        m_avrgResidual += tmpResidual;

        if(tmpResidual > maxResidual1) {
          maxResIndex[0] = i - m_noGhostLayers;
          maxResIndex[1] = j - m_noGhostLayers;
          maxResIndex[2] = k - m_noGhostLayers;
          maxResidual1 = tmpResidual;
        }
      }
    }
  }

  // localCounter = counter;
  localMaxResidual = maxResidual1;
  localAvrgResidual = m_avrgResidual;
  // reset average Residual
  m_avrgResidual = F0;

  MFloat localTotalEnergy = F0;
  MFloat globalTotalEnergy = F0;
  MFloat globalPressure = F0;
  MFloat localPressure = F0;
  if(m_initialCondition == 101) {
    localTotalEnergy = computeTotalKineticEngergy();
    MPI_Allreduce(&localTotalEnergy, &globalTotalEnergy, 1, MPI_DOUBLE, MPI_SUM, m_StructuredComm, AT_,
                  "localTotalEnergy", "globalTotalEnergy");
    globalTotalEnergy /= ((16.0 * pow(4 * atan(1), 3.0))); // divided by the overall volume
    localPressure = computeTotalPressure();
    MPI_Allreduce(&localPressure, &globalPressure, 1, MPI_DOUBLE, MPI_SUM, m_StructuredComm, AT_, "localPressure",
                  "globalPressure");
    globalPressure /= ((16.0 * pow(4 * atan(1), 3.0))); // divided by the overall volume
  }

  // if( noDomains()>1 )
  //  {
  // MPI_Allreduce(m_residualSnd, &m_residualRcv, 1, m_mpiStruct, m_resOp, m_StructuredComm, AT_, "m_residualSnd",
  // "m_residualRcv" );
  MPI_Allreduce(&localAvrgResidual, &accumAvrgResidual, 1, MPI_DOUBLE, MPI_SUM, m_StructuredComm, AT_,
                "localAvrgResidual", "accumAvrgResidual");
  MPI_Allreduce(&localMaxResidual, &globalMaxResidual, 1, MPI_DOUBLE, MPI_MAX, m_StructuredComm, AT_,
                "localMaxResidual", "globalMaxResidual");

  m_avrgResidual = accumAvrgResidual; // m_residualRcv.avrgRes;
  maxResidualOrg = globalMaxResidual;
  // globalMaxResidual=globalMaxResidual;//m_residualRcv.maxRes;
  // for(MInt i=0; i<3; i++)
  //{
  //  maxResIndex[i]=m_residualRcv.maxCellIndex[i];
  //}

  // }

  // cout << "m_avrgResidual = " << m_avrgResidual<< " | totalCells " << m_totalGridCells <<endl;
  m_avrgResidual = m_avrgResidual / m_totalNoCells;
  // write first residuals;
  if(ABS(m_firstMaxResidual) < epsilon) {
    m_firstMaxResidual = mMax(epsilon, globalMaxResidual);
    m_firstAvrgResidual = mMax(epsilon, m_avrgResidual);
    if(m_initialCondition != 101) { // we need an extra treatment of TGV because of symmetie
      if(approx(localMaxResidual, maxResidualOrg,
                m_eps)) { // so only domainId with the max writes out ==> no need to communicate the max index[i]

        // write out values into residual file
        FILE* f_residual;
        f_residual = fopen("./Residual", "a+");
        fprintf(f_residual, "#MaxRes_1: %1.10e \n", m_firstMaxResidual);
        fprintf(f_residual, "#MaxAvgRes_1: %1.10e \n", m_firstAvrgResidual);
        fprintf(f_residual, "#iter, physTime, time, dT, wLoad, avrgRes, maxRes, blockId, i, j, k ");
        fclose(f_residual);
      }
    } else {
      if(domainId() == 0) {
        // write out values into residual file
        FILE* f_residual;
        f_residual = fopen("./Residual", "a+");
        fprintf(f_residual, "#MaxRes_1: %1.10e \n", m_firstMaxResidual);
        fprintf(f_residual, "#MaxAvgRes_1: %1.10e \n", m_firstAvrgResidual);
        fprintf(f_residual,
                "#iter, physTime, time, dT, wLoad, avrgRes, maxRes, blockId, i, j, k, k_mean, p_mean, pProbe ");
      }
    }
  }

  // normalize residuals
  globalMaxResidual = globalMaxResidual / m_firstMaxResidual;
  m_avrgResidual = (m_avrgResidual / m_firstAvrgResidual);

  // question if "( m_avrgResidual >= F0 || m_avrgResidual < F0 ) {} else {"
  // is better to capture the also inf???

  if(std::isnan(m_avrgResidual)) {
    cerr << "Solution diverged, average residual is nan " << endl;
    m_log << "Solution diverged, average residual is nan " << endl;
    saveSolverSolution(true);
    mTerm(1, AT_, "Solution diverged, average residual is nan ");
  }

  // convergence Check

  m_convergence = false;
  if(maxResidual1 < m_convergenceCriterion) {
    m_convergence = true;
  }

  // need again special treatment for TGV due to symmetrie many processors would write out (!!!!IMPORTANT NO CORRECT
  // INDEX WILL BE WRITTEN OUT )
  if(m_initialCondition != 101) {
    // processor with the highest Residual writes out!!! saves communication;
    if(approx(localMaxResidual, maxResidualOrg,
              m_eps)) { // so only domainId with the max writes out ==> no need to communicate the max index[i]
      // write out values into residual file
      FILE* f_residual;
      f_residual = fopen("./Residual", "a+");
      fprintf(f_residual, "%d", globalTimeStep);
      fprintf(f_residual, " %f", m_physicalTime);
      fprintf(f_residual, " %f", m_time);
      fprintf(f_residual, " %f", m_timeStep);
      fprintf(f_residual, " %f", m_workload);
      fprintf(f_residual, " %1.10e", m_avrgResidual);
      fprintf(f_residual, " %1.10e", globalMaxResidual);
      fprintf(f_residual, " %d", m_blockId);
      fprintf(f_residual, " %d", m_nOffsetCells[2] + maxResIndex[0]); // i
      fprintf(f_residual, " %d", m_nOffsetCells[1] + maxResIndex[1]); // j
      fprintf(f_residual, " %d", m_nOffsetCells[0] + maxResIndex[2]); // k
      fprintf(f_residual, "\n");
      fclose(f_residual);
    }
  } else {
    if(domainId() == 0) {
      MFloat dissip = F0;
      if(globalTimeStep == 1) {
        m_kineticEOld = 1.2474901617e-03;
      } else {
        dissip = (globalTotalEnergy - m_kineticEOld) / m_timeStep;
        m_kineticEOld = globalTotalEnergy;
      }
      // write out values into residual file
      // compute the dissipation rate

      FILE* f_residual;
      f_residual = fopen("./Residual", "a+");
      fprintf(f_residual, "%d", globalTimeStep);
      fprintf(f_residual, " %f", m_physicalTime);
      fprintf(f_residual, " %f", m_time);
      fprintf(f_residual, " %f", m_timeStep);
      fprintf(f_residual, " %f", m_workload);
      fprintf(f_residual, " %1.10e", m_avrgResidual);
      fprintf(f_residual, " %1.10e", globalMaxResidual);
      fprintf(f_residual, " %d", m_blockId);
      fprintf(f_residual, " %d", m_nOffsetCells[2] + maxResIndex[0]); // i Will be wrong
      fprintf(f_residual, " %d", m_nOffsetCells[1] + maxResIndex[1]); // j Will be wrong
      fprintf(f_residual, " %d", m_nOffsetCells[0] + maxResIndex[2]); // k Will be wrong
      fprintf(f_residual, " %1.10e", globalTotalEnergy);              // kinetic Energy
      fprintf(f_residual, " %1.10e", globalPressure);                 // averaged pressure
      fprintf(f_residual, " %1.10e", dissip);                         // dissipation rate
      fprintf(f_residual, "\n");
      fclose(f_residual);
    }
  }

  if(maxResidual1 < m_convergenceCriterion) {
    return true;
  } else {
    return false;
  }
}

MFloat FvStructuredSolver3D::computeTotalKineticEngergy() {
  TRACE();
  MFloat localEnergy = F0;
  for(MInt k = m_noGhostLayers; k < m_nCells[0] - m_noGhostLayers; k++) {
    for(MInt j = m_noGhostLayers; j < m_nCells[1] - m_noGhostLayers; j++) {
      for(MInt i = m_noGhostLayers; i < m_nCells[2] - m_noGhostLayers; i++) {
        const MInt cellId = cellIndex(i, j, k);
        localEnergy += (POW2(m_cells->pvariables[PV->U][cellId]) + POW2(m_cells->pvariables[PV->V][cellId])
                        + POW2(m_cells->pvariables[PV->W][cellId]))
                       * m_cells->cellJac[cellId];
      }
    }
  }
  return localEnergy;
}

MFloat FvStructuredSolver3D::computeTotalPressure() {
  TRACE();
  MFloat localPressure = F0;
  for(MInt k = m_noGhostLayers; k < m_nCells[0] - m_noGhostLayers; k++) {
    for(MInt j = m_noGhostLayers; j < m_nCells[1] - m_noGhostLayers; j++) {
      for(MInt i = m_noGhostLayers; i < m_nCells[2] - m_noGhostLayers; i++) {
        const MInt cellId = cellIndex(i, j, k);
        localPressure += m_cells->pvariables[PV->P][cellId] * m_cells->cellJac[cellId];
      }
    }
  }
  return localPressure;
}

inline MFloat FvStructuredSolver3D::pressure(MInt cellId) { return m_cells->pvariables[PV->P][cellId]; }


inline void FvStructuredSolver3D::crossProduct(MFloat* result, const MFloat* vec1, const MFloat* vec2) {
  result[xsd] = vec1[ysd] * vec2[zsd] - vec1[zsd] * vec2[ysd];
  result[ysd] = vec1[zsd] * vec2[xsd] - vec1[xsd] * vec2[zsd];
  result[zsd] = vec1[xsd] * vec2[ysd] - vec1[ysd] * vec2[xsd];
}

inline MInt FvStructuredSolver3D::getPointIdFromCell(const MInt i, const MInt j, const MInt k) {
  return i + (k * (m_nCells[1] + 1) + j) * (m_nCells[2] + 1);
}

inline MInt FvStructuredSolver3D::getPointIdfromPoint(const MInt origin, const MInt incI, const MInt incJ,
                                                      const MInt incK) {
  return origin + incI + incJ * m_nPoints[2] + incK * m_nPoints[2] * m_nPoints[1];
}


void FvStructuredSolver3D::moveGrid(const MBool isRestart, const MBool zeroPos) {
  TRACE();
  const MFloat pi = 4.0 * atan(1);
  const MFloat t = (isRestart) ? m_time : m_time + m_timeStep * m_RKalpha[m_RKStep];

  switch(m_gridMovingMethod) {
    case 1: // oscillating grid in x-direction
    {
      // we need some relaxation function in wall-normal direction
      // here we take a linear function: ( y_max - y ) / y_max
      const MFloat y_max = 10.0; // y of upper domain boundary
      const MFloat frequency = F1 / (F2 * m_timeStep * 1000.0);
      const MFloat amp = PV->UInfinity / (F2 * pi * frequency);

      for(MInt k = 0; k < m_nPoints[0]; k++) {
        for(MInt j = 0; j < (m_nPoints[1]); j++) {
          for(MInt i = 0; i < m_nPoints[2]; i++) {
            const MInt pointId = pointIndex(i, j, k);
            const MFloat x = m_grid->m_initCoordinates[0][pointId];

            m_grid->m_coordinates[0][pointId] =
                x + (1 - m_grid->m_coordinates[1][pointId] / y_max) * amp * sin(2 * pi * t * frequency);
            m_grid->m_velocity[0][pointId] = (1 - m_grid->m_coordinates[1][pointId] / y_max) * amp * 2 * pi * frequency
                                             * cos(2 * pi * t * frequency);
          }
        }
      }

      break;
    }
    case 2: // channel with moving indentation
    {
      const MFloat beta = 4.14l;
      const MFloat StrNum = m_wallVel; // Strouhal Number, set by wallVel
      const MFloat ver = F0;           // ver can be used to translate x coordinate
      const MFloat y_max = F1;         // can be used to scale (old mesh had y_max = m_mgOrgCoordinates[0][32] = 1/30)

      // values from Ralph and Pedley:
      const MFloat x2 = y_max * (ver - 11.75l);
      const MFloat x3 = y_max * (ver - 9.25l);
      const MFloat x4 = y_max * (ver - 1.25);
      const MFloat x5 = y_max * (ver + 1.25);
      const MFloat omega = PV->UInfinity * StrNum * F2 * pi / y_max;
      const MFloat h = F1B2 * 0.38l * (F1 - cos(omega * t));
      const MFloat hvel = F1B2 * 0.38l * omega * sin(omega * t);

      for(MInt k = 0; k < m_nPoints[0]; k++) {
        for(MInt j = 0; j < m_nPoints[1]; j++) {
          for(MInt i = 0; i < m_nPoints[2]; i++) {
            const MInt pointId = pointIndex(i, j, k);
            const MFloat y = m_grid->m_initCoordinates[1][pointId];

            MFloat g = F0;

            if(m_grid->m_coordinates[0][pointId] > x2 && m_grid->m_coordinates[0][pointId] < x3) {
              g = (1.0 + tanhl(beta * (m_grid->m_coordinates[0][pointId] - (x2 + x3) / 2.0l) / y_max)) / 2.0;
            } else if(m_grid->m_coordinates[0][pointId] >= x3 && m_grid->m_coordinates[0][pointId] < x4) {
              g = F1;
            } else if(m_grid->m_coordinates[0][pointId] >= x4 && m_grid->m_coordinates[0][pointId] < x5) {
              g = (1.0 - tanhl(beta * (m_grid->m_coordinates[0][pointId] - (x4 + x5) / 2.0l) / y_max)) / 2.0;
            }

            m_grid->m_coordinates[1][pointId] = y * (F1 - h * g);
            m_grid->m_velocity[1][pointId] = -y * hvel * g;
          }
        }
      }

      break;
    }
    case 3: // piston moving in x-direction
    {
      for(MInt k = 0; k < m_nPoints[0]; k++) {
        for(MInt j = 0; j < m_nPoints[1]; j++) {
          for(MInt i = 0; i < m_nPoints[2]; i++) {
            const MInt pointId = pointIndex(i, j, k);
            const MFloat x = m_grid->m_initCoordinates[0][pointId];

            m_grid->m_coordinates[0][pointId] = x * (1 + t * m_wallVel);
            m_grid->m_velocity[0][pointId] = x * m_wallVel;
            m_grid->m_acceleration[0][pointId] = F0;
          }
        }
      }

      break;
    }
    case 4: // inner grid movement
    {
      const MFloat beta = 16.0l;
      const MFloat StrNum = m_wallVel; // Strouhal Number, set by wallVel
      const MFloat ver = 0.0l;         // ver can be used to translate x coordinate
      const MFloat y_max = 1.0l;       // can be used to scale

      // for Square:
      const MFloat x2 = y_max * (ver + 0.1l);
      const MFloat x3 = y_max * (ver + 0.5l);
      const MFloat x4 = y_max * (ver + 0.5l);
      const MFloat x5 = y_max * (ver + 0.9l);

      const MFloat omega = m_Ma * sqrt(PV->TInfinity) * StrNum * 2.0l * pi / y_max;
      const MFloat h = F1B2 * 0.35l * (F1 - cos(omega * t));
      const MFloat hvel = F1B2 * 0.35l * (omega * sin(omega * t));

      for(MInt k = 0; k < m_nPoints[0]; k++) {
        for(MInt j = 0; j < m_nPoints[1]; j++) {
          for(MInt i = 0; i < m_nPoints[2]; i++) {
            const MInt pointId = pointIndex(i, j, k);

            const MFloat x = m_grid->m_initCoordinates[0][pointId];
            const MFloat y = m_grid->m_initCoordinates[1][pointId];

            MFloat g = F0;
            if(y > x2 && y < x5 && x > x2 && x < x5) {
              g = ((y < x3) ? (1.0l + tanhl(beta * (y - (x2 + x3) / 2.0l) / y_max)) / 2.0l
                            : ((y < x4) ? 1.0l : (1.0l - tanhl(beta * (y - (x4 + x5) / 2.0l) / y_max)) / 2.0l));
            }
            m_grid->m_coordinates[0][pointId] = x * (F1 - h * g * (F1 - x));
            m_grid->m_velocity[0][pointId] = x * (-hvel * g * (F1 - x));

            g = F0;
            if(x > x2 && x < x5 && y > x2 && y < x5) {
              g = ((x < x3) ? (1.0l + tanhl(beta * (x - (x2 + x3) / 2.0l) / y_max)) / 2.0l
                            : ((x < x4) ? 1.0l : (1.0l - tanhl(beta * (x - (x4 + x5) / 2.0l) / y_max)) / 2.0l));
            }
            m_grid->m_coordinates[1][pointId] = y * (F1 - h * g * (F1 - y));
            m_grid->m_velocity[1][pointId] = y * (-hvel * g * (F1 - y));
          }
        }
      }

      break;
    }
    case 9:
    case 10: {
      // traveling wave case
      const MFloat angle = m_waveAngle;
      const MFloat rotSin = sin(angle);
      MFloat t_offset = t - m_movingGridTimeOffset;
      if(zeroPos) {
        t_offset = F0;
      }
      const MFloat transitionLength = m_waveEndTransition - m_waveBeginTransition;
      const MFloat transitionOutLength = m_waveOutEndTransition - m_waveOutBeginTransition;
      const MFloat yTransitionLength = m_waveYEndTransition - m_waveYBeginTransition;

      MFloat fadeInFactor = 0.0;
      MFloat fadeInFactorPrime = 0.0;
      MFloat fadeInFactorPrimePrime = 0.0;
      if(t_offset < m_waveTemporalTransition) {
        fadeInFactor = (1.0 - cos(t_offset / m_waveTemporalTransition * pi)) * F1B2;
        fadeInFactorPrime = (pi / m_waveTemporalTransition) * sin(t_offset / m_waveTemporalTransition * pi) * F1B2;
        fadeInFactorPrimePrime =
            POW2(pi / m_waveTemporalTransition) * cos(t_offset / m_waveTemporalTransition * pi) * F1B2;
      } else {
        fadeInFactor = 1.0;
        fadeInFactorPrime = 0.0;
        fadeInFactorPrimePrime = 0.0;
      }

      if(zeroPos) {
        fadeInFactor = F1;
      }


      for(MInt k = 0; k < m_nPoints[0]; k++) {
        for(MInt j = 0; j < m_nPoints[1]; j++) {
          for(MInt i = 0; i < m_nPoints[2]; i++) {
            const MInt pointId = pointIndex(i, j, k);
            const MFloat xInit = m_grid->m_initCoordinates[0][pointId];
            const MFloat zInit = m_grid->m_initCoordinates[2][pointId];
            const MFloat yInit = m_grid->m_initCoordinates[1][pointId];

            MFloat transitionFactor = F0;
            if(xInit <= m_waveBeginTransition) {
              transitionFactor = F0;
            } else if(xInit > m_waveBeginTransition && xInit < m_waveEndTransition) {
              transitionFactor = (1 - cos((xInit - m_waveBeginTransition) / transitionLength * pi)) * F1B2;
            } else if(m_waveEndTransition <= xInit && xInit <= m_waveOutBeginTransition) {
              transitionFactor = F1;
            } else if(xInit > m_waveOutBeginTransition && xInit < m_waveOutEndTransition) {
              transitionFactor = (1 + cos((xInit - m_waveOutBeginTransition) / transitionOutLength * pi)) * F1B2;
            } else {
              transitionFactor = F0;
            }

            MFloat yTransitionFactor = F1;
            if(yInit <= m_waveYBeginTransition) {
              yTransitionFactor = F1;
            } else if(yInit > m_waveYBeginTransition && yInit < m_waveYEndTransition) {
              yTransitionFactor = (1 + cos((yInit - m_waveYBeginTransition) / yTransitionLength * pi)) * F1B2;
            } else {
              yTransitionFactor = F0;
            }

            const MFloat zPrime = zInit - rotSin * (xInit - m_waveBeginTransition);

            const MFloat func = transitionFactor * yTransitionFactor
                                * (m_waveAmplitude * cos((F2 * pi) / m_waveLength * (zPrime - m_waveSpeed * t_offset)));
            const MFloat funcPrime = transitionFactor * yTransitionFactor * (2 * PI * m_waveSpeed / m_waveLength)
                                     * m_waveAmplitude
                                     * sin((F2 * pi) / m_waveLength * (zPrime - m_waveSpeed * t_offset));
            const MFloat funcPrimePrime = -transitionFactor * yTransitionFactor
                                          * POW2(2 * PI * m_waveSpeed / m_waveLength) * m_waveAmplitude
                                          * cos((F2 * pi) / m_waveLength * (zPrime - m_waveSpeed * t_offset));

            m_grid->m_coordinates[1][pointId] = func * fadeInFactor + yInit;
            m_grid->m_velocity[1][pointId] = func * fadeInFactorPrime + funcPrime * fadeInFactor;
            m_grid->m_acceleration[1][pointId] =
                funcPrimePrime * fadeInFactor + 2 * funcPrime * fadeInFactorPrime + func * fadeInFactorPrimePrime;
          }
        }
      }

      break;
    }
    case 11: {
      // traveling wave channel (Tomiyama & Fukagata 2013)
      MFloat t_offset = t - m_movingGridTimeOffset;
      if(zeroPos) {
        t_offset = F0;
      }

      MFloat fadeInFactor = 0;
      const MFloat timeRelaxation = 50.0;

      if(t_offset < timeRelaxation) {
        fadeInFactor = (1.0 - cos(t_offset / timeRelaxation * pi)) * F1B2;
      } else {
        fadeInFactor = 1.0;
      }

      if(zeroPos) {
        fadeInFactor = F1;
      }

      for(MInt k = 0; k < m_nPoints[0]; k++) {
        for(MInt j = 0; j < m_nPoints[1]; j++) {
          for(MInt i = 0; i < m_nPoints[2]; i++) {
            const MInt pointId = pointIndex(i, j, k);
            const MFloat yInit = m_grid->m_initCoordinates[1][pointId];
            const MFloat zInit = m_grid->m_initCoordinates[2][pointId];
            MFloat yRelaxation = F0;

            if(yInit <= F0) {
              yRelaxation = F1;
            } else if(yInit > F0 && yInit < F1) {
              yRelaxation = F1 - yInit;
            } else if(yInit > F1 && yInit < F2) {
              yRelaxation = yInit - F1;
            } else {
              yRelaxation = F1;
            }

            if(yInit <= F1) {
              m_grid->m_coordinates[1][pointId] =
                  yInit
                  + fadeInFactor
                        * (m_waveAmplitude * yRelaxation
                           * cos((F2 * pi) / m_waveLength * (zInit - m_waveSpeed * t_offset)));
              m_grid->m_velocity[1][pointId] = fadeInFactor
                                               * (m_waveAmplitude * yRelaxation * F2 * pi / m_waveLength * m_waveSpeed
                                                  * sin((F2 * pi) / m_waveLength * (zInit - m_waveSpeed * t_offset)));
            } else {
              m_grid->m_coordinates[1][pointId] =
                  yInit
                  - fadeInFactor
                        * (m_waveAmplitude * yRelaxation
                           * cos((F2 * pi) / m_waveLength * (zInit - m_waveSpeed * t_offset)));
              m_grid->m_velocity[1][pointId] = -fadeInFactor
                                               * (m_waveAmplitude * yRelaxation * F2 * pi / m_waveLength * m_waveSpeed
                                                  * sin((F2 * pi) / m_waveLength * (zInit - m_waveSpeed * t_offset)));
            }
          }
        }
      }

      break;
    }
    case 12: {
      // streamwise traveling wave case
      MFloat t_offset = t - m_movingGridTimeOffset;
      if(zeroPos) {
        t_offset = F0;
      }
      const MFloat transitionLength = m_waveEndTransition - m_waveBeginTransition;
      const MFloat transitionOutLength = m_waveOutEndTransition - m_waveOutBeginTransition;
      const MFloat yTransitionLength = m_waveYEndTransition - m_waveYBeginTransition;

      MFloat fadeInFactor = 0.0;
      MFloat fadeInFactorPrime = 0.0;
      MFloat fadeInFactorPrimePrime = 0.0;
      if(t_offset < m_waveTemporalTransition) {
        fadeInFactor = (1.0 - cos(t_offset / m_waveTemporalTransition * pi)) * F1B2;
        fadeInFactorPrime = (pi / m_waveTemporalTransition) * sin(t_offset / m_waveTemporalTransition * pi) * F1B2;
        fadeInFactorPrimePrime =
            POW2(pi / m_waveTemporalTransition) * cos(t_offset / m_waveTemporalTransition * pi) * F1B2;
      } else {
        fadeInFactor = 1.0;
        fadeInFactorPrime = 0.0;
        fadeInFactorPrimePrime = 0.0;
      }

      if(zeroPos) {
        fadeInFactor = F1;
      }


      for(MInt k = 0; k < m_nPoints[0]; k++) {
        for(MInt j = 0; j < m_nPoints[1]; j++) {
          for(MInt i = 0; i < m_nPoints[2]; i++) {
            const MInt pointId = pointIndex(i, j, k);
            const MFloat xInit = m_grid->m_initCoordinates[0][pointId];
            const MFloat yInit = m_grid->m_initCoordinates[1][pointId];

            MFloat transitionFactor = F0;
            if(xInit <= m_waveBeginTransition) {
              transitionFactor = F0;
            } else if(xInit > m_waveBeginTransition && xInit < m_waveEndTransition) {
              transitionFactor = (1 - cos((xInit - m_waveBeginTransition) / transitionLength * pi)) * F1B2;
            } else if(m_waveEndTransition <= xInit && xInit <= m_waveOutBeginTransition) {
              transitionFactor = F1;
            } else if(xInit > m_waveOutBeginTransition && xInit < m_waveOutEndTransition) {
              transitionFactor = (1 + cos((xInit - m_waveOutBeginTransition) / transitionOutLength * pi)) * F1B2;
            } else {
              transitionFactor = F0;
            }

            MFloat yTransitionFactor = F1;
            if(yInit <= m_waveYBeginTransition) {
              yTransitionFactor = F1;
            } else if(yInit > m_waveYBeginTransition && yInit < m_waveYEndTransition) {
              yTransitionFactor = (1 + cos((yInit - m_waveYBeginTransition) / yTransitionLength * pi)) * F1B2;
            } else {
              yTransitionFactor = F0;
            }

            const MFloat func = transitionFactor * yTransitionFactor
                                * (m_waveAmplitude * cos((F2 * pi) / m_waveLength * (xInit - m_waveSpeed * t_offset)));
            const MFloat funcPrime = transitionFactor * yTransitionFactor * (2 * PI * m_waveSpeed / m_waveLength)
                                     * m_waveAmplitude
                                     * sin((F2 * pi) / m_waveLength * (xInit - m_waveSpeed * t_offset));
            const MFloat funcPrimePrime = -transitionFactor * yTransitionFactor
                                          * POW2(2 * PI * m_waveSpeed / m_waveLength) * m_waveAmplitude
                                          * cos((F2 * pi) / m_waveLength * (xInit - m_waveSpeed * t_offset));

            m_grid->m_coordinates[1][pointId] = func * fadeInFactor + yInit;
            m_grid->m_velocity[1][pointId] = func * fadeInFactorPrime + funcPrime * fadeInFactor;
            m_grid->m_acceleration[1][pointId] =
                funcPrimePrime * fadeInFactor + 2 * funcPrime * fadeInFactorPrime + func * fadeInFactorPrimePrime;
          }
        }
      }

      break;
    }
    case 13: {
      // traveling wave on airfoil
      MFloat t_offset = t - m_movingGridTimeOffset;
      if(zeroPos) {
        t_offset = F0;
      }
      const MFloat transitionLength = m_waveEndTransition - m_waveBeginTransition;
      const MFloat transitionOutLength = m_waveOutEndTransition - m_waveOutBeginTransition;
      const MFloat yTransitionLength = m_waveYEndTransition - m_waveYBeginTransition;

      MFloat fadeInFactor = 0;
      const MFloat timeRelaxation = 1.0; // 80.0

      if(t_offset < timeRelaxation) {
        fadeInFactor = (1.0 - cos(t_offset / timeRelaxation * pi)) * F1B2;
      } else {
        fadeInFactor = 1.0;
      }

      if(zeroPos) {
        fadeInFactor = F1;
      }

      const MInt myBlockId = m_grid->getMyBlockId();

      if(myBlockId == 0) {
        for(MInt i = 0; i < m_nPoints[2]; i++) {
          MFloat transitionFactor = F0;
          const MInt offset = m_nOffsetCells[2];
          const MInt globalPointId = offset + i;
          const MFloat xWall = m_airfoilCoords[0 + globalPointId];
          const MFloat yWall = m_airfoilCoords[m_airfoilNoWallPoints + globalPointId];

          const MFloat wallNormalVec[3] = {m_airfoilNormalVec[0 + globalPointId],
                                           m_airfoilNormalVec[m_airfoilNoWallPoints + globalPointId],
                                           m_airfoilNormalVec[2 * m_airfoilNoWallPoints + globalPointId]};

          if(xWall <= m_waveBeginTransition) {
            transitionFactor = F0;
          } else if(xWall > m_waveBeginTransition && xWall < m_waveEndTransition) {
            transitionFactor = (1 - cos((xWall - m_waveBeginTransition) / transitionLength * pi)) * F1B2;
          } else if(m_waveEndTransition <= xWall && xWall <= m_waveOutBeginTransition) {
            transitionFactor = F1;
          } else if(xWall > m_waveOutBeginTransition && xWall < m_waveOutEndTransition) {
            transitionFactor = (1 + cos((xWall - m_waveOutBeginTransition) / transitionOutLength * pi)) * F1B2;
          } else {
            transitionFactor = F0;
          }

          // linear increase of amplitude
          MFloat heightGrad = m_waveGradientSuction / (m_waveOutBeginTransition - m_waveEndTransition);
          MFloat inc = 1.0 + heightGrad * (xWall - m_waveEndTransition);
          MFloat amp0 = m_waveAmplitudeSuction;

          if(wallNormalVec[1] < 0.0) {
            heightGrad = m_waveGradientPressure / (m_waveOutBeginTransition - m_waveEndTransition);
            inc = 1.0 + heightGrad * (xWall - m_waveEndTransition);
            amp0 = m_waveAmplitudePressure;
          }

          for(MInt k = 0; k < m_nPoints[0]; k++) {
            for(MInt j = 0; j < m_nPoints[1] - 1; j++) {
              const MInt pIJK = pointIndex(i, j, k);
              const MFloat xInit = m_grid->m_initCoordinates[0][pIJK];
              const MFloat yInit = m_grid->m_initCoordinates[1][pIJK];
              const MFloat zInit = m_grid->m_initCoordinates[2][pIJK];


              MFloat yTransitionFactor = F1;
              const MFloat normalDist = sqrt(POW2(xInit - xWall) + POW2(yInit - yWall));
              if(normalDist <= m_waveYBeginTransition) {
                yTransitionFactor = F1;
              } else if(normalDist > m_waveYBeginTransition && normalDist < m_waveYEndTransition) {
                yTransitionFactor = (1 + cos((normalDist - m_waveYBeginTransition) / yTransitionLength * pi)) * F1B2;
              } else {
                yTransitionFactor = F0;
              }

              const MFloat normalDisplacement =
                  inc * fadeInFactor
                  * (amp0 * yTransitionFactor * transitionFactor
                     * (cos((F2 * pi) / m_waveLength * (zInit - m_waveSpeed * t_offset))));
              const MFloat normalVel =
                  inc * fadeInFactor
                  * (amp0 * yTransitionFactor * transitionFactor * F2 * pi / m_waveLength * m_waveSpeed
                     * (sin((F2 * pi) / m_waveLength * (zInit - m_waveSpeed * t_offset))));

              m_grid->m_coordinates[0][pIJK] = xInit + normalDisplacement * wallNormalVec[0];
              m_grid->m_coordinates[1][pIJK] = yInit + normalDisplacement * wallNormalVec[1];

              m_grid->m_velocity[0][pIJK] = normalVel * wallNormalVec[0];
              m_grid->m_velocity[1][pIJK] = normalVel * wallNormalVec[1];
            }
          }
        }
      }
      break;
    }
    case 14: {
      // oscillating cylinder
      MFloat t_offset = t - m_movingGridTimeOffset;
      if(zeroPos) {
        t_offset = F0;
      }

      MFloat fadeInFactor = 0;
      const MFloat timeRelaxation = 1.0;

      if(t_offset < timeRelaxation) {
        fadeInFactor = (1.0 - cos(t_offset / timeRelaxation * pi)) * F1B2;
      } else {
        fadeInFactor = 1.0;
      }

      if(zeroPos) {
        fadeInFactor = F1;
      }

      for(MInt i = 0; i < m_nPoints[2]; i++) {
        for(MInt k = 0; k < m_nPoints[0]; k++) {
          for(MInt j = 0; j < m_nPoints[1] - 1; j++) {
            const MInt pIJK = pointIndex(i, j, k);
            const MFloat x = m_grid->m_initCoordinates[0][pIJK];
            const MFloat y = m_grid->m_initCoordinates[1][pIJK];
            // const MFloat z = m_grid->m_initCoordinates[2][pIJK];

            const MFloat r = sqrt(POW2(x) + POW2(y));

            MFloat spaceTransition = F1;

            if(r < 1.0) {
              spaceTransition = F0;
            } else if(r >= 1.0 && r <= 41.0) {
              spaceTransition = fabs(r - 1.0) / 40.0;
            } else {
              spaceTransition = F1;
            }

            m_grid->m_coordinates[1][pIJK] =
                y * spaceTransition
                + (1.0 - spaceTransition) * (y + fadeInFactor * m_oscAmplitude * sin(2 * PI * m_oscFreq * t_offset));
            m_grid->m_velocity[1][pIJK] =
                (1.0 - spaceTransition)
                * (fadeInFactor * m_oscAmplitude * 2 * PI * m_oscFreq * cos(2 * PI * m_oscFreq * t_offset));
            m_grid->m_acceleration[1][pIJK] =
                (1.0 - spaceTransition)
                * (-fadeInFactor * m_oscAmplitude * POW2(2 * PI * m_oscFreq) * sin(2 * PI * m_oscFreq * t_offset));
          }
        }
      }

      break;
    }
    case 15: {
      // streamwise traveling wave on airfoil
      MFloat t_offset = t - m_movingGridTimeOffset;
      if(zeroPos) {
        t_offset = F0;
      }
      const MFloat transitionLength = m_waveEndTransition - m_waveBeginTransition;
      const MFloat transitionOutLength = m_waveOutEndTransition - m_waveOutBeginTransition;
      const MFloat pressureTransitionLength = m_wavePressureEndTransition - m_wavePressureBeginTransition;
      const MFloat pressureTransitionOutLength = m_wavePressureOutEndTransition - m_wavePressureOutBeginTransition;
      const MFloat yTransitionLength = m_waveYEndTransition - m_waveYBeginTransition;

      MFloat fadeInFactor = 0;
      const MFloat timeRelaxation = 1.0; // 80.0

      if(t_offset < timeRelaxation) {
        fadeInFactor = (1.0 - cos(t_offset / timeRelaxation * pi)) * F1B2;
      } else {
        fadeInFactor = 1.0;
      }

      if(zeroPos) {
        fadeInFactor = F1;
      }

      const MInt myBlockId = m_grid->getMyBlockId();

      if(myBlockId == 0) {
        for(MInt i = 0; i < m_nPoints[2]; i++) {
          MFloat transitionFactor = F0;
          const MInt offset = m_nOffsetCells[2];
          const MInt globalPointId = offset + i;
          const MFloat xWall = m_airfoilCoords[0 + globalPointId];
          const MFloat yWall = m_airfoilCoords[m_airfoilNoWallPoints + globalPointId];

          const MFloat wallNormalVec[3] = {m_airfoilNormalVec[0 + globalPointId],
                                           m_airfoilNormalVec[m_airfoilNoWallPoints + globalPointId],
                                           m_airfoilNormalVec[2 * m_airfoilNoWallPoints + globalPointId]};

          // upper surface
          MFloat amp0 = m_waveAmplitudeSuction;
          if(xWall <= m_waveBeginTransition) {
            transitionFactor = F0;
          } else if(xWall > m_waveBeginTransition && xWall < m_waveEndTransition) {
            transitionFactor = (1 - cos((xWall - m_waveBeginTransition) / transitionLength * pi)) * F1B2;
          } else if(m_waveEndTransition <= xWall && xWall <= m_waveOutBeginTransition) {
            transitionFactor = F1;
          } else if(xWall > m_waveOutBeginTransition && xWall < m_waveOutEndTransition) {
            transitionFactor = (1 + cos((xWall - m_waveOutBeginTransition) / transitionOutLength * pi)) * F1B2;
          } else {
            transitionFactor = F0;
          }

          // lower surface
          if(wallNormalVec[1] < 0.0) {
            amp0 = m_waveAmplitudePressure;
            if(xWall <= m_wavePressureBeginTransition) {
              transitionFactor = F0;
            } else if(xWall > m_wavePressureBeginTransition && xWall < m_wavePressureEndTransition) {
              transitionFactor =
                  (1 - cos((xWall - m_wavePressureBeginTransition) / pressureTransitionLength * pi)) * F1B2;
            } else if(m_wavePressureEndTransition <= xWall && xWall <= m_wavePressureOutBeginTransition) {
              transitionFactor = F1;
            } else if(xWall > m_wavePressureOutBeginTransition && xWall < m_wavePressureOutEndTransition) {
              transitionFactor =
                  (1 + cos((xWall - m_wavePressureOutBeginTransition) / pressureTransitionOutLength * pi)) * F1B2;
            } else {
              transitionFactor = F0;
            }
          }

          for(MInt k = 0; k < m_nPoints[0]; k++) {
            for(MInt j = 0; j < m_nPoints[1] - 1; j++) {
              const MInt pIJK = pointIndex(i, j, k);
              const MFloat xInit = m_grid->m_initCoordinates[0][pIJK];
              const MFloat yInit = m_grid->m_initCoordinates[1][pIJK];

              MFloat yTransitionFactor = F1;
              const MFloat normalDist = sqrt(POW2(xInit - xWall) + POW2(yInit - yWall));
              if(normalDist <= m_waveYBeginTransition) {
                yTransitionFactor = F1;
              } else if(normalDist > m_waveYBeginTransition && normalDist < m_waveYEndTransition) {
                yTransitionFactor = (1 + cos((normalDist - m_waveYBeginTransition) / yTransitionLength * pi)) * F1B2;
              } else {
                yTransitionFactor = F0;
              }

              const MFloat normalDisplacement =
                  fadeInFactor
                  * (amp0 * yTransitionFactor * transitionFactor
                     * (cos((F2 * pi) / m_waveLength * (xInit - m_waveSpeed * t_offset))));
              const MFloat normalVel =
                  fadeInFactor
                  * (amp0 * yTransitionFactor * transitionFactor * F2 * pi / m_waveLength * m_waveSpeed
                     * (sin((F2 * pi) / m_waveLength * (xInit - m_waveSpeed * t_offset))));

              m_grid->m_coordinates[0][pIJK] = xInit + normalDisplacement * wallNormalVec[0];
              m_grid->m_coordinates[1][pIJK] = yInit + normalDisplacement * wallNormalVec[1];

              m_grid->m_velocity[0][pIJK] = normalVel * wallNormalVec[0];
              m_grid->m_velocity[1][pIJK] = normalVel * wallNormalVec[1];
            }
          }
        }
      }
      break;
    }
    default: {
      mTerm(1, AT_, "Grid Moving Method not implemented!");
    }
  }

  RECORD_TIMER_START(m_timers[Timers::MGExchange]);

  if(m_gridMovingMethod != 1) {
    if(m_mgExchangeCoordinates) {
      m_grid->extrapolateGhostPointCoordinates();
    }
    extrapolateGhostPointCoordinatesBC();
  }

  if(noDomains() > 1) {
    if(m_mgExchangeCoordinates) {
      m_grid->exchangePoints(m_sndComm, m_rcvComm, PARTITION_NORMAL);
    }
  }
  RECORD_TIMER_STOP(m_timers[Timers::MGExchange]);

  RECORD_TIMER_START(m_timers[Timers::MGCellCenterCoordinates]);
  m_grid->computeCellCenterCoordinates();
  RECORD_TIMER_STOP(m_timers[Timers::MGCellCenterCoordinates]);

  RECORD_TIMER_START(m_timers[Timers::MGMetrics]);
  m_grid->computeMetrics();
  RECORD_TIMER_STOP(m_timers[Timers::MGMetrics]);

  RECORD_TIMER_START(m_timers[Timers::MGJacobian]);
  m_grid->computeJacobian();

  RECORD_TIMER_STOP(m_timers[Timers::MGJacobian]);
}

void FvStructuredSolver3D::initMovingGrid() {
  TRACE();
  m_log << "Initializing moving grid methods..." << endl;

  // First approach: save whole mesh in m_grid->m_initCoordinates
  for(MInt k = 0; k < m_nPoints[0]; ++k) {
    for(MInt j = 0; j < m_nPoints[1]; ++j) {
      for(MInt i = 0; i < m_nPoints[2]; ++i) {
        const MInt pointId = pointIndex(i, j, k);
        for(MInt isd = xsd; isd < nDim; ++isd) {
          m_grid->m_initCoordinates[isd][pointId] = m_grid->m_coordinates[isd][pointId];
        }
      }
    }
  }

  // Second approach: save only parts of the mesh depending on moving grid case
  switch(m_gridMovingMethod) {
    case 1:
    case 2:
    case 3:
    case 4: {
      break;
    }
    case 9:
      // traveling wave defined by viscous units
      // used Smits formula to compute friction velocity
      {
        m_travelingWave = true;
        if(!m_restart) {
          m_waveTimeStepComputed = false;
        }
        m_waveSpeed = 0.0;
        m_waveLength = 0.0;
        m_waveAmplitude = 0.0;
        m_waveCellsPerWaveLength = 1;
        if(!m_restart) {
          m_waveNoStepsPerCell = 1;
        }

        // time needs to be constant for traveling wave
        m_constantTimeStep = true;

        /*! \property
          \page propertiesFVSTRCTRD
                \section waveLengthPlus
                <code>MInt FvStructuredSolver::m_waveLengthPlus </code>\n
                default = <code> 1.0 </code>\n \n
                Wavelength of the traveling wave in inner units.\n
                Possible values are:\n
                <ul>
                <li>Float > 0.0</li>
                </ul>
                Keywords: <i>WAVE, MOVING, STRUCTURED</i>
              */
        m_waveLengthPlus = 1.0;
        if(Context::propertyExists("waveLengthPlus", m_solverId)) {
          m_waveLengthPlus = Context::getSolverProperty<MFloat>("waveLengthPlus", m_solverId, AT_, &m_waveLengthPlus);
        } else {
          mTerm(1, AT_, "Property waveLengthPlus not specified in property file");
        }

        /*! \property
          \page propertiesFVSTRCTRD
                \section waveAmplitudePlus
                <code>MInt FvStructuredSolver::m_waveAmplitudePlus </code>\n
                default = <code> 1.0 </code>\n \n
                Amplitude of the traveling wave in inner units.\n
                Possible values are:\n
                <ul>
                <li>Float > 0.0</li>
                </ul>
                Keywords: <i>WAVE, MOVING, STRUCTURED</i>
              */
        m_waveAmplitudePlus = 0.0;
        if(Context::propertyExists("waveAmplitudePlus", m_solverId)) {
          m_waveAmplitudePlus =
              Context::getSolverProperty<MFloat>("waveAmplitudePlus", m_solverId, AT_, &m_waveAmplitudePlus);
        } else {
          mTerm(1, AT_, "Property waveAmplitudePlus not specified in property file");
        }

        /*! \property
          \page propertiesFVSTRCTRD
                \section waveTimePlus
                <code>MInt FvStructuredSolver::m_waveTimePlus </code>\n
                default = <code> 1.0 </code>\n \n
                Period time of the traveling wave in inner units.\n
                Possible values are:\n
                <ul>
                <li>Float > 0.0</li>
                </ul>
                Keywords: <i>WAVE, MOVING, STRUCTURED</i>
              */
        m_waveTimePlus = 0.0;
        if(Context::propertyExists("waveTimePlus", m_solverId)) {
          m_waveTimePlus = Context::getSolverProperty<MFloat>("waveTimePlus", m_solverId, AT_, &m_waveTimePlus);
        } else {
          mTerm(1, AT_, "Property waveTimePlus not specified in property file");
        }

        /*! \property
          \page propertiesFVSTRCTRD
                \section waveBeginTransition
                <code>MInt FvStructuredSolver::m_waveBeginTransition </code>\n
                default = <code> 1.0 </code>\n \n
                Start of the transition from flat to wave in x-dir.\n
                Possible values are:\n
                <ul>
                <li>Float</li>
                </ul>
                Keywords: <i>WAVE, MOVING, STRUCTURED</i>
              */
        m_waveBeginTransition = 0.0;
        m_waveBeginTransition =
            Context::getSolverProperty<MFloat>("waveBeginTransition", m_solverId, AT_, &m_waveBeginTransition);

        /*! \property
          \page propertiesFVSTRCTRD
                \section waveEndTransition
                <code>MInt FvStructuredSolver::m_waveEndTransition </code>\n
                default = <code> 1.0 </code>\n \n
                End of the transition from flat to wave in x-dir.\n
                Possible values are:\n
                <ul>
                <li>Float</li>
                </ul>
                Keywords: <i>WAVE, MOVING, STRUCTURED</i>
              */
        m_waveEndTransition = 0.0;
        m_waveEndTransition =
            Context::getSolverProperty<MFloat>("waveEndTransition", m_solverId, AT_, &m_waveEndTransition);

        /*! \property
          \page propertiesFVSTRCTRD
                \section waveOutBeginTransition
                <code>MInt FvStructuredSolver::m_waveOutBeginTransition </code>\n
                default = <code> 1.0 </code>\n \n
                Start of the transition from wave to flat in x-dir.\n
                Possible values are:\n
                <ul>
                <li>Float</li>
                </ul>
                Keywords: <i>WAVE, MOVING, STRUCTURED</i>
              */
        m_waveOutBeginTransition = 1000000.0;
        m_waveOutBeginTransition =
            Context::getSolverProperty<MFloat>("waveOutBeginTransition", m_solverId, AT_, &m_waveOutBeginTransition);

        /*! \property
          \page propertiesFVSTRCTRD
                \section waveOutEndTransition
                <code>MInt FvStructuredSolver::m_waveOutEndTransition </code>\n
                default = <code> 1.0 </code>\n \n
                End of the transition from wave to flat in x-dir.\n
                Possible values are:\n
                <ul>
                <li>Float</li>
                </ul>
                Keywords: <i>WAVE, MOVING, STRUCTURED</i>
              */
        m_waveOutEndTransition = 2000000.0;
        m_waveOutEndTransition =
            Context::getSolverProperty<MFloat>("waveOutEndTransition", m_solverId, AT_, &m_waveOutEndTransition);

        /*! \property
          \page propertiesFVSTRCTRD
                \section waveAngle
                <code>MInt FvStructuredSolver::m_waveAngle </code>\n
                default = <code> 1.0 </code>\n \n
                Angle of the wave.\n
                Possible values are:\n
                <ul>
                <li>Degrees of angle as a Float</li>
                </ul>
                Keywords: <i>WAVE, MOVING, STRUCTURED</i>
              */
        m_waveAngle = 0.0;
        if(Context::propertyExists("waveAngle", m_solverId)) {
          m_waveAngle = Context::getSolverProperty<MFloat>("waveAngle", m_solverId, AT_, &m_waveAngle);
          m_waveAngle *= PI / 180.0;
        }

        /*! \property
          \page propertiesFVSTRCTRD
                \section waveYBeginTransition
                <code>MInt FvStructuredSolver::m_waveYBeginTransition </code>\n
                default = <code> 1.0 </code>\n \n
                End of the transition from wave to flat in x-dir.\n
                Possible values are:\n
                <ul>
                <li>Float</li>
                </ul>
                Keywords: <i>WAVE, MOVING, STRUCTURED</i>
              */
        m_waveYBeginTransition = 2000000.0;
        if(Context::propertyExists("waveYBeginTransition", m_solverId)) {
          m_waveYBeginTransition =
              Context::getSolverProperty<MFloat>("waveYBeginTransition", m_solverId, AT_, &m_waveYBeginTransition);
        }


        /*! \property
          \page propertiesFVSTRCTRD
                \section waveYEndTransition
                <code>MInt FvStructuredSolver::m_waveYEndTransition </code>\n
                default = <code> 1.0 </code>\n \n
                End of the transition from wave to flat in x-dir.\n
                Possible values are:\n
                <ul>
                <li>Float</li>
                </ul>
                Keywords: <i>WAVE, MOVING, STRUCTURED</i>
              */
        m_waveYEndTransition = 2000000.0;
        if(Context::propertyExists("waveYEndTransition", m_solverId)) {
          m_waveYEndTransition =
              Context::getSolverProperty<MFloat>("waveYEndTransition", m_solverId, AT_, &m_waveYEndTransition);
        }

        /*! \property
          \page propertiesFVSTRCTRD
                \section waveTemporalTransition
                <code>MInt FvStructuredSolver::m_waveOutEndTransition </code>\n
                default = <code> 500.0 </code>\n \n
                Acoustic time for wave actuation to transiate from flat plate\n
                to fully extended.\n
                Possible values are:\n
                <ul>
                <li>Float</li>
                </ul>
                Keywords: <i>WAVE, MOVING, STRUCTURED</i>
              */
        m_waveTemporalTransition = 500.0;
        if(Context::propertyExists("waveTemporalTransition", m_solverId)) {
          m_waveTemporalTransition =
              Context::getSolverProperty<MFloat>("waveTemporalTransition", m_solverId, AT_, &m_waveTemporalTransition);
        }

        // compute Wave parameters
        MFloat deltaS = -1.0;
        MFloat cf = -1.0;
        if(m_Re < 5000) {
          // use Smits formula for smaller Reynolds numbers
          deltaS = pow(m_Re, -7.0 / 8.0) * sqrt(2.0 / 0.024);
          cf = 0.024 * pow(m_Re, -F1B4);
        } else {
          // Coles-Fernholz is better for larger Re
          deltaS = ((log(m_Re)) / 0.384 + 4.127) / m_Re;
          cf = 2.0 * pow((log(m_Re) / 0.384 + 4.127), -2.0);
        }
        const MFloat uTau = sqrt(POW2(PV->UInfinity) * CV->rhoInfinity * cf * F1B2);

        m_waveLength = m_waveLengthPlus * deltaS;
        m_waveAmplitude = m_waveAmplitudePlus * deltaS;
        m_waveSpeedPlus = m_waveLengthPlus / m_waveTimePlus;
        m_waveSpeed = m_waveSpeedPlus * uTau;

        // assume equidistant grid in z-direction
        const MFloat deltaZ = abs(m_grid->m_coordinates[2][0] - m_grid->m_coordinates[2][m_nPoints[2] * m_nPoints[1]]);
        m_waveCellsPerWaveLength = round(m_waveLength / deltaZ);


        const MFloat wavePeriod = m_waveLength / m_waveSpeed;
        const MFloat speedAmplitude = 2 * PI * m_waveAmplitude / wavePeriod;

        m_log << "/////////////////// TRAVELING WAVE /////////////////////////////" << endl;
        m_log << "Re: " << m_Re << " c_f: " << cf << " u_tau: " << uTau << " deltaZ: " << deltaZ
              << " viscousUnit delta_v: " << deltaS << endl;
        m_log << "wavelengthPlus: " << m_waveLengthPlus << " AmplitudePlus: " << m_waveAmplitudePlus
              << " SpeedPlus: " << m_waveSpeedPlus << endl;
        m_log << "Wavelength: " << m_waveLength << " Amplitude: " << m_waveAmplitude << " Period: " << wavePeriod
              << " Speed: " << m_waveSpeed << endl;
        m_log << "Max up/down speed: " << speedAmplitude << endl;
        m_log << "Max up/down speed: " << m_waveSpeed * m_waveAmplitude << endl;
        m_log << "////////////////////////////////////////////////////////////////" << endl;
        fixTimeStepTravelingWave();
        break;
      }
    case 10:
      // travelling wave defined by non-plus units
      {
        m_travelingWave = true;
        if(!m_restart) {
          m_waveTimeStepComputed = false;
        }
        m_waveSpeed = 0.0;
        m_waveLength = 0.0;
        m_waveAmplitude = 0.0;
        m_waveCellsPerWaveLength = 1;
        if(!m_restart) {
          m_waveNoStepsPerCell = 1;
        }

        // time needs to be constant for traveling wave
        m_constantTimeStep = true;

        /*! \property
          \page propertiesFVSTRCTRD
                \section waveLength
                <code>MInt FvStructuredSolver::m_waveLength </code>\n
                default = <code> 1.0 </code>\n \n
                Wavelength of the traveling wave.\n
                Possible values are:\n
                <ul>
                <li>Float > 0.0</li>
                </ul>
                Keywords: <i>WAVE, MOVING, STRUCTURED</i>
              */
        m_waveLength = 0.0;
        if(Context::propertyExists("waveLength", m_solverId)) {
          m_waveLength = Context::getSolverProperty<MFloat>("waveLength", m_solverId, AT_, &m_waveLengthPlus);
        } else {
          mTerm(1, AT_, "Property waveLength not specified in property file");
        }

        /*! \property
          \page propertiesFVSTRCTRD
                \section waveAmplitude
                <code>MInt FvStructuredSolver::m_waveAmplitude </code>\n
                default = <code> 1.0 </code>\n \n
                Amplitude of the traveling wave.\n
                Possible values are:\n
                <ul>
                <li>Float > 0.0</li>
                </ul>
                Keywords: <i>WAVE, MOVING, STRUCTURED</i>
              */
        m_waveAmplitude = 0.0;
        if(Context::propertyExists("waveAmplitude", m_solverId)) {
          m_waveAmplitude = Context::getSolverProperty<MFloat>("waveAmplitude", m_solverId, AT_, &m_waveAmplitudePlus);
        } else {
          mTerm(1, AT_, "Property waveAmplitude not specified in property file");
        }

        /*! \property
          \page propertiesFVSTRCTRD
                \section waveTime
                <code>MInt FvStructuredSolver::m_waveTime </code>\n
                default = <code> 1.0 </code>\n \n
                Period time of the traveling wave.\n
                Possible values are:\n
                <ul>
                <li>Float > 0.0</li>
                </ul>
                Keywords: <i>WAVE, MOVING, STRUCTURED</i>
              */
        m_waveTime = 0.0;
        if(Context::propertyExists("waveTime", m_solverId)) {
          m_waveTime = Context::getSolverProperty<MFloat>("waveTime", m_solverId, AT_, &m_waveTimePlus);
        } else {
          mTerm(1, AT_, "Property waveTime not specified in property file");
        }


        /*! \property
          \page propertiesFVSTRCTRD
                \section waveBeginTransition
                <code>MInt FvStructuredSolver::m_waveBeginTransition </code>\n
                default = <code> 1.0 </code>\n \n
                Start of the transition from flat to wave in x-dir.\n
                Possible values are:\n
                <ul>
                <li>Float</li>
                </ul>
                Keywords: <i>WAVE, MOVING, STRUCTURED</i>
              */
        m_waveBeginTransition = 0.0;
        m_waveBeginTransition =
            Context::getSolverProperty<MFloat>("waveBeginTransition", m_solverId, AT_, &m_waveBeginTransition);

        /*! \property
          \page propertiesFVSTRCTRD
                \section waveEndTransition
                <code>MInt FvStructuredSolver::m_waveEndTransition </code>\n
                default = <code> 1.0 </code>\n \n
                End of the transition from flat to wave in x-dir.\n
                Possible values are:\n
                <ul>
                <li>Float</li>
                </ul>
                Keywords: <i>WAVE, MOVING, STRUCTURED</i>
              */
        m_waveEndTransition = 0.0;
        m_waveEndTransition =
            Context::getSolverProperty<MFloat>("waveEndTransition", m_solverId, AT_, &m_waveEndTransition);

        /*! \property
          \page propertiesFVSTRCTRD
                \section waveOutBeginTransition
                <code>MInt FvStructuredSolver::m_waveOutBeginTransition </code>\n
                default = <code> 1.0 </code>\n \n
                Start of the transition from wave to flat in x-dir.\n
                Possible values are:\n
                <ul>
                <li>Float</li>
                </ul>
                Keywords: <i>WAVE, MOVING, STRUCTURED</i>
              */
        m_waveOutBeginTransition = 1000000.0;
        m_waveOutBeginTransition =
            Context::getSolverProperty<MFloat>("waveOutBeginTransition", m_solverId, AT_, &m_waveOutBeginTransition);

        /*! \property
          \page propertiesFVSTRCTRD
                \section waveOutEndTransition
                <code>MInt FvStructuredSolver::m_waveOutEndTransition </code>\n
                default = <code> 1.0 </code>\n \n
                End of the transition from wave to flat in x-dir.\n
                Possible values are:\n
                <ul>
                <li>Float</li>
                </ul>
                Keywords: <i>WAVE, MOVING, STRUCTURED</i>
              */
        m_waveOutEndTransition = 2000000.0;
        m_waveOutEndTransition =
            Context::getSolverProperty<MFloat>("waveOutEndTransition", m_solverId, AT_, &m_waveOutEndTransition);

        /*! \property
          \page propertiesFVSTRCTRD
                \section waveYBeginTransition
                <code>MInt FvStructuredSolver::m_waveYBeginTransition </code>\n
                default = <code> 1.0 </code>\n \n
                End of the transition from wave to flat in x-dir.\n
                Possible values are:\n
                <ul>
                <li>Float</li>
                </ul>
                Keywords: <i>WAVE, MOVING, STRUCTURED</i>
              */
        m_waveYBeginTransition = 2000000.0;
        if(Context::propertyExists("waveYBeginTransition", m_solverId)) {
          m_waveYBeginTransition =
              Context::getSolverProperty<MFloat>("waveYBeginTransition", m_solverId, AT_, &m_waveYBeginTransition);
        }


        /*! \property
          \page propertiesFVSTRCTRD
                \section waveYEndTransition
                <code>MInt FvStructuredSolver::m_waveYEndTransition </code>\n
                default = <code> 1.0 </code>\n \n
                End of the transition from wave to flat in x-dir.\n
                Possible values are:\n
                <ul>
                <li>Float</li>
                </ul>
                Keywords: <i>WAVE, MOVING, STRUCTURED</i>
              */
        m_waveYEndTransition = 2000000.0;
        if(Context::propertyExists("waveYEndTransition", m_solverId)) {
          m_waveYEndTransition =
              Context::getSolverProperty<MFloat>("waveYEndTransition", m_solverId, AT_, &m_waveYEndTransition);
        }

        /*! \property
          \page propertiesFVSTRCTRD
                \section waveAngle
                <code>MInt FvStructuredSolver::m_waveAngle </code>\n
                default = <code> 1.0 </code>\n \n
                Angle of the wave.\n
                Possible values are:\n
                <ul>
                <li>Degrees of angle as a Float</li>
                </ul>
                Keywords: <i>WAVE, MOVING, STRUCTURED</i>
              */
        m_waveAngle = 0.0;
        if(Context::propertyExists("waveAngle", m_solverId)) {
          m_waveAngle = Context::getSolverProperty<MFloat>("waveAngle", m_solverId, AT_, &m_waveAngle);
          m_waveAngle *= PI / 180.0;
        }

        /*! \property
          \page propertiesFVSTRCTRD
                \section waveTemporalTransition
                <code>MInt FvStructuredSolver::m_waveOutEndTransition </code>\n
                default = <code> 500.0 </code>\n \n
                Acoustic time for wave actuation to transiate from flat plate\n
                to fully extended.\n
                Possible values are:\n
                <ul>
                <li>Float</li>
                </ul>
                Keywords: <i>WAVE, MOVING, STRUCTURED</i>
              */
        m_waveTemporalTransition = 500.0;
        if(Context::propertyExists("waveTemporalTransition", m_solverId)) {
          m_waveTemporalTransition =
              Context::getSolverProperty<MFloat>("waveTemporalTransition", m_solverId, AT_, &m_waveTemporalTransition);
        }

        m_waveSpeed = m_waveLength / m_waveTime;
        const MFloat deltaZ = abs(m_grid->m_coordinates[2][0] - m_grid->m_coordinates[2][m_nPoints[2] * m_nPoints[1]]);
        m_waveCellsPerWaveLength = round(m_waveLength / deltaZ);

        const MFloat speedAmplitude = 2 * PI * m_waveAmplitude / m_waveTime;

        m_log << "/////////////////// TRAVELING WAVE /////////////////////////////" << endl;
        m_log << "Wavelength: " << m_waveLength << " Amplitude: " << m_waveAmplitude << " Period: " << m_waveTime
              << " Speed: " << m_waveSpeed << endl;
        m_log << "Max up/down speed: " << speedAmplitude << endl;
        m_log << "Max up/down speed: " << m_waveSpeed * m_waveAmplitude << endl;
        m_log << "////////////////////////////////////////////////////////////////" << endl;

        fixTimeStepTravelingWave();
        break;
      }
    case 11:
      // traveling wave channel (Tomiyama & Fukagata 2013)
      {
        m_travelingWave = true;
        m_waveSpeed = 0.0;
        m_waveLength = 0.0;
        m_waveAmplitude = 0.0;
        m_waveLengthPlus = 0.0;
        m_waveAmplitudePlus = 0.0;
        m_waveTimePlus = 0.0;

        // time needs to be constant for traveling wave
        m_constantTimeStep = true;
        m_waveLengthPlus = Context::getSolverProperty<MFloat>("waveLengthPlus", m_solverId, AT_, &m_waveLengthPlus);
        m_waveAmplitudePlus =
            Context::getSolverProperty<MFloat>("waveAmplitudePlus", m_solverId, AT_, &m_waveAmplitudePlus);
        m_waveTimePlus = Context::getSolverProperty<MFloat>("waveTimePlus", m_solverId, AT_, &m_waveTimePlus);

        // compute Wave parameters
        // const MFloat cf = 0.008185; //0.024*pow(m_Re, -F1B4);
        // const MFloat uTau = sqrt(POW2(PV->UInfinity)*CV->rhoInfinity*cf*F1B2);
        const MFloat uTau = m_ReTau * m_Ma * sqrt(PV->TInfinity) / m_Re;
        const MFloat cf = 2 * POW2(uTau / PV->UInfinity);

        const MFloat mu8 = SUTHERLANDLAW(PV->TInfinity);
        m_waveLength = m_waveLengthPlus / (sqrt(cf / 2.0) * PV->UInfinity * m_Re0 * CV->rhoInfinity / mu8);
        m_waveAmplitude = m_waveAmplitudePlus / (sqrt(cf / 2.0) * PV->UInfinity * m_Re0 * CV->rhoInfinity / mu8);
        m_waveSpeed = (m_waveLengthPlus / m_waveTimePlus) * uTau; /// PV->UInfinity;
        const MFloat deltaZ = abs(m_grid->m_coordinates[2][0] - m_grid->m_coordinates[2][m_nPoints[2] * m_nPoints[1]]);
        m_waveCellsPerWaveLength = round(m_waveLength / deltaZ);

        m_log << "WaveLength: " << m_waveLength << " waveAmplitude: " << m_waveAmplitude
              << " waveSpeed: " << m_waveSpeed << endl;
        fixTimeStepTravelingWave();
        break;
      }
    case 12:
      // streamwise travelling wave defined by non-plus units
      {
        m_streamwiseTravelingWave = true;
        if(!m_restart) {
          m_waveTimeStepComputed = false;
        }
        m_waveSpeed = 0.0;
        m_waveLength = 0.0;
        m_waveAmplitude = 0.0;
        m_waveCellsPerWaveLength = 1;
        if(!m_restart) {
          m_waveNoStepsPerCell = 1;
        }
        // time needs to be constant for traveling wave
        m_constantTimeStep = true;

        /*! \property
          \page propertiesFVSTRCTRD
                \section waveLength
                <code>MInt FvStructuredSolver::m_waveLength </code>\n
                default = <code> 1.0 </code>\n \n
                Wavelength of the traveling wave.\n
                Possible values are:\n
                <ul>
                <li>Float > 0.0</li>
                </ul>
                Keywords: <i>WAVE, MOVING, STRUCTURED</i>
              */
        m_waveLength = 0.0;
        if(Context::propertyExists("waveLength", m_solverId)) {
          m_waveLength = Context::getSolverProperty<MFloat>("waveLength", m_solverId, AT_, &m_waveLengthPlus);
        } else {
          mTerm(1, AT_, "Property waveLength not specified in property file");
        }

        /*! \property
          \page propertiesFVSTRCTRD
                \section waveAmplitudeSuction
                <code>MInt FvStructuredSolver::m_waveAmplitudeSuction </code>\n
                default = <code> 1.0 </code>\n \n
                Amplitude of the traveling wave on the airfoil suction side.\n
                Possible values are:\n
                <ul>
                <li>Float > 0.0</li>
                </ul>
                Keywords: <i>WAVE, MOVING, STRCTRD</i>
              */
        m_waveAmplitudeSuction = 0.0;
        m_waveAmplitudeSuction =
            Context::getSolverProperty<MFloat>("waveAmplitudeSuction", m_blockId, AT_, &m_waveAmplitudeSuction);

        /*! \property
          \page propertiesFVSTRCTRD
                \section waveAmplitudePressure
                <code>MInt FvStructuredSolver::m_waveAmplitudePressure </code>\n
      >>>>>>> origin/improvedPartitioning
                default = <code> 1.0 </code>\n \n
                Amplitude of the traveling wave.\n
                Possible values are:\n
                <ul>
                <li>Float > 0.0</li>
                </ul>
                Keywords: <i>WAVE, MOVING, STRUCTURED</i>
              */
        m_waveAmplitude = 0.0;
        if(Context::propertyExists("waveAmplitude", m_solverId)) {
          m_waveAmplitude = Context::getSolverProperty<MFloat>("waveAmplitude", m_solverId, AT_, &m_waveAmplitudePlus);
        } else {
          mTerm(1, AT_, "Property waveAmplitude not specified in property file");
        }

        /*! \property
          \page propertiesFVSTRCTRD
                \section waveTime
                <code>MInt FvStructuredSolver::m_waveTime </code>\n
                default = <code> 1.0 </code>\n \n
                Period time of the traveling wave.\n
                Possible values are:\n
                <ul>
                <li>Float > 0.0</li>
                </ul>
                Keywords: <i>WAVE, MOVING, STRUCTURED</i>
              */
        m_waveTime = 0.0;
        if(Context::propertyExists("waveTime", m_solverId)) {
          m_waveTime = Context::getSolverProperty<MFloat>("waveTime", m_solverId, AT_, &m_waveTimePlus);
        } else {
          mTerm(1, AT_, "Property waveTime not specified in property file");
        }


        /*! \property
          \page propertiesFVSTRCTRD
                \section waveBeginTransition
                <code>MInt FvStructuredSolver::m_waveBeginTransition </code>\n
                default = <code> 1.0 </code>\n \n
                Start of the transition from flat to wave in x-dir.\n
                Possible values are:\n
                <ul>
                <li>Float</li>
                </ul>
                Keywords: <i>WAVE, MOVING, STRUCTURED</i>
              */
        m_waveBeginTransition = 0.0;
        m_waveBeginTransition =
            Context::getSolverProperty<MFloat>("waveBeginTransition", m_solverId, AT_, &m_waveBeginTransition);

        /*! \property
          \page propertiesFVSTRCTRD
                \section waveEndTransition
                <code>MInt FvStructuredSolver::m_waveEndTransition </code>\n
                default = <code> 1.0 </code>\n \n
                End of the transition from flat to wave in x-dir.\n
                Possible values are:\n
                <ul>
                <li>Float</li>
                </ul>
                Keywords: <i>WAVE, MOVING, STRUCTURED</i>
              */
        m_waveEndTransition = 0.0;
        m_waveEndTransition =
            Context::getSolverProperty<MFloat>("waveEndTransition", m_solverId, AT_, &m_waveEndTransition);

        /*! \property
          \page propertiesFVSTRCTRD
                \section waveOutBeginTransition
                <code>MInt FvStructuredSolver::m_waveOutBeginTransition </code>\n
                default = <code> 1.0 </code>\n \n
                Start of the transition from wave to flat in x-dir.\n
                Possible values are:\n
                <ul>
                <li>Float</li>
                </ul>
                Keywords: <i>WAVE, MOVING, STRUCTURED</i>
              */
        m_waveOutBeginTransition = 1000000.0;
        m_waveOutBeginTransition =
            Context::getSolverProperty<MFloat>("waveOutBeginTransition", m_solverId, AT_, &m_waveOutBeginTransition);

        /*! \property
          \page propertiesFVSTRCTRD
                \section waveOutEndTransition
                <code>MInt FvStructuredSolver::m_waveOutEndTransition </code>\n
                default = <code> 1.0 </code>\n \n
                End of the transition from wave to flat in x-dir.\n
                Possible values are:\n
                <ul>
                <li>Float</li>
                </ul>
                Keywords: <i>WAVE, MOVING, STRUCTURED</i>
              */
        m_waveOutEndTransition = 2000000.0;
        m_waveOutEndTransition =
            Context::getSolverProperty<MFloat>("waveOutEndTransition", m_solverId, AT_, &m_waveOutEndTransition);

        /*! \property
          \page propertiesFVSTRCTRD
                \section waveYBeginTransition
                <code>MInt FvStructuredSolver::m_waveYBeginTransition </code>\n
                default = <code> 1.0 </code>\n \n
                End of the transition from wave to flat in x-dir.\n
                Possible values are:\n
                <ul>
                <li>Float</li>
                </ul>
                Keywords: <i>WAVE, MOVING, STRUCTURED</i>
              */
        m_waveYBeginTransition = 2000000.0;
        if(Context::propertyExists("waveYBeginTransition", m_solverId)) {
          m_waveYBeginTransition =
              Context::getSolverProperty<MFloat>("waveYBeginTransition", m_solverId, AT_, &m_waveYBeginTransition);
        }


        /*! \property
          \page propertiesFVSTRCTRD
                \section waveYEndTransition
                <code>MInt FvStructuredSolver::m_waveYEndTransition </code>\n
                default = <code> 1.0 </code>\n \n
                End of the transition from wave to flat in x-dir.\n
                Possible values are:\n
                <ul>
                <li>Float</li>
                </ul>
                Keywords: <i>WAVE, MOVING, STRUCTURED</i>
              */
        m_waveYEndTransition = 2000000.0;
        if(Context::propertyExists("waveYEndTransition", m_solverId)) {
          m_waveYEndTransition =
              Context::getSolverProperty<MFloat>("waveYEndTransition", m_solverId, AT_, &m_waveYEndTransition);
        }

        /*! \page propertyPage1
          \section waveTemporalTransition
          <code>MInt FvStructuredSolver::m_waveOutEndTransition </code>\n
          default = <code> 500.0 </code>\n \n
          Acoustic time for wave actuation to transiate from flat plate\n
          to fully extended.\n
          Possible values are:\n
          <ul>
          <li>Float</li>
          </ul>
          Keywords: <i>WAVE, MOVING, STRUCTURED</i>
        */
        m_waveTemporalTransition = 500.0;
        if(Context::propertyExists("waveTemporalTransition", m_solverId)) {
          m_waveTemporalTransition =
              Context::getSolverProperty<MFloat>("waveTemporalTransition", m_solverId, AT_, &m_waveTemporalTransition);
        }

        m_waveSpeed = m_waveLength / m_waveTime;
        const MFloat speedAmplitude = 2 * PI * m_waveAmplitude / m_waveTime;

        m_log << "/////////////////// TRAVELING WAVE /////////////////////////////" << endl;
        m_log << "Wavelength: " << m_waveLength << " Amplitude: " << m_waveAmplitude << " Period: " << m_waveTime
              << " Speed: " << m_waveSpeed << endl;
        m_log << "Speed amplitude: " << speedAmplitude << endl;
        m_log << "Phase speed: " << m_waveSpeed << endl;
        m_log << "////////////////////////////////////////////////////////////////" << endl;
        fixTimeStepTravelingWave();
        break;
      }
    case 13:
      // traveling wave on airfoil
      {
        m_travelingWave = true;
        m_waveSpeed = 0.0;
        m_waveCellsPerWaveLength = 1;
        if(!m_restart) {
          m_waveTimeStepComputed = false;
        }
        // time needs to be constant for traveling wave
        m_constantTimeStep = true;
        if(!m_restart) {
          m_waveNoStepsPerCell = 1;
        }

        /*! \property
          \page propertiesFVSTRCTRD
                \section waveLength
                <code>MInt FvStructuredSolver::m_waveLength </code>\n
                default = <code> 1.0 </code>\n \n
                Wavelength of the traveling wave.\n
                Possible values are:\n
                <ul>
                <li>Float > 0.0</li>
                </ul>
                Keywords: <i>WAVE, MOVING, STRUCTURED</i>
              */
        m_waveLength = 0.0;
        m_waveLength = Context::getSolverProperty<MFloat>("waveLength", m_solverId, AT_, &m_waveLength);

        /*! \property
          \page propertiesFVSTRCTRD
                \section waveAmplitudeSuction
                <code>MInt FvStructuredSolver::m_waveAmplitudeSuction </code>\n
                default = <code> 1.0 </code>\n \n
                Amplitude of the traveling wave on the airfoil suction side.\n
                Possible values are:\n
                <ul>
                <li>Float > 0.0</li>
                </ul>
                Keywords: <i>WAVE, MOVING, STRCTRD</i>
              */
        m_waveAmplitudeSuction = 0.0;
        m_waveAmplitudeSuction =
            Context::getSolverProperty<MFloat>("waveAmplitudeSuction", m_blockId, AT_, &m_waveAmplitudeSuction);

        /*! \property
          \page propertiesFVSTRCTRD
                \section waveAmplitudePressure
                <code>MInt FvStructuredSolver::m_waveAmplitudePressure </code>\n
                default = <code> 1.0 </code>\n \n
                Amplitude of the traveling wave on the airfoil suction side.\n
                Possible values are:\n
                <ul>
                <li>Float > 0.0</li>
                </ul>
                Keywords: <i>WAVE, MOVING, STRCTRD</i>
              */
        m_waveAmplitudePressure = 0.0;
        m_waveAmplitudePressure =
            Context::getSolverProperty<MFloat>("waveAmplitudePressure", m_blockId, AT_, &m_waveAmplitudePressure);

        m_waveGradientSuction = 0.5;
        m_waveGradientSuction =
            Context::getSolverProperty<MFloat>("waveGradientSuction", m_blockId, AT_, &m_waveGradientSuction);

        m_waveGradientPressure = 0.0;
        m_waveGradientPressure =
            Context::getSolverProperty<MFloat>("waveGradientPressure", m_blockId, AT_, &m_waveGradientPressure);

        /*! \property
          \page propertiesFVSTRCTRD
                \section waveTime
                <code>MInt FvStructuredSolver::m_waveTime </code>\n
                default = <code> 1.0 </code>\n \n
                Period time of the traveling wave.\n
                Possible values are:\n
                <ul>
                <li>Float > 0.0</li>
                </ul>
                Keywords: <i>WAVE, MOVING, STRUCTURED</i>
              */
        m_waveTime = 0.0;
        m_waveTime = Context::getSolverProperty<MFloat>("waveTime", m_solverId, AT_, &m_waveTimePlus);

        /*! \property
          \page propertiesFVSTRCTRD
                \section waveBeginTransition
                <code>MInt FvStructuredSolver::m_waveBeginTransition </code>\n
                default = <code> 1.0 </code>\n \n
                Start of the transition from flat to wave in x-dir.\n
                Possible values are:\n
                <ul>
                <li>Float</li>
                </ul>
                Keywords: <i>WAVE, MOVING, STRUCTURED</i>
              */
        m_waveBeginTransition = 0.0;
        m_waveBeginTransition =
            Context::getSolverProperty<MFloat>("waveBeginTransition", m_solverId, AT_, &m_waveBeginTransition);

        /*! \property
          \page propertiesFVSTRCTRD
                \section waveEndTransition
                <code>MInt FvStructuredSolver::m_waveEndTransition </code>\n
                default = <code> 1.0 </code>\n \n
                End of the transition from flat to wave in x-dir.\n
                Possible values are:\n
                <ul>
                <li>Float</li>
                </ul>
                Keywords: <i>WAVE, MOVING, STRUCTURED</i>
              */
        m_waveEndTransition = 0.0;
        m_waveEndTransition =
            Context::getSolverProperty<MFloat>("waveEndTransition", m_solverId, AT_, &m_waveEndTransition);

        /*! \property
          \page propertiesFVSTRCTRD
                \section waveOutBeginTransition
                <code>MInt FvStructuredSolver::m_waveOutBeginTransition </code>\n
                default = <code> 1.0 </code>\n \n
                Start of the transition from wave to flat in x-dir.\n
                Possible values are:\n
                <ul>
                <li>Float</li>
                </ul>
                Keywords: <i>WAVE, MOVING, STRUCTURED</i>
              */
        m_waveOutBeginTransition = 1000000.0;
        m_waveOutBeginTransition =
            Context::getSolverProperty<MFloat>("waveOutBeginTransition", m_solverId, AT_, &m_waveOutBeginTransition);

        /*! \property
          \page propertiesFVSTRCTRD
                \section waveOutEndTransition
                <code>MInt FvStructuredSolver::m_waveOutEndTransition </code>\n
                default = <code> 1.0 </code>\n \n
                End of the transition from wave to flat in x-dir.\n
                Possible values are:\n
                <ul>
                <li>Float</li>
                </ul>
                Keywords: <i>WAVE, MOVING, STRUCTURED</i>
              */
        m_waveOutEndTransition = 2000000.0;
        m_waveOutEndTransition =
            Context::getSolverProperty<MFloat>("waveOutEndTransition", m_solverId, AT_, &m_waveOutEndTransition);

        /*! \property
          \page propertiesFVSTRCTRD
                \section waveYBeginTransition
                <code>MInt FvStructuredSolver::m_waveYBeginTransition </code>\n
                default = <code> 1.0 </code>\n \n
                End of the transition from wave to flat in x-dir.\n
                Possible values are:\n
                <ul>
                <li>Float</li>
                </ul>
                Keywords: <i>WAVE, MOVING, STRUCTURED</i>
              */
        m_waveYBeginTransition = 2000000.0;
        if(Context::propertyExists("waveYBeginTransition", m_solverId)) {
          m_waveYBeginTransition =
              Context::getSolverProperty<MFloat>("waveYBeginTransition", m_solverId, AT_, &m_waveYBeginTransition);
        }


        /*! \property
          \page propertiesFVSTRCTRD
                \section waveYEndTransition
                <code>MInt FvStructuredSolver::m_waveYEndTransition </code>\n
                default = <code> 1.0 </code>\n \n
                End of the transition from wave to flat in x-dir.\n
                Possible values are:\n
                <ul>
                <li>Float</li>
                </ul>
                Keywords: <i>WAVE, MOVING, STRUCTURED</i>
              */
        m_waveYEndTransition = 2000000.0;
        if(Context::propertyExists("waveYEndTransition", m_solverId)) {
          m_waveYEndTransition =
              Context::getSolverProperty<MFloat>("waveYEndTransition", m_solverId, AT_, &m_waveYEndTransition);
        }

        const MInt myBlockId = m_grid->getMyBlockId();
        const MInt blockId = 0; // use information from block 0
        m_airfoilNoWallPoints = m_grid->getBlockNoPoints(blockId, 2) + 2 * m_noGhostLayers;

        if(domainId() == 0) {
          cout << "NoWallPoints: " << m_airfoilNoWallPoints << endl;
        }

        MFloatScratchSpace localCoords(nDim, m_airfoilNoWallPoints, AT_, "localCoords");
        MFloatScratchSpace localNormalVec(nDim, m_airfoilNoWallPoints, AT_, "localNormalVec");
        localCoords.fill(-99999.0);
        localNormalVec.fill(-99999.0);
        mAlloc(m_airfoilCoords, nDim * m_airfoilNoWallPoints, "m_airfoilCoords", -99999.0, AT_);
        mAlloc(m_airfoilNormalVec, nDim * m_airfoilNoWallPoints, "m_airfoilNormalVec", -99999.0, AT_);

        if(myBlockId == blockId) {
          if(m_nOffsetCells[1] == 0 && m_nOffsetCells[0] == 0) {
            // only collect variables from points at the wall
            for(MInt i = 0; i < m_nPoints[2]; i++) {
              const MInt offset = m_nOffsetCells[2];
              const MInt localPointId = i; // new postion in Array
              const MInt globalPointId = offset + i;

              const MInt pIJK = pointIndex(i, m_noGhostLayers, m_noGhostLayers);
              const MInt pIJPK = pointIndex(i, m_noGhostLayers + 1, m_noGhostLayers);

              // compute the normal vector
              MFloat normalVec[3] = {F0, F0, F0};
              for(MInt dim = 0; dim < nDim; dim++) {
                normalVec[dim] = m_grid->m_coordinates[dim][pIJPK] - m_grid->m_coordinates[dim][pIJK];
              }

              // normalize normalVec
              const MFloat normalLength = sqrt(POW2(normalVec[0]) + POW2(normalVec[1]) + POW2(normalVec[2]));
              for(MInt dim = 0; dim < nDim; dim++) {
                normalVec[dim] /= normalLength;
              }

              localCoords(globalPointId) = m_grid->m_coordinates[0][localPointId];
              localCoords(m_airfoilNoWallPoints + globalPointId) = m_grid->m_coordinates[1][localPointId];
              localCoords(2 * m_airfoilNoWallPoints + globalPointId) = m_grid->m_coordinates[2][localPointId];

              localNormalVec(globalPointId) = normalVec[0];
              localNormalVec(m_airfoilNoWallPoints + globalPointId) = normalVec[1];
              localNormalVec(2 * m_airfoilNoWallPoints + globalPointId) = normalVec[2];
            }
          }
        }

        MPI_Allreduce(&localCoords[0], &m_airfoilCoords[0], nDim * m_airfoilNoWallPoints, MPI_DOUBLE, MPI_MAX,
                      m_StructuredComm, AT_, "localCoords[0]", "m_airfoilCoords[0]");
        MPI_Allreduce(&localNormalVec[0], &m_airfoilNormalVec[0], nDim * m_airfoilNoWallPoints, MPI_DOUBLE, MPI_MAX,
                      m_StructuredComm, AT_, "localNormalVec[0]", "m_airfoilNormalVec[0]");

        m_waveSpeed = m_waveLength / m_waveTime;
        const MFloat deltaZ = abs(m_grid->m_coordinates[2][0] - m_grid->m_coordinates[2][m_nPoints[2] * m_nPoints[1]]);
        m_waveCellsPerWaveLength = round(m_waveLength / deltaZ);

        const MFloat wavePeriod = m_waveLength / m_waveSpeed;
        const MFloat speedAmplitude = 2 * PI * m_waveAmplitude / wavePeriod;

        m_log << "/////////////////// TRAVELING WAVE /////////////////////////////" << endl;
        m_log << "Wavelength: " << m_waveLength << " Period: " << wavePeriod << endl;
        m_log << "Amplitude Suction: " << m_waveAmplitudeSuction << " Amplitude Pressure: " << m_waveAmplitudePressure
              << endl;
        m_log << "Amplitude Gradient Suction: " << m_waveGradientSuction
              << " Amplitude Gradient Pressure: " << m_waveGradientPressure << endl;
        m_log << "Velocity Amplitude Suction: " << 2 * PI * m_waveAmplitudeSuction / wavePeriod
              << " Velocity Amplitude Pressure: " << 2 * PI * m_waveAmplitudePressure / wavePeriod << endl;
        m_log << "Speed: " << m_waveSpeed << endl;
        m_log << "Max up/down speed: " << speedAmplitude << endl;
        m_log << " Begin InTransition: " << m_waveBeginTransition << " End InTransition: " << m_waveEndTransition
              << endl;
        m_log << " Begin OutTransition: " << m_waveOutBeginTransition
              << " End OutTransition: " << m_waveOutEndTransition << endl;
        m_log << " Begin WallNormalTransition: " << m_waveYBeginTransition
              << " End WallNormalTransition: " << m_waveYEndTransition << endl;
        m_log << "////////////////////////////////////////////////////////////////" << endl;
        fixTimeStepTravelingWave();
        break;
      }
    case 14:
      // oscillating cylinder
      {
        // time needs to be constant for traveling wave
        m_constantTimeStep = true;

        /*! \property
          \page propertiesFVSTRCTRD
                \section oscAmplitude
                <code>MInt FvStructuredSolver::m_oscAmplitude </code>\n
                default = <code> 1.0 </code>\n \n
                Amplitude of the oscillating cylinder motion.\n
                Possible values are:\n
                <ul>
                <li>Float > 0.0</li>
                </ul>
                Keywords: <i>MOVING, STRUCTURED</i>
              */
        m_oscAmplitude = 0.2;
        m_oscAmplitude = Context::getSolverProperty<MFloat>("oscAmplitude", m_solverId, AT_, &m_oscAmplitude);

        m_oscSr = 0.195;
        m_oscSr = Context::getSolverProperty<MFloat>("oscSr", m_solverId, AT_, &m_oscSr);

        MFloat freqFactor = 0.8;
        freqFactor = Context::getSolverProperty<MFloat>("oscFreqFactor", m_solverId, AT_, &freqFactor);

        const MFloat freq0 = m_oscSr * PV->UInfinity / m_referenceLength;
        m_oscFreq = freq0 * freqFactor;
        break;
      }
    case 15:
      // streamwise traveling wave on airfoil
      {
        m_streamwiseTravelingWave = true;
        m_waveSpeed = 0.0;
        m_waveCellsPerWaveLength = 1;
        if(!m_restart) {
          m_waveTimeStepComputed = false;
        }
        // time needs to be constant for traveling wave
        m_constantTimeStep = true;
        if(!m_restart) {
          m_waveNoStepsPerCell = 1;
        }

        /*! \page propertyPage1
          \section waveLength
          <code>MInt FvStructuredSolver::m_waveLength </code>\n
          default = <code> 1.0 </code>\n \n
          Wavelength of the traveling wave.\n
          Possible values are:\n
          <ul>
          <li>Float > 0.0</li>
          </ul>
          Keywords: <i>WAVE, MOVING, STRUCTURED</i>
        */
        m_waveLength = 0.0;
        m_waveLength = Context::getSolverProperty<MFloat>("waveLength", m_solverId, AT_, &m_waveLength);

        /*! \page propertyPage1
          \section waveAmplitudeSuction
          <code>MInt StrctrdBlck::m_waveAmplitudeSuction </code>\n
          default = <code> 1.0 </code>\n \n
          Amplitude of the traveling wave on the airfoil suction side.\n
          Possible values are:\n
          <ul>
          <li>Float > 0.0</li>
          </ul>
          Keywords: <i>WAVE, MOVING, STRCTRD</i>
        */
        m_waveAmplitudeSuction = 0.0;
        m_waveAmplitudeSuction =
            Context::getSolverProperty<MFloat>("waveAmplitudeSuction", m_blockId, AT_, &m_waveAmplitudeSuction);

        /*! \page propertyPage1
          \section waveAmplitudePressure
          <code>MInt StrctrdBlck::m_waveAmplitudePressure </code>\n
          default = <code> 1.0 </code>\n \n
          Amplitude of the traveling wave on the airfoil suction side.\n
          Possible values are:\n
          <ul>
          <li>Float > 0.0</li>
          </ul>
          Keywords: <i>WAVE, MOVING, STRCTRD</i>
        */
        m_waveAmplitudePressure = 0.0;
        m_waveAmplitudePressure =
            Context::getSolverProperty<MFloat>("waveAmplitudePressure", m_blockId, AT_, &m_waveAmplitudePressure);

        /*! \page propertyPage1
          \section waveTime
          <code>MInt FvStructuredSolver::m_waveTime </code>\n
          default = <code> 1.0 </code>\n \n
          Period time of the traveling wave.\n
          Possible values are:\n
          <ul>
          <li>Float > 0.0</li>
          </ul>
          Keywords: <i>WAVE, MOVING, STRUCTURED</i>
        */
        m_waveTime = 0.0;
        m_waveTime = Context::getSolverProperty<MFloat>("waveTime", m_solverId, AT_, &m_waveTimePlus);

        /*! \page propertyPage1
          \section waveBeginTransition
          <code>MInt FvStructuredSolver::m_waveBeginTransition </code>\n
          default = <code> 1.0 </code>\n \n
          Start of the transition from flat to wave in x-dir.\n
          Possible values are:\n
          <ul>
          <li>Float</li>
          </ul>
          Keywords: <i>WAVE, MOVING, STRUCTURED</i>
        */
        m_waveBeginTransition = 0.0;
        m_waveBeginTransition =
            Context::getSolverProperty<MFloat>("waveBeginTransition", m_solverId, AT_, &m_waveBeginTransition);

        /*! \page propertyPage1
          \section waveEndTransition
          <code>MInt FvStructuredSolver::m_waveEndTransition </code>\n
          default = <code> 1.0 </code>\n \n
          End of the transition from flat to wave in x-dir.\n
          Possible values are:\n
          <ul>
          <li>Float</li>
          </ul>
          Keywords: <i>WAVE, MOVING, STRUCTURED</i>
        */
        m_waveEndTransition = 0.0;
        m_waveEndTransition =
            Context::getSolverProperty<MFloat>("waveEndTransition", m_solverId, AT_, &m_waveEndTransition);

        /*! \page propertyPage1
          \section waveOutBeginTransition
          <code>MInt FvStructuredSolver::m_waveOutBeginTransition </code>\n
          default = <code> 1.0 </code>\n \n
          Start of the transition from wave to flat in x-dir.\n
          Possible values are:\n
          <ul>
          <li>Float</li>
          </ul>
          Keywords: <i>WAVE, MOVING, STRUCTURED</i>
        */
        m_waveOutBeginTransition = 1000000.0;
        m_waveOutBeginTransition =
            Context::getSolverProperty<MFloat>("waveOutBeginTransition", m_solverId, AT_, &m_waveOutBeginTransition);

        /*! \page propertyPage1
          \section waveOutEndTransition
          <code>MInt FvStructuredSolver::m_waveOutEndTransition </code>\n
          default = <code> 1.0 </code>\n \n
          End of the transition from wave to flat in x-dir.\n
          Possible values are:\n
          <ul>
          <li>Float</li>
          </ul>
          Keywords: <i>WAVE, MOVING, STRUCTURED</i>
        */
        m_waveOutEndTransition = 2000000.0;
        m_waveOutEndTransition =
            Context::getSolverProperty<MFloat>("waveOutEndTransition", m_solverId, AT_, &m_waveOutEndTransition);


        /*! \page propertyPage1
          \section wavePressureBeginTransition
          <code>MInt FvStructuredSolver::m_wavePressureBeginTransition </code>\n
          default = <code> 1.0 </code>\n \n
          Start of the transition from flat to wave in x-dir.\n
          Possible values are:\n
          <ul>
          <li>Float</li>
          </ul>
          Keywords: <i>WAVE, MOVING, STRUCTURED</i>
        */
        m_wavePressureBeginTransition = 0.0;
        m_wavePressureBeginTransition = Context::getSolverProperty<MFloat>("wavePressureBeginTransition", m_solverId,
                                                                           AT_, &m_wavePressureBeginTransition);

        /*! \page propertyPage1
          \section wavePressureEndTransition
          <code>MInt FvStructuredSolver::m_wavePressureEndTransition </code>\n
          default = <code> 1.0 </code>\n \n
          End of the transition from flat to wave in x-dir.\n
          Possible values are:\n
          <ul>
          <li>Float</li>
          </ul>
          Keywords: <i>WAVE, MOVING, STRUCTURED</i>
        */
        m_wavePressureEndTransition = 0.0;
        m_wavePressureEndTransition = Context::getSolverProperty<MFloat>("wavePressureEndTransition", m_solverId, AT_,
                                                                         &m_wavePressureEndTransition);

        /*! \page propertyPage1
          \section wavePressureOutBeginTransition
          <code>MInt FvStructuredSolver::m_wavePressureOutBeginTransition </code>\n
          default = <code> 1.0 </code>\n \n
          Start of the transition from wave to flat in x-dir.\n
          Possible values are:\n
          <ul>
          <li>Float</li>
          </ul>
          Keywords: <i>WAVE, MOVING, STRUCTURED</i>
        */
        m_wavePressureOutBeginTransition = 1000000.0;
        m_wavePressureOutBeginTransition = Context::getSolverProperty<MFloat>(
            "wavePressureOutBeginTransition", m_solverId, AT_, &m_wavePressureOutBeginTransition);

        /*! \page propertyPage1
          \section wavePressureOutEndTransition
          <code>MInt FvStructuredSolver::m_wavePressureOutEndTransition </code>\n
          default = <code> 1.0 </code>\n \n
          End of the transition from wave to flat in x-dir.\n
          Possible values are:\n
          <ul>
          <li>Float</li>
          </ul>
          Keywords: <i>WAVE, MOVING, STRUCTURED</i>
        */
        m_wavePressureOutEndTransition = 2000000.0;
        m_wavePressureOutEndTransition = Context::getSolverProperty<MFloat>("wavePressureOutEndTransition", m_solverId,
                                                                            AT_, &m_wavePressureOutEndTransition);

        /*! \page propertyPage1
          \section waveYBeginTransition
          <code>MInt FvStructuredSolver::m_waveYBeginTransition </code>\n
          default = <code> 1.0 </code>\n \n
          End of the transition from wave to flat in x-dir.\n
          Possible values are:\n
          <ul>
          <li>Float</li>
          </ul>
          Keywords: <i>WAVE, MOVING, STRUCTURED</i>
        */
        m_waveYBeginTransition = 2000000.0;
        if(Context::propertyExists("waveYBeginTransition", m_solverId)) {
          m_waveYBeginTransition =
              Context::getSolverProperty<MFloat>("waveYBeginTransition", m_solverId, AT_, &m_waveYBeginTransition);
        }


        /*! \page propertyPage1
          \section waveYEndTransition
          <code>MInt FvStructuredSolver::m_waveYEndTransition </code>\n
          default = <code> 1.0 </code>\n \n
          End of the transition from wave to flat in x-dir.\n
          Possible values are:\n
          <ul>
          <li>Float</li>
          </ul>
          Keywords: <i>WAVE, MOVING, STRUCTURED</i>
        */
        m_waveYEndTransition = 2000000.0;
        if(Context::propertyExists("waveYEndTransition", m_solverId)) {
          m_waveYEndTransition =
              Context::getSolverProperty<MFloat>("waveYEndTransition", m_solverId, AT_, &m_waveYEndTransition);
        }

        const MInt myBlockId = m_grid->getMyBlockId();
        const MInt blockId = 0; // use information from block 0
        m_airfoilNoWallPoints = m_grid->getBlockNoPoints(blockId, 2) + 2 * m_noGhostLayers;

        if(domainId() == 0) {
          cout << "NoWallPoints: " << m_airfoilNoWallPoints << endl;
        }

        MFloatScratchSpace localCoords(nDim, m_airfoilNoWallPoints, AT_, "localCoords");
        MFloatScratchSpace localNormalVec(nDim, m_airfoilNoWallPoints, AT_, "localNormalVec");
        localCoords.fill(-99999.0);
        localNormalVec.fill(-99999.0);
        mAlloc(m_airfoilCoords, nDim * m_airfoilNoWallPoints, "m_airfoilCoords", -99999.0, AT_);
        mAlloc(m_airfoilNormalVec, nDim * m_airfoilNoWallPoints, "m_airfoilNormalVec", -99999.0, AT_);

        if(myBlockId == blockId) {
          if(m_nOffsetCells[1] == 0 && m_nOffsetCells[0] == 0) {
            // only collect variables from points at the wall
            for(MInt i = 0; i < m_nPoints[2]; i++) {
              const MInt offset = m_nOffsetCells[2];
              const MInt localPointId = i; // new postion in Array
              const MInt globalPointId = offset + i;

              const MInt pIJK = pointIndex(i, m_noGhostLayers, m_noGhostLayers);
              const MInt pIJPK = pointIndex(i, m_noGhostLayers + 1, m_noGhostLayers);

              // compute the normal vector
              MFloat normalVec[3] = {F0, F0, F0};
              for(MInt dim = 0; dim < nDim; dim++) {
                normalVec[dim] = m_grid->m_coordinates[dim][pIJPK] - m_grid->m_coordinates[dim][pIJK];
              }

              // normalize normalVec
              const MFloat normalLength = sqrt(POW2(normalVec[0]) + POW2(normalVec[1]) + POW2(normalVec[2]));
              for(MInt dim = 0; dim < nDim; dim++) {
                normalVec[dim] /= normalLength;
              }

              localCoords(globalPointId) = m_grid->m_coordinates[0][localPointId];
              localCoords(m_airfoilNoWallPoints + globalPointId) = m_grid->m_coordinates[1][localPointId];
              localCoords(2 * m_airfoilNoWallPoints + globalPointId) = m_grid->m_coordinates[2][localPointId];

              localNormalVec(globalPointId) = normalVec[0];
              localNormalVec(m_airfoilNoWallPoints + globalPointId) = normalVec[1];
              localNormalVec(2 * m_airfoilNoWallPoints + globalPointId) = normalVec[2];
            }
          }
        }

        MPI_Allreduce(&localCoords[0], &m_airfoilCoords[0], nDim * m_airfoilNoWallPoints, MPI_DOUBLE, MPI_MAX,
                      m_StructuredComm, AT_, "localCoords[0]", "m_airfoilCoords[0]");
        MPI_Allreduce(&localNormalVec[0], &m_airfoilNormalVec[0], nDim * m_airfoilNoWallPoints, MPI_DOUBLE, MPI_MAX,
                      m_StructuredComm, AT_, "localNormalVec[0]", "m_airfoilNormalVec[0]");

        m_waveSpeed = m_waveLength / m_waveTime;
        const MFloat deltaZ = abs(m_grid->m_coordinates[2][0] - m_grid->m_coordinates[2][m_nPoints[2] * m_nPoints[1]]);
        m_waveCellsPerWaveLength = round(m_waveLength / deltaZ);

        const MFloat wavePeriod = m_waveLength / m_waveSpeed;
        const MFloat speedAmplitude = 2 * PI * m_waveAmplitude / wavePeriod;

        m_log << "/////////////////// TRAVELING WAVE /////////////////////////////" << endl;
        m_log << "Wavelength: " << m_waveLength << " Period: " << wavePeriod << endl;
        m_log << "Amplitude Suction: " << m_waveAmplitudeSuction << " Amplitude Pressure: " << m_waveAmplitudePressure
              << endl;
        m_log << "Velocity Amplitude Suction: " << 2 * PI * m_waveAmplitudeSuction / wavePeriod
              << " Velocity Amplitude Pressure: " << 2 * PI * m_waveAmplitudePressure / wavePeriod << endl;
        m_log << "Speed: " << m_waveSpeed << endl;
        m_log << "Max up/down speed: " << speedAmplitude << endl;
        m_log << " Begin InTransition: " << m_waveBeginTransition << " End InTransition: " << m_waveEndTransition
              << endl;
        m_log << " Begin OutTransition: " << m_waveOutBeginTransition
              << " End OutTransition: " << m_waveOutEndTransition << endl;
        m_log << " Begin PressureInTransition: " << m_wavePressureBeginTransition
              << " End PressureInTransition: " << m_wavePressureEndTransition << endl;
        m_log << " Begin PressureOutTransition: " << m_wavePressureOutBeginTransition
              << " End PressureOutTransition: " << m_wavePressureOutEndTransition << endl;
        m_log << " Begin WallNormalTransition: " << m_waveYBeginTransition
              << " End WallNormalTransition: " << m_waveYEndTransition << endl;
        m_log << "////////////////////////////////////////////////////////////////" << endl;
        break;
      }
    default: {
      stringstream errorMessage;
      errorMessage << "Grid moving method " << m_gridMovingMethod << " not implemented!" << endl;
      mTerm(1, AT_, errorMessage.str());
    }
  }

  m_grid->saveGrid();

  // now move the grid to the correct position
  if(m_restart) {
    if(m_movingGrid && !m_stgInitialStartup) {
      if(m_movingGridInitialStart) {
        // if this is an initial start of the
        // grid movement, just move to initial pos
        moveGrid(true, false);
        m_grid->saveGrid();
      } else {
        // move to last pos before restart,
        // save and move to current pos again
        // this way the grid velocity is computed
        // correctly in the BC
        m_time -= m_timeStep;
        moveGrid(true, false);
        m_grid->saveGrid();
        m_time += m_timeStep;
        moveGrid(true, false);
      }
    }
  } else {
    moveGrid(true, false);
    m_grid->saveGrid();
  }

  m_log << "Initializing moving grid methods... DONE!" << endl;
}

void FvStructuredSolver3D::applyBodyForce(const MBool isRestart, const MBool zeroPos) {
  TRACE();
  const MFloat pi = 4.0 * atan(1);
  const MFloat t = (isRestart) ? m_time : m_time + m_timeStep * m_RKalpha[m_RKStep];

  switch(m_bodyForceMethod) {
    case 10: {
      // traveling wave case
      MFloat t_offset = t - m_movingGridTimeOffset;
      if(zeroPos) {
        t_offset = F0;
      }
      const MFloat transitionLength = m_waveEndTransition - m_waveBeginTransition;
      const MFloat transitionOutLength = m_waveOutEndTransition - m_waveOutBeginTransition;

      MFloat fadeInFactor = 0.0;
      if(t_offset < m_waveTemporalTransition) {
        fadeInFactor = (1.0 - cos(t_offset / m_waveTemporalTransition * pi)) * F1B2;
      } else {
        fadeInFactor = 1.0;
      }

      if(zeroPos) {
        fadeInFactor = F1;
      }

      for(MInt k = 0; k < m_nCells[0]; k++) {
        for(MInt j = 0; j < m_nCells[1]; j++) {
          for(MInt i = 0; i < m_nCells[2]; i++) {
            const MInt cellId = cellIndex(i, j, k);
            const MFloat x = m_cells->coordinates[0][cellId];
            const MFloat z = m_cells->coordinates[2][cellId];
            const MFloat y = m_cells->coordinates[1][cellId];

            MFloat transitionFactor = F0;
            if(x <= m_waveBeginTransition) {
              transitionFactor = F0;
            } else if(x > m_waveBeginTransition && x < m_waveEndTransition) {
              transitionFactor = (1 - cos((x - m_waveBeginTransition) / transitionLength * pi)) * F1B2;
            } else if(m_waveEndTransition <= x && x <= m_waveOutBeginTransition) {
              transitionFactor = F1;
            } else if(x > m_waveOutBeginTransition && x < m_waveOutEndTransition) {
              transitionFactor = (1 + cos((x - m_waveOutBeginTransition) / transitionOutLength * pi)) * F1B2;
            } else {
              transitionFactor = F0;
            }

            const MFloat force = fadeInFactor * transitionFactor
                                 * (m_waveAmplitude * exp(-y / m_wavePenetrationHeight)
                                    * cos((F2 * pi) / m_waveLength * (z - m_waveSpeed * t_offset)));

            m_cells->rightHandSide[CV->RHO_W][cellId] +=
                m_cells->variables[CV->RHO][cellId] * force * m_cells->cellJac[cellId];
            m_cells->rightHandSide[CV->RHO_E][cellId] +=
                m_cells->variables[CV->RHO_W][cellId] * force * m_cells->cellJac[cellId];
          }
        }
      }

      break;
    }
    case 11: {
      // traveling wave case with predefined force field
      MFloat t_offset = t - m_movingGridTimeOffset;
      if(zeroPos) {
        t_offset = F0;
      }
      const MFloat transitionLength = m_waveEndTransition - m_waveBeginTransition;
      const MFloat transitionOutLength = m_waveOutEndTransition - m_waveOutBeginTransition;

      MFloat fadeInFactor = 0.0;
      if(t_offset < m_waveTemporalTransition) {
        fadeInFactor = (1.0 - cos(t_offset / m_waveTemporalTransition * pi)) * F1B2;
      } else {
        fadeInFactor = 1.0;
      }

      if(zeroPos) {
        fadeInFactor = F1;
      }

      for(MInt k = m_noGhostLayers; k < m_nCells[0] - m_noGhostLayers; k++) {
        for(MInt j = m_noGhostLayers; j < m_nCells[1] - m_noGhostLayers; j++) {
          for(MInt i = m_noGhostLayers; i < m_nCells[2] - m_noGhostLayers; i++) {
            const MInt cellId = cellIndex(i, j, k);

            const MFloat x = m_cells->coordinates[0][cellId];
            const MFloat z = m_cells->coordinates[2][cellId];

            MFloat transitionFactor = F0;
            if(x <= m_waveBeginTransition) {
              transitionFactor = F0;
            } else if(x > m_waveBeginTransition && x < m_waveEndTransition) {
              transitionFactor = (1 - cos((x - m_waveBeginTransition) / transitionLength * pi)) * F1B2;
            } else if(m_waveEndTransition <= x && x <= m_waveOutBeginTransition) {
              transitionFactor = F1;
            } else if(x > m_waveOutBeginTransition && x < m_waveOutEndTransition) {
              transitionFactor = (1 + cos((x - m_waveOutBeginTransition) / transitionOutLength * pi)) * F1B2;
            } else {
              transitionFactor = F0;
            }

            const MFloat pos_abs = z - m_waveSpeed * t_offset;
            const MFloat pos = fmod(pos_abs, m_waveDomainWidth);

            const MInt jj = m_nOffsetCells[1] + (j - m_noGhostLayers);
            MFloat amp = 0.0;

            if(pos < m_waveForceZ[jj]) {
              const MInt p1 = jj;
              const MInt p1end = jj + (m_grid->getMyBlockNoCells(0) - 1) * m_grid->getMyBlockNoCells(1);
              const MFloat zz = m_waveForceZ[p1];
              const MFloat zzstart = m_waveForceZ[p1end] - m_waveDomainWidth;
              const MFloat f = m_waveForceField[p1];
              const MFloat fstart = m_waveForceField[p1end];
              amp = fstart + (pos - zzstart) / (zz - zzstart) * (f - fstart);
            } else if(pos > m_waveForceZ[jj + (m_grid->getMyBlockNoCells(0) - 1) * m_grid->getMyBlockNoCells(1)]) {
              const MInt p1 = jj + (m_grid->getMyBlockNoCells(0) - 1) * m_grid->getMyBlockNoCells(1);
              const MInt p1start = jj + 0 * m_grid->getMyBlockNoCells(1);
              const MFloat zz = m_waveForceZ[p1];
              const MFloat zzend = m_waveForceZ[p1start] + m_waveDomainWidth;
              const MFloat f = m_waveForceField[p1];
              const MFloat fend = m_waveForceField[p1start];
              amp = f + (pos - zz) / (zzend - zz) * (fend - f);
            } else {
              for(MInt kk = 0; kk < m_grid->getMyBlockNoCells(0) - 1; ++kk) {
                const MInt p1 = jj + kk * m_grid->getMyBlockNoCells(1);
                const MInt p1p = jj + (kk + 1) * m_grid->getMyBlockNoCells(1);
                const MFloat zz = m_waveForceZ[p1];
                const MFloat zzp = m_waveForceZ[p1p];
                const MFloat f = m_waveForceField[p1];
                const MFloat fp = m_waveForceField[p1p];
                if(zz <= pos && pos < zzp) {
                  amp = f + (pos - zz) / (zzp - zz) * (fp - f);
                  break;
                }
              }
            }

            const MFloat force = fadeInFactor * transitionFactor * amp * m_waveAmplitude;

            m_cells->rightHandSide[CV->RHO_W][cellId] +=
                m_cells->variables[CV->RHO][cellId] * force * m_cells->cellJac[cellId];
            m_cells->rightHandSide[CV->RHO_E][cellId] +=
                m_cells->variables[CV->RHO_W][cellId] * force * m_cells->cellJac[cellId];
          }
        }
      }

      break;
    }
    default: {
      mTerm(1, AT_, "Body Force Method not implemented!");
    }
  }
}

void FvStructuredSolver3D::initBodyForce() {
  TRACE();
  m_log << "Initializing body force methods..." << endl;

  // Second approach: save only parts of the mesh depending on moving grid case
  switch(m_bodyForceMethod) {
    case 10:
      // travelling wave defined by non-plus units
      {
        m_travelingWave = true;
        if(!m_restart) {
          m_waveTimeStepComputed = false;
        }
        m_waveSpeed = 0.0;
        m_waveLength = 0.0;
        m_waveAmplitude = 0.0;
        m_waveCellsPerWaveLength = 1;
        if(!m_restart) {
          m_waveNoStepsPerCell = 1;
        }

        // time needs to be constant for traveling wave
        m_constantTimeStep = true;

        /*! \property
          \page propertiesFVSTRCTRD
                \section waveLength
                <code>MInt FvStructuredSolver::m_waveLength </code>\n
                default = <code> 1.0 </code>\n \n
                Wavelength of the traveling wave.\n
                Possible values are:\n
                <ul>
                <li>Float > 0.0</li>
                </ul>
                Keywords: <i>WAVE, MOVING, STRUCTURED</i>
              */
        m_waveLength = 0.0;
        m_waveLength = Context::getSolverProperty<MFloat>("waveLength", m_solverId, AT_, &m_waveLength);

        /*! \property
          \page propertiesFVSTRCTRD
                \section waveAmplitude
                <code>MInt FvStructuredSolver::m_waveAmplitude </code>\n
                default = <code> 1.0 </code>\n \n
                Amplitude of the traveling wave.\n
                Possible values are:\n
                <ul>
                <li>Float > 0.0</li>
                </ul>
                Keywords: <i>WAVE, MOVING, STRUCTURED</i>
              */
        m_waveAmplitude = 0.0;
        m_waveAmplitude = Context::getSolverProperty<MFloat>("waveAmplitude", m_solverId, AT_, &m_waveAmplitude);

        m_waveAmplitude = m_waveAmplitude * CV->rhoInfinity * POW2(PV->UInfinity);

        /*! \property
          \page propertiesFVSTRCTRD
                \section wavePenetrationHeight
                <code>MInt FvStructuredSolver::m_wavePenetrationHeight </code>\n
                default = <code> 1.0 </code>\n \n
                Period time of the traveling wave.\n
                Possible values are:\n
                <ul>
                <li>Float > 0.0</li>
                </ul>
                Keywords: <i>WAVE, MOVING, STRUCTURED</i>
              */
        m_wavePenetrationHeight = 0.0;
        m_wavePenetrationHeight =
            Context::getSolverProperty<MFloat>("wavePenetrationHeight", m_solverId, AT_, &m_wavePenetrationHeight);

        /*! \property
          \page propertiesFVSTRCTRD
                \section waveTime
                <code>MInt FvStructuredSolver::m_waveTime </code>\n
                default = <code> 1.0 </code>\n \n
                Period time of the traveling wave.\n
                Possible values are:\n
                <ul>
                <li>Float > 0.0</li>
                </ul>
                Keywords: <i>WAVE, MOVING, STRUCTURED</i>
              */
        m_waveTime = 0.0;
        m_waveTime = Context::getSolverProperty<MFloat>("waveTime", m_solverId, AT_, &m_waveTime);


        /*! \property
          \page propertiesFVSTRCTRD
                \section waveBeginTransition
                <code>MInt FvStructuredSolver::m_waveBeginTransition </code>\n
                default = <code> 1.0 </code>\n \n
                Start of the transition from flat to wave in x-dir.\n
                Possible values are:\n
                <ul>
                <li>Float</li>
                </ul>
                Keywords: <i>WAVE, MOVING, STRUCTURED</i>
              */
        m_waveBeginTransition = 0.0;
        m_waveBeginTransition =
            Context::getSolverProperty<MFloat>("waveBeginTransition", m_solverId, AT_, &m_waveBeginTransition);

        /*! \property
          \page propertiesFVSTRCTRD
                \section waveEndTransition
                <code>MInt FvStructuredSolver::m_waveEndTransition </code>\n
                default = <code> 1.0 </code>\n \n
                End of the transition from flat to wave in x-dir.\n
                Possible values are:\n
                <ul>
                <li>Float</li>
                </ul>
                Keywords: <i>WAVE, MOVING, STRUCTURED</i>
              */
        m_waveEndTransition = 0.0;
        m_waveEndTransition =
            Context::getSolverProperty<MFloat>("waveEndTransition", m_solverId, AT_, &m_waveEndTransition);

        /*! \property
          \page propertiesFVSTRCTRD
                \section waveOutBeginTransition
                <code>MInt FvStructuredSolver::m_waveOutBeginTransition </code>\n
                default = <code> 1.0 </code>\n \n
                Start of the transition from wave to flat in x-dir.\n
                Possible values are:\n
                <ul>
                <li>Float</li>
                </ul>
                Keywords: <i>WAVE, MOVING, STRUCTURED</i>
              */
        m_waveOutBeginTransition = 1000000.0;
        m_waveOutBeginTransition =
            Context::getSolverProperty<MFloat>("waveOutBeginTransition", m_solverId, AT_, &m_waveOutBeginTransition);

        /*! \property
          \page propertiesFVSTRCTRD
                \section waveOutEndTransition
                <code>MInt FvStructuredSolver::m_waveOutEndTransition </code>\n
                default = <code> 1.0 </code>\n \n
                End of the transition from wave to flat in x-dir.\n
                Possible values are:\n
                <ul>
                <li>Float</li>
                </ul>
                Keywords: <i>WAVE, MOVING, STRUCTURED</i>
              */
        m_waveOutEndTransition = 2000000.0;
        m_waveOutEndTransition =
            Context::getSolverProperty<MFloat>("waveOutEndTransition", m_solverId, AT_, &m_waveOutEndTransition);


        /*! \property
          \page propertiesFVSTRCTRD
                \section waveTemporalTransition
                <code>MInt FvStructuredSolver::m_waveOutEndTransition </code>\n
                default = <code> 500.0 </code>\n \n
                Acoustic time for wave actuation to transiate from flat plate\n
                to fully extended.\n
                Possible values are:\n
                <ul>
                <li>Float</li>
                </ul>
                Keywords: <i>WAVE, MOVING, STRUCTURED</i>
              */
        m_waveTemporalTransition = 500.0;
        if(Context::propertyExists("waveTemporalTransition", m_solverId)) {
          m_waveTemporalTransition =
              Context::getSolverProperty<MFloat>("waveTemporalTransition", m_solverId, AT_, &m_waveOutEndTransition);
        }

        m_waveSpeed = m_waveLength / m_waveTime;
        const MFloat deltaZ = abs(m_grid->m_coordinates[2][0] - m_grid->m_coordinates[2][m_nPoints[2] * m_nPoints[1]]);
        m_waveCellsPerWaveLength = round(m_waveLength / deltaZ);

        const MFloat speedAmplitude = 2 * PI * m_waveAmplitude / m_waveTime;

        m_log << "/////////////////// TRAVELING WAVE BODY FORCE //////////////////" << endl;
        m_log << "Wavelength: " << m_waveLength << " Amplitude: " << m_waveAmplitude << " Period: " << m_waveTime
              << " Speed: " << m_waveSpeed << endl;
        m_log << "Max up/down speed: " << speedAmplitude << endl;
        m_log << "Max up/down speed: " << m_waveSpeed * m_waveAmplitude << endl;
        m_log << "////////////////////////////////////////////////////////////////" << endl;

        fixTimeStepTravelingWave();
        break;
      }
    case 11:
      // body forcing strength defined by given distribution read from HDF5 file
      {
        m_travelingWave = true;
        if(!m_restart) {
          m_waveTimeStepComputed = false;
        }
        m_waveSpeed = 0.0;
        m_waveLength = 0.0;
        m_waveAmplitude = 0.0;
        m_waveCellsPerWaveLength = 1;
        if(!m_restart) {
          m_waveNoStepsPerCell = 1;
        }

        // time needs to be constant for traveling wave
        m_constantTimeStep = true;

        /*! \property
          \page propertiesFVSTRCTRD
                \section waveLength
                <code>MInt FvStructuredSolver::m_waveLength </code>\n
                default = <code> 1.0 </code>\n \n
                Wavelength of the traveling wave.\n
                Possible values are:\n
                <ul>
                <li>Float > 0.0</li>
                </ul>
                Keywords: <i>WAVE, MOVING, STRUCTURED</i>
              */
        m_waveLength = 0.0;
        m_waveLength = Context::getSolverProperty<MFloat>("waveLength", m_solverId, AT_, &m_waveLength);

        /*! \property
          \page propertiesFVSTRCTRD
                \section waveAmplitude
                <code>MInt FvStructuredSolver::m_waveAmplitude </code>\n
                default = <code> 1.0 </code>\n \n
                Amplitude of the traveling wave.\n
                Possible values are:\n
                <ul>
                <li>Float > 0.0</li>
                </ul>
                Keywords: <i>WAVE, MOVING, STRUCTURED</i>
              */
        m_waveAmplitude = 0.0;
        m_waveAmplitude = Context::getSolverProperty<MFloat>("waveAmplitude", m_solverId, AT_, &m_waveAmplitude);

        m_waveAmplitude = m_waveAmplitude * CV->rhoInfinity * POW2(PV->UInfinity);

        /*! \property
          \page propertiesFVSTRCTRD
                \section waveTime
                <code>MInt FvStructuredSolver::m_waveTime </code>\n
                default = <code> 1.0 </code>\n \n
                Period time of the traveling wave.\n
                Possible values are:\n
                <ul>
                <li>Float > 0.0</li>
                </ul>
                Keywords: <i>WAVE, MOVING, STRUCTURED</i>
              */
        m_waveTime = 0.0;
        m_waveTime = Context::getSolverProperty<MFloat>("waveTime", m_solverId, AT_, &m_waveTime);


        /*! \property
          \page propertiesFVSTRCTRD
                \section waveBeginTransition
                <code>MInt FvStructuredSolver::m_waveBeginTransition </code>\n
                default = <code> 1.0 </code>\n \n
                Start of the transition from flat to wave in x-dir.\n
                Possible values are:\n
                <ul>
                <li>Float</li>
                </ul>
                Keywords: <i>WAVE, MOVING, STRUCTURED</i>
              */
        m_waveBeginTransition = 0.0;
        m_waveBeginTransition =
            Context::getSolverProperty<MFloat>("waveBeginTransition", m_solverId, AT_, &m_waveBeginTransition);

        /*! \property
          \page propertiesFVSTRCTRD
                \section waveEndTransition
                <code>MInt FvStructuredSolver::m_waveEndTransition </code>\n
                default = <code> 1.0 </code>\n \n
                End of the transition from flat to wave in x-dir.\n
                Possible values are:\n
                <ul>
                <li>Float</li>
                </ul>
                Keywords: <i>WAVE, MOVING, STRUCTURED</i>
              */
        m_waveEndTransition = 0.0;
        m_waveEndTransition =
            Context::getSolverProperty<MFloat>("waveEndTransition", m_solverId, AT_, &m_waveEndTransition);

        /*! \property
          \page propertiesFVSTRCTRD
                \section waveOutBeginTransition
                <code>MInt FvStructuredSolver::m_waveOutBeginTransition </code>\n
                default = <code> 1.0 </code>\n \n
                Start of the transition from wave to flat in x-dir.\n
                Possible values are:\n
                <ul>
                <li>Float</li>
                </ul>
                Keywords: <i>WAVE, MOVING, STRUCTURED</i>
              */
        m_waveOutBeginTransition = 1000000.0;
        m_waveOutBeginTransition =
            Context::getSolverProperty<MFloat>("waveOutBeginTransition", m_solverId, AT_, &m_waveOutBeginTransition);

        /*! \property
          \page propertiesFVSTRCTRD
                \section waveOutEndTransition
                <code>MInt FvStructuredSolver::m_waveOutEndTransition </code>\n
                default = <code> 1.0 </code>\n \n
                End of the transition from wave to flat in x-dir.\n
                Possible values are:\n
                <ul>
                <li>Float</li>
                </ul>
                Keywords: <i>WAVE, MOVING, STRUCTURED</i>
              */
        m_waveOutEndTransition = 2000000.0;
        m_waveOutEndTransition =
            Context::getSolverProperty<MFloat>("waveOutEndTransition", m_solverId, AT_, &m_waveOutEndTransition);


        /*! \property
          \page propertiesFVSTRCTRD
                \section waveTemporalTransition
                <code>MInt FvStructuredSolver::m_waveOutEndTransition </code>\n
                default = <code> 500.0 </code>\n \n
                Acoustic time for wave actuation to transiate from flat plate\n
                to fully extended.\n
                Possible values are:\n
                <ul>
                <li>Float</li>
                </ul>
                Keywords: <i>WAVE, MOVING, STRUCTURED</i>
              */
        m_waveTemporalTransition = 500.0;
        if(Context::propertyExists("waveTemporalTransition", m_solverId)) {
          m_waveTemporalTransition =
              Context::getSolverProperty<MFloat>("waveTemporalTransition", m_solverId, AT_, &m_waveOutEndTransition);
        }

        /*! \property
          \page propertiesFVSTRCTRD
                \section forceField
                <code>MString FvStructuredSolver::m_waveForceField </code>\n
                default = <code>"./out"</code>\n \n
                Distribution of the forcing strenth.\n
                Keywords: <i>SOLUTION, IO, STRUCTURED</i>
              */
        m_waveForceFieldFile = "./slice.hdf5";
        m_waveForceFieldFile = Context::getSolverProperty<MString>("forceFieldFile", m_solverId, AT_);

        ParallelIo::size_type size[2] = {m_grid->getMyBlockNoCells(0), m_grid->getMyBlockNoCells(1)};
        mAlloc(m_waveForceField, size[0] * size[1], "m_waveForceField", 0.0, AT_);
        mAlloc(m_waveForceY, size[0] * size[1], "m_waveForceY", 0.0, AT_);
        mAlloc(m_waveForceZ, size[0] * size[1], "m_waveForceZ", 0.0, AT_);


        if(domainId() == 0) {
          ParallelIoHdf5 pio(m_waveForceFieldFile, maia::parallel_io::PIO_READ, MPI_COMM_SELF);
          stringstream pathStr;
          pathStr << "/";
          std::string path = pathStr.str();

          ParallelIo::size_type offset[2] = {0, 0};

          pio.readArray(&m_waveForceField[0], path, "force_amp", 2, offset, size);
          pio.readArray(&m_waveForceY[0], path, "y", 2, offset, size);
          pio.readArray(&m_waveForceZ[0], path, "z", 2, offset, size);
        }

        MPI_Bcast(&m_waveForceField[0], size[0] * size[1], MPI_DOUBLE, 0, m_StructuredComm, AT_, "m_waveForceField");
        MPI_Bcast(&m_waveForceY[0], size[0] * size[1], MPI_DOUBLE, 0, m_StructuredComm, AT_, "m_waveForceY");
        MPI_Bcast(&m_waveForceZ[0], size[0] * size[1], MPI_DOUBLE, 0, m_StructuredComm, AT_, "m_waveForceZ");


        m_waveSpeed = m_waveLength / m_waveTime;
        const MFloat deltaZ = abs(m_grid->m_coordinates[2][0] - m_grid->m_coordinates[2][m_nPoints[2] * m_nPoints[1]]);
        m_waveCellsPerWaveLength = round(m_waveLength / deltaZ);

        const MFloat speedAmplitude = 2 * PI * m_waveAmplitude / m_waveTime;

        m_log << "/////////////////// TRAVELING WAVE BODY FORCE //////////////////" << endl;
        m_log << "Wavelength: " << m_waveLength << " Amplitude: " << m_waveAmplitude << " Period: " << m_waveTime
              << " Speed: " << m_waveSpeed << endl;
        m_log << "Max up/down speed: " << speedAmplitude << endl;
        m_log << "Max up/down speed: " << m_waveSpeed * m_waveAmplitude << endl;
        m_log << "////////////////////////////////////////////////////////////////" << endl;

        MFloat localDomainWidth = m_cells->coordinates[2][cellIndex(0, 0, m_nCells[0] - m_noGhostLayers)];
        MPI_Allreduce(&localDomainWidth, &m_waveDomainWidth, 1, MPI_DOUBLE, MPI_MAX, m_StructuredComm, AT_,
                      "localDomainWidth", "m_waveDomainWidth");

        fixTimeStepTravelingWave();
        break;
      }
    default: {
      stringstream errorMessage;
      errorMessage << "Grid moving method " << m_gridMovingMethod << " not implemented!" << endl;
      mTerm(1, AT_, errorMessage.str());
    }
  }

  m_log << "Initializing body force methods... DONE!" << endl;
}


void FvStructuredSolver3D::saveNodalBoxes() {
  stringstream filename;

  filename << m_nodalBoxOutputDir << "nodalBoxOutput" << m_outputIterationNumber << m_outputFormat;

  ParallelIoHdf5 pio(filename.str(), maia::parallel_io::PIO_REPLACE, m_StructuredComm);

  writeHeaderAttributes(&pio, "boxes");
  writePropertiesAsAttributes(&pio, "");
  pio.setAttribute(m_nodalBoxNoBoxes, "noBoxes", "");

  if(!m_nodalBoxInitialized) {
    m_nodalBoxInterpolation = make_unique<StructuredInterpolation<nDim>>(m_nCells, m_cells->coordinates,
                                                                         m_cells->pvariables, m_StructuredComm);

    mAlloc(m_nodalBoxLocalPoints, m_nodalBoxNoBoxes, nDim, "m_nodalBoxLocalPoints", 0, AT_);
    mAlloc(m_nodalBoxLocalOffset, m_nodalBoxNoBoxes, nDim, "m_nodalBoxLocalOffset", 0, AT_);
    mAlloc(m_nodalBoxLocalDomainOffset, m_nodalBoxNoBoxes, nDim, "m_nodalBoxLocalDomainOffset", 0, AT_);
    mAlloc(m_nodalBoxLocalSize, m_nodalBoxNoBoxes, "m_nodalBoxLocalSize", 0, AT_);
    mAlloc(m_nodalBoxTotalLocalOffset, m_nodalBoxNoBoxes, "m_nodalBoxLocalSize", 0, AT_);

    m_nodalBoxTotalLocalSize = 0;
    for(MInt b = 0; b < m_nodalBoxNoBoxes; ++b) {
      if(m_nodalBoxBlock[b] == m_blockId
         && ((m_nOffsetPoints[2] <= m_nodalBoxOffset[b][2]
              && m_nodalBoxOffset[b][2] < m_nOffsetPoints[2] + m_nActivePoints[2])
             || (m_nodalBoxOffset[b][2] <= m_nOffsetPoints[2]
                 && m_nOffsetPoints[2] < m_nodalBoxOffset[b][2] + m_nodalBoxPoints[b][2]))
         && ((m_nOffsetPoints[1] <= m_nodalBoxOffset[b][1]
              && m_nodalBoxOffset[b][1] < m_nOffsetPoints[1] + m_nActivePoints[1])
             || (m_nodalBoxOffset[b][1] <= m_nOffsetPoints[1]
                 && m_nOffsetPoints[1] < m_nodalBoxOffset[b][1] + m_nodalBoxPoints[b][1]))
         && ((m_nOffsetPoints[0] <= m_nodalBoxOffset[b][0]
              && m_nodalBoxOffset[b][0] < m_nOffsetPoints[0] + m_nActivePoints[0])
             || (m_nodalBoxOffset[b][0] <= m_nOffsetPoints[0]
                 && m_nOffsetPoints[0]
                        < m_nodalBoxOffset[b][0] + m_nodalBoxPoints[b][0]))) { // the nodalBox is contained

        for(MInt dim = 0; dim < nDim; ++dim) {
          if(m_nOffsetPoints[dim] <= m_nodalBoxOffset[b][dim]
             && m_nodalBoxOffset[b][dim] + m_nodalBoxPoints[b][dim] < m_nOffsetPoints[dim] + m_nActivePoints[dim]) {
            m_nodalBoxLocalPoints[b][dim] = m_nodalBoxPoints[b][dim];
            m_nodalBoxLocalOffset[b][dim] = 0;
            m_nodalBoxLocalDomainOffset[b][dim] = m_nodalBoxOffset[b][dim] - m_nOffsetPoints[dim];
          } else if(m_nOffsetPoints[dim] <= m_nodalBoxOffset[b][dim]) {
            m_nodalBoxLocalPoints[b][dim] = (m_nOffsetPoints[dim] + m_nActivePoints[dim]) - m_nodalBoxOffset[b][dim];
            m_nodalBoxLocalOffset[b][dim] = 0;
            m_nodalBoxLocalDomainOffset[b][dim] = m_nodalBoxOffset[b][dim] - m_nOffsetPoints[dim];
          } else if(m_nodalBoxOffset[b][dim] <= m_nOffsetPoints[dim]
                    && m_nOffsetPoints[dim] + m_nActivePoints[dim]
                           < m_nodalBoxOffset[b][dim] + m_nodalBoxPoints[b][dim]) {
            m_nodalBoxLocalPoints[b][dim] = m_nActivePoints[dim];
            m_nodalBoxLocalOffset[b][dim] = m_nOffsetPoints[dim] - m_nodalBoxOffset[b][dim];
            m_nodalBoxLocalDomainOffset[b][dim] = 0;
          } else {
            m_nodalBoxLocalPoints[b][dim] =
                (m_nodalBoxOffset[b][dim] + m_nodalBoxPoints[b][dim]) - m_nOffsetPoints[dim];
            m_nodalBoxLocalOffset[b][dim] = m_nOffsetPoints[dim] - m_nodalBoxOffset[b][dim];
            m_nodalBoxLocalDomainOffset[b][dim] = 0;
          }
        }

        m_nodalBoxLocalSize[b] = 1;
        for(MInt dim = 0; dim < nDim; ++dim) {
          m_nodalBoxLocalSize[b] *= m_nodalBoxLocalPoints[b][dim];
        }

        m_nodalBoxTotalLocalSize += m_nodalBoxLocalSize[b];
      }

      if(b < m_nodalBoxNoBoxes - 1) {
        m_nodalBoxTotalLocalOffset[b + 1] = m_nodalBoxTotalLocalOffset[b] + m_nodalBoxLocalSize[b];
      }
    }

    if(m_nodalBoxTotalLocalSize > 0) {
      mAlloc(m_nodalBoxCoordinates, nDim, m_nodalBoxTotalLocalSize, "m_nodalBoxCoordinates", 0.0, AT_);
      mAlloc(m_nodalBoxVariables, m_maxNoVariables, m_nodalBoxTotalLocalSize, "m_nodalBoxVariables", 0.0, AT_);
      mAlloc(m_nodalBoxPartnerLocal, m_nodalBoxTotalLocalSize, "m_nodalBoxPartnerLocal", 0, AT_);
    }

    for(MInt b = 0; b < m_nodalBoxNoBoxes; ++b) {
      for(MInt k = m_noGhostLayers + m_nodalBoxLocalDomainOffset[b][0];
          k < m_noGhostLayers + m_nodalBoxLocalDomainOffset[b][0] + m_nodalBoxLocalPoints[b][0];
          ++k) {
        for(MInt j = m_noGhostLayers + m_nodalBoxLocalDomainOffset[b][1];
            j < m_noGhostLayers + m_nodalBoxLocalDomainOffset[b][1] + m_nodalBoxLocalPoints[b][1];
            ++j) {
          for(MInt i = m_noGhostLayers + m_nodalBoxLocalDomainOffset[b][2];
              i < m_noGhostLayers + m_nodalBoxLocalDomainOffset[b][2] + m_nodalBoxLocalPoints[b][2];
              ++i) {
            const MInt pointId = i + (j + k * m_nPoints[1]) * m_nPoints[2];
            const MInt nodalBoxI = i - m_noGhostLayers - m_nodalBoxLocalDomainOffset[b][2];
            const MInt nodalBoxJ = j - m_noGhostLayers - m_nodalBoxLocalDomainOffset[b][1];
            const MInt nodalBoxK = k - m_noGhostLayers - m_nodalBoxLocalDomainOffset[b][0];
            const MInt localTotalId =
                m_nodalBoxTotalLocalOffset[b]
                + (nodalBoxI + (nodalBoxJ + nodalBoxK * m_nodalBoxLocalPoints[b][1]) * m_nodalBoxLocalPoints[b][2]);

            for(MInt dim = 0; dim < nDim; dim++) {
              m_nodalBoxCoordinates[dim][localTotalId] = m_grid->m_coordinates[dim][pointId];
            }
          }
        }
      }
    }

    m_nodalBoxInterpolation->prepareInterpolation(m_nodalBoxTotalLocalSize, m_nodalBoxCoordinates,
                                                  m_nodalBoxPartnerLocal);
    m_nodalBoxInitialized = true;
  }


  for(MInt b = 0; b < m_nodalBoxNoBoxes; ++b) {
    stringstream pathName;
    pathName << "box" << b;

    pio.setAttribute(m_nodalBoxOffset[b][2], "offseti", pathName.str());
    pio.setAttribute(m_nodalBoxOffset[b][1], "offsetj", pathName.str());
    pio.setAttribute(m_nodalBoxOffset[b][0], "offsetk", pathName.str());

    pio.setAttribute(m_nodalBoxPoints[b][2], "sizei", pathName.str());
    pio.setAttribute(m_nodalBoxPoints[b][1], "sizej", pathName.str());
    pio.setAttribute(m_nodalBoxPoints[b][0], "sizek", pathName.str());

    pio.setAttribute(m_nodalBoxBlock[b], "blockId", pathName.str());

    MInt hasCoordinates = 0;

    ParallelIo::size_type size[3] = {m_nodalBoxPoints[b][0], m_nodalBoxPoints[b][1], m_nodalBoxPoints[b][2]};
    for(MInt v = 0; v < m_maxNoVariables; v++) {
      pio.defineArray(maia::parallel_io::PIO_FLOAT, pathName.str(), m_pvariableNames[v], 3, size);
    }

    // create datasets for the variables
    if(m_nodalBoxWriteCoordinates) {
      hasCoordinates = 1;
      pio.defineArray(maia::parallel_io::PIO_FLOAT, pathName.str(), "x", 3, size);
      pio.defineArray(maia::parallel_io::PIO_FLOAT, pathName.str(), "y", 3, size);
      pio.defineArray(maia::parallel_io::PIO_FLOAT, pathName.str(), "z", 3, size);
    }
    // write output to check if coordinates are contained within the variable list
    pio.setAttribute(hasCoordinates, "hasCoordinates", pathName.str());
  }


  if(m_nodalBoxTotalLocalSize > 0) {
    // interpolate primitive variables to nodal grid points
    m_nodalBoxInterpolation->interpolateVariables(m_nodalBoxVariables);
  }


  for(MInt b = 0; b < m_nodalBoxNoBoxes; ++b) {
    // check if the box is contained your inputsolverId
    ParallelIo::size_type localOffset[3] = {m_nodalBoxLocalOffset[b][0], m_nodalBoxLocalOffset[b][1],
                                            m_nodalBoxLocalOffset[b][2]};
    ParallelIo::size_type localSize[3] = {m_nodalBoxLocalPoints[b][0], m_nodalBoxLocalPoints[b][1],
                                          m_nodalBoxLocalPoints[b][2]};
    stringstream pathName;
    pathName << "box" << b;
    if(m_nodalBoxLocalSize[b] > 0) {
      for(MInt v = 0; v < m_maxNoVariables; v++) {
        pio.writeArray(&m_nodalBoxVariables[v][m_nodalBoxTotalLocalOffset[b]], pathName.str(), m_pvariableNames[v],
                       nDim, localOffset, localSize);
      }

      if(m_nodalBoxWriteCoordinates) {
        pio.writeArray(&m_nodalBoxCoordinates[0][m_nodalBoxTotalLocalOffset[b]], pathName.str(), "x", nDim, localOffset,
                       localSize);
        pio.writeArray(&m_nodalBoxCoordinates[1][m_nodalBoxTotalLocalOffset[b]], pathName.str(), "y", nDim, localOffset,
                       localSize);
        pio.writeArray(&m_nodalBoxCoordinates[2][m_nodalBoxTotalLocalOffset[b]], pathName.str(), "z", nDim, localOffset,
                       localSize);
      }
    } else { // write out nothing as box is not contained
      ParallelIo::size_type localBoxPoints[3] = {0, 0, 0};
      ParallelIo::size_type localBoxOffset[3] = {0, 0, 0};
      MFloat empty = 0;

      for(MInt v = 0; v < m_maxNoVariables; ++v) {
        pio.writeArray(&empty, pathName.str(), m_pvariableNames[v], nDim, localBoxOffset, localBoxPoints);
      }

      if(m_nodalBoxWriteCoordinates) {
        pio.writeArray(&empty, pathName.str(), "x", nDim, localBoxOffset, localBoxPoints);
        pio.writeArray(&empty, pathName.str(), "y", nDim, localBoxOffset, localBoxPoints);
        pio.writeArray(&empty, pathName.str(), "z", nDim, localBoxOffset, localBoxPoints);
      }
    }
  }
}


void FvStructuredSolver3D::loadRestartBC2600() {
  if(m_bc2600IsActive && !m_bc2600InitialStartup) {
    if(domainId() == 0) {
      cout << "Loading BC2600 values..." << endl;
    }
    if(m_bc2600) {
      ParallelIo::size_type bcCells[3] = {m_grid->getMyBlockNoCells(0), m_grid->getMyBlockNoCells(1), m_noGhostLayers};
      MInt noCellsBC = bcCells[0] * bcCells[1] * bcCells[2];
      ParallelIo::size_type bcOffset[3] = {0, 0, 0};
      MFloatScratchSpace tmpRestartVars(noCellsBC * m_maxNoVariables, AT_, "m_tmpRestartVars2600");
      if(m_commBC2600MyRank == 0) {
        stringstream restartFileName;
        MString restartFile = Context::getSolverProperty<MString>("restartVariablesFileName", m_solverId, AT_);
        restartFileName << outputDir() << restartFile;

        ParallelIoHdf5 pio(restartFileName.str(), maia::parallel_io::PIO_READ, MPI_COMM_SELF);
        stringstream pathStr;
        pathStr << "/block" << m_blockId << "/bc2600";
        std::string path = pathStr.str();

        for(MInt var = 0; var < m_maxNoVariables; var++) {
          pio.readArray(&tmpRestartVars[var * noCellsBC], path, m_pvariableNames[var], nDim, bcOffset, bcCells);
        }
      }

      MPI_Bcast(&tmpRestartVars(0, 0), noCellsBC * m_maxNoVariables, MPI_DOUBLE, 0, m_commBC2600, AT_,
                "tmpRestartVars(0");

      MInt startGC[3] = {0, 0, 0};
      MInt endGC[3] = {0, 0, 0};

      if(m_bc2600noOffsetCells[1] == 0) {
        startGC[1] = m_noGhostLayers;
      }
      if(m_bc2600noOffsetCells[0] == 0) {
        startGC[0] = m_noGhostLayers;
      }
      if(m_bc2600noOffsetCells[1] + m_bc2600noActiveCells[1] == bcCells[1]) {
        endGC[1] = m_noGhostLayers;
      }
      if(m_bc2600noOffsetCells[0] + m_bc2600noActiveCells[0] == bcCells[0]) {
        endGC[0] = m_noGhostLayers;
      }

      for(MInt i = 0; i < m_noGhostLayers; i++) {
        for(MInt j = startGC[1]; j < m_bc2600noCells[1] - endGC[1]; j++) {
          for(MInt k = startGC[0]; k < m_bc2600noCells[0] - endGC[0]; k++) {
            const MInt cellId = i + (j + k * m_bc2600noCells[1]) * m_bc2600noCells[2];
            MInt globalI = i;
            MInt globalJ = m_bc2600noOffsetCells[1] - m_noGhostLayers + j;
            MInt globalK = m_bc2600noOffsetCells[0] - m_noGhostLayers + k;
            MInt cellIdBC = globalI + (globalJ + globalK * bcCells[1]) * bcCells[2];

            // load values from restart field
            for(MInt var = 0; var < m_maxNoVariables; var++) {
              m_bc2600Variables[var][cellId] = tmpRestartVars[var * noCellsBC + cellIdBC];
            }
          }
        }
      }


      for(MInt i = 0; i < m_bc2600noCells[2]; i++) {
        for(MInt j = 0; j < m_bc2600noCells[1]; j++) {
          MInt cellIdA1 = i + (j + 2 * m_bc2600noCells[1]) * m_bc2600noCells[2];
          MInt cellIdG1 = i + (j + 0 * m_bc2600noCells[1]) * m_bc2600noCells[2];
          MInt cellIdG2 = i + (j + 1 * m_bc2600noCells[1]) * m_bc2600noCells[2];

          for(MInt var = 0; var < m_maxNoVariables; var++) {
            m_bc2600Variables[var][cellIdG1] = m_bc2600Variables[var][cellIdA1];
            m_bc2600Variables[var][cellIdG2] = m_bc2600Variables[var][cellIdA1];
          }

          cellIdA1 = i + (j + (m_bc2600noCells[0] - 3) * m_bc2600noCells[1]) * m_bc2600noCells[2];
          cellIdG1 = i + (j + (m_bc2600noCells[0] - 2) * m_bc2600noCells[1]) * m_bc2600noCells[2];
          cellIdG2 = i + (j + (m_bc2600noCells[0] - 1) * m_bc2600noCells[1]) * m_bc2600noCells[2];

          for(MInt var = 0; var < m_maxNoVariables; var++) {
            m_bc2600Variables[var][cellIdG1] = m_bc2600Variables[var][cellIdA1];
            m_bc2600Variables[var][cellIdG2] = m_bc2600Variables[var][cellIdA1];
          }
        }
      }


      // Fix diagonal cells at end of domain
      if(m_bc2600noOffsetCells[1] + m_bc2600noActiveCells[1] == m_grid->getMyBlockNoCells(1)) {
        for(MInt i = 0; i < m_bc2600noCells[2]; i++) {
          for(MInt k = 0; k < m_bc2600noCells[0]; k++) {
            const MInt cellIdA2 =
                i + (m_noGhostLayers + m_bc2600noActiveCells[1] - 2 + k * m_bc2600noCells[1]) * m_bc2600noCells[2];
            const MInt cellIdA1 =
                i + (m_noGhostLayers + m_bc2600noActiveCells[1] - 1 + k * m_bc2600noCells[1]) * m_bc2600noCells[2];
            const MInt cellIdG1 =
                i + (m_noGhostLayers + m_bc2600noActiveCells[1] + k * m_bc2600noCells[1]) * m_bc2600noCells[2];
            for(MInt var = 0; var < m_maxNoVariables; var++) {
              const MFloat distA1A2 =
                  sqrt(POW2(m_cells->coordinates[0][cellIdA1] - m_cells->coordinates[0][cellIdA2])
                       + POW2(m_cells->coordinates[1][cellIdA1] - m_cells->coordinates[1][cellIdA2])
                       + POW2(m_cells->coordinates[2][cellIdA1] - m_cells->coordinates[2][cellIdA2]));

              const MFloat slope = (m_bc2600Variables[var][cellIdA1] - m_bc2600Variables[var][cellIdA2]) / distA1A2;
              const MFloat distG1A1 =
                  sqrt(POW2(m_cells->coordinates[0][cellIdG1] - m_cells->coordinates[0][cellIdA1])
                       + POW2(m_cells->coordinates[1][cellIdG1] - m_cells->coordinates[1][cellIdA1])
                       + POW2(m_cells->coordinates[2][cellIdG1] - m_cells->coordinates[2][cellIdA1]));
              m_bc2600Variables[var][cellIdG1] = m_bc2600Variables[var][cellIdA1] + distG1A1 * slope;
            }
          }
        }
      }
    }

    if(m_commBC2600MyRank == 0) {
      cout << "Loading BC2600 values... SUCCESSFUL!" << endl;
    }
  }
}


void FvStructuredSolver3D::loadRestartBC2601() {
  // if we have the prescribing boundary,
  // put the values from the restart file into the ghostcells
  if(m_bc2601IsActive && !m_bc2601InitialStartup) {
    if(domainId() == 0) {
      cout << "Loading restart values 2601" << endl;
    }
    ParallelIo::size_type bcCells[3] = {0, 0, 0};
    ParallelIo::size_type bcOffset[3] = {0, 0, 0};

    bcCells[0] = m_grid->getMyBlockNoCells(0);
    bcCells[1] = m_noGhostLayers;
    bcCells[2] = m_grid->getMyBlockNoCells(2);

    MInt noCellsBC = bcCells[0] * bcCells[1] * bcCells[2];
    MFloatScratchSpace tmpRestartVars(noCellsBC * PV->noVariables, AT_, "m_tmpRestartVars2600");

    if(domainId() == 0) {
      stringstream restartFileName;
      MString restartFile = Context::getSolverProperty<MString>("restartVariablesFileName", m_solverId, AT_);
      restartFileName << outputDir() << restartFile;

      ParallelIoHdf5 pio(restartFileName.str(), maia::parallel_io::PIO_READ, MPI_COMM_SELF);
      stringstream pathStr;
      pathStr << "/block" << m_blockId << "/bc2601" << endl;
      std::string path = pathStr.str();

      for(MInt var = 0; var < PV->noVariables; var++) {
        pio.readArray(&tmpRestartVars[var * noCellsBC], path, m_pvariableNames[var].c_str(), nDim, bcOffset, bcCells);
      }
    }

    MPI_Bcast(&tmpRestartVars[0], noCellsBC * PV->noVariables, MPI_DOUBLE, 0, m_StructuredComm, AT_,
              "tmpRestartVars[0]");

    if(m_bc2601) {
      MInt startGC[3] = {0, 0, 0};
      MInt endGC[3] = {0, 0, 0};

      if(m_nOffsetCells[2] == 0) {
        startGC[2] = m_noGhostLayers;
      }
      if(m_nOffsetCells[0] == 0) {
        startGC[0] = m_noGhostLayers;
      }
      if(m_nOffsetCells[2] + m_nActiveCells[2] == m_grid->getMyBlockNoCells(2)) {
        endGC[2] = m_noGhostLayers;
      }
      if(m_nOffsetCells[0] + m_nActiveCells[0] == m_grid->getMyBlockNoCells(0)) {
        endGC[0] = m_noGhostLayers;
      }

      for(MInt i = startGC[2]; i < m_nCells[2] - endGC[2]; i++) {
        for(MInt j = 0; j < m_noGhostLayers; j++) {
          for(MInt k = startGC[0]; k < m_nCells[0] - endGC[0]; k++) {
            MInt cellId = cellIndex(i, j, k);
            MInt globalI = m_nOffsetCells[2] - m_noGhostLayers + i;
            MInt globalJ = j;
            MInt globalK = m_nOffsetCells[0] - m_noGhostLayers + k;
            MInt cellIdBC = globalI + (globalJ + globalK * m_noGhostLayers) * m_grid->getMyBlockNoCells(2);

            // load values from restart field
            for(MInt var = 0; var < PV->noVariables; var++) {
              m_cells->pvariables[var][cellId] = tmpRestartVars[var * noCellsBC + cellIdBC];
            }
          }
        }
      }

      // Fix diagonal cells at start of domain
      if(m_nOffsetCells[2] == 0) {
        for(MInt j = 0; j < m_noGhostLayers; j++) {
          for(MInt k = 0; k < m_nCells[0]; k++) {
            const MInt cellIdA2 = cellIndex(3, j, k);
            const MInt cellIdA1 = cellIndex(2, j, k);
            const MInt cellIdG1 = cellIndex(1, j, k);
            for(MInt var = 0; var < PV->noVariables; var++) {
              const MFloat slope = (m_cells->pvariables[var][cellIdA2] - m_cells->pvariables[var][cellIdA1])
                                   / (m_cells->coordinates[0][cellIdA2] - m_cells->coordinates[0][cellIdA1]);
              m_cells->pvariables[var][cellIdG1] =
                  m_cells->pvariables[var][cellIdA1]
                  + (m_cells->coordinates[0][cellIdG1] - m_cells->coordinates[0][cellIdA1]) * slope;
            }
          }
        }
      }

      // Fix diagonal cells at end of domain
      if(m_nOffsetCells[2] + m_nActiveCells[2] == m_grid->getMyBlockNoCells(2)) {
        for(MInt j = 0; j < m_noGhostLayers; j++) {
          for(MInt k = 0; k < m_nCells[0]; k++) {
            const MInt cellIdA2 = cellIndex(m_noGhostLayers + m_nActiveCells[2] - 2, j, k);
            const MInt cellIdA1 = cellIndex(m_noGhostLayers + m_nActiveCells[2] - 1, j, k);
            const MInt cellIdG1 = cellIndex(m_noGhostLayers + m_nActiveCells[2], j, k);
            for(MInt var = 0; var < PV->noVariables; var++) {
              const MFloat slope = (m_cells->pvariables[var][cellIdA1] - m_cells->pvariables[var][cellIdA2])
                                   / (m_cells->coordinates[0][cellIdA1] - m_cells->coordinates[0][cellIdA2]);
              m_cells->pvariables[var][cellIdG1] =
                  m_cells->pvariables[var][cellIdA1]
                  + (m_cells->coordinates[0][cellIdG1] - m_cells->coordinates[0][cellIdA1]) * slope;
            }
          }
        }
      }
    }
  }
}

void FvStructuredSolver3D::loadRestartSTG(MBool isPrimitiveOutput) {
  if(m_stgIsActive) {
    if(!isPrimitiveOutput && domainId() == 0) {
      cout << "Restart file has conservative variables, converting STG variables to primitive!" << endl;
    }

    stringstream restartFileName;
    MString restartFile = Context::getSolverProperty<MString>("restartVariablesFileName", m_solverId, AT_);
    restartFileName << outputDir() << restartFile;
    ParallelIoHdf5 pio(restartFileName.str(), maia::parallel_io::PIO_READ, m_StructuredComm);

    stringstream blockNumber;
    blockNumber << m_blockId;
    MString blockPathStr = "/block";
    blockPathStr += blockNumber.str();

    MInt restartNoEddies = 0;
    pio.getAttribute(&restartNoEddies, "stgNRAN", "");

    if(restartNoEddies != m_stgMaxNoEddies) {
      m_log << "STG: NRAN in restart file (" << restartNoEddies << ") not the same as given in property file ("
            << m_stgMaxNoEddies << "), creating new random distribution of eddies!" << endl;
      m_stgCreateNewEddies = true;
    } else {
      m_log << "STG: Reading in " << restartNoEddies << " eddies from restart file" << endl;
    }

    // if this is an initialStartup the new
    // eddies need to be created in any case
    if(m_stgInitialStartup) {
      m_stgCreateNewEddies = true;
      ParallelIo::size_type offset[3] = {m_nOffsetCells[0], m_nOffsetCells[1], m_nOffsetCells[2]};
      ParallelIo::size_type size[3] = {m_nActiveCells[0], m_nActiveCells[1], m_nActiveCells[2]};
      // also load nu_t into fq field
      pio.readArray(m_cells->fq[FQ->NU_T], blockPathStr, FQ->fqNames[FQ->NU_T].c_str(), nDim, offset, size);
      FQ->loadedFromRestartFile[FQ->NU_T] = true;
    } else {
      // has to be set manually for restart from RANS profile
      ParallelIo::size_type ninmax = m_stgNoEddieProperties * m_stgMaxNoEddies;
      ParallelIo::size_type VBStart = 0;


      //////////////////////////////////////////////////
      ////////////// LOAD EDDIES ///////////////////////
      //////////////////////////////////////////////////
      if(m_stgCreateNewEddies) {
        if(domainId() == 0) {
          cout << "NRAN in property differs from NRAN in restart file"
               << " NOT READING EDDIES FROM RESTART!" << endl;
        }
      } else {
        MString stgGlobalPathStr = "stgGlobal";

        if(pio.hasDataset("FQeddies", stgGlobalPathStr)) {
          if(domainId() == 0) {
            cout << "FQeddies field is at new position /stgGlobal/FQeddies" << endl;
          }
        } else {
          if(domainId() == 0) {
            cout << "FQeddies field is NOT at new position! Using old path within solver..." << endl;
          }
          stgGlobalPathStr = blockPathStr;
          if(domainId() == 0) {
            cout << "Loading FQeddies from path " << stgGlobalPathStr << endl;
          }
        }

        // do this only serially
        if(globalDomainId() == 0) {
          cout << "Loading STG Eddies..." << endl;
        }
        if(globalDomainId() == 0) {
          ParallelIoHdf5 pioLocal(restartFileName.str(), maia::parallel_io::PIO_READ, MPI_COMM_SELF);
          pioLocal.readArray(m_stgEddies[0], stgGlobalPathStr, "FQeddies", 1, &VBStart, &ninmax);
        }

        MPI_Bcast(m_stgEddies[0], ninmax, MPI_DOUBLE, 0, m_StructuredComm, AT_, "m_stgEddies[0]");
        if(globalDomainId() == 0) {
          cout << "Loading STG Eddies... SUCCESSFUL!" << endl;
        }
      }


      //////////////////////////////////////////////////
      ////////////// LOAD STG VARIABLES ////////////////
      //////////////////////////////////////////////////

      if(m_stgLocal) {
        ParallelIo::size_type bcActiveCells[3] = {m_grid->getMyBlockNoCells(0), m_grid->getMyBlockNoCells(1), 3};
        ParallelIo::size_type bcCells[3] = {m_grid->getMyBlockNoCells(0) + m_noGhostLayers * 2,
                                            m_grid->getMyBlockNoCells(1) + m_noGhostLayers * 2, 3};
        ParallelIo::size_type bcOffset[3] = {0, 0, 0};
        const MInt noCellsBC = bcCells[0] * bcCells[1] * bcCells[2];
        MFloatScratchSpace tmpRestartVars(noCellsBC * m_stgNoVariables, AT_, "tmpRestartVars");

        if(m_commStgMyRank == 0) {
          cout << "Loading STG Datasets..." << endl;
        }
        if(m_commStgMyRank == 0) {
          ParallelIoHdf5 pioLocal(restartFileName.str(), maia::parallel_io::PIO_READ, MPI_COMM_SELF);
          cout << "stg_myRankdomainId:" << domainId() << " m_commStgMyRank:" << m_commStgMyRank << endl;
          for(MInt var = 0; var < m_stgNoVariables; var++) {
            stringstream fieldName;
            stringstream stgPath;
            stgPath << blockPathStr << "/stg";
            fieldName << "stgFQ" << var;
            pioLocal.readArray(&tmpRestartVars[var * noCellsBC], stgPath.str(), fieldName.str(), nDim, bcOffset,
                               bcActiveCells);
          }

          // reorder the cells in the global array
          for(MInt k = (bcActiveCells[0] - 1); k >= 0; k--) {
            for(MInt j = (bcActiveCells[1] - 1); j >= 0; j--) {
              for(MInt i = (bcActiveCells[2] - 1); i >= 0; i--) {
                const MInt cellId_org = i + (j + k * bcActiveCells[1]) * bcActiveCells[2];
                const MInt i_new = i;
                const MInt j_new = j + m_noGhostLayers;
                const MInt k_new = k + m_noGhostLayers;
                const MInt cellId = i_new + (j_new + k_new * bcCells[1]) * bcCells[2];

                for(MInt var = 0; var < m_stgNoVariables; var++) {
                  tmpRestartVars[var * noCellsBC + cellId] = tmpRestartVars[var * noCellsBC + cellId_org];
                  tmpRestartVars[var * noCellsBC + cellId_org] = F0;
                }
              }
            }
          }

          // apply periodic bc
          for(MInt j = 0; j < bcCells[1]; j++) {
            for(MInt i = 0; i < bcCells[2]; i++) {
              const MInt gcId0 = i + (j + (0) * bcCells[1]) * bcCells[2];
              const MInt acId0 = i + (j + (bcCells[0] - 4) * bcCells[1]) * bcCells[2];
              const MInt gcId1 = i + (j + (1) * bcCells[1]) * bcCells[2];
              const MInt acId1 = i + (j + (bcCells[0] - 3) * bcCells[1]) * bcCells[2];

              const MInt gcId2 = i + (j + (bcCells[0] - 2) * bcCells[1]) * bcCells[2];
              const MInt acId2 = i + (j + (2) * bcCells[1]) * bcCells[2];
              const MInt gcId3 = i + (j + (bcCells[0] - 1) * bcCells[1]) * bcCells[2];
              const MInt acId3 = i + (j + (3) * bcCells[1]) * bcCells[2];

              for(MInt var = 0; var < m_stgNoVariables; var++) {
                tmpRestartVars[var * noCellsBC + gcId0] = tmpRestartVars[var * noCellsBC + acId0];
                tmpRestartVars[var * noCellsBC + gcId1] = tmpRestartVars[var * noCellsBC + acId1];
                tmpRestartVars[var * noCellsBC + gcId2] = tmpRestartVars[var * noCellsBC + acId2];
                tmpRestartVars[var * noCellsBC + gcId3] = tmpRestartVars[var * noCellsBC + acId3];
              }
            }
          }

          // extrapolation at the bottoms
          for(MInt k = 0; k < bcCells[0]; k++) {
            for(MInt i = 0; i < bcCells[2]; i++) {
              const MInt gc1 = i + (1 + k * bcCells[1]) * bcCells[2];
              const MInt gc2 = i + (0 + k * bcCells[1]) * bcCells[2];
              const MInt ac1 = i + (2 + k * bcCells[1]) * bcCells[2];
              const MInt ac2 = i + (3 + k * bcCells[1]) * bcCells[2];

              for(MInt var = 0; var < m_stgNoVariables; var++) {
                tmpRestartVars[var * noCellsBC + gc1] =
                    2 * tmpRestartVars[var * noCellsBC + ac1] - tmpRestartVars[var * noCellsBC + ac2];
                tmpRestartVars[var * noCellsBC + gc2] =
                    2 * tmpRestartVars[var * noCellsBC + gc1] - tmpRestartVars[var * noCellsBC + ac1];
              }
            }
          }

          // extrapolation at the top
          for(MInt k = 0; k < bcCells[0]; k++) {
            for(MInt i = 0; i < bcCells[2]; i++) {
              const MInt gc1 = i + (bcCells[1] - 2 + k * bcCells[1]) * bcCells[2];
              const MInt gc2 = i + (bcCells[1] - 1 + k * bcCells[1]) * bcCells[2];
              const MInt ac1 = i + (bcCells[1] - 3 + k * bcCells[1]) * bcCells[2];
              const MInt ac2 = i + (bcCells[1] - 4 + k * bcCells[1]) * bcCells[2];


              for(MInt var = 0; var < m_stgNoVariables; var++) {
                tmpRestartVars[var * noCellsBC + gc1] =
                    2 * tmpRestartVars[var * noCellsBC + ac1] - tmpRestartVars[var * noCellsBC + ac2];
                tmpRestartVars[var * noCellsBC + gc2] =
                    2 * tmpRestartVars[var * noCellsBC + gc1] - tmpRestartVars[var * noCellsBC + ac1];
              }
            }
          }
        }


        // now broadcast to everyone
        MPI_Bcast(&tmpRestartVars[0], noCellsBC * m_stgNoVariables, MPI_DOUBLE, 0, m_commStg, AT_, "tmpRestartVars[0]");
        if(m_commStgMyRank == 0) {
          cout << "Loading STG Datasets... SUCCESSFUL!" << endl;
        }

        //////////////////////////////////////////////////
        ////////// DISTRIBUTE STG VARIABLES //////////////
        //////////////////////////////////////////////////
        for(MInt k = 0; k < m_nCells[0]; k++) {
          for(MInt j = 0; j < m_nCells[1]; j++) {
            for(MInt i = 0; i < 3; i++) {
              MInt cellId = cellIndex(i, j, k);
              MInt globalI = i;
              MInt globalJ = m_nOffsetCells[1] + j;
              MInt globalK = m_nOffsetCells[0] + k;
              MInt cellIdBCGlobal = globalI + (globalJ + globalK * bcCells[1]) * bcCells[2];
              MInt cellIdBC = i + (j + k * m_nCells[1]) * 3;

              // load values from restart field
              for(MInt var = 0; var < m_stgNoVariables; var++) {
                m_cells->stg_fq[var][cellIdBC] = tmpRestartVars[var * noCellsBC + cellIdBCGlobal];
              }

              if(!isPrimitiveOutput) {
                const MFloat rho = m_cells->stg_fq[0][cellIdBC];
                const MFloat rhoU = m_cells->stg_fq[1][cellIdBC];
                const MFloat rhoV = m_cells->stg_fq[2][cellIdBC];
                const MFloat rhoW = m_cells->stg_fq[3][cellIdBC];
                const MFloat rhoE = m_cells->stg_fq[4][cellIdBC];

                const MFloat gammaMinusOne = m_gamma - 1.0;
                const MFloat u = rhoU / rho;
                const MFloat v = rhoV / rho;
                const MFloat w = rhoW / rho;
                const MFloat p = gammaMinusOne * (rhoE - F1B2 * rho * (POW2(u) + POW2(v) + POW2(w)));

                m_cells->stg_fq[PV->RHO][cellIdBC] = rho;
                m_cells->stg_fq[PV->U][cellIdBC] = u;
                m_cells->stg_fq[PV->V][cellIdBC] = v;
                m_cells->stg_fq[PV->W][cellIdBC] = w;
                m_cells->stg_fq[PV->P][cellIdBC] = p;
              }

              if(i < 2) {
                m_cells->pvariables[PV->RHO][cellId] = m_cells->stg_fq[PV->RHO][cellIdBC];
                m_cells->pvariables[PV->U][cellId] = m_cells->stg_fq[PV->U][cellIdBC];
                m_cells->pvariables[PV->V][cellId] = m_cells->stg_fq[PV->V][cellIdBC];
                m_cells->pvariables[PV->W][cellId] = m_cells->stg_fq[PV->W][cellIdBC];
                m_cells->pvariables[PV->P][cellId] = m_cells->stg_fq[PV->P][cellIdBC];
              }
            }
          }
        }
      }
    }

  }
}


void FvStructuredSolver3D::computeVorticity() {
  TRACE();
  MFloat* const RESTRICT u = &m_cells->pvariables[PV->U][0];
  MFloat* const RESTRICT v = &m_cells->pvariables[PV->V][0];
  MFloat* const RESTRICT w = &m_cells->pvariables[PV->W][0];
  MFloat* const RESTRICT vortx = &m_cells->fq[FQ->VORTX][0];
  MFloat* const RESTRICT vorty = &m_cells->fq[FQ->VORTY][0];
  MFloat* const RESTRICT vortz = &m_cells->fq[FQ->VORTZ][0];
  MFloat* const RESTRICT jac = &m_cells->cellJac[0];


  for(MInt k = m_noGhostLayers; k < m_nCells[0] - m_noGhostLayers; k++) {
    for(MInt j = m_noGhostLayers; j < m_nCells[1] - m_noGhostLayers; j++) {
      for(MInt i = m_noGhostLayers; i < m_nCells[2] - m_noGhostLayers; i++) {
        const MInt IJK = i + (j + k * m_nCells[1]) * m_nCells[2];
        const MInt IPJK = (i + 1) + (j + k * m_nCells[1]) * m_nCells[2];
        const MInt IMJK = (i - 1) + (j + k * m_nCells[1]) * m_nCells[2];
        const MInt IJPK = i + ((j + 1) + k * m_nCells[1]) * m_nCells[2];
        const MInt IJMK = i + ((j - 1) + k * m_nCells[1]) * m_nCells[2];
        const MInt IJKP = i + (j + (k + 1) * m_nCells[1]) * m_nCells[2];
        const MInt IJKM = i + (j + (k - 1) * m_nCells[1]) * m_nCells[2];

        const MFloat dudxi = u[IPJK] - u[IMJK];
        const MFloat dudet = u[IJPK] - u[IJMK];
        const MFloat dudze = u[IJKP] - u[IJKM];

        const MFloat dvdxi = v[IPJK] - v[IMJK];
        const MFloat dvdet = v[IJPK] - v[IJMK];
        const MFloat dvdze = v[IJKP] - v[IJKM];

        const MFloat dwdxi = w[IPJK] - w[IMJK];
        const MFloat dwdet = w[IJPK] - w[IJMK];
        const MFloat dwdze = w[IJKP] - w[IJKM];

        const MFloat dvdz = dvdxi * m_cells->cellMetrics[xsd * 3 + zsd][IJK]
                            + dvdet * m_cells->cellMetrics[ysd * 3 + zsd][IJK]
                            + dvdze * m_cells->cellMetrics[zsd * 3 + zsd][IJK];
        const MFloat dwdy = dwdxi * m_cells->cellMetrics[xsd * 3 + ysd][IJK]
                            + dwdet * m_cells->cellMetrics[ysd * 3 + ysd][IJK]
                            + dwdze * m_cells->cellMetrics[zsd * 3 + ysd][IJK];

        const MFloat dudz = dudxi * m_cells->cellMetrics[xsd * 3 + zsd][IJK]
                            + dudet * m_cells->cellMetrics[ysd * 3 + zsd][IJK]
                            + dudze * m_cells->cellMetrics[zsd * 3 + zsd][IJK];
        const MFloat dwdx = dwdxi * m_cells->cellMetrics[xsd * 3 + xsd][IJK]
                            + dwdet * m_cells->cellMetrics[ysd * 3 + xsd][IJK]
                            + dwdze * m_cells->cellMetrics[zsd * 3 + xsd][IJK];

        const MFloat dvdx = dvdxi * m_cells->cellMetrics[xsd * 3 + xsd][IJK]
                            + dvdet * m_cells->cellMetrics[ysd * 3 + xsd][IJK]
                            + dvdze * m_cells->cellMetrics[zsd * 3 + xsd][IJK];
        const MFloat dudy = dudxi * m_cells->cellMetrics[xsd * 3 + ysd][IJK]
                            + dudet * m_cells->cellMetrics[ysd * 3 + ysd][IJK]
                            + dudze * m_cells->cellMetrics[zsd * 3 + ysd][IJK];

        vortx[IJK] = F1B2 * (dwdy - dvdz) / jac[IJK];
        vorty[IJK] = F1B2 * (dudz - dwdx) / jac[IJK];
        vortz[IJK] = F1B2 * (dvdx - dudy) / jac[IJK];
      }
    }
  }
}


/**
 * \brief Function to compute the lambda_2 criterion
 * \author Pascal Meysonnat
 *
 */
void FvStructuredSolver3D::computeLambda2Criterion() {
  TRACE();
  MFloatScratchSpace J(nDim, nDim, AT_, "J");
  MFloat d[3] = {F0, F0, F0};
  MFloat e[3] = {F0, F0, F0};
  J.fill(F0);

  // MInt IMJK, IJMK, IMJMK, IJKM, IJMKM, IMJKM;
  MFloat* const RESTRICT u = &m_cells->pvariables[PV->U][0];
  MFloat* const RESTRICT v = &m_cells->pvariables[PV->V][0];
  MFloat* const RESTRICT w = &m_cells->pvariables[PV->W][0];
  MFloat* const RESTRICT lambda2 = &m_cells->fq[FQ->LAMBDA2][0];

  // compute the lambda2 criterion.
  for(MInt k = m_noGhostLayers - 1; k < m_nCells[0] - m_noGhostLayers + 1; k++) {
    for(MInt j = m_noGhostLayers - 1; j < m_nCells[1] - m_noGhostLayers + 1; j++) {
      for(MInt i = m_noGhostLayers - 1; i < m_nCells[2] - m_noGhostLayers + 1; i++) {
        const MInt cellId = cellIndex(i, j, k);
        const MInt IMJK = cellIndex(i - 1, j, k);
        const MInt IPJK = cellIndex(i + 1, j, k);
        const MInt IJMK = cellIndex(i, j - 1, k);
        const MInt IJPK = cellIndex(i, j + 1, k);
        const MInt IJKM = cellIndex(i, j, k - 1);
        const MInt IJKP = cellIndex(i, j, k + 1);
        const MFloat FcellJac = F1 / m_cells->cellJac[cellId];

        const MFloat dudx = FcellJac
                            * (m_cells->cellMetrics[0][cellId] * (u[IPJK] - u[IMJK])
                               + m_cells->cellMetrics[3][cellId] * (u[IJPK] - u[IJMK])
                               + m_cells->cellMetrics[6][cellId] * (u[IJKP] - u[IJKM]));
        const MFloat dudy = FcellJac
                            * (m_cells->cellMetrics[1][cellId] * (u[IPJK] - u[IMJK])
                               + m_cells->cellMetrics[4][cellId] * (u[IJPK] - u[IJMK])
                               + m_cells->cellMetrics[7][cellId] * (u[IJKP] - u[IJKM]));
        const MFloat dudz = FcellJac
                            * (m_cells->cellMetrics[2][cellId] * (u[IPJK] - u[IMJK])
                               + m_cells->cellMetrics[5][cellId] * (u[IJPK] - u[IJMK])
                               + m_cells->cellMetrics[8][cellId] * (u[IJKP] - u[IJKM]));

        const MFloat dvdx = FcellJac
                            * (m_cells->cellMetrics[0][cellId] * (v[IPJK] - v[IMJK])
                               + m_cells->cellMetrics[3][cellId] * (v[IJPK] - v[IJMK])
                               + m_cells->cellMetrics[6][cellId] * (v[IJKP] - v[IJKM]));
        const MFloat dvdy = FcellJac
                            * (m_cells->cellMetrics[1][cellId] * (v[IPJK] - v[IMJK])
                               + m_cells->cellMetrics[4][cellId] * (v[IJPK] - v[IJMK])
                               + m_cells->cellMetrics[7][cellId] * (v[IJKP] - v[IJKM]));
        const MFloat dvdz = FcellJac
                            * (m_cells->cellMetrics[2][cellId] * (v[IPJK] - v[IMJK])
                               + m_cells->cellMetrics[5][cellId] * (v[IJPK] - v[IJMK])
                               + m_cells->cellMetrics[8][cellId] * (v[IJKP] - v[IJKM]));

        const MFloat dwdx = FcellJac
                            * (m_cells->cellMetrics[0][cellId] * (w[IPJK] - w[IMJK])
                               + m_cells->cellMetrics[3][cellId] * (w[IJPK] - w[IJMK])
                               + m_cells->cellMetrics[6][cellId] * (w[IJKP] - w[IJKM]));
        const MFloat dwdy = FcellJac
                            * (m_cells->cellMetrics[1][cellId] * (w[IPJK] - w[IMJK])
                               + m_cells->cellMetrics[4][cellId] * (w[IJPK] - w[IJMK])
                               + m_cells->cellMetrics[7][cellId] * (w[IJKP] - w[IJKM]));
        const MFloat dwdz = FcellJac
                            * (m_cells->cellMetrics[2][cellId] * (w[IPJK] - w[IMJK])
                               + m_cells->cellMetrics[5][cellId] * (w[IJPK] - w[IJMK])
                               + m_cells->cellMetrics[8][cellId] * (w[IJKP] - w[IJKM]));

        // Compute the matrix (S^2+Omega^2)/2
        J(0, 0) = POW2(dudx) + dvdx * dudy + dwdx * dudz;
        J(0, 1) = F1B2 * (dudx * dudy + dudx * dvdx + dvdy * dudy + dvdy * dvdx + dwdy * dudz + dvdz * dwdx);
        J(0, 2) = F1B2 * (dudx * dudz + dudx * dwdx + dudy * dvdz + dvdx * dwdy + dwdz * dudz + dwdz * dwdx);
        J(1, 0) = J(0, 1);
        J(1, 1) = POW2(dvdy) + dvdx * dudy + dvdz * dwdy;
        J(1, 2) = F1B2 * (dudz * dvdx + dwdx * dudy + dvdy * dvdz + dvdy * dwdy + dwdz * dvdz + dwdz * dwdy);
        J(2, 0) = J(0, 2);
        J(2, 1) = J(1, 2);
        J(2, 2) = POW2(dwdz) + dwdx * dudz + dvdz * dwdy;

        // perform householder tridiagonalization of
        // symmetric real matrix
        tred2(J, nDim, d, e);
        // compute eigenvalues
        tqli2(d, e, nDim);
        // sort eigenvalues
        insertSort(nDim, d);

        lambda2[cellId] = d[1];
      }
    }
  }
}


void FvStructuredSolver3D::saveInterpolatedPoints() {
  TRACE();
  ////////////////////////////////
  ////////// INTERPOLATION ///////
  ////////////////////////////////

  // if it a moving grid the interpolation
  // coefficients need to be computed every time
  if(m_movingGrid) {
    for(MInt pointId = 0; pointId < m_intpPointsNoPointsTotal; pointId++) {
      m_intpPointsHasPartnerLocal[pointId] = 0;
      m_intpPointsHasPartnerGlobal[pointId] = 0;
      for(MInt var = 0; var < PV->noVariables; var++) {
        m_intpPointsVarsLocal[var][pointId] = F0;
        m_intpPointsVarsGlobal[var][pointId] = F0;
      }
    }

    m_pointInterpolation =
        make_unique<StructuredInterpolation<3>>(m_nCells, m_cells->coordinates, m_cells->pvariables, m_StructuredComm);
    // allocate domains points of the lines
    m_pointInterpolation->prepareInterpolation(m_intpPointsNoPointsTotal, m_intpPointsCoordinates,
                                               m_intpPointsHasPartnerLocal);
    MPI_Allreduce(m_intpPointsHasPartnerLocal, m_intpPointsHasPartnerGlobal, m_intpPointsNoPointsTotal, MPI_INT,
                  MPI_SUM, m_StructuredComm, AT_, "m_intpPointsHasPartnerLocal", "m_intpPointsHasPartnerGlobal");
  }

  // calculation of interpolated variables only for the points in domain
  for(MInt pointId = 0; pointId < m_intpPointsNoPointsTotal; pointId++) {
    if(m_intpPointsHasPartnerLocal[pointId]) {
      for(MInt var = 0; var < PV->noVariables; var++) {
        m_intpPointsVarsLocal[var][pointId] = m_pointInterpolation->getInterpolatedVariable(pointId, var);
      }
    }
  }

  MPI_Allreduce(&m_intpPointsVarsLocal[0][0], &m_intpPointsVarsGlobal[0][0],
                m_intpPointsNoPointsTotal * PV->noVariables, MPI_DOUBLE, MPI_SUM, m_StructuredComm, AT_,
                "m_intpPointsVarsLocal[0][0]", "m_intpPointsVarsGlobal[0][0]");

  // calculation of right value, if a point is assigned to more than one domain
  for(MInt pointId = 0; pointId < m_intpPointsNoPointsTotal; pointId++) {
    if(m_intpPointsHasPartnerGlobal[pointId] > 1) {
      for(MInt var = 0; var < PV->noVariables; var++) {
        m_intpPointsVarsGlobal[var][pointId] =
            m_intpPointsVarsGlobal[var][pointId] / (MFloat)m_intpPointsHasPartnerGlobal[pointId];
      }
    }
  }

  if(m_movingGrid) {
    m_pointInterpolation.reset();
  }

  ////////////////////////////////
  ////////// WRITE HDF5 //////////
  ////////////////////////////////
  stringstream fileName;
  fileName << m_intpPointsOutputDir << "interpolatedPoints" << globalTimeStep << m_outputFormat;

  ParallelIoHdf5 pio(fileName.str(), maia::parallel_io::PIO_REPLACE, m_StructuredComm);

  writeHeaderAttributes(&pio, "field");
  writePropertiesAsAttributes(&pio, "");

  pio.setAttribute(m_intpPointsNoLines, "noFields", "");

  for(MInt lineId = 0; lineId < m_intpPointsNoLines; lineId++) {
    ParallelIo::size_type dataOffset[2] = {0, 0};
    ParallelIo::size_type dataSize[2] = {m_intpPointsNoPoints2D[lineId], m_intpPointsNoPoints[lineId]};
    MInt fieldOffset = m_intpPointsOffsets[lineId];

    ParallelIo::size_type noDims = 1;
    dataSize[0] = dataSize[1];
    if(m_intpPoints) {
      noDims = 2;
    }

    stringstream path;
    path << lineId;
    MString solutionpath = "field";
    solutionpath += path.str();
    const char* dsetname = solutionpath.c_str();

    for(MInt v = 0; v < PV->noVariables; v++) {
      pio.defineArray(maia::parallel_io::PIO_FLOAT, dsetname, m_pvariableNames[v], noDims, dataSize);
    }

    pio.defineArray(maia::parallel_io::PIO_FLOAT, dsetname, "x", noDims, dataSize);
    pio.defineArray(maia::parallel_io::PIO_FLOAT, dsetname, "y", noDims, dataSize);
    pio.defineArray(maia::parallel_io::PIO_FLOAT, dsetname, "z", noDims, dataSize);

    if(domainId() == 0) {
      for(MInt v = 0; v < PV->noVariables; v++) {
        pio.writeArray(&m_intpPointsVarsGlobal[v][fieldOffset], dsetname, m_pvariableNames[v], noDims, dataOffset,
                       dataSize);
      }

      pio.writeArray(&m_intpPointsCoordinates[0][fieldOffset], dsetname, "x", noDims, dataOffset, dataSize);
      pio.writeArray(&m_intpPointsCoordinates[1][fieldOffset], dsetname, "y", noDims, dataOffset, dataSize);
      pio.writeArray(&m_intpPointsCoordinates[2][fieldOffset], dsetname, "z", noDims, dataOffset, dataSize);
    } else {
      dataSize[0] = 0;
      dataSize[1] = 0;
      MFloat empty = 0;
      for(MInt v = 0; v < PV->noVariables; v++) {
        pio.writeArray(&empty, dsetname, m_pvariableNames[v], 1, dataOffset, dataSize);
      }

      pio.writeArray(&empty, dsetname, "x", noDims, dataOffset, dataSize);
      pio.writeArray(&empty, dsetname, "y", noDims, dataOffset, dataSize);
      pio.writeArray(&empty, dsetname, "z", noDims, dataOffset, dataSize);
    }
  }
}

void FvStructuredSolver3D::applyInviscidBoundaryCondition() { TRACE(); }

void FvStructuredSolver3D::applyViscousBoundaryCondition() { TRACE(); }

void FvStructuredSolver3D::getSampleVariables(MInt cellId, MFloat* cellVars) {
  cellVars[PV->U] = m_cells->pvariables[PV->U][cellId];
  cellVars[PV->V] = m_cells->pvariables[PV->V][cellId];
  cellVars[PV->W] = m_cells->pvariables[PV->W][cellId];
  cellVars[PV->RHO] = m_cells->pvariables[PV->RHO][cellId];
  cellVars[PV->P] = m_cells->pvariables[PV->P][cellId];
}

MFloat FvStructuredSolver3D::getSampleVorticity(MInt cellId, MInt dim) {
  return m_cells->fq[FQ->VORTICITY[dim]][cellId];
}


/**
 * \brief Loads primitive variables from an HDF5 file
 *
 * \author Frederik Temme
 * \author 14.1.2016
 *
 */
void FvStructuredSolver3D::loadSampleFile(MString fileName) { // loading files for averaging (pre- and postsolve)
  TRACE();

  // open the file
  ParallelIoHdf5 pio(fileName, maia::parallel_io::PIO_READ, m_StructuredComm);

  stringstream blockNumber;
  blockNumber << m_blockId;
  MString blockPathStr = "/block";
  blockPathStr += blockNumber.str();
  ParallelIo::size_type offset[3] = {m_nOffsetCells[0], m_nOffsetCells[1], m_nOffsetCells[2]};
  ParallelIo::size_type size[3] = {m_nActiveCells[0], m_nActiveCells[1], m_nActiveCells[2]};

  for(MInt var = 0; var < PV->noVariables; var++) {
    pio.readArray(m_cells->pvariables[var], blockPathStr, m_pvariableNames[var], nDim, offset, size);
  }

  shiftCellValuesRestart(true);
}


MFloat FvStructuredSolver3D::dvardxyz(MInt IJK, MInt dir, MFloat* var) {
  const MInt INC[3] = {1, m_nCells[2], m_nCells[1] * m_nCells[2]};

  const MInt IPJK = IJK + INC[0];
  const MInt IMJK = IJK - INC[0];
  const MInt IJPK = IJK + INC[1];
  const MInt IJMK = IJK - INC[1];
  const MInt IJKP = IJK + INC[2];
  const MInt IJKM = IJK - INC[2];

  const MFloat dvardxi = F1B2 * (var[IPJK] - var[IMJK]);
  const MFloat dvardet = F1B2 * (var[IJPK] - var[IJMK]);
  const MFloat dvardze = F1B2 * (var[IJKP] - var[IJKM]);

  const MFloat ddxyz =
      (dvardxi * m_cells->cellMetrics[xsd * 3 + dir][IJK] + dvardet * m_cells->cellMetrics[ysd * 3 + dir][IJK]
       + dvardze * m_cells->cellMetrics[zsd * 3 + dir][IJK]);

  return ddxyz / m_cells->cellJac[IJK];
}

MFloat FvStructuredSolver3D::dvardx(MInt IJK, MFloat* var) {
  const MInt INC[3] = {1, m_nCells[2], m_nCells[1] * m_nCells[2]};

  const MInt IPJK = IJK + INC[0];
  const MInt IMJK = IJK - INC[0];
  const MInt IJPK = IJK + INC[1];
  const MInt IJMK = IJK - INC[1];
  const MInt IJKP = IJK + INC[2];
  const MInt IJKM = IJK - INC[2];

  const MFloat dvardxi = var[IPJK] - var[IMJK];
  const MFloat dvardet = var[IJPK] - var[IJMK];
  const MFloat dvardze = var[IJKP] - var[IJKM];

  const MFloat ddx =
      (dvardxi * m_cells->cellMetrics[xsd * 3 + xsd][IJK] + dvardet * m_cells->cellMetrics[ysd * 3 + xsd][IJK]
       + dvardze * m_cells->cellMetrics[zsd * 3 + xsd][IJK]);

  return ddx / m_cells->cellJac[IJK];
}


/**
 * \brief Loads the postprocessing restart file to continue the postprocessing
 *
 * \author Frederik Temme
 * \date 06.01.2016
 */

void FvStructuredSolver3D::loadAverageRestartFile(const MChar* fileName, MFloat** sum, MFloat** square, MFloat** cube,
                                                  MFloat** fourth) {
  TRACE();
  m_log << "loading average restart file ... " << endl;
  if(domainId() == 0) {
    cout << "Opening file: " << fileName << endl;
  }

  // open the file
  MString fileNameStr = fileName;
  ParallelIoHdf5 pio(fileNameStr, maia::parallel_io::PIO_READ, m_StructuredComm);

  // now read in the data!
  m_log << "-> reading in the data ... " << endl;
  pio.getAttribute(&m_noSamples, "noSamples", "");
  m_log << "Current number of postprocessing samples: " << m_noSamples << endl;
  stringstream blockNumber;
  blockNumber << m_blockId;
  MString blockPathStr = "/block";
  blockPathStr += blockNumber.str();
  MInt offset = 0;
  ParallelIo::size_type ioOffset[3] = {m_nOffsetCells[0], m_nOffsetCells[1], m_nOffsetCells[2]};
  ParallelIo::size_type ioSize[3] = {m_nActiveCells[0], m_nActiveCells[1], m_nActiveCells[2]};


  pio.readArray(sum[0], blockPathStr, "u", nDim, ioOffset, ioSize);
  pio.readArray(sum[1], blockPathStr, "v", nDim, ioOffset, ioSize);
  pio.readArray(sum[2], blockPathStr, "w", nDim, ioOffset, ioSize);
  pio.readArray(sum[3], blockPathStr, "rho", nDim, ioOffset, ioSize);
  pio.readArray(sum[4], blockPathStr, "p", nDim, ioOffset, ioSize);
  offset = noVariables();

  if(m_averagingFavre) {
    pio.readArray(m_favre[0], blockPathStr, "um_favre", nDim, ioOffset, ioSize);
    pio.readArray(m_favre[1], blockPathStr, "vm_favre", nDim, ioOffset, ioSize);
    pio.readArray(m_favre[2], blockPathStr, "wm_favre", nDim, ioOffset, ioSize);
    pio.readArray(m_favre[3], blockPathStr, "rhom_favre", nDim, ioOffset, ioSize);
    pio.readArray(m_favre[4], blockPathStr, "pm_favre", nDim, ioOffset, ioSize);
  }

  if(m_averageVorticity) {
    pio.readArray(sum[offset + 0], blockPathStr, "vortx", nDim, ioOffset, ioSize);
    pio.readArray(sum[offset + 1], blockPathStr, "vorty", nDim, ioOffset, ioSize);
    pio.readArray(sum[offset + 2], blockPathStr, "vortz", nDim, ioOffset, ioSize);
  }

  pio.readArray(square[0], blockPathStr, "uu", nDim, ioOffset, ioSize);
  pio.readArray(square[1], blockPathStr, "vv", nDim, ioOffset, ioSize);
  pio.readArray(square[2], blockPathStr, "ww", nDim, ioOffset, ioSize);
  pio.readArray(square[3], blockPathStr, "uv", nDim, ioOffset, ioSize);
  pio.readArray(square[4], blockPathStr, "vw", nDim, ioOffset, ioSize);
  pio.readArray(square[5], blockPathStr, "uw", nDim, ioOffset, ioSize);

  pio.readArray(square[6], blockPathStr, "pp", nDim, ioOffset, ioSize);

  if(m_averageVorticity) {
    pio.readArray(square[7], blockPathStr, "vortxvortx", nDim, ioOffset, ioSize);
    pio.readArray(square[8], blockPathStr, "vortyvorty", nDim, ioOffset, ioSize);
    pio.readArray(square[9], blockPathStr, "vortzvortz", nDim, ioOffset, ioSize);
  }

  if(m_kurtosis || m_skewness) {
    pio.readArray(cube[0], blockPathStr, "uuu", nDim, ioOffset, ioSize);
    pio.readArray(cube[1], blockPathStr, "vvv", nDim, ioOffset, ioSize);
    pio.readArray(cube[2], blockPathStr, "www", nDim, ioOffset, ioSize);
  }

  if(m_kurtosis) {
    pio.readArray(fourth[0], blockPathStr, "uuuu", nDim, ioOffset, ioSize);
    pio.readArray(fourth[1], blockPathStr, "vvvv", nDim, ioOffset, ioSize);
    pio.readArray(fourth[2], blockPathStr, "wwww", nDim, ioOffset, ioSize);
  }

  m_log << "loading Restart file ... SUCCESSFUL " << endl;
  shiftAverageCellValuesRestart();
}

/**
 * \brief Loads the averaged variables again to do further postprocessing
 *
 * \author Marian Albers
 */
void FvStructuredSolver3D::loadAveragedVariables(const MChar* fileName) {
  TRACE();
  m_log << "loading averaged variables file ... " << endl;

  // open the file
  MString fileNameStr = fileName;
  ParallelIoHdf5 pio(fileNameStr, maia::parallel_io::PIO_READ, m_StructuredComm);

  stringstream blockNumber;
  blockNumber << m_blockId;
  MString blockPathStr = "/block";
  blockPathStr += blockNumber.str();
  ParallelIo::size_type ioOffset[3] = {m_nOffsetCells[0], m_nOffsetCells[1], m_nOffsetCells[2]};
  ParallelIo::size_type ioSize[3] = {m_nActiveCells[0], m_nActiveCells[1], m_nActiveCells[2]};

  for(MInt var = 0; var < getNoPPVars(); var++) {
    pio.readArray(m_summedVars[var], blockPathStr, m_avgVariableNames[var], nDim, ioOffset, ioSize);
  }

  if(m_averagingFavre) {
    for(MInt var = 0; var < getNoVars(); var++) {
      pio.readArray(m_favre[var], blockPathStr, m_avgFavreNames[var], nDim, ioOffset, ioSize);
    }
  }

  m_log << "loading Restart file ... SUCCESSFUL " << endl;
  shiftAverageCellValues();
}

/**
 *
 * \author Frederik Temme
 * \date 12.01.2015
 */
void FvStructuredSolver3D::shiftAverageCellValuesRestart() {
  TRACE();
  MInt cellId_org = 0;
  MInt cellId = 0;
  MInt i_new, j_new, k_new;

  // accounting for the ghost layers and shift the values to the right place
  for(MInt k = (m_nActiveCells[0] - 1); k >= 0; k--) {
    for(MInt j = (m_nActiveCells[1] - 1); j >= 0; j--) {
      for(MInt i = (m_nActiveCells[2] - 1); i >= 0; i--) {
        cellId_org = i + (j + k * m_nActiveCells[1]) * m_nActiveCells[2];
        i_new = i + m_noGhostLayers;
        j_new = j + m_noGhostLayers;
        k_new = k + m_noGhostLayers;
        cellId = i_new + (j_new + k_new * m_nCells[1]) * m_nCells[2];

        for(MInt var = 0; var < getNoPPVars(); var++) {
          m_summedVars[var][cellId] = m_summedVars[var][cellId_org];
          m_summedVars[var][cellId_org] = F0;
        }

        if(m_averagingFavre) {
          for(MInt var = 0; var < getNoVars(); var++) {
            m_favre[var][cellId] = m_favre[var][cellId_org];
            m_favre[var][cellId_org] = F0;
          }
        }

        for(MInt var = 0; var < getNoPPSquareVars(); var++) {
          m_square[var][cellId] = m_square[var][cellId_org];
          m_square[var][cellId_org] = F0;
        }

        if(m_kurtosis || m_skewness) {
          for(MInt var = 0; var < nDim; var++) {
            m_cube[var][cellId] = m_cube[var][cellId_org];
            m_cube[var][cellId_org] = F0;
          }
        }

        if(m_kurtosis) {
          for(MInt var = 0; var < nDim; var++) {
            m_fourth[var][cellId] = m_fourth[var][cellId_org];
            m_fourth[var][cellId_org] = F0;
          }
        }
      }
    }
  }
}

/**
 * \brief Shifts the averaged variables
 * \author Marian Albers
 * \date 2016
 */
void FvStructuredSolver3D::shiftAverageCellValues() {
  TRACE();
  MInt cellId_org = 0;
  MInt cellId = 0;
  MInt i_new, j_new, k_new;

  // accounting for the ghost layers and shift the values to the right place
  for(MInt k = (m_nActiveCells[0] - 1); k >= 0; k--) {
    for(MInt j = (m_nActiveCells[1] - 1); j >= 0; j--) {
      for(MInt i = (m_nActiveCells[2] - 1); i >= 0; i--) {
        cellId_org = i + (j + k * m_nActiveCells[1]) * m_nActiveCells[2];
        i_new = i + m_noGhostLayers;
        j_new = j + m_noGhostLayers;
        k_new = k + m_noGhostLayers;
        cellId = i_new + (j_new + k_new * m_nCells[1]) * m_nCells[2];

        for(MInt var = 0; var < getNoPPVars(); var++) {
          m_summedVars[var][cellId] = m_summedVars[var][cellId_org];
          m_summedVars[var][cellId_org] = F0;
        }

        if(m_averagingFavre) {
          for(MInt var = 0; var < getNoVars(); var++) {
            m_favre[var][cellId] = m_favre[var][cellId_org];
            m_favre[var][cellId_org] = F0;
          }
        }
      }
    }
  }
}

template <MFloat (FvStructuredSolver<3>::*pressure_func)(MInt) const>
void FvStructuredSolver3D::computePrimitiveVariables_() {
  const MFloat gammaMinusOne = m_gamma - 1.0;

  MFloat** const RESTRICT cvars = m_cells->variables;
  MFloat** const RESTRICT pvars = m_cells->pvariables;

  for(MInt k = m_noGhostLayers; k < m_nCells[0] - m_noGhostLayers; ++k) {
    for(MInt j = m_noGhostLayers; j < m_nCells[1] - m_noGhostLayers; ++j) {
      for(MInt i = m_noGhostLayers; i < m_nCells[2] - m_noGhostLayers; ++i) {
        const MInt cellId = cellIndex(i, j, k);
        const MFloat fRho = F1 / cvars[CV->RHO][cellId];
        MFloat velPOW2 = F0;
        for(MInt vel = 0; vel < nDim; ++vel) { // compute velocity
          pvars[vel][cellId] = cvars[vel][cellId] * fRho;
          velPOW2 += POW2(pvars[vel][cellId]);
        }

        // density and pressure:
        pvars[PV->RHO][cellId] = cvars[CV->RHO][cellId]; // density
        pvars[PV->P][cellId] =
            gammaMinusOne
            * (cvars[CV->RHO_E][cellId] - F1B2 * pvars[PV->RHO][cellId] * velPOW2 + (this->*pressure_func)(cellId));

        for(MInt ransVar = 0; ransVar < m_noRansEquations; ransVar++) {
          // TODO_SS labels:FV,totest Does it make sense to forbid negative values. BCs sometimes negate the neighbor
          // value
          cvars[CV->RANS_VAR[ransVar]][cellId] = mMax(cvars[CV->RANS_VAR[ransVar]][cellId], F0);
          pvars[PV->RANS_VAR[ransVar]][cellId] = cvars[CV->RANS_VAR[ransVar]][cellId] * fRho;
        }
      }
    }
  }
}

void FvStructuredSolver3D::computePrimitiveVariables() {
  if(noRansEquations(m_ransMethod) == 2)
    computePrimitiveVariables_<&FvStructuredSolver::pressure_twoEqRans>();
  else
    computePrimitiveVariables_();
}


void FvStructuredSolver3D::allocateSingularities() {
  for(MInt i = 0; i < m_hasSingularity; ++i) {
    MInt len[nDim];
    m_singularity[i].totalPoints = 1;
    m_singularity[i].totalCells = 1;

    for(MInt j = 0; j < nDim; j++) {
      len[j] = m_singularity[i].end[j] - m_singularity[i].start[j];
      m_singularity[i].totalPoints *= (len[j] + 1);
      m_singularity[i].totalCells *= len[j];
    }

    // 4 unknowns and 2*Nstar cells
    mAlloc(m_singularity[i].ReconstructionConstants, nDim + 1, m_singularity[i].totalCells * m_singularity[i].Nstar * 2,
           "ReconstructionConstants", 0.0, AT_);
  }
}

// see also the function in lbblckdxqy.cpp for this function
// implemented originally by Georg Eitel Amor for LB

void FvStructuredSolver3D::initFFTW(fftw_complex* uPhysField, fftw_complex* vPhysField, fftw_complex* wPhysField,
                                    MInt lx, MInt ly, MInt lz, MInt noPeakModes) {
  TRACE();
  // first check for odd numbers of cells in each direction
  if(lx % 2 != 0 || ly % 2 != 0 || lz % 2 != 0) {
    stringstream errorMessage;
    errorMessage << " FFTInit: Domain size must NOT be an odd number!: (lx)x(ly)x(lz)-> " << lx << "x" << ly << "x"
                 << lz << endl;
    mTerm(1, AT_, errorMessage.str());
  }

  m_log << " --- initializing FFTW --- " << endl;
  m_log << " domain size = " << lx << "x" << ly << "x" << lz << endl;

  MFloat waveVector[3], k0;

  complex<MFloat>* fourierCoefficient = new complex<MFloat>[3];

  fftw_complex *uHatField, *vHatField, *wHatField;

  fftw_plan planU, planV, planW;

  // 1) Allocation of Fourier coefficients
  uHatField = (fftw_complex*)fftw_malloc(lx * ly * lz * sizeof(fftw_complex));
  vHatField = (fftw_complex*)fftw_malloc(lx * ly * lz * sizeof(fftw_complex));
  wHatField = (fftw_complex*)fftw_malloc(lx * ly * lz * sizeof(fftw_complex));

  // 2) Creation of the plans for the FFTW
  planU = fftw_plan_dft_3d(lx, ly, lz, uHatField, uPhysField, FFTW_BACKWARD, FFTW_MEASURE);
  planV = fftw_plan_dft_3d(lx, ly, lz, vHatField, vPhysField, FFTW_BACKWARD, FFTW_MEASURE);
  planW = fftw_plan_dft_3d(lx, ly, lz, wHatField, wPhysField, FFTW_BACKWARD, FFTW_MEASURE);

  for(MInt p = 0; p < lx; p++) {
    for(MInt q = 0; q < ly; q++) {
      for(MInt r = 0; r < lz; r++) {
        MInt cellId = r + lz * (q + ly * p);
        // u-component
        uHatField[cellId][0] = 0.0;
        uHatField[cellId][1] = 0.0;
        uPhysField[cellId][0] = 0.0;
        uPhysField[cellId][1] = 0.0;

        // v-component
        vHatField[cellId][0] = 0.0;
        vHatField[cellId][1] = 0.0;
        vPhysField[cellId][0] = 0.0;
        vPhysField[cellId][1] = 0.0;

        // w-component
        wHatField[cellId][0] = 0.0;
        wHatField[cellId][1] = 0.0;
        wPhysField[cellId][0] = 0.0;
        wPhysField[cellId][1] = 0.0;
      }
    }
  }

  // IMPORTANT NOTICE ON USE OF FFTW:
  // comment from Georg Eitel Amor:
  // FFTW stores the coefficients for positive wavenumbers in the first half of the array,
  // and those for negative wavenumbers in reverse order in the second half.
  // [0, 1, ... , N/2-1, N/2, ... , N-1]
  //  - the entry for zero-wavenumber is at position 0
  //  - the k-th entry and the (N-k)th entry correspond to wavenumbers with opposite sign
  //  - the entry at position N/2 corresponds to the Nyquist wavenumber and appears only once

  // peak wave number of energy spectrum
  k0 = 2.0 * PI / (lx / noPeakModes);

  for(MInt p = 0; p <= lx / 2; p++) {
    for(MInt q = 0; q <= ly / 2; q++) {
      for(MInt r = 0; r <= lz / 2; r++) {
        // wave-vector: k(p,q,r) = (2 \pi p / lx, 2 \pi q / ly, 2 \pi r / lz)
        waveVector[0] = (p)*2.0 * PI / lx;
        waveVector[1] = (q)*2.0 * PI / ly;
        waveVector[2] = (r)*2.0 * PI / lz;

        getFourierCoefficients(waveVector, k0, fourierCoefficient);

        // 1. Positive frequencies:
        uHatField[r + lz * (q + ly * p)][0] = real(fourierCoefficient[0]);
        uHatField[r + lz * (q + ly * p)][1] = imag(fourierCoefficient[0]);

        vHatField[r + lz * (q + ly * p)][0] = real(fourierCoefficient[1]);
        vHatField[r + lz * (q + ly * p)][1] = imag(fourierCoefficient[1]);

        wHatField[r + lz * (q + ly * p)][0] = real(fourierCoefficient[2]);
        wHatField[r + lz * (q + ly * p)][1] = imag(fourierCoefficient[2]);

        // 2. Negative frequencies:
        if(p > 1 && q > 1 && r > 1) {
          if(p < lx / 2 && q < ly / 2 && r < lz / 2) {
            // since the physical velocity field is real, the coefficients for negative frequencies
            // are the complex conjugate of those for positive frequencies
            uHatField[(lz - r) + lz * ((ly - q) + ly * (lx - p))][0] = uHatField[r + lz * (q + ly * p)][0];
            uHatField[(lz - r) + lz * ((ly - q) + ly * (lx - p))][1] = -uHatField[r + lz * (q + ly * p)][1];

            vHatField[(lz - r) + lz * ((ly - q) + ly * (lx - p))][0] = vHatField[r + lz * (q + ly * p)][0];
            vHatField[(lz - r) + lz * ((ly - q) + ly * (lx - p))][1] = -vHatField[r + lz * (q + ly * p)][1];

            wHatField[(lz - r) + lz * ((ly - q) + ly * (lx - p))][0] = wHatField[r + lz * (q + ly * p)][0];
            wHatField[(lz - r) + lz * ((ly - q) + ly * (lx - p))][1] = -wHatField[r + lz * (q + ly * p)][1];
          }
        }
      }
    }
  }

  // Do Fourier transform (backward, see plan definition)
  // Definition in one dimension:
  // u(x) = \sum_{j=0}^{lx-1} \hat{u}_j exp(i 2 \pi j x / lx)

  fftw_execute(planU);
  fftw_execute(planV);
  fftw_execute(planW);


  // normalize (this preserves the norm of the basis functions)
  for(MInt p = 0; p < lx; p++) {
    for(MInt q = 0; q < ly; q++) {
      for(MInt r = 0; r < lz; r++) {
        uPhysField[r + lz * (q + ly * p)][0] /= sqrt(MFloat(lx * ly * lz));
        vPhysField[r + lz * (q + ly * p)][0] /= sqrt(MFloat(lx * ly * lz));
        wPhysField[r + lz * (q + ly * p)][0] /= sqrt(MFloat(lx * ly * lz));

        uPhysField[r + lz * (q + ly * p)][1] /= sqrt(MFloat(lx * ly * lz));
        vPhysField[r + lz * (q + ly * p)][1] /= sqrt(MFloat(lx * ly * lz));
        wPhysField[r + lz * (q + ly * p)][1] /= sqrt(MFloat(lx * ly * lz));
      }
    }
  }

  fftw_destroy_plan(planU);
  fftw_destroy_plan(planV);
  fftw_destroy_plan(planW);
  fftw_free(uHatField);
  fftw_free(vHatField);
  fftw_free(wHatField);
}

/**
 * \brief Generates a single complex coefficient of Fourier series
 *
 * Original Implementation by Georg Eitel Amor
 * for a given wavenumber k, and a certain energy spectrum
 * (see Appendix of Orszag, 1969)
 *
 */
void FvStructuredSolver3D::getFourierCoefficients(MFloat* k, MFloat k0, complex<MFloat>* fourierCoefficient) {
  TRACE();
  MFloat r[6], s[6], kAbs, energy;
  complex<MFloat> uHat, vHat, wHat;
  // complex<MFloat>* fourierCoefficient;

  kAbs = sqrt(k[0] * k[0] + k[1] * k[1] + k[2] * k[2]);

  // the zero-frequency component is always set to zero, so there is no offset
  if(approx(kAbs, F0, m_eps)) {
    // fourierCoefficient = new complex<MFloat>[3];

    fourierCoefficient[0] = complex<MFloat>(0, 0);
    fourierCoefficient[1] = complex<MFloat>(0, 0);
    fourierCoefficient[2] = complex<MFloat>(0, 0);
  } else {
    // energy = (kAbs/k0)*(kAbs/k0)*(kAbs/k0)*(kAbs/k0) * exp(-2.0*(kAbs/k0)*(kAbs/k0));
    // energy = pow(kAbs/k0,8.0) * exp(-4.0*(kAbs/k0)*(kAbs/k0)); // set spectral distribution
    energy = pow(kAbs / k0, 4.0) * exp(-2.0 * (kAbs / k0) * (kAbs / k0)); // set spectral distribution
    energy *= exp(2.0) * 0.499 * (m_Ma * LBCS)
              * (m_Ma * LBCS); // set maximal fluctuation amplitude to 20% of the freestream velocity (for 128^3: 0.88)

    // determine Fourier coefficients:
    // r and s are Independant random vector fields with independant
    // components (zero mean and rms according to energy spectrum).
    // Each vector has three components for k and another three for -k.

    for(MInt i = 0; i < 6; i++) {
      r[i] = randnormal(0.0, PI * sqrt(energy) / (SQRT2 * kAbs));
      // r[i] = randNumGen.randNorm(0.0, PI*sqrt(energy)/(SQRT2*kAbs));
      s[i] = randnormal(0.0, PI * sqrt(energy) / (SQRT2 * kAbs));
      // s[i] = randNumGen.randNorm(0.0, PI*sqrt(energy)/(SQRT2*kAbs));
    }

    uHat = (1.0 - k[0] * k[0] / (kAbs * kAbs)) * complex<MFloat>(r[0] + r[3], s[0] - s[3])
           - k[0] * k[1] / (kAbs * kAbs) * complex<MFloat>(r[1] + r[4], s[1] - s[4])
           - k[0] * k[2] / (kAbs * kAbs) * complex<MFloat>(r[2] + r[5], s[2] - s[5]);


    vHat = -k[1] * k[0] / (kAbs * kAbs) * complex<MFloat>(r[0] + r[3], s[0] - s[3])
           + (1.0 - k[1] * k[1] / (kAbs * kAbs)) * complex<MFloat>(r[1] + r[4], s[1] - s[4])
           - k[1] * k[2] / (kAbs * kAbs) * complex<MFloat>(r[2] + r[5], s[2] - s[5]);


    wHat = -k[2] * k[0] / (kAbs * kAbs) * complex<MFloat>(r[0] + r[3], s[0] - s[3])
           - k[2] * k[1] / (kAbs * kAbs) * complex<MFloat>(r[1] + r[4], s[1] - s[4])
           + (1.0 - k[2] * k[2] / (kAbs * kAbs)) * complex<MFloat>(r[2] + r[5], s[2] - s[5]);

    // fourierCoefficient = new complex<MFloat>[3];

    // fourierCoefficient[0] = complex<MFloat>(sqrt(2*energy)/SQRT2,sqrt(2*energy)/SQRT2);// uHat;
    // fourierCoefficient[1] = complex<MFloat>(sqrt(2*energy)/SQRT2,sqrt(2*energy)/SQRT2);//vHat;
    // fourierCoefficient[2] = complex<MFloat>(sqrt(2*energy)/SQRT2,sqrt(2*energy)/SQRT2);//wHat;

    fourierCoefficient[0] = uHat;
    fourierCoefficient[1] = vHat;
    fourierCoefficient[2] = wHat;
    // return fourierCoefficient;
  }
}

void FvStructuredSolver3D::computeReconstructionConstantsSVD() {
  MInt nghbr[30], dim = 0;
  MInt start[nDim], end[nDim];
  m_orderOfReconstruction = 1;
  const MInt recDim = (m_orderOfReconstruction == 2) ? (IPOW2(nDim) + 1) : nDim + 1;
  MInt maxNoSingularityRecNghbrIds = 14;
  MFloatScratchSpace tmpA(maxNoSingularityRecNghbrIds, recDim, AT_, "tmpA");
  MFloatScratchSpace tmpC(recDim, maxNoSingularityRecNghbrIds, AT_, "tmpC");
  MFloatScratchSpace weights(maxNoSingularityRecNghbrIds, AT_, "weights");
  MFloat counter = F0;
  MFloat avg = F0;
  MFloat maxc = F0;

  for(MInt i = 0; i < m_hasSingularity; ++i) {
    if(m_singularity[i].BC == -6000) {
      MInt totalCells = 1;
      MInt len1[nDim];

      //(p)reset the reconstruction constants
      for(MInt n = 0; n < nDim + 1; ++n) {
        for(MInt m = 0; m < m_singularity[i].totalCells * m_singularity[i].Nstar * 2; ++m) {
          m_singularity[i].ReconstructionConstants[n][m] = -999;
        }
      }

      for(MInt j = 0; j < nDim; j++) {
        len1[j] = m_singularity[i].end[j] - m_singularity[i].start[j];
        if(len1[j] != 0) totalCells *= len1[j];
      }

      for(MInt n = 0; n < nDim; ++n) {
        if(m_singularity[i].end[n] - m_singularity[i].start[n] > 1) {
          dim = n;
          start[n] = m_singularity[i].start[n] + 1;
          end[n] = m_singularity[i].end[n] - 1;
        } else {
          start[n] = m_singularity[i].start[n];
          end[n] = m_singularity[i].end[n];
        }
      }

      for(MInt kk = start[2]; kk < end[2]; ++kk) {
        for(MInt jj = start[1]; jj < end[1]; ++jj) {
          for(MInt ii = start[0]; ii < end[0]; ++ii) {
            MInt count = 0;
            MInt temp[nDim]{};
            temp[dim] = 1;

            nghbr[count++] = cellIndex(ii, jj, kk);
            nghbr[count++] = cellIndex(ii + temp[0], jj + temp[1], kk + temp[2]);

            // the coordinates of the corner where the viscousflux should be corrected.
            MInt ijk = getPointIdFromCell(ii + m_singularity[i].Viscous[0], jj + m_singularity[i].Viscous[1],
                                          kk + m_singularity[i].Viscous[2]);
            ijk = getPointIdfromPoint(ijk, 1, 1, 1);

            for(MInt m = 0; m < m_singularity[i].Nstar - 1; ++m) {
              MInt* change = m_singularity[i].displacement[m];
              nghbr[count++] = cellIndex(ii + change[0], jj + change[1], kk + change[2]);
              nghbr[count++] = cellIndex(ii + temp[0] + change[0], jj + temp[1] + change[1], kk + temp[2] + change[2]);
            }

            if(count != m_singularity[i].Nstar * 2) {
              cerr << "Something wrong with the singularities in the LS coeffiecient computation" << endl;
            }

            // weighted Least square
            weights.fill(F0);

            // Compute weights with RBF (take mean distance as R0)
            for(MInt n = 0; n < count; n++) {
              MInt nghbrId = nghbr[n];
              MFloat dxdx = F0;
              for(MInt m = 0; m < nDim; ++m) {
                dxdx += POW2(m_cells->coordinates[m][nghbrId] - m_grid->m_coordinates[m][ijk]);
              }

              weights[n] = 1 / dxdx; // RBF( dxdx, POW2( dist) );
            }

            MInt id2 = ii - start[0] + ((jj - start[1]) + (kk - start[2]) * len1[1]) * len1[0];
            MInt ID = id2 * m_singularity[i].Nstar * 2;

            MFloat condNum = computeRecConstSVD(ijk, count, nghbr, ID, i, tmpA, tmpC, weights, recDim);
            avg += condNum;
            maxc = mMax(maxc, condNum);
            counter += F1;
            if(condNum < F0 || condNum > 1e7 || std::isnan(condNum)) {
              cerr << domainId() << " SVD decomposition for pointId " << ijk
                   << " with large condition number: " << condNum << " num of neighbor" << count << "x" << recDim << " "
                   << " coords " << m_grid->m_coordinates[0][ijk] << ", " << m_grid->m_coordinates[1][ijk] << ", "
                   << m_grid->m_coordinates[2][ijk] << endl;
            }
          }
        }
      }
    }
  }
}


#include <numeric>
void FvStructuredSolver3D::exchange6002() {
  mTerm(1, "This 3D version is not tested yet!");
  //  if (m_rans && m_ransMethod!=RANS_KEPSILON)
  //    mTerm(1, "Porous RANS computation is only supported by k-epsilon model!");
  //  if (!m_porous)
  //    mTerm(1, "bc6002 requires the property porous to be set to true!");

  // 0) Check if initBc6002 is called for the first time
  //  for(MInt bcId_=0; bcId_ < bcId; ++bcId_) {
  //    if (m_physicalBCMap[bcId_]->BC==6002) return;
  //  }

  // Determine normal vectors and save for later use
  for(MInt bcId_ = 0; bcId_ < (MInt)m_structuredBndryCnd->m_physicalBCMap.size(); ++bcId_) {
    if(m_structuredBndryCnd->m_physicalBCMap[bcId_]->BC == 6002
       && m_structuredBndryCnd->m_physicalBCMap[bcId_]->Nstar == -1) {
      MInt* start = m_structuredBndryCnd->m_physicalBCMap[bcId_]->start1;
      MInt* end = m_structuredBndryCnd->m_physicalBCMap[bcId_]->end1;

      const MInt IJKP[nDim] = {1, m_nPoints[2], m_nPoints[1] * m_nPoints[2]};
      const MInt IJ[nDim] = {1, m_nCells[2], m_nCells[1] * m_nCells[2]};

      const MInt pp[3][12] = {{0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1},
                              {0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 1},
                              {0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0}};

      const MInt face = m_structuredBndryCnd->m_physicalBCMap[bcId_]->face;

      const MInt normalDir = face / 2;
      const MInt firstTangentialDir = (normalDir + 1) % nDim;
      const MInt secondTangentialDir = (normalDir + 2) % nDim;
      const MInt normalDirStart = start[normalDir];
      const MInt firstTangentialStart = start[firstTangentialDir];
      const MInt firstTangentialEnd = end[firstTangentialDir];
      const MInt secondTangentialStart = start[secondTangentialDir];
      const MInt secondTangentialEnd = end[secondTangentialDir];
      const MInt incp[nDim] = {IJKP[normalDir], IJKP[firstTangentialDir], IJKP[secondTangentialDir]};
      const MInt inc[nDim] = {IJ[normalDir], IJ[firstTangentialDir], IJ[secondTangentialDir]};

      const MInt n = (face % 2) * 2 - 1;                                       //-1,+1
      const MInt g1p = normalDirStart + 2 * ((MInt)(0.5 - (0.5 * (MFloat)n))); //+2,0
      const MInt g1 = normalDirStart + (MInt)(0.5 - (0.5 * (MFloat)n));        //+1,0

      for(MInt t1 = firstTangentialStart; t1 < firstTangentialEnd; t1++) {
        for(MInt t2 = secondTangentialStart; t2 < secondTangentialEnd; t2++) {
          const MInt cellIdG1 = g1 * inc[0] + t1 * inc[1] + t2 * inc[2];

          // compute four surrounding points of surface centroid
          const MInt ijk = g1p * incp[0] + t1 * incp[1] + t2 * incp[2];
          const MInt pp1 = getPointIdfromPoint(ijk, pp[normalDir][0], pp[normalDir][1], pp[normalDir][2]);
          const MInt pp2 = getPointIdfromPoint(ijk, pp[normalDir][3], pp[normalDir][4], pp[normalDir][5]);
          const MInt pp3 = getPointIdfromPoint(ijk, pp[normalDir][6], pp[normalDir][7], pp[normalDir][8]);
          MInt helper[nDim] = {0, 0, 0};
          helper[normalDir] = 1;
          const MInt pp3_ =
              getPointIdfromPoint(ijk, n * helper[0], n * helper[1], n * helper[2]); // point lying outside domain

          // compute the velocity of the surface centroid
          MFloat firstVec[nDim] = {F0, F0, F0};
          MFloat secondVec[nDim] = {F0, F0, F0};
          MFloat normalVec[nDim] = {F0, F0, F0};
          MFloat normalVec_[nDim]{};
          for(MInt dim = 0; dim < nDim; dim++) {
            firstVec[dim] = m_grid->m_coordinates[dim][pp2] - m_grid->m_coordinates[dim][pp1];
            secondVec[dim] = m_grid->m_coordinates[dim][pp3] - m_grid->m_coordinates[dim][pp1];
            normalVec_[dim] = m_grid->m_coordinates[dim][pp3_] - m_grid->m_coordinates[dim][pp1];
          }

          // compute normal vector of surface
          crossProduct(normalVec, firstVec, secondVec);
          const MFloat normalLength = sqrt(POW2(normalVec[0]) + POW2(normalVec[1]) + POW2(normalVec[2]));

          MFloat sgn = (std::inner_product(&normalVec[0], &normalVec[0] + nDim, &normalVec_[0], 0.0) < 0.0) ? -1 : 1;
          if(m_blockType == "fluid") sgn *= -1;

          for(MInt dim = 0; dim < nDim; dim++) {
            normalVec[dim] /= normalLength;
            m_cells->fq[FQ->NORMAL[dim]][cellIdG1] = sgn * normalVec[dim];
          }
        }
      }
    }
  }

  // Determine normal vector at singularities
  for(MInt i = 0; i < m_hasSingularity; ++i) {
    mTerm(1, "Not tested yet!");
    const auto& singularity = m_singularity[i];
    // only correct for bc 6000 not for bc 4000-5000
    if(singularity.BC == -6000) {
      MBool takeIt = false;
      for(MInt n = 0; n < singularity.Nstar; ++n) {
        if(singularity.BCsingular[n] == -6002) takeIt = true;
      }

      if(takeIt) {
        MInt start[nDim], end[nDim];
        for(MInt n = 0; n < nDim; ++n) {
          if(singularity.end[n] - singularity.start[n] > 1) {
            // dim=n;
            // start[n]=singularity.start[n]+1;
            start[n] = singularity.start[n] + 1;
            end[n] = singularity.end[n] - 1;
          } else {
            start[n] = singularity.start[n];
            end[n] = singularity.end[n];
          }
        }

        for(MInt kk = start[2]; kk < end[2]; ++kk) {
          for(MInt jj = start[1]; jj < end[1]; ++jj) {
            for(MInt ii = start[0]; ii < end[0]; ++ii) {
              const MInt IJ = cellIndex(ii, jj, kk);
              if(abs(m_cells->fq[FQ->NORMAL[0]][IJ]) > 1e-8) mTerm(1, "");
              for(MInt m = 0; m < 2; ++m) {
                const MInt* change = singularity.displacement[m];
                const MInt nghbr = cellIndex(ii + change[0], jj + change[1], kk + change[2]);
                for(MInt d = 0; d < nDim; ++d)
                  m_cells->fq[FQ->NORMAL[d]][IJ] += m_cells->fq[FQ->NORMAL[d]][nghbr];
              }
            }
          }
        }
      }
    }
  }

  // TODO_SS labels:FV,toenhance The exchange of all the normals is an overhead, because it is only needed at
  // singularity points
  /////////////////////////////////////////////////////////////////////////////
  /// GATHER & SEND
  /////////////////////////////////////////////////////////////////////////////
  //  MInt sendSizeTotal = 0;
  std::vector<MInt> receiveSizes;
  for(auto& snd : m_sndComm) {
    // TODO_SS labels:FV right now exchange at all 6000er not only 6002 because of singularities
    if(snd->bcId == -6000 || snd->bcId == -6002) {
      // Gather

      MInt size = 1;
      for(MInt dim = 0; dim < nDim; ++dim)
        size *= snd->endInfoCells[dim] - snd->startInfoCells[dim];
      //      std::vector<MFloat> snd->cellBuffer(size*(1+nDim));
      //      sendSizeTotal += size;

      MInt pos = 0;
      for(MInt k = snd->startInfoCells[2]; k < snd->endInfoCells[2]; k++) {
        for(MInt j = snd->startInfoCells[1]; j < snd->endInfoCells[1]; j++) {
          for(MInt i = snd->startInfoCells[0]; i < snd->endInfoCells[0]; i++) {
            const MInt cellId = i + (j + k * m_nCells[1]) * m_nCells[2];
            // TODO_SS labels:FV Latter only allocate FQ->POROSITY for m_blockType==porous
            //        if (m_blockType=="fluid")
            //          snd->cellBuffer[pos++] = 1;
            //        else
            snd->cellBuffer[pos] = m_cells->fq[FQ->POROSITY][cellId];
            for(MInt d = 0; d < nDim; ++d) {
              snd->cellBuffer[(1 + d) * size + pos] = m_cells->fq[FQ->NORMAL[d]][cellId];
            }
            ++pos;
          }
        }
      }

      // Send
      MInt tag = domainId() + (snd->tagHelper) * noDomains();
      MInt err = MPI_Isend((void*)&snd->cellBuffer[0], size * (nDim + 1), MPI_DOUBLE, snd->nghbrId, tag,
                           m_StructuredComm, &snd->mpi_request, AT_, "(void*)&snd->cellBuffer[0]");
      if(err) cout << "rank " << domainId() << " sending throws error " << endl;


      // Determine size of receive buffer
      // TODO_SS labels:FV right now exchange at all 6000er not only 6002 because of singularities
      size = 1;
      for(MInt dim = 0; dim < nDim; ++dim)
        size *= snd->endInfoCells[dim] - snd->startInfoCells[dim];
      receiveSizes.push_back(size);
    }
  }

  /////////////////////////////////////////////////////////////////////////////
  /// RECEIVE
  /////////////////////////////////////////////////////////////////////////////
  std::vector<MFloat> bufferRcv(std::accumulate(receiveSizes.begin(), receiveSizes.end(), 0) * (nDim + 1));
  MInt offset = 0;
  MInt cnt = 0;
  for(auto& rcv : m_rcvComm) {
    // TODO_SS labels:FV right now exchange at all 6000er not only 6002 because of singularities
    if(rcv->bcId == -6002 || rcv->bcId == -6000) {
      const MInt rcvSize = receiveSizes[cnt];
      MInt tag = rcv->nghbrId + (rcv->tagHelper) * noDomains();
      MInt err = MPI_Irecv(&bufferRcv[offset], rcvSize * (1 + nDim), MPI_DOUBLE, rcv->nghbrId, tag, m_StructuredComm,
                           &rcv->mpi_request, AT_, "(void*)&rcvSize");
      if(err) cout << "rank " << domainId() << " sending throws error " << endl;

      offset += rcvSize * (1 + nDim);
      ++cnt;
    }
  }

  /////////////////////////////////////////////////////////////////////////////
  /// WAIT
  /////////////////////////////////////////////////////////////////////////////
  for(auto& snd : m_sndComm) {
    if(snd->bcId == -6002 || snd->bcId == -6000) {
      MPI_Wait(&(snd->mpi_request), &(snd->mpi_status), AT_);
    }
  }


  for(auto& rcv : m_rcvComm) {
    if(rcv->bcId == -6002 || rcv->bcId == -6000) {
      MPI_Wait(&(rcv->mpi_request), &(rcv->mpi_status), AT_);
    }
  }


  /////////////////////////////////////////////////////////////////////////////
  /// SCATTER
  /////////////////////////////////////////////////////////////////////////////
  ScratchSpace<MFloat> normals_temp(m_noCells * nDim, AT_, "normal_temp");
  offset = 0;
  cnt = 0;
  for(auto& rcv : m_rcvComm) {
    if(rcv->bcId == -6002 || rcv->bcId == -6000) {
      // TODO_SS labels:FV right now exchange at all 6000er not only 6002 because of singularities

      MInt k2, j2, i2, id2;
      MInt step2[nDim];
      MInt start1[nDim];
      MInt start2[nDim];
      MInt end2[nDim];
      MInt len2[nDim];
      MInt totalCells = 1;
      MInt len1[nDim];

      for(MInt j = 0; j < nDim; j++) {
        len1[j] = rcv->endInfoCells[j] - rcv->startInfoCells[j];
        if(len1[j] != 0) totalCells *= len1[j];
        // added    check the step for RCV part !!!!!!!!important
        step2[rcv->orderInfo[j]] = rcv->stepInfo[j];
      }

      // Sanity check
      ASSERT(totalCells == receiveSizes[cnt], "");

      for(MInt j = 0; j < nDim; j++) {
        start2[j] = 0;
        end2[j] = len1[j] - 1;
        len2[rcv->orderInfo[j]] = len1[j];
        if(step2[j] < 0) {
          MInt dummy = start2[j];
          start2[j] = end2[j];
          end2[j] = dummy;
        }
      }

      MInt pos = 0;
      k2 = start2[2];
      for(MInt k = rcv->startInfoCells[2]; k < rcv->endInfoCells[2]; k++) {
        j2 = start2[1];
        for(MInt j = rcv->startInfoCells[1]; j < rcv->endInfoCells[1]; j++) {
          i2 = start2[0];
          for(MInt i = rcv->startInfoCells[0]; i < rcv->endInfoCells[0]; i++) {
            start1[rcv->orderInfo[0]] = i2;
            start1[rcv->orderInfo[1]] = j2;
            start1[rcv->orderInfo[2]] = k2;

            id2 = start1[0] + (start1[1] + start1[2] * len2[1]) * len2[0];
            const MInt cellId = i + (j + k * m_nCells[1]) * m_nCells[2];
            m_cells->fq[FQ->POROSITY][cellId] = bufferRcv[offset * (nDim + 1) + id2];
            for(MInt d = 0; d < nDim; ++d) {
              normals_temp[m_noCells * d + cellId] = bufferRcv[offset * (nDim + 1) + (d + 1) * totalCells + id2];
              //            m_cells->fq[FQ->NORMAL[d]][cellId] = bufferRcv[offset*(nDim+1)+(d+1)*totalCells+id2];
            }

            i2 += step2[0];
            pos++;
          }
          j2 += step2[1];
        }
        k2 += step2[2];
      }

      offset += totalCells;
      ++cnt;
    }
  }

  // Determine normal vector at singularities
  for(MInt i = 0; i < m_hasSingularity; ++i) {
    const auto& singularity = m_singularity[i];
    // only correct for bc 6000 not for bc 4000-5000
    if(singularity.BC == -6000) {
      MBool takeIt = false;
      for(MInt n = 0; n < singularity.Nstar; ++n) {
        if(singularity.BCsingular[n] == -6002) takeIt = true;
      }

      if(takeIt) {
        MInt start[nDim], end[nDim];
        for(MInt n = 0; n < nDim; ++n) {
          if(singularity.end[n] - singularity.start[n] > 1) {
            // dim=n;
            // start[n]=singularity.start[n]+1;
            start[n] = singularity.start[n] + 1;
            end[n] = singularity.end[n] - 1;
          } else {
            start[n] = singularity.start[n];
            end[n] = singularity.end[n];
          }
        }

        const MInt nstar = singularity.Nstar;

        for(MInt kk = start[2]; kk < end[2]; ++kk) {
          for(MInt jj = start[1]; jj < end[1]; ++jj) {
            for(MInt ii = start[0]; ii < end[0]; ++ii) {
              const MInt IJ = cellIndex(ii, jj, kk);
              MFloat temp[nDim];
              for(MInt d = 0; d < nDim; ++d) {
                temp[d] = m_cells->fq[FQ->NORMAL[d]][IJ];
              }
              for(MInt m = 0; m < nstar - 1; ++m) {
                const MInt* change = singularity.displacement[m];
                const MInt nghbr = cellIndex(ii + change[0], jj + change[1], kk + change[2]);
                for(MInt d = 0; d < nDim; ++d)
                  temp[d] += normals_temp[m_noCells * d + nghbr];
              }

              MFloat l = 0;
              for(MInt d = 0; d < nDim; ++d) {
                m_cells->fq[FQ->NORMAL[d]][IJ] = temp[d] / (2 * nstar);
                l += POW2(m_cells->fq[FQ->NORMAL[d]][IJ]);
              }
              l = sqrt(l);
              for(MInt d = 0; d < nDim; ++d) {
                m_cells->fq[FQ->NORMAL[d]][IJ] /= l;
              }

              for(MInt m = 0; m < nstar - 1; ++m) {
                const MInt* change = singularity.displacement[m];
                const MInt nghbr = cellIndex(ii + change[0], jj + change[1], kk + change[2]);
                for(MInt d = 0; d < nDim; ++d)
                  m_cells->fq[FQ->NORMAL[d]][nghbr] = m_cells->fq[FQ->NORMAL[d]][IJ];
              }
            }
          }
        }
      }
    }
  }

////////////////////////////////////////////////////////////////////////////////
#if 0
  // 1) Gather
  std::vector<MFloat> sendBuffer;
  std::vector<MInt> sendcounts;
  std::vector<MInt> snghbrs;
  std::vector<MInt> tags;
  for(auto& snd: m_sndComm) {    
    if (snd->bcId==-6002) {
      snghbrs.push_back(snd->nghbrId);
      tags.push_back(domainId()+(snd->tagHelper)*m_solver->noDomains());

      MInt size = 1;
      for (MInt dim = 0; dim < nDim; ++dim)
        size *= snd->endInfoCells[dim] - snd->startInfoCells[dim];
      sendcounts.push_back(size);

      sendBuffer.resize(pos+size);
      for(MInt j=snd->startInfoCells[1]; j<snd->endInfoCells[1]; j++) {
        for(MInt i=snd->startInfoCells[0]; i<snd->endInfoCells[0]; i++) {
          const MInt cellId = cellIndex(i,j);
          // TODO_SS labels:FV Latter only allocate FQ->POROSITY for m_blockType==porous
//          if (m_blockType=="fluid")
//            sendBuffer[pos++] = 1;
//          else
          sendBuffer[pos++] = m_cells->fq[FQ->POROSITY][cellId];
        }
      }
    }
  }

  const MInt noNeighborDomains = snghbrs.size();

  // 2) Send & receive
  std::vector<MInt> recvcounts(noNeighborDomains);
  std::vector<MInt> rdispls(noNeighborDomains);
  std::vector<MFloat> recvBuffer = maia::mpi::mpiExchangePointToPoint(&sendBuffer[0],
                                                                      &snghbrs[0],
                                                                      noNeighborDomains,
                                                                      &sendcounts[0],
                                                                      &snghbrs[0],
                                                                      noNeighborDomains,
                                                                      m_StructuredComm,
                                                                      m_solver->domainId(),
                                                                      1,
                                                                      recvcounts.data(),
                                                                      rdispls.data());
  // 2.5) Send & receive the tags
  std::vector<MInt> sendcounts2(noNeighborDomains, 1);
  std::vector<MInt> recvTags = maia::mpi::mpiExchangePointToPoint(&tags[0],
                                                                  &snghbrs[0],
                                                                  noNeighborDomains,
                                                                  &sendcounts2[0],
                                                                  &snghbrs[0],
                                                                  noNeighborDomains,
                                                                  m_StructuredComm,
                                                                  m_solver->domainId(),
                                                                  1);

  // 3) Scatter
  for(auto& rcv: m_rcvComm) {
    if (rcv->bcId==-6002) {
      const MInt tag = rcv->nghbrId+rcv->tagHelper*m_solver->noDomains();
      MInt n;
      for (n = 0; n < noNeighborDomains; ++n) {
        if (tag==recvTags[n])
          break;
      }
      if (n==noNeighborDomains) mTerm(1, "n == noNeighborDomains");

      const MFloat* const recvBuffer_ = &recvBuffer[rdispls[n]];
      const MInt noReceivedElements = recvcounts[n];

      /////////// following is analoge to  FvStructuredSolver2D::scatter()
      MInt j2, i2, id2;
      MInt  step2[nDim];
      MInt start1[nDim];
      MInt start2[nDim];
      MInt end2[nDim];
      MInt len2[nDim];
      MInt totalCells=1;
      MInt len1[nDim];

      for(MInt j=0; j<nDim; j++) {
        len1[j]=rcv->endInfoCells[j] - rcv->startInfoCells[j];
        if(len1[j]!=0) totalCells*=len1[j];
        //added    check the step for RCV part !!!!!!!!important
        step2[rcv->orderInfo[j]]=rcv->stepInfo[j];
      }

      //TODO_SS labels:FV,totest check if this assert makes sense
      ASSERT(noReceivedElements==totalCells, "noReceivedElements == totalCells");

      for(MInt j=0; j<nDim; j++) {
        start2[j]=0;
        end2[j]=len1[j]-1;
        len2[rcv->orderInfo[j]]=len1[j];
        if(step2[j]<0) {
          MInt dummy=start2[j];
          start2[j]=end2[j];
          end2[j]=dummy;
        }
      }


      MInt pos=0;
      j2=start2[1];
      for(MInt j=rcv->startInfoCells[1]; j<rcv->endInfoCells[1]; j++) {
        i2=start2[0];
        for(MInt i=rcv->startInfoCells[0]; i<rcv->endInfoCells[0]; i++) {
          start1[rcv->orderInfo[0]]=i2;
          start1[rcv->orderInfo[1]]=j2;

          id2=start1[0]+start1[1]*len2[0];
          const MInt cellId = i +(j*m_nCells[1]);
          m_cells->fq[FQ->POROSITY][cellId]= recvBuffer_[id2];

          i2+=step2[0];
          pos++;
        }
        j2+=step2[1];
      }
    }
  }
#endif
  ////////////////////////////////////////////////////////////////////////////////
}


inline MFloat FvStructuredSolver3D::getPSI(MInt I, MInt dim) {
  const MFloat FK = 18.0;
  const MInt IJK[3] = {1, m_nCells[2], m_nCells[1] * m_nCells[2]};
  const MInt IP1 = I + IJK[dim];
  const MInt IM1 = I - IJK[dim];
  const MInt IP2 = I + 2 * IJK[dim];

  const MFloat PIM2 = m_cells->pvariables[PV->P][IM1];
  const MFloat PIM1 = m_cells->pvariables[PV->P][I];
  const MFloat PIP2 = m_cells->pvariables[PV->P][IP2];
  const MFloat PIP1 = m_cells->pvariables[PV->P][IP1];

  const MFloat PSI =
      mMin(F1, FK
                   * mMax(mMax(fabs((PIM2 - PIM1) / mMin(PIM2, PIM1)), fabs((PIM1 - PIP1) / mMin(PIM1, PIP1))),
                          fabs((PIP1 - PIP2) / mMin(PIP1, PIP2))));
  return PSI;
}
