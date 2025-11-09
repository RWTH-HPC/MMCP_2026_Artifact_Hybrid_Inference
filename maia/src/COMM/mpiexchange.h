// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef MPIEXCHANGE_H_
#define MPIEXCHANGE_H_

#include <bitset>
#include <numeric>
#include <set>
#include "INCLUDE/maiamacro.h"
#include "INCLUDE/maiatypes.h"
#include "MEMORY/scratch.h"
#include "UTIL/debug.h"
#include "mpioverride.h"
#include "typetraits.h"

namespace maia {
namespace mpi {


//============================================================================================================
//============================================================================================================
//========================== THIS IS THE BASIC GENERIC EXCHANGE ROUTINE ======================================
//============= PLEASE REFER TO THIS ROUTINE IN MORE SPECIFIC EXCHANGE ROUTINES (see below)===================
//============================================================================================================
//============================================================================================================
/// \brief Generic exchange of data
/// \author Lennart Schneiders
/// \date October 2017
///
/// \tparam U Datatype to be exchanged
/// \param noExDomains Number of exchanging domains
/// \param exDomainId Rank of the exchanging domains
/// \param recvSize Number of values to be received from each exchange partner.
/// \param sendSize Number of values to be send to each exchange partner.
/// \param comm MPI Communicator to be used.
/// \param receiveBuffer Buffer used to receive data.
/// \param sendBuffer Buffer used to send data.
/// \param noDat Number of data elements to be send.
template <typename U>
inline void exchangeBuffer(const MInt noExDomains, const MInt* const exDomainId, const MInt* const recvSize,
                           const MInt* const sendSize, const MPI_Comm comm, U* const receiveBuffer,
                           const U* const sendBuffer, const MInt noDat = 1) {
  TRACE();

  // 0. prepare
  const MPI_Datatype DTYPE = type_traits<U>::mpiType();
  MInt tag = 612;
  ASSERT(tag < MPI_TAG_UB, "");
  ASSERT(noExDomains > 0, "");
  ScratchSpace<MPI_Request> sendRequests(noExDomains, AT_, "sendRequests");
  ScratchSpace<MPI_Request> recvRequests(noExDomains, AT_, "recvRequests");
  sendRequests.fill(MPI_REQUEST_NULL);
  recvRequests.fill(MPI_REQUEST_NULL);

  // Post all receive requests
  for(MInt i = 0, receiveCount = 0; i < noExDomains; i++) {
    MPI_Irecv(&(receiveBuffer[receiveCount]), noDat * recvSize[i], DTYPE, exDomainId[i], tag, comm, &recvRequests[i],
              AT_, "receiveBuffer");
    receiveCount += noDat * recvSize[i];
  }

  // Start sending
  for(MInt i = 0, sendCount = 0; i < noExDomains; i++) {
#if defined(HOST_Klogin)
    MPI_Issend(const_cast<U*>(&sendBuffer[sendCount]), noDat * sendSize[i], DTYPE, exDomainId[i], tag, comm,
               &sendRequests[i]);
#else
    MPI_Isend(&sendBuffer[sendCount], noDat * sendSize[i], DTYPE, exDomainId[i], tag, comm, &sendRequests[i], AT_,
              "&sendBuffer[sendCount]");
#endif
    sendCount += noDat * sendSize[i];
  }

  // Wait for all send and receive requests to finish
  MPI_Waitall(noExDomains, &recvRequests[0], MPI_STATUSES_IGNORE, AT_);
  MPI_Waitall(noExDomains, &sendRequests[0], MPI_STATUSES_IGNORE, AT_);
}
//============================================================================================================

/// \brief Generic exchange of data
/// \author Sven Berger
/// \date September 2019
///
/// \tparam U Datatype to be exchanged
/// \param sendSize Number of values to be send.
/// \param comm MPI Communicator to be used.
/// \param receiveBuffer Buffer used to receive data.
/// \param sendBuffer Buffer used to send data.
/// \param noDat Number of data elements to be send.
template <typename U>
inline MInt exchangeBufferAllToAll(const MInt sendSize, const MPI_Comm comm, U* const receiveBuffer,
                                   const U* const sendBuffer, const MInt noDat = 1) {
  TRACE();

  // 0. prepare
  const MPI_Datatype DTYPE = type_traits<U>::mpiType();
  const MInt tag = 612;
  const MInt noDomains = globalNoDomains();
  ASSERT(tag < MPI_TAG_UB, "");
  ScratchSpace<MPI_Request> mpiRequest(noDomains, AT_, "mpiRequest");
  MIntScratchSpace noToRecv(noDomains, AT_, "noToRecv");
  mpiRequest.fill(MPI_REQUEST_NULL);

  // 1.  gather number of values to be recved
  MPI_Allgather(&sendSize, 1, MPI_INT, &noToRecv[0], 1, MPI_INT, comm, AT_, "sendSize", "noToRecv");
  const MInt totalNo = std::accumulate(noToRecv.begin(), noToRecv.end(), -sendSize);
  ScratchSpace<U> temp(sendSize * noDat + 1, AT_, "temp");
  std::copy(&sendBuffer[0], &sendBuffer[sendSize * noDat], &temp[0]);

  // 2. send
  MInt recvCount = 0;
  for(MInt i = 0; i < noDomains; i++) {
    if(globalDomainId() == i) {
      MPI_Ibcast(&temp[0], noDat * sendSize, DTYPE, i, comm, &mpiRequest[i], AT_, "temp");
    } else {
      MPI_Ibcast(&receiveBuffer[recvCount], noDat * noToRecv[i], DTYPE, i, comm, &mpiRequest[i], AT_, "receiveBuffer");
      recvCount += noDat * noToRecv[i];
    }
  }

  MPI_Waitall(noDomains, &mpiRequest[0], MPI_STATUSES_IGNORE, AT_);
  return totalNo;
}
//============================================================================================================


/// \brief Generic exchange of data
/// \author Lennart Schneiders
/// \date October 2017
///
/// \tparam U
/// \param exDomains Domains with which to exchange Data.
/// \param noValuesToSend Number of values to send (in this type because of m_haloCells)
/// \param noValuesToRecv Number of values to recv (in this type because of m_windowCells)
/// \param comm MPI communicator to be used.
/// \param sendBuffer Buffer used to send data.
/// \param receiveBuffer Buffer used to receive data.
/// \param noDat Number of values per data entry.
template <typename U>
inline void exchangeBuffer(const std::vector<MInt>& exDomains,
                           const std::vector<std::vector<MInt>>& noValuesToRecv,
                           const std::vector<std::vector<MInt>>& noValuesToSend,
                           const MPI_Comm comm,
                           const U* const sendBuffer,
                           U* const receiveBuffer,
                           const MInt noDat = 1) {
  TRACE();

  // 0. prepare
  const auto noNghbrDomains = exDomains.size();
  ASSERT(noNghbrDomains > 0, "");
  ScratchSpace<MInt> noHaloCells(noNghbrDomains, AT_, "noHaloCellsg");
  ScratchSpace<MInt> noWindowCells(noNghbrDomains, AT_, "noWindowCellsg");
  for(MUint i = 0; i < noNghbrDomains; i++) {
    noHaloCells[i] = noValuesToRecv[i].size();
    noWindowCells[i] = noValuesToSend[i].size();
  }

  // 1. exchange
  exchangeBuffer(noNghbrDomains, exDomains.data(), &noHaloCells[0], &noWindowCells[0], comm, receiveBuffer, sendBuffer,
                 noDat);
}

/// \brief Generic exchange of data (std::vector version)
/// \author Lennart Schneiders
/// \date October 2017
///
/// \tparam U
/// \param exDomains Domains with which to exchange Data.
/// \param noValuesToSend Number of values to send.
/// \param noValuesToRecv Number of values to recv.
/// \param comm MPI communicator to be used.
/// \param sendBuffer Buffer used to send data.
/// \param receiveBuffer Buffer used to receive data.
/// \param noDat Number of values per data entry.
template <typename U>
inline void exchangeBuffer(const std::vector<MInt>& exDomains,
                           std::vector<MInt>& noValuesToRecv,
                           std::vector<MInt>& noValuesToSend,
                           const MPI_Comm comm,
                           const U* const sendBuffer,
                           U* const receiveBuffer,
                           const MInt noDat = 1) {
  TRACE();

  // 0. prepare
  const auto noNghbrDomains = exDomains.size();
  ASSERT(noNghbrDomains > 0, "");

  // 1. exchange
  exchangeBuffer(noNghbrDomains, exDomains.data(), &noValuesToRecv[0], &noValuesToSend[0], comm, receiveBuffer,
                 sendBuffer, noDat);
}

/// \brief Generic exchange of data (std::vector version)
/// \author Lennart Schneiders
/// \date October 2017
///
/// \tparam U
/// \param exDomains Domains with which to exchange Data.
/// \param noValuesToSend Number of values to send.
/// \param noValuesToRecv Number of values to recv.
/// \param comm MPI communicator to be used.
/// \param sendBuffer Buffer used to send data.
/// \param receiveBuffer Buffer used to receive data.
/// \param noDat Number of values per data entry.
template <typename U>
inline void exchangeBuffer(const std::vector<MInt>& exDomains,
                           const MInt* const noValuesToSend,
                           const MInt* const noValuesToRecv,
                           const MPI_Comm comm,
                           const U* const sendBuffer,
                           U* const receiveBuffer,
                           const MInt noDat = 1) {
  TRACE();

  // 0. prepare
  const auto noNghbrDomains = exDomains.size();
  ASSERT(noNghbrDomains > 0, "");

  // 1. exchange
  exchangeBuffer(noNghbrDomains, exDomains.data(), &noValuesToRecv[0], &noValuesToSend[0], comm, receiveBuffer,
                 sendBuffer, noDat);
}

/// \brief Generic exchange of data (std::vector version)
/// \author Lennart Schneiders
/// \date October 2017
///
/// \tparam U
/// \param exDomains Domains with which to exchange Data.
/// \param noValues Number of values to send and recv.
/// \param comm MPI communicator to be used.
/// \param sendBuffer Buffer used to send data.
/// \param receiveBuffer Buffer used to receive data.
/// \param noDat Number of values per data entry.
template <typename U>
inline void exchangeBuffer(const std::vector<MInt>& exDomains,
                           std::vector<MInt>& noValues,
                           const MPI_Comm comm,
                           const U* const sendBuffer,
                           U* const receiveBuffer,
                           const MInt noDat = 1) {
  TRACE();

  // 0. prepare
  const auto noNghbrDomains = exDomains.size();
  ASSERT(noNghbrDomains > 0, "");

  // 1. exchange
  exchangeBuffer(noNghbrDomains, exDomains.data(), &noValues[0], &noValues[0], comm, receiveBuffer, sendBuffer, noDat);
}

/// \brief Generic exchange of data (std::vector version)
/// \author Lennart Schneiders
/// \date October 2017
///
/// \tparam U
/// \param exDomains Domains with which to exchange Data.
/// \param noValues Number of values to send and recv.
/// \param comm MPI communicator to be used.
/// \param sendBuffer Buffer used to send data.
/// \param receiveBuffer Buffer used to receive data.
/// \param noDat Number of values per data entry.
template <typename U>
inline void exchangeValues(const std::vector<MInt>& exDomains, MInt noValues, const MPI_Comm comm,
                           const U* const sendBuffer, U* const receiveBuffer, const MInt noDat = 1) {
  TRACE();

  // 0. prepare
  const auto noNghbrDomains = exDomains.size();
  ASSERT(noNghbrDomains > 0, "");

  ScratchSpace<MInt> noV(noNghbrDomains, AT_, "noV");
  noV.fill(noValues);

  // 1. exchange
  exchangeBuffer(noNghbrDomains, exDomains.data(), &noV[0], &noV[0], comm, receiveBuffer, sendBuffer, noDat);
}


//-------------------------------------------------------------------------------------------


/**
 * \brief Generic exchange of data
 *  \author Lennart Schneiders
 *  \date October 2017
 */
template <typename U>
void exchangeData(const MInt noNghbrDomains, const MInt* const nghbrDomains, const MInt* const noHaloCells,
                  const MInt** const /*haloCells*/, const MInt* const noWindowCells, const MInt** const windowCells,
                  const MPI_Comm comm, const U* const data, U* const haloBuffer, const MInt noDat = 1) {
  TRACE();

  // 0. prepare
  MInt sendCount = std::accumulate(noWindowCells, noWindowCells + noNghbrDomains, 0);
  ScratchSpace<U> windowBuffer(mMax(1, noDat * sendCount), AT_, "windowBuffer");

  // 1. gather
  sendCount = 0;
  for(MInt i = 0; i < noNghbrDomains; i++) {
    for(MInt j = 0; j < noWindowCells[i]; j++) {
      for(MInt k = 0; k < noDat; k++) {
        windowBuffer[sendCount] = data[noDat * windowCells[i][j] + k];
        sendCount++;
      }
    }
  }

  // 2. exchange
  exchangeBuffer(noNghbrDomains, nghbrDomains, noHaloCells, noWindowCells, comm, haloBuffer, &windowBuffer[0], noDat);
}


//-------------------------------------------------------------------------------------------


/**
 * \brief Generic exchange of data
 *  \author Lennart Schneiders
 *  \date October 2017
 */
template <typename U>
void exchangeData(const std::vector<MInt>& nghbrDomains, const MInt* const noHaloCells, MInt** const haloCells,
                  MInt* const noWindowCells, MInt** const windowCells, const MPI_Comm comm, const U* const data,
                  U* const haloBuffer, const MInt noDat = 1) {
  TRACE();
  exchangeData(nghbrDomains.size(), nghbrDomains.data(), noHaloCells, haloCells, noWindowCells, windowCells, comm, data,
               haloBuffer, noDat);
}


//-------------------------------------------------------------------------------------------


/**
 * \brief Generic exchange of data
 *  \author Lennart Schneiders
 *  \date October 2017
 */
template <typename U>
void exchangeData(const MInt noNghbrDomains, const MInt* const nghbrDomains, const MInt* const noHaloCells,
                  const MInt** const haloCells, const MInt* const noWindowCells, const MInt** const windowCells,
                  const MPI_Comm comm, U* const data, const MInt noDat = 1) {
  TRACE();

  // 0. prepare
  MInt receiveCount = std::accumulate(noHaloCells, noHaloCells + noNghbrDomains, 0);
  ScratchSpace<U> haloBuffer(mMax(1, noDat * receiveCount), AT_, "haloBuffer");

  // 1. exchange
  exchangeData(noNghbrDomains, nghbrDomains, noHaloCells, haloCells, noWindowCells, windowCells, comm, data,
               &haloBuffer[0], noDat);

  // 2. scatter
  receiveCount = 0;
  for(MInt i = 0; i < noNghbrDomains; i++) {
    for(MInt j = 0; j < noHaloCells[i]; j++) {
      for(MInt k = 0; k < noDat; k++) {
        data[noDat * haloCells[i][j] + k] = haloBuffer[receiveCount];
        receiveCount++;
      }
    }
  }
}


//-------------------------------------------------------------------------------------------


/**
 * \brief Generic exchange of data
 *  \author Lennart Schneiders
 *  \date October 2017
 */
template <typename U>
void exchangeData(const std::vector<MInt>& nghbrDomains, const MInt* const noHaloCells, MInt** const haloCells,
                  MInt* const noWindowCells, MInt** const windowCells, const MPI_Comm comm, U* const data,
                  const MInt noDat = 1) {
  TRACE();
  exchangeData(nghbrDomains.size(), nghbrDomains.data(), noHaloCells, haloCells, noWindowCells, windowCells, comm, data,
               noDat);
}


//-------------------------------------------------------------------------------------------


/**
 * \brief Generic reverse exchange of data
 *  \author Lennart Schneiders
 *  \date October 2017
 */
template <typename U>
void reverseExchangeData(const MInt noNghbrDomains, const MInt* const nghbrDomains, const MInt* const noHaloCells,
                         const MInt** const haloCells, const MInt* const noWindowCells,
                         const MInt** const /*windowCells*/, const MPI_Comm comm, const U* const data,
                         U* const windowBuffer, const MInt noDat = 1) {
  TRACE();

  // 0. prepare
  MInt receiveCount = std::accumulate(noHaloCells, noHaloCells + noNghbrDomains, 0);
  ScratchSpace<U> haloBuffer(mMax(1, receiveCount), AT_, "haloBuffer");

  // 1. gather
  receiveCount = 0;
  for(MInt i = 0; i < noNghbrDomains; i++) {
    for(MInt j = 0; j < noHaloCells[i]; j++) {
      for(MInt k = 0; k < noDat; k++) {
        haloBuffer[receiveCount] = data[noDat * haloCells[i][j] + k];
        receiveCount++;
      }
    }
  }

  // 2. exchange
  exchangeBuffer(noNghbrDomains, nghbrDomains, noWindowCells, noHaloCells, comm, windowBuffer, &haloBuffer[0], noDat);
}


//-------------------------------------------------------------------------------------------


/**
 * \brief Generic exchange of data
 *  \author Lennart Schneiders
 *  \date October 2017
 */
template <typename U>
void exchangeData(const MInt noNghbrDomains, const MInt* const nghbrDomains,
                  const std::vector<std::vector<MInt>>& haloCellVec,
                  const std::vector<std::vector<MInt>>& windowCellVec, const MPI_Comm comm, U* const data,
                  const MInt noDat = 1) {
  TRACE();
  ASSERT(noNghbrDomains == (signed)windowCellVec.size() && noNghbrDomains == (signed)haloCellVec.size(), "");
  ScratchSpace<MInt> noHaloCells(mMax(1, noNghbrDomains), AT_, "noHaloCells");
  ScratchSpace<MInt> noWindowCells(mMax(1, noNghbrDomains), AT_, "noWindowCells");
  ScratchSpace<const MInt*> haloCells(mMax(1, noNghbrDomains), AT_, "haloCells");
  ScratchSpace<const MInt*> windowCells(mMax(1, noNghbrDomains), AT_, "windowCells");
  for(MInt i = 0; i < noNghbrDomains; i++) {
    noHaloCells[i] = (signed)haloCellVec[i].size();
    noWindowCells[i] = (signed)windowCellVec[i].size();
    haloCells[i] = haloCellVec[i].data();
    windowCells[i] = windowCellVec[i].data();
  }

  exchangeData(noNghbrDomains, nghbrDomains, noHaloCells.getPointer(), haloCells.getPointer(),
               noWindowCells.getPointer(), windowCells.getPointer(), comm, data, noDat);
}


//-------------------------------------------------------------------------------------------


/**
 * \brief Generic exchange of data
 *  \author Lennart Schneiders
 *  \date October 2017
 */
template <typename U>
void exchangeData(const MInt noNghbrDomains, const MInt* const nghbrDomains,
                  const std::vector<std::vector<MInt>>& haloCellVec,
                  const std::vector<std::vector<MInt>>& windowCellVec, const MPI_Comm comm, U* const data,
                  U* const haloBuffer, const MInt noDat = 1) {
  TRACE();
  ASSERT(noNghbrDomains == (signed)windowCellVec.size() && noNghbrDomains == (signed)haloCellVec.size(), "");
  ScratchSpace<MInt> noHaloCells(mMax(1, noNghbrDomains), AT_, "noHaloCells");
  ScratchSpace<MInt> noWindowCells(mMax(1, noNghbrDomains), AT_, "noWindowCells");
  ScratchSpace<const MInt*> haloCells(mMax(1, noNghbrDomains), AT_, "haloCells");
  ScratchSpace<const MInt*> windowCells(mMax(1, noNghbrDomains), AT_, "windowCells");
  for(MInt i = 0; i < noNghbrDomains; i++) {
    noHaloCells[i] = (signed)haloCellVec[i].size();
    noWindowCells[i] = (signed)windowCellVec[i].size();
    haloCells[i] = haloCellVec[i].data();
    windowCells[i] = windowCellVec[i].data();
  }

  exchangeData(noNghbrDomains, nghbrDomains, noHaloCells.getPointer(), haloCells.getPointer(),
               noWindowCells.getPointer(), windowCells.getPointer(), comm, data, haloBuffer, noDat);
}


//-------------------------------------------------------------------------------------------


/**
 * \brief Generic exchange of data
 *        NOTE: version called in exchangeData from the solvers!
 *  \author Lennart Schneiders
 *  \date October 2017
 */
template <typename U>
void exchangeData(const std::vector<MInt>& nghbrDomains, const std::vector<std::vector<MInt>>& haloCellVec,
                  const std::vector<std::vector<MInt>>& windowCellVec, const MPI_Comm comm, U* const data,
                  const MInt noDat = 1) {
  TRACE();
  const MInt noNghbrDomains = (signed)nghbrDomains.size();
  exchangeData(noNghbrDomains, nghbrDomains.data(), haloCellVec, windowCellVec, comm, data, noDat);
}


//-------------------------------------------------------------------------------------------


/**
 * \brief Generic exchange of data
 *  \author Lennart Schneiders
 *  \date October 2017
 */
template <std::size_t N>
void exchangeBitset(const std::vector<MInt>& nghbrDomains, const std::vector<std::vector<MInt>>& haloCellVec,
                    const std::vector<std::vector<MInt>>& windowCellVec, const MPI_Comm comm,
                    std::bitset<N>* const data, const MInt noCells, const MInt noDat = 1) {
  TRACE();
  static_assert(N <= 64, "conversion to ulong not appropriate, change to ullong!");
  const auto noNghbrDomains = (signed)nghbrDomains.size();
  ScratchSpace<MUlong> tmp_data(noCells, AT_, "tmp_data");
  for(MUint i = 0; i < nghbrDomains.size(); i++) {
    for(MInt windowCellId : windowCellVec[i]) {
      tmp_data[windowCellId] = data[windowCellId].to_ulong();
    }
  }
  exchangeData(noNghbrDomains, nghbrDomains.data(), haloCellVec, windowCellVec, comm, tmp_data.data(), noDat);
  for(MUint i = 0; i < nghbrDomains.size(); i++) {
    for(MInt haloCellId : haloCellVec[i]) {
      data[haloCellId] = std::bitset<N>(tmp_data[haloCellId]);
    }
  }
}


//-------------------------------------------------------------------------------------------


/**
 * \brief Generic exchange of data
 *  \author Lennart Schneiders
 *  \date October 2017
 */
template <typename U>
void exchangeData(const std::vector<MInt>& nghbrDomains, std::vector<std::vector<MInt>>& haloCellVec,
                  std::vector<std::vector<MInt>>& windowCellVec, const MPI_Comm comm, U* const data,
                  U* const haloBuffer, const MInt noDat = 1) {
  TRACE();
  const auto noNghbrDomains = (signed)nghbrDomains.size();
  exchangeData(noNghbrDomains, nghbrDomains.data(), haloCellVec, windowCellVec, comm, data, haloBuffer, noDat);
}


//-------------------------------------------------------------------------------------------


/**
 * \brief Generic reverse exchange of data
 *  \author Lennart Schneiders
 *  \date October 2017
 */
template <typename U>
void reverseExchangeData(const MInt noNghbrDomains, const MInt* const nghbrDomains,
                         const std::vector<std::vector<MInt>>& haloCellVec,
                         const std::vector<std::vector<MInt>>& windowCellVec, const MPI_Comm comm, const U* const data,
                         U* const windowBuffer, const MInt noDat = 1) {
  TRACE();
  ASSERT(noNghbrDomains == (signed)windowCellVec.size() && noNghbrDomains == (signed)haloCellVec.size(), "");
  ScratchSpace<MInt> noHaloCells(mMax(1, noNghbrDomains), AT_, "noHaloCells");
  ScratchSpace<MInt> noWindowCells(mMax(1, noNghbrDomains), AT_, "noWindowCells");
  ScratchSpace<const MInt*> haloCells(mMax(1, noNghbrDomains), AT_, "haloCells");
  ScratchSpace<const MInt*> windowCells(mMax(1, noNghbrDomains), AT_, "windowCells");
  for(MInt i = 0; i < noNghbrDomains; i++) {
    noHaloCells[i] = (signed)haloCellVec[i].size();
    noWindowCells[i] = (signed)windowCellVec[i].size();
    haloCells[i] = haloCellVec[i].data();
    windowCells[i] = windowCellVec[i].data();
  }

  reverseExchangeData(noNghbrDomains, nghbrDomains, noHaloCells.getPointer(), haloCells.getPointer(),
                      noWindowCells.getPointer(), windowCells.getPointer(), comm, data, windowBuffer, noDat);
}


//-------------------------------------------------------------------------------------------


/**
 * \brief Generic reverse exchange of data
 *  \author Lennart Schneiders
 *  \date October 2017
 */
template <typename U>
void reverseExchangeData(const std::vector<MInt>& nghbrDomains, const std::vector<std::vector<MInt>>& haloCellVec,
                         const std::vector<std::vector<MInt>>& windowCellVec, const MPI_Comm comm, const U* const data,
                         U* const windowBuffer, const MInt noDat = 1) {
  TRACE();
  const auto noNghbrDomains = (signed)nghbrDomains.size();
  reverseExchangeData(noNghbrDomains, nghbrDomains.data(), haloCellVec, windowCellVec, comm, data, windowBuffer, noDat);
}


//-------------------------------------------------------------------------------------------


/**
 * \brief Generic exchange of data
 *  \author Lennart Schneiders
 *  \date October 2017
 */
template <typename U>
void exchangeScattered(const std::vector<MInt>& nghbrDomains, std::vector<MInt>& sendDomainIndices,
                       std::vector<U>& sendData, const MPI_Comm comm, std::vector<MInt>& recvOffsets,
                       std::vector<U>& recvBuffer, const MInt noDat = 1) {
  TRACE();
  const auto noNghbrDomains = (signed)nghbrDomains.size();
  if(sendDomainIndices.size() * ((unsigned)noDat) != sendData.size()) {
    mTerm(1, AT_, "Invalid exchange buffer sizes.");
  }
  ScratchSpace<MInt> recvSize(noNghbrDomains, AT_, "recvSize");
  ScratchSpace<MInt> sendSize(noNghbrDomains, AT_, "sendSize");
  ScratchSpace<MInt> unity(noNghbrDomains, AT_, "unity");
  unity.fill(1);
  recvSize.fill(0);
  sendSize.fill(0);
  for(MInt sendDomainIndice : sendDomainIndices) {
    ASSERT(sendDomainIndice > -1 && sendDomainIndice < noNghbrDomains, "");
    sendSize(sendDomainIndice)++;
  }

  exchangeBuffer(noNghbrDomains, nghbrDomains.data(), &unity[0], &unity[0], comm, &recvSize[0], &sendSize[0]);

  MInt recvCnt = 0;
  MInt sendCnt = 0;
  for(MInt i = 0; i < noNghbrDomains; i++) {
    recvCnt += recvSize[i];
    sendCnt += sendSize[i];
  }

  std::vector<MInt> sendOffsets;
  recvOffsets.clear();
  sendOffsets.resize(noNghbrDomains + 1);
  recvOffsets.resize(noNghbrDomains + 1);

  recvOffsets[0] = 0;
  sendOffsets[0] = 0;
  for(MInt i = 0; i < noNghbrDomains; i++) {
    recvOffsets[i + 1] = recvOffsets[i] + recvSize[i];
    sendOffsets[i + 1] = sendOffsets[i] + sendSize[i];
  }

  std::vector<U> sendBuffer;
  recvBuffer.clear();
  sendBuffer.resize(sendCnt * noDat);
  recvBuffer.resize(recvCnt * noDat);

  std::fill(&sendSize[0], &sendSize[0] + noNghbrDomains, 0);
  for(MUint i = 0; i < sendDomainIndices.size(); i++) {
    MInt idx = sendDomainIndices[i];
    ASSERT(idx > -1 && idx < noNghbrDomains, "");
    for(MInt j = 0; j < noDat; j++) {
      sendBuffer[noDat * (sendOffsets[idx] + sendSize[idx]) + j] = sendData[noDat * i + j];
    }
    sendSize[idx]++;
  }

  exchangeBuffer(noNghbrDomains, nghbrDomains.data(), &recvSize[0], &sendSize[0], comm, recvBuffer.data(),
                 sendBuffer.data(), noDat);
}


//-------------------------------------------------------------------------------------------


/**
 * \brief Generic exchange of data
 *  \author Lennart Schneiders
 *  \date October 2017
 */
inline MUint getBufferSize(const std::vector<std::vector<MInt>>& exchangeCells) {
  return accumulate(exchangeCells.begin(), exchangeCells.end(), 0u,
                    [](const MInt& a, const std::vector<MInt>& b) { return a + static_cast<MInt>(b.size()); });
}


/// \brief Communicate variable amounts of data from each domain to all other domains. Data may also
///        reside on the same domain.
///
/// \author Ansgar Niemoeller (ansgar) <a.niemoeller@aia.rwth-aachen.de>
/// \date 2016-12-22
///
/// \param[in] sendBuffer Pointer to data send buffer.
/// \param[in] dataBlockSize Data solver size, e.g. number of variables per cell.
/// \param[in] dataSizeSend Amount of data to send to each domain, e.g. number of cells.
/// \param[in] dataSizeRecv Amount of data to receive from each domain.
/// \param[out] recvBuffer Pointer to data receive buffer.
template <class dataType>
void exchangeData(const dataType* const sendBuffer, const MInt domainId, const MInt noDomains, const MPI_Comm mpiComm,
                  const MInt dataBlockSize, const MInt* const dataSizeSend, const MInt* const dataSizeRecv,
                  dataType* const recvBuffer) {
  TRACE();

  using namespace maia;

  ScratchSpace<MPI_Request> recvRequests(noDomains, AT_, "recvRequests");
  ScratchSpace<MPI_Request> sendRequests(noDomains, AT_, "sendRequests");

  MIntScratchSpace sendOffsets(noDomains, AT_, "sendOffsets");
  MIntScratchSpace recvOffsets(noDomains, AT_, "recvOffsets");

  sendOffsets[0] = 0;
  recvOffsets[0] = 0;
  // Determine send/receive data offsets
  for(MInt i = 1; i < noDomains; i++) {
    sendOffsets[i] = sendOffsets[i - 1] + dataBlockSize * dataSizeSend[i - 1];
    recvOffsets[i] = recvOffsets[i - 1] + dataBlockSize * dataSizeRecv[i - 1];
  }

  // Start receiving ...
  MInt recvCount = 0;
  for(MInt i = 0; i < noDomains; i++) {
    // ... only if there is any data and this is not the current domain
    if(dataSizeRecv[i] > 0 && i != domainId) {
      MPI_Irecv(&recvBuffer[recvOffsets[i]], dataBlockSize * dataSizeRecv[i], type_traits<dataType>::mpiType(), i, i,
                mpiComm, &recvRequests[recvCount], AT_, "recvBuffer[recvOffsets[i]]");
      recvCount++;
    }
  }

  // Start sending ...
  MInt sendCount = 0;
  for(MInt i = 0; i < noDomains; i++) {
    // ... if there is any data
    if(dataSizeSend[i] > 0) {
      if(i == domainId) {
        // If this is the current domain just copy the data to the receive
        // buffer
        std::copy_n(&sendBuffer[sendOffsets[i]], dataBlockSize * dataSizeSend[i], &recvBuffer[recvOffsets[i]]);
      } else {
#if defined(HOST_Klogin)
        MPI_Isend(const_cast<dataType*>(&sendBuffer[sendOffsets[i]]), dataBlockSize * dataSizeSend[i],
                  type_traits<dataType>::mpiType(), i, domainId, mpiComm, &sendRequests[sendCount], AT_,
                  "const_cast<dataType*>(&sendBuffer[sendOffsets[i]])");
#else
        MPI_Isend(&sendBuffer[sendOffsets[i]], dataBlockSize * dataSizeSend[i], type_traits<dataType>::mpiType(), i,
                  domainId, mpiComm, &sendRequests[sendCount], AT_, "sendBuffer[sendOffsets[i]]");
#endif
        sendCount++;
      }
    }
  }

  // Finish receiving
  if(recvCount > 0) {
    MPI_Waitall(recvCount, &recvRequests[0], MPI_STATUSES_IGNORE, AT_);
  }
  // Finish sending
  if(sendCount > 0) {
    MPI_Waitall(sendCount, &sendRequests[0], MPI_STATUSES_IGNORE, AT_);
  }
}


/// \brief Communicate variable amounts of data from each domain to all neighboring domains.
template <class dataType>
void exchangeData(const dataType* const sendBuffer, const MInt domainId, const MInt noNghbrDomains,
                  const MInt* const nghbrDomainIds, const MPI_Comm mpiComm, const MInt dataBlockSize,
                  const MInt* const dataSizeSend, const MInt* const dataSizeRecv, dataType* const recvBuffer) {
  TRACE();

  using namespace maia;

  if(noNghbrDomains == 0) {
    return;
  }

  ScratchSpace<MPI_Request> recvRequests(noNghbrDomains, AT_, "recvRequests");
  ScratchSpace<MPI_Request> sendRequests(noNghbrDomains, AT_, "sendRequests");

  MIntScratchSpace sendOffsets(noNghbrDomains, AT_, "sendOffsets");
  MIntScratchSpace recvOffsets(noNghbrDomains, AT_, "recvOffsets");

  sendOffsets[0] = 0;
  recvOffsets[0] = 0;
  // Determine send/receive data offsets
  for(MInt i = 1; i < noNghbrDomains; i++) {
    sendOffsets[i] = sendOffsets[i - 1] + dataBlockSize * dataSizeSend[i - 1];
    recvOffsets[i] = recvOffsets[i - 1] + dataBlockSize * dataSizeRecv[i - 1];
  }

  // Start receiving ...
  MInt recvCount = 0;
  for(MInt i = 0; i < noNghbrDomains; i++) {
    const MInt nghbrDomainId = nghbrDomainIds[i];
    // ... only if there is any data and this is not the current domain
    if(dataSizeRecv[i] > 0 && nghbrDomainId != domainId) {
      MPI_Irecv(&recvBuffer[recvOffsets[i]], dataBlockSize * dataSizeRecv[i], type_traits<dataType>::mpiType(),
                nghbrDomainId, nghbrDomainId, mpiComm, &recvRequests[recvCount], AT_, "recvBuffer[recvOffsets[i]]");
      recvCount++;
    }
  }

  // Start sending ...
  MInt sendCount = 0;
  for(MInt i = 0; i < noNghbrDomains; i++) {
    const MInt nghbrDomainId = nghbrDomainIds[i];
    // ... if there is any data
    if(dataSizeSend[i] > 0) {
      if(nghbrDomainId == domainId) {
        // If this is the current domain just copy the data to the receive
        // buffer
        std::copy_n(&sendBuffer[sendOffsets[i]], dataBlockSize * dataSizeSend[i], &recvBuffer[recvOffsets[i]]);
      } else {
#if defined(HOST_Klogin)
        MPI_Isend(const_cast<dataType*>(&sendBuffer[sendOffsets[i]]), dataBlockSize * dataSizeSend[i],
                  type_traits<dataType>::mpiType(), nghbrDomainId, domainId, mpiComm, &sendRequests[sendCount], AT_,
                  "const_cast<dataType*>(&sendBuffer[sendOffsets[i]])");
#else
        MPI_Isend(&sendBuffer[sendOffsets[i]], dataBlockSize * dataSizeSend[i], type_traits<dataType>::mpiType(),
                  nghbrDomainId, domainId, mpiComm, &sendRequests[sendCount], AT_, "sendBuffer[sendOffsets[i]]");
#endif
        sendCount++;
      }
    }
  }

  // Finish receiving
  if(recvCount > 0) {
    MPI_Waitall(recvCount, &recvRequests[0], MPI_STATUSES_IGNORE, AT_);
  }
  // Finish sending
  if(sendCount > 0) {
    MPI_Waitall(sendCount, &sendRequests[0], MPI_STATUSES_IGNORE, AT_);
  }
}


/// Assemble data buffer according to given sorting order.
template <typename DataType>
void assembleDataBuffer(const MInt noCells, const MInt dataBlockSize, const DataType* const data,
                        const MInt* const sortedId, DataType* const buffer) {
  TRACE();

  for(MInt cellId = 0; cellId < noCells; cellId++) {
    // Skip cells without valid index
    if(sortedId[cellId] < 0) {
      continue;
    }

    std::copy_n(&data[dataBlockSize * cellId], dataBlockSize, &buffer[dataBlockSize * sortedId[cellId]]);
  }
}


/// Assemble given data in send buffer and communicate.
template <typename DataType>
void communicateData(const DataType* const data, const MInt noCells, const MInt* const sortedCellId,
                     const MInt noDomains, const MInt domainId, const MPI_Comm mpiComm,
                     const MInt* const noCellsToSendByDomain, const MInt* const noCellsToReceiveByDomain,
                     const MInt dataBlockSize, DataType* const buffer) {
  TRACE();

  ScratchSpace<DataType> sendBuffer(dataBlockSize * std::max(noCellsToSendByDomain[noDomains], 1), FUN_, "sendBuffer");

  maia::mpi::assembleDataBuffer(noCells, dataBlockSize, data, sortedCellId, &sendBuffer[0]);
  maia::mpi::exchangeData(&sendBuffer[0], domainId, noDomains, mpiComm, dataBlockSize, noCellsToSendByDomain,
                          noCellsToReceiveByDomain, buffer);
}


template <std::size_t N>
void communicateBitsetData(const std::bitset<N>* const data, const MInt noCells, const MInt* const sortedCellId,
                           const MInt noDomains, const MInt domainId, const MPI_Comm mpiComm,
                           const MInt* const noCellsToSendByDomain, const MInt* const noCellsToReceiveByDomain,
                           const MInt dataBlockSize, std::bitset<N>* const buffer) {
  TRACE();

  ScratchSpace<MUlong> data2(dataBlockSize * noCells, FUN_, "data2");
  ScratchSpace<MUlong> buffer2(dataBlockSize * noCellsToReceiveByDomain[noDomains], FUN_, "buffer2");

  for(MInt i = 0; i < (signed)data2.size(); i++) {
    data2[i] = data[i].to_ulong();
  }

  communicateData(data2.data(), noCells, sortedCellId, noDomains, domainId, mpiComm, noCellsToSendByDomain,
                  noCellsToReceiveByDomain, dataBlockSize, buffer2.data());

  for(MInt i = 0; i < (signed)buffer2.size(); i++) {
    buffer[i] = std::bitset<N>(buffer2[i]);
  }
}

/// \brief Communicate (optionally) solver-structured data sorted by a global Id
//
//  Sketch:
//   sendData, where locId0...noLocIds is not sorted
//                 locId0            locId1                  locId2               noLocIds
//            |nRows*nCols values|nRows*nCols values|nRows*nCols values|....|nRows*nCols values|
//
//   recvData: data is sorted globally with the provided mapping dataGlobalId(localId).
//             Each process contains data ranging for the globalIds dataOffsets(domainId) to dataOffsets(domainId+1)
//
//            globalId0=dataOffsets(domainId)  globalId1=dataOffsets(domainId)+1 globalIdn=dataOffsets(domainId+1)-1
//            |nRows*nCols values               |nRows*nCols values                  |....|  nRows*nCols values
//
/// \author konstantin
/// \date 2018-09-28
///
/// \param[in] sendData             Pointer to globally and locally unsorted (optionally solverstructured) data
/// \param[in] noLocIds             Number of local dataPoints
/// \param[in] noGlobalIds          Number of global dataPoints
/// \param[in] nRows                Number of rows in solverstructured send data
/// \param[in] nCols                Number of columns in solverstructured send data
/// \param[in] localToGlobal        Mapping of localIds to globalIds
/// \param[in] dataOffsets          The offsets of the global data points after communication
/// \param[in] recvData             Globally ordered (solver-structured) data array
template <typename DataType>
void communicateGlobalyOrderedData(DataType const* const sendData, const MInt noLocIds, const MInt noGlobalIds,
                                   const MInt nRows, const MInt nCols, const MInt noDomains, const MInt domainId,
                                   const MPI_Comm mpiComm, const MIntScratchSpace& localToGlobal,
                                   const MIntScratchSpace& dataOffsets, ScratchSpace<DataType>& recvData) {
  TRACE();

  ScratchSpace<std::vector<MInt>> sendIds(noDomains, AT_, "sendIds");
  ScratchSpace<std::vector<MInt>> recvIds(noDomains, AT_, "recvIds");
  MIntScratchSpace globalToLocal(noGlobalIds, AT_, "noGlobalIds");
  globalToLocal.fill(-1);
  for(MInt localId = 0; localId < noLocIds; localId++) {
    MInt globalId = localToGlobal(localId);
    globalToLocal(globalId) = localId;
    for(MLong dom = 0; dom < noDomains; dom++) {
      if(globalId < dataOffsets(dom + 1) && globalId >= dataOffsets(dom)) {
        sendIds[dom].push_back(globalId);
      }
    }
  }
  MIntScratchSpace sendCount(noDomains, AT_, "sendCount");
  MIntScratchSpace recvCount(noDomains, AT_, "recvCount");
  for(MLong dom = 0; dom < noDomains; dom++) {
    sendCount[dom] = sendIds[dom].size();
    sendIds[dom].shrink_to_fit();
  }
  MPI_Alltoall(&sendCount[0], 1, MPI_INT, &recvCount[0], 1, MPI_INT, mpiComm, AT_, "sendCount[0]", "recvCount[0]");
  for(MInt dom = 0; dom < noDomains; dom++) {
    recvIds[dom].resize(recvCount[dom]);
    recvIds[dom].shrink_to_fit();
  }
  MInt mpiCount = 0;
  ScratchSpace<MPI_Request> sendReq(noDomains, AT_, "sendReq");
  sendReq.fill(MPI_REQUEST_NULL);
  for(MInt dom = 0; dom < noDomains; dom++) {
    if(sendCount[dom] == 0) continue;
    MPI_Issend(&sendIds[dom].front(), sendCount[dom], MPI_INT, dom, 78, mpiComm, &sendReq[mpiCount++], AT_,
               "sendIds[dom].front()");
  }
  for(MInt dom = 0; dom < noDomains; dom++) {
    if(recvCount[dom] == 0) continue;
    MPI_Recv(&recvIds[dom].front(), recvCount[dom], MPI_INT, dom, 78, mpiComm, MPI_STATUS_IGNORE, AT_,
             "recvIds[dom].front()");
  }
  if(mpiCount > 0) MPI_Waitall(mpiCount, &sendReq[0], MPI_STATUSES_IGNORE, AT_);

  ScratchSpace<std::vector<DataType>> sendDataBuf(noDomains, AT_, "sendDataBuf");
  for(MInt dom = 0; dom < noDomains; dom++) {
    sendDataBuf[dom].resize(nRows * nCols * sendCount[dom]);
    sendDataBuf[dom].shrink_to_fit();
    for(MInt cnt = 0; cnt < sendCount[dom]; cnt++) {
      MInt globalId = sendIds[dom][cnt];
      MInt localId = globalToLocal(globalId);
      if(localId < 0) {
        std::cerr << "localId " << localId << " globalId " << globalId << " cnt " << cnt << " sendIds[dom].size() "
                  << sendIds[dom].size() << " sendCount[dom] " << sendCount[dom] << std::endl;
        mTerm(1, AT_, "localId not found -> communication failed.");
      }
      for(MInt row = 0; row < nRows; row++) {
        for(MInt col = 0; col < nCols; col++) {
          sendDataBuf[dom][cnt * nRows * nCols + nCols * row + col] =
              sendData[localId * nRows * nCols + nCols * row + col];
        }
      }
    }
  }

  MIntScratchSpace recvDispl(noDomains + 1, AT_, "recvDispl");
  recvDispl.fill(0);
  recvCount[0] = recvIds[0].size() * nRows * nCols;
  sendCount[0] = sendIds[0].size() * nRows * nCols;
  for(MLong dom = 1; dom < noDomains; dom++) {
    recvCount[dom] = recvIds[dom].size() * nRows * nCols;
    sendCount[dom] = sendIds[dom].size() * nRows * nCols;
    recvDispl[dom] = recvDispl[dom - 1] + recvCount[dom - 1];
  }
  recvDispl[noDomains] = recvDispl[noDomains - 1] + recvCount[noDomains - 1];

  sendReq.fill(MPI_REQUEST_NULL);
  mpiCount = 0;
  ScratchSpace<DataType> tmpData(recvDispl[noDomains], AT_, "tmpData");
  for(MInt dom = 0; dom < noDomains; dom++) {
    if(sendCount[dom] == 0) continue;
    MPI_Issend(&sendDataBuf[dom].front(), sendCount[dom], type_traits<DataType>::mpiType(), dom, 82, mpiComm,
               &sendReq[mpiCount++], AT_, "sendDataBuf[dom].front()");
  }
  for(MInt dom = 0; dom < noDomains; dom++) {
    if(recvCount[dom] == 0) continue;
    MPI_Recv(&tmpData[recvDispl[dom]], recvCount[dom], type_traits<DataType>::mpiType(), dom, 82, mpiComm,
             MPI_STATUS_IGNORE, AT_, "tmpData[recvDispl[dom]]");
  }
  if(mpiCount > 0) MPI_Waitall(mpiCount, &sendReq[0], MPI_STATUSES_IGNORE, AT_);

  recvDispl.fill(0);
  for(MLong dom = 1; dom < noDomains; dom++) {
    recvDispl[dom] = recvDispl[dom - 1] + recvIds[dom - 1].size();
  }

  for(MInt dom = 0; dom < noDomains; dom++) {
    for(MUint cnt = 0; cnt < recvIds[dom].size(); cnt++) {
      MInt globalId = recvIds[dom][cnt];
      MInt recvId = globalId - dataOffsets(domainId);
      if(recvId < 0 || recvId >= (dataOffsets(domainId + 1) - dataOffsets(domainId))) {
        std::cerr << "recvId " << recvId << " globalId " << globalId << "dataOffsets(domainId) "
                  << dataOffsets(domainId) << std::endl;
        mTerm(1, AT_, "recvId exceeds array dimensions.");
      }
      for(MInt row = 0; row < nRows; row++) {
        for(MInt col = 0; col < nCols; col++) {
          recvData[recvId * nCols * nRows + nCols * row + col] =
              tmpData[(recvDispl[dom] + cnt) * nCols * nRows + nCols * row + col];
        }
      }
    }
  }
}


/*
 * \brief Works with non-matching send and receive neighbors. snghbrs and rnghbrs
 *        needs to be consistent (see NDEBUG below). recvcounts and rdispls are
 *        determined, based on sendcounts. Don't need to care about receiveBuffer.
 *
 * \param[in] sendBuffer Pointer to data send buffer.
 * \param[in] snghbrs Neighbor domainIds to which to send data to
 * \param[in] nosnghbrs Number of domains to which to send data to
 * \param[in] sendcounts Number elements in send buffer for each neighbor in snghbrs
 * \param[in] rngbhrs Neighbor domainIds from which to receive data
 * \param[in] nornghbrs Number of domains from which to receive data
 * \param[in] mpi_comm Communicator
 * \param[in] domainId DomainId of current domain
 * \param[in] dataBlockSize
 * \param[out] Any container supporting move semantics is recommendable; you can assign
 *             the output of this function to a variable by auto
 * \param[out] recvcounts_ Entry i specifies number of elemenets to receive from rank rnghbrs[i]
 * \param[out] rdispls_ Entry i specifies the displacement relative to receiveBuffer at which to
 *                      place the data from rnghbrs[i]
 */
template <class T>
std::vector<T> mpiExchangePointToPoint(const T* const sendBuffer,
                                       const MInt* const snghbrs,
                                       const MInt nosnghbrs,
                                       const MInt* const sendcounts,
                                       const MInt* const rnghbrs,
                                       const MInt nornghbrs,
                                       const MPI_Comm& mpi_comm,
                                       const MInt domainId,
                                       const MInt dataBlockSize,
                                       MInt* const recvcounts_ = nullptr,
                                       MInt* const rdispls_ = nullptr) {
  TRACE();

#ifndef NDEBUG
  {
    // TODO labels:COMM,totest uncomment the next two lines, is it necessary to be sorted?
    //  if (!std::is_sorted(&snghbrs[0], &snghbrs[nosnghbrs])) TERMM(-1, "");
    //  if (!std::is_sorted(&rnghbrs[0], &rnghbrs[nornghbrs])) TERMM(-1, "");
    // Check for duplicates
    std::set<MInt> s1(&snghbrs[0], &snghbrs[0] + nosnghbrs);
    ASSERT(s1.size() == static_cast<unsigned>(nosnghbrs), "");
    std::set<MInt> s2(&rnghbrs[0], &rnghbrs[0] + nornghbrs);
    ASSERT(s2.size() == static_cast<unsigned>(nornghbrs), "");

    /*
     * Check if send and receive neighbor domainIds are consistent
     */

    // Determine all neighbor domains of current domain
    MInt noDomains;
    MPI_Comm_size(mpi_comm, &noDomains);
    ScratchSpace<MChar> isReceiveDomain(noDomains, AT_, "isReceiveDomain");
    std::fill(isReceiveDomain.begin(), isReceiveDomain.end(), 0);
    for(MInt d = 0; d < nosnghbrs; d++) {
      isReceiveDomain[snghbrs[d]] = 1;
    }

    // Exchange with all domains
    MPI_Alltoall(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, &isReceiveDomain[0], 1, type_traits<MChar>::mpiType(), mpi_comm,
                 AT_, "MPI_IN_PLACE", "isReceiveDomain[0]");

    // Check if the domains, which will send to the current domain are also included in rnghbrs of current domain
    MInt nornghbrs_ = std::accumulate(isReceiveDomain.begin(), isReceiveDomain.end(), 0);
    if(nornghbrs != nornghbrs_) {
      TERMM(-1, "True number of domains sending to current domain " + std::to_string(domainId)
                    + " differs from expected number " + std::to_string(nornghbrs));
    }
    for(MInt d = 0; d < nornghbrs; d++) {
      if(!isReceiveDomain[rnghbrs[d]]) {
        TERMM(1, "Domain " + std::to_string(domainId) + " has domain " + std::to_string(rnghbrs[d])
                     + " as a neighbor but domain " + std::to_string(rnghbrs[d])
                     + " has nothing to send to current domain.");
      }
    }

    MPI_Barrier(mpi_comm, AT_);
  }
#endif

  // Check if I am a useless domain, which has nothing to send or receive
  if(nosnghbrs == 0 && nornghbrs == 0) {
    std::vector<T> _;
    return _;
  }

  // 1) Determine recvcounts
  ScratchSpace<MPI_Request> recvRequests(std::max(1, nornghbrs), AT_, "recvRequests");
  ScratchSpace<MInt> recvcounts(std::max(1, nornghbrs), AT_, "recvcounts");
  for(MInt d = 0; d < nornghbrs; ++d) {
    recvcounts[d] = 0; // is it necessary?
    MPI_Irecv(&recvcounts[d], 1, type_traits<MInt>::mpiType(), rnghbrs[d], rnghbrs[d], mpi_comm, &recvRequests[d], AT_,
              "recvcounts[d]");
  }

  ScratchSpace<MPI_Request> sendRequests(std::max(1, nosnghbrs), AT_, "sendRequests");
  for(MInt d = 0; d < nosnghbrs; ++d) {
    MPI_Isend(&sendcounts[d], 1, type_traits<MInt>::mpiType(), snghbrs[d], domainId, mpi_comm, &sendRequests[d], AT_,
              "sendcounts[d]");
  }

  MPI_Waitall(nornghbrs, &recvRequests[0], MPI_STATUSES_IGNORE, AT_);
  MPI_Waitall(nosnghbrs, &sendRequests[0], MPI_STATUSES_IGNORE, AT_);
  if(recvcounts_) std::copy_n(&recvcounts[0], nornghbrs, recvcounts_);

  /*
  std::cout << "MYRANK=" << domainId << ": " << "sendinfo : ";
  for (MInt i = 0; i < nosnghbrs; ++i) {
    std::cout << sendcounts[i] << "(" << snghbrs[i] << "),";
  }
  std::cout << std::endl;
  std::cout << "MYRANK=" << domainId << ": " << "recvinfo : ";
  for (MInt i = 0; i < nornghbrs; ++i) {
    std::cout << recvcounts[i] << "(" << rnghbrs[i] << "),";
  }
  std::cout << "   " << std::accumulate(&recvcounts[0], &recvcounts[nornghbrs-1], 0) * dataBlockSize << std::endl <<
  std::flush; MPI_Barrier(mpi_comm, AT_);
  */


  // 2) Alocate memory for receiveBuffer
  const MInt totalBufferSize = std::accumulate(&recvcounts[0], &recvcounts[0] + nornghbrs, 0) * dataBlockSize;
  std::vector<T> receiveBuffer(totalBufferSize);

  // 3) take into account dataBlockSize
  std::vector<MInt> sdispls(nosnghbrs);
  std::vector<MInt> rdispls(nornghbrs);

  //  sdispls[0] = 0;
  //  rdispls[0] = 0;
  for(MInt i = 1; i < nosnghbrs; i++) {
    sdispls[i] = sdispls[i - 1] + dataBlockSize * sendcounts[i - 1];
  }
  for(MInt i = 1; i < nornghbrs; i++) {
    rdispls[i] = rdispls[i - 1] + dataBlockSize * recvcounts[i - 1];
  }
  if(rdispls_) std::copy_n(&rdispls[0], nornghbrs, rdispls_);

  // 4) Send & receive actual data
  std::fill(recvRequests.begin(), recvRequests.end(), MPI_REQUEST_NULL);
  for(MInt d = 0; d < nornghbrs; ++d) {
    if(recvcounts[d] > 0) {
      MPI_Irecv(&receiveBuffer[rdispls[d]], recvcounts[d] * dataBlockSize, type_traits<T>::mpiType(), rnghbrs[d],
                rnghbrs[d], mpi_comm, &recvRequests[d], AT_, "receiveBuffer[rdispls[d]]");
    }
  }

  std::fill(sendRequests.begin(), sendRequests.end(), MPI_REQUEST_NULL);
  for(MInt d = 0; d < nosnghbrs; ++d) {
    if(sendcounts[d] > 0) {
      MPI_Isend(&sendBuffer[sdispls[d]], sendcounts[d] * dataBlockSize, type_traits<T>::mpiType(), snghbrs[d], domainId,
                mpi_comm, &sendRequests[d], AT_, "sendBuffer[sdispls[d]]");
    }
  }

  MPI_Waitall(nornghbrs, &recvRequests[0], MPI_STATUSES_IGNORE, AT_);
  MPI_Waitall(nosnghbrs, &sendRequests[0], MPI_STATUSES_IGNORE, AT_);

  return receiveBuffer;
}

/**
 * \brief  Generic exchange from halo to window cells,
 *         however in this case the value in the halo-cell is added
 *         to the values in the windowcell!
 *  \author Tim Wegmann
 *  \date March 2021
 */
template <typename U>
void reverseExchangeAddData(const std::vector<MInt>& nghbrDomains, const std::vector<std::vector<MInt>>& haloCellVec,
                            const std::vector<std::vector<MInt>>& windowCellVec, const MPI_Comm comm, U** const data,
                            const MInt noDat = 1) {
  TRACE();

  // 0. prepare
  const MInt noNghbrDomains = (signed)nghbrDomains.size();

  ASSERT(noNghbrDomains == (signed)windowCellVec.size() && noNghbrDomains == (signed)haloCellVec.size(), "");

  ScratchSpace<MInt> noHaloCells(mMax(1, noNghbrDomains), AT_, "noHaloCells");
  ScratchSpace<MInt> noWindowCells(mMax(1, noNghbrDomains), AT_, "noWindowCells");
  ScratchSpace<const MInt*> haloCells(mMax(1, noNghbrDomains), AT_, "haloCells");
  ScratchSpace<const MInt*> windowCells(mMax(1, noNghbrDomains), AT_, "windowCells");

  for(MInt i = 0; i < noNghbrDomains; i++) {
    noHaloCells[i] = (signed)haloCellVec[i].size();
    noWindowCells[i] = (signed)windowCellVec[i].size();
    haloCells[i] = haloCellVec[i].data();
    windowCells[i] = windowCellVec[i].data();
  }

  MInt receiveCount = std::accumulate(&noHaloCells[0], &noHaloCells[0] + noNghbrDomains, 0);
  ScratchSpace<U> haloBuffer(mMax(1, noDat * receiveCount), AT_, "haloBuffer");

  MInt receiveSize = std::accumulate(&noWindowCells[0], &noWindowCells[0] + noNghbrDomains, 0);
  ScratchSpace<U> windowBuffer(mMax(1, noDat * receiveSize), AT_, "haloBuffer");

  // 1. gather
  receiveCount = 0;
  for(MInt i = 0; i < noNghbrDomains; i++) {
    for(MInt j = 0; j < noHaloCells[i]; j++) {
      for(MInt k = 0; k < noDat; k++) {
        haloBuffer[receiveCount] = data[haloCells[i][j]][k];
        receiveCount++;
      }
    }
  }

  // 2. exchange
  exchangeBuffer(noNghbrDomains, nghbrDomains.data(), noWindowCells.getPointer(), noHaloCells.getPointer(), comm,
                 &windowBuffer[0], &haloBuffer[0], noDat);

  // 3. scatter
  receiveSize = 0;
  for(MInt i = 0; i < noNghbrDomains; i++) {
    for(MInt j = 0; j < noWindowCells[i]; j++) {
      for(MInt k = 0; k < noDat; k++) {
        data[windowCells[i][j]][k] += windowBuffer[receiveSize];
        receiveSize++;
      }
    }
  }
}

/**
 * \brief  Generic exchange from halo to window cells,
 *         however in this case the value in the halo-cell is added
 *         to the values in the windowcell, simple pointer version
 *  \author Tim Wegmann
 *  \date March 2021
 */
template <typename U>
void reverseExchangeAddData(const std::vector<MInt>& nghbrDomains, const std::vector<std::vector<MInt>>& haloCellVec,
                            const std::vector<std::vector<MInt>>& windowCellVec, const MPI_Comm comm, U* const data,
                            const MInt noDat = 1) {
  TRACE();

  // 0. prepare
  const MInt noNghbrDomains = (signed)nghbrDomains.size();

  ASSERT(noNghbrDomains == (signed)windowCellVec.size() && noNghbrDomains == (signed)haloCellVec.size(), "");

  ScratchSpace<MInt> noHaloCells(mMax(1, noNghbrDomains), AT_, "noHaloCells");
  ScratchSpace<MInt> noWindowCells(mMax(1, noNghbrDomains), AT_, "noWindowCells");
  ScratchSpace<const MInt*> haloCells(mMax(1, noNghbrDomains), AT_, "haloCells");
  ScratchSpace<const MInt*> windowCells(mMax(1, noNghbrDomains), AT_, "windowCells");

  for(MInt i = 0; i < noNghbrDomains; i++) {
    noHaloCells[i] = (signed)haloCellVec[i].size();
    noWindowCells[i] = (signed)windowCellVec[i].size();
    haloCells[i] = haloCellVec[i].data();
    windowCells[i] = windowCellVec[i].data();
  }

  MInt receiveCount = std::accumulate(&noHaloCells[0], &noHaloCells[0] + noNghbrDomains, 0);
  ScratchSpace<U> haloBuffer(mMax(1, noDat * receiveCount), AT_, "haloBuffer");

  MInt receiveSize = std::accumulate(&noWindowCells[0], &noWindowCells[0] + noNghbrDomains, 0);
  ScratchSpace<U> windowBuffer(mMax(1, noDat * receiveSize), AT_, "haloBuffer");

  // 1. gather
  receiveCount = 0;
  for(MInt i = 0; i < noNghbrDomains; i++) {
    for(MInt j = 0; j < noHaloCells[i]; j++) {
      for(MInt k = 0; k < noDat; k++) {
        haloBuffer[receiveCount] = data[noDat * haloCells[i][j] + k];
        receiveCount++;
      }
    }
  }

  // 2. exchange
  exchangeBuffer(noNghbrDomains, nghbrDomains.data(), noWindowCells.getPointer(), noHaloCells.getPointer(), comm,
                 &windowBuffer[0], &haloBuffer[0], noDat);

  // 3. scatter
  receiveSize = 0;
  for(MInt i = 0; i < noNghbrDomains; i++) {
    for(MInt j = 0; j < noWindowCells[i]; j++) {
      for(MInt k = 0; k < noDat; k++) {
        data[noDat * windowCells[i][j] + k] += windowBuffer[receiveSize];
        receiveSize++;
      }
    }
  }
}


} // namespace mpi
} // namespace maia

#endif // MPIEXCHANGE_H_
