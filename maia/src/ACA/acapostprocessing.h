// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef ACAPOSTPROCESSING_H_
#define ACAPOSTPROCESSING_H_

/** \brief Abstract class for ACA post processing
 */
class AcaPostProcessing {
 protected:
  MPI_Comm m_mpiComm;
  MInt m_rank;
  MInt m_noObservers = -1;
  MInt m_noGlobalObservers = -1;
  MInt m_offsetObserver = 0;
  MInt m_noSamples = -1;
  MString m_outPath;
  const MFloat* m_coords = nullptr;
  const MFloat* m_frequencies = nullptr;

 public:
  AcaPostProcessing(const MPI_Comm comm) : m_mpiComm(comm) { MPI_Comm_rank(m_mpiComm, &m_rank); }
  virtual ~AcaPostProcessing() = default;
  virtual void init_() = 0;
  void init(const MInt noObservers, const MInt noGlobalObservers, const MInt offsetObserver, const MInt noSamples,
            const MFloat* const coords, const MFloat* const frequencies, const MString outPath) {
    m_noObservers = noObservers;
    m_noGlobalObservers = noGlobalObservers;
    m_offsetObserver = offsetObserver;
    m_noSamples = noSamples;
    m_coords = coords;
    m_frequencies = frequencies;
    m_outPath = outPath;
    init_();
  }
  virtual void calc(const MInt observerId, const MFloat* const data, const MFloat* const dataComplex = nullptr) = 0;
  virtual void finish() = 0;
};
#endif // ACAPOSTPROCESSING_H_
