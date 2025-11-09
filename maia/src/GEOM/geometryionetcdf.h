// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef GEOMETRYIONETCDF_H
#define GEOMETRYIONETCDF_H

#include "COMM/mpioverride.h"
#include "geometryiobase.h"
#if defined(MAIA_WINDOWS)
#define NC_MAX_NAME 256
#else
#include "IO/maiapnetcdf.h"
#endif

#ifndef MAIA_MS_COMPILER
class PARALLELIO_DEFAULT_BACKEND;
#else
class ParallelIoHdf5;
#endif

class GeometryIONetcdf : public GeometryIOBase {
#ifdef MAIA_MS_COMPILER
  using ParallelIo = ParallelIoHdf5;
#else
  using ParallelIo = PARALLELIO_DEFAULT_BACKEND;
#endif
 public:
  GeometryIONetcdf(const MPI_Comm comm) : GeometryIOBase(comm), m_mpiComm(comm) {}
  geometryAssembly* readPropertyFile(MString fileName);
  void readPropertyFileOldIOMethod(ParallelIo* parallelIo);
  void readPropertyFileNewIOMethod(ParallelIo* parallelIo);

  MInt segmentCount();
  MPI_Comm mpiComm() const { return m_mpiComm; }

 private:
  geometryAssembly* m_geometryAssembly;
  MInt m_noSegments;
  MInt m_noBodies;
  std::pair<geometryPropertyIterator, geometryPropertyIterator> m_pair;
  geometryPropertyMap* m_geometryPropertyMap;
  bodyMap* m_bodyMap;
  MBool m_newIOMethod;

  const MPI_Comm m_mpiComm = MPI_COMM_NULL;
  void readBodiesOldIOMethod(ParallelIo*);
  void readBodiesNewIOMethod(ParallelIo*);

  void buildDefaultBody();
  void makeProperty(GeometryProperty*, MString, ParallelIo*);
  void writeProperties(const MChar* fileName, geometryPropertyMap* pMap);
  MBool checkBodyConsistency();
  MBool checkGeometryPropertyConsistency();

  void distributeBodyProperties();
  void distributeGeometryProperties();
  void receiveBodyProperties();
  void receiveGeometryProperties();
};

#endif
