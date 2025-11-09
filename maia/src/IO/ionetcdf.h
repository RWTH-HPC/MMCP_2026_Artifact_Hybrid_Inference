// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef IONETCDF_H
#define IONETCDF_H

#include "contexttypes.h"
#include "maiapnetcdf.h"

#ifndef MAIA_MS_COMPILER
class PARALLELIO_DEFAULT_BACKEND;
#else
class ParallelIoHdf5;
#endif

class IONetcdf {
#ifdef MAIA_MS_COMPILER
  using ParallelIo = ParallelIoHdf5;
#else
  using ParallelIo = PARALLELIO_DEFAULT_BACKEND;
#endif
 public:
  assembly* readPropertyFile(const MString& fileName);
  MInt solverCount();
  IONetcdf();
  ~IONetcdf();

 private:
  assembly* m_assembly = nullptr;
  MInt m_noSolvers;
  //  std::pair<propertyIterator, propertyIterator> m_pair;
  propertyMap* m_propertyMap = nullptr;
  propertyMap* m_propertyMapLowercase = nullptr;
  zoneMap* m_zoneMap = nullptr;
  void readZones(ParallelIo* bdFile);
  void buildDefaultZone();
  void makeProperty(MProperty*, const MString&, ParallelIo* bdFile);
  MBool checkZoneConsistency();
  MBool checkPropertyConsistency();
};

#endif
