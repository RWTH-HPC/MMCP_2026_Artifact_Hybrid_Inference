// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef GEOMETRYIOBASE_H
#define GEOMETRYIOBASE_H

#include "COMM/mpioverride.h"
#include "geometrycontexttypes.h"

class GeometryIOBase {
 public:
  virtual geometryAssembly* readPropertyFile(MString fileName) = 0;
  virtual MInt segmentCount() = 0;
  virtual void writeProperties(const MChar* /*fileName*/, geometryPropertyMap*){};
  GeometryIOBase(const MPI_Comm comm) {
    MPI_Comm_rank(comm, &m_domainId);
    MPI_Comm_size(comm, &m_noDomains);
  };
  virtual ~GeometryIOBase(){};

  MInt domainId() { return m_domainId; };
  MInt noDomains() { return m_noDomains; };

 private:
  MInt m_domainId;
  MInt m_noDomains;
};

#endif
