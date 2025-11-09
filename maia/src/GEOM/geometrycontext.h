// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef GEOMETRYCONTEXT_H
#define GEOMETRYCONTEXT_H

#include <iostream>
#include "COMM/mpioverride.h"
#include "geometrycontexttypes.h"

class GeometryIOBase;

class GeometryContext {
 public:
  void readPropertyFile(FileType, const MChar* fileName);

  GeometryProperty* getProperty(const MString& name, MInt segment);
  MInt noPropertySegments(const MString& name) {
    if(m_geometryPropertyMap->find(name) == m_geometryPropertyMap->end()) {
      std::cerr << "GeometryContext::noPropertySegments " << name << " not found " << std::endl;
      return 0;
    }
    auto range = m_geometryPropertyMap->equal_range(name);

    MInt no = 0;
    for(auto& i = range.first; i != range.second; i++) {
      no++;
    }
    return no;
  }
  MBool propertyExists(MString name, MInt solver);

  bodyMap getBodies();

  void clear();

  void init();

  void addProperty(GeometryProperty*);

  void writeProperties(MChar* fileName);

  MInt getNoSegments() const { return m_noSegments; };

  ~GeometryContext();
  GeometryContext(const MPI_Comm comm) : m_mpiComm(comm){};
  MPI_Comm mpiComm() const { return m_mpiComm; }

 private:
  GeometryIOBase* m_geometryIOBase;
  MString m_name;
  MInt m_noSegments;
  MInt m_noBodies;
  geometryPropertyMap* m_geometryPropertyMap;
  bodyMap* m_bodyMap;
  geometryAssembly* m_geometryAssembly;
  std::pair<geometryPropertyMap::iterator, geometryPropertyMap::iterator> m_pair;
  const MPI_Comm m_mpiComm = MPI_COMM_NULL;
};

#endif
