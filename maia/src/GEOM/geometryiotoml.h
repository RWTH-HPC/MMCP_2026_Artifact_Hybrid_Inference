// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef GEOMETRYIOTOML_H_
#define GEOMETRYIOTOML_H_

#include <memory>
#include <vector>
#include "COMM/mpioverride.h"
#include "geometryiobase.h"


// Forward declaration to avoid header include
namespace maia {
namespace io {
namespace toml {
class Property;
}
} // namespace io
} // namespace maia

class GeometryIOToml : public GeometryIOBase {
 public:
  GeometryIOToml(const MPI_Comm comm) : GeometryIOBase(comm), m_mpiComm(comm) {}
  geometryAssembly* readPropertyFile(MString fileName);
  void readPropertyFileOldIOMethod(const std::vector<maia::io::toml::Property>& properties);
  void readPropertyFileNewIOMethod(const std::vector<maia::io::toml::Property>& properties);

  MInt segmentCount();
  MPI_Comm mpiComm() const { return m_mpiComm; }
  MString rawText() const { return m_rawText; }

 private:
  geometryAssembly* m_geometryAssembly;
  MInt m_noSegments;
  MInt m_noBodies;
  std::pair<geometryPropertyIterator, geometryPropertyIterator> m_pair;
  geometryPropertyMap* m_geometryPropertyMap;
  bodyMap* m_bodyMap;
  MBool m_newIOMethod;
  MString m_rawText = "";

  const MPI_Comm m_mpiComm = MPI_COMM_NULL;
  void readBodiesOldIOMethod(const std::vector<maia::io::toml::Property>& properties);
  void readBodiesNewIOMethod(const std::vector<maia::io::toml::Property>& properties);

  void buildDefaultBody();
  void makeProperty(GeometryProperty*, const maia::io::toml::Property& prop);
  MBool checkBodyConsistency();
  MBool checkGeometryPropertyConsistency();
};

#endif // #ifndef GEOMETRYIOTOML_H_
