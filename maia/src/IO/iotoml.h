// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef IOTOML_H_
#define IOTOML_H_

#include <memory>
#include "INCLUDE/maiatypes.h"
#include "contexttypes.h"
#include "enums.h"


// Forward declaration to avoid header include
namespace maia {
namespace io {
namespace toml {
class Property;
}
} // namespace io
} // namespace maia

class IOToml {
 public:
  assembly* readPropertyFile(const MString& fileName);
  MInt solverCount();
  MString rawText() const { return m_rawText; }

 private:
  assembly* m_assembly = nullptr;
  MInt m_noSolvers = -1;
  std::pair<propertyIterator, propertyIterator> m_pair;
  propertyMap* m_propertyMap = nullptr;
  propertyMap* m_propertyMapLowercase = nullptr;
  zoneMap* m_zoneMap = nullptr;
  void buildDefaultZone();
  void makeProperty(MProperty*, const maia::io::toml::Property& prop);
  MBool checkZoneConsistency();
  MBool checkPropertyConsistency();
  MString m_rawText = "";
};

#endif // #ifdef IOTOML_H_
