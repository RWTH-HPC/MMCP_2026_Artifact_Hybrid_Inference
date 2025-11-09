// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef GEOMETRYCONTEXTTYPES_H
#define GEOMETRYCONTEXTTYPES_H

#include <map>
#include "geometryproperty.h"

//! define array structures
typedef struct y {
  MString name;
  MInt* segments;
  MInt noSegments;
} Body;
typedef std::multimap<MString, GeometryProperty*> geometryPropertyMap;
typedef geometryPropertyMap::const_iterator geometryPropertyIterator;
typedef std::map<MString, Body*> bodyMap;
typedef bodyMap::const_iterator bodyIterator;
typedef struct b {
  geometryPropertyMap* geometryProperties;
  bodyMap* bodies;
} geometryAssembly;


#endif
