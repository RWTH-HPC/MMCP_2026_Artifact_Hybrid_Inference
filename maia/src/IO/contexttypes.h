// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef CONTEXTTYPES_H
#define CONTEXTTYPES_H

#include <map>
#include "INCLUDE/maiatypes.h"

class MProperty;

typedef std::multimap<MString, MProperty*> propertyMap;
typedef propertyMap::const_iterator propertyIterator;
typedef std::map<MString, MZone*> zoneMap;
typedef zoneMap::const_iterator zoneIterator;
typedef struct a {
  propertyMap* properties;
  propertyMap* propertiesLowercase;
  zoneMap* zones;
} assembly;

#endif
