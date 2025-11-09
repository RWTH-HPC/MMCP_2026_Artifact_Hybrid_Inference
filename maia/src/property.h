// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef PROPERTY_H
#define PROPERTY_H

#include "INCLUDE/maiatypes.h"
#include "enums.h"

class MProperty {
 public:
  MProperty();
  MProperty(MInt size, MString propname, MInt Solver, MString* value);
  MProperty(MInt size, MString propname, MInt Solver, MFloat* value);
  MProperty(MInt size, MString propname, MInt Solver, MInt* value);
  MProperty(MInt size, MString propname, MInt Solver, MBool* value);
  ~MProperty();
  MInt count();
  MInt doesExist();
  MString* asString();
  MInt* asInt();
  MFloat* asFloat();
  MBool* asBool();
  VariableType type();
  void clear();

  MString* asString(MInt index);
  MInt* asInt(MInt index);
  MFloat* asFloat(MInt index);
  MBool* asBool(MInt index);
  VariableType propertyType;
  MInt elements = -1;
  MString name;
  MInt noAccesses = 0;
  MInt noOldAccesses = 0;
  MInt solverId = -1;
  MString* stringField = nullptr;
  MFloat* floatField = nullptr;
  MInt* intField = nullptr;
  MBool* boolField = nullptr;
};


#endif
