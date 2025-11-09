// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef GEOMETRYPROPERTY_H
#define GEOMETRYPROPERTY_H

#include "INCLUDE/maiatypes.h"
#include "enums.h"

class GeometryProperty {
 public:
  GeometryProperty();
  ~GeometryProperty();
  MInt count();
  MString* asString();
  MInt* asInt();
  MFloat* asFloat();
  VariableType type();
  void clear();

  MString* asString(MInt index);
  MInt* asInt(MInt index);
  MFloat* asFloat(MInt index);
  VariableType propertyType;
  MInt elements;
  MString name;
  MInt segmentId;
  MString* stringField;
  MFloat* floatField;
  MInt* intField;
};


#endif
