// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "geometryproperty.h"

#include "globals.h"

using namespace std;

GeometryProperty::~GeometryProperty() {
  DEBUG("GeometryProperty::~GeometryProperty  entry", MAIA_DEBUG_ALLOCATION);

  DEBUG("GeometryProperty::~GeometryProperty  return", MAIA_DEBUG_ALLOCATION);
}

void GeometryProperty::clear() {
  TRACE();
  switch(propertyType) {
    case MINT: {
      delete[] intField;
      break;
    }
    case MFLOAT: {
      delete[] floatField;
      break;
    }
    case MSTRING: {
      delete[] stringField;
      break;
    }
    default: {
      stringstream errorMessage;
      errorMessage << "GeometryProperty::clear(): switch variable 'propertyType' with value " << propertyType
                   << " not matching any case." << endl;
      mTerm(1, AT_, errorMessage.str());
    }
  }
}
GeometryProperty::GeometryProperty() {
  DEBUG("GeometryProperty::GeometryProperty  entry", MAIA_DEBUG_ALLOCATION);
  DEBUG("GeometryProperty::GeometryProperty  return", MAIA_DEBUG_ALLOCATION);
}

VariableType GeometryProperty::type() {
  //  TRACE();

  return propertyType;
}

MInt GeometryProperty::count() {
  //  TRACE();
  return elements;
}

MString* GeometryProperty::asString() {
  //  TRACE();
  return stringField;
}

MInt* GeometryProperty::asInt() {
  //  TRACE();
  return intField;
}

MFloat* GeometryProperty::asFloat() { return floatField; }

MString* GeometryProperty::asString(MInt index) {
  //  TRACE();
  if(index + 1 > elements) {
    cerr << " GeometryProperty::asString()  Requested index out of reach! ";
    // SX8 cannont compile with exception handling
    // throw("  Requested index for property out of reach! ");
  }
  return (stringField + index);
}

MInt* GeometryProperty::asInt(MInt index) {
  //  TRACE();
  if(index + 1 > elements) {
    cerr << " GeometryProperty::asInt()  Requested index out of reach! ";
    // SX8 cannont compile with exception handling
    // throw("  Requested index for property out of reach! ");
  }
  return (intField + index);
}

MFloat* GeometryProperty::asFloat(MInt index) {
  //  TRACE();
  if(index + 1 > elements) {
    cerr << "  GeometryProperty::asFloat()  Requested index out of reach! ";
    // SX8 cannont compile with exception handling
    // throw("  Requested index for property out of reach! ");
  }
  return (floatField + index);
}
