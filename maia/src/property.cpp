// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "property.h"

#include "UTIL/debug.h"
#include "UTIL/functions.h"
#include "UTIL/timer.h"

using namespace std;

MProperty::~MProperty() { clear(); }

void MProperty::clear() {
  TRACE();
  switch(propertyType) {
    case MINT: {
      delete[] intField;

      // Int properties can be requested as float, for which a copy is created. Thus also check the
      // float field if it was allocated.
      if(floatField != nullptr) {
        delete[] floatField;
      }

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
    case MBOOL: {
      delete[] boolField;
      break;
    }
    default: {
      stringstream errorMessage;
      errorMessage << "MProperty::clear(): switch variable 'propertyType' with value " << propertyType
                   << " not matching any case." << endl;
      mTerm(1, AT_, errorMessage.str());
    }
  }
}

MProperty::MProperty() {
  //  TRACE();
  noAccesses = 0;
  noOldAccesses = 0;
}

MProperty::MProperty(MInt size, MString propname, MInt Solver, MString* value) {
  //  TRACE();
  propertyType = MSTRING;
  elements = size;
  name = propname;
  noAccesses = 0;
  noOldAccesses = 0;
  solverId = Solver;
  stringField = new MString[size];
  for(MInt i = 0; i < size; i++)
    stringField[i] = value[i];
}

MProperty::MProperty(MInt size, MString propname, MInt Solver, MFloat* value) {
  //  TRACE();
  propertyType = MFLOAT;
  elements = size;
  name = propname;
  noAccesses = 0;
  noOldAccesses = 0;
  solverId = Solver;
  floatField = new MFloat[size];
  for(MInt i = 0; i < size; i++)
    floatField[i] = value[i];
}

MProperty::MProperty(MInt size, MString propname, MInt Solver, MInt* value) {
  //  TRACE();
  propertyType = MINT;
  elements = size;
  name = propname;
  noAccesses = 0;
  noOldAccesses = 0;
  solverId = Solver;
  intField = new MInt[size];
  for(MInt i = 0; i < size; i++)
    intField[i] = value[i];
}

MProperty::MProperty(MInt size, MString propname, MInt Solver, MBool* value) {
  //  TRACE();
  propertyType = MBOOL;
  elements = size;
  name = propname;
  noAccesses = 0;
  noOldAccesses = 0;
  solverId = Solver;
  boolField = new MBool[size];
  for(MInt i = 0; i < size; i++)
    boolField[i] = value[i];
}

VariableType MProperty::type() { return propertyType; }

MInt MProperty::count() { return elements; }

MString* MProperty::asString() {
  if(propertyType != MSTRING) {
    stringstream errorMessage;
    errorMessage << "Property " << name << " is requested as a String, but it is a ";
    switch(propertyType) {
      case MINT: {
        errorMessage << "INT!!";
        break;
      }
      case MFLOAT: {
        errorMessage << "FLOAT!!";
        break;
      }
      case MSTRING: {
        errorMessage << "STRING!!";
        break;
      }
      case MBOOL: {
        errorMessage << "BOOL!!";
        break;
      }
      default: {
        mTerm(1, AT_, "MProperty::asString(): switch variable 'propertyType' not matching any case");
      }
    }
    mTerm(1, AT_, errorMessage.str());
  }
  return stringField;
}

MBool* MProperty::asBool() {
  if(propertyType != MBOOL) {
    stringstream errorMessage;
    errorMessage << "Property " << name << " is requested as a Bool, but it is a ";
    switch(propertyType) {
      case MINT: {
        errorMessage << "INT!!";
        break;
      }
      case MFLOAT: {
        errorMessage << "FLOAT!!";
        break;
      }
      case MSTRING: {
        errorMessage << "STRING!!";
        break;
      }
      case MBOOL: {
        errorMessage << "BOOL!!";
        break;
      }
      default: {
        mTerm(1, AT_, "MProperty::asBool(): switch variable 'propertyType' not matching any case");
      }
    }
    mTerm(1, AT_, errorMessage.str());
  }
  return boolField;
}

MInt* MProperty::asInt() {
  if(propertyType != MINT) {
    stringstream errorMessage;
    errorMessage << "Property " << name << " is requested as a Int, but it is a ";
    switch(propertyType) {
      case MINT: {
        errorMessage << "INT!!";
        break;
      }
      case MFLOAT: {
        errorMessage << "FLOAT!!";
        break;
      }
      case MSTRING: {
        errorMessage << "STRING!!";
        break;
      }
      case MBOOL: {
        errorMessage << "BOOL!!";
        break;
      }
      default: {
        mTerm(1, AT_, "MProperty::asInt(): switch variable 'propertyType' not matching any case");
      }
    }
    mTerm(1, AT_, errorMessage.str());
  }
  return intField;
}

MFloat* MProperty::asFloat() {
  if(propertyType != MFLOAT && propertyType != MINT) {
    stringstream errorMessage;
    errorMessage << "Property " << name << " is requested as a Float, but it is a ";
    switch(propertyType) {
      case MINT: {
        errorMessage << "INT!!";
        break;
      }
      case MFLOAT: {
        errorMessage << "FLOAT!!";
        break;
      }
      case MSTRING: {
        errorMessage << "STRING!!";
        break;
      }
      case MBOOL: {
        errorMessage << "BOOL!!";
        break;
      }
      default: {
        mTerm(1, AT_, "MProperty::asFloat(): switch variable 'propertyType' not matching any case");
      }
    }
    mTerm(1, AT_, errorMessage.str());
  }

  // If int property is requested as float, create copy of int properties first
  if(propertyType == MINT && floatField == nullptr) {
    floatField = new MFloat[elements];
    copy_n(intField, elements, floatField);
  }

  return floatField;
}

MString* MProperty::asString(MInt index) {
  if(propertyType != MSTRING) {
    stringstream errorMessage;
    errorMessage << "Property " << name << " is requested as a String, but it is a ";
    switch(propertyType) {
      case MINT: {
        errorMessage << "INT!!";
        break;
      }
      case MFLOAT: {
        errorMessage << "FLOAT!!";
        break;
      }
      case MSTRING: {
        errorMessage << "STRING!!";
        break;
      }
      case MBOOL: {
        errorMessage << "BOOL!!";
        break;
      }
      default: {
        mTerm(1, AT_, "MProperty::asString(): switch variable 'propertyType' not matching any case");
      }
    }
    mTerm(1, AT_, errorMessage.str());
  }
  if(index + 1 > elements) {
    stringstream errorMessage;
    errorMessage << " MProperty::asString() for property " << name << " is requested index out of range! It has "
                 << elements << " ,but it is asked for the " << index + 1 << " element! " << endl;
    mTerm(1, AT_, errorMessage.str());
  }
  return (stringField + index);
}

MBool* MProperty::asBool(MInt index) {
  if(propertyType != MBOOL) {
    stringstream errorMessage;
    errorMessage << "Property " << name << " is requested as a Bool, but it is a ";
    switch(propertyType) {
      case MINT: {
        errorMessage << "INT!!";
        break;
      }
      case MFLOAT: {
        errorMessage << "FLOAT!!";
        break;
      }
      case MSTRING: {
        errorMessage << "STRING!!";
        break;
      }
      case MBOOL: {
        errorMessage << "BOOL!!";
        break;
      }
      default: {
        mTerm(1, AT_, "MProperty::asBool(): switch variable 'propertyType' not matching any case");
      }
    }
    mTerm(1, AT_, errorMessage.str());
  }
  if(index + 1 > elements) {
    stringstream errorMessage;
    errorMessage << " MProperty::asBool() for property " << name << " is requested index out of range! It has "
                 << elements << " ,but it is asked for the " << index + 1 << " element! ";
    mTerm(1, AT_, errorMessage.str());
  }
  //  cout << " Using int property : " << name << endl;
  //  cout << " 1st value is         : " << *(intField+index) << endl;
  return (boolField + index);
}
MInt* MProperty::asInt(MInt index) {
  if(propertyType != MINT) {
    stringstream errorMessage;
    errorMessage << "Property " << name << " is requested as a Int, but it is a ";
    switch(propertyType) {
      case MINT: {
        errorMessage << "INT!!";
        break;
      }
      case MFLOAT: {
        errorMessage << "FLOAT!!";
        break;
      }
      case MSTRING: {
        errorMessage << "STRING!!";
        break;
      }
      case MBOOL: {
        errorMessage << "BOOL!!";
        break;
      }
      default: {
        mTerm(1, AT_, "MProperty::asInt(): switch variable 'propertyType' not matching any case");
      }
    }
    mTerm(1, AT_, errorMessage.str());
  }
  if(index + 1 > elements) {
    stringstream errorMessage;
    errorMessage << " MProperty::asInt() for property " << name << " is requested index out of range! It has "
                 << elements << " ,but it is asked for the " << index + 1 << " element! ";
    mTerm(1, AT_, errorMessage.str());
  }
  //  cout << " Using int property : " << name << endl;
  //  cout << " 1st value is         : " << *(intField+index) << endl;
  return (intField + index);
}

MFloat* MProperty::asFloat(MInt index) {
  // Check if property really is a float (or int, as getting int properties as floats is allowed)
  if(propertyType != MFLOAT && propertyType != MINT) {
    stringstream errorMessage;
    errorMessage << "Property " << name << " is requested as a float, but it is a ";
    switch(propertyType) {
      case MINT: {
        errorMessage << "INT!!";
        break;
      }
      case MFLOAT: {
        errorMessage << "FLOAT!!";
        break;
      }
      case MSTRING: {
        errorMessage << "STRING!!";
        break;
      }
      case MBOOL: {
        errorMessage << "BOOL!!";
        break;
      }
      default: {
        mTerm(1, AT_, "MProperty::asInt(): switch variable 'propertyType' not matching any case");
      }
    }
    mTerm(1, AT_, errorMessage.str());
  }

  // Check if size is OK
  if(index + 1 > elements) {
    stringstream errorMessage;
    errorMessage << " MProperty::asFloat() for property " << name << " is requested index out of range! It has "
                 << elements << " ,but it is asked for the " << index + 1 << " element! ";
    mTerm(1, AT_, errorMessage.str());
  }

  // If int property is requested as float, create copy of int properties first
  if(propertyType == MINT && floatField == nullptr) {
    floatField = new MFloat[elements];
    copy_n(intField, elements, floatField);
  }

  return (floatField + index);
}
