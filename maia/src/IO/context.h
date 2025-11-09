// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef CONTEXT_H
#define CONTEXT_H
////////////////////////////////////////////////////////////////////////////////
/// \file \brief This file defines the Context class used for reading properties.
///
/// \date Thu Nov 7 2002
////////////////////////////////////////////////////////////////////////////////

#include <list>
#include "INCLUDE/maiaconstants.h"
#include "UTIL/functions.h"
#include "contexttypes.h"
#include "enums.h"

class IONetcdf;
class IOToml;

/// Type traits for enum type
template <VariableType p>
struct Property2Type;
template <>
struct Property2Type<MSTRING> {
  using type = MString;
};
template <>
struct Property2Type<MINT> {
  using type = MInt;
};
template <>
struct Property2Type<MFLOAT> {
  using type = MFloat;
};
template <>
struct Property2Type<MBOOL> {
  using type = MBool;
};
template <class T>
struct TypeTraits;
template <>
struct TypeTraits<MString> {
  static const VariableType type = MSTRING;
  static const MString name() { return "MString"; }
};
template <>
struct TypeTraits<MInt> {
  static const VariableType type = MINT;
  static const MString name() { return "MInt"; }
};
template <>
struct TypeTraits<MFloat> {
  static const VariableType type = MFLOAT;
  static const MString name() { return "MFloat"; }
};
template <>
struct TypeTraits<MBool> {
  static const VariableType type = MBOOL;
  static const MString name() { return "MBool"; }
};

/** This class manages the properties that are used in the MAIA.
 * The functions and members are statical so that you can acces
 * the properties from every where within the program. The class
 * only stores the pointers to the property structs, and not the
 * structs themselves. For storing those pointers a STL multimap
 * is used with the names of the properties as keys. The class
 * has two template functions for a more convienient useage.
 * To get a basic property, simply call the getBasicProperty() function,
 * with the property name. Furthermore, you can specify a default value and
 * a the position, in case the property is an array.
 * To get a solver specific property, simply call the getSolverProperty() function,
 * with the property name and the solverId. Furthermore, you can specify a default value and
 * a the position, in case the property is an array.
 */
class Context {
 public:
  static void readPropertyFile(FileType, const MString& fileName);

  // returns the length of the property 'name' of solver 'solverId'. If solverId is skipped, the length
  // of the basic property 'name' is returned
  static MInt propertyLength(const MString& name, MInt solverId = m_noSolvers);

  /// Assert that the length of a property matches the given length
  static void assertPropertyLength(const MString& name, const MInt length, const MInt solverId = m_noSolvers);

  // returns the solverId of the property 'name' of solver 'solverId'. If solverId is skipped, the solverId
  // of the basic property 'name' is returned (only used in checkPropertyViolation)
  static MInt propertySolverId(const MString& name, MInt solverId = m_noSolvers);

  // returns true if the property 'name' of solver 'solverId' or the basic property 'name' exists.
  static MBool propertyExists(const MString& name, MInt solver = m_noSolvers);

  // returns true only if the property 'name' of solver 'solverId' exists
  static MBool solverPropertyExists(const MString& name, MInt solver);

  static void clear();

  static void init();

  static void addProperty(MProperty*);

  static void writeProperties(char* fileName);

  static void writePropertiesHumanReadable();

  static void initializationProcessFinished();

  static void communicateProperties();

  static std::list<MProperty*>* receiveProperties(MInt rank);

  static void sendProperties(MInt rank);

  static void checkPropertyViolation(MInt partnerrank, std::list<MProperty*>* prop);

  /// Return name of most recently read property file
  static MString propertyFileName() { return m_name; }

  /// Return unprocessed content of most recently read property file
  static MString propertyFileText() { return m_propertyFileText; }

  static void dump(const MString& fileName);

 private:
  static IONetcdf* m_IONetcdf;
  static IOToml* m_IOToml;
  static MString m_name;
  static MInt m_noSolvers;
  static propertyMap* m_propertyMap;
  static propertyMap* m_propertyMapLowercase;
  static zoneMap* m_zoneMap;
  static assembly* m_assembly;
  static std::pair<propertyMap::iterator, propertyMap::iterator> m_pair;
  static MString m_propertyFileOutputName;
  static MString m_propertyFileText;
  static MInt m_fileType;
  static MBool m_checkingProperties;

  //-------------Get Properties New approach--------------//
 public:
  // Returns the value of a basic property. If the basic property does not
  // exist, the default value will be returned.
  template <typename T>
  static T getBasicProperty(const MString name, const MString& nameOfCallingFunction, const T* default_value,
                            MInt pos = 0) {
    const MBool hasDefault = (default_value != nullptr);
    return getBasicPropertyOverloaded(nameOfCallingFunction, name, hasDefault, default_value, pos);
  }

  // Returns the value of a basic property. Since no default value exists in this case,
  // the requested basic property is required. If it does not exist an error message will
  // be returned.
  template <typename T>
  static T getBasicProperty(const MString name, const MString& nameOfCallingFunction, MInt pos = 0) {
    T* defaultValue = nullptr;
    return getBasicPropertyOverloaded(nameOfCallingFunction, name, false, defaultValue, pos);
  }

  // Returns the value of a solver property. If the solver property of the solver with
  // id solverId does not exist, the basic property will be returned.
  // If the top level property also does not exist, the default value will be returned.
  template <typename T>
  static T getSolverProperty(const MString name, const MInt solverId, const MString& nameOfCallingFunction,
                             const T* default_value, MInt pos = 0) {
    return getSolverPropertyOverloaded(nameOfCallingFunction, name, solverId, true, default_value, pos);
  }

  // Returns the value of a solver property. If the solver property of the solver with
  // id solverId does not exist, the basic property will be returned.
  // Since no default value exists in this case, the requested solver property is
  // required. If it does not exist an error message will be returned.
  template <typename T>
  static T getSolverProperty(const MString name, const MInt solverId, const MString& nameOfCallingFunction,
                             MInt pos = 0) {
    T* defaultValue = nullptr;
    return getSolverPropertyOverloaded(nameOfCallingFunction, name, solverId, false, defaultValue, pos);
  }

 private:
  // overloaded and implicit versions of getBasicProperty() and getSolverProperty()
  template <typename T, typename F>
  static T getBasicPropertyImplicit(const MString& nameOfCallingFunction, const MString& name, MBool has_default,
                                    const T* defaultValue, MInt position, F&& f);

  static MInt getBasicPropertyOverloaded(const MString& nameOfCallingFunction, const MString& name, MBool has_default,
                                         const MInt* default_value, MInt pos);

  static MFloat getBasicPropertyOverloaded(const MString& nameOfCallingFunction, const MString& name, MBool has_default,
                                           const MFloat* default_value, MInt pos);

  static MString getBasicPropertyOverloaded(const MString& nameOfCallingFunction, const MString& name,
                                            MBool has_default, const MString* default_value, MInt pos);

  static MBool getBasicPropertyOverloaded(const MString& nameOfCallingFunction, const MString& name, MBool has_default,
                                          const MBool* default_value, MInt pos);

  template <typename T, typename F>
  static T getSolverPropertyImplicit(const MString& nameOfCallingFunction, const MString& name, const MInt solverId,
                                     MBool has_default, const T* defaultValue, MInt position, F&& f);

  static MInt getSolverPropertyOverloaded(const MString& nameOfCallingFunction, const MString& name,
                                          const MInt solverId, MBool has_default, const MInt* default_value, MInt pos);

  static MFloat getSolverPropertyOverloaded(const MString& nameOfCallingFunction, const MString& name,
                                            const MInt solverId, MBool has_default, const MFloat* default_value,
                                            MInt pos);

  static MString getSolverPropertyOverloaded(const MString& nameOfCallingFunction, const MString& name,
                                             const MInt solverId, MBool has_default, const MString* default_value,
                                             MInt pos);

  static MBool getSolverPropertyOverloaded(const MString& nameOfCallingFunction, const MString& name,
                                           const MInt solverId, MBool has_default, const MBool* default_value,
                                           MInt pos);

  static MBool isSameValue(const MFloat v1, const MFloat v2) { return approx(v1, v2, MFloatEps); }
  static MBool isSameValue(const MBool v1, const MBool v2) { return (v1 == v2); }
  static MBool isSameValue(const MInt v1, const MInt v2) { return (v1 == v2); }
  static MBool isSameValue(const MString v1, const MString v2) { return (v1 == v2); }
};

/// Reads the number of spatial dimensions from the property file.
inline MInt read_nDim() {
  /*! \page propertyPage1
    \section nDim
    <code>MInt MAIA::nDim </code>\n
    default = <code>-1</code>\n \n
    Number of space dimensions, possible values:
    <ul>
    <li>2</li>
    <li>3</li>
    </ul>
    Keywords: <i>GLOBAL</i>
  */
  MInt nDim = -1;
  nDim = Context::getBasicProperty<MInt>("nDim", AT_, &nDim);

  if(nDim == -1) {
    /*! \page propertyPage1
      \section spaceDimensions
      <code>MInt MAIA::nDim </code>\n
      default = <code>-1</code>\n \n
      Old property for number of space dimensions,\n
      DON'T USE ANYMORE, possible values:
      <ul>
      <li>2</li>
      <li>3</li>
      </ul>
      Keywords: <i>GLOBAL</i>
    */
    nDim = Context::getBasicProperty<MInt>("spaceDimensions", AT_);
  }
  if(nDim == -1) {
    mTerm(1, AT_, "nDim / spaceDimensions property not found");
  }
  return nDim;
}

#endif
