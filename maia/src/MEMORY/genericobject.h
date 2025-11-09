// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef MAIA_GENERIC_OBJECT_H_
#define MAIA_GENERIC_OBJECT_H_

#include "IO/infoout.h"
#include "UTIL/functions.h"
#include "globalvariables.h"

#ifdef MAIA_EXTRA_DEBUG
#include <typeinfo>
#endif

/** \brief Prints currently allocated memory
 *
 * \author Lennart Schneiders
 * \date 08.12.2011
 */
void printAllocatedMemory(const MLong oldAllocatedBytes, const MString& solverName, const MPI_Comm comm);

/** \brief Returns memory size in KB, MB or GB
 *
 * \author Lennart Schneiders
 * \date 09.12.2011
 */
MString getMemorySize(MLong noBytes);

/** \brief class containing a generic object
 *
 * \author Lennart Schneiders
 * \date 08.12.2011
 */
class GenericObject {
 public:
  MString objectName;
  MString& getName() { return objectName; }
#ifdef MAIA_EXTRA_DEBUG
  MString callingFunction;
  MString& getCallingFunction() { return callingFunction; }
#endif
  static MInt objectCounter;
  MInt objectId;
  MInt& getObjectId() { return objectId; }
  MLong elementSize;
  MLong& getElementSize() { return elementSize; }
  virtual ~GenericObject() {}
  MBool operator<(GenericObject otherObject) { return (elementSize < otherObject.elementSize); }
};

/** \brief class containing a generic pointer
 *
 * \author Lennart Schneiders
 * \date 08.12.2011
 */
template <typename T>
class GenericPointer : public GenericObject {
 private:
  T* objectPointer;

 public:
  MBool objectIsArray;
  GenericPointer(T* const& t, MString oName, MLong size, MString function, MBool isArray, MBool isNewObject = true);
  T* getObjectPointer() { return objectPointer; }
  ~GenericPointer();
};

template <class T>
GenericPointer<T>::GenericPointer(T* const& t, MString oName, MLong size,
                                  MString
#ifdef MAIA_EXTRA_DEBUG
                                      function
#endif
                                  ,
                                  MBool isArray, MBool isNewObject)
  : objectPointer(t), objectIsArray(isArray) {
  using namespace maia::alloc;
  objectName = oName;
#ifdef MAIA_EXTRA_DEBUG
  callingFunction = function;
#endif
  elementSize = size;
  if(isNewObject) {
    objectId = objectCounter;
    objectCounter++;
  } else
    objectId = objectCounter - 1;
  g_allocatedBytes += size;
  g_maxAllocatedBytes = mMax(g_allocatedBytes, g_maxAllocatedBytes);
}

template <class T>
GenericPointer<T>::~GenericPointer() {
  maia::alloc::g_allocatedBytes -= elementSize;
#ifdef MAIA_EXTRA_DEBUG
  m_log.precision(2);
  m_log << "deleting object " << objectName << " (size: " << getMemorySize(elementSize)
        << ", address: " << objectPointer << ", type: " << typeid(objectPointer).name()
        << ", isArray: " << objectIsArray << ", called by " << callingFunction << ")...";
  m_log.flush();
#endif
  if(objectPointer) {
    if(objectIsArray)
      delete[] objectPointer;
    else
      delete objectPointer;
    objectPointer = nullptr;
  }
#ifdef MAIA_EXTRA_DEBUG
  m_log << "ok" << std::endl;
#endif
}
#endif
