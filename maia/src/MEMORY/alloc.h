// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef MAIA_ALLOC_H
#define MAIA_ALLOC_H
////////////////////////////////////////////////////////////////////////////////
/// \file \brief This file contains the allocation routines mAlloc,
/// mDeallocate, and mDealloc.
////////////////////////////////////////////////////////////////////////////////
#include <vector>
#include "INCLUDE/maiatypes.h"
#include "IO/infoout.h"
#include "collector.h"
#include "genericobject.h"
#include "globalvariables.h"
#include "list.h"
////////////////////////////////////////////////////////////////////////////////
template <class T>
class Collector; // Forward declaration
////////////////////////////////////////////////////////////////////////////////

namespace maia {
namespace alloc {

/// Debug output for mAlloc
#ifdef MAIA_EXTRA_DEBUG

inline void debug(const MString& objectName, const MString& elementsName) {
  MString allocatedObjectName = (objectName == elementsName) ? objectName : objectName + elementsName;
  typedef std::vector<GenericObject*>::iterator GenericObjIt;
  for(GenericObjIt i = g_allocatedObjects.begin(); i != g_allocatedObjects.end(); i++) {
    if((*i)->getName() == allocatedObjectName) {
      std::cerr << "Warning: object " << objectName << " has already been allocated!" << std::endl;
      m_log << "Warning: object " << objectName << " has already been allocated!" << std::endl;
    }
  }
}
inline void debug(const MString& objectName) { debug(objectName, objectName); }
inline void debug_collector(const MString& objectName, const MInt maxSize, const MString& function) {
  debug(objectName);
  m_log << "Creating collector '" << objectName << "' with " << maxSize << " elements as requested by " << function
        << "." << std::endl;
}

inline void debug_list(const MString& objectName, const MInt maxSize, const MString& function) {
  debug(objectName);
  m_log << "Creating list '" << objectName << "' with " << maxSize << " elements as requested by " << function << "."
        << std::endl;
}
#else
inline void debug(const MString&, const MString&) {}
inline void debug(const MString&) {}
inline void debug_collector(const MString&, const MInt, const MString&) {}
inline void debug_list(const MString&, const MInt, const MString&) {}
#endif


/// Find the allocated object for a given pointer
template <class T>
inline MInt findAllocatedObject(const T* const a, std::vector<GenericObject*>::iterator& hit) {
  if(a == nullptr) {
    return 0;
  }

  MInt count = 0;
  for(auto i = g_allocatedObjects.begin(); i != g_allocatedObjects.end(); i++) {
    auto* object = static_cast<GenericPointer<T>*>(*i);
    if(static_cast<T*>((*object).getObjectPointer()) == a) {
      count++;
      hit = i;
      if((*object).objectIsArray) {
        auto j = i + 1;
        while((j != g_allocatedObjects.end()) && ((*j)->getObjectId() == (*i)->getObjectId())) {
          count++;
          j++;
        }
      }
      break;
    }
  }
  return count;
}


/// Check if the given pointer belongs to an allocated object
template <class T>
inline MBool isAllocated(const T* const a) {
  std::vector<GenericObject*>::iterator hit;
  const MInt count = findAllocatedObject(a, hit);
  return (count > 0);
}


#ifdef MAIA_DEBUG_ALLOC
template <class T>
inline void checkPointer(const T* const a, const MString objectName, const MString function) {
  if(isAllocated(a)) {
    std::cerr << "Warning in mAlloc: pointer for object '" << objectName
              << "' already in list of allocated objects, call mDeallocate before reallocating! " << function
              << std::endl;
  } else if(a != nullptr) {
    std::cerr << "Warning in mAlloc: passed pointer != nullptr for object '" << objectName << "'! " << function
              << std::endl;
  }
}
#else
template <class T>
inline void checkPointer(const T* const, const MString&, const MString&) {}
#endif


template <class T>
inline void store_collector(Collector<T>* a, MLong maxSize, const MString& objectName, const MString& functionName) {
  const MLong size = (T::staticElementSize() * (maxSize + 2) * sizeof(char)) + ((maxSize + 1) * sizeof(T));
  g_allocatedObjects.push_back(new GenericPointer<Collector<T>>(a, objectName, size, functionName, false));
}

// This should just be this:
// template<class T, class... Us>  inline Collector<T>* make_collector( const Us... us ) {
//   return new Collector<T>(std::forward<Us>(us)...);
// }
// But in C++03 we need to do this (although it can be done nicer with a trait...):
template <class T>
inline Collector<T>* make_collector(const MLong size, const MInt dummy1, const MInt dummy2, const MInt dummy3,
                                    const MInt dummy4) {
  return new Collector<T>(size, dummy1, dummy2, dummy3, dummy4);
}
template <class T>
inline Collector<T>* make_collector(const MLong size, const MInt dimension, const MFloat dummy, const MInt maxNoSets) {
  return new Collector<T>(size, dimension, dummy, maxNoSets);
}
// Used to allocate:
// - DG: m_surfaces
template <class T>
inline Collector<T>* make_collector(const MLong size, const MInt dimension, const MInt dummy, const MInt dummy1) {
  return new Collector<T>(size, dimension, dummy, dummy1);
}
template <class T>
inline Collector<T>* make_collector(const MLong size, const MInt dimension, const MInt distributions,
                                    const MInt distributions1, const MInt maxNoSurfaces, const MInt dummy1) {
  return new Collector<T>(size, dimension, distributions, distributions1, maxNoSurfaces, dummy1);
}
// Used to allocate:
// - DG: m_elements
template <class T>
inline Collector<T>* make_collector(const MLong size, const MInt dimension, const MInt distributions) {
  return new Collector<T>(size, dimension, distributions);
}
template <class T>
inline Collector<T>* make_collector(const MLong size, const MInt dimension) {
  return new Collector<T>(size, dimension);
}
template <class T>
inline Collector<T>* make_collector(const MLong size) {
  return new Collector<T>(size);
}

} // namespace alloc
} // namespace maia

MLong allocatedBytes();
MLong maxAllocatedBytes();

/** \brief allocates memory for one-dimensional array 'a' of size N
 *
 * \author Lennart Schneiders
 * \date 08.12.2011
 */
template <class T>
void mAlloc(T*& a, const MLong N, const MString& objectName, MString function) {
  using namespace maia::alloc;
  maia::alloc::debug(objectName);
  maia::alloc::checkPointer(a, objectName, function);

#ifdef MAIA_ASSERT_ALLOC
  ASSERT(N > 0, "Error in mAlloc: size must be >0, is " + std::to_string(N) + " for " + objectName + " allocated in "
                    + function);
#endif

  if(N <= 0) {
#ifdef MAIA_DEBUG_ALLOC
    std::cerr << "Warning in mAlloc: size should be >0, is " << std::to_string(N) << " for " << objectName
              << " allocated in " << function << std::endl;
#endif
    return;
  }

  a = new T[N];
  g_allocatedObjects.push_back(new GenericPointer<T>(a, objectName, N * sizeof(T), function, true));
}

/** \brief allocates memory for one-dimensional array 'a' of size N
 *
 * \author Lennart Schneiders
 * \date 08.12.2011
 */
template <class T>
void mAlloc(T*& a, const MLong N, const MString& objectName, const T initialValue, MString function) {
  mAlloc(a, N, objectName, function);
#ifndef SKIP_INITIAL_VALUE_ALLOC
  for(MInt i = 0; i < N; i++)
    a[i] = initialValue;
#endif
}

/** \brief allocates memory for two-dimensional array 'a' of size NxM
 *
 * \author Lennart Schneiders
 * \date 08.12.2011
 */
template <class T>
void mAlloc(T**& a, const MLong N, const MLong M, const MString& objectName, const T initialValue, MString function) {
  using namespace maia::alloc;
  maia::alloc::debug(objectName);
  maia::alloc::checkPointer(a, objectName, function);

#ifdef MAIA_ASSERT_ALLOC
  ASSERT(N > 0 && M > 0, "Error in mAlloc: sizes must be >0, but they are " + std::to_string(N) + " and "
                             + std::to_string(M) + " for " + objectName + " allocated in " + function);
#endif

  if(N <= 0 || M <= 0) {
#ifdef MAIA_DEBUG_ALLOC
    std::cerr << "Error in mAlloc: sizes should be >0, but they are " << std::to_string(N) << " and "
              << std::to_string(M) << " for " << objectName << " allocated in " << function << std::endl;
#endif
    return;
  }

  a = new T*[N];
  g_allocatedObjects.push_back(new GenericPointer<T*>(a, objectName, N * sizeof(T*), function, true));

  const MString elementsName = "_elements";
  maia::alloc::debug(objectName, elementsName);

  T* dummy = new T[N * M];
  g_allocatedObjects.push_back(
      new GenericPointer<T>(dummy, objectName + elementsName, N * M * sizeof(T), function, true, false));

  for(MInt i = 0; i < N; i++) {
    a[i] = (T*)(&(dummy[M * i]));
  }

#ifndef SKIP_INITIAL_VALUE_ALLOC
  for(MInt i = 0; i < N * M; i++)
    dummy[i] = initialValue;
#endif
}


/** \brief allocates memory for two-dimensional array 'a' of size NxM
 *
 * \author Lennart Schneiders
 * \date 08.12.2011
 */
template <class T>
void mAlloc(T**& a, const MLong N, const MLong M, const MString& objectName, MString function) {
  using namespace maia::alloc;
  maia::alloc::debug(objectName);
  maia::alloc::checkPointer(a, objectName, function);

#ifdef MAIA_ASSERT_ALLOC
  ASSERT(N > 0 && M > 0, "Error in mAlloc: sizes must be >0, but they are " + std::to_string(N) + " and "
                             + std::to_string(M) + " for " + objectName + " allocated in " + function);
#endif

  if(N <= 0 || M <= 0) {
#ifdef MAIA_DEBUG_ALLOC
    std::cerr << "Error in mAlloc: sizes should be >0, but they are " << std::to_string(N) << " and "
              << std::to_string(M) << " for " << objectName << " allocated in " << function << std::endl;
#endif
    return;
  }

  a = new T*[N];

  g_allocatedObjects.push_back(new GenericPointer<T*>(a, objectName, N * sizeof(T*), function, true));

  const MString elementsName = "[0]";
  maia::alloc::debug(objectName, elementsName);

  a[0] = new T[N * M];
  g_allocatedObjects.push_back(
      new GenericPointer<T>(a[0], objectName + elementsName, N * M * sizeof(T), function, true, false));

  for(MInt i = 1; i < N; i++) {
    a[i] = (T*)(&(a[0][M * i]));
  }
}


/** \brief allocates memory for two-dimensional array 'a' of size N x arraySizes[0:N-1]
 *
 * \author Lennart Schneiders
 * \date 08.12.2011
 */
template <class T>
void mAlloc(T**& a, const MLong N, const MInt* arraySizes, const MString& objectName, MString function) {
  using namespace maia::alloc;
  maia::alloc::debug(objectName);
  maia::alloc::checkPointer(a, objectName, function);

  size_t size = 0;
  for(MInt k = 0; k < N; k++) {
    size += arraySizes[k];
  }

#ifdef MAIA_ASSERT_ALLOC
  ASSERT(N > 0 && size > 0, "Error in mAlloc: sizes must be >0, but they are " + std::to_string(N) + " and "
                                + std::to_string(size) + " for " + objectName + " allocated in " + function);
#endif

  if(N <= 0 || size <= 0) {
#ifdef MAIA_DEBUG_ALLOC
    std::cerr << "Error in mAlloc: sizes should be >0, but they are " << std::to_string(N) << " and "
              << std::to_string(size) << " for " << objectName << " allocated in " << function << std::endl;
#endif
    return;
  }

  a = new T*[N];
  g_allocatedObjects.push_back(new GenericPointer<T*>(a, objectName, N * sizeof(T*), function, true));

  const MString elementsName = "_elements";
  maia::alloc::debug(objectName, elementsName);

  T* dummy = new T[size];
  g_allocatedObjects.push_back(
      new GenericPointer<T>(dummy, objectName + elementsName, size * sizeof(T), function, true, false));

  size_t cnt = 0;
  for(MInt i = 0; i < N; i++) {
    a[i] = (T*)(&(dummy[cnt]));
    cnt += arraySizes[i];
  }
}


/** \brief allocates memory for two-dimensional array 'a'
 *   of size N x (factor*arraySizes[0:N-1])
 *
 * \author Lennart Schneiders
 * \date 08.12.2011
 */
template <class T>
void mAlloc(T**& a, const MLong N, const MInt* const arraySizes, MInt factor, const MString& objectName,
            MString function) {
  using namespace maia::alloc;
  maia::alloc::debug(objectName);
  maia::alloc::checkPointer(a, objectName, function);

  size_t size = 0;
  for(MInt k = 0; k < N; k++) {
    size += factor * arraySizes[k];
  }

#ifdef MAIA_ASSERT_ALLOC
  ASSERT(N > 0 && size > 0, "Error in mAlloc: sizes must be >0, but they are " + std::to_string(N) + " and "
                                + std::to_string(size) + " for " + objectName + " allocated in " + function);
#endif

  if(N <= 0 || size <= 0) {
#ifdef MAIA_DEBUG_ALLOC
    std::cerr << "Error in mAlloc: sizes should be >0, but they are " << std::to_string(N) << " and "
              << std::to_string(size) << " for " << objectName << " allocated in " << function << std::endl;
#endif
    return;
  }

  a = new T*[N];
  g_allocatedObjects.push_back(new GenericPointer<T*>(a, objectName, N * sizeof(T*), function, true));

  const MString elementsName = "_elements";
  maia::alloc::debug(objectName, elementsName);

  T* dummy = new T[size];
  g_allocatedObjects.push_back(
      new GenericPointer<T>(dummy, objectName + elementsName, size * sizeof(T), function, true, false));

  size_t cnt = 0;
  for(MInt i = 0; i < N; i++) {
    a[i] = (T*)(&(dummy[cnt]));
    cnt += factor * arraySizes[i];
  }
}


/** \brief allocates memory for List container
 *
 * \author Lennart Schneiders
 * \date 08.12.2011
 */
template <class T>
void mAlloc(List<T>*& a, const MLong N, const MString& objectName, const T initialValue, MString function) {
  using namespace maia::alloc;
  maia::alloc::debug_list(objectName, N, function);
  maia::alloc::checkPointer(a, objectName, function);

#ifdef MAIA_ASSERT_ALLOC
  ASSERT(N > 0, "Error in mAlloc: size must be >0, is " + std::to_string(N) + " for " + objectName + " allocated in "
                    + function);
#endif

  if(N <= 0) {
#ifdef MAIA_DEBUG_ALLOC
    std::cerr << "Warning in mAlloc: size should be >0, is " << std::to_string(N) << " for " << objectName
              << " allocated in " << function << std::endl;
#endif
    return;
  }

  a = new List<T>(N);

  g_allocatedObjects.push_back(new GenericPointer<List<T>>(a, objectName, N * sizeof(T), function, false));

#ifndef SKIP_INITIAL_VALUE_ALLOC
  for(MInt i = 0; i < N; i++)
    a->a[i] = initialValue;
#endif
}


/** \brief allocates memory for Collector
 *
 * \author Lennart Schneiders
 * \date 08.12.2011
 */
template <class T>
void mAlloc(Collector<T>*& a, MLong maxSize, MInt dummy1, MInt dummy2, MInt dummy3, MInt dummy4,
            const MString& objectName, MString function) {
  maia::alloc::checkPointer(a, objectName, function);
  maia::alloc::debug_collector(objectName, maxSize, function);
  a = maia::alloc::make_collector<T>(maxSize, dummy1, dummy2, dummy3, dummy4);
  maia::alloc::store_collector(a, maxSize, objectName, function);
}


/** \brief allocates memory for Collector
 *
 * \author Lennart Schneiders
 * \date 08.12.2011
 */
template <class T>
void mAlloc(Collector<T>*& a, MLong maxSize, MInt dimension, MFloat dummy, MInt maxNoSets, const MString& objectName,
            MString function) {
  maia::alloc::checkPointer(a, objectName, function);
  maia::alloc::debug_collector(objectName, maxSize, function);
  a = maia::alloc::make_collector<T>(maxSize, dimension, dummy, maxNoSets);
  maia::alloc::store_collector(a, maxSize, objectName, function);
}


/** \brief allocates memory for Collector
 *
 * \author Lennart Schneiders
 * \date 08.12.2011
 *
 * Used to allocate:
 *  - DG: m_surfaces
 */
template <class T>
void mAlloc(Collector<T>*& a, MLong maxSize, MInt dimension, MInt dummy, MInt dummy1, const MString& objectName,
            MString function) {
  maia::alloc::checkPointer(a, objectName, function);
  maia::alloc::debug_collector(objectName, maxSize, function);
  a = maia::alloc::make_collector<T>(maxSize, dimension, dummy, dummy1);
  maia::alloc::store_collector(a, maxSize, objectName, function);
}


/** \brief Allocates memory for Collector
 *
 * Used to allocate:
 *  - FvBndryCell
 *
 * \author Lennart Schneiders
 * \date 08.12.2011
 */
template <class T>
void mAlloc(Collector<T>*& a, MLong maxSize, MInt dimension, MInt distributions, MInt distributions1,
            MInt maxNoSurfaces, MInt dummy1, const MString& objectName, MString function) {
  maia::alloc::checkPointer(a, objectName, function);
  maia::alloc::debug_collector(objectName, maxSize, function);
  a = maia::alloc::make_collector<T>(maxSize, dimension, distributions, distributions1, maxNoSurfaces, dummy1);
  maia::alloc::store_collector(a, maxSize, objectName, function);
}


/** \brief allocates memory for Collector
 *
 * \author Lennart Schneiders
 * \date 08.12.2011
 *
 * Used to allocate:
 *  - DG: m_elements
 */
template <class T>
void mAlloc(Collector<T>*& a, MLong maxSize, MInt dimension, MInt distributions, const MString& objectName,
            MString function) {
  maia::alloc::checkPointer(a, objectName, function);
  maia::alloc::debug_collector(objectName, maxSize, function);
  a = maia::alloc::make_collector<T>(maxSize, dimension, distributions);
  maia::alloc::store_collector(a, maxSize, objectName, function);
}


/** \brief allocates memory for Collector
 *
 * \author Lennart Schneiders
 * \date 08.12.2011
 */
template <class T>
void mAlloc(Collector<T>*& a, MLong maxSize, MInt dimension, const MString& objectName, MString function) {
  maia::alloc::checkPointer(a, objectName, function);
  maia::alloc::debug_collector(objectName, maxSize, function);
  a = maia::alloc::make_collector<T>(maxSize, dimension);
  maia::alloc::store_collector(a, maxSize, objectName, function);
}


/** \brief allocates memory for Collector
 *
 * \author Lennart Schneiders
 * \date 08.12.2011
 */
template <class T>
void mAlloc(Collector<T>*& a, MLong maxSize, const MString& objectName, MString function) {
  maia::alloc::checkPointer(a, objectName, function);
  maia::alloc::debug_collector(objectName, maxSize, function);
  a = maia::alloc::make_collector<T>(maxSize);
  maia::alloc::store_collector(a, maxSize, objectName, function);
}


/** \brief deallocates the memory previously allocated for element 'a'
 *
 * \author Lennart Schneiders
 * \date 09.12.2011
 */
template <class T>
MBool mDeallocate(T*& a) {
  using namespace maia::alloc;
  if(a == nullptr) {
#ifdef MAIA_EXTRA_DEBUG
    m_log << "cannot deallocate zero pointer " << a << std::endl;
#endif
    return false;
  }

  std::vector<GenericObject*>::iterator hit;
  const MInt count = findAllocatedObject(a, hit);

#ifdef MAIA_EXTRA_DEBUG
  if(count > 0)
    m_log << "Deallocating " << (*hit)->getName() << "."
          << " " << &a << " " << count << std::endl;
  m_log.flush();
  if(count > 0)
    std::cerr << "Deallocating " << (*hit)->getName() << "."
              << " " << &a << " " << count << std::endl;
#endif

  if(count <= 0) {
#ifndef NDEBUG
    std::cerr << "mDeallocate error: " << count << " " << a << std::endl;
    if(count == 0) {
      std::cerr << "Pointer should have been reset to nullptr!" << std::endl;
    }
#endif
    a = nullptr;
    return false;
  }


  for(MInt k = count; k--;) {
    delete *(hit + k);
  }

  if(count == 1)
    g_allocatedObjects.erase(hit);
  else if(count > 1)
    g_allocatedObjects.erase(hit, hit + count);

  if(count > 0) {
    a = nullptr;
  }

  return (count > 0);
}

/** \brief Deallocates all memory allocated previously by mAlloc(...)
 *
 * \author Lennart Schneiders
 * \date 08.12.2011
 */
void mDealloc();

/**\brief Prints the name of all objects which are currently allocated via mAlloc
 *
 * \author Julian Vorspohl
 * \date 20.06.2021
 */
inline void printAllocatedObjects() {
  for(auto&& obj : maia::alloc::g_allocatedObjects) {
#ifdef MAIA_EXTRA_DEBUG
    std::cout << "Name " << obj->getName() << " Function " << obj->getCallingFunction() << std::endl;
#else
    std::cout << "Name " << obj->getName() << std::endl;
#endif
  }
  std::cout << "Total number of allocated objects: " << maia::alloc::g_allocatedObjects.size() << std::endl;
}

/**\brief Prints the name of all objects which are currently allocated via mAlloc and contain the given string
 *
 * \author Julian Vorspohl
 * \date 20.06.2021
 */
#ifdef MAIA_EXTRA_DEBUG
inline void printAllocatedObjects(const MString& filterString) {
  for(auto&& obj : maia::alloc::g_allocatedObjects) {
    const MString fName = obj->getCallingFunction();
    const MBool containsFilterString = (fName.find(filterString) != std::string::npos);
    if(containsFilterString) {
      std::cout << filterString << ": Memory for " << obj->getName() << " from function " << fName
                << " should have been deallocated by now!" << std::endl;
    }
  }
}
#else
inline void printAllocatedObjects(const MString& /*filterString*/) {}
#endif

#endif
