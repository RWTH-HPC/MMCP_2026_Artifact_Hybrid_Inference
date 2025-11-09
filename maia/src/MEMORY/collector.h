// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef COLLECTOR_H
#define COLLECTOR_H

#include "INCLUDE/maiatypes.h"
#include "UTIL/functions.h"

template <typename T>
class Collector {
 public:
  // Note: The "dummy" variables here serve ABSOLUTELY NO OTHER PURPOSE than to differentiate the different constructors
  // (i.e. give them unique signatures so proper overloading can take place)
  Collector(MLong maxSize, MInt dimension, MFloat dummy, MInt maxNoSets); // G cells
  Collector(MLong maxSize, MInt dimension, MInt distributions, MInt distributions1);
  Collector(MLong maxSize, MInt dimension, MInt distributions, MInt distributions1, MInt maxNoSurfaces, MInt dummy1);
  Collector(MLong maxSize, MInt dimension, MInt distributions, MInt dummy, MInt dummy2);
  Collector(MLong maxSize, MInt dimension, MInt distributions);
  Collector(MLong maxSize, MInt dimension);
  Collector(MLong maxSize = 1052000);
  ~Collector() {
    delete[] a;
    delete[] m_rawMemory;
  };
  MInt size() { return (MInt)m_size; };

  MInt maxSize() { return (MInt)m_maxSize; };

  void append();
  MLong memoryUseage() { return (MLong)(m_usedMemory + (m_maxSize + 1) * sizeof(T)); };

  T* operator[](MInt index);
  T* a;
  MLong m_rawMemoryCounter;

  MInt setSize(MInt inputSize) {
    if(inputSize < m_size) m_size = inputSize;
    return (MInt)inputSize;
  }

  MInt resetSize(MInt inputSize) {
    m_size = inputSize;
    return (MInt)inputSize;
  }

  char* getRawPointer();

 private:
  MLong m_maxSize;
  MInt m_size;
  char* m_rawMemory = nullptr;
  MLong m_usedMemory;
  MInt m_distribution;
  MInt m_staticElementSize;
  void init(const MLong maxSize_) {
    m_size = 0;
    m_rawMemoryCounter = 0;
    m_maxSize = maxSize_;
  }

  /// \brief Allocates collector memory
  void allocMemory() {
    m_staticElementSize = T::staticElementSize();
    // An extra element is allocated in m_rawMemory to ensure that all element
    // variables are inside memory even when the all need to be realigned.
    m_usedMemory = m_staticElementSize * (m_maxSize + 16);
    a = new T[m_maxSize + 1];
    if(m_usedMemory > 0) {
      m_rawMemory = new char[m_usedMemory];
    } else {
      m_rawMemory = nullptr;
    }
  }
  /// \brief Allocates collector memory and initializes elements
  void allocMemoryAndInitElements() {
    allocMemory();
    for(MInt i = 0; i < m_maxSize + 1; i++) {
      a[i].allocateElements((void*)(m_rawMemory + m_rawMemoryCounter), (void*)m_rawMemory, i);
      m_rawMemoryCounter += m_staticElementSize;
    }
  }
};

// TODO labels:toenhance make collector constructors a single variadic template:
// template<class T, class.. Us>
// Collector<T>::Collector(MLong maxSize, Us... us) {
//   init(maxSize);
//   T::init(maxSize + 1, std::forward<Us>(us)...); // Why maxSize + 1 ? (dummy cell?)
//   allocMemoryAndInitElements();
// }

template <typename T>
Collector<T>::Collector(MLong inputMaxSize) {
  init(inputMaxSize);
  allocMemory();
  for(MInt i = 0; i < (inputMaxSize + 1); i++) {
    a[i].allocateElements((void*)(m_rawMemory + m_rawMemoryCounter));
    m_rawMemoryCounter += m_staticElementSize;
  }
}

/** Overloading for LB cells, which require the distributions as parameter.
 *
 * Also used for:
 *  - DG: m_elements
 *
 */
template <typename T>
Collector<T>::Collector(MLong inputMaxSize, MInt dimension, MInt dummy) {
  init(inputMaxSize);
  T::init(dimension, dummy, inputMaxSize);
  allocMemoryAndInitElements();
}

// Used for:
// - DG: m_surfaces
template <typename T>
Collector<T>::Collector(MLong inputMaxSize, MInt dimension, MInt distributions, MInt distributions1) {
  init(inputMaxSize);
  T::init(dimension, distributions, distributions1, inputMaxSize + 1);
  allocMemoryAndInitElements();
}

// detailed chemistry
template <typename T>
Collector<T>::Collector(MLong inputMaxSize, MInt dimension, MInt distributions, MInt distributions1, MInt dummy2) {
  init(inputMaxSize);
  T::init(dimension, distributions, distributions1, inputMaxSize + 1, dummy2);
  allocMemoryAndInitElements();
}


template <typename T>
Collector<T>::Collector(MLong inputMaxSize, MInt dimension, MInt distributions, MInt distributions1, MInt maxNoSurfaces,
                        MInt /*dummy1*/) {
  init(inputMaxSize);
  T::init(dimension, distributions, distributions1, inputMaxSize + 1, maxNoSurfaces);
  allocMemoryAndInitElements();
}

template <typename T>
Collector<T>::Collector(MLong inputMaxSize, MInt dimension, MFloat flameSpeed, MInt maxNoSets) {
  init(inputMaxSize);
  T::init(dimension, flameSpeed, inputMaxSize, maxNoSets);
  allocMemoryAndInitElements();
}

template <typename T>
void Collector<T>::append() {
  if(m_size < m_maxSize) {
    m_size++;
  } else {
    std::stringstream errorMessage;
    errorMessage << " Error in collector, maxSize reached ( " << m_maxSize << " elements ).";
    mTerm(1, AT_, errorMessage.str());
  }
}

template <typename T>
char* Collector<T>::getRawPointer() {
  return m_rawMemory;
}


template <typename T>
T* Collector<T>::operator[](MInt index) {
  return &a[index];
}


namespace maia {

/// \brief Helper functions useful for allocating collector memory
namespace collector_memory {

/// \brief Aligns pointer p such that a T stored at its adress is aligned
///
/// \param[in]   p          Pointer to align such that a T stored there is aligned
/// \param[out]  p_aligned  Aligned pointer >= p that can store an aligned T
///
/// \returns an aligned pointer > p such that a T stored there is aligned
template <class T, class U>
inline U* align(U* p) {
  const MLong N = sizeof(T);
  const MLong padding = reinterpret_cast<MLong>(static_cast<void*>(p)) % N;
  if(padding == 0) {
    return p;
  } else {
    return reinterpret_cast<U*>(reinterpret_cast<char*>(p) + padding);
  }
}

/// \brief Stores 1D variables in row-major order
///
/// I.e. the variables are stored in memory as follows:
/// cell_(0)_v(0) ... cell_(0)_v(Nrows) | ... | cell_(cellId)_v(0), cell_(cellId)_v(1), ..., cell_(cellId)_v(Nrows) |
/// cell_(cellId+1)_v(0), ..., | ...
/// ^^^^^^^^^^^^^ base points here!             ^^^^^^^^^^^^^^^^^^^ p is set to point here
///
/// \param[out]     p           Is set to point to the first cell variable.
/// \param[in,out]  base        Pointer to the start of the variable's memory region.
/// \param[in]      cellId      Id of the cell.
/// \param[in]      Nrows       Number of variables to store per cell.
/// \param[in]      maxNoCells  Maximum number of cells to be allocated (used for offset computation).
template <class T>
inline void rowMajor1D(T*& p, void*& base, const MInt cellId, const MInt Nrows, const MInt maxNoCells) {
  base = align<T>(base);
  p = reinterpret_cast<T*>(base) + static_cast<MLong>(1) * Nrows * cellId;
  base = reinterpret_cast<void*>(reinterpret_cast<T*>(base) + static_cast<MLong>(1) * Nrows * maxNoCells);
}

/// \brief Stores 2D variables in row-major order
///
/// In memory, an array of pointers is stored first. The variables itself are allocated behind it.
/// | ptr_0, ..., ptr_(maxNoCells*Nrows) | var_(0,0), ..., var_(maxNoCells*Nrows,maxNoCells*Ncols) |
///  ^^^^^^ base point here
///
/// The array of pointers provides 2D access and is stored as follows:
/// ... | ptr_to_cell_(cellId)_v(0,0), ptr_to_cell_(cellId)_v(1,0), ..., ptr_to_cell_(cellId)_v(Nrows,0) |
/// ptr_to_cell_(cellId+1)_v(0,0), ... | ...
///       ^^^^^^^^^^^^^^^^^^^^^^^^^^^^ p is set to point here!
///
/// The variables are stored as follows:
/// ...| cell_(cellId)_v(0,0), cell_(cellId)_v(0,1), ...,  cell_(cellId)_v(1,0), cell_(cellId)_v(1,1), ..., ....
/// cell_(cellId)_v(Nrows,Ncols) | cell_(cellId+1)_v(0,0), ..., | ...
///     ^^^^^^^^^^^^^^^^^^^^^ ptr_to_cell_(cellId)_v(0,0)  ^^^^^^^^^^^^^ ptr_to_cell_(cellId)_v(1,0) points here!
///
/// \param[out]     p           Is set to point to the first pointer of the cell variables.
/// \param[in,out]  base        Pointer to the start of the variable's memory region.
/// \param[in]      cellId      Id of the cell.
/// \param[in]      Nrows       Number of variables to store per cell.
/// \param[in]      maxNoCells  Maximum number of cells to be allocated (used for offset computation).
template <class T>
ATTRIBUTES1(ATTRIBUTE_NO_AUTOVEC)
inline void rowMajor2D(T**& p, void*& base, const MInt cellId, const MInt Nrows, const MInt Ncols,
                       const MInt maxNoCells) {
  // pointers are stored first in row major order:
  rowMajor1D(p, base, cellId, Nrows, maxNoCells);
  // the data itself is stored after the pointers:
  base = align<T>(base);
  for(MInt i = 0; i < Nrows; ++i) {
    p[i] =
        reinterpret_cast<T*>(base) + static_cast<MLong>(1) * Nrows * Ncols * cellId + static_cast<MLong>(1) * i * Ncols;
  }
  base = reinterpret_cast<T*>(base) + static_cast<MLong>(1) * Nrows * Ncols * maxNoCells;
}

/// \brief Copies 1D cell elements
///
/// \param[in]  to     Elements will be copied here.
/// \param[in]  from   Elements to be copied.
/// \param[in]  Nrows  Number of elements to copy
template <class T>
inline void copyElements1D(T* to, T* from, const MInt Nrows) {
  for(MInt i = 0; i < Nrows; ++i) {
    to[i] = from[i];
  }
}

/// \brief Copies 2D cell elements
///
/// Nrows*Ncols elements will be copied.
///
/// \param[in]  to     Elements will be copied here.
/// \param[in]  from   Elements to be copied.
/// \param[in]  Nrows  Number of rows
/// \param[in]  Ncols  Number of columns
template <class T>
inline void copyElements2D(T* to, T* from, const MInt Nrows, const MInt Ncols) {
  for(MInt i = 0; i < Nrows; ++i) {
    for(MInt j = 0; j < Ncols; ++j) {
      to[i][j] = from[i][j];
    }
  }
}

/// This namespaces includes helpers for managing collector memory "cell-wise"
/// This is _VERY_ wrong and should not be done.
/// The memory is _NOT_ aligned!
namespace unaligned_cell_wise {
/// \brief Store 1D variables in row-majow order relative to the start of
/// each cell's memory solver
///
/// NOTE: DO NOT USE THIS, IS VERY BAD FOR PERFORMANCE, ALIGNMENT, REQUIRED MEMORY....
///
/// Some cells use the collector as follows:
///
/// | ....... collector memory ......... |
///
/// |.....| stuff of the cell i | ..... |
///        ^^^^^^^^^^^^^^^^^^^^^ (expanded below)
/// |.....| array_of_floats, array_of_ints, array_of_floats,... of the cell i | ..... |
///
/// This means that when you loop over some cell pointer, you don't go linearly through memory
/// but you have random access (bad for performance)
///
/// Furthermore: each of the arrays inside each cell needs to be aligned, so the padding
/// you need is not a function of the #of_member_variables but a function of both the
/// #of_cells _AND_ the #of_member_variables -> If a single array is not aligned you
/// need #maxNoCells*padding ! (bad for memory)
///
template <class T>
inline void rowMajor1D(T*& p, void*& base, const MInt Nrows) {
  // rowMajor1D(p,base,0,Nrows,1); // DO NOT ALIGN -> requires too much memory
  p = reinterpret_cast<T*>(base);
  base = reinterpret_cast<void*>(reinterpret_cast<T*>(base) + Nrows);
}

template <class T>
ATTRIBUTES1(ATTRIBUTE_NO_AUTOVEC)
inline void rowMajor2D(T**& p, void*& base, const MInt Nrows, const MInt Ncols) {
  // rowMajor2D(p,base,0,Nrows,Ncols,1); // DO NOT ALIGN -> requires too much memory
  // pointers are stored first in row major order:
  rowMajor1D(p, base, Nrows);
  // the data itself is stored after the pointers:
  for(MInt i = 0; i < Nrows; ++i) {
    p[i] = reinterpret_cast<T*>(base) + static_cast<MLong>(1) * i * Ncols;
  }
  base = reinterpret_cast<T*>(base) + static_cast<MLong>(1) * Nrows * Ncols;
}

} // namespace unaligned_cell_wise

} // namespace collector_memory

} // namespace maia
#endif
