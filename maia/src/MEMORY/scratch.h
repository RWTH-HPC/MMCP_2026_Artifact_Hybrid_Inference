// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef SCRATCH_H
#define SCRATCH_H

#ifndef PVPLUGIN // =============== DO MAIA SCRATCH  ========

#include <cstdint>
#include <iomanip>
#include <list>
#include <new>
#include <ostream>
#include <sstream>
#include <utility>
#include "INCLUDE/maiatypes.h"
#include "IO/infoout.h"
#include "UTIL/functions.h"
#include "config.h"
#include "globalvariables.h"

//#if defined(linux)
//#include <sys/sysinfo.h>
// const uintptr_t ALIGNMENT_BOUNDARY = sysconf(_SC_LEVEL1_DCACHE_LINESIZE); // make alignment boundary equal to
// L1-cache linesize #else const uintptr_t ALIGNMENT_BOUNDARY = 64; // see below #endif


/// \brief type traits for dealing with signed vs unsigned indices
///
/// Note: we should not do this...
/// \todo labels:toenhance use <type_traits> instead of this hack
namespace maia {
struct maia_signed {};
struct maia_unsigned {};
template <class T>
struct is_unsigned {
  static const MBool value = false;
  using type = maia_signed;
};
template <>
struct is_unsigned<MUint> {
  static const MBool value = true;
  using type = maia_unsigned;
};
// used in cartesian grid constructor, comment out to see the errors:
template <>
struct is_unsigned<MUlong> {
  static const MBool value = true;
  using type = maia_unsigned;
};
} // namespace maia

template <class T>
class ScratchSpace;
class ScratchSpaceBase;

using MStringScratchSpace = ScratchSpace<MString>;
using MCharScratchSpace = ScratchSpace<MChar>;
using MIntScratchSpace = ScratchSpace<MInt>;
using MIntScratchSpace = ScratchSpace<MInt>;
using MLongScratchSpace = ScratchSpace<MLong>;
using MFloatScratchSpace = ScratchSpace<MFloat>;
using MBoolScratchSpace = ScratchSpace<MBool>;
using MUshortScratchSpace = ScratchSpace<MUshort>;
using MUcharScratchSpace = ScratchSpace<MUchar>;

using MIntPointerScratchSpace = ScratchSpace<MInt*>;
using MIntPointerScratchSpace = ScratchSpace<MInt*>;
using MLongPointerScratchSpace = ScratchSpace<MLong*>;
using MFloatPointerScratchSpace = ScratchSpace<MFloat*>;
using MBoolPointerScratchSpace = ScratchSpace<MBool*>;

using MFloatPointerPointerScratchSpace = ScratchSpace<MFloat**>;

using MPI_GroupScratchSpace = ScratchSpace<MPI_Group>;

using ScratchList = std::list<ScratchSpaceBase*>;

/** \brief This class holds the complete scratch space
 *
 * \authors Andreas Lintermann, Stephan Schlimpert, Christoph Siewert
 * \date 10.05.2011, 10.06.2011, 21.10.2011
 *
 * This class defines the overall size of the scratch space and tracks
 * the sizes and allocation of newly requested Scratch space arrays.
 *
 **/
class Scratch {
  friend class ScratchSpaceBase;

  static const uintptr_t ALIGNMENT_BOUNDARY = MAIA_SCRATCH_ALIGNMENT_BOUNDARY; // see above

 public:
  /** \brief Constructor
   *
   * \authors Andreas Lintermann, Stephan Schlimpert, Christoph Siewert
   * \date 10.05.2011, 10.06.2011, 21.10.2011
   * \note Base memory now points to an aligned memory address, Lennart Schneiders, 12.12.2012
   *
   * Allocates space of "size * sizeof(MFloat) * Cells" bytes, that is "size" MFloat for every cell.
   *
   * \param[in] size number of MFloat for every cell.
   * \param[in] Cells number of cells.
   *
   **/
  Scratch(MFloat size, MInt Cells);

  /** \brief Destructor
   *
   * \author Andreas Lintermann, Christoph Siewert
   * \date 10.05.2011, 21.10.2011
   *
   **/
  ~Scratch() {
    m_number_of_cells = 0;
    m_number_of_elements = 0;
    m_object_id = 0;
    m_usedmemsize = 0;
    m_nextfree = nullptr;
    m_scratchSpaces.clear();
    m_maxused = nullptr;
    m_report = "n/a";
    delete[] m_totalScratch;
    m_totalScratch = nullptr;
  };

  /** \brief Returns a pointer to the end of the scratch.
   *
   * \authors Andreas Lintermann, Stephan Schlimpert, Christoph Siewert
   * \date 10.05.2011, 10.06.2011, 21.01.2011
   *
   * \return a pointer to the end of the scratch
   **/
  static char* getEndPointer() { return (char*)(size_t)(m_totalScratch + m_number_of_elements - 1); }

  /** \brief Returns the amount of available memory in scratch.
   *
   * \author Andreas Lintermann, Christoph Siewert
   * \date 10.05.2011, 21.10.2011
   *
   * \return the amount of available memory
   **/
  // static size_t getAvailableMemory() {return (getTotalMemory() - m_usedmemsize);}//this leads to errors if the
  // m_nextfree pointer is corrupted for whatever reason, Lennart
  static size_t getAvailableMemory() { return (&m_totalScratch[m_number_of_elements - 1] - m_nextfree); }

  /** \brief Returns the amount of total available memory in scratch.
   *
   * \author Andreas Lintermann. Christoph Siewert
   * \date 10.05.2011, 21.10.2011
   *
   * \return the amount of total available memory
   **/
  static size_t getTotalMemory() { return (getEndPointer() + 1 - m_totalScratch); }

  static MString printSelfScratch();
  static MString printSelf();
  static MString printSelfReport();

  static size_t m_number_of_cells;
  static size_t m_number_of_elements;
  static MInt m_object_id;
  static char* m_totalScratch;
  static size_t m_usedmemsize;
  static char* m_nextfree;
  static ScratchList m_scratchSpaces;
  static char* m_maxused;
  static MString m_report;
};


/** \brief This class is a base class for all ScratchSpaces.
 *
 * \author Andreas Lintermann, Christoph Siewert
 * \date 10.05.2011, 21.10.2011
 *
 * This class is the base of all derived scratch space elements.
 **/
class ScratchSpaceBase {
 public:
  static const uintptr_t ALIGNMENT_BOUNDARY = MAIA_SCRATCH_ALIGNMENT_BOUNDARY; // see above

  const size_t m_memsize;
  const size_t m_memsizePadded;
  const MString m_calling_function;
  const MString m_variable_name;
  MInt m_object_id{};
  MBool m_nonterminal;
  MBool m_destroy;

  /** \brief Constructor
   *
   * \author Andreas Lintermann, Christoph Siewert
   * \date 10.05.2011, 21.10.2011
   * \note Each ScratchSpace now points to an aligned memory address, Lennart Schneiders, 12.12.2012
   *
   * Sets the scratch space elements properties.
   *
   * \param[num] num number of elements in array
   * \param[size] size the size of the array in bytes
   * \param[name] name the name of the calling function
   * \param[varname] varname the name of the array variable
   *
   **/
  ScratchSpaceBase(size_t num, size_t size, MString name, MString varname)
    : m_memsize(num * size),
      m_memsizePadded((num * size) + ((ALIGNMENT_BOUNDARY - ((num * size) % ALIGNMENT_BOUNDARY)) % ALIGNMENT_BOUNDARY)),
      m_calling_function(std::move(name)),
      m_variable_name(std::move(varname)),
      m_nonterminal(false),
      m_destroy(false) {
    ASSERT((((uintptr_t)Scratch::m_nextfree) % ALIGNMENT_BOUNDARY == 0), "Scratch memory is not aligned");
  }
  virtual MString printSelfReport() const = 0;
  virtual MString printSelf() const = 0;
};

/** \brief This class is a ScratchSpace.
 *
 * \author Andreas Lintermann, Christoph Siewert
 * \date 10.05.2011, 21.10.2011
 *
 * This class is a derivative of ScratchSpaceBase and is responsible
 * for the management of the scratch space to be allocated. Scratch objects
 * partially model the container concept and should be usable with some,
 * but not all, STL algorithms.
 **/
template <class T>
class ScratchSpace : public ScratchSpaceBase {
 public:
  using value_type = T;
  using reference = value_type&;
  using const_reference = const value_type&;
  using pointer = value_type*;
  using const_pointer = const value_type*;
  using iterator = pointer;
  using const_iterator = const_pointer;
  using difference_type = long;
  using size_type = std::size_t;

  ScratchSpace(MInt num, const MString& name, const MString& varname); //!< 1D Scratch constructor
  ScratchSpace(MInt num0, MInt num1, const MString& name,
               const MString& varname); //!< 2D Scratch constructor
  ScratchSpace(MInt num0, MInt num1, MInt num2, const MString& name,
               const MString& varname); //!< 3D Scratch constructor
  ~ScratchSpace();

  template <class Integral>
  inline reference operator[](const Integral /*i*/); //!< 1D Scratch access
  template <class Integral>
  inline reference operator()(const Integral /*i*/); //!< 1D Scratch access
  template <class Integral>
  inline reference operator()(const Integral /*i*/, const Integral /*j*/); //!< 2D Scratch access
  template <class Integral>
  inline reference operator()(const Integral /*i*/, const Integral /*j*/, const Integral /*k*/); //!< 3D Scratch access
  template <class Integral>
  inline const_reference operator[](const Integral i) const; //!< 1D Scratch const access
  template <class Integral>
  inline const_reference operator()(const Integral i) const; //!< 1D Scratch const access
  template <class Integral>
  inline const_reference operator()(const Integral /*i*/, const Integral /*j*/) const; //!< 2D Scratch const access
  template <class Integral>
  inline const_reference operator()(const Integral /*i*/, const Integral /*j*/,
                                    const Integral /*k*/) const; //!< 3D Scratch const access
  ScratchSpace<T>&
  operator=(ScratchSpace<T>& /*S*/); //!< copy entries of S to this, given there is enough space available
  template <class>
  friend std::ostream& operator<<(std::ostream& os, const ScratchSpace& s);

  inline iterator begin() {
    checkForEmptyScratch();
    return p;
  }
  inline iterator cbegin() const {
    checkForEmptyScratch();
    return p;
  }
  inline iterator end() {
    checkForEmptyScratch();
    return last;
  }
  inline iterator cend() const {
    checkForEmptyScratch();
    return last;
  }
  pointer data() {
    checkForEmptyScratch();
    return begin();
  }
  const_pointer data() const {
    checkForEmptyScratch();
    return begin();
  }

  inline MInt size0() const { return m_size0; }
  inline MInt size1() const { return m_size1; }
  inline MInt size2() const { return m_size2; }
  inline size_t getMemsize() const { return (last - p) * sizeof(T); }
  inline size_type size() const {
    ASSERT(last >= p, "Scratch cannot have negative size!");
    return last - p;
  }
  inline MBool empty() const {
    ASSERT(last >= p, "Scratch cannot have negative size!");
    return static_cast<bool>(last - p == 0);
  }

  void fill(T val) { std::fill(p, last, val); } //!< fill the scratch with a given value
  MString printSelf() const override;
  MString printSelfReport() const override;

  pointer p;                          ///< Deprecated: use [] instead!
  T* getPointer() const { return p; } ///< Deprecated: use begin() instead!

  // Disable copying and heap allocation
  ScratchSpace() = delete;
  ScratchSpace(const ScratchSpace<T>&) = delete;

  void* operator new(std::size_t) = delete;
  void* operator new(std::size_t, void* p) = delete;
  void* operator new[](std::size_t) = delete;
  void* operator new[](std::size_t, void* p) = delete;
  void operator delete(void*) = delete;
  void operator delete(void* p, void*) = delete;
  void operator delete[](void*) = delete;
  void operator delete[](void* p, void*) = delete;

 private:
  // TODO labels:toenhance when ".p" is removed: const_iterator first, last;
  pointer last;
  // TODO labels:toenhance this should be std::size_t
  const MInt m_size0;
  const MInt m_size1;
  const MInt m_size2;

  /// \brief Assert the bounds of 1D acces to an array
  ///
  /// Note: the condition is tested by checkBounds,
  /// which is overloaded for signed and unsigned indices
  /// using the maia::is_unsigned trait.
  ///
  /// \todo labels:toenhance replace maia::is_unsigned with <type_traits>
  ///
  /// @{
  template <class Integral>
  inline void testBounds(const Integral& i) const;
  template <class Integral>
  inline MBool checkBounds(const Integral& i) const {
    using type = typename maia::is_unsigned<Integral>::type;
    return checkBounds(i, type());
  }
  template <class Integral>
  inline MBool checkBounds(const Integral& i, maia::maia_unsigned /*unused*/) const {
    return static_cast<size_type>(i) < size();
  }
  template <class Integral>
  inline MBool checkBounds(const Integral& i, maia::maia_signed /*unused*/) const {
    return i >= 0 && static_cast<size_type>(i) < size();
  }
  /// @}
  void init();
  inline void checkForEmptyScratch() const;
};


/** \brief Constructor
 *
 * \author Andreas Lintermann, Christoph Siewert
 * \date 10.05.2011, 21.10.2011
 *
 * Allocates a scratch space element of a certain type.
 * If not enough memory is available provided by the Scratch class, it terminates the program run.
 * Otherwise the properties of this management object are set and a pointer to this object is stored in the management
 *list of the Scratch class. If a maximum of memory of the Scratch is used, the state stored to be outputed at
 *the end of the program run to enable memory analysis.
 *
 * \param[num] num number of elements in array
 * \param[name] name the name of the calling function
 * \param[varname] varname the name of the array variable
 *
 * Note: If MAIA_EXTRA_DEBUG is defined writes a SCRATCH_DOUBLE_ALLOCATION WARNING
 * at run-time if allocation of two variables with the same name is performed.
 **/
template <class T>
ScratchSpace<T>::ScratchSpace(MInt num, const MString& name, const MString& varname)
  : ScratchSpaceBase(mMax(1, num), sizeof(T), name, varname),
    p(reinterpret_cast<T*>(Scratch::m_nextfree)),
    last(p + num),
    m_size0(num),
    m_size1(1),
    m_size2(1) {
  this->init();
}

/** \brief Constructor for 2D scratch
 *
 * \author Lennart
 * \date 26.02.2013
 *
 * This is the 2D-extension of the 1D constructor above. Unfortunately, the order of constructor arguments
 * prohibits the use of a constructor default value for the second dimension, so most of the code here is redundant.
 *
 **/
template <class T>
ScratchSpace<T>::ScratchSpace(MInt num0, MInt num1, const MString& name, const MString& varname)
  : ScratchSpaceBase(mMax(1, num0 * num1), sizeof(T), name, varname),
    p(reinterpret_cast<T*>(Scratch::m_nextfree)),
    last(p + (num0 * num1)),
    m_size0(num0),
    m_size1(num1),
    m_size2(1) {
  this->init();
}


/** \brief Constructor for 3D scratch
 *
 * \author Lennart
 * \date 26.02.2013
 *
 * This is the 3D-extension of the 1D constructor above. Unfortunately, the order of constructor arguments
 * prohibits the use of constructor default values for the second and third dimension, so most of the code here is
 *redundant.
 *
 **/
template <class T>
ScratchSpace<T>::ScratchSpace(MInt num0, MInt num1, MInt num2, const MString& name, const MString& varname)
  : ScratchSpaceBase(mMax(1, num0 * num1 * num2), sizeof(T), name, varname),
    p(reinterpret_cast<T*>(Scratch::m_nextfree)),
    last(p + (num0 * num1 * num2)),
    m_size0(num0),
    m_size1(num1),
    m_size2(num2) {
  this->init();
}


/** \brief extended constructor functionality
 *
 * \author Andreas Lintermann, Christoph Siewert
 * \date 10.05.2011, 21.10.2011
 * \brief Content was moved from the constructor body to this function since it is the same for all above constructors
 *(Lennart) \note added 'placement new' operator, as previously the use of non-primitive types resulted in undefined
 *behavior, Lennart 2013
 *
 **/
template <class T>
void ScratchSpace<T>::init() {
  Scratch::m_object_id++;

  m_object_id = Scratch::m_object_id;
  Scratch::m_scratchSpaces.push_back(this);

  /// Check if allocation would exceed the available scratch space:
  if(Scratch::getAvailableMemory() < m_memsizePadded) {
    m_destroy = true;
    std::cerr << Scratch::printSelf();
    m_log << Scratch::printSelf();
    std::stringstream errorMessage, suggestionMessage;
    errorMessage << "Scratch is not big enough - exceeded size in " << m_calling_function << " with variable "
                 << m_variable_name << std::endl;
    errorMessage << "           " << m_memsize / (1024.0 * 1024.0) << " MB required, but there are only "
                 << Scratch::getAvailableMemory() / (1024.0 * 1024.0) << " MB available" << std::endl;
    errorMessage << "           Suggested to raise the property scratchSize from "
                 << MFloat(Scratch::m_number_of_elements) / (sizeof(MFloat) * Scratch::m_number_of_cells) << " to "
                 << MFloat(Scratch::m_number_of_elements + m_memsize - Scratch::getAvailableMemory())
                        / (sizeof(MFloat) * Scratch::m_number_of_cells)
                 << std::endl;
    m_log << errorMessage.str();
    mTerm(1, AT_, errorMessage.str());
  }
  Scratch::m_usedmemsize += m_memsizePadded;
  Scratch::m_nextfree += m_memsizePadded;

  for(MInt i = 0; i < (signed)this->size(); i++) {
    ASSERT(((char*)&p[i] >= &Scratch::m_totalScratch[0])
               && ((char*)&p[i] <= &Scratch::m_totalScratch[Scratch::m_number_of_elements - 1]),
           "ScratchSpace '" + m_variable_name + "' out of bounds during init.");
    new(&p[i]) T(); // 'placement new' operator calling each element's default constructor (no further memory allocated)
  }

  for(MInt i = 0; i < (signed)this->size(); i++) {
    new(&p[i]) T(); // 'placement new' operator calling each element's default constructor
  }

/// Check for double allocation:
#ifdef MAIA_EXTRA_DEBUG
  for(ScratchList::iterator it = Scratch::m_scratchSpaces.begin(); it != Scratch::m_scratchSpaces.end(); ++it) {
    if((*it) == this) continue;
    if((*it)->m_variable_name == this->m_variable_name) {
      m_log << "WARNING: SCRATCH_DOUBLE_ALLOCATION: double usage of variable " << this->m_variable_name
            << " allocated in " << (*it)->m_calling_function << " and " << this->m_calling_function << std::endl;
    }
  }
#endif

  if(Scratch::m_nextfree > Scratch::m_maxused) {
    Scratch::m_maxused = Scratch::m_nextfree;
    Scratch::m_report = Scratch::printSelfScratch();
    ScratchList::iterator iter;
    for(iter = Scratch::m_scratchSpaces.begin(); iter != Scratch::m_scratchSpaces.end(); iter++) {
      Scratch::m_report += (*iter)->printSelfReport();
      Scratch::m_report += "\n";
    }
    Scratch::m_report += "\n\n";
  }
}


/** \brief Destructor
 *
 * \author Andreas Lintermann, Christoph Siewert
 * \date 10.05.2011, 21.10.2011
 *
 * Removes this scratch space elemet from the list of all scratch elements.
 **/
template <class T>
ScratchSpace<T>::~ScratchSpace() {
  if(!m_destroy) {
    auto it = Scratch::m_scratchSpaces.end();
    it--;

    //! Only delete if we are at the last element (guaranteed by C++)
    if((*it)->m_object_id == m_object_id) {
      Scratch::m_nextfree -= m_memsizePadded;
      Scratch::m_usedmemsize -= m_memsizePadded;
      Scratch::m_scratchSpaces.erase(it);
    }

    if(m_nonterminal) {
      Scratch::m_nextfree -= m_memsizePadded;
    }
  }
}

/** \brief Returns a string summing up this scratch space element information.
 *
 * \authors Andreas Lintermann, Stephan Schlimpert, Christoph Siewert
 * \date 10.05.2011, 10.06.2011, 21.10.2011
 *
 * \return a string summing up this scratch space element information
 **/
template <class T>
MString ScratchSpace<T>::printSelf() const {
  std::stringstream id, mem, memsize;
  mem << p;

  id << m_object_id;
  memsize << m_memsize;

  MString message = "ScratchSpace ID: ";
  message += id.str();
  message += "\n-----------------------------\n";
  message += "Calling function:\t";
  message += m_calling_function;
  message += "\nVariable name:\t\t";
  message += m_variable_name;
  message += "\nMemory size:\t\t";
  message += memsize.str();
  message += "\nMemory pointer:\t\t";
  message += mem.str();
  message += "\nSize:\t\t\t";
  message += std::to_string(size());
  message += "\n";

  return message;
}

/** \brief Returns a shortened string summing up this scratch space element information.
 *
 * \author Andreas Lintermann, Christoph Siewert
 * \date 10.05.2011, 21.10.2011
 *
 * This function is used in the report process for tracking the occasion of maximal
 * memory usage during program execution.
 *
 * \return a shortened string summing up this scratch space element information
 **/
template <class T>
MString ScratchSpace<T>::printSelfReport() const {
  std::stringstream id, memsize, funct, var;
  memsize << m_memsize;
  id << m_object_id;

  MString message = "ID: ";
  message += id.str();
  message += "    SIZE: ";
  message += memsize.str();
  message += "    FUNCT: ";
  message += m_calling_function;
  message += "    VAR: ";
  message += m_variable_name;
  return message;
}

template <class T>
template <class Integral>
inline void ScratchSpace<T>::testBounds(const Integral& i) const {
  ASSERT(checkBounds(i), "SCRATCH_OUT_OF_BOUNDS ERROR: accessing "
                             << m_variable_name << " allocated from " << m_calling_function << " with index " << i
                             << " out of range [0," << size() << ")." << std::endl
                             << "Note: if the array range = [0,0) the array is EMPTY. "
                             << "You cannot dereference elements of an empty array!");
}

/// \brief Access the i-th element of the scratch object.
///
/// \return Refernence to the i-th scratch space element
///
/// \param[i] Index of Integral type.
/// \param[Integral] Type that models the Integral concept.
///
/// Note: if NDEBUG is not defined, the array bounds are checked:
/// i >= 0 and i < size.
/// If bound checking fails the program will terminate with a
/// SCRATCH_OUT_OF_BOUNDS ERROR.
///
template <class T>
template <class Integral>
inline typename ScratchSpace<T>::reference ScratchSpace<T>::operator()(const Integral i) {
#if defined(MAIA_ASSERTS)
  testBounds(i);
#endif
  return p[i];
}
template <class T>
template <class Integral>
inline typename ScratchSpace<T>::const_reference ScratchSpace<T>::operator()(const Integral i) const {
#if defined(MAIA_ASSERTS)
  testBounds(i);
#endif
  return p[i];
}

/// For 1D access the [] operator is also provided, the () operator should be preferred.
template <class T>
template <class Integral>
inline typename ScratchSpace<T>::reference ScratchSpace<T>::operator[](const Integral i) {
  return operator()(i);
}
template <class T>
template <class Integral>
inline typename ScratchSpace<T>::const_reference ScratchSpace<T>::operator[](const Integral i) const {
  return operator()(i);
}

/// \brief 2D and 3D Scratch access
///
/// Note: this would probably result in a compilation error
/// if m_sizeX
template <class T>
template <class Integral>
inline typename ScratchSpace<T>::reference ScratchSpace<T>::operator()(const Integral i, const Integral j) {
  ASSERT(i >= 0 && i < m_size0 && j >= 0 && j < m_size1,
         "SCRATCH_OUT_OF_BOUNDS ERROR: accessing " << m_variable_name << " allocated from " << m_calling_function
                                                   << " with indices " << i << "," << j << ". Dimensions are "
                                                   << m_size0 << "," << m_size1 << ".");
  return p[m_size1 * i + j];
}

template <class T>
template <class Integral>
inline typename ScratchSpace<T>::reference ScratchSpace<T>::operator()(const Integral i, const Integral j,
                                                                       const Integral k) {
  ASSERT(i >= 0 && i < m_size0 && j >= 0 && j < m_size1 && k >= 0 && k < m_size2,
         "SCRATCH_OUT_OF_BOUNDS ERROR: accessing "
             << m_variable_name << " allocated from " << m_calling_function << " with indices " << i << "," << j << ","
             << k << ". Dimensions are " << m_size0 << "," << m_size1 << "," << m_size2 << ".");
  return p[m_size1 * m_size2 * i + m_size2 * j + k];
}

template <class T>
template <class Integral>
inline typename ScratchSpace<T>::const_reference ScratchSpace<T>::operator()(const Integral i, const Integral j) const {
  ASSERT(i >= 0 && i < m_size0 && j >= 0 && j < m_size1,
         "SCRATCH_OUT_OF_BOUNDS ERROR: accessing " << m_variable_name << " allocated from " << m_calling_function
                                                   << " with indices " << i << "," << j << ". Dimensions are "
                                                   << m_size0 << "," << m_size1 << ".");
  return p[m_size1 * i + j];
}

template <class T>
template <class Integral>
inline typename ScratchSpace<T>::const_reference ScratchSpace<T>::operator()(const Integral i, const Integral j,
                                                                             const Integral k) const {
  ASSERT(i >= 0 && i < m_size0 && j >= 0 && j < m_size1 && k >= 0 && k < m_size2,
         "SCRATCH_OUT_OF_BOUNDS ERROR: accessing "
             << m_variable_name << " allocated from " << m_calling_function << " with indices " << i << "," << j << ","
             << k << ". Dimensions are " << m_size0 << "," << m_size1 << "," << m_size2 << ".");
  return p[m_size1 * m_size2 * i + m_size2 * j + k];
}

template <class T>
ScratchSpace<T>& ScratchSpace<T>::operator=(ScratchSpace<T>& S) {
  for(MInt i = 0; i < mMin(m_size0, S.m_size0); i++) {
    for(MInt j = 0; j < mMin(m_size1, S.m_size1); j++) {
      for(MInt k = 0; k < mMin(m_size2, S.m_size2); k++) {
        p[m_size1 * m_size2 * i + m_size2 * j + k] = S.p[m_size1 * m_size2 * i + m_size2 * j + k];
      }
    }
  }
  return *this;
}


/// \brief Print contents of scratch space object.
///
/// \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
/// \date 2015-03-13
///
/// \tparam T Element type.
/// \param[in] os Stream object to write to.
/// \param[in] s Scratch space object.
///
/// \return Stream object to allow chained writing
///
/// Note: if T is not convertible to string, the contents are not printed.
template <class T>
std::ostream& operator<<(std::ostream& os, const ScratchSpace<T>& s) {
  os << s.printSelf();
  os << "Content:\n";
  for(typename ScratchSpace<T>::size_type i = 0; i < s.size(); i++) {
    os << std::setw(7) << i << ": " << s[i] << "\n";
  }
  return os;
}


/// \brief Checks if scratch space is empty before returning a pointer to it.
///
/// Note: the check is only performed if NDEBUG is not defined,
/// and MAIA_EXTRA_DEBUG is defined.
///
template <class T>
inline void ScratchSpace<T>::checkForEmptyScratch() const {
#ifdef MAIA_EXTRA_DEBUG
  if(empty()) {
    m_log << "SCRATCH WARNING: passing a pointer to the empty scratch space " << m_variable_name << " allocated from "
          << m_calling_function << "! Make sure this pointer will NEVER be dereferenced!" << std::endl;
  }
#endif
}

#else // =============== DO PARAVIEW PLUGIN SCRATCH  ========

#include <limits>
#include <memory>
#include <vector>
#include "src/INCLUDE/maiatypes.h"
#include "src/UTIL/functions.h"

using std::to_string;

template <typename T>
class ScratchSpace {
  using Storage = std::vector<T>;

 public:
  using value_type = T;
  using reference = typename Storage::reference;
  using const_reference = typename Storage::const_reference;
  using pointer = typename Storage::pointer;
  using const_pointer = typename Storage::const_pointer;
  using iterator = typename Storage::iterator;
  using const_iterator = typename Storage::const_iterator;
  using size_type = MLong;

  ScratchSpace(const size_type size_, const MString&, const MString&)
    : m_size0(size_), m_size1(1), m_data(m_size0 * m_size1), p{*this} {}
  ScratchSpace(const size_type size0, const size_type size1, const MString&, const MString&)
    : m_size0(size0), m_size1(size1), m_data(m_size0 * m_size1), p{*this} {}

  size_type size() const { return m_data.size(); }
  void fill(const T value) { std::fill(begin(), end(), value); }

  reference operator[](const size_type pos) { return m_data.at(pos); }
  reference operator()(const size_type pos) { return m_data.at(pos); }
  const_reference operator[](const size_type pos) const { return m_data.at(pos); }
  const_reference operator()(const size_type pos) const { return m_data.at(pos); }


  reference operator()(const size_type i, const size_type j) {
    if(i >= m_size0) mTerm(1, FUN_, "index i = " + to_string(i) + " out-of-bounds [0, " + to_string(m_size0) + ")");
    if(j >= m_size1) mTerm(1, FUN_, "index j = " + to_string(j) + " out-of-bounds [0, " + to_string(m_size1) + ")");
    return m_data.at(m_size1 * i + j);
  }
  const_reference operator()(const size_type i, const size_type j) const {
    if(i >= m_size0) mTerm(1, FUN_, "index i = " + to_string(i) + " out-of-bounds [0, " + to_string(m_size0) + ")");
    if(j >= m_size1) mTerm(1, FUN_, "index j = " + to_string(j) + " out-of-bounds [0, " + to_string(m_size1) + ")");
    return m_data.at(m_size1 * i + j);
  }

  T* getPointer() { return m_data.data(); }
  pointer data() { return m_data.data(); }
  const_pointer data() const { return m_data.data(); }

  iterator begin() { return m_data.begin(); }
  iterator end() { return m_data.end(); }
  iterator cbegin() const { return m_data.cbegin(); }
  iterator cend() const { return m_data.cend(); }

  struct p_ {
    ScratchSpace<T>& s;
    reference operator[](size_type i) { return s[i]; }
    const_reference operator[](size_type i) const { return s[i]; }
  };

  inline MInt size0() const { return m_size0; }
  inline MInt size1() const { return m_size1; }
  inline MInt size2() const { return m_size2; }

 private:
  size_type m_size0;
  size_type m_size1;
  size_type m_size2;
  Storage m_data;

 public:
  p_ p;
};

using MStringScratchSpace = ScratchSpace<MString>;
using MCharScratchSpace = ScratchSpace<MChar>;
using MIntScratchSpace = ScratchSpace<MInt>;
using MIntScratchSpace = ScratchSpace<MInt>;
using MLongScratchSpace = ScratchSpace<MLong>;
using MFloatScratchSpace = ScratchSpace<MFloat>;
using MBoolScratchSpace = ScratchSpace<MBool>;
using MUshortScratchSpace = ScratchSpace<MUshort>;
using MUcharScratchSpace = ScratchSpace<MUchar>;

class Scratch {
 public:
  static MFloat getAvailableMemory() { return 100; }
};
#endif // PV/MAIA Scratch

#endif // SCRATCH_H
