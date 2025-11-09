// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef CONTAINER_H_
#define CONTAINER_H_

#include <algorithm>
#include <type_traits>
#include <vector>
#include "INCLUDE/maiamacro.h"
#include "INCLUDE/maiatypes.h"
#include "UTIL/functions.h"

// The macro 'CONTAINER_SANITY_CHECKS' enables (potentially expensive) sanity checks for many
// operations. It is enabled for build type "debug"
#ifndef NDEBUG
#define CONTAINER_SANITY_CHECKS
#endif

#define MAIA_CONTAINER_ENSURE(condition, message, at)                                                                  \
  if(!(condition)) {                                                                                                   \
    TERMM(1, std::string("\n\n") + (message) + "\n\n AT: " + at);                                                      \
  }                                                                                                                    \
  do {                                                                                                                 \
  } while(false)

#define MAIA_CONTAINER_ENSURE_VALID_ID(id)                                                                             \
  do {                                                                                                                 \
    MAIA_CONTAINER_ENSURE(this->isValidId((id)),                                                                       \
                          "id = " + std::to_string((id)) + " out-of-bounds [0, " + std::to_string(this->size())        \
                              + ") and is not the dummy cell at \"capacity() - 1\" = "                                 \
                              + std::to_string(this->capacity() - 1),                                                  \
                          AT_);                                                                                        \
  } while(false)


// Sanity-checking macros for normal methods
#ifdef CONTAINER_SANITY_CHECKS
#define ENSURE_VALID_ID(id)                                                                                            \
  do {                                                                                                                 \
    MAIA_CONTAINER_ENSURE_VALID_ID(id);                                                                                \
  } while(false)
#define ENSURE_CONDITION(condition, message)                                                                           \
  do {                                                                                                                 \
    MAIA_CONTAINER_ENSURE(condition, message, AT_);                                                                    \
  } while(false)
#else
#define ENSURE_VALID_ID(id)                                                                                            \
  do {                                                                                                                 \
  } while(false)
#define ENSURE_CONDITION(condition, message)                                                                           \
  do {                                                                                                                 \
  } while(false)
#endif


namespace maia {
namespace container {

/// Base class for collector-like containers.
///
/// To use this as a base class for a specific container type, a derived class must implement the
/// following:
/// - template <class Functor, class T>
//    void rawCopyGeneric(Functor&& c, const T& source, MInt begin, MInt end, MInt to)
/// - void invalidate(MInt begin, MInt end)
/// - void reset()
/// - void resize()
/// Optionally, the following methods may be implemented to maintain internal connectivity on
/// delete and move operations:
/// - void moveConnectivity(MInt begin, MInt end, MInt to)
/// - void deleteConnectivity(MInt begin, MInt end)
template <class Derived, template <class> class Invalid>
class Container {
 public:
  // Constructors
  /// Default c'tor does nothing
  constexpr Container() = default;

  // Size/capacity handling
  /// Return capacity (i.e., maximum number of nodes)
  constexpr MInt capacity() const { return m_capacity; }
  void reset(const MInt capacity);
  void resize(const MInt capacity);
  /// Return size (i.e., currently used number of nodes)
  constexpr MInt size() const { return m_size; }
  void size(const MInt size_);

  // High-level operations to modify the container contents.
  void append(const MInt count);
  void append() { append(1); }
  void shrink(const MInt count);
  void shrink() { shrink(1); }
  template <class T>
  void copy(const T& source, const MInt begin, const MInt end, const MInt to);
  template <class T>
  void copy(const T& source, const MInt from, const MInt to) {
    copy(source, from, from + 1, to);
  }
  void copy(const MInt begin, const MInt end, const MInt to) { copy(derived(), begin, end, to); }
  void copy(const MInt from, const MInt to) { copy(from, from + 1, to); }
  void move(const MInt begin, const MInt end, const MInt to);
  void move(const MInt from, const MInt to) { move(from, from + 1, to); }
  void swap(const MInt a, const MInt b);
  void insert(const MInt begin, const MInt count);
  void insert(const MInt id) { insert(id, 1); }
  void erase(const MInt begin, const MInt end);
  void erase(const MInt id) { erase(id, id + 1); }
  void removeAndShift(const MInt begin, const MInt end);
  void removeAndShift(const MInt id) { removeAndShift(id, id + 1); }
  void removeAndFill(const MInt begin, const MInt end);
  void removeAndFill(const MInt id) { removeAndFill(id, id + 1); }
  void clear();

 protected:
  // Memory management
  template <class T>
  using Storage = std::vector<T>;
  template <class T>
  void resetStorage(const MInt n, Storage<T>& c);
  template <class T>
  void resizeStorage(const MInt n, Storage<T>& c);

  virtual void resize() { TERMM(1, "not implemented"); };

  // Fill data container with invalid values
  template <typename Container_, typename T = typename Container_::value_type>
  void fill_invalid(Container_& c, const MInt begin, const MInt end, const MInt solverSize = 1,
                    const T value = Invalid<T>::value()) {
    std::fill(c.data() + begin * solverSize, c.data() + end * solverSize, value);
  }

  /// Copy [begin, end) range with given solver size from source to dest position of target
  template <typename Container_, typename Functor>
  void copyData(const Container_& source, Container_& target, Functor&& f, const MInt begin, const MInt end,
                const MInt dest, const MInt solverSize = 1) {
    f(source.data() + begin * solverSize, source.data() + end * solverSize, target.data() + dest * solverSize);
  }

  // private:
 public:
  // CRTP accessors
  Derived& derived() { return static_cast<Derived&>(*this); }
  const Derived& derived() const { return static_cast<const Derived&>(*this); }

  // Types
  template <MBool left>
  struct Copy {
    template <class It1, class It2, MBool l = left>
    typename std::enable_if<l, It2>::type operator()(It1 first, It1 last, It2 dest) {
      return std::copy(first, last, dest);
    }
    template <class It1, class It2, MBool l = left>
    typename std::enable_if<!l, It2>::type operator()(It1 first, It1 last, It2 dest) {
      return std::copy_backward(first, last, dest);
    }
  };

  // Low-level modifications
  template <class T>
  void rawCopy(const T& source, const MInt begin, const MInt end, const MInt to);
  template <class T>
  void rawCopy(const T& source, const MInt from, const MInt to) {
    rawCopy(source, from, from + 1, to);
  }
  void deleteConnectivity(const MInt NotUsed(begin), const MInt NotUsed(end)) {}
  void moveConnectivity(const MInt NotUsed(begin), const MInt NotUsed(end), const MInt NotUsed(to)) {}
  void moveConnectivity(const MInt from, const MInt to) { moveConnectivity(from, from + 1, to); }
  constexpr MInt dummy() const { return m_capacity; }

 protected:
  // Sanity checking
  MBool isValidId(const MInt id) const;

 private:
  // Capacity/size
  MInt m_capacity = 0;
  MInt m_size = 0;
};


/// Reset tree, re-create data structures with given capacity, and set size to zero.
template <class Derived, template <class> class Invalid>
void Container<Derived, Invalid>::reset(const MInt capacity_) {
  ENSURE_CONDITION(capacity_ >= 0, "Capacity must be non-negative");

  m_capacity = capacity_;
  m_size = 0;

  derived().reset();
}


/// Resize the container capacity.
template <class Derived, template <class> class Invalid>
void Container<Derived, Invalid>::resize(const MInt capacity_) {
  ENSURE_CONDITION(capacity_ >= 0, "Capacity must be non-negative");
  ENSURE_CONDITION(capacity_ >= m_capacity, "So far only resize to a larger size is tested.");

  if(capacity_ == m_capacity) return;

  m_capacity = capacity_;

  derived().resize();
}


/// Resize tree WITHOUT CONSIDERING ANY NODE CONSISTENCY! Use at own risk and remove ASAP...
template <class Derived, template <class> class Invalid>
void Container<Derived, Invalid>::size(const MInt size_) {
  ENSURE_CONDITION(size_ >= 0, "Size must be non-negative");
  ENSURE_CONDITION(size_ <= capacity(), "Size must not exceed capacity");

  m_size = size_;
}


/// Append nodes to end of tree
template <class Derived, template <class> class Invalid>
void Container<Derived, Invalid>::append(const MInt count) {
  ENSURE_CONDITION(count >= 0, "Count must be non-negative");

  if(size() + count > capacity()) {
    std::cerr << "Container capacity exceeded: is:" + std::to_string(capacity())
                     + ", desired: " + std::to_string(size() + count)
              << std::endl;
    // Allow containers to grow if required -> for now resize with +10% capacity
    const MInt newCapacity = 1.1 * (size() + count);
    std::cerr << "New container capacity is: " << newCapacity << std::endl;
    resize(size() + count);
  }
  derived().invalidate(size(), size() + count);
  m_size += count;
}


/// Remove nodes from end of tree
template <class Derived, template <class> class Invalid>
void Container<Derived, Invalid>::shrink(const MInt count) {
  ENSURE_CONDITION(count >= 0, "Count must be non-negative");
  ENSURE_CONDITION(size() - count >= 0, "New size below zero");

  removeAndShift(size() - count, size());
}


/// Copy nodes to another location without changing any parent/child/neighbor information.
template <class Derived, template <class> class Invalid>
template <class T>
void Container<Derived, Invalid>::copy(const T& source, const MInt begin, const MInt end, const MInt to) {
  ENSURE_CONDITION(begin >= 0 && begin < source.size(), "Begin position outside valid range");
  ENSURE_CONDITION(end - 1 < source.size(), "End position outside valid range");
  ENSURE_VALID_ID(to);
  ENSURE_CONDITION(to + (end - begin) <= size(), "Target range outside valid size");

  // Exit early if there is nothing to do
  if(end <= begin || (&source == this && begin == to)) {
    return;
  }

  rawCopy(source, begin, end, to);
}


/// Move nodes to another location and update parent/child/neighbor information accordingly.
template <class Derived, template <class> class Invalid>
void Container<Derived, Invalid>::move(const MInt begin, const MInt end, const MInt to) {
  ENSURE_VALID_ID(begin);
  ENSURE_VALID_ID(end - 1);
  ENSURE_VALID_ID(to);
  ENSURE_CONDITION(to + (end - begin) - 1 < size(), "Target range outside valid size");

  // Exit early if there is nothing to do
  if(end <= begin || begin == to) {
    return;
  }

  rawCopy(derived(), begin, end, to);
  derived().moveConnectivity(begin, end, to);
  derived().invalidate(begin, end);
}


/// Swap two nodes and update parent/child/neighbor information accordingly.
template <class Derived, template <class> class Invalid>
void Container<Derived, Invalid>::swap(const MInt a, const MInt b) {
  ENSURE_VALID_ID(a);
  ENSURE_VALID_ID(b);

  // Exit early if a and b are the same node
  if(a == b) {
    return;
  }

  // Move a to dummy position
  rawCopy(derived(), a, dummy());
  derived().moveConnectivity(a, dummy());

  // Move b to a
  rawCopy(derived(), b, a);
  derived().moveConnectivity(b, a);

  // Move from dummy position to b
  rawCopy(derived(), dummy(), b);
  derived().moveConnectivity(dummy(), b);

  // Invalidate dummy to be sure
  derived().invalidate(dummy(), dummy() + 1);
}


/// Insert 'count' nodes and push back existing nodes while updating parent/child/neighbor
/// information.
template <class Derived, template <class> class Invalid>
void Container<Derived, Invalid>::insert(const MInt id, const MInt count) {
  ENSURE_CONDITION(isValidId(id) || id == size(), "Id must refer to a valid node or be equal to size()");
  ENSURE_CONDITION(count >= 0, "Count must be non-negative");
  ENSURE_CONDITION(count + size() <= capacity(), "New size exceeds capacity");

  // Exit early if there is nothing to do
  if(count == 0) {
    return;
  }

  m_size += count;
  move(id, size() - count, id + count);
}


/// Erase nodes in range [begin, end) and update parent/child/neighbor information.
template <class Derived, template <class> class Invalid>
void Container<Derived, Invalid>::erase(const MInt begin, const MInt end) {
  ENSURE_VALID_ID(begin);
  ENSURE_VALID_ID(end - 1);

  // Exit early if there is nothing to do
  if(end <= begin) {
    return;
  }

  derived().deleteConnectivity(begin, end);
  derived().invalidate(begin, end);
}


/// Remove nodes in range [begin, end) and pull forward existing nodes while updating
/// parent/child/neighbor information.
template <class Derived, template <class> class Invalid>
void Container<Derived, Invalid>::removeAndShift(const MInt begin, const MInt end) {
  ENSURE_VALID_ID(begin);
  ENSURE_VALID_ID(end - 1);

  // Exit early if there is nothing to do
  if(end <= begin) {
    return;
  }

  derived().deleteConnectivity(begin, end);
  derived().invalidate(begin, end);
  if(end < size()) {
    move(end, size(), begin);
  }
  const MInt count = end - begin;
  m_size -= count;
}


/// Remove nodes in range [begin, end) and fill gap with nodes from end while updating
/// parent/child/neighbor information.
template <class Derived, template <class> class Invalid>
void Container<Derived, Invalid>::removeAndFill(const MInt begin, const MInt end) {
  ENSURE_VALID_ID(begin);
  ENSURE_VALID_ID(end - 1);

  // Exit early if there is nothing to do
  if(end <= begin) {
    return;
  }

  derived().deleteConnectivity(begin, end);
  derived().invalidate(begin, end);
  const MInt count = end - begin;
  if(end < size()) {
    move(std::max(end, size() - count), size(), begin);
  }
  m_size -= count;
}


/// Clear tree by invalidating all nodes and setting size to zero
template <class Derived, template <class> class Invalid>
void Container<Derived, Invalid>::clear() {
  derived().invalidate(0, size());
  m_size = 0;
}


/// Copy range of nodes [begin, end) to range starting at 'to'.
template <class Derived, template <class> class Invalid>
template <class T>
void Container<Derived, Invalid>::rawCopy(const T& source, const MInt begin, const MInt end, const MInt to) {
  // Use different methods for overlapping and non-overlapping ranges, depending on overlap type
  if(to < begin || to >= end) {
    // Non-overlapping ranges or overlap while copying left can be handled by std::copy
    derived().rawCopyGeneric(Copy<true>{}, source, begin, end, to);
  } else {
    // Overlap while copying right can be handled by std::copy_backward
    const MInt count = end - begin;
    derived().rawCopyGeneric(Copy<false>{}, source, begin, end, to + count);
  }
}


/// Create new container with given size and replace original one.
template <class Derived, template <class> class Invalid>
template <class T>
void Container<Derived, Invalid>::resetStorage(const MInt n, Storage<T>& c) {
  // Note: Containers are always initialized with an invalid value as for, e.g., std::vector, all
  // entries are default initialized anyways and not providing an invalid value would zero all int's
  // and probably all float's as well.
  Storage<T>(std::max((m_capacity + 1) * n, 1), Invalid<T>::value()).swap(c);
  // TODO: This step is critical for OpenMP parallelization. In most numa
  // architecture data is placed by a 'first touch'-policy.
  // Storage(==std::vector) does a sequential initialization, i.e., it touches
  // all memory from the core associated with the main thread. That might
  // introduce memory access penalty in OpenMP usage.
}


/// Resize container with given size.
template <class Derived, template <class> class Invalid>
template <class T>
void Container<Derived, Invalid>::resizeStorage(const MInt n, Storage<T>& c) {
  c.resize(std::max((m_capacity + 1) * n, 1), Invalid<T>::value());
}


/// Return whether given id refers to a valid node (auxiliary method).
template <class Derived, template <class> class Invalid>
MBool Container<Derived, Invalid>::isValidId(const MInt id) const {
  // @memSplitGrid Note: allow second last position to be used as a dummy id (needed for FV solver to
  // work atm, changing to the last position requires changing the allocation of solver variables to
  // be of size maxNoCells+1)
  if(!((id >= 0 && id < size()) || id == capacity() - 1)) {
    std::cerr << "id = " << id << std::endl
              << "size() = " << size() << std::endl
              << "capacity() = " << capacity() << std::endl;
  }
  return ((id >= 0 && id < size()) || id == capacity() - 1);
}

} // namespace container
} // namespace maia


// Undefine macros that should not be used outside this file
#undef CONTAINER_SANITY_CHECKS
#undef ENSURE_VALID_ID
#undef ENSURE_CONDITION

#endif // ifndef CONTAINER_H_
