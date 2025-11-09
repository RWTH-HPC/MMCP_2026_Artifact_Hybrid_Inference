// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef PARALLELFOR
#define PARALLELFOR

#include "INCLUDE/maiatypes.h"
#include "compiler_config.h"

#ifdef MAIA_PSTL
#include <algorithm>
#include <execution>
#ifdef MAIA_NVHPC_COMPILER
// WAR: https://nvbugs/3285841
#include <thrust/iterator/counting_iterator.h>
#else
#include <cstddef>
#include <iterator>
#endif
#endif

#if !defined(CHUNK_SIZE)
#define CHUNK_SIZE (4096)
#endif

namespace maia {

#if defined(MAIA_PSTL) && !defined(MAIA_NVHPC_COMPILER)

/** \brief  Dummy iterator class
 *  \author Miro Gondrum
 *  \date   02.12.2021
 *
 * This class is of forward iterator type without really iterating over an STL.
 * It is only holding an integer to be in-/decreased.
 * Purpose: To provide a range for std::for_each without explicitly creating a
 * dummy array only containing a certain integer range.
 *
 * \note The usage for other purpose might be critical as it is not pointing to
 *       a real container.
 */
struct RangeIterator {
  using iterator_category = std::forward_iterator_tag;
  using difference_type = std::ptrdiff_t;
  using value_type = MInt;
  using pointer = MInt*;
  using reference = MInt&;

  RangeIterator(value_type value) : m_value(value) {}

  // FIXME: This might be critical as (pointer const ..) is cast to (pointer ..)
  reference operator*() const {
    pointer p((pointer)&m_value);
    return *p;
  }
  pointer operator->() { return &m_value; }

  // Pre-/Postfix increment
  RangeIterator& operator++() {
    m_value++;
    return *this;
  }
  RangeIterator operator++(MInt) {
    RangeIterator tmp = *this;
    ++(*this);
    return tmp;
  }
  // Pre-/Postfix decrement
  RangeIterator& operator--() {
    m_value--;
    return *this;
  }
  RangeIterator operator--(MInt) {
    RangeIterator tmp = *this;
    --(*this);
    return tmp;
  }

  // some operator
  friend MBool operator==(const RangeIterator& a, const RangeIterator& b) { return a.m_value == b.m_value; };
  friend MBool operator!=(const RangeIterator& a, const RangeIterator& b) { return a.m_value != b.m_value; };
  friend MBool operator>(const RangeIterator& a, const RangeIterator& b) { return a.m_value > b.m_value; };

 private:
  MInt m_value;
};

#endif // defined(MAIA_PSTL) && !defined(MAIA_NVHPC_COMPILER)

/** \brief  Wrapper function for parallel for loop (no PSTL)
 *  \author Miro Gondrum
 *  \date   21.02.2022
 *  \note   PLEASE USE parallelFor(..) function
 */
template <class UnaryFunction>
inline void parallelFor_base(MInt begin, MInt end, UnaryFunction&& f) {
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, CHUNK_SIZE) default(none) shared(begin, end, f)
#endif
  for(MInt i = begin; i < end; i++) {
    f(i);
  }
}

/** \brief  Wrapper function for parallel for loop (PSTL)
 *  \author Miro Gondrum
 *  \date   21.02.2022
 *  \note   PLEASE USE parallelFor(..) function
 */
template <class UnaryFunction>
inline void parallelFor_pstl(MInt begin, MInt end, UnaryFunction&& f) {
#if defined(MAIA_PSTL)
#if defined(MAIA_NVHPC_COMPILER)
  // WAR: https://nvbugs/3285841
  // TODO labels:gpu So far this only works with begin == 0. The
  // RangeIterator is not working for nvhpc's pstl implementation, yet.
  auto begin_ = thrust::counting_iterator(MInt{begin});
  auto end_ = end;
#else
  auto begin_ = RangeIterator(begin);
  auto end_ = RangeIterator(end);
#endif
  // TODO miro: GCC: How to trigger the usage of more threads? Currently in my
  // case it using only 1 thread. Hence, OpenMP is performing better
  std::for_each_n(std::execution::par_unseq, begin_, end_, f);
#else /* defined(MAIA_PSTL) */
  parallelFor_base(begin, end, f);
#endif
}

/** \brief  Wrapper function for parallel for loop
 *  \author Miro Gondrum
 *  \date   02.12.2021
 *
 *  \param[in]  begin start value of iteration
 *  \param[in]  end   end value of iteration
 *  \param[in]  f     function which represents the loops body
 *
 * This function wraps a for loop, such that it is performed parallel depending
 * on the given compiler flags. This loop has to be of following type:
 *  for( int i = begin; i < end ; i++) { f(i); }
 */
template <MBool portedToGpu = false, class UnaryFunction>
inline void parallelFor(MInt begin, MInt end, UnaryFunction&& f) {
  if constexpr(portedToGpu) {
    parallelFor_pstl(begin, end, f);
  } else {
    parallelFor_base(begin, end, f);
  }
}

/** \brief  Wrapper function for parallel for loop (STL container based, no PSTL)
 *  \author Miro Gondrum
 *  \date   21.02.2022
 *  \note   PLEASE USE parallelFor(..) function
 */
template <class UnaryFunction, class T>
inline void parallelFor_base(const std::vector<T>& container, UnaryFunction&& f) {
  const MInt end = container.size();
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, CHUNK_SIZE) default(none) shared(end, f, container)
#endif
  for(MInt i = 0; i < end; i++) {
    f(container[i]);
  }
}

/** \brief  Wrapper function for parallel for loop (STL container based, PSTL)
 *  \author Miro Gondrum
 *  \date   21.02.2022
 *  \note   PLEASE USE parallelFor(..) function
 */
template <class UnaryFunction, class T>
inline void parallelFor_pstl(const std::vector<T>& container, UnaryFunction&& f) {
#if defined(MAIA_PSTL)
  std::for_each_n(std::execution::par_unseq, container.begin(), container.end(), f);
#else
  parallelFor_base(container, f);
#endif
}

/** \brief  Wrapper function for parallel for loop (STL container based)
 *  \author Miro Gondrum
 *  \date   18.02.2022
 *
 *  \param[in]  container STL container over whose entries are looped over
 *  \param[in]  f         function which represents the loops body
 *
 * This function wraps a for loop, such that it is performed parallel depending
 * on the given compiler flags. This loop has to be of following type:
 *  for( auto& item: container) { f(item); }
 */
template <MBool portedToGpu = false, class UnaryFunction, class T>
inline void parallelFor(const std::vector<T>& container, UnaryFunction&& f) {
  if constexpr(portedToGpu) {
    parallelFor_pstl(container, f);
  } else {
    parallelFor_base(container, f);
  }
}


/** \brief  Wrapper function for parallel nested for loops [NON-PSTL]
 *  \author Marian Albers, Miro Gondrum
 *  \date   07.03.2023
 *  \note   PLEASE USE parallelFor(..) function
 */
template <MInt nDim, class UnaryFunction>
inline void parallelFor_base(std::array<MInt, nDim> begin, std::array<MInt, nDim> end, UnaryFunction&& f) {
  if constexpr(nDim == 3) {
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, CHUNK_SIZE) default(none) shared(begin, end, f) collapse(3)
#endif
    for(MInt k = begin[2]; k < end[2]; k++) {
      for(MInt j = begin[1]; j < end[1]; j++) {
        for(MInt i = begin[0]; i < end[0]; i++) {
          f(i, j, k);
        }
      }
    }
  } else if constexpr(nDim == 2) {
#if defined(_OPENMP)
#pragma omp parallel for schedule(static, CHUNK_SIZE) default(none) shared(begin, end, f) collapse(2)
#endif
    for(MInt j = begin[1]; j < end[1]; j++) {
      for(MInt i = begin[0]; i < end[0]; i++) {
        f(i, j);
      }
    }
  } else {
    mTerm(1, AT_, "Only nDim==2 and nDim==3 supported");
  }
}

/** \brief  Wrapper function for parallel nested for loops [PSTL]
 *  \author Marian Albers, Miro Gondrum
 *  \date   07.03.2023
 *  \note   PLEASE USE parallelFor(..) function
 */
template <MInt nDim, class UnaryFunction>
inline void parallelFor_pstl(std::array<MInt, nDim> begin, std::array<MInt, nDim> end, UnaryFunction&& f) {
#if defined(MAIA_PSTL)
  std::array<MInt, nDim> size{};
  const MInt beginI = 0;
  MInt endI = 1;
  for(MInt dim = 0; dim < nDim; ++dim) {
    size[dim] = end[dim] - begin[dim];
    endI *= size[dim];
  }
#if defined(MAIA_NVHPC_COMPILER)
  // WAR: https://nvbugs/3285841
  // TODO labels:gpu So far this only works with begin == 0. The
  // RangeIterator is not working for nvhpc's pstl implementation, yet.
  auto begin_ = thrust::counting_iterator(MInt{beginI});
  auto end_ = endI;
#else
  auto begin_ = RangeIterator(beginI);
  auto end_ = RangeIterator(endI);
#endif

  if constexpr(nDim == 3) {
    std::for_each_n(std::execution::par_unseq, begin_, end_, [=](auto& I) {
      const MInt k = (I / (size[0] * size[1])) + begin[2];
      const MInt j = ((I - k * size[0] * size[1]) / size[0]) + begin[1];
      const MInt i = (I % size[0]) + begin[0];

      f(i, j, k);
    });
  } else if constexpr(nDim == 2) {
    std::for_each_n(std::execution::par_unseq, begin_, end_, [=](auto& I) {
      const MInt j = (I / size[0]) + begin[1];
      const MInt i = (I % size[0]) + begin[0];

      f(i, j);
    });
  }
#else /* defined(MAIA_PSTL) */
  parallelFor_base<nDim>(begin, end, f);
#endif
}


/** \brief  Wrapper function for parallel nested for loops
 *  \author Marian Albers, Miro Gondrum
 *  \date   07.03.2023
 *
 *  \param[in]  begin[nDim] start value of iteration
 *  \param[in]  end[nDim]  end value of iteration
 *  \param[in]  f     function which represents the loops body
 *
 * This function wraps a one or multiple nested for loops,
 * such that it is performed parallel depending
 * on the given compiler flags. This loop has to be of following type:
 *  for( int i = begin[0]; i < end[0] ; i++) { f(i); }
 * or
 *  for( int i = begin[0]; i < end[0] ; i++) {
 *    for( int j = begin[1]; j < end[1] ; j++) {
 *      f(i,j);
 *    }
 *  }
 */
template <MBool portedToGpu = false, MInt nDim, class UnaryFunction>
inline void parallelFor(std::array<MInt, nDim> begin, std::array<MInt, nDim> end, UnaryFunction&& f) {
  if constexpr(portedToGpu) {
    parallelFor_pstl<nDim>(begin, end, f);
  } else {
    parallelFor_base<nDim>(begin, end, f);
  }
}

} // namespace maia
#endif /* PARALLELFOR */
