// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef TENSOR_H_
#define TENSOR_H_

#include <ostream>
#include "INCLUDE/maiatypes.h"
#include "debug.h"
#include "functions.h"

namespace maia {
/**
 * \brief Namespace that contains all classes, functions and constants needed for Tensor.
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
 * \date 2012-11-13
 */
namespace tensor {

// Forward declaration of TensorStorage
namespace detail_ {
template <class T>
class TensorStorage;
}

/**
 * \brief Provides a lightweight and fast class for accessing 1D arrays as multi-dimensional tensors (up to 4D).
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
 * \date 2012-11-09
 *
 * \tparam T Scalar data type of the tensor elements.
 *
 * \details Use with care when doing performance-sensitive work as this has not yet been sufficiently tested for speed.
 */
template <class T>
class Tensor {
 public:
  // Public types
  typedef MLong size_type;
  typedef T value_type;
  typedef value_type* pointer;
  typedef const value_type* const_pointer;
  typedef value_type& reference;
  typedef const value_type& const_reference;
  typedef pointer iterator;
  typedef const_pointer const_iterator;
  typedef MInt dim_type;

  // Constructors & destructor
  Tensor();
  explicit Tensor(size_type n0, size_type n1 = 1, size_type n2 = 1, size_type n3 = 1, size_type n4 = 1);
  Tensor(T* data, size_type n0, size_type n1 = 1, size_type n2 = 1, size_type n3 = 1, size_type n4 = 1);
  ~Tensor();

  size_type size() const;
  size_type resize(size_type n0, size_type n1 = 1, size_type n2 = 1, size_type n3 = 1, size_type n4 = 1);
  void assign(size_type n, const T& newValue);
  void set(const T& value);
  void clear();
  void transpose();

  void swap(Tensor<T>& other);
  template <class TT>
  friend void swap(Tensor<TT>& a, Tensor<TT>& b);

  size_type dim(size_type d) const;
  size_type dim0() const;
  size_type dim1() const;
  size_type dim2() const;
  size_type dim3() const;
  size_type dim4() const;

  Tensor<T>& operator=(Tensor<T> copy);
  Tensor<T>(Tensor<T>&& copy) = default;
  Tensor<T>(const Tensor<T>& copy) = default;


  T& operator[](size_type index);
  const T& operator[](size_type index) const;
  T& operator()(size_type i0);
  const T& operator()(size_type i0) const;
  T& operator()(size_type i0, size_type i1);
  const T& operator()(size_type i0, size_type i1) const;
  T& operator()(size_type i0, size_type i1, size_type i2);
  const T& operator()(size_type i0, size_type i1, size_type i2) const;
  T& operator()(size_type i0, size_type i1, size_type i2, size_type i3);
  const T& operator()(size_type i0, size_type i1, size_type i2, size_type i3) const;
  T& operator()(size_type i0, size_type i1, size_type i2, size_type i3, size_type i4);
  const T& operator()(size_type i0, size_type i1, size_type i2, size_type i3, size_type i4) const;

 private:
  static const dim_type m_maxNoDims = 5;
  typedef detail_::TensorStorage<T> DataStorage;

  DataStorage m_data;
  size_type m_size;
  size_type m_dim[m_maxNoDims]{};
};


/**
 * \brief Default constructor does nothing but setting all internal member variables to zero.
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
 * \date 2012-11-09
 */
template <class T>
inline Tensor<T>::Tensor() : m_data(), m_size(0) {
  m_dim[0] = 0;
  m_dim[1] = 0;
  m_dim[2] = 0;
  m_dim[3] = 0;
  m_dim[4] = 0;
}

/**
 * \brief Initializes the Tensor object by setting the data pointer and the internal dimensions.
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
 * \date 2012-11-09
 *
 * \param[in] data Pointer to previously allocated 1D array.
 * \param[in] n0 Size of first dimension.
 * \param[in] n1 Size of second dimension.
 * \param[in] n2 Size of third dimension.
 * \param[in] n3 Size of fourth dimension.
 * \param[in] n4 Size of fifth dimension.
 */
template <class T>
inline Tensor<T>::Tensor(size_type n0, size_type n1, size_type n2, size_type n3, size_type n4) {
  ASSERT(n0 > 0, "Dimension n0 must be greater than zero");
  ASSERT(n1 > 0, "Dimension n1 must be greater than zero");
  ASSERT(n2 > 0, "Dimension n2 must be greater than zero");
  ASSERT(n3 > 0, "Dimension n3 must be greater than zero");
  ASSERT(n4 > 0, "Dimension n4 must be greater than zero");

  m_dim[0] = n0;
  m_dim[1] = n1;
  m_dim[2] = n2;
  m_dim[3] = n3;
  m_dim[4] = n4;

  m_size = n0 * n1 * n2 * n3 * n4;
  m_data.resize(m_size);
}


/**
 * \brief Initializes the Tensor object by setting the data pointer and the internal dimensions.
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
 * \date 2012-11-09
 *
 * \param[in] data Pointer to previously allocated 1D array.
 * \param[in] n0 Size of first dimension.
 * \param[in] n1 Size of second dimension.
 * \param[in] n2 Size of third dimension.
 * \param[in] n3 Size of fourth dimension.
 * \param[in] n4 Size of fifth dimension.
 */
template <class T>
inline Tensor<T>::Tensor(T* data, size_type n0, size_type n1, size_type n2, size_type n3, size_type n4) {
  ASSERT(n0 > 0, "Dimension n0 must be greater than zero");
  ASSERT(n1 > 0, "Dimension n1 must be greater than zero");
  ASSERT(n2 > 0, "Dimension n2 must be greater than zero");
  ASSERT(n3 > 0, "Dimension n3 must be greater than zero");
  ASSERT(n4 > 0, "Dimension n4 must be greater than zero");

  m_dim[0] = n0;
  m_dim[1] = n1;
  m_dim[2] = n2;
  m_dim[3] = n3;
  m_dim[4] = n4;

  m_size = n0 * n1 * n2 * n3 * n4;
  m_data.resize(m_size, data);
}

/**
 * \brief Resets all internal pointers and variables to zero.
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
 * \date 2012-11-09
 */
template <class T>
inline Tensor<T>::~Tensor() {
  clear();
}


/**
 * \brief Returns the size of the array as product of all five dimensions (i.e. not the actual array size but the size
 * as seen by the Tensor object).
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
 * \date 2012-11-09
 *
 * \return Total size of the data array as seen by the Tensor object.
 */
template <class T>
inline typename Tensor<T>::size_type Tensor<T>::size() const {
  return m_size;
}


/**
 * \brief Deletes the old data structure and creates a new one with the requested dimensions.
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
 * \date 2012-11-12
 *
 * \tparam T
 * \param[in] n0
 * \param[in] n1
 * \param[in] n2
 * \param[in] n3
 * \param[in] n4
 *
 * \return The new size of the data structure, i.e. the number of elements in the tensor.
 *
 * \details Note: After a call to resize the tensor uses an internal storage, even if
 *          it was constructed with an external data pointer.
 */
template <class T>
inline typename Tensor<T>::size_type Tensor<T>::resize(size_type n0, size_type n1, size_type n2, size_type n3,
                                                       size_type n4) {
  m_dim[0] = n0;
  m_dim[1] = n1;
  m_dim[2] = n2;
  m_dim[3] = n3;
  m_dim[4] = n4;

  m_size = n0 * n1 * n2 * n3 * n4;
  m_data.resize(size());

  return size();
}


/**
 * \brief Resizes the tensor to a 1D vector of length 'n' and with initial value 'value'.
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
 * \date 2012-11-13
 *
 * \tparam T
 * \param[in] n
 * \param[in] value
 */
template <class T>
inline void Tensor<T>::assign(size_type n, const T& value) {
  resize(n);
  set(value);
}


/**
 * \brief Initializes tensor to constant value.
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
 * \date 2012-11-13
 *
 * \param[in] value Constant value that is assigned to all elements.
 */
template <class T>
inline void Tensor<T>::set(const T& value) {
  std::fill_n(&m_data[0], size(), value);
}


/**
 * \brief Deletes the data (if internal storage was used) and resets dimensions to zero.
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
 * \date 2012-11-13
 */
template <class T>
inline void Tensor<T>::clear() {
  m_dim[0] = 0;
  m_dim[1] = 0;
  m_dim[2] = 0;
  m_dim[3] = 0;
  m_dim[4] = 0;

  m_size = 0;
  m_data.clear();
}


/**
 * \brief Transposes the first two dimensions (i.e. acting as if the tensor is a matrix).
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
 * \date 2012-11-13
 *
 * \details Note: This creates a temporary copy of the WHOLE tensor, even though only a copy of the
 *          first to dimensions would be needed. Thus this very expensive if transposing a higher-
 *          order tensor.
 */
template <class T>
inline void Tensor<T>::transpose() {
  Tensor<T> tmp(*this);

  std::swap(m_dim[0], m_dim[1]);

  for(MInt i = 0; i < tmp.dim0(); i++) {
    for(MInt j = 0; j < tmp.dim1(); j++) {
      (*this)(j, i) = tmp(i, j);
    }
  }
}


/**
 * \brief Efficiently swap this tensor with another one (uses copying of pointers, no data is moved).
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
 * \date 2012-11-13
 *
 * \tparam T Scalar dat atype of the tensor elements.
 * \param[] other The tensor that should be swapped with.
 */
template <class T>
inline void Tensor<T>::swap(Tensor<T>& other) {
  using std::swap;
  swap(m_data, other.m_data);
  swap(m_size, other.m_size);
  for(dim_type i = 0; i < m_maxNoDims; ++i)
    swap(m_dim[i], other.m_dim[i]);
}


/**
 * \brief Return the size of dimension d.
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
 * \date 2012-11-09
 *
 * \param[in] d The requested dimension.
 *
 * \return Size of dimension d.
 */
template <class T>
inline typename Tensor<T>::size_type Tensor<T>::dim(size_type d) const {
  ASSERT(d >= 0, "Dimension d must be greater than or equal to zero!");
  ASSERT(d < m_maxNoDims, "Dimension d must be less than m_maxNoDims!");

  return m_dim[d];
}


/**
 * \brief Return the size of dimension 0.
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
 * \date 2012-11-09
 *
 * \return Size of dimension 0.
 */
template <class T>
inline typename Tensor<T>::size_type Tensor<T>::dim0() const {
  return m_dim[0];
}


/**
 * \brief Return the size of dimension 1.
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
 * \date 2012-11-09
 *
 * \return Size of dimension 1.
 */
template <class T>
inline typename Tensor<T>::size_type Tensor<T>::dim1() const {
  return m_dim[1];
}


/**
 * \brief Return the size of dimension 2.
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
 * \date 2012-11-09
 *
 * \return Size of dimension 2.
 */
template <class T>
inline typename Tensor<T>::size_type Tensor<T>::dim2() const {
  return m_dim[2];
}


/**
 * \brief Return the size of dimension 3.
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
 * \date 2012-11-09
 *
 * \return Size of dimension 3.
 */
template <class T>
inline typename Tensor<T>::size_type Tensor<T>::dim3() const {
  return m_dim[3];
}


/**
 * \brief Return the size of dimension 4.
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
 * \date 2012-11-09
 *
 * \return Size of dimension 4.
 */
template <class T>
inline typename Tensor<T>::size_type Tensor<T>::dim4() const {
  return m_dim[4];
}


/**
 * \brief Assignment operator.
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
 * \date 2012-11-13
 *
 * \param[in] copy The tensor that is assigned from.
 *
 * \return This tensor with the new data.
 */
template <class T>
inline Tensor<T>& Tensor<T>::operator=(Tensor<T> copy) {
  swap(copy);

  return *this;
}


/**
 * \brief Returns a reference to the element at position index (non-const version).
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
 * \date 2012-11-14
 *
 * \param[in] index Position of the element in the data storage.
 *
 * \return The element at the specified position.
 *
 * \details This provides an array-like access to the tensor, which is especially useful
 *          if the tensor is used as a vector.
 */
template <class T>
inline T& Tensor<T>::operator[](size_type index) {
  ASSERT(index >= 0 && index < size(), "Index out of bounds");

  return m_data[index];
}


/**
 * \brief Returns a reference to the element at position index (const version).
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
 * \date 2012-11-14
 *
 * \param[in] index Position of the element in the data storage.
 *
 * \return The element at the specified position.
 *
 * \details This provides an array-like access to the tensor, which is especially useful
 *          if the tensor is used as a vector.
 */
template <class T>
inline const T& Tensor<T>::operator[](size_type index) const {
  ASSERT(index >= 0 && index < size(), "Index out of bounds");

  return m_data[index];
}


/**
 * \brief Provides tensor access element using index notation (1D, non-const version).
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
 * \date 2012-11-09
 *
 * \param[in] i0 Index of first dimension.
 *
 * \return The array element at the specified index.
 */
template <class T>
inline T& Tensor<T>::operator()(size_type i0) {
  ASSERT(i0 >= 0 && i0 < dim0(), "Index i0 out of bounds");

  return m_data[i0];
}


/**
 * \brief Provides tensor element access using index notation (1D, const version).
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
 * \date 2012-11-09
 *
 * \param[in] i0 Index of first dimension.
 *
 * \return The array element at the specified index.
 */
template <class T>
inline const T& Tensor<T>::operator()(size_type i0) const {
  ASSERT(i0 >= 0 && i0 < dim0(), "Index i0 out of bounds");

  return m_data[i0];
}


/**
 * \brief Provides tensor access element using index notation (2D, non-const version).
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
 * \date 2012-11-09
 *
 * \param[in] i0 Index of first dimension.
 * \param[in] i1 Index of second dimension.
 *
 * \return The array element at the specified index.
 */
template <class T>
inline T& Tensor<T>::operator()(size_type i0, size_type i1) {
  ASSERT(i0 >= 0 && i0 < dim0(), "Index i0 out of bounds");
  ASSERT(i1 >= 0 && i1 < dim1(), "Index i1 out of bounds");

  return m_data[i0 * dim1() + i1];
}


/**
 * \brief Provides tensor element access using index notation (2D, const version).
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
 * \date 2012-11-09
 *
 * \param[in] i0 Index of first dimension.
 * \param[in] i1 Index of second dimension.
 *
 * \return The array element at the specified index.
 */
template <class T>
inline const T& Tensor<T>::operator()(size_type i0, size_type i1) const {
  ASSERT(i0 >= 0 && i0 < dim0(), "Index i0 out of bounds");
  ASSERT(i1 >= 0 && i1 < dim1(), "Index i1 out of bounds");

  return m_data[i0 * dim1() + i1];
}


/**
 * \brief Provides tensor access element using index notation (3D, non-const version).
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
 * \date 2012-11-09
 *
 * \param[in] i0 Index of first dimension.
 * \param[in] i1 Index of second dimension.
 * \param[in] i2 Index of third dimension.
 *
 * \return The array element at the specified index.
 */
template <class T>
inline T& Tensor<T>::operator()(size_type i0, size_type i1, size_type i2) {
  ASSERT(i0 >= 0 && i0 < dim0(), "Index i0 out of bounds");
  ASSERT(i1 >= 0 && i1 < dim1(), "Index i1 out of bounds");
  ASSERT(i2 >= 0 && i2 < dim2(), "Index i2 out of bounds");

  return m_data[i0 * dim1() * dim2() + i1 * dim2() + i2];
}


/**
 * \brief Provides tensor element access using index notation (3D, const version).
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
 * \date 2012-11-09
 *
 * \param[in] i0 Index of first dimension.
 * \param[in] i1 Index of second dimension.
 * \param[in] i2 Index of third dimension.
 *
 * \return The array element at the specified index.
 */
template <class T>
inline const T& Tensor<T>::operator()(size_type i0, size_type i1, size_type i2) const {
  ASSERT(i0 >= 0 && i0 < dim0(), "Index i0 out of bounds");
  ASSERT(i1 >= 0 && i1 < dim1(), "Index i1 out of bounds");
  ASSERT(i2 >= 0 && i2 < dim2(), "Index i2 out of bounds");

  return m_data[i0 * dim1() * dim2() + i1 * dim2() + i2];
}


/**
 * \brief Provides tensor access element using index notation (4D, non-const version).
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
 * \date 2012-11-09
 *
 * \param[in] i0 Index of first dimension.
 * \param[in] i1 Index of second dimension.
 * \param[in] i2 Index of third dimension.
 * \param[in] i3 Index of fourth dimension.
 *
 * \return The array element at the specified index.
 */
template <class T>
inline T& Tensor<T>::operator()(size_type i0, size_type i1, size_type i2, size_type i3) {
  ASSERT(i0 >= 0 && i0 < dim0(), "Index i0 out of bounds");
  ASSERT(i1 >= 0 && i1 < dim1(), "Index i1 out of bounds");
  ASSERT(i2 >= 0 && i2 < dim2(), "Index i2 out of bounds");
  ASSERT(i3 >= 0 && i3 < dim3(), "Index i3 out of bounds");

  return m_data[i0 * dim1() * dim2() * dim3() + i1 * dim2() * dim3() + i2 * dim3() + i3];
}


/**
 * \brief Provides tensor element access using index notation (4D, const version).
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
 * \date 2012-11-09
 *
 * \param[in] i0 Index of first dimension.
 * \param[in] i1 Index of second dimension.
 * \param[in] i2 Index of third dimension.
 * \param[in] i3 Index of fourth dimension.
 *
 * \return The array element at the specified index.
 */
template <class T>
inline const T& Tensor<T>::operator()(size_type i0, size_type i1, size_type i2, size_type i3) const {
  ASSERT(i0 >= 0 && i0 < dim0(), "Index i0 out of bounds");
  ASSERT(i1 >= 0 && i1 < dim1(), "Index i1 out of bounds");
  ASSERT(i2 >= 0 && i2 < dim2(), "Index i2 out of bounds");
  ASSERT(i3 >= 0 && i3 < dim3(), "Index i3 out of bounds");

  return m_data[i0 * dim1() * dim2() * dim3() + i1 * dim2() * dim3() + i2 * dim3() + i3];
}


/**
 * \brief Provides tensor access element using index notation (5D, non-const version).
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
 * \date 2012-11-09
 *
 * \param[in] i0 Index of first dimension.
 * \param[in] i1 Index of second dimension.
 * \param[in] i2 Index of third dimension.
 * \param[in] i3 Index of fourth dimension.
 * \param[in] i4 Index of fifth dimension.
 *
 * \return The array element at the specified index.
 */
template <class T>
inline T& Tensor<T>::operator()(size_type i0, size_type i1, size_type i2, size_type i3, size_type i4) {
  ASSERT(i0 >= 0 && i0 < dim0(), "Index i0 out of bounds");
  ASSERT(i1 >= 0 && i1 < dim1(), "Index i1 out of bounds");
  ASSERT(i2 >= 0 && i2 < dim2(), "Index i2 out of bounds");
  ASSERT(i3 >= 0 && i3 < dim3(), "Index i3 out of bounds");
  ASSERT(i4 >= 0 && i4 < dim4(), "Index i4 out of bounds");

  return m_data[i0 * dim1() * dim2() * dim3() * dim4() + i1 * dim2() * dim3() * dim4() + i2 * dim3() * dim4()
                + i3 * dim4() + i4];
}


/**
 * \brief Provides tensor element access using index notation (5D, const version).
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
 * \date 2012-11-09
 *
 * \param[in] i0 Index of first dimension.
 * \param[in] i1 Index of second dimension.
 * \param[in] i2 Index of third dimension.
 * \param[in] i3 Index of fourth dimension.
 * \param[in] i4 Index of fifth dimension.
 *
 * \return The array element at the specified index.
 */
template <class T>
inline const T& Tensor<T>::operator()(size_type i0, size_type i1, size_type i2, size_type i3, size_type i4) const {
  ASSERT(i0 >= 0 && i0 < dim0(), "Index i0 out of bounds");
  ASSERT(i1 >= 0 && i1 < dim1(), "Index i1 out of bounds");
  ASSERT(i2 >= 0 && i2 < dim2(), "Index i2 out of bounds");
  ASSERT(i3 >= 0 && i3 < dim3(), "Index i3 out of bounds");
  ASSERT(i4 >= 0 && i4 < dim4(), "Index i4 out of bounds");

  return m_data[i0 * dim1() * dim2() * dim3() * dim4() + i1 * dim2() * dim3() * dim4() + i2 * dim3() * dim4()
                + i3 * dim4() + i4];
}


/**
 * \brief Prints the values of the Tensor object together with the respective index.
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
 * \date 2012-11-09
 *
 * \tparam TT Scalar data type of the tensor elements.
 * \param[in] os Output stream that should be written to.
 * \param[in] indexer The indexer object.
 *
 * \return The output stream for further concatenation.
 */
template <class TT>
std::ostream& operator<<(std::ostream& os, const Tensor<TT>& t) {
  for(MInt i = 0; i < t.dim0(); i++) {
    for(MInt j = 0; j < t.dim1(); j++) {
      for(MInt k = 0; k < t.dim2(); k++) {
        for(MInt l = 0; l < t.dim3(); l++) {
          os << "(" << i << "," << j << "," << k << "," << l << ") = ";
          os << t(i, j, k, l) << std::endl;
        }
      }
    }
  }

  os << "Dimensions (i,j,k,l) = (" << t.dim0() << "," << t.dim1() << "," << t.dim2() << "," << t.dim3() << ")"
     << std::endl;

  return os;
}


/**
 * \brief Non-member swap exchanges the contents of two Tensors.
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
 * \date 2012-11-13
 *
 * \tparam TT Scalar data type of the tensor elements.
 * \param[in] a The first tensor to be swapped.
 * \param[in] b The second tensor to be swapped.
 */
template <class TT>
inline void swap(Tensor<TT>& a, Tensor<TT>& b) {
  a.swap(b);
}


////////////////////////////////////////////////////////////////////////////////


// Internal namespace for methods/class that are generally not needed outside of Tensor
namespace detail_ {

/**
 * \brief A vector-like data storage class that supports both internal and external storage mechanisms.
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
 * \date 2012-11-14
 *
 * \tparam T Scalar data type of the tensor elements.
 *
 * \details By default, an internal storage mechanism is used, i.e. memory is allocated and deallocated
 *          automatically. However, there is also the possibility to provide a pointer to an existing
 *          storage location, which allows this class to be used as a "viewport" to data in memory. For
 *          the design of the class, a close look was taken at std::vector, in order to support most
 *          of the dynamic reallocation methods and other operations that std::vector has as well.\n\n
 *
 *          This class is mainly used by Tensor as its internal storage.
 */
template <class T>
class TensorStorage {
 public:
  // Public types
  typedef MLong size_type;
  typedef T value_type;
  typedef value_type* pointer;
  typedef const value_type* const_pointer;
  typedef value_type& reference;
  typedef const value_type& const_reference;
  typedef pointer iterator;
  typedef const_pointer const_iterator;

  // Constructors & destructor
  TensorStorage();
  TensorStorage(const TensorStorage<T>& ts);
  explicit TensorStorage(size_type n);
  TensorStorage(size_type n, T* data);
  ~TensorStorage();

  size_type size() const;
  void resize(size_type n);
  void resize(size_type n, T* data);

  void clear();
  void swap(TensorStorage<T>& other);
  template <class TT>
  friend void swap(TensorStorage<TT>& a, TensorStorage<TT>& b);

  TensorStorage<T>& operator=(TensorStorage<T> copy);
  T& operator[](size_type n);
  const T& operator[](size_type n) const;

 private:
  MBool isAllocated();
  void reallocate(size_type n);

 private:
  T* m_data;
  size_type m_size;
  MBool m_isAllocated;
};


/**
 * \brief Default constructor does not allocate any memory but creates a zero-sized storage container.
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
 * \date 2012-11-14
 */
template <class T>
inline TensorStorage<T>::TensorStorage() : m_data(0), m_size(0), m_isAllocated(false) {
  // Nothing here
}


/**
 * \brief Copy constructor creates a copy of the provided object.
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
 * \date 2012-11-14
 *
 * \param[in] ts Object that is copied.
 *
 * \details Note: The copy constructor will always allocate memory, even if the original object
 *          used external storage.
 */
template <class T>
inline TensorStorage<T>::TensorStorage(const TensorStorage<T>& ts) : m_data(0), m_size(0), m_isAllocated(false) {
  // Allocate new memory and copy the data members
  reallocate(ts.size());
  std::copy(&ts.m_data[0], &ts.m_data[0] + ts.size(), m_data);
}


/**
 * \brief Creates a new TensorStorage object using internal storage.
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
 * \date 2012-11-14
 *
 * \param[in] n Size of the new container.
 */
template <class T>
inline TensorStorage<T>::TensorStorage(size_type n) : m_data(0), m_size(0), m_isAllocated(false) {
  resize(n);
}


/**
 * \brief Creates a new TensorStorage object using external storage, i.e. no memory is allocated.
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
 * \date 2012-11-14
 *
 * \param[in] n Size of the storage in memory.
 * \param[in] data Pointer to existing memory location.
 *
 * \details Since this is supposed to be a very lightweight and fast way to access memory through
 *          TensorStorage, no data is changed to some initial value when this constructor is called.
 *          Furthermore, the user has to ensure that the data location that is pointed to does not get
 *          deallocated before the TensorStorage object is destroyed.
 */
template <class T>
inline TensorStorage<T>::TensorStorage(size_type n, T* data) : m_data(data), m_size(n), m_isAllocated(false) {
  ASSERT(n >= 0, "New storage size must be greater or equal to zero.");
  ASSERT(data != 0, "Data pointer may not be a null pointer.");
}


/**
 * \brief The destructor only calls clear().
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
 * \date 2012-11-14
 */
template <class T>
inline TensorStorage<T>::~TensorStorage() {
  clear();
}


/**
 * \brief Returns the current size of the data container.
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
 * \date 2012-11-14
 *
 * \return Length of the data array.
 */
template <class T>
inline typename TensorStorage<T>::size_type TensorStorage<T>::size() const {
  return m_size;
}


/**
 * \brief Resizes the container by first dropping all existing elements and then reallocating new memory.
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
 * \date 2012-11-14
 *
 * \param[in] n New size of the data container.
 */
template <class T>
inline void TensorStorage<T>::resize(size_type n) {
  ASSERT(n >= 0, "New storage size must be greater or equal to zero.");

  reallocate(n);
}


/**
 * \brief Resizes the container by first dropping all existing elements and then using the external
 *        data pointer for the new storage.
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
 * \date 2012-11-14
 *
 * \param[in] n New size of the data container.
 * \param[in] data Pointer to existing memory location.
 */
template <class T>
inline void TensorStorage<T>::resize(size_type n, T* data) {
  ASSERT(n >= 0, "New storage size must be greater or equal to zero.");
  ASSERT(data != 0, "Data pointer may not be a null pointer.");

  clear();
  m_size = n;
  m_data = data;
}


/**
 * \brief Clears the data storage by deallocating any previously allocated memory and setting the size to zero.
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
 * \date 2012-11-14
 */
template <class T>
inline void TensorStorage<T>::clear() {
  if(isAllocated()) {
    delete[] m_data;
  }

  m_data = 0;
  m_size = 0;
  m_isAllocated = false;
}


/**
 * \brief Swap the contents of the existing TensorStorage object with another one.
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
 * \date 2012-11-14
 *
 * \param[in] other The object with which the contents are exchanged.
 */
template <class T>
inline void TensorStorage<T>::swap(TensorStorage<T>& other) {
  using std::swap;
  swap(m_data, other.m_data);
  swap(m_size, other.m_size);
  swap(m_isAllocated, other.m_isAllocated);
}


/**
 * \brief Swap the contents between two TensorStorage objects.
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
 * \date 2012-11-14
 *
 * \tparam TT
 * \param[in] a First object.
 * \param[in] b Second object.
 */
template <class TT>
inline void swap(TensorStorage<TT>& a, TensorStorage<TT>& b) {
  a.swap(b);
}


/**
 * \brief Assignment operator copies another object to the current one.
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
 * \date 2012-11-14
 *
 * \param[in] copy The object that is copied.
 *
 * \return The assigned object (*this)
 */
template <class T>
inline TensorStorage<T>& TensorStorage<T>::operator=(TensorStorage<T> copy) {
  swap(copy);

  return *this;
}


/**
 * \brief Returns a reference to the element at position n (non-const version).
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
 * \date 2012-11-14
 *
 * \param[in] n Position of the element in the data storage.
 *
 * \return The element at the specified position.
 */
template <class T>
inline T& TensorStorage<T>::operator[](size_type n) {
  ASSERT(n >= 0 && n < size(), "Index out of bounds");

  return m_data[n];
}


/**
 * \brief Returns a reference to the element at position n (const version).
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
 * \date 2012-11-14
 *
 * \param[in] n Position of the element in the data storage.
 *
 * \return The element at the specified position.
 */
template <class T>
inline const T& TensorStorage<T>::operator[](size_type n) const {
  ASSERT(n >= 0 && n < size(), "Index out of bounds");

  return m_data[n];
}


/**
 * \brief Returns true if this class has allocated its own memory, and false if external memory is used.
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
 * \date 2013-05-08
 */
template <class T>
inline MBool TensorStorage<T>::isAllocated() {
  return m_isAllocated;
}


/**
 * \brief Clears any existing internal data storage and then allocates new memory to store as many elements as
 * specified.
 *
 * \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
 * \date 2012-11-14
 *
 * \param[in] n The number of elements that can be stored in the new container.
 */
template <class T>
[[gnu::noinline]] // this is neccesary for flto compilation
void TensorStorage<T>::reallocate(size_type n) {
  clear();

  m_data = (n <= 0) ? nullptr : new T[std::min(n, PTRDIFF_MAX)];
  m_size = n;
  m_isAllocated = true;
}

} // namespace detail_


} // namespace tensor
} // namespace maia

// Define some types for ease of use
typedef maia::tensor::Tensor<MFloat> MFloatVector;
typedef maia::tensor::Tensor<MInt> MIntVector;
typedef maia::tensor::Tensor<MInt> MIntVector;
typedef maia::tensor::Tensor<MFloat> MFloatMatrix;
typedef maia::tensor::Tensor<MInt> MIntMatrix;
typedef maia::tensor::Tensor<MInt> MIntMatrix;
typedef maia::tensor::Tensor<MFloat> MFloatTensor;
typedef maia::tensor::Tensor<MInt> MIntTensor;
typedef maia::tensor::Tensor<MInt> MIntTensor;

#endif /* TENSOR_H_ */
