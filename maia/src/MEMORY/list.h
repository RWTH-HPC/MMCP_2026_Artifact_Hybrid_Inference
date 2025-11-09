// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef LIST_H
#define LIST_H

#include <iostream>
#include <sstream>
#include "INCLUDE/maiatypes.h"
#include "UTIL/functions.h"

template <typename T>
class List {
 public:
  List(MLong maxSize);
  ~List() { delete[] a; };
  MInt size() { return (MInt)m_size; };

  MInt maxSize() { return (MInt)m_maxSize; };

  void append();
  MLong memoryUseage() { return (MLong)((m_maxSize + 1) * sizeof(T)); };

  T* operator[](const MInt index);
  T* a = nullptr;
  MLong m_rawMemoryCounter;

  MInt setSize(MInt inputSize) {
    if(inputSize < m_size) m_size = inputSize;
    return (MInt)inputSize;
  }

  MInt resetSize(MInt inputSize) {
    if(inputSize < 0) {
      mTerm(1, AT_, "Input size is < 0!");
    }
    if(inputSize > m_maxSize) {
      std::stringstream errorMessage;
      errorMessage << " Error in list, maxSize reached ( " << m_maxSize << " cells ).";
      mTerm(1, AT_, errorMessage.str());
    }
    m_size = inputSize;
    return (MInt)inputSize;
  }

  char* getRawPointer();

 private:
  MLong m_maxSize;
  MInt m_size;
};

template <typename T>
List<T>::List(MLong inputMaxSize) {
  m_size = 0;
  m_rawMemoryCounter = 0;
  m_maxSize = inputMaxSize;

  a = new T[(inputMaxSize + 1)];
}

template <typename T>
void List<T>::append() {
  if(m_size < m_maxSize)
    m_size++;
  else {
    std::stringstream errorMessage;
    errorMessage << " Error in list, maxSize reached ( " << m_maxSize << " cells ).";
    mTerm(1, AT_, errorMessage.str());
  }
}


template <typename T>
T* List<T>::operator[](const MInt index) {
  return &a[index];
}
#endif
