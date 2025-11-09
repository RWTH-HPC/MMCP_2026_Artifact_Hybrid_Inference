// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef FACTOR_H_
#define FACTOR_H_

#include <functional>
#include <map>
#include <memory>
#include "INCLUDE/maiatypes.h"
#include "UTIL/functions.h"


/// \brief Class implementing the object factory pattern.
///
/// \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
/// \date 2014-04-05
///
/// \tparam AbstractProduct The (abstract) base class of the objects to be
///                         created.
/// \tparam IdentifierType The type to be used to identify derived types.
/// \tparam ReturnType Return type of the create() method.
/// \tparam ProductCreator The type of the creating function.
/// \tparam Args Additional (optional) arguments that should be passed to the
///              constructor of the created objects.
///
/// This follow the factory pattern as described in Alexandrescu (2001).
///
/// References:
///   Andrei Alexandrescu (2001): Modern C++ Design.
template <class AbstractProduct, typename IdentifierType, typename ReturnType = std::unique_ptr<AbstractProduct>,
          class ProductCreator = std::function<ReturnType()>, class... Args>
class MFactory {
 public:
  /// \brief Function to add new types to the factory.
  ///
  /// \param[in] id Identifier for newly added type.
  /// \param[in] creator Function returning a new instance of the new type.
  ///
  /// \return True if the type was sucessfully added.
  MBool add(const IdentifierType& id, ProductCreator creator) {
    return m_assoc.insert(typename AssocMap::value_type(id, creator)).second;
  }

  /// \brief Function to remove existing types from the factory.
  ///
  /// \param[in] id Identifier for the type to be removed.
  ///
  /// \return True if the type was sucessfully removed.
  MBool remove(const IdentifierType& id) { return m_assoc.erase(id) == true; }

  /// \brief Function to create a new instance of a type.
  ///
  /// \param[in] id Identifier of the type to be created.
  /// \param[in] args Optional arguments to be passed to the type's constructor.
  ///
  /// \return A handle to the newly created object.
  ReturnType create(const IdentifierType& id, Args... args) const {
    auto i = m_assoc.find(id);
    if(i != m_assoc.end()) {
      return (i->second)(args...);
    } else {
      TERMM(1, "Identifier not found.");
    }
  }

 private:
  typedef std::map<IdentifierType, ProductCreator> AssocMap;
  AssocMap m_assoc;
};

#endif // FACTOR_H_
