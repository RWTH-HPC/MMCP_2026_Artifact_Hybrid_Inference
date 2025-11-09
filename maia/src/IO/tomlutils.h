// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef TOMLUTILS_H_
#define TOMLUTILS_H_

#include <memory>
#include <utility>
#include "INCLUDE/maiamacro.h"
#include "INCLUDE/maiatypes.h"
#include "enums.h"

// Disable RTTI and exceptions
#define CPPTOML_NO_RTTI
#define CPPTOML_NO_EXCEPTIONS
#include "cpptoml.h"

namespace maia {
namespace io {
namespace toml {

/// Type traits for enum type
template <VariableType e>
struct Enum2Type;
template <>
struct Enum2Type<MSTRING> {
  using type = MString;
};
template <>
struct Enum2Type<MINT> {
  using type = MInt;
};
template <>
struct Enum2Type<MFLOAT> {
  using type = MFloat;
};
template <>
struct Enum2Type<MBOOL> {
  using type = MBool;
};
template <class T>
struct TypeTraits;
template <>
struct TypeTraits<MString> {
  static const VariableType type = MSTRING;
  static MString name() { return "string"; }
};
template <>
struct TypeTraits<MInt> {
  static const VariableType type = MINT;
  static MString name() { return "int"; }
};
template <>
struct TypeTraits<MFloat> {
  static const VariableType type = MFLOAT;
  static MString name() { return "float"; }
};
template <>
struct TypeTraits<MBool> {
  static const VariableType type = MBOOL;
  static MString name() { return "bool"; }
};


/// Class that represents a single key-value pair for TOML properties
class Property {
 public:
  // C'tors
  Property() = default;
  Property(MString name_, std::vector<MString> data)
    : m_name(std::move(name_)), m_type(MSTRING), m_size(data.size()), m_string(std::move(data)) {}
  Property(MString name_, std::vector<MInt> data)
    : m_name(std::move(name_)), m_type(MINT), m_size(data.size()), m_int(std::move(data)) {}
  Property(MString name_, std::vector<MFloat> data)
    : m_name(std::move(name_)), m_type(MFLOAT), m_size(data.size()), m_float(std::move(data)) {}
  Property(MString name_, std::vector<MBool> data)
    : m_name(std::move(name_)), m_type(MBOOL), m_size(data.size()), m_bool(std::move(data)) {}

  // General member functions
  MString name() const { return m_name; }
  VariableType type() const { return m_type; }
  MString type2string() const;
  MLong size() const { return m_size; }
  MBool valid() const { return m_size != -1; }

  // Accessors
  const std::vector<MString>& asString() const {
    if(type() != MSTRING) {
      TERMM(1, "bad type");
    }
    return m_string;
  }
  const std::vector<MInt>& asInt() const {
    if(type() != MINT) {
      TERMM(1, "bad type");
    }
    return m_int;
  }
  const std::vector<MFloat>& asFloat() const {
    if(type() != MFLOAT) {
      TERMM(1, "bad type");
    }
    return m_float;
  }
  const std::vector<MBool>& asBool() const {
    if(type() != MBOOL) {
      TERMM(1, "bad type");
    }
    return m_bool;
  }

 private:
  MString m_name = "";
  VariableType m_type = MINVALID;
  MLong m_size = -1;
  std::vector<MString> m_string;
  std::vector<MInt> m_int;
  std::vector<MFloat> m_float;
  std::vector<MBool> m_bool;
};


inline MString Property::type2string() const {
  switch(type()) {
    case MSTRING:
      return TypeTraits<Enum2Type<MSTRING>::type>::name();
    case MINT:
      return TypeTraits<Enum2Type<MINT>::type>::name();
    case MFLOAT:
      return TypeTraits<Enum2Type<MFLOAT>::type>::name();
    case MBOOL:
      return TypeTraits<Enum2Type<MBOOL>::type>::name();
    default:
      TERMM(1, "bad value type: " + std::to_string(type()));
  }
}

/// Obtain type information from TOML value
inline VariableType value2type(const std::shared_ptr<cpptoml::base>& value) {
  if(value->as<std::string>()) {
    return MSTRING;
  } else if(value->as<int64_t>()) {
    return MINT;
  } else if(value->as<double>()) {
    return MFLOAT;
  } else if(value->as<bool>()) {
    return MBOOL;
  } else {
    TERMM(1, "bad value type");
  }
}


/// Create property from cpptoml item
inline Property makeProperty(const MString& name, const std::shared_ptr<cpptoml::base>& item) {
  if(item->is_value()) {
    switch(value2type(item)) {
      case MSTRING: {
        std::vector<MString> data(1);
        data[0] = item->as<std::string>()->get();
        return Property(name, data);
      }
      case MINT: {
        std::vector<MInt> data(1);
#if defined(MAIA_GCC_COMPILER)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wnull-dereference"
#endif
        data[0] = static_cast<MInt>(item->as<int64_t>()->get());
#if defined(MAIA_GCC_COMPILER)
#pragma GCC diagnostic pop
#endif
        return Property(name, data);
      }
      case MFLOAT: {
        std::vector<MFloat> data(1);
        data[0] = item->as<double>()->get();
        return Property(name, data);
      }
      case MBOOL: {
        std::vector<MBool> data(1);
        data[0] = item->as<bool>()->get();
        return Property(name, data);
      }
      default:
        TERMM(1, "bad value type");
    }
  } else if(item->is_array() || item->is_table()) {
    if(item->as_array()->get().empty()) {
      TERMM(1, "array of size 0 found, cannot detect type");
    }

    switch(value2type(item->as_array()->at(0))) {
      case MSTRING: {
        auto v = *item->as_array()->get_array_of<std::string>();
        std::vector<MString> data(v.size());
        copy(v.begin(), v.end(), data.begin());
        return Property(name, data);
      }
      case MINT: {
        auto v = *item->as_array()->get_array_of<MLong>();
        std::vector<MInt> data(v.size());
        copy(v.begin(), v.end(), data.begin());
        return Property(name, data);
      }
      case MFLOAT: {
        auto v = *item->as_array()->get_array_of<double>();
        std::vector<MFloat> data(v.size());
        copy(v.begin(), v.end(), data.begin());
        return Property(name, data);
      }
      case MBOOL: {
        auto v = *item->as_array()->get_array_of<bool>();
        std::vector<MBool> data(v.size());
        copy(v.begin(), v.end(), data.begin());
        return Property(name, data);
      }
      default:
        TERMM(1, "bad value type");
    }
  } else {
    TERMM(1, "item is not a value or array");
  }
}

/// Non-recursive search for solver aliases
inline void collectSolverAliases(const std::shared_ptr<cpptoml::table>& tab,
                                 std::map<std::string, int>& solverAliases) {
  // search for key "solverAlias"
  for(auto&& item : *tab) {
    auto&& key = item.first;
    if(key == "solverAlias") {
      auto&& solverAliasesTable = item.second;
      // iterate over all items in solverAlias table
      if(solverAliasesTable->is_table()) {
        for(auto&& alias : *solverAliasesTable->as_table()) {
          // check if solverAlias is given as string
          if(value2type(alias.second) == MSTRING) {
            // store solver alias in map (solverAlias <string> -> solverid <int> )
            solverAliases[alias.second->as<std::string>()->get()] = std::stoi(alias.first);
          } else {
            TERMM(1, "solverAlias need to be defined as string!");
          }
        }
      } else {
        TERMM(1, "solverAlias needs to be in table format: solverAlias.solverId = 'alias'!");
      }
    }
  }
}

/// Recursively traverse TOML table and collect all properties with name, type, and count
inline void collectProperties(const std::shared_ptr<cpptoml::table>& tab,
                              std::vector<std::string>& names,
                              std::vector<Property>& properties,
                              std::map<std::string, int>& solverAliases) {
  //  // TODO labels:IO,toremove remove
  //  for (auto&& item : *tab) {
  //    // Store key, value for convenience
  //    auto&& key = item.first;
  //    auto&& value = item.second;
  //    std::cerr << "key " << key << " value " << value << std::endl;
  //  }

  // Iterate over table
  for(auto&& item : *tab) {
    // Store key, value for convenience
    auto&& key = item.first;
    auto&& value = item.second;
    std::stringstream key_ss;
    key_ss << key;

    if(key_ss.str() == "default") {
      // if default is the key we just add name to properties
      properties.emplace_back(makeProperty(names[0], value));
    } else if(cpptoml::is_number(key_ss.str().c_str()[0]) && !value->is_table()) {
      std::stringstream k;
      k << names[0] << "." << (std::atoi(key.c_str()));
      properties.emplace_back(makeProperty(k.str(), value));
    } else if(std::find_if(solverAliases.begin(), solverAliases.end(),
                           [&key_ss, &solverAliases](std::pair<const std::string, int>& entry) {
                             return !solverAliases.empty() && (entry.first == key_ss.str());
                           })
                  != solverAliases.end()
              && !names.empty()) {
      std::stringstream k;
      k << names[0] << "." << solverAliases[key_ss.str()];
      properties.emplace_back(makeProperty(k.str(), value));
    } else if(value->is_value() || value->is_array()) {
      // Determine fully qualified name
      std::stringstream solverId;
      std::stringstream k;
      MInt counter = 0;
      for(auto&& name : names) {
        solverId << name << ".";
        counter++;
      }

      if(solverId.str().empty()) {
        k << key;
      } else {
        MString subString1 = solverId.str().substr(0, solverId.str().size() - 1);
        MString subString2 = solverId.str().substr(0, solverId.str().find('.'));
        if(subString2 == "solver" || subString1 == "solver") {
          if(strstr(key.c_str(), ".") != nullptr) {
            std::stringstream error;
            error << "Setting a solverId inside the solver table environment '" << subString1 << "' is prohibited!!! ";
            TERMM(1, error.str());
          }
          switch(counter) {
            case 1: {
              k << key;
              break;
            }
            case 2: {
              const char* du = std::strrchr(subString1.c_str(), '.') + 1;
              const MInt singleSolverId = std::atoi(du);
              k << key << "." << singleSolverId;
              break;
            }
            default: {
              std::stringstream error;
              error << "Too many nested tables: " << subString1 << " Please use 'solver.solverId' instead!!!";
              TERMM(1, error.str());
            }
          }
        } else {
          std::stringstream error;
          error << "Unknown table type: " << subString2 << " Please use 'solver' instead!!!";
          TERMM(1, error.str());
        }
      }

      // Add info to collected properties
      properties.emplace_back(makeProperty(k.str(), value));
    } else if(value->is_table()) {
      // If value is a table, add key to list of names and descend one recursion level
      names.push_back(key);
      collectProperties(value->as_table(), names, properties, solverAliases);
      names.pop_back();
    } else {
      TERMM(1, "only values, arrays, and tables are supported");
    }
  }
}


/// Helper function to not have to provide a names array oneself
inline void collectProperties(const std::shared_ptr<cpptoml::table>& tab, std::vector<Property>& properties,
                              std::map<std::string, int>& solverAliases) {
  std::vector<std::string> names;
  collectProperties(tab, names, properties, solverAliases);
}

/// Helper function to not have to provide a names array oneself
inline void collectProperties(const std::shared_ptr<cpptoml::table>& tab, std::vector<Property>& properties) {
  std::vector<std::string> names;
  std::map<std::string, int> aliasesdummy;
  collectProperties(tab, names, properties, aliasesdummy);
}

} // namespace toml
} // namespace io
} // namespace maia

#endif // ifndef TOMLUTILS_H_
