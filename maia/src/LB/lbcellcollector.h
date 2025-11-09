// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef LBCOLLECTOR_H_
#define LBCOLLECTOR_H_

#include <algorithm>
#include <bitset>
#include <type_traits>
#include <vector>
#include "INCLUDE/maiamacro.h"
#include "INCLUDE/maiatypes.h"
#include "MEMORY/container.h"
#include "UTIL/functions.h"
#include "compiler_config.h"
#include "lbcellproperties.h"

// The following macro enables the "Structure-of-Arrays" memory layout for multi-dimensional node
// variables. This might be beneficial for GPU computations. Default is "Array-of-Structures".
// Examples (for nodes nN with four children cM each)
// Array-of-Structures (AOS): n0c0, n0c1, n0c2, n0c3, n1c0, n1c1, n1c2, n1c3, n2c0, n2c1, ...
// Structure-of-Arrays (SOA): n0c0, n1c0, n2c0, n3c0, ..., n0c1, n1c1, n2c1, n3c1, ..., n0c2, ...
#ifdef WAR_NVHPC_PSTL
#define LBCOLLECTOR_SOA_MEMORY_LAYOUT
#endif

// The macro 'LBCOLLECTOR_SANITY_CHECKS_ACCESSORS' enables (potentially very expensive) sanity
// checks for all accessors. It is enabled for build type "extra_debug".
#ifdef MAIA_EXTRA_DEBUG
#define LBCOLLECTOR_SANITY_CHECKS_ACCESSORS
#endif

// Sanity-checking macros for accessors
#if defined(LBCOLLECTOR_SANITY_CHECKS_ACCESSORS) || defined(MAIA_ASSERT_ACCESSORS)
#define ENSURE_VALID_ID_ACCESSOR(id)                                                                                   \
  do {                                                                                                                 \
    MAIA_CONTAINER_ENSURE_VALID_ID(id);                                                                                \
  } while(false)
#define ENSURE_VALID_VARIABLE_ID_ACCESSOR(id)                                                                          \
  do {                                                                                                                 \
    MAIA_CONTAINER_ENSURE(                                                                                             \
        id >= 0 && id < noVariables(),                                                                                 \
        "variable id = " + std::to_string(id) + " is out-of-bounds [0, " + std::to_string(noVariables()) + ")", AT_);  \
  } while(false)
#define ENSURE_VALID_DISTRIBUTION_ID_ACCESSOR(id)                                                                      \
  do {                                                                                                                 \
    MAIA_CONTAINER_ENSURE(id >= 0 && id < noDistributions(),                                                           \
                          "distribution id = " + std::to_string(id) + " is out-of-bounds [0, "                         \
                              + std::to_string(noDistributions()),                                                     \
                          AT_);                                                                                        \
  } while(false)

#define ENSURE_VALID_PROPERTY_ACCESSOR(p)                                                                              \
  do {                                                                                                                 \
    MAIA_CONTAINER_ENSURE(p != LbCell::NumProperties, "Invalid property", AT_);                                        \
  } while(false)
#else
#define ENSURE_VALID_ID_ACCESSOR(id)                                                                                   \
  do {                                                                                                                 \
  } while(false)
#define ENSURE_VALID_VARIABLE_ID_ACCESSOR(id)                                                                          \
  do {                                                                                                                 \
  } while(false)
#define ENSURE_VALID_DISTRIBUTION_ID_ACCESSOR(id)                                                                      \
  do {                                                                                                                 \
  } while(false)
#define ENSURE_VALID_PROPERTY_ACCESSOR(dir)                                                                            \
  do {                                                                                                                 \
  } while(false)
#endif


// Namespace for auxiliary functions/classes
namespace maia {
namespace lb {
namespace collector {

/// Underlying bitset type for property storage
using BitsetType = maia::lb::cell::BitsetType;

// Type traits for invalid values. These values are used to initialize/erase nodes
template <class T>
struct Invalid {};

// Invalid value for ids is 'INT_MIN'
template <>
struct Invalid<MInt> {
  static constexpr MInt value() { return std::numeric_limits<MInt>::min(); }
};

// Invalid value for floats is 'NaN'
template <>
struct Invalid<MFloat> {
  static constexpr MFloat value() {
#ifdef MAIA_PGI_COMPILER
    return std::numeric_limits<MFloat>::quiet_NaN();
#else
    return std::numeric_limits<MFloat>::signaling_NaN();
#endif
  }
};

// Invalid value for BitsetProperties is '0'
template <>
struct Invalid<BitsetType> {
  static constexpr BitsetType value() { return BitsetType(0); }
};

/// Class that represents LB cell collector.
template <MInt nDim>
class LbCellCollector : public maia::container::Container<LbCellCollector<nDim>, Invalid> {
  // Necessary for CRTP
  friend class maia::container::Container<LbCellCollector<nDim>, Invalid>;

  // Make base class functions known to use without this pointer
  using Base = maia::container::Container<LbCellCollector<nDim>, Invalid>;
  using Base::resetStorage;
  template <class T>
  using Storage = typename Base::template Storage<T>;
  using BitsetType = maia::lb::cell::BitsetType;


 public:
  // Types
  template <class T>
  using Invalid = typename maia::lb::collector::Invalid<T>;

  // Constructors
  /// Default c'tor does nothing
  constexpr LbCellCollector() = default;

  // Ensure that base class method is found when called from outside
  using Base::copyData;
  using Base::fill_invalid;
  using Base::reset;

  // Accessors
  MFloat& nu(const MInt id);
  MFloat nu(const MInt id) const;
  MFloat& nuT(const MInt id);
  MFloat nuT(const MInt id) const;
  MFloat& oldNu(const MInt id);
  MFloat oldNu(const MInt id) const;
  MFloat& oldNuT(const MInt id);
  MFloat oldNuT(const MInt id) const;
  MInt& bndId(const MInt id);
  MInt bndId(const MInt id) const;
  MInt& level(const MInt id);
  MInt level(const MInt id) const;
  MFloat& kappa(const MInt id);
  MFloat kappa(const MInt id) const;
  MFloat& diffusivity(const MInt id);
  MFloat diffusivity(const MInt id) const;
  MFloat& spongeFactor(const MInt id);
  MFloat spongeFactor(const MInt id) const;
  MFloat& variables(const MInt id, const MInt eid);
  MFloat variables(const MInt id, const MInt eid) const;
  MFloat* variables_ptr();
  MFloat& oldVariables(const MInt id, const MInt eid);
  MFloat oldVariables(const MInt id, const MInt eid) const;
  MFloat* oldVariables_ptr();
  MFloat& distributions(const MInt id, const MInt eid);
  MFloat distributions(const MInt id, const MInt eid) const;
  MFloat& oldDistributions(const MInt id, const MInt eid);
  MFloat oldDistributions(const MInt id, const MInt eid) const;
  MFloat& distributionsThermal(const MInt id, const MInt eid);
  MFloat distributionsThermal(const MInt id, const MInt eid) const;
  MFloat& oldDistributionsThermal(const MInt id, const MInt eid);
  MFloat oldDistributionsThermal(const MInt id, const MInt eid) const;
  MFloat& distributionsTransport(const MInt id, const MInt eid);
  MFloat distributionsTransport(const MInt id, const MInt eid) const;
  MFloat& oldDistributionsTransport(const MInt id, const MInt eid);
  MFloat oldDistributionsTransport(const MInt id, const MInt eid) const;
  MFloat& externalForces(const MInt id, const MInt eid);
  MFloat externalForces(const MInt id, const MInt eid) const;
  MFloat& previousDistribution(const MInt id, const MInt eid);
  MFloat previousDistribution(const MInt id, const MInt eid) const;
  MFloat& previousVariable(const MInt id, const MInt eid);
  MFloat previousVariable(const MInt id, const MInt eid) const;
  MFloat& uOtherPhase(const MInt id, const MInt dir);
  MFloat uOtherPhase(const MInt id, const MInt dir) const;
  MFloat& invVolumeFraction(const MInt id);
  MFloat invVolumeFraction(const MInt id) const;

  // Property-related accessors
  BitsetType::reference hasProperty(const MInt id, const LbCell p);
  MBool hasProperty(const MInt id, const LbCell p) const;
  void resetProperties(const MInt id);
  BitsetType& allProperties(const MInt id);

  /// Allow setting whether to support thermal computations
  void setThermal(const MBool isThermal_);

  /// Allow setting whether to support transport computations
  void setTransport(const MBool useTransport_);

  /// Update number of variables according to used model
  void setNoVariables();

  /// Sets the number of distributions
  void setNoDistributions(const MInt noDistributions_);

  /// Return number of species
  constexpr MInt isThermal() const { return m_isThermal; }

  /// Return number of species
  constexpr MInt useTransport() const { return m_useTransport; }

  void setSaveUOtherPhase(const MBool saveUOtherPhase_);
  void setSaveVolumeFraction(const MBool saveVolumeFraction_);
  void setSavePrevVars(const MBool savePrevVars_);
  constexpr MInt saveUOtherPhase() const { return m_saveUOtherPhase; }
  constexpr MInt saveVolumeFraction() const { return m_saveVolumeFraction; }
  constexpr MInt savePrevVars() const { return m_savePrevVars; }
  void setSaveNuT(const MBool saveNuT_);
  void setSaveOldNu(const MBool saveOldNu_);
  constexpr MInt saveNuT() const { return m_saveNuT; }
  constexpr MInt saveOldNu() const { return m_saveOldNu; }

  /// Return number of variables
  constexpr MInt noVariables() const { return m_noVariables; }

  /// Return number of cells
  constexpr MInt noCells() const { return this->size(); }

  /// Return number of distributions
  constexpr MInt noDistributions() const { return m_noDistributions; }

  /// Return number of properties defined for each node
  static constexpr MInt noProperties() { return maia::lb::cell::p(LbCell::NumProperties); }

 private:
  // Methods required by base class for CRTP
  void reset();
  void invalidate(const MInt begin, const MInt end);
  template <class Functor, class T>
  void rawCopyGeneric(Functor&& c, const T& source, const MInt begin, const MInt end, const MInt destination);

  /// Number of variables
  MInt m_noVariables = 1 + nDim;

  MInt m_noDistributions = 0;

  MInt m_noCells = 0;

  /// Use thermal model
  MBool m_isThermal = false;

  /// Use transport model
  MBool m_useTransport = false;

  MBool m_saveUOtherPhase = false;
  MBool m_saveVolumeFraction = false;
  MBool m_savePrevVars = false;
  MBool m_saveNuT = false;
  MBool m_saveOldNu = false;

  // Data containers
  Storage<MFloat> m_nu{};
  Storage<MFloat> m_oldNu{};
  Storage<MFloat> m_kappa{};
  Storage<MFloat> m_diffusivity{};
  Storage<MFloat> m_spongeFactor{};
  Storage<MInt> m_bndId{};
  Storage<MInt> m_level{};
  Storage<MFloat> m_variables{};
  Storage<MFloat> m_oldVariables{};
  Storage<MFloat> m_distributions{};
  Storage<MFloat> m_oldDistributions{};
  Storage<MFloat> m_distributionsThermal{};
  Storage<MFloat> m_oldDistributionsThermal{};
  Storage<MFloat> m_distributionsTransport{};
  Storage<MFloat> m_oldDistributionsTransport{};
  Storage<MFloat> m_externalForces{};
  Storage<MFloat> m_previousDistribution{};
  Storage<MFloat> m_previousVariable{};
  Storage<MFloat> m_nuT{};
  Storage<MFloat> m_oldNuT{};
  Storage<MFloat> m_uOtherPhase{};
  Storage<MFloat> m_invVolumeFraction{};
  Storage<BitsetType> m_properties{};
};

/// Reset tree, re-create data structures with given capacity, and set size to zero.
template <MInt nDim>
void LbCellCollector<nDim>::reset() {
  resetStorage(1, m_nu);
  if(saveOldNu()) {
    resetStorage(1, m_oldNu);
  }
  resetStorage(1, m_kappa);
  resetStorage(1, m_diffusivity);
  resetStorage(1, m_spongeFactor);
  resetStorage(1, m_bndId);
  resetStorage(1, m_level);
  resetStorage(noVariables(), m_variables);
  resetStorage(noVariables(), m_oldVariables);
  resetStorage(noDistributions(), m_distributions);
  resetStorage(noDistributions(), m_oldDistributions);
  resetStorage(nDim, m_externalForces);
  if(isThermal()) {
    resetStorage(noDistributions(), m_distributionsThermal);
    resetStorage(noDistributions(), m_oldDistributionsThermal);
  }
  if(useTransport()) {
    resetStorage(noDistributions(), m_distributionsTransport);
    resetStorage(noDistributions(), m_oldDistributionsTransport);
  }
  if(savePrevVars()) {
    resetStorage(noDistributions(), m_previousDistribution);
    resetStorage(noVariables(), m_previousVariable);
  }
  if(saveNuT()) {
    resetStorage(1, m_nuT);
    if(saveOldNu()) {
      resetStorage(1, m_oldNuT);
    }
  }
  if(saveUOtherPhase()) {
    resetStorage(nDim, m_uOtherPhase);
  }
  if(saveVolumeFraction()) {
    resetStorage(1, m_invVolumeFraction);
  }
  resetStorage(1, m_properties);
}


/// Accessor for nu.
template <MInt nDim>
MFloat& LbCellCollector<nDim>::nu(const MInt id) {
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_nu[id];
}
/// Accessor for nu (const version).
template <MInt nDim>
MFloat LbCellCollector<nDim>::nu(const MInt id) const {
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_nu[id];
}

template <MInt nDim>
MFloat& LbCellCollector<nDim>::oldNu(const MInt id) {
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_oldNu[id];
}
template <MInt nDim>
MFloat LbCellCollector<nDim>::oldNu(const MInt id) const {
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_oldNu[id];
}

/// Accessor for kappa.
template <MInt nDim>
MFloat& LbCellCollector<nDim>::kappa(const MInt id) {
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_kappa[id];
}

/// Accessor for kappa (const version).
template <MInt nDim>
MFloat LbCellCollector<nDim>::kappa(const MInt id) const {
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_kappa[id];
}

/// Accessor for the diffusivity.
template <MInt nDim>
MFloat& LbCellCollector<nDim>::diffusivity(const MInt id) {
  // Prevent accidental compilation without support for SoA layout
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_diffusivity[id];
}

/// Accessor for the diffusivity (const version).
template <MInt nDim>
MFloat LbCellCollector<nDim>::diffusivity(const MInt id) const {
  // Prevent accidental compilation without support for SoA layout
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_diffusivity[id];
}

/// Accessor for spongeFactor.
template <MInt nDim>
MFloat& LbCellCollector<nDim>::spongeFactor(const MInt id) {
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_spongeFactor[id];
}
/// Accessor for spongeFactor (const version).
template <MInt nDim>
MFloat LbCellCollector<nDim>::spongeFactor(const MInt id) const {
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_spongeFactor[id];
}

/// Accessor for bndId.
template <MInt nDim>
MInt& LbCellCollector<nDim>::bndId(const MInt id) {
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_bndId[id];
}
/// Accessor for bndId (const version).
template <MInt nDim>
MInt LbCellCollector<nDim>::bndId(const MInt id) const {
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_bndId[id];
}

/// Accessor for level.
template <MInt nDim>
MInt& LbCellCollector<nDim>::level(const MInt id) {
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_level[id];
}
/// Accessor for level (const version).
template <MInt nDim>
MInt LbCellCollector<nDim>::level(const MInt id) const {
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_level[id];
}

template <MInt nDim>
MFloat& LbCellCollector<nDim>::variables(const MInt id, const MInt eid) {
#ifdef LBCOLLECTOR_SOA_MEMORY_LAYOUT
  return m_variables[eid * noCells() + id];
#else
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_VARIABLE_ID_ACCESSOR(eid);
  return m_variables[id * noVariables() + eid];
#endif
}

template <MInt nDim>
MFloat LbCellCollector<nDim>::variables(const MInt id, const MInt eid) const {
#ifdef LBCOLLECTOR_SOA_MEMORY_LAYOUT
  return m_variables[eid * noCells() + id];
#else
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_VARIABLE_ID_ACCESSOR(eid);
  return m_variables[id * noVariables() + eid];
#endif
}

template <MInt nDim>
MFloat* LbCellCollector<nDim>::variables_ptr() {
  return m_variables.data();
}

template <MInt nDim>
MFloat& LbCellCollector<nDim>::oldVariables(const MInt id, const MInt eid) {
#ifdef LBCOLLECTOR_SOA_MEMORY_LAYOUT
  return m_oldVariables[eid * noCells() + id];
#else
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_VARIABLE_ID_ACCESSOR(eid);
  return m_oldVariables[id * noVariables() + eid];
#endif
}

template <MInt nDim>
MFloat LbCellCollector<nDim>::oldVariables(const MInt id, const MInt eid) const {
#ifdef LBCOLLECTOR_SOA_MEMORY_LAYOUT
  return m_oldVariables[eid * noCells() + id];
#else
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_VARIABLE_ID_ACCESSOR(eid);
  return m_oldVariables[id * noVariables() + eid];
#endif
}

template <MInt nDim>
MFloat* LbCellCollector<nDim>::oldVariables_ptr() {
  return m_oldVariables.data();
}

template <MInt nDim>
MFloat& LbCellCollector<nDim>::distributions(const MInt id, const MInt eid) {
#ifdef LBCOLLECTOR_SOA_MEMORY_LAYOUT
  return m_distributions[eid * noCells() + id];
#else
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_DISTRIBUTION_ID_ACCESSOR(eid);
  return m_distributions[id * noDistributions() + eid];
#endif
}

template <MInt nDim>
MFloat LbCellCollector<nDim>::distributions(const MInt id, const MInt eid) const {
#ifdef LBCOLLECTOR_SOA_MEMORY_LAYOUT
  return m_distributions[eid * noCells() + id];
#else
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_DISTRIBUTION_ID_ACCESSOR(eid);
  return m_distributions[id * noDistributions() + eid];
#endif
}

template <MInt nDim>
MFloat& LbCellCollector<nDim>::oldDistributions(const MInt id, const MInt eid) {
#ifdef LBCOLLECTOR_SOA_MEMORY_LAYOUT
  return m_oldDistributions[eid * noCells() + id];
#else
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_DISTRIBUTION_ID_ACCESSOR(eid);
  return m_oldDistributions[id * noDistributions() + eid];
#endif
}

template <MInt nDim>
MFloat LbCellCollector<nDim>::oldDistributions(const MInt id, const MInt eid) const {
#ifdef LBCOLLECTOR_SOA_MEMORY_LAYOUT
  return m_oldDistributions[eid * noCells() + id];
#else
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_DISTRIBUTION_ID_ACCESSOR(eid);
  return m_oldDistributions[id * noDistributions() + eid];
#endif
}

template <MInt nDim>
MFloat& LbCellCollector<nDim>::previousDistribution(const MInt id, const MInt eid) {
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_DISTRIBUTION_ID_ACCESSOR(eid);
  return m_previousDistribution[id * noDistributions() + eid];
}

template <MInt nDim>
MFloat LbCellCollector<nDim>::previousDistribution(const MInt id, const MInt eid) const {
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_DISTRIBUTION_ID_ACCESSOR(eid);
  return m_previousDistribution[id * noDistributions() + eid];
}

template <MInt nDim>
MFloat& LbCellCollector<nDim>::previousVariable(const MInt id, const MInt eid) {
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_DISTRIBUTION_ID_ACCESSOR(eid);
  return m_previousVariable[id * noVariables() + eid];
}

template <MInt nDim>
MFloat LbCellCollector<nDim>::previousVariable(const MInt id, const MInt eid) const {
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_DISTRIBUTION_ID_ACCESSOR(eid);
  return m_previousVariable[id * noVariables() + eid];
}

template <MInt nDim>
MFloat& LbCellCollector<nDim>::nuT(const MInt id) {
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_nuT[id];
}

template <MInt nDim>
MFloat LbCellCollector<nDim>::nuT(const MInt id) const {
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_nuT[id];
}

template <MInt nDim>
MFloat& LbCellCollector<nDim>::uOtherPhase(const MInt id, const MInt dir) {
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_uOtherPhase[id * nDim + dir];
}

template <MInt nDim>
MFloat LbCellCollector<nDim>::uOtherPhase(const MInt id, const MInt dir) const {
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_uOtherPhase[id * nDim + dir];
}

template <MInt nDim>
MFloat& LbCellCollector<nDim>::invVolumeFraction(const MInt id) {
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_invVolumeFraction[id];
}

template <MInt nDim>
MFloat LbCellCollector<nDim>::invVolumeFraction(const MInt id) const {
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_invVolumeFraction[id];
}

template <MInt nDim>
MFloat& LbCellCollector<nDim>::oldNuT(const MInt id) {
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_oldNuT[id];
}

template <MInt nDim>
MFloat LbCellCollector<nDim>::oldNuT(const MInt id) const {
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_oldNuT[id];
}

template <MInt nDim>
MFloat& LbCellCollector<nDim>::distributionsThermal(const MInt id, const MInt eid) {
#ifdef LBCOLLECTOR_SOA_MEMORY_LAYOUT
  return m_distributionsThermal[eid * noCells() + id];
#else
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_DISTRIBUTION_ID_ACCESSOR(eid);
  return m_distributionsThermal[id * noDistributions() + eid];
#endif
}

template <MInt nDim>
MFloat LbCellCollector<nDim>::distributionsThermal(const MInt id, const MInt eid) const {
#ifdef LBCOLLECTOR_SOA_MEMORY_LAYOUT
  return m_distributionsThermal[eid * noCells() + id];
#else
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_DISTRIBUTION_ID_ACCESSOR(eid);
  return m_distributionsThermal[id * noDistributions() + eid];
#endif
}

template <MInt nDim>
MFloat& LbCellCollector<nDim>::oldDistributionsThermal(const MInt id, const MInt eid) {
#ifdef LBCOLLECTOR_SOA_MEMORY_LAYOUT
  return m_oldDistributionsThermal[eid * noCells() + id];
#else
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_DISTRIBUTION_ID_ACCESSOR(eid);
  return m_oldDistributionsThermal[id * noDistributions() + eid];
#endif
}

template <MInt nDim>
MFloat LbCellCollector<nDim>::oldDistributionsThermal(const MInt id, const MInt eid) const {
#ifdef LBCOLLECTOR_SOA_MEMORY_LAYOUT
  return m_oldDistributionsThermal[eid * noCells() + id];
#else
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_DISTRIBUTION_ID_ACCESSOR(eid);
  return m_oldDistributionsThermal[id * noDistributions() + eid];
#endif
}

template <MInt nDim>
MFloat& LbCellCollector<nDim>::distributionsTransport(const MInt id, const MInt eid) {
#ifdef LBCOLLECTOR_SOA_MEMORY_LAYOUT
  return m_distributionsTransport[eid * noCells() + id];
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_DISTRIBUTION_ID_ACCESSOR(eid);
  return m_distributionsTransport[id * noDistributions() + eid];
}

template <MInt nDim>
MFloat LbCellCollector<nDim>::distributionsTransport(const MInt id, const MInt eid) const {
#ifdef LBCOLLECTOR_SOA_MEMORY_LAYOUT
  return m_distributionsTransport[eid * noCells() + id];
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_DISTRIBUTION_ID_ACCESSOR(eid);
  return m_distributionsTransport[id * noDistributions() + eid];
}

template <MInt nDim>
MFloat& LbCellCollector<nDim>::externalForces(const MInt id, const MInt eid) {
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_DISTRIBUTION_ID_ACCESSOR(eid);
  return m_externalForces[id * nDim + eid];
}

template <MInt nDim>
MFloat LbCellCollector<nDim>::externalForces(const MInt id, const MInt eid) const {
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_DISTRIBUTION_ID_ACCESSOR(eid);
  return m_externalForces[id * nDim + eid];
}

template <MInt nDim>
MFloat& LbCellCollector<nDim>::oldDistributionsTransport(const MInt id, const MInt eid) {
#ifdef LBCOLLECTOR_SOA_MEMORY_LAYOUT
  return m_oldDistributionsTransport[eid * noCells() + id];
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_DISTRIBUTION_ID_ACCESSOR(eid);
  return m_oldDistributionsTransport[id * noDistributions() + eid];
}

template <MInt nDim>
MFloat LbCellCollector<nDim>::oldDistributionsTransport(const MInt id, const MInt eid) const {
#ifdef LBCOLLECTOR_SOA_MEMORY_LAYOUT
  return m_oldDistributionsTransport[eid * noCells() + id];
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_DISTRIBUTION_ID_ACCESSOR(eid);
  return m_oldDistributionsTransport[id * noDistributions() + eid];
}

/// Accessor for properties.
template <MInt nDim>
LbCellCollector<nDim>::BitsetType::reference LbCellCollector<nDim>::hasProperty(const MInt id, const LbCell p) {
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_PROPERTY_ACCESSOR(p);
  return m_properties[id][maia::lb::cell::p(p)];
}
/// Accessor for properties (const version).
template <MInt nDim>
MBool LbCellCollector<nDim>::hasProperty(const MInt id, const LbCell p) const {
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_PROPERTY_ACCESSOR(p);
  return m_properties[id][maia::lb::cell::p(p)];
}
/// Reset all properties.
template <MInt nDim>
void LbCellCollector<nDim>::resetProperties(const MInt id) {
  ENSURE_VALID_ID_ACCESSOR(id);
  m_properties[id].reset();
}

/// Accessor for properties.
template <MInt nDim>
BitsetType& LbCellCollector<nDim>::allProperties(const MInt id) {
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_properties[id];
}

/// Set use of thermal model
template <MInt nDim>
void LbCellCollector<nDim>::setThermal(const MBool isThermal_) {
  m_isThermal = isThermal_;
}

/// Set use of transport model
template <MInt nDim>
void LbCellCollector<nDim>::setTransport(const MBool useTransport_) {
  m_useTransport = useTransport_;
}

/// Update number of variables
template <MInt nDim>
void LbCellCollector<nDim>::setNoVariables() {
  if(isThermal() && !useTransport()) {
    // only thermal
    m_noVariables = 1 + nDim + 1;
  } else if(!isThermal() && useTransport()) {
    // only transport
    // TODO: So far both variables must be stored, but this will change
    // when the sysEq formulation is used in the LB
    m_noVariables = 1 + nDim + 2;
  } else if(isThermal() && useTransport()) {
    // coupled thermal and transport
    m_noVariables = 1 + nDim + 2;
  } else {
    // standard model
    m_noVariables = 1 + nDim;
  }
}

template <MInt nDim>
void LbCellCollector<nDim>::setSaveUOtherPhase(const MBool saveUOtherPhase_) {
  m_saveUOtherPhase = saveUOtherPhase_;
}

template <MInt nDim>
void LbCellCollector<nDim>::setSaveVolumeFraction(const MBool saveVolumeFraction_) {
  m_saveVolumeFraction = saveVolumeFraction_;
}

template <MInt nDim>
void LbCellCollector<nDim>::setSavePrevVars(const MBool savePrevVars_) {
  m_savePrevVars = savePrevVars_;
}

template <MInt nDim>
void LbCellCollector<nDim>::setSaveNuT(const MBool saveNuT_) {
  m_saveNuT = saveNuT_;
}

template <MInt nDim>
void LbCellCollector<nDim>::setSaveOldNu(const MBool saveOldNu_) {
  m_saveOldNu = saveOldNu_;
}

/// Sets the number of distributions
template <MInt nDim>
void LbCellCollector<nDim>::setNoDistributions(const MInt noDistributions_) {
  m_noDistributions = noDistributions_;
}

/// Erase range of nodes such that they contain no sensible values anymore.
template <MInt nDim>
void LbCellCollector<nDim>::invalidate(const MInt begin, const MInt end) {

  fill_invalid(m_nu, begin, end);
  if(saveOldNu()) {
    fill_invalid(m_oldNu, begin, end);
  }
  fill_invalid(m_kappa, begin, end);
  fill_invalid(m_diffusivity, begin, end);
  fill_invalid(m_spongeFactor, begin, end);
  fill_invalid(m_bndId, begin, end, 1, -1 /* TODO labels:LB Invalid<MInt>::value() */);
  fill_invalid(m_level, begin, end);
  fill_invalid(m_variables, begin, end, noVariables(), 0.);
  fill_invalid(m_oldVariables, begin, end, noVariables(), 0.);
  fill_invalid(m_distributions, begin, end, noDistributions());
  fill_invalid(m_oldDistributions, begin, end, noDistributions());
  fill_invalid(m_externalForces, begin, end, nDim, 0.0);
  if(isThermal()) {
    fill_invalid(m_distributionsThermal, begin, end, noDistributions());
    fill_invalid(m_oldDistributionsThermal, begin, end, noDistributions());
  }
  if(useTransport()) {
    fill_invalid(m_distributionsTransport, begin, end, noDistributions());
    fill_invalid(m_oldDistributionsTransport, begin, end, noDistributions());
  }
  if(savePrevVars()) {
    fill_invalid(m_previousDistribution, begin, end, noDistributions());
    fill_invalid(m_previousVariable, begin, end, noVariables());
  }
  if(saveNuT()) {
    fill_invalid(m_nuT, begin, end);
    if(saveOldNu()) {
      fill_invalid(m_oldNuT, begin, end);
    }
  }
  if(saveUOtherPhase()) {
    fill_invalid(m_uOtherPhase, begin, end, nDim);
  }
  if(saveVolumeFraction()) {
    fill_invalid(m_invVolumeFraction, begin, end);
  }

  // Properties
  fill_invalid(m_properties, begin, end);
}


/// Helper function for rawCopy(). Destination may refer to beginning or end of target range.
template <MInt nDim>
template <class Functor, class T>
void LbCellCollector<nDim>::rawCopyGeneric(Functor&& c, const T& source, const MInt begin, const MInt end,
                                           const MInt destination) {

  copyData(source.m_nu, m_nu, c, begin, end, destination);
  if(saveOldNu()) {
    copyData(source.m_oldNu, m_oldNu, c, begin, end, destination);
  }
  copyData(source.m_kappa, m_kappa, c, begin, end, destination);
  copyData(source.m_diffusivity, m_diffusivity, c, begin, end, destination);
  copyData(source.m_spongeFactor, m_spongeFactor, c, begin, end, destination);
  copyData(source.m_bndId, m_bndId, c, begin, end, destination);
  copyData(source.m_level, m_level, c, begin, end, destination);
  copyData(source.m_variables, m_variables, c, begin, end, destination, noVariables());
  copyData(source.m_oldVariables, m_oldVariables, c, begin, end, destination, noVariables());
  copyData(source.m_distributions, m_distributions, c, begin, end, destination, noDistributions());
  copyData(source.m_oldDistributions, m_oldDistributions, c, begin, end, destination, noDistributions());
  copyData(source.m_externalForces, m_externalForces, c, begin, end, destination, nDim);
  if(isThermal()) {
    copyData(source.m_distributionsThermal, m_distributionsThermal, c, begin, end, destination, noDistributions());
    copyData(source.m_oldDistributionsThermal, m_oldDistributionsThermal, c, begin, end, destination,
             noDistributions());
  }
  if(useTransport()) {
    copyData(source.m_distributionsTransport, m_distributionsTransport, c, begin, end, destination, noDistributions());
    copyData(source.m_oldDistributionsTransport, m_oldDistributionsTransport, c, begin, end, destination,
             noDistributions());
  }
  if(savePrevVars()) {
    copyData(source.m_previousDistribution, m_previousDistribution, c, begin, end, destination, noDistributions());
    copyData(source.m_previousVariable, m_previousVariable, c, begin, end, destination, noVariables());
  }
  if(saveNuT()) {
    copyData(source.m_nuT, m_nuT, c, begin, end, destination);
    if(saveOldNu()) {
      copyData(source.m_oldNuT, m_oldNuT, c, begin, end, destination);
    }
  }
  if(saveUOtherPhase()) {
    copyData(source.m_uOtherPhase, m_uOtherPhase, c, begin, end, destination, nDim);
  }
  if(saveVolumeFraction()) {
    copyData(source.m_invVolumeFraction, m_invVolumeFraction, c, begin, end, destination);
  }

  // Properties
  copyData(source.m_properties, m_properties, c, begin, end, destination);
}

} // namespace collector
} // namespace lb
} // namespace maia


// Undefine macros that should not be used outside this file
#undef LBCOLLECTOR_SANITY_CHECKS_ACCESSORS
#undef ENSURE_VALID_ID_ACCESSOR
#undef ENSURE_VALID_VARIABLE_ID_ACCESSOR
#undef ENSURE_VALID_DISTRIBUTION_ID_ACCESSOR
#undef ENSURE_VALID_PROPERTY_ACCESSOR

#endif // ifndef LBCOLLECTOR_H_
