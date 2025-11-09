// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef FVCELLCOLLECTOR_H_
#define FVCELLCOLLECTOR_H_

#include <algorithm>
#include <bitset>
#include <type_traits>
#include <vector>
#include "INCLUDE/maiaconstants.h"
#include "INCLUDE/maiamacro.h"
#include "INCLUDE/maiatypes.h"
#include "IO/context.h"
#include "MEMORY/container.h"
#include "UTIL/functions.h"
#include "compiler_config.h"
#include "fvcartesiancellproperties.h"
#include "property.h"

// The following macro enables the "Structure-of-Arrays" memory layout for multi-dimensional node
// variables. This might be beneficial for GPU computations. Default is "Array-of-Structures".
// Examples (for nodes nN with four children cM each)
// Array-of-Structures (AOS): n0c0, n0c1, n0c2, n0c3, n1c0, n1c1, n1c2, n1c3, n2c0, n2c1, ...
// Structure-of-Arrays (SOA): n0c0, n1c0, n2c0, n3c0, ..., n0c1, n1c1, n2c1, n3c1, ..., n0c2, ...
// #define FVCELLCOLLECTOR_SOA_MEMORY_LAYOUT

// Sanity-checking macros for accessors
#if defined(FVCELLCOLLECTOR_SANITY_CHECKS_ACCESSORS) || defined(MAIA_ASSERT_ACCESSORS)
#define ENSURE_VALID_ID_ACCESSOR(id)                                                                                   \
  do {                                                                                                                 \
    MAIA_CONTAINER_ENSURE_VALID_ID(id);                                                                                \
  } while(false)
#define ENSURE_VALID_VARIABLE_ID_ACCESSOR(id)                                                                          \
  do {                                                                                                                 \
    MAIA_CONTAINER_ENSURE(                                                                                             \
        id >= 0 && id < noCVariables(),                                                                                \
        "variable id = " + std::to_string(id) + " out-of-bounds [0, " + std::to_string(noCVariables()) + ")", AT_);    \
  } while(false)
#define ENSURE_VALID_PVARIABLE_ID_ACCESSOR(id)                                                                         \
  do {                                                                                                                 \
    MAIA_CONTAINER_ENSURE(id >= 0 && id < noPVariables(),                                                              \
                          "primitive variable id = " + std::to_string(id) + " out-of-bounds [0, "                      \
                              + std::to_string(noPVariables()) + ")",                                                  \
                          AT_);                                                                                        \
  } while(false)
#define ENSURE_VALID_FVARIABLE_ID_ACCESSOR(id)                                                                         \
  do {                                                                                                                 \
    MAIA_CONTAINER_ENSURE(id >= 0 && id < noFVariables(),                                                              \
                          "flux variable id = " + std::to_string(id) + " out-of-bounds [0, "                           \
                              + std::to_string(noFVariables()) + ")",                                                  \
                          AT_);                                                                                        \
  } while(false)
#define ENSURE_VALID_AVARIABLE_ID_ACCESSOR(id)                                                                         \
  do {                                                                                                                 \
    MAIA_CONTAINER_ENSURE(id >= 0 && id < noAVariables(),                                                              \
                          "additional variable id = " + std::to_string(id) + " out-of-bounds [0, "                     \
                              + std::to_string(noAVariables()) + ")",                                                  \
                          AT_);                                                                                        \
  } while(false)
#define ENSURE_VALID_DIRECTION_ID_ACCESSOR(id)                                                                         \
  do {                                                                                                                 \
    MAIA_CONTAINER_ENSURE(id >= 0 && id < nDim,                                                                        \
                          "direction id = " + std::to_string(id) + " out-of-bounds [0, " + std::to_string(nDim) + ")", \
                          AT_);                                                                                        \
  } while(false)
#define ENSURE_VALID_RECNGHBR_ID_ACCESSOR(id)                                                                          \
  do {                                                                                                                 \
    MAIA_CONTAINER_ENSURE(id >= 0 && id < noRecNghbrs(),                                                               \
                          "reconstruction neighbor id = " + std::to_string(id) + " out-of-bounds [0, "                 \
                              + std::to_string(noRecNghbrs()) + ")",                                                   \
                          AT_);                                                                                        \
  } while(false)
#define ENSURE_VALID_REACTION_ID_ACCESSOR(id)                                                                          \
  do {                                                                                                                 \
    MAIA_CONTAINER_ENSURE(id >= 0 && id < noReactionRates(),                                                           \
                          "reaction rate id = " + std::to_string(id) + " out-of-bounds [0, "                           \
                              + std::to_string(noReactionRates()) + ")",                                               \
                          AT_);                                                                                        \
  } while(false)
#define ENSURE_VALID_IMPLICIT_COEFFICIENT_ID_ACCESSOR(id)                                                              \
  do {                                                                                                                 \
    MAIA_CONTAINER_ENSURE(id >= 0 && id < noImplicitCoefficients(),                                                    \
                          "implict coefficient id = " + std::to_string(id) + " out-of-bounds [0, "                     \
                              + std::to_string(noImplicitCoefficients()) + ")",                                        \
                          AT_);                                                                                        \
  } while(false)
#define ENSURE_VALID_COORDINATE_DIR_ACCESSOR(dir)                                                                      \
  do {                                                                                                                 \
    MAIA_CONTAINER_ENSURE(dir >= 0 && dir < nDim,                                                                      \
                          "coordinate direction dir = " + std::to_string(dir) + " out-of-bounds [0, "                  \
                              + std::to_string(nDim) + ")",                                                            \
                          AT_);                                                                                        \
  } while(false)
#define ENSURE_VALID_PROPERTY_ACCESSOR(p)                                                                              \
  do {                                                                                                                 \
    MAIA_CONTAINER_ENSURE(p != FvCell::NumProperties, "Invalid property", AT_);                                        \
  } while(false)
#else
#define ENSURE_VALID_ID_ACCESSOR(id)                                                                                   \
  do {                                                                                                                 \
  } while(false)
#define ENSURE_VALID_VARIABLE_ID_ACCESSOR(id)                                                                          \
  do {                                                                                                                 \
  } while(false)
#define ENSURE_VALID_PVARIABLE_ID_ACCESSOR(id)                                                                         \
  do {                                                                                                                 \
  } while(false)
#define ENSURE_VALID_FVARIABLE_ID_ACCESSOR(id)                                                                         \
  do {                                                                                                                 \
  } while(false)
#define ENSURE_VALID_AVARIABLE_ID_ACCESSOR(id)                                                                         \
  do {                                                                                                                 \
  } while(false)
#define ENSURE_VALID_DIRECTION_ID_ACCESSOR(id)                                                                         \
  do {                                                                                                                 \
  } while(false)
#define ENSURE_VALID_RECNGHBR_ID_ACCESSOR(id)                                                                          \
  do {                                                                                                                 \
  } while(false)
#define ENSURE_VALID_REACTION_ID_ACCESSOR(id)                                                                          \
  do {                                                                                                                 \
  } while(false)
#define ENSURE_VALID_IMPLICIT_COEFFICIENT_ID_ACCESSOR(id)                                                              \
  do {                                                                                                                 \
  } while(false)
#define ENSURE_VALID_COORDINATE_DIR_ACCESSOR(dir)                                                                      \
  do {                                                                                                                 \
  } while(false)
#define ENSURE_VALID_PROPERTY_ACCESSOR(dir)                                                                            \
  do {                                                                                                                 \
  } while(false)
#endif


/// Namespace for auxiliary functions/classes
namespace maia {
namespace fv {
namespace collector {

/// Underlying bitset type for property storage
using BitsetType = maia::fv::cell::BitsetType;

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

// Invalid value for ... is ...
template <>
struct Invalid<M16X2bit<false>> {
  static constexpr M16X2bit<false>::type value() { return M16X2bit<false>::type{}; }
};


/// Class that represents FV cell collector.
template <MInt nDim>
class FvCellCollector : public maia::container::Container<FvCellCollector<nDim>, Invalid> {
  // Necessary for CRTP
  friend class maia::container::Container<FvCellCollector<nDim>, Invalid>;

  // Make base class functions known to use without this pointer
  using Base = maia::container::Container<FvCellCollector<nDim>, Invalid>;
  using Base::resetStorage;
  using Base::resizeStorage;
  template <class T>
  using Storage = typename Base::template Storage<T>;
  using BitsetType = maia::fv::cell::BitsetType;


 public:
  void setFvCollectorType(const MInt mode);
  void setFvTimeStepType(const MInt mode);
  // accessors to specify FvCell collector type
  constexpr MBool hasCellCenterMeanMolarWeight() const { return m_hasCellCenterMeanMolarWeight; }
  constexpr MBool hasCellCenterGamma() const { return m_hasCellCenterGamma; }
  constexpr MBool hasReactionRates() const { return m_hasReactionRates; }
  constexpr MBool hasReactionRatesBackup() const { return m_hasReactionRatesBackup; }
  constexpr MBool hasPsi() const { return m_hasPsi; }
  constexpr MBool isEEGas() const { return m_isEEGas; }
  constexpr MBool hasLocalTS() const { return m_hasLocalTS; }
  constexpr MBool hasDualTS() const { return m_hasDualTS; }

  // Types
  template <class T>
  using Invalid = typename maia::fv::collector::Invalid<T>;

  // Constructors
  /// Default c'tor does nothing
  constexpr FvCellCollector() = default;

  // Ensure that base class method is found when called from outside
  using Base::copyData;
  using Base::fill_invalid;
  using Base::reset;
  using Base::resize;

  // Accessors
  MFloat& oldVariable(const MInt id, const MInt varId);
  MFloat oldVariable(const MInt id, const MInt varId) const;

  MFloat& variable(const MInt id, const MInt dim);
  MFloat variable(const MInt id, const MInt dim) const;

  MFloat& pvariable(const MInt id, const MInt dim);
  MFloat pvariable(const MInt id, const MInt dim) const;

  MFloat& avariable(const MInt id, const MInt dim);
  MFloat avariable(const MInt id, const MInt dim) const;

  MFloat& rightHandSide(const MInt id, const MInt varId);
  MFloat rightHandSide(const MInt id, const MInt varId) const;

  MFloat& slope(const MInt id, const MInt dimVar, const MInt dimDir);
  MFloat slope(const MInt id, const MInt dimVar, const MInt dimDir) const;

  MInt& bndryCellId(const MInt id);
  MInt bndryCellId(const MInt id) const;

  MInt& noRcnstrctnNghbrIds(const MInt id);
  MInt noRcnstrctnNghbrIds(const MInt id) const;

  MInt& rcnstrctnNghbrId(const MInt id, const MInt dimRecNghbr);
  MInt rcnstrctnNghbrId(const MInt id, const MInt dimRecNghbr) const;

  MInt& reconstructionData(const MInt id);
  MInt reconstructionData(const MInt id) const;

  void nghbrInterface(const MInt id, const MInt dir, const MInt state);
  MInt nghbrInterface(const MInt id, const MInt dir) const;

  MFloat& reactionRate(const MInt id, const MInt dimReaction);
  MFloat reactionRate(const MInt id, const MInt dimReaction) const;

  MFloat& reactionRateBackup(const MInt id, const MInt dimReaction);
  MFloat reactionRateBackup(const MInt id, const MInt dimReaction) const;

  MFloat& psi(const MInt id);
  MFloat psi(const MInt id) const;

  MFloat& implicitCoefficient(const MInt id, const MInt dimCoefficient);
  MFloat implicitCoefficient(const MInt id, const MInt dimCoefficient) const;

  MFloat& speciesReactionRate(const MInt id, const MInt speciesIndex);
  MFloat speciesReactionRate(const MInt id, const MInt speciesIndex) const;

  MFloat& cellCenterMeanMolarWeight(const MInt id);
  MFloat cellCenterMeanMolarWeight(const MInt id) const;

  MFloat& cellCenterGamma(const MInt id);
  MFloat cellCenterGamma(const MInt id) const;

  MFloat& dt1Variable(const MInt id, const MInt dim);
  MFloat dt1Variable(const MInt id, const MInt dim) const;

  MFloat& dt2Variable(const MInt id, const MInt dim);
  MFloat dt2Variable(const MInt id, const MInt dim) const;

  MFloat& localTimeStep(const MInt id);
  MFloat localTimeStep(const MInt id) const;

  MFloat& fluidFraction(const MInt id);
  const MFloat& fluidFraction(const MInt id) const;

  MFloat& spongeFactor(const MInt id);
  MFloat spongeFactor(const MInt id) const;

  MFloat& spongeFactorStart(const MInt id);
  MFloat spongeFactorStart(const MInt id) const;

  MInt& spongeBndryId(const MInt id, const MInt dimDir);
  MInt spongeBndryId(const MInt id, const MInt dimDir) const;

  MFloat& coordinate(const MInt id, const MInt dim);
  MFloat coordinate(const MInt id, const MInt dim) const;

  MInt& level(const MInt cellId);
  MInt level(const MInt cellId) const;

  MFloat& cellVolume(const MInt cellId);
  MFloat cellVolume(const MInt cellId) const;

  MFloat& FcellVolume(const MInt cellId);
  MFloat FcellVolume(const MInt cellId) const;

  // Multilevel-related accessors
  MFloat& tau(const MInt id, const MInt varId);
  MFloat tau(const MInt id, const MInt varId) const;

  MFloat& restrictedRHS(const MInt id, const MInt varId);
  MFloat restrictedRHS(const MInt id, const MInt varId) const;

  MFloat& restrictedVar(const MInt id, const MInt varId);
  MFloat restrictedVar(const MInt id, const MInt varId) const;

  MFloat& storedSlope(const MInt id, const MInt dimVar, const MInt dimDir);
  MFloat storedSlope(const MInt id, const MInt dimVar, const MInt dimDir) const;

  // Property-related accessors
  BitsetType::reference hasProperty(const MInt id, const FvCell p);
  MBool hasProperty(const MInt id, const FvCell p) const;
  void resetProperties(const MInt id);
  BitsetType& properties(const MInt id);

  // Allow setting number of species and rans variables
  void setNoCVariables(const MInt noCVariables_, const MInt noSpecies_);

  // Allow setting number of PVs
  void setNoPVariables(const MInt noPVariables_);

  // Allow setting number of FVs
  void setNoFVariables(const MInt noFVariables_);

  // Allow setting number of AVs
  void setNoAVariables(const MInt noAVariables_);

  /// Return number of species
  constexpr MInt noSpecies() const { return m_noSpecies; }

  /// Allow activation of multilevel calculations
  MInt isMultilevel(const MBool isMultilevel_);

  /// Return whether multilevel is active or not
  constexpr MInt isMultilevel() const { return m_isMultilevel; }

  /// Return number of conservative variables
  constexpr MInt noCVariables() const { return m_noCVariables; }

  /// Return number of primitive variables
  constexpr MInt noPVariables() const { return m_noPVariables; }

  /// Return number of flux variables
  constexpr MInt noFVariables() const { return m_noFVariables; }

  /// Return number of additional variables
  constexpr MInt noAVariables() const { return m_noAVariables; }

  /// Return max number of reconstruction nghbrs
  constexpr MInt noRecNghbrs() const { return m_noRecNghbrs; }

  /// Return max number of reaction rates
  constexpr MInt noReactionRates() const { return m_noReactionRates; }

  /// Return max number of implicit coefficients
  constexpr MInt noImplicitCoefficients() const { return m_noImplicitCoefficients; }

  /// Return number of properties defined for each node
  static constexpr MInt noProperties() { return maia::fv::cell::p(FvCell::NumProperties); }

 private:
  // Bools for fv-collector type
  MBool m_hasCellCenterMeanMolarWeight = false;
  MBool m_hasReactionRates = false;
  MBool m_hasReactionRatesBackup = false;
  MBool m_hasPsi = false;
  MBool m_isEEGas = false;
  MBool m_hasCellCenterGamma = false;
  MBool m_hasLocalTS = false;
  MBool m_hasDualTS = false;


  // Methods required by base class for CRTP
  void reset();
  void resize() override;
  void invalidate(const MInt begin, const MInt end);
  template <class Functor, class T>
  void rawCopyGeneric(Functor&& c, const T& source, const MInt begin, const MInt end, const MInt destination);

  /// Number of variables
  MInt m_noCVariables = 2 + nDim;
  MInt m_noPVariables = 2 + nDim;
  MInt m_noFVariables = 2 + nDim;
  MInt m_noAVariables = 0;

  /// Number of species
  MInt m_noSpecies = 0;

  /// Is multilevel computation activated
  MBool m_isMultilevel = false;

  /// Max number of reconstruction neighbors
  MInt m_noRecNghbrs = IPOW3[nDim] + nDim;

  /// Max number of reaction rates
  MInt m_noReactionRates = 1;

  /// Max number of implicit coefficients
  MInt m_noImplicitCoefficients = 2 * nDim;

  // Data containers
  Storage<MFloat> m_oldVariables{};
  Storage<MFloat> m_variables{};
  Storage<MFloat> m_avariables{};
  Storage<MFloat> m_rightHandSide{};
  Storage<MFloat> m_slopes{};
  Storage<MFloat> m_pvariables{};
  Storage<MInt> m_bndryCellIds{};
  Storage<MInt> m_noRcnstrctnNghbrIds{};
  Storage<MInt> m_rcnstrctnNghbrIds{};
  Storage<MInt> m_reconstructionData{};
  Storage<M16X2bit<false>> m_nghbrInterface{};
  Storage<MFloat> m_reactionRates{};
  Storage<MFloat> m_reactionRatesBackup{};
  Storage<MFloat> m_psi{};
  Storage<MFloat> m_implicitCoefficients{};
  Storage<MFloat> m_speciesReactionRates{};
  Storage<MFloat> m_cellCenterMeanMolarWeight{};
  Storage<MFloat> m_cellCenterGamma{};
  Storage<MFloat> m_dt1Variables{};
  Storage<MFloat> m_dt2Variables{};
  /// Storage of cell-local time-steps for dual time stepping:
  Storage<MFloat> m_localTimeStep_{};
  Storage<MFloat> m_spongeFactor{};
  Storage<MFloat> m_spongeFactorStart{};
  Storage<MInt> m_spongeBndryIds{};
  Storage<MFloat> m_coordinates{};
  Storage<MInt> m_levels{};
  Storage<BitsetType> m_properties{};
  Storage<MFloat> m_cellVolumes{};
  Storage<MFloat> m_FcellVolumes{};

  // Multilevel-related data
  Storage<MFloat> m_tau{};
  Storage<MFloat> m_restrictedRHS{};
  Storage<MFloat> m_restrictedVars{};
  Storage<MFloat> m_storedSlopes{};
};


/// Reset tree, re-create data structures with given capacity, and set size to zero.
template <MInt nDim>
void FvCellCollector<nDim>::reset() {
  resetStorage(noCVariables(), m_oldVariables);
  resetStorage(noCVariables(), m_variables);
  resetStorage(noFVariables(), m_rightHandSide);
  resetStorage(noPVariables() * nDim, m_slopes);
  resetStorage(noPVariables(), m_pvariables);
  resetStorage(1, m_bndryCellIds);
  resetStorage(1, m_noRcnstrctnNghbrIds);
  resetStorage(noRecNghbrs(), m_rcnstrctnNghbrIds);
  resetStorage(1, m_reconstructionData);
  resetStorage(1, m_nghbrInterface);
  if(hasReactionRates()) resetStorage(noReactionRates(), m_reactionRates);
  if(hasReactionRates()) resetStorage(noSpecies(), m_speciesReactionRates);
  if(hasReactionRatesBackup()) resetStorage(noReactionRates(), m_reactionRatesBackup);
  if(hasPsi()) resetStorage(1, m_psi);
  if(hasCellCenterMeanMolarWeight()) resetStorage(1, m_cellCenterMeanMolarWeight);
  if(hasCellCenterGamma()) resetStorage(1, m_cellCenterGamma);
  if(isEEGas()) resetStorage(noImplicitCoefficients(), m_implicitCoefficients);
  resetStorage(noAVariables(), m_avariables);
  if(hasDualTS()) {
    resetStorage(noCVariables(), m_dt1Variables);
    resetStorage(noCVariables(), m_dt2Variables);
    resetStorage(1, m_localTimeStep_);
  } else if(hasLocalTS()) {
    resetStorage(1, m_localTimeStep_);
  }
  resetStorage(1, m_spongeFactor);
  resetStorage(1, m_spongeFactorStart);
  resetStorage(nDim, m_spongeBndryIds);
  resetStorage(nDim, m_coordinates);
  resetStorage(1, m_levels);
  resetStorage(1, m_properties);
  resetStorage(1, m_cellVolumes);
  resetStorage(1, m_FcellVolumes);

  // Initialize multilevel storage only if using more than one grid level
  if(isMultilevel()) {
    resetStorage(noCVariables(), m_tau);
    resetStorage(noCVariables(), m_restrictedRHS);
    resetStorage(noCVariables(), m_restrictedVars);
    resetStorage(noCVariables() * nDim, m_storedSlopes);
  }
}

/// TODO doc
template <MInt nDim>
void FvCellCollector<nDim>::resize() {
  resizeStorage(noCVariables(), m_oldVariables);
  resizeStorage(noCVariables(), m_variables);
  resizeStorage(noFVariables(), m_rightHandSide);
  resizeStorage(noPVariables() * nDim, m_slopes);
  resizeStorage(noPVariables(), m_pvariables);
  resizeStorage(1, m_bndryCellIds);
  resizeStorage(1, m_noRcnstrctnNghbrIds);
  resizeStorage(noRecNghbrs(), m_rcnstrctnNghbrIds);
  resizeStorage(1, m_reconstructionData);
  if(hasReactionRates()) resizeStorage(noReactionRates(), m_reactionRates);
  if(hasReactionRates()) resizeStorage(noSpecies(), m_speciesReactionRates);
  if(hasReactionRatesBackup()) resizeStorage(noReactionRates(), m_reactionRatesBackup);
  if(hasPsi()) resizeStorage(1, m_psi);
  if(hasCellCenterMeanMolarWeight()) resizeStorage(1, m_cellCenterMeanMolarWeight);
  if(hasCellCenterGamma()) resizeStorage(1, m_cellCenterGamma);
  if(isEEGas()) resizeStorage(noImplicitCoefficients(), m_implicitCoefficients);
  resizeStorage(noAVariables(), m_avariables);

  if(hasDualTS()) {
    resizeStorage(noCVariables(), m_dt1Variables);
    resizeStorage(noCVariables(), m_dt2Variables);
    resizeStorage(1, m_localTimeStep_);
  } else if(hasLocalTS()) {
    resizeStorage(1, m_localTimeStep_);
  }

  resizeStorage(1, m_spongeFactor);
  resizeStorage(1, m_spongeFactorStart);
  resizeStorage(nDim, m_spongeBndryIds);
  resizeStorage(nDim, m_coordinates);
  resizeStorage(1, m_levels);
  resizeStorage(1, m_properties);
  resizeStorage(1, m_cellVolumes);
  resizeStorage(1, m_FcellVolumes);

  // Initialize multilevel storage only if using more than one grid level
  if(isMultilevel()) {
    resizeStorage(noCVariables(), m_tau);
    resizeStorage(noCVariables(), m_restrictedRHS);
    resizeStorage(noCVariables(), m_restrictedVars);
    resizeStorage(noCVariables() * nDim, m_storedSlopes);
  }
}


/// Set the fv-collector type
template <MInt nDim>
void FvCellCollector<nDim>::setFvCollectorType(const MInt mode) {
  ASSERT(mode >= 0, "");
  switch(mode) {
    case 0: {
      m_hasCellCenterMeanMolarWeight = false;
      m_hasCellCenterGamma = false;
      m_hasReactionRates = false;
      m_hasReactionRatesBackup = false;
      m_hasPsi = false;
      m_isEEGas = false;
      break;
    }
    case 1: // detailedChemistry:
    {
      m_hasCellCenterMeanMolarWeight = true;
      m_hasCellCenterGamma = true;
      m_hasReactionRates = true;
      m_hasReactionRatesBackup = false;
      m_hasPsi = false;
      m_isEEGas = false;
      break;
    }
    case 2: { // combustion without detailedChemisty
      m_hasCellCenterMeanMolarWeight = false;
      m_hasCellCenterGamma = false;
      m_hasReactionRates = true;
      m_hasReactionRatesBackup = true;
      m_hasPsi = true;
      m_isEEGas = false;
      break;
    }
    case 3: { // EEGas
      m_hasCellCenterMeanMolarWeight = false;
      m_hasCellCenterGamma = false;
      m_hasReactionRates = false;
      m_hasReactionRatesBackup = false;
      m_hasPsi = false;
      m_isEEGas = true;
      break;
    }
    default: {
      mTerm(1, AT_, "Unknown fvCollector Type!");
    }
  }
}
template <MInt nDim>
void FvCellCollector<nDim>::setFvTimeStepType(const MInt mode) {
  ASSERT(mode >= 0, "");
  switch(mode) {
    case 0: {
      // default: global single time-step
      break;
    }
    case 1: {
      // local time-step with standart RK
      m_hasLocalTS = true;
      break;
    }
    case 2: {
      // dual time-stepping with local TS
      m_hasDualTS = true;
      break;
    }
    default: {
      mTerm(1, AT_, "Unknown fvTimeStepType!");
    }
  }
}

/// Accessor for oldVariables.
template <MInt nDim>
MFloat& FvCellCollector<nDim>::oldVariable(const MInt id, const MInt varId) {
// Prevent accidental compilation without support for SoA layout
#ifdef FVCELLCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_VARIABLE_ID_ACCESSOR(varId);
  return m_oldVariables[id * noCVariables() + varId];
}
/// Accessor for oldVariables (const version).
template <MInt nDim>
MFloat FvCellCollector<nDim>::oldVariable(const MInt id, const MInt varId) const {
// Prevent accidental compilation without support for SoA layout
#ifdef FVCELLCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_VARIABLE_ID_ACCESSOR(varId);
  return m_oldVariables[id * noCVariables() + varId];
}

/// Accessor for variables.
template <MInt nDim>
MFloat& FvCellCollector<nDim>::variable(const MInt id, const MInt varId) {
// Prevent accidental compilation without support for SoA layout
#ifdef FVCELLCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_VARIABLE_ID_ACCESSOR(varId);
  return m_variables[id * noCVariables() + varId];
}
/// Accessor for variables (const version).
template <MInt nDim>
MFloat FvCellCollector<nDim>::variable(const MInt id, const MInt varId) const {
// Prevent accidental compilation without support for SoA layout
#ifdef FVCELLCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_VARIABLE_ID_ACCESSOR(varId);
  return m_variables[id * noCVariables() + varId];
}

/// Accessor for pvariable.
template <MInt nDim>
MFloat& FvCellCollector<nDim>::pvariable(const MInt id, const MInt varId) {
// Prevent accidental compilation without support for SoA layout
#ifdef FVCELLCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_PVARIABLE_ID_ACCESSOR(varId);
  return m_pvariables[id * noPVariables() + varId];
}
/// Accessor for pvariable (const version).
template <MInt nDim>
MFloat FvCellCollector<nDim>::pvariable(const MInt id, const MInt varId) const {
// Prevent accidental compilation without support for SoA layout
#ifdef FVCELLCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_PVARIABLE_ID_ACCESSOR(varId);
  return m_pvariables[id * noPVariables() + varId];
}

/// Accessor for additional variables.
template <MInt nDim>
MFloat& FvCellCollector<nDim>::avariable(const MInt id, const MInt varId) {
// Prevent accidental compilation without support for SoA layout
#ifdef FVCELLCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_AVARIABLE_ID_ACCESSOR(varId);
  return m_avariables[id * noAVariables() + varId];
}
/// Accessor for additional variables (const version).
template <MInt nDim>
MFloat FvCellCollector<nDim>::avariable(const MInt id, const MInt varId) const {
// Prevent accidental compilation without support for SoA layout
#ifdef FVCELLCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_AVARIABLE_ID_ACCESSOR(varId);
  return m_avariables[id * noAVariables() + varId];
}

/// Accessor for right hand side.
template <MInt nDim>
MFloat& FvCellCollector<nDim>::rightHandSide(const MInt id, const MInt varId) {
// Prevent accidental compilation without support for SoA layout
#ifdef FVCELLCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_FVARIABLE_ID_ACCESSOR(varId);
  return m_rightHandSide[id * noFVariables() + varId];
}
/// Accessor for right hand side (const version).
template <MInt nDim>
MFloat FvCellCollector<nDim>::rightHandSide(const MInt id, const MInt varId) const {
// Prevent accidental compilation without support for SoA layout
#ifdef FVCELLCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_FVARIABLE_ID_ACCESSOR(varId);
  return m_rightHandSide[id * noFVariables() + varId];
}

/// Accessor for slopes.
template <MInt nDim>
MFloat& FvCellCollector<nDim>::slope(const MInt id, const MInt dimVar, const MInt dimDir) {
// Prevent accidental compilation without support for SoA layout
#ifdef FVCELLCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_PVARIABLE_ID_ACCESSOR(dimVar);
  ENSURE_VALID_DIRECTION_ID_ACCESSOR(dimDir);
  return m_slopes[id * noPVariables() * nDim + dimVar * nDim + dimDir];
}
/// Accessor for slopes (const version).
template <MInt nDim>
MFloat FvCellCollector<nDim>::slope(const MInt id, const MInt dimVar, const MInt dimDir) const {
// Prevent accidental compilation without support for SoA layout
#ifdef FVCELLCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_PVARIABLE_ID_ACCESSOR(dimVar);
  ENSURE_VALID_DIRECTION_ID_ACCESSOR(dimDir);
  return m_slopes[id * noPVariables() * nDim + dimVar * nDim + dimDir];
}

/// Accessor for bndryCellIds.
template <MInt nDim>
MInt& FvCellCollector<nDim>::bndryCellId(const MInt id) {
// Prevent accidental compilation without support for SoA layout
#ifdef FVCELLCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_bndryCellIds[id];
}
/// Accessor for bndryCellIds (const version).
template <MInt nDim>
MInt FvCellCollector<nDim>::bndryCellId(const MInt id) const {
// Prevent accidental compilation without support for SoA layout
#ifdef FVCELLCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_bndryCellIds[id];
}

/// Accessor for noRcnstrctnNghbrIds.
template <MInt nDim>
MInt& FvCellCollector<nDim>::noRcnstrctnNghbrIds(const MInt id) {
// Prevent accidental compilation without support for SoA layout
#ifdef FVCELLCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_noRcnstrctnNghbrIds[id];
}
/// Accessor for noRcnstrctnNghbrIds (const version).
template <MInt nDim>
MInt FvCellCollector<nDim>::noRcnstrctnNghbrIds(const MInt id) const {
// Prevent accidental compilation without support for SoA layout
#ifdef FVCELLCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_noRcnstrctnNghbrIds[id];
}

/// Accessor for rcnstrctnNghbrId.
template <MInt nDim>
MInt& FvCellCollector<nDim>::rcnstrctnNghbrId(const MInt id, const MInt dimRecNghbr) {
// Prevent accidental compilation without support for SoA layout
#ifdef FVCELLCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_RECNGHBR_ID_ACCESSOR(dimRecNghbr);
  return m_rcnstrctnNghbrIds[id * noRecNghbrs() + dimRecNghbr];
}
/// Accessor for rcnstrctnNghbrId (const version).
template <MInt nDim>
MInt FvCellCollector<nDim>::rcnstrctnNghbrId(const MInt id, const MInt dimRecNghbr) const {
// Prevent accidental compilation without support for SoA layout
#ifdef FVCELLCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_RECNGHBR_ID_ACCESSOR(dimRecNghbr);
  return m_rcnstrctnNghbrIds[id * noRecNghbrs() + dimRecNghbr];
}

/// Accessor for reconstructionData.
template <MInt nDim>
MInt& FvCellCollector<nDim>::reconstructionData(const MInt id) {
// Prevent accidental compilation without support for SoA layout
#ifdef FVCELLCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_reconstructionData[id];
}
/// Accessor for reconstructionData (const version).
template <MInt nDim>
MInt FvCellCollector<nDim>::reconstructionData(const MInt id) const {
// Prevent accidental compilation without support for SoA layout
#ifdef FVCELLCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_reconstructionData[id];
}

/// Accessor for nghbrInterface.
template <MInt nDim>
void FvCellCollector<nDim>::nghbrInterface(const MInt id, const MInt dir, const MInt state) {
// Prevent accidental compilation without support for SoA layout
#ifdef FVCELLCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_DIRECTION_ID_ACCESSOR(dir / 2);
  m_nghbrInterface[id].set(dir, state);
}
/// Accessor for rcnstrctnNghbrId (const version).
template <MInt nDim>
MInt FvCellCollector<nDim>::nghbrInterface(const MInt id, const MInt dir) const {
// Prevent accidental compilation without support for SoA layout
#ifdef FVCELLCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_DIRECTION_ID_ACCESSOR(dir / 2);
  return m_nghbrInterface[id].get(dir);
}

/// Accessor for reaction rates.
template <MInt nDim>
MFloat& FvCellCollector<nDim>::reactionRate(const MInt id, const MInt dimReaction) {
// Prevent accidental compilation without support for SoA layout
#ifdef FVCELLCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_REACTION_ID_ACCESSOR(dimReaction);
  return m_reactionRates[id * noReactionRates() + dimReaction];
}
/// Accessor for reaction rates (const version).
template <MInt nDim>
MFloat FvCellCollector<nDim>::reactionRate(const MInt id, const MInt dimReaction) const {
// Prevent accidental compilation without support for SoA layout
#ifdef FVCELLCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_REACTION_ID_ACCESSOR(dimReaction);
  return m_reactionRates[id * noReactionRates() + dimReaction];
}

/// Accessor for reaction rates backup.
template <MInt nDim>
MFloat& FvCellCollector<nDim>::reactionRateBackup(const MInt id, const MInt dimReaction) {
// Prevent accidental compilation without support for SoA layout
#ifdef FVCELLCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_REACTION_ID_ACCESSOR(dimReaction);
  return m_reactionRatesBackup[id * noReactionRates() + dimReaction];
}
/// Accessor for reaction rates backup (const version).
template <MInt nDim>
MFloat FvCellCollector<nDim>::reactionRateBackup(const MInt id, const MInt dimReaction) const {
// Prevent accidental compilation without support for SoA layout
#ifdef FVCELLCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_REACTION_ID_ACCESSOR(dimReaction);
  return m_reactionRatesBackup[id * noReactionRates() + dimReaction];
}

/// Accessor for psi.
template <MInt nDim>
MFloat& FvCellCollector<nDim>::psi(const MInt id) {
// Prevent accidental compilation without support for SoA layout
#ifdef FVCELLCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_psi[id];
}
/// Accessor for psi(const version).
template <MInt nDim>
MFloat FvCellCollector<nDim>::psi(const MInt id) const {
// Prevent accidental compilation without support for SoA layout
#ifdef FVCELLCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_psi[id];
}

/// Accessor for implicit coefficient.
template <MInt nDim>
MFloat& FvCellCollector<nDim>::implicitCoefficient(const MInt id, const MInt dimCoefficient) {
// Prevent accidental compilation without support for SoA layout
#ifdef FVCELLCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_IMPLICIT_COEFFICIENT_ID_ACCESSOR(dimCoefficient);
  return m_implicitCoefficients[id * noImplicitCoefficients() + dimCoefficient];
}
/// Accessor for  implicit coefficient (const version).
template <MInt nDim>
MFloat FvCellCollector<nDim>::implicitCoefficient(const MInt id, const MInt dimCoefficient) const {
// Prevent accidental compilation without support for SoA layout
#ifdef FVCELLCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_IMPLICIT_COEFFICIENT_ID_ACCESSOR(dimCoefficient);
  return m_implicitCoefficients[id * noImplicitCoefficients() + dimCoefficient];
}

/// Accessor for species reaction rates.
template <MInt nDim>
MFloat& FvCellCollector<nDim>::speciesReactionRate(const MInt id, const MInt speciesIndex) {
// Prevent accidental compilation without support for SoA layout
#ifdef FVCELLCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_REACTION_ID_ACCESSOR(speciesIndex);
  return m_speciesReactionRates[id * noSpecies() + speciesIndex];
}

/// Accessor for species reaction rates (const version).
template <MInt nDim>
MFloat FvCellCollector<nDim>::speciesReactionRate(const MInt id, const MInt speciesIndex) const {
// Prevent accidental compilation without support for SoA layout
#ifdef FVCELLCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_REACTION_ID_ACCESSOR(speciesIndex);
  return m_speciesReactionRates[id * noSpecies() + speciesIndex];
}

/// Accessor for species reaction rates.
template <MInt nDim>
MFloat& FvCellCollector<nDim>::cellCenterMeanMolarWeight(const MInt id) {
// Prevent accidental compilation without support for SoA layout
#ifdef FVCELLCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_cellCenterMeanMolarWeight[id];
}

/// Accessor for species reaction rates (const version).
template <MInt nDim>
MFloat FvCellCollector<nDim>::cellCenterMeanMolarWeight(const MInt id) const {
// Prevent accidental compilation without support for SoA layout
#ifdef FVCELLCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_cellCenterMeanMolarWeight[id];
}

/// Accessor for species reaction rates.
template <MInt nDim>
MFloat& FvCellCollector<nDim>::cellCenterGamma(const MInt id) {
// Prevent accidental compilation without support for SoA layout
#ifdef FVCELLCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_cellCenterGamma[id];
}

/// Accessor for species reaction rates (const version).
template <MInt nDim>
MFloat FvCellCollector<nDim>::cellCenterGamma(const MInt id) const {
// Prevent accidental compilation without support for SoA layout
#ifdef FVCELLCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_cellCenterGamma[id];
}

/// Accessor for dt1Variable
template <MInt nDim>
MFloat& FvCellCollector<nDim>::dt1Variable(const MInt id, const MInt varId) {
// Prevent accidental compilation without support for SoA layout
#ifdef FVCELLCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_VARIABLE_ID_ACCESSOR(varId);
  return m_dt1Variables[id * noCVariables() + varId];
}
/// Accessor for dt1Variable (const version).
template <MInt nDim>
MFloat FvCellCollector<nDim>::dt1Variable(const MInt id, const MInt varId) const {
// Prevent accidental compilation without support for SoA layout
#ifdef FVCELLCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_VARIABLE_ID_ACCESSOR(varId);
  return m_dt1Variables[id * noCVariables() + varId];
}

/// Accessor for dt2Variable
template <MInt nDim>
MFloat& FvCellCollector<nDim>::dt2Variable(const MInt id, const MInt varId) {
// Prevent accidental compilation without support for SoA layout
#ifdef FVCELLCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_VARIABLE_ID_ACCESSOR(varId);
  return m_dt2Variables[id * noCVariables() + varId];
}
/// Accessor for dt2Variable (const version).
template <MInt nDim>
MFloat FvCellCollector<nDim>::dt2Variable(const MInt id, const MInt varId) const {
// Prevent accidental compilation without support for SoA layout
#ifdef FVCELLCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_VARIABLE_ID_ACCESSOR(varId);
  return m_dt2Variables[id * noCVariables() + varId];
}

/// Accessor for local timestep (const version).
template <MInt nDim>
MFloat FvCellCollector<nDim>::localTimeStep(const MInt id) const {
#ifdef FVCELLCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_localTimeStep_[id];
}

/// Accessor for local timestep.
template <MInt nDim>
MFloat& FvCellCollector<nDim>::localTimeStep(const MInt id) {
#ifdef FVCELLCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_localTimeStep_[id];
}

/// Accessor for spongeFactor.
template <MInt nDim>
MFloat& FvCellCollector<nDim>::spongeFactor(const MInt id) {
// Prevent accidental compilation without support for SoA layout
#ifdef FVCELLCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_spongeFactor[id];
}
/// Accessor for spongeFactor (const version).
template <MInt nDim>
MFloat FvCellCollector<nDim>::spongeFactor(const MInt id) const {
// Prevent accidental compilation without support for SoA layout
#ifdef FVCELLCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_spongeFactor[id];
}

/// Accessor for spongeFactorStart.
template <MInt nDim>
MFloat& FvCellCollector<nDim>::spongeFactorStart(const MInt id) {
// Prevent accidental compilation without support for SoA layout
#ifdef FVCELLCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_spongeFactorStart[id];
}
/// Accessor for spongeFactorStart (const version).
template <MInt nDim>
MFloat FvCellCollector<nDim>::spongeFactorStart(const MInt id) const {
// Prevent accidental compilation without support for SoA layout
#ifdef FVCELLCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_spongeFactorStart[id];
}

/// Accessor for spongeBndryId.
template <MInt nDim>
MInt& FvCellCollector<nDim>::spongeBndryId(const MInt id, const MInt dimDir) {
// Prevent accidental compilation without support for SoA layout
#ifdef FVCELLCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_DIRECTION_ID_ACCESSOR(dimDir);
  return m_spongeBndryIds[id * nDim + dimDir];
}
/// Accessor for spongeBndryId (const version).
template <MInt nDim>
MInt FvCellCollector<nDim>::spongeBndryId(const MInt id, const MInt dimDir) const {
// Prevent accidental compilation without support for SoA layout
#ifdef FVCELLCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_DIRECTION_ID_ACCESSOR(dimDir);
  return m_spongeBndryIds[id * nDim + dimDir];
}

/// Accessor for coordinates.
template <MInt nDim>
MFloat& FvCellCollector<nDim>::coordinate(const MInt id, const MInt dir) {
// Prevent accidental compilation without support for SoA layout
#ifdef FVCELLCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_COORDINATE_DIR_ACCESSOR(dir);
  return m_coordinates[id * nDim + dir];
}
/// Accessor for coordinates (const version).
template <MInt nDim>
MFloat FvCellCollector<nDim>::coordinate(const MInt id, const MInt dir) const {
// Prevent accidental compilation without support for SoA layout
#ifdef FVCELLCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_COORDINATE_DIR_ACCESSOR(dir);
  return m_coordinates[id * nDim + dir];
}

/// Accessor for level.
template <MInt nDim>
MInt& FvCellCollector<nDim>::level(const MInt id) {
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_levels[id];
}
/// Accessor for level (const version).
template <MInt nDim>
MInt FvCellCollector<nDim>::level(const MInt id) const {
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_levels[id];
}

/// Accessor for cell volume.
template <MInt nDim>
MFloat& FvCellCollector<nDim>::cellVolume(const MInt id) {
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_cellVolumes[id];
}
/// Accessor for cell volume (const version).
template <MInt nDim>
MFloat FvCellCollector<nDim>::cellVolume(const MInt id) const {
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_cellVolumes[id];
}

/// Accessor for inverse cell volume.
template <MInt nDim>
MFloat& FvCellCollector<nDim>::FcellVolume(const MInt id) {
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_FcellVolumes[id];
}
/// Accessor for inverse cell volume (const version).
template <MInt nDim>
MFloat FvCellCollector<nDim>::FcellVolume(const MInt id) const {
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_FcellVolumes[id];
}

/// Accessor for properties.
template <MInt nDim>
FvCellCollector<nDim>::BitsetType::reference FvCellCollector<nDim>::hasProperty(const MInt id, const FvCell p) {
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_PROPERTY_ACCESSOR(p);
  return m_properties[id][maia::fv::cell::p(p)];
}
/// Accessor for properties (const version).
template <MInt nDim>
MBool FvCellCollector<nDim>::hasProperty(const MInt id, const FvCell p) const {
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_PROPERTY_ACCESSOR(p);
  return m_properties[id][maia::fv::cell::p(p)];
}
/// Reset all properties.
template <MInt nDim>
void FvCellCollector<nDim>::resetProperties(const MInt id) {
  ENSURE_VALID_ID_ACCESSOR(id);
  m_properties[id].reset();
}
/// Accessor for properties.
template <MInt nDim>
BitsetType& FvCellCollector<nDim>::properties(const MInt id) {
  ENSURE_VALID_ID_ACCESSOR(id);
  return m_properties[id];
}

/// Accessor for coarse cell correction (tau).
template <MInt nDim>
MFloat& FvCellCollector<nDim>::tau(const MInt id, const MInt varId) {
// Prevent accidental compilation without support for SoA layout
#ifdef FVCELLCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_VARIABLE_ID_ACCESSOR(varId);
  return m_tau[id * noCVariables() + varId];
}
/// Accessor for coarse cell correction (tau) (const version).
template <MInt nDim>
MFloat FvCellCollector<nDim>::tau(const MInt id, const MInt varId) const {
// Prevent accidental compilation without support for SoA layout
#ifdef FVCELLCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_VARIABLE_ID_ACCESSOR(varId);
  return m_tau[id * noCVariables() + varId];
}

/// Accessor for restricted RHS for multilevel computation.
template <MInt nDim>
MFloat& FvCellCollector<nDim>::restrictedRHS(const MInt id, const MInt varId) {
// Prevent accidental compilation without support for SoA layout
#ifdef FVCELLCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_VARIABLE_ID_ACCESSOR(varId);
  return m_restrictedRHS[id * noCVariables() + varId];
}
/// Accessor for restricted RHS for multilevel computation (const version).
template <MInt nDim>
MFloat FvCellCollector<nDim>::restrictedRHS(const MInt id, const MInt varId) const {
// Prevent accidental compilation without support for SoA layout
#ifdef FVCELLCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_VARIABLE_ID_ACCESSOR(varId);
  return m_restrictedRHS[id * noCVariables() + varId];
}

/// Accessor for variables after restriction during multigrid computations.
template <MInt nDim>
MFloat& FvCellCollector<nDim>::restrictedVar(const MInt id, const MInt varId) {
// Prevent accidental compilation without support for SoA layout
#ifdef FVCELLCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_VARIABLE_ID_ACCESSOR(varId);
  return m_restrictedVars[id * noCVariables() + varId];
}
/// Accessor for variables after restriction during multigrid computations (const version).
template <MInt nDim>
MFloat FvCellCollector<nDim>::restrictedVar(const MInt id, const MInt varId) const {
// Prevent accidental compilation without support for SoA layout
#ifdef FVCELLCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_VARIABLE_ID_ACCESSOR(varId);
  return m_restrictedVars[id * noCVariables() + varId];
}

/// Accessor for stored slopes.
template <MInt nDim>
MFloat& FvCellCollector<nDim>::storedSlope(const MInt id, const MInt dimVar, const MInt dimDir) {
// Prevent accidental compilation without support for SoA layout
#ifdef FVCELLCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_PVARIABLE_ID_ACCESSOR(dimVar);
  ENSURE_VALID_DIRECTION_ID_ACCESSOR(dimDir);
  return m_storedSlopes[id * noPVariables() * nDim + dimVar * nDim + dimDir];
}
/// Accessor for stored slopes (const version).
template <MInt nDim>
MFloat FvCellCollector<nDim>::storedSlope(const MInt id, const MInt dimVar, const MInt dimDir) const {
// Prevent accidental compilation without support for SoA layout
#ifdef FVCELLCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif
  ENSURE_VALID_ID_ACCESSOR(id);
  ENSURE_VALID_PVARIABLE_ID_ACCESSOR(dimVar);
  ENSURE_VALID_DIRECTION_ID_ACCESSOR(dimDir);
  return m_storedSlopes[id * noPVariables() * nDim + dimVar * nDim + dimDir];
}

/// Set number of species and update number of variables
template <MInt nDim>
void FvCellCollector<nDim>::setNoCVariables(const MInt noCVariables_, const MInt noSpecies_) {
  m_noSpecies = noSpecies_;
  m_noCVariables = noCVariables_;
  return;
}

/// Update number of primitive variables
template <MInt nDim>
void FvCellCollector<nDim>::setNoPVariables(const MInt noPVariables_) {
  m_noPVariables = noPVariables_;
  return;
}

/// Update number of flux variables
template <MInt nDim>
void FvCellCollector<nDim>::setNoFVariables(const MInt noFVariables_) {
  m_noFVariables = noFVariables_;
  return;
}

/// Update number of additional variables
template <MInt nDim>
void FvCellCollector<nDim>::setNoAVariables(const MInt noAVariables_) {
  m_noAVariables = noAVariables_;
  return;
}

/// Activate or deactivate multilevel-relevant storage
template <MInt nDim>
MInt FvCellCollector<nDim>::isMultilevel(const MBool isMultilevel_) {
  const MInt oldIsMultilevel = m_isMultilevel;
  m_isMultilevel = isMultilevel_;
  return oldIsMultilevel;
}

/// Erase range of nodes such that they contain no sensible values anymore.
template <MInt nDim>
void FvCellCollector<nDim>::invalidate(const MInt begin, const MInt end) {
// Prevent accidental compilation without support for SoA layout
#ifdef FVCELLCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif

  // oldVariables
  fill_invalid(m_oldVariables, begin, end, noCVariables());
  // variables
  fill_invalid(m_variables, begin, end, noCVariables());
  // pvariables
  fill_invalid(m_pvariables, begin, end, noPVariables());
  // right hand side
  fill_invalid(m_rightHandSide, begin, end, noFVariables());
  // slopes
  fill_invalid(m_slopes, begin, end, noPVariables() * nDim);
  // bndry cell ids, -1 as "invalid" value
  // mapping -1: internal, -2: ghost cell
  fill_invalid(m_bndryCellIds, begin, end, 1, -1);
  // noRcnstrctnNghbrIds
  fill_invalid(m_noRcnstrctnNghbrIds, begin, end);
  // rcnstrctnNghbrId
  fill_invalid(m_rcnstrctnNghbrIds, begin, end, noRecNghbrs());
  // m_reconstructionData
  fill_invalid(m_reconstructionData, begin, end);
  // nghbrInterface
  fill_invalid(m_nghbrInterface, begin, end);

  // reactionRates
  if(hasReactionRates()) fill_invalid(m_reactionRates, begin, end, noReactionRates());
  if(hasReactionRatesBackup()) fill_invalid(m_reactionRatesBackup, begin, end, noReactionRates());
  // psi
  if(hasPsi()) fill_invalid(m_psi, begin, end);
  // speciesReactionRates
  if(hasReactionRates()) fill_invalid(m_speciesReactionRates, begin, end, noSpecies());
  if(hasCellCenterMeanMolarWeight()) fill_invalid(m_cellCenterMeanMolarWeight, begin, end);
  if(hasCellCenterGamma()) fill_invalid(m_cellCenterGamma, begin, end);

  // implicit coefficients
  if(isEEGas()) fill_invalid(m_implicitCoefficients, begin, end, noImplicitCoefficients());
  // additional variables
  fill_invalid(m_avariables, begin, end, noAVariables());

  if(hasDualTS()) {
    // dt1Variables
    fill_invalid(m_dt1Variables, begin, end, noCVariables());
    // dt2Variables
    fill_invalid(m_dt2Variables, begin, end, noCVariables());
    fill_invalid(m_localTimeStep_, begin, end);
  } else if(hasLocalTS()) {
    fill_invalid(m_localTimeStep_, begin, end);
  }

  // spongeFactor
  fill_invalid(m_spongeFactor, begin, end);
  // spongeFactorStart
  fill_invalid(m_spongeFactorStart, begin, end);
  // spongeBndryIds
  fill_invalid(m_spongeBndryIds, begin, end, nDim);

  // Properties
  fill_invalid(m_properties, begin, end);
  // Coordinates
  fill_invalid(m_coordinates, begin * nDim, end * nDim);

  // Level
  fill_invalid(m_levels, begin, end);

  // Cell volume
  fill_invalid(m_cellVolumes, begin, end);
  // Inverse cell volume
  fill_invalid(m_FcellVolumes, begin, end);

  // Multilevel-related data
  if(isMultilevel()) {
    // Coarse grid correction
    fill_invalid(m_tau, begin, end, noCVariables());
    // Restricted RHS
    fill_invalid(m_restrictedRHS, begin, end, noCVariables());
    // Variables after restriction
    fill_invalid(m_restrictedVars, begin, end, noCVariables());
    // Temporarily stored slopes
    fill_invalid(m_storedSlopes, begin, end, noPVariables() * nDim);
  }
}


/// Helper function for rawCopy(). Destination may refer to beginning or end of target range.
template <MInt nDim>
template <class Functor, class T>
void FvCellCollector<nDim>::rawCopyGeneric(Functor&& c, const T& source, const MInt begin, const MInt end,
                                           const MInt destination) {
// Prevent accidental compilation without support for SoA layout
#ifdef FVCELLCOLLECTOR_SOA_MEMORY_LAYOUT
#error Missing implementation for structure-of-arrays memory layout.
#endif

  // oldVariables
  copyData(source.m_oldVariables, m_oldVariables, c, begin, end, destination, noCVariables());
  // variables
  copyData(source.m_variables, m_variables, c, begin, end, destination, noCVariables());
  // pvariables
  copyData(source.m_pvariables, m_pvariables, c, begin, end, destination, noPVariables());
  // right hand side
  copyData(source.m_rightHandSide, m_rightHandSide, c, begin, end, destination, noFVariables());
  // slopes
  copyData(source.m_slopes, m_slopes, c, begin, end, destination, noPVariables() * nDim);
  // bndryCellIds
  copyData(source.m_bndryCellIds, m_bndryCellIds, c, begin, end, destination);
  // noRcnstrctnNghbrIds
  copyData(source.m_noRcnstrctnNghbrIds, m_noRcnstrctnNghbrIds, c, begin, end, destination);
  // rcnstrctnNghbrId
  copyData(source.m_rcnstrctnNghbrIds, m_rcnstrctnNghbrIds, c, begin, end, destination, noRecNghbrs());
  // reconstructionData
  copyData(source.m_reconstructionData, m_reconstructionData, c, begin, end, destination);
  // nghbrInterface
  copyData(source.m_nghbrInterface, m_nghbrInterface, c, begin, end, destination);
  // reactionRate
  if(hasReactionRates())
    copyData(source.m_reactionRates, m_reactionRates, c, begin, end, destination, noReactionRates());
  if(hasReactionRatesBackup())
    copyData(source.m_reactionRatesBackup, m_reactionRatesBackup, c, begin, end, destination, noReactionRates());
  // psi
  if(hasPsi()) copyData(source.m_psi, m_psi, c, begin, end, destination);
  // speciesReactionRate
  if(hasReactionRates())
    copyData(source.m_speciesReactionRates, m_speciesReactionRates, c, begin, end, destination, noSpecies());
  if(hasCellCenterMeanMolarWeight())
    copyData(source.m_cellCenterMeanMolarWeight, m_cellCenterMeanMolarWeight, c, begin, end, destination);
  if(hasCellCenterGamma()) copyData(source.m_cellCenterGamma, m_cellCenterGamma, c, begin, end, destination);

  // implicit coefficient
  if(isEEGas())
    copyData(source.m_implicitCoefficients, m_implicitCoefficients, c, begin, end, destination,
             noImplicitCoefficients());

  // additional variables
  copyData(source.m_avariables, m_avariables, c, begin, end, destination, noAVariables());

  if(hasDualTS()) {
    // dt1Variables
    copyData(source.m_dt1Variables, m_dt1Variables, c, begin, end, destination, noCVariables());
    // dt2Variables
    copyData(source.m_dt2Variables, m_dt2Variables, c, begin, end, destination, noCVariables());
    copyData(source.m_localTimeStep_, m_localTimeStep_, c, begin, end, destination);
  } else if(hasLocalTS()) {
    copyData(source.m_localTimeStep_, m_localTimeStep_, c, begin, end, destination);
  }

  // spongeFactor
  copyData(source.m_spongeFactor, m_spongeFactor, c, begin, end, destination);
  // spongeFactorStart
  copyData(source.m_spongeFactorStart, m_spongeFactorStart, c, begin, end, destination);
  // spongeBbndryIds
  copyData(source.m_spongeBndryIds, m_spongeBndryIds, c, begin, end, destination, nDim);
  // Coordinates
  copyData(source.m_coordinates, m_coordinates, c, begin, end, destination, nDim);

  // Properties
  copyData(source.m_properties, m_properties, c, begin, end, destination);

  // Level
  copyData(source.m_levels, m_levels, c, begin, end, destination);

  // Cell volume
  copyData(source.m_cellVolumes, m_cellVolumes, c, begin, end, destination);

  // Inverse cell volume
  copyData(source.m_FcellVolumes, m_FcellVolumes, c, begin, end, destination);

  // Multilevel-related data
  if(isMultilevel()) {
    // Coarse grid correction
    copyData(source.m_tau, m_tau, c, begin, end, destination, noCVariables());
    // Restricted RHS
    copyData(source.m_restrictedRHS, m_restrictedRHS, c, begin, end, destination, noCVariables());
    // Variables after restriction
    copyData(source.m_restrictedVars, m_restrictedVars, c, begin, end, destination, noCVariables());
    // Temporarly stored slopes
    copyData(source.m_storedSlopes, m_storedSlopes, c, begin, end, destination, noPVariables() * nDim);
  }
}

} // namespace collector
} // namespace fv
} // namespace maia


// Undefine macros that should not be used outside this file
#undef FVCELLCOLLECTOR_SANITY_CHECKS_ACCESSORS
#undef ENSURE_VALID_ID_ACCESSOR
#undef ENSURE_VALID_VARIABLE_ID_ACCESSOR
#undef ENSURE_VALID_PVARIABLE_ID_ACCESSOR
#undef ENSURE_VALID_FVARIABLE_ID_ACCESSOR
#undef ENSURE_VALID_AVARIABLE_ID_ACCESSOR
#undef ENSURE_VALID_DIRECTION_ID_ACCESSOR
#undef ENSURE_VALID_RECNGHBR_ID_ACCESSOR
#undef ENSURE_VALID_REACTION_ID_ACCESSOR
#undef ENSURE_VALID_COORDINATE_DIR_ACCESSOR
#undef ENSURE_VALID_PROPERTY_ACCESSOR
#undef ENSURE_VALID_IMPLICIT_COEFFICIENT_ID_ACCESSOR

#endif // ifndef FVCELLCOLLECTOR_H_
