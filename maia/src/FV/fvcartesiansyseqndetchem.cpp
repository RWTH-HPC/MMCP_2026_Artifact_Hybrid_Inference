// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "fvcartesiansyseqndetchem.h"
#include "MEMORY/alloc.h"

#if defined(WITH_CANTERA)

/**
 * @brief Construct a new FvSysEqnDetChem<nDim>::FvSysEqnDetChem object
 *
 * @tparam nDim Number of dimensions
 * @param solverId Id of the solver
 * @param noSpecies Number of species
 */
template <MInt nDim>
FvSysEqnDetChem<nDim>::FvSysEqnDetChem(const MInt solverId, const MInt noSpecies)
  : FvSysEqnNS<nDim>(solverId, noSpecies) {
  Cantera::addDirectory("/aia/opt/cantera/build/data");

  CV = new ConservativeVariables(noSpecies);
  PV = new PrimitiveVariables(noSpecies);
  FV = new FluxVariables(noSpecies);

  readProperties();

  allocateNoDiffusionCoefficients(noSpecies);

  SC = new SurfaceCoefficients(noSpecies, m_noDiffusionCoefficients, m_noThermalDiffusionCoefficients);

  m_species = new SpeciesProperties(noSpecies, *this);
  m_NASA = new NASACoefficients(noSpecies, *this);

  MFloat m_eps = 1.0e-10;
  ASSERT(m_species->referenceTemp - m_NASA->referenceTemp < m_eps,
         "Unequal reference temperature values in NASA coefficients and species standard enthalpy of formation");
}

/**
 * @brief Construct a new FvSysEqnDetChem<nDim>::ConservativeVariables::ConservativeVariables object
 *
 * @tparam nDim Number of dimensions
 * @param noSpecies Number of species
 */
template <MInt nDim>
FvSysEqnDetChem<nDim>::ConservativeVariables::ConservativeVariables(const MInt noSpecies)
  : m_noSpecies(noSpecies), noVariables(nDim + 2 + noSpecies) {
  if(m_noSpecies > 0) {
    mAlloc(RHO_Y, m_noSpecies, "FvSysEqnNS::ConservativeVariables::RHO_Y", AT_);
    for(MUint i = 0; i < m_noSpecies; ++i) {
      RHO_Y[i] = RHO_C + i;
    }
  }
}

/**
 * @brief Construct a new FvSysEqnDetChem<nDim>::FluxVariables::FluxVariables object
 *
 * @tparam nDim Number of dimensions
 * @param noSpecies Number of species
 */
template <MInt nDim>
FvSysEqnDetChem<nDim>::FluxVariables::FluxVariables(const MInt noSpecies) : ConservativeVariables(noSpecies) {}

/**
 * @brief Destroy the FvSysEqnDetChem<nDim>::ConservativeVariables::ConservativeVariables object
 *
 * @tparam nDim Number of dimensions
 */
template <MInt nDim>
FvSysEqnDetChem<nDim>::ConservativeVariables::~ConservativeVariables() {
  mDeallocate(RHO_Y);
}

/**
 * @brief Construct a new FvSysEqnDetChem<nDim>::PrimitiveVariables::PrimitiveVariables object
 *
 * @tparam nDim Number of dimensions
 * @param noSpecies Number of species
 */
template <MInt nDim>
FvSysEqnDetChem<nDim>::PrimitiveVariables::PrimitiveVariables(const MInt noSpecies)
  : m_noSpecies(noSpecies), noVariables(nDim + 2 + noSpecies) {
  if(m_noSpecies > 0) {
    mAlloc(Y, m_noSpecies, "FvSysEqnNS::PrimitiveVariables::Y", AT_);
    for(MUint i = 0; i < m_noSpecies; ++i) {
      Y[i] = C + i;
    }
  }
}

/**
 * @brief Destroy the FvSysEqnDetChem<nDim>::PrimitiveVariables::PrimitiveVariables object
 *
 * @tparam nDim Number of dimensions
 */
template <MInt nDim>
FvSysEqnDetChem<nDim>::PrimitiveVariables::~PrimitiveVariables() {
  mDeallocate(Y);
}

/**
 * @brief Allocates the correct number of diffusion coefficients depending on the chosen diffusion model. Valid entries
 * are "Multi" for multicomponent diffusion model and "Mix" for mixture-averaged model.
 *
 * @tparam nDim Number of dimensions
 * @param noSpecies Number of species
 */
template <MInt nDim>
void FvSysEqnDetChem<nDim>::allocateNoDiffusionCoefficients(const MInt noSpecies) {
  switch(string2enum(m_transportModel)) {
    case Multi: {
      m_multiDiffusion = true;
      m_noDiffusionCoefficients = noSpecies * noSpecies;
      m_noThermalDiffusionCoefficients = noSpecies;
      m_noSurfaceCoefficients = 4 + m_noDiffusionCoefficients + m_noThermalDiffusionCoefficients;
    } break;
    case Mix: {
      m_multiDiffusion = false;
      m_noDiffusionCoefficients = noSpecies;
      m_noThermalDiffusionCoefficients = 0;
      m_noSurfaceCoefficients = 4 + m_noDiffusionCoefficients + m_noThermalDiffusionCoefficients;
    } break;
    default: {
      mTerm(1, AT_, "Diffusion model in properties file does not match any implemented model. Terminating.");
    } break;
  }
}

/**
 * @brief Construct a new FvSysEqnDetChem<nDim>::SpeciesProperties::SpeciesProperties object
 *
 * @tparam nDim
 * @param noSpecies
 * @param sysEqn
 */
template <MInt nDim>
FvSysEqnDetChem<nDim>::SpeciesProperties::SpeciesProperties(const MInt noSpecies, const FvSysEqnDetChem<nDim>& sysEqn) {
  mAlloc(molarMass, noSpecies, "FvSysEqnDetChem::SpeciesProperties::molarMass", AT_);
  mAlloc(fMolarMass, noSpecies, "FvSysEqnDetChem::SpeciesProperties::fMolarMass", AT_);
  mAlloc(specificGasConstant, noSpecies, "FvSysEqnDetChem::SpeciesProperties::specificGasConstant", AT_);
  mAlloc(fSpecificGasConstant, noSpecies, "FvSysEqnDetChem::SpeciesProperties::fSpecificGasConstant", AT_);
  mAlloc(standardHeatFormation, noSpecies, "FvSysEqnDetChem::SpeciesProperties::standardHeatFormation", AT_);

  this->getSpeciesProperties(noSpecies, sysEqn);
}

/**
 * @brief Destroy the FvSysEqnDetChem<nDim>::SpeciesProperties::SpeciesProperties object
 *
 * @tparam nDim
 */
template <MInt nDim>
FvSysEqnDetChem<nDim>::SpeciesProperties::~SpeciesProperties() {
  mDeallocate(molarMass);
  mDeallocate(fMolarMass);
  mDeallocate(specificGasConstant);
  mDeallocate(fSpecificGasConstant);
  mDeallocate(standardHeatFormation);
}

/**
 * @brief Gets important species information and stores it in the member variables of the species struct.
 *
 * @tparam nDim Number of dimensions
 * @param noSpecies Number of species
 * @param sysEqn Equation system
 */
template <MInt nDim>
void FvSysEqnDetChem<nDim>::SpeciesProperties::getSpeciesProperties(const MInt noSpecies,
                                                                    const FvSysEqnDetChem<nDim>& sysEqn) {
  std::shared_ptr<Cantera::Solution> sol(
      newSolution(sysEqn.m_reactionMechanism, sysEqn.m_phaseName, sysEqn.m_transportModel));
  std::shared_ptr<Cantera::ThermoPhase> gas(sol->thermo());

  const MInt numberMechanismSpecies = gas->nSpecies();

  if(numberMechanismSpecies != noSpecies) {
    std::cerr << "noSpecies in mechanism file: " << numberMechanismSpecies
              << ". noSpecies in simulation properties file: " << noSpecies << std::endl;
    mTerm(1, AT_,
          "Number of species included in the mechanism file does not correspond to the number of species defined "
          "in the properties-file. Terminating.");
  }

  gas->getMolecularWeights(molarMass); // SI: [kJ/kmol]]
  speciesName = gas->speciesNames();

  for(MUint s = 0; s < sysEqn.PV->m_noSpecies; s++) {
    standardHeatFormation[s] = gas->Hf298SS(s);                          // SI: [J/kmol]]
    molarMass[s] /= 1000.0;                                              // SI: [kg/mol]
    fMolarMass[s] = F1 / molarMass[s];                                   // SI: [mol/kg]
    standardHeatFormation[s] /= 1000.0;                                  // SI : [J/mol]
    standardHeatFormation[s] = standardHeatFormation[s] * fMolarMass[s]; // SI : [J/kg]
    specificGasConstant[s] = sysEqn.m_gasConstant * fMolarMass[s];
    fSpecificGasConstant[s] = F1 / specificGasConstant[s];
  }

  auto stringIterator = std::find(speciesName.begin(), speciesName.end(), majorSpecies);
  if(stringIterator == speciesName.end()) {
    mTerm(1, AT_, "Given major species name was not found in the species name vector. Terminating.");
  } else {
    majorSpeciesIndex = std::distance(speciesName.begin(), stringIterator);
  }

  // Map an index position with each species name
  for(MInt s = 0; s < noSpecies; s++) {
    speciesMap.insert(std::pair<std::string, MInt>(speciesName[s], s));
  }

#ifndef NDEBUG
  std::cout << "Detailed chemistry species properties." << std::endl;
  for(MInt s = 0; s < numberMechanismSpecies; s++) {
    std::cout << speciesName[s] << ". Molar mass: " << molarMass[s] << ". Heat formation: " << standardHeatFormation[s]
              << ". Specific gas constant: " << specificGasConstant[s] << std::endl;
  }
#endif
}

/**
 * @brief Construct a new FvSysEqnDetChem<nDim>::NASACoefficients::NASACoefficients object. This object stores all
 * the information about the NASA coefficients
 *
 * @tparam nDim Number of dimensions
 * @param noSpecies Number of species
 * @param sysEqn Equation system
 */
template <MInt nDim>
FvSysEqnDetChem<nDim>::NASACoefficients::NASACoefficients(const MInt noSpecies, const FvSysEqnDetChem<nDim>& sysEqn) {
  mAlloc(lowTemp, noNASACoefficients * noSpecies, "FvSysEqnDetChem::NASACoefficients::lowTemp", AT_);
  mAlloc(highTemp, noNASACoefficients * noSpecies, "FvSysEqnDetChem::NASACoefficients::highTemp", AT_);

  mAlloc(integralLowTemp, noNASACoefficientsCpPolynomial * noSpecies,
         "FvSysEqnDetChem::NASACoefficients::integralLowTemp", AT_);
  mAlloc(integralHighTemp, noNASACoefficientsCpPolynomial * noSpecies,
         "FvSysEqnDetChem::NASACoefficients::integralHighTemp", AT_);

  mAlloc(lowTempIntegrationConstantsEnergy, noSpecies,
         "FvSysEqnDetChem::NASACoefficients::lowTempIntegrationConstantsEnergy", AT_);
  mAlloc(highTempIntegrationConstantsEnergy, noSpecies,
         "FvSysEqnDetChem::NASACoefficients::highTempIntegrationConstantsEnergy", AT_);

  mAlloc(lowTempIntegrationConstantsEnthalpy, noSpecies,
         "FvSysEqnDetChem::NASACoefficients::lowTempIntegrationConstantsEnthalpy", AT_);
  mAlloc(highTempIntegrationConstantsEnthalpy, noSpecies,
         "FvSysEqnDetChem::NASACoefficients::highTempIntegrationConstantsEnthalpy", AT_);

  this->getNASACoefficients(noSpecies, sysEqn);
  this->computeSensibleEnergyIntegrationConstants(sysEqn);
}

/**
 * @brief Destroy the FvSysEqnDetChem<nDim>::NASACoefficients::NASACoefficients object
 *
 * @tparam nDim Number of dimensions
 */
template <MInt nDim>
FvSysEqnDetChem<nDim>::NASACoefficients::~NASACoefficients() {
  mDeallocate(lowTemp);
  mDeallocate(highTemp);

  mDeallocate(integralLowTemp);
  mDeallocate(integralHighTemp);

  mDeallocate(lowTempIntegrationConstantsEnergy);
  mDeallocate(highTempIntegrationConstantsEnergy);

  mDeallocate(lowTempIntegrationConstantsEnthalpy);
  mDeallocate(highTempIntegrationConstantsEnthalpy);
}

/**
 * @brief Gets the species 7-NASA Coefficients from the given mechanism file.
The NASA Coefficients are then used to compute the temperature-dependant specific heat capacity
The integral coefficients are the NASA COEFFICIENTS multiplied by 1/(n + 1) and are used
for the integration of the heat capacities to compute the sensible energy.
 *
 * @tparam nDim Number of dimensions
 * @param noSpecies Number of species
 * @param sysEqn Equation system
 */
template <MInt nDim>
void FvSysEqnDetChem<nDim>::NASACoefficients::getNASACoefficients(const MInt noSpecies,
                                                                  const FvSysEqnDetChem<nDim>& sysEqn) {
  std::shared_ptr<Cantera::Solution> sol(
      Cantera::newSolution(sysEqn.m_reactionMechanism, sysEqn.m_phaseName, sysEqn.m_transportModel));
  std::shared_ptr<Cantera::ThermoPhase> gas(sol->thermo());

  MFloat* totalNASACoefficients = nullptr;
  mAlloc(totalNASACoefficients, totalNumberCoefficientsPerSpecies * noSpecies,
         "FvSysEqnDetChem::NASACoefficients::getNASACoefficients::totalNASACoefficients", AT_);

  MFloat* const NASACoeffs = ALIGNED_F(totalNASACoefficients);

  for(MUint s = 0; s < sysEqn.PV->m_noSpecies; s++) {
    MInt offset = totalNumberCoefficientsPerSpecies * s;
    MFloat* const speciesNASACoeffs = ALIGNED_F(NASACoeffs + offset);

    shared_ptr<Cantera::Species> species = gas->species(s);
    shared_ptr<Cantera::SpeciesThermoInterpType> thermoPhaseInf = species->thermo;

    // Dummy variables for Cantera function
    size_t n;
    doublereal tlow, thigh, pref;
    int type;

    thermoPhaseInf->reportParameters(n, type, tlow, thigh, pref, speciesNASACoeffs);
  }
  // Store the NASA coefficients for the low and high temperature regions
  for(MUint s = 0; s < sysEqn.PV->m_noSpecies; s++) {
    const MInt offset = totalNumberCoefficientsPerSpecies * s;
    const MInt offsetRegion = noNASACoefficients * s;
    for(MInt i = 0; i < totalNumberCoefficientsPerSpecies; i++) {
      if(i == 0) {
        continue;
      } else if((i > 0) && (i < (noNASACoefficients + 1))) {
        const MInt index = i - 1;
        highTemp[offsetRegion + index] = totalNASACoefficients[offset + i]; // seems OK
      } else {
        const MInt index = i - 1 - noNASACoefficients;
        lowTemp[offsetRegion + index] = totalNASACoefficients[offset + i]; // seems OK
      }
    }
  }
  // Compute and store integral coefficients
  for(MUint s = 0; s < sysEqn.PV->m_noSpecies; s++) {
    const MInt offsetRegion = noNASACoefficients * s;
    const MInt offsetIntegral = noNASACoefficientsCpPolynomial * s;
    for(MInt i = 0; i < 5; i++) {
      integralLowTemp[offsetIntegral + i] = lowTemp[offsetRegion + i] / (i + 1);
      integralHighTemp[offsetIntegral + i] = highTemp[offsetRegion + i] / (i + 1);
    }
  }

#ifndef NDEBUG
  for(MUint s = 0; s < sysEqn.PV->m_noSpecies; s++) {
    std::cout << sysEqn.m_species->speciesName[s] << std::endl;

    MInt offset = s * totalNumberCoefficientsPerSpecies;
    MInt offset2 = s * noNASACoefficients;
    MInt offset3 = s * noNASACoefficientsCpPolynomial;

    for(MInt i = 0; i < totalNumberCoefficientsPerSpecies; i++) {
      std::cout << totalNASACoefficients[offset + i] << " ";
    }
    std::cout << std::endl;

    for(MInt i = 0; i < noNASACoefficients; i++) {
      std::cout << highTemp[offset2 + i] << " ";
    }
    std::cout << std::endl;
    for(MInt i = 0; i < noNASACoefficientsCpPolynomial; i++) {
      std::cout << integralHighTemp[offset3 + i] << " ";
    }
    std::cout << std::endl;
    for(MInt i = 0; i < noNASACoefficients; i++) {
      std::cout << lowTemp[offset2 + i] << " ";
    }
    std::cout << std::endl;
    for(MInt i = 0; i < noNASACoefficientsCpPolynomial; i++) {
      std::cout << integralLowTemp[offset3 + i] << " ";
    }
    std::cout << std::endl;
  }
#endif

  mDeallocate(totalNASACoefficients);
}

/**
 * @brief Computes the sensible energy and enthalpy integration constants to reuse later (for low and high temp
regions). The constants are a result from definite integration. The low temp constant is computed with the low temp NASA
polynomial and reference temperature. The high temperature constant is composed of two different terms: one arises from
the integration with the low temperature NASA polynomial from reference temperature to the transition temperature. The
other term is the lower bound of the definite integral at the high temeprature region. These constants appear in the
computation of the sensible energy and are stored and reused to save computational time.
 *
 * @tparam nDim Number of dimensions
 * @param sysEqn Equation system
 */
template <MInt nDim>
void FvSysEqnDetChem<nDim>::NASACoefficients::computeSensibleEnergyIntegrationConstants(
    const FvSysEqnDetChem<nDim>& sysEqn) {
  for(MUint s = 0; s < sysEqn.PV->m_noSpecies; s++) {
    MFloat limitTRefEnergy = F0, limitLowTempRegionEnergy = F0, limitHighTempRegionEnergy = F0;
    MFloat limitTRefEnthalpy = F0, limitLowTempRegionEnthalpy = F0, limitHighTempRegionEnthalpy = F0;

    MFloat speciesSpecificGasConstant = sysEqn.m_species->specificGasConstant[s];
    MInt offsetIntegralCoeffs = noNASACoefficientsCpPolynomial * s;

    // Horner's rule for polynomials
    for(MInt i = 4; (i < 5) && (i >= 0); i--) {
      limitTRefEnergy = limitTRefEnergy * referenceTemp + integralLowTemp[offsetIntegralCoeffs + i];
      limitLowTempRegionEnergy = limitLowTempRegionEnergy * transitionTemp + integralLowTemp[offsetIntegralCoeffs + i];
      limitHighTempRegionEnergy =
          limitHighTempRegionEnergy * transitionTemp + integralHighTemp[offsetIntegralCoeffs + i];

      limitTRefEnthalpy = limitTRefEnthalpy * referenceTemp + integralLowTemp[offsetIntegralCoeffs + i];
      limitLowTempRegionEnthalpy =
          limitLowTempRegionEnthalpy * transitionTemp + integralLowTemp[offsetIntegralCoeffs + i];
      limitHighTempRegionEnthalpy =
          limitHighTempRegionEnthalpy * transitionTemp + integralHighTemp[offsetIntegralCoeffs + i];
    }

    limitTRefEnergy *= speciesSpecificGasConstant;
    limitLowTempRegionEnergy *= speciesSpecificGasConstant;
    limitHighTempRegionEnergy *= speciesSpecificGasConstant;

    limitTRefEnthalpy *= speciesSpecificGasConstant * referenceTemp;
    limitLowTempRegionEnthalpy *= speciesSpecificGasConstant * transitionTemp;
    limitHighTempRegionEnthalpy *= speciesSpecificGasConstant * transitionTemp;

    limitTRefEnergy = (limitTRefEnergy - speciesSpecificGasConstant) * referenceTemp;
    limitLowTempRegionEnergy = (limitLowTempRegionEnergy - speciesSpecificGasConstant) * transitionTemp;
    limitHighTempRegionEnergy = (limitHighTempRegionEnergy - speciesSpecificGasConstant) * transitionTemp;

    // Integration constants as defined are SUBSTRACTED in the corresponding equations
    lowTempIntegrationConstantsEnergy[s] = limitTRefEnergy;
    highTempIntegrationConstantsEnergy[s] = limitTRefEnergy - limitLowTempRegionEnergy + limitHighTempRegionEnergy;

    lowTempIntegrationConstantsEnthalpy[s] = limitTRefEnthalpy;
    highTempIntegrationConstantsEnthalpy[s] =
        limitTRefEnthalpy - limitLowTempRegionEnthalpy + limitHighTempRegionEnthalpy;
  }
}

/**
 * @brief Construct a new FvSysEqnDetChem<nDim>::SurfaceCoefficients::SurfaceCoefficients object
 *
 * @tparam nDim Number of dimensions
 * @param noSpecies Number of species
 * @param noDiffusionCoefficients Number of species diffusion coefficients. Equals N^2 for multicomponent model and N
 * for the mixture-averaged model, which is allocated in FvSysEqnDetChem<nDim>::allocateNoDiffusionCoefficients
 * @param noThermalDiffusionCoefficients Number of thermal diffusion coefficients, also allocated in
 * FvSysEqnDetChem<nDim>::allocateNoDiffusionCoefficients
 */
template <MInt nDim>
FvSysEqnDetChem<nDim>::SurfaceCoefficients::SurfaceCoefficients(const MInt noSpecies,
                                                                const MInt noDiffusionCoefficients,
                                                                const MInt noThermalDiffusionCoefficients)
  : m_noDiffusionCoefficients(noDiffusionCoefficients),
    m_noThermalDiffusionCoefficients(noThermalDiffusionCoefficients),
    m_noSurfaceCoefficients(4 + noDiffusionCoefficients + noThermalDiffusionCoefficients) {
  if(noSpecies > 0) {
    mAlloc(D, m_noDiffusionCoefficients, "FvSysEqnNS::SurfaceCoefficients::D", AT_);
    mAlloc(DT, m_noThermalDiffusionCoefficients, "FvSysEqnNS::SurfaceCoefficients::DT", AT_);
    for(MInt i = 0; i < m_noDiffusionCoefficients; ++i) {
      D[i] = D0 + i;
    }
    for(MInt i = 0; i < m_noThermalDiffusionCoefficients; ++i) {
      DT[i] = D0 + m_noDiffusionCoefficients + i;
    }
  }
}

/**
 * @brief Destroy the FvSysEqnDetChem<nDim>::SurfaceCoefficients::SurfaceCoefficients object
 *
 * @tparam nDim Number of dimensions
 */
template <MInt nDim>
FvSysEqnDetChem<nDim>::SurfaceCoefficients::~SurfaceCoefficients() {
  mDeallocate(D);
  mDeallocate(DT);
}

/**
 * @brief Reads important variable values from the properties.toml file.
 *
 * @tparam nDim Number of dimensions
 */
template <MInt nDim>
void FvSysEqnDetChem<nDim>::readProperties() {
  m_gasConstant = 8.314472;
  m_gasConstant = Context::getSolverProperty<MFloat>("gasConstant", m_solverId, AT_, &m_gasConstant);

  m_fGasConstant = F1 / m_gasConstant;

  m_reactionMechanism = Context::getSolverProperty<MString>("reactionMechanism", m_solverId, AT_, &m_reactionMechanism);

  m_transportModel = "Multi";
  m_transportModel = Context::getSolverProperty<MString>("transportModel", m_solverId, AT_, &m_transportModel);

  m_phaseName = Context::getSolverProperty<MString>("phaseName", m_solverId, AT_, &m_phaseName);

  m_computeSrfcCoeffsEveryRKStep = false;
  m_computeSrfcCoeffsEveryRKStep = Context::getSolverProperty<MBool>("computeSrfcCoeffsEveryRKStep", m_solverId, AT_,
                                                                     &m_computeSrfcCoeffsEveryRKStep);

  m_soretEffect = true;
  m_soretEffect = Context::getSolverProperty<MBool>("soretEffect", m_solverId, AT_, &m_soretEffect);

  m_fuelOxidizerStochiometricRatio = 2.0; // default for hydrogen + oxygen
  m_fuelOxidizerStochiometricRatio = Context::getSolverProperty<MFloat>("fuelOxidizerStochiometricRatio", m_solverId,
                                                                        AT_, &m_fuelOxidizerStochiometricRatio);

  m_oxidizer = "O2";
  m_oxidizer = Context::getSolverProperty<MString>("oxidizer", m_solverId, AT_, &m_oxidizer);

  m_fuel = "H2";
  m_fuel = Context::getSolverProperty<MString>("fuel", m_solverId, AT_, &m_fuel);
}

template class FvSysEqnDetChem<2>;
template class FvSysEqnDetChem<3>;

#else // WITH_CANTERA
// To allow compilation of SysEqnDetChem when CANTERA is not defined
template <MInt nDim>
FvSysEqnDetChem<nDim>::FvSysEqnDetChem(const MInt solverId, const MInt noSpecies)
  : FvSysEqnNS<nDim>(solverId, noSpecies) {
  // Does not allow the sysEqn to be used in a simulation if CANTERA has not been compiled
  mTerm(1, AT_,
        "Using detailed chemistry equation system without CANTERA compilation. Recompile enabling --with-cantera in "
        "configure.py file or select a different equation system. Terminating now.");
}

template class FvSysEqnDetChem<2>;
template class FvSysEqnDetChem<3>;
#endif
