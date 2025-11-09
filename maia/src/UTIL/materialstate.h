// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef MAIA_MATERIALSTATE_H
#define MAIA_MATERIALSTATE_H
#include <functional>
#include "globals.h"
#include "solver.h"


template <MInt nDim>
class MaterialState {
 private:
  friend class LPT<nDim>;

  const MInt solverId;

  MFloat m_airCp = 1005;

  MFloat m_gamma = -1;
  MFloat m_gammaMinusOne = -1;
  MFloat m_gasConstant = -1;

  MFloat m_particleTemperature = 293.15;
  MFloat m_temperatureLower = -MFloatMax;
  MFloat m_temperatureUpper = MFloatMax;

  MFloat m_molarMass = -1;
  MFloat m_molarWeightRatio = 1;

  MFloat m_boilingPoint = 300.0;
  MFloat m_bpRefPressure = 101325;

  std::array<MFloat, 4> m_particleDensity = {-1.0, 0.0, 0.0, 0.0};
  std::array<MFloat, 4> m_cp = {0.0, 0.0, 0.0, 0.0};
  std::array<MFloat, 4> m_my = {0.0, 0.0, 0.0, 0.0};
  std::array<MFloat, 4> m_thCond = {0.0, 0.0, 0.0, 0.0};
  std::array<MFloat, 4> m_diffC = {0.0, 0.0, 0.0, 0.0};
  std::array<MFloat, 4> m_lH_ev = {0.0, 0.0, 0.0, 0.0};
  std::array<MFloat, 4> m_surfaceTension = {1.0, 0.0, 0.0, 0.0};
  std::array<MFloat, 4> m_liquidMy = {0.0, 0.0, 0.0, 0.0};
  std::array<MFloat, 4> m_liquidThCond = {0.0, 0.0, 0.0, 0.0};

  std::array<MFloat, 4> m_fluidThCond = {3.227E-3, 8.3894E-5, -1.9858E-8, 0.0};
  std::array<MFloat, 4> m_Pr = {0.815, -0.0004958, 0.0000004514, 0.0};

  MFloat m_densityRatio = 1;
  MFloat m_ambientDensity = 1;

  std::function<MFloat(const MFloat T)> m_viscosityFunction;
  MFloat m_nu = -1.0;

 public:
  MaterialState(const MInt solverId_) : solverId(solverId_) { readProperties(); }

  MFloat m_viscosityFactor = -1.0;
  MFloat m_temperatureFactor = -1.0;
  MFloat m_sutherlandConstant = -1.0;
  MFloat m_sutherlandPlusOne = -1.0;

  // check given temperature against given temperature range and return valid value
  inline MFloat checkWithTemperatureRange(const MFloat temperature_) {
    if(temperature_ < m_temperatureLower) {
      return m_temperatureLower;
    }
    if(temperature_ > m_temperatureUpper) {
      return m_temperatureUpper;
    }
    return temperature_;
  }

  // thermal conductivity of disperse phase in continous phase
  inline MFloat thermalConductivity(const MFloat temperature_) {
    MFloat temperature = checkWithTemperatureRange(temperature_);
    return m_thCond[0] + m_thCond[1] * temperature + m_thCond[2] * POW2(temperature) + m_thCond[3] * POW3(temperature);
  }

  // thermal conductivity of liquid phase in continous phase
  inline MFloat liquidThermalConductivity(const MFloat temperature_) {
    MFloat temperature = checkWithTemperatureRange(temperature_);
    return m_liquidThCond[0] + m_liquidThCond[1] * temperature + m_liquidThCond[2] * POW2(temperature)
           + m_liquidThCond[3] * POW3(temperature);
  }

  // diffusion Coefficient
  inline MFloat diffusionCoefficient(const MFloat temperature_) {
    MFloat temperature = checkWithTemperatureRange(temperature_);
    return m_diffC[0] + m_diffC[1] * temperature + m_diffC[2] * POW2(temperature) + m_diffC[3] * POW3(temperature);
  }

  // latent heat of evaporation of the disperse phase
  inline MFloat latentHeatEvap(const MFloat temperature_) {
    MFloat temperature = checkWithTemperatureRange(temperature_);
    return m_lH_ev[0] + m_lH_ev[1] * temperature + m_lH_ev[2] * POW2(temperature) + m_lH_ev[3] * POW3(temperature);
  }

  // dynamic viscosity as a function of temperature based on coefficients for dynamic viscosity
  inline MFloat dynamicViscosity(const MFloat temperature_) {
    MFloat temperature = checkWithTemperatureRange(temperature_);
    return m_my[0] + m_my[1] * temperature + m_my[2] * POW2(temperature) + m_my[3] * POW3(temperature);
  }

  // dynamic viscosity of liquid phase as a function of temperature based on coefficients for dynamic viscosity of
  // liquid phase
  inline MFloat liquidDynamicViscosity(const MFloat temperature_) {
    MFloat temperature = checkWithTemperatureRange(temperature_);
    return m_liquidMy[0] + m_liquidMy[1] * temperature + m_liquidMy[2] * POW2(temperature)
           + m_liquidMy[3] * POW3(temperature);
  }

  // thermal conductivity of the continous-phase
  inline MFloat airThermalConductivity(const MFloat temperature_) {
    MFloat temperature = checkWithTemperatureRange(temperature_);
    return m_fluidThCond[0] + m_fluidThCond[1] * temperature + m_fluidThCond[2] * POW2(temperature)
           + m_fluidThCond[3] * POW3(temperature);
  }

  // Pr-number of the continous-phase
  inline MFloat airPrandtl(const MFloat temperature_) {
    MFloat temperature = checkWithTemperatureRange(temperature_);
    return m_Pr[0] + m_Pr[1] * temperature + m_Pr[2] * POW2(temperature) + m_Pr[3] * POW3(temperature);
  }

  // dynamic viscosity of the continous phase
  inline MFloat dynViscosityFun(const MFloat temperature_) {
    MFloat temperature = checkWithTemperatureRange(temperature_);
    return m_viscosityFunction(temperature);
  }

  // NOTE: either particle surface tension (dimensional case)
  // or 1 as always equal to the constant reference value (non-dimensional case)
  // Surface Tension as a function of temperature based on coefficients for Surface Tension
  inline MFloat spraySurfaceTension(const MFloat temperature_ = 1) {
    MFloat temperature = checkWithTemperatureRange(temperature_);
    return m_surfaceTension[0] + m_surfaceTension[1] * temperature + m_surfaceTension[2] * POW2(temperature)
           + m_surfaceTension[3] * POW3(temperature);
  }

  // NOTE: either particle temperture (dimensional case)
  // or temperature ratio (non-dimensional case)
  inline MFloat T() { return m_particleTemperature; }

  // NOTE: either particle density (dimensional case)
  // or density ratio and identically with densityRatio (non-dimensional case)
  // particle density as a function of temperature based on coefficients for particle density
  inline MFloat density(const MFloat temperature_ = 1) {
    MFloat temperature = checkWithTemperatureRange(temperature_);
    return m_particleDensity[0] + m_particleDensity[1] * temperature + m_particleDensity[2] * POW2(temperature)
           + m_particleDensity[3] * POW3(temperature);
  }

  inline MFloat densityRatio() { return m_densityRatio; }

  inline MFloat ambientDensityRatio() { return m_densityRatio / m_ambientDensity; }

  // NOTE: either particle cp (dimensional case) particle/air ratio (non-dimensional case)
  // particle cp as a function of temperature based on coefficients for particle cp
  inline MFloat cp(const MFloat temperature_ = 1) {
    MFloat temperature = checkWithTemperatureRange(temperature_);
    return m_cp[0] + m_cp[1] * temperature + m_cp[2] * POW2(temperature) + m_cp[3] * POW3(temperature);
  }

  // NOTE: either particle surface tension (dimensional case)
  // or 1 as always equal to the constant reference value (non-dimensional case)
  inline MFloat airCp() { return m_airCp; }

  // identically function for non-dimensional and dimensional
  // as purely property based!
  inline MFloat boilingPoint() { return m_boilingPoint; }

  // identically function for non-dimensional and dimensional
  // as purely property based!
  inline MFloat bpRefPressure() { return m_bpRefPressure; }

  // NOTE: either particle specific gas constant (dimensional case)
  // or 1/gamma * molarWeight-Ration (non-dimensional case)
  inline MFloat gasConstant() { return m_gasConstant; }

  // identically function for non-dimensional and dimensional
  // as purely property based!
  inline MFloat molWeightRatio() { return m_molarWeightRatio; }

  // NOTE: either unity (dimensional case)
  // or gamma -1 (non-dimensional case)
  inline MFloat gammaMinusOne() { return m_gammaMinusOne; }


  void readProperties() {
    MBool heatCoupling = false;
    heatCoupling = Context::getSolverProperty<MBool>("particleHeatCoupling", solverId, AT_, &heatCoupling);

    MBool evaporation = false;
    evaporation = Context::getSolverProperty<MBool>("particleEvaporation", solverId, AT_, &evaporation);

    MBool nonDimensional = false;
    nonDimensional = Context::getSolverProperty<MBool>("nonDimensionaliseLPT", solverId, AT_, &nonDimensional);

    const std::array<MString, 4> suffix = {"_a", "_b", "_c", "_d"};

    //----------------------------------------- particle related ------------------------------------

    /*! \page propertyPageMaterial
    \section material_Density
    <code>MFloat MaterialState::m_particleDensity</code>\n
    default = -1 \n \n
    Sets the reference density for the different phases.\n
    Keywords: <i>PARTICLE</i>, <i>FINITE VOLUME</i>, <i>MATERIAL</i>
    */
    // NOTE: can be dimensional or non-dimensionalized by the reference density
    if(Context::propertyExists("particleDensity")) {
      m_particleDensity[0] =
          Context::getSolverProperty<MFloat>("particleDensity", solverId, AT_, &m_particleDensity[0]);
    } else {
      m_particleDensity[0] =
          Context::getSolverProperty<MFloat>("particleDensity_a", solverId, AT_, &m_particleDensity[0]);
    }

    /*! \page propertyPageMaterial
    \section material_Temperature
    <code>MFloat MaterialState::m_particleTemperature</code>\n
    default = 293.15K \n \n
    Sets the reference density for the different phases.\n
    Keywords: <i>PARTICLE</i>, <i>FINITE VOLUME</i>
    */
    // NOTE: can be dimensional or non-dimensionalized by the reference Temperature
    // default temperature is set to 20C(293.15K)
    m_particleTemperature =
        Context::getSolverProperty<MFloat>("particleTemperature", solverId, AT_, &m_particleTemperature);

    /*! \page propertyPageMaterial
    \section material_TemperatureRange
    <code>MFloat MaterialState::m_temperatureLower</code>\n
    default = 0.1 \n \n
    Sets the sets the boundaries of the temperature region.\n
    Keywords: <i>PARTICLE</i>, <i>FINITE VOLUME</i>
    */
    // NOTE: can be dimensional or non-dimensionalized by the reference Temperature
    m_temperatureLower = Context::getSolverProperty<MFloat>("temperatureLowerBnd", solverId, AT_, &m_temperatureLower);
    m_temperatureUpper = Context::getSolverProperty<MFloat>("temperatureUpperBnd", solverId, AT_, &m_temperatureUpper);


    /*! \page propertyPageMaterial
    \section material_CP
    <code>MFloat MaterialState::m_cp</code>\n
    default = 1005 \n \n
    Sets the heat capacity of the different species.\n
    Keywords: <i>PARTICLE</i>, <i>FINITE VOLUME</i>
    */
    static constexpr MFloat defaultAirCP = 1005;
    // NOTE: can be dimensional or non-dimensionalized by the fluid cp in the reference state
    if(Context::propertyExists("particleCP")) {
      m_cp[0] = Context::getSolverProperty<MFloat>("particleCP", solverId, AT_, &defaultAirCP);
    } else {
      m_cp[0] = Context::getSolverProperty<MFloat>("particleCP_a", solverId, AT_, &defaultAirCP);
    }

    /*! \page propertyPageMaterial
    \section sprayLiquidSurfaceTension
    <code>MFloat MaterialState::m_spraySurfaceTension</code>\n
    default = 1 \n \n
    Surface tension of the sprayed liquid \n
    Keywords: <i>PARTICLE</i>, <i>SPRAY</i>
    */
    // NOTE: always equal to unity in the non-dimensional case!
    if(Context::propertyExists("sprayLiquidSurfaceTension")) {
      m_surfaceTension[0] =
          Context::getSolverProperty<MFloat>("sprayLiquidSurfaceTension", solverId, AT_, &m_surfaceTension[0]);
    } else {
      m_surfaceTension[0] = Context::getSolverProperty<MFloat>("liquidSurfTens_a", solverId, AT_, &m_surfaceTension[0]);
    }

    if(heatCoupling || evaporation) {
      /*! \page propertyPageMaterial
        \section material_BPRefPressure
      <code>MFloat MaterialState::m_bpRefPressure</code>\n
      default = 1bar \n \n
      Set the reference pressure for the Clausius-Clapeyron equation\n
      Keywords: <i>PARTICLE</i>, <i>FINITE VOLUME</i>
      */
      if(nonDimensional) {
        // reference pressure
        m_bpRefPressure = -1;
      } else {
        // 1 bar
        m_bpRefPressure = 100000;
      }
      m_bpRefPressure = Context::getSolverProperty<MFloat>("particleBPRefPressure", solverId, AT_, &m_bpRefPressure);

      /*! \page propertyPageMaterial
      \section boilingPoint
      <code>MFloat MaterialState::m_boilingPoint</code>\n
      default = 58K (oxygen) \n \n
      Sets the boiling point of the different species.\n
      Keywords: <i>PARTICLE</i>, <i>FINITE VOLUME</i>
      */
      static constexpr MFloat defaultAirTB = 58;
      m_boilingPoint = Context::getSolverProperty<MFloat>("particleBoilingPoint", solverId, AT_, &defaultAirTB);

      /*! \page propertyPageMaterial
      \section material_particleDensity
      <code>MFloat MaterialState::m_particleDensity</code>\n
      default = {0.0, 0.0, 0.0, 0.0} \n \n
      Sets the coefficients for the 3rd Degree polynom equation for the particle density.\n
      Keywords: <i>PARTICLE</i>, <i>FINITE VOLUME</i>, <i>MATERIAL</i>
      */
      for(MInt i = 0; i < 4; ++i) {
        m_particleDensity[i] =
            Context::getSolverProperty<MFloat>("particleDensity" + suffix[i], solverId, AT_, &m_particleDensity[i]);
      }

      /*! \page propertyPageMaterial
      \section material_cp
      <code>MFloat MaterialState::m_cp</code>\n
      default = {0.0, 0.0, 0.0, 0.0} \n \n
      Sets the coefficients for the 3rd Degree polynom equation for the heat capacity.\n
      Keywords: <i>PARTICLE</i>, <i>FINITE VOLUME</i>, <i>MATERIAL</i>
      */
      for(MInt i = 0; i < 4; ++i) {
        m_cp[i] = Context::getSolverProperty<MFloat>("particleCP" + suffix[i], solverId, AT_, &m_cp[i]);
      }

      /*! \page propertyPageMaterial
      \section material_My
      <code>MFloat MaterialState::m_my</code>\n
      default = {0.0, 0.0, 0.0, 0.0} \n \n
      Sets the coefficients for the 3rd Degree polynom equation for the dynamic viscosity.\n
      Keywords: <i>PARTICLE</i>, <i>FINITE VOLUME</i>, <i>MATERIAL</i>
      */
      if(!Context::propertyExists("particleMy_a")) {
        mTerm(1, AT_, "Property particleMy_a not found, but is required for the Simulation!");
      }
      for(MInt i = 0; i < 4; ++i) {
        m_my[i] = Context::getSolverProperty<MFloat>("particleMy" + suffix[i], solverId, AT_, &m_my[i]);
      }

      /*! \page propertyPageMaterial
      \section material_ThCond
      <code>MFloat MaterialState::m_thCond</code>\n
      default = {0.0, 0.0, 0.0, 0.0} \n \n
      Sets the coefficients for the 3rd Degree polynom equation for the thermal conductivity.\n
      Keywords: <i>PARTICLE</i>, <i>FINITE VOLUME</i>, <i>MATERIAL</i>
      */
      if(!Context::propertyExists("particleThCond_a")) {
        mTerm(1, AT_, "Property particleThCond_a not found, but is required for the Simulation!");
      }
      for(MInt i = 0; i < 4; ++i) {
        m_thCond[i] = Context::getSolverProperty<MFloat>("particleThCond" + suffix[i], solverId, AT_, &m_thCond[i]);
      }

      /*! \page propertyPageMaterial
      \section material_DiffC
      <code>MFloat MaterialState::m_diffC</code>\n
      default = {0.0, 0.0, 0.0, 0.0} \n \n
      Sets the coefficients for the 3rd Degree polynom equation for the diffusion coefficient.\n
      Keywords: <i>PARTICLE</i>, <i>FINITE VOLUME</i>, <i>MATERIAL</i>
      */
      if(!Context::propertyExists("particleDiffC_a")) {
        mTerm(1, AT_, "Property particleDiffC_a not found, but is required for the Simulation!");
      }
      for(MInt i = 0; i < 4; ++i) {
        m_diffC[i] = Context::getSolverProperty<MFloat>("particleDiffC" + suffix[i], solverId, AT_, &m_diffC[i]);
      }

      /*! \page propertyPageMaterial
      \section material_LH_ev
      <code>MFloat MaterialState::m_lH_ev</code>\n
      default = {0.0, 0.0, 0.0, 0.0} \n \n
      Sets the coefficients for the 3rd Degree polynom equation for the latent heat of evaporation.\n
      Keywords: <i>PARTICLE</i>, <i>FINITE VOLUME</i>, <i>MATERIAL</i>
      */
      if(!Context::propertyExists("particleLH_ev_a")) {
        mTerm(1, AT_, "Property particleLH_ev_a not found, but is required for the Simulation!");
      }
      for(MInt i = 0; i < 4; ++i) {
        m_lH_ev[i] = Context::getSolverProperty<MFloat>("particleLH_ev" + suffix[i], solverId, AT_, &m_lH_ev[i]);
      }

      /*! \page propertyPageMaterial
      \section material_SurfaceTension
      <code>MFloat MaterialState::m_surfaceTension</code>\n
      default = {1.0, 0.0, 0.0, 0.0} \n \n
      Sets the coefficients for the 3rd Degree polynom equation for the Surface Tension.\n
      Keywords: <i>PARTICLE</i>, <i>FINITE VOLUME</i>, <i>MATERIAL</i>
      */
      for(MInt i = 0; i < 4; ++i) {
        m_surfaceTension[i] =
            Context::getSolverProperty<MFloat>("liquidSurfTens" + suffix[i], solverId, AT_, &m_surfaceTension[i]);
      }

      /*! \page propertyPageMaterial
      \section material_LiquidMy
      <code>MFloat MaterialState::m_liquidMy</code>\n
      default = {0.0, 0.0, 0.0, 0.0} \n \n
      Sets the coefficients for the 3rd Degree polynom equation for the dynamic viscosity
      of the liquid phase.\n
      Keywords: <i>PARTICLE</i>, <i>FINITE VOLUME</i>, <i>MATERIAL</i>
      */
      // Not yet stricly required -> to be changed !!!
      // if(!Context::propertyExists("liquidMy_a")) {
      //   mTerm(1, AT_, "Property liquidMy_a not found, but is required for the Simulation!");
      // }
      for(MInt i = 0; i < 4; ++i) {
        m_liquidMy[i] = Context::getSolverProperty<MFloat>("liquidMy" + suffix[i], solverId, AT_, &m_liquidMy[i]);
      }

      /*! \page propertyPageMaterial
      \section material_LiquidThCond
      <code>MFloat MaterialState::m_liquidThCond</code>\n
      default = {0.0, 0.0, 0.0, 0.0} \n \n
      Sets the coefficients for the 3rd Degree polynom equation for the thermal conductivity
      of the liquid phase.\n
      Keywords: <i>PARTICLE</i>, <i>FINITE VOLUME</i>, <i>MATERIAL</i>
      */
      // Not yet stricly required -> to be changed !!!
      // if(!Context::propertyExists("liquidThCond_a")) {
      //   mTerm(1, AT_, "Property liquidThCond_a not found, but is required for the Simulation!");
      // }
      for(MInt i = 0; i < 4; ++i) {
        m_liquidThCond[i] =
            Context::getSolverProperty<MFloat>("liquidThCond" + suffix[i], solverId, AT_, &m_liquidThCond[i]);
      }

      /*! \page propertyPageMaterial
      \section fluidThCond
      <code>MFloat MaterialState::m_fluidThCond</code>\n
      default = {3.227E-3, 8.3894E-5, -1.9858E-8, 0.0} \n \n
      Sets the coefficients for the 3rd Degree polynom equation for the thermal conductivity
      of the fluid/continous phase.\n
      Keywords: <i>PARTICLE</i>, <i>FINITE VOLUME</i>, <i>MATERIAL</i>
      */
      for(MInt i = 0; i < 4; ++i) {
        m_fluidThCond[i] =
            Context::getSolverProperty<MFloat>("fluidThCond" + suffix[i], solverId, AT_, &m_fluidThCond[i]);
      }

      /*! \page propertyPageMaterial
      \section Pr
      <code>MFloat MaterialState::m_Pr</code>\n
      default = {0.815, -0.0004958, 0.0000004514, 0.0} \n \n
      Sets the coefficients for the 3rd Degree polynom equation for the Pr-number
      of the fluid/continous phase.\n
      Keywords: <i>PARTICLE</i>, <i>FINITE VOLUME</i>, <i>MATERIAL</i>
      */
      for(MInt i = 0; i < 4; ++i) {
        m_Pr[i] = Context::getSolverProperty<MFloat>("fluidPr" + suffix[i], solverId, AT_, &m_Pr[i]);
      }
    }


    // --------------------------------- continues/ambient related --------------------------------

    /*! \page propertyPageMaterial
    \section m_molarWeightRatio
    <code>MFloat MaterialState::m_molarWeightRatio</code>\n
    default = 1 \n \n
    Non-dimensional molar weight ratio: disperse-phase/continous-phase (i.e. M_p / M_ref) \n
    Keywords: <i>PARTICLE</i>, <i>FINITE VOLUME</i>
    */
    m_molarWeightRatio = Context::getSolverProperty<MFloat>("molarWeightRatio", solverId, AT_, &m_molarWeightRatio);

    m_ambientDensity = Context::getSolverProperty<MFloat>("initialDensity", solverId, AT_, &m_ambientDensity);

    /*! \page propertyPageMaterial
    \section viscosityLaw
    <code>MString viscosityLaw</code>\n
    default = SUTHERLAND \n \n
    other options: CONSTANT
    Viscosity law used for the continuous phase \n
    Keywords: <i>PARTICLE</i>, <i>FINITE VOLUME</i>
    */
    MString viscosityLaw = "SUTHERLAND";
    viscosityLaw = Context::getSolverProperty<MString>("viscosityLaw", solverId, AT_, &viscosityLaw);

    /*! \page propertyPageMaterial
      \section referenceTemperature
      <code>MFloat referenceTemperature </code>\n
      default = <code>273.15</code>\n \n
    Reference temperature \f$ T_{\mathrm{ref}}\f$
    Used to scale the Sutherland's constant as follows: \f$ S/T_{\mathrm{ref}} \f$
      possible values are:
      <ul>
      <li>Non-negative floating point values</li>
      </ul>
      Keywords: <i>LPT, VARIABLES</i>
    */
    MFloat referenceTemperature = 273.15;
    referenceTemperature =
        Context::getSolverProperty<MFloat>("referenceTemperature", solverId, AT_, &referenceTemperature);

    /*! \page propertyPageMaterial
      \section sutherlandConstant
      <code>MFloat m_sutherlandConstant </code>\n
      default = <code>110.4 K</code>\n \n
      Sutherland's constant. Used by Sutherland's law.
      possible values are:
      <ul>
      <li>Non-negative floating point values</li>
      </ul>
      Keywords: <i>LPT, VARIABLES</i>
    */
    m_sutherlandConstant = 110.4;
    m_sutherlandConstant =
        Context::getSolverProperty<MFloat>("sutherlandConstant", solverId, AT_, &m_sutherlandConstant);

    m_sutherlandConstant /= referenceTemperature;
    m_sutherlandPlusOne = m_sutherlandConstant + F1;

    if(string2enum(viscosityLaw) == CONSTANT) {
      /*! \page propertyPageMaterial
      \section ambientViscosity
      <code>MFloat MaterialState::m_nu</code>\n
      default = none \n \n
      Constant viscosity of the continuous phase \n
      Keywords: <i>PARTICLE</i>, <i>FINITE VOLUME</i>
      */
      m_nu = Context::getSolverProperty<MFloat>("ambientViscosity", solverId, AT_);
    }

    if(!nonDimensional) {
      /*! \page propertyPageMaterial
      \section molarMass
      <code>MFloat MaterialState::m_molarMass</code>\n
      default = 0.02896 kg/mol (air) \n \n
      Sets the molar weights of the different species.\n
      Keywords: <i>PARTICLE</i>, <i>FINITE VOLUME</i>
      */
      m_molarMass = 0.02896;
      m_molarMass = Context::getSolverProperty<MFloat>("particleMolarWeight", solverId, AT_, &m_molarMass);
      m_gasConstant = 8.3144598 / m_molarMass;

      m_temperatureFactor = 293.15;
      m_temperatureFactor =
          Context::getSolverProperty<MFloat>("ambientTemperature", solverId, AT_, &m_temperatureFactor);

      m_viscosityFactor =
          0.00001716 * pow(m_temperatureFactor / 273.15, 1.5) * (273.15 + 110.4) / (m_temperatureFactor + 110.4);
      m_viscosityFactor = Context::getSolverProperty<MFloat>("ambientDynViscosity", solverId, AT_, &m_viscosityFactor);

      if(Context::propertyExists("particleDensityRatio", solverId)) {
        m_densityRatio = Context::getSolverProperty<MFloat>("particleDensityRatio", solverId, AT_, &m_densityRatio);

      } else if(Context::propertyExists("particleDensity", solverId)
                && Context::propertyExists("ambientDensity", solverId)) {
        MFloat disperse = 1;
        disperse = Context::getSolverProperty<MFloat>("particleDensity", solverId, AT_, &disperse);
        MFloat continous = 1;
        continous = Context::getSolverProperty<MFloat>("ambientDensity", solverId, AT_, &continous);

        m_densityRatio = disperse / continous;

      } else {
        m_densityRatio = 1;
      }

      if(density(T()) < 0.0) {
        MFloat continous = -1;
        continous = Context::getSolverProperty<MFloat>("ambientDensity", solverId, AT_, &continous);
        if(continous > 0.0) {
          m_particleDensity[0] = m_densityRatio * continous;
          m_particleDensity[1] = 0.0;
          m_particleDensity[2] = 0.0;
          m_particleDensity[3] = 0.0;
        } else {
          // particleDensity should be densityRatio * rho_ref;
          mTerm(1, AT_, "Particle-density not set!");
        }
      }

      /*! \page propertyPageMaterial
        \section CP
        <code>MFloat MaterialState::m_airCp</code>\n
        default = 1005 \n \n
        Sets the heat capacity of the different species.\n
        Keywords: <i>PARTICLE</i>, <i>FINITE VOLUME</i>
      */
      // NOTE: always equal to unity in the non-dimensional case!
      m_airCp = Context::getSolverProperty<MFloat>("ambientCPG", solverId, AT_, &m_airCp);

      m_gammaMinusOne = 1;

      switch(string2enum(viscosityLaw)) {
        case SUTHERLAND:
          m_viscosityFunction = [&](const MFloat temperature) {
            return m_viscosityFactor * SUTHERLANDLAW(temperature / m_temperatureFactor);
          };
          break;
        case CONSTANT:
          m_viscosityFunction = [&](const MFloat) { return m_nu; };
          break;
        default:
          mTerm(1, AT_, "Invalid viscosity law");
          break;
      }
    } else {
      m_gamma = 1.4;
      m_gammaMinusOne = m_gamma - 1;
      m_gasConstant = 1 / (m_gamma * molWeightRatio());

      if(Context::propertyExists("ambientDensity", solverId)) {
        mTerm(1, AT_, "Set particle density instead!");
      }

      m_densityRatio = density(T());
      if(m_densityRatio < 0.0) {
        mTerm(1, AT_, "Set particle density!");
      }

      switch(string2enum(viscosityLaw)) {
        case SUTHERLAND:
          m_viscosityFunction = [&](const MFloat temperature) { return SUTHERLANDLAW(temperature); };
          break;
        case CONSTANT:
          m_viscosityFunction = [&](const MFloat) { return m_nu; };
          break;
        default:
          mTerm(1, AT_, "Invalid viscosity law");
          break;
      }
    }
  }
};

template class MaterialState<3>;

#endif // MAIA_MATERIALSTATE_H
