// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "fvparticle.h"
#include "LPT/lptlib.h"

using namespace std;
using namespace maia::lpt;

template <MInt nDim, class SysEqn>
CouplerFvParticle<nDim, SysEqn>::CouplerFvParticle(const MInt couplingId, LPT<nDim>* particle, FvCartesianSolver* fv)
  : Coupling(couplingId), CouplingLpt<nDim, CouplingFv<nDim, SysEqn>>(couplingId, particle, fv), m_sysEqn(1, 0) {
  TRACE();

  m_particle = particle;
  m_fvSolver = fv;

  // get LPT-solver order
  m_noSolverSteps = 1;
  m_noSolverSteps = Context::getBasicProperty<MInt>("recipeMaxNoSteps", AT_, &m_noSolverSteps);

  const MString propNameLpt = "solverOrder_" + std::to_string(lpt().solverId());
  const MString propNameFv = "solverOrder_" + std::to_string(fvSolver().solverId());

  for(MInt step = 0; step < m_noSolverSteps; step++) {
    if(Context::getBasicProperty<MInt>(propNameLpt, AT_, step) > 0) {
      m_lptSolverOrder = step;
    }
    if(Context::getBasicProperty<MInt>(propNameFv, AT_, step) > 0) {
      m_fvSolverOrder = step;
    }
  }

  m_interLeafed = false;
  m_noSolutionSteps = 1;

  if(m_lptSolverOrder == m_fvSolverOrder) {
    m_interLeafed = true;
    m_noSolutionSteps = 5;
    lpt().m_noSolutionSteps = 5;
    if(lpt().domainId() == 0) {
      cerr << "Interleafed Fv - LPT execution at step " << m_lptSolverOrder << endl;
    }
  }

  initData();
}

/** \brief performs the coupling after solver initialization
 *  \author Tim Wegmann
 */
template <MInt nDim, class SysEqn>
void CouplerFvParticle<nDim, SysEqn>::init() {
  TRACE();

  checkProperties();

  initConversion();

  // check for any in-active ranks
  MInt noInactiveFv = 0;
  MInt noInactiveLPT = 0;
  if(!fvSolver().isActive()) {
    noInactiveFv++;
  }
  if(!lpt().isActive()) {
    noInactiveLPT++;
  }

  MPI_Allreduce(MPI_IN_PLACE, &noInactiveFv, 1, MPI_INT, MPI_SUM, globalMaiaCommWorld(), AT_, "MPI_IN_PLACE", "noInactiveFv");
  MPI_Allreduce(MPI_IN_PLACE, &noInactiveLPT, 1, MPI_INT, MPI_SUM, globalMaiaCommWorld(), AT_, "MPI_IN_PLACE",
                "noInactiveLPT");
  if(noInactiveFv > 0) {
    cerr0 << "Fv solver has " << noInactiveFv << " inactive ranks!" << endl;
  }
  if(noInactiveLPT > 0) {
    cerr0 << "LPT solver has " << noInactiveLPT << " inactive ranks!" << endl;
  }
}


/**
 * \fn void CouplerFvParticle<nDim, SysEqn>::initConversion()
 * \brief Calculates the conversion factor for different non-dimensionalisations
 *         in the FV and LPT solvers
 *  \author Tim Wegmann
 */
template <MInt nDim, class SysEqn>
void CouplerFvParticle<nDim, SysEqn>::initConversion() {
  // 1) read the following conversion factors as properties
  MFloat lengthFactor = 1;
  /*! \page propertiesCoupling
      \section lengthFactor
  <code>MFloat lengthFactor</code>\n
  default = 1 \n \n
  Set d_ref/L_ref the conversion between different reference length used in the non-dimensionalisation and the
  calculation of the solver specific Re-number in the FV (L_ref = L_grid) and LPT solver (L_ref = d_ref). \n Keywords:
  <i>PARTICLE</i>, <i>FINITE VOLUME</i>
  */
  lengthFactor = Context::getSolverProperty<MFloat>("lengthFactor", lpt().solverId(), AT_, &lengthFactor);

  MFloat velocityFactor = 1;
  /*! \page propertiesCoupling
      \section velocityFactor
  <code>MFloat velocityFactor</code>\n
  default = 1 \n \n
  Set a0/u_ref the conversion between different velocities used in the non-dimensionalisation and the calculation of the
  solver specific Re-number in the FV (a0) and LPT solver (u_ref). \n Keywords: <i>PARTICLE</i>, <i>FINITE VOLUME</i>
  */
  velocityFactor = Context::getSolverProperty<MFloat>("velocityFactor", lpt().solverId(), AT_, &velocityFactor);

  MFloat viscosityFactor = 1.0;
  /*! \page propertiesCoupling
   * \section viscosityFactor
  <code>MFloat viscosityFactor</code>\n
  default = 1.0 \n \n
  The conversion-factor for FV reference viscosity (mu_0 at T0) to LPT viscosity (mu_ref, at T_ref) as used in the
  solver-specific Re-numbers. The viscosityFactor itself is non-dimensional!

  \n Keywords: <i>PARTICLE</i>, <i>FINITE VOLUME</i>
  */
  viscosityFactor = Context::getSolverProperty<MFloat>("viscosityFactor", lpt().solverId(), AT_, &viscosityFactor);

  const MFloat particleRe = Context::getSolverProperty<MFloat>("Re", lpt().solverId(), AT_);

  // 2) compute the missing conversion factors
  const MFloat ReFactor = fvSolver().sysEqn().m_Re0 / particleRe;
  /*const */ MFloat densityFactor = ReFactor * lengthFactor / (velocityFactor * viscosityFactor);
  /*const */ MFloat pressureFactor = POW2(velocityFactor) * densityFactor;
  /*const */ MFloat temperatureFactor = sqrt(velocityFactor);


  //-----------------------------------------------------------------------------------
  // Fixes for dimensional LPT
  if(!lpt().m_nonDimensional) {
    MFloat specificGasConstant = 287.00283051433; // changes all particle testcases
    static constexpr MFloat defaultAir = 0.02896;
    MFloat molarMass = Context::getSolverProperty<MFloat>("ambientMolarWeight", lpt().solverId(), AT_, &defaultAir);
    MFloat T0 = 293.15;
    T0 = Context::getSolverProperty<MFloat>("ambientTemperature", lpt().solverId(), AT_, &T0);

    if(Context::propertyExists("gasConstant", lpt().solverId())) {
      specificGasConstant =
          Context::getSolverProperty<MFloat>("gasConstant", lpt().solverId(), AT_, &specificGasConstant);
      specificGasConstant = specificGasConstant / molarMass;
    }

    MFloat Tinfty = T0;
    if(Context::getSolverProperty<MInt>("initialCondition", lpt().solverId(), AT_) != 465) {
      Tinfty = T0 * fvSolver().sysEqn().temperature_IR(fvSolver().m_Ma);
    }

    MFloat density0 = -1;
    density0 = Context::getSolverProperty<MFloat>("ambientDensity", lpt().solverId(), AT_, &density0);
    MFloat dynamicViscosity0 = 0.00001716 * pow(T0 / 273.15, 1.5) * (273.15 + 110.4) / (T0 + 110.4);
    dynamicViscosity0 =
        Context::getSolverProperty<MFloat>("ambientDynViscosity", lpt().solverId(), AT_, &dynamicViscosity0);

    if(density0 < 0) {
      density0 = fvSolver().sysEqn().m_Re0 * dynamicViscosity0
                 / sqrt(fvSolver().sysEqn().gamma_Ref() * specificGasConstant * T0);
    }

    const MFloat pressure = density0 * specificGasConstant * Tinfty;

    velocityFactor = sqrt(fvSolver().sysEqn().gamma_Ref() * specificGasConstant * T0);
    densityFactor = density0;
    pressureFactor = pressure * fvSolver().sysEqn().gamma_Ref();
    temperatureFactor = T0;
    viscosityFactor = dynamicViscosity0;

    lpt().m_material->m_temperatureFactor = temperatureFactor;
    lpt().m_material->m_viscosityFactor = viscosityFactor;
  }

  //---------------------------------------------------------------------------------------

  conversionFvLpt.velocity = velocityFactor;
  conversionLptFv.velocity = 1.0 / conversionFvLpt.velocity;
  conversionFvLpt.density = densityFactor;
  conversionLptFv.density = 1.0 / conversionFvLpt.density;
  conversionFvLpt.pressure = pressureFactor;
  conversionLptFv.pressure = 1.0 / conversionFvLpt.pressure;
  conversionLptFv.length = lengthFactor;
  conversionFvLpt.length = 1.0 / conversionLptFv.length;
  conversionFvLpt.temperature = temperatureFactor;
  conversionLptFv.temperature = 1 / conversionFvLpt.temperature;
  conversionFvLpt.viscosity = viscosityFactor;
  conversionLptFv.viscosity = 1 / conversionFvLpt.viscosity;

  ASSERT(fabs(lpt().m_sutherlandConstant - fvSolver().m_sutherlandConstant) < MFloatEps, "");
  ASSERT(fabs(lpt().m_sutherlandPlusOne - fvSolver().m_sutherlandPlusOne) < MFloatEps, "");

  if(lpt().domainId() == 0) {
    cerr << "Fv-Lpt conversion factors are: " << endl;
    cerr << "length      : " << conversionFvLpt.length << endl;
    cerr << "velocity    : " << conversionFvLpt.velocity << endl;
    cerr << "density     : " << conversionFvLpt.density << endl;
    cerr << "pressure    : " << conversionFvLpt.pressure << endl;
    cerr << "temperature : " << conversionFvLpt.temperature << endl;
    cerr << "viscosity   : " << conversionFvLpt.viscosity << endl;
  }

  conversionFvLpt.mass = conversionFvLpt.velocity * conversionFvLpt.density * POW2(conversionFvLpt.length);
  conversionLptFv.mass = 1 / conversionFvLpt.mass;
  conversionFvLpt.momentum = conversionFvLpt.density * POW2(conversionFvLpt.velocity) * POW2(conversionFvLpt.length);
  conversionLptFv.momentum = 1 / conversionFvLpt.momentum;
  conversionFvLpt.energy = conversionFvLpt.density * conversionFvLpt.velocity * POW2(conversionFvLpt.velocity)
                           * POW2(conversionFvLpt.length);
  conversionLptFv.energy = 1 / conversionFvLpt.energy;

  if(lpt().m_ellipsoids) {
    conversionFvLpt.velocitySlope = conversionFvLpt.velocity / conversionFvLpt.length;
    conversionLptFv.velocitySlope = 1 / conversionFvLpt.velocitySlope;
  }
}

/** \brief Sets the initial particle velocity based on the flow field velocity in that cell.
 *         NOTE: This can not be done in the LPT initialCondition,
 *         as the FV flow field has not been transfered before, its set in the FV initialCondition.
 *  \author Tim Wegmann
 */
template <MInt nDim, class SysEqn>
void CouplerFvParticle<nDim, SysEqn>::initParticleVelocity() {
  if(!lpt().m_restart) {
    // init fluid density and velocity, after the values have been transferred to the LPT solver!
    for(MInt i = 0; i < lpt().a_noParticles(); i++) {
      lpt().m_partList[i].initVelocityAndDensity();
    }
    // init fluid density and velocity for ellipsoidal particles
    for(MInt i = 0; i < lpt().a_noEllipsoidalParticles(); i++) {
      lpt().m_partListEllipsoid[i].initVelocityAndDensity();
    }
  }
}

/** \brief Initialize coupling-class-specific Data
 *    \author Tim Wegmann
 */
template <MInt nDim, class SysEqn>
void CouplerFvParticle<nDim, SysEqn>::initData() {
  TRACE();

  m_fvSolverId = fvSolver().m_solverId;
  m_lptSolverId = lpt().m_solverId;

  readProperties();
}


template <MInt nDim, class SysEqn>
void CouplerFvParticle<nDim, SysEqn>::finalizeCouplerInit() {
  TRACE();

  // Transfer the flow field after the intial condition for FV has been applied
  transferCellStatus();
  transferFlowField();
  transferVelocitySlopes();
  if(lpt().isActive()) {
    lpt().receiveFlowField();
    lpt().receiveVelocitySlopes();
    lpt().waitForSendReqs();
  }

  initParticleVelocity();

  if(!lpt().m_restartFile && fvSolver().m_restartFile) {
    if(lpt().domainId() == 0) {
      cerr << "Restart without LPT restartFile -> setting time from Fv-time!" << endl;
    }
    const MFloat conversionTime = conversionFvLpt.length / conversionFvLpt.velocity;
    const MFloat fvTime = fvSolver().m_levelSetMb ? fvSolver().m_physicalTime : fvSolver().m_time;
    lpt().m_time = fvTime * conversionTime;
  }

  // set the timeStep in the LPT solver
  if(m_forceFvTimeStep) {
    transferTimeStep();
  } else {
    unifyTimeStep();
  }

  MBool writeParticleStats = false;
  if(writeParticleStats) {
    writeParticleCellStats();
  }
}

/** \brief read coupler data
 *  \author Tim Wegmann
 */
template <MInt nDim, class SysEqn>
void CouplerFvParticle<nDim, SysEqn>::readProperties() {
  TRACE();

  m_lptFvInterpolation =
      Context::getSolverProperty<MBool>("allowLPTFvInterpolation", m_fvSolverId, AT_, &m_lptFvInterpolation);

  /*! \page propertiesCoupling
      \section fvLptSpeciesId
    <code>MFloat COUPLER::m_fvLPTSpeciesId </code>\n
    default = <code>0</code>\n \n
    Specify the speciesId which the LPT particle liquid phase has in the FV-solver
    <ul>
    positive integer values
    </ul>
    Keywords: <i> LPT, FV, SPECIES, COUPLING </i>
  */
  m_fvLPTSpeciesId = 0;
  m_fvLPTSpeciesId = Context::getSolverProperty<MInt>("fvLptSpeciesId", m_lptSolverId, AT_, &m_fvLPTSpeciesId);

  /*! \page propertiesCoupling
      \section forceFvTimeStep
    <code>MFloat COUPLER::m_forceFvTimeStep </code>\n
    default = <code>true</code>\n \n
    Specify whether the FV-solver alone should force the time-step in the FV-LPT coupling!
    Otherwise, the minimum timestep of LPT and FV solver is computed and enforce!
    This however requires a sub-coupling!
    <ul>
    true,false
    </ul>
    Keywords: <i> LPT, FV, COUPLING </i>
  */
  m_forceFvTimeStep = true;
  m_forceFvTimeStep = Context::getSolverProperty<MBool>("forceFvTimeStep", m_lptSolverId, AT_, &m_forceFvTimeStep);
  if(!m_forceFvTimeStep && !m_interLeafed) {
    mTerm(1, AT_, "Coupled timesteps only working with interleafed setup!");
  }
}

/** \brief Checks property-data which is read in by both lpt-and Fv-Solver
 *    \author Tim Wegmann
 */
template <MInt nDim, class SysEqn>
void CouplerFvParticle<nDim, SysEqn>::checkProperties() {
  TRACE();

  ASSERT(fabs(fvSolver().m_sutherlandConstant - lpt().m_sutherlandConstant) < MFloatEps, "");
  ASSERT(fabs(lpt().m_sutherlandPlusOne - fvSolver().m_sutherlandPlusOne) < MFloatEps, "");
  ASSERT(lpt().a_timeStepComputationInterval() == fvSolver().a_timeStepComputationInterval(), "");
}

/** \brief preCoupler: reset external source terms before the LPT timestep
 *  \author Tim Wegmann
 */
template <MInt nDim, class SysEqn>
void CouplerFvParticle<nDim, SysEqn>::preCouple(MInt recipeStep) {
  TRACE();

  // PRE LPT
  if(recipeStep == m_lptSolverOrder) {
    m_solutionStep = 0;

    // TODO-timw labels:COUPLER,FV,toenhance,totest check if this can be avoided when transfering variables at some
    // other stage!
    //           when removed, testcase results are changed slightly!
    if(!m_interLeafed) fvSolver().computeConservativeVariables();

    // update LPT bndryCells for wall-collision
    updateLPTBndry();
    transferCellStatus();


    if(m_interLeafed) {
      const MFloat conversionTime = conversionFvLpt.length / conversionFvLpt.velocity;
      const MFloat fvTime = fvSolver().m_levelSetMb ? fvSolver().m_physicalTime : fvSolver().m_time;
      if(lpt().isActive() && fabs(fvTime * conversionTime - lpt().m_time) > MFloatEps) {
        cerr << "LPT-time " << lpt().m_time << " Fv-time " << fvTime << endl;
        mTerm(1, AT_, "Miss-matching time!");
      }
    }
  }
}

/** \brief postCouple: exchange source terms after the LPT timeStep
 *  \author Tim Wegmann
 */
template <MInt nDim, class SysEqn>
void CouplerFvParticle<nDim, SysEqn>::postCouple(MInt recipeStep) {
  // when using inter-leafed application, external sources and flow field are transfered in the subCouple!
  if(m_interLeafed) {
    if(recipeStep == m_fvSolverOrder) {
      transferVelocitySlopes(); // transfer only after postTimestep
    }
    return;
  }

  // Post LPT (as the LPT-solver is called before the Fv-solver!)
  if(recipeStep == m_lptSolverOrder) {
    // reset and transfer external Sources to FV-solver
    if(fvSolver().m_hasExternalSource) {
      fvSolver().applyExternalOldSource();
      fvSolver().resetExternalSources();
      transferExternalSources();
    }
  }

  // POST FV
  if(recipeStep == m_fvSolverOrder) {
    transferFlowField();
    transferVelocitySlopes();
    if(lpt().isActive()) {
      lpt().receiveFlowField();
      lpt().receiveVelocitySlopes();
      lpt().waitForSendReqs();
    }
    transferTimeStep();
  }
}

/** \brief postCouple: exchange source terms after the LPT timeStep
 *  \author Tim Wegmann
 */
template <MInt nDim, class SysEqn>
void CouplerFvParticle<nDim, SysEqn>::finalizeBalance(MInt id) {
  // post LPT solver
  if(id == lpt().solverId()) {
    if(fvSolver().m_hasExternalSource) {
      fvSolver().resetExternalSources();
    }
    if(globalTimeStep < 0) {
      prepareAdaptation();
    } else {
      transferExternalSources();
      fvSolver().advanceExternalSource();
    }

  } else if(id == fvSolver().solverId()) {
    transferCellStatus();
    transferFlowField();
    transferVelocitySlopes();
    if(lpt().isActive()) {
      lpt().receiveFlowField();
      lpt().receiveVelocitySlopes();
      lpt().waitForSendReqs();
    }

    updateLPTBndry();
  }
}

/** \brief transfer flow data from FV to LPT
 *
 *  NOTE: each LPT-cell must have at-least a matching FV-parent!
 *
 *  \author Julian Vorspohl, Tim Wegmann
 */
template <MInt nDim, class SysEqn>
void CouplerFvParticle<nDim, SysEqn>::transferFlowField() {
  TRACE();

#ifndef NDEBUG
  MBool diverged = false;
#endif

  for(MInt lptCellId = 0; lptCellId < lpt().noInternalCells(); lptCellId++) {
    if(!lpt().c_isLeafCell(lptCellId)) continue;
    MInt fvCellId = lpt2fvId(lptCellId);
    if(fvCellId > -1 && fvSolver().c_isLeafCell(fvCellId)) {
      // simple transfer if both cells are leaf-cells
      for(MInt n = 0; n < nDim; n++) {
        lpt().a_fluidVelocity(lptCellId, n) =
            fvSolver().a_pvariable(fvCellId, fvSolver().m_sysEqn.PV->VV[n]) * conversionFvLpt.velocity;
      }
      lpt().a_fluidDensity(lptCellId) =
          fvSolver().a_pvariable(fvCellId, fvSolver().m_sysEqn.PV->RHO) * conversionFvLpt.density;
      lpt().a_fluidPressure(lptCellId) =
          fvSolver().a_pvariable(fvCellId, fvSolver().m_sysEqn.PV->P) * conversionFvLpt.pressure;
#ifndef NDEBUG
      if(std::isnan(lpt().a_fluidPressure(lptCellId)) && lpt().a_isValidCell(lptCellId)) {
        diverged = true;
        cerr << "Nan pressure " << lpt().a_fluidPressure(lptCellId) << " " << lpt().a_coordinate(lptCellId, 0) << " "
             << lpt().a_coordinate(lptCellId, 1) << " " << lpt().a_coordinate(lptCellId, nDim - 1) << endl;
        // mTerm(1, AT_, "Nan pressure detected!");
      }
#endif
      lpt().a_fluidTemperature(lptCellId) =
          fvSolver().sysEqn().temperature_ES(fvSolver().a_pvariable(fvCellId, fvSolver().m_sysEqn.PV->RHO),
                                             fvSolver().a_pvariable(fvCellId, fvSolver().m_sysEqn.PV->P))
          * conversionFvLpt.temperature;
      if(lpt().m_evaporation && fvSolver().m_noSpecies > 0) {
        lpt().a_fluidSpecies(lptCellId) =
            mMax(0.0, mMin(1.0, fvSolver().a_pvariable(fvCellId, fvSolver().m_sysEqn.PV->Y[m_fvLPTSpeciesId])));
      }

    } else if(fvCellId > -1 && !fvSolver().c_isLeafCell(fvCellId)) {
      // interpolate from FV-leaf childs if the matching fv-Cell is not a leaf cell
      if(lpt().a_isHalo(lptCellId)) continue;
      // gather all fv-leaf-cell childs for internal cells
      vector<MInt> leafChildIds;
      fvSolver().grid().getAllLeafChilds(fvCellId, leafChildIds);
      ASSERT(!leafChildIds.empty(), "");

      // reset variables
      for(MInt n = 0; n < nDim; n++) {
        lpt().a_fluidVelocity(lptCellId, n) = 0;
      }
      lpt().a_fluidDensity(lptCellId) = 0;
      lpt().a_fluidPressure(lptCellId) = 0;
      lpt().a_fluidTemperature(lptCellId) = 0;
      if(lpt().m_evaporation) {
        lpt().a_fluidSpecies(lptCellId) = 0;
      }
      // interpolate variables by volume
      MFloat volume = 0;
      for(MUint i = 0; i < leafChildIds.size(); i++) {
        const MInt cellId = leafChildIds[i];
        if(fvSolver().a_isInactive(cellId)) continue;
        const MFloat cellVolume = fvSolver().a_cellVolume(cellId);
        ASSERT(fvSolver().c_isLeafCell(cellId), "");
        // ASSERT(cellVolume > 0, "");
        if(cellVolume < 1e-16) {
          mTerm(1, AT_, "Invalid cell-volume!");
        }
        volume += cellVolume;
        for(MInt n = 0; n < nDim; n++) {
          lpt().a_fluidVelocity(lptCellId, n) +=
              cellVolume * fvSolver().a_pvariable(cellId, fvSolver().m_sysEqn.PV->VV[n]);
        }
        lpt().a_fluidDensity(lptCellId) += cellVolume * fvSolver().a_pvariable(cellId, fvSolver().m_sysEqn.PV->RHO);
        lpt().a_fluidPressure(lptCellId) += cellVolume * fvSolver().a_pvariable(cellId, fvSolver().m_sysEqn.PV->P);
        lpt().a_fluidTemperature(lptCellId) +=
            cellVolume
            * fvSolver().sysEqn().temperature_ES(fvSolver().a_pvariable(fvCellId, fvSolver().m_sysEqn.PV->RHO),
                                                 fvSolver().a_pvariable(fvCellId, fvSolver().m_sysEqn.PV->P));
        if(lpt().m_evaporation && fvSolver().m_noSpecies > 0) {
          lpt().a_fluidSpecies(lptCellId) += mMax(
              0.0,
              mMin(1.0, cellVolume * fvSolver().a_pvariable(fvCellId, fvSolver().m_sysEqn.PV->Y[m_fvLPTSpeciesId])));
        }
      }

      if(volume > 1e-16) {
        // meaning at least one child is active
        for(MInt n = 0; n < nDim; n++) {
          lpt().a_fluidVelocity(lptCellId, n) *= conversionFvLpt.velocity / volume;
        }
        lpt().a_fluidDensity(lptCellId) *= conversionFvLpt.density / volume;
        lpt().a_fluidPressure(lptCellId) *= conversionFvLpt.pressure / volume;
        lpt().a_fluidTemperature(lptCellId) *= conversionFvLpt.temperature / volume;
        if(lpt().m_evaporation && fvSolver().m_noSpecies > 0) {
          lpt().a_fluidSpecies(lptCellId) *= 1 / volume;
        }
      } else {
        // meaning completely outside
        for(MInt n = 0; n < nDim; n++) {
          lpt().a_fluidVelocity(lptCellId, n) = conversionFvLpt.velocity;
        }
        lpt().a_fluidDensity(lptCellId) = conversionFvLpt.density;
        lpt().a_fluidPressure(lptCellId) = conversionFvLpt.pressure;
        lpt().a_fluidTemperature(lptCellId) = conversionFvLpt.temperature;
        if(lpt().m_evaporation && fvSolver().m_noSpecies > 0) {
          lpt().a_fluidSpecies(lptCellId) = 1;
        }
      }
    } else if(fvCellId < 0) {
      if(lpt().a_isHalo(lptCellId)) continue;
      // interpolate linear from the FV parents to the LPT-child
      fvCellId = lpt2fvIdParent(lptCellId);
      // ASSERT(fvCellId > -1, "LPT-cell has no matching FV-cell!");
      if(fvCellId < 0) {
        cerr << "LPT-cell is missing fv-cell parent" << lptCellId << " " << fvCellId << " " << lpt().a_level(lptCellId)
             << endl;
      }
      interpolateFVLPT(fvCellId, lptCellId);
    }
  }

  // NOTE: exchange is necessary as particles on window-cells interpolate the flow variables to their
  //       position and will need values on halo-cells to do so!
  if(!lpt().m_nonBlockingComm) {
    lpt().exchangeData(&lpt().a_fluidVariable(0, 0), lpt().PV.noVars());
    if(lpt().m_evaporation) {
      lpt().exchangeData(&lpt().a_fluidSpecies(0), 1);
    }
#ifndef NDEBUG
    if(lpt().isActive()) {
      MPI_Allreduce(MPI_IN_PLACE, &diverged, 1, MPI_C_BOOL, MPI_LOR, lpt().mpiComm(), AT_, "MPI_IN_PLACE", "divCheck");
      if(diverged) {
        lpt().saveDebugRestartFile();
      }
    }
    if(fvSolver().isActive()) {
      MPI_Allreduce(MPI_IN_PLACE, &diverged, 1, MPI_C_BOOL, MPI_LOR, fvSolver().mpiComm(), AT_, "MPI_IN_PLACE",
                    "divCheck");
      if(diverged) {
        fvSolver().writeVtkXmlFiles("QOUT", "GEOM", false, true);
      }
    }
#endif
  } else {
    if(lpt().isActive()) {
      lpt().sendFlowField();
      // the flow field data is received in the time-step
    }
  }
}


/** \brief transfer all relevant bndryCell-data from FV to LPT solver before the LPT timeStep!
 *  \author Tim Wegmann
 */
template <MInt nDim, class SysEqn>
void CouplerFvParticle<nDim, SysEqn>::updateLPTBndry() {
  TRACE();

  if(!lpt().m_wallCollisions) return;

  if(!lpt().isActive()) return;

  // reset all moving bndry-Cells
  for(MInt lptCell = 0; lptCell < lpt().a_noCells(); lptCell++) {
    if(lpt().a_bndryCellId(lptCell) >= lpt().m_noOuterBndryCells) {
      lpt().a_bndryCellId(lptCell) = -1;
    }
  }

  // reset size to stationary outer bndryCells
  lpt().m_bndryCells->resetSize(lpt().m_noOuterBndryCells);

  MFloatScratchSpace surfaceC(nDim, nDim + 1, AT_, "surfaceC");
  MFloatScratchSpace surfaceN(nDim, AT_, "surfaceN");
  MFloatScratchSpace surfaceV(nDim, AT_, "surfaceV");
  MFloatScratchSpace bndryData(nDim + 1, AT_, "surfaceV");

  // update moving bndryCells from FV
  for(MInt bndryId = fvSolver().m_noOuterBndryCells; bndryId < fvSolver().m_fvBndryCnd->m_bndryCells->size();
      bndryId++) {
    const MInt fvCellId = fvSolver().m_bndryCells->a[bndryId].m_cellId;
    // TODO: use splitChildToSplitcell for bndry-cells of sliptchilds!
    if(fvCellId < 0 || fvCellId >= fvSolver().c_noCells()) continue;
    if(!fvSolver().c_isLeafCell(fvCellId)) continue;
    MInt lptCellId = fv2lptId(fvCellId);
    if(lptCellId < 0) {
      lptCellId = fv2lptIdParent(fvCellId);
      if(lptCellId < 0) continue;
      if(lpt().c_level(lptCellId) < lpt().maxRefinementLevel()) continue;
      ASSERT(m_lptFvInterpolation, "");
    }

    if(!lpt().c_isLeafCell(lptCellId)) {
      ASSERT(m_lptFvInterpolation, "");
    }

    MInt noSurf = fvSolver().m_bndryCells->a[bndryId].m_noSrfcs;

    // bndryCell without a bndrySurface, how is that even possible?
    // if(noSurf < 1 ) continue;

    // ASSERT(noSurf == 1, "Only for single bndry-surfaces!");
    if(noSurf > 1) noSurf = 1;

    bndryData[0] = fvSolver().m_fvBndryCnd->m_bndryCells->a[bndryId].m_volume;
    for(MInt i = 0; i < nDim; i++) {
      bndryData[1 + i] = fvSolver().m_fvBndryCnd->m_bndryCells->a[bndryId].m_coordinates[i];
    }

    for(MInt s = 0; s < noSurf; s++) {
      MInt bodyId = fvSolver().m_levelSetMb ? fvSolver().m_bndryCells->a[bndryId].m_srfcs[s]->m_bodyId[0] : -1;
      for(MInt n = 0; n < nDim; n++) {
        surfaceC(noSurf * s + n, 0) = fvSolver().m_fvBndryCnd->m_bndryCells->a[bndryId].m_srfcs[s]->m_coordinates[n];
        for(MInt i = 0; i < nDim; i++) {
          surfaceC(noSurf * s + n, i + 1) =
              fvSolver().m_fvBndryCnd->m_bndryCells->a[bndryId].m_srfcs[s]->m_cutCoordinates[i][n];
        }
        // surfaceV[noSurf * s + n] = fvSolver()
        //                               .m_fvBndryCnd->m_bndryCells->a[bndryId]
        //                               .m_srfcVariables[s]
        //                               ->m_primVars[fvSolver().m_sysEqn.PV->VV[n]];
        if(bodyId > 0) {
          surfaceV[noSurf * s + n] = fvSolver().m_bodyVelocity[bodyId * nDim + n];
        } else {
          surfaceV[noSurf * s + n] = 0.0;
        }
        surfaceN[noSurf * s + n] = fvSolver().m_fvBndryCnd->m_bndryCells->a[bndryId].m_srfcs[s]->m_normalVector[n];
      }
    }

    // append new LPT bndryCell based on FV data
    lpt().addBndryCell(lptCellId, &bndryData[0], &surfaceN[0], surfaceC, &surfaceV[0]);
  }
}

/** \brief set a_noPart in Fv based on a_noParticlesInCell and a_noEllipsoidsInCell in LPT solver
 *         NOTE: this is only necessary if the particle sensor is used during adaptation
 *  \author Tim Wegmann
 */
template <MInt nDim, class SysEqn>
void CouplerFvParticle<nDim, SysEqn>::transferNoParticlesInCell() {
  TRACE();

  if(!fvSolver().m_sensorParticle) return;

  // reset
  for(MInt fvCell = 0; fvCell < fvSolver().c_noCells(); fvCell++) {
    fvSolver().a_noPart(fvCell) = 0;
  }
  for(MInt fvCell = fvSolver().c_noCells(); fvCell < fvSolver().a_noCells(); fvCell++) {
#ifndef NDEBUG
    if(fvSolver().a_noPart(fvCell) > 0) {
      cerr << " Particles in " << fvCell << " " << fvSolver().a_noCells() << " " << fvSolver().a_isHalo(fvCell)
           << fvSolver().a_isBndryGhostCell(fvCell) << endl;
    }
#endif
    fvSolver().a_noPart(fvCell) = 0;
  }

  // transfer
  for(MInt lptCell = 0; lptCell < lpt().noInternalCells(); lptCell++) {
    if(lpt().a_noParticlesInCell(lptCell) < 1 && lpt().a_noEllipsoidsInCell(lptCell) < 1) continue;
    MInt fvCellId = lpt2fvId(lptCell);
    if(fvCellId < 0) {
      ASSERT(m_lptFvInterpolation, "");
      fvCellId = lpt2fvIdParent(lptCell);
      if(fvCellId < 0) continue;
    }

    if(!fvSolver().c_isLeafCell(fvCellId)) {
      ASSERT(m_lptFvInterpolation, "");
    }

    fvSolver().a_noPart(fvCellId) = lpt().a_noParticlesInCell(lptCell) + lpt().a_noEllipsoidsInCell(lptCell);
  }

  fvSolver().exchangeData(&fvSolver().a_noPart(0));
}


/** \brief set external sources in Fv solver based on values in the LPT solver
 *  \author Tim Wegmann
 */
template <MInt nDim, class SysEqn>
void CouplerFvParticle<nDim, SysEqn>::transferExternalSources() {
  TRACE();

  if(!fvSolver().isActive()) return;

  if(!lpt().m_momentumCoupling && !lpt().m_heatCoupling && !lpt().m_massCoupling) return;

  if(lpt().isActive()) lpt().exchangeSourceTerms();

#ifndef NDEBUG
  // check that externalSource is set to zero!
  for(MInt cellId = 0; cellId < fvSolver().a_noCells(); cellId++) {
    for(MInt var = 0; var < fvSolver().CV->noVariables; var++) {
      // NOTE: only valid for single particle coupling!
      // ASSERT(abs(fvSolver().a_externalSource(cellId, var)) < MFloatEps, "");
      if(abs(fvSolver().a_externalSource(cellId, var)) > MFloatEps) {
        mTerm(1, AT_, "External Source not reset correctly!!");
      }
    }
  }
#endif

  for(MInt lptCell = 0; lptCell < lpt().noInternalCells(); lptCell++) {
    // the source term is only set on LPT leaf-cells
    if(!lpt().c_isLeafCell(lptCell)) continue;
    MInt fvCellId = lpt2fvId(lptCell);
    if(fvCellId > -1 && fvSolver().c_isLeafCell(fvCellId)) {
      // simple transfer
      setExternalSourceInCell(fvCellId, lptCell, -1);

    } else if(fvCellId > -1 && !fvSolver().c_isLeafCell(fvCellId)) {
      // fv-cell is refined further: volume average from LPT to FV
      // gather all fv-leaf-cell children for internal cells
      vector<MInt> leafChildIds;
      fvSolver().grid().getAllLeafChilds(fvCellId, leafChildIds);
      MFloat totalVol = 0;
      for(MUint i = 0; i < leafChildIds.size(); i++) {
        const MInt cellId = leafChildIds[i];
        totalVol += fvSolver().a_cellVolume(cellId);
      }
      for(MUint i = 0; i < leafChildIds.size(); i++) {
        const MInt cellId = leafChildIds[i];
        const MFloat vFrac = fvSolver().a_cellVolume(cellId) / totalVol;
        setExternalSourceInCell(cellId, lptCell, vFrac);
      }
    } else {
      // lpt-Cell is refined further
      fvCellId = lpt2fvIdParent(lptCell);
      if(fvCellId < 0) continue;
      setExternalSourceInCell(fvCellId, lptCell, -1);
    }
  }

  // count source-terms in the fv-solver:
  MInt numVars = 2 + nDim;
  MFloatScratchSpace conservativeSums(numVars, AT_, "conservativeSums");
  conservativeSums.fill(0.0);
  const MFloat dt = fvSolver().timeStep();
  for(MInt cellId = 0; cellId < fvSolver().noInternalCells(); cellId++) {
    if(!fvSolver().c_isLeafCell(cellId)) continue;
    if(fvSolver().a_isInactive(cellId)) continue;

    if(fvSolver().m_noSpecies > 0) {
      conservativeSums[0] -= dt * fvSolver().a_externalSource(cellId, fvSolver().CV->RHO_Y[m_fvLPTSpeciesId]);
    }
    conservativeSums[1] -= dt * fvSolver().a_externalSource(cellId, fvSolver().CV->RHO_E);
    for(MInt i = 0; i < nDim; i++) {
      conservativeSums[2 + i] -= dt * fvSolver().a_externalSource(cellId, fvSolver().CV->RHO_VV[i]);
    }
  }

  if(fvSolver().m_vapourData.size() > 0) {
    auto it = fvSolver().m_vapourData.find(globalTimeStep - 1);
    if(it != fvSolver().m_vapourData.end()) {
      for(MInt i = 0; i < numVars; i++) {
        MFloat lastSource = (it->second)[i];
        conservativeSums[i] += lastSource;
      }
    }
  }

  vector<MFloat> tmp(numVars);
  for(MInt j = 0; j < numVars; j++) {
    tmp[j] = conservativeSums[j];
  }

  fvSolver().m_vapourData.insert(make_pair(globalTimeStep, tmp));
}


/** \brief transfer/enforce Fv timeStep onto LPT solver
 *  \author Tim Wegmann
 */
template <MInt nDim, class SysEqn>
void CouplerFvParticle<nDim, SysEqn>::transferTimeStep() {
  TRACE();

  ASSERT(m_forceFvTimeStep, "");

  // set the timeStep in the LPT solver
  const MFloat conversionTime = conversionFvLpt.length / conversionFvLpt.velocity;
  lpt().forceTimeStep(fvSolver().timeStep(true) * conversionTime);
  const MFloat fvTime = fvSolver().m_levelSetMb ? fvSolver().m_physicalTime : fvSolver().m_time;
  if(lpt().isActive()) {
    if(fabs(fvTime * conversionTime - lpt().m_time) > MFloatEps) {
      cerr << "LPT-time " << lpt().m_time << " Fv-time " << fvTime << " " << conversionTime << " " << globalTimeStep
           << " " << lpt().isActive() << endl;
      mTerm(1, AT_, "Miss-matching time!");
    }
  } else {
    lpt().forceTimeStep(fvSolver().timeStep(true) * fvSolver().m_timeRef);

    cerr0 << "Fv-LPT time step set as " << fvSolver().timeStep(true) * fvSolver().m_timeRef << endl;
  }

  // increasing CFL number after the injection!
  if(lpt().m_sprayModel != nullptr && lpt().m_sprayModel->m_injectionCrankAngle > 0) {
    const MFloat cad = fvSolver().crankAngle(fvTime, 0);
    MFloat injDuration = 10;
    injDuration = Context::getSolverProperty<MFloat>("injectorInjectionTime", lpt().solverId(), AT_, &injDuration);
    MInt injCAD = 40;
    if(injDuration > 15) {
      injCAD = 55;
    }
    if(cad > lpt().m_sprayModel->m_injectionCrankAngle + injCAD && fvSolver().m_cfl < 0.7) {
      cerr0 << "Increasing FV CFL for next computation from " << fvSolver().m_cfl << " to 0.95!" << endl;
      fvSolver().m_cfl = 0.95;
    }
  }
}

/** \brief prepate adaptation
 *  \author Tim Wegmann
 */
template <MInt nDim, class SysEqn>
void CouplerFvParticle<nDim, SysEqn>::prepareAdaptation() {
  TRACE();

  // transfer the number of particles for the particle limit sensor!
  transferNoParticlesInCell();
}

/** \brief set external source in FV-solver from LPT fluxes
 *  \author Tim Wegmann
 */
template <MInt nDim, class SysEqn>
void CouplerFvParticle<nDim, SysEqn>::setExternalSourceInCell(const MInt fvCellId,
                                                              const MInt lptCell,
                                                              const MFloat factor) {
  TRACE();

  if(factor < 0) {
    if(lpt().m_momentumCoupling) {
      for(MInt i = 0; i < nDim; i++) {
        fvSolver().a_externalSource(fvCellId, fvSolver().CV->RHO_VV[i]) +=
            lpt().a_momentumFlux(lptCell, i) / conversionFvLpt.momentum;
      }
      fvSolver().a_externalSource(fvCellId, fvSolver().CV->RHO_E) += lpt().a_workFlux(lptCell) / conversionFvLpt.energy;
    }
    if(lpt().m_heatCoupling) {
      fvSolver().a_externalSource(fvCellId, fvSolver().CV->RHO_E) += lpt().a_heatFlux(lptCell) / conversionFvLpt.energy;
    }
    if(lpt().m_massCoupling) {
      fvSolver().a_externalSource(fvCellId, fvSolver().CV->RHO) += lpt().a_massFlux(lptCell) / conversionFvLpt.mass;
      fvSolver().a_externalSource(fvCellId, fvSolver().CV->RHO_Y[m_fvLPTSpeciesId]) +=
          lpt().a_massFlux(lptCell) / conversionFvLpt.mass;
    }
  } else {
    if(lpt().m_momentumCoupling) {
      for(MInt i = 0; i < nDim; i++) {
        fvSolver().a_externalSource(fvCellId, fvSolver().CV->RHO_VV[i]) +=
            lpt().a_momentumFlux(lptCell, i) * factor / conversionFvLpt.momentum;
      }
      fvSolver().a_externalSource(fvCellId, fvSolver().CV->RHO_E) +=
          lpt().a_workFlux(lptCell) * factor / conversionFvLpt.energy;
    }
    if(lpt().m_heatCoupling) {
      fvSolver().a_externalSource(fvCellId, fvSolver().CV->RHO_E) +=
          lpt().a_heatFlux(lptCell) * factor / conversionFvLpt.energy;
    }
    if(lpt().m_massCoupling) {
      fvSolver().a_externalSource(fvCellId, fvSolver().CV->RHO) +=
          lpt().a_massFlux(lptCell) * factor / conversionFvLpt.mass;

      fvSolver().a_externalSource(fvCellId, fvSolver().CV->RHO_Y[m_fvLPTSpeciesId]) +=
          lpt().a_massFlux(lptCell) * factor / conversionFvLpt.mass;
    }
  }
}

/** \brief prepate adaptation
 *  \author Tim Wegmann
 */
template <MInt nDim, class SysEqn>
void CouplerFvParticle<nDim, SysEqn>::finalizeAdaptation(const MInt solverId) {
  TRACE();

  // Post FV solver
  if(solverId == fvSolver().solverId()) {
    transferCellStatus();
    transferFlowField();
    transferVelocitySlopes();
    if(lpt().isActive()) {
      lpt().receiveFlowField();
      lpt().receiveVelocitySlopes();
      lpt().waitForSendReqs();
    }

    updateLPTBndry();

  } else if(solverId == lpt().solverId()) {
    // set external sources in the FV-solver
    // as they are reset there in finalizeAdaptation
    if(fvSolver().m_hasExternalSource) {
      fvSolver().resetExternalSources();
      transferExternalSources();
    }
  }
}

/** \brief interpolate flow variables from the fv-grid to the LPT-cell position
 *
 * \author Tim Wegmann
 * \date March 2020
 */
template <MInt nDim, class SysEqn>
void CouplerFvParticle<nDim, SysEqn>::interpolateFVLPT(const MInt from, const MInt to) {
  TRACE();

  MInt interpolationCells[8] = {0, 0, 0, 0, 0, 0, 0, 0};
  MFloat point[3] = {lpt().c_coordinate(to, 0), lpt().c_coordinate(to, 1), lpt().c_coordinate(to, 2)};

  std::function<MBool(const MInt, const MInt)> neighborCheck = [&](const MInt cellId, const MInt id) {
    return static_cast<MBool>(fvSolver().checkNeighborActive(cellId, id));
  };

  const MInt position = fvSolver().setUpInterpolationStencil(from, interpolationCells, point, neighborCheck, true);

  if(position < 0) {
    for(MInt n = 0; n < nDim; n++) {
      lpt().a_fluidVelocity(to, n) =
          fvSolver().a_pvariable(from, fvSolver().m_sysEqn.PV->VV[n]) * conversionFvLpt.velocity;
    }
    lpt().a_fluidDensity(to) = fvSolver().a_pvariable(from, fvSolver().m_sysEqn.PV->RHO) * conversionFvLpt.density;
    lpt().a_fluidPressure(to) = fvSolver().a_pvariable(from, fvSolver().m_sysEqn.PV->P) * conversionFvLpt.pressure;
    lpt().a_fluidTemperature(to) =
        fvSolver().sysEqn().temperature_ES(fvSolver().a_pvariable(from, fvSolver().m_sysEqn.PV->RHO),
                                           fvSolver().a_pvariable(from, fvSolver().m_sysEqn.PV->P))
        * conversionFvLpt.temperature;
    if(lpt().m_evaporation && fvSolver().m_noSpecies > 0) {
      lpt().a_fluidSpecies(to) = fvSolver().a_pvariable(from, fvSolver().m_sysEqn.PV->Y[m_fvLPTSpeciesId]);
    }

  } else {
    for(MInt n = 0; n < nDim; n++) {
      lpt().a_fluidVelocity(to, n) =
          interpolateVariable(interpolationCells, point, fvSolver().m_sysEqn.PV->VV[n]) * conversionFvLpt.velocity;
    }
    const MFloat rho = interpolateVariable(interpolationCells, point, fvSolver().m_sysEqn.PV->RHO);
    const MFloat p = interpolateVariable(interpolationCells, point, fvSolver().m_sysEqn.PV->P);
    lpt().a_fluidDensity(to) = rho * conversionFvLpt.density;
    lpt().a_fluidPressure(to) = p * conversionFvLpt.pressure;
    lpt().a_fluidTemperature(to) = fvSolver().sysEqn().temperature_ES(rho, p) * conversionFvLpt.temperature;
    if(lpt().m_evaporation && fvSolver().m_noSpecies > 0) {
      lpt().a_fluidSpecies(to) =
          interpolateVariable(interpolationCells, point, fvSolver().m_sysEqn.PV->Y[m_fvLPTSpeciesId]);
    }
  }
}

/**
 * \brief Interpolates the velocity slopes from the fv-grid to the LPT-cell
 * \author Laurent Andre
 * \date August 2022
 */
template <MInt nDim, class SysEqn>
void CouplerFvParticle<nDim, SysEqn>::interpolateVelocitySlopesFVLPT(const MInt from, const MInt to) {
  TRACE();

  MInt interpolationCells[8] = {0, 0, 0, 0, 0, 0, 0, 0};
  MFloat point[3] = {lpt().c_coordinate(to, 0), lpt().c_coordinate(to, 1), lpt().c_coordinate(to, 2)};

  std::function<MBool(const MInt, const MInt)> neighborCheck = [&](const MInt cellId, const MInt id) {
    return static_cast<MBool>(fvSolver().checkNeighborActive(cellId, id));
  };

  const MInt position = fvSolver().setUpInterpolationStencil(from, interpolationCells, point, neighborCheck, true);

  if(position < 0) {
    for(MInt varId = 0; varId < nDim; varId++) {
      for(MInt dir = 0; dir < nDim; dir++) {
        lpt().a_velocitySlope(to, varId, dir) =
            fvSolver().a_slope(from, fvSolver().m_sysEqn.PV->VV[varId], dir) * conversionFvLpt.velocitySlope;
      }
    }
  } else {
    for(MInt varId = 0; varId < nDim; varId++) {
      for(MInt dir = 0; dir < nDim; dir++) {
        lpt().a_velocitySlope(to, varId, dir) =
            interpolateSlope(interpolationCells, point, fvSolver().m_sysEqn.PV->VV[varId], dir)
            * conversionFvLpt.velocitySlope;
      }
    }
  }
}

/** \brief interpolates the fv-variable
/// \author Tim Wegmann
/// \date 2020-04-01
*/
template <MInt nDim, class SysEqn>
MFloat CouplerFvParticle<nDim, SysEqn>::interpolateVariable(MInt* interpolationCells, MFloat* point, MInt v) {
  TRACE();

  std::function<MFloat(const MInt, const MInt)> scalarField = [&](const MInt cellId, const MInt varId) {
    return static_cast<MFloat>(fvSolver().a_pvariable(cellId, varId));
  };

  std::function<MFloat(const MInt, const MInt)> coordinate = [&](const MInt cellId, const MInt id) {
    return static_cast<MFloat>(fvSolver().a_coordinate(cellId, id));
  };

  return fvSolver().template interpolateFieldData<false>(&interpolationCells[0], &point[0], v, scalarField, coordinate);
}

/**
 * \brief interpolate fv slope
 * \author Laurent Andre
 * \date August 2022
 */
template <MInt nDim, class SysEqn>
MFloat CouplerFvParticle<nDim, SysEqn>::interpolateSlope(MInt* interpolationCells, MFloat* point, MInt v, MInt dir) {
  TRACE();

  std::function<MFloat(const MInt, const MInt)> scalarField = [&](const MInt cellId, const MInt varId) {
    return static_cast<MFloat>(*(&fvSolver().a_slope(cellId, 0, 0) + varId));
  };

  std::function<MFloat(const MInt, const MInt)> coordinate = [&](const MInt cellId, const MInt id) {
    return static_cast<MFloat>(fvSolver().a_coordinate(cellId, id));
  };

  MInt valIndex = v * nDim + dir;
  return fvSolver().template interpolateFieldData<false>(&interpolationCells[0], &point[0], valIndex, scalarField,
                                                         coordinate);
}

/** \brief set the isValid status for LPT cells based on the FV-solver cell properties
 *  \author Tim Wegmann
 */
template <MInt nDim, class SysEqn>
void CouplerFvParticle<nDim, SysEqn>::transferCellStatus() {
  for(MInt lptCellId = 0; lptCellId < lpt().a_noCells(); lptCellId++) {
    MInt fvCellId = lpt2fvIdParent(lptCellId);
    if(fvCellId < 0 || fvCellId > fvSolver().a_noCells()) {
      lpt().a_isValidCell(lptCellId) = false;
    } else if(fvSolver().a_hasProperty(fvCellId, FvCell::IsCutOff) && !lpt().grid().isPeriodic(lptCellId)) {
      lpt().a_isValidCell(lptCellId) = false;
    } else if(fvSolver().a_hasProperty(fvCellId, FvCell::IsInSpongeLayer)) {
      lpt().a_isValidCell(lptCellId) = false;
    } else if(fvSolver().a_isInactive(fvCellId)) {
      lpt().a_isValidCell(lptCellId) = false;
    } else {
      lpt().a_isValidCell(lptCellId) = true;
    }
  }
}

/** \brief transfer the FV velocity slopes to the LPT solver
 *
 *  NOTE: each LPT-cell must have at-least a matching FV-parent!
 *
 *  \author Tim Wegmann, update Laurent Andre
 *  \date August 2022 (edit)
 */
template <MInt nDim, class SysEqn>
void CouplerFvParticle<nDim, SysEqn>::transferVelocitySlopes() {
  if(!lpt().m_ellipsoids) return;

#ifndef NDEBUG
  MBool diverged = false;
#endif

  for(MInt lptCellId = 0; lptCellId < lpt().noInternalCells(); lptCellId++) {
    if(!lpt().c_isLeafCell(lptCellId)) continue;
    MInt fvCellId = lpt2fvId(lptCellId);
    if(fvCellId > -1 && fvSolver().c_isLeafCell(fvCellId)) {
      // simple transfer if both cells are leaf-cells
      for(MInt varId = 0; varId < nDim; varId++) {
        for(MInt dir = 0; dir < nDim; dir++) {
          lpt().a_velocitySlope(lptCellId, varId, dir) =
              fvSolver().a_slope(fvCellId, fvSolver().m_sysEqn.PV->VV[varId], dir) * conversionFvLpt.velocitySlope;
#ifndef NDEBUG
          if(std::isnan(lpt().a_velocitySlope(lptCellId, varId, dir)) && lpt().a_isValidCell(lptCellId)
             && globalTimeStep > 0) {
            diverged = true;
            cerr << "Nan velocity slope: " << lpt().a_velocitySlope(lptCellId, varId, dir) << " "
                 << lpt().a_coordinate(lptCellId, 0) << " " << lpt().a_coordinate(lptCellId, 1) << " "
                 << lpt().a_coordinate(lptCellId, nDim - 1) << endl;
            mTerm(1, AT_, "FVParticle: Transfer of slope that is NAN");
          }
#endif
        }
      }
    } else if(fvCellId > -1 && !fvSolver().c_isLeafCell(fvCellId)) {
      // interpolate from FV-leaf childs if the matching fv-Cell is not a leaf cell
      if(lpt().a_isHalo(lptCellId)) continue;
      // gather all fv-leaf-cell childs for internal cells
      vector<MInt> leafChildIds;
      fvSolver().grid().getAllLeafChilds(fvCellId, leafChildIds);
      ASSERT(!leafChildIds.empty(), "");

      // reset variables
      for(MInt varId = 0; varId < nDim; varId++) {
        for(MInt dir = 0; dir < nDim; dir++) {
          lpt().a_velocitySlope(lptCellId, varId, dir) = F0;
        }
      }
      // interpolate variables by volume
      MFloat volume = 0;
      for(MUint i = 0; i < leafChildIds.size(); i++) {
        const MInt cellId = leafChildIds[i];
        if(fvSolver().a_isInactive(cellId)) continue;
        const MFloat cellVolume = fvSolver().a_cellVolume(cellId);
        ASSERT(fvSolver().c_isLeafCell(cellId), "");
        if(cellVolume < 1e-16) mTerm(1, AT_, "Invalid cell-volume!");
        volume += cellVolume;
        for(MInt varId = 0; varId < nDim; varId++) {
          for(MInt dir = 0; dir < nDim; dir++) {
            lpt().a_velocitySlope(lptCellId, varId, dir) +=
                cellVolume * fvSolver().a_slope(fvCellId, fvSolver().m_sysEqn.PV->VV[varId], dir);
          }
        }
      }
      if(volume > 1e-16) {
        // meaning at least one child is active
        for(MInt varId = 0; varId < nDim; varId++) {
          for(MInt dir = 0; dir < nDim; dir++) {
            lpt().a_velocitySlope(lptCellId, varId, dir) *= conversionFvLpt.velocitySlope / volume;
          }
        }
      } else {
        // meaning completely outside
        for(MInt varId = 0; varId < nDim; varId++) {
          for(MInt dir = 0; dir < nDim; dir++) {
            lpt().a_velocitySlope(lptCellId, varId, dir) = conversionFvLpt.velocitySlope;
          }
        }
      }
    } else if(fvCellId < 0) {
      if(lpt().a_isHalo(lptCellId)) continue;
      // interpolate linear from the FV parents to the LPT-child
      fvCellId = lpt2fvIdParent(lptCellId);
      ASSERT(fvCellId > -1, "LPT-cell has no matching FV-cell!");
      interpolateVelocitySlopesFVLPT(fvCellId, lptCellId);
    }
  }

  // NOTE: exchange is necessary as particles on window-cells interpolate the flow variables to their
  //       position and will need values on halo-cells to do so!
  if(!lpt().m_nonBlockingComm) {
    lpt().exchangeData(&lpt().a_velocitySlope(0, 0, 0), nDim * nDim);
#ifndef NDEBUG
    if(lpt().isActive()) {
      MPI_Allreduce(MPI_IN_PLACE, &diverged, 1, MPI_C_BOOL, MPI_LOR, lpt().mpiComm(), AT_, "MPI_IN_PLACE", "divCheck");
      if(diverged) {
        lpt().saveDebugRestartFile();
      }
    }
    if(fvSolver().isActive()) {
      MPI_Allreduce(MPI_IN_PLACE, &diverged, 1, MPI_C_BOOL, MPI_LOR, fvSolver().mpiComm(), AT_, "MPI_IN_PLACE",
                    "divCheck");
      if(diverged) {
        fvSolver().writeVtkXmlFiles("QOUT", "GEOM", false, true);
      }
    }
#endif
  } else {
    if(lpt().isActive()) {
      lpt().sendVelocitySlopes(); // the velocity slopes are received in the time-step
    }
  }
}

/** \brief Find combinded maximum timeStep for Fv and LPT solvers
 *         only possible for interleafed execution
 *         assuming that both solvers have previously computed their independend timeSteps
 *  \author Tim Wegmann
 */
template <MInt nDim, class SysEqn>
void CouplerFvParticle<nDim, SysEqn>::unifyTimeStep() {
  TRACE();

  ASSERT(m_interLeafed, "");

  cerr0 << "Unifying time-step LPT: " << lpt().m_timeStep << " FV: " << fvSolver().timeStep() << endl;

  const MFloat invalidTimeStep = std::numeric_limits<MFloat>::max();
  const MFloat conversionTime = conversionFvLpt.length / conversionFvLpt.velocity;
  MFloat timeFv = invalidTimeStep;
  MFloat timeLpt = invalidTimeStep;
  MInt maxLevelFv = -1;
  MInt maxLevelLpt = -1;
  MFloat newTimeFv = invalidTimeStep;
  MFloat newTimeLpt = invalidTimeStep;

  // determine the new/combined timeStep on ranks which have both Fv and Lpt solvers!
  if(fvSolver().isActive() && lpt().isActive()) {
    timeFv = fvSolver().timeStep(true);
    maxLevelFv = fvSolver().maxLevel();
    timeLpt = lpt().timeStep() / conversionTime;
    maxLevelLpt = lpt().maxLevel();

    MInt levelDiff = maxLevelLpt - maxLevelFv;
    if(levelDiff == 0) {
      newTimeFv = mMin(timeFv, timeLpt);
      newTimeLpt = newTimeFv * conversionTime;
    } else {
      if(levelDiff == 1) {
        // LPT has higher level => lower timeStep!
        newTimeLpt = mMin(newTimeLpt, newTimeFv / 2.0);
        newTimeFv = 2.0 * newTimeLpt;
      } else if(levelDiff == -1) {
        // Fv has the higher level => lower timeStep!
        newTimeFv = mMin(newTimeFv, newTimeLpt / 2.0);
        newTimeLpt = 2.0 * newTimeFv;
      } else {
        mTerm(1, AT_, "TimeStepping currently only implemented for 1 level difference!");
      }
    }
  }

  if(lpt().isActive()) {
    MPI_Allreduce(MPI_IN_PLACE, &newTimeLpt, 1, MPI_DOUBLE, MPI_MIN, lpt().mpiComm(), AT_, "MPI_IN_PLACE",
                  "newTimeLpt");
  }
  if(fvSolver().isActive()) {
    MPI_Allreduce(MPI_IN_PLACE, &newTimeFv, 1, MPI_DOUBLE, MPI_MIN, fvSolver().mpiComm(), AT_, "MPI_IN_PLACE",
                  "newTimeFv");
  }

  if(fvSolver().domainId() == 0) {
    cerr << "Fv timestep was " << timeFv << " combined TS is " << newTimeFv << endl;
  }

  if(lpt().domainId() == 0) {
    cerr << "LPT timestep was " << timeLpt << " combined TS is " << newTimeLpt << endl;
  }


  if(fvSolver().isActive()) {
    fvSolver().forceTimeStep(newTimeFv);
  }
  if(lpt().isActive()) {
    lpt().forceTimeStep(newTimeLpt);
  }
}

/** \brief transfer the FV velocity slopes to the LPT solver
 *  \author Tim Wegmann
 */
template <MInt nDim, class SysEqn>
void CouplerFvParticle<nDim, SysEqn>::subCouple(MInt recipeStep, MInt solverId,
                                                std::vector<MBool>& /*solverCompleted*/) {
  if(!m_interLeafed) return;

  if(fvSolver().isActive() && lpt().isActive()) {
    if(globalTimeStep % lpt().a_timeStepComputationInterval() != 0
       && fabs(lpt().m_timeStep - fvSolver().timeStep()) > MFloatEps) {
      mTerm(1, AT_, "TimeSteps of FV and LPT differ at TS " + to_string(globalTimeStep));
    }
  }

  // only during FV-Lpt computation
  if(recipeStep == m_lptSolverOrder && (solverId == m_fvSolverId || solverId == m_lptSolverId)) {
    m_solutionStep++;

    // meaning after the 4th LPT solution step!
    if(m_solutionStep == 2 * m_noSolutionSteps - 3) {
      if(solverId == m_lptSolverId) {
        // set the source-terms before the last FV-Solution(=RK) Step
        // and before the 4th RK step for FvMb

        if(fvSolver().m_hasExternalSource && !lpt().m_skipLPT) {
          fvSolver().resetExternalSources();
          transferExternalSources();
        } else if(lpt().m_skipLPT) {
          vector<MFloat> tmp(2 + nDim);
          for(MInt j = 0; j < (2 + nDim); j++) {
            tmp[j] = 0.0;
          }
          fvSolver().m_vapourData.insert(make_pair(globalTimeStep, tmp));
        }
      }
    }

    if(m_solutionStep == 2 * m_noSolutionSteps - 2 && !lpt().m_skipLPT) {
      if(fvSolver().m_hasExternalSource) {
        fvSolver().applyExternalOldSource();
      }
    }

    // meaning after the last RK-step!
    if(m_solutionStep == 2 * m_noSolutionSteps) {
      if(solverId == m_fvSolverId && !lpt().m_skipLPT) {
        // after the last FV-solution step!
        transferFlowField();
      }
    }

    if(m_solutionStep == 2 * m_noSolutionSteps && (fvSolver().m_timeStepUpdated || lpt().m_timeStepUpdated)) {
      if(!m_forceFvTimeStep) {
        unifyTimeStep();
      } else {
        transferTimeStep();
      }
    }
  }
}

/** \brief compute particle-cell binning and write to file
 *  \author Tim Wegmann
 */
template <MInt nDim, class SysEqn>
void CouplerFvParticle<nDim, SysEqn>::writeParticleCellStats() {
  TRACE();

  // count particles in cells
  lpt().countParticlesInCells();

  if(!fvSolver().isActive()) {
    return;
  }

  // transfer to FV-grid
  transferNoParticlesInCell();

  // determine max. No Particles in a FV-cell
  MInt maxNoParticles = 0;
  MInt maxParticlesSum = 0;
  for(MInt cellId = 0; cellId < fvSolver().noInternalCells(); cellId++) {
    if(fvSolver().a_noPart(cellId) > maxNoParticles) {
      maxNoParticles = fvSolver().a_noPart(cellId);
    }
    maxParticlesSum += fvSolver().a_noPart(cellId);
  }
  // create binning ranges
  MPI_Allreduce(MPI_IN_PLACE, &maxNoParticles, 1, MPI_INT, MPI_MAX, fvSolver().mpiComm(), AT_, "MPI_IN_PLACE",
                "maxParticles");

  MPI_Allreduce(MPI_IN_PLACE, &maxParticlesSum, 1, MPI_INT, MPI_MAX, fvSolver().mpiComm(), AT_, "MPI_IN_PLACE",
                "maxParticlesSum");

  {
    MFloatScratchSpace noCells(maxNoParticles + 1, AT_, "noCells");
    noCells.fill(0);

    // perform cell-binning
    for(MInt cellId = 0; cellId < fvSolver().noInternalCells(); cellId++) {
      if(!fvSolver().c_isLeafCell(cellId)) {
        continue;
      }
      if(!fvSolver().a_isActive(cellId)) {
        continue;
      }

      const MInt pos = fvSolver().a_noPart(cellId);
      noCells(pos)++;
    }

    // communicate results
    MPI_Allreduce(MPI_IN_PLACE, &noCells(0), maxNoParticles + 1, MPI_DOUBLE, MPI_SUM, fvSolver().mpiComm(), AT_,
                  "MPI_IN_PLACE", "noCells");

    // write to file
    if(fvSolver().domainId() == 0) {
      struct stat buffer {};
      MString name = "noParticlesInCells_" + to_string(fvSolver().solverId()) + "_" + to_string(globalTimeStep);
      const char* cstr = name.c_str();
      if(stat(cstr, &buffer) == 0) {
        rename(cstr, "noParticlesInCells_BAK");
      }
      ofstream ofl;
      ofl.open(name, ios_base::out | ios_base::app);
      ofl << "# 1:noParticles 2:noCells " << endl;
      for(MInt i = 0; i < maxNoParticles + 1; i++) {
        ofl << i << " " << noCells(i) << endl;
      }
      ofl.close();
    }
  }

  // perform min-cell binning
  MIntScratchSpace noMinCells(maxParticlesSum + 1, AT_, "noMinCells");
  noMinCells.fill(0);

  // reset particles to zero on halo-cells
  for(MInt cellId = fvSolver().noInternalCells(); cellId < fvSolver().a_noCells(); cellId++) {
    fvSolver().a_noPart(cellId) = 0;
  }

  MFloat minCellVolume = -1;
  for(MInt id = 0; id < fvSolver().grid().noMinCells(); id++) {
    const MInt cellId = fvSolver().grid().minCell(id);
    if(cellId < 0) {
      continue;
    }
    minCellVolume = fvSolver().c_cellVolumeAtLevel(fvSolver().a_level(cellId));

    // loop over all children and sum particles
    MInt curCount = childLoop(cellId);

    if(curCount > maxParticlesSum || curCount < 0) {
      cerr << "Strange particle count!" << curCount << endl;
    }

    noMinCells(curCount)++;
  }

  // communicate results
  MPI_Allreduce(MPI_IN_PLACE, &noMinCells(0), maxParticlesSum + 1, MPI_INT, MPI_SUM, fvSolver().mpiComm(), AT_,
                "MPI_IN_PLACE", "noMinCells");

  // write to file
  if(fvSolver().domainId() == 0) {
    struct stat buffer {};
    MString name = "noParticlesInMinCells_" + to_string(fvSolver().solverId()) + "_" + to_string(globalTimeStep);
    const char* cstr = name.c_str();
    if(stat(cstr, &buffer) == 0) {
      rename(cstr, "noParticlesInMinCells_BAK");
    }

    MString header = "# Min-cell volume : " + to_string(minCellVolume) + " \n";
    ofstream ofl;
    ofl.open(name, ios_base::out | ios_base::trunc);
    ofl << header << endl;

    for(MInt i = 0; i < maxParticlesSum + 1; i++) {
      ofl << i << " " << noMinCells(i) << endl;
    }
    ofl.close();
  }
}

/** \brief compute particle-cell binning and write to file
 *  \author Tim Wegmann
 */
template <MInt nDim, class SysEqn>
MInt CouplerFvParticle<nDim, SysEqn>::childLoop(MInt cellId) {
  MInt sum = 0;
  if(fvSolver().c_noChildren(cellId) > 0) {
    for(MInt child = 0; child < fvSolver().c_noChildren(cellId); child++) {
      MInt childId = fvSolver().c_childId(cellId, child);
      if(childId < 0) {
        continue;
      }
      sum += childLoop(childId);
    }
  } else {
    sum += fvSolver().a_noPart(cellId);
  }
  return sum;
}

template class CouplerFvParticle<3, FvSysEqnNS<3>>;
template class CouplerFvParticle<3, FvSysEqnRANS<3, RANSModelConstants<RANS_SA_DV>>>;
template class CouplerFvParticle<3, FvSysEqnRANS<3, RANSModelConstants<RANS_FS>>>;
template class CouplerFvParticle<3, FvSysEqnRANS<3, RANSModelConstants<RANS_KOMEGA>>>;
