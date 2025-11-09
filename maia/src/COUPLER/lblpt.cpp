// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "lblpt.h"
#include "LPT/lpt.h"
#include "LPT/lptlib.h"

using namespace std;
using namespace maia::lpt;

template <MInt nDim, MInt nDist, class SysEqn>
LbLpt<nDim, nDist, SysEqn>::LbLpt(const MInt couplingId, LPT<nDim>* particle, LbSolver* lb)
  : Coupling(couplingId), CouplingLpt<nDim, CouplingLB<nDim, nDist, SysEqn>>(couplingId, particle, lb) {
  TRACE();

  m_particle = particle;
  m_lbSolver = lb;

  // get LPT-solver order
  m_noSolverSteps = 1;
  m_noSolverSteps = Context::getBasicProperty<MInt>("recipeMaxNoSteps", AT_, &m_noSolverSteps);
  const MString propName = "solverOrder_" + std::to_string(lpt().solverId());
  for(MInt step = 0; step < m_noSolverSteps; step++) {
    if(Context::getBasicProperty<MInt>(propName, AT_, step) > 0) {
      m_lptSolverOrder = step;
      break;
    }
  }
}

/** \brief Init coupling class
 *
 *  \author Julian Vorspohl <j.vorspohl@aia.rwth-aachen.de>
 *  \date May 2021
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbLpt<nDim, nDist, SysEqn>::init() {
  TRACE();

  m_lbSolverId = lbSolver().solverId();
  m_lptSolverId = lpt().solverId();

  readProperties();

  checkProperties();

  initConversion();

  calculateGridBoundaries();
}

/** \brief Calculates the conversion factor for different non-dimensionalisations
 *         in the LB and LPT solvers
 *  \author Johannes Grafen, adapted from FV
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbLpt<nDim, nDist, SysEqn>::initConversion() {
  TRACE();
  // TODO adjust:
  // 1) read the following conversion factors as properties
  MFloat lengthFactor = F1;
  /*! \page propertyPage1
  <code>MFloat lengthFactor</code>\n
  default = 1 \n \n
  Set L_refLB/L_refLPT the conversion between different reference length used in the non-dimensionalisation and the
  calculation of the solver specific Re-number in the LB and LPT solver. \n Keywords:
  <i>PARTICLE</i>, <i>LATTICE BOLTZMANN</i>
  */
  lengthFactor = Context::getSolverProperty<MFloat>("lengthFactor", lpt().solverId(), AT_, &lengthFactor);

  MFloat velocityFactor = F1;
  /*! \page propertyPage1
  <code>MFloat velocityFactor</code>\n
  default = 1 \n \n
  Set u_refLB/u_refLPT the conversion between different velocities used in the non-dimensionalisation and the
  calculation of the solver specific Re-number in the LB and LPT solver. \n Keywords: <i>PARTICLE</i>, <i>LATTICE
  BOLTZMANN</i>
  */
  velocityFactor = Context::getSolverProperty<MFloat>("velocityFactor", lpt().solverId(), AT_, &velocityFactor);

  MFloat viscosityFactor = F1;
  /*! \page propertyPage1
  <code>MFloat viscosityFactor</code>\n
  default = 1.0 \n \n
  The conversion-factor for LB reference viscosity (mu_ref at T_ref) to LPT viscosity (mu_ref at T_ref) as used in the
  solver-specific Re-numbers. The viscosityFactor itself is non-dimensional!
  \n Keywords: <i>PARTICLE</i>, <i>LATTICE BOLTZMANN</i>
  */
  viscosityFactor = Context::getSolverProperty<MFloat>("viscosityFactor", lpt().solverId(), AT_, &viscosityFactor);

  const MFloat particleRe = Context::getSolverProperty<MFloat>("Re", lpt().solverId(), AT_);

  // 2) compute the missing conversion factors
  const MFloat ReFactor = lbSolver().m_Re / particleRe;
  /*const */ MFloat densityFactor = ReFactor * lengthFactor / (velocityFactor * viscosityFactor);
  /*const */ // MFloat temperatureFactor = sqrt(velocityFactor);

  const MFloat dx = a_cellLengthAtLevel(lbSolver().maxLevel());

  //---------------------------------------------------------------------------------------

  if(lpt().m_Ma >= 0.3) mTerm(1, AT_, "LB-LPT assumes low Mach number!");

  conversionLbLpt.velocity = velocityFactor * F1BCS / lbSolver().m_Ma;
  conversionLptLb.velocity = F1 / conversionLbLpt.velocity;

  conversionLbLpt.density = densityFactor;
  conversionLptLb.density = F1 / conversionLbLpt.density;

  conversionLbLpt.pressure = densityFactor * POW2(conversionLbLpt.velocity);
  conversionLptLb.pressure = F1 / conversionLbLpt.pressure;

  conversionLbLpt.length = lengthFactor * dx;
  conversionLptLb.length = F1 / conversionLbLpt.length;

  /* conversionLbLpt.temperature = temperatureFactor;
  conversionLptLb.temperature = 1.0 / conversionLbLpt.temperature; */

  conversionLbLpt.viscosity = viscosityFactor;
  conversionLptLb.viscosity = F1 / conversionLbLpt.viscosity;

  conversionLbLpt.time = conversionLbLpt.length / conversionLbLpt.velocity;
  conversionLptLb.time = F1 / conversionLbLpt.time;

  conversionLbLpt.mass = conversionLbLpt.velocity * conversionLbLpt.density * POW2(conversionLbLpt.length);
  conversionLptLb.mass = F1 / conversionLbLpt.mass;

  conversionLbLpt.momentum = conversionLbLpt.density * POW2(conversionLbLpt.velocity) * POW2(conversionLbLpt.length);
  conversionLptLb.momentum = F1 / conversionLbLpt.momentum;

  conversionLbLpt.energy = conversionLbLpt.density * POW3(conversionLbLpt.velocity) * POW2(conversionLbLpt.length);
  conversionLptLb.energy = F1 / conversionLbLpt.energy;

  if(lpt().domainId() == 0) {
    cerr << "Lb-Lpt conversion factors are: " << endl;
    cerr << "length       : " << conversionLbLpt.length << endl;
    cerr << "velocity     : " << conversionLbLpt.velocity << endl;
    cerr << "acceleration : " << conversionLbLpt.velocity / conversionLbLpt.time << endl;
    cerr << "density      : " << conversionLbLpt.density << endl;
    cerr << "pressure     : " << conversionLbLpt.pressure << endl;
    cerr << "temperature  : " << conversionLbLpt.temperature << endl;
    cerr << "viscosity    : " << conversionLbLpt.viscosity << endl;
    cerr << "time         : " << conversionLbLpt.time << endl;
    cerr << "mass         : " << conversionLbLpt.mass << endl;
    cerr << "momentum     : " << conversionLbLpt.momentum << endl;
    cerr << "energy       : " << conversionLbLpt.energy << endl;
  }
}

/** \brief calculate Boundaries of box-shaped Grid for calculation of
 *         boundaryCellData
 *  \author Johannes
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbLpt<nDim, nDist, SysEqn>::calculateGridBoundaries() {
  TRACE();

  MFloat halfCellWidth = 0.0;
  for(MInt cellId = 0; cellId < a_noCells(); cellId++) {
    if(!lbSolver().c_isLeafCell(cellId)) continue;
    halfCellWidth = F1B2 * lbSolver().grid().cellLengthAtCell(cellId);
    for(MInt i = 0; i < nDim; i++) {
      m_gridBoundaries[2 * i] = mMin(m_gridBoundaries[2 * i], lbSolver().a_coordinate(cellId, i) - halfCellWidth);
      m_gridBoundaries[2 * i + 1] =
          mMax(m_gridBoundaries[2 * i + 1], lbSolver().a_coordinate(cellId, i) + halfCellWidth);
    }
  }

  for(MInt i = 0; i < nDim; i++) {
    MPI_Allreduce(MPI_IN_PLACE, &m_gridBoundaries[2 * i], 1, maia::type_traits<MFloat>::mpiType(), MPI_MIN,
                  lbSolver().mpiComm(), AT_, "MPI_IN_PLACE", "m_domainBoundaryMin");
    MPI_Allreduce(MPI_IN_PLACE, &m_gridBoundaries[2 * i + 1], 1, maia::type_traits<MFloat>::mpiType(), MPI_MAX,
                  lbSolver().mpiComm(), AT_, "MPI_IN_PLACE", "m_domainBoundaryMax");
  }
}

/** \brief Sets the initial particle velocity based on the flow field velocity in that cell.
 *         NOTE: This can not be done in the LPT initialCondition,
 *         as the FV flow field has not been transfered before, its set in the Lb initialCondition.
 *  \author Tim Wegmann
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbLpt<nDim, nDist, SysEqn>::initParticleVelocity() {
  TRACE();

  // TODO: inital velocity for ellipsoids
  /* if(!lpt().m_restart) {
    for(MInt i = 0; i < (MInt)lpt().m_partListEllipsoid.size(); i++) {
      const MInt lptCellId = lpt().m_partListEllipsoid[i].m_cellId;
      const MInt fvCellId = lpt2fvIdParent(lptCellId);
      for(MInt j = 0; j < nDim; j++) {
        lpt().m_partListEllipsoid[i].m_velocity[j] =
            lpt().m_partListEllipsoid[i].m_velocity[j] + fvSolver().a_pvariable(fvCellId, j);
      }
    } */

  // init fluid density and velocity, after the values have been transferred to the LPT solver!
  for(MInt i = 0; i < lpt().a_noParticles(); i++) {
    lpt().m_partList[i].initVelocityAndDensity();
  }
}

/** \brief Finalize coupler initialization, coupler is ready after this
 *
 *  \author Julian Vorspohl <j.vorspohl@aia.rwth-aachen.de>
 *  \date May 2021
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbLpt<nDim, nDist, SysEqn>::finalizeCouplerInit() {
  TRACE();

  // Transfer the flow field after the intial condition for LB has been applied
  transferFlowField();
  transferCellStatus();

  initParticleVelocity();

  // TODO: check
  if(!lpt().m_restartFile && lbSolver().m_restartFile) {
    if(lpt().domainId() == 0) {
      cerr << "Restart without LPT restartFile -> setting time from Lb-time!" << endl;
    }
    const MFloat conversionTime = conversionLbLpt.length / conversionLbLpt.velocity;
    const MFloat LbTime = lbSolver().m_time;
    lpt().m_time = LbTime * conversionTime;
  }

  // set the timeStep in the LPT solver
  // at this point the timeStep must be enforced by the lb-solver
  transferTimeStep();
}

/** \brief Read coupler properties
 *  \author Julian Vorspohl <j.vorspohl@aia.rwth-aachen.de>
 *  \date May 2021
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbLpt<nDim, nDist, SysEqn>::readProperties() {
  TRACE();

  m_lptlbInterpolation =
      Context::getSolverProperty<MBool>("allowLPTlbInterpolation", m_lbSolverId, AT_, &m_lptlbInterpolation);

  /*! \page propertyPage1
      \section cfl
    <code>MFloat COUPLER::m_forceLbTimeStep </code>\n
    default = <code>true</code>\n \n
    Specify whether the LB-solver alone should force the time-step in the LB-LPT coupling!
    Otherwise, the minimum timestep of LPT and LB solver is computed and enforce!
    This however requires a sub-coupling!
    <ul>
    true,false
    </ul>
    Keywords: <i> LPT, LB, COUPLING </i>
  */
  // m_forceLbTimeStep = true;
  m_forceLbTimeStep = Context::getSolverProperty<MBool>("forceLbTimeStep", m_lptSolverId, AT_, &m_forceLbTimeStep);

  m_CalcSurface = Context::getSolverProperty<MBool>("couplerCalcSurface", m_lptSolverId, AT_, &m_CalcSurface);

  /* if(!m_forceLbTimeStep && !m_interLeafed) {
    mTerm(1, AT_, "Coupled timesteps only working with interleafed setup!");
  } */
}

/** \brief Check coupler properties for validity
 *
 *  \author Julian Vorspohl <j.vorspohl@aia.rwth-aachen.de>
 *  \date May 2021
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbLpt<nDim, nDist, SysEqn>::checkProperties() {
  TRACE();

  cout << "in lblpt.cpp: m_Re = " << lpt().m_Re << " m_Re= " << lpt().m_Re << " m_Ma = " << lpt().m_Ma << "\n";
  const MFloat lptRe = Context::getSolverProperty<MFloat>("Re", lpt().solverId(), AT_);

  cout << "in lblpt.cpp lptRe = " << lptRe << "\n";

  // ASSERT(fabs(lptRe - lbSolver().m_Re) < MFloatEps, "LPT: " << lptRe << " LB " << lbSolver().m_Re);
  if(fabs(lptRe - lbSolver().m_Re) > MFloatEps) {
    cerr << "Warning: Re0 differs between lb- and LPT-solver!"
         << " lb " << lbSolver().m_Re << " LPT " << lptRe << endl;
  }

  // TODO:
  // ASSERT(lpt().a_timeStepComputationInterval() == fvSolver().a_timeStepComputationInterval(), "");
}

/**
 * \brief Coupling before each solver
 *
 * \author Julian Vorspohl <j.vorspohl@aia.rwth-aachen.de>
 * \date May 2021
 *
 * \param[in] recipeStep
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbLpt<nDim, nDist, SysEqn>::preCouple(MInt recipeStep) {
  TRACE();

  // PRE LPT
  if(recipeStep == m_lptSolverOrder) {
    // TODO m_solutionStep = 0;

    // TODO:
    /* if(fvSolver().m_hasExternalSource) {
      if(!m_interLeafed) lbSolver().resetExternalSources();
    } */
    if(lbSolver().m_particleMomentumCoupling) lbSolver().resetExternalSources();

    // update LPT bndryCells for wall-collision
    updateLPTBndry();
  }

  // TODO: implement for interleaved
  /* if(m_interLeafed) {
      const MFloat conversionTime = conversionFvLpt.length / conversionFvLpt.velocity;
      const MFloat fvTime = fvSolver().m_levelSetMb ? fvSolver().m_physicalTime : fvSolver().m_time;
      if(fabs(fvTime * conversionTime - lpt().m_time) > MFloatEps) {
        cerr << "LPT-time " << lpt().m_time << " Fv-time " << fvTime << endl;
        mTerm(1, AT_, "Miss-matching time!");
      }
    } */
}

/** \brief Coupling after each solver
 *
 *  \author Julian Vorspohl <j.vorspohl@aia.rwth-aachen.de>, Johannes Grafen
 *  \date May 2021
 *
 *  \param[in] recipeStep Current step in the execution recipe
 *
 *  POST LPT: Transfer coupling terms from LPT to LB
 *  POST LB : Transfer flow field from LB to LPT
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbLpt<nDim, nDist, SysEqn>::postCouple(MInt recipeStep) {
  TRACE();

  // recipe-order is switched for ellipsoids
  // if(lpt().m_ellipsoids) {
  //  recipeStep = recipeStep == 1 ? 0 : 1;
  //}

  // Post LPT (as the LPT-solver is called before the lb-solver!)
  if(recipeStep == m_lptSolverOrder) {
    updateLbSolver();
  }

  // POST LB
  MInt postLbStep = m_lptSolverOrder < m_noSolverSteps - 1 ? m_lptSolverOrder + 1 : m_lptSolverOrder - 1;
  if(recipeStep == postLbStep) {
    transferFlowField();

    transferCellStatus();
    transferTimeStep();

    // set particle velocity to fluid velocity in first timestep or in TimeStep in which desired velocitySkewness is
    // reached
    if((!(lpt().m_skewnessTimeStep == -1) && (globalTimeStep == 1 && lpt().m_skewnessTimeStep == 0))
       || (globalTimeStep == lpt().m_skewnessTimeStep)) {
      for(auto& part : lpt().m_partList) {
        for(MInt i = 0; i < nDim; i++) {
          MInt lbCellId = lpt2lbId(part.m_cellId);
          part.m_velocity[i] = part.m_oldVel[i] =
              lbSolver().a_variable(lbCellId, lbSolver().PV->VV[i]) * conversionLbLpt.velocity;
        }
      }
    }
  }
}

/** \brief postCouple: exchange source terms after the LPT timeStep
 *  \author Tim Wegmann
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbLpt<nDim, nDist, SysEqn>::finalizeBalance(const MInt id) {
  TRACE();

  std::cout << "in finalize Balance" << std::endl;
  // post LPT solver
  if(id == lpt().solverId()) {
    lbSolver().resetExternalSources();
    if(globalTimeStep < 0) {
      // prepareAdaptation();
    }

  } else if(id == lbSolver().solverId()) {
    transferFlowField();
    transferCellStatus();

    updateLPTBndry();
  }
}

/** \brief Transfer all relevant data from LPT to LB solver.
 *
 *  \author Tim Wegmann
 *  \date May 2021
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbLpt<nDim, nDist, SysEqn>::updateLbSolver() {
  TRACE();

  // transfer external Sources to LB-solver
  if(!lbSolver().isActive() || !lpt().m_momentumCoupling) return;

    // check if noInternalCells must be used
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(MInt lbCellId = 0; lbCellId < lbSolver().a_noCells(); lbCellId++) {
    if(!lbSolver().a_isActive(lbCellId)) continue;
    // check that externalForces is set to zeros
    for(MInt i = 0; i < nDim; i++)
      ASSERT(abs(lbSolver().a_externalForces(lbCellId, i)) < MFloatEps, "");
  }

  // TODO: interleaved
  /*
  const MInt lbMaxLevel = lbSolver().maxLevel();
  #ifdef _OPENMP
  #pragma omp parallel for
  #endif
    for(MInt lptCell = 0; lptCell < lpt().a_noCells(); lptCell++) {
      // TODO-timw: check lpt-exchange!  use below instead
      // for(MInt lptCell = 0; lptCell < lpt().noInternalCells(); lptCell++) {
      // the source term is only set on LPT leaf-cells
      if(!lpt().c_isLeafCell(lptCell)) continue;
      // TODO: always keep noParticles in cell up to date,
      //      so that cells without particles easily be skipped
      MInt lbCellId = lpt2lbId(lptCell);
      if(lbCellId > -1 && lbSolver().c_isLeafCell(lbCellId)) {
        const MInt lbLevel = lbSolver().a_level(lbCellId);
        const MFloat fLbCellVolume = FFPOW2(lbMaxLevel - lbLevel);
        // simple transfer
        for(MInt dir = 0; dir < nDim; dir++) {
          lbSolver().a_externalForces(lbCellId, dir) -=
              lpt().a_momentumFlux(lptCell, dir) * conversionLptLb.momentum * fLbCellVolume * m_factor;
          if(abs(lpt().a_momentumFlux(lptCell, dir)) > 1e-5)
            std::cerr << "lbforce " << lbSolver().a_externalForces(lbCellId, dir) << " conv " <<
  conversionLptLb.momentum
                      << " fLbCV " << fLbCellVolume << " flux " << lpt().a_momentumFlux(lptCell, dir) << std::endl;
        }
      }
    }
   */

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(MInt lptCellId = 0; lptCellId < lpt().a_noCells(); lptCellId++) {
    if(!lpt().c_isLeafCell(lptCellId)) continue;
    const MInt lbCellId = lpt2lbId(lptCellId);
    if(lbCellId < 0) continue;
    if(lbCellId > -1 && lbSolver().c_isLeafCell(lbCellId)) {
      for(MInt i = 0; i < nDim; i++) {
        // Think about different refinement levels of lb-cells and lpt-cells like in fvparticle.cpp ->
        // transferExternalSources()
        // obtain reaction force density related to cellvolume in lb-Units acting on fluid
        lbSolver().a_externalForces(lbCellId, i) -= lpt().a_momentumFlux(lptCellId, i) * conversionLptLb.momentum;
      }
    }
    // TODO different refinement of LPT an LB
  }
  lbSolver().exchangeExternalSources();
}

// TODO: implement interleaved version
/** \brief Transfer the momentum source (forcing) from LPT to LB.
 *  \author Tim Wegmann, Daniel Lauwers
 */
/*template <MInt nDim, MInt nDist, class SysEqn>
void LbLpt<nDim, nDist, SysEqn>::updateForcing() {
  TRACE();

  if(!lbSolver().isActive()) return;

  if(!lpt().m_momentumCoupling && !lpt().m_heatCoupling && !lpt().m_massCoupling) return;

  // TODO-timw: check lpt-exchange! add
  // lpt().exchangeSourceTerms();

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(MInt cellId = 0; cellId < lbSolver().a_noCells(); cellId++) {
    for(MInt var = 0; var < nDim; var++) {
      lbSolver().a_cellForce(cellId, var) = 0.0;
    }
  }

  const MInt lbMaxLevel = lbSolver().maxLevel();
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(MInt lptCell = 0; lptCell < lpt().a_noCells(); lptCell++) {
    // TODO-timw: check lpt-exchange!  use below instead
    // for(MInt lptCell = 0; lptCell < lpt().noInternalCells(); lptCell++) {
    // the source term is only set on LPT leaf-cells
    if(!lpt().c_isLeafCell(lptCell)) continue;
    // TODO: always keep noParticles in cell up to date,
    //      so that cells without particles easily be skipped
    MInt lbCellId = lpt2lbId(lptCell);
    if(lbCellId > -1 && lbSolver().c_isLeafCell(lbCellId)) {
      const MInt lbLevel = lbSolver().a_level(lbCellId);
      const MFloat fLbCellVolume = FFPOW2(lbMaxLevel - lbLevel);
      // simple transfer
      for(MInt dir = 0; dir < nDim; dir++) {
        lbSolver().a_cellForce(lbCellId, dir) -=
            lpt().a_momentumFlux(lptCell, dir) * conversionLptLb.momentum * fLbCellVolume;
        // if(abs(lpt().a_momentumFlux(lptCell, dir)) > 1e-5)
        //   std::cerr << "lbforce " << lbSolver().a_cellForce(lbCellId, dir) << " conv " << conversionLptLb.momentum
        //             << " fLbCV " << fLbCellVolume << " flux " << lpt().a_momentumFlux(lptCell, dir) << std::endl;
      }

    } else if(lbCellId > -1 && !lbSolver().c_isLeafCell(lbCellId)) {
      // lb-cell is refined further: volume average from LPT to LB
      if(lpt().a_isHalo(lptCell)) continue;
      // gather all lb-leaf-cell childs for internal cells
      vector<MInt> leafChildIds;
      lbSolver().grid().getAllLeafChilds(lbCellId, leafChildIds);
      MFloat totalVol = 0;
      for(MUint i = 0; i < leafChildIds.size(); i++) {
        const MInt cellId = leafChildIds[i];
        totalVol += lbSolver().grid().cellVolumeAtLevel(lbSolver().c_level(cellId));
      }
      for(MUint i = 0; i < leafChildIds.size(); i++) {
        const MInt cellId = leafChildIds[i];
        const MInt lbLevel = lbSolver().a_level(cellId);
        const MFloat fLbCellVolume = FFPOW2(lbMaxLevel - lbLevel);
        const MFloat vFrac = lbSolver().grid().cellVolumeAtLevel(lbSolver().c_level(cellId)) / totalVol;
        for(MInt dir = 0; dir < nDim; dir++) {
          lbSolver().a_cellForce(cellId, dir) -=
              lpt().a_momentumFlux(lptCell, dir) * vFrac * conversionLptLb.momentum * fLbCellVolume;
        }
      }
    } else {
      mTerm(1, AT_, "check if race condition exists here!");
      // lpt-Cell is refined further
      lbCellId = lpt2lbIdParent(lptCell);
      if(lbCellId < 0) continue;
      const MInt lbLevel = lbSolver().a_level(lbCellId);
      const MFloat fLbCellVolume = FFPOW2(lbMaxLevel - lbLevel);
      vector<MInt> leafChildIds;
      lpt().grid().getAllLeafChilds(lptCell, leafChildIds);

      // calculate sum of all source terms from all lpt leaf cells
      for(MUint i = 0; i < leafChildIds.size(); i++) {
        const MInt cellId = leafChildIds[i];
        for(MInt dir = 0; dir < nDim; dir++) {
          lbSolver().a_cellForce(lbCellId, dir) -=
              lpt().a_momentumFlux(cellId, dir) * conversionLptLb.momentum * fLbCellVolume;
        }
      }
    }
  }
}
*/

/** \brief Transfer all flow variables from LB to LPT and apply the conversion
 *
 *  \author Julian Vorspohl <j.vorspohl@aia.rwth-aachen.de>
 *  \date May 2021
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbLpt<nDim, nDist, SysEqn>::transferFlowField() {
  TRACE();
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(MInt lbCellId = 0; lbCellId < lbSolver().noInternalCells(); lbCellId++) {
    if(!lbSolver().a_isActive(lbCellId)) {
      continue;
    }

    const MInt lptCellId = lb2lptId(lbCellId);
    if(lptCellId < 0) continue;

    for(MInt n = 0; n < nDim; n++) {
      lpt().a_fluidVelocity(lptCellId, n) =
          lbSolver().a_variable(lbCellId, lbSolver().PV->VV[n]) * conversionLbLpt.velocity;
    }
    lpt().a_fluidDensity(lptCellId) = conversionLbLpt.density * lbSolver().a_variable(lbCellId, lbSolver().PV->RHO);

    // Assume incompressible LB
    lpt().a_fluidPressure(lptCellId) =
        (lbSolver().a_variable(lbCellId, lbSolver().PV->RHO) - 1.0) * F1B3 * conversionLbLpt.pressure + 1.0;
    // TODO: @Julian reactivate for thermal
    /*lpt().a_fluidTemperature(lptCellId) = lbSolver().a_pvariable(lbCellId, lbSolver().PV->T)
     * conversionLbLpt.temperature;*/
  }

  // NOTE: exchange is necessary as particles on window-cells interpolate the flow variables to their
  //       position and will need values on halo-cells to do so!
  lpt().exchangeData(&lpt().a_fluidVariable(0, 0), lpt().PV.noVars());
  if(lpt().m_evaporation) {
    lpt().exchangeData(&lpt().a_fluidSpecies(0), 1);
  }

  lpt().checkCells();
}

/** \brief Transfer all relevant bndryCell-data from LB to LPT.
 *
 *  \author Tim Wegmann, Julian Vorspohl <j.vorspohl@aia.rwth-aachen.de>, Johannes Grafen
 *  \date Dec 2021
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbLpt<nDim, nDist, SysEqn>::updateLPTBndry() {
  TRACE();

  if(!lpt().m_wallCollisions) return;

  if(!lpt().isActive()) return;

  MFloatScratchSpace bBox(2 * nDim, AT_, "bBox");
  lbSolver().m_geometry->getBoundingBox(bBox.getPointer());

  /*
  using Ld = LbLatticeDescriptor<nDim, nDist>;
  // just number of 2D distributions
  static constexpr MInt dist1 = Ld::distFld(0);
  static constexpr MInt maxDist1Idx = dist1; // for loops
  static constexpr MInt dist2 = Ld::distFld(1);
  static constexpr MInt maxDist2Idx = dist1 + dist2; // for loops
  */
  // Just once if no mesh adaptation during run
  for(MInt bndryId = lbSolver().m_noOuterBndryCells; bndryId < a_noBndCells(); bndryId++) {
    const MInt lbCellId = a_bndCellId(bndryId);

    if(!lbSolver().c_isLeafCell(lbCellId)) continue;
    MInt lptCellId = lb2lptId(lbCellId);
    if(lptCellId < 0) {
      lptCellId = lb2lptIdParent(lbCellId);
      if(lptCellId < 0) continue;
      if(lpt().c_level(lptCellId) < lpt().maxRefinementLevel()) continue;
      ASSERT(m_lptlbInterpolation, "");
    }

    if(!lpt().c_isLeafCell(lptCellId)) {
      ASSERT(m_lptlbInterpolation, "");
    }

    const MInt noSurf = 1; // lbSolver().m_bndCnd->m_bndCells->a[bndryId].m_noSrfcs

    // bndryCell without a bndrySurface, how is that even possible?
    // if(noSurf < 1 ) continue;

    ASSERT(noSurf == 1, "Only for single bndry-surfaces!");

    MFloatScratchSpace surfaceC(nDim, nDim + 1, AT_, "surfaceC");
    std::array<MFloat, nDim> surfaceN;
    std::array<MFloat, nDim> surfaceV;
    std::array<MFloat, nDim + 1> bndryData;

    // TODO: check correctness of addbndryCellMethod before reusing!
    /* if(m_CalcSurface) {
      const MFloat dx = lbSolver().c_cellLengthAtCell(lbCellId);
      std::vector<const MFloat*> cellWallNormals; // store all Wall normals
      std::vector<const MFloat*>
          cellPlaneCoordinateVec; // 2D Distribution from cell center to point in the middle of an surface edge
      const MFloat* SCdist{};
      MInt wallCount = 0;
      std::array<MBool, nDim> gridBoxDifference{};

      for(MInt i = 0; i < nDim; i++) {
        bndryData.at(1 + i) = lbSolver().a_coordinate(lbCellId, i);
      }

      for(MInt j = 0; j < maxDist1Idx; j++) {
        if(lbSolver().a_hasNeighbor(lbCellId, j) == 0) {
          SCdist = Ld::ppdfDir(j); // only for 1D, dist from cell center to SC
          cellWallNormals.push_back(Ld::ppdfDir(Ld::oppositeDist(j)));
          wallCount++;
        }
      }

      // compute surface normals out of cellWallNormals, should be independent from number of distributions
      for(MInt n = 0; n < nDim; n++) {
        surfaceV.at(n) = 0; // for simple simulation not needed, can be set to zero;
        for(size_t a = 0; a < cellWallNormals.size(); a++)
          surfaceN.at(n) += cellWallNormals.at(a)[n]; // Assuming wallnormals have just one vec component unequal 0
      }
      maia::math::normalize(surfaceN);

      switch(wallCount) {
        case 0:
          std::cerr << "Something went wrong, " << lbCellId << " is not a boundaryCell" << std::endl;
          break;

        case 1:
          // set bndry-cell volume:
          bndryData.at(0) = POW3(dx);
          // Iterate over 2D distributions
          for(MInt j = maxDist1Idx; j < maxDist2Idx; j++) {
            for(MInt n = 0; n < nDim; n++) {
              // Check for -1 * -1 or 1 * 1; -1 * 1 should not pass since SCdist[n] and ppdf(j,n) are not equal.
              if(SCdist[n] * Ld::ppdfDir(j, n) + MFloatEps > 1.0) {
                cellPlaneCoordinateVec.push_back(Ld::ppdfDir(j));
                if(cellPlaneCoordinateVec.size() == 3) break; // only 3 points are needed
              }
            }
          }

          for(MInt n = 0; n < nDim; n++) {
            gridBoxDifference.at(n) = false;
            surfaceC.at(0).at(n) = bndryData.at(1 + n) - surfaceN.at(n) * dx / 2;
            if(surfaceN.at(n) > 0.0) {
              if(bBox[n] > m_gridBoundaries[2 * n]) {
                surfaceC.at(0).at(n) = bBox[n];
                gridBoxDifference.at(n) = true;
              }
            } else if(surfaceN.at(n) < 0.0) {
              if(bBox[n + nDim] < m_gridBoundaries[2 * n + 1]) {
                surfaceC.at(0).at(n) = bBox[n + nDim];
                gridBoxDifference.at(n) = true;
              }
            }
          }
          break;

        case 2:
          // set bndry-cell volume:
          bndryData.at(0) = POW3(dx) / 2;
          // Iterate over 1D and 2D distributions
          for(MInt j = 0; j < maxDist2Idx; j++) {
            std::array<MFloat, nDim> tmpDist{};
            for(MInt n = 0; n < nDim; n++) {
              tmpDist[n] = Ld::ppdfDir(j, n);
            }
            // search for 1D and 2D distributions that are perpendicular to surface normal vector
            if(abs(std::inner_product(surfaceN.cbegin(), surfaceN.cend(), tmpDist.cbegin(), 0.0)) < MFloatEps)
              cellPlaneCoordinateVec.push_back(Ld::ppdfDir(j));
            if(cellPlaneCoordinateVec.size() == 3) break;
          }

          for(MInt n = 0; n < nDim; n++) {
            surfaceC.at(0).at(n) = bndryData.at(1 + n);
          }
          break;

        case 3:
          // set bndry-cell volume:
          bndryData.at(0) = POW3(dx) * (1.0 - 0.25 * (0.5 * SQRT3 - sqrt(3 - SQRT6)));
          // Iterate over 1D distributions
          for(MInt j = 0; j < maxDist1Idx; j++) {
            std::array<MFloat, nDim> tmpDist{};
            for(MInt n = 0; n < nDim; n++) {
              tmpDist[n] = Ld::ppdfDir(j, n);
            }
            // search for 1D distributions that are pointing in opposite halfspace of surface normal -> va * vb = abs(a)
            // * abs(b) * cos(phi) < 0 -> yes
            if(std::inner_product(surfaceN.cbegin(), surfaceN.cend(), tmpDist.cbegin(), 0.0) < MFloatEps)
              cellPlaneCoordinateVec.push_back(Ld::ppdfDir(j));
            if(cellPlaneCoordinateVec.size() == 3) break;
          }

          for(MInt n = 0; n < nDim; n++) {
            surfaceC.at(0).at(n) = bndryData.at(1 + n) - surfaceN.at(n) * sqrt(3 - SQRT6) * dx; //? check
          }
          break;

        default:
          std::cerr << "No implementation for cells with more than 3 bounding walls!" << std::endl;
          break;
      }

      // find 3 Distributions to get 3 points on surface edges
      MFloat dL = 0.0;
      // loop over coordinate x,y,z
      for(MInt i = 0; i < nDim; i++) {
        if(gridBoxDifference.at(i))
          dL = bndryData.at(1 + i) - surfaceC.at(0).at(i);
        else
          dL = dx / 2;
        // loop over surface point 1,2,3
        for(MInt n = 0; n < nDim; n++) {
          surfaceC.at(1 + n).at(i) = bndryData.at(1 + i) + cellPlaneCoordinateVec.at(n)[i] * dL;
        }
      } */
    // Printing boundary surface information :
    /* std::cout << "lbcellId: " << lbCellId << std::endl
              << "lptcellId: " << lptCellId << std::endl
              << "wallCount: " << wallCount << std::endl
              << "bndrycell Volume: " << bndryData.at(0) << std::endl;
    for(MInt i = 0; i < nDim; i++) std::cout << "bndryData[" << bndryId << "].m_coordinates "<< i << ": " <<
    bndryData.at(1 + i) << std::endl; for(MInt i = 0; i < nDim; i++) std::cout << "surfaceN[ " << i << "] normalized:
    " << surfaceN.at(i) << std::endl; for(MInt i = 0; i < nDim; i++) std::cout << "surfaceC[0][" << i << "]: " <<
    surfaceC.at(0).at(i) << std::endl;

    for(MInt i = 0; i < nDim; i++)
      for(MInt n = 0; n <nDim; n++)
        std::cout << "cellPlaneCoordinateVec[" << i <<"]["<<n<<"]: " << cellPlaneCoordinateVec.at(i)[n] << std::endl;
    for(MInt i = 0; i < nDim; i++)
      for(MInt n = 0; n <nDim; n++)
        std::cout << "surfaceC[" << 1 + i << "]["<< n << "]: " <<  surfaceC.at(1 + i).at(n) << std::endl;
     */

    // clear all temporary vectors
    /* cellWallNormals.clear();
    cellPlaneCoordinateVec.clear();
  } */
    // append new LPT bndryCell based on LB data
    lpt().addBndryCell(lptCellId, &bndryData[0], &surfaceN[0], surfaceC, &surfaceV[0]);
  }
}

/** \brief transfer/enforce Fv timeStep onto LPT solver
 *  \author Tim Wegmann
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbLpt<nDim, nDist, SysEqn>::transferTimeStep() {
  TRACE();

  lpt().forceTimeStep(/*fvSolver().timeStep(true) */ conversionLbLpt.time);
}

// TODO:
/** \brief set the isValid status for LPT cells based on the FV-solver cell properties
 *  \author Tim Wegmann
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbLpt<nDim, nDist, SysEqn>::transferCellStatus() {
  TRACE();
  for(MInt lptCellId = 0; lptCellId < lpt().a_noCells(); lptCellId++) {
    MInt lbCellId = lpt2lbIdParent(lptCellId);
    if(lbCellId < 0 || lbCellId > lbSolver().a_noCells()) {
      lpt().a_isValidCell(lptCellId) = false;
    } else if(!lbSolver().a_isActive(lbCellId)) {
      lpt().a_isValidCell(lptCellId) = false;
    } else {
      lpt().a_isValidCell(lptCellId) = true;
    }
  }
}

// Explicit instantiations
template class LbLpt<3, 19, maia::lb::LbSysEqnIncompressible<3, 19>>;
template class LbLpt<3, 27, maia::lb::LbSysEqnIncompressible<3, 27>>;
