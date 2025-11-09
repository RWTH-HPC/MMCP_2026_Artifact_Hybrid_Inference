// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "lbsolverdxqy.h"

#include <algorithm>
#include "COMM/mpioverride.h"
#include "UTIL/maiafftw.h"
#include "UTIL/parallelfor.h"
#include "globals.h"
#include "lbbndcnddxqy.h"
#include "lbconstants.h"
#include "lbfunctions.h"
#include "lbgridboundarycell.h"
#include "lbinterfacedxqy.h"
#include "lbsyseqn.h"

using namespace std;
using namespace lbconstants;

template <MInt nDim, MInt nDist, class SysEqn>
LbSolverDxQy<nDim, nDist, SysEqn>::LbSolverDxQy(MInt id, MInt noDistributions, GridProxy& gridProxy_,
                                                Geometry<nDim>& geometry_, const MPI_Comm comm)
  : LbSolver<nDim>(id, noDistributions, gridProxy_, geometry_, comm),
    m_srcTermController(this),
    C1(5.0),
    C2(-3.05),
    C3(2.5),
    C4(5.5) // Constants for Law of the Wall
{
  TRACE();

  if(m_isRefined || this->m_adaptation) m_interface = new LbInterfaceDxQy<nDim, nDist, SysEqn>(this);

  m_bndCnd = new LbBndCnd(this);
  LbSolver<nDim>::m_bndCnd = m_bndCnd;

  if(grid().isActive()) {
    this->prepareCommunication();
    if(m_geometry->m_parallelGeometry) {
      this->receiveWindowTriangles();
    }
  }

  m_log << "noCells(offset) =" << this->noInternalCells() << endl;
  m_log << "noCells(size) =" << this->m_cells.size() << endl;

#ifdef WAR_NVHPC_PSTL
  initArraysForPSTL();
#endif
}

/** \brief
 *
 */
template <MInt nDim, MInt nDist, class SysEqn>
LbSolverDxQy<nDim, nDist, SysEqn>::~LbSolverDxQy() {
  TRACE();
  if(m_isRefined || this->m_adaptation) delete m_interface;
  delete m_bndCnd;
}

/** \brief standard OpenMP propagation step
 * \author Andreas Lintermann
 * \date 03.04.2017
 *
 * This function propagates the locally calculated PPDFs to the neighboring cells after the collision step.
 *
 **/
template <MInt nDim, MInt nDist, class SysEqn>
void LbSolverDxQy<nDim, nDist, SysEqn>::propagation_step() {
  TRACE();

  // MInt distStart, distNeighStart, n1Start, n2Start;
  // inner loop references
  const MInt noCells = a_noCells();

  // Required for ported functions, since GPUs cannot access the global variable globalTimeStep
  const MInt gTS = globalTimeStep;
  const MInt maxLevel_ = maxLevel();

  maia::parallelFor<true>(0, noCells, [=](MInt i) {
    if(gTS % IPOW2(maxLevel_ - this->a_level(i)) != 0) return;
    const MInt lastId = nDist - 1;
    for(MInt j = 0; j < lastId; ++j) {
      if(auto n = c_neighborId(i, j); n > -1) {
        a_oldDistribution(n, j) = a_distribution(i, j);
      }
    }
    a_oldDistribution(i, lastId) = a_distribution(i, lastId);
  });

  postPropagationSrcTerm();
  volumeForces();

  if(m_isRefined) {
    restriction();
    prolongation();
  }
}


/** \brief
 *
 * This function propagates the locally calculated PPDFs to the neighboring cells after the collision step.
 * The volumetric refinement scheme requires prolongation to take place before propagation
 *
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbSolverDxQy<nDim, nDist, SysEqn>::propagation_step_vol() {
  TRACE();

  MInt nghbrId, nghbrId2;

  if(m_isRefined) prolongation();

  for(MInt i = 0; i < a_noCells(); i++) {
    if((globalTimeStep) % IPOW2(maxLevel() - this->a_level(i)) == 0) {
      for(MInt j = 0; j < nDist - 1; j += 2) {
        nghbrId = c_neighborId(i, j);
        nghbrId2 = c_neighborId(i, j + 1);

        if(nghbrId > -1) a_oldDistribution(nghbrId, j) = a_distribution(i, j);
        if(nghbrId2 > -1) a_oldDistribution(nghbrId2, j + 1) = a_distribution(i, j + 1);
      }
      a_oldDistribution(i, nDist - 1) = a_distribution(i, nDist - 1);
    }
  }

  restriction();

  //--

  // #define WRITECELLVALUES
#ifdef WRITECELLVALUES
  ofstream ofl;

  if((globalTimeStep) % 40 == 0) {
    if(domainId() == 7) {
      ofl.open("x0p75.dat", ios_base::app);
      ofl << globalTimeStep << " " << a_variable(212707, 0) << " " << a_variable(212707, 1) << " "
          << a_variable(212707, 2) << " " << a_variable(212707, 3) << endl;
      ofl.close();

      ofl.open("x1p5.dat", ios_base::app);
      ofl << globalTimeStep << " " << a_variable(235107, 0) << " " << a_variable(235107, 1) << " "
          << a_variable(235107, 2) << " " << a_variable(235107, 3) << endl;
      ofl.close();
    }
  }
#endif

  //--
}

template <MInt nDim, MInt nDist, class SysEqn>
void LbSolverDxQy<nDim, nDist, SysEqn>::propagation_step_thermal() {
  TRACE();

  // MInt distStart, distNeighStart, n1Start, n2Start;
  // inner loop references
  const MInt noCells = a_noCells();

  // Required for ported functions, since GPUs cannot access the global variable globalTimeStep
  const MInt gTS = globalTimeStep;
  const MInt maxLevel_ = maxLevel();

  maia::parallelFor<true>(0, noCells, [=](MInt i) {
    if(gTS % IPOW2(maxLevel_ - this->a_level(i)) != 0) return;
    const MInt lastId = nDist - 1;
    for(MInt j = 0; j < lastId; ++j) {
      if(auto n = c_neighborId(i, j); n > -1) {
        a_oldDistribution(n, j) = a_distribution(i, j);
        a_oldDistributionThermal(n, j) = a_distributionThermal(i, j);
      }
    }
    a_oldDistribution(i, lastId) = a_distribution(i, lastId);
    a_oldDistributionThermal(i, lastId) = a_distributionThermal(i, lastId);
  });

  if(m_isRefined) {
    restriction();
    prolongation();
  }
}

/** \brief Propagation step for Transport Lattice-Boltzmann
 * \author Shota Ito
 * \date 07.06.22
 *
 * This function propagates the locally computed PPDFs to the neighboring cells after the collision step.
 *
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbSolverDxQy<nDim, nDist, SysEqn>::propagation_step_transport() {
  TRACE();

  // MInt distStart, distNeighStart, n1Start, n2Start;
  // inner loop references
  const MInt noCells = a_noCells();

  // Required for ported functions, since GPUs cannot access the global variable globalTimeStep
  const MInt gTS = globalTimeStep;
  const MInt maxLevel_ = maxLevel();

  maia::parallelFor<true>(0, noCells, [=](MInt i) {
    if(gTS % IPOW2(maxLevel_ - this->a_level(i)) != 0) return;
    const MInt lastId = nDist - 1;
    for(MInt j = 0; j < lastId; ++j) {
      if(auto n = c_neighborId(i, j); n > -1) {
        a_oldDistribution(n, j) = a_distribution(i, j);
        a_oldDistributionTransport(n, j) = a_distributionTransport(i, j);
      }
    }
    a_oldDistribution(i, lastId) = a_distribution(i, lastId);
    a_oldDistributionTransport(i, lastId) = a_distributionTransport(i, lastId);
  });

  if(m_isRefined) {
    restriction();
    prolongation();
  }
}

/** \brief Propagation step for coupled Thermal Transport Lattice-Boltzmann
 * \author Shota Ito
 * \date 13.06.22
 *
 * This function propagates the locally computed PPDFs to the neighboring cells after the collision step.
 *
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbSolverDxQy<nDim, nDist, SysEqn>::propagation_step_thermaltransport() {
  TRACE();

  // MInt distStart, distNeighStart, n1Start, n2Start;
  // inner loop references
  const MInt noCells = a_noCells();

  // Required for ported functions, since GPUs cannot access the global variable globalTimeStep
  const MInt gTS = globalTimeStep;
  const MInt maxLevel_ = maxLevel();

  maia::parallelFor<true>(0, noCells, [=](MInt i) {
    if(gTS % IPOW2(maxLevel_ - this->a_level(i)) != 0) return;
    const MInt lastId = nDist - 1;
    for(MInt j = 0; j < lastId; ++j) {
      if(auto n = c_neighborId(i, j); n > -1) {
        a_oldDistribution(n, j) = a_distribution(i, j);
        a_oldDistributionThermal(n, j) = a_distributionThermal(i, j);
        a_oldDistributionTransport(n, j) = a_distributionTransport(i, j);
      }
    }
    a_oldDistribution(i, lastId) = a_distribution(i, lastId);
    a_oldDistributionThermal(i, lastId) = a_distributionThermal(i, lastId);
    a_oldDistributionTransport(i, lastId) = a_distributionTransport(i, lastId);
  });

  if(m_isRefined) {
    restriction();
    prolongation();
  }
}

/** \brief Propagation step for Thermal Lattice-Boltzmann
 * \author Andreas Lintermann
 * \date 23.02.2011
 *
 * This function propagates the locally calculated PPDFs to the neighboring cells after the collision step.
 * The volumetric refinement scheme requires prolongation to take place before propagation
 *
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbSolverDxQy<nDim, nDist, SysEqn>::propagation_step_thermal_vol() {
  TRACE();

  MInt nghbrId, nghbrId2;

  if(m_isRefined) prolongation();

  for(MInt i = 0; i < a_noCells(); i++) {
    if((globalTimeStep) % IPOW2(maxLevel() - this->a_level(i)) == 0) {
      for(MInt j = 0; j < nDist - 1; j += 2) {
        nghbrId = c_neighborId(i, j);
        nghbrId2 = c_neighborId(i, j + 1);
        if(nghbrId > -1) {
          a_oldDistribution(nghbrId, j) = a_distribution(i, j);
          a_oldDistributionThermal(nghbrId, j) = a_distributionThermal(i, j);
        }
        if(nghbrId2 > -1) {
          a_oldDistribution(nghbrId2, j + 1) = a_distribution(i, j + 1);
          a_oldDistributionThermal(nghbrId2, j + 1) = a_distributionThermal(i, j + 1);
        }
      }
      a_oldDistribution(i, nDist - 1) = a_distribution(i, nDist - 1);
      a_oldDistributionThermal(i, nDist - 1) = a_distributionThermal(i, nDist - 1);
    }
  }

  restriction();
}

/** \brief Propagation step for Transport Lattice-Boltzmann
 * \author Shota Ito
 * \date 07.06.22
 *
 * This function does the same as the propagation_step_thermal_vol but for Transport Lattice Boltzmann.
 *
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbSolverDxQy<nDim, nDist, SysEqn>::propagation_step_transport_vol() {
  TRACE();

  MInt nghbrId, nghbrId2;

  if(m_isRefined) prolongation();

  for(MInt i = 0; i < a_noCells(); i++) {
    if((globalTimeStep) % IPOW2(maxLevel() - this->a_level(i)) == 0) {
      for(MInt j = 0; j < nDist - 1; j += 2) {
        nghbrId = c_neighborId(i, j);
        nghbrId2 = c_neighborId(i, j + 1);
        if(nghbrId > -1) {
          a_oldDistribution(nghbrId, j) = a_distribution(i, j);
          a_oldDistributionTransport(nghbrId, j) = a_distributionTransport(i, j);
        }
        if(nghbrId2 > -1) {
          a_oldDistribution(nghbrId2, j + 1) = a_distribution(i, j + 1);
          a_oldDistributionTransport(nghbrId2, j + 1) = a_distributionTransport(i, j + 1);
        }
      }
      a_oldDistribution(i, nDist - 1) = a_distribution(i, nDist - 1);
      a_oldDistributionTransport(i, nDist - 1) = a_distributionTransport(i, nDist - 1);
    }
  }

  restriction();
}

/** \brief Propagation step for coupled Thermal Transport Lattice-Boltzmann
 * \author Shota Ito
 * \date 13.06.22
 *
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbSolverDxQy<nDim, nDist, SysEqn>::propagation_step_thermaltransport_vol() {
  TRACE();

  MInt nghbrId, nghbrId2;

  if(m_isRefined) prolongation();

  for(MInt i = 0; i < a_noCells(); i++) {
    if((globalTimeStep) % IPOW2(maxLevel() - this->a_level(i)) == 0) {
      for(MInt j = 0; j < nDist - 1; j += 2) {
        nghbrId = c_neighborId(i, j);
        nghbrId2 = c_neighborId(i, j + 1);
        if(nghbrId > -1) {
          a_oldDistribution(nghbrId, j) = a_distribution(i, j);
          a_oldDistributionThermal(nghbrId, j) = a_distributionThermal(i, j);
          a_oldDistributionTransport(nghbrId, j) = a_distributionTransport(i, j);
        }
        if(nghbrId2 > -1) {
          a_oldDistribution(nghbrId2, j + 1) = a_distribution(i, j + 1);
          a_oldDistributionThermal(nghbrId2, j + 1) = a_distributionThermal(i, j + 1);
          a_oldDistributionTransport(nghbrId2, j + 1) = a_distributionTransport(i, j + 1);
        }
      }
      a_oldDistribution(i, nDist - 1) = a_distribution(i, nDist - 1);
      a_oldDistributionThermal(i, nDist - 1) = a_distributionThermal(i, nDist - 1);
      a_oldDistributionTransport(i, nDist - 1) = a_distributionTransport(i, nDist - 1);
    }
  }

  restriction();
}

template <MInt nDim, MInt nDist, class SysEqn>
template <MInt timeStepOffset>
void LbSolverDxQy<nDim, nDist, SysEqn>::updateVariablesFromOldDist_() {
  // Required for ported functions, since GPUs cannot access the global variable globalTimeStep
  const MInt gTSmOffset = globalTimeStep - timeStepOffset; // Because collision step is performed one TS earlier
  const MInt maxLevel_ = maxLevel();

  if(this->isCompressible()) {
    maia::parallelFor<true>(0, m_currentMaxNoCells, [=](MInt index) {
      const MInt pCellId = m_activeCellList[index];
      const MInt lvlDiff = maxLevel_ - a_level(pCellId);
      if(gTSmOffset % IPOW2(lvlDiff) != 0) return;
      swap_variables(pCellId);
#ifdef WAR_NVHPC_PSTL
      MFloat u[nDim] = {F0};
      for(MInt d = 0; d < nDim; d++) {
        u[d] = a_variable(pCellId, d);
      }
      MFloat l_rho = a_variable(pCellId, PV->RHO);
      calculateMacroscopicVariables<true>(pCellId, l_rho, u);
      for(MInt d = 0; d < nDim; d++) {
        a_variable(pCellId, d) = u[d];
      }
      a_variable(pCellId, PV->RHO) = l_rho;
#else
    calculateMacroscopicVariables<true>(pCellId, a_variable(pCellId, PV->RHO), &a_variable(pCellId, PV->U));
#endif
    });
  } else {
    maia::parallelFor<true>(0, m_currentMaxNoCells, [=](MInt index) {
      const MInt pCellId = m_activeCellList[index];
      const MInt lvlDiff = maxLevel_ - a_level(pCellId);
      if(gTSmOffset % IPOW2(lvlDiff) != 0) return;
      swap_variables(pCellId);
#ifdef WAR_NVHPC_PSTL
      MFloat u[nDim] = {F0};
      for(MInt d = 0; d < nDim; d++) {
        u[d] = a_variable(pCellId, d);
      }
      MFloat l_rho = a_variable(pCellId, PV->RHO);
      calculateMacroscopicVariables(pCellId, l_rho, u);
      for(MInt d = 0; d < nDim; d++) {
        a_variable(pCellId, d) = u[d];
      }
      a_variable(pCellId, PV->RHO) = l_rho;
#else
    calculateMacroscopicVariables(pCellId, a_variable(pCellId, PV->RHO), &a_variable(pCellId, PV->U));
#endif
    });
  }
}

template <MInt nDim, MInt nDist, class SysEqn>
void LbSolverDxQy<nDim, nDist, SysEqn>::updateVariablesFromOldDist() {
  updateVariablesFromOldDist_<0>();
}

template <MInt nDim, MInt nDist, class SysEqn>
void LbSolverDxQy<nDim, nDist, SysEqn>::updateVariablesFromOldDist_preCollision() {
  // TODO labels:LB @johannes/julian: You made this step after propagation. Preference against here?
  // if(this->m_isInitRun) {
  //   this->initRunCorrection();
  // } else {
  updateVariablesFromOldDist_<1>();
  // }
}

/** \brief apply volumeForces to the oldDistributions
 * \author Daniel Lauwers
 * \date June 2021
 *
 * apply volumeForces to the oldDistributions. Should be called after propagation.
 *
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbSolverDxQy<nDim, nDist, SysEqn>::volumeForces() {
  if(string2enum(solverMethod()) == MAIA_LATTICE_BGK_TOTALENERGY)
    return; // in this collision step the forcing is added differently. There seems to be no easy way to unify the
            // methods.
  if(m_controlVelocity) controlVelocity();

  // Required for ported functions, since GPUs cannot access the global variable globalTimeStep
  const MInt gTS = globalTimeStep;
  const MInt maxLevel_ = maxLevel();

  // TODO labels:LB Move to lb source terms?
  if(this->m_externalForcing) {
    // add forcing term
    maia::parallelFor<true>(0, m_currentMaxNoCells, [=](MInt index) {
      const MInt pCellId = m_activeCellList[index];
      const MInt lvlDiff = maxLevel_ - a_level(pCellId);
      if((gTS) % IPOW2(lvlDiff) != 0) return;
      for(MInt j = 0; j < nDist - 1; j++) {
        a_oldDistribution(pCellId, j) += FPOW2(maxLevel_ - this->a_level(pCellId)) * m_Fext[j];
      }
    });
  }

  if(m_isEELiquid && m_EELiquid.gravity) {
    maia::parallelFor<true>(0, m_currentMaxNoCells, [=](MInt index) {
      const MInt pCellId = m_activeCellList[index];
      const MInt lvlDiff = maxLevel_ - a_level(pCellId);
      if((gTS) % IPOW2(lvlDiff) != 0) return;
      const MFloat alpha = a_alphaGasLim(pCellId);

      for(MInt j = 0; j < nDist - 1; j++) {
        // gravity included in m_EELiquid.Fg
        // *3.0 is 1/c_s^2 see Haenel_MolekulareGasdynamik Eq. (9.65)
        // -= because the effects of gravity are introduced through buoyancy
        a_oldDistribution(pCellId, j) -= FPOW2(maxLevel_ - this->a_level(pCellId)) * m_EELiquid.Fg[j] * alpha;
      }
    });
  }
}

/** \brief control velocity of a periodic channel flow via volume forces
 * \author Daniel Lauwers
 * \date June 2021
 *
 * This function determines the necessary volume forces to keep periodic channel flow at a bulk velocity of u_inf.
 *
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbSolverDxQy<nDim, nDist, SysEqn>::controlVelocity() {
  if(globalTimeStep % m_velocityControl.interval != 0) return;

  averageGlobalVelocity(
      m_velocityControl.dir); // m_velocityControl.dir defines the component of the velocity that is averaged
  const MFloat velocityGoal = m_Ma / F1BCS;
  if((m_velocityControl.lastGlobalAvgV > 5.0 * velocityGoal) || (m_velocityControl.lastGlobalAvgV < 0.0)
     || !std::isfinite(m_velocityControl.lastGlobalAvgV)) {
    if(domainId() == 0) cerr << "TS " << globalTimeStep << " avVel " << m_velocityControl.lastGlobalAvgV << endl;
    m_log << "TS " << globalTimeStep << " avVel " << m_velocityControl.lastGlobalAvgV << endl;
    TERMM(1, "Velocity control failed!");
  }
  const MFloat error = (m_velocityControl.lastGlobalAvgV - velocityGoal) / velocityGoal;
  m_velocityControl.integratedError += (error + m_velocityControl.previousError) / 2.0 * m_velocityControl.interval;
  m_velocityControl.derivedError = (error - m_velocityControl.previousError) / m_velocityControl.interval;

  MFloat controlSignal = m_velocityControl.KT
                         * (error + F1 / m_velocityControl.KI * m_velocityControl.integratedError
                            + m_velocityControl.KD * m_velocityControl.derivedError);

  controlSignal = mMax(mMin(controlSignal, 100.0), -100.0);

  m_volumeAccel[m_velocityControl.dir] =
      m_volumeAccelBase[m_velocityControl.dir] - controlSignal * m_volumeAccelBase[m_velocityControl.dir];

  this->initVolumeForces();

  stringstream ss;
  ss << "Average velocity at TS " << globalTimeStep << " is " << m_velocityControl.lastGlobalAvgV
     << " (goal: " << velocityGoal << ") accel set to: [" << m_volumeAccel[0] << ", " << m_volumeAccel[1];
  if(nDim == 3) {
    ss << ", " << m_volumeAccel[2];
  }
  ss << "] (base was [" << m_volumeAccelBase[0] << ", " << m_volumeAccelBase[1];
  if(nDim == 3) {
    ss << ", " << m_volumeAccelBase[2];
  }
  ss << "])" << endl;
  ss << "Error: " << error << "  previousError: " << m_velocityControl.previousError
     << "  integratedError: " << m_velocityControl.integratedError
     << "  derivedError: " << m_velocityControl.derivedError << endl;
  if(domainId() == 0) {
    cerr << ss.str();
  }
  m_log << ss.str();

  m_velocityControl.previousError = error;
}

/** \brief calculate average velocity of velocity component dir
 * \author Daniel Lauwers
 * \date June 2021
 *
 * This function calculates the average velocity of all leaf cells (weighted by their volume) thoughout all domains.
 * The value is saved to m_velocityControl.lastGlobalAvgV
 *
 * \param[in] dir the velocity component to average
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbSolverDxQy<nDim, nDist, SysEqn>::averageGlobalVelocity(const MInt dir) {
  MFloat sumV[2] = {}; // Sum of {velocity * cellVolume, cellVolume}
  MFloat c1 = 0.0;

  // the velocity component of sumV is summed by Kahan summation because of potential for different orders of magnitude
  // for the cell volume normal summation is used
  for(MInt cellId = 0; cellId < this->grid().noInternalCells(); cellId++) {
    if(this->c_noChildren(cellId) != 0) continue;
    const MFloat cellVolume = this->grid().cellVolumeAtLevel(this->c_level(cellId));

    const MFloat y1 = a_variable(cellId, dir) * cellVolume - c1;
    const MFloat t1 = sumV[0] + y1;
    c1 = (t1 - sumV[0]) - y1;
    sumV[0] = t1;
    sumV[1] += cellVolume;
  }
  MPI_Allreduce(MPI_IN_PLACE, &sumV, 2, maia::type_traits<MFloat>::mpiType(), MPI_SUM, mpiComm(), AT_, "MPI_IN_PLACE",
                "sumV");

  m_velocityControl.lastGlobalAvgV = sumV[0] / sumV[1];
}

template <MInt nDim, MInt nDist, class SysEqn>
template <MBool useSmagorinsky>
void LbSolverDxQy<nDim, nDist, SysEqn>::clb_collision_step_base() {
  TRACE();
  std::stringstream ss;
  if constexpr(useSmagorinsky) {
    ss << "CLB_SMAGORINSKY collision step not defined for d" << nDim << "q" << nDist << " !";
  } else {
    ss << "CLB collision step not defined for d" << nDim << "q" << nDist << " !";
  }
  TERMM(1, ss.str());
}

template <>
template <>
void LbSolverDxQy<2, 9, maia::lb::LbSysEqnIncompressible<2, 9>>::clb_collision_step_base<false>() {
  TRACE();
  constexpr MInt nDist = 9;

  constexpr MLongFloat K[9][9] = {
      {1.0, 0.0, 0.0, -4.0, 0.0, 0.0, 0.0, 0.0, 4.0},     {1.0, -1.0, 1.0, 2.0, 0.0, 1.0, -1.0, 1.0, 1.0},
      {1.0, -1.0, 0.0, -1.0, 1.0, 0.0, 0.0, -2.0, -2.0},  {1.0, -1.0, -1.0, 2.0, 0.0, -1.0, 1.0, 1.0, 1.0},
      {1.0, 0.0, -1.0, -1.0, -1.0, 0.0, -2.0, 0.0, -2.0}, {1.0, 1.0, -1.0, 2.0, 0.0, 1.0, 1.0, -1.0, 1.0},
      {1.0, 1.0, 0.0, -1.0, 1.0, 0.0, 0.0, 2.0, -2.0},    {1.0, 1.0, 1.0, 2.0, 0.0, -1.0, -1.0, -1.0, 1.0},
      {1.0, 0.0, 1.0, -1.0, -1.0, 0.0, 2.0, 0.0, -2.0}};

  for(MInt id = 0; id < a_noCells(); id++) {
    if((globalTimeStep - 1) % IPOW2(maxLevel() - this->a_level(id)) == 0) {
      m_nu = m_Ma * LBCS / m_Re * m_referenceLength * FFPOW2(maxLevel() - this->a_level(id));

      MLongFloat omega[9] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
      m_omega = 2.0 / (1.0 + 6.0 * m_nu);
      omega[4] = m_omega;
      omega[5] = m_omega;

      swap_variables(id);

      // Calculation of macroscopic variables + forcing term
      for(MInt j = 0; j < nDist; j++) {
        a_variable(id, PV->RHO) += a_oldDistribution(id, j);
      }
      const MFloat rho = a_variable(id, PV->RHO);

      a_variable(id, PV->U) = (a_oldDistribution(id, 1) + a_oldDistribution(id, 4) + a_oldDistribution(id, 5)
                               - a_oldDistribution(id, 7) - a_oldDistribution(id, 6) - a_oldDistribution(id, 0))
                              / rho;

      a_variable(id, PV->V) = (a_oldDistribution(id, 7) + a_oldDistribution(id, 3) + a_oldDistribution(id, 4)
                               - a_oldDistribution(id, 6) - a_oldDistribution(id, 2) - a_oldDistribution(id, 5))
                              / rho;

      const MFloat u = a_variable(id, PV->U);
      const MFloat v = a_variable(id, PV->V);

      const MFloat r = a_oldDistribution(id, 8);
      const MFloat nw = a_oldDistribution(id, 7);
      const MFloat w = a_oldDistribution(id, 0);
      const MFloat sw = a_oldDistribution(id, 6);
      const MFloat s = a_oldDistribution(id, 2);
      const MFloat se = a_oldDistribution(id, 5);
      const MFloat e = a_oldDistribution(id, 1);
      const MFloat ne = a_oldDistribution(id, 4);
      const MFloat n = a_oldDistribution(id, 3);

      // Calculate equilibrium moments...
      MLongFloat keq[9];
      keq[0] = 0.0;
      keq[1] = 0.0;
      keq[2] = 0.0;
      keq[3] = omega[3] * (rho * (u * u + v * v) - e - n - s - w - 2.0 * (se + sw + ne + nw - F1B3 * rho)) / 12.0;

      keq[4] = omega[4] * (n + s - e - w + rho * (u * u - v * v)) / 4.0;

      keq[5] = omega[5] * ((ne + sw - nw - se) - u * v * rho) / 4.0;

      keq[6] = omega[6]
                   * (-1.0
                      * ((se + sw - ne - nw - 2 * u * u * v * rho + v * (rho - n - s - r)) / 4.0
                         + u / 2.0 * (ne - nw - se + sw)))
               - v / 2.0 * (-3.0 * keq[3] - keq[4]) - 2.0 * u * keq[5];
      keq[7] = omega[7]
                   * (-1.0
                      * ((sw + nw - se - ne - 2 * v * v * u * rho + u * (rho - w - e - r)) / 4.0
                         + v / 2.0 * (ne - nw - se + sw) / 2.0))
               - u / 2.0 * (-3.0 * keq[3] + keq[4]) - 2.0 * v * keq[5];
      keq[8] = omega[8]
                   * (F1B4
                      * (F1B9 * rho - ne - nw - se - sw + 2.0 * (u * (ne - nw + se - sw) + v * (ne + nw - se - sw))
                         + 4.0 * u * v * (nw - ne + se - sw) - u * u * (n + ne + nw + s + se + sw)
                         + v * v * (3.0 * u * u * rho - e - ne - nw - se - sw - w)))
               - 2.0 * keq[3] - 2.0 * u * keq[7] - 2.0 * v * keq[6] + 4 * u * v * keq[5]
               - (1.5 * keq[3] - F1B2 * keq[4]) * (u * u + v * v);

      // Matrix vector multiplication (must be optimized)
      MLongFloat k[9];
      for(MInt l = 0; l < 9; l++) {
        k[l] = std::inner_product(&K[l][0], &K[l][9], &keq[0], F0);
      }

      a_distribution(id, 8) = a_oldDistribution(id, 8) + k[0];
      a_distribution(id, 7) = a_oldDistribution(id, 7) + k[1];
      a_distribution(id, 0) = a_oldDistribution(id, 0) + k[2];
      a_distribution(id, 6) = a_oldDistribution(id, 6) + k[3];
      a_distribution(id, 2) = a_oldDistribution(id, 2) + k[4];
      a_distribution(id, 5) = a_oldDistribution(id, 5) + k[5];
      a_distribution(id, 1) = a_oldDistribution(id, 1) + k[6];
      a_distribution(id, 4) = a_oldDistribution(id, 4) + k[7];
      a_distribution(id, 3) = a_oldDistribution(id, 3) + k[8];
    }
  }
}

template <>
template <>
void LbSolverDxQy<3, 19, maia::lb::LbSysEqnIncompressible<3, 19>>::clb_collision_step_base<false>() {
  TRACE();
  constexpr MInt nDist = 19;

  for(MInt i = 0; i < m_currentMaxNoCells; i++) {
    const MInt pCellId = m_activeCellList[i];

    if((globalTimeStep - 1) % IPOW2(maxLevel() - this->a_level(pCellId)) == 0) {
      m_nu = m_Ma * LBCS / m_Re * m_referenceLength;

      m_omega = 2.0 / (1.0 + 6.0 * m_nu * FFPOW2(maxLevel() - this->a_level(pCellId)));
      MFloat omega = m_omega;

      swap_variables(pCellId);

      // Calculate macroscopic variables
      a_variable(pCellId, PV->RHO) = 0.0;
      for(MInt j = 0; j < nDist; j++) {
        a_variable(pCellId, PV->RHO) += a_oldDistribution(pCellId, j);
      }

      a_variable(pCellId, PV->U) = 0.0;
      for(MInt j = 0; j < Ld::dxQyFld(); j++) {
        a_variable(pCellId, PV->U) += a_oldDistribution(pCellId, Ld::pFld(0, j));
        a_variable(pCellId, PV->U) -= a_oldDistribution(pCellId, Ld::nFld(0, j));
      }
      a_variable(pCellId, PV->U) /= a_variable(pCellId, PV->RHO);

      a_variable(pCellId, PV->V) = 0.0;
      for(MInt j = 0; j < Ld::dxQyFld(); j++) {
        a_variable(pCellId, PV->V) += a_oldDistribution(pCellId, Ld::pFld(1, j));
        a_variable(pCellId, PV->V) -= a_oldDistribution(pCellId, Ld::nFld(1, j));
      }
      a_variable(pCellId, PV->V) /= a_variable(pCellId, PV->RHO);

      a_variable(pCellId, PV->W) = 0.0;
      for(MInt j = 0; j < Ld::dxQyFld(); j++) {
        a_variable(pCellId, PV->W) += a_oldDistribution(pCellId, Ld::pFld(2, j));
        a_variable(pCellId, PV->W) -= a_oldDistribution(pCellId, Ld::nFld(2, j));
      }
      a_variable(pCellId, PV->W) /= a_variable(pCellId, PV->RHO);

      const MFloat R = a_oldDistribution(pCellId, 18);
      const MFloat Nw = a_oldDistribution(pCellId, 7);
      const MFloat W = a_oldDistribution(pCellId, 0);
      const MFloat Sw = a_oldDistribution(pCellId, 6);
      const MFloat S = a_oldDistribution(pCellId, 2);
      const MFloat Se = a_oldDistribution(pCellId, 8);
      const MFloat E = a_oldDistribution(pCellId, 1);
      const MFloat Ne = a_oldDistribution(pCellId, 9);
      const MFloat N = a_oldDistribution(pCellId, 3);
      const MFloat Nf = a_oldDistribution(pCellId, 17);
      const MFloat Nb = a_oldDistribution(pCellId, 16);
      const MFloat Sf = a_oldDistribution(pCellId, 15);
      const MFloat Sb = a_oldDistribution(pCellId, 14);
      const MFloat Ef = a_oldDistribution(pCellId, 13);
      const MFloat Eb = a_oldDistribution(pCellId, 12);
      const MFloat Wf = a_oldDistribution(pCellId, 11);
      const MFloat Wb = a_oldDistribution(pCellId, 10);
      const MFloat F = a_oldDistribution(pCellId, 5);
      const MFloat B = a_oldDistribution(pCellId, 4);

      const MFloat rho = Nw + W + Sw + S + Se + E + Ne + N + R + Nf + Nb + Sf + Sb + Ef + Eb + Wf + Wb + F + B;
      const MFloat pi_x = (Ne + E + Se + Ef + Eb - Nw - W - Sw - Wf - Wb);
      const MFloat pi_y = (Ne + N + Nw + Nf + Nb - Se - S - Sw - Sf - Sb);
      const MFloat pi_z = (Nf + Sf + Wf + Ef + F - Nb - Sb - Wb - Eb - B);
      const MFloat vv_x = pi_x / rho;
      const MFloat vv_y = pi_y / rho;
      const MFloat vv_z = pi_z / rho;
      const MFloat vx2 = vv_x * vv_x;
      const MFloat vy2 = vv_y * vv_y;
      const MFloat vz2 = vv_z * vv_z;
      const MFloat uxy = (omega * (Ne - Nw - Se + Sw + vv_x * vv_y * rho - pi_x * vv_y - pi_y * vv_x) * 0.25);
      const MFloat uxz = (omega * (-Eb + Ef + Wb - Wf + vv_x * vv_z * rho - pi_x * vv_z - pi_z * vv_x) * 0.25);
      const MFloat uyz = (omega * (-Nb + Nf + Sb - Sf + vv_z * vv_y * rho - pi_z * vv_y - pi_y * vv_z) * 0.25);

      const MFloat vxy = (omega * (1. / 6.)
                          * (-B - E - F + Nb + Ne + Nf + Nw + Sb + Se + Sf + Sw - W + 2 * (-Eb - Ef + N + S - Wb - Wf)
                             + 2 * (-2 * pi_y * vv_y + pi_x * vv_x + pi_z * vv_z) - rho * (vz2 - 2 * vy2 + vx2)));
      const MFloat vxz = (omega * (1. / 6.)
                          * (-E + Eb + Ef - N + Nb + Nf - S + Sb + Sf - W + Wb + Wf + 2 * (B + F - Ne - Nw - Se - Sw)
                             + 2 * (-2 * pi_z * vv_z + pi_x * vv_x + pi_y * vv_y) - rho * (vx2 + vy2 - 2 * vz2)));
      const MFloat pe = ((1. / 126.)
                         * (-B - E - F - N - S - W + 2 * (-Eb - Ef - Nb - Ne - Nf - Nw - Sb - Se - Sf - Sw - Wb - Wf)
                            + (rho + vv_y * (2 * pi_y - vv_y * rho) + vv_z * (2 * pi_z - vv_z * rho)
                               + vv_x * (2 * pi_x - vv_x * rho))));

      const MFloat UXY = -Ne + Nw + Se - Sw + 4 * uxy;
      const MFloat UXZ = Eb - Ef - Wb + Wf + 4 * uxz;
      const MFloat UYZ = Nb - Nf - Sb + Sf + 4 * uyz;

      const MFloat x = ((0.125
                             * (Nw + Sw + Wf + Wb - Eb - Ef - Ne - Se
                                + vv_x
                                      * (B + Eb + Ef + 84 * pe + F + N + Ne + Nw + S + Se + Sw + Wb + Wf
                                         + 2 * (Nb + Nf + Sb + Sf - vxy - vxz)))
                         + 0.25 * (-vv_y * UXY - vv_z * UXZ) - 0.25 * vv_x * (vv_y * pi_y + vv_z * pi_z)
                         + 0.125 * (vy2 + vz2) * (vv_x * rho - pi_x)));
      const MFloat y = ((0.125
                             * (Sb + Sf + Se + Sw - Nw - Ne - Nf - Nb
                                + vv_y
                                      * (B + E + 84 * pe + F + Nb + Ne + Nf + Nw + Sb + Se + Sf + Sw + W
                                         + 2 * (Eb + Ef + Wb + Wf + vxy)))
                         + 0.25 * (-vv_x * UXY - vv_z * UYZ) - 0.25 * vv_y * (vv_z * pi_z + vv_x * pi_x)
                         + 0.125 * (vx2 + vz2) * (vv_y * rho - pi_y)));
      const MFloat z = ((0.125
                             * (Eb + Nb + Sb + Wb - Ef - Nf - Sf - Wf
                                + vv_z
                                      * (E + Eb + Ef + N + Nb + Nf + S + Sb + Sf + W + Wb + Wf + 84 * pe
                                         + 2 * (Ne + Nw + Se + Sw + vxz)))
                         + 0.25 * (-vv_x * UXZ - vv_y * UYZ) - 0.25 * vv_z * (vv_x * pi_x + vv_y * pi_y)
                         + 0.125 * (vx2 + vy2) * (vv_z * rho - pi_z)));

      const MFloat xxx = ((0.125
                               * (Eb + Ef - Ne + Nw - Se + Sw - Wb - Wf
                                  + vv_x * (-B - Eb - Ef - F + N + Ne + Nw + S + Se + Sw - Wb - Wf + 2 * (vxz - vxy)))
                           + 0.25 * (vv_z * UXZ - vv_y * UXY + vv_x * (vv_z * pi_z - vv_y * pi_y))
                           + 0.125 * ((pi_x - vv_x * rho) * (vz2 - vy2))));
      const MFloat yyy =
          ((0.125
                * (Ne - Nb - Nf + Nw + Sb - Se + Sf - Sw
                   + vv_y * (-E + B + F + Nb - Ne + Nf - Nw + Sb - Se + Sf - Sw - W + 2 * (-vxy - 2 * vxz)))
            + 0.25 * (vv_x * UXY - vv_z * UYZ + vv_y * (vv_x * pi_x - vv_z * pi_z))
            + 0.125 * ((pi_y - vv_y * rho) * (vx2 - vz2))));
      const MFloat zzz =
          ((0.125
                * (Eb - Ef - Nb + Nf - Sb + Sf + Wb - Wf
                   + vv_z * (E + Eb + Ef - N - Nb - Nf - S - Sb - Sf + W + Wb + Wf + 2 * (vxz + 2 * vxy)))
            + 0.25 * (vv_y * UYZ - vv_x * UXZ + vv_z * (vv_y * pi_y - vv_x * pi_x))
            + 0.125 * ((pi_z - vv_z * rho) * (vy2 - vx2))));

      // xyz=(0.125*(Neb-Nef-Nwb+Nwf-Seb+Sef+Swb-Swf-vv_x*UYZ-vv_y*UXZ-vv_z*UXY+rho*vv_x*vv_y*vv_z-(vv_x*vv_y*pi_z+vv_x*pi_y*vv_z+pi_x*vv_y*vv_z)));

      const MFloat X_ = (Eb + Ef + Ne - Nw + Se - Sw - Wb - Wf + 8 * x);
      const MFloat Y_ = (Nb + Ne + Nf + Nw - Sb - Se - Sf - Sw + 8 * y);
      const MFloat Z_ = (Ef - Eb - Nb + Nf - Sb + Sf - Wb + Wf + 8 * z);

      // XnXXX=Eb+Ef-Wb-Wf+4*(x-xxx);
      // XpXXX=Ne-Nw+Se-Sw+4*(x+xxx);
      // YnYYY=Ne+Nw-Se-Sw+4*(y-yyy);
      // YpYYY=Nb+Nf-Sb-Sf+4*(y+yyy);
      // ZnZZZ=-Nb+Nf-Sb+Sf+4*(z-zzz);
      // ZpZZZ=-Eb+Ef-Wb+Wf+4*(z+zzz);

      // XYZ= 8*xyz;
      // EPa=-42*pe-E-Eb-Ef-Ne-Nw-Se-Sw-W-Wb-Wf-2*(vxy+vxz);
      // EPb=-42*pe-B-Eb-Ef-F-Nb-Nf-Sb-Sf-Wb-Wf+2*vxz;
      // EPc=-42*pe-N-Nb-Ne-Nf-Nw-S-Sb-Se-Sf-Sw+2*vxy;

      const MFloat p = ((1. / 12.)
                        * (rho / (3.) - 96 * pe - Eb - Ef - Nb - Ne - Nf - Nw - Sb - Se - Sf - Sw - Wb - Wf
                           + 2 * (+vv_x * X_ + vv_y * Y_ + vv_z * Z_)
                           + vx2
                                 * (-84 * pe - B - Eb - Ef - F - N - Ne - Nw - S - Se - Sw - Wb - Wf
                                    + 2 * (vxy + vxz - Nb - Nf - Sb - Sf))
                           + vy2
                                 * (-84 * pe - B - E - F - Nb - Ne - Nf - Nw - Sb - Se - Sf - Sw - W
                                    + 2 * (-vxy - Eb - Ef - Wb - Wf))
                           + vz2
                                 * (-84 * pe - E - Eb - Ef - N - Nb - Nf - S - Sb - Sf - W - Wb - Wf
                                    + 2 * (-vxz - Ne - Nw - Se - Sw))
                           + 4 * (vv_x * vv_y * UXY + vv_x * vv_z * UXZ + vv_y * vv_z * UYZ)
                           - rho * (vx2 * vy2 + vx2 * vz2 + vy2 * vz2)
                           + 2 * (pi_y * vv_y * (vx2 + vz2) + pi_x * vv_x * (vy2 + vz2) + pi_z * vv_z * (vx2 + vy2))));

      const MFloat a =
          ((1. / 12.)
           * (-Eb - Ef - Ne - Nw - Se - Sw - Wb - Wf
              + 2
                    * (Nb + Nf + Sb + Sf + vv_x * X_
                       + vv_y * (Ne + Nw - Se - Sw + 2 * (Sf + Sb - Nb - Nf - 2 * y - 6 * yyy))
                       + vv_z * (Ef - Eb - Wb + Wf + 2 * (Nb - Nf + Sb - Sf - 2 * z + 6 * zzz)))
              + vx2
                    * (-B - Eb - Ef - F - N - Ne - Nw - S - Se - Sw - Wb - Wf - 84 * pe
                       - 2 * (Nb + Nf + Sb + Sf - vxy - vxz))
              + vy2
                    * (Eb - E + Ef - Ne - Nw - Se - Sw - W + Wb + Wf + 42 * pe
                       + 2 * (B + F + Nb + Nf + Sb + Sf - vxy - 3 * vxz))
              + vz2
                    * (Ne - E - Eb - Ef + Nw + Se + Sw - W - Wb - Wf + 42 * pe
                       + 2 * (N + Nb + Nf + S + Sb + Sf - vxz - 3 * vxy))
              + 4 * (vv_x * (vv_y * UXY + vv_z * UXZ)) - 8 * vv_y * vv_z * UYZ
              + (2 * vy2 * vz2 - vx2 * (vy2 + vz2)) * rho
              + 2 * (pi_y * vv_y * (vx2 - 2 * vz2) + pi_z * vv_z * (vx2 - 2 * vy2) + pi_x * vv_x * (vy2 + vz2))));
      const MFloat c =
          ((1. / 12.)
           * (-Nb - Ne - Nf - Nw - Sb - Se - Sf - Sw
              + 2
                    * (Eb + Ef + Wb + Wf + vv_x * (Ne - Nw + Se - Sw + 2 * (Wb + Wf - Eb - Ef - 2 * x + 6 * xxx))
                       + vv_y * Y_ + vv_z * (Nf - Nb - Sb + Sf + 2 * (Eb - Ef + Wb - Wf - 2 * z - 6 * zzz)))
              + vx2
                    * (Nb - N - Ne + Nf - Nw - S + Sb - Se + Sf - Sw + 42 * pe
                       + 2 * (B + Eb + Ef + F + Wb + Wf + vxy - 2 * vxz))
              - vy2 * (B + E + F + Nb + Ne + Nf + Nw + Sb + Se + Sf + Sw + W + 84 * pe + 2 * (Eb + Ef + Wb + Wf + vxy))
              + vz2
                    * (Ne - N - Nb - Nf + Nw - S - Sb + Se - Sf + Sw + 42 * pe
                       + 2 * (E + Eb + Ef + W + Wb + Wf + 3 * vxy + 2 * vxz))
              + 4 * vv_y * (vv_x * UXY + vv_z * UYZ) - 8 * vv_x * vv_z * UXZ + (2 * vx2 * vz2 - (vx2 + vz2) * vy2) * rho
              + 2 * (pi_x * vv_x * (vy2 - 2 * vz2) + pi_z * vv_z * (vy2 - 2 * vx2) + pi_y * vv_y * (vx2 + vz2))));

      a_distribution(pCellId, 18) = R + (12 * p - 30 * pe);
      a_distribution(pCellId, 7) = Nw + c + a - xxx - yyy - x + y + uxy + p + 8 * pe;
      a_distribution(pCellId, 0) = W - 2 * a + 4 * x - 11 * pe + vxy + vxz - 4 * p;
      a_distribution(pCellId, 6) = Sw + a + c - xxx + yyy - x - y - uxy + p + 8 * pe;
      a_distribution(pCellId, 2) = S - 2 * c + 4 * y - 11 * pe - vxy - 4 * p;
      a_distribution(pCellId, 8) = Se + a + c + xxx + yyy + x - y + uxy + p + 8 * pe;
      a_distribution(pCellId, 1) = E - 2 * a - 4 * x - 11 * pe + vxy + vxz - 4 * p;
      a_distribution(pCellId, 9) = Ne + c + a + x + y + xxx - yyy - uxy + p + 8 * pe;
      a_distribution(pCellId, 3) = N - 2 * c - 4 * y - 11 * pe - vxy - 4 * p; // s_n;
      a_distribution(pCellId, 17) = Nf - a + yyy - zzz + y + z - uyz + p + 8 * pe;
      a_distribution(pCellId, 16) = Nb - a + yyy + zzz + y - z + uyz + p + 8 * pe;
      a_distribution(pCellId, 15) = Sf - a - yyy - zzz - y + z + uyz + p + 8 * pe;
      a_distribution(pCellId, 14) = Sb - a - yyy + zzz - y - z - uyz + p + 8 * pe;
      a_distribution(pCellId, 13) = Ef - c - xxx + zzz + x + z - uxz + p + 8 * pe;
      a_distribution(pCellId, 12) = Eb - c - xxx - zzz + x - z + uxz + p + 8 * pe;
      a_distribution(pCellId, 11) = Wf - c + xxx + zzz - x + z + uxz + p + 8 * pe;
      a_distribution(pCellId, 10) = Wb - c + xxx - zzz - x - z - uxz + p + 8 * pe;
      a_distribution(pCellId, 5) = F + 2 * a + 2 * c - 4 * z - 11 * pe - vxz - 4 * p;
      a_distribution(pCellId, 4) = B + 2 * a + 2 * c + 4 * z - 11 * pe - vxz - 4 * p;
    }
  }
  if(m_calculateDissipation && globalTimeStep % m_solutionInterval == 0) {
    calculateDissipation();
  }
}

template <>
template <MBool useSmagorinsky>
void LbSolverDxQy<3, 27, maia::lb::LbSysEqnIncompressible<3, 27>>::clb_collision_step_base() {
  TRACE();
  constexpr MInt nDist = 27;
  constexpr MInt nDim = 3;
  constexpr MInt nDimSqr = 9;

  [[maybe_unused]] constexpr MFloat C_s = 0.1; /// Smagorinsky constant

  for(MInt i = 0; i < m_currentMaxNoCells; i++) {
    const MInt pCellId = m_activeCellList[i];

    if((globalTimeStep - 1) % IPOW2(maxLevel() - this->a_level(pCellId)) == 0) {
      m_nu = m_Ma * LBCS / m_Re * m_referenceLength;

      m_omega = 2.0 / (1.0 + 6.0 * m_nu * FFPOW2(maxLevel() - this->a_level(pCellId)));
      swap_variables(pCellId);

      // Calculate macroscopic variables
      a_variable(pCellId, PV->RHO) = 0.0;
      for(MInt j = 0; j < nDist; j++) {
        a_variable(pCellId, PV->RHO) += a_oldDistribution(pCellId, j);
      }

      a_variable(pCellId, PV->U) = 0.0;
      for(MInt j = 0; j < Ld::dxQyFld(); j++) {
        a_variable(pCellId, PV->U) += a_oldDistribution(pCellId, Ld::pFld(0, j));
        a_variable(pCellId, PV->U) -= a_oldDistribution(pCellId, Ld::nFld(0, j));
      }
      a_variable(pCellId, PV->U) /= a_variable(pCellId, PV->RHO);

      a_variable(pCellId, PV->V) = 0.0;
      for(MInt j = 0; j < Ld::dxQyFld(); j++) {
        a_variable(pCellId, PV->V) += a_oldDistribution(pCellId, Ld::pFld(1, j));
        a_variable(pCellId, PV->V) -= a_oldDistribution(pCellId, Ld::nFld(1, j));
      }
      a_variable(pCellId, PV->V) /= a_variable(pCellId, PV->RHO);

      a_variable(pCellId, PV->W) = 0.0;
      for(MInt j = 0; j < Ld::dxQyFld(); j++) {
        a_variable(pCellId, PV->W) += a_oldDistribution(pCellId, Ld::pFld(2, j));
        a_variable(pCellId, PV->W) -= a_oldDistribution(pCellId, Ld::nFld(2, j));
      }
      a_variable(pCellId, PV->W) /= a_variable(pCellId, PV->RHO);

      if constexpr(useSmagorinsky) {
        // Calculation of equilibrium and non-equilibrium distribution
        std::array<MFloat, nDist> eqDist{};
        std::array<MFloat, nDist> d{};
        std::array<MFloat, nDim> u{};
        for(MInt dim = 0; dim < nDim; dim++) {
          u[dim] = a_variable(pCellId, dim);
        }
        eqDist = getEqDists(a_variable(pCellId, PV->RHO), u.data());

        for(MInt j = 0; j < nDist; j++) {
          d[j] = a_oldDistribution(pCellId, j) - eqDist[j];
        }

        // Calculation of strain rate tensor (stress tensor)
        std::array<MFloat, nDimSqr> Q{};
        Q[0] = F2B3
                   * (d[0] + d[1] + d[6] + d[7] + d[8] + d[9] + d[10] + d[11] + d[12] + d[13] + d[18] + d[19] + d[20]
                      + d[21] + d[22] + d[23] + d[24] + d[25])
               - F1B3 * (d[2] + d[3] + d[4] + d[5] + d[14] + d[15] + d[16] + d[17] + d[26]);
        Q[1] = (d[6] - d[7] - d[8] + d[9] + d[18] + d[19] - d[20] - d[21] - d[22] - d[23] + d[24] + d[25]);
        Q[2] = (d[10] - d[11] - d[12] + d[13] + d[18] - d[19] + d[20] - d[21] - d[22] + d[23] - d[24] + d[25]);
        Q[3] = Q[1];
        Q[4] = F2B3
                   * (d[2] + d[3] + d[6] + d[7] + d[8] + d[9] + d[14] + d[15] + d[16] + d[17] + d[18] + d[19] + d[20]
                      + d[21] + d[22] + d[23] + d[24] + d[25])
               - F1B3 * (d[0] + d[1] + d[4] + d[5] + d[10] + d[11] + d[12] + d[13] + d[26]);
        Q[5] = (d[14] - d[15] - d[16] + d[17] + d[18] - d[19] - d[20] + d[21] + d[22] - d[23] - d[24] + d[25]);
        Q[6] = Q[2];
        Q[7] = Q[5];
        Q[8] = F2B3
                   * (d[4] + d[5] + d[10] + d[11] + d[12] + d[13] + d[14] + d[15] + d[16] + d[17] + d[18] + d[19]
                      + d[20] + d[21] + d[22] + d[23] + d[24] + d[25])
               - F1B3 * (d[0] + d[1] + d[2] + d[3] + d[6] + d[7] + d[8] + d[9] + d[26]);

        const MFloat Q_mean = sqrt(2.0 * std::inner_product(&Q[0], &Q[nDimSqr], &Q[0], .0));

        // eddy viscosity
        const MFloat nu_t = C_s * C_s * Q_mean;

        m_nu = m_Ma * LBCS / m_Re * m_referenceLength;
        m_nu += nu_t; // this has to adapted for refined grids !!!
        m_omega = 2.0 / (1.0 + 6.0 * m_nu * FFPOW2(maxLevel() - this->a_level(pCellId)));
      }

      const MFloat R = a_oldDistribution(pCellId, 26);
      const MFloat Nw = a_oldDistribution(pCellId, 7);
      const MFloat W = a_oldDistribution(pCellId, 0);
      const MFloat Sw = a_oldDistribution(pCellId, 6);
      const MFloat S = a_oldDistribution(pCellId, 2);
      const MFloat Se = a_oldDistribution(pCellId, 8);
      const MFloat E = a_oldDistribution(pCellId, 1);
      const MFloat Ne = a_oldDistribution(pCellId, 9);
      const MFloat N = a_oldDistribution(pCellId, 3);
      const MFloat Nf = a_oldDistribution(pCellId, 17);
      const MFloat Nb = a_oldDistribution(pCellId, 16);
      const MFloat Sf = a_oldDistribution(pCellId, 15);
      const MFloat Sb = a_oldDistribution(pCellId, 14);
      const MFloat Ef = a_oldDistribution(pCellId, 13);
      const MFloat Eb = a_oldDistribution(pCellId, 12);
      const MFloat Wf = a_oldDistribution(pCellId, 11);
      const MFloat Wb = a_oldDistribution(pCellId, 10);
      const MFloat F = a_oldDistribution(pCellId, 5);
      const MFloat B = a_oldDistribution(pCellId, 4);
      const MFloat Nwf = a_oldDistribution(pCellId, 21);
      const MFloat Nwb = a_oldDistribution(pCellId, 20);
      const MFloat Nef = a_oldDistribution(pCellId, 25);
      const MFloat Neb = a_oldDistribution(pCellId, 24);
      const MFloat Swf = a_oldDistribution(pCellId, 19);
      const MFloat Swb = a_oldDistribution(pCellId, 18);
      const MFloat Sef = a_oldDistribution(pCellId, 23);
      const MFloat Seb = a_oldDistribution(pCellId, 22);

      const MFloat rho = Nw + W + Sw + S + Se + E + Ne + N + R + Nf + Nb + Sf + Sb + Ef + Eb + Wf + Wb + Nwf + Nwb + Nef
                         + Neb + Swf + Swb + Sef + Seb + F + B;
      const MFloat pi_x =
          (Ne + E + Se + Ef + Eb - Nw - W - Sw - Wf - Wb + Nef + Neb + Sef + Seb - Nwf - Nwb - Swf - Swb);
      const MFloat pi_y =
          (Ne + N + Nw + Nf + Nb - Se - S - Sw - Sf - Sb + Nef + Neb + Nwf + Nwb - Sef - Seb - Swf - Swb);
      const MFloat pi_z =
          (Nf + Sf + Wf + Ef + F - Nb - Sb - Wb - Eb - B + Nef + Nwf + Sef + Swf - Neb - Nwb - Seb - Swb);
      const MFloat vv_x = pi_x / rho;
      const MFloat vv_y = pi_y / rho;
      const MFloat vv_z = pi_z / rho;
      const MFloat vx2 = vv_x * vv_x;
      const MFloat vy2 = vv_y * vv_y;
      const MFloat vz2 = vv_z * vv_z;
      const MFloat uxy = (m_omega
                          * (Ne + Neb + Nef - Nw - Nwb - Nwf - Se - Seb - Sef + Sw + Swb + Swf + vv_x * vv_y * rho
                             - pi_x * vv_y - pi_y * vv_x)
                          * 0.25);
      const MFloat uxz = (m_omega
                          * (-Eb + Ef - Neb + Nef + Nwb - Nwf - Seb + Sef + Swb - Swf + Wb - Wf + vv_x * vv_z * rho
                             - pi_x * vv_z - pi_z * vv_x)
                          * 0.25);
      const MFloat uyz = (m_omega
                          * (-Nb - Neb + Nef + Nf - Nwb + Nwf + Sb + Seb - Sef - Sf + Swb - Swf + vv_z * vv_y * rho
                             - pi_z * vv_y - pi_y * vv_z)
                          * 0.25);

      const MFloat vxy = (m_omega * (1. / 6.)
                          * (-B - E - F + Nb + Ne + Nf + Nw + Sb + Se + Sf + Sw - W + 2 * (-Eb - Ef + N + S - Wb - Wf)
                             + 2 * (-2 * pi_y * vv_y + pi_x * vv_x + pi_z * vv_z) - rho * (vz2 - 2 * vy2 + vx2)));
      const MFloat vxz = (m_omega * (1. / 6.)
                          * (-E + Eb + Ef - N + Nb + Nf - S + Sb + Sf - W + Wb + Wf + 2 * (B + F - Ne - Nw - Se - Sw)
                             + 2 * (-2 * pi_z * vv_z + pi_x * vv_x + pi_y * vv_y) - rho * (vx2 + vy2 - 2 * vz2)));
      const MFloat pe = ((1. / 126.)
                         * (-B - E - F - N - S - W + 2 * (-Eb - Ef - Nb - Ne - Nf - Nw - Sb - Se - Sf - Sw - Wb - Wf)
                            + 3 * (-Neb - Nef - Nwb - Nwf - Seb - Sef - Swb - Swf)
                            + (rho + vv_y * (2 * pi_y - vv_y * rho) + vv_z * (2 * pi_z - vv_z * rho)
                               + vv_x * (2 * pi_x - vv_x * rho))));

      const MFloat UXY = -Ne - Neb - Nef + Nw + Nwb + Nwf + Se + Seb + Sef - Sw - Swb - Swf + 4 * uxy;
      const MFloat UXZ = Eb - Ef + Neb - Nef - Nwb + Nwf + Seb - Sef - Swb + Swf - Wb + Wf + 4 * uxz;
      const MFloat UYZ = Nb + Neb - Nef - Nf + Nwb - Nwf - Sb - Seb + Sef + Sf - Swb + Swf + 4 * uyz;

      const MFloat x =
          ((0.125
                * (Nw + Sw + Wf + Wb - Eb - Ef - Ne - Se
                   + vv_x
                         * (B + Eb + Ef + 84 * pe + F + N + Ne + Nw + S + Se + Sw + Wb + Wf
                            + 2 * (Nb + Neb + Nef + Nf + Nwb + Nwf + Sb + Seb + Sef + Sf + Swb + Swf - vxy - vxz)))
            + 0.25 * (-Neb - Nef + Nwb + Nwf - Seb - Sef + Swb + Swf - vv_y * UXY - vv_z * UXZ)
            - 0.25 * vv_x * (vv_y * pi_y + vv_z * pi_z) + 0.125 * (vy2 + vz2) * (vv_x * rho - pi_x)));
      const MFloat y =
          ((0.125
                * (Sb + Sf + Se + Sw - Nw - Ne - Nf - Nb
                   + vv_y
                         * (B + E + 84 * pe + F + Nb + Ne + Nf + Nw + Sb + Se + Sf + Sw + W
                            + 2 * (Eb + Ef + Neb + Nef + Nwb + Nwf + Seb + Sef + Swb + Swf + Wb + Wf + vxy)))
            + 0.25 * (-Neb - Nef - Nwb - Nwf + Seb + Sef + Swb + Swf - vv_x * UXY - vv_z * UYZ)
            - 0.25 * vv_y * (vv_z * pi_z + vv_x * pi_x) + 0.125 * (vx2 + vz2) * (vv_y * rho - pi_y)));
      const MFloat z =
          ((0.125
                * (Eb + Nb + Sb + Wb - Ef - Nf - Sf - Wf
                   + vv_z
                         * (E + Eb + Ef + N + Nb + Nf + S + Sb + Sf + W + Wb + Wf + 84 * pe
                            + 2 * (Ne + Neb + Nef + Nw + Nwb + Nwf + Se + Seb + Sef + Sw + Swb + Swf + vxz)))
            + 0.25 * (-Nef + Neb + Nwb - Nwf + Seb - Sef + Swb - Swf - vv_x * UXZ - vv_y * UYZ)
            - 0.25 * vv_z * (vv_x * pi_x + vv_y * pi_y) + 0.125 * (vx2 + vy2) * (vv_z * rho - pi_z)));

      const MFloat xxx = ((0.125
                               * (Eb + Ef - Ne + Nw - Se + Sw - Wb - Wf
                                  + vv_x * (-B - Eb - Ef - F + N + Ne + Nw + S + Se + Sw - Wb - Wf + 2 * (vxz - vxy)))
                           + 0.25 * (vv_z * UXZ - vv_y * UXY + vv_x * (vv_z * pi_z - vv_y * pi_y))
                           + 0.125 * ((pi_x - vv_x * rho) * (vz2 - vy2))));
      const MFloat yyy =
          ((0.125
                * (Ne - Nb - Nf + Nw + Sb - Se + Sf - Sw
                   + vv_y * (-E + B + F + Nb - Ne + Nf - Nw + Sb - Se + Sf - Sw - W + 2 * (-vxy - 2 * vxz)))
            + 0.25 * (vv_x * UXY - vv_z * UYZ + vv_y * (vv_x * pi_x - vv_z * pi_z))
            + 0.125 * ((pi_y - vv_y * rho) * (vx2 - vz2))));
      const MFloat zzz =
          ((0.125
                * (Eb - Ef - Nb + Nf - Sb + Sf + Wb - Wf
                   + vv_z * (E + Eb + Ef - N - Nb - Nf - S - Sb - Sf + W + Wb + Wf + 2 * (vxz + 2 * vxy)))
            + 0.25 * (vv_y * UYZ - vv_x * UXZ + vv_z * (vv_y * pi_y - vv_x * pi_x))
            + 0.125 * ((pi_z - vv_z * rho) * (vy2 - vx2))));

      const MFloat xyz =
          (0.125
           * (Neb - Nef - Nwb + Nwf - Seb + Sef + Swb - Swf - vv_x * UYZ - vv_y * UXZ - vv_z * UXY
              + rho * vv_x * vv_y * vv_z - (vv_x * vv_y * pi_z + vv_x * pi_y * vv_z + pi_x * vv_y * vv_z)));

      const MFloat X_ =
          (Eb + Ef + Ne - Nw + Se - Sw - Wb - Wf + 2 * (Neb + Nef - Nwb - Nwf + Seb + Sef - Swb - Swf) + 8 * x);
      const MFloat Y_ =
          (Nb + Ne + Nf + Nw - Sb - Se - Sf - Sw + 2 * (Neb + Nef + Nwb + Nwf - Seb - Sef - Swb - Swf) + 8 * y);
      const MFloat Z_ =
          (Ef - Eb - Nb + Nf - Sb + Sf - Wb + Wf + 2 * (Nef - Neb - Nwb + Nwf - Seb + Sef - Swb + Swf) + 8 * z);

      const MFloat XnXXX = Eb + Ef + Neb + Nef - Nwb - Nwf + Seb + Sef - Swb - Swf - Wb - Wf + 4 * (x - xxx);
      const MFloat XpXXX = Ne + Neb + Nef - Nw - Nwb - Nwf + Se + Seb + Sef - Sw - Swb - Swf + 4 * (x + xxx);
      const MFloat YnYYY = Ne + Neb + Nef + Nw + Nwb + Nwf - Se - Seb - Sef - Sw - Swb - Swf + 4 * (y - yyy);
      const MFloat YpYYY = Nb + Neb + Nef + Nf + Nwb + Nwf - Sb - Seb - Sef - Sf - Swb - Swf + 4 * (y + yyy);
      const MFloat ZnZZZ = Nef - Nb - Neb + Nf - Nwb + Nwf - Sb - Seb + Sef + Sf - Swb + Swf + 4 * (z - zzz);
      const MFloat ZpZZZ = Nef - Eb + Ef - Neb - Nwb + Nwf - Seb + Sef - Swb + Swf - Wb + Wf + 4 * (z + zzz);

      const MFloat XYZ = -Neb + Nef + Nwb - Nwf + Seb - Sef - Swb + Swf + 8 * xyz;
      const MFloat EPa = -42 * pe - E - Eb - Ef - Ne - Neb - Nef - Nw - Nwb - Nwf - Se - Seb - Sef - Sw - Swb - Swf - W
                         - Wb - Wf - 2 * (vxy + vxz);
      const MFloat EPb = -42 * pe - B - Eb - Ef - F - Nb - Neb - Nef - Nf - Nwb - Nwf - Sb - Seb - Sef - Sf - Swb - Swf
                         - Wb - Wf + 2 * vxz;
      const MFloat EPc = -42 * pe - N - Nb - Ne - Neb - Nef - Nf - Nw - Nwb - Nwf - S - Sb - Se - Seb - Sef - Sf - Sw
                         - Swb - Swf + 2 * vxy;

      const MFloat p =
          ((1. / 12.)
           * (rho / (3.) - 96 * pe - Eb - Ef - Nb - Ne - Nf - Nw - Sb - Se - Sf - Sw - Wb - Wf
              - 3 * (Nwf + Nwb + Nef + Neb + Swf + Swb + Sef + Seb) + 2 * (+vv_x * X_ + vv_y * Y_ + vv_z * Z_)
              + vx2
                    * (-84 * pe - B - Eb - Ef - F - N - Ne - Nw - S - Se - Sw - Wb - Wf
                       + 2 * (vxy + vxz - Nb - Neb - Nef - Nf - Nwb - Nwf - Sb - Seb - Sef - Sf - Swb - Swf))
              + vy2
                    * (-84 * pe - B - E - F - Nb - Ne - Nf - Nw - Sb - Se - Sf - Sw - W
                       + 2 * (-vxy - Eb - Ef - Neb - Nef - Nwb - Nwf - Seb - Sef - Swb - Swf - Wb - Wf))
              + vz2
                    * (-84 * pe - E - Eb - Ef - N - Nb - Nf - S - Sb - Sf - W - Wb - Wf
                       + 2 * (-vxz - Ne - Neb - Nef - Nw - Nwb - Nwf - Se - Seb - Sef - Sw - Swb - Swf))
              + 4 * (vv_x * vv_y * UXY + vv_x * vv_z * UXZ + vv_y * vv_z * UYZ)
              - rho * (vx2 * vy2 + vx2 * vz2 + vy2 * vz2)
              + 2 * (pi_y * vv_y * (vx2 + vz2) + pi_x * vv_x * (vy2 + vz2) + pi_z * vv_z * (vx2 + vy2))));
      const MFloat a =
          ((1. / 12.)
           * (-Eb - Ef - Ne - Nw - Se - Sw - Wb - Wf
              + 2
                    * (Nb + Nf + Sb + Sf + vv_x * X_
                       + vv_y
                             * (Ne - Neb - Nef + Nw - Nwb - Nwf - Se + Seb + Sef - Sw + Swb + Swf
                                + 2 * (Sf + Sb - Nb - Nf - 2 * y - 6 * yyy))
                       + vv_z
                             * (Ef - Eb + Neb - Nef + Nwb - Nwf + Seb - Sef + Swb - Swf - Wb + Wf
                                + 2 * (Nb - Nf + Sb - Sf - 2 * z + 6 * zzz)))
              + vx2
                    * (-B - Eb - Ef - F - N - Ne - Nw - S - Se - Sw - Wb - Wf - 84 * pe
                       - 2 * (Nb + Neb + Nef + Nf + Nwb + Nwf + Sb + Seb + Sef + Sf + Swb + Swf - vxy - vxz))
              + vy2
                    * (Eb - E + Ef - Ne + Neb + Nef - Nw + Nwb + Nwf - Se + Seb + Sef - Sw + Swb + Swf - W + Wb + Wf
                       + 42 * pe + 2 * (B + F + Nb + Nf + Sb + Sf - vxy - 3 * vxz))
              + vz2
                    * (Ne - E - Eb - Ef + Neb + Nef + Nw + Nwb + Nwf + Se + Seb + Sef + Sw + Swb + Swf - W - Wb - Wf
                       + 42 * pe + 2 * (N + Nb + Nf + S + Sb + Sf - vxz - 3 * vxy))
              + 4 * (vv_x * (vv_y * UXY + vv_z * UXZ)) - 8 * vv_y * vv_z * UYZ
              + (2 * vy2 * vz2 - vx2 * (vy2 + vz2)) * rho
              + 2 * (pi_y * vv_y * (vx2 - 2 * vz2) + pi_z * vv_z * (vx2 - 2 * vy2) + pi_x * vv_x * (vy2 + vz2))));
      const MFloat c =
          ((1. / 12.)
           * (-Nb - Ne - Nf - Nw - Sb - Se - Sf - Sw
              + 2
                    * (Eb + Ef + Wb + Wf
                       + vv_x
                             * (Ne - Neb - Nef - Nw + Nwb + Nwf + Se - Seb - Sef - Sw + Swb + Swf
                                + 2 * (Wb + Wf - Eb - Ef - 2 * x + 6 * xxx))
                       + vv_y * Y_
                       + vv_z
                             * (Nf - Nb + Neb - Nef + Nwb - Nwf - Sb + Seb - Sef + Sf + Swb - Swf
                                + 2 * (Eb - Ef + Wb - Wf - 2 * z - 6 * zzz)))
              + vx2
                    * (Nb - N - Ne + Neb + Nef + Nf - Nw + Nwb + Nwf - S + Sb - Se + Seb + Sef + Sf - Sw + Swb + Swf
                       + 42 * pe + 2 * (B + Eb + Ef + F + Wb + Wf + vxy - 2 * vxz))
              - vy2
                    * (B + E + F + Nb + Ne + Nf + Nw + Sb + Se + Sf + Sw + W + 84 * pe
                       + 2 * (Eb + Ef + Neb + Nef + Nwb + Nwf + Seb + Sef + Swb + Swf + Wb + Wf + vxy))
              + vz2
                    * (Ne - N - Nb + Neb + Nef - Nf + Nw + Nwb + Nwf - S - Sb + Se + Seb + Sef - Sf + Sw + Swb + Swf
                       + 42 * pe + 2 * (E + Eb + Ef + W + Wb + Wf + 3 * vxy + 2 * vxz))
              + 4 * vv_y * (vv_x * UXY + vv_z * UYZ) - 8 * vv_x * vv_z * UXZ + (2 * vx2 * vz2 - (vx2 + vz2) * vy2) * rho
              + 2 * (pi_x * vv_x * (vy2 - 2 * vz2) + pi_z * vv_z * (vy2 - 2 * vx2) + pi_y * vv_y * (vx2 + vz2))));
      const MFloat uxxyz = ((1. / 8.)
                            * (Neb - Nef + Nwb - Nwf - Seb + Sef - Swb + Swf + vv_y * ZpZZZ + vv_z * YnYYY
                               + 2 * vv_x * (XYZ + vv_y * UXZ + vv_z * UXY) + vx2 * UYZ + vv_y * vv_z * EPa
                               + (2 * pi_x - rho * vv_x) * vv_x * vv_y * vv_z + vx2 * (pi_y * vv_z + pi_z * vv_y)));
      const MFloat uxyyz = ((1. / 8.)
                            * (Neb - Nef - Nwb + Nwf + Seb - Sef - Swb + Swf + vv_x * ZnZZZ /*YnYYY*/ + vv_z * XpXXX
                               + 2 * vv_y * (XYZ + vv_x * UYZ + vv_z * UXY) + vy2 * UXZ
                               + vv_x * vv_z * EPc
                               //+rho*3*vy2*vv_x*vv_z
                               + (2 * pi_y - vv_y * rho) * vv_y * vv_x * vv_z + vy2 * (pi_x * vv_z + pi_z * vv_x)));
      const MFloat uxyzz = ((1. / 8.)
                            * (Nwb - Neb - Nef + Nwf + Seb + Sef - Swb - Swf + vv_x * YpYYY + vv_y * XnXXX
                               + 2 * vv_z * (XYZ + vv_y * UXZ + vv_x * UYZ) + vz2 * UXY + vv_x * vv_y * EPb
                               + (2 * pi_z - vv_z * rho) * vv_x * vv_y * vv_z + vz2 * (pi_y * vv_x + pi_x * vv_y)));

      const MFloat A = 4 * a - 32 * pe - Nb - Neb - Nef - Nf - Nwb - Nwf - 4 * p - Sb - Seb - Sef - Sf - Swb - Swf;
      const MFloat C = 4 * (c - p) - Eb - Ef - 32 * pe - Neb - Nef - Nwb - Nwf - Seb - Sef - Swb - Swf - Wb - Wf;
      const MFloat CA = -4 * (a + c + p) - 32 * pe - Ne - Neb - Nef - Nw - Nwb - Nwf - Se - Seb - Sef - Sw - Swb - Swf;
      const MFloat UXYZZ = Nwb - Neb - Nef + Nwf + Seb + Sef - Swb - Swf - 8 * uxyzz;
      const MFloat UXYYZ = Neb - Nef - Nwb + Nwf + Seb - Sef - Swb + Swf - 8 * uxyyz;
      const MFloat UXXYZ = Neb - Nef + Nwb - Nwf - Seb + Sef - Swb + Swf - 8 * uxxyz;

      const MFloat xyyzz =
          ((1. / 8.)
           * (Nwb - Neb - Nef + Nwf - Seb - Sef + Swb + Swf - vv_x * A - 2 * (vv_y * UXYZZ + vv_z * UXYYZ) - vy2 * XnXXX
              - vz2 * XpXXX + 2 * vv_x * (-vv_y * YpYYY - vv_z * ZnZZZ) - 4 * vv_y * vv_z * (XYZ + vv_x * UYZ)
              - vv_x * (vz2 * EPc + vy2 * EPb)
              - 2 * (vy2 * vv_z * UXZ + vv_y * vz2 * UXY)
              //-4*vv_x*vy2*vz2*rho
              - 2 * vv_x * vv_y * vv_z * (pi_y * vv_z + pi_z * vv_y) + vy2 * vz2 * (vv_x * rho - pi_x)));
      const MFloat xxyzz =
          ((1. / 8.)
           * (Seb - Neb - Nef - Nwb - Nwf + Sef + Swb + Swf - vv_y * C - 2 * (vv_x * UXYZZ + vv_z * UXXYZ) - vx2 * YpYYY
              - vz2 * YnYYY + 2 * vv_y * (-vv_x * XnXXX - vv_z * ZpZZZ) - 4 * vv_x * vv_z * (XYZ + vv_y * UXZ)
              - vv_y * (vx2 * EPb + vz2 * EPa)
              - 2 * (vx2 * vv_z * UYZ + vv_x * vz2 * UXY)
              //-4*vv_y*vx2*vz2*rho
              - 2 * vv_x * vv_y * vv_z * (pi_x * vv_z + pi_z * vv_x) + vx2 * vz2 * (vv_y * rho - pi_y)));
      const MFloat xxyyz =
          ((1. / 8.)
           * (Neb - Nef + Nwb - Nwf + Seb - Sef + Swb - Swf - vv_z * CA - 2 * (vv_x * UXYYZ + vv_y * UXXYZ)
              - vx2 * ZnZZZ - vy2 * ZpZZZ + 2 * vv_z * (-vv_x * XpXXX - vv_y * YnYYY)
              - 4 * vv_x * vv_y * (XYZ + vv_z * UXY) - vv_z * (vx2 * EPc + vy2 * EPa)
              - 2 * (vx2 * vv_y * UYZ + vv_x * vy2 * UXZ)
              //-4*vx2*vy2*vv_z*rho
              - 2 * vv_x * vv_y * vv_z * (pi_x * vv_y + pi_y * vv_x) + vx2 * vy2 * (vv_z * rho - pi_z)));

      const MFloat xxyyzz =
          ((1. / 8.)
           * (rho / 27. - Neb - Nef - Nwb - Nwf - Seb - Sef - Swb - Swf
              + 2
                    * (+vv_x * (Neb + Nef - Nwb - Nwf + Seb + Sef - Swb - Swf + 8 * xyyzz)
                       + vv_y * (Neb + Nef + Nwb + Nwf - Seb - Sef - Swb - Swf + 8 * xxyzz)
                       + vv_z * (-Neb + Nef - Nwb + Nwf - Seb + Sef - Swb + Swf + 8 * xxyyz))
              + vx2 * A + vy2 * C + vz2 * CA + 4 * (vv_x * vv_y * UXYZZ + vv_x * vv_z * UXYYZ + vv_y * vv_z * UXXYZ)
              + 2
                    * (vx2 * vv_y * YpYYY + vx2 * vv_z * ZnZZZ + vy2 * vv_z * ZpZZZ + vv_x * vy2 * XnXXX
                       + vv_x * vz2 * XpXXX + vz2 * vv_y * YnYYY)
              + 8 * vv_x * vv_y * vv_z * XYZ + vx2 * vz2 * EPc + vx2 * vy2 * EPb + vy2 * vz2 * EPa
              + 4 * vv_x * vv_y * vv_z * (vv_z * UXY + vv_y * UXZ + vv_x * UYZ)
              + 2 * (vx2 * vy2 * vv_z * pi_z + vx2 * vv_y * pi_y * vz2 + vv_x * pi_x * vy2 * vz2)
              - rho * vx2 * vy2 * vz2));

      a_distribution(pCellId, 26) = R + (12 * p - 30 * pe - 8 * xxyyzz);
      a_distribution(pCellId, 7) =
          Nw + c + a - xxx - yyy - x + y + uxy + p + 8 * pe + 2 * (-xxyyzz - xxyzz + xyyzz + uxyzz);
      a_distribution(pCellId, 0) = W - 2 * a + 4 * x - 11 * pe + vxy + vxz - 4 * p + 4 * (-xyyzz + xxyyzz);
      a_distribution(pCellId, 6) =
          Sw + a + c - xxx + yyy - x - y - uxy + p + 8 * pe + 2 * (-xxyyzz + xxyzz + xyyzz - uxyzz);
      a_distribution(pCellId, 2) = S - 2 * c + 4 * y - 11 * pe - vxy - 4 * p + 4 * (-xxyzz + xxyyzz);
      a_distribution(pCellId, 8) =
          Se + a + c + xxx + yyy + x - y + uxy + p + 8 * pe + 2 * (-xxyyzz + xxyzz - xyyzz + uxyzz);
      a_distribution(pCellId, 1) = E - 2 * a - 4 * x - 11 * pe + vxy + vxz - 4 * p + 4 * (xyyzz + xxyyzz);
      a_distribution(pCellId, 9) =
          Ne + c + a + x + y + xxx - yyy - uxy + p + 8 * pe + 2 * (-xxyyzz - xxyzz - xyyzz - uxyzz);
      a_distribution(pCellId, 3) = N - 2 * c - 4 * y - 11 * pe - vxy - 4 * p + 4 * (xxyzz + xxyyzz); // s_n;
      a_distribution(pCellId, 17) =
          Nf - a + yyy - zzz + y + z - uyz + p + 8 * pe + 2 * (-xxyyzz - xxyzz - xxyyz - uxxyz);
      a_distribution(pCellId, 16) =
          Nb - a + yyy + zzz + y - z + uyz + p + 8 * pe + 2 * (-xxyyzz - xxyzz + xxyyz + uxxyz);
      a_distribution(pCellId, 15) =
          Sf - a - yyy - zzz - y + z + uyz + p + 8 * pe + 2 * (-xxyyzz + xxyzz - xxyyz + uxxyz);
      a_distribution(pCellId, 14) =
          Sb - a - yyy + zzz - y - z - uyz + p + 8 * pe + 2 * (-xxyyzz + xxyzz + xxyyz - uxxyz);
      a_distribution(pCellId, 13) =
          Ef - c - xxx + zzz + x + z - uxz + p + 8 * pe + 2 * (-xxyyzz - xyyzz - xxyyz - uxyyz);
      a_distribution(pCellId, 12) =
          Eb - c - xxx - zzz + x - z + uxz + p + 8 * pe + 2 * (-xxyyzz - xyyzz + xxyyz + uxyyz);
      a_distribution(pCellId, 11) =
          Wf - c + xxx + zzz - x + z + uxz + p + 8 * pe + 2 * (-xxyyzz + xyyzz - xxyyz + uxyyz);
      a_distribution(pCellId, 10) =
          Wb - c + xxx - zzz - x - z - uxz + p + 8 * pe + 2 * (-xxyyzz + xyyzz + xxyyz - uxyyz);
      a_distribution(pCellId, 5) = F + 2 * a + 2 * c - 4 * z - 11 * pe - vxz - 4 * p + 4 * (xxyyzz + xxyyz);
      a_distribution(pCellId, 4) = B + 2 * a + 2 * c + 4 * z - 11 * pe - vxz - 4 * p + 4 * (-xxyyz + xxyyzz);
      a_distribution(pCellId, 21) = Nwf - xyz + uxxyz - uxyyz - uxyzz - (-xxyyz + xyyzz - xxyzz) + xxyyzz;
      a_distribution(pCellId, 20) = Nwb + xyz - uxxyz + uxyyz - uxyzz - (xxyyz + xyyzz - xxyzz) + xxyyzz;
      a_distribution(pCellId, 25) = Nef + xyz + uxxyz + uxyyz + uxyzz - (-xxyyz - xyyzz - xxyzz) + xxyyzz;
      a_distribution(pCellId, 24) = Neb - xyz - uxxyz - uxyyz + uxyzz - (xxyyz - xyyzz - xxyzz) + xxyyzz;
      a_distribution(pCellId, 19) = Swf + xyz - uxxyz - uxyyz + uxyzz - (-xxyyz + xyyzz + xxyzz) + xxyyzz;
      a_distribution(pCellId, 18) = Swb - xyz + uxxyz + uxyyz + uxyzz - (xxyyz + xyyzz + xxyzz) + xxyyzz;
      a_distribution(pCellId, 23) = Sef - xyz - uxxyz + uxyyz - uxyzz - (-xxyyz - xyyzz + xxyzz) + xxyyzz;
      a_distribution(pCellId, 22) = Seb + xyz + uxxyz - uxyyz - uxyzz - (xxyyz - xyyzz + xxyzz) + xxyyzz;
    }
  }
  if(m_calculateDissipation && globalTimeStep % m_solutionInterval == 0) {
    calculateDissipation();
  }
}

/** \fn LbSolverDxQy::clb_collision_step()
 *  \author Miro Gondrum (Refactored)
 *  \date   03.11.2021
 *  \propVal{solverMethod,MAIA_LATTICE_CLB}
 *
 *  This function wraps clb_collision_step_base to access CLB-Algorithm without
 *  a Smagorisnky turbulence model.
 *  \collstep{LbSolverDxQy::clb_collision_step(), LBCclb, Cascaded Collision Step}
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbSolverDxQy<nDim, nDist, SysEqn>::clb_collision_step() {
  clb_collision_step_base<false>();
}

/** \fn LbSolverDxQy::clb_smagorinsky_collision_step()
 *  \author Miro Gondrum (Refactored)
 *  \date   03.11.2021
 *  \propVal{solverMethod,MAIA_LATTICE_CLB_SMAGORINSKY}
 *
 *  This function wraps clb_collision_step_base to access CLB-Algorithm with
 *  a Smagorisnky turbulence model.
 *  \collstep{LbSolverDxQy::clb_smagorinsky_collision_step(), LBCclb_smago, Cascaded Collision Step (Smagorinsky)}
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbSolverDxQy<nDim, nDist, SysEqn>::clb_smagorinsky_collision_step() {
  clb_collision_step_base<true>();
}

/** \fn LbSolverDxQy::cumulant_collision_step()
 *  \author Max Mustermann
 *  \date 1.1.1870
 *  \propVal{solverMethod,MAIA_LATTICE_CUMULANT}
 *
 *  Collision step based on countable cumulants.<br>
 *  Ref: Geier et al. 2015, https://doi.org/10.1016/j.camwa.2015.05.001<br>
 *  \collstep{LbSolverDxQy::cumulant_collision_step(), LBCcumulant, Cumulant Collision Step}
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbSolverDxQy<nDim, nDist, SysEqn>::cumulant_collision_step() {
  TERMM(1, "Cumulant collision step only available for D3Q27, yet!");
}

template <>
void LbSolverDxQy<3, 27, maia::lb::LbSysEqnCompressible<3, 27>>::cumulant_collision_step() {
  TRACE();

  constexpr MInt nDist = 27;

  constexpr MFloat omega2 = F1;
  constexpr MFloat omega3 = F1;
  constexpr MFloat omega4 = F1;
  constexpr MFloat omega5 = F1;
  constexpr MFloat omega6 = F1;
  constexpr MFloat omega7 = F1;
  constexpr MFloat omega8 = F1;
  constexpr MFloat omega9 = F1;
  constexpr MFloat omega10 = F1;

  // Required for ported functions, since GPUs cannot access the global variable globalTimeStep
  const MInt gTS = globalTimeStep;
  const MInt maxLevel_ = maxLevel();

  maia::parallelFor<true>(0, m_currentMaxNoCells, [=](MInt index) {
    const MInt pCellId = m_activeCellList[index];
    const MInt lvlDiff = maxLevel_ - a_level(pCellId);
    if((gTS - 1) % IPOW2(lvlDiff) != 0) return;

    const MFloat l_u = a_variable(pCellId, PV->U);
    const MFloat l_v = a_variable(pCellId, PV->V);
    const MFloat l_w = a_variable(pCellId, PV->W);
    const MFloat l_rho = a_variable(pCellId, PV->RHO);

    // Calculation of relaxation rate
    //-----------------------------------------------------------------------------
    const MFloat tau = F1B2 + F1BCSsq * a_nu(pCellId) * FFPOW2(lvlDiff); // tau_0
    const MFloat omega1 = 1.0 / tau;

    // Calculation of the central moments
    //-----------------------------------------------------------------------------
    const MFloat F1Brho = F1 / l_rho;
    const MFloat u[3] = {l_u, l_v, l_w};
    const MFloat uSq[3] = {u[0] * u[0], u[1] * u[1], u[2] * u[2]};

    MFloat centralMoments[3][3][3] = {{{F0, F0, F0}, {F0, F0, F0}, {F0, F0, F0}},
                                      {{F0, F0, F0}, {F0, F0, F0}, {F0, F0, F0}},
                                      {{F0, F0, F0}, {F0, F0, F0}, {F0, F0, F0}}};

    // info: these loops are consuming (~50%)
    for(MInt dist = 0; dist < 27; dist++) {
#ifdef WAR_NVHPC_PSTL
      const MFloat xDir = (m_idFld[dist][0] - F1) - u[0];
      const MFloat yDir = (m_idFld[dist][1] - F1) - u[1];
      const MFloat zDir = (m_idFld[dist][2] - F1) - u[2];
#else
      const MFloat xDir = (Ld::idFld(dist, 0) - F1) - u[0];
      const MFloat yDir = (Ld::idFld(dist, 1) - F1) - u[1];
      const MFloat zDir = (Ld::idFld(dist, 2) - F1) - u[2];
#endif
      const MFloat oldDist = a_oldDistribution(pCellId, dist);
      const MFloat pow012[3][3] = {{F1, xDir, xDir * xDir}, {F1, yDir, yDir * yDir}, {F1, zDir, zDir * zDir}};

      for(MInt i = 0; i < 3; i++) {
        for(MInt j = 0; j < 3; j++) {
          const MFloat tmp = pow012[0][i] * pow012[1][j] * oldDist;
          centralMoments[i][j][0] += pow012[2][0] * tmp;
          centralMoments[i][j][1] += pow012[2][1] * tmp;
          centralMoments[i][j][2] += pow012[2][2] * tmp;
        }
      }
    }

    // Transformation of central moments in cumulants
    //-----------------------------------------------------------------------------
    MFloat cumulants[3][3][3];

    cumulants[0][0][0] = centralMoments[0][0][0];
    cumulants[1][0][0] = centralMoments[1][0][0];
    cumulants[0][1][0] = centralMoments[0][1][0];
    cumulants[0][0][1] = centralMoments[0][0][1];
    cumulants[2][0][0] = centralMoments[2][0][0];
    cumulants[0][2][0] = centralMoments[0][2][0];
    cumulants[0][0][2] = centralMoments[0][0][2];
    cumulants[1][1][0] = centralMoments[1][1][0];
    cumulants[1][0][1] = centralMoments[1][0][1];
    cumulants[0][1][1] = centralMoments[0][1][1];
    cumulants[1][2][0] = centralMoments[1][2][0];
    cumulants[1][0][2] = centralMoments[1][0][2];
    cumulants[0][1][2] = centralMoments[0][1][2];
    cumulants[2][0][1] = centralMoments[2][0][1];
    cumulants[2][1][0] = centralMoments[2][1][0];
    cumulants[0][2][1] = centralMoments[0][2][1];
    cumulants[1][1][1] = centralMoments[1][1][1];

    cumulants[2][1][1] =
        (centralMoments[2][1][1]
         - (centralMoments[2][0][0] * centralMoments[0][1][1] + F2 * centralMoments[1][1][0] * centralMoments[1][0][1])
               * F1Brho);
    cumulants[1][2][1] =
        (centralMoments[1][2][1]
         - (centralMoments[0][2][0] * centralMoments[1][0][1] + F2 * centralMoments[1][1][0] * centralMoments[0][1][1])
               * F1Brho);
    cumulants[1][1][2] =
        (centralMoments[1][1][2]
         - (centralMoments[0][0][2] * centralMoments[1][1][0] + F2 * centralMoments[0][1][1] * centralMoments[1][0][1])
               * F1Brho);

    cumulants[2][2][0] =
        (centralMoments[2][2][0]
         - (centralMoments[2][0][0] * centralMoments[0][2][0] + F2 * centralMoments[1][1][0] * centralMoments[1][1][0])
               * F1Brho);
    cumulants[2][0][2] =
        (centralMoments[2][0][2]
         - (centralMoments[2][0][0] * centralMoments[0][0][2] + F2 * centralMoments[1][0][1] * centralMoments[1][0][1])
               * F1Brho);
    cumulants[0][2][2] =
        (centralMoments[0][2][2]
         - (centralMoments[0][0][2] * centralMoments[0][2][0] + F2 * centralMoments[0][1][1] * centralMoments[0][1][1])
               * F1Brho);

    cumulants[1][2][2] =
        (centralMoments[1][2][2]
         - (centralMoments[0][0][2] * centralMoments[1][2][0] + centralMoments[0][2][0] * centralMoments[1][0][2]
            + F4 * centralMoments[0][1][1] * centralMoments[1][1][1]
            + F2
                  * (centralMoments[1][0][1] * centralMoments[0][2][1]
                     + centralMoments[1][1][0] * centralMoments[0][1][2]))
               * F1Brho);
    cumulants[2][1][2] =
        (centralMoments[2][1][2]
         - (centralMoments[2][0][0] * centralMoments[0][1][2] + centralMoments[0][0][2] * centralMoments[2][1][0]
            + F4 * centralMoments[1][0][1] * centralMoments[1][1][1]
            + F2
                  * (centralMoments[1][1][0] * centralMoments[1][0][2]
                     + centralMoments[0][1][1] * centralMoments[2][0][1]))
               * F1Brho);
    cumulants[2][2][1] =
        (centralMoments[2][2][1]
         - (centralMoments[2][0][0] * centralMoments[0][2][1] + centralMoments[0][2][0] * centralMoments[2][0][1]
            + F4 * centralMoments[1][1][0] * centralMoments[1][1][1]
            + F2
                  * (centralMoments[1][0][1] * centralMoments[1][2][0]
                     + centralMoments[0][1][1] * centralMoments[2][1][0]))
               * F1Brho);

    const MFloat B0 =
        F4 * centralMoments[1][1][1] * centralMoments[1][1][1] + centralMoments[2][0][0] * centralMoments[0][2][2]
        + centralMoments[0][2][0] * centralMoments[2][0][2] + centralMoments[0][0][2] * centralMoments[2][2][0];

    const MFloat B1 = centralMoments[0][1][1] * centralMoments[2][1][1]
                      + centralMoments[1][0][1] * centralMoments[1][2][1]
                      + centralMoments[1][1][0] * centralMoments[1][1][2];

    const MFloat B2 = centralMoments[1][2][0] * centralMoments[1][0][2]
                      + centralMoments[2][1][0] * centralMoments[0][1][2]
                      + centralMoments[2][0][1] * centralMoments[0][2][1];

    const MFloat B3 = centralMoments[1][1][0] * centralMoments[1][0][1] * centralMoments[0][1][1];

    const MFloat B4 = centralMoments[1][0][1] * centralMoments[1][0][1] * centralMoments[0][2][0]
                      + centralMoments[1][1][0] * centralMoments[1][1][0] * centralMoments[0][0][2]
                      + centralMoments[0][1][1] * centralMoments[0][1][1] * centralMoments[2][0][0];

    const MFloat B5 = centralMoments[2][0][0] * centralMoments[0][2][0] * centralMoments[0][0][2];

    cumulants[2][2][2] = (centralMoments[2][2][2] - (F1 * B0 + F4 * B1 + F2 * B2) * F1Brho
                          + (16.0 * B3 + F4 * B4 + F2 * B5) * F1Brho * F1Brho);

    // Cumulant collision
    //-----------------------------------------------------------------------------
    const MFloat sum = cumulants[2][0][0] + cumulants[0][2][0] + cumulants[0][0][2];
    const MFloat diff1 = cumulants[2][0][0] - cumulants[0][2][0];
    const MFloat diff2 = cumulants[2][0][0] - cumulants[0][0][2];
    const MFloat Dxu =
        -omega1 * F1Brho * F1B2 * (diff1 + diff2) - omega2 * F1Brho * F1B2 * (sum - centralMoments[0][0][0]);
    const MFloat Dyv = Dxu + F3B2 * omega1 * F1Brho * diff1;
    const MFloat Dzw = Dxu + F3B2 * omega1 * F1Brho * diff2;
    const MFloat A1 = (F1 - omega1) * diff1 - F3 * l_rho * (F1 - F1B2 * omega1) * (uSq[0] * Dxu - uSq[1] * Dyv);
    const MFloat A2 = (F1 - omega1) * diff2 - F3 * l_rho * (F1 - F1B2 * omega1) * (uSq[0] * Dxu - uSq[2] * Dzw);
    const MFloat A3 = (F1 - omega2) * sum
                      - F3 * l_rho * (F1 - F1B2 * omega2) * (uSq[0] * Dxu + uSq[1] * Dyv + uSq[2] * Dzw)
                      + omega2 * centralMoments[0][0][0];

    const MFloat A4 = (F1 - omega3) * (cumulants[1][2][0] + cumulants[1][0][2]);
    const MFloat A5 = (F1 - omega3) * (cumulants[2][1][0] + cumulants[0][1][2]);
    const MFloat A6 = (F1 - omega3) * (cumulants[2][0][1] + cumulants[0][2][1]);
    const MFloat A7 = (F1 - omega4) * (cumulants[1][2][0] - cumulants[1][0][2]);
    const MFloat A8 = (F1 - omega4) * (cumulants[2][1][0] - cumulants[0][1][2]);
    const MFloat A9 = (F1 - omega4) * (cumulants[2][0][1] - cumulants[0][2][1]);

    const MFloat A10 = (F1 - omega6) * (cumulants[2][2][0] - F2 * cumulants[2][0][2] + cumulants[0][2][2]);
    const MFloat A11 = (F1 - omega6) * (cumulants[2][2][0] + cumulants[2][0][2] - F2 * cumulants[0][2][2]);
    const MFloat A12 = (F1 - omega7) * (cumulants[2][2][0] + cumulants[2][0][2] + cumulants[0][2][2]);

    cumulants[1][1][0] *= (F1 - omega1);
    cumulants[0][1][1] *= (F1 - omega1);
    cumulants[1][0][1] *= (F1 - omega1);
    cumulants[2][0][0] = (A1 + A2 + A3) * F1B3;
    cumulants[0][0][2] = (A1 - F2 * A2 + A3) * F1B3;
    cumulants[0][2][0] = (-F2 * A1 + A2 + A3) * F1B3;
    cumulants[1][2][0] = (A4 + A5) * F1B2;
    cumulants[1][0][2] = (A4 - A5) * F1B2;
    cumulants[2][1][0] = (A6 + A7) * F1B2;
    cumulants[0][1][2] = (A6 - A7) * F1B2;
    cumulants[2][0][1] = (A8 + A9) * F1B2;
    cumulants[0][2][1] = (A8 - A9) * F1B2;
    cumulants[1][1][1] *= (F1 - omega5);
    cumulants[2][2][0] = F1B3 * (A10 + A11 + A12);
    cumulants[2][0][2] = F1B3 * (A12 - A10);
    cumulants[0][2][2] = F1B3 * (A12 - A11);
    cumulants[2][1][1] *= (F1 - omega8);
    cumulants[1][2][1] *= (F1 - omega8);
    cumulants[1][1][2] *= (F1 - omega8);
    cumulants[2][2][1] *= (F1 - omega9);
    cumulants[2][1][2] *= (F1 - omega9);
    cumulants[1][2][2] *= (F1 - omega9);
    cumulants[2][2][2] *= (F1 - omega10);

    // Transformation of the cumulants in central moments
    //-----------------------------------------------------------------------------
    centralMoments[0][0][0] = cumulants[0][0][0];
    centralMoments[1][0][0] = -cumulants[1][0][0];
    centralMoments[0][1][0] = -cumulants[0][1][0];
    centralMoments[0][0][1] = -cumulants[0][0][1];
    centralMoments[1][1][0] = cumulants[1][1][0];
    centralMoments[1][0][1] = cumulants[1][0][1];
    centralMoments[0][1][1] = cumulants[0][1][1];
    centralMoments[2][0][0] = cumulants[2][0][0];
    centralMoments[0][2][0] = cumulants[0][2][0];
    centralMoments[0][0][2] = cumulants[0][0][2];
    centralMoments[2][1][0] = cumulants[2][1][0];
    centralMoments[2][0][1] = cumulants[2][0][1];
    centralMoments[0][2][1] = cumulants[0][2][1];
    centralMoments[1][2][0] = cumulants[1][2][0];
    centralMoments[1][0][2] = cumulants[1][0][2];
    centralMoments[0][1][2] = cumulants[0][1][2];
    centralMoments[1][1][1] = cumulants[1][1][1];

    centralMoments[2][1][1] =
        (cumulants[2][1][1]
         + (centralMoments[2][0][0] * centralMoments[0][1][1] + F2 * centralMoments[1][1][0] * centralMoments[1][0][1])
               * F1Brho);
    centralMoments[1][2][1] =
        (cumulants[1][2][1]
         + (centralMoments[0][2][0] * centralMoments[1][0][1] + F2 * centralMoments[1][1][0] * centralMoments[0][1][1])
               * F1Brho);
    centralMoments[1][1][2] =
        (cumulants[1][1][2]
         + (centralMoments[0][0][2] * centralMoments[1][1][0] + F2 * centralMoments[0][1][1] * centralMoments[1][0][1])
               * F1Brho);

    centralMoments[2][2][0] =
        (cumulants[2][2][0]
         + (centralMoments[2][0][0] * centralMoments[0][2][0] + F2 * centralMoments[1][1][0] * centralMoments[1][1][0])
               * F1Brho);
    centralMoments[2][0][2] =
        (cumulants[2][0][2]
         + (centralMoments[2][0][0] * centralMoments[0][0][2] + F2 * centralMoments[1][0][1] * centralMoments[1][0][1])
               * F1Brho);
    centralMoments[0][2][2] =
        (cumulants[0][2][2]
         + (centralMoments[0][0][2] * centralMoments[0][2][0] + F2 * centralMoments[0][1][1] * centralMoments[0][1][1])
               * F1Brho);

    centralMoments[1][2][2] =
        (cumulants[1][2][2]
         + (centralMoments[0][0][2] * centralMoments[1][2][0] + centralMoments[0][2][0] * centralMoments[1][0][2]
            + F4 * centralMoments[0][1][1] * centralMoments[1][1][1]
            + F2
                  * (centralMoments[1][0][1] * centralMoments[0][2][1]
                     + centralMoments[1][1][0] * centralMoments[0][1][2]))
               * F1Brho);
    centralMoments[2][1][2] =
        (cumulants[2][1][2]
         + (centralMoments[2][0][0] * centralMoments[0][1][2] + centralMoments[0][0][2] * centralMoments[2][1][0]
            + F4 * centralMoments[1][0][1] * centralMoments[1][1][1]
            + F2
                  * (centralMoments[1][1][0] * centralMoments[1][0][2]
                     + centralMoments[0][1][1] * centralMoments[2][0][1]))
               * F1Brho);
    centralMoments[2][2][1] =
        (cumulants[2][2][1]
         + (centralMoments[2][0][0] * centralMoments[0][2][1] + centralMoments[0][2][0] * centralMoments[2][0][1]
            + F4 * centralMoments[1][1][0] * centralMoments[1][1][1]
            + F2
                  * (centralMoments[1][0][1] * centralMoments[1][2][0]
                     + centralMoments[0][1][1] * centralMoments[2][1][0]))
               * F1Brho);

    const MFloat D0 =
        F4 * centralMoments[1][1][1] * centralMoments[1][1][1] + centralMoments[2][0][0] * centralMoments[0][2][2]
        + centralMoments[0][2][0] * centralMoments[2][0][2] + centralMoments[0][0][2] * centralMoments[2][2][0];

    const MFloat D1 = centralMoments[0][1][1] * centralMoments[2][1][1]
                      + centralMoments[1][0][1] * centralMoments[1][2][1]
                      + centralMoments[1][1][0] * centralMoments[1][1][2];

    const MFloat D2 = centralMoments[1][2][0] * centralMoments[1][0][2]
                      + centralMoments[2][1][0] * centralMoments[0][1][2]
                      + centralMoments[2][0][1] * centralMoments[0][2][1];

    const MFloat D3 = centralMoments[1][1][0] * centralMoments[1][0][1] * centralMoments[0][1][1];

    const MFloat D4 = centralMoments[1][0][1] * centralMoments[1][0][1] * centralMoments[0][2][0]
                      + centralMoments[1][1][0] * centralMoments[1][1][0] * centralMoments[0][0][2]
                      + centralMoments[0][1][1] * centralMoments[0][1][1] * centralMoments[2][0][0];

    const MFloat D5 = centralMoments[2][0][0] * centralMoments[0][2][0] * centralMoments[0][0][2];

    centralMoments[2][2][2] = (cumulants[2][2][2] + (F1 * D0 + F4 * D1 + F2 * D2) * F1Brho
                               - (16.0 * D3 + F4 * D4 + F2 * D5) * F1Brho * F1Brho);

    const MFloat coeff0[3][3] = {{(uSq[0] - u[0]) * F1B2, u[0] - F1B2, F1B2},
                                 {F1 - uSq[0], -F2 * u[0], -F1},
                                 {(uSq[0] + u[0]) * F1B2, u[0] + F1B2, F1B2}};
    const MFloat coeff1[3][3] = {{(uSq[1] - u[1]) * F1B2, u[1] - F1B2, F1B2},
                                 {F1 - uSq[1], -F2 * u[1], -F1},
                                 {(uSq[1] + u[1]) * F1B2, u[1] + F1B2, F1B2}};
    const MFloat coeff2[3][3] = {{(uSq[2] - u[2]) * F1B2, u[2] - F1B2, F1B2},
                                 {F1 - uSq[2], -F2 * u[2], -F1},
                                 {(uSq[2] + u[2]) * F1B2, u[2] + F1B2, F1B2}};

    // info: this loop is very consuming (~30%)
    for(MInt dist = 0; dist < nDist; dist++) {
#ifdef WAR_NVHPC_PSTL
      const MInt i = m_idFld[dist][0];
      const MInt j = m_idFld[dist][1];
      const MInt k = m_idFld[dist][2];
#else
      const MInt i = Ld::idFld(dist, 0);
      const MInt j = Ld::idFld(dist, 1);
      const MInt k = Ld::idFld(dist, 2);
#endif
      const MFloat ki_00 = coeff0[i][0] * centralMoments[0][0][0] + coeff0[i][1] * centralMoments[1][0][0]
                           + coeff0[i][2] * centralMoments[2][0][0];
      const MFloat ki_01 = coeff0[i][0] * centralMoments[0][0][1] + coeff0[i][1] * centralMoments[1][0][1]
                           + coeff0[i][2] * centralMoments[2][0][1];
      const MFloat ki_02 = coeff0[i][0] * centralMoments[0][0][2] + coeff0[i][1] * centralMoments[1][0][2]
                           + coeff0[i][2] * centralMoments[2][0][2];
      const MFloat ki_10 = coeff0[i][0] * centralMoments[0][1][0] + coeff0[i][1] * centralMoments[1][1][0]
                           + coeff0[i][2] * centralMoments[2][1][0];
      const MFloat ki_11 = coeff0[i][0] * centralMoments[0][1][1] + coeff0[i][1] * centralMoments[1][1][1]
                           + coeff0[i][2] * centralMoments[2][1][1];
      const MFloat ki_12 = coeff0[i][0] * centralMoments[0][1][2] + coeff0[i][1] * centralMoments[1][1][2]
                           + coeff0[i][2] * centralMoments[2][1][2];
      const MFloat ki_20 = coeff0[i][0] * centralMoments[0][2][0] + coeff0[i][1] * centralMoments[1][2][0]
                           + coeff0[i][2] * centralMoments[2][2][0];
      const MFloat ki_21 = coeff0[i][0] * centralMoments[0][2][1] + coeff0[i][1] * centralMoments[1][2][1]
                           + coeff0[i][2] * centralMoments[2][2][1];
      const MFloat ki_22 = coeff0[i][0] * centralMoments[0][2][2] + coeff0[i][1] * centralMoments[1][2][2]
                           + coeff0[i][2] * centralMoments[2][2][2];

      const MFloat kij_0 = coeff1[j][0] * ki_00 + coeff1[j][1] * ki_10 + coeff1[j][2] * ki_20;
      const MFloat kij_1 = coeff1[j][0] * ki_01 + coeff1[j][1] * ki_11 + coeff1[j][2] * ki_21;
      const MFloat kij_2 = coeff1[j][0] * ki_02 + coeff1[j][1] * ki_12 + coeff1[j][2] * ki_22;

      a_distribution(pCellId, dist) = coeff2[k][0] * kij_0 + coeff2[k][1] * kij_1 + coeff2[k][2] * kij_2;
    }
  });
}

/** \fn LbSolverDxQy::bgkc_collision_step()
 *  \author Ludwig Eduard Boltzmann
 *  \date   20.2.1844
 *  \propVal{solverMethod,MAIA_LATTICE_BGKC}
 *
 *  Common BGK collision step.<br>
 *  \collstep{LbSolverDxQy::bgkc_collision_step(), LBCbgkc, BGK Collision Step}
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbSolverDxQy<nDim, nDist, SysEqn>::bgkc_collision_step() {
  TRACE();
  const MInt globalTimeStep_ = globalTimeStep;
  // loop over all cells
  maia::parallelFor<true>(0, m_currentMaxNoCells, [=](MInt i) {
    const MInt pCellId = m_activeCellList[i];
    const MInt lvlDiff = maxLevel() - this->a_level(pCellId);
    if((globalTimeStep_ - 1) % IPOW2(lvlDiff) != 0) return;
    // Load macroscopic variables
    const MFloat rho = a_variable(pCellId, PV->RHO);
    std::array<MFloat, nDim> uu;
    for(MInt j = 0; j < nDim; j++) {
      uu[j] = a_variable(pCellId, PV->U + j);
    }
    // Collision
    const MFloat omega = 1.0 / (0.5 + F1BCSsq * a_nu(pCellId) * FFPOW2(lvlDiff));
    std::array<MFloat, nDist> eqDist;
    eqDist = getEqDists(rho, uu.data());
    for(MInt j = 0; j < nDist; j++) {
      a_distribution(pCellId, j) = a_oldDistribution(pCellId, j) + omega * (eqDist[j] - a_oldDistribution(pCellId, j));
    }
  });
}

/** \fn LbSolverDxQy::bgki_collision_step()
 *  \author Ludwig Eduard Boltzmann
 *  \date   20.2.1844
 *  \propVal{solverMethod,MAIA_LATTICE_BGK}
 *
 *  Collision step of the incompressible LBGK algorithm<br>
 *  \collstep{LbSolverDxQy::bgki_collision_step(), LBCbgki, Incompressible BGK Collision Step}
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbSolverDxQy<nDim, nDist, SysEqn>::bgki_collision_step() {
  TRACE();

  const MInt pvrho = PV->RHO;

  constexpr MInt fldlen = Ld::dxQyFld();

  // Update the global Reynolds number if we use the tanh-adaption
  if(m_tanhInit && (globalTimeStep > m_initStartTime)) {
    const MFloat scale =
        0.5 * (m_tanhScaleFactor * tanh((5.0 * (MFloat)(globalTimeStep - m_initStartTime) / (MFloat)m_initTime) - 2.5))
        + 0.5;

    m_Re = m_initRe + scale * (m_finalRe - m_initRe);

    // final timestep reached
    if(globalTimeStep == m_initStartTime + m_initTime) {
      m_Re = m_finalRe;
      m_tanhInit = 0;
    }
  }

  // Update the according nu
  m_nu = m_Ma * LBCS / m_Re * m_referenceLength;

#ifndef WAR_NVHPC_PSTL
  // inner loop references
  const MInt* const RESTRICT activeCellList = ALIGNED_I(m_activeCellList);
  MFloat* const RESTRICT cnu = ALIGNED_MF(&(a_nu(0)));
  const MInt* const RESTRICT clevel = ALIGNED_I(&(this->m_cells.level(0)));
  MFloat* const RESTRICT oldDistributions = ALIGNED_MF(&(a_oldDistribution(0, 0)));
  MFloat* const RESTRICT distributions = ALIGNED_MF(&(a_distribution(0, 0)));
  MFloat* const RESTRICT variables = ALIGNED_MF(&(a_variable(0, 0)));
  // members and globals
  const MInt maxLevel_ = maxLevel();
  // distances
  const MInt distNu = &a_nu(1) - &a_nu(0);
  const MInt distLevel = &(this->m_cells.level(1)) - &(this->m_cells.level(0));
  const MInt distVars = &(a_variable(1, 0)) - &(a_variable(0, 0));
#endif
  // Required for ported functions, since GPUs cannot access the global variable globalTimeStep
  const MInt gTS = globalTimeStep;

  // loop over all cells
  maia::parallelFor<true>(0, m_currentMaxNoCells, [=](MInt i) {
#ifdef WAR_NVHPC_PSTL
    const MInt pCellId = m_activeCellList[i];
    const MInt index = maxLevel() - this->a_level(pCellId);
#else
    const MInt pCellId = activeCellList[i];
    const MInt index = maxLevel_ - clevel[pCellId * distLevel];
#endif

    // perform collision on finest level in every timestep,
    // on the next coarser level every second timestep, and so on
    if((gTS - 1) % IPOW2(index) == 0) {
      // save nu for bndcnd and interface interpolation
#ifndef WAR_NVHPC_PSTL
      cnu[pCellId * distNu] = m_nu;

      const MInt distStart = pCellId * nDist;
      // varStart = (pCellId * distVars)  + (distVarOldVar *
      // !((MInt)floor((MFloat)localGlobalTimeStep/fpow2[index]) % 2));
      const MInt varStart = pCellId * distVars;

      MFloat* const oldDistributionsStart = &oldDistributions[distStart];
      MFloat* const distributionsStart = &distributions[distStart];
      MFloat* const variablesStart = &variables[varStart];
#else
      a_nu(pCellId) = m_nu;
#endif

      const MFloat tmp_FPOW2 = FPOW2(index);
      const MFloat tmp_FFPOW2 = FFPOW2(index);
      const MFloat l_omega = 2.0 / (1.0 + 6.0 * m_nu * tmp_FFPOW2);

      swap_variables(pCellId);

      // Reinit macroscopic variables before recalculation
      MFloat l_v[nDim] = {0.0};
      MFloat l_rho = 0.0;
      std::array<MFloat, nDist - 1> externalForcing{};

      if(m_particleMomentumCoupling) {
        for(MInt dir = 0; dir < nDim; dir++) {
          for(MInt mi = 0; mi < fldlen; mi++) {
#ifdef WAR_NVHPC_PSTL
            externalForcing[m_nFld[dir * fldlen + mi]] +=
                m_tp[m_distType[m_nFld[dir * fldlen + mi]]] * -1.0 * a_externalForces(pCellId, dir) * F1BCSsq;
            externalForcing[m_pFld[dir * fldlen + mi]] +=
                m_tp[m_distType[m_pFld[dir * fldlen + mi]]] * 1.0 * a_externalForces(pCellId, dir) * F1BCSsq;
#else
            externalForcing.at(Ld::nFld(dir, mi)) +=
                Ld::tp(Ld::distType(Ld::nFld(dir, mi))) * -1.0 * a_externalForces(pCellId, dir) * F1BCSsq;
            externalForcing.at(Ld::pFld(dir, mi)) +=
                Ld::tp(Ld::distType(Ld::pFld(dir, mi))) * 1.0 * a_externalForces(pCellId, dir) * F1BCSsq;
#endif
          }
        }
      }
#ifdef WAR_NVHPC_PSTL
      // Add forcing term and calculate density rho
      for(MInt j = 0; j < nDist - 1; j++) {
        if(m_particleMomentumCoupling) {
          a_oldDistribution(pCellId, j) += tmp_FPOW2 * externalForcing[j];
        }
        // Add forcing term and calculate density rho
        l_rho += a_oldDistribution(pCellId, j);
      }
      // Do not forget the last distribution
      l_rho += a_oldDistribution(pCellId, nDist - 1);
#else
      // Add forcing term and calculate density rho
      for(MInt j = 0; j < nDist - 1; j++) {
        if(m_particleMomentumCoupling) {
          oldDistributionsStart[j] += tmp_FPOW2 * externalForcing.at(j);
        }
        // Add forcing term and calculate density rho
        l_rho += oldDistributionsStart[j];
      }
      // Do not forget the last distribution
      l_rho += oldDistributionsStart[nDist - 1];
#endif


      // Calculate macroscopic variables
#ifndef WAR_NVHPC_PSTL
      if constexpr(nDim == 2) {
        // WARNING!!! The velocity should be calculated as in 3d
        // However, result is correct
        // calculation of u
        l_v[0] = oldDistributionsStart[1] + oldDistributionsStart[4] + oldDistributionsStart[5]
                 - oldDistributionsStart[7] - oldDistributionsStart[6] - oldDistributionsStart[0];
        // calculation of v
        l_v[1] = oldDistributionsStart[7] + oldDistributionsStart[3] + oldDistributionsStart[4]
                 - oldDistributionsStart[6] - oldDistributionsStart[2] - oldDistributionsStart[5];
      } else {
#endif
        for(MInt j = 0; j < fldlen; j++) {
          for(MInt d = 0; d < nDim; d++) {
#ifdef WAR_NVHPC_PSTL
            l_v[d] += a_oldDistribution(pCellId, m_pFld[d * fldlen + j]);
            l_v[d] -= a_oldDistribution(pCellId, m_nFld[d * fldlen + j]);
#else
            l_v[d] += oldDistributionsStart[Ld::pFld(d, j)];
            l_v[d] -= oldDistributionsStart[Ld::nFld(d, j)];
#endif
          }
        }
#ifndef WAR_NVHPC_PSTL
      }
#endif

      // Save new macroscopic variables in cell
#ifdef WAR_NVHPC_PSTL
      for(MInt d = 0; d < nDim; d++) {
        a_variable(pCellId, d) = l_v[d];
      }
      a_variable(pCellId, nDim) = l_rho;
#else
      for(MInt d = 0; d < nDim; d++) {
        variablesStart[d] = l_v[d];
      }
      variablesStart[pvrho] = l_rho;
#endif

      const MFloat sqVel = std::inner_product(&l_v[0], &l_v[nDim], &l_v[0], F0);

      std::array<MFloat, nDist> eqDist;
      eqDist = getEqDists(l_rho, sqVel, l_v);

      // Calculation of new distributions
#ifdef WAR_NVHPC_PSTL
      for(MInt j = 0; j < nDist; j++) {
        a_distribution(pCellId, j) =
            a_oldDistribution(pCellId, j) + l_omega * (eqDist[j] - a_oldDistribution(pCellId, j));
      }
#else
      for(MInt j = 0; j < nDist; j++) {
        distributionsStart[j] = oldDistributionsStart[j] + l_omega * (eqDist[j] - oldDistributionsStart[j]);
      }
#endif
    }
  });
}

/** \fn LbSolverDxQy::bgki_collision_step_Guo_forcing()
 *  \author Johannes Grafen
 *  \propVal{solverMethod,MAIA_LATTICE_BGK_GUO_FORCING}
 *
 *  Collision step of the incompressible LBGK algorithm with forcing method of Guo et al.
 *
 *  adapted from conventional bgki-collision step
 *  Discrete lattice effects on the forcing term in the lattice Boltzmann method
 *  DOI: 0.1103/PhysRevE.65.046308<br>
 *  \collstep{LbSolverDxQy::bgki_collision_step_Guo_forcing(), LBCguo, BGKI Guo forcing collision step}
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbSolverDxQy<nDim, nDist, SysEqn>::bgki_collision_step_Guo_forcing() {
  TRACE();

  if(!m_particleMomentumCoupling)
    std::cerr << "Momentum Coupling is not activated, use MAIA_LATTIC_BGK instead!" << std::endl;

  const MInt pvrho = PV->RHO;

  constexpr MInt fldlen = Ld::dxQyFld();

  // Update the global Reynolds number if we use the tanh-adaption
  if(m_tanhInit && (globalTimeStep > m_initStartTime)) {
    const MFloat scale =
        0.5 * (m_tanhScaleFactor * tanh((5.0 * (MFloat)(globalTimeStep - m_initStartTime) / (MFloat)m_initTime) - 2.5))
        + 0.5;

    m_Re = m_initRe + scale * (m_finalRe - m_initRe);

    // final timestep reached
    if(globalTimeStep == m_initStartTime + m_initTime) {
      m_Re = m_finalRe;
      m_tanhInit = 0;
    }
  }

  // Update the according nu
  m_nu = m_Ma * LBCS / m_Re * m_referenceLength;

  // inner loop references
  const MInt* const RESTRICT activeCellList = ALIGNED_I(m_activeCellList);
  MFloat* const RESTRICT cnu = ALIGNED_MF(&(a_nu(0)));
  const MInt* const RESTRICT clevel = ALIGNED_I(&(this->m_cells.level(0)));
  MFloat* const RESTRICT oldDistributions = ALIGNED_MF(&(a_oldDistribution(0, 0)));
  MFloat* const RESTRICT distributions = ALIGNED_MF(&(a_distribution(0, 0)));
  MFloat* const RESTRICT variables = ALIGNED_MF(&(a_variable(0, 0)));

  // members and globals
  const MFloat nu = m_nu;

  // Required for ported functions, since GPUs cannot access the global variable globalTimeStep
  const MInt gTS = globalTimeStep;
  const MInt maxLevel_ = maxLevel();

  // distances
  const MInt distNu = &a_nu(1) - &a_nu(0);
  const MInt distLevel = &(this->m_cells.level(1)) - &(this->m_cells.level(0));
  const MInt distVars = &(a_variable(1, 0)) - &(a_variable(0, 0));

  // loop over all cells
  maia::parallelFor<false>(0, m_currentMaxNoCells, [=](MInt i) {
    const MInt pCellId = activeCellList[i];
    const MInt index = maxLevel_ - clevel[pCellId * distLevel];

    // perform collision on finest level in every timestep,
    // on the next coarser level every second timestep, and so on
    if((gTS - 1) % IPOW2(index) == 0) {
      // save nu for bndcnd and interface interpolation
      cnu[pCellId * distNu] = nu;

      const MInt distStart = pCellId * nDist;
      // varStart = (pCellId * distVars)  + (distVarOldVar *
      // !((MInt)floor((MFloat)localGlobalTimeStep/fpow2[index]) % 2));
      const MInt varStart = pCellId * distVars;

      MFloat* const oldDistributionsStart = &oldDistributions[distStart];
      MFloat* const distributionsStart = &distributions[distStart];

      MFloat* const variablesStart = &variables[varStart];

      const MFloat tmp_FPOW2 = FPOW2(index);
      const MFloat tmp_FFPOW2 = FFPOW2(index);
      const MFloat l_omega = 2.0 / (1.0 + 6.0 * nu * tmp_FFPOW2);

      swap_variables(pCellId);

      // Reinit macroscopic variables before recalculation
      MFloat l_v[nDim] = {0.0};
      MFloat l_rho = 0.0;
      std::array<MFloat, nDist> externalForcing{0.0};

      // calculate density rho and velocity
      for(MInt j = 0; j < nDist - 1; j++) {
        // Add forcing term and calculate density rho
        l_rho += oldDistributionsStart[j];
      }
      // Do not forget the last distribution
      l_rho += oldDistributionsStart[nDist - 1];

      for(MInt j = 0; j < fldlen; j++) {
        for(MInt d = 0; d < nDim; d++) {
          l_v[d] += oldDistributionsStart[Ld::pFld(d, j)];
          l_v[d] -= oldDistributionsStart[Ld::nFld(d, j)];
        }
      }

      // add forcing term to velocity
      for(MInt d = 0; d < nDim; d++) {
        l_v[d] += tmp_FPOW2 * F1B2 * a_externalForces(pCellId, d);
      }

      // Save new macroscopic variables in cell
      for(MInt d = 0; d < nDim; d++) {
        variablesStart[d] = l_v[d];
      }
      variablesStart[pvrho] = l_rho;

      const MFloat sqVel = std::inner_product(&l_v[0], &l_v[nDim], &l_v[0], F0);

      // Compute Sourceterm:
      for(MInt j = 0; j < nDist - 1; j++) {
        // precompute
        std::array<MFloat, nDim> velDiff{};
        MFloat scalarProductVel{};
        for(MInt dir = 0; dir < nDim; dir++) {
          velDiff.at(dir) = Ld::ppdfDir(j, dir) - l_v[dir];
          scalarProductVel += Ld::ppdfDir(j, dir) * l_v[dir];
        }
        for(MInt dir = 0; dir < nDim; dir++) {
          externalForcing.at(j) +=
              tmp_FPOW2 * (F1 - F1B2 * l_omega) * Ld::tp(Ld::distType(j))
              * (F1BCSsq * velDiff.at(dir) + POW2(F1BCSsq) * scalarProductVel * Ld::ppdfDir(j, dir))
              * a_externalForces(pCellId, dir);
        }
      }

      std::array<MFloat, nDist> eqDist;
      eqDist = getEqDists(l_rho, sqVel, l_v);

      // Calculation of new distributions
      for(MInt j = 0; j < nDist; j++) {
        distributionsStart[j] =
            oldDistributionsStart[j] + l_omega * (eqDist[j] - oldDistributionsStart[j]) + externalForcing.at(j);
      }
    }
  });
}

/** \brief Collision step for Transport Lattice-Boltzmann
 * \author Shota Ito
 * \date 07.06.22
 *
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbSolverDxQy<nDim, nDist, SysEqn>::bgkc_transport_collision_step() {
  TRACE();

  const MInt pvrho = PV->RHO;
  const MInt pvc = PV->C;

  constexpr MInt fldlen = Ld::dxQyFld();

  // Update the global Reynolds number if we use the tanh-adaption
  if(m_tanhInit && (globalTimeStep > m_initStartTime)) {
    const MFloat scale =
        0.5 * (m_tanhScaleFactor * tanh((5.0 * (MFloat)(globalTimeStep - m_initStartTime) / (MFloat)m_initTime) - 2.5))
        + 0.5;

    m_Re = m_initRe + scale * (m_finalRe - m_initRe);

    // final timestep reached
    if(globalTimeStep == m_initStartTime + m_initTime) {
      m_Re = m_finalRe;
      m_tanhInit = 0;
    }
  }

  // Update the according nu
  m_nu = m_Ma * LBCS / m_Re * m_referenceLength;

  // Calculate heat diffusion coefficient
  m_diffusivity = m_nu * (m_Re / m_Pe);

#ifndef WAR_NVHPC_PSTL
  // inner loop references
  const MInt* const RESTRICT activeCellList = ALIGNED_I(m_activeCellList);
  MFloat* const RESTRICT cnu = ALIGNED_MF(&(a_nu(0)));
  MFloat* const RESTRICT cdiffusivity = ALIGNED_MF(&(a_diffusivity(0)));
  const MInt* const RESTRICT clevel = ALIGNED_I(&(this->m_cells.level(0)));
  MFloat* const RESTRICT oldDistributions = ALIGNED_MF(&(a_oldDistribution(0, 0)));
  MFloat* const RESTRICT distributions = ALIGNED_MF(&(a_distribution(0, 0)));
  MFloat* const RESTRICT oldDistributionsTransport = ALIGNED_MF(&(a_oldDistributionTransport(0, 0)));
  MFloat* const RESTRICT distributionsTransport = ALIGNED_MF(&(a_distributionTransport(0, 0)));
  MFloat* const RESTRICT variables = ALIGNED_MF(&(a_variable(0, 0)));

  // members and globals
  const MInt maxLevel_ = maxLevel();
  // distances
  const MInt distNu = &a_nu(1) - &a_nu(0);
  const MInt distDiffusivity = &a_diffusivity(1) - &a_diffusivity(0);
  const MInt distLevel = &(this->m_cells.level(1)) - &(this->m_cells.level(0));
  const MInt distVars = &(a_variable(1, 0)) - &(a_variable(0, 0));
#endif
  // Required for ported functions, since GPUs cannot access the global variable globalTimeStep
  const MInt gTS = globalTimeStep;

  // loop over all cells
  maia::parallelFor<true>(0, m_currentMaxNoCells, [=](MInt i) {
#ifdef WAR_NVHPC_PSTL
    const MInt pCellId = m_activeCellList[i];
    const MInt index = maxLevel() - this->a_level(pCellId);
#else
    const MInt pCellId = activeCellList[i];
    const MInt index = maxLevel_ - clevel[pCellId * distLevel];
#endif

    // perform collision on finest level in every timestep,
    // on the next coarser level every second timestep, and so on
    if((gTS - 1) % IPOW2(index) == 0) {
      // save nu for bndcnd and interface interpolation
#ifndef WAR_NVHPC_PSTL
      cnu[pCellId * distNu] = m_nu;
      cdiffusivity[pCellId * distDiffusivity] = m_diffusivity;

      const MInt distStart = pCellId * nDist;
      // varStart = (pCellId * distVars)  + (distVarOldVar *
      // !((MInt)floor((MFloat)localGlobalTimeStep/fpow2[index]) % 2));
      const MInt varStart = pCellId * distVars;

      MFloat* const oldDistributionsStart = &oldDistributions[distStart];
      MFloat* const distributionsStart = &distributions[distStart];
      MFloat* const oldDistributionsTransportStart = &oldDistributionsTransport[distStart];
      MFloat* const distributionsTransportStart = &distributionsTransport[distStart];

      MFloat* const variablesStart = &variables[varStart];
#else
      a_nu(pCellId) = m_nu;
      a_kappa(pCellId) = m_kappa;
#endif

      const MFloat tmp_FFPOW2 = FFPOW2(index);
      const MFloat l_omega = 2.0 / (1.0 + 6.0 * m_nu * tmp_FFPOW2);
      const MFloat l_omegaD = 2.0 / (1.0 + 6.0 * m_diffusivity * tmp_FFPOW2);

      swap_variables(pCellId);

      // Reinit macroscopic variables before recalculation
      MFloat l_v[nDim] = {0.0};
      MFloat l_rho = 0.0;
      MFloat l_c = 0.0;

#ifdef WAR_NVHPC_PSTL
      // Add forcing term and calculate density rho
      for(MInt j = 0; j < nDist - 1; j++) {
        // Add forcing term and calculate density rho
        l_rho += a_oldDistribution(pCellId, j);
        l_c += a_oldDistributionTransport(pCellId, j);
      }
      // Do not forget the last distribution
      l_rho += a_oldDistribution(pCellId, nDist - 1);
      l_c += a_oldDistributionTransport(pCellId, nDist - 1);
#else
      // Add forcing term and calculate density rho
      for(MInt j = 0; j < nDist - 1; j++) {
        l_rho += oldDistributionsStart[j];
        l_c += oldDistributionsTransportStart[j];
      }
      // Do not forget the last distribution
      l_rho += oldDistributionsStart[nDist - 1];
      l_c += oldDistributionsTransportStart[nDist - 1];
#endif

      // Calculate macroscopic variables
#ifndef WAR_NVHPC_PSTL
      if constexpr(nDim == 2) {
        // WARNING!!! The velocity should be calculated as in 3d
        // However, result is correct
        // calculation of u
        l_v[0] = oldDistributionsStart[1] + oldDistributionsStart[4] + oldDistributionsStart[5]
                 - oldDistributionsStart[7] - oldDistributionsStart[6] - oldDistributionsStart[0];
        // calculation of v
        l_v[1] = oldDistributionsStart[7] + oldDistributionsStart[3] + oldDistributionsStart[4]
                 - oldDistributionsStart[6] - oldDistributionsStart[2] - oldDistributionsStart[5];
      } else {
#endif
        for(MInt j = 0; j < fldlen; j++) {
          for(MInt d = 0; d < nDim; d++) {
#ifdef WAR_NVHPC_PSTL
            l_v[d] += a_oldDistribution(pCellId, m_pFld[d * fldlen + j]);
            l_v[d] -= a_oldDistribution(pCellId, m_nFld[d * fldlen + j]);
#else
            l_v[d] += oldDistributionsStart[Ld::pFld(d, j)];
            l_v[d] -= oldDistributionsStart[Ld::nFld(d, j)];
#endif
          }
        }
#ifndef WAR_NVHPC_PSTL
      }
#endif

      // Save new macroscopic variables in cell
#ifdef WAR_NVHPC_PSTL
      for(MInt d = 0; d < nDim; d++) {
        a_variable(pCellId, d) = l_v[d];
      }
      a_variable(pCellId, pvrho) = l_rho;
      a_variable(pCellId, pvc) = l_c;
#else
      for(MInt d = 0; d < nDim; d++) {
        variablesStart[d] = l_v[d];
      }
      variablesStart[pvrho] = l_rho;
      variablesStart[pvc] = l_c;
#endif

      MFloat u[nDim] = {0.0};
      for(MInt d = 0; d < nDim; d++) {
        u[d] = l_v[d] / l_rho;
      }

      std::array<MFloat, nDist> eqDist;
      eqDist = getEqDists(l_rho, u);
      std::array<MFloat, nDist> eqDistC;
      eqDistC = getEqDistsTransport(l_c, u);

      // Calculation of new distributions
      for(MInt j = 0; j < nDist; j++) {
#ifdef WAR_NVHPC_PSTL
        a_distribution(pCellId, j) =
            a_oldDistribution(pCellId, j) + l_omega * (eqDist[j] - a_oldDistribution(pCellId, j));
        a_distributionTransport(pCellId, j) =
            a_oldDistributionTransport(pCellId, j) + l_omegaD * (eqDistC[j] - a_oldDistributionTransport(pCellId, j));
#else
        distributionsStart[j] = oldDistributionsStart[j] + l_omega * (eqDist[j] - oldDistributionsStart[j]);
        distributionsTransportStart[j] = oldDistributionsTransportStart[j] + l_omegaD * (eqDistC[j] - oldDistributionsTransportStart[j]);
#endif
      }
    }
  });
}

/** \brief Collision step for coupled Thermal Transport Lattice-Boltzmann
 * \author Shota Ito
 * \date 13.06.22
 *
 */
template <MInt nDim, MInt nDist, class SysEqn>
template <MInt thermalMode>
void LbSolverDxQy<nDim, nDist, SysEqn>::bgki_thermal_and_transport_collision_step_base() {
  TRACE();

  const MInt pvrho = PV->RHO;
  const MInt pvt = PV->T;
  const MInt pvc = PV->C;

  constexpr MInt fldlen = Ld::dxQyFld();

  // Update the global Reynolds number if we use the tanh-adaption
  if(m_tanhInit && (globalTimeStep > m_initStartTime)) {
    const MFloat scale =
        0.5 * (m_tanhScaleFactor * tanh((5.0 * (MFloat)(globalTimeStep - m_initStartTime) / (MFloat)m_initTime) - 2.5))
        + 0.5;

    m_Re = m_initRe + scale * (m_finalRe - m_initRe);

    // final timestep reached
    if(globalTimeStep == m_initStartTime + m_initTime) {
      m_Re = m_finalRe;
      m_tanhInit = 0;
    }
  }

  // Update the according nu
  m_nu = m_Ma * LBCS / m_Re * m_referenceLength;

  // Calculate heat diffusion coefficient
  m_kappa = m_nu / m_Pr;
  m_diffusivity = m_nu * (m_Re / m_Pe);

#ifndef WAR_NVHPC_PSTL
  // inner loop references
  const MInt* const RESTRICT activeCellList = ALIGNED_I(m_activeCellList);
  MFloat* const RESTRICT cnu = ALIGNED_MF(&(a_nu(0)));
  MFloat* const RESTRICT ckappa = ALIGNED_MF(&(a_kappa(0)));
  MFloat* const RESTRICT cdiffusivity = ALIGNED_MF(&(a_diffusivity(0)));
  const MInt* const RESTRICT clevel = ALIGNED_I(&(this->m_cells.level(0)));
  MFloat* const RESTRICT oldDistributions = ALIGNED_MF(&(a_oldDistribution(0, 0)));
  MFloat* const RESTRICT distributions = ALIGNED_MF(&(a_distribution(0, 0)));
  MFloat* const RESTRICT oldDistributionsThermal = ALIGNED_MF(&(a_oldDistributionThermal(0, 0)));
  MFloat* const RESTRICT distributionsThermal = ALIGNED_MF(&(a_distributionThermal(0, 0)));
  MFloat* const RESTRICT oldDistributionsTransport = ALIGNED_MF(&(a_oldDistributionTransport(0, 0)));
  MFloat* const RESTRICT distributionsTransport = ALIGNED_MF(&(a_distributionTransport(0, 0)));
  MFloat* const RESTRICT variables = ALIGNED_MF(&(a_variable(0, 0)));

  // members and globals
  const MInt maxLevel_ = maxLevel();
  // distances
  const MInt distNu = &a_nu(1) - &a_nu(0);
  const MInt distKappa = &a_kappa(1) - &a_kappa(0);
  const MInt distDiffusivity = &a_diffusivity(1) - &a_diffusivity(0);
  const MInt distLevel = &(this->m_cells.level(1)) - &(this->m_cells.level(0));
  const MInt distVars = &(a_variable(1, 0)) - &(a_variable(0, 0));
#endif
  // Required for ported functions, since GPUs cannot access the global variable globalTimeStep
  const MInt gTS = globalTimeStep;

  const MFloat factor = (thermalMode == 1) ? ((nDim == 2) ? 3.0 : 18.0 / 5.0) : F6;

  // loop over all cells
  maia::parallelFor<true>(0, m_currentMaxNoCells, [=](MInt i) {
#ifdef WAR_NVHPC_PSTL
    const MInt pCellId = m_activeCellList[i];
    const MInt index = maxLevel() - this->a_level(pCellId);
#else
    const MInt pCellId = activeCellList[i];
    const MInt index = maxLevel_ - clevel[pCellId * distLevel];
#endif

    // perform collision on finest level in every timestep,
    // on the next coarser level every second timestep, and so on
    if((gTS - 1) % IPOW2(index) == 0) {
      // save nu for bndcnd and interface interpolation
#ifndef WAR_NVHPC_PSTL
      cnu[pCellId * distNu] = m_nu;
      ckappa[pCellId * distKappa] = m_kappa;
      cdiffusivity[pCellId * distDiffusivity] = m_diffusivity;

      const MInt distStart = pCellId * nDist;
      // varStart = (pCellId * distVars)  + (distVarOldVar *
      // !((MInt)floor((MFloat)localGlobalTimeStep/fpow2[index]) % 2));
      const MInt varStart = pCellId * distVars;

      MFloat* const oldDistributionsStart = &oldDistributions[distStart];
      MFloat* const distributionsStart = &distributions[distStart];
      MFloat* const oldDistributionsTStart = &oldDistributionsThermal[distStart];
      MFloat* const distributionsTStart = &distributionsThermal[distStart];
      MFloat* const oldDistributionsTransportStart = &oldDistributionsTransport[distStart];
      MFloat* const distributionsTransportStart = &distributionsTransport[distStart];

      MFloat* const variablesStart = &variables[varStart];
#else
      a_nu(pCellId) = m_nu;
      a_kappa(pCellId) = m_kappa;
      a_diffusivity(pCellId) = m_diffusivity;
#endif

      const MFloat tmp_FFPOW2 = FFPOW2(index);
      const MFloat l_omega = 2.0 / (1.0 + 6.0 * m_nu * tmp_FFPOW2);
      const MFloat l_omegaT = 2.0 / (1.0 + factor * m_kappa * tmp_FFPOW2);
      const MFloat l_omegaD = 2.0 / (1.0 + 6.0 * m_diffusivity * tmp_FFPOW2);

      swap_variables(pCellId);

      // Reinit macroscopic variables before recalculation
      MFloat l_v[nDim] = {0.0};
      MFloat l_rho = 0.0;
      MFloat l_t = 0.0;
      MFloat l_c = 0.0;

#ifdef WAR_NVHPC_PSTL
      // Add forcing term and calculate density rho
      for(MInt j = 0; j < nDist - 1; j++) {
        // Add forcing term and calculate density rho
        l_rho += a_oldDistribution(pCellId, j);
        l_t += a_oldDistributionThermal(pCellId, j);
        l_c += a_oldDistributionTransport(pCellId, j);
      }
      // Do not forget the last distribution
      l_rho += a_oldDistribution(pCellId, nDist - 1);
      l_t += a_oldDistributionThermal(pCellId, nDist - 1);
      l_c += a_oldDistributionTransport(pCellId, nDist - 1);
#else
      // Add forcing term and calculate density rho
      for(MInt j = 0; j < nDist - 1; j++) {
        l_rho += oldDistributionsStart[j];
        l_t += oldDistributionsTStart[j];
        l_c += oldDistributionsTransportStart[j];
      }
      // Do not forget the last distribution
      l_rho += oldDistributionsStart[nDist - 1];
      l_t += oldDistributionsTStart[nDist - 1];
      l_c += oldDistributionsTransportStart[nDist - 1];
#endif

      // Calculate macroscopic variables
#ifndef WAR_NVHPC_PSTL
      if constexpr(nDim == 2) {
        // WARNING!!! The velocity should be calculated as in 3d
        // However, result is correct
        // calculation of u
        l_v[0] = oldDistributionsStart[1] + oldDistributionsStart[4] + oldDistributionsStart[5]
                 - oldDistributionsStart[7] - oldDistributionsStart[6] - oldDistributionsStart[0];
        // calculation of v
        l_v[1] = oldDistributionsStart[7] + oldDistributionsStart[3] + oldDistributionsStart[4]
                 - oldDistributionsStart[6] - oldDistributionsStart[2] - oldDistributionsStart[5];
      } else {
#endif
        for(MInt j = 0; j < fldlen; j++) {
          for(MInt d = 0; d < nDim; d++) {
#ifdef WAR_NVHPC_PSTL
            l_v[d] += a_oldDistribution(pCellId, m_pFld[d * fldlen + j]);
            l_v[d] -= a_oldDistribution(pCellId, m_nFld[d * fldlen + j]);
#else
            l_v[d] += oldDistributionsStart[Ld::pFld(d, j)];
            l_v[d] -= oldDistributionsStart[Ld::nFld(d, j)];
#endif
          }
        }
#ifndef WAR_NVHPC_PSTL
      }
#endif

      std::array<MFloat, nDim> u{};
      for(MInt d = 0; d < nDim; d++) {
        u[d] = l_v[d] / l_rho;
      }

      // Save new macroscopic variables in cell
#ifdef WAR_NVHPC_PSTL
      for(MInt d = 0; d < nDim; d++) {
        a_variable(pCellId, d) = l_v[d];
      }
      a_variable(pCellId, pvrho) = l_rho;
      a_variable(pCellId, pvc) = l_c;
      IF_CONSTEXPR(thermalMode == 0) { a_variable(pCellId, pvt) = l_t; }
      IF_CONSTEXPR(thermalMode == 1) {
        const MFloat D = (MFloat)nDim;
        const MFloat l_ie_rho = l_t;
        l_t = l_ie_rho * F2 / (D * l_rho);
        a_variable(pCellId, pvt) = l_t;
      }
      IF_CONSTEXPR(thermalMode == 2) {
        const MFloat D = (MFloat)nDim;
        const MFloat sqVel = std::inner_product(&u[0], &u[nDim], &u[0], F0);
        const MFloat l_te_rho = l_t;
        l_t = (l_te_rho / l_rho - sqVel * F1B2) * F2 / D;
        a_variable(pCellId, pvt) = l_t;
        for(MInt d = 0; d < nDim; d++) {
          a_variable(pCellId, d) = l_v[d];
        }
      }
#else
      for(MInt d = 0; d < nDim; d++) {
        variablesStart[d] = l_v[d];
      }
      variablesStart[pvrho] = l_rho;
      variablesStart[pvc] = l_c;
      IF_CONSTEXPR(thermalMode == 0) { variablesStart[pvt] = l_t; }
      IF_CONSTEXPR(thermalMode == 1) {
        const MFloat D = (MFloat)nDim;
        const MFloat l_ie_rho = l_t;
        l_t = l_ie_rho * F2 / (D * l_rho);
        variablesStart[pvt] = l_t;
      }
      IF_CONSTEXPR(thermalMode == 2) {
        const MFloat D = (MFloat)nDim;
        const MFloat sqVel = std::inner_product(&u[0], &u[nDim], &u[0], F0);
        const MFloat l_te_rho = l_t;
        l_t = (l_te_rho / l_rho - sqVel * F1B2) * F2 / D;
        variablesStart[pvt] = l_t;
        for(MInt d = 0; d < nDim; d++) {
          variablesStart[d] = l_v[d];
        }
      }
#endif

      std::array<MFloat, nDist> eqDist;
      eqDist = getEqDists(l_rho, u.data());
      std::array<MFloat, nDist> eqDistT;
      eqDistT = getEqDistsThermal_<thermalMode>(l_t, l_rho, u.data());
      std::array<MFloat, nDist> eqDistC;
      eqDistC = getEqDistsTransport(l_c, u.data());

      std::array<MFloat, nDist> Sj;
      IF_CONSTEXPR(thermalMode == 2) {
        std::array<MFloat, 2 * nDim> b{};
        for(MInt d = 0; d < nDim; d++) {
          b[2 * d] = -u[d];
          b[2 * d + 1] = u[d];
        }
        const MFloat sqVel = std::inner_product(&u[0], &u[nDim], &u[0], F0);
#ifdef WAR_NVHPC_PSTL
        for(MInt j = 0; j < m_distFld[0]; j++) {
          Sj[j] = (l_omegaT - l_omega) * (b[j] - sqVel * F1B2) * (a_oldDistribution(pCellId, j) - eqDist[j]);
        }
        for(MInt j = 0; j < m_distFld[1]; j++) {
          const MInt p = 2 * j;
          const MInt pos = j + m_distFld[0];
          const MFloat tmp = (b[m_mFld1[p]] + b[m_mFld1[p + 1]]);
          Sj[pos] = (l_omegaT - l_omega) * (tmp - sqVel * F1B2) * (a_oldDistribution(pCellId, pos) - eqDist[pos]);
        }
        for(MInt j = 0; j < m_distFld[2]; j++) {
          const MInt p = 3 * j;
          const MInt pos = j + m_distFld[0] + m_distFld[1];
          const MFloat tmp = (b[m_mFld2[p]] + b[m_mFld2[p + 1]] + b[m_mFld2[p + 2]]);
          Sj[pos] = (l_omegaT - l_omega) * (tmp - sqVel * F1B2) * (a_oldDistribution(pCellId, pos) - eqDist[pos]);
        }
        Sj[Ld::lastId()] = (l_omegaT - l_omega) * (F0 - sqVel * F1B2)
                           * (a_oldDistribution(pCellId, Ld::lastId()) - eqDist[Ld::lastId()]);
#else
        for(MInt j = 0; j < Ld::distFld(0); j++) {
          Sj[j] = (l_omegaT - l_omega) * (b[j] - sqVel * F1B2) * (oldDistributionsStart[j] - eqDist[j]);
        }
        for(MInt j = 0; j < Ld::distFld(1); j++) {
          const MInt p = 2 * j;
          const MInt pos = j + Ld::distFld(0);
          const MFloat tmp = (b[Ld::mFld1(p)] + b[Ld::mFld1(p + 1)]);
          Sj[pos] = (l_omegaT - l_omega) * (tmp - sqVel * F1B2) * (oldDistributionsStart[pos] - eqDist[pos]);
        }
        for(MInt j = 0; j < Ld::distFld(2); j++) {
          const MInt p = 3 * j;
          const MInt pos = j + Ld::distFld(0) + Ld::distFld(1);
          const MFloat tmp = (b[Ld::mFld2(p)] + b[Ld::mFld2(p + 1)] + b[Ld::mFld2(p + 2)]);
          Sj[pos] = (l_omegaT - l_omega) * (tmp - sqVel * F1B2) * (oldDistributionsStart[pos] - eqDist[pos]);
        }
        Sj[Ld::lastId()] =
            (l_omegaT - l_omega) * (F0 - sqVel * F1B2) * (oldDistributionsStart[Ld::lastId()] - eqDist[Ld::lastId()]);
#endif
      }

      // Calculation of new distributions
      for(MInt j = 0; j < nDist; j++) {
#ifdef WAR_NVHPC_PSTL
        a_distribution(pCellId, j) =
            a_oldDistribution(pCellId, j) + l_omega * (eqDist[j] - a_oldDistribution(pCellId, j));
        a_distributionThermal(pCellId, j) =
            a_oldDistributionThermal(pCellId, j) + l_omegaT * (eqDistT[j] - a_oldDistributionThermal(pCellId, j));
        a_distributionTransport(pCellId, j) =
            a_oldDistributionTransport(pCellId, j) + l_omegaD * (eqDistC[j] - a_oldDistributionTransport(pCellId, j));
        IF_CONSTEXPR(thermalMode == 2) { a_distributionThermal(pCellId, j) += Sj[j]; }
#else
        distributionsStart[j] = oldDistributionsStart[j] + l_omega * (eqDist[j] - oldDistributionsStart[j]);
        distributionsTStart[j] = oldDistributionsTStart[j] + l_omegaT * (eqDistT[j] - oldDistributionsTStart[j]);
        IF_CONSTEXPR(thermalMode == 2) { distributionsTStart[j] += Sj[j]; }
        distributionsTransportStart[j] = oldDistributionsTransportStart[j] + l_omegaD * (eqDistC[j] - oldDistributionsTransportStart[j]);
#endif
      }
    }
  });
}

/** \brief Collision step for coupled Thermal Transport Lattice-Boltzmann
 * \author Shota Ito
 * \date 13.06.22
 *
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbSolverDxQy<nDim, nDist, SysEqn>::bgkc_thermal_transport_collision_step() {
  TRACE();

  bgki_thermal_and_transport_collision_step_base<0>();
}

/** \brief Collision step for coupled Thermal Transport Lattice-Boltzmann
 * \author Moritz Waldmann
 * \date 09.10.23
 *
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbSolverDxQy<nDim, nDist, SysEqn>::bgkc_innerenergy_transport_collision_step() {
  TRACE();

  bgki_thermal_and_transport_collision_step_base<1>();
}

/** \brief Collision step for coupled Thermal Transport Lattice-Boltzmann
 * \author Moritz Waldmann
 * \date 09.10.23
 *
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbSolverDxQy<nDim, nDist, SysEqn>::bgkc_totalenergy_transport_collision_step() {
  TRACE();

  bgki_thermal_and_transport_collision_step_base<2>();
}

template <MInt nDim, MInt nDist, class SysEqn>
template <MInt thermalMode>
void LbSolverDxQy<nDim, nDist, SysEqn>::bgki_thermal_collision_step_base() {
  TRACE();

  const MInt pvrho = PV->RHO;
  const MInt pvt = PV->T;

  constexpr MInt fldlen = Ld::dxQyFld();

  // Update the global Reynolds number if we use the tanh-adaption
  if(m_tanhInit && (globalTimeStep > m_initStartTime)) {
    const MFloat scale =
        0.5 * (m_tanhScaleFactor * tanh((5.0 * (MFloat)(globalTimeStep - m_initStartTime) / (MFloat)m_initTime) - 2.5))
        + 0.5;

    m_Re = m_initRe + scale * (m_finalRe - m_initRe);

    // final timestep reached
    if(globalTimeStep == m_initStartTime + m_initTime) {
      m_Re = m_finalRe;
      m_tanhInit = 0;
    }
  }

  // Update the according nu
  m_nu = m_Ma * LBCS / m_Re * m_referenceLength;

  // Calculate heat diffusion coefficient
  m_kappa = m_nu / m_Pr;

#ifndef WAR_NVHPC_PSTL
  // inner loop references
  const MInt* const RESTRICT activeCellList = ALIGNED_I(m_activeCellList);
  MFloat* const RESTRICT cnu = ALIGNED_MF(&(a_nu(0)));
  MFloat* const RESTRICT ckappa = ALIGNED_MF(&(a_kappa(0)));
  const MInt* const RESTRICT clevel = ALIGNED_I(&(this->m_cells.level(0)));
  MFloat* const RESTRICT oldDistributions = ALIGNED_MF(&(a_oldDistribution(0, 0)));
  MFloat* const RESTRICT distributions = ALIGNED_MF(&(a_distribution(0, 0)));
  MFloat* const RESTRICT oldDistributionsThermal = ALIGNED_MF(&(a_oldDistributionThermal(0, 0)));
  MFloat* const RESTRICT distributionsThermal = ALIGNED_MF(&(a_distributionThermal(0, 0)));
  MFloat* const RESTRICT variables = ALIGNED_MF(&(a_variable(0, 0)));
  const MFloat* const Fp = ALIGNED_MF(m_Fext);
  // members and globals
  const MInt maxLevel_ = maxLevel();
  // distances
  const MInt distNu = &a_nu(1) - &a_nu(0);
  const MInt distKappa = &a_kappa(1) - &a_kappa(0);
  const MInt distLevel = &(this->m_cells.level(1)) - &(this->m_cells.level(0));
  const MInt distVars = &(a_variable(1, 0)) - &(a_variable(0, 0));
#endif
  // Required for ported functions, since GPUs cannot access the global variable globalTimeStep
  const MInt gTS = globalTimeStep;

  const MFloat factor = (thermalMode == 1) ? ((nDim == 2) ? 3.0 : 18.0 / 5.0) : F6;

  // loop over all cells
  maia::parallelFor<true>(0, m_currentMaxNoCells, [=](MInt i) {
#ifdef WAR_NVHPC_PSTL
    const MInt pCellId = m_activeCellList[i];
    const MInt index = maxLevel() - this->a_level(pCellId);
#else
    const MInt pCellId = activeCellList[i];
    const MInt index = maxLevel_ - clevel[pCellId * distLevel];
#endif

    // perform collision on finest level in every timestep,
    // on the next coarser level every second timestep, and so on
    if((gTS - 1) % IPOW2(index) == 0) {
      // save nu for bndcnd and interface interpolation
#ifndef WAR_NVHPC_PSTL
      cnu[pCellId * distNu] = m_nu;
      ckappa[pCellId * distKappa] = m_kappa;

      const MInt distStart = pCellId * nDist;
      // varStart = (pCellId * distVars)  + (distVarOldVar *
      // !((MInt)floor((MFloat)localGlobalTimeStep/fpow2[index]) % 2));
      const MInt varStart = pCellId * distVars;

      MFloat* const oldDistributionsStart = &oldDistributions[distStart];
      MFloat* const distributionsStart = &distributions[distStart];
      MFloat* const oldDistributionsTStart = &oldDistributionsThermal[distStart];
      MFloat* const distributionsTStart = &distributionsThermal[distStart];

      MFloat* const variablesStart = &variables[varStart];
#else
      a_nu(pCellId) = m_nu;
      a_kappa(pCellId) = m_kappa;
#endif

      const MFloat tmp_FPOW2 = FPOW2(index);
      const MFloat tmp_FFPOW2 = FFPOW2(index);
      const MFloat l_omega = 2.0 / (1.0 + 6.0 * m_nu * tmp_FFPOW2);
      const MFloat l_omegaT = 2.0 / (1.0 + factor * m_kappa * tmp_FFPOW2);

      swap_variables(pCellId);

      // Reinit macroscopic variables before recalculation
      MFloat l_v[nDim] = {0.0};
      MFloat l_rho = 0.0;
      MFloat l_t = 0.0;

#ifdef WAR_NVHPC_PSTL
      // Add forcing term and calculate density rho
      for(MInt j = 0; j < nDist - 1; j++) {
        IF_CONSTEXPR(thermalMode != 2) a_oldDistribution(pCellId, j) += tmp_FPOW2 * m_Fext[j];
        // Add forcing term and calculate density rho
        l_rho += a_oldDistribution(pCellId, j);
        l_t += a_oldDistributionThermal(pCellId, j);
      }
      // Do not forget the last distribution
      l_rho += a_oldDistribution(pCellId, nDist - 1);
      l_t += a_oldDistributionThermal(pCellId, nDist - 1);
#else
      // Add forcing term and calculate density rho
      for(MInt j = 0; j < nDist - 1; j++) {
        IF_CONSTEXPR(thermalMode != 2) oldDistributionsStart[j] += tmp_FPOW2 * Fp[j];
        l_rho += oldDistributionsStart[j];
        l_t += oldDistributionsTStart[j];
      }
      // Do not forget the last distribution
      l_rho += oldDistributionsStart[nDist - 1];
      l_t += oldDistributionsTStart[nDist - 1];
#endif

      // Calculate macroscopic variables
#ifndef WAR_NVHPC_PSTL
      if constexpr(nDim == 2) {
        // WARNING!!! The velocity should be calculated as in 3d
        // However, result is correct
        // calculation of u
        l_v[0] = oldDistributionsStart[1] + oldDistributionsStart[4] + oldDistributionsStart[5]
                 - oldDistributionsStart[7] - oldDistributionsStart[6] - oldDistributionsStart[0];
        // calculation of v
        l_v[1] = oldDistributionsStart[7] + oldDistributionsStart[3] + oldDistributionsStart[4]
                 - oldDistributionsStart[6] - oldDistributionsStart[2] - oldDistributionsStart[5];
      } else {
#endif
        for(MInt j = 0; j < fldlen; j++) {
          for(MInt d = 0; d < nDim; d++) {
#ifdef WAR_NVHPC_PSTL
            l_v[d] += a_oldDistribution(pCellId, m_pFld[d * fldlen + j]);
            l_v[d] -= a_oldDistribution(pCellId, m_nFld[d * fldlen + j]);
#else
            l_v[d] += oldDistributionsStart[Ld::pFld(d, j)];
            l_v[d] -= oldDistributionsStart[Ld::nFld(d, j)];
#endif
          }
        }
#ifndef WAR_NVHPC_PSTL
      }
#endif

      std::array<MFloat, nDim> u{};
      for(MInt d = 0; d < nDim; d++) {
        u[d] = l_v[d] / l_rho;
      }

      // Save new macroscopic variables in cell
#ifdef WAR_NVHPC_PSTL
      for(MInt d = 0; d < nDim; d++) {
        a_variable(pCellId, d) = l_v[d];
      }
      a_variable(pCellId, pvrho) = l_rho;
      IF_CONSTEXPR(thermalMode == 0) { a_variable(pCellId, pvt) = l_t; }
      IF_CONSTEXPR(thermalMode == 1) {
        const MFloat D = (MFloat)nDim;
        const MFloat l_ie_rho = l_t;
        l_t = l_ie_rho * F2 / (D * l_rho);
        a_variable(pCellId, pvt) = l_t;
      }
      IF_CONSTEXPR(thermalMode == 2) {
        const MFloat accelerationForce = (this->m_externalForcing) ? m_densityGradient * F1B3 : F0;
        const MFloat D = (MFloat)nDim;
        // Set the external acceleratrion vector
        // So far only for a channel in x-direction
        std::array<MFloat, nDim> a{};
        a[0] = accelerationForce;
        for(MInt d = 0; d < nDim; d++) {
          l_v[d] += l_rho * tmp_FPOW2 * a[d] * F1B2;
          u[d] = l_v[d] / l_rho;
        }
        const MFloat sqVel = std::inner_product(&u[0], &u[nDim], &u[0], F0);

        MFloat u_x_a = F0;
        for(MInt d = 0; d < nDim; d++) {
          u_x_a += u[d] * a[d];
        }
        const MFloat l_te_rho = l_t + F1B2 * tmp_FPOW2 * l_rho * u_x_a;
        l_t = (l_te_rho / l_rho - sqVel * F1B2) * F2 / D;
        a_variable(pCellId, pvt) = l_t;
        for(MInt d = 0; d < nDim; d++) {
          a_variable(pCellId, d) = l_v[d];
        }
      }
#else
      for(MInt d = 0; d < nDim; d++) {
        variablesStart[d] = l_v[d];
      }
      variablesStart[pvrho] = l_rho;
      IF_CONSTEXPR(thermalMode == 0) { variablesStart[pvt] = l_t; }
      IF_CONSTEXPR(thermalMode == 1) {
        const MFloat D = (MFloat)nDim;
        const MFloat l_ie_rho = l_t;
        l_t = l_ie_rho * F2 / (D * l_rho);
        variablesStart[pvt] = l_t;
      }
      IF_CONSTEXPR(thermalMode == 2) {
        const MFloat accelerationForce = (this->m_externalForcing) ? m_densityGradient * F1B3 : F0;
        const MFloat D = (MFloat)nDim;
        // Set the external acceleratrion vector
        // So far only for a channel in x-direction
        std::array<MFloat, nDim> a{};
        a[0] = accelerationForce;
        for(MInt d = 0; d < nDim; d++) {
          l_v[d] += l_rho * tmp_FPOW2 * a[d] * F1B2;
          u[d] = l_v[d] / l_rho;
        }
        const MFloat sqVel = std::inner_product(&u[0], &u[nDim], &u[0], F0);

        MFloat u_x_a = F0;
        for(MInt d = 0; d < nDim; d++) {
          u_x_a += u[d] * a[d];
        }
        const MFloat l_te_rho = l_t + F1B2 * tmp_FPOW2 * l_rho * u_x_a;
        l_t = (l_te_rho / l_rho - sqVel * F1B2) * F2 / D;
        variablesStart[pvt] = l_t;
        for(MInt d = 0; d < nDim; d++) {
          variablesStart[d] = l_v[d];
        }
      }
#endif

      std::array<MFloat, nDist> eqDist;
      eqDist = getEqDists(l_rho, u.data());
      std::array<MFloat, nDist> eqDistT;
      eqDistT = getEqDistsThermal_<thermalMode>(l_t, l_rho, u.data());

      std::array<MFloat, nDist> forcing;
      std::array<MFloat, nDist> forcingT;
      std::array<MFloat, nDist> Sj;
      IF_CONSTEXPR(thermalMode == 2) {
        const MFloat accelerationForce = (this->m_externalForcing) ? m_densityGradient * F1B3 : F0;
        std::array<MFloat, 2 * nDim> a{};
        a[0] = -accelerationForce;
        a[1] = accelerationForce;
        MFloat u_x_a = F0;
        std::array<MFloat, 2 * nDim> b;
        for(MInt d = 0; d < nDim; d++) {
          b[2 * d] = -u[d];
          b[2 * d + 1] = u[d];
          u_x_a += b[d * 2 + 1] * a[d * 2 + 1];
        }
        const MFloat D = (MFloat)nDim;
        const MFloat sqVel = std::inner_product(&u[0], &u[nDim], &u[0], F0);
        const MFloat l_te_rho = (l_t * D * F1B2 + sqVel * F1B2) * l_rho;
#ifdef WAR_NVHPC_PSTL
        for(MInt j = 0; j < m_distFld[0]; j++) {
          forcing[j] = m_tp[1] * F1BCSsq * l_rho * (a[j] + a[j] * b[j] * F1BCSsq - u_x_a);
          Sj[j] = (l_omegaT - l_omega) * (b[j] - sqVel * F1B2)
                  * (a_oldDistribution(pCellId, j) - eqDist[j] + F1B2 * tmp_FPOW2 * forcing[j]);
          forcingT[j] = m_tp[1] * F1BCSsq * l_te_rho * a[j] + a_oldDistribution(pCellId, j) * a[j];
          forcing[j] *= tmp_FPOW2 * (F1 - l_omega * F1B2);
          forcingT[j] *= tmp_FPOW2 * (F1 - l_omegaT * F1B2);
        }
        for(MInt j = 0; j < m_distFld[1]; j++) {
          const MInt p = 2 * j;
          const MInt pos = j + m_distFld[0];
          const MFloat tmp = (b[m_mFld1[p]] + b[m_mFld1[p + 1]]);
          const MFloat tmpT = (a[m_mFld1[p]] + a[m_mFld1[p + 1]]);
          forcing[pos] = m_tp[2] * F1BCSsq * l_rho * (tmpT + tmpT * tmp * F1BCSsq - u_x_a);
          Sj[pos] = (l_omegaT - l_omega) * (tmp - sqVel * F1B2)
                    * (a_oldDistribution(pCellId, pos) - eqDist[pos] + F1B2 * tmp_FPOW2 * forcing[pos]);
          forcingT[pos] = m_tp[2] * F1BCSsq * l_te_rho * tmpT + a_oldDistribution(pCellId, pos) * tmpT;
          forcing[pos] *= tmp_FPOW2 * (F1 - l_omega * F1B2);
          forcingT[pos] *= tmp_FPOW2 * (F1 - l_omegaT * F1B2);
        }
        for(MInt j = 0; j < m_distFld[2]; j++) {
          const MInt p = 3 * j;
          const MInt pos = j + m_distFld[0] + m_distFld[1];
          const MFloat tmp = (b[m_mFld2[p]] + b[m_mFld2[p + 1]] + b[m_mFld2[p + 2]]);
          const MFloat tmpT = (a[m_mFld2[p]] + a[m_mFld2[p + 1]] + a[m_mFld2[p + 2]]);
          forcing[pos] = m_tp[3] * F1BCSsq * l_rho * (tmpT + tmpT * tmp * F1BCSsq - u_x_a);
          Sj[pos] = (l_omegaT - l_omega) * (tmp - sqVel * F1B2)
                    * (a_oldDistribution(pCellId, pos) - eqDist[pos] + F1B2 * tmp_FPOW2 * forcing[pos]);
          forcingT[pos] = m_tp[3] * F1BCSsq * l_te_rho * tmpT + a_oldDistribution(pCellId, pos) * tmpT;
          forcing[pos] *= tmp_FPOW2 * (F1 - l_omega * F1B2);
          forcingT[pos] *= tmp_FPOW2 * (F1 - l_omegaT * F1B2);
        }
        forcing[Ld::lastId()] = m_tp[0] * l_rho * (-F1BCSsq * u_x_a);
        Sj[Ld::lastId()] = (l_omegaT - l_omega) * (F0 - sqVel * F1B2)
                           * (a_oldDistribution(pCellId, Ld::lastId()) - eqDist[Ld::lastId()]
                              + F1B2 * tmp_FPOW2 * forcing[Ld::lastId()]);
        forcingT[Ld::lastId()] = F0;
        forcing[Ld::lastId()] *= tmp_FPOW2 * (F1 - l_omega * F1B2);
#else
        for(MInt j = 0; j < Ld::distFld(0); j++) {
          forcing[j] = Ld::tp(1) * F1BCSsq * l_rho * (a[j] + a[j] * b[j] * F1BCSsq - u_x_a);
          Sj[j] = (l_omegaT - l_omega) * (b[j] - sqVel * F1B2)
                  * (oldDistributionsStart[j] - eqDist[j] + F1B2 * tmp_FPOW2 * forcing[j]);
          forcingT[j] = Ld::tp(1) * F1BCSsq * l_te_rho * a[j] + oldDistributionsStart[j] * a[j];
          forcing[j] *= tmp_FPOW2 * (F1 - l_omega * F1B2);
          forcingT[j] *= tmp_FPOW2 * (F1 - l_omegaT * F1B2);
        }
        for(MInt j = 0; j < Ld::distFld(1); j++) {
          const MInt p = 2 * j;
          const MInt pos = j + Ld::distFld(0);
          const MFloat tmp = (b[Ld::mFld1(p)] + b[Ld::mFld1(p + 1)]);
          const MFloat tmpT = (a[Ld::mFld1(p)] + a[Ld::mFld1(p + 1)]);
          forcing[pos] = Ld::tp(2) * F1BCSsq * l_rho * (tmpT + tmpT * tmp * F1BCSsq - u_x_a);
          Sj[pos] = (l_omegaT - l_omega) * (tmp - sqVel * F1B2)
                    * (oldDistributionsStart[pos] - eqDist[pos] + F1B2 * tmp_FPOW2 * forcing[pos]);
          forcingT[pos] = Ld::tp(2) * F1BCSsq * l_te_rho * tmpT + oldDistributionsStart[pos] * tmpT;
          forcing[pos] *= tmp_FPOW2 * (F1 - l_omega * F1B2);
          forcingT[pos] *= tmp_FPOW2 * (F1 - l_omegaT * F1B2);
        }
        for(MInt j = 0; j < Ld::distFld(2); j++) {
          const MInt p = 3 * j;
          const MInt pos = j + Ld::distFld(0) + Ld::distFld(1);
          const MFloat tmp = (b[Ld::mFld2(p)] + b[Ld::mFld2(p + 1)] + b[Ld::mFld2(p + 2)]);
          const MFloat tmpT = (a[Ld::mFld2(p)] + a[Ld::mFld2(p + 1)] + a[Ld::mFld2(p + 2)]);
          forcing[pos] = Ld::tp(3) * F1BCSsq * l_rho * (tmpT + tmpT * tmp * F1BCSsq - u_x_a);
          Sj[pos] = (l_omegaT - l_omega) * (tmp - sqVel * F1B2)
                    * (oldDistributionsStart[pos] - eqDist[pos] + F1B2 * tmp_FPOW2 * forcing[pos]);
          forcingT[pos] = Ld::tp(3) * F1BCSsq * l_te_rho * tmpT + oldDistributionsStart[pos] * tmpT;
          forcing[pos] *= tmp_FPOW2 * (F1 - l_omega * F1B2);
          forcingT[pos] *= tmp_FPOW2 * (F1 - l_omegaT * F1B2);
        }
        forcing[Ld::lastId()] = Ld::tp(0) * l_rho * (-F1BCSsq * u_x_a);
        Sj[Ld::lastId()] =
            (l_omegaT - l_omega) * (F0 - sqVel * F1B2)
            * (oldDistributionsStart[Ld::lastId()] - eqDist[Ld::lastId()] + F1B2 * tmp_FPOW2 * forcing[Ld::lastId()]);
        forcingT[Ld::lastId()] = F0;
        forcing[Ld::lastId()] *= tmp_FPOW2 * (F1 - l_omega * F1B2);
#endif
      }

      // Calculation of new distributions
      for(MInt j = 0; j < nDist; j++) {
#ifdef WAR_NVHPC_PSTL
        a_distribution(pCellId, j) =
            a_oldDistribution(pCellId, j) + l_omega * (eqDist[j] - a_oldDistribution(pCellId, j));
        a_distributionThermal(pCellId, j) =
            a_oldDistributionThermal(pCellId, j) + l_omegaT * (eqDistT[j] - a_oldDistributionThermal(pCellId, j));
        IF_CONSTEXPR(thermalMode == 2) {
          a_distribution(pCellId, j) += forcing[j];
          a_distributionThermal(pCellId, j) += forcingT[j] + Sj[j];
        }
#else
        distributionsStart[j] = oldDistributionsStart[j] + l_omega * (eqDist[j] - oldDistributionsStart[j]);
        distributionsTStart[j] = oldDistributionsTStart[j] + l_omegaT * (eqDistT[j] - oldDistributionsTStart[j]);
        IF_CONSTEXPR(thermalMode == 2) {
          distributionsStart[j] += forcing[j];
          distributionsTStart[j] += forcingT[j] + Sj[j];
        }
#endif
      }
    }
  });
}

/** \fn LbSolverDxQy::bgki_thermal_collision_step()
 *  \author Ludwig Eduard Boltzmann
 *  \date   20.2.1844
 *  \propVal{solverMethod,MAIA_LATTICE_BGK_THERMAL}
 *
 *  Collision step of the incompressible LBGK algorithm for thermal LB<br>
 *  \collstep{LbSolverDxQy::bgki_thermal_collision_step(), LBCbgki_thermal, Thermal Incompressible BGK Collision Step}
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbSolverDxQy<nDim, nDist, SysEqn>::bgki_thermal_collision_step() {
  TRACE();

  bgki_thermal_collision_step_base<0>();
}


/** \fn LbSolverDxQy::bgki_innerEnergy_collision_step()
 * \author Moritz Waldmann
 * \date 23.02.2019
 *  \propVal{solverMethod,MAIA_LATTICE_BGK_INNERENERGY}
 *
 * Collision step of the incompressible LBGK algorithm with TLBGK extension using an inner energy distribution
 * function approach
 *
 * Performs the collision step for the BGK method with TLBGK extension:
 *
 * <ul>
 *  <li>the macroscopic variables are obtained by evaluating the moments of the PPDFs.</li>
 *  <li>seperate local velocities are calculated by dividing by the density \f$\rho\f$.</li>
 *  <li>Recalculation of the PPDFs for all directions for the LBGK and TLBGK. Note, that the density \f$rho\f$ has been
 * moved outside the inner brackets of the equilibrium distribution function and that different collision frequencies
 * \f$\omega\f$ and \f$\omega_T\f$ are used for the LBGK and the TLBGK, respectively. The velocities of the TLBGK are
 * used from that one obtained from the moments of the LBGK.
 *  </li>
 * </ul><br>
 *
 * \todo labels:LB natural convection
 * \collstep{LbSolverDxQy::bgki_innerEnergy_collision_step(), LBCbgki_innerEn, Therm. BGKI Collision Step (innerEnergy)}
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbSolverDxQy<nDim, nDist, SysEqn>::bgki_innerEnergy_collision_step() {
  TRACE();

  bgki_thermal_collision_step_base<1>();
}

/** \fn LbSolverDxQy::bgki_totalEnergy_collision_step()
 * \propVal{solverMethod,MAIA_LATTICE_BGK_TOTALENERGY}
 *
 * Collision step of the incompressible LBGK algorithm with TLBGK extension using a total  energy distribution
 * function approach
 *
 * Performs the collision step for the BGK method with TLBGK extension:
 *
 * <ul>
 *  <li>the macroscopic variables are obtained by evaluating the moments of the PPDFs.</li>
 *  <li>seperate local velocities are calculated by dividing by the density \f$\rho\f$.</li>
 *  <li>Recalculation of the PPDFs for all directions for the LBGK and TLBGK. Note, that the density \f$rho\f$ has been
 * moved outside the inner brackets of the equilibrium distribution function and that different collision frequencies
 * \f$\omega\f$ and \f$\omega_T\f$ are used for the LBGK and the TLBGK, respectively. The velocities of the TLBGK are
 * used from that one obtained from the moments of the LBGK.
 *  </li>
 * </ul><br>
 * \collstep{LbSolverDxQy::bgki_totalEnergy_collision_step(), LBCbgki_totalEn, Therm. BGKI Collision Step (totalEnergy)}
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbSolverDxQy<nDim, nDist, SysEqn>::bgki_totalEnergy_collision_step() {
  TRACE();

  bgki_thermal_collision_step_base<2>();
}

/** \fn LbSolverDxQy::bgki_smagorinsky_collision_step()
 * \propVal{solverMethod,MAIA_LATTICE_BGKI_SMAGORINSKY}
 *
 * Collision step of the incompressible LBGK algorithm with SGS modelling
 *
 * only for D3Q19<br>
 * \collstep{LbSolverDxQy::bgki_smagorinsky_collision_step(), LBCsmago_bgki, Smagorinsky BGKI Collision Step}
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbSolverDxQy<nDim, nDist, SysEqn>::bgki_smagorinsky_collision_step() {
  bgki_smagorinsky_collision_step_base<0>();
}

/** \fn LbSolverDxQy::bgki_smago_wall_collision_step()
 * \propVal{solverMethod,MAIA_LATTICE_BGKI_SMAGO_WALL}
 *
 * Collision step of the incompressible LBGK algorithm with SGS modelling
 * and van Driest damping function
 *
 * only for D3Q19 and channel flow<br>
 * \collstep{LbSolverDxQy::bgki_smago_wall_collision_step(), LBCsmago_wall_bgki, Smagorinsky BGKI Collision Step (Wall)}
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbSolverDxQy<nDim, nDist, SysEqn>::bgki_smago_wall_collision_step() {
  bgki_smagorinsky_collision_step_base<1>();
}

/** \fn LbSolverDxQy::bgki_smagorinsky_collision_step2()
 * \propVal{solverMethod,MAIA_LATTICE_BGKI_SMAGORINSKY2}
 *
 * Collision step of the incompressible LBGK algorithm with SGS modelling and calculation of sgs terms along y
 * axis
 *
 * only for D3Q19<br>
 * \collstep{LbSolverDxQy::bgki_smagorinsky_collision_step2(), LBCsmago_bgki2, Smagorinsky BGKI Collision Step (2)}
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbSolverDxQy<nDim, nDist, SysEqn>::bgki_smagorinsky_collision_step2() {
  bgki_smagorinsky_collision_step_base<2>();
}

template <MInt nDim, MInt nDist, class SysEqn>
template <MInt mode>
void LbSolverDxQy<nDim, nDist, SysEqn>::bgki_smagorinsky_collision_step_base() {
  TRACE();

  // For now testing only the D3Q19 algorithm
  MFloat rho = F0;
  MFloat Q = F0;
  MFloat tau = F0;
  MFloat S = F0;

  MInt index = 0;
  MFloat actualCellLength = F0;
  constexpr MInt direction = 1;

  MChar buf[10];
  MString fileName;
  ofstream ofl;

  ScratchSpace<MFloat> diss(m_arraySize[direction], AT_, "diss");
  ScratchSpace<MFloat> totalDiss(m_arraySize[direction], AT_, "totalDiss");
  ScratchSpace<MFloat> SGSDiss(m_arraySize[direction], AT_, "SGSDiss");
  ScratchSpace<MFloat> totalSGSDiss(m_arraySize[direction], AT_, "totalSGSDiss");
  MFloat energy;

  // reset dissipation and energy values
  for(MInt i = 0; i < m_arraySize[direction]; i++) {
    diss[i] = F0;
    SGSDiss[i] = F0;
  }

  energy = F0;

  if(m_tanhInit && (globalTimeStep > m_initStartTime)) {
    MFloat scale =
        0.5 * (m_tanhScaleFactor * tanh((5.0 * (MFloat)(globalTimeStep - m_initStartTime) / (MFloat)m_initTime) - 2.5))
        + 0.5;

    m_Re = m_initRe + scale * (m_finalRe - m_initRe);

    // final timestep reached
    if(globalTimeStep == m_initStartTime + m_initTime) {
      m_Re = m_finalRe;
      m_tanhInit = 0;
    }
  }
  m_nu = m_Ma * LBCS / m_Re * m_referenceLength;

  for(MInt i = 0; i < m_currentMaxNoCells; i++) {
    const MInt pCellId = m_activeCellList[i];

    if((globalTimeStep - 1) % IPOW2(maxLevel() - this->a_level(pCellId))
       == 0) { // perform collision on finest level in every timestep, // on the next coarser level every second
               // timestep, and so on

      swap_variables(pCellId);

      // Calculate macroscopic variables
      a_variable(pCellId, PV->RHO) = 0.0;
      for(MInt j = 0; j < nDist; j++) {
        a_variable(pCellId, PV->RHO) += a_oldDistribution(pCellId, j);
      }
      rho = a_variable(pCellId, PV->RHO);

      MFloat u[nDim];
      for(MInt j = 0; j < nDim; j++) {
        u[j] = F0;
      }

      if constexpr(nDim == 2) {
        // warning!!! the velocity should be calculated as in 3d
        // however, result is correct
        // calculation of u
        u[0] = a_oldDistribution(pCellId, 1) + a_oldDistribution(pCellId, 4) + a_oldDistribution(pCellId, 5)
               - a_oldDistribution(pCellId, 7) - a_oldDistribution(pCellId, 6) - a_oldDistribution(pCellId, 0);
        // calculation of v
        u[1] = a_oldDistribution(pCellId, 7) + a_oldDistribution(pCellId, 3) + a_oldDistribution(pCellId, 4)
               - a_oldDistribution(pCellId, 6) - a_oldDistribution(pCellId, 2) - a_oldDistribution(pCellId, 5);
      } else {
        for(MInt j = 0; j < Ld::dxQyFld(); j++) {
          for(MInt d = 0; d < nDim; d++) {
            u[d] += a_oldDistribution(pCellId, Ld::pFld(d, j));
            u[d] -= a_oldDistribution(pCellId, Ld::nFld(d, j));
          }
        }
      }

      for(MInt j = 0; j < nDim; j++) {
        a_variable(pCellId, j) = u[j];
      }

      // Calculation of equilibrium and non-equilibrium distribution
      std::array<MFloat, nDist> eqDist;
      eqDist = getEqDists(rho, u);

      MFloat nonEqDist[nDist];
      for(MInt j = 0; j < nDist; j++) {
        nonEqDist[j] = a_oldDistribution(pCellId, j) - eqDist[j];
      }

      //--------------------------------------------
      // Calculation of overall viscosity
      // according to Hu, Sterling, Chen 1994

      // 1. Calculation of original relaxation time (on current level)
      tau = F1B2 + 3.0 * m_nu * FFPOW2(maxLevel() - this->a_level(pCellId));

      // 2. Calculation of momentum flux tensor from non-equilibrium parts
      calculateMomentumFlux(pCellId);

      // 3. Calculation of the filtered mean momentum flux
      Q = 0.0;
      for(MInt k = 0; k < nDim * nDim; k++) {
        Q += m_momentumFlux[pCellId][k] * m_momentumFlux[pCellId][k];
      }
      Q = sqrt(2.0 * Q);

      // 4. Calculation of new relaxation time
      if(mode == 1) {
        const MFloat aPlus = 26.0;
        MFloat tmpWidth = m_referenceLength * m_smallestCellLength;
        MFloat yPlus = m_ReTau * (1.0 - fabs(a_coordinate(pCellId, 1)) / tmpWidth);
        MFloat lambda = m_Cs * m_deltaX * (F1 - exp(-yPlus / aPlus));
        tau += F1B2
               * (sqrt(tau * tau
                       + 2.0 * SQRT2 * lambda * lambda * (F1BCSsq * F1BCSsq) * Q
                             * FPOW2(maxLevel() - this->a_level(pCellId)))
                  - tau);
      } else {
        tau += F1B2
               * (sqrt(tau * tau
                       + 2.0 * SQRT2 * m_Cs * m_Cs * m_deltaX * m_deltaX * (F1BCSsq * F1BCSsq) * Q
                             * FPOW2(maxLevel() - this->a_level(pCellId)))
                  - tau);
      }

      m_omega = 1.0 / tau;

      // save total nu for bndcnd and interface interpolation (save on finest level)
      a_nu(pCellId) = FPOW2(maxLevel() - this->a_level(pCellId)) * (tau - F1B2) / 3.0;

      //--------------------------------------------

      // perform relaxation
      for(MInt j = 0; j < nDist; j++) {
        a_distribution(pCellId, j) = a_oldDistribution(pCellId, j) - m_omega * nonEqDist[j];
      }


      //---------------------------
      // save dissipation values
      if(mode == 2) {
        actualCellLength = this->c_cellLengthAtLevel(this->a_level(pCellId));
        index = floor(
            (F1B2 * m_arraySize[direction]
             + (a_coordinate(pCellId, direction) - F1B2 * (actualCellLength)) / (this->c_cellLengthAtLevel(maxLevel())))
            + 0.1);

        // Calculate characteristic filtered rate of strain
        S = F1B2 * F1BCSsq * m_omega * Q;

        // Calculate molecular- and subgrid dissipation
        if(pCellId < this->grid().noInternalCells()) {
          diss[index] += m_nu * S * S;
          SGSDiss[index] += m_Cs * m_Cs * m_deltaX * m_deltaX * S * S * S * FPOW4(maxLevel() - this->a_level(pCellId));
          for(MInt d = 0; d < nDim; d++) {
            energy += F1B2 * (a_variable(pCellId, d) * a_variable(pCellId, d));
          }
        }
      }
    }
  }

  // write dissipation values
  if(mode == 2 && globalTimeStep % m_solutionInterval == 0) {
    // for (MInt i=0; i<m_arraySize[direction]; i++){
    //   m_log <<"cs["<<i<<"]="<<dynamicCs[i]<<endl;
    // }

    // all domains send their information to domain_0
    if(domainId() != 0) {
      MPI_Send(&(diss[0]), m_arraySize[direction], MPI_DOUBLE, 0, 0, mpiComm(), AT_, "(diss[0])");
      MPI_Send(&(SGSDiss[0]), m_arraySize[direction], MPI_DOUBLE, 0, 0, mpiComm(), AT_, "(SGSDiss[0])");
      MPI_Send(&energy, 1, MPI_DOUBLE, 0, 0, mpiComm(), AT_, "energy");
    }

    if(domainId() == 0) {
      MPI_Status status;

      for(MInt i = 0; i < m_arraySize[direction]; i++) {
        totalDiss[i] = diss[i];
        totalSGSDiss[i] = SGSDiss[i];
      }

      for(MInt j = 1; j < noDomains(); j++) {
        MPI_Recv(&(diss[0]), m_arraySize[direction], MPI_DOUBLE, j, 0, mpiComm(), &status, AT_, "(diss[0])");
        MPI_Recv(&(SGSDiss[0]), m_arraySize[direction], MPI_DOUBLE, j, 0, mpiComm(), &status, AT_, "(SGSDiss[0])");
        MPI_Recv(&energy, 1, MPI_DOUBLE, j, 0, mpiComm(), &status, AT_, "energy");
        for(MInt i = 0; i < m_arraySize[direction]; i++) {
          totalDiss[i] += diss[i];
          totalSGSDiss[i] += SGSDiss[i];
        }
      }

      fileName = "totalDiss_";
      sprintf(buf, "%d", globalTimeStep);
      fileName += buf;
      fileName += ".dat";
      m_log << "writing diss file: " << fileName << endl;
      ofl.open(fileName.c_str(), ios_base::out);
      for(MInt i = 0; i < m_arraySize[direction]; i++)
        ofl << globalTimeStep << " " << i << " " << totalDiss[i] << " " << totalSGSDiss[i] << " " << energy << endl;
      ofl.close();
    }
  }
}

/** \fn LbSolverDxQy::bgki_dynamic_smago_collision_step()
 * \propVal{solverMethod,MAIA_LATTICE_BGKI_DYNAMIC_SMAGO}
 *
 * Collision step of the incompressible LBGK algorithm with SGS modelling
 *
 * only for D3Q19<br>
 * \collstep{LbSolverDxQy::bgki_dynamic_smago_collision_step(), LBCdynsmago_bgki, BGKI dyn. Smagorinsky Collision Step}
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbSolverDxQy<nDim, nDist, SysEqn>::bgki_dynamic_smago_collision_step() {
  TRACE();

  updateMacroscopicVariables();
  for(MInt i = 0; i < this->noInternalCells(); i++) {
    calculateMomentumFlux(i);
  }
  calculateSGSTensors();

  std::vector<MFloat> dynamicCs;
  constexpr MInt direction = 1;
  MInt count = 0;
  averageSGSTensors(direction, count, dynamicCs);

  for(MInt i = 0; i < count; i++) {
    if(dynamicCs[i] < 0.0) dynamicCs[i] = 0.0;
  }


  // For now testing only the D3Q19 algorithm
  MFloat rho = F1;
  MFloat u[nDim] = {F0};
  MFloat tau = F0;

  MFloat Q = F0;
  MInt index = 0;
  MFloat actualCellLength = F0;

  if(m_tanhInit && (globalTimeStep > m_initStartTime)) {
    MFloat scale =
        0.5 * (m_tanhScaleFactor * tanh((5.0 * (MFloat)(globalTimeStep - m_initStartTime) / (MFloat)m_initTime) - 2.5))
        + 0.5;

    m_Re = m_initRe + scale * (m_finalRe - m_initRe);

    // final timestep reached
    if(globalTimeStep == m_initStartTime + m_initTime) {
      m_Re = m_finalRe;
      m_tanhInit = 0;
    }
  }
  m_nu = m_Ma * LBCS / m_Re * m_referenceLength;

  for(MInt i = 0; i < m_currentMaxNoCells; i++) {
    const MInt pCellId = m_activeCellList[i];

    if((globalTimeStep - 1) % IPOW2(maxLevel() - this->a_level(pCellId)) == 0) {
      rho = a_variable(pCellId, PV->RHO);

      for(MInt d = 0; d < nDim; d++) {
        u[d] = a_variable(pCellId, d);
      }

      // Calculation of equilibrium and non-equilibrium distribution
      std::array<MFloat, nDist> eqDist;
      eqDist = getEqDists(rho, u);
      MFloat nonEqDist[nDist];
      for(MInt j = 0; j < nDist; j++) {
        nonEqDist[j] = a_oldDistribution(pCellId, j) - eqDist[j];
      }

      //--------------------------------------------
      // Calculation of overall viscosity
      // according to Hu, Sterling, Chen 1994

      // Calculation of original relaxation time (on current level)
      tau = F1B2 + 3.0 * m_nu * FFPOW2(maxLevel() - this->a_level(pCellId));

      // Calculate square of momentum flux tensor
      Q = 0.0;
      for(MInt k = 0; k < nDim * nDim; k++) {
        Q += m_momentumFlux[pCellId][k] * m_momentumFlux[pCellId][k];
      }
      Q = sqrt(2.0 * Q);

      // Calculation of new relaxation time

      // index = (MInt)((a_coordinate(pCellId, direction) - minCoord) / length * MFloat(count - 1) + 0.5);
      index = floor((F1B2 * count
                     + (a_coordinate(pCellId, 1) - F1B2 * (actualCellLength)) / (this->c_cellLengthAtLevel(maxLevel())))
                    + 0.1);

      tau += F1B2
             * (sqrt(tau * tau
                     + 2.0 * SQRT2 * dynamicCs[index] * m_deltaX * m_deltaX * (F1BCSsq * F1BCSsq) * Q
                           * FPOW2(maxLevel() - this->a_level(pCellId)))
                - tau);

      m_omega = 1.0 / tau;

      // save total nu for bndcnd and interface interpolation (save on finest level)
      a_nu(pCellId) = FPOW2(maxLevel() - this->a_level(pCellId)) * (tau - F1B2) / 3.0;

      // perform relaxation
      for(MInt j = 0; j < nDist; j++) {
        a_distribution(pCellId, j) = a_oldDistribution(pCellId, j) - m_omega * nonEqDist[j];
      }
    }
  }
}


/** \fn LbSolverDxQy::bgki_euler_collision_step()
 * \propVal{solverMethod,MAIA_LATTICE_BGKI_EULER_2D}
 *
 * The collision step for the Euler BGK-Algorithm
 * Why is the equilibrium set if Q < eps???
 * Not tested yet, use with caution!!!!<br>
 * \collstep{LbSolverDxQy::bgki_euler_collision_step(), LBCeuler, Incompressible BGK Euler Collision Step}
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbSolverDxQy<nDim, nDist, SysEqn>::bgki_euler_collision_step() {
  TRACE();

  const MInt pvrho = PV->RHO;

  constexpr MInt fldlen = Ld::dxQyFld();

  MFloat l_rho = 1.0;
  MFloat l_v[nDim] = {0.0};

  MFloat eps = 1e-4;

  // Update the global Reynolds number if we use the tanh-adaption
  if(m_tanhInit && (globalTimeStep > m_initStartTime)) {
    const MFloat scale =
        0.5 * (m_tanhScaleFactor * tanh((5.0 * (MFloat)(globalTimeStep - m_initStartTime) / (MFloat)m_initTime) - 2.5))
        + 0.5;

    m_Re = m_initRe + scale * (m_finalRe - m_initRe);

    // final timestep reached
    if(globalTimeStep == m_initStartTime + m_initTime) {
      m_Re = m_finalRe;
      m_tanhInit = 0;
    }
  }

  // Update the according nu
  m_nu = m_Ma * LBCS / m_Re * m_referenceLength;


  for(MInt i = 0; i < m_currentMaxNoCells; i++) {
    const MInt pCellId = m_activeCellList[i];

    if((globalTimeStep - 1) % IPOW2(maxLevel() - this->a_level(pCellId)) == 0) {
      a_oldVariable(i, PV->RHO) = a_variable(i, PV->RHO);
      a_oldVariable(i, PV->U) = a_variable(i, PV->U);
      a_oldVariable(i, PV->V) = a_variable(i, PV->V);

      // Add forcing term and calculate density rho
      l_rho = F0;
      for(MInt j = 0; j < nDist - 1; j++) {
        l_rho += a_oldDistribution(pCellId, j);
      }
      // Do not forget the last distribution
      l_rho += a_oldDistribution(pCellId, nDist - 1);


      // Calculate macroscopic variables
      l_v[0] = F0;
      l_v[1] = F0;
      if constexpr(nDim == 3) l_v[2] = F0;

      if constexpr(nDim == 2) {
        // WARNING!!! The velocity should be calculated as in 3d
        // However, result is correct
        // calculation of u
        l_v[0] = a_oldDistribution(pCellId, 1) + a_oldDistribution(pCellId, 4) + a_oldDistribution(pCellId, 5)
                 - a_oldDistribution(pCellId, 7) - a_oldDistribution(pCellId, 6) - a_oldDistribution(pCellId, 0);
        // calculation of v
        l_v[1] = a_oldDistribution(pCellId, 7) + a_oldDistribution(pCellId, 3) + a_oldDistribution(pCellId, 4)
                 - a_oldDistribution(pCellId, 6) - a_oldDistribution(pCellId, 2) - a_oldDistribution(pCellId, 5);
      } else {
        for(MInt j = 0; j < fldlen; j++) {
          for(MInt d = 0; d < nDim; d++) {
            l_v[d] += a_oldDistribution(pCellId, Ld::pFld(d, j));
            l_v[d] -= a_oldDistribution(pCellId, Ld::nFld(d, j));
          }
        }
      }

      // Save new macroscopic variables in cell
      for(MInt d = 0; d < nDim; d++) {
        a_variable(pCellId, d) = l_v[d];
      }
      a_variable(pCellId, pvrho) = l_rho;

      calculateMomentumFlux(pCellId);

      // 3. Calculation of the filtered mean momentum flux
      MFloat Q = 0.0;
      for(MInt k = 0; k < nDim * nDim; k++) {
        Q += m_momentumFlux[pCellId][k] * m_momentumFlux[pCellId][k];
      }
      Q = sqrt(2.0 * Q);

      for(MInt d = 0; d < nDim; d++) {
        l_v[d] = a_variable(pCellId, d);
      }
      l_rho = a_variable(pCellId, pvrho);

      const MFloat sqVel = std::inner_product(&l_v[0], &l_v[nDim], &l_v[0], F0);

      std::array<MFloat, nDist> eqDist;
      eqDist = getEqDists(l_rho, sqVel, l_v);

      // Calculation of new distributions
      if(fabs(Q) > eps) {
        const MFloat l_omega = 2.0 / (1.0 + 6.0 * m_nu * FFPOW2(maxLevel() - this->a_level(pCellId)));

        // Calculation of new distributions
        for(MInt j = 0; j < nDist; j++) {
          a_distribution(pCellId, j) =
              a_oldDistribution(pCellId, j) + l_omega * (eqDist[j] - a_oldDistribution(pCellId, j));
        }
      } else {
        // Calculation of new distributions
        for(MInt j = 0; j < nDist; j++) {
          a_distribution(pCellId, j) = eqDist[j];
        }
      }
    }
  }
}

/** \brief  Calculate total energy, dissipation, and subgrid dissipation for Smagorinsky
 *  \note   Only on highest refinement-level !!!
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbSolverDxQy<nDim, nDist, SysEqn>::calculateDissipation() {
  TRACE();

  if(globalTimeStep % 5 != 0) {
    return;
  }

  constexpr MInt nDimSqr = nDim * nDim;
  constexpr MInt direction = 1;

  // reset dissipation and energy values
  ScratchSpace<MFloat> diss(m_arraySize[direction], AT_, "diss");
  ScratchSpace<MFloat> SGSDiss(m_arraySize[direction], AT_, "SGSDiss");
  diss.fill(F0);
  SGSDiss.fill(F0);
  MFloat energy = 0.0;

  m_nu = m_Ma * LBCS / m_Re * m_referenceLength;

  for(MInt i = 0; i < m_currentMaxNoCells; i++) {
    const MInt pCellId = m_activeCellList[i];

    calculateMomentumFlux(pCellId);

    // skip lower levels
    if(this->a_level(pCellId) < maxLevel()) continue;

    // // Calculate gradients
    // for(MInt var = 0; var < nDim; var++) {
    //   for(MInt dir = 0; dir < nDim; dir++) {
    //     L = c;
    //     R = c;
    //     if(a_hasNeighbor(c, 2 * dir) > 0) L = c_neighborId(c, 2 * dir);
    //     if(a_hasNeighbor(c, 2 * dir + 1) > 0) R = c_neighborId(c, 2 * dir + 1);
    //     if(L == c) {
    //       gradV[var + 3 * dir] = (a_variable(R, var) - a_variable(c, var));
    //       continue;
    //     }
    //     if(R == c) {
    //       gradV[var + 3 * dir] = (a_variable(c, var) - a_variable(L, var));
    //       continue;
    //     }
    //     gradV[var + 3 * dir] = (a_variable(R, var) - a_variable(L, var)) / 2.0;
    //   }
    // }

    // Calculate strain rate tensor
    // S_{ij} = \frac{1}{2} \partial_j u_i \partial_i u_j
    m_omega = 2.0 / (1.0 + 6.0 * m_nu * FPOW2(maxLevel() - this->a_level(pCellId)));
    MFloat Sij[nDimSqr];
    for(MInt k = 0; k < nDimSqr; k++) {
      Sij[k] = -F1B2 * F1BCSsq * m_omega * m_momentumFlux[i][k];
    }

    // Calculate characteristic filtered rate of strain
    const MFloat Ssqr = 2.0 * std::inner_product(&Sij[0], &Sij[nDimSqr], &Sij[0], .0);
    const MFloat S = sqrt(Ssqr);

    // Calculate molecular- and subgrid dissipation
    const MInt actualCellLength = this->c_cellLengthAtLevel(this->a_level(pCellId));
    const MInt index = floor(
        (F1B2 * m_arraySize[direction]
         + (a_coordinate(pCellId, direction) - F1B2 * (actualCellLength)) / (this->c_cellLengthAtLevel(maxLevel())))
        + 0.1);
    diss[index] += m_nu * Ssqr;
    SGSDiss[index] += m_Cs * m_Cs * m_deltaX * m_deltaX * S * S * S * FPOW4(maxLevel() - this->a_level(pCellId));

    // energy
    MFloat tmp = 0.0;
    for(MInt k = 0; k < nDim; k++) {
      tmp += a_variable(pCellId, k) * a_variable(pCellId, k);
    }
    energy += F1B2 * tmp;
  }

  // Sum-up on root domain
  constexpr MInt rootId = 0;
  MPI_Reduce(MPI_IN_PLACE, &diss[0], m_arraySize[direction], MPI_DOUBLE, MPI_SUM, rootId, mpiComm(), AT_,
             "MPI_IN_PLACE", "diss");
  MPI_Reduce(MPI_IN_PLACE, &SGSDiss[0], m_arraySize[direction], MPI_DOUBLE, MPI_SUM, rootId, mpiComm(), AT_,
             "MPI_IN_PLACE", "SGSDiss");
  MPI_Reduce(MPI_IN_PLACE, &energy, m_arraySize[direction], MPI_DOUBLE, MPI_SUM, rootId, mpiComm(), AT_, "MPI_IN_PLACE",
             "energy");

  // Write out data from root process
  if(domainId() == rootId) {
    ofstream ofl;
    MString fileName = "totalDiss_";
    MChar buf[10];
    sprintf(buf, "%d", globalTimeStep);
    fileName += buf;
    fileName += ".dat";
    m_log << "writing diss file: " << fileName << endl;
    ofl.open(fileName.c_str(), ios_base::out);
    for(MInt i = 0; i < m_arraySize[direction]; i++)
      ofl << globalTimeStep << " " << i << " " << diss[i] << " " << SGSDiss[i] << " " << energy << endl;
    ofl.close();
  }
}

/** \brief  Collision step for the MRT-Algorithm
 *  \author Miro Gondrum (Refactored)
 *  \date   03.11.2021
 *
 *  This function provides a collision step based on the MRT-Algorithm. It
 *  provides a version based on standard and optimized parameters
 *  Ref. :
 *    - (D3Q19) : d'Humieres, 2002, https://doi.org/10.1098/rsta.2001.0955
 *  \note   Currently only available for D2Q9 and D3Q19
 */
template <MInt nDim, MInt nDist, class SysEqn>
template <MBool optimized, MBool useSmagorinsky>
void LbSolverDxQy<nDim, nDist, SysEqn>::mrt_collision_step_base() {
  TERMM(1, "MRT collision step only available for D2Q9 and D3Q19 !");
}

template <>
template <MBool optimized, MBool useSmagorinsky>
void LbSolverDxQy<2, 9, maia::lb::LbSysEqnIncompressible<2, 9>>::mrt_collision_step_base() {
  TRACE();
  constexpr MInt nDist = 9;
  // constexpr MInt nDim = 2;

  constexpr MFloat c1 = -2.0, alpha2 = -8.0, alpha3 = 4.0, gamma1 = F2B3, gamma3 = F2B3, gamma2 = 18.0, gamma4 = -18.0;

  for(MInt i = 0; i < m_currentMaxNoCells; i++) {
    const MInt pCellId = m_activeCellList[i];
    const MInt lvlDiff = maxLevel() - this->a_level(pCellId);
    if((globalTimeStep - 1) % IPOW2(lvlDiff) == 0) {
      std::array<MFloat, nDist> d{};
      for(MInt j = 0; j < nDist; j++) {
        d[j] = a_oldDistribution(pCellId, j);
      }

      std::array<MFloat, nDist> m{};
      m[0] = d[8] + d[1] + d[3] + d[0] + d[2] + d[4] + d[7] + d[6] + d[5];
      m[1] = -4 * d[8] - d[1] - d[3] - d[0] - d[2] + 2 * (d[4] + d[7] + d[6] + d[5]);
      m[2] = 4 * d[8] + 2 * (-d[1] - d[3] - d[0] - d[2]) + d[4] + d[7] + d[6] + d[5];
      m[3] = d[1] - d[0] + d[4] - d[7] - d[6] + d[5];
      m[4] = -2 * (d[1] - d[0]) + d[4] - d[7] - d[6] + d[5];
      m[5] = d[3] - d[2] + d[4] + d[7] - d[6] - d[5];
      m[6] = -2 * (d[3] - d[2]) + d[4] + d[7] - d[6] - d[5];
      m[7] = d[1] - d[3] + d[0] - d[2];
      m[8] = d[4] - d[7] + d[6] - d[5];

      // Relaxation in moment space
      // 1. Set relaxation parameters
      std::array<MFloat, nDist> omega{};
      omega[0] = 0.0;
      omega[1] = 1.63;
      omega[2] = 1.14;
      omega[3] = 0.0;
      omega[5] = 0.0;
      omega[6] = 1.92;
      omega[7] = (2.0 - c1) / (12.0 * a_nu(pCellId) * FFPOW2(lvlDiff) + 1.0 - c1 / 2.0);
      omega[8] = 1.0 / ((2.0 / omega[7] - 1.0) * ((c1 + 4.0) / (2.0 - c1)) + 0.5);
      omega[4] = 3.0 * (2.0 - omega[7]) / (3.0 - omega[7]);

      // 2. Calculate equilibrium Moments
      std::array<MFloat, nDist> EqM{};
      EqM[0] = m[0];
      EqM[3] = m[3];
      EqM[5] = m[5];
      EqM[1] = 0.25 * alpha2 * m[0] + F1B6 * gamma2 * (m[3] * m[3] + m[5] * m[5]);
      EqM[2] = 0.25 * alpha3 * m[0] + F1B6 * gamma4 * (m[3] * m[3] + m[5] * m[5]);
      EqM[4] = 0.5 * c1 * m[3];
      EqM[6] = 0.5 * c1 * m[5];
      EqM[7] = F3B2 * gamma1 * (m[3] * m[3] - m[5] * m[5]);
      EqM[8] = F3B2 * gamma3 * (m[3] * m[5]);

      // 3. Relax Moments
      for(MInt j = 0; j < nDist; j++) {
        m[j] = m[j] - omega[j] * (m[j] - EqM[j]);
      }

      a_distribution(pCellId, 8) = F1B9 * m[0] - F1B9 * m[1] + F1B9 * m[2];
      a_distribution(pCellId, 1) = F1B9 * m[0] - F1B36 * m[1] - F1B18 * m[2] + F1B6 * m[3] - F1B6 * m[4] + F1B4 * m[7];
      a_distribution(pCellId, 3) = F1B9 * m[0] - F1B36 * m[1] - F1B18 * m[2] + F1B6 * m[5] - F1B6 * m[6] - F1B4 * m[7];
      a_distribution(pCellId, 0) = F1B9 * m[0] - F1B36 * m[1] - F1B18 * m[2] - F1B6 * m[3] + F1B6 * m[4] + F1B4 * m[7];
      a_distribution(pCellId, 2) = F1B9 * m[0] - F1B36 * m[1] - F1B18 * m[2] - F1B6 * m[5] + F1B6 * m[6] - F1B4 * m[7];
      a_distribution(pCellId, 4) = F1B9 * m[0] + F1B18 * m[1] + F1B36 * m[2] + F1B6 * m[3] + F1B12 * m[4] + F1B6 * m[5]
                                   + F1B12 * m[6] + F1B4 * m[8];
      a_distribution(pCellId, 7) = F1B9 * m[0] + F1B18 * m[1] + F1B36 * m[2] - F1B6 * m[3] - F1B12 * m[4] + F1B6 * m[5]
                                   + F1B12 * m[6] - F1B4 * m[8];
      a_distribution(pCellId, 6) = F1B9 * m[0] + F1B18 * m[1] + F1B36 * m[2] - F1B6 * m[3] - F1B12 * m[4] - F1B6 * m[5]
                                   - F1B12 * m[6] + F1B4 * m[8];
      a_distribution(pCellId, 5) = F1B9 * m[0] + F1B18 * m[1] + F1B36 * m[2] + F1B6 * m[3] + F1B12 * m[4] - F1B6 * m[5]
                                   - F1B12 * m[6] - F1B4 * m[8];
    }
  }
}

template <>
template <MBool optimized, MBool useSmagorinsky>
void LbSolverDxQy<3, 19, maia::lb::LbSysEqnIncompressible<3, 19>>::mrt_collision_step_base() {
  TRACE();
  constexpr MInt nDist = 19;
  constexpr MInt nDim = 3;
  constexpr MInt nDimSqr = nDim * nDim;

  const MFloat rho_offset = (m_densityFluctuations) ? 1.0 : 0.0;

  // Only relevant for Smagorinsky model
  [[maybe_unused]] const MFloat tmpWidth = m_referenceLength * m_smallestCellLength;
  [[maybe_unused]] constexpr MFloat aPlus = 26.0;

  for(MInt id = 0; id < m_currentMaxNoCells; id++) {
    const MInt pCellId = m_activeCellList[id];

    if((globalTimeStep - 1) % IPOW2(maxLevel() - this->a_level(pCellId)) == 0) {
      m_nu = m_Ma * LBCS / m_Re * m_referenceLength;

      MFloat d[nDist];
      for(MInt j = 0; j < nDist; j++) {
        d[j] = a_oldDistribution(pCellId, j);
      }

      // Set old variables (for residual calculation)
      for(MInt j = 0; j < nDim + 1; j++) {
        a_oldVariable(pCellId, j) = a_variable(pCellId, j);
      }

      MFloat m[nDist];
      m[0] = d[18] + d[1] + d[0] + d[3] + d[2] + d[5] + d[4] + d[9] + d[7] + d[8] + d[6] + d[13] + d[11] + d[12] + d[10]
             + d[17] + d[15] + d[16] + d[14];
      m[1] = -30 * d[18] - 11 * (d[1] + d[0] + d[3] + d[2] + d[5] + d[4])
             + 8 * (d[9] + d[7] + d[8] + d[6] + d[13] + d[11] + d[12] + d[10] + d[17] + d[15] + d[16] + d[14]);
      m[2] = 12 * d[18] - 4 * (d[1] + d[0] + d[3] + d[2] + d[5] + d[4]) + d[9] + d[7] + d[8] + d[6] + d[13] + d[11]
             + d[12] + d[10] + d[17] + d[15] + d[16] + d[14];
      m[3] = d[1] - d[0] + d[9] - d[7] + d[8] - d[6] + d[13] - d[11] + d[12] - d[10];
      m[4] = -4 * (d[1] - d[0]) + d[9] - d[7] + d[8] - d[6] + d[13] - d[11] + d[12] - d[10];
      m[5] = d[3] - d[2] + d[9] + d[7] - d[8] - d[6] + d[17] - d[15] + d[16] - d[14];
      m[6] = -4 * (d[3] - d[2]) + d[9] + d[7] - d[8] - d[6] + d[17] - d[15] + d[16] - d[14];
      m[7] = d[5] - d[4] + d[13] + d[11] - d[12] - d[10] + d[17] + d[15] - d[16] - d[14];
      m[8] = -4 * (d[5] - d[4]) + d[13] + d[11] - d[12] - d[10] + d[17] + d[15] - d[16] - d[14];
      m[9] = 2 * (d[1] + d[0] - d[15] - d[16] - d[14]) - d[3] - d[2] - d[5] - d[4] + d[9] + d[7] + d[8] + d[6] + d[13]
             + d[11] + d[12] + d[10] - 2 * d[17];
      m[10] = -4 * (d[1] + d[0]) + 2 * (d[3] + d[2] + d[5] + d[4] - d[17] - d[15] - d[16] - d[14]) + d[9] + d[7] + d[8]
              + d[6] + d[13] + d[11] + d[12] + d[10];
      m[11] = d[3] + d[2] - d[5] - d[4] + d[9] + d[7] + d[8] + d[6] - d[13] - d[11] - d[12] - d[10];
      m[12] = -2 * (d[3] + d[2] - d[5] - d[4]) + d[9] + d[7] + d[8] + d[6] - d[13] - d[11] - d[12] - d[10];
      m[13] = d[9] - d[7] - d[8] + d[6];
      m[14] = d[17] - d[15] - d[16] + d[14];
      m[15] = d[13] - d[11] - d[12] + d[10];
      m[16] = d[9] - d[7] + d[8] - d[6] - d[13] + d[11] - d[12] + d[10];
      m[17] = -d[9] - d[7] + d[8] + d[6] + d[17] - d[15] + d[16] - d[14];
      m[18] = d[13] + d[11] - d[12] - d[10] - d[17] - d[15] + d[16] + d[14];

      // Relaxation in moment space
      // 1. Set relaxation parameters
      MFloat P[nDist];
      const MFloat tmp = (optimized) ? 0.0 : 1.0;
      P[0] = tmp;
      P[1] = 1.19;
      P[2] = 1.4;
      P[3] = tmp;
      P[4] = 1.2;
      P[5] = tmp;
      P[6] = P[4];
      P[7] = tmp;
      P[8] = P[4];
      P[9] = 2.0 / (FFPOW2(maxLevel() - this->a_level(pCellId)) * 6.0 * m_nu + 1.0);
      P[10] = P[2];
      P[11] = P[9];
      P[12] = P[2];
      P[13] = P[9];
      P[14] = P[9];
      P[15] = P[9];
      P[16] = 1.98;
      P[17] = P[16];
      P[18] = P[16];

      // 2. Calculate equilibrium Moments
      MFloat EqM[nDist];
      EqM[0] = m[0]; // rho
      EqM[3] = m[3]; // j_x
      EqM[5] = m[5]; // j_y
      EqM[7] = m[7]; // j_z

      //---------------------------------------------------------------------------------
      // d'Humieres et al., 2002 (standard/optimized parameters)
      // compressible form
      EqM[1] = -11.0 * m[0] + 19.0 * (m[3] * m[3] + m[5] * m[5] + m[7] * m[7]) / (m[0] + rho_offset); // e

      if constexpr(optimized) {
        EqM[2] = -F475B63 * (m[3] * m[3] + m[5] * m[5] + m[7] * m[7]) / (m[0] + rho_offset); // epsilon
      } else {
        EqM[2] = 3.0 * m[0] - 5.5 * (m[3] * m[3] + m[5] * m[5] + m[7] * m[7]) / (m[0] + rho_offset); // epsilon
      }

      EqM[4] = -F2B3 * m[3]; // q_x
      EqM[6] = -F2B3 * m[5]; // q_y
      EqM[8] = -F2B3 * m[7]; // q_z

      EqM[9] = (2 * m[3] * m[3] - (m[5] * m[5] + m[7] * m[7])) / (m[0] + rho_offset); // 3*p_xx

      if constexpr(optimized) {
        EqM[10] = 0.0;
      } else {
        EqM[10] = -0.5 * EqM[9]; // 3*pi_xx
      }
      EqM[11] = (m[5] * m[5] - m[7] * m[7]) / (m[0] + rho_offset); // p_ww
      if constexpr(optimized) {
        EqM[12] = 0.0;
      } else {
        EqM[12] = -0.5 * EqM[11]; // pi_ww
      }
      EqM[13] = m[3] * m[5] / (m[0] + rho_offset); // p_xy
      EqM[14] = m[5] * m[7] / (m[0] + rho_offset); // p_yz
      EqM[15] = m[3] * m[7] / (m[0] + rho_offset); // p_xz
      EqM[16] = 0;
      EqM[17] = 0;
      EqM[18] = 0;

      // // incompressible form assuming <rho> = 1
      // EqM[1]  = -11.0 * m[0] + 19.0 * ( m[3] * m[3] + m[5] * m[5] + m[7] * m[7] ); // e

      // EqM[2]  = 3.0 * m[0] - 5.5 * ( m[3] * m[3] + m[5] * m[5] + m[7] * m[7] ); // epsilon

      // EqM[4]  = -F2B3 * m[3] ; // q_x
      // EqM[6]  = -F2B3 * m[5] ; // q_y
      // EqM[8]  = -F2B3 * m[7] ; // q_z

      // EqM[9]  = ( 2*m[3]*m[3] - (m[5] * m[5] + m[7] * m[7]) ); // 3*p_xx

      // EqM[10] = -0.5 * EqM[9] ; // 3*pi_xx
      // EqM[11] = ( m[5] * m[5] - m[7] * m[7] ); // p_ww
      // EqM[12] = -0.5 * EqM[11] ; // pi_ww
      // EqM[13] = m[3] * m[5]; // p_xy
      // EqM[14] = m[5] * m[7]; // p_yz
      // EqM[15] = m[3] * m[7]; // p_xz
      //       EqM[16] = 0 ;
      //       EqM[17] = 0 ;
      //       EqM[18] = 0 ;

      //---------------------------------------------------------------------------------

      // 2b. Determine macroscopic variables and calculate residual
      a_variable(pCellId, PV->RHO) = EqM[0];
      for(MInt i = 0; i < nDim; i++) {
        a_variable(pCellId, PV->VV[i]) = EqM[3 + i * 2] / EqM[0];
      }

      // 2c. Smagorinsky model
      if constexpr(useSmagorinsky) {
        // Calculation of equilibrium and non-equilibrium distribution
        MFloat nonEq[nDist];
#ifdef WAR_NVHPC_PSTL
        MFloat dist[nDist] = {F0};
        MFloat u[nDim] = {F0};
        for(MInt dir = 0; dir < nDist; dir++) {
          dist[dir] = a_oldDistribution(pCellId, dir);
        }
        for(MInt dir = 0; dir < nDim; dir++) {
          u[dir] = a_variable(pCellId, dir);
        }
        MFloat l_rho = a_variable(pCellId, PV->RHO);
        lbfunc::calcNonEqDists<nDim, nDist>(l_rho, u, dist, &nonEq[0]);
#else
        lbfunc::calcNonEqDists<nDim, nDist>(a_variable(pCellId, PV->RHO), &a_variable(pCellId, PV->U),
                                            &a_oldDistribution(pCellId, 0), &nonEq[0]);
#endif

        // Smagorinsky model
        const MFloat yPlus = m_ReTau * (1.0 - fabs(a_coordinate(pCellId, 1)) / tmpWidth);
        const MFloat lambda = m_Cs * m_deltaX * (F1 - exp(-yPlus / aPlus));

        // Calculation of momentum flux tensor from non-equilibrium parts
        MFloat c[nDimSqr];
        c[0] = (nonEq[0] + nonEq[1] + nonEq[6] + nonEq[7] + nonEq[8] + nonEq[9] + nonEq[10] + nonEq[11] + nonEq[12]
                + nonEq[13]);
        c[1] = (nonEq[6] - nonEq[7] - nonEq[8] + nonEq[9]);
        c[2] = (nonEq[10] - nonEq[11] - nonEq[12] + nonEq[13]);
        c[3] = c[1];
        c[4] = (nonEq[2] + nonEq[3] + nonEq[6] + nonEq[7] + nonEq[8] + nonEq[9] + nonEq[14] + nonEq[15] + nonEq[16]
                + nonEq[17]);
        c[5] = (nonEq[14] - nonEq[15] - nonEq[16] + nonEq[17]);
        c[6] = c[2];
        c[7] = c[5];
        c[8] = (nonEq[4] + nonEq[5] + nonEq[10] + nonEq[11] + nonEq[12] + nonEq[13] + nonEq[14] + nonEq[15] + nonEq[16]
                + nonEq[17]);

        // Calculation of the filtered mean momentum flux
        const MFloat Q = sqrt(2.0 * std::inner_product(&c[0], &c[nDimSqr], &c[0], .0));

        // Calculation of new relaxation time
        // ..original relaxation time (on current level)
        MFloat tau = F1B2 + 3.0 * m_nu * FFPOW2(maxLevel() - this->a_level(pCellId));
        // ..new
        tau += F1B2
               * (sqrt(tau * tau
                       + 2.0 * SQRT2 * lambda * lambda * (F1BCSsq * F1BCSsq) * Q
                             * FPOW2(maxLevel() - this->a_level(pCellId)))
                  - tau);

        m_omega = 1.0 / tau;

        // save total nu for bndcnd and interface interpolation (save on finest level)
        a_nu(pCellId) = FPOW2(maxLevel() - this->a_level(pCellId)) * (tau - F1B2) / 3.0;

        // set new relaxation time
        P[9] = m_omega;
      }

      // 3. Relax Moments
      for(MInt i = 0; i < nDist; i++) {
        m[i] = m[i] - P[i] * (m[i] - EqM[i]);
      }

      const MFloat ra = F1B19 * m[0];
      const MFloat rb = F11B2394 * m[1];
      const MFloat rc = F1B63 * m[2];
      const MFloat rd = F4B1197 * m[1];
      const MFloat re = F1B252 * m[2];
      const MFloat rf = F1B12 * m[11];
      const MFloat rg = F1B24 * m[12];
      const MFloat rh = 0.25 * m[13];
      const MFloat ri = 0.25 * m[15];
      const MFloat rj = F1B36 * m[9];
      const MFloat rk = 2.0 * rj;
      const MFloat rl = F1B72 * m[10];
      const MFloat rm = 2.0 * rl;
      const MFloat rn = 0.1 * (m[3] - m[4]);
      const MFloat ro = 0.1 * (m[5] - m[6]);
      const MFloat rp = 0.1 * (m[7] - m[8]);
      const MFloat rq = F1B18 * (m[9] - m[10]);
      const MFloat rr = F1B12 * (m[11] - m[12]);
      const MFloat rs = 0.25e-1 * (m[4] + m[6]);
      const MFloat rt = 0.25e-1 * (m[4] - m[6]);
      const MFloat ru = 0.25e-1 * (m[4] - m[8]);
      const MFloat rv = 0.25e-1 * (m[4] + m[8]);
      const MFloat rw = 0.25e-1 * (m[6] + m[8]);
      const MFloat rx = 0.25e-1 * (m[6] - m[8]);
      const MFloat ry = 0.125 * (m[16] - m[17]);
      const MFloat rz = 0.125 * (m[16] + m[17]);
      const MFloat r0 = 0.125 * (m[16] - m[18]);
      const MFloat r1 = 0.125 * (m[16] + m[18]);
      const MFloat r2 = 0.125 * (m[17] + m[18]);
      const MFloat r3 = 0.125 * (m[17] - m[18]);
      const MFloat r4 = 0.1 * (m[5] + m[7]);
      const MFloat r5 = 0.1 * (m[5] - m[7]);
      const MFloat r6 = 0.1 * (m[3] + m[5]);
      const MFloat r7 = 0.1 * (m[3] - m[5]);
      const MFloat r8 = 0.5 * rq;
      const MFloat r9 = 0.1 * (m[3] + m[7]);
      const MFloat r10 = 0.1 * (m[3] - m[7]);
      const MFloat r11 = 0.25 * m[14];

      //       ra =  0.5263157895e-1*m[0];
      //       rb =  0.4594820384e-2*m[1];
      //       rc =  0.1587301587e-1*m[2];
      //       rd =  0.3341687552e-2*m[1];
      //       re =  0.3968253968e-2*m[2];
      //       rf =  0.8333333333e-1*m[11];
      //       rg =  0.4166666667e-1*m[12];
      //       rh =  0.250*m[13];
      //       ri =  0.250*m[15];
      //       rj =  0.2777777778e-1*m[9];
      //       rk =  2*rj;
      //       rl =  0.1388888889e-1*m[10];
      //       rm =  2*rl;
      //       rn =  0.10*(m[3]-m[4]);
      //       ro =  0.10*(m[5]-m[6]);
      //       rp =  0.10*(m[7]-m[8]);
      //       rq =  0.5555555556e-1*(m[9]-m[10]);
      //       rr =  0.8333333333e-1*(m[11]-m[12]);
      //       rs =  0.250e-1*(m[4]+m[6]);
      //       rt =  0.250e-1*(m[4]-m[6]);
      //       ru =  0.250e-1*(m[4]-m[8]);
      //       rv =  0.250e-1*(m[4]+m[8]);
      //       rw =  0.250e-1*(m[6]+m[8]);
      //       rx =  0.250e-1*(m[6]-m[8]);
      //       ry =  0.1250*(m[16]-m[17]);
      //       rz =  0.1250*(m[16]+m[17]);
      //       r0 =  0.1250*(m[16]-m[18]);
      //       r1 =  0.1250*(m[16]+m[18]);
      //       r2 =  0.1250*(m[17]+m[18]);
      //       r3 =  0.1250*(m[17]-m[18]);
      //       r4 =  0.10*(m[5]+m[7]);
      //       r5 =  0.10*(m[5]-m[7]);
      //       r6 =  0.10*(m[3]+m[5]);
      //       r7 =  0.10*(m[3]-m[5]);
      //       r8 =  0.5*rq;
      //       r9 =  0.10*(m[3]+m[7]);
      //       r10=  0.10*(m[3]-m[7]);
      //       r11=  0.250*m[14];

      a_distribution(pCellId, 18) = ra - F5B399 * m[1] + F1B21 * m[2];
      a_distribution(pCellId, 1) = ra - rb - rc + rn + rq;
      a_distribution(pCellId, 0) = ra - rb - rc - rn + rq;
      a_distribution(pCellId, 3) = ra - rb - rc + ro - r8 + rr;
      a_distribution(pCellId, 2) = ra - rb - rc - ro - r8 + rr;
      a_distribution(pCellId, 5) = ra - rb - rc + rp - r8 - rr;
      a_distribution(pCellId, 4) = ra - rb - rc - rp - r8 - rr;
      a_distribution(pCellId, 9) = ra + rd + re + r6 + rs + rj + rl + rf + rg + rh + ry;
      a_distribution(pCellId, 7) = ra + rd + re - r7 - rt + rj + rl + rf + rg - rh - rz;
      a_distribution(pCellId, 8) = ra + rd + re + r7 + rt + rj + rl + rf + rg - rh + rz;
      a_distribution(pCellId, 6) = ra + rd + re - r6 - rs + rj + rl + rf + rg + rh - ry;
      a_distribution(pCellId, 13) = ra + rd + re + r9 + rv + rj + rl - rf - rg + ri - r0;
      a_distribution(pCellId, 11) = ra + rd + re - r10 - ru + rj + rl - rf - rg - ri + r1;
      a_distribution(pCellId, 12) = ra + rd + re + r10 + ru + rj + rl - rf - rg - ri - r1;
      a_distribution(pCellId, 10) = ra + rd + re - r9 - rv + rj + rl - rf - rg + ri + r0;
      a_distribution(pCellId, 17) = ra + rd + re + r4 + rw - rk - rm + r11 + r3;
      a_distribution(pCellId, 15) = ra + rd + re - r5 - rx - rk - rm - r11 - r2;
      a_distribution(pCellId, 16) = ra + rd + re + r5 + rx - rk - rm - r11 + r2;
      a_distribution(pCellId, 14) = ra + rd + re - r4 - rw - rk - rm + r11 - r3;
    }
  }

  if(m_calculateDissipation && globalTimeStep % m_solutionInterval == 0) {
    calculateDissipation();
  }
}

/** \fn LbSolverDxQy::mrt_collision_step()
 *  \author Miro Gondrum (Refactored)
 *  \date   03.11.2021
 * \propVal{solverMethod,MAIA_LATTICE_MRT}
 *
 * Collision step for the MRT-Algorithm with standard parameters
 * This function wraps mrt_collision_step_base to access MRT-Algorithm with
 * standard parameters.<br>
 * \collstep{LbSolverDxQy::mrt_collision_step(), LBCmrt_std, MRT Collision Step (standard)}
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbSolverDxQy<nDim, nDist, SysEqn>::mrt_collision_step() {
  mrt_collision_step_base<false, false>();
}

/** \fn LbSolverDxQy::mrt2_collision_step()
 *  \author Miro Gondrum (Refactored)
 *  \date   03.11.2021
 * \propVal{solverMethod,MAIA_LATTICE_MRT2}
 *
 * Collision step for the MRT-Algorithm with optimized parameters
 * This function wraps mrt_collision_step_base to access MRT-Algorithm with
 * optimized parameters.<br>
 * \collstep{LbSolverDxQy::mrt2_collision_step(), LBCmrt_opt, MRT Collision Step (optimized)}
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbSolverDxQy<nDim, nDist, SysEqn>::mrt2_collision_step() {
  mrt_collision_step_base<true, false>();
}

/** \fn LbSolverDxQy::mrt_smagorinsky_collision_step()
 *  \author Miro Gondrum (Refactored)
 *  \date   03.11.2021
 * \propVal{solverMethod,MAIA_LATTICE_MRT_SMAGORINSKY}
 *
 * Collision step for the MRT-Algorithm + Smagorinsky turbulence model
 * This function wraps mrt_collision_step_base to access MRT-Algorithm with
 * standard parameters. It features a  Smagorinsky turbulence model for
 * calculating the turbulent viscosity with a van-Driest damping for turbulent
 * channel.<br>
 * \collstep{LbSolverDxQy::mrt_smagorinsky_collision_step(), LBCmrt_smago, MRT Smagorinsky Collision Step}
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbSolverDxQy<nDim, nDist, SysEqn>::mrt_smagorinsky_collision_step() {
  if constexpr(nDim == 3 || nDist == 19) {
    mrt_collision_step_base<false, true>();
  } else {
    TERMM(1, "MRT_SMAGORINSKY collision step only available for D3Q19!");
  }
}

template <MInt nDim, MInt nDist, class SysEqn>
template <MBool useSmagorinsky>
void LbSolverDxQy<nDim, nDist, SysEqn>::rbgk_collision_step_base() {
  TRACE();

  if(m_tanhInit && (globalTimeStep > m_initStartTime)) {
    MFloat scale =
        0.5 * (m_tanhScaleFactor * tanh((5.0 * (MFloat)(globalTimeStep - m_initStartTime) / (MFloat)m_initTime) - 2.5))
        + 0.5;

    m_Re = m_initRe + scale * (m_finalRe - m_initRe);

    // final timestep reached
    if(globalTimeStep == m_initStartTime + m_initTime) {
      m_Re = m_finalRe;
      m_tanhInit = 0;
    }
  }
  m_nu = m_Ma * LBCS / m_Re * m_referenceLength;

  for(MInt i = 0; i < m_currentMaxNoCells; i++) {
    const MInt pCellId = m_activeCellList[i];

    if((globalTimeStep - 1) % IPOW2(maxLevel() - this->a_level(pCellId)) == 0) {
      // save nu for bndcnd and interface interpolation
      a_nu(pCellId) = m_nu;

      m_omega = 2.0 / (1.0 + 6.0 * m_nu * FFPOW2(maxLevel() - this->a_level(pCellId)));
      swap_variables(pCellId);

      // Calculate macroscopic variables
      a_variable(pCellId, PV->RHO) = 0.0;
      for(MInt j = 0; j < nDist; j++) {
        a_variable(pCellId, PV->RHO) += a_oldDistribution(pCellId, j);
      }

      const MFloat rho = a_variable(pCellId, PV->RHO);

      MFloat u[nDim];
      for(MInt j = 0; j < nDim; j++) {
        u[j] = F0;
      }

      if constexpr(nDim == 2) {
        // warning!!! the velocity should be calculated as in 3d
        // however, result is correct
        // calculation of u
        u[0] = a_oldDistribution(pCellId, 1) + a_oldDistribution(pCellId, 4) + a_oldDistribution(pCellId, 5)
               - a_oldDistribution(pCellId, 7) - a_oldDistribution(pCellId, 6) - a_oldDistribution(pCellId, 0);
        // calculation of v
        u[1] = a_oldDistribution(pCellId, 7) + a_oldDistribution(pCellId, 3) + a_oldDistribution(pCellId, 4)
               - a_oldDistribution(pCellId, 6) - a_oldDistribution(pCellId, 2) - a_oldDistribution(pCellId, 5);
      } else {
        for(MInt j = 0; j < Ld::dxQyFld(); j++) {
          for(MInt d = 0; d < nDim; d++) {
            u[d] += a_oldDistribution(pCellId, Ld::pFld(d, j));
            u[d] -= a_oldDistribution(pCellId, Ld::nFld(d, j));
          }
        }
      }

      for(MInt j = 0; j < nDim; j++) {
        a_variable(pCellId, j) = u[j];
      }

      // Calculation of equilibrium and non-equilibrium distribution
      std::array<MFloat, nDist> eqDist;
      eqDist = getEqDists(rho, u);

      // Calculation of momentum flux tensor for non-equilibrium parts
      //--------------------------------------------
      // \Pi_{\alpha,\beta} = \sum_i f_i \cdot c_{i,\alpha}c_{i,\beta}

      calculateMomentumFlux(pCellId);

      if constexpr(useSmagorinsky) {
        // Calculation of overall viscosity
        // according to Hu, Sterling, Chen 1994

        // 1. Calculation of original relaxation time (on current level)
        MFloat tau = F1B2 + 3.0 * m_nu * FFPOW2(maxLevel() - this->a_level(pCellId));

        // 2. Calculation of the filtered mean momentum flux
        MFloat Q = F0;
        for(MInt l = 0; l < nDim * nDim; l++) {
          Q += m_momentumFlux[pCellId][l] * m_momentumFlux[pCellId][l];
        }
        Q = sqrt(2 * Q);

        // 3. Calculation of new relaxation time
        tau += F1B2
               * (sqrt(tau * tau
                       + 2.0 * SQRT2 * m_Cs * m_Cs * m_deltaX * m_deltaX * (F1BCSsq * F1BCSsq) * Q
                             * FPOW2(maxLevel() - this->a_level(pCellId)))
                  - tau);

        m_omega = 1.0 / tau;

        // save total nu for bndcnd and interface interpolation (save on finest level)
        a_nu(pCellId) = FPOW2(maxLevel() - this->a_level(pCellId)) * (tau - F1B2) / 3.0;
      }

      // Calculation of new distributions for directions with at least one component
      MFloat trace = F0;

      if constexpr(nDim == 2) {
        trace = m_momentumFlux[pCellId][0] + m_momentumFlux[pCellId][3];
      } else {
        trace = m_momentumFlux[pCellId][0] + m_momentumFlux[pCellId][4] + m_momentumFlux[pCellId][8];
      }

      for(MInt j = 0; j < nDist - 1; j++) {
        const MInt t = Ld::tp(Ld::distType(j));

        a_distribution(pCellId, j) = eqDist[j] - (1.0 - m_omega) * t * F1BCSsq * F1B2 * trace;

        for(MInt k = 0; k < nDim; k++) {
          for(MInt l = 0; l < nDim; l++) {
            a_distribution(pCellId, j) += (1.0 - m_omega) * t * F1BCSsq * F1BCSsq * F1B2 * (Ld::idFld(j, k) - 1)
                                          * (Ld::idFld(j, l) - 1) * m_momentumFlux[pCellId][l + nDim * k];
          }
        }
      }
      // Calculation of new distribution for rest particle distribution (center)
      a_distribution(pCellId, Ld::lastId()) =
          eqDist[Ld::lastId()] - (1.0 - m_omega) * Ld::tp(0) * F1BCSsq * F1B2 * trace;
    }
  }

  if(m_calculateDissipation && globalTimeStep % m_solutionInterval == 0) {
    calculateDissipation();
  }
}

/** \fn LbSolverDxQy::rbgk_collision_step()
 * \propVal{solverMethod,MAIA_LATTICE_RBGK}
 *
 * Collision step for the regularized LBGK-Algorithm<br>
 * \collstep{LbSolverDxQy::rbgk_collision_step(), LBCrbgk_std, Regularized BGK Collision Step}
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbSolverDxQy<nDim, nDist, SysEqn>::rbgk_collision_step() {
  rbgk_collision_step_base<false>();
}

/** \fn LbSolverDxQy::rbgk_smagorinsky_collision_step()
 * \propVal{solverMethod,MAIA_LATTICE_RBGK_SMAGORINSKY}
 *
 * Collision step for the regularized LBGK-Algorithm + Smagorinsky turbulence model<br>
 * \collstep{LbSolverDxQy::rbgk_smagorinsky_collision_step(), LBCrbgk_smago, Regularized BGK Smagorinsky Collision Step}
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbSolverDxQy<nDim, nDist, SysEqn>::rbgk_smagorinsky_collision_step() {
  rbgk_collision_step_base<true>();
}

/** \fn LbSolverDxQy::rbgk_dynamic_smago_collision_step()
 * \propVal{solverMethod,MAIA_LATTICE_RBGK_DYNAMIC_SMAGO}
 *
 * Collision step for the regularized LBGK-Algorithm + dynamic Smagorinsky turbulence model<br>
 * \collstep{LbSolverDxQy::rbgk_dynamic_smago_collision_step(), LBCrbgk_dynsmago, Reg. BGK dyn. Smago. Collision Step}
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbSolverDxQy<nDim, nDist, SysEqn>::rbgk_dynamic_smago_collision_step() {
  TRACE();

  updateMacroscopicVariables();
  for(MInt i = 0; i < this->noInternalCells(); i++) {
    calculateMomentumFlux(i);
  }
  calculateSGSTensors();

  std::vector<MFloat> dynamicCs;
  constexpr MInt direction = 1;
  MInt count = 0;
  averageSGSTensors(direction, count, dynamicCs);

  for(MInt i = 0; i < count; i++) {
    if(dynamicCs[i] < 0.0) dynamicCs[i] = 0.0;
    // if(dynamicCs[i]>0.04)
    //   dynamicCs[i] = 0.04;
  }

  MFloat rho = F0;
  MFloat u[nDim];
  MFloat Q = F0;
  MFloat tau = F0;
  MFloat S = F0;

  MInt index = 0;
  MFloat actualCellLength = F0;

  MChar buf[10];
  MString fileName;
  ofstream ofl;

  // dissipation values
  ScratchSpace<MFloat> diss(m_arraySize[direction], AT_, "diss");
  ScratchSpace<MFloat> totalDiss(m_arraySize[direction], AT_, "totalDiss");
  ScratchSpace<MFloat> SGSDiss(m_arraySize[direction], AT_, "SGSDiss");
  ScratchSpace<MFloat> totalSGSDiss(m_arraySize[direction], AT_, "totalSGSDiss");
  MFloat energy;

  // reset dissipation and energy values
  for(MInt i = 0; i < count; i++) {
    diss[i] = F0;
    SGSDiss[i] = F0;
  }
  energy = F0;

  if(m_tanhInit && (globalTimeStep > m_initStartTime)) {
    MFloat scale =
        0.5 * (m_tanhScaleFactor * tanh((5.0 * (MFloat)(globalTimeStep - m_initStartTime) / (MFloat)m_initTime) - 2.5))
        + 0.5;

    m_Re = m_initRe + scale * (m_finalRe - m_initRe);

    // final timestep reached
    if(globalTimeStep == m_initStartTime + m_initTime) {
      m_Re = m_finalRe;
      m_tanhInit = 0;
    }
  }
  m_nu = m_Ma * LBCS / m_Re * m_referenceLength;

  for(MInt i = 0; i < m_currentMaxNoCells; i++) {
    const MInt pCellId = m_activeCellList[i];

    if((globalTimeStep - 1) % IPOW2(maxLevel() - this->a_level(pCellId)) == 0) {
      rho = a_variable(pCellId, PV->RHO);

      for(MInt d = 0; d < nDim; d++) {
        u[d] = a_variable(pCellId, d);
      }

      // Calculation of equilibrium and non-equilibrium distribution
      std::array<MFloat, nDist> eqDist;
      eqDist = getEqDists(rho, u);

      //--------------------------------------------
      // Calculation of overall viscosity
      // according to Hu, Sterling, Chen 1994

      // 1. Calculation of original relaxation time (on current level)
      tau = F1B2 + 3.0 * m_nu * FFPOW2(maxLevel() - this->a_level(pCellId));

      // 2. Calculation of the filtered mean momentum flux
      Q = 0.0;
      for(MInt k = 0; k < nDim * nDim; k++) {
        Q += m_momentumFlux[pCellId][k] * m_momentumFlux[pCellId][k];
      }
      Q = sqrt(2.0 * Q);

      // 3. Calculation of new relaxation time
      // index = (MInt)((a_coordinate(pCellId, direction) - minCoord) / length * MFloat(count - 1) + 0.5);
      index = floor(
          (F1B2 * count
           + (a_coordinate(pCellId, direction) - F1B2 * (actualCellLength)) / (this->c_cellLengthAtLevel(maxLevel())))
          + 0.1);

      tau += F1B2
             * (sqrt(tau * tau
                     + 2.0 * SQRT2 * dynamicCs[index] * m_deltaX * m_deltaX * (F1BCSsq * F1BCSsq) * Q
                           * FPOW2(maxLevel() - this->a_level(pCellId)))
                - tau);

      m_omega = 1.0 / tau;

      // save total nu for bndcnd and interface interpolation (save on finest level)
      a_nu(pCellId) = FPOW2(maxLevel() - this->a_level(pCellId)) * (tau - F1B2) / 3.0;

      //--------------------------------------------

      // Calculation of new distributions for directions with at least one component
      MFloat trace = F0;

      if constexpr(nDim == 2) {
        trace = m_momentumFlux[pCellId][0] + m_momentumFlux[pCellId][3];
      } else {
        trace = m_momentumFlux[pCellId][0] + m_momentumFlux[pCellId][4] + m_momentumFlux[pCellId][8];
      }

      for(MInt j = 0; j < (nDist - 1); j++) {
        MFloat t = F0;
        if(j < Ld::distFld(0)) {
          t = Ld::tp(1);
        } else if(j < Ld::distFld(1) + Ld::distFld(0)) {
          t = Ld::tp(2);
        } else {
          t = Ld::tp(3);
        }

        a_distribution(pCellId, j) = eqDist[j] - (1.0 - m_omega) * t * F1BCSsq * F1B2 * trace;

        for(MInt k = 0; k < nDim; k++) {
          for(MInt l = 0; l < nDim; l++) {
            a_distribution(pCellId, j) += (1.0 - m_omega) * t * F1BCSsq * F1BCSsq * F1B2 * (Ld::idFld(j, k) - 1)
                                          * (Ld::idFld(j, l) - 1) * m_momentumFlux[pCellId][l + nDim * k];
          }
        }
      }
      // Calculation of new distribution for rest particle distribution (center)
      a_distribution(pCellId, Ld::lastId()) =
          eqDist[Ld::lastId()] - (1.0 - m_omega) * Ld::tp(0) * F1BCSsq * F1B2 * trace;

      //---------------------------
      // save dissipation values

      // Calculate characteristic filtered rate of strain using original viscosity
      m_omega = 2.0 / (1.0 + 6.0 * m_nu * FFPOW2(maxLevel() - this->a_level(pCellId)));
      S = -F1B2 * F1BCSsq * m_omega * Q;

      // Calculate molecular- and subgrid dissipation
      diss[index] += m_nu * S * S;
      SGSDiss[index] += dynamicCs[index] * m_deltaX * m_deltaX * S * S * S * FPOW4(maxLevel() - this->a_level(pCellId));
      for(MInt d = 0; d < nDim; d++) {
        energy += F1B2 * (a_variable(pCellId, d) * a_variable(pCellId, d));
      }
      //---------------------------
    }
  }

  // write dissipation values
  if(m_calculateDissipation && globalTimeStep % m_solutionInterval == 0) {
    // for (MInt i=0; i<count; i++){
    //   m_log <<"cs["<<i<<"]="<<dynamicCs[i]<<endl;
    // }

    // all domains send their information to domain_0
    if(domainId() != 0) {
      MPI_Send(&(diss[0]), count, MPI_DOUBLE, 0, 0, mpiComm(), AT_, "(diss[0])");
      MPI_Send(&(SGSDiss[0]), count, MPI_DOUBLE, 0, 0, mpiComm(), AT_, "(SGSDiss[0])");
      MPI_Send(&energy, 1, MPI_DOUBLE, 0, 0, mpiComm(), AT_, "energy");
    }

    if(domainId() == 0) {
      MPI_Status status;

      for(MInt i = 0; i < count; i++) {
        totalDiss[i] = diss[i];
        totalSGSDiss[i] = SGSDiss[i];
      }

      for(MInt j = 1; j < noDomains(); j++) {
        MPI_Recv(&(diss[0]), count, MPI_DOUBLE, j, 0, mpiComm(), &status, AT_, "(diss[0])");
        MPI_Recv(&(SGSDiss[0]), count, MPI_DOUBLE, j, 0, mpiComm(), &status, AT_, "(SGSDiss[0])");
        MPI_Recv(&energy, 1, MPI_DOUBLE, j, 0, mpiComm(), &status, AT_, "energy");
        for(MInt i = 0; i < count; i++) {
          totalDiss[i] += diss[i];
          totalSGSDiss[i] += SGSDiss[i];
        }
      }

      fileName = "totalDiss_";
      sprintf(buf, "%d", globalTimeStep);
      fileName += buf;
      fileName += ".dat";
      m_log << "writing diss file: " << fileName << endl;
      ofl.open(fileName.c_str(), ios_base::out);
      for(MInt i = 0; i < count; i++)
        ofl << globalTimeStep << " " << i << " " << totalDiss[i] << " " << totalSGSDiss[i] << " " << energy << endl;
      ofl.close();
    }
  }
}


/** \brief Consistent initialization step of the LBGK algorithm
 *
 * The initial velocity field is held constant. Only the pressure
 * is changed by collision. (Mei et al. 2006)
 *
 * Regularized pre-collision functions are used for sake of stability and accuracy !!!
 *
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbSolverDxQy<nDim, nDist, SysEqn>::bgki_init_collision_step() {
  TRACE();

  // capable of D3Q19 and D3Q27
  for(MInt i = 0; i < m_currentMaxNoCells; i++) {
    const MInt pCellId = m_activeCellList[i];

    if((globalTimeStep - 1) % IPOW2(maxLevel() - this->a_level(pCellId)) == 0) {
      m_nu = m_Ma * LBCS / m_Re * m_referenceLength;

      m_omega = 2.0 / (1.0 + 6.0 * m_nu * FFPOW2(maxLevel() - this->a_level(pCellId)));

      swap_variables(pCellId);

      // Calculate macroscopic density
      a_variable(pCellId, PV->RHO) = 0.0;
      for(MInt j = 0; j < nDist; j++) {
        a_variable(pCellId, PV->RHO) += a_oldDistribution(pCellId, j);
      }

      const MFloat rho = a_variable(pCellId, PV->RHO);
      MFloat u[nDim];
      for(MInt d = 0; d < nDim; d++) {
        u[d] = a_variable(pCellId, d);
      }

      // Calculation of equilibrium and non-equilibrium distribution
      std::array<MFloat, nDist> eqDist;
      eqDist = getEqDists(rho, u);

      // Calculation of momentum flux tensor for non-equilibrium parts
      //--------------------------------------------
      // \Pi_{\alpha,\beta} = \sum_i f_i \cdot c_{i,\alpha}c_{i,\beta}
      calculateMomentumFlux(pCellId);

      // Calculation of new distributions for directions with at least one component
      MFloat trace = F0;

      if constexpr(nDim == 2) {
        trace = m_momentumFlux[pCellId][0] + m_momentumFlux[pCellId][3];
      } else {
        trace = m_momentumFlux[pCellId][0] + m_momentumFlux[pCellId][4] + m_momentumFlux[pCellId][8];
      }

      for(MInt j = 0; j < (nDist - 1); j++) {
        const MFloat t = Ld::tp(Ld::distType(j));

        a_distribution(pCellId, j) = eqDist[j] - (1.0 - m_omega) * t * F1BCSsq * F1B2 * trace;

        for(MInt k = 0; k < nDim; k++) {
          for(MInt l = 0; l < nDim; l++) {
            a_distribution(pCellId, j) += (1.0 - m_omega) * t * F1BCSsq * F1BCSsq * F1B2 * (Ld::idFld(j, k) - 1)
                                          * (Ld::idFld(j, l) - 1) * m_momentumFlux[pCellId][l + nDim * k];
          }
        }
      }
      // Calculation of new distribution for rest particle distribution (center)
      a_distribution(pCellId, Ld::lastId()) =
          eqDist[Ld::lastId()] - (1.0 - m_omega) * Ld::tp(0) * F1BCSsq * F1B2 * trace;
    }
  }
}

/** \brief Calculate average SGS tensor
 *
 *  \param[in]  direction   direction of the line
 *  \param[out] count       number of points in direction
 *  \param[out] meanTensors Mean SGS tensor
 *
 * Calculate average values for Cs along a line in chosen direction and return
 * mean tensor and number of points of the line.
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbSolverDxQy<nDim, nDist, SysEqn>::averageSGSTensors(const MInt direction, MInt& count,
                                                          std::vector<MFloat>& meanTensors) {
  TRACE();

  //  MInt tmpId;
  //  MInt saveId;
  //  MInt tmpCt[3];
  //  MFloat tmpMin[3];
  //  MFloat tmpLength[3];
  //
  // /***********************************************/
  // /* Count all cells in x,y, and z-direction*/
  // for(MInt i=0; i < this->grid().noInternalCells(); i++){
  //   // Ignore halo cells
  //   if(cells[i].m_level < 0 )
  //     continue;
  //   if(this->c_noChildren(i)==0){
  //     tmpId=i;
  //     // Go to lowest level
  //     while(c_parentId(tmpId) > -1){
  // 	tmpId=c_parentId(tmpId);
  //     }
  //     saveId = tmpId;
  //     for(MInt j=0; j < 2*nDim; j+=2){
  // 	tmpId = saveId;
  // 	// Go to the last cell in current direction
  // 	while(a_hasNeighbor(tmpId, j)){
  // 	  if( a_coordinate(c_neighborId(tmpId, j), j/2) < a_coordinate(tmpId, j/2) ){
  // 	    tmpId=c_neighborId(tmpId, j);
  // 	  }
  // 	  else
  // 	    break;
  // 	}
  //
  // 	tmpMin[j/2] = a_coordinate(tmpId, j/2);
  // 	tmpLength[j/2] = a_coordinate(tmpId, j/2);
  // 	// Go to the first cell in current direction
  // 	tmpCt[j/2] = 1;
  // 	while(a_hasNeighbor(tmpId, Ld::oppositeDist(j))){
  // 	  if( a_coordinate(c_neighborId(tmpId, Ld::oppositeDist(j)), j/2) > a_coordinate(tmpId, j/2) ){
  // 	    tmpId=c_neighborId(tmpId, Ld::oppositeDist(j));
  // 	    tmpCt[j/2]++;
  // 	  }
  // 	  else
  // 	    break;
  // 	}
  //
  // 	tmpLength[j/2] = (a_coordinate(tmpId, j/2) - tmpLength[j/2]);
  //     }
  //   }
  //   break;
  // }

  const MInt(&tmpCt)[nDim] = m_arraySize;
  count = tmpCt[direction];

  ScratchSpace<MFloat> ML(tmpCt[direction], AT_, "ML");
  ScratchSpace<MFloat> MM(tmpCt[direction], AT_, "MM");
  meanTensors.resize(tmpCt[direction]);

  for(MInt i = 0; i < tmpCt[direction]; i++) {
    MM[i] = F0;
    ML[i] = F0;
    meanTensors[i] = F0;
  }

  // Sum up values
  for(MInt i = 0; i < this->grid().noInternalCells(); i++) {
    const MFloat actualCellLength = this->grid().cellLengthAtCell(i);

    // const MInt index = (MInt)((a_coordinate(i, direction) - tmpMin[direction]) / (MFloat)tmpLength[direction]
    //                              * MFloat(tmpCt[direction] - 1) + 0.5);

    const MInt index =
        floor((F1B2 * tmpCt[direction]
               + (a_coordinate(i, direction) - F1B2 * (actualCellLength)) / (this->c_cellLengthAtLevel(maxLevel())))
              + 0.1);

    if(index < 0 || index >= tmpCt[direction]) {
      cerr << "error: index=" << index << endl;
      m_log << "error: index=" << index << endl;
    }

    ML[index] += m_MijLij[i];
    MM[index] += m_MijMij[i];
  }

  m_log.precision(12);

  // Sum-Up over all domains
  MPI_Allreduce(MPI_IN_PLACE, MM.data(), count, MPI_DOUBLE, MPI_SUM, mpiComm(), AT_, "MPI_IN_PLACE", "MM");
  MPI_Allreduce(MPI_IN_PLACE, ML.data(), count, MPI_DOUBLE, MPI_SUM, mpiComm(), AT_, "MPI_IN_PLACE", "ML");

  // // Normalize
  // for(MInt i=0; i < tmpCt[direction]; i++){
  // 	MM[i] *= ((MFloat)tmpCt[direction] / (MFloat)(tmpCt[0]*tmpCt[1]*tmpCt[2])) ;
  // 	ML[i] *= ((MFloat)tmpCt[direction] / (MFloat)(tmpCt[0]*tmpCt[1]*tmpCt[2])) ;
  // }

  // Divide
  for(MInt i = 2; i < tmpCt[direction] - 2; i++) {
    if(MM[i] > 1e-12) meanTensors[i] = ML[i] / MM[i];
  }
}

/** \brief  Update viscosity (a_nu and a_oldNu)
 *  \author Miro Gondrum
 *  \date   03.01.2024
 *  Reset viscoisty to m_nu default value or according to tanh ramping
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbSolverDxQy<nDim, nDist, SysEqn>::updateViscosity() {
  TRACE();
  // Update the global Reynolds number if we use the tanh-adaption
  if(m_tanhInit && (globalTimeStep > m_initStartTime)) {
    const MFloat scale =
        0.5 * (m_tanhScaleFactor * tanh((5.0 * (MFloat)(globalTimeStep - m_initStartTime) / (MFloat)m_initTime) - 2.5))
        + 0.5;

    m_Re = m_initRe + scale * (m_finalRe - m_initRe);
    // final timestep reached
    if(globalTimeStep == m_initStartTime + m_initTime) {
      m_Re = m_finalRe;
      m_tanhInit = false;
    }
  }
  // Required for ported functions, since GPUs cannot access the global variable globalTimeStep
  const MInt gTS = globalTimeStep;
  const MInt maxLevel_ = maxLevel();
  // update/reset m_nu and a_nu/a_oldNu
  m_nu = m_Ma * LBCS / m_Re * m_referenceLength;
  if(m_cells.saveOldNu()) {
    maia::parallelFor<true>(0, m_currentMaxNoCells, [=](MInt index) {
      const MInt pCellId = m_activeCellList[index];
      const MInt lvlDiff = maxLevel_ - a_level(pCellId);
      if((gTS - 1) % IPOW2(lvlDiff) != 0) return;
      this->a_oldNu(pCellId) = a_nu(pCellId);
      a_nu(pCellId) = m_nu;
    });
  } else {
    maia::parallelFor<true>(0, m_currentMaxNoCells, [=](MInt index) {
      const MInt pCellId = m_activeCellList[index];
      const MInt lvlDiff = maxLevel_ - a_level(pCellId);
      if((gTS - 1) % IPOW2(lvlDiff) != 0) return;
      a_nu(pCellId) = m_nu;
    });
  }
}

/** \brief Update macroscopic variables according to incoming PPDF
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbSolverDxQy<nDim, nDist, SysEqn>::updateMacroscopicVariables() {
  TRACE();

  for(MInt i = 0; i < m_currentMaxNoCells; i++) {
    const MInt pCellId = m_activeCellList[i];
    swap_variables(pCellId);

    // Update macroscopic variables
    a_variable(pCellId, PV->RHO) = 0.0;
    for(MInt j = 0; j < nDist; j++) {
      a_variable(pCellId, PV->RHO) += a_oldDistribution(pCellId, j);
    }
    for(MInt d = 0; d < nDim; d++) {
      a_variable(pCellId, PV->U + d) = 0.0;
      for(MInt j = 0; j < Ld::dxQyFld(); j++) {
        a_variable(pCellId, PV->U + d) += a_oldDistribution(pCellId, Ld::pFld(d, j));
        a_variable(pCellId, PV->U + d) -= a_oldDistribution(pCellId, Ld::nFld(d, j));
      }
    }
  }
}

/** \brief Calculate tensors for dynamic Smagorinsky constant
 *
 * only for D3Q19 and uniform grids
 *
 * Calculation of MijMij and MijLij according to Pope p.620.
 * The test filter is a box filter with a width 2*m_deltaX.
 *
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbSolverDxQy<nDim, nDist, SysEqn>::calculateSGSTensors() {
  TRACE();

  constexpr MFloat coefficients[4] = {F1B8, F1B16, F1B32, F1B64};
  m_nu = m_Ma * LBCS / m_Re * m_referenceLength;

  for(MInt id = 0; id < this->grid().noInternalCells(); id++) {
    //----------------------------------------
    // 1. calculate amount from actual cell
    // Calculate strain rate tensor
    // S_{ij} = \frac{1}{2} (\partial_j u_i + \partial_i u_j)
    m_omega = 2.0 / (1.0 + 6.0 * m_nu * FPOW2(maxLevel() - this->a_level(id)));
    MFloat Sij[9];
    for(MInt k = 0; k < 9; k++) {
      Sij[k] = -F1B2 * F1BCSsq * m_omega * m_momentumFlux[id][k];
    }
    // Calculate caracteristic filtered rate of strain
    MFloat S = sqrt(2.0 * std::inner_product(&Sij[0], &Sij[9], &Sij[0], .0));

    // add values to sum
    MFloat Mij[9]{};
    MFloat SijTilde[9]{};
    for(MInt j = 0; j < 9; j++) {
      Mij[j] += coefficients[0] * 2.0 * m_deltaX * m_deltaX * S * Sij[j];
      SijTilde[j] += coefficients[0] * Sij[j];
    }
    MFloat uuTilde[9]{};
    for(MInt k = 0; k < 3; k++) {
      for(MInt l = 0; l < 3; l++) {
        uuTilde[l + 3 * k] += coefficients[0] * a_variable(id, k) * a_variable(id, l);
      }
    }
    MFloat uTilde[3]{};
    for(MInt k = 0; k < 3; k++) {
      uTilde[k] += coefficients[0] * a_variable(id, k);
    }
    MFloat STilde = coefficients[0] * S;

    //----------------------------------------
    // 2. add amount from neighbor cells
    // neighbors in direction with one component
    for(MInt j = 0; j < Ld::distFld(0); j++) {
      if(!a_hasNeighbor(id, j)) continue;

      const MInt currentNghbr = c_neighborId(id, j);

      // Calculate strain rate tensor
      m_omega = 2.0 / (1.0 + 6.0 * m_nu * FPOW2(maxLevel() - this->a_level(currentNghbr)));
      for(MInt k = 0; k < 9; k++) {
        Sij[k] = -F1B2 * F1BCSsq * m_omega * m_momentumFlux[currentNghbr][k];
      }
      // Calculate caracteristic filtered rate of strain
      S = sqrt(2.0 * std::inner_product(&Sij[0], &Sij[9], &Sij[0], .0));

      // add values to sum
      for(MInt k = 0; k < 9; k++) {
        Mij[k] += coefficients[1] * 2.0 * m_deltaX * m_deltaX * S * Sij[k];
        SijTilde[k] += coefficients[1] * Sij[k];
      }
      for(MInt k = 0; k < 3; k++) {
        for(MInt l = 0; l < 3; l++) {
          uuTilde[l + 3 * k] += coefficients[1] * a_variable(id, k) * a_variable(id, l);
        }
      }
      for(MInt k = 0; k < 3; k++) {
        uTilde[k] += coefficients[1] * a_variable(id, k);
      }
      STilde += coefficients[1] * S;
    }
    // neighbors in direction with two components
    MInt tmpDir = Ld::distFld(0);
    for(MInt j = 0; j < Ld::distFld(1); j++) {
      if(!a_hasNeighbor(id, tmpDir + j)) continue;

      const MInt currentNghbr = c_neighborId(id, tmpDir + j);

      // Calculate strain rate tensor
      m_omega = 2.0 / (1.0 + 6.0 * m_nu * FPOW2(maxLevel() - this->a_level(currentNghbr)));
      for(MInt k = 0; k < 9; k++) {
        Sij[k] = -F1B2 * F1BCSsq * m_omega * m_momentumFlux[currentNghbr][k];
      }
      // Calculate caracteristic filtered rate of strain
      S = sqrt(2.0 * std::inner_product(&Sij[0], &Sij[9], &Sij[0], .0));

      // add values to sum
      for(MInt i = 0; i < 9; i++) {
        Mij[i] += coefficients[2] * 2.0 * m_deltaX * m_deltaX * S * Sij[i];
        SijTilde[i] += coefficients[2] * Sij[i];
      }
      for(MInt k = 0; k < 3; k++) {
        for(MInt l = 0; l < 3; l++) {
          uuTilde[l + 3 * k] += coefficients[2] * a_variable(id, k) * a_variable(id, l);
        }
      }
      for(MInt k = 0; k < 3; k++) {
        uTilde[k] += coefficients[2] * a_variable(id, k);
      }
      STilde += coefficients[2] * S;
    }
    // neighbors in direction with three components
    tmpDir = Ld::distFld(0) + Ld::distFld(1);
    for(MInt j = 0; j < Ld::distFld(2); j++) {
      if(!a_hasNeighbor(id, tmpDir + j)) continue;

      const MInt currentNghbr = c_neighborId(id, tmpDir + j);

      // Calculate strain rate tensor
      m_omega = 2.0 / (1.0 + 6.0 * m_nu * FPOW2(maxLevel() - this->a_level(currentNghbr)));
      for(MInt k = 0; k < 9; k++) {
        Sij[k] = -F1B2 * F1BCSsq * m_omega * m_momentumFlux[currentNghbr][k];
      }
      // Calculate caracteristic filtered rate of strain
      S = sqrt(2.0 * std::inner_product(&Sij[0], &Sij[9], &Sij[0], .0));

      // add values to sum
      for(MInt k = 0; k < 9; k++) {
        Mij[k] += coefficients[3] * 2.0 * m_deltaX * m_deltaX * S * Sij[k];
        SijTilde[k] += coefficients[3] * Sij[k];
      }
      for(MInt k = 0; k < 3; k++) {
        for(MInt l = 0; l < 3; l++) {
          uuTilde[l + 3 * k] += coefficients[3] * a_variable(id, k) * a_variable(id, l);
        }
      }
      for(MInt k = 0; k < 3; k++) {
        uTilde[k] += coefficients[3] * a_variable(id, k);
      }
      STilde += coefficients[3] * S;
    }
    //----------------------------------------


    // finalize Mij and Lij
    for(MInt j = 0; j < 9; j++) {
      Mij[j] -= 2.0 * 4.0 * m_deltaX * m_deltaX * STilde * SijTilde[j];
    }

    MFloat Lij[9]{};
    for(MInt k = 0; k < 3; k++) {
      for(MInt l = 0; l < 3; l++) {
        Lij[l + 3 * k] += uuTilde[l + 3 * k] - uTilde[k] * uTilde[l];
      }
    }

    // Multiply Mij by itself
    m_MijMij[id] = 0.0;
    for(MInt k = 0; k < 9; k++) {
      m_MijMij[id] += Mij[k] * Mij[k];
    }

    // Multiply Mij by Lij
    m_MijLij[id] = 0.0;
    for(MInt k = 0; k < 9; k++) {
      m_MijLij[id] += Mij[k] * Lij[k];
    }
  }
}

/** \brief  Initialize the source term controller
 *  \author Miro Gondrum
 *  \date   01.02.2022
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbSolverDxQy<nDim, nDist, SysEqn>::initSrcTermController() {
  TRACE();
  m_srcTermController.init();
}

/** \brief  Initialize the source term controller
 *  \author Julian Vorspohl
 *  \date   01.05.2024
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbSolverDxQy<nDim, nDist, SysEqn>::initSrcTerms() {
  TRACE();
  m_srcTermController.initSrcTerms();
}

/** \brief  Calls the pre collision routine of the source term controller
 *  \author Miro Gondrum
 *  \date   01.02.2022
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbSolverDxQy<nDim, nDist, SysEqn>::preCollisionSrcTerm() {
  TRACE();
  m_srcTermController.apply_preCollision();
}

/** \brief  Calls the post collision routine of the source term controller
 *  \author Miro Gondrum
 *  \date   01.02.2022
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbSolverDxQy<nDim, nDist, SysEqn>::postCollisionSrcTerm() {
  TRACE();
  m_srcTermController.apply_postCollision();
}

/** \brief  Calls the post collision routine of the source term controller
 *  \author Julian Vorspohl
 *  \date   01.05.2024
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbSolverDxQy<nDim, nDist, SysEqn>::postPropagationSrcTerm() {
  TRACE();
  m_srcTermController.apply_postPropagation();
}

/** \brief
 *
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbSolverDxQy<nDim, nDist, SysEqn>::postCollisionBc() {
  TRACE();

  m_bndCnd->updateVariables();
}

/** \brief
 *
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbSolverDxQy<nDim, nDist, SysEqn>::postPropagationBc() {
  TRACE();

  m_bndCnd->updateRHS();
}

/** \brief
 *
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbSolverDxQy<nDim, nDist, SysEqn>::restartInitLb() {
  TRACE();
  initEqDistFunctions();

  if(m_isThermal) {
    initThermalEqDistFunctions();
  }
  if(m_isTransport) {
    initTransportEqDistFunctions();
  }
}

/** \brief Initializes standard Lattice BGK
 * \author Andreas Lintermann
 * \date 24.01.2011
 *
 * This function sets a function pointer for the right initailization method and then finally calls that method. The
 * decision is based on the property "initMethod" and decides between an initialization with and without FFT.
 *
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbSolverDxQy<nDim, nDist, SysEqn>::initializeLatticeBgk() {
  TRACE();

  m_log << "  + Initializing flow field..." << endl;
  m_log << "    - init method:  " << m_initMethod << endl;

  m_smallestCellLength = (this->c_cellLengthAtLevel(maxLevel())); // macroscopic cell-Length on highest level
  m_nu = m_Ma * LBCS / m_Re * m_referenceLength;                  // viscosity on highest level
  m_omega = 2.0 / (1.0 + 6.0 * m_nu);

  if(m_FftInit) {
    m_log << "    - flow init:    turbulent" << endl;

    switch(string2enum(m_initMethod)) {
      case LB_TURBULENT_PIPE_INIT:
        m_initMethodPtr = &LbSolverDxQy::initLatticeBgkFftPipe;
        break;
      case LB_TURBULENT_CHANNEL_INIT:
        m_initMethodPtr = &LbSolverDxQy::initLatticeBgkFftChannel;
        break;
      case LB_TURBULENT_MIXING_INIT:
        m_initMethodPtr = &LbSolverDxQy::initLatticeBgkFftMixing;
        break;
      case LB_TURBULENT_MIXING_FILTER_INIT:
        m_initMethodPtr = &LbSolverDxQy::initLatticeBgkFftMixingFilter;
        break;
      case LB_TURBULENCE_ISOTROPIC_INIT:
        m_initMethodPtr = &LbSolverDxQy::initLatticeBgkFftIsotropicTurbulence;
        break;
      default: {
        m_log << "lbsolverdxqy::() unknown Init-Method type in combination with FftInit" << endl;
        cerr << "lbsolverdxqy::() unknown Init-Method type in combination with FftInit" << endl;
        DEBUG("lbsolverdxqy::() unknown Init-Method type in combination with FftInit" << endl, MAIA_DEBUG_ASSERTION);
        TERMM(1, "unknown Init-Method type in combination with FftInit");
      }
    }
  } else {
    m_log << "    - flow init:    laminar" << endl;

    switch(string2enum(m_initMethod)) {
      case LB_LAMINAR_CHANNEL_INIT:
        m_initMethodPtr = &LbSolverDxQy::initLatticeBgkLaminarChannel;
        break;
      case LB_TURBULENT_CHANNEL_INIT:
        m_initMethodPtr = &LbSolverDxQy::initLatticeBgkTurbulentChannel;
        break;
      case LB_TURBULENT_MIXING_INIT:
        m_initMethodPtr = &LbSolverDxQy::initLatticeBgkTurbulentMixing;
        break;
      case LB_TURBULENT_BOUNDARY:
        m_initMethodPtr = &LbSolverDxQy::initLatticeBgkTurbulentBoundary;
        break;
      case LB_TURBULENT_PIPE_INIT:
        m_initMethodPtr = &LbSolverDxQy::initLatticeBgkTurbulentPipe;
        break;
      case LB_TURBULENT_DUCT_INIT:
        m_initMethodPtr = &LbSolverDxQy::initLatticeBgkTurbulentDuct;
        break;
      case LB_FROM_ZERO_INIT:
      case LB_LAMINAR_INIT_PX:
      case LB_LAMINAR_INIT_MX:
      case LB_LAMINAR_INIT_PY:
      case LB_LAMINAR_INIT_MY:
      case LB_LAMINAR_INIT_PZ:
      case LB_LAMINAR_INIT_MZ:
        m_initMethodPtr = &LbSolverDxQy::initLatticeBgkLaminar;
        break;
      case LB_LAMINAR_PIPE_INIT:
        m_initMethodPtr = &LbSolverDxQy::initLatticeBgkLaminarPipe;
        break;
      case LB_TGV_INIT:
        m_initMethodPtr = &LbSolverDxQy::initLatticeBgkTGV;
        break;
      case LB_GAUSS_PULSE_INIT:
        m_initMethodPtr = &LbSolverDxQy::initLatticeBgkGaussPulse;
        break;
      case LB_GAUSS_DIFFUSION_INIT:
        m_initMethodPtr = &LbSolverDxQy::initLatticeBgkGaussDiffusion;
        break;
      case LB_GAUSS_ADVECTION_INIT:
        m_initMethodPtr = &LbSolverDxQy::initLatticeBgkGaussAdvection;
        break;
      case LB_LAMINAR_CYLINDER_INIT:
        m_initMethodPtr = &LbSolverDxQy::initLatticeBgkLaminarCylinder;
        break;
      case LB_SPINNING_VORTICIES_INIT:
        m_initMethodPtr = &LbSolverDxQy::initLatticeBgkSpinningVorticies;
        break;
      case LB_CONVECTING_VORTEX_INIT:
        m_initMethodPtr = &LbSolverDxQy::initLatticeBgkVortex;
        break;
      case LB_STEADY_VORTEX_INIT:
        m_initMethodPtr = &LbSolverDxQy::initLatticeBgkVortex;
        break;
      default: {
        m_log << "lbsolverdxqy::() unknown Init-Method type" << endl;
        cerr << "lbsolverdxqy::() unknown Init-Method type" << endl;
        DEBUG("lbsolverdxqy::() unknown Init-Method type" << endl, MAIA_DEBUG_ASSERTION);
        TERMM(1, "unknown Init-Method type");
      }
    }
  }

  if(m_densityFluctuations)
    m_log << "    - init density: fluctuations enanbled" << endl;
  else
    m_log << "    - init density: absolute" << endl;

  (this->*m_initMethodPtr)();
}

/** \fn LbSolverDxQy::initLatticeBgkLaminarChannel()
 *  \author Andreas Lintermann
 *  \date 24.01.2014
 *
 *  \propVal{initMethod,LB_LAMINAR_CHANNEL_INIT}
 *
 * Initializes standard Lattice BGK Laminar Channel.
 * Set property "initMethod" to "LB_LAMINAR_CHANNEL_INIT"
 *
 * \LBIC{LbSolverDxQy::initLatticeBgkLaminarChannel(), initLatticeBgkLaminarChannel, LaminarChannelInit}
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbSolverDxQy<nDim, nDist, SysEqn>::initLatticeBgkLaminarChannel() {
  TRACE();

  // all cells have the same density
  const MFloat rho = (m_densityFluctuations) ? 0.0 : 1.0;

  for(MInt i = 0; i < a_noCells(); i++) {
    a_variable(i, PV->RHO) = rho;
    a_oldVariable(i, PV->RHO) = rho;

    // the viscosity is saved for each cell
    // this is needed for subgrid modelling
    initNu(i, m_nu);
    if(m_isThermal) {
      a_kappa(i) = m_kappa;
      a_variable(i, PV->T) = 1.0;
      a_oldVariable(i, PV->T) = 1.0;
    }
    if(m_isTransport) {
      a_diffusivity(i) = m_diffusivity;
      a_variable(i, PV->C) = 1.0;
      a_oldVariable(i, PV->C) = 1.0;
    }

    for(MInt dir = 1; dir < nDim; dir++) {
      a_variable(i, PV->VV[dir]) = 0.0;
    }

    std::array<MFloat, 2 * nDim> bBox{};
    this->m_geometry->getBoundingBox(bBox.data());
    if(this->m_nonNewtonian) {
      const MFloat tmpWidth = fabs(bBox[nDim + 1] - bBox[1]) * F1B2;
      const MFloat exp = F1 + F1 / this->m_n;

      const MFloat relPos = fabs(a_coordinate(i, 1) / tmpWidth - F1);
      const MFloat parabola = F1 - pow(relPos, exp);

      a_variable(i, PV->VV[0]) = m_Ma * LBCS * parabola;
    } else {
      const MFloat halfWidth = std::max(fabs(bBox[1 + nDim]), fabs(bBox[1]));
      a_variable(i, PV->VV[0]) =
          m_Ma * LBCS
          * ((1 - this->m_CouettePoiseuilleRatio) * a_coordinate(i, 1) / halfWidth + this->m_CouettePoiseuilleRatio
             - this->m_CouettePoiseuilleRatio * (a_coordinate(i, 1) / halfWidth) * (a_coordinate(i, 1) / halfWidth));
      a_oldVariable(i, PV->VV[0]) = a_variable(i, PV->VV[0]);
    }

    for(MInt d = 0; d < nDim; d++) {
      a_oldVariable(i, PV->VV[d]) = a_variable(i, PV->VV[d]);
    }
  }

  initEqDistFunctions();

  if(m_isThermal) {
    initThermalEqDistFunctions();
  }
  if(m_isTransport) {
    initTransportEqDistFunctions();
  }

  if(this->m_nonNewtonian) {
    initNonEqDistFunctions();
    this->exchange();
    this->exchangeOldDistributions();
  }
}

/** \fn LbSolverDxQy::initLatticeBgkLaminarPipe()
 *  \author Andreas Lintermann
 *  \date 24.01.2014
 *
 *  \propVal{initMethod,LB_LAMINAR_PIPE_INIT}
 *
 * Initializes standard Lattice BGK Laminar Pipe.
 * Set property "initMethod" to "LB_LAMINAR_PIPE_INIT"
 *
 * \LBIC{LbSolverDxQy::initLatticeBgkLaminarPipe(), initLatticeBgkLaminarPipe, LaminarPipeInit}
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbSolverDxQy<nDim, nDist, SysEqn>::initLatticeBgkLaminarPipe() {
  TRACE();

  MFloat yPlus;

  ScratchSpace<MFloat> bBox(nDim * 2, AT_, "bBox");
  MIntScratchSpace o_dirs(nDim - 1, AT_, "o_dirs");
  MFloat* bBoxPtr = &bBox[0];
  m_geometry->getBoundingBox(bBoxPtr);

  MFloat rho;

  MInt mainDir = 0;
  MInt otherDirs = 0;
  MFloat maxLength = F0;
  for(MInt d = 0; d < nDim; d++) {
    if(fabs(bBox[d + nDim] - bBox[d]) > maxLength) {
      mainDir = d;
      maxLength = fabs(bBox[d + nDim] - bBox[d]);
    } else {
      o_dirs[otherDirs] = d;
      otherDirs++;
    }
  }

  // all cells have the same density
  if(m_densityFluctuations)
    rho = 0.0;
  else
    rho = 1.0;

  for(MInt i = 0; i < a_noCells(); i++) {
    MFloat radius = 0.0;
    for(MInt dir = 0; dir < nDim; dir++) {
      if(dir == mainDir) continue;
      radius += a_coordinate(i, dir) * a_coordinate(i, dir);
    }
    radius = sqrt(radius);

    const MFloat tmpWidth = 0.5 * fabs(bBox[o_dirs[0] + nDim] - bBox[o_dirs[0]]);

    yPlus = radius / tmpWidth;

    for(MInt dir = 0; dir < nDim; dir++) {
      a_variable(i, PV->VV[dir]) = 0.0;
    }
    a_variable(i, PV->VV[mainDir]) = 2.0 * (1.0 - yPlus * yPlus) * m_Ma * LBCS; // v_max = 2*v_mean
    if(m_isThermal) {
      a_variable(i, PV->T) =
          this->m_initTemperatureKelvin - (this->m_initTemperatureKelvin - F1) * (1.0 - yPlus * yPlus);
    }
    if(m_isTransport) {
      a_variable(i, PV->C) = this->m_initCon - (this->m_initCon - F1) * (1.0 - yPlus * yPlus);
    }

    a_variable(i, PV->RHO) = rho;

    // the viscosity is saved for each cell
    // this is needed for subgrid modelling
    if(m_isThermal) {
      a_kappa(i) = m_kappa;
      a_variable(i, PV->T) = 1.0;
      a_oldVariable(i, PV->T) = 1.0;
    }
    initNu(i, m_nu);
  }

  initEqDistFunctions();

  if(m_isThermal) {
    initThermalEqDistFunctions();
  }
  if(m_isTransport) {
    initTransportEqDistFunctions();
  }
}

/** \fn LbSolverDxQy::initLatticeBgkTGV()
 *  \author Moritz Waldmann
 *  \date 29.11.2019
 *
 *  \propVal{initMethod,LB_TGV_INIT}
 *
 * Initializes standard Lattice BGK TGV.
 * Set property "initMethod" to "LB_TGV_INIT"
 *
 * \LBIC{LbSolverDxQy::initLatticeBgkTGV(), initLatticeBgkTGV, TGVlInit}
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbSolverDxQy<nDim, nDist, SysEqn>::initLatticeBgkTGV() {
  TRACE();
  IF_CONSTEXPR(nDim != 3) TERMM(1, "Only implemented for 3D");

  MFloat rho, rho0, u[3];
  MFloat p, p0;

  // all cells have the same density
  rho0 = F1;
  p0 = rho0 * F1B3;

  for(MInt i = 0; i < a_noCells(); i++) {
    const MFloat xCoord = a_coordinate(i, 0);
    const MFloat yCoord = a_coordinate(i, 1);
    const MFloat zCoord = a_coordinate(i, 2);

    const MFloat x = (cos(F2 * zCoord) + F2) * (cos(F2 * xCoord) + cos(F2 * yCoord));
    const MFloat y = sin(xCoord) * cos(yCoord) * cos(zCoord);
    const MFloat z = cos(xCoord) * sin(yCoord) * cos(zCoord);

    p = p0 + rho0 * (m_Ma * m_Ma) * F1B3 / 16.0 * x;
    rho = F3 * p;
    u[0] = m_Ma * LBCS * y;
    u[1] = -m_Ma * LBCS * z;
    u[2] = F0;

    a_variable(i, PV->RHO) = rho;

    // the viscosity is saved for each cell
    // this is necessary for subgrid modelling
    initNu(i, m_nu);

    a_variable(i, PV->U) = u[0];
    a_variable(i, PV->V) = u[1];
    a_variable(i, PV->W) = u[2];
  }
  initEqDistFunctions();
}

/** \fn LbSolverDxQy::initLatticeBgkTurbulentBoundary()
 *  \author Andreas Lintermann
 *  \date 24.01.2014
 *
 *  \propVal{initMethod,LB_TURBULENT_BOUNDARY}
 *
 * Initializes standard Lattice BGK Turbulent Boundary.
 * Set property "initMethod" to "LB_TURBULENT_BOUNDARY"
 *
 * \LBIC{LbSolverDxQy::initLatticeBgkTurbulentBoundary(), initLatticeBgkTurbulentBoundary, TurbulentBoundaryInit}
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbSolverDxQy<nDim, nDist, SysEqn>::initLatticeBgkTurbulentBoundary() {
  TRACE();

  MFloat yPlus;
  MFloat amp = 0.1;
  MFloat uTau = (MFloat)m_ReTau * m_nu / m_referenceLength;
  ScratchSpace<MFloat> bBox(nDim * 2, AT_, "bBox");
  MFloat* bBoxPtr = &bBox[0];
  m_geometry->getBoundingBox(bBoxPtr);

  MFloat rho;

  // all cells have the same density
  if(m_densityFluctuations)
    rho = 0.0;
  else
    rho = 1.0;

  const MFloat tmpWidth = fabs(bBox[1 + nDim] - bBox[1]);

  for(MInt i = 0; i < a_noCells(); i++) {
    a_variable(i, PV->RHO) = rho;

    yPlus = m_ReTau * (fabs(a_coordinate(i, 1)) / tmpWidth + 0.5);

    if(yPlus <= 5.0) {
      a_variable(i, PV->U) = uTau * yPlus;
    } else {
      if(yPlus <= 30 && yPlus > 5.0) {
        a_variable(i, PV->U) = uTau * (C1 * log(yPlus) + C2);
      } else {
        if(yPlus > 30) {
          a_variable(i, PV->U) = uTau * (C3 * log(yPlus) + C4);
        }
      }
    }

    //      deltaU =amp * 2.0 * (0.5 - 0.2*(MFloat) (1.0*rand()/(RAND_MAX+1.0))) * u[0];
    a_variable(i, PV->U) += amp * 2.0 * (0.5 - (MFloat)(1.0 * rand() / (RAND_MAX + 1.0))) * a_variable(i, PV->U);

    for(MInt dir = 1; dir < nDim; dir++) {
      a_variable(i, PV->VV[dir]) = amp * 2.0 * (0.5 - (MFloat)(1.0 * rand() / (RAND_MAX + 1.0))) * a_variable(i, PV->U);
    }

    // the viscosity is saved for each cell
    // this is needed for subgrid modelling
    initNu(i, m_nu);
  }
  initEqDistFunctions();
}

/** \fn LbSolverDxQy::initLatticeBgkTurbulentChannel()
 *  \author Andreas Lintermann
 *  \date 24.01.2014
 *
 *  \propVal{initMethod,LB_TURBULENT_CHANNEL_INIT}
 *
 * Initializes standard Lattice BGK Turbulent Channel.
 * Set property "initMethod" to "LB_TURBULENT_CHANNEL_INIT"
 *
 * \LBIC{LbSolverDxQy::initLatticeBgkTurbulentChannel(), initLatticeBgkTurbulentChannel, TurbulentChannelInit}
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbSolverDxQy<nDim, nDist, SysEqn>::initLatticeBgkTurbulentChannel() {
  TRACE();

  MFloat yPlus;
  MFloat amp = 0.1;
  MFloat uTau = (MFloat)m_ReTau * m_nu / m_referenceLength;

  MFloat rho;

  // all cells have the same density
  if(m_densityFluctuations)
    rho = 0.0;
  else
    rho = 1.0;

  // tmpWidth = 0.5 * fabs(bBox[1 + nDim]-bBox[1]);
  // tmpwidth equals diameter * cellLength (on maxLevel)
  const MFloat tmpWidth = m_referenceLength * this->c_cellLengthAtLevel(maxLevel());

  for(MInt i = 0; i < a_noCells(); i++) {
    a_variable(i, PV->RHO) = rho;

    yPlus = m_ReTau * (1.0 - fabs(a_coordinate(i, 1)) / tmpWidth);

    if(yPlus <= 5.0) {
      a_variable(i, PV->U) = uTau * yPlus;
    } else {
      if(yPlus <= 30 && yPlus > 5.0) {
        a_variable(i, PV->U) = uTau * (C1 * log(yPlus) + C2);
      } else {
        if(yPlus > 30) {
          a_variable(i, PV->U) = uTau * (C3 * log(yPlus) + C4);
        }
      }
    }

    //      deltaU =amp * 2.0 * (0.5 - 0.2*(MFloat) (1.0*rand()/(RAND_MAX+1.0))) * u[0];
    a_variable(i, PV->U) += amp * 2.0 * (0.5 - (MFloat)(1.0 * rand() / (RAND_MAX + 1.0))) * a_variable(i, PV->U);

    for(MInt dir = 1; dir < nDim; dir++) {
      a_variable(i, PV->VV[dir]) = amp * 2.0 * (0.5 - (MFloat)(1.0 * rand() / (RAND_MAX + 1.0))) * a_variable(i, PV->U);
    }

    // the viscosity is saved for each cell
    // this is needed for subgrid modelling
    initNu(i, m_nu);
  }
  initEqDistFunctions();
}

/** \fn LbSolverDxQy::initLatticeBgkTurbulentDuct()
 *  \author Andreas Lintermann
 *  \date 24.01.2014
 *
 *  \propVal{initMethod,LB_TURBULENT_DUCT_INIT}
 *
 * Initializes standard Lattice BGK Turbulent Duct.
 * Set property "initMethod" to "LB_TURBULENT_DUCT_INIT"
 *
 * \LBIC{LbSolverDxQy::initLatticeBgkTurbulentDuct(), initLatticeBgkTurbulentDuct, TurbulentDuctInit}
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbSolverDxQy<nDim, nDist, SysEqn>::initLatticeBgkTurbulentDuct() {
  TRACE();
  IF_CONSTEXPR(nDim != 3) TERMM(1, "Only implemented for 3D. For 2D use a channel!");

  MFloat deltaU, deltaW; // deltaV;
  MFloat yPlus;
  MFloat amp = 0.1;
  MFloat uTau = (MFloat)m_ReTau * m_nu / m_referenceLength;
  MFloat radius = 0;

  ScratchSpace<MFloat> bBox(nDim * 2, AT_, "bBox");
  MFloat* bBoxPtr = &bBox[0];
  m_geometry->getBoundingBox(bBoxPtr);

  MFloat rho;

  // all cells have the same density
  if(m_densityFluctuations)
    rho = 0.0;
  else
    rho = 1.0;

  const MFloat tmpWidth = fabs(bBox[1 + nDim] - bBox[1]);

  for(MInt i = 0; i < a_noCells(); i++) {
    radius = sqrt(a_coordinate(i, 1) * a_coordinate(i, 1) + a_coordinate(i, 2) * a_coordinate(i, 2));

    yPlus = (1.0 - radius / tmpWidth * 2.0) * m_ReTau;
    if(yPlus <= 5.0) {
      a_variable(i, PV->U) = uTau * yPlus;
    } else {
      if(yPlus <= 30 && yPlus > 5.0) {
        a_variable(i, PV->U) = uTau * (C1 * log(yPlus) + C2);
      } else {
        if(yPlus > 30) {
          a_variable(i, PV->U) = uTau * (C3 * log(yPlus) + C4);
        }
      }
    }

    deltaU = amp * 2.0 * (0.5 - (MFloat)(1.0 * rand() / (RAND_MAX + 1.0))) * a_variable(i, PV->U);
    a_variable(i, PV->U) += deltaU;

    // deltaV =amp * 2.0 * (0.5 - (MFloat) (1.0*rand()/(RAND_MAX+1.0))) * a_variable(i, PV->U);

    deltaW = amp * 2.0 * (0.5 - (MFloat)(1.0 * rand() / (RAND_MAX + 1.0))) * a_variable(i, PV->U);
    a_variable(i, PV->V) = 0.0;
    a_variable(i, PV->W) = deltaW;
    a_variable(i, PV->RHO) = rho;

    // the viscosity is saved for each cell
    // this is needed for subgrid modelling
    initNu(i, m_nu);
  }
  initEqDistFunctions();
}

/** \fn LbSolverDxQy::initLatticeBgkTurbulentMixing()
 *  \author Andreas Lintermann
 *  \date 24.01.2014
 *
 *  \propVal{initMethod,LB_TURBULENT_MIXING_INIT}
 *
 * Initializes standard Lattice BGK Turbulent Mixing Layer.
 * Set property "initMethod" to "LB_TURBULENT_MIXING_INIT"
 *
 * \LBIC{LbSolverDxQy::initLatticeBgkTurbulentMixing(), initLatticeBgkTurbulentMixing, TurbulentMixingInit}
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbSolverDxQy<nDim, nDist, SysEqn>::initLatticeBgkTurbulentMixing() {
  TRACE();
  IF_CONSTEXPR(nDim != 3) mTerm(1, AT_, "Only implemented for 3D. For 2D use a channel!");

  MFloat amp = 0.2 * m_Ma * LBCS;
  MFloat GaussFactor, sigma, deltaOmega0;
  MFloat rho;

  ScratchSpace<MFloat> bBox(nDim * 2, AT_, "bBox");
  MFloat* bBoxPtr = &bBox[0];
  m_geometry->getBoundingBox(bBoxPtr);

  // all cells have the same density
  if(m_densityFluctuations)
    rho = 0.0;
  else
    rho = 1.0;

  // initial vorticity thickness
  deltaOmega0 = m_referenceLength * m_smallestCellLength;

  m_log << " initial macroscopic vorticity thickness: " << deltaOmega0 << endl;

  // width of the Gaussian distribution for perturbation scaling
  // depends on initial vorticity thickness
  sigma = deltaOmega0;

  for(MInt i = 0; i < a_noCells(); i++) {
    // set the mean velocity profile
    a_variable(i, PV->U) = m_Ma * LBCS * tanh(2.0 * a_coordinate(i, 1) / deltaOmega0);
    for(MInt dir = 1; dir < nDim; dir++) {
      a_variable(i, PV->VV[dir]) = 0.0;
    }
    a_variable(i, PV->RHO) = rho;

    GaussFactor = exp(-(a_coordinate(i, 1) / sigma) * (a_coordinate(i, 1) / sigma) * F1B2);

    a_variable(i, PV->U) += amp * F2 * (distrib(randNumGen) - F1B2) * GaussFactor;

    for(MInt dir = 1; dir < nDim; dir++) {
      a_variable(i, PV->VV[dir]) = amp * F2 * (distrib(randNumGen) - F1B2) * GaussFactor;
    }

    // the viscosity is saved for each cell
    // this is needed for subgrid modelling
    initNu(i, m_nu);
  }


  initNonEqDistFunctions();
}

/** \fn LbSolverDxQy::initLatticeBgkTurbulentPipe()
 *  \author Andreas Lintermann
 *  \date 24.01.2014
 *
 *  \propVal{initMethod,LB_TURBULENT_PIPE_INIT}
 *
 * Initializes standard Lattice BGK Turbulent Pipe.
 * Set property "initMethod" to "LB_TURBULENT_PIPE_INIT"
 *
 * \LBIC{LbSolverDxQy::initLatticeBgkTurbulentPipe(), initLatticeBgkTurbulentPipe, TurbulentPipeInit}
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbSolverDxQy<nDim, nDist, SysEqn>::initLatticeBgkTurbulentPipe() {
  TRACE();

  MFloat yPlus;
  MFloat amp = 0.1;
  MFloat uTau = (MFloat)m_ReTau * m_nu / m_referenceLength;
  MFloat radius = 0;
  ScratchSpace<MFloat> bBox(nDim * 2, AT_, "bBox");
  MFloat* bBoxPtr = &bBox[0];

  m_geometry->getBoundingBox(bBoxPtr);

  MFloat rho;

  // all cells have the same density
  if(m_densityFluctuations)
    rho = 0.0;
  else
    rho = 1.0;

  // since interpolated bounce back is used, this is the real half diameter
  const MFloat tmpWidth = 0.5 * fabs(bBox[1 + nDim] - bBox[1]);

  for(MInt i = 0; i < a_noCells(); i++) {
    for(MInt dir = 0; dir < nDim - 1; dir++) {
      radius += a_coordinate(i, dir) * a_coordinate(i, dir);
    }
    radius = sqrt(radius);

    yPlus = m_ReTau * (1.0 - radius / tmpWidth);

    if(yPlus <= 5.0) {
      a_variable(i, PV->VV[nDim - 1]) = uTau * yPlus;
    } else {
      if(yPlus <= 30 && yPlus > 5.0) {
        a_variable(i, PV->VV[nDim - 1]) = uTau * (C1 * log(yPlus) + C2);
      } else {
        if(yPlus > 30) {
          a_variable(i, PV->VV[nDim - 1]) = uTau * (C3 * log(yPlus) + C4);
        }
      }
    }

    for(MInt dir = 0; dir < nDim - 1; dir++) {
      a_variable(i, PV->VV[dir]) =
          amp * 2.0 * (0.5 - (MFloat)(1.0 * rand() / (RAND_MAX + 1.0))) * a_variable(i, PV->VV[nDim - 1]);
    }
    a_variable(i, PV->VV[nDim - 1]) +=
        amp * 2.0 * (0.5 - (MFloat)(1.0 * rand() / (RAND_MAX + 1.0))) * a_variable(i, PV->VV[nDim - 1]);
    a_variable(i, PV->RHO) = rho;

    // the viscosity is saved for each cell
    // this is needed for subgrid modelling
    initNu(i, m_nu);
  }
  initEqDistFunctions();
}

/** \fn LbSolverDxQy::initLatticeBgkVortex()
 *  \author Miro Gondrum
 *  \date 15.05.2021
 *
 *  \propVal{initMethod,LB_STEADY_VORTEX_INIT}
 *  \propVal{initMethod,LB_CONVECTING_VORTEX_INIT}
 *
 * Initializes standard Lattice BGK Vortex (at rest or convecting).
 * Set property "initMethod" to "LB_SPINNING_VORTICIES_INIT" or
 * "LB_CONVECTING_VORTEX_INIT" to use this function.
 *
 * Vortex is represented by vortex core model
 * A. Najafi-Yazdi, 2012. (http://dx.doi.org/10.1016/j.compfluid.2012.07.017)
 *
 * \LBIC{LbSolverDxQy::initLatticeBgkVortex(), initLatticeBgkVortex, VortexInit}
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbSolverDxQy<nDim, nDist, SysEqn>::initLatticeBgkVortex() {
  TRACE();

  MFloat coreRadius = 0.2;
  std::array<MFloat, 2> position;
  MBool eqInit = false;
  position.fill(0.0);

  coreRadius = Context::getSolverProperty<MFloat>("vortexCoreRadius", this->m_solverId, AT_, &coreRadius);
  const MFloat umax = Context::getSolverProperty<MFloat>("vortexMachMax", this->m_solverId, AT_, &m_Ma) * LBCS;
  eqInit = Context::getSolverProperty<MBool>("vortexEqInit", this->m_solverId, AT_, &eqInit);
  if(Context::propertyExists("vortexPosition", this->m_solverId)) {
    for(MInt d = 0; d < 2; d++) {
      position[d] = Context::getSolverProperty<MFloat>("vortexPosition", this->m_solverId, AT_, d);
    }
  }

  std::array<MFloat, 2> u_b = {0.0, 0.0};
  if(string2enum(m_initMethod) == LB_CONVECTING_VORTEX_INIT) {
    u_b[0] = m_Ma * LBCS;
  }

  if(domainId() == 0) {
    std::stringstream ss;
    ss << "Info: " << string2enum(m_initMethod) << std::endl;
    ss << "    vortexCoreRadius : " << coreRadius << std::endl;
    ss << "    vortexUmax       : " << umax << std::endl;
    ss << "    vortexPosition   : " << position[0] << ", " << position[1] << std::endl;
    std::cout << ss.str();
  }

  maia::parallelFor<false>(0, a_noCells(), [=](MInt cellId) {
    const MFloat x = a_coordinate(cellId, 0) - position[0];
    const MFloat y = a_coordinate(cellId, 1) - position[1];
    const MFloat r_rel = std::sqrt(POW2(x) + POW2(y)) / coreRadius;
    const MFloat theta = std::atan2(y, x);
    const MFloat factor = umax * std::exp(0.5);
    // Determine velocity from potential flow theory (x and y direction)
    MFloat up[2] = {0.0, 0.0};
    {
      const MFloat u_tan = factor * r_rel * std::exp(-0.5 * POW2(r_rel));
      up[0] = -std::sin(theta) * u_tan;
      up[1] = std::cos(theta) * u_tan;
    }
    const MFloat rho0 = 1.0;
    // const MFloat rho = rho0 * std::exp(-POW2(factor) / (2.0 * CSsq) * std::exp(1.0 - POW2(r_rel)));
    const MFloat rho = rho0 * std::exp(-0.5 * POW2(factor * F1BCS) * std::exp(-POW2(r_rel)));
    a_variable(cellId, PV->RHO) = rho;
    a_variable(cellId, PV->U) = (u_b[0] + up[0]);
    a_variable(cellId, PV->V) = (u_b[1] + up[1]);
    if constexpr(nDim == 3) a_variable(cellId, PV->W) = 0.0;
    initNu(cellId, m_nu);
  });
  if(eqInit) {
    initEqDistFunctions();
  } else {
    initNonEqDistFunctions();
  }
}

/** \fn LbSolverDxQy::initLatticeBgkSpinningVorticies()
 *  \author Miro Gondrum
 *  \date 23.06.2021
 *
 *  \propVal{initMethod,LB_SPINNING_VORTICIES_INIT}
 *
 * Initializes standard Lattice BGK spinning vorticies.
 * Set property "initMethod" to "LB_SPINNING_VORTICIES_INIT" to use this function.
 *
 * \LBIC{LbSolverDxQy::initLatticeBgkSpinningVorticies(), initLatticeBgkSpinningVorticies, SpinningVorticiesInit}
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbSolverDxQy<nDim, nDist, SysEqn>::initLatticeBgkSpinningVorticies() {
  TRACE();

  // normal zero init
  // initLatticeBgkLaminarDir(-1);

  // Fixed values for initial time and density
  const MFloat t = 0.0;
  // const MFloat gamma = Context::getSolverProperty<MFloat>("circulation", this->m_solverId, AT_);
  const MFloat r0 = Context::getSolverProperty<MFloat>("r0", this->m_solverId, AT_);
  const MFloat rC = Context::getSolverProperty<MFloat>("coreRadius", this->m_solverId, AT_);
  // Vatistas vortex core model (nC = 1 corresponds to Scully model)
  MInt nC = 1;
  nC = Context::getSolverProperty<MInt>("coreModelExponent", this->m_solverId, AT_, &nC);
  // Calculate rotation frequency and offsets
  const MFloat gamma = m_Ma * LBCS * 4.0 * PI * r0; // units !: [u_lb * l_stl]
  const MFloat dx = (this->c_cellLengthAtLevel(maxLevel()));
  const MFloat omega = m_Ma * LBCS / (r0 / dx);
  const MFloat bx = r0 * cos(omega * t);
  const MFloat by = r0 * sin(omega * t);

  const MFloat periode = 2.0 * PI / omega;
  m_log << "    Info spinning vorticities init: " << std::endl
        << "      * Iterations for one full rotation: " << periode << std::endl;
  const MFloat resolutionFactor = rC / dx;
  m_log << "      * coreRadius / dx_maxLevel = " << resolutionFactor << std::endl;


  // Set initial condition for all cells
  maia::parallelFor<false>(0, a_noCells(), [=](MInt cellId) {
    // Calculate vortex-local coordinates
    const MFloat x = a_coordinate(cellId, 0);
    const MFloat y = a_coordinate(cellId, 1);
    const MFloat rPos = std::sqrt(POW2(x - bx) + POW2(y - by));
    const MFloat thetaPos = std::atan2(y - by, x - bx);
    const MFloat rNeg = std::sqrt(POW2(x + bx) + POW2(y + by));
    const MFloat thetaNeg = std::atan2(y + by, x + bx);

    // Determine velocity from potential flow theory (x and y direction)
    MFloat u[2] = {0.0, 0.0};
    {
      const MFloat factorPos = (nC == 1) ? rPos / (rC * rC + rPos * rPos)
                                         : rPos / std::pow(std::pow(rC, 2 * nC) + std::pow(rPos, 2 * nC), 1.0 / nC);
      const MFloat factorNeg = (nC == 1) ? rNeg / (rC * rC + rNeg * rNeg)
                                         : rNeg / std::pow(std::pow(rC, 2 * nC) + std::pow(rNeg, 2 * nC), 1.0 / nC);
      u[0] = (factorPos * std::sin(thetaPos) + factorNeg * std::sin(thetaNeg)) * (-gamma / 2.0 / PI);
      u[1] = (factorPos * std::cos(thetaPos) + factorNeg * std::cos(thetaNeg)) * (gamma / 2.0 / PI);
    }

    // Set macroscopic variables
    // Determine pressure using steady-state Bernoulli equation
    const MFloat usqrB2 = 0.5 * (u[0] * u[0] + u[1] * u[1]);
    // const MFloat rho = 1.0 - usqrB2 / (CSsq + usqrB2);
    const MFloat rho = 1.0 / (1.0 + usqrB2 * F1BCSsq);
    a_variable(cellId, PV->RHO) = rho;
    a_variable(cellId, PV->U) = rho * u[0];
    a_variable(cellId, PV->V) = rho * u[1];
    a_oldVariable(cellId, PV->U) = a_variable(cellId, PV->U);
    a_oldVariable(cellId, PV->V) = a_variable(cellId, PV->V);
    if constexpr(nDim == 3) {
      a_variable(cellId, PV->W) = 0.0;
      a_oldVariable(cellId, PV->W) = 0.0;
    }
    initNu(cellId, m_nu);
  });

  initEqDistFunctions();
}

/** \fn LbSolverDxQy::initLatticeBgkLaminar()
 *  \author Andreas Lintermann
 *  \date 24.01.2011
 *
 *  \propVal{initMethod,LB_FROM_ZERO_INIT}
 *  \propVal{initMethod,LB_LAMINAR_INIT_[PM][XYZ]}
 *
 * Set property "initMethod" to "LB_FROM_ZERO_INIT" or to "LB_LAMINAR_INIT_[PM][XYZ]"
 * to use this function.
 * Simply calls initLatticeBgkLaminarDir(...) with a velocity direction or -1 for zero
 * velocity.
 *
 * \LBIC{LbSolverDxQy::initLatticeBgkLaminar(), initLatticeBgkLaminar, LaminarInit}
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbSolverDxQy<nDim, nDist, SysEqn>::initLatticeBgkLaminar() {
  TRACE();

  switch(string2enum(m_initMethod)) {
    case LB_FROM_ZERO_INIT:
      if constexpr(nDim == 2) {
        m_log << "    - parameters:    u=0.0, v=0.0, " << (m_densityFluctuations ? "rho=0.0" : "rho=1.0")
              << ((m_isThermal) ? ", T=1.0" : "") << ((m_isTransport) ? ", C=1.0" : "") << endl;
      } else {
        m_log << "    - parameters:    u=0.0, v=0.0, w=0, " << (m_densityFluctuations ? "rho=0.0" : "rho=1.0")
              << ((m_isThermal) ? ", T=1.0" : "") << ((m_isTransport) ? ", C=1.0" : "") << endl;
      }
      initLatticeBgkLaminarDir(-1);
      break;
    case LB_LAMINAR_INIT_PX:
      if constexpr(nDim == 2) {
        m_log << "    - parameters:    u=" << m_Ma * LBCS << ", v=0.0, "
              << (m_densityFluctuations ? "rho=0.0" : "rho=1.0") << ((m_isThermal) ? ", T=1.0" : "")
              << ((m_isTransport) ? ", C=1.0" : "") << endl;
      } else {
        m_log << "    - parameters:    u=" << m_Ma * LBCS << ", v=0.0, w=0, "
              << (m_densityFluctuations ? "rho=0.0" : "rho=1.0") << ((m_isThermal) ? ", T=1.0" : "")
              << ((m_isTransport) ? ", C=1.0" : "") << endl;
      }
      initLatticeBgkLaminarDir(0);
      break;
    case LB_LAMINAR_INIT_MX:
      if constexpr(nDim == 2) {
        m_log << "    - parameters:    u=-" << m_Ma * LBCS << ", v=0.0, "
              << (m_densityFluctuations ? "rho=0.0" : "rho=1.0") << ((m_isThermal) ? ", T=1.0" : "")
              << ((m_isTransport) ? ", C=1.0" : "") << endl;
      } else {
        m_log << "    - parameters:    u=-" << m_Ma * LBCS << ", v=0.0, w=0, "
              << (m_densityFluctuations ? "rho=0.0" : "rho=1.0") << ((m_isThermal) ? ", T=1.0" : "")
              << ((m_isTransport) ? ", C=1.0" : "") << endl;
      }
      initLatticeBgkLaminarDir(1);
      break;
    case LB_LAMINAR_INIT_PY:
      if constexpr(nDim == 2) {
        m_log << "    - parameters:    u=0.0, v=" << m_Ma * LBCS << ", "
              << (m_densityFluctuations ? "rho=0.0" : "rho=1.0") << ((m_isThermal) ? ", T=1.0" : "")
              << ((m_isTransport) ? ", C=1.0" : "") << endl;
      } else {
        m_log << "    - parameters:    u=0.0, v=" << m_Ma * LBCS << ", w=0.0, "
              << (m_densityFluctuations ? "rho=0.0" : "rho=1.0") << ((m_isThermal) ? ", T=1.0" : "")
              << ((m_isTransport) ? ", C=1.0" : "") << endl;
      }
      initLatticeBgkLaminarDir(2);
      break;
    case LB_LAMINAR_INIT_MY:
      if constexpr(nDim == 2) {
        m_log << "    - parameters:    u=0.0, v=-" << m_Ma * LBCS << ", "
              << (m_densityFluctuations ? "rho=0.0" : "rho=1.0") << ((m_isThermal) ? ", T=1.0" : "")
              << ((m_isTransport) ? ", C=1.0" : "") << endl;
      } else {
        m_log << "    - parameters:    u=0.0, v=-" << m_Ma * LBCS << ", w=0.0, "
              << (m_densityFluctuations ? "rho=0.0" : "rho=1.0") << ((m_isThermal) ? ", T=1.0" : "")
              << ((m_isTransport) ? ", C=1.0" : "") << endl;
      }
      initLatticeBgkLaminarDir(3);
      break;
    case LB_LAMINAR_INIT_PZ:
      if constexpr(nDim == 2) {
        m_log << "    - parameters:    u=0.0, v=0.0, " << (m_densityFluctuations ? "rho=0.0" : "rho=1.0")
              << ((m_isThermal) ? ", T=1.0" : "") << ((m_isTransport) ? ", C=1.0" : "") << endl;
        initLatticeBgkLaminarDir(-1);
      } else {
        m_log << "    - parameters:    u=0.0, v=0.0, w=" << m_Ma * LBCS << ", "
              << (m_densityFluctuations ? "rho=0.0" : "rho=1.0") << ((m_isThermal) ? ", T=1.0" : "")
              << ((m_isTransport) ? ", C=1.0" : "") << endl;
        initLatticeBgkLaminarDir(4);
      }
      break;
    case LB_LAMINAR_INIT_MZ:
      if constexpr(nDim == 2) {
        m_log << "    - parameters:    u=0.0, v=0.0, " << (m_densityFluctuations ? "rho=0.0" : "rho=1.0")
              << ((m_isThermal) ? ", T=1.0" : "") << ((m_isTransport) ? ", C=1.0" : "") << endl;
        initLatticeBgkLaminarDir(-1);
      } else {
        m_log << "    - parameters:    u=0.0, v=0.0, w=-" << m_Ma * LBCS << ", "
              << (m_densityFluctuations ? "rho=0.0" : "rho=1.0") << ((m_isThermal) ? ", T=1.0" : "")
              << ((m_isTransport) ? ", C=1.0" : "") << endl;
        initLatticeBgkLaminarDir(5);
      }
      break;
    default:
      break;
  }
}

/** \fn LbSolverDxQy::initLatticeBgkGaussPulse()
 *  \author Miro Gondrum
 *  \date 17.03.2020
 *
 *  \propVal{initMethod,LB_GAUSS_PULSE_INIT}
 *
 * Initializes a Gaussian pulse.
 * Set property "initMethod" to "LB_GAUSS_PULSE_INIT" to use this function.
 *
 * \LBIC{LbSolverDxQy::initLatticeBgkGaussPulse(), initLatticeBgkGaussPulse, GaussPulseInit}
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbSolverDxQy<nDim, nDist, SysEqn>::initLatticeBgkGaussPulse() {
  TRACE();

  const MFloat rho0 = (m_densityFluctuations) ? 0.0 : 1.0;
  const MFloat origin[3] = {F0, F0, F0};
  const MFloat amp = 1e-8;
  const MFloat radius = 0.01;

  initLatticeBgkLaminarDir(-1); // normal zero init
  maia::parallelFor<false>(0, a_noCells(), [&](MInt i) {
    MFloat r[nDim];
    for(MInt j = 0; j < nDim; j++)
      r[j] = a_coordinate(i, j) - origin[j];
    const MFloat distance = std::sqrt(std::inner_product(&r[0], &r[nDim], &r[0], .0));
    const MFloat rhoFluct = amp * std::exp(-(distance * distance) / radius);
    a_variable(i, PV->RHO) = rho0 + rhoFluct;
    a_oldVariable(i, PV->RHO) = a_variable(i, PV->RHO);
  });
  initEqDistFunctions();
}

// Initialize a Gaussian hill for diffusion processes. The velocity set uniform to zero.
// Set property "initMethod" to "LB_GAUSS_HILL_DIFFUSION" to use this function.
// Note this has to be used together with the "solverMethod" = "MAIA_LATTICE_BGK_TRANSPORT"
// as well as that this can only used in "2D".
template <MInt nDim, MInt nDist, class SysEqn>
void LbSolverDxQy<nDim, nDist, SysEqn>::initLatticeBgkGaussDiffusion() {
  TRACE();

  const MFloat origin[2] = {F0, F0}; // center of the hill
  const MFloat width = 20.0 * this->c_cellLengthAtLevel(this->maxLevel());
  const MFloat maxConcentration = this->m_initCon;
  const MFloat maxTemp = this->m_initTemperatureKelvin;

  maia::parallelFor<false>(0, a_noCells(), [=](MInt i) {
    // Initialize zero velocity
    a_variable(i, PV->U) = F0;
    a_variable(i, PV->V) = F0;
    a_oldVariable(i, PV->U) = F0;
    a_oldVariable(i, PV->V) = F0;

    // Inizialize density
    MFloat rho = 1.0;
    a_variable(i, PV->RHO) = rho;
    a_oldVariable(i, PV->RHO) = rho;

    // Initialize a Gaussian hill for the conzentration.
    if((m_isTransport) && nDim == 2) {
      MFloat r[nDim];
      for(MInt d = 0; d < nDim; d++)
        r[d] = a_coordinate(i, d) - origin[d];
      const MFloat distSq = pow(r[0], 2) + pow(r[1], 2);
      const MFloat C = maxConcentration * exp(-(distSq) / (2 * pow(width, 2)));
      a_variable(i, PV->C) = C;
      a_oldVariable(i, PV->C) = C;
    }

    if((m_isThermal) && nDim == 2) {
      MFloat r[nDim];
      for(MInt d = 0; d < nDim; d++)
        r[d] = a_coordinate(i, d) - origin[d];
      const MFloat distSq = pow(r[0], 2) + pow(r[1], 2);
      const MFloat T = maxTemp * exp(-(distSq) / (2 * pow(width, 2)));
      a_variable(i, PV->T) = T;
      a_oldVariable(i, PV->T) = T;
    }

    // the viscosity is saved for each cell
    // this is needed for subgrid modelling
    initNu(i, m_nu);
  });

  // Compute the equilibrium distributions
  initEqDistFunctions();

  if(m_isThermal) {
    initThermalEqDistFunctions();
  }
  if(m_isTransport) {
    initTransportEqDistFunctions();
  }
}

// Initialize a Gaussian hill for diffusion processes. The velocity set uniform to zero.
// Set property "initMethod" to "LB_GAUSS_HILL_ADVECTION" to use this function.
// Note this has to be used together with the "solverMethod" = "MAIA_LATTICE_BGK_TRANSPORT"
// as well as that this can only used in "2D".
template <MInt nDim, MInt nDist, class SysEqn>
void LbSolverDxQy<nDim, nDist, SysEqn>::initLatticeBgkGaussAdvection() {
  TRACE();

  const MFloat origin[2] = {F0, F0}; // center of the hill
  const MFloat width = 20.0 * this->c_cellLengthAtLevel(this->maxLevel());
  const MFloat maxConcentration = this->m_initCon;
  const MFloat maxTemp = this->m_initTemperatureKelvin;

  maia::parallelFor<false>(0, a_noCells(), [&](MInt i) {
    // Initialize zero velocity
    MFloat velocity = m_Ma * LBCS;
    a_variable(i, PV->U) = velocity;
    a_variable(i, PV->V) = velocity;
    a_oldVariable(i, PV->U) = velocity;
    a_oldVariable(i, PV->V) = velocity;

    // Inizialize density
    MFloat rho = 1.0;
    a_variable(i, PV->RHO) = rho;
    a_oldVariable(i, PV->RHO) = rho;

    // Initialize a Gaussian hill for the conzentration.
    if((m_isTransport) && nDim == 2) {
      MFloat r[nDim];
      for(MInt d = 0; d < nDim; d++)
        r[d] = a_coordinate(i, d) - origin[d];
      const MFloat distSq = pow(r[0], 2) + pow(r[1], 2);
      const MFloat C = maxConcentration * exp(-(distSq) / (2 * pow(width, 2)));
      a_variable(i, PV->C) = C;
      a_oldVariable(i, PV->C) = C;
    }

    if((m_isThermal) && nDim == 2) {
      MFloat r[nDim];
      for(MInt d = 0; d < nDim; d++)
        r[d] = a_coordinate(i, d) - origin[d];
      const MFloat distSq = pow(r[0], 2) + pow(r[1], 2);
      const MFloat T = maxTemp * exp(-(distSq) / (2 * pow(width, 2)));
      a_variable(i, PV->T) = T;
      a_oldVariable(i, PV->T) = T;
    }

    // the viscosity is saved for each cell
    // this is needed for subgrid modelling
    initNu(i, m_nu);
  });

  // Compute the equilibrium distributions
  initEqDistFunctions();

  if(m_isThermal) {
    initThermalEqDistFunctions();
  }
  if(m_isTransport) {
    initTransportEqDistFunctions();
  }
}

/** \fn LbSolverDxQy::initLatticeBgkLaminarCylinder()
 *  \author Miro Gondrum
 *  \date 31.03.2022
 *
 *  \propVal{initMethod,LB_LAMINAR_CYLINDER_INIT}
 *
 * Initializes a laminar cylinder with perturbation.
 * Set property "initMethod" to "LB_LAMINAR_CYLINDER_INIT" to use this function.
 * Perturbation forces a shorter transient.
 *
 * \LBIC{LbSolverDxQy::initLatticeBgkLaminarCylinder(), initLatticeBgkLaminarCylinder, LaminarCylinderInit}
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbSolverDxQy<nDim, nDist, SysEqn>::initLatticeBgkLaminarCylinder() {
  // normal init in pos +x direction
  initLatticeBgkLaminarDir(0);
  // now add perturbation to force faster transient
  constexpr MFloat factor = 0.1;
  const MFloat bb = 2.0; // bounding box where perturbation is applied
  maia::parallelFor<false>(0, a_noCells(), [=](MInt i) {
    const MFloat x = a_coordinate(i, 0);
    const MFloat y = a_coordinate(i, 1);
    // (x,y) in bb
    if(-bb < x && x < bb && -bb < y && y < bb) {
      const MFloat pertFactor = (1.0 + factor * std::sin(y / 2.0 * PI));
      a_variable(i, PV->U) *= pertFactor;
    }
  });
  initEqDistFunctions();
}

/** \brief Initializes standard Lattice BGK laminar with or without a direction
 *
 * \author Andreas Lintermann
 * \date 10.04.2017
 *
 * The direction parameter is given by:
 * - -1 : set zero velocity
 * - 0/1: positive / negative x
 * - 2/3: positive / negative y
 * - 4/5: positive / negative z
 *
 * \param[in] dir the direction to set
 *
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbSolverDxQy<nDim, nDist, SysEqn>::initLatticeBgkLaminarDir(MInt dir) {
  TRACE();

  MFloat value = m_Ma * LBCS;
  std::array<MInt, nDim - 1> o_dirs{};
  o_dirs[0] = PV->V;
  if constexpr(nDim == 3) o_dirs[1] = PV->W;

  if(dir < 0) {
    dir = PV->U;
    value = 0.0;
  } else {
    if(dir % 2) value *= -1.0;

    switch(dir) {
      case 0:
      case 1:
        dir = PV->U;
        break;
      case 2:
      case 3:
        dir = PV->V;
        o_dirs[0] = PV->U;
        break;
      case 4:
      case 5:
        dir = PV->W;
        o_dirs[0] = PV->U;
        if constexpr(nDim == 3) o_dirs[1] = PV->V;
        break;
      default:
        break;
    }
  }

  maia::parallelFor<false>(0, a_noCells(), [=](MInt i) {
    a_variable(i, dir) = value;
    a_oldVariable(i, dir) = value;
    for(MInt j = 0; j < nDim - 1; j++) {
      a_variable(i, o_dirs[j]) = 0.0;
      a_oldVariable(i, o_dirs[j]) = 0.0;
    }

    MFloat rho = m_densityFluctuations ? 0.0 : 1.0;
    if(m_initDensityGradient) {
      for(MInt j = 0; j < nDim; j++) {
        rho += a_coordinate(i, j) / this->c_cellLengthAtLevel(maxLevel()) * m_volumeAccel[j] * F1BCSsq;
      }
    }
    a_variable(i, PV->RHO) = rho;
    a_oldVariable(i, PV->RHO) = rho;

    if(m_isThermal) {
      a_kappa(i) = m_kappa;
      a_variable(i, PV->T) = 1.0;
      a_oldVariable(i, PV->T) = 1.0;
    }
    if(m_isTransport) {
      a_variable(i, PV->C) = 0.0;
      a_oldVariable(i, PV->C) = 0.0;
    }

    // the viscosity is saved for each cell
    // this is needed for subgrid modelling
    initNu(i, m_nu);
  });

  initEqDistFunctions();

  if(m_isThermal) {
    initThermalEqDistFunctions();
  }
  if(m_isTransport) {
    initTransportEqDistFunctions();
  }
}

/** \brief Calculates equilibrium distribution functions after initialization of macroscopic variables.
 * \author Andreas Lintermann
 * \date 24.01.2011
 *
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbSolverDxQy<nDim, nDist, SysEqn>::initEqDistFunctions() {
  TRACE();

  for(MInt i = 0; i < a_noCells(); i++) {
    const MFloat rho = a_variable(i, PV->RHO);
    std::array<MFloat, nDim> u;
    for(MInt d = 0; d < nDim; d++) {
      u[d] = a_variable(i, d);
    }
    setEqDists(i, rho, u.data());
  }
}

template <MInt nDim, MInt nDist, class SysEqn>
void LbSolverDxQy<nDim, nDist, SysEqn>::initNonEqDistFunctions() {
  TRACE();
  // first init equilibrium everywhere
  initEqDistFunctions();

  // now add non-equilibrium part
  maia::parallelFor<false>(0, a_noCells(), [&](MInt i) {
    // calculate spatial derivatives of velocity vector in the current cell center
    MFloat c[nDim][nDim]; // dv/dz = c[1][2]
    this->calculateVelocityDerivative(i, c);

    MFloat trace = 0.0;
    for(MInt d = 0; d < nDim; d++) {
      trace += c[d][d];
    }

    // add non-eq parts
    // cf. Eq (5.88), Krueger et al., 2017, https://doi.org/10.1007/978-3-319-44649-3 .
    const MInt lvlDiff = maxLevel() - a_level(i);
    const MFloat tau = F1B2 + F1BCSsq * a_nu(i) * FFPOW2(lvlDiff); // tau_0
    const MFloat tauRho = tau * a_variable(i, PV->RHO);
    for(MInt dist = 0; dist < nDist - 1; dist++) {
      const MFloat t = Ld::tp(Ld::distType(dist));
      a_distribution(i, dist) += t * tauRho * trace;
      for(MInt k = 0; k < nDim; k++) {
        for(MInt l = 0; l < nDim; l++) {
          a_distribution(i, dist) -=
              t * tauRho * F1BCSsq * (Ld::idFld(dist, k) - 1) * (Ld::idFld(dist, l) - 1) * c[l][k];
        }
      }
    }
    a_distribution(i, Ld::lastId()) += Ld::tp(0) * tauRho * trace;

    // oldDist = dist
    for(MInt j = 0; j < Ld::lastId() + 1; j++) {
      a_oldDistribution(i, j) = a_distribution(i, j);
    }
  });
}

/** \brief Calculates equilibrium distribution functions after initialization of macroscopic variables.
 * \author Andreas Lintermann
 * \date 23.02.2011
 *
 * This is pretty much the same function as initEqDistFunctions(), except that the thermal macroscopic variable \f$T\f$
 * (temperature) is initilized as well with \f$T = 1.0\f$. After that, the equilibrium distribution functions are
 * calculated for all the LBGK- and the TLBGK-case.
 *
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbSolverDxQy<nDim, nDist, SysEqn>::initThermalEqDistFunctions() {
  TRACE();

  for(MInt i = 0; i < a_noCells(); i++) {
    const MFloat rho = a_variable(i, PV->RHO);
    const MFloat t = a_variable(i, PV->T);
    std::array<MFloat, nDim> u;
    for(MInt d = 0; d < nDim; d++) {
      u[d] = a_variable(i, d);
    }

    setEqDistsThermal(i, t, rho, u.data());
  }
}

/** \brief Calculates equilibrium distribution functions after initialization of macroscopic variables.
 * \author Shota Ito
 * \date 07.06.2022
 *
 * This is pretty much the same function as initEqDistFunctions(), except that the transport macroscopic variable
 * \f$C\f$ (concentration) is initilized as well with \f$C = 1.0\f$. After that, the equilibrium distribution functions
 * are calculated for all the LBGK- and the Transport-LBGK-case.
 *
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbSolverDxQy<nDim, nDist, SysEqn>::initTransportEqDistFunctions() {
  TRACE();

  for(MInt i = 0; i < a_noCells(); i++) {
    const MFloat c = a_variable(i, PV->C);
    std::array<MFloat, nDim> u;
    for(MInt d = 0; d < nDim; d++) {
      u[d] = a_variable(i, d);
    }

    setEqDistsTransport(i, c, u.data());
  }
}

/** \brief Calculates residuals and prints to file
 * \author Andreas Lintermann
 * \date 13.08.2012
 *
 * Calculates the local max. and averaged max. residual for this process. The maximum over all processes is
 * found by using a MPI_Allreduce with the MAX_LOC option. Only the processor with the max. of all values
 * writes to disk.
 * The averaged max. resiudal is found by averaging over the number of cells per process and then summing up
 * and collecting the information on proc. 0 by using MPI_Reduce with the SUM option. Then this residual is
 * additionally averaged by the number of processes.
 *
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbSolverDxQy<nDim, nDist, SysEqn>::calculateResidual() {
  TRACE();

  if(domainId() == 0) mRes.open(m_resFileName, std::ofstream::app);

  constexpr MInt noVars = 1 + nDim;
  for(MInt v = 0; v < noVars; v++) {
    m_residual[v] = F0;
    m_tmpResidual[v] = F0;
    m_tmpResidualLvl[v] = 0;
    m_maxResId[v] = -1;
  }

  MInt count = 0;
  MInt sumCount = 0;

  for(MInt id = 0; id < m_currentMaxNoCells; id++) {
    MInt l_id = m_activeCellList[id];

    // skip halo cells
    if(l_id >= this->grid().noInternalCells()) continue;

    const MInt pCellId = l_id;

    // skip interfacechildren
    if(a_isInterfaceChild(pCellId)) continue;

    MFloat diff[noVars] = {F0};
    const MFloat* const pcoordinates = &a_coordinate(pCellId, 0);
    MInt level = this->a_level(pCellId);

    for(MInt v = 0; v < noVars; v++) {
      diff[v] = fabs(a_variable(pCellId, v) - a_oldVariable(pCellId, v));
      m_residual[v] += diff[v];
      if(m_tmpResidual[v] * m_tmpResidual[v] <= diff[v] * diff[v]) {
        m_tmpResidual[v] = diff[v];
        m_maxResId[v] = l_id;
        for(MInt d = 0; d < nDim; d++) {
          m_rescoordinates[v][d] = pcoordinates[d];
        }
        m_tmpResidualLvl[v] = level;
      }
    }
    count++;
  }
  cerr.precision(7);

  stringstream res_t;
  res_t << globalTimeStep;
  MString res_t_f = res_t.str();

  stringstream res_Re;
  res_Re << (m_referenceLength * m_Ma * LBCS / m_nu);
  MString res_Re_f = res_Re.str();

  struct {
    MFloat val;
    MInt rank;
  } sendBufRes[noVars], rcvBufRes[noVars];

  MFloat sendBufResAvg[noVars] = {0.0};
  for(MInt v = 0; v < noVars; v++) {
    sendBufResAvg[v] = m_residual[v];
    sendBufRes[v].val = m_tmpResidual[v];
  }
  MFloat rcvBufResAvg[noVars];

  for(MInt i = 0; i < noVars; i++)
    sendBufRes[i].rank = domainId();

  MString strs[noVars] = {};
  strs[0] = "U";
  strs[1] = "V";
  if constexpr(nDim == 3) strs[2] = "W";
  strs[nDim] = "RHO";

  MPI_Allreduce(sendBufRes, rcvBufRes, noVars, MPI_DOUBLE_INT, MPI_MAXLOC, mpiComm(), AT_, "sendBufRes", "rcvBufRes");
  MPI_Reduce(sendBufResAvg, rcvBufResAvg, noVars, MPI_DOUBLE, MPI_SUM, 0, mpiComm(), AT_, "sendBufResAvg",
             "rcvBufResAvg");

  // To calculate the avagare residual, get total number of cells without the halo cells
  if(noDomains() == 1) {
    // When serial, then just use the counted number of cells of the root process
    sumCount = count;
  } else {
    // When parallel, then sum up all counted number of cells of each process
    MPI_Reduce(&count, &sumCount, 1, MPI_INT, MPI_SUM, 0, mpiComm(), AT_, "count", "sumCount");
  }

  // all max. averaged residuals
  if(domainId() == 0) {
    mRes << std::endl;
    mRes << "#------------------------------" << std::endl;
    mRes << "In timestep: " << res_t_f << ", with Re: " << res_Re_f << ";" << std::endl;
    mRes << "\n";
    for(MInt i = 0; i < noVars; i++)
      mRes << "   Max. averaged residual found for " << strs[i] << " is " << (rcvBufResAvg[i] / sumCount) << "\n";
  }

  if(noDomains() > 1) {
    for(MInt i = 0; i < noVars; i++) {
      MInt noElements2Snd = 0;
      std::array<MFloat, nDim> sndBufFloat;
      std::array<MInt, 2> sndBufInt;
      std::array<MFloat, nDim> recvElementsF;
      std::array<MInt, 2> recvElementsI;
      MIntScratchSpace noElements2SndPerCPU(noDomains(), AT_, "noElements2SndPerCPU");
      MIntScratchSpace posElements2SndPerCPU(noDomains(), AT_, "posElements2SndPerCPU");

      for(MInt dom = 0; dom < noDomains(); dom++) {
        noElements2SndPerCPU[dom] = 0;
        posElements2SndPerCPU[dom] = 0;
        if(rcvBufRes[i].rank == dom) {
          noElements2SndPerCPU[dom] = nDim;
        }
        if(dom > 0) {
          posElements2SndPerCPU[dom] = posElements2SndPerCPU[dom - 1] + noElements2SndPerCPU[dom - 1];
        }
      }

      if(rcvBufRes[i].rank == domainId() && m_maxResId[i] != -1) {
        noElements2Snd = nDim;
        for(MInt d = 0; d < nDim; d++) {
          sndBufFloat[d] = m_rescoordinates[i][d];
        }
      }

      MPI_Gatherv(sndBufFloat.data(), noElements2Snd, MPI_DOUBLE, recvElementsF.data(), noElements2SndPerCPU.data(),
                  posElements2SndPerCPU.data(), MPI_DOUBLE, 0, mpiComm(), AT_, "sndBufFloat.data()",
                  "recvElementsF.data()");

      for(MInt dom = 0; dom < noDomains(); dom++) {
        noElements2SndPerCPU[dom] = 0;
        posElements2SndPerCPU[dom] = 0;
        if(rcvBufRes[i].rank == dom) {
          noElements2SndPerCPU[dom] = 2;
        }
        if(dom > 0) {
          posElements2SndPerCPU[dom] = posElements2SndPerCPU[dom - 1] + noElements2SndPerCPU[dom - 1];
        }
      }

      if(rcvBufRes[i].rank == domainId() && m_maxResId[i] != -1) {
        noElements2Snd = 2;
        sndBufInt[0] = m_tmpResidualLvl[i];
        sndBufInt[1] = m_maxResId[i];
      }

      MPI_Gatherv(&sndBufInt, noElements2Snd, MPI_INT, recvElementsI.data(), noElements2SndPerCPU.data(),
                  posElements2SndPerCPU.data(), MPI_INT, 0, mpiComm(), AT_, "sndBufInt", "recvElementsI.data()");

      // all max. local residuals
      if(domainId() == 0) {
        mRes << "\n   Max. local residual found for " << strs[i] << " on proc.: " << rcvBufRes[i].rank << " with "
             << rcvBufRes[i].val << "\n"
             << "      Level: " << recvElementsI[0] << "\n"
             << "      ID: " << recvElementsI[1] << "\n"
             << "      Coordinates: (";
        for(MInt d = 0; d < nDim - 1; d++) {
          mRes << recvElementsF[d] << ", ";
        }
        mRes << recvElementsF[nDim - 1] << ")" << endl;
      }
    }
  } else {
    // Print residual for the only process that is running
    // No communication required
    for(MInt i = 0; i < noVars; i++) {
      mRes << "\n   Max. local residual found for " << strs[i] << " on proc.: " << domainId() << " with "
           << rcvBufRes[i].val << "\n"
           << "      Level: " << m_tmpResidualLvl[i] << "\n"
           << "      ID: " << m_maxResId[i] << "\n"
           << "      Coordinates: (";
      for(MInt d = 0; d < nDim - 1; d++) {
        mRes << m_rescoordinates[i][d] << ", ";
      }
      mRes << m_rescoordinates[i][nDim - 1] << ")" << endl;
    }
  }

  // Close Residual file
  if(domainId() == 0) {
    mRes << std::endl;
    mRes.close();
  }
}

// \brief transfer from coarse to fine
template <MInt nDim, MInt nDist, class SysEqn>
void LbSolverDxQy<nDim, nDist, SysEqn>::prolongation() {
  TRACE();
  m_interface->prolongation();
}

// \brief Transfer from fine to coarse
template <MInt nDim, MInt nDist, class SysEqn>
void LbSolverDxQy<nDim, nDist, SysEqn>::restriction() {
  TRACE();
  m_interface->restriction();
}

/**
 * \brief:     Initialize parent variables from children
 * \details:   Call interface function that performs initialization
 * \param[in]: parentId: Solver cell id of parent that will be coarsen
 * \author: Philipp Brokof, Moritz Waldmann
 **/
template <MInt nDim, MInt nDist, class SysEqn>
void LbSolverDxQy<nDim, nDist, SysEqn>::removeChildsLb(const MInt parentId) {
  m_interface->removeChildren(parentId);
}

/**
 * \brief:     Initialize child variables from parent
 * \details:   Call interface function that performs initialization
 * \param[in]: parentId: Solver id of cell that is refined
 *            childIds: Solver ids of new finde cells
 * \author: Philipp Brokof, Moritz Waldmann
 **/
template <MInt nDim, MInt nDist, class SysEqn>
void LbSolverDxQy<nDim, nDist, SysEqn>::refineCellLb(const MInt parentId, const MInt* childIds) {
  m_interface->refineCell(parentId, childIds);
}

/** \page sensorsLB
 *
 *  \section interface Interface
 *
 *  This sensors ensures a band of refined cells at static walls.<br>
 *  Property: <code>INTERFACE</code>
 */
/**
 * \brief Simple boundary sensor for ensuring boundary refinement during adaptation
 *
 * \author: Julian Vorspohl <j.vorspohl@aia.rwth-aachen.de>
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbSolverDxQy<nDim, nDist, SysEqn>::sensorInterface(std::vector<std::vector<MFloat>>& sensors,
                                                        std::vector<std::bitset<64>>& sensorCellFlag,
                                                        std::vector<MFloat>& sensorWeight, MInt sensorOffset,
                                                        MInt sen) {
  m_log << "   - Sensor preparation for the boundary sensor for " << m_maxNoSets << " starting from " << m_levelSetId
        << endl;
  cerr0 << "   - Sensor preparation for the boundary sensor for " << m_maxNoSets << " starting from " << m_levelSetId
        << endl;

  // inList > 0: Cells are marked for refinement
  // inList = 1: Seeding cells for markSurroundingCells
  MIntScratchSpace inList(a_noCells(), AT_, "inList");
  inList.fill(0);

  // Fill inList with cells wich are to be refined
  MInt interfaceCellCounter = 0;
  for(MInt cellId = 0; cellId < a_noCells(); cellId++) {
    if(a_isHalo(cellId)) continue;
    if(inList[cellId] > 0) continue; // Already marked for refinement

    if(this->a_isBndryCell(cellId)) {
      if(a_level(cellId) <= this->m_maxSensorRefinementLevel[sen]) {
        inList[cellId] = 1;
        interfaceCellCounter++;

        // Propagate information to lower levels...
        MInt currentCellId = cellId;
        for(MInt parentLevel = this->a_level(cellId); parentLevel >= minLevel(); parentLevel--) {
          const MInt parentCellId = this->c_parentId(currentCellId);
          if(parentCellId == -1) {
            break;
          }
          inList[parentCellId] = 1;
          interfaceCellCounter++;
          currentCellId = parentCellId;
        }
      }
    }
  }

  std::cout << "Found " << interfaceCellCounter << " static boundary cells..." << std::endl;

  // Fill inList with cells wich are to be refined
  interfaceCellCounter = 0;
  MInt startSet = m_levelSetId;
  for(MInt cellId = 0; cellId < a_noCells(); cellId++) {
    if(a_isHalo(cellId)) continue;
    for(MInt set = startSet; set < m_maxNoSets; set++) {
      if(inList[cellId] > 0) continue; // Already marked for refinement
      if(fabs(a_levelSetFunctionMB(cellId, set)) < this->c_cellLengthAtCell(cellId)) {
        if(a_level(cellId) < this->m_maxSensorRefinementLevel[sen]) {
          inList[cellId] = 1;
          interfaceCellCounter++;
        }
      } else {
        for(MInt dir = 0; dir < 2 * nDim; dir++) {
          if(a_hasNeighbor(cellId, dir) > 0) {
            MInt nghbrId = c_neighborId(cellId, dir);
            if((a_levelSetFunctionMB(nghbrId, set) * a_levelSetFunctionMB(cellId, set) < F0)) {
              if(a_level(cellId) < this->m_maxSensorRefinementLevel[sen]) {
                inList[cellId] = 1;
                interfaceCellCounter++;
              }
              break;
            }
          }
        }
      }
    }
  }

  std::cout << "Found " << interfaceCellCounter << " moving boundary cells..." << std::endl;

  // Mark surrounding cells according to m_bandWidth
  for(MInt level = minLevel(); level < this->m_maxSensorRefinementLevel[sen]; level++) {
    this->exchangeData(inList.data());
    this->markSurrndCells(inList, m_bandWidth[level], level, this->m_refineDiagonals);
  }

  // Transfer inList to sensor
  for(MInt cellId = 0; cellId < noInternalCells(); cellId++) {
    ASSERT(!a_isHalo(cellId), "");
    if(inList(cellId) == 0) {
      if(a_level(cellId) == minLevel()) continue;
      if(inList(c_parentId(cellId))) continue;
      if(c_isLeafCell(cellId)) {
        const MInt gridCellId = grid().tree().solver2grid(cellId);
        sensors[sensorOffset + sen][gridCellId] = -1.0;
        sensorCellFlag[gridCellId][sensorOffset + sen] = true;
      }
    } else {
      if(a_level(cellId) < this->m_maxSensorRefinementLevel[sen]) {
        if(c_noChildren(cellId) > 0) continue;
        const MInt gridCellId = grid().tree().solver2grid(cellId);
        sensors[sensorOffset + sen][gridCellId] = 1.0;
        sensorCellFlag[gridCellId][sensorOffset + sen] = true;
      }
    }
  }
  sensorWeight[sensorOffset + sen] = this->m_sensorWeight[sen];
}

/** \page sensorsLB
 *
 *  \section vorticity Vorticity
 *
 *  Property: <code>VORTICITY</code>
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbSolverDxQy<nDim, nDist, SysEqn>::sensorVorticity(std::vector<std::vector<MFloat>>& sensors,
                                                        std::vector<std::bitset<64>>& sensorCellFlag,
                                                        std::vector<MFloat>& sensorWeight, MInt sensorOffset,
                                                        MInt sen) {
  m_log << "   - Sensor preparation for the vorticity sensor" << endl;
  cerr0 << "   - Sensor preparation for the vorticity sensor" << endl;

  MFloat dx = 0.0;

  std::array<MFloat, nDim> u{};
  std::array<std::array<MFloat, nDim>, nDim> u_p{};
  std::array<std::array<MFloat, nDim>, nDim> u_m{};

  MFloat rho = 0.0;
  MFloat rho_p[nDim] = {0.0};
  MFloat rho_m[nDim] = {0.0};

  MFloat dudy = 0.0;
  MFloat dvdx = 0.0;

  MFloat dwdx = 0.0;
  MFloat dwdy = 0.0;

  MFloat phi_v = 0.0; // Sensor value vorticity

  sensorWeight[sensorOffset + sen] = this->m_sensorWeight[sen];

  // Set cell seize weihts
  // TODO labels:LB Make available by property?
  MFloat cellSizeWeight_v = 1.5;

  for(MInt cellId = 0; cellId < noInternalCells(); cellId++) {
    ASSERT(!a_isHalo(cellId), "");
    const MInt gridId = this->grid().tree().solver2grid(cellId);

    // Mesh adaptation is only performed on active cells
    // --> untag inactive cells
    if(!a_isActive(cellId)) {
      sensorCellFlag[gridId][sensorOffset + sen] = 0;
      continue;
    }

    // Get cell size for cell size weighting
    dx = grid().cellLengthAtLevel(c_level(cellId));

    // Compute macroscopic variables for sensor evaluation
    // All cells performed propagation step before. Thus, in-coming distributions
    // of all cells are in the same time step and comparable
    // --> Compute macroscopic variables from in-coming distributions.
    //
    // Depending on the neighbors of the cell, variables for
    // central, backwards or forwards differentials are computed
    calculateMacroscopicVariables(cellId, rho, u.data());
    if(a_hasNeighbor(cellId, 0) && a_isActive(c_neighborId(cellId, 0))) {
      MInt nghbrId = c_neighborId(cellId, 0);
      calculateMacroscopicVariables(nghbrId, rho_m[0], u_m[0].data());
    }
    if(a_hasNeighbor(cellId, 1) && a_isActive(c_neighborId(cellId, 1))) {
      MInt nghbrId = c_neighborId(cellId, 1);
      calculateMacroscopicVariables(nghbrId, rho_p[0], u_p[0].data());
    }
    if(a_hasNeighbor(cellId, 2) && a_isActive(c_neighborId(cellId, 2))) {
      MInt nghbrId = c_neighborId(cellId, 2);
      calculateMacroscopicVariables(nghbrId, rho_m[1], u_m[1].data());
    }
    if(a_hasNeighbor(cellId, 3) && a_isActive(c_neighborId(cellId, 3))) {
      MInt nghbrId = c_neighborId(cellId, 3);
      calculateMacroscopicVariables(nghbrId, rho_p[1], u_p[1].data());
    }
    if constexpr(nDim == 3) {
      if(a_hasNeighbor(cellId, 4) && a_isActive(c_neighborId(cellId, 4))) {
        MInt nghbrId = c_neighborId(cellId, 4);
        calculateMacroscopicVariables(nghbrId, rho_m[2], u_m[2].data());
      }
      if(a_hasNeighbor(cellId, 5) && a_isActive(c_neighborId(cellId, 5))) {
        MInt nghbrId = c_neighborId(cellId, 5);
        calculateMacroscopicVariables(nghbrId, rho_p[2], u_p[2].data());
      }
    }

    // (1) Difference in total pressure
    // Incident flow total pressure depends on bc!
    // Here so far only for bc10000 --> Assume rho=1 on boundary.
    // Central differences in x-direction
    if((a_hasNeighbor(cellId, 0) && a_isActive(c_neighborId(cellId, 0)))
       && (a_hasNeighbor(cellId, 1) && a_isActive(c_neighborId(cellId, 1)))) {
      dvdx = (u_p[0][1] - u_m[0][1]) / (2.0 * dx);
      if constexpr(nDim == 3) dwdx = (u_p[0][2] - u_m[0][2]) / (2.0 * dx);
    }
    // Forwards difference in x-direction
    else if((!a_hasNeighbor(cellId, 0) || !a_isActive(c_neighborId(cellId, 0)))
            && (a_hasNeighbor(cellId, 1) && a_isActive(c_neighborId(cellId, 1)))) {
      dvdx = (u_p[0][1] - u[1]) / dx;
      if constexpr(nDim == 3) dwdx = (u_p[0][2] - u[2]) / dx;
    }
    // Backwards differences in x-direction
    else if((a_hasNeighbor(cellId, 0) && a_isActive(c_neighborId(cellId, 0)))
            && (!a_hasNeighbor(cellId, 1) || !a_isActive(c_neighborId(cellId, 1)))) {
      dvdx = (u[1] - u_m[0][1]) / dx;
      if constexpr(nDim == 3) dwdx = (u[2] - u_m[0][2]) / dx;
    }
    // Central differences in y-direction
    if((a_hasNeighbor(cellId, 2) && a_isActive(c_neighborId(cellId, 2)))
       && (a_hasNeighbor(cellId, 3) && a_isActive(c_neighborId(cellId, 3)))) {
      dudy = (u_p[1][0] - u_m[1][0]) / (2.0 * dx);
      if constexpr(nDim == 3) dwdy = (u_p[1][2] - u_m[1][2]) / (2.0 * dx);
    }
    // Forwards differences in y-direction
    else if((!a_hasNeighbor(cellId, 2) || !a_isActive(c_neighborId(cellId, 2)))
            && (a_hasNeighbor(cellId, 3) && a_isActive(c_neighborId(cellId, 3)))) {
      dudy = (u_p[1][0] - u[0]) / dx;
      if constexpr(nDim == 3) dwdy = (u_p[1][2] - u[2]) / dx;
    }
    // Backwards differences in y-direction
    else if((a_hasNeighbor(cellId, 2) && a_isActive(c_neighborId(cellId, 2)))
            && (!a_hasNeighbor(cellId, 3) || !a_isActive(c_neighborId(cellId, 3)))) {
      dudy = (u[0] - u_m[1][0]) / dx;
      if constexpr(nDim == 3) dwdy = (u[2] - u_m[1][2]) / dx;
    }
    if constexpr(nDim == 3) {
      MFloat dudz = 0.0;
      MFloat dvdz = 0.0;
      // Central differences in z-direction
      if((a_hasNeighbor(cellId, 4) && a_isActive(c_neighborId(cellId, 4)))
         && (a_hasNeighbor(cellId, 5) && a_isActive(c_neighborId(cellId, 5)))) {
        dudz = (u_p[2][0] - u_m[2][0]) / (2.0 * dx);
        dvdz = (u_p[2][1] - u_m[2][1]) / (2.0 * dx);
      }
      // Forwards difference in z-direction
      else if((!a_hasNeighbor(cellId, 4) || !a_isActive(c_neighborId(cellId, 4)))
              && (a_hasNeighbor(cellId, 5) && a_isActive(c_neighborId(cellId, 5)))) {
        dudz = (u_p[2][0] - u[0]) / dx;
        dvdz = (u_p[2][1] - u[1]) / dx;
      }
      // Backwards differences in z-direction
      else if((a_hasNeighbor(cellId, 4) && a_isActive(c_neighborId(cellId, 4)))
              && (!a_hasNeighbor(cellId, 5) || !a_isActive(c_neighborId(cellId, 5)))) {
        dudz = (u[0] - u_m[2][0]) / dx;
        dvdz = (u[1] - u_m[2][1]) / dx;
      }

      phi_v = sqrt(pow((dwdy - dvdz), 2) + pow((dudz - dwdx), 2) + pow((dvdx - dudy), 2)) * pow(dx, cellSizeWeight_v);
    }
    if constexpr(nDim == 2) phi_v = abs(dvdx - dudy) * pow(dx, cellSizeWeight_v);

    sensors[sensorOffset + sen][gridId] = phi_v;
    sensorCellFlag[gridId][sensorOffset + sen] = 1;
  }
}

/** \page sensorsLB
 *
 *  \section divergence Divergence
 *
 *  Property: <code>DIVERGENCE</code>
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbSolverDxQy<nDim, nDist, SysEqn>::sensorDivergence(std::vector<std::vector<MFloat>>& sensors,
                                                         std::vector<std::bitset<64>>& sensorCellFlag,
                                                         std::vector<MFloat>& sensorWeight, MInt sensorOffset,
                                                         MInt sen) {
  m_log << "   - Sensor preparation for the divergence sensor" << endl;
  cerr0 << "   - Sensor preparation for the divergence sensor" << endl;

  MFloat dx = 0.0;

  std::array<MFloat, nDim> u{};
  std::array<std::array<MFloat, nDim>, nDim> u_p{};
  std::array<std::array<MFloat, nDim>, nDim> u_m{};

  MFloat rho = 0.0;
  MFloat rho_p[nDim] = {0.0};
  MFloat rho_m[nDim] = {0.0};

  MFloat dudx = 0.0;
  MFloat dvdy = 0.0;

  MFloat phi_m = 0.0; // Sensor value mass conservation

  // Set sensor weights for refinement (K_f)
  sensorWeight[sensorOffset + sen] = this->m_sensorWeight[sen];

  // Set cell seize weihts
  // TODO labels:LB Make available by property?
  MFloat cellSizeWeight_m = 1.0;

  for(MInt cellId = 0; cellId < noInternalCells(); cellId++) {
    ASSERT(!a_isHalo(cellId), "");
    const MInt gridId = this->grid().tree().solver2grid(cellId);

    // Mesh adaptation is only performed on active cells
    // --> untag inactive cells
    if(!a_isActive(cellId)) {
      sensorCellFlag[gridId][sensorOffset + sen] = 0;
      continue;
    }

    // Get cell size for cell size weighting
    dx = grid().cellLengthAtLevel(c_level(cellId));

    // Compute macroscopic variables for sensor evaluation
    // All cells performed propagation step before. Thus, in-coming distributions
    // of all cells are in the same time step and comparable
    // --> Compute macroscopic variables from in-coming distributions.
    //
    // Depending on the neighbors of the cell, variables for
    // central, backwards or forwards differentials are computed
    calculateMacroscopicVariables(cellId, rho, u.data());
    if(a_hasNeighbor(cellId, 0) && a_isActive(c_neighborId(cellId, 0))) {
      MInt nghbrId = c_neighborId(cellId, 0);
      calculateMacroscopicVariables(nghbrId, rho_m[0], u_m[0].data());
    }
    if(a_hasNeighbor(cellId, 1) && a_isActive(c_neighborId(cellId, 1))) {
      MInt nghbrId = c_neighborId(cellId, 1);
      calculateMacroscopicVariables(nghbrId, rho_p[0], u_p[0].data());
    }
    if(a_hasNeighbor(cellId, 2) && a_isActive(c_neighborId(cellId, 2))) {
      MInt nghbrId = c_neighborId(cellId, 2);
      calculateMacroscopicVariables(nghbrId, rho_m[1], u_m[1].data());
    }
    if(a_hasNeighbor(cellId, 3) && a_isActive(c_neighborId(cellId, 3))) {
      MInt nghbrId = c_neighborId(cellId, 3);
      calculateMacroscopicVariables(nghbrId, rho_p[1], u_p[1].data());
    }
    if constexpr(nDim == 3) {
      if(a_hasNeighbor(cellId, 4) && a_isActive(c_neighborId(cellId, 4))) {
        MInt nghbrId = c_neighborId(cellId, 4);
        calculateMacroscopicVariables(nghbrId, rho_m[2], u_m[2].data());
      }
      if(a_hasNeighbor(cellId, 5) && a_isActive(c_neighborId(cellId, 5))) {
        MInt nghbrId = c_neighborId(cellId, 5);
        calculateMacroscopicVariables(nghbrId, rho_p[2], u_p[2].data());
      }
    }

    // Central difference in x-direction
    if((a_hasNeighbor(cellId, 0) && a_isActive(c_neighborId(cellId, 0)))
       && (a_hasNeighbor(cellId, 1) && a_isActive(c_neighborId(cellId, 1)))) {
      dudx = (rho_p[0] * u_p[0][0] - rho_m[0] * u_m[0][0]) / (2.0 * dx);
    }
    // Forwards differnce in x-direction
    else if((!a_hasNeighbor(cellId, 0) || !a_isActive(c_neighborId(cellId, 0)))
            && (a_hasNeighbor(cellId, 1) && a_isActive(c_neighborId(cellId, 1)))) {
      dudx = (rho_p[0] * u_p[0][0] - rho * u[0]) / dx;
    }
    // Backwards difference in x-direction
    else if((a_hasNeighbor(cellId, 0) && a_isActive(c_neighborId(cellId, 0)))
            && (!a_hasNeighbor(cellId, 1) || !a_isActive(c_neighborId(cellId, 1)))) {
      dudx = (rho * u[0] - rho_m[0] * u_m[0][0]) / dx;
    }
    // Central differences in y-direction
    if((a_hasNeighbor(cellId, 2) && a_isActive(c_neighborId(cellId, 2)))
       && (a_hasNeighbor(cellId, 3) && a_isActive(c_neighborId(cellId, 3)))) {
      dvdy = (rho_p[1] * u_p[1][1] - rho_m[1] * u_m[1][1]) / (2.0 * dx);
    }
    // Forwards differences in y-direction
    else if((!a_hasNeighbor(cellId, 2) || !a_isActive(c_neighborId(cellId, 2)))
            && (a_hasNeighbor(cellId, 3) && a_isActive(c_neighborId(cellId, 3)))) {
      dvdy = (rho_m[1] * u_p[1][1] - rho * u[1]) / dx;
    }
    // Backwars differences in y-direction
    else if((a_hasNeighbor(cellId, 2) && a_isActive(c_neighborId(cellId, 2)))
            && (!a_hasNeighbor(cellId, 3) || !a_isActive(c_neighborId(cellId, 3)))) {
      dvdy = (rho * u[1] - rho_m[1] * u_m[1][1]) / dx;
    }
    if constexpr(nDim == 3) {
      MFloat dwdz = 0.0;
      // Central differences in z-direction
      if((a_hasNeighbor(cellId, 4) && a_isActive(c_neighborId(cellId, 4)))
         && (a_hasNeighbor(cellId, 5) && a_isActive(c_neighborId(cellId, 5)))) {
        dwdz = (rho_p[2] * u_p[2][2] - rho_m[2] * u_m[2][2]) / (2.0 * dx);
      }
      // Forwards differences in z-direction
      else if((!a_hasNeighbor(cellId, 4) || !a_isActive(c_neighborId(cellId, 4)))
              && (a_hasNeighbor(cellId, 5) && a_isActive(c_neighborId(cellId, 5)))) {
        dwdz = (rho_p[2] * u_p[2][2] - rho * u[2]) / dx;
      }
      // Backwars differences in z-direction
      else if((a_hasNeighbor(cellId, 4) && a_isActive(c_neighborId(cellId, 4)))
              && (!a_hasNeighbor(cellId, 5) || !a_isActive(c_neighborId(cellId, 5)))) {
        dwdz = (rho * u[2] - rho_m[2] * u_m[2][2]) / dx;
      }

      phi_m = abs(dudx + dvdy + dwdz) * pow(dx, cellSizeWeight_m);
      // Test: Do not use divergence but use gradient of solution!
      // phi_m = (abs(dudx)+abs(dvdy)+abs(dwdz))*pow(dx,2.0);
    }
    if constexpr(nDim == 2) phi_m = abs(dudx + dvdy) * pow(dx, cellSizeWeight_m);

    sensors[sensorOffset + sen][gridId] = phi_m;
    sensorCellFlag[gridId][sensorOffset + sen] = 1;
  }
}

/** \page sensorsLB
 *
 *  \section totalPressure Total pressure
 *
 *  Property: <code>TOTALPRESSURE</code>
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbSolverDxQy<nDim, nDist, SysEqn>::sensorTotalPressure(std::vector<std::vector<MFloat>>& sensors,
                                                            std::vector<std::bitset<64>>& sensorCellFlag,
                                                            std::vector<MFloat>& sensorWeight, MInt sensorOffset,
                                                            MInt sen) {
  m_log << "   - Sensor preparation for the total pressure sensor" << endl;
  cerr0 << "   - Sensor preparation for the total pressure sensor" << endl;

  // Set sensor weights for refinement (K_f)
  sensorWeight[sensorOffset + sen] = this->m_sensorWeight[sen];

  // Set cell seize weihts
  // TODO labels:LB Make available by property?
  MFloat cellSizeWeight_p = 1.0;

  for(MInt cellId = 0; cellId < noInternalCells(); cellId++) {
    ASSERT(!a_isHalo(cellId), "");
    const MInt gridId = this->grid().tree().solver2grid(cellId);

    // Mesh adaptation is only performed on active cells
    // --> untag inactive cells
    if(!a_isActive(cellId)) {
      sensorCellFlag[gridId][sensorOffset + sen] = 0;
      continue;
    }

    // Get cell size for cell size weighting
    const MFloat dx = grid().cellLengthAtLevel(c_level(cellId));

    // Compute macroscopic variables for sensor evaluation
    // All cells performed propagation step before. Thus, in-coming distributions
    // of all cells are in the same time step and comparable
    // --> Compute macroscopic variables from in-coming distributions.
    //
    // Depending on the neighbors of the cell, variables for
    // central, backwards or forwards differentials are computed
    MFloat rho{};
    std::array<MFloat, nDim> u{};
    calculateMacroscopicVariables(cellId, rho, u.data());

    // (1) Difference in total pressure
    // Incident flow total pressure depends on bc!
    // Here so far only for bc10000 --> Assume rho=1 on boundary.
    const MFloat rho_0 = 1.0;
    const MFloat p_0 = rho_0 / 3.0;
    const MFloat p = rho / 3.0;
    const MFloat pTotal_0 = 0.5 * rho_0 * (m_Ma * LBCS) * (m_Ma * LBCS) + p_0;
    const MFloat squaredVelocity = std::inner_product(u.begin(), u.end(), u.begin(), 0.0);
    const MFloat pTotal = 0.5 * rho * squaredVelocity + p;
    const MFloat phi_p = abs(pTotal_0 - pTotal) * pow(dx, cellSizeWeight_p);
    sensors[sensorOffset + sen][gridId] = phi_p;
    sensorCellFlag[gridId][sensorOffset + sen] = 1;
  }
}

/** \page sensorsLB
 *
 *  \section meanStress Mean stress
 *
 *  @warning Experimental and only for 2D right now.
 *  Property: <code>MEANSTRESS</code>
 */
// TODO: EXPERIMENTAL, should be reworked!
template <MInt nDim, MInt nDist, class SysEqn>
void LbSolverDxQy<nDim, nDist, SysEqn>::sensorMeanStress(std::vector<std::vector<MFloat>>& sensors,
                                                         std::vector<std::bitset<64>>& sensorCellFlag,
                                                         std::vector<MFloat>& sensorWeight, MInt sensorOffset,
                                                         MInt sen) {
  if(nDim == 3) {
    TERMM(1, "Not implemented for 3D");
  }


  m_log << "   - Sensor preparation for the mean stress sensor" << endl;
  cerr0 << "   - Sensor preparation for the mean stress sensor" << endl;

  // Set sensor weights for refinement (K_f)
  sensorWeight[sensorOffset + sen] = this->m_sensorWeight[sen];

  // Set cell seize weihts
  // TODO: Make available by property?
  const MFloat cellSizeWeight = 1.0;

  for(MInt cellId = 0; cellId < noInternalCells(); cellId++) {
    ASSERT(!a_isHalo(cellId), "");
    const MInt gridId = this->grid().tree().solver2grid(cellId);

    // Mesh adaptation is only performed on active cells
    // --> untag inactive cells
    if(!a_isActive(cellId)) {
      sensorCellFlag[gridId][sensorOffset + sen] = 0;
      continue;
    }

    // Calculate mean stress
    const MFloat rho = a_variable(cellId, PV->RHO);
    std::array<MFloat, 2> u = {a_variable(cellId, PV->U), a_variable(cellId, PV->V)};

    std::array<MFloat, nDist> eqDist{};
    eqDist = getEqDists(rho, u.data());

    // Calculation of non-eq-parts
    std::array<MFloat, 9> d{};
    for(MInt j = 0; j < 9; j++) {
      d[j] = a_oldDistribution(cellId, j) - eqDist[j];
    }

    // Calculation of momentum flux tensor for non-equilibrium parts
    std::array<MFloat, 4> Q{};
    Q[0] = (d[0] + d[1] + d[4] + d[5] + d[6] + d[7]);
    Q[1] = (d[4] - d[5] + d[6] - d[7]);

    Q[2] = Q[1];
    Q[3] = (d[2] + d[3] + d[4] + d[5] + d[6] + d[7]);

    MFloat Qmean = 0.0;
    for(MInt k = 0; k < 4; k++) {
      Qmean += Q[k] * Q[k];
    }

    Qmean = sqrt(2 * Qmean);

    // Get cell size for cell size weighting
    const MFloat dx = grid().cellLengthAtLevel(c_level(cellId));

    const MFloat phi = Qmean * pow(dx, cellSizeWeight);

    sensors[sensorOffset + sen][gridId] = phi;
    sensorCellFlag[gridId][sensorOffset + sen] = 1;
  }
}

/**
 * \brief Restart bndCnd object
 * \author: Philipp Brokof, Moritz Waldmann
 **/
template <MInt nDim, MInt nDist, class SysEqn>
void LbSolverDxQy<nDim, nDist, SysEqn>::restartBndCnd() {
  if(m_isRefined || this->m_adaptation) {
    delete m_interface;
  }

  if(m_isRefined || this->m_adaptation) m_interface = new LbInterfaceDxQy<nDim, nDist, SysEqn>(this);

  delete m_bndCnd;
  m_bndCnd = new LbBndCnd(this);
  LbSolver<nDim>::m_bndCnd = m_bndCnd;
  m_bndCnd->initializeBcData();
}

template <MInt nDim, MInt nDist, class SysEqn>
void LbSolverDxQy<nDim, nDist, SysEqn>::calcNodalLsValues() {
  TRACE();

  // TODO labels:LB Unify for 2D/3D (as far as possible).

  if constexpr(nDim == 3) {
    MInt startSet = this->m_levelSetId;
    static constexpr MInt maxNoDistributions = IPOW3[nDim];
    static constexpr MInt firstCornerDist = 18;
    static constexpr MInt firstEdgeDist = 6;
    static constexpr MInt noEdges = 12;

    // TODO labels:LB Use lattice descriptor.

    static constexpr MInt lb_CornerMappingForSurfaces[6][4] = {{18, 19, 20, 21}, {22, 23, 24, 25}, {18, 19, 22, 23},
                                                               {20, 21, 24, 25}, {18, 20, 22, 24}, {19, 21, 23, 25}};

    static constexpr MInt lb_CornerMappingForEdges[12][2] = {{18, 19}, {20, 21}, {22, 23}, {24, 25},
                                                             {18, 20}, {19, 21}, {22, 24}, {23, 25},
                                                             {18, 22}, {19, 23}, {20, 24}, {21, 25}};

    static constexpr MInt lb_geomIntersToLb[8] = {0, 4, 2, 6, 1, 5, 3, 7};

    static constexpr MFloat factor[26] = {F1,    F1,    F1,    F1,    F1,    F1,    SQRT2, SQRT2, SQRT2,
                                          SQRT2, SQRT2, SQRT2, SQRT2, SQRT2, SQRT2, SQRT2, SQRT2, SQRT2,
                                          SQRT3, SQRT3, SQRT3, SQRT3, SQRT3, SQRT3, SQRT3, SQRT3};

    // 1. Get the nodal Ls values form the CutCandidates Class
    // 2. Find the Cut Cells and calculate the missing Ls values
    // at the cell surfaces and cell edges for those
    // 3. Calculate the Wall distances
    this->m_currentNoG0Cells = 0;
    m_bndCnd->m_boundaryCellsMb.noDistances(nDist - 1);
    m_bndCnd->m_boundaryCellsMb.reset(this->m_noG0CandidatesTotal);
    m_bndCnd->m_boundaryCellMappingMb.clear();

    for(MInt candidate = 0; candidate < this->m_noG0CandidatesTotal; candidate++) {
      for(MInt set = startSet; set < this->m_maxNoSets; set++) {
        const MInt arrayPosNGV = set * (maxNoDistributions - 1);
        for(MInt i = 0; i < IPOW2(nDim); i++) {
          this->m_nodalGValues[candidate][arrayPosNGV + firstCornerDist + lb_geomIntersToLb[i]] =
              this->m_G0Candidates[candidate].nodalValues[set][i];
        }

        // Step a: Find Edges with a cut, i.e. the nodes of the edge have opposite signs in m_nodalGValues
        MInt edgeCounter = 0;
        for(MInt edge = 0; edge < noEdges; edge++) {
          MFloat edgeGValue = this->m_nodalGValues[candidate][arrayPosNGV + Ld::nodalConnectivityVector(edge, 0)]
                              * this->m_nodalGValues[candidate][arrayPosNGV + Ld::nodalConnectivityVector(edge, 1)];

          if(edgeGValue < 0) edgeCounter++;
        }

        // Step b: If there are less than two cuts for a G0Candidate it is not a G0 Cell. Skip that cell!!
        const MInt pCellId = this->m_G0Candidates[candidate].cellId;
        const MInt actualSet = set - startSet;
        if(edgeCounter < 2) {
          this->a_isG0CandidateOfSet(pCellId, actualSet) = false;
          continue;
        } else {
          this->a_isG0CandidateOfSet(pCellId, actualSet) = true;
          this->m_G0CellMapping[pCellId] = this->m_currentNoG0Cells;
          // Use new mapping in bndCnd
          m_bndCnd->m_boundaryCellMappingMb[pCellId] = this->m_currentNoG0Cells;
        }

        // Add relevant information to G0 boundary cell collector
        // Cell is a G0 boundary cell
        const MInt boundaryCellId = m_bndCnd->m_boundaryCellsMb.size();
        m_bndCnd->m_boundaryCellsMb.append();
        m_bndCnd->m_boundaryCellsMb.cellId(boundaryCellId) = pCellId;

        // Step c: Calculate the missing Ls values at the cell surfaces and cell edges
        this->m_G0CellList[this->m_currentNoG0Cells] = pCellId;
        for(MInt dist = 0; dist < firstEdgeDist; dist++) {
          for(MInt i = 0; i < 4; i++) {
            this->m_nodalGValues[candidate][arrayPosNGV + dist] +=
                this->m_nodalGValues[candidate][arrayPosNGV + lb_CornerMappingForSurfaces[dist][i]];
          }
          this->m_nodalGValues[candidate][arrayPosNGV + dist] /= F4;
        }
        for(MInt dist = firstEdgeDist; dist < firstCornerDist; dist++) {
          for(MInt i = 0; i < 2; i++) {
            this->m_nodalGValues[candidate][arrayPosNGV + dist] +=
                this->m_nodalGValues[candidate][arrayPosNGV + lb_CornerMappingForEdges[dist - firstEdgeDist][i]];
          }
          this->m_nodalGValues[candidate][arrayPosNGV + dist] /= F2;
        }


        MInt arrayPosWallDists = (set - startSet) * (nDist - 1);
        // Step d: Check which distributions are cutted and calculate the Wall distance for those
        for(MInt dir = 0; dir < nDist - 1; dir++) {
          if((this->a_levelSetFunctionMB(pCellId, set) * this->m_nodalGValues[candidate][arrayPosNGV + dir]) > 0) {
            // TODO labels:LB Remove offset and use set as third index ...
            m_bndCnd->m_boundaryCellsMb.distance(boundaryCellId, arrayPosWallDists + dir) = F2;
          } else {
            MFloat levelSetRatio =
                abs(this->m_nodalGValues[candidate][arrayPosNGV + dir] / this->a_levelSetFunctionMB(pCellId, set));
            MFloat cellLength = this->c_cellLengthAtLevel(this->a_level(pCellId)) * F1B2 * factor[dir];
            m_bndCnd->m_boundaryCellsMb.distance(boundaryCellId, arrayPosWallDists + dir) =
                cellLength / (levelSetRatio + F1);
          }

          // Determine surface center
          // TODO labels:LB Fix very bad estimate by boundary cell center using approx normal (minWallDist)
          for(MInt n = 0; n < nDim; n++) {
            m_bndCnd->m_boundaryCellsMb.surfaceCenter(boundaryCellId, dir, n) =
                a_coordinate(pCellId, n)
                + m_bndCnd->m_boundaryCellsMb.distance(boundaryCellId, dir) * Ld::ppdfDir(dir, n) * 1 / factor[dir];
          }
        }
        for(MInt n = 0; n < nDim; n++) {
          m_bndCnd->m_boundaryCellsMb.cellCenter(boundaryCellId, n) = a_coordinate(pCellId, n);
        }

        this->m_currentNoG0Cells++;
      }
    }
  } else {
    // D2Q9 Stencil for Cut Cells including the edges:
    //
    //                     14   15
    //                   7 -- 3 -- 4
    //                   |\   |   /|
    //                 9 | \  |  / | 11
    //                   |  \ | /  |
    //                   0 -- * -- 1
    //                   |  / | \  |
    //                 8 | /  |  \ | 10
    //                   |/   |   \|
    //                   6 -- 2 -- 5
    //                     12   13
    //
    // 1. Get the nodal Ls values form the CutCandidates Class
    //   m_nodelGValues holds the LS values for each candidate in the order of the distributions for each set

    // 2. Find the Cut Cells by multipling the nodal LS values of the 4 corners
    //   If the two LS values have opposite signs the edge has a cut if the product is exaktly 0 the cell has to
    //   be considered, since it is still a boundary cell.
    //   Furthermore, the cutPoints are calculated for all edges with a cut. The eges are seperated in two part
    //   which is needed for the calculation of the cut cell volume fraction

    // 3. Calculate the missing Ls values at the cell surfaces and cell edges for the cut cells
    //   by calculating the mean values of the surrounding corners.

    // 4. Calculate the Wall distances using the formular x2 = cellLength / (1 + LS2/LS1)

    const MInt startSet = this->m_levelSetId;
    constexpr MInt maxNoDistributions = 9;
    constexpr MInt firstCornerDist = 4;
    constexpr MInt noEdges = 4;
    // 1. Get the nodal Ls values form the CutCandidates Class
    // 2. Find the Cut Cells and calculate the missing Ls values
    // at the cell surfaces and cell edges for those
    // 3. Calculate the Wall distances
    m_bndCnd->m_boundaryCellsMb.noDistances(nDist - 1);
    m_bndCnd->m_boundaryCellsMb.reset(this->m_noG0CandidatesTotal);
    constexpr MInt lb_CornerMappingForEdges[4][2] = {{6, 7}, {4, 5}, {5, 6}, {7, 4}};
    constexpr MFloat factor[8] = {F1, F1, F1, F1, SQRT2, SQRT2, SQRT2, SQRT2};
    // 1. Get the nodal Ls values form the CutCandidates Class
    // 2. Find the Cut Cells and calculate the missing Ls values
    // at the cell surfaces and cell edges for those
    // 3. Calculate the Wall distances

    // Reset data structure
    this->m_currentNoG0Cells = 0;
    m_bndCnd->m_boundaryCellsMb.noDistances(nDist - 1);
    m_bndCnd->m_boundaryCellsMb.reset(this->m_noG0CandidatesTotal);
    m_bndCnd->m_boundaryCellMappingMb.clear();

    for(MInt candidate = 0; candidate < this->m_noG0CandidatesTotal; candidate++) {
      for(MInt set = startSet; set < this->m_maxNoSets; set++) {
        MInt arrayPosNGV = set * (maxNoDistributions - 1);
        this->m_nodalGValues[candidate][arrayPosNGV + firstCornerDist + 2] =
            this->m_G0Candidates[candidate].nodalValues[set][0];
        this->m_nodalGValues[candidate][arrayPosNGV + firstCornerDist + 1] =
            this->m_G0Candidates[candidate].nodalValues[set][1];
        this->m_nodalGValues[candidate][arrayPosNGV + firstCornerDist + 3] =
            this->m_G0Candidates[candidate].nodalValues[set][2];
        this->m_nodalGValues[candidate][arrayPosNGV + firstCornerDist + 0] =
            this->m_G0Candidates[candidate].nodalValues[set][3];

        // Step a: Find Edges with a cut,
        // i.e. the nodes of the edge have opposite signs in m_nodalGValues
        MInt edgeCounter = 0;
        for(MInt edge = 0; edge < noEdges; edge++) {
          MFloat edgeGValue = this->m_nodalGValues[candidate][arrayPosNGV + Ld::nodalConnectivityVector(edge, 0)]
                              * this->m_nodalGValues[candidate][arrayPosNGV + Ld::nodalConnectivityVector(edge, 1)];
          if(edgeGValue <= 0) edgeCounter++;
        }

        // Step b: If there are less than two cuts for a G0Candidate it is not a G0 Cell.
        // Skip that cell!!
        MInt pCellId = this->m_G0Candidates[candidate].cellId;
        MInt actualSet = set - startSet;
        if(edgeCounter < 2) {
          this->a_isG0CandidateOfSet(pCellId, actualSet) = false;
          continue;
        } else {
          this->a_isG0CandidateOfSet(pCellId, actualSet) = true;
          this->m_G0CellMapping[pCellId] = this->m_currentNoG0Cells;
          // Use new mapping in bndCnd
          m_bndCnd->m_boundaryCellMappingMb[pCellId] = this->m_currentNoG0Cells;
        }

        // Add relevant information to G0 boundary cell collector
        // Cell is a G0 boundary cell
        const MInt boundaryCellId = m_bndCnd->m_boundaryCellsMb.size();
        m_bndCnd->m_boundaryCellsMb.append();
        m_bndCnd->m_boundaryCellsMb.cellId(boundaryCellId) = pCellId;

        // Step c: Calculate the missing Ls values at the cell edges
        this->m_G0CellList[this->m_currentNoG0Cells] = pCellId;
        for(MInt dist = 0; dist < firstCornerDist; dist++) {
          for(MInt i = 0; i < 2; i++) {
            this->m_nodalGValues[candidate][arrayPosNGV + dist] +=
                this->m_nodalGValues[candidate][arrayPosNGV + lb_CornerMappingForEdges[dist][i]];
          }
          this->m_nodalGValues[candidate][arrayPosNGV + dist] /= F2;
        }

        MInt arrayPosWallDists = (set - startSet) * (nDist - 1);
        // Step d: Check which distributions are cutted and calculate the Wall distance for those
        for(MInt dir = 0; dir < nDist - 1; dir++) {
          if((this->a_levelSetFunctionMB(pCellId, set) * this->m_nodalGValues[candidate][arrayPosNGV + dir]) > 0) {
            // TODO labels:LB Remove offset and use set as third index ...
            m_bndCnd->m_boundaryCellsMb.distance(boundaryCellId, arrayPosWallDists + dir) = F2;
          } else {
            const MFloat levelSetRatio =
                abs(this->m_nodalGValues[candidate][arrayPosNGV + dir] / this->a_levelSetFunctionMB(pCellId, set));
            const MFloat cellLength = this->c_cellLengthAtLevel(this->a_level(pCellId)) * F1B2 * factor[dir];
            // TODO labels:LB Remove offset and use set as third index ...
            m_bndCnd->m_boundaryCellsMb.distance(boundaryCellId, arrayPosWallDists + dir) =
                cellLength / (levelSetRatio + F1);
          }

          // Determine surface center
          // TODO labels:LB Fix very bad estimate by boundary cell center using approx normal (minWallDist)
          // OR: Store cross-over point of every cut-direction to make the
          // angular momentum exchange even better!
          for(MInt n = 0; n < nDim; n++) {
            m_bndCnd->m_boundaryCellsMb.surfaceCenter(boundaryCellId, dir, n) =
                a_coordinate(pCellId, n)
                + m_bndCnd->m_boundaryCellsMb.distance(boundaryCellId, n) * Ld::ppdfDir(dir, n) * 1 / factor[dir];
          }
        }
        for(MInt n = 0; n < nDim; n++) {
          m_bndCnd->m_boundaryCellsMb.cellCenter(boundaryCellId, n) = a_coordinate(pCellId, n);
        }

        this->m_currentNoG0Cells++;
      }
    }
  }
}

template <MInt nDim, MInt nDist, class SysEqn>
MBool LbSolverDxQy<nDim, nDist, SysEqn>::maxResidual() {
  TRACE();

  if(globalTimeStep % m_residualInterval == 0) {
    calculateResidual();
  }

  return true;
}

/** brief activate all cells, but the halo cells
 *
 *  The pressure force is calculated on the lowest level of Refinement.
 *  One half of the forcing term is added after the collision, and the other half is
 *  added after streaming. This ensures correct forcing at half-way bounce-back boundaries.
 *
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbSolverDxQy<nDim, nDist, SysEqn>::initPressureForce() {
  TRACE();

  /*********************************************/
  /* 2D LAMINAR CHANNEL INITIALIZATION          */
  /* m_referenceLength must be channel half width !!! */
  /* Ma refers to velocity on centerline !!!   */
  /*********************************************/
  if(string2enum(m_initMethod) == LB_LAMINAR_CHANNEL_INIT) {
    // Re must be defined for channel half height

    m_nu = m_Ma * LBCS / m_Re * m_referenceLength;

    m_densityGradient = m_Ma * LBCS * 9.0 * m_nu / (m_referenceLength * m_referenceLength);

    constexpr MInt dir = 0;
    for(MInt mi = 0; mi < Ld::dxQyFld(); mi++) {
      m_Fext[Ld::nFld(dir, mi)] = Ld::tp(Ld::distType(mi)) * -1.0 * m_densityGradient;
      m_Fext[Ld::pFld(dir, mi)] = Ld::tp(Ld::distType(mi)) * 1.0 * m_densityGradient;
    }
  }

  /*********************************************/
  /* TURBULENT CHANNEL INITIALIZATION          */
  /* m_referenceLength must be channel half width !!! */
  /* Ma refers to velocity on centerline !!!   */
  /*********************************************/
  if(string2enum(m_initMethod) == LB_TURBULENT_CHANNEL_INIT) {
    // For now testing only the D3Q19 algorithm
    MFloat uTau;

    // Re must be defined for channel half width
    m_nu = m_Ma * LBCS / m_Re * m_referenceLength;
    uTau = (MFloat)m_ReTau * m_nu / m_referenceLength;

    // The ratio between the turb. channel length and channel half width : ratio = L/H
    // MFloat ratio = m_domainLength/m_referenceLength;

    m_densityGradient = 3.0 * uTau * uTau / m_referenceLength;

    constexpr MInt dir = 0;
    for(MInt mi = 0; mi < Ld::dxQyFld(); mi++) {
      m_Fext[Ld::nFld(dir, mi)] = Ld::tp(Ld::distType(mi)) * -1.0 * m_densityGradient;
      m_Fext[Ld::pFld(dir, mi)] = Ld::tp(Ld::distType(mi)) * 1.0 * m_densityGradient;
    }
  }

  /**************************************/
  /* TURBULENT DUCT INITIALIZATION      */
  /**************************************/
  if(string2enum(m_initMethod) == LB_TURBULENT_DUCT_INIT) {
    // For now testing only the D3Q19 algorithm
    MFloat uTau;

    uTau = (MFloat)m_ReTau / m_Re * m_Ma * LBCS;
    MFloat ratio = m_domainLength / (m_referenceLength * F1B2);

    m_densityGradient = 3.0 * uTau * uTau * ratio / m_domainLength;

    constexpr MInt dir = 0;
    for(MInt mi = 0; mi < Ld::dxQyFld(); mi++) {
      m_Fext[Ld::nFld(dir, mi)] = Ld::tp(Ld::distType(mi)) * -1.0 * m_densityGradient;
      m_Fext[Ld::pFld(dir, mi)] = Ld::tp(Ld::distType(mi)) * 1.0 * m_densityGradient;
    }
  }

  /*************************************************/
  /* TURBULENT PIPE INITIALIZATION                 */
  /* m_referenceLength must be pipe diameter !!!          */
  /* Ma refers to the maximum streamwise velocity !!! */
  /* u_mean = 0.816 * u_max                        */
  /*************************************************/
  if(string2enum(m_initMethod) == LB_TURBULENT_PIPE_INIT) {
    // For now testing only the D3Q19 algorithm
    MFloat uTau;

    // uTau = (MFloat)m_ReTau/m_Re * m_Ma * LBCS ;
    m_nu = m_Ma * LBCS / m_Re * m_referenceLength;
    uTau = (MFloat)m_ReTau * m_nu / m_referenceLength;

    /*
    MFloat ratio = m_domainLength/m_referenceLength;
    // pressure drop over entire length according to Blasius law (IN THIS CASE Ma DEFINES THE MEAN VELOCITY!)
    m_gradient = 0.3164 * pow(m_Re,-0.25) * ratio * (m_Ma*LBCS) * (m_Ma*LBCS) / 2;//rho=1
    */

    // density drop per cell
    m_densityGradient = 3.0 * 2.0 * uTau * uTau / (F1B2 * m_referenceLength); // Pope, p.293

    if constexpr(nDim != 3) {
      std::stringstream ss;
      ss << "Init method " << m_initMethod << " is only available for 3D." << std::endl;
      TERMM(1, ss.str());
    }
    constexpr MInt dir = 2;
    for(MInt mi = 0; mi < Ld::dxQyFld(); mi++) {
      m_Fext[Ld::nFld(dir, mi)] = Ld::tp(Ld::distType(mi)) * -1.0 * m_densityGradient;
      m_Fext[Ld::pFld(dir, mi)] = Ld::tp(Ld::distType(mi)) * 1.0 * m_densityGradient;
    }
  }

  /*******************************/
  /* LAMINAR PIPE INITIALIZATION */
  /*******************************/
  if(string2enum(m_initMethod) == LB_LAMINAR_PIPE_INIT) {
    m_nu = m_Ma * LBCS / m_Re * m_referenceLength;

    // law of Hagen Poiseuille
    // The pressure gradient is one third of the density gradient !
    //    m_gradient = 3*F1B2*(m_Ma*LBCS)*m_domainLength/m_referenceLength/m_referenceLength*64.0*m_nu;
    m_densityGradient =
        3.0 * 8.0 * (m_Ma * LBCS) * m_nu / ((F1B2 * m_referenceLength) * (F1B2 * m_referenceLength)); // rho=1

    if constexpr(nDim != 3) {
      std::stringstream ss;
      ss << "Init method " << m_initMethod << " is only available for 3D." << std::endl;
      TERMM(1, ss.str());
    }
    constexpr MInt dir = 2;
    for(MInt mi = 0; mi < Ld::dxQyFld(); mi++) {
      m_Fext[Ld::nFld(dir, mi)] = Ld::tp(Ld::distType(mi)) * -1.0 * m_densityGradient;
      m_Fext[Ld::pFld(dir, mi)] = Ld::tp(Ld::distType(mi)) * 1.0 * m_densityGradient;
    }
  }
}

template <MInt nDim, MInt nDist, class SysEqn>
void LbSolverDxQy<nDim, nDist, SysEqn>::initVolumeForces() {
  static MBool firstCall = true;
  static MFloat origFext[nDist];
  if(firstCall) {
    // save the incoming m_Fext to recall them at re-initializations triggered by velocity control
    firstCall = false;
    for(MInt i = 0; i < nDist; i++) {
      origFext[i] = m_Fext[i];
    }
  }
  for(MInt dir = 0; dir < nDim; dir++) {
    for(MInt mi = 0; mi < Ld::dxQyFld(); mi++) {
      m_Fext[Ld::nFld(dir, mi)] =
          origFext[Ld::nFld(dir, mi)] + Ld::tp(Ld::distType(Ld::nFld(dir, mi))) * -1.0 * m_volumeAccel[dir] * 3.0;
      m_Fext[Ld::pFld(dir, mi)] =
          origFext[Ld::pFld(dir, mi)] + Ld::tp(Ld::distType(Ld::pFld(dir, mi))) * 1.0 * m_volumeAccel[dir] * 3.0;
    }
  }

  if(m_isEELiquid) {
    for(MInt dir = 0; dir < nDim; dir++) {
      for(MInt mi = 0; mi < Ld::dxQyFld(); mi++) {
        m_EELiquid.Fg[Ld::nFld(dir, mi)] =
            Ld::tp(Ld::distType(Ld::nFld(dir, mi))) * -1.0 * m_EELiquid.gravityAccelM[dir] * 3.0;
        m_EELiquid.Fg[Ld::pFld(dir, mi)] =
            Ld::tp(Ld::distType(Ld::pFld(dir, mi))) * 1.0 * m_EELiquid.gravityAccelM[dir] * 3.0;
      }
    }
  }
}

/** \brief  Iterative initialize routine to obtained a valid density and non-eq field
 *  \author Miro Gondrum
 *  \date   19.01.2022
 *
 * For the init method following Mei et al. the initial velocity field is kept
 * constant and only density/pressure and the non-equilibrium field is updated.
 * This function does the correction independent of chosen methods. This is done
 * by Resetting
 * 1) .. macroscopic state from (rho, u) to (rho, u0)
 * 2) .. mesoscopic state from f=f_neq + f_eq(rho,u) to f=f_neq + f_eq(rho,u0)
 * Thus, LBM solves only a Poisson equation for the pressure.
 * ref.: Mei et al. 2006: https://doi.org/10.1016/j.compfluid.2005.08.008
 */
template <MInt nDim, MInt nDist, class SysEqn>
template <MBool compressible>
void LbSolverDxQy<nDim, nDist, SysEqn>::initRunCorrection_() {
  TRACE();
  // Required for ported functions, since GPUs cannot access the global variable globalTimeStep
  const MInt gTS = globalTimeStep;

  maia::parallelFor<true>(0, m_currentMaxNoCells, [=](MInt i) {
    const MInt pCellId = m_activeCellList[i];
    const MInt lvlDiff = maxLevel() - c_level(pCellId);
    if((gTS) % IPOW2(lvlDiff) != 0) return;
      // In the following step only rho is updated while u is kept
      // constant with its initial state u0.
      // 1) Calculate rho, rho*u and the associated f_eq(rho, u)
#ifdef WAR_NVHPC_PSTL
    MFloat u[nDim] = {F0};
    for(MInt d = 0; d < nDim; d++) {
      u[d] = a_variable(pCellId, d);
    }
    calculateMacroscopicVariables<compressible>(pCellId, a_variable(pCellId, PV->RHO), u);
    for(MInt d = 0; d < nDim; d++) {
      a_variable(pCellId, d) = u[d];
    }
#else
    calculateMacroscopicVariables<compressible>(pCellId, a_variable(pCellId, PV->RHO), &a_variable(pCellId, PV->U));
#endif
    MFloat eqDist[nDist], trgEqDist[nDist];
#ifdef WAR_NVHPC_PSTL
    for(MInt d = 0; d < nDim; d++) {
      u[d] = a_variable(pCellId, d);
    }
    sysEqn().calcEqDists(a_variable(pCellId, PV->RHO), u, &eqDist[0], m_mFld1.data(), m_mFld2.data(), m_tp.data(),
                         m_distFld.data());
#else
    sysEqn().calcEqDists(a_variable(pCellId, PV->RHO), &a_variable(pCellId, PV->U), &eqDist[0]);
#endif
    // 2) Set state to (rho, u0) and calculate associated f_eq(rho, u0)
    for(MInt d = 0; d < nDim; d++) {
      a_variable(pCellId, PV->U + d) = a_oldVariable(pCellId, PV->U + d);
    }
#ifdef WAR_NVHPC_PSTL
    for(MInt d = 0; d < nDim; d++) {
      u[d] = a_variable(pCellId, d);
    }
    sysEqn().calcEqDists(a_variable(pCellId, PV->RHO), u, &trgEqDist[0], m_mFld1.data(), m_mFld2.data(), m_tp.data(),
                         m_distFld.data());
#else
    sysEqn().calcEqDists(a_variable(pCellId, PV->RHO), &a_variable(pCellId, PV->U), &trgEqDist[0]);
#endif
    // 3) Reset f such that: f_new = f_eq(rho, u0) + f_neq(rho, u)
    for(MInt j = 0; j < nDist; j++) {
      a_oldDistribution(pCellId, j) = a_oldDistribution(pCellId, j) + trgEqDist[j] - eqDist[j];
    }
  });
}

template <MInt nDim, MInt nDist, class SysEqn>
void LbSolverDxQy<nDim, nDist, SysEqn>::initRunCorrection() {
  if(this->isCompressible()) {
    initRunCorrection_<true>();
  } else {
    initRunCorrection_<false>();
  }
}

/** \fn LbSolverDxQy::initLatticeBgkFftChannel()
 * \author Andreas Lintermann
 * \date 24.01.2011
 * \propVal{initMethod,LB_TURBULENT_CHANNEL_INIT}
 *
 * Initializes standard Lattice BGK Turbulent channel with FFT initialization.
 * Set property "initMethod" to "LB_TURBULENT_CHANNEL_INIT" and the property "FFTInit" to "1" to use this function.
 *
 * \LBIC{LbSolverDxQy::initLatticeBgkFftChannel(), initLatticeBgkFftChannel, FftChannelInit}
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbSolverDxQy<nDim, nDist, SysEqn>::initLatticeBgkFftChannel() {
  TRACE();

  MInt xPos, yPos, zPos;
  MInt lx = 0, ly = 0, lz = 0;

  // fft-domain dimensions
  // this holds the size of the domain in number of cells on lowest level
  lx = m_arraySize[0] / FPOW2(maxLevel() - minLevel());
  ly = m_arraySize[1] / FPOW2(maxLevel() - minLevel());
  if(nDim == 3) lz = m_arraySize[2] / FPOW2(maxLevel() - minLevel());

  fftw_complex *uPhysField, *vPhysField, *wPhysField;

  // field of velocities from positve frequencies
  uPhysField = (fftw_complex*)fftw_malloc(lx * ly * lz * sizeof(fftw_complex));
  vPhysField = (fftw_complex*)fftw_malloc(lx * ly * lz * sizeof(fftw_complex));
  wPhysField = (fftw_complex*)fftw_malloc(lx * ly * lz * sizeof(fftw_complex));

  // in multisolver mode the velocity field is created in solver 0 and broadcasted to all others
  ScratchSpace<MFloat> sendRecvBufferU(lx * ly * lz, AT_, "sendRecvBufferU");
  ScratchSpace<MFloat> sendRecvBufferV(lx * ly * lz, AT_, "sendRecvBufferV");
  ScratchSpace<MFloat> sendRecvBufferW(lx * ly * lz, AT_, "sendRecvBufferW");

  if(domainId() == 0) {
    // create velocity field
    maia::math::initFft(uPhysField, vPhysField, wPhysField, lx, ly, lz, m_noPeakModes, m_Ma);

    // copy values to buffer
    for(MInt p = 0; p < lx; p++) {
      for(MInt q = 0; q < ly; q++) {
        for(MInt r = 0; r < lz; r++) {
          sendRecvBufferU[r + lz * (q + ly * p)] = uPhysField[r + lz * (q + ly * p)][0];
          sendRecvBufferV[r + lz * (q + ly * p)] = vPhysField[r + lz * (q + ly * p)][0];
          sendRecvBufferW[r + lz * (q + ly * p)] = wPhysField[r + lz * (q + ly * p)][0];
        }
      }
    }
  }

  MPI_Bcast(&sendRecvBufferU[0], lx * ly * lz, MPI_DOUBLE, 0, mpiComm(), AT_, "sendRecvBufferU[0]");
  MPI_Bcast(&sendRecvBufferV[0], lx * ly * lz, MPI_DOUBLE, 0, mpiComm(), AT_, "sendRecvBufferV[0]");
  MPI_Bcast(&sendRecvBufferW[0], lx * ly * lz, MPI_DOUBLE, 0, mpiComm(), AT_, "sendRecvBufferW[0]");

  if(domainId() != 0) {
    // copy values from buffer
    for(MInt p = 0; p < lx; p++) {
      for(MInt q = 0; q < ly; q++) {
        for(MInt r = 0; r < lz; r++) {
          uPhysField[r + lz * (q + ly * p)][0] = sendRecvBufferU[r + lz * (q + ly * p)];
          vPhysField[r + lz * (q + ly * p)][0] = sendRecvBufferV[r + lz * (q + ly * p)];
          wPhysField[r + lz * (q + ly * p)][0] = sendRecvBufferW[r + lz * (q + ly * p)];
        }
      }
    }
  }

  MFloat actualCellLength;
  MFloat rho, u[3] = {0.0, 0.0, 0.0}, uTau, yPlus;

  // all cells have the same density
  if(m_densityFluctuations)
    rho = 0.0;
  else
    rho = 1.0;

  uTau = (MFloat)m_ReTau * m_nu / m_referenceLength;

  // tmpwidth equals diameter(in cells) * cellLength on maxLevel
  // diameter is half the channel height!
  MFloat tmpWidth = m_referenceLength * m_smallestCellLength;
  // tmpWidth = 0.5 * fabs(bBox[4]-bBox[1]);

  cerr << " --- prescribing mean velocity profile for turbulent channel ---" << endl;
  for(MInt i = 0; i < a_noCells(); i++) {
    initNu(i, m_nu);

    yPlus = m_ReTau * (1.0 - fabs(a_coordinate(i, 1)) / tmpWidth);

    // debug:
    if(yPlus < 0 || yPlus > m_ReTau) {
      stringstream errorMessage;
      errorMessage << "incorrect value for yPlus: " << yPlus << "   exiting...";
      TERMM(1, errorMessage.str());
    }

    // set the mean velocity profile
    if(yPlus <= 5.0) {
      u[0] = uTau * yPlus;
    } else {
      if(yPlus <= 30 && yPlus > 5.0) {
        u[0] = uTau * (C1 * log(yPlus) + C2);
      } else {
        if(yPlus > 30) {
          u[0] = uTau * (C3 * log(yPlus) + C4);
        }
      }
    }

    u[1] = 0.0;
    u[2] = 0.0;

    actualCellLength = this->grid().cellLengthAtCell(i);

    xPos = floor(
        (F1B2 * lx + (a_coordinate(i, 0) - F1B2 * (actualCellLength)) / (this->c_cellLengthAtLevel(minLevel()))) + 0.1);
    yPos = floor(
        (F1B2 * ly + (a_coordinate(i, 1) - F1B2 * (actualCellLength)) / (this->c_cellLengthAtLevel(minLevel()))) + 0.1);
    zPos = floor(
        (F1B2 * lz + (a_coordinate(i, 2) - F1B2 * (actualCellLength)) / (this->c_cellLengthAtLevel(minLevel()))) + 0.1);

    if(xPos > lx - 1 || xPos < 0 || yPos > ly - 1 || yPos < 0 || zPos > lz - 1 || zPos < 0) {
      cerr << "ERROR: wrong array position!" << endl;
      cerr << "pos=" << xPos << ", " << yPos << ", " << zPos << endl;
      cerr << "coorda = (" << a_coordinate(i, 2) << ", " << a_coordinate(i, 2) << ", " << a_coordinate(i, 2) << ")"
           << endl;
      cerr << "actuallength=" << actualCellLength << ", smallestlength=" << m_smallestCellLength << endl;
      cerr << "minlevel=" << minLevel() << ", maxLevel=" << maxLevel() << endl;
      cerr << "lenght on level0=" << this->c_cellLengthAtLevel(0) << endl;
      TERMM(1, "Wrong array position");
    }

    // add normalized real parts of fluctuations
    u[0] += uPhysField[zPos + lz * (yPos + ly * xPos)][0];
    u[1] += vPhysField[zPos + lz * (yPos + ly * xPos)][0];
    u[2] += wPhysField[zPos + lz * (yPos + ly * xPos)][0];

    // // add normalized real parts of fluctuations factorized by u_mean/u_max
    // u[0] += uPhysField[zPos+lz*(yPos+ly*xPos)][0] * u[0]/(uTau * (C3*log(m_ReTau)+C4));
    // u[1] += vPhysField[zPos+lz*(yPos+ly*xPos)][0] * u[0]/(uTau * (C3*log(m_ReTau)+C4));
    // u[2] += wPhysField[zPos+lz*(yPos+ly*xPos)][0] * u[0]/(uTau * (C3*log(m_ReTau)+C4));

    a_variable(i, PV->RHO) = rho;
    a_oldVariable(i, PV->RHO) = rho;
    //    cells[i].m_oldVariablesT1[PV->RHO] = rho;
    a_variable(i, PV->U) = u[0];
    a_oldVariable(i, PV->U) = u[0];
    //    cells[i].m_oldVariablesT1[PV->U] = u[0];
    a_variable(i, PV->V) = u[1];
    a_oldVariable(i, PV->V) = u[1];
    //    cells[i].m_oldVariablesT1[PV->V] = u[1];
    a_variable(i, PV->W) = u[2];
    a_oldVariable(i, PV->W) = u[2];
    //    cells[i].m_oldVariablesT1[PV->W] = u[2];
  }

  // initialize distr. functions
  initEqDistFunctions();

  // Free the FFTW memory
  fftw_free(uPhysField);
  fftw_free(vPhysField);
  fftw_free(wPhysField);
}
template <>
void LbSolverDxQy<2, 9, maia::lb::LbSysEqnIncompressible<2, 9>>::initLatticeBgkFftChannel() {
  TERMM(1, "Init method only available for D3Qy, yet!");
}

/** \fn LbSolverDxQy::initLatticeBgkFftPipe()
 * \author Andreas Lintermann
 * \date 24.01.2011
 * \propVal{initMethod,LB_TURBULENT_PIPE_INIT}
 *
 * Initializes standard Lattice BGK Turbulent Pipe with FFT initialization.
 * Set property "initMethod" to "LB_TURBULENT_PIPE_INIT" and the property "FFTInit" to "1" to use this function.
 *
 * \LBIC{LbSolverDxQy::initLatticeBgkFftPipe(), initLatticeBgkFftPipe, FftPipeInit}
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbSolverDxQy<nDim, nDist, SysEqn>::initLatticeBgkFftPipe() {
  TRACE();
  MInt xPos, yPos, zPos;
  MInt lx = 0, ly = 0, lz = 0;

  // fft-domain dimensions
  // this holds the size of the domain in number of cells on lowest level
  lx = m_arraySize[0] / FPOW2(maxLevel() - minLevel());
  ly = m_arraySize[1] / FPOW2(maxLevel() - minLevel());
  if(nDim == 3) lz = m_arraySize[2] / FPOW2(maxLevel() - minLevel());

  fftw_complex *uPhysField, *vPhysField, *wPhysField;

  // field of velocities from positve frequencies
  uPhysField = (fftw_complex*)fftw_malloc(lx * ly * lz * sizeof(fftw_complex));
  vPhysField = (fftw_complex*)fftw_malloc(lx * ly * lz * sizeof(fftw_complex));
  wPhysField = (fftw_complex*)fftw_malloc(lx * ly * lz * sizeof(fftw_complex));

  // in multisolver mode the velocity field is created in solver 0 and broadcasted to all others
  ScratchSpace<MFloat> sendRecvBufferU(lx * ly * lz, AT_, "sendRecvBufferU");
  ScratchSpace<MFloat> sendRecvBufferV(lx * ly * lz, AT_, "sendRecvBufferV");
  ScratchSpace<MFloat> sendRecvBufferW(lx * ly * lz, AT_, "sendRecvBufferW");


  if(domainId() == 0) {
    // create velocity field
    maia::math::initFft(uPhysField, vPhysField, wPhysField, lx, ly, lz, m_noPeakModes, m_Ma);

    // copy values to buffer
    for(MInt p = 0; p < lx; p++) {
      for(MInt q = 0; q < ly; q++) {
        for(MInt r = 0; r < lz; r++) {
          sendRecvBufferU[r + lz * (q + ly * p)] = uPhysField[r + lz * (q + ly * p)][0];
          sendRecvBufferV[r + lz * (q + ly * p)] = vPhysField[r + lz * (q + ly * p)][0];
          sendRecvBufferW[r + lz * (q + ly * p)] = wPhysField[r + lz * (q + ly * p)][0];
        }
      }
    }
  }

  MPI_Bcast(&sendRecvBufferU[0], lx * ly * lz, MPI_DOUBLE, 0, mpiComm(), AT_, "sendRecvBufferU[0]");
  MPI_Bcast(&sendRecvBufferV[0], lx * ly * lz, MPI_DOUBLE, 0, mpiComm(), AT_, "sendRecvBufferV[0]");
  MPI_Bcast(&sendRecvBufferW[0], lx * ly * lz, MPI_DOUBLE, 0, mpiComm(), AT_, "sendRecvBufferW[0]");

  if(domainId() != 0) {
    // copy values from buffer
    for(MInt p = 0; p < lx; p++) {
      for(MInt q = 0; q < ly; q++) {
        for(MInt r = 0; r < lz; r++) {
          uPhysField[r + lz * (q + ly * p)][0] = sendRecvBufferU[r + lz * (q + ly * p)];
          vPhysField[r + lz * (q + ly * p)][0] = sendRecvBufferV[r + lz * (q + ly * p)];
          wPhysField[r + lz * (q + ly * p)][0] = sendRecvBufferW[r + lz * (q + ly * p)];
        }
      }
    }
  }

  MFloat actualCellLength;

  MFloat rho, u[3] = {0.0, 0.0, 0.0}, uTau, yPlus;
  MFloat radius = 0.0;

  // all cells have the same density
  if(m_densityFluctuations)
    rho = 0.0;
  else
    rho = 1.0;

  uTau = (MFloat)m_ReTau * m_nu / m_referenceLength;

  // tmpwidth equals pipe-radius (in cells) * cellLength on maxLevel
  MFloat tmpWidth = F1B2 * m_referenceLength * m_smallestCellLength;

  cerr << " --- prescribing mean velocity profile for turbulent pipe ---" << endl;
  for(MInt i = 0; i < a_noCells(); i++) {
    initNu(i, m_nu);

    radius = sqrt(a_coordinate(i, 0) * a_coordinate(i, 0) + a_coordinate(i, 1) * a_coordinate(i, 1));

    yPlus = m_ReTau * (1.0 - radius / tmpWidth);

    u[0] = 0.0;
    u[1] = 0.0;

    if(yPlus <= 5.0) {
      u[2] = uTau * yPlus;
    } else {
      if(yPlus <= 30 && yPlus > 5.0) {
        u[2] = uTau * (C1 * log(yPlus) + C2);
      } else {
        if(yPlus > 30) {
          u[2] = uTau * (C3 * log(yPlus) + C4);
        }
      }
    }

    actualCellLength = this->grid().cellLengthAtCell(i);

    xPos = floor(
        (F1B2 * lx + (a_coordinate(i, 0) - F1B2 * (actualCellLength)) / (this->c_cellLengthAtLevel(minLevel()))) + 0.1);
    yPos = floor(
        (F1B2 * ly + (a_coordinate(i, 1) - F1B2 * (actualCellLength)) / (this->c_cellLengthAtLevel(minLevel()))) + 0.1);

    // zero z at the inlet
    zPos = floor(((a_coordinate(i, 2) - F1B2 * (actualCellLength)) / (this->c_cellLengthAtLevel(minLevel()))) + 0.1);

    if(xPos > lx - 1 || xPos < 0 || yPos > ly - 1 || yPos < 0 || zPos > lz - 1 || zPos < 0) {
      stringstream errorMessage;
      errorMessage << "ERROR: wrong array position!" << endl;
      errorMessage << "pos=" << xPos << ", " << yPos << ", " << zPos << endl;
      errorMessage << "coorda = (" << a_coordinate(i, 2) << ", " << a_coordinate(i, 2) << ", " << a_coordinate(i, 2)
                   << ")" << endl;
      errorMessage << "actuallength=" << actualCellLength << ", smallestlength=" << m_smallestCellLength << endl;
      errorMessage << "minlevel=" << minLevel() << ", maxLevel=" << maxLevel() << endl;
      errorMessage << "length on level0=" << this->c_cellLengthAtLevel(0) << endl;
      TERMM(1, errorMessage.str());
    }


    //     m_log <<"pos="<<xPos<<", "<<yPos<<", "<<zPos<<endl;
    //     m_log <<"uvw="<<uPhysField[zPos+lz*(yPos+ly*xPos)][0]<<", "<<vPhysField[zPos+lz*(yPos+ly*xPos)][0]<<",
    //     "<<wPhysField[zPos+lz*(yPos+ly*xPos)][0]<<endl;

    // add normalized real parts of fluctuations factorized by u_mean/u_max
    u[0] += 2.0 * uPhysField[zPos + lz * (yPos + ly * xPos)][0] * u[2] / (uTau * (C3 * log(m_ReTau) + C4));
    u[1] += 2.0 * vPhysField[zPos + lz * (yPos + ly * xPos)][0] * u[2] / (uTau * (C3 * log(m_ReTau) + C4));
    u[2] += 2.0 * wPhysField[zPos + lz * (yPos + ly * xPos)][0] * u[2] / (uTau * (C3 * log(m_ReTau) + C4));

    a_variable(i, PV->RHO) = rho;
    a_oldVariable(i, PV->RHO) = rho;
    //    cells[i].m_oldVariablesT1[PV->RHO] = rho;
    a_variable(i, PV->U) = u[0];
    a_oldVariable(i, PV->U) = u[0];
    //    cells[i].m_oldVariablesT1[PV->U] = u[0];
    a_variable(i, PV->V) = u[1];
    a_oldVariable(i, PV->V) = u[1];
    //    cells[i].m_oldVariablesT1[PV->V] = u[1];
    a_variable(i, PV->W) = u[2];
    a_oldVariable(i, PV->W) = u[2];
    //    cells[i].m_oldVariablesT1[PV->W] = u[2];
  }

  // initialize distr. functions
  initEqDistFunctions();

  // Free the FftW memory
  fftw_free(uPhysField);
  fftw_free(vPhysField);
  fftw_free(wPhysField);
}
template <>
void LbSolverDxQy<2, 9, maia::lb::LbSysEqnIncompressible<2, 9>>::initLatticeBgkFftPipe() {
  TERMM(1, "Init method only available for D3Qy, yet!");
}

/** \fn LbSolverDxQy::initLatticeBgkFftMixing()
 * \author Georg
 * \date 08.12.2011
 * \propVal{initMethod,LB_TURBULENT_MIXING_FILTER_INIT}
 *
 * Initializes standard Lattice BGK Turbulent Boundary with FFT initialization.
 * Set property "initMethod" to "LB_TURBULENT_MIXING_INIT" and the property "FFTInit" to "1" to use this
 *function.
 *
 * \LBIC{LbSolverDxQy::initLatticeBgkFftMixing(), initLatticeBgkFftMixing, FftMixingInit}
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbSolverDxQy<nDim, nDist, SysEqn>::initLatticeBgkFftMixing() {
  TRACE();
  MInt xPos, yPos, zPos;
  MInt lx = 0, ly = 0, lz = 0;

  MFloat actualCellLength;
  MFloat GaussFactor, sigma, deltaOmega0;

  MFloat rho = 1.0, u[3] = {0.0, 0.0, 0.0};

  // fft-domain dimensions
  // this holds the size of the domain in number of cells on lowest level
  lx = m_arraySize[0] / FPOW2(maxLevel() - minLevel());
  ly = m_arraySize[1] / FPOW2(maxLevel() - minLevel());
  if(nDim == 3) lz = m_arraySize[2] / FPOW2(maxLevel() - minLevel());

  fftw_complex *uPhysField, *vPhysField, *wPhysField;

  // field of velocities from positve frequencies (is deleted at the end of the method)
  uPhysField = (fftw_complex*)fftw_malloc(lx * ly * lz * sizeof(fftw_complex));
  vPhysField = (fftw_complex*)fftw_malloc(lx * ly * lz * sizeof(fftw_complex));
  wPhysField = (fftw_complex*)fftw_malloc(lx * ly * lz * sizeof(fftw_complex));

  // in multisolver mode the velocity field is created in solver 0 and broadcasted to all others
  ScratchSpace<MFloat> sendRecvBufferU(lx * ly * lz, AT_, "sendRecvBufferU");
  ScratchSpace<MFloat> sendRecvBufferV(lx * ly * lz, AT_, "sendRecvBufferV");
  ScratchSpace<MFloat> sendRecvBufferW(lx * ly * lz, AT_, "sendRecvBufferW");
  if(domainId() == 0) {
    // create velocity field
    maia::math::initFft(uPhysField, vPhysField, wPhysField, lx, ly, lz, m_noPeakModes, m_Ma);

    // copy values to buffer
    for(MInt p = 0; p < lx; p++) {
      for(MInt q = 0; q < ly; q++) {
        for(MInt r = 0; r < lz; r++) {
          sendRecvBufferU[r + lz * (q + ly * p)] = uPhysField[r + lz * (q + ly * p)][0];
          sendRecvBufferV[r + lz * (q + ly * p)] = vPhysField[r + lz * (q + ly * p)][0];
          sendRecvBufferW[r + lz * (q + ly * p)] = wPhysField[r + lz * (q + ly * p)][0];
        }
      }
    }
  }

  MPI_Bcast(&sendRecvBufferU[0], lx * ly * lz, MPI_DOUBLE, 0, mpiComm(), AT_, "sendRecvBufferU[0]");
  MPI_Bcast(&sendRecvBufferV[0], lx * ly * lz, MPI_DOUBLE, 0, mpiComm(), AT_, "sendRecvBufferV[0]");
  MPI_Bcast(&sendRecvBufferW[0], lx * ly * lz, MPI_DOUBLE, 0, mpiComm(), AT_, "sendRecvBufferW[0]");

  if(domainId() != 0) {
    // copy values from buffer
    for(MInt p = 0; p < lx; p++) {
      for(MInt q = 0; q < ly; q++) {
        for(MInt r = 0; r < lz; r++) {
          uPhysField[r + lz * (q + ly * p)][0] = sendRecvBufferU[r + lz * (q + ly * p)];
          vPhysField[r + lz * (q + ly * p)][0] = sendRecvBufferV[r + lz * (q + ly * p)];
          wPhysField[r + lz * (q + ly * p)][0] = sendRecvBufferW[r + lz * (q + ly * p)];
        }
      }
    }
  }


  ScratchSpace<MFloat> bBox(6, AT_, "bBox");
  MFloat* bBoxPtr = &bBox[0];

  m_geometry->getBoundingBox(bBoxPtr);

  // all cells have the same density
  if(m_densityFluctuations)
    rho = 0.0;
  else
    rho = 1.0;

  // initial vorticity thickness [cells]
  deltaOmega0 = m_referenceLength * m_smallestCellLength;

  // width of the Gaussian distribution for perturbation scaling
  // depends on initial vorticity thickness
  sigma = F1B2 * deltaOmega0;

  m_log << " initial macroscopic vorticity thickness: " << deltaOmega0 << endl;

  cerr << " --- prescribing mean velocity profile for turbulent mixing layer ---" << endl;
  for(MInt i = 0; i < a_noCells(); i++) {
    initNu(i, m_nu);

    GaussFactor = exp(-(a_coordinate(i, 1) / sigma) * (a_coordinate(i, 1) / sigma) * F1B2);

    u[0] = m_Ma * LBCS * tanh(2.0 * a_coordinate(i, 1) / deltaOmega0);
    u[1] = 0.0;
    u[2] = 0.0;

    actualCellLength = this->grid().cellLengthAtCell(i);

    xPos = floor(
        (F1B2 * lx + (a_coordinate(i, 0) - F1B2 * (actualCellLength)) / (this->c_cellLengthAtLevel(minLevel()))) + 0.1);
    yPos = floor(
        (F1B2 * ly + (a_coordinate(i, 1) - F1B2 * (actualCellLength)) / (this->c_cellLengthAtLevel(minLevel()))) + 0.1);

    zPos = floor(
        (F1B2 * lz + (a_coordinate(i, 2) - F1B2 * (actualCellLength)) / (this->c_cellLengthAtLevel(minLevel()))) + 0.1);

    if(xPos > lx - 1 || xPos < 0 || yPos > ly - 1 || yPos < 0 || zPos > lz - 1 || zPos < 0) {
      stringstream errorMessage;
      errorMessage << "ERROR: wrong array position!" << endl;
      errorMessage << "pos=" << xPos << ", " << yPos << ", " << zPos << endl;
      errorMessage << "coorda = (" << a_coordinate(i, 2) << ", " << a_coordinate(i, 2) << ", " << a_coordinate(i, 2)
                   << ")" << endl;
      errorMessage << "actuallength=" << actualCellLength << ", smallestlength=" << m_smallestCellLength << endl;
      errorMessage << "minlevel=" << minLevel() << ", maxLevel=" << maxLevel() << endl;
      errorMessage << "lenght on level0=" << this->c_cellLengthAtLevel(0) << endl;
      TERMM(1, errorMessage.str());
    }

    // m_log <<"pos="<<xPos<<", "<<yPos<<", "<<zPos<<endl;
    // m_log <<"uvw="<<uPhysField[zPos+lz*(yPos+ly*xPos)][0]<<", "<<vPhysField[zPos+lz*(yPos+ly*xPos)][0]<<",
    // "<<wPhysField[zPos+lz*(yPos+ly*xPos)][0]<<endl;

    // add normalized real parts of fluctuations factorized by Gaussian
    u[0] += uPhysField[zPos + lz * (yPos + ly * xPos)][0] * GaussFactor;
    u[1] += vPhysField[zPos + lz * (yPos + ly * xPos)][0] * GaussFactor;
    u[2] += wPhysField[zPos + lz * (yPos + ly * xPos)][0] * GaussFactor;

    a_variable(i, PV->RHO) = rho;
    a_variable(i, PV->U) = u[0];
    a_variable(i, PV->V) = u[1];
    a_variable(i, PV->W) = u[2];

    a_oldVariable(i, PV->RHO) = rho;
    //    cells[i].m_oldVariablesT1[PV->RHO] = rho;

    a_oldVariable(i, PV->U) = u[0];
    //    cells[i].m_oldVariablesT1[PV->U] = u[0];

    a_oldVariable(i, PV->V) = u[1];
    //    cells[i].m_oldVariablesT1[PV->V] = u[1];

    a_oldVariable(i, PV->W) = u[2];
    //  cells[i].m_oldVariablesT1[PV->W] = u[2];
  }

  // takeCurl();

  // initialize distr. functions
  initNonEqDistFunctions();

  // Free the FFTW memory
  fftw_free(uPhysField);
  fftw_free(vPhysField);
  fftw_free(wPhysField);
}
template <>
void LbSolverDxQy<2, 9, maia::lb::LbSysEqnIncompressible<2, 9>>::initLatticeBgkFftMixing() {
  TERMM(1, "Init method only available for D3Qy, yet!");
}

/** \fn LbSolverDxQy::initLatticeBgkFftMixingFilter()
 * \author Georg
 * \date 08.12.2011
 * \propVal{initMethod,LB_TURBULENT_MIXING_FILTER_INIT}
 *
 * Initializes standard Lattice BGK Turbulent Boundary with FFT initialization.
 * The disturbances are filtered to a resolution 4 times lower than the resolution of the disturbance field.
 * Set property "initMethod" to "LB_TURBULENT_MIXING_FILTER_INIT" and the property "FFTInit" to "1" to use this
 * function.
 *
 * \LBIC{LbSolverDxQy::initLatticeBgkFftMixingFilter(), initLatticeBgkFftMixingFilter, FftMixingFilterInit}
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbSolverDxQy<nDim, nDist, SysEqn>::initLatticeBgkFftMixingFilter() {
  TRACE();

  MInt xPos, yPos, zPos;
  MInt lx = 0, ly = 0, lz = 0;
  MInt lxCoarse = 0, lyCoarse = 0, lzCoarse = 0;

  MFloat actualCellLength;
  MFloat GaussFactor, sigma, deltaOmega0;

  MFloat rho, u[3] = {0.0, 0.0, 0.0};

  // all cells have the same density
  if(m_densityFluctuations)
    rho = 0.0;
  else
    rho = 1.0;

  // fft-domain dimensions
  // this holds the size of the domain in number of cells on highest level
  lx = m_arraySize[0];
  ly = m_arraySize[1];
  if(nDim == 3) lz = m_arraySize[2];

  // size of the filtered domain
  lxCoarse = lx / 2;
  lyCoarse = ly / 2;
  lzCoarse = lz / 2;

  // field of velocities filtered to coarse grid (is deleted at the end of the method)
  ScratchSpace<MFloat> uPhysFieldCoarse(lxCoarse * lyCoarse * lzCoarse, AT_, "uPhysFieldCoarse");
  ScratchSpace<MFloat> vPhysFieldCoarse(lxCoarse * lyCoarse * lzCoarse, AT_, "vPhysFieldCoarse");
  ScratchSpace<MFloat> wPhysFieldCoarse(lxCoarse * lyCoarse * lzCoarse, AT_, "wPhysFieldCoarse");

  // in multisolver mode the velocity field is created in solver 0 and broadcasted to all others
  ScratchSpace<MFloat> sendRecvBufferU(lxCoarse * lyCoarse * lzCoarse, AT_, "sendRecvBufferU");
  ScratchSpace<MFloat> sendRecvBufferV(lxCoarse * lyCoarse * lzCoarse, AT_, "sendRecvBufferV");
  ScratchSpace<MFloat> sendRecvBufferW(lxCoarse * lyCoarse * lzCoarse, AT_, "sendRecvBufferW");

  if(domainId() == 0) {
    // create velocity field
    MFloat* uPhysFieldCoarsePtr = &uPhysFieldCoarse[0];
    MFloat* vPhysFieldCoarsePtr = &vPhysFieldCoarse[0];
    MFloat* wPhysFieldCoarsePtr = &wPhysFieldCoarse[0];
    maia::math::initFftFilter(uPhysFieldCoarsePtr, vPhysFieldCoarsePtr, wPhysFieldCoarsePtr, lx, ly, lz, lxCoarse,
                              lyCoarse, lzCoarse, m_noPeakModes, m_Ma);

    // copy values to buffer
    for(MInt p = 0; p < lxCoarse; p++) {
      for(MInt q = 0; q < lyCoarse; q++) {
        for(MInt r = 0; r < lzCoarse; r++) {
          sendRecvBufferU[r + lzCoarse * (q + lyCoarse * p)] = uPhysFieldCoarse[r + lzCoarse * (q + lyCoarse * p)];
          sendRecvBufferV[r + lzCoarse * (q + lyCoarse * p)] = vPhysFieldCoarse[r + lzCoarse * (q + lyCoarse * p)];
          sendRecvBufferW[r + lzCoarse * (q + lyCoarse * p)] = wPhysFieldCoarse[r + lzCoarse * (q + lyCoarse * p)];
        }
      }
    }
  }

  MPI_Bcast(&sendRecvBufferU[0], lxCoarse * lyCoarse * lzCoarse, MPI_DOUBLE, 0, mpiComm(), AT_, "sendRecvBufferU[0]");
  MPI_Bcast(&sendRecvBufferV[0], lxCoarse * lyCoarse * lzCoarse, MPI_DOUBLE, 0, mpiComm(), AT_, "sendRecvBufferV[0]");
  MPI_Bcast(&sendRecvBufferW[0], lxCoarse * lyCoarse * lzCoarse, MPI_DOUBLE, 0, mpiComm(), AT_, "sendRecvBufferW[0]");

  if(domainId() != 0) {
    // copy values from buffer
    for(MInt p = 0; p < lxCoarse; p++) {
      for(MInt q = 0; q < lyCoarse; q++) {
        for(MInt r = 0; r < lzCoarse; r++) {
          uPhysFieldCoarse[r + lzCoarse * (q + lyCoarse * p)] = sendRecvBufferU[r + lzCoarse * (q + lyCoarse * p)];
          vPhysFieldCoarse[r + lzCoarse * (q + lyCoarse * p)] = sendRecvBufferV[r + lzCoarse * (q + lyCoarse * p)];
          wPhysFieldCoarse[r + lzCoarse * (q + lyCoarse * p)] = sendRecvBufferW[r + lzCoarse * (q + lyCoarse * p)];
        }
      }
    }
  }

  ScratchSpace<MFloat> bBox(6, AT_, "bBox");
  MFloat* bBoxPtr = &bBox[0];
  m_geometry->getBoundingBox(bBoxPtr);

  // initial vorticity thickness [cells]
  deltaOmega0 = m_referenceLength * m_smallestCellLength;

  // width of the Gaussian distribution for perturbation scaling
  // depends on initial vorticity thickness
  sigma = deltaOmega0;

  m_log << " initial macroscopic vorticity thickness: " << deltaOmega0 << endl;

  cerr << " --- prescribing mean velocity profile for turbulent mixing layer ---" << endl;
  for(MInt i = 0; i < a_noCells(); i++) {
    initNu(i, m_nu);

    GaussFactor = exp(-(a_coordinate(i, 1) / sigma) * (a_coordinate(i, 1) / sigma) * F1B2);

    u[0] = m_Ma * LBCS * tanh(2.0 * a_coordinate(i, 1) / deltaOmega0);
    u[1] = 0.0;
    u[2] = 0.0;

    actualCellLength = this->grid().cellLengthAtCell(i);

    xPos = floor(
        (F1B2 * lxCoarse + (a_coordinate(i, 0) - F1B2 * (actualCellLength)) / (this->c_cellLengthAtLevel(maxLevel())))
        + 0.1);
    yPos = floor(
        (F1B2 * lyCoarse + (a_coordinate(i, 1) - F1B2 * (actualCellLength)) / (this->c_cellLengthAtLevel(maxLevel())))
        + 0.1);

    zPos = floor(
        (F1B2 * lzCoarse + (a_coordinate(i, 2) - F1B2 * (actualCellLength)) / (this->c_cellLengthAtLevel(maxLevel())))
        + 0.1);

    if(xPos > lxCoarse - 1 || xPos < 0 || yPos > lyCoarse - 1 || yPos < 0 || zPos > lzCoarse - 1 || zPos < 0) {
      cerr << "ERROR: wrong array position!" << endl;
      cerr << "pos=" << xPos << ", " << yPos << ", " << zPos << endl;
      cerr << "coorda = (" << a_coordinate(i, 2) << ", " << a_coordinate(i, 2) << ", " << a_coordinate(i, 2) << ")"
           << endl;
      cerr << "actuallength=" << actualCellLength << ", smallestlength=" << m_smallestCellLength << endl;
      cerr << "minlevel=" << minLevel() << ", maxLevel=" << maxLevel() << endl;
      cerr << "lenght on level0=" << this->c_cellLengthAtLevel(0) << endl;
      TERMM(1, "Wrong array position");
    }

    // m_log <<"pos="<<xPos<<", "<<yPos<<", "<<zPos<<endl;
    // m_log <<"uvw="<<uPhysField[zPos+lz*(yPos+ly*xPos)][0]<<", "<<vPhysField[zPos+lz*(yPos+ly*xPos)][0]<<",
    // "<<wPhysField[zPos+lz*(yPos+ly*xPos)][0]<<endl;

    // add normalized real parts of fluctuations factorized by Gaussian
    u[0] += uPhysFieldCoarse[zPos + lzCoarse * (yPos + lyCoarse * xPos)] * GaussFactor;
    u[1] += vPhysFieldCoarse[zPos + lzCoarse * (yPos + lyCoarse * xPos)] * GaussFactor;
    u[2] += wPhysFieldCoarse[zPos + lzCoarse * (yPos + lyCoarse * xPos)] * GaussFactor;

    a_variable(i, PV->RHO) = rho;
    a_variable(i, PV->U) = u[0];
    a_variable(i, PV->V) = u[1];
    a_variable(i, PV->W) = u[2];

    a_oldVariable(i, PV->RHO) = rho;
    //    cells[i].m_oldVariablesT1[PV->RHO] = rho;

    a_oldVariable(i, PV->U) = u[0];
    //    cells[i].m_oldVariablesT1[PV->U] = u[0];

    a_oldVariable(i, PV->V) = u[1];
    //    cells[i].m_oldVariablesT1[PV->V] = u[1];

    a_oldVariable(i, PV->W) = u[2];
    //  cells[i].m_oldVariablesT1[PV->W] = u[2];
  }

  // takeCurl();

  // initialize distr. functions
  initNonEqDistFunctions();
  // initEqDistFunctions();
}
template <>
void LbSolverDxQy<2, 9, maia::lb::LbSysEqnIncompressible<2, 9>>::initLatticeBgkFftMixingFilter() {
  TERMM(1, "Init method only available for D3Qy, yet!");
}

/** \fn LbSolverDxQy::initLatticeBgkFftIsotropicTurbulence()
 * \author Johannes Grafen adapeted from version of Lennart Schneiders
 * \date 10.02.2022
 * \propVal{initMethod,LB_TURBULENCE_ISOTROPIC_INIT}
 *
 * Creates isotropic Turbulence in box shaped domain, adapted from FV
 * In comparison to FV, where the pressure field is calulated by solving a
 * poisson equation for the pressure in fourier-transformed wavespace,
 * an iterative initialization run for the initialization of the pressure/density field
 * is necessary. See function: initRunCorrection()
 *
 * \LBIC{LbSolverDxQy::initLatticeBgkFftIsotropicTurbulence(), initLatticeBgkFftIsotropicTurbulence,
 * FftIsotropicTurbulenceInit}
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbSolverDxQy<nDim, nDist, SysEqn>::initLatticeBgkFftIsotropicTurbulence() {
  TRACE();
  // adapted from FV:
  IF_CONSTEXPR(nDim == 2) TERMM(-1, "Only works for 3D!");
  const MFloat time0 = MPI_Wtime();
  const MInt fftLevel = grid().maxUniformRefinementLevel();
  const MFloat DX = c_cellLengthAtLevel(fftLevel);
  if(fftLevel > grid().maxUniformRefinementLevel()) {
    mTerm(1, AT_, "Isotropic mesh expected (0).");
  }
  MInt noVars = nDim + 1; // U, V, W, RHO

  std::array<MFloat, nDim * 2> bBox;
  m_geometry->getBoundingBox(bBox.data());

  /* if(Context::propertyExists("cutOffCoordinates", m_solverId)) {
    for(MInt i = 0; i < nDim; i++) {
      bBox[i] = numeric_limits<MInt>::max();
      bBox[nDim + i] = numeric_limits<MInt>::min();
    }
    for(MInt cellId = 0; cellId < a_noCells(); cellId++) {
      if(a_isHalo(cellId)) continue;
      //if(a_isPeriodic(cellId)) continue;
      if(a_isBndryGhostCell(cellId)) continue;
      for(MInt i = 0; i < nDim; i++) {
        bBox[i] = mMin(bBox[i], a_coordinate(cellId, i) - F1B2 * c_cellLengthAtCell(cellId));
        bBox[nDim + i] = mMax(bBox[nDim + i], a_coordinate(cellId, i) + F1B2 * c_cellLengthAtCell(cellId));
      }
    }
    MPI_Allreduce(MPI_IN_PLACE, &bBox[0], nDim, MPI_DOUBLE, MPI_MIN, mpiComm(), AT_, "MPI_IN_PLACE", "bBox[0]");
    MPI_Allreduce(MPI_IN_PLACE, &bBox[nDim], nDim, MPI_DOUBLE, MPI_MAX, mpiComm(), AT_, "MPI_IN_PLACE",
                  "bBox[nDim]");
  } */

  // fft-domain dimensions
  // this holds the size of the domain in number of cells on lowest level

  const MFloat dxeps = 0.1 * DX;
  MInt nx = (bBox[3] - bBox[0] + dxeps) / DX;
  MInt ny = (bBox[4] - bBox[1] + dxeps) / DX;
  MInt nz = (bBox[5] - bBox[2] + dxeps) / DX;
  const MInt size = nx * ny * nz;

  getReLambdaAndUrmsInit();
  if(domainId() == 0) cerr << "UrmsInit: " << m_UrmsInit << endl << "ReLambdaInit: " << m_ReLambda << endl;

  // consistency check
  const MFloat nuCheck = m_referenceLength * m_UrmsInit / m_Re;
  if(abs(nuCheck - m_nu) > 0.01 * m_nu) // 1% difference is tolerated between nuCheck and nuLb
    mTerm(1, AT_, "nu_LB and nu_check are not consistent! Check Re, Ma and UrmsInit");

  MFloat rhoInfinity = F1;
  MFloat UInfinity = m_UrmsInit; // equal to m_Ma * LBCS

  MBool goodSpectrum = false;
  MUlong seed = 0;

  MInt spectrumId = 2; // prescribed energy spectrum, see maiamath.cpp
  spectrumId = Context::getSolverProperty<MInt>("spectrumId", this->m_solverId, AT_, &spectrumId);
  MFloat kpRatio = F4; // peak wave number of prescribed spectrum
  kpRatio = Context::getSolverProperty<MFloat>("kpRatio", this->m_solverId, AT_, &kpRatio);

  while(!goodSpectrum) {
    fftw_complex* uPhysField;
    fftw_complex* nabla2P = nullptr;
    MIntScratchSpace fftInfo(4, AT_, "fftInfo");

    seed = maia::math::initFft(uPhysField, nabla2P, nx, ny, nz, kpRatio, spectrumId, fftInfo, mpiComm(),
                               false); // Using method of Lennart Schneiders, not Ma-number dependent

    if(domainId() == 0) {
      std::cerr << "Theoretical SPECTRUM: kinetic energy: " << 1.5 * POW2(m_UrmsInit)
                << ", dissipation rate * Box length / Urms,0^3: "
                << 72.0 * POW2(kpRatio * PI) * m_nu / (m_UrmsInit * nx) << ")"
                << ", integral length/unit length: " << F3B8 / kpRatio << ")" << std::endl;
    }

    MInt maxRank = fftInfo[0];
    MInt local_n0 = fftInfo[1];
    MInt local_0_start = fftInfo[2];
    MInt alloc_local = fftInfo[3];

    MIntScratchSpace noSendIds(globalNoDomains(), AT_, "noSendIds");
    MInt sendIdsSize = 0;
    MIntScratchSpace offsetsIds(globalNoDomains() + 1, AT_, "offsetsIds");
    MInt locOffsetIds = (globalDomainId() < maxRank) ? ((MInt)local_0_start) * ny * nz : size;
    MPI_Allgather(&locOffsetIds, 1, MPI_INT, &offsetsIds[0], 1, MPI_INT, mpiComm(), AT_, "locOffsetIds",
                  "offsetsIds[0]");
    offsetsIds(globalNoDomains()) = size;
    noSendIds.fill(0);

    // Count no of sendIds
    for(MInt cellId = 0; cellId < a_noCells(); cellId++) {
      if(a_isHalo(cellId)) continue;
      if(a_isBndryGhostCell(cellId)) continue;
      if(a_level(cellId) != fftLevel) continue;
      MFloat actualCellLength = grid().cellLengthAtLevel(a_level(cellId));
      MInt xPos = floor((F1B2 * nx + (a_coordinate(cellId, 0) - F1B2 * (actualCellLength)) / DX) + 0.1);
      MInt yPos = floor((F1B2 * ny + (a_coordinate(cellId, 1) - F1B2 * (actualCellLength)) / DX) + 0.1);
      MInt zPos = floor((F1B2 * nz + (a_coordinate(cellId, 2) - F1B2 * (actualCellLength)) / DX) + 0.1);
      MInt pos = maia::math::getGlobalPosFFTW(xPos, yPos, zPos, ny, nz);
      MInt nghbrDomain = mMin(maxRank - 1, pos / (size / maxRank));
      while(pos < offsetsIds(nghbrDomain) || pos >= offsetsIds(nghbrDomain + 1)) {
        if(pos < offsetsIds(nghbrDomain)) nghbrDomain--;
        if(pos >= offsetsIds(nghbrDomain + 1)) nghbrDomain++;
      }
      if(nghbrDomain >= maxRank) mTerm(1, AT_, "wrong domain");
      noSendIds(nghbrDomain)++;
      sendIdsSize++;
    }
    MIntScratchSpace sendIdsOffsets(globalNoDomains(), AT_, "sendIdsOffsets");
    MIntScratchSpace sendIdsOffsetsTmp(globalNoDomains(), AT_, "sendIdsOffsetsTmp");
    sendIdsOffsets[0] = 0;
    sendIdsOffsetsTmp[0] = 0;
    for(MInt nghbrDomain = 0; nghbrDomain < globalNoDomains() - 1; nghbrDomain++) {
      sendIdsOffsets[nghbrDomain + 1] = sendIdsOffsets[nghbrDomain] + noSendIds[nghbrDomain];
      sendIdsOffsetsTmp[nghbrDomain + 1] = sendIdsOffsets[nghbrDomain] + noSendIds[nghbrDomain];
    }
    MIntScratchSpace sendIds(sendIdsSize, AT_, "sendIdsSize");

    // Fill sendIds
    for(MInt cellId = 0; cellId < a_noCells(); cellId++) {
      if(a_isHalo(cellId)) continue;
      if(a_isBndryGhostCell(cellId)) continue;
      if(a_level(cellId) != fftLevel) continue;
      MFloat actualCellLength = grid().cellLengthAtLevel(a_level(cellId));
      MInt xPos = floor((F1B2 * nx + (a_coordinate(cellId, 0) - F1B2 * (actualCellLength)) / DX) + 0.1);
      MInt yPos = floor((F1B2 * ny + (a_coordinate(cellId, 1) - F1B2 * (actualCellLength)) / DX) + 0.1);
      MInt zPos = floor((F1B2 * nz + (a_coordinate(cellId, 2) - F1B2 * (actualCellLength)) / DX) + 0.1);
      MInt pos = maia::math::getGlobalPosFFTW(xPos, yPos, zPos, nx, ny);
      MInt nghbrDomain = mMin(maxRank - 1, pos / (size / maxRank));
      while(pos < offsetsIds(nghbrDomain) || pos >= offsetsIds(nghbrDomain + 1)) {
        if(pos < offsetsIds(nghbrDomain)) nghbrDomain--;
        if(pos >= offsetsIds(nghbrDomain + 1)) nghbrDomain++;
      }
      if(nghbrDomain >= maxRank) mTerm(1, AT_, "wrong domain");
      sendIds[sendIdsOffsetsTmp[nghbrDomain]++] = pos;
    }
    MIntScratchSpace noRecvIds(globalNoDomains(), AT_, "noRecvIds");
    noRecvIds.fill(0);
    MPI_Alltoall(&noSendIds[0], 1, MPI_INT, &noRecvIds[0], 1, MPI_INT, mpiComm(), AT_, "noSendIds[0]", "noRecvIds[0]");

    MIntScratchSpace recvIdsOffsets(globalNoDomains(), AT_, "recvIdsOffsets");
    recvIdsOffsets[0] = 0;
    MInt recvIdsSize = 0;
    for(MInt nghbrDomain = 0; nghbrDomain < globalNoDomains() - 1; nghbrDomain++) {
      recvIdsOffsets[nghbrDomain + 1] = recvIdsOffsets[nghbrDomain] + noRecvIds[nghbrDomain];
      recvIdsSize += noRecvIds[nghbrDomain];
    }
    recvIdsSize += noRecvIds[globalNoDomains() - 1];

    // Exchange Ids
    MIntScratchSpace recvIds(mMax(1, recvIdsSize), AT_, "recvIds");
    MPI_Alltoallv(&sendIds[0], &noSendIds[0], &sendIdsOffsets[0], MPI_INT, &recvIds[0], &noRecvIds[0],
                  &recvIdsOffsets[0], MPI_INT, mpiComm(), AT_, "sendIds[0]", "recvIds[0]");

    MIntScratchSpace sendVarsOffsets(globalNoDomains(), AT_, "sendVarsOffsets");
    MIntScratchSpace sendVarsOffsetsTmp(globalNoDomains(), AT_, "sendVarsOffsetsTmp");
    MIntScratchSpace noSendVars(globalNoDomains(), AT_, "noSendVars");
    sendVarsOffsets[0] = 0;
    sendVarsOffsetsTmp[0] = 0;
    MInt sendVarsSize = 0;
    for(MInt i = 0; i < globalNoDomains() - 1; i++) {
      sendVarsOffsets[i + 1] = sendVarsOffsets[i] + noRecvIds[i] * 3;
      sendVarsOffsetsTmp[i + 1] = sendVarsOffsetsTmp[i] + noRecvIds[i] * 3;
      noSendVars[i] = noRecvIds[i] * 3;
      sendVarsSize += noRecvIds[i] * 3;
    }
    if(recvIdsSize != 0) {
      noSendVars[globalNoDomains() - 1] = noRecvIds[globalNoDomains() - 1] * 3;
      sendVarsSize += noRecvIds[globalNoDomains() - 1] * 3;
    }

    // Fill sendVars: for each Id received 3 variables have to be send
    MFloatScratchSpace sendVars(mMax(1, sendVarsSize), AT_, "sendVars");
    sendVars.fill(0);
    MFloat cnt = F0;
    MFloat upr[3] = {F0, F0, F0};
    for(MInt i = 0; i < recvIdsSize; i++) {
      MInt pos = recvIds[i];
      ASSERT(pos > -1 && pos < size, "");
      if(!(pos >= ((MInt)local_0_start) * nx * ny && pos < ((MInt)(local_0_start + local_n0) * ny * nx))) {
        mTerm(1, AT_, "Position not available on this domain(1).");
      }
      MInt localPos = pos - (((MInt)local_0_start) * ny * nz);
      if(3 * localPos + 2 > alloc_local) {
        mTerm(1, AT_, "index exceeds array(1)");
      }
      sendVars[i * 3] = uPhysField[3 * localPos][0];
      sendVars[i * 3 + 1] = uPhysField[3 * localPos + 1][0];
      sendVars[i * 3 + 2] = uPhysField[3 * localPos + 2][0];

      MFloat u = uPhysField[3 * localPos][0];
      MFloat v = uPhysField[3 * localPos + 1][0];
      MFloat w = uPhysField[3 * localPos + 2][0];
      upr[0] += POW2(u);
      upr[1] += POW2(v);
      upr[2] += POW2(w);
      cnt++;
    }

    MPI_Allreduce(MPI_IN_PLACE, &cnt, 1, MPI_DOUBLE, MPI_SUM, mpiComm(), AT_, "MPI_IN_PLACE", "cnt");
    MPI_Allreduce(MPI_IN_PLACE, &upr, 3, MPI_DOUBLE, MPI_SUM, mpiComm(), AT_, "MPI_IN_PLACE", "upr");

    MFloat uprime0 = sqrt(F1B3 * (upr[0] + upr[1] + upr[2]) / cnt);

    fftw_free(uPhysField);

    MIntScratchSpace recvVarsOffsets(globalNoDomains(), AT_, "recvVarsOffsets");
    MIntScratchSpace noRecvVars(globalNoDomains(), AT_, "noRecvVars");
    recvVarsOffsets[0] = 0;
    MInt recvVarsSize = 0;
    for(MInt i = 0; i < globalNoDomains() - 1; i++) {
      recvVarsOffsets[i + 1] = recvVarsOffsets[i] + noSendIds[i] * 3;
      noRecvVars[i] = noSendIds[i] * 3;
      recvVarsSize += noSendIds[i] * 3;
    }
    recvVarsSize += noSendIds[globalNoDomains() - 1] * 3;
    noRecvVars[globalNoDomains() - 1] = noSendIds[globalNoDomains() - 1] * 3;

    // Exchange Variables
    MFloatScratchSpace recvVars(recvVarsSize, AT_, "recvVars");
    recvVars.fill(0);
    MPI_Alltoallv(&sendVars[0], &noSendVars[0], &sendVarsOffsets[0], MPI_DOUBLE, &recvVars[0], &noRecvVars[0],
                  &recvVarsOffsets[0], MPI_DOUBLE, mpiComm(), AT_, "sendVars[0]", "recvVars[0]");

    cnt = F0;
    upr[0] = 0;
    upr[1] = 0;
    upr[2] = 0;
    for(MInt cellId = 0; cellId < a_noCells(); cellId++) {
      if(a_isHalo(cellId)) continue;
      if(a_isBndryGhostCell(cellId)) continue;
      if(a_level(cellId) != fftLevel) continue;
      MFloat actualCellLength = grid().cellLengthAtLevel(a_level(cellId));
      MInt xPos = floor((F1B2 * nx + (a_coordinate(cellId, 0) - F1B2 * (actualCellLength)) / DX) + 0.1);
      MInt yPos = floor((F1B2 * ny + (a_coordinate(cellId, 1) - F1B2 * (actualCellLength)) / DX) + 0.1);
      MInt zPos = floor((F1B2 * nz + (a_coordinate(cellId, 2) - F1B2 * (actualCellLength)) / DX) + 0.1);
      MInt pos = maia::math::getGlobalPosFFTW(xPos, yPos, zPos, nx, ny);
      MInt nghbrDomain = mMin(maxRank - 1, pos / (size / maxRank));
      while(pos < offsetsIds(nghbrDomain) || pos >= offsetsIds(nghbrDomain + 1)) {
        if(pos < offsetsIds(nghbrDomain)) nghbrDomain--;
        if(pos >= offsetsIds(nghbrDomain + 1)) nghbrDomain++;
      }
      if(nghbrDomain >= maxRank) mTerm(1, AT_, "wrong domain");
      if(sendIds[sendIdsOffsets[nghbrDomain]] != pos) mTerm(1, AT_, "pos mismatch");
      MFloat u = recvVars[sendIdsOffsets[nghbrDomain] * 3];
      MFloat v = recvVars[sendIdsOffsets[nghbrDomain] * 3 + 1];
      MFloat w = recvVars[sendIdsOffsets[nghbrDomain] * 3 + 2];

      sendIdsOffsets[nghbrDomain]++;
      upr[0] += POW2(u);
      upr[1] += POW2(v);
      upr[2] += POW2(w);
      a_variable(cellId, PV->U) = u;
      a_variable(cellId, PV->V) = v;
      a_variable(cellId, PV->W) = w;
      a_variable(cellId, PV->RHO) = rhoInfinity; // is going to be iterated
      cnt++;
    }


    MPI_Allreduce(MPI_IN_PLACE, &cnt, 1, MPI_DOUBLE, MPI_SUM, mpiComm(), AT_, "MPI_IN_PLACE", "cnt");
    MPI_Allreduce(MPI_IN_PLACE, &upr, 3, MPI_DOUBLE, MPI_SUM, mpiComm(), AT_, "MPI_IN_PLACE", "upr");

    MFloat uprime1 = sqrt(F1B3 * (upr[0] + upr[1] + upr[2]) / cnt);
    upr[0] = sqrt(upr[0] / cnt);
    upr[1] = sqrt(upr[1] / cnt);
    upr[2] = sqrt(upr[2] / cnt);

    if(fabs(uprime0 - uprime1) > exp(-10)) mTerm(1, AT_, "Communication went wrong.");

    if(domainId() == 0)
      cerr << "initial urpime: " << upr[0] << " " << upr[1] << " " << upr[2] << " (" << uprime1 << ")" << endl;

    // for LES: The initial spectrum is cut off for coarse LES grids which
    // results into a lower energy level. If uprime != 1, this energy
    // will be modified such that all resolved scales have the DNS-energy,
    // i.e. large scales get the energy which have been cut off at the small
    // scales.

    // uprime1 = F1;

    for(MInt cellId = 0; cellId < a_noCells(); cellId++) {
      if(a_isHalo(cellId)) continue;
      if(a_isBndryGhostCell(cellId)) continue;
      if(a_level(cellId) != fftLevel) continue;
      // scale velocities such that mean urms is UInfinity, scaling with velocity magnitude of UrmsInit
      MFloat u = a_variable(cellId, PV->U) * UInfinity / uprime1;
      MFloat v = a_variable(cellId, PV->V) * UInfinity / uprime1;
      MFloat w = a_variable(cellId, PV->W) * UInfinity / uprime1;

      a_variable(cellId, PV->U) = u;
      a_variable(cellId, PV->V) = v;
      a_variable(cellId, PV->W) = w;
    }

    exchangeData(&a_variable(0, 0), noVars);

    MFloat umean[3] = {F0, F0, F0};
    MFloat urms[6] = {F0, F0, F0, F0, F0, F0};
    MFloat reyn[6] = {F0, F0, F0, F0, F0, F0};
    MFloat aniso[6] = {F0, F0, F0, F0, F0, F0};
    MFloat dudx[3] = {F0, F0, F0};
    MFloat skew[4] = {F0, F0, F0, F0};
    for(MInt cellId = 0; cellId < a_noCells(); cellId++) {
      if(a_isHalo(cellId)) continue;
      if(a_isBndryGhostCell(cellId)) continue;
      if(a_level(cellId) != fftLevel) continue;
      MFloat u = a_variable(cellId, PV->U);
      MFloat v = a_variable(cellId, PV->V);
      MFloat w = a_variable(cellId, PV->W);
      for(MInt i = 0; i < nDim; i++) {
        MInt n0 = (a_hasNeighbor(cellId, 2 * i) > 0) ? c_neighborId(cellId, 2 * i) : cellId;
        MInt n1 = (a_hasNeighbor(cellId, 2 * i + 1) > 0) ? c_neighborId(cellId, 2 * i + 1) : cellId;
        if(n0 == n1) continue;
        dudx[i] += POW2((a_variable(n1, PV->VV[i]) - a_variable(n0, PV->VV[i]))
                        / ((a_coordinate(n1, i) - a_coordinate(n0, i)) / DX));
        skew[i] += POW3((a_variable(n1, PV->VV[i]) - a_variable(n0, PV->VV[i]))
                        / ((a_coordinate(n1, i) - a_coordinate(n0, i)) / DX));
      }
      urms[0] += POW2(u);
      urms[1] += POW2(v);
      urms[2] += POW2(w);
      urms[3] += u * v;
      urms[4] += u * w;
      urms[5] += v * w;
      umean[0] += u;
      umean[1] += v;
      umean[2] += w;
    }
    MPI_Allreduce(MPI_IN_PLACE, &urms[0], 6, MPI_DOUBLE, MPI_SUM, mpiComm(), AT_, "MPI_IN_PLACE", "urms[0]");
    MPI_Allreduce(MPI_IN_PLACE, &umean[0], 3, MPI_DOUBLE, MPI_SUM, mpiComm(), AT_, "MPI_IN_PLACE", "umean[0]");
    MPI_Allreduce(MPI_IN_PLACE, &dudx[0], 3, MPI_DOUBLE, MPI_SUM, mpiComm(), AT_, "MPI_IN_PLACE", "dudx[0]");
    MPI_Allreduce(MPI_IN_PLACE, &skew[0], 3, MPI_DOUBLE, MPI_SUM, mpiComm(), AT_, "MPI_IN_PLACE", "skew[0]");
    skew[3] = ((skew[0] + skew[1] + skew[2]) / (F3 * cnt)) / pow((dudx[0] + dudx[1] + dudx[2]) / (F3 * cnt), 1.5);
    for(MInt i = 0; i < 3; i++)
      skew[i] = (skew[i] / cnt) / pow(dudx[i] / cnt, 1.5);
    for(MInt i = 0; i < 6; i++)
      reyn[i] = urms[i] / cnt;
    for(MInt i = 0; i < 6; i++)
      aniso[i] = urms[i] / (urms[0] + urms[1] + urms[2]);
    for(MInt i = 0; i < 3; i++)
      aniso[i] -= F1B3;
    for(MInt i = 0; i < 6; i++)
      urms[i] = sqrt(fabs(urms[i]) / cnt);
    for(MInt i = 0; i < 3; i++)
      umean[i] /= cnt;
    for(MInt i = 0; i < 3; i++)
      dudx[i] /= cnt;
    MFloat lambda[3];
    MFloat Rel[4];
    MFloat eps[3];
    for(MInt i = 0; i < 3; i++)
      lambda[i] = urms[i] / sqrt(dudx[i]); // Eq. 6.56 Turbulent Flows Pope two-point correlation

    if(domainId() == 0)
      cerr << "u_mean: " << umean[0] << " " << umean[1] << " " << umean[2] << " (" << (MInt)cnt << ")" << endl;
    if(domainId() == 0)
      cerr << "u_rms/u_inf: " << urms[0] / UInfinity << " " << urms[1] / UInfinity << " " << urms[2] / UInfinity << " "
           << urms[3] / UInfinity << " " << urms[4] / UInfinity << " " << urms[5] / UInfinity << " ("
           << F1B3 * (urms[0] + urms[1] + urms[2]) / UInfinity << ")" << endl;
    if(domainId() == 0)
      cerr << "skewness: " << skew[0] << " " << skew[1] << " " << skew[2] << " (" << skew[3] << ")" << endl;
    if(domainId() == 0)
      cerr << "Anisotropy: " << aniso[0] << " " << aniso[1] << " " << aniso[2] << " " << aniso[3] << " " << aniso[4]
           << " " << aniso[5] << endl;
    if(domainId() == 0)
      cerr << "Reynolds stress: " << reyn[0] << " " << reyn[1] << " " << reyn[2] << " " << reyn[3] << " " << reyn[4]
           << " " << reyn[5] << endl;
    if(domainId() == 0)
      cerr << "Realizability 1: " << reyn[3] * reyn[3] << " <= " << reyn[0] * reyn[1] << ", " << reyn[4] * reyn[4]
           << " <= " << reyn[0] * reyn[2] << ", " << reyn[5] * reyn[5] << " <= " << reyn[1] * reyn[2] << endl;
    if(domainId() == 0)
      cerr << "Realizability 2: det = "
           << reyn[0] * reyn[1] * reyn[2] + F2 * reyn[3] * reyn[4] * reyn[5] - reyn[0] * reyn[5] * reyn[5]
                  - reyn[1] * reyn[4] * reyn[4] - reyn[2] * reyn[4] * reyn[4]
           << " >= 0" << endl;

    // prescribe Taylor microscale Reynolds number
    // for LES: ReLambda is unkown in LES and has to be defined to get the
    // same viscosity as in the DNS.
    const MFloat time1 = MPI_Wtime();
    m_log << "initFFT time " << time1 - time0 << endl;

    /*if(nx < 64) {  // nx < 256 -> for LES, DNS requires refinement level of 8 or higher?
      ! \page propertyPage1
        \section ReLambdaOverwrite
        <code>MInt FvCartesianSolverXD::m_Re</code>\n
        default = <code>none</code>\n \n
        Overwrite Reynolds number for FV isotropic turbulence spectrum initial condition (case
        16). \n \n
        Possible values are:
        <ul>
          <li>Any positive floating point value</li>
        </ul>
        Keywords: <i>FINITE_VOLUME, INITIAL_CONDITION, ISOTROPIC, TURBULENCE</i>

      if(!Context::propertyExists("ReLambdaOverwrite", this->m_solverId)) {
        mTerm(1, AT_, "Undefined ReLambda for LES grid.");
      }
      m_Re = Context::getSolverProperty<MFloat>("ReLambdaOverwrite", this->m_solverId, AT_) / m_referenceLength;
      Re0 = m_Re * m_nu * rhoInfinity / (rhoInfinity * UInfinity);
      m_log << "Overwritten Reynolds number: " << setprecision(15) << m_Re
                << " (Re_L=" << m_Re * (bBox[3] - bBox[0]) << ")"
                << " " << Re0 << endl;
      break;
    }*/

    for(MInt i = 0; i < 3; i++)
      Rel[i] = urms[i] * lambda[i] / m_nu;
    Rel[3] = F1B3 * (urms[0] + urms[1] + urms[2]) * F1B3 * (lambda[0] + lambda[1] + lambda[2]) / m_nu;

    if(domainId() == 0) {
      cerr << "Re_lambda: " << Rel[0] << " " << Rel[1] << " " << Rel[2] << " (" << Rel[3] << ", "
           << F1B3 * (Rel[0] + Rel[1] + Rel[2]) << ")" << endl;

      const MFloat lambdaAvg = F1B3 * (lambda[0] + lambda[1] + lambda[2]);
      cerr << "lambda: " << lambda[0] << " " << lambda[1] << " " << lambda[2] << " (" << lambdaAvg << ") "
           << "lambda / Lb: " << lambdaAvg / m_referenceLength << endl;
    }

    MFloat eps2[3];
    for(MInt i = 0; i < 3; i++) {
      eps[i] = (15.0 * m_nu * dudx[i]) * m_referenceLength / POW3(urms[i]);
      eps2[i] = 15.0 * m_nu * POW2(urms[i] / lambda[i]); // Eq. 6.58 Turbulent Flows Pope two-point correlation
    }

    if(domainId() == 0) {
      cerr << "eps: " << eps[0] << " " << eps[1] << " " << eps[2] << " (" << F1B3 * (eps[0] + eps[1] + eps[2]) << ")"
           << endl;

      const MFloat eps2Avg = F1B3 * (eps2[0] + eps2[1] + eps2[2]);
      cerr << "eps2: " << eps2[0] << " " << eps2[1] << " " << eps2[2] << " (" << eps2Avg << ")"
           << "  eps2 * Box length / Urms,0^3: " << eps2Avg * m_referenceLength / POW3(m_UrmsInit) << endl;

      const MFloat eta =
          pow(m_nu, 0.75)
          / pow(F1B3 * (eps2[0] + eps2[1] + eps2[2]), 0.25); // Eq. 6.1 Turbulent Flows Pope two-point correlation
      cerr << "Kolmogorov length: " << eta << "  Kolmogorov length / Lb: " << eta / m_referenceLength << endl;
    }
    MFloat maxd = mMax(fabs(Rel[0] - Rel[1]), mMax(fabs(Rel[0] - Rel[2]), fabs(Rel[1] - Rel[2])));
    if(maxd < F1) m_log << "spectrum " << maxd << " " << F1B3 * (Rel[0] + Rel[1] + Rel[2]) << " " << seed << endl;
    // if ( maxd < 0.3 && F1B3*(Rel[0]+Rel[1]+Rel[2]) > 79.1 ) {
    // if ( maxd < 0.3 ) {
    goodSpectrum = true;
    if(domainId() == 0) m_log << "spectrum " << seed << endl;
    //}
    //}
  }
  if(domainId() == 0) cerr << "seed " << seed << endl;
  MPI_Barrier(mpiComm(), AT_);

  // TODO: grid().maxUniformRefinementLevel() has to be replaced with maxRefinementLevel() for adaptive meshes with
  // local ref
  if(fftLevel < grid().maxUniformRefinementLevel()) {
    for(MInt cellId = 0; cellId < a_noCells(); cellId++) {
      if(a_isHalo(cellId)) continue;
      if(a_isBndryGhostCell(cellId)) continue;
      if(a_level(cellId) <= fftLevel) continue;
      if(c_noChildren(cellId) > 0) continue;
      MInt parentId = c_parentId(cellId);
      while(a_level(parentId) != fftLevel) {
        parentId = c_parentId(parentId);
      }
      if(a_level(parentId) == fftLevel) {
        for(MInt varId = 0; varId < noVars; varId++) {
          a_variable(cellId, varId) = a_variable(parentId, varId);
          std::cout << "a_variable(" << cellId << "," << varId << "): " << a_variable(parentId, varId) << std::endl;
        }
      }
    }
  }
  exchangeData(&a_variable(0, 0), noVars);

  for(MInt cellId = 0; cellId < a_noCells(); cellId++) {
    for(MInt varId = 0; varId < noVars; varId++) {
      a_oldVariable(cellId, varId) = a_variable(cellId, varId);
    }
  }

  // initialize distr. functions
  initEqDistFunctions();

  m_log << "**************************" << endl;
  m_log << "Initial Condition summary" << endl;
  m_log << "**************************" << endl;
  m_log << "Re = " << m_Re << endl;
  m_log << "Ma = " << m_Ma << endl;
  m_log << "UInfinity = " << UInfinity << endl;
  m_log << "rhoInfinity = " << rhoInfinity << endl;

  // check initial condition
  for(MInt cellId = 0; cellId < noInternalCells(); cellId++) {
    for(MInt v = 0; v < PV->noVariables; v++) {
      if(std::isnan(a_variable(cellId, v))) {
        cerr << "Variable " << v << " cellId : " << cellId << endl;
        mTerm(1, AT_, "Invalid initialcondition formulation!");
      }
    }
  }

  // compute initial values for spectrum
  computeFFTStatistics();
}
template <>
void LbSolverDxQy<2, 9, maia::lb::LbSysEqnIncompressible<2, 9>>::initLatticeBgkFftIsotropicTurbulence() {
  mTerm(1, AT_, "Init method only available for D3Qy, yet!");
}

#ifdef WAR_NVHPC_PSTL
template <MInt nDim, MInt nDist, class SysEqn>
void LbSolverDxQy<nDim, nDist, SysEqn>::initArraysForPSTL() {
  TRACE();

  m_faculty.fill(-1);
  m_nFld.fill(-1);
  m_pFld.fill(-1);
  m_mFld1.fill(-1);
  m_mFld2.fill(-1);
  m_oppositeDist.fill(-1);

  for(MInt i = 0; i < 50; i++) {
    m_faculty[i] = faculty[i];
  }

  MInt fldlen = Ld::dxQyFld();
  for(MInt i = 0; i < fldlen; i++) {
    for(MInt d = 0; d < nDim; d++) {
      m_nFld[d * fldlen + i] = Ld::nFld(d, i);
      m_pFld[d * fldlen + i] = Ld::pFld(d, i);
    }
  }

  for(MInt i = 0; i < POWX(3, nDim); i++) {
    for(MInt j = 0; j < POWX(3, nDim); j++) {
      m_idFld[i][j] = Ld::idFld(i, j);
    }
  }

  for(MInt i = 0; i < 24; i++) {
    m_mFld1[i] = Ld::mFld1(i);
    m_mFld2[i] = Ld::mFld2(i);
  }

  for(MInt i = 0; i < nDist; i++) {
    m_oppositeDist[i] = Ld::oppositeDist(i);
  }

  for(MInt i = 0; i < 4; i++) {
    m_tp[i] = Ld::tp(i);
  }

  for(MInt i = 0; i < 3; i++) {
    m_distFld[i] = Ld::distFld(i);
  }

  for(MInt i = 0; i < nDist; i++) {
    for(MInt j = 0; j < nDim; j++) {
      m_ppdfDir[i * nDim + j] = Ld::ppdfDir(i, j);
    }
  }

  for(MInt i = 0; i < nDist; i++) {
    m_distType[i] = Ld::distType(i);
  }
}
#endif
