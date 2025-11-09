// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef DGCOUPLINGCONDITIONACOUSTICPERTURB_H_
#define DGCOUPLINGCONDITIONACOUSTICPERTURB_H_

#include <algorithm>
#include <iostream>
#include <limits>
#include <sstream>
#include <type_traits>
#include "COMM/mpioverride.h"
#include "DG/dgcartesianelementcollector.h"
#include "DG/dgcartesiangalerkinprojection.h"
#include "DG/dgcartesiansyseqnacousticperturb.h"
#include "DG/dgcartesiantimers.h"
#include "IO/context.h"
#include "IO/parallelio.h"
#include "UTIL/timer.h"
#include "coupling.h"
#include "couplingdgape.h"
#include "typetraits.h"

class Solver;
template <MInt nDim, class SysEqn>
class DgCartesianSolver;

template <MInt nDim, class SysEqn>
class FvCartesianSolverXD;
/// \brief Coupling condition for direct-hybrid LES-CAA computations.
///
/// \author Michael Schlottke-Lakemper (mic) <mic@aia.rwth-aachen.de>
/// \date 2016-07-08
///
/// TODO labels:COUPLER,toremove Re-enable offline coupling or remove it entirely
template <MInt nDim, class FvSysEqn>
class DgCcAcousticPerturb final : public CouplingDgApe<nDim, CouplingFv<nDim, FvSysEqn>> {
  // Typedefs
 public:
  using Base = Coupling;

  using BaseFv = CouplingFv<nDim, FvSysEqn>;
  using FvCartesianSolverType = typename BaseFv::solverType;
  using BaseFv::fvSolver;

  using BaseDg = CouplingDgApe<nDim, BaseFv>;
  using BaseDg::dgSolver;
  using BaseDg::m_calcProjectionError;
  using BaseDg::m_calcSourceDonorCells;
  using BaseDg::m_checkConservation;
  using BaseDg::m_fixedTimeStep;
  using BaseDg::m_hasDgCartesianSolver;
  using BaseDg::m_hasDonorCartesianSolver;
  using BaseDg::m_isFirstSourceCalculation;
  using BaseDg::m_isRestart;
  using BaseDg::m_maxConservationError;
  using BaseDg::m_maxL2Error;
  using BaseDg::m_meanVars;
  using BaseDg::m_meanVarsIndex;
  using BaseDg::m_noActiveDonorCells;
  using BaseDg::m_noDonorCells;
  using BaseDg::m_previousTime;
  using BaseDg::m_timers;
  using BaseDg::noMeanVars;
  using BaseDg::noVelocities;
  using BaseDg::noVorticities;
  using BaseDg::outputDir;
  using BaseDg::sysEqn;
  using DgCartesianSolverType = typename BaseDg::DgCartesianSolverType;
  using ProjectionType = typename BaseDg::ProjectionType;
  using CV = typename BaseDg::CV;
  using MV = typename BaseDg::MV;
  using ST = typename BaseDg::ST;
  using SysEqn = typename BaseDg::SysEqn;
  using Timers = typename BaseDg::Timers;

  // TODO labels:COUPLER,DG use couplingId instead of solverId to request coupling related properties?
  using BaseDg::solverId;

  // Methods
  DgCcAcousticPerturb(const MInt couplingId, FvCartesianSolverType* fvSolver_, DgCartesianSolverType* const dgSolver_)
    : Base(couplingId), BaseDg(couplingId, dgSolver_, fvSolver_) {}
  virtual ~DgCcAcousticPerturb();

  // Main coupling functions
  void init() override {
    sanityCheck();
    BaseDg::init();
  };

  void finalizeCouplerInit() override { calcTimeStep(); };

  /// Load balancing
  void balancePre() override {
    m_loadBalancingReinitStage = 0; // Set reinitialization stage
    this->initCoupler();            // Perform main (re-)initialization of coupler
  };
  void balancePost() override {
    // Nothing to be done here
    m_loadBalancingReinitStage = -1; // Reset reinitialization stage
  };

  /// Methods to inquire coupler data during balancing
  MInt noCellDataDlb() const override {
    return 1; // One (float) data field to communicate
  };
  MInt cellDataTypeDlb(const MInt dataId) const override {
    TERMM_IF_NOT_COND(dataId == 0, "invalid data id");
    return MFLOAT;
  };
  MInt cellDataSizeDlb(const MInt NotUsed(dataId), const MInt NotUsed(cellId)) override;
  // Note: not-templated, only MFloat version required here
  void getCellDataDlb(const MInt dataId, const MInt oldNoCells, const MInt* const bufferIdToCellId,
                      MFloat* const data) override;
  void setCellDataDlb(const MInt NotUsed(dataId), const MFloat* const NotUsed(data)) override;

  /// Return true if time step calculation should be overruled by coupling
  MBool forceTimeStep() const { return true; }

 private:
  void sanityCheck();

  FvCartesianSolverType& donorSolver(const MInt xSolverId = 0) const override { return fvSolver(xSolverId); };
  void getDonorVelocityAndVorticity(const std::vector<MInt>& donorCellIds, MFloatScratchSpace& p_velocity,
                                    MFloatScratchSpace& p_vorticity) override;

  MFloat calcTimeStep();
  void calcSourceLambLinearized(const MFloat* const velocity, const MFloat* const vorticity,
                                MFloat* sourceTerms) override;
  void calcSourceLambNonlinear(const MFloat* const velocity, const MFloat* const vorticity,
                               MFloat* const sourceTerms) override;
  void calcSourceQmII(const MFloat* const velocity, MFloat* const sourceTerms) override;
  void calcSourceQmIII(const MFloat* const velocity, MFloat* sourceTerms) override;
  void calcSourceQe(const MFloat* const velocity, const MFloat time, MFloat* const sourceTerms) override;
  void calcSourceQc(const MFloat* const velocity, MFloat* const sourceTerms, const MFloat time,
                    const MInt timeStep) override;

  // Properties: file and variable names
  // File name pattern for instantaneous data files
  // MString m_instantaneousDataFileNamePattern{};
  // Field width for time step in file name pattern for instantaneous data files
  // MInt m_instantaneousDataFileNamePatternWidth;
  // Names of instantaneous velocity variables in solution files
  // labels:COUPLER Note: the default initializer '{}' was omitted due to a bug in GCC 4.9.x
  // std::array<MString, noVelocities()> m_instantaneousVelocitiesVarNames;
  // Names of instantaneous vorticity variables in solution files
  // labels:COUPLER Note: the default initializer '{}' was omitted due to a bug in GCC 4.9.x
  // std::array<MString, noVorticities()> m_instantaneousVorticitiesVarNames;


  MInt m_loadBalancingReinitStage = -1; ///< DLB reinitialization stage

  std::vector<MFloat> m_oldPressure{}; ///< Store previous pressure term for dp/dt calculation
};

/// \brief Stop the object lifetime timer and calculate global maximum errors.
///
/// \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
/// \date 2015-10-01
template <MInt nDim, class FvSysEqn>
DgCcAcousticPerturb<nDim, FvSysEqn>::~DgCcAcousticPerturb() {
  TRACE();

  MFloat maxGlobalConservationError = -1.0;
  MFloat maxGlobalL2Error = -1.0;

  // Compute global maximum conservation error if activated
  if(m_checkConservation && m_hasDgCartesianSolver) {
    MPI_Allreduce(&m_maxConservationError, &maxGlobalConservationError, 1, maia::type_traits<MFloat>::mpiType(),
                  MPI_MAX, dgSolver().mpiComm(), AT_, "m_maxConservationError", "maxGlobalConservationError");

    m_log << "Global maximum Galerkin projection conservation error: " << std::scientific << maxGlobalConservationError
          << std::endl;
    m_log << "Local maximum Galerkin projection conservation error: " << std::scientific << m_maxConservationError
          << std::endl;
  }

  // Compute global maximum L2 projection error if activated
  if(m_calcProjectionError && m_hasDgCartesianSolver) {
    MPI_Allreduce(&m_maxL2Error, &maxGlobalL2Error, 1, maia::type_traits<MFloat>::mpiType(), MPI_MAX,
                  dgSolver().mpiComm(), AT_, "m_maxL2Error", "maxGlobalL2Error");

    m_log << "Global maximum Galerking projection L2 error: " << std::scientific << maxGlobalL2Error << std::endl;
    m_log << "Local maximum Galerking projection L2 error: " << std::scientific << m_maxL2Error << std::endl;
  }

  // Write the maximum global errors to file
  if(m_hasDgCartesianSolver && dgSolver().domainId() == 0 && (m_checkConservation || m_calcProjectionError)) {
    std::ofstream errorFile(outputDir() + "errors.txt");

    if(m_checkConservation) {
      errorFile << "maxGlobalConservationError: " << std::scientific << maxGlobalConservationError << std::endl;
    }
    if(m_calcProjectionError) {
      errorFile << "maxGlobalProjectionError:   " << std::scientific << maxGlobalL2Error << std::endl;
    }

    errorFile.close();
  }

  RECORD_TIMER_STOP(m_timers[Timers::Class]);
}


/** \brief 	Perform sanity checks
 *  \author Miro Gondrum (refactored)
 *  \date 	28.04.2021
 */
template <MInt nDim, class FvSysEqn>
void DgCcAcousticPerturb<nDim, FvSysEqn>::sanityCheck() {
  TRACE();

  if(!fvSolver().calcSlopesAfterStep()) {
    TERMM(1,
          "Enable calcSlopesAfterStep for FV-solver for coupled simulation (else the slopes are not "
          "initialized prior to the first source term calculation and will contain garbage)");
  }
}

/// \brief Calculate the nonlinear Lamb vector coupling source terms on all
///        donor cells where it is needed.
///
/// \author Ansgar Niemoeller (ansgar) <a.niemoeller@aia.rwth-aachen.de>
/// \date 2016-08-18
///
/// \param[in] velocity Velocities of all donor cells.
/// \param[in] vorticity Vorticities of all donor cells.
/// \param[out] sourceTerms Storage for computed source terms.
template <MInt nDim, class FvSysEqn>
void DgCcAcousticPerturb<nDim, FvSysEqn>::calcSourceLambNonlinear(const MFloat* const velocity,
                                                                  const MFloat* const vorticity,
                                                                  MFloat* const sourceTerms) {
  TRACE();
  RECORD_TIMER_START(m_timers[Timers::CalcSourceLamb]);

  // Index of first lamb vector component in m_meanVars
  const MInt lambIndex = m_meanVarsIndex[MV::LAMB0[0]];

  std::array<MFloat, MAX_SPACE_DIMENSIONS> lambVector{};

  // Calculate source terms on all 'active' donor leaf cells, i.e. that are
  // mapped to elements with nonzero filter values
  for(MInt donorIndex = 0; donorIndex < m_noActiveDonorCells; donorIndex++) {
    const MInt donorId = m_calcSourceDonorCells[donorIndex];
    // Create convenience pointers
    const MFloat* meanVars = &m_meanVars[donorIndex * noMeanVars()];
    const MFloat* velo = &velocity[donorId * noVelocities()];
    const MFloat* vort = &vorticity[donorId * noVorticities()];

    // Compute instantaneous Lamb vector
    if constexpr(nDim == 2) {
      lambVector[0] = -vort[0] * velo[1];
      lambVector[1] = vort[0] * velo[0];
    } else {
      lambVector[0] = vort[1] * velo[2] - vort[2] * velo[1];
      lambVector[1] = vort[2] * velo[0] - vort[0] * velo[2];
      lambVector[2] = vort[0] * velo[1] - vort[1] * velo[0];
    }

    // The complete source term is:
    // S = -L' = -(w x u - bar(w x u))
    // Since the Lamb vector is defined with vorticity as 'nabla x u' but MAIA
    // calculates it as '0.5 nabla x u', we multiply the resulting Lamb vector
    // by 2.0
    for(MInt dim = 0; dim < nDim; dim++) {
      // Note: source terms are accumulated, thus "+=" and not "="
      sourceTerms[donorId * SysEqn::noVars() + dim] += -2.0 * (lambVector[dim] - meanVars[lambIndex + dim]);
    }
  }

  RECORD_TIMER_STOP(m_timers[Timers::CalcSourceLamb]);
}


/// \brief Calculate the q_mII source terms on all donor cells where it is needed.
///
/// \author Bjoern Peeters (bjoern) <b.peeters@aia.rwth-aachen.de>
/// \date 2017-05-01
///
/// \param[in] velocity Velocities of all donor cells.
/// \param[in] time Current simulation time.
/// \param[out] sourceTerms Storage for computed source terms.
template <MInt nDim, class FvSysEqn>
void DgCcAcousticPerturb<nDim, FvSysEqn>::calcSourceQmII(const MFloat* const NotUsed(velocity),
                                                         MFloat* const sourceTerms) {
  TRACE();
  RECORD_TIMER_START(m_timers[Timers::CalcSourceQmII]);

  const MInt rhoIndex = m_meanVarsIndex[MV::RHO0];
  const MInt pIndex = m_meanVarsIndex[MV::P0];
  const MInt drhoIndex = m_meanVarsIndex[MV::DRHO[0]];
  const MInt dpIndex = m_meanVarsIndex[MV::DP[0]];
  const MInt dprhoIndex = m_meanVarsIndex[MV::GRADPRHO[0]];

  // Calculate source terms on all donor leaf cells that are mapped to elements
  // with nonzero filter values
  for(MInt donorIndex = 0; donorIndex < m_noActiveDonorCells; donorIndex++) {
    const MInt donorId = m_calcSourceDonorCells[donorIndex];
    // Create convenience pointers
    const MFloat* meanVars = &m_meanVars[donorIndex * noMeanVars()];
    const MFloat curRho = fvSolver().a_variable(donorId, fvSolver().CV->RHO);
    const MFloat curP = fvSolver().a_pvariable(donorId, fvSolver().PV->P);
    const MFloat rhoM = meanVars[rhoIndex];
    const MFloat pM = meanVars[pIndex];

    // bar(grad(rho))
    std::array<MFloat, nDim> drhoM;
    for(MInt dim = 0; dim < nDim; dim++) {
      drhoM[dim] = meanVars[drhoIndex + dim];
    }

    // bar(grad(p))
    std::array<MFloat, nDim> dpM;
    for(MInt dim = 0; dim < nDim; dim++) {
      dpM[dim] = meanVars[dpIndex + dim];
    }

    // bar(grad(p)/rho)
    std::array<MFloat, nDim> dprhoM;
    for(MInt dim = 0; dim < nDim; dim++) {
      dprhoM[dim] = meanVars[dprhoIndex + dim];
    }

    std::array<MFloat, nDim> qmIISource;
    for(MInt dim = 0; dim < nDim; dim++) {
      // Calculate source term, starting with term 1
      qmIISource[dim] = ((fvSolver().a_slope(donorId, fvSolver().PV->P, dim)) - dpM[dim]) / rhoM
                        // term 2
                        - ((curP - pM) / (rhoM * rhoM)) * drhoM[dim]
                        // term 3
                        - (fvSolver().a_slope(donorId, fvSolver().PV->P, dim)) / curRho
                        // term 4
                        + dprhoM[dim];
    }

    // Store source term at output memory location
    for(MInt dim = 0; dim < nDim; dim++) {
      // Note: source terms are accumulated, thus "+=" and not "="
      sourceTerms[donorId * SysEqn::noVars() + CV::UU[0] + dim] += qmIISource[dim];
    }
  }

  RECORD_TIMER_STOP(m_timers[Timers::CalcSourceQmII]);
}


/// \brief Calculate the q_mIII source terms on all donor cells where it is needed.
///
/// TODO labels:COUPLER,totest requires testing!
template <MInt nDim, class FvSysEqn>
void DgCcAcousticPerturb<nDim, FvSysEqn>::calcSourceQmIII(const MFloat* const velocity, MFloat* sourceTerms) {
  TRACE();
  RECORD_TIMER_START(m_timers[Timers::CalcSourceQmIII]);
  TERMM(1, "requires testing");

  const MInt umIndex = m_meanVarsIndex[MV::UU0[0]];

  // Calculate source terms on all donor leaf cells that are mapped to elements
  // with nonzero filter values
  for(MInt donorIndex = 0; donorIndex < m_noActiveDonorCells; donorIndex++) {
    const MInt donorId = m_calcSourceDonorCells[donorIndex];
    // Create convenience pointers
    const MFloat* meanVars = &m_meanVars[donorIndex * noMeanVars()];

    std::array<MFloat, nDim> qmIIISource;
    std::fill_n(&qmIIISource[0], nDim, 0.0);

    // Calculate source term
    for(MInt dim = 0; dim < nDim; dim++) {
      const MFloat velDiff = velocity[donorId * nDim + dim] - meanVars[umIndex + dim];

      for(MInt i = 0; i < nDim; i++) {
        const MFloat velGrad = fvSolver().a_slope(donorId, fvSolver().PV->VV[dim], i);
        const MFloat velGradM = meanVars[m_meanVarsIndex[MV::GRADU[dim * nDim + i]]];

        qmIIISource[i] += velDiff * (velGrad - velGradM);
      }

      // Substract mean term
      qmIIISource[dim] -= meanVars[m_meanVarsIndex[MV::UGRADU[dim]]];

      for(MInt i = 0; i < nDim; i++) {
        const MInt gradUPos = i * nDim + dim;
        qmIIISource[dim] += meanVars[m_meanVarsIndex[MV::UU0[i]]] * meanVars[m_meanVarsIndex[MV::GRADU[gradUPos]]];
      }
    }

    // Store source term at output memory location
    for(MInt dim = 0; dim < nDim; dim++) {
      // Note: source terms are accumulated, thus "-=" and not "="
      sourceTerms[donorId * SysEqn::noVars() + CV::UU[0] + dim] -= qmIIISource[dim];
    }
  }

  RECORD_TIMER_STOP(m_timers[Timers::CalcSourceQmIII]);
}


/// \brief Calculate the linearized Lamb vector coupling source terms on all
///        donor cells where it is needed.
///
/// \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
/// \date 2015-10-01
///
/// \param[in] velocity Velocities of all donor cells.
/// \param[in] vorticity Vorticities of all donor cells.
/// \param[out] sourceTerms Storage for computed source terms.
template <MInt nDim, class FvSysEqn>
void DgCcAcousticPerturb<nDim, FvSysEqn>::calcSourceLambLinearized(const MFloat* const velocity,
                                                                   const MFloat* const vorticity,
                                                                   MFloat* sourceTerms) {
  TRACE();
  RECORD_TIMER_START(m_timers[Timers::CalcSourceLamb]);

  // Indices of first velocity/vorticity component in m_meanVars
  const MInt veloIndex = m_meanVarsIndex[MV::UU0[0]];
  const MInt vortIndex = m_meanVarsIndex[MV::VORT0[0]];

  std::array<MFloat, MAX_SPACE_DIMENSIONS> veloP{};
  std::array<MFloat, MAX_SPACE_DIMENSIONS> vortP{};

  // Calculate source terms on all 'active' donor leaf cells, i.e. that are
  // mapped to elements with nonzero filter values
  for(MInt donorIndex = 0; donorIndex < m_noActiveDonorCells; donorIndex++) {
    const MInt donorId = m_calcSourceDonorCells[donorIndex];
    // Create convenience pointers
    const MFloat* meanVars = &m_meanVars[donorIndex * noMeanVars()];
    const MFloat* velo = &velocity[donorId * noVelocities()];
    const MFloat* vort = &vorticity[donorId * noVorticities()];

    // Compute velocity fluctuations
    for(MInt dim = 0; dim < nDim; dim++) {
      veloP[dim] = velo[dim] - meanVars[veloIndex + dim];
    }
    // Compute vorticity fluctuations
    for(MInt i = 0; i < noVorticities(); i++) {
      vortP[i] = vort[i] - meanVars[vortIndex + i];
    }

    // The complete source term is:
    // S = -L' = -(w' X u0 + w0 X u')
    // Since the Lamb vector is defined with vorticity as 'nabla x u' but MAIA
    // calculates it as '0.5 nabla x u', we multiply the resulting Lamb vector
    // by 2.0
    const MInt sourcesOffset = donorId * SysEqn::noVars() + CV::UU[0];
    IF_CONSTEXPR(nDim == 2) { // 2D
      // x-component: -(-wz' * u0y - w0z * uy')
      sourceTerms[sourcesOffset] += -2.0 * (-vortP[0] * meanVars[veloIndex + 1] - meanVars[vortIndex] * veloP[1]);

      // y-component: -(wz' * u0x + w0z * ux')
      sourceTerms[sourcesOffset + 1] += -2.0 * (vortP[0] * meanVars[veloIndex] + meanVars[vortIndex] * veloP[0]);
    }
    else { // 3D
      // x-component: -(wy' * u0z - wz' * u0y + w0y * uz' - w0z * uy')
      sourceTerms[sourcesOffset] += -2.0
                                    * (vortP[1] * meanVars[veloIndex + 2] - vortP[2] * meanVars[veloIndex + 1]
                                       + meanVars[vortIndex + 1] * veloP[2] - meanVars[vortIndex + 2] * veloP[1]);

      // y-component: -(wz' * u0x - wx' * u0z + w0z * ux' - w0x * uz')
      sourceTerms[sourcesOffset + 1] += -2.0
                                        * (vortP[2] * meanVars[veloIndex] - vortP[0] * meanVars[veloIndex + 2]
                                           + meanVars[vortIndex + 2] * veloP[0] - meanVars[vortIndex] * veloP[2]);

      // z-component: -(wx' * u0y - wy' * u0x + w0x * uy' - w0y * ux')
      sourceTerms[sourcesOffset + 2] += -2.0
                                        * (vortP[0] * meanVars[veloIndex + 1] - vortP[1] * meanVars[veloIndex]
                                           + meanVars[vortIndex] * veloP[1] - meanVars[vortIndex + 1] * veloP[0]);
    }
  }

  RECORD_TIMER_STOP(m_timers[Timers::CalcSourceLamb]);
}


/// \brief Calculate the q_e source terms on all donor cells where it is needed.
///
/// \author Bjoern Peeters (bjoern) <b.peeters@aia.rwth-aachen.de>
/// \date 2017-05-01
///
/// \param[in] velocity Velocities of all donor cells.
/// \param[in] time Current simulation time.
/// \param[out] sourceTerms Storage for computed source terms.
template <MInt nDim, class FvSysEqn>
void DgCcAcousticPerturb<nDim, FvSysEqn>::calcSourceQe(const MFloat* const velocity,
                                                       const MFloat time,
                                                       MFloat* const sourceTerms) {
  TRACE();
  RECORD_TIMER_START(m_timers[Timers::CalcSourceQe]);
  TERMM(1, "untested in new coupling condition");

  // If this is the first invocation, we need to perform some initializations.
  if(m_isFirstSourceCalculation) {
    // Create persistent storage for old pressure values
    m_oldPressure.resize(m_noActiveDonorCells);

    if(m_isRestart) {
      // If this is a restart, initialize old pressure values from oldVariables, which should have
      // been loaded by setting the property "restartOldVariables" to true.
      // @ansgar TODO labels:COUPLER,FV support for balance, redistribute oldPressure/oldVariables in FV solver
      MFloatScratchSpace buffer(m_noDonorCells, AT_, "buffer");
      fvSolver().oldPressure(&buffer[0]);
      for(MInt donorIndex = 0; donorIndex < m_noActiveDonorCells; donorIndex++) {
        const MInt donorId = m_calcSourceDonorCells[donorIndex];
        m_oldPressure[donorIndex] = buffer[donorId];
      }
    } else {
      // If this is *not* a restart, assume that this is the first time step. In this case, the old
      // pressure values are initialized with the current values to have a gradient of zero
      for(MInt donorIndex = 0; donorIndex < m_noActiveDonorCells; donorIndex++) {
        const MInt donorId = m_calcSourceDonorCells[donorIndex];
        m_oldPressure[donorIndex] = fvSolver().a_pvariable(donorId, fvSolver().PV->P);
      }
    }

    // Calculate slopes for the calculation of the velocity divergence
    fvSolver().LSReconstructCellCenter();
  }

  // Calculate dt (set to infinity if this is the first time step of a non-restarted simulation)
  const MFloat dt =
      (m_isFirstSourceCalculation && !m_isRestart) ? std::numeric_limits<MFloat>::infinity() : time - m_previousTime;

  // Calculate source terms on all donor leaf cells that are mapped to elements with nonzero filter
  // values
  const MInt rhoIndex = m_meanVarsIndex[MV::RHO0];
  const MInt pIndex = m_meanVarsIndex[MV::P0];
  const MInt cIndex = m_meanVarsIndex[MV::C0];
  const MInt veloIndex = m_meanVarsIndex[MV::UU0[0]];
  const MInt dcIndex = m_meanVarsIndex[MV::DC0[0]];
  std::array<MInt, nDim> divIndex;
  for(MInt i = 0; i < nDim; i++) {
    divIndex[i] = m_meanVarsIndex[MV::DU[i]];
  }
  const MInt drhoIndex = m_meanVarsIndex[MV::DRHO[0]];
  const MInt dpIndex = m_meanVarsIndex[MV::DP[0]];

  for(MInt donorIndex = 0; donorIndex < m_noActiveDonorCells; donorIndex++) {
    const MInt donorId = m_calcSourceDonorCells[donorIndex];

    // Calculate dp/dt and store current pressure for the next time step
    const MFloat curP = fvSolver().a_pvariable(donorId, fvSolver().PV->P);
    const MFloat dpdt = (curP - m_oldPressure[donorIndex]) / dt;
    m_oldPressure[donorIndex] = curP;

    // Create convenience pointers
    const MFloat* meanVars = &m_meanVars[donorIndex * noMeanVars()];

    // bar(u)
    std::array<MFloat, nDim> veloM;
    for(MInt dim = 0; dim < nDim; dim++) {
      veloM[dim] = meanVars[veloIndex + dim];
    }

    // bar(grad(c))
    std::array<MFloat, nDim> dcM;
    for(MInt dim = 0; dim < nDim; dim++) {
      dcM[dim] = meanVars[dcIndex + dim];
    }

    // bar(grad(rho))
    std::array<MFloat, nDim> drhoM;
    for(MInt dim = 0; dim < nDim; dim++) {
      drhoM[dim] = meanVars[drhoIndex + dim];
    }

    // bar(grad(p))
    std::array<MFloat, nDim> dpM;
    for(MInt dim = 0; dim < nDim; dim++) {
      dpM[dim] = meanVars[dpIndex + dim];
    }

    // Mean of du/dx, dv/dy, dw/dz for the divergence u
    MFloat divUm = 0.0;
    for(MInt dim = 0; dim < nDim; dim++) {
      divUm += meanVars[divIndex[dim]];
    }

    MFloat divU = 0.0;
    for(MInt dim = 0; dim < nDim; dim++) {
      divU += fvSolver().a_slope(donorId, fvSolver().PV->VV[dim], dim);
    }

    // Compute drho/dt = -(rho* div(u)  + u * grad(rho))
    const MFloat* const velo = &velocity[donorId * noVelocities()];
    MFloat drhodttmp = 0.0;
    for(MInt dim = 0; dim < nDim; dim++) {
      drhodttmp += velo[dim] * fvSolver().a_slope(donorId, fvSolver().PV->RHO, dim);
    }
    const MFloat curRho = fvSolver().a_variable(donorId, fvSolver().CV->RHO);
    const MFloat drhodt = -curRho * divU - drhodttmp;

    // Compute q_e
    const MFloat rhoM = meanVars[rhoIndex];
    const MFloat pM = meanVars[pIndex];
    const MFloat cM = meanVars[cIndex];
    const MFloat cTerm = cM * cM;
    MFloat qeSource = 0.0;
    qeSource += drhodt - (1.0 / cTerm) * dpdt;
    qeSource += (curRho - rhoM) * divUm;
    qeSource -= ((curP - pM) / cTerm) * divUm;
    for(MInt dim = 0; dim < nDim; dim++) {
      // bar(u) * (grad(rho) - bar(grad(rho)))
      qeSource += veloM[dim] * ((fvSolver().a_slope(donorId, fvSolver().PV->RHO, dim)) - drhoM[dim]);
      qeSource -= (1.0 / cTerm) * veloM[dim] * ((fvSolver().a_slope(donorId, fvSolver().PV->P, dim)) - dpM[dim]);
      qeSource += 2.0 * ((curP - pM) / (cM * cM * cM)) * veloM[dim] * dcM[dim];
    }

    // Note: the multiplication with c^2 is technically not part of q_e
    qeSource *= -cTerm;

    // Store source term at output memory location
    // Note: source terms are accumulated, thus "+=" and not "="
    sourceTerms[donorId * SysEqn::noVars() + CV::P] += qeSource;
  }

  RECORD_TIMER_STOP(m_timers[Timers::CalcSourceQe]);
}


/// \brief Calculate the q_c source terms on all donor cells where it is needed.
///
/// \author Bjoern Peeters (bjoern) <b.peeters@aia.rwth-aachen.de>
/// \date 2017-05-01
///
/// \param[in] velocity Velocities of all donor cells.
/// \param[in] time Current simulation time.
/// \param[in] timeStep Current simulation time step.
/// \param[out] sourceTerms Storage for computed source terms.
template <MInt nDim, class FvSysEqn>
void DgCcAcousticPerturb<nDim, FvSysEqn>::calcSourceQc(const MFloat* const velocity,
                                                       MFloat* const sourceTerms,
                                                       const MFloat NotUsed(time),
                                                       const MInt NotUsed(timeStep)) {
  TRACE();
  RECORD_TIMER_START(m_timers[Timers::CalcSourceQc]);

  const MInt veloIndex = m_meanVarsIndex[MV::UU0[0]];
  const MInt drhoIndex = m_meanVarsIndex[MV::DRHO[0]];
  const MInt rhoDivUIndex = m_meanVarsIndex[MV::RHODIVU[0]];
  const MInt uGradRhoIndex = m_meanVarsIndex[MV::UGRADRHO[0]];
  std::array<MInt, nDim> divIndex;
  for(MInt i = 0; i < nDim; i++) {
    divIndex[i] = m_meanVarsIndex[MV::DU[i]];
  }
  const MInt cIndex = m_meanVarsIndex[MV::C0];
  const MInt rhoIndex = m_meanVarsIndex[MV::RHO0];

  // Calculate source terms on all donor leaf cells that are mapped to
  // elements
  // with nonzero filter values
  for(MInt donorIndex = 0; donorIndex < m_noActiveDonorCells; donorIndex++) {
    const MInt donorId = m_calcSourceDonorCells[donorIndex];

    // Create convenience pointer
    const MFloat* const meanVars = &m_meanVars[donorIndex * noMeanVars()];

    // bar(u)
    std::array<MFloat, nDim> veloM;
    for(MInt dim = 0; dim < nDim; dim++) {
      veloM[dim] = meanVars[veloIndex + dim];
    }

    // bar(grad(rho))
    std::array<MFloat, nDim> drhoM;
    for(MInt dim = 0; dim < nDim; dim++) {
      drhoM[dim] = meanVars[drhoIndex + dim];
    }

    // bar(rho * div(u))
    MFloat rhodivuM = 0.0;
    for(MInt dim = 0; dim < nDim; dim++) {
      rhodivuM += meanVars[rhoDivUIndex + dim];
    }

    // bar(u * grad(rho))
    MFloat ugradrhoM = 0.0;
    for(MInt dim = 0; dim < nDim; dim++) {
      ugradrhoM += meanVars[uGradRhoIndex + dim];
    }

    MFloat divUm = 0.0;
    for(MInt dim = 0; dim < nDim; dim++) {
      divUm += meanVars[divIndex[dim]];
    }

    MFloat divU = 0.0;
    divU += fvSolver().a_slope(donorId, fvSolver().PV->U, 0);
    divU += fvSolver().a_slope(donorId, fvSolver().PV->V, 1);
    IF_CONSTEXPR(nDim == 3) { divU += fvSolver().a_slope(donorId, fvSolver().PV->W, 2); }

    // Compute q_c
    const MFloat rhoM = meanVars[rhoIndex];
    const MFloat curRho = fvSolver().a_variable(donorId, fvSolver().CV->RHO);
    const MFloat* const velo = &velocity[donorId * noVelocities()];
    MFloat qcSource = 0.0;
    qcSource += (curRho - rhoM) * (divU - divUm);
    qcSource -= rhodivuM;
    qcSource -= ugradrhoM;
    qcSource += rhoM * divUm;
    for(MInt dim = 0; dim < nDim; dim++) {
      qcSource += (velo[dim] - veloM[dim]) * ((fvSolver().a_slope(donorId, fvSolver().PV->RHO, dim)) - drhoM[dim]);
      qcSource += veloM[dim] * drhoM[dim];
    }

    // Note: the multiplication with c^2 is technically not part of q_c
    qcSource *= -meanVars[cIndex] * meanVars[cIndex];

    // Store source term at output memory location
    // Note: source terms are accumulated, thus "+=" and not "="
    sourceTerms[donorId * SysEqn::noVars() + CV::P] += qcSource;
  }

  RECORD_TIMER_STOP(m_timers[Timers::CalcSourceQc]);
}

template <MInt nDim, class FvSysEqn>
void DgCcAcousticPerturb<nDim, FvSysEqn>::getDonorVelocityAndVorticity(const std::vector<MInt>& donorCellIds,
                                                                       MFloatScratchSpace& p_velocity,
                                                                       MFloatScratchSpace& p_vorticity) {
  // Load velocities of all donor leaf cells (if there are any on this domain)
  for(auto donorId : donorCellIds) {
    for(MInt varId = 0; varId < noVelocities(); varId++) {
      p_velocity(donorId, varId) = fvSolver().a_pvariable(donorId, fvSolver().PV->VV[varId]);
    }
  }

  // Only calculate vorticity if this domain has FV cells at all, otherwise
  // the FV solver pointer is null and a segmentation fault occurs
  // TODO labels:COUPLER,FV only compute vorticity on donor cells where it is needed?
  if(m_hasDonorCartesianSolver) {
    fvSolver().getVorticityT(&p_vorticity[0]);
  }
}

/// \brief Return the time step size dt as determined by the coupling condition.
///
/// \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
/// \date 2015-10-01
///
/// \return The next time step size (dt).
template <MInt nDim, class FvSysEqn>
MFloat DgCcAcousticPerturb<nDim, FvSysEqn>::calcTimeStep() {
  TRACE();

  if(m_fixedTimeStep < 0.0) {
    TERMM(1, "Fixed time step less than zero.");
  }

  fvSolver().forceTimeStep(m_fixedTimeStep);
  dgSolver().forceTimeStep(m_fixedTimeStep);

  // FIXME labels:COUPLER Return proper time step determined by both solvers
  return m_fixedTimeStep;
}

template <MInt nDim, class FvSysEqn>
MInt DgCcAcousticPerturb<nDim, FvSysEqn>::cellDataSizeDlb(const MInt dataId, const MInt gridCellId) {
  TRACE();
  TERMM_IF_NOT_COND(dataId == 0, "invalid data id");

  // Convert to solver cell id and check
  const MInt fvCellId = fvSolver().grid().tree().grid2solver(gridCellId);
  if(fvCellId < 0) {
    return 0;
  }

  MInt dataSize = 0;
  // TODO labels:COUPLER Check that ids in m_calcSourceDonorCells are sorted in get/setCellDataDlb
  const MBool foundDonorCell =
      std::binary_search(m_calcSourceDonorCells.begin(), m_calcSourceDonorCells.end(), fvCellId);
  if(foundDonorCell) {
    dataSize = noMeanVars();
  }
  return dataSize;
}


template <MInt nDim, class FvSysEqn>
void DgCcAcousticPerturb<nDim, FvSysEqn>::getCellDataDlb(const MInt dataId,
                                                         const MInt NotUsed(oldNoCells),
                                                         const MInt* const NotUsed(bufferIdToCellId),
                                                         MFloat* const data) {
  TRACE();
  TERMM_IF_NOT_COND(dataId == 0, "invalid data id");

  const MBool donorIdsSorted = std::is_sorted(m_calcSourceDonorCells.begin(), m_calcSourceDonorCells.end());
  TERMM_IF_NOT_COND(donorIdsSorted, "donor ids not sorted");

  if(m_noActiveDonorCells > 0) {
    std::copy_n(&m_meanVars[0], m_noActiveDonorCells * noMeanVars(), data);
  }
}


template <MInt nDim, class FvSysEqn>
void DgCcAcousticPerturb<nDim, FvSysEqn>::setCellDataDlb(const MInt dataId, const MFloat* const data) {
  TRACE();
  TERMM_IF_NOT_COND(dataId == 0, "invalid data id");

  // Return if not on the right reinitialization stage
  if(m_loadBalancingReinitStage != 0) {
    return;
  }

  const MBool donorIdsSorted = std::is_sorted(m_calcSourceDonorCells.begin(), m_calcSourceDonorCells.end());
  TERMM_IF_NOT_COND(donorIdsSorted, "donor ids not sorted");

  if(m_noActiveDonorCells > 0) {
    std::copy_n(data, m_noActiveDonorCells * noMeanVars(), &m_meanVars[0]);
  }
}

#endif // DGCOUPLINGCONDITIONACOUSTICPERTURB_H_
