// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "postprocessinglblpt.h"
#include "postprocessing.cpp"

using namespace std;

template <MInt nDim>
PostProcessingLbLPT<nDim>::PostProcessingLbLPT(MInt postprocessingId_,
                                               PostData<nDim>* data,
                                               LbSolver<nDim>* ppLbSolver_,
                                               LPT<nDim>* ppLptSolver_)
  : PostProcessingInterface(postprocessingId_),
    PostProcessingLPT<nDim>(postprocessingId_, data, ppLptSolver_),
    PostProcessingLb<nDim>(postprocessingId_, data, ppLbSolver_) {
  m_ppSolverLpt = ppLptSolver_;
  m_ppSolverLb = ppLbSolver_;
}

template <MInt nDim>
void PostProcessingLbLPT<nDim>::getConversionFactors() {
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
  lengthFactor = Context::getSolverProperty<MFloat>("lengthFactor", m_postprocessingId, AT_, &lengthFactor);

  const MFloat dx = lbSolver().c_cellLengthAtLevel(lbSolver().maxLevel());

  m_conversionLbLptLength = lengthFactor * dx;
  m_conversionLptLbLength = F1 / m_conversionLbLptLength;
}

/** \brief init function for particle-laden isotropic turbulence
 * \author Johannes Grafen
 * \date Jul 2022
 **/
template <MInt nDim>
void PostProcessingLbLPT<nDim>::initPLIsoTurbulenceStatistics() {
  TRACE();

  getConversionFactors();
  if(lbSolver().domainId() == 0) {
    std::ofstream PLlog;
    PLlog.open("PLIsoTurbStatistics.log", ios::app);
    PLlog << std::setw(5) << "t"
          << "\t" << std::setw(10) << "tau_p"
          << "\t" << std::setw(10) << "St_eta" << std::endl;
    PLlog.close();
  }
}

/**
 * \brief compute average quantites of particle-laden isotropic turbulence
 *    average particle response time tau_p
 *    average Stokes number in respect of the Kolmogorov length
 * \author Johannes Grafen
 * \date Jul 2022
 *
 **/
template <MInt nDim>
void PostProcessingLbLPT<nDim>::computePLIsoTurbulenceStatistics() {
  TRACE();

  if((lbSolver().a_FFTInterval() > 0 && globalTimeStep % lbSolver().a_FFTInterval() == 0) || m_finalTimeStep) {
    const MInt globalNoPart = lptSolver().globalNoParticles();
    MFloat tau_P = F0;
    MFloat St_eta = F0;

    for(const auto& part : lptSolver().m_partList) {
      tau_P += 1.0 / 18.0 * part.m_densityRatio * POW2(part.m_diameter * m_conversionLptLbLength) / lbSolver().a_Nu();
      St_eta += tau_P / m_tau_eta; // Stokes number
    }

    MPI_Allreduce(MPI_IN_PLACE, &tau_P, 1, MPI_DOUBLE, MPI_SUM, lbSolver().mpiComm(), AT_, "INPLACE", "tau_P");

    MPI_Allreduce(MPI_IN_PLACE, &St_eta, 1, MPI_DOUBLE, MPI_SUM, lbSolver().mpiComm(), AT_, "INPLACE", "St_eta");

    MFloat tau_P_Avg = tau_P / globalNoPart;
    MFloat St_eta_Avg = St_eta / globalNoPart;

    if(lbSolver().domainId() == 0) {
      std::ofstream PLlog;
      PLlog.open("PLIsoTurbStatistics.log", ios::app);
      PLlog << std::setw(5) << globalTimeStep << "\t" << std::setw(10) << tau_P_Avg << "\t" << std::setw(10)
            << St_eta_Avg << std::endl;
      PLlog.close();
    }
  }
}


template class PostProcessingLbLPT<3>;
