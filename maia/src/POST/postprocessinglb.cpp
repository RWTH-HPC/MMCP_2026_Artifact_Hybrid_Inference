// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "postprocessinglb.h"
//#include "postprocessing.h"
#include "postdata.h"
#include "postprocessing.cpp"

using namespace std;

template <MInt nDim>
PostProcessingLb<nDim>::PostProcessingLb(MInt postprocessingId_, PostData<nDim>* data, LbSolver<nDim>* ppSolver_)
  : PostProcessingInterface(postprocessingId_), PostProcessing<nDim, PostProcessingLb<nDim>>(postprocessingId_, data) {
  m_ppSolver = ppSolver_;
}

template <MInt nDim>
PostProcessingLb<nDim>::~PostProcessingLb() {}

template <MInt nDim>
void PostProcessingLb<nDim>::initPointSamplingData() {
  TRACE();

  m_pointData.reset(new PointData<nDim, SolverType>{*m_ppSolver});
  m_pointData->setInputOutputProperties();
  m_pointData->init();
}

template <MInt nDim>
void PostProcessingLb<nDim>::savePointSamplingData() {
  m_pointData->save(m_finalTimeStep);
}

template <MInt nDim>
void PostProcessingLb<nDim>::initSurfaceSamplingData() {
  TRACE();
  m_surfaceData.reset(new SurfaceData<nDim, SolverType>{*m_ppSolver});
  m_surfaceData->setInputOutputProperties();
  m_surfaceData->init();
}

template <MInt nDim>
void PostProcessingLb<nDim>::saveSurfaceSamplingData() {
  m_surfaceData->save(m_finalTimeStep);
}

template <MInt nDim>
void PostProcessingLb<nDim>::initVolumeSamplingData() {
  TRACE();
  m_volumeData.reset(new VolumeData<nDim, SolverType>{*m_ppSolver});
  m_volumeData->setInputOutputProperties();
  m_volumeData->init();
}

template <MInt nDim>
void PostProcessingLb<nDim>::saveVolumeSamplingData() {
  m_volumeData->save(m_finalTimeStep);
}

/** \brief init function for Isotropic Turbulence Statistics
 * \author Johannes Grafen
 * \date Jul 2022
 **/
template <MInt nDim>
void PostProcessingLb<nDim>::initIsoTurbulenceStatistics() {
  TRACE();

  if(solver().domainId() == 0) {
    std::ofstream ofl;
    ofl.open("IsoTurbulenceStatistics.log", ios::app);
    ofl << std::setw(5) << "t"
        << "\t" << std::setw(8) << "ReLambda"
        << "\t" << std::setw(8) << "lambda"
        << "\t" << std::setw(8) << "Eta"
        << "\t" << std::setw(8) << "tau_eta"
        << "\t" << std::setw(8) << "skewness" << endl;
    ofl.close();
  }
}

/**
 * \brief write data for isotropic Turbulence (single phase and particle-laden)
 *
 * \author Johannes Grafen
 * \date 28.02.2022
 *  writes the following parameters to std::cout:
 *  lambda
 *  ReLambda
 *  skewness
 *  eta (kolomogorovlength)
 *  tau_eta (kolomogorov time-scale)
 **/
template <MInt nDim>
void PostProcessingLb<nDim>::computeIsoTurbulenceStatistics() {
  TRACE();

  if((solver().m_fftInterval > 0 && globalTimeStep % solver().m_fftInterval == 0) || m_finalTimeStep) {
    const MInt fftLevel = solver().grid().maxUniformRefinementLevel();
    const MFloat DX = solver().c_cellLengthAtLevel(fftLevel);

    MFloat urms[6] = {F0, F0, F0, F0, F0, F0};
    MFloat dudx[3] = {F0, F0, F0};
    MFloat skew[4] = {F0, F0, F0, F0};
    MFloat lambda[3] = {F0, F0, F0};
    MFloat cnt = F0;

    MFloat lambdaAvg = F0;
    MFloat ReLambda = F0;
    MFloat skewness = F0;
    MFloat eta = F0;

    for(MInt cellId = 0; cellId < solver().a_noCells(); cellId++) {
      if(solver().a_isHalo(cellId)) continue;
      if(solver().a_isBndryGhostCell(cellId)) continue;
      if(solver().a_level(cellId) != fftLevel) continue;
      MFloat u = solver().a_variable(cellId, solver().PV->U);
      MFloat v = solver().a_variable(cellId, solver().PV->V);
      MFloat w = solver().a_variable(cellId, solver().PV->W);
      for(MInt i = 0; i < nDim; i++) {
        MInt n0 = (solver().a_hasNeighbor(cellId, 2 * i) > 0) ? solver().c_neighborId(cellId, 2 * i) : cellId;
        MInt n1 = (solver().a_hasNeighbor(cellId, 2 * i + 1) > 0) ? solver().c_neighborId(cellId, 2 * i + 1) : cellId;
        if(n0 == n1) continue;
        dudx[i] += POW2((solver().a_variable(n1, solver().PV->VV[i]) - solver().a_variable(n0, solver().PV->VV[i]))
                        / ((solver().a_coordinate(n1, i) - solver().a_coordinate(n0, i)) / DX));
        skew[i] += POW3((solver().a_variable(n1, solver().PV->VV[i]) - solver().a_variable(n0, solver().PV->VV[i]))
                        / ((solver().a_coordinate(n1, i) - solver().a_coordinate(n0, i)) / DX));
      }
      urms[0] += POW2(u);
      urms[1] += POW2(v);
      urms[2] += POW2(w);
      urms[3] += u * v;
      urms[4] += u * w;
      urms[5] += v * w;
      cnt++;
    }

    MPI_Allreduce(MPI_IN_PLACE, &cnt, 1, MPI_DOUBLE, MPI_SUM, solver().mpiComm(), AT_, "MPI_IN_PLACE", "cnt");
    MPI_Allreduce(MPI_IN_PLACE, &urms[0], 6, MPI_DOUBLE, MPI_SUM, solver().mpiComm(), AT_, "MPI_IN_PLACE", "urms[0]");
    MPI_Allreduce(MPI_IN_PLACE, &dudx[0], 3, MPI_DOUBLE, MPI_SUM, solver().mpiComm(), AT_, "MPI_IN_PLACE", "dudx[0]");
    MPI_Allreduce(MPI_IN_PLACE, &skew[0], 3, MPI_DOUBLE, MPI_SUM, solver().mpiComm(), AT_, "MPI_IN_PLACE", "skew[0]");
    skew[3] = ((skew[0] + skew[1] + skew[2]) / (F3 * cnt)) / pow((dudx[0] + dudx[1] + dudx[2]) / (F3 * cnt), 1.5);
    for(MInt i = 0; i < 3; i++)
      skew[i] = (skew[i] / cnt) / pow(dudx[i] / cnt, 1.5);
    for(MInt i = 0; i < 6; i++)
      urms[i] = sqrt(fabs(urms[i]) / cnt);
    for(MInt i = 0; i < 3; i++)
      dudx[i] /= cnt;
    for(MInt i = 0; i < 3; i++)
      lambda[i] = urms[i] / sqrt(dudx[i]); // Eq. 6.56 Turbulent Flows Pope two-point correlation

    lambdaAvg = F1B3 * (lambda[0] + lambda[1] + lambda[2]);
    ReLambda = F1B3 * (urms[0] + urms[1] + urms[2]) * lambdaAvg / solver().m_nu;
    skewness = skew[3];

    MFloat eps[3];
    for(MInt i = 0; i < 3; i++) {
      eps[i] = 15.0 * solver().m_nu * POW2(urms[i] / lambda[i]); // Eq. 6.58 Turbulent Flows Pope two-point correlation
    }
    eta = pow(solver().m_nu, 0.75)
          / pow(F1B3 * (eps[0] + eps[1] + eps[2]), 0.25); // Eq. 6.1 Turbulent Flows Pope two-point correlation
    m_tau_eta = pow(solver().m_nu / (F1B3 * (eps[0] + eps[1] + eps[2])), 0.5);

    if(solver().domainId() == 0) {
      std::ofstream ofl;
      ofl.open("IsoTurbulenceStatistics.log", ios::app);
      ofl << std::setw(5) << globalTimeStep << "\t" << std::setw(8) << ReLambda << "\t" << std::setw(8) << lambdaAvg
          << "\t" << std::setw(8) << eta << "\t" << std::setw(8) << m_tau_eta << "\t" << std::setw(8) << skewness
          << endl;
      ofl.close();
    }
  }
}

template class PostProcessingLb<2>;
template class PostProcessingLb<3>;
