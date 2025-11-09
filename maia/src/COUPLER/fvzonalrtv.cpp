// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "fvzonalrtv.h"

#include <algorithm>
#include <iostream>
#include <stack>
#include <vector>
#include "FV/fvcartesiancellcollector.h"
#include "FV/fvransmodelconstants.h"
#include "MEMORY/alloc.h"
#include "UTIL/functions.h"
#include "UTIL/kdtree.h"
#include "globals.h"
#include "globalvariables.h"

using namespace std;

template <MInt nDim, class SysEqn>
FvZonalRTV<nDim, SysEqn>::FvZonalRTV(const MInt couplingId, RANS* R, LES* L)
  : Coupling(couplingId), FvZonal<nDim, SysEqn>(couplingId, R, L) /*, m_RANSSolver(R), m_LESSolver(L)*/ {
  m_log << "Create Zonal coupler for RANS Solver (" << RANSSolver().m_solverId << ") and LES Solver ("
        << LESSolver().m_solverId << ")" << endl;


  /*! \page propertyPage1
    \section zonal
    <code>MInt FvZonalRTV::m_rntStartTimeStep</code>\n
    default = <code>False</code>\n
    Triggers reconstruction of LES Average values form nut after nonZonalRestart
    Keywords: <i>FINITE_VOLUME, FV_ZONAL</i>
  */
  m_rntStartTimeStep = Context::getBasicProperty<MInt>("rntStartTimeStep", AT_, &m_rntStartTimeStep);

  /*! \page propertyPage1
    \section zonal
    <code>MInt FvZonalRTV::m_reconstructAverageFromNut</code>\n
    default = <code>False</code>\n
    Triggers reconstruction of LES Average values form nut after nonZonalRestart
    Keywords: <i>FINITE_VOLUME, FV_ZONAL</i>
  */
  m_reconstructAverageFromNut = false;
  if(Context::propertyExists("reconstructAverageFromNut")) {
    m_reconstructAverageFromNut =
        Context::getBasicProperty<MBool>("reconstructAverageFromNut", AT_, &m_reconstructAverageFromNut);
  }

  /*! \page propertyPage1
  \section zonal
  <code>MInt FvZonalRTV::m_reconstructNut</code>\n
  default = <code>False</code>\n
  Triggers reconstruction of LES Average values form nut after nonZonalRestart
  Keywords: <i>FINITE_VOLUME, FV_ZONAL</i>
*/
  m_reconstructNut = false;
  if(Context::propertyExists("reconstructNut")) {
    m_reconstructNut = Context::getBasicProperty<MBool>("reconstructNut", AT_, &m_reconstructNut);
  }

  /*
    turbulent intensity of the free stream
  */
  /*! \page propertyPage1
    \section zonal
    <code>MFloat FvZonalRTV::m_turbulentIntensity</code>\n
    default = <code>0</code>\n \n
    Turbulent intensity of the free stream.
    Keywords: <i>FINITE_VOLUME, FV_ZONAL</i>
  */
  m_turbulentIntensity = -1;
  if(Context::propertyExists("turbulentIntensity")) {
    m_turbulentIntensity = Context::getBasicProperty<MFloat>("turbulentIntensity", AT_, &m_turbulentIntensity);
  }

  /*
    Length scale correction when turbulent intensity of the free stream
    is reconstructed from data from literature based on the Reynolds number
    (activated if m_turbulentIntensity is not specified in property file)
  */
  m_tuLengthScaleCorrection = F1;
  if(Context::propertyExists("tuLengthScaleCorrection")) {
    m_tuLengthScaleCorrection =
        Context::getBasicProperty<MFloat>("tuLengthScaleCorrection", AT_, &m_tuLengthScaleCorrection);
  }
}


template <MInt nDim, class SysEqn>
void FvZonalRTV<nDim, SysEqn>::init() {
  TRACE();

  m_RANSSolverId = RANSSolver().m_solverId;
  m_LESSolverId = LESSolver().m_solverId;

  LESSolver().m_LESNoVarAverage = LESVarAverageData::noAvgVars;

  initRANSValues();
  initLESValues();

  determineZonalPositions();
}


template <MInt nDim, class SysEqn>
void FvZonalRTV<nDim, SysEqn>::determineZonalPositions() {
  TRACE();

  m_7902faceNormalDir = -1;
  m_7902Position = std::numeric_limits<MFloat>::infinity();

  // read cutoff positions to determine average and sponge cells
  if(Context::propertyExists("cutOffBndryIds")) {
    MInt propertyLength = Context::propertyLength("cutOffBndryIds", m_RANSSolverId);
    for(MInt i = 0; i < propertyLength; i++) {
      MInt bcId = Context::getSolverProperty<MFloat>("cutOffBndryIds", m_RANSSolverId, AT_, i);
      if(bcId == 7902 || bcId == 7905) {
        m_7902faceNormalDir = Context::getSolverProperty<MFloat>("cutOffDirections", m_RANSSolverId, AT_, i);
      }
    }
  }
  if(m_7902faceNormalDir != -1) {
    m_7902faceNormalDir = (m_7902faceNormalDir % 2 == 0) ? m_7902faceNormalDir + 1 : m_7902faceNormalDir - 1;
  }

  m_bcId7902 = -1;
  m_7902wallDir = 1;
  m_7902periodicDir = 2;

  if(Context::propertyExists("bc7902wallDir")) {
    m_7902wallDir = Context::getSolverProperty<MInt>("bc7902wallDir", m_RANSSolverId, AT_);
  }
  if(Context::propertyExists("bc7902periodicDir")) {
    m_7902periodicDir = Context::getSolverProperty<MInt>("bc7902periodicDir", m_RANSSolverId, AT_);
  }

  for(MInt bcId = 0; bcId < RANSSolver().m_fvBndryCnd->m_noCutOffBndryCndIds; bcId++) {
    if(RANSSolver().m_fvBndryCnd->m_cutOffBndryCndIds[bcId] == 7902) {
      m_bcId7902 = bcId;
    }
  }


  // determine real m_7909Position based on cell centers
  MFloat tmpPos7902 = ((m_7902faceNormalDir + 1) % 2 > 0) ? -std::numeric_limits<MFloat>::infinity()
                                                          : std::numeric_limits<MFloat>::infinity();
  MFloat compFactor7902 = ((m_7902faceNormalDir + 1) % 2 > 0) ? F1 : -F1;

  if(LESSolver().grid().isActive()) {
    for(MInt cellId = 0; cellId < a_noFvGridCellsLES(); cellId++) {
      MInt RANSId = convertIdParent(LESSolver(), RANSSolver(), cellId);
      if(RANSId != -1) {
        MFloat pos = LESSolver().a_coordinate(cellId, 0);
        if(compFactor7902 * pos > compFactor7902 * tmpPos7902) {
          tmpPos7902 = pos;
        }
      }
    }

    if(compFactor7902 > 0) {
      MPI_Allreduce(MPI_IN_PLACE, &tmpPos7902, 1, MPI_DOUBLE, MPI_MAX, LESSolver().mpiComm(), AT_, "MPI_IN_PLACE",
                    "tmpPos7902");
    } else {
      MPI_Allreduce(MPI_IN_PLACE, &tmpPos7902, 1, MPI_DOUBLE, MPI_MIN, LESSolver().mpiComm(), AT_, "MPI_IN_PLACE",
                    "tmpPos7902");
    }
    m_7902Position = tmpPos7902;
  }

  if(RANSSolver().grid().isActive()) {
    if(compFactor7902 > 0) {
      MPI_Allreduce(MPI_IN_PLACE, &tmpPos7902, 1, MPI_DOUBLE, MPI_MAX, RANSSolver().mpiComm(), AT_, "MPI_IN_PLACE",
                    "tmpPos7902");
    } else {
      MPI_Allreduce(MPI_IN_PLACE, &tmpPos7902, 1, MPI_DOUBLE, MPI_MIN, RANSSolver().mpiComm(), AT_, "MPI_IN_PLACE",
                    "tmpPos7902");
    }

    if(m_7902faceNormalDir != -1) {
      m_7902Position = tmpPos7902;
    }
  }

  m_averagePos = m_7902Position;
  m_averageDir = abs(m_7902wallDir + m_7902periodicDir - nDim);
  m_averageFaceDir = m_7902faceNormalDir;

  LESSolver().m_averagePos.push_back(m_averagePos);
  LESSolver().m_averageDir.push_back(m_averageDir);
  LESSolver().m_averageReconstructNut.push_back(m_reconstructNut);
}


template <MInt nDim, class SysEqn>
void FvZonalRTV<nDim, SysEqn>::finalizeCouplerInit() {
  TRACE();

  // finalizeLESAverage();
  determineNutReconstructionCells();

  if(globalTimeStep > m_rntStartTimeStep) {
    reconstructNutilde();
  }

  transferSolverData();
}


template <MInt nDim, class SysEqn>
void FvZonalRTV<nDim, SysEqn>::determineNutReconstructionCells() {
  TRACE();

  if(LESSolver().grid().isActive()) {
    m_noRntBcCells = 0;
    for(MInt LESId = 0; LESId < a_noFvGridCellsLES(); LESId++) {
      if(LESSolver().a_isHalo(LESId)) continue;

      MFloat halfCellLength = LESSolver().grid().halfCellLength(LESId);
      MFloat pos = LESSolver().a_coordinate(LESId, m_averageDir);

      if(approx(m_averagePos, pos, halfCellLength)) {
        m_rntBcCells.push_back(LESId);
        m_noRntBcCells++;
      }
    }
  }
}


template <MInt nDim, class SysEqn>
void FvZonalRTV<nDim, SysEqn>::postAdaptation() {
  TRACE();

  if(LESSolver().grid().isActive()) {
    // copy local vectors
    std::vector<MFloat> LESAverageCells_OLD_L;
    std::vector<std::vector<MFloat>> LESVarAverage_OLD;
    std::vector<MInt> rntBcCells_OLD;
    std::vector<MInt> LESAverageCells_OLD;

    for(MInt c = 0; c < (MInt)LESSolver().m_LESAverageCells.size(); c++) {
      LESAverageCells_OLD_L.push_back(LESSolver().m_LESAverageCells[c]);
    }

    for(MInt cc = 0; cc < m_noRntBcCells; cc++) {
      rntBcCells_OLD.push_back(m_rntBcCells[cc]);
    }
  }

  if(globalTimeStep > m_rntStartTimeStep) {
    reconstructNutilde();
  }

  transferSolverData();
}


template <MInt nDim, class SysEqn>
void FvZonalRTV<nDim, SysEqn>::preCouple(MInt /*step*/) {
  TRACE();

  if(LESSolver().grid().isActive()) {
    if(LESSolver().noNeighborDomains() > 0) {
      LESSolver().exchangeData(&LESSolver().a_variable(0, 0), noLESVariables());
      LESSolver().exchangeData(&LESSolver().a_oldVariable(0, 0), noLESVariables());
      LESSolver().exchangeData(&LESSolver().a_pvariable(0, 0), noLESVariables());
    }
  }

  if(RANSSolver().grid().isActive()) {
    if(RANSSolver().noNeighborDomains() > 0) {
      RANSSolver().exchangeData(&RANSSolver().a_variable(0, 0), noRANSVariables());
      RANSSolver().exchangeData(&RANSSolver().a_oldVariable(0, 0), noRANSVariables());
      RANSSolver().exchangeData(&RANSSolver().a_pvariable(0, 0), noRANSVariables());
    }
  }

  if(globalTimeStep % m_zonalTransferInterval == 0) {
    if(globalTimeStep > m_rntStartTimeStep) {
      reconstructNutilde();
    }
    transferSolverData();
  }
}


template <MInt nDim, class SysEqn>
void FvZonalRTV<nDim, SysEqn>::transferSolverData() {
  TRACE();

  // write LES data to RANS Solver and RANS data to LES Solver
  if(RANSSolver().grid().isActive()) {
    for(MInt var = 0; var < noRANSVariables(); var++) {
      for(MInt cellId = 0; cellId < a_noFvGridCellsRANS(); cellId++) {
        MInt LESId = convertIdParent(RANSSolver(), LESSolver(), cellId);

        ASSERT(cellId < (MInt)RANSSolver().m_LESValues[var].size(),
               "Trying to access data [" + to_string(var) + "][" + to_string(cellId)
                   + "] in LESSolver().m_LESVarAverage(RANS) with length "
                   + to_string(RANSSolver().m_LESValues[var].size())
                   + ", domainId: " + to_string(RANSSolver().domainId()));

        if(LESId != -1) {
          vector<MInt>::iterator findAvgId =
              find(LESSolver().m_LESAverageCells.begin(), LESSolver().m_LESAverageCells.end(), LESId);
          if(findAvgId != LESSolver().m_LESAverageCells.end()) {
            MInt avgId = distance(LESSolver().m_LESAverageCells.begin(), findAvgId);
            MInt var_ = var;
            // if(var == nDim + 2) var_ = nDim + 2 + 1;
            RANSSolver().m_LESValues[var][cellId] = LESSolver().m_LESVarAverage[var_][avgId];
          }
        }
      }
    }
  }

  // write RANS data to LES Solver
  if(LESSolver().grid().isActive()) {
    for(MInt cellId = 0; cellId < a_noFvGridCellsLES(); cellId++) {
      MInt RANSId = convertIdParent(LESSolver(), RANSSolver(), cellId);

      if(RANSId != -1) {
        for(MInt var = 0; var < noRANSVariables(); var++) {
          ASSERT(cellId < (MInt)LESSolver().m_RANSValues[var].size(),
                 "Trying to access data [" + to_string(var) + "][" + to_string(cellId)
                     + "] in m_RANSValues with length " + to_string(LESSolver().m_RANSValues[var].size())
                     + ", domainId: " + to_string(LESSolver().domainId()));

          LESSolver().m_RANSValues[var][cellId] = RANSSolver().a_pvariable(RANSId, var);
        }
      }
    }
  }
}


template <MInt nDim, class SysEqn>
void FvZonalRTV<nDim, SysEqn>::reconstructNutilde() {
  TRACE();

  const MInt N = RANSSolver().m_sysEqn.PV->N;

  if(LESSolver().grid().isActive()) {
    // reconstruct nuTilde
    //---------------------------------------------
    const MFloat rRe = F1 / LESSolver().m_sysEqn.m_Re0;

    // Turbulent intensity is calculated and used as a limiter for the turbulent kinetic energy
    /*data taken from:
      Brandt, Luca & Schlatter, Philipp. (2004). Transition in boundary layers subject to free-stream turbulence.
      Journal of Fluid Mechanics. 517. 167 - 198. 10.1017/S0022112004000941.
      ------------
      used to determine the expected freestream turbulence intensity
    */

    MFloat tuFreeStream = 0.05;

    if(abs(m_turbulentIntensity + 1) < 1e-16) {
      std::vector<std::pair<MFloat, MFloat>> tuFreeStreamData = {
          std::make_pair(31171.1711712, 0.0451546391753), std::make_pair(34774.7747748, 0.0428865979381),
          std::make_pair(38378.3783784, 0.0407216494845), std::make_pair(42702.7027027, 0.0385567010309),
          std::make_pair(51711.7117117, 0.0349484536082), std::make_pair(58558.5585586, 0.0327835051546),
          std::make_pair(66846.8468468, 0.0305154639175), std::make_pair(78738.7387387, 0.0278350515464),
          std::make_pair(93513.5135135, 0.0250515463918), std::make_pair(111891.891892, 0.0223711340206),
          std::make_pair(132432.432432, 0.0198969072165), std::make_pair(158738.738739, 0.0175257731959),
          std::make_pair(186126.126126, 0.0154639175258), std::make_pair(215675.675676, 0.0137113402062),
          std::make_pair(243783.783784, 0.0124742268041), std::make_pair(273693.693694, 0.0113402061856),
          std::make_pair(301801.801802, 0.0105154639175)};

      // interpolant freestream turbulence intensity from data
      MInt length = (MInt)tuFreeStreamData.size();
      for(MInt i = 0; i < length; i++) {
        if(tuFreeStreamData[i].first > LESSolver().m_Re * m_tuLengthScaleCorrection) {
          tuFreeStream = tuFreeStreamData[i - 1].second
                         + (tuFreeStreamData[i].second - tuFreeStreamData[i - 1].second)
                               / (tuFreeStreamData[i].first - tuFreeStreamData[i - 1].first)
                               * (LESSolver().m_Re * m_tuLengthScaleCorrection - tuFreeStreamData[i - 1].first);
          break;
        }
      }
      // extrapolat freestream turbulence intensity from data
      if(LESSolver().m_Re * m_tuLengthScaleCorrection < tuFreeStreamData[0].first) {
        tuFreeStream = tuFreeStreamData[0].second
                       + (tuFreeStreamData[1].second - tuFreeStreamData[0].second)
                             / (tuFreeStreamData[1].first - tuFreeStreamData[0].first)
                             * (LESSolver().m_Re * m_tuLengthScaleCorrection - tuFreeStreamData[0].first);
      }
      if(LESSolver().m_Re * m_tuLengthScaleCorrection > tuFreeStreamData[length - 1].first) {
        tuFreeStream = tuFreeStreamData[length - 2].second
                       + (tuFreeStreamData[length - 1].second - tuFreeStreamData[length - 2].second)
                             / (tuFreeStreamData[length - 1].first - tuFreeStreamData[length - 2].first)
                             * (LESSolver().m_Re * m_tuLengthScaleCorrection - tuFreeStreamData[length - 2].first);
      }

      tuFreeStream *= 0.95;

      m_log << "interpolation of freestream turbulence intensity: " << tuFreeStream << endl;

    } else {
      tuFreeStream = m_turbulentIntensity;
    }

    for(MInt c = 0; c < m_noRntBcCells; c++) {
      MInt cellId = m_rntBcCells[c];
      vector<MInt>::iterator findLESId =
          find(LESSolver().m_LESAverageCells.begin(), LESSolver().m_LESAverageCells.end(), cellId);
      if(findLESId == LESSolver().m_LESAverageCells.end()) TERMM(1, "Something went wrong " + to_string(cellId));
      MInt LESAvgId = distance(LESSolver().m_LESAverageCells.begin(), findLESId);

      // calculating turbulent kinetic energy k = 0.5(ui'ui')
      MFloat tuarg = 0.0;
      MFloat k = F0;
      MInt count_ = 0;
      MInt index_ = 0;
      MInt index1 = 0;
      MInt index2 = 0;
      for(MInt i = 0; i < nDim; i++) {
        for(MInt j = count_; j < nDim; j++) {
          index1 = index_ + noRANSVariables();
          if(i == j) {
            ASSERT(c < (MInt)LESSolver().m_LESVarAverage[index1].size(),
                   "Trying to access data [" + to_string(index1) + "][" + to_string(LESAvgId)
                       + "] in LESSolver().m_LESVarAverage with length "
                       + to_string(LESSolver().m_LESVarAverage[index1].size()));

            k += 0.5
                 * (LESSolver().m_LESVarAverage[index1][LESAvgId]
                    - pow(LESSolver().m_LESVarAverage[index2][LESAvgId], F2));
            tuarg += (LESSolver().m_LESVarAverage[index1][LESAvgId]
                      - pow(LESSolver().m_LESVarAverage[index2][LESAvgId], F2));
            index2++;
          }
          index_++;
        }
        count_++;
      }

      tuarg = max(tuarg, epss);

      MFloat tu = sqrt(tuarg / 3.0) / LESSolver().m_Ma;
      MFloat tufac = tanh((tu - tuFreeStream) / tuFreeStream);
      MFloat tulim = max(tu, tuFreeStream);
      tulim = tu / tulim;

      k = k * tufac * tulim;

      // calculating strain rate tensor Sij = 0.5(dui/dxj + duj/dxi)
      MFloat du[3][3]{{F0, F0, F0}, {F0, F0, F0}, {F0, F0, F0}};
      const MInt recData = LESSolver().a_reconstructionData(cellId);

      const MFloat u[3] = {LESSolver().m_LESVarAverage[0][LESAvgId], LESSolver().m_LESVarAverage[1][LESAvgId],
                           LESSolver().m_LESVarAverage[2][LESAvgId]};

      for(MInt nghbr = 0; nghbr < LESSolver().a_noReconstructionNeighbors(cellId); nghbr++) {
        const MInt recNghbrId = LESSolver().a_reconstructionNeighborId(cellId, nghbr);
        const MFloat recConst_x = LESSolver().m_reconstructionConstants[nDim * (recData + nghbr) + 0];
        const MFloat recConst_y = LESSolver().m_reconstructionConstants[nDim * (recData + nghbr) + 1];
        const MFloat recConst_z = LESSolver().m_reconstructionConstants[nDim * (recData + nghbr) + 2];

        vector<MInt>::iterator findRecNghbrId =
            find(LESSolver().m_LESAverageCells.begin(), LESSolver().m_LESAverageCells.end(), recNghbrId);

        if(findRecNghbrId == LESSolver().m_LESAverageCells.end())
          TERMM(1, "Something went wrong in nghbr " + to_string(recNghbrId));

        MInt LESRecNghbrId = distance(LESSolver().m_LESAverageCells.begin(), findRecNghbrId);

        for(MInt dim = 0; dim < nDim; dim++) {
          ASSERT(c < (MInt)LESSolver().m_LESVarAverage[dim].size(),
                 "Trying to access data [" + to_string(dim) + "][" + to_string(LESRecNghbrId)
                     + "] in LESSolver().m_LESVarAverage with length "
                     + to_string(LESSolver().m_LESVarAverage[dim].size()));

          const MFloat delta_u = LESSolver().m_LESVarAverage[dim][LESRecNghbrId] - u[dim];
          if(abs(LESSolver().m_LESVarAverage[dim][LESRecNghbrId]) > epss
             || abs(u[dim] / LESSolver().m_LESVarAverage[dim][LESRecNghbrId]) > 100) {
            du[dim][0] += recConst_x * delta_u;
            du[dim][1] += recConst_y * delta_u;
            du[dim][2] += recConst_z * delta_u;
          }
        }
      }

      MFloat SijSij = F0;
      for(MInt d1 = 0; d1 < nDim; d1++) {
        for(MInt d2 = 0; d2 < nDim; d2++) {
          MFloat sij = 0.5 * (du[d1][d2] + du[d2][d1]);
          SijSij += sij * sij;
        }
      }

      const MFloat SijSijLim = (4.0 * pow(LESSolver().m_Ma, F2)) * 0.0001;
      SijSij = max(SijSij, SijSijLim);
      const MFloat c_mu = 0.09;
      MFloat omega = 1 / sqrt(c_mu) * sqrt(2 * SijSij);
      omega = max(omega, epss);

      // calculating nut
      MFloat nuTurb = max(k / (omega * rRe), LESSolver().m_nuTildeInfinity);
      MFloat nuTilde = LESSolver().m_nuTildeInfinity;
      if(std::isnan(omega) || tulim < 1) {
        nuTurb = LESSolver().m_nuTildeInfinity;
      }

      if(approx(tulim, F1, eps)) {
        const MFloat rho = LESSolver().a_pvariable(cellId, LESSolver().m_sysEqn.PV->RHO);
        const MFloat p = LESSolver().a_pvariable(cellId, LESSolver().m_sysEqn.PV->P);
        const MFloat T = LESSolver().m_gamma * p / rho;
        const MFloat mue = (T * sqrt(T) * LESSolver().m_sutherlandPlusOne) / (T + LESSolver().m_sutherlandConstant);
        const MFloat nuLaminar = mue / rho;

        // Iteration of nuTilde (using Newton, nu_(n+1) = nu_n - f(nu_n)/f'(nu_n))
        MInt NewtonIter = 50;
        MFloat nut = nuTurb;
        MFloat nu_new = nut;
        MFloat nu_old = 0;
        MFloat func_n = 0.0;
        MFloat derv_n = 0.0;
        MInt icount = 0;

        while(abs(nu_new - nu_old) > eps && icount < NewtonIter) {
          nu_old = nu_new;

          MFloat chi = nu_old / nuLaminar;
          MFloat chi3 = pow(chi, 3.0);
          MFloat cv13 = RM_SA_DV::cv1to3;
          MFloat fv1 = chi3 / (chi3 + cv13);

          // current function value
          func_n = nu_old * fv1 - nut;
          // derivative of function
          derv_n = fv1 + 3.0 * fv1 - 3.0 * fv1 * fv1;
          // newton iteration formula
          nu_new = nu_old - func_n / derv_n;

          icount++;
        }

        nuTilde = nu_new;
      }

      LESSolver().m_LESVarAverage[N][LESAvgId] = nuTilde;
    }
  }
}


template <MInt nDim, class SysEqn>
void FvZonalRTV<nDim, SysEqn>::reconstructAverageFromNut() {
  TRACE();

  if(LESSolver().grid().isActive()) {
    for(MInt c = 0; c < (MInt)LESSolver().m_LESAverageCells.size(); c++) {
      MInt cellId = LESSolver().m_LESAverageCells[c];

      MInt RANSId = convertIdParent(LESSolver(), RANSSolver(), cellId);
      if(RANSId > -1) {
        MFloat nuTilde = RANSSolver().a_pvariable(RANSId, RANSSolver().m_sysEqn.PV->N);
        MFloat add = F0;
        for(MInt i = 0; i < nDim; i++) {
          add += F2B3 * RANSSolver().a_slope(RANSId, i, i);
        }

        MInt count = 0;
        MInt index = 0;
        for(MInt i = 0; i < nDim; i++) {
          for(MInt j = count; j < nDim; j++) {
            MInt indexAvg = index + noRANSVariables();

            ASSERT(c < (MInt)LESSolver().m_LESVarAverage[indexAvg].size(),
                   "Trying to access data [" + to_string(indexAvg) + "][" + to_string(c)
                       + "] in LESSolver().m_LESVarAverage with length "
                       + to_string(LESSolver().m_LESVarAverage[indexAvg].size()));

            if(i == j) {
              LESSolver().m_LESVarAverage[indexAvg][c] =
                  LESSolver().m_LESVarAverage[i][c] * LESSolver().m_LESVarAverage[j][c]
                  - nuTilde * (RANSSolver().a_slope(RANSId, i, j) + RANSSolver().a_slope(RANSId, j, i) - add);
            } else {
              LESSolver().m_LESVarAverage[indexAvg][c] =
                  LESSolver().m_LESVarAverage[i][c] * LESSolver().m_LESVarAverage[j][c]
                  - nuTilde * (RANSSolver().a_slope(RANSId, i, j) + RANSSolver().a_slope(RANSId, j, i));
            }

            index++;
          }
          count++;
        }
      }
    }
  }
}


template class FvZonalRTV<2, FvSysEqnRANS<2, RANSModelConstants<RANS_SA_DV>>>;
template class FvZonalRTV<2, FvSysEqnRANS<2, RANSModelConstants<RANS_FS>>>;
template class FvZonalRTV<3, FvSysEqnRANS<3, RANSModelConstants<RANS_SA_DV>>>;
template class FvZonalRTV<3, FvSysEqnRANS<3, RANSModelConstants<RANS_FS>>>;
