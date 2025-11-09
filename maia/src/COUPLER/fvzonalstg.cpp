// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "fvzonalstg.h"

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
FvZonalSTG<nDim, SysEqn>::FvZonalSTG(const MInt couplingId, RANS* R, LES* L)
  : Coupling(couplingId), FvZonal<nDim, SysEqn>(couplingId, R, L) {
  //  Coupling(couplingId), m_RANSSolver(R), m_LESSolver(L) {

  m_log << "Create Zonal coupler for RANS Solver (" << RANSSolver().m_solverId << ") and LES Solver ("
        << LESSolver().m_solverId << ")" << endl;

  /*! \page propertyPage1
    \section zonal
    <code>MFloat FvZonalSTG::m_uvRANSFactor</code>\n
    default = <code>0</code>\n
    scale uv of RANS solution
    Keywords: <i>FINITE_VOLUME, FV_ZONA_STGL</i>
  */
  m_uvRANSFactor = F1;
  if(Context::propertyExists("uvRANSFactor")) {
    m_uvRANSFactor = Context::getBasicProperty<MFloat>("uvRANSFactor", AT_, &m_uvRANSFactor);
  }
}


template <MInt nDim, class SysEqn>
void FvZonalSTG<nDim, SysEqn>::init() {
  TRACE();

  m_RANSSolverId = RANSSolver().m_solverId;
  m_LESSolverId = LESSolver().m_solverId;

  LESSolver().m_LESNoVarAverage = LESVarAverageData::noAvgVars;

  // necessary for STG Method
  initRANSValues();
  initLESValues();

  determineZonalPositions();

  if(m_cylindricCommunication) {
    initCylinderExchange();
    // calcLESSectorAverage();
    calcRANSSectorValues();
  }

  // necessary for boundary condition which is called in finalizeInitSolver
  transferSolverData();
}


template <MInt nDim, class SysEqn>
void FvZonalSTG<nDim, SysEqn>::finalizeCouplerInit() {
  TRACE();

  if(m_cylindricCommunication) {
    initCylinderExchange();
  }

  if(m_cylindricCommunication) {
    // calcLESSectorAverage();
    calcRANSSectorValues();
  }

  transferSolverData();
}


template <MInt nDim, class SysEqn>
void FvZonalSTG<nDim, SysEqn>::finalizeAdaptation(const MInt solverId) {
  TRACE();

  if(solverId != m_RANSSolverId || solverId != m_LESSolverId) return;

  transferSolverData();
}


template <MInt nDim, class SysEqn>
void FvZonalSTG<nDim, SysEqn>::balancePost() {
  TRACE();

  transferSolverData();
}


template <MInt nDim, class SysEqn>
void FvZonalSTG<nDim, SysEqn>::preCouple(MInt /*step*/) {
  TRACE();

  if(globalTimeStep % m_zonalTransferInterval == 0 || globalTimeStep == LESSolver().m_restartTimeStep + 1) {
    if(m_cylindricCommunication) {
      // initCylinderExchange();
      // calcLESSectorAverage();
      calcRANSSectorValues();
    }

    transferSolverData();
  }
}


template <MInt nDim, class SysEqn>
void FvZonalSTG<nDim, SysEqn>::determineZonalPositions() {
  TRACE();

  LESSolver().m_averagePos.clear();
  LESSolver().m_averageDir.clear();

  // if(m_preliminarySTGSponge) {
  //   m_7901Position = m_preliminarySTGSpongePosition;
  //   m_averagePos = m_7901Position;
  //   m_averageDir = m_preliminarySTGSpongeDirection;
  // } else {

  m_7901faceNormalDir = -1;
  m_7909faceNormalDir = -1;
  m_7901Position = std::numeric_limits<MFloat>::infinity();
  m_7909Position = std::numeric_limits<MFloat>::infinity();

  // read cutoff positions to determine average and sponge cells
  if(Context::propertyExists("cutOffBndryIds")) {
    MInt propertyLength = Context::propertyLength("cutOffBndryIds", m_RANSSolverId);
    for(MInt i = 0; i < propertyLength; i++) {
      MInt bcId = Context::getSolverProperty<MFloat>("cutOffBndryIds", m_RANSSolverId, AT_, i);
      if(bcId == 7901) {
        m_7901faceNormalDir = Context::getSolverProperty<MFloat>("cutOffDirections", m_RANSSolverId, AT_, i);
      }
    }
  }
  if(m_7901faceNormalDir != -1) {
    m_7901faceNormalDir = (m_7901faceNormalDir % 2 == 0) ? m_7901faceNormalDir + 1 : m_7901faceNormalDir - 1;
  }

  if(Context::propertyExists("cutOffBndryIds")) {
    MInt propertyLength = Context::propertyLength("cutOffBndryIds", m_LESSolverId);
    for(MInt i = 0; i < propertyLength; i++) {
      MInt bcId = Context::getSolverProperty<MFloat>("cutOffBndryIds", m_LESSolverId, AT_, i);
      if(bcId == 7909) {
        m_7909faceNormalDir = Context::getSolverProperty<MFloat>("cutOffDirections", m_LESSolverId, AT_, i);
      }
    }
  }
  if(m_7909faceNormalDir != -1) {
    m_7909faceNormalDir = (m_7909faceNormalDir % 2 == 0) ? m_7909faceNormalDir + 1 : m_7909faceNormalDir - 1;
  }

  if(m_STGSponge && m_7901faceNormalDir == -1) {
    TERMM(1, "m_STGSponge true, bc7901 has to be set to determine spongeCells!");
  }

  m_7901wallDir = 1;
  m_7901periodicDir = 2;

  if(Context::propertyExists("bc7901wallDir")) {
    m_7901wallDir = Context::getSolverProperty<MInt>("bc7901wallDir", m_RANSSolverId, AT_);
  }
  if(Context::propertyExists("bc7901periodicDir")) {
    m_7901periodicDir = Context::getSolverProperty<MInt>("bc7901periodicDir", m_RANSSolverId, AT_);
  }

  // determine real m_7909Position based on cell centers
  MFloat tmpPos7901 = ((m_7901faceNormalDir + 1) % 2 > 0) ? -std::numeric_limits<MFloat>::infinity()
                                                          : std::numeric_limits<MFloat>::infinity();
  MFloat tmpPos7909 = ((m_7909faceNormalDir + 1) % 2 > 0) ? -std::numeric_limits<MFloat>::infinity()
                                                          : std::numeric_limits<MFloat>::infinity();
  MFloat compFactor7901 = ((m_7901faceNormalDir + 1) % 2 > 0) ? F1 : -F1;
  MFloat compFactor7909 = ((m_7909faceNormalDir + 1) % 2 > 0) ? F1 : -F1;

  if(LESSolver().grid().isActive()) {
    for(MInt cellId = 0; cellId < a_noFvGridCellsLES(); cellId++) {
      MInt RANSId = convertIdParent(LESSolver(), RANSSolver(), cellId);
      if(RANSId != -1) {
        MFloat pos = LESSolver().a_coordinate(cellId, 0);
        if(compFactor7901 * pos > compFactor7901 * tmpPos7901) {
          tmpPos7901 = pos;
        }
        if(compFactor7909 * pos > compFactor7909 * tmpPos7909) {
          tmpPos7909 = pos;
        }
      }
    }

    if(compFactor7901 > 0) {
      MPI_Allreduce(MPI_IN_PLACE, &tmpPos7901, 1, MPI_DOUBLE, MPI_MAX, LESSolver().mpiComm(), AT_, "MPI_IN_PLACE",
                    "tmpPos7901");
    } else {
      MPI_Allreduce(MPI_IN_PLACE, &tmpPos7901, 1, MPI_DOUBLE, MPI_MIN, LESSolver().mpiComm(), AT_, "MPI_IN_PLACE",
                    "tmpPos7901");
    }
    if(compFactor7909 > 0) {
      MPI_Allreduce(MPI_IN_PLACE, &tmpPos7909, 1, MPI_DOUBLE, MPI_MAX, LESSolver().mpiComm(), AT_, "MPI_IN_PLACE",
                    "tmpPos7909");
    } else {
      MPI_Allreduce(MPI_IN_PLACE, &tmpPos7909, 1, MPI_DOUBLE, MPI_MIN, LESSolver().mpiComm(), AT_, "MPI_IN_PLACE",
                    "tmpPos7909");
    }

    if(m_7901faceNormalDir != -1) {
      m_7901Position = tmpPos7901;
    }
    if(m_7909faceNormalDir != -1) {
      m_7909Position = tmpPos7909;
    }
  }

  if(RANSSolver().grid().isActive()) {
    if(compFactor7901 > 0) {
      MPI_Allreduce(MPI_IN_PLACE, &tmpPos7901, 1, MPI_DOUBLE, MPI_MAX, RANSSolver().mpiComm(), AT_, "MPI_IN_PLACE",
                    "tmpPos7901");
    } else {
      MPI_Allreduce(MPI_IN_PLACE, &tmpPos7901, 1, MPI_DOUBLE, MPI_MIN, RANSSolver().mpiComm(), AT_, "MPI_IN_PLACE",
                    "tmpPos7901");
    }
    if(compFactor7909 > 0) {
      MPI_Allreduce(MPI_IN_PLACE, &tmpPos7909, 1, MPI_DOUBLE, MPI_MAX, RANSSolver().mpiComm(), AT_, "MPI_IN_PLACE",
                    "tmpPos7909");
    } else {
      MPI_Allreduce(MPI_IN_PLACE, &tmpPos7909, 1, MPI_DOUBLE, MPI_MIN, RANSSolver().mpiComm(), AT_, "MPI_IN_PLACE",
                    "tmpPos7909");
    }
    MPI_Allreduce(MPI_IN_PLACE, &tmpPos7909, 1, MPI_DOUBLE, MPI_MAX, RANSSolver().mpiComm(), AT_, "MPI_IN_PLACE",
                  "tmpPos7909");

    if(m_7901faceNormalDir != -1) {
      m_7901Position = tmpPos7901;
    }
    if(m_7909faceNormalDir != -1) {
      m_7909Position = tmpPos7909;
    }
  }

  m_averagePos = m_7901Position;
  m_averageDir = abs(m_7901wallDir + m_7901periodicDir - nDim);

  if(!std::isinf(m_7901Position)) {
    // positions to find cells for LES average
    LESSolver().m_7901Position = m_7901Position;
    LESSolver().m_averagePos.push_back(m_averagePos);
    LESSolver().m_averageDir.push_back(m_averageDir);
    LESSolver().m_averageReconstructNut.push_back(false);
  }
}


template <MInt nDim, class SysEqn>
void FvZonalSTG<nDim, SysEqn>::transferSolverData() {
  TRACE();

  // write LES data to RANS Solver
  if(RANSSolver().grid().isActive()) {
    if(!m_cylindricCommunication) {
      for(MInt cellId = 0; cellId < a_noFvGridCellsRANS(); cellId++) {
        MInt LESId = convertIdParent(RANSSolver(), LESSolver(), cellId);
        if(LESId != -1) {
          vector<MInt>::iterator findAvgId =
              find(LESSolver().m_LESAverageCells.begin(), LESSolver().m_LESAverageCells.end(), LESId);
          if(findAvgId != LESSolver().m_LESAverageCells.end()) {
            MInt avgId = distance(LESSolver().m_LESAverageCells.begin(), findAvgId);
            for(MInt var = 0; var < noRANSVariables(); var++) {
              ASSERT(cellId < (MInt)RANSSolver().m_LESValues[var].size(),
                     "Trying to access data [" + to_string(var) + "][" + to_string(cellId)
                         + "] in LESSolver().m_LESVarAverage(RANS) with length "
                         + to_string(RANSSolver().m_LESValues[var].size())
                         + ", domainId: " + to_string(RANSSolver().domainId()));

              RANSSolver().m_LESValues[var][cellId] = LESSolver().m_LESVarAverage[var][avgId];
            }
          }
        }
      }
      //}
    } else {
      if(m_cylindricCommunication && m_cylinderExchangeIdsOffset > -1) {
        for(MInt exchangeIndex = 0; exchangeIndex < m_noRANSExchangeCells; exchangeIndex++) {
          MInt exchangeIndexOffset = exchangeIndex + m_cylinderExchangeIdsOffset;
          MInt RANSId = m_globalCylinderExchangeIds[exchangeIndexOffset];
          if(RANSSolver().a_isBndryGhostCell(RANSId)) continue;
          for(MInt var = 0; var < noRANSVariables(); var++) {
            ASSERT(RANSId < (MInt)RANSSolver().m_LESValues[var].size(),
                   "Trying to access data [" + to_string(var) + "][" + to_string(RANSId)
                           + "] in LESSolver().m_LESVarAverage(RANS) with length "
                           + to_string(RANSSolver().m_LESValues[var].size())
                           + ", domainId: " + to_string(RANSSolver().domainId()) + " | "
                       << RANSSolver().c_noCells());

            RANSSolver().m_LESValues[var][RANSId] =
                m_globalCylinderLESExchangeValues[exchangeIndexOffset * LESSolver().m_LESNoVarAverage + var];
          }
        }
      }
    }
  }

  if(!m_cylindricCommunication) {
    // write RANS data to LES Solver
    if(LESSolver().grid().isActive()) {
      for(MInt LESId = 0; LESId < a_noFvGridCellsLES(); LESId++) {
        MInt RANSId = convertIdParent(LESSolver(), RANSSolver(), LESId);

        if(RANSId != -1) {
          for(MInt var = 0; var < noRANSVariables(); var++) {
            ASSERT(LESId < (MInt)LESSolver().m_RANSValues[var].size(),
                   "Trying to access data [" + to_string(var) + "][" + to_string(LESId)
                           + "] in m_RANSValues with length " + to_string(LESSolver().m_RANSValues[var].size())
                           + ", domainId: " + to_string(LESSolver().domainId()) + " | "
                       << LESSolver().c_noCells());

            LESSolver().m_RANSValues[var][LESId] = RANSSolver().a_pvariable(RANSId, var);
          }

          if(m_STGSponge) {
            // calculate m_uvRans
            const MFloat fre = 1.0 / RANSSolver().sysEqn().m_Re0;
            const MFloat nu_t = RANSSolver().a_pvariable(RANSId, RANSSolver().m_sysEqn.PV->N);

            MFloat du[3][3]{{F0, F0, F0}, {F0, F0, F0}, {F0, F0, F0}};
            const MInt recData = RANSSolver().a_reconstructionData(RANSId);
            const MFloat u[3] = {RANSSolver().a_pvariable(RANSId, RANSSolver().m_sysEqn.PV->U),
                                 RANSSolver().a_pvariable(RANSId, RANSSolver().m_sysEqn.PV->V),
                                 RANSSolver().a_pvariable(RANSId, RANSSolver().m_sysEqn.PV->W)};

            for(MInt nghbr = 0; nghbr < RANSSolver().a_noReconstructionNeighbors(RANSId); nghbr++) {
              const MInt recNghbrId = RANSSolver().a_reconstructionNeighborId(RANSId, nghbr);
              if(recNghbrId > -1) {
                const MFloat recConst_x = RANSSolver().m_reconstructionConstants[nDim * (recData + nghbr) + 0];
                const MFloat recConst_y = RANSSolver().m_reconstructionConstants[nDim * (recData + nghbr) + 1];
                const MFloat recConst_z = RANSSolver().m_reconstructionConstants[nDim * (recData + nghbr) + 2];
                for(MInt dim = 0; dim < nDim; ++dim) {
                  const MFloat delta_u =
                      RANSSolver().a_pvariable(recNghbrId, RANSSolver().m_sysEqn.PV->VV[dim]) - u[dim];
                  du[dim][0] += recConst_x * delta_u;
                  du[dim][1] += recConst_y * delta_u;
                  du[dim][2] += recConst_z * delta_u;
                }
              }
            }
            MFloat sij[3][3]{{F0, F0, F0}, {F0, F0, F0}, {F0, F0, F0}};
            MFloat SijSij = F0;
            for(MInt d1 = 0; d1 < nDim; d1++) {
              for(MInt d2 = 0; d2 < nDim; d2++) {
                sij[d1][d2] = 0.5 * (du[d1][d2] + du[d2][d1]);
                SijSij += sij[d1][d2] * sij[d1][d2];
              }
            }
            const MFloat sr1 = (sij[0][1] + sij[1][0]) * (sij[0][1] + sij[1][0]);
            const MFloat sr2 = (sij[1][2] + sij[2][1]) * (sij[1][2] + sij[2][1]);
            const MFloat sr3 = (sij[0][2] + sij[2][0]) * (sij[0][2] + sij[2][0]);
            const MFloat srt = std::max(sqrt(sr1 + sr2 + sr3), epss);
            const MFloat rr1 = sqrt(sr1) / srt;

            MFloat uvRans = -sqrt(2.0 * SijSij) * rr1 * nu_t * fre;
            LESSolver().m_RANSValues[noRANSVariables()][LESId] = uvRans;
          }
        }
      }
    }
  } else {
    cylinderExchange();
  }
}

template <MInt nDim, class SysEqn>
void FvZonalSTG<nDim, SysEqn>::transferSpongeData() {
  TRACE();

  // write STG Sponge data to LES block
  if(LESSolver().grid().isActive()) {
    if(m_cylindricCommunication) {
      for(MInt cellId = 0; cellId < a_noFvGridCellsLES(); cellId++) {
        // write value to cell
        if(LESSolver().a_isHalo(cellId)) continue;
        MFloat y = LESSolver().a_coordinate(cellId, 1);
        MFloat z = LESSolver().a_coordinate(cellId, 2);
        MFloat sector = maia::math::getSector(y, z, m_azimuthalAngle);

        for(MInt var = 0; var < nDim + 3; var++) {
          ASSERT(cellId < (MInt)LESSolver().m_STGSpongeFactor[var].size(),
                 "Trying to access data [" + to_string(var) + "][" + to_string(cellId)
                     + "] in m_STGSpongeFactor with length " + to_string(LESSolver().m_STGSpongeFactor[var].size())
                     + ", domainId: " + to_string(LESSolver().domainId()));

          if(m_periodicSpongeInterpolationIndex[cellId].size() == 1) {
            MInt interpolationIndex = m_periodicSpongeInterpolationIndex[cellId][0];
            MInt exchangeIndex = m_periodicSpongeCylinderExchangeIndex[interpolationIndex];
            MFloat data = F0;
            if(var < nDim) {
              data = LESSolver().a_pvariable(cellId, var)
                     - m_globalCylinderLESExchangeValues[exchangeIndex * LESSolver().m_LESNoVarAverage + var];
            }
            if(var == nDim) {
              data = m_uvRans[interpolationIndex];
            }
            if(var == nDim + 1) {
              data = alpha * m_uvErr[interpolationIndex] + beta * m_uvInt[interpolationIndex];
            }
            if(var == nDim + 2) {
              MFloat uv_LES =
                  m_globalCylinderLESExchangeValues[exchangeIndex * LESSolver().m_LESNoVarAverage + noLESVariables()
                                                    + 1]
                  - m_globalCylinderLESExchangeValues[exchangeIndex * LESSolver().m_LESNoVarAverage + 0]
                        * m_globalCylinderLESExchangeValues[exchangeIndex * LESSolver().m_LESNoVarAverage + 1];
              data = uv_LES;
            }
            LESSolver().m_STGSpongeFactor[var][cellId] = data;
          }

          if(m_periodicSpongeInterpolationIndex[cellId].size() > 1) {
            // Interpolation - Least square method
            MFloat n = F0;
            MFloat Eyi = F0;
            MFloat Eyi2 = F0;
            MFloat Ezi = F0;
            MFloat Ezi2 = F0;
            MFloat Eyizi = F0;
            MFloat Eui = F0;
            MFloat Euiyi = F0;
            MFloat Euizi = F0;

            for(MInt i = 0; i < (MInt)m_periodicSpongeInterpolationIndex[cellId].size(); i++) {
              // MBool tmp = false;
              MInt interpolationIndex = m_periodicSpongeInterpolationIndex[cellId][i];
              MInt exchangeIndex = m_periodicSpongeCylinderExchangeIndex[interpolationIndex];

              MFloat yExchange = m_globalCylinderExchangeLocations[exchangeIndex * (nDim + 1) + 1];
              MFloat zExchange = m_globalCylinderExchangeLocations[exchangeIndex * (nDim + 1) + 2];
              MFloat rotationAngle =
                  -((sector - m_globalCylinderExchangeLocations[exchangeIndex * (nDim + 1) + 3]) * m_azimuthalAngle)
                  / 180.0 * PI;
              MFloat yRotated = cos(rotationAngle) * yExchange + sin(rotationAngle) * zExchange;
              MFloat zRotated = -sin(rotationAngle) * yExchange + cos(rotationAngle) * zExchange;

              MFloat delta_y = LESSolver().a_coordinate(cellId, 1) - yRotated;
              MFloat delta_z = LESSolver().a_coordinate(cellId, 2) - zRotated;


              n += 1;
              Eyi += delta_y;
              Eyi2 += POW2(delta_y);
              Ezi += delta_z;
              Ezi2 += POW2(delta_z);
              Eyizi += delta_y * delta_z;

              MFloat data = F0;
              if(var < nDim) {
                data = LESSolver().a_pvariable(cellId, var)
                       - m_globalCylinderLESExchangeValues[exchangeIndex * LESSolver().m_LESNoVarAverage + var];
              }
              if(var == nDim) {
                data = m_uvRans[interpolationIndex];
              }
              if(var == nDim + 1) {
                data = alpha * m_uvErr[interpolationIndex] + beta * m_uvInt[interpolationIndex];
              }
              if(var == nDim + 2) {
                MFloat uv_LES =
                    m_globalCylinderLESExchangeValues[exchangeIndex * LESSolver().m_LESNoVarAverage + noLESVariables()
                                                      + 1]
                    - m_globalCylinderLESExchangeValues[exchangeIndex * LESSolver().m_LESNoVarAverage + 0]
                          * m_globalCylinderLESExchangeValues[exchangeIndex * LESSolver().m_LESNoVarAverage + 1];
                data = uv_LES;
              }
              Eui += data;
              Euiyi += data * delta_y;
              Euizi += data * delta_z;
            }

            MFloat detA =
                n * Eyi2 * Ezi2 + 2 * Eyi * Eyizi * Ezi - (Ezi * Eyi2 * Ezi + n * Eyizi * Eyizi + Ezi2 * Eyi * Eyi);

            if(detA < 1e-12) {
              MInt interpolationIndex = m_periodicSpongeInterpolationIndex[cellId][0];
              LESSolver().m_STGSpongeFactor[var][cellId] =
                  m_globalCylinderRANSExchangeValues[interpolationIndex * m_noRANSCylinderExchangeVariables + var];
            } else {
              LESSolver().m_STGSpongeFactor[var][cellId] =
                  1 / detA
                  * (Eui * (Eyi2 * Ezi2 - Eyizi * Eyizi) + Euiyi * (Ezi * Eyizi - Eyi * Ezi2)
                     + Euizi * (Eyi * Eyizi - Eyi2 * Ezi));
            }
          }
        }
      }
    } else {
      for(MInt LESId = 0; LESId < a_noFvGridCellsLES(); LESId++) {
        MInt RANSId = convertIdParent(LESSolver(), RANSSolver(), LESId);
        if(RANSId != -1) {
          for(MInt i = 0; i < LESSolver().m_globalNoSpongeLocations; i++) {
            if(abs(LESSolver().m_globalSpongeLocations[i].first - RANSSolver().a_coordinate(RANSId, m_7901wallDir))
               < eps) {
              LESSolver().m_uvRans[i] = m_uvRans[i];
            }
          }
        }
      }
    }
  }
}


/** \brief Create mapping for LES cells to equivalent RANS cells
 *    \author Jannik Borgelt
 */
template <MInt nDim, class SysEqn>
void FvZonalSTG<nDim, SysEqn>::initCylinderExchange() {
  TRACE();

  if(LESSolver().grid().isActive()) {
    std::vector<MFloat> RANSSectors;
    RANSSectors.push_back(std::numeric_limits<MFloat>::max());
    RANSSectors.push_back(-std::numeric_limits<MFloat>::max());
    std::vector<MFloat> RANSSectorLimits;
    RANSSectorLimits.push_back(std::numeric_limits<MFloat>::max());
    RANSSectorLimits.push_back(-std::numeric_limits<MFloat>::max());
    std::vector<MInt> cylinderExchangeIds;
    std::vector<MFloat> cylinderExchangeLocations;

    mAlloc(m_RANSSectors, 2, "m_RANSSectors", F0, FUN_);
    mAlloc(m_RANSSectorLimits, 2, "m_RANSSectorLimits", F0, FUN_);

    MInt noRANSExchangeCells = 0;
    for(MInt RANSId = 0; RANSId < a_noFvGridCellsRANS(); RANSId++) {
      MFloat halfCellLength = RANSSolver().grid().halfCellLength(RANSId);
      MFloat pos = RANSSolver().a_coordinate(RANSId, 0);
      if(approx(m_7901Position, pos, halfCellLength) || approx(m_7909Position, pos, halfCellLength)) {
        // MInt RANSId = convertIdParent(RANSSolver(), RANSSolver(), RANSId);
        // if(RANSId > -1){
        if(RANSSolver().c_isLeafCell(RANSId) && !RANSSolver().a_isHalo(RANSId)) {
          MFloat x = RANSSolver().a_coordinate(RANSId, 0);
          MFloat y = RANSSolver().a_coordinate(RANSId, 1);
          MFloat z = RANSSolver().a_coordinate(RANSId, 2);
          MFloat sector = maia::math::getSector(y, z, m_azimuthalAngle);
          MFloat angle = maia::math::getAngle(y, z);

          if(sector < RANSSectors[0]) {
            RANSSectors[0] = sector;
          }
          if(sector > RANSSectors[1]) {
            RANSSectors[1] = sector;
          }
          if(angle < RANSSectorLimits[0]) {
            RANSSectorLimits[0] = angle;
          }
          if(angle > RANSSectorLimits[1]) {
            RANSSectorLimits[1] = angle;
          }

          noRANSExchangeCells++;
          cylinderExchangeIds.push_back(RANSId);
          cylinderExchangeLocations.push_back(x);
          cylinderExchangeLocations.push_back(y);
          cylinderExchangeLocations.push_back(z);
          cylinderExchangeLocations.push_back(maia::math::getSector(y, z, m_azimuthalAngle));
          //   }
        }
      }
    }

    if(noRANSExchangeCells > 0) {
      MInt noRANSExchangeCells_ = noRANSExchangeCells;
      // add ghostCells for better boundary interpolation
      for(MInt i = 0; i < noRANSExchangeCells_; i++) {
        MInt RANSId = cylinderExchangeIds[i];
        MInt noNghbrIds = RANSSolver().a_noReconstructionNeighbors(RANSId);
        for(MInt n = 0; n < noNghbrIds; n++) {
          MInt nghbrId = RANSSolver().a_reconstructionNeighborId(RANSId, n);
          if(RANSSolver().a_isBndryGhostCell(nghbrId)) {
            MFloat x = RANSSolver().a_coordinate(nghbrId, 0);
            MFloat y = RANSSolver().a_coordinate(nghbrId, 1);
            MFloat z = RANSSolver().a_coordinate(nghbrId, 2);
            noRANSExchangeCells++;
            cylinderExchangeIds.push_back(nghbrId);
            cylinderExchangeLocations.push_back(x);
            cylinderExchangeLocations.push_back(y);
            cylinderExchangeLocations.push_back(z);
            cylinderExchangeLocations.push_back(maia::math::getSector(y, z, m_azimuthalAngle));
          }
        }
      }
    }

    //--------------------------------------------------------------------------
#ifdef _MB_DEBUG_
    if(noRANSExchangeCells > 0) {
      stringstream f;
      f.clear();
      f << "globalCylinderExchange_" << to_string(RANSSolver().domainId()) << ".txt";
      MString fname = f.str();
      MString header = "#id isHalo x1 x2 x3 sector";
      ofstream data;
      data.precision(16);
      data.open(fname);
      data << header << endl;
      for(MInt i = 0; i < noRANSExchangeCells; i++) {
        MString line = "";
        line.append(to_string(cylinderExchangeIds[i]) + " "
                    + to_string(RANSSolver().a_isBndryGhostCell(cylinderExchangeIds[i])) + " "
                    + to_string(cylinderExchangeLocations[i * (nDim + 1) + 0]) + " "
                    + to_string(cylinderExchangeLocations[i * (nDim + 1) + 1]) + " "
                    + to_string(cylinderExchangeLocations[i * (nDim + 1) + 2]) + " "
                    + to_string(cylinderExchangeLocations[i * (nDim + 1) + 3]));
        data << line << endl;
      }
      data.close();
    }
#endif
    //--------------------------------------------------------------------------

    // determine global limiting angles of RANS sector, necessary for finding equivalent solver cells for interpolation
    MPI_Allreduce(MPI_IN_PLACE, &RANSSectors[0], 1, MPI_DOUBLE, MPI_MIN, LESSolver().mpiComm(), AT_, "MPI_IN_PLACE",
                  "RANSSectors");
    MPI_Allreduce(MPI_IN_PLACE, &RANSSectors[1], 1, MPI_DOUBLE, MPI_MAX, LESSolver().mpiComm(), AT_, "MPI_IN_PLACE",
                  "RANSSectors");
    MPI_Allreduce(MPI_IN_PLACE, &RANSSectorLimits[0], 1, MPI_DOUBLE, MPI_MIN, LESSolver().mpiComm(), AT_,
                  "MPI_IN_PLACE", "RANSSectorLimits");
    MPI_Allreduce(MPI_IN_PLACE, &RANSSectorLimits[1], 1, MPI_DOUBLE, MPI_MAX, LESSolver().mpiComm(), AT_,
                  "MPI_IN_PLACE", "RANSSectorLimits");

    m_RANSSectors[0] = RANSSectors[0];
    m_RANSSectors[1] = RANSSectors[1];
    m_RANSSectorLimits[0] = RANSSectorLimits[0];
    m_RANSSectorLimits[1] = RANSSectorLimits[1];

    m_noCylindricalGlobalExchangeLocations = 0;
    m_noCylindricalGlobalRANSExchangeValues = 0;
    m_noCylindricalGlobalLESExchangeValues = 0;
    m_noCylindricalGlobalExchangeIds = 0;
    m_noGlobalRANSExchangeCells = 0;

    // set up global arrays for exchange
    MPI_Allreduce(&noRANSExchangeCells, &m_noGlobalRANSExchangeCells, 1, MPI_INT, MPI_SUM, LESSolver().mpiComm(), AT_,
                  "noRANSExchangeCells", "m_noGlobalRANSExchangeCells");

    m_noRANSCylinderExchangeVariables = noRANSVariables() + (MInt)m_STGSponge;
    m_noCylindricalGlobalExchangeLocations = m_noGlobalRANSExchangeCells * (nDim + 1);
    m_noCylindricalGlobalRANSExchangeValues = m_noGlobalRANSExchangeCells * m_noRANSCylinderExchangeVariables;
    m_noCylindricalGlobalLESExchangeValues = m_noGlobalRANSExchangeCells * LESSolver().m_LESNoVarAverage;
    m_noCylindricalGlobalExchangeIds = m_noGlobalRANSExchangeCells;

    m_cylinderExchangeIdsOffset = -1;
    m_noRANSExchangeCells = noRANSExchangeCells;
    mAlloc(m_globalCylinderExchangeLocations, m_noCylindricalGlobalExchangeLocations,
           "m_globalCylinderExchangeLocations", F0, FUN_);
    mAlloc(m_globalCylinderRANSExchangeValues, m_noCylindricalGlobalRANSExchangeValues,
           "m_globalCylinderRANSExchangeValues", F0, FUN_);
    mAlloc(m_globalCylinderLESExchangeValues, m_noCylindricalGlobalLESExchangeValues,
           "m_globalCylinderLESExchangeValues", F0, FUN_);
    mAlloc(m_globalCylinderExchangeIds, m_noCylindricalGlobalExchangeIds, "m_globalCylinderExchangeIds", 0, FUN_);

    // Create communicator for cylindic exchange communication
    //--------------------------------------------------------------------------
    MPI_Comm_size(LESSolver().mpiComm(), &m_commSizeCylExchange);
    std::vector<MInt> cellsPerDomain(m_commSizeCylExchange);
    std::vector<MInt> cylRanks(m_commSizeCylExchange);
    MInt myCylRank = LESSolver().domainId();

    // if(m_commSizeCylExchange > 1) {
    MPI_Allgather(&noRANSExchangeCells, 1, MPI_INT, &cellsPerDomain[0], 1, MPI_INT, LESSolver().mpiComm(), AT_,
                  "noBc7909Cells ", "bc7909CellsperDomain");
    MPI_Allgather(&myCylRank, 1, MPI_INT, &cylRanks[0], 1, MPI_INT, LESSolver().mpiComm(), AT_, "myCylRank",
                  "cylRanks");

    MInt noInvolvedRanks = 0;
    std::vector<MInt> involvedRanks(m_commSizeCylExchange);
    for(MInt i = 0; i < m_commSizeCylExchange; i++) {
      if(cellsPerDomain[i]) {
        involvedRanks[noInvolvedRanks] = i;
        ++noInvolvedRanks;
      }
    }

    MPI_Comm commCyl;
    MPI_Group group;
    MPI_Group groupCyl;
    MPI_Comm_group(LESSolver().mpiComm(), &group, AT_, "group");
    MPI_Group_incl(group, noInvolvedRanks, &involvedRanks[0], &groupCyl, AT_);

    MPI_Comm_create(LESSolver().mpiComm(), groupCyl, &commCyl, AT_, "commCyl");
    m_commCyl = commCyl;

    MInt myCommRank = -1;
    for(MInt r = 0; r < noInvolvedRanks; r++) {
      if(myCylRank == involvedRanks[r]) {
        myCommRank = r;
        break;
      }
    }

    m_cylRoot = involvedRanks[0];

    // defined here for global exchange for to determine m_cylinderExchangeIdsOffset
    ScratchSpace<MInt> displsIds(noInvolvedRanks, "displsIds", FUN_);
    //--------------------------------------------------------------------------

    if(noRANSExchangeCells > 0) {
      m_cylCommActive = true;

      // exchange arrays
      MInt noLocations = noRANSExchangeCells * (nDim + 1);
      MInt noIds = noRANSExchangeCells;
      if(noInvolvedRanks > 0) {
        ScratchSpace<MInt> recvbufLocations(noInvolvedRanks, "recvbuf", FUN_);
        ScratchSpace<MInt> recvbufIds(noInvolvedRanks, "recvbuf", FUN_);
        recvbufLocations.fill(0);
        // recvbufValues.fill(0);
        recvbufIds.fill(0);

        MPI_Gather(&noLocations, 1, MPI_INT, &recvbufLocations[0], 1, MPI_INT, 0, commCyl, AT_, "noLocations",
                   "recvbufLocations");
        // MPI_Gather(&noValues, 1, MPI_INT, &recvbufValues[0], 1, MPI_INT, 0, commCyl,
        //           AT_, "noRANSExchangeCells", "recvbufValues");
        MPI_Gather(&noIds, 1, MPI_INT, &recvbufIds[0], 1, MPI_INT, 0, commCyl, AT_, "noIds", "recvbufIds");

        ScratchSpace<MInt> displsLocations(noInvolvedRanks, "displsLocations", FUN_);
        if(myCommRank == 0) {
          MInt offsetLocations = 0;
          // MInt offsetValues = 0;
          MInt offsetIds = 0;
          for(MInt dom = 0; dom < noInvolvedRanks; dom++) {
            displsLocations[dom] = offsetLocations;
            displsIds[dom] = offsetIds;
            offsetLocations += recvbufLocations[dom];
            offsetIds += recvbufIds[dom];
          }
        }

        MPI_Gatherv(&cylinderExchangeLocations[0], noLocations, MPI_DOUBLE, &m_globalCylinderExchangeLocations[0],
                    &recvbufLocations[myCommRank], &displsLocations[myCommRank], MPI_DOUBLE, 0, commCyl, AT_,
                    "cylinderExchangeLocations", "m_globalCylinderExchangeLocations");
        MPI_Gatherv(&cylinderExchangeIds[0], noIds, MPI_INT, &m_globalCylinderExchangeIds[0], &recvbufIds[myCommRank],
                    &displsIds[myCommRank], MPI_INT, 0, commCyl, AT_, "cylinderExchangeIds",
                    "m_globalCylinderExchangeIds");

      } else {
        // JANNIK:is this correct, check this
        // MInt index = 0;
        // for(MInt RANSId = 0; RANSId < a_noFvGridCellsRANS(); RANSId++) {
        //   if(!RANSSolver().a_isHalo(RANSId) && RANSSolver().c_isLeafCell(RANSId)) {
        //     MInt LESId = convertIdParent(RANSSolver(), LESSolver(), RANSId);
        //     if(LESId != -1) {
        //       for(MInt pos = 0; pos < (nDim + 1); pos++) {
        //         m_globalCylinderExchangeLocations[index * (nDim + 1) + pos] =
        //             cylinderExchangeLocations[index * (nDim + 1) + pos];
        //       }
        //       m_globalCylinderExchangeIds[index] = cylinderExchangeIds[index];
        //       index++;
        //     }
        //   }
        // }
      }
    }

    MPI_Bcast(&displsIds[0], noInvolvedRanks, MPI_INT, m_cylRoot, LESSolver().mpiComm(), AT_, "displsIds");
    for(MInt dom = 0; dom < noInvolvedRanks; dom++) {
      if(dom == myCommRank) {
        m_cylinderExchangeIdsOffset = displsIds[dom];
      }
    }

    MPI_Allreduce(MPI_IN_PLACE, &m_cylRoot, 1, MPI_INT, MPI_MAX, LESSolver().mpiComm(), AT_, "MPI_IN_PLACE",
                  "m_cylRoot");
    MPI_Bcast(&m_globalCylinderExchangeLocations[0], m_noCylindricalGlobalExchangeLocations, MPI_DOUBLE, m_cylRoot,
              LESSolver().mpiComm(), AT_, "m_globalCylinderExchangeLocations");
    MPI_Bcast(&m_globalCylinderExchangeIds[0], m_noCylindricalGlobalExchangeIds, MPI_INT, m_cylRoot,
              LESSolver().mpiComm(), AT_, "m_globalCylinderExchangeIds");

    // debug output
    //--------------------------------------------------------------------------
#ifdef _MB_DEBUG_
    if(LESSolver().domainId() == 0) {
      stringstream fn2;
      fn2.clear();
      fn2 << "globalCylinderExchange.txt";
      MString fname2 = fn2.str();
      MString header2 = "#exchangeIndex id x1 x2 x3 sector";
      ofstream data2;
      data2.precision(16);
      data2.open(fname2);
      data2 << header2 << endl;
      for(MInt exchangeIndex = 0; exchangeIndex < m_noGlobalRANSExchangeCells; exchangeIndex++) {
        MString line = "";
        line.append(to_string(exchangeIndex) + " " + to_string(m_globalCylinderExchangeIds[exchangeIndex]) + " "
                    + to_string(m_globalCylinderExchangeLocations[exchangeIndex * (nDim + 1) + 0]) + " "
                    + to_string(m_globalCylinderExchangeLocations[exchangeIndex * (nDim + 1) + 1]) + " "
                    + to_string(m_globalCylinderExchangeLocations[exchangeIndex * (nDim + 1) + 2]) + " "
                    + to_string(m_globalCylinderExchangeLocations[exchangeIndex * (nDim + 1) + 3]));
        data2 << line << endl;
      }
      data2.close();
    }
#endif
    //--------------------------------------------------------------------------

    //
    // Find equivalent LES cells to RANSExchange cells
    //
    mAlloc(m_globalCylinderInterpolationIndex, a_noFvGridCellsLES(), "m_globalCylinderInterpolationIndex", FUN_);
    mAlloc(m_cylinderInterpolationAngle, a_noFvGridCellsLES(), "m_cylinderInterpolationAngle", FUN_);
    mAlloc(m_globalCylinderInterpolationCell, m_noGlobalRANSExchangeCells, "m_globalCylinderInterpolationCell", FUN_);

    for(MInt exchangeIndex = 0; exchangeIndex < m_noCylindricalGlobalLESExchangeValues; exchangeIndex++) {
      m_globalCylinderInterpolationNumber.push_back(0);
    }

    for(MInt exchangeIndex = 0; exchangeIndex < m_noGlobalRANSExchangeCells; exchangeIndex++) {
      m_globalCylinderInterpolationCell[exchangeIndex].clear();
    }

    MInt maxInterpCells = 4;
    vector<vector<pair<MFloat, MInt>>> minDistIndex(
        a_noFvGridCellsLES(),
        vector<pair<MFloat, MInt>>(maxInterpCells, make_pair(std::numeric_limits<MFloat>::infinity(), -1)));
    vector<vector<pair<MFloat, MInt>>> minDistCell(
        m_noGlobalRANSExchangeCells,
        vector<pair<MFloat, MInt>>(maxInterpCells, make_pair(std::numeric_limits<MFloat>::infinity(), -1)));

    // find corresponding cell in exchangeValues
    for(MInt LESId = 0; LESId < a_noFvGridCellsLES(); LESId++) {
      m_globalCylinderInterpolationIndex[LESId].clear();
      MFloat x = LESSolver().a_coordinate(LESId, 0);
      MFloat y = LESSolver().a_coordinate(LESId, 1);
      MFloat z = LESSolver().a_coordinate(LESId, 2);
      MFloat cellHalfLength = F1B2 * LESSolver().c_cellLengthAtLevel(LESSolver().a_level(LESId));
      MFloat sector = maia::math::getSector(y, z, m_azimuthalAngle);
      MFloat angle = maia::math::getAngle(y, z);

      for(MInt exchangeIndex = 0; exchangeIndex < m_noGlobalRANSExchangeCells; exchangeIndex++) {
        MFloat xExchange = m_globalCylinderExchangeLocations[exchangeIndex * (nDim + 1) + 0];
        MFloat yExchange = m_globalCylinderExchangeLocations[exchangeIndex * (nDim + 1) + 1];
        MFloat zExchange = m_globalCylinderExchangeLocations[exchangeIndex * (nDim + 1) + 2];

        // fill data structure for calculation of LES sector average
        for(MInt r = 0; r < 2; r++) { // JANNIK: generalize this
          MFloat rotationAngle_ = -((sector - m_RANSSectors[r]) * m_azimuthalAngle) / 180.0 * PI;
          MFloat yRotated_ = cos(rotationAngle_) * y - sin(rotationAngle_) * z;
          MFloat zRotated_ = sin(rotationAngle_) * y + cos(rotationAngle_) * z;
          if(approx(x, xExchange, F3 * cellHalfLength) && approx(yRotated_, yExchange, F3 * cellHalfLength)
             && approx(zRotated_, zExchange, F3 * cellHalfLength)) {
            // rotated cell is in vicinity of current cell
            // save as interpolation cell
            MInt RANSId_ = m_globalCylinderExchangeIds[exchangeIndex];
            if(!RANSSolver().a_isBndryGhostCell(RANSId_) && LESSolver().c_isLeafCell(LESId)) {
              m_globalCylinderInterpolationCell[exchangeIndex].push_back(
                  make_pair<MInt, MFloat>((MInt)LESId, (MFloat)rotationAngle_));
            }
          }
        }

        // fill data structure for exchange of sector RANS data to LES domain
        MFloat rotationAngle = F0;
        MFloat rotationAngle1 = -((sector - m_RANSSectors[0]) * m_azimuthalAngle) / 180.0 * PI;
        MFloat rotationAngle2 = -((sector - m_RANSSectors[1]) * m_azimuthalAngle) / 180.0 * PI;
        MFloat minAngleDist1 = F0;
        MFloat minAngleDist2 = F0;
        if(angle + rotationAngle1 > m_RANSSectorLimits[0] || angle + rotationAngle1 < m_RANSSectorLimits[1]) {
          if(abs(angle + rotationAngle1 - m_RANSSectorLimits[0])
             < abs(angle + rotationAngle1 - m_RANSSectorLimits[1])) {
            minAngleDist1 = abs(angle + rotationAngle1 - m_RANSSectorLimits[0]);
          } else {
            minAngleDist1 = abs(angle + rotationAngle1 - m_RANSSectorLimits[1]);
          }
        }
        if(angle + rotationAngle2 > m_RANSSectorLimits[0] || angle + rotationAngle2 < m_RANSSectorLimits[1]) {
          if(abs(angle + rotationAngle2 - m_RANSSectorLimits[0])
             < abs(angle + rotationAngle2 - m_RANSSectorLimits[1])) {
            minAngleDist2 = abs(angle + rotationAngle2 - m_RANSSectorLimits[0]);
          } else {
            minAngleDist2 = abs(angle + rotationAngle2 - m_RANSSectorLimits[1]);
          }
        }
        if(minAngleDist2 > minAngleDist1) {
          rotationAngle = rotationAngle2;
        } else {
          rotationAngle = rotationAngle1;
        }
        MInt RANSId = convertId(LESSolver(), RANSSolver(), LESId);
        // equivalent RANS cell exists, no rotation necessary
        if(RANSId != -1) rotationAngle = 0;
        m_cylinderInterpolationAngle[LESId] = rotationAngle;
        MFloat yRotated_ = cos(rotationAngle) * y - sin(rotationAngle) * z;
        MFloat zRotated_ = sin(rotationAngle) * y + cos(rotationAngle) * z;
        if(approx(x, xExchange, F3 * cellHalfLength) && approx(yRotated_, yExchange, F4 * cellHalfLength)
           && approx(zRotated_, zExchange, F4 * cellHalfLength)) {
          MFloat delta_y = yRotated_ - yExchange;
          MFloat delta_z = zRotated_ - zExchange;
          MFloat dist = sqrt(POW2(delta_y) + POW2(delta_z));

          MBool checkWindowHalo = false;
          if(!checkWindowHalo) {
            if(dist < minDistIndex[LESId][maxInterpCells - 1].first) {
              minDistIndex[LESId][maxInterpCells - 1] = make_pair(dist, exchangeIndex);
              // sort by dist
              sort(minDistIndex[LESId].begin(), minDistIndex[LESId].end());
            }
          }
        }
      }
    }

    for(MInt LESId = 0; LESId < a_noFvGridCellsLES(); LESId++) {
      MFloat minY = std::numeric_limits<MFloat>::infinity();
      MFloat maxY = -std::numeric_limits<MFloat>::infinity();
      MFloat minZ = std::numeric_limits<MFloat>::infinity();
      MFloat maxZ = -std::numeric_limits<MFloat>::infinity();
      MFloat rotationAngle = m_cylinderInterpolationAngle[LESId];
      MFloat* coordLES = &LESSolver().a_coordinate(LESId, 0);
      MFloat yRotated = cos(rotationAngle) * coordLES[1] - sin(rotationAngle) * coordLES[2];
      MFloat zRotated = sin(rotationAngle) * coordLES[1] + cos(rotationAngle) * coordLES[2];
      // if(minDistIndex[LESId][0].second != -1 && minDistIndex[LESId][0].first < 1e-16){
      //   MInt exchangeIndex = minDistIndex[LESId][0].second;
      //   m_globalCylinderInterpolationIndex[LESId].push_back(exchangeIndex);
      // } else {
      for(MInt i = 0; i < maxInterpCells; i++) {
        if(minDistIndex[LESId][i].second != -1) {
          MInt exchangeIndex = minDistIndex[LESId][i].second;
          MBool usePoint = true;
          MFloat* coordRANS = &m_globalCylinderExchangeLocations[exchangeIndex * (nDim + 1)];
          if(i == maxInterpCells - 1) {
            if(yRotated > minY && yRotated < maxY && zRotated > minZ && zRotated < maxZ) usePoint = false;
          }
          if(coordRANS[1] > yRotated && coordRANS[1] > maxY) maxY = coordRANS[1];
          if(coordRANS[1] < yRotated && coordRANS[1] < minY) minY = coordRANS[1];
          if(coordRANS[2] > zRotated && coordRANS[2] > maxZ) maxZ = coordRANS[2];
          if(coordRANS[2] < zRotated && coordRANS[2] < minZ) minZ = coordRANS[2];
          if(usePoint) m_globalCylinderInterpolationIndex[LESId].push_back(exchangeIndex);
        }
      }

      // no proper 2D interpolation can be performed, use nearest neighbor extrapolation
      if(yRotated < minY || yRotated > maxY || zRotated < minZ || zRotated > maxZ) {
        if(m_globalCylinderInterpolationIndex[LESId].size() > 1) {
          m_globalCylinderInterpolationIndex[LESId].erase(m_globalCylinderInterpolationIndex[LESId].begin() + 1,
                                                          m_globalCylinderInterpolationIndex[LESId].end());
        }
      }
      // }
    }

    // delete interpolation cells if none is not a halo cell
    for(MInt exchangeIndex = 0; exchangeIndex < m_noGlobalRANSExchangeCells; exchangeIndex++) {
      MInt testForOnlyHalo = 0;
      for(MInt i = 0; i < (MInt)m_globalCylinderInterpolationCell[exchangeIndex].size(); i++) {
        MInt cellId = m_globalCylinderInterpolationCell[exchangeIndex][i].first;
        testForOnlyHalo += (MInt)LESSolver().a_isHalo(cellId);
      }

      if(testForOnlyHalo == (MInt)m_globalCylinderInterpolationCell[exchangeIndex].size()) {
        m_globalCylinderInterpolationCell[exchangeIndex].clear();
      }
    }

    // debug
    //--------------------------------------------------------------------------
#ifdef _MB_DEBUG_
    stringstream fn3;
    fn3.clear();
    fn3 << "globalCylinderLESCell_" << to_string(LESSolver().domainId()) << ".txt";
    MString fname3 = fn3.str();
    MString header3 = "#exchangeIndex ...";
    ofstream data3;
    data3.precision(16);
    data3.open(fname3);
    data3 << header3 << endl;
    for(MInt exchangeIndex = 0; exchangeIndex < m_noGlobalRANSExchangeCells; exchangeIndex++) {
      MString line = to_string(exchangeIndex) + " | ";
      for(MInt i = 0; i < (MInt)m_globalCylinderInterpolationCell[exchangeIndex].size(); i++) {
        MInt cellId = m_globalCylinderInterpolationCell[exchangeIndex][i].first;
        MInt rotationAngle = m_globalCylinderInterpolationCell[exchangeIndex][i].second;
        line.append("(" + to_string(cellId) + " " + to_string(rotationAngle) + " "
                    + to_string(LESSolver().a_isHalo(cellId)) + " ," + to_string(LESSolver().a_coordinate(cellId, 0))
                    + " " + to_string(LESSolver().a_coordinate(cellId, 1)) + " "
                    + to_string(LESSolver().a_coordinate(cellId, 2)) + ") ");
      }
      data3 << line << endl;
    }
    data3.close();

    stringstream fn4;
    fn4.clear();
    fn4 << "globalCylinderLESIndex_" << to_string(LESSolver().domainId()) << ".txt";
    MString fname4 = fn4.str();
    MString header4 = "#LESId exchangeIndex ";
    ofstream data4;
    data4.precision(16);
    data4.open(fname4);
    data4 << header4 << endl;
    for(MInt LESId = 0; LESId < a_noFvGridCellsLES(); LESId++) {
      if(m_globalCylinderInterpolationIndex[LESId].empty()) continue;
      MString line =
          to_string(LESId) + " " + to_string(LESSolver().a_isHalo(LESId)) + " "
          + to_string(LESSolver().c_isLeafCell(LESId)) + " | " + to_string(LESSolver().a_coordinate(LESId, 0)) + " "
          + to_string(LESSolver().a_coordinate(LESId, 1)) + " " + to_string(LESSolver().a_coordinate(LESId, 2)) + " ";
      for(MInt i = 0; i < (MInt)m_globalCylinderInterpolationIndex[LESId].size(); i++) {
        MInt exchangeIndex = m_globalCylinderInterpolationIndex[LESId][i];
        line.append("(" + to_string(exchangeIndex) + ", "
                    + to_string(m_globalCylinderExchangeLocations[exchangeIndex * (nDim + 1) + 0]) + " "
                    + to_string(m_globalCylinderExchangeLocations[exchangeIndex * (nDim + 1) + 1]) + " "
                    + to_string(m_globalCylinderExchangeLocations[exchangeIndex * (nDim + 1) + 2]) + " "
                    + to_string(m_globalCylinderExchangeLocations[exchangeIndex * (nDim + 1) + 3]) + ") ");
      }
      data4 << line << endl;
    }
    data4.close();
#endif
    //--------------------------------------------------------------------------
  }
}


/** \brief Interpolate RANS sector data to LES domain using mapping created in initCylinderExchange
 *    \author Jannik Borgelt
 */
template <MInt nDim, class SysEqn>
void FvZonalSTG<nDim, SysEqn>::cylinderExchange() {
  TRACE();

  if(LESSolver().grid().isActive()) {
    // interpolation
    for(MInt LESId = 0; LESId < a_noFvGridCellsLES(); LESId++) {
      // write value to cell
      if(LESSolver().a_isHalo(LESId)) continue;
      MFloat y = LESSolver().a_coordinate(LESId, 1);
      MFloat z = LESSolver().a_coordinate(LESId, 2);
      MFloat rotationAngle = m_cylinderInterpolationAngle[LESId];
      MFloat yRotated_ = cos(rotationAngle) * y - sin(rotationAngle) * z;
      MFloat zRotated_ = sin(rotationAngle) * y + cos(rotationAngle) * z;

      if(m_globalCylinderInterpolationIndex[LESId].empty()) continue;

      for(MInt var = 0; var < noRANSVariables(); var++) {
        ASSERT(LESId < (MInt)LESSolver().m_RANSValues[var].size(),
               "Trying to access data [" + to_string(var) + "][" + to_string(LESId) + "] in m_RANSValues with length "
                   + to_string(LESSolver().m_RANSValues[var].size())
                   + ", domainId: " + to_string(LESSolver().domainId()));

        if(m_globalCylinderInterpolationIndex[LESId].size() == 1) {
          MInt interpolationIndex = m_globalCylinderInterpolationIndex[LESId][0];
          LESSolver().m_RANSValues[var][LESId] =
              m_globalCylinderRANSExchangeValues[interpolationIndex * m_noRANSCylinderExchangeVariables + var];
        }

        if(m_globalCylinderInterpolationIndex[LESId].size() > 1) {
          MInt interpolationIndex = m_globalCylinderInterpolationIndex[LESId][0];
          MFloat yExchange = m_globalCylinderExchangeLocations[interpolationIndex * (nDim + 1) + 1];
          MFloat zExchange = m_globalCylinderExchangeLocations[interpolationIndex * (nDim + 1) + 2];
          MFloat delta_y = yRotated_ - yExchange;
          MFloat delta_z = zRotated_ - zExchange;

          // no interpolation necessary, since equivalent RANS cell exists
          if(sqrt(pow(delta_y, 2) + pow(delta_z, 2)) < 1e-16 /*halfCellLength/20.0*/) {
            LESSolver().m_RANSValues[var][LESId] =
                m_globalCylinderRANSExchangeValues[interpolationIndex * m_noRANSCylinderExchangeVariables + var];

            // Interpolation - Least square method
          } else {
            MFloat n = F0;
            MFloat Eyi = F0;
            MFloat Eyi2 = F0;
            MFloat Ezi = F0;
            MFloat Ezi2 = F0;
            MFloat Eyizi = F0;
            MFloat Eui = F0;
            MFloat Euiyi = F0;
            MFloat Euizi = F0;

            for(MInt i = 0; i < (MInt)m_globalCylinderInterpolationIndex[LESId].size(); i++) {
              interpolationIndex = m_globalCylinderInterpolationIndex[LESId][i];
              yExchange = m_globalCylinderExchangeLocations[interpolationIndex * (nDim + 1) + 1];
              zExchange = m_globalCylinderExchangeLocations[interpolationIndex * (nDim + 1) + 2];
              delta_y = yRotated_ - yExchange;
              delta_z = zRotated_ - zExchange;
              n += 1;
              Eyi += delta_y;
              Eyi2 += POW2(delta_y);
              Ezi += delta_z;
              Ezi2 += POW2(delta_z);
              Eyizi += delta_y * delta_z;
              Eui += m_globalCylinderRANSExchangeValues[interpolationIndex * m_noRANSCylinderExchangeVariables + var];
              Euiyi += m_globalCylinderRANSExchangeValues[interpolationIndex * m_noRANSCylinderExchangeVariables + var]
                       * delta_y;
              Euizi += m_globalCylinderRANSExchangeValues[interpolationIndex * m_noRANSCylinderExchangeVariables + var]
                       * delta_z;
            }

            MFloat detA =
                n * Eyi2 * Ezi2 + 2 * Eyi * Eyizi * Ezi - (Ezi * Eyi2 * Ezi + n * Eyizi * Eyizi + Ezi2 * Eyi * Eyi);

            if(detA > 1e-16) {
              LESSolver().m_RANSValues[var][LESId] =
                  1 / detA
                  * (Eui * (Eyi2 * Ezi2 - Eyizi * Eyizi) + Euiyi * (Ezi * Eyizi - Eyi * Ezi2)
                     + Euizi * (Eyi * Eyizi - Eyi2 * Ezi));
            }
          }
        }
      }

      // rotate veloctiy components
      MInt vId = LESSolver().sysEqn().PV->V;
      MInt wId = LESSolver().sysEqn().PV->W;
      MFloat v_ = LESSolver().m_RANSValues[vId][LESId];
      MFloat w_ = LESSolver().m_RANSValues[wId][LESId];
      LESSolver().m_RANSValues[vId][LESId] = v_ * cos(-rotationAngle) - w_ * sin(-rotationAngle);
      LESSolver().m_RANSValues[wId][LESId] = v_ * sin(-rotationAngle) + w_ * cos(-rotationAngle);
    }
  }
}


/** \brief fill m_globalCylinderRANSExchangeValues with RANS values for exchange with LES Solver
 *    \author Jannik Borgelt
 */
template <MInt nDim, class SysEqn>
void FvZonalSTG<nDim, SysEqn>::calcRANSSectorValues() {
  TRACE();

  if(LESSolver().grid().isActive()) {
    ScratchSpace<MInt> numberCylinderExchangeValues(m_noCylindricalGlobalRANSExchangeValues,
                                                    "numberCylinderExchangeValues", FUN_);
    numberCylinderExchangeValues.fill(0);

    if(m_cylCommActive) {
      // reset vector
      for(MInt exchangeIndex = 0; exchangeIndex < m_noGlobalRANSExchangeCells; exchangeIndex++) {
        for(MInt var = 0; var < m_noRANSCylinderExchangeVariables; var++) {
          m_globalCylinderRANSExchangeValues[exchangeIndex * m_noRANSCylinderExchangeVariables + var] = F0;
        }
      }

      // fill vector
      for(MInt exchangeIndex_ = 0; exchangeIndex_ < m_noRANSExchangeCells; exchangeIndex_++) {
        MInt exchangeIndex = exchangeIndex_ + m_cylinderExchangeIdsOffset;
        MInt RANSId = m_globalCylinderExchangeIds[exchangeIndex];
        if(RANSId > -1) {
          for(MInt var = 0; var < noRANSVariables(); var++) {
            m_globalCylinderRANSExchangeValues[exchangeIndex * m_noRANSCylinderExchangeVariables + var] =
                RANSSolver().a_pvariable(RANSId, var);
            numberCylinderExchangeValues[exchangeIndex * m_noRANSCylinderExchangeVariables + var] = 1;
          }
        }
      }

      MPI_Allreduce(MPI_IN_PLACE, &m_globalCylinderRANSExchangeValues[0], m_noCylindricalGlobalRANSExchangeValues,
                    MPI_DOUBLE, MPI_SUM, m_commCyl, AT_, "MPI_IN_PLACE", "m_globalCylinderRANSExchangeValues");
      MPI_Allreduce(MPI_IN_PLACE, &numberCylinderExchangeValues[0], m_noCylindricalGlobalRANSExchangeValues, MPI_INT,
                    MPI_SUM, m_commCyl, AT_, "MPI_IN_PLACE", "numberCylinderExchangeValues");
    }

    MPI_Bcast(&m_globalCylinderRANSExchangeValues[0], m_noCylindricalGlobalRANSExchangeValues, MPI_DOUBLE, m_cylRoot,
              LESSolver().mpiComm(), AT_, "m_globalCylinderRANSExchangeValues");
    MPI_Bcast(&numberCylinderExchangeValues[0], m_noCylindricalGlobalRANSExchangeValues, MPI_INT, m_cylRoot,
              LESSolver().mpiComm(), AT_, "numberCylinderExchangeValues");


    for(MInt exchangeIndex = 0; exchangeIndex < m_noGlobalRANSExchangeCells; exchangeIndex++) {
      for(MInt var = 0; var < noRANSVariables(); var++) {
        if(numberCylinderExchangeValues[exchangeIndex * m_noRANSCylinderExchangeVariables + var] > 0) {
          m_globalCylinderRANSExchangeValues[exchangeIndex * m_noRANSCylinderExchangeVariables + var] /=
              numberCylinderExchangeValues[exchangeIndex * m_noRANSCylinderExchangeVariables + var];
        }
      }
    }
  }
}


/** \brief fill m_globalCylinderLESExchangeValues with sector averaged LES values for exchange with RANS Solver
 *    \author Jannik Borgelt
 */
template <MInt nDim, class SysEqn>
void FvZonalSTG<nDim, SysEqn>::calcLESSectorAverage() {
  TRACE();

  if(LESSolver().grid().isActive()) {
    MInt vId = LESSolver().sysEqn().PV->V;
    MInt wId = LESSolver().sysEqn().PV->W;

    for(MInt exchangeIndex = 0; exchangeIndex < m_noGlobalRANSExchangeCells; exchangeIndex++) {
      if((MInt)m_globalCylinderInterpolationCell[exchangeIndex].size() == 0) continue;

      // No interpolation necessary
      if((MInt)m_globalCylinderInterpolationCell[exchangeIndex].size() == 1) {
        MInt cellId = m_globalCylinderInterpolationCell[exchangeIndex][0].first;

        MFloat rotationAngle = m_globalCylinderInterpolationCell[exchangeIndex][0].second;

        vector<MInt>::iterator findAvgId =
            find(LESSolver().m_LESAverageCells.begin(), LESSolver().m_LESAverageCells.end(), cellId);
        if(findAvgId != LESSolver().m_LESAverageCells.end()) {
          MInt avgId = distance(LESSolver().m_LESAverageCells.begin(), findAvgId);
          for(MInt var = 0; var < LESSolver().m_LESNoVarAverage; var++) {
            m_globalCylinderLESExchangeValues[exchangeIndex * LESSolver().m_LESNoVarAverage + var] =
                LESSolver().m_LESVarAverage[var][avgId];
            m_globalCylinderInterpolationNumber[exchangeIndex * LESSolver().m_LESNoVarAverage + var] = 1;
          }

          // rotate veloctiy components
          MFloat v_ = LESSolver().m_LESVarAverage[vId][avgId];
          MFloat w_ = LESSolver().m_LESVarAverage[wId][avgId];
          MFloat v_rot = v_ * cos(rotationAngle) - w_ * sin(rotationAngle);
          MFloat w_rot = v_ * sin(rotationAngle) + w_ * cos(rotationAngle);
          m_globalCylinderLESExchangeValues[exchangeIndex * LESSolver().m_LESNoVarAverage + vId] = v_rot;
          m_globalCylinderLESExchangeValues[exchangeIndex * LESSolver().m_LESNoVarAverage + wId] = w_rot;
        }
        // continue;
      }

      // Interpolation - Least square method
      if(m_globalCylinderInterpolationCell[exchangeIndex].size() > 1) {
        for(MInt var = 0; var < LESSolver().m_LESNoVarAverage; var++) {
          m_globalCylinderLESExchangeValues[exchangeIndex * LESSolver().m_LESNoVarAverage + var] = F0;

          MFloat n = F0;
          MFloat Eyi = F0;
          MFloat Eyi2 = F0;
          MFloat Ezi = F0;
          MFloat Ezi2 = F0;
          MFloat Eyizi = F0;
          MFloat Eui = F0;
          MFloat Euiyi = F0;
          MFloat Euizi = F0;

          for(MInt i = 0; i < (MInt)m_globalCylinderInterpolationCell[exchangeIndex].size(); i++) {
            MInt cellId = m_globalCylinderInterpolationCell[exchangeIndex][i].first;
            MFloat y = LESSolver().a_coordinate(cellId, 1);
            MFloat z = LESSolver().a_coordinate(cellId, 2);
            MFloat rotationAngle = m_globalCylinderInterpolationCell[exchangeIndex][i].second;
            MFloat yRotated_ = cos(rotationAngle) * y - sin(rotationAngle) * z;
            MFloat zRotated_ = sin(rotationAngle) * y + cos(rotationAngle) * z;
            MFloat yExchange = m_globalCylinderExchangeLocations[exchangeIndex * (nDim + 1) + 1];
            MFloat zExchange = m_globalCylinderExchangeLocations[exchangeIndex * (nDim + 1) + 2];
            MFloat delta_y = yExchange - yRotated_;
            MFloat delta_z = zExchange - zRotated_;

            vector<MInt>::iterator findAvgId =
                find(LESSolver().m_LESAverageCells.begin(), LESSolver().m_LESAverageCells.end(), cellId);
            if(findAvgId != LESSolver().m_LESAverageCells.end()) {
              MInt avgId = distance(LESSolver().m_LESAverageCells.begin(), findAvgId);

              // rotate veloctiy components
              MFloat v_ = LESSolver().m_LESVarAverage[vId][avgId];
              MFloat w_ = LESSolver().m_LESVarAverage[wId][avgId];
              MFloat v_rot = v_ * cos(rotationAngle) - w_ * sin(rotationAngle);
              MFloat w_rot = v_ * sin(rotationAngle) + w_ * cos(rotationAngle);

              n += 1;
              Eyi += delta_y;
              Eyi2 += POW2(delta_y);
              Ezi += delta_z;
              Ezi2 += POW2(delta_z);
              Eyizi += delta_y * delta_z;

              if(var == vId) {
                Eui += v_rot;
                Euiyi += v_rot * delta_y;
                Euizi += v_rot * delta_z;
              } else if(var == wId) {
                Eui += w_rot;
                Euiyi += w_rot * delta_y;
                Euizi += w_rot * delta_z;
              } else {
                Eui += LESSolver().m_LESVarAverage[var][avgId];
                Euiyi += LESSolver().m_LESVarAverage[var][avgId] * delta_y;
                Euizi += LESSolver().m_LESVarAverage[var][avgId] * delta_z;
              }
            }
          }

          MFloat detA =
              n * Eyi2 * Ezi2 + 2 * Eyi * Eyizi * Ezi - (Ezi * Eyi2 * Ezi + n * Eyizi * Eyizi + Ezi2 * Eyi * Eyi);

          if(detA >= 1e-12) {
            m_globalCylinderLESExchangeValues[exchangeIndex * LESSolver().m_LESNoVarAverage + var] =
                F1 / detA
                * (Eui * (Eyi2 * Ezi2 - Eyizi * Eyizi) + Euiyi * (Ezi * Eyizi - Eyi * Ezi2)
                   + Euizi * (Eyi * Eyizi - Eyi2 * Ezi));
          }

          if(!approx(m_globalCylinderLESExchangeValues[exchangeIndex * LESSolver().m_LESNoVarAverage + var], F0,
                     1e-12)) {
            m_globalCylinderInterpolationNumber[exchangeIndex * LESSolver().m_LESNoVarAverage + var] = 1;
          }
        }
      }
    }

    MPI_Allreduce(MPI_IN_PLACE, &m_globalCylinderInterpolationNumber[0], m_noCylindricalGlobalLESExchangeValues,
                  MPI_INT, MPI_SUM, LESSolver().mpiComm(), AT_, "MPI_IN_PLACE", "m_globalCylinderInterpolationNumber");

    MPI_Allreduce(MPI_IN_PLACE, &m_globalCylinderLESExchangeValues[0], m_noCylindricalGlobalLESExchangeValues,
                  MPI_DOUBLE, MPI_SUM, LESSolver().mpiComm(), AT_, "MPI_IN_PLACE", "m_globalCylinderLESExchangeValues");


    for(MInt exchangeIndex = 0; exchangeIndex < m_noGlobalRANSExchangeCells; exchangeIndex++) {
      for(MInt var = 0; var < LESSolver().m_LESNoVarAverage; var++) {
        if(m_globalCylinderInterpolationNumber[exchangeIndex * LESSolver().m_LESNoVarAverage + var] > 0) {
          m_globalCylinderLESExchangeValues[exchangeIndex * LESSolver().m_LESNoVarAverage + var] /=
              m_globalCylinderInterpolationNumber[exchangeIndex * LESSolver().m_LESNoVarAverage + var];
        }
      }
    }
  }
}


template class FvZonalSTG<2, FvSysEqnRANS<2, RANSModelConstants<RANS_SA_DV>>>;
template class FvZonalSTG<2, FvSysEqnRANS<2, RANSModelConstants<RANS_FS>>>;
template class FvZonalSTG<3, FvSysEqnRANS<3, RANSModelConstants<RANS_SA_DV>>>;
template class FvZonalSTG<3, FvSysEqnRANS<3, RANSModelConstants<RANS_FS>>>;
