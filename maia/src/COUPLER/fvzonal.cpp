// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "fvzonal.h"

#include <algorithm>
#include <iostream>
#include <stack>
#include <vector>
#include "COUPLER/coupling.h"
#include "FV/fvcartesiancellcollector.h"
#include "MEMORY/alloc.h"
#include "UTIL/functions.h"
#include "UTIL/kdtree.h"
#include "globals.h"
#include "globalvariables.h"

using namespace std;


template <MInt nDim, class SysEqn>
FvZonal<nDim, SysEqn>::FvZonal(const MInt couplingId, RANS* R, LES* L)
  : Coupling(couplingId), CouplingFv<nDim, SysEqn>(couplingId, R), CouplingFv<nDim, FvSysEqnNS<nDim>>(couplingId, L) {
  /*! \page propertyPage1
    \section zonal
    <code>MInt FvZonal::m_zonalAveragingTimeStep</code>\n
    default = <code>NONE</code>\n \n
    Time step at which the zonal averaging starts.
    Keywords: <i>FINITE_VOLUME, FV_ZONAL_STG</i>
  */
  m_zonalAveragingTimeStep = Context::getBasicProperty<MInt>("zonalAveragingTimeStep", AT_, &m_zonalAveragingTimeStep);

  /*! \page propertyPage1
    \section zonal
    <code>MInt FvZonalSTG::m_zonalTransferInterval</code>\n
    default = <code>0</code>\n \n
    Time step interval at which the zonal values are exchanged.
    Keywords: <i>FINITE_VOLUME, FV_ZONAL_STG</i>
  */
  m_zonalTransferInterval = 50;
  if(Context::propertyExists("zonalTransferInterval")) {
    m_zonalTransferInterval = Context::getBasicProperty<MInt>("zonalTransferInterval", AT_, &m_zonalTransferInterval);
  }

  /*! \page propertyPage1
    \section zonal
    <code>MInt FvZonalRTV::m_restartLESAverage</code>\n
    default = <code>False</code>\n
    Triggers loading of LES Average
    Keywords: <i>FINITE_VOLUME, FV_ZONAL</i>
  */
  m_restartLESAverage = false;
  if(Context::propertyExists("restartLESAverage")) {
    m_restartLESAverage = Context::getBasicProperty<MBool>("restartLESAverage", AT_, &m_restartLESAverage);
  }

  /*! \page propertyPage1
    \section zonal
    <code>MBool FvZonal::m_cylindricCommunication</code>\n
    default = <code>0</code>\n
    Triggers communication from sector RANS to full-360 degree LES
    Keywords: <i>FINITE_VOLUME, FV_ZONAL_STG</i>
  */
  m_cylindricCommunication = false;
  if(Context::propertyExists("cylindricCommunication")) {
    m_cylindricCommunication =
        Context::getBasicProperty<MBool>("cylindricCommunication", AT_, &m_cylindricCommunication);
    m_azimuthalAngle = Context::getBasicProperty<MFloat>("azimuthalAngle", AT_, &m_azimuthalAngle);
  }

  if(Context::propertyExists("STGSponge")) {
    m_STGSponge = Context::getBasicProperty<MBool>("STGSponge", AT_, &m_STGSponge);
  }

  mAlloc(LESSolver().m_RANSValues, noRANSVariables() + (MInt)m_STGSponge, "LESSolver().m_RANSValues", AT_);
  LESSolver().m_noRANSVariables = noRANSVariables() + (MInt)m_STGSponge;
  for(MInt i = 0; i < noRANSVariables() + (MInt)m_STGSponge; i++) {
    LESSolver().m_RANSValues[i].clear();
  }

  mAlloc(RANSSolver().m_LESValues, noRANSVariables(), "RANSSolver().m_LESValues", AT_);
  RANSSolver().m_noLESVariables = noRANSVariables();
  for(MInt i = 0; i < noRANSVariables(); i++) {
    RANSSolver().m_LESValues[i].clear();
  }
}


/** \brief Initialize RANSValues for LES Solver
 *    \author Jannik Borgelt
 */
template <MInt nDim, class SysEqn>
void FvZonal<nDim, SysEqn>::initRANSValues() {
  TRACE();

  // write RANS data to LES Solver
  if(LESSolver().grid().isActive()) {
    for(MInt cellId = 0; cellId < a_noFvGridCellsLES(); cellId++) {
      MInt RANSId = convertIdParent(LESSolver(), RANSSolver(), cellId);
      if(RANSId > -1) {
        for(MInt var = 0; var < noRANSVariables(); var++) {
          ASSERT(cellId < (MInt)LESSolver().m_RANSValues[var].size(),
                 "Trying to access data [" + to_string(var) + "][" + to_string(cellId)
                     + "] in m_RANSValues with length " + to_string(LESSolver().m_RANSValues[var].size())
                     + ", domainId: " + to_string(LESSolver().domainId()));

          LESSolver().m_RANSValues[var][cellId] = RANSSolver().a_pvariable(RANSId, var);
        }
      } else {
        for(MInt var = 0; var < noRANSVariables(); var++) {
          if(m_cylindricCommunication) {
            if(var < noLESVariables()) {
              LESSolver().m_RANSValues[var][cellId] = LESSolver().a_pvariable(cellId, var);
            } else {
              LESSolver().m_RANSValues[var][cellId] = F0;
            }
          }
        }
      }
      if(m_STGSponge) {
        LESSolver().m_RANSValues[noRANSVariables()][cellId] = F0;
      }
    }
  }
}


/** \brief Initialize LESValues for RANS Solver
 *    \author Jannik Borgelt
 */
template <MInt nDim, class SysEqn>
void FvZonal<nDim, SysEqn>::initLESValues() {
  TRACE();

  // write RANS data to LES Solver
  if(RANSSolver().grid().isActive()) {
    for(MInt var = 0; var < noRANSVariables(); var++) {
      for(MInt cellId = 0; cellId < a_noFvGridCellsRANS(); cellId++) {
        MInt LESId = convertIdParent(RANSSolver(), LESSolver(), cellId);

        ASSERT(cellId < (MInt)RANSSolver().m_LESValues[var].size(),
               "Trying to access data [" + to_string(var) + "][" + to_string(cellId) + "] in m_RANSValues with length "
                   + to_string(RANSSolver().m_LESValues[var].size())
                   + ", domainId: " + to_string(RANSSolver().domainId()));

        if(LESId != -1) {
          RANSSolver().m_LESValues[var][cellId] = LESSolver().a_pvariable(LESId, var);
          // init nu_t
          if(var == RANSSolver().sysEqn().PV->N) {
            RANSSolver().m_LESValues[var][cellId] = RANSSolver().a_pvariable(cellId, var);
          }
        }
      }
    }
  }
}


template class FvZonal<2, FvSysEqnRANS<2, RANSModelConstants<RANS_SA_DV>>>;
template class FvZonal<2, FvSysEqnRANS<2, RANSModelConstants<RANS_FS>>>;
template class FvZonal<3, FvSysEqnRANS<3, RANSModelConstants<RANS_SA_DV>>>;
template class FvZonal<3, FvSysEqnRANS<3, RANSModelConstants<RANS_FS>>>;
