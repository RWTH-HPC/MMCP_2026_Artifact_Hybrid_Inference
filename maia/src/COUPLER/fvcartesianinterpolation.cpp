// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "fvcartesianinterpolation.h"

#include <algorithm>
#include <stack>
#include "FV/fvcartesiansolverxd.h"
#include "GRID/cartesiannetcdf.h"
#include "IO/parallelio.h"
#include "MEMORY/alloc.h"
#include "UTIL/functions.h"
#include "coupling.h"
#include "globals.h"


using namespace std;


template <MInt nDim, class SysEqnOld, class SysEqnNew>
FvCartesianInterpolation<nDim, SysEqnOld, SysEqnNew>::FvCartesianInterpolation(const MInt couplingId,
                                                                               OldFvSolver* oldS,
                                                                               NewFvSolver* newS)
  : Coupling(couplingId), m_oldSolver(oldS), m_newSolver(newS) {
  TRACE();

  m_newSolverId = newSolver().m_solverId;
  m_oldSolverId = oldSolver().m_solverId;

  /*! \page propertyPage1
    \section zonal
    <code>MBool FvZonal::m_nonZonalRestart</code>\n
    default = <code>0</code>\n
    Triggers loading of restartFile of nonZonalRestartSolver
      Keywords: <i>FINITE_VOLUME, FV_ZONAL_STG</i>
  */
  m_nonZonalRestart = false;
  if(Context::propertyExists("nonZonalRestart")) {
    m_nonZonalRestart = Context::getBasicProperty<MBool>("nonZonalRestart", AT_, &m_nonZonalRestart);
  }

  // Allocate all necessary arrays for zonal exchange (necessary for zonal BC)
  if(!newSolver().m_rans) {
    newSolver().m_noRANSVariables = oldSolver().noVariables();
    mAlloc(newSolver().m_RANSValues, newSolver().m_noRANSVariables, "newSolver().m_RANSValues", AT_);
    for(MInt i = 0; i < newSolver().m_noRANSVariables; i++) {
      newSolver().m_RANSValues[i].clear();
    }
  } else {
    newSolver().m_noLESVariables = oldSolver().noVariables();
    mAlloc(newSolver().m_LESValues, newSolver().m_noLESVariables, "newSolver().m_LESValues", AT_);
    for(MInt i = 0; i < newSolver().m_noLESVariables; i++) {
      newSolver().m_LESValues[i].clear();
    }
  }
}


template <MInt nDim, class SysEqnOld, class SysEqnNew>
void FvCartesianInterpolation<nDim, SysEqnOld, SysEqnNew>::init() {
  TRACE();

  if(m_nonZonalRestart) {
    transferSolverData();
  }
}

template <MInt nDim, class SysEqnOld, class SysEqnNew>
void FvCartesianInterpolation<nDim, SysEqnOld, SysEqnNew>::preCouple(MInt) {
  if(m_nonZonalRestart) {
    for(MInt oldId = 0; oldId < a_noFvGridCellsOld(); oldId++) {
      oldSolver().reduceData(oldId, &oldSolver().a_pvariable(0, 0), oldSolver().noVariables());
    }
    transferSolverData();
  }
}

template <MInt nDim, class SysEqnOld, class SysEqnNew>
void FvCartesianInterpolation<nDim, SysEqnOld, SysEqnNew>::postCouple(MInt) {
  if(m_nonZonalRestart) {
    for(MInt oldId = 0; oldId < a_noFvGridCellsOld(); oldId++) {
      oldSolver().reduceData(oldId, &oldSolver().a_pvariable(0, 0), oldSolver().noVariables());
    }
    transferSolverData();
  }
}

/** \brief interpolate variables from old to new
 *    \author Jannik Borgelt
 */
template <MInt nDim, class SysEqnOld, class SysEqnNew>
void FvCartesianInterpolation<nDim, SysEqnOld, SysEqnNew>::transferSolverData() {
  TRACE();

  m_log << "starting interpolate (FV)" << endl;

  if(newSolver().grid().isActive()) {
    // start from min level and go to max level
    for(MInt newId = 0; newId < a_noFvGridCellsNew(); newId++) {
      MInt oldId = convertIdParent(newSolver(), oldSolver(), newId);

      if(oldId != -1) {
        for(MInt varId = 0; varId < newSolver().noVariables(); varId++) {
          newSolver().a_pvariable(newId, varId) = oldSolver().a_pvariable(oldId, varId);
        }
        IF_CONSTEXPR(SysEqnNew::m_noRansEquations == 0) {
          for(MInt s = 0; s < newSolver().m_noSpecies; ++s) {
            newSolver().a_pvariable(newId, newSolver().PV->Y[s]) =
                (oldSolver().m_noSpecies > 0) ? oldSolver().a_pvariable(oldId, oldSolver().PV->Y[s]) : F0;
          }
        }
      }
    }

    newSolver().computeConservativeVariables();
    newSolver().exchangeAll();

    m_log << "finished interpolate (FV)" << endl;
  }
}


template class FvCartesianInterpolation<2, FvSysEqnRANS<2, RANSModelConstants<RANS_SA_DV>>, FvSysEqnNS<2>>;
template class FvCartesianInterpolation<3, FvSysEqnRANS<3, RANSModelConstants<RANS_SA_DV>>, FvSysEqnNS<3>>;
template class FvCartesianInterpolation<2,
                                        FvSysEqnRANS<2, RANSModelConstants<RANS_SA_DV>>,
                                        FvSysEqnRANS<2, RANSModelConstants<RANS_SA_DV>>>;
template class FvCartesianInterpolation<3,
                                        FvSysEqnRANS<3, RANSModelConstants<RANS_SA_DV>>,
                                        FvSysEqnRANS<3, RANSModelConstants<RANS_SA_DV>>>;
