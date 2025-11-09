// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "lbinterface.h"

#include "IO/context.h"
#include "IO/infoout.h"
#include "lbinterfacecell.h"
#include "lbsolver.h"
#include "property.h"

using namespace std;

/**
 * \brief Base class for the concrete interface treatment
 *
 * \author Moritz Waldmann
 *
 * \param[in] LB solver pointer
 */
template <MInt nDim>
LbInterface<nDim>::LbInterface(LbSolver<nDim>* solver)
  : m_noDistributions(solver->m_noDistributions),
    m_interfaceChildren(solver->m_interfaceChildren),
    m_interfaceParents(solver->m_interfaceParents),
    m_externalForcing(solver->m_externalForcing),
    m_methodId(solver->m_methodId),
    m_nu(solver->m_nu),
    m_isThermal(solver->m_isThermal),
    m_innerEnergy(solver->m_innerEnergy) {
  TRACE();
  m_solver = solver;
  PV = solver->PV;

  m_Fext = solver->m_Fext;
  if(m_externalForcing) {
    cerr << " Using external forcing at Interfaces!" << endl;
    m_log << " Using external forcing at Interfaces!" << endl;
  }
  m_isEELiquid = solver->m_isEELiquid;
  m_cellDependentForcing = m_isEELiquid;
  m_Fg = solver->m_EELiquid.Fg;
}

/**
 * \brief Destructor
 *
 * \author Julian Vorspohl <j.vorspohl@aia.rwth-aachen.de>
 */
template <MInt nDim>
LbInterface<nDim>::~LbInterface() {
  TRACE();
}

// A debug function to print out all interface cells to the log
template <MInt nDim>
void LbInterface<nDim>::printInterfaceCells() {
  TRACE();
  MInt noNghbrs = 0, noParentNghbrs = 0;
  for(MInt i = m_solver->minLevel(); i < m_solver->maxLevel(); i++) {
    m_log << m_interfaceChildren[i - m_solver->minLevel()]->size() << " Interface cells on level: " << i
          << " ----------" << endl;
    for(MInt j = 0; j < m_interfaceChildren[i - m_solver->minLevel()]->size(); j++) {
      for(MInt x = 0; x < m_noDistributions - 1; x++) {
        noNghbrs += m_solver->a_hasNeighbor(m_interfaceChildren[i - m_solver->minLevel()]->a[j].m_cellId, x);
        noParentNghbrs += m_solver->a_hasNeighbor(
            m_solver->c_parentId(m_interfaceChildren[i - m_solver->minLevel()]->a[j].m_cellId), x);
      }
      m_log << " cellID " << m_interfaceChildren[i - m_solver->minLevel()]->a[j].m_cellId
            << " position: " << m_interfaceChildren[i - m_solver->minLevel()]->a[j].m_position
            << " level: " << m_solver->a_level(m_interfaceChildren[i - m_solver->minLevel()]->a[j].m_cellId)
            << " neighbors: " << noNghbrs << " parentNghbrs: " << noParentNghbrs << endl;
      m_log << " interpolation Neighbors: ";
      for(MInt m = 0; m < IPOW2(nDim); m++) {
        m_log << m_interfaceChildren[i - m_solver->minLevel()]->a[j].m_interpolationNeighbors[m] << " ("
              << m_interfaceChildren[i - m_solver->minLevel()]->a[j].m_interpolationCoefficients[m] << "), ";
      }
      m_log << endl << " - - - - -" << endl;
      noNghbrs = 0;
      noParentNghbrs = 0;
    }
  }
  TERMM(0, AT_);
}

// A debug function for making the interface cells visible
template <MInt nDim>
void LbInterface<nDim>::colorInterface() {
  TRACE();
  MBool interfaceParent;
  MBool nonOverlapping;
  for(MInt parentId = 0; parentId < m_solver->grid().noCells(); parentId++) {
    interfaceParent = false;
    nonOverlapping = false;
    for(MInt k = 0; k < IPOW2(nDim); k++) {
      if(m_solver->c_childId(parentId, k) == -1) continue;
      if(m_solver->a_isInterfaceChild(m_solver->c_childId(parentId, k))) {
        nonOverlapping = true;
        break;
      }
    }
    if(nonOverlapping) continue;
    // Check if neighbor is interface parent
    for(MInt n = 0; n < m_noDistributions - 1; n++) {
      for(MInt nk = 0; nk < IPOW2(nDim); nk++) {
        if(m_solver->c_childId(parentId, nk) == -1) continue;
        if(m_solver->a_isInterfaceChild(m_solver->c_childId(m_solver->c_neighborId(parentId, n), nk))) {
          interfaceParent = true;
          break;
        }
      }
      if(interfaceParent) break;
    }
    if(!interfaceParent) continue;
    for(MInt nk = 0; nk < IPOW2(nDim); nk++) {
      if(m_solver->c_childId(parentId, nk) == -1) continue;
      m_solver->a_variable(m_solver->c_childId(parentId, nk), 1) = -1.1;
    }
  }
}

// Explicit instantiations for 2D and 3D
template class LbInterface<2>;
template class LbInterface<3>;
