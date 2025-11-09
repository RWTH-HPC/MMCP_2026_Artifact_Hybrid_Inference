// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "couplerlblb.h"
#include "coupling.h"

template <MInt nDim, MInt nDist, class SysEqn>
CouplerLbLb<nDim, nDist, SysEqn>::CouplerLbLb(const MInt couplingId, std::vector<LbSolver*> solvers)
  : Coupling(couplingId), CouplingLB<nDim, nDist, SysEqn>(couplingId, solvers) {
  initData();
  readProperties();
}

template <MInt nDim, MInt nDist, class SysEqn>
CouplerLbLb<nDim, nDist, SysEqn>::~CouplerLbLb() {
  mDeallocate(m_sourceBaseAddresses);
  mDeallocate(m_destBaseAddresses);
  mDeallocate(m_dataBlockSizes);
}

template <MInt nDim, MInt nDist, class SysEqn>
void CouplerLbLb<nDim, nDist, SysEqn>::initData() {
  m_sourceSolverId = sourceSolver().solverId();
  m_destSolverId = destSolver().solverId();
}

template <MInt nDim, MInt nDist, class SysEqn>
void CouplerLbLb<nDim, nDist, SysEqn>::init() {
  // check if solvers are compatible
  if(sourceSolver().m_updateMacroscopicLocation != destSolver().m_updateMacroscopicLocation) {
    mTerm(1, AT_, "updateAfterPropagation has to be consistent between solvers");
  }
  if(sourceSolver().m_isThermal != destSolver().m_isThermal) {
    mTerm(1, AT_, "isThermal has to be consistent between solvers");
  }

  // make a list of cells that are located within the transfer box. Does not work with adaptation!
  for(MInt sourceId = 0; sourceId < sourceSolver().a_noCells(); sourceId++) {
    if(source2DestId(sourceId) < 0) continue;
    MBool isInside = true;
    for(MInt i = 0; i < nDim; i++) {
      MFloat coordinate = sourceSolver().a_coordinate(sourceId, i);
      if(coordinate < m_transferBox[i] || coordinate > m_transferBox[i + nDim]) {
        isInside = false;
        break;
      }
    }
    if(isInside) {
      m_transferCellIds.push_back(sourceId);
    }
  }

  if(sourceSolver().m_updateMacroscopicLocation == POSTPROPAGATION) {
    if(sourceSolver().m_isThermal) {
      m_noVarsTransfer = 4;
      mAlloc(m_dataBlockSizes, m_noVarsTransfer, "m_dataBlockSizes", 0, AT_);
      mAlloc(m_sourceBaseAddresses, m_noVarsTransfer, "m_sourceBaseAddresses", AT_);
      mAlloc(m_destBaseAddresses, m_noVarsTransfer, "m_destBaseAddresses", AT_);

      m_dataBlockSizes[0] = nDim + 1;
      m_dataBlockSizes[1] = nDim + 1;
      m_dataBlockSizes[2] = nDist;
      m_dataBlockSizes[3] = nDist;

      m_sourceBaseAddresses[0] = &sourceSolver().a_variable(0, 0);
      m_sourceBaseAddresses[1] = &sourceSolver().a_oldVariable(0, 0);
      m_sourceBaseAddresses[2] = &sourceSolver().a_oldDistribution(0, 0);
      m_sourceBaseAddresses[3] = &sourceSolver().a_oldDistributionThermal(0, 0);
      m_destBaseAddresses[0] = &destSolver().a_variable(0, 0);
      m_destBaseAddresses[1] = &destSolver().a_oldVariable(0, 0);
      m_destBaseAddresses[2] = &destSolver().a_oldDistribution(0, 0);
      m_destBaseAddresses[3] = &destSolver().a_oldDistributionThermal(0, 0);
    } else {
      m_noVarsTransfer = 3;
      mAlloc(m_dataBlockSizes, m_noVarsTransfer, "m_noDataBlockSizes", 0, AT_);
      mAlloc(m_sourceBaseAddresses, m_noVarsTransfer, "m_sourceBaseAddresses", AT_);
      mAlloc(m_destBaseAddresses, m_noVarsTransfer, "m_destBaseAddresses", AT_);
      m_dataBlockSizes[0] = nDim + 1;
      m_dataBlockSizes[1] = nDim + 1;
      m_dataBlockSizes[2] = nDist;

      m_sourceBaseAddresses[0] = &sourceSolver().a_variable(0, 0);
      m_sourceBaseAddresses[1] = &sourceSolver().a_oldVariable(0, 0);
      m_sourceBaseAddresses[2] = &sourceSolver().a_oldDistribution(0, 0);
      m_destBaseAddresses[0] = &destSolver().a_variable(0, 0);
      m_destBaseAddresses[1] = &destSolver().a_oldVariable(0, 0);
      m_destBaseAddresses[2] = &destSolver().a_oldDistribution(0, 0);
    }
  } else {
    if(sourceSolver().m_isThermal) {
      m_noVarsTransfer = 3;
      mAlloc(m_dataBlockSizes, m_noVarsTransfer, "m_noDataBlockSizes", 0, AT_);
      mAlloc(m_sourceBaseAddresses, m_noVarsTransfer, "m_sourceBaseAddresses", AT_);
      mAlloc(m_destBaseAddresses, m_noVarsTransfer, "m_destBaseAddresses", AT_);

      m_dataBlockSizes[0] = nDim + 1;
      m_dataBlockSizes[1] = nDist;
      m_dataBlockSizes[2] = nDist;

      m_sourceBaseAddresses[0] = &sourceSolver().a_variable(0, 0);
      m_sourceBaseAddresses[1] = &sourceSolver().a_oldDistribution(0, 0);
      m_sourceBaseAddresses[2] = &sourceSolver().a_oldDistributionThermal(0, 0);
      m_destBaseAddresses[0] = &destSolver().a_variable(0, 0);
      m_destBaseAddresses[1] = &destSolver().a_oldDistribution(0, 0);
      m_destBaseAddresses[2] = &destSolver().a_oldDistributionThermal(0, 0);
    } else {
      m_noVarsTransfer = 2;
      mAlloc(m_dataBlockSizes, m_noVarsTransfer, "m_noDataBlockSizes", 0, AT_);
      mAlloc(m_sourceBaseAddresses, m_noVarsTransfer, "m_sourceBaseAddresses", AT_);
      mAlloc(m_destBaseAddresses, m_noVarsTransfer, "m_destBaseAddresses", AT_);

      m_dataBlockSizes[0] = nDim + 1;
      m_dataBlockSizes[1] = nDist;

      m_sourceBaseAddresses[0] = &sourceSolver().a_variable(0, 0);
      m_sourceBaseAddresses[1] = &sourceSolver().a_oldDistribution(0, 0);
      m_destBaseAddresses[0] = &destSolver().a_variable(0, 0);
      m_destBaseAddresses[1] = &destSolver().a_oldDistribution(0, 0);
    }
  }

#ifndef NDEBUG
  // debug checks
  for(const auto& sourceId : m_transferCellIds) {
    const MInt destId = source2DestId(sourceId);
    MBool isCorrect = true;
    for(MInt i = 0; i < nDim; i++) {
      if(fabs(sourceSolver().a_coordinate(sourceId, i) - destSolver().a_coordinate(destId, i)) > 1e-10) {
        isCorrect = false;
        break;
      }
      if(!isCorrect) {
        std::cerr << "source:" << sourceId << " " << sourceSolver().a_coordinate(sourceId, 0) << " "
                  << sourceSolver().a_coordinate(sourceId, 1) << " " << sourceSolver().a_coordinate(sourceId, 2) << " "
                  << "dest:" << destId << " " << destSolver().a_coordinate(destId, 0) << " "
                  << destSolver().a_coordinate(destId, 1) << " " << destSolver().a_coordinate(destId, 2) << std::endl;
      }
    }
  }
  for(MInt destId = 0; destId < destSolver().a_noCells(); destId++) {
    MBool isInside = true;
    for(MInt i = 0; i < nDim; i++) {
      MFloat coordinate = destSolver().a_coordinate(destId, i);
      if(coordinate < m_transferBox[i] || coordinate > m_transferBox[i + nDim]) {
        isInside = false;
        break;
      }
    }
    if(isInside) {
      const MInt sourceId = dest2SourceId(destId);
      if(std::find(m_transferCellIds.begin(), m_transferCellIds.end(), sourceId) == m_transferCellIds.end()) {
        std::cerr << "missed cell dest:" << destId << " " << destSolver().a_coordinate(destId, 0) << " "
                  << destSolver().a_coordinate(destId, 1) << " " << destSolver().a_coordinate(destId, 2) << " "
                  << "source:" << sourceId << " " << sourceSolver().a_coordinate(sourceId, 0) << " "
                  << sourceSolver().a_coordinate(sourceId, 1) << " " << sourceSolver().a_coordinate(sourceId, 2)
                  << std::endl;
      }
    }
  }
#endif
}

template <MInt nDim, MInt nDist, class SysEqn>
void CouplerLbLb<nDim, nDist, SysEqn>::readProperties() {
  if(Context::propertyExists("lbBox", 0)) {
    for(MInt i = 0; i < nDim * 2; i++) {
      /*! \page propertyPage1
        \section lbBox
        <code>MFloat CouplerLbLb::m_transferBox </code>\n
        This box defines the volume in which the cell values are transfered from one LB solver to the other.\n
        Keywords: <i>LbLb</i>
      */
      m_transferBox[i] = Context::getSolverProperty<MFloat>("lbBox", 0, AT_, i);
    }
  } else {
    mTerm(1, AT_, "Property lbBox required for coupler LBLb!");
  }
}

template <MInt nDim, MInt nDist, class SysEqn>
void CouplerLbLb<nDim, nDist, SysEqn>::postCouple(MInt) {
  for(const auto& sourceId : m_transferCellIds) {
    const MInt destId = source2DestId(sourceId);
    for(MInt var = 0; var < m_noVarsTransfer; var++) {
      memcpy(m_destBaseAddresses[var] + (m_dataBlockSizes[var] * destId),
             m_sourceBaseAddresses[var] + (m_dataBlockSizes[var] * sourceId),
             m_dataBlockSizes[var] * sizeof(MFloat));
    }
  }
}

template class CouplerLbLb<2, 9, maia::lb::LbSysEqnIncompressible<2, 9>>;
template class CouplerLbLb<3, 19, maia::lb::LbSysEqnIncompressible<3, 19>>;
template class CouplerLbLb<3, 27, maia::lb::LbSysEqnIncompressible<3, 27>>;
