// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "couplerlbfveemultiphase.h"
#include "LB/lbconstants.h"
#include "couplerlbfv.h"
#include "coupling.h"

template <MInt nDim, MInt nDist, class SysEqnLb, class SysEqnFv>
class CouplerLbFv;

template <MInt nDim, MInt nDist, class SysEqnLb, class SysEqnFv>
CouplerLbFvEEMultiphase<nDim, nDist, SysEqnLb, SysEqnFv>::CouplerLbFvEEMultiphase(const MInt couplerId, LbSolver* lb,
                                                                                  FvCartesianSolver* fv)
  : Coupling(couplerId), CouplerLbFv<nDim, nDist, SysEqnLb, SysEqnFv>(couplerId, lb, fv) {
  initData();
  readProperties();
  initAlpha();
}

template <MInt nDim, MInt nDist, class SysEqnLb, class SysEqnFv>
void CouplerLbFvEEMultiphase<nDim, nDist, SysEqnLb, SysEqnFv>::preCouple(MInt /*step*/) {
  if(m_alphaConvergenceCheck > 0) {
    lbSolver().storeOldDistributions();
    if(m_updateAfterPropagation) lbSolver().storeOldVariables();
  }
}

template <MInt nDim, MInt nDist, class SysEqnLb, class SysEqnFv>
void CouplerLbFvEEMultiphase<nDim, nDist, SysEqnLb, SysEqnFv>::init() {
  this->initConversionFactors();
  if(!fvSolver().m_EEGas.depthCorrection) {
    mAlloc(m_depthCorrectionValues, fvSolver().maxNoGridCells(), "m_depthCorrectionValues", F0, AT_);
    initDepthcorrection();
  }
}

template <MInt nDim, MInt nDist, class SysEqnLb, class SysEqnFv>
void CouplerLbFvEEMultiphase<nDim, nDist, SysEqnLb, SysEqnFv>::finalizeCouplerInit() {
  if(!fvSolver().m_EEGas.depthCorrection) {
    initDepthcorrection();
  }
  transferUFv2Lb();
}

template <MInt nDim, MInt nDist, class SysEqnLb, class SysEqnFv>
void CouplerLbFvEEMultiphase<nDim, nDist, SysEqnLb, SysEqnFv>::balancePre() {
  if(!fvSolver().m_EEGas.depthCorrection) {
    TERMM(1, "not implemented for in-coupler depthCorrection!");
  }
}

template <MInt nDim, MInt nDist, class SysEqnLb, class SysEqnFv>
void CouplerLbFvEEMultiphase<nDim, nDist, SysEqnLb, SysEqnFv>::subCouple(MInt, MInt,
                                                                         std::vector<MBool>& solverCompleted) {
  static MInt noAlphaIterations = 0;

  if(m_redistributeAlpha && !m_disableSubstepAlphaRedist
     && !(solverCompleted[a_lbSolverId()] && solverCompleted[a_fvSolverId()]))
    correctInvalidAlpha();

  if(solverCompleted[a_lbSolverId()] && solverCompleted[a_fvSolverId()]) {
    if(m_redistributeAlpha) correctInvalidAlpha();
    // both Solvers completed -> check for convergence of alpha
    MBool alphaConverged = true;
    if(noAlphaIterations < m_maxNoAlphaIterations) alphaConverged = checkAlphaConverged();
    if(alphaConverged) {
      transferAlphaFv2Lb();
      transferUFv2Lb();

#ifndef NDEBUG
      if(fvSolver().domainId() == 0 && m_alphaConvergenceCheck != 0 && noAlphaIterations > 3) {
        if(noAlphaIterations < m_maxNoAlphaIterations)
          std::cerr << "Globaltimestep: " << globalTimeStep << "   alpha Converged after " << noAlphaIterations
                    << " Iterations!" << std::endl;
        else
          std::cerr << "Globaltimestep: " << globalTimeStep
                    << "   max number of alphaIterations reached: " << noAlphaIterations << std::endl;
      }
#endif

      noAlphaIterations = 0;
      return;
    } else {
      transferAlphaFv2Lb();
      revertLb();
      revertFv();
      noAlphaIterations++;
      solverCompleted[a_lbSolverId()] = false;
      solverCompleted[a_fvSolverId()] = false;
      return;
    }
  } else if(solverCompleted[a_lbSolverId()]) {
    // LB Solver completed -> transfer new liquid variables
    transferPressureLb2Fv(fvSolver().m_RKalpha[fvSolver().m_RKStep], m_updateFVBC);
    if(fvSolver().m_RKStep == 0)
      transferULb2Fv<true>(fvSolver().m_RKalpha[fvSolver().m_RKStep]);
    else
      transferULb2Fv<false>(fvSolver().m_RKalpha[fvSolver().m_RKStep]);
    transferNuLb2Fv(fvSolver().m_RKalpha[fvSolver().m_RKStep]);
    fvSolver().exchange();
    return;
  } else if(solverCompleted[a_fvSolverId()]) {
    mTerm(1, AT_, "Fv Solver cant be completed before Lb Solver is completed in LB-FV-Euler-Euler-Multiphase!");
  }
}

template <MInt nDim, MInt nDist, class SysEqnLb, class SysEqnFv>
void CouplerLbFvEEMultiphase<nDim, nDist, SysEqnLb, SysEqnFv>::finalizeSubCoupleInit(MInt) {
  transferPressureLb2Fv(1.0, false);
}

template <MInt nDim, MInt nDist, class SysEqnLb, class SysEqnFv>
void CouplerLbFvEEMultiphase<nDim, nDist, SysEqnLb, SysEqnFv>::initData() {}

template <MInt nDim, MInt nDist, class SysEqnLb, class SysEqnFv>
void CouplerLbFvEEMultiphase<nDim, nDist, SysEqnLb, SysEqnFv>::initAlpha() {
  switch(m_initAlphaMethod) {
    case 0: // initialize alphaGas to be zero everywhere
    {
      for(MInt lbCellId = 0; lbCellId < a_noLbCells(); lbCellId++) {
        lbSolver().a_alphaGas(lbCellId) = 0.0;
      }
      for(MInt fvCellId = 0; fvCellId < a_noFvCells(); fvCellId++) {
        fvSolver().a_alphaGas(fvCellId) = 0.0;
      }
      break;
    } // case 0

    case 1: // initialize alphaGas to a constant value given in m_initialAlpha
    {
      const MFloat alpha = m_initialAlpha;
      for(MInt lbCellId = 0; lbCellId < a_noLbCells(); lbCellId++) {
        lbSolver().a_alphaGas(lbCellId) = alpha;
      }
      for(MInt fvCellId = 0; fvCellId < a_noFvCells(); fvCellId++) {
        fvSolver().a_alphaGas(fvCellId) = alpha;
      }
      break;
    } // case 1

    case 2: // initialize alphaGas to a constant value given in m_initialAlpha in a sphere with radius 0.5 around the
            // origin, else m_alphaInf
    {
      const MFloat alpha = m_initialAlpha;
      for(MInt i = 0; i < a_noLbCells(); i++) {
        MFloat xcoord = lbSolver().a_coordinate(i, 0);
        MFloat ycoord = lbSolver().a_coordinate(i, 1);
        MFloat zcoord = lbSolver().a_coordinate(i, 2);
        if(xcoord * xcoord + ycoord * ycoord + zcoord * zcoord < 0.5) {
          lbSolver().a_alphaGas(i) = alpha;
        } else {
          lbSolver().a_alphaGas(i) = m_alphaInf;
        }
      }
      for(MInt i = 0; i < a_noFvCells(); i++) {
        MFloat xcoord = fvSolver().a_coordinate(i, 0);
        MFloat ycoord = fvSolver().a_coordinate(i, 1);
        MFloat zcoord = fvSolver().a_coordinate(i, 2);
        if(xcoord * xcoord + ycoord * ycoord + zcoord * zcoord < 0.5) {
          fvSolver().a_alphaGas(i) = alpha;
        } else {
          fvSolver().a_alphaGas(i) = m_alphaInf;
        }
      }
      break;
    } // case 2

    case 3: // initialize alphaGas to a constant value given in m_initialAlpha below the z value of 0, else m_alphaInf
    {
      const MFloat alpha = m_initialAlpha;
      for(MInt i = 0; i < a_noLbCells(); i++) {
        MFloat zcoord = lbSolver().a_coordinate(i, 2);
        if(zcoord < 0.0) {
          lbSolver().a_alphaGas(i) = alpha;
        } else {
          lbSolver().a_alphaGas(i) = m_alphaInf;
        }
      }
      for(MInt i = 0; i < a_noFvCells(); i++) {
        MFloat zcoord = fvSolver().a_coordinate(i, 2);
        if(zcoord < 0.0) {
          fvSolver().a_alphaGas(i) = alpha;
        } else {
          fvSolver().a_alphaGas(i) = m_alphaInf;
        }
      }
      break;
    }

    case 4: // alpha Gradient in z direction with initial alpha at 0, alphaInf at 1 and (2 initial alpha - alphaInf) at
            // -1
    {
      const MFloat delAlpha = m_initialAlpha - m_alphaInf;
      for(MInt i = 0; i < a_noLbCells(); i++) {
        MFloat zcoord = lbSolver().a_coordinate(i, 2);
        if(zcoord < -1.0) {
          lbSolver().a_alphaGas(i) = m_initialAlpha + delAlpha;
        } else if(zcoord > 1.0) {
          lbSolver().a_alphaGas(i) = m_initialAlpha - delAlpha;
        } else {
          lbSolver().a_alphaGas(i) = m_initialAlpha - zcoord * delAlpha;
        }
      }
      for(MInt i = 0; i < a_noFvCells(); i++) {
        MFloat zcoord = fvSolver().a_coordinate(i, 2);
        if(zcoord < -1.0) {
          fvSolver().a_alphaGas(i) = m_initialAlpha + delAlpha;
        } else if(zcoord > 1.0) {
          fvSolver().a_alphaGas(i) = m_initialAlpha - delAlpha;
        } else {
          fvSolver().a_alphaGas(i) = m_initialAlpha - zcoord * delAlpha;
        }
      }
      break;
    }

    default: {
      mTerm(1, AT_, "initAlphaMethod not matching any case!");
      break;
    } // default

  } // switch (m_initAlphaMethod)
}

template <MInt nDim, MInt nDist, class SysEqnLb, class SysEqnFv>
void CouplerLbFvEEMultiphase<nDim, nDist, SysEqnLb, SysEqnFv>::readProperties() {
  /*! \page propertyPage1
    \section alphaConvergenceCheck
    <code>MInt CouplerLbFvEEMultiphase::m_alphaConvergenceCheck </code>\n
    default = <code>0</code>\n \n
    Method to check for alpha convergence.\n
    <li>0 : don't check for convergence</li>
    <li>1 : check if max diff between alphaOld and alphaNew is smaller the epsAlpha</li>
    Keywords: <i>EEMultiphase</i>
  */
  m_alphaConvergenceCheck = 0;
  if(Context::propertyExists("alphaConvergenceCheck", 0)) {
    m_alphaConvergenceCheck = Context::getBasicProperty<MInt>("alphaConvergenceCheck", AT_, &m_alphaConvergenceCheck);
  }

  /*! \page propertyPage1
    \section maxNoAlphaIterations
    <code>MInt CouplerLbFvEEMultiphase::m_maxNoAlphaIterations </code>\n
    default = <code>4</code>\n \n
    Maximum number of iterations for alpha convergence.\n
    Keywords: <i>EEMultiphase</i>
  */
  m_maxNoAlphaIterations = 4;
  if(Context::propertyExists("maxNoAlphaIterations", 0)) {
    m_maxNoAlphaIterations = Context::getBasicProperty<MInt>("maxNoAlphaIterations", AT_, &m_maxNoAlphaIterations);
  }

  /*! \page propertyPage1
    \section epsAlpha
    <code>MFloat CouplerLbFvEEMultiphase::m_epsAlpha </code>\n
    default = <code>1.0e-6</code>\n \n
    The eps for the alpha convergence check.\n
    Keywords: <i>EEMultiphase</i>
  */
  m_epsAlpha = 1.0e-6;
  if(Context::propertyExists("epsAlpha", 0)) {
    m_epsAlpha = Context::getBasicProperty<MFloat>("epsAlpha", AT_, &m_epsAlpha);
  }

  /*! \page propertyPage1
    \section alphaFloor
    <code>MFloat CouplerLbFvEEMultiphase::m_alphaFloor </code>\n
    default = <code>0.0</code>\n \n
    Minimum acceptable value of alpha.\n
    Keywords: <i>EEMultiphase</i>
  */
  m_alphaFloor = 0.0;
  m_alphaFloor = Context::getBasicProperty<MFloat>("alphaFloor", AT_, &m_alphaFloor);

  /*! \page propertyPage1
    \section alphaCeil
    <code>MFloat CouplerLbFvEEMultiphase::m_alphaCeil </code>\n
    default = <code>1.0</code>\n \n
    Maximum acceptable value of alpha.\n
    Keywords: <i>EEMultiphase</i>
  */
  m_alphaCeil = 1.0;
  m_alphaCeil = Context::getBasicProperty<MFloat>("alphaCeil", AT_, &m_alphaCeil);

  /*! \page propertyPage1
    \section initAlphaMethod
    <code>MInt CouplerLbFvEEMultiphase::m_initAlphaMethod </code>\n
    default = <code>0</code>\n \n
    Method for the initial distribution of alphaGas.\n
    <li>0 : initialize alphaGas to be zero everywhere</li>
    <li>1 : initialize alphaGas to a constant value given in m_initialAlpha</li>
    <li>2 : initialize alphaGas to a constant value given in m_initialAlpha in a sphere with radius 0.5
            around the origin, else alphaInf</li>
    Keywords: <i>EEMultiphase</i>
  */
  m_initAlphaMethod = 0;
  if(Context::propertyExists("initAlphaMethod", 0)) {
    m_initAlphaMethod = Context::getBasicProperty<MInt>("initAlphaMethod", AT_, &m_initAlphaMethod);
  }

  /*! \page propertyPage1
    \section initialAlpha
    <code>MFloat CouplerLbFvEEMultiphase::m_initialAlpha </code>\n
    default = <code>0.0</code>\n \n
    Initial value of alpha in the domain.\n
    Keywords: <i>EEMultiphase</i>
  */
  m_initialAlpha = 0.0;
  if(Context::propertyExists("initialAlpha", 0)) {
    m_initialAlpha = Context::getBasicProperty<MFloat>("initialAlpha", AT_, &m_initialAlpha);
  }

  /*! \page propertyPage1
    \section alphaInf
    <code>MFloat CouplerLbFvEEMultiphase::m_alphaInf </code>\n
    default = <code>m_initialAlpha</code>\n \n
    Infinity value of alpha.\n
    Keywords: <i>EEMultiphase</i>
  */
  m_alphaInf = m_initialAlpha;
  m_alphaInf = Context::getBasicProperty<MFloat>("alphaInf", AT_, &m_alphaInf);

  /*! \page propertyPage1
    \section redistributeAlpha
    <code>MBool CouplerLbFvEEMultiphase::m_redistributeAlpha </code>\n
    default = <code>true</code>\n \n
    Redistribute invalid values for alpha.\n
    Keywords: <i>EEMultiphase</i>
  */
  m_redistributeAlpha = true;
  m_redistributeAlpha = Context::getBasicProperty<MBool>("redistributeAlpha", AT_, &m_redistributeAlpha);

  /*! \page propertyPage1
    \section disableSubstepAlphaRedist
    <code>MBool CouplerLbFvEEMultiphase::m_disableSubstepAlphaRedist </code>\n
    default = <code>true</code>\n \n
    Disable the redistribution of alpha in the RK-substeps.\n
    Keywords: <i>EEMultiphase</i>
  */
  m_disableSubstepAlphaRedist = true;
  m_disableSubstepAlphaRedist =
      Context::getBasicProperty<MBool>("disableSubstepAlphaRedist", AT_, &m_disableSubstepAlphaRedist);

  /*! \page propertyPage1
    \section updateAfterPropagation
    <code>MBool CouplerLbFvEEMultiphase::m_updateAfterPropagation </code>\n
    default = <code>true</code>\n \n
    Update the LB solver after propagation.\n
    Keywords: <i>EEMultiphase</i>
  */
  m_updateAfterPropagation = true;
  m_updateAfterPropagation =
      Context::getSolverProperty<MBool>("updateAfterPropagation", 1, AT_, &m_updateAfterPropagation);

  /*! \page propertyPage1
    \section EEMultiphaseInterpolationFactor
    <code>MFloat CouplerLbFvEEMultiphase::m_interpolationFactor </code>\n
    default = <code>0.5</code>\n \n
    This factor determines to which point in time the values of the other phase are inter/extrapolated.\n
    0.0 is the level of the old timestep, 1.0 is the level of the new timestep.
    Keywords: <i>EEMultiphase</i>
  */
  m_interpolationFactor = 0.5;
  m_interpolationFactor =
      Context::getSolverProperty<MFloat>("EEMultiphaseInterpolationFactor", 1, AT_, &m_interpolationFactor);

  /*! \page propertyPage1
    \section EEMultiphaseUpdateFVBC
    <code>MBool CouplerLbFvEEMultiphase::m_updateFVBC </code>\n
    default = <code>false</code>\n \n
    Update the FV BCs after transfer of the variables.\n
    Keywords: <i>EEMultiphase</i>
  */
  m_updateFVBC = false;
  m_updateFVBC = Context::getSolverProperty<MBool>("EEMultiphaseUpdateFVBC", 1, AT_, &m_updateFVBC);

  if(Context::propertyExists("gravityRefCoords", 0) && Context::propertyExists("depthCorrectionCoefficients", 0)) {
    for(MInt i = 0; i < 3; i++) {
      /*! \page propertyPage1
        \section gravityRefCoords
        <code>MFloat CouplerLbFvEEMultiphase::m_gravityRefCoords[nDim] </code>\n
        Reference Coordinates for density correction as a function of depth below the surface\n
        Keywords: <i>EEMultiphase</i>
      */
      m_gravityRefCoords[i] = Context::getSolverProperty<MFloat>("gravityRefCoords", 0, AT_, i);

      /*! \page propertyPage1
        \section depthCorrectionCoefficients
        <code>MFloat CouplerLbFvEEMultiphase::m_depthCorrectionCoefficients[nDim] </code>\n
        The depthCorrectioncoefficients are dimensionless coefficients for the change in density as a function of depth
        They are defined as (g L_ref)/(R T) with the compontents of the gravity vector g,
        the specific gas constant R (=287.06 J/kgK for air) and
        the Reference Temperature\n
        Keywords: <i>EEMultiphase</i>
      */
      m_depthCorrectionCoefficients[i] = Context::getSolverProperty<MFloat>("depthCorrectionCoefficients", 0, AT_, i);
    }
  } else {
    mTerm(1, AT_, "gravityRefCoords and depthCorrectionCoefficients are required for Euler-Euler multiphase!");
  }
}

template <MInt nDim, MInt nDist, class SysEqnLb, class SysEqnFv>
void CouplerLbFvEEMultiphase<nDim, nDist, SysEqnLb, SysEqnFv>::initDepthcorrection() {
  MFloat TInf = fvSolver().m_TInfinity;
  if(std::isnan(TInf)) {
    TInf = fvSolver().m_initialCondition == 465 || fvSolver().m_initialCondition == 9465
               ? 1.0
               : 1.0 / (1.0 + 0.5 * (fvSolver().m_gamma - 1.0) * POW2(fvSolver().m_Ma));
  }
  for(MInt fvCellId = 0; fvCellId < a_noFvCells(); fvCellId++) {
    MFloat depthCorrectionValue = 0.0;
    for(MInt i = 0; i < nDim; i++) {
      const MFloat deltaH = fvSolver().a_coordinate(fvCellId, i) - m_gravityRefCoords[i];
      depthCorrectionValue +=
          fvSolver().m_EEGas.liquidDensity * m_depthCorrectionCoefficients[i] * deltaH * TInf / fvSolver().m_gamma;
    }
    m_depthCorrectionValues[fvCellId] = depthCorrectionValue;
  }
}

template <MInt nDim, MInt nDist, class SysEqnLb, class SysEqnFv>
void CouplerLbFvEEMultiphase<nDim, nDist, SysEqnLb, SysEqnFv>::transferPressureLb2Fv(const MFloat rkAlpha,
                                                                                     const MBool update) {
  // conversion factor for pressure from LB to FV: 3 / (1 + (gamma-1)/2 * Ma_inf^2)^(gamma/(gamma-1))
  //                                             = 3 * T_inf ^ (gamma/(gamma-1))
  // conversion factor LB RHO to P: p = rho / (3 * gamma)
  // total conversion factor: T_inf ^ (gamma/(gamma-1)) / gamma

  for(MInt lbCellId = 0; lbCellId < a_noLbCells(); lbCellId++) {
    const MInt lvlDiff = lbSolver().maxLevel() - lbSolver().a_level(lbCellId);
    const MInt tsSinceUpdate = globalTimeStep % IPOW2(lvlDiff);
    const MFloat dt =
        m_interpolationFactor > 0.0 ? (rkAlpha * m_interpolationFactor - 1.0 + tsSinceUpdate) * FFPOW2(lvlDiff) : 0.0;
    MInt fvCellId = lb2fvId(lbCellId);
    if(fvCellId < 0) continue;

    if(!fvSolver().m_EEGas.depthCorrection) {
      fvSolver().a_pvariable(fvCellId, fvSolver().PV->P) =
          conversionLbFv.pressure * lbSolver().a_interpolatedVariable(lbCellId, lbSolver().PV->RHO, dt)
              / fvSolver().m_gamma
          + m_depthCorrectionValues[fvCellId];
    } else {
      fvSolver().a_pvariable(fvCellId, fvSolver().PV->P) =
          conversionLbFv.pressure * lbSolver().a_interpolatedVariable(lbCellId, lbSolver().PV->RHO, dt)
          / fvSolver().m_gamma;
    }
  }
  for(MInt bndIndex = 0; bndIndex < (MInt)(lbBndCnd().m_bndCndIds.size()); bndIndex++) {
    for(MInt bndId = lbBndCnd().m_bndCndOffsets[bndIndex]; bndId < lbBndCnd().m_bndCndOffsets[bndIndex + 1]; bndId++) {
      if(lbBndCnd().m_bndCells[bndId].m_isFluid) continue; // extrapolate only solid boundary cells from flow field
      const MInt lbCellId = lbBndCnd().m_bndCells[bndId].m_cellId;
      const MInt fvCellId = lb2fvId(lbCellId);
      if(fvCellId < 0) continue;
      for(MInt i = 0; i < 26; i++) { // look in all directions, including diagonals
        if(!lbSolver().a_hasNeighbor(lbCellId, i)) {
          MInt lbOppNId = lbSolver().c_neighborId(lbCellId, oppositeDirGrid[i]);
          if(lbOppNId >= 0 && !lbSolver().a_isBndryCell(lbOppNId)) {
            MInt fvOppNId = lb2fvId(lbOppNId);
            if(fvOppNId >= a_noFvCells() || fvOppNId < 0) {
              break;
            }
            if(!fvSolver().m_EEGas.depthCorrection) {
              fvSolver().a_pvariable(fvCellId, fvSolver().PV->P) = fvSolver().a_pvariable(fvOppNId, fvSolver().PV->P)
                                                                   - m_depthCorrectionValues[fvOppNId]
                                                                   + m_depthCorrectionValues[fvCellId];
            } else {
              fvSolver().a_pvariable(fvCellId, fvSolver().PV->P) = fvSolver().a_pvariable(fvOppNId, fvSolver().PV->P);
            }
            break;
          }
        }
      }
    }
  }

  if(update) {
    fvSolver().computePV();
    fvSolver().exchange();
    fvSolver().m_fvBndryCnd->updateGhostCellVariables();

    // apply von Neumann boundary conditions
    fvSolver().m_fvBndryCnd->applyNeumannBoundaryCondition();
  }
}

template <MInt nDim, MInt nDist, class SysEqnLb, class SysEqnFv>
template <MBool updateGradients>
void CouplerLbFvEEMultiphase<nDim, nDist, SysEqnLb, SysEqnFv>::transferULb2Fv(const MFloat rkAlpha) {
  static MBool firstCall = true;

  // conversion factor for velocity from LB to FV: sqrt(3 / (1 + (gamma-1)/2 * Ma_inf^2))
  const MFloat ufactor = conversionLbFv.velocity;

  // derivatives need an additional conversion factor
  const MFloat dufactor = conversionLbFv.length;

  // swap the pointers instead of copying all values from m_EEGas.uOtherPhase to m_EEGas.uOtherPhaseOld
  MFloat** tempPointer = fvSolver().m_EEGas.uOtherPhaseOld;
  fvSolver().m_EEGas.uOtherPhaseOld = fvSolver().m_EEGas.uOtherPhase;
  fvSolver().m_EEGas.uOtherPhase = tempPointer;

  lbSolver().exchange();

  // update the velocity arrays
  const MBool interpolate = m_interpolationFactor > 0.0;
  for(MInt lbCellId = 0; lbCellId < a_noLbCells(); lbCellId++) {
    MInt fvCellId = lb2fvId(lbCellId);
    if(fvCellId < 0) continue;

    if(interpolate) {
      const MInt lvlDiff = lbSolver().maxLevel() - lbSolver().a_level(lbCellId);
      const MInt tsSinceUpdate = globalTimeStep % IPOW2(lvlDiff);
      const MFloat dt = (rkAlpha * m_interpolationFactor - 1.0 + tsSinceUpdate) * FFPOW2(lvlDiff);
      // since the gradients are not updated for each RK step, they are interpolated for the whole timestep
      for(MInt i = 0; i < nDim; i++) {
        fvSolver().a_uOtherPhase(fvCellId, i) = ufactor * lbSolver().a_interpolatedVariable(lbCellId, i, dt);
        IF_CONSTEXPR(updateGradients) {
          const MFloat dtGrad = (m_interpolationFactor - 1.0 + tsSinceUpdate) * FFPOW2(lvlDiff);
          for(MInt j = 0; j < nDim; j++) {
            fvSolver().a_gradUOtherPhase(fvCellId, i, j) =
                ufactor * dufactor * lbSolver().calculateInterpolatedDerivative(lbCellId, i, j, dtGrad);
          }
        }
      }
    } else {
      for(MInt i = 0; i < nDim; i++) {
        fvSolver().a_uOtherPhase(fvCellId, i) = ufactor * lbSolver().a_variable(lbCellId, i);
        IF_CONSTEXPR(updateGradients) {
          for(MInt j = 0; j < nDim; j++) {
            fvSolver().a_gradUOtherPhase(fvCellId, i, j) =
                ufactor * dufactor * lbSolver().calculateDerivative(lbCellId, i, j);
          }
        }
      }
    }
    IF_CONSTEXPR(updateGradients) {
      fvSolver().a_vortOtherPhase(fvCellId, 0) =
          fvSolver().a_gradUOtherPhase(fvCellId, 2, 1) - fvSolver().a_gradUOtherPhase(fvCellId, 1, 2);
      fvSolver().a_vortOtherPhase(fvCellId, 1) =
          fvSolver().a_gradUOtherPhase(fvCellId, 0, 2) - fvSolver().a_gradUOtherPhase(fvCellId, 2, 0);
      fvSolver().a_vortOtherPhase(fvCellId, 2) =
          fvSolver().a_gradUOtherPhase(fvCellId, 1, 0) - fvSolver().a_gradUOtherPhase(fvCellId, 0, 1);
    }
  }
  if(firstCall) {
    for(MInt lbCellId = 0; lbCellId < a_noLbCells(); lbCellId++) {
      MInt fvCellId = lb2fvId(lbCellId);
      if(fvCellId < 0) continue;
      for(MInt i = 0; i < 3; i++) {
        fvSolver().a_uOtherPhaseOld(fvCellId, i) = fvSolver().a_uOtherPhase(fvCellId, i);
      }
    }
    firstCall = false;
  }
}

template <MInt nDim, MInt nDist, class SysEqnLb, class SysEqnFv>
void CouplerLbFvEEMultiphase<nDim, nDist, SysEqnLb, SysEqnFv>::transferUFv2Lb() {
  // conversion factor for velocity from LB to FV: sqrt((1 + (gamma-1)/2 * Ma_inf^2) / 3)

  for(MInt fvCellId = 0; fvCellId < a_noFvCells(); fvCellId++) {
    MInt lbCellId = fv2lbId(fvCellId);
    if(lbCellId < 0) continue;

    lbSolver().a_uOtherPhase(lbCellId, 0) =
        conversionFvLb.velocity * fvSolver().a_pvariable(fvCellId, fvSolver().PV->U);
    lbSolver().a_uOtherPhase(lbCellId, 1) =
        conversionFvLb.velocity * fvSolver().a_pvariable(fvCellId, fvSolver().PV->V);
    lbSolver().a_uOtherPhase(lbCellId, 2) =
        conversionFvLb.velocity * fvSolver().a_pvariable(fvCellId, fvSolver().PV->W);
  }
}

template <MInt nDim, MInt nDist, class SysEqnLb, class SysEqnFv>
void CouplerLbFvEEMultiphase<nDim, nDist, SysEqnLb, SysEqnFv>::transferNuLb2Fv(const MFloat rkAlpha) {
  for(MInt lbCellId = 0; lbCellId < a_noLbCells(); lbCellId++) {
    const MInt lvlDiff = lbSolver().maxLevel() - lbSolver().a_level(lbCellId);
    const MInt tsSinceUpdate = globalTimeStep % IPOW2(lvlDiff);
    const MFloat dt =
        m_interpolationFactor > 0.0 ? (rkAlpha * m_interpolationFactor - 1.0 + tsSinceUpdate) * FFPOW2(lvlDiff) : 0.0;
    const MInt fvCellId = lb2fvId(lbCellId);
    if(fvCellId < 0) continue;
    if(fvSolver().a_isBndryGhostCell(fvCellId)) continue;

    fvSolver().a_nuEffOtherPhase(fvCellId) =
        ((1.0 + dt) * lbSolver().a_nu(lbCellId) - dt * lbSolver().a_oldNu(lbCellId)) * conversionLbFv.viscosity;
    fvSolver().a_nuTOtherPhase(fvCellId) =
        ((1.0 + dt) * lbSolver().a_nuT(lbCellId) - dt * lbSolver().a_oldNuT(lbCellId)) * conversionLbFv.viscosity;
  }
}

template <MInt nDim, MInt nDist, class SysEqnLb, class SysEqnFv>
void CouplerLbFvEEMultiphase<nDim, nDist, SysEqnLb, SysEqnFv>::transferAlphaFv2Lb() {
  for(MInt fvCellId = 0; fvCellId < a_noFvCells(); fvCellId++) {
    const MInt lbCellId = fv2lbId(fvCellId);
    if(lbCellId < 0) continue;

    lbSolver().a_alphaGas(lbCellId) = fvSolver().a_alphaGas(fvCellId);
  }
}

template <MInt nDim, MInt nDist, class SysEqnLb, class SysEqnFv>
void CouplerLbFvEEMultiphase<nDim, nDist, SysEqnLb, SysEqnFv>::revertLbVariables() {
  for(MInt cellId = 0; cellId < a_noLbCells(); cellId++) {
    lbSolver().a_variable(cellId, lbSolver().PV->RHO) = lbSolver().a_oldVariable(cellId, lbSolver().PV->RHO);
    lbSolver().a_variable(cellId, lbSolver().PV->U) = lbSolver().a_oldVariable(cellId, lbSolver().PV->U);
    lbSolver().a_variable(cellId, lbSolver().PV->V) = lbSolver().a_oldVariable(cellId, lbSolver().PV->V);
    lbSolver().a_variable(cellId, lbSolver().PV->W) = lbSolver().a_oldVariable(cellId, lbSolver().PV->W);
  }
}

template <MInt nDim, MInt nDist, class SysEqnLb, class SysEqnFv>
void CouplerLbFvEEMultiphase<nDim, nDist, SysEqnLb, SysEqnFv>::revertLbDistributions() {
  for(MInt cellId = 0; cellId < a_noLbCells(); cellId++) {
    for(MInt distr = 0; distr < lbSolver().m_noDistributions; distr++) {
      lbSolver().a_oldDistribution(cellId, distr) = lbSolver().a_previousDistribution(cellId, distr);
    }
  }
}

template <MInt nDim, MInt nDist, class SysEqnLb, class SysEqnFv>
void CouplerLbFvEEMultiphase<nDim, nDist, SysEqnLb, SysEqnFv>::revertLbOldVariables() {
  for(MInt cellId = 0; cellId < a_noLbCells(); cellId++) {
    lbSolver().a_oldVariable(cellId, lbSolver().PV->RHO) = lbSolver().a_previousVariable(cellId, lbSolver().PV->RHO);
    lbSolver().a_oldVariable(cellId, lbSolver().PV->U) = lbSolver().a_previousVariable(cellId, lbSolver().PV->U);
    lbSolver().a_oldVariable(cellId, lbSolver().PV->V) = lbSolver().a_previousVariable(cellId, lbSolver().PV->V);
    lbSolver().a_oldVariable(cellId, lbSolver().PV->W) = lbSolver().a_previousVariable(cellId, lbSolver().PV->W);
  }
}

template <MInt nDim, MInt nDist, class SysEqnLb, class SysEqnFv>
void CouplerLbFvEEMultiphase<nDim, nDist, SysEqnLb, SysEqnFv>::revertLb() {
  if(m_alphaConvergenceCheck <= 0) {
    mTerm(1, AT_, "Didn't store the variables to revert the LB solver!");
  }
  revertLbVariables();
  revertLbDistributions();
  if(m_updateAfterPropagation) revertLbOldVariables();
  if(lbSolver().noDomains() > 1) {
    lbSolver().exchange();
    lbSolver().exchangeOldDistributions();
  }
}

template <MInt nDim, MInt nDist, class SysEqnLb, class SysEqnFv>
void CouplerLbFvEEMultiphase<nDim, nDist, SysEqnLb, SysEqnFv>::revertFvVariables() {
  for(MInt cellId = 0; cellId < a_noFvCells(); cellId++) {
    for(MInt varId = 0; varId < fvSolver().noVariables(); varId++) {
      fvSolver().a_variable(cellId, varId) = fvSolver().a_oldVariable(cellId, varId);
    }
  }
  fvSolver().computePV();
}

template <MInt nDim, MInt nDist, class SysEqnLb, class SysEqnFv>
void CouplerLbFvEEMultiphase<nDim, nDist, SysEqnLb, SysEqnFv>::revertFv() {
  revertFvVariables();
  fvSolver().revertTimestep();
  fvSolver().exchange();
}

template <MInt nDim, MInt nDist, class SysEqnLb, class SysEqnFv>
MBool CouplerLbFvEEMultiphase<nDim, nDist, SysEqnLb, SysEqnFv>::checkAlphaConverged() {
  MBool converged = false;
  switch(m_alphaConvergenceCheck) {
    case 0:
      converged = true;
      break;

    case 1:
      converged = true;
      for(MInt cellId = 0; cellId < a_noFvCells(); cellId++) {
        if(fv2lbId(cellId) < 0) continue;
        if(m_epsAlpha < std::abs(fvSolver().a_alphaGas(cellId) - lbSolver().a_alphaGas(fv2lbId(cellId)))) {
          converged = false;
          break;
        }
      }
      MPI_Allreduce(MPI_IN_PLACE, &converged, 1, maia::type_traits<MBool>::mpiType(), MPI_LAND, fvSolver().mpiComm(),
                    AT_, "MPI_IN_PLACE", "converged");
      break;

    default:
      mTerm(1, AT_, "not a valid alphaConvergenceCheck!");
      break;
  }
  return converged;
}

/// \brief find cells with invalid alpha values and redistribute mass from/to neighboring cells
///        to raise/lower alpha value
///
/// \author Daniel Lauwers
/// \date 2020-11-30
template <MInt nDim, MInt nDist, class SysEqnLb, class SysEqnFv>
void CouplerLbFvEEMultiphase<nDim, nDist, SysEqnLb, SysEqnFv>::correctInvalidAlpha() {
  const MInt noCVars = fvSolver().CV->noVariables;

  fvSolver().computePV();
  fvSolver().exchange();

  // Request conservative variables of halo cells
  // modified from smallcellcorrection()
  {
    MIntScratchSpace sendBufferCnts(mMax(1, fvSolver().noNeighborDomains()), AT_, "sendBufferCnts");
    MIntScratchSpace recvBufferCnts(mMax(1, fvSolver().noNeighborDomains()), AT_, "recvBufferCnts");
    ScratchSpace<MPI_Request> sendReq(mMax(1, fvSolver().noNeighborDomains()), AT_, "sendReq");
    ScratchSpace<MPI_Request> recvReq(mMax(1, fvSolver().noNeighborDomains()), AT_, "recvReq");
    sendReq.fill(MPI_REQUEST_NULL);
    recvReq.fill(MPI_REQUEST_NULL);
    sendBufferCnts.fill(0);
    recvBufferCnts.fill(0);
    for(MInt i = 0; i < fvSolver().noNeighborDomains(); i++) {
      sendBufferCnts(i) = noCVars * fvSolver().m_noMaxLevelWindowCells[i];
      recvBufferCnts(i) = noCVars * fvSolver().m_noMaxLevelHaloCells[i];
      MInt sendBufferCounter = 0;
      for(MInt j = 0; j < fvSolver().m_noMaxLevelWindowCells[i]; j++) {
        MInt cellId = fvSolver().m_maxLevelWindowCells[i][j];
        for(MInt v = 0; v < noCVars; v++) {
          fvSolver().m_sendBuffers[i][sendBufferCounter] = fvSolver().a_variable(cellId, v);
          sendBufferCounter++;
        }
      }
    }
    MInt sendCnt = 0;
    MInt recvCnt = 0;
    if(fvSolver().m_nonBlockingComm) {
      for(MInt i = 0; i < fvSolver().noNeighborDomains(); i++) {
        if(sendBufferCnts(i) == 0) continue;
        ASSERT(sendBufferCnts(i) <= fvSolver().m_noMaxLevelWindowCells[i] * fvSolver().m_dataBlockSize, "");
        MPI_Isend(fvSolver().m_sendBuffers[i], sendBufferCnts(i), MPI_DOUBLE, fvSolver().neighborDomain(i), 12,
                  fvSolver().mpiComm(), &sendReq[sendCnt], AT_, "m_sendBuffers[" + std::to_string(i) + "],0");
        sendCnt++;
      }
      for(MInt i = 0; i < fvSolver().noNeighborDomains(); i++) {
        if(recvBufferCnts(i) == 0) continue;
        ASSERT(recvBufferCnts(i) <= fvSolver().m_noMaxLevelHaloCells[i] * fvSolver().m_dataBlockSize, "");
        MPI_Irecv(fvSolver().m_receiveBuffers[i], recvBufferCnts(i), MPI_DOUBLE, fvSolver().neighborDomain(i), 12,
                  fvSolver().mpiComm(), &recvReq[recvCnt], AT_, "m_receiveBuffers[" + std::to_string(i) + "],0");
        recvCnt++;
      }
      if(recvCnt > 0) MPI_Waitall(recvCnt, &recvReq[0], MPI_STATUSES_IGNORE, AT_);
    } else {
      for(MInt i = 0; i < fvSolver().noNeighborDomains(); i++) {
        if(sendBufferCnts(i) == 0) continue;
        ASSERT(sendBufferCnts(i) <= fvSolver().m_noMaxLevelWindowCells[i] * fvSolver().m_dataBlockSize, "");
        MPI_Issend(fvSolver().m_sendBuffers[i], sendBufferCnts(i), MPI_DOUBLE, fvSolver().neighborDomain(i), 12,
                   fvSolver().mpiComm(), &sendReq[sendCnt], AT_, "m_sendBuffers[" + std::to_string(i) + "],0");
        sendCnt++;
      }
      for(MInt i = 0; i < fvSolver().noNeighborDomains(); i++) {
        if(recvBufferCnts(i) == 0) continue;
        ASSERT(recvBufferCnts(i) <= fvSolver().m_noMaxLevelHaloCells[i] * fvSolver().m_dataBlockSize, "");
        MPI_Recv(fvSolver().m_receiveBuffers[i], recvBufferCnts(i), MPI_DOUBLE, fvSolver().neighborDomain(i), 12,
                 fvSolver().mpiComm(), MPI_STATUS_IGNORE, AT_, "m_receiveBuffers[" + std::to_string(i) + "],0");
        recvCnt++;
      }
    }
    for(MInt i = 0; i < fvSolver().noNeighborDomains(); i++) {
      MInt recvBufferCounter = 0;
      for(MInt j = 0; j < fvSolver().m_noMaxLevelHaloCells[i]; j++) {
        MInt cellId = fvSolver().m_maxLevelHaloCells[i][j];
        for(MInt v = 0; v < noCVars; v++) {
          fvSolver().a_variable(cellId, v) = fvSolver().m_receiveBuffers[i][recvBufferCounter];
          recvBufferCounter++;
        }
      }
    }
    if(sendCnt > 0) MPI_Waitall(sendCnt, &sendReq[0], MPI_STATUSES_IGNORE, AT_);
  }

  // find cells with invalid alpha values and redistribute mass and momentum from/to them
  // TODO labels:COUPLER,toenhance replace struct by std::vector or something else
  struct transferCV {
    MFloat rhoAlpha;
    MFloat rhoUAlpha;
    MFloat rhoVAlpha;
    MFloat rhoWAlpha;
  };
  constexpr MInt noTransferCV = 4;
  std::map<MInt, transferCV>
      haloCorrections; // save which Halo cells are altered to communicate the changes to their domains
  for(MInt cellId = 0; cellId < fvSolver().m_bndryGhostCellsOffset; cellId++) {
    if(fvSolver().a_isHalo(cellId)) continue;
    if(!fvSolver().c_isLeafCell(cellId)) continue;

    const MFloat origAlpha = fvSolver().a_pvariable(cellId, fvSolver().PV->A);

    if(origAlpha < (m_alphaFloor - 1.0e-8)) {
      const MFloat origMass = fvSolver().a_variable(cellId, fvSolver().CV->A_RHO) * fvSolver().a_cellVolume(cellId);
      const MFloat deltaMass = 0.0 - origMass;
      const std::vector<MInt> neighbors = findRedistCells(cellId, true, 0.0);
      std::vector<MInt> redistNeighbors = neighbors;
      MInt noRedistNbs = 0;
      MFloat massNeedPerCell = 0.0;
      while(true) {
        MBool change = false;
        noRedistNbs = redistNeighbors.size();
        massNeedPerCell = deltaMass / noRedistNbs;
        for(auto it = redistNeighbors.begin(); it != redistNeighbors.end(); ++it) {
          const MInt candidateId = *it;
          if(fvSolver().a_variable(candidateId, fvSolver().CV->A_RHO) * fvSolver().a_cellVolume(candidateId)
             < massNeedPerCell) {
            change = true;
            redistNeighbors.erase(it);
            break;
          }
        }
        if(change) continue;
        break;
      }
      MFloat momentumTransfer[3] = {0.0, 0.0, 0.0};
      if(!redistNeighbors.empty()) {
        for(auto redistId : redistNeighbors) {
          const MFloat delRho = massNeedPerCell / fvSolver().a_cellVolume(redistId);
          const MFloat delRhoVV[3] = {delRho * fvSolver().a_pvariable(redistId, fvSolver().PV->VV[0]),
                                      delRho * fvSolver().a_pvariable(redistId, fvSolver().PV->VV[1]),
                                      delRho * fvSolver().a_pvariable(redistId, fvSolver().PV->VV[2])};
          fvSolver().a_variable(redistId, fvSolver().CV->A_RHO) -= delRho;
          for(MInt dir = 0; dir < 3; dir++) {
            fvSolver().a_variable(redistId, fvSolver().CV->A_RHO_VV[dir]) -= delRhoVV[dir];
            momentumTransfer[dir] += delRhoVV[dir] * fvSolver().a_cellVolume(redistId);
          }
          fvSolver().setPrimitiveVariables(redistId);

          if(fvSolver().a_isHalo(redistId)) {
            const transferCV delCV = {
                -delRho,
                -delRhoVV[0],
                -delRhoVV[1],
                -delRhoVV[2],
            };
            auto success = haloCorrections.insert(std::pair<MInt, transferCV>(redistId, delCV));
            if(!success.second) { // redistId already in map
              success.first->second.rhoAlpha += delCV.rhoAlpha;
              success.first->second.rhoUAlpha += delCV.rhoUAlpha;
              success.first->second.rhoVAlpha += delCV.rhoVAlpha;
              success.first->second.rhoWAlpha += delCV.rhoWAlpha;
            }
          }
        }
        fvSolver().a_variable(cellId, fvSolver().CV->A_RHO) +=
            massNeedPerCell * noRedistNbs / fvSolver().a_cellVolume(cellId);
        for(MInt dir = 0; dir < 3; dir++) {
          fvSolver().a_variable(cellId, fvSolver().CV->A_RHO_VV[dir]) +=
              momentumTransfer[dir] / fvSolver().a_cellVolume(cellId);
        }
        fvSolver().setPrimitiveVariables(cellId);
      } else {
        // if there is not enough mass in any cell to get the negative cell to 0, take as much, as you can from cells
        redistNeighbors = neighbors;

        MFloat massTransfer = 0.0;
        for(auto redistId : redistNeighbors) {
          const MFloat massAvail =
              fvSolver().a_variable(redistId, fvSolver().CV->A_RHO) * fvSolver().a_cellVolume(redistId);
          const MFloat delRho = massAvail / fvSolver().a_cellVolume(redistId);
          const MFloat delRhoVV[3] = {delRho * fvSolver().a_pvariable(redistId, fvSolver().PV->VV[0]),
                                      delRho * fvSolver().a_pvariable(redistId, fvSolver().PV->VV[1]),
                                      delRho * fvSolver().a_pvariable(redistId, fvSolver().PV->VV[2])};
          massTransfer += massAvail;
          fvSolver().a_variable(redistId, fvSolver().CV->A_RHO) -= delRho;
          for(MInt dir = 0; dir < 3; dir++) {
            fvSolver().a_variable(redistId, fvSolver().CV->A_RHO_VV[dir]) -= delRhoVV[dir];
            momentumTransfer[dir] += delRhoVV[dir] * fvSolver().a_cellVolume(redistId);
          }
          fvSolver().setPrimitiveVariables(redistId);

          if(fvSolver().a_isHalo(redistId)) {
            const transferCV delCV = {
                -delRho,
                -delRhoVV[0],
                -delRhoVV[1],
                -delRhoVV[2],
            };
            auto success = haloCorrections.insert(std::pair<MInt, transferCV>(redistId, delCV));
            if(!success.second) { // redistId already in map
              success.first->second.rhoAlpha += delCV.rhoAlpha;
              success.first->second.rhoUAlpha += delCV.rhoUAlpha;
              success.first->second.rhoVAlpha += delCV.rhoVAlpha;
              success.first->second.rhoWAlpha += delCV.rhoWAlpha;
            }
          }
        }
        fvSolver().a_variable(cellId, fvSolver().CV->A_RHO) += massTransfer / fvSolver().a_cellVolume(cellId);
        for(MInt dir = 0; dir < 3; dir++) {
          fvSolver().a_variable(cellId, fvSolver().CV->A_RHO_VV[dir]) +=
              momentumTransfer[dir] / fvSolver().a_cellVolume(cellId);
        }
        fvSolver().setPrimitiveVariables(cellId);
      }
    } else if(origAlpha > (m_alphaCeil + 1.0e-8)) {
      const MFloat origMass = fvSolver().a_variable(cellId, fvSolver().CV->A_RHO) * fvSolver().a_cellVolume(cellId);
      const MFloat deltaMass =
          m_alphaCeil * fvSolver().a_pvariable(cellId, fvSolver().PV->RHO) * fvSolver().a_cellVolume(cellId) - origMass;
      std::vector<MInt> neighbors = findRedistCells(cellId, false, m_alphaCeil);
      MBool isOutlet = false;

      if(fvSolver().a_bndryId(cellId) >= 0) {
        const MInt bndryId = fvSolver().a_bndryId(cellId);
        for(MInt srfc = 0; srfc < fvSolver().m_bndryCells->a[bndryId].m_noSrfcs; srfc++) {
          const MInt bndrCndId = fvSolver().m_fvBndryCnd->m_bndryCell[bndryId].m_srfcs[srfc]->m_bndryCndId;
          if(bndrCndId == 1002) {
            isOutlet = true;
            break;
          }
        }
      }
      if(isOutlet) {
        for(MUint dir = 0; dir < 26; dir++) { // loop over ALL the directions!
          const MInt neighborId = fvSolver().grid().neighborList(cellId, dir);
          if(fvSolver().a_isBndryGhostCell(neighborId)) {
            neighbors.push_back(neighborId);
          }
        }
      }
      std::vector<MInt> redistNeighbors = neighbors;

      MInt noRedistNbs = 0;
      MFloat massExcessPerCell = 0.0;
      while(true) {
        MBool change = false;
        noRedistNbs = redistNeighbors.size();
        massExcessPerCell = -deltaMass / noRedistNbs;
        for(auto it = redistNeighbors.begin(); it != redistNeighbors.end(); ++it) {
          const MInt candidateId = *it;
          if(!fvSolver().a_isBndryGhostCell(candidateId)
             && ((fvSolver().a_variable(candidateId, fvSolver().CV->A_RHO)
                  + massExcessPerCell / fvSolver().a_cellVolume(candidateId))
                 > fvSolver().a_pvariable(candidateId, fvSolver().PV->RHO) * m_alphaCeil)) {
            change = true;
            redistNeighbors.erase(it);
            break;
          }
        }
        if(change) continue;
        break;
      }
      MFloat momentumTransfer[3] = {0.0, 0.0, 0.0};
      if(!redistNeighbors.empty() || isOutlet) {
        if(redistNeighbors.empty()) {
          noRedistNbs = 1;
          massExcessPerCell = -deltaMass / noRedistNbs;
          for(MInt dir = 0; dir < 3; dir++) {
            momentumTransfer[dir] += massExcessPerCell * fvSolver().a_pvariable(cellId, fvSolver().PV->VV[dir]);
          }
        }
        for(auto redistId : redistNeighbors) {
          if(fvSolver().a_isBndryGhostCell(redistId)) continue;
          const MFloat delRho = massExcessPerCell / fvSolver().a_cellVolume(redistId);
          const MFloat delRhoVV[3] = {delRho * fvSolver().a_pvariable(redistId, fvSolver().PV->VV[0]),
                                      delRho * fvSolver().a_pvariable(redistId, fvSolver().PV->VV[1]),
                                      delRho * fvSolver().a_pvariable(redistId, fvSolver().PV->VV[2])};
          fvSolver().a_variable(redistId, fvSolver().CV->A_RHO) += delRho;
          for(MInt dir = 0; dir < 3; dir++) {
            fvSolver().a_variable(redistId, fvSolver().CV->A_RHO_VV[dir]) += delRhoVV[dir];
            momentumTransfer[dir] += delRhoVV[dir] * fvSolver().a_cellVolume(redistId);
          }
          fvSolver().setPrimitiveVariables(redistId);

          if(fvSolver().a_isHalo(redistId)) {
            const transferCV delCV = {
                delRho,
                delRhoVV[0],
                delRhoVV[1],
                delRhoVV[2],
            };
            auto success = haloCorrections.insert(std::pair<MInt, transferCV>(redistId, delCV));
            if(!success.second) { // redistId already in map
              success.first->second.rhoAlpha += delCV.rhoAlpha;
              success.first->second.rhoUAlpha += delCV.rhoUAlpha;
              success.first->second.rhoVAlpha += delCV.rhoVAlpha;
              success.first->second.rhoWAlpha += delCV.rhoWAlpha;
            }
          }
        }
        fvSolver().a_variable(cellId, fvSolver().CV->A_RHO) -=
            massExcessPerCell / fvSolver().a_cellVolume(cellId) * noRedistNbs;
        for(MInt dir = 0; dir < 3; dir++) {
          fvSolver().a_variable(cellId, fvSolver().CV->A_RHO_VV[dir]) -=
              momentumTransfer[dir] / fvSolver().a_cellVolume(cellId);
        }
        fvSolver().setPrimitiveVariables(cellId);
      } else {
        // if no neighbor can take enough mass to bring our cell under alphaCeil, give them as much as they can take
        redistNeighbors = neighbors;
        if(redistNeighbors.empty()) {
          redistNeighbors = findRedistCells(cellId, false, std::numeric_limits<MFloat>::max());
        }
        MFloat massTransfer = 0.0;
        for(auto redistId : redistNeighbors) {
          const MFloat massCapacityAvail = (fvSolver().a_pvariable(redistId, fvSolver().PV->RHO) * m_alphaCeil
                                            - fvSolver().a_variable(redistId, fvSolver().CV->A_RHO))
                                           * fvSolver().a_cellVolume(redistId);
          if(massCapacityAvail < 0.0) continue;
          const MFloat delRho = massCapacityAvail / fvSolver().a_cellVolume(redistId);
          const MFloat delRhoVV[3] = {delRho * fvSolver().a_pvariable(redistId, fvSolver().PV->VV[0]),
                                      delRho * fvSolver().a_pvariable(redistId, fvSolver().PV->VV[1]),
                                      delRho * fvSolver().a_pvariable(redistId, fvSolver().PV->VV[2])};
          massTransfer += massCapacityAvail;
          fvSolver().a_variable(redistId, fvSolver().CV->A_RHO) += delRho;
          for(MInt dir = 0; dir < 3; dir++) {
            fvSolver().a_variable(redistId, fvSolver().CV->A_RHO_VV[dir]) += delRhoVV[dir];
            momentumTransfer[dir] += delRhoVV[dir] * fvSolver().a_cellVolume(redistId);
          }
          fvSolver().setPrimitiveVariables(redistId);

          if(fvSolver().a_isHalo(redistId)) {
            const transferCV delCV = {
                delRho,
                delRhoVV[0],
                delRhoVV[1],
                delRhoVV[2],
            };
            auto success = haloCorrections.insert(std::pair<MInt, transferCV>(redistId, delCV));
            if(!success.second) { // redistId already in map
              success.first->second.rhoAlpha += delCV.rhoAlpha;
              success.first->second.rhoUAlpha += delCV.rhoUAlpha;
              success.first->second.rhoVAlpha += delCV.rhoVAlpha;
              success.first->second.rhoWAlpha += delCV.rhoWAlpha;
            }
          }
        }
        fvSolver().a_variable(cellId, fvSolver().CV->A_RHO) -= massTransfer / fvSolver().a_cellVolume(cellId);
        for(MInt dir = 0; dir < 3; dir++) {
          fvSolver().a_variable(cellId, fvSolver().CV->A_RHO_VV[dir]) -=
              momentumTransfer[dir] / fvSolver().a_cellVolume(cellId);
        }
        fvSolver().setPrimitiveVariables(cellId);
      }
    }
  }

  // Communicate the changes in the CVs to the neighboring domains
  // Since this function changed the CV in the halo cells, m_send- and m_receiveBuffer are swapped in this communication
  {
    if(noTransferCV > fvSolver().m_dataBlockSize)
      mTerm(1, AT_, "The send- and receive Buffers are too small for this communication!");
    // gather
    for(MInt i = 0; i < fvSolver().noNeighborDomains(); i++) {
      MInt receiveBufferCounter = 0;
      for(MInt j = 0; j < fvSolver().m_noMaxLevelHaloCells[i]; j++) {
        auto it = haloCorrections.find(fvSolver().m_maxLevelHaloCells[i][j]);
        if(it != haloCorrections.end()) {
          fvSolver().m_receiveBuffers[i][receiveBufferCounter + fvSolver().CV->A_RHO] = it->second.rhoAlpha;
          fvSolver().m_receiveBuffers[i][receiveBufferCounter + fvSolver().CV->A_RHO_U] = it->second.rhoUAlpha;
          fvSolver().m_receiveBuffers[i][receiveBufferCounter + fvSolver().CV->A_RHO_V] = it->second.rhoVAlpha;
          fvSolver().m_receiveBuffers[i][receiveBufferCounter + fvSolver().CV->A_RHO_W] = it->second.rhoWAlpha;
        } else {
          fvSolver().m_receiveBuffers[i][receiveBufferCounter + fvSolver().CV->A_RHO] = 0.0;
          fvSolver().m_receiveBuffers[i][receiveBufferCounter + fvSolver().CV->A_RHO_U] = 0.0;
          fvSolver().m_receiveBuffers[i][receiveBufferCounter + fvSolver().CV->A_RHO_V] = 0.0;
          fvSolver().m_receiveBuffers[i][receiveBufferCounter + fvSolver().CV->A_RHO_W] = 0.0;
        }
        receiveBufferCounter += noTransferCV;
      }
    }

    // send
    MInt bufSize = 0;
    for(MInt i = 0; i < fvSolver().noNeighborDomains(); i++) {
      bufSize = fvSolver().m_noMaxLevelHaloCells[i] * noTransferCV;
      MPI_Issend(fvSolver().m_receiveBuffers[i], bufSize, MPI_DOUBLE, fvSolver().neighborDomain(i), 0,
                 fvSolver().mpiComm(), &fvSolver().m_mpi_request[i], AT_,
                 "m_receiveBuffers[" + std::to_string(i) + "],1");
    }

    // receive
    MPI_Status status;
    for(MInt i = 0; i < fvSolver().noNeighborDomains(); i++) {
      bufSize = fvSolver().m_noMaxLevelWindowCells[i] * noTransferCV;
      MPI_Recv(fvSolver().m_sendBuffers[i], bufSize, MPI_DOUBLE, fvSolver().neighborDomain(i), 0, fvSolver().mpiComm(),
               &status, AT_, "m_sendBuffers[" + std::to_string(i) + "],1");
    }
    for(MInt i = 0; i < fvSolver().noNeighborDomains(); i++) {
      MPI_Wait(&fvSolver().m_mpi_request[i], &status, AT_);
    }

    // scatter
    MInt sendBufferCounter = 0;
    for(MInt i = 0; i < fvSolver().noNeighborDomains(); i++) {
      sendBufferCounter = 0;
      for(MInt j = 0; j < fvSolver().m_noMaxLevelWindowCells[i]; j++) {
        fvSolver().a_variable(fvSolver().m_maxLevelWindowCells[i][j], fvSolver().CV->A_RHO) +=
            fvSolver().m_sendBuffers[i][sendBufferCounter + fvSolver().CV->A_RHO];
        for(MInt dir = 0; dir < nDim; dir++) {
          fvSolver().a_variable(fvSolver().m_maxLevelWindowCells[i][j], fvSolver().CV->A_RHO_VV[dir]) +=
              fvSolver().m_sendBuffers[i][sendBufferCounter + fvSolver().CV->A_RHO_VV[dir]];
        }
        sendBufferCounter += noTransferCV;
      }
    }
  }
  fvSolver().computePV();
  fvSolver().exchange();
}

/** \brief find neighbor Cells to cellId, that are candidates for alpha corrections
 *
 * \author Daniel Lauwers
 *
 * \param[in] cellId the cell id that needs to be corrected
 * \param[in] searchUp if true, search cells with a gas mass > limit, else gas mass < limit
 * \param[in] limit the limit of alpha to be considered
 * \return vector of Ids of candidates for redistribution
 **/
template <MInt nDim, MInt nDist, class SysEqnLb, class SysEqnFv>
std::vector<MInt> CouplerLbFvEEMultiphase<nDim, nDist, SysEqnLb, SysEqnFv>::findRedistCells(const MInt cellId,
                                                                                            const MBool searchUp,
                                                                                            const MFloat limit) {
  std::vector<MInt> redistNeighbors;
  redistNeighbors.reserve(56);

  MIntScratchSpace adjacentLeafCells(56, AT_, "adjacentLeafCells");
  MIntScratchSpace layers(56, AT_, "layers");
  const MInt noLeafCells = fvSolver().getAdjacentLeafCells_d2(cellId, 1, adjacentLeafCells, layers);

  for(MInt i = 0; i < noLeafCells; i++) {
    const MInt neighborId = adjacentLeafCells[i];
    const MFloat alpha = fvSolver().a_pvariable(neighborId, fvSolver().PV->A);
    if(searchUp && alpha < limit) continue;
    if(!searchUp && alpha > limit) continue;
    redistNeighbors.push_back(neighborId);
  }
  return redistNeighbors;
}

template class CouplerLbFvEEMultiphase<3, 27, maia::lb::LbSysEqnIncompressible<3, 27>, FvSysEqnEEGas<3>>;
template class CouplerLbFvEEMultiphase<3, 19, maia::lb::LbSysEqnIncompressible<3, 19>, FvSysEqnEEGas<3>>;
