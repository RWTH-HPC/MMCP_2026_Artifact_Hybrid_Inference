// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "lbinterfacedxqy.h"

#include <algorithm>
#include "IO/infoout.h"
#include "MEMORY/collector.h"
#include "UTIL/parallelfor.h"
#include "lbconstants.h"
#include "lbfunctions.h"
#include "lbinterfacecell.h"
#include "lbparentcell.h"
#include "lbsolver.h"

using namespace std;
using namespace lbconstants;

/**
 * \brief C'tor for the interface class.
 *
 * Storing some pointers to the solvers interface cell structures
 *
 * \param solver
 */
template <MInt nDim, MInt nDist, class SysEqn>
LbInterfaceDxQy<nDim, nDist, SysEqn>::LbInterfaceDxQy(LbSolver<nDim>* solver)
  : LbInterface<nDim>(solver)

{
  m_solver = static_cast<LbSolverDxQy<nDim, nDist, SysEqn>*>(solver);

  /*! \page propertyPage1
    \section interfaceMethod
    <code>MInt LbInterface::m_interfaceMethod</code>\n
    default = <code>"FILIPPOVA"</code>\n\n

    Selects interface method for locally refined meshes.\n
    - FILIPPOVA: Spatial and temporal interpolation of variables, rescaling of non-eq components\n
    - ROHDE: Volumetric formulation (no interpolation, no rescaling)\n\n

    Possible values are:\n
    <ul>
    <li><code>"FILIPPOVA"</code> (off)</li>
    <li><code>"ROHDE"</code> (on)</li>
    </ul>\n
    Keywords: <i>LATTICE_BOLTZMANN, NUMERICAL_SETUP, REFINEMENT</i>
  */
  const MInt solverId = m_solver->solverId();
  m_interfaceMethod = "FILIPPOVA";
  m_interfaceMethod = Context::getSolverProperty<MString>("interfaceMethod", solverId, AT_, &m_interfaceMethod);

  setInterfaceFunctions();

  m_adaptationInitMethod = "INIT_DUPUIS_FILIPPOVA";
  m_adaptationInitMethod =

      Context::getSolverProperty<MString>("interfaceMethod", solverId, AT_, &m_adaptationInitMethod);
  setAdaptationFunctions();
}

/**
 * \brief D'tor for the interface class
 *
 */
template <MInt nDim, MInt nDist, class SysEqn>
LbInterfaceDxQy<nDim, nDist, SysEqn>::~LbInterfaceDxQy() = default;

/**
 * \brief Setting function pointer to the chosen interface treatment method.
 *
 * \author Moritz Waldmann
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbInterfaceDxQy<nDim, nDist, SysEqn>::setInterfaceFunctions() {
  TRACE();

  switch(string2enum(m_interfaceMethod)) {
    // Default property value
    case FILIPPOVA: {
      if(m_isThermal) {
        fProlongation = &LbInterfaceDxQy::prolongationThermalDupuis;
        fRestriction = &LbInterfaceDxQy::restrictionThermalDupuis;
      } else if(m_solver->isCompressible()) {
        fProlongation = &LbInterfaceDxQy::prolongationDupuisCompressible;
        fRestriction = &LbInterfaceDxQy::restrictionDupuisCompressible;
      } else {
        fProlongation = &LbInterfaceDxQy::prolongationDupuis;
        fRestriction = &LbInterfaceDxQy::restrictionDupuis;
      }
      break;
    }
    case ROHDE: {
      if(m_isThermal) {
        fProlongation = &LbInterfaceDxQy::prolongationThermalRohde;
        fRestriction = &LbInterfaceDxQy::restrictionThermalRohde;
      } else {
        fProlongation = &LbInterfaceDxQy::prolongationRohde;
        fRestriction = &LbInterfaceDxQy::restrictionRohde;
      }
      break;
    }

    default: {
      stringstream errorMessage;
      errorMessage << " LbInterface::setInterfaceFunctions: Specified interface condition: " << m_interfaceMethod
                   << " does not exist. Exiting!";
      mTerm(1, AT_, errorMessage.str());
    }
  }
}

/**
 * \brief Performing the chosen prolongation method
 *
 * \author Moritz Waldmann
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbInterfaceDxQy<nDim, nDist, SysEqn>::prolongation() {
  (this->*fProlongation)();
}

/**
 * \brief Performing the chosen restriction method
 *
 * \author Moritz Waldmann
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbInterfaceDxQy<nDim, nDist, SysEqn>::restriction() {
  (this->*fRestriction)();
}

/** coarse to fine grid
 *
 * trilinear spatial interpolation + linear interpolation in time
 * and transformation of non-eq parts
 *
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbInterfaceDxQy<nDim, nDist, SysEqn>::prolongation10() {
  MFloat& tmp2 = m_static_prolongation10_tmp2;
  MFloat(&b)[2 * nDim] = m_static_prolongation10_b;
  MFloat(&c)[nDim * nDim] = m_static_prolongation10_c;
  MFloat& trace = m_static_prolongation10_trace;

  for(MInt level = 0; level < (m_solver->maxLevel() - m_solver->minLevel());
      level++) { // minLevel+level is the actual coarse resolution
                 // minLevel+level+1 is the actual fine resolution

    if((globalTimeStep) % IPOW2(m_solver->maxLevel() - m_solver->minLevel() - level)
       == 0) { // perform prolongation WITHOUT temporal interpolation in every second timestep on finest level,
               // in every fourth timestep on the next coarser level, and so on; starting from t = 1

      for(MInt id = 0; id < m_interfaceChildren[level]->size(); id++) {
        const MInt currentId = m_interfaceChildren[level]->a[id].m_cellId;
        if(currentId >= m_solver->noInternalCells()) continue;
        const MInt parentId = m_solver->c_parentId(m_interfaceChildren[level]->a[id].m_cellId);
        if(parentId >= m_solver->noInternalCells()) continue;

        const MFloat omega_F =
            2.0
            / (1.0 + 6.0 * m_solver->a_nu(currentId) * FFPOW2(m_solver->maxLevel() - m_solver->minLevel() - level - 1));
        // omega_C = 2.0 / ( 1.0 + 6.0 * m_solver->a_nu(currentId) * FFPOW2(m_solver->maxLevel() - m_solver->minLevel()
        // - level));

        //--------------------------------------------------------------------------------------------
        // Spatial interpolation and temporal extrapolation of macroscopic variables and stress tensor
        // f(t+1) = 2*f(t) - f(t-1)

        MFloat rho = 0.0;
        std::array<MFloat, nDim> u{};
        for(MInt m = 0; m < IPOW2(nDim); m++) {
          rho += 2.0 * m_interfaceChildren[level]->a[id].m_interpolationCoefficients[m]
                 * m_solver->a_variable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[m], PV->RHO);
          u[0] += 2.0 * m_interfaceChildren[level]->a[id].m_interpolationCoefficients[m]
                  * m_solver->a_variable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[m], PV->U);
          u[1] += 2.0 * m_interfaceChildren[level]->a[id].m_interpolationCoefficients[m]
                  * m_solver->a_variable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[m], PV->V);
          IF_CONSTEXPR(nDim == 3)
          u[2] += 2.0 * m_interfaceChildren[level]->a[id].m_interpolationCoefficients[m]
                  * m_solver->a_variable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[m], PV->W);
        }

        for(MInt m = 0; m < IPOW2(nDim); m++) {
          rho -= m_interfaceChildren[level]->a[id].m_interpolationCoefficients[m]
                 * m_solver->a_oldVariable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[m], PV->RHO);
          u[0] -= m_interfaceChildren[level]->a[id].m_interpolationCoefficients[m]
                  * m_solver->a_oldVariable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[m], PV->U);
          u[1] -= m_interfaceChildren[level]->a[id].m_interpolationCoefficients[m]
                  * m_solver->a_oldVariable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[m], PV->V);
          IF_CONSTEXPR(nDim == 3)
          u[2] -= m_interfaceChildren[level]->a[id].m_interpolationCoefficients[m]
                  * m_solver->a_oldVariable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[m], PV->W);
        }

        // compute strain rate tensor on fine grid: (du/dx)_fine = 1/2*(du/dx)_coarse
        IF_CONSTEXPR(nDim == 3) {
          for(MInt j = 0; j < 3; j++) {
            // du[j]/dx
            c[j + 0 * 3] =
                2.0 * rho * F1B2
                * ((m_interfaceChildren[level]->a[id].m_interpolationCoefficients[0]
                    + m_interfaceChildren[level]->a[id].m_interpolationCoefficients[1])
                       * (m_solver->a_variable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[1], j)
                          - m_solver->a_variable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[0], j))
                   + (m_interfaceChildren[level]->a[id].m_interpolationCoefficients[2]
                      + m_interfaceChildren[level]->a[id].m_interpolationCoefficients[3])
                         * (m_solver->a_variable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[3], j)
                            - m_solver->a_variable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[2], j))
                   + (m_interfaceChildren[level]->a[id].m_interpolationCoefficients[4]
                      + m_interfaceChildren[level]->a[id].m_interpolationCoefficients[5])
                         * (m_solver->a_variable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[5], j)
                            - m_solver->a_variable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[4], j))
                   + (m_interfaceChildren[level]->a[id].m_interpolationCoefficients[6]
                      + m_interfaceChildren[level]->a[id].m_interpolationCoefficients[7])
                         * (m_solver->a_variable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[7], j)
                            - m_solver->a_variable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[6], j)));
            // du[j]/dy
            c[j + 1 * 3] =
                2.0 * rho * F1B2
                * ((m_interfaceChildren[level]->a[id].m_interpolationCoefficients[2]
                    + m_interfaceChildren[level]->a[id].m_interpolationCoefficients[0])
                       * (m_solver->a_variable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[2], j)
                          - m_solver->a_variable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[0], j))
                   + (m_interfaceChildren[level]->a[id].m_interpolationCoefficients[3]
                      + m_interfaceChildren[level]->a[id].m_interpolationCoefficients[1])
                         * (m_solver->a_variable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[3], j)
                            - m_solver->a_variable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[1], j))
                   + (m_interfaceChildren[level]->a[id].m_interpolationCoefficients[6]
                      + m_interfaceChildren[level]->a[id].m_interpolationCoefficients[4])
                         * (m_solver->a_variable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[6], j)
                            - m_solver->a_variable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[4], j))
                   + (m_interfaceChildren[level]->a[id].m_interpolationCoefficients[7]
                      + m_interfaceChildren[level]->a[id].m_interpolationCoefficients[5])
                         * (m_solver->a_variable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[7], j)
                            - m_solver->a_variable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[5], j)));
            // du[j]/dz
            c[j + 2 * 3] =
                2.0 * rho * F1B2
                * ((m_interfaceChildren[level]->a[id].m_interpolationCoefficients[4]
                    + m_interfaceChildren[level]->a[id].m_interpolationCoefficients[0])
                       * (m_solver->a_variable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[4], j)
                          - m_solver->a_variable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[0], j))
                   + (m_interfaceChildren[level]->a[id].m_interpolationCoefficients[5]
                      + m_interfaceChildren[level]->a[id].m_interpolationCoefficients[1])
                         * (m_solver->a_variable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[5], j)
                            - m_solver->a_variable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[1], j))
                   + (m_interfaceChildren[level]->a[id].m_interpolationCoefficients[6]
                      + m_interfaceChildren[level]->a[id].m_interpolationCoefficients[2])
                         * (m_solver->a_variable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[6], j)
                            - m_solver->a_variable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[2], j))
                   + (m_interfaceChildren[level]->a[id].m_interpolationCoefficients[7]
                      + m_interfaceChildren[level]->a[id].m_interpolationCoefficients[3])
                         * (m_solver->a_variable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[7], j)
                            - m_solver->a_variable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[3], j)));
          }

          for(MInt j = 0; j < 3; j++) {
            // du[j]/dx
            c[j + 0 * 3] -=
                1.0 * rho * F1B2
                * ((m_interfaceChildren[level]->a[id].m_interpolationCoefficients[0]
                    + m_interfaceChildren[level]->a[id].m_interpolationCoefficients[1])
                       * (m_solver->a_oldVariable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[1], j)
                          - m_solver->a_oldVariable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[0], j))
                   + (m_interfaceChildren[level]->a[id].m_interpolationCoefficients[2]
                      + m_interfaceChildren[level]->a[id].m_interpolationCoefficients[3])
                         * (m_solver->a_oldVariable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[3], j)
                            - m_solver->a_oldVariable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[2], j))
                   + (m_interfaceChildren[level]->a[id].m_interpolationCoefficients[4]
                      + m_interfaceChildren[level]->a[id].m_interpolationCoefficients[5])
                         * (m_solver->a_oldVariable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[5], j)
                            - m_solver->a_oldVariable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[4], j))
                   + (m_interfaceChildren[level]->a[id].m_interpolationCoefficients[6]
                      + m_interfaceChildren[level]->a[id].m_interpolationCoefficients[7])
                         * (m_solver->a_oldVariable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[7], j)
                            - m_solver->a_oldVariable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[6],
                                                      j)));
            // du[j]/dy
            c[j + 1 * 3] -=
                1.0 * rho * F1B2
                * ((m_interfaceChildren[level]->a[id].m_interpolationCoefficients[2]
                    + m_interfaceChildren[level]->a[id].m_interpolationCoefficients[0])
                       * (m_solver->a_oldVariable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[2], j)
                          - m_solver->a_oldVariable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[0], j))
                   + (m_interfaceChildren[level]->a[id].m_interpolationCoefficients[3]
                      + m_interfaceChildren[level]->a[id].m_interpolationCoefficients[1])
                         * (m_solver->a_oldVariable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[3], j)
                            - m_solver->a_oldVariable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[1], j))
                   + (m_interfaceChildren[level]->a[id].m_interpolationCoefficients[6]
                      + m_interfaceChildren[level]->a[id].m_interpolationCoefficients[4])
                         * (m_solver->a_oldVariable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[6], j)
                            - m_solver->a_oldVariable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[4], j))
                   + (m_interfaceChildren[level]->a[id].m_interpolationCoefficients[7]
                      + m_interfaceChildren[level]->a[id].m_interpolationCoefficients[5])
                         * (m_solver->a_oldVariable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[7], j)
                            - m_solver->a_oldVariable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[5],
                                                      j)));
            // du[j]/dz
            c[j + 2 * 3] -=
                1.0 * rho * F1B2
                * ((m_interfaceChildren[level]->a[id].m_interpolationCoefficients[4]
                    + m_interfaceChildren[level]->a[id].m_interpolationCoefficients[0])
                       * (m_solver->a_oldVariable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[4], j)
                          - m_solver->a_oldVariable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[0], j))
                   + (m_interfaceChildren[level]->a[id].m_interpolationCoefficients[5]
                      + m_interfaceChildren[level]->a[id].m_interpolationCoefficients[1])
                         * (m_solver->a_oldVariable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[5], j)
                            - m_solver->a_oldVariable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[1], j))
                   + (m_interfaceChildren[level]->a[id].m_interpolationCoefficients[6]
                      + m_interfaceChildren[level]->a[id].m_interpolationCoefficients[2])
                         * (m_solver->a_oldVariable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[6], j)
                            - m_solver->a_oldVariable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[2], j))
                   + (m_interfaceChildren[level]->a[id].m_interpolationCoefficients[7]
                      + m_interfaceChildren[level]->a[id].m_interpolationCoefficients[3])
                         * (m_solver->a_oldVariable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[7], j)
                            - m_solver->a_oldVariable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[3],
                                                      j)));
          }
          trace = c[0] + c[4] + c[8];
        }
        else {
          cout << "Strain rate tensor not implemented for 2D yet" << endl;
          trace = F0;
        }

        //--------------------------------------------------------------------------------------------


        //  	      m_solver->a_variable(currentId, PV->RHO) = rho;
        //  	      m_solver->a_variable(currentId, PV->U) = u[0];
        //  	      m_solver->a_variable(currentId, PV->V) = u[1];
        //  	      m_solver->a_variable(currentId, PV->W) = u[2];

        tmp2 = F0;
        for(MInt d = 0; d < nDim; d++) {
          tmp2 += (u[d] * u[d]);
          b[2 * d] = -u[d];
          b[2 * d + 1] = u[d];
        }

        // Calculate all equilibrium distributions
        std::array<MFloat, nDist> eqDist;
        eqDist = m_solver->getEqDists(rho, tmp2, u.data());

        // interpolate and transform
        for(MInt dist = 0; dist < nDist - 1; dist++) {
          // only the missing distributions
          if(m_solver->a_hasNeighbor(currentId, Ld::oppositeDist(dist)) == 0) {
            const MFloat tp = Ld::tp(Ld::distType(dist));

            // put together eq-part and rescaled non-eq part; scale factor = (omega_C/omega_F) * F1B2

            m_solver->a_oldDistribution(currentId, dist) = eqDist[dist] + tp * trace / omega_F;

            for(MInt k = 0; k < nDim; k++) {
              for(MInt l = 0; l < nDim; l++) {
                m_solver->a_oldDistribution(currentId, dist) -=
                    tp * F1BCSsq * (Ld::idFld(dist, k) - 1) * (Ld::idFld(dist, l) - 1) * c[l + nDim * k] / omega_F;
              }
            }
          }
        }
      } // end of loop over all interface cells

    }
    // If NOT in sync with parent cells apply linear time interpolation
    else {
      if((globalTimeStep - 1) % IPOW2(m_solver->maxLevel() - m_solver->minLevel() - level)
         == 0) { // perform prolongation WITH temporal interpolation in every second timestep on finest level,
                 // in every fourth timestep on the next coarser level, and so on; starting from t = 2

        for(MInt id = 0; id < m_interfaceChildren[level]->size(); id++) {
          const MInt currentId = m_interfaceChildren[level]->a[id].m_cellId;
          if(currentId >= m_solver->noInternalCells()) continue;
          const MInt parentId = m_solver->c_parentId(m_interfaceChildren[level]->a[id].m_cellId);
          if(parentId >= m_solver->noInternalCells()) continue;

          const MFloat omega_F =
              2.0
              / (1.0
                 + 6.0 * m_solver->a_nu(currentId) * FFPOW2(m_solver->maxLevel() - m_solver->minLevel() - level - 1));

          //--------------------------------------------------------------------------------------------
          // Spatial interpolation and temporal extrapolation of macroscopic variables and stress tensor
          // f(t+1/2) = 1.5*f(t) - 0.5*f(t-1)

          MFloat rho = 0;
          std::array<MFloat, nDim> u{};

          for(MInt m = 0; m < IPOW2(nDim); m++) {
            rho += 1.5 * m_interfaceChildren[level]->a[id].m_interpolationCoefficients[m]
                   * m_solver->a_variable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[m], PV->RHO);
            u[0] += 1.5 * m_interfaceChildren[level]->a[id].m_interpolationCoefficients[m]
                    * m_solver->a_variable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[m], PV->U);
            u[1] += 1.5 * m_interfaceChildren[level]->a[id].m_interpolationCoefficients[m]
                    * m_solver->a_variable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[m], PV->V);
            IF_CONSTEXPR(nDim == 3)
            u[2] += 1.5 * m_interfaceChildren[level]->a[id].m_interpolationCoefficients[m]
                    * m_solver->a_variable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[m], PV->W);
          }

          for(MInt m = 0; m < IPOW2(nDim); m++) {
            rho -= 0.5 * m_interfaceChildren[level]->a[id].m_interpolationCoefficients[m]
                   * m_solver->a_oldVariable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[m], PV->RHO);
            u[0] -= 0.5 * m_interfaceChildren[level]->a[id].m_interpolationCoefficients[m]
                    * m_solver->a_oldVariable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[m], PV->U);
            u[1] -= 0.5 * m_interfaceChildren[level]->a[id].m_interpolationCoefficients[m]
                    * m_solver->a_oldVariable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[m], PV->V);
            IF_CONSTEXPR(nDim == 3)
            u[2] -= 0.5 * m_interfaceChildren[level]->a[id].m_interpolationCoefficients[m]
                    * m_solver->a_oldVariable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[m], PV->W);
          }

          // compute strain rate tensor on fine grid: (du/dx)_fine = 1/2*(du/dx)_coarse
          IF_CONSTEXPR(nDim == 3) {
            for(MInt j = 0; j < 3; j++) {
              // du[j]/dx
              c[j + 0 * 3] =
                  1.5 * rho * F1B2
                  * ((m_interfaceChildren[level]->a[id].m_interpolationCoefficients[0]
                      + m_interfaceChildren[level]->a[id].m_interpolationCoefficients[1])
                         * (m_solver->a_variable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[1], j)
                            - m_solver->a_variable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[0], j))
                     + (m_interfaceChildren[level]->a[id].m_interpolationCoefficients[2]
                        + m_interfaceChildren[level]->a[id].m_interpolationCoefficients[3])
                           * (m_solver->a_variable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[3], j)
                              - m_solver->a_variable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[2], j))
                     + (m_interfaceChildren[level]->a[id].m_interpolationCoefficients[4]
                        + m_interfaceChildren[level]->a[id].m_interpolationCoefficients[5])
                           * (m_solver->a_variable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[5], j)
                              - m_solver->a_variable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[4], j))
                     + (m_interfaceChildren[level]->a[id].m_interpolationCoefficients[6]
                        + m_interfaceChildren[level]->a[id].m_interpolationCoefficients[7])
                           * (m_solver->a_variable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[7], j)
                              - m_solver->a_variable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[6],
                                                     j)));
              // du[j]/dy
              c[j + 1 * 3] =
                  1.5 * rho * F1B2
                  * ((m_interfaceChildren[level]->a[id].m_interpolationCoefficients[2]
                      + m_interfaceChildren[level]->a[id].m_interpolationCoefficients[0])
                         * (m_solver->a_variable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[2], j)
                            - m_solver->a_variable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[0], j))
                     + (m_interfaceChildren[level]->a[id].m_interpolationCoefficients[3]
                        + m_interfaceChildren[level]->a[id].m_interpolationCoefficients[1])
                           * (m_solver->a_variable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[3], j)
                              - m_solver->a_variable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[1], j))
                     + (m_interfaceChildren[level]->a[id].m_interpolationCoefficients[6]
                        + m_interfaceChildren[level]->a[id].m_interpolationCoefficients[4])
                           * (m_solver->a_variable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[6], j)
                              - m_solver->a_variable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[4], j))
                     + (m_interfaceChildren[level]->a[id].m_interpolationCoefficients[7]
                        + m_interfaceChildren[level]->a[id].m_interpolationCoefficients[5])
                           * (m_solver->a_variable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[7], j)
                              - m_solver->a_variable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[5],
                                                     j)));
              // du[j]/dz
              c[j + 2 * 3] =
                  1.5 * rho * F1B2
                  * ((m_interfaceChildren[level]->a[id].m_interpolationCoefficients[4]
                      + m_interfaceChildren[level]->a[id].m_interpolationCoefficients[0])
                         * (m_solver->a_variable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[4], j)
                            - m_solver->a_variable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[0], j))
                     + (m_interfaceChildren[level]->a[id].m_interpolationCoefficients[5]
                        + m_interfaceChildren[level]->a[id].m_interpolationCoefficients[1])
                           * (m_solver->a_variable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[5], j)
                              - m_solver->a_variable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[1], j))
                     + (m_interfaceChildren[level]->a[id].m_interpolationCoefficients[6]
                        + m_interfaceChildren[level]->a[id].m_interpolationCoefficients[2])
                           * (m_solver->a_variable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[6], j)
                              - m_solver->a_variable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[2], j))
                     + (m_interfaceChildren[level]->a[id].m_interpolationCoefficients[7]
                        + m_interfaceChildren[level]->a[id].m_interpolationCoefficients[3])
                           * (m_solver->a_variable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[7], j)
                              - m_solver->a_variable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[3],
                                                     j)));
            }

            for(MInt j = 0; j < 3; j++) {
              // du[j]/dx
              c[j + 0 * 3] -=
                  F1B2 * rho * F1B2
                  * ((m_interfaceChildren[level]->a[id].m_interpolationCoefficients[0]
                      + m_interfaceChildren[level]->a[id].m_interpolationCoefficients[1])
                         * (m_solver->a_oldVariable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[1], j)
                            - m_solver->a_oldVariable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[0], j))
                     + (m_interfaceChildren[level]->a[id].m_interpolationCoefficients[2]
                        + m_interfaceChildren[level]->a[id].m_interpolationCoefficients[3])
                           * (m_solver->a_oldVariable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[3], j)
                              - m_solver->a_oldVariable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[2],
                                                        j))
                     + (m_interfaceChildren[level]->a[id].m_interpolationCoefficients[4]
                        + m_interfaceChildren[level]->a[id].m_interpolationCoefficients[5])
                           * (m_solver->a_oldVariable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[5], j)
                              - m_solver->a_oldVariable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[4],
                                                        j))
                     + (m_interfaceChildren[level]->a[id].m_interpolationCoefficients[6]
                        + m_interfaceChildren[level]->a[id].m_interpolationCoefficients[7])
                           * (m_solver->a_oldVariable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[7], j)
                              - m_solver->a_oldVariable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[6],
                                                        j)));
              // du[j]/dy
              c[j + 1 * 3] -=
                  F1B2 * rho * F1B2
                  * ((m_interfaceChildren[level]->a[id].m_interpolationCoefficients[2]
                      + m_interfaceChildren[level]->a[id].m_interpolationCoefficients[0])
                         * (m_solver->a_oldVariable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[2], j)
                            - m_solver->a_oldVariable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[0], j))
                     + (m_interfaceChildren[level]->a[id].m_interpolationCoefficients[3]
                        + m_interfaceChildren[level]->a[id].m_interpolationCoefficients[1])
                           * (m_solver->a_oldVariable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[3], j)
                              - m_solver->a_oldVariable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[1],
                                                        j))
                     + (m_interfaceChildren[level]->a[id].m_interpolationCoefficients[6]
                        + m_interfaceChildren[level]->a[id].m_interpolationCoefficients[4])
                           * (m_solver->a_oldVariable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[6], j)
                              - m_solver->a_oldVariable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[4],
                                                        j))
                     + (m_interfaceChildren[level]->a[id].m_interpolationCoefficients[7]
                        + m_interfaceChildren[level]->a[id].m_interpolationCoefficients[5])
                           * (m_solver->a_oldVariable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[7], j)
                              - m_solver->a_oldVariable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[5],
                                                        j)));
              // du[j]/dz
              c[j + 2 * 3] -=
                  F1B2 * rho * F1B2
                  * ((m_interfaceChildren[level]->a[id].m_interpolationCoefficients[4]
                      + m_interfaceChildren[level]->a[id].m_interpolationCoefficients[0])
                         * (m_solver->a_oldVariable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[4], j)
                            - m_solver->a_oldVariable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[0], j))
                     + (m_interfaceChildren[level]->a[id].m_interpolationCoefficients[5]
                        + m_interfaceChildren[level]->a[id].m_interpolationCoefficients[1])
                           * (m_solver->a_oldVariable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[5], j)
                              - m_solver->a_oldVariable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[1],
                                                        j))
                     + (m_interfaceChildren[level]->a[id].m_interpolationCoefficients[6]
                        + m_interfaceChildren[level]->a[id].m_interpolationCoefficients[2])
                           * (m_solver->a_oldVariable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[6], j)
                              - m_solver->a_oldVariable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[2],
                                                        j))
                     + (m_interfaceChildren[level]->a[id].m_interpolationCoefficients[7]
                        + m_interfaceChildren[level]->a[id].m_interpolationCoefficients[3])
                           * (m_solver->a_oldVariable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[7], j)
                              - m_solver->a_oldVariable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[3],
                                                        j)));
            }
            trace = c[0] + c[4] + c[8];
          }
          else {
            cout << "Strain rate tensor not implemented for 2D yet" << endl;
            trace = F0;
          }

          //--------------------------------------------------------------------------------------------

          tmp2 = F0;
          for(MInt d = 0; d < nDim; d++) {
            tmp2 += (u[d] * u[d]);
            b[2 * d] = -u[d];
            b[2 * d + 1] = u[d];
          }

          // Calculation of new distributions
          std::array<MFloat, nDist> eqDist;
          eqDist = m_solver->getEqDists(rho, tmp2, u.data());

          // create new incoming distributions
          for(MInt dist = 0; dist < nDist - 1; dist++) {
            // only the missing distributions
            if(m_solver->a_hasNeighbor(currentId, Ld::oppositeDist(dist)) == 0) {
              const MFloat tp = Ld::tp(Ld::distType(dist));

              // put together eq-part and rescaled non-eq part; scale factor = (omega_C/omega_F) * F1B2

              m_solver->a_oldDistribution(currentId, dist) = eqDist[dist] + tp * trace / omega_F;

              for(MInt k = 0; k < nDim; k++) {
                for(MInt l = 0; l < nDim; l++) {
                  m_solver->a_oldDistribution(currentId, dist) -=
                      tp * F1BCSsq * (Ld::idFld(dist, k) - 1) * (Ld::idFld(dist, l) - 1) * c[l + nDim * k] / omega_F;
                }
              }
            }
          }
        } // end of loop over all interface cells
      }
    }
  }
}

/** fine to coarse grid
 *
 * average all child cell values
 * and add appropriate non-eq parts
 *
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbInterfaceDxQy<nDim, nDist, SysEqn>::restriction10() {
  MFloat& tmp2 = m_static_restriction10_tmp2;
  MFloat(&b)[2 * nDim] = m_static_restriction10_b;
  MFloat(&c)[nDim * nDim] = m_static_restriction10_c;
  MFloat& trace = m_static_restriction10_trace;

  for(MInt level = 0; level < (m_solver->maxLevel() - m_solver->minLevel());
      level++) { // level is the actual coarse resolution
                 // level+1 is the actual fine resolution

    if((globalTimeStep) % IPOW2(m_solver->maxLevel() - m_solver->minLevel() - level)
       == 0) { // perform restriction from finest level in every second timestep,
               // on the next coarser level every fourth timestep, and so on; starting from t = 1

      for(MInt id = 0; id < m_interfaceParents[level]->size(); id++) {
        const MInt currentId = m_interfaceParents[level]->a[id].m_cellId;
        if(currentId >= m_solver->noInternalCells()) continue;

        // omega_F = 2.0 / ( 1.0 + 6.0 * m_solver->a_nu(currentId) * FFPOW2(m_solver->maxLevel() - m_solver->minLevel()
        // - level - 1));
        const MFloat omega_C =
            2.0 / (1.0 + 6.0 * m_solver->a_nu(currentId) * FFPOW2(m_solver->maxLevel() - m_solver->minLevel() - level));

        MFloat rho = 0;
        std::array<MFloat, nDim> u{};

        //------------------------------------------------------------------------------------
        // temporal extrapolation of the mean macroscopic variables and the mean stress tensor of all child cells
        // f(t+1) = 2*f(t) - f(t-1)

        for(MInt m = 0; m < IPOW2(nDim); m++) {
          rho += 2.0 * FFPOW2(nDim) * m_solver->a_variable(m_solver->c_childId(currentId, m), PV->RHO);
          u[0] += 2.0 * FFPOW2(nDim) * m_solver->a_variable(m_solver->c_childId(currentId, m), PV->U);
          u[1] += 2.0 * FFPOW2(nDim) * m_solver->a_variable(m_solver->c_childId(currentId, m), PV->V);
          IF_CONSTEXPR(nDim == 3)
          u[2] += 2.0 * FFPOW2(nDim) * m_solver->a_variable(m_solver->c_childId(currentId, m), PV->W);
        }
        for(MInt m = 0; m < IPOW2(nDim); m++) {
          rho -= FFPOW2(nDim) * m_solver->a_oldVariable(m_solver->c_childId(currentId, m), PV->RHO);
          u[0] -= FFPOW2(nDim) * m_solver->a_oldVariable(m_solver->c_childId(currentId, m), PV->U);
          u[1] -= FFPOW2(nDim) * m_solver->a_oldVariable(m_solver->c_childId(currentId, m), PV->V);
          IF_CONSTEXPR(nDim == 3)
          u[2] -= FFPOW2(nDim) * m_solver->a_oldVariable(m_solver->c_childId(currentId, m), PV->W);
        }

        // compute strain rate tensor on coarse grid: (du/dx)_coarse = 2*(du/dx)_fine
        IF_CONSTEXPR(nDim == 3) {
          for(MInt j = 0; j < 3; j++) {
            // du[j]/dx
            c[j] = 2.0 * F1B4 * rho * 2.0
                   * ((m_solver->a_variable(m_solver->c_childId(currentId, 1), j)
                       - m_solver->a_variable(m_solver->c_childId(currentId, 0), j))
                      + (m_solver->a_variable(m_solver->c_childId(currentId, 3), j)
                         - m_solver->a_variable(m_solver->c_childId(currentId, 2), j))
                      + (m_solver->a_variable(m_solver->c_childId(currentId, 5), j)
                         - m_solver->a_variable(m_solver->c_childId(currentId, 4), j))
                      + (m_solver->a_variable(m_solver->c_childId(currentId, 7), j)
                         - m_solver->a_variable(m_solver->c_childId(currentId, 6), j)));
            // du[j]/dy
            c[j + 3] = 2.0 * F1B4 * rho * 2.0
                       * ((m_solver->a_variable(m_solver->c_childId(currentId, 2), j)
                           - m_solver->a_variable(m_solver->c_childId(currentId, 0), j))
                          + (m_solver->a_variable(m_solver->c_childId(currentId, 3), j)
                             - m_solver->a_variable(m_solver->c_childId(currentId, 1), j))
                          + (m_solver->a_variable(m_solver->c_childId(currentId, 6), j)
                             - m_solver->a_variable(m_solver->c_childId(currentId, 4), j))
                          + (m_solver->a_variable(m_solver->c_childId(currentId, 7), j)
                             - m_solver->a_variable(m_solver->c_childId(currentId, 5), j)));
            // du[j]/dz
            c[j + 6] = 2.0 * F1B4 * rho * 2.0
                       * ((m_solver->a_variable(m_solver->c_childId(currentId, 4), j)
                           - m_solver->a_variable(m_solver->c_childId(currentId, 0), j))
                          + (m_solver->a_variable(m_solver->c_childId(currentId, 5), j)
                             - m_solver->a_variable(m_solver->c_childId(currentId, 1), j))
                          + (m_solver->a_variable(m_solver->c_childId(currentId, 6), j)
                             - m_solver->a_variable(m_solver->c_childId(currentId, 2), j))
                          + (m_solver->a_variable(m_solver->c_childId(currentId, 7), j)
                             - m_solver->a_variable(m_solver->c_childId(currentId, 3), j)));
          }

          for(MInt j = 0; j < 3; j++) {
            // du[j]/dx
            c[j] -= 1.0 * F1B4 * rho * 2.0
                    * ((m_solver->a_oldVariable(m_solver->c_childId(currentId, 1), j)
                        - m_solver->a_oldVariable(m_solver->c_childId(currentId, 0), j))
                       + (m_solver->a_oldVariable(m_solver->c_childId(currentId, 3), j)
                          - m_solver->a_oldVariable(m_solver->c_childId(currentId, 2), j))
                       + (m_solver->a_oldVariable(m_solver->c_childId(currentId, 5), j)
                          - m_solver->a_oldVariable(m_solver->c_childId(currentId, 4), j))
                       + (m_solver->a_oldVariable(m_solver->c_childId(currentId, 7), j)
                          - m_solver->a_oldVariable(m_solver->c_childId(currentId, 6), j)));
            // du[j]/dy
            c[j + 3] -= 1.0 * F1B4 * rho * 2.0
                        * ((m_solver->a_oldVariable(m_solver->c_childId(currentId, 2), j)
                            - m_solver->a_oldVariable(m_solver->c_childId(currentId, 0), j))
                           + (m_solver->a_oldVariable(m_solver->c_childId(currentId, 3), j)
                              - m_solver->a_oldVariable(m_solver->c_childId(currentId, 1), j))
                           + (m_solver->a_oldVariable(m_solver->c_childId(currentId, 6), j)
                              - m_solver->a_oldVariable(m_solver->c_childId(currentId, 4), j))
                           + (m_solver->a_oldVariable(m_solver->c_childId(currentId, 7), j)
                              - m_solver->a_oldVariable(m_solver->c_childId(currentId, 5), j)));
            // du[j]/dz
            c[j + 6] -= 1.0 * F1B4 * rho * 2.0
                        * ((m_solver->a_oldVariable(m_solver->c_childId(currentId, 4), j)
                            - m_solver->a_oldVariable(m_solver->c_childId(currentId, 0), j))
                           + (m_solver->a_oldVariable(m_solver->c_childId(currentId, 5), j)
                              - m_solver->a_oldVariable(m_solver->c_childId(currentId, 1), j))
                           + (m_solver->a_oldVariable(m_solver->c_childId(currentId, 6), j)
                              - m_solver->a_oldVariable(m_solver->c_childId(currentId, 2), j))
                           + (m_solver->a_oldVariable(m_solver->c_childId(currentId, 7), j)
                              - m_solver->a_oldVariable(m_solver->c_childId(currentId, 3), j)));
          }
          trace = c[0] + c[4] + c[8];
        }
        else {
          cout << "Strain rate tensor not implemented for 2D yet" << endl;
          trace = F0;
        }
        //------------------------------------------------------------------------------------


        // set parent values to averaged child values for subsequent prolongation
        m_solver->a_variable(currentId, PV->RHO) = rho;
        m_solver->a_variable(currentId, PV->U) = u[0];
        m_solver->a_variable(currentId, PV->V) = u[1];
        IF_CONSTEXPR(nDim == 3) m_solver->a_variable(currentId, PV->W) = u[2];

        tmp2 = F0;
        for(MInt d = 0; d < nDim; d++) {
          tmp2 += (u[d] * u[d]);
          b[2 * d] = -u[d];
          b[2 * d + 1] = u[d];
        }

        std::array<MFloat, nDist> eqDist;
        eqDist = m_solver->getEqDists(rho, tmp2, u.data());

        // create new incoming distributions
        for(MInt dist = 0; dist < nDist - 1; dist++) {
          // only the missing distributions
          if(!m_solver->a_hasNeighbor(currentId, Ld::oppositeDist(dist))
             || (!m_solver->a_isInterfaceParent(m_solver->c_neighborId(currentId, Ld::oppositeDist(dist)))
                 && !m_solver->c_isLeafCell(m_solver->c_neighborId(currentId, Ld::oppositeDist(dist))))) {
            const MFloat tp = Ld::tp(Ld::distType(dist));

            // put together eq-part and rescaled non-eq part; scale factor = (omega_C/omega_F) * F1B2

            m_solver->a_oldDistribution(currentId, dist) = eqDist[dist] + tp * trace / omega_C;

            for(MInt k = 0; k < nDim; k++) {
              for(MInt l = 0; l < nDim; l++) {
                m_solver->a_oldDistribution(currentId, dist) -=
                    tp * F1BCSsq * (Ld::idFld(dist, k) - 1) * (Ld::idFld(dist, l) - 1) * c[l + nDim * k] / omega_C;
              }
            }
          }
        }
      }
    }
  }
}

/** \brief Coarse to fine
 *
 * trilinear spatial interpolation + linear interpolation in time
 * and transformation of non-eq parts according to Dupuis et al.
 *
 */
template <MInt nDim, MInt nDist, class SysEqn>
template <MBool compressible>
void LbInterfaceDxQy<nDim, nDist, SysEqn>::prolongationDupuis_() {
#ifdef WAR_NVHPC_PSTL
  const MInt noInternalCells = m_solver->noInternalCells();
#endif

  for(MInt level = 0; level < (m_solver->maxLevel() - m_solver->minLevel());
      level++) { // minLevel+level is the actual coarse resolution
                 // minLevel+level+1 is the actual fine resolution

    if((globalTimeStep) % IPOW2(m_solver->maxLevel() - m_solver->minLevel() - level - 1)
       == 0) { // perform prolongation WITH temporal interpolation in every second timestep on finest level,
               // in every fourth timestep on the next coarser level, and so on; starting from t = 2

      //-----------------------------
      // if in sync with parent:
      if((globalTimeStep) % IPOW2(m_solver->maxLevel() - m_solver->minLevel() - level)
         == 0) { // perform prolongation WITHOUT temporal interpolation in every second timestep on finest level,
                 // in every fourth timestep on the next coarser level, and so on; starting from t = 1

        maia::parallelFor<true>(0, m_interfaceChildren[level]->size(), [=](MInt id) {
          const MInt currentId = m_interfaceChildren[level]->a[id].m_cellId;
#ifdef WAR_NVHPC_PSTL
          if(currentId >= noInternalCells) return;
#else
          if(currentId >= m_solver->noInternalCells()) return;
#endif

          const MFloat omega_F =
              2.0
              / (1.0
                 + 6.0 * m_solver->a_nu(currentId) * FFPOW2(m_solver->maxLevel() - m_solver->minLevel() - level - 1));
          const MFloat omega_C =
              2.0
              / (1.0 + 6.0 * m_solver->a_nu(currentId) * FFPOW2(m_solver->maxLevel() - m_solver->minLevel() - level));

          // Spatial interpolation and temporal extrapolation of macroscopic variables
          // f(t+1) = 2*f(t) - f(t-1)
          MFloat rho = 0.0;
          std::array<MFloat, nDim> u{};

          for(MInt m = 0; m < IPOW2(nDim); m++) {
            rho += 2.0 * m_interfaceChildren[level]->a[id].m_interpolationCoefficients[m]
                   * m_solver->a_variable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[m], PV->RHO);
            for(MInt d = 0; d < nDim; d++) {
              u[d] += 2.0 * m_interfaceChildren[level]->a[id].m_interpolationCoefficients[m]
                      * m_solver->a_variable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[m], d);
            }
          }

          for(MInt m = 0; m < IPOW2(nDim); m++) {
            rho -= m_interfaceChildren[level]->a[id].m_interpolationCoefficients[m]
                   * m_solver->a_oldVariable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[m], PV->RHO);
            for(MInt d = 0; d < nDim; d++) {
              u[d] -= m_interfaceChildren[level]->a[id].m_interpolationCoefficients[m]
                      * m_solver->a_oldVariable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[m], d);
            }
          }

          // Calculate all equilibrium distributions
          std::array<MFloat, nDist> eqDist;
          eqDist = m_solver->getEqDists(rho, u.data());

          // interpolate and transform
          for(MInt dist = 0; dist < nDist - 1; dist++) {
            // only the missing distributions
#ifdef WAR_NVHPC_PSTL
            MInt oppositeDist = m_solver->m_oppositeDist[dist];
#else
            MInt oppositeDist = Ld::oppositeDist(dist);
#endif
            if(m_solver->a_hasNeighbor(currentId, dist) == 0) {
              MFloat tmpDist1 = 0.0;

              for(MInt m = 0; m < IPOW2(nDim); m++) {
                tmpDist1 += m_interfaceChildren[level]->a[id].m_interpolationCoefficients[m]
                            * m_solver->a_oldDistribution(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[m],
                                                          oppositeDist);
              }

              m_solver->a_oldDistribution(currentId, oppositeDist) =
                  (tmpDist1 - eqDist[oppositeDist]) * (omega_C / omega_F) * F1B2 + eqDist[oppositeDist];
            }
          }
        });
      }

      //---------------------------
      // if not in sync with parent: apply temporal interpolation
      else {
        maia::parallelFor<true>(0, m_interfaceChildren[level]->size(), [=](MInt id) {
          const MInt currentId = m_interfaceChildren[level]->a[id].m_cellId;
#ifdef WAR_NVHPC_PSTL
          if(currentId >= noInternalCells) return;
#else
          if(currentId >= m_solver->noInternalCells()) return;
#endif

          const MFloat omega_F =
              2.0
              / (1.0
                 + 6.0 * m_solver->a_nu(currentId) * FFPOW2(m_solver->maxLevel() - m_solver->minLevel() - level - 1));
          const MFloat omega_C =
              2.0
              / (1.0 + 6.0 * m_solver->a_nu(currentId) * FFPOW2(m_solver->maxLevel() - m_solver->minLevel() - level));

          // spatial and temporal interpolation of variables
          // f(t+1/2) = 1.5*f(t) - 0.5*f(t-1)
          MFloat rho = 0;
          std::array<MFloat, nDim> u{};

          for(MInt m = 0; m < IPOW2(nDim); m++) {
            rho += 1.5 * m_interfaceChildren[level]->a[id].m_interpolationCoefficients[m]
                   * m_solver->a_variable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[m], PV->RHO);
            for(MInt d = 0; d < nDim; d++) {
              u[d] += 1.5 * m_interfaceChildren[level]->a[id].m_interpolationCoefficients[m]
                      * m_solver->a_variable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[m], d);
            }
          }

          for(MInt m = 0; m < IPOW2(nDim); m++) {
            rho -= 0.5 * m_interfaceChildren[level]->a[id].m_interpolationCoefficients[m]
                   * m_solver->a_oldVariable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[m], PV->RHO);
            for(MInt d = 0; d < nDim; d++) {
              u[d] -= 0.5 * m_interfaceChildren[level]->a[id].m_interpolationCoefficients[m]
                      * m_solver->a_oldVariable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[m], d);
            }
          }

          // Calculation of new distributions
          std::array<MFloat, nDist> eqDist;
          eqDist = m_solver->getEqDists(rho, u.data());
          std::array<MFloat, nDist> forcing{};
          getCellForcing(forcing, currentId);

          // Transform and interpolate
          for(MInt dist = 0; dist < nDist - 1; dist++) {
            // only for missing distributions
#ifdef WAR_NVHPC_PSTL
            MInt oppositeDist = m_solver->m_oppositeDist[dist];
#else
            MInt oppositeDist = Ld::oppositeDist(dist);
#endif
            if(m_solver->a_hasNeighbor(currentId, dist) == 0) {
              MFloat tmpDist1 = 0.0;
              MFloat tmpDist2 = 0.0;

              for(MInt m = 0; m < IPOW2(nDim); m++) {
                // add the weighted parts to the non-equilibrium distributions using linear temporal interpolation
                // f(t+1/2) = 1/2*f(t) + 1/2*f(t+1)
                const MFloat tmpDist = m_interfaceChildren[level]->a[id].m_interpolationCoefficients[m]
                                       * m_solver->a_oldDistribution(
                                           m_interfaceChildren[level]->a[id].m_interpolationNeighbors[m], oppositeDist);
                if(m_solver->a_hasNeighbor(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[m], dist)) {
                  // the propagation on the coarse grid has not taken place yet,
                  // thus the incoming distribution for the coarse cell can be interpolated as follows 0.5*(tmpDist1 +
                  // tmpDist2)
                  const MInt neighborId =
                      m_solver->c_neighborId(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[m], dist);
                  tmpDist1 += m_interfaceChildren[level]->a[id].m_interpolationCoefficients[m]
                              * m_solver->a_distribution(neighborId, oppositeDist);
                  tmpDist2 += tmpDist;

                  // add forcing terms for coarse grid before transformation
                  // for tmpDist2 the forcing terms have already been updated
                  tmpDist1 += m_interfaceChildren[level]->a[id].m_interpolationCoefficients[m]
                              * FPOW2(m_solver->maxLevel() - m_solver->minLevel() - level) * forcing[oppositeDist];
                } else {
                  // if the interpolation neighbor has no neighbor in actual direction, no temporal interpolation can be
                  // applied this may only occur at a non-periodic boundary
                  tmpDist2 += tmpDist;

                  // forcing term is already included
                  tmpDist1 += tmpDist;
                }
              }

              m_solver->a_oldDistribution(currentId, oppositeDist) =
                  (0.5 * (tmpDist1 + tmpDist2) - eqDist[oppositeDist]) * (omega_C / omega_F) * F1B2
                  + eqDist[oppositeDist];
            }
          }
        });
      } // end of else
    }
  }
}

template <MInt nDim, MInt nDist, class SysEqn>
void LbInterfaceDxQy<nDim, nDist, SysEqn>::prolongationDupuis() {
  prolongationDupuis_<false>();
}

template <MInt nDim, MInt nDist, class SysEqn>
void LbInterfaceDxQy<nDim, nDist, SysEqn>::prolongationDupuisCompressible() {
  prolongationDupuis_<true>();
}

/** \brief Coarse to fine grid for thermal LB
 *
 * \date 31.07.2012
 * \author Andreas Lintermann
 *
 * \todo labels:LB This is currently not working properly, please use ROHDE method.
 *
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbInterfaceDxQy<nDim, nDist, SysEqn>::prolongationThermalDupuis() {
  MFloat* ptr_m_interpolationCoefficients;
  const MFloat factor = (m_solver->m_innerEnergy) ? ((nDim == 2) ? 3.0 : 18.0 / 5.0) : 6.0;

  for(MInt level = 0; level < (m_solver->maxLevel() - m_solver->minLevel());
      level++) { // minLevel+level is the actual coarse resolution
                 // minLevel+level+1 is the actual fine resolution

    if((globalTimeStep) % IPOW2(m_solver->maxLevel() - m_solver->minLevel() - level - 1)
       == 0) { // perform prolongation WITH temporal interpolation in every second timestep on finest level,
               // in every fourth timestep on the next coarser level, and so on; starting from t = 2

      //-----------------------------
      // if in sync with parent:
      if((globalTimeStep) % IPOW2(m_solver->maxLevel() - m_solver->minLevel() - level)
         == 0) { // perform prolongation WITHOUT temporal interpolation in every second timestep on finest level,
                 // in every fourth timestep on the next coarser level, and so on; starting from t = 1

        for(MInt id = 0; id < m_interfaceChildren[level]->size(); id++) {
          const MInt currentId = m_interfaceChildren[level]->a[id].m_cellId;
          if(currentId >= m_solver->noInternalCells()) continue;

          ptr_m_interpolationCoefficients = m_interfaceChildren[level]->a[id].m_interpolationCoefficients;

          const MFloat omega_F =
              2.0
              / (1.0
                 + 6.0 * m_solver->a_nu(currentId) * FFPOW2(m_solver->maxLevel() - m_solver->minLevel() - level - 1));
          const MFloat omega_C =
              2.0
              / (1.0 + 6.0 * m_solver->a_nu(currentId) * FFPOW2(m_solver->maxLevel() - m_solver->minLevel() - level));
          const MFloat omegaT_F = 2.0
                                  / (1.0
                                     + factor * m_solver->a_kappa(currentId)
                                           * FFPOW2(m_solver->maxLevel() - m_solver->minLevel() - level - 1));
          const MFloat omegaT_C =
              2.0
              / (1.0
                 + factor * m_solver->a_kappa(currentId) * FFPOW2(m_solver->maxLevel() - m_solver->minLevel() - level));

          // Spatial interpolation and temporal extrapolation of macroscopic variables
          // f(t+1) = 2*f(t) - f(t-1)
          MFloat rho = 0.0;
          std::array<MFloat, nDim> u{};
          MFloat T = 0.0;

          for(MInt m = 0; m < IPOW2(nDim); m++) {
            rho += 2.0 * ptr_m_interpolationCoefficients[m]
                   * m_solver->a_variable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[m], PV->RHO);
            u[0] += 2.0 * ptr_m_interpolationCoefficients[m]
                    * m_solver->a_variable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[m], PV->U);
            u[1] += 2.0 * ptr_m_interpolationCoefficients[m]
                    * m_solver->a_variable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[m], PV->V);
            IF_CONSTEXPR(nDim == 3)
            u[2] += 2.0 * ptr_m_interpolationCoefficients[m]
                    * m_solver->a_variable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[m], PV->W);
            T += 2.0 * ptr_m_interpolationCoefficients[m]
                 * m_solver->a_variable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[m], PV->T);
            rho -= ptr_m_interpolationCoefficients[m]
                   * m_solver->a_oldVariable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[m], PV->RHO);
            u[0] -= ptr_m_interpolationCoefficients[m]
                    * m_solver->a_oldVariable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[m], PV->U);
            u[1] -= ptr_m_interpolationCoefficients[m]
                    * m_solver->a_oldVariable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[m], PV->V);
            IF_CONSTEXPR(nDim == 3)
            u[2] -= ptr_m_interpolationCoefficients[m]
                    * m_solver->a_oldVariable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[m], PV->W);
            T -= ptr_m_interpolationCoefficients[m]
                 * m_solver->a_oldVariable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[m], PV->T);
          }

          // TODO labels:LB sysEqn: make PV primitive again (not storing rho*u)
          u[0] /= rho;
          u[1] /= rho;
          IF_CONSTEXPR(nDim == 3) u[2] /= rho;

          std::array<MFloat, nDist> eqDist{};
          eqDist = m_solver->getEqDists(rho, u.data());
          std::array<MFloat, nDist> eqDistThermal{};
          eqDistThermal = m_solver->getEqDistsThermal(T, rho, u.data());

          // interpolate and transform
          for(MInt dist = 0; dist < nDist - 1; dist++) {
            // only the missing distributions
            if(m_solver->a_hasNeighbor(currentId, dist) == 0) {
              MFloat tmpDist1 = 0.0;
              MFloat tmpDist1T = 0.0;

              for(MInt m = 0; m < IPOW2(nDim); m++) {
                tmpDist1 += ptr_m_interpolationCoefficients[m]
                            * m_solver->a_oldDistribution(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[m],
                                                          Ld::oppositeDist(dist));
                tmpDist1T += ptr_m_interpolationCoefficients[m]
                             * m_solver->a_oldDistributionThermal(
                                 m_interfaceChildren[level]->a[id].m_interpolationNeighbors[m], Ld::oppositeDist(dist));
              }

              m_solver->a_oldDistribution(currentId, Ld::oppositeDist(dist)) =
                  (tmpDist1 - eqDist[Ld::oppositeDist(dist)]) * (omega_C / omega_F) * F1B2
                  + eqDist[Ld::oppositeDist(dist)];
              m_solver->a_oldDistributionThermal(currentId, Ld::oppositeDist(dist)) =
                  (tmpDist1T - eqDistThermal[Ld::oppositeDist(dist)]) * (omegaT_C / omegaT_F) * F1B2
                  + eqDistThermal[Ld::oppositeDist(dist)];
            }
          }
        }
      }

      //---------------------------
      // if not in sync with parent: apply temporal interpolation
      else {
        for(MInt id = 0; id < m_interfaceChildren[level]->size(); id++) {
          const MInt currentId = m_interfaceChildren[level]->a[id].m_cellId;

          if(currentId >= m_solver->noInternalCells()) continue;

          ptr_m_interpolationCoefficients = m_interfaceChildren[level]->a[id].m_interpolationCoefficients;

          const MFloat omega_F =
              2.0
              / (1.0
                 + 6.0 * m_solver->a_nu(currentId) * FFPOW2(m_solver->maxLevel() - m_solver->minLevel() - level - 1));
          const MFloat omega_C =
              2.0
              / (1.0 + 6.0 * m_solver->a_nu(currentId) * FFPOW2(m_solver->maxLevel() - m_solver->minLevel() - level));
          const MFloat omegaT_F = 2.0
                                  / (1.0
                                     + 6.0 * m_solver->a_kappa(currentId)
                                           * FFPOW2(m_solver->maxLevel() - m_solver->minLevel() - level - 1));
          const MFloat omegaT_C =
              2.0
              / (1.0
                 + 6.0 * m_solver->a_kappa(currentId) * FFPOW2(m_solver->maxLevel() - m_solver->minLevel() - level));

          // spatial and temporal interpolation of variables
          // f(t+1/2) = 1.5*f(t) - 0.5*f(t-1)
          MFloat rho = 0;
          std::array<MFloat, nDim> u{};
          MFloat T = 0.0;

          for(MInt m = 0; m < IPOW2(nDim); m++) {
            rho += 1.5 * ptr_m_interpolationCoefficients[m]
                   * m_solver->a_variable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[m], PV->RHO);
            u[0] += 1.5 * ptr_m_interpolationCoefficients[m]
                    * m_solver->a_variable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[m], PV->U);
            u[1] += 1.5 * ptr_m_interpolationCoefficients[m]
                    * m_solver->a_variable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[m], PV->V);
            IF_CONSTEXPR(nDim == 3)
            u[2] += 1.5 * ptr_m_interpolationCoefficients[m]
                    * m_solver->a_variable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[m], PV->W);
            T += 1.5 * ptr_m_interpolationCoefficients[m]
                 * m_solver->a_variable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[m], PV->T);
            rho -= 0.5 * ptr_m_interpolationCoefficients[m]
                   * m_solver->a_oldVariable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[m], PV->RHO);
            u[0] -= 0.5 * ptr_m_interpolationCoefficients[m]
                    * m_solver->a_oldVariable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[m], PV->U);
            u[1] -= 0.5 * ptr_m_interpolationCoefficients[m]
                    * m_solver->a_oldVariable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[m], PV->V);
            IF_CONSTEXPR(nDim == 3)
            u[2] -= 0.5 * ptr_m_interpolationCoefficients[m]
                    * m_solver->a_oldVariable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[m], PV->W);
            T -= 0.5 * ptr_m_interpolationCoefficients[m]
                 * m_solver->a_oldVariable(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[m], PV->T);
          }

          u[0] /= rho;
          u[1] /= rho;
          IF_CONSTEXPR(nDim == 3) u[2] /= rho;

          std::array<MFloat, nDist> eqDist{};
          eqDist = m_solver->getEqDists(rho, u.data());
          std::array<MFloat, nDist> eqDistThermal{};
          eqDistThermal = m_solver->getEqDistsThermal(T, rho, u.data());

          std::array<MFloat, nDist> forcing{};
          getCellForcing(forcing, currentId);

          // Transform and interpolate
          for(MInt dist = 0; dist < nDist - 1; dist++) {
            // only for missing distributions
            if(m_solver->a_hasNeighbor(currentId, dist) == 0) {
              MFloat tmpDist1 = 0.0;
              MFloat tmpDist2 = 0.0;
              MFloat tmpDist1T = 0.0;
              MFloat tmpDist2T = 0.0;

              for(MInt m = 0; m < IPOW2(nDim); m++) {
                // add the weighted parts to the non-equilibrium distributions using linear temporal interpolation
                // f(t+1/2) = 1/2*f(t) + 1/2*f(t+1)
                if(m_solver->a_hasNeighbor(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[m], dist)) {
                  const MInt neighborId =
                      m_solver->c_neighborId(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[m], dist);

                  // the propagation on the coarse grid has not taken place yet,
                  // thus the incoming distribution for the coarse cell can be interpolated as follows 0.5*(tmpDist1 +
                  // tmpDist2)
                  tmpDist1 +=
                      ptr_m_interpolationCoefficients[m] * m_solver->a_distribution(neighborId, Ld::oppositeDist(dist));
                  tmpDist2 +=
                      ptr_m_interpolationCoefficients[m]
                      * m_solver->a_oldDistribution(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[m],
                                                    Ld::oppositeDist(dist));
                  tmpDist1T += ptr_m_interpolationCoefficients[m]
                               * m_solver->a_distributionThermal(neighborId, Ld::oppositeDist(dist));
                  tmpDist2T +=
                      ptr_m_interpolationCoefficients[m]
                      * m_solver->a_oldDistributionThermal(
                          m_interfaceChildren[level]->a[id].m_interpolationNeighbors[m], Ld::oppositeDist(dist));

                  // add forcing terms for coarse grid before transformation
                  // for tmpDist2 the forcing terms have already been updated
                  tmpDist1 += m_interfaceChildren[level]->a[id].m_interpolationCoefficients[m]
                              * FPOW2(m_solver->maxLevel() - m_solver->minLevel() - level)
                              * forcing[Ld::oppositeDist(dist)];
                }

                // if the interpolation neighbor has no neighbor in actual direction, no temporal interpolation can be
                // applied this may only occur at a non-periodic boundary
                else {
                  tmpDist2 +=
                      ptr_m_interpolationCoefficients[m]
                      * m_solver->a_oldDistribution(m_interfaceChildren[level]->a[id].m_interpolationNeighbors[m],
                                                    Ld::oppositeDist(dist));
                  tmpDist2T +=
                      ptr_m_interpolationCoefficients[m]
                      * m_solver->a_oldDistributionThermal(
                          m_interfaceChildren[level]->a[id].m_interpolationNeighbors[m], Ld::oppositeDist(dist));

                  // set tmpDist2=tmpDist1
                  // forcing term is already included
                  tmpDist1 = tmpDist2;
                  tmpDist1T = tmpDist2T;
                }
              }

              m_solver->a_oldDistribution(currentId, Ld::oppositeDist(dist)) =
                  (0.5 * (tmpDist1 + tmpDist2) - eqDist[Ld::oppositeDist(dist)]) * (omega_C / omega_F) * F1B2
                  + eqDist[Ld::oppositeDist(dist)];
              m_solver->a_oldDistributionThermal(currentId, Ld::oppositeDist(dist)) =
                  (0.5 * (tmpDist1T + tmpDist2T) - eqDistThermal[Ld::oppositeDist(dist)]) * (omegaT_C / omegaT_F) * F1B2
                  + eqDistThermal[Ld::oppositeDist(dist)];
            }
          }
        }
      } // end of else
    }
  }
}

/** \brief Fine to coarse grid
 *
 * average all child cell values
 * and transform the non-eq parts according to Dupuis et al.
 *
 */
template <MInt nDim, MInt nDist, class SysEqn>
template <MBool compressible>
void LbInterfaceDxQy<nDim, nDist, SysEqn>::restrictionDupuis_() {
  TRACE();

#ifdef WAR_NVHPC_PSTL
  const MInt noInternalCells = m_solver->noInternalCells();
#endif

  for(MInt level = 0; level < (m_solver->maxLevel() - m_solver->minLevel());
      level++) { // level is the actual coarse resolution
                 // level+1 is the actual fine resolution

    const MFloat nu_F = m_solver->m_Ma * LBCS / m_solver->m_Re * m_solver->m_referenceLength
                        * FFPOW2(m_solver->maxLevel() - m_solver->minLevel() - level - 1);
    const MFloat omega_F = 2.0 / (1.0 + 6.0 * nu_F);
    const MFloat nu_C = m_solver->m_Ma * LBCS / m_solver->m_Re * m_solver->m_referenceLength
                        * FFPOW2(m_solver->maxLevel() - m_solver->minLevel() - level);
    const MFloat omega_C = 2.0 / (1.0 + 6.0 * nu_C);


    if((globalTimeStep) % IPOW2(m_solver->maxLevel() - m_solver->minLevel() - level)
       == 0) { // perform restriction from finest level in every second timestep,
               // on the next coarser level every fourth timestep, and so on; starting from t = 1

      maia::parallelFor<true>(0, m_interfaceParents[level]->size(), [=](MInt id) {
        const MInt currentId = m_interfaceParents[level]->a[id].m_cellId;
#ifdef WAR_NVHPC_PSTL
        if(currentId >= noInternalCells) return;
#else
        if(currentId >= m_solver->noInternalCells()) return;
#endif

        MFloat rho = 0;
        std::array<MFloat, nDim> u{};

        // temporal extrapolated mean of the macroscopic variables of the child cells
        // f(t+1) = 2*f(t) - f(t-1)
        for(MInt m = 0; m < IPOW2(nDim); m++) {
          const MInt childCellId = m_solver->c_childId(currentId, m);
          rho += FFPOW2(nDim)
                 * (2.0 * m_solver->a_variable(childCellId, PV->RHO) - m_solver->a_oldVariable(childCellId, PV->RHO));
          for(MInt n = 0; n < nDim; n++) {
            u[n] +=
                FFPOW2(nDim) * (2.0 * m_solver->a_variable(childCellId, n) - m_solver->a_oldVariable(childCellId, n));
          }
        }

        std::array<MFloat, nDist> eqDist;
        eqDist = m_solver->getEqDists(rho, u.data());

        // Interpolate missing distributions
        for(MInt dist = 0; dist < (nDist - 1); dist++) {
#ifdef WAR_NVHPC_PSTL
          MInt oppositeDist = m_solver->m_oppositeDist[dist];
#else
          MInt oppositeDist = Ld::oppositeDist(dist);
#endif
          if(!m_solver->a_isInterfaceParent(m_solver->c_neighborId(currentId, oppositeDist))
             && !m_solver->c_isLeafCell(m_solver->c_neighborId(currentId, oppositeDist))) {
            // Calculate mean distribution from all children
            MFloat tmpDist = 0.0;
            for(MInt m = 0; m < IPOW2(nDim); m++) {
              tmpDist += FFPOW2(nDim) * m_solver->a_oldDistribution(m_solver->c_childId(currentId, m), dist);
            }

            // transform distribution
            m_solver->a_oldDistribution(currentId, dist) =
                (tmpDist - eqDist[dist]) * (omega_F / omega_C) * 2.0 + eqDist[dist];
          }
        }
      });
    }
  }
}

template <MInt nDim, MInt nDist, class SysEqn>
void LbInterfaceDxQy<nDim, nDist, SysEqn>::restrictionDupuis() {
  restrictionDupuis_<false>();
}

template <MInt nDim, MInt nDist, class SysEqn>
void LbInterfaceDxQy<nDim, nDist, SysEqn>::restrictionDupuisCompressible() {
  restrictionDupuis_<true>();
}

/** \brief Fine to coarse grid for thermal LB
 *
 * \date 31.07.2012
 * \author Andreas Lintermann
 *
 * \todo labels:LB This is currently not working properly, please use ROHDE method instead.
 *
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbInterfaceDxQy<nDim, nDist, SysEqn>::restrictionThermalDupuis() {
  const MFloat factor = (m_solver->m_innerEnergy) ? ((nDim == 2) ? 3.0 : 18.0 / 5.0) : 6.0;
  for(MInt level = 0; level < (m_solver->maxLevel() - m_solver->minLevel());
      level++) { // level is the actual coarse resolution
                 // level+1 is the actual fine resolution

    if((globalTimeStep) % IPOW2(m_solver->maxLevel() - m_solver->minLevel() - level)
       == 0) { // perform restriction from finest level in every second timestep,
               // on the next coarser level every fourth timestep, and so on; starting from t = 1

      for(MInt id = 0; id < m_interfaceParents[level]->size(); id++) {
        const MInt currentId = m_interfaceParents[level]->a[id].m_cellId;
        if(currentId >= m_solver->noInternalCells()) continue;

        const MFloat omega_F =
            2.0
            / (1.0 + 6.0 * m_solver->a_nu(currentId) * FFPOW2(m_solver->maxLevel() - m_solver->minLevel() - level - 1));
        const MFloat omega_C =
            2.0 / (1.0 + 6.0 * m_solver->a_nu(currentId) * FFPOW2(m_solver->maxLevel() - m_solver->minLevel() - level));
        const MFloat omegaT_F = 2.0
                                / (1.0
                                   + factor * m_solver->a_kappa(currentId)
                                         * FFPOW2(m_solver->maxLevel() - m_solver->minLevel() - level - 1));
        const MFloat omegaT_C =
            2.0
            / (1.0
               + factor * m_solver->a_kappa(currentId) * FFPOW2(m_solver->maxLevel() - m_solver->minLevel() - level));

        MFloat rho = F0;
        std::array<MFloat, nDim> u{};
        MFloat T = F0;

        // temporal extrapolated mean of the macroscopic variables of the child cells
        // f(t+1) = 2*f(t) - f(t-1)
        for(MInt m = 0; m < IPOW2(nDim); m++) {
          rho += F2 * FFPOW2(nDim) * m_solver->a_variable(m_solver->c_childId(currentId, m), PV->RHO);
          u[0] += F2 * FFPOW2(nDim) * m_solver->a_variable(m_solver->c_childId(currentId, m), PV->U);
          u[1] += F2 * FFPOW2(nDim) * m_solver->a_variable(m_solver->c_childId(currentId, m), PV->V);
          IF_CONSTEXPR(nDim == 3)
          u[2] += F2 * FFPOW2(nDim) * m_solver->a_variable(m_solver->c_childId(currentId, m), PV->W);
          T += F2 * FFPOW2(nDim) * m_solver->a_variable(m_solver->c_childId(currentId, m), PV->T);
          rho -= FFPOW2(nDim) * m_solver->a_oldVariable(m_solver->c_childId(currentId, m), PV->RHO);
          u[0] -= FFPOW2(nDim) * m_solver->a_oldVariable(m_solver->c_childId(currentId, m), PV->U);
          u[1] -= FFPOW2(nDim) * m_solver->a_oldVariable(m_solver->c_childId(currentId, m), PV->V);
          IF_CONSTEXPR(nDim == 3)
          u[2] -= FFPOW2(nDim) * m_solver->a_oldVariable(m_solver->c_childId(currentId, m), PV->W);
          T -= FFPOW2(nDim) * m_solver->a_oldVariable(m_solver->c_childId(currentId, m), PV->T);
        }

        for(MInt d = 0; d < nDim; d++) {
          u[d] /= rho;
        }

        std::array<MFloat, nDist> eqDist{};
        eqDist = m_solver->getEqDists(rho, u.data());
        std::array<MFloat, nDist> eqDistThermal{};
        eqDistThermal = m_solver->getEqDistsThermal(T, rho, u.data());

        // Interpolate missing distributions
        for(MInt dist = 0; dist < (nDist - 1); dist++) {
          if(!m_solver->a_isInterfaceParent(m_solver->c_neighborId(currentId, Ld::oppositeDist(dist)))
             && !m_solver->c_isLeafCell(m_solver->c_neighborId(currentId, Ld::oppositeDist(dist)))) {
            // Calculate mean distribution from all children
            MFloat tmpDist = 0.0;
            MFloat tmpDistT = 0.0;
            for(MInt m = 0; m < IPOW2(nDim); m++) {
              tmpDist += FFPOW2(nDim) * m_solver->a_oldDistribution(m_solver->c_childId(currentId, m), dist);
              tmpDistT += FFPOW2(nDim) * m_solver->a_oldDistributionThermal(m_solver->c_childId(currentId, m), dist);
            }

            // transform distribution
            m_solver->a_oldDistribution(currentId, dist) =
                (tmpDist - eqDist[dist]) * (omega_F / omega_C) * 2.0 + eqDist[dist];
            m_solver->a_oldDistributionThermal(currentId, dist) =
                (tmpDistT - eqDistThermal[dist]) * (omegaT_F / omegaT_C) * 2.0 + eqDistThermal[dist];
          }
        }
      }
    }
  }
}

/** coarse to fine grid
 *
 * outgoing distributions are handed over from parent to children
 *
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbInterfaceDxQy<nDim, nDist, SysEqn>::prolongationRohde() {
  std::array<MFloat, nDist> forcing{};
  getCellForcing(forcing);

  for(MInt level = 0; level < (m_solver->maxLevel() - m_solver->minLevel());
      level++) { // minLevel+level is the actual coarse resolution
                 // minLevel+level+1 is the actual fine resolution

    if((globalTimeStep) % IPOW2(m_solver->maxLevel() - m_solver->minLevel() - level - 1) == 0) {
      //-----------------------
      // intermittant prolongation
      if((globalTimeStep) % IPOW2(m_solver->maxLevel() - m_solver->minLevel() - level) == 0) {
        for(MInt id = 0; id < m_interfaceParents[level]->size(); id++) {
          const MInt currentId = m_interfaceParents[level]->a[id].m_cellId;
          if(currentId >= m_solver->noInternalCells()) continue;
          for(MInt m = 0; m < IPOW2(nDim); m++) {
            // compensate collision on fine grid
            if(m_cellDependentForcing) getCellForcing(forcing, m_solver->c_childId(currentId, m));
            for(MInt dist = 0; dist < nDist; dist++) {
              // remove forcing term
              m_solver->a_oldDistribution(m_solver->c_childId(currentId, m), dist) -=
                  FPOW2(m_solver->maxLevel() - m_solver->minLevel() - level - 1) * forcing[dist];

              // reset outgoing distributions
              m_solver->a_distribution(m_solver->c_childId(currentId, m), dist) =
                  m_solver->a_oldDistribution(m_solver->c_childId(currentId, m), dist);
            }
            // copy variables from parent
            for(MInt var = 0; var < (nDim + 1); var++) {
              m_solver->a_variable(m_solver->c_childId(currentId, m), var) = m_solver->a_variable(currentId, var);
            }
          }
        }

      }
      //--------------------------
      // regular prolongation
      else {
        for(MInt id = 0; id < m_interfaceParents[level]->size(); id++) {
          const MInt currentId = m_interfaceParents[level]->a[id].m_cellId;
          if(currentId >= m_solver->noInternalCells()) continue;

          // overwrite variables and distributions
          for(MInt m = 0; m < IPOW2(nDim); m++) {
            // hand over distributions
            if(m_cellDependentForcing) getCellForcing(forcing, m_solver->c_childId(currentId, m));
            for(MInt dist = 0; dist < nDist; dist++) {
              m_solver->a_distribution(m_solver->c_childId(currentId, m), dist) =
                  m_solver->a_distribution(currentId, dist);

              // add half of forcing term for coarse grid
              m_solver->a_distribution(m_solver->c_childId(currentId, m), dist) +=
                  F1B2 * FPOW2(m_solver->maxLevel() - m_solver->minLevel() - level) * forcing[dist];
              // remove half of forcing term for fine grid
              m_solver->a_distribution(m_solver->c_childId(currentId, m), dist) -=
                  F1B2 * FPOW2(m_solver->maxLevel() - m_solver->minLevel() - level - 1) * forcing[dist];
            }

            // copy variables from parent
            for(MInt var = 0; var < (nDim + 1); var++) {
              m_solver->a_variable(m_solver->c_childId(currentId, m), var) = m_solver->a_variable(currentId, var);
            }
          }
        }
      }
    }
  }
}

/** \brief Coarse to fine grid for thermal LB
 *
 * \date 31.07.2012
 * \author Andreas Linermann
 *
 * Outgoing distributions are handed over from parent to children.
 *
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbInterfaceDxQy<nDim, nDist, SysEqn>::prolongationThermalRohde() {
  std::array<MFloat, nDist> forcing{};
  getCellForcing(forcing);

  for(MInt level = 0; level < (m_solver->maxLevel() - m_solver->minLevel());
      level++) { // minLevel+level is the actual coarse resolution
                 // minLevel+level+1 is the actual fine resolution

    if((globalTimeStep) % IPOW2(m_solver->maxLevel() - m_solver->minLevel() - level - 1) == 0) {
      //-----------------------
      // intermittant prolongation
      if((globalTimeStep) % IPOW2(m_solver->maxLevel() - m_solver->minLevel() - level) == 0) {
        for(MInt id = 0; id < m_interfaceParents[level]->size(); id++) {
          const MInt currentId = m_interfaceParents[level]->a[id].m_cellId;
          if(currentId >= m_solver->noInternalCells()) continue;

          for(MInt m = 0; m < IPOW2(nDim); m++) {
            // compensate collision on fine grid
            if(m_cellDependentForcing) getCellForcing(forcing, m_solver->c_childId(currentId, m));
            for(MInt dist = 0; dist < nDist; dist++) {
              // remove forcing term
              m_solver->a_oldDistribution(m_solver->c_childId(currentId, m), dist) -=
                  FPOW2(m_solver->maxLevel() - m_solver->minLevel() - level - 1) * forcing[dist];

              // reset outgoing distributions
              m_solver->a_distribution(m_solver->c_childId(currentId, m), dist) =
                  m_solver->a_oldDistribution(m_solver->c_childId(currentId, m), dist);
              m_solver->a_distributionThermal(m_solver->c_childId(currentId, m), dist) =
                  m_solver->a_oldDistributionThermal(m_solver->c_childId(currentId, m), dist);
            }
            // copy variables from parent
            for(MInt var = 0; var < (nDim + 2); var++) {
              m_solver->a_variable(m_solver->c_childId(currentId, m), var) = m_solver->a_variable(currentId, var);
            }
          }
        }

      }
      //--------------------------
      // regular prolongation
      else {
        for(MInt id = 0; id < m_interfaceParents[level]->size(); id++) {
          const MInt currentId = m_interfaceParents[level]->a[id].m_cellId;
          if(currentId >= m_solver->noInternalCells()) continue;

          // overwrite variables and distributions
          for(MInt m = 0; m < IPOW2(nDim); m++) {
            // hand over distributions
            if(m_cellDependentForcing) getCellForcing(forcing, m_solver->c_childId(currentId, m));
            for(MInt dist = 0; dist < nDist; dist++) {
              m_solver->a_distribution(m_solver->c_childId(currentId, m), dist) =
                  m_solver->a_distribution(currentId, dist);
              m_solver->a_distributionThermal(m_solver->c_childId(currentId, m), dist) =
                  m_solver->a_distributionThermal(currentId, dist);

              // add half of forcing term for coarse grid
              m_solver->a_distribution(m_solver->c_childId(currentId, m), dist) +=
                  F1B2 * FPOW2(m_solver->maxLevel() - m_solver->minLevel() - level) * forcing[dist];
              // remove half of forcing term for fine grid
              m_solver->a_distribution(m_solver->c_childId(currentId, m), dist) -=
                  F1B2 * FPOW2(m_solver->maxLevel() - m_solver->minLevel() - level - 1) * forcing[dist];
            }

            // copy variables from parent
            for(MInt var = 0; var < (nDim + 2); var++) {
              m_solver->a_variable(m_solver->c_childId(currentId, m), var) = m_solver->a_variable(currentId, var);
            }
          }
        }
      }
    }
  }
}


/** fine to coarse grid
 *
 * missing incoming distributions of parent are filled with average child cell values
 *
 */

template <MInt nDim, MInt nDist, class SysEqn>
void LbInterfaceDxQy<nDim, nDist, SysEqn>::restrictionRohde() {
  std::array<MFloat, nDist> forcing{};
  getCellForcing(forcing);

  for(MInt level = 0; level < (m_solver->maxLevel() - m_solver->minLevel());
      level++) { // minLevel+level is the actual coarse resolution
    // minLevel+level+1 is the actual fine resolution

    if((globalTimeStep) % IPOW2(m_solver->maxLevel() - m_solver->minLevel() - level) == 0) {
      for(MInt id = 0; id < m_interfaceParents[level]->size(); id++) {
        const MInt currentId = m_interfaceParents[level]->a[id].m_cellId;
        if(currentId >= m_solver->noInternalCells()) continue;
        if(m_cellDependentForcing) getCellForcing(forcing, currentId);

        for(MInt dist = 0; dist < (nDist - 1); dist++) {
          if(!m_solver->a_isInterfaceParent(m_solver->c_neighborId(currentId, Ld::oppositeDist(dist)))
             && !m_solver->c_isLeafCell(m_solver->c_neighborId(currentId, Ld::oppositeDist(dist)))) {
            m_solver->a_oldDistribution(currentId, dist) = F0;

            for(MInt m = 0; m < IPOW2(nDim); m++) {
              m_solver->a_oldDistribution(currentId, dist) +=
                  FFPOW2(nDim) * m_solver->a_oldDistribution(m_solver->c_childId(currentId, m), dist);
            }

            // add half of forcing term for fine grid
            m_solver->a_oldDistribution(currentId, dist) +=
                F1B2 * FPOW2(m_solver->maxLevel() - m_solver->minLevel() - level - 1) * forcing[dist];
            // compensate half of forcing term for coarse grid
            m_solver->a_oldDistribution(currentId, dist) -=
                F1B2 * FPOW2(m_solver->maxLevel() - m_solver->minLevel() - level) * forcing[dist];
          }
        }
      }
    }
  }
}

/** \brief Fine to coarse grid fot thermal LB
 *
 * \date 31.07.2012
 * \author Andreas Linermann
 *
 * Missing incoming distributions of parent are filled with average child cell values.
 *
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbInterfaceDxQy<nDim, nDist, SysEqn>::restrictionThermalRohde() {
  std::array<MFloat, nDist> forcing{};
  getCellForcing(forcing);

  for(MInt level = 0; level < (m_solver->maxLevel() - m_solver->minLevel());
      level++) { // minLevel+level is the actual coarse resolution
    // minLevel+level+1 is the actual fine resolution

    if((globalTimeStep) % IPOW2(m_solver->maxLevel() - m_solver->minLevel() - level) == 0) {
      for(MInt id = 0; id < m_interfaceParents[level]->size(); id++) {
        const MInt currentId = m_interfaceParents[level]->a[id].m_cellId;
        if(currentId >= m_solver->noInternalCells()) continue;
        if(m_cellDependentForcing) getCellForcing(forcing, currentId);

        for(MInt dist = 0; dist < (nDist - 1); dist++) {
          if(!m_solver->a_isInterfaceParent(m_solver->c_neighborId(currentId, Ld::oppositeDist(dist)))
             && !m_solver->c_isLeafCell(m_solver->c_neighborId(currentId, Ld::oppositeDist(dist)))) {
            m_solver->a_oldDistribution(currentId, dist) = F0;
            m_solver->a_oldDistributionThermal(currentId, dist) = F0;

            for(MInt m = 0; m < IPOW2(nDim); m++) {
              m_solver->a_oldDistribution(currentId, dist) +=
                  FFPOW2(nDim) * m_solver->a_oldDistribution(m_solver->c_childId(currentId, m), dist);
              m_solver->a_oldDistributionThermal(currentId, dist) +=
                  FFPOW2(nDim) * m_solver->a_oldDistributionThermal(m_solver->c_childId(currentId, m), dist);
            }

            // add half of forcing term for fine grid
            m_solver->a_oldDistribution(currentId, dist) +=
                F1B2 * FPOW2(m_solver->maxLevel() - m_solver->minLevel() - level - 1) * forcing[dist];
            // compensate half of forcing term for coarse grid
            m_solver->a_oldDistribution(currentId, dist) -=
                F1B2 * FPOW2(m_solver->maxLevel() - m_solver->minLevel() - level) * forcing[dist];
          }
        }
      }
    }
  }
}

//------------------ ADAPTATION -----------------------//

/**
 * \brief Setting the adaptation functions chosen via property
 *
 * \author Julian Vorspohl <j.vorspohl@aia.rwth-aachen.de>
 */
template <MInt nDim, MInt nDist, class SysEqn>
void LbInterfaceDxQy<nDim, nDist, SysEqn>::setAdaptationFunctions() {
  TRACE();

  switch(string2enum(m_adaptationInitMethod)) {
    // Default property value
    case INIT_DUPUIS_FILIPPOVA: {
      fRefineCell = &LbInterfaceDxQy::refineCellDupuis;
      fRemoveChildren = &LbInterfaceDxQy::removeChildsDupuisFilippova;
      break;
    }
    case INIT_COPYPASTE: {
      fRefineCell = &LbInterfaceDxQy::refineCellCopyPaste;
      fRemoveChildren = &LbInterfaceDxQy::removeChildsCopyPaste;
      break;
    }
    default: {
      stringstream errorMessage;
      errorMessage << " LbInterface::setAdaptationFunctions: Specified interface condition: " << m_adaptationInitMethod
                   << " does not exist. Exiting!";
      mTerm(1, AT_, errorMessage.str());
    }
  }
}

template <MInt nDim, MInt nDist, class SysEqn>
void LbInterfaceDxQy<nDim, nDist, SysEqn>::refineCell(const MInt parentId, const MInt* childIds) {
  (this->*fRefineCell)(parentId, childIds);
}

template <MInt nDim, MInt nDist, class SysEqn>
void LbInterfaceDxQy<nDim, nDist, SysEqn>::removeChildren(const MInt parentId) {
  (this->*fRemoveChildren)(parentId);
}

/**
 * \brief     Initialize child variables from parent
 * \details   - interpolation
 *            - temporal extrapolation
 *            - Dupuis
 * \param[in] parentId: Solver cell id of cell that will be refined
 * \param[in] childIds: Solver cell ids of new children
 * \author    Philipp Brokof, Moritz Waldmann
 **/
template <MInt nDim, MInt nDist, class SysEqn>
void LbInterfaceDxQy<nDim, nDist, SysEqn>::refineCellDupuis(const MInt parentId, const MInt* childIds) {
  // Halo cells are initialized by exchange
  // ToDo: This check is unnecessary by now
  //       because this function should not be called for Halos anymore.
  if(!m_solver->a_isHalo(parentId)) {
    // Find parent neighbors
    // ToDo: Check if c_neighbor of solver can be used. Should be possible now because
    // this function is calles at the end of reinitAfterAdaptation and c_neighbor is
    // updated before.
    MInt parentNeighbors[IPOW3[nDim] - 1];
    for(MInt n = 0; n < IPOW3[nDim] - 1; n++) {
      parentNeighbors[n] = m_solver->c_neighborId(parentId, n);
    }

    const MInt parentLevel = m_solver->c_level(parentId);
    const MInt dParentLevel = m_solver->maxLevel() - parentLevel;
    const MFloat omega_F = 2.0 / (1.0 + 6.0 * m_solver->a_nu(parentId) * FFPOW2(dParentLevel - 1));
    const MFloat omega_C = 2.0 / (1.0 + 6.0 * m_solver->a_nu(parentId) * FFPOW2(dParentLevel));
    std::array<MFloat, nDist> forcing{};
    getCellForcing(forcing);

    // Parent and child level have done propagation in the last timestep and will do collision after this...
    if((globalTimeStep - 1) % IPOW2(dParentLevel) == 0) {
      // Loop over all children that are created by refining the cell
      for(MInt child = 0; child < IPOW2(nDim); child++) {
        if(childIds[child] < 0) continue;
        const MInt childId = childIds[child];

        m_solver->a_nu(childId) = m_solver->a_nu(parentId);

        // Get position of child in parent
        MInt position = 0;
        for(MInt dim = 0; dim < nDim; dim++) {
          if(m_solver->a_coordinate(childId, dim) < m_solver->a_coordinate(parentId, dim)) {
            position += 0;
          } else {
            position += IPOW2(dim);
          }
        }

        // Set interpolation neighbors
        // All interpolation neighbors are on parent level
        MInt interpolationNeighbors[IPOW2(nDim)];
        for(MInt k = 0; k < IPOW2(nDim); k++) {
          const MInt dir = Ld::intNghbrArray(position, k);
          if(dir < IPOW3[nDim] - 1) {
            interpolationNeighbors[k] = parentNeighbors[dir];
          } else {
            interpolationNeighbors[k] = parentId;
          }
        }
        // Set interpolation coefficients
        MFloat interpolationCoefficients[IPOW2(nDim)];
        for(MInt k = 0; k < IPOW2(nDim); k++) {
          interpolationCoefficients[k] = Ld::linearInterpolationCoefficients(position, k);
        }

        MFloat rho = 0;
        std::array<MFloat, nDim> u{};
        MFloat rhoOld = 0;
        std::array<MFloat, nDim> uOld{};

        // Spatial interpolation
        for(MInt m = 0; m < IPOW2(nDim); m++) {
          const MInt interpNeighbor = interpolationNeighbors[m];
          const MFloat interpCoeff = interpolationCoefficients[m];

          // Macroscopic vars at t
          rho += interpCoeff * m_solver->a_variable(interpNeighbor, PV->RHO);
          for(MInt n = 0; n < nDim; n++) {
            u[n] += interpCoeff * m_solver->a_variable(interpNeighbor, PV->VV[n]);
          }

          // Macroscopic vars at t-1
          rhoOld += interpCoeff * m_solver->a_variable(interpNeighbor, PV->RHO);
          for(MInt n = 0; n < nDim; n++) {
            uOld[n] += interpCoeff * m_solver->a_variable(interpNeighbor, PV->VV[n]);
          }
        }

        // Temporal extrapolation: f(t+1) = 2*f(t) - f(t-1)
        // Needed for initialization of oldDistribution
        const MFloat rho_tp1 = 2 * rho - rhoOld;
        std::array<MFloat, nDim> u_tp1{};
        for(MInt d = 0; d < nDim; d++) {
          u_tp1[d] = 2 * u[d] - uOld[d];
        }

        // Temporal extrapolation: f(t+1/2) = 1.5*f(t) - 0.5*f(t-1)
        // Needed for initialization of macroscopic variables
        MFloat rho_tp12 = 1.5 * rho - 0.5 * rhoOld;
        std::array<MFloat, nDim> u_tp12{};
        for(MInt d = 0; d < nDim; d++) {
          u_tp12[d] = 1.5 * u[d] - 0.5 * uOld[d];
        }

        // Setting macroscopic variables in child cell
        m_solver->a_variable(childId, PV->RHO) = rho_tp12;
        for(MInt n = 0; n < nDim; n++) {
          m_solver->a_variable(childId, PV->VV[n]) = u_tp12[n];
        }
        m_solver->a_oldVariable(childId, PV->RHO) = rho;
        for(MInt n = 0; n < nDim; n++) {
          m_solver->a_oldVariable(childId, PV->VV[n]) = u[n];
        }

        // Calculate all equilibrium distribution for tp1
        std::array<MFloat, nDist> eqDist_tp1{};
        eqDist_tp1 = m_solver->getEqDists(rho_tp1, u_tp1.data());

        // Transform and interpolate all distributions
        for(MInt dist = 0; dist < nDist; dist++) {
          MFloat tmpDist = 0.0;
          for(MInt m = 0; m < IPOW2(nDim); m++) {
            tmpDist += interpolationCoefficients[m] * m_solver->a_oldDistribution(interpolationNeighbors[m], dist);
          }

          // Distribution before next collision
          m_solver->a_oldDistribution(childId, dist) =
              (tmpDist - eqDist_tp1[dist]) * (omega_C / omega_F) * F1B2 + eqDist_tp1[dist];

          // All cells need to perform collision step next to initialize this variable
          m_solver->a_distribution(childId, dist) = m_solver->a_oldDistribution(childId, dist);
        }
      }
    } else {
      // Loop over all children that are created by refining the cell
      for(MInt child = 0; child < IPOW2(nDim); child++) {
        if(childIds[child] < 0) continue;
        const MInt childId = childIds[child];

        m_solver->a_nu(childId) = m_solver->a_nu(parentId);

        // Get position of child in parent
        MInt position = 0;
        for(MInt dim = 0; dim < nDim; dim++) {
          if(m_solver->a_coordinate(childId, dim) < m_solver->a_coordinate(parentId, dim)) {
            position += 0;
          } else {
            position += IPOW2(dim);
          }
        }

        // Set interpolation neighbors
        // All interpolation neighbors are on parent level
        MInt interpolationNeighbors[IPOW2(nDim)];
        for(MInt k = 0; k < IPOW2(nDim); k++) {
          const MInt dir = Ld::intNghbrArray(position, k);
          if(dir < IPOW3[nDim] - 1) {
            interpolationNeighbors[k] = parentNeighbors[dir];
          } else {
            interpolationNeighbors[k] = parentId;
          }
        }
        // Set interpolation coefficients
        MFloat interpolationCoefficients[IPOW2(nDim)];
        for(MInt k = 0; k < IPOW2(nDim); k++) {
          interpolationCoefficients[k] = Ld::linearInterpolationCoefficients(position, k);
        }

        MFloat rho = 0;
        std::array<MFloat, nDim> u{};
        MFloat rhoOld = 0;
        std::array<MFloat, nDim> uOld{};

        // Spatial interpolation
        for(MInt m = 0; m < IPOW2(nDim); m++) {
          const MInt interpNeighbor = interpolationNeighbors[m];
          const MFloat interpCoeff = interpolationCoefficients[m];

          // Macroscopic vars at t
          rho += interpCoeff * m_solver->a_variable(interpNeighbor, PV->RHO);
          for(MInt n = 0; n < nDim; n++) {
            u[n] += interpCoeff * m_solver->a_variable(interpNeighbor, PV->VV[n]);
          }

          // Macroscopic vars at t-1
          rhoOld += interpCoeff * m_solver->a_variable(interpNeighbor, PV->RHO);
          for(MInt n = 0; n < nDim; n++) {
            uOld[n] += interpCoeff * m_solver->a_variable(interpNeighbor, PV->VV[n]);
          }
        }

        // Temporal extrapolation: f(t+1/2) = 1.5*f(t) - 0.5*f(t-1)
        // Needed for initialization of oldDistribution
        const MFloat rho_tp1 = 1.5 * rho - 0.5 * rhoOld;
        std::array<MFloat, nDim> u_tp1{};
        for(MInt d = 0; d < nDim; d++) {
          u_tp1[d] = 1.5 * u[d] - 0.5 * uOld[d];
        }

        // Temporal extrapolation: f(t-1/2) = 0.5*f(t) + 0.5*f(t-1)
        // Needed for initialization of macroscopic variables
        MFloat rho_tp12 = 0.5 * rho + 0.5 * rhoOld;
        std::array<MFloat, nDim> u_tp12{};
        for(MInt d = 0; d < nDim; d++) {
          u_tp12[d] = 0.5 * u[d] + 0.5 * uOld[d];
        }

        // Setting macroscopic variables in child cell
        m_solver->a_variable(childId, PV->RHO) = rho;
        for(MInt n = 0; n < nDim; n++) {
          m_solver->a_variable(childId, PV->VV[n]) = u[n];
        }
        m_solver->a_oldVariable(childId, PV->RHO) = rho_tp12;
        for(MInt n = 0; n < nDim; n++) {
          m_solver->a_oldVariable(childId, PV->VV[n]) = u_tp12[n];
        }

        // Calculate all equilibrium distribution for tp1
        std::array<MFloat, nDist> eqDist_tp1{};
        eqDist_tp1 = m_solver->getEqDists(rho_tp1, u_tp1.data());

        // Transform and interpolate all distributions
        MFloat sumOfOldDist = F0;
        for(MInt dist = 0; dist < nDist - 1; dist++) {
          MBool propagateNeighborValue = false;

          if(m_solver->a_hasNeighbor(childId, Ld::oppositeDist(dist))) {
            if(!m_solver->a_hasProperty(childId, LbCell::WasNewlyCreated)) {
              propagateNeighborValue = true;
            }
          }

          if(propagateNeighborValue) {
            m_solver->a_oldDistribution(childId, dist) =
                m_solver->a_distribution(m_solver->c_neighborId(childId, Ld::oppositeDist(dist)), dist)
                + FPOW2(m_solver->maxLevel() - m_solver->a_level(childId)) * forcing[dist];
          } else {
            MFloat tmpDist1 = 0.0;
            MFloat tmpDist2 = 0.0;
            for(MInt m = 0; m < IPOW2(nDim); m++) {
              MFloat tmpDist =
                  interpolationCoefficients[m] * m_solver->a_oldDistribution(interpolationNeighbors[m], dist);
              if(m_solver->a_hasNeighbor(interpolationNeighbors[m], Ld::oppositeDist(dist))) {
                const MInt neighborId = m_solver->c_neighborId(interpolationNeighbors[m], Ld::oppositeDist(dist));
                tmpDist1 += interpolationCoefficients[m] * m_solver->a_distribution(neighborId, dist);
                tmpDist1 += interpolationCoefficients[m] * FPOW2(m_solver->maxLevel() - m_solver->a_level(childId) + 1)
                            * forcing[dist];
                tmpDist2 += tmpDist;
              } else {
                // if the interpolation neighbor has no neighbor in actual direction, no temporal interpolation can be
                // applied this may only occur at a non-periodic boundary
                tmpDist2 += tmpDist;
                tmpDist1 += tmpDist;
              }
            }

            // Distribution before next collision
            m_solver->a_oldDistribution(childId, dist) =
                (0.5 * (tmpDist1 + tmpDist2) - eqDist_tp1[dist]) * (omega_C / omega_F) * F1B2 + eqDist_tp1[dist];
          }

          // All cells need to perform collision step next to initialize this variable
          m_solver->a_distribution(childId, dist) = m_solver->a_oldDistribution(childId, dist);

          sumOfOldDist += m_solver->a_oldDistribution(childId, dist);
        }
        // Dont forget the rest distribution
        m_solver->a_oldDistribution(childId, Ld::lastId()) = rho - sumOfOldDist;
        m_solver->a_distribution(childId, Ld::lastId()) = rho - sumOfOldDist;
      }
    }
  }
}

/**
 * \brief     Initialize child variables from parent
 * \details   Copy parent --> children
 * \param[in] parentId: Solver cell id of parent that will be coarsen
 *            childIds: Solver ids of new children
 * \author: Philipp Brokof, Moritz Waldmann
 **/
template <MInt nDim, MInt nDist, class SysEqn>
void LbInterfaceDxQy<nDim, nDist, SysEqn>::refineCellCopyPaste(const MInt parentId, const MInt* childIds) {
  if(!m_solver->a_isHalo(parentId)) {
    for(MInt child = 0; child < IPOW2(nDim); child++) {
      if(childIds[child] < 0) continue;
      for(MInt v = 0; v < m_solver->noVariables(); v++) {
        m_solver->a_variable(childIds[child], v) = m_solver->a_variable(parentId, v);
        m_solver->a_oldVariable(childIds[child], v) = m_solver->a_oldVariable(parentId, v);
      }
      for(MInt dir = 0; dir < nDist; dir++) {
        m_solver->a_distribution(childIds[child], dir) = m_solver->a_distribution(parentId, dir);
        m_solver->a_oldDistribution(childIds[child], dir) = m_solver->a_oldDistribution(parentId, dir);
      }
    }
  }
}

/**
 * \brief     Initialize parent variables from children
 * \details   Average children --> parent
 * \param[in] parentId: Solver cell id of parent that will be coarsen
 * \author: Philipp Brokof, Moritz Waldmann
 **/
template <MInt nDim, MInt nDist, class SysEqn>
void LbInterfaceDxQy<nDim, nDist, SysEqn>::removeChildsCopyPaste(const MInt parentId) {
  if(!m_solver->a_isHalo(parentId)) {
    std::array<MFloat, nDist> forcing{};
    getCellForcing(forcing);

    // Get children that are to remove
    MInt children[IPOW2(nDim)];
    for(MInt m = 0; m < IPOW2(nDim); m++) {
      children[m] = m_solver->c_childId(parentId, m);
    }

    m_solver->a_variable(parentId, PV->RHO) = 0.0;
    m_solver->a_variable(parentId, PV->U) = 0.0;
    m_solver->a_variable(parentId, PV->V) = 0.0;
    IF_CONSTEXPR(nDim == 3) m_solver->a_variable(parentId, PV->W) = 0.0;
    m_solver->a_oldVariable(parentId, PV->RHO) = 0.0;
    m_solver->a_oldVariable(parentId, PV->U) = 0.0;
    m_solver->a_oldVariable(parentId, PV->V) = 0.0;
    IF_CONSTEXPR(nDim == 3) m_solver->a_oldVariable(parentId, PV->W) = 0.0;
    for(MInt i = 0; i < nDist; i++) {
      m_solver->a_distribution(parentId, i) = 0.0;
      m_solver->a_oldDistribution(parentId, i) = 0.0;
    }

    // Spatial interpolation: Average all child values
    for(MInt m = 0; m < IPOW2(nDim); m++) {
      m_solver->a_variable(parentId, PV->RHO) += FFPOW2(nDim) * m_solver->a_variable(children[m], PV->RHO);
      m_solver->a_variable(parentId, PV->U) += FFPOW2(nDim) * m_solver->a_variable(children[m], PV->U);
      m_solver->a_variable(parentId, PV->V) += FFPOW2(nDim) * m_solver->a_variable(children[m], PV->V);
      IF_CONSTEXPR(nDim == 3)
      m_solver->a_variable(parentId, PV->W) += FFPOW2(nDim) * m_solver->a_variable(children[m], PV->W);
      m_solver->a_oldVariable(parentId, PV->RHO) += FFPOW2(nDim) * m_solver->a_oldVariable(children[m], PV->RHO);
      m_solver->a_oldVariable(parentId, PV->U) += FFPOW2(nDim) * m_solver->a_oldVariable(children[m], PV->U);
      m_solver->a_oldVariable(parentId, PV->V) += FFPOW2(nDim) * m_solver->a_oldVariable(children[m], PV->V);
      IF_CONSTEXPR(nDim == 3)
      m_solver->a_oldVariable(parentId, PV->W) += FFPOW2(nDim) * m_solver->a_oldVariable(children[m], PV->W);
      for(MInt i = 0; i < nDist; i++) {
        m_solver->a_distribution(parentId, i) += FFPOW2(nDim) * m_solver->a_distribution(children[m], i);
        m_solver->a_oldDistribution(parentId, i) += FFPOW2(nDim) * m_solver->a_oldDistribution(children[m], i);
      }
    }
  }
}

/**
 * \brief     Initialize parent variables from children
 * \details   - average
 *            - temporal extrapolation
 *            - Dupuis if propagation is in-sync
 * \param[in] parentId: Solver cell id of parent that will be coarsen
 *
 * \author: Philipp Brokof, Moritz Waldmann
 **/
template <MInt nDim, MInt nDist, class SysEqn>
void LbInterfaceDxQy<nDim, nDist, SysEqn>::removeChildsDupuisFilippova(const MInt parentId) {
  if(!m_solver->a_isHalo(parentId)) {
    const MInt parentLevel = m_solver->c_level(parentId);
    const MInt dParentLevel = m_solver->maxLevel() - parentLevel;

    // Relaxation paramter for parent and child refinement level
    const MFloat omega_F = F2 / (F1 + 6.0 * m_solver->a_nu(parentId) * FFPOW2(dParentLevel - 1));
    const MFloat omega_C = F2 / (F1 + 6.0 * m_solver->a_nu(parentId) * FFPOW2(dParentLevel));
    std::array<MFloat, nDist> forcing{};
    getCellForcing(forcing);

    // Get children that are to remove
    MInt children[IPOW2(nDim)];
    for(MInt m = 0; m < IPOW2(nDim); m++) {
      children[m] = m_solver->c_childId(parentId, m);
    }

    // Parent and child level have done propagation in the last timestep and will do collision after this...
    if((globalTimeStep - 1) % IPOW2(dParentLevel) == 0) {
      MFloat rho = F0;
      MFloat u[nDim] = {F0};
      MFloat rhoOld = F0;
      MFloat uOld[nDim] = {F0};
      MFloat rho_tp1 = F0;
      MFloat u_tp1[nDim] = {F0};
      MFloat rho_tm3 = F0;
      MFloat u_tm3[nDim] = {F0};

      // Spatial interpolation: Sum up all child values
      for(MInt m = 0; m < IPOW2(nDim); m++) {
        // Macroscopic variables at time t
        rho += FFPOW2(nDim) * m_solver->a_variable(children[m], PV->RHO);
        for(MInt n = 0; n < nDim; n++) {
          u[n] += FFPOW2(nDim) * m_solver->a_variable(children[m], PV->VV[n]);
        }
        // Macroscopic variables at time t-1
        rhoOld += FFPOW2(nDim) * m_solver->a_oldVariable(children[m], PV->RHO);
        for(MInt n = 0; n < nDim; n++) {
          uOld[n] += FFPOW2(nDim) * m_solver->a_oldVariable(children[m], PV->VV[n]);
        }
      }

      // Temporal extrapolation: f(t+1) = 2*f(t) - f(t-1)
      // Needed for initialization of old distribution
      rho_tp1 = 2 * rho - rhoOld;
      for(MInt d = 0; d < nDim; d++) {
        u_tp1[d] = 2 * u[d] - uOld[d];
      }

      // Temporal extrapolation: f(t-3) = 3*f(t-1) - 2*f(t)
      // Needed for initialization of macroscopic variables
      rho_tm3 = 3 * rhoOld - 2 * rho;
      for(MInt d = 0; d < nDim; d++) {
        u_tm3[d] = 3 * uOld[d] - 2 * u[d];
      }

      // Setting macroscopic variables on parent level
      m_solver->a_variable(parentId, PV->RHO) = rhoOld;
      for(MInt n = 0; n < nDim; n++) {
        m_solver->a_variable(parentId, PV->VV[n]) = uOld[n];
      }
      m_solver->a_oldVariable(parentId, PV->RHO) = rho_tm3;
      for(MInt n = 0; n < nDim; n++) {
        m_solver->a_oldVariable(parentId, PV->VV[n]) = u_tm3[n];
      }

      // Calculate equilibrium for time t+1
      std::array<MFloat, nDist> eqDistTp1{};
      eqDistTp1 = m_solver->getEqDists(rho_tp1, u_tp1);

      // Transform and interpolate post-propagation distributions according to DUPUIS
      for(MInt dist = 0; dist < nDist; dist++) {
        // Calculate mean distribution from all children
        MFloat tmpDist = 0.0;
        for(MInt m = 0; m < IPOW2(nDim); m++) {
          tmpDist += FFPOW2(nDim) * m_solver->a_oldDistribution(children[m], dist);
        }

        m_solver->a_oldDistribution(parentId, dist) =
            (tmpDist - eqDistTp1[dist]) * (omega_F / omega_C) * 2.0 + eqDistTp1[dist];

        // Cell need to perform collision step next to initialize this variable
        m_solver->a_distribution(parentId, dist) = m_solver->a_oldDistribution(parentId, dist);
      }
    } else {
      MFloat rho{};
      std::array<MFloat, nDim> u{};
      MFloat rho_tm1{};
      std::array<MFloat, nDim> u_tm1{};

      // Spatial interpolation: Sum up all child values
      for(MInt m = 0; m < IPOW2(nDim); m++) {
        // Macroscopic variables at time t
        rho += FFPOW2(nDim) * m_solver->a_variable(children[m], PV->RHO);
        for(MInt n = 0; n < nDim; n++) {
          u[n] += FFPOW2(nDim) * m_solver->a_variable(children[m], PV->VV[n]);
        }
        // Macroscopic variables at time t-1
        rho_tm1 += FFPOW2(nDim) * m_solver->a_oldVariable(children[m], PV->RHO);
        for(MInt n = 0; n < nDim; n++) {
          u_tm1[n] += FFPOW2(nDim) * m_solver->a_oldVariable(children[m], PV->VV[n]);
        }
      }

      // Temporal (backward) extrapolation f(t-2) = 2*f(t-1) - f(t)
      MFloat rho_tm2{};
      std::array<MFloat, nDim> u_tm2{};
      rho_tm2 = 2 * rho_tm1 - 1 * rho;
      for(MInt d = 0; d < nDim; d++) {
        u_tm2[d] = 2 * u_tm1[d] - 1 * u[d];
      }

      // Setting macroscopic variables for parent cell
      m_solver->a_variable(parentId, PV->RHO) = rho;
      for(MInt d = 0; d < nDim; d++) {
        m_solver->a_variable(parentId, PV->VV[d]) = u[d];
      }
      m_solver->a_oldVariable(parentId, PV->RHO) = rho_tm2;
      for(MInt d = 0; d < nDim; d++) {
        m_solver->a_oldVariable(parentId, PV->VV[d]) = u_tm2[d];
      }

      // Calculate equilibirum for time t
      std::array<MFloat, nDist> eqDist{};
      eqDist = m_solver->getEqDists(rho, u.data());

      // Transform and interpolate post-collision distributions according to FILIPPOVA
      const MFloat tau_F = 1 / omega_F;
      const MFloat tau_C = 1 / omega_C;
      for(MInt dist = 0; dist < nDist; dist++) {
        // Calculate mean distribution from all children
        MFloat tmpDist = 0.0;
        for(MInt m = 0; m < IPOW2(nDim); m++) {
          tmpDist += FFPOW2(nDim) * m_solver->a_distribution(children[m], dist);
        }

        m_solver->a_distribution(parentId, dist) =
            (tmpDist - eqDist[dist]) * ((tau_C - 1) / (tau_F - 1)) * 2.0 + eqDist[dist];

        m_solver->a_oldDistribution(parentId, dist) = m_solver->a_distribution(parentId, dist);
      }
    }
  }
}
