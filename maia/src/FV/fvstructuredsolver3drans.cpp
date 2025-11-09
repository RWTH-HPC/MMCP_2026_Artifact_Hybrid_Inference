// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "fvstructuredsolver3drans.h"
#include "globals.h"

#include <cstdlib>
#include "GRID/structuredpartition.h"
#include "fvstructuredsolverwindowinfo.h"
#if not defined(MAIA_MS_COMPILER)
#include <unistd.h>
#endif
// temporaray
#include <vector>

using namespace std;


FvStructuredSolver3DRans::FvStructuredSolver3DRans(FvStructuredSolver3D* solver)
  : m_StructuredComm(solver->m_StructuredComm),
    m_solverId(solver->m_solverId),
    m_nCells(solver->m_nCells),
    m_nPoints(solver->m_nPoints),
    m_noCells(solver->m_noCells),
    m_cells(solver->m_cells),
    CV(solver->CV),
    PV(solver->PV),
    FQ(solver->FQ),
    m_noGhostLayers(solver->m_noGhostLayers),
    m_eps(solver->m_eps),
    m_chi(solver->m_chi),
    m_sutherlandConstant(solver->m_sutherlandConstant),
    m_sutherlandPlusOne(solver->m_sutherlandPlusOne) {
  m_solver = solver;

  initFluxMethod();
}

FvStructuredSolver3DRans::~FvStructuredSolver3DRans() {}

void FvStructuredSolver3DRans::initFluxMethod() {
  // m_structuredBndryCndRans = new StructuredBndryCnd3DRans(m_solver, m_solver->m_noSpecies);
  // set pointer to the right AUSM for the RANS equations;

  /*! \property
    \page propertiesFVSTRCTRD
    \section ransMethod
    <code>MInt FvStructuredSolver3DRans::m_ransMethod </code>\n
    default = <code> 1.0 </code>\n \n
    Name of the RANS method to be used.\n
    Possible values are:\n
    <ul>
    <li>RANS_SA_DV</li>
    </ul>
    Keywords: <i>RANS, STRUCTURED</i>
  */
  m_ransMethod =
      static_cast<RansMethod>(string2enum(Context::getSolverProperty<MString>("ransMethod", m_solverId, AT_)));
  m_dsIsComputed = false;
  switch(m_ransMethod) {
    case RANS_SA_DV: {
      mAlloc(m_cells->saFlux1, nDim, m_noCells, "m_cells->saFlux1", -999999.9, AT_);
      mAlloc(m_cells->saFlux2, nDim, m_noCells, "m_cells->saFlux2", -999999.9, AT_);
      mAlloc(m_cells->prodDest, m_noCells, "m_cells->prodDest", -999999.9, AT_);
      compTurbVisc = &FvStructuredSolver3DRans::computeTurbViscosity_SA;
      viscFluxMethod = &FvStructuredSolver3DRans::viscousFlux_SA;
      switch(m_solver->CV->noVariables) {
        case 5: {
          reconstructSurfaceData = &FvStructuredSolver3DRans::Muscl_AusmDV<5>;
          break;
        }
        case 6: {
          reconstructSurfaceData = &FvStructuredSolver3DRans::Muscl_AusmDV<6>;
          break;
        }
        case 7: {
          reconstructSurfaceData = &FvStructuredSolver3DRans::Muscl_AusmDV<7>;
          break;
        }
        default: {
          stringstream errorMessage;
          errorMessage << "Number of Variables " << m_solver->CV->noVariables
                       << " not implemented! in temlate Rans AUSM " << endl;
          mTerm(1, AT_, errorMessage.str());
        }
      }
      break;
    }
    case RANS_FS: {
      viscFluxMethod = &FvStructuredSolver3DRans::viscousFlux_FS;
      compTurbVisc = &FvStructuredSolver3DRans::computeTurbViscosity_FS;

      if(m_solver->m_limiter) {
        mAlloc(m_cells->ql, PV->noVariables, m_noCells, "m_cells->ql", F0, AT_);
        mAlloc(m_cells->qr, PV->noVariables, m_noCells, "m_cells->qr", F0, AT_);
        cout << "Using RANS with Fares-Schroeder model and limited AusmDV" << endl;
        switch(m_solver->CV->noVariables) {
          case 5: {
            reconstructSurfaceData = &FvStructuredSolver3DRans::Muscl_AusmDV_Limited<5>;
            break;
          }
          case 6: {
            reconstructSurfaceData = &FvStructuredSolver3DRans::Muscl_AusmDV_Limited<6>;
            break;
          }
          case 7: {
            reconstructSurfaceData = &FvStructuredSolver3DRans::Muscl_AusmDV_Limited<7>;
            break;
          }
          default: {
            stringstream errorMessage;
            errorMessage << "Number of Variables " << m_solver->CV->noVariables
                         << " not implemented! in temlate Rans AUSM " << endl;
            mTerm(1, AT_, errorMessage.str());
          }
        }
      } else {
        cout << "Using RANS with Fares-Schroeder model and AusmDV" << endl;
        switch(m_solver->CV->noVariables) {
          case 5: {
            reconstructSurfaceData = &FvStructuredSolver3DRans::Muscl_AusmDV<5>;
            break;
          }
          case 6: {
            reconstructSurfaceData = &FvStructuredSolver3DRans::Muscl_AusmDV<6>;
            break;
          }
          case 7: {
            reconstructSurfaceData = &FvStructuredSolver3DRans::Muscl_AusmDV<7>;
            break;
          }
          default: {
            stringstream errorMessage;
            errorMessage << "Number of Variables " << m_solver->CV->noVariables
                         << " not implemented! in temlate Rans AUSM " << endl;
            mTerm(1, AT_, errorMessage.str());
          }
        }
      }
      break;
    }
    default: {
      mTerm(1, AT_, "RANS METHOD wsa not specified in properties");
      break;
    }
  }
}

void FvStructuredSolver3DRans::Muscl(MInt) { (this->*reconstructSurfaceData)(); }


void FvStructuredSolver3DRans::viscousFluxRANS() { (this->*viscFluxMethod)(); }


void FvStructuredSolver3DRans::viscousFlux_SA() {
  computeTurbViscosity_SA();

  // call the standard LES viscous flux
  m_solver->viscousFluxLES<>();

  // OTHER variables required to calculate the laminar viscous fluxes
  const MFloat rRe = F1 / m_solver->m_Re0;

  MFloat* const RESTRICT u = ALIGNED_F(&m_cells->pvariables[PV->U][0]);
  MFloat* const RESTRICT v = ALIGNED_F(&m_cells->pvariables[PV->V][0]);
  MFloat* const RESTRICT w = ALIGNED_F(&m_cells->pvariables[PV->W][0]);
  MFloat* const RESTRICT rho = ALIGNED_F(&m_cells->pvariables[PV->RHO][0]);
  MFloat* const RESTRICT nuTilde = ALIGNED_F(&m_cells->pvariables[PV->RANS_VAR[0]][0]);
  MFloat* const RESTRICT muLam = ALIGNED_F(&m_cells->fq[FQ->MU_L][0]);
  MFloat* const RESTRICT T = ALIGNED_F(&m_cells->temperature[0]);

  MFloat* const* const RESTRICT eflux = ALIGNED_F(m_cells->eFlux);
  MFloat* const* const RESTRICT fflux = ALIGNED_F(m_cells->fFlux);
  MFloat* const* const RESTRICT gflux = ALIGNED_F(m_cells->gFlux);
  MFloat* const* const RESTRICT sa_1flux = ALIGNED_F(m_cells->saFlux1);
  MFloat* const* const RESTRICT sa_2flux = ALIGNED_F(m_cells->saFlux2);

  for(MInt k = m_noGhostLayers - 1; k < m_nCells[0] - m_noGhostLayers + 1; k++) {
    for(MInt j = m_noGhostLayers - 1; j < m_nCells[1] - m_noGhostLayers + 1; j++) {
      for(MInt i = m_noGhostLayers - 1; i < m_nCells[2] - m_noGhostLayers + 1; i++) {
        // get the adjacent cells;
        const MInt IJK = cellIndex(i, j, k);
        const MInt IPJK = cellIndex((i + 1), j, k);
        const MInt IPJPK = cellIndex((i + 1), (j + 1), k);
        const MInt IJPK = cellIndex(i, (j + 1), k);
        const MInt IJKP = cellIndex(i, j, (k + 1));
        const MInt IPJKP = cellIndex((i + 1), j, (k + 1));
        const MInt IPJPKP = cellIndex((i + 1), (j + 1), (k + 1));
        const MInt IJPKP = cellIndex(i, (j + 1), (k + 1));

        const MInt IMJK = cellIndex((i - 1), j, k);
        const MInt IJMK = cellIndex(i, (j - 1), k);
        const MInt IJKM = cellIndex(i, j, (k - 1));

        const MFloat cornerMetrics[9] = {
            m_cells->cornerMetrics[0][IJK], m_cells->cornerMetrics[1][IJK], m_cells->cornerMetrics[2][IJK],
            m_cells->cornerMetrics[3][IJK], m_cells->cornerMetrics[4][IJK], m_cells->cornerMetrics[5][IJK],
            m_cells->cornerMetrics[6][IJK], m_cells->cornerMetrics[7][IJK], m_cells->cornerMetrics[8][IJK]};


        const MFloat dnutldxi = F1B4
                                * (nuTilde[IPJPKP] + nuTilde[IPJPK] + nuTilde[IPJKP] + nuTilde[IPJK] - nuTilde[IJPKP]
                                   - nuTilde[IJPK] - nuTilde[IJKP] - nuTilde[IJK]);
        const MFloat dnutldet = F1B4
                                * (nuTilde[IPJPKP] + nuTilde[IJPKP] + nuTilde[IPJPK] + nuTilde[IJPK] - nuTilde[IPJKP]
                                   - nuTilde[IJKP] - nuTilde[IPJK] - nuTilde[IJK]);
        const MFloat dnutldze = F1B4
                                * (nuTilde[IPJPKP] + nuTilde[IJPKP] + nuTilde[IPJKP] + nuTilde[IJKP] - nuTilde[IPJPK]
                                   - nuTilde[IJPK] - nuTilde[IPJK] - nuTilde[IJK]);


        const MFloat nutldAvg = F1B8
                                * (nuTilde[IPJPKP] + nuTilde[IJPKP] + nuTilde[IJPK] + nuTilde[IPJPK] + nuTilde[IPJKP]
                                   + nuTilde[IJKP] + nuTilde[IJK] + nuTilde[IPJK]);

        const MFloat nuLamAvg = F1B8
                                * (muLam[IPJPKP] / rho[IPJPKP] + muLam[IJPKP] / rho[IJPKP] + muLam[IJPK] / rho[IJPK]
                                   + muLam[IPJPK] / rho[IPJPK] + muLam[IPJKP] / rho[IPJKP] + muLam[IJKP] / rho[IJKP]
                                   + muLam[IJK] / rho[IJK] + muLam[IPJK] / rho[IPJK]);

        const MFloat dnutldx = dnutldxi * cornerMetrics[xsd * 3 + xsd] + dnutldet * cornerMetrics[ysd * 3 + xsd]
                               + dnutldze * cornerMetrics[zsd * 3 + xsd];

        const MFloat dnutldy = dnutldxi * cornerMetrics[xsd * 3 + ysd] + dnutldet * cornerMetrics[ysd * 3 + ysd]
                               + dnutldze * cornerMetrics[zsd * 3 + ysd];

        const MFloat dnutldz = dnutldxi * cornerMetrics[xsd * 3 + zsd] + dnutldet * cornerMetrics[ysd * 3 + zsd]
                               + dnutldze * cornerMetrics[zsd * 3 + zsd];

        const MFloat Frj = rRe / m_cells->cornerJac[IJK];

        const MFloat sax1 = Frj * (nuLamAvg + (1.0 + RM_SA_DV::cb2) * nutldAvg)
                            * (dnutldx * cornerMetrics[xsd * 3 + xsd] + dnutldy * cornerMetrics[xsd * 3 + ysd]
                               + dnutldz * cornerMetrics[xsd * 3 + zsd]);

        const MFloat say1 = Frj * (nuLamAvg + (1.0 + RM_SA_DV::cb2) * nutldAvg)
                            * (dnutldx * cornerMetrics[ysd * 3 + xsd] + dnutldy * cornerMetrics[ysd * 3 + ysd]
                               + dnutldz * cornerMetrics[ysd * 3 + zsd]);

        const MFloat saz1 = Frj * (nuLamAvg + (1.0 + RM_SA_DV::cb2) * nutldAvg)
                            * (dnutldx * cornerMetrics[zsd * 3 + xsd] + dnutldy * cornerMetrics[zsd * 3 + ysd]
                               + dnutldz * cornerMetrics[zsd * 3 + zsd]);

        const MFloat sax2 = -Frj * RM_SA_DV::cb2
                            * (dnutldx * cornerMetrics[xsd * 3 + xsd] + dnutldy * cornerMetrics[xsd * 3 + ysd]
                               + dnutldz * cornerMetrics[xsd * 3 + zsd]);
        const MFloat say2 = -Frj * RM_SA_DV::cb2
                            * (dnutldx * cornerMetrics[ysd * 3 + xsd] + dnutldy * cornerMetrics[ysd * 3 + ysd]
                               + dnutldz * cornerMetrics[ysd * 3 + zsd]);

        const MFloat saz2 = -Frj * RM_SA_DV::cb2
                            * (dnutldx * cornerMetrics[zsd * 3 + xsd] + dnutldy * cornerMetrics[zsd * 3 + ysd]
                               + dnutldz * cornerMetrics[zsd * 3 + zsd]);

        // for dwdy
        const MFloat dw1 = w[IPJK] - w[IMJK];
        const MFloat dw2 = w[IJPK] - w[IJMK];
        const MFloat dw3 = w[IJKP] - w[IJKM];

        // for dvdz
        const MFloat dv1 = v[IPJK] - v[IMJK];
        const MFloat dv2 = v[IJPK] - v[IJMK];
        const MFloat dv3 = v[IJKP] - v[IJKM];

        // for dudz
        const MFloat du1 = u[IPJK] - u[IMJK];
        const MFloat du2 = u[IJPK] - u[IJMK];
        const MFloat du3 = u[IJKP] - u[IJKM];

        const MFloat vorti =
            (m_cells->cellMetrics[xsd * 3 + ysd][IJK] * dw1) + (m_cells->cellMetrics[ysd * 3 + ysd][IJK] * dw2)
            + (m_cells->cellMetrics[zsd * 3 + ysd][IJK] * dw3) -

            (m_cells->cellMetrics[xsd * 3 + zsd][IJK] * dv1) - (m_cells->cellMetrics[ysd * 3 + zsd][IJK] * dv2)
            - (m_cells->cellMetrics[zsd * 3 + zsd][IJK] * dv3);

        const MFloat vortj =
            (m_cells->cellMetrics[xsd * 3 + zsd][IJK] * du1) + (m_cells->cellMetrics[ysd * 3 + zsd][IJK] * du2)
            + (m_cells->cellMetrics[zsd * 3 + zsd][IJK] * du3) -

            (m_cells->cellMetrics[xsd * 3 + xsd][IJK] * dw1) - (m_cells->cellMetrics[ysd * 3 + xsd][IJK] * dw2)
            - (m_cells->cellMetrics[zsd * 3 + xsd][IJK] * dw3);

        const MFloat vortk =
            (m_cells->cellMetrics[xsd * 3 + xsd][IJK] * dv1) + (m_cells->cellMetrics[ysd * 3 + xsd][IJK] * dv2)
            + (m_cells->cellMetrics[zsd * 3 + xsd][IJK] * dv3) -

            (m_cells->cellMetrics[xsd * 3 + ysd][IJK] * du1) - (m_cells->cellMetrics[ysd * 3 + ysd][IJK] * du2)
            - (m_cells->cellMetrics[zsd * 3 + ysd][IJK] * du3);

        MFloat s = (vorti * vorti) + (vortj * vortj) + (vortk * vortk);
        s = F1B2 * sqrt(s) / m_cells->cellJac[IJK];

        // assuming wall distance function
        const MFloat distance = m_cells->fq[FQ->WALLDISTANCE][IJK];
        const MFloat Fdist2 = 1.0 / (distance * distance);
        const MFloat chi = nuTilde[IJK] * rho[IJK] / (SUTHERLANDLAW(T[IJK]) / rho[IJK]);
        const MFloat chip3 = chi * chi * chi;
        const MFloat Fv1 = chip3 / (chip3 + RM_SA_DV::cv1to3);
        const MFloat Fv2 = F1 - (chi / (F1 + chi * Fv1));

        const MFloat term = nuTilde[IJK] * Fdist2 * RM_SA_DV::Fkap2;
        const MFloat stilde = s + term * Fv2 * rRe;
        const MFloat r = min(10.0, rRe * term / stilde);

        const MFloat g = r + RM_SA_DV::cw2 * (pow(r, 6) - r);
        const MFloat Fwterm = (1 + RM_SA_DV::cw3to6) / (pow(g, 6) + RM_SA_DV::cw3to6);
        const MFloat Fw = g * pow(Fwterm, (1.0 / 6.0));
        const MFloat prodValue = rho[IJK] * RM_SA_DV::cb1 * (F1 - RM_SA_DV::Ft2) * stilde * nuTilde[IJK];
        const MFloat destValue = rRe * rho[IJK] * (RM_SA_DV::cw1 * Fw - RM_SA_DV::cb1 * RM_SA_DV::Fkap2 * RM_SA_DV::Ft2)
                                 * pow(nuTilde[IJK], 2.0) * Fdist2;

        m_cells->prodDest[IJK] = (prodValue - destValue) * m_cells->cellJac[IJK];

        eflux[0][IJK] = sax1; // diffusion of nutilde for every cell
        eflux[1][IJK] = sax2; // diffusion of nutilde for every cell

        fflux[0][IJK] = say1; // diffusion of nutilde for every cell
        fflux[1][IJK] = say2; // diffusion of nutilde for every cell

        gflux[0][IJK] = saz1; // diffusion of nutilde for every cell
        gflux[1][IJK] = saz2; // diffusion of nutilde for every cell
      }
    }
  }

  for(MInt k = m_noGhostLayers; k < m_nCells[0] - m_noGhostLayers; k++) {
    for(MInt j = m_noGhostLayers; j < m_nCells[1] - m_noGhostLayers; j++) {
      for(MInt i = m_noGhostLayers - 1; i < m_nCells[2] - m_noGhostLayers; i++) {
        const MInt IJK = cellIndex(i, j, k);
        const MInt IJMK = cellIndex(i, (j - 1), k);
        const MInt IJKM = cellIndex(i, j, (k - 1));
        const MInt IJMKM = cellIndex(i, (j - 1), (k - 1));

        sa_1flux[0][IJK] = F1B4 * (eflux[0][IJK] + eflux[0][IJKM] + eflux[0][IJMK] + eflux[0][IJMKM]);
        sa_2flux[0][IJK] = F1B4 * (eflux[1][IJK] + eflux[1][IJKM] + eflux[1][IJMK] + eflux[1][IJMKM]);
      }
    }
  }


  for(MInt k = m_noGhostLayers; k < m_nCells[0] - m_noGhostLayers; k++) {
    for(MInt j = m_noGhostLayers - 1; j < m_nCells[1] - m_noGhostLayers; j++) {
      for(MInt i = m_noGhostLayers; i < m_nCells[2] - m_noGhostLayers; i++) {
        const MInt IJK = cellIndex(i, j, k);
        const MInt IMJK = cellIndex((i - 1), j, k);
        const MInt IJKM = cellIndex(i, j, (k - 1));
        const MInt IMJKM = cellIndex((i - 1), j, (k - 1));

        sa_1flux[1][IJK] = F1B4 * (fflux[0][IJK] + fflux[0][IJKM] + fflux[0][IMJK] + fflux[0][IMJKM]);
        sa_2flux[1][IJK] = F1B4 * (fflux[1][IJK] + fflux[1][IJKM] + fflux[1][IMJK] + fflux[1][IMJKM]);
      }
    }
  }

  for(MInt k = m_noGhostLayers - 1; k < m_nCells[0] - m_noGhostLayers; k++) {
    for(MInt j = m_noGhostLayers; j < m_nCells[1] - m_noGhostLayers; j++) {
      for(MInt i = m_noGhostLayers; i < m_nCells[2] - m_noGhostLayers; i++) {
        const MInt IJK = cellIndex(i, j, k);
        const MInt IMJK = cellIndex((i - 1), j, k);
        const MInt IJMK = cellIndex(i, (j - 1), k);
        const MInt IMJMK = cellIndex((i - 1), (j - 1), k);

        sa_1flux[2][IJK] = F1B4 * (gflux[0][IJK] + gflux[0][IMJK] + gflux[0][IJMK] + gflux[0][IMJMK]);
        sa_2flux[2][IJK] = F1B4 * (gflux[1][IJK] + gflux[1][IMJK] + gflux[1][IJMK] + gflux[1][IMJMK]);
      }
    }
  }

  // separate loop for adding the prodn nad destrn terms for tur kin viscosity transport variable
  for(MInt k = m_noGhostLayers; k < m_nCells[0] - m_noGhostLayers; k++) {
    for(MInt j = m_noGhostLayers; j < m_nCells[1] - m_noGhostLayers; j++) {
      for(MInt i = m_noGhostLayers; i < m_nCells[2] - m_noGhostLayers; i++) {
        const MInt IJK = cellIndex(i, j, k);
        const MInt IMJK = cellIndex(i - 1, j, k);
        const MInt IJMK = cellIndex(i, j - 1, k);
        const MInt IJKM = cellIndex(i, j, k - 1);
        const MFloat dissipation_term =
            (((sa_1flux[0][IJK] - sa_1flux[0][IMJK]) + ((sa_2flux[0][IJK] - sa_2flux[0][IMJK]) * nuTilde[IJK]))
                 * rho[IJK] * RM_SA_DV::Fsigma
             + ((sa_1flux[1][IJK] - sa_1flux[1][IJMK]) + ((sa_2flux[1][IJK] - sa_2flux[1][IJMK]) * nuTilde[IJK]))
                   * rho[IJK] * RM_SA_DV::Fsigma
             + ((sa_1flux[2][IJK] - sa_1flux[2][IJKM]) + ((sa_2flux[2][IJK] - sa_2flux[2][IJKM]) * nuTilde[IJK]))
                   * rho[IJK] * RM_SA_DV::Fsigma);

        m_cells->rightHandSide[CV->RANS_FIRST][IJK] += dissipation_term;
        m_cells->rightHandSide[CV->RANS_FIRST][IJK] += m_cells->prodDest[IJK];
      }
    }
  }
}


void FvStructuredSolver3DRans::viscousFlux_FS() {
  const MFloat eps = 1e-16;
  const MFloat epss = 1e-34;

  computeTurbViscosity_FS();

  // call the standard LES viscous flux
  m_solver->viscousFluxLES<>();

  // OTHER variables required to calculate the laminar viscous fluxes
  const MInt noCells = m_noCells;
  const MFloat rRe = F1 / m_solver->m_Re0;

  MFloat* const RESTRICT u = ALIGNED_F(&m_cells->pvariables[PV->U][0]);
  MFloat* const RESTRICT v = ALIGNED_F(&m_cells->pvariables[PV->V][0]);
  MFloat* const RESTRICT w = ALIGNED_F(&m_cells->pvariables[PV->W][0]);
  MFloat* const RESTRICT rho = ALIGNED_F(&m_cells->pvariables[PV->RHO][0]);
  MFloat* const RESTRICT nuTilde = ALIGNED_F(&m_cells->pvariables[PV->RANS_VAR[0]][0]);
  MFloat* const RESTRICT muLam = ALIGNED_F(&m_cells->fq[FQ->MU_L][0]);

  MFloat* const* const RESTRICT eflux = ALIGNED_F(m_cells->eFlux);
  MFloat* const* const RESTRICT fflux = ALIGNED_F(m_cells->fFlux);
  MFloat* const* const RESTRICT gflux = ALIGNED_F(m_cells->gFlux);

  MFloatScratchSpace uvwn(noCells, 4, nDim, AT_, "uvwn");

  for(MInt k = m_noGhostLayers; k < m_nCells[0] - m_noGhostLayers; k++) {
    for(MInt j = m_noGhostLayers; j < m_nCells[1] - m_noGhostLayers; j++) {
      for(MInt i = m_noGhostLayers; i < m_nCells[2] - m_noGhostLayers; i++) {
        // get the adjacent cells;
        const MInt IJK = cellIndex(i, j, k);
        const MInt IPJK = cellIndex((i + 1), j, k);
        const MInt IJPK = cellIndex(i, (j + 1), k);
        const MInt IJKP = cellIndex(i, j, (k + 1));

        const MInt IMJK = cellIndex((i - 1), j, k);
        const MInt IJMK = cellIndex(i, (j - 1), k);
        const MInt IJKM = cellIndex(i, j, (k - 1));

        const MFloat cellMetrics[9] = {
            m_cells->cellMetrics[0][IJK], m_cells->cellMetrics[1][IJK], m_cells->cellMetrics[2][IJK],
            m_cells->cellMetrics[3][IJK], m_cells->cellMetrics[4][IJK], m_cells->cellMetrics[5][IJK],
            m_cells->cellMetrics[6][IJK], m_cells->cellMetrics[7][IJK], m_cells->cellMetrics[8][IJK]};

        const MFloat fjac = F1B2 / m_cells->cellJac[IJK];

        const MFloat dudxi = fjac * (u[IPJK] - u[IMJK]);
        const MFloat dudet = fjac * (u[IJPK] - u[IJMK]);
        const MFloat dudze = fjac * (u[IJKP] - u[IJKM]);

        const MFloat dvdxi = fjac * (v[IPJK] - v[IMJK]);
        const MFloat dvdet = fjac * (v[IJPK] - v[IJMK]);
        const MFloat dvdze = fjac * (v[IJKP] - v[IJKM]);

        const MFloat dwdxi = fjac * (w[IPJK] - w[IMJK]);
        const MFloat dwdet = fjac * (w[IJPK] - w[IJMK]);
        const MFloat dwdze = fjac * (w[IJKP] - w[IJKM]);

        const MFloat dnutdxi = fjac * (nuTilde[IPJK] - nuTilde[IMJK]);
        const MFloat dnutdet = fjac * (nuTilde[IJPK] - nuTilde[IJMK]);
        const MFloat dnutdze = fjac * (nuTilde[IJKP] - nuTilde[IJKM]);


        uvwn(IJK, 0, 0) = cellMetrics[0] * dudxi + cellMetrics[3] * dudet + cellMetrics[6] * dudze;
        uvwn(IJK, 0, 1) = cellMetrics[1] * dudxi + cellMetrics[4] * dudet + cellMetrics[7] * dudze;
        uvwn(IJK, 0, 2) = cellMetrics[2] * dudxi + cellMetrics[5] * dudet + cellMetrics[8] * dudze;

        uvwn(IJK, 1, 0) = cellMetrics[0] * dvdxi + cellMetrics[3] * dvdet + cellMetrics[6] * dvdze;
        uvwn(IJK, 1, 1) = cellMetrics[1] * dvdxi + cellMetrics[4] * dvdet + cellMetrics[7] * dvdze;
        uvwn(IJK, 1, 2) = cellMetrics[2] * dvdxi + cellMetrics[5] * dvdet + cellMetrics[8] * dvdze;

        uvwn(IJK, 2, 0) = cellMetrics[0] * dwdxi + cellMetrics[3] * dwdet + cellMetrics[6] * dwdze;
        uvwn(IJK, 2, 1) = cellMetrics[1] * dwdxi + cellMetrics[4] * dwdet + cellMetrics[7] * dwdze;
        uvwn(IJK, 2, 2) = cellMetrics[2] * dwdxi + cellMetrics[5] * dwdet + cellMetrics[8] * dwdze;

        uvwn(IJK, 3, 0) = cellMetrics[0] * dnutdxi + cellMetrics[3] * dnutdet + cellMetrics[6] * dnutdze;
        uvwn(IJK, 3, 1) = cellMetrics[1] * dnutdxi + cellMetrics[4] * dnutdet + cellMetrics[7] * dnutdze;
        uvwn(IJK, 3, 2) = cellMetrics[2] * dnutdxi + cellMetrics[5] * dnutdet + cellMetrics[8] * dnutdze;
      }
    }
  }


  for(MInt var = 0; var < 4; var++) {
    // i start
    for(MInt k = m_noGhostLayers - 1; k < m_nCells[0] - m_noGhostLayers + 1; k++) {
      for(MInt j = m_noGhostLayers - 1; j < m_nCells[1] - m_noGhostLayers + 1; j++) {
        const MInt II1 = cellIndex(m_noGhostLayers - 1, j, k);
        const MInt II2 = cellIndex(m_noGhostLayers, j, k);
        const MInt II3 = cellIndex(m_noGhostLayers + 1, j, k);

        uvwn(II1, var, 0) = 2.0 * uvwn(II2, var, 0) - uvwn(II3, var, 0);
        uvwn(II1, var, 1) = 2.0 * uvwn(II2, var, 1) - uvwn(II3, var, 1);
        uvwn(II1, var, 2) = 2.0 * uvwn(II2, var, 2) - uvwn(II3, var, 2);
      }
    }

    // i end
    for(MInt k = m_noGhostLayers - 1; k < m_nCells[0] - m_noGhostLayers + 1; k++) {
      for(MInt j = m_noGhostLayers - 1; j < m_nCells[1] - m_noGhostLayers + 1; j++) {
        const MInt II1 = cellIndex(m_nCells[2] - m_noGhostLayers + 1, j, k);
        const MInt II2 = cellIndex(m_nCells[2] - m_noGhostLayers, j, k);
        const MInt II3 = cellIndex(m_nCells[2] - m_noGhostLayers - 1, j, k);

        uvwn(II1, var, 0) = 2.0 * uvwn(II2, var, 0) - uvwn(II3, var, 0);
        uvwn(II1, var, 1) = 2.0 * uvwn(II2, var, 1) - uvwn(II3, var, 1);
        uvwn(II1, var, 2) = 2.0 * uvwn(II2, var, 2) - uvwn(II3, var, 2);
      }
    }

    // j start
    for(MInt k = m_noGhostLayers - 1; k < m_nCells[0] - m_noGhostLayers + 1; k++) {
      for(MInt i = m_noGhostLayers - 1; i < m_nCells[2] - m_noGhostLayers + 1; i++) {
        const MInt II1 = cellIndex(i, m_noGhostLayers - 1, k);
        const MInt II2 = cellIndex(i, m_noGhostLayers, k);
        const MInt II3 = cellIndex(i, m_noGhostLayers + 1, k);

        uvwn(II1, var, 0) = 2.0 * uvwn(II2, var, 0) - uvwn(II3, var, 0);
        uvwn(II1, var, 1) = 2.0 * uvwn(II2, var, 1) - uvwn(II3, var, 1);
        uvwn(II1, var, 2) = 2.0 * uvwn(II2, var, 2) - uvwn(II3, var, 2);
      }
    }

    // j end
    for(MInt k = m_noGhostLayers - 1; k < m_nCells[0] - m_noGhostLayers + 1; k++) {
      for(MInt i = m_noGhostLayers - 1; i < m_nCells[2] - m_noGhostLayers + 1; i++) {
        const MInt II1 = cellIndex(i, m_nCells[1] - m_noGhostLayers + 1, k);
        const MInt II2 = cellIndex(i, m_nCells[1] - m_noGhostLayers, k);
        const MInt II3 = cellIndex(i, m_nCells[1] - m_noGhostLayers - 1, k);

        uvwn(II1, var, 0) = 2.0 * uvwn(II2, var, 0) - uvwn(II3, var, 0);
        uvwn(II1, var, 1) = 2.0 * uvwn(II2, var, 1) - uvwn(II3, var, 1);
        uvwn(II1, var, 2) = 2.0 * uvwn(II2, var, 2) - uvwn(II3, var, 2);
      }
    }


    // k start
    for(MInt j = m_noGhostLayers - 1; j < m_nCells[1] - m_noGhostLayers + 1; j++) {
      for(MInt i = m_noGhostLayers - 1; i < m_nCells[2] - m_noGhostLayers + 1; i++) {
        const MInt II1 = cellIndex(i, j, m_noGhostLayers - 1);
        const MInt II2 = cellIndex(i, j, m_noGhostLayers);
        const MInt II3 = cellIndex(i, j, m_noGhostLayers + 1);

        uvwn(II1, var, 0) = 2.0 * uvwn(II2, var, 0) - uvwn(II3, var, 0);
        uvwn(II1, var, 1) = 2.0 * uvwn(II2, var, 1) - uvwn(II3, var, 1);
        uvwn(II1, var, 2) = 2.0 * uvwn(II2, var, 2) - uvwn(II3, var, 2);
      }
    }

    // k end
    for(MInt j = m_noGhostLayers - 1; j < m_nCells[1] - m_noGhostLayers + 1; j++) {
      for(MInt i = m_noGhostLayers - 1; i < m_nCells[2] - m_noGhostLayers + 1; i++) {
        const MInt II1 = cellIndex(i, j, m_nCells[0] - m_noGhostLayers + 1);
        const MInt II2 = cellIndex(i, j, m_nCells[0] - m_noGhostLayers);
        const MInt II3 = cellIndex(i, j, m_nCells[0] - m_noGhostLayers - 1);

        uvwn(II1, var, 0) = 2.0 * uvwn(II2, var, 0) - uvwn(II3, var, 0);
        uvwn(II1, var, 1) = 2.0 * uvwn(II2, var, 1) - uvwn(II3, var, 1);
        uvwn(II1, var, 2) = 2.0 * uvwn(II2, var, 2) - uvwn(II3, var, 2);
      }
    }
  }

  MFloatScratchSpace omega(noCells, AT_, "omega");
  MFloatScratchSpace dumm(noCells, AT_, "dumm");
  omega.fill(0.0);
  dumm.fill(0.0);

  for(MInt j = 0; j < nDim; j++) {
    for(MInt i = 0; i < nDim; i++) {
      for(MInt kk = m_noGhostLayers - 1; kk < m_nCells[0] - m_noGhostLayers + 1; kk++) {
        for(MInt jj = m_noGhostLayers - 1; jj < m_nCells[1] - m_noGhostLayers + 1; jj++) {
          for(MInt ii = m_noGhostLayers - 1; ii < m_nCells[2] - m_noGhostLayers + 1; ii++) {
            const MInt IJK = cellIndex(ii, jj, kk);
            const MFloat Sij = F1B2 * (uvwn(IJK, i, j) + uvwn(IJK, j, i));
            dumm(IJK) += Sij * Sij;
          }
        }
      }
    }
  }

  MFloat fac = 1.0 / sqrt(RM_FS::fabetcs);

  for(MInt kk = m_noGhostLayers - 1; kk < m_nCells[0] - m_noGhostLayers + 1; kk++) {
    for(MInt jj = m_noGhostLayers - 1; jj < m_nCells[1] - m_noGhostLayers + 1; jj++) {
      for(MInt ii = m_noGhostLayers - 1; ii < m_nCells[2] - m_noGhostLayers + 1; ii++) {
        const MInt IJK = cellIndex(ii, jj, kk);
        omega(IJK) = mMax(eps, fac * sqrt(2.0 * dumm(IJK)));
      }
    }
  }

  dumm.fill(0.0);

  for(MInt k = 0; k < nDim; k++) {
    for(MInt j = 0; j < nDim; j++) {
      for(MInt i = 0; i < nDim; i++) {
        for(MInt kk = m_noGhostLayers; kk < m_nCells[0] - m_noGhostLayers; kk++) {
          for(MInt jj = m_noGhostLayers; jj < m_nCells[1] - m_noGhostLayers; jj++) {
            for(MInt ii = m_noGhostLayers; ii < m_nCells[2] - m_noGhostLayers; ii++) {
              const MInt IJK = cellIndex(ii, jj, kk);
              const MFloat Oij = 0.5 * (uvwn(IJK, i, j) - uvwn(IJK, j, i));
              const MFloat Ojk = 0.5 * (uvwn(IJK, j, k) - uvwn(IJK, k, j));
              const MFloat Ski = 0.5 * (uvwn(IJK, k, i) + uvwn(IJK, i, k));
              dumm(IJK) += Oij * Ojk * Ski;
            }
          }
        }
      }
    }
  }

  fac = 1.0 / pow(RM_FS::fabetcs, 3.0);

  for(MInt k = m_noGhostLayers; k < m_nCells[0] - m_noGhostLayers; k++) {
    for(MInt j = m_noGhostLayers; j < m_nCells[1] - m_noGhostLayers; j++) {
      for(MInt i = m_noGhostLayers; i < m_nCells[2] - m_noGhostLayers; i++) {
        const MInt IJK = cellIndex(i, j, k);
        const MInt IPJK = cellIndex((i + 1), j, k);
        const MInt IJPK = cellIndex(i, (j + 1), k);
        const MInt IJKP = cellIndex(i, j, (k + 1));

        const MInt IMJK = cellIndex((i - 1), j, k);
        const MInt IJMK = cellIndex(i, (j - 1), k);
        const MInt IJKM = cellIndex(i, j, (k - 1));

        const MFloat fjac = 0.5 / m_cells->cellJac[IJK];

        const MFloat cellMetrics[9] = {
            m_cells->cellMetrics[0][IJK], m_cells->cellMetrics[1][IJK], m_cells->cellMetrics[2][IJK],
            m_cells->cellMetrics[3][IJK], m_cells->cellMetrics[4][IJK], m_cells->cellMetrics[5][IJK],
            m_cells->cellMetrics[6][IJK], m_cells->cellMetrics[7][IJK], m_cells->cellMetrics[8][IJK]};

        const MFloat ddxi = fjac * (omega(IPJK) - omega(IMJK));
        const MFloat ddeta = fjac * (omega(IJPK) - omega(IJMK));
        const MFloat ddzeta = fjac * (omega(IJKP) - omega(IJKM));

        const MFloat domdx = cellMetrics[0] * ddxi + cellMetrics[3] * ddeta + cellMetrics[6] * ddzeta;
        const MFloat domdy = cellMetrics[1] * ddxi + cellMetrics[4] * ddeta + cellMetrics[7] * ddzeta;
        const MFloat domdz = cellMetrics[2] * ddxi + cellMetrics[5] * ddeta + cellMetrics[8] * ddzeta;


        const MFloat dudx = uvwn(IJK, 0, 0);
        const MFloat dudy = uvwn(IJK, 0, 1);
        const MFloat dudz = uvwn(IJK, 0, 2);

        const MFloat dvdx = uvwn(IJK, 1, 0);
        const MFloat dvdy = uvwn(IJK, 1, 1);
        const MFloat dvdz = uvwn(IJK, 1, 2);

        const MFloat dwdx = uvwn(IJK, 2, 0);
        const MFloat dwdy = uvwn(IJK, 2, 1);
        const MFloat dwdz = uvwn(IJK, 2, 2);

        const MFloat dndx = uvwn(IJK, 3, 0);
        const MFloat dndy = uvwn(IJK, 3, 1);
        const MFloat dndz = uvwn(IJK, 3, 2);

        const MFloat CDNOM = dndx * domdx + dndy * domdy + dndz * domdz;
        const MFloat fpsik = 0.0;
        const MFloat crdif = 0.0;

        const MFloat betas = RM_FS::fabetcs * (1.0 + RM_FS::fapsik1 * fpsik) / (1.0 + RM_FS::fapsik2 * fpsik);
        const MFloat beta = RM_FS::fabetc;

        const MFloat P =
            (dudy * (dudy + dvdx) + dudz * (dudz + dwdx) + dvdx * (dvdx + dudy) + dvdz * (dvdz + dwdy)
             + dwdx * (dwdx + dudz) + dwdy * (dwdy + dvdz) + 2.0 * dudx * dudx + 2.0 * dvdy * dvdy + 2.0 * dwdz * dwdz);

        const MFloat fv2t = 1.0;
        const MFloat prod1 = (1.0 - RM_FS::faalpha * fv2t) * nuTilde[IJK] * rho[IJK] / omega(IJK) * P;
        MFloat prod = prod1 - (1.0 - RM_FS::faalpha * fv2t) * F2B3 * nuTilde[IJK] * (dudx + dvdy + dwdz) * rho[IJK];
        const MFloat dest = (betas - beta) * nuTilde[IJK] * omega(IJK) * rho[IJK] + rho[IJK] * crdif * 0.125;
        const MFloat prodwall = mMin(
            0.0, 2.0 * rRe * rho[IJK] / omega(IJK) * (muLam[IJK] / rho[IJK] + RM_FS::fasigma * nuTilde[IJK]) * CDNOM);

        prod += prodwall;

        m_cells->rightHandSide[CV->RANS_FIRST][IJK] += m_cells->cellJac[IJK] * (prod - dest);
      }
    }
  }


  for(MInt k = m_noGhostLayers - 1; k < m_nCells[0] - m_noGhostLayers + 1; k++) {
    for(MInt j = m_noGhostLayers - 1; j < m_nCells[1] - m_noGhostLayers + 1; j++) {
      for(MInt i = m_noGhostLayers - 1; i < m_nCells[2] - m_noGhostLayers + 1; i++) {
        const MInt IJK = cellIndex(i, j, k);
        const MInt IPJK = cellIndex((i + 1), j, k);
        const MInt IPJPK = cellIndex((i + 1), (j + 1), k);
        const MInt IJPK = cellIndex(i, (j + 1), k);
        const MInt IJKP = cellIndex(i, j, (k + 1));
        const MInt IPJKP = cellIndex((i + 1), j, (k + 1));
        const MInt IPJPKP = cellIndex((i + 1), (j + 1), (k + 1));
        const MInt IJPKP = cellIndex(i, (j + 1), (k + 1));

        const MFloat cornerMetrics[9] = {
            m_cells->cornerMetrics[0][IJK], m_cells->cornerMetrics[1][IJK], m_cells->cornerMetrics[2][IJK],
            m_cells->cornerMetrics[3][IJK], m_cells->cornerMetrics[4][IJK], m_cells->cornerMetrics[5][IJK],
            m_cells->cornerMetrics[6][IJK], m_cells->cornerMetrics[7][IJK], m_cells->cornerMetrics[8][IJK]};

        const MFloat nutldAvg = F1B8
                                * (nuTilde[IPJPKP] + nuTilde[IJPKP] + nuTilde[IJPK] + nuTilde[IPJPK] + nuTilde[IPJKP]
                                   + nuTilde[IJKP] + nuTilde[IJK] + nuTilde[IPJK]);

        const MFloat nuLamAvg = F1B8
                                * (muLam[IPJPKP] / rho[IPJPKP] + muLam[IJPKP] / rho[IJPKP] + muLam[IJPK] / rho[IJPK]
                                   + muLam[IPJPK] / rho[IPJPK] + muLam[IPJKP] / rho[IPJKP] + muLam[IJKP] / rho[IJKP]
                                   + muLam[IJK] / rho[IJK] + muLam[IPJK] / rho[IPJK]);

        const MFloat dnutldxi = F1B4
                                * (nuTilde[IPJPKP] + nuTilde[IPJPK] + nuTilde[IPJKP] + nuTilde[IPJK] - nuTilde[IJPKP]
                                   - nuTilde[IJPK] - nuTilde[IJKP] - nuTilde[IJK]);
        const MFloat dnutldet = F1B4
                                * (nuTilde[IPJPKP] + nuTilde[IJPKP] + nuTilde[IPJPK] + nuTilde[IJPK] - nuTilde[IPJKP]
                                   - nuTilde[IJKP] - nuTilde[IPJK] - nuTilde[IJK]);
        const MFloat dnutldze = F1B4
                                * (nuTilde[IPJPKP] + nuTilde[IJPKP] + nuTilde[IPJKP] + nuTilde[IJKP] - nuTilde[IPJPK]
                                   - nuTilde[IJPK] - nuTilde[IPJK] - nuTilde[IJK]);

        const MFloat dnutldx = (dnutldxi * cornerMetrics[xsd * 3 + xsd] + dnutldet * cornerMetrics[ysd * 3 + xsd]
                                + dnutldze * cornerMetrics[zsd * 3 + xsd]);

        const MFloat dnutldy = (dnutldxi * cornerMetrics[xsd * 3 + ysd] + dnutldet * cornerMetrics[ysd * 3 + ysd]
                                + dnutldze * cornerMetrics[zsd * 3 + ysd]);

        const MFloat dnutldz = (dnutldxi * cornerMetrics[xsd * 3 + zsd] + dnutldet * cornerMetrics[ysd * 3 + zsd]
                                + dnutldze * cornerMetrics[zsd * 3 + zsd]);

        const MFloat Frj = rRe / mMax(fabs(m_cells->cornerJac[IJK]), epss);

        eflux[6][IJK] = (Frj * (nuLamAvg + RM_FS::fasigma * nutldAvg)
                         * (dnutldx * cornerMetrics[xsd * 3 + xsd] + dnutldy * cornerMetrics[xsd * 3 + ysd]
                            + dnutldz * cornerMetrics[xsd * 3 + zsd]));

        fflux[6][IJK] = (Frj * (nuLamAvg + RM_FS::fasigma * nutldAvg)
                         * (dnutldx * cornerMetrics[ysd * 3 + xsd] + dnutldy * cornerMetrics[ysd * 3 + ysd]
                            + dnutldz * cornerMetrics[ysd * 3 + zsd]));

        gflux[6][IJK] = (Frj * (nuLamAvg + RM_FS::fasigma * nutldAvg)
                         * (dnutldx * cornerMetrics[zsd * 3 + xsd] + dnutldy * cornerMetrics[zsd * 3 + ysd]
                            + dnutldz * cornerMetrics[zsd * 3 + zsd]));
      }
    }
  }

  for(MInt k = m_noGhostLayers; k < m_nCells[0] - m_noGhostLayers; k++) {
    for(MInt j = m_noGhostLayers; j < m_nCells[1] - m_noGhostLayers; j++) {
      for(MInt i = m_noGhostLayers - 1; i < m_nCells[2] - m_noGhostLayers; i++) {
        const MInt IJK = cellIndex(i, j, k);
        const MInt IJMK = cellIndex(i, (j - 1), k);
        const MInt IJKM = cellIndex(i, j, (k - 1));
        const MInt IJMKM = cellIndex(i, (j - 1), (k - 1));

        eflux[3][IJK] = F1B4 * (eflux[6][IJK] + eflux[6][IJKM] + eflux[6][IJMK] + eflux[6][IJMKM]);
      }
    }
  }


  for(MInt k = m_noGhostLayers; k < m_nCells[0] - m_noGhostLayers; k++) {
    for(MInt j = m_noGhostLayers - 1; j < m_nCells[1] - m_noGhostLayers; j++) {
      for(MInt i = m_noGhostLayers; i < m_nCells[2] - m_noGhostLayers; i++) {
        const MInt IJK = cellIndex(i, j, k);
        const MInt IMJK = cellIndex((i - 1), j, k);
        const MInt IJKM = cellIndex(i, j, (k - 1));
        const MInt IMJKM = cellIndex((i - 1), j, (k - 1));

        fflux[3][IJK] = F1B4 * (fflux[6][IJK] + fflux[6][IJKM] + fflux[6][IMJK] + fflux[6][IMJKM]);
      }
    }
  }

  for(MInt k = m_noGhostLayers - 1; k < m_nCells[0] - m_noGhostLayers; k++) {
    for(MInt j = m_noGhostLayers; j < m_nCells[1] - m_noGhostLayers; j++) {
      for(MInt i = m_noGhostLayers; i < m_nCells[2] - m_noGhostLayers; i++) {
        const MInt IJK = cellIndex(i, j, k);
        const MInt IMJK = cellIndex((i - 1), j, k);
        const MInt IJMK = cellIndex(i, (j - 1), k);
        const MInt IMJMK = cellIndex((i - 1), (j - 1), k);

        gflux[3][IJK] = F1B4 * (gflux[6][IJK] + gflux[6][IMJK] + gflux[6][IJMK] + gflux[6][IMJMK]);
      }
    }
  }

  // separate loop for adding the dissipation term for tur kin viscosity transport variable
  for(MInt k = m_noGhostLayers; k < m_nCells[0] - m_noGhostLayers; k++) {
    for(MInt j = m_noGhostLayers; j < m_nCells[1] - m_noGhostLayers; j++) {
      for(MInt i = m_noGhostLayers; i < m_nCells[2] - m_noGhostLayers; i++) {
        const MInt IJK = cellIndex(i, j, k);
        const MInt IMJK = cellIndex(i - 1, j, k);
        const MInt IJMK = cellIndex(i, j - 1, k);
        const MInt IJKM = cellIndex(i, j, k - 1);
        const MFloat dissipation_term =
            ((eflux[3][IJK] - eflux[3][IMJK]) + (fflux[3][IJK] - fflux[3][IJMK]) + (gflux[3][IJK] - gflux[3][IJKM]));

        m_cells->rightHandSide[CV->RANS_FIRST][IJK] += dissipation_term * rho[IJK];
      }
    }
  }
}


void FvStructuredSolver3DRans::computeTurbViscosity() { (this->*compTurbVisc)(); }

void FvStructuredSolver3DRans::computeTurbViscosity_SA() {
  // OTHER variables required to calculate the laminar viscous fluxes
  const MFloat transPos = m_solver->m_ransTransPos;
  MFloat* const RESTRICT p = &m_cells->pvariables[PV->P][0];
  MFloat* const RESTRICT rho = &m_cells->pvariables[PV->RHO][0];
  MFloat* const RESTRICT nuTilde = &m_cells->pvariables[PV->RANS_VAR[0]][0];
  MFloat* const RESTRICT T = &m_cells->temperature[0];

  for(MInt i = 0; i < m_noCells; i++) {
    if(m_cells->coordinates[0][i] <= transPos) {
      nuTilde[i] = 0.0;
    } else {
      nuTilde[i] = mMin(nuTilde[i], 3000.0);
    }
    T[i] = m_solver->m_gamma * p[i] / rho[i];
    // decode the kinematic turbulent viscosity from the turb dynamic visc arrays
    const MFloat nuLaminar = SUTHERLANDLAW(T[i]) / rho[i];
    const MFloat chi = nuTilde[i] / (nuLaminar);
    const MFloat fv1 = pow(chi, 3) / (pow(chi, 3) + RM_SA_DV::cv1to3);
    const MFloat nuTurb = fv1 * nuTilde[i];
    m_cells->fq[FQ->NU_T][i] = nuTurb;
    m_cells->fq[FQ->MU_T][i] = rho[i] * nuTurb;
  }
}

void FvStructuredSolver3DRans::computeTurbViscosity_FS() {
  // OTHER variables required to calculate the laminar viscous fluxes
  const MFloat transPos = m_solver->m_ransTransPos;
  MFloat* const RESTRICT p = &m_cells->pvariables[PV->P][0];
  MFloat* const RESTRICT rho = &m_cells->pvariables[PV->RHO][0];
  MFloat* const RESTRICT nuTilde = &m_cells->pvariables[PV->RANS_VAR[0]][0];
  MFloat* const RESTRICT T = &m_cells->temperature[0];

  for(MInt i = 0; i < m_noCells; i++) {
    if(m_cells->coordinates[0][i] <= transPos) {
      nuTilde[i] = 0.0;
    } else {
      nuTilde[i] = mMin(nuTilde[i], 1200.0);
    }
    T[i] = m_solver->m_gamma * p[i] / rho[i];
    const MFloat nuLaminar = SUTHERLANDLAW(T[i]) / rho[i];
    const MFloat chi = nuTilde[i] / nuLaminar;
    const MFloat nuTurb = nuTilde[i] * pow(chi, 3.0) / (pow(chi, 3.0) + RM_FS::facv1to3);
    m_cells->fq[FQ->NU_T][i] = nuTurb;
    m_cells->fq[FQ->MU_T][i] = nuTurb * rho[i];
  }
}

template <MInt noVars>
void FvStructuredSolver3DRans::Muscl_AusmDV() {
  TRACE();

  // stencil identifier
  const MInt IJK[3] = {1, m_nCells[2], m_nCells[1] * m_nCells[2]};

  const MFloat* const RESTRICT x = ALIGNED_F(m_cells->coordinates[0]);
  const MFloat* const RESTRICT y = ALIGNED_F(m_cells->coordinates[1]);
  const MFloat* const RESTRICT z = ALIGNED_F(m_cells->coordinates[2]);
  const MFloat* const* const RESTRICT vars = ALIGNED_F(m_cells->pvariables);
  MFloat* const* const RESTRICT dss = ALIGNED_F(m_cells->dss);
  MFloat* const* const RESTRICT flux = ALIGNED_F(m_cells->flux);
  MFloat* const RESTRICT cellRhs = ALIGNED_MF(m_cells->rightHandSide[0]);

  const MUint noCells = m_noCells;
  const MFloat gamma = m_solver->m_gamma;
  const MFloat gammaMinusOne = gamma - 1.0;
  const MFloat FgammaMinusOne = F1 / gammaMinusOne;

  const MInt noCellsI = m_nCells[2] - 2;
  const MInt noCellsJ = m_nCells[1] - 2;
  const MInt noCellsK = m_nCells[0] - 2;

  const MInt noCellsIP1 = m_nCells[2] - 1;
  const MInt noCellsJP1 = m_nCells[1] - 1;
  const MInt noCellsKP1 = m_nCells[0] - 1;


  if(!m_dsIsComputed) {
    for(MInt dim = 0; dim < nDim; dim++) {
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for(MInt k = 0; k < noCellsKP1; k++) {
        for(MInt j = 0; j < noCellsJP1; j++) {
          for(MInt i = 0; i < noCellsIP1; i++) {
            const MInt I = cellIndex(i, j, k);
            const MInt IP1 = I + IJK[dim];
            dss[dim][I] = sqrt(POW2(x[IP1] - x[I]) + POW2(y[IP1] - y[I]) + POW2(z[IP1] - z[I]));
          }
        }
      }
    }

    m_dsIsComputed = true;
  }

  for(MInt dim = 0; dim < nDim; dim++) {
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(MInt k = 1; k < noCellsK; k++) {
      for(MInt j = 1; j < noCellsJ; j++) {
#if defined(MAIA_INTEL_COMPILER)
#pragma ivdep
#pragma vector always
#endif
        for(MInt i = 1; i < noCellsI; i++) {
          const MInt I = cellIndex(i, j, k);
          const MInt IP1 = I + IJK[dim];
          const MInt IM1 = I - IJK[dim];
          const MInt IP2 = I + 2 * IJK[dim];

          const MFloat DS = dss[dim][I];
          const MFloat DSM1 = dss[dim][IM1];
          const MFloat DSP1 = dss[dim][IP1];

          const MFloat DSP = DS / POW2(DSP1 + DS);
          const MFloat DSM = DS / POW2(DSM1 + DS);

          // unrolled the loop so the compiler
          // can optimize better
          const MFloat DQU = vars[PV->U][IP1] - vars[PV->U][I];
          const MFloat DQPU = vars[PV->U][IP2] - vars[PV->U][IP1];
          const MFloat DQMU = vars[PV->U][I] - vars[PV->U][IM1];
          MFloat UL = vars[PV->U][I] + DSM * (DSM1 * DQU + DS * DQMU);
          MFloat UR = vars[PV->U][IP1] - DSP * (DS * DQPU + DSP1 * DQU);

          const MFloat DQV = vars[PV->V][IP1] - vars[PV->V][I];
          const MFloat DQPV = vars[PV->V][IP2] - vars[PV->V][IP1];
          const MFloat DQMV = vars[PV->V][I] - vars[PV->V][IM1];
          MFloat VL = vars[PV->V][I] + DSM * (DSM1 * DQV + DS * DQMV);
          MFloat VR = vars[PV->V][IP1] - DSP * (DS * DQPV + DSP1 * DQV);

          const MFloat DQW = vars[PV->W][IP1] - vars[PV->W][I];
          const MFloat DQPW = vars[PV->W][IP2] - vars[PV->W][IP1];
          const MFloat DQMW = vars[PV->W][I] - vars[PV->W][IM1];
          MFloat WL = vars[PV->W][I] + DSM * (DSM1 * DQW + DS * DQMW);
          MFloat WR = vars[PV->W][IP1] - DSP * (DS * DQPW + DSP1 * DQW);

          const MFloat DQP = vars[PV->P][IP1] - vars[PV->P][I];
          const MFloat DQPP = vars[PV->P][IP2] - vars[PV->P][IP1];
          const MFloat DQMP = vars[PV->P][I] - vars[PV->P][IM1];
          const MFloat PL = vars[PV->P][I] + DSM * (DSM1 * DQP + DS * DQMP);
          const MFloat PR = vars[PV->P][IP1] - DSP * (DS * DQPP + DSP1 * DQP);

          const MFloat DQRHO = vars[PV->RHO][IP1] - vars[PV->RHO][I];
          const MFloat DQPRHO = vars[PV->RHO][IP2] - vars[PV->RHO][IP1];
          const MFloat DQMRHO = vars[PV->RHO][I] - vars[PV->RHO][IM1];
          const MFloat RHOL = vars[PV->RHO][I] + DSM * (DSM1 * DQRHO + DS * DQMRHO);
          const MFloat RHOR = vars[PV->RHO][IP1] - DSP * (DS * DQPRHO + DSP1 * DQRHO);

          const MFloat DQNUTILDE = vars[PV->RANS_VAR[0]][IP1] - vars[PV->RANS_VAR[0]][I];
          const MFloat DQPNUTILDE = vars[PV->RANS_VAR[0]][IP2] - vars[PV->RANS_VAR[0]][IP1];
          const MFloat DQMNUTILDE = vars[PV->RANS_VAR[0]][I] - vars[PV->RANS_VAR[0]][IM1];
          const MFloat NUTILDEL = vars[PV->RANS_VAR[0]][I] + DSM * (DSM1 * DQNUTILDE + DS * DQMNUTILDE);
          const MFloat NUTILDER = vars[PV->RANS_VAR[0]][IP1] - DSP * (DS * DQPNUTILDE + DSP1 * DQNUTILDE);
          const MFloat surf0 = m_cells->surfaceMetrics[dim * 3 + 0][I];
          const MFloat surf1 = m_cells->surfaceMetrics[dim * 3 + 1][I];
          const MFloat surf2 = m_cells->surfaceMetrics[dim * 3 + 2][I];
          const MFloat dxdtau = m_cells->dxt[dim][I];

          const MFloat FRHOL = F1 / RHOL;
          const MFloat FRHOR = F1 / RHOR;

          const MFloat PLfRHOL = PL / RHOL;
          const MFloat PRfRHOR = PR / RHOR;
          const MFloat e0 = PLfRHOL * FgammaMinusOne + 0.5 * (POW2(UL) + POW2(VL) + POW2(WL)) + PLfRHOL;
          const MFloat e1 = PRfRHOR * FgammaMinusOne + 0.5 * (POW2(UR) + POW2(VR) + POW2(WR)) + PRfRHOR;


          // compute lenght of metric vector for normalization
          const MFloat DGRAD = sqrt(POW2(surf0) + POW2(surf1) + POW2(surf2));
          const MFloat FDGRAD = F1 / DGRAD;

          // scale by metric length to get velocity in the new basis (get normalized basis vectors)
          const MFloat UUL = ((UL * surf0 + VL * surf1 + WL * surf2) - dxdtau) * FDGRAD;


          const MFloat UUR = ((UR * surf0 + VR * surf1 + WR * surf2) - dxdtau) * FDGRAD;

          MFloat AL = FRHOL * PL;
          MFloat AR = FRHOR * PR;

          const MFloat FALR = 2.0 / (AL + AR);
          const MFloat ALPHAL = AL * FALR;
          const MFloat ALPHAR = AR * FALR;

          AL = sqrt(gamma * AL);
          AR = sqrt(gamma * AR);
          AL = mMax(AL, AR);
          AR = AL;

          const MFloat XMAL = UUL / AL;
          const MFloat XMAR = UUR / AR;

          AL = AL * DGRAD;
          AR = AR * DGRAD;

          const MFloat RHOAL = AL * RHOL;
          const MFloat RHOAR = AR * RHOR;

          const MFloat FDV = 0.3;
          const MFloat DXDXEZ = m_cells->coordinates[0][IP1] - m_cells->coordinates[0][I];
          const MFloat DYDXEZ = m_cells->coordinates[1][IP1] - m_cells->coordinates[1][I];
          const MFloat DZDXEZ = m_cells->coordinates[2][IP1] - m_cells->coordinates[2][I];
          MFloat SV = 2.0 * DGRAD / (m_cells->cellJac[I] + m_cells->cellJac[IP1]) * (FDV + (F1 - FDV) * getPSI(I, dim));
          const MFloat SV1 = F1 * SV * DXDXEZ;
          const MFloat SV2 = F1 * SV * DYDXEZ;
          const MFloat SV3 = F1 * SV * DZDXEZ;

          const MFloat XMAL1 = mMin(F1, mMax(-F1, XMAL));
          const MFloat XMAR1 = mMin(F1, mMax(-F1, XMAR));

          MFloat FXMA = F1B2 * (XMAL1 + fabs(XMAL1));
          const MFloat XMALP = ALPHAL * (F1B4 * POW2(XMAL1 + F1) - FXMA) + FXMA + (mMax(F1, XMAL) - F1);
          FXMA = F1B2 * (XMAR1 - fabs(XMAR1));
          const MFloat XMARM = ALPHAR * (-F1B4 * POW2(XMAR1 - F1) - FXMA) + FXMA + (mMin(-F1, XMAR) + F1);

          const MFloat FLP = PL * ((F2 - XMAL1) * POW2(F1 + XMAL1));
          const MFloat FRP = PR * ((F2 + XMAR1) * POW2(F1 - XMAR1));
          const MFloat PLR = F1B4 * (FLP + FRP);

          const MFloat RHOUL = XMALP * RHOAL;
          const MFloat RHOUR = XMARM * RHOAR;
          const MFloat RHOU = RHOUL + RHOUR;
          const MFloat RHOU2 = F1B2 * RHOU;
          const MFloat ARHOU2 = fabs(RHOU2);

          const MFloat UUL2 = SV1 * UUL;
          const MFloat UUR2 = SV1 * UUR;
          UL = UL - UUL2;
          UR = UR - UUR2;
          const MFloat UUL3 = SV2 * UUL;
          const MFloat UUR3 = SV2 * UUR;
          VL = VL - UUL3;
          VR = VR - UUR3;
          const MFloat UUL4 = SV3 * UUL;
          const MFloat UUR4 = SV3 * UUR;
          WL = WL - UUL4;
          WR = WR - UUR4;

          flux[CV->RHO_U][I] = RHOU2 * (UL + UR) + ARHOU2 * (UL - UR) + PLR * surf0 + RHOUL * UUL2 + RHOUR * UUR2;
          flux[CV->RHO_V][I] = RHOU2 * (VL + VR) + ARHOU2 * (VL - VR) + PLR * surf1 + RHOUL * UUL3 + RHOUR * UUR3;
          flux[CV->RHO_W][I] = RHOU2 * (WL + WR) + ARHOU2 * (WL - WR) + PLR * surf2 + RHOUL * UUL4 + RHOUR * UUR4;
          flux[CV->RHO_E][I] = RHOU2 * (e0 + e1) + ARHOU2 * (e0 - e1) + PLR * dxdtau;
          flux[CV->RHO][I] = RHOU;
          flux[CV->RANS_VAR[0]][I] = RHOU2 * (NUTILDEL + NUTILDER) + ARHOU2 * (NUTILDEL - NUTILDER);
        }
      }
    }

    // FLUX BALANCE
    for(MUint v = 0; v < noVars; v++) {
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for(MInt k = m_noGhostLayers; k < m_nCells[0] - m_noGhostLayers; k++) {
        for(MInt j = m_noGhostLayers; j < m_nCells[1] - m_noGhostLayers; j++) {
#if defined(MAIA_INTEL_COMPILER)
#pragma ivdep
#pragma vector always
#endif
          for(MInt i = m_noGhostLayers; i < m_nCells[2] - m_noGhostLayers; i++) {
            const MInt I = cellIndex(i, j, k);
            const MInt IM1 = I - IJK[dim];
            const MUint offset = v * noCells;
            MFloat* const RESTRICT rhs = ALIGNED_F(cellRhs + offset);
            rhs[I] += flux[v][IM1] - flux[v][I];
          }
        }
      }
    }
  }
}

template void FvStructuredSolver3DRans::Muscl_AusmDV<5>();
template void FvStructuredSolver3DRans::Muscl_AusmDV<6>();
template void FvStructuredSolver3DRans::Muscl_AusmDV<7>();


template <MInt noVars>
void FvStructuredSolver3DRans::Muscl_AusmDV_Limited() {
  TRACE();

  // stencil identifier
  const MInt IJK[3] = {1, m_nCells[2], m_nCells[1] * m_nCells[2]};

  const MFloat* const RESTRICT x = ALIGNED_F(m_cells->coordinates[0]);
  const MFloat* const RESTRICT y = ALIGNED_F(m_cells->coordinates[1]);
  const MFloat* const RESTRICT z = ALIGNED_F(m_cells->coordinates[2]);
  const MFloat* const* const RESTRICT vars = ALIGNED_F(m_cells->pvariables);
  MFloat* const* const RESTRICT ql = ALIGNED_F(m_cells->ql);
  MFloat* const* const RESTRICT qr = ALIGNED_F(m_cells->qr);
  MFloat* const* const RESTRICT dss = ALIGNED_F(m_cells->dss);

  MFloat* const* const RESTRICT flux = ALIGNED_F(m_cells->flux);
  MFloat* const RESTRICT cellRhs = ALIGNED_MF(m_cells->rightHandSide[0]);

  const MFloat EPSLIM = 1e-14;
  const MFloat EPSS = 1e-34;
  const MUint noCells = m_noCells;
  const MFloat gamma = m_solver->m_gamma;
  const MFloat gammaMinusOne = gamma - 1.0;
  const MFloat FgammaMinusOne = F1 / gammaMinusOne;

  const MInt noCellsI = m_nCells[2] - 2;
  const MInt noCellsJ = m_nCells[1] - 2;
  const MInt noCellsK = m_nCells[0] - 2;

  const MInt noCellsIP1 = m_nCells[2] - 1;
  const MInt noCellsJP1 = m_nCells[1] - 1;
  const MInt noCellsKP1 = m_nCells[0] - 1;

  if(!m_dsIsComputed) {
    for(MInt dim = 0; dim < nDim; dim++) {
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for(MInt k = 0; k < noCellsKP1; k++) {
        for(MInt j = 0; j < noCellsJP1; j++) {
          for(MInt i = 0; i < noCellsIP1; i++) {
            const MInt I = cellIndex(i, j, k);
            const MInt IP1 = I + IJK[dim];
            dss[dim][I] = sqrt(POW2(x[IP1] - x[I]) + POW2(y[IP1] - y[I]) + POW2(z[IP1] - z[I]));
          }
        }
      }
    }

    m_dsIsComputed = true;
  }

  for(MInt dim = 0; dim < nDim; dim++) {
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(MInt k = 1; k < noCellsKP1; k++) {
      for(MInt j = 1; j < noCellsJP1; j++) {
        for(MInt i = 1; i < noCellsIP1; i++) {
          const MInt I = cellIndex(i, j, k);
          const MInt IP1 = I + IJK[dim];
          const MInt IM1 = I - IJK[dim];

          const MFloat DS = dss[dim][I];
          const MFloat DSL = dss[dim][IM1];
          const MFloat DSR = dss[dim][I];
          const MFloat FDS = DS / POW2(DSL + DSR) * 2.0;

          // unrolled the loop so the compiler
          // can optimize better
          const MFloat DQLU = DSL * (vars[PV->U][IP1] - vars[PV->U][I]);
          const MFloat DQRU = DSR * (vars[PV->U][I] - vars[PV->U][IM1]);
          const MFloat PHIU = mMax(F0, DQLU * DQRU / (POW2(DQLU) + POW2(DQRU) + EPSLIM * mMax(EPSS, DSL * DSR)));
          ql[PV->U][I] = vars[PV->U][I] + FDS * (DQLU + DQRU) * PHIU;
          qr[PV->U][IM1] = vars[PV->U][I] - FDS * (DQLU + DQRU) * PHIU;

          const MFloat DQLV = DSL * (vars[PV->V][IP1] - vars[PV->V][I]);
          const MFloat DQRV = DSR * (vars[PV->V][I] - vars[PV->V][IM1]);
          const MFloat PHIV = mMax(F0, DQLV * DQRV / (POW2(DQLV) + POW2(DQRV) + EPSLIM * mMax(EPSS, DSL * DSR)));
          ql[PV->V][I] = vars[PV->V][I] + FDS * (DQLV + DQRV) * PHIV;
          qr[PV->V][IM1] = vars[PV->V][I] - FDS * (DQLV + DQRV) * PHIV;

          const MFloat DQLW = DSL * (vars[PV->W][IP1] - vars[PV->W][I]);
          const MFloat DQRW = DSR * (vars[PV->W][I] - vars[PV->W][IM1]);
          const MFloat PHIW = mMax(F0, DQLW * DQRW / (POW2(DQLW) + POW2(DQRW) + EPSLIM * mMax(EPSS, DSL * DSR)));
          ql[PV->W][I] = vars[PV->W][I] + FDS * (DQLW + DQRW) * PHIW;
          qr[PV->W][IM1] = vars[PV->W][I] - FDS * (DQLW + DQRW) * PHIW;

          const MFloat DQLP = DSL * (vars[PV->P][IP1] - vars[PV->P][I]);
          const MFloat DQRP = DSR * (vars[PV->P][I] - vars[PV->P][IM1]);
          const MFloat PHIP = mMax(F0, DQLP * DQRP / (POW2(DQLP) + POW2(DQRP) + EPSLIM * mMax(EPSS, DSL * DSR)));
          ql[PV->P][I] = vars[PV->P][I] + FDS * (DQLP + DQRP) * PHIP;
          qr[PV->P][IM1] = vars[PV->P][I] - FDS * (DQLP + DQRP) * PHIP;

          const MFloat DQLRHO = DSL * (vars[PV->RHO][IP1] - vars[PV->RHO][I]);
          const MFloat DQRRHO = DSR * (vars[PV->RHO][I] - vars[PV->RHO][IM1]);
          const MFloat PHIRHO =
              mMax(F0, DQLRHO * DQRRHO / (POW2(DQLRHO) + POW2(DQRRHO) + EPSLIM * mMax(EPSS, DSL * DSR)));
          ql[PV->RHO][I] = vars[PV->RHO][I] + FDS * (DQLRHO + DQRRHO) * PHIRHO;
          qr[PV->RHO][IM1] = vars[PV->RHO][I] - FDS * (DQLRHO + DQRRHO) * PHIRHO;

          const MFloat DQLNUTILDE = DSL * (vars[PV->RANS_VAR[0]][IP1] - vars[PV->RANS_VAR[0]][I]);
          const MFloat DQRNUTILDE = DSR * (vars[PV->RANS_VAR[0]][I] - vars[PV->RANS_VAR[0]][IM1]);
          const MFloat PHINUTILDE = mMax(
              F0, DQLNUTILDE * DQRNUTILDE / (POW2(DQLNUTILDE) + POW2(DQRNUTILDE) + EPSLIM * mMax(EPSS, DSL * DSR)));
          ql[PV->RANS_VAR[0]][I] = vars[PV->RANS_VAR[0]][I] + FDS * (DQLNUTILDE + DQRNUTILDE) * PHINUTILDE;
          qr[PV->RANS_VAR[0]][IM1] = vars[PV->RANS_VAR[0]][I] - FDS * (DQLNUTILDE + DQRNUTILDE) * PHINUTILDE;
        }
      }
    }


    for(MInt k = 1; k < noCellsK; k++) {
      for(MInt j = 1; j < noCellsJ; j++) {
        for(MInt i = 1; i < noCellsI; i++) {
          const MInt I = cellIndex(i, j, k);
          const MInt IP1 = I + IJK[dim];

          MFloat UL = ql[PV->U][I];
          MFloat UR = qr[PV->U][I];

          MFloat VL = ql[PV->V][I];
          MFloat VR = qr[PV->V][I];

          MFloat WL = ql[PV->W][I];
          MFloat WR = qr[PV->W][I];

          const MFloat PL = ql[PV->P][I];
          const MFloat PR = qr[PV->P][I];

          const MFloat RHOL = ql[PV->RHO][I];
          const MFloat RHOR = qr[PV->RHO][I];

          const MFloat NUTILDEL = ql[PV->RANS_VAR[0]][I];
          const MFloat NUTILDER = qr[PV->RANS_VAR[0]][I];

          const MFloat surf0 = m_cells->surfaceMetrics[dim * 3 + 0][I];
          const MFloat surf1 = m_cells->surfaceMetrics[dim * 3 + 1][I];
          const MFloat surf2 = m_cells->surfaceMetrics[dim * 3 + 2][I];
          const MFloat dxdtau = m_cells->dxt[dim][I];

          const MFloat FRHOL = F1 / RHOL;
          const MFloat FRHOR = F1 / RHOR;

          const MFloat PLfRHOL = PL / RHOL;
          const MFloat PRfRHOR = PR / RHOR;
          const MFloat e0 = PLfRHOL * FgammaMinusOne + 0.5 * (POW2(UL) + POW2(VL) + POW2(WL)) + PLfRHOL;
          const MFloat e1 = PRfRHOR * FgammaMinusOne + 0.5 * (POW2(UR) + POW2(VR) + POW2(WR)) + PRfRHOR;


          // compute lenght of metric vector for normalization
          const MFloat DGRAD = sqrt(POW2(surf0) + POW2(surf1) + POW2(surf2));
          const MFloat FDGRAD = F1 / DGRAD;

          // scale by metric length to get velocity in the new basis (get normalized basis vectors)
          const MFloat UUL = ((UL * surf0 + VL * surf1 + WL * surf2) - dxdtau) * FDGRAD;


          const MFloat UUR = ((UR * surf0 + VR * surf1 + WR * surf2) - dxdtau) * FDGRAD;

          MFloat AL = FRHOL * PL;
          MFloat AR = FRHOR * PR;

          const MFloat FALR = 2.0 / (AL + AR);
          const MFloat ALPHAL = AL * FALR;
          const MFloat ALPHAR = AR * FALR;

          AL = sqrt(gamma * AL);
          AR = sqrt(gamma * AR);
          AL = mMax(AL, AR);
          AR = AL;

          const MFloat XMAL = UUL / AL;
          const MFloat XMAR = UUR / AR;

          AL = AL * DGRAD;
          AR = AR * DGRAD;

          const MFloat RHOAL = AL * RHOL;
          const MFloat RHOAR = AR * RHOR;

          const MFloat FDV = 0.3;
          const MFloat DXDXEZ = m_cells->coordinates[0][IP1] - m_cells->coordinates[0][I];
          const MFloat DYDXEZ = m_cells->coordinates[1][IP1] - m_cells->coordinates[1][I];
          const MFloat DZDXEZ = m_cells->coordinates[2][IP1] - m_cells->coordinates[2][I];
          MFloat SV = 2.0 * DGRAD / (m_cells->cellJac[I] + m_cells->cellJac[IP1]) * (FDV + (F1 - FDV) * getPSI(I, dim));
          const MFloat SV1 = F1 * SV * DXDXEZ;
          const MFloat SV2 = F1 * SV * DYDXEZ;
          const MFloat SV3 = F1 * SV * DZDXEZ;

          const MFloat XMAL1 = mMin(F1, mMax(-F1, XMAL));
          const MFloat XMAR1 = mMin(F1, mMax(-F1, XMAR));

          MFloat FXMA = F1B2 * (XMAL1 + fabs(XMAL1));
          const MFloat XMALP = ALPHAL * (F1B4 * POW2(XMAL1 + F1) - FXMA) + FXMA + (mMax(F1, XMAL) - F1);
          FXMA = F1B2 * (XMAR1 - fabs(XMAR1));
          const MFloat XMARM = ALPHAR * (-F1B4 * POW2(XMAR1 - F1) - FXMA) + FXMA + (mMin(-F1, XMAR) + F1);

          const MFloat FLP = PL * ((F2 - XMAL1) * POW2(F1 + XMAL1));
          const MFloat FRP = PR * ((F2 + XMAR1) * POW2(F1 - XMAR1));
          const MFloat PLR = F1B4 * (FLP + FRP);

          const MFloat RHOUL = XMALP * RHOAL;
          const MFloat RHOUR = XMARM * RHOAR;
          const MFloat RHOU = RHOUL + RHOUR;
          const MFloat RHOU2 = F1B2 * RHOU;
          const MFloat ARHOU2 = fabs(RHOU2);

          const MFloat UUL2 = SV1 * UUL;
          const MFloat UUR2 = SV1 * UUR;
          UL = UL - UUL2;
          UR = UR - UUR2;
          const MFloat UUL3 = SV2 * UUL;
          const MFloat UUR3 = SV2 * UUR;
          VL = VL - UUL3;
          VR = VR - UUR3;
          const MFloat UUL4 = SV3 * UUL;
          const MFloat UUR4 = SV3 * UUR;
          WL = WL - UUL4;
          WR = WR - UUR4;

          flux[CV->RHO_U][I] = RHOU2 * (UL + UR) + ARHOU2 * (UL - UR) + PLR * surf0 + RHOUL * UUL2 + RHOUR * UUR2;
          flux[CV->RHO_V][I] = RHOU2 * (VL + VR) + ARHOU2 * (VL - VR) + PLR * surf1 + RHOUL * UUL3 + RHOUR * UUR3;
          flux[CV->RHO_W][I] = RHOU2 * (WL + WR) + ARHOU2 * (WL - WR) + PLR * surf2 + RHOUL * UUL4 + RHOUR * UUR4;
          flux[CV->RHO_E][I] = RHOU2 * (e0 + e1) + ARHOU2 * (e0 - e1) + PLR * dxdtau;
          flux[CV->RHO][I] = RHOU;
          flux[CV->RANS_VAR[0]][I] = RHOU2 * (NUTILDEL + NUTILDER) + ARHOU2 * (NUTILDEL - NUTILDER);
        }
      }
    }

    // FLUX BALANCE
    for(MUint v = 0; v < noVars; v++) {
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for(MInt k = m_noGhostLayers; k < m_nCells[0] - m_noGhostLayers; k++) {
        for(MInt j = m_noGhostLayers; j < m_nCells[1] - m_noGhostLayers; j++) {
#if defined(MAIA_INTEL_COMPILER)
#pragma ivdep
#pragma vector always
#endif
          for(MInt i = m_noGhostLayers; i < m_nCells[2] - m_noGhostLayers; i++) {
            const MInt I = cellIndex(i, j, k);
            const MInt IM1 = I - IJK[dim];
            const MUint offset = v * noCells;
            MFloat* const RESTRICT rhs = ALIGNED_F(cellRhs + offset);
            rhs[I] += flux[v][IM1] - flux[v][I];
          }
        }
      }
    }
  }
}

template void FvStructuredSolver3DRans::Muscl_AusmDV_Limited<5>();
template void FvStructuredSolver3DRans::Muscl_AusmDV_Limited<6>();
template void FvStructuredSolver3DRans::Muscl_AusmDV_Limited<7>();

inline MInt FvStructuredSolver3DRans::cellIndex(MInt i, MInt j, MInt k) {
  return i + (j + k * m_nCells[1]) * m_nCells[2];
}

inline MInt FvStructuredSolver3DRans::getCellIdfromCell(MInt origin, MInt incI, MInt incJ, MInt incK) {
  return origin + incI + incJ * m_nCells[2] + incK * m_nCells[2] * m_nCells[1];
}

inline MFloat FvStructuredSolver3DRans::getPSI(MInt I, MInt dim) {
  const MFloat FK = 18.0;
  const MInt IJK[3] = {1, m_nCells[2], m_nCells[1] * m_nCells[2]};
  const MInt IP1 = I + IJK[dim];
  const MInt IM1 = I - IJK[dim];
  const MInt IP2 = I + 2 * IJK[dim];

  const MFloat PIM2 = m_cells->pvariables[PV->P][IM1];
  const MFloat PIM1 = m_cells->pvariables[PV->P][I];
  const MFloat PIP2 = m_cells->pvariables[PV->P][IP2];
  const MFloat PIP1 = m_cells->pvariables[PV->P][IP1];

  const MFloat PSI =
      mMin(F1,
           FK
               * mMax(mMax(fabs((PIM2 - PIM1) / mMin(PIM2, PIM1)), fabs((PIM1 - PIP1) / mMin(PIM1, PIP1))),
                      fabs((PIP1 - PIP2) / mMin(PIP1, PIP2))));
  return PSI;
}
