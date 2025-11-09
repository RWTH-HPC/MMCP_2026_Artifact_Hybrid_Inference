// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "fvstructuredsolver2d.h"
#include "COMM/mpioverride.h"
#include "INCLUDE/maiaconstants.h"
#include "IO/parallelio_hdf5.h"
#include "globals.h"

using namespace std;


/**
 * \brief Constructor for the 2D structured solver
 */
FvStructuredSolver2D::FvStructuredSolver2D(MInt solverId, StructuredGrid<2>* grid_, MBool* propertiesGroups,
                                           const MPI_Comm comm)
  : FvStructuredSolver<2>(solverId, grid_, propertiesGroups, comm) {
  TRACE();
  const MLong oldAllocatedBytes = allocatedBytes();

  // count the no of necessary FQ fields and allocate
  initializeFQField();

  // compute the cell center coordinates from point coordinates
  computeCellCentreCoordinates();

  if(m_rans) {
    m_structuredBndryCnd = new StructuredBndryCnd2D<true>(this, m_grid);
  } else {
    m_structuredBndryCnd = new StructuredBndryCnd2D<false>(this, m_grid);
  }

  allocateSingularities();

  // assign coordinates to all ghost points
  addGhostPointCoordinateValues();

  // allocate memory for aux data maps (cf,cp)
  allocateAuxDataMaps();

  // if we are Rans we should allocate a new RANS solver
  if(m_rans == true) {
    m_ransSolver = new FvStructuredSolver2DRans(this);
  }

  computeCellCentreCoordinates();
  RECORD_TIMER_START(m_timers[Timers::ComputeMetrics]);
  m_grid->computeMetrics();
  RECORD_TIMER_STOP(m_timers[Timers::ComputeMetrics]);
  RECORD_TIMER_START(m_timers[Timers::ComputeJacobian]);
  m_grid->computeJacobian();
  RECORD_TIMER_STOP(m_timers[Timers::ComputeJacobian]);

  // TODO_SS labels:FV By now I am not sure if Code performs correctly for wrongly oriented meshes
  for(MInt j = m_noGhostLayers - 1; j < m_nCells[0] - 1; j++) {
    for(MInt i = m_noGhostLayers - 1; i < m_nCells[1] - 1; i++) {
      const MInt cellId = cellIndex(i, j);
      if(m_cells->cellJac[cellId] < 0.0) {
        mTerm(1, "Negative Jacobian found! Check first if code can cope this!");
      }
    }
  }

  initFluxMethod();

  m_convergence = false;

  // Assign handlers to the correct boundary conditions
  assignBndryCells();

  // Computation of modified wall distance in porous computation requires to set the porosity,
  // Da-number etc. first; on the other hand we need to wait for initializeFQField to be called
  if(m_porous) {
    initPorous();
    // exchange6002();
  }

  if(m_rans) {
    if(m_ransMethod == RANS_SA_DV || m_ransMethod == RANS_KEPSILON) {
      if(m_setLocalWallDistance) m_structuredBndryCnd->computeLocalExtendedDistancesAndSetComm();
      // m_structuredBndryCnd->computeLocalWallDistances();
      else
        m_structuredBndryCnd->computeWallDistances();
    }

    // utau is required in RANS models for the computation of y+; In case the wall distance is very
    // large, the utau computation is skipped; In such situations the default value of utau=0 would
    // yield a y+=0 even very far away from walls
    // TODO_SS labels:FV currently in UTAU we save the inverse of the turbulent length scale
    if(FQ->neededFQVariables[FQ->UTAU]) {
      for(MInt cellId = 0; cellId < m_noCells; ++cellId) {
        m_cells->fq[FQ->UTAU][cellId] = 1.0;
      }
    }
  }

  printAllocatedMemory(oldAllocatedBytes, "FvStructuredSolver2D", m_StructuredComm);
  RECORD_TIMER_STOP(m_timers[Timers::Constructor]);
}

void FvStructuredSolver2D::initFluxMethod() {
  // set the MUSCL-scheme to the right function
  if(m_rans == true) {
    reconstructSurfaceData = &FvStructuredSolver2D::MusclRANS;
    viscFluxMethod = &FvStructuredSolver2D::viscousFluxRANS;
  } else {
    if(m_viscCompact)
      viscFluxMethod = &FvStructuredSolver2D::viscousFluxLESCompact<>;
    else
      viscFluxMethod = &FvStructuredSolver2D::viscousFluxLES<>;
    if(m_limiter) {
      switch(string2enum(m_limiterMethod)) {
        case ALBADA: {
          m_log << "Using VAN ALBADA limiter!" << endl;
          reconstructSurfaceData = &FvStructuredSolver2D::MusclAlbada;
          break;
        }
        default: {
          stringstream errorMessage;
          errorMessage << "Limiter function " << m_limiterMethod << " not implemented!" << endl;
          mTerm(1, AT_, errorMessage.str());
        }
      }
    } else {
      if(m_musclScheme == "Standard") {
        m_log << "Using unlimited MUSCL! (standard Formulation)" << endl;
        if(m_ausmScheme == "Standard") {
          switch(CV->noVariables) {
            case 4: {
              reconstructSurfaceData = &FvStructuredSolver2D::Muscl_<&FvStructuredSolver2D::AusmLES, 4>;
              break;
            }
            case 5: {
              reconstructSurfaceData = &FvStructuredSolver2D::Muscl_<&FvStructuredSolver2D::AusmLES, 5>;
              break;
            }
            case 6: {
              reconstructSurfaceData = &FvStructuredSolver2D::Muscl_<&FvStructuredSolver2D::AusmLES, 6>;
              break;
            }
            default: {
              stringstream errorMessage;
              errorMessage << "Number of Variables " << CV->noVariables << " not implemented in template AUSM!" << endl;
              mTerm(1, AT_);
            }
          }
        } else if(m_ausmScheme == "PTHRC") {
          switch(CV->noVariables) {
            case 4: {
              reconstructSurfaceData = &FvStructuredSolver2D::Muscl_<&FvStructuredSolver2D::AusmLES_PTHRC, 4>;
              break;
            }
            case 5: {
              reconstructSurfaceData = &FvStructuredSolver2D::Muscl_<&FvStructuredSolver2D::AusmLES_PTHRC, 5>;
              break;
            }
            case 6: {
              reconstructSurfaceData = &FvStructuredSolver2D::Muscl_<&FvStructuredSolver2D::AusmLES_PTHRC, 6>;
              break;
            }
            default: {
              stringstream errorMessage;
              errorMessage << "Number of Variables " << CV->noVariables << " not implemented in template AUSM!" << endl;
              mTerm(1, AT_);
            }
          }
        } else if(m_ausmScheme == "AUSMDV") {
          switch(CV->noVariables) {
            case 4: {
              reconstructSurfaceData = &FvStructuredSolver2D::Muscl_<&FvStructuredSolver2D::AusmDV, 4>;
              break;
            }
            case 5: {
              reconstructSurfaceData = &FvStructuredSolver2D::Muscl_<&FvStructuredSolver2D::AusmDV, 5>;
              break;
            }
            case 6: {
              reconstructSurfaceData = &FvStructuredSolver2D::Muscl_<&FvStructuredSolver2D::AusmDV, 6>;
              break;
            }
            default: {
              stringstream errorMessage;
              errorMessage << "Number of Variables " << CV->noVariables << " not implemented in template AUSM!" << endl;
              mTerm(1, AT_);
            }
          }
        }
      } else if(m_musclScheme == "Stretched") {
        m_log << "Using unlimited MUSCL (streched Grids)";
        mAlloc(m_cells->cellLength, nDim, m_noCells, "m_cells->cellLength", -F1, AT_);
        computeCellLength();
        if(m_ausmScheme == "Standard") {
          switch(CV->noVariables) {
            case 4: {
              reconstructSurfaceData = &FvStructuredSolver2D::MusclStretched_<&FvStructuredSolver2D::AusmLES, 4>;
              break;
            }
            case 5: {
              reconstructSurfaceData = &FvStructuredSolver2D::MusclStretched_<&FvStructuredSolver2D::AusmLES, 5>;
              break;
            }
            case 6: {
              reconstructSurfaceData = &FvStructuredSolver2D::MusclStretched_<&FvStructuredSolver2D::AusmLES, 6>;
              break;
            }
            default: {
              stringstream errorMessage;
              errorMessage << "Number of Variables " << CV->noVariables << " not implemented in template AUSM!" << endl;
              mTerm(1, AT_);
            }
          }
        } else if(m_ausmScheme == "PTHRC") {
          switch(CV->noVariables) {
            case 4: {
              reconstructSurfaceData = &FvStructuredSolver2D::MusclStretched_<&FvStructuredSolver2D::AusmLES_PTHRC, 4>;
              break;
            }
            case 5: {
              reconstructSurfaceData = &FvStructuredSolver2D::MusclStretched_<&FvStructuredSolver2D::AusmLES_PTHRC, 5>;
              break;
            }
            case 6: {
              reconstructSurfaceData = &FvStructuredSolver2D::MusclStretched_<&FvStructuredSolver2D::AusmLES_PTHRC, 6>;
              break;
            }
            default: {
              stringstream errorMessage;
              errorMessage << "Number of Variables " << CV->noVariables << " not implemented in template AUSM!" << endl;
              mTerm(1, AT_);
            }
          }
        }
      }
    }
  }
}

FvStructuredSolver2D::~FvStructuredSolver2D() {
  TRACE();

  delete m_structuredBndryCnd;
  if(m_rans) {
    delete m_ransSolver;
  }

  if(m_hasSingularity) {
    delete[] m_singularity;
  }
}


void FvStructuredSolver2D::computeCellLength() {
  // this function can be moved into the MusclSchemeStreched later but for testing it is easier
  // REMEMBER: FOR MOVINg GRIDS THIS NEEDS TO BE CALLED EACH TIME

  for(MInt j = 0; j < m_nCells[0]; j++) {
    for(MInt i = 0; i < m_nCells[1]; i++) {
      const MInt cellId = cellIndex(i, j);
      const MInt P1 = getPointIdFromCell(i, j);
      const MInt P2 = getPointIdFromPoint(P1, 1, 0);
      const MInt P3 = getPointIdFromPoint(P1, 1, 1);
      const MInt P4 = getPointIdFromPoint(P1, 0, 1);
      //----------Idirection
      // face 1
      const MFloat f1x = F1B2 * (m_grid->m_coordinates[0][P1] + m_grid->m_coordinates[0][P4]);
      const MFloat f1y = F1B2 * (m_grid->m_coordinates[1][P1] + m_grid->m_coordinates[1][P4]);
      // face 2
      const MFloat f2x = F1B2 * (m_grid->m_coordinates[0][P2] + m_grid->m_coordinates[0][P3]);
      const MFloat f2y = F1B2 * (m_grid->m_coordinates[1][P2] + m_grid->m_coordinates[1][P3]);
      m_cells->cellLength[0][cellId] = sqrt(POW2(f2x - f1x) + POW2(f2y - f1y));
      //----------Jdirection
      // face 3
      const MFloat f3x = F1B2 * (m_grid->m_coordinates[0][P1] + m_grid->m_coordinates[0][P2]);
      const MFloat f3y = F1B2 * (m_grid->m_coordinates[1][P1] + m_grid->m_coordinates[1][P2]);
      // face 4
      const MFloat f4x = F1B4 * (m_grid->m_coordinates[0][P3] + m_grid->m_coordinates[0][P4]);
      const MFloat f4y = F1B4 * (m_grid->m_coordinates[1][P3] + m_grid->m_coordinates[1][P4]);
      m_cells->cellLength[1][cellId] = sqrt(POW2(f4x - f3x) + POW2(f4y - f3y));
    }
  }
}


/**
 * \brief initalize the solution step
 *
 * \author Pascal Meysonnat
 *
 */

void FvStructuredSolver2D::initSolutionStep(MInt mode) {
  TRACE();

  std::ignore = mode;

  // Compute infinity values from property file
  // and (if no restart) fill cells according
  // to the initialCondition property
  initialCondition();

  if(m_restart) {
    loadRestartFile();
  }

  setTimeStep();

  // initialize moving grid
  // functions and move grid
  // to correct position
  if(m_movingGrid) {
    RECORD_TIMER_START(m_timers[Timers::MovingGrid]);
    RECORD_TIMER_START(m_timers[Timers::MGMoveGrid]);
    initMovingGrid();
    RECORD_TIMER_STOP(m_timers[Timers::MGMoveGrid]);
    RECORD_TIMER_STOP(m_timers[Timers::MovingGrid]);
  }

  // Get the correct values
  // in the exchange ghostcells
  exchange();

  // Call the init function of each BC
  initBndryCnds();

  // Apply boundary conditions
  // and fill the non-exchange ghostcells
  applyBoundaryCondition();

  // Convert SA turb. quantities to k-epsilon when restarting from SA solution
  if(m_restart) convertSA2KEPS();

  // Check for NaNs
  checkNans();

  computeConservativeVariables();
}


void FvStructuredSolver2D::convertSA2KEPS() {
  TRACE();

  MBool restartFromSA = false;
  restartFromSA = Context::getSolverProperty<MBool>("restartFromSA", m_solverId, AT_, &restartFromSA);
  if(!restartFromSA) return;

  if(domainId() == 0) cout << "\033[1;31m !!!Converting SA turbulence quantities to k & epsilon!!!\033[0m\n" << endl;

  MFloat* const RESTRICT p = &m_cells->pvariables[PV->P][0];
  MFloat* const RESTRICT rho = &m_cells->pvariables[PV->RHO][0];
  MFloat* const RESTRICT nuTilde = &m_cells->pvariables[PV->RANS_VAR[0]][0];

  for(MInt i = 0; i < m_noCells; ++i) {
    const MFloat T = m_gamma * p[i] / rho[i];
    // decode the kinematic turbulent viscosity from the turb dynamic visc arrays
    const MFloat nuLaminar = SUTHERLANDLAW(T) / rho[i];
    const MFloat chi = nuTilde[i] / (nuLaminar);
    const MFloat fv1 = pow(chi, 3) / (pow(chi, 3) + RM_SA_DV::cv1to3);
    const MFloat nuTurb = fv1 * nuTilde[i];
    m_cells->fq[FQ->NU_T][i] = nuTurb;
    m_cells->fq[FQ->MU_T][i] = rho[i] * nuTurb;
  }

  // compute friction velocity at wall
  m_structuredBndryCnd->computeFrictionCoef();
  // communicate wall properties to all cells
  m_structuredBndryCnd->distributeWallAndFPProperties();

  const MFloat rRe0 = 1.0 / m_Re0;
  const MFloat fac_nonDim = m_keps_nonDimType ? 1.0 : PV->UInfinity;
  MFloat* const RESTRICT u = &m_cells->pvariables[PV->U][0];
  MFloat* const RESTRICT v = &m_cells->pvariables[PV->V][0];
  const MFloat* const RESTRICT utau = &m_cells->fq[FQ->UTAU][0];
  const MFloat* const RESTRICT wallDist = &m_cells->fq[FQ->WALLDISTANCE][0];
  const MFloat* const RESTRICT muTurb = &m_cells->fq[FQ->MU_T][0];
  for(MInt j = m_noGhostLayers; j < m_nCells[0] - m_noGhostLayers; j++) {
    for(MInt i = m_noGhostLayers; i < m_nCells[1] - m_noGhostLayers; i++) {
      const MInt IJ = cellIndex(i, j);
      const MInt IPJ = cellIndex((i + 1), j);
      const MInt IMJ = cellIndex(i - 1, j);
      const MInt IJM = cellIndex(i, j - 1);
      const MInt IJP = cellIndex(i, (j + 1));

      const MFloat invCellJac = 1.0 / m_cells->cellJac[IJ];

      const MFloat dudxi = 0.5 * (u[IPJ] - u[IMJ]);
      const MFloat dudet = 0.5 * (u[IJP] - u[IJM]);
      const MFloat dvdxi = 0.5 * (v[IPJ] - v[IMJ]);
      const MFloat dvdet = 0.5 * (v[IJP] - v[IJM]);

      const MFloat dudx =
          dudxi * m_cells->cellMetrics[xsd * 2 + xsd][IJ] + dudet * m_cells->cellMetrics[ysd * 2 + xsd][IJ];
      const MFloat dudy =
          dudxi * m_cells->cellMetrics[xsd * 2 + ysd][IJ] + dudet * m_cells->cellMetrics[ysd * 2 + ysd][IJ];
      const MFloat dvdx =
          dvdxi * m_cells->cellMetrics[xsd * 2 + xsd][IJ] + dvdet * m_cells->cellMetrics[ysd * 2 + xsd][IJ];
      const MFloat dvdy =
          dvdxi * m_cells->cellMetrics[xsd * 2 + ysd][IJ] + dvdet * m_cells->cellMetrics[ysd * 2 + ysd][IJ];

      const MFloat P_ = (rRe0 / POW2(fac_nonDim) * invCellJac * muTurb[IJ]
                         * (2.0 * POW2(dudx) + 2.0 * POW2(dvdy) + POW2(dudy + dvdx)))
                        * invCellJac;
      const MFloat rhoEps = P_ / fac_nonDim;
      m_cells->pvariables[PV->RANS_VAR[1]][IJ] = rhoEps / rho[IJ];
      m_cells->variables[PV->RANS_VAR[1]][IJ] = rhoEps;

      //
      //      const MFloat T = m_gamma*p[i]/rho[i];
      //      const MFloat nuLaminar = SUTHERLANDLAW(T)/rho[IJ];
      // TODO_SS labels:FV
      const MFloat yp = utau[IJ] * wallDist[IJ]; // utau[IJ]*wallDist[IJ]/nuLaminar;
      const MFloat f_mu = 1 - exp(-RM_KEPS::C3 * m_Re0 * yp);
      const MFloat rhoTKE = sqrt(rhoEps * muTurb[IJ] / (RM_KEPS::C_mu * f_mu * m_Re0 * fac_nonDim));
      m_cells->pvariables[PV->RANS_VAR[0]][IJ] = rhoTKE / rho[IJ];
      m_cells->variables[PV->RANS_VAR[0]][IJ] = rhoTKE;
    }
  }

  //
  exchange();
  applyBoundaryCondition();
}

/**
 * \brief Computation of infinity values for the conservative and primitive variables
 *        Initialization ot the entire flow field
 */

void FvStructuredSolver2D::initialCondition() {
  TRACE();

  const MFloat gammaMinusOne = m_gamma - 1.0;
  MFloat UT;
  // MFloat pressure=F0;
  // MFloat Frho;         switched off cause of compiler

  PV->TInfinity = 1.0 / (1.0 + F1B2 * gammaMinusOne * POW2(m_Ma));
  UT = m_Ma * sqrt(PV->TInfinity);
  PV->UInfinity = UT * cos(m_angle[0]) * cos(m_angle[1]);
  PV->VInfinity = UT * sin(m_angle[0]) * cos(m_angle[1]);
  PV->VVInfinity[0] = PV->UInfinity;
  PV->VVInfinity[1] = PV->VInfinity;
  PV->PInfinity = pow(PV->TInfinity, (m_gamma / gammaMinusOne)) / m_gamma;

  // compute conservative variables
  CV->rhoInfinity = pow(PV->TInfinity, (1.0 / gammaMinusOne));
  CV->rhoUInfinity = CV->rhoInfinity * PV->UInfinity;
  CV->rhoVInfinity = CV->rhoInfinity * PV->VInfinity;
  CV->rhoVVInfinity[0] = CV->rhoUInfinity;
  CV->rhoVVInfinity[1] = CV->rhoVInfinity;
  CV->rhoEInfinity = PV->PInfinity / gammaMinusOne + CV->rhoInfinity * (F1B2 * POW2(UT));

  // internal Reynolds number Re0 = Re / ( rho8*M*sqrt(T8)/T8^F072)
  m_Re0 = m_Re * SUTHERLANDLAW(PV->TInfinity) / (CV->rhoInfinity * m_Ma * sqrt(PV->TInfinity));

  // reference enthalpies (needed for combustion computations)
  m_hInfinity = PV->PInfinity / CV->rhoInfinity * m_gamma / gammaMinusOne;

  // reference time (convection time)
  m_timeRef = UT / m_referenceLength;

  m_deltaP = F0;
  // pressure loss per unit length dp = rho_00 u_tau^2 L / D ) here: D=1.0, L=1;
  // channel: dp = rho_00 u_tau^2 L / D )
  // m_deltaP = POW2( m_Ma * m_ReTau * sqrt(PV->TInfinity) / m_Re  ) * CV->rhoInfinity / m_referenceLength;x
  // result is obtained by making deltap dimensionless with a_0^2 and rho_0
  // pipe: dp = lambda * L/D * rho/2 * u^2, lambda = 0.3164 Re^(-1/4) (Blasius)

  if(m_rans) {
    if(m_ransMethod == RANS_KEPSILON) {
      const MFloat rho = CV->rhoInfinity;
      const MFloat lamVisc = SUTHERLANDLAW(PV->TInfinity);
      const MFloat k8 = m_keps_nonDimType ? 1.5 * POW2(UT * m_I) : 1.5 * POW2(m_I);
      PV->ransInfinity[0] = k8;
      if(m_kepsICMethod == 1) {
        PV->ransInfinity[1] = RM_KEPS::C_mu * pow(k8, 1.5) / m_epsScale;
      } else if(m_kepsICMethod == 2) {
        PV->ransInfinity[1] = m_keps_nonDimType ? RM_KEPS::C_mu * rho * POW2(k8) * m_Re0 / (lamVisc * m_epsScale)
                                                : RM_KEPS::C_mu * rho * POW2(k8) * m_Re0 * UT / (lamVisc * m_epsScale);
      } else {
        PV->ransInfinity[1] = m_keps_nonDimType ? m_epsScale * UT * UT * UT : m_epsScale;
      }
      CV->ransInfinity[0] = CV->rhoInfinity * PV->ransInfinity[0];
      CV->ransInfinity[1] = CV->rhoInfinity * PV->ransInfinity[1];
    } else {
      MFloat chi = 0.1;
      chi = Context::getSolverProperty<MFloat>("chi", m_solverId, AT_, &chi);
      const MFloat lamVisc = SUTHERLANDLAW(PV->TInfinity);
      CV->ransInfinity[0] = chi * (lamVisc);
      PV->ransInfinity[0] = chi * (lamVisc / CV->rhoInfinity);
    }
  }

  m_log << "**************************" << endl;
  m_log << "Initial Condition summary" << endl;
  m_log << "**************************" << endl;
  m_log << "Re = " << m_Re << endl;
  m_log << "Re0 = " << m_Re0 << endl;
  m_log << "Ma = " << m_Ma << endl;
  m_log << "TInfinity = " << PV->TInfinity << endl;
  m_log << "UInfinity = " << PV->UInfinity << endl;
  m_log << "VInfinity = " << PV->VInfinity << endl;
  m_log << "PInfinity = " << PV->PInfinity << endl;
  m_log << "rhoInfinity = " << CV->rhoInfinity << endl;
  m_log << "rhoEInfinity = " << CV->rhoEInfinity << endl;
  for(MInt ransVarId = 0; ransVarId < m_noRansEquations; ++ransVarId)
    m_log << "Rans" << ransVarId << "Infinity = " << PV->ransInfinity[ransVarId] * 1e8 << "e-8" << endl;
  m_log << "referenceTime = " << m_timeRef << endl;

  if(domainId() == 0) {
    cout << "**************************" << endl;
    cout << "Initial Condition summary" << endl;
    cout << "**************************" << endl;
    cout << "Re = " << m_Re << endl;
    cout << "Re0 = " << m_Re0 << endl;
    cout << "Ma = " << m_Ma << endl;
    cout << "TInfinity = " << PV->TInfinity << endl;
    cout << "UInfinity = " << PV->UInfinity << endl;
    cout << "VInfinity = " << PV->VInfinity << endl;
    cout << "PInfinity = " << PV->PInfinity << endl;
    cout << "rhoInfinity = " << CV->rhoInfinity << endl;
    cout << "rhoEInfinity = " << CV->rhoEInfinity << endl;
    for(MInt ransVarId = 0; ransVarId < m_noRansEquations; ++ransVarId)
      cout << "Rans" << ransVarId << "Infinity = " << PV->ransInfinity[ransVarId] * 1e8 << "e-8" << endl;
    cout << "referenceTime = " << m_timeRef << endl;
  }


  if(!m_restart) {
    // inflow condition
    // ----------------

    switch(m_initialCondition) {
      case 0:
      case 465: // quiscient state
      {
        MFloat u_infty[nDim]{};
        if(m_initialCondition == 0) {
          u_infty[0] = PV->VVInfinity[0];
          u_infty[1] = PV->VVInfinity[1];
        }
        // parallel inflow field
        for(MInt cellId = 0; cellId < m_noCells; cellId++) {
          // go through every cell
          m_cells->pvariables[PV->RHO][cellId] = CV->rhoInfinity;
          for(MInt i = 0; i < nDim; i++) {
            m_cells->pvariables[PV->VV[i]][cellId] = u_infty[i]; // PV->VVInfinity[i];
          }

          m_cells->pvariables[PV->P][cellId] = PV->PInfinity;

          if(m_rans) {
            for(MInt ransVarId = 0; ransVarId < m_noRansEquations; ++ransVarId) {
              m_cells->pvariables[PV->RANS_VAR[ransVarId]][cellId] = PV->ransInfinity[ransVarId];
            }
          }
        }
        break;
      }
      case 43: {
        // parallel inflow field
        for(MInt cellId = 0; cellId < m_noCells; cellId++) {
          // go through every cell
          m_cells->pvariables[PV->RHO][cellId] = CV->rhoInfinity;
          for(MInt i = 0; i < nDim; i++) {
            m_cells->pvariables[PV->VV[i]][cellId] = F0;
          }

          m_cells->pvariables[PV->P][cellId] = PV->PInfinity;

          MFloat radius =
              sqrt(POW2(m_cells->coordinates[0][cellId] - 0.5) + POW2(m_cells->coordinates[1][cellId] - 0.5));
          // impose pressure peak in the middle of the domain
          if(radius <= 0.05) {
            MFloat pAmp = 0.005;
            MFloat pressureSignal = sin(radius / 0.05 * PI) * pAmp + PV->PInfinity;
            m_cells->pvariables[PV->P][cellId] = pressureSignal;
          }
        }
        break;
      }
      case 314: // stagnating flow field
      {
        for(MInt cellid = 0; cellid < m_noCells; cellid++) {
          m_cells->pvariables[PV->RHO][cellid] = CV->rhoInfinity;
          for(MInt i = 0; i < nDim; i++) {
            m_cells->pvariables[PV->VV[i]][cellid] = F0;
          }

          m_cells->pvariables[PV->P][cellid] = PV->PInfinity;
        }
        cout << "I.C. stagnating flow field was applied! " << endl;
        break;
      }
      case 333: {
        // parallel inflow field
        for(MInt cellId = 0; cellId < m_noCells; cellId++) {
          // go through every cell
          m_cells->pvariables[PV->RHO][cellId] = CV->rhoInfinity;
          for(MInt i = 0; i < nDim; i++) {
            m_cells->pvariables[PV->VV[i]][cellId] = PV->VVInfinity[i];
          }

          m_cells->pvariables[PV->P][cellId] = PV->PInfinity;

          // impose pressure peak in the middle of the domain
          if(m_cells->coordinates[0][cellId] > 0.4 && m_cells->coordinates[0][cellId] < 0.5) {
            MFloat pAmp = 0.005;
            MFloat xCoordinate = m_cells->coordinates[0][cellId] - 0.4;
            MFloat pressureSignal = sin(xCoordinate / 0.1 * PI) * pAmp + PV->PInfinity;
            m_cells->pvariables[PV->P][cellId] = pressureSignal / (m_gamma - F1);
          }
        }
        break;
      }
      case 79091: // Turbulent plate
      {
        const MFloat epss = 1e-10;
        const MFloat reTheta = 1000.0;
        const MFloat theta = 1.0;
        const MFloat delta0 = 72.0 / 7.0 * theta;
        const MFloat K = 0.4;
        const MFloat C1 = 3.573244189003983; // With coles
        const MFloat PI1 = 0.55;
        const MFloat cf = 0.024 / pow(reTheta, 0.25);

        for(MInt j = 0; j < m_nCells[0]; j++) {
          for(MInt i = 0; i < m_nCells[1]; i++) {
            const MInt cellId = cellIndex(i, j);
            const MFloat mu = SUTHERLANDLAW(PV->TInfinity);
            const MFloat utau = sqrt(cf / 2.0) * m_Ma * sqrt(PV->TInfinity);
            const MFloat yplus = m_cells->coordinates[1][cellId] * sqrt(cf / 2.) * CV->rhoUInfinity / mu * m_Re0;
            const MFloat eta = m_cells->coordinates[1][cellId] / delta0; // y/delta

            // 1-7th profile
            // log-law + wake
            if(m_cells->coordinates[1][cellId] > delta0) {
              m_cells->pvariables[PV->U][cellId] = PV->UInfinity; // Outside BL
            } else if(yplus < 10) {
              m_cells->pvariables[PV->U][cellId] = utau * yplus;
            } else {
              m_cells->pvariables[PV->U][cellId] = mMin(
                  utau * ((1. / K) * log(max(yplus, epss)) + C1 + 2 * PI1 / K * (3 * eta * eta - 2 * eta * eta * eta)),
                  PV->UInfinity);
            }

            m_cells->pvariables[PV->V][cellId] = PV->VInfinity;
            m_cells->pvariables[PV->RHO][cellId] = CV->rhoInfinity;
            m_cells->pvariables[PV->P][cellId] = PV->PInfinity;

            if(m_rans) {
              for(MInt ransVarId = 0; ransVarId < m_noRansEquations; ++ransVarId) {
                m_cells->pvariables[PV->RANS_VAR[ransVarId]][cellId] = PV->ransInfinity[ransVarId];
              }
            }
          }
        }

        break;
      }
      case 999: {
        // blasius laminar bl
        m_log << "Blasius initial condition (incompressible)" << endl;
        if(!m_useBlasius) mTerm(1, "property Blasius not set. Refer to the description of the property");
        if(domainId() == 0) {
          ofstream blasiusf;
          // write velocity to file
          blasiusf.open("velocity_x0.dat", ios::trunc);
          if(blasiusf) {
            blasiusf << "#y eta u v" << endl;
            MFloat d0 = 0.0;
            MFloat d1 = 0.0;
            MFloat d2 = 0.0;
            MBool d0Set = false;

            for(MInt i = 0; i < m_blasius_noPoints; i++) {
              // coord
              const MFloat y = i * 0.05;
              blasiusf << y << " " << getBlasiusEta(F0, y);
              // velocity
              MFloat vel[nDim];
              getBlasiusVelocity(F0, y, vel);
              for(MInt dim = 0; dim < nDim; dim++)
                blasiusf << " " << vel[dim] / PV->UInfinity;
              blasiusf << endl;

              if(!d0Set && vel[0] >= 0.99 * PV->UInfinity) {
                d0 = y;
                d0Set = true;
              }

              if(y < 10.0) {
                d1 += (1 - vel[0] / PV->UInfinity) * 0.05;
                d2 += vel[0] / PV->UInfinity * (1 - vel[0] / PV->UInfinity) * 0.05;
              }
            }
            blasiusf.close();

            cout << "x0: " << m_blasius_x0 << endl;
            cout << "d0: " << d0 << " d1: " << d1 << " d2: " << d2 << endl;
          }
        }

        for(MInt cellid = 0; cellid < m_noCells; cellid++) {
          MFloat vel[nDim];
          getBlasiusVelocity(cellid, vel);
          for(MInt i = 0; i < nDim; i++) {
            m_cells->pvariables[PV->VV[i]][cellid] = vel[i];
          }
          m_cells->pvariables[PV->P][cellid] = PV->PInfinity;
          m_cells->pvariables[PV->RHO][cellid] = CV->rhoInfinity;
        }

        break;
      }

      default: {
        // put the parallel flow field input in here
        // force output that no specific initial condition was chosen
        m_log << "No (correct) initial Condition is given! Used initial Condtion of parallel inflow!!!!!" << endl;
        break;
      }
    }
  }
}

void FvStructuredSolver2D::initMovingGrid() {
  TRACE();

  MInt pointId = 0;

  // First approach: save whole mesh in m_mgInitCoordinates (for analytical channel with indentation)

  for(MInt j = 0; j < m_nPoints[0]; ++j) {
    for(MInt i = 0; i < m_nPoints[1]; ++i) {
      pointId = pointIndex(i, j);
      for(MInt isd = xsd; isd < nDim; ++isd) {
        m_grid->m_initCoordinates[isd][pointId] = m_grid->m_coordinates[isd][pointId];
      }
    }
  }

  // Second approach: save only parts of the mesh depending on moving grid case
  switch(m_gridMovingMethod) {
    case 4: {
      // do nothing
      break;
    }
    case 10:
      // travelling wave defined by non-plus units
      {
        m_travelingWave = true;
        if(!m_restart) {
          m_waveTimeStepComputed = false;
        }
        m_waveSpeed = 0.0;
        m_waveLength = 0.0;
        m_waveAmplitude = 0.0;
        m_waveCellsPerWaveLength = 1;
        if(!m_restart) {
          m_waveNoStepsPerCell = 1;
        }

        // time needs to be constant for traveling wave
        m_constantTimeStep = true;

        /*! \property
          \page propertiesFVSTRCTRD
          \section waveLength
          <code>MInt FvStructuredSolver::m_waveLength </code>\n
          default = <code> 1.0 </code>\n \n
          Wavelength of the traveling wave.\n
          Possible values are:\n
          <ul>
          <li>Float > 0.0</li>
          </ul>
          Keywords: <i>WAVE, MOVING, STRUCTURED</i>
        */
        m_waveLength = 0.0;
        m_waveLength = Context::getSolverProperty<MFloat>("waveLength", m_solverId, AT_, &m_waveLengthPlus);

        /*! \property
          \page propertiesFVSTRCTRD
          \section waveAmplitude
          <code>MInt FvStructuredSolver::m_waveAmplitude </code>\n
          default = <code> 1.0 </code>\n \n
          Amplitude of the traveling wave.\n
          Possible values are:\n
          <ul>
          <li>Float > 0.0</li>
          </ul>
          Keywords: <i>WAVE, MOVING, STRUCTURED</i>
        */
        m_waveAmplitude = 0.0;
        m_waveAmplitude = Context::getSolverProperty<MFloat>("waveAmplitude", m_solverId, AT_, &m_waveAmplitudePlus);

        /*! \property
          \page propertiesFVSTRCTRD
          \section waveTime
          <code>MInt FvStructuredSolver::m_waveTime </code>\n
          default = <code> 1.0 </code>\n \n
          Period time of the traveling wave.\n
          Possible values are:\n
          <ul>
          <li>Float > 0.0</li>
          </ul>
          Keywords: <i>WAVE, MOVING, STRUCTURED</i>
        */
        m_waveTime = 0.0;
        m_waveTime = Context::getSolverProperty<MFloat>("waveTime", m_solverId, AT_, &m_waveTimePlus);

        /*! \property
          \page propertiesFVSTRCTRD
          \section waveYBeginTransition
          <code>MInt FvStructuredSolver::m_waveYBeginTransition </code>\n
          default = <code> 1.0 </code>\n \n
          End of the transition from wave to flat in x-dir.\n
          Possible values are:\n
          <ul>
          <li>Float</li>
          </ul>
          Keywords: <i>WAVE, MOVING, STRUCTURED</i>
        */
        m_waveYBeginTransition = 2000000.0;
        if(Context::propertyExists("waveYBeginTransition", m_solverId)) {
          m_waveYBeginTransition =
              Context::getSolverProperty<MFloat>("waveYBeginTransition", m_solverId, AT_, &m_waveYBeginTransition);
        }

        /*! \property
          \page propertiesFVSTRCTRD
          \section waveYEndTransition
          <code>MInt FvStructuredSolver::m_waveYEndTransition </code>\n
          default = <code> 1.0 </code>\n \n
          End of the transition from wave to flat in x-dir.\n
          Possible values are:\n
          <ul>
          <li>Float</li>
          </ul>
          Keywords: <i>WAVE, MOVING, STRUCTURED</i>
        */
        m_waveYEndTransition = 2000000.0;
        if(Context::propertyExists("waveYEndTransition", m_solverId)) {
          m_waveYEndTransition =
              Context::getSolverProperty<MFloat>("waveYEndTransition", m_solverId, AT_, &m_waveYEndTransition);
        }


        /*! \property
          \page propertiesFVSTRCTRD
          \section waveTemporalTransition
          <code>MInt FvStructuredSolver::m_waveOutEndTransition </code>\n
          default = <code> 500.0 </code>\n \n
          Acoustic time for wave actuation to transiate from flat plate\n
          to fully extended.\n
          Possible values are:\n
          <ul>
          <li>Float</li>
          </ul>
          Keywords: <i>WAVE, MOVING, STRUCTURED</i>
        */
        m_waveTemporalTransition = 500.0;
        if(Context::propertyExists("waveTemporalTransition", m_solverId)) {
          m_waveTemporalTransition =
              Context::getSolverProperty<MFloat>("waveTemporalTransition", m_solverId, AT_, &m_waveOutEndTransition);
        }

        m_waveSpeed = m_waveLength / m_waveTime;
        const MFloat deltaX = abs(m_grid->m_coordinates[0][0] - m_grid->m_coordinates[0][m_nPoints[1]]);
        m_waveCellsPerWaveLength = round(m_waveLength / deltaX);

        const MFloat speedAmplitude = 2 * PI * m_waveAmplitude / m_waveTime;

        m_log << "/////////////////// TRAVELING WAVE /////////////////////////////" << endl;
        m_log << "Wavelength: " << m_waveLength << " Amplitude: " << m_waveAmplitude << " Period: " << m_waveTime
              << " Speed: " << m_waveSpeed << endl;
        m_log << "Max up/down speed: " << speedAmplitude << endl;
        m_log << "Max up/down speed: " << m_waveSpeed * m_waveAmplitude << endl;
        m_log << "////////////////////////////////////////////////////////////////" << endl;

        fixTimeStepTravelingWave();
        break;
      }
    case 12:
      // streamwise travelling wave defined by non-plus units
      {
        m_streamwiseTravelingWave = true;
        if(!m_restart) {
          m_waveTimeStepComputed = false;
        }
        m_waveSpeed = 0.0;
        m_waveLength = 0.0;
        m_waveAmplitude = 0.0;
        m_waveCellsPerWaveLength = 1;
        if(!m_restart) {
          m_waveNoStepsPerCell = 1;
        }
        // time needs to be constant for traveling wave
        m_constantTimeStep = true;

        /*! \property
          \page propertiesFVSTRCTRD
                \section waveLength
                <code>MInt FvStructuredSolver::m_waveLength </code>\n
                default = <code> 1.0 </code>\n \n
                Wavelength of the traveling wave.\n
                Possible values are:\n
                <ul>
                <li>Float > 0.0</li>
                </ul>
                Keywords: <i>WAVE, MOVING, STRUCTURED</i>
              */
        m_waveLength = 0.0;
        if(Context::propertyExists("waveLength", m_solverId)) {
          m_waveLength = Context::getSolverProperty<MFloat>("waveLength", m_solverId, AT_, &m_waveLengthPlus);
        } else {
          mTerm(1, AT_, "Property waveLength not specified in property file");
        }

        /*! \property
          \page propertiesFVSTRCTRD
                \section waveAmplitudePressure
                <code>MInt FvStructuredSolver::m_waveAmplitudePressure </code>\n
      >>>>>>> origin/improvedPartitioning
                default = <code> 1.0 </code>\n \n
                Amplitude of the traveling wave.\n
                Possible values are:\n
                <ul>
                <li>Float > 0.0</li>
                </ul>
                Keywords: <i>WAVE, MOVING, STRUCTURED</i>
              */
        m_waveAmplitude = 0.0;
        if(Context::propertyExists("waveAmplitude", m_solverId)) {
          m_waveAmplitude = Context::getSolverProperty<MFloat>("waveAmplitude", m_solverId, AT_, &m_waveAmplitudePlus);
        } else {
          mTerm(1, AT_, "Property waveAmplitude not specified in property file");
        }

        /*! \property
          \page propertiesFVSTRCTRD
                \section waveTime
                <code>MInt FvStructuredSolver::m_waveTime </code>\n
                default = <code> 1.0 </code>\n \n
                Period time of the traveling wave.\n
                Possible values are:\n
                <ul>
                <li>Float > 0.0</li>
                </ul>
                Keywords: <i>WAVE, MOVING, STRUCTURED</i>
              */
        m_waveTime = 0.0;
        if(Context::propertyExists("waveTime", m_solverId)) {
          m_waveTime = Context::getSolverProperty<MFloat>("waveTime", m_solverId, AT_, &m_waveTimePlus);
        } else {
          mTerm(1, AT_, "Property waveTime not specified in property file");
        }


        /*! \property
          \page propertiesFVSTRCTRD
                \section waveBeginTransition
                <code>MInt FvStructuredSolver::m_waveBeginTransition </code>\n
                default = <code> 1.0 </code>\n \n
                Start of the transition from flat to wave in x-dir.\n
                Possible values are:\n
                <ul>
                <li>Float</li>
                </ul>
                Keywords: <i>WAVE, MOVING, STRUCTURED</i>
              */
        m_waveBeginTransition = 0.0;
        m_waveBeginTransition =
            Context::getSolverProperty<MFloat>("waveBeginTransition", m_solverId, AT_, &m_waveBeginTransition);

        /*! \property
          \page propertiesFVSTRCTRD
                \section waveEndTransition
                <code>MInt FvStructuredSolver::m_waveEndTransition </code>\n
                default = <code> 1.0 </code>\n \n
                End of the transition from flat to wave in x-dir.\n
                Possible values are:\n
                <ul>
                <li>Float</li>
                </ul>
                Keywords: <i>WAVE, MOVING, STRUCTURED</i>
              */
        m_waveEndTransition = 0.0;
        m_waveEndTransition =
            Context::getSolverProperty<MFloat>("waveEndTransition", m_solverId, AT_, &m_waveEndTransition);

        /*! \property
          \page propertiesFVSTRCTRD
                \section waveOutBeginTransition
                <code>MInt FvStructuredSolver::m_waveOutBeginTransition </code>\n
                default = <code> 1.0 </code>\n \n
                Start of the transition from wave to flat in x-dir.\n
                Possible values are:\n
                <ul>
                <li>Float</li>
                </ul>
                Keywords: <i>WAVE, MOVING, STRUCTURED</i>
              */
        m_waveOutBeginTransition = 1000000.0;
        m_waveOutBeginTransition =
            Context::getSolverProperty<MFloat>("waveOutBeginTransition", m_solverId, AT_, &m_waveOutBeginTransition);

        /*! \property
          \page propertiesFVSTRCTRD
                \section waveOutEndTransition
                <code>MInt FvStructuredSolver::m_waveOutEndTransition </code>\n
                default = <code> 1.0 </code>\n \n
                End of the transition from wave to flat in x-dir.\n
                Possible values are:\n
                <ul>
                <li>Float</li>
                </ul>
                Keywords: <i>WAVE, MOVING, STRUCTURED</i>
              */
        m_waveOutEndTransition = 2000000.0;
        m_waveOutEndTransition =
            Context::getSolverProperty<MFloat>("waveOutEndTransition", m_solverId, AT_, &m_waveOutEndTransition);

        /*! \property
          \page propertiesFVSTRCTRD
                \section waveYBeginTransition
                <code>MInt FvStructuredSolver::m_waveYBeginTransition </code>\n
                default = <code> 1.0 </code>\n \n
                End of the transition from wave to flat in x-dir.\n
                Possible values are:\n
                <ul>
                <li>Float</li>
                </ul>
                Keywords: <i>WAVE, MOVING, STRUCTURED</i>
              */
        m_waveYBeginTransition = 2000000.0;
        if(Context::propertyExists("waveYBeginTransition", m_solverId)) {
          m_waveYBeginTransition =
              Context::getSolverProperty<MFloat>("waveYBeginTransition", m_solverId, AT_, &m_waveYBeginTransition);
        }


        /*! \property
          \page propertiesFVSTRCTRD
                \section waveYEndTransition
                <code>MInt FvStructuredSolver::m_waveYEndTransition </code>\n
                default = <code> 1.0 </code>\n \n
                End of the transition from wave to flat in x-dir.\n
                Possible values are:\n
                <ul>
                <li>Float</li>
                </ul>
                Keywords: <i>WAVE, MOVING, STRUCTURED</i>
              */
        m_waveYEndTransition = 2000000.0;
        if(Context::propertyExists("waveYEndTransition", m_solverId)) {
          m_waveYEndTransition =
              Context::getSolverProperty<MFloat>("waveYEndTransition", m_solverId, AT_, &m_waveYEndTransition);
        }

        /*! \page propertyPage1
          \section waveTemporalTransition
          <code>MInt FvStructuredSolver::m_waveOutEndTransition </code>\n
          default = <code> 500.0 </code>\n \n
          Acoustic time for wave actuation to transiate from flat plate\n
          to fully extended.\n
          Possible values are:\n
          <ul>
          <li>Float</li>
          </ul>
          Keywords: <i>WAVE, MOVING, STRUCTURED</i>
        */
        m_waveTemporalTransition = 500.0;
        if(Context::propertyExists("waveTemporalTransition", m_solverId)) {
          m_waveTemporalTransition =
              Context::getSolverProperty<MFloat>("waveTemporalTransition", m_solverId, AT_, &m_waveTemporalTransition);
        }

        m_waveSpeed = m_waveLength / m_waveTime;
        const MFloat speedAmplitude = 2 * PI * m_waveAmplitude / m_waveTime;

        m_log << "/////////////////// TRAVELING WAVE /////////////////////////////" << endl;
        m_log << "Wavelength: " << m_waveLength << " Amplitude: " << m_waveAmplitude << " Period: " << m_waveTime
              << " Speed: " << m_waveSpeed << endl;
        m_log << "Speed amplitude: " << speedAmplitude << endl;
        m_log << "Phase speed: " << m_waveSpeed << endl;
        m_log << "////////////////////////////////////////////////////////////////" << endl;
        fixTimeStepTravelingWave();
        break;
      }
    case 14:
      // oscillating cylinder
      {
        // time needs to be constant for traveling wave
        m_constantTimeStep = true;

        /*! \property
          \page propertiesFVSTRCTRD
          \section oscAmplitude
          <code>MInt FvStructuredSolver::m_oscAmplitude </code>\n
          default = <code> 1.0 </code>\n \n
          Amplitude of the oscillating cylinder motion.\n
          Possible values are:\n
          <ul>
          <li>Float > 0.0</li>
          </ul>
          Keywords: <i>MOVING, STRUCTURED</i>
        */
        m_oscAmplitude = 0.2;
        m_oscAmplitude = Context::getSolverProperty<MFloat>("oscAmplitude", m_solverId, AT_, &m_oscAmplitude);

        m_oscSr = 0.195;
        m_oscSr = Context::getSolverProperty<MFloat>("oscSr", m_solverId, AT_, &m_oscSr);

        MFloat freqFactor = 0.8;
        freqFactor = Context::getSolverProperty<MFloat>("oscFreqFactor", m_solverId, AT_, &freqFactor);

        const MFloat freq0 = m_oscSr * PV->UInfinity / m_referenceLength;
        m_oscFreq = freq0 * freqFactor;
        break;
      }
    default: {
      stringstream errorMessage;
      errorMessage << "Grid moving method " << m_gridMovingMethod << " not implemented!" << endl;
      mTerm(1, AT_, errorMessage.str());
    }
  }

  m_grid->saveGrid();

  // now move the grid to the correct position
  if(m_restart) {
    if(m_movingGrid) {
      if(m_movingGridInitialStart) {
        // if this is an initial start of the
        // grid movement, just move to initial pos
        moveGrid(true, false);
        m_grid->saveGrid();
      } else {
        // move to last pos before restart,
        // save and move to current pos again
        // this way the grid velocity is computed
        // correctly in the BC
        m_time -= m_timeStep;
        moveGrid(true, false);
        m_grid->saveGrid();
        m_time += m_timeStep;
        moveGrid(true, false);
      }
    }
  } else {
    moveGrid(true, false);
    m_grid->saveGrid();
  }

  m_log << "Initializing moving grid methods... DONE!" << endl;
}


void FvStructuredSolver2D::moveGrid(const MBool isRestart, const MBool zeroPos) {
  TRACE();
  const MFloat pi = 4.0 * atan(1);
  const MFloat t = (isRestart) ? m_time : m_time + m_timeStep * m_RKalpha[m_RKStep];

  switch(m_gridMovingMethod) {
    case 4: // inner grid movement
    {
      const MFloat beta = 16.0l;
      const MFloat StrNum = m_wallVel; // Strouhal Number, set by wallVel
      const MFloat ver = 0.0l;         // ver can be used to translate x coordinate
      const MFloat y_max = 1.0l;       // can be used to scale

      // for Square:
      const MFloat x2 = y_max * (ver + 0.1l);
      const MFloat x3 = y_max * (ver + 0.5l);
      const MFloat x4 = y_max * (ver + 0.5l);
      const MFloat x5 = y_max * (ver + 0.9l);

      const MFloat omega = m_Ma * sqrt(PV->TInfinity) * StrNum * 2.0l * pi / y_max;
      const MFloat h = F1B2 * 0.35l * (F1 - cos(omega * t));
      const MFloat hvel = F1B2 * 0.35l * (omega * sin(omega * t));

      for(MInt j = 0; j < m_nPoints[0]; j++) {
        for(MInt i = 0; i < m_nPoints[1]; i++) {
          const MInt pointId = pointIndex(i, j);

          const MFloat x = m_grid->m_initCoordinates[0][pointId];
          const MFloat y = m_grid->m_initCoordinates[1][pointId];

          MFloat g = F0;
          if(y > x2 && y < x5 && x > x2 && x < x5) {
            g = ((y < x3) ? (1.0l + tanhl(beta * (y - (x2 + x3) / 2.0l) / y_max)) / 2.0l
                          : ((y < x4) ? 1.0l : (1.0l - tanhl(beta * (y - (x4 + x5) / 2.0l) / y_max)) / 2.0l));
          }
          m_grid->m_coordinates[0][pointId] = x * (F1 - h * g * (F1 - x));
          m_grid->m_velocity[0][pointId] = x * (-hvel * g * (F1 - x));

          g = F0;
          if(x > x2 && x < x5 && y > x2 && y < x5) {
            g = ((x < x3) ? (1.0l + tanhl(beta * (x - (x2 + x3) / 2.0l) / y_max)) / 2.0l
                          : ((x < x4) ? 1.0l : (1.0l - tanhl(beta * (x - (x4 + x5) / 2.0l) / y_max)) / 2.0l));
          }
          m_grid->m_coordinates[1][pointId] = y * (F1 - h * g * (F1 - y));
          m_grid->m_velocity[1][pointId] = y * (-hvel * g * (F1 - y));
        }
      }

      break;
    }
    case 10: {
      // traveling wave case
      MFloat t_offset = t - m_movingGridTimeOffset;
      if(zeroPos) {
        t_offset = F0;
      }

      const MFloat ytransitionLength = m_waveYEndTransition - m_waveYBeginTransition;

      MFloat fadeInFactor = 0.0;
      MFloat fadeInFactorPrime = 0.0;
      MFloat fadeInFactorPrimePrime = 0.0;
      if(t_offset < m_waveTemporalTransition) {
        fadeInFactor = (1.0 - cos(t_offset / m_waveTemporalTransition * pi)) * F1B2;
        fadeInFactorPrime = (pi / m_waveTemporalTransition) * sin(t_offset / m_waveTemporalTransition * pi) * F1B2;
        fadeInFactorPrimePrime =
            POW2(pi / m_waveTemporalTransition) * cos(t_offset / m_waveTemporalTransition * pi) * F1B2;
      } else {
        fadeInFactor = 1.0;
        fadeInFactorPrime = 0.0;
        fadeInFactorPrimePrime = 0.0;
      }

      if(zeroPos) {
        fadeInFactor = F1;
      }

      for(MInt j = 0; j < m_nPoints[0]; j++) {
        for(MInt i = 0; i < m_nPoints[1]; i++) {
          const MInt pointId = pointIndex(i, j);
          const MFloat xInit = m_grid->m_initCoordinates[0][pointId];
          const MFloat yInit = m_grid->m_initCoordinates[1][pointId];


          // To modify the activation function along the wall normal direction
          MFloat ytransitionFactor = F0;
          if(yInit <= m_waveYBeginTransition) {
            ytransitionFactor = F1;
          } else if(yInit > m_waveYBeginTransition && yInit < m_waveYEndTransition) {
            ytransitionFactor = (1 + cos((yInit - m_waveYBeginTransition) / ytransitionLength * pi)) * F1B2;
          } else {
            ytransitionFactor = F0;
          }

          // Modified activation function
          const MFloat func =
              ytransitionFactor * (m_waveAmplitude * cos((F2 * pi) / m_waveLength * (xInit - m_waveSpeed * t_offset)));
          const MFloat funcPrime = ytransitionFactor * (2 * PI * m_waveSpeed / m_waveLength) * m_waveAmplitude
                                   * sin((F2 * pi) / m_waveLength * (xInit - m_waveSpeed * t_offset));
          const MFloat funcPrimePrime = -ytransitionFactor * POW2(2 * PI * m_waveSpeed / m_waveLength) * m_waveAmplitude
                                        * cos((F2 * pi) / m_waveLength * (xInit - m_waveSpeed * t_offset));

          // Base activation function
          /* const MFloat func = transitionFactor
           * (m_waveAmplitude * cos((F2 * pi) / m_waveLength * (zPrime - m_waveSpeed * t_offset)));
           const MFloat funcPrime = transitionFactor * (2 * PI * m_waveSpeed / m_waveLength) * m_waveAmplitude
           * sin((F2 * pi) / m_waveLength * (zPrime - m_waveSpeed * t_offset));
           const MFloat funcPrimePrime = -transitionFactor * POW2(2 * PI * m_waveSpeed / m_waveLength)
           * m_waveAmplitude
           * cos((F2 * pi) / m_waveLength * (zPrime - m_waveSpeed * t_offset));
           */

          m_grid->m_coordinates[1][pointId] = func * fadeInFactor + yInit;
          m_grid->m_velocity[1][pointId] = func * fadeInFactorPrime + funcPrime * fadeInFactor;
          m_grid->m_acceleration[1][pointId] =
              funcPrimePrime * fadeInFactor + 2 * funcPrime * fadeInFactorPrime + func * fadeInFactorPrimePrime;
        }
      }

      break;
    }
    case 12: {
      // streamwise traveling wave case
      MFloat t_offset = t - m_movingGridTimeOffset;
      if(zeroPos) {
        t_offset = F0;
      }
      const MFloat transitionLength = m_waveEndTransition - m_waveBeginTransition;
      const MFloat transitionOutLength = m_waveOutEndTransition - m_waveOutBeginTransition;
      const MFloat yTransitionLength = m_waveYEndTransition - m_waveYBeginTransition;

      MFloat fadeInFactor = 0.0;
      MFloat fadeInFactorPrime = 0.0;
      MFloat fadeInFactorPrimePrime = 0.0;
      if(t_offset < m_waveTemporalTransition) {
        fadeInFactor = (1.0 - cos(t_offset / m_waveTemporalTransition * pi)) * F1B2;
        fadeInFactorPrime = (pi / m_waveTemporalTransition) * sin(t_offset / m_waveTemporalTransition * pi) * F1B2;
        fadeInFactorPrimePrime =
            POW2(pi / m_waveTemporalTransition) * cos(t_offset / m_waveTemporalTransition * pi) * F1B2;
      } else {
        fadeInFactor = 1.0;
        fadeInFactorPrime = 0.0;
        fadeInFactorPrimePrime = 0.0;
      }

      if(zeroPos) {
        fadeInFactor = F1;
      }


      for(MInt j = 0; j < m_nPoints[0]; j++) {
        for(MInt i = 0; i < m_nPoints[1]; i++) {
          const MInt pointId = pointIndex(i, j);
          const MFloat xInit = m_grid->m_initCoordinates[0][pointId];
          const MFloat yInit = m_grid->m_initCoordinates[1][pointId];

          MFloat transitionFactor = F0;
          if(xInit <= m_waveBeginTransition) {
            transitionFactor = F0;
          } else if(xInit > m_waveBeginTransition && xInit < m_waveEndTransition) {
            transitionFactor = (1 - cos((xInit - m_waveBeginTransition) / transitionLength * pi)) * F1B2;
          } else if(m_waveEndTransition <= xInit && xInit <= m_waveOutBeginTransition) {
            transitionFactor = F1;
          } else if(xInit > m_waveOutBeginTransition && xInit < m_waveOutEndTransition) {
            transitionFactor = (1 + cos((xInit - m_waveOutBeginTransition) / transitionOutLength * pi)) * F1B2;
          } else {
            transitionFactor = F0;
          }

          MFloat yTransitionFactor = F1;
          if(yInit <= m_waveYBeginTransition) {
            yTransitionFactor = F1;
          } else if(yInit > m_waveYBeginTransition && yInit < m_waveYEndTransition) {
            yTransitionFactor = (1 + cos((yInit - m_waveYBeginTransition) / yTransitionLength * pi)) * F1B2;
          } else {
            yTransitionFactor = F0;
          }

          const MFloat func = transitionFactor * yTransitionFactor
                              * (m_waveAmplitude * cos((F2 * pi) / m_waveLength * (xInit - m_waveSpeed * t_offset)));
          const MFloat funcPrime = transitionFactor * yTransitionFactor * (2 * PI * m_waveSpeed / m_waveLength)
                                   * m_waveAmplitude * sin((F2 * pi) / m_waveLength * (xInit - m_waveSpeed * t_offset));
          const MFloat funcPrimePrime = -transitionFactor * yTransitionFactor
                                        * POW2(2 * PI * m_waveSpeed / m_waveLength) * m_waveAmplitude
                                        * cos((F2 * pi) / m_waveLength * (xInit - m_waveSpeed * t_offset));

          m_grid->m_coordinates[1][pointId] = func * fadeInFactor + yInit;
          m_grid->m_velocity[1][pointId] = func * fadeInFactorPrime + funcPrime * fadeInFactor;
          m_grid->m_acceleration[1][pointId] =
              funcPrimePrime * fadeInFactor + 2 * funcPrime * fadeInFactorPrime + func * fadeInFactorPrimePrime;
        }
      }

      break;
    }
    case 14: {
      // oscillating cylinder
      MFloat t_offset = t - m_movingGridTimeOffset;
      if(zeroPos) {
        t_offset = F0;
      }

      MFloat fadeInFactor = 0;
      const MFloat timeRelaxation = 1.0;

      if(t_offset < timeRelaxation) {
        fadeInFactor = (1.0 - cos(t_offset / timeRelaxation * pi)) * F1B2;
      } else {
        fadeInFactor = 1.0;
      }

      if(zeroPos) {
        fadeInFactor = F1;
      }

      for(MInt i = 0; i < m_nPoints[1]; i++) {
        for(MInt j = 0; j < m_nPoints[0] - 1; j++) {
          const MInt pIJ = pointIndex(i, j);
          const MFloat x = m_grid->m_initCoordinates[0][pIJ];
          const MFloat y = m_grid->m_initCoordinates[1][pIJ];

          const MFloat r = sqrt(POW2(x) + POW2(y));

          MFloat spaceTransition = F1;

          if(r < 1.0) {
            spaceTransition = F0;
          } else if(r >= 1.0 && r <= 41.0) {
            spaceTransition = fabs(r - 1.0) / 40.0;
          } else {
            spaceTransition = F1;
          }

          m_grid->m_coordinates[1][pIJ] =
              y * spaceTransition
              + (1.0 - spaceTransition) * (y + fadeInFactor * m_oscAmplitude * sin(2 * PI * m_oscFreq * t_offset));
          m_grid->m_velocity[1][pIJ] =
              (1.0 - spaceTransition)
              * (fadeInFactor * m_oscAmplitude * 2 * PI * m_oscFreq * cos(2 * PI * m_oscFreq * t_offset));
          m_grid->m_acceleration[1][pIJ] =
              (1.0 - spaceTransition)
              * (-fadeInFactor * m_oscAmplitude * POW2(2 * PI * m_oscFreq) * sin(2 * PI * m_oscFreq * t_offset));
        }
      }

      break;
    }
    default: {
      mTerm(1, AT_, "Grid Moving Method not implemented!");
    }
  }

  RECORD_TIMER_START(m_timers[Timers::MGExchange]);

  if(m_gridMovingMethod != 1) {
    if(m_mgExchangeCoordinates) {
      m_grid->extrapolateGhostPointCoordinates();
    }
    extrapolateGhostPointCoordinatesBC();
  }

  if(noDomains() > 1) {
    if(m_mgExchangeCoordinates) {
      m_grid->exchangePoints(m_sndComm, m_rcvComm, PARTITION_NORMAL);
    }
  }
  RECORD_TIMER_STOP(m_timers[Timers::MGExchange]);

  RECORD_TIMER_START(m_timers[Timers::MGCellCenterCoordinates]);
  m_grid->computeCellCenterCoordinates();
  RECORD_TIMER_STOP(m_timers[Timers::MGCellCenterCoordinates]);

  RECORD_TIMER_START(m_timers[Timers::MGMetrics]);
  m_grid->computeMetrics();
  RECORD_TIMER_STOP(m_timers[Timers::MGMetrics]);

  RECORD_TIMER_START(m_timers[Timers::MGJacobian]);
  m_grid->computeJacobian();

  RECORD_TIMER_STOP(m_timers[Timers::MGJacobian]);
}

void FvStructuredSolver2D::assignBndryCells() {
  TRACE();
  m_structuredBndryCnd->assignBndryCnds();
}

void FvStructuredSolver2D::initBndryCnds() {
  TRACE();
  m_structuredBndryCnd->correctBndryCndIndices();
}

void FvStructuredSolver2D::applyBoundaryCondition() {
  TRACE();
  // treat Dirichlet and Neumann BC in one go!!!
  m_structuredBndryCnd->applyDirichletNeumannBC();
}

void FvStructuredSolver2D::computeCellCentreCoordinates() {
  // function to compute the coordinates at cell centre
  // calculated over I, J loop but changed to one array
  for(MInt j = 0; j < m_nCells[0]; ++j) {
    for(MInt i = 0; i < m_nCells[1]; ++i) {
      const MInt IJ = getPointIdFromCell(i, j);
      const MInt IP1J = getPointIdFromPoint(IJ, 1, 0);
      const MInt IJP1 = getPointIdFromPoint(IJ, 0, 1);
      const MInt IP1JP1 = getPointIdFromPoint(IJ, 1, 1);
      const MInt cellId = cellIndex(i, j);

      for(MInt dim = 0; dim < nDim; dim++) {
        // average the coordinates for cell centre data
        m_cells->coordinates[dim][cellId] = F1B4
                                            * (m_grid->m_coordinates[dim][IJ] + m_grid->m_coordinates[dim][IP1J]
                                               + m_grid->m_coordinates[dim][IJP1] + m_grid->m_coordinates[dim][IP1JP1]);
      }
    }
  }
}

MBool FvStructuredSolver2D::maxResidual() {
  TRACE();

  if(globalTimeStep % m_residualInterval != 0) return true;
  MFloat epsilon = pow(10.0, -10.0);
  m_avrgResidual = F0;
  MInt cellId = F0;
  MFloat tmpResidual = F0;
  MFloat maxResidual1 = F0;
  MInt maxResIndex[3];
  // MInt localCounter=F0;
  MFloat maxResidualOrg = F0;
  MFloat localMaxResidual = F0;
  MFloat localAvrgResidual = F0;
  MFloat accumAvrgResidual = F0;
  MFloat globalMaxResidual = F0;
  // MInt accumCounter=0;
  for(MInt dim = 0; dim < nDim; dim++) {
    maxResIndex[dim] = F0;
  }

  if(!m_localTimeStep) {
    for(MInt j = m_noGhostLayers; j < m_nCells[0] - m_noGhostLayers; j++) {
      for(MInt i = m_noGhostLayers; i < m_nCells[1] - m_noGhostLayers; i++) {
        cellId = cellIndex(i, j);
        tmpResidual = m_timeStep / (m_cfl * m_cells->cellJac[cellId]) * fabs(m_cells->rightHandSide[CV->RHO][cellId]);
        m_avrgResidual += tmpResidual;
        if(tmpResidual > maxResidual1) {
          maxResIndex[0] = i - m_noGhostLayers;
          maxResIndex[1] = j - m_noGhostLayers;
          maxResidual1 = tmpResidual;
        }
      }
    }
  } else {
    for(MInt j = m_noGhostLayers; j < m_nCells[0] - m_noGhostLayers; j++) {
      for(MInt i = m_noGhostLayers; i < m_nCells[1] - m_noGhostLayers; i++) {
        cellId = cellIndex(i, j);
        tmpResidual = m_cells->localTimeStep[cellId] / (m_cfl * m_cells->cellJac[cellId])
                      * fabs(m_cells->rightHandSide[CV->RHO][cellId]);
        m_avrgResidual += tmpResidual;
        if(tmpResidual > maxResidual1) {
          maxResIndex[0] = i - m_noGhostLayers;
          maxResIndex[1] = j - m_noGhostLayers;
          maxResidual1 = tmpResidual;
        }
      }
    }
  }

  // localCounter = counter;
  localMaxResidual = maxResidual1;
  localAvrgResidual = m_avrgResidual;
  // reset average Residual
  m_avrgResidual = F0;

  MPI_Allreduce(&localAvrgResidual, &accumAvrgResidual, 1, MPI_DOUBLE, MPI_SUM, m_StructuredComm, AT_,
                "localAvrgResidual", "accumAvrgResidual");
  MPI_Allreduce(&localMaxResidual, &globalMaxResidual, 1, MPI_DOUBLE, MPI_MAX, m_StructuredComm, AT_,
                "localMaxResidual", "globalMaxResidual");
  m_avrgResidual = accumAvrgResidual; // m_residualRcv.avrgRes;
  maxResidualOrg = globalMaxResidual;
  m_avrgResidual = m_avrgResidual / m_totalNoCells;

  // write first residuals;
  if(fabs(m_firstMaxResidual) < epsilon) {
    m_firstMaxResidual = mMax(epsilon, globalMaxResidual);
    m_firstAvrgResidual = mMax(epsilon, m_avrgResidual);
    if(approx(localMaxResidual, maxResidualOrg, m_eps)) {
      // write out values into residual file
      FILE* f_residual;
      f_residual = fopen("./Residual", "a+");
      fprintf(f_residual, "#MaxRes_1: %1.10e \n", m_firstMaxResidual);
      fprintf(f_residual, "#MaxAvgRes_1: %1.10e \n", m_firstAvrgResidual);
      fprintf(f_residual, "#iter, physTime, time, dT, wLoad, avrgRes, maxRes, blockId, i, j");
      fclose(f_residual);
    }
  }

  // normalize residuals
  globalMaxResidual = globalMaxResidual / m_firstMaxResidual;
  m_avrgResidual = (m_avrgResidual / m_firstAvrgResidual);

  if(std::isnan(m_avrgResidual)) {
    cerr << "Solution diverged, average residual is nan " << endl;
    m_log << "Solution diverged, average residual is nan " << endl;
    RECORD_TIMER_STOP(m_timers[Timers::MainLoop]);
    RECORD_TIMER_STOP(m_timers[Timers::Run]);
    saveSolverSolution(true);
    savePartitions();
    mTerm(1, AT_, "Solution diverged, average residual is nan ");
  }

  // convergence Check
  m_convergence = false;
  if(maxResidual1 < m_convergenceCriterion) {
    m_convergence = true;
  }

  // processor with the highest Residual writes out!!! saves communication;
  if(approx(localMaxResidual, maxResidualOrg,
            m_eps)) { // so only cpu with the max writes out ==> no need to communicate the max index[i]
    // write out values into residual file
    FILE* f_residual;
    f_residual = fopen("./Residual", "a+");
    fprintf(f_residual, "%d", globalTimeStep);
    fprintf(f_residual, " %f", m_physicalTime);
    fprintf(f_residual, " %f", m_time);
    fprintf(f_residual, " %f", m_timeStep);
    fprintf(f_residual, " %f", m_workload);
    fprintf(f_residual, " %1.10e", m_avrgResidual);
    fprintf(f_residual, " %1.10e", globalMaxResidual);
    fprintf(f_residual, " %d", m_blockId);
    fprintf(f_residual, " %d", m_nOffsetCells[1] + maxResIndex[0]); // i
    fprintf(f_residual, " %d", m_nOffsetCells[0] + maxResIndex[1]); // j
    fprintf(f_residual, "\n");
    fclose(f_residual);
  }

  if(maxResidual1 < m_convergenceCriterion) {
    return true;
  } else {
    return false;
  }
}

inline MFloat FvStructuredSolver2D::crossProduct(MFloat vec1[2], MFloat vec2[2]) {
  MFloat result = vec1[xsd] * vec2[ysd] - vec1[ysd] * vec2[xsd];
  return result;
}


inline MInt FvStructuredSolver2D::getPointIdFromCell(MInt i, MInt j) { return i + (j * (m_nCells[1] + 1)); }

inline MInt FvStructuredSolver2D::getPointIdFromPoint(MInt origin, MInt incI, MInt incJ) {
  return origin + incI + incJ * m_nPoints[1];
}

inline MInt FvStructuredSolver2D::getCellIdFromCell(MInt origin, MInt incI, MInt incJ) {
  return origin + incI + incJ * m_nCells[1];
}

inline MInt FvStructuredSolver2D::cellIndex(MInt i, MInt j) { return i + (j * m_nCells[1]); }

inline MInt FvStructuredSolver2D::pointIndex(MInt i, MInt j) { return i + (j * m_nPoints[1]); }

MFloat FvStructuredSolver2D::pressure(MInt cellId) { return m_cells->pvariables[PV->P][cellId]; }

void FvStructuredSolver2D::addGhostPointCoordinateValues() {
  TRACE();

  if(m_debugOutput) {
    for(MInt j = 0; j < (m_nCells[0]); j++) {
      for(MInt i = 0; i < (m_nCells[1]); i++) {
        MInt cellId = cellIndex(i, j);
        m_cells->fq[FQ->CELLID][cellId] = cellId;
        m_cells->fq[FQ->BLOCKID][cellId] = domainId();
      }
    }
  }

  // 1) extrapolate GhostPoints
  m_grid->extrapolateGhostPointCoordinates();
  // 2) communicate GhostPoints

  m_grid->exchangePoints(m_sndComm, m_rcvComm, PARTITION_NORMAL);
  m_grid->exchangePoints(m_sndComm, m_rcvComm, PERIODIC_BC);

  extrapolateGhostPointCoordinatesBC();

  computeCellCentreCoordinates();

  // MUST be done after cell center computation!!!
  m_grid->exchangePoints(m_sndComm, m_rcvComm, SINGULAR);
  m_grid->exchangePoints(m_sndComm, m_rcvComm, PERIODIC_BC_SINGULAR);

  if(m_hasSingularity > 0) {
    computeReconstructionConstantsSVD();
  }

  // 3) write the totalGridFile with GhostPoints
  if(m_savePartitionOutput) {
    m_grid->writePartitionedGrid();
  }
}


void FvStructuredSolver2D::extrapolateGhostPointCoordinatesBC() {
  for(MInt bcId = 0; bcId < (MInt)m_structuredBndryCnd->m_physicalBCMap.size(); ++bcId) {
    // all the periodic BCs are NOT included.
    // also skip the channel bc
    if(m_structuredBndryCnd->m_physicalBCMap[bcId]->BC == 2401
       || m_structuredBndryCnd->m_physicalBCMap[bcId]->BC == 2402
       || (m_structuredBndryCnd->m_physicalBCMap[bcId]->BC >= 6000
           && m_structuredBndryCnd->m_physicalBCMap[bcId]->BC < 6010)) {
      continue;
    }

    MInt* start = m_structuredBndryCnd->m_physicalBCMap[bcId]->start1;
    MInt* end = m_structuredBndryCnd->m_physicalBCMap[bcId]->end1;
    MInt index = m_structuredBndryCnd->m_physicalBCMap[bcId]->face / 2;
    MInt step = m_structuredBndryCnd->m_physicalBCMap[bcId]->face % 2;
    MInt pos[2], fix[2], mirror[2], ij[2], extendij[2];
    MInt pointId, FixPointId, MirrorPointId;

    extendij[0] = 1;
    extendij[1] = 1;
    extendij[index] = 0;

    for(ij[1] = start[1]; ij[1] < end[1] + extendij[1]; ++ij[1]) {
      for(ij[0] = start[0]; ij[0] < end[0] + extendij[0]; ++ij[0]) {
        for(MInt m = 0; m < 2; ++m) {
          if(index == m) {
            if(step == 1) {
              pos[m] = ij[m] + 1;
              fix[m] = start[m];
              mirror[m] = 2 * fix[m] - pos[m];
            } else {
              pos[m] = ij[m];
              fix[m] = end[m];
              mirror[m] = 2 * fix[m] - pos[m];
            }
          } else {
            pos[m] = ij[m];
            fix[m] = ij[m];
            mirror[m] = ij[m];
          }
        } // m

        pointId = pointIndex(pos[0], pos[1]);
        FixPointId = pointIndex(fix[0], fix[1]);
        MirrorPointId = pointIndex(mirror[0], mirror[1]);

        for(MInt dim = 0; dim < nDim; dim++) {
          m_grid->m_coordinates[dim][pointId] =
              (2 * m_grid->m_coordinates[dim][FixPointId] - m_grid->m_coordinates[dim][MirrorPointId]);
        }
      } // ij
    }

  } // bcid
}


void FvStructuredSolver2D::Muscl(MInt NotUsed(timerId)) {
  TRACE();

  if(m_movingGrid) {
    RECORD_TIMER_START(m_timers[Timers::MovingGrid]);
    if(m_RKStep == 0) {
      RECORD_TIMER_START(m_timers[Timers::MGSaveGrid]);
      m_grid->saveGrid();
      m_grid->saveCellJacobian();
      RECORD_TIMER_STOP(m_timers[Timers::MGSaveGrid]);
    }

    RECORD_TIMER_START(m_timers[Timers::MGMoveGrid]);
    moveGrid(false, false);
    RECORD_TIMER_STOP(m_timers[Timers::MGMoveGrid]);

    // compute the volume fluxes
    RECORD_TIMER_START(m_timers[Timers::MGVolumeFlux]);
    m_grid->computeDxt(m_timeStep, m_RKalpha, m_RKStep);
    RECORD_TIMER_STOP(m_timers[Timers::MGVolumeFlux]);
    RECORD_TIMER_STOP(m_timers[Timers::MovingGrid]);
  }

  RECORD_TIMER_START(m_timers[Timers::ConvectiveFlux]);
  (this->*reconstructSurfaceData)();
  RECORD_TIMER_STOP(m_timers[Timers::ConvectiveFlux]);
}


// Muscl reconstruction with Albada limiter
void FvStructuredSolver2D::MusclAlbada() {
  TRACE();

  MFloat* const* const RESTRICT flux = ALIGNED_F(m_cells->flux);
  MInt cellId = 0, IP1 = 0, IM1 = 0, IP2 = 0;
  // grid stretching factors
  MFloat DS = F0, DSP1 = F0, DSM1 = F0, DSP = F0, DSM = F0;
  // stencil identifier
  MInt IJ[2] = {1, m_nCells[1]};

  // flow variables differences
  // MFloat DQ=F0, DQP1=F0, DQM1=F0;
  MFloatScratchSpace DQ(CV->noVariables, AT_, "DQ");
  MFloatScratchSpace DQP1(CV->noVariables, AT_, "DQP1");
  MFloatScratchSpace DQM1(CV->noVariables, AT_, "DQM1");

  // left and right state
  MFloat epsLim = m_eps;
  MFloat smps = F0;
  MFloat dummy = F0, dummy1 = F0;

  MFloat pIM2 = F0, pIM1 = F0, pIP1 = F0, pIP2 = F0;
  for(MInt i = 0; i < CV->noVariables; i++) {
    DQ[i] = F0;
    DQP1[i] = F0;
    DQM1[i] = F0;
  }
  // reduce to onedimensional arrays
  MFloat* __restrict x = &m_cells->coordinates[0][0];
  MFloat* __restrict y = &m_cells->coordinates[1][0];

  MFloat phi = F0, psi = F0, vel = F0;

  /////////IMPORTANT PARAMETER
  // MFloat epsi=F1;
  // MFloat kappa=F1B3;
  /////////END IMPORTANT PARAMETER
  for(MInt dim = 0; dim < nDim; ++dim) {
    for(MInt j = m_noGhostLayers - 1; j < m_nCells[0] - m_noGhostLayers; ++j) {
      for(MInt i = m_noGhostLayers - 1; i < m_nCells[1] - m_noGhostLayers; ++i) {
        // cell ids
        cellId = cellIndex(i, j);

        IP1 = cellId + IJ[dim];
        IM1 = cellId - IJ[dim];
        IP2 = cellId + 2 * IJ[dim];

        // distances q_i+1 - q_i
        DS = sqrt(POW2(x[IP1] - x[cellId]) + POW2(y[IP1] - y[cellId]));
        // distances q_i - q_i-1
        DSM1 = sqrt(POW2(x[cellId] - x[IM1]) + POW2(y[cellId] - y[IM1]));
        DSP1 = sqrt(POW2(x[IP2] - x[IP1]) + POW2(y[IP2] - y[IP1]));
        // account for grid stretching
        // like tfs comment
        // DSP=F2*DS/POW2(DSP1+DS);
        // DSM=F2*DS/POW2(DSM1+DS);
        // like tfs
        DSP = DS / POW2(DSP1 + DS);
        DSM = DS / POW2(DSM1 + DS);


        for(MInt var = 0; var < CV->noVariables; ++var) {
          DQ[var] = m_cells->variables[var][IP1] - m_cells->variables[var][cellId];

          DQP1[var] = m_cells->variables[var][IP2] - m_cells->variables[var][IP1];

          DQM1[var] = m_cells->variables[var][cellId] - m_cells->variables[var][IM1];
          // limiter
        }
        vel = F0;
        for(MInt dim1 = 0; dim1 < nDim; ++dim1) {
          vel += POW2(m_cells->variables[CV->RHO_VV[dim1]][IM1] / m_cells->variables[CV->RHO][IM1]);
        }
        pIM2 = m_cells->variables[CV->RHO_E][IM1] - F1B2 * m_cells->variables[CV->RHO][IM1] * vel;
        vel = F0;
        for(MInt dim1 = 0; dim1 < nDim; ++dim1) {
          vel += POW2(m_cells->variables[CV->RHO_VV[dim1]][cellId] / m_cells->variables[CV->RHO][cellId]);
        }
        pIM1 = m_cells->variables[CV->RHO_E][cellId] - F1B2 * m_cells->variables[CV->RHO][cellId] * vel;

        vel = F0;
        for(MInt dim1 = 0; dim1 < nDim; ++dim1) {
          vel += POW2(m_cells->variables[CV->RHO_VV[dim1]][IP2] / m_cells->variables[CV->RHO][IP2]);
        }
        pIP2 = m_cells->variables[CV->RHO_E][IP2] - F1B2 * m_cells->variables[CV->RHO][IP2] * vel;
        vel = F0;
        for(MInt dim1 = 0; dim1 < nDim; ++dim1) {
          vel += POW2(m_cells->variables[CV->RHO_VV[dim1]][IP1] / m_cells->variables[CV->RHO][IP1]);
        }
        pIP1 = m_cells->variables[CV->RHO_E][IP1] - F1B2 * m_cells->variables[CV->RHO][IP1] * vel;

        smps = DS * DSP1;

        dummy = fabs(pIM2 - F2 * pIM1 + pIP1) / (pIM2 + F2 * pIM1 + pIP1);
        dummy1 = fabs(pIM1 - F2 * pIP1 + pIP2) / (pIM1 + F2 * pIP1 + pIP2);

        psi = mMin(F1, F6 * mMax(dummy, dummy1));
        epsLim = mMax(m_eps, pow(F1B2 * smps, F5));


        for(MInt var = 0; var < CV->noVariables; ++var) {
          phi = F1B2
                - (F1B2
                   - mMax(F0, (DQP1[var] * DQM1[var] * smps + F1B2 * epsLim)
                                  / (POW2(DQP1[var] * DS) + POW2(DQM1[var] * DSP1) + epsLim)))
                      * psi;

          m_QLeft[var] = m_cells->variables[var][cellId] + DSM * (DSM1 * DQ[var] + DS * DQM1[var]) * phi;
          m_QRight[var] = m_cells->variables[var][IP1] - DSP * (DS * DQP1[var] + DSP1 * DQ[var]) * phi;

          // PHI(IP,IM,ID)=F2-(F2-MAX(F0,(DQ(IP,ID)*DQ(IM,ID)*SMSP+EPSMP2)
          //     &     /((DQ(IP,ID)*SM)**2+(DQ(IM,ID)*SP)**2+EPSMP)))*PSI
        }


        AusmLES(m_QLeft, m_QRight, dim, cellId); // Flux balance in AUSM
      }
    }


    // FLUX BALANCE
    for(MInt v = 0; v < CV->noVariables; ++v) {
      for(MInt j = m_noGhostLayers; j < m_nCells[0] - m_noGhostLayers; ++j) {
        for(MInt i = m_noGhostLayers; i < m_nCells[1] - m_noGhostLayers; ++i) {
          const MInt I = cellIndex(i, j);
          IM1 = I - IJ[dim];
          m_cells->rightHandSide[v][I] += flux[v][IM1] - flux[v][I];
        }
      }
    }
  }
}

void FvStructuredSolver2D::MusclRANS() { m_ransSolver->Muscl(); }

template <FvStructuredSolver2D::fluxmethod ausm, MInt noVars>
void FvStructuredSolver2D::Muscl_() {
  TRACE();
  // stencil identifier
  const MUint noCells = m_noCells;
  const MInt IJK[2] = {1, m_nCells[1]};

  const MFloat* const RESTRICT x = ALIGNED_F(m_cells->coordinates[0]);
  const MFloat* const RESTRICT y = ALIGNED_F(m_cells->coordinates[1]);

  const MFloat* const RESTRICT cellVariables = ALIGNED_F(m_cells->pvariables[0]);
  MFloat* const RESTRICT cellRhs = ALIGNED_MF(m_cells->rightHandSide[0]);
  MFloat* const RESTRICT qleft = ALIGNED_MF(m_QLeft);
  MFloat* const RESTRICT qright = ALIGNED_MF(m_QRight);
  MFloat* const* const RESTRICT flux = ALIGNED_F(m_cells->flux);

  for(MInt dim = 0; dim < nDim; dim++) {
    for(MInt j = m_noGhostLayers - 1; j < m_nCells[0] - m_noGhostLayers; j++) {
      for(MInt i = m_noGhostLayers - 1; i < m_nCells[1] - m_noGhostLayers; i++) {
        const MInt I = cellIndex(i, j);
        const MInt IP1 = I + IJK[dim];
        const MInt IM1 = I - IJK[dim];
        const MInt IP2 = I + 2 * IJK[dim];

        // distances q_i+1 - q_i
        const MFloat DS = sqrt(POW2(x[IP1] - x[I]) + POW2(y[IP1] - y[I]));
        // distances q_i - q_i-1
        const MFloat DSM1 = sqrt(POW2(x[I] - x[IM1]) + POW2(y[I] - y[IM1]));
        const MFloat DSP1 = sqrt(POW2(x[IP2] - x[IP1]) + POW2(y[IP2] - y[IP1]));
        const MFloat DSP = DS / POW2(DSP1 + DS);
        const MFloat DSM = DS / POW2(DSM1 + DS);

        for(MUint v = 0; v < noVars; ++v) {
          const MUint offset = v * noCells;
          const MFloat* const RESTRICT vars = ALIGNED_F(cellVariables + offset);
          const MFloat DQ = vars[IP1] - vars[I];
          const MFloat DQP1 = vars[IP2] - vars[IP1];
          const MFloat DQM1 = vars[I] - vars[IM1];
          qleft[v] = vars[I] + DSM * (DSM1 * DQ + DS * DQM1);
          qright[v] = vars[IP1] - DSP * (DS * DQP1 + DSP1 * DQ);
        }

        (this->*ausm)(m_QLeft, m_QRight, dim, I); // Flux balance in AUSM
      }
    }

    // FLUX BALANCE
    for(MUint v = 0; v < noVars; v++) {
      for(MInt j = m_noGhostLayers; j < m_nCells[0] - m_noGhostLayers; j++) {
        for(MInt i = m_noGhostLayers; i < m_nCells[1] - m_noGhostLayers; i++) {
          const MInt I = cellIndex(i, j);
          const MInt IM1 = I - IJK[dim];
          const MUint offset = v * noCells;
          MFloat* const RESTRICT rhs = ALIGNED_F(cellRhs + offset);
          rhs[I] += flux[v][IM1] - flux[v][I];
        }
      }
    }
  }
}

// standard Ausm
template void FvStructuredSolver2D::Muscl_<&FvStructuredSolver2D::AusmLES, 5>();
template void FvStructuredSolver2D::Muscl_<&FvStructuredSolver2D::AusmLES, 6>();
template void FvStructuredSolver2D::Muscl_<&FvStructuredSolver2D::AusmLES, 7>();
// pthrc Ausm
template void FvStructuredSolver2D::Muscl_<&FvStructuredSolver2D::AusmLES_PTHRC, 5>();
template void FvStructuredSolver2D::Muscl_<&FvStructuredSolver2D::AusmLES_PTHRC, 6>();
template void FvStructuredSolver2D::Muscl_<&FvStructuredSolver2D::AusmLES_PTHRC, 7>();
// ausm dv
template void FvStructuredSolver2D::Muscl_<&FvStructuredSolver2D::AusmDV, 5>();
template void FvStructuredSolver2D::Muscl_<&FvStructuredSolver2D::AusmDV, 6>();
template void FvStructuredSolver2D::Muscl_<&FvStructuredSolver2D::AusmDV, 7>();


template <FvStructuredSolver2D::fluxmethod ausm, MInt noVars>
void FvStructuredSolver2D::MusclStretched_() {
  TRACE();

  // stencil identifier
  const MInt IJK[2] = {1, m_nCells[1]};
  const MFloat* const RESTRICT cellVariables = ALIGNED_F(m_cells->pvariables[0]);
  const MFloat* const RESTRICT cellLength = ALIGNED_F(m_cells->cellLength[0]);
  MFloat* const RESTRICT cellRhs = ALIGNED_MF(m_cells->rightHandSide[0]);
  MFloat* const RESTRICT qleft = ALIGNED_MF(m_QLeft);
  MFloat* const RESTRICT qright = ALIGNED_MF(m_QRight);
  MFloat* const* const RESTRICT flux = ALIGNED_F(m_cells->flux);
  /////////IMPORTANT PARAMETER
  const MFloat phi = F1;
  const MFloat kappa = F0; // F1B3;
  /////////END IMPORTANT PARAMETER
  for(MInt dim = 0; dim < nDim; dim++) {
    const MUint dimOffset = dim * m_noCells;
    const MFloat* const RESTRICT length = ALIGNED_F(cellLength + dimOffset);

    for(MInt j = m_noGhostLayers - 1; j < m_nCells[0] - m_noGhostLayers; j++) {
      for(MInt i = m_noGhostLayers - 1; i < m_nCells[1] - m_noGhostLayers; i++) {
        const MInt I = cellIndex(i, j);
        const MInt IP1 = I + IJK[dim];
        const MInt IM1 = I - IJK[dim];
        const MInt IP2 = I + 2 * IJK[dim];

        const MFloat rp = (length[I] + length[IP1]) / (F2 * length[I]);
        const MFloat rm = (length[I] + length[IM1]) / (F2 * length[I]);
        const MFloat f = phi / (F2 * (rp + rm));
        const MFloat f1 = (rm + kappa * phi) / rp;
        const MFloat f2 = (rp - kappa * phi) / rm;

        const MFloat rp1 = (length[IP1] + length[IP2]) / (F2 * length[IP1]);
        const MFloat rm1 = (length[IP1] + length[I]) / (F2 * length[IP1]);
        const MFloat fa = phi / (F2 * (rp1 + rm1));
        const MFloat fb = (rm1 - kappa * phi) / rp1;
        const MFloat fc = (rp1 + kappa * phi) / rm1;

        for(MUint v = 0; v < noVars; v++) {
          const MUint offset = v * m_noCells;
          const MFloat* const RESTRICT vars = ALIGNED_F(cellVariables + offset);
          // left variables
          const MFloat DQ = (vars[IP1] - vars[I]);
          const MFloat DQM1 = (vars[I] - vars[IM1]);
          qleft[v] = vars[I] + f * (f1 * DQ + f2 * DQM1);

          // right variables
          const MFloat DQP1 = (vars[IP2] - vars[IP1]);
          const MFloat DQ1 = (vars[IP1] - vars[I]);
          qright[v] = vars[IP1] - fa * (fb * DQP1 + fc * DQ1);
        }

        (this->*ausm)(m_QLeft, m_QRight, dim, I); // Flux balance in AUSM
      }
    }

    // FLUX BALANCE
    for(MUint v = 0; v < noVars; v++) {
      for(MInt j = m_noGhostLayers; j < m_nCells[0] - m_noGhostLayers; j++) {
        for(MInt i = m_noGhostLayers; i < m_nCells[1] - m_noGhostLayers; i++) {
          const MInt I = cellIndex(i, j);
          const MInt IM1 = I - IJK[dim];
          const MUint offset = v * m_noCells;
          MFloat* const RESTRICT rhs = ALIGNED_F(cellRhs + offset);
          rhs[I] += flux[v][IM1] - flux[v][I];
        }
      }
    }
  }
}
// standard Ausm
template void FvStructuredSolver2D::MusclStretched_<&FvStructuredSolver2D::AusmLES, 5>();
template void FvStructuredSolver2D::MusclStretched_<&FvStructuredSolver2D::AusmLES, 6>();
template void FvStructuredSolver2D::MusclStretched_<&FvStructuredSolver2D::AusmLES, 7>();
// pthrc
template void FvStructuredSolver2D::MusclStretched_<&FvStructuredSolver2D::AusmLES_PTHRC, 5>();
template void FvStructuredSolver2D::MusclStretched_<&FvStructuredSolver2D::AusmLES_PTHRC, 6>();
template void FvStructuredSolver2D::MusclStretched_<&FvStructuredSolver2D::AusmLES_PTHRC, 7>();


void FvStructuredSolver2D::Ausm() {
  // Ausm routines have been moved and are called from inside Muscl (better performance)
}


/**
 * \brief AUSM Central
 *
 *  Can be used for moving grids, dxt term is included
 */
// inline void FvStructuredSolver2D::AusmNew(MFloat* QLeft, MFloat* QRight, const MInt dim, const MInt I)
inline void FvStructuredSolver2D::AusmLES(MFloat* QLeft, MFloat* QRight, const MInt dim, const MInt I) {
  const MFloat gamma = m_gamma;
  const MFloat gammaMinusOne = gamma - 1.0;
  const MFloat FgammaMinusOne = F1 / gammaMinusOne;

  const MFloat surf0 = m_cells->surfaceMetrics[dim * 2 + 0][I];
  const MFloat surf1 = m_cells->surfaceMetrics[dim * 2 + 1][I];

  const MFloat dxdtau = m_cells->dxt[dim][I];

  // calculate pressure
  const MFloat PL = QLeft[PV->P];
  const MFloat UL = QLeft[PV->U];
  const MFloat VL = QLeft[PV->V];
  const MFloat RHOL = QLeft[PV->RHO];

  const MFloat PR = QRight[PV->P];
  const MFloat UR = QRight[PV->U];
  const MFloat VR = QRight[PV->V];
  const MFloat RHOR = QRight[PV->RHO];

  // compute lenght of metric vector for normalization
  const MFloat metricLength = sqrt(POW2(surf0) + POW2(surf1));
  const MFloat fMetricLength = F1 / metricLength;

  // scale by metric length to get velocity in the new basis (get normalized basis vectors)
  const MFloat UUL = ((UL * surf0 + VL * surf1) - dxdtau) * fMetricLength;


  const MFloat UUR = ((UR * surf0 + VR * surf1) - dxdtau) * fMetricLength;


  // speed of sound
  const MFloat AL = sqrt(gamma * mMax(m_eps, PL / mMax(m_eps, RHOL)));
  const MFloat AR = sqrt(gamma * mMax(m_eps, PR / mMax(m_eps, RHOR)));

  const MFloat MAL = UUL / AL;
  const MFloat MAR = UUR / AR;

  const MFloat MALR = F1B2 * (MAL + MAR);
  const MFloat PLR = F1B2 * (PL + PR);

  const MFloat RHO_AL = RHOL * AL;
  const MFloat RHO_AR = RHOR * AR;

  const MFloat PLfRHOL = PL / RHOL;
  const MFloat PRfRHOR = PR / RHOR;

  const MFloat e0 = PLfRHOL * FgammaMinusOne + 0.5 * (POW2(UL) + POW2(VL)) + PLfRHOL;
  const MFloat e1 = PRfRHOR * FgammaMinusOne + 0.5 * (POW2(UR) + POW2(VR)) + PRfRHOR;

  const MFloat RHOU = F1B2 * (MALR * (RHO_AL + RHO_AR) + fabs(MALR) * (RHO_AL - RHO_AR)) * metricLength;
  const MFloat RHOU2 = F1B2 * RHOU;
  // multiply by metric length to take surface area into account
  const MFloat AbsRHO_U2 = fabs(RHOU2);

  MFloat* const* const RESTRICT flux = ALIGNED_MF(m_cells->flux);

  flux[CV->RHO_U][I] = RHOU2 * (UL + UR) + AbsRHO_U2 * (UL - UR) + PLR * surf0;
  flux[CV->RHO_V][I] = RHOU2 * (VL + VR) + AbsRHO_U2 * (VL - VR) + PLR * surf1;
  flux[CV->RHO_E][I] = RHOU2 * (e0 + e1) + AbsRHO_U2 * (e0 - e1) + PLR * dxdtau;
  flux[CV->RHO][I] = RHOU;
}


/**
 *  \brief Same AUSM scheme as AusmLES with additional damping controlled
 *  by the 4th order pressure derivative. Pressure needs to computed
 *  beforehand.
 *
 */
inline void FvStructuredSolver2D::AusmLES_PTHRC(MFloat* QLeft, MFloat* QRight, MInt dim, MInt I) {
  const MFloat gamma = m_gamma;
  const MFloat gammaMinusOne = gamma - 1.0;
  const MFloat FgammaMinusOne = F1 / gammaMinusOne;

  const MFloat surf0 = m_cells->surfaceMetrics[dim * 2 + 0][I];
  const MFloat surf1 = m_cells->surfaceMetrics[dim * 2 + 1][I];

  const MFloat* const RESTRICT p = ALIGNED_F(m_cells->pvariables[PV->P]);

  const MFloat dxdtau = m_cells->dxt[dim][I];

  // calculate pressure
  const MFloat PL = QLeft[PV->P];
  const MFloat UL = QLeft[PV->U];
  const MFloat VL = QLeft[PV->V];
  const MFloat RHOL = QLeft[PV->RHO];

  const MFloat PR = QRight[PV->P];
  const MFloat UR = QRight[PV->U];
  const MFloat VR = QRight[PV->V];
  const MFloat RHOR = QRight[PV->RHO];

  // compute lenght of metric vector for normalization
  const MFloat metricLength = sqrt(POW2(surf0) + POW2(surf1));
  const MFloat fMetricLength = F1 / metricLength;

  // scale by metric length to get velocity in the new basis (get normalized basis vectors)
  const MFloat UUL = ((UL * surf0 + VL * surf1) - dxdtau) * fMetricLength;


  const MFloat UUR = ((UR * surf0 + VR * surf1) - dxdtau) * fMetricLength;


  // speed of sound
  const MFloat AL = sqrt(gamma * mMax(m_eps, PL / mMax(m_eps, RHOL)));
  const MFloat AR = sqrt(gamma * mMax(m_eps, PR / mMax(m_eps, RHOR)));

  const MFloat MAL = UUL / AL;
  const MFloat MAR = UUR / AR;

  const MFloat MALR = F1B2 * (MAL + MAR);

  // compute splitting pressure
  const MInt IPJK = getCellIdFromCell(I, 1, 0);
  const MInt IMJK = getCellIdFromCell(I, -1, 0);
  const MInt IP2JK = getCellIdFromCell(I, 2, 0);
  const MInt IM2JK = getCellIdFromCell(I, -2, 0);

  const MInt IJPK = getCellIdFromCell(I, 0, 1);
  const MInt IJMK = getCellIdFromCell(I, 0, -1);
  const MInt IJP2K = getCellIdFromCell(I, 0, 2);
  const MInt IJM2K = getCellIdFromCell(I, 0, -2);

  const MFloat p4I4 = F4 * (p[IPJK] + p[IMJK]) - F6 * (p[I]) - p[IP2JK] - p[IM2JK];
  const MFloat p4J4 = F4 * (p[IJPK] + p[IJMK]) - F6 * (p[I]) - p[IJP2K] - p[IJM2K];

  const MFloat cfac = 1.0 / 1.3;
  const MFloat pfac = fabs(p4I4) + fabs(p4J4);
  MFloat fac = cfac * pfac;
  fac = min(1 / 64.0, fac * 5.0);

  const MFloat PLR = PL * (F1B2 + fac * MAL) + PR * (F1B2 - fac * MAR);

  const MFloat RHO_AL = RHOL * AL;
  const MFloat RHO_AR = RHOR * AR;

  const MFloat PLfRHOL = PL / RHOL;
  const MFloat PRfRHOR = PR / RHOR;

  const MFloat e0 = PLfRHOL * FgammaMinusOne + 0.5 * (POW2(UL) + POW2(VL)) + PLfRHOL;
  const MFloat e1 = PRfRHOR * FgammaMinusOne + 0.5 * (POW2(UR) + POW2(VR)) + PRfRHOR;

  const MFloat RHOU = F1B2 * (MALR * (RHO_AL + RHO_AR) + fabs(MALR) * (RHO_AL - RHO_AR)) * metricLength;
  const MFloat RHOU2 = F1B2 * RHOU;
  // multiply by metric length to take surface area into account
  const MFloat AbsRHO_U2 = fabs(RHOU2);

  MFloat* const* const RESTRICT flux = m_cells->flux;

  flux[CV->RHO_U][I] = RHOU2 * (UL + UR) + AbsRHO_U2 * (UL - UR) + PLR * surf0;
  flux[CV->RHO_V][I] = RHOU2 * (VL + VR) + AbsRHO_U2 * (VL - VR) + PLR * surf1;
  flux[CV->RHO_E][I] = RHOU2 * (e0 + e1) + AbsRHO_U2 * (e0 - e1) + PLR * dxdtau;
  flux[CV->RHO][I] = RHOU;
}

void FvStructuredSolver2D::AusmDV(MFloat* QLeft, MFloat* QRight, const MInt dim, const MInt I) {
  const MFloat gamma = m_gamma;
  const MFloat gammaMinusOne = gamma - 1.0;
  const MFloat FgammaMinusOne = F1 / gammaMinusOne;

  const MFloat surf0 = m_cells->surfaceMetrics[dim * 2 + 0][I];
  const MFloat surf1 = m_cells->surfaceMetrics[dim * 2 + 1][I];

  // left side
  const MFloat RHOL = QLeft[PV->RHO];
  const MFloat FRHOL = F1 / RHOL;
  MFloat UL = QLeft[PV->U];
  MFloat VL = QLeft[PV->V];
  const MFloat PL = QLeft[PV->P];

  // right side
  const MFloat RHOR = QRight[PV->RHO];
  const MFloat FRHOR = F1 / RHOR;
  MFloat UR = QRight[PV->U];
  MFloat VR = QRight[PV->V];
  const MFloat PR = QRight[PV->P];

  const MFloat PLfRHOL = PL / RHOL;
  const MFloat PRfRHOR = PR / RHOR;
  const MFloat e0 = PLfRHOL * FgammaMinusOne + 0.5 * (POW2(UL) + POW2(VL)) + PLfRHOL;
  const MFloat e1 = PRfRHOR * FgammaMinusOne + 0.5 * (POW2(UR) + POW2(VR)) + PRfRHOR;


  // compute lenght of metric vector for normalization
  const MFloat DGRAD = sqrt(POW2(surf0) + POW2(surf1));
  const MFloat FDGRAD = F1 / DGRAD;

  // scale by metric length to get velocity in the new basis (get normalized basis vectors)
  const MFloat UUL = ((UL * surf0 + VL * surf1)) * FDGRAD;


  const MFloat UUR = ((UR * surf0 + VR * surf1)) * FDGRAD;

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

  const MInt IJK[2] = {1, m_nCells[1]};
  const MInt IP1 = I + IJK[dim];

  const MFloat FDV = 0.3;
  const MFloat DXDXEZ = m_cells->coordinates[0][IP1] - m_cells->coordinates[0][I];
  const MFloat DYDXEZ = m_cells->coordinates[1][IP1] - m_cells->coordinates[1][I];
  MFloat SV = 2.0 * DGRAD / (m_cells->cellJac[I] + m_cells->cellJac[IP1]) * (FDV + (F1 - FDV) * getPSI(I, dim));
  const MFloat SV1 = SV * DXDXEZ;
  const MFloat SV2 = SV * DYDXEZ;

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


  MFloat* const* const RESTRICT flux = m_cells->flux;

  flux[CV->RHO_U][I] = RHOU2 * (UL + UR) + ARHOU2 * (UL - UR) + PLR * surf0 + RHOUL * UUL2 + RHOUR * UUR2;
  flux[CV->RHO_V][I] = RHOU2 * (VL + VR) + ARHOU2 * (VL - VR) + PLR * surf1 + RHOUL * UUL3 + RHOUR * UUR3;
  flux[CV->RHO_E][I] = RHOU2 * (e0 + e1) + ARHOU2 * (e0 - e1);
  flux[CV->RHO][I] = RHOU;
}

void FvStructuredSolver2D::computeTimeStep() {
  TRACE();
  m_timeStep = 1000.0;
  const MFloat* const RESTRICT dxtx = ALIGNED_F(m_cells->dxt[0]);
  const MFloat* const RESTRICT dxty = ALIGNED_F(m_cells->dxt[1]);
  const MFloat* const* const RESTRICT metric = m_cells->cellMetrics;

  for(MInt j = m_noGhostLayers; j < m_nCells[0] - m_noGhostLayers; j++) {
    for(MInt i = m_noGhostLayers; i < m_nCells[1] - m_noGhostLayers; i++) {
      const MInt cellId = cellIndex(i, j);
      const MFloat Frho = F1 / m_cells->pvariables[PV->RHO][cellId];

      // compute the speed of sound
      const MFloat speedOfSound = sqrt(m_gamma * pressure(cellId) * Frho);

      // no need for simplified metrics, since information is already contained
      // in cell metrics
      const MFloat lenXi = sqrt(POW2(metric[0][cellId]) + POW2(metric[1][cellId]));

      const MFloat lenEt = sqrt(POW2(metric[2][cellId]) + POW2(metric[3][cellId]));

      // contravariant velocities
      MFloat U_c = F0;
      MFloat V_c = F0;

      for(MInt isd = xsd; isd < nDim; isd++) {
        U_c += m_cells->pvariables[PV->VV[isd]][cellId] * metric[xsd * nDim + isd][cellId];
        V_c += m_cells->pvariables[PV->VV[isd]][cellId] * metric[ysd * nDim + isd][cellId];
      }

      // subtract grid velocity
      U_c -= dxtx[cellId];
      V_c -= dxty[cellId];

      U_c = fabs(U_c);
      V_c = fabs(V_c);

      // has area information in it due to metric terms
      const MFloat eigenvalue = U_c + V_c + speedOfSound * (lenXi + lenEt);

      // divide volume information (jacobian) through area to get characteristic length for CFL
      const MFloat deltaT = m_cfl * m_cells->cellJac[cellId] / eigenvalue;
      if(m_localTimeStep) {
        m_cells->localTimeStep[cellId] = deltaT;
        m_timeStep = F1;
        m_timeRef = F1;
      } else {
        // TODO_SS labels:FV if jacobian is negative, then timeStep will be negative and the following mMin won't
        // work
        m_timeStep = mMin(m_timeStep, deltaT);
      }
    }
  }
}

void FvStructuredSolver2D::updateSpongeLayer() {
  TRACE();
  if(m_useSponge) m_structuredBndryCnd->updateSpongeLayer();
}

MBool FvStructuredSolver2D::rungeKuttaStep() {
  TRACE();
  const MInt noVars = CV->noVariables;
  const MUint noCells = m_noCells;
  const MFloat rkAlpha = m_RKalpha[m_RKStep];
  const MFloat rkFactor = rkAlpha * m_timeStep;

#ifdef MAIA_EXTRA_DEBUG
  savePartitions(); // testing only
#endif

  MFloat* const RESTRICT oldVars = ALIGNED_F(m_cells->oldVariables[0]);
  MFloat* const RESTRICT vars = ALIGNED_F(m_cells->variables[0]);
  MFloat* const RESTRICT oldCellJac = ALIGNED_MF(m_cells->oldCellJac);
  const MFloat* const RESTRICT cellJac = ALIGNED_MF(m_cells->cellJac);
  const MFloat* const RESTRICT rhs = ALIGNED_MF(m_cells->rightHandSide[0]);

  // set old variables
  if(m_RKStep == 0) {
    for(MInt v = 0; v < noVars; v++) {
      const MUint offset = v * noCells;
      MFloat* const RESTRICT oldCellVars = ALIGNED_F(oldVars + offset);
      const MFloat* const RESTRICT cellVars = ALIGNED_F(vars + offset);
      for(MInt cellId = 0; cellId < m_noCells; cellId++) {
        oldCellVars[cellId] = cellVars[cellId];
      }
    }
  }


  switch(m_rungeKuttaOrder) {
    case 2: {
      if(m_localTimeStep) {
        for(MInt v = 0; v < noVars; v++) {
          const MUint cellOffset = v * noCells;
          MFloat* const RESTRICT cellVars = vars + cellOffset;
          const MFloat* const RESTRICT oldCellVars = oldVars + cellOffset;
          const MFloat* const RESTRICT cellRhs = rhs + cellOffset;

          for(MInt j = m_noGhostLayers; j < m_nCells[0] - m_noGhostLayers; j++) {
            for(MInt i = m_noGhostLayers; i < m_nCells[1] - m_noGhostLayers; i++) {
              const MInt cellId = cellIndex(i, j);
              const MFloat factor = (rkAlpha * m_cells->localTimeStep[cellId]) / m_cells->cellJac[cellId];
              cellVars[cellId] = oldCellVars[cellId] + factor * cellRhs[cellId];
            }
          }
        }
      } else if(m_movingGrid) {
        for(MInt v = 0; v < noVars; v++) {
          const MUint cellOffset = v * noCells;
          MFloat* const RESTRICT cellVars = vars + cellOffset;
          const MFloat* const RESTRICT oldCellVars = oldVars + cellOffset;
          const MFloat* const RESTRICT cellRhs = rhs + cellOffset;

          for(MInt j = m_noGhostLayers; j < m_nCells[0] - m_noGhostLayers; j++) {
            for(MInt i = m_noGhostLayers; i < m_nCells[1] - m_noGhostLayers; i++) {
              const MInt cellId = cellIndex(i, j);
              cellVars[cellId] =
                  (oldCellVars[cellId] * oldCellJac[cellId] + rkFactor * cellRhs[cellId]) / cellJac[cellId];
            }
          }
        }
      } else {
        for(MInt v = 0; v < noVars; v++) {
          const MUint cellOffset = v * noCells;
          MFloat* const RESTRICT cellVars = vars + cellOffset;
          const MFloat* const RESTRICT oldCellVars = oldVars + cellOffset;
          const MFloat* const RESTRICT cellRhs = rhs + cellOffset;

          for(MInt j = m_noGhostLayers; j < m_nCells[0] - m_noGhostLayers; j++) {
            for(MInt i = m_noGhostLayers; i < m_nCells[1] - m_noGhostLayers; i++) {
              const MInt cellId = cellIndex(i, j);
              const MFloat factor = (rkAlpha * m_timeStep) / m_cells->cellJac[cellId];
              cellVars[cellId] = oldCellVars[cellId] + factor * cellRhs[cellId];
            }
          }
        }
      }
      break;
    }
    case 3: {
      for(MInt v = 0; v < noVars; v++) {
        const MUint cellOffset = v * noCells;
        MFloat* const RESTRICT cellVars = vars + cellOffset;
        const MFloat* const RESTRICT oldCellVars = oldVars + cellOffset;
        const MFloat* const RESTRICT cellRhs = rhs + cellOffset;

        for(MInt j = m_noGhostLayers; j < m_nCells[0] - m_noGhostLayers; j++) {
          for(MInt i = m_noGhostLayers; i < m_nCells[1] - m_noGhostLayers; i++) {
            const MInt cellId = cellIndex(i, j);
            const MFloat factor = (rkAlpha * m_timeStep) / m_cells->cellJac[cellId];
            cellVars[cellId] =
                rkAlpha * cellVars[cellId] + (F1 - rkAlpha) * oldCellVars[cellId] - factor * cellRhs[cellId];
          }
        }
      }
      break;
    }
    default: {
      stringstream errorMessage;
      errorMessage << "Given RungeKutta Order " << m_rungeKuttaOrder << " not implemented! " << endl;
      mTerm(1, AT_, errorMessage.str());
    }
  }

  ++m_RKStep;

  if(m_RKStep == m_noRKSteps) {
    m_physicalTime += m_timeStep * m_timeRef;
    m_time += m_timeStep;

    m_RKStep = 0;

    return true;
  } else {
    return false;
  }
}

void FvStructuredSolver2D::viscousFlux() { (this->*viscFluxMethod)(); }

void FvStructuredSolver2D::viscousFluxRANS() { m_ransSolver->viscousFluxRANS(); }


template <MBool twoEqRans>
void FvStructuredSolver2D::viscousFluxLES() {
  TRACE();
  const MFloat rPrLam = F1 / m_Pr;
  const MFloat rPrTurb = F1 / 0.9;
  const MFloat rRe = F1 / m_Re0;
  const MFloat gammaMinusOne = m_gamma - 1.0;
  const MFloat FgammaMinusOne = F1 / gammaMinusOne;

  const MFloat* const RESTRICT u = &m_cells->pvariables[PV->U][0];
  const MFloat* const RESTRICT v = &m_cells->pvariables[PV->V][0];
  const MFloat* const RESTRICT p = &m_cells->pvariables[PV->P][0];
  const MFloat* const RESTRICT rho = &m_cells->pvariables[PV->RHO][0];
  MFloat* const RESTRICT T = &m_cells->temperature[0];
  MFloat* const RESTRICT muLam = &m_cells->fq[FQ->MU_L][0];
  MFloat* const RESTRICT muTurb = &m_cells->fq[FQ->MU_T][0]; // this is zero for LES

  MFloat* const* const RESTRICT eflux = ALIGNED_MF(m_cells->eFlux);
  MFloat* const* const RESTRICT fflux = ALIGNED_MF(m_cells->fFlux);
  MFloat* const* const RESTRICT vflux = ALIGNED_MF(m_cells->viscousFlux);

  // only relevant for 2-eq Rans model (Waiting for if constexpr)
  const MFloat* const RESTRICT TKE = (twoEqRans) ? ALIGNED_F(m_cells->pvariables[PV->RANS_VAR[0]]) : nullptr;
  MFloat TKEcorner = 0;

  for(MInt j = m_noGhostLayers - 1; j < m_nCells[0] - m_noGhostLayers + 1; j++) {
    for(MInt ii = m_noGhostLayers - 1; ii < m_nCells[1] - m_noGhostLayers + 1; ii++) {
      const MInt I = cellIndex(ii, j);
      T[I] = m_gamma * p[I] / rho[I];
      muLam[I] = SUTHERLANDLAW(T[I]);
    }
  }

  MFloat tau1, tau2, tau4;
  MFloat dTdx, dTdy;

  for(MInt j = m_noGhostLayers - 1; j < m_nCells[0] - m_noGhostLayers; j++) {
    for(MInt i = m_noGhostLayers - 1; i < m_nCells[1] - m_noGhostLayers; i++) {
      // get the adjacent cells;
      const MInt IJ = cellIndex(i, j);
      const MInt IPJ = cellIndex((i + 1), j);
      const MInt IPJP = cellIndex((i + 1), (j + 1));
      const MInt IJP = cellIndex(i, (j + 1));

      const MFloat dudxi = F1B2 * (u[IPJP] + u[IPJ] - u[IJP] - u[IJ]);
      const MFloat dudet = F1B2 * (u[IPJP] + u[IJP] - u[IPJ] - u[IJ]);

      const MFloat dvdxi = F1B2 * (v[IPJP] + v[IPJ] - v[IJP] - v[IJ]);
      const MFloat dvdet = F1B2 * (v[IPJP] + v[IJP] - v[IPJ] - v[IJ]);

      const MFloat dTdxi = F1B2 * (T[IPJP] + T[IPJ] - T[IJP] - T[IJ]);
      const MFloat dTdet = F1B2 * (T[IPJP] + T[IJP] - T[IPJ] - T[IJ]);

      const MFloat uAvg = F1B4 * (u[IJP] + u[IPJP] + u[IJ] + u[IPJ]);
      const MFloat vAvg = F1B4 * (v[IJP] + v[IPJP] + v[IJ] + v[IPJ]);

      const MFloat muLamAvg = F1B4 * (muLam[IJP] + muLam[IPJP] + muLam[IJ] + muLam[IPJ]);
      const MFloat muTurbAvg = F1B4 * (muTurb[IJP] + muTurb[IPJP] + muTurb[IJ] + muTurb[IPJ]);

      if(twoEqRans) {
        if(m_rans2eq_mode == "production")
          TKEcorner =
              -2 / 3 * F1B4 * (rho[IJP] * TKE[IJP] + rho[IPJP] * TKE[IPJP] + rho[IJ] * TKE[IJ] + rho[IPJ] * TKE[IPJ]);
      }

      const MFloat cornerMetrics[4] = {m_cells->cornerMetrics[0][IJ], m_cells->cornerMetrics[1][IJ],
                                       m_cells->cornerMetrics[2][IJ], m_cells->cornerMetrics[3][IJ]};

      // compute tau1 = 2 du/dx - 2/3 ( du/dx + dv/dy + dw/dz )

      // tau_xx = 4/3*( du/dxi * dxi/dx + du/deta * deta/dx + du/dzeta * dzeta/dx )
      //            - 2/3*( dv/dxi * dxi/dy + dv/deta * deta/dy + dv/dzeta * dzeta/dy)
      //            - 2/3*( dw/dxi * dxi/dz + dw/deta * deta/dz + dw/dzeta * dzeta/dz )
      tau1 = F4B3 * (dudxi * cornerMetrics[xsd * 2 + xsd] + dudet * cornerMetrics[ysd * 2 + xsd]) -

             F2B3 * (dvdxi * cornerMetrics[xsd * 2 + ysd] + dvdet * cornerMetrics[ysd * 2 + ysd]);

      // compute tau2 = du/dy + dv/dx

      // tau_xy = du/dxi * dxi/dy + du/deta * deta/dy + du/dzeta * dzeta/dy
      //        + dv/dxi * dxi/dx + dv/deta * deta/dx + dv/dzeta * dzeta/dx
      tau2 = dudxi * cornerMetrics[xsd * 2 + ysd] + dudet * cornerMetrics[ysd * 2 + ysd] +

             dvdxi * cornerMetrics[xsd * 2 + xsd] + dvdet * cornerMetrics[ysd * 2 + xsd];

      // compute tau4 = 2 dv/dy - 2/3 ( du/dx + dv/dy + dw/dz )

      // tau_yy = 4/3*( dv/dxi * dxi/dy + dv/deta * deta/dy + dv/dzeta * dzeta/dy )
      //            - 2/3*( du/dxi * dxi/dx + du/deta * deta/dx + du/dzeta * dzeta/dx)
      //            - 2/3*( dw/dxi * dxi/dz + dw/deta * deta/dz + dw/dzeta * dzeta/dz )
      tau4 = F4B3 * (dvdxi * cornerMetrics[xsd * 2 + ysd] + dvdet * cornerMetrics[ysd * 2 + ysd]) -

             F2B3 * (dudxi * cornerMetrics[xsd * 2 + xsd] + dudet * cornerMetrics[ysd * 2 + xsd]);


      dTdx = dTdxi * cornerMetrics[xsd * 2 + xsd] + dTdet * cornerMetrics[ysd * 2 + xsd];

      dTdy = dTdxi * cornerMetrics[xsd * 2 + ysd] + dTdet * cornerMetrics[ysd * 2 + ysd];

      const MFloat fJac = 1.0 / m_cells->cornerJac[IJ];
      const MFloat mueOverRe = rRe * fJac * (muLamAvg + muTurbAvg);
      tau1 = mueOverRe * tau1 + TKEcorner;
      tau2 *= mueOverRe;
      tau4 = mueOverRe * tau4 + TKEcorner;

      const MFloat muCombined = FgammaMinusOne * rRe * fJac * (rPrLam * muLamAvg + rPrTurb * muTurbAvg);
      const MFloat qx = muCombined * dTdx + uAvg * tau1 + vAvg * tau2;
      const MFloat qy = muCombined * dTdy + uAvg * tau2 + vAvg * tau4;

      // efluxes
      eflux[0][IJ] = tau1 * cornerMetrics[xsd * 2 + xsd] + tau2 * cornerMetrics[xsd * 2 + ysd];

      eflux[1][IJ] = tau2 * cornerMetrics[xsd * 2 + xsd] + tau4 * cornerMetrics[xsd * 2 + ysd];

      eflux[2][IJ] = qx * cornerMetrics[xsd * 2 + xsd] + qy * cornerMetrics[xsd * 2 + ysd];

      // ffluxes
      fflux[0][IJ] = tau1 * cornerMetrics[ysd * 2 + xsd] + tau2 * cornerMetrics[ysd * 2 + ysd];

      fflux[1][IJ] = tau2 * cornerMetrics[ysd * 2 + xsd] + tau4 * cornerMetrics[ysd * 2 + ysd];

      fflux[2][IJ] = qx * cornerMetrics[ysd * 2 + xsd] + qy * cornerMetrics[ysd * 2 + ysd];
    }
  }


  // viscous flux correction for the singular points
  // m_hasSingularity=0 means no singular points in this block, otherwise do flux correction
  if(m_hasSingularity > 0) {
    viscousFluxCorrection();
  }

  for(MInt var = 0; var < 3; ++var) {
    for(MInt j = m_noGhostLayers; j < m_nCells[0] - m_noGhostLayers; ++j) {
      for(MInt i = m_noGhostLayers - 1; i < m_nCells[1] - m_noGhostLayers; ++i) {
        const MInt IJ = cellIndex(i, j);
        const MInt IJM = cellIndex(i, (j - 1));

        vflux[0][IJ] = F1B2 * (eflux[var][IJ] + eflux[var][IJM]);
      }
    }

    for(MInt j = m_noGhostLayers - 1; j < m_nCells[0] - m_noGhostLayers; ++j) {
      for(MInt i = m_noGhostLayers; i < m_nCells[1] - m_noGhostLayers; ++i) {
        const MInt IJ = cellIndex(i, j);
        const MInt IMJ = cellIndex((i - 1), j);

        vflux[1][IJ] = F1B2 * (fflux[var][IJ] + fflux[var][IMJ]);
      }
    }

    for(MInt j = m_noGhostLayers; j < m_nCells[0] - m_noGhostLayers; ++j) {
      for(MInt i = m_noGhostLayers; i < m_nCells[1] - m_noGhostLayers; ++i) {
        const MInt IJ = cellIndex(i, j);
        const MInt IMJ = cellIndex((i - 1), j);
        const MInt IJM = cellIndex(i, (j - 1));
        m_cells->rightHandSide[var][IJ] += (vflux[0][IJ] - vflux[0][IMJ] + vflux[1][IJ] - vflux[1][IJM]);
      }
    }
  }
}
template void FvStructuredSolver2D::viscousFluxLES<true>();


/*template <MBool twoEqRans>
  void FvStructuredSolver2D::viscousFluxLESCompact()
  {
  TRACE();
  const MInt noCells = m_noCells;
  const MFloat rPrLam = F1/m_Pr;
  const MFloat rPrTurb = F1/0.9;
  const MFloat rRe = F1/m_Re0;
  const MFloat gammaMinusOne = m_gamma - 1.0;
  const MFloat FgammaMinusOne = F1 / gammaMinusOne;

  MFloat* const RESTRICT u = &m_cells->pvariables[PV->U][0];
  MFloat* const RESTRICT v = &m_cells->pvariables[PV->V][0];
  MFloat* const RESTRICT p = &m_cells->pvariables[PV->P][0];
  MFloat* const RESTRICT rho = &m_cells->pvariables[PV->RHO][0];
  MFloat* const RESTRICT T = &m_cells->temperature[0];
  MFloat* const RESTRICT muLam = &m_cells->fq[FQ->MU_L][0];
  MFloat* const RESTRICT muTurb = &m_cells->fq[FQ->MU_T][0]; //this is zero for LES

  MFloat *const RESTRICT eflux= ALIGNED_MF(m_cells->eFlux);
  MFloat *const RESTRICT fflux= ALIGNED_MF(m_cells->fFlux);

  for(MInt j=m_noGhostLayers-1; j<m_nCells[0]-m_noGhostLayers+1; j++) {
  for(MInt ii=m_noGhostLayers-1; ii<m_nCells[1]-m_noGhostLayers+1; ii++) {
  const MInt I=cellIndex(ii,j);
  T[I] = m_gamma*p[I]/rho[I];
  muLam[I] = SUTHERLANDLAW(T[I]);
  }
  }

  // Save some aux. variables at the corners
  for(MInt j=m_noGhostLayers-1; j<m_nCells[0]-m_noGhostLayers; j++) {
  for(MInt i=m_noGhostLayers-1; i<m_nCells[1]-m_noGhostLayers; i++) {
  const MInt IJ   = cellIndex(i,j);
  const MInt IPJ  = cellIndex((i+1),j);
  const MInt IPJP = cellIndex((i+1),(j+1));
  const MInt IJP  = cellIndex(i,(j+1));

  const MFloat dudxi=F1B2*(u[IPJP]+u[IPJ]-u[IJP]-u[IJ]);
  const MFloat dudet=F1B2*(u[IPJP]+u[IJP]-u[IPJ]-u[IJ]);

  const MFloat dvdxi=F1B2*(v[IPJP]+v[IPJ]-v[IJP]-v[IJ]);
  const MFloat dvdet=F1B2*(v[IPJP]+v[IJP]-v[IPJ]-v[IJ]);

  const MFloat dTdxi=F1B2*(T[IPJP]+T[IPJ]-T[IJP]-T[IJ]);
  const MFloat dTdet=F1B2*(T[IPJP]+T[IJP]-T[IPJ]-T[IJ]);

  // Temporarily save corner derivatives
  fflux[ 0*noCells+IJ ] = dudxi;
  fflux[ 1*noCells+IJ ] = dvdxi;
  fflux[ 2*noCells+IJ ] = dTdxi;

  eflux[ 0*noCells+IJ ] = dudet;
  eflux[ 1*noCells+IJ ] = dvdet;
  eflux[ 2*noCells+IJ ] = dTdet;
  }
  }

  // xi-Flux
  for(MInt j=m_noGhostLayers; j<m_nCells[0]-m_noGhostLayers; j++) {
  for(MInt i=m_noGhostLayers-1; i<m_nCells[1]-m_noGhostLayers; i++) {
  const MInt IJ   = cellIndex(i,j);
  const MInt IPJ  = cellIndex(i+1,j);
  const MInt IJM  = cellIndex(i,j-1);
  //      const MFloat surf[nDim*nDim] = {0.5*(m_cells->cornerMetrics[0*noCells+IJ] +
  m_cells->cornerMetrics[0*noCells+IJM]),
  //                                        0.5*(m_cells->cornerMetrics[1*noCells+IJ] +
  m_cells->cornerMetrics[1*noCells+IJM]),
  //                                        0.5*(m_cells->cornerMetrics[2*noCells+IJ] +
  m_cells->cornerMetrics[2*noCells+IJM]),
  //                                        0.5*(m_cells->cornerMetrics[3*noCells+IJ] +
  m_cells->cornerMetrics[3*noCells+IJM])}; const MFloat *const surf = m_cells->surfaceMetrics[IJ]; const MFloat
  invSurfJac = 2.0/(m_cells->cornerJac[IJ]+m_cells->cornerJac[IJM]); //TODO_SS labels:FV Is this the right way to go?

  const MFloat dudxi = u[IPJ]-u[IJ];
  const MFloat dvdxi = v[IPJ]-v[IJ];
  const MFloat dTdxi = T[IPJ]-T[IJ];
  const MFloat dudet = 0.5*(eflux[0*noCells+IJ]+eflux[0*noCells+IJM]);
  const MFloat dvdet = 0.5*(eflux[1*noCells+IJ]+eflux[1*noCells+IJM]);
  const MFloat dTdet = 0.5*(eflux[2*noCells+IJ]+eflux[2*noCells+IJM]);

  const MFloat uAvg = F1B2*(u[IJ]+u[IPJ]);
  const MFloat vAvg = F1B2*(v[IJ]+v[IPJ]);

  const MFloat muLamAvg = F1B2*(muLam[IJ]+muLam[IPJ]);
  const MFloat muTurbAvg = F1B2*(muTurb[IJ]+muTurb[IPJ]);

  // compute tau1 = 2 du/dx - 2/3 ( du/dx + dv/dy + dw/dz )

  // tau_xx = 4/3*( du/dxi * dxi/dx + du/deta * deta/dx + du/dzeta * dzeta/dx )
  //            - 2/3*( dv/dxi * dxi/dy + dv/deta * deta/dy + dv/dzeta * dzeta/dy)
  //            - 2/3*( dw/dxi * dxi/dz + dw/deta * deta/dz + dw/dzeta * dzeta/dz )
  MFloat tau1 = F4B3 * ( dudxi * surf[ xsd * 2 + xsd ]   +
  dudet * surf[ ysd * 2 + xsd ] ) -

  F2B3 * ( dvdxi * surf[ xsd * 2 + ysd ]   +
  dvdet * surf[ ysd * 2 + ysd ] );

  // compute tau2 = du/dy + dv/dx

  // tau_xy = du/dxi * dxi/dy + du/deta * deta/dy + du/dzeta * dzeta/dy
  //        + dv/dxi * dxi/dx + dv/deta * deta/dx + dv/dzeta * dzeta/dx
  MFloat tau2 = dudxi * surf[ xsd * 2 + ysd ] +
  dudet * surf[ ysd * 2 + ysd ] +

  dvdxi * surf[ xsd * 2 + xsd ] +
  dvdet * surf[ ysd * 2 + xsd ];

  // compute tau4 = 2 dv/dy - 2/3 ( du/dx + dv/dy + dw/dz )

  // tau_yy = 4/3*( dv/dxi * dxi/dy + dv/deta * deta/dy + dv/dzeta * dzeta/dy )
  //            - 2/3*( du/dxi * dxi/dx + du/deta * deta/dx + du/dzeta * dzeta/dx)
  //            - 2/3*( dw/dxi * dxi/dz + dw/deta * deta/dz + dw/dzeta * dzeta/dz )
  MFloat tau4 = F4B3 * ( dvdxi * surf[ xsd * 2 + ysd ] +
  dvdet * surf[ ysd * 2 + ysd ] ) -

  F2B3 * ( dudxi * surf[ xsd * 2 + xsd ] +
  dudet * surf[ ysd * 2 + xsd ]);

  const MFloat dTdx = dTdxi * surf[ xsd * 2 + xsd ] +
  dTdet * surf[ ysd * 2 + xsd ];

  const MFloat dTdy = dTdxi * surf[ xsd * 2 + ysd ] +
  dTdet * surf[ ysd * 2 + ysd ];

  const MFloat mueOverRe = rRe*invSurfJac*(muLamAvg + muTurbAvg);
  tau1=mueOverRe*tau1;
  tau2*=mueOverRe;
  tau4=mueOverRe*tau4;

  const MFloat muCombined=FgammaMinusOne*rRe*invSurfJac*(rPrLam*muLamAvg + rPrTurb*muTurbAvg);
  const MFloat qx=muCombined*dTdx+uAvg*tau1+vAvg*tau2;
  const MFloat qy=muCombined*dTdy+uAvg*tau2+vAvg*tau4;

  //efluxes
  eflux[ 0*noCells+IJM ]     = tau1 * surf[ xsd * 2 + xsd ] +
  tau2 * surf[ xsd * 2 + ysd ];

  eflux[ 1*noCells+IJM ] = tau2 * surf[ xsd * 2 + xsd ] +
  tau4 * surf[ xsd * 2 + ysd ];

  eflux[ 2*noCells+IJM ] = qx * surf[ xsd * 2 + xsd ] +
  qy * surf[ xsd * 2 + ysd ];
  }
  }

  // eta-Flux
  for(MInt j=m_noGhostLayers-1; j<m_nCells[0]-m_noGhostLayers; j++) {
    for(MInt i=m_noGhostLayers; i<m_nCells[1]-m_noGhostLayers; i++) {
      const MInt IJ   = cellIndex(i,j);
      const MInt IJP  = cellIndex(i,j+1);
      const MInt IMJ  = cellIndex(i-1,j);
      //      const MFloat surf[nDim*nDim] = {0.5*(m_cells->cornerMetrics[0*noCells+IJ] +
      m_cells->cornerMetrics[0*noCells+IMJ]),
                      //                                        0.5*(m_cells->cornerMetrics[1*noCells+IJ] +
                      m_cells->cornerMetrics[1*noCells+IMJ]),
                    //                                        0.5*(m_cells->cornerMetrics[2*noCells+IJ] +
    m_cells->cornerMetrics[2*noCells+IMJ]),
//                                        0.5*(m_cells->cornerMetrics[3*noCells+IJ] +
 m_cells->cornerMetrics[3*noCells+IMJ])}; const MFloat *const surf = m_cells->surfaceMetrics[IJ]; const MFloat
 invSurfJac = 2.0/(m_cells->cornerJac[IJ]+m_cells->cornerJac[IMJ]); //TODO_SS labels:FV Is this the right way to go?

const MFloat dudet = (u[IJP]-u[IJ]);
const MFloat dvdet = (v[IJP]-v[IJ]);
const MFloat dTdet = (T[IJP]-T[IJ]);
const MFloat dudxi = 0.5*(fflux[0*noCells+IJ]+fflux[0*noCells+IMJ]);
const MFloat dvdxi = 0.5*(fflux[1*noCells+IJ]+fflux[1*noCells+IMJ]);
const MFloat dTdxi = 0.5*(fflux[2*noCells+IJ]+fflux[2*noCells+IMJ]);

const MFloat uAvg = F1B2*(u[IJ]+u[IJP]);
const MFloat vAvg = F1B2*(v[IJ]+v[IJP]);

const MFloat muLamAvg = F1B2*(muLam[IJ]+muLam[IJP]);
const MFloat muTurbAvg = F1B2*(muTurb[IJ]+muTurb[IJP]);

// compute tau1 = 2 du/dx - 2/3 ( du/dx + dv/dy + dw/dz )

// tau_xx = 4/3*( du/dxi * dxi/dx + du/deta * deta/dx + du/dzeta * dzeta/dx )
//            - 2/3*( dv/dxi * dxi/dy + dv/deta * deta/dy + dv/dzeta * dzeta/dy)
//            - 2/3*( dw/dxi * dxi/dz + dw/deta * deta/dz + dw/dzeta * dzeta/dz )
MFloat tau1 = F4B3 * ( dudxi * surf[ xsd * 2 + xsd ]   +
dudet * surf[ ysd * 2 + xsd ] ) -

F2B3 * ( dvdxi * surf[ xsd * 2 + ysd ]   +
dvdet * surf[ ysd * 2 + ysd ] );

// compute tau2 = du/dy + dv/dx

// tau_xy = du/dxi * dxi/dy + du/deta * deta/dy + du/dzeta * dzeta/dy
//        + dv/dxi * dxi/dx + dv/deta * deta/dx + dv/dzeta * dzeta/dx
MFloat tau2 = dudxi * surf[ xsd * 2 + ysd ] +
dudet * surf[ ysd * 2 + ysd ] +

dvdxi * surf[ xsd * 2 + xsd ] +
dvdet * surf[ ysd * 2 + xsd ];

// compute tau4 = 2 dv/dy - 2/3 ( du/dx + dv/dy + dw/dz )

// tau_yy = 4/3*( dv/dxi * dxi/dy + dv/deta * deta/dy + dv/dzeta * dzeta/dy )
//            - 2/3*( du/dxi * dxi/dx + du/deta * deta/dx + du/dzeta * dzeta/dx)
//            - 2/3*( dw/dxi * dxi/dz + dw/deta * deta/dz + dw/dzeta * dzeta/dz )
MFloat tau4 = F4B3 * ( dvdxi * surf[ xsd * 2 + ysd ] +
dvdet * surf[ ysd * 2 + ysd ] ) -

F2B3 * ( dudxi * surf[ xsd * 2 + xsd ] +
dudet * surf[ ysd * 2 + xsd ]);

const MFloat dTdx = dTdxi * surf[ xsd * 2 + xsd ] +
dTdet * surf[ ysd * 2 + xsd ];

const MFloat dTdy = dTdxi * surf[ xsd * 2 + ysd ] +
dTdet * surf[ ysd * 2 + ysd ];

const MFloat mueOverRe = rRe*invSurfJac*(muLamAvg + muTurbAvg);
tau1=mueOverRe*tau1;
tau2*=mueOverRe;
tau4=mueOverRe*tau4;

const MFloat muCombined=FgammaMinusOne*rRe*invSurfJac*(rPrLam*muLamAvg + rPrTurb*muTurbAvg);
const MFloat qx=muCombined*dTdx+uAvg*tau1+vAvg*tau2;
const MFloat qy=muCombined*dTdy+uAvg*tau2+vAvg*tau4;

//ffluxes
fflux[ 0*noCells+IMJ ] = tau1 * surf[ ysd * 2 + xsd ] + tau2 * surf[ ysd * 2 + ysd ];

fflux[ 1*noCells+IMJ ] = tau2 * surf[ ysd * 2 + xsd ] + tau4 * surf[ ysd * 2 + ysd ];

fflux[ 2*noCells+IMJ ] = qx * surf[ ysd * 2 + xsd ] +  qy * surf[ ysd * 2 + ysd ];

}
}

//TODO_SS labels:FV
for(MInt var = 0; var < 3; ++var) {
for(MInt j=m_noGhostLayers; j<m_nCells[0]-m_noGhostLayers; ++j) {
for(MInt i=m_noGhostLayers; i<m_nCells[1]-m_noGhostLayers; ++i) {
const MInt IJ    = cellIndex(i,j);
const MInt IMJ   = cellIndex(i-1,j);
const MInt IJM   = cellIndex(i,j-1);
const MInt IMJM  = cellIndex(i-1,j-1);

m_cells->rightHandSide[var][IJ]+= (eflux[var*noCells+IJM] - eflux[var*noCells+IMJM] +
fflux[var*noCells+IMJ] - fflux[var*noCells+IMJM]);

}
}
}
}
template void FvStructuredSolver2D::viscousFluxLESCompact<true>();
*/

template <MBool twoEqRans>
void FvStructuredSolver2D::viscousFluxLESCompact() {
  TRACE();
  const MFloat rPrLam = F1 / m_Pr;
  const MFloat rPrTurb = F1 / 0.9;
  const MFloat rRe = F1 / m_Re0;
  const MFloat gammaMinusOne = m_gamma - 1.0;
  const MFloat FgammaMinusOne = F1 / gammaMinusOne;

  MFloat* const RESTRICT u = &m_cells->pvariables[PV->U][0];
  MFloat* const RESTRICT v = &m_cells->pvariables[PV->V][0];
  MFloat* const RESTRICT p = &m_cells->pvariables[PV->P][0];
  MFloat* const RESTRICT rho = &m_cells->pvariables[PV->RHO][0];
  MFloat* const RESTRICT T = &m_cells->temperature[0];
  MFloat* const RESTRICT muLam = &m_cells->fq[FQ->MU_L][0];
  MFloat* const RESTRICT muTurb = &m_cells->fq[FQ->MU_T][0]; // this is zero for LES

  MFloat* const* const RESTRICT eflux = ALIGNED_MF(m_cells->eFlux);
  MFloat* const* const RESTRICT fflux = ALIGNED_MF(m_cells->fFlux);
  MFloat* const* const RESTRICT dT = ALIGNED_MF(m_cells->dT);

  // only relevant for 2-eq Rans model (Waiting for if constexpr)
  const MFloat* const RESTRICT TKE = (twoEqRans) ? ALIGNED_F(m_cells->pvariables[PV->RANS_VAR[0]]) : nullptr;
  MFloat TKEcorner = 0;

  for(MInt j = m_noGhostLayers - 1; j < m_nCells[0] - m_noGhostLayers + 1; j++) {
    for(MInt ii = m_noGhostLayers - 1; ii < m_nCells[1] - m_noGhostLayers + 1; ii++) {
      const MInt I = cellIndex(ii, j);
      T[I] = m_gamma * p[I] / rho[I];
      muLam[I] = SUTHERLANDLAW(T[I]);
    }
  }

  // Save some aux. variables at the corners
  for(MInt j = m_noGhostLayers - 1; j < m_nCells[0] - m_noGhostLayers; j++) {
    for(MInt i = m_noGhostLayers - 1; i < m_nCells[1] - m_noGhostLayers; i++) {
      const MInt IJ = cellIndex(i, j);
      const MInt IPJ = cellIndex((i + 1), j);
      const MInt IPJP = cellIndex((i + 1), (j + 1));
      const MInt IJP = cellIndex(i, (j + 1));

      const MFloat invCornerJac = F1 / m_cells->cornerJac[IJ];
      const MFloat cornerMetrics[nDim * nDim] = {
          m_cells->cornerMetrics[0][IJ] * invCornerJac, m_cells->cornerMetrics[1][IJ] * invCornerJac,
          m_cells->cornerMetrics[2][IJ] * invCornerJac, m_cells->cornerMetrics[3][IJ] * invCornerJac};

      const MFloat dudxi = F1B2 * (u[IPJP] + u[IPJ] - u[IJP] - u[IJ]);
      const MFloat dudet = F1B2 * (u[IPJP] + u[IJP] - u[IPJ] - u[IJ]);

      const MFloat dvdxi = F1B2 * (v[IPJP] + v[IPJ] - v[IJP] - v[IJ]);
      const MFloat dvdet = F1B2 * (v[IPJP] + v[IJP] - v[IPJ] - v[IJ]);

      const MFloat dTdxi = F1B2 * (T[IPJP] + T[IPJ] - T[IJP] - T[IJ]);
      const MFloat dTdet = F1B2 * (T[IPJP] + T[IJP] - T[IPJ] - T[IJ]);

      eflux[0][IJ] = dudet * cornerMetrics[ysd * 2 + ysd] + dvdet * cornerMetrics[ysd * 2 + xsd];
      eflux[1][IJ] = dudet * cornerMetrics[ysd * 2 + xsd];
      eflux[2][IJ] = dvdet * cornerMetrics[ysd * 2 + ysd];

      fflux[0][IJ] = dudxi * cornerMetrics[xsd * 2 + ysd] + dvdxi * cornerMetrics[xsd * 2 + xsd];
      fflux[1][IJ] = dudxi * cornerMetrics[xsd * 2 + xsd];
      fflux[2][IJ] = dvdxi * cornerMetrics[xsd * 2 + ysd];

      dT[0][IJ] = dTdet * cornerMetrics[ysd * 2 + xsd];
      dT[1][IJ] = dTdet * cornerMetrics[ysd * 2 + ysd];
      dT[2][IJ] = dTdxi * cornerMetrics[xsd * 2 + xsd];
      dT[3][IJ] = dTdxi * cornerMetrics[xsd * 2 + ysd];
    }
  }

  if(m_hasSingularity) {
    viscousFluxCompactCorrection();
  }

  // xi-Flux
  for(MInt j = m_noGhostLayers; j < m_nCells[0] - m_noGhostLayers; j++) {
    for(MInt i = m_noGhostLayers - 1; i < m_nCells[1] - m_noGhostLayers; i++) {
      const MInt IJ = cellIndex(i, j);
      const MInt IPJ = cellIndex(i + 1, j);
      const MInt IJM = cellIndex(i, j - 1);
      const MFloat surf0 = m_cells->surfaceMetrics[0][IJ];
      const MFloat surf1 = m_cells->surfaceMetrics[1][IJ];
      // TODO_SS labels:FV here surface jac (m_cells->surfJac) verwenden!!!
      const MFloat invSurfJac =
          2.0 / (m_cells->cornerJac[IJ] + m_cells->cornerJac[IJM]); // TODO_SS labels:FV Is this the right way to go?

      const MFloat dudxi = (u[IPJ] - u[IJ]) * invSurfJac;
      const MFloat dvdxi = (v[IPJ] - v[IJ]) * invSurfJac;
      const MFloat dTdxi = (T[IPJ] - T[IJ]) * invSurfJac;

      const MFloat uAvg = F1B2 * (u[IJ] + u[IPJ]);
      const MFloat vAvg = F1B2 * (v[IJ] + v[IPJ]);

      const MFloat muLamAvg = F1B2 * (muLam[IJ] + muLam[IPJ]);
      const MFloat muTurbAvg = F1B2 * (muTurb[IJ] + muTurb[IPJ]);

      // The underscores in the variable names should indicate, that these are just
      // one contribution to the full, e.g. dudx
      const MFloat dudx_ = F1B2 * (eflux[1][IJ] + eflux[1][IJM]);
      const MFloat dvdy_ = F1B2 * (eflux[2][IJ] + eflux[2][IJM]);
      const MFloat S12_ = F1B2 * (eflux[0][IJ] + eflux[0][IJM]);
      const MFloat dTdx_ = F1B2 * (dT[0][IJ] + dT[0][IJM]);
      const MFloat dTdy_ = F1B2 * (dT[1][IJ] + dT[1][IJM]);

      if(twoEqRans) {
        if(m_rans2eq_mode == "production") TKEcorner = -2 / 3 * F1B2 * (rho[IJ] * TKE[IJ] + rho[IPJ] * TKE[IPJ]);
      }

      // compute tau1 = 2 du/dx - 2/3 ( du/dx + dv/dy + dw/dz )

      // tau_xx = 4/3*( du/dxi * dxi/dx + du/deta * deta/dx + du/dzeta * dzeta/dx )
      //            - 2/3*( dv/dxi * dxi/dy + dv/deta * deta/dy + dv/dzeta * dzeta/dy)
      //            - 2/3*( dw/dxi * dxi/dz + dw/deta * deta/dz + dw/dzeta * dzeta/dz )
      MFloat tau1 = F4B3 * (dudxi * surf0 + dudx_) - F2B3 * (dvdxi * surf1 + dvdy_);

      // compute tau2 = du/dy + dv/dx

      // tau_xy = du/dxi * dxi/dy + du/deta * deta/dy + du/dzeta * dzeta/dy
      //        + dv/dxi * dxi/dx + dv/deta * deta/dx + dv/dzeta * dzeta/dx
      MFloat tau2 = dudxi * surf1 + dvdxi * surf0 + S12_;

      // compute tau4 = 2 dv/dy - 2/3 ( du/dx + dv/dy + dw/dz )

      // tau_yy = 4/3*( dv/dxi * dxi/dy + dv/deta * deta/dy + dv/dzeta * dzeta/dy )
      //            - 2/3*( du/dxi * dxi/dx + du/deta * deta/dx + du/dzeta * dzeta/dx)
      //            - 2/3*( dw/dxi * dxi/dz + dw/deta * deta/dz + dw/dzeta * dzeta/dz )
      MFloat tau4 = F4B3 * (dvdxi * surf1 + dvdy_) - F2B3 * (dudxi * surf0 + dudx_);

      const MFloat dTdx = dTdxi * surf0 + dTdx_;
      const MFloat dTdy = dTdxi * surf1 + dTdy_;

      const MFloat mueOverRe = rRe * (muLamAvg + muTurbAvg);
      tau1 = mueOverRe * tau1 + TKEcorner;
      tau2 *= mueOverRe;
      tau4 = mueOverRe * tau4 + TKEcorner;

      const MFloat muCombined = FgammaMinusOne * rRe * (rPrLam * muLamAvg + rPrTurb * muTurbAvg);
      const MFloat qx = muCombined * dTdx + uAvg * tau1 + vAvg * tau2;
      const MFloat qy = muCombined * dTdy + uAvg * tau2 + vAvg * tau4;

      MFloat limiterVisc = 1.0;
      if(m_limiterVisc) {
        // TODO_SS labels:FV limiter should preserve conservativity
        const MFloat mu_ref = F1B2 * (m_cells->fq[FQ->LIMITERVISC][IJ] + m_cells->fq[FQ->LIMITERVISC][IPJ]);
        limiterVisc = std::min(1.0, mu_ref / std::abs(mueOverRe)); // TODO_SS labels:FV evtl: scale with Re
      }

      // efluxes
      eflux[0][IJM] = (tau1 * surf0 + tau2 * surf1) * limiterVisc;
      eflux[1][IJM] = (tau2 * surf0 + tau4 * surf1) * limiterVisc;
      eflux[2][IJM] = (qx * surf0 + qy * surf1) * limiterVisc; // TODO_SS labels:FV aufpassen mit limiter
    }
  }

  // eta-Flux
  for(MInt j = m_noGhostLayers - 1; j < m_nCells[0] - m_noGhostLayers; j++) {
    for(MInt i = m_noGhostLayers; i < m_nCells[1] - m_noGhostLayers; i++) {
      const MInt IJ = cellIndex(i, j);
      const MInt IJP = cellIndex(i, j + 1);
      const MInt IMJ = cellIndex(i - 1, j);
      const MFloat surf0 = m_cells->surfaceMetrics[2 + 0][IJ];
      const MFloat surf1 = m_cells->surfaceMetrics[2 + 1][IJ];
      // TODO_SS labels:FV here surface jac (m_cells->surfJac) verwenden!!!
      const MFloat invSurfJac =
          2.0 / (m_cells->cornerJac[IJ] + m_cells->cornerJac[IMJ]); // TODO_SS labels:FV Is this the right way to go?

      const MFloat dudet = (u[IJP] - u[IJ]) * invSurfJac;
      const MFloat dvdet = (v[IJP] - v[IJ]) * invSurfJac;
      const MFloat dTdet = (T[IJP] - T[IJ]) * invSurfJac;

      const MFloat uAvg = F1B2 * (u[IJ] + u[IJP]);
      const MFloat vAvg = F1B2 * (v[IJ] + v[IJP]);

      const MFloat muLamAvg = F1B2 * (muLam[IJ] + muLam[IJP]);
      const MFloat muTurbAvg = F1B2 * (muTurb[IJ] + muTurb[IJP]);

      // The underscores in the variable names should indicate, that these are just
      // one contribution to the full, e.g. dudx
      const MFloat dudx_ = F1B2 * (fflux[1][IJ] + fflux[1][IMJ]);
      const MFloat dvdy_ = F1B2 * (fflux[2][IJ] + fflux[2][IMJ]);
      const MFloat S12_ = F1B2 * (fflux[0][IJ] + fflux[0][IMJ]);
      const MFloat dTdx_ = F1B2 * (dT[2][IJ] + dT[2][IMJ]);
      const MFloat dTdy_ = F1B2 * (dT[3][IJ] + dT[3][IMJ]);

      if(twoEqRans) {
        if(m_rans2eq_mode == "production") TKEcorner = -2 / 3 * F1B2 * (rho[IJ] * TKE[IJ] + rho[IJP] * TKE[IJP]);
      }

      // compute tau1 = 2 du/dx - 2/3 ( du/dx + dv/dy + dw/dz )

      // tau_xx = 4/3*( du/dxi * dxi/dx + du/deta * deta/dx + du/dzeta * dzeta/dx )
      //            - 2/3*( dv/dxi * dxi/dy + dv/deta * deta/dy + dv/dzeta * dzeta/dy)
      //            - 2/3*( dw/dxi * dxi/dz + dw/deta * deta/dz + dw/dzeta * dzeta/dz )
      MFloat tau1 = F4B3 * (dudx_ + dudet * surf0) - F2B3 * (dvdy_ + dvdet * surf1);

      // compute tau2 = du/dy + dv/dx

      // tau_xy = du/dxi * dxi/dy + du/deta * deta/dy + du/dzeta * dzeta/dy
      //        + dv/dxi * dxi/dx + dv/deta * deta/dx + dv/dzeta * dzeta/dx
      MFloat tau2 = dudet * surf1 + dvdet * surf0 + S12_;

      // compute tau4 = 2 dv/dy - 2/3 ( du/dx + dv/dy + dw/dz )

      // tau_yy = 4/3*( dv/dxi * dxi/dy + dv/deta * deta/dy + dv/dzeta * dzeta/dy )
      //            - 2/3*( du/dxi * dxi/dx + du/deta * deta/dx + du/dzeta * dzeta/dx)
      //            - 2/3*( dw/dxi * dxi/dz + dw/deta * deta/dz + dw/dzeta * dzeta/dz )
      MFloat tau4 = F4B3 * (dvdy_ + dvdet * surf1) - F2B3 * (dudx_ + dudet * surf0);

      const MFloat dTdx = dTdet * surf0 + dTdx_;
      const MFloat dTdy = dTdet * surf1 + dTdy_;

      const MFloat mueOverRe = rRe * (muLamAvg + muTurbAvg);
      tau1 = mueOverRe * tau1 + TKEcorner;
      tau2 *= mueOverRe;
      tau4 = mueOverRe * tau4 + TKEcorner;

      const MFloat muCombined = FgammaMinusOne * rRe * (rPrLam * muLamAvg + rPrTurb * muTurbAvg);
      const MFloat qx = muCombined * dTdx + uAvg * tau1 + vAvg * tau2;
      const MFloat qy = muCombined * dTdy + uAvg * tau2 + vAvg * tau4;

      MFloat limiterVisc = 1.0;
      if(m_limiterVisc) {
        // TODO_SS labels:FV limiter should preserve conservativity
        const MFloat mu_ref = F1B2 * (m_cells->fq[FQ->LIMITERVISC][IJ] + m_cells->fq[FQ->LIMITERVISC][IJP]);
        limiterVisc = std::min(1.0, mu_ref / std::abs(mueOverRe)); // TODO_SS labels:FV evtl: scale with Re
      }

      // ffluxes
      fflux[0][IMJ] = (tau1 * surf0 + tau2 * surf1) * limiterVisc;
      fflux[1][IMJ] = (tau2 * surf0 + tau4 * surf1) * limiterVisc;
      fflux[2][IMJ] = (qx * surf0 + qy * surf1) * limiterVisc; // TODO_SS labels:FV aufpassen mit limiter
    }
  }

  for(MInt var = 0; var < 3; ++var) {
    for(MInt j = m_noGhostLayers; j < m_nCells[0] - m_noGhostLayers; ++j) {
      for(MInt i = m_noGhostLayers; i < m_nCells[1] - m_noGhostLayers; ++i) {
        const MInt IJ = cellIndex(i, j);
        const MInt IMJ = cellIndex(i - 1, j);
        const MInt IJM = cellIndex(i, j - 1);
        const MInt IMJM = cellIndex(i - 1, j - 1);

        MFloat limiterVisc = 1.0;
        /*        if (m_limiterVisc) {
  // TODO_SS labels:FV limiter should preserve conservativity
  const MFloat mu_ref = m_cells->fq[FQ->LIMITERVISC][IJ];
  // TODO_SS labels:FV evtl: scale with Re
  limiterVisc = std::min(1.0, mu_ref/std::abs(rRe*(muLam[IJ]+muTurb[IJ])));
  }*/
        m_cells->rightHandSide[var][IJ] +=
            (eflux[var][IJM] - eflux[var][IMJM] + (fflux[var][IMJ] - fflux[var][IMJM])) * limiterVisc;
      }
    }
  }
}
template void FvStructuredSolver2D::viscousFluxLESCompact<true>();


template <MBool twoEqRans>
void FvStructuredSolver2D::viscousFluxCorrection() {
  mTerm(1, "I don't think that this correction is ok, becuase of the cornerMetrics at the singular point!");
  const MFloat rPr = F1 / m_Pr;
  const MFloat rPrTurb = F1 / 0.9;
  const MFloat rRe = F1 / m_Re0;
  const MFloat gammaMinusOne = m_gamma - 1.0;
  const MFloat FgammaMinusOne = F1 / gammaMinusOne;

  const MFloat* const RESTRICT u_ = &m_cells->variables[PV->U][0];
  const MFloat* const RESTRICT v_ = &m_cells->variables[PV->V][0];
  const MFloat* const RESTRICT p_ = &m_cells->variables[PV->P][0];
  const MFloat* const RESTRICT rho_ = &m_cells->variables[PV->RHO][0];
  const MFloat* const RESTRICT muTurb_ = &m_cells->fq[FQ->MU_T][0]; // this is zero for LES

  MFloat* const* const RESTRICT eflux = &m_cells->eFlux[0];
  MFloat* const* const RESTRICT fflux = &m_cells->fFlux[0];

  // only relevant for 2-eq Rans model (Waiting for if constexpr)
  const MFloat* const RESTRICT TKE_ = (twoEqRans) ? ALIGNED_F(m_cells->pvariables[PV->RANS_VAR[0]]) : nullptr;

  //  MInt dim = 0;
  MInt start[nDim], end[nDim], nghbr[10];
  MInt len1[nDim];
  //  MInt totalCells;

  for(MInt i = 0; i < m_hasSingularity; ++i) {
    // only correct for bc 6000 not for bc 4000-5000
    if(m_singularity[i].BC == -6000) {
      //      totalCells=1;
      for(MInt j = 0; j < nDim; j++) {
        len1[j] = m_singularity[i].end[j] - m_singularity[i].start[j];
        //        if(len1[j]!=0)  totalCells*=len1[j];
      }

      for(MInt n = 0; n < nDim; ++n) {
        if(m_singularity[i].end[n] - m_singularity[i].start[n] > 1) {
          mTerm(1, "In 2D not possible!");
          //          dim=n;
          // start[n]=m_singularity[i].start[n]+1;
          start[n] = m_singularity[i].start[n] + 1;
          end[n] = m_singularity[i].end[n] - 1;
        } else {
          start[n] = m_singularity[i].start[n];
          end[n] = m_singularity[i].end[n];
        }
      }

      MFloat u[10], v[10], T[10], muTurb[10], TKE[10];
      MFloat U, V, t, muTurb_corner, TKEcorner, dudx, dudy, dvdx, dvdy, dTdx, dTdy;
      MFloat tau1, tau2, tau4;
      MFloat mueOverRe, mue, mueH;
      for(MInt jj = start[1]; jj < end[1]; ++jj) {
        for(MInt ii = start[0]; ii < end[0]; ++ii) {
          MInt count = 0;
          //            MInt temp[nDim]{};
          MInt IJK = cellIndex(ii + m_singularity[i].Viscous[0], jj + m_singularity[i].Viscous[1]);

          const MFloat cornerMetrics[nDim * nDim] = {m_cells->cornerMetrics[0][IJK], m_cells->cornerMetrics[1][IJK],
                                                     m_cells->cornerMetrics[2][IJK], m_cells->cornerMetrics[3][IJK]};

          //            temp[dim]=1;
          nghbr[count++] = cellIndex(ii, jj);
          //            nghbr[count++]=cellIndex(ii+temp[0],jj+temp[1],kk+temp[2]);

          for(MInt m = 0; m < m_singularity[i].Nstar - 1; ++m) {
            MInt* change = m_singularity[i].displacement[m];
            nghbr[count++] = cellIndex(ii + change[0], jj + change[1]);
            //              nghbr[count++]=cellIndex(ii+temp[0]+change[0],jj+temp[1]+change[1],kk+temp[2]+change[2]);
          }

          if(count != m_singularity[i].Nstar) {
            cout << "what the hell! it is wrong!!! count=" << count << " Nstar=" << m_singularity[i].Nstar << endl;
          }

          for(MInt m = 0; m < m_singularity[i].Nstar; ++m) {
            u[m] = u_[nghbr[m]];
            v[m] = v_[nghbr[m]];
            //              T[m]=(m_gamma*gammaMinusOne*(rhoE[nghbr[m]]-F1B2*rho[nghbr[m]]*(POW2(u[m])+POW2(v[m]))))/rho[nghbr[m]];
            T[m] = m_gamma * p_[nghbr[m]] / rho_[nghbr[m]];
            muTurb[m] = muTurb_[nghbr[m]];
            if(twoEqRans && m_rans2eq_mode == "production") TKE[m] = rho_[nghbr[m]] * TKE_[nghbr[m]];
          }

          U = F0;
          V = F0;
          t = F0;
          muTurb_corner = F0;
          TKEcorner = F0;
          dudx = F0;
          dudy = F0;
          dvdx = F0;
          dvdy = F0;
          dTdx = F0;
          dTdy = F0;

          MInt id2 = ii - start[0] + (jj - start[1]) * len1[0];

          for(MInt n = 0; n < count; n++) {
            MInt ID = id2 * count + n;
            U += m_singularity[i].ReconstructionConstants[0][ID] * u[n];
            dudx += m_singularity[i].ReconstructionConstants[1][ID] * u[n];
            dudy += m_singularity[i].ReconstructionConstants[2][ID] * u[n];

            V += m_singularity[i].ReconstructionConstants[0][ID] * v[n];
            dvdx += m_singularity[i].ReconstructionConstants[1][ID] * v[n];
            dvdy += m_singularity[i].ReconstructionConstants[2][ID] * v[n];

            t += m_singularity[i].ReconstructionConstants[0][ID] * T[n];
            dTdx += m_singularity[i].ReconstructionConstants[1][ID] * T[n];
            dTdy += m_singularity[i].ReconstructionConstants[2][ID] * T[n];

            muTurb_corner += m_singularity[i].ReconstructionConstants[0][ID] * muTurb[n];

            if(twoEqRans && m_rans2eq_mode == "production")
              TKEcorner += m_singularity[i].ReconstructionConstants[0][ID] * TKE[n];
          }

          //            cout << globalTimeStep << "(" << m_RKStep << ") dom=" << domainId() << " x|y=" <<
          //            m_cells->coordinates[0][IJK] << "|" << m_cells->coordinates[1][IJK] << " U=" << U << " V=" << V
          //            << " t=" << t << " cornerMetrics=" << cornerMetrics[0] << "|" << cornerMetrics[1] << "|" <<
          //            cornerMetrics[2] << "|" << cornerMetrics[3] << endl;

          tau1 = 2 * dudx - 2 / 3 * (dudx + dvdy);
          tau2 = dudy + dvdx;
          tau4 = 2 * dvdy - 2 / 3 * (dudx + dvdy);

          mue = SUTHERLANDLAW(t);
          TKEcorner *= -2 / 3;
          mueOverRe = (mue + muTurb_corner) * rRe;
          tau1 = mueOverRe * tau1 + TKEcorner;
          tau2 *= mueOverRe;
          tau4 = mueOverRe * tau4 + TKEcorner;
          mueH = FgammaMinusOne * rRe * (rPr * mue + rPrTurb * muTurb_corner);

          const MFloat qx = mueH * dTdx + U * tau1 + V * tau2;
          const MFloat qy = mueH * dTdy + U * tau2 + V * tau4;

          // efluxes
          eflux[0][IJK] = tau1 * cornerMetrics[xsd * 2 + xsd] + tau2 * cornerMetrics[xsd * 2 + ysd];
          eflux[1][IJK] = tau2 * cornerMetrics[xsd * 2 + xsd] + tau4 * cornerMetrics[xsd * 2 + ysd];
          eflux[2][IJK] = qx * cornerMetrics[xsd * 2 + xsd] + qy * cornerMetrics[xsd * 2 + ysd];

          // ffluxes
          fflux[0][IJK] = tau1 * cornerMetrics[ysd * 2 + xsd] + tau2 * cornerMetrics[ysd * 2 + ysd];
          fflux[1][IJK] = tau2 * cornerMetrics[ysd * 2 + xsd] + tau4 * cornerMetrics[ysd * 2 + ysd];
          fflux[2][IJK] = qx * cornerMetrics[ysd * 2 + xsd] + qy * cornerMetrics[ysd * 2 + ysd];
        }
      }
    }
  }
}


template <MBool twoEqRans>
void FvStructuredSolver2D::viscousFluxCompactCorrection() {
  MFloat* RESTRICT u_ = &m_cells->pvariables[PV->U][0];
  MFloat* RESTRICT v_ = &m_cells->pvariables[PV->V][0];
  MFloat* RESTRICT p_ = &m_cells->pvariables[PV->P][0];
  MFloat* RESTRICT rho_ = &m_cells->pvariables[PV->RHO][0];

  MFloat* const* const RESTRICT eflux = ALIGNED_MF(m_cells->eFlux);
  MFloat* const* const RESTRICT fflux = ALIGNED_MF(m_cells->fFlux);
  MFloat* const* const RESTRICT dT = ALIGNED_MF(m_cells->dT);

  //  MInt dim = 0;
  MInt start[nDim], end[nDim], nghbr[10];
  MInt len1[nDim];
  //  MInt totalCells;

  for(MInt i = 0; i < m_hasSingularity; ++i) {
    // only correct for bc 6000 not for bc 4000-5000
    if(m_singularity[i].BC == -6000) {
      //      totalCells=1;
      for(MInt j = 0; j < nDim; j++) {
        len1[j] = m_singularity[i].end[j] - m_singularity[i].start[j];
        //        if(len1[j]!=0)  totalCells*=len1[j];
      }

      for(MInt n = 0; n < nDim; ++n) {
        if(m_singularity[i].end[n] - m_singularity[i].start[n] > 1) {
          mTerm(1, "In 2D not possible!");
          //          dim=n;
          // start[n]=m_singularity[i].start[n]+1;
          start[n] = m_singularity[i].start[n] + 1;
          end[n] = m_singularity[i].end[n] - 1;
        } else {
          start[n] = m_singularity[i].start[n];
          end[n] = m_singularity[i].end[n];
        }
      }

      MFloat u[10], v[10], T[10];
      MFloat U, V, t;
      for(MInt jj = start[1]; jj < end[1]; ++jj) {
        for(MInt ii = start[0]; ii < end[0]; ++ii) {
          MInt count = 0;
          //            MInt temp[nDim]{};
          MInt IJ_ = cellIndex(ii + m_singularity[i].Viscous[0], jj + m_singularity[i].Viscous[1]);

          // dxidx, dxidy, detadx, detady
          const MFloat surfaceMetrics[nDim * nDim] = {
              m_cells->surfaceMetricsSingularity[0][i], m_cells->surfaceMetricsSingularity[1][i],
              m_cells->surfaceMetricsSingularity[2][i], m_cells->surfaceMetricsSingularity[2][i]};

          //            temp[dim]=1;
          const MInt IJ = cellIndex(ii, jj);
          nghbr[count++] = IJ;
          //            nghbr[count++]=cellIndex(ii+temp[0],jj+temp[1],kk+temp[2]);

          for(MInt m = 0; m < m_singularity[i].Nstar - 1; ++m) {
            MInt* change = m_singularity[i].displacement[m];
            nghbr[count++] = cellIndex(ii + change[0], jj + change[1]);
            //              nghbr[count++]=cellIndex(ii+temp[0]+change[0],jj+temp[1]+change[1],kk+temp[2]+change[2]);
          }

          if(count != m_singularity[i].Nstar) {
            cout << "what the hell! it is wrong!!! count=" << count << " Nstar=" << m_singularity[i].Nstar << endl;
          }

          for(MInt m = 0; m < m_singularity[i].Nstar; ++m) {
            u[m] = u_[nghbr[m]];
            v[m] = v_[nghbr[m]];
            //              T[m]=(m_gamma*gammaMinusOne*(rhoE[nghbr[m]]-F1B2*rho[nghbr[m]]*(POW2(u[m])+POW2(v[m]))))/rho[nghbr[m]];
            T[m] = m_gamma * p_[nghbr[m]] / rho_[nghbr[m]];
          }

          U = F0;
          V = F0;
          t = F0;

          MInt id2 = ii - start[0] + (jj - start[1]) * len1[0];

          for(MInt n = 0; n < count; n++) {
            MInt ID = id2 * count + n;
            U += m_singularity[i].ReconstructionConstants[0][ID] * u[n];
            V += m_singularity[i].ReconstructionConstants[0][ID] * v[n];
            t += m_singularity[i].ReconstructionConstants[0][ID] * T[n];
          }

          //            cout << globalTimeStep << "(" << m_RKStep << ") dom=" << domainId() << " x|y=" <<
          //            m_cells->coordinates[0][IJK] << "|" << m_cells->coordinates[1][IJK] << " U=" << U << " V=" << V
          //            << " t=" << t << " cornerMetrics=" << cornerMetrics[0] << "|" << cornerMetrics[1] << "|" <<
          //            cornerMetrics[2] << "|" << cornerMetrics[3] << endl;

          const MInt sign_xi = 2 * m_singularity[i].Viscous[0] + 1;  //-1,+1
          const MInt sign_eta = 2 * m_singularity[i].Viscous[1] + 1; //-1,+1

          const MInt IPMJ = getCellIdFromCell(IJ, sign_xi, 0);
          const MInt IJPM = getCellIdFromCell(IJ, 0, sign_eta);

          const MFloat dudet = sign_eta * 0.5 * (U - 0.5 * (u[0] + u_[IPMJ]));
          const MFloat dvdet = sign_eta * 0.5 * (V - 0.5 * (v[0] + v_[IPMJ]));
          const MFloat T1 = m_gamma * p_[IPMJ] / rho_[IPMJ];
          const MFloat dTdet = sign_eta * 0.5 * (t - 0.5 * (T[0] + T1));

          const MFloat dudxi = sign_xi * 0.5 * (U - 0.5 * (u[0] + u_[IJPM]));
          const MFloat dvdxi = sign_xi * 0.5 * (V - 0.5 * (v[0] + v_[IJPM]));
          const MFloat T2 = m_gamma * p_[IJPM] / rho_[IJPM];
          const MFloat dTdxi = sign_xi * 0.5 * (t - 0.5 * (T[0] + T2));

          eflux[0][IJ_] = dudet * surfaceMetrics[ysd * 2 + ysd] + dvdet * surfaceMetrics[ysd * 2 + xsd];
          eflux[1][IJ_] = dudet * surfaceMetrics[ysd * 2 + xsd];
          eflux[2][IJ_] = dvdet * surfaceMetrics[ysd * 2 + ysd];

          fflux[0][IJ_] = dudxi * surfaceMetrics[xsd * 2 + ysd] + dvdxi * surfaceMetrics[xsd * 2 + xsd];
          fflux[1][IJ_] = dudxi * surfaceMetrics[xsd * 2 + xsd];
          fflux[2][IJ_] = dvdxi * surfaceMetrics[xsd * 2 + ysd];

          dT[0][IJ_] = dTdet * surfaceMetrics[ysd * 2 + xsd];
          dT[1][IJ_] = dTdet * surfaceMetrics[ysd * 2 + ysd];
          dT[2][IJ_] = dTdxi * surfaceMetrics[xsd * 2 + xsd];
          dT[3][IJ_] = dTdxi * surfaceMetrics[xsd * 2 + ysd];
        }
      }
    }
  }
}


void FvStructuredSolver2D::computePorousRHS(MBool isRans) {
  TRACE();

  static constexpr MFloat minMFloat =
      1e-16; // std::min(std::numeric_limits<MFloat>::min(), 1.0/std::numeric_limits<MFloat>::max());
  const MFloat rRe0 = 1.0 / m_Re0;

  const MFloat* const* const RESTRICT pvars = m_cells->pvariables;
  const MFloat* const RESTRICT muLam = &m_cells->fq[FQ->MU_L][0];
  const MFloat* const RESTRICT por = &m_cells->fq[FQ->POROSITY][0];
  const MFloat* const RESTRICT Da = &m_cells->fq[FQ->DARCY][0];
  const MFloat* const RESTRICT cf = &m_cells->fq[FQ->FORCH][0];
  const MFloat* const RESTRICT utau = &m_cells->fq[FQ->UTAU][0];
  const MFloat* const RESTRICT wallDist = &m_cells->fq[FQ->WALLDISTANCE][0];

  if(isRans) {
    const MFloat* const RESTRICT u_ = &pvars[PV->U][0];
    const MFloat* const RESTRICT v_ = &pvars[PV->V][0];
    const MFloat* const RESTRICT rho = &pvars[PV->RHO][0];
    const MFloat* const RESTRICT p = &pvars[PV->P][0];
    const MFloat* const RESTRICT TKE = &pvars[PV->RANS_VAR[0]][0];
    const MFloat* const RESTRICT eps = &pvars[PV->RANS_VAR[1]][0];
    const MFloat* const RESTRICT muTurb = &m_cells->fq[FQ->MU_T][0];

    MFloat* const* const RESTRICT eflux = ALIGNED_F(m_cells->eFlux);
    MFloat* const* const RESTRICT fflux = ALIGNED_F(m_cells->fFlux);
    MFloat* const* const RESTRICT vflux = ALIGNED_F(m_cells->viscousFlux);
    MFloat* const* const RESTRICT sa_1flux = ALIGNED_F(m_cells->saFlux1);
    MFloat* const* const RESTRICT sa_2flux = ALIGNED_F(m_cells->saFlux2);

    for(MInt j = m_noGhostLayers - 1; j < m_nCells[0] - m_noGhostLayers; ++j) {
      for(MInt i = m_noGhostLayers - 1; i < m_nCells[1] - m_noGhostLayers; ++i) {
        // get the adjacent cells;
        const MInt IJ = cellIndex(i, j);
        const MInt IPJ = cellIndex((i + 1), j);
        const MInt IPJP = cellIndex((i + 1), (j + 1));
        const MInt IJP = cellIndex(i, (j + 1));

        // Compute Reynolds stresses at the corners
        const MFloat dudxi = F1B2 * (u_[IPJP] + u_[IPJ] - u_[IJP] - u_[IJ]);
        const MFloat dudet = F1B2 * (u_[IPJP] + u_[IJP] - u_[IPJ] - u_[IJ]);

        const MFloat dvdxi = F1B2 * (v_[IPJP] + v_[IPJ] - v_[IJP] - v_[IJ]);
        const MFloat dvdet = F1B2 * (v_[IPJP] + v_[IJP] - v_[IPJ] - v_[IJ]);

        const MFloat rhoAvg = F1B4 * (rho[IJP] + rho[IPJP] + rho[IJ] + rho[IPJ]);
        const MFloat muTurbAvg = F1B4 * (muTurb[IJP] + muTurb[IPJP] + muTurb[IJ] + muTurb[IPJ]);
        const MFloat nuTurbAvg = muTurbAvg / rhoAvg;

        const MFloat TKEcorner =
            (m_rans2eq_mode == "production") ? -2 / 3 * F1B4 * (TKE[IJP] + TKE[IPJP] + TKE[IJ] + TKE[IPJ]) : 0.0;

        const MFloat cornerMetrics[4] = {m_cells->cornerMetrics[0][IJ], m_cells->cornerMetrics[1][IJ],
                                         m_cells->cornerMetrics[2][IJ], m_cells->cornerMetrics[3][IJ]};

        // compute tau1 = 2 du/dx - 2/3 ( du/dx + dv/dy + dw/dz )

        // tau_xx = 4/3*( du/dxi * dxi/dx + du/deta * deta/dx + du/dzeta * dzeta/dx )
        //            - 2/3*( dv/dxi * dxi/dy + dv/deta * deta/dy + dv/dzeta * dzeta/dy)
        //            - 2/3*( dw/dxi * dxi/dz + dw/deta * deta/dz + dw/dzeta * dzeta/dz )
        MFloat tau1 = F4B3 * (dudxi * cornerMetrics[xsd * 2 + xsd] + dudet * cornerMetrics[ysd * 2 + xsd]) -

                      F2B3 * (dvdxi * cornerMetrics[xsd * 2 + ysd] + dvdet * cornerMetrics[ysd * 2 + ysd]);

        // compute tau2 = du/dy + dv/dx

        // tau_xy = du/dxi * dxi/dy + du/deta * deta/dy + du/dzeta * dzeta/dy
        //        + dv/dxi * dxi/dx + dv/deta * deta/dx + dv/dzeta * dzeta/dx
        MFloat tau2 = dudxi * cornerMetrics[xsd * 2 + ysd] + dudet * cornerMetrics[ysd * 2 + ysd] +

                      dvdxi * cornerMetrics[xsd * 2 + xsd] + dvdet * cornerMetrics[ysd * 2 + xsd];

        // compute tau4 = 2 dv/dy - 2/3 ( du/dx + dv/dy + dw/dz )

        // tau_yy = 4/3*( dv/dxi * dxi/dy + dv/deta * deta/dy + dv/dzeta * dzeta/dy )
        //            - 2/3*( du/dxi * dxi/dx + du/deta * deta/dx + du/dzeta * dzeta/dx)
        //            - 2/3*( dw/dxi * dxi/dz + dw/deta * deta/dz + dw/dzeta * dzeta/dz )
        MFloat tau4 = F4B3 * (dvdxi * cornerMetrics[xsd * 2 + ysd] + dvdet * cornerMetrics[ysd * 2 + ysd]) -

                      F2B3 * (dudxi * cornerMetrics[xsd * 2 + xsd] + dudet * cornerMetrics[ysd * 2 + xsd]);

        const MFloat fJac = 1.0 / m_cells->cornerJac[IJ];
        const MFloat mueOverRe = rRe0 * fJac * nuTurbAvg;
        tau1 = mueOverRe * tau1 + TKEcorner;
        tau2 *= mueOverRe;
        tau4 = mueOverRe * tau4 + TKEcorner;

        vflux[0][IJ] = -tau1;
        vflux[1][IJ] = -tau2;
        vflux[2][IJ] = -tau4;

        // c_Dp-term
        const MFloat dTKEdxi = F1B2 * (TKE[IPJP] + TKE[IPJ] - TKE[IJP] - TKE[IJ]);
        const MFloat dTKEdet = F1B2 * (TKE[IPJP] + TKE[IJP] - TKE[IPJ] - TKE[IJ]);

        const MFloat depsdxi = F1B2 * (eps[IPJP] + eps[IPJ] - eps[IJP] - eps[IJ]);
        const MFloat depsdet = F1B2 * (eps[IPJP] + eps[IJP] - eps[IPJ] - eps[IJ]);

        //      const MFloat TKEAvg = F1B4*(TKE[IJP]+TKE[IPJP]+TKE[IJ]+TKE[IPJ]);
        //      const MFloat epsAvg = F1B4*(eps[IJP]+eps[IPJP]+eps[IJ]+eps[IPJ]);
        //      const MFloat temp = POW2(TKEAvg)/std::max(epsAvg, minMFloat);

        const MFloat dTKEdx = dTKEdxi * cornerMetrics[xsd * 2 + xsd] + dTKEdet * cornerMetrics[ysd * 2 + xsd];

        const MFloat dTKEdy = dTKEdxi * cornerMetrics[xsd * 2 + ysd] + dTKEdet * cornerMetrics[ysd * 2 + ysd];

        const MFloat depsdx = depsdxi * cornerMetrics[xsd * 2 + xsd] + depsdet * cornerMetrics[ysd * 2 + xsd];

        const MFloat depsdy = depsdxi * cornerMetrics[xsd * 2 + ysd] + depsdet * cornerMetrics[ysd * 2 + ysd];

        // Compute indicator function
        const MFloat dpdxi = F1B2 * (p[IPJP] + p[IPJ] - p[IJP] - p[IJ]);
        const MFloat dpdet = F1B2 * (p[IPJP] + p[IJP] - p[IPJ] - p[IJ]);

        const MFloat dpdx = dpdxi * cornerMetrics[xsd * 2 + xsd] + dpdet * cornerMetrics[ysd * 2 + xsd];

        const MFloat dpdy = dpdxi * cornerMetrics[xsd * 2 + ysd] + dpdet * cornerMetrics[ysd * 2 + ysd];

        const MFloat uAvg = F1B4 * (u_[IJP] + u_[IPJP] + u_[IJ] + u_[IPJ]);
        const MFloat vAvg = F1B4 * (v_[IJP] + v_[IPJP] + v_[IJ] + v_[IPJ]);
        const MFloat velAbs = sqrt(POW2(uAvg) + POW2(vAvg));
        const MFloat muLamAvg = F1B4 * (muLam[IJP] + muLam[IPJP] + muLam[IJ] + muLam[IPJ]);
        const MFloat porAvg = F1B4 * (por[IJP] + por[IPJP] + por[IJ] + por[IPJ]);
        // TODO_SS labels:FV,toenhance For now I take the Da and cf of the current cell and not at the corner;
        //       To do it correctly we need to exchange Da and cf similar to porosity
        const MFloat rDaAvg = 1.0 / Da[IJ];
        const MFloat cfAvg = cf[IJ];
        const MFloat indicator = (dpdx * uAvg + dpdy * vAvg)
                                 / std::max(POW2(velAbs)
                                                * (rRe0 * porAvg * muLamAvg * rDaAvg
                                                   + rhoAvg * POW2(porAvg) * cfAvg * sqrt(rDaAvg) * velAbs),
                                            minMFloat);
        const MFloat c_Dp_eff = m_c_Dp * (-0.5 * tanh(indicator - 5.0) + 0.5);
        const MFloat c_Dp_eps_eff = m_c_Dp_eps * (-0.5 * tanh(indicator - 5.0) + 0.5);
        m_cells->fq[FQ->POROUS_INDICATOR][IJ] =
            (-0.5 * tanh(indicator - 5.0) + 0.5); // it's not exactly the cell center indicator

        // TODO_SS labels:FV evtl. das f_mu wieder loeschen, wenn der Code auch ohne laueft
        const MFloat nuLaminarAvg = muLamAvg / rhoAvg;
        const MFloat utauAvg = F1B4 * (utau[IJP] + utau[IPJP] + utau[IJ] + utau[IPJ]);
        const MFloat wallDistAvg = F1B4 * (wallDist[IJP] + wallDist[IPJP] + wallDist[IJ] + wallDist[IPJ]);
        const MFloat yp = utauAvg * wallDistAvg / nuLaminarAvg;
        const MFloat f_mu = 1.0 - exp(-RM_KEPS::C3 * 40 * m_Re0 * yp); // TODO_SS labels:FV

        MFloat limiterVisc1 = 1.0;
        MFloat limiterVisc2 = 1.0;
        if(m_limiterVisc) {
          // TODO_SS labels:FV limiter should preserve conservativity; auch in computePorousRHSCorrection einbauen
          const MFloat mu_ref = F1B4
                                * (m_cells->fq[FQ->LIMITERVISC][IPJP] + m_cells->fq[FQ->LIMITERVISC][IPJ]
                                   + m_cells->fq[FQ->LIMITERVISC][IJP] + m_cells->fq[FQ->LIMITERVISC][IJ]);
          //        limiterVisc1 = std::min(1.0, mu_ref/std::abs( c_Dp_eff*temp ));
          //        limiterVisc2 = std::min(1.0, mu_ref/std::abs( c_Dp_eps_eff*temp ));
          limiterVisc1 = std::min(1.0,
                                  mu_ref
                                      / std::abs(rRe0 * c_Dp_eff * RM_KEPS::rsigma_k
                                                 * muTurbAvg)); // TODO_SS labels:FV scale with Re-number?
          limiterVisc2 = std::min(1.0, mu_ref / std::abs(rRe0 * c_Dp_eps_eff * RM_KEPS::rsigma_eps * muTurbAvg));
        }

        const MFloat Frj = rRe0 * fJac;

        const MFloat sax1 = Frj * c_Dp_eff * RM_KEPS::rsigma_k * muTurbAvg * f_mu * limiterVisc1
                            * (dTKEdx * cornerMetrics[xsd * 2 + xsd] + dTKEdy * cornerMetrics[xsd * 2 + ysd]);
        const MFloat say1 = Frj * c_Dp_eff * RM_KEPS::rsigma_k * muTurbAvg * f_mu * limiterVisc1
                            * (dTKEdx * cornerMetrics[ysd * 2 + xsd] + dTKEdy * cornerMetrics[ysd * 2 + ysd]);

        const MFloat sax2 = Frj * c_Dp_eps_eff * muTurbAvg * RM_KEPS::rsigma_eps * f_mu * limiterVisc2
                            * (depsdx * cornerMetrics[xsd * 2 + xsd] + depsdy * cornerMetrics[xsd * 2 + ysd]);
        const MFloat say2 = Frj * c_Dp_eps_eff * muTurbAvg * RM_KEPS::rsigma_eps * f_mu * limiterVisc2
                            * (depsdx * cornerMetrics[ysd * 2 + xsd] + depsdy * cornerMetrics[ysd * 2 + ysd]);

        //      const MFloat sax1 = fJac*c_Dp_eff*temp*
        //        (dTKEdx * cornerMetrics[ xsd * 2 + xsd ]+
        //         dTKEdy * cornerMetrics[ xsd * 2 + ysd ]);
        //      const MFloat say1 = fJac*c_Dp_eff*temp*
        //        (dTKEdx * cornerMetrics[ ysd * 2 + xsd ]+
        //         dTKEdy * cornerMetrics[ ysd * 2 + ysd ]);

        // TODO_SS labels:FV for now I assume that whenever rans-porous, then sa_1flux array is allocated
        sa_1flux[0][IJ] = sax1;
        sa_1flux[1][IJ] = say1;
        sa_2flux[0][IJ] = sax2;
        sa_2flux[1][IJ] = say2;
      }
    }

    if(m_hasSingularity) computePorousRHSCorrection();

    for(MInt j = m_noGhostLayers; j < m_nCells[0] - m_noGhostLayers; ++j) {
      for(MInt i = m_noGhostLayers - 1; i < m_nCells[1] - m_noGhostLayers; ++i) {
        const MInt IJ = cellIndex(i, j);
        const MInt IJM = cellIndex(i, (j - 1));

        sa_1flux[0][IJM] = F1B2 * (sa_1flux[0][IJ] + sa_1flux[0][IJM]);
        sa_2flux[0][IJM] = F1B2 * (sa_2flux[0][IJ] + sa_2flux[0][IJM]);
      }
    }
    for(MInt j = m_noGhostLayers - 1; j < m_nCells[0] - m_noGhostLayers; ++j) {
      for(MInt i = m_noGhostLayers; i < m_nCells[1] - m_noGhostLayers; ++i) {
        const MInt IJ = cellIndex(i, j);
        const MInt IMJ = cellIndex((i - 1), j);

        sa_1flux[1][IMJ] = F1B2 * (sa_1flux[1][IJ] + sa_1flux[1][IMJ]);
        sa_2flux[1][IMJ] = F1B2 * (sa_2flux[1][IJ] + sa_2flux[1][IMJ]);
      }
    }

    for(MInt var = 0; var < 3; ++var) {
      // intermediate saving of Reynolds stresses in the surface centers
      for(MInt j = m_noGhostLayers; j < m_nCells[0] - m_noGhostLayers; ++j) {
        for(MInt i = m_noGhostLayers - 1; i < m_nCells[1] - m_noGhostLayers; ++i) {
          const MInt IJ = cellIndex(i, j);
          const MInt IJM = cellIndex(i, j - 1);

          eflux[var][IJ] = 0.5 * (vflux[var][IJ] + vflux[var][IJM]);
        }
      }
      for(MInt j = m_noGhostLayers - 1; j < m_nCells[0] - m_noGhostLayers; ++j) {
        for(MInt i = m_noGhostLayers; i < m_nCells[1] - m_noGhostLayers; ++i) {
          const MInt IJ = cellIndex(i, j);
          const MInt IMJ = cellIndex((i - 1), j);

          fflux[var][IJ] = 0.5 * (vflux[var][IJ] + vflux[var][IMJ]);
        }
      }
    }

    // intermediate saving of velocity magnitude in the cell centers
    // TODO_SS labels:FV this intermediate could have done in the first for loop already
    for(MInt j = m_noGhostLayers - 1; j < m_nCells[0] - m_noGhostLayers + 1; ++j) {
      for(MInt i = m_noGhostLayers - 1; i < m_nCells[1] - m_noGhostLayers + 1; ++i) {
        const MInt IJ = cellIndex(i, j);

        vflux[0][IJ] = sqrt(POW2(u_[IJ]) + POW2(v_[IJ]));
        vflux[1][IJ] = 1.0 / std::max(vflux[0][IJ], minMFloat);
      }
    }

    for(MInt j = m_noGhostLayers; j < m_nCells[0] - m_noGhostLayers; ++j) {
      for(MInt i = m_noGhostLayers; i < m_nCells[1] - m_noGhostLayers; ++i) {
        const MInt IJ = cellIndex(i, j);
        const MInt IMJ = cellIndex(i - 1, j);
        const MInt IJM = cellIndex(i, j - 1);
        const MInt IPJ = cellIndex(i + 1, j);
        const MInt IJP = cellIndex(i, j + 1);

        const MFloat velAbs = vflux[0][IJ]; // sqrt(POW2(u_[IJ]) + POW2(v_[IJ]));
        const MFloat rDa = 1.0 / Da[IJ];
        const MFloat porPOW2 = POW2(por[IJ]);

        const MFloat invCellJac = 1.0 / m_cells->cellJac[IJ];

        // Average Reynolds stresses in the surface centers into cell center
        const MFloat tau1 = 0.25 * (eflux[0][IJ] + eflux[0][IMJ] + fflux[0][IJ] + fflux[0][IJM]);
        const MFloat tau2 = 0.25 * (eflux[1][IJ] + eflux[1][IMJ] + fflux[1][IJ] + fflux[1][IJM]);
        const MFloat tau4 = 0.25 * (eflux[2][IJ] + eflux[2][IMJ] + fflux[2][IJ] + fflux[2][IJM]);

        const MFloat dtau1dxi = eflux[0][IJ] - eflux[0][IMJ];
        const MFloat dtau2dxi = eflux[1][IJ] - eflux[1][IMJ];
        const MFloat dtau4dxi = eflux[2][IJ] - eflux[2][IMJ];
        const MFloat dtau1det = fflux[0][IJ] - fflux[0][IJM];
        const MFloat dtau2det = fflux[1][IJ] - fflux[1][IJM];
        const MFloat dtau4det = fflux[2][IJ] - fflux[2][IJM];

        const MFloat dtau1dx =
            (dtau1dxi * m_cells->cellMetrics[xsd * 2 + xsd][IJ] + dtau1det * m_cells->cellMetrics[ysd * 2 + xsd][IJ])
            * invCellJac;
        const MFloat dtau1dy =
            (dtau1dxi * m_cells->cellMetrics[xsd * 2 + ysd][IJ] + dtau1det * m_cells->cellMetrics[ysd * 2 + ysd][IJ])
            * invCellJac;
        const MFloat dtau2dx =
            (dtau2dxi * m_cells->cellMetrics[xsd * 2 + xsd][IJ] + dtau2det * m_cells->cellMetrics[ysd * 2 + xsd][IJ])
            * invCellJac;
        const MFloat dtau2dy =
            (dtau2dxi * m_cells->cellMetrics[xsd * 2 + ysd][IJ] + dtau2det * m_cells->cellMetrics[ysd * 2 + ysd][IJ])
            * invCellJac;
        const MFloat dtau4dx =
            (dtau4dxi * m_cells->cellMetrics[xsd * 2 + xsd][IJ] + dtau4det * m_cells->cellMetrics[ysd * 2 + xsd][IJ])
            * invCellJac;
        const MFloat dtau4dy =
            (dtau4dxi * m_cells->cellMetrics[xsd * 2 + ysd][IJ] + dtau4det * m_cells->cellMetrics[ysd * 2 + ysd][IJ])
            * invCellJac;

        // central discretization
        const MFloat dTKEdxi = 0.5 * (TKE[IPJ] - TKE[IMJ]);
        const MFloat dTKEdet = 0.5 * (TKE[IJP] - TKE[IJM]);
        const MFloat depsdxi = 0.5 * (eps[IPJ] - eps[IMJ]);
        const MFloat depsdet = 0.5 * (eps[IJP] - eps[IJM]);
        const MFloat dTKEdx =
            (dTKEdxi * m_cells->cellMetrics[xsd * 2 + xsd][IJ] + dTKEdet * m_cells->cellMetrics[ysd * 2 + xsd][IJ])
            * invCellJac;
        const MFloat dTKEdy =
            (dTKEdxi * m_cells->cellMetrics[xsd * 2 + ysd][IJ] + dTKEdet * m_cells->cellMetrics[ysd * 2 + ysd][IJ])
            * invCellJac;
        const MFloat depsdx =
            (depsdxi * m_cells->cellMetrics[xsd * 2 + xsd][IJ] + depsdet * m_cells->cellMetrics[ysd * 2 + xsd][IJ])
            * invCellJac;
        const MFloat depsdy =
            (depsdxi * m_cells->cellMetrics[xsd * 2 + ysd][IJ] + depsdet * m_cells->cellMetrics[ysd * 2 + ysd][IJ])
            * invCellJac;

        // central discretization
        const MFloat dvelAbsdxi = 0.5 * (vflux[0][IPJ] - vflux[0][IMJ]);
        const MFloat dvelAbsdet = 0.5 * (vflux[0][IJP] - vflux[0][IJM]);
        const MFloat dvelAbsdx = (dvelAbsdxi * m_cells->cellMetrics[xsd * 2 + xsd][IJ]
                                  + dvelAbsdet * m_cells->cellMetrics[ysd * 2 + xsd][IJ])
                                 * invCellJac;
        const MFloat dvelAbsdy = (dvelAbsdxi * m_cells->cellMetrics[xsd * 2 + ysd][IJ]
                                  + dvelAbsdet * m_cells->cellMetrics[ysd * 2 + ysd][IJ])
                                 * invCellJac;

        // central discretization
        const MFloat dvel11dxi = 0.5 * (u_[IPJ] * u_[IPJ] * vflux[1][IPJ] - u_[IMJ] * u_[IMJ] * vflux[1][IMJ]);
        const MFloat dvel11det = 0.5 * (u_[IJP] * u_[IJP] * vflux[1][IJP] - u_[IJM] * u_[IJM] * vflux[1][IJM]);
        const MFloat dvel12dxi = 0.5 * (u_[IPJ] * v_[IPJ] * vflux[1][IPJ] - u_[IMJ] * v_[IMJ] * vflux[1][IMJ]);
        const MFloat dvel12det = 0.5 * (u_[IJP] * v_[IJP] * vflux[1][IJP] - u_[IJM] * v_[IJM] * vflux[1][IJM]);
        const MFloat dvel22dxi = 0.5 * (v_[IPJ] * v_[IPJ] * vflux[1][IPJ] - v_[IMJ] * v_[IMJ] * vflux[1][IMJ]);
        const MFloat dvel22det = 0.5 * (v_[IJP] * v_[IJP] * vflux[1][IJP] - v_[IJM] * v_[IJM] * vflux[1][IJM]);
        const MFloat dvel11dx =
            (dvel11dxi * m_cells->cellMetrics[xsd * 2 + xsd][IJ] + dvel11det * m_cells->cellMetrics[ysd * 2 + xsd][IJ])
            * invCellJac;
        const MFloat dvel11dy =
            (dvel11dxi * m_cells->cellMetrics[xsd * 2 + ysd][IJ] + dvel11det * m_cells->cellMetrics[ysd * 2 + ysd][IJ])
            * invCellJac;
        const MFloat dvel12dx =
            (dvel12dxi * m_cells->cellMetrics[xsd * 2 + xsd][IJ] + dvel12det * m_cells->cellMetrics[ysd * 2 + xsd][IJ])
            * invCellJac;
        const MFloat dvel12dy =
            (dvel12dxi * m_cells->cellMetrics[xsd * 2 + ysd][IJ] + dvel12det * m_cells->cellMetrics[ysd * 2 + ysd][IJ])
            * invCellJac;
        const MFloat dvel22dx =
            (dvel22dxi * m_cells->cellMetrics[xsd * 2 + xsd][IJ] + dvel22det * m_cells->cellMetrics[ysd * 2 + xsd][IJ])
            * invCellJac;
        const MFloat dvel22dy =
            (dvel22dxi * m_cells->cellMetrics[xsd * 2 + ysd][IJ] + dvel22det * m_cells->cellMetrics[ysd * 2 + ysd][IJ])
            * invCellJac;

        const MFloat u = u_[IJ];
        const MFloat v = v_[IJ];

        const MFloat temp = u * u * tau1 + 2.0 * u * v * tau2 + v * v * tau4;
        const MFloat kovereps = TKE[IJ] / std::max(eps[IJ], minMFloat);
        const MFloat velAbsInv = vflux[1][IJ];
        const MFloat tau111 = 3.0 * (tau1 * dtau1dx + tau2 * dtau1dy);
        const MFloat tau112 = 2.0 * (tau1 * dtau2dx + tau2 * dtau2dy) + tau2 * dtau1dx + tau4 * dtau1dy;
        const MFloat tau122 = tau1 * dtau4dx + tau2 * dtau4dy + 2 * (tau2 * dtau2dx + tau4 * dtau2dy);
        const MFloat tau222 = 3.0 * (tau2 * dtau4dx + tau4 * dtau4dy);

        // c_Dp-term
        const MInt IMJM = cellIndex(i - 1, j - 1);
        /*        MFloat limiterVisc = 1.0;
      if (m_limiterVisc) {
                  //TODO_SS labels:FV limiter should preserve conservativity
                  const MFloat mu_ref = m_cells->fq[FQ->LIMITERVISC][IJ];
                  limiterVisc = std::min(1.0, mu_ref/std::abs( m_c_Dp*TKE[IJ]*kovereps ));
      }*/
        const MFloat diffusion_TKE =
            ((sa_1flux[0][IJM] - sa_1flux[0][IMJM]) + (sa_1flux[1][IMJ] - sa_1flux[1][IMJM])); // * limiterVisc;
        const MFloat diffusion_eps =
            ((sa_2flux[0][IJM] - sa_2flux[0][IMJM]) + (sa_2flux[1][IMJ] - sa_2flux[1][IMJM])); // * limiterVisc;

        /*        m_cells->fq[FQ->VAR5][IJ] = diffusion_TKE;
      m_cells->fq[FQ->VAR6][IJ] = (- 2.0*rRe0*rDa*por[IJ]*muLam[IJ]*TKE[IJ]
      - 0.5*sqrt(rDa)*porPOW2*cf[IJ]*rho[IJ]*(
      4.0*TKE[IJ]*velAbs
      + 2.0*velAbsInv*temp))*m_cells->cellJac[IJ];*/

        // RHS
        m_cells->rightHandSide[CV->RHO_U][IJ] +=
            (-rRe0 * rDa * por[IJ] * muLam[IJ] * u
             - sqrt(rDa) * porPOW2 * cf[IJ] * rho[IJ]
                   * ((velAbs + TKE[IJ] * velAbsInv - 0.5 * pow(velAbsInv, 3) * temp) * u
                      + velAbsInv * (u * tau1 + v * tau2)))
            * m_cells->cellJac[IJ];
        m_cells->rightHandSide[CV->RHO_V][IJ] +=
            (-rRe0 * rDa * por[IJ] * muLam[IJ] * v
             - sqrt(rDa) * porPOW2 * cf[IJ] * rho[IJ]
                   * ((velAbs + TKE[IJ] * velAbsInv - 0.5 * pow(velAbsInv, 3) * temp) * v
                      + velAbsInv * (u * tau2 + v * tau4)))
            * m_cells->cellJac[IJ];
        m_cells->rightHandSide[CV->RANS_VAR[0]][IJ] +=
            diffusion_TKE
            + (-2.0 * rRe0 * rDa * por[IJ] * muLam[IJ] * TKE[IJ]
               - 0.5 * sqrt(rDa) * porPOW2 * cf[IJ] * rho[IJ]
                     * (4.0 * TKE[IJ] * velAbs + 2.0 * velAbsInv * temp
                        + 3.0 * velAbsInv * m_c_t * kovereps * (u * tau111 + u * tau122 + v * tau112 + v * tau222)
                        - pow(velAbsInv, 3) * m_c_t * kovereps
                              * (u * u * u * tau111 + 3.0 * u * u * v * tau112 + 3.0 * u * v * v * tau122
                                 + v * v * v * tau222)))
                  * m_cells->cellJac[IJ];
        m_cells->rightHandSide[CV->RANS_VAR[1]][IJ] +=
            diffusion_eps
            + (-2.0 * rRe0 * rDa * por[IJ] * muLam[IJ] * eps[IJ]
               - sqrt(rDa) * porPOW2 * cf[IJ] * rho[IJ]
                     * (F8B3 * velAbs * eps[IJ]
                        + 2.0 * rRe0 * (muLam[IJ] / rho[IJ]) * (dvelAbsdx * dTKEdx + dvelAbsdy * dTKEdy)
                        - 2.0 * m_c_eps * kovereps * velAbsInv
                              * (u * tau1 * depsdx + v * tau2 * depsdx + u * tau2 * depsdy + v * tau4 * depsdy)
                        + rRe0 * (muLam[IJ] / rho[IJ])
                              * (dvel11dx * dtau1dx + dvel11dy * dtau1dy + 2.0 * dvel12dx * dtau2dx
                                 + 2.0 * dvel12dy * dtau2dy + dvel22dx * dtau4dx + dvel22dy * dtau4dy)))
                  * m_cells->cellJac[IJ];

        //        m_cells->fq[FQ->VAR3][IJ] = rho[IJ]*diffusion_TKE;
        //        m_cells->fq[FQ->VAR4][IJ] = - 2.0*rRe0*rDa*por[IJ]*muLam[IJ]*TKE[IJ]*m_cells->cellJac[IJ];
        //        m_cells->fq[FQ->VAR5][IJ] = -
        //        0.5*sqrt(rDa)*porPOW2*cf[IJ]*rho[IJ]*(4.0*TKE[IJ]*velAbs)*m_cells->cellJac[IJ];
        //        m_cells->fq[FQ->VAR6][IJ] = -
        //        0.5*sqrt(rDa)*porPOW2*cf[IJ]*rho[IJ]*(2.0*velAbsInv*temp)*m_cells->cellJac[IJ];

        //        const MFloat dx0 = m_cells->coordinates[1][IJ] - m_cells->coordinates[1][IJM];
        //        const MFloat dx1 = m_cells->coordinates[1][IJP] - m_cells->coordinates[1][IJ];
        //        const MFloat ddTKE = 2.0*(TKE[IJP]*dx0 - TKE[IJ]*(dx0+dx1) +
        //        TKE[IJM]*dx1)/(POW2(dx0)*dx1+POW2(dx1)*dx0); m_cells->fq[FQ->VAR7][IJ] =
        //        m_c_Dp*POW2(TKE[IJ])/eps[IJ]*rho[IJ]*ddTKE;///(2.0*rRe0*rDa*por[IJ]*muLam[IJ]*TKE[IJ]);
        /*        std::ignore = (rho[IJ]*diffusion_TKE
      - 2.0*rRe0*rDa*por[IJ]*muLam[IJ]*TKE[IJ]
      - 0.5*sqrt(rDa)*porPOW2*cf[IJ]*rho[IJ]*(
      4.0*TKE[IJ]*velAbs
      + 2.0*velAbsInv*temp
      + 3.0*velAbsInv*m_c_t*kovereps*(u*tau111+u*tau122+v*tau112+v*tau222)
      -
      pow(velAbsInv,3)*m_c_t*kovereps*(u*u*u*tau111+3.0*u*u*v*tau112+3.0*u*v*v*tau122+v*v*v*tau222)))*m_cells->cellJac[IJ];
      std::ignore = (-2.0*rRe0*rDa*por[IJ]*muLam[IJ]*eps[IJ]
      - sqrt(rDa)*porPOW2*cf[IJ]*rho[IJ]*(
      F8B3*velAbs*eps[IJ]
      + 2.0*rRe0*(muLam[IJ]/rho[IJ])*(dvelAbsdx*dTKEdx+dvelAbsdy*dTKEdy)
      - 2.0*m_c_eps*kovereps*velAbsInv*(u*tau1*depsdx+v*tau2*depsdx+u*tau2*depsdy+v*tau4*depsdy)
      +
      rRe0*(muLam[IJ]/rho[IJ])*(dvel11dx*dtau1dx+dvel11dy*dtau1dy+2.0*dvel12dx*dtau2dx+2.0*dvel12dy*dtau2dy+dvel22dx*dtau4dx+dvel22dy*dtau4dy)))*m_cells->cellJac[IJ];*/
      }
    }
  } else { // LES
    for(MInt j = m_noGhostLayers; j < m_nCells[0] - m_noGhostLayers; ++j) {
      for(MInt i = m_noGhostLayers; i < m_nCells[1] - m_noGhostLayers; ++i) {
        const MInt IJK = cellIndex(i, j);
        MFloat velAbs = F0;
        for(MInt dim = 0; dim < nDim; ++dim) {
          velAbs += POW2(pvars[PV->VV[dim]][IJK]);
        }
        velAbs = sqrt(velAbs);
        const MFloat rDa = 1.0 / Da[IJK];
        const MFloat porPOW2 = POW2(por[IJK]);
        for(MInt dim = 0; dim < nDim; ++dim) {
          m_cells->rightHandSide[CV->RHO_VV[dim]][IJK] +=
              -(rRe0 * rDa * por[IJK] * muLam[IJK] + porPOW2 * sqrt(rDa) * cf[IJK] * velAbs * pvars[PV->RHO][IJK])
              * pvars[PV->VV[dim]][IJK] * m_cells->cellJac[IJK];
        }
      }
    }
  }
}


void FvStructuredSolver2D::computePorousRHSCorrection() {
  static constexpr MFloat minMFloat =
      1e-16; // std::min(std::numeric_limits<MFloat>::min(), 1.0/std::numeric_limits<MFloat>::max());
  const MFloat rRe0 = 1.0 / m_Re0;

  MFloat* RESTRICT u_ = &m_cells->pvariables[PV->U][0];
  MFloat* RESTRICT v_ = &m_cells->pvariables[PV->V][0];
  MFloat* RESTRICT p_ = &m_cells->pvariables[PV->P][0];
  MFloat* RESTRICT rho_ = &m_cells->pvariables[PV->RHO][0];
  const MFloat* const RESTRICT muLam_ = &m_cells->fq[FQ->MU_L][0];
  MFloat* RESTRICT muTurb_ = &m_cells->fq[FQ->MU_T][0];
  MFloat* const RESTRICT TKE_ = &m_cells->pvariables[PV->RANS_VAR[0]][0];
  MFloat* const RESTRICT EPS_ = &m_cells->pvariables[PV->RANS_VAR[1]][0];
  const MFloat* const RESTRICT por_ = &m_cells->fq[FQ->POROSITY][0];
  const MFloat* const RESTRICT Da = &m_cells->fq[FQ->DARCY][0];
  const MFloat* const RESTRICT cf = &m_cells->fq[FQ->FORCH][0];

  //  MFloat *const RESTRICT eflux= ALIGNED_MF(m_cells->eFlux);
  //  MFloat *const RESTRICT fflux= ALIGNED_MF(m_cells->fFlux);
  MFloat* const* const RESTRICT vflux = ALIGNED_MF(m_cells->viscousFlux);
  MFloat* const* const RESTRICT sa_1flux = ALIGNED_F(m_cells->saFlux1);
  MFloat* const* const RESTRICT sa_2flux = ALIGNED_F(m_cells->saFlux2);

  MInt start[nDim], end[nDim], nghbr[6];
  MInt len1[nDim];

  for(MInt i = 0; i < m_hasSingularity; ++i) {
    // only correct for bc 6000 not for bc 4000-5000
    if(m_singularity[i].BC == -6000) {
      //      totalCells=1;
      for(MInt j = 0; j < nDim; j++) {
        len1[j] = m_singularity[i].end[j] - m_singularity[i].start[j];
        //        if(len1[j]!=0)  totalCells*=len1[j];
      }

      for(MInt n = 0; n < nDim; ++n) {
        if(m_singularity[i].end[n] - m_singularity[i].start[n] > 1) {
          mTerm(1, "In 2D not possible!");
          //          dim=n;
          // start[n]=m_singularity[i].start[n]+1;
          start[n] = m_singularity[i].start[n] + 1;
          end[n] = m_singularity[i].end[n] - 1;
        } else {
          start[n] = m_singularity[i].start[n];
          end[n] = m_singularity[i].end[n];
        }
      }

      MFloat u[6], v[6], p[6], rho[6], muTurb[6], muLam[6], TKE[6], EPS[6], POR[6];
      MFloat U, V, RHO, muTurb_corner, muLam_corner, TKEcorner, EPScorner, PORcorner, dudx, dudy, dvdx, dvdy, dTKEdx,
          dTKEdy, depsdx, depsdy, dpdx, dpdy;
      MFloat tau1, tau2, tau4;
      MFloat mueOverRe;
      for(MInt jj = start[1]; jj < end[1]; ++jj) {
        for(MInt ii = start[0]; ii < end[0]; ++ii) {
          MInt count = 0;
          //            MInt temp[nDim]{};
          MInt IJ_ = cellIndex(ii + m_singularity[i].Viscous[0], jj + m_singularity[i].Viscous[1]);

          //            temp[dim]=1;
          const MInt IJ = cellIndex(ii, jj);
          nghbr[count++] = IJ;
          //            nghbr[count++]=cellIndex(ii+temp[0],jj+temp[1],kk+temp[2]);

          const MInt IPMJ_ = getCellIdFromCell(IJ, m_singularity[i].Viscous[0], 0);
          const MInt IJPM_ = getCellIdFromCell(IJ, 0, m_singularity[i].Viscous[1]);

          const MFloat surfaceMetrics[nDim * nDim] = {
              m_cells->surfaceMetrics[0][IPMJ_], m_cells->surfaceMetrics[1][IPMJ_], m_cells->surfaceMetrics[2][IJPM_],
              m_cells->surfaceMetrics[3][IJPM_]};

          for(MInt m = 0; m < m_singularity[i].Nstar - 1; ++m) {
            MInt* change = m_singularity[i].displacement[m];
            nghbr[count++] = cellIndex(ii + change[0], jj + change[1]);
            //              nghbr[count++]=cellIndex(ii+temp[0]+change[0],jj+temp[1]+change[1],kk+temp[2]+change[2]);
          }

          if(count != m_singularity[i].Nstar) {
            cout << "what the hell! it is wrong!!! count=" << count << " Nstar=" << m_singularity[i].Nstar << endl;
          }

          for(MInt m = 0; m < m_singularity[i].Nstar; ++m) {
            u[m] = u_[nghbr[m]];
            v[m] = v_[nghbr[m]];
            muTurb[m] = muTurb_[nghbr[m]] / rho_[nghbr[m]]; // nuTurb
            muLam[m] = muLam_[nghbr[m]];
            TKE[m] = TKE_[nghbr[m]];
            EPS[m] = EPS_[nghbr[m]];
            p[m] = p_[nghbr[m]];
            rho[m] = rho_[nghbr[m]];
            POR[m] = por_[nghbr[m]];
          }

          U = F0;
          V = F0;
          RHO = F0;
          muTurb_corner = F0;
          muLam_corner = F0;
          TKEcorner = F0, EPScorner = F0, PORcorner = F0;
          dudx = F0;
          dudy = F0;
          dvdx = F0;
          dvdy = F0;
          dTKEdx = F0;
          dTKEdy = F0;
          depsdx = F0;
          depsdy = F0;
          dpdx = F0;
          dpdy = F0;

          MInt id2 = ii - start[0] + (jj - start[1]) * len1[0];

          for(MInt n = 0; n < count; n++) {
            MInt ID = id2 * count + n;
            U += m_singularity[i].ReconstructionConstants[0][ID] * u[n];
            dudx += m_singularity[i].ReconstructionConstants[1][ID] * u[n];
            dudy += m_singularity[i].ReconstructionConstants[2][ID] * u[n];

            V += m_singularity[i].ReconstructionConstants[0][ID] * v[n];
            dvdx += m_singularity[i].ReconstructionConstants[1][ID] * v[n];
            dvdy += m_singularity[i].ReconstructionConstants[2][ID] * v[n];

            TKEcorner += m_singularity[i].ReconstructionConstants[0][ID] * TKE[n];
            dTKEdx += m_singularity[i].ReconstructionConstants[1][ID] * TKE[n];
            dTKEdy += m_singularity[i].ReconstructionConstants[2][ID] * TKE[n];

            dpdx += m_singularity[i].ReconstructionConstants[1][ID] * p[n];
            dpdy += m_singularity[i].ReconstructionConstants[2][ID] * p[n];

            RHO += m_singularity[i].ReconstructionConstants[0][ID] * rho[n];

            EPScorner += m_singularity[i].ReconstructionConstants[0][ID] * EPS[n];
            depsdx += m_singularity[i].ReconstructionConstants[1][ID] * EPS[n];
            depsdy += m_singularity[i].ReconstructionConstants[2][ID] * EPS[n];

            muTurb_corner += m_singularity[i].ReconstructionConstants[0][ID] * muTurb[n];
            muLam_corner += m_singularity[i].ReconstructionConstants[0][ID] * muLam[n];
            PORcorner += m_singularity[i].ReconstructionConstants[0][ID] * POR[n];
          }

          //            cout << globalTimeStep << "(" << m_RKStep << ") dom=" << domainId() << " x|y=" <<
          //            m_cells->coordinates[0][IJK] << "|" << m_cells->coordinates[1][IJK] << " U=" << U << " V=" << V
          //            << " t=" << t << " cornerMetrics=" << cornerMetrics[0] << "|" << cornerMetrics[1] << "|" <<
          //            cornerMetrics[2] << "|" << cornerMetrics[3] << endl;


          tau1 = 2 * dudx - 2 / 3 * (dudx + dvdy);
          tau2 = dudy + dvdx;
          tau4 = 2 * dvdy - 2 / 3 * (dudx + dvdy);

          TKEcorner *= -2 / 3;
          mueOverRe = muTurb_corner * rRe0;
          tau1 = mueOverRe * tau1 + TKEcorner;
          tau2 *= mueOverRe;
          tau4 = mueOverRe * tau4 + TKEcorner;

          vflux[0][IJ_] = -tau1;
          vflux[1][IJ_] = -tau2;
          vflux[2][IJ_] = -tau4;


          const MFloat temp = POW2(TKEcorner) / std::max(EPScorner, minMFloat);
          const MFloat velAbs = sqrt(POW2(U) + POW2(V));
          // TODO_SS labels:FV,toenhance For now I take the Da and cf of the current cell and not at the corner;
          //       To do it correctly we need to exchange Da and cf similar to porosity
          const MFloat rDaAvg = 1.0 / Da[IJ];
          const MFloat cfAvg = cf[IJ];
          const MFloat indicator = (dpdx * U + dpdy * V)
                                   / std::max(POW2(velAbs)
                                                  * (rRe0 * PORcorner * muLam_corner * rDaAvg
                                                     + RHO * POW2(PORcorner) * cfAvg * sqrt(rDaAvg) * velAbs),
                                              minMFloat);
          const MFloat c_Dp_eff = m_c_Dp * (-0.5 * tanh(indicator - 5.0) + 0.5);
          const MFloat c_Dp_eps_eff = m_c_Dp_eps * (-0.5 * tanh(indicator - 5.0) + 0.5);
          m_cells->fq[FQ->POROUS_INDICATOR][IJ_] =
              (-0.5 * tanh(indicator - 5.0) + 0.5); // it's not exactly the cell center indicator

          // TODO_SS labels:FV f_mu & limiterVisc !!!!!!!!
          const MFloat sax1 = c_Dp_eff * RM_KEPS::rsigma_k * temp
                              * (dTKEdx * surfaceMetrics[xsd * 2 + xsd] + dTKEdy * surfaceMetrics[xsd * 2 + ysd]);
          const MFloat say1 = c_Dp_eff * RM_KEPS::rsigma_k * temp
                              * (dTKEdx * surfaceMetrics[ysd * 2 + xsd] + dTKEdy * surfaceMetrics[ysd * 2 + ysd]);
          const MFloat sax2 = c_Dp_eps_eff * RM_KEPS::rsigma_eps * temp
                              * (depsdx * surfaceMetrics[xsd * 2 + xsd] + depsdy * surfaceMetrics[xsd * 2 + ysd]);
          const MFloat say2 = c_Dp_eps_eff * RM_KEPS::rsigma_eps * temp
                              * (depsdx * surfaceMetrics[ysd * 2 + xsd] + depsdy * surfaceMetrics[ysd * 2 + ysd]);

          // TODO_SS labels:FV for now I assume that whenever rans-porous, then sa_1flux array is allocated
          sa_1flux[0][IJ_] = sax1;
          sa_1flux[1][IJ_] = say1;
          sa_2flux[0][IJ_] = sax2;
          sa_2flux[1][IJ_] = say2;
        }
      }
    }
  }
}


void FvStructuredSolver2D::loadRestartBC2600() {
  if(m_bc2600IsActive && !m_bc2600InitialStartup) {
    if(domainId() == 0) {
      cout << "Loading BC2600 values..." << endl;
    }

    ParallelIo::size_type bcCells[2] = {m_grid->getMyBlockNoCells(0), m_noGhostLayers};
    MInt noCellsBC = bcCells[0] * bcCells[1];
    ParallelIo::size_type bcOffset[2] = {0, 0};
    MFloatScratchSpace tmpRestartVars(noCellsBC * CV->noVariables, AT_, "m_tmpRestartVars2600");

    if(domainId() == 0) {
      stringstream restartFileName;
      MString restartFile = Context::getSolverProperty<MString>("restartVariablesFileName", m_solverId, AT_);
      restartFileName << outputDir() << restartFile;

      ParallelIoHdf5 pio(restartFileName.str(), maia::parallel_io::PIO_READ, MPI_COMM_SELF);
      stringstream pathStr;
      pathStr << "/block" << m_blockId << "/bc2600" << endl;
      std::string path = pathStr.str();

      for(MInt var = 0; var < CV->noVariables; var++) {
        cout << "Loading " << m_variableNames[var] << " offset: " << var * noCellsBC << endl;
        pio.readArray(&tmpRestartVars[var * noCellsBC], path, m_variableNames[var], nDim, bcOffset, bcCells);
      }

    }

    MPI_Bcast(&tmpRestartVars[0], noCellsBC * CV->noVariables, MPI_DOUBLE, 0, m_StructuredComm, AT_,
              "tmpRestartVars[0]");

    if(domainId() == 0) {
      cout << "Loading BC2600 values... SUCCESSFUL!" << endl;
    }

    if(m_bc2600) {
      MInt startGC[2] = {0, 0};
      MInt endGC[2] = {0, 0};

      if(m_bc2600noOffsetCells[0] == 0) {
        startGC[0] = m_noGhostLayers;
      }
      if(m_bc2600noOffsetCells[0] + m_bc2600noActiveCells[0] == bcCells[0]) {
        endGC[0] = m_noGhostLayers;
      }

      MInt startI = 0, endI = 0;
      if(m_bc2600Face == 0) {
        startI = 0;
        endI = m_noGhostLayers;
      } else if(m_bc2600Face == 1) {
        startI = m_nCells[2] - m_noGhostLayers - 1;
        endI = m_nCells[2];
      }

      for(MInt i = startI; i < endI; i++) {
        for(MInt j = startGC[0]; j < m_bc2600noCells[0] - endGC[0]; j++) {
          MInt cellId = cellIndex(i, j);
          MInt globalI = i;
          MInt globalJ = m_bc2600noOffsetCells[0] - m_noGhostLayers + j;
          MInt cellIdBC = globalI + globalJ * bcCells[1];

          // load values from restart field
          for(MInt var = 0; var < CV->noVariables; var++) {
            m_cells->variables[var][cellId] = tmpRestartVars[var * noCellsBC + cellIdBC];
          }
        }
      }


      // Fix diagonal cells at end of domain
      if(m_bc2600noOffsetCells[0] + m_bc2600noActiveCells[0] == m_grid->getMyBlockNoCells(0)) {
        for(MInt i = startI; i < endI; i++) {
          const MInt cellIdA2 = cellIndex(i, m_noGhostLayers + m_bc2600noActiveCells[0] - 2);
          const MInt cellIdA1 = cellIndex(i, m_noGhostLayers + m_bc2600noActiveCells[0] - 1);
          const MInt cellIdG1 = cellIndex(i, m_noGhostLayers + m_bc2600noActiveCells[0]);
          for(MInt var = 0; var < CV->noVariables; var++) {
            const MFloat distA1A2 = sqrt(POW2(m_cells->coordinates[0][cellIdA1] - m_cells->coordinates[0][cellIdA2])
                                         + POW2(m_cells->coordinates[2][cellIdA1] - m_cells->coordinates[1][cellIdA2]));
            const MFloat slope = (m_cells->variables[var][cellIdA1] - m_cells->variables[var][cellIdA2]) / distA1A2;
            const MFloat distG1A1 = sqrt(POW2(m_cells->coordinates[0][cellIdG1] - m_cells->coordinates[0][cellIdA1])
                                         + POW2(m_cells->coordinates[2][cellIdG1] - m_cells->coordinates[1][cellIdA1]));
            m_cells->variables[var][cellIdG1] = m_cells->variables[var][cellIdA1] + distG1A1 * slope;
          }
        }
      }
    }
  }
}

template <MFloat (FvStructuredSolver<2>::*pressure_func)(MInt) const>
void FvStructuredSolver2D::computePrimitiveVariables_() {
  const MFloat gammaMinusOne = m_gamma - 1.0;

  MFloat** const RESTRICT cvars = m_cells->variables;
  MFloat** const RESTRICT pvars = m_cells->pvariables;

  for(MInt j = m_noGhostLayers; j < m_nCells[0] - m_noGhostLayers; ++j) {
    for(MInt i = m_noGhostLayers; i < m_nCells[1] - m_noGhostLayers; ++i) {
      const MInt cellId = cellIndex(i, j);
      const MFloat fRho = F1 / cvars[CV->RHO][cellId];
      MFloat velPOW2 = F0;
      for(MInt vel = 0; vel < nDim; ++vel) { // compute velocity
        pvars[vel][cellId] = cvars[vel][cellId] * fRho;
        velPOW2 += POW2(pvars[vel][cellId]);
      }

      // density and pressure:
      pvars[PV->RHO][cellId] = cvars[CV->RHO][cellId]; // density
      pvars[PV->P][cellId] =
          gammaMinusOne
          * (cvars[CV->RHO_E][cellId] - F1B2 * pvars[PV->RHO][cellId] * velPOW2 + (this->*pressure_func)(cellId));

      for(MInt ransVar = 0; ransVar < m_noRansEquations; ransVar++) {
        cvars[CV->RANS_VAR[ransVar]][cellId] =
            mMax(cvars[CV->RANS_VAR[ransVar]][cellId], 1e-16 /*F0*/); // TODO_SS labels:FV maybe activate this
        pvars[PV->RANS_VAR[ransVar]][cellId] = cvars[CV->RANS_VAR[ransVar]][cellId] * fRho;
      }
    }
  }
}

void FvStructuredSolver2D::computePrimitiveVariables() {
  if(noRansEquations(m_ransMethod) == 2) {
    if(m_rans2eq_mode == "production")
      computePrimitiveVariables_<&FvStructuredSolver::pressure_twoEqRans>();
    else
      computePrimitiveVariables_();
  } else
    computePrimitiveVariables_();
}


void FvStructuredSolver2D::allocateSingularities() {
  for(MInt i = 0; i < m_hasSingularity; ++i) {
    MInt len[nDim];
    m_singularity[i].totalPoints = 1;
    m_singularity[i].totalCells = 1;

    for(MInt j = 0; j < nDim; j++) {
      len[j] = m_singularity[i].end[j] - m_singularity[i].start[j];
      m_singularity[i].totalPoints *= (len[j] + 1);
      m_singularity[i].totalCells *= len[j];
    }

    // 4 unknowns and Nstar cells
    mAlloc(m_singularity[i].ReconstructionConstants, nDim + 1, m_singularity[i].totalCells * m_singularity[i].Nstar,
           "ReconstructionConstants", 0.0, AT_);
  }
}


void FvStructuredSolver2D::gather(const MBool periodicExchange,
                                  std::vector<std::unique_ptr<StructuredComm<nDim>>>& sndComm) {
  for(auto& snd : sndComm) {
    if(isPeriodicComm(snd) && !periodicExchange) continue;
    if(periodicExchange && skipPeriodicDirection(snd)) continue;

    MInt pos = 0;
    for(MInt var = 0; var < snd->noVars; var++) {
      for(MInt j = snd->startInfoCells[1]; j < snd->endInfoCells[1]; j++) {
        for(MInt i = snd->startInfoCells[0]; i < snd->endInfoCells[0]; i++) {
          const MInt cellId = i + j * m_nCells[1];
          snd->cellBuffer[pos] = snd->variables[var][cellId];
          pos++;
        }
      }
    }
  }
}


void FvStructuredSolver2D::scatter(const MBool periodicExchange,
                                   std::vector<std::unique_ptr<StructuredComm<nDim>>>& rcvComm) {
  // the ordering of the grid points can be different from
  // sending instance ==> reorder it and copy it to the
  // right place

  for(auto& rcv : rcvComm) {
    if(isPeriodicComm(rcv) && !periodicExchange) continue;
    if(periodicExchange && skipPeriodicDirection(rcv)) continue;

    MInt step2[2];
    MInt start1[2];
    MInt start2[2];
    MInt end2[2];
    MInt len2[2];
    MInt totalCells = 1;
    MInt len1[2];

    for(MInt j = 0; j < nDim; j++) {
      len1[j] = rcv->endInfoCells[j] - rcv->startInfoCells[j];
      if(len1[j] != 0) totalCells *= len1[j];
      step2[rcv->orderInfo[j]] = rcv->stepInfo[j];
    }

    for(MInt j = 0; j < nDim; j++) {
      start2[j] = 0;
      end2[j] = len1[j] - 1;
      len2[rcv->orderInfo[j]] = len1[j];
      if(step2[j] < 0) {
        MInt dummy = start2[j];
        start2[j] = end2[j];
        end2[j] = dummy;
      }
    }

    MInt pos = 0;
    for(MInt var = 0; var < rcv->noVars; var++) {
      MInt j2 = start2[1];
      for(MInt j = rcv->startInfoCells[1]; j < rcv->endInfoCells[1]; j++) {
        MInt i2 = start2[0];
        for(MInt i = rcv->startInfoCells[0]; i < rcv->endInfoCells[0]; i++) {
          start1[rcv->orderInfo[0]] = i2;
          start1[rcv->orderInfo[1]] = j2;

          const MInt id2 = var * totalCells + start1[0] + (start1[1]) * len2[0];
          const MInt cellId = i + j * m_nCells[1];
          rcv->variables[var][cellId] = rcv->cellBuffer[id2];

          i2 += step2[0];
          pos++;
        }
        j2 += step2[1];
      }
    }
  }
}


void FvStructuredSolver2D::computeReconstructionConstantsSVD() {
  MInt nghbr[15 /*30*/]; //,dim = 0;
  MInt start[nDim], end[nDim];
  m_orderOfReconstruction = 1;
  const MInt recDim = (m_orderOfReconstruction == 2) ? (IPOW2(nDim) + 1) : nDim + 1;
  MInt maxNoSingularityRecNghbrIds = 7; // 14;
  MFloatScratchSpace tmpA(maxNoSingularityRecNghbrIds, recDim, AT_, "tmpA");
  MFloatScratchSpace tmpC(recDim, maxNoSingularityRecNghbrIds, AT_, "tmpC");
  MFloatScratchSpace weights(maxNoSingularityRecNghbrIds, AT_, "weights");
  MFloat counter = F0;
  MFloat avg = F0;
  MFloat maxc = F0;

  for(MInt i = 0; i < m_hasSingularity; ++i) {
    if(m_singularity[i].BC == -6000) {
      //      MInt totalCells=1;
      MInt len1[nDim];

      //(p)reset the reconstruction constants
      for(MInt n = 0; n < nDim + 1; ++n) {
        for(MInt m = 0; m < m_singularity[i].totalCells * m_singularity[i].Nstar; ++m) {
          m_singularity[i].ReconstructionConstants[n][m] = -999;
        }
      }

      for(MInt j = 0; j < nDim; j++) {
        len1[j] = m_singularity[i].end[j] - m_singularity[i].start[j];
        //        if(len1[j]!=0)  totalCells*=len1[j];
        ASSERT(len1[j] == 1, "");
      }

      for(MInt n = 0; n < nDim; ++n) {
        //        if(m_singularity[i].end[n]-m_singularity[i].start[n]>1) {
        //          dim=n;
        //          start[n]=m_singularity[i].start[n]+1;
        //          end[n]=m_singularity[i].end[n]-1;
        //        } else {
        start[n] = m_singularity[i].start[n];
        end[n] = m_singularity[i].end[n];
        //        }
      }

      //      for( MInt kk = start[2]; kk <end[2]; ++kk ) {
      for(MInt jj = start[1]; jj < end[1]; ++jj) {
        for(MInt ii = start[0]; ii < end[0]; ++ii) {
          MInt count = 0;
          //            MInt temp[nDim]{};
          //            temp[dim]=1;

          nghbr[count++] = cellIndex(ii, jj);
          //            nghbr[count++]=cellIndex(ii+temp[0],jj+temp[1],kk+temp[2]);

          // the coordinates of the corner where the viscousflux should be corrected.
          MInt ijk = getPointIdFromCell(ii + m_singularity[i].Viscous[0], jj + m_singularity[i].Viscous[1]);
          ijk = getPointIdFromPoint(ijk, 1, 1);

          for(MInt m = 0; m < m_singularity[i].Nstar - 1; ++m) {
            MInt* change = m_singularity[i].displacement[m];
            nghbr[count++] = cellIndex(ii + change[0], jj + change[1]);
            //              nghbr[count++]=cellIndex(ii+temp[0]+change[0],jj+temp[1]+change[1],kk+temp[2]+change[2]);
          }

          if(count != m_singularity[i].Nstar) {
            cerr << "Something wrong with the singularities in the LS coeffiecient computation" << endl;
          }

          // weighted Least square
          weights.fill(F0);

          // Compute weights with RBF (take mean distance as R0)
          for(MInt n = 0; n < count; n++) {
            MInt nghbrId = nghbr[n];
            MFloat dxdx = F0;
            for(MInt m = 0; m < nDim; ++m) {
              dxdx += POW2(m_cells->coordinates[m][nghbrId] - m_grid->m_coordinates[m][ijk]);
            }

            weights[n] = 1 / dxdx; // RBF( dxdx, POW2( dist) );
          }

          MInt id2 = ii - start[0] + (jj - start[1]) * len1[0];
          ASSERT(id2 == 0, "");
          MInt ID = id2 * m_singularity[i].Nstar;

          MFloat condNum = computeRecConstSVD(ijk, count, nghbr, ID, i, tmpA, tmpC, weights, recDim);
          avg += condNum;
          maxc = mMax(maxc, condNum);
          counter += F1;
          if(condNum < F0 || condNum > 1e7 || std::isnan(condNum)) {
            cerr << domainId() << " SVD decomposition for pointId " << ijk
                 << " with large condition number: " << condNum << " num of neighbor" << count << "x" << recDim << " "
                 << " coords " << m_grid->m_coordinates[0][ijk] << ", " << m_grid->m_coordinates[1][ijk] << endl;
          }
        }
      }
      //      }
    }
  }
}


#include <numeric>
void FvStructuredSolver2D::exchange6002() {
  if(m_rans && m_ransMethod != RANS_KEPSILON) mTerm(1, "Porous RANS computation is only supported by k-epsilon model!");
  //  if (!m_porous)
  //    mTerm(1, "bc6002 requires the property porous to be set to true!");

  // 0) Check if initBc6002 is called for the first time
  //  for(MInt bcId_=0; bcId_ < bcId; ++bcId_) {
  //    if (m_physicalBCMap[bcId_]->BC==6002) return;
  //  }

  // Determine normal vectors and save for later use
  for(MInt bcId_ = 0; bcId_ < (MInt)m_structuredBndryCnd->m_physicalBCMap.size(); ++bcId_) {
    if(m_structuredBndryCnd->m_physicalBCMap[bcId_]->BC == 6002
       && m_structuredBndryCnd->m_physicalBCMap[bcId_]->Nstar == -1) {
      MInt* start = m_structuredBndryCnd->m_physicalBCMap[bcId_]->start1;
      MInt* end = m_structuredBndryCnd->m_physicalBCMap[bcId_]->end1;

      const MInt IJKP[nDim] = {1, m_nPoints[1]};
      const MInt IJ[nDim] = {1, m_nCells[1]};

      const MInt pp[2][4] = {{0, 0, 0, 1}, {0, 0, 1, 0}};

      const MInt face = m_structuredBndryCnd->m_physicalBCMap[bcId_]->face;

      const MInt normalDir = face / 2;
      const MInt firstTangentialDir = (normalDir + 1) % nDim;
      const MInt normalDirStart = start[normalDir];
      const MInt firstTangentialStart = start[firstTangentialDir];
      const MInt firstTangentialEnd = end[firstTangentialDir];
      const MInt incp[nDim] = {IJKP[normalDir], IJKP[firstTangentialDir]};
      const MInt inc[nDim] = {IJ[normalDir], IJ[firstTangentialDir]};

      const MInt n = (face % 2) * 2 - 1;                                       //-1,+1
      const MInt g1p = normalDirStart + 2 * ((MInt)(0.5 - (0.5 * (MFloat)n))); //+2,0
      const MInt g1 = normalDirStart + (MInt)(0.5 - (0.5 * (MFloat)n));        //+1,0

      for(MInt t1 = firstTangentialStart; t1 < firstTangentialEnd; t1++) {
        const MInt cellIdG1 = g1 * inc[0] + t1 * inc[1];

        // compute four surrounding points of surface centroid
        const MInt ij = g1p * incp[0] + t1 * incp[1];
        const MInt pp1 = getPointIdFromPoint(ij, pp[normalDir][0], pp[normalDir][1]);
        const MInt pp2 = getPointIdFromPoint(ij, pp[normalDir][2], pp[normalDir][3]);
        const MInt pp3 = getPointIdFromPoint(ij, n * (1 - normalDir), n * normalDir); // point lying outside domain

        // compute the velocity of the surface centroid
        MFloat firstVec[nDim] = {F0, F0};
        MFloat normalVec[nDim] = {F0, F0};
        MFloat normalVec_[nDim]{};
        for(MInt dim = 0; dim < nDim; dim++) {
          firstVec[dim] = m_grid->m_coordinates[dim][pp2] - m_grid->m_coordinates[dim][pp1];
          normalVec_[dim] = m_grid->m_coordinates[dim][pp3] - m_grid->m_coordinates[dim][pp1];
        }

        // compute normal vector of surface
        normalVec[0] = -firstVec[1];
        normalVec[1] = firstVec[0];
        const MFloat normalLength = sqrt(POW2(normalVec[0]) + POW2(normalVec[1]));

        MFloat sgn = (std::inner_product(&normalVec[0], &normalVec[0] + nDim, &normalVec_[0], 0.0) < 0.0) ? -1 : 1;
        if(m_blockType == "fluid") sgn *= -1;

        for(MInt dim = 0; dim < nDim; dim++) {
          normalVec[dim] /= normalLength;
          m_cells->fq[FQ->NORMAL[dim]][cellIdG1] = sgn * normalVec[dim];
        }
      }
    }
  }

  // Determine normal vector at singularities
  for(MInt i = 0; i < m_hasSingularity; ++i) {
    const auto& singularity = m_singularity[i];
    // only correct for bc 6000 not for bc 4000-5000
    if(singularity.BC == -6000) {
      MBool takeIt = false;
      for(MInt n = 0; n < singularity.Nstar; ++n) {
        if(singularity.BCsingular[n] == -6002) takeIt = true;
      }

      if(takeIt) {
        MInt start[nDim], end[nDim];
        for(MInt n = 0; n < nDim; ++n) {
          if(singularity.end[n] - singularity.start[n] > 1) {
            mTerm(1, "In 2D not possible!");
            // dim=n;
            // start[n]=singularity.start[n]+1;
            start[n] = singularity.start[n] + 1;
            end[n] = singularity.end[n] - 1;
          } else {
            start[n] = singularity.start[n];
            end[n] = singularity.end[n];
          }
        }

        for(MInt jj = start[1]; jj < end[1]; ++jj) {
          for(MInt ii = start[0]; ii < end[0]; ++ii) {
            const MInt IJ = cellIndex(ii, jj);
            if(abs(m_cells->fq[FQ->NORMAL[0]][IJ]) > 1e-8) mTerm(1, "");
            for(MInt m = 0; m < 2; ++m) {
              const MInt* change = singularity.displacement[m];
              const MInt nghbr = cellIndex(ii + change[0], jj + change[1]);
              for(MInt d = 0; d < nDim; ++d)
                m_cells->fq[FQ->NORMAL[d]][IJ] += m_cells->fq[FQ->NORMAL[d]][nghbr];
            }
          }
        }
      }
    }
  }

  // TODO_SS labels:FV,toenhance The exchange of all the normals is an overhead, because it is only needed at
  // singularity points
  /////////////////////////////////////////////////////////////////////////////
  /// GATHER & SEND
  /////////////////////////////////////////////////////////////////////////////
  //  MInt sendSizeTotal = 0;
  std::vector<MInt> receiveSizes;
  for(auto& snd : m_sndComm) {
    // TODO_SS labels:FV right now exchange at all 6000er not only 6002 because of singularities
    if(snd->bcId == -6000 || snd->bcId == -6002) {
      // Gather

      MInt size = 1;
      for(MInt dim = 0; dim < nDim; ++dim)
        size *= snd->endInfoCells[dim] - snd->startInfoCells[dim];
      //      std::vector<MFloat> bufferSnd(size*(1+nDim));
      //      sendSizeTotal += size;

      MInt pos = 0;
      for(MInt j = snd->startInfoCells[1]; j < snd->endInfoCells[1]; j++) {
        for(MInt i = snd->startInfoCells[0]; i < snd->endInfoCells[0]; i++) {
          const MInt cellId = i + (j * m_nCells[1]);
          // TODO_SS labels:FV Latter only allocate FQ->POROSITY for m_blockType==porous
          //        if (m_blockType=="fluid")
          //          bufferSnd[pos++] = 1;
          //        else
          snd->cellBuffer[pos] = m_cells->fq[FQ->POROSITY][cellId];
          for(MInt d = 0; d < nDim; ++d) {
            snd->cellBuffer[(1 + d) * size + pos] = m_cells->fq[FQ->NORMAL[d]][cellId];
          }
          ++pos;
        }
      }

      // Send
      MInt tag = domainId() + (snd->tagHelper) * noDomains();
      MInt err = MPI_Isend((void*)&snd->cellBuffer[0], size * (nDim + 1), MPI_DOUBLE, snd->nghbrId, tag,
                           m_StructuredComm, &snd->mpi_request, AT_, "(void*)&snd->cellBuffer[0]");
      if(err) cout << "rank " << domainId() << " sending throws error " << endl;


      // Determine size of receive buffer
      // TODO_SS labels:FV right now exchange at all 6000er not only 6002 because of singularities
      size = 1;
      for(MInt dim = 0; dim < nDim; ++dim)
        size *= snd->endInfoCells[dim] - snd->startInfoCells[dim];
      receiveSizes.push_back(size);
    }
  }

  /////////////////////////////////////////////////////////////////////////////
  /// RECEIVE
  /////////////////////////////////////////////////////////////////////////////
  std::vector<MFloat> bufferRcv(std::accumulate(receiveSizes.begin(), receiveSizes.end(), 0) * (nDim + 1));
  MInt offset = 0;
  MInt cnt = 0;
  for(auto& rcv : m_rcvComm) {
    // TODO_SS labels:FV right now exchange at all 6000er not only 6002 because of singularities
    if(rcv->bcId == -6002 || rcv->bcId == -6000) {
      const MInt rcvSize = receiveSizes[cnt];
      MInt tag = rcv->nghbrId + (rcv->tagHelper) * noDomains();
      MInt err = MPI_Irecv(&bufferRcv[offset], rcvSize * (1 + nDim), MPI_DOUBLE, rcv->nghbrId, tag, m_StructuredComm,
                           &rcv->mpi_request, AT_, "(void*)&rcvSize");
      if(err) cout << "rank " << domainId() << " sending throws error " << endl;

      offset += rcvSize * (1 + nDim);
      ++cnt;
    }
  }

  /////////////////////////////////////////////////////////////////////////////
  /// WAIT
  /////////////////////////////////////////////////////////////////////////////
  for(auto& snd : m_sndComm) {
    if(snd->bcId == -6002 || snd->bcId == -6000) {
      MPI_Wait(&(snd->mpi_request), &(snd->mpi_status), AT_);
    }
  }


  for(auto& rcv : m_rcvComm) {
    if(rcv->bcId == -6002 || rcv->bcId == -6000) {
      MPI_Wait(&(rcv->mpi_request), &(rcv->mpi_status), AT_);
    }
  }

  /////////////////////////////////////////////////////////////////////////////
  /// SCATTER
  /////////////////////////////////////////////////////////////////////////////
  // TODO_SS labels:FV,toremove
  ScratchSpace<MFloat> normals_temp(m_noCells * nDim, AT_, "normal_temp");
  offset = 0;
  cnt = 0;
  for(auto& rcv : m_rcvComm) {
    if(rcv->bcId == -6002 || rcv->bcId == -6000) {
      // TODO_SS labels:FV right now exchange at all 6000er not only 6002 because of singularities

      MInt j2, i2, id2;
      MInt step2[nDim];
      MInt start1[nDim];
      MInt start2[nDim];
      MInt end2[nDim];
      MInt len2[nDim];
      MInt totalCells = 1;
      MInt len1[nDim];

      for(MInt j = 0; j < nDim; j++) {
        len1[j] = rcv->endInfoCells[j] - rcv->startInfoCells[j];
        if(len1[j] != 0) totalCells *= len1[j];
        // added    check the step for RCV part !!!!!!!!important
        step2[rcv->orderInfo[j]] = rcv->stepInfo[j];
      }

      // Sanity check
      ASSERT(totalCells == receiveSizes[cnt], "");

      for(MInt j = 0; j < nDim; j++) {
        start2[j] = 0;
        end2[j] = len1[j] - 1;
        len2[rcv->orderInfo[j]] = len1[j];
        if(step2[j] < 0) {
          MInt dummy = start2[j];
          start2[j] = end2[j];
          end2[j] = dummy;
        }
      }

      MInt pos = 0;
      j2 = start2[1];
      for(MInt j = rcv->startInfoCells[1]; j < rcv->endInfoCells[1]; j++) {
        i2 = start2[0];
        for(MInt i = rcv->startInfoCells[0]; i < rcv->endInfoCells[0]; i++) {
          start1[rcv->orderInfo[0]] = i2;
          start1[rcv->orderInfo[1]] = j2;

          id2 = start1[0] + start1[1] * len2[0];
          const MInt cellId = i + (j * m_nCells[1]);
          m_cells->fq[FQ->POROSITY][cellId] = bufferRcv[offset * (nDim + 1) + id2];
          for(MInt d = 0; d < nDim; ++d) {
            normals_temp[m_noCells * d + cellId] = bufferRcv[offset * (nDim + 1) + (d + 1) * totalCells + id2];
            //            m_cells->fq[FQ->NORMAL[d]][cellId] = bufferRcv[offset*(nDim+1)+(d+1)*totalCells+id2];
          }

          i2 += step2[0];
          pos++;
        }
        j2 += step2[1];
      }

      offset += totalCells;
      ++cnt;
    }
  }

  // Determine normal vector at singularities
  for(MInt i = 0; i < m_hasSingularity; ++i) {
    const auto& singularity = m_singularity[i];
    // only correct for bc 6000 not for bc 4000-5000
    if(singularity.BC == -6000) {
      MBool takeIt = false;
      for(MInt n = 0; n < singularity.Nstar; ++n) {
        if(singularity.BCsingular[n] == -6002) takeIt = true;
      }

      if(takeIt) {
        MInt start[nDim], end[nDim];
        for(MInt n = 0; n < nDim; ++n) {
          if(singularity.end[n] - singularity.start[n] > 1) {
            mTerm(1, "In 2D not possible!");
            // dim=n;
            // start[n]=singularity.start[n]+1;
            start[n] = singularity.start[n] + 1;
            end[n] = singularity.end[n] - 1;
          } else {
            start[n] = singularity.start[n];
            end[n] = singularity.end[n];
          }
        }

        const MInt nstar = singularity.Nstar;

        for(MInt jj = start[1]; jj < end[1]; ++jj) {
          for(MInt ii = start[0]; ii < end[0]; ++ii) {
            const MInt IJ = cellIndex(ii, jj);
            MFloat temp[nDim];
            for(MInt d = 0; d < nDim; ++d) {
              temp[d] = m_cells->fq[FQ->NORMAL[d]][IJ];
            }
            for(MInt m = 0; m < nstar - 1; ++m) {
              const MInt* change = singularity.displacement[m];
              const MInt nghbr = cellIndex(ii + change[0], jj + change[1]);
              for(MInt d = 0; d < nDim; ++d)
                temp[d] += normals_temp[m_noCells * d + nghbr];
            }

            MFloat l = 0;
            for(MInt d = 0; d < nDim; ++d) {
              m_cells->fq[FQ->NORMAL[d]][IJ] = temp[d] / (2 * nstar);
              l += POW2(m_cells->fq[FQ->NORMAL[d]][IJ]);
            }
            l = sqrt(l);
            for(MInt d = 0; d < nDim; ++d) {
              m_cells->fq[FQ->NORMAL[d]][IJ] /= l;
            }

            for(MInt m = 0; m < nstar - 1; ++m) {
              const MInt* change = singularity.displacement[m];
              const MInt nghbr = cellIndex(ii + change[0], jj + change[1]);
              for(MInt d = 0; d < nDim; ++d)
                m_cells->fq[FQ->NORMAL[d]][nghbr] = m_cells->fq[FQ->NORMAL[d]][IJ];
            }
          }
        }
      }
    }
  }

  ////////////////////////////////////////////////////////////////////////////////
#if 0
  // 1) Gather
  std::vector<MFloat> sendBuffer;
  std::vector<MInt> sendcounts;
  std::vector<MInt> snghbrs;
  std::vector<MInt> tags;
  for(auto& snd: m_sndComm) {    
    if (snd->bcId==-6002) {
      snghbrs.push_back(snd->nghbrId);
      tags.push_back(domainId()+(snd->tagHelper)*m_solver->noDomains());

      MInt size = 1;
      for (MInt dim = 0; dim < nDim; ++dim)
        size *= snd->endInfoCells[dim] - snd->startInfoCells[dim];
      sendcounts.push_back(size);

      sendBuffer.resize(pos+size);
      for(MInt j=snd->startInfoCells[1]; j<snd->endInfoCells[1]; j++) {
        for(MInt i=snd->startInfoCells[0]; i<snd->endInfoCells[0]; i++) {
          const MInt cellId = cellIndex(i,j);
          // TODO_SS labels:FV Latter only allocate FQ->POROSITY for m_blockType==porous
	  //          if (m_blockType=="fluid")
	  //            sendBuffer[pos++] = 1;
	  //          else
          sendBuffer[pos++] = m_cells->fq[FQ->POROSITY][cellId];
        }
      }
    }
  }

  const MInt noNeighborDomains = snghbrs.size();

  // 2) Send & receive
  std::vector<MInt> recvcounts(noNeighborDomains);
  std::vector<MInt> rdispls(noNeighborDomains);
  std::vector<MFloat> recvBuffer = maia::mpi::mpiExchangePointToPoint(&sendBuffer[0],
                                                                      &snghbrs[0],
                                                                      noNeighborDomains,
                                                                      &sendcounts[0],
                                                                      &snghbrs[0],
                                                                      noNeighborDomains,
                                                                      m_StructuredComm,
                                                                      m_solver->domainId(),
                                                                      1,
                                                                      recvcounts.data(),
                                                                      rdispls.data());
  // 2.5) Send & receive the tags
  std::vector<MInt> sendcounts2(noNeighborDomains, 1);
  std::vector<MInt> recvTags = maia::mpi::mpiExchangePointToPoint(&tags[0],
                                                                  &snghbrs[0],
                                                                  noNeighborDomains,
                                                                  &sendcounts2[0],
                                                                  &snghbrs[0],
                                                                  noNeighborDomains,
                                                                  m_StructuredComm,
                                                                  m_solver->domainId(),
                                                                  1);

  // 3) Scatter
  for(auto& rcv: m_rcvComm) {
    if (rcv->bcId==-6002) {
      const MInt tag = rcv->nghbrId+rcv->tagHelper*m_solver->noDomains();
      MInt n;
      for (n = 0; n < noNeighborDomains; ++n) {
        if (tag==recvTags[n])
          break;
      }
      if (n==noNeighborDomains) mTerm(1, "n == noNeighborDomains");

      const MFloat* const recvBuffer_ = &recvBuffer[rdispls[n]];
      const MInt noReceivedElements = recvcounts[n];

      /////////// following is analoge to  FvStructuredSolver2D::scatter()
      MInt j2, i2, id2;
      MInt  step2[nDim];
      MInt start1[nDim];
      MInt start2[nDim];
      MInt end2[nDim];
      MInt len2[nDim];
      MInt totalCells=1;
      MInt len1[nDim];

      for(MInt j=0; j<nDim; j++) {
        len1[j]=rcv->endInfoCells[j] - rcv->startInfoCells[j];
        if(len1[j]!=0) totalCells*=len1[j];
        //added    check the step for RCV part !!!!!!!!important
        step2[rcv->orderInfo[j]]=rcv->stepInfo[j];
      }

      //TODO_SS labels:FV check if this assert makes sense
      ASSERT(noReceivedElements==totalCells, "noReceivedElements==totalCells");

      for(MInt j=0; j<nDim; j++) {
        start2[j]=0;
        end2[j]=len1[j]-1;
        len2[rcv->orderInfo[j]]=len1[j];
        if(step2[j]<0) {
          MInt dummy=start2[j];
          start2[j]=end2[j];
          end2[j]=dummy;
        }
      }
      
      MInt* startInfo=rcv->startInfoCells;
      MInt* endInfo=rcv->endInfoCells;      

      MInt pos=0;
      j2=start2[1];
      for(MInt j=startInfo[1]; j<endInfo[1]; j++) {
        i2=start2[0];
        for(MInt i=startInfo[0]; i<endInfo[0]; i++) {
          start1[rcv->orderInfo[0]]=i2;
          start1[rcv->orderInfo[1]]=j2;

          id2=start1[0]+start1[1]*len2[0];
          const MInt cellId = i +(j*m_nCells[1]);
          m_cells->fq[FQ->POROSITY][cellId]= recvBuffer_[id2];

          i2+=step2[0];
          pos++;
        }
        j2+=step2[1];
      }
    }
  }
#endif
  ////////////////////////////////////////////////////////////////////////////////
}


inline MFloat FvStructuredSolver2D::getPSI(MInt I, MInt dim) {
  const MFloat FK = 18.0;
  const MInt IJK[2] = {1, m_nCells[1]};
  const MInt IP1 = I + IJK[dim];
  const MInt IM1 = I - IJK[dim];
  const MInt IP2 = I + 2 * IJK[dim];

  const MFloat PIM2 = m_cells->pvariables[PV->P][IM1];
  const MFloat PIM1 = m_cells->pvariables[PV->P][I];
  const MFloat PIP2 = m_cells->pvariables[PV->P][IP2];
  const MFloat PIP1 = m_cells->pvariables[PV->P][IP1];

  const MFloat PSI =
      mMin(F1, FK
                   * mMax(mMax(fabs((PIM2 - PIM1) / mMin(PIM2, PIM1)), fabs((PIM1 - PIP1) / mMin(PIM1, PIP1))),
                          fabs((PIP1 - PIP2) / mMin(PIP1, PIP2))));
  return PSI;
}
