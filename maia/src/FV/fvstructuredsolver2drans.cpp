// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "fvstructuredsolver2drans.h"
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

class FvStructuredSolver2D;

FvStructuredSolver2DRans::FvStructuredSolver2DRans(FvStructuredSolver2D* solver)
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
  // set all the methods for the 2d RANS code here

  initFluxMethod();
}

FvStructuredSolver2DRans::~FvStructuredSolver2DRans() {}

void FvStructuredSolver2DRans::initFluxMethod() {
  // set pointer to the right AUSM for the RANS equations;

  /*! \property
    \page propertiesFVSTRCTRD
    \section ransMethod
    <code>RansMethod FvStructuredSolver2DRans::m_ransMethod </code>\n
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
    case RANS_SA:
    case RANS_SA_DV: {
      mAlloc(m_cells->saFlux1, nDim, m_noCells, "m_cells->saFlux1", -999999.9, AT_);
      mAlloc(m_cells->saFlux2, nDim, m_noCells, "m_cells->saFlux2", -999999.9, AT_);
      mAlloc(m_cells->prodDest, m_noCells, "m_cells->prodDest", -999999.9, AT_);
      compTurbVisc = &FvStructuredSolver2DRans::computeTurbViscosity_SA;
      viscFluxMethod = &FvStructuredSolver2DRans::viscousFlux_SA;
      switch(m_solver->CV->noVariables) {
        case 5: {
          reconstructSurfaceData = &FvStructuredSolver2DRans::Muscl_AusmDV<5, 1>;
          break;
        }
        case 6: {
          reconstructSurfaceData = &FvStructuredSolver2DRans::Muscl_AusmDV<6, 1>;
          break;
        }
        case 7: {
          reconstructSurfaceData = &FvStructuredSolver2DRans::Muscl_AusmDV<7, 1>;
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
    case RANS_KEPSILON: {
      setAndAllocate_KEPSILON();

      compTurbVisc = &FvStructuredSolver2DRans::computeTurbViscosity_KEPSILON;
      viscFluxMethod = &FvStructuredSolver2DRans::viscousFlux_KEPSILON2;

      // TODO_SS labels:FV
      MBool AusmDV = true;
      AusmDV = Context::getSolverProperty<MBool>("AusmDV", m_solverId, AT_, &AusmDV);

      if(m_solver->m_limiter) {
        mAlloc(m_cells->ql, PV->noVariables, m_noCells, "m_cells->ql", F0, AT_);
        mAlloc(m_cells->qr, PV->noVariables, m_noCells, "m_cells->qr", F0, AT_);
        m_log << "Using RANS with K-Epsilon model and limited AusmDV" << endl;
        switch(m_solver->CV->noVariables) {
          case 6: {
            if(m_solver->m_rans2eq_mode == "production")
              reconstructSurfaceData = (AusmDV) ? &FvStructuredSolver2DRans::Muscl_AusmDV_Limited<6, 2, true>
                                                : &FvStructuredSolver2DRans::Muscl_Ausm_Limited<6, 2, true>;
            else
              reconstructSurfaceData = (AusmDV) ? &FvStructuredSolver2DRans::Muscl_AusmDV_Limited<6, 2>
                                                : &FvStructuredSolver2DRans::Muscl_Ausm_Limited<6, 2>;
            break;
          }
          case 7: {
            if(m_solver->m_rans2eq_mode == "production")
              reconstructSurfaceData = (AusmDV) ? &FvStructuredSolver2DRans::Muscl_AusmDV_Limited<7, 2, true>
                                                : &FvStructuredSolver2DRans::Muscl_Ausm_Limited<7, 2, true>;
            else
              reconstructSurfaceData = (AusmDV) ? &FvStructuredSolver2DRans::Muscl_AusmDV_Limited<7, 2>
                                                : &FvStructuredSolver2DRans::Muscl_Ausm_Limited<7, 2>;
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
        if(!AusmDV) mTerm(1, "Not implemented yet!");
        m_log << "Using RANS with K-Epsilon model and AusmDV" << endl;
        switch(m_solver->CV->noVariables) {
          case 6: {
            if(m_solver->m_rans2eq_mode == "production")
              reconstructSurfaceData = &FvStructuredSolver2DRans::Muscl_AusmDV<6, 2, true>;
            else
              reconstructSurfaceData = &FvStructuredSolver2DRans::Muscl_AusmDV<6, 2>;
            break;
          }
          case 7: {
            if(m_solver->m_rans2eq_mode == "production")
              reconstructSurfaceData = &FvStructuredSolver2DRans::Muscl_AusmDV<7, 2, true>;
            else
              reconstructSurfaceData = &FvStructuredSolver2DRans::Muscl_AusmDV<7, 2>;
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

void FvStructuredSolver2DRans::Muscl(MInt) { (this->*reconstructSurfaceData)(); }

template <MInt noVars, MInt noRANS, MBool rans2eq_production_mode>
void FvStructuredSolver2DRans::Muscl_AusmDV() {
  TRACE();

  // stencil identifier
  const MInt IJ[2] = {1, m_nCells[1]};

  const MFloat* const RESTRICT x = ALIGNED_F(m_cells->coordinates[0]);
  const MFloat* const RESTRICT y = ALIGNED_F(m_cells->coordinates[1]);
  const MFloat* const* const RESTRICT vars = ALIGNED_F(m_cells->pvariables);
  MFloat* const* const RESTRICT dss = ALIGNED_F(m_cells->dss);
  MFloat* const* const RESTRICT flux = ALIGNED_F(m_cells->flux);
  MFloat* const RESTRICT cellRhs = ALIGNED_MF(m_cells->rightHandSide[0]);

  const MUint noCells = m_noCells;
  const MFloat gamma = m_solver->m_gamma;
  const MFloat gammaMinusOne = gamma - 1.0;
  const MFloat FgammaMinusOne = F1 / gammaMinusOne;

  const MInt noCellsI = m_nCells[1] - 2;
  const MInt noCellsJ = m_nCells[0] - 2;

  const MInt noCellsIP1 = m_nCells[1] - 1;
  const MInt noCellsJP1 = m_nCells[0] - 1;

  if(!m_dsIsComputed) {
    for(MInt dim = 0; dim < nDim; dim++) {
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for(MInt j = 0; j < noCellsJP1; j++) {
        for(MInt i = 0; i < noCellsIP1; i++) {
          const MInt I = cellIndex(i, j);
          const MInt IP1 = I + IJ[dim];
          dss[dim][I] = sqrt(POW2(x[IP1] - x[I]) + POW2(y[IP1] - y[I]));
        }
      }
    }

    m_dsIsComputed = true;
  }


  for(MInt dim = 0; dim < nDim; dim++) {
    for(MInt j = 1; j < noCellsJ; j++) {
      for(MInt i = 1; i < noCellsI; i++) {
        const MInt I = cellIndex(i, j);
        const MInt IP1 = I + IJ[dim];
        const MInt IM1 = I - IJ[dim];
        const MInt IP2 = I + 2 * IJ[dim];

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

        // in C++17 and later use "if constexpr"
        // for SA-model: nutilde
        // for k-epsilon-model: k and epsilon
        MFloat RANSL[noRANS];
        MFloat RANSR[noRANS];
        for(MInt v = 0; v < noRANS; ++v) {
          const MFloat DQRANSVAR = vars[PV->RANS_VAR[v]][IP1] - vars[PV->RANS_VAR[v]][I];
          const MFloat DQPRANSVAR = vars[PV->RANS_VAR[v]][IP2] - vars[PV->RANS_VAR[v]][IP1];
          const MFloat DQMRANSVAR = vars[PV->RANS_VAR[v]][I] - vars[PV->RANS_VAR[v]][IM1];
          RANSL[v] = vars[PV->RANS_VAR[v]][I] + DSM * (DSM1 * DQRANSVAR + DS * DQMRANSVAR);
          RANSR[v] = vars[PV->RANS_VAR[v]][IP1] - DSP * (DS * DQPRANSVAR + DSP1 * DQRANSVAR);
        }
        //        const MFloat DQNUTILDE   = vars[PV->RANS_VAR[0]][IP1]-vars[PV->RANS_VAR[0]][I];
        //        const MFloat DQPNUTILDE = vars[PV->RANS_VAR[0]][IP2]-vars[PV->RANS_VAR[0]][IP1];
        //        const MFloat DQMNUTILDE = vars[PV->RANS_VAR[0]][I]-vars[PV->RANS_VAR[0]][IM1];
        //        const MFloat NUTILDEL = vars[PV->RANS_VAR[0]][I]+DSM*(DSM1*DQNUTILDE+DS*DQMNUTILDE);
        //        const MFloat NUTILDER = vars[PV->RANS_VAR[0]][IP1]-DSP*(DS*DQPNUTILDE+DSP1*DQNUTILDE);

        const MFloat surf0 = m_cells->surfaceMetrics[dim * 2 + 0][I];
        const MFloat surf1 = m_cells->surfaceMetrics[dim * 2 + 1][I];

        const MFloat FRHOL = F1 / RHOL;
        const MFloat FRHOR = F1 / RHOR;

        const MFloat PLfRHOL = PL / RHOL;
        const MFloat PRfRHOR = PR / RHOR;
        const MFloat e0 = (rans2eq_production_mode)
                              ? PLfRHOL * FgammaMinusOne + 0.5 * (POW2(UL) + POW2(VL)) + RANSL[0] + PLfRHOL
                              : PLfRHOL * FgammaMinusOne + 0.5 * (POW2(UL) + POW2(VL)) + PLfRHOL;
        const MFloat e1 = (rans2eq_production_mode)
                              ? PRfRHOR * FgammaMinusOne + 0.5 * (POW2(UR) + POW2(VR)) + RANSR[0] + PRfRHOR
                              : PRfRHOR * FgammaMinusOne + 0.5 * (POW2(UR) + POW2(VR)) + PRfRHOR;


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

        const MFloat FDV = 0.3;
        const MFloat DXDXEZ = m_cells->coordinates[0][IP1] - m_cells->coordinates[0][I];
        const MFloat DYDXEZ = m_cells->coordinates[1][IP1] - m_cells->coordinates[1][I];
        MFloat SV = 2.0 * DGRAD / (m_cells->cellJac[I] + m_cells->cellJac[IP1]) * (FDV + (F1 - FDV) * getPSI(I, dim));
        const MFloat SV1 = F1 * SV * DXDXEZ;
        const MFloat SV2 = F1 * SV * DYDXEZ;

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

        flux[CV->RHO_U][I] = RHOU2 * (UL + UR) + ARHOU2 * (UL - UR) + PLR * surf0 + RHOUL * UUL2 + RHOUR * UUR2;
        flux[CV->RHO_V][I] = RHOU2 * (VL + VR) + ARHOU2 * (VL - VR) + PLR * surf1 + RHOUL * UUL3 + RHOUR * UUR3;
        flux[CV->RHO_E][I] = RHOU2 * (e0 + e1) + ARHOU2 * (e0 - e1);
        flux[CV->RHO][I] = RHOU;
        for(MInt v = 0; v < noRANS; ++v) {
          flux[CV->RANS_VAR[v]][I] = RHOU2 * (RANSL[v] + RANSR[v]) + ARHOU2 * (RANSL[v] - RANSR[v]);
          //        flux[I+noCells*CV->RANS_VAR[0]] = RHOU2*(NUTILDEL+NUTILDER) + ARHOU2*(NUTILDEL-NUTILDER);
        }
      }
    }

    // FLUX BALANCE
    for(MUint v = 0; v < noVars; v++) {
      for(MInt j = m_noGhostLayers; j < m_nCells[0] - m_noGhostLayers; j++) {
        for(MInt i = m_noGhostLayers; i < m_nCells[1] - m_noGhostLayers; i++) {
          const MInt I = cellIndex(i, j);
          const MInt IM1 = I - IJ[dim];
          const MUint offset = v * noCells;
          MFloat* const RESTRICT rhs = ALIGNED_F(cellRhs + offset);
          rhs[I] += flux[v][IM1] - flux[v][I];
        }
      }
    }
  }
}
template void FvStructuredSolver2DRans::Muscl_AusmDV<5, 1>();
template void FvStructuredSolver2DRans::Muscl_AusmDV<6, 1>();
template void FvStructuredSolver2DRans::Muscl_AusmDV<7, 1>();
template void FvStructuredSolver2DRans::Muscl_AusmDV<6, 2>();
template void FvStructuredSolver2DRans::Muscl_AusmDV<7, 2>();


template <MInt noVars, MInt noRANS, MBool rans2eq_production_mode>
void FvStructuredSolver2DRans::Muscl_AusmDV_Limited() {
  TRACE();

  // stencil identifier
  const MInt IJ[2] = {1, m_nCells[1]};

  const MFloat* const RESTRICT x = ALIGNED_F(m_cells->coordinates[0]);
  const MFloat* const RESTRICT y = ALIGNED_F(m_cells->coordinates[1]);
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

  const MInt noCellsI = m_nCells[1] - 2;
  const MInt noCellsJ = m_nCells[0] - 2;

  const MInt noCellsIP1 = m_nCells[1] - 1;
  const MInt noCellsJP1 = m_nCells[0] - 1;

  if(!m_dsIsComputed) {
    for(MInt dim = 0; dim < nDim; dim++) {
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for(MInt j = 0; j < noCellsJP1; j++) {
        for(MInt i = 0; i < noCellsIP1; i++) {
          const MInt I = cellIndex(i, j);
          const MInt IP1 = I + IJ[dim];
          dss[dim][I] = sqrt(POW2(x[IP1] - x[I]) + POW2(y[IP1] - y[I]));
        }
      }
    }

    m_dsIsComputed = true;
  }

  for(MInt dim = 0; dim < nDim; dim++) {
    for(MInt j = 1; j < noCellsJP1; j++) {
      for(MInt i = 1; i < noCellsIP1; i++) {
        const MInt I = cellIndex(i, j);
        const MInt IP1 = I + IJ[dim];
        const MInt IM1 = I - IJ[dim];

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

        for(MInt v = 0; v < noRANS; ++v) {
          const MFloat DQLRANSVAR = DSL * (vars[PV->RANS_VAR[v]][IP1] - vars[PV->RANS_VAR[v]][I]);
          const MFloat DQRRANSVAR = DSR * (vars[PV->RANS_VAR[v]][I] - vars[PV->RANS_VAR[v]][IM1]);
          const MFloat PHIRANSVAR = mMax(
              F0, DQLRANSVAR * DQRRANSVAR / (POW2(DQLRANSVAR) + POW2(DQRRANSVAR) + EPSLIM * mMax(EPSS, DSL * DSR)));
          ql[PV->RANS_VAR[v]][I] = vars[PV->RANS_VAR[v]][I] + FDS * (DQLRANSVAR + DQRRANSVAR) * PHIRANSVAR;
          qr[PV->RANS_VAR[v]][IM1] = vars[PV->RANS_VAR[v]][I] - FDS * (DQLRANSVAR + DQRRANSVAR) * PHIRANSVAR;
        }
      }
    }


    for(MInt j = 1; j < noCellsJ; j++) {
      for(MInt i = 1; i < noCellsI; i++) {
        const MInt I = cellIndex(i, j);
        const MInt IP1 = I + IJ[dim];

        MFloat UL = ql[PV->U][I];
        MFloat UR = qr[PV->U][I];

        MFloat VL = ql[PV->V][I];
        MFloat VR = qr[PV->V][I];

        const MFloat PL = ql[PV->P][I];
        const MFloat PR = qr[PV->P][I];

        const MFloat RHOL = ql[PV->RHO][I];
        const MFloat RHOR = qr[PV->RHO][I];

        MFloat RANSL[noRANS];
        MFloat RANSR[noRANS];
        for(MInt v = 0; v < noRANS; ++v) {
          RANSL[v] = ql[PV->RANS_VAR[v]][I];
          RANSR[v] = qr[PV->RANS_VAR[v]][I];
        }

        const MFloat surf0 = m_cells->surfaceMetrics[dim * 2 + 0][I];
        const MFloat surf1 = m_cells->surfaceMetrics[dim * 2 + 1][I];

        const MFloat FRHOL = F1 / RHOL;
        const MFloat FRHOR = F1 / RHOR;

        const MFloat PLfRHOL = PL / RHOL;
        const MFloat PRfRHOR = PR / RHOR;
        const MFloat e0 = (rans2eq_production_mode)
                              ? PLfRHOL * FgammaMinusOne + 0.5 * (POW2(UL) + POW2(VL)) + RANSL[0] + PLfRHOL
                              : PLfRHOL * FgammaMinusOne + 0.5 * (POW2(UL) + POW2(VL)) + PLfRHOL;
        const MFloat e1 = (rans2eq_production_mode)
                              ? PRfRHOR * FgammaMinusOne + 0.5 * (POW2(UR) + POW2(VR)) + RANSR[0] + PRfRHOR
                              : PRfRHOR * FgammaMinusOne + 0.5 * (POW2(UR) + POW2(VR)) + PRfRHOR;


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

        const MFloat FDV = 0.3;
        const MFloat DXDXEZ = m_cells->coordinates[0][IP1] - m_cells->coordinates[0][I];
        const MFloat DYDXEZ = m_cells->coordinates[1][IP1] - m_cells->coordinates[1][I];
        MFloat SV = 2.0 * DGRAD / (m_cells->cellJac[I] + m_cells->cellJac[IP1]) * (FDV + (F1 - FDV) * getPSI(I, dim));
        const MFloat SV1 = F1 * SV * DXDXEZ;
        const MFloat SV2 = F1 * SV * DYDXEZ;

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

        flux[CV->RHO_U][I] = RHOU2 * (UL + UR) + ARHOU2 * (UL - UR) + PLR * surf0 + RHOUL * UUL2 + RHOUR * UUR2;
        flux[CV->RHO_V][I] = RHOU2 * (VL + VR) + ARHOU2 * (VL - VR) + PLR * surf1 + RHOUL * UUL3 + RHOUR * UUR3;
        flux[CV->RHO_E][I] = RHOU2 * (e0 + e1) + ARHOU2 * (e0 - e1);
        flux[CV->RHO][I] = RHOU;
        for(MInt v = 0; v < noRANS; ++v) {
          flux[CV->RANS_VAR[v]][I] = RHOU2 * (RANSL[v] + RANSR[v]) + ARHOU2 * (RANSL[v] - RANSR[v]);
        }
      }
    }

    // FLUX BALANCE
    for(MUint v = 0; v < noVars; v++) {
      for(MInt j = m_noGhostLayers; j < m_nCells[0] - m_noGhostLayers; j++) {
        for(MInt i = m_noGhostLayers; i < m_nCells[1] - m_noGhostLayers; i++) {
          const MInt I = cellIndex(i, j);
          const MInt IM1 = I - IJ[dim];
          const MUint offset = v * noCells;
          MFloat* const RESTRICT rhs = ALIGNED_F(cellRhs + offset);
          rhs[I] += flux[v][IM1] - flux[v][I];
        }
      }
    }
  }
}
template void FvStructuredSolver2DRans::Muscl_AusmDV_Limited<6, 2>();
template void FvStructuredSolver2DRans::Muscl_AusmDV_Limited<7, 2>();


template <MInt noVars, MInt noRANS, MBool rans2eq_production_mode>
void FvStructuredSolver2DRans::Muscl_Ausm_Limited() {
  TRACE();

  // stencil identifier
  const MInt IJ[2] = {1, m_nCells[1]};

  const MFloat* const RESTRICT x = ALIGNED_F(m_cells->coordinates[0]);
  const MFloat* const RESTRICT y = ALIGNED_F(m_cells->coordinates[1]);
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

  const MInt noCellsI = m_nCells[1] - 2;
  const MInt noCellsJ = m_nCells[0] - 2;

  const MInt noCellsIP1 = m_nCells[1] - 1;
  const MInt noCellsJP1 = m_nCells[0] - 1;

  if(!m_dsIsComputed) {
    for(MInt dim = 0; dim < nDim; dim++) {
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for(MInt j = 0; j < noCellsJP1; j++) {
        for(MInt i = 0; i < noCellsIP1; i++) {
          const MInt I = cellIndex(i, j);
          const MInt IP1 = I + IJ[dim];
          dss[dim][I] = sqrt(POW2(x[IP1] - x[I]) + POW2(y[IP1] - y[I]));
        }
      }
    }

    m_dsIsComputed = true;
  }

  for(MInt dim = 0; dim < nDim; dim++) {
    for(MInt j = 1; j < noCellsJP1; j++) {
      for(MInt i = 1; i < noCellsIP1; i++) {
        const MInt I = cellIndex(i, j);
        const MInt IP1 = I + IJ[dim];
        const MInt IM1 = I - IJ[dim];

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

        for(MInt v = 0; v < noRANS; ++v) {
          const MFloat DQLRANSVAR = DSL * (vars[PV->RANS_VAR[v]][IP1] - vars[PV->RANS_VAR[v]][I]);
          const MFloat DQRRANSVAR = DSR * (vars[PV->RANS_VAR[v]][I] - vars[PV->RANS_VAR[v]][IM1]);
          const MFloat PHIRANSVAR = mMax(
              F0, DQLRANSVAR * DQRRANSVAR / (POW2(DQLRANSVAR) + POW2(DQRRANSVAR) + EPSLIM * mMax(EPSS, DSL * DSR)));
          ql[PV->RANS_VAR[v]][I] = vars[PV->RANS_VAR[v]][I] + FDS * (DQLRANSVAR + DQRRANSVAR) * PHIRANSVAR;
          qr[PV->RANS_VAR[v]][IM1] = vars[PV->RANS_VAR[v]][I] - FDS * (DQLRANSVAR + DQRRANSVAR) * PHIRANSVAR;
        }
      }
    }


    for(MInt j = 1; j < noCellsJ; j++) {
      for(MInt i = 1; i < noCellsI; i++) {
        const MInt I = cellIndex(i, j);

        MFloat UL = ql[PV->U][I];
        MFloat UR = qr[PV->U][I];

        MFloat VL = ql[PV->V][I];
        MFloat VR = qr[PV->V][I];

        const MFloat PL = ql[PV->P][I];
        const MFloat PR = qr[PV->P][I];

        const MFloat RHOL = ql[PV->RHO][I];
        const MFloat RHOR = qr[PV->RHO][I];

        MFloat RANSL[noRANS];
        MFloat RANSR[noRANS];
        for(MInt v = 0; v < noRANS; ++v) {
          RANSL[v] = ql[PV->RANS_VAR[v]][I];
          RANSR[v] = qr[PV->RANS_VAR[v]][I];
        }

        const MFloat surf0 = m_cells->surfaceMetrics[dim * 2 + 0][I];
        const MFloat surf1 = m_cells->surfaceMetrics[dim * 2 + 1][I];

        // compute lenght of metric vector for normalization
        const MFloat DGRAD = sqrt(POW2(surf0) + POW2(surf1));
        const MFloat FDGRAD = F1 / DGRAD;

        // scale by metric length to get velocity in the new basis (get normalized basis vectors)
        const MFloat UUL = ((UL * surf0 + VL * surf1)) * FDGRAD;


        const MFloat UUR = ((UR * surf0 + VR * surf1)) * FDGRAD;

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
        const MFloat e0 = (rans2eq_production_mode)
                              ? PLfRHOL * FgammaMinusOne + 0.5 * (POW2(UL) + POW2(VL)) + RANSL[0] + PLfRHOL
                              : PLfRHOL * FgammaMinusOne + 0.5 * (POW2(UL) + POW2(VL)) + PLfRHOL;
        const MFloat e1 = (rans2eq_production_mode)
                              ? PRfRHOR * FgammaMinusOne + 0.5 * (POW2(UR) + POW2(VR)) + RANSR[0] + PRfRHOR
                              : PRfRHOR * FgammaMinusOne + 0.5 * (POW2(UR) + POW2(VR)) + PRfRHOR;

        const MFloat RHOU = F1B2 * (MALR * (RHO_AL + RHO_AR) + fabs(MALR) * (RHO_AL - RHO_AR)) * DGRAD;
        const MFloat RHOU2 = F1B2 * RHOU;
        const MFloat ARHOU2 = fabs(RHOU2);

        flux[CV->RHO_U][I] = RHOU2 * (UL + UR) + ARHOU2 * (UL - UR) + PLR * surf0;
        flux[CV->RHO_V][I] = RHOU2 * (VL + VR) + ARHOU2 * (VL - VR) + PLR * surf1;
        flux[CV->RHO_E][I] = RHOU2 * (e0 + e1) + ARHOU2 * (e0 - e1);
        flux[CV->RHO][I] = RHOU;
        for(MInt v = 0; v < noRANS; ++v) {
          flux[CV->RANS_VAR[v]][I] = RHOU2 * (RANSL[v] + RANSR[v]) + ARHOU2 * (RANSL[v] - RANSR[v]);
        }
      }
    }

    // FLUX BALANCE
    for(MUint v = 0; v < noVars; v++) {
      for(MInt j = m_noGhostLayers; j < m_nCells[0] - m_noGhostLayers; j++) {
        for(MInt i = m_noGhostLayers; i < m_nCells[1] - m_noGhostLayers; i++) {
          const MInt I = cellIndex(i, j);
          const MInt IM1 = I - IJ[dim];
          const MUint offset = v * noCells;
          MFloat* const RESTRICT rhs = ALIGNED_F(cellRhs + offset);
          rhs[I] += flux[v][IM1] - flux[v][I];
        }
      }
    }
  }
}
template void FvStructuredSolver2DRans::Muscl_Ausm_Limited<6, 2>();
template void FvStructuredSolver2DRans::Muscl_Ausm_Limited<7, 2>();


void FvStructuredSolver2DRans::viscousFluxRANS() { (this->*viscFluxMethod)(); }


void FvStructuredSolver2DRans::viscousFlux_SA() {
  computeTurbViscosity_SA();

  // call the standard LES viscous flux
  m_solver->viscousFluxLES<>();

  // OTHER variables required to calculate the laminar viscous fluxes
  const MFloat rRe = F1 / m_solver->m_Re0;

  MFloat* const* const RESTRICT eflux = ALIGNED_F(m_cells->eFlux);
  MFloat* const* const RESTRICT fflux = ALIGNED_F(m_cells->fFlux);
  MFloat* const* const RESTRICT sa_1flux = ALIGNED_F(m_cells->saFlux1);
  MFloat* const* const RESTRICT sa_2flux = ALIGNED_F(m_cells->saFlux2);

  MFloat* const RESTRICT u = ALIGNED_F(m_cells->pvariables[PV->U]);
  MFloat* const RESTRICT v = ALIGNED_F(m_cells->pvariables[PV->V]);
  MFloat* const RESTRICT rho = ALIGNED_F(m_cells->pvariables[PV->RHO]);
  MFloat* const RESTRICT T = ALIGNED_F(m_cells->temperature);
  MFloat* const RESTRICT nuTilde = ALIGNED_F(m_cells->pvariables[PV->RANS_VAR[0]]);
  MFloat* const RESTRICT muLam = ALIGNED_F(m_cells->fq[FQ->MU_L]);

  for(MInt j = m_noGhostLayers - 1; j < m_nCells[0] - m_noGhostLayers + 1; j++) {
    for(MInt i = m_noGhostLayers - 1; i < m_nCells[1] - m_noGhostLayers + 1; i++) {
      // get the adjacent cells;

      const MInt IJ = cellIndex(i, j);
      const MInt IPJ = cellIndex((i + 1), j);
      const MInt IPJP = cellIndex((i + 1), (j + 1));
      const MInt IJP = cellIndex(i, (j + 1));
      const MInt IMJ = cellIndex((i - 1), j);
      const MInt IJM = cellIndex(i, (j - 1));

      const MFloat cornerMetrics[9] = {m_cells->cornerMetrics[0][IJ], m_cells->cornerMetrics[1][IJ],
                                       m_cells->cornerMetrics[2][IJ], m_cells->cornerMetrics[3][IJ]};

      const MFloat dnutdxi = F1B2 * (nuTilde[IPJP] + nuTilde[IPJ] - nuTilde[IJP] - nuTilde[IJ]);
      const MFloat dnutdet = F1B2 * (nuTilde[IPJP] + nuTilde[IJP] - nuTilde[IPJ] - nuTilde[IJ]);

      const MFloat nutldAvg = F1B4 * (nuTilde[IJP] + nuTilde[IPJP] + nuTilde[IJ] + nuTilde[IPJ]);

      const MFloat nuLamAvg =
          F1B4 * (muLam[IJP] / rho[IJP] + muLam[IPJP] / rho[IPJP] + muLam[IJ] / rho[IJ] + muLam[IPJ] / rho[IPJ]);


      const MFloat dnutldx =
          dnutdxi * m_cells->cornerMetrics[xsd * 2 + xsd][IJ] + dnutdet * m_cells->cornerMetrics[ysd * 2 + xsd][IJ];

      const MFloat dnutldy =
          dnutdxi * m_cells->cornerMetrics[xsd * 2 + ysd][IJ] + dnutdet * m_cells->cornerMetrics[ysd * 2 + ysd][IJ];


      const MFloat Frj = rRe / m_cells->cornerJac[IJ];

      const MFloat sax1 = Frj * (nuLamAvg + (1.0 + RM_SA_DV::cb2) * nutldAvg)
                          * (dnutldx * cornerMetrics[xsd * 2 + xsd] + dnutldy * cornerMetrics[xsd * 2 + ysd]);
      const MFloat sax2 =
          -Frj * RM_SA_DV::cb2 * (dnutldx * cornerMetrics[xsd * 2 + xsd] + dnutldy * cornerMetrics[xsd * 2 + ysd]);
      const MFloat say1 = Frj * (nuLamAvg + (1.0 + RM_SA_DV::cb2) * nutldAvg)
                          * (dnutldx * cornerMetrics[ysd * 2 + xsd] + dnutldy * cornerMetrics[ysd * 2 + ysd]);
      const MFloat say2 =
          -Frj * RM_SA_DV::cb2 * (dnutldx * cornerMetrics[ysd * 2 + xsd] + dnutldy * cornerMetrics[ysd * 2 + ysd]);


      // compute vorticity
      const MFloat du1 = u[IPJ] - u[IMJ];
      const MFloat du2 = u[IJP] - u[IJM];
      const MFloat dv1 = v[IPJ] - v[IMJ];
      const MFloat dv2 = v[IJP] - v[IJM];
      const MFloat vortk =
          (m_cells->cellMetrics[xsd * 2 + xsd][IJ] * dv1) + (m_cells->cellMetrics[ysd * 2 + xsd][IJ] * dv2)
          - (m_cells->cellMetrics[xsd * 2 + ysd][IJ] * du1) - (m_cells->cellMetrics[ysd * 2 + ysd][IJ] * du2);
      const MFloat s = F1B2 * fabs(vortk) / m_cells->cellJac[IJ];


      const MFloat distance = m_cells->fq[FQ->WALLDISTANCE][IJ];
      const MFloat Fdist2 = 1.0 / (distance * distance);
      const MFloat chi = nuTilde[IJ] * rho[IJ] / (SUTHERLANDLAW(T[IJ]));
      const MFloat chip3 = chi * chi * chi;
      const MFloat Fv1 = chip3 / (chip3 + RM_SA_DV::cv1to3);
      const MFloat Fv2 = F1 - (chi / (F1 + chi * Fv1));

      const MFloat term = nuTilde[IJ] * Fdist2 * RM_SA_DV::Fkap2;
      const MFloat stilde = s + term * Fv2 * rRe;
      const MFloat r = min(10.0, rRe * term / stilde);

      const MFloat g = r + RM_SA_DV::cw2 * (pow(r, 6) - r);
      const MFloat Fwterm = (1 + RM_SA_DV::cw3to6) / (pow(g, 6) + RM_SA_DV::cw3to6);
      const MFloat Fw = g * pow(Fwterm, (1.0 / 6.0));
      const MFloat prodValue = rho[IJ] * RM_SA_DV::cb1 * (F1 - RM_SA_DV::Ft2) * stilde * nuTilde[IJ];
      const MFloat destValue = rRe * rho[IJ] * (RM_SA_DV::cw1 * Fw - RM_SA_DV::cb1 * RM_SA_DV::Fkap2 * RM_SA_DV::Ft2)
                               * pow(nuTilde[IJ], 2.0) * Fdist2;

      m_cells->prodDest[IJ] = (prodValue - destValue) * m_cells->cellJac[IJ];

      eflux[0][IJ] = sax1; // diffusion of nutilde for every cell
      eflux[1][IJ] = sax2; // diffusion of nutilde for every cell

      fflux[0][IJ] = say1; // diffusion of nutilde for every cell
      fflux[1][IJ] = say2; // diffusion of nutilde for every cell
    }
  }

  ///////////////////////////////////////////
  //////////// SA1/SA2 //////////////////////
  ///////////////////////////////////////////
  for(MInt j = m_noGhostLayers; j < m_nCells[0] - m_noGhostLayers; ++j) {
    for(MInt i = m_noGhostLayers - 1; i < m_nCells[1] - m_noGhostLayers; ++i) {
      const MInt IJ = cellIndex(i, j);
      const MInt IJM = cellIndex(i, (j - 1));

      sa_1flux[0][IJ] = F1B2 * (eflux[0][IJ] + eflux[0][IJM]);
      sa_2flux[0][IJ] = F1B2 * (eflux[1][IJ] + eflux[1][IJM]);
    }
  }

  for(MInt j = m_noGhostLayers - 1; j < m_nCells[0] - m_noGhostLayers; ++j) {
    for(MInt i = m_noGhostLayers; i < m_nCells[1] - m_noGhostLayers; ++i) {
      const MInt IJ = cellIndex(i, j);
      const MInt IMJ = cellIndex((i - 1), j);

      sa_1flux[1][IJ] = F1B2 * (fflux[0][IJ] + fflux[0][IMJ]);
      sa_2flux[1][IJ] = F1B2 * (fflux[1][IJ] + fflux[1][IMJ]);
    }
  }


  // separate loop for adding the prodn nad destrn terms for tur kin viscosity transport variable
  for(MInt j = m_noGhostLayers; j < m_nCells[0] - m_noGhostLayers; j++) {
    for(MInt i = m_noGhostLayers; i < m_nCells[1] - m_noGhostLayers; i++) {
      const MInt IJ = cellIndex(i, j);
      const MInt IMJ = cellIndex((i - 1), j);
      const MInt IJM = cellIndex(i, (j - 1));

      const MFloat dissipation_term =
          (((sa_1flux[0][IJ] - sa_1flux[0][IMJ]) + ((sa_2flux[0][IJ] - sa_2flux[0][IMJ]) * nuTilde[IJ])) * rho[IJ]
               * RM_SA_DV::Fsigma
           + ((sa_1flux[1][IJ] - sa_1flux[1][IJM]) + ((sa_2flux[1][IJ] - sa_2flux[1][IJM]) * nuTilde[IJ])) * rho[IJ]
                 * RM_SA_DV::Fsigma);

      m_cells->rightHandSide[CV->RANS_VAR[0]][IJ] += dissipation_term;
      m_cells->rightHandSide[CV->RANS_VAR[0]][IJ] += m_cells->prodDest[IJ];
    }
  }
}

void FvStructuredSolver2DRans::viscousFlux_KEPSILON() {
  ASSERT(!m_solver->m_hasSingularity, "Get your ass out of this fucntion!!!");

  // compute friction velocity at wall
  m_solver->m_structuredBndryCnd->computeFrictionCoef();
  // communicate wall properties to all cells
  m_solver->m_structuredBndryCnd->distributeWallAndFPProperties();

  computeTurbViscosity_KEPSILON();

  // call the standard LES viscous flux
  if(m_solver->m_viscCompact)
    m_solver->viscousFluxLESCompact<true>();
  else
    m_solver->viscousFluxLES<true>();

  // OTHER variables required to calculate the laminar viscous fluxes
  static constexpr MFloat minMFloat =
      1e-16; // std::min(std::numeric_limits<MFloat>::min(), 1.0/std::numeric_limits<MFloat>::max());
  const MFloat rRe0 = F1 / m_solver->m_Re0;
  const MFloat gamma = m_solver->m_gamma;
  const MFloat fac = (m_solver->m_rans2eq_mode == "production") ? 1.0 : 0.0;
  const MFloat fac_nonDim = m_solver->m_keps_nonDimType ? 1.0 : PV->UInfinity;

  MFloat* const* const RESTRICT eflux = ALIGNED_F(m_cells->eFlux);
  MFloat* const* const RESTRICT fflux = ALIGNED_F(m_cells->fFlux);
  MFloat* const* const RESTRICT sa_1flux = ALIGNED_F(m_cells->saFlux1);
  MFloat* const* const RESTRICT sa_2flux = ALIGNED_F(m_cells->saFlux2);

  MFloat* const RESTRICT u = ALIGNED_F(m_cells->pvariables[PV->U]);
  MFloat* const RESTRICT v = ALIGNED_F(m_cells->pvariables[PV->V]);
  MFloat* const RESTRICT rho = ALIGNED_F(m_cells->pvariables[PV->RHO]);
  MFloat* const RESTRICT p = ALIGNED_F(m_cells->pvariables[PV->P]);
  //  MFloat* const RESTRICT T = ALIGNED_F(m_cells->temperature);
  MFloat* const RESTRICT TKE = ALIGNED_F(m_cells->pvariables[PV->RANS_VAR[0]]);
  MFloat* const RESTRICT eps = ALIGNED_F(m_cells->pvariables[PV->RANS_VAR[1]]);
  MFloat* const RESTRICT muLam = ALIGNED_F(m_cells->fq[FQ->MU_L]);
  MFloat* const RESTRICT muTurb = ALIGNED_F(m_cells->fq[FQ->MU_T]);
  const MFloat* const RESTRICT utau = &m_cells->fq[FQ->UTAU][0];
  const MFloat* const RESTRICT wallDist = &m_cells->fq[FQ->WALLDISTANCE][0];

  for(MInt j = m_noGhostLayers - 1; j < m_nCells[0] - m_noGhostLayers; j++) {
    for(MInt i = m_noGhostLayers - 1; i < m_nCells[1] - m_noGhostLayers; i++) {
      // get the adjacent cells;

      const MInt IJ = cellIndex(i, j);
      const MInt IPJ = cellIndex((i + 1), j);
      const MInt IPJP = cellIndex((i + 1), (j + 1));
      const MInt IJP = cellIndex(i, (j + 1));

      const MFloat cornerMetrics[nDim * nDim] = {m_cells->cornerMetrics[0][IJ], m_cells->cornerMetrics[1][IJ],
                                                 m_cells->cornerMetrics[2][IJ], m_cells->cornerMetrics[3][IJ]};

      // TODO_SS labels:FV if wall is along eta=constant, then at wall dTKEdxi and depsdxi must be zero! Check!!!

      const MFloat dTKEdxi = F1B2 * (TKE[IPJP] + TKE[IPJ] - TKE[IJP] - TKE[IJ]);
      const MFloat dTKEdet = F1B2 * (TKE[IPJP] + TKE[IJP] - TKE[IPJ] - TKE[IJ]);

      const MFloat dTKEdx = dTKEdxi * cornerMetrics[xsd * 2 + xsd] + dTKEdet * cornerMetrics[ysd * 2 + xsd];

      const MFloat dTKEdy = dTKEdxi * cornerMetrics[xsd * 2 + ysd] + dTKEdet * cornerMetrics[ysd * 2 + ysd];


      const MFloat depsdxi = F1B2 * (eps[IPJP] + eps[IPJ] - eps[IJP] - eps[IJ]);
      const MFloat depsdet = F1B2 * (eps[IPJP] + eps[IJP] - eps[IPJ] - eps[IJ]);

      const MFloat depsdx = depsdxi * cornerMetrics[xsd * 2 + xsd] + depsdet * cornerMetrics[ysd * 2 + xsd];

      const MFloat depsdy = depsdxi * cornerMetrics[xsd * 2 + ysd] + depsdet * cornerMetrics[ysd * 2 + ysd];

      // TODO_SS labels:FV muTurbAvg at wall surface must be zero! Check!
      const MFloat muLamAvg = F1B4 * (muLam[IJP] + muLam[IPJP] + muLam[IJ] + muLam[IPJ]);
      const MFloat muTurbAvg = F1B4 * (muTurb[IJP] + muTurb[IPJP] + muTurb[IJ] + muTurb[IPJ]);

      const MFloat Frj = rRe0 / m_cells->cornerJac[IJ];
      const MFloat mu_k = muLamAvg + muTurbAvg * RM_KEPS::rsigma_k;
      const MFloat mu_eps = muLamAvg + muTurbAvg * RM_KEPS::rsigma_eps;

      const MFloat sax1 = Frj * mu_k * (dTKEdx * cornerMetrics[xsd * 2 + xsd] + dTKEdy * cornerMetrics[xsd * 2 + ysd]);
      const MFloat sax2 =
          Frj * mu_eps * (depsdx * cornerMetrics[xsd * 2 + xsd] + depsdy * cornerMetrics[xsd * 2 + ysd]);
      const MFloat say1 = Frj * mu_k * (dTKEdx * cornerMetrics[ysd * 2 + xsd] + dTKEdy * cornerMetrics[ysd * 2 + ysd]);
      const MFloat say2 =
          Frj * mu_eps * (depsdx * cornerMetrics[ysd * 2 + xsd] + depsdy * cornerMetrics[ysd * 2 + ysd]);

      eflux[0][IJ] = sax1; // diffusion of nutilde for every cell
      eflux[1][IJ] = sax2; // diffusion of nutilde for every cell

      fflux[0][IJ] = say1; // diffusion of nutilde for every cell
      fflux[1][IJ] = say2; // diffusion of nutilde for every cell
    }
  }

  ///////////////////////////////////////////
  //////////// SA1/SA2 //////////////////////
  ///////////////////////////////////////////
  for(MInt j = m_noGhostLayers; j < m_nCells[0] - m_noGhostLayers; ++j) {
    for(MInt i = m_noGhostLayers - 1; i < m_nCells[1] - m_noGhostLayers; ++i) {
      const MInt IJ = cellIndex(i, j);
      const MInt IJM = cellIndex(i, (j - 1));

      sa_1flux[0][IJ] = F1B2 * (eflux[0][IJ] + eflux[0][IJM]);
      sa_2flux[0][IJ] = F1B2 * (eflux[1][IJ] + eflux[1][IJM]);
    }
  }

  for(MInt j = m_noGhostLayers - 1; j < m_nCells[0] - m_noGhostLayers; ++j) {
    for(MInt i = m_noGhostLayers; i < m_nCells[1] - m_noGhostLayers; ++i) {
      const MInt IJ = cellIndex(i, j);
      const MInt IMJ = cellIndex((i - 1), j);

      sa_1flux[1][IJ] = F1B2 * (fflux[0][IJ] + fflux[0][IMJ]);
      sa_2flux[1][IJ] = F1B2 * (fflux[1][IJ] + fflux[1][IMJ]);
    }
  }


  for(MInt j = m_noGhostLayers; j < m_nCells[0] - m_noGhostLayers; j++) {
    for(MInt i = m_noGhostLayers; i < m_nCells[1] - m_noGhostLayers; i++) {
      const MInt IJ = cellIndex(i, j);
      const MInt IMJ = cellIndex((i - 1), j);
      const MInt IJM = cellIndex(i, (j - 1));
      const MInt IPJ = cellIndex((i + 1), j);
      const MInt IJP = cellIndex(i, (j + 1));

      /////////////////////////////////////////////////////////////////////////
      /// Assemble diffusion terms from previously computed surface fluxes
      /////////////////////////////////////////////////////////////////////////
      const MFloat diffusion_TKE = (sa_1flux[0][IJ] - sa_1flux[0][IMJ]) + (sa_1flux[1][IJ] - sa_1flux[1][IJM]);

      const MFloat diffusion_eps = (sa_2flux[0][IJ] - sa_2flux[0][IMJ]) + (sa_2flux[1][IJ] - sa_2flux[1][IJM]);

      /////////////////////////////////////////////////////////////////////////
      ///  Production/Destruction terms + wall damping terms
      /////////////////////////////////////////////////////////////////////////
      const MFloat invCellJac = 1.0 / m_cells->cellJac[IJ];

      // compute
      const MFloat dudxi = 0.5 * (u[IPJ] - u[IMJ]); // TODO_SS labels:FV,totest check the 0.5; it should be correct
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

      // Production P = rho*nu_t/Re * ( 2*dudx^2 + 2*dvdy^2 + (dudy+dvdx)^2 - (2/3)*(dudx+dvdy)^2)
      //                - 2/3*rho*k * ( dudx + dvdy )
      const MFloat P = ((rRe0 / POW2(fac_nonDim)) * invCellJac * muTurb[IJ]
                            * (2.0 * POW2(dudx) + 2.0 * POW2(dvdy) + POW2(dudy + dvdx) - fac * F2B3 * POW2(dudx + dvdy))
                        - fac * F2B3 * rho[IJ] * TKE[IJ] * (dudx + dvdy))
                       * invCellJac;

      // Ret is only required for inner cells, for which eps & TKE must be positive
      const MFloat Ret =
          fac_nonDim * m_solver->m_Re0 * rho[IJ] * POW2(TKE[IJ]) / muLam[IJ] / std::max(eps[IJ], minMFloat);
      const MFloat f2 = 1 - 0.222222222222222 /*=0.4/1.8*/ * exp(-POW2(Ret) * 0.027777777777777777 /*=1/36*/);

      const MFloat Mt2 = 1.5 * TKE[IJ] * rho[IJ] / (gamma * p[IJ]);
      const MFloat Fdist2 = 1.0 / POW2(wallDist[IJ]);
      const MFloat k_diss1 = fac_nonDim * rho[IJ] * eps[IJ] * (1 + Mt2);
      const MFloat k_diss2 = 2.0 * rRe0 * muLam[IJ] * TKE[IJ] * Fdist2;

      const MFloat tau = m_cells->turbTimeScale[IJ]; // eps[IJ]/std::max(TKE[IJ], minMFloat); // time scale?!
      const MFloat eps_prod = RM_KEPS::C1 * tau * P;
      const MFloat yp = utau[IJ] * wallDist[IJ]; // utau[IJ]*wallDist[IJ]*rho[IJ]/muLam[IJ]; // TODO_SS labels:FV
      const MFloat exp_wall = exp(-RM_KEPS::C4 * m_solver->m_Re0 * yp);

      // limiter
      MFloat f2_ = f2;
      MFloat f3 = 1.0;
      if(m_solver->m_solutionAnomaly && m_cells->isAnomCell[IJ]) {
        const MFloat invRhoEps = 1.0 / (rho[IJ] * eps[IJ]);
        const MFloat RHS1 = invRhoEps / tau * diffusion_eps / m_cells->cellJac[IJ]
                            - invRhoEps * diffusion_TKE / m_cells->cellJac[IJ] + (RM_KEPS::C1 - 1.0) * P * invRhoEps
                            + fac_nonDim * (1 + Mt2);
        const MFloat f2_limit = (RHS1 + k_diss2 * invRhoEps * (1 - exp_wall)) / (fac_nonDim * RM_KEPS::C2);
        f2_ = std::min(1.0, std::max(f2, f2_limit));
        const MFloat diff = std::max(0.0, f2_limit - 1.0);
        f3 = std::max(0.0, 1.0 - fac_nonDim * RM_KEPS::C2 * diff / (k_diss2 * invRhoEps * (1 - exp_wall)));
      }

      const MFloat k_diss = k_diss1 + f3 * k_diss2;
      const MFloat eps_diss = tau * (fac_nonDim * RM_KEPS::C2 * f2_ * rho[IJ] * eps[IJ] + f3 * k_diss2 * exp_wall);

      const MFloat prodDest_TKE = (P - k_diss) * m_cells->cellJac[IJ];
      const MFloat prodDest_eps = (eps_prod - eps_diss) * m_cells->cellJac[IJ];

      /////////////////////////////////////////////////////////////////////////
      /// Adding up RHS
      /////////////////////////////////////////////////////////////////////////
      m_cells->rightHandSide[CV->RANS_VAR[0]][IJ] += diffusion_TKE;
      m_cells->rightHandSide[CV->RANS_VAR[0]][IJ] += prodDest_TKE;
      m_cells->rightHandSide[CV->RANS_VAR[1]][IJ] += diffusion_eps;
      m_cells->rightHandSide[CV->RANS_VAR[1]][IJ] += prodDest_eps;
    }
  }
}

void FvStructuredSolver2DRans::viscousFlux_KEPSILON2() {
  // compute friction velocity at wall
  m_solver->m_structuredBndryCnd->computeFrictionCoef();
  // communicate wall properties to all cells
  m_solver->m_structuredBndryCnd->distributeWallAndFPProperties();

  computeTurbViscosity_KEPSILON();

  // OTHER variables required to calculate the laminar viscous fluxes
  static constexpr MFloat minMFloat =
      1e-16; // std::min(std::numeric_limits<MFloat>::min(), 1.0/std::numeric_limits<MFloat>::max());
  const MFloat rRe0 = F1 / m_solver->m_Re0;
  const MFloat gamma = m_solver->m_gamma;
  const MFloat fac = (m_solver->m_rans2eq_mode == "production") ? 1.0 : 0.0;
  const MFloat fac_nonDim = m_solver->m_keps_nonDimType ? 1.0 : PV->UInfinity;
  const MFloat transPos = m_solver->m_ransTransPos;

  MFloat* const* const RESTRICT eflux = ALIGNED_F(m_cells->eFlux);
  MFloat* const* const RESTRICT fflux = ALIGNED_F(m_cells->fFlux);
  MFloat* const* const RESTRICT sa_1flux = ALIGNED_F(m_cells->saFlux1);
  MFloat* const* const RESTRICT sa_2flux = ALIGNED_F(m_cells->saFlux2);

  MFloat* const RESTRICT u = ALIGNED_F(m_cells->pvariables[PV->U]);
  MFloat* const RESTRICT v = ALIGNED_F(m_cells->pvariables[PV->V]);
  MFloat* const RESTRICT rho = ALIGNED_F(m_cells->pvariables[PV->RHO]);
  MFloat* const RESTRICT p = ALIGNED_F(m_cells->pvariables[PV->P]);
  //  MFloat* const RESTRICT T = ALIGNED_F(m_cells->temperature);
  MFloat* const RESTRICT TKE = ALIGNED_F(m_cells->pvariables[PV->RANS_VAR[0]]);
  MFloat* const RESTRICT eps = ALIGNED_F(m_cells->pvariables[PV->RANS_VAR[1]]);
  MFloat* const RESTRICT muLam = ALIGNED_F(m_cells->fq[FQ->MU_L]);
  MFloat* const RESTRICT muTurb = ALIGNED_F(m_cells->fq[FQ->MU_T]);
  const MFloat* const RESTRICT utau = &m_cells->fq[FQ->UTAU][0];
  const MFloat* const RESTRICT wallDist = &m_cells->fq[FQ->WALLDISTANCE][0];

  for(MInt j = m_noGhostLayers - 1; j < m_nCells[0] - m_noGhostLayers; j++) {
    for(MInt i = m_noGhostLayers - 1; i < m_nCells[1] - m_noGhostLayers; i++) {
      // get the adjacent cells;

      const MInt IJ = cellIndex(i, j);
      const MInt IPJ = cellIndex((i + 1), j);
      const MInt IPJP = cellIndex((i + 1), (j + 1));
      const MInt IJP = cellIndex(i, (j + 1));

      const MFloat cornerMetrics[nDim * nDim] = {m_cells->cornerMetrics[0][IJ], m_cells->cornerMetrics[1][IJ],
                                                 m_cells->cornerMetrics[2][IJ], m_cells->cornerMetrics[3][IJ]};

      // TODO_SS labels:FV,totest if wall is along eta=constant, then at wall dTKEdxi and depsdxi must be zero! Check!!!

      const MFloat dTKEdxi = F1B2 * (TKE[IPJP] + TKE[IPJ] - TKE[IJP] - TKE[IJ]);
      const MFloat dTKEdet = F1B2 * (TKE[IPJP] + TKE[IJP] - TKE[IPJ] - TKE[IJ]);

      const MFloat depsdxi = F1B2 * (eps[IPJP] + eps[IPJ] - eps[IJP] - eps[IJ]);
      const MFloat depsdet = F1B2 * (eps[IPJP] + eps[IJP] - eps[IPJ] - eps[IJ]);

      const MFloat metricTerms = cornerMetrics[xsd * 2 + xsd] * cornerMetrics[ysd * 2 + xsd]
                                 + cornerMetrics[xsd * 2 + ysd] * cornerMetrics[ysd * 2 + ysd];

      const MFloat invCornerJac = F1 / POW2(m_cells->cornerJac[IJ]);

      const MFloat sax1 = invCornerJac * dTKEdet * metricTerms;
      const MFloat sax2 = invCornerJac * depsdet * metricTerms;
      const MFloat say1 = invCornerJac * dTKEdxi * metricTerms;
      const MFloat say2 = invCornerJac * depsdxi * metricTerms;

      eflux[0][IJ] = sax1; // diffusion of nutilde for every cell
      eflux[1][IJ] = sax2; // diffusion of nutilde for every cell

      fflux[0][IJ] = say1; // diffusion of nutilde for every cell
      fflux[1][IJ] = say2; // diffusion of nutilde for every cell

      /////////////////////////////////////////////////////////////////////////
      /// Special production term treatment
      /////////////////////////////////////////////////////////////////////////
      if(m_P_keps) {
        const MFloat dudxi = F1B2 * (u[IPJP] + u[IPJ] - u[IJP] - u[IJ]);
        const MFloat dudet = F1B2 * (u[IPJP] + u[IJP] - u[IPJ] - u[IJ]);

        const MFloat dvdxi = F1B2 * (v[IPJP] + v[IPJ] - v[IJP] - v[IJ]);
        const MFloat dvdet = F1B2 * (v[IPJP] + v[IJP] - v[IPJ] - v[IJ]);

        const MFloat S11 = dudxi * cornerMetrics[xsd * 2 + xsd] + dudet * cornerMetrics[ysd * 2 + xsd];
        const MFloat S12 = 0.5
                           * (dudxi * cornerMetrics[xsd * 2 + ysd] + dudet * cornerMetrics[ysd * 2 + ysd]
                              + dvdxi * cornerMetrics[xsd * 2 + xsd] + dvdet * cornerMetrics[ysd * 2 + xsd]);
        const MFloat S22 = dvdxi * cornerMetrics[xsd * 2 + ysd] + dvdet * cornerMetrics[ysd * 2 + ysd];
        const MFloat S = S11 * S11 + 2 * S12 * S12 + S22 * S22;
        // TODO_SS labels:FV,totest muTurbAvg at wall surface should be zero! Check!
        const MFloat muTurbAvg = F1B4 * (muTurb[IJP] + muTurb[IPJP] + muTurb[IJ] + muTurb[IPJ]);

        //        if ((m_solver->domainId()==7 && IJ==83) || (m_solver->domainId()==8 && IJ==203))
        //          muTurbAvg = 0;

        const MFloat mueOverRe = rRe0 * invCornerJac * muTurbAvg / POW2(fac_nonDim);
        // The simplified (incompressible) version reads: 2*rho*nu_t/Re * S_ij*S_ij
        m_cells->P_keps[IJ] = 2.0 * mueOverRe * S;
      }
    }
  }


  // diffusive flux correction for the singular points
  // m_hasSingularity=0 means no singular points in this solver, otherwise do flux correction
  if(m_solver->m_hasSingularity > 0) {
    diffusiveFluxCorrection();
  }

  ///////////////////////////////////////////
  //////////// SA1/SA2 //////////////////////
  ///////////////////////////////////////////
  for(MInt j = m_noGhostLayers; j < m_nCells[0] - m_noGhostLayers; ++j) {
    for(MInt i = m_noGhostLayers - 1; i < m_nCells[1] - m_noGhostLayers; ++i) {
      const MInt IJ = cellIndex(i, j);
      const MInt IJM = cellIndex(i, (j - 1));
      const MInt IPJ = cellIndex(i + 1, j);
      const MFloat surf0 = m_cells->surfaceMetrics[0][IJ];
      const MFloat surf1 = m_cells->surfaceMetrics[1][IJ];

      const MFloat dTKEdxi = TKE[IPJ] - TKE[IJ];
      const MFloat depsdxi = eps[IPJ] - eps[IJ];

      const MFloat muLamAvg_xi = F1B2 * (muLam[IJ] + muLam[IPJ]);
      const MFloat muTurbAvg_xi = F1B2 * (muTurb[IJ] + muTurb[IPJ]);

      // TODO_SS labels:FV here surface jac (m_cells->surfJac) verwenden!!!
      const MFloat surfJac =
          0.5
          * (m_cells->cornerJac[IJ] + m_cells->cornerJac[IJM]); // TODO_SS labels:FV,totest Is this the right way to go?
      const MFloat Frj = rRe0 * surfJac;
      const MFloat mu_k_xi = muLamAvg_xi + muTurbAvg_xi * RM_KEPS::rsigma_k;
      const MFloat mu_eps_xi = muLamAvg_xi + muTurbAvg_xi * RM_KEPS::rsigma_eps;

      const MFloat temp_TKE = dTKEdxi * (POW2(surf0) + POW2(surf1)) * POW2(1.0 / surfJac);
      const MFloat temp_eps = depsdxi * (POW2(surf0) + POW2(surf1)) * POW2(1.0 / surfJac);

      MFloat limiterVisc1 = 1.0;
      MFloat limiterVisc2 = 1.0;
      if(m_solver->m_limiterVisc) {
        // TODO_SS labels:FV nochmal richtig anschauen
        const MFloat mu_ref = F1B2 * (m_cells->fq[FQ->LIMITERVISC][IJ] + m_cells->fq[FQ->LIMITERVISC][IPJ]);
        limiterVisc1 = std::min(1.0, mu_ref / std::abs(rRe0 * mu_k_xi));
        limiterVisc2 = std::min(1.0, mu_ref / std::abs(rRe0 * mu_eps_xi));
      }

      sa_1flux[0][IJ] = Frj * mu_k_xi * (temp_TKE + F1B2 * (eflux[0][IJ] + eflux[0][IJM])) * limiterVisc1;
      sa_2flux[0][IJ] = Frj * mu_eps_xi * (temp_eps + F1B2 * (eflux[1][IJ] + eflux[1][IJM])) * limiterVisc2;
    }
  }

  for(MInt j = m_noGhostLayers - 1; j < m_nCells[0] - m_noGhostLayers; ++j) {
    for(MInt i = m_noGhostLayers; i < m_nCells[1] - m_noGhostLayers; ++i) {
      const MInt IJ = cellIndex(i, j);
      const MInt IMJ = cellIndex((i - 1), j);
      const MInt IJP = cellIndex(i, j + 1);
      const MFloat surf0 = m_cells->surfaceMetrics[2 + 0][IJ];
      const MFloat surf1 = m_cells->surfaceMetrics[2 + 1][IJ];

      const MFloat dTKEdet = TKE[IJP] - TKE[IJ];
      const MFloat depsdet = eps[IJP] - eps[IJ];

      const MFloat muLamAvg_eta = F1B2 * (muLam[IJ] + muLam[IJP]);
      const MFloat muTurbAvg_eta = F1B2 * (muTurb[IJ] + muTurb[IJP]);

      //      const MFloat limit = 1e12;
      //      MFloat mu = muLamAvg_eta+muTurbAvg_eta;
      //      m_cells->fq[FQ->VAR1][IJ] = 1.0;
      //      if (m_cells->coordinates[0][IJ]>0.0) {
      //        const MInt IJM = cellIndex(i,j-1);
      //        const MFloat mu_ref =
      //        limit*POW2(0.5*(m_cells->coordinates[1][IJP]-m_cells->coordinates[1][IJM]))/m_solver->m_timeRef;
      //        m_cells->fq[FQ->VAR1][IJ] = std::min(mu_ref/mu, 1.0);
      //      }

      // TODO_SS labels:FV here surface jac (m_cells->surfJac) verwenden!!!
      const MFloat surfJac =
          0.5
          * (m_cells->cornerJac[IJ] + m_cells->cornerJac[IMJ]); // TODO_SS labels:FV,totest Is this the right way to go?
      const MFloat Frj = rRe0 * surfJac;
      const MFloat mu_k_eta = muLamAvg_eta + muTurbAvg_eta * RM_KEPS::rsigma_k;
      const MFloat mu_eps_eta = muLamAvg_eta + muTurbAvg_eta * RM_KEPS::rsigma_eps;

      const MFloat temp_TKE = dTKEdet * (POW2(surf0) + POW2(surf1)) * POW2(1.0 / surfJac);
      const MFloat temp_eps = depsdet * (POW2(surf0) + POW2(surf1)) * POW2(1.0 / surfJac);

      MFloat limiterVisc1 = 1.0;
      MFloat limiterVisc2 = 1.0;
      if(m_solver->m_limiterVisc) {
        // TODO_SS labels:FV nochmal richtig anschauen
        const MFloat mu_ref = F1B2 * (m_cells->fq[FQ->LIMITERVISC][IJ] + m_cells->fq[FQ->LIMITERVISC][IJP]);
        limiterVisc1 = std::min(1.0, mu_ref / std::abs(rRe0 * mu_k_eta));
        limiterVisc2 = std::min(1.0, mu_ref / std::abs(rRe0 * mu_eps_eta));
      }

      sa_1flux[1][IJ] = Frj * mu_k_eta * (temp_TKE + F1B2 * (fflux[0][IJ] + fflux[0][IMJ])) * limiterVisc1;
      sa_2flux[1][IJ] = Frj * mu_eps_eta * (temp_eps + F1B2 * (fflux[1][IJ] + fflux[1][IMJ])) * limiterVisc2;
    }
  }


  for(MInt j = m_noGhostLayers; j < m_nCells[0] - m_noGhostLayers; j++) {
    for(MInt i = m_noGhostLayers; i < m_nCells[1] - m_noGhostLayers; i++) {
      const MInt IJ = cellIndex(i, j);
      const MInt IMJ = cellIndex((i - 1), j);
      const MInt IJM = cellIndex(i, (j - 1));
      const MInt IPJ = cellIndex((i + 1), j);
      const MInt IJP = cellIndex(i, (j + 1));

      // Actualy skipping should not be necessary, but somehow for localTimeStepping, the code might react
      // to the non-meaningful values even though these are not used
      if(m_cells->coordinates[0][IJ] <= transPos) continue;

      /////////////////////////////////////////////////////////////////////////
      /// Assemble diffusion terms from previously computed surface fluxes
      /////////////////////////////////////////////////////////////////////////
      MFloat limiterVisc1 = 1.0;
      MFloat limiterVisc2 = 1.0;
      /*      if (m_solver->m_limiterVisc) {
              //TODO_SS labels:FV nochmal richtig anschauen
              const MFloat mu_ref = m_cells->fq[FQ->LIMITERVISC][IJ];
              limiterVisc1 = std::min(1.0, mu_ref/std::abs(rRe0*(muLam[IJ]+muTurb[IJ]*RM_KEPS::rsigma_k)));
              limiterVisc2 = std::min(1.0, mu_ref/std::abs(rRe0*(muLam[IJ]+muTurb[IJ]*RM_KEPS::rsigma_eps)));
            }*/

      const MFloat diffusion_TKE =
          (sa_1flux[0][IJ] - sa_1flux[0][IMJ]) + (sa_1flux[1][IJ] - sa_1flux[1][IJM]) * limiterVisc1;

      const MFloat diffusion_eps =
          (sa_2flux[0][IJ] - sa_2flux[0][IMJ]) + (sa_2flux[1][IJ] - sa_2flux[1][IJM]) * limiterVisc2;

      /////////////////////////////////////////////////////////////////////////
      ///  Production/Destruction terms + wall damping terms
      /////////////////////////////////////////////////////////////////////////
      const MFloat invCellJac = 1.0 / m_cells->cellJac[IJ];

      MFloat P;
      if(m_P_keps) {
        const MInt IMJM = cellIndex(i - 1, j - 1);
        P = F1B4 * (m_cells->P_keps[IJ] + m_cells->P_keps[IMJ] + m_cells->P_keps[IJM] + m_cells->P_keps[IMJM]);
      } else {
        // compute
        const MFloat dudxi = 0.5 * (u[IPJ] - u[IMJ]); // TODO_SS labels:FV,totest check the 0.5; it should be correct
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

        // Production P = rho*nu_t/Re * ( 2*dudx^2 + 2*dvdy^2 + (dudy+dvdx)^2 - (2/3)*(dudx+dvdy)^2)
        //                - 2/3*rho*k * ( dudx + dvdy )
        P = ((rRe0 / POW2(fac_nonDim)) * invCellJac * muTurb[IJ]
                 * (2.0 * POW2(dudx) + 2.0 * POW2(dvdy) + POW2(dudy + dvdx) - fac * F2B3 * POW2(dudx + dvdy))
             - fac * F2B3 * rho[IJ] * TKE[IJ] * (dudx + dvdy))
            * invCellJac;
      }

      // Ret is only required for inner cells, for which eps & TKE must be positive
      const MFloat Ret =
          fac_nonDim * m_solver->m_Re0 * rho[IJ] * POW2(TKE[IJ]) / muLam[IJ] / std::max(eps[IJ], minMFloat);
      const MFloat f2 = 1 - 0.222222222222222 /*=0.4/1.8*/ * exp(-POW2(Ret) * 0.027777777777777777 /*=1/36*/);

      const MFloat Mt2 = 1.5 * TKE[IJ] * rho[IJ] / (gamma * p[IJ]);
      const MFloat Fdist2 = 1.0 / POW2(wallDist[IJ]);
      const MFloat k_diss1 = fac_nonDim * rho[IJ] * eps[IJ] * (1 + Mt2);
      const MFloat k_diss2 = 2.0 * rRe0 * muLam[IJ] * TKE[IJ] * Fdist2;

      const MFloat tau = m_cells->turbTimeScale[IJ]; // eps[IJ]/std::max(TKE[IJ], minMFloat); // time scale?!
      const MFloat eps_prod = RM_KEPS::C1 * tau * P;
      const MFloat yp = utau[IJ] * wallDist[IJ]; // utau[IJ]*wallDist[IJ]*rho[IJ]/muLam[IJ]; // TODO_SS labels:FV
      const MFloat exp_wall = exp(-RM_KEPS::C4 * m_solver->m_Re0 * yp);

      // limiter
      MFloat f2_ = f2;
      MFloat f3 = 1.0;
      if(m_solver->m_solutionAnomaly && m_cells->isAnomCell[IJ]) {
        const MFloat invRhoEps = 1.0 / std::max(rho[IJ] * eps[IJ], minMFloat);
        // The additional porous terms which might be part of the rhs are omitted
        const MFloat RHS1 = invRhoEps / tau * diffusion_eps / m_cells->cellJac[IJ]
                            - invRhoEps * diffusion_TKE / m_cells->cellJac[IJ] + (RM_KEPS::C1 - 1.0) * P * invRhoEps
                            + fac_nonDim * (1 + Mt2);
        const MFloat f2_limit = (RHS1 + k_diss2 * invRhoEps * (1 - exp_wall)) / (fac_nonDim * RM_KEPS::C2);
        f2_ = std::min(1.0, std::max(f2, f2_limit));
        const MFloat diff = std::max(0.0, f2_limit - 1.0);
        f3 = std::max(0.0, 1.0 - fac_nonDim * RM_KEPS::C2 * diff / (k_diss2 * invRhoEps * (1 - exp_wall)));
      }

      const MFloat k_diss = k_diss1 + f3 * k_diss2;
      const MFloat eps_diss = tau * (fac_nonDim * RM_KEPS::C2 * f2_ * rho[IJ] * eps[IJ] + f3 * k_diss2 * exp_wall);

      const MFloat prodDest_TKE = (P - k_diss) * m_cells->cellJac[IJ];
      const MFloat prodDest_eps = (eps_prod - eps_diss) * m_cells->cellJac[IJ];

      /////////////////////////////////////////////////////////////////////////
      /// Adding up RHS
      /////////////////////////////////////////////////////////////////////////
      //      m_cells->fq[FQ->VAR1][IJ] = 2.0*POW2(dudx);//m_cells->rightHandSide[CV->RANS_VAR[0]][IJ];
      //      m_cells->fq[FQ->VAR2][IJ] = dudy;//diffusion_TKE;
      //      m_cells->fq[FQ->VAR3][IJ] = prodDest_TKE;
      m_cells->rightHandSide[CV->RANS_VAR[0]][IJ] += diffusion_TKE;
      m_cells->rightHandSide[CV->RANS_VAR[0]][IJ] += prodDest_TKE;
      m_cells->rightHandSide[CV->RANS_VAR[1]][IJ] += diffusion_eps;
      m_cells->rightHandSide[CV->RANS_VAR[1]][IJ] += prodDest_eps;
      //      m_cells->fq[FQ->VAR4][IJ] = dvdx;//m_cells->rightHandSide[CV->RANS_VAR[0]][IJ];
      //      m_cells->fq[FQ->VAR1][IJ] = rho[IJ]*v[IJ];
      //      m_cells->fq[FQ->VAR2][IJ] = v[IJ];
    }
  }

  // call the standard LES viscous flux
  if(m_solver->m_viscCompact)
    m_solver->viscousFluxLESCompact<true>();
  else
    m_solver->viscousFluxLES<true>();
}


void FvStructuredSolver2DRans::computeTurbViscosity() { (this->*compTurbVisc)(); }

void FvStructuredSolver2DRans::computeTurbViscosity_SA() {
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
      nuTilde[i] = mMin(mMax(nuTilde[i], 0.0), 3000.0);
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

void FvStructuredSolver2DRans::computeTurbViscosity_KEPSILON() {
  // nu_t is also required in the ghost cells at boundaries because of e.g. viscousFluxLES
  // OTHER variables required to calculate the laminar viscous fluxes
  static constexpr MFloat minMFloat =
      1e-16; // std::min(std::numeric_limits<MFloat>::min(), 1.0/std::numeric_limits<MFloat>::max());
  const MFloat transPos = m_solver->m_ransTransPos;
  const MFloat Re0 = m_solver->m_Re0;
  const MFloat fac_nonDim = m_solver->m_keps_nonDimType ? 1.0 : PV->UInfinity;
  static const MFloat sqrtF2B3 = sqrt(2.0) / 3.0;
  const MFloat c_wd_eff = m_solver->m_c_wd;
  MFloat* const RESTRICT p = &m_cells->pvariables[PV->P][0];
  MFloat* const RESTRICT rho = &m_cells->pvariables[PV->RHO][0];
  MFloat* const RESTRICT TKE = &m_cells->pvariables[PV->RANS_VAR[0]][0];
  MFloat* const RESTRICT eps = &m_cells->pvariables[PV->RANS_VAR[1]][0];
  MFloat* const RESTRICT T = &m_cells->temperature[0];
  MFloat* const RESTRICT u = &m_cells->pvariables[PV->U][0];
  MFloat* const RESTRICT v = &m_cells->pvariables[PV->V][0];
  const MFloat* const RESTRICT utau = &m_cells->fq[FQ->UTAU][0];
  const MFloat* const RESTRICT utau2 =
      (m_solver->m_blockType != "porous") ? &m_cells->fq[FQ->UTAU][0] : &m_cells->fq[FQ->UTAU2][0];
  const MFloat* const RESTRICT wallDist = &m_cells->fq[FQ->WALLDISTANCE][0];
  MFloat* const RESTRICT muLam = &m_cells->fq[FQ->MU_L][0];
  const MFloat* const RESTRICT por = &m_cells->fq[FQ->POROSITY][0];
  const MFloat* const RESTRICT Da = &m_cells->fq[FQ->DARCY][0];

  for(MInt jj = m_noGhostLayers - 1; jj < m_nCells[0] - m_noGhostLayers + 1; jj++) {
    for(MInt ii = m_noGhostLayers - 1; ii < m_nCells[1] - m_noGhostLayers + 1; ii++) {
      const MInt i = cellIndex(ii, jj);
      if(m_cells->coordinates[0][i] <= transPos) {
        TKE[i] = PV->ransInfinity[0];
        eps[i] = PV->ransInfinity[1];
      } else {
        // TODO_SS labels:FV maybe else not necessary
      }
      T[i] = m_solver->m_gamma * p[i] / rho[i];
      // decode the kinematic turbulent viscosity from the turb dynamic visc arrays
      muLam[i] = SUTHERLANDLAW(T[i]); // TODO_SS labels:FV put this maybe somewhere else
      //    const MFloat nuLaminar = muLam[i]/rho[i];
      const MFloat yp = utau[i] * wallDist[i]; /// nuLaminar;
      // TODO_SS labels:FV,totest Using Da & por in halo cells at fluid-porous interfaces might be not correct
      // TODO_SS labels:FV,toenhance In the following think if it makes sense to take the minmum of yp or wallDistance
      //    const MFloat yp_ = (m_solver->m_blockType=="porous" && c_wd_eff*sqrt(Da[i]/por[i])<wallDist[i]) ?
      //    utau2[i]*c_wd_eff*sqrt(Da[i]/por[i]) : yp;
      const MFloat yp_ = (m_solver->m_blockType == "porous" && utau2[i] * c_wd_eff * sqrt(Da[i] / por[i]) < yp)
                             ? utau2[i] * c_wd_eff * sqrt(Da[i] / por[i])
                             : yp;
      const MFloat f_mu = 1.0 - exp(-RM_KEPS::C3 * Re0 * yp_);
      const MFloat sgn = eps[i] < 0 ? -1.0 : 1.0;
      const MFloat eps_i = sgn * mMax(fabs(eps[i]), minMFloat); // eps & TKE can be negative in BC treatment

      //
      // Note: actually the computation of the derivatives here, needs to be in the same way as in the
      //      computation of the production term
      // TODO_SS labels:FV determine neighbor cells by i+i, etc.
      const MInt IPJ = getCellIdfromCell(i, 1, 0);
      const MInt IMJ = getCellIdfromCell(i, -1, 0);
      const MInt IJP = getCellIdfromCell(i, 0, 1);
      const MInt IJM = getCellIdfromCell(i, 0, -1);
      const MFloat dudxi = 0.5 * (u[IPJ] - u[IMJ]); // TODO_SS labels:FV,totest check the 0.5; it should be correct
      const MFloat dudet = 0.5 * (u[IJP] - u[IJM]);
      const MFloat dvdxi = 0.5 * (v[IPJ] - v[IMJ]);
      const MFloat dvdet = 0.5 * (v[IJP] - v[IJM]);

      const MFloat dudx =
          dudxi * m_cells->cellMetrics[xsd * 2 + xsd][i] + dudet * m_cells->cellMetrics[ysd * 2 + xsd][i];
      const MFloat dudy =
          dudxi * m_cells->cellMetrics[xsd * 2 + ysd][i] + dudet * m_cells->cellMetrics[ysd * 2 + ysd][i];
      const MFloat dvdx =
          dvdxi * m_cells->cellMetrics[xsd * 2 + xsd][i] + dvdet * m_cells->cellMetrics[ysd * 2 + xsd][i];
      const MFloat dvdy =
          dvdxi * m_cells->cellMetrics[xsd * 2 + ysd][i] + dvdet * m_cells->cellMetrics[ysd * 2 + ysd][i];

      const MFloat invCellJac = 1.0 / m_cells->cellJac[i];
      const MFloat S11 = F2B3 * dudx - F1B3 * dvdy;
      const MFloat S12 = 0.5 * (dudy + dvdx);
      const MFloat S22 = F2B3 * dvdy - F1B3 * dudx;
      const MFloat S = sqrt(S11 * S11 + 2 * S12 * S12 + S22 * S22) * invCellJac;
      const MFloat temp = sqrtF2B3 / std::max(S, minMFloat);
      const MFloat nuTurbMax = temp * fabs(TKE[i]) * Re0 * POW2(fac_nonDim);

      const MFloat nuTurb_true = fabs(Re0 * RM_KEPS::C_mu * f_mu * POW2(TKE[i]) / eps_i * fac_nonDim);
      //    m_cells->fq[FQ->VAR2][i] = nuTurbMax;
      //    if (nuTurbMax<nuTurb_true)
      //      m_cells->fq[FQ->VAR1][i] = true;
      //    else
      //      m_cells->fq[FQ->VAR1][i] = false;
      const MFloat nuTurb = sgn * std::min(nuTurbMax, nuTurb_true);
      m_cells->fq[FQ->NU_T][i] = nuTurb;
      m_cells->fq[FQ->MU_T][i] = rho[i] * nuTurb;
      // turbTimeScale is always positive, because, if TKE negative in boundary ghost cell, then also eps
      m_cells->turbTimeScale[i] =
          RM_KEPS::C_mu * f_mu * fabs(TKE[i]) * Re0 * fac_nonDim / std::max(fabs(nuTurb), minMFloat);

      //
      /*    MFloat tau1 = F4B3 * dudx - F2B3 * dvdy;
          MFloat tau2 = dudy + dvdx;
          const MFloat mueOverRe = 1.0/(Re0*m_cells->cellJac[i])*nuTurb;
          m_cells->fq[FQ->VAR3][i] = -mueOverRe*tau1+F2B3*TKE[i];
          m_cells->fq[FQ->VAR4][i] = -mueOverRe*tau2;*/
    }
  }
}


void FvStructuredSolver2DRans::setAndAllocate_KEPSILON() {
  mAlloc(m_cells->saFlux1, nDim, m_noCells, "m_cells->saFlux1", -999999.9, AT_);
  mAlloc(m_cells->saFlux2, nDim, m_noCells, "m_cells->saFlux2", -999999.9, AT_);
  mAlloc(m_cells->turbTimeScale, m_noCells, "m_cells->turbTimeScale", -999999.9, AT_);

  // Modified production term computation required temporary buffer
  m_P_keps = Context::getSolverProperty<MBool>("P_keps", m_solverId, AT_, &m_P_keps);
  if(m_P_keps) {
    if(m_solver->m_hasSingularity) mTerm(1, "P_keps=true for grid with singularities not implemented yet!");
    mAlloc(m_cells->P_keps, m_noCells, "m_cells->P_keps", -999999.9, AT_);
  }

  //
  m_solver->m_solutionAnomaly = false;
  m_solver->m_solutionAnomaly =
      Context::getSolverProperty<MBool>("solutionAnomaly", m_solverId, AT_, &m_solver->m_solutionAnomaly);
  if(m_solver->m_solutionAnomaly) {
    // TODO_SS labels:FV smoothly blend the anomaly area from the rest
    // By default if solution anomaly treatment is enabled all cells will be taken into account
    mAlloc(m_cells->isAnomCell, m_noCells, "m_cells->isAnomCell", true, AT_);
    if(Context::propertyExists("anomPoint", m_solverId) + Context::propertyExists("anomRadius") == 1)
      mTerm(1, "Either you specify 'anomPoint' and 'anomRadius' or neither of them!");
    if(Context::propertyExists("anomPoint", m_solverId)) {
      for(MInt i = 0; i < m_noCells; ++i) {
        m_cells->isAnomCell[i] = false;
      }
      const MInt noAnomPoints = Context::propertyLength("anomPoint", m_solverId);
      const MInt noAnomRadius = Context::propertyLength("anomRadius", m_solverId);
      if(noAnomPoints % nDim != 0 || noAnomPoints != nDim * noAnomRadius)
        mTerm(1, "IQ>10 mandatory for usage of k-epsilon model!!!");
      for(MInt n = 0; n < noAnomRadius; ++n) {
        const MFloat rPOW2 = POW2(Context::getSolverProperty<MFloat>("anomRadius", m_solverId, AT_, n));
        const MFloat c[nDim] = {Context::getSolverProperty<MFloat>("anomPoint", m_solverId, AT_, n * nDim),
                                Context::getSolverProperty<MFloat>("anomPoint", m_solverId, AT_, n * nDim + 1)};
        for(MInt i = 0; i < m_noCells; ++i) {
          if(POW2(m_cells->coordinates[0][i] - c[0]) + POW2(m_cells->coordinates[1][i] - c[1]) < rPOW2)
            m_cells->isAnomCell[i] = true;
        }
      }
    }
  } // if(m_solver->m_solutionAnomaly)
}


void FvStructuredSolver2DRans::diffusiveFluxCorrection() {
  MFloat* const* const RESTRICT eflux = ALIGNED_F(m_cells->eFlux);
  MFloat* const* const RESTRICT fflux = ALIGNED_F(m_cells->fFlux);
  MFloat* const RESTRICT TKE_ = ALIGNED_F(m_cells->pvariables[PV->RANS_VAR[0]]);
  MFloat* const RESTRICT EPS_ = ALIGNED_F(m_cells->pvariables[PV->RANS_VAR[1]]);

  const auto& singularity = m_solver->m_singularity;

  //  MInt dim = 0;
  MInt start[nDim], end[nDim], nghbr[10];
  MInt len1[nDim];
  //  MInt totalCells;

  for(MInt i = 0; i < m_solver->m_hasSingularity; ++i) {
    // only correct for bc 6000 not for bc 4000-5000
    if(singularity[i].BC == -6000) {
      //      totalCells=1;
      for(MInt j = 0; j < nDim; j++) {
        len1[j] = singularity[i].end[j] - singularity[i].start[j];
        //        if(len1[j]!=0)  totalCells*=len1[j];
      }

      for(MInt n = 0; n < nDim; ++n) {
        if(singularity[i].end[n] - singularity[i].start[n] > 1) {
          mTerm(1, "In 2D not possible!");
          //          dim=n;
          // start[n]=singularity[i].start[n]+1;
          start[n] = singularity[i].start[n] + 1;
          end[n] = singularity[i].end[n] - 1;
        } else {
          start[n] = singularity[i].start[n];
          end[n] = singularity[i].end[n];
        }
      }

      MFloat TKE[10], EPS[10];
      MFloat TKEcorner, EPScorner;
      for(MInt jj = start[1]; jj < end[1]; ++jj) {
        for(MInt ii = start[0]; ii < end[0]; ++ii) {
          MInt count = 0;
          //            MInt temp[nDim]{};

          const MInt IJ_ = cellIndex(ii + singularity[i].Viscous[0], jj + singularity[i].Viscous[1]);
          //            const MFloat* const surf = m_cells->surfaceMetrics[IJ_];
          // dxidx, dxidy, detadx, detady
          const MFloat surfaceMetricsS[nDim * nDim] = {
              m_cells->surfaceMetricsSingularity[0][i], m_cells->surfaceMetricsSingularity[1][i],
              m_cells->surfaceMetricsSingularity[2][i], m_cells->surfaceMetricsSingularity[3][i]};

          //            temp[dim]=1;
          const MInt IJ = cellIndex(ii, jj);
          nghbr[count++] = IJ;
          //            nghbr[count++]=cellIndex(ii+temp[0],jj+temp[1],kk+temp[2]);

          const MInt IPMJ_ = getCellIdfromCell(IJ, singularity[i].Viscous[0], 0);
          const MInt IJPM_ = getCellIdfromCell(IJ, 0, singularity[i].Viscous[1]);

          const MFloat surfaceMetrics[nDim * nDim] = {
              m_cells->surfaceMetrics[0][IPMJ_], m_cells->surfaceMetrics[1][IPMJ_], m_cells->surfaceMetrics[2][IJPM_],
              m_cells->surfaceMetrics[3][IJPM_]};

          for(MInt m = 0; m < singularity[i].Nstar - 1; ++m) {
            MInt* change = singularity[i].displacement[m];
            nghbr[count++] = cellIndex(ii + change[0], jj + change[1]);
            //              nghbr[count++]=cellIndex(ii+temp[0]+change[0],jj+temp[1]+change[1],kk+temp[2]+change[2]);
          }

          if(count != singularity[i].Nstar) {
            cout << "what the hell! it is wrong!!!" << endl;
          }

          for(MInt m = 0; m < singularity[i].Nstar; ++m) {
            TKE[m] = TKE_[nghbr[m]];
            EPS[m] = EPS_[nghbr[m]];
          }

          TKEcorner = F0;
          EPScorner = F0;

          MInt id2 = ii - start[0] + (jj - start[1]) * len1[0];

          for(MInt n = 0; n < count; n++) {
            MInt ID = id2 * count + n;
            TKEcorner += singularity[i].ReconstructionConstants[0][ID] * TKE[n];
            EPScorner += singularity[i].ReconstructionConstants[0][ID] * EPS[n];
          }

          const MInt sign_xi = 2 * singularity[i].Viscous[0] + 1;
          const MInt sign_eta = 2 * singularity[i].Viscous[1] + 1;

          const MInt IPMJ = getCellIdfromCell(IJ, sign_xi, 0);
          const MInt IJPM = getCellIdfromCell(IJ, 0, sign_eta);

          const MFloat dTKEdet = sign_eta * 0.5 * (TKEcorner - 0.5 * (TKE[0] + TKE_[IPMJ]));
          const MFloat dEPSdet = sign_eta * 0.5 * (EPScorner - 0.5 * (EPS[0] + EPS_[IPMJ]));

          const MFloat dTKEdxi = sign_xi * 0.5 * (TKEcorner - 0.5 * (TKE[0] + TKE_[IJPM]));
          const MFloat dEPSdxi = sign_xi * 0.5 * (EPScorner - 0.5 * (EPS[0] + EPS_[IJPM]));

          //            const MFloat metricTerms = cornerMetrics[ xsd * 2 + xsd ]*cornerMetrics[ ysd * 2 + xsd ]
          //                                       + cornerMetrics[ xsd * 2 + ysd ]*cornerMetrics[ ysd * 2 + ysd ];

          const MFloat metricTerms_xi = surfaceMetrics[xsd * 2 + xsd] * surfaceMetricsS[ysd * 2 + xsd]
                                        + surfaceMetrics[xsd * 2 + ysd] * surfaceMetricsS[ysd * 2 + ysd];
          const MFloat metricTerms_eta = surfaceMetricsS[xsd * 2 + xsd] * surfaceMetrics[ysd * 2 + xsd]
                                         + surfaceMetricsS[xsd * 2 + ysd] * surfaceMetrics[ysd * 2 + ysd];

          const MFloat invSurfJac_xi = F1 / POW2(m_cells->surfJac[IPMJ_]);
          const MFloat invSurfJac_eta = F1 / POW2(m_cells->surfJac[m_noCells + IJPM_]);

          const MFloat sax1 = invSurfJac_xi * dTKEdet * metricTerms_xi;
          const MFloat sax2 = invSurfJac_xi * dEPSdet * metricTerms_xi;
          const MFloat say1 = invSurfJac_eta * dTKEdxi * metricTerms_eta;
          const MFloat say2 = invSurfJac_eta * dEPSdxi * metricTerms_eta;

          eflux[0][IJ_] = sax1; // diffusion of nutilde for every cell
          eflux[1][IJ_] = sax2; // diffusion of nutilde for every cell

          fflux[0][IJ_] = say1; // diffusion of nutilde for every cell
          fflux[1][IJ_] = say2; // diffusion of nutilde for every cell

          //            cout << globalTimeStep << "(" << m_solver->m_RKStep << ") dom=" << m_solver->domainId() <<
          //            setprecision(10) << " x|y=" << m_cells->coordinates[0][IJ] << "|" << m_cells->coordinates[1][IJ]
          //            << " eflux=" << eflux[ 0*noCells+IJ_ ] << "|" << eflux[ 1*noCells+IJ_ ] << " fflux=" << fflux[
          //            0*noCells+IJ_ ] << "|" << fflux[ 1*noCells+IJ_ ] << " metricTerms_xi=" << metricTerms_xi << "
          //            metricTerms_eta=" << metricTerms_eta << " " << m_cells->surfaceMetricsSingularity[i][0] << "|"
          //            << m_cells->surfaceMetricsSingularity[i][1] << "|" << m_cells->surfaceMetricsSingularity[i][2]
          //            << "|" << m_cells->surfaceMetricsSingularity[i][3]
          //              << " " << surfaceMetrics[0] << "|" << surfaceMetrics[1] << "|" << surfaceMetrics[2] << " " <<
          //              surfaceMetrics[3] << endl;
        }
      }
    }
  }
}


inline MInt FvStructuredSolver2DRans::cellIndex(MInt i, MInt j) { return i + j * m_nCells[1]; }

inline MInt FvStructuredSolver2DRans::getCellIdfromCell(MInt origin, MInt incI, MInt incJ) {
  return origin + incI + incJ * m_nCells[1];
}

inline MFloat FvStructuredSolver2DRans::getPSI(MInt I, MInt dim) {
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
      mMin(F1,
           FK
               * mMax(mMax(fabs((PIM2 - PIM1) / mMin(PIM2, PIM1)), fabs((PIM1 - PIP1) / mMin(PIM1, PIP1))),
                      fabs((PIP1 - PIP2) / mMin(PIP1, PIP2))));
  return PSI;
}
