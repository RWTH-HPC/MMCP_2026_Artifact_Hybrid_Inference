// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "lsfvcombustion.h"

#include <algorithm>
#include <stack>
#include "IO/parallelio.h"
#include "MEMORY/alloc.h"
#include "UTIL/functions.h"
#include "COMM/mpioverride.h"
#include "globals.h"
#include "UTIL/kdtree.h"
#include "coupling.h"
#include "globalvariables.h"

using namespace std;

template <MInt nDim, class SysEqn>
LsFvCombustion<nDim, SysEqn>::LsFvCombustion(const MInt couplingId, LsSolver* ls, FvCartesianSolver* fv)
  : Coupling(couplingId),
    // CouplingLS<nDim>(couplingId, ls),
    // CouplingFv<nDim, SysEqn>(couplingId, dynamic_cast<Solver*>(fv)),
    m_lsSolver(ls),
    m_fvSolver(fv) {
  TRACE();
  // fvSolver().setCombustionCouplingClassPointer(this);
  // lsSolver().setCombustionCouplingClassPointer(this);
  m_maxNoSets = lsSolver().m_maxNoSets;
}

template <MInt nDim, class SysEqn>
void LsFvCombustion<nDim, SysEqn>::init() {
  mAlloc(fvSolver().m_levelSetValues, lsSolver().m_noSets, "fvSolver().m_levelSetValues", AT_);
  mAlloc(fvSolver().m_curvatureG, lsSolver().m_noSets, "fvSolver().m_curvatureG", AT_);
  mAlloc(fvSolver().m_flameSpeedG, lsSolver().m_noSets, "fvSolver().m_flameSpeedG", AT_);
  computeGCellTimeStep();
}
template <MInt nDim, class SysEqn>
void LsFvCombustion<nDim, SysEqn>::preCouple(MInt) {
  computeGCellTimeStep();
  exchangeCouplingData();
}
template <MInt nDim, class SysEqn>
void LsFvCombustion<nDim, SysEqn>::postCouple(MInt) {}

template <MInt nDim, class SysEqn>
void LsFvCombustion<nDim, SysEqn>::finalizeCouplerInit() {
  for(MInt s = 0; s < lsSolver().m_noSets; s++) {
    fvSolver().m_levelSetValues[s].resize(fvSolver().a_noCells());
    fvSolver().m_curvatureG[s].resize(fvSolver().a_noCells());
    fvSolver().m_flameSpeedG[s].resize(fvSolver().a_noCells());
  }
  computeGCellTimeStep();
  setRhoFlameTubeInLs();
  setRhoInfinityInLs();
}

template <MInt nDim, class SysEqn>
void LsFvCombustion<nDim, SysEqn>::constructExtensionVelocity() {
  if(lsSolver().m_gRKStep != 0) return;
  if(!lsSolver().m_computeExtVel) return;

  if(!lsSolver().m_smoothExtVel) {
    fastInterfaceExtensionVelocity();
    return;
  }
  switch(lsSolver().m_computeExtVel) {
    case 700: {
      for(MInt set = lsSolver().m_startSet; set < lsSolver().m_noSets; set++) {
        MFloatScratchSpace fluidDensity(lsSolver().a_noG0Cells(set), AT_, "fluidDensity");
        fluidDensity.fill(-F1);
        if(!lsSolver().m_interpolateFlowFieldToFlameFront) {
          collectGEquationModelDataOpt(fluidDensity.getPointer(), set);
        } else {
          collectGEquationModelDataOptInterpolate(fluidDensity.getPointer(), set);
        }
        lsSolver().computeExtensionVelocityGEQUPVMarksteinOpt(fluidDensity.getPointer(), set);
      }
      break;
    }
    default: {
      break;
    }
  }
}

template <MInt nDim, class SysEqn>
void LsFvCombustion<nDim, SysEqn>::collectGEquationModelDataOpt(MFloat* fluidDensity, MInt set) {
  TRACE();

  MInt baseGCellId, cellId;
  for(MInt id = 0; id < lsSolver().a_noBandCells(set); id++) {
    cellId = lsSolver().a_bandCellId(id, set);
    for(MInt i = 0; i < lsSolver().nDim; i++) {
      lsSolver().a_extensionVelocityG(cellId, i, set) = F0;
    }
    lsSolver().a_correctedBurningVelocity(cellId, set) = F0;
  }
  //}

  // take only non-ghost cells!
  for(MInt id = 0; id < lsSolver().a_noG0Cells(set); id++) {
    if(lsSolver().a_isHalo(lsSolver().a_G0CellId(id, set))
       || lsSolver().a_level(lsSolver().a_G0CellId(id, set)) != lsSolver().a_maxGCellLevel()) {
      for(MInt i = 0; i < lsSolver().nDim; i++)
        lsSolver().a_extensionVelocityG(lsSolver().a_G0CellId(id, set), i, set) = F0;
      fluidDensity[id] = -F1;

    } else {
      // find the parent cell which is connected to the flow grid
      baseGCellId = lsSolver().a_G0CellId(id, set);

      if(lsSolver().a_isGBoundaryCellG(baseGCellId, 0)) continue;
      if(lsSolver().a_isBndryCellG(baseGCellId)) continue;

      while(ls2fvId(baseGCellId) == -1) {
        baseGCellId = lsSolver().c_parentId(baseGCellId);
      }
      if(baseGCellId == -1) {
        cerr << "ERROR: no parent cell found for connecting" << endl;
      }
      MInt flowCell = baseGCellId;

      fluidDensity[id] = collectFvDataForCollectGEquationModelDataOpt(flowCell, id);
    }
  }
}


template <MInt nDim, class SysEqn>
void LsFvCombustion<nDim, SysEqn>::exchangeCouplingData() {
  // set levelset function in the fv solver
  for(MInt c = 0; c < fvSolver().m_bndryGhostCellsOffset; c++) {
    if(fvSolver().grid().tree().hasProperty(c, Cell::IsHalo)) {
      continue;
    }
    if(!fvSolver().a_hasProperty(c, SolverCell::IsOnCurrentMGLevel)) {
      continue;
    }
    MInt cLs = fv2lsId(c);
    // TODO labels:COUPLER,FV,LS,toenhance some cells need to be skipped in the FV solver.
    //                              uses a very small dummy value for now.
    if(cLs < 0) {
      fvSolver().a_levelSetFunction(c, 0) = -999999;
      continue;
    }
    MInt gc = cLs;

    if(ABS(lsSolver().m_outsideGValue - ABS(lsSolver().a_levelSetFunctionG(gc, 0))) < 0.02) {
      fvSolver().a_levelSetFunction(c, 0) = -999999;
      continue;
    }
    if(lsSolver().a_level(gc) != fvSolver().maxRefinementLevel()) {
      fvSolver().a_levelSetFunction(c, 0) = -999999;
      continue;
    }

    fvSolver().a_levelSetFunction(c, 0) = lsSolver().a_levelSetFunctionG(gc, 0);
    fvSolver().a_flameSpeed(c, 0) = lsSolver().a_flameSpeedG(gc, 0);
    fvSolver().a_curvatureG(c, 0) = lsSolver().a_curvatureG(gc, 0);
  }

  // set flow field data in the LS solver
  constructExtensionVelocity();
}

template <MInt nDim, class SysEqn>
void LsFvCombustion<nDim, SysEqn>::computeGCellTimeStep() {
  lsSolver().m_timeStep = fvSolver().timeStep(true);
}

template <MInt nDim, class SysEqn>
void LsFvCombustion<nDim, SysEqn>::setLsTimeStep(MFloat timeStep) {
  lsSolver().m_timeStep = timeStep;
}


template <MInt nDim, class SysEqn>
MFloat LsFvCombustion<nDim, SysEqn>::collectFvDataForCollectGEquationModelDataOpt(MInt flowCell, MInt id) {
  MInt flowCellFv = ls2fvId(flowCell);

  if(flowCellFv == -1) {
    return -99;
  } else {
    for(MInt i = 0; i < nDim; i++) {
      lsSolver().a_extensionVelocityG(lsSolver().a_G0CellId(id, 0), i, 0) =
          fvSolver().a_pvariable(flowCellFv, fvSolver().PV->VV[i]);
    }

    return F1 / fvSolver().a_pvariable(flowCellFv, fvSolver().PV->RHO);
  }
}

template <MInt nDim, class SysEqn>
void LsFvCombustion<nDim, SysEqn>::setRhoFlameTubeInLs() {
  lsSolver().m_rhoFlameTube = fvSolver().m_rhoFlameTube;
}
template <MInt nDim, class SysEqn>
void LsFvCombustion<nDim, SysEqn>::setRhoInfinityInLs() {
  lsSolver().m_rhoInfinity = fvSolver().m_rhoInfinity;
}


template <MInt nDim, class SysEqn>
MFloat LsFvCombustion<nDim, SysEqn>::collectGFromCouplingClass(MInt i) {
  if(i + 1 > fvSolver().grid().tree().size()) {
    return -99;
  }
  if(fvSolver().grid().tree().hasProperty(i, Cell::IsHalo)) {
    return -99;
  }
  MInt iLs = fv2lsId(i);
  // sohel: for some reason iLs can be a big number in parallel runs.. maybe value is broken for halo cells?
  if(iLs < 0 || iLs > lsSolver().grid().tree().size()) {
    return -99;
  } else {
    return (lsSolver().a_levelSetFunctionG(iLs, 0)); // numeric_limits<MFloat>::infinity() ;
  }
}
template <MInt nDim, class SysEqn>
MFloat LsFvCombustion<nDim, SysEqn>::collectCurvFromCouplingClass(MInt i) {
  if(i + 1 > fvSolver().grid().tree().size()) {
    return -99;
  }
  if(fvSolver().grid().tree().hasProperty(i, Cell::IsHalo)) {
    return -99;
  }
  MInt iLs = fv2lsId(i);
  if(iLs < 0 || iLs > lsSolver().grid().tree().size()) {
    return -99;
  } else {
    return (lsSolver().a_curvatureG(iLs, 0)); // numeric_limits<MFloat>::infinity() ;
  }
}
template <MInt nDim, class SysEqn>
MInt LsFvCombustion<nDim, SysEqn>::noLevelSetFieldData() {
  MInt noLevelSetFieldData = 0;
  if(fvSolver().m_levelSet) {
    if(lsSolver().m_writeOutAllLevelSetFunctions) {
      noLevelSetFieldData += lsSolver().m_noSets;
      noLevelSetFieldData += lsSolver().m_noSets;
    } else
      noLevelSetFieldData += 1;
    if(lsSolver().m_writeOutAllCurvatures)
      noLevelSetFieldData += lsSolver().m_noSets;
    else
      noLevelSetFieldData += 1;
  }

  return noLevelSetFieldData;
}


template <MInt nDim, class SysEqn>
void LsFvCombustion<nDim, SysEqn>::saveOutputLS() {
  if(fvSolver().m_combustion && fvSolver().m_recordPressure && globalTimeStep % 10 == 0) {
    lsSolver().computeZeroLevelSetArcLength();
    if(fvSolver().domainId() == 0) {
      FILE* datei;
      datei = fopen("pressureSensor", "a+");
      fprintf(datei, " %d", globalTimeStep);
      fprintf(datei, " %f", fvSolver().m_time);
      fprintf(datei, " %-10.10f", fvSolver().a_pvariable(fvSolver().m_cellToRecordData, fvSolver().PV->P));
      fprintf(datei, " %-10.10f", fvSolver().m_meanPressure);
      fprintf(datei, " %-10.10f", fvSolver().a_pvariable(fvSolver().m_cellToRecordData, fvSolver().PV->V));
      fprintf(datei, " %-10.10f", lsSolver().m_arcLength);
      fprintf(datei, " %-10.10f", fvSolver().m_totalHeatReleaseRate);
      fprintf(datei, " %-10.10f", lsSolver().m_minFlameFrontPosition[1]);
      fprintf(datei, " %-10.10f", lsSolver().m_maxFlameFrontPosition[1]);
      fprintf(datei, " %-10.10f", lsSolver().m_meanFlameFrontPosition[1]);
      fprintf(datei, "\n");
      fclose(datei);
    }
  }

  if(fvSolver().m_combustion && fvSolver().m_recordFlameFrontPosition && globalTimeStep % 10 == 0) {
    lsSolver().computeZeroLevelSetArcLength();
    if(fvSolver().domainId() == 0) {
      FILE* datei;
      datei = fopen("flameFrontData", "a+");
      fprintf(datei, " %d", globalTimeStep);
      fprintf(datei, " %f", fvSolver().m_time);
      fprintf(datei, " %f", lsSolver().m_arcLength);
      fprintf(datei, " %-10.15f", lsSolver().m_minFlameFrontPosition[1]);
      fprintf(datei, " %-10.15f", lsSolver().m_maxFlameFrontPosition[1]);
      fprintf(datei, "\n");
      fclose(datei);
    }
  }

  if(fvSolver().m_combustion && !fvSolver().m_forcing && !fvSolver().m_structuredFlameOutput
     && fvSolver().domainId() == 0) {
    FILE* datei00;
    datei00 = fopen("flameSurfaceAreaROH_steady", "a+");
    fprintf(datei00, " %d", globalTimeStep);
    fprintf(datei00, " %f", fvSolver().m_time);
    fprintf(datei00, " %-10.10f", lsSolver().m_massConsumption);
    fprintf(datei00, " %-10.10f", fvSolver().m_totalHeatReleaseRate);
    fprintf(datei00, " %-10.10f", lsSolver().m_arcLength);
    fprintf(datei00, "\n");
    fclose(datei00);
  }


  if(fvSolver().m_combustion && (fvSolver().m_structuredFlameOutput || fvSolver().m_structuredFlameOutputLevel == 7)) {
    MInt sweepStart = 0;
    MInt sweepStartRow = 0;
    MInt current;
    MInt sample, currentCycle;
    MFloat forcingAmplitude;
    MInt samplingStartCycle = fvSolver().m_samplingStartCycle;
    MInt samplingEndCycle = fvSolver().m_samplingEndCycle;
    MFloat St = fvSolver().m_flameStrouhal;
    //---

    // compute the cycle time and the cycle
    sample = (MInt)(globalTimeStep / fvSolver().m_noTimeStepsBetweenSamples);

    currentCycle = (MInt)(sample / fvSolver().m_samplesPerCycle);

    // compute the flame surface area (unfiltered signal)
    lsSolver().computeZeroLevelSetArcLength();
    forcingAmplitude = fvSolver().m_forcingAmplitude * sin(St * fvSolver().m_time);

    if(fvSolver().domainId() == 0) {
      FILE* datei0;
      datei0 = fopen("flameSurfaceAreaROH", "a+");
      fprintf(datei0, " %d", globalTimeStep);
      fprintf(datei0, " %d", currentCycle);
      fprintf(datei0, " %d", sample);
      fprintf(datei0, " %f", fvSolver().m_time);
      fprintf(datei0, " %f", forcingAmplitude);
      fprintf(datei0, " %-10.10f", lsSolver().m_massConsumption);
      fprintf(datei0, " %-10.10f", fvSolver().m_totalHeatReleaseRate);
      fprintf(datei0, " %-10.10f", lsSolver().m_arcLength);
      fprintf(datei0, "\n");
      fclose(datei0);
    }

    // really?...
    switch(fvSolver().m_structuredFlameOutputLevel) {
      default: {
        if(globalTimeStep % fvSolver().m_noTimeStepsBetweenSamples != 0) {
          return;
        }
        break;
      }
    }


    if(fvSolver().domainId() == 0) {
      // output heat release over velocity forcing (filtered signal)
      FILE* datei;
      datei = fopen("flameSurfaceArea", "a+");
      fprintf(datei, " %d", globalTimeStep);
      fprintf(datei, " %d", currentCycle);
      fprintf(datei, " %d", sample);
      fprintf(datei, " %f", fvSolver().m_time);
      fprintf(datei, " %-10.10f", forcingAmplitude);
      fprintf(datei, " %-10.10f", lsSolver().m_massConsumption);
      fprintf(datei, " %-10.10f", fvSolver().m_totalHeatReleaseRate);
      fprintf(datei, " %-10.10f", lsSolver().m_arcLength);
      fprintf(datei, "\n");
      fclose(datei);
    }

    // really?...
    switch(fvSolver().m_structuredFlameOutputLevel) {
      default: {
        if(currentCycle < samplingStartCycle || currentCycle > samplingEndCycle) {
          return;
        }
        break;
      }
    }


    /** \brief Modified structured grid output (ASCII) of the calculated variables rho,c,v,u,... for the forced response
     * cases
     *
     * \author Stephan Schlimpert
     * \date 15.11.2010, June 2011
     *
     *  Structured Output written for a single flame (case 1) or a flame with a plenum above.
     *  Used in Matlab in order to compute the transfer functions of the flames via the velocity perturbations of the
     * centerline and the flame surface perturbations computed, see flameSurfaceArea.
     */
    switch(fvSolver().m_structuredFlameOutputLevel) {
        ///* flame with tube above (new version) + pd/dt information, correctded flame front amplitude computation
      case 400: {
        MFloat FtimeStep = 1 / fvSolver().timeStep();
        MFloat dPdT = -9999999.0;
        const MFloat gammaMinusOne = fvSolver().m_gamma - 1.0;
        MFloat fRho = -9999999.0;
        MFloat vel = -9999999.0;
        MFloat velPOW2 = -9999999.0;
        MFloat pold = -9999999.0;

        /** \brief
         * find cell on r_min in the left lower corner of the flame tube (= inlet of tube from plenum), only in the
         * first iteration order of output: rows: x; columns: y
         *
         */

        if((MInt)sample
               <= (MInt)(samplingStartCycle * fvSolver().m_samplesPerCycle + fvSolver().m_noTimeStepsBetweenSamples)
           || globalTimeStep <= (MInt)(fvSolver().m_restartTimeStep + fvSolver().m_noTimeStepsBetweenSamples)) {
          // for plenum simulation
          if(fvSolver().m_plenum) {
            for(MInt ac = 0; ac < fvSolver().m_noActiveCells; ac++) {
              if(fvSolver().a_coordinate(fvSolver().m_activeCellIds[ac], 0) > -(fvSolver().m_radiusFlameTube + 0.1))
                continue;
              if(fvSolver().a_coordinate(fvSolver().m_activeCellIds[ac], 1) > -1.1) continue;
              if(fvSolver().a_hasNeighbor(fvSolver().m_activeCellIds[ac], 0) > 0
                 && fvSolver().a_hasNeighbor(fvSolver().m_activeCellIds[ac], 2)
                        > 0) // number of neighbors in -x and -y direction > 0
                continue;
              fvSolver().m_sweepStartFirstCell =
                  fvSolver().m_activeCellIds[ac]; /// Structured Output starts from this founded Cell in every timestep
              break;
            }
            // only flame version
          } else {
            for(MInt ac = 0; ac < fvSolver().m_noActiveCells; ac++) {
              if(fvSolver().a_hasNeighbor(fvSolver().m_activeCellIds[ac], 0) > 0) continue;
              if(fvSolver().a_hasNeighbor(fvSolver().m_activeCellIds[ac], 2) > 0) continue;
              if(fvSolver().a_coordinate(fvSolver().m_activeCellIds[ac], 1) > -0.8) continue;
              if(fvSolver().a_level(fvSolver().m_activeCellIds[ac]) != fvSolver().maxRefinementLevel()) continue;
              fvSolver().m_sweepStartFirstCell = fvSolver().m_activeCellIds[ac];
              break;
            }
          }
          m_log << "Structured Output starting from Cell ... :" << endl;
          m_log << "cellId: " << fvSolver().m_sweepStartFirstCell << endl;
          m_log << "x: " << fvSolver().a_coordinate(fvSolver().m_sweepStartFirstCell, 0) << endl;
          m_log << "y: " << fvSolver().a_coordinate(fvSolver().m_sweepStartFirstCell, 1) << endl;

          cerr << "Structured Output starting from Cell ... :" << endl;
          cerr << "cellId: " << fvSolver().m_sweepStartFirstCell << endl;
          cerr << "x: " << fvSolver().a_coordinate(fvSolver().m_sweepStartFirstCell, 0) << endl;
          cerr << "y: " << fvSolver().a_coordinate(fvSolver().m_sweepStartFirstCell, 1) << endl;

          FILE* pv6;
          stringstream StartEndCoord;
          StartEndCoord << "out/StartEndCoord_" << currentCycle << "_" << (MInt)sample;
          pv6 = fopen((StartEndCoord.str()).c_str(), "w");
          current = fvSolver().m_sweepStartFirstCell;
          fprintf(pv6, "%f ", fvSolver().a_coordinate(current, 0));
          fprintf(pv6, "%f ", fvSolver().a_coordinate(current, 1));
          // writing next cell values in x -direction = row direction
          while(fvSolver().a_hasNeighbor(current, 1) > 0) { // number of Neighbors in +x direction > 0
            current = fvSolver().c_neighborId(current, 1);  // get next neighbor as current cell Id
            fprintf(pv6, "%f ", fvSolver().a_coordinate(current, 0));
            fprintf(pv6, "%f ", fvSolver().a_coordinate(current, 1));
          }
          fprintf(pv6, "\n");

          // get next cell in +y direction from first cell = fvSolver().m_sweepStartFirstCell
          sweepStartRow = fvSolver().m_sweepStartFirstCell;
          while(fvSolver().a_hasNeighbor(sweepStartRow, 3) > 0) {
            sweepStartRow = fvSolver().c_neighborId(sweepStartRow, 3); // get next cell in +y direction from first cell
                                                                       // = fvSolver().m_sweepStartFirstCell
            // check if the geometry continues in -x direction
            while(fvSolver().a_hasNeighbor(sweepStartRow, 0) > 0
                  && fvSolver().a_hasNeighbor(fvSolver().c_neighborId(sweepStartRow, 0), 3) > 0) {
              sweepStartRow = fvSolver().c_neighborId(sweepStartRow, 0); // get next cell in -x direction
            }
            current = sweepStartRow;
            fprintf(pv6, "%f ", fvSolver().a_coordinate(current, 0));
            fprintf(pv6, "%f ", fvSolver().a_coordinate(current, 1));
            // next element in +x direction
            while(fvSolver().a_hasNeighbor(current, 1) > 0) {
              current = fvSolver().c_neighborId(current, 1);
              fprintf(pv6, "%f ", fvSolver().a_coordinate(current, 0));
              fprintf(pv6, "%f ", fvSolver().a_coordinate(current, 1));
            }
            // check if the geometry continues in +y direction by searching cells in +x direction
            while(fvSolver().a_hasNeighbor(sweepStartRow, 3) == 0 && fvSolver().a_hasNeighbor(sweepStartRow, 1) > 0) {
              sweepStartRow = fvSolver().c_neighborId(sweepStartRow, 1); // get next cell in +x direction
            }
            fprintf(pv6, "\n");
          }
          fclose(pv6);
        }
        // skip output if samplingStartCycle not reached
        if(((MInt)sample < (MInt)(samplingStartCycle * fvSolver().m_samplesPerCycle))
           || globalTimeStep <= (MInt)(fvSolver().m_restartTimeStep + fvSolver().m_noTimeStepsBetweenSamples))
          break;
        if(lsSolver().m_noSets > 1) {
          mTerm(1, AT_,
                "This routine should only be called if lsSolver().m_noSets = 1, which is not the case here... "
                "Please "
                "check!");
        }
        // compute local flame front amplitude
        // computeFlameFrontAmplitude(currentCycle,sample,0);

        if(((MInt)sample < (MInt)(samplingStartCycle * fvSolver().m_samplesPerCycle))) break;

        // write out primitive variables and G
        FILE* pv0;
        FILE* pv1;
        FILE* pv2;
        FILE* pv3;
        FILE* pv4;
        FILE* pv5;
        FILE* pv6;
        stringstream u, v, p, dpdt, rho, c, G;
        u << "out/u_" << currentCycle << "_" << (MInt)sample;
        v << "out/v_" << currentCycle << "_" << (MInt)sample;
        p << "out/p_" << currentCycle << "_" << (MInt)sample;
        dpdt << "out/dpdt_" << currentCycle << "_" << (MInt)sample;
        rho << "out/rho_" << currentCycle << "_" << (MInt)sample;
        c << "out/c_" << currentCycle << "_" << (MInt)sample;
        G << "out/G_" << currentCycle << "_" << (MInt)sample;

        pv0 = fopen((u.str()).c_str(), "w");
        pv1 = fopen((v.str()).c_str(), "w");
        pv2 = fopen((rho.str()).c_str(), "w");
        pv3 = fopen((p.str()).c_str(), "w");
        pv4 = fopen((c.str()).c_str(), "w");
        pv5 = fopen((G.str()).c_str(), "w");
        pv6 = fopen((dpdt.str()).c_str(), "w");
        if(fvSolver().m_sweepStartFirstCell < 0) {
          MString errorMessage = "sweep start cell is negative";
          mTerm(1, AT_, errorMessage);
        }
        // write the first row from first founded cell in the flame tube
        current = fvSolver().m_sweepStartFirstCell;

        // writing first element on row
        fprintf(pv0, "%f ", fvSolver().a_pvariable(current, 0));
        fprintf(pv1, "%f ", fvSolver().a_pvariable(current, 1));
        fprintf(pv2, "%f ", fvSolver().a_pvariable(current, 2));
        fprintf(pv3, "%f ", fvSolver().a_pvariable(current, 3));
        fprintf(pv4, "%f ", fvSolver().a_pvariable(current, 4));
        fprintf(pv5, "%f ", lsSolver().a_levelSetFunctionG(fv2lsId(current), 0));
        // compute the velocities
        fRho = F1 / fvSolver().a_oldVariable(current, fvSolver().CV->RHO);
        velPOW2 = F0;
        for(MInt i = 0; i < nDim; i++) {
          vel = fvSolver().a_oldVariable(current, fvSolver().CV->RHO_VV[i]) * fRho;
          velPOW2 += POW2(vel);
        }

        dPdT = fvSolver().a_pvariable(current, 3);
        // compute primitive old pressure from conservative variables
        pold = gammaMinusOne
               * (fvSolver().a_oldVariable(current, fvSolver().CV->RHO_E)
                  - F1B2 * fvSolver().a_oldVariable(current, fvSolver().CV->RHO) * velPOW2);
        dPdT -= pold;
        dPdT *= FtimeStep;
        fprintf(pv6, "%f ", dPdT);
        // writing next cell values in x -direction = row direction
        while(fvSolver().a_hasNeighbor(current, 1) > 0) { // number of Neighbors in +x direction > 0
          current = fvSolver().c_neighborId(current, 1);  // get next neighbor as current cell Id
          fprintf(pv0, "%f ", fvSolver().a_pvariable(current, 0));
          fprintf(pv1, "%f ", fvSolver().a_pvariable(current, 1));
          fprintf(pv2, "%f ", fvSolver().a_pvariable(current, 2));
          fprintf(pv3, "%f ", fvSolver().a_pvariable(current, 3));
          fprintf(pv4, "%f ", fvSolver().a_pvariable(current, 4));
          fprintf(pv5, "%f ", lsSolver().a_levelSetFunctionG(fv2lsId(current), 0));
          // compute the velocities
          fRho = F1 / fvSolver().a_oldVariable(current, fvSolver().CV->RHO);
          velPOW2 = F0;
          for(MInt i = 0; i < nDim; i++) {
            vel = fvSolver().a_oldVariable(current, fvSolver().CV->RHO_VV[i]) * fRho;
            velPOW2 += POW2(vel);
          }

          dPdT = fvSolver().a_pvariable(current, 3);
          // compute primitive old pressure from conservative variables
          pold = gammaMinusOne
                 * (fvSolver().a_oldVariable(current, fvSolver().CV->RHO_E)
                    - F1B2 * fvSolver().a_oldVariable(current, fvSolver().CV->RHO) * velPOW2);
          dPdT -= pold;
          dPdT *= FtimeStep;
          fprintf(pv6, "%f ", dPdT);
        }
        // space between this row and next row -> column direction is y direction
        fprintf(pv0, "\n");
        fprintf(pv1, "\n");
        fprintf(pv2, "\n");
        fprintf(pv3, "\n");
        fprintf(pv4, "\n");
        fprintf(pv5, "\n");
        fprintf(pv6, "\n");
        // get next cell in +y direction from first cell = fvSolver().m_sweepStartFirstCell
        sweepStartRow = fvSolver().m_sweepStartFirstCell;
        while(fvSolver().a_hasNeighbor(sweepStartRow, 3) > 0) {
          sweepStartRow = fvSolver().c_neighborId(sweepStartRow, 3);
          // check if the geometry continues in -x direction
          while(fvSolver().a_hasNeighbor(sweepStartRow, 0) > 0
                && fvSolver().a_hasNeighbor(fvSolver().c_neighborId(sweepStartRow, 0), 3) > 0) {
            sweepStartRow = fvSolver().c_neighborId(sweepStartRow, 0); // get next cell in -x direction
          }
          current = sweepStartRow;
          // write the next row
          fprintf(pv0, "%f ", fvSolver().a_pvariable(current, 0));
          fprintf(pv1, "%f ", fvSolver().a_pvariable(current, 1));
          fprintf(pv2, "%f ", fvSolver().a_pvariable(current, 2));
          fprintf(pv3, "%f ", fvSolver().a_pvariable(current, 3));
          fprintf(pv4, "%f ", fvSolver().a_pvariable(current, 4));
          fprintf(pv5, "%f ", lsSolver().a_levelSetFunctionG(fv2lsId(current), 0));
          // compute the velocities
          fRho = F1 / fvSolver().a_oldVariable(current, fvSolver().CV->RHO);
          velPOW2 = F0;
          for(MInt i = 0; i < nDim; i++) {
            vel = fvSolver().a_oldVariable(current, fvSolver().CV->RHO_VV[i]) * fRho;
            velPOW2 += POW2(vel);
          }

          dPdT = fvSolver().a_pvariable(current, 3);
          // compute primitive old pressure from conservative variables
          pold = gammaMinusOne
                 * (fvSolver().a_oldVariable(current, fvSolver().CV->RHO_E)
                    - F1B2 * fvSolver().a_oldVariable(current, fvSolver().CV->RHO) * velPOW2);
          dPdT -= pold;
          dPdT *= FtimeStep;
          fprintf(pv6, "%f ", dPdT);

          // next element in +x direction
          while(fvSolver().a_hasNeighbor(current, 1) > 0) {
            current = fvSolver().c_neighborId(current, 1);
            fprintf(pv0, "%f ", fvSolver().a_pvariable(current, 0));
            fprintf(pv1, "%f ", fvSolver().a_pvariable(current, 1));
            fprintf(pv2, "%f ", fvSolver().a_pvariable(current, 2));
            fprintf(pv3, "%f ", fvSolver().a_pvariable(current, 3));
            fprintf(pv4, "%f ", fvSolver().a_pvariable(current, 4));
            fprintf(pv5, "%f ", lsSolver().a_levelSetFunctionG(fv2lsId(current), 0));

            // compute the velocities
            fRho = F1 / fvSolver().a_oldVariable(current, fvSolver().CV->RHO);
            velPOW2 = F0;
            for(MInt i = 0; i < nDim; i++) {
              vel = fvSolver().a_oldVariable(current, fvSolver().CV->RHO_VV[i]) * fRho;
              velPOW2 += POW2(vel);
            }
            dPdT = fvSolver().a_pvariable(current, 3);
            // compute primitive old pressure from conservative variables
            pold = gammaMinusOne
                   * (fvSolver().a_oldVariable(current, fvSolver().CV->RHO_E)
                      - F1B2 * fvSolver().a_oldVariable(current, fvSolver().CV->RHO) * velPOW2);
            dPdT -= pold;
            dPdT *= FtimeStep;
            fprintf(pv6, "%f ", dPdT);
          }

          fprintf(pv0, "\n");
          fprintf(pv1, "\n");
          fprintf(pv2, "\n");
          fprintf(pv3, "\n");
          fprintf(pv4, "\n");
          fprintf(pv5, "\n");
          fprintf(pv6, "\n");

          // check if the geometry continues in +y direction by searching cells in +x direction
          while(fvSolver().a_hasNeighbor(sweepStartRow, 3) == 0 && fvSolver().a_hasNeighbor(sweepStartRow, 1) > 0) {
            sweepStartRow = fvSolver().c_neighborId(sweepStartRow, 1); // get next cell in +x direction
          }
        }
        fclose(pv0);
        fclose(pv1);
        fclose(pv2);
        fclose(pv3);
        fclose(pv4);
        fclose(pv5);
        fclose(pv6);

        break;
      }
        // netcdf output
      case 5: {
        {
          const MFloat gammaMinusOne = fvSolver().m_gamma - F1;
          const MFloat FgammaMinusOne = F1 / gammaMinusOne;
          MFloat reactionEnthalpy = (fvSolver().m_burntUnburntTemperatureRatio - F1) * FgammaMinusOne;
          MFloat halfWidth = 9999999.0;
          const MFloat eps = fvSolver().c_cellLengthAtLevel(fvSolver().maxRefinementLevel()) * 0.00000000001;
          MInt cnt = 0;
          MInt nghbrId = -1;
          MFloat check = 9999.0;
          MFloatScratchSpace heatRelease(fvSolver().noInternalCells(), AT_, "heatRelease");
          MFloatScratchSpace pressure(fvSolver().noInternalCells(), AT_, "pressure");
          MFloatScratchSpace coords(fvSolver().noInternalCells(), AT_, "coords");

          // here comes our line output
          for(MInt cell = 0; cell < fvSolver().noInternalCells(); cell++) {
            if(!fvSolver().a_hasProperty(cell, SolverCell::IsOnCurrentMGLevel)) continue;

            halfWidth = F1B2 * fvSolver().c_cellLengthAtCell(cell) + eps;

            if(fvSolver().a_coordinate(cell, 0) > halfWidth || fvSolver().a_coordinate(cell, 0) < F0) continue;

            // -x neighbor
            nghbrId = fvSolver().c_neighborId(cell, 0);

            check = fvSolver().a_coordinate(cell, 0) + fvSolver().a_coordinate(nghbrId, 0);

            if(fabs(check) > eps) mTerm(1, AT_, "ERROR: coord check failed");

            // interpolate to center line
            pressure[cnt] = fvSolver().a_pvariable(cell, 3);
            pressure[cnt] += fvSolver().a_pvariable(nghbrId, 3);
            pressure[cnt] *= F1B2;

            heatRelease[cnt] = fvSolver().a_reactionRate(cell, 0) * fvSolver().a_cellVolume(cell) * reactionEnthalpy;
            heatRelease[cnt] +=
                fvSolver().a_reactionRate(nghbrId, 0) * fvSolver().a_cellVolume(nghbrId) * reactionEnthalpy;
            heatRelease[cnt] *= F1B2;

            coords[cnt] = fvSolver().a_coordinate(cell, 1);

            cnt++;
          }
          // communicate center line data if available

          MInt root = 0;
          MFloat totalCnt = 0;
          MIntScratchSpace globalCnt(fvSolver().noDomains(), AT_, "globalCnt");
          MIntScratchSpace offsetIO(fvSolver().noDomains(), AT_, "offsetIO");

          MPI_Gather(&cnt, 1, MPI_INT, globalCnt.getPointer(), 1, MPI_INT, root, fvSolver().mpiComm(), AT_, "cnt",
                     "globalCnt.getPointer()");

          for(MInt i = 0; i < fvSolver().noDomains(); i++) {
            offsetIO[i] = totalCnt;
            totalCnt += globalCnt[i];
          }
          if(fvSolver().domainId() != root) totalCnt = 0;

          MFloatScratchSpace tmp(totalCnt, AT_, "tmp");
          MFloatScratchSpace tmp1(totalCnt, AT_, "tmp1");
          MFloatScratchSpace tmp2(totalCnt, AT_, "tmp2");

          MPI_Gatherv(pressure.getPointer(), cnt, MPI_DOUBLE, tmp.getPointer(), globalCnt.getPointer(),
                      offsetIO.getPointer(), MPI_DOUBLE, root, fvSolver().mpiComm(), AT_, "pressure.getPointer()",
                      "tmp.getPointer()");
          MPI_Gatherv(heatRelease.getPointer(), cnt, MPI_DOUBLE, tmp1.getPointer(), globalCnt.getPointer(),
                      offsetIO.getPointer(), MPI_DOUBLE, root, fvSolver().mpiComm(), AT_, "heatRelease.getPointer()",
                      "tmp1.getPointer()");
          MPI_Gatherv(coords.getPointer(), cnt, MPI_DOUBLE, tmp2.getPointer(), globalCnt.getPointer(),
                      offsetIO.getPointer(), MPI_DOUBLE, root, fvSolver().mpiComm(), AT_, "coords.getPointer()",
                      "tmp2.getPointer()");

          FILE* pv0;
          FILE* pv1;
          FILE* pv2;

          stringstream p;
          stringstream h;
          stringstream y;

          p << "out/centerline_p";
          h << "out/centerline_h";
          y << "out/centerline_y";
          if(fvSolver().domainId() == root) {
            pv0 = fopen((p.str()).c_str(), "a+");
            pv1 = fopen((h.str()).c_str(), "a+");
            pv2 = fopen((y.str()).c_str(), "w");

            for(MInt c = 0; c < totalCnt; c++) {
              fprintf(pv0, "%f ", tmp[c]);
              fprintf(pv0, "\n");
              fprintf(pv1, "%f ", tmp1[c]);
              fprintf(pv1, "\n");
              fprintf(pv2, "%f ", tmp2[c]);
              fprintf(pv2, "\n");
            }

            fclose(pv0);
            fclose(pv1);
            fclose(pv2);
          }
        }
        // output when time is equal to a sample time
        if(globalTimeStep % fvSolver().m_noTimeStepsBetweenSamples != 0) {
          return;
        }
        // skip output if samplingStartCycle not reached
        if(((MInt)sample < (MInt)(samplingStartCycle * fvSolver().m_samplesPerCycle))) break;

        stringstream varFileName;
        varFileName << fvSolver().outputDir() << "Q_" << currentCycle << "_" << (MInt)sample << ParallelIo::fileExt();

        MFloatScratchSpace dbVariables(fvSolver().a_noCells() * (fvSolver().CV->noVariables + 5), AT_, "dbVariables");
        MIntScratchSpace idVariables(0, AT_, "idVariables");
        MFloatScratchSpace dbParameters(5, AT_, "dbParameters");
        MIntScratchSpace idParameters(4, AT_, "idParameters");
        vector<MString> dbVariablesName;
        vector<MString> idVariablesName;
        vector<MString> dbParametersName;
        vector<MString> idParametersName;
        vector<MString> name;

        if(fvSolver().m_levelSet) {
          MFloatScratchSpace levelSetFunction(fvSolver().a_noCells(), AT_, "levelSetFunction");
          MFloatScratchSpace curvature(fvSolver().a_noCells(), AT_, "curvature");

          for(MInt cell = 0; cell < fvSolver().a_noCells(); cell++) {
            const MInt gCellId = fv2lsId(cell);
            levelSetFunction[cell] = lsSolver().a_levelSetFunctionG(gCellId, 0);
            curvature[cell] = lsSolver().a_curvatureG(gCellId, 0);
          }
          stringstream gName;
          stringstream gCurv;
          gName << "G";
          gCurv << "curv";

          name.push_back(gName.str());
          fvSolver().collectVariables(levelSetFunction.begin(), dbVariables, name, dbVariablesName, 1,
                                      fvSolver().a_noCells());
          name.clear();
          name.push_back(gCurv.str());
          fvSolver().collectVariables(curvature.begin(), dbVariables, name, dbVariablesName, 1, fvSolver().a_noCells());
        }
        {
          const MFloat gammaMinusOne = fvSolver().m_gamma - F1;
          const MFloat FgammaMinusOne = F1 / gammaMinusOne;
          MFloat reactionEnthalpy = (fvSolver().m_burntUnburntTemperatureRatio - F1) * FgammaMinusOne;
          MFloat FtimeStep = 1 / fvSolver().timeStep();
          MFloat fRho;
          MFloat vel = -9999999.0;
          MFloat velPOW2 = -9999999.0;
          MFloat pold;

          MFloatScratchSpace dPdT(fvSolver().a_noCells(), AT_, "dPdT");
          MFloatScratchSpace dHdT(fvSolver().a_noCells(), AT_, "dHdT");
          MFloatScratchSpace h(fvSolver().a_noCells(), AT_, "h");

          for(MInt cell = 0; cell < fvSolver().a_noCells(); cell++) {
            dPdT[cell] = F0;
            dHdT[cell] = F0;
            h[cell] = F0;

            if(!fvSolver().a_hasProperty(cell, SolverCell::IsOnCurrentMGLevel)) continue;

            // compute the velocities
            fRho = F1 / fvSolver().a_oldVariable(cell, fvSolver().CV->RHO);
            velPOW2 = F0;
            for(MInt i = 0; i < nDim; i++) {
              vel = fvSolver().a_oldVariable(cell, fvSolver().CV->RHO_VV[i]) * fRho;
              velPOW2 += POW2(vel);
            }

            dPdT[cell] = fvSolver().a_pvariable(cell, 3);
            // compute primitive old pressure from conservative variables
            pold = gammaMinusOne
                   * (fvSolver().a_oldVariable(cell, fvSolver().CV->RHO_E)
                      - F1B2 * fvSolver().a_oldVariable(cell, fvSolver().CV->RHO) * velPOW2);
            dPdT[cell] -= pold;
            dPdT[cell] *= FtimeStep;

            h[cell] = fvSolver().a_reactionRate(cell, 0);
            h[cell] *= fvSolver().a_cellVolume(cell) * reactionEnthalpy;
            dHdT[cell] = fvSolver().a_reactionRate(cell, 0);
            dHdT[cell] -= fvSolver().a_reactionRateBackup(cell, 0);
            dHdT[cell] *= fvSolver().a_cellVolume(cell) * reactionEnthalpy;
            dHdT[cell] *= FtimeStep;
          }
          stringstream dPdTname;
          stringstream dHdTname;
          stringstream hName;
          dPdTname << "dPdT";
          dHdTname << "dHdT";
          hName << "h";
          name.clear();
          name.push_back(dPdTname.str());
          fvSolver().collectVariables(dPdT.begin(), dbVariables, name, dbVariablesName, 1, fvSolver().a_noCells());
          name.clear();
          name.push_back(dHdTname.str());
          fvSolver().collectVariables(dHdT.begin(), dbVariables, name, dbVariablesName, 1, fvSolver().a_noCells());
          name.clear();
          name.push_back(hName.str());
          fvSolver().collectVariables(h.begin(), dbVariables, name, dbVariablesName, 1, fvSolver().a_noCells());
        }

        name.clear();
        for(MInt v = 0; v < fvSolver().PV->noVariables; v++) {
          name.push_back(fvSolver().m_variablesName[v]);
        }
        fvSolver().collectVariables(&fvSolver().a_pvariable(0, 0), dbVariables, name, dbVariablesName,
                                    fvSolver().PV->noVariables, fvSolver().a_noCells());

        fvSolver().setRestartFileOutputTimeStep();
        fvSolver().collectParameters(fvSolver().m_noSamples, idParameters, "noSamples", idParametersName);
        fvSolver().collectParameters(globalTimeStep, idParameters, "globalTimeStep", idParametersName);
        fvSolver().collectParameters(fvSolver().m_time, dbParameters, "time", dbParametersName);
        fvSolver().collectParameters(fvSolver().m_restartFileOutputTimeStep, dbParameters, "timeStep",
                                     dbParametersName);
        fvSolver().collectParameters(fvSolver().m_noTimeStepsBetweenSamples, idParameters, "noTimeStepsBetweenSamples",
                                     idParametersName);
        fvSolver().collectParameters(fvSolver().m_physicalTime, dbParameters, "physicalTime", dbParametersName);
        fvSolver().collectParameters((MInt)fvSolver().m_forcing, idParameters, "forcing", idParametersName);


        m_log << "Writing structured output at time setp " << globalTimeStep << endl;
        switch(string2enum(fvSolver().m_outputFormat)) {
          case NETCDF: {
            fvSolver().saveGridFlowVarsPar((varFileName.str()).c_str(), fvSolver().a_noCells(),
                                           fvSolver().noInternalCells(), dbVariables, dbVariablesName, 0, idVariables,
                                           idVariablesName, 0, dbParameters, dbParametersName, idParameters,
                                           idParametersName, fvSolver().m_recalcIds);
            break;
          }
          default: {
            mTerm(1, AT_, "change solution output format to NETCDF or change code");
            break;
          }
        }
        break;
      }
      default: {
        // find cell on r_max in the left lower corner
        // order of output: rows: x; columns: y
        for(MInt ac = 0; ac < fvSolver().m_noActiveCells; ac++) {
          if(fvSolver().a_hasNeighbor(fvSolver().m_activeCellIds[ac], 0) > 0) continue;
          if(fvSolver().a_hasNeighbor(fvSolver().m_activeCellIds[ac], 2) > 0) continue;
          if(ABS(fvSolver().a_coordinate(fvSolver().m_activeCellIds[ac], 0)) > 0.5) continue;
          if(fvSolver().a_level(fvSolver().m_activeCellIds[ac]) != fvSolver().maxRefinementLevel()) continue;
          sweepStart = fvSolver().m_activeCellIds[ac];
          break;
        }

        // write out primitive variables and G
        FILE* pv0;
        FILE* pv1;
        FILE* pv2;
        FILE* pv3;
        FILE* pv4;
        FILE* pv5;
        stringstream u, v, p, rho, c, G;
        u << "out/u_" << currentCycle << "_" << (MInt)sample;
        v << "out/v_" << currentCycle << "_" << (MInt)sample;
        p << "out/p_" << currentCycle << "_" << (MInt)sample;
        rho << "out/rho_" << currentCycle << "_" << (MInt)sample;
        c << "out/c_" << currentCycle << "_" << (MInt)sample;
        G << "out/G_" << currentCycle << "_" << (MInt)sample;
        pv0 = fopen((u.str()).c_str(), "w");
        pv1 = fopen((v.str()).c_str(), "w");
        pv2 = fopen((rho.str()).c_str(), "w");
        pv3 = fopen((p.str()).c_str(), "w");
        pv4 = fopen((c.str()).c_str(), "w");
        pv5 = fopen((G.str()).c_str(), "w");
        // write the first row
        current = sweepStart;
        fprintf(pv0, "%f ", fvSolver().a_pvariable(current, 0));
        fprintf(pv1, "%f ", fvSolver().a_pvariable(current, 1));
        fprintf(pv2, "%f ", fvSolver().a_pvariable(current, 2));
        fprintf(pv3, "%f ", fvSolver().a_pvariable(current, 3));
        fprintf(pv4, "%f ", fvSolver().a_pvariable(current, 4));
        fprintf(pv5, "%f ", lsSolver().a_levelSetFunctionG(fv2lsId(current), 0));
        while(fvSolver().a_hasNeighbor(current, 1) > 0) {
          current = fvSolver().c_neighborId(current, 1);
          if(ABS(fvSolver().a_coordinate(current, 0)) > 0.5) break;
          fprintf(pv0, "%f ", fvSolver().a_pvariable(current, 0));
          fprintf(pv1, "%f ", fvSolver().a_pvariable(current, 1));
          fprintf(pv2, "%f ", fvSolver().a_pvariable(current, 2));
          fprintf(pv3, "%f ", fvSolver().a_pvariable(current, 3));
          fprintf(pv4, "%f ", fvSolver().a_pvariable(current, 4));
          fprintf(pv5, "%f ", lsSolver().a_levelSetFunctionG(fv2lsId(current), 0));
        }
        fprintf(pv0, "\n");
        fprintf(pv1, "\n");
        fprintf(pv2, "\n");
        fprintf(pv3, "\n");
        fprintf(pv4, "\n");
        fprintf(pv5, "\n");

        while(fvSolver().a_hasNeighbor(sweepStart, 3) > 0) {
          sweepStart = fvSolver().c_neighborId(sweepStart, 3);
          current = sweepStart;
          // write the next row
          fprintf(pv0, "%f ", fvSolver().a_pvariable(current, 0));
          fprintf(pv1, "%f ", fvSolver().a_pvariable(current, 1));
          fprintf(pv2, "%f ", fvSolver().a_pvariable(current, 2));
          fprintf(pv3, "%f ", fvSolver().a_pvariable(current, 3));
          fprintf(pv4, "%f ", fvSolver().a_pvariable(current, 4));
          fprintf(pv5, "%f ", lsSolver().a_levelSetFunctionG(fv2lsId(current), 0));
          while(fvSolver().a_hasNeighbor(current, 1) > 0) {
            current = fvSolver().c_neighborId(current, 1);
            if(ABS(fvSolver().a_coordinate(current, 0)) > 0.5) break;
            fprintf(pv0, "%f ", fvSolver().a_pvariable(current, 0));
            fprintf(pv1, "%f ", fvSolver().a_pvariable(current, 1));
            fprintf(pv2, "%f ", fvSolver().a_pvariable(current, 2));
            fprintf(pv3, "%f ", fvSolver().a_pvariable(current, 3));
            fprintf(pv4, "%f ", fvSolver().a_pvariable(current, 4));
            fprintf(pv5, "%f ", lsSolver().a_levelSetFunctionG(fv2lsId(current), 0));
          }
          fprintf(pv0, "\n");
          fprintf(pv1, "\n");
          fprintf(pv2, "\n");
          fprintf(pv3, "\n");
          fprintf(pv4, "\n");
          fprintf(pv5, "\n");
        }
        fclose(pv0);
        fclose(pv1);
        fclose(pv2);
        fclose(pv3);
        fclose(pv4);
        fclose(pv5);
        break;
      }
    }
  }
}


// works only with zeroth level-set function!
template <MInt nDim, class SysEqn>
void LsFvCombustion<nDim, SysEqn>::computeSourceTerms() {
  TRACE();


  // get all the parameters from the fv solver...
  MInt m_constantFlameSpeed = fvSolver().m_constantFlameSpeed;
  MFloat m_subfilterVariance = fvSolver().m_subfilterVariance;
  MFloat m_rhoEInfinity = fvSolver().m_rhoEInfinity;
  MInt m_useCorrectedBurningVelocity = fvSolver().m_useCorrectedBurningVelocity;
  MFloat m_integralLengthScale = fvSolver().m_integralLengthScale;
  MFloat m_marksteinLength = fvSolver().m_marksteinLength;
  MFloat m_dampingDistanceFlameBase = fvSolver().m_dampingDistanceFlameBase;
  MFloat m_sutherlandConstant = fvSolver().m_sutherlandConstant;
  MFloat m_sutherlandPlusOne = fvSolver().m_sutherlandPlusOne;
  MFloat m_TInfinity = fvSolver().m_TInfinity;
  MFloat m_flameSpeed = fvSolver().m_flameSpeed;
  MFloat m_turbFlameSpeed = fvSolver().m_turbFlameSpeed;
  MFloat m_noReactionCells = fvSolver().m_noReactionCells;
  MFloat m_NuT = fvSolver().m_NuT;
  MFloat m_ScT = fvSolver().m_ScT;
  MFloat m_Pr = fvSolver().m_Pr;
  MFloat m_integralAmplitude = fvSolver().m_integralAmplitude;
  MFloat m_Re0 = fvSolver().m_sysEqn.m_Re0;
  MFloat m_yOffsetFlameTube = fvSolver().m_yOffsetFlameTube;
  MInt m_temperatureChange = fvSolver().m_temperatureChange;
  MInt m_heatReleaseDamp = fvSolver().m_heatReleaseDamp;
  MInt m_bndryGhostCellsOffset = fvSolver().m_bndryGhostCellsOffset;
  MInt m_initialCondition = fvSolver().m_initialCondition;
  MFloat m_rhoUnburnt = fvSolver().m_rhoUnburnt;
  MFloat m_rhoFlameTube = fvSolver().m_rhoFlameTube;
  MInt m_restartTimeStep = fvSolver().m_restartTimeStep;
  //  #define debugOutput

  const MInt noCells = fvSolver().a_noCells();
  MInt gc;
  MFloat factor, Psi, xi = F0, cbar, a = F0, a1 = F0, a2 = F0, b = F0, b2 = F0, b1 = F0, c1 = F0, FD, Rr,
                      sigma; //,d,omega
  MFloat reactionEnthalpy, rhoBar;
  const MFloat gammaMinusOne = fvSolver().m_gamma - F1;
  const MFloat FgammaMinusOne = F1 / gammaMinusOne;
  MFloat denominator;
  MFloat FlaminarFlameThickness = F1 / fvSolver().m_laminarFlameThickness;
  MFloat rhoU, Dth;
  MFloat rhoUFrhoB = fvSolver().m_burntUnburntTemperatureRatio;
  MFloat rhoBurnt = fvSolver().m_rhoInfinity / fvSolver().m_burntUnburntTemperatureRatio;
  MFloat rhoJump = F1 / fvSolver().m_burntUnburntTemperatureRatio - F1;
  MFloat FrhoBurnt = fvSolver().m_burntUnburntTemperatureRatio;
  MFloat factor1 = fvSolver().m_c0 / (F1 - fvSolver().m_c0);
  MFloat echekkiFerzigerPrefactor = F1 / (F1 - fvSolver().m_c0) * (F1 / (F1 - fvSolver().m_c0) - F1);
  const MFloat sq1F2 = sqrt(F1B2);
  const MFloat eps = fvSolver().c_cellLengthAtLevel(fvSolver().maxRefinementLevel()) * 0.00000000001;
  MFloat levelSetNegative = F0, levelSetPlus = F0;

  if(lsSolver().m_noSets > 1) {
    mTerm(1, AT_,
          "This routine should only be called if lsSolver().m_noSets = 1, which is not the case here... Please "
          "check!");
  }

  //---
  // reset
  fvSolver().m_maxReactionRate = F0;
  fvSolver().m_totalHeatReleaseRate = F0;
  // compute the heat release
  reactionEnthalpy = (fvSolver().m_burntUnburntTemperatureRatio - F1) * FgammaMinusOne;
  // save old reaction rate
  if(fvSolver().m_RKStep == 0) {
    // if(hasReactionRates() && hasReactionRatesBackup())
    memcpy(&fvSolver().a_reactionRateBackup(0, 0), &fvSolver().a_reactionRate(0, 0),
           sizeof(MFloat) * fvSolver().noInternalCells());
  }
  // reset
  for(MInt c = 0; c < noCells; c++) {
    fvSolver().a_reactionRate(c, 0) = F0;
  }
  // compute the subfilter variance
  sigma = m_subfilterVariance * fvSolver().c_cellLengthAtLevel(fvSolver().maxLevel());
  factor = F1 / (sqrt(F2) * sigma);
  // return for no-heat release combustion
  if(reactionEnthalpy < m_rhoEInfinity * 0.00001) {
    return;
  }
  switch(m_initialCondition) {
    case 17516: {
      /** brief important correction of jump condition (implemented like presented in Dissertation D. Hartmann)
       *
       * \author Stephan Schlimpert
       * \date November 2011
       */
      MInt diverged = 0;
      MFloat diffTimeStep = 50000;
      if(m_temperatureChange && (globalTimeStep - m_restartTimeStep) <= diffTimeStep) {
        MFloat diffTemp =
            (fvSolver().m_burntUnburntTemperatureRatioEnd - fvSolver().m_burntUnburntTemperatureRatioStart);
        fvSolver().m_burntUnburntTemperatureRatio = (diffTemp / diffTimeStep) * (globalTimeStep - m_restartTimeStep)
                                                    + fvSolver().m_burntUnburntTemperatureRatioStart;
      }
      rhoBurnt = m_rhoFlameTube / fvSolver().m_burntUnburntTemperatureRatio;
      FrhoBurnt = fvSolver().m_burntUnburntTemperatureRatio * F1 / m_rhoFlameTube;
      rhoJump = F1 - fvSolver().m_burntUnburntTemperatureRatio;
      for(MInt c = 0; c < m_bndryGhostCellsOffset; c++) {
        if(fvSolver().grid().tree().hasProperty(c, Cell::IsHalo)) {
          continue;
        }
        if(!fvSolver().a_hasProperty(c, SolverCell::IsOnCurrentMGLevel)) {
          continue;
        }
        MInt cLs = fv2lsId(c);
        if(cLs < 0) {
          continue;
        }
        gc = cLs;
        // MInt gcFv= ls2fvId(gc);
        if(gc == -1) {
          continue;
        }
        // continue if the g cell is out of the band
        // the source term is in this case zero
        if(ABS(lsSolver().m_outsideGValue - ABS(lsSolver().a_levelSetFunctionG(gc, 0))) < 0.02) {
          continue;
        }
        if(lsSolver().a_level(gc) != fvSolver().maxRefinementLevel()) {
          continue;
        }
        // compute the Echekki-Ferziger constant
        // - compute the inverse eddy viscosity (s/m^2)
        // - - density of the unburnt gas rho^u = fvSolver().m_rhoInfinity
        // - - assuming m_TInfinity is the temperature of the unburnt gas -> muInf=SUTHERLANDLAW(m_TInfinity)
        // - - DthInf = muInf^u / ( rho^u *Pr )
        FD = F1 / fvSolver().m_DthInfinity;
        // - compute Rr
        levelSetNegative = lsSolver().a_levelSetFunctionG(gc, 0) - m_noReactionCells;
        levelSetPlus = lsSolver().a_levelSetFunctionG(gc, 0) + m_noReactionCells;

        if(m_useCorrectedBurningVelocity) {
          if(lsSolver().m_sharpDamp) {
            Rr = echekkiFerzigerPrefactor * POW2(lsSolver().a_correctedBurningVelocity(gc, 0)) * FD * F1B4
                 * (1 + tanh(levelSetPlus * 100.0)) * (1 - tanh(levelSetNegative * 100.0));
          } else {
            Rr = echekkiFerzigerPrefactor * POW2(lsSolver().a_correctedBurningVelocity(gc, 0)) * FD;
          }
        } else {
          if(lsSolver().m_sharpDamp) {
            Rr = echekkiFerzigerPrefactor
                 * POW2(lsSolver().a_flameSpeedG(gc, 0) * (F1 - lsSolver().a_curvatureG(gc, 0) * m_marksteinLength))
                 * FD * F1B4 * (1 + tanh(levelSetPlus * 100.0)) * (1 - tanh(levelSetNegative * 100.0));
          } else {
            Rr = echekkiFerzigerPrefactor
                 * POW2(lsSolver().a_flameSpeedG(gc, 0) * (F1 - lsSolver().a_curvatureG(gc, 0) * m_marksteinLength))
                 * FD;
          }
        }

        // damp out the reaction rate to the wall w(y=0.0341) = 0 by damping the model constant Rr(y=0.0341)=0
        if(m_heatReleaseDamp) {
          if(lsSolver().c_coordinate(gc, 1) < (0.0234 + m_dampingDistanceFlameBase)) {
            if(lsSolver().c_coordinate(gc, 1) > 0.0234) {
              Rr *= 1. / m_dampingDistanceFlameBase * (lsSolver().c_coordinate(gc, 1) - 0.0234);
            }
          }
          if(lsSolver().c_coordinate(gc, 1) < 0.0234) {
            Rr = F0;
          }
        }
        // compute Psi
        // - compute xi
        xi = lsSolver().a_levelSetFunctionG(gc, 0) * factor; //< xi = G(x,t)/(sigma*sqrt(2))

        // - compute c bar
        a = xi - sq1F2 * factor1 * sigma * FlaminarFlameThickness;
        a1 = F1B2 * (POW2(sigma * factor1 * FlaminarFlameThickness))
             - SQRT2 * xi * sigma * factor1 * FlaminarFlameThickness;
        if(a > F3) {
          a2 = (F1 - fvSolver().m_c0) * exp(a1) * F1B2
               * (-erfc(a) + F2); // erfc computes 1-erf and is more accurate for large a
        } else {
          a2 = (F1 - fvSolver().m_c0) * exp(a1) * F1B2 * (erf(a) + F1);
        }
        b = xi + sq1F2 * sigma * FlaminarFlameThickness;
        b1 = F1B2 * POW2(sigma * FlaminarFlameThickness) + SQRT2 * xi * sigma * FlaminarFlameThickness;
        if(b > F3) {
          b2 = fvSolver().m_c0 * exp(b1) * F1B2 * (-erfc(b)); // erfc computes 1-erf and is more accurate for large b
        } else {
          b2 = fvSolver().m_c0 * exp(b1) * F1B2 * (erf(b) - F1);
        }
        c1 = F1B2 * (erf(xi) - F1);
        cbar = F1 - (a2 + b2 - c1);

        // - compute Psi
        if(ABS(F1 - cbar) < eps) {
          Psi = F1;
        } else {
          Psi = a2 / ((F1 - cbar) * POW2((rhoUFrhoB + cbar * rhoJump)));
        }

        fvSolver().a_psi(c) = Psi;

        // compute rhoBar
        rhoBar = m_rhoUnburnt / (1 - fvSolver().a_pvariable(c, fvSolver().PV->C) * rhoJump);

        // compute the source term
        fvSolver().a_reactionRate(c, 0) =
            m_Re0 * rhoBar * rhoUFrhoB * Rr * (F1 - fvSolver().a_pvariable(c, fvSolver().PV->C)) * Psi;

        // catch nan reaction rate
        if(!(fvSolver().a_reactionRate(c, 0) >= F0) && !(fvSolver().a_reactionRate(c, 0) <= F0)) {
          diverged = 1;
          cerr << "reaction rate is nan!!!" << endl;
          // cerr << "b=" << b << " " << "b1=" << b1 << " " << "exp(b1)=" << exp(b1) << " " << "b2=" << b2 << " " <<
          // "c1="
          //   << c1 << endl;
          // cerr << "a=" << a << " " << "a1=" << a1 << " " << "a2=" << a2 << endl;
          // cerr << "g cell info for " << gc << endl;
          // cerr << "gc_coordx=" << lsSolver().c_coordinate(gc, 0) << " " << "gc_coordy=" <<
          // lsSolver().c_coordinate(gc, 1) << " glevel "
          //   << lsSolver().a_level(gc) << " inBand " << lsSolver().a_inBandG(gc, 0) << " in shadow layer " <<
          //   lsSolver().a_inShadowLayerG(gc)
          //   << " no childs " << lsSolver().a_noChildrenG(gc) << endl;
          // cerr << "is GWindow Cell " << lsSolver().a_isGWindowCellG(gc) << endl;
          // cerr << "is GHalo Cell " << lsSolver().a_isGHaloCellG(gc) << endl;
          // cerr << "levelsetfunction=" << lsSolver().a_levelSetFunctionG(IDX_LSSET(gc, 0)] << " " << "xi=" << xi <<
          // endl; cerr << "cbar=" << cbar << " " << "Psi=" << Psi << " " << "rhoBar=" << rhoBar << " " << "c="
          ///    << fvSolver().a_pvariable(c, fvSolver().PV->C) << endl;
          // cerr << "rho=" << fvSolver().a_pvariable(c, fvSolver().PV->RHO) << endl;
          // cerr << "rhoUnburnt=" << m_rhoUnburnt << endl;
          // cerr << "global cell info for " << fvSolver().c_globalId(c) << endl;
          // cerr << "cell level " << a_level(c) << endl;
          // cerr << "window cell " << fvSolver().a_hasProperty(c, Cell::IsWindow) << endl;
          // cerr << "halo cell " << fvSolver().a_hasProperty(c, Cell::IsHalo) << endl;
          // cerr << "x=" << fvSolver().a_coordinate(c, 0) << " "
          //<< "y=" << fvSolver().a_coordinate(c, 1) << " "
          //<< "z=" << fvSolver().a_coordinate(c, 2) << endl;
          // fvSolver().a_reactionRate(c, 0) = 1000.0;
        }

        // compute the source terms
        fvSolver().a_rightHandSide(c, fvSolver().CV->RHO_C) -=
            fvSolver().a_reactionRate(c, 0) * fvSolver().a_cellVolume(c);
        fvSolver().a_rightHandSide(c, fvSolver().CV->RHO_E) -=
            fvSolver().a_reactionRate(c, 0) * fvSolver().a_cellVolume(c) * reactionEnthalpy;

        // compute the maximum reaction rate and the total heat release
        if(!fvSolver().a_hasProperty(c, Cell::IsHalo)) {
          fvSolver().m_maxReactionRate = mMax(fvSolver().a_reactionRate(c, 0), fvSolver().m_maxReactionRate);
          fvSolver().m_totalHeatReleaseRate +=
              fvSolver().a_reactionRate(c, 0) * fvSolver().a_cellVolume(c) * reactionEnthalpy;
        }
      }
      MPI_Allreduce(MPI_IN_PLACE, &diverged, 1, MPI_INT, MPI_MAX, fvSolver().mpiComm(), AT_, "MPI_IN_PLACE",
                    "diverged");
      if(diverged == 1) {
        cerr << "Solution diverged, writing NETCDF file for debugging" << endl;
        MString errorMessage = "reaction rate is nan";
        fvSolver().saveSolverSolution(1);
        stringstream fileNameVtk;
        stringstream levelSetFileName;
        levelSetFileName << "restartLevelSetNewGCellGrid_"
                         << "debug";

        stringstream levelSetFileNameSol;
        levelSetFileNameSol << "restartLevelSetNew_"
                            << "debug";

        lsSolver().writeRestartLevelSetFileCG(1, (levelSetFileName.str()).c_str(), (levelSetFileNameSol.str()).c_str());

        mTerm(1, AT_, errorMessage);
      }

      break;
    }
    case 1751600: {
      /** brief important correction of jump condition (implemented like presented in Dissertation D. Hartmann)
       *
       * \author Stephan Schlimpert
       * \date November 2011
       */
      //      MFloat
      //      minRr=10000,maxRr=-10000,minPsi=10000,maxPsi=-10000,minVol=10000,maxVol=-10000,minA2=10000,maxA2=-10000,maxCurvature=-10000,minCurvature=10000;
      MInt diverged = 0;
      MFloat diffTimeStep = 50000;
      if(m_temperatureChange && (globalTimeStep - m_restartTimeStep) <= diffTimeStep) {
        MFloat diffTemp =
            (fvSolver().m_burntUnburntTemperatureRatioEnd - fvSolver().m_burntUnburntTemperatureRatioStart);
        fvSolver().m_burntUnburntTemperatureRatio = (diffTemp / diffTimeStep) * (globalTimeStep - m_restartTimeStep)
                                                    + fvSolver().m_burntUnburntTemperatureRatioStart;
      }
      rhoBurnt = m_rhoFlameTube / fvSolver().m_burntUnburntTemperatureRatio;
      FrhoBurnt = fvSolver().m_burntUnburntTemperatureRatio * F1 / m_rhoFlameTube;
      rhoJump = F1 - fvSolver().m_burntUnburntTemperatureRatio;
      for(MInt c = 0; c < m_bndryGhostCellsOffset; c++) {
        if(fvSolver().grid().tree().hasProperty(c, Cell::IsHalo)) {
          continue;
        }
        if(!fvSolver().a_hasProperty(c, SolverCell::IsOnCurrentMGLevel)) {
          continue;
        }
        MInt cLs = fv2lsId(c);
        if(cLs < 0) {
          continue;
        }
        gc = cLs;
        if(gc == -1) {
          continue;
        }
        // continue if the g cell is out of the band
        // the source term is in this case zero
        if(ABS(lsSolver().m_outsideGValue - ABS(lsSolver().a_levelSetFunctionG(gc, 0))) < eps) {
          continue;
        }
        // compute the Echekki-Ferziger constant
        // - compute the inverse eddy viscosity (s/m^2)
        // - - density of the unburnt gas rho^u = fvSolver().m_rhoInfinity
        // - - assuming m_TInfinity is the temperature of the unburnt gas -> muInf=SUTHERLANDLAW(m_TInfinity)
        // - - DthInf = muInf^u / ( rho^u *Pr )
        FD = F1 / fvSolver().m_DthInfinity;
        // - compute Rr
        levelSetNegative = lsSolver().a_levelSetFunctionG(gc, 0) - m_noReactionCells;
        levelSetPlus = lsSolver().a_levelSetFunctionG(gc, 0) + m_noReactionCells;

        if(!m_constantFlameSpeed) {
          if(lsSolver().m_sharpDamp) {
            Rr = echekkiFerzigerPrefactor
                 * POW2(lsSolver().a_flameSpeedG(gc, 0) * (F1 - lsSolver().a_curvatureG(gc, 0) * m_marksteinLength))
                 * FD * F1B4 * (1 + tanh(levelSetPlus * 100.0)) * (1 - tanh(levelSetNegative * 100.0));
          } else {
            Rr = echekkiFerzigerPrefactor
                 * POW2(lsSolver().a_flameSpeedG(gc, 0) * (F1 - lsSolver().a_curvatureG(gc, 0) * m_marksteinLength))
                 * FD;
          }
        } else {
          if(lsSolver().m_sharpDamp) {
            Rr = echekkiFerzigerPrefactor * POW2(lsSolver().a_flameSpeedG(gc, 0)) * FD * F1B4
                 * (1 + tanh(levelSetPlus * 100.0)) * (1 - tanh(levelSetNegative * 100.0));
          } else {
            Rr = echekkiFerzigerPrefactor * POW2(lsSolver().a_flameSpeedG(gc, 0)) * FD;
          }
        }
        // damp out the reaction rate to the wall w(y=0.0341) = 0 by damping the model constant Rr(y=0.0341)=0
        if(m_heatReleaseDamp) {
          if(lsSolver().c_coordinate(gc, 1) < (m_yOffsetFlameTube + m_dampingDistanceFlameBase)) {
            if(lsSolver().c_coordinate(gc, 1) > m_yOffsetFlameTube) {
              Rr *= 1. / m_dampingDistanceFlameBase * (lsSolver().c_coordinate(gc, 1) - m_yOffsetFlameTube);
            }
          }
          if(lsSolver().c_coordinate(gc, 1) < m_yOffsetFlameTube) {
            Rr = F0;
          }
        }

        xi = lsSolver().a_levelSetFunctionG(gc, 0) * factor; //< xi = G(x,t)/(sigma*sqrt(2))

        // - compute c bar
        a = xi - sq1F2 * factor1 * sigma * FlaminarFlameThickness;
        a1 = F1B2 * (POW2(sigma * factor1 * FlaminarFlameThickness))
             - SQRT2 * xi * sigma * factor1 * FlaminarFlameThickness;
        if(a > F3) {
          a2 = (F1 - fvSolver().m_c0) * exp(a1) * F1B2
               * (-erfc(a) + F2); // erfc computes 1-erf and is more accurate for large a
        } else {
          a2 = (F1 - fvSolver().m_c0) * exp(a1) * F1B2 * (erf(a) + F1);
        }
        b = xi + sq1F2 * sigma * FlaminarFlameThickness;
        b1 = F1B2 * POW2(sigma * FlaminarFlameThickness) + SQRT2 * xi * sigma * FlaminarFlameThickness;
        if(b > F3) {
          b2 = fvSolver().m_c0 * exp(b1) * F1B2 * (-erfc(b)); // erfc computes 1-erf and is more accurate for large b
        } else {
          b2 = fvSolver().m_c0 * exp(b1) * F1B2 * (erf(b) - F1);
        }
        c1 = F1B2 * (erf(xi) - F1);
        cbar = F1 - (a2 + b2 - c1);

        // - compute Psi
        if(ABS(F1 - cbar) < eps) {
          Psi = F1;
        } else {
          Psi = a2 / ((F1 - cbar) * POW2((rhoUFrhoB + cbar * rhoJump)));
        }

        fvSolver().a_psi(c) = Psi;

        // compute rhoBar
        rhoBar = m_rhoUnburnt / (1 - fvSolver().a_pvariable(c, fvSolver().PV->C) * rhoJump);

        // compute the source term
        fvSolver().a_reactionRate(c, 0) =
            m_Re0 * rhoBar * rhoUFrhoB * Rr * (F1 - fvSolver().a_pvariable(c, fvSolver().PV->C)) * Psi;

        // catch nan reaction rate
        if(!(fvSolver().a_reactionRate(c, 0) >= F0) && !(fvSolver().a_reactionRate(c, 0) <= F0)) {
          diverged = 1;
          // m_log << "reaction rate is nan!!!" << endl;
          // m_log << "b=" << b << " " << "b1=" << b1 << " " << "exp(b1)=" << exp(b1) << " " << "b2=" << b2 << " "
          //<< "c1=" << c1 << endl;
          // m_log << "a=" << a << " " << "a1=" << a1 << " " << "a2=" << a2 << endl;
          // m_log << "g cell info for " << gc << endl;
          // m_log << "gc_coordx=" << lsSolver().c_coordinate(gc, 0) << " " << "gc_coordy=" <<
          // lsSolver().c_coordinate(gc, 1) << " glevel "
          //<< lsSolver().a_level(gc) << " inBand " << lsSolver().a_inBandG(gc, 0) << " in shadow layer " <<
          // lsSolver().a_inShadowLayerG(gc)
          //<< " no childs " << lsSolver().a_noChildrenG(gc) << endl;
          // m_log << "is GWindow Cell " << lsSolver().a_isGWindowCellG(gc) << endl;
          // m_log << "is GHalo Cell " << lsSolver().a_isGHaloCellG(gc) << endl;
          // m_log << "levelsetfunction=" << lsSolver().a_levelSetFunctionG(IDX_LSSET(gc, 0)] << " " << "xi=" << xi
          // << endl; m_log << "cbar=" << cbar << " " << "Psi=" << Psi << " " << "rhoBar=" << rhoBar << " " << "c="
          //<< fvSolver().a_pvariable(c, fvSolver().PV->C) << endl;
          // m_log << "rho=" << fvSolver().a_pvariable(c, fvSolver().PV->RHO) << endl;
          // m_log << "rhoUnburnt=" << m_rhoUnburnt << endl;
          // m_log << "global cell info for " << fvSolver().c_globalId(c) << endl;
          // m_log << "cell level " << a_level(c) << endl;
          // m_log << "window cell " << fvSolver().a_hasProperty(c, Cell::IsWindow) << endl;
          // m_log << "halo cell " << fvSolver().a_hasProperty(c, Cell::IsHalo) << endl;
          // m_log << "x=" << fvSolver().a_coordinate(c, 0) << " "
          //<< "y=" << fvSolver().a_coordinate(c, 1) << " "
          //<< "z=" << fvSolver().a_coordinate(c, 2) << endl;
          // fvSolver().a_reactionRate(c, 0) = 1000.0;
        }

        // compute the source terms
        fvSolver().a_rightHandSide(c, fvSolver().CV->RHO_C) -=
            fvSolver().a_reactionRate(c, 0) * fvSolver().a_cellVolume(c);
        fvSolver().a_rightHandSide(c, fvSolver().CV->RHO_E) -=
            fvSolver().a_reactionRate(c, 0) * fvSolver().a_cellVolume(c) * reactionEnthalpy;

        // compute the maximum reaction rate and the total heat release
        if(!fvSolver().a_hasProperty(c, Cell::IsHalo)) {
          fvSolver().m_maxReactionRate = mMax(fvSolver().a_reactionRate(c, 0), fvSolver().m_maxReactionRate);
          fvSolver().m_totalHeatReleaseRate +=
              fvSolver().a_reactionRate(c, 0) * fvSolver().a_cellVolume(c) * reactionEnthalpy;
        }
      }
      MPI_Allreduce(MPI_IN_PLACE, &diverged, 1, MPI_INT, MPI_MAX, fvSolver().mpiComm(), AT_, "MPI_IN_PLACE",
                    "diverged");
      if(diverged == 1) {
        // m_log << "Solution diverged, writing NETCDF file for debugging" << endl;
        MString errorMessage = "reaction rate is nan";
        // saveSolverSolution(1);
        // stringstream fileNameVtk;
        // stringstream levelSetFileName;
        // levelSetFileName << "restartLevelSetNewGCellGrid_" << "debug";
        //
        // stringstream levelSetFileNameSol;
        // levelSetFileNameSol << "restartLevelSetNew_" << "debug";
        //

        mTerm(1, AT_, errorMessage);
      }

      break;
    }
      // new flame speed model for turbulent flames, see Pitsch et al. 2005
    case 5401000:
    case 2751600: {
      MInt diverged = 0;
      MFloat diffTimeStep = 50000;
      if(m_temperatureChange && (globalTimeStep - m_restartTimeStep) <= diffTimeStep) {
        MFloat diffTemp =
            (fvSolver().m_burntUnburntTemperatureRatioEnd - fvSolver().m_burntUnburntTemperatureRatioStart);
        fvSolver().m_burntUnburntTemperatureRatio = (diffTemp / diffTimeStep) * (globalTimeStep - m_restartTimeStep)
                                                    + fvSolver().m_burntUnburntTemperatureRatioStart;
      }
      rhoBurnt = m_rhoFlameTube / fvSolver().m_burntUnburntTemperatureRatio;

      FrhoBurnt = fvSolver().m_burntUnburntTemperatureRatio * F1 / m_rhoFlameTube;

      rhoJump = F1 - fvSolver().m_burntUnburntTemperatureRatio;

      for(MInt c = 0; c < m_bndryGhostCellsOffset; c++) {
        if(fvSolver().grid().tree().hasProperty(c, Cell::IsHalo)) {
          continue;
        }
        if(!fvSolver().a_hasProperty(c, SolverCell::IsOnCurrentMGLevel)) {
          continue;
        }
        MInt cLs = fv2lsId(c);
        if(cLs < 0) {
          continue;
        }
        gc = cLs;

        if(gc == -1) {
          continue;
        }

        // continue if the g cell is out of the band
        // the source term is in this case zero
        if(ABS(lsSolver().m_outsideGValue - ABS(lsSolver().a_levelSetFunctionG(gc, 0))) < eps) {
          continue;
        }

        // compute the Echekki-Ferziger constant
        // - compute the inverse eddy viscosity (s/m^2)
        // - - density of the unburnt gas rho^u = fvSolver().m_rhoInfinity
        // - - assuming m_TInfinity is the temperature of the unburnt gas -> muInf=SUTHERLANDLAW(m_TInfinity)
        // - - DthInf = muInf^u / ( rho^u *Pr )
        FD = F1 / fvSolver().m_DthInfinity;
        // - compute Rr
        levelSetNegative = lsSolver().a_levelSetFunctionG(gc, 0) - m_noReactionCells;
        levelSetPlus = lsSolver().a_levelSetFunctionG(gc, 0) + m_noReactionCells;

        // turb. flame speed model, see Pitsch et al. 2005
        MFloat turbFlameSpeed = m_turbFlameSpeed;                             // once calculated
        MFloat delta = fvSolver().c_cellLengthAtLevel(fvSolver().maxLevel()); // LES filter width equals grid siz
        MFloat uAmpl = pow(m_integralAmplitude, 3);                           // filtered velocity
        uAmpl *= delta;
        uAmpl /= m_integralLengthScale;
        uAmpl = pow(uAmpl, F1B3); // filtered velocity Pitsch et al. 2005
        MFloat flameSpeed = m_flameSpeed;

        MFloat Dt = m_NuT / m_ScT;

        MFloat FDa = Dt;
        if(fvSolver().m_Da > 1) {
          FDa *= F1 / pow(fvSolver().m_Da, 2);
        }
        // flame stretch effects (Freitag et al. 2007)
        flameSpeed *= (F1
                       - lsSolver().a_curvatureG(gc, 0)
                             * m_marksteinLength); //(+=fvSolver().m_DthInfinity * lsSolver().a_curvatureG( gc , 0)   );
        // turb. flame speed
        flameSpeed += turbFlameSpeed;
        // turb. stretch effects
        flameSpeed -= FDa * lsSolver().a_curvatureG(gc, 0);
        // take the power of 2
        flameSpeed = pow(flameSpeed, 2);


        Rr = echekkiFerzigerPrefactor * flameSpeed * FD;

        // damp out the reaction rate to the wall w(y=0.0341) = 0 by damping the model constant Rr(y=0.0341)=0
        if(m_heatReleaseDamp) {
          if(lsSolver().c_coordinate(gc, 1) < (m_yOffsetFlameTube + m_dampingDistanceFlameBase)) {
            if(lsSolver().c_coordinate(gc, 1) > m_yOffsetFlameTube) {
              Rr *= 1. / m_dampingDistanceFlameBase * (lsSolver().c_coordinate(gc, 1) - m_yOffsetFlameTube);
            }
          }
          if(lsSolver().c_coordinate(gc, 1) < m_yOffsetFlameTube) {
            Rr = F0;
          }
        }

        // compute Psi
        // - compute xi

        xi = lsSolver().a_levelSetFunctionG(gc, 0) * factor; //< xi = G(x,t)/(sigma*sqrt(2))

        // - compute c bar
        a = xi - sq1F2 * factor1 * sigma * FlaminarFlameThickness;
        a1 = F1B2 * (POW2(sigma * factor1 * FlaminarFlameThickness))
             - SQRT2 * xi * sigma * factor1 * FlaminarFlameThickness;
        if(a > F3) {
          a2 = (F1 - fvSolver().m_c0) * exp(a1) * F1B2
               * (-erfc(a) + F2); // erfc computes 1-erf and is more accurate for large a
        } else {
          a2 = (F1 - fvSolver().m_c0) * exp(a1) * F1B2 * (erf(a) + F1);
        }
        b = xi + sq1F2 * sigma * FlaminarFlameThickness;
        b1 = F1B2 * POW2(sigma * FlaminarFlameThickness) + SQRT2 * xi * sigma * FlaminarFlameThickness;
        if(b > F3) {
          b2 = fvSolver().m_c0 * exp(b1) * F1B2 * (-erfc(b)); // erfc computes 1-erf and is more accurate for large b
        } else {
          b2 = fvSolver().m_c0 * exp(b1) * F1B2 * (erf(b) - F1);
        }
        c1 = F1B2 * (erf(xi) - F1);
        cbar = F1 - (a2 + b2 - c1);

        // - compute Psi
        if(ABS(F1 - cbar) < eps) {
          Psi = F1;
        } else {
          Psi = a2 / ((F1 - cbar) * POW2((rhoUFrhoB + cbar * rhoJump)));
        }

        fvSolver().a_psi(c) = Psi;

        // compute rhoBar
        rhoBar = m_rhoUnburnt / (1 - fvSolver().a_pvariable(c, fvSolver().PV->C) * rhoJump);

        // compute the source term
        fvSolver().a_reactionRate(c, 0) =
            m_Re0 * rhoBar * rhoUFrhoB * Rr * (F1 - fvSolver().a_pvariable(c, fvSolver().PV->C)) * Psi;

        // catch nan reaction rate
        if(!(fvSolver().a_reactionRate(c, 0) >= F0) && !(fvSolver().a_reactionRate(c, 0) <= F0)) {
          diverged = 1;
          // m_log << "reaction rate is nan!!!" << endl;
          // m_log << "b=" << b << " " << "b1=" << b1 << " " << "exp(b1)=" << exp(b1) << " " << "b2=" << b2 << " "
          //<< "c1=" << c1 << endl;
          // m_log << "a=" << a << " " << "a1=" << a1 << " " << "a2=" << a2 << endl;
          // m_log << "g cell info for " << gc << endl;
          // m_log << "gc_coordx=" << lsSolver().c_coordinate(gc, 0) << " " << "gc_coordy=" <<
          // lsSolver().c_coordinate(gc, 1) << " glevel "
          //<< lsSolver().a_level(gc) << " inBand " << lsSolver().a_inBandG(gc, 0) << " in shadow layer " <<
          // lsSolver().a_inShadowLayerG(gc)
          //<< " no childs " << lsSolver().a_noChildrenG(gc) << endl;
          // m_log << "is GWindow Cell " << lsSolver().a_isGWindowCellG(gc) << endl;
          // m_log << "is GHalo Cell " << lsSolver().a_isGHaloCellG(gc) << endl;
          // m_log << "levelsetfunction=" << lsSolver().a_levelSetFunctionG(IDX_LSSET(gc, 0)] << " " << "xi=" << xi
          // << endl; m_log << "cbar=" << cbar << " " << "Psi=" << Psi << " " << "rhoBar=" << rhoBar << " " << "c="
          //<< fvSolver().a_pvariable(c, fvSolver().PV->C) << endl;
          // m_log << "rho=" << fvSolver().a_pvariable(c, fvSolver().PV->RHO) << endl;
          // m_log << "rhoUnburnt=" << m_rhoUnburnt << endl;
          // m_log << "global cell info for " << fvSolver().c_globalId(c) << endl;
          // m_log << "cell level " << a_level(c) << endl;
          // m_log << "window cell " << fvSolver().a_hasProperty(c, Cell::IsWindow) << endl;
          // m_log << "halo cell " << fvSolver().a_hasProperty(c, Cell::IsHalo) << endl;
          // m_log << "x=" << fvSolver().a_coordinate(c, 0) << " "
          //<< "y=" << fvSolver().a_coordinate(c, 1) << " "
          //<< "z=" << fvSolver().a_coordinate(c, 2) << endl;
          // fvSolver().a_reactionRate(c, 0) = 1000.0;
        }

        // compute the source terms
        fvSolver().a_rightHandSide(c, fvSolver().CV->RHO_C) -=
            fvSolver().a_reactionRate(c, 0) * fvSolver().a_cellVolume(c);
        fvSolver().a_rightHandSide(c, fvSolver().CV->RHO_E) -=
            fvSolver().a_reactionRate(c, 0) * fvSolver().a_cellVolume(c) * reactionEnthalpy;

        // compute the maximum reaction rate and the total heat release
        if(!fvSolver().a_hasProperty(c, Cell::IsHalo)) {
          fvSolver().m_maxReactionRate = mMax(fvSolver().a_reactionRate(c, 0), fvSolver().m_maxReactionRate);
          fvSolver().m_totalHeatReleaseRate +=
              fvSolver().a_reactionRate(c, 0) * fvSolver().a_cellVolume(c) * reactionEnthalpy;
        }
      }
      MPI_Allreduce(MPI_IN_PLACE, &diverged, 1, MPI_INT, MPI_MAX, fvSolver().mpiComm(), AT_, "MPI_IN_PLACE",
                    "diverged");
      if(diverged == 1) {
        // m_log << "Solution diverged, writing NETCDF file for debugging" << endl;
        MString errorMessage = "reaction rate is nan";
        // saveSolverSolution(1);
        // stringstream fileNameVtk;
        mTerm(1, AT_, errorMessage);
      }
      break;
    }
    default: {
      for(MInt c = 0; c < m_bndryGhostCellsOffset; c++) {
        if(fvSolver().grid().tree().hasProperty(c, Cell::IsHalo)) {
          continue;
        }
        if(!fvSolver().a_hasProperty(c, SolverCell::IsOnCurrentMGLevel)) {
          continue;
        }
        MInt cLs = fv2lsId(c);
        if(cLs < 0) {
          continue;
        }
        gc = cLs;
        if(gc == -1) {
          continue;
        }

        // continue if the g cell is out of the band
        // the source term is in this case zero
        if(ABS(lsSolver().m_outsideGValue - ABS(lsSolver().a_levelSetFunctionG(gc, 0))) < eps) {
          continue;
        }

        // compute the Echekki-Ferziger constant
        // - compute the inverse eddy viscosity (s/m^2)
        // - - density of the unburnt gas
        rhoU = fvSolver().m_rhoInfinity;
        // - - assuming m_TInfinity is the temperature of the unburnt gas
        // - - Dth = mue^u / ( rho^u *Pr )
        Dth = SUTHERLANDLAW(m_TInfinity) / (rhoU * m_Pr);
        FD = F1 / Dth;
        // - compute Rr
        Rr = echekkiFerzigerPrefactor * POW2(lsSolver().a_correctedBurningVelocity(gc, 0)) * FD;

        // compute Psi
        // - compute xi
        xi = lsSolver().a_levelSetFunctionG(gc, 0) * factor;
        // - compute c bar
        a = xi - sq1F2 * factor1 * sigma * FlaminarFlameThickness;

        a1 = F1B2 * (POW2(sigma * factor1 * FlaminarFlameThickness))
             - SQRT2 * xi * sigma * factor1 * FlaminarFlameThickness;
        if(a > F3) {
          a2 = (F1 - fvSolver().m_c0) * exp(a1) * F1B2
               * (-erfc(a) + F2); // erfc computes 1-erf and is more accurate for large a
        } else {
          a2 = (F1 - fvSolver().m_c0) * exp(a1) * F1B2 * (erf(a) + F1);
        }
        b = xi + sq1F2 * sigma * FlaminarFlameThickness;
        b1 = F1B2 * POW2(sigma * FlaminarFlameThickness) + SQRT2 * xi * sigma * FlaminarFlameThickness;
        if(b > F3) {
          b2 = fvSolver().m_c0 * exp(b1) * F1B2 * (-erfc(b)); // erfc computes 1-erf and is more accurate for large b
        } else {
          b2 = fvSolver().m_c0 * exp(b1) * F1B2 * (erf(b) - F1);
        }
        c1 = F1B2 * (erf(xi) - F1);
        cbar = F1 - (a2 + b2 - c1);
        // - compute Psi
        denominator = (F1 - cbar) * POW2((rhoUFrhoB + cbar * rhoJump * FrhoBurnt));
        if(approx(denominator, 0.0, MFloatEps)) {
          Psi = F1;
        } else {
          Psi = a2 / denominator;
        }

        // compute rhoBar
        rhoBar = fvSolver().m_rhoInfinity / (1 - fvSolver().a_pvariable(c, fvSolver().PV->C) * rhoJump / rhoBurnt);

        // compute the source term
        fvSolver().a_reactionRate(c, 0) = m_Re0 * fvSolver().a_pvariable(c, fvSolver().PV->RHO) * rhoUFrhoB * Rr
                                          * (F1 - fvSolver().a_pvariable(c, fvSolver().PV->C)) * Psi;

        // catch nan reaction rate
        if(!(fvSolver().a_reactionRate(c, 0) >= F0) && !(fvSolver().a_reactionRate(c, 0) <= F0)) {
          // cerr << "reaction rate is nan!!!" << endl;
          // cerr << b << " " << b1 << " " << exp(b1) << " " << b2 << " " << c1 << endl;
          // cerr << a << " " << a1 << " " << a2 << endl;
          // cerr << gc << " " << lsSolver().c_coordinate(gc, 0) << " " << lsSolver().c_coordinate(gc, 1) << " " <<
          // lsSolver().c_coordinate(gc, 2)
          //<< endl;
          // cerr << lsSolver().a_levelSetFunctionG(IDX_LSSET(gc, 0)] << " " << xi << endl;
          // cerr << cbar << " " << Psi << " " << rhoBar << " " << fvSolver().a_pvariable(c, fvSolver().PV->C) << endl;
          // cerr << "cell info for " << c << endl;
          // cerr << a_level(c) << endl;
          // cerr << fvSolver().a_coordinate(c, 0) << " "
          //<< fvSolver().a_coordinate(c, 1) << " "
          //<< fvSolver().a_coordinate(c, 2) << endl;
          // fvSolver().a_reactionRate(c, 0) = 1000.0;
          // saveSolverSolution(1);
          mTerm(1, AT_, "reaction rate is nan!!!");
        }

        // compute the source terms
        fvSolver().a_rightHandSide(c, fvSolver().CV->RHO_C) -=
            fvSolver().a_reactionRate(c, 0) * fvSolver().a_cellVolume(c);
        fvSolver().a_rightHandSide(c, fvSolver().CV->RHO_E) -=
            fvSolver().a_reactionRate(c, 0) * fvSolver().a_cellVolume(c) * reactionEnthalpy;

        // compute the maximum reaction rate and the total heat release
        fvSolver().m_maxReactionRate = mMax(fvSolver().a_reactionRate(c, 0), fvSolver().m_maxReactionRate);
        fvSolver().m_totalHeatReleaseRate +=
            fvSolver().a_reactionRate(c, 0) * fvSolver().a_cellVolume(c) * reactionEnthalpy;
      }
      break;
    }
  }
  if(fvSolver().noDomains() > 1) {
    MPI_Allreduce(MPI_IN_PLACE, &fvSolver().m_totalHeatReleaseRate, 1, MPI_DOUBLE, MPI_SUM, fvSolver().mpiComm(), AT_,
                  "MPI_IN_PLACE", "fvSolver().m_totalHeatReleaseRate");
    MPI_Allreduce(MPI_IN_PLACE, &fvSolver().m_maxReactionRate, 1, MPI_DOUBLE, MPI_MAX, fvSolver().mpiComm(), AT_,
                  "MPI_IN_PLACE", "fvSolver().m_maxReactionRate");
  }
}

// ----------------------------------------------------------------------------------------

/** \brief transfers v from the flow to the G-grid (highest level) via interpolation
 *
 * \author Stephan Schlimpert
 * \date Feb 2012
 *
 * guarantees pocket formation at the flame tip
 *
 * currently out of use!
 *
 */
template <MInt nDim, class SysEqn>
void LsFvCombustion<nDim, SysEqn>::collectGEquationModelDataOptInterpolate(MFloat* fluidDensity, MInt set) {
  TRACE();

  std::ignore = fluidDensity;
  std::ignore = set;

  /*
    MInt baseGCell,flowCell,cellId;
    MFloat Fdx;
    MFloat FlengthLevel0 = F1 / grid().cellLengthAtLevel(0);
    //---

    // reset extension velocity
    //  if(globalTimeStep == m_restartTimeStep){
    for( MInt id = 0; id < lsSolver().a_noBandCells(set); id++ ){
      cellId = lsSolver().a_bandCellId(id, set);
      a_extensionVelocityG( cellId, 0 , set)  = F0;
      a_extensionVelocityG( cellId, 1 , set)  = F0;
      a_correctedBurningVelocity(cellId, set) = F0;
    }
    //}


    // take only non-ghost cells!
    for( MInt id=0; id < lsSolver().a_noG0Cells(set); id++ ) {
      if( a_isHalo( a_G0CellId(id, set))) {
        for( MInt i=0; i<nDim; i++ )
    a_extensionVelocityG( a_G0CellId(id, set), i, set)  = F0;
        fluidDensity[ gc ] = -F1; //FrhoInf;

      } else {

        // find the parent cell which is connected to the flow grid
        baseGCell = a_G0CellId(id, set);
        while( combustion().ls2fvId( baseGCell ) == -1 ) {
          baseGCell = a_parentGId( baseGCell );
        }
        if(baseGCell==-1) {
          cerr << "ERROR: no parent cell found for connecting" << endl;
        }
        flowCell = combustion().ls2fvId( baseGCell );
        if(flowCell==-1)
        cerr << "ERROR: no parent cell found for connecting" << endl;



        // compute the gradient on the flow cell (rewrite, just used for a proof of concept!!!)
        Fdx =  FlengthLevel0 * FPOW2( grid().tree().level( flowCell ) );
        for( MInt v=0; v<fvSolverD().PV->noVariables; v++ ) {
    for( MInt i=0; i<nDim; i++ ) {
      if( grid().tree().hasNeighbor( flowCell ,  2*i ) > 0 ){
        if( grid().tree().hasNeighbor( flowCell ,  2*i+1 ) > 0 ){
          fvSolverD().a_slope( flowCell ,  v ,  i ) = F1B2 * Fdx *
      (fvSolverD().a_pvariable( grid().tree().neighbor( flowCell ,  2*i+1 ) ,  v ) -
    fvSolverD().a_pvariable( grid().tree().neighbor( flowCell ,  2*i ) ,  v ) );
        } else {
          fvSolverD().a_slope( flowCell ,  v ,  i ) = Fdx *
      (fvSolverD().a_pvariable( flowCell ,  v ) -
    fvSolverD().a_pvariable( grid().tree().neighbor( flowCell ,  2*i ) ,  v ) );
        }
      }else{
        if( grid().tree().hasNeighbor( flowCell ,  2*i+1 ) > 0 ) {
          fvSolverD().a_slope( flowCell ,  v ,  i ) = Fdx *
      (fvSolverD().a_pvariable( grid().tree().neighbor( flowCell ,  2*i+1 ) ,  v ) -
    fvSolverD().a_pvariable( flowCell ,  v ) );
        } else {
          fvSolverD().a_slope( flowCell ,  v ,  i ) = F0;
          fvSolverD().a_slope( flowCell ,  v ,  i ) = F0;
        }
      }
    }
        }

        // store the flow velocity
        for( MInt i=0; i<nDim; i++ ){
    a_extensionVelocityG( a_G0CellId(id, set), i, set)  =fvSolverD().a_pvariable( flowCell ,
    fvSolverD().PV->VV[ i ] ); for( MInt j=0; j<nDim; j++ ){ a_extensionVelocityG( a_G0CellId(id, set), i
    , set)  += fvSolverD().a_slope( flowCell ,  fvSolverD().PV->VV[ j ]  ,  j ) * ( c_coordinate( a_G0CellId(
    id, set),  j ) - grid().tree().coordinate( flowCell ,  j ) );
    }
        }
        fluidDensity[ gc ] =fvSolverD().a_pvariable( flowCell ,  fvSolverD().PV->RHO );

        for( MInt i=0; i<nDim; i++ ){
    fluidDensity[ gc ] +=
      fvSolverD().a_slope( flowCell ,  fvSolverD().PV->RHO  ,  i ) * ( c_coordinate( a_G0CellId(id, set),
    i ) - grid().tree().coordinate( flowCell ,  i ) );

        }
        // store the density
        fluidDensity[ gc ] = F1 / fluidDensity[ gc ]; //fvSolverD().a_pvariable( flowCell ,  fvSolverD().PV->RHO );
      }
    }
  */
}

//-----------------------------------------------------------------------------------------


/** \brief
 *
 * \author Daniel Hartmann
 * \date 2007
 *
 * currently out of use!
 *
 */
template <MInt nDim, class SysEqn>
void LsFvCombustion<nDim, SysEqn>::fastInterfaceExtensionVelocity() {
  TRACE();

  for(MInt set = lsSolver().m_startSet; set < lsSolver().m_noSets; set++) {
    // compute the extension velocity of G0 cut cells
    for(MInt id = 0; id < lsSolver().a_bandLayer(0, set); id++) {
      MFloat correctedBurningVelocity =
          lsSolver().a_flameSpeedG(lsSolver().a_bandCellId(id, set), set)
          * (F1 - lsSolver().a_curvatureG(lsSolver().a_bandCellId(id, set), set) * lsSolver().m_marksteinLength);

      for(MInt i = 0; i < nDim; i++)
        lsSolver().a_extensionVelocityG(lsSolver().a_bandCellId(id, set), i, set) =
            fvSolver().a_variable(lsSolver().a_bandCellId(id, set), fvSolver().PV->VV[i])
            + lsSolver().a_normalVectorG(lsSolver().a_bandCellId(id, set), i, set) * correctedBurningVelocity;
    }
  }
}

template class LsFvCombustion<2, FvSysEqnNS<2>>;
template class LsFvCombustion<3, FvSysEqnNS<3>>;
template class LsFvCombustion<2, FvSysEqnRANS<2, RANSModelConstants<RANS_SA_DV>>>;
template class LsFvCombustion<3, FvSysEqnRANS<3, RANSModelConstants<RANS_SA_DV>>>;
template class LsFvCombustion<2, FvSysEqnRANS<2, RANSModelConstants<RANS_FS>>>;
template class LsFvCombustion<3, FvSysEqnRANS<3, RANSModelConstants<RANS_FS>>>;
template class LsFvCombustion<2, FvSysEqnRANS<2, RANSModelConstants<RANS_KOMEGA>>>;
template class LsFvCombustion<3, FvSysEqnRANS<3, RANSModelConstants<RANS_KOMEGA>>>;
