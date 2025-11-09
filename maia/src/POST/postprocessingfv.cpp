// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "postprocessingfv.h"
//#include "postprocessing.h"
#include "postprocessing.cpp"

using namespace std;

template <MInt nDim, class SysEqn>
PostProcessingFv<nDim, SysEqn>::PostProcessingFv(MInt postprocessingId_,
                                                 PostData<nDim>* data,
                                                 FvCartesianSolverXD<nDim, SysEqn>* ppSolver_)
  : PostProcessingInterface(postprocessingId_),
    PostProcessing<nDim, PostProcessingFv<nDim, SysEqn>>(postprocessingId_, data) {
  m_ppSolver = ppSolver_;
}


template <MInt nDim, class SysEqn>
PostProcessingFv<nDim, SysEqn>::~PostProcessingFv() {}

template <MInt nDim, class SysEqn>
void PostProcessingFv<nDim, SysEqn>::initPostProcessing() {
  TRACE();

  PostProcessing<nDim, PostProcessingFv<nDim, SysEqn>>::initPostProcessing();

  solver().m_skewness = this->m_skewness;
  solver().m_kurtosis = this->m_kurtosis;
  solver().m_activeMeanVars = this->m_activeMeanVars;

  // for(MInt i = 0; i < (signed)this->m_activeMeanVars.size(); i++){
  //   solver().m_activeMeanVars.insert(this->m_activeMeanVars[i]);
  // }
}


template <MInt nDim, class SysEqn>
void PostProcessingFv<nDim, SysEqn>::initPointSamplingData() {
  TRACE();
  m_pointData.reset(new PointData<nDim, SolverType>{*m_ppSolver});
  m_pointData->setInputOutputProperties();
  m_pointData->init();
}

template <MInt nDim, class SysEqn>
void PostProcessingFv<nDim, SysEqn>::savePointSamplingData() {
  if(!isActive()) return;
  m_pointData->save(m_finalTimeStep);
}

template <MInt nDim, class SysEqn>
void PostProcessingFv<nDim, SysEqn>::initSurfaceSamplingData() {
  TRACE();
  m_surfaceData.reset(new SurfaceData<nDim, SolverType>{*m_ppSolver});
  m_surfaceData->setInputOutputProperties();
  m_surfaceData->init();
}

template <MInt nDim, class SysEqn>
void PostProcessingFv<nDim, SysEqn>::saveSurfaceSamplingData() {
  if(!isActive()) return;
  m_surfaceData->save(m_finalTimeStep);
}

template <MInt nDim, class SysEqn>
void PostProcessingFv<nDim, SysEqn>::initVolumeSamplingData() {
  TRACE();
  m_volumeData.reset(new VolumeData<nDim, SolverType>{*m_ppSolver});
  m_volumeData->setInputOutputProperties();
  m_volumeData->init();
}

template <MInt nDim, class SysEqn>
void PostProcessingFv<nDim, SysEqn>::saveVolumeSamplingData() {
  if(!isActive()) return;
  m_volumeData->save(m_finalTimeStep);
}

template <MInt nDim, class SysEqn>
void PostProcessingFv<nDim, SysEqn>::initMovingAverage() {
  TRACE();

  PostProcessing<nDim, PostProcessingFv<nDim, SysEqn>>::initMovingAverage();

  solver().m_movingAvgInterval = this->m_movingAverageInterval;
}

template <MInt nDim, class SysEqn>
void PostProcessingFv<nDim, SysEqn>::initAveragingProperties() {
  TRACE();

  PostProcessing<nDim, PostProcessingFv<nDim, SysEqn>>::initAveragingProperties();

  solver().m_averageVorticity = this->m_averageVorticity;
  solver().m_averageSpeedOfSound = this->m_averageSpeedOfSound;
}


template <MInt nDim, class SysEqn>
void PostProcessingFv<nDim, SysEqn>::probeLinePeriodicPost() {
  TRACE();

  if(!solver().grid().isActive()) return;

  if(m_postprocessFileName != "") {
    m_log << "    ^        * probe line for file " << m_postprocessFileName << endl;
    // TERMM(1, "FIXME untested");
    IF_CONSTEXPR(SysEqn::m_noRansEquations == 0) postData().loadMeanFile(m_postprocessFileName);
    probeLinePeriodic();
  } else {
    for(MInt t = m_averageStartTimestep; t <= m_averageStopTimestep; t += m_averageInterval) {
      solver().loadSampleVariables(t);
      probeLinePeriodic();
    }
  }
}

template <MInt nDim, class SysEqn>
void PostProcessingFv<nDim, SysEqn>::probeLinePeriodic() {
  TRACE();

  using namespace maia::parallel_io;

  // check for average timestep or postprocessFileName
  if((globalTimeStep >= m_averageStartTimestep && (globalTimeStep - m_averageStartTimestep) % m_averageInterval == 0
      && globalTimeStep <= m_averageStopTimestep)
     || globalTimeStep == 0) {
    MInt noVars = m_noVariables;
    if(m_postData->isMeanFile()) {
      noVars = m_postData->fileNoVars();
    }

    MInt step = (m_postData->isMeanFile()) ? 0 : globalTimeStep;
    stringstream fileName;
    fileName << solver().outputDir() << "probeLines_" << step << ParallelIo::fileExt();
    ParallelIo parallelIo(fileName.str(), maia::parallel_io::PIO_REPLACE, solver().mpiComm());

    // define all arrays in output file
    for(MInt probeLineId = 0; probeLineId < m_noProbeLines; probeLineId++) {
      stringstream varNameBase;
      varNameBase << "line_" << probeLineId;
      string coordName = varNameBase.str() + "_coordinates";

      parallelIo.defineArray(PIO_FLOAT, coordName, m_globalNoProbeLineIds[probeLineId]);
      parallelIo.setAttribute(probeLineId, "lineId", coordName);

      varNameBase << "_var_";
      for(MInt varId = 0; varId < noVars; varId++) {
        stringstream varName;
        varName << varNameBase.str() << varId;
        parallelIo.defineArray(PIO_FLOAT, varName.str(), m_globalNoProbeLineIds[probeLineId]);
        parallelIo.setAttribute(probeLineId, "lineId", varName.str());
      }
    }

    for(MInt probeLineId = 0; probeLineId < m_noProbeLines; probeLineId++) {
      if(m_probeLineDirection[probeLineId] < 0 || m_probeLineDirection[probeLineId] >= nDim) {
        continue;
      }

      m_log << "    ^        * probe line timestep " << globalTimeStep << " for line #" << probeLineId << endl;

      MInt noIds = m_noProbeLineIds[probeLineId];
      if(noIds == 0) {
        noIds = 1;
      } // avoid dereferencing array with length 0 in writeArray(...)
      ScratchSpace<MFloat> vars(noVars * noIds, "vars", FUN_);

      MInt probeId;
      const MFloat* cellVars = 0;

      // collect local variables
      for(MInt i = 0; i < m_noProbeLineIds[probeLineId]; i++) {
        probeId = m_probeLineIds[probeLineId][i];
        PostProcessing<nDim, PostProcessingFv<nDim, SysEqn>>::getSampleVariables(probeId, cellVars, false);
        for(MInt varId = 0; varId < noVars; varId++) {
          vars[i * noVars + varId] = cellVars[varId];
        }
      }

      parallelIo.setOffset(m_noProbeLineIds[probeLineId], m_probeLineOffsets[probeLineId]);

      // write to file
      stringstream varNameBase;
      varNameBase << "line_" << probeLineId;
      string coordName = varNameBase.str() + "_coordinates";
      parallelIo.writeArray(m_probeLinePositions[probeLineId], coordName);

      varNameBase << "_var_";
      for(MInt varId = 0; varId < noVars; varId++) { // write all variables to file
        stringstream varName;
        varName << varNameBase.str() << varId;
        parallelIo.writeArray(&(vars[varId]), varName.str(), noVars);
      }
    }
  } else {
    // PROBE_LINE_POST
    //--------------------------------------------------------------------------
    //--------------------------------------------------------------------------

    for(MInt probeLineId = 0; probeLineId < m_noProbeLines; probeLineId++) {
      if(m_probeLineAverageDirection[probeLineId] < 0 || m_probeLineAverageDirection[probeLineId] >= nDim) {
        continue;
      }

      MInt probeId;

      if(m_probeLineDirection[probeLineId] == 1 || m_probeLineDirection[probeLineId] == 0) {
        MInt noVars;
        if(solver().m_rans) {
          noVars = m_noVariables + 1; // k
        } else {
          noVars = postData().noVariables();
        }

        MInt noIds = m_globalNoProbeLineIds[probeLineId];
        if(noIds == 0) {
          noIds = 1;
        }

        // ScratchSpace<MFloat> vars( noVars*ids, "vars", FUN_ );
        ScratchSpace<MFloat> avgVars(noVars * noIds, "avgVars", FUN_);
        ScratchSpace<MFloat> globalAvgVars(noVars * noIds, "globalAvgVars", FUN_);
        avgVars.fill(F0);
        globalAvgVars.fill(F0);

        if(solver().m_rans) {
          // RANS SOLVER
          for(MInt p = 0; p < m_globalNoProbeLineIds[probeLineId]; p++) {
            // collect local variables
            for(MInt i = 0; i < m_noProbeLineAverageIds[probeLineId][p]; i++) {
              probeId = m_probeLineAverageIds[probeLineId][p][i];
              const MFloat* cellVars = 0;
              PostProcessing<nDim, PostProcessingFv<nDim, SysEqn>>::getSampleVariables(probeId, cellVars, false);

              // RANS Variables
              for(MInt varId = 0; varId < m_noVariables; varId++) {
                avgVars[p * noVars + varId] += cellVars[varId];
              }

              // turbulent kinetic energy
              MFloat k = F0;
              const MFloat* cellVarsDeriv1;
              getSampleVarsDerivatives(probeId, cellVarsDeriv1);
              const MFloatTensor deriv1(const_cast<MFloat*>(cellVarsDeriv1), m_noVariables, nDim);
              MFloat SijSij = F0;
              for(MInt d1 = 0; d1 < nDim; d1++) {
                for(MInt d2 = 0; d2 < nDim; d2++) {
                  MFloat sij = 0.5 * (deriv1(d1, d2) + deriv1(d2, d1));
                  SijSij += sij * sij;
                }
              }

              k = sqrt(2.0 * SijSij / 0.09) * cellVars[5] / solver().sysEqn().m_Re0;

              avgVars[p * noVars + m_noVariables] += k;
            }
          }
        } else {
          // LES SOLVER
          for(MInt p = 0; p < m_globalNoProbeLineIds[probeLineId]; p++) {
            // collect local variables
            for(MInt i = 0; i < m_noProbeLineAverageIds[probeLineId][p]; i++) {
              // const MFloat* cellVars = 0;
              probeId = m_probeLineAverageIds[probeLineId][p][i];
              // PostProcessingBlock<nDim, Block>::getSampleVariables( probeId, cellVars );
              MInt dataId = convertIdParent(solver(), postData(), probeId);
              if(dataId != -1) {
                for(MInt varId = 0; varId < noVars; varId++) {
                  avgVars[p * noVars + varId] += postData().m_averagedVars[dataId][varId];
                }

                // SijSij
                // calculating strain rate tensor Sij = 0.5(dui/dxj + duj/dxi)
                std::vector<std::vector<MFloat>> du(nDim, std::vector<MFloat>(nDim, F0));
                const MInt recData = solver().a_reconstructionData(probeId);
                std::vector<MFloat> u{postData().m_averagedVars[dataId][0], postData().m_averagedVars[dataId][1],
                                      postData().m_averagedVars[dataId][2]};

                for(MInt nghbr = 0; nghbr < solver().a_noReconstructionNeighbors(probeId); nghbr++) {
                  const MInt recNghbrId = solver().a_reconstructionNeighborId(probeId, nghbr);
                  if(recNghbrId > -1 && recNghbrId < solver().noInternalCells()) {
                    MInt dataNgbhrId = convertIdParent(solver(), postData(), recNghbrId);
                    if(dataNgbhrId > -1) {
                      // std::ignore = recData;
                      // std::ignore = u;
                      const MFloat recConst_x = solver().m_reconstructionConstants[nDim * (recData + nghbr) + 0];
                      const MFloat recConst_y = solver().m_reconstructionConstants[nDim * (recData + nghbr) + 1];
                      const MFloat recConst_z = solver().m_reconstructionConstants[nDim * (recData + nghbr) + 2];
                      for(MInt dim = 0; dim < nDim; ++dim) {
                        MFloat delta_u = postData().m_averagedVars[dataNgbhrId][dim] - u[dim];
                        du[dim][0] += recConst_x * delta_u;
                        du[dim][1] += recConst_y * delta_u;
                        du[dim][2] += recConst_z * delta_u;
                      }
                    }
                  }
                }
                std::vector<std::vector<MFloat>> sij(nDim, std::vector<MFloat>(nDim, F0));
                MFloat SijSij = F0;
                for(MInt d1 = 0; d1 < nDim; d1++) {
                  for(MInt d2 = 0; d2 < nDim; d2++) {
                    sij[d1][d2] = 0.5 * (du[d1][d2] + du[d2][d1]);
                    SijSij += sij[d1][d2] * sij[d1][d2];
                  }
                }
                avgVars[p * noVars + (noVars - 1)] += SijSij;
              }
            }
          }
        }

        MPI_Allreduce(&avgVars[0], &globalAvgVars[0], noVars * noIds, MPI_DOUBLE, MPI_SUM, solver().mpiComm(), AT_,
                      "avgVars", "globalAvgVars");


        if(solver().domainId() == 0) {
          // Calculate spanwise average
          //--------------------------------------------------------------------
          MInt tempId = 0;
          for(MInt p = 0; p < m_globalNoProbeLineIds[probeLineId]; p++) {
            for(MInt varId = 0; varId < noVars; varId++) {
              m_globalProbeLineAverageVars[probeLineId][p][varId] =
                  globalAvgVars[tempId] / m_globalNoProbeLineAverageIds[probeLineId][p];
              tempId++;
            }
          }

          // vector<MInt> outputIds;
          // if(solver().m_rans) {
          //   outputIds = {0, 1, 2, 3, 4, 5, 6};
          // } else {
          //   outputIds = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
          // }
          //--------------------------------------------------------------------


          // OUTPUT
          //--------------------------------------------------------------------
          // Sorting
          ScratchSpace<MFloat> sorting(m_globalNoProbeLineIds[probeLineId], "sorting", FUN_);
          vector<int> sortingIndex(m_globalNoProbeLineIds[probeLineId]);
          for(MInt p = 0; p < m_globalNoProbeLineIds[probeLineId]; p++) {
            sorting[p] = m_probeLineAverageCoordinates[probeLineId][p];
          }
          int x = 0;
          iota(sortingIndex.begin(), sortingIndex.end(), x++); // Initializing
          sort(sortingIndex.begin(), sortingIndex.end(), [&](int i, int j) { return sorting[i] < sorting[j]; });

          // Write text file
          MString dir;
          if(solver().m_rans) {
            dir = "#y um vm wm rhom pm num ym k";
          } else {
            if(solver().m_solverId == 1) {
              dir = "#y um vm wm rhom pm num ym u'u' u'v' u'w' v'v' v'w' w'w' p'";
            } else {
              dir = "#y um vm wm rhom pm u'u' u'v' u'w' v'v' v'w' w'w' p' SijSij";
            }
          }
          ofstream lineprob;
          lineprob.precision(8);
          MString fname = "probeLine" + to_string(solver().m_solverId) + "_" + to_string(probeLineId) + "_"
                          + to_string(globalTimeStep) + ".txt";
          cerr << "Writing " << fname << endl;
          lineprob.open(fname);
          for(MInt p = 0; p < m_globalNoProbeLineIds[probeLineId]; p++) {
            MInt index = sortingIndex[p];
            MString line = "";
            line.append(to_string(m_probeLineAverageCoordinates[probeLineId][index]));
            if(p == 0) lineprob << dir << endl;
            for(MInt k = 0; k < noVars; k++) {
              MInt varId = k; // outputIds[k];
              line.append(" " + to_string(m_globalProbeLineAverageVars[probeLineId][index][varId]));
            }
            lineprob << line << endl;
          }
          lineprob.close();
          //--------------------------------------------------------------------
        }
      } else if(m_probeLineDirection[probeLineId] == 2) {
        //
      }
    }
  }
}

template <MInt nDim, class SysEqn>
void PostProcessingFv<nDim, SysEqn>::initSprayData() {
  TRACE();

  PostProcessing<nDim, PostProcessingFv<nDim, SysEqn>>::initSprayData();

  if(solver().domainId() == 0) {
    cerr << "Allocating Post-processing Fv data!" << endl;
  }

  const MInt vapourPenSize = nDim == 3 ? 16 : 12;
  mAlloc(m_vapourPen, m_sprayDataSize, vapourPenSize, "m_vapourPen", F0, AT_);
  mAlloc(m_vapourCV, m_sprayDataSize, 5 + 2 * nDim, "m_vapourCV", F0, AT_);
}

template <MInt nDim, class SysEqn>
void PostProcessingFv<nDim, SysEqn>::vapourMass(const MInt speciesId) {
  TRACE();

  // calculate dimensionless vapour information from finite volume solver
  const MBool hasSpecies = solver().m_noSpecies > 0 ? true : false;

  // reset
  for(MInt i = 0; i < 5 + 2 * nDim; i++) {
    m_vapourCV[m_sprayDataStep][i] = 0;
  }

  if(solver().isActive()) {
    for(MInt cellId = 0; cellId < solver().noInternalCells(); cellId++) {
      if(!solver().c_isLeafCell(cellId)) continue;
      if(solver().a_isInactive(cellId)) continue;
      const MFloat cellV = 1 / solver().a_FcellVolume(cellId);
      const MFloat rho = solver().a_variable(cellId, solver().m_sysEqn.CV->RHO);
      const MFloat rhoY = hasSpecies ? solver().a_variable(cellId, solver().m_sysEqn.CV->RHO_Y[speciesId]) : 0;

      m_vapourCV[m_sprayDataStep][0] += cellV * rhoY;
      m_vapourCV[m_sprayDataStep][1] += cellV * (rho - rhoY);
      m_vapourCV[m_sprayDataStep][2] += cellV * solver().a_variable(cellId, solver().m_sysEqn.CV->RHO_E);
      m_vapourCV[m_sprayDataStep][3] += cellV * solver().a_variable(cellId, solver().m_sysEqn.CV->RHO_U);
      m_vapourCV[m_sprayDataStep][4] += cellV * solver().a_variable(cellId, solver().m_sysEqn.CV->RHO_V);
      IF_CONSTEXPR(nDim == 3) {
        m_vapourCV[m_sprayDataStep][5] += cellV * solver().a_variable(cellId, solver().m_sysEqn.CV->RHO_W);
      }
    }
    if(solver().m_hasExternalSource) {
      MInt numVars = 2 + nDim;
      for(MInt i = 0; i < numVars; i++) {
        MFloat source = (solver().m_vapourData.find(globalTimeStep)->second)[i];
        m_vapourCV[m_sprayDataStep][3 + nDim + i] = source;
      }
    }
  }
}

template <MInt nDim, class SysEqn>
void PostProcessingFv<nDim, SysEqn>::vapourPenetration(MFloat spawnCoord[nDim]) {
  TRACE();

  const MInt noYLimits = 3;
  const MInt refLvl = 10;
  const MFloat yLimitRefLvl = 0.001;

  // concentration limits for differecent penetration lengths
  array<MFloat, noYLimits> yLimits{};
  yLimits[0] = 0.001;
  yLimits[1] = 0.01;
  yLimits[2] = 0.05;

  // x,y,z, axial penetration and distance
  array<array<MFloat, noYLimits>, nDim + 1> maxVapPen{};
  array<MFloat, nDim> maxVapPenRefLvl{};
  for(MInt n = 0; n < noYLimits; n++) {
    for(MInt i = 0; i < nDim + 1; i++) {
      maxVapPen[i][n] = 0.0;
    }
  }
  for(MInt i = 0; i < nDim; i++) {
    maxVapPenRefLvl[i] = 0.0;
  }

  auto assignLarger = [&](MFloat& A, MFloat b) {
    if(b > A) {
      A = b;
    }
  };

  const MBool hasSpecies = solver().m_noSpecies > 0 ? true : false;

  if(solver().isActive() && hasSpecies) {
    for(MInt cellId = 0; cellId < solver().noInternalCells(); cellId++) {
      if(solver().a_isInactive(cellId)) continue;
      if(solver().c_isLeafCell(cellId)) {
        const MFloat rho = solver().a_variable(cellId, nDim + 1);
        const MFloat Y = hasSpecies ? solver().a_variable(cellId, nDim + 2) / rho : 0;
        array<MFloat, nDim> centerCoords{};

        for(MInt n = 0; n < noYLimits; n++) {
          if(Y > yLimits[n]) {
            for(MInt i = 0; i < nDim; i++) {
              centerCoords[i] = solver().a_coordinate(cellId, i);
            }
            array<MFloat, nDim + 1> distance{};
            for(int i = 0; i < nDim; i++) {
              distance[i] = abs(centerCoords[i] - spawnCoord[i]);
            }
            distance[nDim] = maia::math::distance(centerCoords, spawnCoord);

            for(int i = 0; i < nDim + 1; i++) {
              assignLarger(maxVapPen[i][n], distance[i]);
            }
          } else {
            break;
          }
        }
      } else {
        if(solver().a_level(cellId) == refLvl) {
          solver().reduceData(cellId, &solver().a_variable(0, 0), solver().CV->noVariables);
          const MFloat rho = solver().a_variable(cellId, nDim + 1);
          const MFloat Y = hasSpecies ? solver().a_variable(cellId, nDim + 2) / rho : 0;
          if(Y > yLimitRefLvl) {
            for(MInt i = 0; i < nDim; i++) {
              assignLarger(maxVapPenRefLvl[i], solver().a_coordinate(cellId, i));
            }
          }
        }
      }
    }
  }

  // save data:
  MInt it = 0;
  for(MInt n = 0; n < noYLimits; n++) {
    // absolute distance
    m_vapourPen[m_sprayDataStep][it++] = maxVapPen[nDim][n];
    // x/y/z-distance
    for(MInt i = 0; i < nDim; i++) {
      m_vapourPen[m_sprayDataStep][it++] = maxVapPen[i][n];
    }
  }
  for(MInt i = 0; i < nDim; i++) {
    m_vapourPen[m_sprayDataStep][it++] = maxVapPenRefLvl[i];
  }
}
