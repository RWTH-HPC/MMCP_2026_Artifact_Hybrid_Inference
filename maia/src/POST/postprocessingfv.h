// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef POSTPROCESSINGFV_H_
#define POSTPROCESSINGFV_H_

#include <array>
#include <map>
#include <set>
#include <vector>

#include "globals.h"
#include "postdata.h"
#include "postprocessing.h"
#include "samplingdata.h"

#include "FV/fvcartesiansyseqntraits.h"

template <MInt nDim, class SysEqn>
class FvCartesianSolverXD;

template <MInt nDim, class SysEqn>
class PostProcessingFv : public PostProcessing<nDim, PostProcessingFv<nDim, SysEqn>> {
 private:
  template <MInt nDim_, class ppType>
  friend class PostProcessing;

 public:
  using SolverType = FvCartesianSolverXD<nDim, SysEqn>;

  using Base = PostProcessing<nDim, PostProcessingFv<nDim, SysEqn>>;
  using Base::isActive;
  using Base::m_averageInterval;
  using Base::m_averageStartTimestep;
  using Base::m_averageStopTimestep;
  using Base::m_averageVorticity;
  using Base::m_finalTimeStep;
  using Base::m_globalNoProbeLineAverageIds;
  using Base::m_globalNoProbeLineIds;
  using Base::m_globalProbeLineAverageVars;
  using Base::m_movAvgNoVariables;
  using Base::m_movAvgVariables;
  using Base::m_movAvgVarNames;
  using Base::m_movingAverageCounter;
  using Base::m_movingAverageDataPoints;
  using Base::m_movingAverageInterval;
  using Base::m_noProbeLineAverageIds;
  using Base::m_noProbeLineIds;
  using Base::m_noProbeLines;
  using Base::m_noVariables;
  using Base::m_postData;
  using Base::m_postprocessFileName;
  using Base::m_postprocessingId;
  using Base::m_probeLineAverageCoordinates;
  using Base::m_probeLineAverageDirection;
  using Base::m_probeLineAverageIds;
  using Base::m_probeLineDirection;
  using Base::m_probeLineIds;
  using Base::m_probeLineOffsets;
  using Base::m_probeLinePositions;
  using Base::m_sprayComputeInterval;
  using Base::m_sprayDataSize;
  using Base::m_sprayDataStep;
  using Base::m_sprayWriteInterval;
  using Base::m_vapourCV;
  using Base::m_vapourPen;
  using Base::postData;

  // Constructor
  PostProcessingFv(MInt postprocessingId_, PostData<nDim>* data, SolverType* ppSolver_);

  virtual ~PostProcessingFv();

  void initPostProcessing() override;

  SolverType& solver() const { return *m_ppSolver; }

 private:
  SolverType* m_ppSolver;

 protected:
  void probeLinePeriodicPost() override;
  void probeLinePeriodic() override;
  void initMovingAverage() override;

  void initAveragingProperties() override;

  void initPointSamplingData() override;
  void savePointSamplingData() override;

  void initSurfaceSamplingData() override;
  void saveSurfaceSamplingData() override;

  void initVolumeSamplingData() override;
  void saveVolumeSamplingData() override;

  void initSprayData() override;

  // TODO labels:FV Has to be overriden for FV because it calculates vorticity = 0.5 * nabla x u
  void calcVorticity(const MFloatTensor& deriv, MFloat vorticity[nDim * 2 - 3]) override {
    if constexpr(nDim == 2) {
      vorticity[0] = 0.5 * (deriv(1, 0) - deriv(0, 1));
    } else { //(nDim ==3)
      vorticity[0] = 0.5 * (deriv(2, 1) - deriv(1, 2));
      vorticity[1] = 0.5 * (deriv(0, 2) - deriv(2, 0));
      vorticity[2] = 0.5 * (deriv(1, 0) - deriv(0, 1));
    }
  }
  void getVorticity(MFloat* const vorticity) override { solver().getVorticity(&vorticity[0]); };
  void getVorticityT(MFloat* const vorticity) override { solver().getVorticityT(&vorticity[0]); };
  void getSampleVarsDerivatives(MInt cellId, const MFloat*& vars) { solver().getSampleVarsDerivatives(cellId, vars); };
  MBool getSampleVarsDerivatives(const MInt cellId, std::vector<MFloat>& vars) {
    return solver().getSampleVarsDerivatives(cellId, vars);
  };
  MFloat& vorticityAtCell(const MInt cellId, const MInt dir) override { return solver().vorticityAtCell(cellId, dir); };
  void getPrimitiveVariables(MInt cellId, MFloat* Xp, MFloat* vars, MInt order) override {
    solver().getPrimitiveVariables(cellId, Xp, vars, order);
  };

  void computeAcousticSourceTermQe(MFloatScratchSpace& QeI, MFloatScratchSpace& QeIII, MFloatScratchSpace& cSquared,
                                   MFloatScratchSpace& drhodt) override {
    solver().computeAcousticSourceTermQe(QeI, QeIII, cSquared, drhodt);
  };

  void getPrimitiveVarNames(MString* names) const override { solver().sysEqn().PV->getPrimitiveVariableNames(names); }
  MFloat getBoundaryHeatFlux(const MInt cellId) const override { return solver().getBoundaryHeatFlux(cellId); };
  void vapourPenetration(MFloat spawnCoord[nDim]);
  void vapourMass(const MInt);

  void advanceDataStep() { m_sprayDataStep++; };
  void resetDataStep() { m_sprayDataStep = 0; };

  std::unique_ptr<PointData<nDim, SolverType>> m_pointData;
  std::unique_ptr<SurfaceData<nDim, SolverType>> m_surfaceData;
  std::unique_ptr<VolumeData<nDim, SolverType>> m_volumeData;
};

#endif // POSTPROCESSINGFV_H_
