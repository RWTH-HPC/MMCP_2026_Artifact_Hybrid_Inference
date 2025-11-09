// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef FVZONAL_H_
#define FVZONAL_H_

#include <vector>
#include "FV/fvcartesiansolverxd.h"
#include "FV/fvcartesiansyseqnns.h"
#include "FV/fvcartesiansyseqnrans.h"
#include "INCLUDE/maiaconstants.h"
#include "INCLUDE/maiatypes.h"
#include "MEMORY/scratch.h"
#include "UTIL/functions.h"
#include "coupling.h"
#include "couplingutils.h"

template <MInt nDim, class SysEqn>
class FvZonal : public CouplingFv<nDim, SysEqn>, public CouplingFv<nDim, FvSysEqnNS<nDim>> {
 private:
  friend class CouplingFv<nDim, FvSysEqnNS<nDim>>;
  friend class CouplingFv<nDim, SysEqn>;

 public:
  using RANS = FvCartesianSolverXD<nDim, SysEqn>;
  using LES = FvCartesianSolverXD<nDim, FvSysEqnNS<nDim>>;

 public:
  // Constructor
  FvZonal(const MInt couplingId, RANS* R, LES* L);
  virtual ~FvZonal() = default;

  // empty main functions
  void init() override{};
  void finalizeSubCoupleInit(MInt){};
  void postCouple(MInt) override{};
  void cleanUp(){};
  void checkProperties() override{};

 protected:
  RANS& RANSSolver() const { return *CouplingFv<nDim, SysEqn>::m_fvSolvers[0]; }
  LES& LESSolver() const { return *CouplingFv<nDim, FvSysEqnNS<nDim>>::m_fvSolvers[0]; }

  MInt noExchangeVariables() { return nDim + 2; };
  MInt noLESVariables() { return LESSolver().noVariables(); };
  MInt noRANSVariables() { return RANSSolver().noVariables(); };
  MFloat getAveragingFactor() {
    return F1B8 * LESSolver().m_Ma * sqrt(LESSolver().m_TInfinity) * LESSolver().timeStep();
  };

  MInt a_noFvCellsLES() const { return LESSolver().a_noCells(); }
  MInt a_noFvGridCellsLES() const { return LESSolver().c_noCells(); }
  MInt a_noFvCellsRANS() const { return RANSSolver().a_noCells(); }
  MInt a_noFvGridCellsRANS() const { return RANSSolver().c_noCells(); }

  void initRANSValues();
  void initLESValues();

  MInt m_RANSSolverId;
  MInt m_LESSolverId;

  const MInt m_noReconstructNutVars = 6;

  MInt m_zonalAveragingTimeStep;
  MInt m_zonalTransferInterval;
  MInt m_averageTimeSteps;
  MBool m_restartLESAverage;

  MInt m_LESNoVarAverage;

  MBool m_cylindricCommunication;
  MFloat m_azimuthalAngle;

  MBool m_STGSponge = false;
};


#endif // ifndef FVZONAL_H_
