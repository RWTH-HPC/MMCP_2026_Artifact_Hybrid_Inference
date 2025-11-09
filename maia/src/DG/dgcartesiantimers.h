// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef DGTIMERS_H_
#define DGTIMERS_H_

namespace maia {
namespace dg {

// Create struct for easy timer identification
struct Timers_ {
  // Enum to store timer "names"
  enum {
    SolverType,
    Constructor,
    Run,
    RunInit,
    InitSolverObjects,
    InitMortarProjection,
    InitialCondition,
    InitData,
    InitMainLoop,
    MainLoop,
    CalcTimeStep,
    ResetExternalSources,
    AdaptiveRefinement,
    RungeKuttaStep,
    TimeDeriv,
    ResetRHS,
    Coupling,
    Prolong,
    ForwardProjection,
    SurfExchange,
    SurfExchangeComm,
    SECommSend,
    SECommRecv,
    SurfExchangeCopy,
    SECopySend,
    SECopyRecv,
    SurfExchangeWait,
    SEWaitSend,
    SEWaitRecv,
    VolInt,
    Flux,
    FluxBndry,
    FluxInner,
    FluxMPI,
    SurfInt,
    Jacobian,
    Sources,
    ExternalSources,
    Sponge,
    TimeInt,
    MainLoopIO,
    Analysis,
    CleanUp,
    Destructor,
    Accumulated,
    IO,
    SaveSolutionFile,
    SaveRestartFile,
    AnalyzeTimeStep,
    MPI,
    MPIComm,
    MPICopy,
    MPIWait,

    // Special enum value used to initialize timer array
    _count
  };
};

} // namespace dg
} // namespace maia

#endif // ifndef DGTIMERS_H_
