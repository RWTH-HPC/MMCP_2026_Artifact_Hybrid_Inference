// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef STRUCTUREDTIMERS_H_
#define STRUCTUREDTIMERS_H_

namespace maia {
namespace structured {

// Create struct for easy timer identification
struct Timers_ {
  // Enum to store timer "names"
  enum {
    Structured,
    Run,
    RunInit,
    MainLoop,
    Constructor,
    GridDecomposition,
    GridReading,
    Init,
    LoadRestart,
    LoadVariables,
    LoadSponge,
    LoadSTG,
    BuildUpSponge,
    ComputeMetrics,
    ComputeJacobian,
    ConvectiveFlux,
    ViscousFlux,
    SandpaperTrip,
    MovingGrid,
    MGVolumeFlux,
    MGMoveGrid,
    MGExchange,
    MGSaveGrid,
    MGCellCenterCoordinates,
    MGMetrics,
    MGSurfaceMetrics,
    MGCellMetrics,
    MGCornerMetrics,
    MGJacobian,
    WaveSpanwiseReordering,
    Exchange,
    Gather,
    Send,
    SendWait,
    Receive,
    ReceiveWait,
    Scatter,
    BoundaryCondition,
    RungeKutta,
    SaveOutput,
    SaveSolution,
    SaveForces,
    SavePlanes,
    SaveAuxdata,
    SaveBoxes,
    SaveIntpPoints,
    SetTimeStep,
    UpdateSponge,

    // Special enum value used to initialize timer array
    _count
  };
};

} // namespace structured
} // namespace maia

#endif // ifndef STRUCTUREDTIMERS_H_
