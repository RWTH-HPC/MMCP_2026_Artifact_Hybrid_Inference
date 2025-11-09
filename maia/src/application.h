// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef APPLICATION_H
#define APPLICATION_H

#include <memory>
#include <vector>
#include "INCLUDE/maiatypes.h"
#include "solver.h"

template <MInt nDim>
class LPT;
class ExecutionRecipe;
class Solver;
template <MInt nDim>
class CartesianGrid;
class Coupling;
class PostProcessingInterface;

template <MInt nDim>
class Geometry;

template <MInt nDim>
class StructuredGrid;

/** \class Application
 *
 * \brief Manages the initialisation of the solvers, the methods
 * depending to the solvers and the start of the solving steps
 *
 *
 */
class Application {
 public:
  Application();
  ~Application();

  template <MInt nDim>
  void run();

  Application(const Application&) = delete;
  Application(Application&&) = delete;
  Application& operator=(const Application&) = delete;
  Application& operator=(Application&&) = delete;

 private:
  MInt gridType(const MString solverType);

  MInt getInitialTimeStep(const MBool restartFile, const MPI_Comm comm);
  void initTimings();
  void collectTimingsAndSolverInformation(const MBool finalTimeStep);
  void storeTimingsAndSolverInformation(const MBool finalTimeStep);
  void cleanUp();

  template <MInt nDim>
  void createSolver(CartesianGrid<nDim>* grid, const MInt solverId, Geometry<nDim>* geometry, MBool* propertiesGroup,
                    MBool& isActive);
  template <MInt nDim>
  void createSolver(StructuredGrid<nDim>* grid, const MInt solverId, MBool* propertiesGroups, const MPI_Comm comm);
  template <MInt nDim>
  void createSolver(const MInt solverId, const MPI_Comm comm);

  template <MInt nDim>
  void createCoupler(MInt);
  template <MInt nDim>
  PostProcessingInterface* createPostProcessing(MInt);
  ExecutionRecipe* createRecipe();

  static MBool* readPropertiesGroups();

  // variables used by the FV method
  MBool m_dualTimeStepping;
  /// Maximum number of iterations
  MInt m_maxIterations;
  /// The number of timesteps before executing a restart-backup
  MInt m_restartBackupInterval;
  /// The number of solvers
  MInt m_noSolvers = -1;
  /// The list of solvers
  std::vector<std::unique_ptr<Solver>> m_solvers;
  /// Initial adaptation
  MBool m_initialAdaptation = true;
  /// The number of couplers
  MInt m_noCouplers = 0;
  /// The list of couplers
  std::vector<std::unique_ptr<Coupling>> m_couplers;
  /// The number of postprocessing solvers
  MInt m_noPostProcessing;
  // Post-processing
  MBool m_postProcessing;
  // Post-data solverId
  MInt m_postDataSolverId = -1;
  /// Memory statistics controller
  MBool m_displayMemoryStatistics = false;
  /// post-processing in-solve position
  MBool m_ppAfterTS = false;

  // auto saving
  //  MInt m_noMinutesEndAutoSave = 0;

  /// Solver/coupler timings for performance evaluations
  MBool m_writeSolverTimings = false;
  /// Switch timings mode between ALL timings (default) and a reduced/essential timings mode
  MBool m_writeAllSolverTimings = true;
  /// Write interval for timings
  MInt m_solverTimingsWriteInterval = -1;
  /// Sampling interval for timings
  MInt m_solverTimingsSampleInterval = 1;
  MInt m_noGlobalSolverTimers = -1;
  const MInt m_maxNoSolverTimings = 100000;
  std::vector<std::vector<MFloat>> m_solverTimings{};
  std::vector<MFloat> m_solverTimingsPrevTime{};
  std::vector<MInt> m_solverTimingsTimeStep{};
  std::vector<MString> m_solverTimingsNames{};
  std::vector<std::pair<MString, MInt>> m_domainInfo{};
};

#endif
