// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef STRUCTUREDPOSTPROCESSING_H
#define STRUCTUREDPOSTPROCESSING_H

#include <map>
#include <vector>
#include "INCLUDE/maiatypes.h"

#include "INCLUDE/maiaconstants.h"
#include "UTIL/functions.h"

// Forard declarations
class StructuredCell;

template <MInt nDim, class SolverType>
class StructuredPostprocessing {
  typedef StructuredCell PPCell;

 public:
  StructuredPostprocessing();
  ~StructuredPostprocessing();

  void postprocessPreInit();
  void postprocessPreSolve();
  void postprocessPostSolve();
  void postprocessInSolve();

 protected:
  void initStructuredPostprocessing();
  void initAverageIn();
  void initAverageVariables();
  void initTimeStepProperties();
  void initMovingAverage();
  void initProductionVariables();
  void initDissipationVariables();
  void averageSolutionsInSolve();
  void averageSolutions();
  void addAveragingSample();
  void addTempWaveSample();
  void saveAveragedSolution(MInt);
  void computeAveragedSolution();
  void computeAverageSkinFriction();
  void subtractPeriodicFluctuations();
  void subtractMean();
  void movingAverage();
  void movingAveragePost();
  void computeProductionTerms();
  void computeDissipationTerms();
  void decomposeCf();
  void decomposeCfDouble();
  void writeGradients();
  void loadAveragedSolution();

  void saveAverageRestart();

  void loadMeanFile(const MChar* fileName);
  void getSampleVariables(MInt cellId, const MFloat*& vars);
  MInt getNoPPVars();
  MInt getNoVars();
  MInt getNoPPSquareVars();

  MBool m_postprocessing;
  MInt m_noPostprocessingOps;
  MString* m_postprocessingOps = nullptr;

  MInt m_dissFileStart;
  MInt m_dissFileEnd;
  MInt m_dissFileStep;
  MString m_dissFileDir;
  MString m_dissFilePrefix;
  MInt m_dissFileBoxNr;

  // this vector holds all the functions to be called
  // it is initialized with three entries (each position of postprocessing)
  typedef void (StructuredPostprocessing::*tpost)();
  typedef std::vector<tpost> tvecpost;
  std::vector<tvecpost> m_postprocessingMethods;
  std::vector<std::vector<MString>> m_postprocessingMethodsDesc;

 public:
  MInt m_restartTimeStep = -1;

 protected:
  // Number of variables of the corresponding solver
  MInt m_noVariables;

  // Averaging variables
  MFloat** m_summedVars = nullptr;
  MFloat** m_square = nullptr;
  MFloat** m_cube = nullptr;
  MFloat** m_fourth = nullptr;
  MFloat** m_tempWaveSample = nullptr;
  MFloat** m_favre = nullptr;

  // Kahan summation
  MBool m_useKahan;
  MFloat** m_cSum = nullptr;
  MFloat** m_ySum = nullptr;
  MFloat** m_tSum = nullptr;
  MFloat** m_cSquare = nullptr;
  MFloat** m_ySquare = nullptr;
  MFloat** m_tSquare = nullptr;
  MFloat** m_cCube = nullptr;
  MFloat** m_yCube = nullptr;
  MFloat** m_tCube = nullptr;
  MFloat** m_cFourth = nullptr;
  MFloat** m_yFourth = nullptr;
  MFloat** m_tFourth = nullptr;

  MBool m_twoPass;
  MBool m_skewness;
  MBool m_kurtosis;
  MBool m_averageVorticity = false;
  MBool m_averagingFavre = false;
  MFloat** m_production = nullptr;
  MString* m_avgVariableNames = nullptr;
  MString* m_avgFavreNames = nullptr;
  MFloat* m_dissipation;
  MFloat** m_gradients;
  MString* m_gradientNames;

  MInt m_movingAvgInterval;
  MInt m_movingAvgDataPoints;
  MInt m_movingAvgCounter;
  MInt m_movAvgNoVariables;
  MFloat** m_movAvgVariables = nullptr;
  MString* m_movAvgVarNames = nullptr;

  // periodic flucs
  MFloat** m_spanAvg;

  // averaging
  MInt m_averageInterval = 0;
  MInt m_averageStartTimestep = 0;
  MInt m_averageStopTimestep = 0;
  MInt m_averageRestartInterval = 0;
  MInt m_averageRestart = 0;
  MInt m_noSamples = 0;

  // moving Grid
  MBool m_movingGrid;

  static const MInt xsd = 0;
  static const MInt ysd = 1;
  static const MInt zsd = 2;

  MBool m_computeProductionTerms;
  MBool m_computeDissipationTerms;

  MString m_postprocessFileName;

  MFloat m_sutherlandConstant;
  MFloat m_sutherlandPlusOne;

 private:
  MString m_solverType;
  SolverType* ppsolver() { return static_cast<SolverType*>(this); } ///< CRTP
};
#endif
