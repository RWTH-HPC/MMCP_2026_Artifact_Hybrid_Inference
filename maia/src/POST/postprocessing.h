// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef POSTPROCESSING_H_
#define POSTPROCESSING_H_

#include <array>
#include <map>
#include <set>
#include <vector>

#include "COUPLER/couplingutils.h"
#include "GRID/cartesiangrid.h"
#include "GRID/cartesiangridproxy.h"
#include "globals.h"
#include "postdata.h"
#include "solver.h"

// Forward declarations
template <MInt nDim>
class CartesianGrid;


/** \brief compare struct for 1D coordinates
 *
 * \author A. Niemoeller
 * \date 09.12.2013
 **/
struct coord_comp_1d_ {
  MBool operator()(const MFloat& left, const MFloat& right) const {
    return ((left < right) && (std::fabs(left - right) > 1000 * MFloatEps));
  }
};

/** \brief compare struct for 2D coordinates
 *
 * \author A. Niemoeller
 * \date 09.12.2013
 **/
struct coord_comp_2d_ {
  MBool operator()(const std::pair<MFloat, MFloat>& left, const std::pair<MFloat, MFloat>& right) const {
    if(left.first < right.first) {
      return true;
    } else if(approx(left.first, right.first, 100 * MFloatEps)) {
      if(left.second < right.second) {
        return true;
      } else {
        return false;
      }
    } else {
      return false;
    }
  }
};


class PostProcessingInterface {
 public:
  PostProcessingInterface(const MInt postprocessingId_) : m_postprocessingId(postprocessingId_){};
  virtual ~PostProcessingInterface() = default;

  PostProcessingInterface(const PostProcessingInterface&) = delete;
  PostProcessingInterface& operator=(const PostProcessingInterface&) = delete;

  // virtual void postprocessPreInit() = 0;
  virtual void initPostProcessing() = 0;
  virtual void postprocessPreSolve() = 0;
  virtual void postprocessPostSolve() = 0;
  virtual void postprocessInSolve(MBool finalTimeStep) = 0;
  virtual void postprocessSolution() = 0;
  // to access timer functions in PostProcessingController
  virtual Solver* mSolver() const = 0;

  MInt a_postprocessingId() const { return m_postprocessingId; };
  MBool m_finalTimeStep = false;

 protected:
  MInt m_postprocessingId;
};

template <MInt nDim, class ppType>
class PostProcessing : virtual public PostProcessingInterface {
 public:
  using PPGrid = CartesianGrid<nDim>;
  using Cell = typename PPGrid::Cell;
  using GridProxy = typename maia::grid::Proxy<nDim>;
  using Data = PostData<nDim>;

  // Constructor
  PostProcessing(MInt postprocessingId_, PostData<nDim>* data) : PostProcessingInterface(postprocessingId_) {
    m_postData = data;
  }

  ~PostProcessing();

  // void postprocessPreInit() override;
  void initPostProcessing() override;
  void postprocessPreSolve() override;
  void postprocessPostSolve() override;
  void postprocessInSolve(MBool finalTimeStep) override;
  void postprocessSolution() override;
  // to access timer functions in PostProcessingController
  Solver* mSolver() const override { return static_cast<Solver*>(&pp->solver()); };

  Data& postData() const { return *m_postData; }

 protected:
  MBool isActive() const { return pp->solver().grid().isActive(); };

  void transferSensorWeights();

  void initReduceToLevelAvg();
  void initTimeStepProperties();
  void initTimeStepPropertiesSlice();
  void initTimeStepPropertiesLine();
  virtual void initAveragingProperties();
  void initProbePoint();
  void initProbeLine();
  void initProbeSlice();
  void initSpatialAveraging();

  void initProbeArbitraryLine();
  void initProbeArbitrarySlice();

  void initProbeLinePeriodic();
  virtual void probeLinePeriodic() { TERMM(1, "Not implemented for this solverType"); };
  virtual void probeLinePeriodicPost();

  virtual void initSprayData();
  virtual void computeSprayData() { TERMM(1, "Not implemented for this solverType"); };
  virtual void writeSprayData() { TERMM(1, "Not implemented for this solverType"); };

  // implemented in PostProcessingLpt
  virtual void initLPTSolutionFile() { TERMM(1, "Not implemented for this solverType"); };
  virtual void writeLPTSolutionFile() { TERMM(1, "Not implemented for this solverType"); };
  virtual void initParticleStatistics() { TERMM(1, "Not implemented for this solverType"); };
  virtual void computeParticleStatistics() { TERMM(1, "Not implemented for this solverType"); };

  // implemented in PostProcessingLbLpt
  virtual void initPLIsoTurbulenceStatistics() { TERMM(1, "Not implemented for this solverType"); };
  virtual void computePLIsoTurbulenceStatistics() { TERMM(1, "Not implemented for this solverType"); };

  // allocates memory for average pre/in/post
  virtual void initAverageVariables();

  void pp_saveCoarseSolution();
  void averageSolutions();
  virtual void averageSolutionsInSolve();
  template <MBool inSolution = false>
  void computeAndSaveDivergence();
  void computeAndSaveMean();
  void initAverageSolutionsSlice();
  void averageSolutionsSlice();
  void probeLocation();
  void findClosestProbePointsGridCells();

  // implemented in PostProcessingLb
  virtual void initPointSamplingData() { TERMM(1, "Not implemented for this solverType"); };
  virtual void savePointSamplingData() { TERMM(1, "Not implemented for this solverType"); };
  virtual void initSurfaceSamplingData() { TERMM(1, "Not implemented for this solverType"); };
  virtual void saveSurfaceSamplingData() { TERMM(1, "Not implemented for this solverType"); };
  virtual void initVolumeSamplingData() { TERMM(1, "Not implemented for this solverType"); };
  virtual void saveVolumeSamplingData() { TERMM(1, "Not implemented for this solverType"); };
  virtual void initIsoTurbulenceStatistics() { TERMM(1, "Not implemented for this solverType"); };
  virtual void computeIsoTurbulenceStatistics() { TERMM(1, "Not implemented for this solverType"); };

  void initWritePointData();
  void writePointData();

  void spatialAveraging();
  void spatialAveragingPost();
  void createCellToMap1D(std::map<MFloat, MInt, coord_comp_1d_>& coordinates,
                         std::map<MFloat, MFloat, coord_comp_1d_>& cell_to_map);
  void createCellToMap2D(std::map<std::pair<MFloat, MFloat>, MInt, coord_comp_2d_>& coordinates,
                         std::map<std::pair<MFloat, MFloat>, std::pair<MFloat, MFloat>, coord_comp_2d_>& cell_to_map);

  void probeLine();
  void probeLinePre();
  void probeLinePost();
  void probeSlice();
  void probeSliceIn();
  void probeSlicePre();
  void probeSlicePost();

  void collectMinLvlCells();
  void findContainingCell(const MFloat* coord, MInt& cellId);
  void probeArbitraryLine();
  void probeArbitraryLinePost();
  void probeArbitrarySlice();
  void probeArbitrarySlicePost();

  virtual void initMovingAverage();
  void movingAverage();
  void movingAveragePost();

  void initCorrelation();
  void initPeriodicSliceAverage();
  void periodicSliceAverage();

  MInt findNearestGridCell(const MFloat* coord);
  void movePointsToGrid(MFloat* in_points, MInt in_noPoints, MInt in_moveType);

  void getSampleVariables(MInt cellId, const MFloat*& vars, MBool mode);
  void getSampleVariables(MInt const cellId, std::vector<MFloat>& vars);
  void calcSamplingVar(const MInt cellId, const MInt sampleVarId, MFloat* const vars);

  void saveSliceAiaFileFormat(const MInt step, const MInt noVars, MFloatScratchSpace& vars, const MInt sliceId);

  // TODO labels:POST This is a very inefficient way of doing. Here, the derivative of ALL vars in ALL directions are
  // calculated. Only trace of velocity gradient is needed !
  MFloat calcDivergence(const MInt cellIdSolver) {
    const MInt noVariables = pp->solver().noVariables();
    std::vector<MFloat> cellVarsDeriv(noVariables * nDim);
    MFloat divergence = 0.0;
    if(pp->getSampleVarsDerivatives(cellIdSolver, cellVarsDeriv)) {
      const MFloatTensor deriv(cellVarsDeriv.data(), noVariables, nDim);
      for(MInt d = 0; d < nDim; d++) {
        divergence += deriv(d, d);
      }
    }
    return divergence;
  }

  virtual void calcVorticity(const MFloatTensor& deriv, MFloat vorticity[nDim * 2 - 3]) {
    if constexpr(nDim == 2) {
      vorticity[0] = deriv(1, 0) - deriv(0, 1);
    } else { //(nDim ==3)
      vorticity[0] = deriv(2, 1) - deriv(1, 2);
      vorticity[1] = deriv(0, 2) - deriv(2, 0);
      vorticity[2] = deriv(1, 0) - deriv(0, 1);
    }
  }
  // not implemented for LB
  virtual void getVorticity(MFloat* const vorticity) {
    std::ignore = vorticity;
    TERMM(1, "Not implemented for this solverType");
  };
  virtual void getVorticityT(MFloat* const vorticity) {
    std::ignore = vorticity;
    TERMM(1, "Not implemented for this solverType");
  };
  virtual void getPrimitiveVariables(MInt, MFloat*, MFloat*, MInt) { TERMM(1, "Not implemented for this solverType"); };

  virtual void
  computeAcousticSourceTermQe(MFloatScratchSpace&, MFloatScratchSpace&, MFloatScratchSpace&, MFloatScratchSpace&) {
    TERMM(1, "Not implemented for this solverType");
  }

  void getSampleVarsDerivatives(const MInt cellId, const MFloat*& vars) {
    std::ignore = cellId;
    vars = nullptr;
    // TERMM(1, "Not implemented for this solverType");
  };
  MBool getSampleVarsDerivatives(const MInt /*cellId*/, std::vector<MFloat>& /*vars*/) {
    return false; // = not implemented for this solver
    // TERMM(1, "Not implemented for this solverType");
  };
  virtual MFloat& vorticityAtCell(const MInt cellId, const MInt dir) {
    std::ignore = cellId;
    std::ignore = dir;
    TERMM(1, "Not implemented for this solverType");
  };
  virtual MFloat getBoundaryHeatFlux(const MInt cellId) const {
    std::ignore = cellId;
    TERMM(1, "Not implemented for this solver");
  };
  virtual void getPrimitiveVarNames(MString* names) const {
    std::ignore = names;
    TERMM(1, "Not implemented for this solverType");
  };

  MInt m_noPostprocessingOps{};
  MString* m_postprocessingOps = nullptr;

  // this vector holds all the functions to be called
  // it is initialized with three entries (each position of postprocessing)
  typedef void (PostProcessing::*tpost)();
  typedef std::vector<tpost> tvecpost;
  tvecpost m_postprocessingSolution;
  std::vector<tvecpost> m_postprocessingMethods;
  std::vector<std::vector<MString>> m_postprocessingMethodsDesc;

  // Relevant for CAA simulations
  void neededMeanVarsForSourceTerm(const MInt sourceTerm, std::vector<MInt>& meanVars) const;

 public:
  MInt m_noLocalVars{};
  MInt m_restartTimeStep{};
  MBool m_statisticCombustionAnalysis{};
  MBool m_acousticAnalysis{};

 private:
  // Status of postprocessPreInit()
  /* MBool m_isPreInitialized = false; */

  MString m_solverType;
  ppType* pp = static_cast<ppType*>(this); ///< CRTP

  // PP_ToDo-GridProxy:change this to proxy of postData and check that it is still working
  //                  -> create new testcases
  GridProxy* m_gridProxy;

 protected:
  // Pointer to PostData class
  Data* m_postData = nullptr;

  MBool isMeanFile() {
    if(m_postData == nullptr ||
       // necessary for postprocessing without postdata
       postData().m_solverId == pp->solver().m_solverId) {
      return false;
    } else {
      return postData().isMeanFile();
    }
  }

  // Number of variables of the corresponding solver
  MInt m_noVariables{};

  // Averaging variables
  // MFloat** m_summedVars = nullptr;
  // MFloat** m_square = nullptr;
  // MFloat** m_cube = nullptr;
  // MFloat** m_fourth = nullptr;
  MBool m_square = false;
  MBool m_skewness = false;
  MBool m_kurtosis = false;
  MBool m_computeAndSaveMean = false;
  // not used in any testcase
  MBool m_twoPass = false;

  // Kahan summation
  // not used in any testcase
  MBool m_useKahan{};
  // MFloat** m_cSum = nullptr;
  // MFloat** m_ySum = nullptr;
  // MFloat** m_tSum = nullptr;
  // MFloat** m_cSquare = nullptr;
  // MFloat** m_ySquare = nullptr;
  // MFloat** m_tSquare = nullptr;
  // MFloat** m_cCube = nullptr;
  // MFloat** m_yCube = nullptr;
  // MFloat** m_tCube = nullptr;
  // MFloat** m_cFourth = nullptr;
  // MFloat** m_yFourth = nullptr;
  // MFloat** m_tFourth = nullptr;

  MBool m_correlation = false;

  // spanwise average
  std::vector<MInt>* m_cell2globalIndex = nullptr;
  MInt* m_noPeriodicSliceCells = nullptr;
  MFloat** m_sliceAverage = nullptr;
  MInt m_globalnoSlicePositions = 0;
  // MFloat* m_globalPosX = nullptr;
  // MFloat* m_globalPosY = nullptr;
  // MInt* m_globalCountZ = nullptr;
  std::multimap<MFloat, MFloat> m_sliceGlobalPositions;

  MBool m_averageVorticity = false;
  MBool m_averageSpeedOfSound = false;
  MFloat m_gamma = -1.0;

  // moving average
  MInt m_movingAverageStartTimestep = 0;
  MInt m_movingAverageStopTimestep = 0;
  MInt m_movingAverageInterval = 1;
  MInt m_movingAverageDataPoints{};
  MInt m_movingAverageCounter{};
  MFloat** m_movAvgVariables = nullptr;
  MInt m_movAvgNoVariables{};
  MString* m_movAvgVarNames = nullptr;

  // spatial averaging
  MInt m_spatialDirection1{};
  MInt m_spatialDirection2{};
  MFloat m_spatialWeightSinglePoint{};
  MFloat m_spatialLvlWeight[20]{};
  MInt* m_spatialDispls = nullptr;
  MInt* m_spatialVarsDispls = nullptr;
  MInt* m_spatialRecvcnts = nullptr;
  MInt* m_spatialVarsRecvcnts = nullptr;
  MInt m_spatialCoordSum{};

  std::map<MFloat, MInt, coord_comp_1d_> m_spatialLineCoordinates;
  std::map<MFloat, MFloat, coord_comp_1d_> m_spatialLineCellToMap;

  std::map<MFloat, MInt, coord_comp_1d_> m_spatialGlobalLineCoordinates;
  std::map<MFloat, MFloat, coord_comp_1d_> m_spatialGlobalLineCellToMap;

  MInt m_spatialLineNoCells = -1;
  MFloat* m_spatialLineAllVars = nullptr;
  MFloat* m_spatialLineAllCoord = nullptr;

  std::map<std::pair<MFloat, MFloat>, MInt, coord_comp_2d_> m_spatialPlaneCoordinates;
  std::map<std::pair<MFloat, MFloat>, std::pair<MFloat, MFloat>, coord_comp_2d_> m_spatialPlaneCellToMap;

  std::map<std::pair<MFloat, MFloat>, MInt, coord_comp_2d_> m_spatialGlobalPlaneCoordinates;
  std::map<std::pair<MFloat, MFloat>, std::pair<MFloat, MFloat>, coord_comp_2d_> m_spatialGlobalPlaneCellToMap;

  MInt m_spatialPlaneNoCells = -1;
  MFloat* m_spatialPlaneAllVars = nullptr;
  MFloat* m_spatialPlaneAllCoord = nullptr;

  // line probing
  MInt m_noProbeLines = -1;
  MInt* m_probeLineDirection = nullptr;
  MInt* m_probeLinePeriodic = nullptr;
  MFloat** m_probeLineCoordinates = nullptr;
  MInt** m_probeLineIds = nullptr;
  MInt* m_noProbeLineIds = nullptr;
  MFloat** m_probeLinePositions = nullptr;
  MInt* m_globalNoProbeLineIds = nullptr;
  MInt* m_probeLineOffsets = nullptr;
  MInt** m_noGlobalProbeLineIds = nullptr;
  // line time step properties
  MInt m_probeLineInterval = -1;
  MInt m_probeLineStartTimestep = 0;
  MInt m_probeLineStopTimestep = -1;

  // line probe correlation
  MInt* m_correlationDirection = nullptr;
  MInt* m_correlationVariableIndex = nullptr;
  MFloat** m_correlationCoordinates = nullptr;
  MInt m_noCorrelationLines = -1;
  MInt* m_noCorrelationIds = nullptr;
  MInt** m_correlationIds = nullptr;
  MInt** m_globalCorrelationIds = nullptr;
  MFloat** m_globalCorrelationPositions = nullptr;
  MInt* m_globalNoCorrelationIds = nullptr;
  MFloat** m_correlationPositions = nullptr;
  MInt** m_correlationIndexMapping = nullptr;
  MFloat** m_correlationExchangeVar = nullptr;
  MFloat** m_correlationExchangeVarMean = nullptr;
  MFloat** m_correlationExchangeVarRMS = nullptr;
  MFloat** m_globalCorrelationExchangeVar = nullptr;
  MFloat** m_globalCorrelationExchangeVarMean = nullptr;
  MFloat** m_globalCorrelationExchangeVarRMS = nullptr;

  // line probe averaging
  MInt* m_probeLineAverageDirection = nullptr;
  MFloat** m_probeLineAverageCoordinates = nullptr;
  MInt** m_probeLineAverageCoordinatesSign = nullptr;
  MInt** m_noProbeLineAverageIds = nullptr;
  MFloat*** m_globalProbeLineAverageVars = nullptr;
  // MFloat*** m_globalProbeLineAverageRMS = nullptr;
  MInt m_noProbeLineAverageSteps;
  MInt** m_globalNoProbeLineAverageIds = nullptr;
  MInt*** m_probeLineAverageIds = nullptr;
  MFloat*** m_probeLineAveragePositions = nullptr;

  // slice probing
  MInt m_noProbeSlices = -1;
  MInt* m_probeSliceDir = nullptr;
  MFloat* m_probeSliceCoordinate = nullptr;
  MInt** m_probeSliceIds = nullptr;
  MInt* m_noProbeSliceIds = nullptr;
  MFloat** m_probeSlicePositions = nullptr;
  MInt* m_globalNoProbeSliceIds = nullptr;
  MInt* m_probeSliceOffsets = nullptr;
  MInt** m_noGlobalProbeSliceIds = nullptr;
  MString* m_probeSliceGridNames = nullptr;
  MString* m_sliceAxis = nullptr;
  MFloat* m_sliceIntercept = nullptr;
  MBool m_sliceAiaFileFormat = false;
  MBool m_optimizedSliceIo = true;
  MInt* m_noProbeSliceNoHilbertIds = nullptr;
  MInt* m_noProbeSliceMaxNoHilbertIds = nullptr;
  MInt** m_noProbeSliceHilbertInfo = nullptr;
  MInt* m_noProbeSliceNoContHilbertIds = nullptr;
  MInt* m_noProbeSliceMaxNoContHilbertIds = nullptr;
  MInt** m_noProbeSliceContHilbertInfo = nullptr;
  // slice time step properties
  MInt m_probeSliceInterval = -1;
  MInt m_probeSliceStartTimestep = 0;
  MInt m_probeSliceStopTimestep = -1;

  /// List of slice variables
  std::vector<MInt> m_sliceVarIds{};
  /// Number of variables for each slice variable
  std::vector<MInt> m_noSliceVars{};
  /// List of variable names for each slice variable
  std::vector<std::vector<MString>> m_sliceVarNames{};

  MInt m_noMinLvlIds{};
  MInt* m_minLvlIds = nullptr;

  // arbitrary probing
  MInt m_movePointsToGrid{};
  MBool m_spatialAveraging{};

  // arbitrary line probing
  MFloat* m_arbLinePoints = nullptr;
  MInt m_noArbLines = 0;
  MInt* m_noArbLineIds = nullptr;
  MFloat* m_arbLinePointsDistribution = nullptr;
  MInt* m_arbLineHasNewDistribution = nullptr;
  MInt* m_arbLineFull = nullptr;
  MInt* m_moveLinePointsToGrid = nullptr;
  MInt* m_globalNoArbLineIds = nullptr;
  MInt** m_arbLineIds = nullptr;
  MInt* m_arbLineOffsets = nullptr;
  MFloat** m_arbLineCoordinates = nullptr;

  // arbitrary slice probing
  MFloat* m_arbSlicePoints = nullptr;
  MInt m_noArbSlices{};
  MInt* m_noArbSlicePoints = nullptr;
  MFloat* m_arbSlicePointsDistribution = nullptr;
  MInt* m_globalNoArbSlicePoints = nullptr;
  MInt** m_arbSliceIds = nullptr;
  MInt* m_arbSliceOffsets = nullptr;
  MFloat** m_arbSliceCoordinates = nullptr;

  // loading file for averaging
  MString m_postprocessFileName = "";

  // genereal purpose
  // PP_ToDo-LB:delete this
  MFloat** m_localVars = nullptr;

  // averaging
  MBool m_averageRestart{};
  MInt m_averageInterval{};
  MInt m_averageStartTimestep{};
  MInt m_averageStopTimestep{};

  MInt m_noAveragedVorticities{};
  MInt m_noSpeedOfSoundVars{};
  MInt m_noSourceVars{};
  MBool m_needVorticity = false;

  // Reynolds stresses
  MString m_ReStressesAverageFileName;

  // probing
  MFloat** m_probeCoordinates = nullptr;
  MString m_probePath;
  MInt* m_probeCellIds = nullptr;
  std::ofstream* m_probeFileStreams = nullptr;
  MInt m_noProbePoints{};
  MInt m_probeInterval;
  MInt m_probeWriteInterval;

  // Point data
  MString m_pdFileName;
  MInt m_pdStartTimestep{};
  MInt m_pdStopTimestep{};
  MInt m_pdRestartInterval{};
  std::vector<MFloat> m_pdPoints;
  std::vector<MInt> m_pdCells;
  std::vector<MFloat> m_pdVars;
  MInt m_pdNoPoints{};

  // spray data
  MInt m_sprayComputeInterval = 50;
  MInt m_sprayWriteInterval = 50;
  MInt m_sprayDataSize = 50;
  MInt m_sprayDataStep = 0;

  MFloat** m_particleCV = nullptr;
  MFloat** m_particlePen = nullptr;
  MFloat** m_sprayStat = nullptr;
  MFloat** m_sprayTimes = nullptr;
  MFloat** m_injectionData = nullptr;

  MFloat** m_vapourCV = nullptr;
  MFloat** m_vapourPen = nullptr;

  // LPT solution file
  MInt m_LPTSolutionInterval = 50;
  MBool m_forceOutput = false;

  //////////////////////////////////////////////////////////////////////////////

  // Relevant for CAA simulations
  // List of active source terms
  std::vector<MInt> m_activeSourceTerms{};
  // Hold indices for source terms
  struct ST {
    static constexpr const MInt Q_mI = 0;
    static constexpr const MInt Q_mI_linear = 1;
    static constexpr const MInt Q_mII = 2;
    static constexpr const MInt Q_mIII = 3;
    static constexpr const MInt Q_mIV = 4;
    static constexpr const MInt Q_c = 5;
    static constexpr const MInt Q_e = 6;

    static constexpr const MInt totalNoSourceTerms = 7;
  };

  // Source term names corresponding to ST
  std::array<MString, ST::totalNoSourceTerms> s_sourceTermNames = {
      {"q_mI", "q_mI_linear", "q_mII", "q_mIII", "q_mIV", "q_c", "q_e"}};

  // List of active mean variables for all active source terms
  std::set<MInt> m_activeMeanVars{};
  // Hold indices for mean variables
  struct MV {
    // Mean Lamb vector
    static constexpr const MInt LAMB0 = 0;

    // Mean vorticity
    static constexpr const MInt VORT0 = 1;

    // du/dx, du/dy, dw/dz for the divergence
    static constexpr const MInt DU = 2;

    // Mean gradient of rho
    static constexpr const MInt DRHO = 3;

    // Mean gradient of p
    static constexpr const MInt DP = 4;

    // Mean gradient of rho*div(u)
    static constexpr const MInt RHODIVU = 5;

    // Mean gradient of u*grad(rho)
    static constexpr const MInt UGRADRHO = 6;

    // Mean of (gradient of p divided by rho)
    static constexpr const MInt GRADPRHO = 7;

    // Mean gradients of velocity components (contains MV::DU)
    static constexpr const MInt GRADU = 8;

    // Sum of products of velocity and velocity gradients:
    // u * grad(u) + v * grad(v) + w * grad(w)
    static constexpr const MInt UGRADU = 9;
  };
};


#endif // POSTPROCESSING_H_
