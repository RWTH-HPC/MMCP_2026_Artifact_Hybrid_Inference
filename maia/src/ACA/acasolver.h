// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef ACASOLVER_H_
#define ACASOLVER_H_

#include "COMM/mpioverride.h"
#include "INCLUDE/maiamacro.h"
#include "INCLUDE/maiatypes.h"
#include "IO/parallelio.h"
#include "UTIL/debug.h"
#include "UTIL/functions.h"
#include "acaobserverdatacollector.h"
#include "acapostprocessing.h"
#include "acasurfacedatacollector.h"
#include "compiler_config.h"
#include "solver.h"


namespace maia {
namespace acoustic {

// Create struct for easy timer identification
struct Timers_ {
  // Enum to store timer "names"
  enum {
    SolverType,
    Run,
    InitParallelization,
    ReadSurfaceData,
    InitObservers,
    CalcSourceTerms,
    WindowAndTransformSources,
    CalcSurfaceIntegrals,
    CombineSurfaceIntegrals,
    BacktransformObservers,
    SaveObserverSignals,
    LoadObserverSignals,
    Postprocessing,

    // Special enum value used to initialize timer array
    _count
  };
};

/// Parameter for analytical acoustic source terms
struct SourceParameters {
  MBool perturbed;
  MFloat omega;
  MFloat rho0;
  MFloat beta;
  MFloat amplitude;
  MFloat Ma;
};

struct SourceVars {
  MFloat* rho = nullptr;
  MFloat* u = nullptr;
  MFloat* v = nullptr;
  MFloat* w = nullptr;
  MFloat* p = nullptr;
};

} // namespace acoustic
} // namespace maia


/// \brief Acoustic extrapolation of near/near-far field data to the true far-field using an
///        integral method (currently only FWH, but should allow for other methods to be implemented
///        as well)
///
/// Reference for implementation/validation of FWH solver: 'Predition of far-field noise using the
/// Ffowcs Williams-Hawkings method', Maximilian Kerstan, Master thesis, 2020.
template <MInt nDim>
class AcaSolver : public Solver {
  // Types
 private:
  using SurfaceDataCollector = maia::acoustic_analogy::collector::SurfaceDataCollector<nDim>;
  using ObserverDataCollector = maia::acoustic_analogy::observer_collector::ObserverDataCollector<nDim>;
  using SourceParameters = maia::acoustic::SourceParameters;
  using SourceVars = maia::acoustic::SourceVars;
  using Timers = maia::acoustic::Timers_;

  enum NONDIMBASIS {
    NONDIMACA = 0,
    NONDIMSTAGNATION = 1,
    NONDIMLB = 2,
  };

  // Methods
 public:
  AcaSolver(const MInt solverId, const MPI_Comm comm);
  ~AcaSolver() override {
    averageTimer();
    RECORD_TIMER_STOP(m_timers[Timers::SolverType]);
  };

  MBool solutionStep() override;
  MBool solverConverged() { return m_isFinished; };
  MBool isActive() const override { return true; } // All ranks are active

  void initSolver() override;
  void finalizeInitSolver() override{};

  // Required to override but should not be used
  MFloat time() const override {
    m_log << "WARNING: Empty body of AcaSolver::time() called" << std::endl;
    return NAN;
  }
  MInt noVariables() const override { return 0; }
  MInt noInternalCells() const override { return 0; }
  void preTimeStep() override {}
  void postTimeStep() {}
  virtual void saveSolverSolution(const MBool, const MBool){};
  virtual void cleanUp(){};

  // Solver internal methods
 private:
  // Constructor methods
  void initTimers();
  void averageTimer();
  void setInputOutputProperties();
  void setNumericalProperties();

  /// Return number of surfaces
  MInt noSurfaces() { return m_noSurfaces; }

  /// Return if this rank has data of this surface
  MBool hasSurface(const MInt surfaceId) { return m_noSurfaceElements.at(surfaceId) > 0; }

  /// Return number of (local) surfaces elements for a surface
  MInt noSurfaceElements(const MInt surfaceId) { return m_noSurfaceElements.at(surfaceId); }

  /// Return total number of (local) surfaces elements
  MInt totalNoSurfaceElements() { return std::accumulate(m_noSurfaceElements.begin(), m_noSurfaceElements.end(), 0); }

  /// Return local offset of a surface in m_surfaceData
  MInt localSurfaceOffset(const MInt sId) {
    return std::accumulate(&m_noSurfaceElements.at(0), &m_noSurfaceElements.at(sId), 0);
  }

  /// Return number of samples (i.e. number of time steps)
  MInt noSamples() { return m_noSamples; }

  /// Return local number of observer points
  MInt noObservers() { return m_noObservers; }
  /// Return local number of observer points
  MInt observerOffset() { return m_offsetObserver; }
  /// Return number of observer points
  MInt globalNoObservers() { return m_noGlobalObservers; }
  /// Return local observer id
  MInt a_localObserverId(const MInt globalObserverId) {
    if(globalObserverId < m_offsetObserver || (m_offsetObserver + m_noObservers) < globalObserverId) {
      return -1;
    } else {
      return globalObserverId - m_offsetObserver;
    }
  };
  /// Return global observer id
  MInt a_globalObserverId(const MInt localObserverId) { return localObserverId + m_offsetObserver; };
  /// Accessor for global observer coordinates
  MFloat a_globalObserverCoordinate(const MInt globalObserverId, const MInt dir) {
    return m_globalObserverCoordinates[globalObserverId * nDim + dir];
  };
  /// Accessor for global observer coordinates
  std::array<MFloat, nDim> getGlobalObserverCoordinates(const MInt globalObserverId) {
    std::array<MFloat, nDim> coord;
    for(MInt dir = 0; dir < nDim; dir++) {
      coord[dir] = a_globalObserverCoordinate(globalObserverId, dir);
    }
    return coord;
  };
  /// Get domain id of (global) observer
  MInt getObserverDomainId(const MInt globalObserverId) {
    const MInt minNoObsPerDomain = globalNoObservers() / noDomains();
    const MInt remainderObs = globalNoObservers() % noDomains();
    if(globalObserverId < remainderObs * (minNoObsPerDomain + 1)) {
      return globalObserverId / (minNoObsPerDomain + 1);
    } else {
      return (globalObserverId - remainderObs) / minNoObsPerDomain;
    }
  };

  /// Initialize observer points
  void initObserverPoints();

  /// Initialize parallelization (distribution of surface segments on all domains)
  void initParallelization();
  /// Distribute observers on all ranks
  void distributeObservers();

  /// IO methods
  MString getSurfaceDataFileName(const MInt surfaceId, const MInt fileNo);
  void readSurfaceData();
  void readSurfaceDataFromFile(const MInt surfaceId);
  void storeSurfaceData(const MInt surfaceId);
  void readNoSegmentsAndSamples(const MInt surfaceId, MInt& noSegments, MInt& noSamples, MString& inputFileName);

  /// Info functions
  inline void printMessage(const MString msg) {
    m_log << msg << std::endl;
    if(domainId() == 0) std::cout << msg << std::endl;
  };

  /// Main compute functions
  void calcSourceTerms();
  void calcWindowFactors(MFloat* const p_window);
  void windowAndTransformSources();
  void windowAndTransformObservers();
  void calcFrequencies();
  void calcSurfaceIntegralsForObserver(const std::array<MFloat, nDim> coord, MFloat* const p_complexVariables);
  void calcSurfaceIntegralsAtSourceTime(const MInt globalObserverId);
  MFloat calcTimeDerivative(const MInt segmentId, const MInt varId, const MInt tau);
  void interpolateFromSourceToObserverTime(const MFloat distMinObservers, MFloat* const p_perturbedFWH);
  void combineSurfaceIntegrals(const MInt globalObserverId, MFloat* const p_complexVariables);
  void combineSurfaceIntegralsInTime(const MInt globalObserverId, MFloat* const p_perturbedFWH);
  void calcMinDistanceForObservers(MFloat* const distMinObservers);
  void calculateObserverInFrequency();
  void calculateObserverInTime();
  void backtransformObservers();
  void saveObserverSignals();
  void loadObserverSignals();
  void postprocessObserverSignals();

  /// Functions for Fourier Transformation
  void FastFourierTransform(MFloat* dataArray, const MInt size, const MInt dir, MFloat* result, const MBool inPlace);
  void DiscreteFourierTransform(MFloat* dataArray, const MInt size, const MInt dir, MFloat* result,
                                const MBool inPlace);
  void checkNoSamplesPotencyOfTwo();

  /// Source term compute kernels
  void calcSourceTermsFwhFreq(const MInt segmentId);
  void calcSourceTermsFwhTime(const MInt segmentId);

  /// Coordinate transformation
  void transformCoordinatesToFlowDirection2D(MFloat* const fX, MFloat* const fY);
  void transformCoordinatesToFlowDirection3D(MFloat* const fX, MFloat* const fY, MFloat* const fZ);
  void transformBackToOriginalCoordinate2D(MFloat* const fX, MFloat* const fY);
  void transformBackToOriginalCoordinate3D(MFloat* const fX, MFloat* const fY, MFloat* const fZ);

  /// Generation of surfaces
  void generateSurfaces();
  void generateSurfacePlane(const MInt sId);
  void generateSurfaceCircle(const MInt sId);
  MInt readNoSegmentsSurfaceGeneration(const MInt sId, std::vector<MInt>& noSegments, const MInt lengthCheck = -1);

  /// Generation of analytical input data
  void generateSurfaceData();
  /// Monopole
  void genMonopoleAnalytic2D(const MFloat* coord, const SourceParameters& param, SourceVars& vars);
  void genMonopoleAnalytic3D(const MFloat* coord, const SourceParameters& param, SourceVars& vars);
  /// Dipole
  void genDipoleAnalytic2D(const MFloat* coord, const SourceParameters& param, SourceVars& vars);
  void genDipoleAnalytic3D(const MFloat* coord, const SourceParameters& param, SourceVars& vars);
  /// Quadrupole
  void genQuadrupoleAnalytic2D(const MFloat* coord, const SourceParameters& param, SourceVars& vars);
  void genQuadrupoleAnalytic3D(const MFloat* coord, const SourceParameters& param, SourceVars& vars);
  /// Inviscid vortex convection
  void genVortexConvectionAnalytic2D(const MFloat* obsCoord, const SourceParameters& m_sourceParameters,
                                     SourceVars& sourceVars);

  /// Symmetry boundary condition
  void applySymmetryBc(const std::array<MFloat, nDim> coord, MFloat* const p_complexVariables);

  // void genSurfaceDataNumeric2D(const MInt segmentId, const SourceParameters& param);
  // void genSurfaceDataNumeric3D(const MInt segmentId, const SourceParameters& &param);

  // Postprocessing operations
  std::vector<std::unique_ptr<AcaPostProcessing>> m_post;
  void calcObsPressureAnalytic();

  // Change the non-dimensionalizion of input data from a previous simulation
  void initConversionFactors();
  void changeDimensionsSurfaceData();
  void changeDimensionsObserverData();
  void computeDt();

  // Accuracy computation
  void computeAccuracy();

  /// Member variables
  /// Solver status
  MBool m_isFinished = false;

  MInt m_acaResultsNondimMode = 0;

  /// Conversion-Factors for the input conversion
  struct ConversionFactors {
    MFloat length = 1.0;
    MFloat density = 1.0;
    MFloat velocity = 1.0;
    MFloat time = 1.0;
    MFloat pressure = 1.0;
    MFloat frequency = 1.0;
  } m_input2Aca, m_aca2Output;

  /// Number of surfaces
  MInt m_noSurfaces = -1;
  /// Ids of the surface input file
  std::vector<MInt> m_surfaceIds{};

  MBool m_useMergedInputFile;
  MInt m_inputFileIndexStart;
  MInt m_inputFileIndexEnd;

  /// Number of (local) surface elements for each surface
  std::vector<MInt> m_noSurfaceElements{};

  /// Allow multiple surfaces per rank
  MBool m_allowMultipleSurfacesPerRank = false;

  /// Total number of surface elements for each surface
  std::vector<MInt> m_noSurfaceElementsGlobal{};

  /// Local surface element offsets for each surface
  std::vector<MInt> m_surfaceElementOffsets{};

  /// Surface local MPI communicators
  std::vector<MPI_Comm> m_mpiCommSurface{};
  /// Surface local number of domains
  std::vector<MInt> m_noDomainsSurface{};
  /// Surface local domain id
  std::vector<MInt> m_domainIdSurface{};

  /// Weighting factor for each surface (for surface averaging)
  std::vector<MFloat> m_surfaceWeightingFactor{};

  /// Surface input file names that were used for sampling
  std::vector<MString> m_surfaceInputFileName{};

  /// Number of samples
  MInt m_noSamples = -1;
  /// Total number of samples in input data files
  MInt m_totalNoSamples = -1;
  /// Stride of used samples in input data files
  MInt m_sampleStride = 1;
  /// Offset in input data files from which to start reading samples
  MInt m_sampleOffset = 0;

  /// Collector for surface data
  SurfaceDataCollector m_surfaceData;

  /// Acoustic extrapolation method to use, see enum ExtrapolationMethods in enums.h
  MInt m_acousticMethod = -1;
  MInt m_fwhTimeShiftType = 0;
  MInt m_fwhMeanStateType = 0;

  MInt m_noVariables = -1;
  MInt m_noComplexVariables = -1;

  /// Set of input variables
  std::map<MString, MInt> m_inputVars{};

  /// Window type
  MInt m_windowType = 0;
  /// Normalization factor of window
  MFloat m_windowNormFactor = -1.0;

  /// Fourier type
  MInt m_transformationType = 0;

  /// Number of samples
  MInt m_noGlobalObservers = -1;

  std::vector<MFloat> m_globalObserverCoordinates; /// Coordinates of (global) observer points
  ObserverDataCollector m_observerData;            /// Collector for (local) observer data

  MInt m_offsetObserver = 0; /// Offset of local observer
  MInt m_noObservers = 0;    /// Number of local observer

  /// Coordinate system transformation
  std::array<MFloat, nDim * nDim> m_transformationMatrix;
  MBool m_zeroFlowAngle = false;

  /// I/O
  MString m_inputDir;              /// Input directory
  MString m_outputFilePrefix = ""; /// Prefix for all output files

  /// Enable storing of surface data (testing or output of generated data)
  MBool m_storeSurfaceData = false;

  /// File containing observer points
  MString m_observerFileName;
  /// Enable observer point generation (instead of reading from file)
  MBool m_generateObservers = false;


  /// Enable generation of analytical surface data
  MBool m_generateSurfaceData = false;

  /// Source Type
  MInt m_sourceType = -1;

  /// Parameter vector
  maia::acoustic::SourceParameters m_sourceParameters;

  /// Omega factor
  MFloat m_omega_factor = -1.0;

  /// Wave length lambda
  MFloat m_lambdaZero = -1.0;
  MFloat m_lambdaMach = -1.0;

  /// Non-dimensionalization variables
  MFloat m_cInfty = -1.0;
  MFloat m_rhoInfty = -1.0;

  /// Analytic Fourier Transformation, also used for Validation
  std::vector<MFloat> m_Analyticfreq{};

  /// Analytic rms pressure for accuracy function
  std::vector<MFloat> m_rmsP_Analytic{};

  /// Constant time step size
  MFloat m_dt = -1.0;

  /// Mach vector
  MFloat m_MaVec[3];
  /// Mach number
  MFloat m_Ma = -1.0;
  /// Mach number used for change of nondimensionalization
  MFloat m_MaDim = -1.0;

  static constexpr MFloat m_gamma = 1.4; // isentropic exponent

  /// Can FFT be used?
  MBool m_FastFourier = false;

  /// Storage for times and frequencies
  std::vector<MFloat> m_times{};
  std::vector<MFloat> m_frequencies{};
  MBool m_hasTimes = false;

  /// Post Processing of observer signals
  MInt m_noPostprocessingOps = 0;
  std::vector<MInt> m_postprocessingOps{};

  MBool m_acaPostprocessingMode = false;
  MString m_acaPostprocessingFile = "";

  // Timers
  // Timer group which holds all solver-wide timers
  MInt m_timerGroup = -1;
  // Stores all solver-wide timers
  std::array<MInt, Timers::_count> m_timers{};

  //////////////////////////////////////////////////////////////////////////////
  // Hold indices for variables for different formulations
  struct FWH {
    // Flow variables: u, v, w, rho, p
    static constexpr const MInt U = 0;
    static constexpr const MInt V = 1;
    static constexpr const MInt W = nDim - 1;
    static constexpr const MInt RHO = nDim;
    static constexpr const MInt P = nDim + 1;

    // Source terms
    static constexpr const MInt FX = nDim + 2;
    static constexpr const MInt FY = nDim + 3;
    static constexpr const MInt FZ = 2 * nDim + 1;
    static constexpr const MInt Q = 2 * nDim + 2;
    // Number of variables
    static constexpr const MInt noVars = 2 * nDim + 3;

    // Complex variables
    static constexpr const MInt FX_C = 0;
    static constexpr const MInt FY_C = 1;
    static constexpr const MInt FZ_C = nDim - 1;
    static constexpr const MInt Q_C = nDim;
    // Number of complex variables
    static constexpr const MInt noComplexVars = nDim + 1;
  };

  struct FWH_TIME {
    // Flow variables: u, v, w, rho, p
    static constexpr const MInt U = 0;
    static constexpr const MInt V = 1;
    static constexpr const MInt W = nDim - 1;
    static constexpr const MInt RHO = nDim;
    static constexpr const MInt P = nDim + 1;

    // Surface velocities
    static constexpr const MInt surfMX = nDim + 2;
    static constexpr const MInt surfMY = nDim + 3;
    static constexpr const MInt surfMZ = 2 * nDim + 1;

    // Source terms
    static constexpr const MInt srcUn = 2 * nDim + 2;
    static constexpr const MInt srcLX = 2 * nDim + 3;
    static constexpr const MInt srcLY = 2 * nDim + 4;
    static constexpr const MInt srcLZ = 3 * nDim + 2;

    // Source-Time FW-H varaiables
    static constexpr const MInt timeObs = 3 * nDim + 3;
    static constexpr const MInt fwhPP = 3 * nDim + 4;

    // Number of variables
    static constexpr const MInt noVars = 3 * nDim + 5;
  };

  struct FWH_APE {
    // Perturbed APE variables: u', v', w', p'
    static constexpr const MInt UP = 0;
    static constexpr const MInt VP = 1;
    static constexpr const MInt WP = nDim - 1;
    static constexpr const MInt PP = nDim;
    // Source terms
    static constexpr const MInt FX = nDim + 1;
    static constexpr const MInt FY = nDim + 2;
    static constexpr const MInt FZ = 2 * nDim;
    static constexpr const MInt Q = 2 * nDim + 1;
    // Number of variables
    static constexpr const MInt noVars = 2 * nDim + 2;

    // Complex variables
    static constexpr const MInt FX_C = 0;
    static constexpr const MInt FY_C = 1;
    static constexpr const MInt FZ_C = nDim - 1;
    static constexpr const MInt Q_C = nDim;
    // Number of complex variables
    static constexpr const MInt noComplexVars = nDim + 1;
  };
  //////////////////////////////////////////////////////////////////////////////

  // Hold indices for different postprocessing operations
  struct PP {
    static constexpr const MInt RMS_PRESSURE = 0;
    static constexpr const MInt OASPL = 1;
    static constexpr const MInt SPL = 2;
    static constexpr const MInt pABS = 3;
    static constexpr const MInt calcObsPress = 4;
  };

  // Hold indices for different window types
  struct WINDOW {
    static constexpr const MInt NONE = 0;
    static constexpr const MInt HANNING = 1;
    static constexpr const MInt HAMMING = 2;
    static constexpr const MInt MODHANNING = 3;

    static const std::array<MString, 4> windowNames;
  };

  // Hold indices for different transformation types
  struct FOURIER {
    static constexpr const MInt DFT = 0;
    static constexpr const MInt FFT = 1;
  };

  /// Hold data for symmetry boundary condition
  struct AcaSymBc {
    std::array<MFloat, nDim> origin;
    std::array<MFloat, nDim> normal;
  };
  std::vector<AcaSymBc> m_symBc;
};

#endif // ACASOLVER_H_
