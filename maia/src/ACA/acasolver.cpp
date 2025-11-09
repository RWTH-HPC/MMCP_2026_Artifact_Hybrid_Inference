// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "acasolver.h"
#include "INCLUDE/maiaconstants.h"
#include "INCLUDE/maiatypes.h"
#include "IO/context.h"
#include "MEMORY/scratch.h"
#include "UTIL/functions.h"
#include "UTIL/parallelfor.h"
#include "UTIL/tensor.h"
#include "acapostprocessing.hpp"
#include "typetraits.h"

using namespace maia;

template <MInt nDim>
const std::array<MString, 4> AcaSolver<nDim>::WINDOW::windowNames = {"none", "Hanning", "Hamming", "modified Hanning"};

template <MInt nDim>
AcaSolver<nDim>::AcaSolver(const MInt solverId, const MPI_Comm comm) : Solver(solverId, comm, true) {}

/// \brief Initialize solver.
template <MInt nDim>
void AcaSolver<nDim>::initSolver() {
  TRACE();

  // Initialize all solver-wide timers
  initTimers();
  // Input/output
  setInputOutputProperties();
  // Numerical method
  setNumericalProperties();
}

/// \brief Initialize solver timers.
template <MInt nDim>
void AcaSolver<nDim>::initTimers() {
  TRACE();

  // Invalidate timer ids
  std::fill(m_timers.begin(), m_timers.end(), -1);

  // Create timer group & timer for solver, and start the timer
  NEW_TIMER_GROUP_NOCREATE(m_timerGroup, "AcaSolver (solverId = " + std::to_string(m_solverId) + ")");
  NEW_TIMER_NOCREATE(m_timers[Timers::SolverType], "total object lifetime", m_timerGroup);
  RECORD_TIMER_START(m_timers[Timers::SolverType]);

  // Create regular solver-wide timers
  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::Run], "run", m_timers[Timers::SolverType]);

  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::InitParallelization], "initParallelization", m_timers[Timers::Run]);
  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::ReadSurfaceData], "readSurfaceData", m_timers[Timers::Run]);

  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::InitObservers], "initObservers", m_timers[Timers::Run]);
  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::CalcSourceTerms], "calcSourceTerms", m_timers[Timers::Run]);
  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::WindowAndTransformSources], "windowAndTransformSources",
                         m_timers[Timers::Run]);
  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::CalcSurfaceIntegrals], "calcSurfaceIntegrals", m_timers[Timers::Run]);
  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::CombineSurfaceIntegrals], "combineSurfaceIntegrals", m_timers[Timers::Run]);
  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::BacktransformObservers], "backtransformObservers", m_timers[Timers::Run]);
  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::SaveObserverSignals], "saveObserverSignals", m_timers[Timers::Run]);
  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::Postprocessing], "postprocessing", m_timers[Timers::Run]);
}

/** \brief  Average all FWH timer over participating ranks
 *  \author Miro Gondrum
 *  \date   21.02.2024
 */
template <MInt nDim>
void AcaSolver<nDim>::averageTimer() {
  // 0) map timer ids for safety
  std::vector<MInt> timerIds_;
  timerIds_.reserve(11);
  timerIds_.emplace_back(m_timers[Timers::Run]);
  timerIds_.emplace_back(m_timers[Timers::InitParallelization]);
  timerIds_.emplace_back(m_timers[Timers::ReadSurfaceData]);
  timerIds_.emplace_back(m_timers[Timers::InitObservers]);
  timerIds_.emplace_back(m_timers[Timers::CalcSourceTerms]);
  timerIds_.emplace_back(m_timers[Timers::WindowAndTransformSources]);
  timerIds_.emplace_back(m_timers[Timers::CalcSurfaceIntegrals]);
  timerIds_.emplace_back(m_timers[Timers::CombineSurfaceIntegrals]);
  timerIds_.emplace_back(m_timers[Timers::BacktransformObservers]);
  timerIds_.emplace_back(m_timers[Timers::SaveObserverSignals]);
  timerIds_.emplace_back(m_timers[Timers::Postprocessing]);
  const MInt noTimers = timerIds_.size();
  // 1) fill buffer with local timer values
  std::vector<MFloat> timerValues_;
  timerValues_.reserve(noTimers);
  for(MInt i = 0; i < noTimers; i++) {
    timerValues_.emplace_back(RETURN_TIMER_TIME(timerIds_[i]));
  }
  // 2) collect values from all ranks and use max value
  MPI_Allreduce(MPI_IN_PLACE, timerValues_.data(), noTimers, maia::type_traits<MFloat>::mpiType(), MPI_MAX, mpiComm(),
                AT_, "MPI_IN_PLACE", "timerValues_");
  for(MInt i = 0; i < noTimers; i++) {
    SET_RECORD(timerIds_[i], timerValues_[i]);
  }
}


/// \brief Read/set all I/O related properties
template <MInt nDim>
void AcaSolver<nDim>::setInputOutputProperties() {
  TRACE();

  /*! \property
    \page propertyPageACA ACA
    \section inputDir
    <code>MString AcaSolver::m_inputDir</code>\n
    Input directory for reading surface data.\n\n
    Keywords: <i>ACA, input</i>
  */
  m_inputDir = Context::getSolverProperty<MString>("inputDir", m_solverId, AT_);
  m_inputDir = testcaseDir() + m_inputDir;

  /*! \property
    \page propertyPageACA ACA
    \section observerFileName
    <code>MString AcaSolver::m_observerFileName</code>\n
    default = <code>""</code>\n\n
    File name which specify observer location in a space separator ed ASCII file.
    Keywords: <i>ACA, observer</i>
  */
  m_observerFileName = "";
  m_observerFileName = Context::getSolverProperty<MString>("observerFileName", m_solverId, AT_, &m_observerFileName);

  /*! \property
    \page propertyPageACA ACA
    \section generateObservers
    <code>MBool AcaSolver::m_generateObservers</code>\n
    default = <code>false</code>\n\n
    Flag whether observer location are generated or read from oberverFileName.\n
    Keywords: <i>ACA, observer</i>
  */
  m_generateObservers = false;
  m_generateObservers = Context::getSolverProperty<MBool>("generateObservers", m_solverId, AT_, &m_generateObservers);

  TERMM_IF_COND(m_observerFileName == "" && m_generateObservers == false,
                "No observer file provided and generation of observers disabled.");

  // Prefix added to all output files
  /*! \property
    \page propertyPageACA ACA
    \section outputFilePrefix
    <code>MString AcaSolver::m_outputFilePrefix</code>\n
    default = <code>""</code>\n\n
    Prefix for output data.\n
    Keywords: <i>ACA, output</i>
  */
  m_outputFilePrefix = "";
  m_outputFilePrefix = Context::getSolverProperty<MString>("outputFilePrefix", m_solverId, AT_, &m_outputFilePrefix);

  // Number of surfaces
  /*! \property
    \page propertyPageACA ACA
    \section noSurfaces
    <code>MInt AcaSolver::m_noSurfaces</code>\n
    default = <code>0</code>\n\n
    Number of permeable surface data files.\n
    Keywords: <i>ACA, input</i>
  */
  m_noSurfaces = 0;
  m_noSurfaces = Context::getSolverProperty<MInt>("noSurfaces", m_solverId, AT_, &m_noSurfaces);

  // Surface ids
  if(Context::propertyExists("surfaceIds", m_solverId)) {
    const MInt noGivenSrfcIds = Context::propertyLength("surfaceIds", m_solverId);
    if(noSurfaces() != 0 && noSurfaces() != noGivenSrfcIds) {
      std::cerr << "WARNING: Mismatch between properties 'noSurfaces' and size of surfaceIds" << std::endl;
    }
    m_noSurfaces = noGivenSrfcIds;
    TERMM_IF_COND(m_noSurfaces <= 0, "Number of surfaces (noSurfaces) needs to be > 0");
    m_surfaceIds.resize(noSurfaces());
    for(MInt id = 0; id < noSurfaces(); id++) {
      /*! \property
        \page propertyPageACA ACA
        \section surfaceIds
        <code>MBool AcaSolver::m_surfaceIds</code>\n
        default = <code>[0,1,2,..,n-1]</code>\n\n
        Surface ids to be used for reading surface data file. Default is a
        consecutive ordering.\n
        Keywords: <i>ACA, input</i>
      */
      m_surfaceIds[id] = Context::getSolverProperty<MInt>("surfaceIds", m_solverId, AT_, id);
    }
  } else {
    TERMM_IF_COND(m_noSurfaces <= 0, "Number of surfaces (noSurfaces) needs to be > 0");
    m_surfaceIds.resize(noSurfaces());
    for(MInt id = 0; id < noSurfaces(); id++) {
      m_surfaceIds[id] = id;
    }
  }

  // Use merged input file or multiple input files
  /*! \property
    \page propertyPageACA ACA
    \section useMergedInputFile
    <code>MBool AcaSolver::m_useMergedInputFile</code>\n
    default = <code>true</code>\n\n
    Flag whether input surface data is in a merged file (true) or in the non
    merged format as written by the post processing class (false).\n
    Keywords: <i>ACA, input</i>
  */
  m_useMergedInputFile = true;
  m_useMergedInputFile =
      Context::getSolverProperty<MBool>("useMergedInputFile", m_solverId, AT_, &m_useMergedInputFile);
  if(m_useMergedInputFile == false) {
    // inputFileIndexStart and inputFileIndexEnd
    /*! \property
      \page propertyPageACA ACA
      \section inputFileIndexStart
      <code>MInt AcaSolver::m_inputFileIndexStart</code>\n\n
      For useMergedInputFile=true provides the start file index of the surface file to be read.\n
      Keywords: <i>ACA, input</i>
    */
    m_inputFileIndexStart = Context::getSolverProperty<MInt>("inputFileIndexStart", m_solverId, AT_);
    /*! \property
      \page propertyPageACA ACA
      \section inputFileIndexEnd
      <code>MInt AcaSolver::m_inputFileIndexEnd</code>\n\n
      For useMergedInputFile=true provides the end file index of the surface file to be read.\n
      Keywords: <i>ACA, input</i>
    */
    m_inputFileIndexEnd = Context::getSolverProperty<MInt>("inputFileIndexEnd", m_solverId, AT_);
    TERMM_IF_COND(m_inputFileIndexStart > m_inputFileIndexEnd,
                  "'inputFileIndexStart' must not be greater than 'm_inputFileIndexEnd'.");
  }

  /*! \property
    \page propertyPageACA ACA
    \section acaNoSamples
    <code>MInt AcaSolver::m_noSamples</code>\n
    default = <code>-1</code>\n\n
    User specified number of samples to use for the computation (`-1` means
    all samples in data file).\n
    Keywords: <i>ACA, input</i>
  */
  m_noSamples = -1;
  m_noSamples = Context::getSolverProperty<MInt>("acaNoSamples", m_solverId, AT_, &m_noSamples);

  /*! \property
    \page propertyPageACA ACA
    \section acaSampleStride
    <code>MInt AcaSolver::m_sampleStride</code>\n
    default = <code>1</code>\n\n
    Stride used to read samples, i.e., read only every n-th sample.\n
    Keywords: <i>ACA, input</i>
  */
  m_sampleStride = 1;
  m_sampleStride = Context::getSolverProperty<MInt>("acaSampleStride", m_solverId, AT_, &m_sampleStride);

  /*! \property
    \page propertyPageACA ACA
    \section acaSampleOffset
    <code>MInt AcaSolver::m_sampleOffset</code>\n
    default = <code>0</code>\n\n
    Offset used to read samples.\n
    Keywords: <i>ACA, input</i>
  */
  m_sampleOffset = 0;
  m_sampleOffset = Context::getSolverProperty<MInt>("acaSampleOffset", m_solverId, AT_, &m_sampleOffset);

  /*! \property
    \page propertyPageACA ACA
    \section m_allowMultipleSurfacesPerRank
    <code>MBool m_allowMultipleSurfacesPerRank</code>\n
    default = <code>false</code>\n\n
    By default each rank hold data of one surface. This does not make sense if
    noSurfaces > noRanks, which is the purpose of this property.\n
    Keywords: <i>ACA</i>
  */
  m_allowMultipleSurfacesPerRank = Context::getSolverProperty<MBool>("allowMultipleSurfacesPerRank", m_solverId, AT_,
                                                                     &m_allowMultipleSurfacesPerRank);

  /*! \property
    \page propertyPageACA ACA
    \section generateSurfaceData
    <code>MBool AcaSolver::m_generateSurfaceData</code>\n
    default = <code>false</code>\n\n
    Flag whether surface data is generated from analytical function.\n
    Keywords: <i>ACA, input</i>
  */
  m_generateSurfaceData = false;
  m_generateSurfaceData =
      Context::getSolverProperty<MBool>("generateSurfaceData", m_solverId, AT_, &m_generateSurfaceData);

  /*! \property
    \page propertyPageACA ACA
    \section acaResultsNondimMode
    <code>MBool AcaSolver::m_acaResultsNondimMode</code>\n
    default = <code>0</code>\n\n
    Possible values are: \n
    <ul>
    <li><code>0</code> non-dim ACA</li>
    <li><code>1</code> non-dim stagnation</li>
    <li><code>2</code> non-dim LB</li>
    </ul>
    Keywords: <i>ACA, output</i>
  */
  m_acaResultsNondimMode = 0;
  m_acaResultsNondimMode =
      Context::getSolverProperty<MInt>("acaResultsNondimMode", m_solverId, AT_, &m_acaResultsNondimMode);

  /*! \property
    \page propertyPageACA ACA
    \section storeSurfaceData
    <code>MBool AcaSolver::m_storeSurfaceData</code>\n
    default = <code>false</code>\n\n
    Check if (generated/loaded) surface should be written to a file (for testcases) .\n
    Keywords: <i>ACA, output</i>
  */
  m_storeSurfaceData = false;
  m_storeSurfaceData = Context::getSolverProperty<MBool>("storeSurfaceData", m_solverId, AT_, &m_storeSurfaceData);

  // Read in the postprocessing operations
  m_noPostprocessingOps = 0;
  if(Context::propertyExists("postprocessingOps", m_solverId)) {
    m_noPostprocessingOps = Context::propertyLength("postprocessingOps", m_solverId);
    m_postprocessingOps.resize(m_noPostprocessingOps, -1);
  }
  for(MInt i = 0; i < m_noPostprocessingOps; i++) {
    /*! \property
      \page propertyPageACA ACA
      \section postprocessingOps
      <code>MInt AcaSolver::m_postprocessingOps</code>\n
      default = <code>false</code>\n\n
      Id of the used post processing operations.\n
      Possible values are: \n
      <ul>
      <li><code>0</code> calculate RMS pressure</li>
      <li><code>1</code> calculate OASPL</li>
      <li><code>2</code> calculate SPL</li>
      <li><code>3</code> calculate absoulute pressure</li>
      <li><code>4</code> calculate observer pressure</li>
      </ul>
      Keywords: <i>ACA, postprocessing</i>
    */
    m_postprocessingOps[i] = Context::getSolverProperty<MInt>("postprocessingOps", m_solverId, AT_, i);
  }

  /*! \property
    \page propertyPageACA ACA
    \section acaPostprocessingMode
    <code>MBool AcaSolver::m_acaPostprocessingMode</code>\n
    default = <code>false</code>\n\n
    Flag if only observer data from a previous run is read instead of performing
    a new simulation run, i.e., running in a only post processing mode.\n
    Keywords: <i>ACA, postprocessing</i>
  */
  m_acaPostprocessingMode = false;
  m_acaPostprocessingMode =
      Context::getSolverProperty<MBool>("acaPostprocessingMode", m_solverId, AT_, &m_acaPostprocessingMode);

  if(m_acaPostprocessingMode) {
    TERMM_IF_COND(m_noPostprocessingOps == 0, "acaPostprocessingMode is enabled but noPostprocessingOps is zero.");

    /*! \property
      \page propertyPageACA ACA
      \section acaPostprocessingFile
      <code>MString AcaSolver::m_acaPostprocessingFile</code>\n\n
      File name where input observer data are given.\n
      Keywords: <i>ACA, postprocessing</i>
    */
    m_acaPostprocessingFile = Context::getSolverProperty<MString>("acaPostprocessingFile", m_solverId, AT_);
  }

  /*! \property
    \page propertyPageACA ACA
    \section symmetryOrigin
    This property defines the origins of symmetry planes.
    <ul>
    <li><code>Any float array with length of multiple nDim</code> </li>
    </ul>
    Keywords: <i>ACA, BC, symmetry</i>
  */
  const MString propNameSymOrigin = "symmetryOrigin";
  /*! \property
    \page propertyPageACA ACA
    \section symmetryNormal
    This property defines the normals of symmetry planes.
    <ul>
    <li><code>Any float array with length of multiple nDim</code> </li>
    </ul>
    Keywords: <i>ACA, BC, symmetry</i>
  */
  const MString propNameSymNormal = "symmetryNormal";
  if(!(Context::propertyExists(propNameSymOrigin, m_solverId)
       || Context::propertyExists(propNameSymNormal, m_solverId)))
    return;
  const MInt propLength_symOrign = Context::propertyLength(propNameSymOrigin, m_solverId);
  const MInt propLength_symNormal = Context::propertyLength(propNameSymNormal, m_solverId);
  TERMM_IF_COND(propLength_symOrign % nDim != 0, propNameSymOrigin + " has to be of a length of multiple nDim.");
  TERMM_IF_COND(propLength_symNormal % nDim != 0, propNameSymNormal + " has to be of a length of multiple nDim.");
  const MInt noSym = propLength_symOrign / nDim;
  TERMM_IF_COND(noSym != propLength_symNormal / nDim,
                propNameSymOrigin + " and " + propNameSymNormal + " have to be of the same length.");
  // TODO labels:ACA Make more than one symmetry possible?
  TERMM_IF_COND(noSym > 1, "So far not more than 1 symmetry BC possible.");
  m_symBc.resize(noSym);
  for(MInt i = 0; i < noSym; i++) {
    MFloat lengthSymNormal = 0.0;
    for(MInt d = 0; d < nDim; d++) {
      m_symBc[i].origin[d] = Context::getSolverProperty<MFloat>(propNameSymOrigin, m_solverId, AT_, i * nDim + d);
      m_symBc[i].normal[d] = Context::getSolverProperty<MFloat>(propNameSymNormal, m_solverId, AT_, i * nDim + d);
      lengthSymNormal += POW2(m_symBc[i].normal[d]);
    }
    lengthSymNormal = sqrt(lengthSymNormal);
    for(MInt d = 0; d < nDim; d++) {
      m_symBc[i].normal[d] /= lengthSymNormal; // normalized normal vector
    }
  }
}


/// \brief Read/set all numerical properties
template <MInt nDim>
void AcaSolver<nDim>::setNumericalProperties() {
  TRACE();

  /*! \property
    \page propertyPageACA ACA
    \section acousticMethod
    <code>MInt AcaSolver::acousticMethod</code>\n
    Set method to be used for acoustic far-field prediction.\n
    Possible values are:\n
    <ul>
    <li><code>"FWH_METHOD"</code> Input data describes flow field \f$(\mv{u},\rho,p)\f$ </li>
    <li><code>"FWH_APE_METHOD"</code> Input data describes acoustic field \f$(\mv{u}^\prime,p^\prime)\f$ </li>
    </ul>
    Keywords: <i>ACA, numerical method</i>
  */
  const MString methodName = Context::getSolverProperty<MString>("acousticMethod", solverId(), AT_);
  m_acousticMethod = -1;
  m_acousticMethod = string2enum(methodName);
  m_inputVars.clear();

  switch(m_acousticMethod) {
    case FWH_METHOD:
      // Set number of (complex) variables
      m_noVariables = FWH::noVars;
      m_noComplexVariables = FWH::noComplexVars;

      // Store names of input variables with corresponding storage position in m_surfaceData
      m_inputVars["u"] = FWH::U;
      m_inputVars["v"] = FWH::V;
      if constexpr(nDim == 3) {
        m_inputVars["w"] = FWH::W;
      }
      m_inputVars["rho"] = FWH::RHO;
      m_inputVars["p"] = FWH::P;

      break;
    case FWH_APE_METHOD:
      // Set number (complex) variables
      m_noVariables = FWH_APE::noVars;
      m_noComplexVariables = FWH_APE::noComplexVars;

      // Store names of input variables with corresponding storage position in m_surfaceData
      m_inputVars["u"] = FWH_APE::UP;
      m_inputVars["v"] = FWH_APE::VP;
      if constexpr(nDim == 3) {
        m_inputVars["w"] = FWH_APE::WP;
      }
      m_inputVars["p"] = FWH_APE::PP;
      // APE = Acoustic Perturbed Equation

      break;
    default:
      TERMM(1, "Invalid acoustic extrapolation method: '" + std::to_string(m_acousticMethod) + "'");
      break;
  }

  // FWH_Farassat_1A time domain formulation
  if(string2enum(solverMethod()) == MAIA_FWH_TIME) m_noVariables = FWH_TIME::noVars;

  // Use observer time shift
  /*! \property
    \page propertyPageACA ACA
    \section
    <code>MInt AcaSolver::m_fwhTimeShiftType</code>\n
    Possible values are: \n
    <ul>
    <li><code>0</code> Zero time shift  </li>
    <li><code>1</code> Local time shift for each observer </li>
    <li><code>2</code> Global time shift for all observers </li>
    </ul>
    Define the constant time shift between source time and observer time,
    to make good use of the data. timeObserver = timeSource + delta_t.
    To plot the p' contour at a certain time, Zero or Global time shift
    for all observers is required.\n
    Keywords: <i>ACA, time-domain formulation</i>
  */
  m_fwhTimeShiftType = 0;
  m_fwhTimeShiftType = Context::getSolverProperty<MInt>("fwhTimeShiftType", m_solverId, AT_, &m_fwhTimeShiftType);

  // Set mean variables
  /*! \property
    \page propertyPageACA ACA
    \section
    <code>MInt AcaSolver::m_fwhMeanStateType</code>\n
    Possible values are: \n
    <ul>
    <li><code>0</code> use far-field flow states as mean properties, e.g. uMean = U_inf </li>
    <li><code>1</code> use local time-averaged flow properties as mean properties, e.g. uMean = <u> </li>
    </ul>
    Keywords: <i>ACA, numerical method</i>
  */
  m_fwhMeanStateType = 0;
  m_fwhMeanStateType = Context::getSolverProperty<MInt>("fwhMeanStateType", m_solverId, AT_, &m_fwhMeanStateType);

  // Type used at windowing function: default setting is 0 = NONE
  /*! \property
    \page propertyPageACA ACA
    \section windowType
    <code>MInt AcaSolver::m_windowType</code>\n
    default = <code>0</code>\n\n
    The windowing function used to weighting the data before applying Fourier
    transformation. If the input signal is not an exact integer period this
    helps reducing spectral leakage.\n
    Possible values are: \n
    <ul>
    <li><code>0</code> None             </li>
    <li><code>1</code> Hanning          </li>
    <li><code>2</code> Hamming          </li>
    <li><code>3</code> Modified Hamming </li>
    </ul>
    Keywords: <i>ACA, numerical method</i>
  */
  m_windowType = 0;
  m_windowType = Context::getSolverProperty<MInt>("windowType", m_solverId, AT_, &m_windowType);

  /*! \property
    \page propertyPageACA ACA
    \section transformationType
    <code>MInt AcaSolver::m_transformationType</code>\n
    default = <code>1</code>\n\n
    Flag if Fast Fourier transformation (FFT) shall be used. Later it will be
    checked if this possible as FFT is only possible for a number samples being
    a potency of 2.\n
    Possible values are: \n
    <ul>
    <li><code>0</code> discrete </li>
    <li><code>1</code> fast     </li>
    </ul>
    Keywords: <i>ACA, numerical method</i>
  */
  m_transformationType = 1;
  m_transformationType = Context::getSolverProperty<MInt>("transformationType", m_solverId, AT_, &m_transformationType);

  // Surface weighting factor for surface averaging (default 1.0)
  m_surfaceWeightingFactor.resize(noSurfaces(), 1.0);
  if(Context::propertyExists("surfaceWeightingFactor")) {
    for(MInt sId = 0; sId < noSurfaces(); sId++) {
      /*! \property
        \page propertyPageACA ACA
        \section surfaceWeightingFactor
        <code>MFloat weightingFactor</code>\n
        \n
        Weight applied to the source terms of the corresponding surface.\n
        Keywords: <i>ACA, numerical method</i>
      */
      const MFloat weightingFactor = Context::getSolverProperty<MFloat>("surfaceWeightingFactor", solverId(), AT_, sId);
      TERMM_IF_COND(weightingFactor < 0.0 || weightingFactor > 1.0,
                    "Error: surface weighting factor should be in range (0.0,1.0], is: "
                        + std::to_string(weightingFactor));
      m_surfaceWeightingFactor[sId] = weightingFactor;

      m_log << "Surface weighting factor for surface #" << m_surfaceIds[sId] << ": " << m_surfaceWeightingFactor[sId]
            << std::endl;
    }
  } else {
    m_log << "Surface weighting not enabled." << std::endl;
  }

  m_log << "Mach vector: ";
  MFloat machSq = 0.0;
  for(MInt i = 0; i < nDim; i++) {
    m_MaVec[i] = 0.0;
    /*! \property
      \page propertyPageACA ACA
      \section Ma
      <code>MFloat AcaSolver::m_MaVec[3]</code>\n
      default = <code>[0.0, 0.0, 0.0]</code>\n\n
      Free stream Mach number in vector notation.\n
      Keywords: <i>ACA, flow variables</i>
    */
    m_MaVec[i] = Context::getSolverProperty<MFloat>("Ma", m_solverId, AT_, &m_MaVec[i], i);
    m_log << m_MaVec[i] << " ";
    machSq += POW2(m_MaVec[i]);
  }
  const MFloat Mach = sqrt(machSq);
  m_Ma = Mach;

  /*! \property
    \page propertyPageACA ACA
    \section MaDim
    <code>MFloat AcaSolver::m_MaDim</code>\n
    default = Taking absolute value of Ma in vector notation.\n\n
    Mach number for change of non-dimensionalization (e.g. for a jet).\n
    Keywords: <i>ACA, flow variables</i>
  */
  if(Context::propertyExists("MaDim", m_solverId)) {
    m_MaDim = Context::getSolverProperty<MFloat>("MaDim", m_solverId, AT_);
  } else {
    m_MaDim = m_Ma;
  }
  m_log << "; Ma = " << m_Ma << "; Ma_dim = " << m_MaDim << std::endl;

  // Coordinate system transformation
  /* Project the vectors to the flow direction U_inf = (U_inf, 0, 0) to apply FWH formulations.
   * Used in functions: transformCoordinatesToFlowDirection2D, transformCoordinatesToFlowDirection3D
   * transformBackToOriginalCoordinate2D, transformBackToOriginalCoordinate3D
   */
  MFloat theta = 0.0;
  MFloat alpha = 0.0;
  if constexpr(nDim == 3) {
    if(!approx(m_MaVec[2], 0.0, MFloatEps) && !approx(m_MaVec[0], 0.0, MFloatEps)) {
      theta = atan2(m_MaVec[2], m_MaVec[0]);
    }
    if(!approx(m_MaVec[1], 0.0, MFloatEps)) {
      alpha = asin(m_MaVec[1] / m_Ma);
    }
    m_transformationMatrix[0] = cos(alpha) * cos(theta);
    m_transformationMatrix[1] = sin(alpha);
    m_transformationMatrix[2] = cos(alpha) * sin(theta);
    m_transformationMatrix[3] = -sin(alpha) * cos(theta);
    m_transformationMatrix[4] = cos(alpha);
    m_transformationMatrix[5] = -sin(alpha) * sin(theta); // "+" in Locard's formulation
    m_transformationMatrix[6] = -sin(theta);
    m_transformationMatrix[7] = 0.0;
    m_transformationMatrix[8] = cos(theta);
  } else if constexpr(nDim == 2) {
    if(!approx(m_MaVec[1], 0.0, MFloatEps) && !approx(m_MaVec[0], 0.0, MFloatEps)) {
      theta = atan2(m_MaVec[1], m_MaVec[0]);
    }
    m_transformationMatrix[0] = cos(theta);
    m_transformationMatrix[1] = sin(theta);
    m_transformationMatrix[2] = -sin(theta);
    m_transformationMatrix[3] = cos(theta);
  }
  MFloat flowAngle = std::sqrt(theta * theta + alpha * alpha);
  if(approx(flowAngle, 0.0, MFloatEps)) m_zeroFlowAngle = true;

  if(m_generateSurfaceData) {
    /*! \property
      \page propertyPageACA ACA
      \section sourceType
      <code>MInt AcaSolver::m_sourceType</code>\n
      default = <code>false</code>\n\n
      Type of generated source, e.g. Monopole, Dipole and so on. Used in
      `generateSurfaceData()` and `calcObsPressAnalytic()`.\n
      Possible values are: \n
      <ul>
      <li><code>0</code> Analytic Monopole    </li>
      <li><code>1</code> Analytic Dipole      </li>
      <li><code>2</code> Analytic Quadrupole  </li>
      <li><code>3</code> Inviscid vortex convection </li>
      </ul>
      Keywords: <i>ACA, input</i>
    */
    m_sourceType = -1;
    m_sourceType = Context::getSolverProperty<MInt>("sourceType", m_solverId, AT_, &m_sourceType);

    m_sourceParameters.perturbed = (m_acousticMethod == FWH_METHOD) ? false : true;

    // Mach number from -x direction
    m_sourceParameters.Ma = Mach;
    constexpr MFloat c0 = 1.0; // per definition

    /*! \property
      \page propertyPageACA ACA
      \section rho0
      <code>MFloat AcaSolver::m_sourceParameters.rho0</code>\n
      default = <code>1.0</code>\n\n
      Ambient density for source parameters.\n
      Keywords: <i>ACA, input</i>
    */
    m_sourceParameters.rho0 = 1.0;
    m_sourceParameters.rho0 = Context::getSolverProperty<MFloat>("rho0", m_solverId, AT_, &m_sourceParameters.rho0);

    // values
    /*! \property
      \page propertyPageACA ACA
      \section omega_factor
      <code>MFloat AcaSolver::m_omega_factor</code>\n
      default = <code>2D: 4.0, 3D: 4.0/46.0</code>\n\n
      Omega factor used for source parameters.\n
      Keywords: <i>ACA, input</i>
    */
    m_omega_factor = (nDim == 2) ? 4.0 : 4.0 / 46.0;
    m_omega_factor = Context::getSolverProperty<MFloat>("omega_factor", m_solverId, AT_, &m_omega_factor);

    /*! \property
      \page propertyPageACA ACA
      \section omega
      <code>MFloat AcaSolver::m_sourceParameters.omega</code>\n
      default = <code> m_omega_factor * PI * c0;
      Omega for source parameters.\n
      Keywords: <i>ACA, input</i>
    */
    m_sourceParameters.omega = m_omega_factor * PI * c0;
    m_sourceParameters.omega = Context::getSolverProperty<MFloat>("omega", m_solverId, AT_, &m_sourceParameters.omega);

    /*! \property
      \page propertyPageACA ACA
      \section beta
      <code>MFloat AcaSolver::m_sourceParameters.beta</code>\n
      default = <code>sqrt(1 - Mach * Mach)</code>\n\n
      Beta factor for source parameters.\n
      Keywords: <i>ACA, input</i>
    */
    m_sourceParameters.beta = sqrt(1 - Mach * Mach);
    m_sourceParameters.beta = Context::getSolverProperty<MFloat>("beta", m_solverId, AT_, &m_sourceParameters.beta);

    // Amplitude
    m_sourceParameters.amplitude = (nDim == 2 ? 0.0001 : 0.01) * c0; // [m^2/s]

    {
      std::stringstream ss;
      ss << "Source parameters" << std::endl;
      ss << "   perturbed: " << m_sourceParameters.perturbed << std::endl;
      ss << "   omega: " << m_sourceParameters.omega << std::endl;
      ss << "   rho0: " << m_sourceParameters.rho0 << std::endl;
      ss << "   beta: " << m_sourceParameters.beta << std::endl;
      ss << "   amplitude: " << m_sourceParameters.amplitude << std::endl;
      ss << "   Ma: " << m_sourceParameters.Ma << std::endl;
      m_log << ss.str();
    }

    /*! \property
      \page propertyPageACA ACA
      \section generateSurfaceDataNoSamples
      <code>MInt AcaSolver::m_noSamples</code>\n
      default = <code>false</code>\n\n
      Set number of samples for generation of surface data.\n
      Keywords: <i>ACA, input</i>
    */
    m_noSamples = -1;
    m_noSamples = Context::getSolverProperty<MInt>("generateSurfaceDataNoSamples", m_solverId, AT_);
    m_totalNoSamples = m_noSamples;
    m_sampleStride = 1;
    m_sampleOffset = 0;

    // If uneven number of samples, take on less
    if(m_noSamples % 2 > 0) {
      m_noSamples = m_noSamples - 1;
      if(domainId() == 0) {
        std::cout << "Note: Number of Samples was uneven which could lead to mistakes. This is why the "
                     "new number of Samples is one less."
                  << std::endl;
      }
    }

    /*! \property
      \page propertyPageACA ACA
      \section noPeriods
      <code>MFloat noPeriods</code>\n
      default = <code>1.0</code>\n\n
      Read number of periods.\n
      Keywords: <i>ACA, input</i>
    */
    MFloat noPeriods = 1.0;
    noPeriods = Context::getSolverProperty<MFloat>("noPeriods", m_solverId, AT_, &noPeriods);
    TERMM_IF_COND(noPeriods < 1.0, "The number of periods is below 1!");

    /*! \property
      \page propertyPageACA ACA
      \section timeStepSize
      <code>MFloat AcaSolver::m_dt</code>\n
      default = <code>noPeriods * 2.0 * PI / (noSamples * m_sourceParameters.omega)</code>\n\n
      Computation of timeStep or read in if given.\n
      Keywords: <i>ACA, input</i>
    */
    m_dt = noPeriods * 2.0 * PI / (noSamples() * m_sourceParameters.omega);
    m_dt = Context::getSolverProperty<MFloat>("timeStepSize", m_solverId, AT_, &m_dt);
    const MFloat freq_max = m_sourceParameters.omega / (2 * PI);
    // Check Nyquist criterion afterwards
    TERMM_IF_COND(m_dt >= 1.0 / (2.0 * freq_max), "Nyquist criterions is not met!");

    // COUT the max frequency and the sampling rate to show that the Nyquist criterion is met
    if(domainId() == 0) {
      std::cout << " >> NOTE: The maximum frequency is f_max = omega / (2*pi) = " << freq_max
                << ". The Nyquist criterion is ensured ( m_dt < 1 /(2*freq_max) ), since..." << std::endl;
      std::cout << " >> m_dt = " << m_dt << " < " << 1.0 / (2.0 * freq_max) << " = 1 / (2*freq_max) <<" << std::endl;
    }
  } // End of "if m_generateSurfaceData"

  if(domainId() == 0) {
    std::cout << "Acoustic extrapolation method: " << methodName << std::endl;
    std::cout << "  * noVars = " << m_noVariables << std::endl;
    std::cout << "  * noComplexVars = " << m_noComplexVariables << std::endl;

    std::cout << "  * Input variables:" << std::endl;
    for(auto& var : m_inputVars) {
      std::cout << "    * " << var.first << ", position " << var.second << std::endl;
    }
  }
}


/// \brief Main compute method to perform acoustic far-field prediction
template <MInt nDim>
MBool AcaSolver<nDim>::solutionStep() {
  TRACE();
  RECORD_TIMER_START(m_timers[Timers::Run]);

  // Distribute surface elements among ranks and create separate communicators etc.
  initParallelization();

  initConversionFactors();

  // read in data of all surfaces
  // Note: for validating/testing instead of loading the data from disk you can generate the
  // analytical input data
  readSurfaceData();

  // read in observer points or generate the according to some properties if specified
  initObserverPoints();

  if(!m_acaPostprocessingMode) {
    // Compute source terms depending on method
    calcSourceTerms();

    // Compute observer pressure signal by using time or frequency FW-H formulation
    if(string2enum(solverMethod()) == MAIA_FWH_FREQUENCY) {
      // apply window to source terms and transform
      windowAndTransformSources();

      calcFrequencies();

      // Calculate the observer signals in frequency domain
      calculateObserverInFrequency();

      // backtransformation of local observer signal from frequency into time domain
      backtransformObservers();
    } else if(string2enum(solverMethod()) == MAIA_FWH_TIME) {
      // Surface Integration
      calculateObserverInTime();

      calcFrequencies();

      // Windowing and FFT of the observer data
      windowAndTransformObservers();
    }

    // change dimensions from ACA to desired output format
    changeDimensionsObserverData();

    // output of observer signals (frequency and time domain)
    saveObserverSignals();

    // perform additional post-processing operations with the observer signals, e.g., compute
    // RMS-pressure, average over observers, ..., and output the data
    postprocessObserverSignals();

    // Calculates the accuracy, the difference between analytic and computed solution. Can be only be
    // used with analytical data at the moment.
    // computeAccuracy();
  } else {
    loadObserverSignals();

    // Postprocessing with loaded observer data
    postprocessObserverSignals();
  }

  RECORD_TIMER_STOP(m_timers[Timers::Run]);

  m_isFinished = true;
  return true;
}


/// \brief Compute source terms for all surface elements depending on the selected method
template <MInt nDim>
void AcaSolver<nDim>::calcSourceTerms() {
  TRACE();
  RECORD_TIMER_START(m_timers[Timers::CalcSourceTerms]);
  printMessage("- Calculate source terms");

  // Note: totalNoSurfaceElements() returns the local, total number of elements for a surface
  // (may differ from domain to domain)

  if(string2enum(solverMethod()) == MAIA_FWH_FREQUENCY) {
    m_log << "Using Frequency-Domain FW-H Formulation" << std::endl;
    for(MInt i = 0; i < totalNoSurfaceElements(); i++) {
      calcSourceTermsFwhFreq(i);
    }
  } else if(string2enum(solverMethod()) == MAIA_FWH_TIME) {
    m_log << "Using Farassat 1A Time-Domain FW-H Formulation" << std::endl;
    for(MInt i = 0; i < totalNoSurfaceElements(); i++) {
      calcSourceTermsFwhTime(i);
    }
  } else {
    TERMM(1, "Unknown solverMethod");
  }

  RECORD_TIMER_STOP(m_timers[Timers::CalcSourceTerms]);
}

/// \brief Calculate FWH source terms for the frequency FWH formulation
/// Reference:
/// [1] Lockard, D.P., 2000. An efficient, two-dimensional ... doi.org/10.1006/jsvi.1999.2522
/// [2] Lockard, D.P., 2002, A comparison of Ffowcs ... doi.org/10.2514/6.2002-2580
template <MInt nDim>
void AcaSolver<nDim>::calcSourceTermsFwhFreq(const MInt segmentId) {
  TRACE();

  // Get pointer to normal vector
  const MFloat* const normal = &m_surfaceData.surfaceNormal(segmentId);
  // First two normal vector components (2D)
  const MFloat dfdx = normal[0];
  const MFloat dfdy = normal[1];
  // Mach numbers (2D)
  const MFloat Max = m_MaVec[0];
  const MFloat May = m_MaVec[1];

  const MFloat c_inf = 1.0;   // by definition
  const MFloat rho_inf = 1.0; // by definition
  const MFloat p_inf = rho_inf * (c_inf * c_inf) / m_gamma;
  const MFloat noSamplesInv = 1.0 / noSamples();

  if constexpr(nDim == 3) {
    // Third normal vector component for 3D
    const MFloat dfdz = normal[2];
    // Third Mach number for 3D
    const MFloat Maz = m_MaVec[2];

    // Compute mean velocities
    MFloat uMean = 0.0;
    MFloat vMean = 0.0;
    MFloat wMean = 0.0;
    if(m_acousticMethod == FWH_METHOD) {
      switch(m_fwhMeanStateType) {
        case 0: { // Default: use infnity properties to compute u'
          uMean = m_MaVec[0];
          vMean = m_MaVec[1];
          wMean = m_MaVec[2];
          break;
        }
        case 1: { // Use mean velocities to compute u'
          for(MInt t = 0; t < noSamples(); t++) {
            uMean += m_surfaceData.variables(segmentId, FWH::U, t);
            vMean += m_surfaceData.variables(segmentId, FWH::V, t);
            wMean += m_surfaceData.variables(segmentId, FWH::W, t);
          }
          uMean *= noSamplesInv;
          vMean *= noSamplesInv;
          wMean *= noSamplesInv;
          break;
        }
        default: {
          TERMM(1, "fwh mean state type not implemented");
          break;
        }
      }
    }

    // Time loop over all samples for 3D
    maia::parallelFor<false>(0, noSamples(), [=](MInt t) {
      // Get required variables (u', v', w', p, rho)
      std::vector<MInt> VarDim;
      MFloat up = 0.0;
      MFloat vp = 0.0;
      MFloat wp = 0.0;
      MFloat p = 0.0;
      MFloat rho = 0.0;
      if(m_acousticMethod == FWH_METHOD) {
        VarDim = {FWH::FX, FWH::FY, FWH::FZ, FWH::Q};
        up = m_surfaceData.variables(segmentId, FWH::U, t) - uMean;
        vp = m_surfaceData.variables(segmentId, FWH::V, t) - vMean;
        wp = m_surfaceData.variables(segmentId, FWH::W, t) - wMean;
        p = m_surfaceData.variables(segmentId, FWH::P, t);
        rho = m_surfaceData.variables(segmentId, FWH::RHO, t);
      } else if(m_acousticMethod == FWH_APE_METHOD) {
        VarDim = {FWH_APE::FX, FWH_APE::FY, FWH_APE::FZ, FWH_APE::Q};
        up = m_surfaceData.variables(segmentId, FWH_APE::UP, t);
        vp = m_surfaceData.variables(segmentId, FWH_APE::VP, t);
        wp = m_surfaceData.variables(segmentId, FWH_APE::WP, t);
        p = m_surfaceData.variables(segmentId, FWH_APE::PP, t) + p_inf;
        rho = m_surfaceData.variables(segmentId, FWH_APE::PP, t) + rho_inf;
      }

      const MFloat f1 = (up + Max) * dfdx + (vp + May) * dfdy + (wp + Maz) * dfdz;
      const MFloat f2 = Max * dfdx + May * dfdy + Maz * dfdz;

      // Calculation of source terms (3D) Fx, Fy, Fz and Q for each sample.
      m_surfaceData.variables(segmentId, VarDim[0], t) = p * dfdx + rho * (up - Max) * f1 + Max * f2;
      m_surfaceData.variables(segmentId, VarDim[1], t) = p * dfdy + rho * (vp - May) * f1 + May * f2;
      m_surfaceData.variables(segmentId, VarDim[2], t) = p * dfdz + rho * (wp - Maz) * f1 + Maz * f2;
      m_surfaceData.variables(segmentId, VarDim[3], t) = rho * f1 - f2;
    });
  } else {
    // Compute mean velocities
    MFloat uMean = 0.0;
    MFloat vMean = 0.0;
    if(m_acousticMethod == FWH_METHOD) {
      switch(m_fwhMeanStateType) {
        case 0: { // Default: use infnity properties to compute u'
          uMean = m_MaVec[0];
          vMean = m_MaVec[1];
          break;
        }
        case 1: { // Use mean velocities to compute u'
          for(MInt t = 0; t < noSamples(); t++) {
            uMean += m_surfaceData.variables(segmentId, FWH::U, t);
            vMean += m_surfaceData.variables(segmentId, FWH::V, t);
          }
          uMean = uMean * noSamplesInv;
          vMean = vMean * noSamplesInv;
          break;
        }
        default: {
          TERMM(1, "fwh mean state type not implemented");
          break;
        }
      }
    }

    // Time loop over all samples for 2D
    maia::parallelFor<false>(0, noSamples(), [=](MInt t) {
      // Get required variables (u', v', p, rho)
      std::vector<MInt> VarDim;
      MFloat up = 0.0;
      MFloat vp = 0.0;
      MFloat p = 0.0;
      MFloat rho = 0.0;
      if(m_acousticMethod == FWH_METHOD) {
        VarDim = {FWH::FX, FWH::FY, FWH::Q};
        up = m_surfaceData.variables(segmentId, FWH::U, t) - uMean;
        vp = m_surfaceData.variables(segmentId, FWH::V, t) - vMean;
        p = m_surfaceData.variables(segmentId, FWH::P, t);
        rho = m_surfaceData.variables(segmentId, FWH::RHO, t);
      } else if(m_acousticMethod == FWH_APE_METHOD) {
        VarDim = {FWH_APE::FX, FWH_APE::FY, FWH_APE::Q};
        up = m_surfaceData.variables(segmentId, FWH_APE::UP, t);
        vp = m_surfaceData.variables(segmentId, FWH_APE::VP, t);
        p = m_surfaceData.variables(segmentId, FWH_APE::PP, t) + p_inf;
        rho = m_surfaceData.variables(segmentId, FWH_APE::PP, t) + rho_inf;
      }

      const MFloat f1 = (up + Max) * dfdx + (vp + May) * dfdy;
      const MFloat f2 = Max * dfdx + May * dfdy;

      // Calculation of source terms (2D) Fx, Fy and Q for each sample.
      m_surfaceData.variables(segmentId, VarDim[0], t) = p * dfdx + rho * (up - Max) * f1 + Max * f2;
      m_surfaceData.variables(segmentId, VarDim[1], t) = p * dfdy + rho * (vp - May) * f1 + May * f2;
      m_surfaceData.variables(segmentId, VarDim[2], t) = rho * f1 - f2;
    });
  }
}

/// \brief Calculate FWH source terms for the time domain FWH formulation
/// Reference:
/// [1] Brentner, K. S. and F. Farassat, 1998, Analytical Comparison ... doi.org/10.2514/2.558
/// [2] Casalino, D., 2003, An advanced time approach ... doi.org/10.1016/s0022-460x(02)00986-0
/// [3] Brès, G., 2010, A Ffowcs Williams - Hawkings ... doi.org/10.2514/6.2010-3711
template <MInt nDim>
void AcaSolver<nDim>::calcSourceTermsFwhTime(const MInt segmentId) {
  TRACE();

  // Get pointer to normal vector
  const MFloat* const normal = &m_surfaceData.surfaceNormal(segmentId);
  const MFloat c_inf = 1.0;   // by definition
  const MFloat rho_inf = 1.0; // by definition
  const MFloat p_inf = rho_inf * (c_inf * c_inf) / m_gamma;
  const MFloat noSamplesInv = 1.0 / noSamples();

  if constexpr(nDim == 3) {
    // Compute mean variables
    MFloat uMean = 0.0;
    MFloat vMean = 0.0;
    MFloat wMean = 0.0;
    MFloat pMean = 0.0;
    if(m_acousticMethod == FWH_METHOD) {
      switch(m_fwhMeanStateType) {
        case 0: { // Default: use infnity properties to compute u'
          uMean = m_MaVec[0];
          vMean = m_MaVec[1];
          wMean = m_MaVec[2];
          pMean = p_inf;
          break;
        }
        case 1: { // Use mean velocities to compute u'
          for(MInt t = 0; t < noSamples(); t++) {
            uMean += m_surfaceData.variables(segmentId, FWH::U, t);
            vMean += m_surfaceData.variables(segmentId, FWH::V, t);
            wMean += m_surfaceData.variables(segmentId, FWH::W, t);
            pMean += m_surfaceData.variables(segmentId, FWH::P, t);
          }
          uMean *= noSamplesInv;
          vMean *= noSamplesInv;
          wMean *= noSamplesInv;
          pMean *= noSamplesInv;
          break;
        }
        default: {
          TERMM(1, "fwh mean state type not implemented");
          break;
        }
      }
    }

    // Set surface geometry properties
    // TODO: This assumes a stationary sampling surface! Impliment the version for moving surface...
    MFloat nx = normal[0];
    MFloat ny = normal[1];
    MFloat nz = normal[2];
    MFloat surfVX = -m_MaVec[0];
    MFloat surfVY = -m_MaVec[1];
    MFloat surfVZ = -m_MaVec[2];

    maia::parallelFor<false>(0, noSamples(), [=](MInt tau) {
      // Get surface (f=0) velocity --> based on the ABSOLUTE frame of reference
      m_surfaceData.variables(segmentId, FWH_TIME::surfMX, tau) = surfVX;
      m_surfaceData.variables(segmentId, FWH_TIME::surfMY, tau) = surfVY;
      m_surfaceData.variables(segmentId, FWH_TIME::surfMZ, tau) = surfVZ;

      // Get primative variables --> based on the ABSOLUTE frame of reference
      MFloat up = 0.0;
      MFloat vp = 0.0;
      MFloat wp = 0.0;
      MFloat pp = 0.0;
      MFloat rho = 0.0;
      if(m_acousticMethod == FWH_METHOD) {
        up = m_surfaceData.variables(segmentId, FWH::U, tau) - uMean;
        vp = m_surfaceData.variables(segmentId, FWH::V, tau) - vMean;
        wp = m_surfaceData.variables(segmentId, FWH::W, tau) - wMean;
        pp = m_surfaceData.variables(segmentId, FWH::P, tau) - pMean;
        rho = m_surfaceData.variables(segmentId, FWH::RHO, tau);
      } else if(m_acousticMethod == FWH_APE_METHOD) {
        up = m_surfaceData.variables(segmentId, FWH_APE::UP, tau);
        vp = m_surfaceData.variables(segmentId, FWH_APE::VP, tau);
        wp = m_surfaceData.variables(segmentId, FWH_APE::WP, tau);
        pp = m_surfaceData.variables(segmentId, FWH_APE::PP, tau);
        rho = rho_inf + pp; // = rho_inf + p'/(c_inf * c_inf)
      }

      // Compute sources
      const MFloat un = up * nx + vp * ny + wp * nz;
      const MFloat vn = surfVX * nx + surfVY * ny + surfVZ * nz;
      const MFloat deltaUn = un - vn;

      m_surfaceData.variables(segmentId, FWH_TIME::srcUn, tau) = vn + rho / rho_inf * deltaUn;
      m_surfaceData.variables(segmentId, FWH_TIME::srcLX, tau) = pp * nx + rho * up * deltaUn;
      m_surfaceData.variables(segmentId, FWH_TIME::srcLY, tau) = pp * ny + rho * vp * deltaUn;
      m_surfaceData.variables(segmentId, FWH_TIME::srcLZ, tau) = pp * nz + rho * wp * deltaUn;

      // Transform the source terms to the flow direction coordinate
      transformCoordinatesToFlowDirection3D(&m_surfaceData.variables(segmentId, FWH_TIME::surfMX, tau),
                                            &m_surfaceData.variables(segmentId, FWH_TIME::surfMY, tau),
                                            &m_surfaceData.variables(segmentId, FWH_TIME::surfMZ, tau));
      transformCoordinatesToFlowDirection3D(&m_surfaceData.variables(segmentId, FWH_TIME::srcLX, tau),
                                            &m_surfaceData.variables(segmentId, FWH_TIME::srcLY, tau),
                                            &m_surfaceData.variables(segmentId, FWH_TIME::srcLZ, tau));

      /*// Debugging
      MFloat* surfCoordinates = &m_surfaceData.surfaceCoordinates(segmentId);
      if(segmentId == 0 && tau == 0) {
        m_log << "tau = " << tau
              << ", surfCoordinates[0] = " << surfCoordinates[0]
              << ", surfCoordinates[1] = " << surfCoordinates[1]
              << ", surfCoordinates[2] = " << surfCoordinates[2]
              << ", nx = " << normal[0]
              << ", ny = " << normal[1]
              << ", nz = " << normal[2]
              << ", surfVX = " << surfVX
              << ", surfVY = " << surfVY
              << ", surfVZ = " << surfVZ
              << ", up = " << up
              << ", vp = " << vp
              << ", wp = " << wp
              << ", rho = " << rho
              << ", pp = " << pp
              << ", un = " << un
              << ", vn = " << vn
              << ", deltaUn = " << deltaUn
              << ", srcUn = " << std::setprecision(16) << m_surfaceData.variables(segmentId, FWH_TIME::srcUn, tau)
              << ", srcLX = " << std::setprecision(16) << m_surfaceData.variables(segmentId, FWH_TIME::srcLX, tau)
              << ", srcLY = " << std::setprecision(16) << m_surfaceData.variables(segmentId, FWH_TIME::srcLY, tau)
              << ", srcLZ = " << std::setprecision(16) << m_surfaceData.variables(segmentId, FWH_TIME::srcLZ, tau)
              << std::endl;
      }*/
    });
  } else if constexpr(nDim == 2) {
    TERMM(1, "Time Domain FW-H not implimented in 2D: Requires free space Green function in 2D");
  } // End of dimension switch
}

/// \brief Compute the window function factors and the normalization factor
template <MInt nDim>
void AcaSolver<nDim>::calcWindowFactors(MFloat* const p_window) {
  const MFloat N = static_cast<MFloat>(noSamples());
  MFloat normFactor = 1.0;

  switch(m_windowType) {
    case WINDOW::NONE: {
      m_log << " - using no window function" << std::endl;
      break;
    }
    case WINDOW::HANNING: {
      m_log << " - using Hanning window function" << std::endl;
      MFloat sum_w_sq = 0.0;
      for(MInt t = 0; t < noSamples(); t++) {
        const MFloat w = 0.5 * (1.0 - cos(2.0 * PI * t / N));
        sum_w_sq += w * w;
        p_window[t] = w;
      }
      normFactor = 1.0 / sqrt(sum_w_sq / N);
      break;
    }
    case WINDOW::HAMMING: {
      m_log << " - using Hamming window function" << std::endl;
      MFloat sum_w_sq = 0.0;
      for(MInt t = 0; t < noSamples(); t++) {
        const MFloat w = 0.54 - 0.46 * cos(2.0 * PI * t / N);
        sum_w_sq += w * w;
        p_window[t] = w;
      }
      normFactor = 1.0 / sqrt(sum_w_sq / N);
      break;
    }
    case WINDOW::MODHANNING: {
      m_log << " - using modified-Hanning window function" << std::endl;
      MFloat sum_w_sq = 0.0;
      for(MInt t = 0; t < noSamples(); t++) {
        MFloat w = 0.0;
        if(t < N / 8.0 || t > 7.0 * N / 8.0) {
          w = 0.5 * (1.0 - cos(8.0 * PI * t / N));
        } else {
          w = 1.0;
        }
        sum_w_sq += w * w;
        p_window[t] = w;
      }
      normFactor = 1.0 / sqrt(sum_w_sq / N);
      break;
    }
    default: {
      TERMM(1, "Invalid window type: '" + std::to_string(m_windowType) + "'");
      break;
    }
  }
  m_log << " - normalization factor: " << normFactor << std::endl;
  m_windowNormFactor = normFactor;
}

/// \brief Apply window function to source terms and transform into frequency domain
template <MInt nDim>
void AcaSolver<nDim>::windowAndTransformSources() {
  TRACE();
  RECORD_TIMER_START(m_timers[Timers::WindowAndTransformSources]);
  printMessage("- Window and transform source terms");

  // Compute the window function factors and the normalization factor
  std::vector<MFloat> window(noSamples());
  std::fill(window.begin(), window.end(), 1.0);
  calcWindowFactors(window.data());

  // Set source term indices depending on used method
  std::vector<MInt> VarDim;
  std::vector<MInt> VarDimC;
  if(m_acousticMethod == FWH_METHOD) {
    if constexpr(nDim == 3) {
      VarDim = {FWH::FX, FWH::FY, FWH::FZ, FWH::Q};
      VarDimC = {FWH::FX_C, FWH::FY_C, FWH::FZ_C, FWH::Q_C};
    } else {
      VarDim = {FWH::FX, FWH::FY, FWH::Q};
      VarDimC = {FWH::FX_C, FWH::FY_C, FWH::Q_C};
    }
  } else if(m_acousticMethod == FWH_APE_METHOD) {
    if constexpr(nDim == 3) {
      VarDim = {FWH_APE::FX, FWH_APE::FY, FWH_APE::FZ, FWH_APE::Q};
      VarDimC = {FWH_APE::FX_C, FWH_APE::FY_C, FWH_APE::FZ_C, FWH_APE::Q_C};
    } else {
      VarDim = {FWH_APE::FX, FWH_APE::FY, FWH_APE::Q};
      VarDimC = {FWH_APE::FX_C, FWH_APE::FY_C, FWH_APE::Q_C};
    }
  } else {
    TERMM(1, "Unknown acoustic method.");
  }

  MInt noSourceTerms = 0;
  if constexpr(nDim == 3) {
    noSourceTerms = 4;
  } else if constexpr(nDim == 2) {
    noSourceTerms = 3;
  }

  const MFloat noSamplesF = static_cast<MFloat>(noSamples());
  for(MInt segmentId = 0; segmentId < totalNoSurfaceElements(); segmentId++) {
    for(MInt i = 0; i < noSourceTerms; i++) {
      // Sum up source term over all samples
      MFloat sumSource = 0.0;
      for(MInt t = 0; t < noSamples(); t++) {
        sumSource += m_surfaceData.variables(segmentId, VarDim[i], t);
      }
      const MFloat avgSource = sumSource / noSamplesF; // Compute average

      // Substract average from all samples to get the fluctuations
      // After substracting the mean apply the window function [see e.g. Mendez et al. or Lockard]
      for(MInt t = 0; t < noSamples(); t++) {
        m_surfaceData.variables(segmentId, VarDim[i], t) -= avgSource;
        m_surfaceData.variables(segmentId, VarDim[i], t) *= window[t];
      }
    }

    // Copy variables into complex variables array: Re1, Im1, Re2, Im2, Re3, Im3, ... with Im_n = 0
    for(MInt i = 0; i < noSourceTerms; i++) {
      for(MInt t = 0; t < noSamples(); t++) {
        m_surfaceData.complexVariables(segmentId, VarDimC[i], t, 0) =
            m_surfaceData.variables(segmentId, VarDim[i], t);              // Real part
        m_surfaceData.complexVariables(segmentId, VarDimC[i], t, 1) = 0.0; // Imaginary part
      }
    }

    // Use FFT or DFT to transform source terms (depending on possibility and the set property)
    checkNoSamplesPotencyOfTwo();
    if(m_FastFourier == true && m_transformationType == 1) {
      for(MInt i = 0; i < noSourceTerms; i++) {
        FastFourierTransform(&m_surfaceData.complexVariables(segmentId, VarDimC[i]), noSamples(), 1, nullptr, true);
      }
    } else if(m_FastFourier == false || (m_FastFourier == true && m_transformationType == 0)) {
      for(MInt i = 0; i < noSourceTerms; i++) {
        // 1 in function means forward; -1 is backward
        DiscreteFourierTransform(&m_surfaceData.complexVariables(segmentId, VarDimC[i]), noSamples(), 1, nullptr, true);
      }
    } else {
      TERMM(1, "Error Fourier-Transform: Check transformationType and noSamples.");
    }

    // Normalize for conservation of energy
    for(MInt t = 0; t < noSamples(); t++) {
      for(MInt i = 0; i < noSourceTerms; i++) {
        m_surfaceData.complexVariables(segmentId, VarDimC[i], t, 0) *= m_windowNormFactor;
        m_surfaceData.complexVariables(segmentId, VarDimC[i], t, 1) *= m_windowNormFactor;
      }
    }

    // Calculate frequencies is moved out of segmentId Loop because it is the same for every
    // segmentId. Moved to calcSurfaceIntegral Calculate time: It is calculated at
    // "generateSurfaceData" m_times = dt * noSamples() Initialize observerdata array. Function
    // initObserverPoints
  } // End of loop over all segments

  RECORD_TIMER_STOP(m_timers[Timers::WindowAndTransformSources]);
}

/// \brief Apply window function to the observer time series data and transform it into frequency domain
template <MInt nDim>
void AcaSolver<nDim>::windowAndTransformObservers() {
  TRACE();

  // Compute the window function factors and the normalization factor
  std::vector<MFloat> window(noSamples());
  std::fill(window.begin(), window.end(), 1.0);
  calcWindowFactors(window.data());

  const MFloat noSamplesF = static_cast<MFloat>(noSamples());
  for(MInt observerId = 0; observerId < noObservers(); observerId++) {
    // Sum up observer pressure over all samples
    MFloat sumObsPressure = 0.0;
    for(MInt t = 0; t < noSamples(); t++) {
      sumObsPressure += m_observerData.variables(observerId, 0, t);
    }
    const MFloat aveObsPressure = sumObsPressure / noSamplesF; // Compute average

    // Substract average from all samples to get the fluctuations
    // After substracting the mean apply the window function [see e.g. Mendez et al. or Lockard]
    // copy variables into complex variables array: re1, im1, re2, im2, re3, im3, ... with im_n = 0
    for(MInt t = 0; t < noSamples(); t++) {
      // Do nothing to "m_observerData.variables" output the original data
      m_observerData.complexVariables(observerId, 0, t, 0) =
          window[t] * (m_observerData.variables(observerId, 0, t) - aveObsPressure); // Real part
      m_observerData.complexVariables(observerId, 0, t, 1) = 0.0;                    // Imaginary part
    }

    // Use FFT or DFT to transform source terms (depending on possibility and the set property)
    checkNoSamplesPotencyOfTwo();
    if(m_FastFourier == true && m_transformationType == 1) {
      FastFourierTransform(&m_observerData.complexVariables(observerId, 0), noSamples(), 1, nullptr, true);
    } else if(m_FastFourier == false || (m_FastFourier == true && m_transformationType == 0)) {
      // 1 in function means forward; -1 is backward
      DiscreteFourierTransform(&m_observerData.complexVariables(observerId, 0), noSamples(), 1, nullptr, true);
    } else {
      TERMM(1, "Error Fourier-Transform: Check transformationType and noSamples.");
    }

    // Normalize for conservation of energy
    for(MInt t = 0; t < noSamples(); t++) {
      m_observerData.complexVariables(observerId, 0, t, 0) *= m_windowNormFactor;
      m_observerData.complexVariables(observerId, 0, t, 1) *= m_windowNormFactor;
    }
  } // End of observer loop
}


/// \brief Fast Fourier transformation.
///
/// Reference: 'Numerical Recipes in C', Cambridge, p. 507 ff.
///
/// \param dataArray Pointer to complex data array.
/// \param size Size of the complex data array (number of samples).
/// \param dir Direction of transformation: 1=forward; -1=backward.
/// \param result Optional storage for transformed data, used if inPlace==false.
/// \param inPlace Determines if result overwrites the input data or stored in result array.
template <MInt nDim>
void AcaSolver<nDim>::FastFourierTransform(MFloat* dataArray, const MInt size, const MInt dir, MFloat* result,
                                           const MBool inPlace) {
  TRACE();
  TERMM_IF_NOT_COND(dir == 1 || dir == -1, "Invalid direction");
  TERMM_IF_NOT_COND(size % 2 == 0, "FFT: not a multiple of 2 samples.");

  const MFloat isign = (dir == 1) ? -1.0 : 1.0;
  const MBool dim = (dir == 1) ? true : false;
  const MInt size2 = 2 * size; // full array size (real+imag)

  // Set pointer to result storage
  MFloat* const data = (inPlace) ? dataArray : result;
  // Copy input data if not performing transform in-place (overwriting input)
  if(!inPlace) {
    ASSERT(result != nullptr, "Pointer to result storage is nullptr.");
    std::copy_n(&dataArray[0], size2, result);
  }

#define SWAP(a, b)                                                                                                     \
  rtemp = (a);                                                                                                         \
  (a) = (b);                                                                                                           \
  (b) = rtemp
  // Bit reversing: Reorder by inverse bit notation, e.g. first number: 0000 --> 0000, second
  // number: 0001 -> 1000
  const MInt nn = size;
  MInt n = nn << 1;
  MInt m;
  MFloat rtemp;
  MInt j = 1;
  for(MInt i = 1; i < n; i += 2) {
    if(j > i) {
      SWAP(data[j - 1], data[i - 1]);
      SWAP(data[j], data[i]);
    }
    m = nn;
    while(m >= 2 && j > m) {
      j -= m;
      m >>= 1;
    }
    j += m;
  }

  // Fourier Transformation (Danielson-Lanczos)
  MUlong istep, p, k, i;
  MFloat wtemp, wr, wpr, wpi, wi, theta, tempr, tempi;
  const MUlong nLong = static_cast<MUlong>(size2);
  MUlong mmax = 2;
  while(nLong > mmax) {
    istep = mmax << 1;
    theta = isign * (6.28318530717959 / mmax);
    wtemp = sin(0.5 * theta);
    wpr = -2.0 * wtemp * wtemp;
    wpi = sin(theta);
    wr = 1.0;
    wi = 0.0;
    for(p = 1; p < mmax; p += 2) {
      for(i = p; i <= nLong; i += istep) {
        k = i + mmax;
        tempr = wr * data[k - 1] - wi * data[k];
        tempi = wr * data[k] + wi * data[k - 1];
        data[k - 1] = data[i - 1] - tempr;
        data[k] = data[i] - tempi;
        data[i - 1] += tempr;
        data[i] += tempi;
      }
      wr = (wtemp = wr) * wpr - wi * wpi + wr;
      wi = wi * wpr + wtemp * wpi + wi;
    }
    mmax = istep;
  }

  // Dimensionierung
  if(dim) {
    for(MInt z = 0; z < size2; z++) {
      data[z] /= size;
    }
  }
}


/// \brief Discrete Fourier transform (computational expensive O(n^2) instead of O(n*logn) compared
///        to FFT)
///
/// \param dataArray Pointer to complex data array.
/// \param size Size of the complex data array (number of samples).
/// \param dir Direction of transformation: 1=forward; -1=backward.
/// \param result Optional storage for transformed data, used if inPlace==false.
/// \param inPlace Determines if result overwrites the input data or stored in result array.
template <MInt nDim>
void AcaSolver<nDim>::DiscreteFourierTransform(MFloat* dataArray, const MInt size, const MInt dir, MFloat* result,
                                               const MBool inPlace) {
  TRACE();
  const MFloat sizeF = static_cast<MFloat>(size);
  const MInt size2 = 2 * size;
  MFloat* resultPtr = nullptr;
  std::vector<MFloat> resultArray;

  if(inPlace) {
    resultArray.resize(size2); // Temporary solution storage
    resultPtr = &resultArray[0];
  } else {
    ASSERT(result != nullptr, "Pointer to result storage is nullptr.");
    resultPtr = result; // Store solution in provided storage
  }

  if(dir == 1) { // FOURIER FORWARD
    const MFloat con = 2.0 * PI / sizeF;
    for(MInt k = 0; k < size; k++) {
      for(MInt i = 0; i < size; i++) {
        // cos(-a) = cos(a) and sin(-a) = -sin(a)
        // Actual Discrete Fourier Transform
        resultPtr[2 * k] += dataArray[2 * i] * cos(k * con * i) + dataArray[2 * i + 1] * sin(k * con * i);
        resultPtr[2 * k + 1] += -1.0 * dataArray[2 * i] * sin(k * con * i) + dataArray[2 * i + 1] * cos(k * con * i);
      }
      resultPtr[2 * k] /= sizeF;
      resultPtr[2 * k + 1] /= sizeF;
    }
  } else if(dir == -1) { // FOURIER BACKWARD
    for(MInt t = 0; t < size; t++) {
      const MFloat arg = 2.0 * PI * t / sizeF;
      for(MInt k = 0; k < size; k++) {
        resultPtr[2 * t] += dataArray[2 * k] * cos(arg * k) - dataArray[2 * k + 1] * sin(arg * k);
        resultPtr[2 * t + 1] += dataArray[2 * k] * sin(arg * k) + dataArray[2 * k + 1] * cos(arg * k);
      }
    }
  } else {
    TERMM(1, "Please insert 1 for forward or -1 for backward");
  }

  if(inPlace) {
    // Copy results back into dataArray overwriting the input data
    std::copy_n(&resultPtr[0], size2, &dataArray[0]);
  }
}


/// \brief Calculate the frequencies
template <MInt nDim>
void AcaSolver<nDim>::calcFrequencies() {
  TRACE();
  printMessage("- Calculate frequencies");
  // Calculate frequencies, m_frequencies is already resized with noSamples
  // Numerical Recipes in C by Cambridge page 503
  // Structure: 0 +1 +2 +3 ... +noSamples/2 ... -3 -2 -1
  for(MInt i = 0, k = noSamples() - 1; i <= noSamples() / 2; k--, i++) {
    if(k > noSamples() / 2) {
      m_frequencies[k] = -1.0 * static_cast<MFloat>(i + 1) / (noSamples() * m_dt);
    }
    m_frequencies[i] = static_cast<MFloat>(i) / (noSamples() * m_dt);
  }
  // Frequencies from Hz --> rad/s: now omega
  for(MInt i = 0; i < noSamples(); i++) {
    m_frequencies[i] *= 2.0 * PI;
  }
}

/// \brief  Calculate the FWH surface integrals for a certain observer location
/// \author Miro Gondrum
/// \date   15.06.2023
/// \param[in]  coord               Location of the observer
/// \param[out] p_complexVariables  Output of complex p', length of 2*noSamples()
template <MInt nDim>
void AcaSolver<nDim>::calcSurfaceIntegralsForObserver(const std::array<MFloat, nDim> coord,
                                                      MFloat* const p_complexVariables) {
  RECORD_TIMER_START(m_timers[Timers::CalcSurfaceIntegrals]);
  std::vector<MInt> VarDimC;
  if(m_acousticMethod == FWH_METHOD) {
    if constexpr(nDim == 3) {
      VarDimC = {FWH::FX_C, FWH::FY_C, FWH::FZ_C, FWH::Q_C};
    } else {
      VarDimC = {FWH::FX_C, FWH::FY_C, FWH::Q_C};
    }
  } else if(m_acousticMethod == FWH_APE_METHOD) {
    if constexpr(nDim == 3) {
      VarDimC = {FWH_APE::FX_C, FWH_APE::FY_C, FWH_APE::FZ_C, FWH_APE::Q_C};
    } else {
      VarDimC = {FWH_APE::FX_C, FWH_APE::FY_C, FWH_APE::Q_C};
    }
  } else {
    TERMM(1, "Unknown acoustic method.");
  }

  // Calculate Mach number
  const MFloat mach = sqrt(std::inner_product(&m_MaVec[0], &m_MaVec[nDim], &m_MaVec[0], 0.0));
  // Prandtl-Glauert-Factor
  const MFloat betaSq = 1.0 - mach * mach;
  const MFloat beta = sqrt(betaSq);
  const MFloat fbeta2 = 1 / betaSq;

  if constexpr(nDim == 2) {
    const MFloat Max = m_MaVec[0];
    const MFloat May = m_MaVec[1];

    MFloat theta = 0;
    if(!approx(Max, 0.0, MFloatEps) && !approx(May, 0.0, MFloatEps)) {
      theta = atan2(May, Max);
    }

    // Pre-computing sine and cosine terms
    const MFloat sinTheta = sin(theta);
    const MFloat cosTheta = cos(theta);

    // INNER LOOP OVER ALL FREQUENCIES (with only positive frequencies because bessel is not valid
    // for negativ nor zero). Since nw = 0 is always a 0, it is left out in the loop --> start at
    // nw = 1
    maia::parallelFor<false>(1, noSamples() / 2 + 1, [=](MInt nw) {
      // non-dimensionalized omega is already in m_frequencies
      const MFloat waveNumber = m_frequencies[nw];

      // Set everything to 0
      MFloat integralFxRe = 0.0;
      MFloat integralFxIm = 0.0;
      MFloat integralFyRe = 0.0;
      MFloat integralFyIm = 0.0;
      MFloat integralQRe = 0.0;
      MFloat integralQIm = 0.0;

      // LOOP OVER ALL LOCAL SURFACE ELEMENTS
      for(MInt segmentId = 0; segmentId < totalNoSurfaceElements(); segmentId++) {
        // Calculation of Prandtl-Glauert-Coordinates
        // Division by reference length is leaved out bc it would cancel out later
        const MFloat* const surfCoordinates = &m_surfaceData.surfaceCoordinates(segmentId);
        const MFloat surfX = surfCoordinates[0];
        const MFloat surfY = surfCoordinates[1];
        const MFloat X = (coord[0] - surfX) * cosTheta + (coord[1] - surfY) * sinTheta;
        const MFloat Y = -1.0 * (coord[0] - surfX) * sinTheta + (coord[1] - surfY) * cosTheta;

        const MFloat d = sqrt(X * X + betaSq * Y * Y);
        const MFloat fd = 1 / d;

        // Argument in Hankelfunction argH and its derivatives
        const MFloat argH = waveNumber * fbeta2 * d;
        const MFloat dargHdxi = waveNumber * fbeta2 * fd * (-1.0 * X * cosTheta + betaSq * Y * sinTheta);
        const MFloat dargHdeta = waveNumber * fbeta2 * fd * (-1.0 * X * sinTheta - betaSq * Y * cosTheta);

        // Argument for exp function and its derivaties
        const MFloat argExp = mach * waveNumber * X * fbeta2;
        const MFloat dargEdxi = -1.0 * mach * waveNumber * cosTheta * fbeta2;
        const MFloat dargEdeta = -1.0 * mach * waveNumber * sinTheta * fbeta2;
        const MFloat C1 = 1.0 / (4.0 * beta);

        const MFloat J0 = j0(argH);
        const MFloat J1 = j1(argH);
        const MFloat Y0 = y0(argH);
        const MFloat Y1 = y1(argH);

        const MFloat cosargE = cos(argExp);
        const MFloat sinargE = sin(argExp);

        // Calculation of the Green's function. y0, j0, y1 and j1 are used for Hankel function
        const MFloat greenValue_Re = C1 * (cosargE * Y0 - sinargE * J0);
        const MFloat greenValue_Im = C1 * (cosargE * J0 + sinargE * Y0);

        // Calculation of the derivation of the Greenfunction resp. xi and eta
        const MFloat dGreendXi_Re =
            C1 * (dargEdxi * (-1.0 * cosargE * J0 - sinargE * Y0) + dargHdxi * (-1.0 * cosargE * Y1 + sinargE * J1));
        const MFloat dGreendXi_Im =
            C1 * (dargEdxi * (cosargE * Y0 - sinargE * J0) + dargHdxi * (-1.0 * cosargE * J1 - sinargE * Y1));
        const MFloat dGreendEta_Re =
            C1 * (dargEdeta * (-1.0 * cosargE * J0 - sinargE * Y0) + dargHdeta * (-1.0 * cosargE * Y1 + sinargE * J1));
        const MFloat dGreendEta_Im =
            C1 * (dargEdeta * (cosargE * Y0 - sinargE * J0) + dargHdeta * (-1.0 * cosargE * J1 - sinargE * Y1));

        // Get Grid size dL
        const MFloat dL = m_surfaceData.surfaceArea(segmentId);
        // Get local values
        const MFloat FxRe = m_surfaceData.complexVariables(segmentId, VarDimC[0], nw, 0);
        const MFloat FxIm = m_surfaceData.complexVariables(segmentId, VarDimC[0], nw, 1);

        const MFloat FyRe = m_surfaceData.complexVariables(segmentId, VarDimC[1], nw, 0);
        const MFloat FyIm = m_surfaceData.complexVariables(segmentId, VarDimC[1], nw, 1);

        const MFloat QRe = m_surfaceData.complexVariables(segmentId, VarDimC[2], nw, 0);
        const MFloat QIm = m_surfaceData.complexVariables(segmentId, VarDimC[2], nw, 1);

        // Calculate integral terms
        integralFxRe += dL * (FxRe * dGreendXi_Re - FxIm * dGreendXi_Im);
        integralFxIm += dL * (FxRe * dGreendXi_Im + FxIm * dGreendXi_Re);

        integralFyRe += dL * (FyRe * dGreendEta_Re - FyIm * dGreendEta_Im);
        integralFyIm += dL * (FyRe * dGreendEta_Im + FyIm * dGreendEta_Re);

        integralQRe += -1.0 * waveNumber * dL * (QRe * greenValue_Im + QIm * greenValue_Re);
        integralQIm += waveNumber * dL * (QRe * greenValue_Re - QIm * greenValue_Im);
      } // end segmentId Loop

      // sum up all integral terms
      const MFloat integralRe = -1.0 * (integralFxRe + integralFyRe + integralQRe);
      const MFloat integralIm = -1.0 * (integralFxIm + integralFyIm + integralQIm);

      // Store pressure values in observerData
      p_complexVariables[2 * nw + 0] = integralRe;
      p_complexVariables[2 * nw + 1] = integralIm;
    });
  } else if constexpr(nDim == 3) {
    const MFloat Max = m_MaVec[0];
    const MFloat May = m_MaVec[1];
    const MFloat Maz = m_MaVec[2];

    MFloat theta = 0.0;
    MFloat alpha = 0.0;
    // Determination of the flow direction
    // Use approx(..., 0.0, MFloatEps) instead of Maz == 0.0 || May == 0.0. Compares Ma with a
    // really low number epsilon
    if(!approx(Maz, 0.0, MFloatEps) && !approx(Max, 0.0, MFloatEps)) {
      theta = atan2(Maz, Max);
    }
    if(!approx(May, 0.0, MFloatEps)) {
      alpha = asin(May / mach);
    }

    // Precompute trigonometric values
    const MFloat sinAlpha = sin(alpha);
    const MFloat cosAlpha = cos(alpha);
    const MFloat sinTheta = sin(theta);
    const MFloat cosTheta = cos(theta);

    // FREQUENCIE LOOP: only for positive frequencies without nw=0
    maia::parallelFor<false>(1, noSamples() / 2 + 1, [=](MInt nw) {
      // Calculate wave number: Non-dimensionalized omega is already saved in m_frequencies
      const MFloat waveNumber = m_frequencies[nw];

      // Set values to zero
      MFloat integralFxRe = 0.0;
      MFloat integralFxIm = 0.0;
      MFloat integralFyRe = 0.0;
      MFloat integralFyIm = 0.0;
      MFloat integralFzRe = 0.0;
      MFloat integralFzIm = 0.0;
      MFloat integralQRe = 0.0;
      MFloat integralQIm = 0.0;

      // SURFACE LOOP
      for(MInt segmentId = 0; segmentId < totalNoSurfaceElements(); segmentId++) {
        MFloat* surfCoordinates = &m_surfaceData.surfaceCoordinates(segmentId);

        const MFloat x = coord[0];
        const MFloat y = coord[1];
        const MFloat z = coord[2];
        const MFloat xi = surfCoordinates[0];
        const MFloat eta = surfCoordinates[1];
        const MFloat zeta = surfCoordinates[2];

        /*
        // Version without Prandtl-Glauert
        const MFloat X = x - xi;
        const MFloat Y = y - eta;
        const MFloat Z = z - zeta;
        const MFloat R = sqrt(POW2(X) + POW2(Y) + POW2(Z));

        const MFloat arg = waveNumber * R;
        const MFloat sinarg = sin(arg);
        const MFloat cosarg = cos(arg);
        const MFloat factor1 = waveNumber / (4.0 * PI * POW2(R));
        const MFloat factor2 = 1.0 / (4.0 * PI * POW3(R));
        const MFloat factorRe = -factor1 * sinarg - factor2 * cosarg;
        const MFloat factorIm = -factor1 * cosarg + factor2 * sinarg;

        const MFloat greenValue_Re = -cosarg / (4.0 * PI * R);
        const MFloat greenValue_Im = sinarg / (4.0 * PI * R);
        const MFloat dGreendXi_Re = X * factorRe;
        const MFloat dGreendXi_Im = X * factorIm;
        const MFloat dGreendEta_Re = Y * factorRe;
        const MFloat dGreendEta_Im = Y * factorIm;
        const MFloat dGreendZeta_Re = Z * factorRe;
        const MFloat dGreendZeta_Im = Z * factorIm;
        */

        // Version with Prandtl-Glauert, cf. Lockard2002
        const MFloat X = (x - xi) * cosAlpha * cosTheta + (y - eta) * sinAlpha + (z - zeta) * cosAlpha * sinTheta;
        // The Y coordinate differs from the paper, as there is a typo of the sign in front of the (z-zeta) term
        const MFloat Y = -(x - xi) * sinAlpha * cosTheta + (y - eta) * cosAlpha - (z - zeta) * sinAlpha * sinTheta;
        const MFloat Z = -(x - xi) * sinTheta + (z - zeta) * cosTheta;
        const MFloat d = std::sqrt(X * X + betaSq * (Y * Y + Z * Z));
        const MFloat k = waveNumber;

        // partial d / partial xi
        const MFloat F1B4pid = 1.0 / (4.0 * PI * d);
        const MFloat d_xi = (-X * cosAlpha * cosTheta + betaSq * (Y * sinAlpha * cosTheta + Z * sinTheta)) / d;
        const MFloat d_eta = (-X * sinAlpha - betaSq * Y * cosAlpha) / d;
        const MFloat d_zeta = (-X * cosAlpha * sinTheta + betaSq * (-Y * sinAlpha * sinTheta - Z * cosTheta)) / d;

        const MFloat arg = k * (mach * X - d) / betaSq;
        const MFloat cosArg = cos(arg);
        const MFloat sinArg = sin(arg);

        // Pre-terms in the derivatives
        const MFloat P0 = F1B4pid * mach * k / betaSq;
        const MFloat P1 = F1B4pid / d;
        const MFloat P2 = F1B4pid * k / betaSq;

        // clang-format off
        const MFloat greenValue_Re  = - F1B4pid * cosArg;
        const MFloat greenValue_Im  = - F1B4pid * sinArg;
        const MFloat dGreendXi_Re   =-P0 * sinArg * cosAlpha * cosTheta  + P1 * d_xi   * cosArg - P2 * d_xi   * sinArg;
        const MFloat dGreendXi_Im   = P0 * cosArg * cosAlpha * cosTheta  + P1 * d_xi   * sinArg + P2 * d_xi   * cosArg;
        const MFloat dGreendEta_Re  =-P0 * sinArg * sinAlpha             + P1 * d_eta  * cosArg - P2 * d_eta  * sinArg;
        const MFloat dGreendEta_Im  = P0 * cosArg * sinAlpha             + P1 * d_eta  * sinArg + P2 * d_eta  * cosArg;
        const MFloat dGreendZeta_Re =-P0 * sinArg * cosAlpha * sinTheta  + P1 * d_zeta * cosArg - P2 * d_zeta * sinArg;
        const MFloat dGreendZeta_Im = P0 * cosArg * cosAlpha * sinTheta  + P1 * d_zeta * sinArg + P2 * d_zeta * cosArg;
        // clang-format on

        const MFloat dL = m_surfaceData.surfaceArea(segmentId);

        const MFloat FxRe = m_surfaceData.complexVariables(segmentId, VarDimC[0], nw, 0);
        const MFloat FxIm = m_surfaceData.complexVariables(segmentId, VarDimC[0], nw, 1);

        const MFloat FyRe = m_surfaceData.complexVariables(segmentId, VarDimC[1], nw, 0);
        const MFloat FyIm = m_surfaceData.complexVariables(segmentId, VarDimC[1], nw, 1);

        const MFloat FzRe = m_surfaceData.complexVariables(segmentId, VarDimC[2], nw, 0);
        const MFloat FzIm = m_surfaceData.complexVariables(segmentId, VarDimC[2], nw, 1);

        const MFloat QRe = m_surfaceData.complexVariables(segmentId, VarDimC[3], nw, 0);
        const MFloat QIm = m_surfaceData.complexVariables(segmentId, VarDimC[3], nw, 1);

        // Calculation of integrals
        integralFxRe += dL * (FxRe * dGreendXi_Re - FxIm * dGreendXi_Im);
        integralFxIm += dL * (FxRe * dGreendXi_Im + FxIm * dGreendXi_Re);

        integralFyRe += dL * (FyRe * dGreendEta_Re - FyIm * dGreendEta_Im);
        integralFyIm += dL * (FyRe * dGreendEta_Im + FyIm * dGreendEta_Re);

        integralFzRe += dL * (FzRe * dGreendZeta_Re - FzIm * dGreendZeta_Im);
        integralFzIm += dL * (FzRe * dGreendZeta_Im + FzIm * dGreendZeta_Re);

        integralQRe += -1.0 * waveNumber * dL * (QRe * greenValue_Im + QIm * greenValue_Re);
        integralQIm += waveNumber * dL * (QRe * greenValue_Re - QIm * greenValue_Im);
      } // End of segmentId loop

      // sum up all integral terms
      const MFloat integralRe = -1.0 * (integralFxRe + integralFyRe + integralFzRe + integralQRe);
      const MFloat integralIm = -1.0 * (integralFxIm + integralFyIm + integralFzIm + integralQIm);

      // Saving Pressure in Observerdata variable p for each positive frequencie nw between 1 and
      // size/2+1. 0 = RE; 1 = IM
      p_complexVariables[2 * nw + 0] = integralRe;
      p_complexVariables[2 * nw + 1] = integralIm;
    });
  }
  // nw = 0 is not defined and therefore set to zero to avoid "NaN-errors"
  const MInt nw = 0;
  p_complexVariables[2 * nw + 0] = 0.0;
  p_complexVariables[2 * nw + 1] = 0.0;
  RECORD_TIMER_STOP(m_timers[Timers::CalcSurfaceIntegrals]);
}

/// \brief  Calculate the FWH surface integrals for a certain observer using the time-domain formulation
/// \author Zhe Yang
/// \date   15.05.2024
/// \Reference: Brentner, K. S. and F. Farassat, 1998, Analytical Comparison ... doi.org/10.2514/2.558
/// \param[in]  globalObserverId
/// \param[out] m_surfaceData.variables(segmentId, FWH_TIME::timeObs, tau)  observer time
/// \param[out] m_surfaceData.variables(segmentId, FWH_TIME::fwhPP, tau) observer pressure signal
template <MInt nDim>
void AcaSolver<nDim>::calcSurfaceIntegralsAtSourceTime(const MInt globalObserverId) {
  // Get observer coordinates
  const std::array<MFloat, nDim> coord = getGlobalObserverCoordinates(globalObserverId);

  // Calculate Mach number and the Prandtl-Glauert-Factor
  const MFloat mach = sqrt(std::inner_product(&m_MaVec[0], &m_MaVec[nDim], &m_MaVec[0], 0.0));
  const MFloat betaSq = 1.0 - mach * mach;

  // Loop over all the source timeSteps
  maia::parallelFor<false>(0, noSamples(), [=](MInt tau) {
    if constexpr(nDim == 3) {
      // SURFACE LOOP
      for(MInt segmentId = 0; segmentId < totalNoSurfaceElements(); segmentId++) {
        MFloat* surfCoordinates = &m_surfaceData.surfaceCoordinates(segmentId);

        // Geometric Distance & Coordinate Transformation
        MFloat distX = coord[0] - surfCoordinates[0];
        MFloat distY = coord[1] - surfCoordinates[1];
        MFloat distZ = coord[2] - surfCoordinates[2];
        transformCoordinatesToFlowDirection3D(&distX, &distY, &distZ);

        // Effective Acoustic Distance
        // Reference: Brès, G., 2010, A Ffowcs Williams - Hawkings ... doi.org/10.2514/6.2010-3711
        const MFloat R = sqrt(distX * distX + betaSq * (distY * distY + distZ * distZ));
        const MFloat dist = (R - mach * distX) / betaSq;

        const MFloat invsersDist = 1.0 / dist;
        const MFloat invsersDistSq = invsersDist * invsersDist;
        const MFloat distUnitX = (distX - mach * R) * invsersDist / betaSq;
        const MFloat distUnitY = distY * invsersDist;
        const MFloat distUnitZ = distZ * invsersDist;

        // compute observer time
        const MFloat timeSource = tau * m_dt;
        const MFloat timeObserver = timeSource + dist;

        // Get surface Mach number and source terms
        const MFloat Mx = m_surfaceData.variables(segmentId, FWH_TIME::surfMX, tau);
        const MFloat My = m_surfaceData.variables(segmentId, FWH_TIME::surfMY, tau);
        const MFloat Mz = m_surfaceData.variables(segmentId, FWH_TIME::surfMZ, tau);
        const MFloat srcUn = m_surfaceData.variables(segmentId, FWH_TIME::srcUn, tau);
        const MFloat srcLX = m_surfaceData.variables(segmentId, FWH_TIME::srcLX, tau);
        const MFloat srcLY = m_surfaceData.variables(segmentId, FWH_TIME::srcLY, tau);
        const MFloat srcLZ = m_surfaceData.variables(segmentId, FWH_TIME::srcLZ, tau);

        // Project variables on the radiation direction
        const MFloat Mr = distUnitX * Mx + distUnitY * My + distUnitZ * Mz;
        const MFloat inverseDopplerFactor = 1.0 / (1.0 - Mr);
        const MFloat inverseDopplerFactorSq = inverseDopplerFactor * inverseDopplerFactor;
        const MFloat Lr = distUnitX * srcLX + distUnitY * srcLY + distUnitZ * srcLZ;
        const MFloat LM = srcLX * Mx + srcLY * My + srcLZ * Mz;

        // Compute time derivative in second order accuracy
        const MFloat dotMX = calcTimeDerivative(segmentId, FWH_TIME::surfMX, tau);
        const MFloat dotMY = calcTimeDerivative(segmentId, FWH_TIME::surfMY, tau);
        const MFloat dotMZ = calcTimeDerivative(segmentId, FWH_TIME::surfMZ, tau);
        const MFloat dotUn = calcTimeDerivative(segmentId, FWH_TIME::srcUn, tau);
        const MFloat dotLX = calcTimeDerivative(segmentId, FWH_TIME::srcLX, tau);
        const MFloat dotLY = calcTimeDerivative(segmentId, FWH_TIME::srcLY, tau);
        const MFloat dotLZ = calcTimeDerivative(segmentId, FWH_TIME::srcLZ, tau);

        // Project time derivative on the radiation direction
        const MFloat dotMr = distUnitX * dotMX + distUnitY * dotMY + distUnitZ * dotMZ;
        const MFloat dotLr = distUnitX * dotLX + distUnitY * dotLY + distUnitZ * dotLZ;

        const MFloat MrFactor0 = invsersDist * inverseDopplerFactorSq;
        const MFloat MrFactor1 =
            (dist * dotMr + Mr - mach * mach) * invsersDistSq * inverseDopplerFactor * inverseDopplerFactorSq;

        // Compute Thickness Source
        const MFloat surfIntegralPPT_0 = dotUn * MrFactor0;
        const MFloat surfIntegralPPT_1 = srcUn * MrFactor1;

        // Compute Loading Source
        const MFloat surfIntegralPPL_0 = dotLr * MrFactor0;
        const MFloat surfIntegralPPL_1 = (Lr - LM) * MrFactor0 * invsersDist;
        const MFloat surfIntegralPPL_2 = Lr * MrFactor1;

        // Save the computed source value
        m_surfaceData.variables(segmentId, FWH_TIME::timeObs, tau) = timeObserver;
        m_surfaceData.variables(segmentId, FWH_TIME::fwhPP, tau) =
            surfIntegralPPT_0 + surfIntegralPPT_1 + surfIntegralPPL_0 + surfIntegralPPL_1 + surfIntegralPPL_2;
        m_surfaceData.variables(segmentId, FWH_TIME::fwhPP, tau) *= m_surfaceData.surfaceArea(segmentId);

        /*// Debugging
        if(segmentId == 0 && tau == 0) {
          m_log << "tau = " << tau
                << ", timeSource = " << timeSource
                << ", timeObserver = " << timeObserver
                << ", coord[0] = " << coord[0]
                << ", coord[1] = " << coord[1]
                << ", coord[2] = " << coord[2]
                << ", surfCoordinates[0] = " << surfCoordinates[0]
                << ", surfCoordinates[1] = " << surfCoordinates[1]
                << ", surfCoordinates[2] = " << surfCoordinates[2]
                << ", distX = " << distX
                << ", distY = " << distY
                << ", distZ = " << distZ
                << ", dist = " << dist
                << ", distUnitX = " << distUnitX
                << ", distUnitY = " << distUnitY
                << ", distUnitZ = " << distUnitZ
                << ", dL = " << m_surfaceData.surfaceArea(segmentId)
                << ", srcUn = " << srcUn
                << ", dotUn = " << dotUn
                << ", Lr = " << Lr
                << ", dotLr = " << dotLr
                << ", LM = " << LM
                << ", Mr = " << Mr
                << ", dotMr = " << dotMr
                << ", MSq = " << MSq
                << ", MrFactor0 = " << MrFactor0
                << ", MrFactor1 = " << MrFactor1
                << ", surfIntegralPPT_0= " << surfIntegralPPT_0
                << ", surfIntegralPPT_1= " << surfIntegralPPT_1
                << ", surfIntegralPPL_0= " << surfIntegralPPL_0
                << ", surfIntegralPPL_1= " << surfIntegralPPL_1
                << ", surfIntegralPPL_2= " << surfIntegralPPL_2
                << std::endl;
        }*/
      } // End of surface loop
    } else if constexpr(nDim == 2) {
      TERMM(1, "Time Domain FW-H not implimented in 2D: Requires free space Green function in 2D");
    } // End of the dimension switch
  }); // End of the time loop
}

template <MInt nDim>
void AcaSolver<nDim>::transformCoordinatesToFlowDirection2D(MFloat* const fX, MFloat* const fY) {
  if(m_zeroFlowAngle) return;

  // Project the vector on the flow direction U_inf = (U_inf, 0)
  const MFloat fX_old = *fX;
  const MFloat fY_old = *fY;

  *fX = fX_old * m_transformationMatrix[0] + fY_old * m_transformationMatrix[1];
  *fY = fX_old * m_transformationMatrix[2] + fY_old * m_transformationMatrix[3];
}

template <MInt nDim>
void AcaSolver<nDim>::transformCoordinatesToFlowDirection3D(MFloat* const fX, MFloat* const fY, MFloat* const fZ) {
  if(m_zeroFlowAngle) return;

  // Project the vector on the flow direction U_inf = (U_inf, 0, 0)
  const MFloat fX_old = *fX;
  const MFloat fY_old = *fY;
  const MFloat fZ_old = *fZ;

  *fX = fX_old * m_transformationMatrix[0] + fY_old * m_transformationMatrix[1] + fZ_old * m_transformationMatrix[2];
  *fY = fX_old * m_transformationMatrix[3] + fY_old * m_transformationMatrix[4] + fZ_old * m_transformationMatrix[5];
  *fZ = fX_old * m_transformationMatrix[6] + fY_old * m_transformationMatrix[7] + fZ_old * m_transformationMatrix[8];
}

template <MInt nDim>
void AcaSolver<nDim>::transformBackToOriginalCoordinate2D(MFloat* const fX, MFloat* const fY) {
  if(m_zeroFlowAngle) return;

  // Back transform to the original coordinate system
  // Use negative angular position (-theta) in m_transformationMatrix calculation
  const MFloat fX_old = *fX;
  const MFloat fY_old = *fY;

  *fX = fX_old * m_transformationMatrix[0] - fY_old * m_transformationMatrix[1];
  *fY = -fX_old * m_transformationMatrix[2] + fY_old * m_transformationMatrix[3];
}

template <MInt nDim>
void AcaSolver<nDim>::transformBackToOriginalCoordinate3D(MFloat* const fX, MFloat* const fY, MFloat* const fZ) {
  if(m_zeroFlowAngle) return;

  // Back transform to the original coordinate system
  // Use negative angular position (-alpha) and (-theta) in m_transformationMatrix calculation
  const MFloat fX_old = *fX;
  const MFloat fY_old = *fY;
  const MFloat fZ_old = *fZ;

  *fX = fX_old * m_transformationMatrix[0] - fY_old * m_transformationMatrix[1] - fZ_old * m_transformationMatrix[2];
  *fY = -fX_old * m_transformationMatrix[3] + fY_old * m_transformationMatrix[4] + fZ_old * m_transformationMatrix[5];
  *fZ = -fX_old * m_transformationMatrix[6] + fY_old * m_transformationMatrix[7] + fZ_old * m_transformationMatrix[8];
}

template <MInt nDim>
MFloat AcaSolver<nDim>::calcTimeDerivative(const MInt segmentId, const MInt varId, const MInt tau) {
  // Compute time derivative in second order accuarate
  MFloat dotVar = 0.0;

  if(tau == 0) {
    dotVar = -1.0 * m_surfaceData.variables(segmentId, varId, tau + 2)
             + 4.0 * m_surfaceData.variables(segmentId, varId, tau + 1)
             - 3.0 * m_surfaceData.variables(segmentId, varId, tau);
  } else if(tau == noSamples() - 1) {
    dotVar = +3.0 * m_surfaceData.variables(segmentId, varId, tau)
             - 4.0 * m_surfaceData.variables(segmentId, varId, tau - 1)
             + 1.0 * m_surfaceData.variables(segmentId, varId, tau - 2);
  } else {
    dotVar = +1.0 * m_surfaceData.variables(segmentId, varId, tau + 1)
             - 1.0 * m_surfaceData.variables(segmentId, varId, tau - 1);
  }
  return (dotVar / 2.0 / m_dt);
}

/// \brief Find the minimum distance from the surface elements to the observers
template <MInt nDim>
void AcaSolver<nDim>::calcMinDistanceForObservers(MFloat* const distMinObservers) {
  if(m_fwhTimeShiftType == 0) return;

  // Calculate Mach number and the Prandtl-Glauert-Factor
  const MFloat mach = sqrt(std::inner_product(&m_MaVec[0], &m_MaVec[nDim], &m_MaVec[0], 0.0));
  const MFloat betaSq = 1.0 - mach * mach;

  printMessage("- Loop over observers to find minimum distance from the source to the observers");
  for(MInt k = 0; k < globalNoObservers(); k++) {
    MFloat distMinLocal = std::numeric_limits<MFloat>::max();

    // Get observer coordinates
    const std::array<MFloat, nDim> coord = getGlobalObserverCoordinates(k);

    // SURFACE LOOP
    for(MInt segmentId = 0; segmentId < totalNoSurfaceElements(); segmentId++) {
      MFloat* surfCoordinates = &m_surfaceData.surfaceCoordinates(segmentId);

      if constexpr(nDim == 3) {
        // Geometric Distance & Coordinate Transformation
        MFloat distX = coord[0] - surfCoordinates[0];
        MFloat distY = coord[1] - surfCoordinates[1];
        MFloat distZ = coord[2] - surfCoordinates[2];
        transformCoordinatesToFlowDirection3D(&distX, &distY, &distZ);

        // Effective Acoustic Distance
        const MFloat R = sqrt(distX * distX + betaSq * (distY * distY + distZ * distZ));
        const MFloat dist = (R - mach * distX) / betaSq;

        // Find the minimum distance of all segments
        distMinLocal = (dist < distMinLocal) ? dist : distMinLocal;
      } else if constexpr(nDim == 2) {
        TERMM(1, "Time Domain FW-H not implimented in 2D: Requires free space Green function in 2D");
      } // End of the dimension switch
    }   // End of surface loop

    distMinObservers[k] = distMinLocal;
  } // End of observer loop

  MPI_Allreduce(MPI_IN_PLACE, &distMinObservers[0], globalNoObservers(), type_traits<MFloat>::mpiType(), MPI_MIN,
                mpiComm(), AT_, "MPI_IN_PLACE", "distMinObservers");

  if(m_fwhTimeShiftType == 2) {
    const MFloat disMinGlobal = *std::min_element(distMinObservers, distMinObservers + globalNoObservers());
    std::fill(distMinObservers, distMinObservers + globalNoObservers(), disMinGlobal);
    if(domainId() == 0) {
      std::cout << "Apply universal time shift for all observers = " << disMinGlobal << std::endl;
    }
  }
}

/// \brief Interpolate the observer pressure signal from source time to observer time using the advanced time approach
/// by Casalino. Reference: Casalino, D., 2003, An advanced time approach ... doi.org/10.1016/s0022-460x(02)00986-0
template <MInt nDim>
void AcaSolver<nDim>::interpolateFromSourceToObserverTime(MFloat distMinObservers, MFloat* const p_perturbedFWH) {
  const MFloat inverseFourPI = 1.0 / 4.0 / PI;

  // Loop over all the observer timeSteps
  maia::parallelFor<false>(0, noSamples(), [=](MInt t) {
    // Set the observer time
    MFloat timeObserver = t * m_dt + distMinObservers;

    // initialize surface integral
    MFloat sumFwhPP = 0.0;

    // Surface element loop
    for(MInt segmentId = 0; segmentId < totalNoSurfaceElements(); segmentId++) {
      // initialize interpolation value
      MFloat interpolateFwhPP = 0.0;
      MInt t_lower = 0;
      MInt t_upper = noSamples() - 1;

      // Read the observer lower/upper bound
      MFloat timeObserver_lower = m_surfaceData.variables(segmentId, FWH_TIME::timeObs, t_lower);
      MFloat timeObserver_upper = m_surfaceData.variables(segmentId, FWH_TIME::timeObs, t_upper);

      // Search in FWH_TIME::timeObs
      if(timeObserver < timeObserver_lower || timeObserver > timeObserver_upper) {
        // In case the timeObserver is out of bound
        interpolateFwhPP = 0.0;
      } else if(approx(timeObserver, timeObserver_lower, MFloatEps)) {
        // In case the timeObserver is at the lower bound
        interpolateFwhPP = m_surfaceData.variables(segmentId, FWH_TIME::fwhPP, t_lower);
      } else if(approx(timeObserver, timeObserver_upper, MFloatEps)) {
        // In case the timeObserver is at the upper bound
        interpolateFwhPP = m_surfaceData.variables(segmentId, FWH_TIME::fwhPP, t_upper);
      } else {
        // In the case of: timeObserver_lower < timeObserver < timeObserver_upper
        // Use binary search algorithm to search

        MInt t_range = t_upper - t_lower;
        do {
          // Compute t_mid
          const MFloat t_mid = (t_range % 2 == 0) ? ((t_lower + t_upper) / 2) : ((t_lower + t_upper - 1) / 2);
          const MFloat timeObserver_mid = m_surfaceData.variables(segmentId, FWH_TIME::timeObs, t_mid);

          if(approx(timeObserver, timeObserver_mid, MFloatEps)) {
            // In case the timeObserver is at the middle
            interpolateFwhPP = m_surfaceData.variables(segmentId, FWH_TIME::fwhPP, t_mid);
            break;
          }

          // Update t_upper and t_lower
          if(timeObserver < timeObserver_mid) {
            t_upper = t_mid;
          } else if(timeObserver > timeObserver_mid) {
            t_lower = t_mid;
          }
          t_range = t_upper - t_lower;

          // For the smallest un-dividable region, do linear interpolation to get the observer fwhPP
          if(t_range == 1) {
            timeObserver_lower = m_surfaceData.variables(segmentId, FWH_TIME::timeObs, t_lower);
            timeObserver_upper = m_surfaceData.variables(segmentId, FWH_TIME::timeObs, t_upper);
            const MFloat fwhPP_lower = m_surfaceData.variables(segmentId, FWH_TIME::fwhPP, t_lower);
            const MFloat fwhPP_upper = m_surfaceData.variables(segmentId, FWH_TIME::fwhPP, t_upper);
            const MFloat fwhPP_slop = (fwhPP_upper - fwhPP_lower) / (timeObserver_upper - timeObserver_lower);
            interpolateFwhPP = fwhPP_lower + (timeObserver - timeObserver_lower) * fwhPP_slop;
            break;
          }
        } while(t_range > 1);
      } // End of the search in FWH_TIME::timeObs

      sumFwhPP += interpolateFwhPP;
    } // End of the surface loop

    p_perturbedFWH[t] = inverseFourPI * sumFwhPP;
  }); // End of the observer time loop (parallel)
}

/// \brief Combine the FWH surface integrals of all ranks
/// \param[in] globalObserverId   Global observer id
/// \param[in] p_complexVariables Output of complex p', length of 2*noSamples()
template <MInt nDim>
void AcaSolver<nDim>::combineSurfaceIntegrals(const MInt globalObserverId, MFloat* const p_complexVariables) {
  TRACE();
  RECORD_TIMER_START(m_timers[Timers::CombineSurfaceIntegrals]);

  // Data is stored in unconvenient way (globalObserverId, 0) so that Allreduce won't work. Put Data in linear
  // sequence --> Allreduce --> Back in the convenient data structure
  // (globalObserverId, 0)

  // First: Copy all Observerdata to a databuffer
  // obsDataSize is noSamples, because we have imaginary and real part but only for half of samples
  // (all positive frequencies including zero)
  constexpr MInt noRealAndIm = 2;
  const MInt halfNoSamplesRealAndIm = (noSamples() / 2 + 1) * noRealAndIm;
  const MInt dataBuffSize = halfNoSamplesRealAndIm * 1.0;
  std::vector<MFloat> obsDataBuff(dataBuffSize, 0.0);

  const MInt trgDomain = getObserverDomainId(globalObserverId);
  const MInt observerId = a_localObserverId(globalObserverId);

  TERMM_IF_COND(trgDomain == domainId() && (observerId < 0 || observerId > noObservers()),
                "Error in observer partitioning!");

  // The actual copying: Whole observer Data from observer "globalObserverId" to the position globalObserverId *
  // obsDataSize in databuffer
  std::copy_n(p_complexVariables, halfNoSamplesRealAndIm, obsDataBuff.data());
  // obsDataBuff now contains all complex pressure values of all observers in linear sequence

  if(domainId() == trgDomain) {
    // Sums up all parts of the surface integral and distributes the results back to all processes.
    MPI_Reduce(MPI_IN_PLACE, obsDataBuff.data(), dataBuffSize, type_traits<MFloat>::mpiType(), MPI_SUM, trgDomain,
               mpiComm(), AT_, "MPI_IN_PLACE", "obsDataBuff");
    // obsDataBuff now contains the FWH-Integral value of each observer and each frequency in linear
    // sequence.

    // Split the data buffer into the convenient data structure
    std::copy_n(obsDataBuff.data(), halfNoSamplesRealAndIm, &m_observerData.complexVariables(observerId, 0));
  } else {
    MPI_Reduce(obsDataBuff.data(), obsDataBuff.data(), dataBuffSize, type_traits<MFloat>::mpiType(), MPI_SUM, trgDomain,
               mpiComm(), AT_, "MPI_IN_PLACE", "obsDataBuff");
  }

  // At first the segments were split up on processers. Each processes computed
  // the integral for all observers for his segments. Then all integral were
  // summed up. Now, all processes have the correct (local) observer data.

  RECORD_TIMER_STOP(m_timers[Timers::CombineSurfaceIntegrals]);
}

/// \brief Combine the FWH surface integrals of the time domain method of all ranks
template <MInt nDim>
void AcaSolver<nDim>::combineSurfaceIntegralsInTime(const MInt globalObserverId, MFloat* const p_perturbedFWH) {
  TRACE();

  // Data is stored in unconvenient way (globalObserverId, 0) so that Allreduce won't work. Put Data in linear
  // sequence --> Allreduce --> Back in the convenient data structure (globalObserverId, 0)

  // First: Copy all Observerdata to a databuffer
  const MInt dataBuffSize = noSamples();
  std::vector<MFloat> obsDataBuff(dataBuffSize, 0.0);

  const MInt trgDomain = getObserverDomainId(globalObserverId);
  const MInt observerId = a_localObserverId(globalObserverId);

  TERMM_IF_COND(trgDomain == domainId() && (observerId < 0 || observerId > noObservers()),
                "Error in observer partitioning!");

  // The actual copying: Whole observer Data from observer "globalObserverId" to
  // the position globalObserverId * obsDataSize in databuffer
  std::copy_n(p_perturbedFWH, dataBuffSize, obsDataBuff.data());
  // obsDataBuff now contains all complex pressure values of all observers in linear sequence

  if(domainId() == trgDomain) {
    // Sums up all parts of the surface integral and distributes the results back to all processes.
    MPI_Reduce(MPI_IN_PLACE, obsDataBuff.data(), dataBuffSize, type_traits<MFloat>::mpiType(), MPI_SUM, trgDomain,
               mpiComm(), AT_, "MPI_IN_PLACE", "obsDataBuff");
    // obsDataBuff now contains the FWH-Integral value of each observer and each frequency in linear sequence.

    // Split the data buffer into the convenient data structure
    std::copy_n(obsDataBuff.data(), dataBuffSize, &m_observerData.variables(observerId, 0));
  } else {
    MPI_Reduce(obsDataBuff.data(), obsDataBuff.data(), dataBuffSize, type_traits<MFloat>::mpiType(), MPI_SUM, trgDomain,
               mpiComm(), AT_, "MPI_IN_PLACE", "obsDataBuff");
  }

  // At first the segments were split up on processers. Each processes computed
  // the integral for all observers for his segments. Then all integral were
  // summed up. Now, all processes have the correct (local) observer data.
}

/// \brief Calculate the observer signals in frequency domain
template <MInt nDim>
void AcaSolver<nDim>::calculateObserverInFrequency() {
  printMessage("- Calculate observer signals in frequency domain");
  ScratchSpace<MFloat> complexVariablesBuffer(2 * noSamples(), AT_, "complexVariablesBuffer");
  complexVariablesBuffer.fill(0.0);
  MFloat* const p_complexVariables = complexVariablesBuffer.data();
  for(MInt globalObserverId = 0; globalObserverId < globalNoObservers(); globalObserverId++) {
    std::array<MFloat, nDim> coord = getGlobalObserverCoordinates(globalObserverId);
    calcSurfaceIntegralsForObserver(coord, p_complexVariables);
    applySymmetryBc(coord, p_complexVariables);
    combineSurfaceIntegrals(globalObserverId, p_complexVariables);
  }
}

/// \brief Calculate the observer signals in time domain
template <MInt nDim>
void AcaSolver<nDim>::calculateObserverInTime() {
  // Set memory for min observer distance
  ScratchSpace<MFloat> distMinBuffer(globalNoObservers(), AT_, "distMinBuffer");
  distMinBuffer.fill(0.0);
  MFloat* const distMinObservers = distMinBuffer.data();

  // Set memory for observer pressure signal
  ScratchSpace<MFloat> perturbedPressureBuffer(noSamples(), AT_, "perturbedPressureBuffer");
  perturbedPressureBuffer.fill(0.0);
  MFloat* const p_perturbedFWH = perturbedPressureBuffer.data();

  printMessage("- Calculate observer signals in time domain");
  calcMinDistanceForObservers(distMinObservers);

  for(MInt globalObserverId = 0; globalObserverId < globalNoObservers(); globalObserverId++) {
    // m_log << "Compute for observer = " << globalObserverId << "/" << globalNoObservers() << std::endl;
    calcSurfaceIntegralsAtSourceTime(globalObserverId);
    interpolateFromSourceToObserverTime(distMinObservers[globalObserverId], p_perturbedFWH);
    combineSurfaceIntegralsInTime(globalObserverId, p_perturbedFWH);
  }
}

/// \brief Evenly distribute observers on all ranks
template <MInt nDim>
void AcaSolver<nDim>::distributeObservers() {
  TRACE();
  const MInt minNoObsPerDomain = globalNoObservers() / noDomains();
  const MInt remainderObs = globalNoObservers() % noDomains();
  if(domainId() < remainderObs) {
    m_noObservers = minNoObsPerDomain + 1;
    m_offsetObserver = domainId() * (minNoObsPerDomain + 1);
  } else {
    m_noObservers = minNoObsPerDomain;
    m_offsetObserver = domainId() * minNoObsPerDomain + remainderObs;
  }
}

/// \brief Backtransform the observer frequency signals into the time domain.
template <MInt nDim>
void AcaSolver<nDim>::backtransformObservers() {
  TRACE();
  RECORD_TIMER_START(m_timers[Timers::BacktransformObservers]);
  printMessage("- Backtransforming the observer signals");

  // Data buffer for backtransformed pressure data
  MFloatScratchSpace backtransform(noSamples() * 2, AT_, "backtransform");

  // BACKTRANSFORM
  for(MInt obs = 0; obs < noObservers(); obs++) {
    // Copy the first calculated half to the other half starting from behind going to the middle
    // (since it is a two-sided frequency plot)
    for(MInt i = 1; i < noSamples() / 2; i++) {
      m_observerData.complexVariables(obs, 0, noSamples() - i, 0) = m_observerData.complexVariables(obs, 0, i, 0);
      m_observerData.complexVariables(obs, 0, noSamples() - i, 1) =
          -1.0 * m_observerData.complexVariables(obs, 0, i, 1);
    }

    std::fill_n(&backtransform[0], 2 * noSamples(), 0.0);
    // Backtransformation
    if(m_FastFourier == true && m_transformationType == 1) {
      FastFourierTransform(&m_observerData.complexVariables(obs, 0), noSamples(), -1, &backtransform[0], false);
    } else if(m_FastFourier == false || (m_FastFourier == true && m_transformationType == 0)) {
      DiscreteFourierTransform(&m_observerData.complexVariables(obs, 0), noSamples(), -1, &backtransform[0], false);
    }

    // Storing the results in observer.variables since it is back in time domain (only real part)
    for(MInt i = 0; i < noSamples(); i++) {
      m_observerData.variables(obs, 0, i) = backtransform[2 * i];
    }
  }

  RECORD_TIMER_STOP(m_timers[Timers::BacktransformObservers]);
}


/// \brief Save the observer signals.
/// Save observer data: 1. pressure data in time domain and 2. pressure data in
/// frequency domain.
template <MInt nDim>
void AcaSolver<nDim>::saveObserverSignals() {
  TRACE();
  RECORD_TIMER_START(m_timers[Timers::SaveObserverSignals]);

  printMessage("- Save observer signals");

  using namespace maia::parallel_io;
  const MString fileName = outputDir() + m_outputFilePrefix + "observerData" + ParallelIo::fileExt();

  ParallelIo file(fileName, PIO_REPLACE, mpiComm());

  ParallelIo::size_type dimSizesCoordinates[] = {globalNoObservers(), nDim};
  file.defineArray(PIO_FLOAT, "coordinates", 2, &dimSizesCoordinates[0]);

  file.defineArray(PIO_FLOAT, "time", noSamples());
  file.defineArray(PIO_FLOAT, "frequency", noSamples());

  ParallelIo::size_type dimSizesVarTime[] = {globalNoObservers(), noSamples()};
  file.defineArray(PIO_FLOAT, "pressure_time", 2, &dimSizesVarTime[0]);

  ParallelIo::size_type dimSizesVarFreq[] = {globalNoObservers(), 2 * noSamples()};
  file.defineArray(PIO_FLOAT, "pressure_frequency", 2, &dimSizesVarFreq[0]);

  // Store main parameters/configuration as attributes
  const MString acousticMethod = (m_acousticMethod == FWH_METHOD) ? "FWH" : "FWH_APE";
  file.setAttribute(acousticMethod, "acousticMethod");
  file.setAttribute(m_Ma, "Ma");
  if(!approx(m_Ma, m_MaDim, MFloatEps)) {
    file.setAttribute(m_MaDim, "Ma_dim");
  }
  file.setAttributes(m_MaVec, "Ma_i", nDim);
  // Sample attributes
  file.setAttribute(m_noSamples, "noSamples");
  file.setAttribute(m_sampleStride, "sampleStride");
  file.setAttribute(m_sampleOffset, "sampleOffset");
  file.setAttribute(m_dt, "dt");
  // Windowing
  file.setAttribute(WINDOW::windowNames[m_windowType], "windowType");
  file.setAttribute(m_windowNormFactor, "windowNormFactor");
  // Surface attributes
  for(MInt sId = 0; sId < noSurfaces(); sId++) {
    file.setAttribute(m_surfaceInputFileName[sId], "inputFileName_" + std::to_string(m_surfaceIds[sId]));
    file.setAttribute(m_noSurfaceElementsGlobal[sId], "noSurfaceElements_" + std::to_string(m_surfaceIds[sId]));
    file.setAttribute(m_surfaceWeightingFactor[sId], "surfaceWeight_" + std::to_string(m_surfaceIds[sId]));
  }
  file.setAttribute(m_observerFileName, "observerFileName");

  // Store also observer coordinates to have everything in one file
  if(domainId() == 0) {
    file.setOffset(globalNoObservers(), 0, 2);
  } else {
    file.setOffset(0, 0, 2);
  }
  file.writeArray(m_globalObserverCoordinates.data(), "coordinates");

  // write times and frequencies
  if(domainId() == 0) {
    file.setOffset(noSamples(), 0);
  } else {
    file.setOffset(0, 0);
  }
  file.writeArray(&m_times[0], "time");
  file.writeArray(&m_frequencies[0], "frequency");

  file.setOffset(noObservers(), observerOffset(), 2);
  ScratchSpace<MFloat> pressure_time(std::max(noObservers(), 1), noSamples(), AT_, "pressure_time");
  for(MInt obs = 0; obs < noObservers(); obs++) {
    for(MInt t = 0; t < noSamples(); t++) {
      pressure_time(obs, t) = m_observerData.variables(obs, 0, t);
    }
  }
  file.writeArray(&pressure_time[0], "pressure_time");

  file.setOffset(noObservers(), observerOffset(), 2);
  ScratchSpace<MFloat> pressure_frequency(std::max(noObservers(), 1), 2 * noSamples(), AT_, "pressure_frequency");
  for(MInt obs = 0; obs < noObservers(); obs++) {
    for(MInt f = 0; f < noSamples(); f++) {
      pressure_frequency(obs, 2 * f) = m_observerData.complexVariables(obs, 0, f, 0);
      pressure_frequency(obs, 2 * f + 1) = m_observerData.complexVariables(obs, 0, f, 1);
    }
  }
  file.writeArray(&pressure_frequency[0], "pressure_frequency");

  RECORD_TIMER_STOP(m_timers[Timers::SaveObserverSignals]);
}


/// \brief Load the observer signals from a file for further postprocessing.
template <MInt nDim>
void AcaSolver<nDim>::loadObserverSignals() {
  TRACE();
  RECORD_TIMER_START(m_timers[Timers::LoadObserverSignals]);
  using namespace maia::parallel_io;

  const MInt obsOffset = observerOffset();
  const MInt obsCount = noObservers();

  std::stringstream ss;
  ss << "- Load observer signals from file: " << m_acaPostprocessingFile << std::endl;
  printMessage(ss.str());
  ParallelIo file(m_acaPostprocessingFile, PIO_READ, mpiComm());

  // read times and frequencies on rank 0 and broadcast
  if(domainId() == 0) {
    file.setOffset(noSamples(), 0);
  } else {
    file.setOffset(0, 0);
  }
  file.readArray(&m_times[0], "time");
  file.readArray(&m_frequencies[0], "frequency");

  MPI_Bcast(&m_times[0], noSamples(), type_traits<MFloat>::mpiType(), 0, mpiComm(), AT_, "&m_times[0]");
  MPI_Bcast(&m_frequencies[0], noSamples(), type_traits<MFloat>::mpiType(), 0, mpiComm(), AT_, "&m_frequencies[0]");

  // Read observer data in time domain
  file.setOffset(obsCount, obsOffset, 2);
  ScratchSpace<MFloat> pressure_time(std::max(obsCount, 1), noSamples(), AT_, "pressure_time");

  file.readArray(&pressure_time[0], "pressure_time");
  for(MInt obs = 0; obs < obsCount; obs++) {
    for(MInt t = 0; t < noSamples(); t++) {
      m_observerData.variables(obs, 0, t) = pressure_time(obs, t);
    }
  }

  // Read observer data in frequency domain
  file.setOffset(obsCount, obsOffset, 2);
  ScratchSpace<MFloat> pressure_frequency(std::max(obsCount, 1), 2 * noSamples(), AT_, "pressure_frequency");
  file.readArray(&pressure_frequency[0], "pressure_frequency");

  for(MInt obs = 0; obs < obsCount; obs++) {
    for(MInt f = 0; f < noSamples(); f++) {
      m_observerData.complexVariables(obs, 0, f, 0) = pressure_frequency(obs, 2 * f);
      m_observerData.complexVariables(obs, 0, f, 1) = pressure_frequency(obs, 2 * f + 1);
    }
  }

  RECORD_TIMER_STOP(m_timers[Timers::LoadObserverSignals]);
}


/// \brief Postprocessing of observer signals.
template <MInt nDim>
void AcaSolver<nDim>::postprocessObserverSignals() {
  TRACE();
  RECORD_TIMER_START(m_timers[Timers::Postprocessing]);
  printMessage("- Post process observer signals");

  // Loop over all postprocessing operations
  for(MInt op = 0; op < m_noPostprocessingOps; op++) {
    const MInt operation = m_postprocessingOps[op];
    switch(operation) {
      case PP::RMS_PRESSURE: {
        m_post.push_back(std::make_unique<AcaPostProcessingRMS<nDim>>(mpiComm()));
        break;
      }
      case PP::OASPL: {
        m_post.push_back(std::make_unique<AcaPostProcessingOASPL<nDim>>(mpiComm()));
        break;
      }
      case PP::SPL: {
        m_post.push_back(std::make_unique<AcaPostProcessingSPL<nDim>>(mpiComm()));
        break;
      }
      case PP::pABS: {
        m_post.push_back(std::make_unique<AcaPostProcessingPABS<nDim>>(mpiComm()));
        break;
      }
      case PP::calcObsPress: {
        TERMM_IF_NOT_COND(m_generateSurfaceData, "Only useful for analytical testcases with generated surface data.");
        calcObsPressureAnalytic();
        break;
      }
      default: {
        std::cerr << "  "
                  << "skippinginvalid postprocessing operation: " << operation << std::endl;
        break;
      }
    }
  }

  for(auto&& post : m_post) {
    post->init(noObservers(), globalNoObservers(), observerOffset(), noSamples(), m_globalObserverCoordinates.data(),
               m_frequencies.data(), outputDir() + m_outputFilePrefix);
    for(MInt obs = 0; obs < noObservers(); obs++) {
      post->calc(obs, &m_observerData.variables(obs, 0), &m_observerData.complexVariables(obs, 0));
    }
    post->finish();
  }

  RECORD_TIMER_STOP(m_timers[Timers::Postprocessing]);
}

// Only useful for analytical test cases
template <MInt nDim>
void AcaSolver<nDim>::calcObsPressureAnalytic() {
  TRACE();

  const MInt noVars = (m_acousticMethod == FWH_METHOD) ? nDim + 2 : nDim + 1;
  MFloatScratchSpace vars(noVars, noObservers(), noSamples(), AT_, "varP");
  vars.fill(0.0);

  for(MInt obsId = 0; obsId < noObservers(); obsId++) {
    // Pointers to pass to functions
    const MFloat* obsCoord = &m_observerData.observerCoordinates(obsId);
    SourceVars sourceVars;
    if(m_acousticMethod == FWH_METHOD) {
      sourceVars.p = &vars(FWH::P, obsId, 0);
      sourceVars.u = &vars(FWH::U, obsId, 0);
      sourceVars.v = &vars(FWH::V, obsId, 0);
      sourceVars.w = &vars(FWH::W, obsId, 0);
      sourceVars.rho = &vars(FWH::RHO, obsId, 0);
    } else {
      sourceVars.p = &vars(FWH_APE::PP, obsId, 0);
      sourceVars.u = &vars(FWH_APE::UP, obsId, 0);
      sourceVars.v = &vars(FWH_APE::VP, obsId, 0);
      sourceVars.w = &vars(FWH_APE::WP, obsId, 0);
    }

    if constexpr(nDim == 2) {
      switch(m_sourceType) {
        case 0: {
          genMonopoleAnalytic2D(obsCoord, m_sourceParameters, sourceVars);
          break;
        }
        case 1: {
          genDipoleAnalytic2D(obsCoord, m_sourceParameters, sourceVars);
          break;
        }
        case 2: {
          genQuadrupoleAnalytic2D(obsCoord, m_sourceParameters, sourceVars);
          break;
        }
        case 3: {
          genVortexConvectionAnalytic2D(obsCoord, m_sourceParameters, sourceVars);
          break;
        }
        default: {
          TERMM(1, "source type not implemented");
          break;
        }
      }
    } else {
      switch(m_sourceType) {
        case 0: {
          genMonopoleAnalytic3D(obsCoord, m_sourceParameters, sourceVars);
          break;
        }
        case 1: {
          genDipoleAnalytic3D(obsCoord, m_sourceParameters, sourceVars);
          break;
        }
        case 2: {
          genQuadrupoleAnalytic3D(obsCoord, m_sourceParameters, sourceVars);
          break;
        }
        default: {
          TERMM(1, "source type not implemented");
          break;
        }
      }
    }
  }

  // OUTPUT
  using namespace maia::parallel_io;
  const MString fileName = outputDir() + m_outputFilePrefix + "obsPressAnalytic" + ParallelIo::fileExt();

  ParallelIo file(fileName, PIO_REPLACE, mpiComm());

  ParallelIo::size_type dimSizes[] = {globalNoObservers(), noSamples()};
  file.defineArray(PIO_FLOAT, "obsPressAnalytic", 2, &dimSizes[0]);

  file.setOffset(noObservers(), observerOffset(), 2);
  if(m_acousticMethod == FWH_METHOD) {
    file.writeArray(&vars(FWH::P, 0, 0), "obsPressAnalytic");
  } else {
    file.writeArray(&vars(FWH_APE::PP, 0, 0), "obsPressAnalytic");
  }
}

/// \brief Initialize the parallelization, i.e., distribute all surface segments among ranks
///        create communicators etc.
template <MInt nDim>
void AcaSolver<nDim>::initParallelization() {
  TRACE();
  printMessage("- Init parallelization");
  RECORD_TIMER_START(m_timers[Timers::InitParallelization]);

  m_noSurfaceElements.resize(noSurfaces());
  m_noSurfaceElementsGlobal.resize(noSurfaces());
  m_surfaceElementOffsets.resize(noSurfaces());
  m_surfaceInputFileName.resize(noSurfaces(), "");

  std::vector<MInt> totalNoSegments(noSurfaces());

  // TODO labels:ACA maybe do this only on rank 0
  if(m_generateSurfaceData) {
    for(MInt sId = 0; sId < noSurfaces(); sId++) {
      std::vector<MInt> dummyVec{};
      // Determine total number of segments for each surface
      // const MInt totalNoSegments = readNoSegmentsSurfaceGeneration(sId, dummyVec);
      totalNoSegments[sId] = readNoSegmentsSurfaceGeneration(sId, dummyVec);
    }
  } else {
    m_totalNoSamples = -1;

    // Determine total number of samples and number of segments for each surface from input data files
    for(MInt sId = 0; sId < noSurfaces(); sId++) {
      MInt noSamples = -1;
      readNoSegmentsAndSamples(sId, totalNoSegments[sId], noSamples, m_surfaceInputFileName[sId]);

      // Check that number of samples matches for all surfaces
      TERMM_IF_NOT_COND(m_totalNoSamples < 0 || m_totalNoSamples == noSamples,
                        "number of samples mismatch for different surfaces");
      m_totalNoSamples = noSamples;
    }

    m_noSamples = (m_noSamples == -1) ? m_totalNoSamples : m_noSamples;

    // Check noSamples/stride/offset
    const MInt usedSampleRange = m_noSamples * m_sampleStride + m_sampleOffset;
    TERMM_IF_COND(m_noSamples < 1 || m_sampleStride < 1 || m_sampleOffset < 0 || usedSampleRange > m_totalNoSamples,
                  "Error: number of samples and stride need to be > 0, offset >= 0 and samples*stride+offset cannot "
                  "exceed the total number of samples (noSamples="
                      + std::to_string(m_noSamples) + ", stride=" + std::to_string(m_sampleStride) + ", offset="
                      + std::to_string(m_sampleOffset) + ", samples*stride+offset=" + std::to_string(usedSampleRange)
                      + ", totalNoSamples=" + std::to_string(m_totalNoSamples) + ").");

    // Only even number of samples
    if(m_noSamples % 2 > 0) {
      m_noSamples = m_noSamples - 1;
      if(domainId() == 0) {
        std::cout << "Note: Number of Samples was uneven which is not feasible. The new number of samples is: "
                  << m_noSamples << std::endl;
      }
    }
  }

  // Store global number of segments for all surfaces
  std::copy(totalNoSegments.begin(), totalNoSegments.end(), &m_noSurfaceElementsGlobal[0]);

  // Determine distribution of surface segments among all processes
  const MInt globalNoSegments = std::accumulate(totalNoSegments.begin(), totalNoSegments.end(), 0);

  for(MInt sId = 0; sId < noSurfaces(); sId++) {
    m_log << "Surface #" << m_surfaceIds[sId] << " has " << totalNoSegments[sId]
          << " segments; surface input file name: " << m_surfaceInputFileName[sId]
          << "; weighting factor: " << m_surfaceWeightingFactor[sId] << std::endl;
  }
  m_log << "Global number of segments: " << globalNoSegments << std::endl;
  m_log << "Number of samples: " << m_noSamples << " (totalNoSamples = " << m_totalNoSamples
        << ", stride = " << m_sampleStride << ", offset = " << m_sampleOffset << ")" << std::endl;

  MInt localNoSegments = -1;

  // Distribute all segments equally among all domains. This means that a process can have segments
  // of different surfaces, e.g., when using only a small number of processes. However, for larger
  // production cases it makes sense to limit each process to have only segments of a single surface
  // if the I/O of data of different surfaces should be fully parallel.
  if(m_allowMultipleSurfacesPerRank) {
    // Average number of segments (rounded down)
    const MInt avgNoSegments = floor(static_cast<MFloat>(globalNoSegments) / static_cast<MFloat>(noDomains()));
    localNoSegments = avgNoSegments;
    // Equally distribute missing segments
    if(domainId() < globalNoSegments - noDomains() * avgNoSegments) {
      localNoSegments += 1;
    }
  } else { // Restrict each rank to have only segments of one surface to improve performance in production cases
    std::vector<MInt> noDomainsPerSurface(noSurfaces());
    // Average number of segments
    const MFloat avgNoSegments = static_cast<MFloat>(globalNoSegments) / static_cast<MFloat>(noDomains());
    // Distribute ranks on surfaces
    for(MInt i = 0; i < noSurfaces(); i++) {
      noDomainsPerSurface[i] = std::max(1, static_cast<MInt>(std::floor(totalNoSegments[i] / avgNoSegments)));
    }
    // Assign remaining ranks to the surfaces with the highest number of segments per rank
    while(noDomains() - std::accumulate(noDomainsPerSurface.begin(), noDomainsPerSurface.end(), 0) > 0) {
      MFloat maxWeight = 0.0;
      MInt maxIndex = -1;
      for(MInt i = 0; i < noSurfaces(); i++) {
        const MFloat weight = totalNoSegments[i] / noDomainsPerSurface[i];
        if(weight > maxWeight) {
          maxWeight = weight;
          maxIndex = i;
        }
      }
      noDomainsPerSurface[maxIndex]++;
    }
    const MInt noDomainsSum = std::accumulate(noDomainsPerSurface.begin(), noDomainsPerSurface.end(), 0);
    TERMM_IF_NOT_COND(noDomainsSum == noDomains(),
                      "number of domains mismatch " + std::to_string(noDomainsSum) + " " + std::to_string(noDomains()));

    // Determine the local surface
    MInt localSurfaceId = -1;
    MInt localSurfaceDomainId = -1;
    for(MInt i = 0, sumDomains = 0; i < noSurfaces(); i++) {
      if(sumDomains + noDomainsPerSurface[i] > domainId()) {
        localSurfaceId = i;
        localSurfaceDomainId = domainId() - sumDomains;
        break;
      }
      sumDomains += noDomainsPerSurface[i];
    }

    const MInt avgNoSegmentsSurface = floor(totalNoSegments[localSurfaceId] / noDomainsPerSurface[localSurfaceId]);
    // Determine local number of segments of this surface
    localNoSegments = avgNoSegmentsSurface;
    if(localSurfaceDomainId
       < totalNoSegments[localSurfaceId] - noDomainsPerSurface[localSurfaceId] * avgNoSegmentsSurface) {
      localNoSegments += 1;
    }
  }

  // Communicate number of segments per domain
  std::vector<MInt> noSegmentsPerDomain(noDomains());
  MPI_Allgather(&localNoSegments, 1, type_traits<MInt>::mpiType(), &noSegmentsPerDomain[0], 1,
                type_traits<MInt>::mpiType(), mpiComm(), AT_, "localNoSegments", "noSegmentsPerDomain[0]");

  // Sanity check
  const MInt globalCount = std::accumulate(noSegmentsPerDomain.begin(), noSegmentsPerDomain.end(), 0);
  TERMM_IF_NOT_COND(globalCount == globalNoSegments, "globalNoSegments mismatch");

  // Compute domain offsets
  std::vector<MInt> domainSegmentOffsets(std::max(noDomains() + 1, 1));
  domainSegmentOffsets[0] = 0;
  for(MInt i = 1; i < noDomains() + 1; i++) {
    domainSegmentOffsets[i] = domainSegmentOffsets[i - 1] + noSegmentsPerDomain[i - 1];
  }

  MIntScratchSpace noSegments_(noDomains(), noSurfaces(), AT_, "noSegments_");

  MInt surfaceOffset = 0;
  // Check which surfaces are present on the local rank
  for(MInt sId = 0; sId < noSurfaces(); sId++) {
    const MInt nextSurfaceOffset = surfaceOffset + totalNoSegments[sId];
    const MInt start = domainSegmentOffsets[domainId()];
    const MInt end = domainSegmentOffsets[domainId()] + localNoSegments;

    MInt noLocalSegments = 0;
    MInt localSegmentsOffset = 0;

    if(start < surfaceOffset && end >= surfaceOffset) {
      // First local segment belongs to a previous surface but last local segment belongs to current
      // or a subsequent surface
      noLocalSegments = std::min(totalNoSegments[sId], end - surfaceOffset);
      localSegmentsOffset = 0;
    } else if(start >= surfaceOffset && start < nextSurfaceOffset) {
      // First local segment belongs to this surface
      noLocalSegments = std::min(nextSurfaceOffset - start, end - start);
      localSegmentsOffset = start - surfaceOffset;
    }

    m_noSurfaceElements[sId] = noLocalSegments;
    m_surfaceElementOffsets[sId] = localSegmentsOffset;

    // Store in communication buffer
    noSegments_(domainId(), sId) = noLocalSegments;

    // Set new offset for next surface
    surfaceOffset = nextSurfaceOffset;
  }

  MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, &noSegments_[0], noSurfaces(), type_traits<MInt>::mpiType(),
                mpiComm(), AT_, "MPI_IN_PLACE", "noSegments_[0]");

  m_mpiCommSurface.resize(noSurfaces());
  m_noDomainsSurface.resize(noSurfaces());
  m_domainIdSurface.resize(noSurfaces());

  // Create MPI communicators for each surface
  for(MInt sId = 0; sId < noSurfaces(); sId++) {
    MIntScratchSpace domains(noDomains(), AT_, "domains");
    MInt noDomains_ = 0;

    // Determine list of subdomains with segments of this surface
    for(MInt d = 0; d < noDomains(); d++) {
      if(noSegments_(d, sId) > 0) {
        domains[noDomains_] = d;
        noDomains_++;
      }
    }

    m_log << "Create MPI communicator for surface #" << m_surfaceIds[sId] << " with " << noDomains_
          << " ranks:" << std::endl;
#ifndef NDEBUG
    for(MInt d = 0; d < noDomains_; d++) {
      m_log << "  * local rank #" << d << " global rank #" << domains[d] << " with " << noSegments_(domains[d], sId)
            << " segments" << std::endl;
    }
#endif

    // Obtain new group
    MPI_Group globalGroup, localGroup;
    MPI_Comm_group(mpiComm(), &globalGroup, AT_, "globalGroup");
    MPI_Group_incl(globalGroup, noDomains_, &domains[0], &localGroup, AT_);

    // Create new communicator and clean up
    MPI_Comm_create(mpiComm(), localGroup, &m_mpiCommSurface[sId], AT_,
                    "m_mpiCommSurface[" + std::to_string(sId) + "]");

    MPI_Group_free(&globalGroup, AT_);
    MPI_Group_free(&localGroup, AT_);

    if(noSurfaceElements(sId) > 0) {
      // If surface is active on current domain, determine local rank and number of domains
      MPI_Comm_rank(m_mpiCommSurface[sId], &m_domainIdSurface[sId]);
      MPI_Comm_size(m_mpiCommSurface[sId], &m_noDomainsSurface[sId]);

      TERMM_IF_NOT_COND(m_noDomainsSurface[sId] == noDomains_, "number of domains mismatch");
    } else {
      // Reset for inactive surface
      m_domainIdSurface[sId] = -1;
      m_noDomainsSurface[sId] = -1;
    }
  }

  RECORD_TIMER_STOP(m_timers[Timers::InitParallelization]);
}


/// \brief Return the name of a surface data
/// \input[in]  surfaceId id of the surface element
/// \input[in]  fileNo    index of the file
template <MInt nDim>
inline MString AcaSolver<nDim>::getSurfaceDataFileName(const MInt surfaceId, const MInt fileNo) {
  TRACE();
  const MString baseFileName = "surface_data_";
  const MString baseFileNameSuffix = ".Netcdf";
  std::ostringstream os;
  os << m_inputDir << baseFileName << std::to_string(m_surfaceIds[surfaceId]) << "_" << std::setw(8)
     << std::setfill('0') << fileNo << baseFileNameSuffix;
  return (MString)os.str();
}

/// \brief Read the surface data (or generate for an analytical case)
template <MInt nDim>
void AcaSolver<nDim>::readSurfaceData() {
  TRACE();
  printMessage("- Read surface data");
  RECORD_TIMER_START(m_timers[Timers::ReadSurfaceData]);

  // Set data sizes in surface data collector
  m_surfaceData.setNoVars(m_noVariables);
  m_surfaceData.setNoComplexVars(m_noComplexVariables);
  m_surfaceData.setNoSamples(m_noSamples);

  // Reset surface data collector to total number of local surface elements
  const MInt totalNoSurfaceElements_ = totalNoSurfaceElements();
  TERMM_IF_COND(totalNoSurfaceElements_ < 1, "Less surface elements than computational ranks is critical.");
  m_surfaceData.reset(totalNoSurfaceElements_);

  // Set storage size for times and frequencies
  m_times.resize(noSamples());
  m_frequencies.resize(noSamples());

  // Read in the surface data and store in m_surfaceData
  // or generate the surface elements and the corresponding analytical input data
  if(m_generateSurfaceData) {
    generateSurfaces();
    generateSurfaceData();

    // Testing: store generated data
    for(MInt sId = 0; sId < noSurfaces(); sId++) {
      if(hasSurface(sId) && m_storeSurfaceData) {
        storeSurfaceData(sId);
      }
    }
  } else if(m_acaPostprocessingMode) {
    // Nothing to do if just postprocessing is enabled
  } else {
    for(MInt sId = 0; sId < noSurfaces(); sId++) {
      // Read data if surface present on local rank
      if(hasSurface(sId)) {
        readSurfaceDataFromFile(sId);

        // Testing: store data for validating the I/O
        if(m_storeSurfaceData) {
          storeSurfaceData(sId);
        }
      }
    }

    { // Compare the loaded times on all ranks
      ScratchSpace<MFloat> timesBuffer(noSamples(), AT_, "timesBuffer");
      std::copy_n(&m_times[0], noSamples(), &timesBuffer[0]);

      MPI_Bcast(&timesBuffer[0], noSamples(), type_traits<MFloat>::mpiType(), 0, mpiComm(), AT_, "&timesBuffer[0]");

      for(MInt i = 0; i < noSamples(); i++) {
        TERMM_IF_COND(!approx(timesBuffer[i], m_times[i], MFloatEps), "Error: times for surfaces do not match.");
      }
    }

    // The simulation non-dimensionalises with different variables! For the FWH simulation one needs
    // to change how the numbers are non-dimensionalized
    changeDimensionsSurfaceData();
    computeDt();
  }

  RECORD_TIMER_STOP(m_timers[Timers::ReadSurfaceData]);
}


/// \brief Read the number of segments and samples for the given surface id.
template <MInt nDim>
void AcaSolver<nDim>::readNoSegmentsAndSamples(const MInt surfaceId, MInt& noSegments, MInt& noSamples,
                                               MString& inputFileName) {
  TRACE();
  using namespace maia::parallel_io;

  if(m_useMergedInputFile) {
    /*! \property
      \page propertyPageACA ACA
      \section surfaceDataFileName_
      <code>MString fileName</code>\n
      Name of the surface data files to be read in.\n
      Here, `surfaceDataFileName_0` is used for the first,
      `surfaceDataFileName_1` for the second, and so on.\n
      Keywords: <i>ACA, input</i>
    */
    const MString fileName = m_inputDir
                             + Context::getSolverProperty<MString>(
                                 "surfaceDataFileName_" + std::to_string(m_surfaceIds[surfaceId]), m_solverId, AT_);
    ParallelIo file(fileName, PIO_READ, mpiComm());

    noSegments = -1;
    noSamples = -1;

    for(auto& var : m_inputVars) {
      const MString varName = var.first;

      std::vector<MLong> dims = file.getArrayDims(varName);

      TERMM_IF_NOT_COND(noSegments < 0 || dims[0] == noSegments, "number of segments mismatch between variables");
      TERMM_IF_NOT_COND(noSamples < 0 || dims[1] == noSamples, "number of samples mismatch between variables");

      noSegments = dims[0];
      noSamples = dims[1];
    }

    // Read surface input file name if stored as attribute in file
    if(file.hasAttribute("inputFileName")) {
      file.getAttribute(&inputFileName, "inputFileName");
    }
  } else {
    noSamples = 0;
    noSegments = 0;
    const MString baseFileName = "surface_data_";
    const MString baseFileNameSuffix = ".Netcdf";
    const MString varName = "pointStates";
    for(MInt fileNo = m_inputFileIndexStart; fileNo < m_inputFileIndexEnd + 1; fileNo++) {
      const MString fileName = getSurfaceDataFileName(surfaceId, fileNo);
      ParallelIo file(fileName, PIO_READ, mpiComm());

      std::vector<MLong> dims = file.getArrayDims(varName);

      TERMM_IF_COND(noSegments != dims[0] && noSegments != 0, "number of segments mismatch between files");
      TERMM_IF_COND(dims[2] < (MInt)m_inputVars.size(), "number of input vars > number of variables in file");

      noSegments = dims[0];
      noSamples += dims[1];

      // TODO labels:ACA sanity check: if variables in input file are the same as requested

      // Read surface input file name if stored as attribute in file
      if(file.hasAttribute("inputFileName", varName)) {
        file.getAttribute(&inputFileName, "inputFileName", varName);
      }
    }
  }
}


/// \brief Read the data for one surface.
template <MInt nDim>
void AcaSolver<nDim>::readSurfaceDataFromFile(const MInt surfaceId) {
  TRACE();
  using namespace maia::parallel_io;

  const MInt count = m_noSurfaceElements[surfaceId];
  const MInt offset = m_surfaceElementOffsets[surfaceId];
  const MInt surfaceOffset = localSurfaceOffset(surfaceId);

  const MString measureName = (nDim == 3) ? "area" : "length";

  // Create surface elements and store its data in m_surfaceData
  auto createSurfaceElements = [&](const MInt l_surfaceId,
                                   const MFloatScratchSpace& l_measure,
                                   const MFloatScratchSpace& l_coordinates,
                                   const MFloatScratchSpace& l_normalVector) {
    MFloat totalSurfaceMeasure = 0.0;
    MFloat totalSurfaceMeasureWeighted = 0.0;

    const MFloat weight = m_surfaceWeightingFactor[l_surfaceId];
    for(MInt i = 0; i < count; i++) {
      const MInt id = surfaceOffset + i;
      // Add a new surface element
      m_surfaceData.append();
      // Check that the total size of m_surfaceData is correct
      TERMM_IF_NOT_COND(m_surfaceData.size() == id + 1, "m_surfaceData size mismatch.");

      MFloat* surfCoordinates = &m_surfaceData.surfaceCoordinates(id);
      MFloat* surfNormal = &m_surfaceData.surfaceNormal(id);

      // Store segment area/length (weighted with surface averaging factor)
      m_surfaceData.surfaceArea(id) = weight * l_measure[i];
      totalSurfaceMeasure += l_measure[i];
      totalSurfaceMeasureWeighted += (weight * l_measure[i]);

      // Set coordinates and normal vector
      for(MInt j = 0; j < nDim; j++) {
        surfCoordinates[j] = l_coordinates[i * nDim + j];
        surfNormal[j] = l_normalVector[i * nDim + j];
      }
    }

    // Compute total area/length of surface
    MPI_Allreduce(MPI_IN_PLACE, &totalSurfaceMeasure, 1, type_traits<MFloat>::mpiType(), MPI_SUM,
                  m_mpiCommSurface[surfaceId], AT_, "MPI_IN_PLACE", "totalSurfaceMeasure");
    MPI_Allreduce(MPI_IN_PLACE, &totalSurfaceMeasureWeighted, 1, type_traits<MFloat>::mpiType(), MPI_SUM,
                  m_mpiCommSurface[surfaceId], AT_, "MPI_IN_PLACE", "totalSurfaceMeasureWeighted");
    if(m_domainIdSurface[surfaceId] == 0) {
      std::cout << "Total " << measureName << " of surface #" << surfaceId << ": " << totalSurfaceMeasure
                << " (weighted: " << totalSurfaceMeasureWeighted << ")" << std::endl;
    }
  };

  if(m_useMergedInputFile) {
    const MString fileName = m_inputDir
                             + Context::getSolverProperty<MString>(
                                 "surfaceDataFileName_" + std::to_string(m_surfaceIds[surfaceId]), m_solverId, AT_);
    ParallelIo file(fileName, PIO_READ, m_mpiCommSurface[surfaceId]);

    MFloatScratchSpace coordinates(count * nDim, AT_, "coordinates");
    MFloatScratchSpace normalVector(count * nDim, AT_, "normalVector");
    MFloatScratchSpace measure(count, AT_, "measure");

    // Read surface area/length
    file.setOffset(count, offset);
    file.readArray(&measure[0], measureName);

    // Read coordinates and normal vectors
    file.setOffset(count, offset, 2);
    file.readArray(&coordinates[0], "coordinates");
    file.readArray(&normalVector[0], "normalVector");

    // Create surface elements and store its data in m_surfaceData
    createSurfaceElements(surfaceId, measure, coordinates, normalVector);

    // Read times only on surface local rank 0 and broadcast to other ranks in this communicator
    ScratchSpace<MFloat> timesBuffer(m_totalNoSamples, AT_, "timesBuffer");
    if(m_domainIdSurface[surfaceId] == 0) {
      // Read full array, then copy selected range to m_times
      file.setOffset(m_totalNoSamples, 0);
      file.readArray(&timesBuffer[0], "time");

      // Check the times with those of other surfaces already loaded on this rank (before overwriting them), later
      // cross-check the times with all surfaces
      if(m_hasTimes) {
        for(MInt i = 0; i < noSamples(); i++) {
          TERMM_IF_COND(!approx(timesBuffer[m_sampleOffset + i * m_sampleStride], m_times[i], MFloatEps),
                        "Error: local times for surfaces do not match.");
        }
      }

      // Store times in m_times
      if(m_sampleStride == 1) {
        std::copy_n(&timesBuffer[m_sampleOffset], noSamples(), &m_times[0]);
      } else {
        for(MInt i = 0; i < noSamples(); i++) {
          m_times[i] = timesBuffer[m_sampleOffset + i * m_sampleStride];
        }
      }
    } else {
      file.setOffset(0, 0);
      file.readArray(&timesBuffer[0], "time");
    }

    MPI_Bcast(&m_times[0], noSamples(), type_traits<MFloat>::mpiType(), 0, m_mpiCommSurface[surfaceId], AT_,
              "&m_times[0]");
    m_hasTimes = true;

    // Set 2D offset (second dimension is the time axis which is set internally to read all samples)
    // TODO labels:ACA,IO directly read only selected sample range once this is supported by the ParallelIO class
    file.setOffset(count, offset, 2);

    // Data buffer
    ScratchSpace<MFloat> varData(count, m_totalNoSamples, AT_, "varData");
    // Read datasets and store in m_surfaceData
    for(auto& var : m_inputVars) {
      const MString varName = var.first;

      if(m_domainIdSurface[surfaceId] == 0) {
        std::cout << "Surface #" << m_surfaceIds[surfaceId] << ": Loading variable: " << varName << std::endl;
      }

      file.readArray(&varData[0], varName);

      const MInt varIndex = var.second;
      for(MInt i = 0; i < count; i++) {
        const MInt id = surfaceOffset + i;
        for(MInt t = 0; t < noSamples(); t++) {
          m_surfaceData.variables(id, varIndex, t) = varData(i, m_sampleOffset + t * m_sampleStride);
        }
      }
    }
  } else { // ! m_useMergedInputFile -------------------------------------------
    TERMM_IF_COND(m_sampleOffset > 0, "Sample offset for multiple input files not supported yet");
    TERMM_IF_COND(m_sampleStride > 1, "Sample stride for multiple input files not supported yet");

    const MInt noDom = m_noDomainsSurface[surfaceId];
    std::vector<MInt> sendCount(noDom);
    std::vector<MInt> sendCountNdim(noDom);
    std::vector<MInt> displs(noDom);
    std::vector<MInt> displsNdim(noDom);
    {
      std::fill(sendCount.begin(), sendCount.end(), 0);
      std::fill(displs.begin(), displs.end(), 0);
      sendCount[m_domainIdSurface[surfaceId]] = count;
      displs[m_domainIdSurface[surfaceId]] = offset;
      MPI_Allreduce(MPI_IN_PLACE, sendCount.data(), sendCount.size(), MPI_INT, MPI_SUM, m_mpiCommSurface[surfaceId],
                    AT_, "MPI_IN_PLACE", "sendCount");
      MPI_Allreduce(MPI_IN_PLACE, displs.data(), displs.size(), MPI_INT, MPI_SUM, m_mpiCommSurface[surfaceId], AT_,
                    "MPI_IN_PLACE", "displs");
      for(MInt i = 0; i < (MInt)sendCount.size(); i++) {
        sendCountNdim[i] = sendCount[i] * nDim;
        displsNdim[i] = displs[i] * nDim;
      }
    }
    const MInt noTotalCount = std::accumulate(sendCount.begin(), sendCount.end(), 0);

    //--Read sortIndex from first input file
    std::vector<MInt> sortIndex;
    {
      auto readSortIndex = [&](const MInt fileNo, std::vector<MInt>& sortIndex_) {
        //-Open the correct sample file
        const MString fileName = getSurfaceDataFileName(surfaceId, fileNo);
        ParallelIo file(fileName, PIO_READ, m_mpiCommSurface[surfaceId]);
        //-First read the famous sortIndex_ from first file
        const MString varName = "sortIndex";
        if(m_domainIdSurface[surfaceId] == 0) {
          sortIndex_.resize(noTotalCount);
          file.setOffset(noTotalCount, 0);
          file.readArray(sortIndex_.data(), varName);
        } else {
          sortIndex_.resize(0);
          file.setOffset(0, 0);
          file.readArray(sortIndex_.data(), varName);
        }
      };
      readSortIndex(m_inputFileIndexStart, sortIndex);
      // Sanity check: Is sort index the same for all input files?
      std::vector<MInt> tmpSortIndex(sortIndex.size());
      for(MInt fileNo = 1 + m_inputFileIndexStart; fileNo < m_inputFileIndexEnd + 1; fileNo++) {
        readSortIndex(fileNo, tmpSortIndex);
        if(!std::equal(sortIndex.begin(), sortIndex.end(), tmpSortIndex.begin())) {
          std::stringstream ss;
          ss << "SortIndex varies between ";
          ss << getSurfaceDataFileName(surfaceId, m_inputFileIndexStart);
          ss << " and ";
          ss << getSurfaceDataFileName(surfaceId, fileNo);
          TERMM(1, ss.str());
        }
      }
    }

    //--Read geometry data from separate file
    {
      MFloatScratchSpace measure(count, AT_, "measure");
      MFloatScratchSpace coordinates(count * nDim, AT_, "coordinates");
      MFloatScratchSpace normalVector(count * nDim, AT_, "normalVector");

      //-Read everything from first rank and sort it
      const MString varName = "geomMeasure";
      const MString baseFileNameSuffix = ".Netcdf";
      const MString baseGeometryFileName = "surface_geometryInfo_";
      std::ostringstream fileName;
      fileName << m_inputDir << baseGeometryFileName << std::to_string(m_surfaceIds[surfaceId]) << baseFileNameSuffix;
      ParallelIo file(fileName.str(), PIO_READ, m_mpiCommSurface[surfaceId]);
      if(m_domainIdSurface[surfaceId] == 0) {
        MFloatScratchSpace sendBuffer(noTotalCount * nDim, AT_, "sendBuffer");
        MFloatScratchSpace readBuffer(noTotalCount * nDim, AT_, "readBuffer");
        // Read surface area/length
        file.setOffset(noTotalCount, 0);
        file.readArray(readBuffer.data(), varName);
        for(MInt i = 0; i < noTotalCount; i++) {
          sendBuffer(i) = readBuffer(sortIndex[i]);
        }
        MPI_Scatterv(sendBuffer.data(), sendCount.data(), displs.data(), MPI_DOUBLE, measure.data(), count, MPI_DOUBLE,
                     0, m_mpiCommSurface[surfaceId], AT_, "send", "recv");

        // Read coordinates and normal vectors
        file.setOffset(noTotalCount * nDim, 0);
        file.readArray(readBuffer.data(), "centroid");
        for(MInt i = 0; i < noTotalCount; i++) {
          for(MInt j = 0; j < nDim; j++) {
            sendBuffer(i * nDim + j) = readBuffer(sortIndex[i] * nDim + j);
          }
        }
        MPI_Scatterv(sendBuffer.data(), sendCountNdim.data(), displsNdim.data(), MPI_DOUBLE, coordinates.data(),
                     count * nDim, MPI_DOUBLE, 0, m_mpiCommSurface[surfaceId], AT_, "send", "recv");
        file.readArray(readBuffer.data(), "normalVector");
        for(MInt i = 0; i < noTotalCount; i++) {
          for(MInt j = 0; j < nDim; j++) {
            sendBuffer(i * nDim + j) = readBuffer(sortIndex[i] * nDim + j);
          }
        }
        MPI_Scatterv(sendBuffer.data(), sendCountNdim.data(), displsNdim.data(), MPI_DOUBLE, normalVector.data(),
                     count * nDim, MPI_DOUBLE, 0, m_mpiCommSurface[surfaceId], AT_, "send", "recv");
        //-Scatter data (send)
      } else { // not root rank -> receive sorted data
        //-Scatter data (receive)
        MInt* dummySendBuffer = nullptr;
        file.setOffset(0, 0);
        file.readArray(dummySendBuffer, varName);
        MPI_Scatterv(dummySendBuffer, sendCount.data(), displs.data(), MPI_DOUBLE, measure.data(), count, MPI_DOUBLE, 0,
                     m_mpiCommSurface[surfaceId], AT_, "send", "recv");
        file.readArray(dummySendBuffer, "centroid");
        MPI_Scatterv(dummySendBuffer, sendCountNdim.data(), displsNdim.data(), MPI_DOUBLE, coordinates.data(),
                     count * nDim, MPI_DOUBLE, 0, m_mpiCommSurface[surfaceId], AT_, "send", "recv");
        file.readArray(dummySendBuffer, "normalVector");
        MPI_Scatterv(dummySendBuffer, sendCountNdim.data(), displsNdim.data(), MPI_DOUBLE, normalVector.data(),
                     count * nDim, MPI_DOUBLE, 0, m_mpiCommSurface[surfaceId], AT_, "send", "recv");
      }

      // Create surface elements and store its data in m_surfaceData
      createSurfaceElements(surfaceId, measure, coordinates, normalVector);
    }

    //--Read times only on surface local rank 0 and broadcast to other ranks in this communicator
    {
      const MString varName = "time";
      // const MString varName = "timeStep";
      ScratchSpace<MFloat> timesBuffer(m_totalNoSamples, AT_, "timesBuffer");
      if(m_domainIdSurface[surfaceId] == 0) {
        MInt offsetInBuffer = 0;
        for(MInt fileNo = m_inputFileIndexStart; fileNo < m_inputFileIndexEnd + 1; fileNo++) {
          // get fileName of fileNo th sample
          const MString fileName = getSurfaceDataFileName(surfaceId, fileNo);
          ParallelIo file(fileName, PIO_READ, m_mpiCommSurface[surfaceId]);

          std::vector<MLong> dims = file.getArrayDims(varName);
          const MInt samplesInFile = dims[0];

          file.setOffset(samplesInFile, 0);
          file.readArray(&timesBuffer[offsetInBuffer], varName);

          offsetInBuffer += samplesInFile;
        }
        // Check the times with those of other surfaces already loaded on this rank (before overwriting them), later
        // cross-check the times with all surfaces
        if(m_hasTimes) {
          for(MInt i = 0; i < noSamples(); i++) {
            TERMM_IF_COND(!approx(timesBuffer[m_sampleOffset + i * m_sampleStride], m_times[i], MFloatEps),
                          "Error: local times for surfaces do not match.");
          }
        }

        // Store times in m_times
        if(m_sampleStride == 1) {
          std::copy_n(&timesBuffer[m_sampleOffset], noSamples(), &m_times[0]);
        } else {
          for(MInt i = 0; i < noSamples(); i++) {
            m_times[i] = timesBuffer[m_sampleOffset + i * m_sampleStride];
          }
        }
      } else { // not root process
        for(MInt fileNo = m_inputFileIndexStart; fileNo < m_inputFileIndexEnd + 1; fileNo++) {
          // get fileName of fileNo th sample
          const MString fileName = getSurfaceDataFileName(surfaceId, fileNo);
          // read only dummy to ensure valid 'parallel' reading
          ParallelIo file(fileName, PIO_READ, m_mpiCommSurface[surfaceId]);
          file.setOffset(0, 0);
          file.readArray(&timesBuffer[0], varName);
        }
      }
      MPI_Bcast(&m_times[0], noSamples(), type_traits<MFloat>::mpiType(), 0, m_mpiCommSurface[surfaceId], AT_,
                "&m_times[0]");
      m_hasTimes = true;
    }

    //--Read data from each file and sort it correctly into my buffers
    {
      const MInt noVars = m_inputVars.size();
      MInt offsetInBuffer = 0;
      for(MInt fileNo = m_inputFileIndexStart; fileNo < m_inputFileIndexEnd + 1; fileNo++) {
        //-Open the correct sample file
        const MString fileName = getSurfaceDataFileName(surfaceId, fileNo);
        ParallelIo file(fileName, PIO_READ, m_mpiCommSurface[surfaceId]);

        //-Read variable data
        {
          // Read data
          const MString varName = "pointStates";
          file.setOffset(count, offset, 3);

          std::vector<MLong> dims = file.getArrayDims(varName);
          const MInt samplesInFile = dims[1];

          ScratchSpace<MFloat> varData(count, samplesInFile, noVars, AT_, "varData");
          file.readArray(varData.data(), varName);

          // Create mapping from input variables in file to order of stored variable
          std::vector<MInt> varOrder;
          {
            std::map<MInt, MString> fileVars{};
            for(MInt varIndex = 0; varIndex < noVars; varIndex++) {
              MString var;
              file.getAttribute(&var, "var_" + std::to_string(varIndex), varName);
              fileVars[varIndex] = var;
            }
            for(MInt varIndex = 0; varIndex < noVars; varIndex++) {
              varOrder.push_back(m_inputVars[fileVars[varIndex]]);
            }
          }

          // Store in local data
          for(MInt i = 0; i < count; i++) {
            const MInt id = surfaceOffset + i;
            for(MInt t = 0; t < samplesInFile && offsetInBuffer + t < noSamples(); t++) {
              for(MInt varIndex = 0; varIndex < noVars; varIndex++) {
                const MInt var = varOrder[varIndex];
                // TODO labels:ACA  add possibility to set a strideA
                // TODO labels:ACA add offset, e.g. skip first X samples in files
                m_surfaceData.variables(id, var, offsetInBuffer + t) = varData(id, t, varIndex);
              }
            }
          }
          offsetInBuffer += samplesInFile;
        }
      }
    }
  }
}


/// \brief Store the surface data.
/// Note: used for validating the data input and for the output of generated data.
template <MInt nDim>
void AcaSolver<nDim>::storeSurfaceData(const MInt surfaceId) {
  TRACE();
  using namespace maia::parallel_io;
  const MString fileName = outputDir() + m_outputFilePrefix + "out_surfaceData_"
                           + std::to_string(m_surfaceIds[surfaceId]) + ParallelIo::fileExt();
  std::stringstream ss;
  ss << "  "
     << "Testing: Storing surface data in " << fileName;
  printMessage(ss.str());

  ParallelIo file(fileName, PIO_REPLACE, m_mpiCommSurface[surfaceId]);

  const MInt totalCount = m_noSurfaceElementsGlobal[surfaceId];
  const MInt count = m_noSurfaceElements[surfaceId];
  const MInt offset = m_surfaceElementOffsets[surfaceId];

  const MString measureName = (nDim == 3) ? "area" : "length";
  file.defineArray(PIO_FLOAT, measureName, totalCount);

  file.defineArray(PIO_FLOAT, "time", noSamples());

  ParallelIo::size_type dimSizesXD[] = {totalCount, nDim};

  file.defineArray(PIO_FLOAT, "coordinates", 2, &dimSizesXD[0]);
  file.defineArray(PIO_FLOAT, "normalVector", 2, &dimSizesXD[0]);

  ParallelIo::size_type dimSizesVar[] = {totalCount, noSamples()};

  for(auto& var : m_inputVars) {
    const MString varName = var.first;
    file.defineArray(PIO_FLOAT, varName, 2, &dimSizesVar[0]);
  }

  const MInt surfaceOffset = localSurfaceOffset(surfaceId);

  // Write surface area/length
  file.setOffset(count, offset);
  file.writeArray(&m_surfaceData.surfaceArea(surfaceOffset), measureName);

  // Write surface area/length
  file.setOffset(((m_domainIdSurface[surfaceId] == 0) ? noSamples() : 0), 0);
  file.writeArray(&m_times[0], "time");

  // Write coordinates and normal vectors
  file.setOffset(count, offset, 2);
  file.writeArray(&m_surfaceData.surfaceCoordinates(surfaceOffset), "coordinates");
  file.writeArray(&m_surfaceData.surfaceNormal(surfaceOffset), "normalVector");

  file.setOffset(count, offset, 2);
  ScratchSpace<MFloat> varBuffer(count, noSamples(), AT_, "varBuffer");
  for(auto& var : m_inputVars) {
    const MInt varIndex = var.second;
    for(MInt i = 0; i < count; i++) {
      for(MInt t = 0; t < noSamples(); t++) {
        varBuffer(i, t) = m_surfaceData.variables(surfaceOffset + i, varIndex, t);
      }
    }

    file.writeArray(&varBuffer[0], var.first);
  }
}


/// \brief Generation of surfaces (for analytical cases).
template <MInt nDim>
void AcaSolver<nDim>::generateSurfaces() {
  TRACE();

  // Create all surfaces according to their specified type
  for(MInt sId = 0; sId < noSurfaces(); sId++) {
    /*! \property
      \page propertyPageACA ACA
      \section surfaceType_
      <code>MInt surfaceType</code>\n\n
      Type of the generated surfaces.\n
      Here, `surfaceType_0` is used for the first, `surfaceType_1` for the
      second, and so on.\n
      Possible values are: \n
      <ul>
      <li><code>0</code> plane </li>
      <li><code>1</code> circle </li>
      </ul>
      Keywords: <i>ACA, input</i>
    */
    const MInt surfaceType =
        Context::getSolverProperty<MInt>("surfaceType_" + std::to_string(m_surfaceIds[sId]), solverId(), AT_);

    m_surfaceInputFileName[sId] = (surfaceType == 0) ? "generated plane" : "generated circle";

    // Nothing to do if the surface is not present
    if(!hasSurface(sId)) {
      continue;
    }

    switch(surfaceType) {
      case 0:
        generateSurfacePlane(sId);
        break;
      case 1:
        generateSurfaceCircle(sId);
        break;
      default:
        TERMM(1, "Unknown surface type: " + std::to_string(surfaceType) + "; possible values: 0 - plane.");
        break;
    }
  }
}


/// \brief Determine the number of segments to be generated for a surface, returns the total number
///        of segments and a vector with the number of segments in each direction for the specific
///        surface type
template <MInt nDim>
MInt AcaSolver<nDim>::readNoSegmentsSurfaceGeneration(const MInt sId, std::vector<MInt>& noSegments,
                                                      const MInt lengthCheck) {
  /*! \property
    \page propertyPageACA ACA
    \section noSegments_
    <code>MInt noSegments</code>\n\n
    Number of surface segments to be generated for a surface file.\n
    Here, `noSegments_0` is used for the first, `noSegments_1` for the
    second, and so on.\n
    Keywords: <i>ACA, input</i>
  */
  const MString noSegmentsName = "noSegments_" + std::to_string(m_surfaceIds[sId]);

  // If a length is given, check if it matches the property length
  if(lengthCheck > 0) {
    Context::assertPropertyLength(noSegmentsName, lengthCheck, solverId());
  }

  MInt totalNoSegments = 1;
  const MInt propLength = Context::propertyLength(noSegmentsName);
  noSegments.resize(propLength);

  // Read all property values
  for(MInt i = 0; i < propLength; i++) {
    noSegments[i] = Context::getSolverProperty<MInt>(noSegmentsName, m_solverId, AT_, i);
    TERMM_IF_NOT_COND(noSegments[i] >= 0, "Number of segments needs to be >= 0");

    // Determine total number of segments by product of number of segments for each direction
    if(noSegments[i] > 0) {
      totalNoSegments *= noSegments[i];
    }
  }

  TERMM_IF_NOT_COND(totalNoSegments > 0, "Total number of segments needs to be > 0");
  return totalNoSegments;
}


/// \brief Generate a surface plane (Cartesian).
template <MInt nDim>
void AcaSolver<nDim>::generateSurfacePlane(const MInt sId) {
  TRACE();

  // Weighting factor for surface averaging
  const MFloat weight = m_surfaceWeightingFactor[sId];

  if constexpr(nDim == 3) {
    // Determine number of segments for each Cartesian direction
    std::vector<MInt> noSegmentsPerPlane{};
    // using nDim because there are all 3 directions in the properties, e.g. [0, 4, 4].
    MInt totalNoSeg = readNoSegmentsSurfaceGeneration(sId, noSegmentsPerPlane, nDim);

    /*! \property
      \page propertyPageACA ACA
      \section surfaceNormalDir_
      <code>MInt surfaceNormalDir</code>\n\n
      Get surface outer normal direction (0:-x; 1:x; 2:-y; 3:y; 4:-z; 5:z).\n
      Here, `surfaceNormalDir_0` is used for the first, `surfaceNormalDir_1` for
      the second, and so on.\n
      Keywords: <i>ACA, input</i>
    */
    const MInt surfaceNormalDir =
        Context::getSolverProperty<MInt>("surfaceNormalDir_" + std::to_string(m_surfaceIds[sId]), m_solverId, AT_);
    TERMM_IF_NOT_COND(surfaceNormalDir >= 0 && surfaceNormalDir < 2 * nDim,
                      "Invalid surface outer normal direction: " + std::to_string(surfaceNormalDir));

    // Determine surface normal vector (axis aligned)
    std::vector<MFloat> surfaceNormal(nDim, 0.0);
    for(MInt i = 0; i < nDim; i++) {
      if(surfaceNormalDir == 2 * i) { // 0 == -x; 2 == -y; 4 == -z
        surfaceNormal[i] = -1.0;
      } else if(surfaceNormalDir == 2 * i + 1) { // 1 == x; 3 == y; 5 == z
        surfaceNormal[i] = 1.0;
      } else {
        surfaceNormal[i] = 0.0;
      }
    }

    // Three points define the plane in 3D, make sure it is axis aligned or fix surface normal
    // vector
    std::vector<MFloat> planeCoordinates(2 * nDim, 0.0);

    const MInt noCoords = 2 * nDim;
    /*! \property
      \page propertyPageACA ACA
      \section surfaceCoordinates_
      <code>std::vector<MFloat> planeCoordinates(2*nDim)</code>\n\n
      default = <code>[0.0, 0.0, 0.0]</code>\n\n
      Number of surface segments to be generated for a surface file.\n
      Here, `surfaceCoordinates_0` is used for the first, `surfaceCoordinates_1` for the
      second, and so on.\n
      Keywords: <i>ACA, input</i>
    */
    const MString coordName = "surfaceCoordinates_" + std::to_string(m_surfaceIds[sId]);
    Context::assertPropertyLength(coordName, noCoords, solverId());

    // Read coordinates
    for(MInt i = 0; i < noCoords; i++) {
      planeCoordinates[i] = Context::getSolverProperty<MFloat>(coordName, m_solverId, AT_, i);
    }

    // Get length of each side
    const MFloat length_x = std::abs(planeCoordinates[nDim] - planeCoordinates[0]);
    const MFloat length_y = std::abs(planeCoordinates[nDim + 1] - planeCoordinates[1]);
    const MFloat length_z = std::abs(planeCoordinates[nDim + 2] - planeCoordinates[2]);

    MFloat start_x, start_y, start_z, noSeg_x, noSeg_y, noSeg_z, stepSize_x, stepSize_y, stepSize_z;
    // Determine coordinates, number of Segments per direction and stepSize per direction
    if(approx(length_x, 0.0, MFloatEps)) {
      start_x = planeCoordinates[0];
      noSeg_x = 0;
      stepSize_x = 0;
    } else {
      noSeg_x =
          Context::getSolverProperty<MFloat>("noSegments_" + std::to_string(m_surfaceIds[sId]), m_solverId, AT_, 0);
      stepSize_x = length_x / noSeg_x;
      start_x = planeCoordinates[0] + stepSize_x / 2.0;
    }

    if(approx(length_y, 0.0, MFloatEps)) {
      start_y = planeCoordinates[1];
      noSeg_y = 0;
      stepSize_y = 0;
    } else {
      noSeg_y =
          Context::getSolverProperty<MFloat>("noSegments_" + std::to_string(m_surfaceIds[sId]), m_solverId, AT_, 1);
      stepSize_y = length_y / noSeg_y;
      start_y = planeCoordinates[1] + stepSize_y / 2.0;
    }

    if(approx(length_z, 0.0, MFloatEps)) {
      start_z = planeCoordinates[2];
      noSeg_z = 0;
      stepSize_z = 0;
    } else {
      noSeg_z =
          Context::getSolverProperty<MFloat>("noSegments_" + std::to_string(m_surfaceIds[sId]), m_solverId, AT_, 2);
      stepSize_z = length_z / noSeg_z;
      start_z = planeCoordinates[2] + stepSize_z / 2.0;
    }

    // Check for axis aligned surface
    TERMM_IF_NOT_COND(
        (approx(length_x, 0.0, MFloatEps) || approx(length_y, 0.0, MFloatEps) || approx(length_z, 0.0, MFloatEps)),
        "Plane not axis aligned.");

    std::vector<MFloat> coordinates_buff(nDim * totalNoSeg, 0.0);
    MFloat delta_segment = 0.0;
    // Determine segment area and save coordinates in data buffer
    if(approx(stepSize_x, 0.0, MFloatEps)) {
      delta_segment = stepSize_y * stepSize_z;

      for(MInt j = 0; j < noSeg_z; j++) {
        for(MInt i = 0; i < noSeg_y; i++) {
          coordinates_buff[j * 3 * noSeg_y + i * 3] = start_x;
          coordinates_buff[j * 3 * noSeg_y + i * 3 + 1] = start_y + i * stepSize_y;
          coordinates_buff[j * 3 * noSeg_y + i * 3 + 2] = start_z + j * stepSize_z;
        }
      }
    } else if(approx(stepSize_y, 0.0, MFloatEps)) {
      delta_segment = stepSize_x * stepSize_z;

      for(MInt j = 0; j < noSeg_z; j++) {
        for(MInt i = 0; i < noSeg_x; i++) {
          coordinates_buff[j * 3 * noSeg_x + i * 3] = start_x + i * stepSize_x;
          coordinates_buff[j * 3 * noSeg_x + i * 3 + 1] = start_y;
          coordinates_buff[j * 3 * noSeg_x + i * 3 + 2] = start_z + j * stepSize_z;
        }
      }

    } else if(approx(stepSize_z, 0.0, MFloatEps)) {
      delta_segment = stepSize_x * stepSize_y;

      for(MInt j = 0; j < noSeg_x; j++) {
        for(MInt i = 0; i < noSeg_y; i++) {
          coordinates_buff[j * 3 * noSeg_y + i * 3] = start_x + j * stepSize_x;
          coordinates_buff[j * 3 * noSeg_y + i * 3 + 1] = start_y + i * stepSize_y;
          coordinates_buff[j * 3 * noSeg_y + i * 3 + 2] = start_z;
        }
      }
    }

    // Surface Offset in the case that one processor needs to compute points on different surfaces.
    const MInt surfaceOffset = localSurfaceOffset(sId);

    // Loop over all local elements
    for(MInt i = 0; i < m_noSurfaceElements[sId]; i++) {
      const MInt id = surfaceOffset + i;
      // Add a new surface element
      m_surfaceData.append();
      // Check that the total size of m_surfaceData is correct
      TERMM_IF_NOT_COND(m_surfaceData.size() == id + 1, "m_surfaceData size mismatch.");

      // Id of segment on the surface. In case processor starts computing points somewhere on the
      // surface.
      const MInt segmentId = i + m_surfaceElementOffsets[sId];
      // Check for valid segment id
      TERMM_IF_NOT_COND(segmentId >= 0 && segmentId < totalNoSeg, "invalid segment id");

      // Read out coordinates from buffer
      const MFloat xcenter = coordinates_buff[3 * segmentId];
      const MFloat ycenter = coordinates_buff[3 * segmentId + 1];
      const MFloat zcenter = coordinates_buff[3 * segmentId + 2];

      // Store segment center coordinates
      MFloat* surfCoordinates = &m_surfaceData.surfaceCoordinates(id);
      surfCoordinates[0] = xcenter;
      surfCoordinates[1] = ycenter;
      surfCoordinates[2] = zcenter;

      // Store segment area (weighted with surface averaging factor)
      m_surfaceData.surfaceArea(id) = weight * delta_segment;

      // Store segment normal vector
      MFloat* surfNormal = &m_surfaceData.surfaceNormal(id);
      surfNormal[0] = surfaceNormal[0];
      surfNormal[1] = surfaceNormal[1];
      surfNormal[2] = surfaceNormal[2];
    }
  } else if constexpr(nDim == 2) {
    // Determine number of segments for each Cartesian direction
    std::vector<MInt> noSegmentsPerDir{};
    readNoSegmentsSurfaceGeneration(sId, noSegmentsPerDir, nDim - 1);

    // Get surface outer normal direction (0:-x; 1:x; 2:-y; 3:y; 4:-z; 5:z)
    const MInt surfaceNormalDir =
        Context::getSolverProperty<MInt>("surfaceNormalDir_" + std::to_string(m_surfaceIds[sId]), m_solverId, AT_);
    TERMM_IF_NOT_COND(surfaceNormalDir >= 0 && surfaceNormalDir < 2 * nDim,
                      "Invalid surface outer normal direction: " + std::to_string(surfaceNormalDir));

    // Determine surface normal vector (axis aligned)
    std::vector<MFloat> surfaceNormal(nDim, 0.0);
    for(MInt i = 0; i < nDim; i++) {
      if(surfaceNormalDir == 2 * i) { // 0 == -x; 2 == -y; 4 == -z
        surfaceNormal[i] = -1.0;
      } else if(surfaceNormalDir == 2 * i + 1) { // 1 == x; 3 == y; 5 == z
        surfaceNormal[i] = 1.0;
      } else {
        surfaceNormal[i] = 0.0;
      }
    }

    // Two points define the plane in 2D, make sure it is axis aligned or fix surface normal vector
    std::vector<MFloat> planeCoordinates(2 * nDim, 0.0);

    const MInt noCoords = 2 * nDim;
    const MString coordName = "surfaceCoordinates_" + std::to_string(m_surfaceIds[sId]);
    Context::assertPropertyLength(coordName, noCoords, solverId());

    // Read coordinates
    for(MInt i = 0; i < noCoords; i++) {
      planeCoordinates[i] = Context::getSolverProperty<MFloat>(coordName, m_solverId, AT_, i);
    }

    const MFloat noSegmentsInv = 1.0 / (noSegmentsPerDir[0]);
    const MFloat Delta_x = (planeCoordinates[nDim] - planeCoordinates[0]);
    const MFloat Delta_y = (planeCoordinates[nDim + 1] - planeCoordinates[1]);

    // Check for axis aligned surface
    TERMM_IF_NOT_COND(approx(Delta_x, 0.0, MFloatEps) || approx(Delta_y, 0.0, MFloatEps), "Plane not axis aligned.");

    const MInt surfaceOffset = localSurfaceOffset(sId);

    // Loop over local number of surface elememts
    for(MInt i = 0; i < m_noSurfaceElements[sId]; i++) {
      const MInt id = surfaceOffset + i;
      // Add a new surface element
      m_surfaceData.append();
      // Check that the total size of m_surfaceData is correct
      TERMM_IF_NOT_COND(m_surfaceData.size() == id + 1, "m_surfaceData size mismatch.");

      // Id of segment on the surface
      const MInt segmentId = i + m_surfaceElementOffsets[sId];
      // Check for valid segment id
      TERMM_IF_NOT_COND(segmentId >= 0 && segmentId < noSegmentsPerDir[0], "invalid segment id");

      // First segment point
      const MFloat p1_x = planeCoordinates[0] + segmentId * noSegmentsInv * Delta_x;
      const MFloat p1_y = planeCoordinates[1] + segmentId * noSegmentsInv * Delta_y;

      // Second segment point
      const MFloat p2_x = planeCoordinates[0] + (segmentId + 1) * noSegmentsInv * Delta_x;
      const MFloat p2_y = planeCoordinates[1] + (segmentId + 1) * noSegmentsInv * Delta_y;

      // Segment center
      const MFloat xcenter = 0.5 * (p2_x + p1_x);
      const MFloat ycenter = 0.5 * (p2_y + p1_y);

      // Determine segment length
      const MFloat dx = (p2_x - p1_x);
      const MFloat dy = (p2_y - p1_y);
      const MFloat delta_segment = sqrt(dx * dx + dy * dy);

      // Store segment center coordinates
      MFloat* surfCoordinates = &m_surfaceData.surfaceCoordinates(id);
      surfCoordinates[0] = xcenter;
      surfCoordinates[1] = ycenter;

      // Store segment area (weighted with surface averaging factor)
      m_surfaceData.surfaceArea(id) = weight * delta_segment;

      // Store segment normal vector
      MFloat* surfNormal = &m_surfaceData.surfaceNormal(id);
      surfNormal[0] = surfaceNormal[0];
      surfNormal[1] = surfaceNormal[1];
    }
  }
}


/// \brief Generate a circular surface in 2D.
template <MInt nDim>
void AcaSolver<nDim>::generateSurfaceCircle(const MInt sId) {
  TRACE();

  // Weighting factor for surface averaging
  const MFloat weight = m_surfaceWeightingFactor[sId];

  if constexpr(nDim != 2) {
    TERMM(1, "Only implemented in 2D.");
  }

  // TODO labels:ACA Fix this
  TERMM(1, "The 2D implementation is not tested and seems not to work properly.");

  // Creating Circle with radius and number of Elements
  MInt totalNoSurfaceElements = 100;
  totalNoSurfaceElements = Context::getSolverProperty<MInt>("noSegments_" + std::to_string(m_surfaceIds[sId]),
                                                            m_solverId, AT_, &totalNoSurfaceElements);

  MFloat radius = 1.0;
  radius = Context::getSolverProperty<MFloat>("radius", m_solverId, AT_, &radius);

  const MFloat step_rad = 2 * PI / static_cast<MFloat>(totalNoSurfaceElements);

  std::vector<MFloat> surfaceNormal(2, 0.0);
  const MInt surfaceOffset = localSurfaceOffset(sId);
  for(MInt i = 0; i < m_noSurfaceElements[sId]; i++) {
    const MInt id = surfaceOffset + i;
    // Add a new surface element
    m_surfaceData.append();
    // Check that the total size of m_surfaceData is correct
    TERMM_IF_NOT_COND(m_surfaceData.size() == id + 1, "m_surfaceData size mismatch.");

    // First surface point
    const MFloat p1_x = radius * cos(i * step_rad);
    const MFloat p1_y = radius * sin(i * step_rad);

    // Second surface point
    const MFloat p2_x = radius * cos((i + 1) * step_rad);
    const MFloat p2_y = radius * sin((i + 1) * step_rad);

    // Segment center
    const MFloat xcenter = 0.5 * (p2_x + p1_x);
    const MFloat ycenter = 0.5 * (p2_y + p1_y);

    // Determine segment length
    const MFloat dx = (p2_x - p1_x);
    const MFloat dy = (p2_y - p1_y);
    const MFloat delta_segment = sqrt(dx * dx + dy * dy);

    // Store segment center coordinates
    MFloat* surfCoordinates = &m_surfaceData.surfaceCoordinates(id);
    surfCoordinates[0] = xcenter;
    surfCoordinates[1] = ycenter;

    // Store segment area (weighted with surface averaging factor)
    m_surfaceData.surfaceArea(i) = weight * delta_segment;

    // Calculate surface normal vector. Simply the gradient and normed
    surfaceNormal[0] = 2 * xcenter / (sqrt(4 * xcenter * xcenter + 4 * ycenter * ycenter));
    surfaceNormal[1] = 2 * ycenter / (sqrt(4 * xcenter * xcenter + 4 * ycenter * ycenter));

    // Store segment normal vector
    MFloat* surfNormal = &m_surfaceData.surfaceNormal(id);
    surfNormal[0] = surfaceNormal[0];
    surfNormal[1] = surfaceNormal[1];
  }
}


/// \brief Generate surface data for an analytical test case.
template <MInt nDim>
void AcaSolver<nDim>::generateSurfaceData() {
  TRACE();

  const MInt noVars = m_surfaceData.noVars();
  TERMM_IF_NOT_COND(noVars == m_noVariables, "number of variables mismatch");

  // TODO labels:ACA Read in analytic or semi-analytic
  const MBool analytic = true;
  /* analytic = Context::getSolverProperty<MBool>("analyticSourceGeneration", m_solverId, AT_,
   * &analytic); */

  // Initialize every variable with zeros
  std::fill_n(&m_surfaceData.variables(0, 0, 0), totalNoSurfaceElements() * noSamples() * noVars, 0.0);

  TERMM_IF_NOT_COND(m_dt > 0.0, "Invalid constant time step size: " + std::to_string(m_dt));
  // Initialize times (constant timestep size)

  for(MInt t = 0; t < noSamples(); t++) {
    m_times[t] = m_dt * t;
  }

  // TODO labels:ACA,toenhance add possibility for superposition of multiple sources with different source parameters
  const MInt noSources = 1;

  for(MInt sourceId = 0; sourceId < noSources; sourceId++) {
    // TODO labels:ACA,toenhance change structure (performance)?
    // Iterate over all surface segments and all observers.
    for(MInt segmentId = 0; segmentId < totalNoSurfaceElements(); segmentId++) {
      const MFloat* coord = &m_surfaceData.surfaceCoordinates(segmentId);

      SourceVars sourceVars;
      if(m_acousticMethod == FWH_METHOD) {
        sourceVars.p = &m_surfaceData.variables(segmentId, FWH::P);
        sourceVars.u = &m_surfaceData.variables(segmentId, FWH::U);
        sourceVars.v = &m_surfaceData.variables(segmentId, FWH::V);
        sourceVars.w = &m_surfaceData.variables(segmentId, FWH::W);
        sourceVars.rho = &m_surfaceData.variables(segmentId, FWH::RHO);
      } else {
        sourceVars.p = &m_surfaceData.variables(segmentId, FWH_APE::PP);
        sourceVars.u = &m_surfaceData.variables(segmentId, FWH_APE::UP);
        sourceVars.v = &m_surfaceData.variables(segmentId, FWH_APE::VP);
        sourceVars.w = &m_surfaceData.variables(segmentId, FWH_APE::WP);
      }

      if(analytic) {
        // Full analytic formulation
        if constexpr(nDim == 2) {
          switch(m_sourceType) {
            case 0: {
              genMonopoleAnalytic2D(coord, m_sourceParameters, sourceVars);
              break;
            }
            case 1: {
              genDipoleAnalytic2D(coord, m_sourceParameters, sourceVars);
              break;
            }
            case 2: {
              genQuadrupoleAnalytic2D(coord, m_sourceParameters, sourceVars);
              break;
            }
            case 3: {
              genVortexConvectionAnalytic2D(coord, m_sourceParameters, sourceVars);
              break;
            }
            default: {
              TERMM(1, "source type not implemented");
              break;
            }
          }
        } else {
          switch(m_sourceType) {
            case 0: {
              genMonopoleAnalytic3D(coord, m_sourceParameters, sourceVars);
              break;
            }
            case 1: {
              genDipoleAnalytic3D(coord, m_sourceParameters, sourceVars);
              break;
            }
            case 2: {
              genQuadrupoleAnalytic3D(coord, m_sourceParameters, sourceVars);
              break;
            }

            default: {
              TERMM(1, "source type not implemented");
              break;
            }
          }
        }
      } else {
        TERMM(1, "Error: semi-analytic source generation not implemented.");
        // TODO labels:ACA Semi-analytic formulation (numeric differentiation of acoustic potential)
        // 1. generate the acoustic potential (without derivatives) for the current
        // segment and a surrounding stencil of points
        // 2. use finite differences to:
        // 2.1 compute the derivatives to get the correct acoustic potential (eg. dipole/quadrupole)
        // 2.2 compute the derivatives to get the variables p', u', ...
        /* if constexpr(nDim == 2) { */
        /*   genSurfaceDataNumeric2D(segmentId, m_sourceParameters); */
        /* } else { */
        /*   genSurfaceDataNumeric3D(segmentId, m_sourceParameters); */
        /* } */
      }
    } // end of segmentId for loop
  }   // end of sourceType loop
}


/// \brief 2D analytical monopole generation
template <MInt nDim>
void AcaSolver<nDim>::genMonopoleAnalytic2D(const MFloat* coord, const SourceParameters& param, SourceVars& vars) {
  MFloat x = coord[0];
  MFloat y = coord[1];
  transformCoordinatesToFlowDirection2D(&x, &y);

  const MFloat omega = param.omega;
  const MFloat c0 = 1.0; // by definition
  const MFloat rho0 = param.rho0;
  const MFloat p0 = rho0 * c0 * c0 / m_gamma;
  const MFloat beta = param.beta;
  const MFloat amplitude = param.amplitude;
  const MFloat u0 = param.Ma * c0;
  const MFloat mach = param.Ma;

  // Precompute time constant values/factors
  const MFloat k = omega / c0;

  const MFloat d = sqrt(x * x + beta * beta * y * y);
  const MFloat xfd = x / d;
  const MFloat yfd = y / d;
  const MFloat beta2 = (beta * beta);
  const MFloat fbeta2 = 1 / beta2;
  const MFloat kfbeta2 = k / (beta * beta);
  const MFloat amplwf4beta = amplitude * omega / (4.0 * beta);
  const MFloat amplkf4beta3 = amplitude * k / (4.0 * beta * beta2);
  const MFloat argHankel = kfbeta2 * d;
  const MFloat J0 = j0(argHankel);
  const MFloat J1 = j1(argHankel);
  const MFloat Y0 = y0(argHankel);
  const MFloat Y1 = y1(argHankel);

  for(MInt sample = 0; sample < noSamples(); sample++) {
    // Note: adapted from acoupot3d/acou_pot.cpp::calc_p/ux/uy/uz_analytic()
    const MFloat time = sample * m_dt;
    const MFloat arg = omega * time + mach * k * x * fbeta2;

    const MFloat cosarg = cos(arg);
    const MFloat sinarg = sin(arg);

    const MFloat dPhidt = -1.0 * amplwf4beta * (cosarg * J0 + sinarg * Y0);
    const MFloat dPhidx =
        -1.0 * mach * amplkf4beta3 * (cosarg * J0 + sinarg * Y0) + amplkf4beta3 * xfd * (sinarg * J1 - cosarg * Y1);
    const MFloat dPhidy = amplkf4beta3 * yfd * beta2 * (-1.0 * cosarg * Y1 + sinarg * J1);
    const MFloat pp = -1.0 * rho0 * (dPhidt + u0 * dPhidx);

    if(param.perturbed) {
      vars.u[sample] = dPhidx;
      vars.v[sample] = dPhidy;
      vars.p[sample] = pp;
    } else {
      vars.u[sample] = dPhidx + u0;
      vars.v[sample] = dPhidy;
      vars.p[sample] = pp + p0;
      vars.rho[sample] = rho0 + pp / (c0 * c0);
    }
    transformBackToOriginalCoordinate2D(&vars.u[sample], &vars.v[sample]);
  }
}


/// \brief 3D analytic monopole generation
template <MInt nDim>
void AcaSolver<nDim>::genMonopoleAnalytic3D(const MFloat* coord, const SourceParameters& param, SourceVars& vars) {
  MFloat x = coord[0];
  MFloat y = coord[1];
  MFloat z = coord[2];
  transformCoordinatesToFlowDirection3D(&x, &y, &z);

  const MFloat omega = param.omega;
  const MFloat c0 = 1.0; // by definition
  const MFloat rho0 = param.rho0;
  const MFloat p0 = rho0 * c0 * c0 / m_gamma;
  const MFloat beta = param.beta;
  const MFloat amplitude = param.amplitude;
  const MFloat u0 = param.Ma * c0;
  const MFloat mach = param.Ma;

  // Precompute time constant values/factors
  // Wavenumber k needs to be dimensionless
  MFloat k = omega / c0;

  const MFloat d = sqrt(x * x + beta * beta * (y * y + z * z));
  const MFloat xfd = x / d;
  const MFloat xfd2 = x / (d * d);
  const MFloat yfd = y / d;
  const MFloat yfd2 = y / (d * d);
  const MFloat zfd = z / d;
  const MFloat zfd2 = z / (d * d);
  const MFloat ampf4pid = amplitude / (4.0 * PI * d);
  const MFloat beta2 = (beta * beta);
  const MFloat fbeta2 = 1.0 / (beta * beta);

  for(MInt sample = 0; sample < noSamples(); sample++) {
    // Note: adapted from acoupot3d/acou_pot.cpp::calc_p/ux/uy/uz_analytic()
    const MFloat time = sample * m_dt;
    const MFloat arg = omega * time - k * (d - mach * x) * fbeta2;

    const MFloat cosarg = cos(arg);
    const MFloat sinarg = sin(arg);

    const MFloat dPhidt = -1.0 * ampf4pid * omega * sinarg;
    const MFloat dPhidx = -1.0 * ampf4pid * xfd2 * cosarg + ampf4pid * k * fbeta2 * sinarg * (xfd - mach);
    const MFloat dPhidy = -1.0 * ampf4pid * beta2 * yfd2 * cosarg + ampf4pid * k * yfd * sinarg;
    const MFloat dPhidz = -1.0 * ampf4pid * beta2 * zfd2 * cosarg + ampf4pid * k * zfd * sinarg;
    const MFloat pp = -1.0 * rho0 * (dPhidt + u0 * dPhidx);

    if(param.perturbed) {
      vars.u[sample] = dPhidx;
      vars.v[sample] = dPhidy;
      vars.w[sample] = dPhidz;
      vars.p[sample] = pp;
    } else {
      vars.u[sample] = dPhidx + u0;
      vars.v[sample] = dPhidy;
      vars.w[sample] = dPhidz;
      vars.p[sample] = pp + p0;
      vars.rho[sample] = rho0 + pp / (c0 * c0);
    }
    transformBackToOriginalCoordinate3D(&vars.u[sample], &vars.v[sample], &vars.w[sample]);
  }
}


/// \brief 2D analytical dipole generation
template <MInt nDim>
void AcaSolver<nDim>::genDipoleAnalytic2D(const MFloat* coord, const SourceParameters& param, SourceVars& vars) {
  MFloat x = coord[0];
  MFloat y = coord[1];
  transformCoordinatesToFlowDirection2D(&x, &y);

  const MFloat omega = param.omega;
  const MFloat c0 = 1.0; // by definition
  const MFloat rho0 = param.rho0;
  const MFloat p0 = rho0 * c0 * c0 / m_gamma;
  const MFloat beta = param.beta;
  const MFloat G = param.amplitude;
  const MFloat u0 = param.Ma * c0;
  const MFloat mach = param.Ma;

  // Precompute time constant values/factors
  const MFloat beta2 = (beta * beta);
  const MFloat fbeta2 = 1 / beta2;

  MFloat k = omega / c0;

  const MFloat d = sqrt(x * x + beta2 * y * y);
  const MFloat dargdx = mach * k * fbeta2;
  const MFloat argHankel = k * d * fbeta2;
  const MFloat dHdx = k * x * fbeta2 / d;
  const MFloat dHdxdx = k * y * y / (d * d * d);
  const MFloat dHdxdy = -1.0 * k * x * y / (d * d * d);
  const MFloat dHdy = k * y / d;
  const MFloat C1 = -1.0 * mach * G * k / (4.0 * beta2 * beta);
  const MFloat C2 = G / (4.0 * beta);

  const MFloat J0 = j0(argHankel);
  const MFloat J1 = j1(argHankel);
  const MFloat dJ1 = J0 - (1 / argHankel) * J1;

  const MFloat Y0 = y0(argHankel);
  const MFloat Y1 = y1(argHankel);
  const MFloat dY1 = Y0 - (1 / argHankel) * Y1;

  for(MInt sample = 0; sample < noSamples(); sample++) {
    const MFloat time = sample * m_dt;
    const MFloat arg = omega * time + mach * k * x * fbeta2;

    const MFloat cosarg = cos(arg);
    const MFloat sinarg = sin(arg);

    const MFloat dPhidt =
        C1 * omega * (-1.0 * sinarg * J0 + cosarg * Y0) + C2 * dHdx * omega * (cosarg * J1 + sinarg * Y1);
    const MFloat dPhidx = -1.0 * C1 * (dargdx * (sinarg * J0 - cosarg * Y0) + dHdx * (cosarg * J1 + sinarg * Y1))
                          + C2
                                * (dHdxdx * (sinarg * J1 - cosarg * Y1) + dHdx * dargdx * (cosarg * J1 + sinarg * Y1)
                                   + dHdx * dHdx * (sinarg * dJ1 - cosarg * dY1));

    const MFloat dPhidy = -1.0 * C1 * dHdy * (cosarg * J1 + sinarg * Y1)
                          + C2 * (dHdxdy * (sinarg * J1 - cosarg * Y1) + dHdx * dHdy * (sinarg * dJ1 - cosarg * dY1));
    const MFloat pp = -1.0 * rho0 * (dPhidt + u0 * dPhidx);

    if(param.perturbed) {
      vars.u[sample] = dPhidx;
      vars.v[sample] = dPhidy;
      vars.p[sample] = pp;
    } else {
      vars.u[sample] = dPhidx + u0;
      vars.v[sample] = dPhidy;
      vars.p[sample] = pp + p0;
      vars.rho[sample] = rho0 + pp / (c0 * c0);
    }
    transformBackToOriginalCoordinate2D(&vars.u[sample], &vars.v[sample]);
  }
}


/// \brief 3D dipole analytic
template <MInt nDim>
void AcaSolver<nDim>::genDipoleAnalytic3D(const MFloat* coord, const SourceParameters& param, SourceVars& vars) {
  MFloat x = coord[0];
  MFloat y = coord[1];
  MFloat z = coord[2];
  transformCoordinatesToFlowDirection3D(&x, &y, &z);

  const MFloat omega = param.omega;
  const MFloat c0 = 1.0; // by definition
  const MFloat rho0 = param.rho0;
  const MFloat p0 = rho0 * c0 * c0 / m_gamma;
  const MFloat beta = param.beta;
  const MFloat amplitude = param.amplitude;
  const MFloat u0 = param.Ma * c0;
  const MFloat mach = param.Ma;

  MFloat k = omega / c0;

  const MFloat beta2 = beta * beta;
  const MFloat fbeta2 = 1 / beta2;
  const MFloat d = sqrt(x * x + beta * beta * (y * y + z * z));
  const MFloat d2 = d * d;
  const MFloat fd2 = 1 / d2;
  const MFloat fd3 = 1 / (d2 * d);
  const MFloat C1 = amplitude / (4 * PI);
  const MFloat C2 = k * fbeta2;

  const MFloat P1 = fd3 - 3 * x * x * fd2 * fd3;
  const MFloat P2 = 1 * fd2 - 2 * x * x * fd2 * fd2;
  const MFloat P3 = -1.0 * x * fd3;

  const MFloat dargdx = -1.0 * C2 * (x / d - mach);
  const MFloat dargdy = -1.0 * k * y / d;
  const MFloat dargdz = -1.0 * k * z / d;

  for(MInt sample = 0; sample < noSamples(); sample++) {
    // Note: adapted from acoupot3d/acou_pot.cpp::calc_p/ux/uy/uz_analytic()
    const MFloat time = sample * m_dt;
    const MFloat arg = omega * time - k * (d - mach * x) * fbeta2;

    const MFloat cosarg = cos(arg);
    const MFloat sinarg = sin(arg);

    const MFloat dPhidt =
        C1 * omega * fd3 * x * sinarg + C1 * omega * fd2 * x * C2 * cosarg - C1 * omega * C2 * (mach / d) * cosarg;
    const MFloat dPhidx = -1.0 * C1 * P1 * cosarg + C1 * x * fd3 * dargdx * sinarg + C1 * P2 * C2 * sinarg
                          + C1 * x * fd2 * C2 * dargdx * cosarg - C1 * C2 * P3 * mach * sinarg
                          - C1 * C2 * (mach / d) * dargdx * cosarg;
    const MFloat dPhidy = C1 * x * 3 * y * beta2 * fd3 * fd2 * cosarg + C1 * x * fd3 * dargdy * sinarg
                          - C1 * x * C2 * 2 * y * beta2 * fd2 * fd2 * sinarg + C1 * x * C2 * fd2 * dargdy * cosarg
                          + C1 * C2 * y * beta2 * mach * fd3 * sinarg - C1 * C2 * (mach / d) * dargdy * cosarg;
    const MFloat dPhidz = C1 * x * 3 * z * beta2 * fd3 * fd2 * cosarg + C1 * x * fd3 * dargdz * sinarg
                          - C1 * x * C2 * 2 * z * beta2 * fd2 * fd2 * sinarg + C1 * x * C2 * fd2 * dargdz * cosarg
                          + C1 * C2 * z * beta2 * mach * fd3 * sinarg - C1 * C2 * (mach / d) * dargdz * cosarg;
    const MFloat pp = -1.0 * rho0 * (dPhidt + u0 * dPhidx);

    if(param.perturbed) {
      vars.u[sample] = dPhidx;
      vars.v[sample] = dPhidy;
      vars.w[sample] = dPhidz;
      vars.p[sample] = pp;
    } else {
      vars.u[sample] = dPhidx + u0;
      vars.v[sample] = dPhidy;
      vars.w[sample] = dPhidz;
      vars.p[sample] = pp + p0;
      vars.rho[sample] = rho0 + pp / (c0 * c0);
    }
    transformBackToOriginalCoordinate3D(&vars.u[sample], &vars.v[sample], &vars.w[sample]);
  }
}


/// \brief 2D quadrupole analytic
template <MInt nDim>
void AcaSolver<nDim>::genQuadrupoleAnalytic2D(const MFloat* coord, const SourceParameters& param, SourceVars& vars) {
  MFloat x = coord[0];
  MFloat y = coord[1];
  transformCoordinatesToFlowDirection2D(&x, &y);

  const MFloat omega = param.omega;
  const MFloat c0 = 1.0; // by definition
  const MFloat rho0 = param.rho0;
  const MFloat p0 = rho0 * c0 * c0 / m_gamma;
  const MFloat beta = param.beta;
  const MFloat amplitude = param.amplitude;
  const MFloat u0 = param.Ma * c0;
  const MFloat mach = param.Ma;

  MFloat k = omega / c0;

  const MFloat beta2 = beta * beta;
  const MFloat fbeta2 = 1 / beta2;
  const MFloat d = sqrt(x * x + beta2 * y * y);
  const MFloat dargdx = mach * k * fbeta2;
  const MFloat argHankel = k * d * fbeta2;
  const MFloat K1 = -1.0 * mach * amplitude * k / (4.0 * beta2 * beta);
  const MFloat K2 = amplitude / (4.0 * beta);

  const MFloat dHdx = k * x * fbeta2 / d;
  const MFloat dHdy = k * y / d;

  const MFloat dHdxdx = k * y * y / (d * d * d);
  const MFloat dHdxdy = -1.0 * k * x * y / (d * d * d);
  const MFloat dHdydx = dHdxdy;
  const MFloat dHdydy = k * x * x / (d * d * d);

  const MFloat dHdxdydx = k * y * (2 * x * x - beta2 * y * y) / (d * d * d * d * d);
  const MFloat dHdxdydy = k * x * (2 * beta2 * y * y - x * x) / (d * d * d * d * d);

  const MFloat dHdxfH = x / (d * d);
  const MFloat dHdxfHdx = -1.0 * (x * x - beta2 * y * y) / (d * d * d * d);
  const MFloat dHdxfHdy = -1.0 * 2 * beta2 * x * y / (d * d * d * d);


  const MFloat J0 = j0(argHankel);
  const MFloat J1 = j1(argHankel);
  const MFloat dJ1 = J0 - (1 / argHankel) * J1;

  const MFloat Y0 = y0(argHankel);
  const MFloat Y1 = y1(argHankel);
  const MFloat dY1 = Y0 - (1 / argHankel) * Y1;

  for(MInt sample = 0; sample < noSamples(); sample++) {
    const MFloat time = sample * m_dt;
    const MFloat arg = omega * time + mach * k * x * fbeta2;

    const MFloat cosarg = cos(arg);
    const MFloat sinarg = sin(arg);

    const MFloat dPhidt = K1 * dHdy * sinarg * omega * J1 - K1 * dHdy * cosarg * omega * Y1
                          + K2 * dHdxdy * cosarg * omega * J1 + K2 * dHdxdy * sinarg * omega * Y1
                          + K2 * dHdx * dHdy * cosarg * omega * J0 - K2 * dHdxfH * dHdy * cosarg * omega * J1
                          + K2 * dHdx * dHdy * sinarg * omega * Y0 - K2 * dHdxfH * dHdy * sinarg * omega * Y1;

    const MFloat dPhidx =
        -1.0 * K1 * dHdydx * cosarg * J1 + K1 * dHdy * sinarg * dargdx * J1 - K1 * dHdy * cosarg * dJ1 * dHdx
        - K1 * dHdydx * sinarg * Y1 - K1 * dHdy * cosarg * dargdx * Y1 - K1 * dHdy * sinarg * dY1 * dHdx
        + K2 * dHdxdydx * sinarg * J1 + K2 * dHdxdy * cosarg * dargdx * J1 + K2 * dHdxdy * sinarg * dJ1 * dHdx
        - K2 * dHdxdydx * cosarg * Y1 + K2 * dHdxdy * sinarg * dargdx * Y1 - K2 * dHdxdy * cosarg * dY1 * dHdx
        + K2 * dHdxdx * dHdy * sinarg * J0 + K2 * dHdx * dHdydx * sinarg * J0 + K2 * dHdx * dHdy * cosarg * dargdx * J0
        - K2 * dHdx * dHdy * sinarg * J1 * dHdx - K2 * dHdxfHdx * dHdy * sinarg * J1
        - K2 * dHdxfH * dHdydx * sinarg * J1 - K2 * dHdxfH * dHdy * cosarg * dargdx * J1
        - K2 * dHdxfH * dHdy * sinarg * dJ1 * dHdx - K2 * dHdxdx * dHdy * cosarg * Y0 - K2 * dHdx * dHdydx * cosarg * Y0
        + K2 * dHdx * dHdy * sinarg * dargdx * Y0 + K2 * dHdx * dHdy * cosarg * Y1 * dHdx
        + K2 * dHdxfHdx * dHdy * cosarg * Y1 + K2 * dHdxfH * dHdydx * cosarg * Y1
        - K2 * dHdxfH * dHdy * sinarg * dargdx * Y1 + K2 * dHdxfH * dHdy * cosarg * dY1 * dHdx;

    const MFloat dPhidy =
        -1.0 * K1 * dHdydy * cosarg * J1 - K1 * dHdy * cosarg * dJ1 * dHdy - K1 * dHdydy * sinarg * Y1
        - K1 * dHdy * sinarg * dY1 * dHdy + K2 * dHdxdydy * sinarg * J1 + K2 * dHdxdy * sinarg * dJ1 * dHdy
        - K2 * dHdxdydy * cosarg * Y1 - K2 * dHdxdy * cosarg * dY1 * dHdy + K2 * dHdxdy * dHdy * sinarg * J0
        + K2 * dHdx * dHdydy * sinarg * J0 - K2 * dHdx * dHdy * sinarg * J1 * dHdy - K2 * dHdxfHdy * dHdy * sinarg * J1
        - K2 * dHdxfH * dHdydy * sinarg * J1 - K2 * dHdxfH * dHdy * sinarg * dJ1 * dHdy
        - K2 * dHdxdy * dHdy * cosarg * Y0 - K2 * dHdx * dHdydy * cosarg * Y0 + K2 * dHdx * dHdy * cosarg * Y1 * dHdy
        + K2 * dHdxfHdy * dHdy * cosarg * Y1 + K2 * dHdxfH * dHdydy * cosarg * Y1
        + K2 * dHdxfH * dHdy * cosarg * dY1 * dHdy;

    const MFloat pp = -1.0 * rho0 * (dPhidt + u0 * dPhidx);

    if(param.perturbed) {
      vars.u[sample] = dPhidx;
      vars.v[sample] = dPhidy;
      vars.p[sample] = pp;
    } else {
      vars.u[sample] = dPhidx + u0;
      vars.v[sample] = dPhidy;
      vars.p[sample] = pp + p0;
      vars.rho[sample] = rho0 + pp / (c0 * c0);
    }
    transformBackToOriginalCoordinate2D(&vars.u[sample], &vars.v[sample]);
  }
}


/// \brief 3D quadrupole analytic
template <MInt nDim>
void AcaSolver<nDim>::genQuadrupoleAnalytic3D(const MFloat* coord, const SourceParameters& param, SourceVars& vars) {
  MFloat x = coord[0];
  MFloat y = coord[1];
  MFloat z = coord[2];
  transformCoordinatesToFlowDirection3D(&x, &y, &z);

  const MFloat omega = param.omega;
  const MFloat c0 = 1.0; // by definition
  const MFloat rho0 = param.rho0;
  const MFloat p0 = rho0 * c0 * c0 / m_gamma;
  const MFloat beta = param.beta;
  const MFloat amplitude = param.amplitude;
  const MFloat u0 = param.Ma * c0;
  const MFloat mach = param.Ma;

  MFloat k = omega / c0;

  const MFloat beta2 = beta * beta;
  const MFloat fbeta2 = 1 / beta2;
  const MFloat d = sqrt(x * x + beta * beta * (y * y + z * z));
  const MFloat d2 = d * d;
  const MFloat d3 = d2 * d;
  const MFloat fd2 = 1 / d2;
  const MFloat fd3 = 1 / d3;
  const MFloat C1 = amplitude / (4 * PI);
  const MFloat C2 = k * fbeta2;

  const MFloat dargdx = -1.0 * C2 * (x / d - mach);
  const MFloat dargdy = -1.0 * k * y / d;
  const MFloat dargdz = -1.0 * k * z / d;

  for(MInt sample = 0; sample < noSamples(); sample++) {
    // Note: adapted from acoupot3d/acou_pot.cpp::calc_p/ux/uy/uz_analytic()
    const MFloat time = sample * m_dt;
    const MFloat arg = omega * time - k * (d - mach * x) * fbeta2;

    const MFloat cosarg = cos(arg);
    const MFloat sinarg = sin(arg);

    const MFloat dPhidt = -1.0 * C1 * 3 * x * y * beta2 * fd2 * fd3 * sinarg * omega
                          - C1 * 3 * k * x * y * fd2 * fd2 * cosarg * omega + C1 * k * C2 * x * y * fd3 * sinarg * omega
                          + C1 * mach * k * y * fd3 * cosarg * omega - C1 * mach * k * C2 * y * fd2 * sinarg * omega;

    const MFloat dPhidx =
        C1 * 3 * y * beta2 * (fd2 * fd3 - 5 * x * x * fd2 * fd2 * fd3) * cosarg
        - C1 * 3 * x * y * beta2 * fd2 * fd3 * sinarg * dargdx
        - C1 * 3 * k * y * (fd2 * fd2 - 4 * x * x * fd3 * fd3) * sinarg
        - C1 * 3 * k * x * y * fd2 * fd2 * cosarg * dargdx - C1 * k * C2 * y * (fd3 - 3 * x * x * fd3 * fd2) * cosarg
        + C1 * k * C2 * x * y * fd3 * sinarg * dargdx - C1 * mach * k * y * 3 * x * fd3 * fd2 * sinarg
        + C1 * mach * k * y * fd3 * cosarg * dargdx - C1 * mach * k * C2 * y * 2 * x * fd2 * fd2 * cosarg
        - C1 * mach * k * C2 * y * fd2 * sinarg * dargdx;

    const MFloat dPhidy =
        C1 * 3 * x * beta2 * (fd2 * fd3 - 5 * beta2 * y * y * fd2 * fd2 * fd3) * cosarg
        - C1 * 3 * x * y * beta * beta * fd3 * fd2 * sinarg * dargdy
        - C1 * 3 * k * x * (fd2 * fd2 - 4 * beta2 * y * y * fd3 * fd3) * sinarg
        - C1 * 3 * k * x * y * fd2 * fd2 * cosarg * dargdy
        - C1 * k * C2 * x * (fd3 - 3 * beta2 * y * y * fd2 * fd3) * cosarg + C1 * k * C2 * x * y * fd3 * sinarg * dargdy
        + C1 * mach * k * (fd3 - 3 * beta2 * y * y * fd2 * fd3) * sinarg + C1 * mach * k * y * fd3 * cosarg * dargdy
        + C1 * mach * k * C2 * (fd2 - 2 * beta2 * y * y * fd2 * fd2) * cosarg
        - C1 * mach * k * C2 * y * fd2 * sinarg * dargdy;

    const MFloat dPhidz =
        -1.0 * C1 * 15 * x * y * z * beta2 * beta2 * fd2 * fd2 * fd3 * cosarg
        - C1 * 3 * x * y * beta2 * fd2 * fd3 * sinarg * dargdz + C1 * 12 * k * x * y * z * beta2 * fd3 * fd3 * sinarg
        - C1 * 3 * k * x * y * fd2 * fd2 * cosarg * dargdz + C1 * 3 * k * k * x * y * z * fd2 * fd3 * cosarg
        + C1 * k * C2 * x * y * fd3 * sinarg * dargdz - C1 * 3 * mach * k * y * z * beta2 * fd2 * fd3 * sinarg
        + C1 * mach * k * y * fd3 * cosarg * dargdz - C1 * 2 * mach * k * k * y * z * fd2 * fd2 * cosarg
        - C1 * mach * k * C2 * y * fd2 * sinarg * dargdz;

    const MFloat pp = -1.0 * rho0 * (dPhidt + u0 * dPhidx);

    if(param.perturbed) {
      vars.u[sample] = dPhidx;
      vars.v[sample] = dPhidy;
      vars.w[sample] = dPhidz;
      vars.p[sample] = pp;
    } else {
      vars.u[sample] = dPhidx + u0;
      vars.v[sample] = dPhidy;
      vars.w[sample] = dPhidz;
      vars.p[sample] = pp + p0;
      vars.rho[sample] = rho0 + pp / (c0 * c0);
    }
    transformBackToOriginalCoordinate3D(&vars.u[sample], &vars.v[sample], &vars.w[sample]);
  }
}

/** \brief  Calculate solution vars for analytical inviscid convecting vortex
 *  \author Miro Gondrum
 *  \date   11.05.2023
 *
 *  A vortex is convected in an inviscid flow with M_infty in x-direction. The
 *  analytical solution is known to not generate noise. Hence, all noise
 *  generate is of erroneous nature.
 */
template <MInt nDim>
void AcaSolver<nDim>::genVortexConvectionAnalytic2D(const MFloat* coord, const SourceParameters& param,
                                                    SourceVars& vars) {
  const MFloat x = coord[0];
  const MFloat y = coord[1];
  constexpr MFloat cs = 1.0; // by definition
  const MFloat rho0 = param.rho0;

  const MFloat fac = cs / (2.0 * PI);
  for(MInt sample = 0; sample < noSamples(); sample++) {
    const MFloat t = sample * m_dt;
    const MFloat xRel = x - param.Ma * cs * t;
    const MFloat yRel = y;
    const MFloat r = sqrt(POW2(xRel) + POW2(yRel));
    const MFloat theta = atan2(yRel, xRel);
    const MFloat utheta = fac * r * exp(0.5 * (1 - r * r));
    const MFloat up = utheta * sin(theta);
    const MFloat vp = -utheta * cos(theta);
    const MFloat pp = pow(1 + (m_gamma - 1) * exp(1 - r * r) / (8.0 * PI * PI), m_gamma / (m_gamma - 1.0)) - 1.0;
    if(param.perturbed) {
      vars.u[sample] = up;
      vars.v[sample] = vp;
      vars.p[sample] = pp;
    } else {
      vars.u[sample] = param.Ma * cs + up;
      vars.v[sample] = vp;
      vars.p[sample] = pp;
      if(vars.rho != nullptr) vars.rho[sample] = rho0 * pow((1.0 + pp), 1.0 / m_gamma);
    }
  }
}

/** \brief  Apply symmetry boundary condition
 *  \author Miro Gondrum
 *  \date   15.06.2023
 *  \param[in]  coord               Location of the observer
 *  \param[out] p_complexVariables  Output of complex p', length of 2*noSamples()
 *
 *  Calculates surface integral for mirrored observer points. The resulting
 *  complex pressure is then added to the image observer points.
 */
template <MInt nDim>
void AcaSolver<nDim>::applySymmetryBc(const std::array<MFloat, nDim> coord, MFloat* const p_complexVariables) {
  TRACE();
  const MInt noSym = m_symBc.size();
  if(noSym == 0) return;
  TERMM_IF_COND(noSym > 1, "WARNING: only implemented for one symmetry BC so far");
  auto transform = [](std::array<MFloat, nDim>& coord_, const AcaSymBc& bc) {
    const MFloat coordTimesNormal = std::inner_product(&coord_[0], &coord_[nDim], &bc.normal[0], .0);
    for(MInt d = 0; d < nDim; d++) {
      coord_[d] += -2 * coordTimesNormal * bc.normal[d] + 2 * bc.origin[d];
    }
  };
  MFloatScratchSpace complexVariables_(2 * noSamples(), AT_, "complexVariables_");
  // transform observer position
  std::array<MFloat, nDim> symCoord = coord;
  transform(symCoord, m_symBc[0]);
  // Calculate pressure for transformed point before adding to its 'image' observer
  calcSurfaceIntegralsForObserver(symCoord, complexVariables_.data());
  for(MInt i = 0; i < noSamples(); i++) {
    p_complexVariables[2 * i + 0] += complexVariables_[2 * i + 0];
    p_complexVariables[2 * i + 1] += complexVariables_[2 * i + 1];
  }
}

/// \brief Initialize the observer points (load/generate coordinates and init collector)
template <MInt nDim>
void AcaSolver<nDim>::initObserverPoints() {
  TRACE();
  RECORD_TIMER_START(m_timers[Timers::InitObservers]);

  MInt noGlobalObservers = 0;
  printMessage("- Init observer points");
  std::stringstream ss;

  if(m_generateObservers) { // TODO labels:ACA,toenhance make more general/flexible
    MString observerGenerationType = "CIRCLE";
    observerGenerationType =
        Context::getSolverProperty<MString>("observerGenerationType", m_solverId, AT_, &observerGenerationType);
    /*! \property
      \page propertyPageACA ACA
      \section noObservers
      <code>MInt noGlobalObservers</code>\n
      default = <code>100</code>\n\n
      Number of observer points that are created when `generateObservers` is
      `True`.\n
      Keywords: <i>ACA, observer</i>
    */
    noGlobalObservers = 100;
    noGlobalObservers = Context::getSolverProperty<MInt>("noObservers", m_solverId, AT_, &noGlobalObservers);
    if(noGlobalObservers < 1) {
      TERMM(1, "The number of Observers cannot be lower than 1!");
    }
    m_globalObserverCoordinates.resize(noGlobalObservers * nDim);

    ss << "  "
       << "Generate " << noGlobalObservers << " observer with mode " << observerGenerationType;

    if(observerGenerationType == MString("CIRCLE")) {
      // Default: 100 points on a circle with R=50

      const MFloat step_rad = 2 * PI / static_cast<MFloat>(noGlobalObservers);
      /*! \property
        \page propertyPageACA ACA
        \section observerRadius
        <code>MFloat radius</code>\n
        default = <code>50.0</code>\n\n
        Radius of t he circle on which observer points are distributed.\n
        Keywords: <i>ACA, observer</i>
      */
      MFloat radius = 50.0;
      radius = Context::getSolverProperty<MFloat>("observerRadius", m_solverId, AT_, &radius);
      if(radius <= 0.0) {
        TERMM(1, "The radius cannot be equal or lower than 0!");
      }

      for(MInt i = 0; i < noGlobalObservers; i++) {
        // TODO labels:ACA properties for center
        const MFloat x = radius * cos(i * step_rad);
        const MFloat y = radius * sin(i * step_rad);

        m_globalObserverCoordinates[i * nDim] = x;
        m_globalObserverCoordinates[i * nDim + 1] = y;
        if constexpr(nDim == 3) {
          m_globalObserverCoordinates[i * nDim + 2] = 0.0; // TODO labels:ACA make more flexible for 3D
        }
      }

    } else if(observerGenerationType == MString("PLANE")) {
      // input
      MInt normalDir = 2;
      if constexpr(nDim > 2) {
        normalDir = Context::getSolverProperty<MInt>("observerGenerationPlaneNormal", m_solverId, AT_, &normalDir);
      }
      MFloat dim[4];
      for(MInt i = 0; i < 4; i++) {
        dim[i] = Context::getSolverProperty<MFloat>("observerGenerationPlaneDimension", m_solverId, AT_, i);
      }
      // Specify spacing dx[] in the plane after determining number of points
      // per direction nx[]
      const MFloat l[2] = {dim[2] - dim[0], dim[3] - dim[1]};
      MInt nx[2];
      nx[0] = sqrt(noGlobalObservers * l[0] / l[1]);
      nx[1] = noGlobalObservers / nx[0];
      const MFloat dx[2] = {l[0] / nx[0], l[1] / nx[1]};
      // generate points
      for(MInt i = 0; i < noGlobalObservers; i++) {
        const MInt j = i % nx[0];
        const MInt k = i / nx[0];
        m_globalObserverCoordinates[i * nDim + 0] = dim[0] + dx[0] * j;
        m_globalObserverCoordinates[i * nDim + 1] = dim[1] + dx[1] * k;
      }
    } else {
      ss << " is not possible.";
    }
  } else if(m_observerFileName != "") {
    // Read number of observers and m_globalObserverCoordinates from given file on rank 0
    if(domainId() == 0) {
      noGlobalObservers = loadPointCoordinatesFromFile(m_observerFileName, nDim, m_globalObserverCoordinates);

      if(noGlobalObservers < 1) {
        TERMM(1, "Error: no observer points found in file " + m_observerFileName);
      }
      ss << "  "
         << "Read " << noGlobalObservers << " observer points from file " << m_observerFileName;
    }

    // Broadcast number of observers to all ranks
    MPI_Bcast(&noGlobalObservers, 1, type_traits<MInt>::mpiType(), 0, mpiComm(), AT_, "noObservers");

    // Resize storage for coordinates on other ranks
    if(domainId() != 0) {
      m_globalObserverCoordinates.resize(noGlobalObservers * nDim);
    }

    // Broadcast observer coordinates from rank 0
    MPI_Bcast(m_globalObserverCoordinates.data(), noGlobalObservers * nDim, type_traits<MFloat>::mpiType(), 0,
              mpiComm(), AT_, "coordinates");
  } else {
    TERMM(1, "Error: no observer file name specified and generation of observers not enabled.");
  }

  m_noGlobalObservers = noGlobalObservers;

  // number of observer (complex) variables, for now just the acoustic pressure
  const MInt noObsVars = 1;
  const MInt noObsComplexVars = 1;

  // Set data sizes in observer collector
  m_observerData.setNoVars(noObsVars);
  m_observerData.setNoComplexVars(noObsComplexVars);
  m_observerData.setNoSamples(m_noSamples);

  // Calculate offsets and number of local observers
  distributeObservers();

  // Resize observer collector only for local observer data !
  m_observerData.reset(noObservers());

  // Add (local) observers to collector and set data
  for(MInt localObserverId = 0; localObserverId < noObservers(); localObserverId++) {
    m_observerData.append();
    const MInt globalObserverId = a_globalObserverId(localObserverId);
    MFloat* obsCoord = &m_observerData.observerCoordinates(localObserverId);
    for(MInt d = 0; d < nDim; d++) {
      obsCoord[d] = m_globalObserverCoordinates[globalObserverId * nDim + d];
    }
  }

  printMessage(ss.str());
  RECORD_TIMER_STOP(m_timers[Timers::InitObservers]);
}

/// \brief Initialize the conversion factors required to convert between
//         in-/output dimension and the dimensions used in ACA
template <MInt nDim>
void AcaSolver<nDim>::initConversionFactors() {
  // TODO labels:ACA IMHO it does not make sense to have a different dimensions
  // or even quantities (wave number or frequency) in the output. The output
  // should have the same dimensions for a certain solver to avoid confusion
  // when switching between different input solver.
  std::stringstream ss;
  printMessage("- Init conversion factors");
  switch(m_acaResultsNondimMode) {
    case NONDIMACA:
    case NONDIMSTAGNATION: {
      const MFloat constant = 1.0 + (m_gamma - 1.0) * 0.5 * m_MaDim * m_MaDim;
      m_input2Aca.length = 1.0;
      m_input2Aca.density = pow(constant, (1.0 / (m_gamma - 1.0)));
      m_input2Aca.velocity = sqrt(constant);
      m_input2Aca.time = m_input2Aca.length / m_input2Aca.velocity;
      m_input2Aca.pressure = pow(constant, (m_gamma / (m_gamma - 1.0)));
      constexpr MFloat HzToRad = 2.0 * PI;
      m_input2Aca.frequency = HzToRad / m_input2Aca.time;
      // .. and corresponding inverse factors
      if(m_acaResultsNondimMode == NONDIMACA) {
        // Default definition is correct -> dump out in ACA dimension
      } else if(m_acaResultsNondimMode == NONDIMSTAGNATION) {
        m_aca2Output.length = 1.0 / m_input2Aca.length;
        m_aca2Output.density = 1.0 / m_input2Aca.density;
        m_aca2Output.velocity = 1.0 / m_input2Aca.velocity;
        m_aca2Output.time = 1.0 / m_input2Aca.time;
        m_aca2Output.pressure = 1.0 / m_input2Aca.pressure;
        m_aca2Output.frequency = 1.0 / m_input2Aca.frequency;
      }
      break;
    }
    case NONDIMLB: {
      // To reduce some confusion of the LB user. In LB all quantities are made
      // non-dimensional with dxLb (cell length at max level), density at
      // infinity, and xi_0=dxLb/dtLb. But, the time in the surface sampling
      // file is given by t_infty = t * u_infty / L_FV , L_FV=1.
      m_input2Aca.density = 1.0;
      m_input2Aca.time = 1.0 / m_MaDim; // not LB to ACA, as in sampling file t_infty is dumped
      // TODO labels:ACA How to choose dxLb (question of input only)
      // const MFloat dxLb = Context::getSolverProperty<MFloat>("dxLb", m_solverId, AT_, 0);
      // m_input2Aca.length = dxLb;
      m_input2Aca.length = std::numeric_limits<double>::quiet_NaN();
      m_input2Aca.velocity = F1BCS;
      m_input2Aca.pressure = m_input2Aca.density * POW2(m_input2Aca.velocity);
      break;
    }
    default: {
      TERMM(1, "Invalid acaResultsNondimMode: '" + std::to_string(m_acaResultsNondimMode) + "'");
    }
  }
  ss << "  MaDim = " << m_MaDim << std::endl;
  ss << "  Conversion factors input2Aca | aca2Output:" << std::endl;
  ss << "    velocity = " << m_input2Aca.velocity << " | " << m_aca2Output.velocity << std::endl;
  ss << "    pressure = " << m_input2Aca.pressure << " | " << m_aca2Output.pressure << std::endl;
  ss << "    density  = " << m_input2Aca.density << " | " << m_aca2Output.density << std::endl;
  ss << "    time     = " << m_input2Aca.time << " | " << m_aca2Output.time << std::endl;
  ss << "    frequency= " << m_input2Aca.frequency << " | " << m_aca2Output.frequency;
  printMessage(ss.str());
}

/// \brief Change of the non-dimensionalization of input data from the stagnation state (0) non-dim.
///        to the undisturbed freestream state (inf) non-dim.
template <MInt nDim>
void AcaSolver<nDim>::changeDimensionsSurfaceData() {
  // All variables which are read in need to be non-dimensionalized with the freestream c and rho.
  // Index 0 in this (and only in this) function descripes the stagnation state.
  m_log << "Change dimensions from stagnation to freestream state" << std::endl;

  // APE or not APE? That is the question
  const MInt dummy_u = (m_acousticMethod == FWH_METHOD) ? FWH::U : FWH_APE::UP;
  const MInt dummy_v = (m_acousticMethod == FWH_METHOD) ? FWH::V : FWH_APE::VP;
  const MInt dummy_w = (m_acousticMethod == FWH_METHOD) ? FWH::W : FWH_APE::WP;
  const MInt dummy_p = (m_acousticMethod == FWH_METHOD) ? FWH::P : FWH_APE::PP;

  // Apply norm factor to the variables
  for(MInt seg = 0; seg < totalNoSurfaceElements(); seg++) {
    for(MInt t = 0; t < noSamples(); t++) {
      m_surfaceData.variables(seg, dummy_u, t) *= m_input2Aca.velocity;
      m_surfaceData.variables(seg, dummy_v, t) *= m_input2Aca.velocity;
      if constexpr(nDim == 3) m_surfaceData.variables(seg, dummy_w, t) *= m_input2Aca.velocity;
      m_surfaceData.variables(seg, dummy_p, t) *= m_input2Aca.pressure;
      if(m_acousticMethod == FWH_METHOD) {
        m_surfaceData.variables(seg, FWH::RHO, t) *= m_input2Aca.density;
      }
    } // end of sample
  }   // end of seg

  // Norm factor for time:
  for(MInt i = 0; i < noSamples(); i++) {
    m_times[i] *= m_input2Aca.time;
  }
}


template <MInt nDim>
void AcaSolver<nDim>::changeDimensionsObserverData() {
  printMessage("- Change dimension of observer data");
  if(approx(m_aca2Output.frequency, 1.0, MFloatEps)) {
    for(MInt i = 0; i < noSamples(); i++) {
      m_frequencies[i] *= m_aca2Output.frequency;
    }
  }
  if(approx(m_aca2Output.pressure, 1.0, MFloatEps)) {
    for(MInt obs = 0; obs < noObservers(); obs++) {
      for(MInt t = 0; t < noSamples(); t++) {
        m_observerData.variables(obs, 0, t) *= m_aca2Output.pressure;
      }
    }
    for(MInt obs = 0; obs < noObservers(); obs++) {
      for(MInt f = 0; f < noSamples(); f++) {
        m_observerData.complexVariables(obs, 0, f, 0) *= m_aca2Output.pressure;
        m_observerData.complexVariables(obs, 0, f, 1) *= m_aca2Output.pressure;
      }
    }
  }
  if(approx(m_aca2Output.time, 1.0, MFloatEps)) {
    m_dt *= m_aca2Output.time;
    for(MInt i = 0; i < noSamples(); i++) {
      m_times[i] *= m_aca2Output.time;
    }
  }
}

template <MInt nDim>
void AcaSolver<nDim>::computeDt() {
  // Compute time step size
  m_dt = (m_times[noSamples() - 1] - m_times[0]) / (static_cast<MFloat>(noSamples() - 1));

  // Check for equidistant timesteps - necessary for the computation!
  for(MInt i = 1; i < noSamples(); i++) {
    const MFloat delta = m_times[i] - m_times[i - 1];
    const MFloat relErr = std::fabs(m_dt - delta) / m_dt;
    if(relErr > 1e-10) {
      std::ostringstream err;
      err << std::setprecision(12) << std::scientific << "Error: time step size mismatch; m_dt = " << m_dt
          << "; delta = " << delta << "; relErr = " << relErr << "; i = " << i << "; times[i] = " << m_times[i]
          << "; times[i-1] = " << m_times[i - 1];
      TERMM(1, err.str());
    }
  }
}

template <MInt nDim>
void AcaSolver<nDim>::checkNoSamplesPotencyOfTwo() {
  // Check if noSamples is a potency of 2 to see if FFT can be used
  MInt numSampl = noSamples();
  MInt prove = 1;
  m_FastFourier = false;

  // if potency of 2: numSampl == 0; if not: numSampl == 1
  while(numSampl > 1) {
    numSampl /= 2;
    prove *= 2;
  }
  if(prove == noSamples()) {
    m_FastFourier = true;
    m_log << " >> NOTE: Number of samples is a potency of 2. Thus Fast Fourier Transform can "
             "be used, if not defined otherwise (transformationType = 0 for DFT)."
          << std::endl;
  } else { // prove != noSamples()
    m_FastFourier = false;
    m_transformationType = 0;
    m_log << " >> NOTE: Number of samples is not a potency of 2. Thus Fast Fourier Transform "
             "cannot be used, no matter what."
          << std::endl;
  }
}

// Only works with analytical test cases
// TODO labels:ACA,cleanup check/cleanup
template <MInt nDim>
void AcaSolver<nDim>::computeAccuracy() {
  TERMM(1, "TODO check this function");
  if(domainId() == 0) {
    std::cout << "8. Starting accuarcy calculation." << std::endl;

    // These parameters have to be set
    // -------------------
    const MInt observer = globalNoObservers();
    MInt noSegs = 100;
    noSegs = Context::getSolverProperty<MFloat>("noSegments_" + std::to_string(0), m_solverId, AT_, 0);
    if(noSegs == 0) {
      noSegs = Context::getSolverProperty<MFloat>("noSegments_" + std::to_string(0), m_solverId, AT_, 1);
    }

    MFloat box_size = 2.0 * 2;
    box_size = Context::getSolverProperty<MFloat>("surfaceCoordinates_" + std::to_string(0), m_solverId, AT_, 0);
    box_size = 2.0 * std::fabs(box_size);
    // --------------------

    // Alternative expression: SamplePoints per Period even with flow
    MFloat noPeriods = 1.0;
    noPeriods = Context::getSolverProperty<MFloat>("noPeriods", m_solverId, AT_, &noPeriods);

    MFloat omega_factor = 4.0;
    omega_factor = Context::getSolverProperty<MFloat>("omega_factor", m_solverId, AT_, &omega_factor);

    MFloat Delta_Seg = box_size / noSegs;
    if constexpr(nDim == 3) {
      Delta_Seg *= Delta_Seg;
    }

    // Wave length lambda
    m_lambdaZero = 2.0 / m_omega_factor;
    m_lambdaMach = 2.0 * (1.0 - m_Ma) / m_omega_factor;

    // minimal resolution of one period according to the thesis.
    const MFloat R_min = noSamples() * m_lambdaMach / (noPeriods * m_lambdaZero);
    const MFloat Z = Delta_Seg / R_min;
    const MFloat Y = Delta_Seg / m_lambdaMach;

    // Read in the analytic solution for even periods to be able to compare with uneven periods
    std::vector<MFloat> Analytic_EvenPeriods(globalNoObservers() * noSamples() / 2 * 5, 0.0);
    std::vector<MFloat> Analytic_EvenPeriods_Amplitude(globalNoObservers() * noSamples() / 2, 0.0);
    std::vector<MFloat> Analytic_EvenPeriods_Frequencies(globalNoObservers() * noSamples() / 2, 0.0);
    std::vector<MFloat> RMS_EvenPeriods(globalNoObservers(), 0.0);

    const MInt EP_Samples = 64;
    if(noSamples() == EP_Samples) {
      std::ifstream file;
      MString line;
      MString Name = "./test/AnalyticSolution_ForEvenPeriods.dat";
      file.open(Name, std::ios::in);
      MInt count = 0;
      MFloat i1, i2, i3, i4, i5;
      if(file.is_open()) {
        while(getline(file, line)) {
          std::stringstream linebuffer(line);
          linebuffer >> i1 >> i2 >> i3 >> i4 >> i5;
          Analytic_EvenPeriods[5 * count] = i1;
          Analytic_EvenPeriods[5 * count + 1] = i2;
          Analytic_EvenPeriods[5 * count + 2] = i3;
          Analytic_EvenPeriods[5 * count + 3] = i4;
          Analytic_EvenPeriods[5 * count + 4] = i5;

          Analytic_EvenPeriods_Amplitude[count] = i4;
          Analytic_EvenPeriods_Frequencies[count] = i1;
          count += 1;
        }
      } else {
        std::cerr << "Fehler beim Öffnen der Datei: " << Name << "." << std::endl;
      }

      std::ifstream file1;
      MString line1;
      MString Name1 = "./test/RMSEvenPeriods.dat";
      file1.open(Name1, std::ios::in);
      MInt count1 = 0;
      MFloat i_1, i_2;
      if(file1.is_open()) {
        while(getline(file1, line1)) {
          std::stringstream linebuffer(line1);
          linebuffer >> i_1 >> i_2;
          RMS_EvenPeriods[count1] = i_2;
          count1 += 1;
        }
      } else {
        std::cerr << "Fehler beim Öffnen der Datei: " << Name1 << "." << std::endl;
      }
    }
    // RMS PRESSURE
    std::vector<MFloat> p_rms(globalNoObservers(), 0.0);

    for(MInt obs = 0; obs < globalNoObservers(); obs++) {
      for(MInt t = 0; t < noSamples(); t++) {
        p_rms[obs] += m_observerData.variables(obs, 0, t) * m_observerData.variables(obs, 0, t);
      }
      p_rms[obs] = sqrt(p_rms[obs] / noSamples());
    }

    // Quality (Variance) check
    MFloat max_diff = 0.0;
    MFloat max_diff_2 = 0.0;
    MFloat diff_normalized = 0.0;
    MFloat diff_2_normalized = 0.0;
    MFloat sum_diff_sq = 0.0;
    MFloat sum_diff_2_sq = 0.0;
    MFloat sum_analytic_sq = 0.0;
    MFloat sum_a_ep_sq = 0.0;
    std::vector<MFloat> diff(observer, 0.0);
    std::vector<MFloat> diff_2(observer, 0.0);
    for(MInt i = 0; i < observer; i++) {
      diff[i] = std::fabs(m_rmsP_Analytic[i] - p_rms[i]);
      sum_diff_sq += diff[i] * diff[i];
      sum_analytic_sq += m_rmsP_Analytic[i];
      diff_normalized = diff[i] / m_rmsP_Analytic[i] * 100;

      // In case of uneven periods
      if(noSamples() == EP_Samples) {
        diff_2[i] = std::fabs(RMS_EvenPeriods[i] - p_rms[i]);
        sum_diff_2_sq += diff_2[i] * diff_2[i];
        sum_a_ep_sq += RMS_EvenPeriods[i];
        diff_2_normalized = diff_2[i] / RMS_EvenPeriods[i] * 100;
      }
      if(diff_normalized >= max_diff) {
        if(m_rmsP_Analytic[i] > 10e-13) {
          max_diff = diff_normalized; // [%]
        }
        if(EP_Samples == 64) {
          max_diff_2 = diff_2_normalized; //[%]
        }
      }
    }
    const MFloat value_top = sqrt(sum_diff_sq / globalNoObservers());
    const MFloat value_bot = sum_analytic_sq / globalNoObservers();
    const MFloat Q_C = 100 * value_top / value_bot; // [%]

    // In case of uneven periods
    MFloat Q_C_UnevenPeriods = 0.0;
    if(noSamples() == EP_Samples) {
      const MFloat value_top_UP = sqrt(sum_diff_2_sq / globalNoObservers());
      const MFloat value_bot_UP = sum_a_ep_sq / globalNoObservers();
      Q_C_UnevenPeriods = 100 * value_top_UP / value_bot_UP;
    }

    // Frequency, Amplitude comparison for all observers
    std::vector<MFloat> diff_Ampl(globalNoObservers() * noSamples(), 0.0);
    std::vector<MFloat> max_Ampl(globalNoObservers(), 0.0);
    std::vector<MFloat> max_Ampl_2(globalNoObservers(), 0.0);
    std::vector<MFloat> max_Ampl_a(globalNoObservers(), 0.0);
    std::vector<MFloat> max_Ampl_A_EvenPeriods(globalNoObservers(), 0.0);
    std::vector<MFloat> max_Ampl_Diff(globalNoObservers(), 0.0);
    std::vector<MFloat> max_Ampl_Diff_2(globalNoObservers(), 0.0);
    std::vector<MFloat> max_Ampl_Diff_normalized(globalNoObservers(), 0.0);
    std::vector<MFloat> max_Ampl_Diff_normalized_2(globalNoObservers(), 0.0);
    std::vector<MFloat> freq_At_Max_Ampl(globalNoObservers(), 0.0);
    std::vector<MFloat> freq_At_Max_Ampl_2(globalNoObservers(), 0.0);
    std::vector<MFloat> freq_At_Max_Ampl_a(globalNoObservers(), 0.0);
    std::vector<MFloat> freq_At_Max_Ampl_A_EvenPeriods(globalNoObservers(), 0.0);
    std::vector<MFloat> freq_Diff(globalNoObservers(), 0.0);
    std::vector<MFloat> freq_Diff_2(globalNoObservers(), 0.0);
    std::vector<MFloat> freq_Diff_normalized(globalNoObservers(), 0.0);
    std::vector<MFloat> freq_Diff_normalized_2(globalNoObservers(), 0.0);
    std::vector<MFloat> mean_Ampl_analytic(globalNoObservers(), 0.0);
    std::vector<MFloat> sum_Diff_Amplitudes_UnevenPeriods(globalNoObservers(), 0.0);
    std::vector<MFloat> SideLobes(globalNoObservers(), 0.0);
    MFloat sum_max_Ampl_Diff_normalized = 0.0;
    MFloat sum_freq_Diff_normalized = 0.0;
    MFloat sum_max_Ampl_Diff_normalized_2 = 0.0;
    MFloat sum_freq_Diff_normalized_2 = 0.0;
    MFloat sum_SideLobes = 0.0;
    MFloat Amplitude_EvenPeriods = 0.0;
    for(MInt obs = 0; obs < globalNoObservers(); obs++) {
      for(MInt i = 1; i < noSamples() / 2 + 1; i++) {
        const MFloat real_c = m_observerData.complexVariables(obs, 0, i, 0);
        const MFloat imag_c = m_observerData.complexVariables(obs, 0, i, 1);
        const MFloat real_a = m_Analyticfreq[obs * noSamples() * 2 + 2 * i];
        const MFloat imag_a = m_Analyticfreq[obs * noSamples() * 2 + 2 * i + 1];
        const MFloat Amplitude_c = 2 * sqrt(real_c * real_c + imag_c * imag_c);
        const MFloat Amplitude_a = 2 * sqrt(real_a * real_a + imag_a * imag_a);
        // Comparison with even periods in case uneven periods are used
        if(noSamples() == EP_Samples) {
          Amplitude_EvenPeriods = Analytic_EvenPeriods_Amplitude[obs * noSamples() / 2 + (i - 1)];

          // Get Frequency at which the analytic source for even periods oscillates
          if(Amplitude_EvenPeriods > max_Ampl_A_EvenPeriods[obs]) {
            freq_At_Max_Ampl_A_EvenPeriods[obs] = Analytic_EvenPeriods_Frequencies[obs * noSamples() / 2 + (i - 1)];
            max_Ampl_A_EvenPeriods[obs] = Amplitude_EvenPeriods;
          }
        }

        // Get Frequency at which analytic source oscillates
        if(Amplitude_a > max_Ampl_a[obs]) {
          freq_At_Max_Ampl_a[obs] = m_frequencies[i];
          max_Ampl_a[obs] = Amplitude_a;
        }

        // Get max Amplitude computed and its frequency
        if(Amplitude_c > max_Ampl[obs]) {
          freq_At_Max_Ampl[obs] = m_frequencies[i];
          max_Ampl[obs] = Amplitude_c;
        }

        // Side lobes
        SideLobes[obs] += std::fabs(Amplitude_EvenPeriods - Amplitude_c);
      }
      // Normalize the maximum difference in amplitude and frequency with the maximum analytical
      // ampl. and freq. Index 2 is for uneven Periods comparison
      max_Ampl_Diff[obs] = std::fabs(max_Ampl_a[obs] - max_Ampl[obs]);
      freq_Diff[obs] = std::fabs(freq_At_Max_Ampl[obs] - freq_At_Max_Ampl_a[obs]);

      max_Ampl_Diff_normalized[obs] = max_Ampl_Diff[obs] / max_Ampl_a[obs];
      freq_Diff_normalized[obs] = freq_Diff[obs] / freq_At_Max_Ampl_a[obs];

      sum_max_Ampl_Diff_normalized += max_Ampl_Diff_normalized[obs];
      sum_freq_Diff_normalized += freq_Diff_normalized[obs];

      if(noSamples() == EP_Samples) {
        max_Ampl_Diff_2[obs] = std::fabs(max_Ampl_A_EvenPeriods[obs] - max_Ampl[obs]);
        freq_Diff_2[obs] = std::fabs(freq_At_Max_Ampl[obs] - freq_At_Max_Ampl_A_EvenPeriods[obs]);
        sum_SideLobes += (SideLobes[obs] / (noSamples() / 2));

        max_Ampl_Diff_normalized_2[obs] = max_Ampl_Diff_2[obs] / max_Ampl_A_EvenPeriods[obs];
        freq_Diff_normalized_2[obs] = freq_Diff_2[obs] / freq_At_Max_Ampl_A_EvenPeriods[obs];

        // In case of uneven and even periods comparison
        sum_max_Ampl_Diff_normalized_2 += max_Ampl_Diff_normalized_2[obs];
        sum_freq_Diff_normalized_2 += freq_Diff_normalized_2[obs];
      }
    }
    const MFloat Ampl_aver_Diff = 100 * sum_max_Ampl_Diff_normalized / globalNoObservers(); // [%]
    const MFloat Freq_aver_Diff = 100 * sum_freq_Diff_normalized / globalNoObservers();     // [%]

    // In case of uneven and even periods comparison
    MFloat Side_Lobes = 0.0;
    MFloat Freq_aver_Diff_UnevenPeriods = 0.0;
    MFloat Ampl_aver_Diff_UnevenPeriods = 0.0;
    if(noSamples() == EP_Samples) {
      Ampl_aver_Diff_UnevenPeriods = 100 * sum_max_Ampl_Diff_normalized_2 / globalNoObservers(); //[%]
      Freq_aver_Diff_UnevenPeriods = 100 * sum_freq_Diff_normalized_2 / globalNoObservers();     //[%]
      Side_Lobes = 100 * sum_SideLobes / globalNoObservers();                                    //[%]
      // Those minor deviations are due to the writing out in a file and re-reading it again
      // (Analytic Even Periods). Thus, they match perfectly but a deviation is falsly shown
      if(Freq_aver_Diff_UnevenPeriods < 1e-10) {
        Freq_aver_Diff_UnevenPeriods = 0.0;
      }
    }

    // OUTPUT OF max Amplitude at given frequencie
    MString fileFreq = "./FreqComparison/maxFreqCompare_";
    fileFreq += std::to_string(noSegs);
    fileFreq += ".dat";
    std::ofstream outfile;
    outfile.open(fileFreq);
    outfile << Y << " " << Freq_aver_Diff << " " << Ampl_aver_Diff << " " << Side_Lobes << std::endl;
    outfile.close();

    // Freq for samples
    MString fileF = "./FreqComparison/maxFreqCompareSamples_";
    fileF += std::to_string(noSamples());
    fileF += ".dat";
    outfile.open(fileF);
    outfile << R_min << " " << Freq_aver_Diff << " " << Ampl_aver_Diff << " " << Side_Lobes << std::endl;

    outfile.close();

    if(noSamples() == EP_Samples) {
      // Freq for periods
      MString fileFr = "./FreqComparison/maxFreqComparePeriods_";
      fileFr += std::to_string(noPeriods);
      fileFr += ".dat";
      outfile.open(fileFr);
      outfile << noPeriods << " " << Freq_aver_Diff_UnevenPeriods << " " << Ampl_aver_Diff_UnevenPeriods << " "
              << Side_Lobes << " " << Ampl_aver_Diff << std::endl;
      outfile.close();
    }

    // OUTPUT IF ONE VARIES NUMBER OF SEGMENTS
    MString ofile = "./Diff_seg/accuracy_";
    ofile += std::to_string(noSegs);
    ofile += ".dat";
    outfile.open(ofile);
    outfile << Y << " " << Q_C << " " << max_diff << " " << Z << " " << Side_Lobes << std::endl;
    outfile.close();

    // OUTPUT IF ONE VARIES NUMBER OF SAMPLES
    MString fileName = "./Diff_Samples/accuracy_";
    fileName += std::to_string(noSamples());
    fileName += ".dat";
    outfile.open(fileName);
    outfile << R_min << " " << Q_C << " " << max_diff << " " << Side_Lobes << std::endl;
    outfile.close();

    if(noSamples() == EP_Samples) {
      // OUTPUT IF ONE VARIES NUMBER OF SAMPLES
      MString fileN = "./Diff_Periods/accuracy_";
      fileN += std::to_string(noPeriods);
      fileN += ".dat";
      outfile.open(fileN);
      outfile << noPeriods << " " << Q_C_UnevenPeriods << " " << max_diff_2 << " " << Side_Lobes << std::endl;
      outfile.close();
    }

    // OUTPUT IF ONE VARIES NUMBER OF SAMPLES
    MString fileNa = "./Diff_omega/accuracy_";
    fileNa += std::to_string(omega_factor);
    fileNa += ".dat";
    outfile.open(fileNa);
    outfile << omega_factor << " " << Q_C << " " << max_diff << " " << Side_Lobes << std::endl;
    outfile.close();

    std::cout << "8. Accuarcy calculation finished." << std::endl;
  }
}


template class AcaSolver<2>;
template class AcaSolver<3>;
