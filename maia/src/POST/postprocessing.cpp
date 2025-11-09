// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "postprocessing.h"
#include <cmath>
#include "COMM/mpioverride.h"

#include "globals.h"
// Required for FV instantiation:
#include "FV/fvcartesiansolverxd.h"
// Required for LB instantiation:
#if not defined(MAIA_DISABLE_LB)
#include "LB/lbsolver.h"
#endif
// Required for DG instantiation:
#if not defined(MAIA_DISABLE_DG)
#include "DG/dgcartesiansolver.h"
#include "DG/dgcartesiansyseqnacousticperturb.h"
#include "DG/dgcartesiansyseqnlinearscalaradv.h"
#endif
#include "UTIL/parallelfor.h"

#include "IO/parallelio.h"

template <MInt nDim, class ppType>
const constexpr MInt PostProcessing<nDim, ppType>::MV::LAMB0;
template <MInt nDim, class ppType>
const constexpr MInt PostProcessing<nDim, ppType>::MV::VORT0;
template <MInt nDim, class ppType>
const constexpr MInt PostProcessing<nDim, ppType>::MV::DU;
template <MInt nDim, class ppType>
const constexpr MInt PostProcessing<nDim, ppType>::MV::DRHO;
template <MInt nDim, class ppType>
const constexpr MInt PostProcessing<nDim, ppType>::MV::DP;
template <MInt nDim, class ppType>
const constexpr MInt PostProcessing<nDim, ppType>::MV::RHODIVU;
template <MInt nDim, class ppType>
const constexpr MInt PostProcessing<nDim, ppType>::MV::UGRADRHO;
template <MInt nDim, class ppType>
const constexpr MInt PostProcessing<nDim, ppType>::MV::GRADPRHO;
template <MInt nDim, class ppType>
const constexpr MInt PostProcessing<nDim, ppType>::MV::GRADU;
template <MInt nDim, class ppType>
const constexpr MInt PostProcessing<nDim, ppType>::MV::UGRADU;

template <MInt nDim>
class PostProcessingLb;

using namespace std;


/** \brief Destructor for the massive paralle postprocessing solver
 *
 * \author Andreas Lintermann
 * \date 26.08.2012
 *
 * \tparam[in] T celltype
 *
 **/
template <MInt nDim, class ppType>
PostProcessing<nDim, ppType>::~PostProcessing() {
  TRACE();

  if(m_noPostprocessingOps > 0)
    for(MInt op = 0; op < m_noPostprocessingOps; op++) {
      switch(string2enum(m_postprocessingOps[op])) {
        case PP_PROBE_POINT_PRE:
        case PP_PROBE_POINT_POST:
        case PP_PROBE_POINT_IN: {
          for(MInt np = 0; np < m_noProbePoints; np++)
            m_probeFileStreams[np].close();

          delete[] m_probeFileStreams;
          break;
        }
        case PP_AVERAGE_PRE:
        case PP_AVERAGE_IN:
        case PP_AVERAGE_POST: {
          // mDeallocate(m_summedVars);
          // mDeallocate(m_square);
          if(m_useKahan) {
            // mDeallocate(m_cSum);
            // mDeallocate(m_ySum);
            // mDeallocate(m_tSum);

            // mDeallocate(m_cSquare);
            // mDeallocate(m_ySquare);
            // mDeallocate(m_tSquare);

            if(m_kurtosis) {
              // mDeallocate(m_cCube);
              // mDeallocate(m_yCube);
              // mDeallocate(m_tCube);

              // mDeallocate(m_cFourth);
              // mDeallocate(m_yFourth);
              // mDeallocate(m_tFourth);
            } else if(m_skewness) {
              // mDeallocate(m_cCube);
              // mDeallocate(m_yCube);
              // mDeallocate(m_tCube);
            }
          }
          if(m_kurtosis) {
            // mDeallocate(m_cube);
            // mDeallocate(m_fourth);
          } else if(m_skewness) {
            // mDeallocate(m_cube);
          }
          break;
        }
        case PP_MOVING_AVERAGE_IN:
        case PP_MOVING_AVERAGE_PRE:
        case PP_MOVING_AVERAGE_POST: {
          break;
        }
        case PP_SPATIAL_AVERAGE_PRE:
        case PP_SPATIAL_AVERAGE_POST:
        case PP_SPATIAL_AVERAGE_IN: {
          break;
        }
        case PP_PROBE_LINE_PRE:
        case PP_PROBE_LINE_POST:
        case PP_PROBE_LINE_IN: {
          break;
        }
        case PP_PROBE_ARB_LINE_PRE:
        case PP_PROBE_ARB_LINE_POST:
        case PP_PROBE_ARB_LINE_IN: {
          break;
        }
        case PP_PROBE_SLICE_PRE:
        case PP_PROBE_SLICE_POST:
        case PP_PROBE_SLICE_IN: {
          break;
        }
        case PP_PROBE_ARB_SLICE_PRE:
        case PP_PROBE_ARB_SLICE_POST:
        case PP_PROBE_ARB_SLICE_IN: {
          break;
        }
        case PP_WRITEPOINTS_IN: {
          break;
        }
        case PP_REDUCE_TO_LEVEL_PRE:
        case PP_REDUCE_TO_LEVEL_POST: {
          break;
        }
        case PP_AVERAGE_SLICE_PRE: {
          break;
        }
        case PP_SPRAY_STATS:
        case PP_PARTICLE_STATISTICS:
        case PP_PARTICLE_SOLUTION: {
          break;
        }
        case PP_PL_ISO_TURBULENCE_STATISTICS:
        case PP_ISO_TURBULENCE_STATISTICS: {
          break;
        }
        default: {
          mTerm(1, AT_, "Unknown postprocessing operation");
        }
      }
    }
  delete[] m_postprocessingOps;
  //}
}


template <MInt nDim, class ppType>
void PostProcessing<nDim, ppType>::postprocessPreSolve() {
  TRACE();

  for(MInt op = 0; op < (signed)m_postprocessingMethods[0].size(); op++) {
    (this->*(m_postprocessingMethods[0][op]))();
  }
}


template <MInt nDim, class ppType>
void PostProcessing<nDim, ppType>::postprocessInSolve(MBool finalTimeStep) {
  TRACE();

  m_finalTimeStep = finalTimeStep;
  for(MInt op = 0; op < (signed)m_postprocessingMethods[1].size(); op++) {
    (this->*(m_postprocessingMethods[1][op]))();
  }
}


template <MInt nDim, class ppType>
void PostProcessing<nDim, ppType>::postprocessSolution() {
  TRACE();

  for(MInt op = 0; op < (signed)m_postprocessingSolution.size(); op++) {
    (this->*(m_postprocessingSolution[op]))();
  }
}


template <MInt nDim, class ppType>
void PostProcessing<nDim, ppType>::postprocessPostSolve() {
  TRACE();

  m_finalTimeStep = true;
  for(MInt op = 0; op < (signed)m_postprocessingMethods[2].size(); op++) {
    (this->*(m_postprocessingMethods[2][op]))();
  }
}


/**
 * \fn void PostProcessing<nDim,ppType>::initPostProcessing()
 * \brief Reads all required properties in and prepares for postprocessing
 *
 * \author Andreas Lintermann
 * \date 12.09.2012_" + to_string(m_postprocessingId)
 *
 **/
template <MInt nDim, class ppType>
void PostProcessing<nDim, ppType>::initPostProcessing() {
  TRACE();

  m_gridProxy = &pp->solver().grid();

  m_log << "  + Initializing postprocessing" << endl << endl;

  // Read & store solver type to allow multiple solvers to coexist
  m_solverType = Context::getSolverProperty<MString>("solvertype", m_postprocessingId, AT_);

  // global properties already required here
  /*! \page propertiesPP
    \section restartTimeStep
    <code>MInt PostprocesssingSolver::m_restartTimeStep</code>\n
    default = <code>0</code>\n\n
    This property determines the timestep to restart from. Is only active, if
    restartFile is set.
    <ul>
    <li><code>timestep</code> </li>
    </ul>\n
    Keywords: <i>GENERAL, GLOBAL</i>
  */
  m_restartTimeStep = Context::getSolverProperty<MInt>("restartTimeStep", m_postprocessingId, AT_);

  // Determine number of variables in solver
  m_noVariables = pp->solver().noVariables();

  /*! \page propertiesPP
    \section pp_fileName
    <code>MInt PostprocesssingSolver::m_postprocessFileName</code>\n
    default = <code>""</code>\n\n
    Provide a file that is used for postprocessing (e.g. spatial averaging or probing).\n
    <ul>
    <li><code>filename</code></li>
    </ul>\n
    Keywords: <i>POSTPROCESSING</i>
  */
  m_postprocessFileName = "";
  m_postprocessFileName =
      Context::getSolverProperty<MString>("pp_fileName", m_postprocessingId, AT_, &m_postprocessFileName);

  m_noPostprocessingOps = 0;

  m_log << "  + Postprocessing is active" << endl;
  for(MInt i = 0; i < 3; i++) {
    tvecpost tmp;
    vector<MString> st;
    m_postprocessingMethods.push_back(tmp);
    m_postprocessingMethodsDesc.push_back(st);
  }

  if(Context::propertyExists("apeSourceTermsAverages", m_postprocessingId)) {
    /*! \page propertiesPP
      \section apeSourceTermsAverages
      <code>vector<MInt> PostProcessing::m_activeSourceTerms</code>\n
      default = <code>none</code>\n \n
      Names of the coupling source terms to use.\n
      The possible source term names are stored in s_sourceTermNames.\n
      Keywords: <i>COUPLING, PHYSICS, SOURCE_TERM</i>
    */
    // Get the number of source terms
    const MInt noSourceTerms = Context::propertyLength("apeSourceTermsAverages", m_postprocessingId);

    // Ensure that at least one source term is loaded
    if(noSourceTerms == 0) {
      TERMM(1, "No APE source term average was specified.");
    }

    // Loop over all given source terms
    for(MInt i = 0; i < noSourceTerms; i++) {
      const MString sourceTermName =
          Context::getSolverProperty<MString>("apeSourceTermsAverages", m_postprocessingId, AT_, i);


      // Find the source term name in the list of all source terms
      auto sourceTermIt = std::find(s_sourceTermNames.begin(), s_sourceTermNames.end(), sourceTermName);

      // Exit if source term name not found
      // if(postData().getVariablePropertyIndex(sourceTermName) = -1){
      if(sourceTermIt == s_sourceTermNames.end()) {
        TERMM(1, "The given APE source term average '" + sourceTermName + "' was not found.");
      }

      // Determine source term index
      const MInt sourceTerm = std::distance(s_sourceTermNames.begin(), sourceTermIt);

      // Check if this source term is already active
      if(std::find(m_activeSourceTerms.begin(), m_activeSourceTerms.end(), sourceTerm) != m_activeSourceTerms.end()) {
        TERMM(1, "Given APE source term averages '" + sourceTermName + "' already present. Check your property file.");
      }

      // Valid source term, store in list of active source terms
      m_activeSourceTerms.push_back(sourceTerm);
    }
    // Sort all active source terms by id
    std::sort(m_activeSourceTerms.begin(), m_activeSourceTerms.end());

    // Determine unique list of needed mean variables for all active source terms
    // and store in m_activeMeanVars
    for(auto&& sourceTerm : m_activeSourceTerms) {
      std::vector<MInt> sourceTermMeanVars;
      // Get the needed mean variables for the current source term
      neededMeanVarsForSourceTerm(sourceTerm, sourceTermMeanVars);

      // Insert into set of needed mean variables
      m_activeMeanVars.insert(sourceTermMeanVars.begin(), sourceTermMeanVars.end());
    }
  }

  /*! \page ugPostprocessing
    \section postprocessingOp
    <code>MString* PostProcessing::m_postprocessingOps</code>\n
    default = <code>empty</code>\n\n
    This property is a list of postprocessing operations to be performed
    <ul>
    <li><code>PP_AVERAGE_PRE</code> </li>
    <li><code>PP_AVERAGE_IN</code> </li>
    <li><code>PP_AVERAGE_POST</code> </li>
    <li><code>PP_COMPUTE_DIVERGENCEVELOCITY_PRE</code> </li>
    <li><code>PP_COMPUTE_DIVERGENCEVELOCITY_IN</code> </li>
    <li><code>PP_COMPUTE_DIVERGENCEVELOCITY_POST</code> </li>
    <li><code>PP_MOVING_AVERAGE_PRE</code> </li>
    <li><code>PP_MOVING_AVERAGE_IN</code> </li>
    <li><code>PP_MOVING_AVERAGE_POST</code> </li>
    <li><code>PP_PROBE_POINT_PRE</code> </li>
    <li><code>PP_PROBE_POINT_IN</code> </li>
    <li><code>PP_PROBE_POINT_POST</code> </li>
    <li><code>PP_PROBE_LINE_PRE</code> </li>
    <li><code>PP_PROBE_LINE_IN</code> </li>
    <li><code>PP_PROBE_LINE_POST</code> </li>
    <li><code>PP_PROBE_LINE_PERIODIC_IN</code> </li>
    <li><code>PP_PROBE_LINE_PERIODIC_POST</code> </li>
    <li><code>PP_PROBE_ARB_LINE_PRE</code>  </li>
    <li><code>PP_PROBE_ARB_LINE_IN</code>  </li>
    <li><code>PP_PROBE_ARB_LINE_POST</code>  </li>
    <li><code>PP_PROBE_SLICE_PRE</code> </li>
    <li><code>PP_PROBE_SLICE_POST</code> </li>
    <li><code>PP_PROBE_SLICE_IN</code> </li>
    <li><code>PP_PROBE_ARB_SLICE_PRE</code> </li>
    <li><code>PP_PROBE_ARB_SLICE_IN</code> </li>
    <li><code>PP_PROBE_ARB_SLICE_POST</code> </li>
    <li><code>PP_AVERAGE_SLICE_PRE</code> </li>
    <li><code>PP_REDUCE_TO_LEVEL_PRE</code> </li>
    <li><code>PP_REDUCE_TO_LEVEL_POST</code> </li>
    <li><code>PP_REDUCE_TO_LEVEL_AVERAGES_PRE</code></li>
    <li><code>PP_SPATIAL_AVERAGE_PRE</code> </li>
    <li><code>PP_SPATIAL_AVERAGE_POST</code> </li>
    <li><code>PP_SPATIAL_AVERAGE_IN</code> </li>
    <li><code>PP_AVERAGE_REYNOLDS_STRESSES_AND_AVERAGES_IN</code> </li>
    </ul>\n
    Keywords: <i>GENERAL, GLOBAL, POSTPROCESSING</i>
  */
  m_noPostprocessingOps = Context::propertyLength("postProcessingOps_" + to_string(m_postprocessingId));

  m_log << "    - number of operations: " << m_noPostprocessingOps << endl;
  m_log << "    - requested operations: " << endl;
  m_postprocessingOps = nullptr;

  if(m_noPostprocessingOps > 0) {
    m_postprocessingOps = new MString[m_noPostprocessingOps];

    for(MInt op = 0; op < m_noPostprocessingOps; op++) {
      m_postprocessingOps[op] = Context::getBasicProperty<MString>("postProcessingOps_" + to_string(m_postprocessingId),
                                                                   AT_, &m_postprocessingOps[op], op);
      m_log << "      * " << m_postprocessingOps[op] << " (enum: " << string2enum(m_postprocessingOps[op]) << ")"
            << endl;

      switch(string2enum(m_postprocessingOps[op])) {
        case PP_AVERAGE_PRE: {
          m_postprocessingMethodsDesc[0].push_back(m_postprocessingOps[op]);
          m_postprocessingMethods[0].push_back(&PostProcessing::averageSolutions);
          initTimeStepProperties();
          initAveragingProperties();
          initAverageVariables();
          break;
        }
        case PP_AVERAGE_IN: {
          m_postprocessingMethodsDesc[1].push_back(m_postprocessingOps[op]);
          m_postprocessingMethods[1].push_back(&PostProcessing::averageSolutionsInSolve);
          m_postprocessingSolution.push_back(&PostProcessing::computeAndSaveMean);

          // Add call to average-in-solve to preSolve pp-calls such that at a restart with
          // restartTimeStep == averageStartTimestep that time step is also taken into account!
          m_postprocessingMethodsDesc[0].push_back(m_postprocessingOps[op]);
          m_postprocessingMethods[0].push_back(&PostProcessing::averageSolutionsInSolve);
          initTimeStepProperties();
          initAveragingProperties();
          if(m_correlation) {
            initCorrelation();
          }
          initAverageVariables();
          break;
        }
        case PP_AVERAGE_POST: {
          m_postprocessingMethodsDesc[2].push_back(m_postprocessingOps[op]);
          m_postprocessingMethods[2].push_back(&PostProcessing::averageSolutions);
          initTimeStepProperties();
          initAveragingProperties();
          initAverageVariables();
          break;
        }
        case PP_COMPUTE_DIVERGENCEVELOCITY_PRE: {
          m_postprocessingMethodsDesc[0].push_back(m_postprocessingOps[op]);
          m_postprocessingMethods[0].push_back(&PostProcessing::computeAndSaveDivergence<false>);
          break;
        }
        case PP_COMPUTE_DIVERGENCEVELOCITY_IN: {
          m_postprocessingMethodsDesc[1].push_back(m_postprocessingOps[op]);
          m_postprocessingMethods[1].push_back(&PostProcessing::computeAndSaveDivergence<true>);
          break;
        }
        case PP_COMPUTE_DIVERGENCEVELOCITY_POST: {
          m_postprocessingMethodsDesc[2].push_back(m_postprocessingOps[op]);
          m_postprocessingMethods[2].push_back(&PostProcessing::computeAndSaveDivergence<false>);
          break;
        }
        case PP_MOVING_AVERAGE_PRE: {
          m_postprocessingMethodsDesc[0].push_back(m_postprocessingOps[op]);
          m_postprocessingMethods[0].push_back(&PostProcessing::movingAveragePost);
          initAveragingProperties();
          initMovingAverage();
          break;
        }
        case PP_MOVING_AVERAGE_IN: {
          m_postprocessingMethodsDesc[1].push_back(m_postprocessingOps[op]);
          m_postprocessingMethods[1].push_back(&PostProcessing::movingAverage);
          initAveragingProperties();
          initMovingAverage();
          break;
        }
        case PP_MOVING_AVERAGE_POST: {
          m_postprocessingMethodsDesc[2].push_back(m_postprocessingOps[op]);
          m_postprocessingMethods[2].push_back(&PostProcessing::movingAveragePost);
          initAveragingProperties();
          initMovingAverage();
          break;
        }
        case PP_PROBE_POINT_PRE: {
          m_postprocessingMethodsDesc[0].push_back(m_postprocessingOps[op]);
          m_postprocessingMethods[0].push_back(&PostProcessing::probeLocation);
          initProbePoint();
          break;
        }
        case PP_PROBE_POINT_IN: {
          m_postprocessingMethodsDesc[1].push_back(m_postprocessingOps[op]);
          m_postprocessingMethods[1].push_back(&PostProcessing::probeLocation);
          initProbePoint();

          break;
        }
        case PP_PROBE_POINT_POST: {
          m_postprocessingMethodsDesc[2].push_back(m_postprocessingOps[op]);
          m_postprocessingMethods[2].push_back(&PostProcessing::probeLocation);
          initProbePoint();
          break;
        }
        case PP_PROBE_LINE_PRE: {
          m_postprocessingMethodsDesc[0].push_back(m_postprocessingOps[op]);
          m_postprocessingMethods[0].push_back(&PostProcessing::probeLinePre);
          initProbeLine();
          initTimeStepPropertiesLine();
          break;
        }
        case PP_PROBE_LINE_IN: {
          m_postprocessingMethodsDesc[1].push_back(m_postprocessingOps[op]);
          m_postprocessingMethods[1].push_back(&PostProcessing::probeLine);
          initProbeLine();
          initTimeStepPropertiesLine();
          break;
        }
        case PP_PROBE_LINE_POST: {
          m_postprocessingMethodsDesc[2].push_back(m_postprocessingOps[op]);
          m_postprocessingMethods[2].push_back(&PostProcessing::probeLinePost);
          initAverageVariables();
          initProbeLine();
          initTimeStepPropertiesLine();
          initAveragingProperties();
          // initTimeStepProperties();
          break;
        }
        case PP_PROBE_LINE_PERIODIC_IN: {
          m_postprocessingMethodsDesc[1].push_back(m_postprocessingOps[op]);
          m_postprocessingMethods[1].push_back(&PostProcessing::probeLinePeriodic);
          initProbeLine();
          initProbeLinePeriodic();
          initTimeStepPropertiesLine();
          break;
        }
        case PP_PROBE_LINE_PERIODIC_POST: {
          m_postprocessingMethodsDesc[2].push_back(m_postprocessingOps[op]);
          m_postprocessingMethods[2].push_back(&PostProcessing::probeLinePeriodicPost);
          initTimeStepPropertiesLine();
          initAveragingProperties();
          if(m_correlation) {
            initCorrelation();
          }
          initAverageVariables();
          initProbeLine();
          initProbeLinePeriodic();
          break;
        }
        case PP_SLICE_AVERAGE: {
          m_postprocessingMethodsDesc[1].push_back(m_postprocessingOps[op]);
          m_postprocessingMethods[1].push_back(&PostProcessing::periodicSliceAverage);
          initTimeStepProperties();
          initAveragingProperties();
          if(m_correlation) {
            initCorrelation();
          }
          initAverageVariables();
          initPeriodicSliceAverage();
          break;
        }
        case PP_PROBE_ARB_LINE_PRE: {
          m_postprocessingMethodsDesc[0].push_back(m_postprocessingOps[op]);
          m_postprocessingMethods[0].push_back(&PostProcessing::probeArbitraryLinePost);
          initProbeArbitraryLine();
          initTimeStepPropertiesLine();
          break;
        }
        case PP_PROBE_ARB_LINE_IN: {
          m_postprocessingMethodsDesc[1].push_back(m_postprocessingOps[op]);
          m_postprocessingMethods[1].push_back(&PostProcessing::probeArbitraryLine);
          initProbeArbitraryLine();
          initTimeStepPropertiesLine();
          break;
        }
        case PP_PROBE_ARB_LINE_POST: {
          m_postprocessingMethodsDesc[2].push_back(m_postprocessingOps[op]);
          m_postprocessingMethods[2].push_back(&PostProcessing::probeArbitraryLinePost);
          initAveragingProperties();
          initAverageVariables();
          initProbeArbitraryLine();
          initTimeStepPropertiesLine();
          break;
        }
        case PP_PROBE_SLICE_PRE: {
          m_postprocessingMethodsDesc[0].push_back(m_postprocessingOps[op]);
          m_postprocessingMethods[0].push_back(&PostProcessing::probeSlicePre);
          initProbeSlice();
          initTimeStepPropertiesSlice();
          break;
        }
        case PP_PROBE_SLICE_POST: {
          m_postprocessingMethodsDesc[2].push_back(m_postprocessingOps[op]);
          m_postprocessingMethods[2].push_back(&PostProcessing::probeSlicePost);
          initAveragingProperties();
          initAverageVariables();
          initProbeSlice();
          initTimeStepPropertiesSlice();
          break;
        }
        case PP_PROBE_SLICE_IN: {
          m_postprocessingMethodsDesc[1].push_back(m_postprocessingOps[op]);
          m_postprocessingMethods[1].push_back(&PostProcessing::probeSliceIn);
          initProbeSlice();
          initTimeStepPropertiesSlice();
          break;
        }
        case PP_PROBE_ARB_SLICE_PRE: {
          m_postprocessingMethodsDesc[0].push_back(m_postprocessingOps[op]);
          m_postprocessingMethods[0].push_back(&PostProcessing::probeArbitrarySlicePost);
          initProbeArbitrarySlice();
          initTimeStepPropertiesSlice();
          break;
        }
        case PP_PROBE_ARB_SLICE_IN: {
          m_postprocessingMethodsDesc[1].push_back(m_postprocessingOps[op]);
          m_postprocessingMethods[1].push_back(&PostProcessing::probeArbitrarySlice);
          initProbeArbitrarySlice();
          initTimeStepPropertiesSlice();
          break;
        }
        case PP_PROBE_ARB_SLICE_POST: {
          m_postprocessingMethodsDesc[2].push_back(m_postprocessingOps[op]);
          m_postprocessingMethods[2].push_back(&PostProcessing::probeArbitrarySlicePost);
          initProbeArbitrarySlice();
          initTimeStepPropertiesSlice();
          break;
        }
        case PP_AVERAGE_SLICE_PRE: {
          m_postprocessingMethodsDesc[0].push_back(m_postprocessingOps[op]);
          m_postprocessingMethods[0].push_back(&PostProcessing::averageSolutionsSlice);
          initTimeStepProperties();
          initTimeStepPropertiesSlice();
          break;
        }
        case PP_REDUCE_TO_LEVEL_PRE: {
          m_postprocessingMethodsDesc[0].push_back(m_postprocessingOps[op]);
          m_postprocessingMethods[0].push_back(&PostProcessing::pp_saveCoarseSolution);
          break;
        }
        case PP_REDUCE_TO_LEVEL_POST: {
          m_postprocessingMethodsDesc[2].push_back(m_postprocessingOps[op]);
          m_postprocessingMethods[2].push_back(&PostProcessing::pp_saveCoarseSolution);
          break;
        }
          // PP_ToDo-LB: pp_saveCoarseSolution has to be done in postdata
          //            delete m_localVars
        case PP_REDUCE_TO_LEVEL_AVERAGES_PRE: {
          m_postprocessingMethodsDesc[0].push_back(m_postprocessingOps[op]);
          m_postprocessingMethods[0].push_back(&PostProcessing::pp_saveCoarseSolution);
          /*this array now holds the averages and the entries of the Reynolds stress tensor and all required values

            thermal        non-thermal
            ---------------------------------
            0:  u_mean         u_mean
            1:  v_mean         v_mean
            2:  w_mean         w_mean
            3:  rho_mean       rho_mean
            4:  t_mean         u^2
            5:  u^2            uv
            6:  uv             uw
            7:  uw             v^2
            8:  v^2            vw
            9:  vw             w^2
            10: w^2

            5/4:  (u'u')_mean = (u^2)_mean - (u_mean)^2
            6/5:  (u'v')_mean = (uv)_mean - u_mean * v_mean
            7/6:  (u'w')_mean = (uw)_mean - u_mean * w_mean
            8/7:  (v'v')_mean = (v^2)_mean - (v_mean)^2
            9/8:  (v'w')_mean = (vw)_mean - v_mean * v_mean
            10/9: (w'w')_mean = (w^2)_mean - (w_mean)^2

          */

          m_noLocalVars = m_noVariables + 2 * nDim;
          mAlloc(m_localVars, pp->solver().grid().noInternalCells(), m_noLocalVars, "m_localVars", F0, AT_);
          initReduceToLevelAvg();
          break;
        }
        case PP_SPATIAL_AVERAGE_PRE: {
          m_postprocessingMethodsDesc[0].push_back(m_postprocessingOps[op]);
          m_postprocessingMethods[0].push_back(&PostProcessing::spatialAveragingPost);
          initTimeStepProperties();
          initSpatialAveraging();
          break;
        }
        case PP_SPATIAL_AVERAGE_POST: {
          m_postprocessingMethodsDesc[2].push_back(m_postprocessingOps[op]);
          m_postprocessingMethods[2].push_back(&PostProcessing::spatialAveragingPost);
          initTimeStepProperties();
          initSpatialAveraging();
          break;
        }
        case PP_SPATIAL_AVERAGE_IN: {
          m_postprocessingMethodsDesc[1].push_back(m_postprocessingOps[op]);
          m_postprocessingMethods[1].push_back(&PostProcessing::spatialAveraging);
          initTimeStepProperties();
          initSpatialAveraging();
          break;
        }
        case PP_WRITEPOINTS_IN: {
          m_postprocessingMethodsDesc[1].push_back(m_postprocessingOps[op]);
          m_postprocessingMethods[1].push_back(&PostProcessing::writePointData);
          initWritePointData();
          break;
        }
        case PP_SPRAY_STATS: {
          m_postprocessingMethodsDesc[1].push_back(m_postprocessingOps[op]);
          m_postprocessingMethods[1].push_back(&PostProcessing::computeSprayData);
          m_postprocessingMethods[2].push_back(&PostProcessing::writeSprayData);
          initSprayData();
          break;
        }
        case PP_PARTICLE_SOLUTION: {
          m_postprocessingMethodsDesc[1].push_back(m_postprocessingOps[op]);
          m_postprocessingMethods[1].push_back(&PostProcessing::writeLPTSolutionFile);
          initLPTSolutionFile();
          break;
        }
        case PP_POINT_SAMPLING_IN: {
          m_postprocessingMethodsDesc[1].push_back(m_postprocessingOps[op]);
          m_postprocessingMethods[1].push_back(&PostProcessing::savePointSamplingData);
          initPointSamplingData();
          break;
        }
        case PP_SURFACE_SAMPLING_IN: {
          m_postprocessingMethodsDesc[1].push_back(m_postprocessingOps[op]);
          m_postprocessingMethods[1].push_back(&PostProcessing::saveSurfaceSamplingData);
          initSurfaceSamplingData();
          break;
        }
        case PP_VOLUME_SAMPLING_IN: {
          m_postprocessingMethodsDesc[1].push_back(m_postprocessingOps[op]);
          m_postprocessingMethods[1].push_back(&PostProcessing::saveVolumeSamplingData);
          initVolumeSamplingData();
          break;
        }
        case PP_PARTICLE_STATISTICS: {
          m_postprocessingMethodsDesc[1].push_back(m_postprocessingOps[op]);
          m_postprocessingMethods[1].push_back(&PostProcessing::computeParticleStatistics);
          initParticleStatistics();
          break;
        }
        case PP_ISO_TURBULENCE_STATISTICS: {
          m_postprocessingMethodsDesc[1].push_back(m_postprocessingOps[op]);
          m_postprocessingMethods[1].push_back(&PostProcessing::computeIsoTurbulenceStatistics);
          initIsoTurbulenceStatistics();
          break;
        }
        case PP_PL_ISO_TURBULENCE_STATISTICS: {
          m_postprocessingMethodsDesc[1].push_back(m_postprocessingOps[op]);
          m_postprocessingMethods[1].push_back(&PostProcessing::computePLIsoTurbulenceStatistics);
          initPLIsoTurbulenceStatistics();
          break;
        }
        default: {
          mTerm(1, AT_, "Unknown postprocessing operation");
        }
      }
    }
  } else {
    mTerm(1, AT_, "noProcessingOps for PostProcessing_" + to_string(m_postprocessingId) + " equal 0");
  }
}

/**
 * \fn void PostProcessing<nDim,ppType>::initTimeStepProperties()
 * \brief Initializes timestep properties for postprocessing
 *
 * \author Andreas Lintermann (last modified Ansgar Niemoeller, 07/14)
 * \date 14.09.2012
 *
 * reads properties pp_averageStartTimestep, pp_averageStopTimestep and pp_averageInterval
 *
 * \tparam[in] Solver solvertype
 *
 **/
template <MInt nDim, class ppType>
void PostProcessing<nDim, ppType>::initTimeStepProperties() {
  TRACE();

  // Init vars
  m_averageInterval = 1;
  m_averageStartTimestep = 0;
  m_averageStopTimestep = 0;
  m_averageRestart = false;

  /*! \page propertiesPP
    \section pp_averageInterval
    <code>MInt PostProcessing::m_averageInterval</code>\n
    default = <code>1</code>\n\n
    This property determines the time step interval used for averaging solutions.
    <ul>
    <li><code>interval</code> </li>
    </ul>
    Keywords: <i>POSTPROCESSING</i>
  */
  m_averageInterval =
      Context::getSolverProperty<MInt>("pp_averageInterval", m_postprocessingId, AT_, &m_averageInterval);

  if(m_averageInterval <= 0) {
    TERMM(1, "Please specify the property 'averageInterval' (with a value > 0) ...");
  }

  /*! \page propertiesPP
    \section pp_averageStartTimestep
    <code>MInt PostProcessing::m_averageStartTimestep</code>\n
    default = <code>0</code>\n\n
    This property determines the start timestep used for averaging.
    <ul>
    <li><code>timestep</code> </li>
    </ul>
    Keywords: <i>POSTPROCESSING</i>
  */
  m_averageStartTimestep =
      Context::getSolverProperty<MInt>("pp_averageStartTimestep", m_postprocessingId, AT_, &m_averageStartTimestep);

  if(m_averageStartTimestep < 0) {
    TERMM(1, "Please specify the property 'averageStartTimestep' (with a value >= 0) ...");
  }

  /*! \page propertiesPP
    \section pp_averageStopTimestep
    <code>MInt PostProcessing::m_averageStopTimestep</code>\n
    default = <code>0</code>\n\n
    This property determines the stop timestep used for averaging.
    <ul>
    <li><code>timestep</code> </li>
    </ul>
    Keywords: <i>POSTPROCESSING</i>
  */
  m_averageStopTimestep =
      Context::getSolverProperty<MInt>("pp_averageStopTimestep", m_postprocessingId, AT_, &m_averageStopTimestep);

  if(m_averageStopTimestep <= 0 || m_averageStopTimestep <= m_averageStartTimestep) {
    TERMM(1, "Please specify the property 'averageStopTimestep' (> averageStartTimestep) ...");
  }

  /*! \page propertiesPP
    \section pp_averageRestart
    <code>MInt PostProcessing::m_averageRestart</code>\n
    default = <code>0</code>\n\n
    This property determines if we should restart from our last averaging.
    <ul>
    <li><code>0</code> disabled</li>
    <li><code>1</code> enabled</li>
    </ul>
    Keywords: <i>GENERAL, POSTPROCESSING</i>
  */
  if(Context::propertyExists("pp_averageRestart", m_postprocessingId)) {
    m_averageRestart =
        Context::getSolverProperty<MBool>("pp_averageRestart", m_postprocessingId, AT_, &m_averageRestart);
  }
}

/// \fn void PostProcessing<nDim,ppType>::initAveragingProperties()
/// \brief Initialize properties relevant for temporal averaging
///
/// \author Michael Schlottke-Lakemper (mic) <mic@aia.rwth-aachen.de>
/// \date 2015-11-19
template <MInt nDim, class ppType>
void PostProcessing<nDim, ppType>::initAveragingProperties() {
  TRACE();

  m_correlation = Context::getSolverProperty<MBool>("pp_correlation", m_postprocessingId, AT_, &m_correlation);

  /*! \page propertiesPP
    \section pp_skewness
    <code>MInt PostprocesssingSolver::m_square</code>\n
    default = <code>0</code>\n\n
    This property determines if squares are computed for RMS values when PP_AVERAGE_PRE/IN/POST is activated.
    <ul>
    <li><code>0</code> deactivated</li>
    <li><code>1</code> active</li>
    </ul>\n
    Keywords: <i>POSTPROCESSING</i>
  */
  m_square = Context::getSolverProperty<MBool>("pp_square", m_postprocessingId, AT_, &m_square);

  /*! \page propertiesPP
    \section pp_skewness
    <code>MInt PostprocesssingSolver::m_skewness</code>\n
    default = <code>0</code>\n\n
    This property determines if skewness is computed when PP_AVERAGE_PRE/IN/POST is activated.
    <ul>
    <li><code>0</code> deactivated</li>
    <li><code>1</code> active</li>
    </ul>\n
    Keywords: <i>POSTPROCESSING</i>
  */
  m_skewness = Context::getSolverProperty<MBool>("pp_skewness", m_postprocessingId, AT_, &m_skewness);
  if(m_skewness) {
    m_square = true;
  }

  /*! \page propertiesPP
    \section pp_kurtosis
    <code>MInt PostprocesssingSolver::m_kurtosis</code>\n
    default = <code>0</code>\n\n
    This property determines if kurtosis (and skewness) is computed when PP_AVERAGE_PRE/IN/POST is activated.
    <ul>
    <li><code>0</code> deactivated</li>
    <li><code>1</code> active</li>
    </ul>\n
    Keywords: <i>POSTPROCESSING</i>
  */
  m_kurtosis = Context::getSolverProperty<MBool>("pp_kurtosis", m_postprocessingId, AT_, &m_kurtosis);
  if(m_kurtosis) {
    m_skewness = true;
    m_square = true;
  }

  /*! \page propertiesPP
    \section pp_twoPass
    <code>MInt PostprocesssingSolver::m_twoPass</code>\n
    default = <code>0</code>\n\n
    This property determines if two-pass averaging is performed in PP_AVERAGE_PRE/POST.\n
    Either m_twoPass or m_useKahan should be activated.
    <ul>
    <li><code>0</code> deactivated</li>
    <li><code>1</code> active</li>
    </ul>\n
    Keywords: <i>POSTPROCESSING</i>
  */
  m_twoPass = false;
  m_twoPass = Context::getSolverProperty<MBool>("pp_twoPass", m_postprocessingId, AT_, &m_twoPass);

  /*! \page propertiesPP
    \section pp_useKahan
    <code>MInt PostprocesssingSolver::m_useKahan</code>\n
    default = <code>0</code>\n\n
    This property determines if kahan summation is performed in PP_AVERAGE_PRE/IN/POST. \n
    Either m_twoPass or m_useKahan should be activated.
    <ul>
    <li><code>0</code> deactivated</li>
    <li><code>1</code> active</li>
    </ul>\n
    Keywords: <i>POSTPROCESSING</i>
  */
  m_useKahan = false;
  m_useKahan = Context::getSolverProperty<MBool>("pp_useKahan", m_postprocessingId, AT_, &m_useKahan);

  /*! \page propertiesPP
    \section pp_averageVorticity
    <code>MInt PostProcessing::m_averageVorticity</code>\n
    default = <code>0</code>\n
    This property determines if the vorticity vector is considered in the computation of averages.
    <ul>
    <li><code>0</code> disabled</li>
    <li><code>1</code> enabled</li>
    </ul>
    Keywords: <i>POSTPROCESSING</i>
  */

  m_averageVorticity = false;
  m_averageVorticity =
      Context::getSolverProperty<MBool>("pp_averageVorticity", m_postprocessingId, AT_, &m_averageVorticity);

  /*! \page propertiesPP
    \section pp_averageSpeedOfSound
    <code>MInt PostprocesssingSolver::m_averageSpeedOfSound</code>\n
    default = <code>0</code>\n\n
    This property determines if the speed of sound and its derivatives are considered in the
    computation of averages.
    <ul>
    <li><code>0</code> disabled</li>
    <li><code>1</code> enabled</li>
    </ul>
    Keywords: <i>POSTPROCESSING</i>
  */
  m_averageSpeedOfSound = false;
  m_averageSpeedOfSound =
      Context::getSolverProperty<MBool>("pp_averageSpeedOfSound", m_postprocessingId, AT_, &m_averageSpeedOfSound);

  // TODO labels:PP Rename this property as it is refering to specific thermal acoustic stuff
  m_acousticAnalysis = false;
  m_acousticAnalysis =
      Context::getSolverProperty<MBool>("acousticAnalysis", m_postprocessingId, AT_, &m_acousticAnalysis);

  /*! \page propertiesPP
      \section statisticCombustionAnalysis
      <code>MBool FvCartesianSolver::m_statisticCombustionAnalysis </code>\n
      default = <code>false</code>\n \n
      Sets the acoustic analysis /n
      possible values are:
      <ul>
      <li> true - on  </li>
      <li> false - off </li>
      </ul>
      Keywords: COMBUSTION, ACOUSTIC, ANALYSIS
  */
  m_statisticCombustionAnalysis = 0;
  m_statisticCombustionAnalysis = Context::getSolverProperty<MBool>("statisticCombustionAnalysis", m_postprocessingId,
                                                                    AT_, &m_statisticCombustionAnalysis);

  // Number of variable that need to be averaged for the source terms
  m_noSourceVars = m_activeMeanVars.size() * nDim; // NOTE: assumes each mean var has nDim values

  if(m_activeMeanVars.find(MV::VORT0) != m_activeMeanVars.end()) {
    // If VORT0 is part of active mean vars, the equation for the number of source variables changes
    // because of the vorticity. In 2D there is only one vorticity variable.
    // Number of active mean vars without voriticity
    m_noSourceVars -= nDim;
    // Add correct number of vorticites
    m_noSourceVars += 2 * nDim - 3;
  }

  // Correct number of source variables if GRADU is averaged (ndim*ndim)
  if(m_activeMeanVars.find(MV::GRADU) != m_activeMeanVars.end()) {
    if(m_activeMeanVars.find(MV::DU) != m_activeMeanVars.end()) {
      m_noSourceVars -= nDim; // averaging of DU is contained in GRADU
    }
    m_noSourceVars -= nDim;        // GRADU counted above with nDim
    m_noSourceVars += nDim * nDim; // add correct size of GRADU
  }

  // For some source terms additional properties are set
  for(auto&& sourceTerm : m_activeSourceTerms) {
    if(sourceTerm == ST::Q_mI_linear) {
      m_averageVorticity = true;
    }
    if(sourceTerm == ST::Q_c || sourceTerm == ST::Q_e) {
      m_averageSpeedOfSound = true;
    }
  }
}


/** \brief allocates memory for averageSolutions() and averageSolutionsInSolve()
 *
 * \author A. Niemoeller
 * \date 09.12.2013
 *
 **/
template <MInt nDim, class ppType>
void PostProcessing<nDim, ppType>::initAverageVariables() {
  // TODO:move to postData!

#if not defined(MAIA_DISABLE_LB)
  constexpr MBool isothermal = (std::is_same<ppType, PostProcessingLb<nDim>>::value);
  m_averageSpeedOfSound = (m_averageSpeedOfSound && !isothermal);
#else
  constexpr MBool isothermal = false;
#endif

  // Allocate memory for summed variables
  m_noAveragedVorticities = (m_averageVorticity) * (nDim * 2 - 3);
  m_noSpeedOfSoundVars = (m_averageSpeedOfSound) * (1 + nDim);
  m_needVorticity = (m_averageVorticity || (m_activeMeanVars.find(MV::LAMB0) != m_activeMeanVars.end()));

  const MInt noSummedVars =
      m_noVariables                                          // primitive variables
      + (3 * (nDim - 1) + 1) * m_square                      // Reynolds stresses +  pressure amplitude
      + nDim * m_skewness                                    // sknewness
      + nDim * m_kurtosis                                    // kurtosis
      + static_cast<MInt>(m_statisticCombustionAnalysis) * 3 // combustion variables (hm,c',h': cm in m_noVariables)
      + m_noAveragedVorticities                              // vorticities
      + m_noSpeedOfSoundVars                                 // speed of sound and derivatives
      + m_noSourceVars                                       // APE source terms
      + m_correlation * m_noCorrelationLines;

  if(noSummedVars > postData().noVariables()) {
    std::stringstream ss;
    ss << "PostData: noVariables in postprocessing not correct: postData().noVariables(): " << postData().noVariables()
       << " m_varNames.size(): " << noSummedVars << std::endl;
    TERMM(1, ss.str());
  }

  MInt count = 0;
  postData().setPropertyVariableOffset("primitive", count, m_noVariables);
  count += m_noVariables;

  if(m_square) {
    const MInt offset = 3 * (nDim - 1) + 1;
    postData().setPropertyVariableOffset("square", count, offset);
    count += offset;
  }

  if(m_kurtosis) {
    postData().setPropertyVariableOffset("skewness", count, nDim);
    count += nDim;
    postData().setPropertyVariableOffset("kurtosis", count, nDim);
    count += nDim;
    ASSERT(m_skewness, "required as well");
  } else if(m_skewness) {
    postData().setPropertyVariableOffset("skewness", count, nDim);
    count += nDim;
  }

  if(m_statisticCombustionAnalysis) {
    postData().setPropertyVariableOffset("statisticCombustionAnalysis", count, 3);
    count += 3;
  }

  if(m_averageVorticity) {
    postData().setPropertyVariableOffset("averageVorticity", count, m_noAveragedVorticities);
    count += m_noAveragedVorticities;
  }

  if(m_averageSpeedOfSound) {
    postData().setPropertyVariableOffset("averageSpeedOfSound", count, m_noSpeedOfSoundVars);
    count += m_noSpeedOfSoundVars;
  }

  if(m_noSourceVars > 0) {
    postData().m_sourceVarsIndex = count;
    if(m_activeMeanVars.find(MV::LAMB0) != m_activeMeanVars.end()) {
      postData().setPropertyVariableOffset("lamb0", count, nDim);
      count += nDim;
    }

    if(m_activeMeanVars.find(MV::GRADU) != m_activeMeanVars.end()) {
      postData().setPropertyVariableOffset("gradu", count, nDim * nDim);
      count += nDim * nDim;
    } else if(m_activeMeanVars.find(MV::DU) != m_activeMeanVars.end()) {
      postData().setPropertyVariableOffset("du", count, nDim);
      count += nDim;
    }

    if(m_activeMeanVars.find(MV::UGRADU) != m_activeMeanVars.end()) {
      postData().setPropertyVariableOffset("ugradu", count, nDim);
      count += nDim;
    }

    if(m_activeMeanVars.find(MV::DRHO) != m_activeMeanVars.end()) {
      postData().setPropertyVariableOffset("drho", count, nDim);
      count += nDim;
    }

    if(m_activeMeanVars.find(MV::DP) != m_activeMeanVars.end()) {
      postData().setPropertyVariableOffset("dp", count, nDim);
      count += nDim;
    }

    if(m_activeMeanVars.find(MV::RHODIVU) != m_activeMeanVars.end()) {
      postData().setPropertyVariableOffset("rhodivu", count, nDim);
      count += nDim;
    }

    if(m_activeMeanVars.find(MV::UGRADRHO) != m_activeMeanVars.end()) {
      postData().setPropertyVariableOffset("ugradrho", count, nDim);
      count += nDim;
    }

    if(m_activeMeanVars.find(MV::GRADPRHO) != m_activeMeanVars.end()) {
      postData().setPropertyVariableOffset("gradprho", count, nDim);
      count += nDim;
    }
  }

  if(m_correlation) {
    postData().setPropertyVariableOffset("correlation", count, m_noCorrelationLines);
    count += m_noCorrelationLines;
  }

  postData().setVariableNames();

  if(count > postData().noVariables()) {
    cerr << "Missmatching post data " << count << " " << postData().noVariables() << endl;
    mTerm(1, AT_, "Variables size not matching!");
  }

  // Allocate memory for squared variables
  // mAlloc(m_square, noInternalCells, 3 * (nDim - 1) + 1 + (MInt)(m_statisticCombustionAnalysis)*2, "m_square", F0,
  //           FUN_);
  // + 1 for pressure ampl + 2 ampl for combustion variables, c',h'

  // if(m_kurtosis != 0 /*&& !m_twoPass*/) {
  //   mAlloc(m_cube, noInternalCells, nDim, "m_cube", F0, FUN_);
  //   mAlloc(m_fourth, noInternalCells, nDim, "m_fourth", F0, FUN_);
  // } else if(m_skewness != 0 /*&& !m_twoPass*/) {
  //   mAlloc(m_cube, noInternalCells, nDim, "m_cube", F0, FUN_);
  // }

  if(m_useKahan) // allocate memory for kahan summation
  {
    m_log << "m_useKahan is activated" << endl;
    TERMM(1, "FIXME Kahan summation is untested and probably broken");
    // mAlloc(m_cSum, noInternalCells, m_noVariables + 3 * (nDim - 1), "m_cSum", F0, AT_);
    // mAlloc(m_tSum, noInternalCells, m_noVariables + 3 * (nDim - 1), "m_tSum", F0, AT_);
    // mAlloc(m_ySum, noInternalCells, m_noVariables + 3 * (nDim - 1), "m_ySum", F0, AT_);
    // mAlloc(m_cSquare, noInternalCells, 3 * (nDim - 1), "m_cSquare", F0, AT_);
    // mAlloc(m_tSquare, noInternalCells, 3 * (nDim - 1), "m_tSquare", F0, AT_);
    // mAlloc(m_ySquare, noInternalCells, 3 * (nDim - 1), "m_ySquare", F0, AT_);
    if(m_kurtosis) {
      // mAlloc(m_cCube, noInternalCells, nDim, "m_cCube", F0, AT_);
      // mAlloc(m_tCube, noInternalCells, nDim, "m_tCube", F0, AT_);
      // mAlloc(m_yCube, noInternalCells, nDim, "m_yCube", F0, AT_);
      // mAlloc(m_cFourth, noInternalCells, nDim, "m_cFourth", F0, AT_);
      // mAlloc(m_tFourth, noInternalCells, nDim, "m_tFourth", F0, AT_);
      // mAlloc(m_yFourth, noInternalCells, nDim, "m_yFourth", F0, AT_);
    } else if(m_skewness) {
      // mAlloc(m_cCube, noInternalCells, nDim, "m_cCube", F0, AT_);
      // mAlloc(m_tCube, noInternalCells, nDim, "m_tCube", F0, AT_);
      // mAlloc(m_yCube, noInternalCells, nDim, "m_yCube", F0, AT_);
    }
  }
}


template <MInt nDim, class ppType>
void PostProcessing<nDim, ppType>::initCorrelation() {
  TRACE();

  if(!postData().grid().isActive()) return;

  GridProxy* gridProxy = &postData().grid();

  // Read properties
  m_noCorrelationLines = Context::propertyLength("pp_correlationDirection", m_postprocessingId);
  if(m_noCorrelationLines <= 0) {
    TERMM(1, "Error: line probing activated but number of probe line directions is "
                 + std::to_string(m_noCorrelationLines) + "!");
  }

  mAlloc(m_correlationDirection, m_noCorrelationLines, "m_correlationDirection", AT_);
  mAlloc(m_correlationVariableIndex, m_noCorrelationLines, "m_correlationVariableIndex", AT_);

  for(MInt correlationId = 0; correlationId < m_noCorrelationLines; correlationId++) {
    m_correlationDirection[correlationId] =
        Context::getSolverProperty<MInt>("pp_correlationDirection", m_postprocessingId, AT_, correlationId);
    if(m_correlationDirection[correlationId] < 0 || m_correlationDirection[correlationId] >= nDim) {
      m_log << "invalid value for property correlationDirection: " << m_correlationDirection << " at position "
            << correlationId << endl;
    }
    m_correlationVariableIndex[correlationId] =
        Context::getSolverProperty<MInt>("pp_correlationVariableIndex", m_postprocessingId, AT_, correlationId);
  }

  ScratchSpace<MInt> correlationDir(m_noCorrelationLines, (nDim - 1), AT_, "correlationDir");
  mAlloc(m_correlationCoordinates, m_noCorrelationLines, (nDim - 1), "m_correlationCoordinates", AT_);

  for(MInt correlationId = 0; correlationId < m_noCorrelationLines; correlationId++) {
    if(m_correlationDirection[correlationId] == 0) {
      correlationDir(correlationId, 0) = 1;
    } else {
      correlationDir(correlationId, 0) = 0;
    }
    m_correlationCoordinates[correlationId][0] = Context::getSolverProperty<MFloat>(
        "pp_correlationCoordinates", m_postprocessingId, AT_, (nDim - 1) * correlationId);
    IF_CONSTEXPR(nDim == 3) {
      if(m_correlationDirection[correlationId] == 2) {
        correlationDir(correlationId, 1) = 1;
      } else {
        correlationDir(correlationId, 1) = 2;
      }
      m_correlationCoordinates[correlationId][1] = Context::getSolverProperty<MFloat>(
          "pp_correlationCoordinates", m_postprocessingId, AT_, 2 * correlationId + 1);
    }
  }

  for(MInt correlationId = 0; correlationId < m_noCorrelationLines; correlationId++) {
    m_log << "correlation in correlation function direction " << m_correlationDirection[correlationId]
          << " with correlationCorrelationDirection "
          << "m_noCorrelationLines is " << m_noCorrelationLines;
    m_log << std::endl;
  }
  mAlloc(m_noCorrelationIds, m_noCorrelationLines, "m_noCorrelationIds", AT_);
  mAlloc(m_correlationIds, m_noCorrelationLines, "m_correlationIds", AT_);
  mAlloc(m_correlationExchangeVar, m_noCorrelationLines, "m_correlationExchangeVar", AT_);
  mAlloc(m_correlationExchangeVarMean, m_noCorrelationLines, "m_correlationExchangeVarMean", AT_);
  mAlloc(m_correlationExchangeVarRMS, m_noCorrelationLines, "m_correlationExchangeVarRMS", AT_);
  mAlloc(m_correlationPositions, m_noCorrelationLines, "m_correlationPositions", AT_);
  mAlloc(m_globalNoCorrelationIds, m_noCorrelationLines, "m_globalNoCorrelationIds", AT_);

  // suppress valgrind error
  for(MInt i = 0; i < m_noCorrelationLines; i++) {
    m_correlationIds[i] = nullptr;
    m_correlationPositions[i] = nullptr;
  }

  for(MInt correlationId = 0; correlationId < m_noCorrelationLines; correlationId++) {
    map<MFloat, MInt, coord_comp_1d_> coordinates;
    std::vector<MFloat> coord(nDim, std::numeric_limits<MFloat>::max());
    MFloat halfCellLength = NAN;
    MBool ok;
    MBool ok2;
    MBool nghbrFound;
    for(MInt cellId = 0; cellId < postData().grid().noInternalCells(); cellId++) {
      if(gridProxy->tree().noChildren(cellId) > 0) continue; // skip non leaf cells
      for(MInt dim = 0; dim < nDim; dim++) {
        coord[dim] = postData().c_coordinate(cellId, dim);
      }
      ok = true;
      ok2 = true;
      nghbrFound = false;

      halfCellLength = gridProxy->halfCellLength(cellId);
      for(MInt i = 0; i < nDim - 1; i++) {
        // check if distance is greater than half of the cell length (in direction i)
        if(abs(m_correlationCoordinates[correlationId][i] - coord[correlationDir(correlationId, i)])
           >= halfCellLength) {
          ok = false;
        }

        // check if distance is a bit more than half of the cell length (line between cells)
        if(abs(m_correlationCoordinates[correlationId][i] - coord[correlationDir(correlationId, i)])
           < (halfCellLength + MFloatEps)) {
          for(MInt dirId = 0; dirId < 2 * nDim; dirId++) {
            if(dirId == 2 * m_correlationDirection[correlationId]
               || dirId == 2 * m_correlationDirection[correlationId] + 1) {
              continue;
            }
            if(gridProxy->tree().hasNeighbor(cellId, dirId) != 1) {
              continue;
            }
            MInt nId = gridProxy->tree().neighbor(cellId, dirId);
            // check if neighbor cell is already in the probe line
            if(coordinates.find(postData().c_coordinate(nId, m_correlationDirection[correlationId]))
               != coordinates.end()) {
              nghbrFound = true;
            }
          }
        } else {
          ok2 = false;
        }
      }
      if(ok || (ok2 && !nghbrFound)) {
        if(postData().a_isBndryGhostCell(cellId)) continue;
        coordinates[coord[m_correlationDirection[correlationId]]] = cellId;
      }
    }

    m_noCorrelationIds[correlationId] = coordinates.size(); // local number of cells

    if(m_noCorrelationIds[correlationId] > 0) {
      mAlloc(m_correlationIds[correlationId], m_noCorrelationIds[correlationId], "m_correlationIds", AT_);
      mAlloc(m_correlationPositions[correlationId], m_noCorrelationIds[correlationId], "m_correlationPositions", AT_);
    } else {
      // Allocate dummy array to avoid errors when dereferencing during writing
      mAlloc(m_correlationIds[correlationId], 1, "m_correlationIds", AT_);
      mAlloc(m_correlationPositions[correlationId], 1, "m_correlationPositions", AT_);
    }
    auto it = coordinates.begin();
    for(MInt i = 0; it != coordinates.end(); it++, i++) {
      m_correlationIds[correlationId][i] = (*it).second;
      m_correlationPositions[correlationId][i] = (*it).first;
    }
  }

  // fill global arrays
  MPI_Allreduce(m_noCorrelationIds, m_globalNoCorrelationIds, m_noCorrelationLines, MPI_INT, MPI_SUM,
                postData().mpiComm(), AT_, "m_noCorrelationIds", "m_globalNoCorrelationIds");

  MInt maxGlobalNoCorrelationIds = 0;
  for(MInt correlationId = 0; correlationId < m_noCorrelationLines; correlationId++) {
    if(m_globalNoCorrelationIds[correlationId] >= maxGlobalNoCorrelationIds) {
      maxGlobalNoCorrelationIds = m_globalNoCorrelationIds[correlationId];
    }
  }

  mAlloc(m_globalCorrelationIds, m_noCorrelationLines, maxGlobalNoCorrelationIds, "m_globalCorrelationIds", 0, AT_);
  mAlloc(m_globalCorrelationPositions, m_noCorrelationLines, maxGlobalNoCorrelationIds, "m_globalCorrelationPositions",
         F0, AT_);
  mAlloc(m_globalCorrelationExchangeVar, m_noCorrelationLines, maxGlobalNoCorrelationIds,
         "m_globalCorrelationExchangeVar", F0, AT_);
  mAlloc(m_globalCorrelationExchangeVarMean, m_noCorrelationLines, maxGlobalNoCorrelationIds,
         "m_globalCorrelationExchangeVarMean", F0, AT_);
  mAlloc(m_globalCorrelationExchangeVarRMS, m_noCorrelationLines, maxGlobalNoCorrelationIds,
         "m_globalCorrelationExchangeVarRMS", F0, AT_);

  // create mapping for cells
  mAlloc(m_correlationIndexMapping, m_noCorrelationLines, postData().noInternalCells(), "m_correlationIndexMapping", -1,
         AT_);


  for(MInt correlationId = 0; correlationId < m_noCorrelationLines; correlationId++) {
    MInt noDomain = postData().noDomains();
    ScratchSpace<MInt> recvbuf(noDomain, "recvbuf", AT_);
    recvbuf.fill(0);

    MPI_Gather(&m_noCorrelationIds[correlationId], 1, MPI_INT, &recvbuf[0], 1, MPI_INT, 0, postData().mpiComm(), AT_,
               "m_noCorrelationIds[correlationId]", "postData().noDomain()");

    ScratchSpace<MInt> displs(noDomain, "displs", AT_);
    if(postData().domainId() == 0) {
      MInt offset = 0;
      for(MInt dom = 0; dom < noDomain; dom++) {
        displs[dom] = offset;
        offset += recvbuf[dom];
      }
    }

    MPI_Gatherv(&m_correlationIds[correlationId][0], m_noCorrelationIds[correlationId], MPI_INT,
                &m_globalCorrelationIds[correlationId][0], &recvbuf[postData().domainId()],
                &displs[postData().domainId()], MPI_INT, 0, postData().mpiComm(), AT_, "m_correlationIds",
                "m_globalCorrelationIds");

    MPI_Gatherv(&m_correlationPositions[correlationId][0], m_noCorrelationIds[correlationId], MPI_DOUBLE,
                &m_globalCorrelationPositions[correlationId][0], &recvbuf[postData().domainId()],
                &displs[postData().domainId()], MPI_DOUBLE, 0, postData().mpiComm(), AT_, "m_correlationPositions",
                "m_globalCorrelationPositions");

    MPI_Bcast(m_globalCorrelationIds[correlationId], m_globalNoCorrelationIds[correlationId], MPI_INT, 0,
              postData().mpiComm(), AT_, "m_probeLineAverageCoordinates");

    MPI_Bcast(m_globalCorrelationPositions[correlationId], m_globalNoCorrelationIds[correlationId], MPI_DOUBLE, 0,
              postData().mpiComm(), AT_, "m_globalCorrelationPositions");

    mAlloc(m_correlationExchangeVar[correlationId], postData().grid().noInternalCells(),
           "m_correlationExchangeVar[correlationId]", AT_);
    mAlloc(m_correlationExchangeVarMean[correlationId], postData().grid().noInternalCells(),
           "m_correlationExchangeVarMean[correlationId]", AT_);
    mAlloc(m_correlationExchangeVarRMS[correlationId], postData().grid().noInternalCells(),
           "m_correlationExchangeVarRMS[correlationId]", AT_);

    for(MInt dataId = 0; dataId < postData().noInternalCells(); dataId++) {
      MFloat dataPosition = postData().c_coordinate(dataId, m_correlationDirection[correlationId]);
      MFloat smallerDist = std::numeric_limits<MFloat>::max();
      MFloat biggerDist = std::numeric_limits<MFloat>::max();
      MInt index1 = -1;
      MInt index2 = -1;
      for(MInt i = 0; i < m_globalNoCorrelationIds[correlationId]; i++) {
        MFloat corrPosition = m_globalCorrelationPositions[correlationId][i];
        // change this
        if(m_correlationIndexMapping[correlationId][dataId] == -1) {
          MFloat dist = abs(corrPosition - dataPosition);
          if(corrPosition <= dataPosition && dist < smallerDist) {
            smallerDist = dist;
            index1 = i;
          }
          if(corrPosition > dataPosition && dist < biggerDist) {
            biggerDist = dist;
            index2 = i;
          }
        }
      }

      if(smallerDist <= biggerDist) {
        m_correlationIndexMapping[correlationId][dataId] = index1;
      } else {
        m_correlationIndexMapping[correlationId][dataId] = index2;
      }
    }
  }
}


/**
 * \fn void PostProcessing<nDim,ppType>::initMovingAverage()
 * \brief Initializes properties and allocates memory for moving averaging
 *
 * \author Ansgar Niemoeller
 * \date 09.07.2014
 *
 * reads properties pp_movingAvgInterval, pp_movingAvgDataPoints
 * (pp_averageVorticity already read in initAveragingVariables)
 *
 * \param[in] grid pointer to the grid
 *
 **/
template <MInt nDim, class ppType>
void PostProcessing<nDim, ppType>::initMovingAverage() {
  TRACE();

  /*! \page propertiesPP
    \section pp_movingAverageInterval
    <code>MInt PostprocesssingSolver::m_movingAverageInterval</code>\n
    default = <code>1</code>\n\n
    This property determines the interval between timesteps considered for moving average,
    e.g. if set to 2 at an averaging timestep n (see pp_averagStartTimestep, pp_averageStopTimestep,
    pp_averageInterval) the timesteps n, n-2, n-4, ... are used to compute the moving average.\n
    Note: this has to be a factor of m_averageInterval\n
    See also pp_movingAverageDataPoints, pp_averageVorticity.\n
    Keywords: <i>POSTPROCESSING</i>
  */
  m_movingAverageInterval =
      Context::getSolverProperty<MInt>("pp_movingAverageInterval", m_postprocessingId, AT_, &m_movingAverageInterval);

  /*! \page propertiesPP
  \section pp_movingAverageStartTimestep
  <code>MInt PostProcessing::m_movingAverageStartTimestep</code>\n
  default = <code>0</code>\n\n
  This property determines the start timestep used for averaging.\n
  Keywords: <i>POSTPROCESSING, MOVINGAVERAGE</i>
*/
  m_movingAverageStartTimestep = Context::getSolverProperty<MInt>("pp_movingAverageStartTimestep", m_postprocessingId,
                                                                  AT_, &m_averageStartTimestep);

  if(m_movingAverageStartTimestep <= 0) {
    TERMM(1, "Please specify the property 'averageStartTimestep' (with a value > 0) ...");
  }

  /*! \page propertiesPP
    \section pp_movingAverageStopTimestep
    <code>MInt PostProcessing::m_movingAverageStopTimestep</code>\n
    default = <code>0</code>\n\n
    This property determines the stop timestep used for averaging.\n
    Keywords: <i>POSTPROCESSING, MOVINGAVERAGE</i>
  */
  m_movingAverageStopTimestep =
      Context::getSolverProperty<MInt>("pp_movingAverageStopTimestep", m_postprocessingId, AT_, &m_averageStopTimestep);

  if(m_movingAverageInterval < 1) {
    TERMM(1, "m_movingAverageInterval has to be >=1");
  }

  /*! \page propertiesPP
    \section pp_movingAverageDataPoints
    <code>MInt PostprocesssingSolver::m_movingAverageDataPoints</code>\n
    default = <code>empty</code>\n\n
    This property determines the number of timesteps (data points) used for moving average
    computation.\n See also pp_movingAverageInterval, pp_averageVorticity.\n
    <ul>
    <li> integer <code>&ge 2</code> </li>
    </ul>\n
    Keywords: <i>POSTPROCESSING</i>
  */
  m_movingAverageDataPoints = Context::getSolverProperty<MInt>("pp_movingAverageDataPoints", m_postprocessingId, AT_,
                                                               &m_movingAverageDataPoints);
  if(m_movingAverageDataPoints < 2) {
    TERMM(1, "m_movingAverageDataPoints has to be at least 2");
  }

  m_movingAverageCounter = 0;

  m_movAvgNoVariables = m_noVariables;
  if(m_averageVorticity) {
    m_movAvgNoVariables += (2 * nDim - 3);
  }
  mAlloc(m_movAvgVariables, pp->solver().grid().noInternalCells(), m_movAvgNoVariables * m_movingAverageDataPoints,
         "m_movAvgVariables", F0, AT_);
  mAlloc(m_movAvgVarNames, 2 * m_movAvgNoVariables, "m_movAvgVarNames", MString("default"), AT_);

  pp->getPrimitiveVarNames(m_movAvgVarNames);
  if(m_averageVorticity == 1) {
    if(nDim == 2) {
      m_movAvgVarNames[m_noVariables] = "w";
    } else if(nDim == 3) {
      m_movAvgVarNames[m_noVariables] = "wx";
      m_movAvgVarNames[m_noVariables + 1] = "wy";
      m_movAvgVarNames[m_noVariables + 2] = "wz";
    }
  }

  for(MInt varId = 0; varId < m_movAvgNoVariables; varId++) {
    m_movAvgVarNames[m_movAvgNoVariables + varId] = m_movAvgVarNames[varId] + "m";
  }
}

/**
 * \fn void PostProcessing<nDim,ppType>::initProbePoint()
 * \brief Initializes properties for point probing
 *
 * \author Andreas Lintermann
 * \date 14.09.2012
 *
 * \tparam[in] T celltype
 * \param[in]  grid pointer to the grid
 **/
template <MInt nDim, class ppType>
void PostProcessing<nDim, ppType>::initProbePoint() {
  TRACE();

  m_noProbePoints = Context::propertyLength("pp_probeCoordinates", m_postprocessingId) / nDim;

  m_log << "        = number of probe points: " << m_noProbePoints << endl;

  mAlloc(m_probeCoordinates, m_noProbePoints, nDim, "m_probeCoordinates", F0, AT_);
  mAlloc(m_probeCellIds, m_noProbePoints, "m_probeCellIds", AT_);

  m_probeFileStreams = new ofstream[m_noProbePoints];

  /*! \page propertiesPP
    \section pp_probeCoordinates
    <code>MFloat** Solver::m_probeCoordinates</code>\n
    default = <code>(0,0,0), ... ,(0,0,0)</code>\n\n
    This property defines the probe coordinates.
    <ul>
    <li><code>coordinates</code> </li>
    </ul>
    Keywords: <i>POSTPROCESSING</i>
  */

  for(MInt np = 0, t = 0; np < m_noProbePoints; np++) {
    m_probeCellIds[np] = -1;
    for(MInt i = 0; i < nDim; i++, t++)
      m_probeCoordinates[np][i] = Context::getSolverProperty<MFloat>("pp_probeCoordinates", m_postprocessingId, AT_, t);
  }

  /*! \page propertiesPP
    \section pp_probePath
    <code>MString PostProcessing::m_probePath</code>\n
    default = <code>"./probes/"</code>\n\n
    This property defines the path for storing probe files.
    <ul>
    <li><code>relative path</code> </li>
    </ul>
    Keywords: <i>GENERAL, POSTPROCESSING</i>
  */

  m_probePath = "./probes/";
  m_probePath = Context::getSolverProperty<MString>("pp_probePath", m_postprocessingId, AT_, &m_probePath);

  m_log << "        = probe path: " << m_probePath << endl;

  findClosestProbePointsGridCells();

  m_log << "        = local probe points: " << endl;
  for(MInt np = 0; np < m_noProbePoints; np++) {
    if(m_probeCellIds[np] != -1) {
      MString filename = m_probePath + "probe_";
      MChar buf[100];
      sprintf(buf, "%d", np);

      filename.append(buf);
      filename += ".dat";

      m_probeFileStreams[np].open(filename.c_str(), ios_base::app);
      m_probeFileStreams[np].precision(16);
      m_probeFileStreams[np] << "# Coordinates provided (really used) - local cellId:\n#";
      for(MInt i = 0; i < nDim; i++)
        m_probeFileStreams[np] << "   " << m_probeCoordinates[np][i] << " ("
                               << pp->solver().a_coordinate(m_probeCellIds[np], i) << ")  ";
      m_probeFileStreams[np] << " - " << m_probeCellIds[np] << endl;

      m_log << "          # number " << np << endl;
      m_log << "            # coordinates: ";
      for(MInt d = 0; d < nDim; d++)
        m_log << m_probeCoordinates[np][d] << " ";
      m_log << endl;
      m_log << "            # according cell id: " << m_probeCellIds[np] << endl;
      m_log << "            # real coordinates: ";
      for(MInt d = 0; d < nDim; d++)
        m_log << pp->solver().a_coordinate(m_probeCellIds[np], d) << " ";
      m_log << endl;
    }
  }
}


/**
 * \fn void PostProcessing<nDim,ppType>::initProbeLine()
 * \brief Initializes properties and memory for line probing
 *
 * \author A. Niemoeller
 * \date 09.12.2013
 *
 * reads properties probeLineDirection and probeLineCoordinates which define a line (or multiple lines) in the specified
 *direction with fixed coordinate(s) in the other dimensions \n assembles all cell Ids for probing and communicates data
 *for later line probing \n
 *
 * \param[in] grid pointer to the grid
 *
 **/
template <MInt nDim, class ppType>
void PostProcessing<nDim, ppType>::initProbeLine() {
  TRACE();

  if(!pp->solver().grid().isActive()) return;

  /*! \page propertiesPP
    \section pp_probeLineDirection
    <code>MInt* PostprocesssingSolver::m_probeLineDirection</code>\n
    default = <code>empty</code>\n
    This property determines the coordinate direction for each probe line (size: number of probe
    lines).
    <ul>
    <li><code>0</code> x</li>
    <li><code>1</code> y</li>
    <li><code>2</code> z (3D)</li>
    </ul>
    Keywords: <i>POSTPROCESSING</i>
  */
  m_noProbeLines = Context::propertyLength("pp_probeLineDirection", m_postprocessingId);
  if(m_noProbeLines <= 0) {
    TERMM(1, "Error: line probing activated but number of probe line directions is " + std::to_string(m_noProbeLines)
                 + "!");
  }

  mAlloc(m_probeLineDirection, m_noProbeLines, "m_probeLineDirection", AT_);

  for(MInt probeLineId = 0; probeLineId < m_noProbeLines; probeLineId++) { // read all probe line directions
    m_probeLineDirection[probeLineId] =
        Context::getSolverProperty<MInt>("pp_probeLineDirection", m_postprocessingId, AT_, probeLineId);
    if(m_probeLineDirection[probeLineId] < 0 || m_probeLineDirection[probeLineId] >= nDim) {
      m_log << "invalid value for property probeLineDirection: " << m_probeLineDirection << " at position "
            << probeLineId << endl;
    }
  }

  /*! \page propertiesPP
    \section pp_probeLineCoordinates
    <code>MInt** PostprocesssingSolver::m_probeLineCoordinates</code>\n
    default = <code>empty</code>\n
    This property determines the coordinate(s) of the probe line(s).\n
    In 3D the first provided value corresponds to the lower of the remaining directions (e.g. x if
    probeLineDirection=1(=y)).
    <ul>
    <li>coordinate array of size <code>(nDim-1)*(number of probe lines)</code></li>
    </ul>
    Keywords: <i>POSTPROCESSING</i>
  */

  ScratchSpace<MInt> lineDir(m_noProbeLines, (nDim - 1), AT_, "lineDir");
  mAlloc(m_probeLineCoordinates, m_noProbeLines, (nDim - 1), "m_probeLineCoordinates", AT_);

  for(MInt probeLineId = 0; probeLineId < m_noProbeLines;
      probeLineId++) { // read line coordinates and set corresponding directions
    if(m_probeLineDirection[probeLineId] == 0) {
      lineDir(probeLineId, 0) = 1;
    } else {
      lineDir(probeLineId, 0) = 0;
    }
    m_probeLineCoordinates[probeLineId][0] = Context::getSolverProperty<MFloat>(
        "pp_probeLineCoordinates", m_postprocessingId, AT_, (nDim - 1) * probeLineId);
    IF_CONSTEXPR(nDim == 3) {
      if(m_probeLineDirection[probeLineId] == 2)
        lineDir(probeLineId, 1) = 1;
      else
        lineDir(probeLineId, 1) = 2;
      m_probeLineCoordinates[probeLineId][1] =
          Context::getSolverProperty<MFloat>("pp_probeLineCoordinates", m_postprocessingId, AT_, 2 * probeLineId + 1);
    }
  }

  for(MInt probeLineId = 0; probeLineId < m_noProbeLines; probeLineId++) { // write all lines to log
    m_log << "probeLine in direction " << m_probeLineDirection[probeLineId] << " with probeLineCoordinates ";
    for(MInt i = 0; i < nDim - 1; i++) {
      m_log << lineDir(probeLineId, i) << " " << m_probeLineCoordinates[probeLineId][i] << " ";
    }
    m_log << endl;
  }

  mAlloc(m_noProbeLineIds, m_noProbeLines, "m_noProbeLineIds", AT_); // local number of cells for every line
  mAlloc(m_probeLineIds, m_noProbeLines, "m_probeLineIds", AT_);     // local ids of cells for every line
  mAlloc(m_probeLinePositions, m_noProbeLines, "m_probeLinePositions", AT_);
  mAlloc(m_probeLineOffsets, m_noProbeLines, "m_probeLineOffsets", AT_);
  mAlloc(m_globalNoProbeLineIds, m_noProbeLines, "m_globalNoProbeLineIds", AT_); // global #cells on every line
  // mAlloc(m_noGlobalProbeLineIds, m_noProbeLines, pp->solver().m_noDomains, "m_noGlobalProbeLineIds", AT_);

  // suppress valgrind error
  for(MInt i = 0; i < m_noProbeLines; i++) {
    m_probeLineIds[i] = nullptr;
    m_probeLinePositions[i] = nullptr;
  }

  for(MInt probeLineId = 0; probeLineId < m_noProbeLines; probeLineId++) { // loop over all probe lines

    map<MFloat, MInt, coord_comp_1d_> coordinates;
    const MFloat* coord = nullptr;
    MFloat halfCellLength = NAN;
    MBool ok;
    MBool ok2;
    MBool nghbrFound;
    for(MInt cellId = 0; cellId < pp->solver().grid().noInternalCells(); cellId++) {
      if(m_gridProxy->tree().noChildren(cellId) > 0) continue; // skip non leaf cells
      coord = &pp->solver().a_coordinate(cellId, 0);
      ok = true;
      ok2 = true;
      nghbrFound = false;

      halfCellLength = m_gridProxy->halfCellLength(cellId);
      for(MInt i = 0; i < nDim - 1; i++) {
        // check if distance is greater than half of the cell length (in direction i)
        if(abs(m_probeLineCoordinates[probeLineId][i] - coord[lineDir(probeLineId, i)]) >= halfCellLength) ok = false;

        // check if distance is a bit more than half of the cell length (line between cells)
        if(abs(m_probeLineCoordinates[probeLineId][i] - coord[lineDir(probeLineId, i)])
           < (halfCellLength + MFloatEps)) {
          for(MInt dirId = 0; dirId < 2 * nDim; dirId++) { // check all directions (x-,x+,...)
            if(dirId == 2 * m_probeLineDirection[probeLineId]
               || dirId == 2 * m_probeLineDirection[probeLineId] + 1) { // skip line directions
              continue;
            }
            if(m_gridProxy->tree().hasNeighbor(cellId, dirId) != 1) { // skip if cells has no equal level neighbor in
                                                                      // direction dirId
              continue;
            }
            MInt nId = m_gridProxy->tree().neighbor(cellId, dirId);
            // check if neighbor cell is already in the probe line
            if(coordinates.find(pp->solver().a_coordinate(nId, m_probeLineDirection[probeLineId]))
               != coordinates.end()) {
              nghbrFound = true;
            }
          }
        } else {
          ok2 = false;
        }
      }
      if(ok || (ok2 && !nghbrFound)) {
        if(pp->solver().a_isBndryGhostCell(cellId)) continue;
        coordinates[coord[m_probeLineDirection[probeLineId]]] = cellId;
        // m_log << "adding cellId to probeLine " << cellId << " " << coord[0] << " " << coord[1] << " " << coord[2]
        // << "\t" << m_probeLineCoordinates[probeLineId][0] << " " << m_probeLineCoordinates[probeLineId][1] << " " <<
        // " " << halfCellLength << endl;
      }
    }

    m_noProbeLineIds[probeLineId] = coordinates.size(); // local number of cells

    if(m_noProbeLineIds[probeLineId] > 0) {
      mAlloc(m_probeLineIds[probeLineId], m_noProbeLineIds[probeLineId], "m_probeLineIds",
             AT_); // collector for cellIds
      mAlloc(m_probeLinePositions[probeLineId], m_noProbeLineIds[probeLineId], "m_probeLinePositions", AT_);
    } else {
      // Allocate dummy array to avoid errors when dereferencing during writing
      mAlloc(m_probeLinePositions[probeLineId], 1, "m_probeLinePositions", AT_);
    }

    // testing
    // MPI_Gather(&(m_noProbeLineIds[probeLineId]), 1, MPI_INT,
    // 	       m_noGlobalProbeLineIds[probeLineId], 1, MPI_INT, 0,
    //            pp->solver().mpiComm(), AT_,
    // 	       "(m_noProbeLineIds[probeLineId])", "m_noGlobalProbeLineIds[probeLineId]");
    // m_log << "line #" << probeLineId << ": ";
    // for(MInt i=0; i<pp->solver().grid().m_noDomains; i++){
    //  m_log << m_noGlobalProbeLineIds[probeLineId][i] << " ";
    //}
    // m_log << endl;
    //

    auto it = coordinates.begin();
    for(MInt i = 0; it != coordinates.end();
        it++, i++) { // sorted local cell ids for probing and their respective coordinates
      m_probeLineIds[probeLineId][i] = (*it).second;
      m_probeLinePositions[probeLineId][i] = (*it).first;
    }

    ParallelIo::size_type localCount, offset, totalCount;
    localCount = m_noProbeLineIds[probeLineId];
    ParallelIo::calcOffset(localCount, &offset, &totalCount, pp->solver().mpiComm());
    m_probeLineOffsets[probeLineId] = offset;
  }


  MPI_Allreduce(m_noProbeLineIds, m_globalNoProbeLineIds, m_noProbeLines, MPI_INT, MPI_SUM, pp->solver().mpiComm(), AT_,
                "m_noProbeLineIds", "m_globalNoProbeLineIds");
}


template <MInt nDim, class ppType>
void PostProcessing<nDim, ppType>::initProbeLinePeriodic() {
  TRACE();

  if(!pp->solver().grid().isActive()) return;

  mAlloc(m_probeLinePeriodic, m_noProbeLines, "m_probeLinePeriodic", AT_);

  for(MInt probeLineId = 0; probeLineId < m_noProbeLines; probeLineId++) { // read all probe line directions
    m_probeLinePeriodic[probeLineId] =
        Context::getSolverProperty<MInt>("pp_probeLinePeriodic", m_postprocessingId, AT_, probeLineId);
  }

  for(MInt probeLineId = 0; probeLineId < m_noProbeLines; probeLineId++) { // write all lines to log
    m_log << "probeLine in periodic direction " << m_probeLinePeriodic[probeLineId] << " with probeLinePeriodic ";
    m_log << endl;
  }

  const MInt noVars = m_noVariables    // primitive variables
                      + 3 * (nDim - 1) // Reynolds stress components
                      + 1;             // pressure amplitude p'

  ScratchSpace<MInt> lineDir(m_noProbeLines, (nDim - 1), AT_, "lineDir");

  m_noProbeLineAverageSteps = 0;

  mAlloc(m_probeLineAverageDirection, m_noProbeLines, "m_probeLineAverageDirection", AT_);

  vector<vector<MFloat>> probeLineAverageCoordinates_temp{};

  MInt id = -1;
  for(MInt probeLineId = 0; probeLineId < m_noProbeLines; probeLineId++) { // loop over all probe lines
    m_probeLineAverageDirection[probeLineId] = m_probeLinePeriodic[probeLineId];
    vector<MFloat> temp{};
    vector<MInt> sign{};
    for(MInt i = 0; i < m_noProbeLineIds[probeLineId]; i++) {
      MInt probeId = m_probeLineIds[probeLineId][i];
      id++;
      // periodic direction has to be z-direction
      MFloat pos = m_gridProxy->tree().coordinate(probeId, m_probeLineDirection[probeLineId]);
      temp.emplace_back(pos);
    }

    probeLineAverageCoordinates_temp.emplace_back(temp);
  }

  // Assumes constant mesh in periodic direction
  mAlloc(m_probeLineAverageCoordinates, m_noProbeLines, m_globalNoProbeLineIds[0], "m_probeLineAverageCoordinates", F0,
         AT_);

  for(MInt probeLineId = 0; probeLineId < m_noProbeLines; probeLineId++) {
    MInt noDomain = pp->solver().noDomains();
    ScratchSpace<MInt> recvbuf(noDomain, "recvbuf", AT_);
    recvbuf.fill(0);

    MPI_Gather(&m_noProbeLineIds[probeLineId], 1, MPI_INT, &recvbuf[0], 1, MPI_INT, 0, pp->solver().mpiComm(), AT_,
               "m_noProbeLineIds[probeLineId]", "ppblock()->noDomain()");

    ScratchSpace<MInt> displs(noDomain, "displs", AT_);
    if(pp->solver().domainId() == 0) {
      MInt offset = 0;
      for(MInt dom = 0; dom < noDomain; dom++) {
        displs[dom] = offset;
        offset += recvbuf[dom];
      }
    }

    MPI_Gatherv(&probeLineAverageCoordinates_temp[probeLineId][0], m_noProbeLineIds[probeLineId], MPI_DOUBLE,
                &m_probeLineAverageCoordinates[probeLineId][0], &recvbuf[pp->solver().domainId()],
                &displs[pp->solver().domainId()], MPI_DOUBLE, 0, pp->solver().mpiComm(), AT_,
                "probeLineAverageCoordinates_temp", "m_probeLineAverageCoordinates");

    MPI_Bcast(m_probeLineAverageCoordinates[probeLineId], m_globalNoProbeLineIds[probeLineId], MPI_DOUBLE, 0,
              pp->solver().mpiComm(), AT_, "m_probeLineAverageCoordinates");
  }

  mAlloc(m_noProbeLineAverageIds, m_noProbeLines, m_globalNoProbeLineIds[0], "m_noProbeLineAverageIds", 0, AT_);
  mAlloc(m_probeLineAverageIds, m_noProbeLines, m_globalNoProbeLineIds[0], "m_probeLineAverageIds", AT_);
  mAlloc(m_probeLineAveragePositions, m_noProbeLines, m_globalNoProbeLineIds[0], "m_probeLineAveragePositions", AT_);
  mAlloc(m_globalNoProbeLineAverageIds, m_noProbeLines, m_globalNoProbeLineIds[0], "m_globalNoProbeLineAverageIds", 0,
         AT_);
  mAlloc(m_globalProbeLineAverageVars, m_noProbeLines, m_globalNoProbeLineIds[0], "m_globalProbeLineAverageVars", AT_);


  for(MInt probeLineId = 0; probeLineId < m_noProbeLines; probeLineId++) {
    m_probeLineAverageDirection[probeLineId] = m_probeLinePeriodic[probeLineId];

    for(MInt p = 0; p < m_globalNoProbeLineIds[probeLineId]; p++) {
      TERMM_IF_COND(nDim == 2, "FIXME: the size of periodicCoordiantes was nDim-1, for 2D the below code was writing "
                               "out of array bounds. Check if the code works as intended.");
      vector<MFloat> periodicCoordiantes(nDim);

      periodicCoordiantes[0] = m_probeLineCoordinates[probeLineId][0];
      periodicCoordiantes[1] = m_probeLineAverageCoordinates[probeLineId][p];

      // currently only possible for periodic in z-direction
      lineDir(probeLineId, 0) = m_probeLineDirection[probeLineId] == 0 ? 1 : 0;
      lineDir(probeLineId, 1) = m_probeLineDirection[probeLineId] == 0 ? 0 : 1;

      map<MFloat, MInt, coord_comp_1d_> coordinates;
      const MFloat* coord;
      MFloat halfCellLength;
      MBool ok, ok2, nghbrFound;

      for(MInt cellId = 0; cellId < pp->solver().grid().noInternalCells(); cellId++) {
        if(m_gridProxy->tree().noChildren(cellId) > 0) continue; // skip non leaf cells
        coord = &pp->solver().a_coordinate(cellId, 0);
        ok = true;
        ok2 = true;
        nghbrFound = false;

        halfCellLength = m_gridProxy->halfCellLength(cellId);

        for(MInt i = 0; i < nDim - 1; i++) {
          // check if distance is greater than half of the cell length (in direction i)
          if(abs(periodicCoordiantes[i] - coord[lineDir(probeLineId, i)]) >= halfCellLength) ok = 0;

          // check if distance is a bit more than half of the cell length (line between cells)
          if(abs(periodicCoordiantes[i] - coord[lineDir(probeLineId, i)]) < (halfCellLength + MFloatEps)) {
            for(MInt dirId = 0; dirId < 2 * nDim; dirId++) { // check all directions (x-,x+,...)
              if(dirId == 2 * m_probeLineAverageDirection[probeLineId]
                 || dirId == 2 * m_probeLineAverageDirection[probeLineId] + 1) { // skip line directions
                continue;
              }
              if(m_gridProxy->tree().hasNeighbor(cellId, dirId)
                 != 1) { // skip if cells has no equal level neighbor in direction dirId
                continue;
              }
              MInt nId = m_gridProxy->tree().neighbor(cellId, dirId);
              // check if neighbor cell is already in the probe line
              if(coordinates.find(pp->solver().a_coordinate(nId, m_probeLineAverageDirection[probeLineId]))
                 != coordinates.end()) {
                nghbrFound = true;
              }
            }
          } else {
            ok2 = false;
          }
        }
        if(ok || (ok2 && !nghbrFound)) {
          if(pp->solver().a_isBndryGhostCell(cellId)) continue;
          coordinates[coord[m_probeLineAverageDirection[probeLineId]]] = cellId;
        }
      }

      m_noProbeLineAverageIds[probeLineId][p] = coordinates.size(); // local number of cells

      if(m_noProbeLineAverageIds[probeLineId][p] > 0) {
        mAlloc(m_probeLineAverageIds[probeLineId][p], m_noProbeLineAverageIds[probeLineId][p],
               "m_probeLineAverageIds[probeLineId][p]", AT_);
        mAlloc(m_probeLineAveragePositions[probeLineId][p], m_noProbeLineAverageIds[probeLineId][p],
               "m_probeLineAveragePositions[probeLineId][p]", AT_);
      } else {
        // Allocate dummy array to avoid errors when dereferencing during writing
        mAlloc(m_probeLineAveragePositions[probeLineId][p], 1, "m_probeLineAveragePositions", AT_);
      }

      auto it = coordinates.begin();
      for(MInt i = 0; it != coordinates.end(); it++, i++) {
        m_probeLineAverageIds[probeLineId][p][i] = (*it).second;
        m_probeLineAveragePositions[probeLineId][p][i] = (*it).first;
      }

      MPI_Allreduce(&m_noProbeLineAverageIds[probeLineId][p], &m_globalNoProbeLineAverageIds[probeLineId][p], 1,
                    MPI_INT, MPI_SUM, pp->solver().mpiComm(), AT_, "m_noProbeLineAverageIds",
                    "m_globalNoProbeLineAverageIds");

      mAlloc(m_globalProbeLineAverageVars[probeLineId][p], noVars, "m_globalProbeLineAverageVars[probeLineId][p]", F0,
             AT_);

      for(MInt varId = 0; varId < m_noVariables; varId++) {
        m_globalProbeLineAverageVars[probeLineId][p][varId] = F0;
      }
    }
  }
}

/// \brief Read the time-step properties for slice probing.
template <MInt nDim, class ppType>
void PostProcessing<nDim, ppType>::initTimeStepPropertiesLine() {
  TRACE();

  /*! \page propertiesPP
    \section pp_probeLineInterval
    <code>MInt PostProcessing::m_probeLineInterval</code>\n
    default = <code>none</code>\n\n
    This property determines the time step interval used for line probing.\n
    Keywords: <i>POSTPROCESSING, PROBELINE</i>
  */
  m_probeLineInterval = Context::getSolverProperty<MInt>("pp_probeLineInterval", m_postprocessingId, AT_);
  TERMM_IF_COND(m_probeLineInterval <= 0, "Error: property 'pp_probeLineInterval' needs to be > 0");

  /*! \page propertiesPP
    <code>MInt PostProcessing::m_probeLineStartTimestep</code>\n\n
    This property determines the start timestep used for line probing.\n
    default = <code>0</code>\n
    Keywords: <i>POSTPROCESSING, PROBELINE</i>
  */
  m_probeLineStartTimestep = 0;
  m_probeLineStartTimestep =
      Context::getSolverProperty<MInt>("pp_probeLineStartTimestep", m_postprocessingId, AT_, &m_probeLineStartTimestep);

  /*! \page propertiesPP
    \section pp_probeLineStopTimestep
    <code>MInt PostProcessing::m_probeLineStopTimestep</code>\n
    default = <code>numeric_limits<MInt>::max()</code>\n\n
    This property determines the stop timestep used for line probing.\n
    Keywords: <i>POSTPROCESSING, PROBELINE</i>
  */
  m_probeLineStopTimestep = std::numeric_limits<MInt>::max();
  m_probeLineStopTimestep =
      Context::getSolverProperty<MInt>("pp_probeLineStopTimestep", m_postprocessingId, AT_, &m_probeLineStopTimestep);

  TERMM_IF_COND(m_probeLineStopTimestep <= 0 || m_probeLineStopTimestep <= m_probeLineStartTimestep,
                "Error: pp_probeLineStopTimestep needs to be > pp_probeLineStartTimestep.");

  m_log << "Line probing: interval = " << m_probeLineInterval << "; startTimeStep = " << m_probeLineStartTimestep
        << "; stopTimeStep = " << m_probeLineStopTimestep << std::endl;
}

/**
 * \fn void PostProcessing<nDim,ppType>::initProbeSlice()
 * \brief Initializes properties and memory for slice probing
 *
 * \author A. Niemoeller
 * \date 09.12.2013
 *
 * reads properties probeSliceDir1 / probeSliceDir2 and probeSliceCoordinate which define a slice (or multiple slices)
 *in the specified directions with fixed coordinate in the remaining dimension \n assembles all cell Ids for probing and
 *communicates data for later slice probing \n
 *
 * \param[in] grid pointer to the grid
 *
 **/
template <MInt nDim, class ppType>
void PostProcessing<nDim, ppType>::initProbeSlice() {
  TRACE();

  /*! \page propertiesPP
    \section pp_sliceAiaFileFormat
    <code>MInt PostprocesssingSolver::m_sliceAiaFileFormat</code>\n
    default = <code>0</code>\n\n
    This property determines if the slice files are written in the Aia File Format and therefore readable by the Aia
    Paraview Reader. \n <ul> <li><code>0</code> deactivated</li> <li><code>1</code> active</li>
    </ul>\n
    Keywords: <i>GLOBAL, POSTPROCESSING</i>
  */
  m_sliceAiaFileFormat = false;
  m_sliceAiaFileFormat =
      Context::getSolverProperty<MBool>("pp_sliceAiaFileFormat", m_postprocessingId, AT_, &m_sliceAiaFileFormat);

  /*! \page propertiesPP
    \section _optimizedSliceIo
    <code>MBool PostprocesssingSolver::m_optimizedSliceIo</code>\n
    default = <code>0</code>\n\n
    This property accelerates the slice IO for big files.\n
    If this option is not working for a specific case, set it to false.\n
    </li> <li><code>1</code> active</li>
    </ul>\n
    Keywords: <i>GLOBAL, POSTPROCESSING</i>
  */
  m_optimizedSliceIo = true;
  m_optimizedSliceIo = Context::getBasicProperty<MBool>("optimizedSliceIo", AT_, &m_optimizedSliceIo);

  const MBool activeRank = pp->solver().grid().isActive();
  if(!activeRank && !m_sliceAiaFileFormat) return; // Note: createGridSlice needs to be called from all ranks!

  IF_CONSTEXPR(nDim == 2) {
    m_log << "slice probing for 2D not useful" << endl;
    return;
  }
  // Implementation of CartesianGrid::createGridSlice()
  if(m_sliceAiaFileFormat) {
    /*! \page propertiesPP
      \section pp_sliceAxis
      <code>MString * PostprocesssingSolver::m_sliceAxis</code>\n\n
      This property determines the axis of the probe slice(s).
      <ul>
      <li><code>x</code> x axis</li>
      <li><code>y</code> y axis</li>
      <li><code>z</code> z axis</li>
      </ul>\n
      Keywords: <i>GLOBAL, POSTPROCESSING</i>
    */
    m_noProbeSlices = Context::propertyLength("sliceAxis", m_postprocessingId);
    TERMM_IF_COND(m_noProbeSlices < 1, "Error: slice probing enabled but property sliceAxis not defined.");

    mAlloc(m_sliceAxis, m_noProbeSlices, "m_sliceAxis", AT_);
    mAlloc(m_sliceIntercept, m_noProbeSlices, "m_sliceIntercept", AT_);

    for(MInt sliceId = 0; sliceId < m_noProbeSlices; sliceId++) {
      m_sliceAxis[sliceId] = Context::getSolverProperty<MString>("sliceAxis", m_postprocessingId, AT_, sliceId);
    }
    /*! \page propertiesPP
      \section pp_probeSliceDir
      <code>MFloat * PostprocesssingSolver::m_probeSliceDir</code>\n\n
      This property determines the x, y or z-coordinates of the axes. \n
      Possible values are: all coordinates inside the domain \n
      Values that are exactly between two cell centers don't work (e.g. 0)
      Keywords: <i>GLOBAL, POSTPROCESSING</i>
    */
    for(MInt sliceId = 0; sliceId < m_noProbeSlices; sliceId++) {
      m_sliceIntercept[sliceId] =
          Context::getSolverProperty<MFloat>("sliceIntercept", m_postprocessingId, AT_, sliceId);
    }

    m_sliceVarIds.clear();
    m_noSliceVars.clear();
    m_sliceVarNames.clear();
    pp->solver().getSolverSamplingProperties(m_sliceVarIds, m_noSliceVars, m_sliceVarNames, "Slice");
    // Allocate additional storage in the solver for the sampling variables if required
    pp->solver().initSolverSamplingVariables(m_sliceVarIds, m_noSliceVars);

    if(!activeRank) {
      for(MInt sliceId = 0; sliceId < m_noProbeSlices; sliceId++) {
        stringstream gridName;
        gridName << "sliceGrid_" << m_sliceAxis[sliceId] << "_" << m_sliceIntercept[sliceId];

        const MString filePath = pp->solver().outputDir() + gridName.str();
        pp->solver().grid().raw().createGridSlice(m_sliceAxis[sliceId], m_sliceIntercept[sliceId], filePath, -1,
                                                  nullptr, nullptr, nullptr, nullptr, nullptr, nullptr);
      }
      // Return for inactive ranks after participating in creating the grid-slice of the whole grid (all ranks)
      return;
    }

    mAlloc(m_probeSliceIds, m_noProbeSlices, "m_probeSliceIds", AT_);
    mAlloc(m_noProbeSliceIds, m_noProbeSlices, "m_noProbeSliceIds", AT_);
    mAlloc(m_probeSliceGridNames, m_noProbeSlices, "m_probeSliceGridNames", AT_);
    mAlloc(m_globalNoProbeSliceIds, m_noProbeSlices, "m_globalNoProbeSliceIds", AT_); // global #ids in each slice
    mAlloc(m_noGlobalProbeSliceIds, m_noProbeSlices, pp->solver().m_noDomains, "m_noGlobalProbeSliceIds", AT_);
    // Variables for every slice
    // Number of HilbertIds on this domain determines number of functions calls for data writing
    mAlloc(m_noProbeSliceNoHilbertIds, m_noProbeSlices, "m_noProbeSliceNoHilbertIds", AT_);
    // Store hilbertInfo, which is hilbertId, noCellsPerHilbertId (for offset in local data array) and the offset for
    // file writing
    mAlloc(m_noProbeSliceHilbertInfo, m_noProbeSlices, "m_noProbeSliceHilbertInfo", AT_);
    // Maximum number of HilbertIds is used to call the function for writing data on every domain equally often
    mAlloc(m_noProbeSliceMaxNoHilbertIds, m_noProbeSlices, "m_noProbeSliceMaxNoHilbertIds", AT_);

    if(m_optimizedSliceIo) {
      mAlloc(m_noProbeSliceNoContHilbertIds, m_noProbeSlices, "m_noProbeSliceNoContHilbertIds", AT_);
      mAlloc(m_noProbeSliceContHilbertInfo, m_noProbeSlices, "m_noProbeSliceContHilbertInfo", AT_);
      mAlloc(m_noProbeSliceMaxNoContHilbertIds, m_noProbeSlices, "m_noProbeSliceMaxNoContHilbertIds", AT_);
    }

    MIntScratchSpace tmpSliceCellIds(pp->solver().grid().noInternalCells(), AT_, "tmpSliceCellIds");
    MIntScratchSpace tmpSliceHilbertInfo(pp->solver().grid().noInternalCells() * 3, AT_, "tmpSliceHilbertInfo");

    const MInt contHInfoSize = (m_optimizedSliceIo) ? pp->solver().grid().noInternalCells() * 3 : 0;
    MIntScratchSpace tmpSliceContHilbertInfo(contHInfoSize, AT_, "tmpSliceContHilbertInfo");

    for(MInt i = 0; i < m_noProbeSlices; i++) {
      m_probeSliceIds[i] = nullptr;
      m_noProbeSliceHilbertInfo[i] = nullptr;
      m_noProbeSliceNoHilbertIds[i] = 0;

      if(m_optimizedSliceIo) {
        m_noProbeSliceContHilbertInfo[i] = nullptr;
        m_noProbeSliceNoContHilbertIds[i] = 0;
      }
    }

    for(MInt sliceId = 0; sliceId < m_noProbeSlices; sliceId++) {
      stringstream gridName;

      gridName << "sliceGrid_" << m_sliceAxis[sliceId] << "_" << m_sliceIntercept[sliceId];
      m_probeSliceGridNames[sliceId] = gridName.str();

      MInt tmpNoSliceCellIds = 0;
      MInt tmpNoHilbertIds = 0;
      MInt tmpNoContHilbertIds = 0;

      MInt* tmpNoContHilbertIdsP = nullptr;
      MInt* tmpSliceContHilbertInfoP = nullptr;
      if(m_optimizedSliceIo) {
        tmpNoContHilbertIdsP = &tmpNoContHilbertIds;
        tmpSliceContHilbertInfoP = &tmpSliceContHilbertInfo[0];
      }

      const MString filePath = pp->solver().outputDir() + m_probeSliceGridNames[sliceId];
      pp->solver().grid().raw().createGridSlice(
          m_sliceAxis[sliceId], m_sliceIntercept[sliceId], filePath, pp->solver().solverId(), &tmpNoSliceCellIds,
          tmpSliceCellIds.getPointer(), &tmpNoHilbertIds, tmpSliceHilbertInfo.getPointer(), tmpNoContHilbertIdsP,
          tmpSliceContHilbertInfoP);

      gridName << ParallelIo::fileExt();
      m_probeSliceGridNames[sliceId] = gridName.str();

      m_noProbeSliceIds[sliceId] = tmpNoSliceCellIds;
      m_noProbeSliceNoHilbertIds[sliceId] = tmpNoHilbertIds;

      if(tmpNoSliceCellIds > 0) mAlloc(m_probeSliceIds[sliceId], tmpNoSliceCellIds, "m_probeSliceIds", AT_);

      for(MInt slicedCellId = 0; slicedCellId < tmpNoSliceCellIds; slicedCellId++) {
        // Convert grid cellId into block cellId
        m_probeSliceIds[sliceId][slicedCellId] = pp->solver().grid().tree().grid2solver(tmpSliceCellIds[slicedCellId]);
      }

      if(tmpNoHilbertIds > 0) {
        mAlloc(m_noProbeSliceHilbertInfo[sliceId], tmpNoHilbertIds * 3, "m_noProbeSliceHilbertInfo", AT_);
        copy_n(&tmpSliceHilbertInfo[0], tmpNoHilbertIds * 3, &m_noProbeSliceHilbertInfo[sliceId][0]);
      }

      if(m_optimizedSliceIo) {
        m_noProbeSliceNoContHilbertIds[sliceId] = tmpNoContHilbertIds;
        if(tmpNoContHilbertIds > 0) {
          mAlloc(m_noProbeSliceContHilbertInfo[sliceId], tmpNoContHilbertIds * 3, "m_noProbeSliceContHilbertInfo", AT_);
        }
        copy_n(&tmpSliceContHilbertInfo[0], tmpNoContHilbertIds * 3, &m_noProbeSliceContHilbertInfo[sliceId][0]);
      }
    }
    MPI_Allreduce(m_noProbeSliceIds, m_globalNoProbeSliceIds, m_noProbeSlices, MPI_INT, MPI_SUM, pp->solver().mpiComm(),
                  AT_, "m_noProbeSliceIds", "m_globalNoProbeSliceIds");
    MPI_Allreduce(m_noProbeSliceNoHilbertIds, m_noProbeSliceMaxNoHilbertIds, m_noProbeSlices, MPI_INT, MPI_MAX,
                  pp->solver().mpiComm(), AT_, "m_noProbeSliceNoHilbertIds", "m_noProbeSliceMaxNoHilbertIds");
    if(m_optimizedSliceIo) {
      MPI_Allreduce(m_noProbeSliceNoContHilbertIds, m_noProbeSliceMaxNoContHilbertIds, m_noProbeSlices, MPI_INT,
                    MPI_MAX, pp->solver().mpiComm(), AT_, "m_noProbeSliceNoContHilbertIds",
                    "m_noProbeSliceMaxNoContHilbertIds");
    }
  } else {
    /*! \page propertiesPP
      \section pp_probeSliceDir
      <code>MInt PostprocesssingSolver::m_probeSliceDir</code>\n\n
      This property determines the directions of the probe slice(s).
      <ul>
      <li><code>0</code> x</li>
      <li><code>1</code> y</li>
      <li><code>2</code> z</li>
      </ul>\n
      Keywords: <i>GLOBAL, POSTPROCESSING</i>
    */
    if(Context::propertyLength("pp_probeSliceDir", m_postprocessingId) % 2 != 0) {
      m_log << "odd number of probe slice directions provided by probeSliceDir, should be even!" << endl;
    }
    m_noProbeSlices = Context::propertyLength("pp_probeSliceDir", m_postprocessingId) / 2;
    if(m_noProbeSlices == 0) {
      m_log << "not enough probe slice directions provided by probeSliceDir!" << endl;
      return;
    }

    mAlloc(m_probeSliceDir, 2 * m_noProbeSlices, "m_probeSliceDir", AT_);
    mAlloc(m_probeSliceCoordinate, m_noProbeSlices, "m_probeSliceDir", AT_);

    ScratchSpace<MInt> dir(m_noProbeSlices, AT_, "dir");
    for(MInt dirId = 0; dirId < 2 * m_noProbeSlices; dirId++) {
      m_probeSliceDir[dirId] = Context::getSolverProperty<MInt>("pp_probeSliceDir", m_postprocessingId, AT_, dirId);
    }
    for(MInt sliceId = 0; sliceId < m_noProbeSlices; sliceId++) {
      if((m_probeSliceDir[2 * sliceId] < 0 || m_probeSliceDir[2 * sliceId] >= nDim)
         || (m_probeSliceDir[2 * sliceId + 1] < 0 || m_probeSliceDir[2 * sliceId + 1] >= nDim)
         || (m_probeSliceDir[2 * sliceId] == m_probeSliceDir[2 * sliceId + 1])) {
        m_log << "invalid value for property probeSliceDir for slice #" << sliceId << ": "
              << m_probeSliceDir[2 * sliceId] << " " << m_probeSliceDir[2 * sliceId + 1] << endl;
        return;
      }
      if(m_probeSliceDir[2 * sliceId] > m_probeSliceDir[2 * sliceId + 1]) {
        MInt tmp = m_probeSliceDir[2 * sliceId];
        m_probeSliceDir[2 * sliceId] = m_probeSliceDir[2 * sliceId + 1];
        m_probeSliceDir[2 * sliceId + 1] = tmp;
      }
      dir[sliceId] = 3 - (m_probeSliceDir[2 * sliceId] + m_probeSliceDir[2 * sliceId + 1]);
    }

    /*! \page propertiesPP
    \section pp_probeSliceCoordinate
    <code>MInt PostprocesssingSolver::m_probeSliceCoordinate</code>\n\n
    default = <code>empty</code>\n
    This property determines the coordinate of the probe slice(s) in the remaining direction (size:
    number of probe slices).
    <ul>
    <li><code>coordinate</code></li>
    </ul>
    Keywords: <i>POSTPROCESSING</i>
    */
    if(Context::propertyLength("pp_probeSliceCoordinate", m_postprocessingId) < m_noProbeSlices) {
      m_log << "to few #probeSliceCoordinate provided, found "
            << Context::propertyLength("pp_probeSliceCoordinate", m_postprocessingId) << " needed " << m_noProbeSlices
            << endl
            << "using default value 0.0 for missing values" << endl;
    }
    for(MInt sliceId = 0; sliceId < m_noProbeSlices; sliceId++) {
      if(sliceId < Context::propertyLength("pp_probeSliceCoordinate", m_postprocessingId)) {
        m_probeSliceCoordinate[sliceId] =
            Context::getSolverProperty<MFloat>("pp_probeSliceCoordinate", m_postprocessingId, AT_, sliceId);
      } else {
        m_probeSliceCoordinate[sliceId] = 0.0;
      }

      m_log << "probeSlice #" << sliceId << " in directions " << m_probeSliceDir[2 * sliceId] << ", "
            << m_probeSliceDir[2 * sliceId + 1] << " with probeSliceCoordinate " << m_probeSliceCoordinate[sliceId]
            << endl;
    }

    mAlloc(m_noProbeSliceIds, m_noProbeSlices, "m_noProbeSliceIds", AT_); // local #cells in each slice
    mAlloc(m_probeSliceIds, m_noProbeSlices, "m_probeSliceIds", AT_);     // collector for cellIds for each slice
    mAlloc(m_probeSlicePositions, m_noProbeSlices, "m_probeSlicePositions",
           AT_); // collector for cell positions in slice
    mAlloc(m_probeSliceOffsets, m_noProbeSlices, "m_probeSliceOffsets", AT_);
    mAlloc(m_globalNoProbeSliceIds, m_noProbeSlices, "m_globalNoProbeSliceIds", AT_); // global #ids in each slice
    mAlloc(m_noGlobalProbeSliceIds, m_noProbeSlices, pp->solver().m_noDomains, "m_noGlobalProbeSliceIds", AT_);

    // suppress valgrind error
    for(MInt i = 0; i < m_noProbeSlices; i++) {
      m_probeSliceIds[i] = nullptr;
      m_probeSlicePositions[i] = nullptr;
    }

    const MFloat* coord = nullptr;
    const MFloat* coord2 = nullptr;
    MFloat halfCellLength = NAN;
    MInt dirId = 0;
    MInt nId = 0;
    MBool nghbrFound = false;
    for(MInt sliceId = 0; sliceId < m_noProbeSlices; sliceId++) { // loop over all slices
      map<pair<MFloat, MFloat>, MInt, coord_comp_2d_> coordinates{};
      for(MInt cellId = 0; cellId < pp->solver().grid().noInternalCells(); cellId++) { // collect cellIds on slice
        if(m_gridProxy->tree().noChildren(cellId) > 0) continue;                       // skip non leaf cells
        coord = &pp->solver().a_coordinate(cellId, 0);

        halfCellLength = m_gridProxy->halfCellLength(cellId);
        if(abs(m_probeSliceCoordinate[sliceId] - coord[dir[sliceId]]) < halfCellLength) {
          coordinates[pair<MFloat, MFloat>(coord[m_probeSliceDir[2 * sliceId]],
                                           coord[m_probeSliceDir[2 * sliceId + 1]])] = cellId;
        } else {
          nghbrFound = false;
          // check if distance is a bit more than half of the cell length (slice between cells)
          if(abs(m_probeSliceCoordinate[sliceId] - coord[dir[sliceId]]) < (halfCellLength + MFloatEps)) {
            for(MInt i = 0; i < 2; i++) {
              dirId = 2 * dir[sliceId] + i; // directions orthogonal to slice
              if(m_gridProxy->tree().hasNeighbor(cellId, dirId)
                 != 1) { // skip if cells has no equal level neighbor in direction dirId
                continue;
              }
              nId = m_gridProxy->tree().neighbor(cellId, dirId);
              coord2 = &pp->solver().a_coordinate(nId, 0);
              // check if neighbor cell is already in the probe slice
              if(coordinates.find(pair<MFloat, MFloat>(coord2[m_probeSliceDir[2 * sliceId]],
                                                       coord2[m_probeSliceDir[2 * sliceId + 1]]))
                 != coordinates.end()) {
                nghbrFound = true;
              }
            }
            if(!nghbrFound) {
              coordinates[pair<MFloat, MFloat>(coord[m_probeSliceDir[2 * sliceId]],
                                               coord[m_probeSliceDir[2 * sliceId + 1]])] = cellId;
            }
          }
        }
      }

      m_noProbeSliceIds[sliceId] = coordinates.size();

      if(m_noProbeSliceIds[sliceId] > 0) {
        mAlloc(m_probeSliceIds[sliceId], m_noProbeSliceIds[sliceId], "m_probeSliceIds", AT_);
        mAlloc(m_probeSlicePositions[sliceId], 2 * m_noProbeSliceIds[sliceId], "m_probeSlicePositions", AT_);
      } else {
        // Allocate dummy array to avoid errors when dereferencing during writing
        mAlloc(m_probeSlicePositions[sliceId], 2, "m_probeSlicePositions", AT_);
      }

      // testing
      MPI_Gather(&(m_noProbeSliceIds[sliceId]), 1, MPI_INT, m_noGlobalProbeSliceIds[sliceId], 1, MPI_INT, 0,
                 pp->solver().mpiComm(), AT_, "(m_noProbeSliceIds[sliceId])", "m_noGlobalProbeSliceIds[sliceId]");
      // m_log << "slice #" << sliceId << ": ";
      // for(MInt i=0; i<pp->solver().m_noDomains; i++){
      //  m_log << m_noGlobalProbeSliceIds[sliceId][i] << " ";
      //}
      // m_log << endl;
      //

      auto it = coordinates.begin();
      for(MInt i = 0; it != coordinates.end(); it++, i++) {
        m_probeSliceIds[sliceId][i] = (*it).second;
        m_probeSlicePositions[sliceId][2 * i] = (*it).first.first;
        m_probeSlicePositions[sliceId][2 * i + 1] = (*it).first.second;
      }

      ParallelIo::size_type localCount, offset, totalCount;
      localCount = m_noProbeSliceIds[sliceId];
      ParallelIo::calcOffset(localCount, &offset, &totalCount, pp->solver().mpiComm());
      m_probeSliceOffsets[sliceId] = offset;
    }

    MPI_Allreduce(m_noProbeSliceIds, m_globalNoProbeSliceIds, m_noProbeSlices, MPI_INT, MPI_SUM, pp->solver().mpiComm(),
                  AT_, "m_noProbeSliceIds", "m_globalNoProbeSliceIds");
  }
}


/// \brief Read the time-step properties for slice probing.
template <MInt nDim, class ppType>
void PostProcessing<nDim, ppType>::initTimeStepPropertiesSlice() {
  TRACE();

  /*! \page propertiesPP
    \section pp_probeSliceInterval
    <code>MInt PostProcessing::m_probeSliceInterval</code>\n
    default = <code>none</code>\n\n
    This property determines the time step interval used for averaging solutions.
    <ul>
    <li><code>interval</code> </li>
    </ul>
    Keywords: <i>POSTPROCESSING, SLICE</i>
  */
  m_probeSliceInterval = Context::getSolverProperty<MInt>("pp_probeSliceInterval", m_postprocessingId, AT_);
  TERMM_IF_COND(m_probeSliceInterval <= 0, "Error: property 'pp_probeSliceInterval' needs to be > 0");

  /*! \page propertiesPP
    \section pp_probeSliceStartTimestep
    <code>MInt PostProcessing::m_probeSliceStartTimestep</code>\n
    default = <code>0</code>\n\n
    This property determines the start timestep used for slice probing.
    <ul>
    <li><code>timestep</code> </li>
    </ul>
    Keywords: <i>POSTPROCESSING, SLICE</i>
  */
  m_probeSliceStartTimestep = 0;
  m_probeSliceStartTimestep = Context::getSolverProperty<MInt>("pp_probeSliceStartTimestep", m_postprocessingId, AT_,
                                                               &m_probeSliceStartTimestep);

  /*! \page propertiesPP
    \section pp_probeSliceStopTimestep
    <code>MInt PostProcessing::m_probeSliceStopTimestep</code>\n
    default = <code>numeric_limits<MInt>::max()</code>\n\n
    This property determines the stop timestep used for slice probing.
    <ul>
    <li><code>timestep</code> </li>
    </ul>
    Keywords: <i>POSTPROCESSING</i>
  */
  m_probeSliceStopTimestep = std::numeric_limits<MInt>::max();
  m_probeSliceStopTimestep =
      Context::getSolverProperty<MInt>("pp_probeSliceStopTimestep", m_postprocessingId, AT_, &m_probeSliceStopTimestep);
  TERMM_IF_COND(m_probeSliceStopTimestep <= 0 || m_probeSliceStopTimestep <= m_probeSliceStartTimestep,
                "Error: pp_probeSliceStopTimestep needs to be > pp_probeSliceStartTimestep.");

  m_log << "Slice probing: interval = " << m_probeSliceInterval << "; startTimeStep = " << m_probeSliceStartTimestep
        << "; stopTimeStep = " << m_probeSliceStopTimestep << std::endl;
}


/** \brief Initializes properties for slice averaging.
 *
 * \author D. Zilles
 * \date 07.09.2017
 *
 * loads sliceAxes and sliceIntercepts. \n
 *
 **/
template <MInt nDim, class ppType>
void PostProcessing<nDim, ppType>::initAverageSolutionsSlice() {
  TRACE();

  m_noProbeSlices = Context::propertyLength("sliceAxis", m_postprocessingId);

  mAlloc(m_sliceAxis, m_noProbeSlices, "m_sliceAxis", AT_);
  mAlloc(m_sliceIntercept, m_noProbeSlices, "m_sliceIntercept", AT_);

  for(MInt sliceId = 0; sliceId < m_noProbeSlices; sliceId++) {
    m_sliceAxis[sliceId] = Context::getSolverProperty<MString>("sliceAxis", m_postprocessingId, AT_, sliceId);
  }
  for(MInt sliceId = 0; sliceId < m_noProbeSlices; sliceId++) {
    m_sliceIntercept[sliceId] = Context::getSolverProperty<MFloat>("sliceIntercept", m_postprocessingId, AT_, sliceId);
  }
}


// TODO labels:PP
/** \brief initializes data for spatial averaging
 *
 * \author A. Niemoeller
 * \date 23.04.14
 *
 * \param[in] grid pointer to the grid
 *
 **/
template <MInt nDim, class ppType>
void PostProcessing<nDim, ppType>::initSpatialAveraging() {
  TRACE();

  /*! \page propertiesPP
    \section pp_spatialDirection1
    <code>MInt PostprocesssingSolver::m_spatialDirection1</code>\n
    default = <code>-1</code>\n\n
    This property determines the first direction for spatial averaging.\n
    If pp_spatialDirection1 and pp_spatialDirection2 are both -1, spatial averaging of the whole
    flow field on a single point is performed.\n
    If one of the two directions is unequal -1, spatial averaging on a line in the specified
    direction is performed.\n
    If both directions are unequal -1 and not equal to each other, same is done for a slice in these
    two directions (3D only).\n
    <ul>
    <li><code>-1</code> no direction</li>
    <li><code>0</code> x</li>
    <li><code>1</code> y</li>
    <li><code>2</code> z (3D only)</li>
    </ul>\n
    Keywords: <i>POSTPROCESSING</i>
  */
  m_spatialDirection1 = -1;
  m_spatialDirection1 =
      Context::getSolverProperty<MInt>("pp_spatialDirection1", m_postprocessingId, AT_, &m_spatialDirection1);

  /*! \page propertiesPP
    \section pp_spatialDirection2
    <code>MInt PostprocesssingSolver::m_spatialDirection2</code>\n
    default = <code>-1</code>\n\n
    This property determines the second direction for spatial averaging.\n
    See pp_spatialDirection1 for more information.\n
    Keywords: <i>POSTPROCESSING</i>
  */
  m_spatialDirection2 = -1;
  m_spatialDirection2 =
      Context::getSolverProperty<MInt>("pp_spatialDirection2", m_postprocessingId, AT_, &m_spatialDirection2);

  MInt noCells = pp->solver().grid().noInternalCells();
  MInt noVars = m_noVariables;
  MInt dir1 = m_spatialDirection1;
  MInt dir2 = m_spatialDirection2;
  MInt lvlCell;

  // precalculate weights
  for(MInt lvl = 0; lvl < 20; lvl++) {
    m_spatialLvlWeight[lvl] = 1 / pow(2, nDim * lvl);
  }

  if(dir1 == -1 && dir2 == -1) { // single point
    m_spatialWeightSinglePoint = 0;
    for(MInt cellId = 0; cellId < noCells; cellId++) {         // calculate local total weight
      if(m_gridProxy->tree().noChildren(cellId) > 0) continue; // check if leaf cell
      lvlCell = m_gridProxy->tree().level(cellId);
      m_spatialWeightSinglePoint += m_spatialLvlWeight[lvlCell];
    }
  } // point end

  else if((dir1 == -1 && dir2 != -1) || (dir2 == -1 && dir1 != -1)) { // line
    dir1 = (dir2 == -1) ? dir1 : dir2;

    if(dir1 < 0 || dir1 > nDim - 1) {
      m_log << "invalid spatial direction" << endl;
      return;
    }

    // single logfiles for testing
    stringstream logName;
    logName << "log_rank_" << pp->solver().grid().domainId();
    ofstream log;
    log.open((logName.str()).c_str());
    log.precision(15);
    //

    map<MFloat, MInt, coord_comp_1d_>& coordinates =
        m_spatialLineCoordinates; // TODO labels:PP,toremove remove reference

    MInt cellLvl;
    for(MInt cellId = 0; cellId < noCells; cellId++) { // collect cell coordinates in spatial direction
      cellLvl = m_gridProxy->tree().level(cellId);
      if(m_gridProxy->tree().noChildren(cellId) == 0) {
        coordinates[pp->solver().a_coordinate(cellId, dir1)] = cellLvl;
      }
    }

    createCellToMap1D(coordinates, m_spatialLineCellToMap); // reduce coordinates to positions where to average and
                                                            // store a mapping of positions to average positions

    log << "local coordinates" << endl;
    for(auto it = coordinates.begin(); it != coordinates.end(); it++) {
      log << (*it).first << "\t" << (*it).second << endl;
    }
    log << endl;

    MInt noSolvers = pp->solver().m_noDomains; // rename

    m_spatialLineNoCells = coordinates.size();
    m_spatialCoordSum = 0;
    MPI_Reduce(&m_spatialLineNoCells, &m_spatialCoordSum, 1, MPI_INT, MPI_SUM, 0, pp->solver().mpiComm(), AT_,
               "m_spatialLineNoCells", "m_spatialCoordSum");

    ScratchSpace<MFloat> coords(m_spatialLineNoCells, AT_, "coords");
    ScratchSpace<MInt> levels(m_spatialLineNoCells, AT_, "levels");
    //    ScratchSpace<MInt> no_coords( noSolvers, AT_, "no_coords" );
    ScratchSpace<MInt> all_levels(m_spatialCoordSum, AT_, "all_levels");

    mAlloc(m_spatialDispls, noSolvers, "m_spatialDispls", AT_);
    mAlloc(m_spatialVarsDispls, noSolvers, "m_spatialVarsDispls", AT_);
    mAlloc(m_spatialRecvcnts, noSolvers, "m_spatialRecvcnts", AT_);
    mAlloc(m_spatialVarsRecvcnts, noSolvers, "m_spatialVarsRecvcnts", AT_);
    mAlloc(m_spatialLineAllCoord, m_spatialCoordSum, "m_spatialLineAllCoord", AT_); // TODO labels:PP root only

    MInt pos = 0;
    for(auto it = coordinates.begin(); it != coordinates.end(); it++, pos++) {
      coords[pos] = (*it).first;
      levels[pos] = (*it).second;
    }

    MPI_Gather(&m_spatialLineNoCells, 1, MPI_INT, m_spatialRecvcnts, 1, MPI_INT, 0, pp->solver().mpiComm(), AT_,
               "m_spatialLineNoCells", "m_spatialRecvcnts");

    log << "m_spatialRecvcnts";
    for(MInt i = 0; i < noSolvers; i++) {
      log << m_spatialRecvcnts[i] << " ";
    }
    log << endl;

    m_spatialDispls[0] = 0;
    log << "displs 0 ";
    for(MInt solverId = 1; solverId < noSolvers; solverId++) {
      m_spatialDispls[solverId] = m_spatialDispls[solverId - 1] + m_spatialRecvcnts[solverId - 1];
      log << m_spatialDispls[solverId] << " ";
    }
    log << endl;

    MPI_Gatherv(coords.begin(), m_spatialLineNoCells, MPI_DOUBLE, m_spatialLineAllCoord, m_spatialRecvcnts,
                m_spatialDispls, MPI_DOUBLE, 0, pp->solver().mpiComm(), AT_, "coords.begin()", "m_spatialLineAllCoord");
    MPI_Gatherv(levels.begin(), m_spatialLineNoCells, MPI_INT, all_levels.begin(), m_spatialRecvcnts, m_spatialDispls,
                MPI_INT, 0, pp->solver().mpiComm(), AT_, "levels.begin()", "all_levels.begin()");

    if(pp->solver().domainId() == 0) {
      for(MInt cellId = 0; cellId < m_spatialCoordSum; cellId++) {
        m_spatialGlobalLineCoordinates[m_spatialLineAllCoord[cellId]] = all_levels[cellId];
      }
      createCellToMap1D(m_spatialGlobalLineCoordinates, m_spatialGlobalLineCellToMap);


      log << "global coordinates " << m_spatialGlobalLineCoordinates.size() << endl;
      for(auto it = m_spatialGlobalLineCoordinates.begin(); it != m_spatialGlobalLineCoordinates.end(); it++) {
        log << (*it).first << "\t" << (*it).second << endl;
      }
      log << endl;
    }

    //    mAlloc( m_spatialLineAllVars, m_spatialCoordSum*(noVars+1), "m_spatialLineAllVars", AT_); // TODO labels:PP
    //    noVars can change when file is loaded
    m_spatialLineAllVars = 0;

    if(pp->solver().domainId() == 0) {
      log << "recvcnts     displs" << endl;
      for(MInt solverId = 0; solverId < noSolvers; solverId++) {
        m_spatialVarsDispls[solverId] = (noVars + 1) * m_spatialDispls[solverId];
        m_spatialVarsRecvcnts[solverId] = (noVars + 1) * m_spatialRecvcnts[solverId];
        log << m_spatialVarsRecvcnts[solverId] << " " << m_spatialVarsDispls[solverId] << endl;
      }
    }

  }      // line end
  else { // plane
    if(dir1 < 0 || dir1 > nDim - 1 || dir2 < 0 || dir2 > nDim - 1 || dir1 == dir2) {
      m_log << "invalid spatial directions " << dir1 << " " << dir2 << endl;
      return;
    }
    if(dir2 < dir1) {
      MInt tmpDir = dir1;
      dir1 = dir2;
      dir2 = tmpDir;
    }

    // single logfiles for testing
    stringstream logName;
    logName << "log_rank_" << pp->solver().grid().domainId();
    ofstream log;
    log.open((logName.str()).c_str());
    log.precision(15);
    //

    map<pair<MFloat, MFloat>, MInt, coord_comp_2d_>& coordinates =
        m_spatialPlaneCoordinates; // TODO labels:PP remove reference

    MInt cellLvl;
    for(MInt cellId = 0; cellId < noCells; cellId++) { // collect cell coordinates in spatial direction
      cellLvl = m_gridProxy->tree().level(cellId);
      if(m_gridProxy->tree().noChildren(cellId) == 0) {
        coordinates[pair<MFloat, MFloat>(pp->solver().a_coordinate(cellId, dir1),
                                         pp->solver().a_coordinate(cellId, dir2))] = cellLvl;
      }
    }

    createCellToMap2D(coordinates, m_spatialPlaneCellToMap); // reduce coordinates to positions where to average and
                                                             // store a mapping of positions to average positions

    /*
        log << "local coordinates" << endl;
        for( map< MFloat,MInt >::const_iterator it=coordinates.begin(); it!=coordinates.end(); it++){
          log << (*it).first << "\t" << (*it).second << endl;
        }
        log << endl;
    */

    MInt noSolvers = pp->solver().m_noDomains; // rename

    m_spatialPlaneNoCells = coordinates.size();
    m_spatialCoordSum = 0;
    MPI_Reduce(&m_spatialPlaneNoCells, &m_spatialCoordSum, 1, MPI_INT, MPI_SUM, 0, pp->solver().mpiComm(), AT_,
               "m_spatialPlaneNoCells", "m_spatialCoordSum");

    ScratchSpace<MFloat> coords(2 * m_spatialPlaneNoCells, AT_, "coords");
    ScratchSpace<MInt> levels(m_spatialPlaneNoCells, AT_, "levels");
    ScratchSpace<MInt> no_coords(noSolvers, AT_, "no_coords");
    ScratchSpace<MInt> coords_displs(noSolvers, AT_, "coords_displs");
    //    ScratchSpace<MInt> no_ids( noSolvers, AT_, "no_ids" );
    ScratchSpace<MInt> all_levels(m_spatialCoordSum, AT_, "all_levels");

    mAlloc(m_spatialDispls, noSolvers, "m_spatialDispls", AT_);
    mAlloc(m_spatialRecvcnts, noSolvers, "m_spatialRecvcnts", AT_);
    mAlloc(m_spatialVarsDispls, noSolvers, "m_spatialVarsDispls", AT_);
    mAlloc(m_spatialVarsRecvcnts, noSolvers, "m_spatialVarsRecvcnts", AT_);
    mAlloc(m_spatialPlaneAllCoord, 2 * m_spatialCoordSum, "m_spatialPlaneAllCoord",
           AT_); // TODO labels:PP root only?

    MInt pos = 0;
    for(auto it = coordinates.begin(); it != coordinates.end(); it++, pos++) {
      coords[2 * pos] = (*it).first.first;
      coords[2 * pos + 1] = (*it).first.second;
      levels[pos] = (*it).second;
    }

    MPI_Gather(&m_spatialPlaneNoCells, 1, MPI_INT, m_spatialRecvcnts, 1, MPI_INT, 0, pp->solver().mpiComm(), AT_,
               "m_spatialPlaneNoCells", "m_spatialRecvcnts");

    for(MInt i = 0; i < noSolvers; i++) {
      no_coords[i] = 2 * m_spatialRecvcnts[i];
    }

    log << "no_coords ";
    for(MInt i = 0; i < noSolvers; i++) {
      log << no_coords[i] << " ";
    }
    log << endl;

    m_spatialDispls[0] = 0;
    coords_displs[0] = 0;
    log << "displs 0 ";
    for(MInt solverId = 1; solverId < noSolvers; solverId++) {
      m_spatialDispls[solverId] = m_spatialDispls[solverId - 1] + no_coords[solverId - 1];
      coords_displs[solverId] = 2 * m_spatialDispls[solverId];
      log << m_spatialDispls[solverId] << " ";
    }
    log << endl;

    MPI_Gatherv(coords.begin(), 2 * m_spatialPlaneNoCells, MPI_DOUBLE, m_spatialPlaneAllCoord, no_coords.begin(),
                coords_displs.begin(), MPI_DOUBLE, 0, pp->solver().mpiComm(), AT_, "coords.begin()",
                "m_spatialPlaneAllCoord");
    MPI_Gatherv(levels.begin(), m_spatialPlaneNoCells, MPI_INT, all_levels.begin(), m_spatialRecvcnts, m_spatialDispls,
                MPI_INT, 0, pp->solver().mpiComm(), AT_, "levels.begin()", "all_levels.begin()");

    if(pp->solver().domainId() == 0) {
      for(MInt cellId = 0; cellId < m_spatialCoordSum; cellId++) {
        m_spatialGlobalPlaneCoordinates[pair<MFloat, MFloat>(
            m_spatialPlaneAllCoord[2 * cellId], m_spatialPlaneAllCoord[2 * cellId + 1])] = all_levels[cellId];
      }
      createCellToMap2D(m_spatialGlobalPlaneCoordinates, m_spatialGlobalPlaneCellToMap);

      /*
            log << "global coordinates " << m_spatialGlobalLineCoordinates.size() << endl;
            for( map< MFloat,MInt >::const_iterator it=m_spatialGlobalLineCoordinates.begin();
         it!=m_spatialGlobalLineCoordinates.end(); it++){ log << (*it).first << "\t" << (*it).second << endl;
            }
            log << endl;
      */
    }

    //    mAlloc( m_spatialPlaneAllVars, m_spatialCoordSum*(noVars+1), "m_spatialPlaneAllVars", AT_ );
    m_spatialPlaneAllVars = 0;

    if(pp->solver().domainId() == 0) {
      log << "recvcnts     displs" << endl;
      for(MInt solverId = 0; solverId < noSolvers; solverId++) {
        m_spatialVarsDispls[solverId] = (noVars + 1) * m_spatialDispls[solverId];
        m_spatialVarsRecvcnts[solverId] = (noVars + 1) * m_spatialRecvcnts[solverId];
        log << m_spatialRecvcnts[solverId] << " " << m_spatialVarsDispls[solverId] << endl;
      }
    }

  } // plane end
}


/**
 * \fn void PostProcessing<nDim,ppType>::initProbeArbitraryLine()
 * \brief initializes properties and data for arbitrary line probing
 *
 * \author Ansgar Niemoeller
 * \date 07.06.2014
 *
 * \param[in] grid pointer to the grid
 *
 **/
template <MInt nDim, class ppType>
void PostProcessing<nDim, ppType>::initProbeArbitraryLine() {
  TRACE();

  MInt noPoints = 0;

  /*! \page propertiesPP
    \section pp_arbLinePoints
    <code>MFloat* PostProcessing::m_arbLinePoints</code>\n\n
    default = <code>0</code>\n
    This property determines the start/end points of arbitrary probe lines.\n
    Specify each point with an x, y and z (3D only) component.\n
    A line is defined by two points and the number of lines is determined by the given total number
    of points.
    <ul>
    <li><code>coordinate</code></li>
    </ul>
    Keywords: <i>POSTPROCESSING</i>
  */
  noPoints = Context::propertyLength("pp_arbLinePoints", m_postprocessingId);
  if(noPoints == 0 || noPoints % (2 * nDim) != 0) {
    m_log << "property pp_arbLinePoints has " << Context::propertyLength("pp_arbLinePoints", m_postprocessingId)
          << " entries, multiple of " << 2 * nDim << " needed";
    m_noArbLines = 0;
    return;
  }
  mAlloc(m_arbLinePoints, noPoints, "m_arbLinePoints", AT_);
  m_noArbLines = noPoints / (2 * nDim);
  for(MInt i = 0; i < noPoints; i++) {
    m_arbLinePoints[i] = Context::getSolverProperty<MFloat>("pp_arbLinePoints", m_postprocessingId, AT_, i);
  }


  mAlloc(m_noArbLineIds, m_noArbLines, "m_noArbLineIds", AT_);

  /*! \page propertiesPP
    \section pp_noPointsArbLine
    <code>MInt* PostProcessing::m_noArbLineIds</code>\n\n
    default = <code>empty</code>\n
    This property determines the number of probing points for each arbitrary line.\n Note: the
    number of points written to disk can be less than the provided number of points, e.g. when the
    line starts/ends outside the computational domain.
    <ul>
    <li><code>integer &ge; 2</code></li>
    </ul>
    Keywords: <i>POSTPROCESSING</i>
  */
  if(Context::propertyLength("pp_noPointsArbLine", m_postprocessingId) < m_noArbLines) {
    m_log << "property pp_noPointsArbLine has " << Context::propertyLength("pp_noPointsArbLine", m_postprocessingId)
          << " entries, " << m_noArbLines << " needed";
    m_noArbLines = 0;
    return;
  }

  for(MInt i = 0; i < m_noArbLines; i++) {
    m_noArbLineIds[i] = Context::getSolverProperty<MInt>("pp_noPointsArbLine", m_postprocessingId, AT_, i);
    if(m_noArbLineIds[i] < 2) {
      m_log << "property pp_noPointsArbLine[ " << i << " ] needs to be at least 2" << endl;
      m_noArbLineIds[i] = 2;
    }
  }


  mAlloc(m_arbLineIds, m_noArbLines, "m_arbLineIds", AT_);
  mAlloc(m_arbLineCoordinates, m_noArbLines, "m_arbLineCoordinates", AT_);
  mAlloc(m_arbLineOffsets, m_noArbLines, "m_arbLineOffsets", AT_);
  mAlloc(m_globalNoArbLineIds, m_noArbLines, "m_globalNoArbLineIds", AT_);

  // suppress valgrind error
  for(MInt i = 0; i < m_noArbLines; i++) {
    m_arbLineIds[i] = nullptr;
    m_arbLineCoordinates[i] = nullptr;
  }


  MFloat point[3]{};

  for(MInt lineId = 0; lineId < m_noArbLines; lineId++) { // loop over all lines

    mAlloc(m_arbLineIds[lineId], m_noArbLineIds[lineId], "m_arbLineIds", AT_);
    mAlloc(m_arbLineCoordinates[lineId], m_noArbLineIds[lineId] * nDim, "m_arbLineCoordinates", AT_);

    MInt count = 0;
    // find the cellIds for all specified positions on the line
    for(MInt pId = 0; pId < m_noArbLineIds[lineId]; pId++) {
      for(MInt dimId = 0; dimId < nDim; dimId++) {
        point[dimId] =
            m_arbLinePoints[2 * nDim * lineId + dimId]
            + (MFloat)pId / (m_noArbLineIds[lineId] - 1)
                  * (m_arbLinePoints[2 * nDim * lineId + dimId + nDim] - m_arbLinePoints[2 * nDim * lineId + dimId]);
      }

      const MInt cellId = pp->solver().grid().raw().findContainingLeafCell(point);

      MInt tmpCellId = std::numeric_limits<MInt>::max();

      if(cellId > -1) tmpCellId = m_gridProxy->tree().globalId(cellId);

      MPI_Allreduce(MPI_IN_PLACE, &tmpCellId, 1, MPI_INT, MPI_MIN, pp->solver().mpiComm(), AT_, "MPI_IN_PLACE",
                    "tmpCellId");

      if(tmpCellId == std::numeric_limits<MInt>::max())
        m_log << "no cell contains line " << lineId << ", point " << pId << endl;

      if(cellId == -1) {
        continue;
      }

      if(m_gridProxy->tree().globalId(cellId) > tmpCellId) {
        cout << " deleting double identified containing cell " << lineId << " " << pId << endl;
        continue;
      }

      m_arbLineIds[lineId][count] = cellId;
      for(MInt dimId = 0; dimId < nDim; dimId++) {
        m_arbLineCoordinates[lineId][nDim * count + dimId] = point[dimId];
      }
      count++;
    }

    //#cells found on the line
    m_noArbLineIds[lineId] = count;

    ParallelIo::size_type localCount, offset, totalCount;
    localCount = m_noArbLineIds[lineId];
    ParallelIo::calcOffset(localCount, &offset, &totalCount, pp->solver().mpiComm());
    m_arbLineOffsets[lineId] = offset;
  }

  MPI_Allreduce(m_noArbLineIds, m_globalNoArbLineIds, m_noArbLines, MPI_INT, MPI_SUM, pp->solver().mpiComm(), AT_,
                "m_noArbLineIds", "m_globalNoArbLineIds");
}

/** \fn void PostProcessing<nDim,ppType>::initProbeArbitrarySlice()
 * \brief initializes properties and data for arbitrary slice probing
 *
 * \author Ansgar Niemoeller, Sven Berger
 * \date 07.06.2014
 *
 * \param[in] grid pointer to the grid
 *
 **/

template <MInt nDim, class ppType>
void PostProcessing<nDim, ppType>::initProbeArbitrarySlice() {
  TRACE();


  /*! \page propertiesPP
    \section pp_arbSlicePoints
    <code>MFloat* PostProcessing::m_arbSlicePoints</code>\n\n
    default = <code>empty</code>\n
    This property defines the points which span arbitrary probe slices.
    Specify each point with an x, y and z component.\n
    A slice is defined by three points, the number of slices is determined by the given total number
    of points.\n Given points a, b and c for a slice (with vectors b-a and c-a nonparallel), the
    positions on the slice are determined by\n
    p_ij = a + 1/(i-1)*(b-a) + 1/(j-1)*(c-a)\n
    The parameters i and j (number of points in direction b-a and c-a) are determined by the
    property pp_noPointsArbSlice.\n
    <ul>
    <li><code>coordinates</code></li>
    </ul>
    Keywords: <i>POSTPROCESSING</i>
  */

  const MInt no_pp_arbSlicePoints = Context::propertyLength("pp_arbSlicePoints", m_postprocessingId);

  if(no_pp_arbSlicePoints == 0 || no_pp_arbSlicePoints % (3 * nDim) != 0) {
    m_log << "property pp_arbSlicePoints has " << no_pp_arbSlicePoints << " entries, multiple of " << 3 * nDim
          << " needed";
    TERMM(1, "ERROR: Invalid number of pp_arbSlicePoints");
  } else {
    mAlloc(m_arbSlicePoints, no_pp_arbSlicePoints, AT_, "m_arbSlicePoints");
    m_noArbSlices = no_pp_arbSlicePoints / (3 * nDim);
    for(MInt i = 0; i < no_pp_arbSlicePoints; i++) {
      m_arbSlicePoints[i] = Context::getSolverProperty<MFloat>("pp_arbSlicePoints", m_postprocessingId, AT_, i);
    }

    // TODO labels:PP,DOC check if valid points
  }

  mAlloc(m_noArbSlicePoints, 2 * m_noArbSlices, "m_noArbSlicePoints", AT_);

  /*! \page propertiesPP
    \section pp_noPointsArbSlice
    <code>MInt* PostProcessing::m_noArbSlicePoints</code>\n\n
    default = <code>empty</code>\n
    Number of points on the two directions of arbitrary probe slices, give a count for each
    direction of each slice. See property pp_arbSlicePoints for more information.
    <ul>
    <li><code>integer &ge; 2</code></li>
    </ul>
    Keywords: <i>POSTPROCESSING</i>
  */
  const MInt no_pp_noPointsArbSlice = Context::propertyLength("pp_noPointsArbSlice", m_postprocessingId);
  if(no_pp_noPointsArbSlice != 2 * m_noArbSlices) {
    m_log << "property pp_noPointsArbSlice has " << no_pp_noPointsArbSlice << " entries, " << m_noArbSlices
          << " needed";
    TERMM(1, "ERROR: Invalid number of pp_noPointsArbSlice");
  } else {
    for(MInt i = 0; i < 2 * m_noArbSlices; i++) {
      m_noArbSlicePoints[i] = Context::getSolverProperty<MInt>("pp_noPointsArbSlice", m_postprocessingId, AT_, i);
      if(m_noArbSlicePoints[i] < 2) {
        m_log << "property pp_noPointsArbSlice[ " << i << " ] needs to be at least 2" << endl;
        m_noArbSlicePoints[i] = 2;
      }
    }
  }

  /*! \page propertiesPP
    \section pp_arbSlicePointsDistribution
    <code>MFloat* PostProcessing::m_arbSlicePointsDistribution</code>\n\n
    Distribution of the points in instead of the typical equidistant distribution we normaly uses.\n
    Keywords: <i>GLOBAL, POSTPROCESSING</i>
  */

  MBool hasNewDistribution = false;
  MFloat stdVal = 0.0;
  m_log << "NYRAE DEBUG 2 Context::propertyLength is "
        << Context::propertyLength("pp_arbSlicePointsDistribution", m_postprocessingId) << endl;

  if(Context::propertyLength("pp_arbSlicePointsDistribution", m_postprocessingId) > 1) {
    mAlloc(m_arbSlicePointsDistribution, Context::propertyLength("pp_arbSlicePointsDistribution", m_postprocessingId),
           "m_arbSlicePointsDistribution", AT_);
    for(MInt i = 0; i < Context::propertyLength("pp_arbSlicePointsDistribution", m_postprocessingId); ++i) {
      m_arbSlicePointsDistribution[i] =
          Context::getSolverProperty<MFloat>("pp_arbSlicePointsDistribution", m_postprocessingId, AT_, &stdVal, i);
    }

    hasNewDistribution = true;
  }

  /*! \page propertiesPP
    \section pp_movePointsToGrid
    <code>MInt* PostProcessing::m_movePointsToGrid</code>\n\n
    This property is used to activate the moving of the source points for arbitrary lines or slices onto the grid.\n
    0: Points are not moved (default).\n
    1: Points outside of the grid are move to the nearest cell center.\n
    2: All points are moved to the nearest cell center.
    Keywords: <i>GLOBAL, POSTPROCESSING</i>
  */

  m_movePointsToGrid = 0;
  if(Context::propertyLength("pp_movePointsToGrid", m_postprocessingId) > 0) {
    m_movePointsToGrid = Context::getSolverProperty<MInt>("pp_movePointsToGrid", m_postprocessingId, AT_);
  }
  m_log << "NYRAE DEBUG 4 " << endl;

  /*! \page propertiesPP
    \section pp_spatialAveraging
    <code>MInt* PostProcessing::m_spatialAveraging</code>\n\n
    This property is used to activate the averaging of all variables in the slice/line probed.\n
    0: Deactivated (default).\n
    1: Activated.\n
    Keywords: <i>GLOBAL, POSTPROCESSING</i>
  */

  m_spatialAveraging = false;
  if(Context::propertyLength("pp_spartialAveraging", m_postprocessingId) > 0) {
    m_spatialAveraging = Context::getSolverProperty<MBool>("pp_spartialAveraging", m_postprocessingId, AT_);
  }

  collectMinLvlCells();

  if(0 < m_movePointsToGrid) {
    movePointsToGrid(m_arbSlicePoints, no_pp_arbSlicePoints, m_movePointsToGrid);
  }

  mAlloc(m_arbSliceIds, m_noArbSlices, "m_arbSliceIds", AT_);
  mAlloc(m_arbSliceCoordinates, m_noArbSlices, "m_arbSliceCoordinates", AT_);
  mAlloc(m_arbSliceOffsets, m_noArbSlices, "m_arbSliceOffsets", AT_);
  mAlloc(m_globalNoArbSlicePoints, m_noArbSlices, "m_globalNoArbSlicePoints", AT_);

  MFloat point[3];
  MInt cellId, count, noIds;

  MInt sliceIdOffset = 0;
  for(MInt sliceId = 0; sliceId < m_noArbSlices; sliceId++) { // loop over all slices

    noIds = m_noArbSlicePoints[2 * sliceId] * m_noArbSlicePoints[2 * sliceId + 1];

    mAlloc(m_arbSliceIds[sliceId], noIds, "m_arbSliceIds", AT_);
    mAlloc(m_arbSliceCoordinates[sliceId], noIds * nDim, "m_arbSliceCoordinates", AT_);

    count = 0;
    // find the cellIds for all specified positions on the slice

    for(MInt pId1 = 0; pId1 < m_noArbSlicePoints[2 * sliceId]; pId1++) {
      for(MInt pId2 = 0; pId2 < m_noArbSlicePoints[2 * sliceId + 1]; pId2++) {
        for(MInt dimId = 0; dimId < nDim; dimId++) {
          if(hasNewDistribution) {
            point[dimId] = m_arbSlicePoints[3 * nDim * sliceId + dimId]
                           + m_arbSlicePointsDistribution[pId1 + (sliceIdOffset)]
                                 * (m_arbSlicePoints[3 * nDim * sliceId + dimId + nDim]
                                    - m_arbSlicePoints[3 * nDim * sliceId + dimId])
                           + m_arbSlicePointsDistribution[pId2 + (sliceIdOffset + m_noArbSlicePoints[2 * sliceId])]
                                 * (m_arbSlicePoints[3 * nDim * sliceId + 2 * nDim + dimId]
                                    - m_arbSlicePoints[3 * nDim * sliceId + dimId]);
          } else {
            point[dimId] = m_arbSlicePoints[3 * nDim * sliceId + dimId]
                           + (MFloat)pId1 / (m_noArbSlicePoints[2 * sliceId] - 1)
                                 * (m_arbSlicePoints[3 * nDim * sliceId + dimId + nDim]
                                    - m_arbSlicePoints[3 * nDim * sliceId + dimId])
                           + (MFloat)pId2 / (m_noArbSlicePoints[2 * sliceId + 1] - 1)
                                 * (m_arbSlicePoints[3 * nDim * sliceId + 2 * nDim + dimId]
                                    - m_arbSlicePoints[3 * nDim * sliceId + dimId]);
          }
        }


        findContainingCell(point, cellId);
        if(cellId == -1) {
          continue;
        }
        m_arbSliceIds[sliceId][count] = cellId;
        for(MInt dimId = 0; dimId < nDim; dimId++) {
          m_arbSliceCoordinates[sliceId][nDim * count + dimId] = point[dimId];
        }
        count++;
      }
    }

    //#cells found on the slice
    m_noArbSlicePoints[sliceId] = count;

    ParallelIo::size_type localCount, offset, totalCount;
    localCount = m_noArbSlicePoints[sliceId];
    ParallelIo::calcOffset(localCount, &offset, &totalCount, pp->solver().mpiComm());
    m_arbSliceOffsets[sliceId] = offset;

    sliceIdOffset += m_noArbSlicePoints[2 * sliceId] + m_noArbSlicePoints[(2 * sliceId) + 1];
  }

  MPI_Allreduce(m_noArbSlicePoints, m_globalNoArbSlicePoints, m_noArbSlices, MPI_INT, MPI_SUM, pp->solver().mpiComm(),
                AT_, "m_noArbSlicePoints", "m_globalNoArbSlicePoints");
}

/**
 * \fn void PostProcessing<nDim,ppType>::initReduceToLevelAvg()
 * \brief Initializes properties for grid level reduction
 *
 * \author Andreas Lintermann
 * \date 14.09.2012
 *
 * \tparam[in] T celltype
 *
 **/
template <MInt nDim, class ppType>
void PostProcessing<nDim, ppType>::initReduceToLevelAvg() {
  TRACE();

  m_ReStressesAverageFileName = "";
  m_ReStressesAverageFileName = Context::getSolverProperty<MString>("pp_ReStressesAverageFileName", m_postprocessingId,
                                                                    AT_, &m_ReStressesAverageFileName);
  if(m_ReStressesAverageFileName.empty()) {
    mTerm(1, AT_, "Please specify the property 'ReStressesAverageFileName' ...");
  }
}

/**
 * \fn void PostProcessing<nDim,ppType>::initPeriodicSliceAverage()
 * \brief Initializes the periodic averaging on a slice
 *
 * \author Jannik Borgelt
 * \date 10.2023
 *
 * \tparam[in] T celltype
 *
 **/
template <MInt nDim, class ppType>
void PostProcessing<nDim, ppType>::initPeriodicSliceAverage() {
  TRACE();

  if(!postData().grid().isActive()) return;

  MInt noSlicePositions = 0;
  std::multimap<MFloat, MFloat> slicePositions;

  std::vector<MFloat> coord(2, F0);
  for(MInt cellId = 0; cellId < postData().c_noCells(); cellId++) {
    if(postData().c_noChildren(cellId) > 0) continue; // skip non leaf cells
    if(postData().a_isHalo(cellId)) continue;         // skip halo cells
    coord[0] = postData().c_coordinate(cellId, 0);
    coord[1] = postData().c_coordinate(cellId, 1);

    // see if coord is already present
    MInt count = slicePositions.count(coord[0]);

    if(count > 0) {
      MBool alreadyInserted = false;
      for(auto it = slicePositions.equal_range(coord[0]).first; it != slicePositions.equal_range(coord[0]).second;
          ++it) {
        if(approx((*it).second, coord[1], 1e-16)) {
          alreadyInserted = true;
        }
      }
      if(!alreadyInserted) {
        slicePositions.insert(std::make_pair(coord[0], coord[1]));
        noSlicePositions++;
      }
    } else {
      // insert
      slicePositions.insert(std::make_pair(coord[0], coord[1]));
      noSlicePositions++;
    }
  }

  // debug output slicePositions
  //----------------------------------------------------------------------------
  // stringstream fn;
  // fn.clear();
  // fn << "slicePositions_" << postData().domainId() << ".txt";
  // MString fname = fn.str();
  // ofstream slicePositionsOutput;
  // slicePositionsOutput.precision(16);
  // slicePositionsOutput.open(fname);
  // for(auto it = slicePositions.begin(); it != slicePositions.end(); ++it){
  //   MString line = "";
  //   line.append(to_string((*it).first) + " " + to_string((*it).second) + " " + to_string(postData().domainId()));
  //   slicePositionsOutput << line << endl;
  // }
  // slicePositionsOutput.close();
  //----------------------------------------------------------------------------

  // exchange multimap
  MFloatScratchSpace posX(noSlicePositions, AT_, "posX");
  MFloatScratchSpace posY(noSlicePositions, AT_, "posY");
  MInt count_ = 0;
  for(auto it = slicePositions.begin(); it != slicePositions.end(); ++it) {
    posX[count_] = (*it).first;
    posY[count_] = (*it).second;
    count_++;
  }

  MPI_Allreduce(&noSlicePositions, &m_globalnoSlicePositions, 1, MPI_INT, MPI_SUM, postData().mpiComm(), AT_,
                "noSlicePositions", "m_globalnoSlicePositions");

  MFloatScratchSpace globalPosX(m_globalnoSlicePositions, AT_, "globalPosX");
  MFloatScratchSpace globalPosY(m_globalnoSlicePositions, AT_, "globalPosY");

  ScratchSpace<MInt> recvbuf(postData().noDomains(), "recvbuf", FUN_);
  recvbuf.fill(0);
  MPI_Gather(&noSlicePositions, 1, MPI_INT, &recvbuf[0], 1, MPI_INT, 0, postData().mpiComm(), AT_, "noSlicePositions",
             "recvbuf");

  ScratchSpace<MInt> displs(postData().noDomains(), "displspos", FUN_);
  if(postData().domainId() == 0) {
    MInt offset = 0;
    for(MInt dom = 0; dom < postData().noDomains(); dom++) {
      displs[dom] = offset;
      offset += recvbuf[dom];
    }
  }

  MPI_Gatherv(&posX[0], noSlicePositions, MPI_DOUBLE, &globalPosX[0], &recvbuf[0], &displs[postData().domainId()],
              MPI_DOUBLE, 0, postData().mpiComm(), AT_, "posX", "globalPosX");
  MPI_Gatherv(&posY[0], noSlicePositions, MPI_DOUBLE, &globalPosY[0], &recvbuf[0], &displs[postData().domainId()],
              MPI_DOUBLE, 0, postData().mpiComm(), AT_, "posY", "globalCountZ");
  MPI_Bcast(&globalPosX[0], m_globalnoSlicePositions, MPI_DOUBLE, 0, postData().mpiComm(), AT_, "globalPosX");
  MPI_Bcast(&globalPosY[0], m_globalnoSlicePositions, MPI_DOUBLE, 0, postData().mpiComm(), AT_, "globalPosY");

  // find double entries
  for(MInt i = 0; i < m_globalnoSlicePositions; i++) {
    coord[0] = globalPosX[i];
    coord[1] = globalPosY[i];

    // see if coord is already present
    MInt count = m_sliceGlobalPositions.count(coord[0]);
    if(count > 0) {
      MBool alreadyInserted = false;
      for(auto it = m_sliceGlobalPositions.equal_range(coord[0]).first;
          it != m_sliceGlobalPositions.equal_range(coord[0]).second;
          ++it) {
        if(approx((*it).second, coord[1], 1e-16)) {
          alreadyInserted = true;
          break;
        }
      }
      if(!alreadyInserted) {
        m_sliceGlobalPositions.insert(std::make_pair(coord[0], coord[1]));
      }
    } else {
      // insert
      m_sliceGlobalPositions.insert(std::make_pair(coord[0], coord[1]));
    }
  }

  m_globalnoSlicePositions = m_sliceGlobalPositions.size();

  mAlloc(m_sliceAverage, postData().noVariables(), m_globalnoSlicePositions, "m_sliceAverage", AT_);
  for(MInt i = 0; i < m_globalnoSlicePositions; i++) {
    for(MInt varId = 0; varId < postData().noVariables(); varId++) {
      m_sliceAverage[varId][i] = F0;
    }
  }

  // determine index of cells in globalPos
  mAlloc(m_cell2globalIndex, m_globalnoSlicePositions, "m_cell2globalIndex", AT_);
  mAlloc(m_noPeriodicSliceCells, m_globalnoSlicePositions, "m_noPeriodicSliceCells", 0, AT_);
  for(MInt cellId = 0; cellId < postData().c_noCells(); cellId++) {
    if(postData().c_noChildren(cellId) > 0) continue; // skip non leaf cells
    if(postData().a_isHalo(cellId)) continue;         // skip halo cells
    MFloat x = postData().c_coordinate(cellId, 0);
    MFloat y = postData().c_coordinate(cellId, 1);
    MInt index = 0;
    for(auto it = m_sliceGlobalPositions.begin(); it != m_sliceGlobalPositions.end(); ++it) {
      if(approx((*it).first, x, 1e-16) && approx((*it).second, y, 1e-16)) {
        m_cell2globalIndex[index].push_back(cellId);
        m_noPeriodicSliceCells[index]++;
        break;
      }
      index++;
    }
  }

  MPI_Allreduce(MPI_IN_PLACE, &m_noPeriodicSliceCells[0], m_globalnoSlicePositions, MPI_INT, MPI_SUM,
                postData().mpiComm(), AT_, "MPI_IN_PLACE,", "m_noPeriodicSliceCells");

#ifndef NDEBUG
  // debug output slicePositions
  //----------------------------------------------------------------------------
  if(postData().domainId() == 0) {
    stringstream fng;
    fng.clear();
    fng << "sliceGlobalPositions.txt";
    MString fnameg = fng.str();
    ofstream sliceGlobalPositionsOutput;
    sliceGlobalPositionsOutput.precision(16);
    sliceGlobalPositionsOutput.open(fnameg);
    MInt index = 0;
    for(auto it = m_sliceGlobalPositions.begin(); it != m_sliceGlobalPositions.end(); ++it) {
      MString line = "";
      line.append(to_string((*it).first) + " " + to_string((*it).second) + " "
                  + to_string(m_noPeriodicSliceCells[index]));
      index++;
      sliceGlobalPositionsOutput << line << endl;
    }
    sliceGlobalPositionsOutput.close();
  }
  //----------------------------------------------------------------------------
#endif
}


template <MInt nDim, class ppType>
void PostProcessing<nDim, ppType>::periodicSliceAverage() {
  TRACE();

  if(!postData().grid().isActive()) return;

  postData().loadMeanFile(m_postprocessFileName);

  if(postData().domainId() == 0) cerr << "Calculating periodic average" << endl;

  for(MInt i = 0; i < m_globalnoSlicePositions; i++) {
    for(MInt c = 0; c < (MInt)m_cell2globalIndex[i].size(); c++) {
      MInt cellId = m_cell2globalIndex[i][c];
      for(MInt varId = 0; varId < postData().noVariables(); varId++) {
        m_sliceAverage[varId][i] += postData().m_averagedVars[cellId][varId] / m_noPeriodicSliceCells[i];
      }
    }
  }

  // exchange sliceAverage
  for(MInt varId = 0; varId < postData().noVariables(); varId++) {
    MPI_Allreduce(MPI_IN_PLACE, &m_sliceAverage[varId][0], m_globalnoSlicePositions, MPI_DOUBLE, MPI_SUM,
                  postData().mpiComm(), AT_, "MPI_IN_PLACE", "m_sliceAverage");
  }

  for(MInt i = 0; i < m_globalnoSlicePositions; i++) {
    for(MInt c = 0; c < (MInt)m_cell2globalIndex[i].size(); c++) {
      MInt cellId = m_cell2globalIndex[i][c];
      for(MInt varId = 0; varId < postData().noVariables(); varId++) {
        postData().a_variable(cellId, varId) = m_sliceAverage[varId][i];
      }
    }
  }

  postData().saveRestartFile(false);


#ifndef NDEBUG
  // debug output sliceAverage
  //----------------------------------------------------------------------------
  if(postData().domainId() == 0) {
    stringstream fng;
    fng.clear();
    fng << "sliceAverage_" << postData().domainId() << ".txt";
    MString fnameg = fng.str();
    ofstream sliceGlobalPositionsOutput;
    sliceGlobalPositionsOutput.precision(16);
    sliceGlobalPositionsOutput.open(fnameg);
    MInt i = 0;
    for(auto it = m_sliceGlobalPositions.begin(); it != m_sliceGlobalPositions.end(); ++it) {
      MString line = "";
      line.append(to_string((*it).first) + " " + to_string((*it).second) + " " + to_string(m_noPeriodicSliceCells[i])
                  + " " + to_string(m_sliceAverage[0][i]));
      i++;
      sliceGlobalPositionsOutput << line << endl;
    }
    sliceGlobalPositionsOutput.close();
  }
  //----------------------------------------------------------------------------
#endif
}

template <MInt nDim, class ppType>
void PostProcessing<nDim, ppType>::neededMeanVarsForSourceTerm(const MInt sourceTerm,
                                                               std::vector<MInt>& meanVars) const {
  TRACE();

  meanVars.clear();
  switch(sourceTerm) {
    case ST::Q_mI:
      // Mean Lamb vector
      meanVars.push_back(MV::LAMB0);
      break;
    case ST::Q_mI_linear:
      // Mean vorticities
      meanVars.push_back(MV::VORT0);
      break;
    case ST::Q_mII:
      // Mean gradient of rho
      meanVars.push_back(MV::DRHO);

      // Mean gradient of p
      meanVars.push_back(MV::DP);

      // Mean of (gradient of p divided by rho)
      meanVars.push_back(MV::GRADPRHO);
      break;
    case ST::Q_mIII:
      // Mean gradients of velocity components (contains MV::DU)
      meanVars.push_back(MV::GRADU);

      // Sum of products of velocity and velocity gradients:
      // u * grad(u) + v * grad(v) + w * grad(w)
      meanVars.push_back(MV::UGRADU);
      break;
    case ST::Q_e:
      // Components of divergence of u
      meanVars.push_back(MV::DU);

      // Mean gradient of rho
      meanVars.push_back(MV::DRHO);

      // Mean gradient of p
      meanVars.push_back(MV::DP);
      break;
    case ST::Q_c:
      // Components of divergence of u
      meanVars.push_back(MV::DU);

      // Mean gradient of rho
      meanVars.push_back(MV::DRHO);

      // Mean gradient of rho*div(u)
      meanVars.push_back(MV::RHODIVU);

      // Mean gradient of u*grad(rho)
      meanVars.push_back(MV::UGRADRHO);
      break;
    default:
      TERMM(1, "Source term '" + s_sourceTermNames[sourceTerm] + "' not implemented yet.");
      break;
  }
}


/** \brief Averages solutions.
 *
 * \author A. Niemoeller
 * \date 09.12.2013
 *
 * loads variables for all averaging timesteps and performs averaging (mean values and reynolds stress components) \n
 * additionally computation of skewness and kurtosis of the velocity components is performed if activated (see
 *properties "skewness"/"kurtosis") \n summation of variables can be done alternatively using kahan summation or two
 *pass in contrast to normal summation (see properties "useKahan"/"twoPass") \n solution output filename will be
 *"Mean_[averageStartTimestep]-[averageStopTimestep]"
 *
 **/
template <MInt nDim, class ppType>
void PostProcessing<nDim, ppType>::averageSolutions() {
  TRACE();

  if(!pp->solver().grid().isActive()) return;

  // Find out how many solutions are required and how the weighting will be
  const MFloat weight = 1.0 / (((m_averageStopTimestep - m_averageStartTimestep) / m_averageInterval) + 1);

  const MInt noCells = postData().grid().noInternalCells();

  const MFloat* cellVars = nullptr;
  MFloat* heatRelease = nullptr;

  if(m_twoPass) { // two pass activated, compute mean values first
    TERMM(1, "FIXME two-pass averaging is untested and probably broken");
    for(MInt t = m_averageStartTimestep; t <= m_averageStopTimestep; t += m_averageInterval) {
      m_log << "    ^        * Mean computation timestep " << t << " with weight " << weight << "\n";

      pp->solver().loadSampleVariables(t); // load variables for timestep t

      if(m_statisticCombustionAnalysis) pp->solver().getHeatRelease(heatRelease);

      for(MInt dataId = 0; dataId < noCells; dataId++) {
        MInt cellId = convertIdParent(postData(), pp->solver(), dataId);
        if(cellId != -1) {
          getSampleVariables(cellId, cellVars, true);

          if(cellVars == nullptr) {
            return;
          }

          for(MInt varId = 0; varId < m_noVariables; varId++) { // sum up all variables
            // m_summedVars[cellId][varId] += cellVars[varId];
            postData().a_variable(dataId, varId) += cellVars[varId];
          }
        }
      }
    }
    // weighting of summed primitive variables
    for(MInt dataId = 0; dataId < noCells; dataId++) {
      MInt cellId = convertIdParent(postData(), pp->solver(), dataId);
      if(cellId != -1) {
        for(MInt varId = 0; varId < m_noVariables; varId++) {
          // m_summedVars[cellId][varId] *= weight;
          postData().a_variable(dataId, varId) *= weight;
        }
      }
    }
  } // end of mean calculation for twoPass

  // Start averaging, timestep loop
  // PP_ToDo-Average:is this correct
  // MInt restartT = m_restartTimeStep;
  // if(m_averageRestart) {
  //   restartT += m_averageInterval;
  // } else {
  MInt restartT = m_averageStartTimestep;
  // }
  // TODO labels:PP move to separate method e.g. addAveragingSample() (see fv-structured PP)
  for(MInt t = restartT; t <= m_averageStopTimestep; t += m_averageInterval) {
    m_log << "    ^        * Averaging timestep " << t << " with weight " << weight << "\n";
    pp->solver().loadSampleVariables(t); // load variables for timestep t

    if(m_statisticCombustionAnalysis) pp->solver().getHeatRelease(heatRelease);

    for(MInt dataId = 0; dataId < noCells; dataId++) { // cell loop
      const MInt cellId = convertIdParent(postData(), pp->solver(), dataId);
      if(cellId != -1) {
        getSampleVariables(cellId, cellVars, true);

        if(cellVars == nullptr) {
          return;
        }

        if(m_twoPass) // two pass, second part
        {
          for(MInt varId = 0; varId < nDim; varId++) {
            // m_summedVars[cellId][varId + m_noVariables] +=
            //     (cellVars[varId] - m_summedVars[cellId][varId]) * (cellVars[varId] - m_summedVars[cellId][varId]);
            // postData().a_variable(dataId, varId) +=
            //   (cellVars[varId] - postData().a_variable(dataId, varId)) * (cellVars[varId] -
            //   m_summedVars[cellId][varId]);
          }
          for(MInt varId = 0; varId < 2 * nDim - 3; varId++) {
            // m_summedVars[cellId][m_noVariables + nDim + varId] +=
            //     (cellVars[varId % nDim] - m_summedVars[cellId][varId])
            //     * (cellVars[(varId + 1) % nDim] - m_summedVars[cellId][(varId + 1) % nDim]);
            // postData().a_variable(dataId, m_noVariables + nDim + varId) +=
            //   (cellVars[varId % nDim] - postData().a_variable(dataId, varId))
            //   * (cellVars[(varId + 1) % nDim] - m_summedVars[cellId][(varId + 1) % nDim]);
          }
          if(m_kurtosis) {
            for(MInt varId = 0; varId < nDim; varId++) {
              // m_summedVars[cellId][m_noVariables + 3 * (nDim - 1) + varId] +=
              //     pow(cellVars[varId] - m_summedVars[cellId][varId], 3);
              // m_summedVars[cellId][m_noVariables + 3 * (nDim - 1) + nDim + varId] +=
              //     pow(cellVars[varId] - m_summedVars[cellId][varId], 4);
              // postData().a_variable(dataId, m_noVariables + 3 * (nDim - 1) + varId) +=
              //   pow(cellVars[varId] - postData().a_variable(dataId,varId), 3);
              // postData().a_variable(dataId, m_noVariables + 3 * (nDim - 1) + nDim + varId) +=
              //   pow(cellVars[varId] - postData().a_variable(dataId, varId), 4);
            }
          } else if(m_skewness) {
            for(MInt varId = 0; varId < nDim; varId++) {
              // m_summedVars[cellId][m_noVariables + 3 * (nDim - 1) + varId] +=
              //     pow(cellVars[varId] - m_summedVars[cellId][varId], 3);
              // postData().a_variable(dataId, m_noVariables + 3 * (nDim - 1) + varId) +=
              //   pow(cellVars[varId] - postData().a_variable(dataId, varId), 3);
            }
          }
        }                   // two pass end
        else if(m_useKahan) // Kahan summation activated, reduced error in summation of variables
        {
          if(dataId == 0) m_log << "start kahan summation" << endl;
          /* Kahan summation pseudocode

             sum=0; c=0;
             for i=0 to input.length
             y = input[i] -c;
             t = sum + y;
             c = (t-sum) - y;
             sum = t;

          */
          for(MInt varId = 0; varId < m_noVariables; varId++) { // sum up all variables
            // m_ySum[cellId][varId] = cellVars[varId] - m_cSum[cellId][varId];
            // m_tSum[cellId][varId] = m_summedVars[cellId][varId] + m_ySum[cellId][varId];
            // m_cSum[cellId][varId] = (m_tSum[cellId][varId] - m_summedVars[cellId][varId]) - m_ySum[cellId][varId];
            // m_summedVars[cellId][varId] = m_tSum[cellId][varId];
          }
          for(MInt varId = 0; varId < nDim; varId++) { // squares of velocity components (u*u,v*v(,w*w))
            // m_ySquare[cellId][varId] = (cellVars[varId] * cellVars[varId]) - m_cSquare[cellId][varId];
            // m_tSquare[cellId][varId] = m_square[cellId][varId] + m_ySquare[cellId][varId];
            // m_cSquare[cellId][varId] = (m_tSquare[cellId][varId] - m_square[cellId][varId]) -
            // m_ySquare[cellId][varId]; m_square[cellId][varId] = m_tSquare[cellId][varId];
          }
          for(MInt varId = 0; varId < 2 * nDim - 3;
              varId++) { // products of different velocity components (u*v(,v*w,w*u))
            // m_ySquare[cellId][nDim + varId] =
            //(cellVars[varId]) * (cellVars[(varId + 1) % nDim]) - m_cSquare[cellId][nDim + varId];
            // m_tSquare[cellId][nDim + varId] = m_square[cellId][nDim + varId] + m_ySquare[cellId][nDim + varId];
            // m_cSquare[cellId][nDim + varId] =
            //(m_tSquare[cellId][nDim + varId] - m_square[cellId][nDim + varId]) - m_ySquare[cellId][nDim + varId];
            // m_square[cellId][nDim + varId] = m_tSquare[cellId][nDim + varId];
          }
          if(m_kurtosis) { // compute third and fourth power of velocity components for skewness and kurtosis
            for(MInt varId = 0; varId < nDim; varId++) {
              // m_yCube[cellId][varId] = pow(cellVars[varId], 3) - m_cCube[cellId][varId];
              // m_tCube[cellId][varId] = m_cube[cellId][varId] + m_yCube[cellId][varId];
              // m_cCube[cellId][varId] = (m_tCube[cellId][varId] - m_cube[cellId][varId]) - m_yCube[cellId][varId];
              // m_cube[cellId][varId] = m_tCube[cellId][varId];

              // m_yFourth[cellId][varId] = pow(cellVars[varId], 4) - m_cFourth[cellId][varId];
              // m_tFourth[cellId][varId] = m_fourth[cellId][varId] + m_yFourth[cellId][varId];
              // m_cFourth[cellId][varId] = (m_tFourth[cellId][varId] - m_fourth[cellId][varId]) -
              // m_yFourth[cellId][varId]; m_fourth[cellId][varId] = m_tFourth[cellId][varId];
            }
          } else if(m_skewness) { // compute only third power of velocity components for skewness
            for(MInt varId = 0; varId < nDim; varId++) {
              // m_yCube[cellId][varId] = pow(cellVars[varId], 3) - m_cCube[cellId][varId];
              // m_tCube[cellId][varId] = m_cube[cellId][varId] + m_yCube[cellId][varId];
              // m_cCube[cellId][varId] = (m_tCube[cellId][varId] - m_cube[cellId][varId]) - m_yCube[cellId][varId];
              // m_cube[cellId][varId] = m_tCube[cellId][varId];
            }
          }
        }    // kahan summation end
        else // normal summation
        {
          // Primitive variables
          MInt varIndex = postData().getPropertyVariableOffset("primitive").first;
          for(MInt varId = 0; varId < m_noVariables; varId++) {
            postData().a_variable(dataId, varIndex + varId) += cellVars[varId];
          }

          if(m_square) {
            MInt varSquare = postData().getPropertyVariableOffset("square").first;
            for(MInt varId = 0; varId < nDim; varId++) {
              // squares of velocity components (u*u,v*v(,w*w))
              postData().a_variable(dataId, varSquare + varId) += cellVars[varId] * cellVars[varId];
            }
            for(MInt varId = 0; varId < 2 * nDim - 3; varId++) {
              // products of different velocity components (u*v(,v*w,w*u))
              postData().a_variable(dataId, varSquare + nDim + varId) +=
                  (cellVars[varId % nDim]) * (cellVars[(varId + 1) % nDim]);
            }
            // squares of pressure  p*p
            postData().a_variable(dataId, varSquare + 3 * nDim - 3) += cellVars[nDim + 1] * cellVars[nDim + 1];
          }

          // third and fourth powers of velocity components (skewness and kurtosis)
          if(m_kurtosis) {
            ASSERT(m_square, "");
            ASSERT(m_skewness, "");
            MInt varSkew = postData().getPropertyVariableOffset("skewness").first;
            MInt varKurt = postData().getPropertyVariableOffset("kurtosis").first;
            for(MInt varId = 0; varId < nDim; varId++) {
              postData().a_variable(dataId, varSkew + varId) += pow(cellVars[varId], 3);
              postData().a_variable(dataId, varKurt + varId) += pow(cellVars[varId], 4);
            }
            // third powers of velocity components (skewness)
          } else if(m_skewness) {
            ASSERT(m_square, "");
            MInt varSkew = postData().getPropertyVariableOffset("skewness").first;
            for(MInt varId = 0; varId < nDim; varId++) {
              postData().a_variable(dataId, varSkew + varId) += pow(cellVars[varId], 3);
            }
          }

          // hm,c',h'
          if(m_statisticCombustionAnalysis) {
            MInt varComb = postData().getPropertyVariableOffset("statisticCombustionAnalysis").first;
            postData().a_variable(dataId, varComb) += heatRelease[cellId];
            postData().a_variable(dataId, varComb + 1) += cellVars[nDim + 2] * cellVars[nDim + 2];
            postData().a_variable(dataId, varComb + 2) += heatRelease[cellId] * heatRelease[cellId];
          }

          // Vorticities
          if(m_averageVorticity) {
            MInt varVort = postData().getPropertyVariableOffset("averageVorticity").first;
            for(MInt varId = 0; varId < m_noAveragedVorticities; varId++) {
              postData().a_variable(dataId, varVort + varId) += pp->vorticityAtCell(cellId, varId);
            }
          }
        }
      }
    } // normal summation end
  }   // cell loop end

  if(m_twoPass) {
    TERMM(1, "FIXME");
    // weighting for twoPass
    for(MInt dataId = 0; dataId < noCells; dataId++) {
      MInt cellId = convertIdParent(postData(), pp->solver(), dataId);
      if(cellId != -1) {
        for(MInt varId = 0; varId < 3 * nDim - 3 + nDim * (m_skewness + m_kurtosis); varId++) {
          // m_summedVars[cellId][m_noVariables + varId] *= weight;
        }
      }
    }
  } else {
    m_computeAndSaveMean = true;
    computeAndSaveMean();
  }

  m_log << "    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ \n" << endl;
}

/** \brief Averages solutions during solving.
 *
 * \author A. Niemoeller
 * \date 09.12.2013
 *
 * performs averaging during solver run (mean values and reynolds stress components) \n
 * additionally computation of skewness and kurtosis of the velocity components is performed if activated (see
 *properties "skewness"/"kurtosis") \n summation of variables can be done alternatively using kahan summation in
 *contrast to normal summation (see property "useKahan") \n solution output filename will be
 *"Mean_[averageStartTimestep]-[averageStopTimestep]"
 *
 **/
template <MInt nDim, class ppType>
void PostProcessing<nDim, ppType>::averageSolutionsInSolve() {
  TRACE();

  // Prepare to write restart file for all ranks
  if(globalTimeStep == m_averageStopTimestep) {
    postData().m_forceWriteRestart = true;
  }

  if(!pp->solver().grid().isActive()) return;

  const MInt noCells = postData().grid().noInternalCells();

  // Following bool is used to disable direct addressing of [nDim+1], which
  // corresponds to the pressure in the FV solver. For LB, pressure is not a
  // solution variable as it is solving for an isothermal equation of state.
  // This might be done more beautiful in the future? Or we decide that each
  // solver has to provide this value even if meaningless?
#if not defined(MAIA_DISABLE_LB)
  constexpr MBool isothermal = (std::is_same<ppType, PostProcessingLb<nDim>>::value);
  m_averageSpeedOfSound = (m_averageSpeedOfSound && !isothermal);
#else
  constexpr MBool isothermal = false;
#endif

  // Determine the value of "gamma" from solver if needed
  // This check is placed here to be called during the first loop iteration, as
  // getDimensionalizationParams() does not work (yet) while being in the c'tor
  // of the PP solver
  if(m_averageSpeedOfSound && m_gamma < 0.0) {
    vector<pair<MFloat, MString>> dimParams;
    pp->solver().getDimensionalizationParams(dimParams);

    m_gamma = -1.0;
    for(auto&& p : dimParams) {
      if(p.second == "gamma") {
        m_gamma = p.first;
        break;
      }
    }

    if(m_gamma < 0.0) {
      TERMM(1, "Dimensionalization parameter for 'gamma' not found but needed "
               "for averaging the speed of sound.");
    }
  }

  const MBool skipRestartStep = m_averageRestart && (globalTimeStep == m_restartTimeStep);

  if(globalTimeStep == m_averageStartTimestep) {
    // Reset everything at averaging start timestep
    std::fill_n(&postData().a_variable(0, 0), noCells * postData().noVariables(), 0.0);
  }

  // Lets do the averaging only if we are on the right timestep
  if(globalTimeStep >= m_averageStartTimestep && (globalTimeStep - m_averageStartTimestep) % m_averageInterval == 0
     && globalTimeStep <= m_averageStopTimestep && !skipRestartStep) {
    m_log << "    ^        * Averaging timestep " << globalTimeStep << "\n";

    MFloat* heatRelease = 0;

    if(m_statisticCombustionAnalysis) {
      pp->solver().calculateHeatRelease();
      pp->solver().getHeatRelease(heatRelease);
    }

    if(m_acousticAnalysis) {
      MFloatScratchSpace QeI(pp->solver().a_noCells(), AT_, "QeI");
      MFloatScratchSpace QeIII(pp->solver().a_noCells(), AT_, "QeIII");
      MFloatScratchSpace cSquared(pp->solver().a_noCells(), AT_, "cSquared");
      MFloatScratchSpace drhodt(pp->solver().a_noCells(), AT_, "drhodt");
      pp->computeAcousticSourceTermQe(QeI, QeIII, cSquared, drhodt);
    }

    // To trigger calculating derivative of primitive variables
    const MBool needDerivative = m_averageSpeedOfSound || m_needVorticity
                                 || m_activeMeanVars.find(MV::GRADU) != m_activeMeanVars.end()
                                 || m_activeMeanVars.find(MV::DU) != m_activeMeanVars.end()
                                 || m_activeMeanVars.find(MV::UGRADU) != m_activeMeanVars.end()
                                 || m_activeMeanVars.find(MV::DRHO) != m_activeMeanVars.end()
                                 || m_activeMeanVars.find(MV::DP) != m_activeMeanVars.end()
                                 || m_activeMeanVars.find(MV::RHODIVU) != m_activeMeanVars.end()
                                 || m_activeMeanVars.find(MV::UGRADRHO) != m_activeMeanVars.end()
                                 || m_activeMeanVars.find(MV::GRADPRHO) != m_activeMeanVars.end();

    // exchange for correlation
    if(m_correlation) {
      const MFloat* cellVarsCorr = nullptr;
      for(MInt correlationId = 0; correlationId < m_noCorrelationLines; correlationId++) {
        for(MInt i = 0; i < m_noCorrelationIds[correlationId]; i++) {
          MInt id = m_correlationIds[correlationId][i];
          MInt cellId = convertIdParent(postData(), pp->solver(), id);
          getSampleVariables(cellId, cellVarsCorr, true);
          MInt corrVarIndex = m_correlationVariableIndex[correlationId];
          m_correlationExchangeVar[correlationId][i] = cellVarsCorr[corrVarIndex];
        }

        MInt noDomain = postData().noDomains();
        ScratchSpace<MInt> recvbuf(noDomain, "recvbuf", AT_);
        recvbuf.fill(0);

        MPI_Gather(&m_noCorrelationIds[correlationId], 1, MPI_INT, &recvbuf[0], 1, MPI_INT, 0, postData().mpiComm(),
                   AT_, "m_noCorrelationIds[correlationId]", "postData().noDomain()");

        ScratchSpace<MInt> displs(noDomain, "displs", AT_);
        if(postData().domainId() == 0) {
          MInt offset = 0;
          for(MInt dom = 0; dom < noDomain; dom++) {
            displs[dom] = offset;
            offset += recvbuf[dom];
          }
        }

        MPI_Gatherv(&m_correlationExchangeVar[correlationId][0], m_noCorrelationIds[correlationId], MPI_DOUBLE,
                    &m_globalCorrelationExchangeVar[correlationId][0], &recvbuf[postData().domainId()],
                    &displs[postData().domainId()], MPI_DOUBLE, 0, postData().mpiComm(), AT_,
                    "m_correlationExchangeVar", "m_globalCorrelationExchangeVar");

        MPI_Bcast(m_globalCorrelationExchangeVar[correlationId], m_globalNoCorrelationIds[correlationId], MPI_DOUBLE, 0,
                  postData().mpiComm(), AT_, "m_globalCorrelationExchangeVar");
      }
    }

#if defined(MAIA_GCC_COMPILER) && GNU < 10
    // NOTE: This is a bug in GCC 9 we're using on HAWK. With this fix clang
    // (correctly) complains.
    maia::parallelFor(0, noCells, [&, isothermal](MInt dataId) {
#else
    maia::parallelFor(0, noCells, [&](MInt dataId) {
#endif
      const MInt cellId = convertIdParent(postData(), pp->solver(), dataId);
      if(cellId != -1) {
        std::vector<MFloat> cellVars(m_noVariables);
        getSampleVariables(cellId, cellVars);

        if(m_useKahan) { // kahan summation activated
          /* Kahan summation pseudocode

             sum=0; c=0;
             for i=0 to input.length
             y = input[i] -c;
             t = sum + y;
             c = (t-sum) - y;
             sum = t;

          */
          for(MInt varId = 0; varId < m_noVariables; varId++) { // sum up all variables
            // m_ySum[cellId][varId] = cellVars[varId] - m_cSum[cellId][varId];
            // m_tSum[cellId][varId] = m_summedVars[cellId][varId] + m_ySum[cellId][varId];
            // m_cSum[cellId][varId] = (m_tSum[cellId][varId] - m_summedVars[cellId][varId]) - m_ySum[cellId][varId];
            // m_summedVars[cellId][varId] = m_tSum[cellId][varId];
          }
          for(MInt varId = 0; varId < nDim; varId++) { // squares of velocity components (u*u,v*v(,w*w))
            // m_ySquare[cellId][varId] = (cellVars[varId] * cellVars[varId]) - m_cSquare[cellId][varId];
            // m_tSquare[cellId][varId] = m_square[cellId][varId] + m_ySquare[cellId][varId];
            // m_cSquare[cellId][varId] = (m_tSquare[cellId][varId] - m_square[cellId][varId]) -
            // m_ySquare[cellId][varId]; m_square[cellId][varId] = m_tSquare[cellId][varId];
          }
          for(MInt varId = 0; varId < 2 * nDim - 3;
              varId++) { // products of different velocity components (u*v(,v*w,w*u))
            // m_ySquare[cellId][nDim + varId] =
            //     (cellVars[varId]) * (cellVars[(varId + 1) % nDim]) - m_cSquare[cellId][nDim + varId];
            // m_tSquare[cellId][nDim + varId] = m_square[cellId][nDim + varId] + m_ySquare[cellId][nDim + varId];
            // m_cSquare[cellId][nDim + varId] =
            //     (m_tSquare[cellId][nDim + varId] - m_square[cellId][nDim + varId]) - m_ySquare[cellId][nDim + varId];
            // m_square[cellId][nDim + varId] = m_tSquare[cellId][nDim + varId];
          }
          if(m_kurtosis) { // compute third and fourth power of velocity components for skewness and kurtosis
            for(MInt varId = 0; varId < nDim; varId++) {
              // m_yCube[cellId][varId] = pow(cellVars[varId], 3) - m_cCube[cellId][varId];
              // m_tCube[cellId][varId] = m_cube[cellId][varId] + m_yCube[cellId][varId];
              // m_cCube[cellId][varId] = (m_tCube[cellId][varId] - m_cube[cellId][varId]) - m_yCube[cellId][varId];
              // m_cube[cellId][varId] = m_tCube[cellId][varId];

              // m_yFourth[cellId][varId] = pow(cellVars[varId], 4) - m_cFourth[cellId][varId];
              // m_tFourth[cellId][varId] = m_fourth[cellId][varId] + m_yFourth[cellId][varId];
              // m_cFourth[cellId][varId] = (m_tFourth[cellId][varId] - m_fourth[cellId][varId]) -
              // m_yFourth[cellId][varId]; m_fourth[cellId][varId] = m_tFourth[cellId][varId];
            }
          } else if(m_skewness) { // compute only third power of velocity components for skewness
            for(MInt varId = 0; varId < nDim; varId++) {
              // m_yCube[cellId][varId] = pow(cellVars[varId], 3) - m_cCube[cellId][varId];
              // m_tCube[cellId][varId] = m_cube[cellId][varId] + m_yCube[cellId][varId];
              // m_cCube[cellId][varId] = (m_tCube[cellId][varId] - m_cube[cellId][varId]) - m_yCube[cellId][varId];
              // m_cube[cellId][varId] = m_tCube[cellId][varId];
            }
          }
        }      // kahan summation end
        else { // normal summation

          // Primitive variables
          for(MInt varId = 0; varId < m_noVariables; varId++) {
            postData().a_variable(dataId, varId) += cellVars[varId];
          }

          if(m_square) {
            // squares of velocity components (u*u,v*v(,w*w))
            const MInt varIndex = postData().getPropertyVariableOffset("square").first;
            constexpr MInt noSquareDiag = nDim;
            for(MInt varId = 0; varId < noSquareDiag; varId++) {
              postData().a_variable(dataId, varIndex + varId) += cellVars[varId] * cellVars[varId];
            }
            // products of different velocity components (u*v(,v*w,w*u))
            constexpr MInt noSquareMixed = 2 * nDim - 3;
            for(MInt varId = 0; varId < noSquareMixed; varId++) {
              postData().a_variable(dataId, varIndex + nDim + varId) +=
                  (cellVars[varId % nDim]) * (cellVars[(varId + 1) % nDim]);
            }
            // squares of pressure  p*p
            if constexpr(isothermal) {
              postData().a_variable(dataId, varIndex + noSquareDiag + noSquareMixed) += cellVars[nDim] * cellVars[nDim];
            } else {
              postData().a_variable(dataId, varIndex + noSquareDiag + noSquareMixed) +=
                  cellVars[nDim + 1] * cellVars[nDim + 1];
            }
          }

          if(m_correlation) {
            MInt varIndex = postData().getPropertyVariableOffset("correlation").first;
            for(MInt correlationId = 0; correlationId < m_noCorrelationLines; correlationId++) {
              MInt index = m_correlationIndexMapping[correlationId][dataId];
              if(index != -1) {
                MInt corrVarIndex = m_correlationVariableIndex[correlationId];
                MFloat refCorrVar = m_globalCorrelationExchangeVar[correlationId][index];
                postData().a_variable(dataId, varIndex + correlationId) += refCorrVar * cellVars[corrVarIndex];
              }
            }
          }

          // third and fourth powers of velocity components (skewness and kurtosis)
          if(m_kurtosis) {
            ASSERT(m_square, "");
            ASSERT(m_skewness, "");
            const MInt varSkew = postData().getPropertyVariableOffset("skewness").first;
            const MInt varKurt = postData().getPropertyVariableOffset("kurtosis").first;
            ASSERT(varSkew + (nDim - 1) < postData().noVariables(), "");
            ASSERT(varKurt + (nDim - 1) < postData().noVariables(), "");
            for(MInt varId = 0; varId < nDim; varId++) {
              postData().a_variable(dataId, varSkew + varId) += pow(cellVars[varId], 3);
              postData().a_variable(dataId, varKurt + varId) += pow(cellVars[varId], 4);
            }
            // third powers of velocity components (skewness)
          } else if(m_skewness) {
            ASSERT(m_square, "");
            MInt varSkew_ = postData().getPropertyVariableOffset("skewness").first;
            ASSERT(varSkew_ + (nDim - 1) < postData().noVariables(), "");
            for(MInt varId = 0; varId < nDim; varId++) {
              postData().a_variable(dataId, varSkew_ + varId) += pow(cellVars[varId], 3);
            }
          }

          // hm,c',h'
          if(m_statisticCombustionAnalysis) {
            MInt varIndex = postData().getPropertyVariableOffset("statisticCombustionAnalysis").first;
            postData().a_variable(dataId, varIndex) += heatRelease[cellId];
            postData().a_variable(dataId, varIndex + 1) += cellVars[nDim + 2] * cellVars[nDim + 2];
            postData().a_variable(dataId, varIndex + 2) += heatRelease[cellId] * heatRelease[cellId];
          }

          if(needDerivative) {
            // Get pointer to derivatives of primitive variables
            // Instead of returning a pointer an array is filled for two reasons:
            //  1. solver can do some kind of calculation (certain gradient might not be calculated during solution
            //  step)
            //  2. solver may use SoA or AoS memory layout (flexibility)
            std::vector<MFloat> cellVarsDeriv(m_noVariables * nDim);
            if(pp->getSampleVarsDerivatives(cellId, cellVarsDeriv)) {
              // if(cellVarsDeriv != nullptr) {
              // const MFloatTensor deriv(const_cast<MFloat*>(cellVarsDeriv), m_noVariables, nDim);
              const MFloatTensor deriv(const_cast<MFloat*>(cellVarsDeriv.data()), m_noVariables, nDim);
              // Speed of sound and derivatives
              if(m_averageSpeedOfSound) {
                // Store convenience variables
                const MFloat rho = cellVars[nDim];
                const MFloat p = cellVars[nDim + 1];
                const MFloat pByRho = p / rho;

                // Speed of sound: c = sqrt(gamma * p / rho)
                MInt varIndex = postData().getPropertyVariableOffset("averageSpeedOfSound").first;
                postData().a_variable(dataId, varIndex) += sqrt(m_gamma * pByRho);

                // Derivatives of speed of sound
                const MFloat factor = 0.5 * sqrt(m_gamma / (rho * p));
                for(MInt dimId = 0; dimId < nDim; dimId++) {
                  postData().a_variable(dataId, varIndex + dimId + 1) +=
                      factor * (deriv(nDim + 1, dimId) - pByRho * deriv(nDim, dimId));
                }
              }

              if(m_needVorticity) {
                constexpr MInt noVorticities = nDim * 2 - 3;
                MFloat vort[nDim * 2 - 3];
                calcVorticity(deriv, vort);

                // Vorticities
                if(m_averageVorticity) {
                  MInt varIndex = postData().getPropertyVariableOffset("averageVorticity").first;
                  for(MInt varId = 0; varId < noVorticities; varId++) {
                    postData().a_variable(dataId, varIndex + varId) += vort[varId];
                  }
                }
                // Lamb vector (vorticity x velocity)
                if(m_activeMeanVars.find(MV::LAMB0) != m_activeMeanVars.end()) {
                  MInt varIndex = postData().getPropertyVariableOffset("lamb0").first;
                  ASSERT(varIndex > -1 && varIndex < postData().noVariables(), "");
                  if constexpr(nDim == 2) {
                    postData().a_variable(dataId, varIndex) -= vort[0] * cellVars[1];
                    postData().a_variable(dataId, varIndex + 1) += vort[0] * cellVars[0];
                  } else {
                    postData().a_variable(dataId, varIndex) += vort[1] * cellVars[2] - vort[2] * cellVars[1];
                    postData().a_variable(dataId, varIndex + 1) += vort[2] * cellVars[0] - vort[0] * cellVars[2];
                    postData().a_variable(dataId, varIndex + 2) += vort[0] * cellVars[1] - vort[1] * cellVars[0];
                  }
                }
              }

              // Velocity gradients
              if(m_activeMeanVars.find(MV::GRADU) != m_activeMeanVars.end()) {
                MInt varIndex = postData().getPropertyVariableOffset("gradu").first;
                ASSERT(varIndex > -1 && varIndex < postData().noVariables(), "");
                MInt indexGradU = 0;
                for(MInt veloId = 0; veloId < nDim; veloId++) {
                  for(MInt dimId = 0; dimId < nDim; dimId++) {
                    postData().a_variable(dataId, varIndex + indexGradU) += deriv(veloId, dimId);
                    indexGradU++;
                  }
                }
                // offset += nDim * nDim;
              } else if(m_activeMeanVars.find(MV::DU) != m_activeMeanVars.end()) {
                MInt varIndex = postData().getPropertyVariableOffset("du").first;
                ASSERT(varIndex > -1 && varIndex < postData().noVariables(), "");
                // du/dx, dv/dy, dw/dz for the divergence of the velocity field
                for(MInt dimId = 0; dimId < nDim; dimId++) {
                  postData().a_variable(dataId, varIndex + dimId) += deriv(dimId, dimId);
                }
              }

              // u * grad(u) + v * grad(v) + w * grad(w)
              if(m_activeMeanVars.find(MV::UGRADU) != m_activeMeanVars.end()) {
                MInt varIndex = postData().getPropertyVariableOffset("ugradu").first;
                ASSERT(varIndex > -1 && varIndex < postData().noVariables(), "");
                for(MInt dimId = 0; dimId < nDim; dimId++) {
                  for(MInt veloId = 0; veloId < nDim; veloId++) {
                    postData().a_variable(dataId, varIndex + dimId /*indexUGradU*/) +=
                        cellVars[veloId] * deriv(veloId, dimId);
                  }
                }
              }

              // grad(rho)
              if(m_activeMeanVars.find(MV::DRHO) != m_activeMeanVars.end()) {
                MInt varIndex = postData().getPropertyVariableOffset("drho").first;
                ASSERT(varIndex > -1 && varIndex < postData().noVariables(), "");
                for(MInt dimId = 0; dimId < nDim; dimId++) {
                  postData().a_variable(dataId, varIndex + dimId) += deriv(nDim, dimId);
                }
              }

              // grad(p)
              if constexpr(!isothermal) {
                if(m_activeMeanVars.find(MV::DP) != m_activeMeanVars.end()) {
                  MInt varIndex = postData().getPropertyVariableOffset("dp").first;
                  ASSERT(varIndex > -1 && varIndex < postData().noVariables(), "");
                  for(MInt dimId = 0; dimId < nDim; dimId++) {
                    postData().a_variable(dataId, varIndex + dimId) += deriv(nDim + 1, dimId);
                  }
                }
              }

              // rho*div(u)
              if(m_activeMeanVars.find(MV::RHODIVU) != m_activeMeanVars.end()) {
                MInt varIndex = postData().getPropertyVariableOffset("rhodivu").first;
                ASSERT(varIndex > -1 && varIndex < postData().noVariables(), "");
                const MFloat rho = cellVars[nDim];
                for(MInt dimId = 0; dimId < nDim; dimId++) {
                  postData().a_variable(dataId, varIndex + dimId) += rho * deriv(dimId, dimId);
                }
              }

              // u*grad(rho)
              if(m_activeMeanVars.find(MV::UGRADRHO) != m_activeMeanVars.end()) {
                MInt varIndex = postData().getPropertyVariableOffset("ugradrho").first;
                ASSERT(varIndex > -1 && varIndex < postData().noVariables(), "");
                for(MInt dimId = 0; dimId < nDim; dimId++) {
                  postData().a_variable(dataId, varIndex + dimId) += cellVars[dimId] * deriv(nDim, dimId);
                }
              }

              // grad(p)/rho
              if constexpr(!isothermal) {
                if(m_activeMeanVars.find(MV::GRADPRHO) != m_activeMeanVars.end()) {
                  MInt varIndex = postData().getPropertyVariableOffset("gradprho").first;
                  ASSERT(varIndex > -1 && varIndex < postData().noVariables(), "");
                  const MFloat rho = cellVars[nDim];
                  for(MInt dimId = 0; dimId < nDim; dimId++) {
                    postData().a_variable(dataId, varIndex + dimId) += (deriv(nDim + 1, dimId) / rho);
                  }
                }
              }
            }
          }
        }
      }
    });
  }

  // Write mean file if we are finished
  if(globalTimeStep == m_averageStopTimestep) {
    m_computeAndSaveMean = true;
  }
}

template <MInt nDim, class ppType>
template <MBool inSolution>
void PostProcessing<nDim, ppType>::computeAndSaveDivergence() {
  TRACE();

  if(!(postData().isActive())) return;

  const MInt timeStep = pp->solver().getCurrentTimeStep();
  if constexpr(inSolution) {
    const MInt solutionInterval = pp->solver().m_solutionInterval;
    const MInt solutionOffset = pp->solver().m_solutionOffset;
    if(!(timeStep % solutionInterval == 0 && timeStep >= solutionOffset)) return;
  }

  // Calculate for each cell the divergence
  constexpr MInt noVars = 1;
  const MInt noCellsPost = postData().grid().noInternalCells();
  MFloatScratchSpace divergence(noCellsPost, noVars, AT_, "divergence");
  divergence.fill(0.0);
  for(MInt cellIdPost = 0; cellIdPost < noCellsPost; cellIdPost++) {
    const MInt cellIdSolver = convertIdParent(postData(), pp->solver(), cellIdPost);
    if(cellIdSolver < 0) continue;
    divergence(cellIdSolver) = calcDivergence(cellIdSolver);
  }

  // Write output file
  stringstream ss;
  ss << pp->solver().outputDir() << "Divergence_s" << to_string(postData().solverId()) << "_" << setw(8) << setfill('0')
     << timeStep << ParallelIo::fileExt();
  const MString fileName = ss.str();

  m_log << "Saving divergence velocity field" << endl;

  std::vector<MString> variableNames = {"DIVERGENCE_U"};
  postData().saveDataFile(false, fileName, noVars, variableNames, divergence.data());

  ParallelIo parallelIo(fileName, maia::parallel_io::PIO_APPEND, postData().mpiComm());
  parallelIo.setAttribute(pp->solver().time(), "time");
}

template <MInt nDim, class ppType>
void PostProcessing<nDim, ppType>::computeAndSaveMean() {
  TRACE();

  if(!(postData().isActive() && m_computeAndSaveMean)) return;

  // Find out how many solutions are averaged for weighting the summed variables
  const MFloat weight = 1.0 / (((m_averageStopTimestep - m_averageStartTimestep) / m_averageInterval) + 1);
  m_log << "Computing averages with weight: " << setprecision(12) << weight << "\n";

  const MInt noCells = postData().grid().noInternalCells();
  const MInt noVars = postData().noVariables();
  MFloatScratchSpace averagedVars(noCells, noVars, AT_, "averagedVars");
  averagedVars.fill(0.0);

  const MInt varPrim = postData().getPropertyVariableOffset("primitive").first;
  const MInt varSquare = postData().getPropertyVariableOffset("square").first;
  const MInt varSkew = postData().getPropertyVariableOffset("skewness").first;
  const MInt varKurt = postData().getPropertyVariableOffset("kurtosis").first;
  const MInt varCorrelation = postData().getPropertyVariableOffset("correlation").first;
  ASSERT(varPrim < postData().noVariables(), "postData noVariables too low");
  ASSERT(!m_square || varSquare < postData().noVariables(), "postData noVariables too low");
  ASSERT(!m_kurtosis || varKurt < postData().noVariables(), "postData noVariables too low");
  ASSERT(!m_skewness || varSkew < postData().noVariables(), "postData noVariables too low");
  ASSERT(!m_correlation || varCorrelation < postData().noVariables(), "postData noVariables too low");

#if not defined(MAIA_DISABLE_LB)
  constexpr MBool isothermal = (std::is_same<ppType, PostProcessingLb<nDim>>::value);
  m_averageSpeedOfSound = (m_averageSpeedOfSound && !isothermal);
#else
  constexpr MBool isothermal = false;
#endif

  // exchange for correlation
  if(m_correlation) {
    for(MInt correlationId = 0; correlationId < m_noCorrelationLines; correlationId++) {
      for(MInt i = 0; i < m_noCorrelationIds[correlationId]; i++) {
        MInt dataId = m_correlationIds[correlationId][i];
        MInt corrVarIndex = m_correlationVariableIndex[correlationId];
        // sum_vref
        m_correlationExchangeVarMean[correlationId][i] = postData().a_variable(dataId, varPrim + corrVarIndex);
        // sum_vref_vref
        m_correlationExchangeVarRMS[correlationId][i] = postData().a_variable(dataId, varSquare + corrVarIndex);
      }

      MInt noDomain = postData().noDomains();
      ScratchSpace<MInt> recvbuf(noDomain, "recvbuf", AT_);
      recvbuf.fill(0);

      MPI_Gather(&m_noCorrelationIds[correlationId], 1, MPI_INT, &recvbuf[0], 1, MPI_INT, 0, postData().mpiComm(), AT_,
                 "m_noCorrelationIds[correlationId]", "postData().noDomain()");

      ScratchSpace<MInt> displs(noDomain, "displs", AT_);
      if(postData().domainId() == 0) {
        MInt offset = 0;
        for(MInt dom = 0; dom < noDomain; dom++) {
          displs[dom] = offset;
          offset += recvbuf[dom];
        }
      }

      MPI_Gatherv(&m_correlationExchangeVarMean[correlationId][0], m_noCorrelationIds[correlationId], MPI_DOUBLE,
                  &m_globalCorrelationExchangeVarMean[correlationId][0], &recvbuf[postData().domainId()],
                  &displs[postData().domainId()], MPI_DOUBLE, 0, postData().mpiComm(), AT_,
                  "m_correlationExchangeVarMean", "m_globalCorrelationExchangeVarMean");
      MPI_Gatherv(&m_correlationExchangeVarRMS[correlationId][0], m_noCorrelationIds[correlationId], MPI_DOUBLE,
                  &m_globalCorrelationExchangeVarRMS[correlationId][0], &recvbuf[postData().domainId()],
                  &displs[postData().domainId()], MPI_DOUBLE, 0, postData().mpiComm(), AT_,
                  "m_correlationExchangeVarRMS", "m_globalCorrelationExchangeVarRMS");

      MPI_Bcast(m_globalCorrelationExchangeVarMean[correlationId], m_globalNoCorrelationIds[correlationId], MPI_DOUBLE,
                0, postData().mpiComm(), AT_, "m_globalCorrelationExchangeVarMean");
      MPI_Bcast(m_globalCorrelationExchangeVarRMS[correlationId], m_globalNoCorrelationIds[correlationId], MPI_DOUBLE,
                0, postData().mpiComm(), AT_, "m_globalCorrelationExchangeVarRMS");
    }
  }


  for(MInt dataId = 0; dataId < noCells; dataId++) {
    MInt cellId = convertIdParent(postData(), pp->solver(), dataId);
    if(cellId != -1) {
      std::copy_n(&postData().a_variable(dataId, 0), noVars, &averagedVars(dataId, 0));

      // Weighting of all summed variables -> mean
      for(MInt varId = 0; varId < m_noVariables; varId++) {
        averagedVars(dataId, varPrim + varId) *= weight;
      }

      if(m_kurtosis) {
        // compute skewness and kurtosis of velocity components
        // e.g. skewness(u) = mean(u^3) - 3*u_mean*mean(u^2) + 2*u_mean^3
        // e.g. kurtosis(u) = mean(u^4) - 4*u_mean*mean(u^3) + 6*u_mean^2*mean(u^2) - 3*u_mean^4
        for(MInt varId = 0; varId < nDim; varId++) {
          averagedVars(dataId, varKurt + varId) =
              weight * averagedVars(dataId, varKurt + varId)
              - 4 * weight * averagedVars(dataId, varPrim + varId) * averagedVars(dataId, varSkew + varId)
              + 6 * weight * averagedVars(dataId, varSquare + varId) * pow(averagedVars(dataId, varPrim + varId), 2)
              - 3 * pow(averagedVars(dataId, varPrim + varId), 4);
          averagedVars(dataId, varSkew + varId) =
              weight * averagedVars(dataId, varSkew + varId)
              - 3 * weight * averagedVars(dataId, varPrim + varId) * averagedVars(dataId, varSquare + varId)
              + 2 * pow(averagedVars(dataId, varPrim + varId), 3);
        }

      } else if(m_skewness) {
        // compute skewness of velocity components
        for(MInt varId = 0; varId < nDim; varId++) {
          averagedVars(dataId, varSkew + varId) =
              weight * averagedVars(dataId, varSkew + varId)
              - 3 * weight * averagedVars(dataId, varPrim + varId) * averagedVars(dataId, varSquare + varId)
              + 2 * pow(averagedVars(dataId, varPrim + varId), 3);
        }
      }

      if(m_square) {
        // compute u',v'(,w') ( e.g. u'=mean(u^2)-(u_mean))^2 )
        for(MInt varId = 0; varId < nDim; varId++) {
          averagedVars(dataId, varSquare + varId) =
              weight * averagedVars(dataId, varSquare + varId) - pow(averagedVars(dataId, varPrim + varId), 2);
        }
        // compute u'v'(,v'w',w'u')  ( e.g. u'v'=mean(u*v)-u_mean*v_mean )
        for(MInt varId = 0; varId < 2 * nDim - 3; varId++) {
          averagedVars(dataId, varSquare + varId + nDim) =
              weight * averagedVars(dataId, varSquare + varId + nDim)
              - averagedVars(dataId, varPrim + varId % nDim) * averagedVars(dataId, varPrim + (varId + 1) % nDim);
        }
        // compute p'*p'
        if constexpr(isothermal) {
          averagedVars(dataId, varSquare + 3 * nDim - 3) =
              CSsq * CSsq * // So far this factor is only valid for LB: <p'p'> = CSsq^2 * <rho'rho'>
              (weight * averagedVars(dataId, varSquare + 3 * nDim - 3)
               - averagedVars(dataId, varPrim + nDim) * averagedVars(dataId, varPrim + nDim));
        } else {
          averagedVars(dataId, varSquare + 3 * nDim - 3) =
              weight * averagedVars(dataId, varSquare + 3 * nDim - 3)
              - averagedVars(dataId, varPrim + nDim + 1) * averagedVars(dataId, varPrim + nDim + 1);
        }
      }

      if(m_correlation) {
        MInt varIndex = postData().getPropertyVariableOffset("correlation").first;
        for(MInt correlationId = 0; correlationId < m_noCorrelationLines; correlationId++) {
          MInt index = m_correlationIndexMapping[correlationId][dataId];

          MFloat sum_varRef = m_globalCorrelationExchangeVarMean[correlationId][index];
          MFloat sum_varRefSquare = m_globalCorrelationExchangeVarRMS[correlationId][index];

          MInt corrVarIndex = m_correlationVariableIndex[correlationId];
          MFloat mean_varRef = weight * averagedVars(dataId, varIndex + correlationId)
                               - weight * sum_varRef * averagedVars(dataId, varPrim + corrVarIndex);

          MFloat mean_varRefSquare = weight * sum_varRefSquare - std::pow(weight * sum_varRef, 2);

          MFloat mean_varSquare = averagedVars(dataId, varSquare + corrVarIndex);

          averagedVars(dataId, varIndex + correlationId) = mean_varRef / std::sqrt(mean_varRefSquare * mean_varSquare);
        }
      }

      // compute hm,c',h'
      if(m_statisticCombustionAnalysis) {
        const MInt varComb = postData().getPropertyVariableOffset("statisticCombustionAnalysis").first;
        averagedVars(dataId, varComb) *= weight;
        averagedVars(dataId, varComb + 1) =
            averagedVars(dataId, varComb + 1) * weight - pow(averagedVars(dataId, varPrim + nDim + 2), 2);
        averagedVars(dataId, varComb + 2) =
            averagedVars(dataId, varComb + 2) * weight - pow(averagedVars(dataId, varComb), 2);
      }

      if(m_averageVorticity) {
        const MInt varVort = postData().getPropertyVariableOffset("averageVorticity").first;
        for(MInt varId = 0; varId < m_noAveragedVorticities; varId++) {
          averagedVars(dataId, varVort + varId) *= weight;
        }
      }

      if(m_averageSpeedOfSound) {
        const MInt varSos = postData().getPropertyVariableOffset("averageSpeedOfSound").first;
        for(MInt varId = 0; varId < m_noSpeedOfSoundVars; varId++) {
          averagedVars(dataId, varSos + varId) *= weight;
        }
      }

      if(m_noSourceVars > 0) {
        const MInt varSourceVars = postData().m_sourceVarsIndex;
        for(MInt varId = 0; varId < m_noSourceVars; varId++) {
          averagedVars(dataId, varSourceVars + varId) *= weight;
        }
      }
    }
  }

  stringstream ss;
  ss << pp->solver().outputDir() << "Mean_s" << to_string(postData().solverId()) << "_" << setw(8) << setfill('0')
     << m_averageStartTimestep << "-" << setw(8) << setfill('0') << m_averageStopTimestep << ParallelIo::fileExt();
  const MString name = ss.str();

  m_log << "Saving averaged variables " << name << endl;

  postData().saveDataFile(false, name, noVars, postData().m_variablesName, &averagedVars(0, 0));

  ParallelIo parallelIo(name, maia::parallel_io::PIO_APPEND, postData().mpiComm());
  parallelIo.setAttribute(1, "isMeanFile");
  parallelIo.setAttribute(m_averageStartTimestep, "averageStartTimestep");
  parallelIo.setAttribute(m_averageStopTimestep, "averageStopTimestep");
  parallelIo.setAttribute(m_averageInterval, "averageInterval");

  // reset the value
  m_computeAndSaveMean = false;
}

/** \brief moving average for all averaging timesteps
 *
 * \author Ansgar Niemoeller
 * \date 09.07.2014
 *
 * loads restart files and performs moving average
 *
 **/
template <MInt nDim, class ppType>
void PostProcessing<nDim, ppType>::movingAveragePost() {
  TRACE();

  MInt toNextOut;
  for(MInt t = m_movingAverageStartTimestep; t <= m_movingAverageStopTimestep; t += m_movingAverageInterval) {
    toNextOut = m_movingAverageInterval - (t - m_movingAverageStartTimestep) % m_movingAverageInterval;
    if(toNextOut < m_movingAverageDataPoints * m_movingAverageInterval) {
      pp->solver().loadSampleVariables(t);
      movingAverage();
    }
  }
}

/** \brief perform moving average
 *
 * \author Ansgar Niemoeller
 * \date 09.07.2014
 *
 * stores all required variables on timesteps which are considered for moving average
 * (see properties pp_movingAverageInterval, pp_movingAverageDataPoints)
 * and computes the moving average on averaging timesteps
 * (see properties pp_average[Start/Stop]Timestep, pp_averageInterval)\n
 * computes moving averages of primitive variables, these and the current variables
 * are written using ParallelIo\n
 * can additionaly be performed on the vorticity vector (see pp_averageVorticity)
 *
 **/
template <MInt nDim, class ppType>
void PostProcessing<nDim, ppType>::movingAverage() {
  TRACE();

  if(!pp->solver().grid().isActive()) return;

  MInt noCells = pp->solver().grid().noInternalCells();
  MInt noVars = m_noVariables;
  if(m_averageVorticity == 1) {
    noVars += (2 * nDim - 3);
  }
  const MFloat* vars = 0;

  // store data for averaging if we are on the right timestep
  if(globalTimeStep >= m_movingAverageStartTimestep - (m_movingAverageDataPoints - 1) * m_movingAverageInterval
     && (globalTimeStep - m_movingAverageStartTimestep) % m_movingAverageInterval == 0
     && globalTimeStep <= m_movingAverageStopTimestep) {
    MInt offset = (m_movingAverageCounter % m_movingAverageDataPoints) * noVars;

    for(MInt cellId = 0; cellId < noCells; cellId++) {
      // pp->solver().getSampleVariables(cellId, vars);
      getSampleVariables(cellId, vars, true);
      for(MInt varId = 0; varId < m_noVariables; varId++) {
        m_movAvgVariables[cellId][offset + varId] = vars[varId];
      }
    }

    if(m_averageVorticity) {
      ScratchSpace<MFloat> vorticity(pp->solver().a_noCells() * (2 * nDim - 3), AT_, "vorticity");
      pp->getVorticity(&vorticity[0]);
      for(MInt cellId = 0; cellId < noCells; cellId++) {
        for(MInt varId = m_noVariables; varId < noVars; varId++) {
          m_movAvgVariables[cellId][offset + varId] =
              vorticity[pp->solver().a_noCells() * (varId - m_noVariables) + cellId];
        }
      }
    }

    m_movingAverageCounter++;
  }

  // calculate moving average and write to file
  if(globalTimeStep >= m_movingAverageStartTimestep
     && (globalTimeStep - m_movingAverageStartTimestep) % m_movingAverageInterval == 0
     && globalTimeStep <= m_movingAverageStopTimestep) {
    m_log << "    ^        * Moving average timestep " << globalTimeStep << "\n";

    MInt noPoints = min(m_movingAverageCounter, m_movingAverageDataPoints);
    MInt offset = ((m_movingAverageCounter - 1) % m_movingAverageDataPoints) * noVars;
    ScratchSpace<MFloat> averagedVars(noCells * noVars, AT_, "averagedVars");
    ScratchSpace<MFloat> currentVars(noCells * noVars, AT_, "currentVars");

    for(MInt i = 0; i < noCells * noVars; i++) {
      averagedVars[i] = 0;
    }

    // collect all variables for current timeStep
    for(MInt cellId = 0; cellId < noCells; cellId++) {
      for(MInt varId = 0; varId < noVars; varId++) {
        currentVars[cellId * noVars + varId] = m_movAvgVariables[cellId][offset + varId];
      }
    }

    // sum up variables
    for(MInt cellId = 0; cellId < noCells; cellId++) {
      for(MInt pId = 0; pId < noPoints; pId++) {
        for(MInt varId = 0; varId < noVars; varId++) {
          averagedVars[cellId * noVars + varId] += m_movAvgVariables[cellId][pId * noVars + varId];
        }
      }
    }
    // normalize
    for(MInt i = 0; i < noCells * noVars; i++) {
      averagedVars[i] /= noPoints;
    }

    // output to parallelIo file
    ParallelIo::size_type cellsOffset, totalNoCells;
    MFloat currentTime = pp->solver().time();
    stringstream fileName;
    fileName << pp->solver().outputDir() << "movingAverage_" << globalTimeStep << ParallelIo::fileExt();
    ParallelIo parallelIo(fileName.str(), maia::parallel_io::PIO_REPLACE, pp->solver().mpiComm());
    parallelIo.calcOffset(noCells, &cellsOffset, &totalNoCells, pp->solver().mpiComm());

    // write global attributes
    parallelIo.setAttribute(pp->solver().grid().gridInputFileName(), "gridFile");
    parallelIo.setAttribute(globalTimeStep, "timeStep");
    parallelIo.setAttribute(currentTime, "time");
    parallelIo.setAttribute(m_movingAverageInterval, "movingAvgInterval");
    parallelIo.setAttribute(m_movingAverageDataPoints, "movingAvgDataPoints");

    // write parameters for dimensionalization
    vector<pair<MFloat, MString>> dimParams;
    pp->solver().getDimensionalizationParams(dimParams);

    for(MInt i = 0; i < (MInt)dimParams.size(); i++) {
      parallelIo.setAttribute(dimParams[i].first, dimParams[i].second);
    }

    // define all arrays in output file
    for(MInt varId = 0; varId < noVars; varId++) {
      stringstream varName;
      varName << "variables" << varId;
      parallelIo.defineArray(maia::parallel_io::PIO_FLOAT, varName.str(), totalNoCells);
      parallelIo.setAttribute(m_movAvgVarNames[varId], "name", varName.str());

      varName << "_mean";
      parallelIo.defineArray(maia::parallel_io::PIO_FLOAT, varName.str(), totalNoCells);
      parallelIo.setAttribute(m_movAvgVarNames[m_movAvgNoVariables + varId], "name", varName.str());
    }

    parallelIo.setOffset(noCells, cellsOffset);
    // write all variables
    for(MInt varId = 0; varId < noVars; varId++) {
      stringstream varName;
      varName << "variables" << varId;
      parallelIo.writeArray(&(currentVars[varId]), varName.str(), noVars);
      varName << "_mean";
      parallelIo.writeArray(&(averagedVars[varId]), varName.str(), noVars);
    }
  }
}

/** \brief Probes values at a location
 *
 * \author Andreas Lintermann
 * \date 22.08.2012
 *
 * \tparam[in] T celltype
 *
 **/
template <MInt nDim, class ppType>
void PostProcessing<nDim, ppType>::probeLocation() {
  TRACE();

  if(!pp->solver().grid().isActive()) return;

  const MFloat* cellVars = 0;
  for(MInt np = 0; np < m_noProbePoints; np++) {
    if(m_probeCellIds[np] == -1) continue;

    // pp->solver().getSampleVariables(m_probeCellIds[np], cellVars); // get variables for single cell
    getSampleVariables(m_probeCellIds[np], cellVars, true);

    // just to be sure, also write the domainId, in case the point lies exacly between 2 domains
    m_probeFileStreams[np] << pp->solver().domainId() << " " << globalTimeStep << " ";
    for(MInt i = 0; i < m_noVariables; i++) {
      m_probeFileStreams[np] << cellVars[i] << " ";
    }
    m_probeFileStreams[np] << endl;
  }
}

/** \brief line probing for all average timesteps or for a specified file
 *
 * \author Jannik Borgelt
 * \date 01.11.2020
 *
 * loads variables for all average timesteps or from specified file and performs periodic line probing
 *
 **/
template <MInt nDim, class ppType>
void PostProcessing<nDim, ppType>::probeLinePeriodicPost() {
  TRACE();

  if(!pp->solver().grid().isActive()) return;

  if(m_postprocessFileName != "") {
    m_log << "    ^        * probe line for file " << m_postprocessFileName << endl;
    TERMM(1, "FIXME untested");
    postData().loadMeanFile(m_postprocessFileName);
    probeLinePeriodic();
  } else {
    for(MInt t = m_averageStartTimestep; t <= m_averageStopTimestep; t += m_averageInterval) {
      pp->solver().loadSampleVariables(t);
      probeLinePeriodic();
    }
  }
}

/** \brief line probing for all average timesteps or for a specified file
 *
 * \author Jannik Borgelt
 * \date 14.03.2022
 *
 * loads variables for all average timesteps or from specified file and performs line probing
 *
 **/
template <MInt nDim, class ppType>
void PostProcessing<nDim, ppType>::probeLinePre() {
  TRACE();

  if(!pp->solver().grid().isActive()) return;

  for(MInt t = m_probeLineStartTimestep; t <= m_probeLineStopTimestep; t += m_probeLineInterval) {
    pp->solver().loadSampleVariables(t);
    probeLine();
  }
}

/** \brief line probing for all average timesteps or for a specified file
 *
 * \author A. Niemoeller
 * \date 09.12.2013
 *
 * loads variables for all average timesteps or from specified file and performs line probing
 *
 **/
template <MInt nDim, class ppType>
void PostProcessing<nDim, ppType>::probeLinePost() {
  TRACE();

  if(!pp->solver().grid().isActive()) return;

  m_log << "    ^        * probe line for file " << m_postprocessFileName << endl;
  postData().loadMeanFile(m_postprocessFileName);
  probeLine();
}

/** \brief Probes values on all probe lines
 *
 * \author A. Niemoeller
 * \date 09.12.2013
 *
 * writes single output file "probeLines_[globalTimeStep].Netcdf" containing all lines
 *
 **/
template <MInt nDim, class ppType>
void PostProcessing<nDim, ppType>::probeLine() {
  TRACE();

  if(!pp->solver().grid().isActive()) return;

  using namespace maia::parallel_io;

  // check for average timestep or postprocessFileName (globalTimeStep==0 if file is loaded with
  if((globalTimeStep >= m_probeLineStartTimestep
      && (globalTimeStep - m_probeLineStartTimestep) % m_probeLineInterval == 0
      && globalTimeStep <= m_probeLineStopTimestep)
     || globalTimeStep == 0) {
    MInt noVars = m_noVariables;
    vector<MString> datasetNames;
    datasetNames.clear();
    if(isMeanFile()) {
      noVars = postData().fileNoVars();
      datasetNames = postData().fileVarNames();
      TERMM_IF_NOT_COND((MInt)datasetNames.size() == noVars, "file var names size mismatch");
    } else {
      pp->solver().getSampleVariableNames(datasetNames);
    }

    const MInt step = (isMeanFile()) ? 0 : globalTimeStep;
    stringstream fileName;
    fileName << pp->solver().outputDir() << "probeLines_" << step << ParallelIo::fileExt();
    ParallelIo parallelIo(fileName.str(), maia::parallel_io::PIO_REPLACE, pp->solver().mpiComm());

    // define all arrays in output file
    // TODO labels:PP @ansgar_pp add full set of coordinates to file, or store line position as attributes
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
        if((MInt)datasetNames.size() == noVars) {
          parallelIo.setAttribute(datasetNames[varId], "name", varName.str());
        }
      }
    }

    for(MInt probeLineId = 0; probeLineId < m_noProbeLines; probeLineId++) { // loop over all probe lines
      if(m_probeLineDirection[probeLineId] < 0 || m_probeLineDirection[probeLineId] >= nDim) {
        continue;
      }

      m_log << "    ^        * probe line timestep " << globalTimeStep << " for line #" << probeLineId << endl;

      ScratchSpace<MFloat> vars(noVars * std::max(m_noProbeLineIds[probeLineId], 1), "vars", AT_);
      const MFloat* cellVars = nullptr;

      // collect local variables
      for(MInt i = 0; i < m_noProbeLineIds[probeLineId]; i++) {
        const MInt probeId = m_probeLineIds[probeLineId][i];
        getSampleVariables(probeId, cellVars, false);
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
  }
}

/** \brief arbitrary line probing for all average timesteps or for a specified file
 *
 * \author Ansgar Niemoeller
 * \date 07.06.2014
 *
 * loads variables for all average timesteps or from specified file and performs arbitrary line probing
 *
 **/
template <MInt nDim, class ppType>
void PostProcessing<nDim, ppType>::probeArbitraryLinePost() {
  TRACE();

  if(m_postprocessFileName != "") {
    m_log << "    ^        * probe arbitrary line for file " << m_postprocessFileName << endl;
    TERMM(1, "FIXME untested");
    postData().loadMeanFile(m_postprocessFileName);
    probeArbitraryLine();
  } else {
    for(MInt t = m_probeLineStartTimestep; t <= m_probeLineStopTimestep; t += m_probeLineInterval) {
      pp->solver().loadSampleVariables(t);
      probeArbitraryLine();
    }
  }
}

/** \brief Probes values on all arbitrary probe lines
 *
 * \author Ansgar Niemoeller
 * \date 07.06.2014
 *
 * writes single output file "probeArbitraryLines_[globalTimeStep].Netcdf" containing all lines
 *
 **/
template <MInt nDim, class ppType>
void PostProcessing<nDim, ppType>::probeArbitraryLine() {
  TRACE();

  using namespace maia::parallel_io;

  // check for average timestep or postprocessFileName (globalTimeStep==0 if file is loaded with
  if((globalTimeStep >= m_probeLineStartTimestep
      && (globalTimeStep - m_probeLineStartTimestep) % m_probeLineInterval == 0
      && globalTimeStep <= m_probeLineStopTimestep)
     || globalTimeStep == 0) {
    MInt noVars = m_noVariables;
    if(isMeanFile()) {
      noVars = postData().fileNoVars();
    }

    MInt step = (isMeanFile()) ? 0 : globalTimeStep;
    stringstream fileName;
    fileName << pp->solver().outputDir() << "probeArbitraryLines_" << step << ParallelIo::fileExt();
    ParallelIo parallelIo(fileName.str(), maia::parallel_io::PIO_REPLACE, pp->solver().mpiComm());

    // define all arrays in output file
    for(MInt probeLineId = 0; probeLineId < m_noArbLines; probeLineId++) {
      stringstream varNameBase;
      varNameBase << "line_" << probeLineId;

      for(MInt dimId = 0; dimId < nDim; dimId++) {
        stringstream coordName;
        coordName << varNameBase.str() << "_coordinate" << dimId;
        parallelIo.defineArray(PIO_FLOAT, coordName.str(), m_globalNoArbLineIds[probeLineId]);
        parallelIo.setAttribute(probeLineId, "lineId", coordName.str());
      }

      varNameBase << "_var_";
      for(MInt varId = 0; varId < noVars; varId++) {
        stringstream varName;
        varName << varNameBase.str() << varId;
        parallelIo.defineArray(PIO_FLOAT, varName.str(), m_globalNoArbLineIds[probeLineId]);
        parallelIo.setAttribute(probeLineId, "lineId", varName.str());
      }

      stringstream varName;
      varName << varNameBase.str() << noVars;
      parallelIo.defineArray(PIO_FLOAT, varName.str(), m_globalNoArbLineIds[probeLineId]);
      parallelIo.setAttribute(probeLineId, "lineId", varName.str());
    }

    parallelIo.defineScalar(PIO_FLOAT, "time");
    parallelIo.defineScalar(PIO_INT, "noLines");
    parallelIo.defineScalar(PIO_INT, "noVars");
    parallelIo.defineScalar(PIO_INT, "nDim");

    for(MInt probeLineId = 0; probeLineId < m_noArbLines; probeLineId++) { // loop over all probe lines
      m_log << "    ^        * probe arbitrary line timestep " << globalTimeStep << " for line #" << probeLineId
            << endl;

      MInt noIds = m_noArbLineIds[probeLineId];
      ScratchSpace<MFloat> vars((noVars + 1) * mMax(noIds, 1), "vars", AT_); // +1 for heat flux

      MInt probeId;
      MFloatScratchSpace pointVars(noVars, AT_, "pointVars");
      MFloat position[3];

      // collect local variables
      for(MInt i = 0; i < m_noArbLineIds[probeLineId]; i++) {
        probeId = m_arbLineIds[probeLineId][i];

        for(MInt j = 0; j < nDim; j++) {
          position[j] = m_arbLineCoordinates[probeLineId][nDim * i + j];
        }

        if(isMeanFile()) {
          // for(MInt varId = 0; varId < noVars; varId++) {
          //   pointVars[varId] = m_averagedVars[probeId][varId];
          // }
        } else {
          if(string2enum(pp->solver().solverType()) == MAIA_FINITE_VOLUME
             || string2enum(m_solverType) == MAIA_FV_GEQU_PV || string2enum(m_solverType) == MAIA_FV_LEVELSET
             || string2enum(m_solverType) == MAIA_FV_MB || string2enum(m_solverType) == MAIA_FV_MB_NEW_RK) {
            // reinterpret_cast<FvCartesianSolver<nDim>*>(pp->solver()).getPrimitiveVariables(probeId, position,
            // &pointVars[0], 1);
            pp->getPrimitiveVariables(probeId, position, &pointVars[0], 1);
          } else {
            TERMM(-1, "ERROR: Not implemented!");
          }
        }

        for(MInt varId = 0; varId < noVars; varId++) {
          vars[i * (noVars + 1) + varId] = pointVars[varId];
        }
        // only relevant for fv-solver types...
        if(string2enum(pp->solver().solverType()) == MAIA_FINITE_VOLUME || string2enum(m_solverType) == MAIA_FV_GEQU_PV
           || string2enum(m_solverType) == MAIA_FV_LEVELSET || string2enum(m_solverType) == MAIA_FV_MB
           || string2enum(m_solverType) == MAIA_FV_MB_NEW_RK) {
          vars[i * (noVars + 1) + noVars] =
              // reinterpret_cast<FvCartesianSolver<nDim>*>(pp->solver()).getBoundaryHeatFlux(probeId);
              pp->getBoundaryHeatFlux(probeId);
        }
      }

      parallelIo.setOffset(m_noArbLineIds[probeLineId], m_arbLineOffsets[probeLineId]);

      // write to file
      parallelIo.writeScalar(pp->solver().time(), "time");
      parallelIo.writeScalar(nDim, "nDim");
      parallelIo.writeScalar(m_noArbLines, "noLines");
      parallelIo.writeScalar(noVars + 1, "noVars");

      stringstream varNameBase;
      varNameBase << "line_" << probeLineId;
      for(MInt dimId = 0; dimId < nDim; dimId++) {
        stringstream coordName;
        coordName << varNameBase.str() << "_coordinate" << dimId;
        parallelIo.writeArray(&(m_arbLineCoordinates[probeLineId][dimId]), coordName.str(), nDim);
      }

      varNameBase << "_var_";
      for(MInt varId = 0; varId < noVars + 1; varId++) { // write all variables to file
        stringstream varName;
        varName << varNameBase.str() << varId;
        parallelIo.writeArray(&(vars[varId]), varName.str(), noVars + 1);
      }
    }
  }
}

/** \brief slice probing for all average timesteps or a specified file
 *
 * \author j. Borgelt
 * \date 09.03.2022
 *
 * loads variables for all average timesteps or from specified file and performs slice probing
 *
 **/
template <MInt nDim, class ppType>
void PostProcessing<nDim, ppType>::probeSlicePre() {
  TRACE();

  if(!pp->solver().grid().isActive()) return;

  // TODO labels:PP recombine probeSlicePre/Post
  for(MInt t = m_probeSliceStartTimestep; t <= m_probeSliceStopTimestep; t += m_probeSliceInterval) {
    m_log << "    ^        * probe slices for timestep " << t << endl;
    pp->solver().loadSampleVariables(t);
    probeSlice();
  }
}

/** \brief slice probing for all average timesteps or a specified file
 *
 * \author A. Niemoeller
 * \date 09.12.2013
 *
 * loads variables for all average timesteps or from specified file and performs slice probing
 *
 **/
template <MInt nDim, class ppType>
void PostProcessing<nDim, ppType>::probeSlicePost() {
  TRACE();

  if(!pp->solver().grid().isActive()) return;

  m_log << "    ^        * probe slice for file " << m_postprocessFileName << endl;
  postData().loadMeanFile(m_postprocessFileName);
  probeSlice();
}

template <MInt nDim, class ppType>
void PostProcessing<nDim, ppType>::probeSliceIn() {
  TRACE();

  if(!pp->solver().grid().isActive()) return;

  if(m_statisticCombustionAnalysis) {
    pp->solver().calculateHeatRelease();
  }

  probeSlice();
}

/** \brief Probes values on all probe slices
 *
 * \author A. Niemoeller
 * \date 09.12.2013
 *
 * writes single output file "probeSlices_[globalTimeStep].Netcdf" containing all slices
 * writes multiple files in the Aia File Format (Paraview) with m_sliceAiaFileFormat=1
 **/
template <MInt nDim, class ppType>
void PostProcessing<nDim, ppType>::probeSlice() {
  TRACE();

  if(!pp->solver().grid().isActive()) return;

  using namespace maia::parallel_io;

  IF_CONSTEXPR(nDim == 2) return;

  const MBool isSliceStep = (globalTimeStep >= m_probeSliceStartTimestep
                             && (globalTimeStep - m_probeSliceStartTimestep) % m_probeSliceInterval == 0
                             && globalTimeStep <= m_probeSliceStopTimestep)
                            || globalTimeStep == 0 || isMeanFile();

  cerr << "isSliceStep: " << pp->solver().domainId() << " " << isSliceStep << endl;
  cerr << pp->solver().domainId() << " " << m_probeSliceStartTimestep << " " << m_probeSliceStopTimestep << " "
       << m_probeSliceInterval << endl;
  cerr << (MBool)(globalTimeStep >= m_probeSliceStartTimestep) << " "
       << (MBool)((globalTimeStep - m_probeSliceStartTimestep) % m_probeSliceInterval == 0) << " "
       << (MBool)(globalTimeStep <= m_probeSliceStopTimestep) << " " << (MBool)(globalTimeStep == 0) << " "
       << (MBool)(isMeanFile()) << endl;

  // Implementation of CartesianGrid::createGridSlice()
  if(m_sliceAiaFileFormat) {
    if(isSliceStep) {
      MInt noVars = std::accumulate(m_noSliceVars.begin(), m_noSliceVars.end(), 0);
      MUint noVarIds = m_noSliceVars.size();

      if(isMeanFile()) {
        noVars = postData().fileNoVars();
        noVarIds = 1;
      } else {
        ASSERT(!m_sliceVarIds.empty(), "");
        pp->solver().calcSamplingVariables(m_sliceVarIds, false);
      }

      const MInt step = (isMeanFile()) ? 0 : globalTimeStep;

      for(MInt probeSliceId = 0; probeSliceId < m_noProbeSlices; probeSliceId++) {
        const MInt noIds = m_noProbeSliceIds[probeSliceId];
        MFloatScratchSpace vars(max(noIds, 1), noVars, AT_, "vars");

        // collect local variables
        for(MInt i = 0; i < m_noProbeSliceIds[probeSliceId]; i++) {
          const MInt probeId = m_probeSliceIds[probeSliceId][i];

          MInt varOffset = 0;
          for(MUint v = 0; v < noVarIds; v++) {
            ASSERT(!isMeanFile() || v < 1, "Error for probeSlice of mean file");
            const MInt varId = (isMeanFile()) ? -1 : m_sliceVarIds[v];
            calcSamplingVar(probeId, varId, &vars(i, varOffset));
            varOffset += m_noSliceVars[v];
          }
        }
        saveSliceAiaFileFormat(step, noVars, vars, probeSliceId);
      }
    }
  } else {
    // check for average timestep or postprocessFileName (globalTimeStep==0 if file is loaded with
    if(isSliceStep) {
      MInt noVars = m_noVariables;
      if(isMeanFile()) {
        noVars = postData().fileNoVars();
      }

      MInt step = (isMeanFile()) ? 0 : globalTimeStep;
      stringstream fileName;
      fileName << pp->solver().outputDir() << "probeSlices_" << step << ParallelIo::fileExt();
      ParallelIo parallelIo(fileName.str(), maia::parallel_io::PIO_REPLACE, pp->solver().mpiComm());

      // define all arrays in output file
      for(MInt probeSliceId = 0; probeSliceId < m_noProbeSlices; probeSliceId++) {
        stringstream varNameBase;
        varNameBase << "slice_" << probeSliceId;
        string coordName = varNameBase.str() + "_coordinates0";

        parallelIo.defineArray(PIO_FLOAT, coordName, m_globalNoProbeSliceIds[probeSliceId]);
        parallelIo.setAttribute(probeSliceId, "sliceId", coordName);

        coordName = varNameBase.str() + "_coordinates1";

        parallelIo.defineArray(PIO_FLOAT, coordName, m_globalNoProbeSliceIds[probeSliceId]);
        parallelIo.setAttribute(probeSliceId, "sliceId", coordName);

        varNameBase << "_var_";
        for(MInt varId = 0; varId < noVars; varId++) {
          stringstream varName;
          varName << varNameBase.str() << varId;
          parallelIo.defineArray(PIO_FLOAT, varName.str(), m_globalNoProbeSliceIds[probeSliceId]);
          parallelIo.setAttribute(probeSliceId, "sliceId", varName.str());
        }
      }

      for(MInt sliceId = 0; sliceId < m_noProbeSlices; sliceId++) { // loop over all slices
        if((m_probeSliceDir[2 * sliceId] < 0 || m_probeSliceDir[2 * sliceId] >= nDim)
           || (m_probeSliceDir[2 * sliceId + 1] < 0 || m_probeSliceDir[2 * sliceId + 1] >= nDim)
           || (m_probeSliceDir[2 * sliceId] == m_probeSliceDir[2 * sliceId + 1])) {
          return;
        }

        m_log << "    ^        * probe slice timestep " << globalTimeStep << " for probe slice #" << sliceId << endl;

        MInt noIds = m_noProbeSliceIds[sliceId];
        if(noIds == 0) noIds = 1; // avoid dereferencing array with length 0 in writeArray(...)
        ScratchSpace<MFloat> vars(noVars * noIds, "vars", AT_);

        const MFloat* cellVars = nullptr;

        // collect local variables
        for(MInt i = 0; i < m_noProbeSliceIds[sliceId]; i++) {
          const MInt probeId = m_probeSliceIds[sliceId][i];
          getSampleVariables(probeId, cellVars, false);
          for(MInt varId = 0; varId < noVars; varId++) {
            vars[i * noVars + varId] = cellVars[varId];
          }
        }

        parallelIo.setOffset(m_noProbeSliceIds[sliceId], m_probeSliceOffsets[sliceId]);

        // write to file
        stringstream varNameBase;
        varNameBase << "slice_" << sliceId;
        string coordName = varNameBase.str() + "_coordinates0";
        parallelIo.writeArray(&(m_probeSlicePositions[sliceId][0]), coordName, 2);
        coordName = varNameBase.str() + "_coordinates1";
        parallelIo.writeArray(&(m_probeSlicePositions[sliceId][1]), coordName, 2);

        varNameBase << "_var_";
        for(MInt varId = 0; varId < noVars; varId++) { // write all variables to file
          stringstream varName;
          varName << varNameBase.str() << varId;
          parallelIo.writeArray(&(vars[varId]), varName.str(), noVars);
        }
      }
    }
  }
}

/** \brief arbitrary slice probing for all average timesteps or for a specified file
 *
 * \author Ansgar Niemoeller
 * \date 07.06.2014
 *
 * loads variables for all average timesteps or from specified file and performs arbitrary slice probing
 *
 **/
template <MInt nDim, class ppType>
void PostProcessing<nDim, ppType>::probeArbitrarySlicePost() {
  TRACE();

  if(m_postprocessFileName != "") {
    m_log << "    ^        * probe arbitrary slice for file " << m_postprocessFileName << endl;
    TERMM(1, "FIXME untested");
    postData().loadMeanFile(m_postprocessFileName);
    probeArbitrarySlice();
  } else {
    for(MInt t = m_averageStartTimestep; t <= m_averageStopTimestep; t += m_averageInterval) {
      pp->solver().loadSampleVariables(t);
      probeArbitrarySlice();
    }
  }
}

/** \brief Probes values on all arbitrary probe slices
 *
 * \author Ansgar Niemoeller
 * \date 07.06.2014
 *
 * writes single output file "probeArbitrarySlices_[globalTimeStep].Netcdf" containing all slices
 *
 **/
template <MInt nDim, class ppType>
void PostProcessing<nDim, ppType>::probeArbitrarySlice() {
  TRACE();

  using namespace maia::parallel_io;

  // check for average timestep or postprocessFileName (globalTimeStep==0 if file is loaded with
  if((globalTimeStep >= m_averageStartTimestep && (globalTimeStep - m_averageStartTimestep) % m_averageInterval == 0
      && globalTimeStep <= m_averageStopTimestep)
     || globalTimeStep == 0) {
    MInt noVars = m_noVariables;
    if(isMeanFile()) {
      noVars = postData().fileNoVars();
    }

    const MInt step = (isMeanFile()) ? 0 : globalTimeStep;
    stringstream fileName;
    fileName << pp->solver().outputDir() << "probeArbitrarySlices_" << step << ParallelIo::fileExt();
    ParallelIo parallelIo(fileName.str(), maia::parallel_io::PIO_REPLACE, pp->solver().mpiComm());

    if(!m_spatialAveraging) {
      cout << globalDomainId() << ":Creating " << m_noArbSlices << " arbitary slices....";
      m_log << "Creating " << m_noArbSlices << " arbitary slices...";
    } else {
      cout << globalDomainId() << ":Creating " << m_noArbSlices << " arbitary slices with spatial averaging....";
      m_log << "Creating " << m_noArbSlices << " arbitary slices with spatial averaging....";
    }

    // define all arrays in output file
    for(MInt probeSliceId = 0; probeSliceId < m_noArbSlices; probeSliceId++) {
      stringstream varNameBase;
      varNameBase << "slice_" << probeSliceId;

      for(MInt dimId = 0; dimId < nDim; dimId++) {
        stringstream coordName;
        coordName << varNameBase.str() << "_coordinate" << dimId;
        parallelIo.defineArray(PIO_FLOAT, coordName.str(), m_globalNoArbSlicePoints[probeSliceId]);
        parallelIo.setAttribute(probeSliceId, "sliceId", coordName.str());
      }

      varNameBase << "_var_";
      for(MInt varId = 0; varId < noVars; varId++) {
        stringstream varName;
        varName << varNameBase.str() << varId;
        parallelIo.defineArray(PIO_FLOAT, varName.str(), m_globalNoArbSlicePoints[probeSliceId]);
        parallelIo.setAttribute(probeSliceId, "sliceId", varName.str());
      }

      //      if (m_spatialAveraging) {
      //        stringstream spatialAveragingName;
      //        spatialAveragingName << varNameBase.str() << "spatialAveraging";
      //        parallelIo.defineArray(PIO_FLOAT, spatialAveragingName.str(), noVars);
      //        parallelIo.setAttribute(probeSliceId, "sliceId", spatialAveragingName.str());
      //      }
    }

    for(MInt probeSliceId = 0; probeSliceId < m_noArbSlices; probeSliceId++) { // loop over all probe slices
      m_log << "    ^        * probe arbitrary slice timestep " << globalTimeStep << " for slice #" << probeSliceId
            << endl;

      MInt noPoints = m_noArbSlicePoints[probeSliceId];

      if(noPoints == 0) {
        noPoints = 1;
      } // avoid dereferencing array with length 0 in writeArray(...)
      MFloatScratchSpace vars(noVars * noPoints, AT_, "vars");
      MFloatScratchSpace pointVars(noVars, AT_, "pointVars");

      MInt probeId;
      // MFloatScratchSpace pointVars(noVars, AT_, "pointVars");
      MFloat position[3];

      // collect local variables
      for(MInt i = 0; i < m_noArbSlicePoints[probeSliceId]; i++) {
        probeId = m_arbSliceIds[probeSliceId][i];

        for(MInt j = 0; j < nDim; j++) {
          position[j] = m_arbSliceCoordinates[probeSliceId][nDim * i + j];
        }

        if(isMeanFile()) {
          TERMM(1, "FIXME");
          // for(MInt varId = 0; varId < noVars; varId++) {
          //   pointVars[varId] = m_averagedVars[probeId][varId];
          // }
        } else {
          pp->solver().getInterpolatedVariables(probeId, position, &pointVars[0]);
        }

        for(MInt varId = 0; varId < noVars; varId++) {
          vars[i * noVars + varId] = pointVars[varId];
        }
      }

      parallelIo.setOffset(m_noArbSlicePoints[probeSliceId], m_arbSliceOffsets[probeSliceId]);

      // write to file
      stringstream varNameBase;
      varNameBase << "slice_" << probeSliceId;
      for(MInt dimId = 0; dimId < nDim; dimId++) {
        stringstream coordName;
        coordName << varNameBase.str() << "_coordinate" << dimId;
        parallelIo.writeArray(&(m_arbSliceCoordinates[probeSliceId][dimId]), coordName.str(), nDim);
      }

      varNameBase << "_var_";
      for(MInt varId = 0; varId < noVars; varId++) { // write all variables to file
        stringstream varName;
        varName << varNameBase.str() << varId;
        parallelIo.writeArray(&(vars[varId]), varName.str(), noVars);
      }


      //      if (m_spatialAveraging) {
      //        MFloatScratchSpace spatialAveraging(noVars, AT_, "spatialAveraging");
      //        spatialAveraging.fill(0.0);
      //
      //        MInt globalNoIds = m_noArbSlicePoints[probeSliceId];
      //        if (1 < pp->solver().noDomains()) {
      //          MPI_Allreduce(MPI_IN_PLACE, &globalNoIds, 1, MPI_INT, MPI_SUM, pp->solver().mpiComm(), AT_,
      //          "MPI_IN_PLACE", "globalNoIds" );
      //        }
      //
      //        for (MInt varId = 0; varId < noVars; ++varId) {
      //          for (MInt i = 0; i < m_noArbSlicePoints[probeSliceId]; ++i) {
      //            spatialAveraging[varId] += vars[(i * noVars) + varId];
      //          }
      //
      //          if (0 < m_noArbSlicePoints[probeSliceId]) {
      //            spatialAveraging[varId] /= globalNoIds;
      //          }
      //
      //          if (1 < pp->solver().noDomains()) {
      //            MPI_Allreduce(MPI_IN_PLACE, &spatialAveraging[varId], 1, MPI_DOUBLE, MPI_SUM,
      //            pp->solver().mpiComm(), AT_, "MPI_IN_PLACE", "spatialAveraging[varId]" );
      //          }
      //        }
      //
      //
      //        if (0 == pp->solver().domainId()) {
      //          parallelIo.setOffset(noVars, 0);
      //        } else {
      //          parallelIo.setOffset(0, 0);
      //        }
      //        stringstream spatialAveragingName;
      //        spatialAveragingName << varNameBase.str() << "spatialAveraging";
      //        parallelIo.writeArray(&spatialAveraging[0], spatialAveragingName.str());
      //      }
    }
  }
}


/** \brief Averages slices.
 *
 * \author D. Zilles
 * \date 09.07.2017
 *
 * performs averaging before solver run (mean values) \n
 * currently only for slices created with pp_sliceAiaFileFormat=1 \n
 * reads slices stored in the output directory (probeSlice_[sliceAxis]_[sliceIntercept]_[timestep])" \n
 * solution output filename will be
 *"Mean_probeSlice_[sliceAxis]_[sliceIntercept]_[averageStartTimestep]-[averageStopTimestep]"
 *
 **/

template <MInt nDim, class ppType>
void PostProcessing<nDim, ppType>::averageSolutionsSlice() {
  TRACE();

  if(!pp->solver().grid().isActive()) return;

  using namespace maia::parallel_io;

  // loop through slices
  for(MInt sliceId = 0; sliceId < m_noProbeSlices; sliceId++) {
    stringstream initFileName;
    initFileName << pp->solver().outputDir() << "probeSlice_" << m_sliceAxis[sliceId] << "_"
                 << m_sliceIntercept[sliceId] << "_" << m_averageStartTimestep << ParallelIo::fileExt();

    ParallelIo initParallelIo(initFileName.str(), PIO_READ, pp->solver().mpiComm());

    // domain decomposition
    ParallelIo::size_type offset, totalCount;
    initParallelIo.readScalar(&totalCount, "noCells");
    MInt localCount = totalCount / pp->solver().noDomains();
    if(pp->solver().domainId() == pp->solver().noDomains() - 1) {
      localCount += (totalCount - pp->solver().noDomains() * localCount);
    }

    initParallelIo.calcOffset(localCount, &offset, &totalCount, pp->solver().mpiComm());

    MInt datasetSize = 0;
    vector<MString> datasetNames;

    // get dataset names
    vector<MString> variableNames = initParallelIo.getDatasetNames();

    for(MUint var = 0; var < variableNames.size(); var++) { // search for all datasets named variables0, variables1...

      if(variableNames[var].find("variables") != MString::npos) datasetSize++;
    }

    for(MInt var = 0; var < datasetSize; var++) {
      string tmpDatasetName;
      stringstream variableName;
      variableName << "variables" << var;
      initParallelIo.getAttribute(&tmpDatasetName, "name", variableName.str());
      datasetNames.push_back(tmpDatasetName);
    }

    // Initialize arrays
    MFloatScratchSpace tmpVar(localCount, AT_, "tmpVar");
    MFloatScratchSpace tmpSum(datasetSize, localCount, AT_, "tmpSum");
    MFloatScratchSpace tmpSquare(datasetSize, localCount, AT_, "tmpSquare");

    for(MInt var = 0; var < datasetSize; var++) {
      for(MInt cellId = 0; cellId < localCount; cellId++) {
        tmpSum(var, cellId) = 0;
        tmpSquare(var, cellId) = 0;
      }
    }

    // loop through timesteps
    for(MInt t = m_averageStartTimestep; t <= m_averageStopTimestep; t += m_averageInterval) {
      stringstream fileName;
      fileName << pp->solver().outputDir() << "probeSlice_" << m_sliceAxis[sliceId] << "_" << m_sliceIntercept[sliceId]
               << "_" << t << ParallelIo::fileExt();
      ParallelIo parallelIo(fileName.str(), PIO_READ, pp->solver().mpiComm());
      parallelIo.setOffset(localCount, offset);

      // loop through datasets
      for(MInt var = 0; var < datasetSize; var++) {
        stringstream variableName;
        variableName << "variables" << var;
        parallelIo.readArray(tmpVar.getPointer(), variableName.str());

        for(MInt cellId = 0; cellId < localCount; cellId++) {
          tmpSum[cellId * datasetSize + var] += tmpVar[cellId];
          tmpSquare[cellId * datasetSize + var] += tmpVar[cellId] * tmpVar[cellId];
        }
      }
    }
    MFloat weight = 1.0 / (((m_averageStopTimestep - m_averageStartTimestep) / m_averageInterval) + 1);
    for(MInt var = 0; var < datasetSize; var++) {
      for(MInt cellId = 0; cellId < localCount; cellId++) {
        tmpSum[cellId * datasetSize + var] = tmpSum[cellId * datasetSize + var] * weight;
        tmpSquare[cellId * datasetSize + var] =
            tmpSquare[cellId * datasetSize + var] * weight
            - tmpSum[cellId * datasetSize + var] * tmpSum[cellId * datasetSize + var];
      }
    }

    // save meanfile
    stringstream solutionFileName;
    solutionFileName << pp->solver().outputDir() << "Mean_probeSlice_" << m_sliceAxis[sliceId] << "_"
                     << m_sliceIntercept[sliceId] << "_" << m_averageStartTimestep << "-" << m_averageStopTimestep
                     << ParallelIo::fileExt();
    MString gridName;
    initParallelIo.getAttribute(&gridName, "gridFile", "");
    ParallelIo solutionParallelIo(solutionFileName.str(), PIO_REPLACE, pp->solver().mpiComm());
    solutionParallelIo.setAttribute(gridName, "gridFile");
    solutionParallelIo.defineScalar(PIO_INT, "noCells");

    for(MInt var = 0; var < datasetSize; var++) {
      stringstream varName;
      varName << "variables" << var;
      solutionParallelIo.defineArray(PIO_FLOAT, varName.str(), totalCount);
      solutionParallelIo.setAttribute(datasetNames[var] + "m", "name", varName.str());
    }
    for(MInt var = 0; var < datasetSize; var++) {
      stringstream varName;
      varName << "variables" << datasetSize + var;
      solutionParallelIo.defineArray(PIO_FLOAT, varName.str(), totalCount);
      solutionParallelIo.setAttribute(datasetNames[var] + "'", "name", varName.str());
    }

    solutionParallelIo.setOffset(localCount, offset);
    solutionParallelIo.writeScalar(totalCount, "noCells");

    for(MInt var = 0; var < datasetSize; var++) {
      stringstream varName;
      varName << "variables" << var;
      solutionParallelIo.writeArray(&(tmpSum[var]), varName.str(), datasetSize);
    }
    for(MInt var = 0; var < datasetSize; var++) {
      stringstream varName;
      varName << "variables" << datasetSize + var;
      solutionParallelIo.writeArray(&(tmpSquare[var]), varName.str(), datasetSize);
    }
  }
}

/** \brief access to variables for averaging
 *
 * \author A. Niemoeller
 * \date 23.04.14
 *
 * \param[in] cellId id of requested cell
 * \param[inout] vars pointer which is set to the start of the cell variables
 *
 **/
template <MInt nDim, class ppType>
void PostProcessing<nDim, ppType>::getSampleVariables(MInt cellId, const MFloat*& vars, MBool mode) {
  TRACE();

  if(mode) {
    pp->solver().getSampleVariables(cellId, vars);
  } else {
    if(isMeanFile()) {
      // MInt dataId = solver2DataIdParent(cellId);
      // vars = &postData().a_variable(dataId, 0);
      // TODO labels:PP is this correct, or use dataId?
      vars = &postData().a_averagedVariable(cellId, 0);
    } else {
      pp->solver().getSampleVariables(cellId, vars);
    }
  }
}
template <MInt nDim, class ppType>
void PostProcessing<nDim, ppType>::getSampleVariables(MInt const cellId, std::vector<MFloat>& vars) {
  TRACE();
  // TODO labels:PP Make this equivalent to getSampleVariables(MInt cellId, const MFloat*& vars, MBool mode) ?
  pp->solver().getSampleVariables(cellId, vars);
}

/// \brief Get the sampling variables associated with the given sampleVarId for the given cell
template <MInt nDim, class ppType>
void PostProcessing<nDim, ppType>::calcSamplingVar(const MInt cellId, const MInt sampleVarId, MFloat* const vars) {
  TRACE();

  if(isMeanFile()) {
    TERMM_IF_NOT_COND(sampleVarId == -1, "Error: sampleVarId for mean file needs to be -1.");
    MInt dataId = convertIdParent(pp->solver(), postData(), cellId);
    std::copy_n(&postData().a_variable(dataId, 0), postData().fileNoVars(), vars);
  } else {
    const MBool interpolate = false;
    pp->solver().calcSamplingVarAtPoint(nullptr, cellId, sampleVarId, vars, interpolate);
  }
}


/** \brief Spatial averaging for all average timesteps or for a specified file
 *
 * \author A. Niemoeller
 * \date 23.04.14
 *
 * loads variables for all averaging timesteps or from specified file and performes spatial averaging
 *
 **/
template <MInt nDim, class ppType>
void PostProcessing<nDim, ppType>::spatialAveragingPost() {
  TRACE();

  if(m_postprocessFileName != "") {
    m_log << "    ^        * spatial averaging for file " << m_postprocessFileName << endl;
    TERMM(1, "FIXME untested");
    postData().loadMeanFile(m_postprocessFileName);

    if(isMeanFile()) {
      if(pp->solver().domainId() == 0) {
        for(MInt i = 0; i < pp->solver().noDomains(); i++) {
          m_spatialVarsDispls[i] = (postData().fileNoVars() + 1) * m_spatialDispls[i];
          m_spatialVarsRecvcnts[i] = (postData().fileNoVars() + 1) * m_spatialRecvcnts[i];
        }
      }
    }

    spatialAveraging();
  } else {
    for(MInt t = m_averageStartTimestep; t <= m_averageStopTimestep; t += m_averageInterval) {
      pp->solver().loadSampleVariables(t);
      spatialAveraging();
    }
  }
}

// TODO labels:PP
/** \brief performs spatial averaging
 *
 * \author A. Niemoeller
 * \date 23.04.14
 *
 * depending on the set properties for the spatial directions, the flow field is averaged on a point/line/slice
 *
 **/
template <MInt nDim, class ppType>
void PostProcessing<nDim, ppType>::spatialAveraging() {
  TRACE();

  if((globalTimeStep >= m_averageStartTimestep && (globalTimeStep - m_averageStartTimestep) % m_averageInterval == 0
      && globalTimeStep <= m_averageStopTimestep)
     || globalTimeStep == 0) {
    m_log << "    ^        * spatial averaging timestep " << globalTimeStep << "\n";

    MInt noCells = pp->solver().grid().noInternalCells();
    MInt noVars = m_noVariables;
    if(isMeanFile()) {
      noVars = postData().fileNoVars();
    }
    MInt lvlCell;

    const MFloat* cellVars = 0;

    MInt dir1 = m_spatialDirection1;
    MInt dir2 = m_spatialDirection2;

    if(dir1 == -1 && dir2 == -1) { // single point

      ScratchSpace<MFloat> spatial_sum(noVars + 1, AT_, "spatial_sum");
      for(MInt i = 0; i < (noVars + 1); i++) {
        spatial_sum[i] = 0;
      }

      // average
      for(MInt cellId = 0; cellId < noCells; cellId++) {
        if(pp->solver().grid().tree().noChildren(cellId) > 0) continue; // check if leaf cell

        getSampleVariables(cellId, cellVars, false);

        lvlCell = pp->solver().grid().tree().level(cellId);

        for(MInt varId = 0; varId < noVars; varId++) {
          spatial_sum[varId] += 1 / m_spatialWeightSinglePoint * (m_spatialLvlWeight[lvlCell] * cellVars[varId]);
        }
      }
      spatial_sum[noVars] = m_spatialWeightSinglePoint;

      // exchange results
      MFloat globalWeight = m_spatialWeightSinglePoint;
      MPI_Allreduce(MPI_IN_PLACE, &globalWeight, 1, MPI_DOUBLE, MPI_SUM, pp->solver().mpiComm(), AT_, "MPI_IN_PLACE",
                    "globalWeight");

      for(MInt varId = 0; varId < noVars; varId++) {
        spatial_sum[varId] = m_spatialWeightSinglePoint / globalWeight * spatial_sum[varId];
      }

      MPI_Allreduce(MPI_IN_PLACE, spatial_sum.begin(), noVars, MPI_DOUBLE, MPI_SUM, pp->solver().mpiComm(), AT_,
                    "MPI_IN_PLACE", "spatial_sum.begin()"); // allreduce not nececcary

      spatial_sum[noVars] = globalWeight;
      //

      // TODO labels:PP,IO
      // write to file
      if(pp->solver().domainId() == 0) {
        stringstream fileName;
        fileName << "spatialAveraging_singlePoint.txt";
        ofstream file;
        file.open((fileName.str()).c_str());
        file.precision(12);
        for(MInt varId = 0; varId < noVars; varId++) {
          file << spatial_sum[varId] << "\t";
        }
        file << endl;
        file.close();
      }
    } else if((dir1 == -1 && dir2 != -1) || (dir2 == -1 && dir1 != -1)) {                        // line
      ScratchSpace<MFloat> spatial_sum(m_spatialLineNoCells * (noVars + 1), AT_, "spatial_sum"); // member?
      if(m_spatialLineAllVars == 0) {
        mAlloc(m_spatialLineAllVars, m_spatialCoordSum * (noVars + 1), "m_spatialLineAllVars", AT_);
      }

      for(MInt i = 0; i < m_spatialLineNoCells * (noVars + 1); i++) {
        spatial_sum[i] = 0;
      }

      MInt coordId, index;
      MFloat mapCoord, cellCoord, weight;
      for(MInt cellId = 0; cellId < noCells; cellId++) {
        if(pp->solver().grid().tree().noChildren(cellId) > 0) continue; // check if leaf cell

        getSampleVariables(cellId, cellVars, false);

        cellCoord = pp->solver().a_coordinate(cellId, dir1);
        map<MFloat, MFloat, coord_comp_1d_>::const_iterator it3 = m_spatialLineCellToMap.find(cellCoord);

        if(it3 != m_spatialLineCellToMap.end()) {
          mapCoord = (*it3).second;
        } else {
          mapCoord = cellCoord;
        }

        map<MFloat, MInt, coord_comp_1d_>::iterator it_coord = m_spatialLineCoordinates.find(mapCoord);
        coordId = distance(m_spatialLineCoordinates.begin(), it_coord);

        if(coordId >= m_spatialLineNoCells) {
          continue;
        }

        lvlCell = pp->solver().grid().tree().level(cellId);
        weight = spatial_sum[coordId * (noVars + 1) + noVars] + m_spatialLvlWeight[lvlCell]; // new summed weight

        for(MInt varId = 0; varId < noVars; varId++) {
          index = coordId * (noVars + 1) + varId;

          spatial_sum[index] = 1 / weight
                               * (m_spatialLvlWeight[lvlCell] * cellVars[varId]
                                  + spatial_sum[coordId * (noVars + 1) + noVars] * spatial_sum[index]);
        }
        spatial_sum[coordId * (noVars + 1) + noVars] = weight;
      }

      MPI_Gatherv(spatial_sum.begin(), m_spatialLineNoCells * (noVars + 1), MPI_DOUBLE, m_spatialLineAllVars,
                  m_spatialVarsRecvcnts, m_spatialVarsDispls, MPI_DOUBLE, 0, pp->solver().mpiComm(), AT_,
                  "spatial_sum.begin()", "m_spatialLineAllVars");

      ScratchSpace<MFloat> final_vars(m_spatialGlobalLineCoordinates.size() * (noVars + 1), AT_,
                                      "final_vars"); // make member?
      for(MUint varId = 0; varId < m_spatialGlobalLineCoordinates.size() * (noVars + 1); varId++) {
        final_vars[varId] = 0;
      }
      if(pp->solver().domainId() == 0) {
        for(MInt cellId = 0; cellId < m_spatialCoordSum; cellId++) { // assemble gathered variables
          mapCoord = m_spatialLineAllCoord[cellId];
          map<MFloat, MInt, coord_comp_1d_>::iterator it_coord = m_spatialGlobalLineCoordinates.find(mapCoord);

          if(it_coord == m_spatialGlobalLineCoordinates.end()) {
            it_coord = m_spatialGlobalLineCoordinates.find((*m_spatialGlobalLineCellToMap.find(mapCoord)).second);
          }
          if(it_coord == m_spatialGlobalLineCoordinates.end()) {
            continue;
          }

          coordId = distance(m_spatialGlobalLineCoordinates.begin(), it_coord);

          index = coordId * (noVars + 1);
          weight = m_spatialLineAllVars[cellId * (noVars + 1) + noVars];

          for(MInt varId = 0; varId < noVars; varId++) {
            final_vars[index + varId] = 1 / (final_vars[index + noVars] + weight)
                                        * (weight * m_spatialLineAllVars[cellId * (noVars + 1) + varId]
                                           + final_vars[index + noVars] * final_vars[index + varId]);
          }
          final_vars[index + noVars] += weight;
        }

        // test
        MFloat summed_weight = 0;
        for(MUint i = 0; i < m_spatialGlobalLineCoordinates.size(); i++) {
          summed_weight += final_vars[i * (noVars + 1) + noVars];
        }
        m_log << "spatial averaging: summed_weight " << summed_weight << endl;
        //
      }

      // write to file
      if(pp->solver().domainId() == 0) {
        stringstream fileName;
        fileName << "spatialAveraging_line_" << dir1 << "_" << globalTimeStep << ".txt";
        ofstream file;
        file.open((fileName.str()).c_str());
        file.precision(12);

        map<MFloat, MInt, coord_comp_1d_>::const_iterator it_ = m_spatialGlobalLineCoordinates.begin();
        for(MUint i = 0; i < m_spatialGlobalLineCoordinates.size(); i++) {
          file << (*it_).first << " ";
          for(MInt varId = 0; varId < noVars + 1; varId++) {
            file << final_vars[i * (noVars + 1) + varId] << " ";
          }
          it_++;
          file << endl;
        }
        file.close();
      }
    } else { // plane
      if(dir1 < 0 || dir1 > nDim - 1 || dir2 < 0 || dir2 > nDim - 1 || dir1 == dir2) {
        m_log << "invalid spatial directions " << dir1 << " " << dir2 << endl;
        return;
      }

      ScratchSpace<MFloat> spatial_sum(m_spatialPlaneNoCells * (noVars + 1), AT_, "spatial_sum"); // member?

      if(m_spatialPlaneAllVars == 0) {
        mAlloc(m_spatialPlaneAllVars, m_spatialCoordSum * (noVars + 1), "m_spatialPlaneAllVars", AT_);
      }

      for(MInt i = 0; i < m_spatialPlaneNoCells * (noVars + 1); i++) {
        spatial_sum[i] = 0;
      }

      MInt coordId, index;
      MFloat weight;
      pair<MFloat, MFloat> cellCoord, mapCoord;
      for(MInt cellId = 0; cellId < noCells; cellId++) {
        if(pp->solver().grid().tree().noChildren(cellId) > 0) continue; // check if leaf cell

        getSampleVariables(cellId, cellVars, false);

        cellCoord = make_pair(pp->solver().a_coordinate(cellId, dir1), pp->solver().a_coordinate(cellId, dir2));
        map<pair<MFloat, MFloat>, pair<MFloat, MFloat>, coord_comp_2d_>::const_iterator it3 =
            m_spatialPlaneCellToMap.find(cellCoord);

        if(it3 != m_spatialPlaneCellToMap.end()) {
          mapCoord = (*it3).second;
        } else {
          mapCoord = cellCoord;
        }

        map<pair<MFloat, MFloat>, MInt, coord_comp_2d_>::iterator it_coord = m_spatialPlaneCoordinates.find(mapCoord);
        coordId = distance(m_spatialPlaneCoordinates.begin(), it_coord);

        if(coordId >= m_spatialPlaneNoCells) {
          continue;
        }

        lvlCell = pp->solver().grid().tree().level(cellId);
        weight = spatial_sum[coordId * (noVars + 1) + noVars] + m_spatialLvlWeight[lvlCell]; // new summed weight

        for(MInt varId = 0; varId < noVars; varId++) {
          index = coordId * (noVars + 1) + varId;

          spatial_sum[index] = 1 / weight
                               * (m_spatialLvlWeight[lvlCell] * cellVars[varId]
                                  + spatial_sum[coordId * (noVars + 1) + noVars] * spatial_sum[index]);
        }
        spatial_sum[coordId * (noVars + 1) + noVars] = weight;
      }

      MPI_Gatherv(spatial_sum.begin(), m_spatialPlaneNoCells * (noVars + 1), MPI_DOUBLE, m_spatialPlaneAllVars,
                  m_spatialVarsRecvcnts, m_spatialVarsDispls, MPI_DOUBLE, 0, pp->solver().mpiComm(), AT_,
                  "spatial_sum.begin()", "m_spatialPlaneAllVars");

      ScratchSpace<MFloat> final_vars(m_spatialGlobalPlaneCoordinates.size() * (noVars + 1), AT_,
                                      "final_vars"); // make member?
      for(MUint varId = 0; varId < m_spatialGlobalPlaneCoordinates.size() * (noVars + 1); varId++) {
        final_vars[varId] = 0;
      }
      if(pp->solver().domainId() == 0) {
        for(MInt cellId = 0; cellId < m_spatialCoordSum; cellId++) { // assemble gathered variables
          mapCoord = make_pair(m_spatialPlaneAllCoord[2 * cellId], m_spatialPlaneAllCoord[2 * cellId + 1]);
          map<pair<MFloat, MFloat>, MInt, coord_comp_2d_>::iterator it_coord =
              m_spatialGlobalPlaneCoordinates.find(mapCoord);

          if(it_coord == m_spatialGlobalPlaneCoordinates.end()) {
            it_coord = m_spatialGlobalPlaneCoordinates.find((*m_spatialGlobalPlaneCellToMap.find(mapCoord)).second);
          }
          if(it_coord == m_spatialGlobalPlaneCoordinates.end()) {
            continue;
          }

          coordId = distance(m_spatialGlobalPlaneCoordinates.begin(), it_coord);

          index = coordId * (noVars + 1);
          weight = m_spatialPlaneAllVars[cellId * (noVars + 1) + noVars];

          for(MInt varId = 0; varId < noVars; varId++) {
            final_vars[index + varId] = 1 / (final_vars[index + noVars] + weight)
                                        * (weight * m_spatialPlaneAllVars[cellId * (noVars + 1) + varId]
                                           + final_vars[index + noVars] * final_vars[index + varId]);
          }
          final_vars[index + noVars] += weight;
        }

        // test
        MFloat summed_weight = 0;
        for(MUint i = 0; i < m_spatialGlobalPlaneCoordinates.size(); i++) {
          summed_weight += final_vars[i * (noVars + 1) + noVars];
        }
        m_log << "spatial averaging: summed_weight " << summed_weight << endl;
        //
      }

      // write to file
      if(pp->solver().domainId() == 0) {
        stringstream fileName;
        fileName << "spatialAveraging_plane_" << dir1 << "_" << dir2 << "_" << globalTimeStep << ".txt";
        ofstream file;
        file.open((fileName.str()).c_str());
        file.precision(12);

        /*
        file << "coordinate\t" << "um\t\t" << "vm\t\t";
        IF_CONSTEXPR(nDim == 3) {
          file << "wm\t\t";
        }
        file << "rhom\t\t" << "pm\t\t" << "summed weight" << endl;
        */

        auto it_ = m_spatialGlobalPlaneCoordinates.begin();
        for(MUint i = 0; i < m_spatialGlobalPlaneCoordinates.size(); i++) {
          file << (*it_).first.first << " " << (*it_).first.second << " ";
          for(MInt varId = 0; varId < noVars + 1; varId++) {
            file << final_vars[i * (noVars + 1) + varId] << " ";
          }
          it_++;
          file << endl;
        }
        file.close();
      }
    }
  }
}


// TODO labels:PP
/** \brief creates mapping for spatial averaging on a line
 *
 * \author A. Niemoeller
 * \date 23.04.14
 *
 * \param[inout] coordinates positions and levels of all cells projected on the line
 * \param[out] cell_to_map mapping from position to averaging position on the line
 * \param[in] grid pointer to the grid
 *
 **/
template <MInt nDim, class ppType>
void PostProcessing<nDim, ppType>::createCellToMap1D(map<MFloat, MInt, coord_comp_1d_>& coordinates,
                                                     map<MFloat, MFloat, coord_comp_1d_>& cell_to_map) {
  TRACE();

  MFloat length0 = m_gridProxy->cellLengthAtLevel(0);
  const MInt minLvl = m_gridProxy->minLevel();

  multimap<MFloat, MFloat>
      map_to_cell; // store small cells with key equal to the large cell coordinate in which to average
  pair<multimap<MFloat, MFloat>::iterator, multimap<MFloat, MFloat>::iterator> range;
  multimap<MFloat, MFloat>::iterator it_multi, it_multi2;
  MFloat coord;

  MInt c = 0;

  // find largest cells in the cut planes
  auto it = coordinates.begin();
  auto it2 = coordinates.begin();
  auto itTmp = coordinates.begin();
  auto itTmp2 = coordinates.begin();
  MFloat c1, c2;
  for(it = it2 = coordinates.begin(); it != coordinates.end();) {
    if(it2 == coordinates.end()) {
      it++;
      if(it == coordinates.end()) break;
      it2 = it;
      it2++;
      continue;
    } else if((*it).second == (*it2).second) { // same level
      it2++;
      continue;
    } else {
      itTmp = it;
      itTmp2 = it2;
      if((*it).second > (*it2).second) {
        it = it2;
        it2 = itTmp;
      }
      c1 = (*it).first;
      c2 = (*it2).first;
      if(fabs(c1 - c2) > length0 / FPOW2(minLvl)) {
        it = itTmp;
        it++;
        it2 = it;
        it2++;
        continue;
      } else if(((c1 < c2) && (c1 + 0.5 * length0 / FPOW2((*it).second) > c2))
                || ((c1 > c2) && (c1 - 0.5 * length0 / FPOW2((*it).second) < c2))) {
        range = map_to_cell.equal_range((*it2).first);
        c = distance(range.first, range.second);

        MInt r_count = 0;
        for(it_multi = range.first; c != 0 && r_count != c; r_count++) {
          it_multi2 = it_multi;
          it_multi2++; // access to next element before erase of it_multi
          coord = (*it_multi).second;
          map_to_cell.erase(it_multi);
          map_to_cell.insert(pair<MFloat, MFloat>((*it).first, coord));
          it_multi = it_multi2;
        }

        map_to_cell.insert(pair<MFloat, MFloat>((*it).first, (*it2).first));
        if(itTmp == it) {
          itTmp2++;
          coordinates.erase(it2);
          it2 = itTmp2;
        } else {
          itTmp++;
          coordinates.erase(it2);
          it = itTmp;
          it2 = itTmp2;
        }
      } else {
        it = itTmp;
        it2 = itTmp2;
        it2++;
      }
    }
  }

  // cell to map: cell/position -> average cell/position
  for(it_multi = map_to_cell.begin(); it_multi != map_to_cell.end(); it_multi++) {
    cell_to_map[(*it_multi).second] = (*it_multi).first;
  }
}

// TODO labels:PP
/** \brief creates mapping for spatial averaging on a slice
 *
 * \author A. Niemoeller
 * \date 23.04.14
 *
 * \param[inout] coordinates positions and levels of all cells projected on the slice
 * \param[out] cell_to_map mapping from position to averaging position on the slice
 * \param[in] grid pointer to the grid
 *
 **/
template <MInt nDim, class ppType>
void PostProcessing<nDim, ppType>::createCellToMap2D(
    map<pair<MFloat, MFloat>, MInt, coord_comp_2d_>& coordinates,
    map<pair<MFloat, MFloat>, pair<MFloat, MFloat>, coord_comp_2d_>& cell_to_map) {
  TRACE();

  MFloat length0 = m_gridProxy->cellLengthAtLevel(0);
  const MInt minLvl = m_gridProxy->minLevel();

  multimap<pair<MFloat, MFloat>, pair<MFloat, MFloat>>
      map_to_cell; // store small cells with key equal to the large cell coordinate in which to average
  pair<multimap<pair<MFloat, MFloat>, pair<MFloat, MFloat>>::iterator,
       multimap<pair<MFloat, MFloat>, pair<MFloat, MFloat>>::iterator>
      range;
  multimap<pair<MFloat, MFloat>, pair<MFloat, MFloat>>::iterator it_multi, it_multi2;
  pair<MFloat, MFloat> coord;

  MInt c = 0;

  // find largest cells in the cut planes
  auto it = coordinates.begin();
  auto it2 = coordinates.begin();
  auto itTmp = coordinates.begin();
  auto itTmp2 = coordinates.begin();
  pair<MFloat, MFloat> c1, c2;
  for(it = it2 = coordinates.begin(); it != coordinates.end();) {
    if(it2 == coordinates.end()) {
      it++;
      if(it == coordinates.end()) break;
      it2 = it;
      it2++;
      continue;
    } else if((*it).second == (*it2).second) { // same level
      it2++;
      continue;
    } else {
      itTmp = it;
      itTmp2 = it2;
      if((*it).second > (*it2).second) {
        it = it2;
        it2 = itTmp;
      }
      c1 = (*it).first;
      c2 = (*it2).first;
      if(fabs(c1.first - c2.first) > length0 / FPOW2(minLvl) && fabs(c1.second - c2.second) > length0 / FPOW2(minLvl)) {
        it2++;
        continue;
      }

      else if(((c1.first < c2.first) && (c1.first + 0.5 * length0 / FPOW2((*it).second) > c2.first))
              || ((c1.first > c2.first) && (c1.first - 0.5 * length0 / FPOW2((*it).second) < c2.first))) {
        if(((c1.second < c2.second) && (c1.second + 0.5 * length0 / FPOW2((*it).second) > c2.second))
           || ((c1.second > c2.second) && (c1.second - 0.5 * length0 / FPOW2((*it).second) < c2.second))) {
          range = map_to_cell.equal_range((*it2).first);
          c = distance(range.first, range.second);

          MInt r_count = 0;
          for(it_multi = range.first; c != 0 && r_count != c; r_count++) {
            it_multi2 = it_multi;
            it_multi2++; // access to next element before erase of it_multi
            coord = (*it_multi).second;
            map_to_cell.erase(it_multi);
            map_to_cell.insert(pair<pair<MFloat, MFloat>, pair<MFloat, MFloat>>((*it).first, coord));
            it_multi = it_multi2;
          }

          map_to_cell.insert(pair<pair<MFloat, MFloat>, pair<MFloat, MFloat>>((*it).first, (*it2).first));
          if(itTmp == it) {
            itTmp2++;
            coordinates.erase(it2);
            it2 = itTmp2;
          } else {
            itTmp++;
            coordinates.erase(it2);
            it = itTmp;
            it2 = itTmp2;
          }
        }
      } else {
        it = itTmp;
        it2 = itTmp2;
        it2++;
      }
    }
  }

  // cell to map: cell/position -> average cell/position
  for(it_multi = map_to_cell.begin(); it_multi != map_to_cell.end(); it_multi++) {
    cell_to_map[(*it_multi).second] = (*it_multi).first;
  }
}

/** \brief Finds the closest cell to a probe point
 *
 * \author Andreas Lintermann
 * \date 22.08.2012
 *
 * \tparam[in] T celltype
 *
 **/
template <MInt nDim, class ppType>
void PostProcessing<nDim, ppType>::findClosestProbePointsGridCells() {
  TRACE();

  m_log << "        = finding closest cells for all probe points on this domain" << endl;
  /*
    MFloatScratchSpace sqleveldiags(pp->solver().grid().maxLevel(), AT_, "leveldiags");
    for(MInt lev = 0; lev < pp->solver().grid().maxLevel(); lev++)
    {
    MFloat len = pp->solver().grid().cellLengthAtLevel(lev+1);
    sqleveldiags.p[lev] = 3*len*len;
    }
  */

  for(MInt np = 0; np < m_noProbePoints; np++) {
    // Initial length
    MFloat dist = 0.0;
    for(MInt d = 0; d < nDim; d++)
      dist += (pp->solver().a_coordinate(0, d) - m_probeCoordinates[np][d])
              * (pp->solver().a_coordinate(0, d) - m_probeCoordinates[np][d]);

    m_probeCellIds[np] = 0;

    for(MInt i = 0; i < pp->solver().grid().noInternalCells(); i++) {
      // skip cells without children
      if(m_gridProxy->tree().noChildren(i) != 0) continue;

      MFloat tmpdist = 0;
      for(MInt d = 0; d < nDim; d++)
        tmpdist += (pp->solver().a_coordinate(i, d) - m_probeCoordinates[np][d])
                   * (pp->solver().a_coordinate(i, d) - m_probeCoordinates[np][d]);

      if(tmpdist <= dist) {
        m_probeCellIds[np] = i;
        dist = tmpdist;
      }
    }

    MBool is_inside = true;
    for(MInt d = 0; d < nDim; d++)
      is_inside = is_inside
                  && (fabs(pp->solver().a_coordinate(m_probeCellIds[np], d) - m_probeCoordinates[np][d])
                      <= m_gridProxy->cellLengthAtLevel(m_gridProxy->tree().level(m_probeCellIds[np]) + 1));

    if(!is_inside) m_probeCellIds[np] = -1;
  }
}

template <MInt nDim, class ppType>
void PostProcessing<nDim, ppType>::saveSliceAiaFileFormat(const MInt step, const MInt noVars, MFloatScratchSpace& vars,
                                                          const MInt sliceId) {
  TRACE();

  stringstream fileName;
  vector<MString> datasetNames;

  if(!isMeanFile()) {
    // TODO labels:PP @ansgar use a unique name for each solver
    fileName << pp->solver().outputDir() << "probeSlice_" << m_sliceAxis[sliceId] << "_" << m_sliceIntercept[sliceId]
             << "_" << step << ParallelIo::fileExt();

    // Assemble all slice variable names
    datasetNames.clear();
    for(MUint v = 0; v < m_sliceVarIds.size(); v++) {
      datasetNames.insert(datasetNames.end(), m_sliceVarNames[v].begin(), m_sliceVarNames[v].end());
    }
  } else {
    const MUint slashPos = m_postprocessFileName.find_last_of("/");
    const MUint pos = (slashPos < m_postprocessFileName.length()) ? slashPos + 1 : 0;
    const MString fileNameRaw = m_postprocessFileName.substr(pos, m_postprocessFileName.find_last_of(".") - pos);
    // TODO labels:PP @ansgar use a unique name for each solver
    fileName << pp->solver().outputDir() << fileNameRaw << "_slice_" << m_sliceAxis[sliceId] << "_"
             << m_sliceIntercept[sliceId] << ParallelIo::fileExt();

    datasetNames = postData().fileVarNames();

    // if(m_fileVarNames.size() > 0) {
    //   datasetNames = m_fileVarNames;
    // } else {
    //   datasetNames = {"um", "vm", "wm",   "rhom", "pm",   "cm", "hm", "u'",
    //                   "v'", "w'", "u'v'", "v'w'", "w'u'", "p'", "c'", "h'"};
    // }
  }


  ParallelIo parallelIo(fileName.str(), maia::parallel_io::PIO_REPLACE, pp->solver().mpiComm());
  const MString varNameBase = "variables";

  // Write attributes
  parallelIo.setAttribute(m_globalNoProbeSliceIds[sliceId], "noCells");
  parallelIo.setAttribute(pp->solver().solverId(), "solverId");
  parallelIo.setAttribute(pp->solver().getCurrentTimeStep(), "globalTimeStep");
  parallelIo.setAttribute(pp->solver().time(), "time");
  parallelIo.setAttribute(m_probeSliceGridNames[sliceId], "gridFile");

  for(MInt varId = 0; varId < noVars; varId++) {
    stringstream varName;
    varName << varNameBase << varId;
    parallelIo.defineArray(maia::parallel_io::PIO_FLOAT, varName.str(), m_globalNoProbeSliceIds[sliceId]);
    stringstream name;
    if(datasetNames.size() != (MUint)noVars) {
      name << "var" << varId;
    } else {
      name << datasetNames[varId];
    }
    parallelIo.setAttribute(name.str(), "name", varName.str());
  }


  // Set variables for writing depening on the selected IO method
  MInt maxNoIds = -1;
  MInt noIds = -1;
  MInt* hilbertInfo = nullptr;
  if(m_optimizedSliceIo) {
    maxNoIds = m_noProbeSliceMaxNoContHilbertIds[sliceId];
    noIds = m_noProbeSliceNoContHilbertIds[sliceId];
    hilbertInfo = m_noProbeSliceContHilbertInfo[sliceId];
  } else {
    maxNoIds = m_noProbeSliceMaxNoHilbertIds[sliceId];
    noIds = m_noProbeSliceNoHilbertIds[sliceId];
    hilbertInfo = m_noProbeSliceHilbertInfo[sliceId];
  }

  const MFloat writeTimeStart = wallTime();
  for(MInt varId = 0; varId < noVars; varId++) { // write all variables to file
    stringstream varName;
    varName << varNameBase << varId;

    // Call the writeArray function on every domain equally often, since data is written in chunks with the same
    // hilbertId or chunks of contiguous hilbert ids
    for(MInt h = 0, pOffset = 0; h < maxNoIds; h++) {
      if(h < noIds) {
        // write data with same hilbertId or contigous chunk of hilbertIds, since in the slices cells are sorted by
        // their hilbertId of the slice
        parallelIo.setOffset(hilbertInfo[h * 3 + 1], hilbertInfo[h * 3 + 2]);
        parallelIo.writeArray(&(vars[varId + pOffset * noVars]), varName.str(), noVars);
        // local offset in data array
        pOffset += hilbertInfo[h * 3 + 1];
      } else {
        // dummy calls if current domain has finished writing data, but other domains dont
        parallelIo.setOffset(0, 0);
        parallelIo.writeArray(&(vars[varId]), varName.str(), noVars);
      }
    }
  }
  const MFloat writeTimeTotal = wallTime() - writeTimeStart;
  if(pp->solver().domainId() == 0) {
    std::cerr << "Slice #" << sliceId << " " << fileName.str() << " write time: " << writeTimeTotal << std::endl;
  }
}

/** \brief stores all min level cell ids, TODO labels:PP wrong for partition level shift
 *
 * \author Ansgar Niemoeller
 * \date 07.06.2014
 *
 * \param[in] grid pointer to the grid
 *
 **/
template <MInt nDim, class ppType>
void PostProcessing<nDim, ppType>::collectMinLvlCells() {
  TRACE();

  MInt noCells = pp->solver().grid().noInternalCells();
  ScratchSpace<MInt> minLvlCellIds(noCells, AT_, "minLvlCellIds");
  MInt noMinLvlCells = 0;

  for(MInt cellId = 0; cellId < noCells; cellId++) {
    if(m_gridProxy->tree().parent(cellId) == -1 || m_gridProxy->tree().parent(cellId) >= noCells) { // ask jerry!
      minLvlCellIds[noMinLvlCells++] = cellId;
    }
  }

  m_noMinLvlIds = noMinLvlCells;
  mAlloc(m_minLvlIds, m_noMinLvlIds, "m_minLvlIds", AT_);

  for(MInt cellId = 0; cellId < noMinLvlCells; cellId++) {
    m_minLvlIds[cellId] = minLvlCellIds[cellId];
  }
}

/** \brief get id of leaf cell containing a given point
 *
 * \author Ansgar Niemoeller
 * \date 07.06.2014
 *
 * Important: this function requires a previous call to collectMinLvlCells(...) \n
 * starts with the min level cells and searches down to the leaf cell containing a given point \n
 * returns -1 if no cell is found \n
 * if the given point lies on a cell boundary and between different domains it is considered
 * to belong only to the domain with the lowest id (cellId==-1 for all other domains)
 *
 * \param[in] coord point of interest
 * \param[out] cellId id of leaf cell containing coord
 * \param[in] grid pointer to the grid
 *
 **/
template <MInt nDim, class ppType>
void PostProcessing<nDim, ppType>::findContainingCell(const MFloat* coord, MInt& cellId) {
  TRACE();

  MInt noMinLvlCells = m_noMinLvlIds;
  MInt curId, childCellId;
  MFloat halfCellLength;
  const MFloat* cellCoords;

  cellId = -1;
  // find min level cell containing coord
  for(MInt i = 0; i < noMinLvlCells; i++) {
    MInt cnt = 0;
    curId = m_minLvlIds[i];
    cellCoords = &pp->solver().a_coordinate(curId, 0);
    halfCellLength = m_gridProxy->halfCellLength(curId);

    for(MInt dimId = 0; dimId < nDim; dimId++) {
      if(abs(coord[dimId] - cellCoords[dimId]) <= halfCellLength) cnt++;
    }
    if(cnt == nDim) {
      cellId = curId;
      break;
    }
  }

  if(cellId == -1) return;

  // loop down to leaf cell containing coord
  while(m_gridProxy->tree().noChildren(cellId) > 0) {
    MInt childId;
    for(childId = 0; childId < IPOW2(nDim); childId++) {
      MInt cnt = 0;
      childCellId = m_gridProxy->tree().child(cellId, childId);
      if(childCellId < 0) continue;
      cellCoords = &pp->solver().a_coordinate(childCellId, 0);
      halfCellLength = m_gridProxy->halfCellLength(childCellId);

      for(MInt dimId = 0; dimId < nDim; dimId++) {
        if(abs(coord[dimId] - cellCoords[dimId]) <= halfCellLength) cnt++;
      }

      if(cnt == nDim) {
        cellId = childCellId;
        break;
      }
    }
    if(childId == IPOW2(nDim)) {
      cellId = -1;
      return;
    }
  }

  // check if coord is on cell face and determine if neighboring cell is in another domain
  // if coord lies between domains it is considered to belong to the domain with the lowest id
  MBool cellInOtherDomain = 0;
  MFloat distance;
  MInt direction;
  MInt nghbrId, globalId, parentId;
  if(cellId != -1) {
    cellCoords = &pp->solver().a_coordinate(cellId, 0);
    halfCellLength = m_gridProxy->halfCellLength(cellId);
    for(MInt dimId = 0; dimId < nDim; dimId++) {
      distance = abs(coord[dimId] - cellCoords[dimId]);
      if(distance >= halfCellLength && distance < halfCellLength + MFloatEps) {
        direction = (coord[dimId] - cellCoords[dimId] < 0) ? 0 : 1;

        globalId = -1;
        // check for neighbor
        if(m_gridProxy->tree().hasNeighbor(cellId, 2 * dimId + direction)) { // same level neighbor
          nghbrId = m_gridProxy->tree().neighbor(cellId, 2 * dimId + direction);
          globalId = m_gridProxy->tree().globalId(nghbrId);
        } else { // check for neighbor on parent level
          parentId = m_gridProxy->tree().parent(cellId);
          if(parentId != -1) {
            if(m_gridProxy->tree().hasNeighbor(parentId, 2 * dimId + direction)) {
              nghbrId = m_gridProxy->tree().neighbor(parentId, 2 * dimId + direction);
              globalId = m_gridProxy->tree().globalId(nghbrId);
            }
          }
        }

        if(globalId != -1) {
          // check if neighboring cell is in a "lower" domain
          if(pp->solver().grid().raw().findNeighborDomainId(globalId) < pp->solver().domainId()) {
            cellInOtherDomain = 1;
          }
        }
      }
    }
    if(cellInOtherDomain) {
      cellId = -1;
    }
  }
}


/// \brief Reads the coordinates of points for which the velocities and pressure should be saved
///
/// \author Bjoern Peeters (bjoern) <b.peeters@aia.rwth-aachen.de>
/// \date 2017-05-01
///
/// \param[in] grid Pointer to the grid
///
template <MInt nDim, class ppType>
void PostProcessing<nDim, ppType>::initWritePointData() {
  TRACE();

  /*! \page propertiesPP
    \section pp_pdFileName
    <code>MInt PostProcesssingSolver::m_pdFileName</code>\n
    default = <code>""</code>\n\n
    This property determines the filename for the input coordinates.\n
    <ul>
    <li><code>filename</code></li>
    </ul>\n
    Keywords: <i>GLOBAL, POSTPROCESSING</i>
  */
  m_pdFileName = "";
  m_pdFileName = Context::getSolverProperty<MString>("pp_pdFileName", m_postprocessingId, AT_, &m_pdFileName);

  /*! \page propertiesPP
    \section pp_pdStartTimestep
    <code>MInt PostProcesssingSolver::m_pdStartTimestep</code>\n
    default = <code>0</code>\n\n
    This property determines the timestep after which the storing of the point data starts.\n
    Keywords: <i>GLOBAL, POSTPROCESSING</i>
  */
  m_pdStartTimestep = 0;
  m_pdStartTimestep = Context::getSolverProperty<MInt>("pp_pdStartTimestep", m_postprocessingId, AT_);

  /*! \page propertiesPP
    \section pp_pdStopTimestep
    <code>MInt PostProcesssingSolver::m_pdStopTimestep</code>\n
    default = <code>0</code>\n\n
    This property determines the timestep after which the storing of the point data stops.\n
    Keywords: <i>GLOBAL, POSTPROCESSING</i>
  */
  m_pdStopTimestep = 0;
  m_pdStopTimestep = Context::getSolverProperty<MInt>("pp_pdStopTimestep", m_postprocessingId, AT_);

  /*! \page propertiesPP
    \section pp_pdRestartInterval
    <code>MInt PostProcesssingSolver::m_pdRestartInterval</code>\n
    default = <code>0</code>\n\n
    This property determines the interval after which the storing of the point data writes a restart
    file.\n
    Keywords: <i>GLOBAL, POSTPROCESSING</i>
  */
  m_pdRestartInterval = 0;
  m_pdRestartInterval = Context::getSolverProperty<MInt>("pp_pdRestartInterval", m_postprocessingId, AT_);

  IF_CONSTEXPR(nDim == 3) {
    TERMM(1, "This function is untested for 3D. First make sure it works, then delete this warning...");
  }
  // Store u,v, (w) and p
  MInt noVars = nDim + 1;
  MInt noCoords = 0;
  MInt noPoints = 0;
  MInt noTimesteps = m_pdStopTimestep - m_pdStartTimestep;
  if(noTimesteps < 0) {
    TERMM(1, "The stop timestep for the writing of point data must be greater than the start timestep.");
  }

  collectMinLvlCells();

  vector<MFloat> tmpCoordinates{};
  if(pp->solver().domainId() == 0) {
    MString line;
    MFloat curFloat;
    istringstream iss;
    ifstream csvFile(m_pdFileName);
    // Read all Lines and get coordniates
    while(getline(csvFile, line)) {
      iss.str(line);
      iss.clear();

      for(MInt i = 0; i < nDim; i++) {
        iss >> curFloat;
        if(iss.fail()) {
          ostringstream err;
          // Start line count at one
          err << "Error at line " << noPoints + 1 << ": " << line << "\n"
              << "Either wrong dimension (nDim = " << nDim << ") or otherwise wrong format."
              << "Format should be nDim floats seperated by spaces per line.";
          TERMM(1, err.str());
        }
        tmpCoordinates.push_back(curFloat);
      }
      noCoords++;
    }
  }
  MPI_Bcast(&noCoords, 1, MPI_INT, 0, pp->solver().mpiComm(), AT_, "noCoords");

  vector<MFloat> points(noCoords * nDim);
  if(pp->solver().domainId() == 0) {
    copy(tmpCoordinates.begin(), tmpCoordinates.end(), points.begin());
  }
  MPI_Bcast(&points[0], noCoords * nDim, MPI_DOUBLE, 0, pp->solver().mpiComm(), AT_, "points[0]");

  // Find the containing cells for the coordinates and save their cell id
  MInt cellId = 0;
  for(MInt pointId = 0; pointId < noCoords; pointId++) {
    findContainingCell(&points[nDim * pointId], cellId);

    if(cellId != -1) {
      noPoints++;
      m_pdPoints.insert(m_pdPoints.end(), &points[nDim * pointId], &points[nDim * (pointId + 1)]);
      m_pdCells.push_back(cellId);
    }
  }

  m_pdNoPoints = noPoints;
  m_pdVars.resize(noPoints * noTimesteps * noVars, 0.0);
}


/// \brief Stores the velocities and pressure for all points and writes them to a Netcdf file with
/// the prefix "microphones_"
///
/// \author Bjoern Peeters (bjoern) <b.peeters@aia.rwth-aachen.de>
/// \date 2017-05-01
///
template <MInt nDim, class ppType>
void PostProcessing<nDim, ppType>::writePointData() {
  TRACE();

  IF_CONSTEXPR(nDim == 3) {
    TERMM(1, "This function is untested for 3D. First make sure it works, then delete this warning...");
  }
  // Store u,v, (w) and p
  MInt noVars = nDim + 1;
  MInt noPoints = (m_pdNoPoints > 0) ? m_pdNoPoints : 1;
  MInt noTimesteps = m_pdStopTimestep - m_pdStartTimestep;

  if(globalTimeStep > m_pdStartTimestep && globalTimeStep <= m_pdStopTimestep) {
    if(m_pdNoPoints > 0) {
      MInt curTimestep = globalTimeStep - m_pdStartTimestep;
      for(MInt pointId = 0; pointId < m_pdNoPoints; pointId++) {
        const MFloat* cellVars = 0;
        // pp->solver().getSampleVariables(m_pdCells[pointId], cellVars);
        getSampleVariables(m_pdCells[pointId], cellVars, true);
        MInt startPos = pointId * (noTimesteps * noVars) + (curTimestep - 1) * noVars;
        if(cellVars == nullptr) continue;
        // Store the velocities
        for(MInt dimId = 0; dimId < nDim; dimId++) {
          m_pdVars[startPos + dimId] = cellVars[dimId];
        }
        // Store p
        m_pdVars[startPos + nDim] = cellVars[nDim + 1];
      }
    }
  }

  if(globalTimeStep == m_pdStopTimestep || (globalTimeStep % m_pdRestartInterval == 0)) {
    // Save to Netcdf file
    MString outputDir;
    outputDir = Context::getSolverProperty<MString>("outputDir", m_postprocessingId, AT_);

    ParallelIo::size_type pointsOffset, totalNoPoints;
    stringstream fileName;
    if(globalTimeStep != m_pdStopTimestep) {
      fileName << outputDir << "microphones_" << m_pdStartTimestep << "_" << m_pdStopTimestep << "_" << globalTimeStep
               << ".Netcdf";
    } else {
      fileName << outputDir << "microphones_" << m_pdStartTimestep << "_" << m_pdStopTimestep << ".Netcdf";
    }
    ParallelIo parallelIo(fileName.str(), maia::parallel_io::PIO_REPLACE, pp->solver().mpiComm());
    parallelIo.calcOffset(m_pdNoPoints, &pointsOffset, &totalNoPoints, pp->solver().mpiComm());


    MString dimNames[3] = {"x", "y", "z"};

    for(MInt dimId = 0; dimId < nDim; dimId++) {
      stringstream dimName;
      dimName << "coordinates" << dimId;
      parallelIo.defineArray(maia::parallel_io::PIO_FLOAT, dimName.str(), totalNoPoints);
      parallelIo.setAttribute(dimNames[dimId], "name", dimName.str());
    }
    parallelIo.defineArray(maia::parallel_io::PIO_FLOAT, "cellIds", totalNoPoints);

    ParallelIo::size_type dimSizes[] = {totalNoPoints, noTimesteps, noVars};
    parallelIo.defineArray(maia::parallel_io::PIO_FLOAT, "microphones", 3, &dimSizes[0]);
    parallelIo.setAttribute("desc", "description", "microphones");
    parallelIo.setAttribute("timesteps", "dim_0", "microphones");
    parallelIo.setAttribute("points", "dim_1", "microphones");
    parallelIo.setAttribute("vars", "dim_2", "microphones");
    parallelIo.setAttribute("u", "var_0", "microphones");
    parallelIo.setAttribute("v", "var_1", "microphones");
    IF_CONSTEXPR(nDim == 2) { parallelIo.setAttribute("p", "var_2", "microphones"); }
    else {
      parallelIo.setAttribute("w", "var_2", "microphones");
      parallelIo.setAttribute("p", "var_3", "microphones");
    }

    parallelIo.setOffset(m_pdNoPoints, pointsOffset);


    // Write coordinates as x and y
    for(MInt dimId = 0; dimId < nDim; dimId++) {
      ScratchSpace<MFloat> coord(noPoints, AT_, "coordinates");
      for(MInt i = 0; i < m_pdNoPoints; i++) {
        coord[i] = m_pdPoints[nDim * i + dimId];
      }
      stringstream dimName;
      dimName << "coordinates" << dimId;
      parallelIo.writeArray(&(coord[0]), dimName.str());
    }

    // Write cellIds
    ScratchSpace<MFloat> cellIds(noPoints, AT_, "cellIds");
    for(MInt i = 0; i < m_pdNoPoints; i++) {
      cellIds[i] = m_pdCells[i];
    }
    parallelIo.writeArray(&(cellIds[0]), "cellIds");

    // Write values at microphones
    ScratchSpace<MFloat> microphones(noPoints, noTimesteps, noVars, AT_, "microphones_data");
    if(m_pdNoPoints > 0) {
      for(MInt i = 0; i < noPoints; i++) {
        for(MInt j = 0; j < noTimesteps; j++) {
          for(MInt k = 0; k < noVars; k++) {
            microphones(i, j, k) = m_pdVars[i * (noTimesteps * noVars) + j * noVars + k];
          }
        }
      }
    }

    parallelIo.setOffset(m_pdNoPoints, pointsOffset, 3);
    parallelIo.writeArray(&microphones[0], "microphones");
  }
}


/** \brief Searches for the nearest boundary cell from the given coordinate.
 *
 * \author Jerry Grimmen
 * \date 08.08.2016
 *
 * \param[out] id of the cell
 * \param[in] coord to check
 * \param[in] grid pointer to the grid
 *
 **/
template <MInt nDim, class ppType>
MInt PostProcessing<nDim, ppType>::findNearestGridCell(const MFloat* coord) {
  TRACE();

  MInt out_cellId = -1;
  MFloat t_distance = 0.0;

  m_log << "grid size is " << m_gridProxy->tree().size() << endl;

  // Loop over all (boundary) cells
  for(MInt i = 0; i < m_gridProxy->tree().size(); ++i) {
    if(m_gridProxy->tree().noChildren(i) > 0) {
      continue;
    }
    if(m_gridProxy->tree().hasProperty(i, Cell::IsHalo)) {
      continue;
    }
    if(m_gridProxy->tree().globalId(i) < 0) {
      continue;
    }
    // Get Distance to the cell
    MFloat i_distance = 0.0;

    i_distance += (coord[0] - pp->solver().a_coordinate(i, 0)) * (coord[0] - pp->solver().a_coordinate(i, 0));
    i_distance += (coord[1] - pp->solver().a_coordinate(i, 1)) * (coord[1] - pp->solver().a_coordinate(i, 1));

    IF_CONSTEXPR(nDim == 3) {
      i_distance += (coord[2] - pp->solver().a_coordinate(i, 2)) * (coord[2] - pp->solver().a_coordinate(i, 2));
    }

    if((i_distance < t_distance) || (-1 == out_cellId)) {
      t_distance = i_distance;
      out_cellId = i;
    }
  }

  //  if (m_gridProxy->tree().hasProperty(out_cellId, Cell::IsHalo)) {
  //    out_cellId = -1;
  //  }

  return out_cellId;
}

/** \brief Move the points onto the grid cells.
 * Property m_movePointsToGrid should be set before calling this.
 *
 * \author Jerry Grimmen
 * \date 15.08.2016
 *
 * \param[in] coordinates of the points to move
 * \param[in] number of coordinates
 * \param[in] grid pointer to the grid
 *
 **/
template <MInt nDim, class ppType>
void PostProcessing<nDim, ppType>::movePointsToGrid(MFloat* in_points, MInt in_noPoints, MInt in_moveType) {
  TRACE();

  struct MFloatInt {
    MFloat val;
    MInt rank;
  };

  MFloatInt* distRank = new MFloatInt[in_noPoints / nDim];

  MIntScratchSpace gridCellId(in_noPoints / nDim, AT_, "gridPointId");
  gridCellId.fill(-1);

  for(MInt i = 0; i < in_noPoints; i += nDim) {
    MFloat checkPoint[3];

    checkPoint[0] = in_points[i];
    checkPoint[1] = in_points[i + 1];
    IF_CONSTEXPR(nDim == 3) { checkPoint[2] = in_points[i + 2]; }

    m_log << "checkPoint is " << checkPoint[0] << ", " << checkPoint[1] << ", " << checkPoint[2] << endl;

    if(1 == in_moveType) {
      gridCellId[i / nDim] = pp->solver().grid().raw().findContainingLeafCell(checkPoint);
    }

    if(-1 == gridCellId[i / nDim]) {
      gridCellId[i / nDim] = findNearestGridCell(checkPoint);

      m_log << "gridCellId is " << gridCellId[i / nDim] << endl;

      distRank[i / nDim].val = 0.0;
      distRank[i / nDim].val += (in_points[i] - pp->solver().a_coordinate(gridCellId[i / nDim], 0))
                                * (in_points[i] - pp->solver().a_coordinate(gridCellId[i / nDim], 0));
      distRank[i / nDim].val += (in_points[i + 1] - pp->solver().a_coordinate(gridCellId[i / nDim], 1))
                                * (in_points[i + 1] - pp->solver().a_coordinate(gridCellId[i / nDim], 1));
      IF_CONSTEXPR(nDim == 3) {
        distRank[i / nDim].val += (in_points[i + 2] - pp->solver().a_coordinate(gridCellId[i / nDim], 2))
                                  * (in_points[i + 2] - pp->solver().a_coordinate(gridCellId[i / nDim], 2));
      }
      distRank[i / nDim].rank = pp->solver().domainId();
    }
  }

  m_log << "Before correction old coordinates for each points are..." << endl;
  for(MInt i = 0; i < in_noPoints; i += nDim) {
    m_log << "Coordinate for point " << i / nDim << " is " << in_points[i] << ", " << in_points[i + 1] << ", "
          << in_points[i + 2] << endl;
  }

  if(1 < pp->solver().noDomains()) {
    m_log << "Before Communication" << endl;
    for(MInt i = 0; i < in_noPoints; i += nDim) {
      m_log << "Distance for point " << i / nDim << " is " << distRank[i / nDim].val << " on rank "
            << distRank[i / nDim].rank << endl;
    }

    // Check which rank has the lowest distance to each point
    MPI_Allreduce(MPI_IN_PLACE, distRank, in_noPoints / nDim, MPI_DOUBLE_INT, MPI_MINLOC, pp->solver().mpiComm(), AT_,
                  "MPI_IN_PLACE", "distRank");

    m_log << "After Communication" << endl;
    for(MInt i = 0; i < in_noPoints; i += nDim) {
      m_log << "Smallest distance for point " << i / nDim << " is " << distRank[i / nDim].val << " on rank "
            << distRank[i / nDim].rank << endl;
    }
  }

  for(MInt i = 0; i < in_noPoints; i += nDim) {
    // If I'm the right rank, move the point and broadcast the new coordinates
    if(distRank[i / nDim].rank == pp->solver().domainId()) {
      in_points[i] = pp->solver().a_coordinate(gridCellId[i / nDim], 0);
      in_points[i + 1] = pp->solver().a_coordinate(gridCellId[i / nDim], 1);
      IF_CONSTEXPR(nDim == 3) { in_points[i + 2] = pp->solver().a_coordinate(gridCellId[i / nDim], 2); }

      if(1 < pp->solver().noDomains()) {
        MPI_Bcast(&in_points[i], nDim, MPI_DOUBLE, distRank[i / nDim].rank, pp->solver().mpiComm(), AT_,
                  "in_points[i]");
      }
    } else {
      // If I'm not the right rank, receive the new coordinates
      if(1 < pp->solver().noDomains()) {
        MPI_Bcast(&in_points[i], nDim, MPI_DOUBLE, distRank[i / nDim].rank, pp->solver().mpiComm(), AT_,
                  "in_points[i]");
      }
    }
  }

  m_log << "After correction new coordinates for each points are..." << endl;
  for(MInt i = 0; i < in_noPoints; i += nDim) {
    m_log << "Coordinate for point " << i / nDim << " is " << in_points[i] << ", " << in_points[i + 1] << ", "
          << in_points[i + 2] << endl;
  }

  delete[] distRank;
  distRank = nullptr;
}

/** \brief init function for sray data
 * \author Tim Wegmann
 * \date Jan 2022
 **/
template <MInt nDim, class ppType>
void PostProcessing<nDim, ppType>::initSprayData() {
  m_sprayComputeInterval =
      Context::getSolverProperty<MInt>("postProcTimeInterval", m_postprocessingId, AT_, &m_sprayComputeInterval);

  /*! \page propertiesPP
    \section particleOutputStep
    <code>MInt LPT::m_outputStep</code>\n
    default = <code>50</code> \n \n
    Number of time steps between output of particle data. \n
    Keywords: <i>PARTICLE</i>
 */
  m_sprayWriteInterval =
      Context::getSolverProperty<MInt>("particleOutputStep", m_postprocessingId, AT_, &m_sprayWriteInterval);

  m_sprayDataSize = m_sprayWriteInterval / m_sprayComputeInterval;

  if(m_sprayDataSize * m_sprayComputeInterval != m_sprayWriteInterval) {
    cerr << m_sprayDataSize << " " << m_sprayWriteInterval << " " << m_sprayComputeInterval << endl;
    mTerm(1, AT_, "Missmatching spray intervals");
  }
}

/** \brief Reduces the current grid and saves the data
 *
 * \author Andreas Lintermann, Thomas Schilden
 * \date 26.08.2012
 *
 * The current grid is reduced in the solver and the IO routine in the solver is used to save the data.
 *
 *
 **/
template <MInt nDim, class ppType>
void PostProcessing<nDim, ppType>::pp_saveCoarseSolution() {
  TRACE();
#if not defined(MAIA_DISABLE_LB)
  if(string2enum(pp->solver().solverType()) == MAIA_LATTICE_BOLTZMANN) {
    pp->solver().saveCoarseSolution();
  }
#else
  TERM(-1);
#endif
}
