// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "fvstructuredpostprocessing.h"

#include "fvstructuredcell.h"
#include "fvstructuredsolver.h"
#include "globals.h"

using namespace std;

/**
 * \brief Constructor for the postprocessing solver
 *
 * \author Andreas Lintermann
 * \date 12.09.2012
 *
 * Reads in options for postprocessing by calling initProcessingSolver().
 *
 *
 **/
template <MInt nDim, class SolverType>
StructuredPostprocessing<nDim, SolverType>::StructuredPostprocessing() : m_postprocessing(false) {
  TRACE();
}

/**
 * \brief Destructor for the postprocessing solver
 *
 * \author Andreas Lintermann
 * \date 26.08.2012
 *
 * \tparam[in] T celltype
 *
 **/
template <MInt nDim, class SolverType>
StructuredPostprocessing<nDim, SolverType>::~StructuredPostprocessing() {
  TRACE();

  if(m_postprocessing) {
    if(m_noPostprocessingOps > 0)
      for(MInt op = 0; op < m_noPostprocessingOps; op++) {
        switch(string2enum(m_postprocessingOps[op])) {
          case PP_AVERAGE_PRE:
          case PP_AVERAGE_IN:
          case PP_AVERAGE_POST:
          case PP_TAUW_PRE:
          case PP_LOAD_AVERAGED_SOLUTION_PRE:
          case PP_SUBTRACT_MEAN:
          case PP_WRITE_GRADIENTS:
          case PP_DECOMPOSE_CF: {
            mDeallocate(m_summedVars);
            mDeallocate(m_square);
            if(m_useKahan) {
              mDeallocate(m_cSum);
              mDeallocate(m_ySum);
              mDeallocate(m_tSum);

              mDeallocate(m_cSquare);
              mDeallocate(m_ySquare);
              mDeallocate(m_tSquare);

              if(m_kurtosis) {
                mDeallocate(m_cCube);
                mDeallocate(m_yCube);
                mDeallocate(m_tCube);

                mDeallocate(m_cFourth);
                mDeallocate(m_yFourth);
                mDeallocate(m_tFourth);
              } else if(m_skewness) {
                mDeallocate(m_cCube);
                mDeallocate(m_yCube);
                mDeallocate(m_tCube);
              }
            }
            if(m_kurtosis) {
              mDeallocate(m_cube);
              mDeallocate(m_fourth);
            } else if(m_skewness) {
              mDeallocate(m_cube);
            }
            break;
          }
          case PP_COMPUTE_PRODUCTION_TERMS_PRE: {
            mDeallocate(m_production);
            break;
          }
          case PP_COMPUTE_DISSIPATION_TERMS_PRE: {
            mDeallocate(m_dissipation);
            break;
          }
          case PP_MOVING_AVERAGE_IN:
          case PP_MOVING_AVERAGE_PRE:
          case PP_MOVING_AVERAGE_POST: {
            break;
          }
          default: {
            mTerm(1, AT_, "Unknown postprocessing operation");
          }
        }
      }
    delete[] m_postprocessingOps;
  }
}

template <MInt nDim, class SolverType>
void StructuredPostprocessing<nDim, SolverType>::initStructuredPostprocessing() {
  TRACE();


  m_movingGrid = false;
  m_movingGrid = Context::getSolverProperty<MBool>("movingGrid", ppsolver()->solverId(), AT_, &m_movingGrid);


  /*! \property
    \page propertiesFVSTRCTRD
    \section postprocessing
    <code>MInt StructuredPostprocessing::postprocessing</code>\n
    default = <code>0</code>\n\n
    This property determines the postrprocessing.
    <ul>
    <li><code>0</code> deactivated</li>
    <li><code>1</code> active</li>
    </ul>\n
    Keywords: <i>GENERAL, GLOBAL, POSTPROCESSING</i>
  */
  m_postprocessing = false;
  m_postprocessing =
      Context::getSolverProperty<MBool>("postprocessing", ppsolver()->solverId(), AT_, &m_postprocessing);

  /*! \property
    \page propertiesFVSTRCTRD
    \section pp_skewness
    <code>MInt PostprocesssingSolver::m_skewness</code>\n
    default = <code>0</code>\n\n
    This propertpy determines if skewness is computed when PP_AVERAGE_PRE/IN/POST is activated.
    <ul>
    <li><code>0</code> deactivated</li>
    <li><code>1</code> active</li>
    </ul>\n
    Keywords: <i>GLOBAL, POSTPROCESSING</i>
  */

  m_averageVorticity = 0;
  m_averageVorticity =
      Context::getSolverProperty<MBool>("pp_averageVorticity", ppsolver()->solverId(), AT_, &m_averageVorticity);

  m_skewness = false;
  m_skewness = Context::getSolverProperty<MBool>("pp_skewness", ppsolver()->solverId(), AT_, &m_skewness);

  /*! \property
    \page propertiesFVSTRCTRD
    \section pp_kurtosis
    <code>MInt PostprocesssingSolver::m_kurtosis</code>\n
    default = <code>0</code>\n\n
    This property determines if kurtosis (and skewness) is computed when PP_AVERAGE_PRE/IN/POST is activated.
    <ul>
    <li><code>0</code> deactivated</li>
    <li><code>1</code> active</li>
    </ul>\n
    Keywords: <i>GLOBAL, POSTPROCESSING</i>
  */
  m_kurtosis = false;
  m_kurtosis = Context::getSolverProperty<MBool>("pp_kurtosis", ppsolver()->solverId(), AT_, &m_kurtosis);
  if(m_kurtosis) m_skewness = 1;

  /*! \property
    \page propertiesFVSTRCTRD
    \section pp_turbulentProduction
    <code>MInt PostprocesssingSolver::m_computeProductionTerms</code>\n
    default = <code>0</code>\n\n
    Determines if the turbulent production terms should be computed after the normal averaging
    <ul>
    <li><code>0</code> deactivated</li>
    <li><code>1</code> active</li>
    </ul>\n
    Keywords: <i>GLOBAL, POSTPROCESSING</i>
  */
  m_computeProductionTerms = false;
  m_computeProductionTerms = Context::getSolverProperty<MBool>(
      "pp_turbulentProduction", ppsolver()->solverId(), AT_, &m_computeProductionTerms);

  //  if(m_kurtosis==1 || m_skewness==1)
  //    mTerm(1,AT_,"check fvsolver and all computations because pressure ampl was introduced at position 11");

  /*! \property
    \page propertiesFVSTRCTRD
    \section pp_twoPass
    <code>MInt PostprocesssingSolver::m_twoPass</code>\n
    default = <code>0</code>\n\n
    This property determines if two-pass averaging is performed in PP_AVERAGE_PRE/POST.\n
    Either m_twoPass or m_useKahan should be activated.
    <ul>
    <li><code>0</code> deactivated</li>
    <li><code>1</code> active</li>
    </ul>\n
    Keywords: <i>GLOBAL, POSTPROCESSING</i>
  */
  m_twoPass = false;
  m_twoPass = Context::getSolverProperty<MBool>("pp_twoPass", ppsolver()->solverId(), AT_, &m_twoPass);

  /*! \property
    \page propertiesFVSTRCTRD
    \section pp_useKahan
    <code>MInt PostprocesssingSolver::m_useKahan</code>\n
    default = <code>0</code>\n\n
    This property determines if kahan summation is performed in PP_AVERAGE_PRE/IN/POST. \n
    Either m_twoPass or m_useKahan should be activated.
    <ul>
    <li><code>0</code> deactivated</li>
    <li><code>1</code> active</li>
    </ul>\n
    Keywords: <i>GLOBAL, POSTPROCESSING</i>
  */
  m_useKahan = false;
  m_useKahan = Context::getSolverProperty<MBool>("pp_useKahan", ppsolver()->solverId(), AT_, &m_useKahan);

  /*! \property
    \page propertiesFVSTRCTRD
    \section pp_fileName
    <code>MInt PostprocesssingSolver::m_postprocessFileName</code>\n
    default = <code>""</code>\n\n
    This property determines a filename for averaging.\n
    <ul>
    <li><code>filename</code></li>
    </ul>\n
    Keywords: <i>GLOBAL, POSTPROCESSING</i>
  */


  m_postprocessFileName = "";
  m_postprocessFileName =
      Context::getSolverProperty<MString>("pp_fileName", ppsolver()->solverId(), AT_, &m_postprocessFileName);

  // Init vars
  m_averageInterval = 0;
  m_averageRestart = 0;
  m_averageRestartInterval = 0;
  m_averageStartTimestep = 0;
  m_averageStopTimestep = 0;
  m_noPostprocessingOps = 0;
  m_noSamples = 0;

  m_noVariables = ppsolver()->noVariables();

  if(m_postprocessing) {
    for(MInt i = 0; i < 3; i++) {
      tvecpost tmp;
      vector<MString> st;
      m_postprocessingMethods.push_back(tmp);
      m_postprocessingMethodsDesc.push_back(st);
    }

    /*! \property
      \page propertiesFVSTRCTRD
      \section postprocessingOps
      <code>MString* Solver::m_postprocessingOps</code>\n
      default = <code>empty</code>\n\n
      This property is a list of postprocessing operations to be performed
      <ul>
      <li><code>PP_AVERAGE_PRE</code> </li>
      <li><code>PP_AVERAGE_POST</code> </li>
      <li><code>PP_AVERAGE_IN</code> </li>
      <li><code>PP_MOVING_AVERAGE_PRE</code> </li>
      <li><code>PP_MOVING_AVERAGE_POST</code> </li>
      <li><code>PP_MOVING_AVERAGE_IN</code> </li>
      </ul>\n
      Keywords: <i>GENERAL, GLOBAL, POSTPROCESSING</i>
    */
    m_noPostprocessingOps = Context::propertyLength("postprocessingOps", ppsolver()->solverId());
    m_postprocessingOps = nullptr;
    if(m_noPostprocessingOps > 0) {
      // mAlloc( m_postprocessingOps, m_noPostprocessingOps, "m_postprocessingOps",  AT_ );
      m_postprocessingOps = new MString[m_noPostprocessingOps];
      for(MInt op = 0; op < m_noPostprocessingOps; op++) {
        m_postprocessingOps[op] = Context::getSolverProperty<MString>(
            "postprocessingOps", ppsolver()->solverId(), AT_, &m_postprocessingOps[op], op);

        switch(string2enum(m_postprocessingOps[op])) {
          case PP_AVERAGE_IN: {
            m_postprocessingMethodsDesc[1].push_back(m_postprocessingOps[op]);
            m_postprocessingMethods[1].push_back(&StructuredPostprocessing::averageSolutionsInSolve);
            initAverageIn();
            initAverageVariables();
            break;
          }
          case PP_AVERAGE_PRE: {
            m_postprocessingMethodsDesc[0].push_back(m_postprocessingOps[op]);
            m_postprocessingMethods[0].push_back(&StructuredPostprocessing::averageSolutions);
            initTimeStepProperties();
            initAverageVariables();
            break;
          }
          case PP_AVERAGE_POST: {
            m_postprocessingMethodsDesc[2].push_back(m_postprocessingOps[op]);
            m_postprocessingMethods[2].push_back(&StructuredPostprocessing::averageSolutions);
            initTimeStepProperties();
            initAverageVariables();
            break;
          }
          case PP_LOAD_AVERAGED_SOLUTION_PRE: {
            m_postprocessingMethodsDesc[0].push_back(m_postprocessingOps[op]);
            m_postprocessingMethods[0].push_back(&StructuredPostprocessing::loadAveragedSolution);
            initTimeStepProperties();
            initAverageVariables();
            break;
          }
          case PP_COMPUTE_PRODUCTION_TERMS_PRE: {
            m_postprocessingMethodsDesc[0].push_back(m_postprocessingOps[op]);
            m_postprocessingMethods[0].push_back(&StructuredPostprocessing::computeProductionTerms);
            initProductionVariables();
            break;
          }
          case PP_COMPUTE_DISSIPATION_TERMS_PRE: {
            m_postprocessingMethodsDesc[0].push_back(m_postprocessingOps[op]);
            m_postprocessingMethods[0].push_back(&StructuredPostprocessing::computeDissipationTerms);
            initDissipationVariables();
            break;
          }
          case PP_TAUW_PRE: {
            m_postprocessingMethodsDesc[0].push_back(m_postprocessingOps[op]);
            m_postprocessingMethods[0].push_back(&StructuredPostprocessing::computeAverageSkinFriction);
            initTimeStepProperties();
            initAverageVariables();
            break;
          }
          case PP_SUBTRACT_PERIODIC_FLUCTUATIONS: {
            m_postprocessingMethodsDesc[0].push_back(m_postprocessingOps[op]);
            m_postprocessingMethods[0].push_back(&StructuredPostprocessing::subtractPeriodicFluctuations);

            mAlloc(m_spanAvg, m_noVariables, ppsolver()->getNoCells(), "m_spanAvg", F0, FUN_);
            initTimeStepProperties();
            initAverageVariables();
            break;
          }
          case PP_SUBTRACT_MEAN: {
            m_postprocessingMethodsDesc[0].push_back(m_postprocessingOps[op]);
            m_postprocessingMethods[0].push_back(&StructuredPostprocessing::subtractMean);
            initTimeStepProperties();
            initAverageVariables();
            break;
          }
          case PP_WRITE_GRADIENTS: {
            m_postprocessingMethodsDesc[0].push_back(m_postprocessingOps[op]);
            m_postprocessingMethods[0].push_back(&StructuredPostprocessing::writeGradients);
            // for 3 dimensions and 9 variables (u,v,w,uu,vv,ww,uv,uw,vw)
            mAlloc(m_gradients, 3 * 9, ppsolver()->getNoCells(), "m_gradients", F0, FUN_);
            initTimeStepProperties();
            initAverageVariables();
            break;
          }
          case PP_DECOMPOSE_CF: {
            m_postprocessingMethodsDesc[0].push_back(m_postprocessingOps[op]);
            m_postprocessingMethods[0].push_back(&StructuredPostprocessing::decomposeCfDouble);
            initTimeStepProperties();
            initAverageVariables();
            break;
          }
          case PP_MOVING_AVERAGE_IN: {
            m_postprocessingMethodsDesc[1].push_back(m_postprocessingOps[op]);
            m_postprocessingMethods[1].push_back(&StructuredPostprocessing::movingAverage);
            initTimeStepProperties();
            initMovingAverage();
            break;
          }
          case PP_MOVING_AVERAGE_PRE: {
            m_postprocessingMethodsDesc[0].push_back(m_postprocessingOps[op]);
            m_postprocessingMethods[0].push_back(&StructuredPostprocessing::movingAveragePost);
            initTimeStepProperties();
            initMovingAverage();
            break;
          }
          case PP_MOVING_AVERAGE_POST: {
            m_postprocessingMethodsDesc[2].push_back(m_postprocessingOps[op]);
            m_postprocessingMethods[2].push_back(&StructuredPostprocessing::movingAveragePost);
            initTimeStepProperties();
            initMovingAverage();
            break;
          }
          default: {
            mTerm(1, AT_, "Unknown postprocessing operation");
          }
        }
      }
    }
  }
}

/**
 *
 * \author Frederik Temme
 * \modified 22.12.2015
 */
template <MInt nDim, class SolverType>
void StructuredPostprocessing<nDim, SolverType>::postprocessPreInit() {
  TRACE();

  if(m_averageRestart) {
    m_restartTimeStep = ppsolver()->m_restartTimeStep;
  }

  if(m_noPostprocessingOps != 0 && m_postprocessing == 1) {
    for(MInt op = 0; op < m_noPostprocessingOps; op++) {
      switch(string2enum(m_postprocessingOps[op])) {
        case PP_AVERAGE_POST:
        case PP_AVERAGE_PRE:
        case PP_AVERAGE_IN: {
          // load the averaging restart
          if(m_restartTimeStep > m_averageStartTimestep && m_restartTimeStep <= m_averageStopTimestep) {
            MString name = ppsolver()->outputDir() + "PostprocessingRestart_";
            MChar buf1[10];
            MChar buf2[10];
            sprintf(buf1, "%d", m_averageStartTimestep);
            sprintf(buf2, "%d", m_restartTimeStep);
            name.append(buf1);
            name += "-";
            name.append(buf2);
            name.append(ppsolver()->m_outputFormat);

            m_log << "\n\n"
                  << "    ^^^^^^^^^^^^^^^ Entering postprocessing mode PreInit ^^^^^^^^^^^^^^^^ \n"
                  << "    ^   - Loading restart for mean flow calculation: " << name << "\n"
                  << "    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ \n\n";

            ppsolver()->loadAverageRestartFile(name.c_str(), m_summedVars, m_square, m_cube, m_fourth);
          }
          break;
        }
        case PP_SUBTRACT_MEAN: {
          if(ppsolver()->domainId() == 0) {
            cout << "Subtracting mean from restart file" << endl;
          }

          break;
        }
        default: {
          // has to be filled
        }
      }
    }
  }
}

template <MInt nDim, class SolverType>
void StructuredPostprocessing<nDim, SolverType>::postprocessPreSolve() {
  TRACE();

  if(m_noPostprocessingOps != 0 && m_postprocessing == 1) {
    m_log << "\n\n"
          << "    ^^^^^^^^^^^^^^^ Entering postprocessing mode PreSolve ^^^^^^^^^^^^^^^ \n"
          << "    ^   - Activated operations are:\n";
    for(MInt op = 0; op < m_noPostprocessingOps; op++)
      m_log << "    ^      + " << m_postprocessingOps[op] << "\n";
    m_log << "    ^   - Running:\n";

    for(MInt op = 0; op < (signed)m_postprocessingMethods[0].size(); op++) {
      m_log << "    ^      + " << m_postprocessingMethodsDesc[0][op] << "\n";
      (this->*(m_postprocessingMethods[0][op]))();
    }
    m_log << "    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ \n" << endl;
  }
}

template <MInt nDim, class SolverType>
void StructuredPostprocessing<nDim, SolverType>::postprocessInSolve() {
  TRACE();

  if(m_noPostprocessingOps != 0 && m_postprocessing == 1) {
    for(MInt op = 0; op < (signed)m_postprocessingMethods[1].size(); op++) {
      (this->*(m_postprocessingMethods[1][op]))();
    }
  }
}

template <MInt nDim, class SolverType>
void StructuredPostprocessing<nDim, SolverType>::postprocessPostSolve() {
  TRACE();

  if(m_noPostprocessingOps != 0 && m_postprocessing == 1) {
    m_log << "\n\n"
          << "    ^^^^^^^^^^^^^^^ Entering postprocessing mode PostSolve ^^^^^^^^^^^^^^^ \n"
          << "    ^   - Activated operations are:\n\n";
    for(MInt op = 0; op < m_noPostprocessingOps; op++)
      m_log << "    ^      + " << m_postprocessingOps[op] << "\n";
    m_log << "    ^   - Running:\n";

    for(MInt op = 0; op < (signed)m_postprocessingMethods[2].size(); op++) {
      m_log << "    ^      + " << m_postprocessingMethodsDesc[2][op] << "\n";
      (this->*(m_postprocessingMethods[2][op]))();
    }
    m_log << "    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ \n" << endl;
  }
}


/**
 *
 * \author Frederik Temme
 * \modified 17.12.2015
 *
 */
template <MInt nDim, class SolverType>
void StructuredPostprocessing<nDim, SolverType>::averageSolutionsInSolve() {
  // Average only at the right timestep
  if(globalTimeStep >= m_averageStartTimestep && globalTimeStep <= m_averageStopTimestep) {
    if((globalTimeStep - m_averageStartTimestep) % m_averageInterval == 0) {
      if(!ppsolver()->isTravelingWave()) {
        addAveragingSample();
      } else {
        ppsolver()->spanwiseWaveReorder();
        addTempWaveSample();
        if(ppsolver()->domainId() == 0) {
          cout << ">>>>> wave sample with interval  " << m_averageInterval
               << " time steps at time step: " << globalTimeStep << " has been added  <<<<<" << endl;
        }
      }
    }
  }

  // Compute the averaged solution and write to file
  if(globalTimeStep == m_averageStopTimestep) {
    saveAveragedSolution(globalTimeStep);
  }
}

/**
 * \brief Adds for the travelling wave setups
 *
 * This method is similar to addAveragingSample but
 * was adapted for the phase-averaging of the travelling
 * wave setups
 *
 * \author Marian Albers (original by Ansgar Niemoeller)
 * \date 15.12.2016
 *
 *
 **/
template <MInt nDim, class SolverType>
void StructuredPostprocessing<nDim, SolverType>::addTempWaveSample() {
  TRACE();
  const MInt noCells = ppsolver()->getNoCells();
  const MInt noAveragedVorticities = (m_averageVorticity != 0) * (nDim * 2 - 3);

  for(MInt cellId = 0; cellId < noCells; cellId++) {
    MInt offset = 0;
    MInt offsetSquare = 0;
    // Primitive variables
    for(MInt varId = 0; varId < m_noVariables; varId++) {
      m_summedVars[varId + offset][cellId] += m_tempWaveSample[varId][cellId];
    }
    offset += m_noVariables;

    // Favre-averaging
    if(m_averagingFavre) {
      const MFloat rho = m_tempWaveSample[nDim][cellId];
      for(MInt varId = 0; varId < m_noVariables; varId++) {
        m_favre[varId][cellId] += m_tempWaveSample[varId][cellId] * rho;
      }
    }

    // Vorticities
    if(m_averageVorticity) {
      for(MInt varId = 0; varId < noAveragedVorticities; varId++) {
        m_summedVars[varId + offset][cellId] += m_tempWaveSample[varId + offset][cellId];
      }
      offset += noAveragedVorticities;
    }

    // squares of velocity components ( uu, vv, ww )
    for(MInt varId = 0; varId < nDim; varId++) {
      m_square[varId + offsetSquare][cellId] += m_tempWaveSample[varId][cellId] * m_tempWaveSample[varId][cellId];
    }
    offsetSquare += nDim;

    // product of different velocity componets ( uv, vw, wu )
    for(MInt varId = 0; varId < 2 * nDim - 3; varId++) {
      m_square[varId + offsetSquare][cellId] +=
          m_tempWaveSample[varId % nDim][cellId] * m_tempWaveSample[(varId + 1) % nDim][cellId];
    }
    offsetSquare += 2 * nDim - 3;

    // square of pressure (pp)
    m_square[offsetSquare][cellId] += m_tempWaveSample[nDim + 1][cellId] * m_tempWaveSample[nDim + 1][cellId];
    offsetSquare += 1;

    // squares of the vorticity
    if(m_averageVorticity) {
      for(MInt varId = 0; varId < noAveragedVorticities; varId++) {
        m_square[offsetSquare + varId][cellId] +=
            m_tempWaveSample[varId + m_noVariables][cellId] * m_tempWaveSample[varId + m_noVariables][cellId];
      }
      offsetSquare += noAveragedVorticities;
    }

    // third and fouth powers of velocity components (skewness and kurtosis)
    if(m_kurtosis) {
      for(MInt varId = 0; varId < nDim; varId++) {
        m_cube[varId][cellId] += pow(m_tempWaveSample[varId][cellId], 3);
        m_fourth[varId][cellId] += pow(m_tempWaveSample[varId][cellId], 4);
      }
    }
    // third power if velocity components (skewness)
    if(m_skewness) {
      for(MInt varId = 0; varId < nDim; varId++) {
        m_cube[varId][cellId] += pow(m_tempWaveSample[varId][cellId], 3);
      }
    }
  }
}


/**
 *
 * \author Frederik Temme
 * \modified 19.1.2016
 *
 */
template <MInt nDim, class SolverType>
void StructuredPostprocessing<nDim, SolverType>::averageSolutions() {
  TRACE();

  /***********************************************************
  in properties: solutionInterval has to be equal to pp_averageInterval
                 if using this method.
  ***********************************************************/

  m_log << "    ^          * Averaging solutions ";


  // set the summationstart for averaging
  MFloat avgStart = m_restartTimeStep;

  for(MInt avgTimestep = avgStart; avgTimestep <= m_averageStopTimestep; avgTimestep += m_averageInterval) {
    // load samples
    stringstream filename;
    filename << ppsolver()->outputDir() << avgTimestep << ppsolver()->m_outputFormat;
    ppsolver()->loadSampleFile(filename.str());
    ppsolver()->exchange();
    ppsolver()->applyBoundaryCondition();
    addAveragingSample();

    // Write average restart file
    if(m_averageRestartInterval != 0
       && (avgTimestep >= m_averageStartTimestep && avgTimestep % m_averageRestartInterval == 0
           && avgTimestep <= m_averageStopTimestep)) {
      saveAverageRestart();
    }
  }

  saveAveragedSolution(m_averageStopTimestep);
}

template <MInt nDim, class SolverType>
void StructuredPostprocessing<nDim, SolverType>::saveAveragedSolution(MInt endTimeStep) {
  TRACE();

  computeAveragedSolution();

  // output filename
  MString name = ppsolver()->outputDir() + "Mean_";
  MChar buf1[10];
  MChar buf2[10];
  sprintf(buf1, "%d", m_averageStartTimestep);
  sprintf(buf2, "%d", endTimeStep);
  name.append(buf1);
  name += "-";
  name.append(buf2);

  m_log << "         ^   saving averaged variables " << name << endl;

  ppsolver()->saveAveragedVariables(name, getNoPPVars(), m_summedVars);
}


/**
 * \brief Computes the mean variables from summed vars
 *
 * Computes the correct averaged solution from all added samples
 * Also computes the RMS components of the velocities, the pressure
 * and the vorticities (if desired)
 *
 * \author Marian Albers (original by Ansgar Niemoeller)
 * \date 01.02.2016
 *
 *
 **/
template <MInt nDim, class SolverType>
void StructuredPostprocessing<nDim, SolverType>::computeAveragedSolution() {
  TRACE();

  const MInt noCells = ppsolver()->getNoCells();
  const MInt noAveragedVorticities = (m_averageVorticity != 0) * (nDim * 2 - 3);

  const MFloat weight = F1 / m_noSamples; // F1/(((m_averageStopTimestep-m_averageStartTimestep)/m_averageInterval)+1);
  MInt offset = 0;
  MInt offsetSquare = 0;

  // mean of summed primitive variables
  for(MInt cellId = 0; cellId < noCells; cellId++) {
    for(MInt varId = 0; varId < m_noVariables; varId++) {
      m_summedVars[varId + offset][cellId] *= weight;
    }
  }
  offset += m_noVariables;

  // mean of summed primitive variables
  if(m_averagingFavre) {
    for(MInt cellId = 0; cellId < noCells; cellId++) {
      const MFloat frhom = F1 / m_summedVars[nDim][cellId];
      for(MInt varId = 0; varId < m_noVariables; varId++) {
        m_favre[varId][cellId] = m_favre[varId][cellId] * weight * frhom;
      }
    }
  }

  // Weighting of summed vorticity variables -> mean
  if(m_averageVorticity) {
    for(MInt cellId = 0; cellId < noCells; cellId++) {
      for(MInt varId = 0; varId < noAveragedVorticities; varId++) {
        m_summedVars[varId + offset][cellId] *= weight;
      }
    }

    offset += noAveragedVorticities;
  }

  // compute mean(u'u'),mean(v'v'),mean(w'w') ( e.g.mean(u'u')=mean(u^2)-(u_mean))^2 )
  for(MInt cellId = 0; cellId < noCells; cellId++) {
    for(MInt varId = 0; varId < nDim; varId++) {
      m_summedVars[varId + offset][cellId] =
          weight * m_square[varId][cellId] - m_summedVars[varId][cellId] * m_summedVars[varId][cellId];
    }
  }
  offset += nDim;
  offsetSquare += nDim;


  // compute mean(u'v'),mean(v'w'),mean(w'u') ( e.g. mean(u'v')=mean(u*v)-u_mean*v_mean )
  for(MInt cellId = 0; cellId < noCells; cellId++) {
    for(MInt varId = 0; varId < 2 * nDim - 3; varId++) {
      m_summedVars[varId + offset][cellId] =
          weight * m_square[varId + offsetSquare][cellId]
          - m_summedVars[varId % nDim][cellId] * m_summedVars[(varId + 1) % nDim][cellId];
    }
  }
  offset += 2 * nDim - 3;
  offsetSquare += 2 * nDim - 3;

  if(m_kurtosis) {
    // compute skewness and kurtosis of velocity components
    // e.g. skewness(u) = mean(u^3) - 3*u_mean*mean(u^2) + 2*u_mean^3
    // e.g. kurtosis(u) = mean(u^4) - 4*u_mean*mean(u^3) + 6*u_mean^2*mean(u^2) - 3*u_mean^4
    for(MInt cellId = 0; cellId < noCells; cellId++) {
      for(MInt varId = 0; varId < nDim; varId++) {
        m_summedVars[varId + offset][cellId] = weight * m_cube[varId][cellId]
                                               - 3 * weight * m_summedVars[varId][cellId] * m_square[varId][cellId]
                                               + 2 * pow(m_summedVars[varId][cellId], 3);

        m_summedVars[varId + offset + nDim][cellId] =
            weight * m_fourth[varId][cellId] - 4 * weight * m_cube[varId][cellId] * m_summedVars[varId][cellId]
            + 6 * weight * m_square[varId][cellId] * m_summedVars[varId][cellId] * m_summedVars[varId][cellId]
            - 3 * pow(m_summedVars[varId][cellId], 4);
      }
    }
    offset += 2 * nDim;

  } else if(m_skewness) {
    // compute skewness of velocity components
    for(MInt cellId = 0; cellId < noCells; cellId++) {
      for(MInt varId = 0; varId < nDim; varId++) {
        m_summedVars[varId + offset][cellId] = weight * m_cube[varId][cellId]
                                               - 3 * weight * m_summedVars[varId][cellId] * m_square[varId][cellId]
                                               + 2 * pow(m_summedVars[varId][cellId], 3);
      }
    }
    offset += nDim;
  }

  // compute p'*p'
  for(MInt cellId = 0; cellId < noCells; cellId++) {
    m_summedVars[offset][cellId] = weight * m_square[offsetSquare][cellId]
                                   - m_summedVars[nDim + 1][cellId] * m_summedVars[nDim + 1][cellId]; // pressure
  }
  offset += 1;
  offsetSquare += 1;

  if(m_averageVorticity) {
    // compute vorticity symmetric rms
    for(MInt cellId = 0; cellId < noCells; cellId++) {
      for(MInt varId = 0; varId < nDim; varId++) {
        m_summedVars[varId + offset][cellId] =
            weight * m_square[varId + offsetSquare][cellId]
            - m_summedVars[varId + m_noVariables][cellId] * m_summedVars[varId + m_noVariables][cellId];
      }
    }

    offset += noAveragedVorticities;
    offsetSquare += noAveragedVorticities;
  }


  // add aditional variables here -> otherwise the code for kurtosis will overwrite the summed vars of the new
  // introduced variables
}

/**
 * \brief Adds one sample to the summedVars
 *
 * Adds a sample of the current time step to the
 * summedVars (and other) fields
 *
 * \author Marian Albers (original by Ansgar Niemoeller)
 * \date 12.12.2016
 *
 *
 **/
template <MInt nDim, class SolverType>
void StructuredPostprocessing<nDim, SolverType>::addAveragingSample() {
  TRACE();

  m_log << "     ^        * Averaging timestep " << globalTimeStep << "\n";
  m_noSamples++;

  const MInt noCells = ppsolver()->getNoCells();
  MFloatScratchSpace cellVars(m_noVariables, AT_, "cellVars");

  const MInt noAveragedVorticities = (m_averageVorticity != 0) * (nDim * 2 - 3);

  if(m_averageVorticity) {
    ppsolver()->computeVorticity();
  }

  for(MInt cellId = 0; cellId < noCells; cellId++) {
    /* List of Variables in m_summedVars after following computation
       0=mean(u)      1=mean(v)      2=mean(w)
       3=mean(rho)    4=mean(p)
       5=mean(u'u')   6=mean(v'v')  7=mean(w'w')
       8=mean(u'v')   9=mean(v'w') 10=mean(w'u')
       11=skew u      12=skew v     13=skew w
       14=kurt u      15=kurt v     16=kurt w
    */

    // Calculation of primitive variables
    ppsolver()->getSampleVariables(cellId, cellVars.begin());

    if(m_useKahan) { // Kahan summation

      /* Kahan summation pseudocode
         sum=0; c=0;
         for i=0 to input.length
         y = input[i] -c;
         t = sum + y;
         c = (t-sum) - y;
         sum = t;
      */

      MInt offsetSquare = 0;
      for(MInt varId = 0; varId < m_noVariables; varId++) { // sum up all primitive variables
        m_ySum[varId][cellId] = cellVars[varId] - m_cSum[varId][cellId];
        m_tSum[varId][cellId] = m_summedVars[varId][cellId] + m_ySum[varId][cellId];
        m_cSum[varId][cellId] = (m_tSum[varId][cellId] - m_summedVars[varId][cellId]) - m_ySum[varId][cellId];
        m_summedVars[varId][cellId] = m_tSum[varId][cellId];
      }
      for(MInt varId = 0; varId < nDim; varId++) { // squares of velocity components (u*u,v*v,w*w)
        m_ySquare[varId][cellId] = (cellVars[varId] * cellVars[varId]) - m_cSquare[varId][cellId];
        m_tSquare[varId][cellId] = m_square[varId][cellId] + m_ySquare[varId][cellId];
        m_cSquare[varId][cellId] = (m_tSquare[varId][cellId] - m_square[varId][cellId]) - m_ySquare[varId][cellId];
        m_square[varId][cellId] = m_tSquare[varId][cellId];
      }
      offsetSquare += 3;
      for(MInt varId = 0; varId < nDim; varId++) { // products of different velocity components in order (u*v,v*w,w*)
        m_ySquare[varId + offsetSquare][cellId] =
            (cellVars[varId % 3] * cellVars[(varId + 1) % 3]) - m_cSquare[varId + offsetSquare][cellId];
        m_tSquare[varId + offsetSquare][cellId] =
            m_square[varId + offsetSquare][cellId] + m_ySquare[varId + offsetSquare][cellId];
        m_cSquare[varId + offsetSquare][cellId] =
            (m_tSquare[varId + offsetSquare][cellId] - m_square[varId + offsetSquare][cellId])
            - m_ySquare[varId + offsetSquare][cellId];
        m_square[varId + offsetSquare][cellId] = m_tSquare[varId + offsetSquare][cellId];
      }
      if(m_kurtosis) { // compute third and fourth power of velocity components for skewness and kurtosis
        for(MInt varId = 0; varId < nDim; varId++) {
          m_yCube[varId][cellId] = pow(cellVars[varId], 3) - m_cCube[varId][cellId];
          m_tCube[varId][cellId] = m_cube[varId][cellId] + m_yCube[varId][cellId];
          m_cCube[varId][cellId] = (m_tCube[varId][cellId] - m_cube[varId][cellId]) - m_yCube[varId][cellId];
          m_cube[varId][cellId] = m_tCube[varId][cellId];

          m_yFourth[varId][cellId] = pow(cellVars[varId], 4) - m_cFourth[varId][cellId];
          m_tFourth[varId][cellId] = m_fourth[varId][cellId] + m_yFourth[varId][cellId];
          m_cFourth[varId][cellId] = (m_tFourth[varId][cellId] - m_fourth[varId][cellId]) - m_yFourth[varId][cellId];
          m_fourth[varId][cellId] = m_tFourth[varId][cellId];
        }
      } else if(m_skewness) { // compute only third power of velocity components for skewness
        for(MInt varId = 0; varId < nDim; varId++) {
          m_yCube[varId][cellId] = pow(cellVars[varId], 3) - m_cCube[varId][cellId];
          m_tCube[varId][cellId] = m_cube[varId][cellId] + m_yCube[varId][cellId];
          m_cCube[varId][cellId] = (m_tCube[varId][cellId] - m_cube[varId][cellId]) - m_yCube[varId][cellId];
          m_cube[varId][cellId] = m_tCube[varId][cellId];
        }
      }

    } else { // normal summation
      // Reset offsets
      MInt offset = 0;
      MInt offsetSquare = 0;

      // Primitive variables
      for(MInt varId = 0; varId < m_noVariables; varId++) {
        m_summedVars[varId + offset][cellId] += cellVars[varId];
      }
      offset += m_noVariables;

      // Favre-averaging
      if(m_averagingFavre) {
        const MFloat rho = cellVars[nDim];
        for(MInt varId = 0; varId < m_noVariables; varId++) {
          m_favre[varId][cellId] += cellVars[varId] * rho;
        }
      }

      // Vorticities
      if(m_averageVorticity) {
        for(MInt varId = 0; varId < noAveragedVorticities; varId++) {
          m_summedVars[varId + offset][cellId] += ppsolver()->getSampleVorticity(cellId, varId);
        }
        offset += noAveragedVorticities;
      }

      for(MInt varId = 0; varId < nDim; varId++) { // squares of velocity components (u*u,v*v(,w*w))
        m_square[varId][cellId] += cellVars[varId] * cellVars[varId];
      }
      offsetSquare += nDim;
      for(MInt varId = 0; varId < 2 * nDim - 3; varId++) { // products of different velocity components
                                                           // (u*v(,v*w,w*u))
        m_square[offsetSquare + varId][cellId] += (cellVars[varId % nDim]) * (cellVars[(varId + 1) % nDim]);
      }
      offsetSquare += 2 * nDim - 3;
      m_square[offsetSquare][cellId] += cellVars[nDim + 1] * cellVars[nDim + 1]; // squares of pressure  p*p
      offsetSquare += 1;

      // add aditional variables here -> otherwise the code for kurtosis will overwrite the summed vars of the new
      // introduced variables
      // squares of the vorticity
      if(m_averageVorticity) {
        for(MInt varId = 0; varId < noAveragedVorticities; varId++) {
          m_square[offsetSquare + varId][cellId] +=
              ppsolver()->getSampleVorticity(cellId, varId) * ppsolver()->getSampleVorticity(cellId, varId);
        }
        offsetSquare += noAveragedVorticities;
      }


      if(m_kurtosis) { // third and fourth powers of velocity components (skewness and kurtosis)
        for(MInt varId = 0; varId < nDim; varId++) {
          m_cube[varId][cellId] += pow(cellVars[varId], 3);
          m_fourth[varId][cellId] += pow(cellVars[varId], 4);
        }
      } else if(m_skewness) { // third powers of velocity components (skewness)
        for(MInt varId = 0; varId < nDim; varId++) {
          m_cube[varId][cellId] += pow(cellVars[varId], 3);
        }
      }
    }
  }
}

/**
 * \brief Loads an averaged file again
 *
 * Loads an Postprocessing Mean file into the corresponding field
 * such that further pp actions can be performed
 *
 * Specify file name in 'pp_fileName' in the property file and
 * set 'pp_averageRestart' if this is an PostprocessingRestart and not
 * a Mean file
 *
 * \author Marian Albers
 * \date 01.02.2016
 *
 *
 **/
template <MInt nDim, class SolverType>
void StructuredPostprocessing<nDim, SolverType>::loadAveragedSolution() {
  if(ppsolver()->domainId() == 0) {
    cout << "Loading postprocessing file " << m_postprocessFileName << endl;
  }
  if(m_averageRestart) {
    ppsolver()->loadAverageRestartFile(m_postprocessFileName.c_str(), m_summedVars, m_square, m_cube, m_fourth);
    computeAveragedSolution();
  } else {
    if(ppsolver()->domainId() == 0) {
      cout << "Loading file " << m_postprocessFileName << endl;
    }
    ppsolver()->loadAveragedVariables(m_postprocessFileName.c_str());
  }

  if(ppsolver()->domainId() == 0) {
    cout << "Filling ghost-cells..." << endl;
  }
  vector<MFloat*> ppVariables;
  for(MInt var = 0; var < getNoPPVars(); var++) {
    ppVariables.push_back(m_summedVars[var]);
  }
  ppsolver()->gcFillGhostCells(ppVariables);

  if(ppsolver()->domainId() == 0) {
    cout << "Filling ghost-cells... FINISHED!" << endl;
  }

  const MInt noCells = ppsolver()->getNoCells();
  for(MInt cellId = 0; cellId < noCells; cellId++) {
    for(MInt var = 0; var < m_noVariables; var++) {
      ppsolver()->saveVarToPrimitive(cellId, var, m_summedVars[var][cellId]);
    }
  }

  ppsolver()->exchange();

  if(ppsolver()->isMovingGrid()) {
    cout << "Moving grid to correct position!" << endl;
    ppsolver()->moveGrid(true, true);
  }

  ppsolver()->applyBoundaryCondition();

  cout << "Computing conservative variables" << endl;
  ppsolver()->computeConservativeVariables();
}

/**
 * \brief Computes skin friction of an averaged field
 *
 * \author Marian Albers
 * \date 12.12.2016
 *
 *
 **/
template <MInt nDim, class SolverType>
void StructuredPostprocessing<nDim, SolverType>::computeAverageSkinFriction() {
  ppsolver()->saveAuxData();
}


template <MInt nDim, class SolverType>
void StructuredPostprocessing<nDim, SolverType>::subtractPeriodicFluctuations() {
  if(m_movingGrid) {
    if(ppsolver()->domainId() == 0) {
      cout << "Moving grid" << endl;
    }
    ppsolver()->moveGrid(true, true);

    if(ppsolver()->domainId() == 0) {
      cout << "Subtracting mean" << endl;
    }

    vector<MFloat*> spanAvgList;

    // save summed vars as preparation for spanwise average
    for(MInt var = 0; var < m_noVariables; var++) {
      for(MInt I = 0; I < ppsolver()->getNoCells(); I++) {
        m_spanAvg[var][I] = m_summedVars[var][I];
      }
    }

    for(MInt var = 0; var < m_noVariables; var++) {
      spanAvgList.push_back(m_spanAvg[var]);
    }

    ppsolver()->spanwiseAvgZonal(spanAvgList);

    for(MInt var = 0; var < m_noVariables; var++) {
      for(MInt I = 0; I < ppsolver()->getNoCells(); I++) {
        const MFloat varMean = m_summedVars[var][I];
        const MFloat varSpanAvg = m_spanAvg[var][I];
        const MFloat varPerFluc = varMean - varSpanAvg;
        const MFloat varMeanWOPerFluc = varMean - varPerFluc;
        ppsolver()->saveVarToPrimitive(I, var, varMeanWOPerFluc);
      }
    }
  } else {
    MFloatScratchSpace cellVars(m_noVariables, AT_, "cellVars");

    for(MInt I = 0; I < ppsolver()->getNoCells(); I++) {
      ppsolver()->getSampleVariables(I, cellVars.begin());
      for(MInt var = 0; var < m_noVariables; var++) {
        const MFloat varMean = m_summedVars[var][I];
        ppsolver()->saveVarToPrimitive(I, var, varMean);
      }
    }
  }

  ppsolver()->exchange();
  ppsolver()->applyBoundaryCondition();
  ppsolver()->computeConservativeVariables();
}


template <MInt nDim, class SolverType>
void StructuredPostprocessing<nDim, SolverType>::subtractMean() {
  if(m_movingGrid) {
    if(ppsolver()->domainId() == 0) {
      cout << "Moving grid" << endl;
    }
    ppsolver()->moveGrid(true, true);
    if(ppsolver()->domainId() == 0) {
      cout << "Reordering cells" << endl;
    }
    ppsolver()->spanwiseWaveReorder();
    if(ppsolver()->domainId() == 0) {
      cout << "Subtracting mean" << endl;
    }

    for(MInt var = 0; var < m_noVariables; var++) {
      for(MInt I = 0; I < ppsolver()->getNoCells(); I++) {
        const MFloat varMean = m_summedVars[var][I];
        const MFloat varInst = m_tempWaveSample[var][I];
        const MFloat varFluc = varInst - varMean;
        ppsolver()->saveVarToPrimitive(I, var, varFluc);
      }
    }
  } else {
    MFloatScratchSpace cellVars(m_noVariables, AT_, "cellVars");

    if(ppsolver()->domainId() == 0) {
      cout << "Subtracting mean" << endl;
    }
    for(MInt I = 0; I < ppsolver()->getNoCells(); I++) {
      ppsolver()->getSampleVariables(I, cellVars.begin());
      for(MInt var = 0; var < m_noVariables; var++) {
        const MFloat varMean = m_summedVars[var][I];
        const MFloat varInst = cellVars[var];
        const MFloat varFluc = varInst - varMean;
        ppsolver()->saveVarToPrimitive(I, var, varFluc);
      }
    }
  }

  ppsolver()->exchange();
  ppsolver()->applyBoundaryCondition();
  ppsolver()->computeConservativeVariables();

  ppsolver()->computeLambda2Criterion();
  ppsolver()->saveBoxes();
}


/**
 * \brief Computes the production terms from an averaged field
 *
 * \author Marian Albers
 * \date 10.01.2017
 *
 *
 **/
template <MInt nDim, class SolverType>
void StructuredPostprocessing<nDim, SolverType>::computeProductionTerms() {
  TRACE();

  if(ppsolver()->isMovingGrid()) {
    cout << "Moving grid to correct position!" << endl;
    ppsolver()->moveGrid(true, true);
  }

  const MInt noAveragedVorticities = (m_averageVorticity != 0) * (nDim * 2 - 3);
  const MInt offset = m_noVariables + noAveragedVorticities;

  MFloat* ubar = &m_summedVars[0][0];
  MFloat* vbar = &m_summedVars[1][0];
  MFloat* wbar = &m_summedVars[2][0];

  MFloat* uu = &m_summedVars[offset + 0][0];
  MFloat* vv = &m_summedVars[offset + 1][0];
  MFloat* ww = &m_summedVars[offset + 2][0];
  MFloat* uv = &m_summedVars[offset + 3][0];
  MFloat* vw = &m_summedVars[offset + 4][0];
  MFloat* uw = &m_summedVars[offset + 5][0];

  for(MInt cellId = 0; cellId < ppsolver()->getNoCells(); cellId++) {
    m_production[0][cellId] = -uu[cellId] * ppsolver()->dvardxyz(cellId, 0, ubar)
                              - uv[cellId] * ppsolver()->dvardxyz(cellId, 1, ubar)
                              - uw[cellId] * ppsolver()->dvardxyz(cellId, 2, ubar);
    m_production[1][cellId] = -uv[cellId] * ppsolver()->dvardxyz(cellId, 0, vbar)
                              - vv[cellId] * ppsolver()->dvardxyz(cellId, 1, vbar)
                              - vw[cellId] * ppsolver()->dvardxyz(cellId, 2, vbar);
    m_production[2][cellId] = -uw[cellId] * ppsolver()->dvardxyz(cellId, 0, wbar)
                              - vw[cellId] * ppsolver()->dvardxyz(cellId, 1, wbar)
                              - ww[cellId] * ppsolver()->dvardxyz(cellId, 2, wbar);
  }

  ppsolver()->saveProductionTerms(m_postprocessFileName.c_str(), m_production);
}

/**
 * \brief Computes the production terms from an averaged field
 *
 * \author Marian Albers
 * \date 10.01.2017
 *
 *
 **/
template <MInt nDim, class SolverType>
void StructuredPostprocessing<nDim, SolverType>::computeDissipationTerms() {
  TRACE();

  if(ppsolver()->isMovingGrid()) {
    cout << "Moving grid to correct position!" << endl;
    ppsolver()->moveGrid(true, true);
  }

  MInt noFiles = (m_dissFileEnd - m_dissFileStart) / m_dissFileStep;
  MInt noSamples = 0;

  if(ppsolver()->domainId() == 0) {
    cout << "Computing dissipation..." << endl;
  }

  MFloatScratchSpace diss1(ppsolver()->getNoCells(), AT_, "diss1");
  MFloatScratchSpace diss2(ppsolver()->getNoCells(), AT_, "diss2");

  for(MInt n = 0; n < noFiles; n++) {
    MInt currentStep = m_dissFileStart + n * m_dissFileStep;
    MBool result = ppsolver()->loadBoxFile(m_dissFileDir, m_dissFilePrefix, currentStep, m_dissFileBoxNr);

    if(result == false) {
      continue;
    }
    noSamples++;

    if(ppsolver()->isMovingGrid()) {
      ppsolver()->spanwiseWaveReorder();
      if(ppsolver()->domainId() == 0) {
        cout << "After spanwise wave reorder" << endl;
      }


      for(MInt I = 0; I < ppsolver()->getNoCells(); I++) {
        for(MInt var = 0; var < 5; var++) {
          ppsolver()->setPV(var, I, m_tempWaveSample[var][I]);
        }
      }
    }

    ppsolver()->exchange();
    ppsolver()->applyBoundaryCondition();

    MFloatScratchSpace velFluc(3, ppsolver()->getNoCells(), AT_, "uFluc");

    if(ppsolver()->domainId() == 0) {
      cout << "Computing fluctuations..." << endl;
    }

    for(MInt I = 0; I < ppsolver()->getNoCells(); I++) {
      for(MInt var = 0; var < 5; var++) {
        const MFloat varMean = m_summedVars[var][I];
        const MFloat varInst = ppsolver()->getPV(var, I);
        const MFloat fluc = varInst - varMean;
        ppsolver()->setPV(var, I, fluc);
      }
    }

    ppsolver()->exchange();
    ppsolver()->applyBoundaryCondition();

    for(MInt I = 0; I < ppsolver()->getNoCells(); I++) {
      for(MInt var = 0; var < 3; var++) {
        velFluc(var, I) = ppsolver()->getPV(var, I);
      }
    }

    if(ppsolver()->domainId() == 0) {
      cout << "Computing strain tensor..." << endl;
    }

    for(MInt I = 0; I < ppsolver()->getNoCells(); I++) {
      for(MInt ii = 0; ii < 3; ii++) {
        for(MInt jj = 0; jj < 3; jj++) {
          const MFloat sij = ppsolver()->dvardxyz(I, jj, &velFluc(ii, 0));
          const MFloat sji = ppsolver()->dvardxyz(I, ii, &velFluc(jj, 0));

          diss1[I] += sij * sij;
          diss2[I] += sij * sji;
        }
      }
    }
  }

  for(MInt I = 0; I < ppsolver()->getNoCells(); I++) {
    diss1[I] /= noSamples;
    diss2[I] /= noSamples;

    m_dissipation[I] = diss1[I] + diss2[I];
  }

  if(ppsolver()->domainId() == 0) {
    cout << "Computing dissipation..." << endl;
  }

  ppsolver()->saveDissipation(m_postprocessFileName.c_str(), m_dissipation);
}

template <MInt nDim, class SolverType>
void StructuredPostprocessing<nDim, SolverType>::movingAverage() {
  TRACE();
}

template <MInt nDim, class SolverType>
void StructuredPostprocessing<nDim, SolverType>::movingAveragePost() {
  TRACE();
}


template <MInt nDim, class SolverType>
void StructuredPostprocessing<nDim, SolverType>::writeGradients() {
  TRACE();
  if(m_movingGrid) {
    if(ppsolver()->domainId() == 0) {
      cout << "Moving grid" << endl;
    }
    ppsolver()->moveGrid(true, true);
  }

  const MInt noVars = 9;
  const MInt noCells = ppsolver()->getNoCells();
  const MInt noAveragedVorticities = (m_averageVorticity != 0) * (nDim * 2 - 3);
  const MInt offset = m_noVariables + noAveragedVorticities;


  // spanwise average for all
  if(!m_movingGrid) {
    vector<MFloat*> spanAvgList;
    for(MInt var = 0; var < getNoPPVars(); var++) {
      spanAvgList.push_back(m_summedVars[var]);
    }

    ppsolver()->spanwiseAvgZonal(spanAvgList);
    vector<MFloat*> ppVariables;
    for(MInt var = 0; var < getNoPPVars(); var++) {
      ppVariables.push_back(m_summedVars[var]);
    }
    ppsolver()->gcFillGhostCells(ppVariables);
  }


  // first for u,v,w
  for(MInt dim = 0; dim < nDim; dim++) {
    for(MInt var = 0; var < 3; var++) {
      for(MInt cellId = 0; cellId < noCells; cellId++) {
        m_gradients[var + dim * noVars][cellId] = ppsolver()->dvardxyz(cellId, dim, m_summedVars[var]);
      }
    }

    // then for uu,vv,ww,uv,uw,vw
    for(MInt var = 0; var < 6; var++) {
      for(MInt cellId = 0; cellId < noCells; cellId++) {
        m_gradients[3 + var + dim * noVars][cellId] = ppsolver()->dvardxyz(cellId, dim, m_summedVars[var + offset]);
      }
    }
  }

  MStringScratchSpace gradientNames(27, AT_, "gradientNames");

  gradientNames[0] = "dudx";
  gradientNames[1] = "dvdx";
  gradientNames[2] = "dwdx";
  gradientNames[3] = "duudx";
  gradientNames[4] = "dvvdx";
  gradientNames[5] = "dwwdx";
  gradientNames[6] = "duvdx";
  gradientNames[7] = "duwdx";
  gradientNames[8] = "dvwdx";

  gradientNames[9] = "dudy";
  gradientNames[10] = "dvdy";
  gradientNames[11] = "dwdy";
  gradientNames[12] = "duudy";
  gradientNames[13] = "dvvdy";
  gradientNames[14] = "dwwdy";
  gradientNames[15] = "duvdy";
  gradientNames[16] = "duwdy";
  gradientNames[17] = "dvwdy";

  gradientNames[18] = "dudz";
  gradientNames[19] = "dvdz";
  gradientNames[20] = "dwdz";
  gradientNames[21] = "duudz";
  gradientNames[22] = "dvvdz";
  gradientNames[23] = "dwwdz";
  gradientNames[24] = "duvdz";
  gradientNames[25] = "duwdz";
  gradientNames[26] = "dvwdz";

  MString gradientFileName = "meanGradients.hdf5";
  ppsolver()->saveGradients(gradientFileName.c_str(), m_gradients, &gradientNames[0]);
}

template <MInt nDim, class SolverType>
void StructuredPostprocessing<nDim, SolverType>::decomposeCf() {
  TRACE();
  if(m_movingGrid) {
    if(ppsolver()->domainId() == 0) {
      cout << "Moving grid" << endl;
    }
    ppsolver()->moveGrid(true, true);
  }

  const MInt noCells = ppsolver()->getNoCells();
  const MInt noAveragedVorticities = (m_averageVorticity != 0) * (nDim * 2 - 3);
  const MInt offset = m_noVariables + noAveragedVorticities;

  MFloatScratchSpace uTilde(noCells, AT_, "uTilde");
  MFloatScratchSpace vTilde(noCells, AT_, "vTilde");
  MFloatScratchSpace wTilde(noCells, AT_, "wTilde");
  MFloatScratchSpace rhoTilde(noCells, AT_, "rhoTilde");
  MFloatScratchSpace pTilde(noCells, AT_, "pTilde");

  MFloatScratchSpace uuTilde(noCells, AT_, "uuTilde");
  MFloatScratchSpace uvTilde(noCells, AT_, "uvTilde");
  MFloatScratchSpace uwTilde(noCells, AT_, "uwTilde");
  MFloatScratchSpace mue(noCells, AT_, "mue");

  MFloat* const u = &m_summedVars[0][0];
  MFloat* const v = &m_summedVars[1][0];
  MFloat* const w = &m_summedVars[2][0];
  MFloat* const rho = &m_summedVars[3][0];
  MFloat* const p = &m_summedVars[4][0];
  MFloat* const uu = &m_summedVars[offset][0];
  MFloat* const uv = &m_summedVars[offset + 3][0];
  MFloat* const uw = &m_summedVars[offset + 4][0];

  m_sutherlandConstant = ppsolver()->getSutherlandConstant();
  m_sutherlandPlusOne = m_sutherlandConstant + F1;

  for(MInt I = 0; I < ppsolver()->getNoCells(); I++) {
    const MFloat T = ppsolver()->getGamma() * p[I] / rho[I];
    mue[I] = SUTHERLANDLAW(T);
  }

  uuTilde.fill(F0);
  uvTilde.fill(F0);
  uwTilde.fill(F0);

  if(ppsolver()->domainId() == 0) {
    cout << "Computing spanwise average/periodic fluctuations" << endl;
  }

  if(m_movingGrid) {
    MFloatScratchSpace spanAvg(5, noCells, AT_, "spanAvg");

    // for moving grid compute spanwise average to compute periodic fluctuations
    for(MInt var = 0; var < m_noVariables; var++) {
      for(MInt I = 0; I < ppsolver()->getNoCells(); I++) {
        spanAvg(var, I) = m_summedVars[var][I];
      }
    }

    vector<MFloat*> spanVariables;
    for(MInt var = 0; var < m_noVariables; var++) {
      spanVariables.push_back(&spanAvg(var, 0));
    }
    ppsolver()->spanwiseAvgZonal(spanVariables);
    ppsolver()->gcFillGhostCells(spanVariables);

    for(MInt I = 0; I < ppsolver()->getNoCells(); I++) {
      uTilde[I] = u[I] - spanAvg(0, I);
      vTilde[I] = v[I] - spanAvg(1, I);
      wTilde[I] = w[I] - spanAvg(2, I);
      rhoTilde[I] = rho[I] - spanAvg(3, I);
      pTilde[I] = p[I] - spanAvg(4, I);

      uuTilde[I] = uTilde[I] * uTilde[I];
      uvTilde[I] = uTilde[I] * vTilde[I];
      uwTilde[I] = uTilde[I] * wTilde[I];

      u[I] = spanAvg(0, I);
      v[I] = spanAvg(1, I);
      w[I] = spanAvg(2, I);
      rho[I] = spanAvg(3, I);
      p[I] = spanAvg(4, I);
    }
  } else {
    // for reference case only compute the spanwise average of the summed vars
    vector<MFloat*> ppVariables;
    for(MInt var = 0; var < getNoPPVars(); var++) {
      ppVariables.push_back(m_summedVars[var]);
    }
    ppsolver()->spanwiseAvgZonal(ppVariables);
    ppsolver()->gcFillGhostCells(ppVariables);
  }

  MFloatScratchSpace dudx(noCells, AT_, "dudx");
  MFloatScratchSpace dudy(noCells, AT_, "dudy");
  MFloatScratchSpace dudz(noCells, AT_, "dudz");
  MFloatScratchSpace dpdx(noCells, AT_, "dpdx");
  MFloatScratchSpace duTildedx(noCells, AT_, "duTildedx");
  MFloatScratchSpace duTildedy(noCells, AT_, "duTildedy");
  MFloatScratchSpace duTildedz(noCells, AT_, "duTildedz");
  MFloatScratchSpace dpTildedx(noCells, AT_, "dpTildedx");
  MFloatScratchSpace duuTildedx(noCells, AT_, "duuTildedx");
  MFloatScratchSpace duwTildedz(noCells, AT_, "duuTildedz");
  MFloatScratchSpace duudx(noCells, AT_, "dudx");
  MFloatScratchSpace duwdz(noCells, AT_, "dudx");

  dudx.fill(F0);
  dudy.fill(F0);
  dudz.fill(F0);
  dpdx.fill(F0);
  duTildedx.fill(F0);
  duTildedy.fill(F0);
  duTildedz.fill(F0);
  dpTildedx.fill(F0);
  duuTildedx.fill(F0);
  duwTildedz.fill(F0);
  duudx.fill(F0);
  duwdz.fill(F0);

  MFloatScratchSpace nududx(noCells, AT_, "nududx");
  MFloatScratchSpace nududz(noCells, AT_, "nududz");
  MFloatScratchSpace nududxx(noCells, AT_, "nududxx");
  MFloatScratchSpace nududzz(noCells, AT_, "nududzz");

  nududx.fill(F0);
  nududz.fill(F0);
  nududxx.fill(F0);
  nududzz.fill(F0);

  MFloatScratchSpace nuduTildedx(noCells, AT_, "nuduTildedx");
  MFloatScratchSpace nuduTildedz(noCells, AT_, "nuduTildedz");
  MFloatScratchSpace nuduTildedxx(noCells, AT_, "nuduTildedxx");
  MFloatScratchSpace nuduTildedzz(noCells, AT_, "nuduTildedzz");

  nuduTildedx.fill(F0);
  nuduTildedz.fill(F0);
  nuduTildedxx.fill(F0);
  nuduTildedzz.fill(F0);

  if(ppsolver()->domainId() == 0) {
    cout << "Computing gradients" << endl;
  }

  const MInt noGC = ppsolver()->getNoGhostLayers();

  const MInt* nCells = ppsolver()->getCellGrid();
  const MInt* nActiveCells = ppsolver()->getActiveCells();
  const MInt* nOffsetCells = ppsolver()->getOffsetCells();

  for(MInt i = 1; i < nCells[2] - 1; i++) {
    for(MInt k = 1; k < nCells[0] - 1; k++) {
      for(MInt j = 1; j < nCells[1] - 1; j++) {
        const MInt I = i + (j + k * nCells[1]) * nCells[2];
        dudx[I] = ppsolver()->dvardxyz(I, 0, u);
        dudy[I] = ppsolver()->dvardxyz(I, 1, u);
        dudz[I] = ppsolver()->dvardxyz(I, 2, u);
        dpdx[I] = ppsolver()->dvardxyz(I, 0, p);

        duTildedx[I] = ppsolver()->dvardxyz(I, 0, &uTilde[0]);
        duTildedy[I] = ppsolver()->dvardxyz(I, 1, &uTilde[0]);
        duTildedz[I] = ppsolver()->dvardxyz(I, 2, &uTilde[0]);
        dpTildedx[I] = ppsolver()->dvardxyz(I, 0, &pTilde[0]);

        duudx[I] = ppsolver()->dvardxyz(I, 0, uu);
        duwdz[I] = ppsolver()->dvardxyz(I, 2, uw);

        duuTildedx[I] = ppsolver()->dvardxyz(I, 0, &uuTilde[0]);
        duwTildedz[I] = ppsolver()->dvardxyz(I, 2, &uwTilde[0]);

        nududx[I] = mue[I] / rho[I] * dudx[I];
        nududz[I] = mue[I] / rho[I] * dudz[I];

        nuduTildedx[I] = mue[I] / rho[I] * duTildedx[I];
        nuduTildedz[I] = mue[I] / rho[I] * duTildedz[I];
      }
    }
  }

  if(ppsolver()->domainId() == 0) {
    cout << "Computing second order gradients" << endl;
  }

  for(MInt i = 1; i < nCells[2] - 1; i++) {
    for(MInt k = 1; k < nCells[0] - 1; k++) {
      for(MInt j = 1; j < nCells[1] - 1; j++) {
        const MInt I = i + (j + k * nCells[1]) * nCells[2];
        nududxx[I] = ppsolver()->dvardxyz(I, 0, &nududx[0]);
        nududzz[I] = ppsolver()->dvardxyz(I, 2, &nududz[0]);

        nuduTildedxx[I] = ppsolver()->dvardxyz(I, 0, &nuduTildedx[0]);
        nuduTildedzz[I] = ppsolver()->dvardxyz(I, 2, &nuduTildedz[0]);
      }
    }
  }

  const MInt globalNoCellsI = ppsolver()->getGrid()->getMyBlockNoCells(2);
  const MInt globalNoCellsK = ppsolver()->getGrid()->getMyBlockNoCells(0);
  const MInt totalNoCellsIK = globalNoCellsI * globalNoCellsK;

  const MInt noDecompVars = 11;
  MFloatScratchSpace cfDecomposedLocal(noDecompVars, totalNoCellsIK, AT_, "cfDecomposedLocal");
  MFloatScratchSpace cfDecomposedGlobal(noDecompVars, totalNoCellsIK, AT_, "cfDecomposedGlobal");

  cfDecomposedLocal.fill(F0);

  const MFloat gammaMinusOne = ppsolver()->getGamma() - 1.0;
  const MFloat t8 = 1.0 / (1.0 + F1B2 * gammaMinusOne * POW2(ppsolver()->getMa()));
  const MFloat u8 = ppsolver()->getMa() * sqrt(t8);

  const MFloat fac = 2.0 / (POW3(u8));
  const MFloat fre0 = 1.0 / ppsolver()->getRe0();

  if(ppsolver()->domainId() == 0) {
    cout << "Computing decomposition activeCells[0]: " << nActiveCells[0] << " activeCells[1]: " << nActiveCells[1]
         << " activeCells[2]: " << nActiveCells[2] << endl;
    cout << "GlobalCellsI: " << globalNoCellsI << endl;
  }

  for(MInt i = 0; i < nActiveCells[2]; i++) {
    for(MInt k = 0; k < nActiveCells[0]; k++) {
      for(MInt j = 0; j < nActiveCells[1]; j++) {
        const MInt globalId = (i + nOffsetCells[2]) + (k + nOffsetCells[1]) * (globalNoCellsI);
        const MInt I = i + noGC + ((j + noGC) + (k + noGC) * nCells[1]) * nCells[2];

        const MFloat dy = ppsolver()->getCellLengthY(i + noGC, j + noGC, k + noGC);

        cfDecomposedLocal(0, globalId) += fac * fre0 * mue[I] / rho[I] * dudy[I] * dudy[I] * dy;

        cfDecomposedLocal(1, globalId) += fac * fre0 * mue[I] / rho[I] * dudy[I] * duTildedy[I] * dy;

        cfDecomposedLocal(2, globalId) += fac * (-uv[I] * dudy[I] * dy);

        cfDecomposedLocal(3, globalId) += fac * (-uvTilde[I] * dudy[I] * dy);

        cfDecomposedLocal(4, globalId) += fac * ((u[I] - u8) * (u[I] * dudx[I] + v[I] * dudy[I] + w[I] * dudz[I]) * dy);

        cfDecomposedLocal(5, globalId) +=
            fac * ((u[I] - u8) * (u[I] * duTildedx[I] + v[I] * duTildedy[I] + w[I] * duTildedz[I]) * dy);

        cfDecomposedLocal(6, globalId) +=
            fac * ((u[I] - u8) * (uTilde[I] * dudx[I] + vTilde[I] * dudy[I] + wTilde[I] * dudz[I]) * dy);

        cfDecomposedLocal(7, globalId) += fac * (u[I] - u8) * dpdx[I] * dy;

        cfDecomposedLocal(8, globalId) += fac * (u[I] - u8) * dpTildedx[I] * dy;

        cfDecomposedLocal(9, globalId) -=
            fac * (u[I] - u8) * (fre0 * nududxx[I] + fre0 * nuduTildedxx[I] - duuTildedx[I] - duudx[I]) * dy;
        cfDecomposedLocal(10, globalId) -=
            fac * (u[I] - u8) * (fre0 * nududzz[I] + fre0 * nuduTildedzz[I] - duwTildedz[I] - duwdz[I]) * dy;
      }
    }
  }

  if(ppsolver()->domainId() == 0) {
    cout << "Before MPI Allreduce" << endl;
  }

  MPI_Allreduce(&cfDecomposedLocal(0, 0),
                &cfDecomposedGlobal(0, 0),
                noDecompVars * totalNoCellsIK,
                MPI_DOUBLE,
                MPI_SUM,
                ppsolver()->getCommunicator(),
                AT_,
                "cfDecomposedLocal",
                "cfDecomposedGlobal");

  MFloatScratchSpace cfDecomposedLine(noDecompVars, globalNoCellsI, AT_, "cfDecomposedLine");
  cfDecomposedLine.fill(F0);

  if(ppsolver()->domainId() == 0) {
    cout << "Now averaging in spanwise direction, no cells spanwise: " << globalNoCellsK << endl;
  }

  for(MInt var = 0; var < noDecompVars; var++) {
    for(MInt i = 0; i < globalNoCellsI; i++) {
      for(MInt k = 0; k < globalNoCellsK; k++) {
        const MInt globalId = i + k * globalNoCellsI;
        cfDecomposedLine(var, i) += cfDecomposedGlobal(var, globalId);
      }
      cfDecomposedLine(var, i) /= (MFloat)globalNoCellsK;
    }
  }

  // now write to file
  if(ppsolver()->domainId() == 0) {
    MString filename = "./cf_decomposed.dat";
    FILE* f_forces;
    f_forces = fopen(filename.c_str(), "w");
    for(MInt i = 0; i < globalNoCellsI; i++) {
      for(MInt var = 0; var < noDecompVars; var++) {
        fprintf(f_forces, "%f ", cfDecomposedLine(var, i));
      }
      fprintf(f_forces, "\n");
    }
    fclose(f_forces);
  }

  if(ppsolver()->domainId() == 0) {
    cout << "Cf decomposition finished!" << endl;
  }
}


template <MInt nDim, class SolverType>
void StructuredPostprocessing<nDim, SolverType>::decomposeCfDouble() {
  TRACE();
  if(m_movingGrid) {
    if(ppsolver()->domainId() == 0) {
      cout << "Moving grid" << endl;
    }
    ppsolver()->moveGrid(true, true);
  }

  const MInt noCells = ppsolver()->getNoCells();
  const MInt noAveragedVorticities = (m_averageVorticity != 0) * (nDim * 2 - 3);
  const MInt offset = m_noVariables + noAveragedVorticities;

  MFloatScratchSpace mue(noCells, AT_, "mue");

  MFloat* const u = &m_summedVars[0][0];
  MFloat* const v = &m_summedVars[1][0];
  MFloat* const w = &m_summedVars[2][0];
  MFloat* const rho = &m_summedVars[3][0];
  MFloat* const p = &m_summedVars[4][0];
  MFloat* const uu = &m_summedVars[offset][0];
  MFloat* const uv = &m_summedVars[offset + 3][0];
  MFloat* const uw = &m_summedVars[offset + 4][0];

  m_sutherlandConstant = ppsolver()->getSutherlandConstant();
  m_sutherlandPlusOne = m_sutherlandConstant + F1;

  for(MInt I = 0; I < ppsolver()->getNoCells(); I++) {
    const MFloat T = ppsolver()->getGamma() * p[I] / rho[I];
    mue[I] = SUTHERLANDLAW(T);
  }

  if(!m_movingGrid) {
    // for reference case only compute the spanwise average of the summed vars
    vector<MFloat*> ppVariables;
    for(MInt var = 0; var < getNoPPVars(); var++) {
      ppVariables.push_back(m_summedVars[var]);
    }
    ppsolver()->spanwiseAvgZonal(ppVariables);
    ppsolver()->gcFillGhostCells(ppVariables);
  }

  MFloatScratchSpace dudx(noCells, AT_, "dudx");
  MFloatScratchSpace dudy(noCells, AT_, "dudy");
  MFloatScratchSpace dudz(noCells, AT_, "dudz");
  MFloatScratchSpace dpdx(noCells, AT_, "dpdx");
  MFloatScratchSpace duudx(noCells, AT_, "dudx");
  MFloatScratchSpace duwdz(noCells, AT_, "dudx");

  dudx.fill(F0);
  dudy.fill(F0);
  dudz.fill(F0);
  dpdx.fill(F0);
  duudx.fill(F0);
  duwdz.fill(F0);

  MFloatScratchSpace nududx(noCells, AT_, "nududx");
  MFloatScratchSpace nududz(noCells, AT_, "nududz");
  MFloatScratchSpace nududxx(noCells, AT_, "nududxx");
  MFloatScratchSpace nududzz(noCells, AT_, "nududzz");

  nududx.fill(F0);
  nududz.fill(F0);
  nududxx.fill(F0);
  nududzz.fill(F0);

  if(ppsolver()->domainId() == 0) {
    cout << "Computing gradients" << endl;
  }

  const MInt noGC = ppsolver()->getNoGhostLayers();

  const MInt* nCells = ppsolver()->getCellGrid();
  const MInt* nActiveCells = ppsolver()->getActiveCells();
  const MInt* nOffsetCells = ppsolver()->getOffsetCells();

  for(MInt i = 1; i < nCells[2] - 1; i++) {
    for(MInt k = 1; k < nCells[0] - 1; k++) {
      for(MInt j = 1; j < nCells[1] - 1; j++) {
        const MInt I = i + (j + k * nCells[1]) * nCells[2];
        dudx[I] = ppsolver()->dvardxyz(I, 0, u);
        dudy[I] = ppsolver()->dvardxyz(I, 1, u);
        dudz[I] = ppsolver()->dvardxyz(I, 2, u);
        dpdx[I] = ppsolver()->dvardxyz(I, 0, p);

        duudx[I] = ppsolver()->dvardxyz(I, 0, uu);
        duwdz[I] = ppsolver()->dvardxyz(I, 2, uw);

        nududx[I] = mue[I] / rho[I] * dudx[I];
        nududz[I] = mue[I] / rho[I] * dudz[I];
      }
    }
  }

  if(ppsolver()->domainId() == 0) {
    cout << "Computing second order gradients" << endl;
  }

  for(MInt i = 1; i < nCells[2] - 1; i++) {
    for(MInt k = 1; k < nCells[0] - 1; k++) {
      for(MInt j = 1; j < nCells[1] - 1; j++) {
        const MInt I = i + (j + k * nCells[1]) * nCells[2];
        nududxx[I] = ppsolver()->dvardxyz(I, 0, &nududx[0]);
        nududzz[I] = ppsolver()->dvardxyz(I, 2, &nududz[0]);
      }
    }
  }

  const MInt globalNoCellsI = ppsolver()->getGrid()->getMyBlockNoCells(2);
  const MInt globalNoCellsJ = ppsolver()->getGrid()->getMyBlockNoCells(1);
  const MInt globalNoCellsK = ppsolver()->getGrid()->getMyBlockNoCells(0);
  const MInt totalNoCellsIK = globalNoCellsI * globalNoCellsK;

  const MInt noDecompVars = 7;
  MFloatScratchSpace cfDecomposedLocal(noDecompVars, totalNoCellsIK, AT_, "cfDecomposedLocal");
  MFloatScratchSpace cfDecomposedGlobal(noDecompVars, totalNoCellsIK, AT_, "cfDecomposedGlobal");

  cfDecomposedLocal.fill(F0);

  const MFloat gammaMinusOne = ppsolver()->getGamma() - 1.0;
  const MFloat t8 = 1.0 / (1.0 + F1B2 * gammaMinusOne * POW2(ppsolver()->getMa()));
  const MFloat u8 = ppsolver()->getMa() * sqrt(t8);

  const MFloat fac = 2.0 / (POW3(u8));
  const MFloat fre0 = 1.0 / ppsolver()->getRe0();

  if(ppsolver()->domainId() == 0) {
    cout << "Computing decomposition activeCells[0]: " << nActiveCells[0] << " activeCells[1]: " << nActiveCells[1]
         << " activeCells[2]: " << nActiveCells[2] << endl;
    cout << "GlobalCellsI: " << globalNoCellsI << endl;
  }

  for(MInt i = 0; i < nActiveCells[2]; i++) {
    for(MInt k = 0; k < nActiveCells[0]; k++) {
      for(MInt j = 0; j < nActiveCells[1]; j++) {
        const MInt I = i + noGC + ((j + noGC) + (k + noGC) * nCells[1]) * nCells[2];
        const MInt globalId = (i + nOffsetCells[2]) + (k + nOffsetCells[1]) * (globalNoCellsI);
        const MFloat dy = ppsolver()->getCellLengthY(i + noGC, j + noGC, k + noGC);

        cfDecomposedLocal(0, globalId) += ppsolver()->getCellCoordinate(I, 0);

        cfDecomposedLocal(1, globalId) += fac * fre0 * mue[I] / rho[I] * dudy[I] * dudy[I] * dy;

        cfDecomposedLocal(2, globalId) += fac * (-uv[I] * dudy[I] * dy);

        cfDecomposedLocal(3, globalId) += fac * ((u[I] - u8) * (u[I] * dudx[I] + v[I] * dudy[I] + w[I] * dudz[I]) * dy);

        cfDecomposedLocal(4, globalId) += fac * (u[I] - u8) * dpdx[I] * dy;

        cfDecomposedLocal(5, globalId) -= fac * (u[I] - u8) * (fre0 * nududxx[I] - duudx[I]) * dy;
        cfDecomposedLocal(6, globalId) -= fac * (u[I] - u8) * (fre0 * nududzz[I] - duwdz[I]) * dy;
      }
    }
  }

  if(ppsolver()->domainId() == 0) {
    cout << "Before MPI Allreduce" << endl;
  }

  MPI_Allreduce(&cfDecomposedLocal(0, 0),
                &cfDecomposedGlobal(0, 0),
                noDecompVars * totalNoCellsIK,
                MPI_DOUBLE,
                MPI_SUM,
                ppsolver()->getCommunicator(),
                AT_,
                "cfDecomposedLocal",
                "cfDecomposedGlobal");

  MFloatScratchSpace cfDecomposedLine(noDecompVars, globalNoCellsI, AT_, "cfDecomposedLine");
  cfDecomposedLine.fill(F0);

  if(ppsolver()->domainId() == 0) {
    cout << "Now averaging in spanwise direction, no cells spanwise: " << globalNoCellsK << endl;
  }

  for(MInt var = 0; var < noDecompVars; var++) {
    for(MInt i = 0; i < globalNoCellsI; i++) {
      for(MInt k = 0; k < globalNoCellsK; k++) {
        const MInt globalId = i + k * globalNoCellsI;
        cfDecomposedLine(var, i) += cfDecomposedGlobal(var, globalId);
      }
      cfDecomposedLine(var, i) /= (MFloat)globalNoCellsK;

      if(var == 0) {
        cfDecomposedLine(var, i) /= (MFloat)globalNoCellsJ;
      }
    }
  }

  // now write to file
  if(ppsolver()->domainId() == 0) {
    MString filename = "./cf_decomposed.dat";
    FILE* f_forces;
    f_forces = fopen(filename.c_str(), "w");
    for(MInt i = 0; i < globalNoCellsI; i++) {
      for(MInt var = 0; var < noDecompVars; var++) {
        fprintf(f_forces, "%f ", cfDecomposedLine(var, i));
      }
      fprintf(f_forces, "\n");
    }
    fclose(f_forces);
  }

  if(ppsolver()->domainId() == 0) {
    cout << "Cf decomposition finished!" << endl;
  }
}

/**
 * \brief Initializes properties for averaging during solver run
 *
 * \author Andreas Lintermann (last modified Ansgar Niemoeller 07/14)
 * \date 14.09.2012
 *
 * \tparam[in] SolverType solvertype
 *
 **/
template <MInt nDim, class SolverType>
void StructuredPostprocessing<nDim, SolverType>::initAverageIn() {
  TRACE();

  initTimeStepProperties();

  /*! \property
    \page propertiesFVSTRCTRD
    \section pp_averageRestart
    <code>MString* Solver::m_averageRestart</code>\n
    default = <code>0</code>\n\n
    This property determines if we should restart from our last averaging.
    <ul>
    <li><code>0</code> turned off </li>
    <li><code>1</code> turned on</li>
    </ul>\n
    Keywords: <i>GENERAL, GLOBAL, POSTPROCESSING</i>
  */
  m_averageRestart =
      Context::getSolverProperty<MInt>("pp_averageRestart", ppsolver()->solverId(), AT_, &m_averageRestart);

  /*! \property
    \page propertiesFVSTRCTRD
    \section pp_averageRestartInterval
    <code>MString* Solver::m_averageRestartInterval</code>\n
    default = <code>0</code>\n\n
    This property determines the interval to write averaging restart files. Has to be a multiple of averageInterval and
    of restartInterval. <ul> <li><code>interval</code> </li>
    </ul>\n
    Keywords: <i>GENERAL, GLOBAL, POSTPROCESSING</i>
  */
  m_averageRestartInterval = Context::getSolverProperty<MInt>(
      "pp_averageRestartInterval", ppsolver()->solverId(), AT_, &m_averageRestartInterval);

  /*! \property
    \page propertiesFVSTRCTRD
    \section pp_averagingFavre
    <code>MInt PostprocesssingSolver::m_averagingFavre</code>\n
    default = <code>0</code>\n\n
    Computes additional density-correlated averages
    <ul>c
    <li><code>0</code> disabled</li>
    <li><code>1</code> enabled</li>
    </ul>\n
    Keywords: <i>GLOBAL, POSTPROCESSING</i>
  */
  m_averagingFavre = false;
  m_averagingFavre =
      Context::getSolverProperty<MBool>("pp_averagingFavre", ppsolver()->solverId(), AT_, &m_averagingFavre);


  if(m_averageRestartInterval % ppsolver()->restartInterval() != 0) {
    mTerm(1, AT_, "The property 'averageRestartInterval' has to be a multiple of the property 'restartInterval'...");
  }
}

/**
 * \brief allocates memory for averageSolutions() and averageSolutionsInSolve()
 *
 * \author A. Niemoeller
 * \date 09.12.2013
 *
 * \param[in] noInternalCells the number of internal cells on which averaging is performed
 *
 **/
template <MInt nDim, class SolverType>
void StructuredPostprocessing<nDim, SolverType>::initAverageVariables() {
  TRACE();
  const MInt noCells = ppsolver()->getNoCells();
  const MInt noVars = getNoPPVars();
  const MInt noSquareVars = getNoPPSquareVars();

  mAlloc(m_summedVars, noVars, noCells, "m_summedVars", F0, FUN_); // +1 for pressure ampl
  mAlloc(m_square, noSquareVars, noCells, "m_square", F0, FUN_);   // + 1 for pressure ampl

  if(m_averagingFavre) {
    mAlloc(m_favre, getNoVars(), noCells, "m_favre", F0, FUN_); // for Favre averaging
  }

  if(m_kurtosis) {
    m_log << "Allocating cube and fourth field for kurtosis computation for " << noCells << " cells" << endl;
    mAlloc(m_cube, nDim, noCells, "m_cube", F0, FUN_);
    mAlloc(m_fourth, nDim, noCells, "m_fourth", F0, FUN_);
  } else if(m_skewness /*&& !m_twoPass*/) {
    mAlloc(m_cube, nDim, noCells, "m_cube", F0, FUN_);
  }

  if(m_useKahan) // allocate memory for kahan summation
  {
    m_log << "m_useKahan is activated" << endl;
    mAlloc(m_cSum, m_noVariables + 3 * (nDim - 1), noCells, "m_cSum", F0, FUN_);
    mAlloc(m_tSum, m_noVariables + 3 * (nDim - 1), noCells, "m_tSum", F0, FUN_);
    mAlloc(m_ySum, m_noVariables + 3 * (nDim - 1), noCells, "m_ySum", F0, FUN_);
    mAlloc(m_cSquare, 3 * (nDim - 1), noCells, "m_cSquare", F0, FUN_);
    mAlloc(m_tSquare, 3 * (nDim - 1), noCells, "m_tSquare", F0, FUN_);
    mAlloc(m_ySquare, 3 * (nDim - 1), noCells, "m_ySquare", F0, FUN_);
    if(m_kurtosis) {
      mAlloc(m_cCube, nDim, noCells, "m_cCube", F0, FUN_);
      mAlloc(m_tCube, nDim, noCells, "m_tCube", F0, FUN_);
      mAlloc(m_yCube, nDim, noCells, "m_yCube", F0, FUN_);
      mAlloc(m_cFourth, nDim, noCells, "m_cFourth", F0, FUN_);
      mAlloc(m_tFourth, nDim, noCells, "m_tFourth", F0, FUN_);
      mAlloc(m_yFourth, nDim, noCells, "m_yFourth", F0, FUN_);
    } else if(m_skewness) {
      mAlloc(m_cCube, nDim, noCells, "m_cCube", F0, FUN_);
      mAlloc(m_tCube, nDim, noCells, "m_tCube", F0, FUN_);
      mAlloc(m_yCube, nDim, noCells, "m_yCube", F0, FUN_);
    }
  }

  mAlloc(m_avgVariableNames, noVars, "m_avgVariableNames", AT_);
  mAlloc(m_avgFavreNames, getNoVars(), "m_avgFavreNames", AT_);


  // Mean values
  m_avgVariableNames[0] = "um";
  m_avgVariableNames[1] = "vm";
  IF_CONSTEXPR(nDim == 3) { m_avgVariableNames[2] = "wm"; }
  m_avgVariableNames[nDim] = "rhom";
  m_avgVariableNames[nDim + 1] = "pm";
  MInt offset = m_noVariables;

  // Vorticities
  if(m_averageVorticity) {
    IF_CONSTEXPR(nDim == 3) {
      m_avgVariableNames[offset + 0] = "vortxm";
      m_avgVariableNames[offset + 1] = "vortym";
      m_avgVariableNames[offset + 2] = "vortzm";
      offset += 3;
    }
    else {
      m_avgVariableNames[offset + 0] = "vortzm";
      offset += 1;
    }
  }

  // reynolds stress components
  m_avgVariableNames[offset + 0] = "uu";
  m_avgVariableNames[offset + 1] = "vv";
  offset += 2;
  IF_CONSTEXPR(nDim == 3) {
    m_avgVariableNames[offset + 0] = "ww";
    m_avgVariableNames[offset + 1] = "uv";
    m_avgVariableNames[offset + 2] = "vw";
    m_avgVariableNames[offset + 3] = "uw";
    offset += 4;
  }
  else {
    m_avgVariableNames[offset + 0] = "uv";
    offset += 1;
  }

  // Skewness variables
  if(noVars > m_noVariables + 3 * (nDim - 1) + 1 && m_skewness) {
    IF_CONSTEXPR(nDim == 3) {
      m_avgVariableNames[offset + 0] = "uuu";
      m_avgVariableNames[offset + 1] = "vvv";
      m_avgVariableNames[offset + 2] = "www";
      offset += 3;
    }
    else {
      m_avgVariableNames[offset + 0] = "uuu";
      m_avgVariableNames[offset + 1] = "vvv";
      offset += 2;
    }
  }

  // Kurtosis variables
  if((noVars > m_noVariables + 3 * (nDim - 1) + nDim + 1) && m_kurtosis) {
    IF_CONSTEXPR(nDim == 3) {
      m_avgVariableNames[offset + 0] = "uuuu";
      m_avgVariableNames[offset + 1] = "vvvv";
      m_avgVariableNames[offset + 2] = "wwww";
      offset += 3;
    }
    else {
      m_avgVariableNames[offset + 0] = "uuuu";
      m_avgVariableNames[offset + 1] = "vvvv";
      offset += 2;
    }
  }

  // pressure fluctuation
  m_avgVariableNames[offset + 0] = "pp";
  offset += 1;

  // rms of the vorticities
  if(m_averageVorticity) {
    IF_CONSTEXPR(nDim == 3) {
      m_avgVariableNames[offset + 0] = "vortrmsx";
      m_avgVariableNames[offset + 1] = "vortrmsy";
      m_avgVariableNames[offset + 2] = "vortrmsz";
      offset += 3;
    }
    else {
      m_avgVariableNames[offset + 0] = "vortrmsz";
      offset += 1;
    }
  }

  if(m_averagingFavre) {
    // Mean values
    m_avgFavreNames[0] = "um_favre";
    m_avgFavreNames[1] = "vm_favre";
    IF_CONSTEXPR(nDim == 3) { m_avgFavreNames[2] = "wm_favre"; }
    m_avgFavreNames[nDim] = "rhom_favre";
    m_avgFavreNames[nDim + 1] = "pm_favre";
  }
}

template <MInt nDim, class SolverType>
void StructuredPostprocessing<nDim, SolverType>::initProductionVariables() {
  const MInt noCells = ppsolver()->getNoCells();
  mAlloc(m_production, nDim, noCells, "m_production", F0, FUN_);
}

template <MInt nDim, class SolverType>
void StructuredPostprocessing<nDim, SolverType>::initDissipationVariables() {
  const MInt noCells = ppsolver()->getNoCells();
  mAlloc(m_dissipation, noCells, "m_dissipation", F0, FUN_);


  m_dissFileDir = "";
  m_dissFileDir = Context::getSolverProperty<MString>("pp_dissFileDir", ppsolver()->solverId(), AT_, &m_dissFileDir);

  m_dissFilePrefix = "";
  m_dissFilePrefix =
      Context::getSolverProperty<MString>("pp_dissFilePrefix", ppsolver()->solverId(), AT_, &m_dissFilePrefix);

  m_dissFileBoxNr = 0;
  m_dissFileBoxNr = Context::getSolverProperty<MInt>("pp_dissFileBoxNr", ppsolver()->solverId(), AT_, &m_dissFileBoxNr);

  m_dissFileStart = -1;
  m_dissFileStart = Context::getSolverProperty<MInt>("pp_dissFileStart", ppsolver()->solverId(), AT_, &m_dissFileStart);

  m_dissFileStep = -1;
  m_dissFileStep = Context::getSolverProperty<MInt>("pp_dissFileStep", ppsolver()->solverId(), AT_, &m_dissFileStep);

  m_dissFileEnd = -1;
  m_dissFileEnd = Context::getSolverProperty<MInt>("pp_dissFileEnd", ppsolver()->solverId(), AT_, &m_dissFileEnd);
}


/**
 * \brief Initializes timestep properties for postprocessing
 *
 * \author Andreas Lintermann (last modified Ansgar Niemoeller, 07/14)
 * \date 14.09.2012
 *
 * reads properties pp_averageStartTimestep, pp_averageStopTimestep and pp_averageInterval
 *
 * \tparam[in] SolverType solvertype
 *
 **/
template <MInt nDim, class SolverType>
void StructuredPostprocessing<nDim, SolverType>::initTimeStepProperties() {
  TRACE();

  /*! \property
    \page propertiesFVSTRCTRD
    \section pp_averageInterval
    <code>MString* Solver::m_averageInterval</code>\n
    default = <code>0</code>\n\n
    This property determines the interval of the solutions used for averaging.
    <ul>
    <li><code>interval</code> </li>
    </ul>\n
    Keywords: <i>GENERAL, GLOBAL, POSTPROCESSING</i>
  */
  m_averageInterval = 0;
  m_averageInterval =
      Context::getSolverProperty<MInt>("pp_averageInterval", ppsolver()->solverId(), AT_, &m_averageInterval);


  if(m_averageInterval == 0) {
    mTerm(1, AT_, "Please specify the property 'averageInterval' ...");
  }

  /*! \property
    \page propertiesFVSTRCTRD
    \section pp_averageStartTimestep
    <code>MString* Solver::m_averageStartTimestep</code>\n
    default = <code>0</code>\n\n
    This property determines the start timestep used for averaging.
    <ul>
    <li><code>timestep</code> </li>
    </ul>\n
    Keywords: <i>GENERAL, GLOBAL, POSTPROCESSING</i>
  */
  m_averageStartTimestep = -1;
  m_averageStartTimestep =
      Context::getSolverProperty<MInt>("pp_averageStartTimestep", ppsolver()->solverId(), AT_, &m_averageStartTimestep);

  if(m_averageStartTimestep == 0) {
    mTerm(1, AT_, "Please specify the property 'averageStartTimestep' ...");
  }

  /*! \property
    \page propertiesFVSTRCTRD
    \section pp_averageStopTimestep
    <code>MString* Solver::m_averageStopTimestep</code>\n
    default = <code>0</code>\n\n
    This property determines the stop timestep used for averaging.
    <ul>
    <li><code>timestep</code> </li>
    </ul>\n
    Keywords: <i>GENERAL, GLOBAL, POSTPROCESSING</i>
  */
  m_averageStopTimestep = -1;
  m_averageStopTimestep =
      Context::getSolverProperty<MInt>("pp_averageStopTimestep", ppsolver()->solverId(), AT_, &m_averageStopTimestep);

  if(m_averageStopTimestep == 0) {
    mTerm(1, AT_, "Please specify the property 'averageStopTimestep' ...");
  }

  /*! \property
    \page propertiesFVSTRCTRD
    \section pp_averageInterval
    <code>MString* Solver::m_averageRestart</code>\n
    default = <code>0</code>\n\n
    This property determines if we should restart from our last averaging.
    <ul>
    <li><code>0</code> turned off </li>
    <li><code>1</code> turned on</li>
    </ul>\n
    Keywords: <i>GENERAL, GLOBAL, POSTPROCESSING</i>
  */
  m_averageRestart = 0;
  m_averageRestart =
      Context::getSolverProperty<MInt>("pp_averageRestart", ppsolver()->solverId(), AT_, &m_averageRestart);


  /*! \property
    \page propertiesFVSTRCTRD
    \section pp_averageRestartInterval
    <code>MString* Solver::m_averageRestartInterval</code>\n
    default = <code>0</code>\n\n
    This property determines the interval to write averaging restart files. Has to be a multiple of averageInterval and
    of restartInterval. <ul> <li><code>interval</code> </li>
    </ul>\n
    Keywords: <i>GENERAL, GLOBAL, POSTPROCESSING</i>
  */
  m_averageRestartInterval = 1000000;
  m_averageRestartInterval = Context::getSolverProperty<MInt>(
      "pp_averageRestartInterval", ppsolver()->solverId(), AT_, &m_averageRestartInterval);

  /**
   * \brief initializes properties and allocates memory for moving averaging
   *
   * \author Ansgar Niemoeller
   * \date 09.07.2014
   *
   * reads properties pp_movingAvgInterval, pp_movingAvgDataPoints and pp_averageVorticity
   *
   * \param[in] grid pointer to the grid
   *
   **/
}
template <MInt nDim, class SolverType>
void StructuredPostprocessing<nDim, SolverType>::initMovingAverage() {
  TRACE();

  const MInt noCells = ppsolver()->getNoCells();

  /*! \property
    \page propertiesFVSTRCTRD
    \section pp_movingAvgInterval
    <code>MInt PostprocesssingSolver::m_movingAvgInterval</code>\n
    default = <code>1</code>\n\n
    This property determines the interval between timesteps considered for moving average,
    e.g. if set to 2 at an averaging timestep n (see pp_averagStartTimestep, pp_averageStopTimestep, pp_averageInterval)
    the timesteps n, n-2, n-4, ... are used to compute the moving average\n
    Note: this has to be a factor of m_averageInterval\n
    see also pp_movingAvgDataPoints, pp_averageVorticity
    Keywords: <i>GLOBAL, POSTPROCESSING</i>
  */
  m_movingAvgInterval = 1;
  m_movingAvgInterval =
      Context::getSolverProperty<MInt>("pp_movingAvgInterval", ppsolver()->solverId(), AT_, &m_movingAvgInterval);
  if(m_movingAvgInterval < 1) {
    TERMM(1, "m_movingAvgInterval has to be >=1");
  }
  if(m_averageInterval % m_movingAvgInterval != 0 || m_movingAvgInterval > m_averageInterval) {
    TERMM(1, "m_movingAvgInterval has to be a factor of m_averageInterval");
  }

  /*! \property
    \page propertiesFVSTRCTRD
    \section pp_movingAvgDataPoints
    <code>MInt PostprocesssingSolver::m_movingAvgDataPoints</code>\n
    This property determines the number of timesteps (data points) used for moving average computation\n
    see also pp_movingAvgInterval, pp_averageVorticity
    Keywords: <i>GLOBAL, POSTPROCESSING</i>
  */
  m_movingAvgDataPoints =
      Context::getSolverProperty<MInt>("pp_movingAvgDataPoints", ppsolver()->solverId(), AT_, &m_movingAvgDataPoints);
  if(m_movingAvgDataPoints < 2) {
    TERMM(1, "m_movingAvgDataPoints has to be at least 2");
  }

  /*! \property
    \page propertiesFVSTRCTRD
    \section pp_averageVorticity
    <code>MInt PostprocesssingSolver::m_pp.m_averageVorticity</code>\n
    default = <code>0</code>\n\n
    This property determines if the vorticity vector is considered in computation of averages
    <ul>c
    <li><code>0</code> disabled</li>
    <li><code>1</code> enabled</li>
    </ul>\n
    Keywords: <i>GLOBAL, POSTPROCESSING</i>
  */
  m_averageVorticity = 0;
  m_averageVorticity =
      Context::getSolverProperty<MBool>("pp_averageVorticity", ppsolver()->solverId(), AT_, &m_averageVorticity);

  m_movingAvgCounter = 0;

  m_movAvgNoVariables = m_noVariables;
  if(m_averageVorticity == 1) {
    m_movAvgNoVariables += (2 * nDim - 3);
  }
  mAlloc(m_movAvgVariables, noCells, m_movAvgNoVariables * m_movingAvgDataPoints, "m_movAvgVariables", F0, FUN_);
  mAlloc(m_movAvgVarNames, 2 * m_movAvgNoVariables, "m_movAvgVarNames", MString("default"), FUN_);
}


template <MInt nDim, class SolverType>
void StructuredPostprocessing<nDim, SolverType>::saveAverageRestart() {
  TRACE();

  if(m_postprocessing && (globalTimeStep >= m_averageStartTimestep && globalTimeStep <= m_averageStopTimestep)) {
    MString name = ppsolver()->outputDir() + "PostprocessingRestart_";
    MChar buf1[10];
    MChar buf2[10];
    sprintf(buf1, "%d", m_averageStartTimestep);
    sprintf(buf2, "%d", globalTimeStep);
    name.append(buf1);
    name += "-";
    name.append(buf2);
    m_log << "       ^     saving average restart " << name << endl;

    ppsolver()->saveAverageRestart(name, getNoPPVars(), m_summedVars, m_square, m_cube, m_fourth);
  }
}

/**
 * \brief Returns number of postprocessing variables
 *
 * \author Marian Albers
 * \date 01.08.2016
 *
 *
 **/
template <MInt nDim, class SolverType>
MInt StructuredPostprocessing<nDim, SolverType>::getNoPPVars() {
  TRACE();
  MInt c = 0;
  if(m_kurtosis)
    c = 2;
  else if(m_skewness)
    c = 1;


  // Determine number of averaged variables
  // Mean vorticities and symmetric rms components (6 for 3D, 4 for 2D)
  const MInt noAveragedVorticities = (m_averageVorticity != 0) * (nDim * 2 - 3);
  const MInt noVars = m_noVariables            // primitive variables
                      + noAveragedVorticities  // mean vorticities
                      + 3 * (nDim - 1)         // Reynolds stress components
                      + nDim * c               // skewness/kurtosis
                      + 1                      // pressure amplitude p'
                      + noAveragedVorticities; // rms vorticities

  return noVars;
}

/**
 * \brief Returns number of pp Square variables
 *
 * \author Marian Albers
 * \date 01.08.2016
 *
 *
 **/
template <MInt nDim, class SolverType>
MInt StructuredPostprocessing<nDim, SolverType>::getNoPPSquareVars() {
  TRACE();

  // Determine number of averaged variables
  const MInt noAveragedVorticities = (m_averageVorticity != 0) * (nDim * 2 - 3);
  const MInt noVars = 3 * (nDim - 1)           // uu,vv,ww,uv,uw,vw
                      + 1                      // pressure amplitude p'
                      + noAveragedVorticities; // vorticity rms
  return noVars;
}

template <MInt nDim, class SolverType>
MInt StructuredPostprocessing<nDim, SolverType>::getNoVars() {
  TRACE();
  const MInt noVars = nDim + 2;
  return noVars;
}

// Explicit instantiations for 2D and 3D
template class StructuredPostprocessing<2, FvStructuredSolver<2>>;
template class StructuredPostprocessing<3, FvStructuredSolver<3>>;
