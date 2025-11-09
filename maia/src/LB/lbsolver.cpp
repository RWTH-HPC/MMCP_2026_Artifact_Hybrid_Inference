// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "lbsolver.h"

#include <cstring>
#include <functional>
#include <unordered_map>
#include <utility>
#include "COMM/mpioverride.h"
#include "GEOM/geometryelement.h"
#include "IO/logtable.h"
#include "IO/parallelio.h"
#include "UTIL/maiamath.h"
#include "UTIL/parallelfor.h"
#include "globals.h"
#include "lbbndcnd.h"
#include "lbconstants.h"
#include "lbfunctions.h"
#include "lbgridboundarycell.h"
#include "lbinterfacecell.h"
#include "lblatticedescriptor.h"
#include "lbparentcell.h"

using namespace std::placeholders;
using namespace std;

// the swap endian directive uses changes the endian for the binary VTK output
#define SWAP_ENDIAN

using namespace lbconstants;

template <MInt nDim>
LbSolver<nDim>::LbSolver(MInt id, MInt dist, GridProxy& gridProxy_, Geometry<nDim>& geometry_, const MPI_Comm comm)
  : maia::CartesianSolver<nDim, LbSolver<nDim>>(id, gridProxy_, comm, true)

    /// [Splitt] The following is part of a first step to splitt
    /// CartesianGrid from the inheritance hierarchy:
    ///
    /// - in order to avoid renaming a lot of access to CartesianGrid data
    ///   members, references are introduced. These references must be
    ///   initialized in the solver's constructor.
    ///
    /// \todo labels:LB this references will be removed in future commits
    ,
    m_geometry(&geometry_),
    /// the references to be initialized end here

    m_noDistributions(dist) {
  TRACE();

  initTimer();
  RECORD_TIMER_START(m_t.solver);

  m_geometryIntersection = new GeometryIntersection<nDim>(&grid(), m_geometry);

  auto solverMethodEnum = string2enum(solverMethod());
  if(solverMethodEnum == MAIA_LATTICE_BGK_THERMAL || solverMethodEnum == MAIA_LATTICE_BGK_INNERENERGY
     || solverMethodEnum == MAIA_LATTICE_BGK_TOTALENERGY || solverMethodEnum == MAIA_LATTICE_BGK_THERMAL_TRANSPORT
     || solverMethodEnum == MAIA_LATTICE_BGK_INNERENERGY_TRANSPORT
     || solverMethodEnum == MAIA_LATTICE_BGK_TOTALENERGY_TRANSPORT) {
    m_cells.setThermal(true);
    m_isThermal = 1;
  } else {
    m_cells.setThermal(false);
  }
  if(solverMethodEnum == MAIA_LATTICE_BGK_TRANSPORT || solverMethodEnum == MAIA_LATTICE_BGK_THERMAL_TRANSPORT
     || solverMethodEnum == MAIA_LATTICE_BGK_INNERENERGY_TRANSPORT
     || solverMethodEnum == MAIA_LATTICE_BGK_TOTALENERGY_TRANSPORT) {
    m_cells.setTransport(true);
    m_isTransport = 1;
  } else {
    m_cells.setTransport(false);
  }
  if(solverMethodEnum == MAIA_LATTICE_BGK_INNERENERGY || solverMethodEnum == MAIA_LATTICE_BGK_INNERENERGY_TRANSPORT) {
    m_innerEnergy = 1;
  } else if(solverMethodEnum == MAIA_LATTICE_BGK_TOTALENERGY
            || solverMethodEnum == MAIA_LATTICE_BGK_TOTALENERGY_TRANSPORT) {
    m_totalEnergy = 1;
  }

  // update number Variables
  m_cells.setNoVariables();

  /*! \page propertyPage1
    \section EELiquid
    <code>MBool LbSolver::m_isEELiquid</code>\n
    default = <code>false</code>\n\n
    Enable the LB-component of the coupled LB-FV Euler-Euler method for bubbly flows.\n
    Keywords: <i>LATTICE_BOLTZMANN, EEMultiphase</i>
  */
  m_isEELiquid = false;
  m_isEELiquid = Context::getSolverProperty<MBool>("EELiquid", m_solverId, AT_, &m_isEELiquid);
  if(m_isEELiquid) {
    MInt l_savePrevVars = 0;
    l_savePrevVars = Context::getBasicProperty<MInt>("alphaConvergenceCheck", AT_, &l_savePrevVars);
    m_cells.setSaveUOtherPhase(true);
    m_cells.setSaveVolumeFraction(true);
    m_cells.setSavePrevVars(l_savePrevVars > 0);
    m_cells.setSaveNuT(true);
    m_cells.setSaveOldNu(true);

    /*! \page propertyPage1
      \section LBRestartWithoutAlpha
      <code>MBool LbSolver::m_EELiquid.restartWithoutAlpha</code>\n
      default = <code>false</code>\n\n
      Restart an E-E simulation from a restart file without alpha.\n
      Keywords: <i>LATTICE_BOLTZMANN, EEMultiphase</i>
    */
    m_EELiquid.restartWithoutAlpha = false;
    m_EELiquid.restartWithoutAlpha =
        Context::getSolverProperty<MBool>("LBRestartWithoutAlpha", m_solverId, AT_, &m_EELiquid.restartWithoutAlpha);
  }

  if(Context::propertyExists("nonNewtonianModel", m_solverId)) {
    m_cells.setSaveOldNu(true);
  }

  {
    MBool updateAfterPropagation = false;
    if(m_isEELiquid) {
      updateAfterPropagation = true;
    }

    /*! \page propertyPage1
      \section updateAfterPropagation
      <code>MBool LbSolver::updateAfterPropagation</code>\n
      default = <code>false</code>\n\n
      Update the macroscopic variables at the end of the timestep.\n
      Keywords: <i>LATTICE_BOLTZMANN</i>
    */
    updateAfterPropagation =
        Context::getSolverProperty<MBool>("updateAfterPropagation", m_solverId, AT_, &updateAfterPropagation);
    if(updateAfterPropagation) {
      if(solverMethodEnum != MAIA_LATTICE_CUMULANT) {
        TERMM(1, "For now only implemented for the cumulant collision step!"); // TODO labels:LB,DOC daniell
      } else {
        m_updateMacroscopicLocation = POSTPROPAGATION;
      }
    }
  }

  m_cells.setNoDistributions(m_noDistributions);

  // Instead of noCells() grid().tree().capacity() was used in m_cells.reset() before;
  // it was changed because there was no difference at that time and to avoid direct access of the gird()
  //  m_cells.reset(grid().noCells());
  m_cells.reset(grid().raw().treeb().capacity());
  m_cells.append(grid().noCells());

  // Update cell collector with info from grid
  updateCellCollectorFromGrid();

  grid().findEqualLevelNeighborsParDiagonal(false);

  // Initialize status flags for adaptation
  m_adaptationSinceLastRestart = false;
  m_adaptationSinceLastSolution = false;

  // Allocate memory for recalcIds if adaptation is used
  /*! \page propertyPage1
    \section adaptation
    <code>MBool LbSolver::adaptation</code>\n
    default = <code>false</code>\n\n
    Switch for using adaptation.\n
    Keywords: <i>LATTICE_BOLTZMANN, ADAPTATION</i>
    */
  m_adaptation = false;
  m_adaptation = Context::getSolverProperty<MBool>("adaptation", m_solverId, AT_, &m_adaptation);

  if(m_adaptation) {
    mAlloc(this->m_recalcIds, maxNoGridCells(), "m_recalcIds", -1, AT_);
    for(MInt i = 0; i < maxNoGridCells(); i++) {
      this->m_recalcIds[i] = i;
    }
  } else {
    this->m_recalcIds = nullptr;
  }

  // Read property cards for lb adaptation
  if(m_adaptation) {
    /*! \page propertyPage1
      \section singleAdaptation
      <code>MBool LbSolver::singleAdaptation</code>\n
      default = <code>true</code>\n\n
      Switch for single adaptation in one adaptation run.\n
      Keywords: <i>LATTICE_BOLTZMANN, ADAPTATION</i>
      */
    this->m_singleAdaptation = true;
    this->m_singleAdaptation =
        Context::getSolverProperty<MBool>("singleAdaptation", m_solverId, AT_, &this->m_singleAdaptation);

    /*! \page propertyPage1
      \section adaptationInitMethod
      <code>MString LbSolver::adaptationInitMethod</code>\n
      default = <code></code>\n\n
      Defines the method for initializing newly created cells after adaptation.\n
      Keywords: <i>LATTICE_BOLTZMANN</i>
      */
    m_adaptationInitMethod = "INIT_FILIPPOVA";
    m_adaptationInitMethod =
        Context::getSolverProperty<MString>("adaptationInitMethod", m_solverId, AT_, &m_adaptationInitMethod);

    if(!domainId()) {
      cerr << "adaptation initialization method is: " << m_adaptationInitMethod << endl;
    }
  }

  // Initialize grid file name and path
  m_reinitFileName = this->grid().gridInputFileName().c_str();
  m_reinitFilePath = outputDir() + m_reinitFileName;

  /*! \page propertyPage1
    \section noSpecies
    <code>MInt LbSolver::noSpecies</code>\n
    default = <code>0</code>\n\n
    Defines the number of species.\n
    Keywords: <i>LATTICE_BOLTZMANN</i>
    */
  MInt noSpecies = 0;
  noSpecies = Context::getSolverProperty<MInt>("noSpecies", m_solverId, AT_, &noSpecies);

  CV = new MConservativeVariables<nDim>(noSpecies);
  PV = new MPrimitiveVariables<nDim>(noSpecies);

  // Set number of variables
  // TODO: Do avoid a linear combination of different implemented methods
  // we should introduce sys equations at some point. In a first step
  // maybe it makes sense to transform m_istTransport into an enum which
  // holds value for 'thermal', 'humidity', or 'whateverInFutureMayAppear'
  // (pre-step for sysEqn like class)
  if(m_isTransport && !m_isThermal) {
    m_noVariables = nDim + 1 + 2 * m_isTransport;
  } else {
    m_noVariables = nDim + 1 + m_isThermal + m_isTransport;
  }

  m_isRefined = (grid().maxUniformRefinementLevel() < maxLevel());

  MString fileName;

  RECORD_TIMER_START(m_t.initSolver);
  m_log << endl;
  m_log << "#########################################################################################################"
           "############"
        << endl;
  m_log << "##                                             Initializing LB solver                                    "
           "          ##"
        << endl;
  m_log << "#########################################################################################################"
           "############"
        << endl;

  /*! \page propertyPage1
    \section nonBlockingComm
    <code>MInt LbSolver::m_nonBlockingComm</code>\n
    default = <code>0</code>\n\n
    This property is a switch for using non-blocking commuication
    <ul>
    <li><code>0</code> (off)</li>
    <li><code>1</code> (on)</li>
    </ul>\n
    Keywords: <i>LATTICE_BOLTZMANN</i>
  */
  m_nonBlockingComm = false;

  m_nonBlockingComm = Context::getSolverProperty<MBool>("nonBlockingComm", m_solverId, AT_, &m_nonBlockingComm);

  // Starting with flow solution
  switch(dist) {
    case 9:
      break;
    case 15:
      break;
    case 19:
      break;
    case 27:
      break;
    default:
      stringstream errorMessage;
      errorMessage << " LbSolver::() Invalid no of distributions : " << dist << " Exiting...";
      TERMM(1, errorMessage.str());
  }

  if(grid().isActive()) {
    // list of active cells
    //  m_activeCellList = new MInt[grid().noCells()];
    mAlloc(m_activeCellList, grid().noCells(), "m_activeCellList", AT_);
  }

  /*! \page propertyPage1
  \section solutionOffset
  <code>MInt solver::m_solutionOffset</code>\n
  default = <code>0</code>\n \n
  which time step to start writing out solution
  Possible values are:
  <ul>
    <li> Int </li>
  </ul>
  Keywords: <i>output</i>
  */
  m_solutionOffset = 0;
  m_solutionOffset = Context::getSolverProperty<MInt>("solutionOffset", m_solverId, AT_, &m_solutionOffset);

  /*! \page propertyPage1
    \section LBinitMethod initMethod
    <code>MString LbSolver::m_initMethod</code>\n
    default = <code>0</code>\n\n
    This property describes the initialization procedure.\n
    Initial variables and distribution functions are set accordingly (in 3d m_initMethodPtr is set)\n
    Additionally the pressure force depends on this setting.\n\n
    Possible values are:\n
    <ul>
    <li><code>"LB_FROM_ZERO_INIT"</code> </li>
    <li><code>"LB_LAMINAR_INIT_PX"</code> </li>
    <li><code>"LB_LAMINAR_INIT_MX"</code> </li>
    <li><code>"LB_LAMINAR_INIT_PY"</code> </li>
    <li><code>"LB_LAMINAR_INIT_MY"</code> </li>
    <li><code>"LB_LAMINAR_INIT_PZ"</code> </li>
    <li><code>"LB_LAMINAR_INIT_MZ"</code> </li>
    <li><code>"LB_LAMINAR_PIPE_INIT"</code> </li>
    <li><code>"LB_VORTEX_INIT"</code> </li>
    <li><code>"LB_LAMINAR_CHANNEL_INIT"</code> </li>
    <li><code>"LB_TURBULENT_CHANNEL_INIT"</code> </li>
    <li><code>"LB_TURBULENT_MIXING_INIT"</code> </li>
    <li><code>"LB_TURBULENT_MIXING_FILTER_INIT"</code> </li>
    <li><code>"LB_TURBULENT_BOUNDARY"</code> </li>
    <li><code>"LB_TURBULENT_PIPE_INIT"</code> </li>
    <li><code>"LB_TURBULENT_DUCT_INIT"</code> </li>
    <li><code>"LB_SOUND_PULSE_INIT"</code> </li>
    <li><code>"LB_SPINNING_VORTICIES_INIT"</code> </li>
    <li><code>"LB_VORTEX_INIT"</code> </li>
    </ul>\n
    Keywords: <i>LATTICE_BOLTZMANN, INITILIZATION</i>
  */
  m_initMethod = Context::getSolverProperty<MString>("initMethod", m_solverId, AT_);
  /*! \page propertyPage1
     \section interpolationType
     <code>MString LbSolver::m_interpolationType</code>\n
     default = <code>"LINEAR_INTERPOLATION"</code>\n\n
     This property sets the interpolation type at grid-refinement interfaces.\n\n

     Possible values are:\n
     <ul>
     <li><code>"LINEAR_INTERPOLATION"</code> </li>
     <li><code>"QUADRATIC_INTERPOLATION"</code> </li>
     <li><code>"CUBIC_INTERPOLATION"</code> </li>
     </ul>\n
     Keywords: <i>LATTICE_BOLTZMANN, INITILIZATION, NUMERICAL_SETUP</i>
   */
  // 1. read as string
  MString interpolationType = "LINEAR_INTERPOLATION";

  interpolationType = Context::getSolverProperty<MString>("interpolationType", m_solverId, AT_, &interpolationType);

  // 2. recast as LbInterpolationType
  m_interpolationType = (LbInterpolationType)string2enum(interpolationType);

  /*! \page propertyPage1
    \section CouettePoiseuilleRatio
    <code>MFloat LbSolver::m_CouettePoiseuilleRatio</code>\n
    default = <code>0</code>\n\n
    Sets the ratio of Couette and Poiseuille flow in a channel.\n
    If m_CouettePoiseuilleRatio = 0 the flow is a Couette flow.\n
    If m_CouettePoiseuilleRatio = inf the flow is a Poiseuille flow.\n
    The equations are implemented in MAIA_LAMINAR_CHANNEL_INIT and bc1002 and bc1022.\n
    Possible values are:\n
    <ul>
    <li><code>any positive float</li>
    </ul>\n
    Keywords: <i>LATTICE_BOLTZMANN, NUMERICAL_SETUP, INITILIZATION</i>
  */
  m_CouettePoiseuilleRatio = F0;
  m_CouettePoiseuilleRatio =
      Context::getSolverProperty<MFloat>("CouettePoiseuilleRatio", m_solverId, AT_, &m_CouettePoiseuilleRatio);

  /*! \page propertyPage1
    \section calcTotalPressureGradient
    <code>MFloat LbSolver::m_calcTotalPressureGradient</code>\n
    default = <code>0</code>\n\n
    Derivative of rho is calculated by default. If calcTotalPressureGradient is on
    the derivative of the total pressure is calculated
    Possible values are:\n
    <ul>
    <li><code>0, 1</li>
    </ul>\n
    Keywords: <i>LATTICE_BOLTZMANN, NUMERICAL_SETUP, INITILIZATION</i>
  */
  m_calcTotalPressureGradient = 0;
  m_calcTotalPressureGradient =
      Context::getSolverProperty<MInt>("calcTotalPressureGradient", m_solverId, AT_, &m_calcTotalPressureGradient);

  /*! \page propertyPage1
    \section densityFluctuations
    <code>MBool LbSolver::m_densityFluctuations</code>\n
    default = <code>0</code>\n\n
    Switch for use of density fluctuations.\n
    If turned on the mean density is zero. This reduces the round-off error.\n
    Some collision steps and boundary conditions are not prepared for this yet.\n\n
    Possible values are:\n
    <ul>
    <li><code>0</code> (off)</li>
    <li><code>1</code> (on)</li>
    </ul>\n
    Keywords: <i>LATTICE_BOLTZMANN, NUMERICAL_SETUP, INITILIZATION</i>
  */
  m_densityFluctuations = false;
  m_densityFluctuations =
      Context::getSolverProperty<MBool>("densityFluctuations", m_solverId, AT_, &m_densityFluctuations);

  /*! \page propertyPage1
    \section calculateDissipation
    <code>MBool LbSolver::m_calculateDissipation</code>\n
    default = <code>0</code>\n\n
    Switch for calculation of dissipation.\n
    Values of dissipation and energy are written to dissipation.dat in each timestep.\n
    This may slow down the solver severely!\n\n
    Possible values are:\n
    <ul>
    <li><code>0</code> (off)</li>
    <li><code>1</code> (on)</li>
    </ul>\n
    Keywords: <i>LATTICE_BOLTZMANN, INPUT_OUTPUT</i>
  */
  m_calculateDissipation = false;
  m_calculateDissipation =
      Context::getSolverProperty<MBool>("calculateDissipation", m_solverId, AT_, &m_calculateDissipation);

  /*! \page propertyPage1
    \section FFTInit
    <code>MInt LbSolver::m_FftInit</code>\n
    default = <code>false</code>\n\n
    This property is a switch for Fourier-transform initialization.\n\n
    Possible values are:
    <ul>
    <li><code>0</code> (off)</li>
    <li><code>1</code> (on)</li>
    </ul>\n
    Keywords: <i>LATTICE_BOLTZMANN, INITILIZATION</i>
  */
  m_FftInit = false;
  m_FftInit = Context::getSolverProperty<MBool>("FFTInit", m_solverId, AT_, &m_FftInit);

  m_fftInterval =
      Context::propertyExists("fftInterval", 0) ? Context::getSolverProperty<MInt>("fftInterval", m_solverId, AT_) : 0;

  /*! \page propertyPage1
    \section domainSize
    <code>MInt LbSolver::m_arraySize</code>\n
    default = <code>-</code>\n\n
    Array size for FFTInit given on the highest level of refinement.\n
    Values must not be odd numbers!\n\n
    Possible values in 3d are:
    <ul>
    <li><code>24,64,30</code></li>
    </ul>\n
    Keywords: <i>LATTICE_BOLTZMANN, NUMERICAL_SETUP, INITILIZATION</i>
  */
  for(MInt dir = 0; dir < nDim; dir++) {
    m_arraySize[dir] = 2;
    m_arraySize[dir] = Context::getSolverProperty<MInt>("arraySize", m_solverId, AT_, &m_arraySize[dir], dir);
  }

  /*! \page propertyPage1
    \section noPeakModes
    <code>MInt LbSolver::m_noPeakModes</code>\n
    default = <code>1</code>\n\n

    Defines how often the FFT-mode with the highest energy fits into the x-length of the domain.\n\n
    Keywords: <i>LATTICE_BOLTZMANN, INITILIZATION</i>
  */
  m_noPeakModes = 1;
  m_noPeakModes = Context::getSolverProperty<MInt>("m_noPeakModes", m_solverId, AT_, &m_noPeakModes);

  if(m_FftInit && m_noPeakModes == 0) {
    stringstream errorMessage;
    errorMessage << " noPeakModes was not defined for FFTInit, exiting";
    TERMM(1, errorMessage.str());
  }

  // ------------------------------------
  // subgrid properties
  // ------------------------------------

  /*! \page propertyPage1
    \section smagorinskyConstant
    <code>MFloat LbSolver::m_Cs</code>\n
    default = <code>0.1</code>\n\n
    This property defines the Smagorinsky constant for LES computations.\n
    Keywords: <i>LATTICE_BOLTZMANN</i>
  */
  m_Cs = 0.1;
  m_Cs = Context::getSolverProperty<MFloat>("smagorinskyConstant", m_solverId, AT_, &m_Cs);

  /*! \page propertyPage1
    \section filterWidth
    <code>MFloat LbSolver::m_deltaX</code>\n
    default = <code>1.0</code>\n\n
    This property defines the Smagorinsky filter width for LES computations.\n
    Keywords: <i>LATTICE_BOLTZMANN</i>
  */
  m_deltaX = 1.0;
  m_deltaX = Context::getSolverProperty<MFloat>("filterWidth", m_solverId, AT_, &m_deltaX);

  // LES related fields
  if(solverMethodEnum == MAIA_LATTICE_BGKI_SMAGORINSKY || solverMethodEnum == MAIA_LATTICE_BGKI_SMAGORINSKY2
     || solverMethodEnum == MAIA_LATTICE_BGKI_SMAGO_WALL || solverMethodEnum == MAIA_LATTICE_BGKI_DYNAMIC_SMAGO
     || solverMethodEnum == MAIA_LATTICE_BGK_INIT || solverMethodEnum == MAIA_LATTICE_RBGK_SMAGORINSKY
     || solverMethodEnum == MAIA_LATTICE_RBGK_DYNAMIC_SMAGO || solverMethodEnum == MAIA_LATTICE_BGK) {
    // Allocate space for momentum flux tensor
    mAlloc(m_momentumFlux, a_noCells(), nDim * nDim, "m_momentumFlux", F0, AT_);
    // Allocate space for SGS tensors
    mAlloc(m_MijMij, a_noCells(), "m_MijMij", F0, AT_);
    mAlloc(m_MijLij, a_noCells(), "m_MijLij", F0, AT_);
  }

  //------------------------------------
  // Thermal properties
  //------------------------------------

  /*! \page propertyPage1
   \section Pr
   <code>MFloat LbSolver::m_Pr</code>\n
   default = <code>0.72</code>\n\n
   This property defines the Prandtl number, which is required to calculate the heat coefiicient k.
   Example numbers are:
   <ul>
   <li>Air: <code>0.7-0.8</code></li>
   <li>Water: <code>7</code></li>
   </ul>
   Keywords: <i>THERMAL_LATTICE_BOLTZMANN</i>
 */
  m_Pr = 0.72;
  m_Pr = Context::getSolverProperty<MFloat>("Pr", m_solverId, AT_, &m_Pr);

  /*! \page propertyPage1
    \section initTemperatureKelvin
    <code>MFloat LbSolver::m_initTemperatureKelvin</code>\n
    default = <code>1.0546</code>\n\n
    This property defines the dimensionless temperature in Kelvin.
    \n\n
    Keywords: <i>THERMAL_LATTICE_BOLTZMANN</i>
  */
  m_initTemperatureKelvin = 1.0546;
  m_initTemperatureKelvin =
      Context::getSolverProperty<MFloat>("initTemperatureKelvin", m_solverId, AT_, &m_initTemperatureKelvin);

  /*! \page propertyPage1
    \section blasiusPos
    <code>MFloat LbSolver::m_blasiusPos</code>\n
    default = <code>0.5</code>\n\n
    This property defines the position on a flat plate where to evaluate a Blasius solution
    for the inflow condition.
    \n\n
    Keywords: <i>THERMAL_LATTICE_BOLTZMANN</i>
  */
  m_blasiusPos = 0.1;
  m_blasiusPos = Context::getSolverProperty<MFloat>("blasiusPos", m_solverId, AT_, &m_blasiusPos);

  //------------------------------------
  // Transport properties
  //------------------------------------

  /*! \page propertyPage1
   \section Pe
   <code>MFloat LbSolver::m_Pe</code>\n
   default = <code>100</code>\n\n
   This property defines the Peclet number, which is used to compute the diffusivity.
   </ul>
   Keywords: <i>TRANSPORT_LATTICE_BOLTZMANN</i>
 */
  m_Pe = 100;
  m_Pe = Context::getSolverProperty<MFloat>("Pe", m_solverId, AT_, &m_Pe);

  /*! \page propertyPage1
    \section initCon
    <code>MFloat LbSolver::m_initCon</code>\n
    default = <code>1.0</code>\n\n
    This property defines the dimensionless concentration of the passively transported media.
    \n\n
    Keywords: <i>TRANSPORT_LATTICE_BOLTZMANN</i>
  */
  m_initCon = 1.0;
  m_initCon = Context::getSolverProperty<MFloat>("initCon", m_solverId, AT_, &m_initCon);

  //------------------------------------

  /*! \page propertyPage1
    \section alpha
    <code>MFloat LbSolver::m_alpha</code>\n
    default = <code>0.0</code>\n\n
    This property defines the Womerleynumber for unsteady flow.
    \n\n
    Keywords: <i>LATTICE_BOLTZMANN</i>
  */
  m_alpha = 0.0;
  m_alpha = Context::getSolverProperty<MFloat>("alpha", m_solverId, AT_, &m_alpha);

  /*! \page propertyPage1
    \section saveDerivatives
    <code>MBool LbSolver::m_saveDerivatives</code>\n
    default = <code>0</code>\n\n
    This property defines if derivatives should be additionally stored.
    \n\n
    Keywords: <i>LATTICE_BOLTZMANN</i>
  */
  m_saveDerivatives = false;
  if(Context::propertyExists("saveDerivatives", m_solverId))
    m_saveDerivatives = Context::getSolverProperty<MBool>("saveDerivatives", m_solverId, AT_);

  // ------------------------------------

  /*! \page propertyPage1
    \section tanhInit
    <code>MBool LbSolver::m_tanhInit</code>\n
    default = <code>0</code>\n\n
    This property defines if the Reynolds number should be increased according to a tanh-function.\n
    Keywords: <i>LATTICE_BOLTZMANN</i>
  */
  m_tanhInit = false;
  m_tanhInit = Context::getSolverProperty<MBool>("tanhInit", m_solverId, AT_, &m_tanhInit);

  if(m_tanhInit) {
    // the follwoing must be defined, otherwise exit with error

    /*! \page propertyPage1
      \section initRe
      <code>MFloat LbSolver::m_initRe</code>\n
      This property defines if the Reynolds number to start with for tanh-function increase
      Keywords: <i>LATTICE_BOLTZMANN</i>
    */
    m_initRe = Context::getSolverProperty<MFloat>("initRe", m_solverId, AT_);

    /*! \page propertyPage1
      \section initTime
      <code>MInt LbSolver::m_initTime</code>\n
      This property defines the number of LB iterations for Reynolds number tanh-function increase.\n
      Keywords: <i>LATTICE_BOLTZMANN</i>
    */
    m_initTime = Context::getSolverProperty<MInt>("initTime", m_solverId, AT_);

    /*! \page propertyPage1
      \section initStartTime
      <code>MInt LbSolver::m_initStartTime</code>\n
      This property defines the LB iteration number to start a Reynolds number tanh-function increase.\n
      Keywords: <i>LATTICE_BOLTZMANN</i>
    */
    m_initStartTime = Context::getSolverProperty<MInt>("initStartTime", m_solverId, AT_);
  }

  // ------------------------------------

  /*! \page propertyPage1
    \section Ma
    <code>MFloat LbSolver::m_Ma</code>\n
    default = <code>0.1</code>\n\n
    This property defines the Mach number of the LB simulation.
    Note that due to the LB's quadi-incompressibility, only Mach numbers of up to Ma=0.3 should be used.
    Keywords: <i>LATTICE_BOLTZMANN</i>
  */
  m_Ma = 0.1;
  m_Ma = Context::getSolverProperty<MFloat>("Ma", m_solverId, AT_, &m_Ma);

  /*! \page propertyPage1
    \section Re
    <code>MFloat LbSolver::m_Re</code>\n
    default = <code>100.0</code>\n\n
    This property defines the Reynolds number of the LB simulation.
    The Reynolds number, together with the velocity that is calculated from the Mach number, and the
    reference length are used to calculate the viscosity of the fluid in the LB.
    Keywords: <i>LATTICE_BOLTZMANN</i>
  */
  m_Re = 100.0;
  m_Re = Context::getSolverProperty<MFloat>("Re", m_solverId, AT_, &m_Re);

  if(Context::propertyExists("ReTau", m_solverId)) {
    m_ReTau = Context::getSolverProperty<MFloat>("ReTau", m_solverId, AT_);
  } else {
    m_Re = Context::getSolverProperty<MFloat>("Re", m_solverId, AT_);
  }

  /*! \page propertyPage1
    \section rho1
    <code>MFloat LbSolver::m_rho1</code>\n
    default = <code>1.0</code>\n\n
    This property defines the reference density rho1.\n
    Keywords: <i>LATTICE_BOLTZMANN</i>
  */
  m_rho1 = 1.0;
  m_rho1 = Context::getSolverProperty<MFloat>("rho1", m_solverId, AT_, &m_rho1);

  /*! \page propertyPage1
    \section rho2
    <code>MFloat LbSolver::m_rho2</code>\n
    default = <code>1.0</code>\n\n
    This property defines the reference density rho2.\n
    Keywords: <i>LATTICE_BOLTZMANN</i>
  */
  m_rho2 = 1.0;
  m_rho2 = Context::getSolverProperty<MFloat>("rho2", m_solverId, AT_, &m_rho2);

  // ------------------------------------

  // Calculate cell lengths
  if(grid().isActive()) {
    mAlloc(m_cellLength, maxLevel() + 1, "m_cellLength", -F1, AT_);
    for(MInt level = 0; level < maxLevel() + 1; level++) {
      m_cellLength[level] = c_cellLengthAtLevel(level);
    }
  }

  {
    // Sanity check for current implementation
    /*if(isActive() && (grid().maxRefinementLevel() != maxLevel())) {
      TERMM(1,
            "WARNING: Local max level and possible "
            "maxRefinementLevel from the property file are not consistent.");
    }*/
    // reference lengths
    MInt chk = 0;
    stringstream ss1;
    if(Context::propertyExists("referenceLength", m_solverId)) {
      /*! \page propertyPage1
        \section referenceLength
        <code>MFloat LbSolver::m_referenceLength </code>\n
        This property sets the reference length in stl units, whereas internally
        the reference length is converted to LB units.
        Keywords: <i>LATTICE_BOLTZMANN</i>
      */
      const MFloat referenceLengthSTL =
          Context::getSolverProperty<MFloat>("referenceLength", m_solverId, AT_, &referenceLengthSTL);
      const MFloat lengthMaxRefinementLvl = c_cellLengthAtLevel(grid().maxRefinementLevel());
      m_referenceLengthSTL = referenceLengthSTL;
      m_referenceLength = referenceLengthSTL / lengthMaxRefinementLvl;
      ss1 << "referenceLength ";
      chk++;
    }
    if(Context::propertyExists("referenceLengthLB", m_solverId)) {
      /*! \page propertyPage1
        \section referenceLengthLB
        <code>MFloat LbSolver::m_referenceLength</code>\n
        default = <code>is calculated by calculateReferenceLength based on the reference length</code>\n\n
        This property overrides the reference length in cell units.
        Keywords: <i>LATTICE_BOLTZMANN</i>
      */
      m_referenceLength = Context::getSolverProperty<MFloat>("referenceLengthLB", m_solverId, AT_);
      ss1 << "referenceLengthLB ";
      chk++;
    }
    if(Context::propertyExists("referenceLengthSegId", m_solverId)) {
      /*! \page propertyPage1
        \section referenceLengthSegId
        <code>MFloat LbSolver::m_referenceLengthSegId</code>\n
        This property is mandatory if property referenceLengthLB is not given. Defines the segment id to be used for the
        calculation of the characteristic length in cell units.
        Keywords: <i>LATTICE_BOLTZMANN</i>
      */
      m_referenceLengthSegId = Context::getSolverProperty<MInt>("referenceLengthSegId", m_solverId, AT_);
      m_referenceLength = calculateReferenceLength(m_referenceLengthSegId);
      ss1 << "referenceLengthSegId ";
      chk++;
    }
    if(chk == 0) {
      stringstream ss2;
      ss2 << "One of the following properties must be given: ";
      ss2 << "referenceLength referenceLengthLB referenceLengthSegId";
      TERMM(1, ss2.str());
    } else if(chk != 1) {
      stringstream ss2;
      ss2 << "Redundant conflicting input properties: ";
      ss2 << ss1.str();
      TERMM(1, ss2.str());
    }

    if(grid().isActive()) {
      chk = 0;
      m_domainLength = FPOW2(maxLevel()) / reductionFactor();
      if(Context::propertyExists("domainLength", m_solverId)) {
        /*! \page propertyPage1
        \section domainLength
        <code>MFloat LbSolver::m_domainLength</code>\n\n
        This property overrides the maximum length in stl units.
        Internally the domain length is converted to LB units.
        Keywords: <i>LATTICE_BOLTZMANN</i>
        */
        const MFloat domainLengthSTL =
            Context::getSolverProperty<MFloat>("domainLength", m_solverId, AT_, &domainLengthSTL);
        const MFloat lengthMaxRefinementLvl = c_cellLengthAtLevel(grid().maxRefinementLevel());
        m_domainLength = domainLengthSTL / lengthMaxRefinementLvl;
        chk++;
      }
      if(Context::propertyExists("domainLengthLB", m_solverId)) {
        /*! \page propertyPage1
        \section domainLengthLB
        <code>MFloat LbSolver::m_domainLength</code>\n
        default = <code>2^maxLevel() / reductionFactor()</code>\n\n
        This property overrides the maximum length in cell units.
        Keywords: <i>LATTICE_BOLTZMANN</i>
        */
        m_domainLength = Context::getSolverProperty<MFloat>("domainLengthLB", m_solverId, AT_, &m_domainLength);
        chk++;
      }
      if(chk > 1) {
        TERMM(1, "Redundant conflicting input properties: domainLength, domainLengthLB");
      }
    }
  }

  // if external forcing is turned off, the forcing terms are set to zero
  mAlloc(m_Fext, m_noDistributions, "m_Fext", F0, AT_);

  mAlloc(m_EELiquid.Fg, m_noDistributions, "m_EELiquid.Fg", F0, AT_);

  for(MInt mi = 0; mi < m_noDistributions; mi++) {
    m_EELiquid.Fg[mi] = F0;
  }

  if(m_isEELiquid) {
    if(!Context::propertyExists("EELiquidGravity", m_solverId)
       || !Context::propertyExists("EELiquidGravityAccel", m_solverId))
      TERMM(1, "Missing property EELiquidGravity or EELiquidGravityAccel!");

    /*! \page propertyPage1
      \section initialAlpha
      <code>MFloat LbSolver::m_EELiquid.initialAlpha</code>\n
      default = <code>0.0</code>\n\n
      Initial value of alpha in the domain.\n
      Keywords: <i>LATTICE_BOLTZMANN, EEMultiphase</i>
    */
    m_EELiquid.initialAlpha = 0.0;
    m_EELiquid.initialAlpha =
        Context::getSolverProperty<MFloat>("initialAlpha", m_solverId, AT_, &m_EELiquid.initialAlpha);

    /*! \page propertyPage1
      \section alphaInf
      <code>MFloat LbSolver::m_EELiquid.alphaInf</code>\n
      default = <code>m_EELiquid.initialAlpha</code>\n\n
      Infinity value of alpha.\n
      Keywords: <i>LATTICE_BOLTZMANN, EEMultiphase</i>
    */
    m_EELiquid.alphaInf = m_EELiquid.initialAlpha;
    m_EELiquid.alphaInf = Context::getSolverProperty<MFloat>("alphaInf", m_solverId, AT_, &m_EELiquid.alphaInf);

    if(domainId() == 0) cerr << "LB Solver EELiquid!" << endl;
  }

  /*! \page propertyPage1
    \section EELiquidGravity
    <code>MBool LbSolver::m_EELiquid.gravity</code>\n
    default = <code>false</code>\n\n
    Enable the influence of buoyancy on the liquid phase.\n
    Keywords: <i>LATTICE_BOLTZMANN, EEMultiphase</i>
  */
  m_EELiquid.gravity = false;
  m_EELiquid.gravity = Context::getSolverProperty<MBool>("EELiquidGravity", m_solverId, AT_, &m_EELiquid.gravity);

  /*! \page propertyPage1
    \section initDensityGradient
    <code>MBool LbSolver::m_initDensityGradient</code>\n
    default = <code>false</code>\n\n
    Initialize the density with a gradient according to m_volumeAccel.\n
    Keywords: <i>LATTICE_BOLTZMANN</i>
  */
  m_initDensityGradient = false;
  m_initDensityGradient =
      Context::getSolverProperty<MBool>("initDensityGradient", m_solverId, AT_, &m_initDensityGradient);

  if(m_EELiquid.gravity) {
    if(!Context::propertyExists("EELiquidGravityAccel", m_solverId)) {
      TERMM(1, "Missing property EELiquidGravityAccel!");
    }
    for(MInt i = 0; i < nDim; i++) {
      /*! \page propertyPage1
        \section EELiquidGravityAccel
        <code>MFloat LbSolver::m_EELiquid.gravityAccelM[nDim]</code>\n
        Gravity acceleration used for buoyancy effects.\n
        Gravity has to be specified in non-dimensional form by multiplying by (\delta t ^ 2) / (\delta x)
                                                                            = (\delta x) / (3 * a_inf^2) (see wiki)
        Keywords: <i>LATTICE_BOLTZMANN, EEMultiphase</i>
      */
      m_EELiquid.gravityAccelM[i] = 0.0;
      m_EELiquid.gravityAccelM[i] = Context::getSolverProperty<MFloat>("EELiquidGravityAccel", m_solverId, AT_, i);
    }
  }

  m_externalForcing = false;
  if(Context::propertyExists("volumeAcceleration", m_solverId) || Context::propertyExists("Ga", m_solverId)) {
    m_externalForcing = true;

    if(Context::propertyExists("volumeAcceleration", m_solverId)) {
      for(MInt i = 0; i < nDim; i++) {
        /*! \page propertyPage1
          \section volumeAcceleration
          <code>MFloat LbSolver::m_volumeAccel</code>\n
          default = <code>[0.0, 0.0, 0.0]</code>\n\n
          This property defines the amount of acceleration applied in each Cartesian direction
          The acceleration has to be specified in LB non-dimensional form
          by multiplying by (\delta t ^ 2) / (\delta x) = (\delta x) / (3 * a_inf^2) (see wiki)\n
          Keywords: <i>LATTICE_BOLTZMANN</i>
        */
        m_volumeAccel[i] += Context::getSolverProperty<MFloat>("volumeAcceleration", m_solverId, AT_, i);
      }
    }
    /*! \page propertyPage1
      \section Ga
      <code>MFloat LbSolver::m_Ga</code>\n
      default = <code>0.0</code>\n\n
      Galileo number\n
      Gravity is applied in negative y-direction.\n
      Keywords: <i>LATTICE_BOLTZMANN</i>
    */
    if(Context::propertyExists("Ga", m_solverId)) {
      m_Ga = Context::getSolverProperty<MFloat>("Ga", m_solverId, AT_, &m_Ga);
      const MFloat nu = m_Ma / F1BCS / m_Re * m_referenceLength;
      const MFloat gravity = POW2(m_Ga) * POW2(nu) / POW3(m_referenceLength);
      m_volumeAccel[1] -= gravity;
    }
  }

  /*! \page propertyPage1
    \section externalForcing
    <code>MBool LbBndSolver::m_externalForcing</code>\n
    This property defines if external forcing should be activated.\n
    <ul>
    <li><code>0</code> (off)</li>
    <li><code>1</code> (on)</li>
    </ul>\n
    Keywords: <i>LATTICE_BOLTZMANN</i>
  */
  m_externalForcing = Context::getSolverProperty<MBool>("externalForcing", m_solverId, AT_, &m_externalForcing);

  m_particleMomentumCoupling = false;
  m_particleMomentumCoupling =
      Context::getSolverProperty<MBool>("particleMomentumCoupling", m_solverId, AT_, &m_particleMomentumCoupling);

  m_saveExternalForces = false;
  m_saveExternalForces =
      Context::getSolverProperty<MBool>("saveExternalForces", m_solverId, AT_, &m_saveExternalForces);

  /*m_gravity = POW2(m_Ma * LBCS) / POW2(m_Fr) / m_referenceLength;
  m_Ga = Context::getSolverProperty<MFloat>("Ga", m_solverId, AT_, &m_Ga);*/

  /*! \page propertyPage1
    \section velocityControl
    <code>MInt LbSolver::m_velocityControl.dir</code>\n
    default = <code>-1</code>\n\n
    This property defines if velocityControl should be activated and what the main direction of flow is.
    The velocity control algorithm averages the velocity in the corresponding direction over the whole domain.
    The volumeAcceleration is subsequently controlled using a PID-controller.
    The target for the averaged velocity is u_infinity.
    <ul>
    <li><code>-1</code> (off)</li>
    <li><code>0</code> (on), (positive) x-direction</li>
    <li><code>1</code> (on), (positive) y-direction</li>
    <li><code>2</code> (on), (positive) z-direction</li>
    </ul>\n
    Keywords: <i>LATTICE_BOLTZMANN</i>
  */
  m_velocityControl.dir = -1;
  m_velocityControl.dir = Context::getSolverProperty<MInt>("velocityControl", m_solverId, AT_, &m_velocityControl.dir);
  if(m_velocityControl.dir < -1 || m_velocityControl.dir > 2) {
    TERMM(1, "Invalid Cartesian direction for velocityControl!");
  }
  if(m_velocityControl.dir >= 0) {
    m_controlVelocity = true;
  }

  /*! \page propertyPage1
      \section velocityControlRestart
      <code>MBool LbSolver::m_velocityControl.restart</code>\n

      Set this property to true if you restart a run with velocityControl to read values like previousError
      and integratedError fromm the restart file. Else these values are set to zero.
  */
  m_velocityControl.restart = false;
  m_velocityControl.restart =
      Context::getSolverProperty<MBool>("velocityControlRestart", m_solverId, AT_, &m_velocityControl.restart);

  /*! \page propertyPage1
    \section velocityControlInterval
    <code>MInt LbSolver::m_velocityControlInterval</code>\n
    default = <code>100</code>\n\n
    This property defines the interval in which the velocity is averaged
    and the volumeforce is adjusted for velocityControl\n
    Keywords: <i>LATTICE_BOLTZMANN</i>
  */
  m_velocityControl.interval = 100;
  m_velocityControl.interval =
      Context::getSolverProperty<MInt>("velocityControlInterval", m_solverId, AT_, &m_velocityControl.interval);
  if(m_velocityControl.interval < 1) {
    TERMM(1, "Invalid velocityControlInterval!");
  }

  /*! \page propertyPage1
    \section velocityControlKT
    <code>MFloat LbSolver::m_velocityControlKT</code>\n
    default = <code>1.0</code>\n\n
    velocityControl works as a PID-controller following the formula:
    controlSignal = KT * (err + 1/KI * integral(err)dt + KD * (d err / dt))\n
    m_volumeAccel = m_volumeAccelBase - controlSinal * m_volumeAccelBase (in the dimension of m_velocityControl)
    Keywords: <i>LATTICE_BOLTZMANN</i>
  */
  m_velocityControl.KT = 1.0;
  m_velocityControl.KT =
      Context::getSolverProperty<MFloat>("velocityControlKT", m_solverId, AT_, &m_velocityControl.KT);

  /*! \page propertyPage1
    \section velocityControlKI
    <code>MFloat LbSolver::m_velocityControlKI</code>\n
    default = <code>10000.0</code>\n\n
    velocityControl works as a PID-controller following the formula:
    controlSignal = KT * (err + 1/KI * integral(err)dt + KD * (d err / dt))\n
    m_volumeAccel = m_volumeAccelBase - controlSinal * m_volumeAccelBase (in the dimension of m_velocityControl)
    Keywords: <i>LATTICE_BOLTZMANN</i>
  */
  m_velocityControl.KI = 10000.0;
  m_velocityControl.KI =
      Context::getSolverProperty<MFloat>("velocityControlKI", m_solverId, AT_, &m_velocityControl.KI);

  /*! \page propertyPage1
    \section velocityControlKD
    <code>MFloat LbSolver::m_velocityControlKD</code>\n
    default = <code>10.0</code>\n\n
    velocityControl works as a PID-controller following the formula:
    controlSignal = KT * (err + 1/KI * integral(err)dt + KD * (d err / dt))\n
    m_volumeAccel = m_volumeAccelBase - controlSinal * m_volumeAccelBase (in the dimension of m_velocityControl)
    Keywords: <i>LATTICE_BOLTZMANN</i>
  */
  m_velocityControl.KD = 10.0;
  m_velocityControl.KD =
      Context::getSolverProperty<MFloat>("velocityControlKD", m_solverId, AT_, &m_velocityControl.KD);

  if(m_controlVelocity) {
    for(MInt i = 0; i < nDim; i++) {
      m_volumeAccelBase[i] = m_volumeAccel[i];
    }
  }

  /*! \page propertyPage1
    \section solidLayer
    <code>MBool LbSolver::m_solidLayerExtension</code>\n
    default = <code>false</code>\n\n
    enables special treatment if solid cells outside of the domain are not deleted during grid generation
    Keywords: <i>LATTICE_BOLTZMANN</i>
  */
  m_solidLayerExtension = false;
  m_solidLayerExtension = Context::getSolverProperty<MBool>("solidLayer", m_solverId, AT_, &m_solidLayerExtension);

  /*! \page propertyPage1
    \section writeLsData
    <code>MBool LbSolver::m_writeLsData</code>\n
    default = <code>false</code>\n\n
    If enabled, the Level-Set data stored in the LB solver is also written out, this includes
    the Level-Set, the Body Id and the isActive state
    Keywords: <i>LATTICE_BOLTZMANN</i>
  */
  m_writeLsData = false;
  m_writeLsData = Context::getSolverProperty<MBool>("writeLsData", m_solverId, AT_, &m_writeLsData);

  /*! \page propertyPage1
    \section useOnlyCollectedLS
    <code>MBool LbSolver::m_useOnlyCollectedLS</code>\n
    default = <code>false</code>\n\n
    If enabled, only the 0th level-set is transferred to the LB solver.
    Keywords: <i>LATTICE_BOLTZMANN</i>
  */
  m_useOnlyCollectedLS = false;
  m_useOnlyCollectedLS =
      Context::getSolverProperty<MBool>("useOnlyCollectedLs", m_solverId, AT_, &m_useOnlyCollectedLS);

  /*! \page propertyPage1
    \section allowBndryAsG0
    <code>MBool LbSolver::m_allowBndryAsG0</code>\n
    default = <code>false</code>\n\n
    If enabled, bndry cells can be G0 candidates / G0 boundary cells.
    Keywords: <i>LATTICE_BOLTZMANN</i>
  */
  m_allowBndryAsG0 = false;
  m_allowBndryAsG0 = Context::getSolverProperty<MBool>("allowBndryAsG0", m_solverId, AT_, &m_allowBndryAsG0);

  // allocate and initialize temporary residual and variable arrays :
  mAlloc(m_rescoordinates, nDim + 1, nDim, "m_rescoordinates", F0, AT_);
  mAlloc(m_residual, nDim + 1, "m_residual", F0, AT_);
  mAlloc(m_tmpResidual, nDim + 1, "m_tmpResidual", F0, AT_);
  mAlloc(m_tmpResidualLvl, nDim + 1, "m_tmpResidualLvl", 0, AT_);
  mAlloc(m_maxResId, nDim + 1, "m_ResId", 0, AT_);

  // Calculate simulation setup variables
  m_finalRe = m_Re;
  if(m_tanhInit) {
    m_Re = m_initRe;
    m_tanhScaleFactor = 1.0 / (1.0 - (tanh(-2.5) + 1.0));
  }

  m_nu = m_Ma * LBCS / m_Re * m_referenceLength;
  m_omega = 2.0 / (1.0 + 6.0 * m_nu);


  if(Context::propertyExists("nonNewtonianModel", m_solverId)) {
    m_nonNewtonian = true;
    const MString modelString = Context::getSolverProperty<MString>("nonNewtonianModel", m_solverId, AT_);
    const auto model = string2enum(modelString);
    if(model == CARREAU) {
      m_n = Context::getSolverProperty<MFloat>("nonNewtonian_n_bndry", m_solverId, AT_);
    } else if(model == POWERLAW) {
      m_n = Context::getSolverProperty<MFloat>("nonNewtonian_n", m_solverId, AT_);
    } else {
      std::stringstream ss;
      ss << "nonNewtonianModel : " << modelString << " is not implemented!" << std::endl;
      TERMM(1, ss.str());
    }
  }

  if(m_isThermal != 0) {
    m_kappa = m_nu / m_Pr;
    m_omegaT = 2.0 / (1.0 + 6.0 * m_kappa);
  }

  if(m_isTransport != 0) {
    m_diffusivity = m_nu * (m_Re / m_Pe);
    m_omegaD = 2.0 / (1.0 + 6.0 * m_diffusivity);
  }

  // pulsatile frequency
  m_pulsatileFrequency = m_nu * ((m_alpha * m_alpha) / (m_referenceLength * m_referenceLength));

  RECORD_TIMER_STOP(m_t.initSolver);

  // interface cell treatment
  if(m_isRefined || m_adaptation) {
    /*! \page propertyPage1
    \section correctInterfaceBcCells
    <code>MBool LbSolver::m_correctInterfaceBcCells</code>\n
    default = <code>false</code>\n\n
    Activating this properties toggle interface cells on boundaries including
    shifting and correction of the interpolation stencil and its coefficients
    <ul>
    <li><code>false</code>li>
    <li><code>true</code></li>
    </ul>\n
    Keywords: <i>LATTICE BOLTZMANN</i>
    */
    m_correctInterfaceBcCells =
        Context::getSolverProperty<MBool>("correctInterfaceBcCells", m_solverId, AT_, &m_correctInterfaceBcCells);
    if(grid().isActive()) {
      treatInterfaceCells();
    }
  }

  // set active cells
  if(grid().isActive()) {
    setActiveCellList();
  }

  // load restart file if necessary
  m_initRestart = false;
  m_initFromCoarse = false;
  if(m_restartFile) {
    /*! \page propertyPage1
    \section initRestart
    <code>MBool LbSolver::m_initRestart</code>\n
    default = <code>0</code>\n\n
    This property defines if the PPDFs should be initialized
    with the eq. PPDFs of the given PVs loaded from restart file
    <ul>
    <li><code>0</code> (off)</li>
    <li><code>1</code> (on)</li>
    </ul>\n
    Keywords: <i>LATTICE BOLTZMANN</i>
    */
    m_initRestart = Context::getSolverProperty<MBool>("initRestart", m_solverId, AT_);

    /*! \page propertyPage1
    \section initFromCoarse
    <code>MBool LbSolver::m_initFromCoarse</code>\n
    default = <code>0</code>\n\n
    This property defines if a restart should be initialized
    from a solution on a coarser mesh.
    <ul>
    <li><code>0</code> (off)</li>
    <li><code>1</code> (on)</li>
    </ul>\n
    Keywords: <i>LATTICE BOLTZMANN</i>
    */
    m_initFromCoarse = Context::getSolverProperty<MBool>("initFromCoarse", m_solverId, AT_, &m_initFromCoarse);
  }

  /*! \page propertyPage1
  \section isInitRun
  <code>MBool LbSolver::m_isInitRun</code>\n
  default = <code>false</code>\n\n
  This property defines if an initRun is performed. That is, the initialized
  velocity field is kept constant, while the density is iterated (Analogue to
  solving Poisson equation). Hereby, a correct density and non-equilibrium field
  is obtained.
  Ref.: Mei et al. 2006: https://doi.org/10.1016/j.compfluid.2005.08.008
  <ul>
  <li><code>true</code> (off)</li>
  <li><code>false</code> (on)</li>
  </ul>\n
  Keywords: <i>LATTICE BOLTZMANN</i>
  */
  m_isInitRun = false;
  m_isInitRun = Context::getSolverProperty<MBool>("isInitRun", m_solverId, AT_, &m_isInitRun);

  MString lbInterfaceMethod = "FILIPPOVA";
  lbInterfaceMethod = Context::getSolverProperty<MString>("interfaceMethod", m_solverId, AT_, &lbInterfaceMethod);

  m_log << endl;
  m_log << "#########################################################################################################"
           "#############"
        << endl;
  m_log << "##                                             Methods and flow field init                               "
           "           ##"
        << endl;
  m_log << "#########################################################################################################"
           "#############"
        << endl
        << endl;
  m_log << "  + Initializing methods..." << endl;
  m_log << "    - solver method:     " << solverMethod() << endl;
  m_log << "    - interface method: " << lbInterfaceMethod << endl << endl;

  switch(string2enum(lbInterfaceMethod)) {
    case FILIPPOVA:
      if(m_isThermal && !m_isTransport) {
        m_propagationStepMethod = &LbSolver::propagation_step_thermal;
      } else if(!m_isThermal && m_isTransport) {
        m_propagationStepMethod = &LbSolver::propagation_step_transport;
      } else if(m_isThermal && m_isTransport) {
        m_propagationStepMethod = &LbSolver::propagation_step_thermaltransport;
      } else {
        m_propagationStepMethod = &LbSolver::propagation_step;
      }
      break;
    case ROHDE:
      if(m_isThermal && !m_isTransport) {
        m_propagationStepMethod = &LbSolver::propagation_step_thermal_vol;
      } else if(!m_isThermal && m_isTransport) {
        m_propagationStepMethod = &LbSolver::propagation_step_transport_vol;
      } else if(m_isThermal && m_isTransport) {
        m_propagationStepMethod = &LbSolver::propagation_step_thermaltransport_vol;
      } else {
        m_propagationStepMethod = &LbSolver::propagation_step_vol;
      }
      break;
    default: {
      TERMM(1, "Unknown lb interface method: Exiting!");
    }
  }

  MInt nonBlockingComm = static_cast<MInt>(m_nonBlockingComm);

  if(nonBlockingComm != 0) {
    TERMM(1, "Due to the solver reconstruction, this exchange method is not available at the moment");
    m_exchangeMethod = &LbSolver::exchangeLbNB;
  } else {
    m_exchangeMethod = &LbSolver::exchangeLb;
  }

  switch(solverMethodEnum) {
    case MAIA_LATTICE_BGK:
      m_initializeMethod = &LbSolver::initializeLatticeBgk;
      m_solutionStepMethod = &LbSolver::bgki_collision_step;
      break;
    case MAIA_LATTICE_BGK_GUO_FORCING:
      m_initializeMethod = &LbSolver::initializeLatticeBgk;
      m_solutionStepMethod = &LbSolver::bgki_collision_step_Guo_forcing;
      break;
    case MAIA_LATTICE_BGK_INIT:
      m_initializeMethod = &LbSolver::initializeLatticeBgk;
      m_solutionStepMethod = &LbSolver::bgki_init_collision_step;
      break;
    case MAIA_LATTICE_BGKC:
      m_initializeMethod = &LbSolver::initializeLatticeBgk;
      m_solutionStepMethod = &LbSolver::bgkc_collision_step;
      m_updateMacroscopicLocation =
          (m_updateMacroscopicLocation == INCOLLISION) ? PRECOLLISION : m_updateMacroscopicLocation;
      m_isCompressible = true; // TODO labels:LB move to sysEqn-like class
      break;
    case MAIA_LATTICE_BGKI_SMAGORINSKY:
      m_initializeMethod = &LbSolver::initializeLatticeBgk;
      m_solutionStepMethod = &LbSolver::bgki_smagorinsky_collision_step;
      break;
    case MAIA_LATTICE_BGKI_SMAGORINSKY2:
      m_initializeMethod = &LbSolver::initializeLatticeBgk;
      m_solutionStepMethod = &LbSolver::bgki_smagorinsky_collision_step2;
      break;
    case MAIA_LATTICE_BGKI_SMAGO_WALL:
      m_initializeMethod = &LbSolver::initializeLatticeBgk;
      m_solutionStepMethod = &LbSolver::bgki_smago_wall_collision_step;
      break;
    case MAIA_LATTICE_BGKI_DYNAMIC_SMAGO:
      m_initializeMethod = &LbSolver::initializeLatticeBgk;
      m_solutionStepMethod = &LbSolver::bgki_dynamic_smago_collision_step;
      break;
    case MAIA_LATTICE_RBGK_DYNAMIC_SMAGO:
      m_initializeMethod = &LbSolver::initializeLatticeBgk;
      m_solutionStepMethod = &LbSolver::rbgk_dynamic_smago_collision_step;
      break;
    case MAIA_LATTICE_BGKI_EULER_2D:
      m_initializeMethod = &LbSolver::initializeLatticeBgk;
      m_solutionStepMethod = &LbSolver::bgki_euler_collision_step;
      break;
    case MAIA_LATTICE_BGK_THERMAL:
      m_initializeMethod = &LbSolver::initializeLatticeBgk;
      m_solutionStepMethod = &LbSolver::bgki_thermal_collision_step;
      m_isCompressible = true; // TODO labels:LB move to sysEqn-like class
      break;
    case MAIA_LATTICE_BGK_INNERENERGY:
      m_initializeMethod = &LbSolver::initializeLatticeBgk;
      m_solutionStepMethod = &LbSolver::bgki_innerEnergy_collision_step;
      m_isCompressible = true; // TODO labels:LB move to sysEqn-like class
      break;
    case MAIA_LATTICE_BGK_TOTALENERGY:
      m_initializeMethod = &LbSolver::initializeLatticeBgk;
      m_solutionStepMethod = &LbSolver::bgki_totalEnergy_collision_step;
      m_isCompressible = true; // TODO labels:LB move to sysEqn-like class
      break;
    case MAIA_LATTICE_BGK_TRANSPORT:
      m_initializeMethod = &LbSolver::initializeLatticeBgk;
      m_solutionStepMethod = &LbSolver::bgkc_transport_collision_step;
      m_isCompressible = true; // TODO labels:LB move to sysEqn-like class
      break;
    case MAIA_LATTICE_BGK_THERMAL_TRANSPORT:
      m_initializeMethod = &LbSolver::initializeLatticeBgk;
      m_solutionStepMethod = &LbSolver::bgkc_thermal_transport_collision_step;
      m_isCompressible = true; // TODO labels:LB move to sysEqn-like class
      break;
    case MAIA_LATTICE_BGK_INNERENERGY_TRANSPORT:
      m_initializeMethod = &LbSolver::initializeLatticeBgk;
      m_solutionStepMethod = &LbSolver::bgkc_innerenergy_transport_collision_step;
      m_isCompressible = true; // TODO labels:LB move to sysEqn-like class
      break;
    case MAIA_LATTICE_BGK_TOTALENERGY_TRANSPORT:
      m_initializeMethod = &LbSolver::initializeLatticeBgk;
      m_solutionStepMethod = &LbSolver::bgkc_totalenergy_transport_collision_step;
      m_isCompressible = true; // TODO labels:LB move to sysEqn-like class
      break;
    case MAIA_LATTICE_RBGK:
      m_initializeMethod = &LbSolver::initializeLatticeBgk;
      m_solutionStepMethod = &LbSolver::rbgk_collision_step;
      break;
    case MAIA_LATTICE_RBGK_SMAGORINSKY:
      m_initializeMethod = &LbSolver::initializeLatticeBgk;
      m_solutionStepMethod = &LbSolver::rbgk_smagorinsky_collision_step;
      break;
    case MAIA_LATTICE_MRT:
      m_initializeMethod = &LbSolver::initializeLatticeBgk;
      m_solutionStepMethod = &LbSolver::mrt_collision_step;
      m_updateMacroscopicLocation =
          (m_updateMacroscopicLocation == INCOLLISION) ? PRECOLLISION : m_updateMacroscopicLocation;
      m_isCompressible = true; // TODO labels:LB move to sysEqn-like class
      break;
    case MAIA_LATTICE_MRT2:
      m_initializeMethod = &LbSolver::initializeLatticeBgk;
      m_solutionStepMethod = &LbSolver::mrt2_collision_step;
      break;
    case MAIA_LATTICE_CLB:
      m_initializeMethod = &LbSolver::initializeLatticeBgk;
      m_solutionStepMethod = &LbSolver::clb_collision_step;
      break;
    case MAIA_LATTICE_CLB_SMAGORINSKY:
      m_initializeMethod = &LbSolver::initializeLatticeBgk;
      m_solutionStepMethod = &LbSolver::clb_smagorinsky_collision_step;
      break;
    case MAIA_LATTICE_MRT_SMAGORINSKY:
      m_initializeMethod = &LbSolver::initializeLatticeBgk;
      m_solutionStepMethod = &LbSolver::mrt_smagorinsky_collision_step;
      break;
    case MAIA_LATTICE_CUMULANT:
      m_initializeMethod = &LbSolver::initializeLatticeBgk;
      m_solutionStepMethod = &LbSolver::cumulant_collision_step;
      m_updateMacroscopicLocation =
          (m_updateMacroscopicLocation == INCOLLISION) ? PRECOLLISION : m_updateMacroscopicLocation;
      m_isCompressible = true; // TODO labels:LB move to sysEqn-like class
      break;
    default: {
      TERMM(1, "Unknown LB method: Exiting!");
    }
  }

  // If m_restartFile is set, the function pointer m_initializeMethod will be overriden
  // with the loadRestartFile, because the distributions and variables will
  // be taken from the restart file. If one of the properties initRestart or initFromCoarse
  // is set, the funktion restartInitLb will be used instead of loadRestartFile, since the
  // distributions have to be calculated.
  if(m_restartFile || m_initFromRestartFile) {
    m_initMethod = "FROM_RESTART_FILE";
    m_initializeMethod = &LbSolver::loadRestartFile;
    if(m_initRestart || m_initFromCoarse) {
      m_initializeMethod = &LbSolver::restartInitLb;
    }
  }
  if(Context::propertyExists("UrmsInit") || Context::propertyExists("ReLambda")) {
    getReLambdaAndUrmsInit();

    if(m_innerBandWidth != nullptr) mDeallocate(m_innerBandWidth);
    mAlloc(m_innerBandWidth, grid().maxRefinementLevel(), "m_innerBandWidth", F0, AT_);
    if(m_outerBandWidth != nullptr) mDeallocate(m_outerBandWidth);
    mAlloc(m_outerBandWidth, grid().maxRefinementLevel(), "m_outerBandWidth", F0, AT_);
    if(m_bandWidth != nullptr) mDeallocate(m_bandWidth);
    mAlloc(m_bandWidth, grid().maxRefinementLevel(), "m_bandWidth", 0, AT_);

    /*! \page propertyPage1
    \section mbBandwidth
    <code>MFloat FvMbSolverXD::distFac </code>\n
    default = {18.0, 9.0}\n
    Sets the distance factor which is used to calculate the inner (distFac[0]) and outer (distFac[1]) bandwidth.\n
    Possible values are:
    <ul>
      <li>Any floating point values.</li>
    </ul>
    Keywords: <i>MOVING BOUNDARY</i>
  */

    MFloat distFac[2] = {18.0, 9.0};
    for(MInt i = 0; i < 2; i++) {
      distFac[i] = Context::getSolverProperty<MFloat>("mbBandWidth", m_solverId, AT_, &distFac[i], i);
    }
    m_outerBandWidth[grid().maxRefinementLevel() - 1] = distFac[0] * c_cellLengthAtLevel(grid().maxRefinementLevel());
    m_bandWidth[grid().maxRefinementLevel() - 1] = distFac[0];
    for(MInt i = grid().maxRefinementLevel() - 2; i >= 0; i--) {
      m_outerBandWidth[i] = m_outerBandWidth[i + 1] + (distFac[1] * c_cellLengthAtLevel(i + 1));
      m_bandWidth[i] = (m_bandWidth[i + 1] / 2) + 1 + distFac[1];
    }
    for(MInt i = 0; i < grid().maxRefinementLevel(); i++) {
      m_innerBandWidth[i] = -m_outerBandWidth[i];
      m_log << "bandwidth level " << i << ": " << m_innerBandWidth[i] << " " << m_outerBandWidth[i] << endl;
    }
  }
  // Initialise a txt-file for the residual
  // Only process 0 should create a txt.file
  if(domainId() == 0) {
    // Create a txt file
    m_resFileName = "maia_res";
    m_resFileName += ".txt";
    mRes.open(m_resFileName, std::ofstream::app);

    // Print header
    mRes << "This is the residual file of the LB solver." << std::endl;
    mRes << std::endl;
    mRes << "#-----------------------------------------------------------------" << std::endl;
    mRes << "#---SOLVER INFORMATION " << std::endl;
    mRes << std::endl;
    mRes << "Number of ranks:			" << noDomains() << std::endl;
    mRes << "Used solver method:			" << solverMethod() << std::endl;
    mRes << "Used initialization method: 		" << m_initMethod << std::endl;
    mRes << "Number dimensions: 			" << nDim << std::endl;
    mRes << "Number of distributions: 		" << m_noDistributions << std::endl;
    mRes << "Minimal spatial refinement level:	" << minLevel() << std::endl;
    mRes << "Maximal spatial refinement level:	" << maxLevel() << std::endl;
    mRes << std::endl;
    mRes << "#-----------------------------------------------------------------" << std::endl;
    mRes << "#---RESIDUAL INFORMATION " << std::endl;
    mRes << std::endl;
    mRes << "Residual intervall:			" << m_residualInterval << std::endl;
    mRes << std::endl;
    mRes << "#-----------------------------------------------------------------" << std::endl;
    mRes << "#---FLOW PARAMETERS " << std::endl;
    mRes << std::endl;
    mRes << "Mach number in the simulation:		" << m_Ma << std::endl;
    mRes << "Reynolds number in the simulation:	" << m_Re << std::endl;
    mRes << std::endl;
    mRes << "#-----------------------------------------------------------------" << std::endl;
    mRes << std::endl << std::endl;
    mRes << "Calculated resudial for each interval: " << std::endl;

    mRes.close();
  }
}

template <MInt nDim>
LbSolver<nDim>::~LbSolver() {
  TRACE();

  resetComm();

  mDeallocate(m_residual);
  mDeallocate(m_tmpResidual);
  mDeallocate(m_tmpResidualLvl);
  mDeallocate(m_maxResId);
  mDeallocate(m_rescoordinates);

  mDeallocate(m_Fext);

  mDeallocate(m_EELiquid.Fg);

  mDeallocate(m_momentumFlux);
  mDeallocate(m_MijMij);
  mDeallocate(m_MijLij);

  mDeallocate(m_cellLength);

  if(m_interfaceChildren.size() != 0) {
    m_interfaceChildren.clear();
  }
  if(m_interfaceParents.size() != 0) {
    m_interfaceParents.clear();
  }
  mDeallocate(m_activeCellList);

  RECORD_TIMER_STOP(m_t.solver);

  averageTimer();
}

/** \brief  Initialize all timer groups, timer, and sub timer required by the LB solver
 *  \author Miro Gondrum
 *  \date   18.02.2022
 */
template <MInt nDim>
void LbSolver<nDim>::initTimer() {
  TRACE();

  NEW_TIMER_GROUP(tg_solver, "LB Solver (solverId=" + std::to_string(m_solverId) + ")");
  NEW_TIMER_NOCREATE(m_t.solver, "complete solver", tg_solver);

  NEW_SUB_TIMER_NOCREATE(m_t.initSolver, "init solver", m_t.solver);
  NEW_SUB_TIMER_NOCREATE(m_t.solutionStep, "SolutionStep", m_t.solver);
  NEW_SUB_TIMER_NOCREATE(m_t.collision, "Collision", m_t.solutionStep);
  NEW_SUB_TIMER_NOCREATE(m_t.collisionBC, "CollisionBC", m_t.solutionStep);
  NEW_SUB_TIMER_NOCREATE(m_t.propagation, "Propagation", m_t.solutionStep);
  NEW_SUB_TIMER_NOCREATE(m_t.propagationBC, "PropagationBC", m_t.solutionStep);
  NEW_SUB_TIMER_NOCREATE(m_t.exchange, "Exchange", m_t.solutionStep);
  NEW_SUB_TIMER_NOCREATE(m_t.packing, "Packing", m_t.exchange);
  NEW_SUB_TIMER_NOCREATE(m_t.unpacking, "Unpacking", m_t.exchange);
  NEW_SUB_TIMER_NOCREATE(m_t.communication, "Communication", m_t.exchange);
  NEW_SUB_TIMER_NOCREATE(m_t.residual, "Residuum", m_t.solutionStep);
  NEW_SUB_TIMER_NOCREATE(m_t.srcTerms, "SourceTerms", m_t.solutionStep);

  NEW_SUB_TIMER_NOCREATE(m_t.findG0Cells, "Find G0 cells", m_t.solver);
  NEW_SUB_TIMER_NOCREATE(m_t.resetListsMb, "Reset lists MB", m_t.findG0Cells);
  NEW_SUB_TIMER_NOCREATE(m_t.findG0Candidates, "Find G0 candidates", m_t.findG0Cells);
  NEW_SUB_TIMER_NOCREATE(m_t.geomNodal, "Geometry intersection compute Nodal", m_t.findG0Cells);
  NEW_SUB_TIMER_NOCREATE(m_t.geomExchange, "Geometry intersection exchange", m_t.findG0Cells);
  NEW_SUB_TIMER_NOCREATE(m_t.calcNodalValues, "Calc nodal values", m_t.findG0Cells);

  NEW_SUB_TIMER_NOCREATE(m_t.prepComm, "prepare communication", m_t.solver);

  NEW_SUB_TIMER_NOCREATE(m_t.fft, "FFT", m_t.solver);
}

/** \brief  Average all LB timer over participating ranks
 *  \author Miro Gondrum
 *  \date   18.02.2022
 */
template <MInt nDim>
void LbSolver<nDim>::averageTimer() {
  TRACE();
  if(!grid().isActive()) return;

  // Determine timer operation
  MString m_timerType = "max";
  m_timerType = Context::getSolverProperty<MString>("timerType", m_solverId, AT_, &m_timerType);

  // 0) map timer ids for safety
  std::vector<MInt> timerIds_;
  timerIds_.reserve(17);
  timerIds_.emplace_back(m_t.solver);
  timerIds_.emplace_back(m_t.initSolver);
  timerIds_.emplace_back(m_t.solutionStep);
  timerIds_.emplace_back(m_t.collision);
  timerIds_.emplace_back(m_t.collisionBC);
  timerIds_.emplace_back(m_t.propagation);
  timerIds_.emplace_back(m_t.propagationBC);
  timerIds_.emplace_back(m_t.srcTerms);
  timerIds_.emplace_back(m_t.exchange);
  timerIds_.emplace_back(m_t.residual);
  timerIds_.emplace_back(m_t.findG0Cells);
  timerIds_.emplace_back(m_t.resetListsMb);
  timerIds_.emplace_back(m_t.findG0Candidates);
  timerIds_.emplace_back(m_t.geomNodal);
  timerIds_.emplace_back(m_t.geomExchange);
  timerIds_.emplace_back(m_t.calcNodalValues);
  timerIds_.emplace_back(m_t.prepComm);
  timerIds_.emplace_back(m_t.fft);
  const MInt noTimers = timerIds_.size();

  // 1) fill buffer with local timer values
  std::vector<MFloat> timerValues_;
  timerValues_.reserve(noTimers);
  for(MInt i = 0; i < noTimers; i++) {
    timerValues_.emplace_back(RETURN_TIMER_TIME(timerIds_[i]));
  }

  // 2) collect values from all ranks
  if(m_timerType == "average") {
    MPI_Allreduce(MPI_IN_PLACE, timerValues_.data(), noTimers, maia::type_traits<MFloat>::mpiType(), MPI_SUM, mpiComm(),
                  AT_, "MPI_IN_PLACE", "timerValues_");
  } else {
    MPI_Allreduce(MPI_IN_PLACE, timerValues_.data(), noTimers, maia::type_traits<MFloat>::mpiType(), MPI_MAX, mpiComm(),
                  AT_, "MPI_IN_PLACE", "timerValues_");
  }

  // 3) perform averaging on timer and4) set new timer values
  if(m_timerType == "average") {
    const MInt noDomains_ = noDomains();
    for(MInt i = 0; i < noTimers; i++) {
      const MFloat meanValue = timerValues_[i] / noDomains_;
      SET_RECORD(timerIds_[i], meanValue);
    }
  } else {
    for(MInt i = 0; i < noTimers; i++) {
      SET_RECORD(timerIds_[i], timerValues_[i]);
    }
  }
}

template <MInt nDim>
void LbSolver<nDim>::printScalingVariables() {
  std::map<MString, MInt> scalingVars;

  scalingVars["noCells"] = a_noCells();
  scalingVars["noInternalCells"] = noInternalCells();
  scalingVars["noHaloCells"] = a_noCells() - grid().noInternalCells();
  scalingVars["noActiveCells"] = m_currentMaxNoCells;

  std::vector<MInt> recvBuffer;
  recvBuffer.resize(noDomains());
  for(const auto& [name, value] : scalingVars) {
    MPI_Gather(&value, 1, maia::type_traits<MInt>::mpiType(), recvBuffer.data(), 1, maia::type_traits<MInt>::mpiType(),
               0, mpiComm(), AT_, "send", "recv");

    if(!domainId()) {
      // average
      MFloat avgValue = 0;
      avgValue = std::accumulate(recvBuffer.begin(), recvBuffer.end(), avgValue);
      avgValue /= noDomains();
      std::cout << "Average " + name + " per Domain are: " << avgValue << std::endl;

      // maximum
      const MInt maxValue = *std::max_element(recvBuffer.begin(), recvBuffer.end());
      std::cout << "Max " + name + " per Domain are: " << maxValue << std::endl;

      // minimum
      const MInt minValue = *std::min_element(recvBuffer.begin(), recvBuffer.end());
      std::cout << "Min " + name + " per Domain are: " << minValue << std::endl;
    }
  }
}

/// Copy selected information from grid to cell collector to allow modifications and/or faster access.
///
/// \author Michael Schlottke-Lakemper (mic) <mic@aia.rwth-aachen.de>
/// \date 2018-09-18
template <MInt nDim>
void LbSolver<nDim>::updateCellCollectorFromGrid() {
  // Store level information in cell collector for faster access
  for(MInt i = 0; i < a_noCells(); i++) {
    m_cells.level(i) = c_level(i);
    ASSERT(a_level(i) = c_level(i), "");
    a_isHalo(i) = false;
    a_isWindow(i) = false;
  }
  for(MInt i = 0; i < noNeighborDomains(); i++) {
    for(MInt c = 0; c < noHaloCells(i); c++) {
      a_isHalo(haloCell(i, c)) = true;
    }
  }
  for(MInt i = 0; i < noNeighborDomains(); i++) {
    for(MInt c = 0; c < noWindowCells(i); c++) {
      a_isWindow(windowCell(i, c)) = true;
    }
  }
}


/** \brief Prepares the communication
 *
 * \author Andreas Lintermann
 * \date 05.04.2013, 03.08.2015
 *
 * The following is done in this function:
 *
 * Check if we want to do a reduced communication
 *  - if no, then run the previous code (run prepareCommunicationNormal)
 *  - if yes, the run the reduced code ( run prepareCommunicationReduced)
 *
 **/
template <MInt nDim>
void LbSolver<nDim>::prepareCommunication() {
  TRACE();

  RECORD_TIMER_START(m_t.prepComm);

  m_log << endl;
  m_log << "#########################################################################################################"
           "#############"
        << endl;
  m_log << "##                                                Communication                                          "
           "           ##"
        << endl;
  m_log << "#########################################################################################################"
           "#############"
        << endl;

  /*! \page propertyPage1
    \section reducedComm
    <code>MInt LbSolver::m_reducedComm</code>\n
    default = <code>0</code>\n\n
    This property defines if reduced communication should be activated.
    <ul>
    <li><code>0</code> (off)</li>
    <li><code>1</code> (on)</li>
    </ul>\n
    Keywords: <i>LATTICE BOLTZMANN</i>
  */
  m_reducedComm = false;
  m_reducedComm = Context::getSolverProperty<MBool>("reducedComm", m_solverId, AT_, &m_reducedComm);

  if(m_reducedComm) {
    // 1. set the right function pointers
    m_gatherMethod = &LbSolver::gatherReduced;
    m_scatterMethod = &LbSolver::scatterReduced;
    m_sendMethod = &LbSolver::sendReduced;
    m_receiveMethod = &LbSolver::receiveReduced;
    prepareCommunicationReduced();
  } else {
    // 1. Set the right function pointers
    m_gatherMethod = &LbSolver::gatherNormal;
    m_scatterMethod = &LbSolver::scatterNormal;
    m_sendMethod = &LbSolver::sendNormal;
    m_receiveMethod = &LbSolver::receiveNormal;
    prepareCommunicationNormal();
  }
  m_log << " lbsolver: domain " << domainId() << ", has " << a_noCells() << " cells" << endl;

  RECORD_TIMER_STOP(m_t.prepComm);
}

/** \brief Prepares the communication by allocating buffers etc. for normal communication
 *
 * \author Andreas Lintermann
 * \date 03.08.2015
 *
 *  1. set the right function pointers
 *
 *  2. set the data solver sizes and allocate space for the offsets:
 *     + depending on the LB method only a certain amount of data has to be transferred:
 *       - normal LBGK            : PPDFs
 *       - thermal LBGK           : PPDFs, thermal PPDFs
 *       - normal, refined LBGK   : PPDFs, PVs, old PVs
 *       - thermal, refined LBGK  : PPDFs, thermal PPDFs, PVs, oldPVs
 *     + the array m_dataBlockSizes holds the size of one dataset per cell
 *     + the array m_baseAddresses holds the base pointer of the datasets
 *     + the variable m_dataBlockSizeTotal holds the total size of data per cell to be transferred
 *
 *   3. allocate the offsets:
 *     + m_nghbrOffsetsWindow holds the offsets in the total
 *
 *   4. prints communication information
 **/
template <MInt nDim>
void LbSolver<nDim>::prepareCommunicationNormal() {
  TRACE();

  m_log << "  + Preparing normal communication ..." << endl << endl;

  // 2. Set the data solver sizes and allocate space for the offsets
  if(m_isRefined || m_isEELiquid) {
    // Transfer PVs, old PVs, thermal and normal PPDFs  : sum = 4 offsets
    if(m_isThermal && !m_isTransport) {
      m_noElementsTransfer = 4;
      mAlloc(m_dataBlockSizes, m_noElementsTransfer, "m_noDataBlockSizes", 0, AT_);
      mAlloc(m_baseAddresses, m_noElementsTransfer, "m_baseAddresses", AT_);

      m_dataBlockSizes[0] = m_noDistributions;
      m_dataBlockSizes[1] = m_noDistributions;
      m_dataBlockSizes[2] = nDim + 2;
      m_dataBlockSizes[3] = nDim + 2;

      m_baseAddresses[0] = &a_distribution(0, 0);
      m_baseAddresses[1] = &a_distributionThermal(0, 0);
      m_baseAddresses[2] = &a_variable(0, 0);
      m_baseAddresses[3] = &a_oldVariable(0, 0);
    } else if(m_isTransport && !m_isThermal) {
      m_noElementsTransfer = 4;
      mAlloc(m_dataBlockSizes, m_noElementsTransfer, "m_noDataBlockSizes", 0, AT_);
      mAlloc(m_baseAddresses, m_noElementsTransfer, "m_baseAddresses", AT_);

      m_dataBlockSizes[0] = m_noDistributions;
      m_dataBlockSizes[1] = m_noDistributions;
      m_dataBlockSizes[2] = nDim + 2;
      m_dataBlockSizes[3] = nDim + 2;

      m_baseAddresses[0] = &a_distribution(0, 0);
      m_baseAddresses[1] = &a_distributionTransport(0, 0);
      m_baseAddresses[2] = &a_variable(0, 0);
      m_baseAddresses[3] = &a_oldVariable(0, 0);
    } else if(m_isThermal && m_isTransport) {
      m_noElementsTransfer = 5;
      mAlloc(m_dataBlockSizes, m_noElementsTransfer, "m_noDataBlockSizes", 0, AT_);
      mAlloc(m_baseAddresses, m_noElementsTransfer, "m_baseAddresses", AT_);

      m_dataBlockSizes[0] = m_noDistributions;
      m_dataBlockSizes[1] = m_noDistributions;
      m_dataBlockSizes[2] = m_noDistributions;
      m_dataBlockSizes[3] = nDim + 3;
      m_dataBlockSizes[4] = nDim + 3;

      m_baseAddresses[0] = &a_distribution(0, 0);
      m_baseAddresses[1] = &a_distributionThermal(0, 0);
      m_baseAddresses[2] = &a_distributionTransport(0, 0);
      m_baseAddresses[3] = &a_variable(0, 0);
      m_baseAddresses[4] = &a_oldVariable(0, 0);
    }
    // Transfer PVs, oldPVs and normal PPDFs           : sum = 3 offsets
    else {
      m_noElementsTransfer = 3;
      mAlloc(m_dataBlockSizes, m_noElementsTransfer, "m_noDataBlockSizes", 0, AT_);
      mAlloc(m_baseAddresses, m_noElementsTransfer, "m_baseAddresses", AT_);

      m_dataBlockSizes[0] = m_noDistributions;
      m_dataBlockSizes[1] = nDim + 1;
      m_dataBlockSizes[2] = nDim + 1;

      m_baseAddresses[0] = &a_distribution(0, 0);
      m_baseAddresses[1] = &a_variable(0, 0);
      m_baseAddresses[2] = &a_oldVariable(0, 0);
    }
  } else {
    // Transfer thermal and normal PPDFs        : sum = 2 offsets
    if(m_isThermal && !m_isTransport) {
      m_noElementsTransfer = 2;
      mAlloc(m_dataBlockSizes, m_noElementsTransfer, "m_noDataBlockSizes", 0, AT_);
      mAlloc(m_baseAddresses, m_noElementsTransfer, "m_baseAddresses", AT_);

      m_dataBlockSizes[0] = m_noDistributions;
      m_dataBlockSizes[1] = m_noDistributions;

      m_baseAddresses[0] = &a_distribution(0, 0);
      m_baseAddresses[1] = &a_distributionThermal(0, 0);
    } else if(m_isTransport && !m_isThermal) {
      m_noElementsTransfer = 2;
      mAlloc(m_dataBlockSizes, m_noElementsTransfer, "m_noDataBlockSizes", 0, AT_);
      mAlloc(m_baseAddresses, m_noElementsTransfer, "m_baseAddresses", AT_);

      m_dataBlockSizes[0] = m_noDistributions;
      m_dataBlockSizes[1] = m_noDistributions;

      m_baseAddresses[0] = &a_distribution(0, 0);
      m_baseAddresses[1] = &a_distributionTransport(0, 0);

    } else if(m_isTransport && m_isThermal) {
      m_noElementsTransfer = 3;
      mAlloc(m_dataBlockSizes, m_noElementsTransfer, "m_noDataBlockSizes", 0, AT_);
      mAlloc(m_baseAddresses, m_noElementsTransfer, "m_baseAddresses", AT_);

      m_dataBlockSizes[0] = m_noDistributions;
      m_dataBlockSizes[1] = m_noDistributions;
      m_dataBlockSizes[2] = m_noDistributions;

      m_baseAddresses[0] = &a_distribution(0, 0);
      m_baseAddresses[1] = &a_distributionThermal(0, 0);
      m_baseAddresses[2] = &a_distributionTransport(0, 0);

    }
    // Transfer normal PPDFs                    : sum = 1 offset
    else {
      // This is needed if gradients are calculated (e.g. for LODI equations)
      m_noElementsTransfer = 2;
      mAlloc(m_dataBlockSizes, m_noElementsTransfer, "m_noDataBlockSizes", 0, AT_);
      mAlloc(m_baseAddresses, m_noElementsTransfer, "m_baseAddresses", AT_);

      m_dataBlockSizes[0] = m_noDistributions;
      m_dataBlockSizes[1] = nDim + 1;

      m_baseAddresses[0] = &a_distribution(0, 0);
      m_baseAddresses[1] = &a_variable(0, 0);
    }
  }

  // calculate the total size of the data to be transferred per cell
  m_dataBlockSizeTotal = 0;
  for(MInt i = 0; i < m_noElementsTransfer; i++)
    m_dataBlockSizeTotal += m_dataBlockSizes[i];

  // 3. Allocate the offsets
  if(noNeighborDomains() > 0) {
    mAlloc(m_nghbrOffsetsWindow, noNeighborDomains(), m_noElementsTransfer, "m_nghbrOffsetsWindow", 0, AT_);
    mAlloc(m_nghbrOffsetsHalo, noNeighborDomains(), m_noElementsTransfer, "m_nghbrOffsetsHalo", 0, AT_);
  }
  for(MInt n = 0; n < noNeighborDomains(); n++) {
    for(MInt v = 1; v < m_noElementsTransfer; v++) {
      m_nghbrOffsetsWindow[n][v] = m_nghbrOffsetsWindow[n][v - 1] + noWindowCells(n) * m_dataBlockSizes[v - 1];
      m_nghbrOffsetsHalo[n][v] = m_nghbrOffsetsHalo[n][v - 1] + noHaloCells(n) * m_dataBlockSizes[v - 1];
    }
  }

  if(noNeighborDomains() > 0) {
    // Allocate the buffers
    ScratchSpace<MInt> haloCellsCnt(noNeighborDomains(), AT_, "noHaloCells");
    ScratchSpace<MInt> windowCellsCnt(noNeighborDomains(), AT_, "noWindowCells");
    for(MInt d = 0; d < noNeighborDomains(); d++) {
      haloCellsCnt[d] = noHaloCells(d);
      windowCellsCnt[d] = noWindowCells(d);
    }
    mAlloc(m_sendBuffers, noNeighborDomains(), &windowCellsCnt[0], m_dataBlockSizeTotal, "m_sendBuffers", AT_);
    mAlloc(m_receiveBuffers, noNeighborDomains(), &haloCellsCnt[0], m_dataBlockSizeTotal, "m_receiveBuffers", AT_);

    // blocking communication
    if(!m_nonBlockingComm) {
      mAlloc(mpi_request, noNeighborDomains(), "mpi_request", AT_);
    }
    // non-blocking communication
    else {
      mAlloc(mpi_requestS, noNeighborDomains(), "mpi_requestS", AT_);
      mAlloc(mpi_requestR, noNeighborDomains(), "mpi_requestR", AT_);
    }
  }

  // 4. Print the information
  printCommunicationMethod();
}

/** \brief Prepares the communication by allocating buffers etc. for normal communication
 *
 * \author Andreas Lintermann, Moritz Waldmann
 * \date 20.11.2021
 *
 *   1. set the right function pointers
 *   2. count the elements per neighbor domain that we need to transfer
 *   3. fill the arrays depending on the number of elements
 *   4. allocate the buffers
 *   5. print the communication information
 *
 **/
template <MInt nDim>
void LbSolver<nDim>::prepareCommunicationReduced() {
  TRACE();

  m_log << "  + Preparing reduced communication ..." << endl << endl;

  // blocking communication
  if(noNeighborDomains() > 0) {
    if(!m_nonBlockingComm) {
      mAlloc(mpi_request, noNeighborDomains(), "mpi_request", AT_);
    } else {
      mAlloc(mpi_requestS, noNeighborDomains(), "mpi_requestS", AT_);
      mAlloc(mpi_requestR, noNeighborDomains(), "mpi_requestR", AT_);
    }

    ScratchSpace<MInt> noHalosPerDomain(noNeighborDomains(), AT_, "noHalosPerDomain");
    ScratchSpace<MInt> noWindowsPerDomain(noNeighborDomains(), AT_, "noWindowsPerDomain");
    for(MInt n = 0; n < noNeighborDomains(); n++) {
      noHalosPerDomain[n] = noHaloCells(n) * (m_noDistributions - 1) + 1;
      noWindowsPerDomain[n] = noWindowCells(n) * (m_noDistributions - 1) + 1;
    }

    mAlloc(m_noWindowDistDataPerDomain, noNeighborDomains(), "m_noWindowDistDataPerDomain", 0, AT_);
    mAlloc(m_noHaloDistDataPerDomain, noNeighborDomains(), "m_noHaloDistDataPerDomain", 0, AT_);

    mAlloc(m_windowDistsForExchange, noNeighborDomains(), &noWindowsPerDomain[0], "m_windowDistsForExchange", AT_);
    mAlloc(m_haloDistsForExchange, noNeighborDomains(), &noHalosPerDomain[0], "m_haloDistsForExchange", AT_);
  }

  mAlloc(m_needsFurtherExchange, a_noCells(), "needsFurtherExchange", 0, AT_);
  markCellsForAdditionalComm();

  if(m_isRefined) {
    // 2. count the elements per neighbor domain that we need to transfer
    m_log << "    - counting PPDFs per neighbor" << endl;

    for(MInt n = 0; n < noNeighborDomains(); n++) {
      MInt haloStart = domainOffset(neighborDomain(n));
      MInt haloEnd = domainOffset(neighborDomain(n) + 1);

      MInt cnt = 0;
      for(MInt j = 0; j < noWindowCells(n); j++) {
        if(m_needsFurtherExchange[windowCell(n, j)] > 0) {
          for(MInt d = 0; d < m_noDistributions - 1; d++) {
            m_noWindowDistDataPerDomain[n]++;
            m_windowDistsForExchange[n][j * (m_noDistributions - 1) + d] = cnt;
            cnt++;
          }
        } else {
          for(MInt d = 0; d < m_noDistributions - 1; d++) {
            m_windowDistsForExchange[n][j * (m_noDistributions - 1) + d] = -1;
            MInt ngh = c_neighborId(windowCell(n, j), d);
            if(ngh != -1 && a_isHalo(ngh)) {
              MInt global = c_globalId(ngh);
              if(global >= haloStart && global < haloEnd) {
                m_noWindowDistDataPerDomain[n]++;
                m_windowDistsForExchange[n][j * (m_noDistributions - 1) + d] = cnt;
                cnt++;
              }
            }
          }
        }
      }
      cnt = 0;
      for(MInt j = 0; j < noHaloCells(n); j++) {
        if(m_needsFurtherExchange[haloCell(n, j)] > 0) {
          for(MInt d = 0; d < m_noDistributions - 1; d++) {
            m_noHaloDistDataPerDomain[n]++;
            m_haloDistsForExchange[n][j * (m_noDistributions - 1) + d] = cnt;
            cnt++;
          }
        } else {
          for(MInt d = 0; d < m_noDistributions - 1; d++) {
            m_haloDistsForExchange[n][j * (m_noDistributions - 1) + d] = -1;
            MInt ngh = c_neighborId(haloCell(n, j), d);
            if(ngh != -1 && a_isWindow(ngh)) {
              m_noHaloDistDataPerDomain[n]++;
              m_haloDistsForExchange[n][j * (m_noDistributions - 1) + d] = cnt;
              cnt++;
            }
          }
        }
      }

      m_log << "      * neighbor " << neighborDomain(n) << ":" << endl;
      m_log << "        = windows :" << m_noWindowDistDataPerDomain[n] << endl;
      m_log << "        = halos :" << m_noHaloDistDataPerDomain[n] << endl;
    }

    // 3. fill the arrays depending on the number of elements
    if(m_isThermal && !m_isTransport) {
      m_noDistsTransfer = 2;
      m_noVarsTransfer = 2;
      m_noElementsTransfer = m_noDistsTransfer + m_noVarsTransfer;
    } else if(m_isTransport && !m_isThermal) {
      m_noDistsTransfer = 2;
      m_noVarsTransfer = 2;
      m_noElementsTransfer = m_noDistsTransfer + m_noVarsTransfer;
    } else if(m_isTransport && m_isThermal) {
      m_noDistsTransfer = 3;
      m_noVarsTransfer = 2;
      m_noElementsTransfer = m_noDistsTransfer + m_noVarsTransfer;
    } else {
      m_noDistsTransfer = 1;
      m_noVarsTransfer = 2;
      m_noElementsTransfer = m_noDistsTransfer + m_noVarsTransfer;
    }
  } else {
    // 2. count the elements per neighbor domain that we need to transfer
    m_log << "    - counting PPDFs per neighbor" << endl;

    for(MInt n = 0; n < noNeighborDomains(); n++) {
      MInt haloStart = domainOffset(neighborDomain(n));
      MInt haloEnd = domainOffset(neighborDomain(n) + 1);

      MInt cnt = 0;
      for(MInt j = 0; j < noWindowCells(n); j++) {
        if(!c_isLeafCell(windowCell(n, j))) {
          for(MInt d = 0; d < m_noDistributions - 1; d++) {
            m_windowDistsForExchange[n][j * (m_noDistributions - 1) + d] = -1;
          }
        } else {
          if(m_needsFurtherExchange[windowCell(n, j)] > 0) {
            for(MInt d = 0; d < m_noDistributions; d++) {
              m_noWindowDistDataPerDomain[n]++;
              m_windowDistsForExchange[n][j * (m_noDistributions - 1) + d] = cnt;
              cnt++;
            }
          } else {
            for(MInt d = 0; d < m_noDistributions - 1; d++) {
              m_windowDistsForExchange[n][j * (m_noDistributions - 1) + d] = -1;
              MInt ngh = c_neighborId(windowCell(n, j), d);
              if(ngh != -1 && a_isHalo(ngh)) {
                MInt global = c_globalId(ngh);
                if(global >= haloStart && global < haloEnd) {
                  m_noWindowDistDataPerDomain[n]++;
                  m_windowDistsForExchange[n][j * (m_noDistributions - 1) + d] = cnt;
                  cnt++;
                }
              }
            }
          }
        }
      }
      cnt = 0;
      for(MInt j = 0; j < noHaloCells(n); j++) {
        if(!c_isLeafCell(haloCell(n, j))) {
          for(MInt d = 0; d < m_noDistributions - 1; d++) {
            m_haloDistsForExchange[n][j * (m_noDistributions - 1) + d] = -1;
          }
        } else {
          if(m_needsFurtherExchange[haloCell(n, j)] > 0) {
            for(MInt d = 0; d < m_noDistributions; d++) {
              m_noHaloDistDataPerDomain[n]++;
              m_haloDistsForExchange[n][j * (m_noDistributions - 1) + d] = cnt;
              cnt++;
            }
          } else {
            for(MInt d = 0; d < m_noDistributions - 1; d++) {
              m_haloDistsForExchange[n][j * (m_noDistributions - 1) + d] = -1;
              MInt ngh = c_neighborId(haloCell(n, j), d);
              if(ngh != -1 && a_isWindow(ngh)) {
                m_noHaloDistDataPerDomain[n]++;
                m_haloDistsForExchange[n][j * (m_noDistributions - 1) + d] = cnt;
                cnt++;
              }
            }
          }
        }
      }

      m_log << "      * neighbor " << neighborDomain(n) << ":" << endl;
      m_log << "        = windows :" << m_noWindowDistDataPerDomain[n] << endl;
      m_log << "        = halos :" << m_noHaloDistDataPerDomain[n] << endl;
    }

    // 3. fill the arrays depending on the number of elements
    if(m_isThermal && !m_isTransport) {
      m_noDistsTransfer = 2;
      m_noVarsTransfer = 0;
      m_noElementsTransfer = m_noDistsTransfer + m_noVarsTransfer;
    } else if(m_isTransport && !m_isThermal) {
      m_noDistsTransfer = 2;
      m_noVarsTransfer = 0;
      m_noElementsTransfer = m_noDistsTransfer + m_noVarsTransfer;
    } else if(m_isTransport && m_isThermal) {
      m_noDistsTransfer = 3;
      m_noVarsTransfer = 0;
      m_noElementsTransfer = m_noDistsTransfer + m_noVarsTransfer;
    } else {
      m_noDistsTransfer = 1;
      m_noVarsTransfer = 1;
      m_noElementsTransfer = m_noDistsTransfer + m_noVarsTransfer;
    }
  }

  // 4. allocate the buffers
  if(noNeighborDomains() > 0) {
    ScratchSpace<MInt> haloDataPerDomain(noNeighborDomains(), AT_, "noHaloCells");
    ScratchSpace<MInt> windowDataPerDomain(noNeighborDomains(), AT_, "noWindowCells");

    for(MInt n = 0; n < noNeighborDomains(); n++) {
      windowDataPerDomain(n) =
          m_noWindowDistDataPerDomain[n] * m_noDistsTransfer + noWindowCells(n) * m_noVariables * m_noVarsTransfer;
      haloDataPerDomain(n) =
          m_noHaloDistDataPerDomain[n] * m_noDistsTransfer + noHaloCells(n) * m_noVariables * m_noVarsTransfer;
    }

    mAlloc(m_sendBuffers, noNeighborDomains(), &windowDataPerDomain[0], "m_sendBuffers", AT_);
    mAlloc(m_receiveBuffers, noNeighborDomains(), &haloDataPerDomain[0], "m_receiveBuffers", AT_);
  }

  // 5. print the communication information
  printCommunicationMethod();
}

/** \brief Marks cells which need additional communication compared to regular cells, for which the reduced
 *  communication is sufficient
 *
 * \author Moritz Waldmann
 * \date 20.11.2021
 *
 *  1. find interface cells, find boundary cells, find G0-cells
 *  2. mark the neighbors of these cells
 *  3. mark the neighbors neighbors of these cells
 *
 **/
template <MInt nDim>
void LbSolver<nDim>::markCellsForAdditionalComm() {
  TRACE();

  // All interface cells need further exchange, i.e. all distributions need to be exchanged
  MInt periodicSimulation = 0;
  for(MInt i = 0; i < a_noCells(); i++) {
    if(grid().isPeriodic(i)) {
      periodicSimulation = 1;
      break;
    }
  }
  MPI_Allreduce(MPI_IN_PLACE, &periodicSimulation, 1, MPI_INT, MPI_MAX, mpiComm(), AT_, "MPI_IN_PLACE",
                "periodicSimulation");

  if(periodicSimulation) {
    MInt noPeriodicHalos = 0;
    MIntScratchSpace noPeriodicHalosPerDomain(noNeighborDomains(), AT_, "noDomains");
    std::vector<MInt> myGlobIdsPeriHalos;
    for(MInt i = 0; i < a_noCells(); i++) {
      if(!a_isHalo(i)) continue;
      if(grid().isPeriodic(i)) {
        myGlobIdsPeriHalos.push_back(c_globalId(i));
        noPeriodicHalos++;
      }
    }

    for(MInt i = 0; i < noNeighborDomains(); i++) {
      MPI_Issend(&noPeriodicHalos, 1, MPI_INT, neighborDomain(i), 0, mpiComm(), &mpi_request[i], AT_,
                 "noPeriodicHalos");
    }
    MPI_Status status;
    for(MInt i = 0; i < noNeighborDomains(); i++) {
      MPI_Recv(&noPeriodicHalosPerDomain(i), 1, MPI_INT, neighborDomain(i), 0, mpiComm(), &status, AT_,
               "noPeriodicHalosPerDomain(i)");
    }
    for(MInt i = 0; i < noNeighborDomains(); i++) {
      MPI_Wait(&mpi_request[i], &status, AT_);
    }
    MInt** allGlobIdsPeriHalos{};
    mAlloc(allGlobIdsPeriHalos, noNeighborDomains(), &noPeriodicHalosPerDomain[0], "allGlobalIdsPeriHalos", AT_);

    for(MInt i = 0; i < noNeighborDomains(); i++) {
      MPI_Issend(myGlobIdsPeriHalos.data(), noPeriodicHalos, MPI_INT, neighborDomain(i), 0, mpiComm(), &mpi_request[i],
                 AT_, "noPeriodicHalos");
    }
    for(MInt i = 0; i < noNeighborDomains(); i++) {
      MPI_Recv(allGlobIdsPeriHalos[i], noPeriodicHalosPerDomain(i), MPI_INT, neighborDomain(i), 0, mpiComm(), &status,
               AT_, "noPeriodicHalosPerDomain(i)");
    }
    for(MInt i = 0; i < noNeighborDomains(); i++) {
      MPI_Wait(&mpi_request[i], &status, AT_);
    }

    for(MInt i = 0; i < a_noCells(); i++) {
      if(a_isHalo(i)) continue;
      MInt globalId = c_globalId(i);
      MBool found = false;
      for(MInt d = 0; d < noNeighborDomains(); d++) {
        for(MInt halo = 0; halo < noPeriodicHalosPerDomain(d); halo++) {
          if(allGlobIdsPeriHalos[d][halo] == globalId) {
            m_needsFurtherExchange[i] = 1;
            break;
          }
        }
        if(found) break;
      }
    }
    mDeallocate(allGlobIdsPeriHalos);
  }

  for(MInt i = 0; i < a_noCells(); i++) {
    if(a_isHalo(i)) continue;
    if(a_isInterfaceParent(i)) m_needsFurtherExchange[i] = 1;
    if(a_isInterfaceChild(i)) m_needsFurtherExchange[i] = 1;
    if(a_isBndryCell(i)) m_needsFurtherExchange[i] = 1;
    if(c_parentId(i) > 0) {
      if(a_isInterfaceParent(c_parentId(i))) m_needsFurtherExchange[i] = 1;
    }
    for(MInt n = 0; n < IPOW3[nDim] - 1; n++) {
      if(c_neighborId(i, n) > -1) {
        if(a_isActive(i) != a_isActive(c_neighborId(i, n))) {
          m_needsFurtherExchange[i] = 1;
          break;
        }
      }
    }
  }
  this->exchangeData(&(m_needsFurtherExchange[0]), 1);

  // Also the neighbors of the interface cells need special treatment
  for(MInt i = 0; i < a_noCells(); i++) {
    if(a_isHalo(i)) continue;
    if(m_needsFurtherExchange[i]) continue;
    for(MInt d = 0; d < m_noDistributions - 1; d++) {
      MInt nghbr = c_neighborId(i, d);
      if(nghbr > -1) {
        if(m_needsFurtherExchange[nghbr] == 1) m_needsFurtherExchange[i] = 2;
      }
    }
  }

  this->exchangeData(&(m_needsFurtherExchange[0]), 1);
  // Also the neighbors neighbors of the interface cells need special treatment
  for(MInt i = 0; i < a_noCells(); i++) {
    if(a_isHalo(i)) continue;
    if(m_needsFurtherExchange[i]) continue;
    for(MInt d = 0; d < m_noDistributions - 1; d++) {
      MInt nghbr = c_neighborId(i, d);
      if(nghbr > -1) {
        if(m_needsFurtherExchange[nghbr] == 2) m_needsFurtherExchange[i] = 3;
      }
    }
  }
  this->exchangeData(&(m_needsFurtherExchange[0]), 1);
}


/** \brief This function prints the communication setup to the logfile
 *
 * \author Andreas Lintermann
 * \date 08.04.2013
 *
 **/
template <MInt nDim>
void LbSolver<nDim>::printCommunicationMethod() {
  TRACE();

  if(!m_nonBlockingComm)
    m_log << "  + Communication method: blocking " << (m_reducedComm ? " (reduced)" : " (normal)") << endl << endl;
  else
    m_log << "  + Communication method: non-blocking" << (m_reducedComm ? " (reduced)" : " (normal)") << endl << endl;

  m_log << "  + Elements to be transferred: " << m_noElementsTransfer << " ( ";

  if(m_isRefined) {
    m_log << "PVs, oldPVs, PPDFs";
    if(m_isThermal)
      m_log << ", thPPDFs )" << endl;
    else
      m_log << " )" << endl;
  } else {
    m_log << "PPDFs";
    if(m_isThermal)
      m_log << ", thPPDFs )" << endl;
    else
      m_log << " )" << endl;
  }
  if(!m_reducedComm)
    for(MInt i = 0; i < m_noElementsTransfer; i++)
      m_log << "    - size [" << i << "]  : " << m_dataBlockSizes[i] << endl;

  if(!m_reducedComm) m_log << "    - total size: " << m_dataBlockSizeTotal << endl << endl;

  m_log << "  + Neighboring domains are:" << endl;

  for(MInt i = 0; i < noNeighborDomains(); i++) {
    m_log << "    - domain id " << neighborDomain(i) << ":" << endl;
    m_log << "      * window:              " << noWindowCells(i) << endl;
    m_log << "      * halo:                " << noHaloCells(i) << endl;
    m_log << "      * window data offsets: ";

    if(!m_reducedComm) {
      for(MInt v = 0; v < m_noElementsTransfer; v++)
        m_log << m_nghbrOffsetsWindow[i][v] << " ";
      m_log << endl;

      m_log << "      * halo data offsets:   ";
      for(MInt v = 0; v < m_noElementsTransfer; v++)
        m_log << m_nghbrOffsetsHalo[i][v] << " ";
      m_log << endl;
    }

    if(!m_reducedComm) {
      m_log << "      * total window buffer: " << (noWindowCells(i) * m_dataBlockSizeTotal * sizeof(MFloat)) << " bytes"
            << endl;
      m_log << "      * total halo buffer  : " << (noHaloCells(i) * m_dataBlockSizeTotal * sizeof(MFloat)) << " bytes"
            << endl;
    } else {
      m_log << "      * total window buffer: " << m_noWindowDistDataPerDomain[i] * sizeof(MFloat) << " bytes" << endl;
      m_log << "      * total halo buffer  : " << m_noHaloDistDataPerDomain[i] * sizeof(MFloat) << " bytes" << endl;
    }
  }

  m_log << endl << endl;
}

/** \brief Creates a subcommuincator
 *
 * \author Andreas Lintermann
 * \date 18.09.2015
 *
 * \param[in] ranks an array indicating which domains participate in the communicator
 * \param[in] noRanks size of the array
 * \param[in] comm the communicator
 *
 **/
template <MInt nDim>
void LbSolver<nDim>::createMPIComm(const MInt* ranks, MInt noRanks, MPI_Comm* comm) {
  TRACE();

  MIntScratchSpace ownerList(noRanks, AT_, "ownerList");
  for(MInt d = 0, l = 0; d < noDomains(); d++) {
    if(ranks[d] > 0) {
      ownerList[l] = d;
      l++;
    }
  }

  MPI_Group tmp_group;
  MPI_Group group;
  MPI_Comm_group(mpiComm(), &tmp_group, AT_, "tmp_group");
  MPI_Group_incl(tmp_group, noRanks, ownerList.getPointer(), &group, AT_);

  MPI_Comm_create(mpiComm(), group, comm, AT_, "comm");
}

/** \brief Write the boundary of a segment to a VTK
 *
 * \author Andreas Lintermann
 * \date 18.09.2015
 *
 * \param[in] bndVs array containing the points
 * \param[in] num size of the array
 *
 **/
template <MInt nDim>
inline void LbSolver<nDim>::writeSegmentBoundaryVTK(MFloat** bndVs, MInt num) {
  TRACE();

  stringstream name;
  name << "out/" << domainId() << "_boundary.vtk";
  ofstream st;
  st.open(name.str());
  st << "# vtk DataFile Version 3.0" << endl;
  st << "vtk output" << endl;
  st << "ASCII" << endl;
  st << "DATASET POLYDATA" << endl;
  st << "POINTS " << num << " float" << endl;

  for(MInt i = 0; i < num; i++) {
    for(MInt d = 0; d < 3; d++)
      st << bndVs[i][d] << " ";
    st << endl;
  }
  st << endl;
  st << "LINES " << 1 << " " << num + 2 << endl;
  st << num + 1 << " ";
  for(MInt i = 0; i < num; i++) {
    st << i << " ";
  }
  st << "0" << endl;

  st.close();
}

/** \brief  Update m_time.
 *  \author Miro Gondrum
 *  \date   19.05.2020
 *  dt corresponds to the real physical time step dt_phys, a reference length
 *  L_ref=1, and the real reference velocity u_inf such that:
 *  dt= dt_phys * u_inf / L_ref
 */
template <MInt nDim>
void LbSolver<nDim>::updateTime() {
  TRACE();
  constexpr MFloat L_ref = 1.0; // since it is used this way in FvSolver
  const MFloat dx = c_cellLengthAtLevel(maxLevel());
  m_dt = dx * m_Ma * LBCS / L_ref;
  m_time = getCurrentTimeStep() * m_dt;
}

/** \brief This function calculates the diamteter in cell units for parallel geometries
 *
 * \author Andreas Lintermann
 * \date 16.09.2015
 *
 * \param[in] segmentId the id of the segment to use for the calculation
 *
 **/
template <MInt nDim>
inline MFloat LbSolver<nDim>::calcCharLenParallelSplit(MInt segmentId) {
  TRACE();

  MInt own;
  MInt sumowners;
  MInt firstOwner;
  MIntScratchSpace owners(noDomains(), AT_, "owners");

  m_log << "      * segment owned by:          ";
  m_geometry->determineSegmentOwnership(segmentId, &own, &sumowners, &firstOwner, owners.getPointer());
  for(MInt d = 0; d < noDomains(); d++)
    if(owners[d] > 0) m_log << d << " ";

  m_log << endl;
  m_log << "      * sum of owners:             " << sumowners << endl;
  m_log << "      * root of communication:     " << firstOwner << endl;


  // my domain owns the whole segment
  if(sumowners == 1) {
    MFloat result = 0.0;
    if(own) result = calcCharLenAll(segmentId);

    MPI_Bcast(&result, 1, MPI_DOUBLE, firstOwner, mpiComm(), AT_, "result");

    return result;
  }

  // we need to share the geometry
  else {
    MFloat diameter;

    // build communicator
    MPI_Comm charComm;
    createMPIComm(owners.getPointer(), sumowners, &charComm);

    if(own) {
      // collect the triangles for testing first
      MInt offStart = m_geometry->m_segmentOffsets[segmentId];
      MInt offEnd = m_geometry->m_segmentOffsets[segmentId + 1];
      MInt numElements = offEnd - offStart;

      m_log << "      * number of local triangles: " << numElements << endl;
      m_log << "      * segment offsets:           " << m_geometry->m_segmentOffsets[segmentId] << " - "
            << m_geometry->m_segmentOffsets[segmentId + 1] << endl;

      // (normals + vertices)
      MInt noTriInfo = nDim * nDim;
      MIntScratchSpace myOriginalIds(numElements, AT_, "myOriginalIds");
      MFloatScratchSpace segTriangles(noTriInfo * numElements, AT_, "segTriangles");

      for(MInt t = offStart, i = 0, j = 0; t < offEnd; t++, i++) {
        myOriginalIds[i] = m_geometry->elements[t].m_originalId;
        for(MInt v = 0; v < nDim; v++)
          for(MInt d = 0; d < nDim; d++, j++)
            segTriangles[j] = m_geometry->elements[t].m_vertices[v][d];
      }

      MIntScratchSpace numElemPerCPU(sumowners, AT_, "numElemPerCPU");
      MPI_Allgather(&numElements, 1, MPI_INT, numElemPerCPU.getPointer(), 1, MPI_INT, charComm, AT_, "numElements",
                    "numElemPerCPU.getPointer()");
      MIntScratchSpace numTriInfoPerCPU(sumowners, AT_, "numTriInfoPerCPU");

      m_log << "      * triangles per domain:      ";
      MInt sumallelem = 0;
      for(MInt i = 0; i < sumowners; i++) {
        m_log << numElemPerCPU[i] << " ";
        sumallelem += numElemPerCPU[i];
        numTriInfoPerCPU[i] = numElemPerCPU[i] * noTriInfo;
      }
      m_log << endl;
      m_log << "      * sum of global triangles:   " << sumallelem << endl;


      MIntScratchSpace displOrig(sumowners, AT_, "displOrig");
      MIntScratchSpace displTris(sumowners, AT_, "displTris");
      displOrig[0] = 0;
      displTris[0] = 0;
      for(MInt d = 1; d < sumowners; d++) {
        displOrig[d] = displOrig[d - 1] + numElemPerCPU[d - 1];
        displTris[d] = displTris[d - 1] + numTriInfoPerCPU[d - 1];
      }

      MIntScratchSpace allOriginalIds(sumallelem, AT_, "allOriginalIds");
      MPI_Gatherv(myOriginalIds.getPointer(), numElements, MPI_INT, allOriginalIds.getPointer(),
                  numElemPerCPU.getPointer(), displOrig.getPointer(), MPI_INT, 0, charComm, AT_,
                  "myOriginalIds.getPointer()", "allOriginalIds.getPointer()");

      MFloatScratchSpace allSegTriangles(noTriInfo * sumallelem, AT_, "allSegTriangles");
      MPI_Gatherv(segTriangles.getPointer(), noTriInfo * numElements, MPI_DOUBLE, allSegTriangles.getPointer(),
                  numTriInfoPerCPU.getPointer(), displTris.getPointer(), MPI_DOUBLE, 0, charComm, AT_,
                  "segTriangles.getPointer()", "allSegTriangles.getPointer()");

      // now the firstOwner has received all relevant triangle information, calculate the characteristic length
      if(domainId() == firstOwner) {
        set<MInt> uniqueTriangles;
        for(MInt i = 0; i < sumallelem; i++)
          uniqueTriangles.insert(allOriginalIds[i]);

        MInt noUniqueTris = uniqueTriangles.size();
        m_log << "      * sum of unique triangles:   " << noUniqueTris << endl;

        MIntScratchSpace dbl(noUniqueTris, AT_, "dbl");
        for(MInt i = 0; i < noUniqueTris; i++)
          dbl[i] = 0;

        // this contains the offsets in the list of triangles which we want to use for the calculation
        MIntScratchSpace keepOffsets(noUniqueTris, AT_, "keepOffsets");
        for(MInt i = 0, j = 0; i < sumallelem; i++) {
          MInt dist = distance(uniqueTriangles.begin(), uniqueTriangles.find(allOriginalIds[i]));
          if(dist != noUniqueTris && dbl[dist] == 0) {
            keepOffsets[j] = i * noTriInfo;
            dbl[dist] = 1;
            j++;
          }
        }

        MInt num;
        MFloat** bndVs = m_geometry->GetBoundaryVertices(segmentId, allSegTriangles.getPointer(),
                                                         keepOffsets.getPointer(), noUniqueTris, &num);

        if(m_geometry->m_debugParGeom) writeSegmentBoundaryVTK(bndVs, num);

        // Calculate circumference
        MFloat circ = m_geometry->calcCircumference(bndVs, num);
        m_log << "      * circumference:             " << circ << endl;

        // Get size of surface
        MFloat size = m_geometry->GetBoundarySize(allSegTriangles.getPointer(), keepOffsets.getPointer(), noUniqueTris);
        m_log << "      * area:                      " << size << endl;

        // Get the hydraulic diameter
        MFloat hydraulic_diam = 4 * size / circ;
        m_log << "      * hydraulic diameter:        " << hydraulic_diam << endl;

        // Get bounding box
        std::array<MFloat, nDim * 2> bBox;
        m_geometry->getBoundingBox(bBox.data());

        // Find min and max in bounding box
        MFloat maxlength = 0.0;

        for(MInt i = 0; i < nDim; i++)
          if(fabs(bBox[i + nDim] - bBox[i]) > maxlength) maxlength = fabs(bBox[i + nDim] - bBox[i]);

        diameter = hydraulic_diam / (maxlength / FPOW2(maxLevel())) / reductionFactor();
      }
    }
    MPI_Bcast(&diameter, 1, MPI_DOUBLE, firstOwner, mpiComm(), AT_, "diameter");

    return diameter;
  }
  return 0.0;
}

/** \brief This function calculates the diamteter in cell units
 * \author Andreas Lintermann
 * \date 16.09.2015
 *
 * The diameter in cell units (on highest refinement level) is calculated by using the hydraulic diameter
 *
 * \f[h = \frac{4\cdot A}{C},\f]
 *
 * where \f$A\f$ is the area of a reference geometry element and \f$C\f$ is
 * its circumference.
 *
 * Then, the diameter is defined as follows:
 *
 * \f[d = \frac{h\cdot 2^{l_{max}}}{max\cdot r},\f]
 *
 * where \f$l_{max}\f$ is the maximal level in the grid, \f$max\f$ the maximal
 * length of the geometric bounding box and \f$r\f$ the reduction factor.
 *
 * \param[in] segmentId the id of the segment to use for the calculation
 *
 **/
template <MInt nDim>
inline MFloat LbSolver<nDim>::calcCharLenAll(MInt segmentId) {
  TRACE();

  // Get the boundary points of the segment
  MInt num = 0;
  MFloat** bndVs = m_geometry->GetBoundaryVertices(segmentId, nullptr, nullptr, 0, &num);

  // Calculate circumference
  MFloat circ = m_geometry->calcCircumference(bndVs, num);

  // Get size of surface
  MFloat size = m_geometry->GetBoundarySize(segmentId);

  // Get the hydraulic diameter
  MFloat hydraulic_diam = 4 * size / circ;

  // Get bounding box
  std::array<MFloat, nDim * 2> bBox;
  m_geometry->getBoundingBox(bBox.data());

  // Find min and max in bounding box
  MFloat maxlength = 0.0;

  for(MInt i = 0; i < nDim; i++)
    if(fabs(bBox[i + nDim] - bBox[i]) > maxlength) maxlength = fabs(bBox[i + nDim] - bBox[i]);

  m_log << "      * circumference:        " << circ << endl;
  m_log << "      * area:                 " << size << endl;
  m_log << "      * hydraulic diameter:   " << hydraulic_diam << endl;

  return hydraulic_diam / (maxlength / FPOW2(maxLevel())) / reductionFactor();
}

/** \brief Prepares the calculation of the characteristic length
 *
 * \author Andreas Lintermaann
 * \date 16.09.2015
 *
 * \param[in] segmentId the id of the segment to use for the calculation
 *
 **/
template <MInt nDim>
MFloat LbSolver<nDim>::calculateReferenceLength(MInt segmentId) {
  TRACE();

  m_log << "  + Calculating characteristic length:" << endl;
  m_log << "    - type:       " << (m_geometry->m_parallelGeometry ? "parallel" : "serial") << endl;
  m_log << "    - segment id: " << segmentId << endl;

  if constexpr(nDim == 2) {
    stringstream errorMsg;
    errorMsg << "ERROR: no implementation for the calulation of the characteristic length in 2D!" << endl;
    m_log << errorMsg.str();
    TERMM(1, errorMsg.str());
  }

  MFloat ret = 0.0;
  if(m_geometry->m_parallelGeometry)
    ret = calcCharLenParallelSplit(segmentId);
  else
    ret = calcCharLenAll(segmentId);

  m_log << "      * charcteristic length: " << ret << " [cell units]" << endl;

  return ret;
}

template <MInt nDim>
void LbSolver<nDim>::returnCellInfo(MInt cellId) {
  TRACE();

  static ofstream ofl;
  MInt nghbrId = 0;

  // write only the first timesteps
  if(globalTimeStep < FPOW2(maxLevel())) {
    ofl.open("cellInfo.dat", ios_base::app);

    ofl << " cell informations: " << endl;
    ofl << " -------------------" << endl;
    ofl << " id=" << cellId << ", domain=" << domainId() << ", level=" << a_level(cellId) << endl;
    ofl << " (x,y,z)=(" << a_coordinate(cellId, 0) << ", " << a_coordinate(cellId, 1) << ", " << a_coordinate(cellId, 2)
        << ")" << endl;
    ofl << "   neighbors: " << endl;
    for(MInt k = 0; k < 26; k++) {
      nghbrId = c_neighborId(cellId, k);
      ofl << "     " << k << ": noNghbrs=" << a_hasNeighbor(cellId, k) << ", id=" << nghbrId << endl;
      ofl << " (x,y,z)=(" << a_coordinate(nghbrId, 0) << ", " << a_coordinate(nghbrId, 1) << ", "
          << a_coordinate(nghbrId, 2) << ")" << endl;
    }
    ofl << "   children: " << endl;
    ofl << "   noChildren: " << c_noChildren(cellId) << endl;
    for(MInt k = 0; k < IPOW2(nDim); k++) {
      ofl << "     " << k << ": " << c_childId(cellId, k) << endl;
    }
    ofl << "   parent: " << c_parentId(cellId) << endl;

    if(a_isInterface(cellId)) {
      for(MInt k = 0; k < 26; k++) {
        nghbrId = c_neighborId(c_parentId(cellId), k);
        ofl << "     " << k << ": parent noNghbrs=" << a_hasNeighbor(c_parentId(cellId), k) << ", id=" << nghbrId
            << endl;
        ofl << " (x,y,z)=(" << a_coordinate(nghbrId, 0) << ", " << a_coordinate(nghbrId, 1) << ", "
            << a_coordinate(nghbrId, 2) << ")" << endl;
      }
    }

    ofl.close();
  }
}

template <MInt nDim>
void LbSolver<nDim>::writeGridToTecFile(const MChar* fileName) {
  TRACE();
  // open a file
  MInt tmpCounter = 0;
  ofstream ofl;
  ofl.open(fileName);


  if(ofl) {
    if constexpr(nDim == 3) {
      // write header
      ofl << " TITLE = \" AMAZONAS Tecfile \" " << endl;
      ofl << R"( VARIABLES = "X", "Y", "Z", "U", "V", "W", "RHO")" << endl;
      ofl << " ZONE N=" << m_gridPoints->size() << ", ";
      ofl << " E=" << m_extractedCells->size() << ", ";
      ofl << " ZONETYPE="
          << "FEBRICK"
          << ", ";
      ofl << " DATAPACKING="
          << "SOLVER" << endl;
      ofl << " VARLOCATION=([4, 5, 6, 7]=CELLCENTERED)" << endl;
      // write the points to the datafile
    } else {
      // write header
      ofl << " TITLE = \" AMAZONAS Tecfile \" " << endl;
      ofl << R"( VARIABLES = "X", "Y", "U", "V", "RHO")" << endl;
      ofl << " ZONE N=" << m_gridPoints->size() << ", ";
      ofl << " E=" << m_extractedCells->size() << ", ";
      ofl << " ZONETYPE="
          << "FEQUADRILATERAL"
          << ", ";
      ofl << " DATAPACKING="
          << "SOLVER" << endl;
      ofl << " VARLOCATION=([3, 4, 5]=CELLCENTERED)" << endl;
      // write the points to the datafile
    }

    for(MInt axisId = 0; axisId < nDim; axisId++) {
      for(MInt pointId = 0; pointId < m_gridPoints->size(); pointId++) {
        ofl << m_gridPoints->a[pointId].m_coordinates[axisId] << " ";
        tmpCounter++;
        if(tmpCounter == 10) {
          ofl << endl; // don't write everything in one line
          tmpCounter = 0;
        }
      }
    }
    ofl << endl;

    tmpCounter = 0;
    for(MInt varId = 0; varId < nDim + 1; varId++) {
      for(MInt cellId = 0; cellId < m_extractedCells->size(); cellId++) {
        ofl << a_variable(m_extractedCells->a[cellId].m_cellId, varId) << " ";

        tmpCounter++;
        if(tmpCounter == 8) {
          ofl << endl; // don't write everything in one line
          tmpCounter = 0;
        }
      }
      ofl << endl;
    }
    ofl << endl;
    ofl << endl;

    tmpCounter = 0;
    // Write the cells to the datafile, i.e. write all the points of which
    // the cells consist. (TECPLOT needs FORTRAN COUNTING i.e. starting from 1, not 0
    // and TECPLOT needs a special point order (see manual))
    for(MInt cellId = 0; cellId < m_extractedCells->size(); cellId++) {
      ofl << m_extractedCells->a[cellId].m_pointIds[2] + 1 << " ";
      ofl << m_extractedCells->a[cellId].m_pointIds[3] + 1 << " ";
      ofl << m_extractedCells->a[cellId].m_pointIds[1] + 1 << " ";
      ofl << m_extractedCells->a[cellId].m_pointIds[0] + 1 << " ";
      if constexpr(nDim == 3) {
        ofl << m_extractedCells->a[cellId].m_pointIds[6] + 1 << " ";
        ofl << m_extractedCells->a[cellId].m_pointIds[7] + 1 << " ";
        ofl << m_extractedCells->a[cellId].m_pointIds[5] + 1 << " ";
        ofl << m_extractedCells->a[cellId].m_pointIds[4] + 1 << " ";
      }
      ofl << endl;
    }
  } // end if the file could be opened
  else {
    cerr << "ERROR! COULD NOT OPEN FILE " << fileName << " for writing! " << endl;
  }
}

template <MInt nDim>
void LbSolver<nDim>::writeGridToVtkFile(const MChar* fileName) {
  TRACE();
  // open a file

  ofstream ofl;
  MInt noFieldElements = nDim + 1;

#define WRITE_CELLIDS
#define WRITE_OLD_VARIABLES
#define WRITE_DISTRIBUTIONS
#ifdef WRITE_CELLIDS
  noFieldElements++;
#endif
#ifdef WRITE_OLD_VARIABLES
  noFieldElements += nDim + 1;
#endif
#ifdef WRITE_DISTRIBUTIONS
  noFieldElements += m_noDistributions;
#endif

  ofl.open(fileName);
  ofl.precision(12);

  if(ofl) {
    // write header
    ofl << "# vtk DataFile Version 3.0 " << endl;
    ofl << "Amazonas output file" << endl;
    ofl << "BINARY " << endl;
    ofl << "DATASET UNSTRUCTURED_GRID " << endl;
    ofl << "POINTS " << m_gridPoints->size() << " float" << endl;

    // 	ofl.close();
    // 	ofl.open(fileName, ios_base::out|ios_base::app|ios_base::binary);
    float tmpCoordinate = 0;
    for(MInt pointId = 0; pointId < m_gridPoints->size(); pointId++) {
      for(MInt axisId = 0; axisId < nDim; axisId++) {
        tmpCoordinate = (float)m_gridPoints->a[pointId].m_coordinates[axisId];
#ifdef SWAP_ENDIAN
        tmpCoordinate = floatSwap(tmpCoordinate);
#endif
        ofl.write(reinterpret_cast<char*>(&tmpCoordinate), sizeof(float));
        // 	    ofl.write(reinterpret_cast<char *> (m_gridPoints->a[pointId].m_coordinates), sizeof(MFloat)*nDim);
      }

      // Add dummy dimension for 2D
      if constexpr(nDim == 2) {
#ifdef SWAP_ENDIAN
        tmpCoordinate = floatSwap(0.0);
#endif
        ofl.write(reinterpret_cast<char*>(&tmpCoordinate), sizeof(float));
      }
    }
    ofl.close();
    ofl.open(fileName, ios_base::app);
    ofl << endl;
    // ------------------------------------------
    MInt number = IPOW2(nDim);
    ofl << "CELLS " << m_extractedCells->size() << " " << m_extractedCells->size() * (number + 1) << endl;
    // Write the cells to the datafile, i.e. write all the points of which
    // the cells consist.
    ofl.close();
    ofl.open(fileName, ios_base::out | ios_base::app | ios_base::binary);

    MInt pointId = 0;
#ifdef SWAP_ENDIAN
    number = intSwap(number);
#endif
    for(MInt cellId = 0; cellId < m_extractedCells->size(); cellId++) {
      ofl.write(reinterpret_cast<char*>(&number), sizeof(MInt));
      //  	  ofl.write(reinterpret_cast<char *> (m_extractedCells->a[ cellId ].m_pointIds), sizeof(MInt)*number);
      for(MInt points = 0; points < IPOW2(nDim); points++) {
        pointId = m_extractedCells->a[cellId].m_pointIds[points];
#ifdef SWAP_ENDIAN
        pointId = intSwap(pointId);
#endif
        ofl.write(reinterpret_cast<char*>(&pointId), sizeof(MInt));
      }
    }
    ofl.close();
    ofl.open(fileName, ios_base::app);
    ofl << endl;
    ofl << "CELL_TYPES " << m_extractedCells->size() << endl;
    ofl.close();
    ofl.open(fileName, ios_base::out | ios_base::app | ios_base::binary);
    if constexpr(nDim == 3) {
      pointId = 11; // use VTK_VOXEL in 3D
    } else {
      pointId = 8; // use VTK_PIXEL in 2D
    }

#ifdef SWAP_ENDIAN
    pointId = intSwap(pointId);
#endif
    for(MInt cellId = 0; cellId < m_extractedCells->size(); cellId++) {
      ofl.write(reinterpret_cast<char*>(&pointId), sizeof(MInt));
    }

    ofl.close();
    ofl.open(fileName, ios_base::app);
    ofl << endl;
    ofl << "CELL_DATA " << m_extractedCells->size() << endl;
    ofl << "SCALARS scalars double" << endl;
    ofl << "LOOKUP_TABLE default" << endl;

    ofl.close();
    ofl.open(fileName, ios_base::out | ios_base::app | ios_base::binary);
    MFloat fnumber = 1.0;
#ifdef SWAP_ENDIAN
    fnumber = doubleSwap(fnumber);
#endif
    for(MInt cellId = 0; cellId < m_extractedCells->size(); cellId++) {
      ofl.write(reinterpret_cast<char*>(&fnumber), sizeof(MFloat));
    }
    ofl.close();

#ifdef WRITE_VECTOR_DATA
    // Write velocities as vector
    ofl.open(fileName, ios_base::app);
    ofl << endl;
    ofl << "VECTORS vectors double" << endl;
    ofl.close();

    ofl.open(fileName, ios_base::out | ios_base::app | ios_base::binary);
    for(MInt cellId = 0; cellId < m_extractedCells->size(); cellId++) {
      for(MInt varId = 0; varId < nDim; varId++) {
        fnumber = a_variable(m_extractedCells->a[cellId].m_cellId, varId);
#ifdef SWAP_ENDIAN
        fnumber = doubleSwap(fnumber);
#endif
        ofl.write(reinterpret_cast<char*>(&fnumber), sizeof(MFloat));
      }
    }

    ofl.close();
#endif

    // Write all variables as field data
    ofl.open(fileName, ios_base::app);
    ofl << endl;
    ofl << "FIELD FieldData " << noFieldElements << endl;

    for(MInt varId = 0; varId < nDim + 1; varId++) {
      ofl.close();
      ofl.open(fileName, ios_base::app);
      ofl << endl;
      ofl << "variables" << varId << " 1 " << m_extractedCells->size() << " double" << endl;
      ofl.close();
      ofl.open(fileName, ios_base::out | ios_base::app | ios_base::binary);
      for(MInt cellId = 0; cellId < m_extractedCells->size(); cellId++) {
        // 	    if(varId==nDim)//subtract 1 from density
        // 	      fnumber = a_variable(m_extractedCells->a[cellId].m_cellId, varId);
        // 	    else
        //  	    if(varId==nDim)//write actual domainId
        //  	      fnumber = domainId();
        //  	    else
        fnumber = a_variable(m_extractedCells->a[cellId].m_cellId, varId);
#ifdef SWAP_ENDIAN
        fnumber = doubleSwap(fnumber);
#endif
        ofl.write(reinterpret_cast<char*>(&fnumber), sizeof(MFloat));
      }
    }

#ifdef WRITE_OLD_VARIABLES
    for(MInt varId = 0; varId < nDim + 1; varId++) {
      ofl.close();
      ofl.open(fileName, ios_base::app);
      ofl << endl;
      ofl << "residual" << varId << " 1 " << m_extractedCells->size() << " double" << endl;
      ofl.close();
      ofl.open(fileName, ios_base::out | ios_base::app | ios_base::binary);
      for(MInt cellId = 0; cellId < m_extractedCells->size(); cellId++) {
        fnumber = a_variable(m_extractedCells->a[cellId].m_cellId, varId)
                  - a_oldVariable(m_extractedCells->a[cellId].m_cellId, varId);
#ifdef SWAP_ENDIAN
        fnumber = doubleSwap(fnumber);
#endif
        ofl.write(reinterpret_cast<char*>(&fnumber), sizeof(MFloat));
      }
    }
#endif
#ifdef WRITE_DISTRIBUTIONS
    /***********WRITE OUT CELL DISTRIBUTIONS***********************/
    for(MInt distId = 0; distId < m_noDistributions; distId++) {
      ofl.close();
      ofl.open(fileName, ios_base::app);
      ofl << endl;
      ofl << "distributions" << distId << " 1 " << m_extractedCells->size() << " double" << endl;
      ofl.close();
      ofl.open(fileName, ios_base::out | ios_base::app | ios_base::binary);
      for(MInt cellId = 0; cellId < m_extractedCells->size(); cellId++) {
        fnumber = a_distribution(m_extractedCells->a[cellId].m_cellId, distId);
#ifdef SWAP_ENDIAN
        fnumber = doubleSwap(fnumber);
#endif
        ofl.write(reinterpret_cast<char*>(&fnumber), sizeof(MFloat));
      }
    }
#endif

#ifdef WRITE_CELLIDS
    /***********WRITE OUT CELL IDS***********************/
    MInt inumber = 0;
    ofl.close();
    ofl.open(fileName, ios_base::app);
    ofl << endl;
    ofl << "cellIds "
        << " 1 " << m_extractedCells->size() << " int" << endl;
    ofl.close();
    ofl.open(fileName, ios_base::out | ios_base::app | ios_base::binary);
    for(MInt cellId = 0; cellId < m_extractedCells->size(); cellId++) {
      inumber = m_extractedCells->a[cellId].m_cellId;
#ifdef SWAP_ENDIAN
      inumber = intSwap(inumber);
#endif
      ofl.write(reinterpret_cast<char*>(&inumber), sizeof(MInt));
    }
#endif
    ofl.close();
  } // end if the file could be opened
  else {
    cerr << "ERROR! COULD NOT OPEN FILE " << fileName << " for writing! " << endl;
  }
}

template <MInt nDim>
void LbSolver<nDim>::writeGridToVtkFileAscii(const MChar* fileName) {
  TRACE();
  // open a file
  MInt tmpCounter = 0;
  ofstream ofl;
  ofl.open(fileName);

  if(ofl) {
    if constexpr(nDim == 3) {
      // write header
      ofl << "# vtk DataFile Version 2.0 " << endl;
      ofl << "Amazonas output file" << endl;
      ofl << "ASCII " << endl;
      ofl << "DATASET UNSTRUCTURED_GRID " << endl;
      ofl << "POINTS " << m_gridPoints->size() << " float" << endl;

      for(MInt pointId = 0; pointId < m_gridPoints->size(); pointId++) {
        for(MInt axisId = 0; axisId < nDim; axisId++) {
          ofl << m_gridPoints->a[pointId].m_coordinates[axisId] << " ";
        }
        ofl << "   ";
        tmpCounter++;
        if(tmpCounter == 10) {
          ofl << endl; // don't write everything in one line
          tmpCounter = 0;
        }
      }
      ofl << endl;
      ofl << endl;

      ofl << "CELLS " << m_extractedCells->size() << " " << m_extractedCells->size() * 9 << endl;
      // Write the cells to the datafile, i.e. write all the points of which
      // the cells consist.
      for(MInt cellId = 0; cellId < m_extractedCells->size(); cellId++) {
        ofl << "8 ";
        ofl << m_extractedCells->a[cellId].m_pointIds[0] << " ";
        ofl << m_extractedCells->a[cellId].m_pointIds[1] << " ";
        ofl << m_extractedCells->a[cellId].m_pointIds[2] << " ";
        ofl << m_extractedCells->a[cellId].m_pointIds[3] << " ";
        if constexpr(nDim == 3) {
          ofl << m_extractedCells->a[cellId].m_pointIds[4] << " ";
          ofl << m_extractedCells->a[cellId].m_pointIds[5] << " ";
          ofl << m_extractedCells->a[cellId].m_pointIds[6] << " ";
          ofl << m_extractedCells->a[cellId].m_pointIds[7] << " ";
        }
        ofl << endl;
      }

      ofl << "CELL_TYPES " << m_extractedCells->size() << endl;
      for(MInt cellId = 0; cellId < m_extractedCells->size(); cellId++) {
        ofl << "11 " << endl; // VTK_VOXEL
      }
      ofl << "CELL_DATA " << m_extractedCells->size() << endl;
      ofl << "SCALARS scalars float" << endl;
      ofl << "LOOKUP_TABLE default" << endl;
      for(MInt cellId = 0; cellId < m_extractedCells->size(); cellId++) {
        ofl << "1.0" << endl; // VTK_VOXEL
      }
      ofl << "FIELD FieldData " << nDim + 1 << endl;
      // Write variables
      for(MInt varId = 0; varId < nDim + 1; varId++) {
        ofl << "variables" << varId << " 1 " << m_extractedCells->size() << " float" << endl;
        tmpCounter = 0;
        for(MInt cellId = 0; cellId < m_extractedCells->size(); cellId++) {
          ofl << a_variable(m_extractedCells->a[cellId].m_cellId, varId) << " ";
          tmpCounter++;
        }
        ofl << endl; // don't write everything in one line
      }
      ofl << endl;

      // Write variables old variables
      ofl << "FIELD FieldData " << nDim + 1 << endl;
      for(MInt varId = 0; varId < nDim + 1; varId++) {
        ofl << "oldVariables" << varId << " 1 " << m_extractedCells->size() << " float" << endl;
        tmpCounter = 0;
        for(MInt cellId = 0; cellId < m_extractedCells->size(); cellId++) {
          ofl << a_oldVariable(m_extractedCells->a[cellId].m_cellId, varId) << " ";
          tmpCounter++;
        }
        ofl << endl; // don't write everything in one line
      }
      ofl << endl;
    }
  } // end if the file could be opened
  else {
    cerr << "ERROR! COULD NOT OPEN FILE " << fileName << " for writing! " << endl;
  }
}

/** \brief loads a restart from a coarser solution and interpolates to a finer mesh
 *
 * \author Andreas Lintermann
 * \date 29.01.2016
 *
 * The algorithm does the following:
 *
 *  1.  count all non-leaf cells and determine offsets for reading the file
 *  2.  determine ids of the coarser cells
 *  3.  read variables
 *  4.  exchange coarse windows
 *  4.1 collect windows and halos on a coarse level
 *  4.2 construct send buffer and send
 *  4.3 construct receive buffer and receive
 *  4.4 distribute received information
 *  5.  fill finer cells
 *  5.1 fill all cells surrounded by neighbors by interpolation all others by the restart values
 *
 * \param[in] fileName the name of the file to read from
 *
 **/
template <MInt nDim>
void LbSolver<nDim>::loadRestartWithoutDistributionsParFromCoarse(const MChar* fileName) {
  TRACE();

  if constexpr(nDim != 2 && nDim != 3) {
    cerr << " In global function loadGridFlowVariablesThermalPar: wrong number of dimensions !" << endl;
    exit(0);
    return;
  }

  using namespace maia::parallel_io;
  ParallelIo parallelIo(fileName, PIO_READ, mpiComm());

  // read time
  if(!m_initFromRestartFile) {
    parallelIo.getAttribute(&m_time, "time");
    parallelIo.getAttribute(&globalTimeStep, "globalTimeStep");
  }

  if(m_velocityControl.restart) {
    parallelIo.getAttribute(&m_velocityControl.lastGlobalAvgV, "lastGlobalAvgV");
    parallelIo.getAttribute(&m_velocityControl.previousError, "previousError");
    parallelIo.getAttribute(&m_velocityControl.integratedError, "integratedError");
    parallelIo.getAttribute(&m_velocityControl.derivedError, "derivedError");
    for(MInt i = 0; i < nDim; i++) {
      parallelIo.getAttribute(&m_volumeAccel[i], "volumeAcceleration_" + std::to_string(i));
    }
  }

  MInt coarse_count = 0;
  parallelIo.readScalar(&coarse_count, "noCells");

  // 1. count all non-leaf cells and determine offsets for reading the file
  MInt local_count = 0;
  for(MInt i = 0; i < noInternalCells(); i++) {
    if(c_noChildren(i) > 0) local_count++;
  }

  MIntScratchSpace nonFine_count(noDomains(), AT_, "nonFineCount");
  for(MInt d = 0; d < noDomains(); d++)
    if(d == domainId())
      nonFine_count[d] = local_count;
    else
      nonFine_count[d] = 0;

  MPI_Allreduce(MPI_IN_PLACE, nonFine_count.getPointer(), noDomains(), MPI_INT, MPI_SUM, mpiComm(), AT_, "MPI_IN_PLACE",
                "nonFine_count.getPointer()");

  MInt off = 0;
  for(MInt d = 0; d < domainId(); d++)
    off += nonFine_count[d];

  MInt all = 0;
  for(MInt d = 0; d < noDomains(); d++)
    all += nonFine_count[d];

  m_log << "    - loading information: " << endl;
  m_log << "      * number of all cells:   " << all << endl;
  m_log << "      * number of local cells: " << local_count << endl;
  m_log << "      * domain offset:         " << off << endl;

  const MPI_Offset firstGlobalId = off;
  const MPI_Offset localNoCells = local_count;
  parallelIo.setOffset(localNoCells, firstGlobalId);

  m_log << "    - determining ids of the coarse cells" << endl;
  // 2. determine ids of the coarser cells
  MIntScratchSpace coarse_ids(local_count, AT_, "coarse_ids");
  for(MInt i = 0, j = 0; i < noInternalCells(); ++i)
    if(c_noChildren(i) > 0) {
      coarse_ids[j] = i;
      j++;
    }

  m_log << "    - reading variables " << endl;
  // 3. read variables
  {
    MFloatScratchSpace tmparray(local_count, AT_, "tmparray");

    // first read all variables
    for(MInt v = 0; v < m_noVariables; v++) {
      MString name = "variables" + to_string(v);

      parallelIo.readArray(tmparray.getPointer(), name);
      for(MInt i = 0; i < local_count; ++i)
        a_variable(coarse_ids[i], v) = tmparray.p[i];
    }

    // now read all oldVariables
    for(MInt v = 0; v < m_noVariables; v++) {
      MString name = "oldVariables" + to_string(v);

      parallelIo.readArray(tmparray.getPointer(), name);
      for(MInt i = 0; i < local_count; ++i)
        a_variable(coarse_ids[i], v) = tmparray.p[i];
    }
  }

  m_log << "    - exchanging coarse window cells" << endl;
  // 4. exchange coarse windows
  // 4.1 collect windows and halos on a coarse level
  m_log << "      * collect windows and halos on a coarse level " << endl;
  vector<vector<MInt>> coarseWindowPD;
  vector<vector<MInt>> coarseHaloPD;
  for(MInt d = 0; d < noNeighborDomains(); d++) {
    vector<MInt> cW;
    for(MInt w = 0; w < noWindowCells(d); w++)
      if(a_level(windowCell(d, w)) == maxLevel() - 1) cW.push_back(windowCell(d, w));

    vector<MInt> cH;
    for(MInt h = 0; h < noHaloCells(d); h++)
      if(a_level(haloCell(d, h)) == maxLevel() - 1) cH.push_back(haloCell(d, h));

    coarseWindowPD.push_back(cW);
    coarseHaloPD.push_back(cH);
  }

  m_log << "        = number of coarse window cells for " << noNeighborDomains() << " neighboring domains" << endl;
  for(MInt d = 0; d < noNeighborDomains(); d++)
    m_log << "          # domain " << neighborDomain(d) << ":" << coarseWindowPD[d].size() << endl;
  m_log << "        = number of coarse halo cells for " << noNeighborDomains() << " neighboring domains" << endl;
  for(MInt d = 0; d < noNeighborDomains(); d++)
    m_log << "          # domain " << neighborDomain(d) << ":" << coarseHaloPD[d].size() << endl;

  // 4.2 construct send buffer and send
  m_log << "      * constructing send buffer with size: ";
  MInt totalWindowBuf = 0;
  MIntScratchSpace windowOff(noNeighborDomains(), AT_, "windowOff");
  for(MInt d = 0; d < noNeighborDomains(); d++) {
    windowOff[d] = totalWindowBuf;
    totalWindowBuf += coarseWindowPD[d].size() * 2 * m_noVariables;
  }

  m_log << totalWindowBuf << endl;

  MFloatScratchSpace sendBuf(totalWindowBuf, AT_, "sendBuf");
  for(MInt d = 0; d < noNeighborDomains(); d++) {
    MInt pos = windowOff[d];
    for(MInt c = 0; c < (MInt)coarseWindowPD[d].size(); c++) {
      for(MInt v = 0; v < m_noVariables; v++, pos++)
        sendBuf[pos] = a_variable(coarseWindowPD[d][c], v);
      for(MInt v = 0; v < m_noVariables; v++, pos++)
        sendBuf[pos] = a_oldVariable(coarseWindowPD[d][c], v);
    }
    MPI_Issend(&sendBuf[windowOff[d]], coarseWindowPD[d].size() * 2 * m_noVariables, MPI_DOUBLE, neighborDomain(d), 0,
               mpiComm(), &mpi_request[d], AT_, "sendBuf[windowOff[d]]");
  }

  // 4.3 construct receive buffer and receive
  m_log << "      * constructing receive buffer and receive" << endl;
  MInt totalHaloBuf = 0;
  MIntScratchSpace haloOff(noNeighborDomains(), AT_, "haloOff");
  for(MInt d = 0; d < noNeighborDomains(); d++) {
    haloOff[d] = totalHaloBuf;
    totalHaloBuf += coarseHaloPD[d].size() * 2 * m_noVariables;
  }

  MFloatScratchSpace receiveBuf(totalHaloBuf, AT_, "receiveBuf");

  MPI_Status status;
  for(MInt d = 0; d < noNeighborDomains(); d++)
    MPI_Recv(&receiveBuf[haloOff[d]], coarseHaloPD[d].size() * 2 * m_noVariables, MPI_DOUBLE, neighborDomain(d), 0,
             mpiComm(), &status, AT_, "receiveBuf[haloOff[d]]");

  for(MInt d = 0; d < noNeighborDomains(); d++)
    MPI_Wait(&mpi_request[d], &status, AT_);

  // 4.4 distribute received information
  m_log << "      * distributing received information" << endl;
  for(MInt d = 0; d < noNeighborDomains(); d++) {
    MInt pos = haloOff[d];
    for(MInt c = 0; c < (MInt)coarseHaloPD[d].size(); c++) {
      for(MInt v = 0; v < m_noVariables; v++, pos++)
        a_variable(coarseHaloPD[d][c], v) = receiveBuf[pos];
      for(MInt v = 0; v < m_noVariables; v++, pos++)
        a_oldVariable(coarseHaloPD[d][c], v) = receiveBuf[pos];
    }
  }

  // 5. fill finer cells
  m_log << "    - filling all finer cells" << endl;
  MInt elems = IPOW2(nDim);
  MIntScratchSpace dirs(elems, AT_, "elems");

  if constexpr(nDim == 2) {
    for(MInt d = 0; d < elems; d++) {
      dirs[d] = childPos2D[d];
    }
  } else {
    for(MInt d = 0; d < elems; d++) {
      dirs[d] = childPos3D[d];
    }
  }

  // 5.1 fill all cells surrounded by neighbors by interpolation all others by the restart values
  m_log << "      * filling all cells with nieghbors" << endl;
  for(MInt i = 0; i < noInternalCells(); i++)
    if(a_level(i) == maxLevel()) {
      MInt parent = c_parentId(i);

      // find myself under my parent
      MInt l = 0;
      for(MInt c = 0; c < elems; c++)
        if(c_childId(parent, c) == i) {
          l = c;
          break;
        }

      MInt dir = dirs[l];

      // great neighbor in direction
      if(a_hasNeighbor(parent, dir)) {
        MInt parentNeighbor = c_neighborId(parent, dir);
        for(MInt v = 0; v < m_noVariables; v++) {
          a_variable(i, v) = F3B4 * a_variable(parent, v) + F1B4 * a_variable(parentNeighbor, v);
          a_oldVariable(i, v) = F3B4 * a_oldVariable(parent, v) + F1B4 * a_oldVariable(parentNeighbor, v);
        }
      }
      // otherwise fill with value from restart
      else {
        for(MInt v = 0; v < m_noVariables; v++) {
          a_variable(i, v) = a_variable(parent, v);
          a_oldVariable(i, v) = a_oldVariable(parent, v);
        }
      }
    }

  // 5.2 fill finer halo cells
  m_log << "      * filling halos" << endl;

  MInt noNeigh = pow(nDim, nDim);
  for(MInt d = 0; d < noNeighborDomains(); d++)
    for(MInt h = 0; h < noHaloCells(d); h++) {
      MInt haloId = haloCell(d, h);

      if(a_level(haloId) == maxLevel()) {
        // find a neighbor to the halo cell which is a window cell
        MInt haloWinNeigh = -1;
        MInt dir = 0;
        for(; dir < noNeigh; dir++)
          if(a_hasNeighbor(haloId, dir) && a_isWindow(c_neighborId(haloId, dir))) {
            haloWinNeigh = c_neighborId(haloId, dir);
            break;
          }

        // found a window cell
        if(haloWinNeigh >= 0) {
          MInt winParent = c_parentId(haloWinNeigh);
          const MInt opp =
              LbLatticeDescriptorBase<3>::oppositeDist(dir); // TODO labels:LB dxqy: 3 replaceable by nDim ?

          // check neighbor of the parent of the window cell in opoosite direction
          MInt parent = -1;
          if(winParent >= 0 && a_hasNeighbor(winParent, opp)) parent = c_neighborId(winParent, opp);

          if(parent < 0) parent = winParent;

          for(MInt v = 0; v < m_noVariables; v++) {
            a_variable(haloId, v) = a_variable(parent, v);
            a_oldVariable(haloId, v) = a_oldVariable(parent, v);
          }
        }
      }
    }
}

/** \brief This function loads the flow information of the cells
 *        such as variables and attributes like u_velocity,density,etc. for restart purpose.
 *
 * \author Andreas Lintermann, Konstantin Froehlich
 * \date 13.02.2012
 * In contrast to loadGridFlowVariables(const MChar* fileName) this function loads the
 * information for restart in parallel.
 * The attribute 'name' of the variables are set to their according meaning.
 *
 * \param[in] fileName the name of the file to load the data from
 *
 **/
template <MInt nDim>
void LbSolver<nDim>::loadRestartWithDistributionsPar(const MChar* fileName) {
  TRACE();
  if constexpr(nDim != 2 && nDim != 3) {
    cerr << " In global function loadGridFlowVariablesThermalPar: wrong number of dimensions !" << endl;
    exit(0);
    return;
  }

  using namespace maia::parallel_io;
  ParallelIo parallelIo(fileName, PIO_READ, mpiComm());

  // read time
  if(!m_initFromRestartFile) {
    parallelIo.getAttribute(&m_time, "time");
    parallelIo.getAttribute(&globalTimeStep, "globalTimeStep");
  }

  if(m_velocityControl.restart) {
    parallelIo.getAttribute(&m_velocityControl.lastGlobalAvgV, "lastGlobalAvgV");
    parallelIo.getAttribute(&m_velocityControl.previousError, "previousError");
    parallelIo.getAttribute(&m_velocityControl.integratedError, "integratedError");
    parallelIo.getAttribute(&m_velocityControl.derivedError, "derivedError");
    for(MInt i = 0; i < nDim; i++) {
      parallelIo.getAttribute(&m_volumeAccel[i], "volumeAcceleration_" + std::to_string(i));
    }
  }

  MInt pv_veloc[nDim];
  if constexpr(nDim == 2) {
    pv_veloc[0] = PV->U;
    pv_veloc[1] = PV->V;
  } else {
    pv_veloc[0] = PV->U;
    pv_veloc[1] = PV->V;
    pv_veloc[2] = PV->W;
  }
  MInt pvrho = PV->RHO;
  MInt pvt = PV->T;
  MInt pvc = PV->C;
  MString rhoName;
  MString oldrhoName;
  MString temperatureName;
  MString oldtemperatureName;
  MString concentrationName;
  MString oldconcentrationName;
  MString EELiquidName;
  MString distributions;
  MString oldDistributions;
  MString distributionsThermal;
  MString oldDistributionsThermal;
  MString distributionsTransport;
  MString oldDistributionsTransport;
  if constexpr(nDim == 3) {
    rhoName = "variables3";
    oldrhoName = "oldVariables3";
    temperatureName = "variables4";
    oldtemperatureName = "oldVariables4";
    concentrationName = "variables5";
    oldconcentrationName = "oldVariables5";
    EELiquidName = "oldVariables4";
  } else {
    rhoName = "variables2";
    oldrhoName = "oldVariables2";
    temperatureName = "variables3";
    oldtemperatureName = "oldVariables3";
    concentrationName = "variables4";
    oldconcentrationName = "oldVariables4";
  }


  // Set file offsets (first globalId and #of cells to be read by this process):
  const MPI_Offset firstGlobalId = domainOffset(domainId());
  const MPI_Offset localNoCells = noInternalCells();
  parallelIo.setOffset(localNoCells, firstGlobalId);


  // first read 2D variables
  {
    MFloatScratchSpace tmparray(noInternalCells(), AT_, "tmparray");

    // Velocity u
    parallelIo.readArray(tmparray.getPointer(), "variables0");
    for(MInt i = 0; i < noInternalCells(); ++i)
      a_variable(i, pv_veloc[0]) = tmparray[i];

    // Velocity old_u
    parallelIo.readArray(tmparray.getPointer(), "oldVariables0");
    for(MInt i = 0; i < noInternalCells(); ++i)
      a_oldVariable(i, pv_veloc[0]) = tmparray[i];

    // Velocity v
    parallelIo.readArray(tmparray.getPointer(), "variables1");
    for(MInt i = 0; i < noInternalCells(); ++i)
      a_variable(i, pv_veloc[1]) = tmparray.p[i];

    // Velocity old_v
    parallelIo.readArray(tmparray.getPointer(), "oldVariables1");
    for(MInt i = 0; i < noInternalCells(); ++i)
      a_oldVariable(i, pv_veloc[1]) = tmparray.p[i];

    // Density rho
    parallelIo.readArray(tmparray.getPointer(), rhoName);
    for(MInt i = 0; i < noInternalCells(); ++i)
      a_variable(i, pvrho) = tmparray.p[i];

    // Density old_rho
    parallelIo.readArray(tmparray.getPointer(), oldrhoName);
    for(MInt i = 0; i < noInternalCells(); ++i)
      a_oldVariable(i, pvrho) = tmparray.p[i];

    if(m_isThermal) {
      // Temperature t
      parallelIo.readArray(tmparray.getPointer(), temperatureName);
      for(MInt i = 0; i < noInternalCells(); ++i)
        a_variable(i, pvt) = tmparray.p[i];

      // Temperature old_t
      parallelIo.readArray(tmparray.getPointer(), oldtemperatureName);
      for(MInt i = 0; i < noInternalCells(); ++i)
        a_oldVariable(i, pvt) = tmparray.p[i];
    }
    if(m_isTransport) {
      // Concentration c
      parallelIo.readArray(tmparray.getPointer(), concentrationName);
      for(MInt i = 0; i < noInternalCells(); ++i)
        a_variable(i, pvc) = tmparray.p[i];

      // Concentration old_c
      parallelIo.readArray(tmparray.getPointer(), oldconcentrationName);
      for(MInt i = 0; i < noInternalCells(); ++i)
        a_oldVariable(i, pvc) = tmparray.p[i];
    }
    if(m_isEELiquid) {
      if(!m_EELiquid.restartWithoutAlpha) {
        parallelIo.readArray(tmparray.getPointer(), EELiquidName);
        for(MInt i = 0; i < noInternalCells(); ++i)
          a_alphaGas(i) = tmparray.p[i];
      } else {
        for(MInt i = 0; i < noInternalCells(); ++i)
          a_alphaGas(i) = 0.0;
      }
    }
  }

  // now read 3D variables
  if constexpr(nDim == 3) {
    MFloatScratchSpace tmparray(noInternalCells(), AT_, "tmparray");

    // Velocity w
    parallelIo.readArray(tmparray.getPointer(), "variables2");
    for(MInt i = 0; i < noInternalCells(); ++i)
      a_variable(i, pv_veloc[2]) = tmparray.p[i];

    // Velocity old_w
    parallelIo.readArray(tmparray.getPointer(), "oldVariables2");
    for(MInt i = 0; i < noInternalCells(); ++i)
      a_oldVariable(i, pv_veloc[2]) = tmparray.p[i];
  }


  // Distributions
  for(MInt j = 0; j < m_noDistributions; j++) {
    distributions = "distributions" + to_string(j);

    MFloatScratchSpace tmp_distributions(noInternalCells(), AT_, "tmp_distributions");
    parallelIo.readArray(tmp_distributions.getPointer(), distributions);
    for(MInt i = 0; i < noInternalCells(); ++i)
      a_distribution(i, j) = tmp_distributions.p[i];
  }

  // old Distributions
  for(MInt j = 0; j < m_noDistributions; j++) {
    oldDistributions = "oldDistributions" + to_string(j);

    MFloatScratchSpace tmp_oldDistributions(noInternalCells(), AT_, "tmp_oldDistributions");
    parallelIo.readArray(tmp_oldDistributions.getPointer(), oldDistributions);
    for(MInt i = 0; i < noInternalCells(); ++i)
      a_oldDistribution(i, j) = tmp_oldDistributions.p[i];
  }

  if(m_isThermal) {
    // Thermal Distributions
    for(MInt j = 0; j < m_noDistributions; j++) {
      distributionsThermal = "distributionsThermal" + to_string(j);

      MFloatScratchSpace tmp_distributionsThermal(noInternalCells(), AT_, "tmp_distributionsThermal");
      parallelIo.readArray(tmp_distributionsThermal.getPointer(), distributionsThermal);
      for(MInt i = 0; i < noInternalCells(); ++i)
        a_distributionThermal(i, j) = tmp_distributionsThermal.p[i];
    }

    // old Thermal Distributions
    for(MInt j = 0; j < m_noDistributions; j++) {
      oldDistributionsThermal = "oldDistributionsThermal" + to_string(j);

      MFloatScratchSpace tmp_oldDistributionsThermal(noInternalCells(), AT_, "tmp_oldDistributionsThermal");
      parallelIo.readArray(tmp_oldDistributionsThermal.getPointer(), oldDistributionsThermal);
      for(MInt i = 0; i < noInternalCells(); ++i)
        a_oldDistributionThermal(i, j) = tmp_oldDistributionsThermal.p[i];
    }

    for(MInt i = 0; i < a_noCells(); i++) {
      a_kappa(i) = m_kappa;
    }
  }

  if(m_isTransport) {
    // Transport Distributions
    for(MInt j = 0; j < m_noDistributions; j++) {
      distributionsTransport = "distributionsTransport" + to_string(j);

      MFloatScratchSpace tmp_distributionsTransport(noInternalCells(), AT_, "tmp_distributionsTransport");
      parallelIo.readArray(tmp_distributionsTransport.getPointer(), distributionsTransport);
      for(MInt i = 0; i < noInternalCells(); ++i)
        a_distributionTransport(i, j) = tmp_distributionsTransport.p[i];
    }

    // old Transport Distributions
    for(MInt j = 0; j < m_noDistributions; j++) {
      oldDistributionsTransport = "oldDistributionsTransport" + to_string(j);

      MFloatScratchSpace tmp_oldDistributionsTransport(noInternalCells(), AT_, "tmp_oldDistributionsTransport");
      parallelIo.readArray(tmp_oldDistributionsTransport.getPointer(), oldDistributionsTransport);
      for(MInt i = 0; i < noInternalCells(); ++i)
        a_oldDistributionTransport(i, j) = tmp_oldDistributionsTransport.p[i];
    }

    for(MInt i = 0; i < a_noCells(); i++) {
      a_diffusivity(i) = m_diffusivity;
    }
  }

  if(!m_nonNewtonian) {
    for(MInt i = 0; i < a_noCells(); i++) {
      initNu(i, m_nu);
    }
  } else {
    MFloatScratchSpace tmp_nu(noInternalCells(), AT_, "tmp_nu");
    parallelIo.readArray(tmp_nu.getPointer(), "viscosity");
    for(MInt i = 0; i < noInternalCells(); i++) {
      initNu(i, tmp_nu(i));
    }
    this->exchangeData(&a_nu(0));
  }

  m_bndCnd->bcDataReadRestartData(parallelIo);

  // do an exchange to have the right values in the halo cells
  if(noDomains() > 1) {
    exchange(1);
    exchangeOldDistributions();
  } else if(noNeighborDomains() > 0) {
    // TODO labels:LB Periodic domain on noDomain == 1 also creates halos that
    //                demand for an exchange. This has to be treaten somehow
    std::stringstream ss;
    ss << "WARNING: Serial run and periodic (?) is problematic as halo cells are"
          " not exchanged even though periodic BC make usage of it!"
       << std::endl;
    m_log << ss.str();
    std::cerr << ss.str();
  }
}

/** \brief Exchanges the old distributions which is required for a restart
 *
 * \author Andreas Lintermann
 * \date 21.012.2015
 *
 **/
template <MInt nDim>
void LbSolver<nDim>::exchangeOldDistributions() {
  TRACE();

  MInt noTransfer = 0;

  if(m_isThermal && !m_isTransport) {
    noTransfer = 2;
  } else if(m_isTransport && !m_isThermal) {
    noTransfer = 2;
  } else if(m_isTransport && m_isThermal) {
    noTransfer = 3;
  } else {
    noTransfer = 1;
  }

  MInt sumwin = 0;
  MInt sumhalo = 0;
  for(MInt n = 0; n < noNeighborDomains(); n++) {
    sumwin += noWindowCells(n);
    sumhalo += noHaloCells(n);
  }

  MFloatScratchSpace sendBuf(noTransfer * m_noDistributions * sumwin, AT_, "sendBuf");
  MFloatScratchSpace recBuf(noTransfer * m_noDistributions * sumhalo, AT_, "recBuf");

  MIntScratchSpace sendNeighOffset(noNeighborDomains(), AT_, "sendNeighOffset");
  MIntScratchSpace recNeighOffset(noNeighborDomains(), AT_, "recNeighOffset");
  sendNeighOffset[0] = 0;
  recNeighOffset[0] = 0;

  for(MInt n = 1; n < noNeighborDomains(); n++) {
    sendNeighOffset[n] = sendNeighOffset[n - 1] + noWindowCells(n - 1) * noTransfer * m_noDistributions;
    recNeighOffset[n] = recNeighOffset[n - 1] + noHaloCells(n - 1) * noTransfer * m_noDistributions;
  }

  for(MInt n = 0; n < noNeighborDomains(); n++) {
    MInt sndBuf = sendNeighOffset[n];
    for(MInt j = 0; j < noWindowCells(n); j++) {
      for(MInt dist = 0; dist < m_noDistributions; dist++) {
        MInt offset = j * m_noDistributions + dist;
        sendBuf[sndBuf + offset] = a_oldDistribution(windowCell(n, j), dist);
      }
    }
    if(m_isThermal) {
      sndBuf += m_noDistributions * noWindowCells(n);
      for(MInt j = 0; j < noWindowCells(n); j++) {
        for(MInt dist = 0; dist < m_noDistributions; dist++) {
          MInt offset = j * m_noDistributions + dist;
          sendBuf[sndBuf + offset] = a_oldDistributionThermal(windowCell(n, j), dist);
        }
      }
    }
    if(m_isTransport) {
      sndBuf += m_noDistributions * noWindowCells(n);
      for(MInt j = 0; j < noWindowCells(n); j++) {
        for(MInt dist = 0; dist < m_noDistributions; dist++) {
          MInt offset = j * m_noDistributions + dist;
          sendBuf[sndBuf + offset] = a_oldDistributionTransport(windowCell(n, j), dist);
        }
      }
    }
  }

  // 2. Send the gathered information.
  for(MInt n = 0; n < noNeighborDomains(); n++) {
    MInt bufsize = noWindowCells(n) * noTransfer * m_noDistributions;
    MPI_Issend(&sendBuf[sendNeighOffset[n]], bufsize, MPI_DOUBLE, neighborDomain(n), 0, mpiComm(), &mpi_request[n], AT_,
               "sendBuf[sendNeighOffset[n]]");
  }

  // 3. Receive data from neighbors.
  MPI_Status status;

  for(MInt n = 0; n < noNeighborDomains(); n++) {
    MInt bufsize = noHaloCells(n) * noTransfer * m_noDistributions;
    MPI_Recv(&recBuf[recNeighOffset[n]], bufsize, MPI_DOUBLE, neighborDomain(n), 0, mpiComm(), &status, AT_,
             "recBuf[recNeighOffset[n]]");
  }

  for(MInt n = 0; n < noNeighborDomains(); n++)
    MPI_Wait(&mpi_request[n], &status, AT_);

  // 4. scatter the received data
  for(MInt n = 0; n < noNeighborDomains(); n++) {
    MInt recBuffer = recNeighOffset[n];
    for(MInt j = 0; j < noHaloCells(n); j++) {
      for(MInt dist = 0; dist < m_noDistributions; dist++) {
        MInt offset = j * m_noDistributions + dist;
        a_oldDistribution(haloCell(n, j), dist) = recBuf[recBuffer + offset];
      }
    }
    if(m_isThermal) {
      recBuffer += m_noDistributions * noHaloCells(n);
      for(MInt j = 0; j < noHaloCells(n); j++) {
        for(MInt dist = 0; dist < m_noDistributions; dist++) {
          MInt offset = j * m_noDistributions + dist;
          a_oldDistributionThermal(haloCell(n, j), dist) = recBuf[recBuffer + offset];
        }
      }
    }
    if(m_isTransport) {
      recBuffer += m_noDistributions * noHaloCells(n);
      for(MInt j = 0; j < noHaloCells(n); j++) {
        for(MInt dist = 0; dist < m_noDistributions; dist++) {
          MInt offset = j * m_noDistributions + dist;
          a_oldDistributionTransport(haloCell(n, j), dist) = recBuf[recBuffer + offset];
        }
      }
    }
  }
}

/** \brief This function stores the flow information of the cells
 *        such as variables and attributes like u_velocity,density,etc.,
 *
 * \author Andreas Lintermann, Konstantin Froehlich
 * \date 13.02.2012
 * In contrast to saveUvwOnly(const MChar* fileName) this function also stores the
 * information for restart in parallel.
 * The attribute 'name' of the variables are set to their according meaning.
 *
 * \param[in] fileName the name of the file to write to
 **/
template <MInt nDim>
void LbSolver<nDim>::saveRestartWithDistributionsPar(const MChar* fileName, const MChar* gridInputFileName,
                                                     MInt* recalcIdTree) {
  TRACE();
  if constexpr(nDim != 2 && nDim != 3) {
    cerr << " In global function saveRestartWithDistributionsPar: wrong number of dimensions !" << endl;
    exit(0);
    return;
  }

  using namespace maia::parallel_io;
  ParallelIo parallelIo(fileName, PIO_REPLACE, mpiComm());
  parallelIo.defineScalar(PIO_INT, "noCells");

  // total #of cells in the grid == global cell index of the last cell in the last process:
  const MInt totalNoCells = domainOffset(noDomains());

  //--specify helper functions
  auto defFloatArray = [&](const MString arrayName, const MString varName, const MInt length) {
    parallelIo.defineArray(PIO_FLOAT, arrayName, length);
    parallelIo.setAttribute(varName, "name", arrayName);
  };
  auto defIntArray = [&](const MString arrayName, const MString varName, const MInt length) {
    parallelIo.defineArray(PIO_INT, arrayName, length);
    parallelIo.setAttribute(varName, "name", arrayName);
  };
  // function: define macroscopic variables
  auto defineMacroscopicVariables = [&](const MString name, const MString prefix, const MBool old) {
    // velocities
    const MString velNames[3] = {"u", "v", "w"};
    for(MInt d = 0; d != nDim; ++d) {
      defFloatArray(name + std::to_string(d), prefix + velNames[d], totalNoCells);
    }
    // density
    defFloatArray(name + std::to_string(nDim), prefix + "rho", totalNoCells);
    // temperature
    if(!m_isTransport && m_isThermal) {
      defFloatArray(name + std::to_string(nDim + 1), prefix + "t", totalNoCells);
    }
    if(m_isTransport && !m_isThermal) {
      defFloatArray(name + std::to_string(nDim + 1), prefix + "c", totalNoCells);
    }
    if(m_isTransport && m_isThermal) {
      defFloatArray(name + std::to_string(nDim + 1), prefix + "t", totalNoCells);
      defFloatArray(name + std::to_string(nDim + 2), prefix + "c", totalNoCells);
    }
    if(m_writeLsData && !old) {
      if(m_isThermal) {
        defIntArray(name + std::to_string(nDim + 2), prefix + "isActive", totalNoCells);
        for(MInt set = 0; set < m_maxNoSets; set++) {
          defFloatArray(name + std::to_string(nDim + 3 + set), prefix + "G_" + std::to_string(set), totalNoCells);
        }
        for(MInt set = 0; set < m_maxNoSets; set++) {
          defIntArray(name + std::to_string(nDim + 3 + m_maxNoSets + set), prefix + "Body_" + std::to_string(set),
                      totalNoCells);
        }
      } else {
        defIntArray(name + std::to_string(nDim + 1), prefix + "isActive", totalNoCells);
        for(MInt set = 0; set < m_maxNoSets; set++) {
          defFloatArray(name + std::to_string(nDim + 2 + set), prefix + "G_" + std::to_string(set), totalNoCells);
        }
        for(MInt set = 0; set < m_maxNoSets; set++) {
          defIntArray(name + std::to_string(nDim + 2 + m_maxNoSets + set), prefix + "Body_" + std::to_string(set),
                      totalNoCells);
        }
      }
    }
    if(m_isEELiquid) {
      defFloatArray(name + std::to_string(nDim + 1), prefix + "alphaG", totalNoCells);
    }
  };

  auto defineLocalViscosity = [&](const MString& name) { defFloatArray(name, name, totalNoCells); };

  // function: define distributions
  auto defineDistributions = [&, totalNoCells](const MString name) {
    for(MInt j = 0; j != m_noDistributions; ++j) {
      defFloatArray(name + std::to_string(j), name + std::to_string(j), totalNoCells);
      if(m_isThermal) {
        const MString thermalName = name + "Thermal" + std::to_string(j);
        defFloatArray(thermalName, thermalName, totalNoCells);
      }
      if(m_isTransport) {
        const MString transportName = name + "Transport" + std::to_string(j);
        defFloatArray(transportName, transportName, totalNoCells);
      }
    }
  };
  // function: write macroscopic variable
  auto writeMacroscopicVariable = [&](const MInt index, const MInt suffix) {
    MFloatScratchSpace tmp(noInternalCells(), AT_, "tmp");
    for(MInt i = 0; i < noInternalCells(); ++i) {
      tmp[i] = a_variable(recalcIdTree != nullptr ? recalcIdTree[i] : i, index);
    }
    parallelIo.writeArray(tmp.getPointer(), "variables" + std::to_string(suffix));
  };
  // function: write old macroscopic variable
  auto writeMacroscopicOldVariable = [&](const MInt index, const MInt suffix) {
    MFloatScratchSpace tmp(noInternalCells(), AT_, "tmp");
    for(MInt i = 0; i < noInternalCells(); ++i) {
      tmp[i] = a_oldVariable(recalcIdTree != nullptr ? recalcIdTree[i] : i, index);
    }
    parallelIo.writeArray(tmp.getPointer(), "oldVariables" + std::to_string(suffix));
  };
  // function: write old macroscopic variable
  auto writeMacroscopicActiveState = [&](const MInt suffix) {
    MIntScratchSpace tmp(noInternalCells(), AT_, "tmp");
    for(MInt i = 0; i < noInternalCells(); ++i) {
      tmp[i] = (MInt)a_isActive(recalcIdTree != nullptr ? recalcIdTree[i] : i);
    }
    parallelIo.writeArray(tmp.getPointer(), "variables" + std::to_string(suffix));
  };
  auto writeLevelSet = [&](const MInt set, const MInt suffix) {
    MFloatScratchSpace tmp(noInternalCells(), AT_, "tmp");
    for(MInt i = 0; i < noInternalCells(); i++) {
      tmp[i] = a_levelSetFunctionMB(recalcIdTree != nullptr ? recalcIdTree[i] : i, set);
    }
    parallelIo.writeArray(tmp.getPointer(), "variables" + std::to_string(suffix));
  };
  auto writeBodyId = [&](const MInt set, const MInt suffix) {
    MIntScratchSpace tmp(noInternalCells(), AT_, "tmp");
    for(MInt i = 0; i < noInternalCells(); i++) {
      tmp[i] = a_associatedBodyIds(recalcIdTree != nullptr ? recalcIdTree[i] : i, set);
    }
    parallelIo.writeArray(tmp.getPointer(), "variables" + std::to_string(suffix));
  };
  auto writeViscosity = [&]() {
    MFloatScratchSpace tmp(noInternalCells(), AT_, "tmp");
    for(MInt i = 0; i < noInternalCells(); ++i) {
      tmp[i] = a_nu(recalcIdTree != nullptr ? recalcIdTree[i] : i);
    }
    parallelIo.writeArray(tmp.getPointer(), "viscosity");
  };

  //--define arrays
  defineMacroscopicVariables("variables", "", false);
  defineMacroscopicVariables("oldVariables", "old_", true);
  if(m_nonNewtonian) defineLocalViscosity("viscosity");
  defineDistributions("distributions");
  defineDistributions("oldDistributions");
  m_bndCnd->bcDataWriteRestartHeader(parallelIo);

  //--define global attributes
  parallelIo.setAttribute(solverId(), "solverId");
  parallelIo.setAttribute(gridInputFileName, "gridFile", "");
  parallelIo.setAttribute(m_time, "time");
  parallelIo.setAttribute(globalTimeStep, "globalTimeStep");
  parallelIo.setAttribute(m_Re, "ReynoldsNumber");
  parallelIo.setAttribute(m_Ma, "MachNumber");

  if(m_controlVelocity) {
    parallelIo.setAttribute(m_velocityControl.lastGlobalAvgV, "lastGlobalAvgV");
    parallelIo.setAttribute(m_velocityControl.previousError, "previousError");
    parallelIo.setAttribute(m_velocityControl.integratedError, "integratedError");
    parallelIo.setAttribute(m_velocityControl.derivedError, "derivedError");
    for(MInt i = 0; i < nDim; i++) {
      parallelIo.setAttribute(m_volumeAccel[i], "volumeAcceleration_" + std::to_string(i));
    }
  }

  //--write scalars
  parallelIo.writeScalar(totalNoCells, "noCells");

  //--write arrays
  // Set file offsets (first globalId and #of cells to be written by this process):
  const MPI_Offset firstGlobalId = domainOffset(domainId());
  const MPI_Offset localNoCells = noInternalCells();
  parallelIo.setOffset(localNoCells, firstGlobalId);

  // Macroscopic Variables
  MInt suffix = 0;
  writeMacroscopicVariable(PV->U, suffix++);
  writeMacroscopicVariable(PV->V, suffix++);
  if constexpr(nDim == 3) writeMacroscopicVariable(PV->W, suffix++);
  writeMacroscopicVariable(PV->RHO, suffix++);
  if(m_isThermal) writeMacroscopicVariable(PV->T, suffix++);
  if(m_isTransport) writeMacroscopicVariable(PV->C, suffix++);
  if(m_writeLsData) {
    writeMacroscopicActiveState(suffix++);
    for(MInt set = 0; set < m_maxNoSets; set++) {
      writeLevelSet(set, suffix++);
    }
    for(MInt set = 0; set < m_maxNoSets; set++) {
      writeBodyId(set, suffix++);
    }
  }
  if(m_isEELiquid) {
    MFloatScratchSpace tmp_alphaG(noInternalCells(), AT_, "tmp_alphaG");
    for(MInt i = 0; i < noInternalCells(); ++i) {
      const MInt id = recalcIdTree != nullptr ? recalcIdTree[i] : i;
      tmp_alphaG[i] = a_alphaGas(id);
    }
    parallelIo.writeArray(tmp_alphaG.getPointer(), "variables" + std::to_string(suffix++));
  }
  if(m_nonNewtonian) writeViscosity();

  suffix = 0;
  writeMacroscopicOldVariable(PV->U, suffix++);
  writeMacroscopicOldVariable(PV->V, suffix++);
  if constexpr(nDim == 3) writeMacroscopicOldVariable(PV->W, suffix++);
  writeMacroscopicOldVariable(PV->RHO, suffix++);
  if(m_isThermal) writeMacroscopicOldVariable(PV->T, suffix++);
  if(m_isTransport) writeMacroscopicOldVariable(PV->C, suffix++);
  // Distributions
  for(MInt j = 0; j < m_noDistributions; j++) {
    const MString distributions = "distributions" + to_string(j);

    MFloatScratchSpace tmp_distributions(noInternalCells(), AT_, "tmp_distributions");
    for(MInt i = 0; i < noInternalCells(); ++i)
      tmp_distributions[i] = a_distribution(recalcIdTree != nullptr ? recalcIdTree[i] : i, j);
    parallelIo.writeArray(tmp_distributions.getPointer(), distributions);
  }

  // old Distributions
  for(MInt j = 0; j < m_noDistributions; j++) {
    const MString oldDistributions = "oldDistributions" + to_string(j);

    MFloatScratchSpace tmp_oldDistributions(noInternalCells(), AT_, "tmp_oldDistributions");
    for(MInt i = 0; i < noInternalCells(); ++i)
      tmp_oldDistributions[i] = a_oldDistribution(recalcIdTree != nullptr ? recalcIdTree[i] : i, j);
    parallelIo.writeArray(tmp_oldDistributions.getPointer(), oldDistributions);
  }

  if(m_isThermal) {
    // Thermal Distributions
    for(MInt j = 0; j < m_noDistributions; j++) {
      const MString distributionsThermal = "distributionsThermal" + to_string(j);

      MFloatScratchSpace tmp_distributionsThermal(noInternalCells(), AT_, "tmp_distributionsThermal");
      for(MInt i = 0; i < noInternalCells(); ++i)
        tmp_distributionsThermal[i] = a_distributionThermal(recalcIdTree != nullptr ? recalcIdTree[i] : i, j);
      parallelIo.writeArray(tmp_distributionsThermal.getPointer(), distributionsThermal);
    }

    // old Thermal Distributions
    for(MInt j = 0; j < m_noDistributions; j++) {
      const MString oldDistributionsThermal = "oldDistributionsThermal" + to_string(j);

      MFloatScratchSpace tmp_oldDistributionsThermal(noInternalCells(), AT_, "tmp_oldDistributionsThermal");
      for(MInt i = 0; i < noInternalCells(); ++i)
        tmp_oldDistributionsThermal[i] = a_oldDistributionThermal(recalcIdTree != nullptr ? recalcIdTree[i] : i, j);
      parallelIo.writeArray(tmp_oldDistributionsThermal.getPointer(), oldDistributionsThermal);
    }
  }

  if(m_isTransport) {
    // Transport Distributions
    for(MInt j = 0; j < m_noDistributions; j++) {
      const MString distributionsTransport = "distributionsTransport" + to_string(j);

      MFloatScratchSpace tmp_distributionsTransport(noInternalCells(), AT_, "tmp_distributionsTransport");
      for(MInt i = 0; i < noInternalCells(); ++i)
        tmp_distributionsTransport[i] = a_distributionTransport(recalcIdTree != nullptr ? recalcIdTree[i] : i, j);
      parallelIo.writeArray(tmp_distributionsTransport.getPointer(), distributionsTransport);
    }

    // old Thermal Distributions
    for(MInt j = 0; j < m_noDistributions; j++) {
      const MString oldDistributionsTransport = "oldDistributionsTransport" + to_string(j);

      MFloatScratchSpace tmp_oldDistributionsTransport(noInternalCells(), AT_, "tmp_oldDistributionsTransport");
      for(MInt i = 0; i < noInternalCells(); ++i)
        tmp_oldDistributionsTransport[i] = a_oldDistributionTransport(recalcIdTree != nullptr ? recalcIdTree[i] : i, j);
      parallelIo.writeArray(tmp_oldDistributionsTransport.getPointer(), oldDistributionsTransport);
    }
  }
  // Former variables
  m_bndCnd->bcDataWriteRestartData(parallelIo);
}

/**
 * \brief Saving the current solution and write new grid file if needed
 *
 * \author Julian Vorspohl <j.vorspohl@aia.rwth-aachen.de>
 */
template <MInt nDim>
void LbSolver<nDim>::saveOutput() {
  TRACE();

  stringstream PvFileName;
  PvFileName << outputDir() << "PV_" << getIdentifier() << globalTimeStep << ParallelIo::fileExt();
  if(domainId() == 0)
    std::cerr << "Writing output file: PV_" << getIdentifier() << globalTimeStep << " at m_time " << this->m_time
              << std::endl;


  // If grid is adapted, write new grid file and compute recalcIdTree.
  // Different ids in solver/grid and new grid file are mapped by recalcIdTree.
  if(m_adaptationSinceLastSolution) {
    m_reinitFileName = "grid" + std::to_string(globalTimeStep) + ".Netcdf";
    m_reinitFilePath = outputDir() + m_reinitFileName;
    saveAdaptedGridFile(this->m_recalcIds);
  }
  std::vector<MInt> recalcCellIdsSolver(0);
  MInt noCells;
  MInt noInternalCellIds;
  std::vector<MInt> reorderedCellIds(0);
  this->calcRecalcCellIdsSolver(this->m_recalcIds, noCells, noInternalCellIds, recalcCellIdsSolver, reorderedCellIds);
  m_adaptationSinceLastSolution = false;
  saveUVWRhoTOnlyPar(PvFileName.str().c_str(), m_reinitFileName.c_str(), recalcCellIdsSolver.data());
}

/** \brief Save adapted gridFile
 *
 * \details: Move from save output in lbsolverdxqy to lbsolver because
 *           this class is friend of cartesian grid.
 *           New grid file is written and change in cellIds is stored in m_recalcIds
 *
 * \author Philipp Brokof, Moritz Waldmann
 **/
template <MInt nDim>
void LbSolver<nDim>::saveAdaptedGridFile(MInt* const p_recalcCellIds) {
  // Reset p_recalcCellIds
  for(MInt i = 0; i < this->maxNoGridCells(); i++) {
    p_recalcCellIds[i] = -1;
  }

  // Write new Grid and get new recalcIdTree
  grid().raw().saveGrid(this->m_reinitFilePath.c_str(), p_recalcCellIds);
}

/** \brief This function stores all necessary variables and derivative information
 *
 * \author Andreas Lintermann
 * \date 10.03.2016
 *
 * This method uses central differences for the calulation of the derivatives.
 *
 * \param[in] fileName name of the file to write
 * \param[in] gridInputFileName name of the according grid file, will be written as an attribute
 * \param[in] recalcIdTree the new ordering, save only these cells
 *
 **/
template <MInt nDim>
void LbSolver<nDim>::saveUVWRhoTOnlyPar(const MChar* fileName, const MChar* gridInputFileName, MInt* recalcIdTree) {
  TRACE();

  if(m_saveDerivatives) {
    exchange(1);
  }

  using namespace maia::parallel_io;
  ASSERT(nDim == 2 || nDim == 3, "wrong number of dimensions!");

  ParallelIo parallelIo(fileName, PIO_REPLACE, mpiComm());
  parallelIo.defineScalar(PIO_INT, "noCells");
  parallelIo.setAttribute(solverId(), "solverId");
  // total #of cells in the grid == global cell index of the last cell in the last process:
  const MInt totalNoCells = domainOffset(noDomains());

  // defines the array and creates the corresponding attributes
  auto defFloatArray = [&](const MString& arrayName, const MString& varName) {
    parallelIo.defineArray(PIO_FLOAT, arrayName, totalNoCells);
    parallelIo.setAttribute(varName, "name", arrayName);
  };

  auto defineMacroscopicVariables = [&](const MString& name, const MString& prefix) {
    MInt counter = 0;
    // define variables
    if constexpr(nDim == 2) {
      const MString namesStd[3] = {"u", "v", "rho"};

      for(; counter < 3; counter++)
        defFloatArray(name + std::to_string(counter), prefix + namesStd[counter]);
    } else {
      const MString namesStd[4] = {"u", "v", "w", "rho"};

      for(; counter < 4; counter++)
        defFloatArray(name + std::to_string(counter), prefix + namesStd[counter]);
    }

    if(m_isThermal) {
      const MString temperature = "t";
      defFloatArray(name + std::to_string(counter), prefix + temperature);
      counter++;
    }

    if(m_isTransport) {
      const MString concentration = "c";
      defFloatArray(name + std::to_string(counter), prefix + concentration);
      counter++;
    }
    if(m_particleMomentumCoupling && m_saveExternalForces) {
      const MString externalForcing[3] = {"F_ext_x", "F_ext_y", "F_ext_z"};
      for(MInt d = 0; d < nDim; d++, counter++)
        defFloatArray(name + std::to_string(counter), prefix + externalForcing[d]);
    }

    if(m_isEELiquid) {
      const MString EELiquid = "alphaG";
      defFloatArray(name + std::to_string(counter), prefix + EELiquid);
      counter++;
    }

    if(m_saveDerivatives) {
      if constexpr(nDim == 2) {
        // define derivatives:
        const MString names[10] = {"du/dx",   "du/dy",       "dv/dx",       "dv/dy",           "dp_t/dx",
                                   "dp_t/dy", "vorticity_z", "Q-criterion", "Delta-criterion", "Lambda2"};

        for(MInt d = 0; d < 10; d++, counter++)
          defFloatArray(name + std::to_string(counter), prefix + names[d]);
      }

      else {
        // define derivatives:
        const MString names[19] = {"du/dx",         "du/dy",       "du/dz",           "dv/dx",       "dv/dy",
                                   "dv/dz",         "dw/dx",       "dw/dy",           "dw/dz",       "dp_t/dx",
                                   "dp_t/dy",       "dp_t/dz",     "vorticity_x",     "vorticity_y", "vorticity_z",
                                   "vorticity_mag", "Q-criterion", "Delta-criterion", "Lambda2"};

        for(MInt d = 0; d < 19; d++, counter++)
          defFloatArray(name + std::to_string(counter), prefix + names[d]);

        if(m_isThermal) {
          const MString namesT[3] = {"dT/dx", "dT/dy", "dT/dz"};
          for(MInt d = 0; d < nDim; d++, counter++)
            defFloatArray(name + std::to_string(counter), prefix + namesT[d]);
        }

        if(m_isTransport) {
          const MString namesC[3] = {"dC/dx", "dC/dy", "dC/dz"};
          for(MInt d = 0; d < nDim; d++, counter++)
            defFloatArray(name + std::to_string(counter), prefix + namesC[d]);
        }
      }
    }
  };

  defineMacroscopicVariables("variables", "");

  // Creating global attributes
  MString gridFileName = gridInputFileName;
  parallelIo.setAttribute(gridFileName, "gridFile", "");
  parallelIo.setAttribute(m_time, "time");
  parallelIo.setAttribute(globalTimeStep, "globalTimeStep");
  parallelIo.setAttribute(m_Re, "ReynoldsNumber");
  parallelIo.setAttribute(m_Ma, "MachNumber");

  if(m_controlVelocity) {
    parallelIo.setAttribute(m_velocityControl.lastGlobalAvgV, "lastGlobalAvgV");
    parallelIo.setAttribute(m_velocityControl.previousError, "previousError");
    parallelIo.setAttribute(m_velocityControl.integratedError, "integratedError");
    parallelIo.setAttribute(m_velocityControl.derivedError, "derivedError");
    for(MInt i = 0; i < nDim; i++) {
      parallelIo.setAttribute(m_volumeAccel[i], "volumeAcceleration_" + std::to_string(i));
    }
  }

  // number of Cells.
  parallelIo.writeScalar(totalNoCells, "noCells");

  // Set file offsets (first globalId and #of cells to be written by this process):
  const MPI_Offset firstGlobalId = domainOffset(domainId());
  const MPI_Offset localNoCells = noInternalCells();
  parallelIo.setOffset(localNoCells, firstGlobalId);

  MString name = "variables";

  MInt start_derivatives = m_noVariables;
  // Write macroscopic variables:
  {
    MInt counter = 0;
    ScratchSpace<MFloat> buffer(localNoCells, FUN_, "buffer");
    // write velocities and density:
    for(MInt d = 0; d != nDim + 1; ++d) {
      for(MInt cellId = 0; cellId < localNoCells; cellId++) {
        buffer[cellId] = a_variable(recalcIdTree != nullptr ? recalcIdTree[cellId] : cellId, d);
      }
      parallelIo.writeArray(&buffer[0], name + std::to_string(counter++));
    }

    // write temperature:
    if(m_isThermal && !m_isTransport) {
      for(MInt cellId = 0; cellId < localNoCells; cellId++) {
        buffer[cellId] = a_variable(recalcIdTree != nullptr ? recalcIdTree[cellId] : cellId, PV->T);
      }
      parallelIo.writeArray(&buffer[0], name + std::to_string(counter++));
    }
    if(m_isTransport && !m_isThermal) {
      for(MInt cellId = 0; cellId < localNoCells; cellId++) {
        buffer[cellId] = a_variable(recalcIdTree != nullptr ? recalcIdTree[cellId] : cellId, PV->C);
      }
      parallelIo.writeArray(&buffer[0], name + std::to_string(nDim + 1));
    }
    if(m_isTransport && m_isThermal) {
      for(MInt cellId = 0; cellId < localNoCells; cellId++) {
        buffer[cellId] = a_variable(recalcIdTree != nullptr ? recalcIdTree[cellId] : cellId, PV->T);
      }
      parallelIo.writeArray(&buffer[0], name + std::to_string(nDim + 1));

      for(MInt cellId = 0; cellId < localNoCells; cellId++) {
        buffer[cellId] = a_variable(recalcIdTree != nullptr ? recalcIdTree[cellId] : cellId, PV->C);
      }
      parallelIo.writeArray(&buffer[0], name + std::to_string(nDim + 2));
    }

    // write external Force:
    if(m_particleMomentumCoupling && m_saveExternalForces) {
      for(MInt d = 0; d < nDim; d++) {
        for(MInt cellId = 0; cellId < localNoCells; cellId++) {
          buffer[cellId] = a_externalForces(recalcIdTree != nullptr ? recalcIdTree[cellId] : cellId, d);
        }
        parallelIo.writeArray(&buffer[0], name + std::to_string(counter++));
      }
      start_derivatives += nDim;
    }

    if(m_isEELiquid) {
      for(MInt cellId = 0; cellId < localNoCells; cellId++) {
        buffer[cellId] = a_alphaGas(recalcIdTree != nullptr ? recalcIdTree[cellId] : cellId);
      }
      parallelIo.writeArray(&buffer[0], name + std::to_string(counter++));
      start_derivatives += 1;
    }
  }

  // write derivatives
  if(m_saveDerivatives) {
    // d1 = du, dv, dw;    d2 = dx, dy, dz
    ScratchSpace<MFloat> derivatives(nDim * nDim + nDim, localNoCells, FUN_, "derivatives");
    for(MInt d1 = 0; d1 < nDim; ++d1)
      for(MInt d2 = 0; d2 < nDim; ++d2)
        for(MInt cellId = 0; cellId < localNoCells; cellId++)
          derivatives(d1 * nDim + d2, cellId) = calculateDerivative(cellId, d1, d2);

    for(MInt d = 0; d < nDim; d++) {
      for(MInt cellId = 0; cellId < localNoCells; cellId++) {
        if(m_calcTotalPressureGradient != 0) {
          derivatives(nDim * nDim + d, cellId) = calculatePressureDerivative(cellId, d);
        } else {
          derivatives(nDim * nDim + d, cellId) = calculateDerivative(cellId, PV->RHO, d);
        }
      }
    }

    // write derivatives to disk
    for(MInt d1 = 0; d1 < nDim; ++d1)
      for(MInt d2 = 0; d2 < nDim; ++d2)
        parallelIo.writeArray(&derivatives(d1 * nDim + d2, 0),
                              name + std::to_string(start_derivatives + d1 * nDim + d2));

    for(MInt d = 0; d < nDim; d++) {
      parallelIo.writeArray(&derivatives(nDim * nDim + d, 0),
                            name + std::to_string(start_derivatives + nDim * nDim + d));
    }

    MInt start_vor;
    MInt start_vor_mag;
    MInt start_q;
    MInt start_d;
    MInt start_lambda2;

    if constexpr(nDim == 2) {
      start_vor = start_derivatives + nDim * nDim + nDim;
      start_q = start_vor + 1;
      start_d = start_q + 1;
      start_lambda2 = start_d + 1;
    } else {
      start_vor = start_derivatives + nDim * nDim + nDim;
      start_vor_mag = start_vor + nDim;
      start_q = start_vor_mag + 1;
      start_d = start_q + 1;
      start_lambda2 = start_d + 1;
    }
    // all vorticity components
    {// create vorticity elements and write to disk
     // also calculate the magnitude of the vorticity
     if constexpr(nDim == 2) {
       ScratchSpace<MFloat> vorticity_c(localNoCells, FUN_, "vorticity_c");
       for(MInt cellId = 0; cellId < localNoCells; cellId++) {
         vorticity_c[cellId] = derivatives(2, cellId) - derivatives(1, cellId);
       }
       parallelIo.writeArray(&vorticity_c[0], name + std::to_string(start_vor));
     }
     if constexpr(nDim == 3) {
       ScratchSpace<MFloat> vorticity_mag(localNoCells, FUN_, "vorticity_mag");
       for(MInt cellId = 0; cellId < localNoCells; cellId++) {
         vorticity_mag[cellId] = 0.0;
       }

       for(MInt d = 0; d < nDim; ++d) {
         MInt v_c1 = (d + 1) % (nDim);
         MInt v_c2 = (d + 2) % (nDim);
         ScratchSpace<MFloat> vorticity_c(localNoCells, FUN_, "vorticity_c");
         for(MInt cellId = 0; cellId < localNoCells; cellId++) {
           vorticity_c[cellId] = derivatives(v_c2 * nDim + v_c1, cellId) - derivatives(v_c1 * nDim + v_c2, cellId);
           vorticity_mag[cellId] += vorticity_c[cellId] * vorticity_c[cellId];
         }

         parallelIo.writeArray(&vorticity_c[0], name + std::to_string(start_vor + d));
       }

       // write vorticity magnitude
       for(MInt cellId = 0; cellId < localNoCells; cellId++) {
         vorticity_mag[cellId] = sqrt(vorticity_mag[cellId]);
       }

       parallelIo.writeArray(&vorticity_mag[0], name + std::to_string(start_vor_mag));
     }
}
    // q-criterion and delta-criterion
    {
      // first q-criterion
      ScratchSpace<MFloat> q_criterion(localNoCells, FUN_, "q_criterion");
      MFloatScratchSpace strain(nDim, nDim, FUN_, "strain");
      MFloatScratchSpace vorticity(nDim, nDim, FUN_, "vorticity");

      for(MInt cellId = 0; cellId < localNoCells; cellId++) {
        getStrainTensor(derivatives, cellId, strain);
        getVorticityTensor(derivatives, cellId, vorticity);

        MFloat strainFrobSq = maia::math::frobeniusMatrixNormSquared(strain, nDim, nDim);
        MFloat vorticityFrobSq = maia::math::frobeniusMatrixNormSquared(vorticity, nDim, nDim);

        q_criterion[cellId] = 0.5 * (vorticityFrobSq - strainFrobSq);
      }

      parallelIo.writeArray(&q_criterion[0], name + std::to_string(start_q));

      // second delta-criterion (requires q-criterion)
      {
        ScratchSpace<MFloat> d_criterion(localNoCells, FUN_, "d_criterion");

        for(MInt cellId = 0; cellId < localNoCells; cellId++) {
          std::array<std::array<MFloat, nDim>, nDim> dyad;

          for(MInt i = 0; i < nDim; i++) {
            for(MInt j = 0; j < nDim; j++) {
              dyad[i][j] = derivatives(j * nDim + i, cellId);
            }
          }
          MFloat det = maia::math::determinant(dyad);
          d_criterion[cellId] = pow((q_criterion[cellId] / 3.0), 3.0) + (0.5 * det) * (0.5 * det);
        }
        parallelIo.writeArray(&d_criterion[0], name + std::to_string(start_d));
      }

      // third lambda 2
      {
        ScratchSpace<MFloat> lambda_2(localNoCells, FUN_, "lambda_2");
        MFloatScratchSpace strain_sq(nDim, nDim, FUN_, "strain_sq");
        MFloatScratchSpace vorticity_sq(nDim, nDim, FUN_, "vorticity_sq");
        MFloatScratchSpace sv(nDim, nDim, FUN_, "sv");

        for(MInt cellId = 0; cellId < localNoCells; cellId++) {
          getStrainTensor(derivatives, cellId, strain);
          getVorticityTensor(derivatives, cellId, vorticity);

          maia::math::multiplyMatricesSq(strain, strain, strain_sq, nDim);
          maia::math::multiplyMatricesSq(vorticity, vorticity, vorticity_sq, nDim);

          maia::math::addMatrices(strain_sq, vorticity_sq, sv, nDim, nDim);

          if constexpr(nDim == 3) {
            MFloat matrix[3][3];
            for(MInt i = 0; i < nDim; i++) {
              for(MInt j = 0; j < nDim; j++) {
                matrix[i][j] = sv(i, j);
              }
            }

            MFloat evs[3]{};
            maia::math::calcEigenValues(matrix, evs);
            MInt l2_index = 0;
            for(MInt i = 0; i < nDim; i++) {
              if((evs[i] > evs[((i == 0) ? 0 : (i - 1)) % nDim] && evs[i] < evs[(i + 1) % nDim])
                 || (evs[i] > evs[(i + 1) % nDim] && evs[i] < evs[((i == 0) ? 0 : (i - 1)) % nDim])) {
                l2_index = i;
                break;
              }
            }
            lambda_2[cellId] = evs[l2_index];
          } else {
            TERMM(-1, "Code is not correct for 2D applications");
          }
        }
        parallelIo.writeArray(&lambda_2[0], name + std::to_string(start_lambda2));
      }
      if(m_isThermal && !m_isTransport) {
        ScratchSpace<MFloat> derivativesT(nDim, localNoCells, FUN_, "derivativesT");

        for(MInt d = 0; d < nDim; d++) {
          for(MInt cellId = 0; cellId < localNoCells; cellId++)
            derivativesT(d, cellId) = calculateDerivative(cellId, PV->T, d);
        }
        for(MInt d = 0; d < nDim; d++) {
          parallelIo.writeArray(&derivativesT(d, 0), name + std::to_string(start_lambda2 + 1 + d));
        }
      }
      if(m_isTransport && !m_isThermal) {
        ScratchSpace<MFloat> derivativesC(nDim, localNoCells, FUN_, "derivativesC");

        for(MInt d = 0; d < nDim; d++) {
          for(MInt cellId = 0; cellId < localNoCells; cellId++)
            derivativesC(d, cellId) = calculateDerivative(cellId, PV->C, d);
        }
        for(MInt d = 0; d < nDim; d++) {
          parallelIo.writeArray(&derivativesC(d, 0), name + std::to_string(start_lambda2 + 1 + d));
        }
      }
      if(m_isTransport && m_isThermal) {
        ScratchSpace<MFloat> derivativesT(nDim, localNoCells, FUN_, "derivativesT");

        for(MInt d = 0; d < nDim; d++) {
          for(MInt cellId = 0; cellId < localNoCells; cellId++)
            derivativesT(d, cellId) = calculateDerivative(cellId, PV->T, d);
        }
        for(MInt d = 0; d < nDim; d++) {
          parallelIo.writeArray(&derivativesT(d, 0), name + std::to_string(start_lambda2 + 1 + d));
        }

        ScratchSpace<MFloat> derivativesC(nDim, localNoCells, FUN_, "derivativesC");

        for(MInt d = 0; d < nDim; d++) {
          for(MInt cellId = 0; cellId < localNoCells; cellId++)
            derivativesC(d, cellId) = calculateDerivative(cellId, PV->C, d);
        }
        for(MInt d = 0; d < nDim; d++) {
          parallelIo.writeArray(&derivativesC(d, 0), name + std::to_string(start_lambda2 + 2 + d));
        }
      }
    }
  }
}

/**
 * \brief calculate ReLambda and UrmsInit for isotropic Turbulence init
 *
 * \author Johannes Grafen
 * \date 20.02.2022
 *
 **/
template <MInt nDim>
void LbSolver<nDim>::getReLambdaAndUrmsInit() {
  TRACE();
  const MFloat Lb = m_referenceLength;
  // cout << "Lb: " << Lb << endl;
  MFloat kpRatio = F4; // peak wave number of prescribed spectrum
  kpRatio = Context::getSolverProperty<MFloat>("kpRatio", this->m_solverId, AT_, &kpRatio);
  if(Context::propertyExists("UrmsInit")) {
    m_UrmsInit = Context::getSolverProperty<MFloat>("UrmsInit", this->m_solverId, AT_);
    m_ReLambda = F1 / (F2 * PI * kpRatio) * sqrt(F5 / F6) * m_Re; // assuming that Re = Lb * m_UrmsInit / nu
    /*see literature: Lattice Boltzmann simulation of turbulent flow laden with finite-size particles
    Hui Gao, Hui Li, Lian-Ping Wang*/
  } else if(Context::propertyExists("ReLambda")) {
    m_ReLambda = Context::getSolverProperty<MFloat>("ReLambda", this->m_solverId, AT_);
    m_UrmsInit = m_ReLambda * 2 * PI * kpRatio * m_nu * sqrt(F6 / F5) / Lb;
  } else {
    mTerm(1, AT_, "UrmsInit or ReLambda and nu must be specified in property file");
  }
}

template <MInt nDim>
void LbSolver<nDim>::computeFFTStatistics() {
  mTerm(1, AT_, "computeFFTStatistics not implemented for 2D");
}

/**
 * \brief adapted from FV_MB, compute statistics for isotropic Turbulence
 *
 * \author Johannes Grafen
 * \date 20.02.2022
 *
 **/
template <>
void LbSolver<3>::computeFFTStatistics() {
  TRACE();
  constexpr MInt nDim = 3;

  std::array<MFloat, nDim * 2> bBox;
  m_geometry->getBoundingBox(bBox.data());

  MInt noVel = nDim;

  // fft-domain dimensions
  // this holds the size of the domain in number of cells on lowest level
  const MInt fftLevel = grid().maxUniformRefinementLevel();
  if(fftLevel > grid().maxUniformRefinementLevel()) mTerm(1, AT_, "Non-isotropic mesh is not implemented, yet");
  if(fftLevel < minLevel()) mTerm(1, AT_, "Parents missing for fftLevel < minLevel()");

  const MFloat dx = c_cellLengthAtLevel(fftLevel);
  const MFloat dxeps = 0.1 * dx;
  MInt nx = (bBox[3] - bBox[0] + dxeps) / dx;
  MInt ny = (bBox[4] - bBox[1] + dxeps) / dx;
  MInt nz = (bBox[5] - bBox[2] + dxeps) / dx;

  // MFloatScratchSpace coupling(a_noCells(), nDim, AT_, "coupling");
  MFloatScratchSpace pVariables(a_noCells(), noVel, AT_, "pVariables");
  // MFloatScratchSpace cVariables(a_noCells(), m_noVars, AT_, "cVariables");
  MFloatScratchSpace oldVariables(a_noCells(), noVel, AT_, "oldVariables");
  /*MInt filterLvlLaplace = fftLevel;
  filterLvlLaplace = Context::getSolverProperty<MInt>("filterLvlLaplace", m_solverId, AT_, &filterLvlLaplace);
  */
  // MFloat couplingCheck = determineCoupling(coupling);
  for(MInt cellId = 0; cellId < a_noCells(); cellId++) {
    for(MInt i = 0; i < noVel; i++) {
      pVariables(cellId, i) = a_variable(cellId, i);
      // cVariables(cellId, i) = a_variable(cellId, i);
      oldVariables(cellId, i) = a_oldVariable(cellId, i);
    }
  }
  /* if(domainId() == 0)
    cerr << "coupling check " << couplingCheck << " " << m_couplingRate << " " << couplingCheck / (DX * DX * DX)
          << " " << m_couplingRate / (DX * DX * DX) << " " << DX * m_couplingRate / POW3(DX * m_UInfinity) << endl; */

  // Apply volumetric filter only for particle resolved simulations with local grid refinement
  /* if(maxRefinementLevel() != fftLevel) {
    for(MInt cellId = 0; cellId < a_noCells(); cellId++) {
      if(a_hasProperty(cellId, SolverCell::IsSplitCell)) continue;
      if(a_isPeriodic(cellId)) continue;
      if(m_bndryCellIds->a[cellId] < -1) continue;
      if(a_isHalo(cellId)) continue;
      if(a_level(cellId) != fftLevel) continue;
      if(!a_hasProperty(cellId, SolverCell::IsSplitChild) && c_noChildren(cellId) > 0) {
        reduceData(cellId, &coupling(0), nDim, false);
        reduceData(cellId, &(pVariables(0, 0)), m_noVars);
        reduceData(cellId, &(cVariables(0, 0)), m_noVars);
        reduceData(cellId, &(oldVariables(0, 0)), m_noVars);
      }
    }
  } */

  MInt noLocDat = 0;
  for(MInt cellId = 0; cellId < a_noCells(); cellId++) {
    if(a_isHalo(cellId)) continue;
    if(a_level(cellId) != fftLevel) continue;
    noLocDat++;
  }
  MFloatScratchSpace velDat(noLocDat, 3 * nDim, AT_, "velDat");
  MIntScratchSpace velPos(noLocDat, AT_, "velPos");

  noLocDat = 0;

  // set UrmsInit and ReLambda once after iterative init of rho
  if(globalTimeStep == 0 || globalTimeStep / m_fftInterval == 1) {
    getReLambdaAndUrmsInit();
    if(domainId() == 0) cerr << "UrmsInit: " << m_UrmsInit << endl << "ReLambdaInit: " << m_ReLambda << endl;
  }

  // const MFloat dt = dx * m_Ma * LBCS / m_UrmsInit;
  // const MFloat velocityFactor = m_UrmsInit * F1BCS / m_Ma;
  // const MFloat nu =  m_UrmsInit * nx / m_Re;
  /* if(domainId() == 0)
    cout << "viscosity for computeEnergySpectrum: " << m_nu << endl; */

  for(MInt cellId = 0; cellId < a_noCells(); cellId++) {
    if(a_isHalo(cellId)) continue;
    if(a_level(cellId) != fftLevel) continue;

    MFloat actualCellLength = c_cellLengthAtCell(cellId);
    MInt xPos = floor((F1B2 * nx + (a_coordinate(cellId, 0) - F1B2 * (actualCellLength)) / dx) + 0.1);
    MInt yPos = floor((F1B2 * ny + (a_coordinate(cellId, 1) - F1B2 * (actualCellLength)) / dx) + 0.1);
    MInt zPos = floor((F1B2 * nz + (a_coordinate(cellId, 2) - F1B2 * (actualCellLength)) / dx) + 0.1);
    if(xPos > nx - 1 || xPos < 0 || yPos > ny - 1 || yPos < 0 || zPos > nz - 1 || zPos < 0) {
      cerr << "ERROR: wrong array position!" << endl;
      cerr << "pos=" << xPos << ", " << yPos << ", " << zPos << endl;
      cerr << "coorda = (" << a_coordinate(cellId, 0) << ", " << a_coordinate(cellId, 1) << ", "
           << a_coordinate(cellId, 2) << ")" << endl;
      cerr << "actuallength=" << actualCellLength << endl;
      cerr << "minlevel=" << minLevel() << ", maxLevel=" << maxLevel() << endl;
      cerr << "lenght on level0=" << c_cellLengthAtLevel(0) << endl;
      mTerm(1, AT_, "Wrong array position");
    }

    // TODO: adjust for lb usage
    MInt pos = zPos + nz * (yPos + ny * xPos);
    for(MInt i = 0; i < nDim; i++) {
      // Transform velocity data from LB to physical space, because spectrum is computed in physical space
      velDat(noLocDat, i) = pVariables(cellId, PV->VV[i]);
      velDat(noLocDat, nDim + i) = a_externalForces(cellId, i); // cellVolume in lb units is 1
      // MFloat u1 = cVariables(cellId, CV->RHO_VV[i]) / cVariables(cellId, CV->RHO);
      // MFloat u0 = oldVariables(cellId, CV->RHO_VV[i]) / oldVariables(cellId, CV->RHO);
      MFloat u1 = pVariables(cellId, PV->VV[i]);
      MFloat u0 = oldVariables(cellId, PV->VV[i]);
      velDat(noLocDat, 2 * nDim + i) = ((u1 - u0) / F1); // rate of change beweeten two timeteps
    }
    velPos(noLocDat) = pos;
    noLocDat++;
  }

  maia::math::computeEnergySpectrum(velDat, velPos, 3 * nDim, noLocDat, nx, ny, nz, m_nu, mpiComm());
}

template <MInt nDim>
inline void LbSolver<nDim>::getStrainTensor(MFloatScratchSpace& derivatives, MInt cellId, MFloatScratchSpace& strain) {
  // TRACE();

  // diagonal
  for(MInt d1 = 0; d1 < nDim; d1++)
    strain(d1, d1) = derivatives(d1 * nDim + d1, cellId);

  // off-diagonal
  for(MInt d1 = 0; d1 < nDim - 1; d1++)
    for(MInt d2 = d1 + 1; d2 < nDim; d2++) {
      MInt a = d1 * nDim + d2;
      MInt b = d2 * nDim + d1;
      strain(d1, d2) = 0.5 * (derivatives(a, cellId) + derivatives(b, cellId));
      strain(d2, d1) = strain(d1, d2);
    }
}


template <MInt nDim>
inline void LbSolver<nDim>::getVorticityTensor(MFloatScratchSpace& derivatives, MInt cellId,
                                               MFloatScratchSpace& vorticity) {
  // TRACE();

  // diagonal
  for(MInt d1 = 0; d1 < nDim; d1++)
    vorticity(d1, d1) = 0.0;

  // off-diagonal
  for(MInt d1 = 0; d1 < nDim - 1; d1++)
    for(MInt d2 = d1 + 1; d2 < nDim; d2++) {
      MInt a = d1 * nDim + d2;
      MInt b = d2 * nDim + d1;
      vorticity(d1, d2) = 0.5 * (derivatives(a, cellId) - derivatives(b, cellId));
      vorticity(d2, d1) = 0.5 * (derivatives(b, cellId) - derivatives(a, cellId));
    }
}

/** \brief Calculates a spatial derivative for the total pressure
 *
 * \author Moritz Waldmann
 * \date 28.11.2019
 *
 * This function uses:
 *  - 4th order central differences if both neighbors and their neighbors in spaceDir are available
 *  - 3rd order forward or backward differences if both neighbors and one neighbors neighbor in spaceDir are available
 *  - 2nd order central differences if both neighbors in spaceDir are available
 *  - 1st order forward or backward differences if only one neighbor is available
 *  - zero gradient if no neighbor is available
 *
 * \param[in] cellId the cell id
 * \param[in] spaceDir the spatial direction
 * \return the derivative
 *
 **/
template <MInt nDim>
inline MFloat LbSolver<nDim>::calculatePressureDerivative(MInt cellId, MInt spaceDir) {
  // TRACE();
  MInt lbdir1 = 2 * spaceDir;
  MInt lbdir2 = lbdir1 + 1;

  MInt left = c_neighborId(cellId, lbdir1);
  MInt right = c_neighborId(cellId, lbdir2);
  MInt leftleft = (c_neighborId(cellId, lbdir1) > -1) ? c_neighborId(left, lbdir1) : -1;
  MInt rightright = (c_neighborId(cellId, lbdir2) > -1) ? c_neighborId(right, lbdir2) : -1;
  MFloat gradient = F0;

  MFloat pt = a_variable(cellId, PV->RHO) * F1B3
              + F1B2 / a_variable(cellId, PV->RHO)
                    * (a_variable(cellId, PV->U) * a_variable(cellId, PV->U)
                       + a_variable(cellId, PV->V) * a_variable(cellId, PV->V)
                       + a_variable(cellId, PV->W) * a_variable(cellId, PV->W));

  if((left > -1) && (right > -1)) {
    MFloat ptLeft =
        a_variable(left, PV->RHO) * F1B3
        + F1B2 / a_variable(left, PV->RHO)
              * (a_variable(left, PV->U) * a_variable(left, PV->U) + a_variable(left, PV->V) * a_variable(left, PV->V)
                 + a_variable(left, PV->W) * a_variable(left, PV->W));

    MFloat ptRight = a_variable(right, PV->RHO) * F1B3
                     + F1B2 / a_variable(right, PV->RHO)
                           * (a_variable(right, PV->U) * a_variable(right, PV->U)
                              + a_variable(right, PV->V) * a_variable(right, PV->V)
                              + a_variable(right, PV->W) * a_variable(right, PV->W));

    if(leftleft > -1 && rightright > -1) { // use central differences 4th order

      MFloat ptLeftLeft = a_variable(leftleft, PV->RHO) * F1B3
                          + F1B2 / a_variable(leftleft, PV->RHO)
                                * (a_variable(leftleft, PV->U) * a_variable(leftleft, PV->U)
                                   + a_variable(leftleft, PV->V) * a_variable(leftleft, PV->V)
                                   + a_variable(leftleft, PV->W) * a_variable(leftleft, PV->W));

      MFloat ptRightRight = a_variable(rightright, PV->RHO) * F1B3
                            + F1B2 / a_variable(rightright, PV->RHO)
                                  * (a_variable(rightright, PV->U) * a_variable(rightright, PV->U)
                                     + a_variable(rightright, PV->V) * a_variable(rightright, PV->V)
                                     + a_variable(rightright, PV->W) * a_variable(rightright, PV->W));

      gradient = (-ptRightRight + 8.0 * ptRight - 8.0 * ptLeft + ptLeftLeft) / 12.0;

    } else {
      if(leftleft > -1) { // backward differences 3nd order

        MFloat ptLeftLeft = a_variable(leftleft, PV->RHO) * F1B3
                            + F1B2 / a_variable(leftleft, PV->RHO)
                                  * (a_variable(leftleft, PV->U) * a_variable(leftleft, PV->U)
                                     + a_variable(leftleft, PV->V) * a_variable(leftleft, PV->V)
                                     + a_variable(leftleft, PV->W) * a_variable(leftleft, PV->W));

        gradient = (2.0 * ptRight + 3.0 * pt - 6.0 * ptLeft + ptLeftLeft) / 6.0;
      } else if(rightright > -1) { // forward differences 3nd order

        MFloat ptRightRight = a_variable(rightright, PV->RHO) * F1B3
                              + F1B2 / a_variable(rightright, PV->RHO)
                                    * (a_variable(rightright, PV->U) * a_variable(rightright, PV->U)
                                       + a_variable(rightright, PV->V) * a_variable(rightright, PV->V)
                                       + a_variable(rightright, PV->W) * a_variable(rightright, PV->W));

        gradient = (-ptRightRight + 6.0 * ptRight - 3.0 * pt - 2.0 * ptLeft) / 6.0;
      } else { // use central differences 2nd order
        gradient = (ptRight - ptLeft) / 2.0;
      }
    }
  } else { // use forward or backward differences 1st order
    // backward
    if(left > -1) {
      MFloat ptLeft =
          a_variable(left, PV->RHO) * F1B3
          + F1B2 / a_variable(left, PV->RHO)
                * (a_variable(left, PV->U) * a_variable(left, PV->U) + a_variable(left, PV->V) * a_variable(left, PV->V)
                   + a_variable(left, PV->W) * a_variable(left, PV->W));

      gradient = pt - ptLeft;
    }
    // forward
    else if(right > -1) {
      MFloat ptRight = a_variable(right, PV->RHO) * F1B3
                       + F1B2 / a_variable(right, PV->RHO)
                             * (a_variable(right, PV->U) * a_variable(right, PV->U)
                                + a_variable(right, PV->V) * a_variable(right, PV->V)
                                + a_variable(right, PV->W) * a_variable(right, PV->W));

      gradient = ptRight - pt;
    }
  }
  return gradient;
}

template <MInt nDim>
void LbSolver<nDim>::saveCoarseSolution() {
  TRACE();

  if(g_multiSolverGrid) {
    TERMM(1, "Does not work with multi-solver grid concept");
  }

  /*! \page propertyPage1
    \section pp_reductionLevel
    <code>MString* PostProcessingSolver::m_reductionLevel</code>\n
    default = <code>maxLevel()</code>\n\n
    This property determines the level to reduce the grid and the solution to.
    <ul>
    <li><code>level</code> </li>
    </ul>\n
    Keywords: <i>POSTPROCESSING</i>
  */
  MInt reductionLevel = Context::getSolverProperty<MInt>("pp_reductionLevel", m_solverId, AT_);

  // 1. Reduce the grid
  // In C++11 prefer lambdas over bind
#if defined(MAIA_MS_COMPILER)
  MInt recalcIdTreeSize =
      grid().raw().reduceToLevel(reductionLevel, [&](const MInt cellId) { this->interpolateToParentCells(cellId); });
#else
  // In C++03: prefer boost::bind, std::bind1st is deprecated.
  MInt recalcIdTreeSize =
      grid().raw().reduceToLevel(reductionLevel, std::bind(&LbSolver<nDim>::interpolateToParentCells, this, _1));
#endif

  MChar buf1[16];
  MChar buf2[16];
  MString preName = "";
  MString preName2 = "";
  MString preName3 = "";
  sprintf(buf1, "%d", (globalTimeStep - 1));
  sprintf(buf2, "%d", reductionLevel);

  preName.append("Lvl");
  preName.append(buf2);
  preName.append("_");
  preName2 = preName;
  preName2.append(buf1);
  preName2.append("_");
  preName3 = preName + "PV_";
  preName3.append(buf1);
  preName3 = preName3 + ParallelIo::fileExt();

  MString tmpFileNameGrid = preName2 + grid().gridInputFileName();
  MString tmpFileNameGridPath = outputDir() + tmpFileNameGrid;
  MString tmpOutputFileName = outputDir() + preName3;

  // This contains the recalculated Ids sorted by the tree after the save-call
  MIntScratchSpace recalcIdTree(recalcIdTreeSize, AT_, "recalcIdTree");
  grid().raw().saveGrid(tmpFileNameGridPath.c_str(), grid().raw().m_haloCells, grid().raw().m_windowCells,
                        grid().raw().m_azimuthalHaloCells, grid().raw().m_azimuthalUnmappedHaloCells,
                        grid().raw().m_azimuthalWindowCells, recalcIdTree.begin());

  // 3. Call function to write interpolated solution
  saveUVWRhoTOnlyPar(tmpOutputFileName.c_str(), tmpFileNameGrid.c_str(), recalcIdTree.begin());

  // 4. Exit MAIA
  TERMM(0, AT_);
}

/** \brief Interpolates variables from one level to its parent level
 *
 * \author Thomas Schilden
 * \date 31.08.2017
 *
 * \param[in] cellId
 * \param[out] the vars
 *
 **/
template <MInt nDim>
void LbSolver<nDim>::accessSampleVariables(MInt cellId, MFloat*& vars) {
  TRACE();
  vars = &(a_variable(cellId, 0));
}

/** \brief Interpolates variables from one level to its parent level
 *
 * \author Andreas Lintermann
 * \date 26.08.2012
 *
 * The interpolation is done by averaging over the children of the parent and considering the volume change.
 *
 * \param[in] parentlevel the level to interpolate to, must be smaller than maxLevel()
 *
 **/
template <MInt nDim>
void LbSolver<nDim>::interpolateToParentCells(MInt parentlevel) {
  TRACE();

  if(parentlevel == grid().maxLevel()) {
    TERMM(1, "Interpolation to parent cells residing on the finest level makes no sense, exiting ... ");
  }

  for(MInt id = 0; id < noInternalCells(); id++) {
    // skip cells not on parentlevel
    if(a_level(id) != parentlevel || c_isLeafCell(id)) continue;

    for(MInt v = 0; v < m_noVariables; v++) {
      a_variable(id, v) = .0;
    }

    MFloat avg_factor = 1.0 / c_noChildren(id);

    for(MInt ch = 0; ch < IPOW2(nDim); ch++) {
      // skip
      if(c_childId(id, ch) == -1) continue;

      // sum up
      for(MInt v = 0; v < m_noVariables; v++) {
        a_variable(id, v) += avg_factor * (a_variable(c_childId(id, ch), v));
      }
    }
  }
}

/** \brief This function loads the flow information of the cells
 *        such as variables and attributes like u_velocity,density,etc. for restart purpose.
 *
 * \author Andreas Lintermann, Konstantin Froehlich
 * \date 14.02.2012
 * In contrast to loadGridFlowVariables(const MChar* fileName) this function loads the
 * information for restart in parallel only for the primitive variables.
 * The attribute 'name' of the variables are set to their according meaning.
 *
 * \param[in] fileName the name of the file to load the data from
 *
 */
template <MInt nDim>
void LbSolver<nDim>::loadRestartWithoutDistributionsPar(const MChar* fileName) {
  TRACE();
  if constexpr(nDim != 2 && nDim != 3) {
    cerr << " In global function loadRestartWithoutDistributionsPar: wrong number of dimensions !" << endl;
    exit(0);
    return;
  }

  // File loading.
  using namespace maia::parallel_io;
  ParallelIo parallelIo(fileName, PIO_READ, mpiComm());

  // read time
  parallelIo.getAttribute(&m_time, "time");
  parallelIo.getAttribute(&globalTimeStep, "globalTimeStep");

  if(m_velocityControl.restart) {
    parallelIo.getAttribute(&m_velocityControl.lastGlobalAvgV, "lastGlobalAvgV");
    parallelIo.getAttribute(&m_velocityControl.previousError, "previousError");
    parallelIo.getAttribute(&m_velocityControl.integratedError, "integratedError");
    parallelIo.getAttribute(&m_velocityControl.derivedError, "derivedError");
    for(MInt i = 0; i < nDim; i++) {
      parallelIo.getAttribute(&m_volumeAccel[i], "volumeAcceleration_" + std::to_string(i));
    }
  }

  // This should be the same for all the variables
  MPI_Offset dimLen = noInternalCells();
  MPI_Offset start = domainOffset(domainId());
  parallelIo.setOffset(dimLen, start);
  // Load our variables

  { // Velocity u
    MFloatScratchSpace tmp_velocityU((MInt)dimLen, AT_, "tmp_velocityU");
    parallelIo.readArray(tmp_velocityU.getPointer(), "variables0");
    for(MInt i = 0; i < (MInt)dimLen; ++i)
      a_variable(i, PV->U) = tmp_velocityU.p[i];
  }

  { // Velocity v
    MFloatScratchSpace tmp_velocityV((MInt)dimLen, AT_, "tmp_velocityV");
    parallelIo.readArray(tmp_velocityV.getPointer(), "variables1");
    for(MInt i = 0; i < (MInt)dimLen; ++i)
      a_variable(i, PV->V) = tmp_velocityV.p[i];
  }

  if constexpr(nDim == 3) {
    { // Velocity w
      MFloatScratchSpace tmp_velocityW((MInt)dimLen, AT_, "tmp_velocityW");
      parallelIo.readArray(tmp_velocityW.getPointer(), "variables2");
      for(MInt i = 0; i < (MInt)dimLen; ++i)
        a_variable(i, PV->W) = tmp_velocityW.p[i];
    }

    { // Density rho
      MFloatScratchSpace tmp_rho((MInt)dimLen, AT_, "tmp_rho");
      parallelIo.readArray(tmp_rho.getPointer(), "variables3");
      for(MInt i = 0; i < (MInt)dimLen; ++i)
        a_variable(i, PV->RHO) = tmp_rho.p[i];
    }

    if(m_isThermal && !m_isTransport) {
      // Temperature t
      MFloatScratchSpace tmp_t((MInt)dimLen, AT_, "tmp_t");
      parallelIo.readArray(tmp_t.getPointer(), "variables4");
      for(MInt i = 0; i < (MInt)dimLen; ++i)
        a_variable(i, PV->T) = tmp_t.p[i];
    }
    if(m_isTransport && !m_isThermal) {
      // Concentration c
      MFloatScratchSpace tmp_c((MInt)dimLen, AT_, "tmp_c");
      parallelIo.readArray(tmp_c.getPointer(), "variables4");
      for(MInt i = 0; i < (MInt)dimLen; ++i)
        a_variable(i, PV->C) = tmp_c.p[i];
    }
    if(m_isTransport && m_isThermal) {
      // Temperature t
      MFloatScratchSpace tmp_t((MInt)dimLen, AT_, "tmp_t");
      parallelIo.readArray(tmp_t.getPointer(), "variables4");
      for(MInt i = 0; i < (MInt)dimLen; ++i)
        a_variable(i, PV->T) = tmp_t.p[i];
      // Concentration c
      MFloatScratchSpace tmp_c((MInt)dimLen, AT_, "tmp_c");
      parallelIo.readArray(tmp_c.getPointer(), "variables5");
      for(MInt i = 0; i < (MInt)dimLen; ++i)
        a_variable(i, PV->C) = tmp_c.p[i];
    }
    if(m_isEELiquid) {
      MFloatScratchSpace tmp_alphaG((MInt)dimLen, AT_, "tmp_alphaG");
      parallelIo.readArray(tmp_alphaG.getPointer(), "variables4");
      for(MInt i = 0; i < (MInt)dimLen; ++i)
        a_alphaGas(i) = tmp_alphaG.p[i];
    }
  } else {
    { // Density rho
      MFloatScratchSpace tmp_rho((MInt)dimLen, AT_, "tmp_rho");
      parallelIo.readArray(tmp_rho.getPointer(), "variables2");
      for(MInt i = 0; i < (MInt)dimLen; ++i)
        a_variable(i, PV->RHO) = tmp_rho.p[i];
    }

    if(m_isThermal && !m_isTransport) {
      // Temperature t
      MFloatScratchSpace tmp_t((MInt)dimLen, AT_, "tmp_t");
      parallelIo.readArray(tmp_t.getPointer(), "variables3");
      for(MInt i = 0; i < (MInt)dimLen; ++i)
        a_variable(i, PV->T) = tmp_t.p[i];
    }
    if(m_isTransport && !m_isThermal) {
      // Concentration c
      MFloatScratchSpace tmp_c((MInt)dimLen, AT_, "tmp_c");
      parallelIo.readArray(tmp_c.getPointer(), "variables3");
      for(MInt i = 0; i < (MInt)dimLen; ++i)
        a_variable(i, PV->C) = tmp_c.p[i];
    }
    if(m_isTransport && m_isThermal) {
      // Temperature t
      MFloatScratchSpace tmp_t((MInt)dimLen, AT_, "tmp_t");
      parallelIo.readArray(tmp_t.getPointer(), "variables3");
      for(MInt i = 0; i < (MInt)dimLen; ++i)
        a_variable(i, PV->T) = tmp_t.p[i];
      // Concentration c
      MFloatScratchSpace tmp_c((MInt)dimLen, AT_, "tmp_c");
      parallelIo.readArray(tmp_c.getPointer(), "variables4");
      for(MInt i = 0; i < (MInt)dimLen; ++i)
        a_variable(i, PV->C) = tmp_c.p[i];
    }
  }

  { // old Velocity u
    MFloatScratchSpace tmp_oldVelocityU((MInt)dimLen, AT_, "tmp_oldVelocityU");
    parallelIo.readArray(tmp_oldVelocityU.getPointer(), "oldVariables0");
    for(MInt i = 0; i < (MInt)dimLen; ++i) {
      a_oldVariable(i, PV->U) = tmp_oldVelocityU.p[i];
    }
  }

  { // old Velocity v
    MFloatScratchSpace tmp_oldVelocityV((MInt)dimLen, AT_, "tmp_oldVelocityV");
    parallelIo.readArray(tmp_oldVelocityV.getPointer(), "oldVariables1");
    for(MInt i = 0; i < (MInt)dimLen; ++i)
      a_oldVariable(i, PV->V) = tmp_oldVelocityV.p[i];
  }

  if constexpr(nDim == 3) {
    { // old Velocity w
      MFloatScratchSpace tmp_oldVelocityW((MInt)dimLen, AT_, "tmp_oldVelocityW");
      parallelIo.readArray(tmp_oldVelocityW.getPointer(), "oldVariables2");
      for(MInt i = 0; i < (MInt)dimLen; ++i)
        a_oldVariable(i, PV->W) = tmp_oldVelocityW.p[i];
    }

    { // old Density rho
      MFloatScratchSpace tmp_oldRho((MInt)dimLen, AT_, "tmp_oldRho");
      parallelIo.readArray(tmp_oldRho.getPointer(), "oldVariables3");
      for(MInt i = 0; i < (MInt)dimLen; ++i)
        a_oldVariable(i, PV->RHO) = tmp_oldRho.p[i];
    }

    if(m_isThermal && !m_isTransport) {
      { // old Temperature t
        MFloatScratchSpace tmp_oldT((MInt)dimLen, AT_, "tmp_oldT");
        parallelIo.readArray(tmp_oldT.getPointer(), "oldVariables4");
        for(MInt i = 0; i < (MInt)dimLen; ++i)
          a_oldVariable(i, PV->T) = tmp_oldT.p[i];
      }
    }
    if(m_isTransport && !m_isThermal) {
      { // old Concentration c
        MFloatScratchSpace tmp_oldC((MInt)dimLen, AT_, "tmp_oldC");
        parallelIo.readArray(tmp_oldC.getPointer(), "oldVariables4");
        for(MInt i = 0; i < (MInt)dimLen; ++i)
          a_oldVariable(i, PV->C) = tmp_oldC.p[i];
      }
    }
    if(m_isTransport && m_isThermal) {
      { // old Temperature t
        MFloatScratchSpace tmp_oldT((MInt)dimLen, AT_, "tmp_oldT");
        parallelIo.readArray(tmp_oldT.getPointer(), "oldVariables4");
        for(MInt i = 0; i < (MInt)dimLen; ++i)
          a_oldVariable(i, PV->T) = tmp_oldT.p[i];
      }
      { // old Concentration c
        MFloatScratchSpace tmp_oldC((MInt)dimLen, AT_, "tmp_oldC");
        parallelIo.readArray(tmp_oldC.getPointer(), "oldVariables5");
        for(MInt i = 0; i < (MInt)dimLen; ++i)
          a_oldVariable(i, PV->C) = tmp_oldC.p[i];
      }
    }
  } else {
    { // old Density rho
      MFloatScratchSpace tmp_oldRho((MInt)dimLen, AT_, "tmp_oldRho");
      parallelIo.readArray(tmp_oldRho.getPointer(), "oldVariables2");
      for(MInt i = 0; i < (MInt)dimLen; ++i)
        a_oldVariable(i, PV->RHO) = tmp_oldRho.p[i];
    }

    if(m_isThermal && !m_isTransport) {
      // old Temperature t
      MFloatScratchSpace tmp_oldT((MInt)dimLen, AT_, "tmp_oldT");
      parallelIo.readArray(tmp_oldT.getPointer(), "oldVariables3");
      for(MInt i = 0; i < (MInt)dimLen; ++i)
        a_oldVariable(i, PV->T) = tmp_oldT.p[i];
    }
    if(m_isTransport && !m_isThermal) {
      // old Concentration c
      MFloatScratchSpace tmp_oldC((MInt)dimLen, AT_, "tmp_oldC");
      parallelIo.readArray(tmp_oldC.getPointer(), "oldVariables3");
      for(MInt i = 0; i < (MInt)dimLen; ++i)
        a_oldVariable(i, PV->C) = tmp_oldC.p[i];
    }
    if(m_isTransport && m_isThermal) {
      // old Temperature t
      MFloatScratchSpace tmp_oldT((MInt)dimLen, AT_, "tmp_oldT");
      parallelIo.readArray(tmp_oldT.getPointer(), "oldVariables3");
      for(MInt i = 0; i < (MInt)dimLen; ++i)
        a_oldVariable(i, PV->T) = tmp_oldT.p[i];
      // old Concentration c
      MFloatScratchSpace tmp_oldC((MInt)dimLen, AT_, "tmp_oldC");
      parallelIo.readArray(tmp_oldC.getPointer(), "oldVariables4");
      for(MInt i = 0; i < (MInt)dimLen; ++i)
        a_oldVariable(i, PV->C) = tmp_oldC.p[i];
    }
  }

  m_bndCnd->bcDataReadRestartData(parallelIo);
}

/** \brief
 *
 *  Sets the neighbor information for all directions including diagonals.
 *  Uses the information stored in m_neighborList.
 *
 */
template <MInt nDim>
void LbSolver<nDim>::determineEqualDiagonals() {
  TRACE();
  // do nothing, all diagonals are created in createMultiSolverInformationPar(), see cartesiangrid.h
}

/** \brief writes the flow variables and distributions to disc
 *
 *
 */
template <MInt nDim>
void LbSolver<nDim>::saveRestartFile() {
  TRACE();

  stringstream fileName;
  stringstream grid_;

  fileName << outputDir() << "restart_" << getIdentifier() << globalTimeStep << ParallelIo::fileExt();
  grid_ << grid().gridInputFileName();

  saveRestartWithDistributionsPar(fileName.str().c_str(), grid().gridInputFileName().c_str());
}

/** \brief loads the flow variables and distributions from disc
 *
 *
 */
template <MInt nDim>
void LbSolver<nDim>::loadRestartFile() {
  TRACE();

  NEW_TIMER_GROUP(t_restart, "Restart");
  NEW_TIMER(t_restartAll, "restart file loading", t_restart);
  RECORD_TIMER_START(t_restartAll);

  m_log << "#########################################################################################################"
           "#############"
        << endl;
  m_log << "##                                              Loading restart file                                     "
           "           ##"
        << endl;
  m_log << "#########################################################################################################"
           "#############"
        << endl
        << endl;


  MBool initializeRestart = false;
  initializeRestart = Context::getSolverProperty<MBool>("initRestart", m_solverId, AT_, &initializeRestart);
  stringstream varFileName;

  m_log << "  + Restart file information:" << endl;
  m_log << "    - restart time step: " << m_restartTimeStep << endl;
  if(m_initFromCoarse)
    m_log << "    - type:              "
          << "from coarse mac. variables " << endl;
  else
    m_log << "    - type:              " << (initializeRestart ? "from mac. variables" : "from distributions") << endl;

  varFileName << restartDir() << "restart_" << getIdentifier() << m_restartTimeStep << ParallelIo::fileExt();

  m_log << "    - file name:         " << varFileName.str() << endl << endl;
  m_log << "  + Loading the file ..." << endl;

  if(m_initFromCoarse)
    loadRestartWithoutDistributionsParFromCoarse((varFileName.str()).c_str());
  else {
    // restart with distributions
    if(!initializeRestart) loadRestartWithDistributionsPar((varFileName.str()).c_str());
    // restart only from macroscopic variables
    else
      loadRestartWithoutDistributionsPar((varFileName.str()).c_str());
  }
  m_log << endl;

  if(m_velocityControl.restart) initVolumeForces();

  RECORD_TIMER_STOP(t_restartAll);
}

//----------------------------------------------------------------------------

/** \brief Exchanges all required data between neighboring processes with non-blocking communication.
 *
 * \author Andreas Lintermann
 * \date 08.04.2013
 *
 * Algorithm is as follows:
 *
 *  1. Wait for unfinished sending operations from last iteration.
 *  2. Gather all information from the cells in a send buffer and send per neighbor domain.
 *  3. Wait for unfinished receiving operations from last iteration.
 *  4. Scatter all data from the receive buffer to the cells per finished neighbor domain.
 *  5. Open new receive.
 *
 **/
template <MInt nDim>
void LbSolver<nDim>::exchangeLbNB(MInt mode) {
  // 1. Wait for unfinished sending operations from last iteration.
  if(noNeighborDomains() <= 0) return;

  if(mode == 0) {
    // TODO labels:LB,COMM As soon as the steps of the collision step are separated
    // the implementation of the non-blocking communication can go on.
    // At the moment the communication cannot be hidden behind anything.
    // 0. Calculate the macroscopic variables for all cells first.
    // They are needed for the boundary conditons
    //    RECORD_TIMER_START(m_t.collision);
    //    for(MInt i = 0; i < m_currentMaxNoCells; i++) {
    //      MInt pCellId = m_activeCellList[i];
    //      MFloat rho = 0;
    //      Mfloat u[nDim] = {F0};
    //      this->calculateMacroscopicVariables(pCellId, &rho, u);
    //    a_variables(pCellId, PV->RHO) = rho;
    //    for(MInt d = 0; d < nDim; d++) {
    //      a_variables(pCellId, d) = u[d];
    //    RECORD_TIMER_STOP(m_t.collision);

    // 1. Start the receive of data from neighbors.
    RECORD_TIMER_START(m_t.communication);
    (this->*m_receiveMethod)();
    RECORD_TIMER_STOP(m_t.communication);
    // 2. Perform the collision of the window cells
    //    RECORD_TIMER_START(m_t.collision);
    //    (this->*m_solutionStepMethod)(m_activeWinCellList, m_currentMaxNoWinCells);
    //    RECORD_TIMER_STOP(m_t.collision);

    // 3. Gather all information from the cells in a send buffer.
    RECORD_TIMER_START(m_t.packing);
    (this->*m_gatherMethod)();
    RECORD_TIMER_STOP(m_t.packing);

    // 4. Send the gathered information. In the meantime the collision for all cells except
    // window cells is done
    RECORD_TIMER_START(m_t.communication);
    (this->*m_sendMethod)();
    RECORD_TIMER_STOP(m_t.communication);

  } else {
    RECORD_TIMER_START(m_t.communication);
    MPI_Waitall(noNeighborDomains(), &mpi_requestR[0], MPI_STATUS_IGNORE, AT_);
    RECORD_TIMER_STOP(m_t.communication);

    // 5. Scatter all data from the receive buffer to the cells.
    RECORD_TIMER_START(m_t.unpacking);
    (this->*m_scatterMethod)();
    RECORD_TIMER_STOP(m_t.unpacking);

    RECORD_TIMER_START(m_t.communication);
    MPI_Waitall(noNeighborDomains(), &mpi_requestS[0], MPI_STATUS_IGNORE, AT_);
    RECORD_TIMER_STOP(m_t.communication);
  }
}

/** \brief Exchanges all required data between neighboring processes with blocking communication.
 *
 * \author Andreas Lintermann
 * \date 05.04.2013
 *
 * Algorithm is as follows:
 *
 *  1. Gather all information from the cells in a send buffer.
 *  2. Send the gathered information.
 *  3. Receive data from neighbors.
 *  4. Scatter all data from the receive buffer to the cells.
 *
 **/
template <MInt nDim>
void LbSolver<nDim>::exchangeLb(MInt mode) {
  TRACE();

  if(mode == 0) return; // only executed for non-blocking

  // 1. Gather all information from the cells in a send buffer.
  RECORD_TIMER_START(m_t.packing);
  (this->*m_gatherMethod)();
  RECORD_TIMER_STOP(m_t.packing);

  RECORD_TIMER_START(m_t.communication);
  // 2. Send the gathered information.
  (this->*m_sendMethod)();
  // 3. Receive data from neighbors.
  (this->*m_receiveMethod)();
  RECORD_TIMER_STOP(m_t.communication);

  // 4. Scatter all data from the receive buffer to the cells.
  RECORD_TIMER_START(m_t.unpacking);
  (this->*m_scatterMethod)();
  RECORD_TIMER_STOP(m_t.unpacking);
}

template <MInt nDim>
void LbSolver<nDim>::exchange(MInt mode) {
  //  TRACE();
  RECORD_TIMER_START(m_t.solutionStep);
  RECORD_TIMER_START(m_t.exchange);
  (this->*m_exchangeMethod)(mode);
  RECORD_TIMER_STOP(m_t.exchange);
  RECORD_TIMER_STOP(m_t.solutionStep);
}


//----------------------------------------------------------------------------


/** \brief Collects all necessary data from the cells into the send buffer.
 *
 * \author Andreas Lintermann
 * \date 05.04.2013
 *
 * This is a highly efficient gather function which collect data in the send buffer
 * as previously determined by the method.
 *
 * \param[in] var the id of the data element to be transferred
 *
 **/
template <MInt nDim>
inline void LbSolver<nDim>::gatherNormal() {
  TRACE();

  for(MInt n = 0; n < noNeighborDomains(); n++) {
    for(MInt var = 0; var < m_noElementsTransfer; var++) {
      MInt sendBufferCounter = m_nghbrOffsetsWindow[n][var];
      if(var == 0) {
        maia::parallelFor(0, noWindowCells(n), [=](MInt j) {
          for(MInt d = 0; d < m_noDistributions; d++) {
            MInt id = windowCell(n, j);
            MInt offset = sendBufferCounter + (j * m_noDistributions + d);
            m_sendBuffers[n][offset] = a_distribution(id, d);
          }
        });
      } else if(var == 1 && m_isThermal) {
        maia::parallelFor(0, noWindowCells(n), [=](MInt j) {
          for(MInt d = 0; d < m_noDistributions; d++) {
            MInt id = windowCell(n, j);
            MInt offset = sendBufferCounter + (j * m_noDistributions + d);
            m_sendBuffers[n][offset] = a_distributionThermal(id, d);
          }
        });
      } else if(var == 1 && m_isTransport) {
        maia::parallelFor(0, noWindowCells(n), [=](MInt j) {
          for(MInt d = 0; d < m_noDistributions; d++) {
            MInt id = windowCell(n, j);
            MInt offset = sendBufferCounter + (j * m_noDistributions + d);
            m_sendBuffers[n][offset] = a_distributionTransport(id, d);
          }
        });
      } else if(var == 2 && m_isThermal && m_isTransport) {
        maia::parallelFor(0, noWindowCells(n), [=](MInt j) {
          for(MInt d = 0; d < m_noDistributions; d++) {
            MInt id = windowCell(n, j);
            MInt offset = sendBufferCounter + (j * m_noDistributions + d);
            m_sendBuffers[n][offset] = a_distributionTransport(id, d);
          }
        });
      } else if(var == 1 || (var == 2 && (m_isThermal || m_isTransport))
                || (var == 3 && m_isThermal && m_isTransport)) {
        maia::parallelFor(0, noWindowCells(n), [=](MInt j) {
          for(MInt v = 0; v < m_noVariables; v++) {
            MInt id = windowCell(n, j);
            MInt offset = sendBufferCounter + (j * m_noVariables + v);
            m_sendBuffers[n][offset] = a_variable(id, v);
          }
        });
      } else {
        maia::parallelFor(0, noWindowCells(n), [=](MInt j) {
          for(MInt v = 0; v < m_noVariables; v++) {
            MInt id = windowCell(n, j);
            MInt offset = sendBufferCounter + (j * m_noVariables + v);
            m_sendBuffers[n][offset] = a_oldVariable(id, v);
          }
        });
      }
      //      for(MInt j = 0; j < noWindowCells(n); j++) {
      //        memcpy((void*)&m_sendBuffers[n][sendBufferCounter],
      //               m_baseAddresses[var] + (m_dataBlockSizes[var] * windowCell(n, j)),
      //               m_dataBlockSizes[var] * sizeof(MFloat));
      //        sendBufferCounter += m_dataBlockSizes[var];
      //      }
    }
  }
}

/** \brief Collects all necessary data from the cells into the send buffer.
 *
 * \author Andreas Lintermann
 * \date 05.04.2013
 *
 * This is a highly efficient gather function which collect data in the send buffer
 * as previously determined by the method.
 *
 * \param[in] var the id of the data element to be transferred
 *
 **/
template <MInt nDim>
inline void LbSolver<nDim>::gatherReduced() {
  TRACE();

  for(MInt n = 0; n < noNeighborDomains(); n++) {
    for(MInt var = 0; var < m_noElementsTransfer; var++) {
      if(var == 0) {
        maia::parallelFor<true>(0, noWindowCells(n), [=](MInt j) {
          for(MInt d = 0; d < (m_noDistributions - 1); d++) {
            if(m_windowDistsForExchange[n][j * (m_noDistributions - 1) + d] > -1) {
              MInt offset = m_windowDistsForExchange[n][j * (m_noDistributions - 1) + d];
              m_sendBuffers[n][offset] = a_distribution(windowCell(n, j), d);
            }
          }
        });
      } else if(var == 1 && m_isThermal) {
        maia::parallelFor<true>(0, noWindowCells(n), [=](MInt j) {
          for(MInt d = 0; d < (m_noDistributions - 1); d++) {
            if(m_windowDistsForExchange[n][j * (m_noDistributions - 1) + d] > -1) {
              MInt offset =
                  m_windowDistsForExchange[n][j * (m_noDistributions - 1) + d] + m_noWindowDistDataPerDomain[n];
              m_sendBuffers[n][offset] = a_distributionThermal(windowCell(n, j), d);
            }
          }
        });
      } else if(var == 1 && m_isTransport) {
        maia::parallelFor<true>(0, noWindowCells(n), [=](MInt j) {
          for(MInt d = 0; d < (m_noDistributions - 1); d++) {
            if(m_windowDistsForExchange[n][j * (m_noDistributions - 1) + d] > -1) {
              MInt offset =
                  m_windowDistsForExchange[n][j * (m_noDistributions - 1) + d] + m_noWindowDistDataPerDomain[n];
              m_sendBuffers[n][offset] = a_distributionTransport(windowCell(n, j), d);
            }
          }
        });
      } else if(var == 2 && m_isThermal && m_isTransport) {
        maia::parallelFor<true>(0, noWindowCells(n), [=](MInt j) {
          for(MInt d = 0; d < (m_noDistributions - 1); d++) {
            if(m_windowDistsForExchange[n][j * (m_noDistributions - 1) + d] > -1) {
              MInt offset =
                  m_windowDistsForExchange[n][j * (m_noDistributions - 1) + d] + 2 * m_noWindowDistDataPerDomain[n];
              m_sendBuffers[n][offset] = a_distributionTransport(windowCell(n, j), d);
            }
          }
        });
      } else if(var == 1 || (var == 2 && (m_isThermal || m_isTransport))
                || (var == 3 && m_isThermal && m_isTransport)) {
        maia::parallelFor<true>(0, noWindowCells(n), [=](MInt j) {
          for(MInt v = 0; v < m_noVariables; v++) {
            MInt id = windowCell(n, j);
            MInt offset = (j * m_noVariables + v) + (m_noDistsTransfer * m_noWindowDistDataPerDomain[n]);
            m_sendBuffers[n][offset] = a_variable(id, v);
          }
        });
      } else {
        maia::parallelFor<true>(0, noWindowCells(n), [=](MInt j) {
          for(MInt v = 0; v < m_noVariables; v++) {
            MInt id = windowCell(n, j);
            MInt offset = (j * m_noVariables + v) + (m_noDistsTransfer * m_noWindowDistDataPerDomain[n])
                          + (noWindowCells(n) * m_noVariables);
            m_sendBuffers[n][offset] = a_oldVariable(id, v);
          }
        });
      }
    }
    //    MInt totalData =
    //        m_noDistsTransfer * m_noWindowDistDataPerDomain[n] + noWindowCells(n) * m_noVariables * m_noVarsTransfer;
    //    for(MInt j = 0; j < totalData; j++) {
    //      memcpy((void*)&m_sendBuffers[n][j], m_commPtWindow[n][j], sizeof(MFloat));
    //    }
  }
}

/** \brief Scatteres all necessary data in from the receive buffer to the cells.
 *
 * \author Andreas Lintermann
 * \date 05.04.2013
 *
 * This is a highly efficient scatter function which distributes data from the receive buffer
 * as previously determined by the method.
 *
 * \param[in] var the id of the data element to be scattered
 *
 **/
template <MInt nDim>
inline void LbSolver<nDim>::scatterNormal() {
  TRACE();

  for(MInt n = 0; n < noNeighborDomains(); n++) {
    for(MInt var = 0; var < m_noElementsTransfer; var++) {
      MInt receiveBufferCounter = m_nghbrOffsetsHalo[n][var];
      if(var == 0) {
        maia::parallelFor(0, noHaloCells(n), [=](MInt j) {
          for(MInt d = 0; d < m_noDistributions; d++) {
            MInt id = haloCell(n, j);
            MInt offset = receiveBufferCounter + (j * m_noDistributions + d);
            a_distribution(id, d) = m_receiveBuffers[n][offset];
          }
        });
      } else if(var == 1 && m_isThermal) {
        maia::parallelFor(0, noHaloCells(n), [=](MInt j) {
          for(MInt d = 0; d < m_noDistributions; d++) {
            MInt id = haloCell(n, j);
            MInt offset = receiveBufferCounter + (j * m_noDistributions + d);
            a_distributionThermal(id, d) = m_receiveBuffers[n][offset];
          }
        });
      } else if(var == 1 && m_isTransport) {
        maia::parallelFor(0, noHaloCells(n), [=](MInt j) {
          for(MInt d = 0; d < m_noDistributions; d++) {
            MInt id = haloCell(n, j);
            MInt offset = receiveBufferCounter + (j * m_noDistributions + d);
            a_distributionTransport(id, d) = m_receiveBuffers[n][offset];
          }
        });
      } else if(var == 2 && m_isThermal && m_isTransport) {
        maia::parallelFor(0, noHaloCells(n), [=](MInt j) {
          for(MInt d = 0; d < m_noDistributions; d++) {
            MInt id = haloCell(n, j);
            MInt offset = receiveBufferCounter + (j * m_noDistributions + d);
            a_distributionTransport(id, d) = m_receiveBuffers[n][offset];
          }
        });
      } else if(var == 1 || (var == 2 && (m_isThermal || m_isTransport))
                || (var == 3 && m_isThermal && m_isTransport)) {
        maia::parallelFor(0, noHaloCells(n), [=](MInt j) {
          for(MInt v = 0; v < m_noVariables; v++) {
            MInt id = haloCell(n, j);
            MInt offset = receiveBufferCounter + (j * m_noVariables + v);
            a_variable(id, v) = m_receiveBuffers[n][offset];
          }
        });
      } else {
        maia::parallelFor(0, noHaloCells(n), [=](MInt j) {
          for(MInt v = 0; v < m_noVariables; v++) {
            MInt id = haloCell(n, j);
            MInt offset = receiveBufferCounter + (j * m_noVariables + v);
            a_oldVariable(id, v) = m_receiveBuffers[n][offset];
          }
        });
      }
      //      for(MInt j = 0; j < noHaloCells(n); j++) {
      //        memcpy(m_baseAddresses[var] + (m_dataBlockSizes[var] * haloCell(n, j)),
      //               (void*)&m_receiveBuffers[n][receiveBufferCounter],
      //               m_dataBlockSizes[var] * sizeof(MFloat));
      //        receiveBufferCounter += m_dataBlockSizes[var];
      //      }
    }
  }
}

/** \brief Scatteres all necessary data in from the receive buffer to the cells.
 *
 * \author Andreas Lintermann
 * \date 30.07.2015
 *
 * This is a highly efficient scatter function which distributes data from the receive buffer
 * as previously determined by the method.
 *
 * \param[in] var the id of the data element to be scattered
 *
 **/
template <MInt nDim>
inline void LbSolver<nDim>::scatterReduced() {
  TRACE();

  for(MInt n = 0; n < noNeighborDomains(); n++) {
    for(MInt var = 0; var < m_noElementsTransfer; var++) {
      if(var == 0) {
        maia::parallelFor<true>(0, noHaloCells(n), [=](MInt j) {
          for(MInt d = 0; d < (m_noDistributions - 1); d++) {
            if(m_haloDistsForExchange[n][j * (m_noDistributions - 1) + d] > -1) {
              MInt offset = m_haloDistsForExchange[n][j * (m_noDistributions - 1) + d];
              a_distribution(haloCell(n, j), d) = m_receiveBuffers[n][offset];
            }
          }
        });
      } else if(var == 1 && m_isThermal) {
        maia::parallelFor<true>(0, noHaloCells(n), [=](MInt j) {
          for(MInt d = 0; d < (m_noDistributions - 1); d++) {
            if(m_haloDistsForExchange[n][j * (m_noDistributions - 1) + d] > -1) {
              MInt offset = m_haloDistsForExchange[n][j * (m_noDistributions - 1) + d] + m_noHaloDistDataPerDomain[n];
              a_distributionThermal(haloCell(n, j), d) = m_receiveBuffers[n][offset];
            }
          }
        });
      } else if(var == 1 && m_isTransport) {
        maia::parallelFor<true>(0, noHaloCells(n), [=](MInt j) {
          for(MInt d = 0; d < (m_noDistributions - 1); d++) {
            if(m_haloDistsForExchange[n][j * (m_noDistributions - 1) + d] > -1) {
              MInt offset = m_haloDistsForExchange[n][j * (m_noDistributions - 1) + d] + m_noHaloDistDataPerDomain[n];
              a_distributionTransport(haloCell(n, j), d) = m_receiveBuffers[n][offset];
            }
          }
        });
      } else if(var == 2 && m_isThermal && m_isTransport) {
        maia::parallelFor<true>(0, noHaloCells(n), [=](MInt j) {
          for(MInt d = 0; d < (m_noDistributions - 1); d++) {
            if(m_haloDistsForExchange[n][j * (m_noDistributions - 1) + d] > -1) {
              MInt offset =
                  m_haloDistsForExchange[n][j * (m_noDistributions - 1) + d] + 2 * m_noHaloDistDataPerDomain[n];
              a_distributionTransport(haloCell(n, j), d) = m_receiveBuffers[n][offset];
            }
          }
        });
      } else if(var == 1 || (var == 2 && (m_isThermal || m_isTransport))
                || (var == 3 && m_isThermal && m_isTransport)) {
        maia::parallelFor<true>(0, noHaloCells(n), [=](MInt j) {
          for(MInt v = 0; v < m_noVariables; v++) {
            MInt id = haloCell(n, j);
            MInt offset = (j * m_noVariables + v) + m_noDistsTransfer * m_noHaloDistDataPerDomain[n];
            a_variable(id, v) = m_receiveBuffers[n][offset];
          }
        });
      } else {
        maia::parallelFor<true>(0, noHaloCells(n), [=](MInt j) {
          for(MInt v = 0; v < m_noVariables; v++) {
            MInt id = haloCell(n, j);
            MInt offset = (j * m_noVariables + v) + m_noDistsTransfer * m_noHaloDistDataPerDomain[n]
                          + (noHaloCells(n) * m_noVariables);
            a_oldVariable(id, v) = m_receiveBuffers[n][offset];
          }
        });
      }
    }
    //    MInt totalData =
    //        m_noDistsTransfer * m_noHaloDistDataPerDomain[n] + noHaloCells(n) * m_noVariables * m_noVarsTransfer;
    //    for(MInt j = 0; j < totalData; j++) {
    //      memcpy(m_commPtHalo[n][j], (void*)&m_receiveBuffers[n][j], sizeof(MFloat));
    //    }
  }
}

/** \brief Sends the send buffer to corresponding processes.
 *
 * \author Andreas Lintermann
 * \date 05.04.2013
 *
 * This method uses blocking communication.
 *
 **/
template <MInt nDim>
void LbSolver<nDim>::sendNormal() {
  MInt bufSize = 0;
  for(MInt i = 0; i < noNeighborDomains(); i++) {
    bufSize = noWindowCells(i) * m_dataBlockSizeTotal;
    MPI_Issend(m_sendBuffers[i], bufSize, MPI_DOUBLE, neighborDomain(i), 0, mpiComm(), &mpi_request[i], AT_,
               "m_sendBuffers[i]");
  }
}

/** \brief Sends the send buffer to corresponding processes.
 *
 * \author Andreas Lintermann
 * \date 30.07.2015
 *
 * This method uses blocking communication.
 *
 **/
template <MInt nDim>
void LbSolver<nDim>::sendReduced() {
  TRACE();

  for(MInt i = 0; i < noNeighborDomains(); i++) {
    MInt totalData =
        m_noDistsTransfer * m_noWindowDistDataPerDomain[i] + noWindowCells(i) * m_noVariables * m_noVarsTransfer;
    if(totalData == 0) continue;
    MPI_Issend(m_sendBuffers[i], totalData, MPI_DOUBLE, neighborDomain(i), 0, mpiComm(), &mpi_request[i], AT_,
               "m_sendBuffers[i]");
  }
}

/** \brief Receives the data from the corresponding processes.
 *
 * \author Andreas Lintermann
 * \date 05.04.2013
 *
 * This method uses blocking communication.
 *
 **/
template <MInt nDim>
void LbSolver<nDim>::receiveNormal() {
  TRACE();

  // Test if send has been performed
  MPI_Status status;

  for(MInt i = 0; i < noNeighborDomains(); i++) {
    const MInt bufSize = noHaloCells(i) * m_dataBlockSizeTotal;
    MPI_Recv(m_receiveBuffers[i], bufSize, MPI_DOUBLE, neighborDomain(i), 0, mpiComm(), &status, AT_,
             "m_receiveBuffers[i]");
  }
  for(MInt i = 0; i < noNeighborDomains(); i++)
    MPI_Wait(&mpi_request[i], &status, AT_);
}

/** \brief Receives the data from the corresponding processes.
 *
 * \author Andreas Lintermann
 * \date 30.07.2015
 *
 * This method uses blocking communication.
 *
 **/
template <MInt nDim>
void LbSolver<nDim>::receiveReduced() {
  TRACE();

  // Test if send has been performed
  MPI_Status status;

  for(MInt i = 0; i < noNeighborDomains(); i++) {
    MInt totalData =
        m_noDistsTransfer * m_noHaloDistDataPerDomain[i] + noHaloCells(i) * m_noVariables * m_noVarsTransfer;
    if(totalData == 0) continue;
    MPI_Recv(m_receiveBuffers[i], totalData, MPI_DOUBLE, neighborDomain(i), 0, mpiComm(), &status, AT_,
             "m_receiveBuffers[i]");
  }
  for(MInt i = 0; i < noNeighborDomains(); i++) {
    MInt totalData =
        m_noDistsTransfer * m_noHaloDistDataPerDomain[i] + noHaloCells(i) * m_noVariables * m_noVarsTransfer;
    if(totalData == 0) continue;
    MPI_Wait(&mpi_request[i], MPI_STATUS_IGNORE, AT_);
  }
}

/** brief update list of active cells
 *
 * \author Georg Eitel-Amor, Andreas Lintermann
 * \date unknown, 19.12.2017
 *
 * Runs over all cells and deactivates those that are refined. Two cases have to be considered:
 *
 *  1. the mesh is uniformly refined (only leaf cells stay active)
 *  2. the mesh is locally refined (leaf cells and interface parents stay active)
 *
 * In a final stage the active cell list is set.
 *
 **/
template <MInt nDim>
void LbSolver<nDim>::setActiveCellList() {
  TRACE();

  m_log << endl;
  m_log << "  + Setting active cell list" << endl;

  setIsActiveDefaultStates();
  fillActiveCellList();
}

/** brief activate inner cells
 *
 *
 */
template <MInt nDim>
void LbSolver<nDim>::activateInnerCells() {
  TRACE();

  // m_activeCellList = m_allInnerCells;
  m_activeCellList = nullptr;
  m_currentMaxNoCells = (m_cells.size()
                         //			 - m_overallNoWindowCells
                         //			 - m_overallNoHaloCells
  );
}

/** brief activate only the window cells
 *
 *
 */
template <MInt nDim>
void LbSolver<nDim>::activateWindowCells() {
  TRACE();

  // m_activeCellList = m_allWindowCells;
  m_activeCellList = nullptr;
  // m_currentMaxNoCells = m_overallNoWindowCells;
  m_currentMaxNoCells = 0;
}

/** brief activate all cells, but the halo cells
 *
 *  This function is useful for the calculation of the residual
 *  since it should not be evaluated for halo cells
 *  This function makes use of the fact, that all halo cells
 *  are stored at the end of the collector. Therefore there is
 *  no need for an specific cell list.
 */
template <MInt nDim>
void LbSolver<nDim>::activateAllButHaloCells() {
  TRACE();

  // m_currentMaxNoCells = m_cells.size() - m_overallNoHaloCells;
  m_currentMaxNoCells = m_cells.size();
}

template <MInt nDim>
void LbSolver<nDim>::treatInterfaceCells() {
  TRACE();

  // determine interface cells and interpolation information
  NEW_SUB_TIMER(t_ifaceCells, "interface cells", m_t.solver);
  RECORD_TIMER_START(t_ifaceCells);

  m_log << endl;
  m_log << "  + Interface cell generation" << endl;

  resetInterfaceCells();
  initializeInterfaceCells();
  buildInterfaceCells();
  setInterpolationNeighbors();
  setInterpolationCoefficients();

  RECORD_TIMER_STOP(t_ifaceCells);
}

template <MInt nDim>
void LbSolver<nDim>::correctInterfaceBcCells() {
  TRACE();
  if(m_isRefined && m_correctInterfaceBcCells) {
    setInterpolationNeighborsBC();
    setInterpolationCoefficientsBC();
  }
}

/** \brief Allocate memory for cells located at an interface
 *
 * children and parents are stored seperately
 * \author Georg Eitel-Amor
 */
template <MInt nDim>
void LbSolver<nDim>::initializeInterfaceCells() {
  TRACE();

  const MInt interfaceLevels = maxLevel() - minLevel();

  if(interfaceLevels == 0) return;

  /*! \page propertyPage1
    \section interfaceCellSize
    <code>MFloat LbSolver::m_interfaceCellSize</code>\n
    default = <code>0.5</code>\n\n
    Ratio of number of interface cells to overall number of cells.
    <ul>
    <li><code>0.0 - 1.0</code> (it makes no sense to use more than 1.0)</li>
    </ul>\n
    Keywords: <i>LATTICE_BOLTZMANN, NUMERICAL_SETUP, INITILIZATION</i>
  */
  m_interfaceCellSize = 0.5;
  m_interfaceCellSize = Context::getSolverProperty<MFloat>("interfaceCellSize", m_solverId, AT_, &m_interfaceCellSize);

  if(m_interfaceCellSize > noInternalCells())
    TERMM(1,
          " LbSolver::initializeInterfaceCells: You have requested more space for interface cells than for "
          "normal "
          "internal cells. What is the sense of that? Exiting!");

  MInt nointerfaceCellSize = (MInt)(m_interfaceCellSize * noInternalCells());

  m_log << "    - initializing interface cells for levels " << minLevel() << " to " << maxLevel() << endl;
  m_log << "      * currently using a size of " << nointerfaceCellSize << " for interface child and parent cells"
        << endl;
  m_log << "      * ratio: " << m_interfaceCellSize << endl;
  m_log << "      * collector size: "
        << (MFloat)(IPOW2(nDim) * sizeof(MFloat) + IPOW2(nDim) * sizeof(MInt)) * (MFloat)nointerfaceCellSize
               / (1024.0 * 1024.0)
        << " MB" << endl;


  // ----------- build collector for interface child cells
  for(MInt l = 0; l < interfaceLevels; l++) {
    m_interfaceChildren.emplace_back();
  }

  for(MInt i = 0; i < interfaceLevels; i++) {
    m_interfaceChildren[i] = new Collector<LbInterfaceCell>(nointerfaceCellSize, nDim, m_noDistributions);
  }

  // ----------- build collector for interface parent cells
  for(MInt l = 0; l < interfaceLevels; l++) {
    m_interfaceParents.emplace_back();
  }

  for(MInt i = 0; i < interfaceLevels; i++) {
    m_interfaceParents[i] = new Collector<LbParentCell>(nointerfaceCellSize, nDim, m_noDistributions);
  }
}

/** \brief Reset collectors for interface cells
 *\author Georg Eitel-Amor
 */
template <MInt nDim>
void LbSolver<nDim>::resetInterfaceCells() {
  TRACE();

  m_log << "    - resetting interface cells" << endl;
  const MInt interfaceLevels = maxLevel() - minLevel();

  if(interfaceLevels == 0) return;

  MInt noCells = m_cells.size();

  // reset interface marker for all cells
  //  for (MInt i = 0; i < noInternalCells(); i++){
  for(MInt i = 0; i < noCells; i++) {
    a_isInterfaceChild(i) = false;
    a_isInterfaceParent(i) = false;
  }

  mDeallocate(m_noInterfaceChildren);
  mDeallocate(m_noInterfaceParents);

  if(m_interfaceChildren.size() != 0) {
    for(auto&& children : m_interfaceChildren) {
      delete children;
    }
    m_interfaceChildren.clear();
  }
  if(m_interfaceParents.size() != 0) {
    for(auto&& parents : m_interfaceParents) {
      delete parents;
    }
    m_interfaceParents.clear();
  }
}

/** \brief Identify the cells located at an interface
 *
 * Diagonals must have been created already from the solver class.
 * \author Georg Eitel-Amor
 */
template <MInt nDim>
void LbSolver<nDim>::buildInterfaceCells() {
  TRACE();

  m_log << "    - building interface cells" << endl;

  const MInt interfaceLevels = maxLevel() - minLevel();

  if(interfaceLevels == 0) return;

  MInt noCells = m_cells.size();

  MInt pointedId = -1;
  MInt parentId = -1;
  MInt maxNoNghbrs;

  MBool emptyParent;

  if constexpr(nDim == 3) {
    maxNoNghbrs = 26;
  } else {
    maxNoNghbrs = 8;
  }

  // 1.a) ----------- Find interface children ----------------

  m_log << "      * finding interface children" << endl;

  mAlloc(m_noInterfaceChildren, interfaceLevels, "m_noInterfaceChildren", 0, AT_);

  // a) mark interface cells
  for(MInt i = 0; i < noInternalCells(); i++) {
    // skip cells on lowest level
    if(a_level(i) == minLevel()) continue;

    // // skip boundary cells
    // if (a_onlyBoundary(i))
    // 	  continue;

    // if cell has a parent cell
    if(c_parentId(i) > -1) {
      parentId = c_parentId(i);
      pointedId = i;

      const MInt currentInterfaceLevel = a_level(pointedId) - minLevel() - 1;

      for(MInt dist = 0; dist < maxNoNghbrs; dist++) {
        // if an equal-level neighbor doesn't exist but a parent-neighbor exists
        if(a_hasNeighbor(pointedId, dist) == 0 && a_hasNeighbor(parentId, dist)) {
          // if parent neighbor has no children an interface was found
          if(c_isLeafCell(c_neighborId(parentId, dist))) {
            // cell is interface cell
            // -> mark cell and increase counter
            a_isInterfaceChild(pointedId) = true;
            m_noInterfaceChildren[currentInterfaceLevel]++;
            break;
          }
        }
      }

      if(!m_correctInterfaceBcCells && a_isInterfaceChild(pointedId)) {
        // unmark cells which lie at a non-periodic boundary of the domain
        for(MInt dist = 0; dist < maxNoNghbrs; dist++) {
          if(a_hasNeighbor(pointedId, dist) == 0 && a_hasNeighbor(parentId, dist) == 0) {
            a_isInterfaceChild(pointedId) = false;
            m_noInterfaceChildren[currentInterfaceLevel]--;
            break;
          }
        }
      }
    }
  }

  for(MInt level = 0; level < interfaceLevels; level++)
    if(m_noInterfaceChildren[level] > 0)
      m_log << "        # found children for level " << minLevel() + level + 1 << ": " << m_noInterfaceChildren[level]
            << endl;


  // 1.b) ------------- write interface children to collector ------------

  m_log << "        # writing children to collector " << endl;

  for(MInt i = 0; i < noCells; i++) { // XXX exclude halo cellls?

    // skip all unmarked cells
    if(!a_isInterfaceChild(i)) continue;

    parentId = c_parentId(i);
    pointedId = i;

    // append cell to collector
    const MInt currentInterfaceLevel = a_level(i) - minLevel() - 1;
    m_interfaceChildren[currentInterfaceLevel]->append();

    auto& currentInterfaceCell =
        m_interfaceChildren[currentInterfaceLevel]->a[m_interfaceChildren[currentInterfaceLevel]->size() - 1];
    currentInterfaceCell.m_cellId = i;

    // Determine position in parent cell
    currentInterfaceCell.m_position = 0;
    for(MInt dim = 0; dim < nDim; dim++) {
      if(a_coordinate(pointedId, dim) < a_coordinate(parentId, dim)) {
        currentInterfaceCell.m_position += 0;
      } else {
        currentInterfaceCell.m_position += IPOW2(dim);
      }
    }
  }

  // 2.a) ----------- Find interface parents ----------------

  m_log << "      * finding interface parents" << endl;

  mAlloc(m_noInterfaceParents, interfaceLevels, "m_noInterfaceParents", 0, AT_);

  // 2.a1) mark interface parent cells
  for(MInt i = 0; i < noCells; i++) { // XXX exclude halo cellls?

    // skip finest cells
    if(a_level(i) == maxLevel()) continue;
    // skip non-parent cells
    if(c_isLeafCell(i)) continue;

    for(MInt k = 0; k < IPOW2(nDim); k++) {
      if(c_childId(i, k) < 0 || c_childId(i, k) >= noCells) {
        if(c_childId(i, k) != -1) cout << "STRANGE ID " << c_childId(i, k) << endl;
        continue;
      }
      if(a_isInterfaceChild(c_childId(i, k))) {
        // Cell is interface parent
        // -> mark cell and increase counter if number of children is correct
        if(c_noChildren(i) != IPOW2(nDim) && i < noInternalCells()) {
          stringstream errorMessage;
          errorMessage << " LbInterface::identifyInterfaceCells: Insufficient number of children in interface "
                          "parent cell with id "
                       << i << ", Exiting!";
          TERMM(1, errorMessage.str());
        }

        a_isInterfaceParent(i) = true;

        const MInt currentInterfaceLevel = a_level(i) - minLevel();
        m_noInterfaceParents[currentInterfaceLevel]++;

        break;
      }
    }
  }

  for(MInt level = 0; level < interfaceLevels; level++)
    if(m_noInterfaceParents[level] > 0)
      m_log << "        # found parents for level " << minLevel() + level + 1 << ": " << m_noInterfaceParents[level]
            << endl;


  m_log << "      * adding halo cells " << endl;

  // 2.a2) add halo cells
  for(MInt i = noInternalCells(); i < noCells; i++) {
    // skip finest cells
    if(a_level(i) == maxLevel()) continue;

    // skip non-parent cells
    if(c_isLeafCell(i)) continue;

    // if children are missing, skip the parent cell
    emptyParent = false;
    for(MInt child = 0; child < IPOW2(nDim); child++) {
      if(c_childId(i, child) < 0) {
        emptyParent = true;
        break;
      }
    }
    if(emptyParent) continue;

    for(MInt k = 0; k < maxNoNghbrs; k++) {
      if(a_hasNeighbor(i, k) > 0 && c_isLeafCell(c_neighborId(i, k))) {
        a_isInterfaceParent(i) = true;

        const MInt currentInterfaceLevel = a_level(i) - minLevel();
        m_noInterfaceParents[currentInterfaceLevel]++;
        break;
      }
    }
  }

  for(MInt level = 0; level < interfaceLevels; level++)
    if(m_noInterfaceParents[level] > 0)
      m_log << "        # new number of parents for level " << minLevel() + level + 1 << ": "
            << m_noInterfaceParents[level] << endl;

  m_log << "        # writing parents to collector " << endl;

  // 2.b) ------------- write interface parents to collector --------------
  //  for (MInt i = 0; i < noInternalCells(); i++){
  for(MInt i = 0; i < noCells; i++) {
    pointedId = i;

    // skip all unmarked cells
    if(!a_isInterfaceParent(i)) continue;

    // append cell to collector
    const MInt currentInterfaceLevel = a_level(i) - minLevel();
    m_interfaceParents[currentInterfaceLevel]->append();

    auto& currentParentCell =
        m_interfaceParents[currentInterfaceLevel]->a[m_interfaceParents[currentInterfaceLevel]->size() - 1];
    currentParentCell.m_cellId = i;
  }
}

/**\brief sets the coefficients for interpolation at grid interfaces
 *
 *\author Georg Eitel-Amor
 */
template <MInt nDim>
void LbSolver<nDim>::setInterpolationCoefficients() {
  TRACE();

  m_log << "    - setting interpolation coefficients" << endl;

  MInt noInterpolationNghbrs;
  noInterpolationNghbrs = IPOW2(nDim);
  ScratchSpace<MFloat> interpolationCoefficients(IPOW2(nDim), noInterpolationNghbrs, AT_, "interpolationCoefficients");

  switch(m_interpolationType) {
    case LINEAR_INTERPOLATION: {
      for(MInt i = 0; i < IPOW2(nDim); i++) {
        for(MInt j = 0; j < noInterpolationNghbrs; j++) {
          interpolationCoefficients(i, j) = LbLatticeDescriptorBase<nDim>::linearInterpolationCoefficients(i, j);
        }
      }
      break;
    }
    case QUADRATIC_INTERPOLATION: {
      // not implemented yet
      TERMM(1,
            " In function LbInterfaceD2Q9::setInterpolationCoefficients() :  quadratic interpolation is not "
            "implemented yet!");
      break;
    }
    case CUBIC_INTERPOLATION: {
      // not implemented yet
      TERMM(1,
            " In function LbInterfaceD2Q9::setInterpolationCoefficients() :  cubic interpolation is not "
            "implemented yet!");
      break;
    }
    default: {
      TERMM(1, " In function LbInterfaceD2Q9::setInterpolationCoefficients() :  Unknown interpolation type!");
    }
  }

  MInt position = 0;
  for(MInt i = 0; i < (maxLevel() - minLevel()); i++) {
    for(MInt j = 0; j < m_interfaceChildren[i]->size(); j++) {
      position = m_interfaceChildren[i]->a[j].m_position;
      for(MInt k = 0; k < noInterpolationNghbrs; k++) {
        m_interfaceChildren[i]->a[j].m_interpolationCoefficients[k] = interpolationCoefficients(position, k);
      }
    }
  }
}

/*\brief sets the neighbors for interpolation at grid interfaces
 *
 *\author Georg Eitel-Amor
 */
template <MInt nDim>
void LbSolver<nDim>::setInterpolationNeighbors() {
  TRACE();

  m_log << "    - setting interpolation neighbors" << endl;

  MInt position = 0;
  MInt currentId = 0;

  // ------- 2D ---------
  if constexpr(nDim == 2) {
    // by now only linear interpolation
    switch(m_interpolationType) {
      case LINEAR_INTERPOLATION: {
        for(MInt i = 0; i < maxLevel() - minLevel(); i++) {
          for(MInt j = 0; j < m_noInterfaceChildren[i]; j++) {
            currentId = m_interfaceChildren[i]->a[j].m_cellId;
            position = m_interfaceChildren[i]->a[j].m_position;
            for(MInt k = 0; k < IPOW2(nDim); k++) {
              const MInt dir = LbLatticeDescriptorBase<nDim>::intNghbrArray(position, k);
              if(dir < IPOW3[nDim] - 1) {
                // Access interpolation neighbor via child or parent cell
                // (depending on which is connected)
                if(a_hasNeighbor(currentId, dir) == 0) {
                  m_interfaceChildren[i]->a[j].m_interpolationNeighbors[k] = c_neighborId(c_parentId(currentId), dir);
                } else {
                  m_interfaceChildren[i]->a[j].m_interpolationNeighbors[k] = c_parentId(c_neighborId(currentId, dir));
                }
              } else {
                m_interfaceChildren[i]->a[j].m_interpolationNeighbors[k] = c_parentId(currentId);
              }
            }
          }
        }
        break;
      }
      case QUADRATIC_INTERPOLATION:
        // not implemented yet
        TERMM(1,
              " In function LbInterfaceD2Q9::setInterpolationCoefficients() :  quadratic interpolation not "
              "implemented yet!");
        break;
      case CUBIC_INTERPOLATION:
        // not implemented yet
        TERMM(1,
              " In function LbInterfaceD2Q9::setInterpolationCoefficients() :  cubic interpolation not "
              "implemented yet!");
        break;
      default: {
        TERMM(1, " In function LbInterfaceD2Q9::setInterpolationCoefficients() :  Unknown interpolation type!");
      }
    }
  }
  // ------- 3D ---------
  else {
    // by now only linear interpolation
    for(MInt i = 0; i < maxLevel() - minLevel(); i++) {
      for(MInt j = 0; j < m_noInterfaceChildren[i]; j++) {
        currentId = m_interfaceChildren[i]->a[j].m_cellId;
        position = m_interfaceChildren[i]->a[j].m_position;
        for(MInt k = 0; k < IPOW2(nDim); k++) {
          const MInt dir = LbLatticeDescriptorBase<nDim>::intNghbrArray(position, k);
          if(dir < IPOW3[nDim] - 1) {
            // interpolation neighbor is neighbor of parent
            const MInt parentId = c_parentId(currentId);
            const MInt nghbrId = c_neighborId(parentId, dir);
            m_interfaceChildren[i]->a[j].m_interpolationNeighbors[k] = nghbrId;
          } else {
            // interpolation neighbor is own parent
            m_interfaceChildren[i]->a[j].m_interpolationNeighbors[k] = c_parentId(currentId);
          }
        }
      }
    }
  }
}

/** \brief  Re-Set interpolation neighbors for interface boundary cells
 *  \author Miro Gondrum
 *  \date   15.10.2020
 */
template <MInt nDim>
void LbSolver<nDim>::setInterpolationNeighborsBC() {
  switch(nDim) {
    case 2: { //--2D------------------------------------------------------------
      // TODO labels:LB,toenhance needs to be extended to 2D
      std::cout << "WARNING: Refined boundary cells are problematic in 2D! "
                << "This is not implemented, yet." << std::endl;
      break;
    }
    case 3: { //--3D------------------------------------------------------------
      for(MInt i = 0; i < (maxLevel() - minLevel()); i++) {
        // for(MInt j = 0; j < m_interfaceChildren[i]->size(); j++) {
        for(MInt j = 0; j < m_noInterfaceChildren[i]; j++) {
          const MInt cellId = m_interfaceChildren[i]->a[j].m_cellId;
          // check: if at least one interpolNghbr is missing or is a BC cell
          //        if this is the case -> adjust stencil
          // adjustment is for now done by simply shifting the complete stencil
          // befor rechecking whether shifted one is working. Shifting is
          // maximal done once for each Cartesian direction. If then no correct
          // stencil is found an error is called !
          MBool performCorrection = false;
          constexpr MInt maxNoNghbrs = IPOW3[nDim] - 1;
          std::array<MBool, maxNoNghbrs> validShiftDirection;
          validShiftDirection.fill(true);
          for(MInt k = 0; k < m_interfaceChildren[i]->a[j].m_noInterpolationNeighbors; k++) {
            const MInt nghbrId = m_interfaceChildren[i]->a[j].m_interpolationNeighbors[k];
            if(nghbrId < 0 || !a_isActive(nghbrId) || a_isBndryCell(nghbrId)) {
              performCorrection = true;
              // check which Cartesian directions are non-valid for shift
              // ..from parent cell's point of view or from original interpolNghbr's point of view
              const MInt srcId = (nghbrId < 0) ? c_parentId(cellId) : nghbrId;
              // check which direction are non-valid for shift
              for(MInt d = 0; d < maxNoNghbrs; d++) {
                if(a_hasNeighbor(srcId, d)) {
                  const MInt srcNghbrId = c_neighborId(srcId, d);
                  if(!a_isActive(srcNghbrId) || a_isBndryCell(srcNghbrId)) {
                    validShiftDirection[d] = false;
                  }
                } else {
                  validShiftDirection[d] = false;
                }
              }
            }
          }
          // now shift in first validShiftDirection
          if(performCorrection) {
            MBool validShift = false;
            for(MInt d = 0; d < maxNoNghbrs; d++) {
              if(validShiftDirection[d]) {
                validShift = true;
                for(MInt k = 0; k < m_interfaceChildren[i]->a[j].m_noInterpolationNeighbors; k++) {
                  const MInt nghbrId = m_interfaceChildren[i]->a[j].m_interpolationNeighbors[k];
                  if(nghbrId > -1) {
                    m_interfaceChildren[i]->a[j].m_interpolationNeighbors[k] = c_neighborId(nghbrId, d);
                  } else {
                    std::array<MInt, nDim> newDirVec;
                    const MInt position = m_interfaceChildren[i]->a[j].m_position;
                    for(MInt l = 0; l < nDim; l++) {
                      // shift original direction vector by shift vector
                      // TODO labels:LB dxqy: get orignal direction less deterministic !!
                      newDirVec[l] = LbLatticeDescriptorBase<nDim>::idFld(
                                         LbLatticeDescriptorBase<nDim>::intNghbrArray(position, k), l)
                                     + (LbLatticeDescriptorBase<nDim>::idFld(d, l) - 1);
                    }
                    const MInt newDir = LbLatticeDescriptorBase<3>::dirFld(newDirVec[0], newDirVec[1],
                                                                           newDirVec[2]); // TODO labels:LB dxqy
                    const MInt parentId = c_parentId(cellId);
                    if(newDir == IPOW3[nDim] - 1) {
                      m_interfaceChildren[i]->a[j].m_interpolationNeighbors[k] = parentId;
                    } else {
                      m_interfaceChildren[i]->a[j].m_interpolationNeighbors[k] = c_neighborId(parentId, newDir);
                    }
                  }
                }
                break;
              }
            }
            if(!validShift) {
              std::stringstream err;
              err << "ERROR: interpolation stencil for BC interface could not be corrected!";
              err << " [globalId: " << c_globalId(cellId) << "]" << std::endl;
              TERMM(1, err.str());
            } else {
              // store to reset interpolation coefficients later
              InterfaceLink tmp = {i, j};
              m_interfaceChildrenBc.push_back(tmp);
            }
          } // end correction
        }
      }
      break;
    }
    default: {
      TERMM(1, "Only nDim=2 and nDim=3 allowed!");
    }
  }
}

/** \brief  Re-set interpolation coefficients
 *  \author Miro Gondrum
 *  \date   15.10.2020
 */
template <MInt nDim>
void LbSolver<nDim>::setInterpolationCoefficientsBC() {
  // if(m_interfaceChildrenBc.size() <= 0) return;
  // constexpr MInt noInterpolationNghbrs = IPOW2(nDim);
  switch(nDim) {
    case 2: { //--2D------------------------------------------------------------
      // TODO labels:LB,toenhance needs to be extended to 2D
      // A warning is already called in setInterpolationNeighborsBC()
      break;
    }
    case 3: { //--3D------------------------------------------------------------
      for(auto link : m_interfaceChildrenBc) {
        const MInt i = link.levelId;
        const MInt j = link.interfaceId;
        // weigths are determined for a 8 interpolation neighbors, which are
        // ordered in a cube based on tri-linear interpolation
        // 1) determine 3 weigths for first nghbr for each Cartesian dir
        std::array<MFloat, nDim> c0;
        const MInt cellId = m_interfaceChildren[i]->a[j].m_cellId;
        const MInt cellId0 = m_interfaceChildren[i]->a[j].m_interpolationNeighbors[0];
        for(MInt d = 0; d < nDim; d++) {
          const MInt cellId1 = m_interfaceChildren[i]->a[j].m_interpolationNeighbors[IPOW2(d)];
          const MFloat x = a_coordinate(cellId, d);
          const MFloat x0 = a_coordinate(cellId0, d);
          const MFloat x1 = a_coordinate(cellId1, d);
          c0[d] = (x1 - x) / (x1 - x0);
        }
        // 2) assemble final weigth by combination of c0
        for(MInt m = 0; m < 8; m++) {
          MFloat weight = 1.0;
          MInt tmp = m;
          for(MInt d = nDim - 1; d > -1; d--) {
            const MInt chk0 = tmp / IPOW2(d);
            tmp = tmp % IPOW2(d);
            const MBool chk = (chk0 == 0);
            weight *= (chk) ? c0[d] : 1 - c0[d];
          }
          m_interfaceChildren[i]->a[j].m_interpolationCoefficients[m] = weight;
        }
      }
      break;
    }
    default: {
      TERMM(1, "Only nDim=2 and nDim=3 allowed!");
    }
  }
}

/*! \brief load variables for the specified timeStep
 *
 * \author A. Niemoeller
 * \date 11.12.2013
 *
 * loads the variables from a restartFile for a given timeStep
 *
 * \param[in] timeStep timestep of restartfile
 *
 */
template <MInt nDim>
void LbSolver<nDim>::loadSampleVariables(MInt timeStep) {
  MString filename = outputDir() + "PV_";
  MChar buf[10];
  sprintf(buf, "%d", timeStep);
  filename.append(buf);
  filename += ParallelIo::fileExt();

  loadRestartWithoutDistributionsPar(filename.c_str());
}

/*! \brief read only access to primitive variables of a single cell
 *
 * \author A. Niemoeller
 * \date 11.12.2013
 *
 * get cell variables of a single cell
 *
 * \param[in] cellId cell that is accessed
 * \param[in,out] vars pointer to the variables
 *
 */
template <MInt nDim>
void LbSolver<nDim>::getSampleVariables(MInt cellId, const MFloat*& vars) {
#ifdef WAR_NVHPC_PSTL
  // TODO labels:LB,GPU resolve following TERMM
  TERMM(1, "getSampleVariables: var runs out of scope after this function call!");
  if(m_isThermal) {
    MFloat var[nDim + 2] = {F0};
    for(MInt d = 0; d < nDim + 2; d++) {
      var[d] = a_variable(cellId, d);
    }
    vars = var;
  } else {
    MFloat var[nDim + 1] = {F0};
    for(MInt d = 0; d < nDim + 1; d++) {
      var[d] = a_variable(cellId, d);
    }
    vars = var;
  }
#else
  vars = a_variables_ptr(cellId);
#endif
}

template <MInt nDim>
void LbSolver<nDim>::getSampleVariables(const MInt cellId, std::vector<MFloat>& vars) {
  const MInt noVars = vars.size();
  ASSERT(noVars <= m_noVariables, "noVars > m_noVariables");
  for(MInt varId = 0; varId < noVars; varId++) {
    vars[varId] = a_variable(cellId, varId);
  }
}

template <MInt nDim>
MBool LbSolver<nDim>::getSampleVarsDerivatives(const MInt cellId, std::vector<MFloat>& vars) {
  auto isValid = [&](const MInt l_cellId, const MInt l_dir) -> MInt {
    if(l_cellId != -1) {
      const MInt nghbrId = c_neighborId(l_cellId, l_dir);
      if(nghbrId > -1 && a_isActive(nghbrId)) {
        return nghbrId;
      }
    }
    return -1;
  };
  const MFloat F1bdx = FFPOW2(maxLevel() - a_level(cellId)); // LB units
  for(MInt dir = 0; dir < nDim; dir++) {
    const MInt dir0 = 2 * dir;
    const MInt dir1 = 2 * dir + 1;
    // 1) set coefficients
    constexpr MInt noCellsInStencil = 5;
    std::array<MFloat, noCellsInStencil> coef;
    std::array<MInt, noCellsInStencil> cellIds;
    // 1.1) get cellIds and coefficients
    cellIds[1] = isValid(cellId, dir0);     // left
    cellIds[0] = isValid(cellIds[1], dir0); // leftleft
    cellIds[2] = cellId;                    // middle
    cellIds[3] = isValid(cellId, dir1);     // right
    cellIds[4] = isValid(cellIds[3], dir1); // rightright
    getDerivativeStencilAndCoefficient(cellIds, coef);
    // 1.2) scale coefficients with correct dx in case of refined
    for(MInt i = 0; i < noCellsInStencil; i++) {
      coef[i] *= F1bdx;
    }
    // 2) calculate gradients
    // 2.1) .. velocity
    for(MInt veloId = 0; veloId < nDim; veloId++) {
      const MInt index = nDim * veloId + dir;
      vars[index] = 0.0;
      for(MInt i = 0; i < noCellsInStencil; i++) {
        vars[index] += cellIds[i] < 0 ? 0.0 : coef[i] * a_variable(cellIds[i], veloId);
      }
    }
    // 2.2) .. density
    for(MInt varId = nDim; varId < m_noVariables; varId++) {
      const MInt index = nDim * varId + dir;
      vars[index] = 0.0;
      for(MInt i = 0; i < noCellsInStencil; i++) {
        vars[index] += cellIds[i] < 0 ? 0.0 : coef[i] * a_variable(cellIds[i], varId);
      }
    }
    // 2.3) .. alpha?
    if(m_isEELiquid) {
      const MInt index = m_noVariables * nDim + dir;
      for(MInt i = 0; i < noCellsInStencil; i++) {
        vars[index] += cellIds[i] < 0 ? 0.0 : coef[i] * a_alphaGas(cellIds[i]);
      }
    }
  }
  return true;
}

template <MInt nDim>
void LbSolver<nDim>::storeOldDistributions() {
  if(!m_cells.savePrevVars()) {
    TERMM(1, "Cell collector not configured to stor previous values!");
  }
  for(MInt cellId = 0; cellId < a_noCells(); cellId++) {
    for(MInt distr = 0; distr < m_noDistributions; distr++) {
      a_previousDistribution(cellId, distr) = a_oldDistribution(cellId, distr);
    }
  }
}

template <MInt nDim>
void LbSolver<nDim>::storeOldVariables() {
  if(!m_cells.savePrevVars()) {
    TERMM(1, "Cell collector not configured to stor previous values!");
  }
  for(MInt cellId = 0; cellId < a_noCells(); cellId++) {
    a_previousVariable(cellId, PV->RHO) = a_oldVariable(cellId, PV->RHO);
    a_previousVariable(cellId, PV->U) = a_oldVariable(cellId, PV->U);
    a_previousVariable(cellId, PV->V) = a_oldVariable(cellId, PV->V);
    a_previousVariable(cellId, PV->W) = a_oldVariable(cellId, PV->W);
  }
}

template <MInt nDim>
void LbSolver<nDim>::initNu(const MInt cellId, const MFloat nu) {
  a_nu(cellId) = nu;
  if(m_cells.saveNuT()) {
    a_nuT(cellId) = 0.0;
    if(m_cells.saveOldNu()) {
      a_oldNuT(cellId) = 0.0;
    }
  }
  if(m_cells.saveOldNu()) {
    a_oldNu(cellId) = nu;
  }
}

/// Perform the solution step of the lattice Boltzmann method
template <MInt nDim>
MBool LbSolver<nDim>::solutionStep() {
  TRACE();

  RECORD_TIMER_START(m_t.solutionStep);

  // calls the exchange method
  RECORD_TIMER_START(m_t.exchange);
  (this->*m_exchangeMethod)(0);
  RECORD_TIMER_STOP(m_t.exchange);

  // update macroscopic variables
  if(m_updateMacroscopicLocation == PRECOLLISION) {
    updateVariablesFromOldDist_preCollision();
  }

  // update/reset viscosity
  updateViscosity();

  // pre-collison source term call
  RECORD_TIMER_START(m_t.srcTerms);
  preCollisionSrcTerm();
  RECORD_TIMER_STOP(m_t.srcTerms);

  // collision step
  RECORD_TIMER_START(m_t.collision);
  (this->*m_solutionStepMethod)();
  RECORD_TIMER_STOP(m_t.collision);

  // post-collison source term call
  RECORD_TIMER_START(m_t.srcTerms);
  postCollisionSrcTerm();
  RECORD_TIMER_STOP(m_t.srcTerms);

  // post-collision boundary conditions
  RECORD_TIMER_START(m_t.collisionBC);
  postCollisionBc();
  RECORD_TIMER_STOP(m_t.collisionBC);

  // calls the exchange method
  RECORD_TIMER_START(m_t.exchange);
  (this->*m_exchangeMethod)(1);
  RECORD_TIMER_STOP(m_t.exchange);

  // propagation step
  RECORD_TIMER_START(m_t.propagation);
  (this->*m_propagationStepMethod)();
  RECORD_TIMER_STOP(m_t.propagation);

  // TODO labels:LB @johannes/julian: I made this step in updateVariablesFromOldDist_preCollison. Preference for here?
  if(m_isInitRun) {
    this->initRunCorrection();
  }

  // set number outer Bnd cells for lblpt coupler
  m_noOuterBndryCells = this->m_bndCnd->m_bndCells.size();
  // post-propagation boundary conditions
  RECORD_TIMER_START(m_t.propagationBC);
  postPropagationBc();
  RECORD_TIMER_STOP(m_t.propagationBC);

  writeInfo();

  if(m_updateMacroscopicLocation == POSTPROPAGATION) {
    updateVariablesFromOldDist();
  }

  RECORD_TIMER_STOP(m_t.solutionStep);

  return true;
}

template <MInt nDim>
void LbSolver<nDim>::writeInfo() {
  TRACE();

  if(!isActive()) return;
  RECORD_TIMER_START(m_t.residual);

  if(globalTimeStep % m_residualInterval == 0) {
    maxResidual();
  }
  RECORD_TIMER_STOP(m_t.residual);
}

/// \brief Initialize the solver
template <MInt nDim>
void LbSolver<nDim>::initSolver() {
  TRACE();

  // set cell properties
  resetActiveCellList();

  // reset interface BC cell interpolation neighbors
  // TODO labels:LB this can not be called by treatInterfaceCell since a_isBndry is not
  // set correctly in advance.
  correctInterfaceBcCells();

  if(m_externalForcing) {
    initPressureForce();
  }
  if(m_EELiquid.gravity || m_externalForcing) {
    initVolumeForces();
  }

  if(grid().isActive()) {
    initSrcTermController();
    (this->*m_initializeMethod)();
    m_bndCnd->initializeBcData();
    initSrcTerms();
  }

  m_log << endl << endl;
}


/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/**
 * \brief
 * \author Moritz Waldmann
 */
template <MInt nDim>
void LbSolver<nDim>::prepareAdaptation() {
  TRACE();
  if(!m_refinedParents.empty()) m_refinedParents.clear();
}

/** \brief: Compute the sensors for adaptation
 *
 * \author: Philipp Brokof, Moritz Waldmann
 */
template <MInt nDim>
void LbSolver<nDim>::setSensors(std::vector<std::vector<MFloat>>& sensors, std::vector<MFloat>& sensorWeight,
                                std::vector<std::bitset<64>>& sensorCellFlag, std::vector<MInt>& sensorSolverId) {
  cerr0 << "Set " << this->m_noSensors << " Sensors for Adaptation" << endl;

  TRACE();
  ASSERT(m_freeIndices.empty(), "");

  // The sensors are added at the end of the previous sensor-vectors!
  // Implement sensor in the for loop below in this order.
  const auto sensorOffset = (signed)sensors.size();
  ASSERT(sensorOffset == 0 || grid().raw().treeb().noSolvers() > 1, "");
  sensors.resize(sensorOffset + this->m_noSensors, vector<MFloat>(grid().raw().m_noInternalCells, F0));
  sensorWeight.resize(sensorOffset + this->m_noSensors, -1);
  sensorCellFlag.resize(grid().raw().m_noInternalCells, sensorOffset + this->m_noSensors);
  sensorSolverId.resize(sensorOffset + this->m_noSensors, solverId());
  ASSERT(sensorOffset + this->m_noSensors < CartesianGrid<nDim>::m_maxNoSensors, "Increase bitset size!");

  if(!isActive()) {
    for(MInt sen = 0; sen < this->m_noSensors; sen++) {
      sensorWeight[sensorOffset + sen] = this->m_sensorWeight[sen];
    }
    return;
  }

  // only set sensors if the adaptation was globally triggered in buildLevelSetTubeCG!
  if(m_adaptation) {
    // Perform exchange to update Halos that are involved in gradient computation
    if(noNeighborDomains() > 0) {
      exchange(1);
      exchangeOldDistributions();
    }

    if(domainId() == 0) {
      cerr << "Setting sensor for the lb-solver adaptation!" << endl;
    }

    for(MInt sen = 0; sen < this->m_noSensors; sen++) {
      (this->*(this->m_sensorFnPtr[sen]))(sensors, sensorCellFlag, sensorWeight, sensorOffset, sen);
    }
  }
}


/**
 * \brief needed for sampling
 * \author Benyamin Krisna
 */
template <MInt nDim>
void LbSolver<nDim>::getSolverSamplingProperties(std::vector<MInt>& samplingVarIds,
                                                 std::vector<MInt>& noSamplingVars,
                                                 std::vector<std::vector<MString>>& samplingVarNames,
                                                 const MString featureName) {
  TRACE();

  // Read sampling variable names
  std::vector<MString> varNamesList;
  MInt noVars = readSolverSamplingVarNames(varNamesList, featureName);

  // Set default sampling variables if none specified
  if(noVars == 0) {
    varNamesList.emplace_back("LB_PV");
    noVars = 1;
  }

  for(MInt i = 0; i < noVars; i++) {
    const MInt samplingVar = string2enum(varNamesList[i]);
    std::vector<MString> varNames;

    auto samplingVarIt = std::find(samplingVarIds.begin(), samplingVarIds.end(), samplingVar);
    if(samplingVarIt != samplingVarIds.end()) {
      TERMM(1, "Sampling variable '" + varNamesList[i] + "' already specified.");
    }

    switch(samplingVar) {
      case LB_PV: {
        ASSERT(m_noVariables == nDim + 1,
               "Error: additional primitive variables not supported for sampling data feature yet");

        samplingVarIds.push_back(LB_PV);
        noSamplingVars.push_back(m_noVariables + 1); //+1 because of addition of p

        varNames.resize(m_noVariables + 1); //+1 because of addition of p
        varNames[PV->U] = "u";
        varNames[PV->V] = "v";
        if constexpr(nDim == 3) {
          varNames[PV->W] = "w";
        }
        varNames[PV->P] = "p"; // add a p variable to be sampled, that can be induced from rho
        varNames[PV->RHO] = "rho";

        samplingVarNames.push_back(varNames);
        break;
      }
      default: {
        TERMM(1, "Unknown sampling variable: " + varNamesList[i]);
        break;
      }
    }
  }
}

/**
 * \brief
 * \author Moritz Waldmann
 */
template <MInt nDim>
void LbSolver<nDim>::finalizeAdaptation() {
  TRACE();

  if(!grid().isActive()) return;
  if(m_maxNoSets > -1) {
    mDeallocate(m_sendBufferMB);
    mDeallocate(m_receiveBufferMB);
    mDeallocate(m_associatedBodyIds);
    mDeallocate(m_levelSetValues);
    mDeallocate(m_G0CellMapping);
    mDeallocate(m_isG0CandidateOfSet);

    // TODO: noNeighborDomains > 0 instead?
    if(noDomains() > 1) {
      MInt sumwin = 0;
      MInt sumhalo = 0;
      for(MInt d = 0; d < noNeighborDomains(); d++) {
        sumwin += noWindowCells(d);
        sumhalo += noHaloCells(d);
      }
      mAlloc(m_sendBufferMB, sumwin, "m_sendBufferMB", 0, AT_);
      mAlloc(m_receiveBufferMB, sumhalo, "m_receiveBufferMB", 0, AT_);
    }

    // body properties
    mAlloc(m_associatedBodyIds, m_noLevelSetsUsedForMb * a_noCells(), "m_associatedBodyIds", -1, AT_);
    mAlloc(m_levelSetValues, m_noLevelSetsUsedForMb * a_noCells(), "m_levelSetValues", F0, AT_);
    mAlloc(m_G0CellMapping, a_noCells(), "m_G0CellMapping", -1, AT_);

    MInt startSet = m_levelSetId; // if bodyToSetMode = 11 the 0th set is the collected set, which should be skiped
    MInt arrayLength = m_maxNoSets - startSet;
    mAlloc(m_isG0CandidateOfSet, a_noCells(), arrayLength, "m_isG0CandidateOfSet", false, AT_);

    // TODO: noNeighborDomains > 0 instead?
    if(noDomains() > 1) {
      grid().updateLeafCellExchange();
    }

    resetActiveCellList();
  }

  // These status flag distinguishes if an adaptation was performed or not.
  // If adaptation was performed, recalcIds need to be passed to saveUVWRhoTOnlyPar
  // because the ids in the newly written grid file differ from the
  // ids used in solver/grid.
  m_adaptationSinceLastRestart = true;
  m_adaptationSinceLastSolution = true;
}

/**
 * \brief: Reinit lb solver after adaptation
 *
 * \details: - update cell lists that are stored in the solver
 *           - restart communication
 *           - restart boundary condition object
 *           - call initialization of refined cells
 *
 * \authors: Sohel Herff
 *           Tim Wegmann (level set)
 *           Philipp Brokof, Moritz Waldmann (lb)
 */
template <MInt nDim>
void LbSolver<nDim>::postAdaptation() {
  TRACE();

  for(MInt i = 0; i < grid().raw().noNeighborDomains(); i++) {
    for(MInt j = 0; j < grid().raw().noHaloCells(i); j++) {
      MInt gridId = grid().raw().m_haloCells[i][j];
      MInt l_solverId = this->grid().tree().grid2solver(gridId);
      if(l_solverId < 0) {
        continue;
      }

      if(!grid().raw().treeb().solver(gridId, m_solverId)) {
        cerr << "removing a Halo-cell " << endl;
        this->removeCellId(l_solverId);
      }
    }
  }

  // In meshAdaptation() the order of window and halo cells can change
  // even if no adaptation happens
  this->compactCells();

  if(!g_multiSolverGrid) {
    for(MInt gridCellId = 0; gridCellId < grid().raw().treeb().size(); gridCellId++) {
      ASSERT(grid().tree().solver2grid(gridCellId) == gridCellId, "");
      ASSERT(grid().tree().grid2solver(gridCellId) == gridCellId, "");
    }
  }

  grid().updateOther();
  updateDomainInfo(grid().domainId(), grid().noDomains(), grid().mpiComm(), AT_);

  if(!grid().isActive()) return;

  m_cells.size(c_noCells());
  m_freeIndices.clear();

  if(!g_multiSolverGrid) ASSERT(a_noCells() == c_noCells() && c_noCells() == grid().raw().treeb().size(), "");

  // FIXME labels:LB
  // set the cell weights.... i dont know why this is done by the solvers right now
  for(MInt i = 0; i < grid().raw().treeb().size(); i++) {
    grid().raw().treeb().weight(i) = 1;
  }

  updateCellCollectorFromGrid();

  grid().findEqualLevelNeighborsParDiagonal(false);

  m_isRefined = (grid().maxUniformRefinementLevel() < maxLevel());

  resetComm();

  if(m_reducedComm) {
    prepareCommunicationReduced();
  } else {
    prepareCommunicationNormal();
  }

  resetCellLists();

  // Do an exchange to have the right values in halo cells
  if(noNeighborDomains() > 0) {
    exchange(1);
    exchangeOldDistributions();
  }

  restartBndCnd();

  // Reset level set

  if(!grid().isActive()) return;
  if(m_maxNoSets > -1) {
    mDeallocate(m_associatedBodyIds);
    mDeallocate(m_levelSetValues);
    mDeallocate(m_G0CellMapping);
    mDeallocate(m_isG0CandidateOfSet);

    mAlloc(m_associatedBodyIds, m_noLevelSetsUsedForMb * a_noCells(), "m_associatedBodyIds", -1, AT_);
    mAlloc(m_levelSetValues, m_noLevelSetsUsedForMb * a_noCells(), "m_levelSetValues", F0, AT_);
    mAlloc(m_G0CellMapping, a_noCells(), "m_G0CellMapping", -1, AT_);

    MInt startSet = m_levelSetId; // if bodyToSetMode = 11 the 0th set is the collected set, which should be skiped
    MInt arrayLength = m_maxNoSets - startSet;
    mAlloc(m_isG0CandidateOfSet, a_noCells(), arrayLength, "m_isG0CandidateOfSet", false, AT_);
  }

  // Initialize refined cells after solver is updated because spatial interpolation
  // accross domain interfaces can become necessary. Also see refineCell.
  if(globalTimeStep < 1) {
    resetActiveCellList(1); // needs to be done here, since now bndry cells are known
    (this->*m_initializeMethod)();
    m_bndCnd->initializeBcData();
  } else {
    initializeRefinedCellsPerLevel();
    resetActiveCellList(1);
    m_bndCnd->initializeBcData();
    initializeNewInterfaceParents();
  }

  // New cells have been initialized, reset flag
  for(MInt i = 0; i < a_noCells(); i++) {
    a_hasProperty(i, SolverCell::WasNewlyCreated) = false;
  }

  ASSERT(m_freeIndices.empty(), "");
}

/**
 * \brief    This functions create solver entries for newly created
             cells due to adaptive mesh refinement.
 * \details  The cells are not initialized becuase spatial interpolation accros
             domains can be necessary and Halo/Window connection is not yet
             reinitialized. Therefore, cells are initialized in reinitAfterAdaptation
             by a call of initializeRefinedCellsPerLevel as soon as the Halo/Window
             connection is established.
 * \author   Philipp Brokof
 */
template <MInt nDim>
void LbSolver<nDim>::refineCell(const MInt gridCellId) {
  // Get solver cell id of the parent
  MInt solverParentId = grid().tree().grid2solver(gridCellId);
  if(!g_multiSolverGrid) ASSERT(solverParentId == gridCellId, "");

    ASSERT(grid().raw().a_hasProperty(gridCellId, Cell::WasRefined), "");

  // Loop over newly created children of the parent cell
  for(MInt child = 0; child < grid().m_maxNoChilds; child++) {
    // Get grid cell id of new child
    MInt gridChildId = grid().raw().a_childId(gridCellId, child);
    if(gridChildId == -1) continue;

    if(!g_multiSolverGrid) ASSERT(grid().raw().a_hasProperty(gridChildId, Cell::WasNewlyCreated), "");

    // If solver is inactive all cells musst be halo cells!
    if(!grid().isActive()) ASSERT(grid().raw().a_isHalo(gridChildId), "");
    // If child exists in grid but is not located inside solver geometry
    if(!grid().solverFlag(gridChildId, solverId())) continue;

    // Create cell id for new child in solver
    const MInt solverChildId = this->createCellId(gridChildId);
    a_hasProperty(solverChildId, SolverCell::WasNewlyCreated) = true;
  }
  m_refinedParents.insert(solverParentId);
}


/**
 * \brief
 * \author Moritz Waldmann
 *
 **/

template <MInt nDim>
void LbSolver<nDim>::removeChilds(const MInt gridCellId) {
  TRACE();

  // Get solver id of the parent that will be coarsen
  MInt solverParentId = grid().tree().grid2solver(gridCellId);
  ASSERT(solverParentId > -1 && solverParentId < m_cells.size(), "solverParentId is: " << solverParentId);
  ASSERT(c_noChildren(solverParentId) > 0, "");
  if(!g_multiSolverGrid) ASSERT(solverParentId == gridCellId, "");

  // SOLVER SPECIFIC PART
  // Initialize variables of parent from its children
  this->removeChildsLb(solverParentId);

  // Remove solver ids of deleted children
  for(MInt c = 0; c < grid().m_maxNoChilds; c++) {
    MInt childId = c_childId(solverParentId, c);
    if(childId < 0) continue;
    this->removeCellId(childId);
  }
  if(!g_multiSolverGrid) {
    ASSERT((grid().raw().treeb().size() - m_cells.size()) <= grid().m_maxNoChilds, "");
  }
}

/**
 * \brief
 **/
template <MInt nDim>
void LbSolver<nDim>::removeCell(const MInt gridCellId) {
  TRACE();

  // Get solver id of the parent that will be coarsen
  MInt cellId = grid().tree().grid2solver(gridCellId);
  if(cellId > -1) {
    ASSERT(cellId < m_cells.size(), "cellId is: " << cellId);
    ASSERT(c_noChildren(cellId) == 0, "");

    this->removeCellId(cellId);
  }
}

// this function should be moved to Solver as soon as cartesiansolver.h has been removed!!!
// this function should be moved to Solver as soon as cartesiansolver.h has been removed!!!
// this function should be moved to Solver as soon as cartesiansolver.h has been removed!!!
template <MInt nDim>
void LbSolver<nDim>::resizeGridMap() {
  grid().resizeGridMap(m_cells.size());
}

/**
 * \brief
 * \author Moritz Waldmann
 */

template <MInt nDim>
void LbSolver<nDim>::swapCells(const MInt cellId0, const MInt cellId1) {
  if(cellId1 == cellId0) return;

  for(MInt v = 0; v < m_noVariables; v++) {
    std::swap(a_variable(cellId1, v), a_variable(cellId0, v));
    std::swap(a_oldVariable(cellId1, v), a_oldVariable(cellId0, v));
  }
  for(MInt dir = 0; dir < m_noDistributions; dir++) {
    std::swap(a_distribution(cellId1, dir), a_distribution(cellId0, dir));
    std::swap(a_oldDistribution(cellId1, dir), a_oldDistribution(cellId0, dir));
  }
  if(m_isThermal) {
    for(MInt dir = 0; dir < m_noDistributions; dir++) {
      std::swap(a_distributionThermal(cellId1, dir), a_distributionThermal(cellId0, dir));
      std::swap(a_oldDistributionThermal(cellId1, dir), a_oldDistributionThermal(cellId0, dir));
    }
  }
  if(m_isTransport) {
    for(MInt dir = 0; dir < m_noDistributions; dir++) {
      std::swap(a_distributionTransport(cellId1, dir), a_distributionTransport(cellId0, dir));
      std::swap(a_oldDistributionTransport(cellId1, dir), a_oldDistributionTransport(cellId0, dir));
    }
  }
  std::swap(m_cells.allProperties(cellId1), m_cells.allProperties(cellId0));
  std::swap(a_nu(cellId1), a_nu(cellId0));
  std::swap(a_kappa(cellId1), a_kappa(cellId0));
  std::swap(a_diffusivity(cellId1), a_diffusivity(cellId0));

  std::swap(a_bndId(cellId1), a_bndId(cellId0));
  std::swap(a_level(cellId1), a_level(cellId0));

  if(!m_refinedParents.empty()) {
    auto it0 = m_refinedParents.find(cellId0);
    auto it1 = m_refinedParents.find(cellId1);
    if(it0 != m_refinedParents.end() && it1 == m_refinedParents.end()) {
      // nothing to be done
    } else if(it0 != m_refinedParents.end()) {
      m_refinedParents.erase(it0);
      m_refinedParents.insert(cellId1);
    } else if(it1 != m_refinedParents.end()) {
      m_refinedParents.erase(it1);
      m_refinedParents.insert(cellId0);
    }
  }
}

/**
 * \brief
 * \author Lennart Schneiders
 */

template <MInt nDim>
void LbSolver<nDim>::swapProxy(const MInt cellId0, const MInt cellId1) {
  grid().swapGridIds(cellId0, cellId1);
}

/**
 * \brief
 * \author Moritz Waldmann
 */

template <MInt nDim>
void LbSolver<nDim>::resetComm() {
  if(m_reducedComm) {
    mDeallocate(m_noWindowDistDataPerDomain);
    mDeallocate(m_noHaloDistDataPerDomain);
    mDeallocate(m_nghbrOffsetsWindow);
    mDeallocate(m_nghbrOffsetsHalo);
    mDeallocate(m_sendBuffers);
    mDeallocate(m_receiveBuffers);
    mDeallocate(m_needsFurtherExchange);
    mDeallocate(m_windowDistsForExchange);
    mDeallocate(m_haloDistsForExchange);
    if(!m_nonBlockingComm) {
      mDeallocate(mpi_request);
    } else {
      mDeallocate(mpi_requestS);
      mDeallocate(mpi_requestR);
    }
  } else {
    mDeallocate(m_dataBlockSizes);
    mDeallocate(m_baseAddresses);
    mDeallocate(m_nghbrOffsetsWindow);
    mDeallocate(m_nghbrOffsetsHalo);
    mDeallocate(m_sendBuffers);
    mDeallocate(m_receiveBuffers);
    if(!m_nonBlockingComm) {
      mDeallocate(mpi_request);
    } else {
      mDeallocate(mpi_requestS);
      mDeallocate(mpi_requestR);
    }
  }
}

/**
 * \brief
 * \author Moritz Waldmann
 */
template <MInt nDim>
void LbSolver<nDim>::resetCellLists(MBool resize /*= true*/) {
  if(resize) {
    // resize list of active cells
    mDeallocate(m_activeCellList);
    mAlloc(m_activeCellList, grid().noCells(), "m_activeCellList", AT_);
  }

  if(m_isRefined) {
    resetInterfaceCells();
    initializeInterfaceCells();
    buildInterfaceCells();
    setInterpolationNeighbors();
    setInterpolationCoefficients();
  }
}

/**
 * \brief This function prepares a restart that is handled by the grid-controller!
 * \author Moritz Waldmann
 */

template <MInt nDim>
MBool LbSolver<nDim>::prepareRestart(MBool writeRestart, MBool& writeGridRestart) {
  TRACE();

  writeGridRestart = false;

  if(((globalTimeStep % m_restartInterval) == 0) || writeRestart) {
    writeRestart = true;

    if(m_adaptationSinceLastRestart) {
      writeGridRestart = true;
    }
  }

  return writeRestart;
}

/**
 * \brief This function resets the grid-trigger after a restart that is handled by the grid-controller!
 * \author Moritz Waldmann
 */

template <MInt nDim>
void LbSolver<nDim>::reIntAfterRestart(MBool doneRestart) {
  TRACE();

  if(doneRestart) {
    m_adaptationSinceLastRestart = false;
  }
}

/**
 * \brief This function writes restart that is handled by the grid-controller!
 * \author Moritz Waldmann
 */
template <MInt nDim>
void LbSolver<nDim>::writeRestartFile(const MBool writeRestart, const MBool /*writeBackup*/, const MString gridFileName,
                                      MInt* recalcIdTree) {
  TRACE();

  if(writeRestart) {
    stringstream fileName;
    fileName << outputDir() << "restart_" << getIdentifier() << globalTimeStep << ParallelIo::fileExt();

    std::vector<MInt> recalcCellIdsSolver(0);
    MInt noCells;
    MInt noInternalCellIds;
    std::vector<MInt> reorderedCellIds(0);
    this->calcRecalcCellIdsSolver(recalcIdTree, noCells, noInternalCellIds, recalcCellIdsSolver, reorderedCellIds);
    MInt* pointerRecalcIds = (recalcIdTree == nullptr) ? nullptr : recalcCellIdsSolver.data();
    saveRestartWithDistributionsPar(fileName.str().c_str(), gridFileName.c_str(), pointerRecalcIds);
  }
}

/**
 * \brief sets the cell-weight for balancing and a restarting
 * \author Moritz Waldmann
 */
template <MInt nDim>
void LbSolver<nDim>::setCellWeights(MFloat* solverCellWeight) {
  TRACE();
  const MInt noCellsGrid = grid().raw().treeb().size();
  const MInt offset = noCellsGrid * solverId();

  for(MInt cellId = 0; cellId < a_noCells(); cellId++) {
    const MInt gridCellId = grid().tree().solver2grid(cellId);
    const MInt id = gridCellId + offset;
    solverCellWeight[id] = F1; // 0.2;
  }
}

template <MInt nDim>
void LbSolver<nDim>::initializeNewInterfaceParents() {
  const MInt interfaceLevels = maxLevel() - minLevel();

  for(MInt level = 0; level < interfaceLevels; level++) {
    for(MInt p = 0; p < m_noInterfaceParents[level]; p++) {
      auto& parent = m_interfaceParents[level]->a[p];
      const MInt parentId = parent.m_cellId;
      if(!a_wasActive(parentId) && a_isActive(parentId)) {
        this->removeChildsLb(parentId);
      }
    }
  }
}

/** \brief Initialize refined cells after adaptation
 *
 * \details: Initialization is performed level-wise from maxUniformRefinementLevel
 *           to maxLevel. Perform an exchange after each level is initialized.
 *           This is necessary, becuase Halo cells of the previous level
 *           are needed to initialize refined cells on the next higher level.
 *
 * \author: Philipp Brokof, Moritz Waldmann
 */
template <MInt nDim>
void LbSolver<nDim>::initializeRefinedCellsPerLevel() {
  for(MInt level = grid().maxUniformRefinementLevel(); level < grid().maxLevel(); level++) {
    // Exchange to update Halo cells that are needed for spatial interpolation
    // TODO labels:LB Find out how to exchange cells only on the current level!
    if(noNeighborDomains() > 0) {
      exchange(1);
      exchangeOldDistributions();
    }

    for(auto&& parentId : m_refinedParents) {
      if(c_level(parentId) != level) continue;
      if(a_isHalo(parentId)) continue;
      if(!g_multiSolverGrid)
        if(!grid().raw().a_hasProperty(parentId, Cell::WasRefined)) continue;

      // Gather children of the refined parent cell and call initialization method
      std::vector<MInt> childIds(grid().m_maxNoChilds, -1);
      for(MInt child = 0; child < grid().m_maxNoChilds; child++) {
        MInt childId = c_childId(parentId, child);
        if(childId == -1) {
          continue;
        }
        if(!g_multiSolverGrid) ASSERT(grid().raw().a_hasProperty(childId, Cell::WasNewlyCreated), "");

        childIds[child] = childId;
      }
      this->refineCellLb(parentId, childIds.data());
    }
  }
}

template <MInt nDim>
void LbSolver<nDim>::initSolverSamplingVariables(const std::vector<MInt>& varIds,
                                                 const std::vector<MInt>& noSamplingVars) {
  TRACE();

  // TODO labels:LB,DLB enable use with balancing
  if(m_isInitSamplingVars) {
    m_log << "Sampling variables already initialized." << std::endl;
    return;
  }

  for(MUint i = 0; i < varIds.size(); i++) {
    MFloat** varPointer = nullptr;
    MInt dataBlockSize = noSamplingVars[i];

    // No additional storage for primitive variables
    if(varIds[i] == LB_PV) {
      dataBlockSize = 0;
    }

    if(dataBlockSize > 0) {
      mAlloc(varPointer, a_noCells(), dataBlockSize, "m_samplingVariables_" + std::to_string(varIds[i]), 0.0, AT_);
    }
    m_samplingVariables.push_back(varPointer);
  }

  // Set flag to avoid multiple initializations when using multiple sampling features
  m_isInitSamplingVars = true;
}

template <MInt nDim>
void LbSolver<nDim>::calcSamplingVariables(const std::vector<MInt>& varIds, const MBool NotUsed(exchange)) {
  for(MUint i = 0; i < varIds.size(); i++) {
    switch(varIds[i]) {
      case LB_PV: {
        // Nothing to do
        break;
      }
      default: {
        TERMM(1, "Invalid variable id");
        break;
      }
    }
  }
}

template <MInt nDim>
void LbSolver<nDim>::interpolateVariablesInCell(const MInt cellId, const MFloat* /*position*/,
                                                MFloat* interpolationResult) {
  // TODO labels:LB Add trilinear or least square interpolation
  interpolationResult[PV->RHO] = a_variable(cellId, PV->RHO);
  for(MInt d = 0; d < nDim; d++) {
    interpolationResult[PV->U + d] = a_variable(cellId, PV->U + d);
  }
  interpolationResult[PV->P] = interpolationResult[PV->RHO] * CSsq;
}

/// \brief Calculate the sampling variables at a given point in a cell.
template <MInt nDim>
void LbSolver<nDim>::calcSamplingVarAtPoint(const MFloat* point, const MInt id, const MInt sampleVarId, MFloat* state,
                                            const MBool NotUsed(interpolate)) {
  // Note: interpolateVariablesInCell() requires variables for all cells
  switch(sampleVarId) {
    case LB_PV: {
      interpolateVariablesInCell(id, point, state);
      break;
    }
    default: {
      TERMM(1, "Invalid variable id");
      break;
    }
  }
}


template <MInt nDim>
void LbSolver<nDim>::saveSolverSolution(MBool /*forceOutput*/, const MBool /*finalTimeStep*/) {
  if(globalTimeStep % m_solutionInterval == 0 && globalTimeStep >= m_solutionOffset) saveOutput();

  if(m_fftInterval > 0 && globalTimeStep % m_fftInterval == 0
     && (m_initMethod == "LB_TURBULENCE_ISOTROPIC_INIT" || m_initMethod == "FROM_RESTART_FILE")) {
    RECORD_TIMER_START(m_t.fft);
    computeFFTStatistics();
    RECORD_TIMER_STOP(m_t.fft);
  }
}

template <MInt nDim>
void LbSolver<nDim>::finalizeInitSolver() {
  if(!isActive()) return;
  if(m_updateMacroscopicLocation == POSTPROPAGATION && !m_restartFile) {
    updateVariablesFromOldDist();
  }

  outputInitSummary();
}

/**
 * \brief Output summary to console and log.
 *
 * \author Julian Vorspohl <j.vorspohl@aia.rwth-aachen.de>
 */
template <MInt nDim>
void LbSolver<nDim>::outputInitSummary() {
  using namespace maia::logtable;

  // TODO labels:LB,IO This should be moved into an appropriate sysEqn-like structure...
  const MInt noVars = nDim + 1;
  std::vector<MString> cvs{};
  cvs.push_back("rho");
  cvs.push_back("rho_u");
  cvs.push_back("rho_v");
  if(nDim == 3) {
    cvs.push_back("rho_w");
  }

  // INIT SUMMARY FRAME
  Frame summary(" SOLVER " + std::to_string(solverId()) + " INITIALIZATION SUMMARY AT TIME STEP "
                + std::to_string(globalTimeStep));

  // PROBLEM SUMMARY GROUP
  Group& problem = summary.addGroup(" PROBLEM SUMMARY");
  problem.addData("System of Equations", solverMethod());
  problem.addData("Number of dimensions", nDim);
  Data& vars = problem.addData("Number of variables", noVars);

  MBool first = true;
  for(auto&& cv : cvs) {
    if(first) {
      vars.addData("Conservative variable name(s)", cv);
    } else {
      vars.addData("", cv);
    }
    first = false;
  }

  problem.addData("Reynolds Number", m_Re);
  problem.addData("Mach Number", m_Ma);
  problem.addData("Reference length (LB)", m_referenceLength);

  problem.addData("Restart", m_restartFile);
  if(m_restartFile) {
    const MString restartFileName =
        restartDir() + "restart_" + getIdentifier() + std::to_string(m_restartTimeStep) + ParallelIo::fileExt();
    problem.addData("Initial condition", restartFileName);
  } else {
    Data& initData = problem.addData("Initial condition", m_initMethod);
    if(m_tanhInit) {
      initData.addData("tanhInit Re", m_initRe);
      initData.addData("tanhInit start time", m_initStartTime);
      initData.addData("tanhInit time", m_initTime);
    }
  }
  problem.addData("Start time (non-dimensionalized)", globalTimeStep);
  problem.addData("Final time (non-dimensionalized)", globalTimeStep + g_timeSteps);

  // DISCRETIZATION SUMMARY GROUP
  Group& discret = summary.addGroup("DISCRETIZATION SUMMARY");
  discret.addData("Number of distributions", m_noDistributions);
  discret.addData("Collision operator", solverMethod());
  discret.addBlank();
  Data& omega = discret.addData("Collision frequency", m_omega);
  omega.addData("Relaxation time", 1.0 / m_omega);
  omega.addData("Num. viscosity", m_nu);
  if(m_isThermal) {
    Data& omegaT = discret.addData("Collision frequency thermal", m_omegaT);
    omegaT.addData("Relaxation time thermal", 1.0 / m_omegaT);
    omegaT.addData("Num. diffusion thermal", m_kappa);
  }
  if(m_isTransport) {
    Data& omegaD = discret.addData("Collision frequency transport", m_omegaD);
    omegaD.addData("Relaxation time transport", 1.0 / m_omegaD);
    omegaD.addData("Num. diffusion transport", m_diffusivity);
  }
  discret.addBlank();
  discret.addData("External forcing", m_externalForcing);

  discret.addBlank();
  MString interfaceMethod = "FILIPPOVA";
  interfaceMethod = Context::getSolverProperty<MString>("interfaceMethod", m_solverId, AT_, &interfaceMethod);
  if(m_isRefined) {
    discret.addData("Interface method", interfaceMethod);
  }
  discret.addData("Adaptation init method", m_adaptationInitMethod);

  // PARALLELIZATION SUMMARY
  Group& parallel = summary.addGroup("PARALLELIZATION SUMMARY");
  parallel.addData("Domain id", domainId());
  parallel.addData("Number of neighbor domains", grid().noNeighborDomains());
  parallel.addData("Total number of domains", noDomains());

  // GRID SUMMARY LOCAL
  Group& gridLocal = summary.addGroup("GRID SUMMARY (LOCAL)");
  Data& minlvl = gridLocal.addData("Minimum grid level", grid().minLevel());
  minlvl.addData("cell length", c_cellLengthAtLevel(minLevel()));

  Data& maxlvl = gridLocal.addData("Maximum grid level", grid().maxLevel());
  maxlvl.addData("cell length", c_cellLengthAtLevel(maxLevel()));

  gridLocal.addBlank();
  Data& noCellsLocal = gridLocal.addData("Number of cells", a_noCells());
  noCellsLocal.addData("internal cells", grid().noInternalCells());
  noCellsLocal.addData("halo cells", a_noCells() - grid().noInternalCells());
  gridLocal.addData("Number of active cells", m_currentMaxNoCells);
  gridLocal.addData("Max. number of cells (collector size)", grid().raw().treeb().capacity());
  gridLocal.addData("Memory utilization", 100.0 * grid().noCells() / grid().raw().treeb().capacity());
  gridLocal.addBlank();
  if(m_isRefined) {
    const MInt levelRange = maxLevel() - minLevel();
    const MInt totalNoInterfaceChildren =
        std::accumulate(&m_noInterfaceChildren[0], &m_noInterfaceChildren[levelRange], 0);
    Data& ic = gridLocal.addData("Number of interface children", totalNoInterfaceChildren);
    for(MInt i = 0; i < levelRange; i++) {
      const MInt noInterfaceChildren = m_noInterfaceChildren[i];
      if(noInterfaceChildren > 0) {
        ic.addData("on level " + std::to_string(minLevel() + i + 2), noInterfaceChildren);
      }
    }
    const MInt totalNoInterfaceParents =
        std::accumulate(&m_noInterfaceParents[0], &m_noInterfaceParents[levelRange], 0);
    Data& ip = gridLocal.addData("Number of interface parents", totalNoInterfaceParents);
    for(MInt i = 0; i < levelRange; i++) {
      const MInt noInterfaceParents = m_noInterfaceParents[i];
      if(noInterfaceParents > 0) {
        ip.addData("on level " + std::to_string(minLevel() + i + 1), noInterfaceParents);
      }
    }
  }

  // Finally build string for summary frame
  MString s = summary.buildString();

  if(domainId() == 0) {
    std::cout << s;
  }
  m_log << s << std::endl;
}


template <MInt nDim>
void LbSolver<nDim>::preTimeStep() {
  TRACE();

  // update m_time for output
  updateTime();
}


//#################################################################################
// MOVING BOUNDARY PART STARTS HERE!!!!!!
//#################################################################################

template <MInt nDim>
void LbSolver<nDim>::preCoupleLs(std::vector<MInt>& maxGCellLevels) {
  TRACE();

  if(!grid().isActive()) {
    return;
  }

  RECORD_TIMER_START(m_t.findG0Cells);

  this->exchangeData(&a_levelSetFunctionMB(0, 0));

  MInt startSet = m_levelSetId; // if bodyToSetMode = 11 the 0th set is the collected set, which should be skiped

  resetActiveCellList(2);

  if(m_reducedComm) {
    resetComm();
    prepareCommunicationReduced();
  }

  RECORD_TIMER_START(m_t.findG0Candidates);
  findG0Candidates(maxGCellLevels); // findes all the G0Candidates of all Sets
  RECORD_TIMER_STOP(m_t.findG0Candidates);

  MBool gapClosure = false;
  m_G0Candidates.clear();

  MInt candidates = 0;
  for(MInt i = 0; i < a_noCells(); i++) {
    MBool isG0Candidate = false;
    for(MInt set = 0; set < (m_maxNoSets - startSet); set++) {
      if(a_isG0CandidateOfSet(i, set)) {
        isG0Candidate = true;
        break;
      }
    }
    m_G0CellMapping[i] = -1;
    if(isG0Candidate) {
      if(a_isHalo(i)) continue;
      m_G0Candidates.emplace_back();
      m_G0Candidates[candidates].cellId = i;
      m_G0CellMapping[i] = candidates;
      candidates++;
    }
  }

  MBoolScratchSpace isGapCell(a_noCells(), AT_, "isGapCell");
  isGapCell.fill(false);

  RECORD_TIMER_START(m_t.geomNodal);
  m_geometryIntersection->computeNodalValues(m_G0Candidates, &m_G0CellMapping[0], &a_levelSetFunctionMB(0, 0),
                                             &a_associatedBodyIds(0, 0), &isGapCell(0), gapClosure);
  RECORD_TIMER_STOP(m_t.geomNodal);

  if(noNeighborDomains() > 0) {
    RECORD_TIMER_START(m_t.geomExchange);

    // Save pointer to leafWindow/leafHaloCells per neighborDomain to get 'const MInt**' argument
    ScratchSpace<const MInt*> p_leafWindowCells(noNeighborDomains(), AT_, "p_leafWindowCells");
    ScratchSpace<const MInt*> p_leafHaloCells(noNeighborDomains(), AT_, "p_leafHaloCells");
    ScratchSpace<MInt> noLeafWindowCells(noNeighborDomains(), AT_, "noLeafWindowCells");
    ScratchSpace<MInt> noLeafHaloCells(noNeighborDomains(), AT_, "noLeafHaloCells");
    for(MInt d = 0; d < noNeighborDomains(); d++) {
      p_leafWindowCells[d] = &grid().leafWindowCell(d, 0);
      p_leafHaloCells[d] = &grid().leafHaloCell(d, 0);
      noLeafWindowCells[d] = grid().noLeafWindowCells(d);
      noLeafHaloCells[d] = grid().noLeafHaloCells(d);
    }

    m_geometryIntersection->exchangeNodalValues(p_leafWindowCells.data(), noLeafWindowCells.data(),
                                                p_leafHaloCells.data(), m_G0Candidates, &m_G0CellMapping[0]);
    RECORD_TIMER_STOP(m_t.geomExchange);
  }

  if(m_G0CellList != nullptr) {
    mDeallocate(m_G0CellList);
  }
  if(m_nodalGValues != nullptr) {
    mDeallocate(m_nodalGValues);
  }

  if(!m_G0Candidates.empty()) {
    MInt noNodalGValues = m_noLevelSetsUsedForMb * (IPOW3[nDim] - 1);
    m_noG0CandidatesTotal = (MInt)m_G0Candidates.size();
    mAlloc(m_G0CellList, m_noG0CandidatesTotal, "m_G0CellList", -1, AT_);
    mAlloc(m_nodalGValues, m_noG0CandidatesTotal, noNodalGValues, "m_nodalGValues", F0, AT_);
  } else {
    m_noG0CandidatesTotal = 0;
  }

  // reset Mapping
  std::fill_n(m_G0CellMapping, a_noCells(), -1);

  RECORD_TIMER_START(m_t.calcNodalValues);
  calcNodalLsValues();
  RECORD_TIMER_STOP(m_t.calcNodalValues);

  RECORD_TIMER_STOP(m_t.findG0Cells);
}

template <MInt nDim>
void LbSolver<nDim>::createBndryToBodyMapping(maia::coupling::Mapping& bndryToBodyMapping,
                                              maia::coupling::Mapping& bodyToBndryMapping) {
  // TODO labels:LB Use collector!
  // Create Mapping bndryCell -> bodyId
  bndryToBodyMapping.clear();
  bodyToBndryMapping.clear();
  for(MInt i = 0; i < m_currentNoG0Cells; i++) {
    const MInt bndryCellId = m_G0CellList[i];
    if(a_associatedBodyIds(bndryCellId, 0) >= 0) {
      bndryToBodyMapping.insert(i) = a_associatedBodyIds(bndryCellId, 0);
      bodyToBndryMapping.insert(a_associatedBodyIds(bndryCellId, 0)) = i;
    }
  }

  // Sanity checks
  for(MInt candidate = 0; candidate < m_currentNoG0Cells; candidate++) {
    const MInt cellId = m_G0CellList[candidate];
    ASSERT(cellId > -1, "G0CellList broken!");
    ASSERT(m_G0CellMapping[cellId] > -1, "G0CellMapping broken!");
    std::stringstream err;
    err << "G0CellMapping broken " << m_G0CellMapping[cellId] << " not eq to " << candidate;
    ASSERT(m_G0CellMapping[cellId] == candidate, err.str());
  }
}


/** \brief  Reset m_activeCellList as well as a_isActive
 *  \author Miro Gondrum
 *  \date   23.10.2020
 */
template <MInt nDim>
void LbSolver<nDim>::resetActiveCellList(MInt mode) {
  TRACE();
  //--store wasActive state-----------------------------------------------------
  for(MInt i = 0; i < a_noCells(); i++) {
    a_wasActive(i) = a_isActive(i);
  }
  //--initialization of a_isActive----------------------------------------------
  if(mode != 2) {
    setIsActiveDefaultStates();
  }
  setInActiveBndryCells();
  if(mode != 1) {
    setInActiveMBCells();
  }
  // TODO: Between the above and the below part, coupler might be called in
  // future. Then the upper and lower part has to be split into several
  // functions. (resetActiveStates and resetActiveCellList)
  //--reset activeCellList------------------------------------------------------
  fillActiveCellList();
}

// Reset external forces vector for each cell
template <MInt nDim>
void LbSolver<nDim>::resetExternalSources() {
  TRACE();
  maia::parallelFor(0, a_noCells(), [&](MInt cellId) {
    for(MInt dim = 0; dim < nDim; dim++) {
      a_externalForces(cellId, dim) = 0.0;
    }
  });
}

template <MInt nDim>
void LbSolver<nDim>::exchangeExternalSources() {
  TRACE();

  ASSERT(m_particleMomentumCoupling, "");

  if(noDomains() > 1 && m_particleMomentumCoupling) {
    maia::mpi::reverseExchangeAddData(grid().neighborDomains(), grid().haloCells(), grid().windowCells(), mpiComm(),
                                      &a_externalForces(0, 0), nDim);
  }
}

/** \brief  Set a_isActive to default values
 *  \author Miro Gondrum
 *  \date   23.10.2020
 *  a_isActive is 1 for all leaf and interfaceParents otherwise it is 0
 */
template <MInt nDim>
inline void LbSolver<nDim>::setIsActiveDefaultStates() {
  TRACE();
  const MInt noCells = m_cells.size();
  if(m_isRefined) {
    for(MInt i = 0; i < noCells; i++) {
      const MBool l_isActive = c_isLeafCell(i) || a_isInterfaceParent(i);
      a_isActive(i) = (l_isActive) ? 1 : 0;
    }
  } else {
    for(MInt i = 0; i < noCells; i++) {
      const MBool l_isActive = c_isLeafCell(i);
      a_isActive(i) = (l_isActive) ? 1 : 0;
    }
  }
  if(m_solidLayerExtension) {
    if(m_initialActiveCells != nullptr) {
      mDeallocate(m_initialActiveCells);
    }
    mAlloc(m_initialActiveCells, noCells, "m_initialActiveCells", 1, AT_);
    for(MInt i = 0; i < noCells; i++) {
      if(a_isActive(i)) {
        a_isActive(i) = !m_geometry->pointIsInside2(&a_coordinate(i, 0));
      }
      m_initialActiveCells[i] = a_isActive(i);
    }
  }
}

/** \brief  Set isActive = 0 for solid boundary cells
 *  \author Miro Gondrum
 *  \date   23.10.2020
 */
template <MInt nDim>
inline void LbSolver<nDim>::setInActiveBndryCells() {
  TRACE();
  // find ids of wall boundaries
  vector<MInt> wallIds;
  for(MInt i = 0; i < (MInt)(m_bndCnd->m_bndCndSegIds.size()); i++) {
    MBool is_periodic = false;
    if(m_bndCnd->m_noPeriodicSegments != 0)
      for(MInt j = 0; j < m_bndCnd->m_noPeriodicSegments; j++)
        if(m_bndCnd->m_bndCndSegIds[i] == m_bndCnd->m_periodicSegmentsIds[j]) {
          is_periodic = true;
          break;
        }
    MBool is_inout = false;
    for(MInt j = 0; j < m_bndCnd->m_noInOutSegments; j++)
      if(m_bndCnd->m_bndCndSegIds[i] == m_bndCnd->m_inOutSegmentsIds[j]) {
        is_inout = true;
        break;
      }
    if(!is_periodic & !is_inout) wallIds.push_back(i);
  }
  // set non-fluid cell to inactive
  for(auto wallId : wallIds) {
    for(MInt i = m_bndCnd->m_bndCndOffsets[wallId]; i < m_bndCnd->m_bndCndOffsets[wallId + 1]; i++) {
      if(m_bndCnd->m_bndCells[i].m_isFluid == false) {
        const MInt id = m_bndCnd->m_bndCells[i].m_cellId;
        a_isActive(id) = false;
      }
    }
  }
}


/** \brief  Set isActive = 0 for cells located within levelset solid to
 *  \author Miro Gondrum
 *  \date   23.10.2020
 */
template <MInt nDim>
inline void LbSolver<nDim>::setInActiveMBCells() {
  TRACE();
  if(m_maxNoSets == -1 || m_levelSetValues == nullptr) return;
  for(MInt i = 0; i < a_noCells(); i++) {
    if(m_solidLayerExtension) {
      if(!m_initialActiveCells[i]) continue;
    }
    for(MInt set = 0; set < m_maxNoSets; set++) {
      if(a_associatedBodyIds(i, set) < 0) continue;
      if(!c_isLeafCell(i)) continue;
      if(a_levelSetFunctionMB(i, set) < 0) {
        a_isActive(i) = 0;
      } else {
        a_isActive(i) = 1;
      }
    }
  }
}

/** \brief  Fills active cell's Ids into activeCellList
 *  \author Miro Gondrum
 *  \date   22.01.2021
 */
template <MInt nDim>
inline void LbSolver<nDim>::fillActiveCellList() {
  if(!isActive()) return;
  const MInt noLevels = maxLevel() - minLevel() + 1;

  // counting active cells - total and per level
  m_currentMaxNoCells = 0;
  ScratchSpace<MInt> noActiveCellsPerLevel(noLevels, AT_, "noActiveCellsPerLevel");
  noActiveCellsPerLevel.fill(0);
  for(MInt i = 0; i < a_noCells(); i++) {
    if(a_isActive(i)) {
      const MInt lvl = a_level(i) - minLevel();
      m_activeCellList[m_currentMaxNoCells] = i;
      noActiveCellsPerLevel[lvl]++;
      m_currentMaxNoCells++;
    }
  }
  // fill left entries with non-valid value
  for(MInt i = m_currentMaxNoCells; i < a_noCells(); i++) {
    m_activeCellList[i] = -1;
  }
  //  m_log << "    + number of active cells: " << m_currentMaxNoCells << endl << endl;
  // fill offset for better accessing later
  if(m_activeCellListLvlOffset != nullptr) {
    mDeallocate(m_activeCellListLvlOffset);
  }
  mAlloc(m_activeCellListLvlOffset, noLevels + 1, "m_activeCellListLvlOffset", 0, AT_);
  m_activeCellListLvlOffset[0] = 0;
  for(MInt i = 1; i < noLevels + 1; i++) {
    m_activeCellListLvlOffset[i] = m_activeCellListLvlOffset[i - 1] + noActiveCellsPerLevel[i - 1];
  }
  // sort active cells per level TODO labels:LB,toenhance sorting algorithm might be worth improving (HEAP-sort?)
  {
    MInt level = minLevel();
    MInt i = 0;
    for(MInt l = 0; l < noLevels - 1; l++) {
      MInt j = i + 1;
      for(; i < m_activeCellListLvlOffset[l + 1]; i++) {
        if(a_level(m_activeCellList[i]) == level) {
          j = i + 1;
          continue; // already at correct position
        }
        for(; j < m_currentMaxNoCells; j++) {
          if(a_level(m_activeCellList[j]) == level) {
            std::swap(m_activeCellList[i], m_activeCellList[j]);
            break; // now with correct level, next position i
          }
        }
      }
      level++;
    }
  }
}

/** \brief Exchanges the information about beeing a G0 Candidate in a Lattice-Boltzmann-Level-Set-Coupling
 *
 * \author Moritz Waldmann
 * \date 18.08.2019
 *
 **/
template <MInt nDim>
void LbSolver<nDim>::findG0Candidates(std::vector<MInt>& maxGCellLevels) {
  TRACE();

  const MInt startSet = m_levelSetId;

  for(MInt set = startSet; set < m_maxNoSets; set++) {
    const MInt positionInArray = set - startSet;
    const MFloat maxDistanceG0 = 4.0 * c_cellLengthAtLevel(maxGCellLevels[set]);

    for(MInt i = 0; i < a_noCells(); i++) {
      a_isG0CandidateOfSet(i, positionInArray) = false;
      if(a_isHalo(i)) continue;
      if(!c_isLeafCell(i)) continue;
      if(a_associatedBodyIds(i, set) == -2) continue;
      // if a_associatedBodyId is set to -2 it means that the corresponding cell is
      // not cell belonging to the LS solver. Hence, no valid LS values and body Ids
      // can be given here. Cells which are located next to cells which are not
      // owned by the LS solver require special treatment. These cells should not
      // be G0 candidates because their nodel G-value cannot be calculated accurately.
      // Hence, these cells are filltered out.
      MBool hasNeighborInPureLbDomain = false;
      for(MInt n = 0; n < IPOW3[nDim] - 1; n++) {
        if(c_neighborId(i, n) < 0) {
          continue;
        } else {
          if(a_associatedBodyIds(c_neighborId(i, n), set) == -2) {
            hasNeighborInPureLbDomain = true;
            a_associatedBodyIds(i, set) = -1;
            break;
          }
        }
      }
      if(hasNeighborInPureLbDomain) {
        a_isG0CandidateOfSet(i, positionInArray) = false;
      } else {
        if(m_allowBndryAsG0) {
          if((fabs(a_levelSetFunctionMB(i, set)) < maxDistanceG0)) {
            a_isG0CandidateOfSet(i, positionInArray) = true;
          } else {
            a_isG0CandidateOfSet(i, positionInArray) = false;
          }
        } else {
          if((fabs(a_levelSetFunctionMB(i, set)) < maxDistanceG0) && !a_isBndryCell(i)) {
            a_isG0CandidateOfSet(i, positionInArray) = true;
          } else {
            a_isG0CandidateOfSet(i, positionInArray) = false;
          }
        }
      }
    }
  }

  if(noDomains() > 1) {
    exchangeG0Candidates();
  }
}


/** \brief Exchanges the information about beeing a G0 Candidate in a Lattice-Boltzmann-Level-Set-Coupling
 *
 * \author Moritz Waldmann
 * \date 18.08.2019
 *
 **/
template <MInt nDim>
void LbSolver<nDim>::exchangeG0Candidates() {
  TRACE();

  MInt startSet = m_levelSetId; // if bodyToSetMode = 11 the 0th set is the collected set, which should be skiped
  for(MInt set = 0; set < (m_maxNoSets - startSet); set++) {
    // 1. Gathered information (to which level set does the cell belong) of the window cells.
    MInt sendNeighOffset = 0;
    for(MInt n = 0; n < noNeighborDomains(); n++) {
      for(MInt j = 0; j < noWindowCells(n); j++) {
        m_sendBufferMB[sendNeighOffset + j] = a_isG0CandidateOfSet(windowCell(n, j), set);
      }
      sendNeighOffset += noWindowCells(n);
    }

    // 2. Send the gathered information.
    sendNeighOffset = 0;
    for(MInt n = 0; n < noNeighborDomains(); n++) {
      MInt bufsize = noWindowCells(n);
      MPI_Issend(&m_sendBufferMB[sendNeighOffset], bufsize, MPI_CXX_BOOL, neighborDomain(n), 0, mpiComm(),
                 &mpi_request[n], AT_, "m_sendBufferMB[sendNeighOffset]");
      sendNeighOffset += noWindowCells(n);
    }

    // 3. Receive data from neighbors.
    MPI_Status status;
    MInt recNeighOffset = 0;
    for(MInt n = 0; n < noNeighborDomains(); n++) {
      MInt bufsize = noHaloCells(n);
      MPI_Recv(&m_receiveBufferMB[recNeighOffset], bufsize, MPI_CXX_BOOL, neighborDomain(n), 0, mpiComm(), &status, AT_,
               "m_receiveBufferMB[recNeighOffset]");
      recNeighOffset += noHaloCells(n);
    }

    for(MInt n = 0; n < noNeighborDomains(); n++)
      MPI_Wait(&mpi_request[n], &status, AT_);

    // 4. scatter the received data
    recNeighOffset = 0;
    for(MInt n = 0; n < noNeighborDomains(); n++) {
      for(MInt j = 0; j < noHaloCells(n); j++) {
        a_isG0CandidateOfSet(haloCell(n, j), set) = m_receiveBufferMB[recNeighOffset + j];
      }
      recNeighOffset += noHaloCells(n);
    }
  }
}

/**
 * \brief sets up interpolation stencil for levelSet interpolation
 * mode = 0 -> trilinear interpolation
 * cellId: cell on the highes refinement level which contains the point

 * \author Claudia Guenther
 * \date 03/2012
 */
template <MInt nDim>
MInt LbSolver<nDim>::setUpLbInterpolationStencil(const MInt cellId, MInt* interpolationCells, MFloat* point) {
  TRACE();

  std::function<MBool(const MInt, const MInt)> alwaysTrue = [&](const MInt, const MInt) { return true; };

  return this->setUpInterpolationStencil(cellId, interpolationCells, point, alwaysTrue, false);
}

template <MInt nDim>
MFloat LbSolver<nDim>::interpolateFieldDataLb(MInt* interpolationCells, MFloat* point, MInt set,
                                              std::function<MFloat(MInt, MInt)> scalarField,
                                              std::function<MFloat(MInt, MInt)> coordinate) {
  TRACE();

  return this->template interpolateFieldData<true>(interpolationCells, point, set, scalarField, coordinate);
}

template <MInt nDim>
void LbSolver<nDim>::initializeMovingBoundaries() {
  TRACE();

  //########################################################################
  // start moving boundary properties
  //########################################################################

  /*! \page propertyPage1
    \section trackMovingBndry
    <code>MBool LbSolver::m_trackMovingBndry </code>\n
    default = <code>true</code>\n
    Also read in fvsolver.h and lssolver.cpp\n
    Triggers the displacement of bodies in the moving boundary solver using the G Field.\n
    Possible values are:
    <ul>
    <li>true: displace bodies</li>
    <li>false: do not displace bodies</li>
    </ul>
    Keywords: <i>MOVING BOUNDARY, BODY DISPLACEMENT</i>
  */
  m_trackMovingBndry = true;
  m_trackMovingBndry = Context::getSolverProperty<MBool>("trackMovingBndry", m_solverId, AT_, &m_trackMovingBndry);

  /*! \page propertyPage1
    \section trackMbStart
    <code>MInt LbSolver::m_trackMbStart </code>\n
    default = <code>numeric_limits<MInt>::max()</code>\n
    Also read in fvsolver.h and lssolver.cpp\n
    For time steps smaller than m_trackMbStart, the bodies are not displaced and the G Field is not updated\n
    Possible values are:
    <ul>
    <li>integer > 0 and < numeric_limits<MInt>::max()</li>
    </ul>
    Keywords: <i>MOVING BOUNDARY, BODY DISPLACEMENT</i>
  */
  m_trackMbStart = 0;
  m_trackMbStart = Context::getSolverProperty<MInt>("trackMbStart", m_solverId, AT_, &m_trackMbStart);

  /*! \page propertyPage1
    \section trackMbEnd
    <code>MInt LbSolver::m_trackMbEnd </code>\n
    default = <code>numeric_limits<MInt>::max()</code>\n
    Also read in fvsolver.h and lssolver.cpp\n
    For time steps larger than m_trackMbEnd, the bodies are not displaced and the G Field is not updated\n
    Possible values are:
    <ul>
    <li>integer > 0 and < numeric_limits<MInt>::max()</li>
    </ul>
    Keywords: <i>MOVING BOUNDARY, BODY DISPLACEMENT</i>
  */
  m_trackMbEnd = numeric_limits<MInt>::max();
  m_trackMbEnd = Context::getSolverProperty<MInt>("trackMbEnd", m_solverId, AT_, &m_trackMbEnd);

  if(m_sendBufferMB != nullptr) mDeallocate(m_sendBufferMB);
  if(m_receiveBufferMB != nullptr) mDeallocate(m_receiveBufferMB);
  if(noDomains() > 1 && grid().isActive()) {
    MInt sumwin = 0;
    MInt sumhalo = 0;
    for(MInt d = 0; d < noNeighborDomains(); d++) {
      sumwin += noWindowCells(d);
      sumhalo += noHaloCells(d);
    }
    mAlloc(m_sendBufferMB, sumwin, "m_sendBufferMB", 0, AT_);
    mAlloc(m_receiveBufferMB, sumhalo, "m_receiveBufferMB", 0, AT_);
  }
  //########################################################################
  // end moving boundary properties
  //########################################################################

  //########################################################################
  // start updating ActiveCellList for Moving Boundaries
  //########################################################################

  // body properties
  if(m_associatedBodyIds != nullptr) mDeallocate(m_associatedBodyIds);
  mAlloc(m_associatedBodyIds, m_noLevelSetsUsedForMb * a_noCells(), "m_associatedBodyIds", -1, AT_);
  if(m_levelSetValues != nullptr) mDeallocate(m_levelSetValues);
  mAlloc(m_levelSetValues, m_noLevelSetsUsedForMb * a_noCells(), "m_levelSetValues", F0, AT_);

  MInt startSet = m_levelSetId; // if bodyToSetMode = 11 the 0th set is the collected set, which should be skiped
  MInt arrayLength = m_maxNoSets - startSet;
  if(m_isG0CandidateOfSet != nullptr) mDeallocate(m_isG0CandidateOfSet);
  mAlloc(m_isG0CandidateOfSet, a_noCells(), arrayLength, "m_isG0CandidateOfSet", false, AT_);

  //########################################################################
  // end updating ActiveCellList for Moving Boundaries
  //########################################################################

  if(m_G0CellMapping != nullptr) mDeallocate(m_G0CellMapping);
  mAlloc(m_G0CellMapping, a_noCells(), "m_G0CellMapping", -1, AT_);
  if(m_innerBandWidth != nullptr) mDeallocate(m_innerBandWidth);
  mAlloc(m_innerBandWidth, grid().maxRefinementLevel(), "m_innerBandWidth", F0, AT_);
  if(m_bandWidth != nullptr) mDeallocate(m_bandWidth);
  mAlloc(m_bandWidth, grid().maxRefinementLevel(), "m_bandWidth", 0, AT_);
  if(m_outerBandWidth != nullptr) mDeallocate(m_outerBandWidth);
  mAlloc(m_outerBandWidth, grid().maxRefinementLevel(), "m_outerBandWidth", F0, AT_);

  /*! \page propertyPage1
  \section mbBandwidth
  <code>MFloat LbSolver::distFac </code>\n
  default = {18.0, 9.0}\n
  Sets the distance factor which is used to calculate the inner (distFac[0]) and outer (distFac[1]) bandwidth.\n
  Possible values are:
  <ul>
    <li>Any floating point values.</li>
  </ul>
  Keywords: <i>MOVING BOUNDARY</i>
*/
  MFloat distFac[2] = {18.0, 9.0};
  for(MInt i = 0; i < 2; i++) {
    distFac[i] = Context::getSolverProperty<MFloat>("mbBandWidth", m_solverId, AT_, &distFac[i], i);
  }

  m_outerBandWidth[grid().maxRefinementLevel() - 1] = distFac[0] * c_cellLengthAtLevel(grid().maxRefinementLevel());
  m_bandWidth[grid().maxRefinementLevel() - 1] = distFac[0];
  for(MInt i = grid().maxRefinementLevel() - 2; i >= 0; i--) {
    m_outerBandWidth[i] = m_outerBandWidth[i + 1] + (distFac[1] * c_cellLengthAtLevel(i + 1));
    m_bandWidth[i] = (m_bandWidth[i + 1] / 2) + 1 + distFac[1];
  }
  for(MInt i = 0; i < grid().maxRefinementLevel(); i++) {
    m_innerBandWidth[i] = -m_outerBandWidth[i];
    m_log << "bandwidth level " << i << ": " << m_innerBandWidth[i] << " " << m_outerBandWidth[i] << endl;
  }

  /*! \page propertyPage1
    \section refineDiagonals
    <code>MBool LbSolver::m_refineDiagonals </code>\n
    default = <code>true</code>\n \n
    Determines whether the diagonal cells for the interface sensor should be refined as well!
    <ul>
    <li>true</li>
    <li>false</li>
    </ul>
    Keywords: <i>SENSOR, ADAPTATION</i>
  */
  m_refineDiagonals = true;
  m_refineDiagonals = Context::getSolverProperty<MBool>("refineDiagonals", m_solverId, AT_, &m_refineDiagonals);

  if(noNeighborDomains() > 0 && grid().isActive()) {
    grid().updateLeafCellExchange();
  }
}

/// \brief Reset solver data before new data is set during load balancing
///
/// \author Julian Vorspohl <j.vorspohl@aia.rwth-aachen.de>
template <MInt nDim>
void LbSolver<nDim>::balancePre() {
  // Update the grid proxy for this solver
  grid().update();

  // Just reset parallelization information if solver is not active and cleanup containers
  if(!isActive()) {
    updateDomainInfo(-1, -1, MPI_COMM_NULL, AT_);
    return;
  }

  // Set new domain info for solver
  updateDomainInfo(grid().domainId(), grid().noDomains(), grid().mpiComm(), AT_);

  // Clear cell collector
  m_cells.clear();

  // Resize cell collector to internal cells
  m_cells.append(grid().noInternalCells());

  // Now we are for the cell collector to be filled
  m_setCellDataFinished = false;
}

/// \brief Reinitialize solver data during load balancing after data has been transferred and set
///
/// \author Julian Vorspohl <j.vorspohl@aia.rwth-aachen.de>
template <MInt nDim>
void LbSolver<nDim>::balancePost() {
  // Cell collector data has been set
  m_setCellDataFinished = true;

  // Nothing to do if solver is not active
  if(!isActive()) {
    return;
  }

  // append the halo-cells
  m_cells.append(c_noCells() - m_cells.size());

  updateCellCollectorFromGrid();

  grid().findEqualLevelNeighborsParDiagonal(false);

  resetCellLists();

  // Restart boundary condition object to update cell lists
  restartBndCnd();
  resetActiveCellList();

  resetComm();

  if(m_reducedComm) {
    prepareCommunicationReduced();
  } else {
    prepareCommunicationNormal();
  }

  // Do an exchange to have the right values in halo cells
  if(noNeighborDomains() > 0) {
    exchange(1);
    exchangeOldDistributions();
  }

  // If grid adaptation is used make sure that a grid file is written with the next restart file
  // since the cells will be sorted during balacing, while the old grid after adaptation will
  // contain the unsorted grid cells
  if(m_adaptation) {
    m_adaptationSinceLastRestart = true;
    m_adaptationSinceLastSolution = true;
  }
}

/**
 * \brief Return the number of LB load types
 * \author Miro Gondrum
 * \date  15.01.2021
 */
template <MInt nDim>
MInt LbSolver<nDim>::noLoadTypes() const {
  const MInt noLevels = maxLevel() - minLevel() + 1;
  const MInt noBndryTypes = 1;

  return noLevels + noBndryTypes;
}

/**
 * \brief Return the default weights for LB load types
 * \author Miro Gondrum
 * \date  15.01.2021
 */
template <MInt nDim>
void LbSolver<nDim>::getDefaultWeights(MFloat* weights, std::vector<MString>& names) const {
  TRACE();

  std::cerr << "WARNING: Using powers of 0.5 as weights. Could lead to large number of cells for cases with large "
               "refinement level range."
            << std::endl;

  MInt count = 0;

  const MInt noLevels = maxLevel() - minLevel() + 1;
  for(MInt level = 0; level < noLevels; level++) {
    const MInt lvlDiff = noLevels - level;
    weights[count] = FFPOW2(lvlDiff);
    std::stringstream ss;
    ss << "lb_leaf_cell_level_" << minLevel() + level;
    names[count] = ss.str();
    count++;
  }

  // Boundary cells
  weights[count] = 1.0;
  names[count] = "lb_bndry_cell";
}

/* \brief Return the load of a single cell
 * \author Miro Gondrum
 * \date  15.01.2021
 * \param[in] cellId  Requested grid cell id
 * \param[in] weights Computational weights for different simulation components
 * \return Cell load
 */
template <MInt nDim>
MFloat LbSolver<nDim>::getCellLoad(const MInt gridCellId, const MFloat* const weights) const {
  TRACE();
  ASSERT(isActive(), "solver is not active");

  // Convert to solver cell id and check
  const MInt cellId = grid().tree().grid2solver(gridCellId);
  if(cellId < 0) {
    return 0;
  }
  // sanity check
  if(cellId < 0 || cellId >= grid().noInternalCells()) {
    TERMM(1, "The given cell id is invalid.");
  }

  // lb specific part
  MFloat cellLoad = 0.0; // Default cell load
  if(a_isActive(cellId)) {
    const MInt relLevel = c_level(cellId) - minLevel();
    cellLoad += weights[relLevel];
  }

  const MInt noLevels = maxLevel() - minLevel() + 1;
  if(a_isBndryCell(cellId)) {
    cellLoad += weights[noLevels];
  }

  return cellLoad;
}

// Get total number of loadQuantities for this domain...
template <MInt nDim>
void LbSolver<nDim>::getLoadQuantities(MInt* const loadQuantities) const {
  TRACE();

  // Reset
  std::fill_n(&loadQuantities[0], noLoadTypes(), 0);

  // Nothing to do if solver is not active
  if(!isActive()) {
    return;
  }

  const MInt noLevels = maxLevel() - minLevel() + 1;

  for(MInt level = 0; level < noLevels; level++) {
    const MInt noActiveCellsOfLevel = m_activeCellListLvlOffset[level + 1] - m_activeCellListLvlOffset[level];
    loadQuantities[level] = noActiveCellsOfLevel;
  }

  MInt noBoundaryCells = 0.0;
  for(MInt i = 0; i < a_noCells(); i++) {
    noBoundaryCells += a_isBndryCell(i);
  }
  loadQuantities[noLevels] = noBoundaryCells;
}
