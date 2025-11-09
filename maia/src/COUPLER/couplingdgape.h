// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef COUPLINGDGAPE_H_
#define COUPLINGDGAPE_H_

#include "DG/dgcartesiangalerkinprojection.h"
#include "DG/dgcartesiansolver.h"
#include "UTIL/parallelfor.h"
#include "coupling.h"

/* idea
                      +----------+
          +---------->| Coupling |<---------+
          |           +----+-----+           |
          |                ^                 |
   +------+-----+    +-----+------+    +-----+------+
   | CouplingFv |    | CouplingDg |    | CouplingLb |
   +------+-----+    +-----+------+    +-----+------+
          ^                ^                 ^
          |                |                 |
          |         +------+--------+        |
          |    +--->| CouplingDgApe |<--+    |
          |    |    +---------------+   |    |
          |    |                        |    |
 +--------+----+-------+              +-+----+----------+
 | DgCcAcousticPerturb |              | CouplingLbDgApe |
 +---------------------+              +-----------------+
*/

template <MInt nDim_>
struct MeanVariables;
template <>
struct MeanVariables<2> {
  static constexpr const MInt nDim = 2;
  static constexpr const MInt noVorticities = 1;

  // Mean Lamb vector
  static constexpr const MInt LAMB0_X = 0;
  static constexpr const MInt LAMB0_Y = 1;
  static constexpr const MInt LAMB0[nDim] = {0, 1};

  // Mean velocities
  static constexpr const MInt U0 = 2;
  static constexpr const MInt V0 = 3;
  static constexpr const MInt UU0[nDim] = {2, 3};

  static constexpr const MInt RHO0 = 4;
  static constexpr const MInt P0 = 5;
  static constexpr const MInt C0 = 6;

  // Mean speed of sound
  static constexpr const MInt DC0_X = 7;
  static constexpr const MInt DC0_Y = 8;
  static constexpr const MInt DC0[nDim] = {7, 8};

  // Mean vorticity
  static constexpr const MInt VORT0_Z = 9;
  static constexpr const MInt VORT0[noVorticities] = {9};

  // Mean velocity gradient
  static constexpr const MInt DU_DX = 10;
  static constexpr const MInt DU_DY = 11;

  static constexpr const MInt DV_DX = 12;
  static constexpr const MInt DV_DY = 13;

  static constexpr const MInt GRADU[nDim * nDim] = {10, 11, 12, 13};

  // du/dx, dv/dy, dw/dz for the divergence
  // NOTE: be careful when accessing data since entries of DU are not consecutive in memory if the
  // full mean velocity gradient is loaded
  static constexpr const MInt DU[nDim] = {10, 13};

  // Mean gradient of rho
  static constexpr const MInt DRHO_DX = 14;
  static constexpr const MInt DRHO_DY = 15;
  static constexpr const MInt DRHO[nDim] = {14, 15};

  // Mean gradient of p
  static constexpr const MInt DP_DX = 16;
  static constexpr const MInt DP_DY = 17;
  static constexpr const MInt DP[nDim] = {16, 17};

  // Mean gradient of rho*div(u)
  static constexpr const MInt RHODIVU_X = 18;
  static constexpr const MInt RHODIVU_Y = 19;
  static constexpr const MInt RHODIVU[nDim] = {18, 19};

  // Mean gradient of u*grad(rho)
  static constexpr const MInt UGRADRHO_X = 20;
  static constexpr const MInt UGRADRHO_Y = 21;
  static constexpr const MInt UGRADRHO[nDim] = {20, 21};

  // Mean of (gradient of p divided by rho)
  static constexpr const MInt GRADPRHO_X = 22;
  static constexpr const MInt GRADPRHO_Y = 23;
  static constexpr const MInt GRADPRHO[nDim] = {22, 23};

  // Sum of products of velocity components with corresponding velocity gradients
  static constexpr const MInt UGRADU_X = 24;
  static constexpr const MInt UGRADU_Y = 25;
  static constexpr const MInt UGRADU[nDim] = {24, 25};

  // Note: If the totalNoMeanVars value is changed here, it must also be changed
  // for s_totalNoMeanVars
  static constexpr const MInt totalNoMeanVars = 26;
  static const std::array<MString, totalNoMeanVars> names;
  static const std::array<MString, 6> nodeVarNames;
};
template <>
struct MeanVariables<3> {
  static constexpr const MInt nDim = 3;
  static constexpr const MInt noVorticities = 3;

  // Mean lamb vector
  static constexpr const MInt LAMB0_X = 0;
  static constexpr const MInt LAMB0_Y = 1;
  static constexpr const MInt LAMB0_Z = 2;
  static constexpr const MInt LAMB0[nDim] = {0, 1, 2};
  // Mean velocities
  static constexpr const MInt U0 = 3;
  static constexpr const MInt V0 = 4;
  static constexpr const MInt W0 = 5;
  static constexpr const MInt UU0[nDim] = {3, 4, 5};

  static constexpr const MInt RHO0 = 6;
  static constexpr const MInt P0 = 7;
  static constexpr const MInt C0 = 8;

  // Mean speed of sound
  static constexpr const MInt DC0_X = 9;
  static constexpr const MInt DC0_Y = 10;
  static constexpr const MInt DC0_Z = 11;
  static constexpr const MInt DC0[nDim] = {9, 10, 11};

  // Mean vorticities
  static constexpr const MInt VORT0_X = 12;
  static constexpr const MInt VORT0_Y = 13;
  static constexpr const MInt VORT0_Z = 14;
  static constexpr const MInt VORT0[noVorticities] = {12, 13, 14};

  // Mean velocity gradient
  static constexpr const MInt DU_DX = 15;
  static constexpr const MInt DU_DY = 16;
  static constexpr const MInt DU_DZ = 17;

  static constexpr const MInt DV_DX = 18;
  static constexpr const MInt DV_DY = 19;
  static constexpr const MInt DV_DZ = 20;

  static constexpr const MInt DW_DX = 21;
  static constexpr const MInt DW_DY = 22;
  static constexpr const MInt DW_DZ = 23;

  static constexpr const MInt GRADU[nDim * nDim] = {15, 16, 17, 18, 19, 20, 21, 22, 23};

  // du/dx, du/dy, dw/dz for the divergence
  // NOTE: be careful when accessing data since entries of DU are not consecutive in memory if the
  // full mean velocity gradient is loaded
  static constexpr const MInt DU[nDim] = {15, 19, 23};

  // Mean gradient of rho
  static constexpr const MInt DRHO_DX = 24;
  static constexpr const MInt DRHO_DY = 25;
  static constexpr const MInt DRHO_DZ = 26;
  static constexpr const MInt DRHO[nDim] = {24, 25, 26};

  // Mean gradient of p
  static constexpr const MInt DP_DX = 27;
  static constexpr const MInt DP_DY = 28;
  static constexpr const MInt DP_DZ = 29;
  static constexpr const MInt DP[nDim] = {27, 28, 29};

  // Mean gradient of rho*div(u)
  static constexpr const MInt RHODIVU_X = 30;
  static constexpr const MInt RHODIVU_Y = 31;
  static constexpr const MInt RHODIVU_Z = 32;
  static constexpr const MInt RHODIVU[nDim] = {30, 31, 32};

  // Mean gradient of u*grad(rho)
  static constexpr const MInt UGRADRHO_X = 33;
  static constexpr const MInt UGRADRHO_Y = 34;
  static constexpr const MInt UGRADRHO_Z = 35;
  static constexpr const MInt UGRADRHO[nDim] = {33, 34, 35};

  // Mean of (gradient of p divided by rho)
  static constexpr const MInt GRADPRHO_X = 36;
  static constexpr const MInt GRADPRHO_Y = 37;
  static constexpr const MInt GRADPRHO_Z = 38;
  static constexpr const MInt GRADPRHO[nDim] = {36, 37, 38};

  // Sum of products of velocity components with corresponding velocity gradients
  static constexpr const MInt UGRADU_X = 39;
  static constexpr const MInt UGRADU_Y = 40;
  static constexpr const MInt UGRADU_Z = 41;
  static constexpr const MInt UGRADU[nDim] = {39, 40, 41};

  // Note: If the totalNoMeanVars value is changed here, it must also be changed
  // for s_totalNoMeanVars
  static constexpr const MInt totalNoMeanVars = 42;
  static const std::array<MString, totalNoMeanVars> names;
  static const std::array<MString, 8> nodeVarNames;
};

/** \brief  Intermediate class for coupling DG solver with APE sysEqn
 *  \author Miro Gondrum (refactored stuff of Ansgar Niemoeller)
 *  \date   22.02.2021
 */
template <MInt nDim, class CouplingDonor>
class CouplingDgApe : public CouplingDonor, public CouplingDg<nDim, DgSysEqnAcousticPerturb<nDim>> {
  // Typedefs
 public:
  using Base = Coupling;
  using Base::startLoadTimer;
  using Base::stopLoadTimer;

  using BaseDonor = CouplingDonor;
  using DonorSolverType = typename BaseDonor::solverType;

  using SysEqn = DgSysEqnAcousticPerturb<nDim>;
  using ProjectionType = DgGalerkinProjection<nDim>;
  using BaseDg = CouplingDg<nDim, DgSysEqnAcousticPerturb<nDim>>;
  using DgCartesianSolverType = typename BaseDg::solverType;

  using BaseDg::dgSolver;
  using BaseDg::elements;
  using BaseDg::externalSource;
  using BaseDg::maxPolyDeg;
  using BaseDg::minPolyDeg;
  using BaseDg::noElements;
  using BaseDg::returnIdleRecord;
  using BaseDg::returnLoadRecord;
  using BaseDg::saveNodalData;
  using BaseDg::solverId; // TODO labels:COUPLER,DG use couplingId instead of solverId to request coupling related
                          // properties?
  using BaseDg::sysEqn;
  using MV = MeanVariables<nDim>; ///< Hold indices for mean variables

  //--Ctor & Dtor
  CouplingDgApe(const MInt couplingId, DgCartesianSolverType* dg, DonorSolverType* ds);
  virtual ~CouplingDgApe() = default;

  //--Methods
  // Virtual methods to override
  void init() override {
    readProperties();
    initCoupler();
    initData();
  };
  void finalizeSubCoupleInit(MInt /*couplingStep*/) override final{};
  // void finalizeCouplerInit() override{}; // defined in child classes
  void preCouple(MInt /*recipestep*/) final;
  void subCouple(MInt /*recipeStep*/, MInt /*solverId*/, std::vector<MBool>& /*solverCompleted*/) override final{};
  void postCouple(MInt /*recipestep*/) override final{};
  void cleanUp() override final{};

  void getCouplingTimings(std::vector<std::pair<MString, MFloat>>& timings, const MBool allTimings) override {
    const MString namePrefix = "c" + std::to_string(this->couplerId()) + "_";

    timings.emplace_back(namePrefix + "loadCouplerApe", returnLoadRecord());
    timings.emplace_back(namePrefix + "idleCouplerApe", returnIdleRecord());

#ifdef MAIA_TIMER_FUNCTION
    timings.emplace_back(namePrefix + "preCouple", RETURN_TIMER_TIME(m_timers[Timers::PreCouple]));

    if(allTimings) {
      // Full set of timings
      timings.emplace_back(namePrefix + "storeSources", RETURN_TIMER_TIME(m_timers[Timers::StoreSources]));
      timings.emplace_back(namePrefix + "calcSources", RETURN_TIMER_TIME(m_timers[Timers::CalcSources]));
      timings.emplace_back(namePrefix + "projectSourceTerms", RETURN_TIMER_TIME(m_timers[Timers::ProjectSourceTerms]));

      timings.emplace_back(namePrefix + "applyStoredSources", RETURN_TIMER_TIME(m_timers[Timers::ApplyStoredSources]));
    } else {
      // Reduced/essential set of timings
      timings.emplace_back(namePrefix + "projectSourceTerms", RETURN_TIMER_TIME(m_timers[Timers::ProjectSourceTerms]));
    }
#endif
  };

  void getDomainDecompositionInformation(std::vector<std::pair<MString, MInt>>& domainInfo) override final {
    const MString namePrefix = "c" + std::to_string(this->couplerId()) + "_";

    const MInt noCoupledDgElements = m_calcSourceElements.size();
    const MInt noCoupledDonorCells = m_calcSourceDonorCells.size();

    domainInfo.emplace_back(namePrefix + "noCoupledDgElements", noCoupledDgElements);
    domainInfo.emplace_back(namePrefix + "noCoupledDonorCells", noCoupledDonorCells);
  };
  MInt noCouplingTimers(const MBool allTimings) const override {
#ifdef MAIA_TIMER_FUNCTION
    if(allTimings) {
      return 7;
    } else {
      return 4;
    }
#else
    return 2;
#endif
  };

  //
  constexpr static MInt noVelocities() { return nDim; }                 ///< Return number of velocity variables
  constexpr static MInt noVorticities() { return (nDim == 2) ? 1 : 3; } ///< Return number of vorticity variables
  MInt noMeanVars() const { return m_activeMeanVars.size(); }           ///< Return number of mean variables

 protected:
  //--Methods
  void calcInitialCondition(const MFloat time);
  void initCoupler();
  void initData();
  void initProjection();
  void initRestart(const MFloat time, const MFloat dt);
  void initSourceFilter();
  void initTimers();

  MBool loadCouplingData(const MString& filename, const MString& name_, const MInt stride, MFloat* data);
  void loadMeanQuantities(const MBool skipNodeVars = false);
  void projectToElement(const MInt elementId, const MInt noVars, const MFloat* data, const MFloat* defaultValues,
                        MFloat* target);
  void saveFilterValues();
  void saveFilterValuesDonor();

  // Source term calculation
  void applyStoredSources(const MFloat time);
  void storeSources(const MFloat time, const MInt timeStep);

  virtual void performUnitConversion(const MString& /*name*/, const MInt /*count*/, const MInt /*stride*/,
                                     MFloat* /*data*/){};
  virtual void calcSourceLambLinearized(const MFloat* const velocity, const MFloat* const vorticity,
                                        MFloat* sourceTerms) = 0;
  virtual void calcSourceLambNonlinear(const MFloat* const velocity, const MFloat* const vorticity,
                                       MFloat* const sourceTerms) = 0;
  virtual void calcSourceQmII(const MFloat* const velocity, MFloat* const sourceTerms) = 0;
  virtual void calcSourceQmIII(const MFloat* const velocity, MFloat* sourceTerms) = 0;
  virtual void calcSourceQe(const MFloat* const velocity, const MFloat time, MFloat* const sourceTerms) = 0;
  virtual void calcSourceQc(const MFloat* const velocity, MFloat* const sourceTerms, const MFloat time,
                            const MInt timeStep) = 0;

  void saveSourceTermsDonorGrid(const MInt timeStep, const MFloat* const data);
  void saveSourceTermsTargetGrid(const MInt timeStep);

  // Virtual methods to be overriden
  void readProperties() override final;
  void checkProperties() override final{};

  // Donor solver specific 'accessors'
  virtual DonorSolverType& donorSolver(const MInt solverId = 0) const = 0;
  virtual void getDonorVelocityAndVorticity(const std::vector<MInt>& donorCellIds, MFloatScratchSpace& p_velocity,
                                            MFloatScratchSpace& p_vorticity) = 0;

  //--Variables
  MBool m_isRestart = false;               ///< Store whether this is a restart (in case special treatment is necessary)
  MBool m_hasDgCartesianSolver = false;    ///< Store whether this domain has DG cells/elements
  MBool m_hasDonorCartesianSolver = false; ///< Store whether this domain has Cartesian donor solver cells

  MInt m_maxNoNodesXD = -1; ///< Maximum number of nodes of an element (corresponding to maxPolyDeg)

  // Time step calculation
  // MInt m_timeStepOffset = -1; ///< Time step offset between CAA simulation and LES simulation
  MFloat m_fixedTimeStep = -1.0; ///< Fixed time step to use
  MFloat m_previousTime =
      std::numeric_limits<MFloat>::quiet_NaN(); ///< Previous time for the calculation of time derivatives
  MInt m_previousTimeStep = -1; ///< Previous time step (to determine whether new sources need to be calculated)

  std::vector<MInt> m_activeMeanVars{};     ///< List of active mean variables for all active source terms
  MString m_meanDataFileName{};             ///< File name for mean quantities
  std::vector<MInt> m_activeSourceTerms{};  ///< List of active source terms
  MBool m_saveSourceTermsDonorGrid = false; ///< Store whether the sources on the donor grid should be saved as well
  MInt m_saveSourceTermsInterval = -1;      ///< Interval at which source term data should be saved to disk
  MBool m_calcProjectionError = false;      ///< Calculate the L2 error of the Galerkin projection

  // Donor projection information
  MInt m_noActiveDonorCells = -1; ///< Number of 'active' donor cells, i.e. those donor leaf cells on
                                  ///< which source  terms need to be computed
  MInt m_noDonorCells = -1;       ///< Total number of donor cells on this domain
  std::vector<std::vector<MInt>> m_elementDonorLevels = {}; ///< Donor cell levels relative to element
  std::vector<std::vector<MInt>> m_elementDonorMap = {};    ///< Mapping from donor cells to elements
  std::vector<std::vector<MInt>> m_elementDonorPos = {};    ///< Donor cell positions on the corresponding cell level

  // Projection
  std::vector<ProjectionType> m_projection{}; ///< Galerkin projection (access by polynomial degree)

  // Source term calculation
  MBool m_isFirstSourceCalculation = true;    ///< Store whether this is the first calculation of the source terms
  std::vector<MFloat> m_sourceTerms{};        ///< Local storage for source terms
  std::vector<MInt> m_calcSourceDonorCells{}; ///< List of all donor cell ids for which source terms need to be computed
  std::vector<MInt> m_calcSourceElements{};   ///< List of all element ids for which source terms need to be computed
  MFloat m_sourceFactor = 1.0;                ///< Factor by which the source terms are multiplied

  // Error checking
  MBool m_checkConservation = false;    ///< Check if each Galerkin projection is conservative
  MFloat m_maxConservationError = -1.0; ///< Maximum conservation error
  MFloat m_maxL2Error = -1.0;           ///< Maximum computed L2 error of the Galerkin projection

  // Mean variables
  std::vector<MFloat> m_meanVars{}; ///< Local storage for mean variables of the donor cells

  // Geometry
  std::vector<MInt> m_elementInsideBody{}; ///< Marker if elements are inside a geometric object

  // Filter
  // MBool m_applySourceTermFilter = true; ///< Store whether to apply filtering to calculated source terms
  Filter<nDim> m_sourceTermFilter;          ///< Auxiliary object that handles source filtering
  std::vector<MFloat> m_filter{};           ///< Local storage for source filter values
  std::vector<MFloat> m_filterDonor{};      ///< Local storage for source filter values on donor cells
  MInt m_noCutModesLowPass = 0;             ///< Number of modes to cut using the low-pass source term
  MBool m_useLowPassFilter = false;         ///< Switch low pass on and off to allow disabling e.g. for node vars
  std::vector<MFloatTensor> m_vdmLowPass{}; ///< Vandermonde matrix/matrices for the low-pass filter
  MBool m_applySourceFilterDonor = false;   ///< Apply source filter on donor cells or on DG elements after projection
  MBool m_saveSourceTermFilter = false; ///< Store whether filter values should be written to a file at initialization

  MBool m_projectionFilter = false; ///< Use spatial projection filter (i.e. exclude elements from the projection)
  std::vector<MFloat> m_projectionFilterBox{}; ///< Spatial projection filter box (excluding)

  //--Structs
  // Total number of mean variables (Note: MV::totalNoMeanVars can't be used).
  // For the list of mean variables, see class definition for MV below.
  static const MInt s_totalNoMeanVars = MV::totalNoMeanVars;
  // List of indices of all mean variables in m_activeMeanVars, i.e. for each
  // active mean variable X, its position in m_activeMeanVars is given by
  // m_meanVarsIndex[MV::X]
  std::array<MInt, s_totalNoMeanVars> m_meanVarsIndex{};
  // static const std::array<MString, s_totalNoMeanVars> s_meanVarNames; ///< Names for all possible mean variables in
  // MV
  struct ST;                                                                  ///< Hold indices for source terms
  static const std::array<MString, ST::totalNoSourceTerms> s_sourceTermNames; ///< Source term names corresponding to ST

  // Timing
  struct Timers {
    enum {
      // Timer group and main timer
      TimerGroup,
      Class,

      // Regular (method-specific timers)
      Constructor,
      Init,
      InitialCondition,
      LoadMeanQuantities,
      PreCouple,
      StoreSources,
      LoadInstantaneous,
      CalcSources,
      CalcSourceLamb,
      CalcSourceQmII,
      CalcSourceQmIII,
      CalcSourceQe,
      CalcSourceQc,
      ProjectSourceTerms,
      ApplyStoredSources,

      // Accumulated timers for certain subsystems
      Accumulated,
      LoadCouplingData,

      // Counter
      _count
    };
  };
  std::array<MInt, Timers::_count> m_timers;

 private:
  void neededMeanVarsForSourceTerm(const MInt sourceTerm, std::vector<MInt>& meanVars) const;
};

// Use to differentiate between different source terms
template <MInt nDim, class CouplingDonor>
struct CouplingDgApe<nDim, CouplingDonor>::ST {
  static constexpr const MInt Q_mI = 0;
  static constexpr const MInt Q_mI_linear = 1;
  static constexpr const MInt Q_mII = 2;
  static constexpr const MInt Q_mIII = 3;
  static constexpr const MInt Q_mIV = 4;
  static constexpr const MInt Q_c = 5;
  static constexpr const MInt Q_e = 6;

  static constexpr const MInt totalNoSourceTerms = 7;
};

// Source term names
template <MInt nDim, class CouplingDonor>
const std::array<MString, CouplingDgApe<nDim, CouplingDonor>::ST::totalNoSourceTerms>
    CouplingDgApe<nDim, CouplingDonor>::s_sourceTermNames = {
        {"q_mI", "q_mI_linear", "q_mII", "q_mIII", "q_mIV", "q_c", "q_e"}};

//--Ctor & Dtor-----------------------------------------------------------------

/// \brief Initialize timers and read properties in c'tor.
///
/// \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
/// \date 2015-10-01
///
/// \param[in] solver_ A reference to the DG solver in which this coupling class
///                   is used.
template <MInt nDim, class CouplingDonor>
CouplingDgApe<nDim, CouplingDonor>::CouplingDgApe(const MInt couplingId, DgCartesianSolverType* dg, DonorSolverType* ds)
  : Coupling(couplingId),
    BaseDonor(couplingId, ds),
    CouplingDg<nDim, SysEqn>(couplingId, dg),
    m_sourceTermFilter(dg->solverId()) {
  TRACE();

  // Ensure consistency in arguments
  static_assert(MV::totalNoMeanVars == s_totalNoMeanVars, "Total number of mean vars is inconsistent");

  initTimers();

  RECORD_TIMER_STOP(m_timers[Timers::Constructor]);
}

//--Methods (public)------------------------------------------------------------

/// \brief Calculate source terms and add to the external source terms of the DG solver.
template <MInt nDim, class CouplingDonor>
void CouplingDgApe<nDim, CouplingDonor>::preCouple(MInt) {
  TRACE();
  RECORD_TIMER_START(m_timers[Timers::PreCouple]);

  const MInt timeStep = (m_hasDgCartesianSolver) ? dgSolver().getCurrentTimeStep() : -1;
  const MFloat time = (m_hasDgCartesianSolver) ? dgSolver().time() : -1.0;

  if(m_hasDgCartesianSolver && timeStep == m_previousTimeStep) {
    TERMM(1, "Error: preCouple already called for this time step.");
  }

  // Calculate source terms
  storeSources(time, timeStep);
  // Add stored source terms to the external source buffer of the DG solver
  applyStoredSources(time);
  // Update time step
  m_previousTimeStep = timeStep;

  RECORD_TIMER_STOP(m_timers[Timers::PreCouple]);
}

//--Methods (protected)---------------------------------------------------------

/// \brief Apply initial conditions for this coupling condition.
///
/// \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
/// \date 2015-10-01
///
/// \param[in] time Current simulation time for time-dependent ICs
template <MInt nDim, class CouplingDonor>
void CouplingDgApe<nDim, CouplingDonor>::calcInitialCondition(const MFloat NotUsed(time)) {
  TRACE();
  RECORD_TIMER_START(m_timers[Timers::InitialCondition]);

  // Load mean quantities from disk into node variables
  loadMeanQuantities();

  RECORD_TIMER_STOP(m_timers[Timers::InitialCondition]);
}

/// \brief Initialize the coupling condition (data structures).
///
/// \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
/// \date 2015-10-01
///
/// This also allocates memory for local storage of mean variables and source terms.
template <MInt nDim, class CouplingDonor>
void CouplingDgApe<nDim, CouplingDonor>::initCoupler() {
  TRACE();
  RECORD_TIMER_START(m_timers[Timers::Init]);

  m_hasDgCartesianSolver = dgSolver().isActive();
  if(!m_hasDgCartesianSolver) {
    TERMM_IF_NOT_COND(noElements() == 0, "inactive DG solver should have no elements");
  }

  // Determine the total number of donor cells on this domain, i.e., total number of donor cells on
  // this domain
  m_noDonorCells = donorSolver().noInternalCells();
  if(m_noDonorCells < 0) { // Reset if donorSolver is inactive
    m_noDonorCells = 0;
  }
  m_hasDonorCartesianSolver = (m_noDonorCells > 0);
  TERMM_IF_NOT_COND(m_hasDonorCartesianSolver == donorSolver().isActive(), "Donor solver status error");

  // Calculate maximum number of nodes
  m_maxNoNodesXD = (m_hasDgCartesianSolver) ? ipow(maxPolyDeg() + 1, nDim) : -1;

  // Initialize the Galerkin projection(s) and determine the information needed
  // to perform the Galerkin projection for all elements that have donor cells
  initProjection();

  // Allocate persistent storage for the source filter values. The maximum
  // polynomial degree is used here to be consisent with the memory allocation
  // in the solver.
  m_filter.resize(std::max(noElements() * m_maxNoNodesXD, 0));
  // Storage for donor cell filter
  m_filterDonor.resize(m_noDonorCells);

  // Allocate storage for element ids for which coupling source terms need to be
  // computed
  m_calcSourceElements.reserve(noElements());

  // Calculate source filter values and store element ids on which source terms
  // need to be computed, i.e. elements that have nonzero filter value(s) AND
  // also have donor cell(s)
  initSourceFilter();

  // Save filter values to file
  if(m_saveSourceTermFilter) {
    saveFilterValues();
    saveFilterValuesDonor();
  }

  // Allocate persistent storage for the source terms of all elements for which
  // coupling source terms need to be computed. The maximum polynomial degree is
  // used here to be consistent with the memory allocation in the solver.
  m_sourceTerms.resize(std::max(static_cast<MInt>(m_calcSourceElements.size()) * m_maxNoNodesXD * SysEqn::noVars(), 0));

  // Determine the set of donor cells on which source terms need to be computed
  std::set<MInt> calcSourceDonorCells;
  for(MInt elementId : m_calcSourceElements) {
    calcSourceDonorCells.insert(m_elementDonorMap[elementId].begin(), m_elementDonorMap[elementId].end());
  }
  // Number of 'active' donor cells
  m_noActiveDonorCells = calcSourceDonorCells.size();
  if(!m_hasDgCartesianSolver) {
    TERMM_IF_COND(m_noActiveDonorCells > 0, "Inactive DG solver cannot have active donor cells.");
  }

  // Store 'active' donor cell ids in m_calcSourceDonorCells
  m_calcSourceDonorCells.resize(m_noActiveDonorCells);
  m_calcSourceDonorCells.assign(calcSourceDonorCells.begin(), calcSourceDonorCells.end());

  // Allocate persistent storage for the mean variables of the donor cells
  m_meanVars.resize(m_noActiveDonorCells * noMeanVars());

  if(m_noCutModesLowPass > 0 && m_hasDgCartesianSolver) {
    ////////////////////////////////////////////////////////////////////////////
    // Temporarily create interpolation objects
    ////////////////////////////////////////////////////////////////////////////

    // Check for SBP Mode
    const MInt defaultSbpMode = 0;
    const MInt sbpMode = Context::getSolverProperty("sbpMode", solverId(), AT_, &defaultSbpMode);

    // Determine SBP Operator
    const MString defaultSbpOperator = "";
    const MString sbpOperator =
        Context::getSolverProperty<MString>("sbpOperator", solverId(), AT_, &defaultSbpOperator);

    // Determine initial number of nodes
    MInt initNoNodes1D = -1;
    initNoNodes1D = Context::getSolverProperty<MInt>("initNoNodes", solverId(), AT_, &initNoNodes1D);

    // Determine DG integration method (same as in dgsolver.cpp)
    const MString defaultDgIntegrationMethod = "DG_INTEGRATE_GAUSS";
    const MInt dgIntegrationMethod = string2enum(
        Context::getSolverProperty<MString>("dgIntegrationMethod", solverId(), AT_, &defaultDgIntegrationMethod));

    // Determine DG polynomial type (same as in dgsolver.cpp)
    const MString defaultDgPolynomialType = "DG_POLY_LEGENDRE";
    const MInt dgPolynomialType =
        string2enum(Context::getSolverProperty<MString>("dgPolynomialType", solverId(), AT_, &defaultDgPolynomialType));

    // Convert integers to enums
    auto polyType = static_cast<DgPolynomialType>(dgPolynomialType);
    auto intMethod = static_cast<DgIntegrationMethod>(dgIntegrationMethod);

    // Create & init interpolation objects
    std::vector<DgInterpolation> interpolation(maxPolyDeg() + 1);
    for(MInt i = std::max(0, minPolyDeg() - m_noCutModesLowPass); i < maxPolyDeg() + 1; i++) {
      const MInt noNodes1D = sbpMode ? initNoNodes1D : (i + 1);
      interpolation[i].init(i, polyType, noNodes1D, intMethod, sbpMode, sbpOperator);
    }

    ////////////////////////////////////////////////////////////////////////////
    // Create & store Vandermonde matrices
    ////////////////////////////////////////////////////////////////////////////

    // Create
    m_vdmLowPass.clear();
    m_vdmLowPass.resize(maxPolyDeg() + 1 - m_noCutModesLowPass);

    // Create matrices
    for(MInt polyDeg = std::max(0, minPolyDeg() - m_noCutModesLowPass);
        polyDeg < maxPolyDeg() + 1 - m_noCutModesLowPass;
        polyDeg++) {
      const MInt noNodesIn = polyDeg + 1;
      const MInt noNodesOut = polyDeg + 1 + m_noCutModesLowPass;
      m_vdmLowPass[polyDeg].resize(noNodesOut, noNodesIn);
      const DgInterpolation& interpIn = interpolation[polyDeg];
      const DgInterpolation& interpOut = interpolation[polyDeg + m_noCutModesLowPass];
      maia::dg::interpolation::calcPolynomialInterpolationMatrix(noNodesIn,
                                                                 &interpIn.m_nodes[0],
                                                                 noNodesOut,
                                                                 &interpOut.m_nodes[0],
                                                                 &interpIn.m_wBary[0],
                                                                 &m_vdmLowPass[polyDeg][0]);
    }
  } else {
    m_vdmLowPass.clear();
  }

  RECORD_TIMER_STOP(m_timers[Timers::Init]);
}

/// \brief Initialize the data (initial conditions or for a restart) of the coupling condition.
///
/// Note: not included anymore in initCoupler() since not required, e.g., during balancing.
template <MInt nDim, class CouplingDonor>
void CouplingDgApe<nDim, CouplingDonor>::initData() {
  TRACE();
  RECORD_TIMER_START(m_timers[Timers::Init]);

  if(dgSolver().m_restart) {
    // Load mean vars for source term calculation at a restart
    const MFloat time = dgSolver().time();
    // TODO labels:COUPLER,DG,toremove @ansgar_unified is dt already set?
    // const MFloat dt = dgSolver().m_dt;
    // initRestart(time, dt);
    // initRestart(time, -1.0);
    initRestart(time, m_fixedTimeStep); // NOTE: use prescribed fixed time step
  } else {
    // Load DG node variables and mean vars for source term calculation
    calcInitialCondition(-1.0);
    RECORD_TIMER_STOP(m_timers[Timers::Init]);

    // Apply the initial conditions of the solver (again) to allow overwriting node variables set
    // by the coupling
    if(m_hasDgCartesianSolver) {
      RECORD_TIMER_START(dgSolver().m_timers[maia::dg::Timers_::RunInit]);
      RECORD_TIMER_START(dgSolver().m_timers[maia::dg::Timers_::InitData]);
      dgSolver().initialCondition();
      RECORD_TIMER_STOP(dgSolver().m_timers[maia::dg::Timers_::InitData]);
      RECORD_TIMER_STOP(dgSolver().m_timers[maia::dg::Timers_::RunInit]);
    }

    RECORD_TIMER_START(m_timers[Timers::Init]);

    // Update node variables at surfaces in the solver and perform extension
    if(m_hasDgCartesianSolver && SysEqn::noNodeVars() > 0) {
      dgSolver().updateNodeVariables();
      dgSolver().extendNodeVariables();
    }
  }

  RECORD_TIMER_STOP(m_timers[Timers::Init]);
}

/// \brief Initialize the projection information for spatial coupling.
///
/// \author Ansgar Niemoeller (ansgar) <a.niemoeller@aia.rwth-aachen.de>
/// \date 2016-07-16
///
/// Initializes the required Galerkin projections for spatial coupling and
/// determines all information needed to perform the Galerkin projection for all
/// elements that have a mapping.
template <MInt nDim, class CouplingDonor>
void CouplingDgApe<nDim, CouplingDonor>::initProjection() {
  TRACE();

  if(!m_hasDgCartesianSolver) {
    // Clear member variables if DG solver is not active
    m_projection.clear();
    m_elementDonorMap.clear();
    m_elementDonorLevels.clear();
    m_elementDonorPos.clear();
    return;
  }

  // Create map from donor-cells to dg-elements
  std::vector<std::vector<MInt>> elementGridMap(noElements());

  for(MInt elementId = 0; elementId < noElements(); elementId++) {
    const MInt cellId = elements().cellId(elementId);
    const MInt gridCellId = dgSolver().grid().tree().solver2grid(cellId);

    elementGridMap[elementId].clear();
    // Get all leaf-cells of the donor-solver that are mapped to the global cell of the current
    // dg-element
    dgSolver().grid().raw().createLeafCellMapping(donorSolver().solverId(), gridCellId, elementGridMap[elementId]);
  }

  // Initialize the required Galerkin projections needed for spatial coupling
  {
    m_projection.clear();
    m_projection.reserve(maxPolyDeg() + 1);

    // Determine the number of refinement levels of the donor grid
    MInt noLevels = 1;
    if(m_hasDonorCartesianSolver) {
      noLevels = donorSolver().grid().maxLevel() - donorSolver().grid().minLevel() + 1;
    }

    m_log << "Initialize Galerkin projections for " << noLevels << " levels." << std::endl;

    // Dummy projections such that m_projection can be accessed by polyDeg
    for(MInt i = 0; i < std::max(0, minPolyDeg() - m_noCutModesLowPass); i++) {
      m_projection.push_back(ProjectionType());
    }

    // Initialize projection for each used polynomial degree
    for(MInt i = std::max(0, minPolyDeg() - m_noCutModesLowPass); i <= maxPolyDeg(); i++) {
      m_log << "Initialize Galerkin projection (p = " << i << ")... ";
      m_projection.push_back(ProjectionType(i, noLevels, dgSolver().grid().lengthLevel0(), dgSolver().solverId()));
      m_log << "done" << std::endl;
    }
  }

  m_elementDonorMap.resize(noElements());
  m_elementDonorLevels.resize(noElements());
  m_elementDonorPos.resize(noElements());

  // Determine information needed to apply the Galerkin projection
  for(MInt elementId = 0; elementId < noElements(); elementId++) {
    const MInt polyDeg = elements().polyDeg(elementId);

    // Clear any previous projection information for this element
    m_elementDonorMap[elementId].clear();
    m_elementDonorLevels[elementId].clear();
    m_elementDonorPos[elementId].clear();

    // Check if element has mapped donor cells and skip iteration otherwise
    if(elementGridMap[elementId].empty()) {
      continue;
    }

    const MInt noMappedLeafCells = elementGridMap[elementId].size();
    MIntScratchSpace donorCellLevels(noMappedLeafCells, AT_, "donorCellLevels");
    MFloatScratchSpace donorCellCoordinates(noMappedLeafCells, nDim, AT_, "donorCellCoordinates");

    // Information on target element
    const MInt targetCellId = elements().cellId(elementId);
    const MInt targetLevel = dgSolver().grid().tree().level(targetCellId);
    const MFloat halfCellLength = 0.5 * dgSolver().grid().cellLengthAtLevel(targetLevel);
    const MFloat* targetElementCenter = &dgSolver().grid().tree().coordinate(targetCellId, 0);

    // Check if element is inside projection filter box (excluding filter)
    if(m_projectionFilter) {
      MBool hasProjection = false;
      for(MInt dimId = 0; dimId < nDim; dimId++) {
        const MFloat coord = targetElementCenter[dimId];
        if(coord + halfCellLength < m_projectionFilterBox[dimId]
           || coord - halfCellLength > m_projectionFilterBox[nDim + dimId]) {
          hasProjection = true;
          break;
        }
      }
      if(!hasProjection) {
        continue;
      }
    }

    MInt noDonorCells = 0;
    // Loop over all mapped donor leaf-cells
    for(MInt i = 0; i < noMappedLeafCells; i++) {
      const MInt donorGridCellId = elementGridMap[elementId][i];
      const MInt donorCellId = donorSolver().grid().tree().grid2solver(donorGridCellId);

      if(donorCellId < 0) {
        TERMM(1, "Error in mapped donor leaf cells.");
      }

      // Store donor cell id in element donor map
      m_elementDonorMap[elementId].push_back(donorCellId);
      // Store the corresponding refinement level of the donor cell
      donorCellLevels[noDonorCells] = donorSolver().grid().tree().level(donorCellId);

      // Get original cell center coordinates
      // NOTE: For all donor cells (i.e. also for boundary and cut cells) the
      // original cell position is needed to establish the projection
      // information. There is no special treatment of cut cells at the
      // moment, for now they are treated as if they were normal cells.
      for(MInt d = 0; d < nDim; d++) {
        donorCellCoordinates(noDonorCells, d) = donorSolver().grid().tree().coordinate(donorCellId, d);
      }

      noDonorCells++;
    }

    if(noDonorCells == 0) {
      TERMM(1, "Donor cells for element (elementId=" + std::to_string(elementId)
                   + ") not found, but elementGridMap[elementId] "
                     "with elementGridMap[elementId].size()="
                   + std::to_string(elementGridMap[elementId].size()) + " contains mapped donor cells.");
    }

    // Allocate storage for the projection information of this element
    m_elementDonorLevels[elementId].resize(noDonorCells);
    m_elementDonorPos[elementId].resize(noDonorCells);

    // Check for one-to-multiple mapping and continue if so (there is no projection information to
    // compute)
    if(noDonorCells == 1) {
      const MInt donorCellId = m_elementDonorMap[elementId][0];
      const MInt donorLevel = donorSolver().grid().tree().level(donorCellId);
      if(targetLevel > donorLevel) {
        m_elementDonorLevels[elementId][0] = donorLevel - targetLevel;
        m_elementDonorPos[elementId][0] = -1;
        continue;
      }
    }

    // Compute and store the projection information for this element (this includes non-mapped
    // volumes)
    m_projection[polyDeg].calcProjectionInformation(noDonorCells, &donorCellLevels[0], &donorCellCoordinates[0],
                                                    targetLevel, targetElementCenter, m_elementDonorPos[elementId],
                                                    m_elementDonorLevels[elementId]);
  }
}

/// \brief Perform initializations for a restart.
///
/// \author Michael Schlottke-Lakemper (mic) <mic@aia.rwth-aachen.de>
///
/// \param[in] time Simulation time at time of restart.
/// \param[in] dt Previous time step size before restart.
template <MInt nDim, class CouplingDonor>
void CouplingDgApe<nDim, CouplingDonor>::initRestart(const MFloat time, const MFloat dt) {
  TRACE();

  // Store that this is a restart
  m_isRestart = true;

  // Calculate previous time for use in Qe source term
  m_previousTime = time - dt;

  RECORD_TIMER_START(m_timers[Timers::InitialCondition]);
  // Only load mean data for source term calculation since node variables are already read from the
  // nodevars-file
  loadMeanQuantities(true);
  RECORD_TIMER_STOP(m_timers[Timers::InitialCondition]);
}

/// \brief Calculate source filter values for all elements and store element ids
///        on which source terms need to be computed.
///
/// \author Ansgar Niemoeller (ansgar) <a.niemoeller@aia.rwth-aachen.de>
/// \date 2016-07-04
template <MInt nDim, class CouplingDonor>
void CouplingDgApe<nDim, CouplingDonor>::initSourceFilter() {
  TRACE();

  // Make sure that the element list is empty
  m_calcSourceElements.clear();

  if(!m_hasDgCartesianSolver) {
    m_filter.clear();
    m_elementInsideBody.clear();
    return;
  }

  const MInt filterOffset = m_maxNoNodesXD;

  // Temporary storage before permanent allocation with exact size
  MBoolScratchSpace hasNonZeroFilter(noElements(), AT_, "hasNonZeroFilter");
  std::fill(hasNonZeroFilter.begin(), hasNonZeroFilter.end(), false);

  // Reset inside body marker for all elements
  m_elementInsideBody.resize(noElements());
  std::fill(m_elementInsideBody.begin(), m_elementInsideBody.end(), 0);

  MInt noInsideElements = 0;

  for(MInt elementId = 0; elementId < noElements(); elementId++) {
    const MInt noNodes1D = elements().polyDeg(elementId) + 1;
    const MInt noNodesXD = ipow(noNodes1D, nDim);

    // Pointer to coordinates of source term
    const MFloat* const coordinates = &elements().nodeCoordinates(elementId);

    // Pointer to filter values of current element
    MFloat* const elementFilter = &m_filter[elementId * filterOffset];

    // Calculate filter values for all node coordinates of the element
    m_sourceTermFilter.filter(coordinates, noNodesXD, elementFilter);

    // Set filter values to zero if element is excluded via projection filter box
    // @ansgar_mb TODO labels:COUPLER,DG make this a function?
    if(m_projectionFilter) {
      MBool hasProjection = false;
      for(MInt dimId = 0; dimId < nDim; dimId++) {
        const MFloat coord = coordinates[dimId];
        if(coord < m_projectionFilterBox[dimId] || coord > m_projectionFilterBox[nDim + dimId]) {
          hasProjection = true;
          break;
        }
      }
      if(!hasProjection) {
        std::fill_n(elementFilter, noNodesXD, 0.0);
      }
    }

    // Check if any filter value is nonzero (above epsilon), i.e. source terms need to be
    // computed for this element
    hasNonZeroFilter[elementId] =
        std::any_of(elementFilter, elementFilter + noNodesXD, [](MFloat filter) { return filter > MFloatEps; });

    // @ansgar TODO labels:COUPLER fix this, what about geometries that are not inside the source region -> mean
    // vars (see loadMeanQuantities)?
    if(hasNonZeroFilter[elementId] && m_elementDonorMap[elementId].empty()) {
      m_elementInsideBody[elementId] = 1;

      // Reset filter
      hasNonZeroFilter[elementId] = false;
      std::fill_n(elementFilter, noNodesXD, 0.0);

      noInsideElements++;
      std::cerr << "Element " << std::to_string(elementId)
                << " has nonzero filter value(s), but there are no donor "
                   "cells mapped to this element. It will be considered to be inside the geometry."
                << std::endl;
      TERMM(1, "TODO check! (element has nonzero filter but no donor cells)");
    }
  }

  MPI_Allreduce(MPI_IN_PLACE, &noInsideElements, 1, maia::type_traits<MInt>::mpiType(), MPI_SUM, dgSolver().mpiComm(),
                AT_, "MPI_IN_PLACE", "noInsideElements");
  if(noInsideElements > 0) {
    m_log << "WARNING: Number of DG elements considered to be inside the geometry: " << noInsideElements << std::endl;
    if(dgSolver().domainId() == 0) {
      std::cerr << "WARNING: Number of DG elements considered to be inside the geometry: " << noInsideElements
                << std::endl;
    }
  }

  // Permanently store elements with non-zero filter
  const MInt noNonZeroElements = std::count(hasNonZeroFilter.begin(), hasNonZeroFilter.end(), true);
  m_calcSourceElements.resize(noNonZeroElements);

  // Store element id if source terms need to be computed, i.e. there is at
  // least one nonzero filter value AND there are donor cells mapped to this
  // element.
  for(MInt elementId = 0, nonZeroElementId = 0; elementId < noElements(); elementId++) {
    if(hasNonZeroFilter[elementId]) {
      if(!m_elementDonorMap[elementId].empty()) {
        m_calcSourceElements[nonZeroElementId] = elementId;
        nonZeroElementId++;
      } else {
        TERMM(1, "Element " + std::to_string(elementId)
                     + " has nonzero filter value(s), but there are no donor "
                       "cells mapped to this element. Go fix the source term "
                       "filter in your property file.");
      }
    }
  }

  // Compute spatial filter values on donor cells
  if(m_hasDonorCartesianSolver) {
    for(MInt i = 0; i < m_noDonorCells; i++) {
      const MFloat* const coordinates = &donorSolver().a_coordinate(i, 0);
      MFloat* const cellFilter = &m_filterDonor[i];

      m_sourceTermFilter.filter(coordinates, 1, cellFilter);
    }
  }
}

/// \brief Create timers.
///
/// \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
/// \date 2015-10-01
template <MInt nDim, class CouplingDonor>
void CouplingDgApe<nDim, CouplingDonor>::initTimers() {
  TRACE();

  // Create timer group & timer for solver, and start the timer
  NEW_TIMER_GROUP_NOCREATE(m_timers[Timers::TimerGroup],
                           "CouplingDgApe (solverId = " + std::to_string(solverId()) + ")");
  NEW_TIMER_NOCREATE(m_timers[Timers::Class], "total object lifetime", m_timers[Timers::TimerGroup]);
  RECORD_TIMER_START(m_timers[Timers::Class]);

  // Create & start constructor timer
  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::Constructor], "Constructor", m_timers[Timers::Class]);
  RECORD_TIMER_START(m_timers[Timers::Constructor]);

  // Create regular solver-wide timers
  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::Init], "init", m_timers[Timers::Class]);
  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::InitialCondition], "calcInitialCondition", m_timers[Timers::Class]);
  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::LoadMeanQuantities], "loadMeanQuantities",
                         m_timers[Timers::InitialCondition]);

  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::PreCouple], "PreCouple", m_timers[Timers::Class]);
  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::StoreSources], "storeSources", m_timers[Timers::PreCouple]);
  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::LoadInstantaneous], "loadInstantaneousQuantities",
                         m_timers[Timers::StoreSources]);

  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::CalcSources], "calcSources", m_timers[Timers::StoreSources]);
  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::CalcSourceLamb], "calcSourceLamb", m_timers[Timers::CalcSources]);
  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::CalcSourceQmII], "calcSourceQmII", m_timers[Timers::CalcSources]);
  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::CalcSourceQmIII], "calcSourceQmIII", m_timers[Timers::CalcSources]);
  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::CalcSourceQe], "calcSourceQe", m_timers[Timers::CalcSources]);
  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::CalcSourceQc], "calcSourceQc", m_timers[Timers::CalcSources]);

  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::ProjectSourceTerms], "projectSourceTerms", m_timers[Timers::StoreSources]);
  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::ApplyStoredSources], "applyStoredSources", m_timers[Timers::PreCouple]);

  // Create accumulated timers that monitor selected subsystems
  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::Accumulated], "selected accumulated timers", m_timers[Timers::Class]);
  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::LoadCouplingData], "loadCouplingData", m_timers[Timers::Accumulated]);
}

/// \brief Auxiliary method to load coupling data (i.e. anything coming from
///        the LES) from a file.
///
/// \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
/// \date 2015-10-01
///
/// \param[in] filename File name to open.
/// \param[in] name_ Attribute 'name' of variable to read.
/// \param[in] stride Store data with given stride.
/// \param[out] data Pointer to memory where read data should be stored.
template <MInt nDim, class CouplingDonor>
MBool CouplingDgApe<nDim, CouplingDonor>::loadCouplingData(const MString& filename,
                                                           const MString& name_,
                                                           const MInt stride,
                                                           MFloat* data) {
  TRACE();
  RECORD_TIMER_START(m_timers[Timers::Accumulated]);
  RECORD_TIMER_START(m_timers[Timers::LoadCouplingData]);

  using namespace maia::parallel_io;
  ParallelIo file(filename, PIO_READ, dgSolver().grid().raw().mpiComm());

  MString datasetName;
  MString variableName;
  MBool foundVariable = false;
  // Get the names of all datasets in the file
  std::vector<MString> datasets = file.getDatasetNames();
  // Loop over all datasets and check if it has the right 'name' attribute
  for(auto&& dataset : datasets) {
    if(file.hasAttribute("name", dataset)) {
      file.getAttribute(&variableName, "name", dataset);
      if(variableName == name_) {
        datasetName = dataset;
        foundVariable = true;
        break;
      }
    }
  }
  // Read if variable was found in file
  if(foundVariable) {
    ParallelIo::size_type count = 0;
    ParallelIo::size_type start = 0;

    if(m_hasDonorCartesianSolver) {
      count = donorSolver().noInternalCells();
      start = donorSolver().domainOffset(donorSolver().domainId());
    }

    file.setOffset(count, start);
    file.readArray(data, datasetName, stride);

    const MInt noInternalCells = donorSolver().noInternalCells();
    performUnitConversion(variableName, noInternalCells, stride, data);
  }

  RECORD_TIMER_STOP(m_timers[Timers::LoadCouplingData]);
  RECORD_TIMER_STOP(m_timers[Timers::Accumulated]);
  return foundVariable;
}

/// \brief Load mean velocities from file and store them to the node variables.
///
/// \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
/// \date 2015-10-01
template <MInt nDim, class CouplingDonor>
void CouplingDgApe<nDim, CouplingDonor>::loadMeanQuantities(const MBool skipNodeVars) {
  TRACE();
  RECORD_TIMER_START(m_timers[Timers::LoadMeanQuantities]);

  // Get the default node variables for regions without LES
  std::array<MFloat, SysEqn::noNodeVars()> defaultNodeVars;
  sysEqn().getDefaultNodeVars(&defaultNodeVars[0]);

  // Get the default node variables for regions inside a LES geometry
  std::array<MFloat, SysEqn::noNodeVars()> defaultNodeVarsBody;
  sysEqn().getDefaultNodeVarsBody(&defaultNodeVarsBody[0]);

  // Skip loading of the mean file (e.g. when doing a scaling and there is no useful mean file to be loaded)
  MBool skipLoadMeanFile = false;
  skipLoadMeanFile = Context::getSolverProperty<MBool>("skipLoadMeanFile", solverId(), AT_, &skipLoadMeanFile);

  // Load mean data for node variables, apply the projection and store results
  // in nodeVars
  if(!skipNodeVars && !skipLoadMeanFile) {
    m_log << "Coupling condition: loading node vars from mean file." << std::endl;
    MFloatScratchSpace data(std::max(m_noDonorCells, 1), SysEqn::noNodeVars(), AT_, "data");
    // Load mean velocity data (data is stored with stride SysEqn::noNodeVars())
    std::vector<MInt> undefinedNodeVarIds; // nodeVarNames that are not available in mean file
    for(MInt varId = 0; varId < SysEqn::noNodeVars(); varId++) {
      if(!loadCouplingData(m_meanDataFileName, MV::nodeVarNames[varId], SysEqn::noNodeVars(), &data(0, varId))) {
        // Variable name was not found in mean file
        undefinedNodeVarIds.push_back(varId);
      }
    }

    // Sanity check whether all data is available
    if(!undefinedNodeVarIds.empty()) {
      std::stringstream ss;
      ss << "Warning: Following dataset/s with attribute 'name' have not been found in file '" << m_meanDataFileName
         << "': ";
      for(MInt i : undefinedNodeVarIds) {
        ss << MV::nodeVarNames[i] << ", ";
      }
      ss << std::endl << "Those dataset/s have to be specified by the initial condition.";
      // Note: Here is only a warning called, since it might make sense not
      // providing some of the data through the mean file. For example the speed
      // of sound is not provided by the LB solver but set by the initial
      // condition as constant.
      // TERMM(1, ss.str());
      m_log << ss.str() << std::endl;
    }

    // Disable low-pass filter for node vars
    const MBool previousLowPassFilterState = m_useLowPassFilter;
    m_useLowPassFilter = false;

#ifdef _OPENMP
#pragma omp parallel for
#endif
    // Project data to elements and store in nodeVars
    for(MInt elementId = 0; elementId < noElements(); elementId++) {
      const MInt noDonorCells = m_elementDonorMap[elementId].size();
      const MInt noNonMappedCells = m_elementDonorPos[elementId].size() - noDonorCells;

      // FIXME labels:COUPLER @ansgar this is not correct for a general case where there are bodies
      // in the donor-solver that are located outside the source region as DG-elements inside these
      // objects are not identified to be inside a body!
      if((m_elementInsideBody[elementId] == 0) && noNonMappedCells == 0) {
        projectToElement(elementId, SysEqn::noNodeVars(), &data[0], &defaultNodeVars[0],
                         &elements().nodeVars(elementId));
      } else {
        projectToElement(elementId, SysEqn::noNodeVars(), &data[0], &defaultNodeVarsBody[0],
                         &elements().nodeVars(elementId));
      }
    }

    // Enable low-pass filter again after interpolation of node vars
    m_useLowPassFilter = previousLowPassFilterState;
  } else {
    m_log << "Coupling condition: skipped loading node vars from mean file." << std::endl;
  }

  // Load mean data for source term calculation on donor grid
  if(!skipLoadMeanFile) {
    MFloatScratchSpace meanVars(std::max(m_noDonorCells, 1), noMeanVars(), AT_, "meanVars");
    m_log << "Coupling condition noMeanVars: " << noMeanVars() << std::endl;
    // Load mean variables data (stored with stride noMeanVars() in meanVars)
    std::vector<MInt> undefinedMeanVarIds; // meanVarNames that are not available in mean file
    for(MInt varId = 0; varId < noMeanVars(); varId++) {
      m_log << "Load coupling data: " << varId << " " << MV::names[m_activeMeanVars[varId]] << std::endl;
      if(!loadCouplingData(m_meanDataFileName, MV::names[m_activeMeanVars[varId]], noMeanVars(), &meanVars[varId])) {
        // Variable name was not found in mean file
        undefinedMeanVarIds.push_back(varId);
      }
    }

    // Sanity check whether all data is available
    if(!undefinedMeanVarIds.empty()) {
      std::stringstream ss;
      ss << "Warning: Following dataset/s with attribute 'name' have not been found in file '" << m_meanDataFileName
         << "': ";
      for(MInt i : undefinedMeanVarIds) {
        ss << MV::names[m_activeMeanVars[i]] << ", ";
      }
      ss << std::endl;
      TERMM(1, ss.str());
    }

    // Store the mean variables data of all donor leaf cells that are mapped to
    // elements that need to calculate coupling source terms in m_meanVars
    for(MInt donorIndex = 0; donorIndex < m_noActiveDonorCells; donorIndex++) {
      const MInt donorId = m_calcSourceDonorCells[donorIndex];
      std::copy_n(&meanVars(donorId, 0), noMeanVars(), &m_meanVars[donorIndex * noMeanVars()]);
    }
  } else {
    m_log << "Coupling condition: skip loading mean vars, set to zero instead." << std::endl;
    std::fill(m_meanVars.begin(), m_meanVars.end(), 0.0);
  }

  RECORD_TIMER_STOP(m_timers[Timers::LoadMeanQuantities]);
}

/// \brief Project the given data fields onto a single element
///
/// \author Ansgar Niemoeller (ansgar) <a.niemoeller@aia.rwth-aachen.de>
/// \date 2016-07-16
///
/// \param[in] elementId The id of the element to perform the projection onto.
/// \param[in] noVar Number of variables to project.
/// \param[in] data Donor field(s).
/// \param[in] defaultValues Default values to use if there is no donor cell.
/// \param[out] target Target field containing the projected variables.
///
/// The variables in data need to be stored as: [u1, v1, w1, u2, v2, w2, ...]
template <MInt nDim, class CouplingDonor>
void CouplingDgApe<nDim, CouplingDonor>::projectToElement(const MInt elementId, const MInt noVars, const MFloat* data,
                                                          const MFloat* defaultValues, MFloat* target) {
  TRACE();
  const MInt polyDeg = elements().polyDeg(elementId);
  const MInt noNodesXD = ipow(polyDeg + 1, nDim);
  const MInt noDonorCells = m_elementDonorMap[elementId].size();
  // Number of non-mapped cells/volumes
  const MInt noNonMappedCells = m_elementDonorPos[elementId].size() - noDonorCells;

  // Check if more than one donor cell
  if(noDonorCells > 1) {
    // Apply the Galerkin projection
    if(!m_useLowPassFilter || m_noCutModesLowPass == 0) {
      m_projection[polyDeg].apply(noDonorCells, noNonMappedCells, &m_elementDonorMap[elementId][0],
                                  &m_elementDonorLevels[elementId][0], &m_elementDonorPos[elementId][0], &data[0],
                                  defaultValues, noVars, target);
    } else {
      MFloatScratchSpace tmp(noNodesXD * noVars, AT_, "tmp");
      m_projection[polyDeg - m_noCutModesLowPass].apply(
          noDonorCells, noNonMappedCells, &m_elementDonorMap[elementId][0], &m_elementDonorLevels[elementId][0],
          &m_elementDonorPos[elementId][0], &data[0], defaultValues, noVars, &tmp[0]);
      maia::dg::interpolation::interpolateNodes<nDim>(&tmp[0],
                                                      &m_vdmLowPass[polyDeg - m_noCutModesLowPass][0],
                                                      polyDeg + 1 - m_noCutModesLowPass,
                                                      polyDeg + 1,
                                                      noVars,
                                                      target);
    }

    if(m_checkConservation || m_calcProjectionError) {
      if(noNonMappedCells > 0) {
        TERMM(1, "Fix conservation/projection error comp. for the case with non-mapped cells");
      }

      const MInt cellId = elements().cellId(elementId);
      const MFloat targetLength = dgSolver().grid().cellLengthAtCell(cellId);
      MFloatScratchSpace errors(noVars, AT_, "errors");

      // Check if the projection is conservative
      if(m_checkConservation) {
        m_projection[polyDeg].calcConservationError(noDonorCells, &m_elementDonorMap[elementId][0],
                                                    &m_elementDonorLevels[elementId][0], &data[0], noVars, targetLength,
                                                    target, &errors[0]);

        const MFloat maxError = *(std::max_element(errors.begin(), errors.end()));
        m_maxConservationError = std::max(m_maxConservationError, maxError);
      }

      // Calculate the projection error
      if(m_calcProjectionError) {
        m_projection[polyDeg].calcL2Error(noDonorCells, &m_elementDonorMap[elementId][0],
                                          &m_elementDonorLevels[elementId][0], &m_elementDonorPos[elementId][0],
                                          &data[0], noVars, targetLength, target, &errors[0]);

        // Calculate maximum error
        const MFloat maxError = *(std::max_element(errors.begin(), errors.end()));
        m_maxL2Error = std::max(m_maxL2Error, maxError);
      }
    }
  } else if(noDonorCells == 1) {
    // One-to-one mapping
    MFloatTensor targetF(target, noNodesXD, noVars);

    // Store donor cell variables in all element nodes
    const MInt donorCellId = m_elementDonorMap[elementId][0];
    const MInt offset = donorCellId * noVars;
    for(MInt i = 0; i < noNodesXD; i++) {
      for(MInt v = 0; v < noVars; v++) {
        targetF(i, v) = data[offset + v];
      }
    }
  } else {
    // No mapped cell: use default value(s)
    MFloatTensor targetF(target, noNodesXD, noVars);
    for(MInt i = 0; i < noNodesXD; i++) {
      for(MInt v = 0; v < noVars; v++) {
        targetF(i, v) = defaultValues[v];
      }
    }
  }
}

/// \brief Save filter values to file.
///
/// \author Ansgar Niemoeller (ansgar) <a.niemoeller@aia.rwth-aachen.de>
/// \date 2016-12-03
template <MInt nDim, class CouplingDonor>
void CouplingDgApe<nDim, CouplingDonor>::saveFilterValues() {
  TRACE();

  m_log << "Saving file with filter values ... ";

  if(!m_hasDgCartesianSolver) {
    return;
  }

  const std::vector<MString> filterVarName{"sourcefilter:" + m_sourceTermFilter.name()};
  saveNodalData("sourcefilter", 1, filterVarName, &m_filter[0]);

  m_log << "done" << std::endl;
}

/// \brief Save the filter values on the donor cells
template <MInt nDim, class CouplingDonor>
void CouplingDgApe<nDim, CouplingDonor>::saveFilterValuesDonor() {
  TRACE();
  m_log << "Saving file with filter values on donor cells... ";

  if(!m_hasDonorCartesianSolver) {
    return;
  }

  std::stringstream fileName;
  fileName << donorSolver().outputDir() << "sourcefilterDonor" << ParallelIo::fileExt();

  using namespace maia::parallel_io;
  ParallelIo parallelIo(fileName.str(), PIO_REPLACE, donorSolver().mpiComm());

  const MInt localNoCells = m_noDonorCells;
  // Determine offset and global number of cells
  ParallelIo::size_type offset = 0;
  ParallelIo::size_type globalNoCells = 0;
  parallelIo.calcOffset(localNoCells, &offset, &globalNoCells, donorSolver().mpiComm());

  // Set attributes
  parallelIo.setAttribute(globalNoCells, "noCells");
  parallelIo.setAttribute(donorSolver().grid().gridInputFileName(), "gridFile");
  parallelIo.setAttribute(donorSolver().solverId(), "solverId");

  const MString name_ = "sourcefilter:" + m_sourceTermFilter.name();
  parallelIo.defineArray(PIO_FLOAT, name_, globalNoCells);
  parallelIo.setOffset(localNoCells, offset);
  parallelIo.writeArray(&m_filterDonor[0], name_);

  m_log << "done" << std::endl;
}

/// \brief Apply the stored source terms to the time derivative.
///
/// \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
/// \date 2015-10-01
template <MInt nDim, class CouplingDonor>
void CouplingDgApe<nDim, CouplingDonor>::applyStoredSources(const MFloat NotUsed(time)) {
  TRACE();
  RECORD_TIMER_START(m_timers[Timers::ApplyStoredSources]);

  const MInt dataSize = m_maxNoNodesXD * SysEqn::noVars();

  // Note: moved this into the DG solver applyExternalSource to remove the dependency
  // on the simulation time in this function, this way the external sources can be computed once in
  // the preCouple step before starting the first RK stage!
  // Source term ramp up factor (linear in time)
  /* const MFloat rampUpFactor */
  /*     = (m_useSourceRampUp) */
  /*           ? maia::filter::slope::linear(startTime(), */
  /*                                        startTime() + m_sourceRampUpTime, time) */
  /*           : 1.0; */

  const MInt noSourceElements = m_calcSourceElements.size();
  // Add source terms from persistent storage to right-hand side of all
  // elements with both nonzero filter values and donor cells.
  maia::parallelFor(0, noSourceElements, [&](MInt elementIndex) {
    const MInt elementId = m_calcSourceElements[elementIndex];
    const MInt noNodes1D = elements().polyDeg(elementId) + 1;
    const MInt noNodes1D3 = (nDim == 3) ? noNodes1D : 1;
    MFloatTensor r(externalSource() + elementId * dataSize, noNodes1D, noNodes1D, noNodes1D3, SysEqn::noVars());
    MFloatTensor s(&m_sourceTerms[0] + elementIndex * dataSize, noNodes1D, noNodes1D, noNodes1D3, SysEqn::noVars());

    // Add source terms to all equations
    for(MInt i = 0; i < noNodes1D; i++) {
      for(MInt j = 0; j < noNodes1D; j++) {
        for(MInt k = 0; k < noNodes1D3; k++) {
          for(MInt l = 0; l < SysEqn::noVars(); l++) {
            r(i, j, k, l) += m_sourceFactor * s(i, j, k, l);
          }
        }
      }
    }
  });

  RECORD_TIMER_STOP(m_timers[Timers::ApplyStoredSources]);
}

/// \brief Load coupling data from LES, then calculate, accumulate, project and
///        store coupling source terms.
///
/// \author Ansgar Niemoeller (ansgar) <a.niemoeller@aia.rwth-aachen.de>
/// \date 2016-08-18
///
/// \param[in] time Current simulation time.
/// \param[in] timeStep Current simulation time step.
///
/// The purpose of this method is to calculate the sources only once per time
/// step, and not once per Runge-Kutta step.
template <MInt nDim, class CouplingDonor>
void CouplingDgApe<nDim, CouplingDonor>::storeSources(const MFloat time, const MInt timeStep) {
  TRACE();
  RECORD_TIMER_START(m_timers[Timers::StoreSources]);

  // Storage for velocities and vorticities of the donor cells
  MFloatScratchSpace velocity(std::max(m_noDonorCells, 1), noVelocities(), AT_, "velocity");
  MFloatScratchSpace vorticity(std::max(m_noDonorCells, 1), noVorticities(), AT_, "vorticity");

  // Load all needed coupling data here
  // NOTE: for now (Q_mI and Q_mI_linear) the velocities and the vorticities are
  // needed for both source terms so there is no check yet which variables are
  // actually needed, but this will change in the future.
  RECORD_TIMER_START(m_timers[Timers::LoadInstantaneous]);
  getDonorVelocityAndVorticity(m_calcSourceDonorCells, velocity, vorticity);
  RECORD_TIMER_STOP(m_timers[Timers::LoadInstantaneous]);

  // Source terms on donor cells (i.e. all source terms summed up)
  MFloatScratchSpace sourceTerms(std::max(m_noDonorCells, 1), SysEqn::noVars(), AT_, "sourceTerms");
  std::fill(sourceTerms.begin(), sourceTerms.end(), 0.0);

  RECORD_TIMER_START(m_timers[Timers::CalcSources]);
  // Compute all source terms and accumulate to get the total source terms
  for(auto&& sourceTerm : m_activeSourceTerms) {
    switch(sourceTerm) {
      case ST::Q_mI:
        // Perturbed Lamb vector
        calcSourceLambNonlinear(&velocity[0], &vorticity[0], &sourceTerms(0, 0));
        break;
      case ST::Q_mI_linear:
        // Linearized Lamb vector
        calcSourceLambLinearized(&velocity[0], &vorticity[0], &sourceTerms(0, 0));
        break;
      case ST::Q_mII:
        // Pressure linearization
        calcSourceQmII(&velocity[0], &sourceTerms(0, 0));
        break;
      case ST::Q_mIII:
        calcSourceQmIII(&velocity[0], &sourceTerms(0, 0));
        break;
      case ST::Q_e:
        // Perturbed energy
        calcSourceQe(&velocity[0], time, &sourceTerms(0, 0));
        break;
      case ST::Q_c:
        calcSourceQc(&velocity[0], &sourceTerms(0, 0), time, timeStep);
        break;
      default:
        TERMM(1, "Source term '" + s_sourceTermNames[sourceTerm] + "' not implemented.");
        break;
    }
  }
  RECORD_TIMER_STOP(m_timers[Timers::CalcSources]);

  // Set flag to indicate that the source terms have already been calculated at least once
  m_isFirstSourceCalculation = false;

  // Store source terms on donor grid
  if(m_saveSourceTermsInterval > 0 && timeStep % m_saveSourceTermsInterval == 0 && m_saveSourceTermsDonorGrid) {
    stopLoadTimer(AT_);
    maia::dlb::g_dlbTimerController.disableAllDlbTimers();

    saveSourceTermsDonorGrid(timeStep, &sourceTerms[0]);

    maia::dlb::g_dlbTimerController.enableAllDlbTimers();
    startLoadTimer(AT_);
  }

  // TODO labels:COUPLER save filtered or unfiltered sources on donor grid?
  if(m_applySourceFilterDonor) {
    // Apply spatial filter to source terms on donor cells
    const MInt noVars = SysEqn::noVars();
    for(MInt i = 0; i < m_noDonorCells; i++) {
      const MFloat filter = m_filterDonor[i];
      for(MInt j = 0; j < noVars; j++) {
        sourceTerms[i * noVars + j] *= filter;
      }
    }
  }

  // Default source terms (NaN): only needed for the call to projectToElement()
  // but this should never be used/needed as an element without mapped LES cells
  // should not do anything here.
  /* const MFloat nan_value = std::numeric_limits<MFloat>::quiet_NaN(); */
  /* const std::array<MFloat, MAX_SPACE_DIMENSIONS + 1> defaultSourceTerms = { */
  /*     {nan_value, nan_value, nan_value, nan_value}}; */
  // Note: changed to 0 since non-mapped cells can now appear in the projection information that
  // will use these values
  const std::array<MFloat, MAX_SPACE_DIMENSIONS + 1> defaultSourceTerms = {{0.0, 0.0, 0.0, 0.0}};

  // Project accumulated source terms on elements and apply the spatial filter
  RECORD_TIMER_START(m_timers[Timers::ProjectSourceTerms]);
  const MInt noSourceElements = m_calcSourceElements.size();
  const MInt sourceOffset = m_maxNoNodesXD * SysEqn::noVars();

  maia::parallelFor(0, noSourceElements, [&](MInt elementIndex) {
    const MInt elementId = m_calcSourceElements[elementIndex];
    const MInt noNodes1D = elements().polyDeg(elementId) + 1;
    const MInt noNodesXD = ipow(noNodes1D, nDim);
    const MInt elementOffset = elementIndex * sourceOffset;

    MFloatTensor elementSourceTerms(&m_sourceTerms[elementOffset], noNodesXD, SysEqn::noVars());
    MFloatTensor filter(&m_filter[elementId * m_maxNoNodesXD], noNodesXD);

    // Project source terms on element
    projectToElement(elementId, SysEqn::noVars(), &sourceTerms[0], &defaultSourceTerms[0], &elementSourceTerms[0]);

    // Apply spatial filter on element (if not already applied on donor cells)
    if(!m_applySourceFilterDonor) {
      for(MInt i = 0; i < noNodesXD; i++) {
        for(MInt v = 0; v < SysEqn::noVars(); v++) {
          elementSourceTerms(i, v) *= filter(i);
        }
      }
    }
  });
  RECORD_TIMER_STOP(m_timers[Timers::ProjectSourceTerms]);

  // Store previous time for calculating dp/dt in the Qe source term
  m_previousTime = time;

  // Store source terms on target grid
  if(m_saveSourceTermsInterval > 0 && timeStep % m_saveSourceTermsInterval == 0) {
    stopLoadTimer(AT_);
    maia::dlb::g_dlbTimerController.disableAllDlbTimers();

    saveSourceTermsTargetGrid(timeStep);

    maia::dlb::g_dlbTimerController.enableAllDlbTimers();
    startLoadTimer(AT_);
  }

  RECORD_TIMER_STOP(m_timers[Timers::StoreSources]);
}

/// \brief Store the source terms on the donor grid.
///
/// \author Ansgar Niemoeller (ansgar) <a.niemoeller@aia.rwth-aachen.de>
/// \date 2016-12-03
///
/// \param[in] timeStep Current time step.
/// \param[in] data Pointer to source term data on donor grid.
template <MInt nDim, class CouplingDonor>
void CouplingDgApe<nDim, CouplingDonor>::saveSourceTermsDonorGrid(const MInt timeStep, const MFloat* const data) {
  TRACE();
  CHECK_TIMERS_IO();

  if(!m_hasDonorCartesianSolver) {
    return;
  }

  const MInt noVariables = SysEqn::noVars();

  std::stringstream fileName;
  fileName << donorSolver().outputDir() << "sourceTermsDonorGrid_" << std::setw(8) << std::setfill('0') << timeStep
           << ParallelIo::fileExt();

  // Define source term names
  std::vector<MString> sourceTermNames(noVariables);
  if constexpr(nDim == 2) {
    sourceTermNames = std::vector<MString>{"source_u", "source_v", "source_p"};
  } else {
    sourceTermNames = std::vector<MString>{"source_u", "source_v", "source_w", "source_p"};
  }

  // Create data out object to save data to disk
  using namespace maia::parallel_io;
  ParallelIo parallelIo(fileName.str(), PIO_REPLACE, donorSolver().mpiComm());

  const MInt localNoCells = m_noDonorCells;
  // Determine offset and global number of cells
  ParallelIo::size_type offset = 0;
  ParallelIo::size_type globalNoCells = 0;
  parallelIo.calcOffset(localNoCells, &offset, &globalNoCells, donorSolver().mpiComm());

  // Set attributes
  parallelIo.setAttribute(globalNoCells, "noCells");
  parallelIo.setAttribute(donorSolver().grid().gridInputFileName(), "gridFile");
  parallelIo.setAttribute(donorSolver().solverId(), "solverId");

  // Define arrays in file
  for(MInt i = 0; i < noVariables; i++) {
    const MString name_ = "variables" + std::to_string(i);
    parallelIo.defineArray(PIO_FLOAT, name_, globalNoCells);
    parallelIo.setAttribute(sourceTermNames[i], "name", name_);
  }

  parallelIo.setOffset(localNoCells, offset);

  // Write source term data to file
  for(MInt i = 0; i < noVariables; i++) {
    const MString name_ = "variables" + std::to_string(i);
    parallelIo.writeArray(&data[i], name_, noVariables);
  }
}

/// \brief Store the source terms on the target grid, i.e. after interpolation.
///
/// \author Ansgar Niemoeller (ansgar) <a.niemoeller@aia.rwth-aachen.de>
/// \date 2016-12-03
///
/// \param[in] timeStep Current time step.
template <MInt nDim, class CouplingDonor>
void CouplingDgApe<nDim, CouplingDonor>::saveSourceTermsTargetGrid(const MInt timeStep) {
  TRACE();
  CHECK_TIMERS_IO();

  const MInt noVariables = SysEqn::noVars();
  const MInt dataSize = m_maxNoNodesXD * noVariables;

  // Buffer for source terms on all elements to write to file
  ScratchSpace<MFloat> sourceTerms(noElements(), dataSize, AT_, "sourceTerms");

  // Copy source terms to buffer
  const MInt noSourceElements = m_calcSourceElements.size();
  for(MInt eId = 0; eId < noSourceElements; eId++) {
    const MInt elementId = m_calcSourceElements[eId];
    const MInt noNodesXD = ipow(elements().polyDeg(elementId) + 1, nDim);

    // Set source and destination pointers
    MFloat* const dest = &sourceTerms[elementId * dataSize];
    const MFloat* const src = &m_sourceTerms[eId * dataSize];

    // Copy all source terms of this element
    std::copy_n(src, noNodesXD * noVariables, dest);
  }

  std::stringstream fileNameBase;
  fileNameBase << "sourceTerms_" << std::setw(8) << std::setfill('0') << timeStep;

  // Define source term names
  std::vector<MString> sourceTermNames(noVariables);
  if constexpr(nDim == 2) {
    sourceTermNames = std::vector<MString>{"source_u", "source_v", "source_p"};
  } else {
    sourceTermNames = std::vector<MString>{"source_u", "source_v", "source_w", "source_p"};
  }

  // Save source terms to file
  saveNodalData(fileNameBase.str(), noVariables, sourceTermNames, &sourceTerms[0]);
}


/// \fn void CouplingDgApe<nDim>::readProperties()
/// \brief Read properties and store in member variables.
///
/// \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
/// \date 2015-10-02
template <MInt nDim, class CouplingDonor>
void CouplingDgApe<nDim, CouplingDonor>::readProperties() {
  TRACE();

  /*! \page propertyPage1
    \section meanDataFileName
    <code>MString CouplingDgApe::m_meanDataFileName </code>\n
    default = <code>none</code>\n \n
    Name of the file containing the mean velocities and vorticities.\n
    Keywords: <i>COUPLING, I/O, MEAN_DATA</i>
  */
  m_meanDataFileName = Context::getSolverProperty<MString>("meanDataFileName", solverId(), AT_);

  /*! \page propertyPage1
    \section sourceTerms
    <code>vector<MInt> CouplingDgApe::m_activeSourceTerms</code>\n
    default = <code>none</code>\n \n
    Names of the coupling source terms to use.\n
    The possible source term names are stored in s_sourceTermNames.\n
    Keywords: <i>COUPLING, PHYSICS, SOURCE_TERM</i>
  */
  // Get the number of source terms
  const MInt noSourceTerms = Context::propertyLength("sourceTerms", solverId());

  // Ensure that at least one source term is loaded
  if(noSourceTerms == 0) {
    TERMM(1, "No source term was specified.");
  }

  MBool q_mI_active = false;
  // Loop over all given source terms
  for(MInt i = 0; i < noSourceTerms; i++) {
    const MString sourceTermName = Context::getSolverProperty<MString>("sourceTerms", solverId(), AT_, i);

    // Find the source term name in the list of all source terms
    auto sourceTermIt = std::find(s_sourceTermNames.begin(), s_sourceTermNames.end(), sourceTermName);

    // Exit if source term name not found
    if(sourceTermIt == s_sourceTermNames.end()) {
      TERMM(1, "The given source term '" + sourceTermName + "' was not found.");
    }

    // Determine source term index
    const MInt sourceTerm = std::distance(s_sourceTermNames.begin(), sourceTermIt);

    // Check if this source term is already active
    if(std::find(m_activeSourceTerms.begin(), m_activeSourceTerms.end(), sourceTerm) != m_activeSourceTerms.end()) {
      TERMM(1, "Given source term '" + sourceTermName + "' already present. Check your property file.");
    }

    // Make sure that only one q_mI source term is active
    if(sourceTerm == ST::Q_mI || sourceTerm == ST::Q_mI_linear) {
      if(q_mI_active) {
        TERMM(1, "Source terms 'q_mI' and 'q_mI_linear' cannot be active at the same time.");
      }
      q_mI_active = true;
    }

    // Valid source term, store in list of active source terms
    m_activeSourceTerms.push_back(sourceTerm);
  }
  // Sort all active source terms by id
  std::sort(m_activeSourceTerms.begin(), m_activeSourceTerms.end());

  m_log << "Direct-hybrid CFD-CAA multi-solver coupling condition" << std::endl;
  m_log << "Activated coupling source terms:" << std::endl;
  for(auto&& sourceTerm : m_activeSourceTerms) {
    m_log << s_sourceTermNames[sourceTerm] << std::endl;
  }
  m_log << "Number of active source terms: " << noSourceTerms << std::endl;

  /*! \page propertyPage1
    \section saveSourceTermsInterval
    <code>MInt CouplingDgApe::m_saveSourceTermsInterval </code>\n
    default = <code>-1</code>\n \n
    Interval at which the source terms are stored to disk. If the value is < 1,
    no source term files are created.\n
    Possible values are:
    <ul>
      <li>any integer</li>
    </ul>
    Keywords: <i>COUPLING, SPATIAL_INTERPOLATION, SOURCE_TERMS, I/O</i>
  */
  m_saveSourceTermsInterval = -1;
  m_saveSourceTermsInterval =
      Context::getSolverProperty<MInt>("saveSourceTermsInterval", solverId(), AT_, &m_saveSourceTermsInterval);

  /*! \page propertyPage1
    \section saveSourceTermsDonorGrid
    <code>MBool CouplingDgApe::m_saveSourceTermsDonorGrid </code>\n
    default = <code>false</code>\n \n
    Enable/disable the saving of the source terms on the donor grid as well.
    Property is only relevant if "saveSourceTermsInterval" is > 0.\n
    Keywords: <i>COUPLING, SPATIAL_INTERPOLATION, FILTERING, I/O</i>
  */
  m_saveSourceTermsDonorGrid = false;
  m_saveSourceTermsDonorGrid =
      Context::getSolverProperty<MBool>("saveSourceTermsDonorGrid", solverId(), AT_, &m_saveSourceTermsDonorGrid);

  // Determine unique list of needed mean variables for all active source terms
  // and store in m_activeMeanVars
  std::set<MInt> neededMeanVars;
  for(auto&& sourceTerm : m_activeSourceTerms) {
    std::vector<MInt> sourceTermMeanVars;
    // Get the needed mean variables for the current source term
    neededMeanVarsForSourceTerm(sourceTerm, sourceTermMeanVars);
    // Insert into set of needed mean variables
    neededMeanVars.insert(sourceTermMeanVars.begin(), sourceTermMeanVars.end());
  }
  m_activeMeanVars.assign(neededMeanVars.begin(), neededMeanVars.end());

  // Store mapping from mean variable to corresponding index in m_activeMeanVars
  // Initialize with -1, since non-active mean variables should never be used
  std::fill_n(m_meanVarsIndex.begin(), s_totalNoMeanVars, -1);
  MInt meanVarPos = 0;
  // Loop over all active mean variables and store their position
  for(MInt meanVar : m_activeMeanVars) {
    m_meanVarsIndex[meanVar] = meanVarPos;
    meanVarPos++;
  }

  /*! \page propertyPage1
    \section fixedTimeStep
    <code>MInt CouplingDgApe::m_fixedTimeStep </code>\n
    default = <code>none</code>\n \n
    Use this fixed time step size when using offline coupling.\n
    Keywords: <i>COUPLING, TIME_STEP</i>
  */
  m_fixedTimeStep = Context::getSolverProperty<MFloat>("fixedTimeStep", solverId(), AT_, &m_fixedTimeStep);

  // Initialize source term filter
  m_sourceTermFilter.init();

  /*! \page propertyPage1
    \section applySourceFilterDonor
    <code>MBool CouplingDgApe::m_applySourceFilterDonor</code>\n
    default = <code>false</code>\n \n
    Enable/disable the application of the source term filter on the donor cells instead of on the DG
    elements after spatial projection of the source terms.\n
    Keywords: <i>COUPLING, SPATIAL_INTERPOLATION, FILTERING, I/O</i>
  */
  m_applySourceFilterDonor = false;
  m_applySourceFilterDonor =
      Context::getSolverProperty<MBool>("applySourceFilterDonor", solverId(), AT_, &m_applySourceFilterDonor);

  /*! \page propertyPage1
    \section saveSourceTermFilter
    <code>MBool CouplingDgApe::m_saveSourceTermFilter </code>\n
    default = <code>false</code>\n \n
    Enable/disable the saving of the spatial source term filter values to a file
    during the initialization phase.\n
    Keywords: <i>COUPLING, SPATIAL_INTERPOLATION, FILTERING, I/O</i>
  */
  m_saveSourceTermFilter = false;
  m_saveSourceTermFilter =
      Context::getSolverProperty<MBool>("saveSourceTermFilter", solverId(), AT_, &m_saveSourceTermFilter);

  /*! \page propertyPage1
    \section checkConservation
    <code>MBool CouplingDgApe::m_checkConservation </code>\n
    default = <code>false</code>\n \n
    Select if for each Galerkin projection the conservation error is calculated.
    If the error exceeds a certain threshold the simulation will abort.\n
    Keywords: <i>COUPLING, SPATIAL_INTERPOLATION, CONSERVATION_ERROR</i>
  */
  m_checkConservation = false;
  m_checkConservation = Context::getSolverProperty<MBool>("checkConservation", solverId(), AT_, &m_checkConservation);

  /*! \page propertyPage1
    \section calcProjectionError
    <code>MBool CouplingDgApe::m_calcProjectionError</code>\n
    default = <code>false</code>\n \n
    Select if for each Galerkin projection the interpolation error in the
    L2-norm is calculated. The maximum global projection error is written to the
    m_log at the end of the simulation.\n
    Keywords: <i>COUPLING, SPATIAL_INTERPOLATION, PROJECTION_ERROR</i>
  */
  m_calcProjectionError = false;
  m_calcProjectionError =
      Context::getSolverProperty<MBool>("calcProjectionError", solverId(), AT_, &m_calcProjectionError);

  /*! \page propertyPage1
    \section sourceFactor
    <code>MFloat CouplingDgApe::m_sourceFactor</code>\n
    default = <code>1.0</code>\n \n
    Set a factor by which all source terms are scaled.\n
    Keywords: <i>COUPLING, SOURCE_TERM, SCALING</i>
  */
  m_sourceFactor = 1.0;
  m_sourceFactor = Context::getSolverProperty<MFloat>("sourceFactor", solverId(), AT_, &m_sourceFactor);

  /*! \page propertyPage1
    \section noCutModesLowPass
    <code>MFloat CouplingDgApe::m_noCutModesLowPass</code>\n
    default = <code>0</code>\n \n
    Allows to specify a low-pass filter for the source terms.\n
    Keywords: <i>COUPLING, SOURCE_TERM, FILTER</i>
  */
  m_noCutModesLowPass = 0;
  m_noCutModesLowPass = Context::getSolverProperty<MInt>("noCutModesLowPass", solverId(), AT_, &m_noCutModesLowPass);

  /*! \page propertyPage1
    \section projectionFilter
    <code>MBool CouplingDgApe::m_projectionFilter</code>\n
    default = <code>false</code>\n \n
    Enables a projection filter box to exclude a region from using mean/node varibles from the mean
    file. Specify box via projectionFilterBox property.\n
    Keywords: <i>COUPLING, MEAN_VARIABLES, FILTER</i>
  */
  m_projectionFilter = false;
  m_projectionFilter = Context::getSolverProperty<MBool>("projectionFilter", solverId(), AT_, &m_projectionFilter);

  if(m_projectionFilter) {
    m_log << "Using projection filter box" << std::endl;
    m_projectionFilterBox.resize(2 * nDim);
    for(MInt i = 0; i < 2 * nDim; i++) {
      /*! \page propertyPage1
        \section projectionFilterBox
        <code>std::vector<MFloat> CouplingDgApe::m_projectionFilterBox</code>\n
        Projection filter box to exclude a region from using mean/node varibles from the mean
        file.\n
        Keywords: <i>COUPLING, MEAN_VARIABLES, FILTER</i>
      */
      m_projectionFilterBox[i] = Context::getSolverProperty<MFloat>("projectionFilterBox", solverId(), AT_, i);
      m_log << "projection filter box " << i << " " << m_projectionFilterBox[i] << std::endl;
    }
  }
}

//--Methods (private)-----------------------------------------------------------

/// \brief Return the needed mean variables for a given source term.
///
/// \author Ansgar Niemoeller (ansgar) <a.niemoeller@aia.rwth-aachen.de>
/// \date 2016-07-30
///
/// \param[in] sourceTerm The id of the requestest source term.
/// \param[out] meanVars Vector of needed mean variable ids.
template <MInt nDim, class CouplingDonor>
void CouplingDgApe<nDim, CouplingDonor>::neededMeanVarsForSourceTerm(const MInt sourceTerm,
                                                                     std::vector<MInt>& meanVars) const {
  TRACE();

  meanVars.clear();
  switch(sourceTerm) {
    case ST::Q_mI:
      // Mean Lamb vector
      for(MInt i = 0; i < nDim; i++) {
        meanVars.push_back(MV::LAMB0[i]);
      }
      break;
    case ST::Q_mI_linear:
      // Mean velocities
      for(MInt i = 0; i < nDim; i++) {
        meanVars.push_back(MV::UU0[i]);
      }
      // Mean vorticities
      for(MInt i = 0; i < noVorticities(); i++) {
        meanVars.push_back(MV::VORT0[i]);
      }
      break;
    case ST::Q_mII:
      // Mean density and pressure
      meanVars.push_back(MV::RHO0);
      meanVars.push_back(MV::P0);
      // Mean gradient of rho
      for(MInt i = 0; i < nDim; i++) {
        meanVars.push_back(MV::DRHO[i]);
      }
      // Mean gradient of p
      for(MInt i = 0; i < nDim; i++) {
        meanVars.push_back(MV::DP[i]);
      }
      // Mean of (gradient of p divided by rho)
      for(MInt i = 0; i < nDim; i++) {
        meanVars.push_back(MV::GRADPRHO[i]);
      }
      break;
    case ST::Q_mIII:
      // Mean velocities
      for(MInt i = 0; i < nDim; i++) {
        meanVars.push_back(MV::UU0[i]);
      }
      // Mean velocity gradients
      for(MInt i = 0; i < nDim * nDim; i++) {
        meanVars.push_back(MV::GRADU[i]);
      }
      // Sum of products of velocity components with corresponding velocity gradients
      for(MInt i = 0; i < nDim; i++) {
        meanVars.push_back(MV::UGRADU[i]);
      }
      break;
    case ST::Q_e:
      meanVars.push_back(MV::RHO0);
      meanVars.push_back(MV::P0);
      meanVars.push_back(MV::C0);
      // Mean velocities
      for(MInt i = 0; i < nDim; i++) {
        meanVars.push_back(MV::UU0[i]);
      }
      // Mean gradient of c
      for(MInt i = 0; i < nDim; i++) {
        meanVars.push_back(MV::DC0[i]);
      }
      // Components of divergence of u
      for(MInt i = 0; i < nDim; i++) {
        meanVars.push_back(MV::DU[i]);
      }
      // Mean gradient of rho
      for(MInt i = 0; i < nDim; i++) {
        meanVars.push_back(MV::DRHO[i]);
      }
      // Mean gradient of p
      for(MInt i = 0; i < nDim; i++) {
        meanVars.push_back(MV::DP[i]);
      }
      break;
    case ST::Q_c:
      // Mean density, pressure, and speed of sound
      meanVars.push_back(MV::RHO0);
      meanVars.push_back(MV::P0);
      meanVars.push_back(MV::C0);
      // Mean velocities
      for(MInt i = 0; i < nDim; i++) {
        meanVars.push_back(MV::UU0[i]);
      }
      // Components of divergence of u
      for(MInt i = 0; i < nDim; i++) {
        meanVars.push_back(MV::DU[i]);
      }
      // Mean gradient of rho
      for(MInt i = 0; i < nDim; i++) {
        meanVars.push_back(MV::DRHO[i]);
      }
      // Mean gradient of rho*div(u)
      for(MInt i = 0; i < nDim; i++) {
        meanVars.push_back(MV::RHODIVU[i]);
      }
      // Mean gradient of u*grad(rho)
      for(MInt i = 0; i < nDim; i++) {
        meanVars.push_back(MV::UGRADRHO[i]);
      }
      break;
    default:
      TERMM(1, "Source term '" + s_sourceTermNames[sourceTerm] + "' not implemented yet.");
      break;
  }
}

#endif // COUPLINGDGAPE_H_
