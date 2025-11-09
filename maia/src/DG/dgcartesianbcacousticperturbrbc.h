// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef DGBOUNDARYCONDITIONACOUSTICPERTURBRBC_H_
#define DGBOUNDARYCONDITIONACOUSTICPERTURBRBC_H_

#include "UTIL/timer.h"
#include "dgcartesianboundarycondition.h"

template <MInt nDim>
class DgBcAcousticPerturbRBC final : public DgBoundaryCondition<nDim, DgSysEqnAcousticPerturb<nDim>> {
  // Typedefs
 public:
  using Base = DgBoundaryCondition<nDim, DgSysEqnAcousticPerturb<nDim>>;
  using Base::applyJacobian;
  using Base::begin;
  using Base::calcRegularSurfaceFlux;
  using Base::calcSourceTerms;
  using Base::calcSurfaceIntegral;
  using Base::calcVolumeIntegral;
  using Base::dt;
  using Base::elements;
  using Base::end;
  using Base::flux;
  using Base::getElementByCellId;
  using Base::getHElementId;
  using Base::helements;
  using Base::id;
  using Base::integrationMethod;
  using Base::interpolation;
  using Base::isMpiSurface;
  using Base::isRestart;
  using Base::maxNoNodes1D;
  using Base::maxPolyDeg;
  using Base::needHElementForCell;
  using Base::resetBuffer;
  using Base::solver;
  using Base::subTimeStepRk;
  using Base::surfaces;
  using Base::sysEqn;
  using Base::timeIntegrationScheme;
  using typename Base::SolverType;
  using typename Base::SysEqn;

  // Methods
  DgBcAcousticPerturbRBC(SolverType& solver_, MInt bcId) : Base(solver_, bcId) {
    // Create timers
    // @ansgar TODO labels:DG,TIMERS this creates new timers every time DLB is performed and new RBCs are created
    initTimers();
  }
  ~DgBcAcousticPerturbRBC() {
    // Free allocated memory
    mDeallocate(m_rbSolution);
    // Stop the class timer
    RECORD_TIMER_STOP(m_timers[Timers::Class]);
  }

  MString name() const override { return "radiation boundary condition"; }

  void init() override;

  void initTimers();

  void apply(const MFloat time) override;

  void applyAtSurface(const MInt surfaceId, const MFloat NotUsed(time));

  void calcFlux(const MFloat* nodeVars, const MFloat* q, const MInt noNodes1D, MFloat* flux_) const;

  void calcRiemann(const MFloat* nodeVarsL,
                   const MFloat* nodeVarsR,
                   const MFloat* stateL,
                   const MFloat* stateR,
                   const MInt noNodes1D,
                   const MInt dirId,
                   MFloat* flux_) const;

  void calcSource(const MFloat* nodeVars, const MFloat* u, const MInt noNodes1D, const MFloat t, const MFloat* x,
                  MFloat* src) const;

  MInt noRestartVars() const override { return Base::SysEqn::noVars(); };
  MInt getLocalNoNodes() const override;
  MString restartVarName(const MInt id_) const override;
  void setRestartVariable(const MInt id_, const MFloat* const data) override;
  void getRestartVariable(const MInt id_, MFloat* const data) const override;

  MInt noBcElements() const override { return m_rbElements.size(); };
  MBool hasBcElement(const MInt elementId) const override;

  /// Dynamic load balancing
  // Number of variables in BC
  MInt noCellDataDlb() const override { return noRestartVars(); }

  // Data types of variables
  MInt cellDataTypeDlb(const MInt NotUsed(dataId)) const override {
    return MFLOAT; // Note: assumes there are not MInt to communicate
  }

  // Data size of variables
  MInt cellDataSizeDlb(const MInt dataId, const MInt cellId) const override;

  // Get BC solution data
  void getCellDataDlb(const MInt dataId, MFloat* const data) const override { getRestartVariable(dataId, data); }

  // Set BC solution data
  void setCellDataDlb(const MInt dataId, const MFloat* const data) override { setRestartVariable(dataId, data); }

 private:
  void calcFlux1D(const MFloat* nodeVars, const MFloat* q, const MInt noNodes1D, const MInt dirId, MFloat* flux_) const;

  void calcElementNodeVars();
  void calcSurfaceNodeVars();
  void calcNodeVars(const MFloat* pos, const MFloat* nodeVarsAPE, MFloat* nodeVars);

  // Member variables
 private:
  using ElementCollector = maia::dg::collector::ElementCollector<nDim, SysEqn>;
  using HElementCollector = maia::dg::collector::HElementCollector<nDim, SysEqn>;
  using SurfaceCollector = maia::dg::collector::SurfaceCollector<nDim, SysEqn>;

  // Element variables
  // Collector of rb-elements
  ElementCollector m_rbElements;

  // Store corresponding rb-element id for each boundary surface
  std::vector<MInt> m_rbElementId;
  // Store corresponding rb-surface id for each boundary surface
  std::vector<MInt> m_rbSurfaceId;

  // Store corresponding element id for each rb-element
  std::vector<MInt> m_elementId;
  // Store corresponding surface id for each rb-surface
  std::vector<MInt> m_surfaceId;

  // Surface variables
  // Collector of rb-surfaces
  SurfaceCollector m_rbSurfaces;

  // Store additional rb-boundary surface ids (with a different bcId) and their
  // corresponding rb-element
  std::vector<std::pair<MInt, MInt>> m_addBndrySrfcIds;

  // h-refinement
  // Collector of elements
  HElementCollector m_hRbElements;

  // Store solution of RBC system
  MFloat** m_rbSolution = nullptr;

  // Total number of values (DOF * noVars * noRbElements)
  MInt m_internalDataSize = -1;

  // Current Runge Kutta stage
  MInt m_rkStage = -1;

  // RBC reference position which should be close to the acoustic sources
  std::array<MFloat, MAX_SPACE_DIMENSIONS> m_rbcReferencePos{};

  // Store if RBC node variables are initialized and initial condition is set
  MBool m_isInitialized = false;

  // Number of node variables
  static const constexpr MInt s_noNodeVars = nDim + 2;

  // Timing
  struct Timers {
    enum {
      // Timer group and main timer
      TimerGroup,
      Class,

      // Regular (method-specific timers)
      Init,
      Apply,
      InitializeVars,
      CopySurfaceVars,
      Prolong,
      TimeDeriv,
      RungeKuttaStep,
      ApplyAtSurface,
      SetRestartVars,
      GetRestartVars,

      // Counter
      _count
    };
  };
  std::array<MInt, Timers::_count> m_timers;
};


/// \brief Initialize RBC.
///
/// \author Ansgar Niemoeller (ansgar) <a.niemoeller@aia.rwth-aachen.de>
/// \date 2015-12-13
///
/// 1. Create rb-elements for all elements along the boundary.
/// 2. Create needed h-rb-elements.
/// 3. Create rb-surfaces for all rb-elements.
/// 4. Initialize h-rb-elements and create missing h-refined surfaces.
template <MInt nDim>
void DgBcAcousticPerturbRBC<nDim>::init() {
  TRACE();
  RECORD_TIMER_START(m_timers[Timers::Init]);

  /*! \property
    \page propertyPageDG DG
    \section rbcReferencePos
    <code>std::array<MFloat, nDim> DgBcAcousticPerturbRBC::m_rbcReferencePos</code>\n
    default = <code>none</code>\n \n
    Set the reference position for the radiation boundary condition(s), which should be close to the
    acoustic sources (e.g. assume a point source in the acoustic source region).\n
    Possible values are:
    <ul>
      <li>point coordinates</li>
    </ul>
    Keywords: <i>DISCONTINUOUS_GALERKIN, BOUNDARY CONDITION</i>
  */
  for(MInt i = 0; i < nDim; i++) {
    m_rbcReferencePos[i] = Context::getSolverProperty<MFloat>("rbcReferencePos", solver().solverId(), AT_, i);
  }

  // Element mapping: elementId -> rbElementId
  std::map<MInt, MInt> elementMap;
  MInt noRbElements = 0;

  // Collect rb-elements to be created
  // Loop over boundary and check internal elements
  for(MInt srfcId = begin(); srfcId < end(); srfcId++) {
    const MInt internalSide = surfaces().internalSideId(srfcId);
    const MInt elementId = surfaces().nghbrElementIds(srfcId, internalSide);

    // Check if rb-element id already exists for current element
    const auto it = elementMap.find(elementId);
    if(it != elementMap.end()) {
      // rb-element already exists, store its id for current boundary surface
      m_rbElementId.push_back(it->second);
      continue;
    }

    // Determine new rb-element id
    const MInt rbElementId = noRbElements++;

    // Store mapping: elementId -> rbElementId
    elementMap[elementId] = rbElementId;

    // Store rb-element id for current boundary surface
    m_rbElementId.push_back(rbElementId);
  }

  MInt maxNoSurfaces = 2 * nDim * noRbElements;
  MInt noHRbElements = 0;
  ScratchSpace<MInt> hRbElementIds(noRbElements, AT_, "hRbElementIds");

  m_elementId.clear();
  // Store element ids in ascending order, create rb-elements in this order and
  // check for h-refinement
  for(const auto& it : elementMap) {
    const MInt elementId = it.first;
    m_elementId.push_back(elementId);

    // Check if h-rb-element is needed, store rb-element id
    const MInt cellId = elements().cellId(elementId);
    if(needHElementForCell(cellId)) {
      hRbElementIds[noHRbElements] = it.second;
      noHRbElements++;

      // Upper bound for additional surfaces needed by the h-refinement
      maxNoSurfaces += 2 * nDim * (2 * (nDim - 1) - 1);
    }
  }

  // Allocate rb-elements
  m_rbElements.maxPolyDeg(maxPolyDeg());
  m_rbElements.maxNoNodes1D(maxNoNodes1D());
  m_rbElements.noNodeVars(s_noNodeVars);
  m_rbElements.reset(noRbElements);

  // Allocate rb-surfaces
  m_rbSurfaces.maxPolyDeg(maxPolyDeg());
  m_rbSurfaces.maxNoNodes1D(maxNoNodes1D());
  m_rbSurfaces.noNodeVars(s_noNodeVars);
  m_rbSurfaces.reset(maxNoSurfaces);

  // Allocate h-rb-elements
  m_hRbElements.reset(noHRbElements);

  // Create rb-elements
  for(MInt rbId = 0; rbId < noRbElements; rbId++) {
    m_rbElements.append();

    // Copy member variables from elements to rb-elements
    m_rbElements.copy(elements(), m_elementId[rbId], rbId);
  }

  // Create h-rb-elements
  for(MInt hRbId = 0; hRbId < noHRbElements; hRbId++) {
    m_hRbElements.append();

    // Set corresponding rb-element id
    m_hRbElements.elementId(hRbId) = hRbElementIds[hRbId];
  }

  // Create rb-surfaces
  m_rbSurfaceId.assign(end() - begin() + 1, -1);
  std::map<MInt, MInt> srfcMap;
  const MInt noDirs = 2 * nDim;
  // Iterate over rb-elements
  for(MInt rbId = 0; rbId < m_rbElements.size(); rbId++) {
    const MInt eId = m_elementId[rbId]; // Element id

    // Loop over all directions
    for(MInt dir = 0; dir < noDirs; dir++) {
      const MInt sId = elements().surfaceIds(eId, dir);

      const MInt nghbrSideId = dir % 2;
      const MInt elementSideId = (nghbrSideId + 1) % 2;

      // Check if rb-surface is already created
      const auto srfcMapIt = srfcMap.find(sId);
      if(srfcMapIt == srfcMap.end()) { // rb-surface does not exist
        const MInt rbSurfaceId = m_rbSurfaces.size();
        m_rbSurfaces.append();

        // Store mapping: surfaceId -> rbSurfaceId
        srfcMap[sId] = rbSurfaceId;

        // Store surface id
        m_surfaceId.push_back(sId);

        // Store rb-surface ids only for boundary surfaces
        if(begin() <= sId && sId < end()) {
          m_rbSurfaceId[sId - begin()] = rbSurfaceId;
        }

        // Store rb-surface id in rb-element
        m_rbElements.surfaceIds(rbId, dir) = rbSurfaceId;

        // Copy variables of surface to rb-surface
        m_rbSurfaces.copy(surfaces(), sId, rbSurfaceId);

        // Store neighboring rb-element ids in rb-surface
        m_rbSurfaces.nghbrElementIds(rbSurfaceId, elementSideId) = rbId;
        m_rbSurfaces.nghbrElementIds(rbSurfaceId, nghbrSideId) = -1;

        const MInt internalSide = surfaces().internalSideId(sId);

        // Check for boundary surface which is not belonging to the RBC boundary
        // Note: Such boundary surfaces are assumed to belong also to a
        // radiation boundary, i.e. the RBC surface flux is computed with the
        // prolonged RBC solution as inner and outer state
        // TODO labels:DG other boundary conditions (e.g. solid wall)
        if(internalSide != -1 && (sId < begin() || sId >= end()) && !isMpiSurface(sId)) {
          // Store the rb-surface id and also the corresponding rb-element id
          // such that later the RBC-solution can also be prolonged to this
          // rb-surface
          m_addBndrySrfcIds.push_back(std::make_pair(rbSurfaceId, rbId));
        }
      } else { // rb-surface exists, store id in rb-element
        const MInt rbSurfaceId = srfcMapIt->second;
        m_rbElements.surfaceIds(rbId, dir) = rbSurfaceId;

        // Store neighboring rb-element id in rb-surface
        m_rbSurfaces.nghbrElementIds(rbSurfaceId, elementSideId) = rbId;
      }
    }
  }

  // Store h-refined surface ids in h-rb-elements, create missing surfaces
  const MInt noSurfs = 2 * (nDim - 1);
  for(MInt hRbId = 0; hRbId < noHRbElements; hRbId++) {
    const MInt rbId = m_hRbElements.elementId(hRbId);
    const MInt eId = m_elementId[rbId];
    const MInt hId = getHElementId(eId);

    for(MInt dir = 0; dir < 2 * nDim; dir++) {
      const MInt side = 1 - dir % 2;
      const MInt oppositeSide = dir % 2;

      for(MInt hSurfId = 0; hSurfId < noSurfs; hSurfId++) {
        // Id of h-refined surface
        const MInt sId = helements().hrefSurfaceIds(hId, dir, hSurfId);

        MInt hSId = -1;
        const auto it = srfcMap.find(sId);

        // For an existing h-refined surface check if a rb-surface exists
        if(sId != -1 && it != srfcMap.end()) {
          hSId = it->second;

          // (h-)rb-surface exists, just store the neighboring rb-element id
          m_rbSurfaces.nghbrElementIds(hSId, side) = rbId;
        } else if(sId != -1 && hSId == -1) {
          // rb-surface does not exist, create it if there is a h-refined
          // surface. This covers the cases of h-mpi surfaces and h-refinement
          // parallel to the boundary in the interior of the domain
          hSId = m_rbSurfaces.size();
          m_rbSurfaces.append();

          // Store surface id
          m_surfaceId.push_back(sId);

          // Copy variables of surface to rb-surface
          m_rbSurfaces.copy(surfaces(), sId, hSId);

          // Store neighboring rb-element ids in rb-surface
          m_rbSurfaces.nghbrElementIds(hSId, side) = rbId;
          m_rbSurfaces.nghbrElementIds(hSId, oppositeSide) = -1;
        }

        // Store h-rb-surface id in h-rb-element
        // -1 if there is no h-refinement in this direction
        m_hRbElements.hrefSurfaceIds(hRbId, dir, hSurfId) = hSId;
      }
    }
  }

  // Allocate space for RBC solution
  const MInt dataSize = Base::SysEqn::noVars() * ipow(maxNoNodes1D(), nDim);
  mDeallocate(m_rbSolution);
  if(noRbElements > 0) {
    mAlloc(m_rbSolution, noRbElements, dataSize, "m_rbSolution", F0, AT_);
  }
  m_internalDataSize = noRbElements * dataSize;

  RECORD_TIMER_STOP(m_timers[Timers::Init]);
}


/// \brief Initialize all timers and start the class timer.
///
/// \author Ansgar Niemoeller (ansgar) <a.niemoeller@aia.rwth-aachen.de>
/// \date 2016-01-24
template <MInt nDim>
void DgBcAcousticPerturbRBC<nDim>::initTimers() {
  TRACE();

  // Create timer group & timer for solver, and start the timer
  NEW_TIMER_GROUP_NOCREATE(m_timers[Timers::TimerGroup],
                           "DgBcAcousticPerturbRBC (solverId = " + std::to_string(solver().solverId())
                               + ", bcId = " + std::to_string(id()) + ")");
  NEW_TIMER_NOCREATE(m_timers[Timers::Class], "total object lifetime", m_timers[Timers::TimerGroup]);
  RECORD_TIMER_START(m_timers[Timers::Class]);

  // Create regular solver-wide timers
  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::Init], "init", m_timers[Timers::Class]);

  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::Apply], "apply", m_timers[Timers::Class]);
  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::InitializeVars], "initializeVars", m_timers[Timers::Apply]);
  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::CopySurfaceVars], "copySurfaceVars", m_timers[Timers::Apply]);
  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::Prolong], "prolong", m_timers[Timers::Apply]);
  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::TimeDeriv], "timeDeriv", m_timers[Timers::Apply]);
  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::RungeKuttaStep], "RungeKuttaStep", m_timers[Timers::Apply]);
  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::ApplyAtSurface], "applyAtSurface", m_timers[Timers::Apply]);

  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::SetRestartVars], "setRestartVars", m_timers[Timers::Class]);
  NEW_SUB_TIMER_NOCREATE(m_timers[Timers::GetRestartVars], "getRestartVars", m_timers[Timers::Class]);
}


/// \brief Apply RBC.
///
/// \author Ansgar Niemoeller (ansgar) <a.niemoeller@aia.rwth-aachen.de>
/// \date 2015-12-13
///
/// 0. Initialize:
///     - Calculate RBC node variables for rb-elements and rb-surfaces.
///     - Copy APE initial condition to rb-elements.
///
/// 1. Copy APE surface variables to rb-surfaces.
/// 2. Prolong RBC solution to boundary surfaces and copy to corresponding
///    rb-surfaces (both sides).
/// 3. Solve RBC system.
/// 4. Compute APE boundary surface flux with prolonged RBC solution as outer
///    state.
template <MInt nDim>
void DgBcAcousticPerturbRBC<nDim>::apply(const MFloat time) {
  using namespace maia::dg::interpolation;
  TRACE();

  // Nothing to be done if there are no RBC elements
  if(m_rbElements.size() == 0) {
    return;
  }

  RECORD_TIMER_START(m_timers[Timers::Apply]);

  // 0. Initialize node variables and set initial condition if not already done
  const MInt* const noNodes = &m_rbElements.noNodes1D(0);
  if(!m_isInitialized) {
    RECORD_TIMER_START(m_timers[Timers::InitializeVars]);

    // Set initial condition (same as for APE)
    if(!isRestart()) {
      for(MInt rbId = 0; rbId < m_rbElements.size(); rbId++) {
        const MInt eId = m_elementId[rbId];
        const MInt noNodes1D = noNodes[rbId];
        const MInt noNodesXD = ipow(noNodes1D, nDim);
        const MInt dataSize = noNodesXD * Base::SysEqn::noVars();

        // Copy element variables to rb-element
        std::copy_n(&elements().variables(eId), dataSize, &m_rbElements.variables(rbId));
      }
    }

    // Copy APE nodeVars to outer state of boundary surfaces
    for(MInt srfcId = begin(); srfcId < end(); srfcId++) {
      const MInt noNodes1D = surfaces().noNodes1D(srfcId);
      const MInt noNodes1D3 = (nDim == 3) ? noNodes1D : 1;
      const MInt noNodesXD = noNodes1D * noNodes1D3;
      const MInt internalSide = surfaces().internalSideId(srfcId);
      const MInt outerSide = (internalSide + 1) % 2;

      // Copy node vars to outer surface state
      std::copy_n(&surfaces().nodeVars(srfcId, internalSide),
                  noNodesXD * Base::SysEqn::noNodeVars(),
                  &surfaces().nodeVars(srfcId, outerSide));
    }

    // Also copy APE node vars to outer state of additional boundary surfaces
    const MInt noAddBcIds = m_addBndrySrfcIds.size();
    for(MInt sId = 0; sId < noAddBcIds; sId++) {
      const MInt rbSrfcId = m_addBndrySrfcIds[sId].first;
      const MInt srfcId = m_surfaceId[rbSrfcId];

      const MInt noNodes1D = surfaces().noNodes1D(srfcId);
      const MInt noNodes1D3 = (nDim == 3) ? noNodes1D : 1;
      const MInt noNodesXD = noNodes1D * noNodes1D3;
      const MInt internalSide = surfaces().internalSideId(srfcId);
      const MInt outerSide = (internalSide + 1) % 2;

      // Copy node vars to outer surface state
      std::copy_n(&surfaces().nodeVars(srfcId, internalSide),
                  noNodesXD * Base::SysEqn::noNodeVars(),
                  &surfaces().nodeVars(srfcId, outerSide));
    }

    // Calculate RBC node variables
    calcElementNodeVars();
    calcSurfaceNodeVars();

    // Set initial Runge Kutta stage
    m_rkStage = 0;

    m_isInitialized = true;

    RECORD_TIMER_STOP(m_timers[Timers::InitializeVars]);
  }

  // 1. Copy surface variables to rb-surfaces for surface flux computation
  // For boundary surfaces: variables are later overwritten with prolonged
  // rb-solutions
  RECORD_TIMER_START(m_timers[Timers::CopySurfaceVars]);

  for(MInt rbSId = 0; rbSId < m_rbSurfaces.size(); rbSId++) {
    const MInt sId = m_surfaceId[rbSId];

    const MInt noNodesXD = ipow(maxNoNodes1D(), nDim - 1);
    const MInt dataSize = noNodesXD * Base::SysEqn::noVars();

    std::copy_n(&surfaces().variables(sId, 0), dataSize, &m_rbSurfaces.variables(rbSId, 0));
    std::copy_n(&surfaces().variables(sId, 1), dataSize, &m_rbSurfaces.variables(rbSId, 1));
  }

  RECORD_TIMER_STOP(m_timers[Timers::CopySurfaceVars]);

  // 2.1 Prolong RBC solution to boundary surfaces
  RECORD_TIMER_START(m_timers[Timers::Prolong]);

  for(MInt srfcId = begin(); srfcId < end(); srfcId++) {
    // Get corresponding rb-element id for boundary surface
    const MInt rbId = m_rbElementId[srfcId - begin()];

    const MInt polyDeg = m_rbElements.polyDeg(rbId);
    const MInt noNodes1D = m_rbElements.noNodes1D(rbId);
    const MInt orientation = surfaces().orientation(srfcId);
    const MInt internalSide = surfaces().internalSideId(srfcId);
    const MInt outerSide = (internalSide + 1) % 2;
    const MInt dir = 2 * orientation + outerSide;

    // Prolong to outer state of boundary surface
    MFloat* src = &m_rbElements.variables(rbId);
    MFloat* dest = &surfaces().variables(srfcId, outerSide);

    if(integrationMethod() == DG_INTEGRATE_GAUSS) {
      prolongToFaceGauss<nDim, Base::SysEqn::noVars()>(src, dir, noNodes1D,
                                                       &interpolation(polyDeg, noNodes1D).m_LFace[0][0],
                                                       &interpolation(polyDeg, noNodes1D).m_LFace[1][0], dest);
    } else if(integrationMethod() == DG_INTEGRATE_GAUSS_LOBATTO) {
      prolongToFaceGaussLobatto<nDim, Base::SysEqn::noVars()>(src, dir, noNodes1D, dest);
    }

    const MInt rbSId = m_rbSurfaceId[srfcId - begin()];
    const MInt noNodes1D3 = (nDim == 3) ? noNodes1D : 1;
    const MInt noNodesXD = noNodes1D * noNodes1D3;
    const MInt dataSize = noNodesXD * Base::SysEqn::noVars();

    // Copy prolonged rbc solution to rb boundary surfaces
    std::copy_n(dest, dataSize, &m_rbSurfaces.variables(rbSId, 0));
    std::copy_n(dest, dataSize, &m_rbSurfaces.variables(rbSId, 1));
  }

  // 2.2 Prolong RBC solution to additional boundary surfaces
  const MInt noAddBcIds = m_addBndrySrfcIds.size();
  for(MInt srfcId = 0; srfcId < noAddBcIds; srfcId++) {
    const MInt rbSrfcId = m_addBndrySrfcIds[srfcId].first;
    const MInt rbId = m_addBndrySrfcIds[srfcId].second;

    const MInt polyDeg = m_rbElements.polyDeg(rbId);
    const MInt noNodes1D = m_rbElements.noNodes1D(rbId);
    const MInt orientation = m_rbSurfaces.orientation(rbSrfcId);
    const MInt internalSide = m_rbSurfaces.internalSideId(rbSrfcId);
    const MInt outerSide = (internalSide + 1) % 2;
    const MInt dir = 2 * orientation + outerSide;

    // Prolong to left state of boundary surface
    MFloat* src = &m_rbElements.variables(rbId);
    MFloat* dest = &m_rbSurfaces.variables(rbSrfcId, 0);

    if(integrationMethod() == DG_INTEGRATE_GAUSS) {
      prolongToFaceGauss<nDim, Base::SysEqn::noVars()>(src, dir, noNodes1D,
                                                       &interpolation(polyDeg, noNodes1D).m_LFace[0][0],
                                                       &interpolation(polyDeg, noNodes1D).m_LFace[1][0], dest);
    } else if(integrationMethod() == DG_INTEGRATE_GAUSS_LOBATTO) {
      prolongToFaceGaussLobatto<nDim, Base::SysEqn::noVars()>(src, dir, noNodes1D, dest);
    }

    const MInt noNodes1D3 = (nDim == 3) ? noNodes1D : 1;
    const MInt dataSize = noNodes1D * noNodes1D3 * Base::SysEqn::noVars();

    // Copy prolonged left state to right state
    std::copy_n(&m_rbSurfaces.variables(rbSrfcId, 0), dataSize, &m_rbSurfaces.variables(rbSrfcId, 1));
  }

  RECORD_TIMER_STOP(m_timers[Timers::Prolong]);

  // 3. Solve RBC system
  RECORD_TIMER_START(m_timers[Timers::TimeDeriv]);

  // Reset rbc-rhs
  resetBuffer(m_internalDataSize, &m_rbElements.rightHandSide(0));

  // Store rbc-solution and copy element variables to rb-elements for flux
  // computation
  for(MInt rbId = 0; rbId < m_rbElements.size(); rbId++) {
    const MInt eId = m_elementId[rbId];
    const MInt noNodes1D = noNodes[rbId];
    const MInt noNodesXD = ipow(noNodes1D, nDim);
    const MInt dataSize = noNodesXD * Base::SysEqn::noVars();

    // Store rb-element variables in m_rbSolution
    std::copy_n(&m_rbElements.variables(rbId), dataSize, m_rbSolution[rbId]);

    // Copy element variables to rb-element
    std::copy_n(&elements().variables(eId), dataSize, &m_rbElements.variables(rbId));
  }

  // Calculate volume integral (with APE variables)
  calcVolumeIntegral(m_rbElements.size(), m_rbElements, *this);

  // Calculate surface fluxes
  // Inner fluxes are computed with the APE variables on the surfaces
  // Boundary fluxes use the prolonged rbc solutions for both the left and right
  // state
  calcRegularSurfaceFlux(0, m_rbSurfaces.size(), m_rbSurfaces, *this);

  // Calculate surface integral
  calcSurfaceIntegral(0, m_rbElements.size(), m_rbElements, m_rbSurfaces, m_hRbElements, m_hRbElements.size());

  // Apply jacobian
  applyJacobian(m_rbElements.size(), m_rbElements);

  // Calculate source terms
  calcSourceTerms(time, m_rbElements.size(), m_rbElements, *this);

  // Copy rbc-solutions back to rb-elements for time integration
  for(MInt rbId = 0; rbId < m_rbElements.size(); rbId++) {
    const MInt noNodes1D = noNodes[rbId];
    const MInt noNodesXD = ipow(noNodes1D, nDim);
    const MInt dataSize = noNodesXD * Base::SysEqn::noVars();

    std::copy_n(m_rbSolution[rbId], dataSize, &m_rbElements.variables(rbId));
  }

  RECORD_TIMER_STOP(m_timers[Timers::TimeDeriv]);

  // Perform rk substep
  RECORD_TIMER_START(m_timers[Timers::RungeKuttaStep]);
  subTimeStepRk(dt(), m_rkStage, m_internalDataSize, &m_rbElements.rightHandSide(0), &m_rbElements.variables(0),
                &m_rbElements.timeIntStorage(0));

  // Define the number of Rk-steps according to the Scheme being used
  if(timeIntegrationScheme() == DG_TIMEINTEGRATION_CARPENTER_4_5) {
    m_rkStage = (m_rkStage + 1) % 5;
  } else if(timeIntegrationScheme() == DG_TIMEINTEGRATION_TOULORGEC_4_8) {
    m_rkStage = (m_rkStage + 1) % 8;
  } else if(timeIntegrationScheme() == DG_TIMEINTEGRATION_NIEGEMANN_4_14) {
    m_rkStage = (m_rkStage + 1) % 14;
  } else if(timeIntegrationScheme() == DG_TIMEINTEGRATION_NIEGEMANN_4_13) {
    m_rkStage = (m_rkStage + 1) % 13;
  } else if(timeIntegrationScheme() == DG_TIMEINTEGRATION_TOULORGEC_3_7) {
    m_rkStage = (m_rkStage + 1) % 7;
  } else if(timeIntegrationScheme() == DG_TIMEINTEGRATION_TOULORGEF_4_8) {
    m_rkStage = (m_rkStage + 1) % 8;
  }

  RECORD_TIMER_STOP(m_timers[Timers::RungeKuttaStep]);

  // 4. Apply boundary condition
  RECORD_TIMER_START(m_timers[Timers::ApplyAtSurface]);
  maia::dg::bc::loop(&DgBcAcousticPerturbRBC::applyAtSurface, this, begin(), end(), time);
  RECORD_TIMER_STOP(m_timers[Timers::ApplyAtSurface]);

  RECORD_TIMER_STOP(m_timers[Timers::Apply]);
}


/// \brief Calculate boundary surface flux.
///
/// \author Ansgar Niemoeller (ansgar) <a.niemoeller@aia.rwth-aachen.de>
/// \date 2015-12-13
template <MInt nDim>
void DgBcAcousticPerturbRBC<nDim>::applyAtSurface(const MInt surfaceId, const MFloat NotUsed(time)) {
  // TRACE();

  MFloat* nodeVarsL = &surfaces().nodeVars(surfaceId, 0);
  MFloat* nodeVarsR = &surfaces().nodeVars(surfaceId, 1);
  MFloat* stateL = &surfaces().variables(surfaceId, 0);
  MFloat* stateR = &surfaces().variables(surfaceId, 1);
  const MInt noNodes1D = surfaces().noNodes1D(surfaceId);
  const MInt dirId = surfaces().orientation(surfaceId);

  // The boundary state contains the prolonged rbc-solution and the same node
  // variables as the inner state
  MFloat* f = flux(surfaceId);
  sysEqn().calcRiemann(nodeVarsL, nodeVarsR, stateL, stateR, noNodes1D, dirId, f);
}


/// \brief Calculates the physical fluxes in all directions for all integration
///        points in an element.
///
/// \author Ansgar Niemoeller (ansgar) <a.niemoeller@aia.rwth-aachen.de>
/// \date 2015-12-13
///
/// Note: for argument description, see DgSysEqnAcousticPerturb::calcFlux.
template <MInt nDim>
void DgBcAcousticPerturbRBC<nDim>::calcFlux(const MFloat* nodeVars,
                                            const MFloat* q,
                                            const MInt noNodes1D,
                                            MFloat* flux_) const {
  // TRACE();

  const MInt noVars = Base::SysEqn::noVars();
  const MInt noNodes1D3 = (nDim == 3) ? noNodes1D : 1;
  const MFloatTensor U(const_cast<MFloat*>(q), noNodes1D, noNodes1D, noNodes1D3, noVars);
  const MFloatTensor coeff(const_cast<MFloat*>(nodeVars), noNodes1D, noNodes1D, noNodes1D3, s_noNodeVars);

  MFloatTensor f(flux_, noNodes1D, noNodes1D, noNodes1D3, nDim, noVars);

  // The following flux equations are calculated here (if 2D, consider all
  // z-components to be zero and sin(theta)=1):
  // vg: speed of wave propagation
  // theta, phi: angles in cylindrical coordinates

  // Flux in x-direction:
  // f1 = vg*sin(theta)*cos(phi)
  // f_x(CV[u]) = f1*u
  // f_x(CV[v]) = f1*v
  // f_x(CV[w]) = f1*w
  // f_x(CV[p]) = f1*p

  // Flux in y-direction:
  // f2 = vg*sin(theta)*sin(phi)
  // f_y(CV[u]) = f2*u
  // f_y(CV[v]) = f2*v
  // f_y(CV[w]) = f2*w
  // f_y(CV[p]) = f2*p

  // Flux in z-direction:
  // f3 = vg*cos(theta)
  // f_z(CV[u]) = f3*u
  // f_z(CV[v]) = f3*v
  // f_z(CV[w]) = f3*w
  // f_z(CV[p]) = f3*p

  for(MInt i = 0; i < noNodes1D; i++) {
    for(MInt j = 0; j < noNodes1D; j++) {
      for(MInt k = 0; k < noNodes1D3; k++) {
        for(MInt d = 0; d < nDim; d++) {
          for(MInt var = 0; var < noVars; var++) {
            f(i, j, k, d, var) = coeff(i, j, k, d) * U(i, j, k, var);
          }
        }
      }
    }
  }
}


/// \brief Calculates the physical fluxes in direction 'dirId' on a surface.
///
/// \author Ansgar Niemoeller (ansgar) <a.niemoeller@aia.rwth-aachen.de>
/// \date 2015-12-13
///
/// The flux is obtained by the projection of the speed of wave propagation on
/// the given direction multiplied with the variables.
/// Note: for argument description, see DgSysEqnAcousticPerturb::calcFlux1D.
template <MInt nDim>
void DgBcAcousticPerturbRBC<nDim>::calcFlux1D(const MFloat* nodeVars, const MFloat* q, const MInt noNodes1D,
                                              const MInt dirId, MFloat* flux_) const {
  // TRACE();

  const MInt noVars = Base::SysEqn::noVars();
  const MInt noNodes1D3 = (nDim == 3) ? noNodes1D : 1;
  const MFloatTensor U(const_cast<MFloat*>(q), noNodes1D, noNodes1D3, noVars);
  const MFloatTensor coeff(const_cast<MFloat*>(nodeVars), noNodes1D, noNodes1D3, s_noNodeVars);
  MFloatTensor f(flux_, noNodes1D, noNodes1D3, noVars);

  for(MInt i = 0; i < noNodes1D; i++) {
    for(MInt j = 0; j < noNodes1D3; j++) {
      for(MInt var = 0; var < noVars; var++) {
        f(i, j, var) = coeff(i, j, dirId) * U(i, j, var);
      }
    }
  }
}


/// \brief Calculates the source terms of the RBC system.
///
/// \author Ansgar Niemoeller (ansgar) <a.niemoeller@aia.rwth-aachen.de>
/// \date 2015-12-13
///
/// The RBC source terms are given by:
///  2D: S = vg/(2*r)*U
///  3D: S = vg/r*U
/// Note: for argument description, see DgSysEqnAcousticPerturb::calcSource.
template <MInt nDim>
void DgBcAcousticPerturbRBC<nDim>::calcSource(const MFloat* nodeVars,
                                              const MFloat* u,
                                              const MInt noNodes1D,
                                              const MFloat NotUsed(t),
                                              const MFloat* NotUsed(x),
                                              MFloat* src) const {
  // TRACE();

  const MInt noVars = Base::SysEqn::noVars();

  const MInt noNodes1D3 = (nDim == 3) ? noNodes1D : 1;
  MFloatTensor s(src, noNodes1D, noNodes1D, noNodes1D3, noVars);

  const MFloatTensor U(const_cast<MFloat*>(u), noNodes1D, noNodes1D, noNodes1D3, noVars);
  const MFloatTensor coeff(const_cast<MFloat*>(nodeVars), noNodes1D, noNodes1D, noNodes1D3, nDim + 2);

  for(MInt i = 0; i < noNodes1D; i++) {
    for(MInt j = 0; j < noNodes1D; j++) {
      for(MInt k = 0; k < noNodes1D3; k++) {
        for(MInt var = 0; var < noVars; var++) {
          s(i, j, k, var) = coeff(i, j, k, nDim + 1) * U(i, j, k, var);
        }
      }
    }
  }
}


/// \brief Calculates the numerical flux at a surface given two states.
///
/// \author Ansgar Niemoeller (ansgar) <a.niemoeller@aia.rwth-aachen.de>
/// \date 2015-12-13
///
/// The numerical flux is calculated using the local Lax-Friedrichs flux scheme.
/// Note: for argument description, see DgSysEqnAcousticPerturb::calcRiemann.
template <MInt nDim>
void DgBcAcousticPerturbRBC<nDim>::calcRiemann(const MFloat* nodeVarsL,
                                               const MFloat* nodeVarsR,
                                               const MFloat* stateL,
                                               const MFloat* stateR,
                                               const MInt noNodes1D,
                                               const MInt dirId,
                                               MFloat* flux_) const {
  // TRACE();

  const MInt noVars = Base::SysEqn::noVars();
  const MInt noNodes1D3 = (nDim == 3) ? noNodes1D : 1;
  const MFloatTensor uL(const_cast<MFloat*>(stateL), noNodes1D, noNodes1D3, noVars);
  const MFloatTensor uR(const_cast<MFloat*>(stateR), noNodes1D, noNodes1D3, noVars);

#ifndef _OPENMP
  MFloatScratchSpace maxLambda(noNodes1D, noNodes1D3, AT_, "maxLambda");
#else
  MFloatTensor maxLambda(noNodes1D, noNodes1D3);
#endif

  // Left and right node variables
  const MFloatTensor cL(const_cast<MFloat*>(nodeVarsL), noNodes1D, noNodes1D3, s_noNodeVars);
  const MFloatTensor cR(const_cast<MFloat*>(nodeVarsR), noNodes1D, noNodes1D3, s_noNodeVars);

  // Calculate maximum eigenvalue from speed of sound propagation
  for(MInt i = 0; i < noNodes1D; i++) {
    for(MInt j = 0; j < noNodes1D3; j++) {
      maxLambda(i, j) = std::max(fabs(cL(i, j, nDim)), fabs(cR(i, j, nDim)));
    }
  }

#ifndef _OPENMP
  MFloatScratchSpace fluxL(noNodes1D, noNodes1D3, noVars, AT_, "fluxL");
  MFloatScratchSpace fluxR(noNodes1D, noNodes1D3, noVars, AT_, "fluxR");
#else
  MFloatTensor fluxL(noNodes1D, noNodes1D3, noVars);
  MFloatTensor fluxR(noNodes1D, noNodes1D3, noVars);
#endif

  // Calculate flux from left and right state
  calcFlux1D(nodeVarsL, stateL, noNodes1D, dirId, &fluxL[0]);
  calcFlux1D(nodeVarsR, stateR, noNodes1D, dirId, &fluxR[0]);

  // Solve Riemann problem
  MFloatTensor riemann(flux_, noNodes1D, noNodes1D3, noVars);
  for(MInt i = 0; i < noNodes1D; i++) {
    for(MInt j = 0; j < noNodes1D3; j++) {
      for(MInt n = 0; n < noVars; n++) {
        riemann(i, j, n) = 0.5 * ((fluxL(i, j, n) + fluxR(i, j, n)) - maxLambda(i, j) * (uR(i, j, n) - uL(i, j, n)));
      }
    }
  }
}


/// \brief Calculate node variables of rb-elements.
///
/// \author Ansgar Niemoeller (ansgar) <a.niemoeller@aia.rwth-aachen.de>
/// \date 2015-12-13
///
/// Calculate all node variables on the rb-elements using the node variables on
/// the 'normal' elements.
template <MInt nDim>
void DgBcAcousticPerturbRBC<nDim>::calcElementNodeVars() {
  TRACE();

  // Loop over rb-elements
  for(MInt rbId = 0; rbId < m_rbElements.size(); rbId++) {
    const MInt eId = m_elementId[rbId]; // Corresponding element id
    const MInt noNodesXD = m_rbElements.noNodesXD(rbId);
    const MFloat* coordinates = &m_rbElements.nodeCoordinates(rbId);

    for(MInt pId = 0; pId < noNodesXD; pId++) {
      calcNodeVars(&coordinates[pId * nDim],
                   &elements().nodeVars(eId) + pId * Base::SysEqn::noNodeVars(),
                   &m_rbElements.nodeVars(rbId) + pId * s_noNodeVars);
    }
  }
}


/// \brief Calculate node variables of rb-surfaces
///
/// \author Ansgar Niemoeller (ansgar) <a.niemoeller@aia.rwth-aachen.de>
/// \date 2015-12-13
///
/// Calculate all node variables on the rb-surfaces using the node variables on
/// the 'normal' surfaces.
template <MInt nDim>
void DgBcAcousticPerturbRBC<nDim>::calcSurfaceNodeVars() {
  TRACE();

  // Loop over rb-surfaces
  for(MInt rbSId = 0; rbSId < m_rbSurfaces.size(); rbSId++) {
    const MInt sId = m_surfaceId[rbSId];
    const MFloat* coordinates = &m_rbSurfaces.nodeCoords(rbSId);
    const MInt noNodesXD = m_rbSurfaces.noNodesXD(rbSId);

    for(MInt pId = 0; pId < noNodesXD; pId++) {
      calcNodeVars(&coordinates[pId * nDim],
                   &((&surfaces().nodeVars(sId, 0))[pId * Base::SysEqn::noNodeVars()]),
                   &((&m_rbSurfaces.nodeVars(rbSId, 0))[pId * s_noNodeVars]));

      calcNodeVars(&coordinates[pId * nDim],
                   &((&surfaces().nodeVars(sId, 1))[pId * Base::SysEqn::noNodeVars()]),
                   &((&m_rbSurfaces.nodeVars(rbSId, 1))[pId * s_noNodeVars]));
    }
  }
}


/// \brief Calculate RBC node variables at a given position.
///
/// \author Ansgar Niemoeller (ansgar) <a.niemoeller@aia.rwth-aachen.de>
/// \date 2015-12-13
///
/// The 'nDim + 2' RBC node variables are computed at a given position with the
/// corresponding mean velocities.
/// The first nDim node variables contain the projection of the speed of wave
/// propagation from the acoustic source region onto the Cartesian coordinate
/// directions. These are the needed coefficients for the flux computations in
/// the corresponding directions.
/// At position 'nDim + 1' the speed of wave propagation is stored, which is
/// needed for the Riemann solver.
/// The last entry contains the source term coefficient which is given by the
/// ratio of the speed of wave propagation and the distance to the reference
/// position (2D: with factor of 1/2).
///
/// \param[in] pos position where to calculate the node variables.
/// \param[in] nodeVarsAPE node variables of the APE-system at given position.
/// \param[out] nodeVars computed RBC node variables.
template <MInt nDim>
void DgBcAcousticPerturbRBC<nDim>::calcNodeVars(const MFloat* pos, const MFloat* nodeVarsAPE, MFloat* nodeVars) {
  // Mean speed of sound
  const MFloat c = nodeVarsAPE[Base::SysEqn::CV::C0];

  IF_CONSTEXPR(nDim == 2) { // 2D
    const MFloat x = pos[0];
    const MFloat y = pos[1];

    const MFloat x0 = m_rbcReferencePos[0];
    const MFloat y0 = m_rbcReferencePos[1];

    // Inverse distance of given position and reference position
    const MFloat rInv = 1.0 / sqrt(POW2(x - x0) + POW2(y - y0));

    // Orientation of given position with respect to the reference position
    const MFloat cosPhi = (x - x0) * rInv;
    const MFloat sinPhi = (y - y0) * rInv;

    // Mean velocities
    const MFloat u0 = nodeVarsAPE[Base::SysEqn::CV::U0];
    const MFloat v0 = nodeVarsAPE[Base::SysEqn::CV::V0];

    // Speed of wave propagation
    //
    // v_g = u0 * cos(phi) + v0 * sin(phi)
    //       + sqrt(c^2 - (-u0 * sin(phi) + v0 * cos(phi))^2)
    //
    // c: mean sound speed
    // u0, v0: mean velocities
    // phi: angle between given position and reference position
    //
    // Reference: see 3D case
    const MFloat vg = u0 * cosPhi + v0 * sinPhi + sqrt(POW2(c) - POW2(-u0 * sinPhi + v0 * cosPhi));

    // Projection of speed of wave propagation in both Cartesian coordinate
    // directions (coefficients for flux computation)
    nodeVars[0] = vg * cosPhi;
    nodeVars[1] = vg * sinPhi;

    // Speed of wave propagation (needed in calcRiemann)
    nodeVars[2] = vg;

    // Source term coefficient
    nodeVars[3] = 0.5 * vg * rInv;
  }
  else { // 3D
    const MFloat x = pos[0];
    const MFloat y = pos[1];
    const MFloat z = pos[2];

    const MFloat x0 = m_rbcReferencePos[0];
    const MFloat y0 = m_rbcReferencePos[1];
    const MFloat z0 = m_rbcReferencePos[2];

    // Inverse distance of given position and reference position
    const MFloat rInv = 1.0 / sqrt(POW2(x - x0) + POW2(y - y0) + POW2(z - z0));

    // Inverse distance in x-y-plane of given position and reference position
    const MFloat r2Inv = 1.0 / sqrt(POW2(x - x0) + POW2(y - y0));

    // Orientation of given position with respect to the reference position
    // (spherical coordinates (r, theta, phi))
    const MFloat cosPhi = (x - x0) * r2Inv;
    const MFloat sinPhi = (y - y0) * r2Inv;
    const MFloat cosTheta = (z - z0) * rInv;
    const MFloat sinTheta = sqrt(1 - POW2(cosTheta));

    // Mean velocities
    const MFloat u0 = nodeVarsAPE[Base::SysEqn::CV::UU0[0]];
    const MFloat v0 = nodeVarsAPE[Base::SysEqn::CV::UU0[1]];
    const MFloat w0 = nodeVarsAPE[Base::SysEqn::CV::UU0[2]];

    // Speed of wave propagation in radial direction
    const MFloat u_er = u0 * sinTheta * cosPhi + v0 * sinTheta * sinPhi + w0 * cosTheta;

    // Speed of wave propagation in theta direction
    const MFloat u_et = u0 * cosTheta * cosPhi + v0 * cosTheta * sinPhi + w0 * (-sinTheta);

    // Speed of wave propagation in phi direction
    const MFloat u_ep = u0 * (-sinPhi) + v0 * cosPhi;

    // Speed of wave propagation
    //
    // v_g = u_m*e_r + sqrt(c^2 - (u_m*e_t)^2 - (u_m*e_p)^2)
    //
    // u_m: mean flow velocity vector u_m = [u0, v0, w0]'
    // e_r: unit vector in radial direction
    //      e_r = [sin(phi)*cos(theta), sin(theta)*sin(phi), cos(theta)]'
    // e_t: unit vector in theta direction
    //      e_t = [cos(theta)*cos(phi), cos(theta)*sin(phi), -sin(theta)]'
    // e_p: unit vector in phi direction
    //      e_p = [-sin(phi), cos(phi), 0]'
    // c: mean sound speed
    //
    // Reference: C .Bogey, C. Bailly; Three-dimensional non-reflective boundary
    // conditions for acoustic simulations: far field formulation and validation
    // test cases, Eqns. (1) and (2)
    const MFloat vg = u_er + sqrt(POW2(c) - POW2(u_et) - POW2(u_ep));

    // Projection of speed of wave propagation in all Cartesian coordinate
    // directions (coefficients for flux computation)
    nodeVars[0] = vg * sinTheta * cosPhi;
    nodeVars[1] = vg * sinTheta * sinPhi;
    nodeVars[2] = vg * cosTheta;

    // Speed of wave propagation (needed in calcRiemann)
    nodeVars[3] = vg;

    // Source term coefficient
    nodeVars[4] = vg * rInv;
  }
}


/// \brief Return local number of nodes.
///
/// \author Ansgar Niemoeller (ansgar) <a.niemoeller@aia.rwth-aachen.de>
/// \date 2016-01-12
template <MInt nDim>
MInt DgBcAcousticPerturbRBC<nDim>::getLocalNoNodes() const {
  TRACE();

  MInt localNoNodes = 0;
  for(MInt rbId = 0; rbId < m_rbElements.size(); rbId++) {
    const MInt noNodesXD = m_rbElements.noNodesXD(rbId);

    localNoNodes += noNodesXD;
  }

  return localNoNodes;
}


/// \brief Return name of restart variable.
///
/// \author Ansgar Niemoeller (ansgar) <a.niemoeller@aia.rwth-aachen.de>
/// \date 2016-01-12
///
/// \param[in] id_ Requested variable id.
template <MInt nDim>
MString DgBcAcousticPerturbRBC<nDim>::restartVarName(const MInt id_) const {
  TRACE();

  std::stringstream varName;
  varName << "bc" << id() << "_" << Base::SysEqn::consVarNames(id_);

  return varName.str();
}


/// \brief Copy restart variable data to boundary conditon.
///
/// \author Ansgar Niemoeller (ansgar) <a.niemoeller@aia.rwth-aachen.de>
/// \date 2016-01-12
///
/// \param[in] id_ Variable id.
/// \param[in] data Pointer to restart variable data.
template <MInt nDim>
void DgBcAcousticPerturbRBC<nDim>::setRestartVariable(const MInt id_, const MFloat* const data) {
  TRACE();
  RECORD_TIMER_START(m_timers[Timers::SetRestartVars]);

  // Copy variables to rb elements (in ascending cell/element id order)
  MInt offset = 0;
  for(MInt rbId = 0; rbId < m_rbElements.size(); rbId++) {
    const MInt noNodesXD = m_rbElements.noNodesXD(rbId);
    MFloat* const vars = &m_rbElements.variables(rbId);
    for(MInt n = 0; n < noNodesXD; n++) {
      vars[n * Base::SysEqn::noVars() + id_] = data[offset + n];
    }

    offset += noNodesXD;
  }

  RECORD_TIMER_STOP(m_timers[Timers::SetRestartVars]);
}


/// \brief Copy restart variable data to pointer.
///
/// \author Ansgar Niemoeller (ansgar) <a.niemoeller@aia.rwth-aachen.de>
/// \date 2016-01-12
///
/// \param[in] id_ Variable id.
/// \param[in] data Pointer to store restart variable data.
template <MInt nDim>
void DgBcAcousticPerturbRBC<nDim>::getRestartVariable(const MInt id_, MFloat* const data) const {
  TRACE();
  RECORD_TIMER_START(m_timers[Timers::GetRestartVars]);

  // Copy variables to data buffer (in ascending cell/element id order)
  MInt offset = 0;
  for(MInt rbId = 0; rbId < m_rbElements.size(); rbId++) {
    const MInt noNodesXD = m_rbElements.noNodesXD(rbId);

    for(MInt n = 0; n < noNodesXD; n++) {
      data[offset + n] = m_rbElements.variables(rbId, n * Base::SysEqn::noVars() + id_);
    }

    offset += noNodesXD;
  }

  RECORD_TIMER_STOP(m_timers[Timers::GetRestartVars]);
}


/// \brief Solution data size for given cell and data id.
///
/// \author Ansgar Niemoeller (ansgar) <a.niemoeller@aia.rwth-aachen.de>
template <MInt nDim>
MInt DgBcAcousticPerturbRBC<nDim>::cellDataSizeDlb(const MInt dataId, const MInt cellId) const {
  TRACE();

  // Only the variables need to be considered
  if(dataId < 0 || dataId > Base::SysEqn::noVars()) {
    TERMM(1, "Invalid dataId.");
  }

  // Get element id and find corresponding RB-element (if there is one)
  const MInt elementId = getElementByCellId(cellId);
  auto rbElemIt = std::find(m_elementId.begin(), m_elementId.end(), elementId);

  MInt dataSize = 0;
  if(rbElemIt != m_elementId.end()) {
    // Compute rb-element id
    const MInt rbId = std::distance(m_elementId.begin(), rbElemIt);
    const MInt noNodesXD = m_rbElements.noNodesXD(rbId);

    dataSize = noNodesXD;
  }

  return dataSize;
}


/// \brief Check if there is a boundary condition element associated with the given element
///
/// \author Ansgar Niemoeller (ansgar) <a.niemoeller@aia.rwth-aachen.de>
template <MInt nDim>
MBool DgBcAcousticPerturbRBC<nDim>::hasBcElement(const MInt elementId) const {
  TRACE();

  return std::find(m_elementId.begin(), m_elementId.end(), elementId) != m_elementId.end();
}


#endif // DGBOUNDARYCONDITIONACOUSTICPERTURBRBC_H_
