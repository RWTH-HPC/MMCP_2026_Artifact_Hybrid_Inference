// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef DGSPONGE_H_
#define DGSPONGE_H_

#include <algorithm>
#include <array>
#include <limits>
#include <memory>
#include <numeric>
#include <set>
#include <utility>
#include <vector>

#include "COMM/mpioverride.h"
#include "GRID/cartesiangridproxy.h"
#include "INCLUDE/maiatypes.h"
#include "IO/context.h"
#include "MEMORY/collector.h"
#include "MEMORY/scratch.h"
#include "UTIL/functions.h"
#include "config.h"
#include "dgcartesianboundarycondition.h"
#include "dgcartesianboundaryconditionfactory.h"
#include "dgcartesianelementcollector.h"
#include "dgcartesianspongeelementcollector.h"
#include "dgcartesiansurfacecollector.h"
#include "property.h"
#include "typetraits.h"

/**
 * \brief Container for attributes characterizing sponge layer.
 *
 * \author Hans Yu <hans.yu@rwth-aachen.de>
 * \date April 16th, 2015
 */
template <MInt nDim>
struct DgSpongeLayer {
  MInt m_bcId;
  MInt m_orientation;
  MFloat m_offset;
  MFloat m_thickness;
  std::array<MFloat, nDim> m_minCoord;
  std::array<MFloat, nDim> m_maxCoord;
};


/**
 * \brief Container for sponge elements.
 *
 * \author Hans Yu <hans.yu@rwth-aachen.de>
 * \date April 8th, 2015
 */
template <MInt nDim, class SysEqn>
class DgSponge {
  using BC = typename DgBoundaryConditionFactory<nDim, SysEqn>::ReturnType;
  using Grid = typename maia::grid::Proxy<nDim>;
  using SpongeElementCollector = maia::dg::collector::SpongeElementCollector<nDim, SysEqn>;
  using SpongeLayerType = DgSpongeLayer<nDim>;
  using ElementCollector = maia::dg::collector::ElementCollector<nDim, SysEqn>;
  using SurfaceCollector = maia::dg::collector::SurfaceCollector<nDim, SysEqn>;

 public:
  // Methods
  /// Return the MPI communicator used by the corresponding solver
  MPI_Comm mpiComm() const { return m_mpiComm; }
  /// Return the domain id of the corresponding solver
  MInt domainId() const { return m_domainId; }
  /// Return the number of domains of the corresponding solver
  MInt noDomains() const { return m_noDomains; }
  DgSponge(const MInt solverId, const MPI_Comm comm);
  void init(const MInt maxPolyDeg, Grid* grid_, ElementCollector* elements_, SurfaceCollector* surfaces_,
            std::vector<BC>* boundaryConditions, SysEqn* sysEqn, const MPI_Comm comm);
  void calcSourceTerms();

  // Accessor methods
  MInt noSpongeElements() const { return m_spongeElements.size(); }
  MInt elementId(const MInt seId) const { return m_spongeElements.elementId(seId); }
  MInt spongeElementId(const MInt eId) const;
  MFloat spongeEta(const MInt seId, const MInt pos) { return m_spongeElements.spongeEta(seId, pos); }

 private:
  // Sponge Layer
  void checkSpongeBoundaryConditions() const;
  void initSpongeElements();
  void initSpongeLayer(SpongeLayerType& spongeLayer, const BC& boundaryCondition);
  void exchangeSpongeLayers(std::vector<SpongeLayerType>& spongeLayers);
  void calcSpongeEtaForAllNodes(const MInt seId, const std::vector<SpongeLayerType>& spongeLayers);
  MFloat distance(const SpongeLayerType& spongeLayer, const MFloat* otherPoint) const;
  MBool pointInsideSpongeLayer(const SpongeLayerType& spongeLayer, const MFloat* otherPoint, MFloat& diff) const;

  // Solver variables
  const MInt m_solverId = -1;
  MPI_Comm m_mpiComm = MPI_COMM_NULL;
  MInt m_domainId = -1;
  MInt m_noDomains = -1;

  // Grid variables
  Grid* m_grid = nullptr;

  // Polynomial degree variables
  MInt m_maxPolyDeg = -1;

  // Element variables
  ElementCollector* m_elements = nullptr;

  // Boundary condition variables
  std::vector<BC>* m_boundaryConditions = nullptr;

  // Surface variables
  SurfaceCollector* m_surfaces = nullptr;

  // Source calculation variables
  SysEqn* m_sysEqn = nullptr;

  // Sponge elements
  SpongeElementCollector m_spongeElements;
};


/// \brief Constructor accepts solver id and MPI communicator.
///
/// \author Michael Schlottke (mic) <mic@aia.rwth-aachen.de>
/// \date 2015-08-31
///
/// \param[in] solverId The id of the corresponding solver.
/// \param[in] comm The MPI communicator to use.
template <MInt nDim, class SysEqn>
DgSponge<nDim, SysEqn>::DgSponge(const MInt solverId, const MPI_Comm comm) : m_solverId(solverId), m_mpiComm(comm) {
  TRACE();
}


/**
 * \brief Sets attributes so that sponge can access necessary components in
 *solver. Checks sponge boundary condition ids. Creates sponge elements.
 *
 * \author Hans Yu <hans.yu@rwth-aachen.de>
 * \date April 8th, 2015
 *
 * \param[in] solverId
 * \param[in] maxPolyDeg
 * \param[in] grid
 * \param[in] elements
 * \param[in] surfaces
 * \param[in] boundaryConditions
 * \param[in] sysEqn
 */
template <MInt nDim, class SysEqn>
void DgSponge<nDim, SysEqn>::init(const MInt maxPolyDeg, Grid* grid_, ElementCollector* elements_,
                                  SurfaceCollector* surfaces_, std::vector<BC>* boundaryConditions, SysEqn* sysEqn,
                                  const MPI_Comm comm) {
  TRACE();

  m_mpiComm = comm;
  // Determine domain id and number of domains
  MPI_Comm_rank(mpiComm(), &m_domainId);
  MPI_Comm_size(mpiComm(), &m_noDomains);

  m_maxPolyDeg = maxPolyDeg;
  m_grid = grid_;
  m_elements = elements_;
  m_surfaces = surfaces_;
  m_boundaryConditions = boundaryConditions;
  m_sysEqn = sysEqn;

  // Make sure that all BCs marked as spongeBCs are actually present in
  // simulation setup (sanity check)
  checkSpongeBoundaryConditions();

  // Intialize a sponge element for each DG element affected by a sponge layer
  initSpongeElements();
}


/**
 * \brief Check whether boundary conditions specified to have a sponge exist.
 *
 * \author Hans Yu <hans.yu@rwth-aachen.de>
 * \date April 20th, 2015
 */
template <MInt nDim, class SysEqn>
void DgSponge<nDim, SysEqn>::checkSpongeBoundaryConditions() const {
  TRACE();

  // Get number of boundary conditions over all domains
  const MInt noBcs = m_boundaryConditions->size();
  const MInt noDomains1 = (domainId() != 0) ? 1 : noDomains();
  MIntScratchSpace noBcsPerDomain(noDomains1, AT_, "noBcsPerDomain");
  MPI_Gather(const_cast<MInt*>(&noBcs), 1, maia::type_traits<MInt>::mpiType(), &noBcsPerDomain[0], 1,
             maia::type_traits<MInt>::mpiType(), 0, mpiComm(), AT_, "const_cast<MInt*>(&noBcs)", "noBcsPerDomain[0]");

  /////////////////////////////////////////////////////////////////////////////
  // Gather information on root
  /////////////////////////////////////////////////////////////////////////////

  // sendbuf
  const MInt noBcs1 = (noBcs == 0) ? 1 : noBcs;
  MIntScratchSpace bcIds(noBcs1, AT_, "bcIds");
  for(MInt i = 0; i < noBcs; i++) {
    bcIds[i] = (*m_boundaryConditions)[i]->id();
  }

  // recvbuf
  const MInt noBcsTotal1 = (domainId() != 0) ? 1 : std::accumulate(noBcsPerDomain.begin(), noBcsPerDomain.end(), 0);
  MIntScratchSpace bcIdsPerDomain(noBcsTotal1, AT_, "bcIdsPerDomain");

  // displs
  MIntScratchSpace displs(noDomains1, AT_, "displs");
  displs[0] = 0;
  for(MInt i = 1; i < noDomains1; i++) {
    displs[i] = displs[i - 1] + noBcsPerDomain[i - 1];
  }

  // Get boundary condition ids for all domains
  MPI_Gatherv(&bcIds[0], noBcs, maia::type_traits<MInt>::mpiType(), &bcIdsPerDomain[0], &noBcsPerDomain[0], &displs[0],
              maia::type_traits<MInt>::mpiType(), 0, mpiComm(), AT_, "bcIds[0]", "bcIdsPerDomain[0]");

  // Copy to set to eliminate repeated entries
  std::set<MInt> bcIdsUnique(bcIdsPerDomain.begin(), bcIdsPerDomain.end());

  /////////////////////////////////////////////////////////////////////////////
  // Read properties and set m_spongeBoundaryConditionIds
  /////////////////////////////////////////////////////////////////////////////

  // get property
  // MProperty* spongeBcProperty = nullptr;
  if(Context::propertyExists("spongeBCs", m_solverId)) {
    /*! \property
      \page propertyPageDG DG
      \section spongeBCs
      <code>MInt spongeBcProperty</code>\n
      default = <code>n/a</code>\n \n
      Specifies at which interval the state at the specified points should be
      sampled.\n
      Possible values are:
      <ul>
        <li>any valid boundary condition id (i.e., a bcId that is actually in use)</li>
      </ul>
      Keywords: <i>DISCONTINUOUS_GALERKIN, SPONGE, BOUNDARY CONDITION</i>
    */
  } else {
    TERMM(1, "Sponge is activated, but spongeBCs property is not defined.");
  }

  // read sponge boundary condition ids
  const MInt noSpongeBcs = Context::propertyLength("spongeBCs", m_solverId);
  for(MInt i = 0; i < noSpongeBcs; i++) {
    const MInt spongeBcId = Context::getSolverProperty<MInt>("spongeBCs", m_solverId, AT_, i);
    // check whether boundary condition is valid
    if(domainId() == 0 && bcIdsUnique.count(spongeBcId) == 0) {
      TERMM(1, "BC " + std::to_string(spongeBcId) + " in spongeBCs is not a valid boundary condition id.");
    }
  }
}


/**
 * \brief Initialisation of SpongeLayer. Calculation of spongeEta.
 *        Sponge is an additional Sourceoperator:
 *        DelL: du/dt + L(u) = DelL(u)
 *
 *        DelL = spongeSigma * spongeEta * DeltaStates
 *        spongeSigma = heuristical value (approx 0.5)
 *        spongeEta = (x_sp / L_sp)^2
 *        L_sp = spongeLayerThickness
 *        x_sp = distance of a point to the inner spongeLayer
 *        DeltaStages = state_infinity - state_current
 *
 *        Literature: Hartmann Diss. eq. 3.58
 *
 * \author Vitali Pauz <v.pauz@aia.rwth-aachen.de>, Sven Berger, Hans Yu
 * <hans.yu@rwth-aachen.de>
 * \date April 20th, 2015
 */
template <MInt nDim, class SysEqn>
void DgSponge<nDim, SysEqn>::initSpongeElements() {
  TRACE();

  // Setting default sponge thickness to nonsense value, so it can be easily
  // discerned later if default sponge thickness applies.
  MFloat defaultThickness = -std::numeric_limits<MFloat>::infinity();
  /*! \property
    \page propertyPageDG DG
    \section defaultSpongeThickness
    <code>MFloat defaultThickness</code>\n
    default = <code>-infinity</code>\n \n
    Default thickness of the sponge layer. Can be overridden for each sponge boundary condition.\n
    Possible values are:
    <ul>
      <li>any non-negative real value</li>
    </ul>
    Keywords: <i>DISCONTINUOUS_GALERKIN, SPONGE, THICKNESS</i>
  */
  defaultThickness = Context::getSolverProperty<MFloat>("defaultSpongeThickness", m_solverId, AT_, &defaultThickness);

  // get property
  const MInt noSpongeBcs = Context::propertyLength("spongeBCs", m_solverId);

  // 1. for sponge bc's calculate a boundary plane
  std::vector<SpongeLayerType> spongeLayers(noSpongeBcs);
  for(MInt i = 0; i < noSpongeBcs; i++) {
    // Store for convenience
    SpongeLayerType& layer = spongeLayers[i];

    // Determine sponge BC id
    const MInt spongeBcId = Context::getSolverProperty<MInt>("spongeBCs", m_solverId, AT_, i);

    layer.m_bcId = spongeBcId;

    // Determine thickness
    const MString propName = "spongeThickness_" + std::to_string(spongeBcId);
    /*! \page propertyPageDG DG
      \section spongeThickness_N
      <code>MFloat layer.m_thickness</code>\n
      default = <code>defaultSpongeThickness</code>\n \n
      Thickness of the sponge layer for N-th sponge boundary condition.\n
      Possible values are:
      <ul>
        <li>any non-negative real value</li>
      </ul>
      Keywords: <i>DISCONTINUOUS_GALERKIN, SPONGE, THICKNESS</i>
    */
    layer.m_thickness = Context::getSolverProperty<MFloat>(propName, m_solverId, AT_, &defaultThickness);

    if(layer.m_thickness < 0.0) {
      TERMM(1,
            "No valid definition for spongeThickness_" + std::to_string(spongeBcId)
                + ". Specify positive value or check your default sponge size.");
    }

    // Initialize sponge layer
    const auto boundaryConditionIt = std::find_if(m_boundaryConditions->begin(),
                                                  m_boundaryConditions->end(),
                                                  [spongeBcId](const BC& bc) { return bc->id() == spongeBcId; });
    if(boundaryConditionIt == m_boundaryConditions->end()) {
      // This should never happen, since all boundary conditions are initialized
      // on all domains even if no actual boundary surfaces exist locally
      TERMM(1, "Boundary condition id " + std::to_string(spongeBcId) + " not found");
    }
    if((*boundaryConditionIt)->count() > 0) {
      // If sponge BC id has surfaces on current domain, initialize it properly
      const BC& boundaryCondition = *boundaryConditionIt;
      initSpongeLayer(layer, boundaryCondition);
    } else {
      // Otherwise use *sensible* default values for minimum and maximum extent
      std::fill_n(layer.m_minCoord.begin(), nDim, std::numeric_limits<MFloat>::max());
      std::fill_n(layer.m_maxCoord.begin(), nDim, -std::numeric_limits<MFloat>::max());
    }
  }

  // 1.A if necessary exchange boundary planes
  exchangeSpongeLayers(spongeLayers);

  // 2. allocate memory for the sponge layer
  MBoolScratchSpace spongeElementFlags(m_elements->size(), AT_, "spongeElementFlags");
  std::fill(spongeElementFlags.begin(), spongeElementFlags.end(), false);

// count the number of elements in the sponge layer
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(MInt eId = 0; eId < m_elements->size(); eId++) {
    const MInt noNodes = m_elements->noNodesXD(eId);

    const MFloat* const nodeCoordinates = &m_elements->nodeCoordinates(eId);
    // check if any node of the element is inside the spongelayer
    for(MInt nodeId = 0; nodeId < noNodes; nodeId++) {
      const MFloat* const point = &nodeCoordinates[nodeId * nDim];

      // check for each boundary plane
      for(size_t i = 0; i < spongeLayers.size(); i++) {
        MFloat diff = 0.0;
        if(pointInsideSpongeLayer(spongeLayers[i], point, diff)) {
          spongeElementFlags[eId] = true;
          break;
        }
      }

      // leave loop since already determined element is inside of sponge layer
      if(spongeElementFlags[eId]) {
        break;
      }
    }
  }

  // Allocate memory for sponge elements
  const MInt spongeElementsCount = std::count(spongeElementFlags.begin(), spongeElementFlags.end(), true);
  m_spongeElements.maxPolyDeg(m_maxPolyDeg);
  m_spongeElements.reset(spongeElementsCount);

  // Determine sponge element id and create sponge element in collector
  for(size_t i = 0; i < spongeElementFlags.size(); i++) {
    if(spongeElementFlags[i]) {
      const MInt seId = noSpongeElements();
      m_spongeElements.append();
      m_spongeElements.elementId(seId) = i;
    }
  }

// 3. calculate eta for all sponge elements
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(MInt seId = 0; seId < noSpongeElements(); seId++) {
    calcSpongeEtaForAllNodes(seId, spongeLayers);
  }
}


/**
 * \brief Determines boundary plane.
 *
 * \author Hans Yu <hans.yu@rwth-aachen.de>
 * \date April 20th, 2015
 */
template <MInt nDim, class SysEqn>
void DgSponge<nDim, SysEqn>::initSpongeLayer(SpongeLayerType& spongeLayer, const BC& boundaryCondition) {
  TRACE();

  const MInt firstSrfcId = boundaryCondition->begin();
  const MFloat* srfcCoord = &m_surfaces->coords(firstSrfcId, 0);

  // Store orientation
  MInt orientation = m_surfaces->orientation(firstSrfcId);
  MInt direction = m_surfaces->internalSideId(firstSrfcId);
  spongeLayer.m_orientation = 2 * orientation + direction;

  // Correct orientation of the boundary condition manually
  // used for the direct-hybrid testcase of a acoustically
  // exicted flame
  if(spongeLayer.m_bcId == 200) {
    orientation = 0;
    direction = 0;
    spongeLayer.m_orientation = 0;
  }
  if(spongeLayer.m_bcId == 201) {
    orientation = 1;
    direction = 1;
    spongeLayer.m_orientation = 3;
  }
  if(spongeLayer.m_bcId == 202) {
    orientation = 1;
    direction = 1;
    spongeLayer.m_orientation = 3;
  }
  if(spongeLayer.m_bcId == 203) {
    orientation = 0;
    direction = 1;
    spongeLayer.m_orientation = 1;
  }

  // Store offset (distance to origin)
  spongeLayer.m_offset = srfcCoord[orientation];

  // Determine minimum and maximum coordinates
  std::fill_n(&spongeLayer.m_minCoord[0], nDim, std::numeric_limits<MFloat>::max());
  spongeLayer.m_minCoord[orientation] = srfcCoord[orientation];
  std::fill_n(&spongeLayer.m_maxCoord[0], nDim, -std::numeric_limits<MFloat>::max());
  spongeLayer.m_maxCoord[orientation] = srfcCoord[orientation];

  // Determine minimum and maximum extent of boundary by checking each surface
  for(MInt i = boundaryCondition->begin(); i < boundaryCondition->end(); i++) {
    const MFloat* srfcCenter = &m_surfaces->coords(i, 0);
    const MInt internalSideId = m_surfaces->internalSideId(i);
    // FIXME labels:DG,toenhance Boundary surfaces with two neighbor elements should not exist. For
    // now, they are simply skipped.
    if(internalSideId == -1) {
      continue;
    }
    const MInt eId = m_surfaces->nghbrElementIds(i, internalSideId);
    const MInt cellId = m_elements->cellId(eId);
    const MFloat halfCellLength = m_grid->halfCellLength(cellId);
    for(MInt j = 0; j < nDim; j++) {
      if(j != orientation) {
        spongeLayer.m_minCoord[j] = std::min(spongeLayer.m_minCoord[j], srfcCenter[j] - halfCellLength);
        spongeLayer.m_maxCoord[j] = std::max(spongeLayer.m_maxCoord[j], srfcCenter[j] + halfCellLength);
      }
    }
  }
}


/**
 * \brief Exchange sponge layers between all domains.
 *
 * \author Sven Berger, Hans Yu <hans.yu@rwth-aachen.de>
 * \date   April 21st, 2015
 */
template <MInt nDim, class SysEqn>
void DgSponge<nDim, SysEqn>::exchangeSpongeLayers(std::vector<SpongeLayerType>& spongeLayers) {
  TRACE();

  /////////////////////////////////////////////////////////////////////////////
  // Exchange sponge layer ids
  /////////////////////////////////////////////////////////////////////////////

  // Set up communication buffer
  const MInt noSpongeBcs = spongeLayers.size();
  ScratchSpace<MChar> whichSpongeLayers(noDomains(), noSpongeBcs, AT_, "whichSpongeLayers");

  // For each sponge BC id, store if it has surfaces on current domain
  for(MInt i = 0; i < noSpongeBcs; i++) {
    const MInt spongeBcId = spongeLayers[i].m_bcId;
    whichSpongeLayers(domainId(), i) =
        std::find_if(m_boundaryConditions->begin(),
                     m_boundaryConditions->end(),
                     [spongeBcId](const BC& bc) { return bc->id() == spongeBcId && bc->count() > 0; })
        != m_boundaryConditions->end();
  }

  // Exchange with all domains
  MPI_Allgather(MPI_IN_PLACE, noSpongeBcs, maia::type_traits<MChar>::mpiType(), &whichSpongeLayers(0, 0), noSpongeBcs,
                maia::type_traits<MChar>::mpiType(), mpiComm(), AT_, "MPI_IN_PLACE", "whichSpongeLayers(0");

  /////////////////////////////////////////////////////////////////////////////
  // Exchange boundary planes and thicknesses
  /////////////////////////////////////////////////////////////////////////////

  // One broadcast per sponge BC:
  // 1.) If whichSpongeLayerRecv is false, nothing has to be done.
  // 2.) If whichSpongeLayerRecv is true, make sure that domains with higher id
  // do not broadcast.
  // 3.) Prepare and perform broadcast.
  // Iterating through whichSpongeLayer is equivalent to iterating through
  // all domains and sponge BCs
  for(MInt i = 0; i < noDomains(); i++) {
    for(MInt j = 0; j < noSpongeBcs; j++) {
      // Step 1
      if(!whichSpongeLayers(i, j)) {
        continue;
      }

      // Step 2
      for(MInt k = i + 1; k < noDomains(); k++) {
        whichSpongeLayers(k, j) = false;
      }

      // Step 3
      const MInt spongeBcId = spongeLayers[j].m_bcId;
      SpongeLayerType& spongeLayer =
          *std::find_if(spongeLayers.begin(), spongeLayers.end(),
                        [spongeBcId](const SpongeLayerType& sl) { return sl.m_bcId == spongeBcId; });
      MPI_Bcast(&spongeLayer.m_orientation, 1, maia::type_traits<MInt>::mpiType(), i, mpiComm(), AT_,
                "spongeLayer.m_orientation");
      MPI_Bcast(&spongeLayer.m_offset, 1, maia::type_traits<MFloat>::mpiType(), i, mpiComm(), AT_,
                "spongeLayer.m_offset");
      MPI_Allreduce(MPI_IN_PLACE, &spongeLayer.m_minCoord[0], nDim, maia::type_traits<MFloat>::mpiType(), MPI_MIN,
                    mpiComm(), AT_, "MPI_IN_PLACE", "spongeLayer.m_minCoord[0]");
      MPI_Allreduce(MPI_IN_PLACE, &spongeLayer.m_maxCoord[0], nDim, maia::type_traits<MFloat>::mpiType(), MPI_MAX,
                    mpiComm(), AT_, "MPI_IN_PLACE", "spongeLayer.m_maxCoord[0]");
    }
  }
}


/**
 * \brief Return sponge element id for given element id.
 *        If none exists, return -1.
 *
 * \author Hans Yu <hans.yu@rwth-aachen.de>
 * \date April 2nd, 2015
 */
template <MInt nDim, class SysEqn>
MInt DgSponge<nDim, SysEqn>::spongeElementId(const MInt eId) const {
  const MInt* const begin = &m_spongeElements.elementId(0);
  const MInt* const end = &m_spongeElements.elementId(0) + noSpongeElements();
  const MInt* const low = std::lower_bound(begin, end, eId);

  // Return not found if std::lower_bound does not find anything
  if(low == end || *low != eId) {
    return -1;
  }

  // Otherwise return found sponge element id
  return std::distance(begin, low);
}


/**
 * \brief Calculates the sponge terms, that can be considered as an addition
 *source terms for each node and adds them to the time derivative of the
 *conservative variables.
 *
 * \author Vitali Pauz <v.pauz@aia.rwth-aachen.de>
 * \date April 4th, 2014
 */
template <MInt nDim, class SysEqn>
void DgSponge<nDim, SysEqn>::calcSourceTerms() {
  TRACE();

  // const MInt* polyDeg = &m_elements->polyDeg(0);
  const MInt* noNodes1D = &m_elements->noNodes1D(0);
  const MInt maxDataBlockSize = m_elements->maxNoNodesXD() * SysEqn::noVars();

  std::vector<MFloat> spongeSources(maxDataBlockSize);

#ifdef _OPENMP
#pragma omp parallel for firstprivate(spongeSources)
#endif
  for(MInt seId = 0; seId < noSpongeElements(); seId++) {
    const MInt eId = m_spongeElements.elementId(seId);
    const MInt dataBlockSize = m_elements->noNodesXD(eId) * SysEqn::noVars();
    m_sysEqn->calcSpongeSource(&m_elements->nodeVars(eId),
                               &m_elements->variables(eId),
                               noNodes1D[eId],
                               &m_spongeElements.spongeEta(seId),
                               &spongeSources[0]);
    MFloat* const rhs = &m_elements->rightHandSide(eId);
    for(MInt dataId = 0; dataId < dataBlockSize; dataId++) {
      rhs[dataId] += spongeSources[dataId];
    }
  }
}


/**
 * \brief Calculation of SpongeEta for all nodes of one DG-Element.
 *
 * \author Vitali Pauz <v.pauz@aia.rwth-aachen.de>
 * \date April 4th, 2014
 */
template <MInt nDim, class SysEqn>
void DgSponge<nDim, SysEqn>::calcSpongeEtaForAllNodes(const MInt seId,
                                                      const std::vector<SpongeLayerType>& spongeLayers) {
  TRACE();

  const MInt eId = m_spongeElements.elementId(seId);
  const MInt noNodes1D = m_elements->noNodes1D(eId);
  const MInt noNodes1D3 = (nDim == 3) ? noNodes1D : 1;

  MFloatTensor nodeCoordinates(&m_elements->nodeCoordinates(eId), noNodes1D, noNodes1D, noNodes1D3, nDim);

  MFloatTensor sEta(&m_spongeElements.spongeEta(seId), noNodes1D, noNodes1D, noNodes1D3);
  std::fill_n(&sEta(0, 0, 0), sEta.size(), 0.0);

  for(MInt i = 0; i < noNodes1D; i++) {
    for(MInt j = 0; j < noNodes1D; j++) {
      for(MInt k = 0; k < noNodes1D3; k++) {
        for(size_t s = 0; s < spongeLayers.size(); s++) {
          MFloat diff = 0.0;
          if(pointInsideSpongeLayer(spongeLayers[s], &nodeCoordinates(i, j, k, 0), diff)) {
            const MFloat eta = POW2(diff / spongeLayers[s].m_thickness);
            sEta(i, j, k) = std::max(sEta(i, j, k), eta);
          }
        }
      }
    }
  }
}


/**
 * \brief Calculates the distance between a point and a boundary plane. The
 *distance is signed, i.e., it is positive if the point is inside.
 *
 * \author Sven Berger, Hans Yu <hans.yu@rwth-aachen.de>
 * \date   April 21st, 2015
 */
template <MInt nDim, class SysEqn>
MFloat DgSponge<nDim, SysEqn>::distance(const SpongeLayerType& spongeLayer, const MFloat* point) const {
  // 0=x, 1=y, 2=z
  const MInt orientation = (spongeLayer.m_orientation - spongeLayer.m_orientation % 2) / 2;
  // 1=inside in positive direction, -1=inside in negative direction
  const MInt direction = 2 * (spongeLayer.m_orientation % 2) - 1;

  return direction * (point[orientation] - spongeLayer.m_offset);
}


/**
 * \brief Returns true if point is inside sponge layer.
 *        Calculates difference between thickness of sponge layer and distance
 *        of point to boundary plane.
 *
 * \author Hans Yu <hans.yu@rwth-aachen.de>
 * \date   April 21st, 2015
 */
template <MInt nDim, class SysEqn>
MBool DgSponge<nDim, SysEqn>::pointInsideSpongeLayer(const SpongeLayerType& spongeLayer, const MFloat* point,
                                                     MFloat& diff) const {
  // distance to boundary plane
  MFloat dist = distance(spongeLayer, point);
  if(dist < 0.0) {
    return false;
  }
  diff = spongeLayer.m_thickness - dist;
  if(diff < 0.0) {
    return false;
  }

  // within min and max coordinates
  for(MInt i = 0; i < nDim; i++) {
    const MInt orientation = (spongeLayer.m_orientation - spongeLayer.m_orientation % 2) / 2;
    if(orientation == i) {
      continue;
    }

    if(point[i] < spongeLayer.m_minCoord[i] || point[i] > spongeLayer.m_maxCoord[i]) {
      return false;
    }
  }

  return true;
}

#endif // DGSPONGE_H_
