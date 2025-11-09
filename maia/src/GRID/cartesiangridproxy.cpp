// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "cartesiangridproxy.h"

#include <algorithm>
#include <numeric>
#include "COMM/mpiexchange.h"
#include "COMM/mpioverride.h"
#include "IO/context.h"

using namespace std;


namespace maia {
namespace grid {

template <MInt nDim>
Proxy<nDim>::Proxy(const MInt solverId, Grid& grid_, Geometry<nDim>& geometry_)
  : m_solverId(solverId), m_grid(grid_), m_tree(solverId, grid_.treeb()), m_geometry(&geometry_) {
  if(raw().addSolverToGrid()) {
    MBool useSolverGeometry = false;
    useSolverGeometry = Context::getSolverProperty<MBool>("useSolverGeometry", solverId, AT_, &useSolverGeometry);
    if(useSolverGeometry && solverId == (raw().treeb().noSolvers() - 1)) {
      setSolverFlagsForAddedSolver();
    }
  }

  update();

  m_maxNoCells = raw().maxNoCells();

  if(!raw().paraViewPlugin()) {
    m_maxRefinementLevel = Context::getSolverProperty<MInt>("maxRfnmntLvl", solverId, AT_);

    if(m_maxRefinementLevel < maxLevel()) {
      TERMM(1, "Error: property maxRfnmntLvl is smaller than the maxLevel in the grid for solver #"
                   + std::to_string(solverId) + " (" + std::to_string(m_maxRefinementLevel) + " < "
                   + std::to_string(m_maxLevel) + ")");
    }
    m_maxNoCells = Context::getSolverProperty<MInt>("maxNoCells", solverId, AT_, &m_maxNoCells);
    if(Context::propertyExists("noDomains")) {
      const MInt testNoDomains = Context::getBasicProperty<MInt>("noDomains", AT_, 0);
      if(globalNoDomains() < testNoDomains) {
        // Here, the number of maxNoCells is scaled. This is useful if a test
        // case is specified to run with a certain number of ranks
        // ('noDomains' property in run.toml) but needs to be run on a lower
        // number of mpiranks, p.e. by running maia on an accelerator.
        m_maxNoCells *= testNoDomains / (MFloat)globalNoDomains();
        cerr0 << "noDomain > number of used mpi ranks! Therefore, increasing maxNoCells: " << m_maxNoCells << std::endl;
      }
    }
  }
}

template <MInt nDim>
Proxy<nDim>::Proxy(const MInt solverId, Grid& grid_)
  : m_solverId(solverId), m_grid(grid_), m_tree(solverId, grid_.treeb()) {
  update();

  if(!raw().paraViewPlugin()) {
    TERMM(1, "This constructor should only be called by the ParaView plugin");
  }
}

template <MInt nDim>
Proxy<nDim>::~Proxy() {
  TRACE();

  // Delete created communicator if existing and if it is not MPI_COMM_WORLD
  if(m_mpiComm != MPI_COMM_NULL && m_mpiComm != MPI_COMM_WORLD) {
    MPI_Comm_free(&m_mpiComm, AT_, "m_mpiComm");
  }
}


/// Update all solver-local information that is stored in the solver grid.
///
/// \author Michael Schlottke-Lakemper (mic) <mic@aia.rwth-aachen.de>
/// \date 2018-02-06
template <MInt nDim>
void Proxy<nDim>::update() {
  TRACE();

  // Update grid map
  updateGridMap();

  // Update information on parallelization including neighbor domains and window/halo cells
  updateParallelizationInfo();

  // Perform the remaining steps only if the solver is active on the current domain
  if(isActive()) {
    // Update general tree data
    updateTreeData();

    // Update tree-internal data
    // empty!
    m_tree.update(*this);

    // Update other grid-related information
    updateGridInfo();

    // update cutOff
    if(!raw().paraViewPlugin()) {
      updateCutOff();
    }
  }
}


template <MInt nDim>
void Proxy<nDim>::updateOther() {
  TRACE();

  // Update information on parallelization including neighbor domains and window/halo cells
  updateParallelizationInfo();

  // Perform the remaining steps only if the solver is active on the current domain
  if(isActive()) {
    // Update general tree data
    updateTreeData();

    // Update tree-internal data
    m_tree.update(*this);

    // Update other grid-related information
    updateGridInfo();

    // update cutOff
    updateCutOff();
  }
}

template <MInt nDim>
void Proxy<nDim>::resizeGridMap(const MInt solverSize) {
  TRACE();
  m_tree.m_solver2grid.resize(solverSize, -1);
  m_tree.m_grid2solver.resize(raw().treeb().size() + 1, -1); // +1 is to map "-1" to "-1"

  // ASSERT( m_tree.m_grid2solver[0] == -1, "" );
  m_tree.m_grid2solver[0] = -1;
}


template <MInt nDim>
void Proxy<nDim>::initGridMap() {
  TRACE();
  // Properly size maps
  if(g_multiSolverGrid) {
    m_tree.m_solver2grid.assign(raw().treeb().noNodesBySolver(solverId()), -1);
  } else {
    // FIXME/TODO labels:GRID Use proper solverId once multiple solvers can use a single grid
    m_tree.m_solver2grid.assign(raw().treeb().noNodesBySolver(0), -1);
  }
  m_tree.m_grid2solver.assign(raw().treeb().size() + 1, -1); // +1 is to map "-1" to "-1"
  m_tree.m_grid2solver[0] = -1;
}

template <MInt nDim>
void Proxy<nDim>::updateGridMap() {
  TRACE();

  // The grid is unaware of the solver geometry. Therefore, the solver flags
  // of the azimuthal halo cells need to be corrected in the proxy.
  if(azimuthalPeriodicity()) correctAzimuthalHaloCells();

  // Properly size maps
  if(g_multiSolverGrid) {
    m_tree.m_solver2grid.assign(raw().treeb().noNodesBySolver(solverId()), -1);
  } else {
    // FIXME/TODO labels:GRID Use proper solverId once multiple solvers can use a single grid
    m_tree.m_solver2grid.assign(raw().treeb().noNodesBySolver(0), -1);
  }
  m_tree.m_grid2solver.assign(raw().treeb().size() + 1, -1); // +1 is to map "-1" to "-1"

  // Create forward map (solver -> grid) (only if there are any cells on the present domain)
  if(m_tree.m_solver2grid.size() > 0) {
    if(g_multiSolverGrid) {
      raw().treeb().nodesBySolver(solverId(), m_tree.m_solver2grid.data());
    } else {
      // FIXME/TODO labels:GRID Use proper solverId once multiple solvers can use a single grid
      raw().treeb().nodesBySolver(0, m_tree.m_solver2grid.data());
    }
  }

  // Create reverse map (grid -> solver)
  for(MInt cellId = 0; cellId < static_cast<MInt>(m_tree.m_solver2grid.size()); cellId++) {
    const MInt gridCellId = m_tree.m_solver2grid[cellId];
    m_tree.m_grid2solver[gridCellId + 1] = cellId; // +1 is to map "-1" to "-1"
  }
}


template <MInt nDim>
void Proxy<nDim>::updateParallelizationInfo() {
  TRACE();

  // Update internal cell count
  m_noInternalCells = 0;
  for(MInt i = 0; i < raw().noInternalCells(); i++) {
    m_noInternalCells += (m_tree.grid2solver(i) > -1);
  }

  // Create new MPI communicator (but first delete existing communicator if present)
  if(m_mpiComm != MPI_COMM_NULL && m_mpiComm != MPI_COMM_WORLD) {
    MPI_Comm_free(&m_mpiComm, AT_, "m_mpiComm");
  }

  // Determine number of internal cells on all domain
  MIntScratchSpace noInternalCells_(raw().noDomains(), AT_, "noInternalCells_");
  noInternalCells_[raw().domainId()] = noInternalCells();
  MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, &noInternalCells_[0], 1, type_traits<MInt>::mpiType(),
                raw().mpiComm(), AT_, "MPI_IN_PLACE", "noInternalCells_[0]");

  // Determine list of subdomains with cells
  MIntScratchSpace domains(raw().noDomains(), AT_, "domains");
  MInt noDomains_ = 0;
  for(MInt d = 0; d < raw().noDomains(); d++) {
    if(noInternalCells_[d] > 0) {
      domains[noDomains_] = d;
      noDomains_++;
    }
  }

  // Obtain new group
  MPI_Group globalGroup, localGroup;
  MPI_Comm_group(raw().mpiComm(), &globalGroup, AT_, "globalGroup");

  // @ansgar_largeJob debug stack memory usage
  // writeMemoryStatistics(MPI_COMM_WORLD, globalNoDomains(), globalDomainId(), AT_, "grid proxy before
  // MPI_Group_incl");

  // NOTE: allocates large memory chunk on stack on O(100000) ranks (at least on Hawk, 12/21) and scales linear!
  MPI_Group_incl(globalGroup, noDomains_, &domains[0], &localGroup, AT_);

  // @ansgar_largeJob debug stack memory usage
  // writeMemoryStatistics(MPI_COMM_WORLD, globalNoDomains(), globalDomainId(), AT_, "grid proxy after MPI_Group_incl");

  // Create new communicator and clean up
  MPI_Comm_create(raw().mpiComm(), localGroup, &m_mpiComm, AT_, "m_mpiComm");

  MPI_Group_free(&globalGroup, AT_);
  MPI_Group_free(&localGroup, AT_);

  if(noDomains_ == raw().noDomains()) {
    m_hasInactiveRanks = false;
  } else {
    m_hasInactiveRanks = true;
  }

  // If solver is active on current domain, update domain info, otherwise return early.
  if(noInternalCells() > 0) {
    // Store new domain id and number of domains in communicator
    MPI_Comm_rank(mpiComm(), &m_domainId);
    MPI_Comm_size(mpiComm(), &m_noDomains);

    if(noDomains() != noDomains_) {
      TERMM(1, "an error occurred while creating the solver-specific MPI communicator");
    }

#ifndef NDEBUG
    if(g_multiSolverGrid && domainId() == 0 && noDomains() < 50) {
      cout << "Solver " << solverId() << ": Global domain " << globalDomainId()
           << " is the root domain of a new solver-specific communicator with the following " << noDomains_
           << " domain(s): " << domains[0];
      for(MInt i = 1; i < noDomains_; i++) {
        cout << ", " << domains[i];
      }
      cout << endl;
    }
#endif

    // Mark domain as active
    m_isActive = true;
  } else {
    // If there are no local internal cells, mark this domain as inactive and reset domain info
    m_isActive = false;
    m_domainId = -1;
    m_noDomains = -1;
    m_noInternalCells = -1;
    std::vector<MLong>().swap(m_domainOffsets);
    std::vector<MInt>().swap(m_neighborDomains);
    std::vector<MInt>().swap(m_neighborDomainIndex);

    std::vector<std::vector<MInt>>().swap(m_windowCells);
    std::vector<std::vector<MInt>>().swap(m_haloCells);
    std::map<MInt, MInt>().swap(m_global2solver);
    m_maxLevel = -1;

    return;
  }

  // Determine domain offsets (last entry contains global number of cells)
  m_domainOffsets.assign(noDomains() + 1, -1);
  // m_domainOffsets[0] = 0;
  m_domainOffsets[0] = bitOffset();
  for(MInt d = 1; d < noDomains() + 1; d++) {
    m_domainOffsets[d] = m_domainOffsets[d - 1] + (MLong)noInternalCells_[domains[d - 1]];
  }

#ifndef NDEBUG
  // Sanity check: domain offsets are the same on all domains
  checkOffsetConsistency();
#endif

  // Create map from global to solver domain ids
  std::map<MInt, MInt>().swap(m_global2solver); // Clear map first
  for(MInt d = 0; d < noDomains(); d++) {
    m_global2solver[domains[d]] = d;
  }

  if(!g_multiSolverGrid) {
    for(MInt d = 0; d < noDomains(); d++) {
      ASSERT(m_global2solver[d] == d, "");
    }
  }

  // Determine neighbor domains
  // Map from old position to new position
  std::map<MInt, MInt> neighborDomains;
  for(MInt d = 0; d < raw().noNeighborDomains(); d++) {
    // Skip domain if it is not among solver-specific domains
    if(m_global2solver.count(raw().neighborDomain(d)) == 0) {
      continue;
    }

    // Check in window cells
    for(MInt c = 0; c < raw().noWindowCells(d); c++) {
      if((raw().haloMode() == 0 && m_tree.grid2solver(raw().windowCell(d, c)) > -1)
         || (raw().haloMode() > 0 && m_tree.grid2solver(raw().windowCell(d, c)) > -1
             && raw().isSolverWindowCell(d, raw().windowCell(d, c), solverId()))) { // WH_old
        const MInt position = neighborDomains.size();
        neighborDomains[d] = position;
        break;
      }
    }

    // Continue early if domain already found
    if(neighborDomains.count(d) > 0) {
      continue;
    }

    // Check in halo cells
    for(MInt c = 0; c < raw().noHaloCells(d); c++) {
      if(m_tree.grid2solver(raw().haloCell(d, c)) > -1) {
        const MInt position = neighborDomains.size();
        neighborDomains[d] = position;
        break;
      }
    }
  }
  m_neighborDomains.assign(neighborDomains.size(), -1);
  m_neighborDomainIndex.assign(noDomains(), -1);
  for(auto&& m : neighborDomains) {
    // m.second holds new position, m.first holds position in the global neighbor domain arrays
    m_neighborDomains[m.second] = raw().neighborDomain(m.first);
  }

  // Update list of neighbor domain ids using a map from global domain id to solver domain id
  for(MInt d = 0; d < noNeighborDomains(); d++) {
    m_neighborDomains[d] = m_global2solver[m_neighborDomains[d]];
    m_neighborDomainIndex[m_neighborDomains[d]] = d;
  }

#ifndef NDEBUG
  // Sanity check: neighbor domains are consistent on all domains
  checkNeighborConsistency();
#endif

  setupWindowHaloConnectivityOnLeafLvl(neighborDomains);

#ifndef NDEBUG
  // Sanity check: window and halo cells are consistent on all domains
  checkWindowHaloConsistency();
#endif


  // Setup azimuthal periodic connection
  m_azimuthalNeighborDomains.clear();
  if(azimuthalPeriodicity()) {
    neighborDomains.clear();

    MInt windowCnt = 0;
    MInt haloCnt = 0;
    MInt tmpCnt = 0;
    for(MInt d = 0; d < raw().noAzimuthalNeighborDomains(); d++) {
      // Skip domain if it is not among solver-specific domains
      if(m_global2solver.count(raw().azimuthalNeighborDomain(d)) == 0) {
        windowCnt += raw().noAzimuthalWindowCells(d);
        haloCnt += raw().noAzimuthalHaloCells(d);
        continue;
      }

      // Check in window cells
      tmpCnt = windowCnt;
      for(MInt c = 0; c < raw().noAzimuthalWindowCells(d); c++) {
        if(m_tree.grid2solver(raw().azimuthalWindowCell(d, c)) > -1 && m_isOutsideWindow[tmpCnt] == false) {
          const MInt position = neighborDomains.size();
          neighborDomains[d] = position;
          break;
        }
        tmpCnt++;
      }
      windowCnt += raw().noAzimuthalWindowCells(d);

      // Continue early if domain already found
      if(neighborDomains.count(d) > 0) {
        haloCnt += raw().noAzimuthalHaloCells(d);
        continue;
      }

      // Check in halo cells
      tmpCnt = haloCnt;
      for(MInt c = 0; c < raw().noAzimuthalHaloCells(d); c++) {
        if(m_tree.grid2solver(raw().azimuthalHaloCell(d, c)) > -1 && m_isOutsideHalo[tmpCnt] == false) {
          const MInt position = neighborDomains.size();
          neighborDomains[d] = position;
          break;
        }
        tmpCnt++;
      }
      haloCnt += raw().noAzimuthalHaloCells(d);
    }
    m_azimuthalNeighborDomains.assign(neighborDomains.size(), -1);
    m_azimuthalNeighborDomainIndex.assign(noDomains(), -1);
    for(auto&& m : neighborDomains) {
      // m.second holds new position, m.first holds position in the global neighbor domain arrays
      m_azimuthalNeighborDomains[m.second] = raw().azimuthalNeighborDomain(m.first);
    }

    // Update list of neighbor domain ids using a map from global domain id to solver domain id
    for(MInt d = 0; d < noAzimuthalNeighborDomains(); d++) {
      m_azimuthalNeighborDomains[d] = m_global2solver[m_azimuthalNeighborDomains[d]];
      m_azimuthalNeighborDomainIndex[m_azimuthalNeighborDomains[d]] = d;
    }

#ifndef NDEBUG
    // Sanity check: neighbor domains are consistent on all domains
    checkNeighborConsistencyAzimuthal();
#endif

    m_azimuthalWindowCells.clear();
    m_azimuthalWindowCells.resize(noAzimuthalNeighborDomains());
    m_azimuthalHaloCells.clear();
    m_azimuthalHaloCells.resize(noAzimuthalNeighborDomains());
    m_azimuthalUnmappedHaloCells.clear();
    m_azimuthalUnmappedHaloDomains.clear();

    windowCnt = 0;
    haloCnt = 0;
    tmpCnt = 0;
    for(MInt d = 0; d < raw().noAzimuthalNeighborDomains(); d++) {
      if(m_global2solver.count(raw().azimuthalNeighborDomain(d)) == 0) {
        windowCnt += raw().noAzimuthalWindowCells(d);
        haloCnt += raw().noAzimuthalHaloCells(d);
        continue;
      }

      // Store old and new position for convenience
      const MInt oldDomainPosition = d;
      const MInt newDomainPosition = m_azimuthalNeighborDomainIndex[m_global2solver[raw().azimuthalNeighborDomain(d)]];
      if(newDomainPosition < 0) {
        windowCnt += raw().noAzimuthalWindowCells(d);
        haloCnt += raw().noAzimuthalHaloCells(d);
        continue;
      }

      // Determine window cells
      tmpCnt = windowCnt;
      std::vector<MInt> windowCells;
      windowCells.reserve(raw().noAzimuthalWindowCells(oldDomainPosition));
      for(MInt c = 0; c < raw().noAzimuthalWindowCells(oldDomainPosition); c++) {
        const MInt windowCellId = m_tree.grid2solver(raw().azimuthalWindowCell(oldDomainPosition, c));
        if(windowCellId > -1 && m_isOutsideWindow[tmpCnt] == false) {
          windowCells.push_back(windowCellId);
        }
        tmpCnt++;
      }
      m_azimuthalWindowCells[newDomainPosition] = windowCells;
      windowCnt += raw().noAzimuthalWindowCells(d);

      // Determine halo cells
      tmpCnt = haloCnt;
      std::vector<MInt> haloCells;
      haloCells.reserve(raw().noAzimuthalHaloCells(oldDomainPosition));
      for(MInt c = 0; c < raw().noAzimuthalHaloCells(oldDomainPosition); c++) {
        const MInt haloCellId = m_tree.grid2solver(raw().azimuthalHaloCell(oldDomainPosition, c));
        if(haloCellId > -1) {
          if(m_isOutsideHalo[tmpCnt] == false) {
            haloCells.push_back(haloCellId);
          } else {
            m_azimuthalUnmappedHaloCells.push_back(haloCellId);
            m_azimuthalUnmappedHaloDomains.push_back(m_global2solver[raw().azimuthalNeighborDomain(d)]);
          }
        }
        tmpCnt++;
      }
      m_azimuthalHaloCells[newDomainPosition] = haloCells;
      haloCnt += raw().noAzimuthalHaloCells(d);
    }

    for(MInt c = 0; c < raw().noAzimuthalUnmappedHaloCells(); c++) {
      const MInt haloCellId = m_tree.grid2solver(raw().azimuthalUnmappedHaloCell(c));
      if(haloCellId > -1) {
        m_azimuthalUnmappedHaloCells.push_back(haloCellId);
        m_azimuthalUnmappedHaloDomains.push_back(m_global2solver[raw().azimuthalUnmappedHaloDomain(c)]);
      }
    }

#ifndef NDEBUG
    // Sanity check: window and halo cells are consistent on all domains
    checkWindowHaloConsistencyAzimuthal();
#endif
  }

  // The following output may be helpful until the multi-solver capabilities are sufficiently tested
  // stringstream ss;
  // ss << "---------------------------------------------------------------------------------------"
  //    << endl;
  // ss << solverId() << " " << domainId() << " Window/halo information" << endl;
  // ss << solverId() << " " << domainId() << " =======================" << endl;
  // ss << solverId() << " " << domainId() << " noNeighborDomains(): " << noNeighborDomains()
  //    << endl;
  // ss << solverId() << " " << domainId() << " neighbor domains:";
  // for (MInt i = 0; i < noNeighborDomains(); i++) {
  //   ss << " " << neighborDomain(i);
  // }
  // ss << endl;
  // for (MInt i = 0; i < noNeighborDomains(); i++) {
  //   // Window cells
  //   ss << solverId() << " " << domainId() << " has " << noWindowCells(i)
  //       << " window cells for neighbor " << i << " (domain id: " << neighborDomain(i) << "):";
  //   for (MInt j = 0; j < noWindowCells(i); j++) {
  //     ss << " " << raw().treeb().globalId(m_tree.solver2grid(windowCell(i, j)));
  //   }
  //   ss << endl;

  //   // Halo cells
  //   ss << solverId() << " " << domainId() << " has " << noHaloCells(i)
  //       << " halo   cells for neighbor " << i << " (domain id: " << neighborDomain(i) << "):";
  //   for (MInt j = 0; j < noHaloCells(i); j++) {
  //     ss << " " << raw().treeb().globalId(m_tree.solver2grid(haloCell(i, j)));
  //   }
  //   ss << endl;
  // }
  // ss << "---------------------------------------------------------------------------------------"
  //    << endl;
  // cout << ss.str() << endl;
  // MInt totalHaloCells = 0;
  // for (MInt i = 0; i < noNeighborDomains(); i++)
  //   totalHaloCells += noHaloCells(i);
  // MPI_Allreduce(MPI_IN_PLACE, &totalHaloCells, 1, type_traits<MInt>::mpiType(), MPI_SUM, mpiComm(), AT_,
  //     "MPI_IN_PLACE", "totalHaloCells");
  // if (domainId()==0)
  //   cout << " totalNoHaloCells of solver=" << solverId() << ": " << totalHaloCells << endl;
}


template <MInt nDim>
void Proxy<nDim>::updateLeafCellExchange() {
  TRACE();

  // reset info
  std::vector<std::vector<MInt>>().swap(m_leafWindowCells);
  m_leafWindowCells.resize(noNeighborDomains());
  std::vector<std::vector<MInt>>().swap(m_leafHaloCells);
  m_leafHaloCells.resize(noNeighborDomains());

  m_leafRecvSize = 0;
  m_leafSendSize = 0;

  std::vector<MInt>().swap(m_leafSendNeighborDomains);
  std::vector<MInt>().swap(m_leafRecvNeighborDomains);

  if(noNeighborDomains() == 0) {
    return;
  }
  // meaning only active-rangs for parallel computation

  ScratchSpace<MPI_Request> sendReq(noNeighborDomains(), AT_, "sendReq");
  ScratchSpace<MPI_Request> recvReq(noNeighborDomains(), AT_, "recvReq");
  sendReq.fill(MPI_REQUEST_NULL);
  recvReq.fill(MPI_REQUEST_NULL);

  MInt noHaloSums = 0;
  MInt noWindowSums = 0;
  for(MInt d = 0; d < noNeighborDomains(); d++) {
    noHaloSums += noHaloCells(d);
    noWindowSums += noWindowCells(d);
  }

  MIntScratchSpace sendBuffer(noHaloSums, AT_, "sendBuffer");
  sendBuffer.fill(-1);
  MIntScratchSpace receiveBuffer(noWindowSums, AT_, "receiveBuffer");
  receiveBuffer.fill(-1);

  // set leaf-cell property for halos
  MInt haloId = 0;
  for(MInt i = 0; i < noNeighborDomains(); i++) {
    std::vector<MInt> leafHaloCells;
    leafHaloCells.reserve(noHaloCells(i));

    for(MInt j = 0; j < noHaloCells(i); j++) {
      const MInt cellId = haloCell(i, j);
      if(!tree().isLeafCell(cellId)) {
        haloId++;
        continue;
      }
      leafHaloCells.push_back(cellId);
      sendBuffer[haloId++] = cellId;
    }
    m_leafHaloCells[i] = leafHaloCells;
    if(!leafHaloCells.empty()) {
      m_leafRecvNeighborDomains.push_back(i);
    }
  }

  // sending from halo -> window
  MInt sendOffset = 0;
  MInt sendCnt = 0;
  for(MInt i = 0; i < noNeighborDomains(); i++) {
    const MInt bufSize = noHaloCells(i);
    if(bufSize == 0) {
      continue;
    }
    MPI_Issend(&sendBuffer[sendOffset], bufSize, MPI_INT, neighborDomain(i), 0, mpiComm(), &sendReq[i], AT_,
               "sendBuffer[sendOffset]");
    sendOffset += bufSize;
    sendCnt++;
  }

  MInt receiveOffset = 0;
  for(MInt i = 0; i < noNeighborDomains(); i++) {
    const MInt bufSize = noWindowCells(i);
    if(bufSize == 0) {
      continue;
    }
    MPI_Recv(&receiveBuffer[receiveOffset], bufSize, MPI_INT, neighborDomain(i), 0, mpiComm(), MPI_STATUS_IGNORE, AT_,
             "receiveBuffer[receiveOffset]");
    receiveOffset += bufSize;
  }
  if(sendCnt > 0) {
    MPI_Waitall(noNeighborDomains(), &sendReq[0], MPI_STATUSES_IGNORE, AT_);
  }

  // write window cell information
  MInt windowId = 0;
  for(MInt i = 0; i < noNeighborDomains(); i++) {
    std::vector<MInt> leafWindowCells;
    leafWindowCells.reserve(noWindowCells(i));

    for(MInt j = 0; j < noWindowCells(i); j++) {
      if(receiveBuffer[windowId++] != -1) {
        ASSERT(windowCell(i, j) > -1, "");
        leafWindowCells.push_back(windowCell(i, j));
      }
    }
    m_leafWindowCells[i] = leafWindowCells;
    if(!leafWindowCells.empty()) {
      m_leafSendNeighborDomains.push_back(i);
    }
  }

  // reduce neighbor-domains for leaf-cell communication
  for(MInt n = 0; n < noLeafRecvNeighborDomains(); n++) {
    const MInt d = leafRecvNeighborDomain(n);
    m_leafRecvSize += noLeafHaloCells(d);
  }

  for(MInt n = 0; n < noLeafSendNeighborDomains(); n++) {
    const MInt d = leafSendNeighborDomain(n);
    m_leafSendSize += noLeafWindowCells(d);
  }
}


/// Update tree-related data.
///
/// \author Michael Schlottke-Lakemper (mic) <mic@aia.rwth-aachen.de>
/// \date 2018-04-13
template <MInt nDim>
void Proxy<nDim>::updateTreeData() {
  TRACE();

  // Determine solver-specific global ids for internal cells
  m_tree.m_globalIds.assign(noCells(), -1);

  MInt localCnt = 0;
  const MInt firstGridMinCell = raw().minCell(0);

  // Partition level shift: check if the first min-cell is a partition level ancestor halo min-cell
  // (i.e. located on another domain and the first local partition cell corresponding to the
  // domain-offset  is an offspring of this min-cell) (similar to
  // cartesiangrid::computeGlobalIds())
  if(raw().a_hasProperty(firstGridMinCell, Cell::IsPartLvlAncestor)
     && raw().a_hasProperty(firstGridMinCell, Cell::IsHalo)) {
    descendStoreGlobalId(firstGridMinCell, localCnt);
  }

  // Loop over all min-level cells, start with min-cell at position 1 if there is a partition level
  // shift handled by the case above
  for(MInt i = std::max(0, std::min(1, localCnt)); i < raw().noMinCells(); i++) {
    const MInt gridCellId = raw().minCell(i);

    // Skip min-level cells not belonging to this solver or if its a halo
    if(!solverFlag(gridCellId, solverId())) continue;
    if(raw().a_hasProperty(gridCellId, Cell::IsHalo)) continue;
    descendStoreGlobalId(gridCellId, localCnt);
  }

  // Determine solver-specific global ids for halo cells
  if(m_neighborDomains.size() > 0) {
    mpi::exchangeData(m_neighborDomains, m_haloCells, m_windowCells, mpiComm(), m_tree.m_globalIds.data());
  }
  if(m_azimuthalNeighborDomains.size() > 0) {
    mpi::exchangeData(m_azimuthalNeighborDomains, m_azimuthalHaloCells, m_azimuthalWindowCells, mpiComm(),
                      m_tree.m_globalIds.data());
  }

  // If this is a single solver case the determined global ids have to match the ones in the grid
  // Note: this check was not suitable for a targetGridMinLevel != minLevel because the min-level
  // cell hilbert ids were not determined correctly such that the ordering was wrong
  if(!g_multiSolverGrid) { // && raw().targetGridMinLevel() == minLevel()) {
    for(MInt i = 0; i < noInternalCells(); i++) {
      ASSERT(m_tree.m_globalIds[i] == raw().a_globalId(i),
             std::to_string(domainId()) + " globalId mismatch: " + std::to_string(i) + ", "
                 + std::to_string(m_tree.m_globalIds[i]) + " != " + std::to_string(raw().a_globalId(i)));
    }
  }

  m_tree.createGlobalToLocalIdMapping();
}


/// Update tree-related data.
///
/// \author Tim Wegmann
/// \date 2018-10-09
template <MInt nDim>
void Proxy<nDim>::descendStoreGlobalId(MInt gridCellId, MInt& localCnt) {
  if(!solverFlag(gridCellId, solverId())) return;

  // Update global id (unless it is a halo cell)
  if(!raw().a_hasProperty(gridCellId, Cell::IsHalo)) {
    ASSERT(solverFlag(gridCellId, solverId()), "");
    ASSERT(localCnt < noInternalCells(), "");

    const MInt solverCellId = m_tree.grid2solver(gridCellId);
    m_tree.m_globalIds[solverCellId] = m_domainOffsets[domainId()] + localCnt++;

    ASSERT(m_tree.m_globalIds[solverCellId] >= m_domainOffsets[domainId()]
               && m_tree.m_globalIds[solverCellId] < m_domainOffsets[domainId() + 1]
               && m_tree.m_globalIds[solverCellId] < m_domainOffsets[domainId()] + noInternalCells(),
           "");
  }

  // Descend tree to all children that belong to this solver
  for(MInt child = 0; child < ipow(2, nDim); child++) {
    if(raw().a_childId(gridCellId, child) < 0) continue;
    if(!solverFlag(gridCellId, solverId())) continue;
    descendStoreGlobalId(raw().a_childId(gridCellId, child), localCnt);
  }
}

/** \brief Find neighbor domain id corresponding to given solver-specific globalId
    \author Tim Wegmann
    \date October 2018
  */
template <MInt nDim>
MInt Proxy<nDim>::findNeighborDomainId(const MLong globalId) {
  const MInt domain = raw().findNeighborDomainId(globalId, noDomains(), &m_domainOffsets[0]);
  if(!g_multiSolverGrid) {
    ASSERT(domain == raw().findNeighborDomainId(globalId), "");
  }
  return domain;
}


/// Update general grid information that may rely on an up-to-date tree.
///
/// \author Michael Schlottke-Lakemper (mic) <mic@aia.rwth-aachen.de>
/// \date 2018-04-13
template <MInt nDim>
void Proxy<nDim>::updateGridInfo() {
  TRACE();

  // Update maxLevel
  m_maxLevel = -1;

  for(MInt i = 0; i < noInternalCells(); i++) {
    m_maxLevel = max(m_maxLevel, m_tree.level(i));
  }
  MPI_Allreduce(MPI_IN_PLACE, &m_maxLevel, 1, type_traits<MInt>::mpiType(), MPI_MAX, mpiComm(), AT_, "MPI_IN_PLACE",
                "m_maxLevel");

  if(minLevel() < 0 || minLevel() > m_maxLevel) {
    mTerm(1, AT_, "Inconsistent min/max levels: " + to_string(minLevel()) + "/" + to_string(m_maxLevel));
  }
}


/// Check domain offsets for consistency.
///
/// \author Michael Schlottke-Lakemper (mic) <mic@aia.rwth-aachen.de>
/// \date 2018-03-15
///
/// This method checks if the following information matches on corresponding domains:
/// - domain offset on domain 0 is 0
/// - offsets are strictly monotonically increasing
/// - offsets are the same on all domains
/// Note: m_neighborDomains must already contain solver-local domain ids!
template <MInt nDim>
void Proxy<nDim>::checkOffsetConsistency() const {
  TRACE();

  // Check if first offset is zero
  if(m_domainOffsets[0] != bitOffset()) {
    TERMM(1, "Domain offset for domain 0 must be equal to the 32-Bit-Offset (0)!");
  }

  // Check if offsets are strictly monotonically increasing
  for(MInt i = 1; i < noDomains() + 1; i++) {
    if(m_domainOffsets[i] <= m_domainOffsets[i - 1]) {
      TERMM(1, "Domain offsets are not strictly monotonically increasing");
    }
  }

  // Check if offsets are the same on all domains
  MLongScratchSpace rootOffsets(noDomains() + 1, AT_, "rootOffsets");
  copy(m_domainOffsets.begin(), m_domainOffsets.end(), rootOffsets.begin());
  MPI_Bcast(&rootOffsets[0], noDomains() + 1, type_traits<MLong>::mpiType(), 0, mpiComm(), AT_, "rootOffsets[0]");
  for(MInt i = 0; i < noDomains() + 1; i++) {
    if(m_domainOffsets[i] != rootOffsets[i]) {
      TERMM(1, "Domain offsets differ from root domain");
    }
  }

  if(!g_multiSolverGrid) {
    for(MInt i = 1; i < noDomains() + 1; i++) {
      ASSERT(m_domainOffsets[i] == raw().domainOffset(i), "");
    }
  }
}


/// Check neighbor domains for consistency.
///
/// \author Michael Schlottke-Lakemper (mic) <mic@aia.rwth-aachen.de>
/// \date 2018-03-15
///
/// This method checks if the following information matches on corresponding domains:
/// - neighbor domains are sorted ascendingly
/// - domains are neighbors
/// Note: m_neighborDomains must already contain solver-local domain ids!
template <MInt nDim>
void Proxy<nDim>::checkNeighborConsistency() const {
  TRACE();

  // Check if neighbor domains are sorted ascendingly
  if(!std::is_sorted(m_neighborDomains.begin(), m_neighborDomains.end())) {
    TERMM(1, "Neighbor domains are not sorted in ascending order.");
  }

  // Determine all neighbor domains of current domain
  MCharScratchSpace isNeighborDomain(noDomains(), AT_, "isNeighborDomain");
  fill(isNeighborDomain.begin(), isNeighborDomain.end(), 0);
  for(MInt d = 0; d < noNeighborDomains(); d++) {
    isNeighborDomain[neighborDomain(d)] = 1;
  }

  // Exchange with all domains in solver-local communicator
  MPI_Alltoall(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, &isNeighborDomain[0], 1, type_traits<MChar>::mpiType(), mpiComm(),
               AT_, "MPI_IN_PLACE", "isNeighborDomain[0]");

  // Check if all domains that the current domain considers as neighbors also reciprocate
  for(MInt d = 0; d < noNeighborDomains(); d++) {
    if(!isNeighborDomain[neighborDomain(d)]) {
      TERMM(1, "Domain " + to_string(domainId()) + " has domain " + to_string(neighborDomain(d))
                   + " as a neighbor but not the other way around.");
    }
  }
}


/// Create window/halo connectivity on leaf cells.
///
/// \author XXX
/// \date 2021-06-23
///
template <MInt nDim>
void Proxy<nDim>::setupWindowHaloConnectivityOnLeafLvl(std::map<MInt, MInt>& neighborDomains) {
  TRACE();

  std::vector<std::vector<MInt>>().swap(m_windowCells);
  std::vector<std::vector<MInt>>().swap(m_haloCells);

  // WH_old
  if(raw().haloMode() > 0) {
    // Determine window/halo cells
    m_windowCells.resize(noNeighborDomains());
    m_haloCells.resize(noNeighborDomains());
    for(auto&& m : neighborDomains) {
      // Store old and new position for convenience
      const MInt oldDomainPosition = m.first;
      const MInt newDomainPosition = m.second;


      // Determine window cells
      std::vector<MInt> windowCells;
      windowCells.reserve(raw().noWindowCells(oldDomainPosition));
      for(MInt c = 0; c < raw().noWindowCells(oldDomainPosition); c++) {
        const MInt windowCellId = m_tree.grid2solver(raw().windowCell(oldDomainPosition, c));
        if(windowCellId > -1) {
          ASSERT(m_tree.solver2grid(windowCellId) == raw().windowCell(oldDomainPosition, c),
                 std::to_string(raw().windowCell(oldDomainPosition, c)) + "|" + std::to_string(windowCellId) + "|"
                     + std::to_string(m_tree.solver2grid(windowCellId)));
          if(raw().isSolverWindowCell(oldDomainPosition, raw().windowCell(oldDomainPosition, c), solverId())) {
            windowCells.push_back(windowCellId);
          }
        }
      }
      m_windowCells[newDomainPosition] = windowCells;
      // Check won't work in case of periodicity
      // TERMM_IF_NOT_COND((signed)m_windowCells[newDomainPosition].size()==raw().noSolverWindowCells(oldDomainPosition,m_solverId),
      // std::to_string((signed)m_windowCells[newDomainPosition].size()) + " " +
      // std::to_string(raw().noSolverWindowCells(oldDomainPosition,m_solverId)));

      // Determine halo cells
      std::vector<MInt> haloCells;
      haloCells.reserve(raw().noHaloCells(oldDomainPosition));
      for(MInt c = 0; c < raw().noHaloCells(oldDomainPosition); c++) {
        const MInt haloCellId = m_tree.grid2solver(raw().haloCell(oldDomainPosition, c));
        ASSERT((haloCellId > -1) == raw().a_solver(raw().haloCell(oldDomainPosition, c), m_solverId),
               "Shut down your PC!");
        // TODO_SS labels:GRID the last check is currently necessary because of tagActiveWindowsOnLeafLvl
        if(haloCellId > -1 && raw().a_solver(raw().haloCell(oldDomainPosition, c), m_solverId)) {
          ASSERT(m_tree.solver2grid(haloCellId) == raw().haloCell(oldDomainPosition, c), "");
          haloCells.push_back(haloCellId);
        }
      }
      m_haloCells[newDomainPosition] = haloCells;
      ASSERT((signed)m_haloCells[newDomainPosition].size()
                 == raw().noSolverHaloCells(oldDomainPosition, m_solverId, true),
             std::to_string(m_haloCells[newDomainPosition].size()) + " "
                 + std::to_string(raw().noSolverHaloCells(oldDomainPosition, m_solverId, true)));
    }
  } else {
    // Determine window/halo cells
    m_windowCells.resize(noNeighborDomains());
    m_haloCells.resize(noNeighborDomains());
    for(auto&& m : neighborDomains) {
      // Store old and new position for convenience
      const MInt oldDomainPosition = m.first;
      const MInt newDomainPosition = m.second;


      // Determine window cells
      std::vector<MInt> windowCells;
      windowCells.reserve(raw().noWindowCells(oldDomainPosition));
      for(MInt c = 0; c < raw().noWindowCells(oldDomainPosition); c++) {
        const MInt windowCellId = m_tree.grid2solver(raw().windowCell(oldDomainPosition, c));
        if(windowCellId > -1) {
          windowCells.push_back(windowCellId);
        }
      }
      m_windowCells[newDomainPosition] = windowCells;

      // Determine halo cells
      std::vector<MInt> haloCells;
      haloCells.reserve(raw().noHaloCells(oldDomainPosition));
      for(MInt c = 0; c < raw().noHaloCells(oldDomainPosition); c++) {
        const MInt haloCellId = m_tree.grid2solver(raw().haloCell(oldDomainPosition, c));
        if(haloCellId > -1) {
          haloCells.push_back(haloCellId);
        }
      }
      m_haloCells[newDomainPosition] = haloCells;
    }
  }
}


/// Check window/halo cells for consistency.
///
/// \author Michael Schlottke-Lakemper (mic) <mic@aia.rwth-aachen.de>
/// \date 2018-03-15
///
/// This method checks if the following information matches on corresponding domains:
/// - for each pair of neighbor domains: number of window/halo cells
/// - for each pair of neighbor domains: grid global ids of window/halo cells
template <MInt nDim>
void Proxy<nDim>::checkWindowHaloConsistency() const {
  TRACE();
  if(globalNoDomains() == 1) {
    return;
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Check #1: number of window/halo cells match
  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Start receiving number of window cells from each neighbor domain
  ScratchSpace<MPI_Request> recvRequests(max(1, noNeighborDomains()), AT_, "recvRequests");
  MIntScratchSpace noWindowCellsRecv(max(1, noNeighborDomains()), AT_, "noWindowCellsRecv");
  for(MInt d = 0; d < noNeighborDomains(); d++) {
    noWindowCellsRecv[d] = -1;
    MPI_Irecv(&noWindowCellsRecv[d], 1, type_traits<MInt>::mpiType(), neighborDomain(d), neighborDomain(d), mpiComm(),
              &recvRequests[d], AT_, "noWindowCellsRecv[d]");
  }

  // Start sending number of window cells to each neighbor domain
  ScratchSpace<MPI_Request> sendRequests(max(1, noNeighborDomains()), AT_, "sendRequests");
  MIntScratchSpace noWindowCellsSend(max(1, noNeighborDomains()), AT_, "noWindowCellsSend");
  for(MInt d = 0; d < noNeighborDomains(); d++) {
    noWindowCellsSend[d] = noWindowCells(d);
    MPI_Isend(&noWindowCellsSend[d], 1, type_traits<MInt>::mpiType(), neighborDomain(d), domainId(), mpiComm(),
              &sendRequests[d], AT_, "noWindowCellsSend[d]");
  }

  // Finish MPI communication
  MPI_Waitall(noNeighborDomains(), &recvRequests[0], MPI_STATUSES_IGNORE, AT_);
  MPI_Waitall(noNeighborDomains(), &sendRequests[0], MPI_STATUSES_IGNORE, AT_);

  // Check if received number of window cells matches the local number of halo cells
  for(MInt d = 0; d < noNeighborDomains(); d++) {
    if(noWindowCellsRecv[d] != noHaloCells(d)) {
      TERMM(1, "Solver " + to_string(m_solverId) + " d=" + to_string(domainId())
                   + ": Number of window cells from domain " + to_string(neighborDomain(d))
                   + " does not match local number of halo cells; window: " + to_string(noWindowCellsRecv[d])
                   + " ,halo: " + to_string(noHaloCells(d)));
    }
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Check #2: grid global ids of window/halo cells match
  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Start receiving window cell global ids from each neighbor domain
  fill(recvRequests.begin(), recvRequests.end(), MPI_REQUEST_NULL);
  const MInt totalNoWindowCellsRecv = accumulate(noWindowCellsRecv.begin(), noWindowCellsRecv.end(), 0);
  MIntScratchSpace windowCellsRecv(max(1, totalNoWindowCellsRecv), AT_, "windowCellsRecv");
  fill(windowCellsRecv.begin(), windowCellsRecv.end(), -1);
  for(MInt d = 0, offset = 0; d < noNeighborDomains(); d++) {
    if(noHaloCells(d) > 0) {
      MPI_Irecv(&windowCellsRecv[offset], noHaloCells(d), type_traits<MInt>::mpiType(), neighborDomain(d),
                neighborDomain(d), mpiComm(), &recvRequests[d], AT_, "windowCellsRecv[offset]");
    }
    offset += noHaloCells(d);
  }

  // Start sending window cell global ids to each neighbor domain
  fill(sendRequests.begin(), sendRequests.end(), MPI_REQUEST_NULL);
  const MInt totalNoWindowCellsSend = accumulate(noWindowCellsSend.begin(), noWindowCellsSend.end(), 0);
  MIntScratchSpace windowCellsSend(max(1, totalNoWindowCellsSend), AT_, "windowCellsSend");
  for(MInt d = 0, offset = 0; d < noNeighborDomains(); d++) {
    for(MInt c = 0; c < noWindowCellsSend[d]; c++) {
      windowCellsSend[offset + c] = raw().treeb().globalId(m_tree.solver2grid(windowCell(d, c)));
    }
    if(noWindowCells(d) > 0) {
      MPI_Isend(&windowCellsSend[offset], noWindowCells(d), type_traits<MInt>::mpiType(), neighborDomain(d), domainId(),
                mpiComm(), &sendRequests[d], AT_, "windowCellsSend[offset]");
    }
    offset += noWindowCells(d);
  }

  // Finish MPI communication
  MPI_Waitall(noNeighborDomains(), &recvRequests[0], MPI_STATUSES_IGNORE, AT_);
  MPI_Waitall(noNeighborDomains(), &sendRequests[0], MPI_STATUSES_IGNORE, AT_);

  // Check if received window cell global ids match the local halo cell global ids
  for(MInt d = 0, offset = 0; d < noNeighborDomains(); d++) {
    for(MInt c = 0; c < noHaloCells(d); c++) {
      const MInt cellId = m_tree.solver2grid(haloCell(d, c));
      const MInt globalId = raw().treeb().globalId(cellId);
      // If halo cell has periodic flag, its global id should be -1
      if(windowCellsRecv[offset + c] != globalId
         && !(raw().treeb().hasProperty(cellId, Cell::IsPeriodic) && globalId == -1)) {
        TERMM(1, "Global id of window cell " + to_string(c) + " from domain " + to_string(neighborDomain(d))
                     + " does not match local halo cell gobal id (" + to_string(windowCellsRecv[offset + c]) + " vs. "
                     + to_string(raw().treeb().globalId(m_tree.solver2grid(haloCell(d, c)))) + ")");
      }
    }
    offset += noHaloCells(d);
  }
}


/** \brief Exchange DG halo/window cell data for visualization with ParaView
 **
 ** @author: Ansgar Niemoeller
 ** @date: Yesterday
 **/
template <MInt nDim>
void Proxy<nDim>::exchangeHaloCellsForVisualizationDG(MFloat* data,
                                                      const MInt* const polyDegs,
                                                      const MInt* const dataOffsets) {
  std::vector<MInt> dataSizeSend(noNeighborDomains(), 0);
  std::vector<MInt> dataSizeRecv(noNeighborDomains(), 0);
  MInt totalSendSize = 0;
  MInt totalRecvSize = 0;

  for(MInt nghbrId = 0; nghbrId < noNeighborDomains(); nghbrId++) {
    const MInt noWindowSend = noWindowCells(nghbrId);
    const MInt noHaloRecv = noHaloCells(nghbrId);

    // Determine number of window cell DOFs
    MInt noDOFsSend = 0;
    for(MInt i = 0; i < noWindowSend; i++) {
      noDOFsSend += ipow(polyDegs[windowCell(nghbrId, i)] + 1, nDim);
    }
    dataSizeSend[nghbrId] = noDOFsSend;
    totalSendSize += noDOFsSend;

    // Determine number of halo cell DOFs
    MInt noDOFsRecv = 0;
    for(MInt i = 0; i < noHaloRecv; i++) {
      noDOFsRecv += ipow(polyDegs[haloCell(nghbrId, i)] + 1, nDim);
    }
    dataSizeRecv[nghbrId] = noDOFsRecv;
    totalRecvSize += noDOFsRecv;
  }

  // Storage for send data
  std::vector<MFloat> windowData(totalSendSize, std::numeric_limits<MFloat>::infinity());

  // Copy data to send buffer
  MInt sendBufCtr = 0;
  for(MInt nghbrId = 0; nghbrId < noNeighborDomains(); nghbrId++) {
    const MInt noWindowSend = noWindowCells(nghbrId);

    for(MInt i = 0; i < noWindowSend; i++) {
      const MInt cellId = windowCell(nghbrId, i);
      const MInt noDOFs = ipow(polyDegs[cellId] + 1, nDim);

      std::copy(&data[dataOffsets[cellId]], &data[dataOffsets[cellId] + noDOFs], &windowData[sendBufCtr]);

      sendBufCtr += noDOFs;
    }
  }

  if(sendBufCtr != totalSendSize) {
    TERMM(1, "Send buffer counter does not match total send size.");
  }

  // Receive buffer for halo data
  std::vector<MFloat> haloData(totalRecvSize, std::numeric_limits<MFloat>::infinity());
  // Communicate data with neighbor domains, data is stored ordered by neighbor domain in the buffer
  maia::mpi::exchangeData(&windowData[0], m_domainId, noNeighborDomains(), &m_neighborDomains[0], mpiComm(), 1,
                          &dataSizeSend[0], &dataSizeRecv[0], &haloData[0]);

  // Reorder data from halo cell buffer
  MInt haloDataOffset = 0;
  for(MInt nghbrId = 0; nghbrId < noNeighborDomains(); nghbrId++) {
    const MInt noHaloRecv = noHaloCells(nghbrId);

    // Loop over all halos from this neighbor domain and copy to corresponding position in data
    for(MInt i = 0; i < noHaloRecv; i++) {
      const MInt haloCellId = haloCell(nghbrId, i);
      const MInt noDOFsHalo = ipow(polyDegs[haloCellId] + 1, nDim);
      std::copy_n(&haloData[haloDataOffset], noDOFsHalo, &data[dataOffsets[haloCellId]]);
      haloDataOffset += noDOFsHalo;
    }
  }
}

/** \brief Exchange DG-SBP halo/window cell data for visualization with ParaView
 **
 ** @author: Ansgar Niemoeller
 ** @date: Yesterday
 **/
template <MInt nDim>
void Proxy<nDim>::exchangeHaloCellsForVisualizationSBP(MFloat* data,
                                                       const MInt* const noNodes1D,
                                                       const MInt* const dataOffsets) {
  std::vector<MInt> dataSizeSend(noNeighborDomains(), 0);
  std::vector<MInt> dataSizeRecv(noNeighborDomains(), 0);
  MInt totalSendSize = 0;
  MInt totalRecvSize = 0;

  for(MInt nghbrId = 0; nghbrId < noNeighborDomains(); nghbrId++) {
    const MInt noWindowSend = noWindowCells(nghbrId);
    const MInt noHaloRecv = noHaloCells(nghbrId);

    // Determine number of window cell DOFs
    MInt noDOFsSend = 0;
    for(MInt i = 0; i < noWindowSend; i++) {
      noDOFsSend += ipow(noNodes1D[windowCell(nghbrId, i)], nDim);
    }
    dataSizeSend[nghbrId] = noDOFsSend;
    totalSendSize += noDOFsSend;

    // Determine number of halo cell DOFs
    MInt noDOFsRecv = 0;
    for(MInt i = 0; i < noHaloRecv; i++) {
      noDOFsRecv += ipow(noNodes1D[haloCell(nghbrId, i)], nDim);
    }
    dataSizeRecv[nghbrId] = noDOFsRecv;
    totalRecvSize += noDOFsRecv;
  }

  // Storage for send data
  std::vector<MFloat> windowData(totalSendSize, std::numeric_limits<MFloat>::infinity());

  // Copy data to send buffer
  MInt sendBufCtr = 0;
  for(MInt nghbrId = 0; nghbrId < noNeighborDomains(); nghbrId++) {
    const MInt noWindowSend = noWindowCells(nghbrId);

    for(MInt i = 0; i < noWindowSend; i++) {
      const MInt cellId = windowCell(nghbrId, i);
      const MInt noDOFs = ipow(noNodes1D[cellId], nDim);

      std::copy(&data[dataOffsets[cellId]], &data[dataOffsets[cellId] + noDOFs], &windowData[sendBufCtr]);

      sendBufCtr += noDOFs;
    }
  }

  if(sendBufCtr != totalSendSize) {
    TERMM(1, "Send buffer counter does not match total send size.");
  }

  // Receive buffer for halo data
  std::vector<MFloat> haloData(totalRecvSize, std::numeric_limits<MFloat>::infinity());
  // Communicate data with neighbor domains, data is stored ordered by neighbor domain in the buffer
  maia::mpi::exchangeData(&windowData[0], m_domainId, noNeighborDomains(), &m_neighborDomains[0], mpiComm(), 1,
                          &dataSizeSend[0], &dataSizeRecv[0], &haloData[0]);

  // Reorder data from halo cell buffer
  MInt haloDataOffset = 0;
  for(MInt nghbrId = 0; nghbrId < noNeighborDomains(); nghbrId++) {
    const MInt noHaloRecv = noHaloCells(nghbrId);

    // Loop over all halos from this neighbor domain and copy to corresponding position in data
    for(MInt i = 0; i < noHaloRecv; i++) {
      const MInt haloCellId = haloCell(nghbrId, i);
      const MInt noDOFsHalo = ipow(noNodes1D[haloCellId], nDim);
      std::copy_n(&haloData[haloDataOffset], noDOFsHalo, &data[dataOffsets[haloCellId]]);
      haloDataOffset += noDOFsHalo;
    }
  }
}

/// Determine on which side of the azimuthal periodic boundary coords are located
/// \author Thomas Hoesgen
/// \date 2022-03-20
template <MInt nDim>
MInt Proxy<nDim>::determineAzimuthalBoundarySide(const MFloat* coords) {
  TRACE();

  MFloat side = 0;
  MFloat coordsCyl[3];

  raw().cartesianToCylindric(coords, coordsCyl);

  side = raw().m_azimuthalBbox.azimuthalSide(/*coords,*/ coordsCyl[1]);

  return side;
}

/** \brief update proxy cutOff information
 **
 ** @author: Tim Wegmann
 ** @date: the day before yesterday
 **/
template <MInt nDim>
void Proxy<nDim>::updateCutOff() {
  TRACE();

  m_tree.m_isCutOffCell.assign(noCells(), -1);

  if(hasCutOff()) {
    /*! \page propertyPage1
    \section cutOffDirections
    <code>cutOffDirections </code>\n
    default = <code>none</code>\n \n
    Create cut off boundary at the respective border of the domain (0: create cut off bndry in -x direction, and so on).
    Only works together with property cutOffBndryIds <ul> <li> 0, 1, 2, 3, 4 (2D), 5 (3D) </li>
    </ul>
    Keywords: <i> FINITE_VOLUME, CUTOFF </i>
  */
    const MInt noCutOffDirections = Context::propertyLength("cutOffDirections", solverId());

    ScratchSpace<MInt> isCutOff(m_tree.size(), AT_, "isCutOff");
    isCutOff.fill(-1);

    // mark internal cutOff cells
    for(MInt i = 0; i < noCutOffDirections; i++) {
      const MInt dir = Context::getSolverProperty<MInt>("cutOffDirections", solverId(), AT_, i);

      for(MInt cellId = 0; cellId < noInternalCells(); cellId++) {
        if(tree().hasNeighbor(cellId, dir) != 0) continue;
        if(!tree().isLeafCell(cellId)) continue;
        if(tree().parent(cellId) > -1) {
          if(tree().hasNeighbor(tree().parent(cellId), dir) != 0) continue;
        }
        if(m_tree.m_isCutOffCell[cellId] < 0) {
          m_tree.m_isCutOffCell[cellId] = dir;
          isCutOff[cellId] = dir;
        } else {
          // NOTE: if instead of the direction prime numbers are used
          //      the addition of the prime number leads to the recalculation of the used directions!
          //      for cells which are part of multiple cutOffs!
          m_tree.m_isCutOffCell[cellId] += dir;
          isCutOff[cellId] += dir;
        }
      }
    }

    // exchange property to mark halo cutOff cells
    if(m_neighborDomains.size() > 0) {
      mpi::exchangeData(m_neighborDomains, m_haloCells, m_windowCells, mpiComm(), &isCutOff[0], 1);
    }

    for(MInt cellId = noInternalCells(); cellId < m_tree.size(); cellId++) {
      m_tree.m_isCutOffCell[cellId] = isCutOff[cellId];
    }
  }
}


template <MInt nDim>
void Proxy<nDim>::findEqualLevelNeighborsParDiagonal(MBool idsAreGlobal) {
  static constexpr MInt maxNoNeighbors = (nDim == 2) ? 8 : 26;

  m_log << "detecting neighbors for " << noCells() << " cells ... ";

  // Allocate neighbor list and neighbor element
  m_log << "allocating new array (" << noCells() << " entries)...";
  mDeallocate(m_neighborList);
  mAlloc(m_neighborList, noCells(), maxNoNeighbors, "m_neighborList", -1, AT_);
  m_log << " allocated... ";

  static constexpr MInt twoPowers[6] = {1, 2, 4, 8, 16, 32};

  IF_CONSTEXPR(nDim == 2) {
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(MInt id = 0; id < noCells(); id++) { // m_noRegularCells;

      // -----------------------------
      // Find all equal level neighbors
      for(MInt i = 0; i < raw().m_counter2D; i++) {
        // First step
        MInt current = id;

        if(tree().hasNeighbor(current, raw().m_paths2D[i][0]) == 1) {
          // take the global Id
          MInt currentGlobal = tree().neighbor(current, raw().m_paths2D[i][0]);

          MInt nCode = raw().m_neighborCode[twoPowers[raw().m_paths2D[i][0]]];
          // save also global neighbor Ids laying in other domains (here: no exception for outer laying domain cells)
          m_neighborList[id][nCode] = currentGlobal;

          if(currentGlobal == -1) {
            continue;
          }

          // continue if cell is laying in other domain
          if(idsAreGlobal && raw().m_globalToLocalId.find(currentGlobal) == raw().m_globalToLocalId.end()) {
            continue;
          }

          // take the local Id
          current = idsAreGlobal ? raw().m_globalToLocalId.at(tree().neighbor(current, raw().m_paths2D[i][0]))
                                 : tree().neighbor(current, raw().m_paths2D[i][0]);

          // Second step
          if(tree().hasNeighbor(current, raw().m_paths2D[i][1]) == 1) {
            // save also global neighbor Ids laying in other domains (here: no exception for outer laying domain cells)
            nCode = raw().m_neighborCode[twoPowers[raw().m_paths2D[i][0]] + twoPowers[raw().m_paths2D[i][1]]];
            m_neighborList[id][nCode] = tree().neighbor(current, raw().m_paths2D[i][1]);
          }
        }
      } // End loop over neighbor directions

      // For LB: Find also diagonal neighbors that are not connected via cartesian neighbors (2D)!
      // Some mapping that is needed:
      // E.g. For child on position 0, child on position 3 in parent is the wanted diagonal neighbor
      if(raw().m_lbGridChecks) {
        MInt posPosMap[4] = {3, 2, 1, 0};

        // E.g. For child 0, the neighbor that can be affected by this exceptional case is neighbor 6
        MInt posDiagNghbrMap[4] = {6, 5, 7, 4};

        // E.g. For child 0, the paths on paren level are 0,2 or 2,0
        MInt parentPaths[4][2][2] = {{{0, 2}, {2, 0}}, {{1, 2}, {2, 1}}, {{0, 3}, {3, 0}}, {{1, 3}, {3, 1}}};
        MInt numberOfPossiblePaths = 2;

        // This case can only happen if cells are not on maxUniformRefinementLevel
        if(tree().level(id) > maxUniformRefinementLevel()) {
          MInt parent = tree().parent(id);

          // Get postion of child in parent
          MInt position = 0;
          for(MInt dim = 0; dim < 2; dim++) {
            if(tree().coordinate(id, dim) < tree().coordinate(parent, dim)) {
              position += 0;
            } else {
              position += IPOW2(dim);
            }
          }

          // Check if diagonal neighbor has not been already found via cartesian neighbor link
          for(MInt j = 0; j < numberOfPossiblePaths; j++) {
            if(m_neighborList[id][posDiagNghbrMap[position]] == -1) {
              // Perform first step on parent level
              if(tree().hasNeighbor(parent, parentPaths[position][j][0]) == 1) {
                MInt firstParentNghbrGlobal = tree().neighbor(parent, parentPaths[position][j][0]);
                if(firstParentNghbrGlobal == -1) continue;
                if(idsAreGlobal
                   && raw().m_globalToLocalId.find(firstParentNghbrGlobal) == raw().m_globalToLocalId.end())
                  continue;
                MInt firstParentNghbr =
                    idsAreGlobal ? raw().m_globalToLocalId.at(tree().neighbor(parent, parentPaths[position][j][0]))
                                 : tree().neighbor(parent, parentPaths[position][j][0]);

                // Perform second step on parent level
                if(tree().hasNeighbor(firstParentNghbr, parentPaths[position][j][1]) == 1) {
                  MInt secondParentNghbrGlobal = tree().neighbor(firstParentNghbr, parentPaths[position][j][1]);
                  if(secondParentNghbrGlobal == -1) continue;
                  if(idsAreGlobal
                     && raw().m_globalToLocalId.find(secondParentNghbrGlobal) == raw().m_globalToLocalId.end())
                    continue;
                  MInt secondParentNghbr =
                      idsAreGlobal
                          ? raw().m_globalToLocalId.at(tree().neighbor(firstParentNghbr, parentPaths[position][j][1]))
                          : tree().neighbor(firstParentNghbr, parentPaths[position][j][1]);
                  if(tree().hasChildren(secondParentNghbr)) {
                    if(tree().child(secondParentNghbr, posPosMap[position]) != -1) {
                      MInt nCode = posDiagNghbrMap[position];
                      m_neighborList[id][nCode] = tree().child(secondParentNghbr, posPosMap[position]);
                      break;
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  else {
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(MInt id = 0; id < noCells(); id++) { // m_noRegularCells

      // -----------------------------
      // Find all equal level neighbors
      for(MInt i = 0; i < raw().m_counter3D; i++) {
        // First step
        MInt current = id;
        if(tree().hasNeighbor(current, raw().m_paths3D[i][0]) == 1) {
          // take the global Id
          MInt currentGlobal = tree().neighbor(current, raw().m_paths3D[i][0]);

          MInt nCode = raw().m_neighborCode[twoPowers[raw().m_paths3D[i][0]]];
          // save also global neighbor Ids laying in other domains (here: no exception for outer laying domain cells)
          m_neighborList[id][nCode] = currentGlobal;

          if(currentGlobal == -1) {
            continue;
          }

          // continue if cell is laying in other domain
          if(idsAreGlobal && raw().m_globalToLocalId.find(currentGlobal) == raw().m_globalToLocalId.end()) {
            continue;
          }

          // take the local Id
          current = idsAreGlobal ? raw().m_globalToLocalId.at(tree().neighbor(current, raw().m_paths3D[i][0]))
                                 : tree().neighbor(current, raw().m_paths3D[i][0]);

          // Second step
          if(tree().hasNeighbor(current, raw().m_paths3D[i][1]) == 1) {
            // take the global Id
            currentGlobal = tree().neighbor(current, raw().m_paths3D[i][1]);

            nCode = raw().m_neighborCode[twoPowers[raw().m_paths3D[i][0]] + twoPowers[raw().m_paths3D[i][1]]];
            // save also global neighbor Ids laying in other domains (here: no exception for outer laying domain cells)
            m_neighborList[id][nCode] = currentGlobal;

            if(currentGlobal == -1) {
              continue;
            }

            // continue if cell is laying in other domain
            if(idsAreGlobal && raw().m_globalToLocalId.find(currentGlobal) == raw().m_globalToLocalId.end()) {
              continue;
            }

            // take the local Id
            current = idsAreGlobal ? raw().m_globalToLocalId.at(tree().neighbor(current, raw().m_paths3D[i][1]))
                                   : tree().neighbor(current, raw().m_paths3D[i][1]);

            // Third step
            if(tree().hasNeighbor(current, raw().m_paths3D[i][2]) == 1) {
              // save also global neighbor Ids laying in other domains (here: no exception for outer laying domain
              // cells)
              nCode = raw().m_neighborCode[twoPowers[raw().m_paths3D[i][0]] + twoPowers[raw().m_paths3D[i][1]]
                                           + twoPowers[raw().m_paths3D[i][2]]];
              m_neighborList[id][nCode] = tree().neighbor(current, raw().m_paths3D[i][2]);
            }
          }
        }
      }

      // For LB: Find also diagonal neighbors that are not connected via cartesian neighbors (3D).
      // Some mapping:
      // E.g. child at pos. 0 has neighbors 6,10,14,18 that can be affected by the exception
      if(raw().m_lbGridChecks) {
        MInt posDiagNghbrMap[8][4] = {{6, 10, 14, 18}, {8, 12, 14, 22}, {7, 10, 16, 20}, {9, 12, 16, 24},
                                      {6, 11, 15, 19}, {8, 13, 15, 23}, {7, 11, 17, 21}, {9, 13, 17, 25}};

        // E.g. the neighbors of child 0 in the parent neighbors 6,10,14,18 are the children 3,5,6,7
        MInt posPosMap[8][4] = {{3, 5, 6, 7}, {2, 4, 7, 6}, {1, 7, 4, 5}, {0, 6, 5, 4},
                                {7, 1, 2, 3}, {6, 0, 3, 2}, {5, 3, 0, 1}, {4, 2, 1, 0}};

        // Index is diag direction id minus 6.
        // E.g. path to diag neighbor 6 can be 2,0 or 0,2
        MInt diagNghbrPathMap[20][4][3] = {
            {{2, 0, -1}, {0, 2, -1}, {-1, -1, -1}, {-1, -1, -1}}, {{3, 0, -1}, {0, 3, -1}, {-1, -1, -1}, {-1, -1, -1}},
            {{2, 1, -1}, {1, 2, -1}, {-1, -1, -1}, {-1, -1, -1}}, {{3, 1, -1}, {1, 3, -1}, {-1, -1, -1}, {-1, -1, -1}},
            {{4, 0, -1}, {0, 4, -1}, {-1, -1, -1}, {-1, -1, -1}}, {{5, 0, -1}, {0, 5, -1}, {-1, -1, -1}, {-1, -1, -1}},
            {{4, 1, -1}, {1, 4, -1}, {-1, -1, -1}, {-1, -1, -1}}, {{5, 1, -1}, {1, 5, -1}, {-1, -1, -1}, {-1, -1, -1}},
            {{4, 2, -1}, {2, 4, -1}, {-1, -1, -1}, {-1, -1, -1}}, {{5, 2, -1}, {2, 5, -1}, {-1, -1, -1}, {-1, -1, -1}},
            {{4, 3, -1}, {3, 4, -1}, {-1, -1, -1}, {-1, -1, -1}}, {{5, 3, -1}, {3, 5, -1}, {-1, -1, -1}, {-1, -1, -1}},
            {{0, 2, 4}, {2, 0, 4}, {4, 2, 0}, {4, 0, 2}},         {{0, 2, 5}, {2, 0, 5}, {5, 0, 2}, {5, 2, 0}},
            {{0, 3, 4}, {3, 0, 4}, {4, 0, 3}, {4, 3, 0}},         {{0, 3, 5}, {3, 0, 5}, {5, 0, 3}, {5, 3, 0}},
            {{1, 2, 4}, {2, 1, 4}, {4, 1, 2}, {4, 2, 1}},         {{1, 2, 5}, {2, 1, 5}, {5, 1, 2}, {5, 2, 1}},
            {{1, 3, 4}, {3, 1, 4}, {4, 1, 3}, {4, 3, 1}},         {{1, 3, 5}, {3, 1, 5}, {5, 1, 3}, {5, 3, 1}}};

        // This case can only happen if cells are not on maxUniformRefinementlevel
        if(tree().level(id) > maxUniformRefinementLevel()) {
          MInt parent = tree().parent(id);

          // Get position of child in parent
          MInt position = 0;
          for(MInt dim = 0; dim < nDim; dim++) {
            if(tree().coordinate(id, dim) < tree().coordinate(parent, dim)) {
              position += 0;
            } else {
              position += IPOW2(dim);
            }
          }

          // Loop over possible affected neighbors
          for(MInt i = 0; i < 4; i++) {
            MInt neighbor = posDiagNghbrMap[position][i];

            // Check if diagonal neighbor has not been already found via Cartesian neighbor link
            if(m_neighborList[id][neighbor] == -1) { // not equal for testing :) --> Overwrites above. Thus, must be the
                                                     // same as without this extension!

              // Loop over paths
              for(MInt j = 0; j < 4; j++) {
                // Continue if there is no path left
                if(diagNghbrPathMap[neighbor - 6][j][0] == -1) continue;

                // Perform first step on parent level
                if(tree().hasNeighbor(parent, diagNghbrPathMap[neighbor - 6][j][0]) == 1) {
                  MInt firstParentNghbrGlobal = tree().neighbor(parent, diagNghbrPathMap[neighbor - 6][j][0]);
                  if(firstParentNghbrGlobal == -1) continue;
                  if(idsAreGlobal
                     && raw().m_globalToLocalId.find(firstParentNghbrGlobal) == raw().m_globalToLocalId.end())
                    continue;
                  MInt firstParentNghbr =
                      idsAreGlobal
                          ? raw().m_globalToLocalId.at(tree().neighbor(parent, diagNghbrPathMap[neighbor - 6][j][0]))
                          : tree().neighbor(parent, diagNghbrPathMap[neighbor - 6][j][0]);

                  // Perform second step on parent level
                  if(tree().hasNeighbor(firstParentNghbr, diagNghbrPathMap[neighbor - 6][j][1]) == 1) {
                    MInt secondParentNghbrGlobal =
                        tree().neighbor(firstParentNghbr, diagNghbrPathMap[neighbor - 6][j][1]);
                    if(secondParentNghbrGlobal == -1) continue;
                    if(idsAreGlobal
                       && raw().m_globalToLocalId.find(secondParentNghbrGlobal) == raw().m_globalToLocalId.end())
                      continue;
                    MInt secondParentNghbr =
                        idsAreGlobal ? raw().m_globalToLocalId.at(
                            tree().neighbor(firstParentNghbr, diagNghbrPathMap[neighbor - 6][j][1]))
                                     : tree().neighbor(firstParentNghbr, diagNghbrPathMap[neighbor - 6][j][1]);

                    // Check if a third step is necessary
                    if(diagNghbrPathMap[neighbor - 6][j][2] != -1) {
                      if(tree().hasNeighbor(secondParentNghbr, diagNghbrPathMap[neighbor - 6][j][2]) == 1) {
                        MInt thirdParentNghbrGlobal =
                            tree().neighbor(secondParentNghbr, diagNghbrPathMap[neighbor - 6][j][2]);
                        if(thirdParentNghbrGlobal == -1) continue;
                        if(idsAreGlobal
                           && raw().m_globalToLocalId.find(thirdParentNghbrGlobal) == raw().m_globalToLocalId.end())
                          continue;
                        MInt thirdParentNghbr =
                            idsAreGlobal ? raw().m_globalToLocalId.at(
                                tree().neighbor(secondParentNghbr, diagNghbrPathMap[neighbor - 6][j][2]))
                                         : tree().neighbor(secondParentNghbr, diagNghbrPathMap[neighbor - 6][j][2]);
                        if(tree().hasChildren(thirdParentNghbr)) {
                          if(tree().child(thirdParentNghbr, posPosMap[position][i]) != -1) {
                            MInt nCode = neighbor;
                            m_neighborList[id][nCode] = tree().child(thirdParentNghbr, posPosMap[position][i]);
                            break;
                          }
                        }
                      }
                    } else {
                      if(tree().hasChildren(secondParentNghbr)) {
                        if(tree().child(secondParentNghbr, posPosMap[position][i]) != -1) {
                          MInt nCode = neighbor;
                          m_neighborList[id][nCode] = tree().child(secondParentNghbr, posPosMap[position][i]);
                          break;
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}


/**
 * identify and store all direct neighbors
 * function sweeps over all cartesian leaf cells
 * and determines all of their direct leaf neighbors
 * neighbor data is stored into two scratch arrays that are allocated by
 * the calling method and transferred
 *
 * @author Daniel Hartmann, multiSolver update Tim Wegmann
 */
template <MInt nDim>
void Proxy<nDim>::findCartesianNghbIds() {
  TRACE();

  if(m_storeNghbrIds != nullptr) {
    mDeallocate(m_storeNghbrIds);
  }
  mAlloc(m_storeNghbrIds, m_noDirs * IPOW2(nDim) * noCells(), "m_storeNghbrIds", -1, AT_);

  if(m_identNghbrIds != nullptr) {
    mDeallocate(m_identNghbrIds);
  }
  mAlloc(m_identNghbrIds, m_noDirs * noCells() + 1, "m_identNghbrIds", -1, AT_);

  MInt pos = 0;

  // --------- loop over all cartesian cells ----------
  for(MInt cellId = 0; cellId < noCells(); cellId++) {
    for(MInt dirId = 0; dirId < m_noDirs; dirId++) {
      m_identNghbrIds[cellId * m_noDirs + dirId] = pos;
    }
    if(!tree().isLeafCell(cellId)) continue;
    for(MInt dirId = 0; dirId < m_noDirs; dirId++) {
      const MInt spaceId = dirId / 2;
      const MInt sideId = dirId % 2;
      const MInt nghbrSideId = (sideId + 1) % 2;

      // set the location in m_storeNghbrIds, where neighbors of cellId are
      // stored
      m_identNghbrIds[cellId * m_noDirs + dirId] = pos;

      // ----------------------------------------------------------------
      // --- identify neighbor in considered direction ------------------
      // *-> this can be a neighbor of the cell itself or a neighbor of *
      // *   its parent, grandparent, and so forth                      *
      // ----------------------------------------------------------------

      // check if cell has a neighbor in the investigated
      // direction
      if(tree().hasNeighbor(cellId, dirId) > 0) {
        // -> cell has a neighbor
        // this neighbor:
        MInt nghbrId = tree().neighbor(cellId, dirId);
        const MInt nghbrLvl = tree().level(nghbrId);

        // --------------------------------------------------------------
        // --- investigate neighbor -------------------------------------
        // *-> check if neighbor has children, grandchildren, and so    *
        // *   forth, which are direct neighbors of the considered cell *
        // --------------------------------------------------------------

        if(tree().isLeafCell(nghbrId)) {
          // neighbor is not a parent and can be stored as neighbor of the cell
          m_storeNghbrIds[pos] = nghbrId;
          pos++;
        } else {
          // neighbor is a parent
          // investigate all children and grandchildren
          // temporary array that holds all parent cells to be
          // investigated
          MInt tempNghbrs[16];
          MInt tempNghbrs1[16];
          tempNghbrs[0] = nghbrId;
          MBool nghbrsFound = false;

          // ----------------------------------------------------------
          // --- investigate all neighbors that are parents -----------
          // ----------------------------------------------------------

          MInt maxIndex = 1;
          while(nghbrsFound == false) {
            // initialize counter for the number of investigated
            // neighboring children that are parents themselves and
            // thus added to a temporary array
            MInt k = 0;

            // ------ loop over all neighbor parents ------------------
            for(MInt index = 0; index < maxIndex; index++) {
              // investigated neighbor parent cell
              ASSERT(index > -1 && index < 16, "");
              nghbrId = tempNghbrs[index];

              // check the location of the children within the parent cell
              // by means of comparison of the grid cell cogs of the
              // parent cell and the children in considered space direction

              // in case none of the children is identified
              // as a direct neighbor, do nothing, since then
              // the considered cell side is a non-fluid side

              // --------- loop over neighbor children ----------------
              for(MInt childId = 0; childId < IPOW2(nDim); childId++) {
                const MInt nghbrChildId = tree().child(nghbrId, childId);
                if(nghbrChildId < 0) continue;
                if((tree().coordinate(nghbrChildId, spaceId) - tree().coordinate(nghbrId, spaceId))
                       * ((MFloat)nghbrSideId - 0.5)
                   > 0) {
                  if(tree().isLeafCell(nghbrChildId)) {
                    // child cell is not a parent and stored as a neighbor
                    m_storeNghbrIds[pos] = nghbrChildId;
                    pos++;
                  } else {
                    // add child cell to temporary array for
                    // further investigation
                    // ASSERT(k > -1 && k < 16,
                    //       "Too many local child neighbors. Something might be wrong with "
                    //       "your mesh!");
                    if(k > 14) {
                      cerr << "Large neighbor count " << k << " " << nghbrLvl << " " << tree().level(nghbrChildId)
                           << endl;
                    }
                    tempNghbrs1[k] = nghbrChildId;
                    k++;
                  }
                }
              }
            }

            maxIndex = k;
            if(k == 0) {
              nghbrsFound = true;
            } else {
              // transfer the cells collected in tempNghbrs1 to
              // tempNghbrs
              for(MInt l = 0; l < k; l++) {
                ASSERT(l > -1 && l < 16, "");
                tempNghbrs[l] = tempNghbrs1[l];
              }
            }
          }
        }
      } else {
        // -> cell does not have a neighbor of equivalent level
        MInt currentCellId = cellId;
        MBool nghbrExist = false;
        // loop until a parent is found who has a neighbor in considered direction
        // break when the parent ID becomes -1 (boundary/cutOff in that direction)
        while(nghbrExist == false) {
          const MInt parentId = tree().parent(currentCellId);
          // if parent does not exist, break
          if(parentId < 0 || parentId >= noCells()) break;


          // break if currentCell does not lie at the parent cell face
          // whose normal points to considered direction
          if((tree().coordinate(currentCellId, spaceId) - tree().coordinate(parentId, spaceId)) * ((float)sideId - 0.5)
             < 0) {
            break;
          }

          // identify parentId of current cell and store it as current cell
          currentCellId = parentId;

          // check if current cell has a neighbor in regarded direction
          if(tree().hasNeighbor(currentCellId, dirId) > 0) {
            nghbrExist = true;
            const MInt nghbrId = tree().neighbor(currentCellId, dirId);
            // store this neighbor if it is regarded - multilevel
            if(tree().isLeafCell(nghbrId)) {
              m_storeNghbrIds[pos] = nghbrId;
              pos++;
            }
          }
        }
      }
    }
  }
  m_identNghbrIds[noCells() * m_noDirs] = pos;
}

/*! \brief Obtain list of direct neighbors of given cell.
 *         Requires m_identNghbrIds, as such findCartesianNghbIds() must have been called
 *         before!
 * \author Sven Berger, Tim Wegmann
 * \date May 2016
 * \param[in] cellId id of the cell for which all neighbors should be obtained
 * \param[out] nghbrList list of neighboring cells
 *
 */
template <MInt nDim>
void Proxy<nDim>::findDirectNghbrs(const MInt cellId, vector<MInt>& nghbrList) {
  ASSERT(cellId >= 0, "ERROR: Invalid cellId!");
  ASSERT(m_storeNghbrIds != nullptr, "");
  ASSERT(m_identNghbrIds != nullptr, "");

  nghbrList.push_back(cellId);

  // this prevents neighbors to be identified multiple times
  auto uniqueAdd = [&](MInt newElement) {
    auto it = find(nghbrList.begin(), nghbrList.end(), newElement);
    if(it == nghbrList.end() && newElement >= 0) {
      nghbrList.push_back(newElement);
      return true;
    }
    return false;
  };

  ///
  /// Gather all neighbouring cells and Diagonal neighbor(s) in the xy-plane (in total 10 additional cells (3D)):
  ///
  static constexpr MInt diagDirs[4]{3, 2, 0, 1};
  for(MInt dir = 0; dir < 2 * nDim; dir++) {
    // iterate overall direct neighbors
    for(MInt j = m_identNghbrIds[2 * cellId * nDim + dir]; j < m_identNghbrIds[2 * cellId * nDim + dir + 1]; j++) {
      const MInt nghbrId = m_storeNghbrIds[j];
      uniqueAdd(nghbrId);
      if(dir < 4) {
        // add diagonal neighbors in x-y plane
        for(MInt t = m_identNghbrIds[2 * nghbrId * nDim + diagDirs[dir]];
            t < m_identNghbrIds[2 * nghbrId * nDim + diagDirs[dir] + 1];
            t++) {
          const MInt diagNghbr = m_storeNghbrIds[t];
          uniqueAdd(diagNghbr);
        }
      }
    }
  }

  ///
  /// In the 3D case look for additional diagonal and tridiagonal neighbors in the z-directions (in total 16 additional
  /// cells (3D)):
  ///
  IF_CONSTEXPR(nDim == 3) {
    for(MInt dirZ = 4; dirZ < 6; dirZ++) {
      for(MInt j = m_identNghbrIds[2 * cellId * nDim + dirZ]; j < m_identNghbrIds[2 * cellId * nDim + dirZ + 1]; j++) {
        const MInt nghbrIdZ = m_storeNghbrIds[j];
        for(MInt dir = 0; dir < 4; dir++) {
          for(MInt t = m_identNghbrIds[2 * nghbrIdZ * nDim + dir]; t < m_identNghbrIds[2 * nghbrIdZ * nDim + dir + 1];
              t++) {
            const MInt nghbrId = m_storeNghbrIds[t];
            uniqueAdd(nghbrId);
            // tridiagonal neighbors
            for(MInt v = m_identNghbrIds[2 * nghbrId * nDim + diagDirs[dir]];
                v < m_identNghbrIds[2 * nghbrId * nDim + diagDirs[dir] + 1];
                v++) {
              const MInt triDiagNghbrId = m_storeNghbrIds[v];
              uniqueAdd(triDiagNghbrId);
            }
          }
        }
      }
    }
  }
}

/*! \brief Obtain list of neighbors for the given extend,
 *         using m_storeNghbrIds and m_identNghbrIds
 * \author Sven Berger, Tim Wegmann
 * \date May 2016
 * \param[in] cellId id of the cell for which all neighbors should be obtained
 * \param[in] layer maximum extend of neighborhood in number of cells
 * \param[out] nghbrList list of neighboring cells
 *
 */
template <MInt nDim>
void Proxy<nDim>::findNeighborHood(const MInt cellId, const MInt layer, vector<MInt>& nghbrList) {
  ASSERT(layer > 0, "ERROR: Invalid layer!");
  ASSERT(cellId >= 0, "ERROR: Invalid cellId!");
  ASSERT(m_storeNghbrIds != nullptr, "");
  ASSERT(m_identNghbrIds != nullptr, "");

  // this prevents neighbors to be identified multiple times
  auto uniqueAdd = [&](MInt newElement) {
    auto it = find(nghbrList.begin(), nghbrList.end(), newElement);
    if(it == nghbrList.end() && newElement >= 0) {
      nghbrList.push_back(newElement);
      return true;
    }
    return false;
  };

  vector<MInt> tempNghbrCollection{};

  ///
  /// 1. direct neighbors of cell
  ///
  findDirectNghbrs(cellId, nghbrList);

  ///
  /// 2. increase neighborhood to given size
  ///
  for(MInt currEx = 1; currEx < layer; currEx++) {
    // for each newly found neighbor find all neighbors
    for(auto& nghbr : nghbrList) {
      findDirectNghbrs(nghbr, tempNghbrCollection);
    }
    nghbrList.clear();
    // add all new unique neighbor ids to vector
    for(auto& nghbr : tempNghbrCollection) {
      if(uniqueAdd(nghbr) && currEx + 1 < layer) {
        // optimisation
        nghbrList.push_back(nghbr);
      }
    }
  }
}


template <MInt nDim>
MBool Proxy<nDim>::checkOutsideGeometry(MInt gridId) {
  TRACE();
  static constexpr MInt cornerIndices[8][3] = {{-1, -1, -1}, {1, -1, -1}, {-1, 1, -1}, {1, 1, -1},
                                               {-1, -1, 1},  {1, -1, 1},  {-1, 1, 1},  {1, 1, 1}};
  MFloat corner[3] = {0, 0, 0};
  MFloat cellHalfLength = raw().halfCellLength(gridId);

  MInt cnt = 0;
#ifndef PMODEC
  for(MInt i = 0; i < IPOW2(nDim); i++) {
    for(MInt dim = 0; dim < nDim; dim++) {
      corner[dim] = raw().a_coordinate(gridId, dim) + cornerIndices[i][dim] * cellHalfLength;
    }
    if(nDim == 2) {
      if(!m_geometry->pointIsInside(corner)) {
        cnt++; // pointIsInside == true if Point is outside fluid domain
        break;
      }
    } else {
      if(!m_geometry->pointIsInside2(corner)) {
        cnt++; // pointIsInside == true if Point is outside fluid domain
        break;
      }
    }
  }
#endif
  MInt outside = false;
  if(cnt == 0) {
    outside = true;
  }
  return outside;
}

/**
 * \brief Retrieves all direct and diagonal neighboring cells of the given cell on the child level if available
 * \author Lennart Schneiders
 */
template <MInt nDim>
MInt Proxy<nDim>::getAdjacentGridCells(MInt cellId, MInt noLayers, MIntScratchSpace& adjacentCells, MInt level,
                                       MInt diagonalNeighbors) {
  MInt noDirs = 2 * nDim;

  std::array<std::array<MInt, 4>, nDim> tmpNghbrs;
  set<MInt> nghbrs;

  nghbrs.insert(cellId);
  for(MInt layer = 1; layer <= noLayers; layer++) {
    set<MInt> nextLayer;
    for(auto& nghbrId : nghbrs) {
      for(MInt dir0 = 0; dir0 < noDirs; dir0++) {
        const MInt cnt0 = getNghbrCells(nghbrId, level, &tmpNghbrs[0][0], dir0);
        for(MInt c0 = 0; c0 < cnt0; c0++) {
          nextLayer.insert(tmpNghbrs[0][c0]);

          if(diagonalNeighbors > 0) {
            for(MInt dir1 = 0; dir1 < noDirs; dir1++) {
              if((dir1 / 2) == (dir0 / 2)) continue;
              const MInt cnt1 = getNghbrCells(tmpNghbrs[0][c0], level, &tmpNghbrs[1][0], dir1, dir0);
              for(MInt c1 = 0; c1 < cnt1; c1++) {
                nextLayer.insert(tmpNghbrs[1][c1]);

                if(nDim == 3 && diagonalNeighbors == 2) {
                  for(MInt dir2 = 0; dir2 < noDirs; dir2++) {
                    if(((dir2 / 2) == (dir0 / 2)) || ((dir2 / 2) == (dir1 / 2))) continue;
                    const MInt cnt2 = getNghbrCells(tmpNghbrs[1][c1], level, &tmpNghbrs[2][0], dir2, dir1, dir0);
                    for(MInt c2 = 0; c2 < cnt2; c2++) {
                      nextLayer.insert(tmpNghbrs[2][c2]);
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    for(MInt it : nextLayer) {
      nghbrs.insert(it);
    }
  }
  nghbrs.erase(cellId);

  MInt cnt = 0;
  for(MInt nghbr : nghbrs) {
    ASSERT(tree().level(nghbr) <= level + 1, "");
    ASSERT(cnt < (signed)adjacentCells.size(), to_string(nghbrs.size()) + " " + to_string(adjacentCells.size()));
    adjacentCells[cnt] = nghbr;
    cnt++;
  }
  return cnt;
}


template <MInt nDim>
MInt Proxy<nDim>::getNghbrCells(MInt cellId, MInt level, MInt* nghbrs, MInt dir0, MInt dir1, MInt dir2) {
  MInt count = 0;
  MInt nghbrId0 = -1;
  if(tree().hasNeighbor(cellId, dir0) > 0) {
    nghbrId0 = tree().neighbor(cellId, dir0);
  } else if(tree().parent(cellId) > -1) {
    if(tree().hasNeighbor(tree().parent(cellId), dir0) > 0) {
      nghbrId0 = tree().neighbor(tree().parent(cellId), dir0);
    }
  }
  if(nghbrId0 < 0) return 0;

  if(tree().noChildren(nghbrId0) > 0 && tree().level(nghbrId0) <= level) {
    for(MInt child = 0; child < m_maxNoChilds; child++) {
      if(!childCode[dir0][child]) continue;
      if(dir1 > -1 && !childCode[dir1][child]) continue;
      if(dir2 > -1 && !childCode[dir2][child]) continue;
      if(tree().child(nghbrId0, child) > -1) {
        nghbrs[count] = tree().child(nghbrId0, child);
        count++;
      }
    }
  } else {
    nghbrs[count] = nghbrId0;
    count++;
  }

  return count;
}

/// Check azimuthal neighbor domains for consistency.
///
/// \author Thomas Hoesgen
/// \date 2020-11-02
///
/// This method checks if the following information matches on corresponding domains:
/// - azimuthal neighbor domains are sorted ascendingly
/// - domains are azimuthal neighbors
/// Note: m_azimuthalNeighborDomains must already contain solver-local domain ids!
template <MInt nDim>
void Proxy<nDim>::checkNeighborConsistencyAzimuthal() const {
  TRACE();

  // Check if neighbor domains are sorted ascendingly
  if(!std::is_sorted(m_azimuthalNeighborDomains.begin(), m_azimuthalNeighborDomains.end())) {
    TERMM(1, "Neighbor domains are not sorted in ascending order.");
  }

  // Determine all neighbor domains of current domain
  MCharScratchSpace isNeighborDomain(noDomains(), AT_, "isNeighborDomain");
  fill(isNeighborDomain.begin(), isNeighborDomain.end(), 0);
  for(MInt d = 0; d < noAzimuthalNeighborDomains(); d++) {
    isNeighborDomain[azimuthalNeighborDomain(d)] = 1;
  }

  // Exchange with all domains in solver-local communicator
  MPI_Alltoall(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, &isNeighborDomain[0], 1, type_traits<MChar>::mpiType(), mpiComm(),
               AT_, "MPI_IN_PLACE", "isNeighborDomain[0]");

  // Check if all domains that the current domain considers as neighbors also reciprocate
  for(MInt d = 0; d < noAzimuthalNeighborDomains(); d++) {
    if(!isNeighborDomain[azimuthalNeighborDomain(d)]) {
      TERMM(1, "Domain " + to_string(domainId()) + " has domain " + to_string(neighborDomain(d))
                   + " as a neighbor but not the other way around.");
    }
  }
}


/// Check azimuthal window/halo cells for consistency.
///
/// \author Thomas Hoesgen
/// \date 2020-11-01
///
/// This method checks if the following information matches on corresponding domains:
/// - for each pair of azimuthal neighbor domains: number of azimuthal window/halo cells
/// - for each pair of azimuthal neighbor domains: grid global ids of azimuthal window/halo cells
template <MInt nDim>
void Proxy<nDim>::checkWindowHaloConsistencyAzimuthal() const {
  TRACE();

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Check #1: number of azimuthal window/halo cells match
  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Start receiving number of azimuthal window cells from each azimuthal neighbor domain
  ScratchSpace<MPI_Request> recvRequests(max(1, noAzimuthalNeighborDomains()), AT_, "recvRequests");
  MIntScratchSpace noWindowCellsRecv(max(1, noAzimuthalNeighborDomains()), AT_, "noWindowCellsRecv");
  for(MInt d = 0; d < noAzimuthalNeighborDomains(); d++) {
    noWindowCellsRecv[d] = -1;
    MPI_Irecv(&noWindowCellsRecv[d], 1, type_traits<MInt>::mpiType(), azimuthalNeighborDomain(d),
              azimuthalNeighborDomain(d), mpiComm(), &recvRequests[d], AT_, "noWindowCellsRecv[d]");
  }

  // Start sending number of azimuthal window cells to each azimuthal neighbor domain
  ScratchSpace<MPI_Request> sendRequests(max(1, noAzimuthalNeighborDomains()), AT_, "sendRequests");
  MIntScratchSpace noWindowCellsSend(max(1, noAzimuthalNeighborDomains()), AT_, "noWindowCellsSend");
  for(MInt d = 0; d < noAzimuthalNeighborDomains(); d++) {
    noWindowCellsSend[d] = noAzimuthalWindowCells(d);
    MPI_Isend(&noWindowCellsSend[d], 1, type_traits<MInt>::mpiType(), azimuthalNeighborDomain(d), domainId(), mpiComm(),
              &sendRequests[d], AT_, "noWindowCellsSend[d]");
  }

  // Finish MPI communication
  MPI_Waitall(noAzimuthalNeighborDomains(), &recvRequests[0], MPI_STATUSES_IGNORE, AT_);
  MPI_Waitall(noAzimuthalNeighborDomains(), &sendRequests[0], MPI_STATUSES_IGNORE, AT_);

  // Check if received number of azimuthal window cells matches the local number of azimuthal halo cells
  for(MInt d = 0; d < noAzimuthalNeighborDomains(); d++) {
    if(noWindowCellsRecv[d] != noAzimuthalHaloCells(d)) {
      TERMM(1, "Solver " + to_string(m_solverId) + ": Number of azimuthal window cells from domain "
                   + to_string(azimuthalNeighborDomain(d))
                   + " does not match local number of azimuthal halo cells; window: " + to_string(noWindowCellsRecv[d])
                   + " ,halo: " + to_string(noAzimuthalHaloCells(d)));
    }
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Check #2: grid global ids of azimuthal window/halo cells match
  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Start receiving azimuthal window cell global ids from each azimuthal neighbor domain
  fill(recvRequests.begin(), recvRequests.end(), MPI_REQUEST_NULL);
  const MInt totalNoWindowCellsRecv = accumulate(noWindowCellsRecv.begin(), noWindowCellsRecv.end(), 0);
  MIntScratchSpace windowCellsRecv(max(1, totalNoWindowCellsRecv), AT_, "windowCellsRecv");
  fill(windowCellsRecv.begin(), windowCellsRecv.end(), -1);
  for(MInt d = 0, offset = 0; d < noAzimuthalNeighborDomains(); d++) {
    if(noAzimuthalHaloCells(d) > 0) {
      MPI_Irecv(&windowCellsRecv[offset], noAzimuthalHaloCells(d), type_traits<MInt>::mpiType(),
                azimuthalNeighborDomain(d), azimuthalNeighborDomain(d), mpiComm(), &recvRequests[d], AT_,
                "windowCellsRecv[offset]");
    }
    offset += noAzimuthalHaloCells(d);
  }

  // Start sending window cell global ids to each neighbor domain
  fill(sendRequests.begin(), sendRequests.end(), MPI_REQUEST_NULL);
  const MInt totalNoWindowCellsSend = accumulate(noWindowCellsSend.begin(), noWindowCellsSend.end(), 0);
  MIntScratchSpace windowCellsSend(max(1, totalNoWindowCellsSend), AT_, "windowCellsSend");
  for(MInt d = 0, offset = 0; d < noAzimuthalNeighborDomains(); d++) {
    for(MInt c = 0; c < noWindowCellsSend[d]; c++) {
      windowCellsSend[offset + c] = raw().treeb().globalId(m_tree.solver2grid(azimuthalWindowCell(d, c)));
    }
    if(noAzimuthalWindowCells(d) > 0) {
      MPI_Isend(&windowCellsSend[offset], noAzimuthalWindowCells(d), type_traits<MInt>::mpiType(),
                azimuthalNeighborDomain(d), domainId(), mpiComm(), &sendRequests[d], AT_, "windowCellsSend[offset]");
    }
    offset += noAzimuthalWindowCells(d);
  }

  // Finish MPI communication
  MPI_Waitall(noAzimuthalNeighborDomains(), &recvRequests[0], MPI_STATUSES_IGNORE, AT_);
  MPI_Waitall(noAzimuthalNeighborDomains(), &sendRequests[0], MPI_STATUSES_IGNORE, AT_);

  // Check if received window cell global ids match the local halo cell global ids
  for(MInt d = 0, offset = 0; d < noAzimuthalNeighborDomains(); d++) {
    for(MInt c = 0; c < noAzimuthalHaloCells(d); c++) {
      const MInt cellId = m_tree.solver2grid(azimuthalHaloCell(d, c));
      const MInt globalId = raw().treeb().globalId(cellId);
      if(windowCellsRecv[offset + c] != globalId) {
        TERMM(1, "Global id of azimuthal window cell " + to_string(c) + " from domain "
                     + to_string(azimuthalNeighborDomain(d)) + " does not match local azimuthal halo cell gobal id ("
                     + to_string(windowCellsRecv[offset + c]) + " vs. "
                     + to_string(raw().treeb().globalId(m_tree.solver2grid(azimuthalHaloCell(d, c)))) + ")");
      }
    }
    offset += noAzimuthalHaloCells(d);
  }
}


template <MInt nDim>
void Proxy<nDim>::correctAzimuthalHaloCells() {
  TRACE();

  if(raw().noAzimuthalNeighborDomains() > 0) {
    MUint haloCnt = maia::mpi::getBufferSize(raw().azimuthalHaloCells());
    MUint windowCnt = maia::mpi::getBufferSize(raw().azimuthalWindowCells());
    m_isOutsideWindow.resize(windowCnt);
    m_isOutsideHalo.resize(haloCnt);

    MInt cnt = 0;
    for(MInt d = 0; d < raw().noAzimuthalNeighborDomains(); d++) {
      for(MInt c = 0; c < raw().noAzimuthalWindowCells(d); c++) {
        const MInt gridWindowId = raw().azimuthalWindowCell(d, c);
        MBool outside = checkOutsideGeometry(gridWindowId);
        if(outside) { // Fully outside solver domain
          m_isOutsideWindow[cnt] = true;
        } else {
          m_isOutsideWindow[cnt] = false;
        }
        cnt++;
      }
    }
    mpi::exchangeBuffer(raw().azimuthalNeighborDomains(), raw().azimuthalHaloCells(), raw().azimuthalWindowCells(),
                        raw().mpiComm(), m_isOutsideWindow.data(), m_isOutsideHalo.data());

    cnt = 0;
    for(MInt d = 0; d < raw().noAzimuthalNeighborDomains(); d++) {
      for(MInt c = 0; c < raw().noAzimuthalHaloCells(d); c++) {
        const MInt gridHaloId = raw().azimuthalHaloCell(d, c);
        ASSERT(raw().treeb().hasProperty(gridHaloId, Cell::IsPeriodic), "");
        MBool outside = checkOutsideGeometry(gridHaloId);
        if(outside) { // Halo is fully outside domain
          raw().treeb().solver(gridHaloId, solverId()) = false;
          m_isOutsideHalo[cnt] = true;
        } else {
          if(!raw().treeb().solver(gridHaloId, solverId())) {
            if(m_isOutsideHalo[cnt]) { // Window is not in solver domain - Make it unmapped
              raw().treeb().solver(gridHaloId, solverId()) = true;
              m_isOutsideHalo[cnt] = true;
            } else { // Window is not refined for solver
              m_isOutsideHalo[cnt] = true;
            }
          } else {
            m_isOutsideHalo[cnt] = false;
          }
        }
        cnt++;
      }
    }
    mpi::exchangeBuffer(raw().azimuthalNeighborDomains(), raw().azimuthalWindowCells(), raw().azimuthalHaloCells(),
                        raw().mpiComm(), m_isOutsideHalo.data(), m_isOutsideWindow.data());
  }

  for(MInt c = 0; c < raw().noAzimuthalUnmappedHaloCells(); c++) {
    const MInt gridHaloId = raw().azimuthalUnmappedHaloCell(c);
    ASSERT(raw().treeb().hasProperty(gridHaloId, Cell::IsPeriodic), "");
    ASSERT(raw().a_globalId(gridHaloId) == -1, "Unmapped grid halo cell " + std::to_string(gridHaloId)
                                                   + " globalId not -1: " + to_string(raw().a_globalId(gridHaloId)));
    MBool outside = checkOutsideGeometry(gridHaloId);
    if(outside) {
      raw().treeb().solver(gridHaloId, solverId()) = false;
    }
  }
  // Correct leaf level
  for(MInt c = 0; c < raw().noAzimuthalUnmappedHaloCells(); c++) {
    const MInt gridHaloId = raw().azimuthalUnmappedHaloCell(c);
    raw().a_isLeafCell(gridHaloId, solverId()) = false;
    if(raw().treeb().solver(gridHaloId, solverId())) {
      raw().a_isLeafCell(gridHaloId, solverId()) =
          !raw().a_hasChildren(gridHaloId, solverId()) && !raw().a_hasProperty(gridHaloId, Cell::IsPartLvlAncestor);
    }
  }
}

/*! \brief smooth the grid-based values over neighboring leaf cell
 *         asure conservation of the total value
 * \author ???, Tim Wegmann
 */
template <MInt nDim>
void Proxy<nDim>::smoothFilter(const MInt level, MFloat** value) {
  vector<MInt> levelCellId;
  // cout << "finding all cells on level " << level - 1 << endl;
  getLocalSameLevelCellIds(levelCellId, level - 1);

  // cerr << "test2 " << value[0][0] << endl;
  // cerr << "test3 " << value[38252][0] << endl;
  // cerr << "test "<< value[38253][0] << endl;

  for(auto& cellId : levelCellId) {
    const MInt noChilds = tree().noChildren(cellId);
    if(noChilds == 0) {
      continue;
    }

    for(MInt childId = 0; childId < noChilds; childId++) {
      const MInt noChilds2 = tree().noChildren(tree().child(cellId, childId));
      if(noChilds2 == 0) {
        continue;
      }

      for(MInt childId2 = 0; childId2 < noChilds2; childId2++) {
        for(MInt dimId = 0; dimId < nDim; dimId++) {
          const MInt childCellId2 = tree().child(tree().child(cellId, childId), childId2);
          // cerr << "cellId " << cellId << endl;
          // cerr << "childId " << tree().child(cellId, childId) << endl;
          // cerr << value[cellId][0] << endl;
          value[tree().child(cellId, childId)][dimId] += value[childCellId2][dimId] * 1.0 / noChilds2;
          //        value[cellId][dimId] += 0;
        }
      }
      for(MInt dimId = 0; dimId < nDim; dimId++) {
        value[cellId][dimId] += value[tree().child(cellId, childId)][dimId] * 1.0 / noChilds;
      }
    }

    for(MInt childId = 0; childId < noChilds; childId++) {
      const MInt noChilds2 = tree().noChildren(tree().child(cellId, childId));
      for(MInt childId2 = 0; childId2 < noChilds2; childId2++) {
        for(MInt dimId = 0; dimId < nDim; dimId++) {
          const MInt childCellId2 = tree().child(tree().child(cellId, childId), childId2);

          value[childCellId2][dimId] = 0.4 * value[cellId][dimId] + 0.4 * value[tree().child(cellId, childId)][dimId]
                                       + 0.2 * value[childCellId2][dimId];
          // value[tree().child(cellId, childId)][dimId] = value[cellId][dimId];
        }
      }
      for(MInt dimId = 0; dimId < nDim; dimId++) {
        value[tree().child(cellId, childId)][dimId] = 0;
      }
    }

    for(MInt dimId = 0; dimId < nDim; dimId++) {
      value[cellId][dimId] = 0;
    }

    //    if(domainId() == 0) cout << " cellId " << cellId << " with level " << tree().level(cellId) << " " <<
    //    value[cellId]
    //    << endl;
  }
}

/*! \brief find all cells on the same level
 *          used by the smoothFilter!
 * \author ???, Tim Wegmann
 */
template <MInt nDim>
void Proxy<nDim>::getLocalSameLevelCellIds(vector<MInt>& levelCellId, const MInt level) {
  const MInt noLocalPartitionCells = localPartitionCellOffsets(1) - localPartitionCellOffsets(0);

  // go down to searched for level using recursion
  std::function<void(MInt)> goDownAndPush = [&](const MInt cellId) {
    const MInt cellLevel = tree().level(cellId);
    if(level == cellLevel) {
      levelCellId.push_back(cellId);
    } else {
      for(MInt childId = 0; childId < tree().noChildren(cellId); childId++) {
        const MInt childCellId = tree().child(cellId, childId);
        goDownAndPush(childCellId);
      }
    }
  };

  // loop over partition cells
  for(MInt i = 0; i < noLocalPartitionCells; i++) {
    const MInt partionCellId = localPartitionCellLocalIds(i);
    const MInt partionLevel = tree().level(partionCellId);
    if(partionLevel == level) {
      levelCellId.push_back(partionCellId);
    } else if(partionLevel >= level) {
      // TODO labels:GRID replace by a more sensible check before the loop (e.g. partitionLevel....)
      TERMM(-1, "ERROR: Invalid level!");
    }

    if(!tree().hasChildren(partionCellId)) {
      continue;
    }

    // start of recursion
    for(MInt childId = 0; childId < tree().noChildren(partionCellId); childId++) {
      const MInt partitionChildId = tree().child(partionCellId, childId);
      goDownAndPush(partitionChildId);
    }
  }
}

/*! \brief gathers all leaf child cells for a given parentId
 * \author Tim Wegmann
 */
template <MInt nDim>
void Proxy<nDim>::getAllLeafChilds(const MInt parentId, vector<MInt>& levelChilds) {
  if(tree().isLeafCell(parentId)) {
    levelChilds.push_back(parentId);
    return;
  }

  for(MInt child = 0; child < m_maxNoChilds; child++) {
    if(!tree().hasChild(parentId, child)) continue;
    const MInt childId = tree().child(parentId, child);
    getAllLeafChilds(childId, levelChilds);
  }
}

/*! \brief Set solver flags for newly added solver
 * \author Thomas Hoesgen
 */
template <MInt nDim>
void Proxy<nDim>::setSolverFlagsForAddedSolver() {
  TRACE();

  MIntScratchSpace solverFlag(raw().maxNoCells(), AT_, "solverFlag");
  solverFlag.fill(-1);

  for(MInt gridId = 0; gridId < raw().noInternalCells(); gridId++) {
    if(raw().a_level(gridId) != raw().minLevel()) {
      continue;
    }
    MBool outside = checkOutsideGeometry(gridId);
    MInt flag = (!outside ? 1 : -1);
    solverFlag[gridId] = flag;
    std::vector<MInt> offsprings(raw().maxNoCells());
    raw().createLeafCellMapping(solverId(), gridId, offsprings, true);
    for(auto it = offsprings.begin(); it != offsprings.end(); it++) {
      solverFlag[*it] = flag;
    }
  }

  for(MInt gridId = 0; gridId < raw().noInternalCells(); gridId++) {
    if(solverFlag[gridId] == 1) {
      raw().treeb().solver(gridId, solverId()) = true;
    } else {
      raw().treeb().solver(gridId, solverId()) = false;
    }
  }

  maia::mpi::exchangeData(raw().neighborDomains(), raw().haloCells(), raw().windowCells(), raw().mpiComm(),
                          solverFlag.data());

  for(MInt gridId = raw().noInternalCells(); gridId < raw().maxNoCells(); gridId++) {
    if(solverFlag[gridId] == 1) {
      raw().treeb().solver(gridId, solverId()) = true;
    }
  }
}

/*! \brief updates solver-grid partition cell offsets and globalIds!
 * \author Tim Wegmann
 */
template <MInt nDim>
void Proxy<nDim>::updatePartitionCellOffsets() {
  TRACE();

  if(m_partitionCellGlobalId != nullptr) {
    mDeallocate(m_partitionCellGlobalId);
  }

  MLong localPartitionCellOffsetsGrid[3] = {static_cast<MLong>(-1), static_cast<MLong>(-1), static_cast<MLong>(-1)};

  MBool useRestartOffsets = !raw().updatedPartitionCells();
  if(raw().m_localPartitionCellOffsetsRestart[2] < 0) useRestartOffsets = false;

  if(useRestartOffsets) {
    for(MInt i = 0; i < 3; i++) {
      localPartitionCellOffsetsGrid[i] = raw().m_localPartitionCellOffsetsRestart[i];
    }
  } else {
    for(MInt i = 0; i < 3; i++) {
      localPartitionCellOffsetsGrid[i] = raw().localPartitionCellOffsets(i);
    }
  }

  mAlloc(m_partitionCellGlobalId, localPartitionCellOffsetsGrid[2], "m_partitionCellGlobalId", static_cast<MLong>(-1),
         AT_);

  MInt noLocalPartitionCellsGrid = localPartitionCellOffsetsGrid[1] - localPartitionCellOffsetsGrid[0];
  MInt noLocalPartitionCells = 0;
  MInt partitionGridCellId;
  for(MInt id = 0; id < noLocalPartitionCellsGrid; id++) {
    if(useRestartOffsets) {
      partitionGridCellId = raw().m_localPartitionCellLocalIdsRestart[id];
    } else {
      partitionGridCellId = raw().localPartitionCellLocalIds(id);
    }
    if(partitionGridCellId < 0) continue;
    if(m_tree.grid2solver(partitionGridCellId) < 0) continue;
    m_partitionCellGlobalId[noLocalPartitionCells] = m_tree.globalId(m_tree.grid2solver(partitionGridCellId));
    noLocalPartitionCells++;
  }

  MIntScratchSpace localPartitionCellCounts(noDomains(), AT_, "localPartitionCellCounts");
  MPI_Allgather(&noLocalPartitionCells, 1, MPI_INT, &localPartitionCellCounts[0], 1, MPI_INT, mpiComm(), AT_,
                "noLocalPartitionCells", "localPartitionCellCounts[0]");
  MLong noPartitionCellsGlobal = noLocalPartitionCells;
  MPI_Allreduce(MPI_IN_PLACE, &noPartitionCellsGlobal, 1, type_traits<MLong>::mpiType(), MPI_SUM, mpiComm(), AT_,
                "MPI_IN_PLACE", "noPartitionCellsGlobal");

  MInt offset = 0;
  for(MInt dId = 0; dId < domainId(); dId++) {
    offset += localPartitionCellCounts[dId];
  }

  m_localPartitionCellOffsets[0] = offset;
  m_localPartitionCellOffsets[1] = offset + noLocalPartitionCells;
  m_localPartitionCellOffsets[2] = noPartitionCellsGlobal;
}

/** \brief domain id containing a given global cell id
 *
 * \author Ansgar Niemoeller
 * \date 23.06.2014
 *
 * \param[in] globalId global cell id
 * \return domain id containing globalId, -1 if not found
 *
 */
template <MInt nDim>
MInt Proxy<nDim>::domainContainingCell(MInt globalId) const {
  MInt domain = -1;
  for(MInt domainId_ = 0; domainId_ < noDomains(); domainId_++) {
    if(globalId >= m_domainOffsets[domainId_] && globalId < m_domainOffsets[domainId_ + 1]) {
      domain = domainId_;
      break;
    }
  }
  return domain;
}

} // namespace grid
} // namespace maia

// Explicit instantiations
template class maia::grid::Proxy<2>;
template class maia::grid::Proxy<3>;
