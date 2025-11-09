// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include <cmath>
#include <fstream>
#include <iostream>
#include <numeric>

#include "COMM/mpiexchange.h"
#include "COMM/mpioverride.h"
#include "INCLUDE/maiatypes.h"
#include "IO/parallelio_hdf5.h"
#include "UTIL/kdtree.h"
#include "UTIL/pointbox.h"
#include "fvstructuredbndrycnd2d.h"
#include "fvstructuredsolver.h"
#include "fvstructuredsolver2d.h"


using namespace std;

template <MBool isRans>
StructuredBndryCnd2D<isRans>::StructuredBndryCnd2D(FvStructuredSolver<2>* solver, StructuredGrid<2>* grid)
  : StructuredBndryCnd<2>(solver, grid) {
  TRACE();

  const MLong oldAllocatedBytes = allocatedBytes();

  m_solver = static_cast<FvStructuredSolver2D*>(solver);

  // read rescaling bc properties
  if(m_solver->m_rescalingCommGrRoot >= 0) {
    m_rescalingBLT = 5.95;
    m_rescalingBLT = Context::getSolverProperty<MFloat>("rescalingBLT", m_solverId, AT_);

    cout << "Rescaling BLT: " << m_rescalingBLT << endl;
  }

  // compute sponge coefficients for each cell
  if(m_solver->m_useSponge) {
    RECORD_TIMER_START(m_solver->timer(Timers::BuildUpSponge));
    readAndDistributeSpongeCoordinates();
    RECORD_TIMER_STOP(m_solver->timer(Timers::BuildUpSponge));
  }

  printAllocatedMemory(oldAllocatedBytes, "StructuredBndryCnd2D", m_StructuredComm);
}

template <MBool isRans>
StructuredBndryCnd2D<isRans>::~StructuredBndryCnd2D() {}

template <MBool isRans>
void StructuredBndryCnd2D<isRans>::readAndDistributeSpongeCoordinates() {
  TRACE();

  vector<unique_ptr<StructuredWindowMap<nDim>>>& spongeInfo = m_solver->m_windowInfo->m_spongeInfoMap;
  MInt noSpongeInfo = spongeInfo.size(); // contains the size of the maps
  MInt memSize = 0;

  // 1)
  // allocate the space for all the coordinates of the sponge in Scratch!
  // memory will not be needed later ==> Scratch!

  // 1.1) determine the storage size
  // determine the size and to store the whole data
  for(MInt i = 0; i < noSpongeInfo; ++i) {
    MInt size = 1;
    for(MInt dim = 0; dim < nDim; ++dim) {
      size *= (spongeInfo[i]->end1[dim] - spongeInfo[i]->start1[dim] + 1);
    }
    memSize += size;
  }

  // 1.2) allocate the memory
  MFloatScratchSpace coordMem(nDim * memSize, AT_, "spongeCoordinates");
  MFloatPointerScratchSpace spongeCoords(noSpongeInfo, AT_, "spongeCoordPointer");

  MInt totMemSize = 0;
  for(MInt i = 0; i < noSpongeInfo; ++i) {
    MInt size = 1;
    for(MInt dim = 0; dim < nDim; ++dim) {
      size *= (spongeInfo[i]->end1[dim] - spongeInfo[i]->start1[dim] + 1);
    }
    spongeCoords[i] = &coordMem[totMemSize];
    totMemSize += nDim * size;
  }


  // 2)
  // we do not need the corner points but the centre coordinates of the face
  // we need to allocate the memory for the faces (==cell size) for the sponge face
  // ->determine the size and allocate the memory (again only scratch)

  // 2.1) calculate the number of cells to store.
  MInt cellmemSize = 0;
  // determine the size and to store the whole data
  for(MInt i = 0; i < noSpongeInfo; ++i) {
    MInt size = 1;
    for(MInt dim = 0; dim < nDim; ++dim) {
      if(spongeInfo[i]->end1[dim] - spongeInfo[i]->start1[dim] == 0) continue;
      size *= (spongeInfo[i]->end1[dim] - spongeInfo[i]->start1[dim]);
    }
    cellmemSize += size;
  }
  // 2.2) allocate the space for all the coordinates in Scratch!
  // memory will not be needed later!
  MFloatScratchSpace coordCellMem(nDim * cellmemSize, AT_, "spongeCellCoordinates");
  MFloatPointerScratchSpace spongeSurfCoords(noSpongeInfo, AT_, "spongeCellCoordPointer");

  MInt totCellMemSize = 0;
  for(MInt i = 0; i < noSpongeInfo; ++i) {
    MInt size = 1;
    for(MInt dim = 0; dim < nDim; ++dim) {
      if(spongeInfo[i]->end1[dim] - spongeInfo[i]->start1[dim] == 0) continue;
      size *= (spongeInfo[i]->end1[dim] - spongeInfo[i]->start1[dim]);
    }
    spongeSurfCoords[i] = &coordCellMem[totCellMemSize];
    totCellMemSize += nDim * size;
  }


  // 3) read in the coordinates of the grid points
  // open file for reading the grid data
  MString gridFileName = m_grid->m_gridInputFileName;

  // unique identifier needed to associate grid and solution in case of restart
  // const char* aUID= new char[18];


  // open file and read number of solvers and uid
  ParallelIoHdf5 pio(gridFileName, maia::parallel_io::PIO_READ, m_StructuredComm);

  // read the data in and distribute ist

  // the split of reading and distributing is done on purpose!
  // this is because we then know what the library is doing
  // and not what the io_library such as hdf or netcdf should do
  // but does not
  // a good library should handle that directly! TEST IT
  memSize = 1;
  // read in the data if  processor zero else read nothing!!!
  ParallelIo::size_type ioOffset[2] = {0, 0};
  ParallelIo::size_type ioSize[2] = {0, 0};
  if(m_solver->domainId() == 0) {
    for(MInt i = 0; i < noSpongeInfo; ++i) {
      memSize = 1;
      for(MInt dim = 1; dim >= 0; --dim) {
        ioSize[dim] = (spongeInfo[i]->end1[1 - dim] - spongeInfo[i]->start1[1 - dim] + 1);
        memSize *= ioSize[dim];
        ioOffset[dim] = spongeInfo[i]->start1[1 - dim];
      }
      // read in the data if  processor zero else read nothing!!!
      // determine the Solver name
      MString bName = "block";
      stringstream number;
      number << spongeInfo[i]->Id1;
      bName += number.str();
      pio.readArray(&spongeCoords[i][0], bName, "x", 2, ioOffset, ioSize);
      pio.readArray(&spongeCoords[i][memSize], bName, "y", 2, ioOffset, ioSize);
    }
  } else {
    for(MInt i = 0; i < noSpongeInfo; ++i) {
      MString bName = "block";
      stringstream number;
      number << spongeInfo[i]->Id1;
      bName += number.str();
      MFloat empty = 0;
      pio.readArray(&empty, bName, "x", 2, ioOffset, ioSize);
      pio.readArray(&empty, bName, "y", 2, ioOffset, ioSize);
    }
  }

  // now broadcast the information to everyone!!!
  MPI_Bcast(&spongeCoords[0][0], totMemSize, MPI_DOUBLE, 0, m_StructuredComm, AT_, "spongeCoords[0][0]");

  // 4) computing the coordinates of surface center from corner points;

  for(MInt ii = 0; ii < noSpongeInfo; ++ii) {
    MInt label, size1, count = 0;
    for(label = 0; label < nDim; label++) { // 2== dimensions
      if(spongeInfo[ii]->end1[label] - spongeInfo[ii]->start1[label] == 0) break;
    }
    switch(label) {
      case 0: {
        size1 = spongeInfo[ii]->end1[1] - spongeInfo[ii]->start1[1] + 1;
        for(MInt i = 0; i < size1 - 1; i++) {
          MInt I = i;
          MInt IP = i + 1;
          spongeSurfCoords[ii][count] = 0.5 * (spongeCoords[ii][I] + spongeCoords[ii][IP]);
          spongeSurfCoords[ii][count + (size1 - 1)] =
              0.5 * (spongeCoords[ii][I + size1] + spongeCoords[ii][IP + size1]);
          count++;
        }
        break;
      }
      case 1: {
        size1 = spongeInfo[ii]->end1[0] - spongeInfo[ii]->start1[0] + 1;
        for(MInt i = 0; i < size1 - 1; i++) {
          MInt I = i;
          MInt IP = i + 1;
          spongeSurfCoords[ii][count] = 0.5 * (spongeCoords[ii][I] + spongeCoords[ii][IP]);
          spongeSurfCoords[ii][count + (size1 - 1)] =
              0.5 * (spongeCoords[ii][I + size1] + spongeCoords[ii][IP + size1]);
          count++;
        }
        break;
      }
      default:
        mTerm(1, AT_, "sponge direction is messed up");
    }
  }

  // now everyone can determine the shortest distance

  m_log << "sponge layer build: searching for the nearest points (building simgma sponge)" << endl;
  // cout << "seraching for the nearest points(building sigma sponge)" << endl;
  // build a k-d-tree for a quick search:
  // 1) rearragne the coordinates into points;
  for(MInt i = 0; i < noSpongeInfo; ++i) {
    MInt noPoints = 1;
    for(MInt dim = 0; dim < nDim; ++dim) {
      if(spongeInfo[i]->end1[dim] - spongeInfo[i]->start1[dim] == 0) continue;
      noPoints *= (spongeInfo[i]->end1[dim] - spongeInfo[i]->start1[dim]);
    }
    // build up the points (surface centres)
    vector<Point<2>> pts;
    for(MInt j = 0; j < noPoints; ++j) {
      Point<2> a(spongeSurfCoords[i][j], spongeSurfCoords[i][j + noPoints]);
      pts.push_back(a);
    }

    // build up the tree
    KDtree<2> tree(pts);
    MFloat distance = -1.0;
    // go through all the cells an determine the closest distance
    MFloat spongeThickness = spongeInfo[i]->spongeThickness;
    for(MInt id = 0; id < m_noCells; ++id) {
      distance = -1.1111111111111111; // to check
      Point<2> pt(m_cells->coordinates[0][id], m_cells->coordinates[1][id]);
      (void)tree.nearest(pt, distance);
      if(distance <= spongeThickness) {
        MFloat spongeFactor =
            spongeInfo[i]->sigma * pow((spongeThickness - distance) / spongeThickness, spongeInfo[i]->beta);
        m_cells->fq[FQ->SPONGE_FACTOR][id] = mMax(m_cells->fq[FQ->SPONGE_FACTOR][id], spongeFactor);
      }
    }
  }

  m_log << "Spone layer SUCESSFUL: finished building sigma sponge " << endl;
}

template <MBool isRans>
void StructuredBndryCnd2D<isRans>::updateSpongeLayer() {
  // damp to infinity values
  // for the moment only rho and rhoE are damped
  const MFloat gammaMinusOne = m_solver->m_gamma - 1.0;
  switch(m_spongeLayerType) {
    case 1: {
      MFloat deltaRhoE = F0, deltaRho = F0;
      for(MInt j = m_noGhostLayers; j < m_nCells[0] - m_noGhostLayers; j++) {
        for(MInt i = m_noGhostLayers; i < m_nCells[1] - m_noGhostLayers; i++) {
          // compute the forcing term
          // for the pressure or engery and the velocity
          MInt cellId = cellIndex(i, j);
          const MFloat rhoE =
              m_cells->pvariables[PV->P][cellId] / gammaMinusOne
              + F1B2 * m_cells->pvariables[PV->RHO][cellId]
                    * (POW2(m_cells->pvariables[PV->U][cellId]) + POW2(m_cells->pvariables[PV->V][cellId]));

          deltaRhoE = rhoE - CV->rhoEInfinity;
          deltaRho = m_cells->pvariables[PV->RHO][cellId] - (CV->rhoInfinity * m_targetDensityFactor);

          m_cells->rightHandSide[CV->RHO_E][cellId] -= m_cells->fq[FQ->SPONGE_FACTOR][cellId] * deltaRhoE
                                                       * m_cells->cellJac[cellId]; // deltaP * m_cells->cellJac[cellId];
          m_cells->rightHandSide[CV->RHO][cellId] -=
              m_cells->fq[FQ->SPONGE_FACTOR][cellId] * deltaRho * m_cells->cellJac[cellId];
        }
      }
      break;
    }
    case 2: {
      // damp to values in the FQ field (set at startup, from RANS etc.)
      // damp to values in the FQ field (set at startup, from RANS etc.)
      MFloat deltaRhoE = F0, deltaRho = F0;
      for(MInt j = m_noGhostLayers; j < m_nCells[0] - m_noGhostLayers; j++) {
        for(MInt i = m_noGhostLayers; i < m_nCells[1] - m_noGhostLayers; i++) {
          MInt cellId = cellIndex(i, j);
          const MFloat rhoE =
              m_cells->pvariables[PV->P][cellId] / gammaMinusOne
              + F1B2 * m_cells->pvariables[PV->RHO][cellId]
                    * (POW2(m_cells->pvariables[PV->U][cellId]) + POW2(m_cells->pvariables[PV->V][cellId]));

          deltaRhoE = rhoE - m_cells->fq[FQ->SPONGE_RHO_E][cellId];
          deltaRho =
              m_cells->pvariables[PV->RHO][cellId] - (m_cells->fq[FQ->SPONGE_RHO][cellId] * m_targetDensityFactor);
          m_cells->rightHandSide[CV->RHO_E][cellId] -=
              m_cells->fq[FQ->SPONGE_FACTOR][cellId] * deltaRhoE * m_cells->cellJac[cellId];
          m_cells->rightHandSide[CV->RHO][cellId] -=
              m_cells->fq[FQ->SPONGE_FACTOR][cellId] * deltaRho * m_cells->cellJac[cellId];
        }
      }

      break;
    }
    case 6: {
      // damp to PInfinity
      const MFloat FgammaMinusOne = F1 / (m_solver->m_gamma - F1);
      for(MInt j = m_noGhostLayers; j < m_nCells[0] - m_noGhostLayers; j++) {
        for(MInt i = m_noGhostLayers; i < m_nCells[1] - m_noGhostLayers; i++) {
          const MInt cellId = cellIndex(i, j);
          const MFloat deltaP = (m_cells->pvariables[PV->P][cellId] - PV->PInfinity) * FgammaMinusOne;
          m_cells->rightHandSide[CV->RHO_E][cellId] -=
              m_cells->fq[FQ->SPONGE_FACTOR][cellId] * deltaP * m_cells->cellJac[cellId];
        }
      }
      break;
    }
    default:
      mTerm(1, AT_, "Sponge type doesn't exist");
  }
}


//>marian: new function to compute distance to nearest wall for each cell
template <MBool isRans>
void StructuredBndryCnd2D<isRans>::computeWallDistances() {
  vector<unique_ptr<StructuredWindowMap<nDim>>>& wallDistInfo = m_solver->m_windowInfo->m_wallDistInfoMap;
  MInt noWallDistInfo = wallDistInfo.size(); // contains the size of the maps
  MInt memSize = 0;

  // initialize array with high numbers
  for(MInt id = 0; id < m_noCells; ++id) {
    m_cells->fq[FQ->WALLDISTANCE][id] = 99999;
  }

  // 1)
  // memory will not be needed later ==> Scratch!

  // 1.1) determine the storage size
  // determine the size and to store the whole data
  for(MInt i = 0; i < noWallDistInfo; ++i) {
    MInt size = 1;
    for(MInt dim = 0; dim < 2; ++dim) {
      size *= (wallDistInfo[i]->end1[dim] - wallDistInfo[i]->start1[dim] + 1);
    }
    memSize += size;
  }

  // 1.2) allocate the memory
  MFloatScratchSpace coordMem(2 * memSize, AT_, "wallCoordinates");
  MFloatPointerScratchSpace wallDistCoords(noWallDistInfo, AT_, "wallCoordsPointer");

  MInt totMemSize = 0;
  for(MInt i = 0; i < noWallDistInfo; ++i) {
    MInt size = 1;
    for(MInt dim = 0; dim < 2; ++dim) {
      size *= (wallDistInfo[i]->end1[dim] - wallDistInfo[i]->start1[dim] + 1);
    }
    wallDistCoords[i] = &coordMem[totMemSize];
    totMemSize += 2 * size;
  }

  // 2)
  // we do not need the corner points but the centre coordinates of the face
  // we need to allocate the memory for the faces (==cell size) for the sponge face
  // ->determine the size and allocate the memory (again only scratch)

  // 2.1) calculate the number of cells to store.
  MInt cellmemSize = 0;
  // determine the size and to store the whole data
  for(MInt i = 0; i < noWallDistInfo; ++i) {
    MInt size = 1;
    for(MInt dim = 0; dim < 2; ++dim) {
      if(wallDistInfo[i]->end1[dim] - wallDistInfo[i]->start1[dim] == 0) continue;
      size *= (wallDistInfo[i]->end1[dim] - wallDistInfo[i]->start1[dim]);
    }
    cellmemSize += size;
  }
  // 2.2) allocate the space for all the coordinates in Scratch!
  // memory will not be needed later!
  MFloatScratchSpace coordCellMem(2 * cellmemSize, AT_, "wallDistCellCoordinates");
  MFloatPointerScratchSpace wallDistSurfCoords(noWallDistInfo, AT_, "wallDistCellCoordPointer");

  MInt totCellMemSize = 0;
  for(MInt i = 0; i < noWallDistInfo; ++i) {
    MInt size = 1;
    for(MInt dim = 0; dim < 2; ++dim) {
      if(wallDistInfo[i]->end1[dim] - wallDistInfo[i]->start1[dim] == 0) continue;
      size *= (wallDistInfo[i]->end1[dim] - wallDistInfo[i]->start1[dim]);
    }
    wallDistSurfCoords[i] = &coordCellMem[totCellMemSize];
    totCellMemSize += 2 * size;
  }

  // 3) read in the coordinates of the grid points
  // open file for reading the grid data
  MString gridFileName = m_grid->m_gridInputFileName;

  // open file and read number of solvers and uid
  ParallelIoHdf5 pio(gridFileName, maia::parallel_io::PIO_READ, m_solver->m_StructuredComm);

  // read the data in and distribute ist

  // the split of reading and distributing is done on purpose!
  // this is because we then know what the library is doing
  // and not what the io_library such as hdf or netcdf should do
  // but does not
  // a good library should handle that directly! TEST IT
  memSize = 1;
  // read in the data if  processor zero else read nothing!!!
  ParallelIo::size_type ioOffset[2] = {0, 0};
  ParallelIo::size_type ioSize[2] = {0, 0};
  if(m_solver->domainId() == 0) {
    for(MInt i = 0; i < noWallDistInfo; ++i) {
      memSize = 1;
      for(MInt dim = 1; dim >= 0; --dim) {
        ioSize[dim] = (wallDistInfo[i]->end1[1 - dim] - wallDistInfo[i]->start1[1 - dim] + 1);
        memSize *= ioSize[dim];
        ioOffset[dim] = wallDistInfo[i]->start1[1 - dim];
      }
      // read in the data if  processor zero else read nothing!!!
      // determine the Solver name
      MString bName = "block";
      stringstream number;
      number << wallDistInfo[i]->Id1;
      bName += number.str();
      pio.readArray(&wallDistCoords[i][0], bName, "x", 2, ioOffset, ioSize);
      pio.readArray(&wallDistCoords[i][memSize], bName, "y", 2, ioOffset, ioSize);
    }
  } else {
    for(MInt i = 0; i < noWallDistInfo; ++i) {
      MString bName = "block";
      stringstream number;
      number << wallDistInfo[i]->Id1;
      bName += number.str();
      MFloat empty = 0;
      pio.readArray(&empty, bName, "x", 2, ioOffset, ioSize);
      pio.readArray(&empty, bName, "y", 2, ioOffset, ioSize);
    }
  }

  // cout << "broadcasting wall-bc information" << endl;
  // now broadcast the information to everyone!!!
  MPI_Bcast(&wallDistCoords[0][0], totMemSize, MPI_DOUBLE, 0, m_solver->m_StructuredComm, AT_, "wallDistCoords[0][0]");
  // cout << "broadcast end" << endl;

  // 4) computing the coordinates of surface center from corner points;

  for(MInt ii = 0; ii < noWallDistInfo; ++ii) {
    MInt label, size1, count = 0;
    for(label = 0; label < 2; label++) { // 2== dimensions
      if(wallDistInfo[ii]->end1[label] - wallDistInfo[ii]->start1[label] == 0) break;
    }
    switch(label) {
      case 0: {
        size1 = wallDistInfo[ii]->end1[1] - wallDistInfo[ii]->start1[1] + 1;
        for(MInt i = 0; i < size1 - 1; i++) {
          MInt I = i;
          MInt IP = i + 1;
          wallDistSurfCoords[ii][count] = 0.5 * (wallDistCoords[ii][I] + wallDistCoords[ii][IP]);
          wallDistSurfCoords[ii][count + (size1 - 1)] =
              0.5 * (wallDistCoords[ii][I + size1] + wallDistCoords[ii][IP + size1]);
          count++;
        }
        break;
      }
      case 1: {
        size1 = wallDistInfo[ii]->end1[0] - wallDistInfo[ii]->start1[0] + 1;
        for(MInt i = 0; i < size1 - 1; i++) {
          MInt I = i;
          MInt IP = i + 1;
          wallDistSurfCoords[ii][count] = 0.5 * (wallDistCoords[ii][I] + wallDistCoords[ii][IP]);
          wallDistSurfCoords[ii][count + (size1 - 1)] =
              0.5 * (wallDistCoords[ii][I + size1] + wallDistCoords[ii][IP + size1]);
          count++;
        }
        break;
      }
      default:
        mTerm(1, AT_, "wall direction is messed up");
    }
  }

  // now everyone can determine the shortest distance

  m_log << "wall distance computation: searching for the nearest wall" << endl;
  // cout << "seraching for the nearest" << endl;
  // build a k-d-tree for a quick search:
  // 1) rearragne the coordinates into points;
  for(MInt i = 0; i < noWallDistInfo; ++i) {
    MInt noPoints = 1;
    for(MInt dim = 0; dim < 2; ++dim) {
      if(wallDistInfo[i]->end1[dim] - wallDistInfo[i]->start1[dim] == 0) continue;
      noPoints *= (wallDistInfo[i]->end1[dim] - wallDistInfo[i]->start1[dim]);
    }
    // build up the points (surface centres)
    vector<Point<2>> pts;
    for(MInt j = 0; j < noPoints; ++j) {
      Point<2> a(wallDistSurfCoords[i][j], wallDistSurfCoords[i][j + noPoints]);
      pts.push_back(a);
    }
    // build up the tree
    KDtree<2> tree(pts);
    MFloat distance = -1.0;

    // go through all the cells an determine the closest distance
    for(MInt id = 0; id < m_noCells; ++id) {
      distance = -1.1111111111111111; // to check
      Point<2> pt(m_cells->coordinates[0][id], m_cells->coordinates[1][id]);
      (void)tree.nearest(pt, distance);

      // take minimum because another bc1000 might be further away than the current but would overwrite the actually
      // closer distance
      m_cells->fq[FQ->WALLDISTANCE][id] = mMin(m_cells->fq[FQ->WALLDISTANCE][id], distance);
    }
  }

  m_log << "Wall Distance Computation SUCESSFUL: Saved minimum distance to next wall for all cells " << endl;
}
//<marian


// This function does the same as computeWallDistance(); But instead of reading the wall surface
// coordinates from the grid file, the wall surface coordinates are gathered from all ranks, i.e.
// in a first step each rank loops over m_auxDataMap, which is supposed to contain all wall maps;
// Hence the call to FvStructuredSolverWindowInfo<nDim>::setWallInformation is unnecessary.
// By this approach we can store for each cell its nearest wall neighbor cell. This information can
// be utilized to compute, e.g., the wall distance in units of y+ (i.e. we need to communicate the
// friction velocity during the simulation)
/*template <MBool isRans>
void StructuredBndryCnd2D<isRans>::computeLocalWallDistances(){

  // Initialize array with high numbers
  MFloat maxWallDistance = numeric_limits<MFloat>::max();
  maxWallDistance = Context::getSolverProperty<MFloat>("maxWallDistance", m_solverId, AT_, &maxWallDistance);
  for(MInt id=0; id<m_noCells; ++id){
    m_cells->fq[FQ->WALLDISTANCE][id] = maxWallDistance;
  }

  const vector<StructuredWindowMap*> wallDistInfo = m_auxDataMap;//m_solver->m_windowInfo->m_wallDistInfoMap;
  const MInt noWallDistInfo = wallDistInfo.size(); //contains the size of the maps
  std::vector<MInt> wallNormalIndex(noWallDistInfo);
  std::vector<MInt> normalDir(noWallDistInfo);
  std::vector<std::array<MInt, 2*nDim>> startEndTangential(noWallDistInfo);
  // Number of cells on current rank to store
  MInt cellmemSize = 0;
  for (MInt nn=0; nn<noWallDistInfo;++nn) {
    const MInt* const start = wallDistInfo[nn]->start1;
    const MInt* const end = wallDistInfo[nn]->end1;
    const MInt face = wallDistInfo[nn]->face;
    normalDir[nn] = face/2;
    wallNormalIndex[nn] = start[normalDir[nn]];

    MInt size = 1;
    for (MInt dim = 0; dim < nDim-1; ++dim) {
      const MInt tangentialDir = (normalDir[nn]+(dim+1))%nDim;
      startEndTangential[nn][dim*2+0] = start[tangentialDir];// + m_noGhostLayers;
      startEndTangential[nn][dim*2+1] = end[tangentialDir];// - m_noGhostLayers;
      size *= startEndTangential[nn][dim*2+1] - startEndTangential[nn][dim*2+0];
    }
    cellmemSize += size;
#ifndef NDEBUG
    cout << "Debug output: wallDistInfo=" << nn << ": face=" << face << "; normalDir="
         << normalDir[nn] << "; wallNormalIndex=" << wallNormalIndex[nn] << "; startEndTangential="
         << startEndTangential[nn][0] << "|" << startEndTangential[nn][1] << "normal start/end=" << start[normalDir[nn]]
<< "|" << end[normalDir[nn]] << endl; #endif
  }

  //2.2) allocate the space for all the coordinates in Scratch!
  //memory will not be needed later!
  MFloatScratchSpace coordCellMem(std::max(1,nDim*cellmemSize), AT_, "coordCellMem" );
  MFloatPointerScratchSpace wallDistSurfCoords(nDim, AT_, "wallDistSurfCoords");
  for (MInt dim = 0; dim < nDim; ++dim)
    wallDistSurfCoords[dim] = &coordCellMem[dim*cellmemSize];

  //4) computing the coordinates of surface center from corner points;
  MInt cnt = 0;
  for (MInt nn=0; nn<noWallDistInfo; ++nn) {
    MInt tIdx;
    MInt& i = (normalDir[nn]==0) ? wallNormalIndex[nn] : tIdx;
    MInt& j = (normalDir[nn]==0) ? tIdx : wallNormalIndex[nn];
    MInt inc[] = {0, 0};
    inc[1-normalDir[nn]] = 1;
    for (tIdx = startEndTangential[nn][0]; tIdx < startEndTangential[nn][1]; ++tIdx) {
      const MInt p1 = getPointIdFromCell(i, j);
      const MInt p2 = getPointIdFromPoint(p1, inc[0], inc[1]);
      wallDistSurfCoords[0][cnt] = 0.5*(m_grid->m_coordinates[0][p1] + m_grid->m_coordinates[0][p2]);
      wallDistSurfCoords[1][cnt] = 0.5*(m_grid->m_coordinates[1][p1] + m_grid->m_coordinates[1][p2]);
      ++cnt;
    }
  }

  // Broadcast all wall surface coordinates to all ranks
  // Note: actually this approach does not determine the absolute shortest
  // distance to a wall, but only the shortest distance among all wall surface centroids
  MInt noRanks;
  MPI_Comm_size(m_solver->m_StructuredComm, &noRanks);
  // recvcounts: How many wall surface coordinates to receive from each domain
  std::vector<MInt> recvcounts(noRanks);
  MPI_Allgather(&cellmemSize, 1, MPI_INT, recvcounts.data(), 1, MPI_INT, m_solver->m_StructuredComm, AT_,
                "cellmemSize", "recvcounts");
  std::vector<MInt> displs(noRanks);
  for (MInt r = 0; r < noRanks-1; ++r) {
    displs[r+1] = displs[r] + recvcounts[r];
  }
  const MInt noTotalWallPoints = displs[noRanks-1] + recvcounts[noRanks-1];
  // wallDistSurfCoordsAll: contains the wall coordinates of all wall surfaces of all ranks
  std::vector<std::vector<MFloat>> wallDistSurfCoordsAll(nDim);
  for (MInt dim = 0; dim < nDim; ++dim) {
    wallDistSurfCoordsAll[dim].resize(noTotalWallPoints);
    MPI_Allgatherv(&wallDistSurfCoords[dim][0], cellmemSize, MPI_DOUBLE,
        wallDistSurfCoordsAll[dim].data(), recvcounts.data(), displs.data(), MPI_DOUBLE,
        m_solver->m_StructuredComm, AT_, "wallDistSurfCoords", "wallDistSurfCoordsAll");
  }

  // --> Now all ranks have all wall surface coordinates

  // globalWallIds: store the positions of all wall surfaces in wallDistSurfCoordsAll, which are relevant
  //                for current rank
  std::set<MInt> globalWallIds;
  // cellId2globalWallId: mapping of cellId to the position of its closest wall surface in wallDistSurfCoordsAll
  std::vector<MInt> cellId2globalWallId(m_noCells, -1);

  // Build k-d-tree for a quick search
  // In the following Point<3> is used, but 3rd dimension is just dummy
  vector < Point<3> > pts;
  for(MInt globalWallId=0; globalWallId<noTotalWallPoints; ++globalWallId){
    Point<3> a(wallDistSurfCoordsAll[0][globalWallId],wallDistSurfCoordsAll[1][globalWallId], 0, globalWallId);
    pts.push_back(a);
  }
  //build up the tree
  KDtree<3> tree(pts);
  MFloat distance = -1.0;

  //go through all the cells an determine the closest distance
  for (MInt cellId = 0; cellId < m_noCells; ++cellId) {
    distance=-1.1111111111111111; //to check
    Point<3> pt(m_cells->coordinates[0][cellId],m_cells->coordinates[1][cellId], 0);
    const MInt globalWallId = tree.nearest(pt, distance);

    if (distance < m_cells->fq[FQ->WALLDISTANCE][cellId]) {
      globalWallIds.insert(globalWallId);
      cellId2globalWallId[cellId] = globalWallId;
      m_cells->fq[FQ->WALLDISTANCE][cellId] = distance;
    }
  }

  //
  auto getDomainId = [&displs, noRanks] (const MInt id, const MInt hint = 0) {
    for (MInt rank = hint; rank < (signed)displs.size()-1; ++rank) {
      if (displs[rank+1]>id)
        return rank;
    }
    return noRanks - 1;
  };

  // globalWallId2localId: mapping of globalWallId in wallDistSurfCoordsAll to local wallId on current rank
  std::vector<MInt> globalWallId2localId((globalWallIds.empty()) ? 0 : *globalWallIds.rbegin()+1);
  // recvWallIds: local wallIds from each rank, requested by current rank
  std::vector<MInt> recvWallIds(globalWallIds.size());
  std::vector<MInt> sendcounts(noRanks);
  MInt domainId_temp = 0;
  MInt localWallId = 0;
  for (auto it = globalWallIds.begin(); it != globalWallIds.end(); ++it, ++localWallId) {
    globalWallId2localId[*it] = localWallId;
    domainId_temp = getDomainId(*it, domainId_temp);
    // transform globalWallId "*it" to local wallId of rank domainId_temp
    recvWallIds[localWallId] = *it-displs[domainId_temp];
    ++sendcounts[domainId_temp];
  }

  // m_wallCellId2recvCell: mapping from local cellId to the position, where its closest wall surface information is
stored for (MInt cellId = 0; cellId < m_noCells; ++cellId) { if (cellId2globalWallId[cellId]==-1) continue;
    m_wallCellId2recvCell.insert({cellId, globalWallId2localId[cellId2globalWallId[cellId]]});
  }

  // Sanity check
  if (m_wallCellId2recvCell.size()<globalWallIds.size())
    mTerm(1, "Something is inconsistent!");

  // Broadcast send/receive infos
  std::vector<MInt> snghbrs(noRanks);
  std::iota(snghbrs.begin(), snghbrs.end(), 0);
  // m_wallSendcounts: # cells to be sent to each rank
  m_wallSendcounts.resize(noRanks);
  // m_wallSendCells: local wallIds to be sent to each rank
  m_wallSendCells = maia::mpi::mpiExchangePointToPoint(&recvWallIds[0],
                                                      &snghbrs[0],
                                                      noRanks,
                                                      &sendcounts[0],
                                                      &snghbrs[0],
                                                      noRanks,
                                                      m_solver->m_StructuredComm,
                                                      m_solver->domainId(),
                                                      1,
                                                      m_wallSendcounts.data());

  // --> m_wallCellId2recvCell, m_wallSendCells and m_wallSendcounts can be used later to get information
  //     from nearest wall neighbor, e.g., friction velocity at nearest wall cell

  // Sanity checks
  //1)
  if ((signed)m_wallSendCells.size()!=std::accumulate(&m_wallSendcounts[0], &m_wallSendcounts[0]+noRanks, 0))
    mTerm(1, "Neeeeeeiiiiiiiiin!");

  //2)
  if (cellmemSize!=0) {
    MInt min = std::numeric_limits<MInt>::max();
    MInt max = -1;
    for (auto id : m_wallSendCells) {
      if (id>=cellmemSize)
        mTerm(1, "Something went wrong!");
      min = mMin(id, min);
      max = mMax(id, max);
    }
    if (min!=0) {
      mTerm(1, "Scheisse: min!=0");
    }
    if (max!=cellmemSize-1)
      mTerm(1, "Scheisse ......");
  }

  // Some debug output
#ifndef NDEBUG
  cout << "::computeLocalWallDistance: rank=" << m_solver->domainId() << " cellmemSize="
       << cellmemSize << " noTotalWallPoints=" << noTotalWallPoints << " #cells near wall="
       << m_wallCellId2recvCell.size() << " #globalWallIds=" << globalWallIds.size() << endl;
#endif
}*/


// Computes shortest distance to wall, which is described by the grid as the sum of linear line elements
// and stores the information required to exchange interpolated data from the nearest neighbors to a cell
// NOTE: A call to this function is now replaced by a call to computeLocalExtendedDistancesAndSetComm() ->
//       a call to computeLocalWallDistances() and a call to computeLocalExtendedDistancesAndSetComm(), where
//       m_auxDataMap contains only wall maps should be equivalent
template <MBool isRans>
void StructuredBndryCnd2D<isRans>::computeLocalWallDistances() {
  // Initialize array with high numbers
  MFloat maxWallDistance = numeric_limits<MFloat>::max();
  maxWallDistance = Context::getSolverProperty<MFloat>("maxWallDistance", m_solverId, AT_, &maxWallDistance);
  for(MInt id = 0; id < m_noCells; ++id) {
    m_cells->fq[FQ->WALLDISTANCE][id] = maxWallDistance;
  }

  const vector<unique_ptr<StructuredWindowMap<nDim>>>& wallDistInfo =
      m_auxDataMap;                                                   // m_solver->m_windowInfo->m_wallDistInfoMap;
  const MInt noWallDistInfo = wallDistInfo.size();                    // contains the size of the maps
  std::vector<MInt> wallNormalIndex(noWallDistInfo);
  std::vector<MInt> normalDir(noWallDistInfo);
  std::vector<std::array<MInt, 2 * nDim>> startEndTangential(noWallDistInfo);
  // Number of grid points on current rank to store
  MInt cellmemSize = 0;
  for(MInt nn = 0; nn < noWallDistInfo; ++nn) {
    // start & end represent exactly the range of grid point indices, which form the map/boundary
    const MInt* const start = wallDistInfo[nn]->start1;
    const MInt* const end = wallDistInfo[nn]->end1;
    const MInt face = wallDistInfo[nn]->face;
    normalDir[nn] = face / 2;
    wallNormalIndex[nn] = start[normalDir[nn]];

    MInt size = 1;
    for(MInt dim = 0; dim < nDim - 1; ++dim) {
      const MInt tangentialDir = (normalDir[nn] + (dim + 1)) % nDim;
      startEndTangential[nn][dim * 2 + 0] = start[tangentialDir]; // + m_noGhostLayers;
      startEndTangential[nn][dim * 2 + 1] =
          end[tangentialDir] + 1; //+1, because grid coords not cell coords// - m_noGhostLayers;
      size *= startEndTangential[nn][dim * 2 + 1] - startEndTangential[nn][dim * 2 + 0];
    }
    cellmemSize += size;
#ifndef NDEBUG
    cout << "Debug output: wallDistInfo=" << nn << ": face=" << face << "; normalDir=" << normalDir[nn]
         << "; wallNormalIndex=" << wallNormalIndex[nn] << "; startEndTangential=" << startEndTangential[nn][0] << "|"
         << startEndTangential[nn][1] << "normal start/end=" << start[normalDir[nn]] << "|" << end[normalDir[nn]]
         << endl;
#endif
  }

  // 2.2) allocate the space for all the wall grid coordinates in Scratch!
  // memory will not be needed later!
  MFloatScratchSpace coordGridMem(std::max(1, nDim * cellmemSize), AT_, "coordGridMem");
  MFloatPointerScratchSpace wallDistSurfCoords(nDim, AT_, "wallDistSurfCoords");
  for(MInt dim = 0; dim < nDim; ++dim)
    wallDistSurfCoords[dim] = &coordGridMem[dim * cellmemSize];
  // Store if a grid point has a neighbor to the left or right on the same domain
  MIntScratchSpace isStartEndPoint(std::max(1, cellmemSize), AT_, "isStartEndPoint");

  // 4) Put wall grid points and information about existence of neighbor grid point in temporary buffers;
  MInt cnt = 0;
  for(MInt nn = 0; nn < noWallDistInfo; ++nn) {
    MInt tIdx;
    MInt& i = (normalDir[nn] == 0) ? wallNormalIndex[nn] : tIdx;
    MInt& j = (normalDir[nn] == 0) ? tIdx : wallNormalIndex[nn];
    for(tIdx = startEndTangential[nn][0]; tIdx < startEndTangential[nn][1]; ++tIdx) {
      const MInt p = getPointIdFromCell(i, j);
      wallDistSurfCoords[0][cnt] = m_grid->m_coordinates[0][p];
      wallDistSurfCoords[1][cnt] = m_grid->m_coordinates[1][p];
      // The case no neighbor to the left and right can not occur
      if(tIdx == startEndTangential[nn][0])
        isStartEndPoint[cnt] = 1; // no neighbor to the left on this rank
      else if(tIdx == startEndTangential[nn][1] - 1)
        isStartEndPoint[cnt] = -1; // no neighbor to the right on this rank
      else
        isStartEndPoint[cnt] =
            0; // TODO_SS labels:FV,GRID,totest is this necessary, or will it be default initialized anyway?
      ++cnt;
    }
  }

  // Broadcast all wall grid coordinates and isStartEndPoint to all ranks
  MInt noRanks;
  MPI_Comm_size(m_solver->m_StructuredComm, &noRanks);
  // recvcounts: How many wall grid coordinates to receive from each domain
  std::vector<MInt> recvcounts(noRanks);
  MPI_Allgather(&cellmemSize, 1, MPI_INT, recvcounts.data(), 1, MPI_INT, m_solver->m_StructuredComm, AT_, "cellmemSize",
                "recvcounts");
  std::vector<MInt> displs(noRanks);
  for(MInt r = 0; r < noRanks - 1; ++r) {
    displs[r + 1] = displs[r] + recvcounts[r];
  }
  const MInt noTotalWallPoints = displs[noRanks - 1] + recvcounts[noRanks - 1];
  // wallDistSurfCoordsAll: contains the wall grid coordinates of all ranks
  std::vector<std::vector<MFloat>> wallDistSurfCoordsAll(nDim);
  for(MInt dim = 0; dim < nDim; ++dim) {
    wallDistSurfCoordsAll[dim].resize(noTotalWallPoints);
    MPI_Allgatherv(&wallDistSurfCoords[dim][0], cellmemSize, MPI_DOUBLE, wallDistSurfCoordsAll[dim].data(),
                   recvcounts.data(), displs.data(), MPI_DOUBLE, m_solver->m_StructuredComm, AT_, "wallDistSurfCoords",
                   "wallDistSurfCoordsAll");
  }
  std::vector<MInt> isStartEndPointAll(noTotalWallPoints);
  MPI_Allgatherv(&isStartEndPoint[0], cellmemSize, MPI_INT, isStartEndPointAll.data(), recvcounts.data(), displs.data(),
                 MPI_INT, m_solver->m_StructuredComm, AT_, "isStartEndPoint", "isStartEndPointAll");

  // --> Now all ranks have all wall grid coordinates and the information if that grid point has any
  //     neighbor to the left or right on same domain

  // globalWallIds: store the positions of all wall grid points in wallDistSurfCoordsAll, which are relevant
  //                for current rank
  std::set<MInt> globalWallIds;
  // cellId2globalWallId: mapping of cellId to the positions of the two cells which encloses the nearest point at the
  // wall
  std::vector<std::pair<MInt, MInt>> cellId2globalWallId(m_noCells, std::make_pair(-1, -1));
  // weighting: weighting of the wall cell which is closest to a given cell; the weighting of its neighbor is
  // 1-weighting
  std::vector<MFloat> weighting(m_noCells, 0.0);

  // Build k-d-tree for a quick search
  // In the following Point<3> is used, but 3rd dimension is just dummy
  std::vector<Point<3>> pts;
  // cellIdShifts is used to retrieve the correct cellIds from the grid ids
  std::vector<MInt> cellIdShifts(noTotalWallPoints + 1, 0);
  for(MInt globalWallId = 0; globalWallId < noTotalWallPoints; ++globalWallId) {
    Point<3> a(wallDistSurfCoordsAll[0][globalWallId], wallDistSurfCoordsAll[1][globalWallId], 0, globalWallId);
    pts.push_back(a);
    cellIdShifts[globalWallId + 1] = cellIdShifts[globalWallId] - 1 * (isStartEndPointAll[globalWallId] == -1);
  }
  // build up the tree
  KDtree<3> tree(pts);
  static constexpr MFloat radius = std::numeric_limits<MFloat>::max();
  static constexpr MFloat eps = 1e-16; // I hope this is a suitable value

  // go through all the cells an determine the closest distance
  MInt results[3];
  MFloat dists[3];
  for(MInt cellId = 0; cellId < m_noCells; ++cellId) {
    Point<3> pt(m_cells->coordinates[0][cellId], m_cells->coordinates[1][cellId], 0);
    MBool overflow = false;
    MInt nfound = tree.locatenearest(pt, radius, results, dists, 3, overflow);
    if(nfound != 3) mTerm(1, "Come on, your mesh should have at least 3 wall grid points!");
    const MInt globalWallId1 = results[0];
    MInt globalWallId1_ = results[1];
    MInt globalWallId1__ = results[2];
    nfound = 1;
    const MFloat r2 = POW2(wallDistSurfCoordsAll[0][globalWallId1] - wallDistSurfCoordsAll[0][globalWallId1_])
                      + POW2(wallDistSurfCoordsAll[1][globalWallId1] - wallDistSurfCoordsAll[1][globalWallId1_]);
    const MFloat r3 = POW2(wallDistSurfCoordsAll[0][globalWallId1] - wallDistSurfCoordsAll[0][globalWallId1__])
                      + POW2(wallDistSurfCoordsAll[1][globalWallId1] - wallDistSurfCoordsAll[1][globalWallId1__]);
    nfound += r2 < eps;
    nfound += r3 < eps;
    // It might happen that all three points have same distance to target point, but point 1 and 3 coincide and not
    // 1 and 2
    if(nfound == 2 && r3 < r2) {
      const MInt temp = globalWallId1_;
      globalWallId1_ = globalWallId1__;
      globalWallId1__ = temp;
    }

    ///////////////////////////////////////////////////////////////////////////
    /// Set distance, nearest neigbor cell identifiers & weighting
    ///////////////////////////////////////////////////////////////////////////
    // Case 1: wall grid point has only one neighbor, because it is the begin or end of a solid wall
    // Case 2: wall grid point has two neighbors, both on one rank (nfound==1)
    // Case 3: wall grid point has two neighbors, but one neighbor is on a different rank (nfound==2)
    MInt globalCellId1 = -1;
    MInt globalCellId2 = -1;
    MFloat distance{}, weight;
    if(nfound == 1 && abs(isStartEndPointAll[globalWallId1]) == 1) {
      // Case 1:
      const MInt globalWallId2 = globalWallId1 + isStartEndPointAll[globalWallId1];
      globalCellId1 = (isStartEndPointAll[globalWallId1] == -1) ? globalWallId1 - 1 : globalWallId1;
      globalCellId1 += cellIdShifts[globalWallId1];
      const MFloat p1[nDim] = {wallDistSurfCoordsAll[0][globalWallId1], wallDistSurfCoordsAll[1][globalWallId1]};
      const MFloat p2[nDim] = {wallDistSurfCoordsAll[0][globalWallId2], wallDistSurfCoordsAll[1][globalWallId2]};
      const MFloat pTarget[nDim] = {m_cells->coordinates[0][cellId], m_cells->coordinates[1][cellId]};
      MFloat fraction;
      shortestDistanceToLineElement(p1, p2, pTarget, distance, fraction);
      // Sanity check
      if(fraction > 0.5 + 1e-8) {
        stringstream msg;
        msg << "ERROR: p1=" << p1[0] << "|" << p1[1] << " p2=" << p2[0] << "|" << p2[1] << " pTarget=" << pTarget[0]
            << "|" << pTarget[1] << "  distance=" << distance << " fraction=" << fraction << endl;
        cout << msg.str();
        mTerm(1, "Fraction is >0.5, so nearest grid point is supposed to be p2 and not p1!");
      }
      weight = 1.0;
    } else { // grid point has neighbors to both sides (Case 1 & Case 2)
      // globalCellId1 lies in between the grid points globalWallId1 and globalWallId2
      // globalCellId2 lies in between the grid points globalWallId1 and globalWallId3 or
      //                       between the grid points globalWallId1_ and globalWallId3
      MInt globalWallId2;
      MInt globalWallId3;
      if(nfound == 1) {
        // Case 2:
        globalWallId2 = globalWallId1 - 1;
        globalWallId3 = globalWallId1 + 1;
        globalCellId1 = globalWallId2;
        globalCellId1 += cellIdShifts[globalWallId2];
        globalCellId2 = globalWallId1;
        globalCellId2 += cellIdShifts[globalWallId1];
      } else if(nfound == 2) {
        // Case 3:
        if(abs(isStartEndPointAll[globalWallId1]) + abs(isStartEndPointAll[globalWallId1_]) != 2) {
          cout << "ERROR: (x|y)_1=" << setprecision(12) << wallDistSurfCoordsAll[0][globalWallId1] << "|"
               << wallDistSurfCoordsAll[1][globalWallId1] << " (x|y)_2=" << wallDistSurfCoordsAll[0][globalWallId1_]
               << "|" << wallDistSurfCoordsAll[1][globalWallId1_] << endl;
          mTerm(1, "One grid point found twice. Both grid points have neighbors to both sides, which is unexpected!");
        }
        globalWallId2 = globalWallId1 + isStartEndPointAll[globalWallId1];
        globalWallId3 = globalWallId1_ + isStartEndPointAll[globalWallId1_];
        globalCellId1 = (isStartEndPointAll[globalWallId1] == -1) ? globalWallId1 - 1 : globalWallId1;
        globalCellId1 += cellIdShifts[globalWallId1];
        globalCellId2 = (isStartEndPointAll[globalWallId1_] == -1) ? globalWallId1_ - 1 : globalWallId1_;
        globalCellId2 += cellIdShifts[globalWallId1_];
      } else { // nfound==3
        stringstream msg;
        msg << "p1=" << wallDistSurfCoordsAll[0][globalWallId1] << "|" << wallDistSurfCoordsAll[1][globalWallId1]
            << ", p2=" << wallDistSurfCoordsAll[0][globalWallId1_] << "|" << wallDistSurfCoordsAll[1][globalWallId1_]
            << ", p3=" << wallDistSurfCoordsAll[0][globalWallId1__] << "|" << wallDistSurfCoordsAll[1][globalWallId1__]
            << endl;
        cout << msg.str();
        mTerm(1, "Three wall grid points coincide. This situation is unknown!");
      }

      // Now compute distances to both line elements
      MFloat distances[2];
      MFloat fractions[2];
      MFloat lineLengths[2];
      const MFloat p1[nDim] = {wallDistSurfCoordsAll[0][globalWallId1], wallDistSurfCoordsAll[1][globalWallId1]};
      const MFloat p2[nDim] = {wallDistSurfCoordsAll[0][globalWallId2], wallDistSurfCoordsAll[1][globalWallId2]};
      const MFloat p3[nDim] = {wallDistSurfCoordsAll[0][globalWallId3], wallDistSurfCoordsAll[1][globalWallId3]};
      const MFloat pTarget[nDim] = {m_cells->coordinates[0][cellId], m_cells->coordinates[1][cellId]};
      lineLengths[0] = shortestDistanceToLineElement(p1, p2, pTarget, distances[0], fractions[0]);
      lineLengths[1] = shortestDistanceToLineElement(p1, p3, pTarget, distances[1], fractions[1]);
      const MInt i = distances[1] < distances[0];
      // Sanity check
      if(fractions[i] > 0.5 + 1e-8) {
        stringstream msg;
        msg << "Error p1=" << p1[0] << "|" << p1[1] << " p2=" << p2[0] << "|" << p2[1] << " p3=" << p3[0] << "|"
            << p3[1] << " pTarget=" << pTarget[0] << "|" << pTarget[1] << "  distance=" << distance
            << " fraction=" << fractions[i] << " " << fractions[1 - i] << endl;
        cout << msg.str();
        mTerm(1, "Fraction >0.5, so nearest grid point can not be p1!");
      }
      // If wall is strongly concave, things might get nasty
      weight = (fractions[i] * lineLengths[i] + 0.5 * lineLengths[1 - i]) / (0.5 * (lineLengths[0] + lineLengths[1]));
      distance = distances[i];
      if(i == 1) {
        // Swap so that weight belongs to globalCellId1; weight of globalCellId2 is 1-weight
        const MInt temp = globalCellId1;
        globalCellId1 = globalCellId2;
        globalCellId2 = temp;
      }
    }

    if(distance < m_cells->fq[FQ->WALLDISTANCE][cellId]) {
      globalWallIds.insert(globalCellId1);
      globalWallIds.insert(globalCellId2);
      cellId2globalWallId[cellId] = std::make_pair(globalCellId1, globalCellId2);
      weighting[cellId] = weight;
      m_cells->fq[FQ->WALLDISTANCE][cellId] = distance;

      // DEBUG
      //      m_cells->fq[FQ->VAR1][cellId] = globalCellId1*weight + (1-weight)*globalCellId2;
    }
  }
  globalWallIds.erase(-1);

  // Convert displs based on grid points to displs based on cells
  for(MInt i = 1; i < (signed)displs.size(); ++i) {
    displs[i] = displs[i] + cellIdShifts[displs[i]];
  }
  auto getDomainId = [&displs, noRanks](const MInt id, const MInt hint = 0) {
    for(MInt rank = hint; rank < (signed)displs.size() - 1; ++rank) {
      if(displs[rank + 1] > id) return rank;
    }
    return noRanks - 1;
  };

  // globalWallId2localId: mapping of globalWallId to local wallId on current rank
  std::vector<MInt> globalWallId2localId((globalWallIds.empty()) ? 0 : *globalWallIds.rbegin() + 1);
  // recvWallIds: local wallIds from each rank, requested by current rank
  std::vector<MInt> recvWallIds(globalWallIds.size());
  std::vector<MInt> sendcounts(noRanks);
  MInt domainId_temp = 0;
  MInt localWallId = 0;
  for(auto it = globalWallIds.begin(); it != globalWallIds.end(); ++it, ++localWallId) {
    globalWallId2localId[*it] = localWallId;
    domainId_temp = getDomainId(*it, domainId_temp);
    // transform globalWallId "*it" to local wallId of rank domainId_temp
    recvWallIds[localWallId] = *it - displs[domainId_temp];
    ++sendcounts[domainId_temp];
  }

  // m_wallCellId2recvCell: mapping from local cellId to the position, where its closest wall surface information is
  // stored
  for(MInt cellId = 0; cellId < m_noCells; ++cellId) {
    if(cellId2globalWallId[cellId].first == -1 && cellId2globalWallId[cellId].second == -1) continue;
    // In case one globalWallId is -1, replace that id with the other one. Since weight of such a cell
    // is zero, it doesn't matter
    const MInt temp = (cellId2globalWallId[cellId].second == -1) ? cellId2globalWallId[cellId].first
                                                                 : cellId2globalWallId[cellId].second;
    auto temp2 = std::make_tuple(globalWallId2localId[cellId2globalWallId[cellId].first], globalWallId2localId[temp],
                                 weighting[cellId]);
    m_wallCellId2recvCell.insert({cellId, temp2});
  }

  // Broadcast send/receive infos
  std::vector<MInt> snghbrs(noRanks);
  std::iota(snghbrs.begin(), snghbrs.end(), 0);
  // m_wallSendcounts: # cells to be sent to each rank
  m_wallSendcounts.resize(noRanks);
  // m_wallSendCells: local wallIds to be sent to each rank
  m_wallSendCells = maia::mpi::mpiExchangePointToPoint(&recvWallIds[0],
                                                       &snghbrs[0],
                                                       noRanks,
                                                       &sendcounts[0],
                                                       &snghbrs[0],
                                                       noRanks,
                                                       m_solver->m_StructuredComm,
                                                       m_solver->domainId(),
                                                       1,
                                                       m_wallSendcounts.data());

  // --> m_wallCellId2recvCell, m_wallSendCells and m_wallSendcounts can be used later to get information
  //     from nearest wall neighbor, e.g., friction velocity at nearest wall cell

  // Sanity checks
  // 1)
  if((signed)m_wallSendCells.size() != std::accumulate(&m_wallSendcounts[0], &m_wallSendcounts[0] + noRanks, 0))
    mTerm(1, "Neeeeeeiiiiiiiiin!");

  // 2)
  if(cellmemSize != 0) {
    MInt min = std::numeric_limits<MInt>::max();
    MInt max = -1;
    for(auto id : m_wallSendCells) {
      // per wallDistInfo, we have one grid point more then grid cell
      if(id >= cellmemSize - noWallDistInfo) {
        cout << "DEBUG id=" << id << " cellmemSize=" << cellmemSize << endl;
        mTerm(1, "Something went wrong!");
      }
      min = mMin(id, min);
      max = mMax(id, max);
    }
    if(min != 0) {
      mTerm(1, "Scheisse: min!=0");
    }
    if(max != cellmemSize - 1 - noWallDistInfo) {
      cout << "SCHEISSE max=" << max << " cellmemSize=" << cellmemSize
           << " cellIdShifts[cellIdShifts.size()-2]=" << cellIdShifts[cellIdShifts.size() - 2] << endl;
      mTerm(1, "Scheisse ......");
    }
  }

  // Some debug output
#ifndef NDEBUG
  cout << "::computeLocalWallDistance: rank=" << m_solver->domainId() << " cellmemSize=" << cellmemSize
       << " noTotalWallPoints=" << noTotalWallPoints << " #cells near wall=" << m_wallCellId2recvCell.size()
       << " #globalWallIds=" << globalWallIds.size() << endl;
#endif
}


// Function returns shortest distance of point pTarget to line element consisting of the points
// p1 & p2 and the fraction at which the closest point on line element lies
template <MBool isRans>
MFloat StructuredBndryCnd2D<isRans>::shortestDistanceToLineElement(const MFloat (&p1)[nDim],
                                                                   const MFloat (&p2)[nDim],
                                                                   const MFloat (&pTarget)[nDim],
                                                                   MFloat& distance,
                                                                   MFloat& fraction) {
  TRACE();

  const MFloat linep1p2[nDim] = {p2[0] - p1[0], p2[1] - p1[1]};
  const MFloat linep1p2LengthPOW2 = POW2(linep1p2[0]) + POW2(linep1p2[1]);
  const MFloat linep1pT[nDim] = {pTarget[0] - p1[0], pTarget[1] - p1[1]};
  const MFloat dotProd = std::inner_product(std::begin(linep1p2), std::end(linep1p2), std::begin(linep1pT), 0.0);
  // If pProjection is required uncomment those lines
  // MFloat pProjection[nDim];
  if(dotProd < 0) {
    // std::copy(std::begin(p1), std::end(p1), std::begin(pProjection));
    distance = sqrt(std::inner_product(std::begin(linep1pT), std::end(linep1pT), std::begin(linep1pT), 0.0));
    fraction = 0.0;
  } else if(dotProd > linep1p2LengthPOW2) {
    // std::copy(std::begin(p2), std::end(p2), std::begin(pProjection));
    const MFloat linep2pT[nDim] = {pTarget[0] - p2[0], pTarget[1] - p2[1]};
    distance = sqrt(std::inner_product(std::begin(linep2pT), std::end(linep2pT), std::begin(linep2pT), 0.0));
    fraction = 1.0;
  } else {
    fraction = dotProd / linep1p2LengthPOW2;
    const MFloat pProjection[nDim] = {p1[0] + fraction * linep1p2[0], p1[1] + fraction * linep1p2[1]};
    const MFloat linepPpT[nDim] = {pTarget[0] - pProjection[0], pTarget[1] - pProjection[1]};
    distance = sqrt(std::inner_product(std::begin(linepPpT), std::end(linepPpT), std::begin(linepPpT), 0.0));
  }
  // fraction >0.5 would indicate that p2 is closer to pTarget then p1
  if(fraction < 0.0 || fraction > 1.0) mTerm(1, "fraction < 0.0 or fraction > 1.0");
  return linep1p2LengthPOW2;
}


///////////////////////////////////////////////////////////////////////////////////////////////////
// In the following the function computeLocalWallDistance() is split into the functions
//  - computeDistance2Map
//  - setUpNearMapComm
// This enables to process multiple groups of maps, e.g., maps of solid walls and maps of
// fluid-porous interfaces. I.e. we compute the distance to the nearest wall and the distance
// to the nearest fluid-porous interface. In between both functions, we can weight the distance
// to a solid wall differently from the distance to a fluid-porous interface and pick the smallest
// weighted distance from both lists (as is done in the function getCloserMap);
// computeLocalWallDistances is equal to:
/*
template <MBool isRans>
void StructuredBndryCnd2D<isRans>::computeLocalWallDistances(){
  std::vector<MInt> cellId2globalWallId;
  std::vector<MInt>  displs;
  computeDistance2Map(m_auxDataMap, m_cells->fq[FQ->WALLDISTANCE], cellId2globalWallId, displs);
  setUpNearMapComm(cellId2globalWallId, displs, m_wallCellId2recvCell, m_wallSendcounts, m_wallSendCells);
}
*/


// TODO_SS labels:FV,COMM Currently the wall distances in halo cells is not taken from the respective window cells
template <MBool isRans>
void StructuredBndryCnd2D<isRans>::computeLocalExtendedDistancesAndSetComm() {
  TRACE();

  // Store solid maps and fluid-porous maps in seperate lists
  vector<unique_ptr<StructuredWindowMap<nDim>>> wallDistInfo;
  vector<unique_ptr<StructuredWindowMap<nDim>>> FPDistInfo;
  for(auto& auxDataMap : m_auxDataMap) {
    if(auxDataMap->isFluidPorousInterface) {
      unique_ptr<StructuredWindowMap<nDim>> temp = make_unique<StructuredWindowMap<nDim>>();
      m_solver->m_windowInfo->mapCpy(auxDataMap, temp);
      FPDistInfo.push_back(move(temp));
    } else if((int)(auxDataMap->BC / 1000.0) == 1) {
      unique_ptr<StructuredWindowMap<nDim>> temp = make_unique<StructuredWindowMap<nDim>>();
      m_solver->m_windowInfo->mapCpy(auxDataMap, temp);
      wallDistInfo.push_back(move(auxDataMap));
    }
  }
#ifndef NDEBUG
  std::cout << "DomainId=" << m_solver->domainId() << ": " << FPDistInfo.size() << " FP maps & " << wallDistInfo.size()
            << " wall maps." << endl;
#endif

  //
  MInt noRanks;
  MPI_Comm_size(m_solver->m_StructuredComm, &noRanks);

  // Allocate temporary space; the distance to solid wall is temporarily saved in m_cells->fq
  ScratchSpace<MFloat> FP_distance(m_noCells, AT_, "FP_distance");

  // Solid
  m_wallSendcounts.assign(noRanks, 0);
  // r_cellId2globalWallId: mapping of cellId to the positions of the two cells which encloses the nearest point at the
  // wall
  std::vector<std::pair<MInt, MInt>> cellId2globalWallId; //(m_noCells, std::make_pair(-1, -1));
  std::vector<MInt> wallDispls;
  // weighting: weighting of the wall cell which is closest to a given cell; the weighting of its neighbor is
  // 1-weighting
  std::vector<MFloat> wallWeighting; //(m_noCells, 0.0);
  computeDistance2Map(wallDistInfo, m_cells->fq[FQ->WALLDISTANCE], cellId2globalWallId, wallDispls, wallWeighting);

  // Porous
  m_FPSendcounts.assign(noRanks, 0);
  // r_cellId2globalWallId: mapping of cellId to the positions of the two cells which encloses the nearest point at the
  // wall
  std::vector<std::pair<MInt, MInt>> cellId2globalFPId; //(m_noCells, std::make_pair(-1,-1));
  std::vector<MInt> FPDispls;
  // weighting: weighting of the wall cell which is closest to a given cell; the weighting of its neighbor is
  // 1-weighting
  std::vector<MFloat> FPWeighting; //(m_noCells, 0.0);
  computeDistance2Map(FPDistInfo, &FP_distance[0], cellId2globalFPId, FPDispls, FPWeighting);

  //////////////////////////
  // Before we reset we reset the communication information of porous solvers about fluid-porous interfaces,
  // we save that information
  std::vector<std::pair<MInt, MInt>> cellId2globalFPId_(cellId2globalFPId.begin(), cellId2globalFPId.end());
  if(m_solver->m_blockType != "porous") {
    for(MInt cellId = 0; cellId < m_noCells; ++cellId) {
      cellId2globalFPId_[cellId] = std::make_pair(-1, -1);
    }
  }
  setUpNearMapComm(cellId2globalFPId_, FPDispls, FPWeighting, m_FPCellId2recvCell_porous, m_FPSendcounts_porous,
                   m_FPSendCells_porous);
  /////////////////////////

  // Here we need to modify the walldistance; if we are inside porous solver set wall FPdistance to max value
  modifyFPDistance(FPDistInfo, &FP_distance[0], cellId2globalFPId, wallWeighting);

  //
  getCloserMap(m_cells->fq[FQ->WALLDISTANCE], cellId2globalWallId, &FP_distance[0], cellId2globalFPId,
               m_cells->fq[FQ->WALLDISTANCE]);

  // Solid
  setUpNearMapComm(cellId2globalWallId, wallDispls, wallWeighting, m_wallCellId2recvCell, m_wallSendcounts,
                   m_wallSendCells);

  // Porous
  setUpNearMapComm(cellId2globalFPId, FPDispls, FPWeighting, m_FPCellId2recvCell, m_FPSendcounts, m_FPSendCells);
}


/**
 * \brief Modifies the fluid-porous distance
 *
 * [Mößner, Michael, and Rolf Radespiel. "Modelling of turbulent flow over porous media using a
 * volume averaging approach and a Reynolds stress model." Computers & Fluids 108 (2015): 25-42.]
 * */
template <MBool isRans>
void StructuredBndryCnd2D<isRans>::modifyFPDistance(const vector<unique_ptr<StructuredWindowMap<nDim>>>& FPDistInfo,
                                                    MFloat* const FP_distance,
                                                    std::vector<std::pair<MInt, MInt>>& cellId2globalFPId,
                                                    const std::vector<MFloat>& wallWeighting) {
  if(!m_solver->m_porous) return;

  MFloat maxWallDistance = numeric_limits<MFloat>::max();
  maxWallDistance = Context::getSolverProperty<MFloat>("maxWallDistance", m_solverId, AT_, &maxWallDistance);

  // If this is a porous solver, then discard the distance to fluid-porous interface
  if(m_solver->m_blockType == "porous") {
    for(MInt cellId = 0; cellId < m_noCells; ++cellId) {
      FP_distance[cellId] = maxWallDistance;
      cellId2globalFPId[cellId] = std::make_pair(-1, -1);
    }
    //    return;
  }

  const MFloat* const RESTRICT por = &m_cells->fq[FQ->POROSITY][0];
  const MFloat* const RESTRICT Da = &m_cells->fq[FQ->DARCY][0];
  const auto& c_wd = m_solver->m_c_wd;

  const MInt noWallDistInfo = FPDistInfo.size(); // contains the size of the maps
  std::vector<MInt> wallNormalIndex(noWallDistInfo);
  std::vector<MInt> normalDir(noWallDistInfo);
  std::vector<std::array<MInt, 2 * nDim>> startEndTangential(noWallDistInfo);
  // Number of cells on current rank to store
  MInt cellmemSize = 0;
  for(MInt nn = 0; nn < noWallDistInfo; ++nn) {
    const MInt* const start = FPDistInfo[nn]->start1;
    const MInt* const end = FPDistInfo[nn]->end1;
    const MInt face = FPDistInfo[nn]->face;
    normalDir[nn] = face / 2;
    wallNormalIndex[nn] = start[normalDir[nn]] - (face % 2);

    MInt size = 1;
    for(MInt dim = 0; dim < nDim - 1; ++dim) {
      const MInt tangentialDir = (normalDir[nn] + (dim + 1)) % nDim;
      startEndTangential[nn][dim * 2 + 0] = start[tangentialDir]; // + m_noGhostLayers;
      startEndTangential[nn][dim * 2 + 1] = end[tangentialDir];   // - m_noGhostLayers;
      size *= startEndTangential[nn][dim * 2 + 1] - startEndTangential[nn][dim * 2 + 0];
    }
    cellmemSize += size;
#ifndef NDEBUG
    cout << "Debug output: FPDistInfo=" << nn << ": face=" << face << "; normalDir=" << normalDir[nn]
         << "; wallNormalIndex=" << wallNormalIndex[nn] << "; startEndTangential=" << startEndTangential[nn][0] << "|"
         << startEndTangential[nn][1] << "normal start/end=" << start[normalDir[nn]] << "|" << end[normalDir[nn]]
         << endl;
#endif
  }

  // 2.2) allocate the space in Scratch!
  // memory will not be needed later!
  MFloatScratchSpace correctionMem(std::max(1, cellmemSize), AT_, "correctionMem");

  // 4) computing the correction terms
  MInt cnt = 0;
  for(MInt nn = 0; nn < noWallDistInfo; ++nn) {
    MInt tIdx;
    MInt& i = (normalDir[nn] == 0) ? wallNormalIndex[nn] : tIdx;
    MInt& j = (normalDir[nn] == 0) ? tIdx : wallNormalIndex[nn];
    for(tIdx = startEndTangential[nn][0]; tIdx < startEndTangential[nn][1]; ++tIdx) {
      const MInt cellId = cellIndex(i, j);
      correctionMem[cnt++] = c_wd * sqrt(Da[cellId] / por[cellId]);
    }
  }

  // Broadcast all to all ranks
  MInt noRanks;
  MPI_Comm_size(m_solver->m_StructuredComm, &noRanks);
  // recvcounts: How many wall surface coordinates to receive from each domain
  std::vector<MInt> recvcounts(noRanks);
  MPI_Allgather(&cellmemSize, 1, MPI_INT, recvcounts.data(), 1, MPI_INT, m_solver->m_StructuredComm, AT_, "cellmemSize",
                "recvcounts");
  std::vector<MInt> displs(noRanks);
  for(MInt r = 0; r < noRanks - 1; ++r) {
    displs[r + 1] = displs[r] + recvcounts[r];
  }
  const MInt noTotalWallPoints = displs[noRanks - 1] + recvcounts[noRanks - 1];
  // wallDistSurfCoordsAll: contains the wall coordinates of all wall surfaces of all ranks
  std::vector<MFloat> correctionMemAll(noTotalWallPoints);
  MPI_Allgatherv(&correctionMem[0], cellmemSize, MPI_DOUBLE, correctionMemAll.data(), recvcounts.data(), displs.data(),
                 MPI_DOUBLE, m_solver->m_StructuredComm, AT_, "correctionMem", "correctionMemAll");

  // IMPORTANT: Here we assume for simplicity that by the correction the nearest fluid-porous point does not change;
  //            This is true if the material parameters along the fluid-porous interface are constant
  for(MInt cellId = 0; cellId < m_noCells; ++cellId) {
    const MInt globalFPId1 = cellId2globalFPId[cellId].first;
    // If globalFPId1==-1 then also globalFPId2==-1 and if globalFPId1!=-1 then also globalFPId2!=-1
    if(globalFPId1 > -1) {
      const MInt globalFPId2 = cellId2globalFPId[cellId].second;
      const MFloat my_FP_distance = FP_distance[cellId];
      const MFloat correction1 = correctionMemAll[globalFPId1];
      const MFloat correction2 = correctionMemAll[globalFPId2];
      const MFloat correction = wallWeighting[cellId] * correction1 + (1 - wallWeighting[cellId]) * correction2;
      const MFloat my_new_FP_distance = my_FP_distance + correction;
      if(my_new_FP_distance < maxWallDistance) {
        FP_distance[cellId] = my_new_FP_distance;
      } else {
        FP_distance[cellId] = maxWallDistance;
        cellId2globalFPId[cellId] = std::make_pair(-1, -1);
      }
    }
  }
}


/**
 * \brief Compute shortest distance to given set of maps
 *
 * Computes shortest distance to wall/fluid-porous interface, which is described by the grid as the sum of linear line
 * elements and stores the information required to exchange interpolated data from the nearest
 * neighbors to a cell.
 *
 * \date 2020-03-09
 *
 * \param[in] wallDistInfo e.g maps of all walls or fluid-porous-interfaces
 * \param[out] r_cells_fq_distance e.g. m_cells->fq[FQ->WALLDISTANCE]
 * \param[out] r_cellId2globalWallId
 * \parma[out] r_displs
 */
template <MBool isRans>
void StructuredBndryCnd2D<isRans>::computeDistance2Map(
    const vector<unique_ptr<StructuredWindowMap<nDim>>>& wallDistInfo,
    MFloat* const r_cells_fq_distance,
    std::vector<std::pair<MInt, MInt>>& r_cellId2globalWallId,
    std::vector<MInt>& r_displs,
    std::vector<MFloat>& r_weighting) {
  // Initialize array with high numbers
  MFloat maxWallDistance = numeric_limits<MFloat>::max();
  maxWallDistance = Context::getSolverProperty<MFloat>("maxWallDistance", m_solverId, AT_, &maxWallDistance);
  for(MInt id = 0; id < m_noCells; ++id) {
    r_cells_fq_distance[id] = maxWallDistance;
  }

  //  const vector<StructuredWindowMap*> wallDistInfo = m_auxDataMap;//m_solver->m_windowInfo->m_wallDistInfoMap;
  const MInt noWallDistInfo = wallDistInfo.size(); // contains the size of the maps
  std::vector<MInt> wallNormalIndex(noWallDistInfo);
  std::vector<MInt> normalDir(noWallDistInfo);
  std::vector<std::array<MInt, 2 * nDim>> startEndTangential(noWallDistInfo);
  // Number of grid points on current rank to store
  MInt cellmemSize = 0;
  for(MInt nn = 0; nn < noWallDistInfo; ++nn) {
    // start & end represent exactly the range of grid point indices, which form the map/boundary
    const MInt* const start = wallDistInfo[nn]->start1;
    const MInt* const end = wallDistInfo[nn]->end1;
    const MInt face = wallDistInfo[nn]->face;
    normalDir[nn] = face / 2;
    wallNormalIndex[nn] = start[normalDir[nn]];

    MInt size = 1;
    for(MInt dim = 0; dim < nDim - 1; ++dim) {
      const MInt tangentialDir = (normalDir[nn] + (dim + 1)) % nDim;
      startEndTangential[nn][dim * 2 + 0] = start[tangentialDir]; // + m_noGhostLayers;
      startEndTangential[nn][dim * 2 + 1] =
          end[tangentialDir] + 1; //+1, because grid coords not cell coords// - m_noGhostLayers;
      size *= startEndTangential[nn][dim * 2 + 1] - startEndTangential[nn][dim * 2 + 0];
    }
    cellmemSize += size;
#ifndef NDEBUG
    cout << "Debug output: wallDistInfo=" << nn << ": face=" << face << "; normalDir=" << normalDir[nn]
         << "; wallNormalIndex=" << wallNormalIndex[nn] << "; startEndTangential=" << startEndTangential[nn][0] << "|"
         << startEndTangential[nn][1] << "normal start/end=" << start[normalDir[nn]] << "|" << end[normalDir[nn]]
         << endl;
#endif
  }

  // 2.2) allocate the space for all the coordinates in Scratch!
  // memory will not be needed later!
  MFloatScratchSpace coordGridMem(std::max(1, nDim * cellmemSize), AT_, "coordGridMem");
  MFloatPointerScratchSpace wallDistSurfCoords(nDim, AT_, "wallDistSurfCoords");
  for(MInt dim = 0; dim < nDim; ++dim)
    wallDistSurfCoords[dim] = &coordGridMem[dim * cellmemSize];
  // Store if a grid point has a neighbor to the left or right on the same domain
  MIntScratchSpace isStartEndPoint(std::max(1, cellmemSize), AT_, "isStartEndPoint");

  // 4) Put wall grid points and information about existence of neighbor grid point in temporary buffers;
  MInt cnt = 0;
  for(MInt nn = 0; nn < noWallDistInfo; ++nn) {
    MInt tIdx;
    MInt& i = (normalDir[nn] == 0) ? wallNormalIndex[nn] : tIdx;
    MInt& j = (normalDir[nn] == 0) ? tIdx : wallNormalIndex[nn];
    for(tIdx = startEndTangential[nn][0]; tIdx < startEndTangential[nn][1]; ++tIdx) {
      const MInt p = getPointIdFromCell(i, j);
      wallDistSurfCoords[0][cnt] = m_grid->m_coordinates[0][p];
      wallDistSurfCoords[1][cnt] = m_grid->m_coordinates[1][p];
      // WRONG: The case no neighbor to the left and right can not occur
      if(tIdx == startEndTangential[nn][0] && tIdx == startEndTangential[nn][1] - 1)
        isStartEndPoint[cnt] = 3;
      else if(tIdx == startEndTangential[nn][0])
        isStartEndPoint[cnt] = 1; // no neighbor to the left on this rank
      else if(tIdx == startEndTangential[nn][1] - 1)
        isStartEndPoint[cnt] = -1; // no neighbor to the right on this rank
      else
        isStartEndPoint[cnt] =
            0; // TODO_SS labels:FV,GRID,totest is this necessary, or will it be default initialized anyway?
      ++cnt;
    }
  }

  // Broadcast all wall grid coordinates and isStartEndPoint to all ranks
  MInt noRanks;
  MPI_Comm_size(m_solver->m_StructuredComm, &noRanks);
  // recvcounts: How many wall grid coordinates to receive from each domain
  std::vector<MInt> recvcounts(noRanks);
  MPI_Allgather(&cellmemSize, 1, MPI_INT, recvcounts.data(), 1, MPI_INT, m_solver->m_StructuredComm, AT_, "cellmemSize",
                "recvcounts");
  r_displs.resize(noRanks);
  for(MInt r = 0; r < noRanks - 1; ++r) {
    r_displs[r + 1] = r_displs[r] + recvcounts[r];
  }
  const MInt noTotalWallPoints = r_displs[noRanks - 1] + recvcounts[noRanks - 1];
  // wallDistSurfCoordsAll: contains the wall grid coordinates of all ranks
  std::vector<std::vector<MFloat>> wallDistSurfCoordsAll(nDim);
  for(MInt dim = 0; dim < nDim; ++dim) {
    wallDistSurfCoordsAll[dim].resize(noTotalWallPoints);
    MPI_Allgatherv(&wallDistSurfCoords[dim][0], cellmemSize, MPI_DOUBLE, wallDistSurfCoordsAll[dim].data(),
                   recvcounts.data(), r_displs.data(), MPI_DOUBLE, m_solver->m_StructuredComm, AT_,
                   "wallDistSurfCoords", "wallDistSurfCoordsAll");
  }
  std::vector<MInt> isStartEndPointAll(noTotalWallPoints);
  MPI_Allgatherv(&isStartEndPoint[0], cellmemSize, MPI_INT, isStartEndPointAll.data(), recvcounts.data(),
                 r_displs.data(), MPI_INT, m_solver->m_StructuredComm, AT_, "isStartEndPoint", "isStartEndPointAll");
  // Also send own inputBoxId
  std::vector<MInt> blockIds(noRanks, -1);
  const MInt myBlockId = m_solver->m_blockId;
  MPI_Allgather(&myBlockId, 1, MPI_INT, blockIds.data(), 1, MPI_INT, m_solver->m_StructuredComm, AT_, "myBlockId",
                "blockId");

  // --> Now all ranks have all wall grid coordinates and the information if that grid point has any
  //     neighbor to the left or right on same domain

  // r_cellId2globalWallId: mapping of cellId to the positions of the two cells which encloses the nearest point at the
  // wall
  r_cellId2globalWallId.assign(m_noCells, std::make_pair(-1, -1));
  // weighting: weighting of the wall cell which is closest to a given cell; the weighting of its neighbor is
  // 1-weighting
  r_weighting.assign(m_noCells, 0.0);

  // Return early: the return is done here because we need to allocate the vector
  // r_cellId2globalWallId, before we can exit this function
  if(noTotalWallPoints == 0) return;

  // Build k-d-tree for a quick search
  // In the following Point<3> is used, but 3rd dimension is just dummy
  vector<Point<3>> pts;
  // cellIdShifts is used to retrieve the correct cellIds from the grid ids
  std::vector<MInt> cellIdShifts(noTotalWallPoints + 1, 0);
  for(MInt globalWallId = 0; globalWallId < noTotalWallPoints; ++globalWallId) {
    Point<3> a(wallDistSurfCoordsAll[0][globalWallId], wallDistSurfCoordsAll[1][globalWallId], 0, globalWallId);
    pts.push_back(a);
    cellIdShifts[globalWallId + 1] = cellIdShifts[globalWallId] - 1 * (isStartEndPointAll[globalWallId] == -1);
  }
  // build up the tree
  KDtree<3> tree(pts);

  auto getDomainId = [&r_displs, noRanks](const MInt id, const MInt hint = 0) {
    for(MInt rank = hint; rank < (signed)r_displs.size() - 1; ++rank) {
      if(r_displs[rank + 1] > id) return rank;
    }
    return noRanks - 1;
  };

  static constexpr MFloat radius = std::numeric_limits<MFloat>::max();
  static constexpr MFloat eps = 1e-16; // I hope this is a suitable value

  // go through all the cells an determine the closest distance
  MInt results[3];
  MFloat dists[3];
  for(MInt cellId = 0; cellId < m_noCells; ++cellId) {
    Point<3> pt(m_cells->coordinates[0][cellId], m_cells->coordinates[1][cellId], 0);
    MBool overflow = false;
    MInt nfound = tree.locatenearest(pt, radius, results, dists, 3, overflow);
    if(nfound != 3) mTerm(1, "Come on, your mesh should have at least 3 wall grid points!");
    MInt globalWallId1 = results[0];
    MInt globalWallId1_ = results[1];
    MInt globalWallId1__ = results[2];


    nfound = 1;
    const MFloat r2 = POW2(wallDistSurfCoordsAll[0][globalWallId1] - wallDistSurfCoordsAll[0][globalWallId1_])
                      + POW2(wallDistSurfCoordsAll[1][globalWallId1] - wallDistSurfCoordsAll[1][globalWallId1_]);
    const MFloat r3 = POW2(wallDistSurfCoordsAll[0][globalWallId1] - wallDistSurfCoordsAll[0][globalWallId1__])
                      + POW2(wallDistSurfCoordsAll[1][globalWallId1] - wallDistSurfCoordsAll[1][globalWallId1__]);
    nfound += r2 < eps;
    nfound += r3 < eps;
    // It might happen that all three points have same distance to target point, but point 1 and 3 coincide and not
    // 1 and 2
    if(nfound == 2 && r3 < r2) {
      const MInt temp = globalWallId1_;
      globalWallId1_ = globalWallId1__;
      globalWallId1__ = temp;
    }
    // Very special case: Because of the domain decomposition the start or end of a wall is on one domain without a
    // neighbor
    if(nfound == 2) {
      if(isStartEndPointAll[globalWallId1] == 3 && isStartEndPointAll[globalWallId1_] == 3) mTerm(1, "");
      if(isStartEndPointAll[globalWallId1] == 3) {
        globalWallId1 = globalWallId1_;
        nfound = 1;
      }
      if(isStartEndPointAll[globalWallId1_] == 3) {
        nfound = 1;
      }
    }

    ///////////////////////////////////////////////////////////////////////////
    /// Set distance, nearest neigbor cell identifiers & weighting
    ///////////////////////////////////////////////////////////////////////////
    // Case 1: wall grid point has only one neighbor, because it is the begin or end of a solid wall
    // Case 2: wall grid point has two neighbors, both on one rank (nfound==1)
    // Case 3: wall grid point has two neighbors, but one neighbor is on a different rank (nfound==2)
    MInt globalCellId1 = -1;
    MInt globalCellId2 = -1;
    MFloat distance = std::numeric_limits<MFloat>::max(), weight;
    if(nfound == 1 && abs(isStartEndPointAll[globalWallId1]) == 1) {
      // Case 1:
      const MInt globalWallId2 = globalWallId1 + isStartEndPointAll[globalWallId1];
      globalCellId1 = (isStartEndPointAll[globalWallId1] == -1) ? globalWallId1 - 1 : globalWallId1;
      globalCellId1 += cellIdShifts[globalWallId1];
      const MFloat p1[nDim] = {wallDistSurfCoordsAll[0][globalWallId1], wallDistSurfCoordsAll[1][globalWallId1]};
      const MFloat p2[nDim] = {wallDistSurfCoordsAll[0][globalWallId2], wallDistSurfCoordsAll[1][globalWallId2]};
      const MFloat pTarget[nDim] = {m_cells->coordinates[0][cellId], m_cells->coordinates[1][cellId]};
      MFloat fraction;
      shortestDistanceToLineElement(p1, p2, pTarget, distance, fraction);
      // Sanity check
      if(fraction > 0.5 + 1e-8) {
        stringstream msg;
        msg << "ERROR: p1=" << p1[0] << "|" << p1[1] << " p2=" << p2[0] << "|" << p2[1] << " pTarget=" << pTarget[0]
            << "|" << pTarget[1] << "  distance=" << distance << " fraction=" << fraction << endl;
        cout << msg.str();
        mTerm(1, "Fraction is >0.5, so nearest grid point is supposed to be p2 and not p1!");
      }
      weight = 1.0;
    } else { // grid point has neighbors to both sides (Case 1 & Case 2)
      // globalCellId1 lies in between the grid points globalWallId1 and globalWallId2
      // globalCellId2 lies in between the grid points globalWallId1 and globalWallId3 or
      //                       between the grid points globalWallId1_ and globalWallId3
      MInt globalWallId2;
      MInt globalWallId3;
      if(nfound == 1 /*&& isStartEndPointAll[globalWallId1]==0*/) {
        // Case 2:
        globalWallId2 = globalWallId1 - 1;
        globalWallId3 = globalWallId1 + 1;
        globalCellId1 = globalWallId2;
        globalCellId1 += cellIdShifts[globalWallId2];
        globalCellId2 = globalWallId1;
        globalCellId2 += cellIdShifts[globalWallId1];
      } else if(nfound == 2) { // globalWallId1 and globalWallId1_ coincide
        // Case 3:
        if(abs(isStartEndPointAll[globalWallId1]) + abs(isStartEndPointAll[globalWallId1_]) != 2) {
          // It can be a interface point which exists on two different solvers
          const MInt blockId1 = blockIds[getDomainId(globalWallId1)];
          const MInt blockId1_ = blockIds[getDomainId(globalWallId1_)];
          if(blockId1 == blockId1_) {
            cout << m_solver->domainId() << "|" << cellId << " " << globalWallId1 << " " << globalWallId1_ << " "
                 << globalWallId1__ << ": cellCoord=" << setprecision(12) << m_cells->coordinates[0][cellId] << "|"
                 << m_cells->coordinates[1][cellId] << "; wallCoord1=" << wallDistSurfCoordsAll[0][globalWallId1] << "|"
                 << wallDistSurfCoordsAll[1][globalWallId1]
                 << "; wallCoords1_=" << wallDistSurfCoordsAll[0][globalWallId1_] << "|"
                 << wallDistSurfCoordsAll[1][globalWallId1_]
                 << "; wallCoords1__=" << wallDistSurfCoordsAll[0][globalWallId1__] << "|"
                 << wallDistSurfCoordsAll[1][globalWallId1__] << endl;

            cout << "ERROR: (x|y)_1=" << setprecision(12) << wallDistSurfCoordsAll[0][globalWallId1] << "|"
                 << wallDistSurfCoordsAll[1][globalWallId1] << " (x|y)_2=" << wallDistSurfCoordsAll[0][globalWallId1_]
                 << "|" << wallDistSurfCoordsAll[1][globalWallId1_] << endl;
            mTerm(1, "One grid point found twice. Both grid points have neighbors to both sides, which is unexpected!");
          }
          if(blockId1 != myBlockId && blockId1_ != myBlockId) mTerm(1, "This case is not implemented yet!");
          if(blockId1_ == myBlockId) globalWallId1 = globalWallId1_;
          // Case 2:
          globalWallId2 = globalWallId1 - 1;
          globalWallId3 = globalWallId1 + 1;
          globalCellId1 = globalWallId2;
          globalCellId1 += cellIdShifts[globalWallId2];
          globalCellId2 = globalWallId1;
          globalCellId2 += cellIdShifts[globalWallId1];
        } else {
          globalWallId2 = globalWallId1 + isStartEndPointAll[globalWallId1];
          globalWallId3 = globalWallId1_ + isStartEndPointAll[globalWallId1_];
          globalCellId1 = (isStartEndPointAll[globalWallId1] == -1) ? globalWallId1 - 1 : globalWallId1;
          globalCellId1 += cellIdShifts[globalWallId1];
          globalCellId2 = (isStartEndPointAll[globalWallId1_] == -1) ? globalWallId1_ - 1 : globalWallId1_;
          globalCellId2 += cellIdShifts[globalWallId1_];
        }
      } else { // nfound==3
        // Here we could also check if the point belongs to different solvers
        stringstream msg;
        msg << setprecision(12) << "p1=" << wallDistSurfCoordsAll[0][globalWallId1] << "|"
            << wallDistSurfCoordsAll[1][globalWallId1] << ", p2=" << wallDistSurfCoordsAll[0][globalWallId1_] << "|"
            << wallDistSurfCoordsAll[1][globalWallId1_] << ", p3=" << wallDistSurfCoordsAll[0][globalWallId1__] << "|"
            << wallDistSurfCoordsAll[1][globalWallId1__] << endl;
        cout << msg.str();
        mTerm(1, "Three wall grid points coincide. This situation is unknown!");
      }

      // Now compute distances to both line elements
      MFloat distances[2];
      MFloat fractions[2];
      MFloat lineLengths[2];
      const MFloat p1[nDim] = {wallDistSurfCoordsAll[0][globalWallId1], wallDistSurfCoordsAll[1][globalWallId1]};
      const MFloat p2[nDim] = {wallDistSurfCoordsAll[0][globalWallId2], wallDistSurfCoordsAll[1][globalWallId2]};
      const MFloat p3[nDim] = {wallDistSurfCoordsAll[0][globalWallId3], wallDistSurfCoordsAll[1][globalWallId3]};
      const MFloat pTarget[nDim] = {m_cells->coordinates[0][cellId], m_cells->coordinates[1][cellId]};
      lineLengths[0] = shortestDistanceToLineElement(p1, p2, pTarget, distances[0], fractions[0]);
      lineLengths[1] = shortestDistanceToLineElement(p1, p3, pTarget, distances[1], fractions[1]);
      const MInt i = distances[1] < distances[0];
      // Sanity check
      if(fractions[i] > 0.5 + 1e-8) {
        stringstream msg;
        msg << "Error p1=" << p1[0] << "|" << p1[1] << " p2=" << p2[0] << "|" << p2[1] << " p3=" << p3[0] << "|"
            << p3[1] << " pTarget=" << pTarget[0] << "|" << pTarget[1] << "  distance=" << distance
            << " fraction=" << fractions[i] << " " << fractions[1 - i] << endl;
        cout << msg.str();
        mTerm(1, "Fraction >0.5, so nearest grid point can not be p1!");
      }
      // If wall is strongly concave, things might get nasty
      weight = (fractions[i] * lineLengths[i] + 0.5 * lineLengths[1 - i]) / (0.5 * (lineLengths[0] + lineLengths[1]));
      distance = distances[i];
      if(i == 1) {
        // Swap so that weight belongs to globalCellId1; weight of globalCellId2 is 1-weight
        const MInt temp = globalCellId1;
        globalCellId1 = globalCellId2;
        globalCellId2 = temp;
      }
    }

    if(distance < r_cells_fq_distance[cellId]) {
      ASSERT(globalCellId1 != -1, "");
      // In case one globalId is -1, replace that id with the other one. Since weight of such a cell
      // is zero, it doesn't matter
      globalCellId2 = (globalCellId2 == -1) ? globalCellId1 : globalCellId2;
      r_cellId2globalWallId[cellId] = std::make_pair(globalCellId1, globalCellId2);
      r_weighting[cellId] = weight;
      r_cells_fq_distance[cellId] = distance;

      // DEBUG
      //      m_cells->fq[FQ->VAR1][cellId] = globalCellId1*weight + (1-weight)*globalCellId2;
    }
  }

  // Convert displs based on grid points to displs based on cells
  for(MInt i = 1; i < (signed)r_displs.size(); ++i) {
    r_displs[i] = r_displs[i] + cellIdShifts[r_displs[i]];
  }

  // Some debug output
#ifndef NDEBUG
  MInt noNearWallCells = 0;
  // globalWallIds: store the positions of all wall surfaces in wallDistSurfCoordsAll, which are relavant
  //                for current rank
  std::set<MInt> globalWallIds;
  for(MInt cellId = 0; cellId < m_noCells; ++cellId) {
    if(r_cellId2globalWallId[cellId].first != -1) {
      ++noNearWallCells;
      globalWallIds.insert(r_cellId2globalWallId[cellId].first);
    }
    if(r_cellId2globalWallId[cellId].second != -1) {
      globalWallIds.insert(r_cellId2globalWallId[cellId].second);
    }
  }

  cout << "::computeDistance2Map: rank=" << m_solver->domainId() << " cellmemSize=" << cellmemSize
       << " noTotalWallPoints=" << noTotalWallPoints << " #cells near wall=" << noNearWallCells
       << " #globalWallIds=" << globalWallIds.size() << endl;
#endif
}

#if 0
/**
 * \brief Compute shortest distance to given set of maps
 *
 * \date 2020-03-09
 *
 * \param[in] wallDistInfo e.g maps of all walls or fluid-porous-interfaces
 * \param[out] r_cells_fq_distance e.g. m_cells->fq[FQ->WALLDISTANCE]
 * \param[out] r_cellId2globalWallId
 * \parma[out] r_displs
 */
template <MBool isRans>
void StructuredBndryCnd2D<isRans>::computeDistance2Map(const vector<unique_ptr<StructuredWindowMap<nDim>>>& wallDistInfo,
                                                       MFloat* const r_cells_fq_distance,
                                                       std::vector<MInt>& r_cellId2globalWallId,
                                                       std::vector<MInt>& r_displs){

  // Initialize array with high numbers
  MFloat maxWallDistance = numeric_limits<MFloat>::max();
  maxWallDistance = Context::getSolverProperty<MFloat>("maxWallDistance", m_solverId, AT_, &maxWallDistance);
  for(MInt id=0; id<m_noCells; ++id){
    r_cells_fq_distance[id] = maxWallDistance;
  }

//  const vector<StructuredWindowMap*> wallDistInfo = m_auxDataMap;//m_solver->m_windowInfo->m_wallDistInfoMap;
  const MInt noWallDistInfo = wallDistInfo.size(); //contains the size of the maps
  std::vector<MInt> wallNormalIndex(noWallDistInfo);
  std::vector<MInt> normalDir(noWallDistInfo);
  std::vector<std::array<MInt, 2*nDim>> startEndTangential(noWallDistInfo);
  // Number of cells on current rank to store
  MInt cellmemSize = 0;
  for (MInt nn=0; nn<noWallDistInfo;++nn) {
    const MInt* const start = wallDistInfo[nn]->start1;
    const MInt* const end = wallDistInfo[nn]->end1;
    const MInt face = wallDistInfo[nn]->face;
    normalDir[nn] = face/2;
    wallNormalIndex[nn] = start[normalDir[nn]];

    MInt size = 1;
    for (MInt dim = 0; dim < nDim-1; ++dim) {
      const MInt tangentialDir = (normalDir[nn]+(dim+1))%nDim;
      startEndTangential[nn][dim*2+0] = start[tangentialDir];// + m_noGhostLayers;
      startEndTangential[nn][dim*2+1] = end[tangentialDir];// - m_noGhostLayers;
      size *= startEndTangential[nn][dim*2+1] - startEndTangential[nn][dim*2+0];
    }
    cellmemSize += size;
//#ifndef NDEBUG
    cout << "Debug output: wallDistInfo=" << nn << ": face=" << face << "; normalDir=" 
         << normalDir[nn] << "; wallNormalIndex=" << wallNormalIndex[nn] << "; startEndTangential="
         << startEndTangential[nn][0] << "|" << startEndTangential[nn][1] << "normal start/end=" << start[normalDir[nn]] << "|" << end[normalDir[nn]] << endl;
//#endif
  }

  //2.2) allocate the space for all the coordinates in Scratch!
  //memory will not be needed later!
  MFloatScratchSpace coordCellMem(std::max(1,nDim*cellmemSize), AT_, "coordCellMem" );
  MFloatPointerScratchSpace wallDistSurfCoords(nDim, AT_, "wallDistSurfCoords");
  for (MInt dim = 0; dim < nDim; ++dim)
    wallDistSurfCoords[dim] = &coordCellMem[dim*cellmemSize];

  //4) computing the coordinates of surface center from corner points;
  MInt cnt = 0;
  for (MInt nn=0; nn<noWallDistInfo; ++nn) {
    MInt tIdx;
    MInt& i = (normalDir[nn]==0) ? wallNormalIndex[nn] : tIdx;
    MInt& j = (normalDir[nn]==0) ? tIdx : wallNormalIndex[nn];
    MInt inc[] = {0, 0};
    inc[1-normalDir[nn]] = 1;
    for (tIdx = startEndTangential[nn][0]; tIdx < startEndTangential[nn][1]; ++tIdx) {
      const MInt p1 = getPointIdFromCell(i, j);
      const MInt p2 = getPointIdFromPoint(p1, inc[0], inc[1]);
      wallDistSurfCoords[0][cnt] = 0.5*(m_grid->m_coordinates[0][p1] + m_grid->m_coordinates[0][p2]);
      wallDistSurfCoords[1][cnt] = 0.5*(m_grid->m_coordinates[1][p1] + m_grid->m_coordinates[1][p2]);
      ++cnt;
    }
  }

  // Broadcast all wall surface coordinates to all ranks
  // Note: actually this approach does not determine the absolute shortest
  // distance to a wall, but only the shortest distance among all wall surface centroids
  MInt noRanks;
  MPI_Comm_size(m_solver->m_StructuredComm, &noRanks);
  // recvcounts: How many wall surface coordinates to receive from each domain
  std::vector<MInt> recvcounts(noRanks);
  MPI_Allgather(&cellmemSize, 1, MPI_INT, recvcounts.data(), 1, MPI_INT, m_solver->m_StructuredComm, AT_,
                "cellmemSize", "recvcounts");
  r_displs.resize(noRanks);
  for (MInt r = 0; r < noRanks-1; ++r) {
    r_displs[r+1] = r_displs[r] + recvcounts[r];
  }
  const MInt noTotalWallPoints = r_displs[noRanks-1] + recvcounts[noRanks-1];
  // wallDistSurfCoordsAll: contains the wall coordinates of all wall surfaces of all ranks
  std::vector<std::vector<MFloat>> wallDistSurfCoordsAll(nDim);
  for (MInt dim = 0; dim < nDim; ++dim) {
    wallDistSurfCoordsAll[dim].resize(noTotalWallPoints);
    MPI_Allgatherv(&wallDistSurfCoords[dim][0], cellmemSize, MPI_DOUBLE,
        wallDistSurfCoordsAll[dim].data(), recvcounts.data(), r_displs.data(), MPI_DOUBLE,
        m_solver->m_StructuredComm, AT_, "wallDistSurfCoords", "wallDistSurfCoordsAll");
  }

  // --> Now all ranks have all wall surface coordinates

  // r_cellId2globalWallId: mapping of cellId to the position of its closest wall surface in wallDistSurfCoordsAll
  r_cellId2globalWallId.assign(m_noCells, -1);

  // Build k-d-tree for a quick search
  // In the following Point<3> is used, but 3rd dimension is just dummy
  vector < Point<3> > pts;
  for(MInt globalWallId=0; globalWallId<noTotalWallPoints; ++globalWallId){
    Point<3> a(wallDistSurfCoordsAll[0][globalWallId],wallDistSurfCoordsAll[1][globalWallId], 0, globalWallId);
    pts.push_back(a);
  }
  //build up the tree
  KDtree<3> tree(pts);
  MFloat distance = -1.0;
 
  //go through all the cells an determine the closest distance
  for (MInt cellId = 0; cellId < m_noCells; ++cellId) {
    distance=-1.1111111111111111; //to check
    Point<3> pt(m_cells->coordinates[0][cellId],m_cells->coordinates[1][cellId], 0);
    const MInt globalWallId = tree.nearest(pt, distance);

    if (distance < r_cells_fq_distance[cellId]) {
      r_cellId2globalWallId[cellId] = globalWallId;
      r_cells_fq_distance[cellId] = distance;
    }
  }

  // Some debug output
#ifndef NDEBUG
  MInt noNearWallCells = 0;
  // globalWallIds: store the positions of all wall surfaces in wallDistSurfCoordsAll, which are relavant
  //                for current rank
  std::set<MInt> globalWallIds;
  for (MInt cellId = 0; cellId < m_noCells; ++cellId) {
    if (r_cellId2globalWallId[cellId]!=-1) {
      ++noNearWallCells;
      globalWallIds.insert(r_cellId2globalWallId[cellId]);
    }
  }

  cout << "::computeDistance2Map: rank=" << m_solver->domainId() << " cellmemSize=" 
       << cellmemSize << " noTotalWallPoints=" << noTotalWallPoints << " #cells near wall=" 
       << noNearWallCells << " #globalWallIds=" << globalWallIds.size() << endl;
#endif
}
#endif


/**
 * \brief
 *
 * \date 2020-03-19
 *
 * \param[in] cells_fq_distance1
 * \param[in/out] cellId2globalWallId1
 * \param[in] cells_fq_distance2
 * \param[in/out] cellId2globalWallId2
 * \param[out] cells_fq_distance
 * \param[in] comparator
 */
template <MBool isRans>
template <typename T>
void StructuredBndryCnd2D<isRans>::getCloserMap(const MFloat* const cells_fq_distance1,
                                                std::vector<std::pair<MInt, MInt>>& cellId2globalWallId1,
                                                const MFloat* const cells_fq_distance2,
                                                std::vector<std::pair<MInt, MInt>>& cellId2globalWallId2,
                                                MFloat* const cells_fq_distance,
                                                T comparator) {
  // Loop over all cells of current domain and decide which one is closer
  for(MInt cellId = 0; cellId < m_noCells; ++cellId) {
    if(cellId2globalWallId1[cellId].first != -1 && cellId2globalWallId2[cellId].first != -1) {
      // In case both are possible candidates, apply comparator
      const MBool swtch = comparator(cellId, cells_fq_distance1[cellId], cells_fq_distance2[cellId]);
      if(swtch) {
        cellId2globalWallId2[cellId] = std::make_pair(-1, -1);
        cells_fq_distance[cellId] = cells_fq_distance1[cellId];
      } else {
        cellId2globalWallId1[cellId] = std::make_pair(-1, -1);
        cells_fq_distance[cellId] = cells_fq_distance2[cellId];
      }
    } else if(cellId2globalWallId1[cellId].first != -1) {
      cells_fq_distance[cellId] = cells_fq_distance1[cellId];
    } else {
      cells_fq_distance[cellId] = cells_fq_distance2[cellId];
    }
  }
}

/**
 * \brief
 *
 * \date 2020-03-19
 *
 * \param[in] cellId2globalWallId
 * \param[in] displs
 * \param[out] r_cellId2recvCell
 * \param[out] r_sendcounts
 * \param[out] r_sendCells
 */
template <MBool isRans>
void StructuredBndryCnd2D<isRans>::setUpNearMapComm(const std::vector<std::pair<MInt, MInt>>& cellId2globalWallId,
                                                    const std::vector<MInt>& displs,
                                                    const std::vector<MFloat>& weighting,
                                                    std::map<MInt, std::tuple<MInt, MInt, MFloat>>& r_cellId2recvCell,
                                                    std::vector<MInt>& r_sendcounts,
                                                    std::vector<MInt>& r_sendCells) {
  // Determine unqiue set of globalWallIds
  std::set<MInt> globalWallIds;
  MInt noNearMapCells = 0;
  for(MInt cellId = 0; cellId < m_noCells; ++cellId, ++noNearMapCells) {
    if(cellId2globalWallId[cellId].first != -1) {
      globalWallIds.insert(cellId2globalWallId[cellId].first);
      // if cellId2globalWallId[cellId].first!=-1, then also cellId2globalWallId[cellId].second!=-1
      globalWallIds.insert(cellId2globalWallId[cellId].second);
    }
  }

  const MInt noRanks = displs.size();
  auto getDomainId = [&displs, noRanks](const MInt id, const MInt hint = 0) {
    for(MInt rank = hint; rank < (signed)displs.size() - 1; ++rank) {
      if(displs[rank + 1] > id) return rank;
    }
    return noRanks - 1;
  };

  // globalWallId2localId: mapping of globalWallId in wallDistSurfCoordsAll to local wallId on current rank
  std::vector<MInt> globalWallId2localId((globalWallIds.empty()) ? 0 : *globalWallIds.rbegin() + 1);
  // recvWallIds: local wallIds from each rank, requested by current rank
  std::vector<MInt> recvWallIds(globalWallIds.size());
  std::vector<MInt> sendcounts(noRanks, 0);
  MInt domainId_temp = 0;
  MInt localWallId = 0;
  for(auto it = globalWallIds.begin(); it != globalWallIds.end(); ++it, ++localWallId) {
    globalWallId2localId[*it] = localWallId;
    domainId_temp = getDomainId(*it, domainId_temp);
    // transform globalWallId "*it" to local wallId of rank domainId_temp
    recvWallIds[localWallId] = *it - displs[domainId_temp];
    ++sendcounts[domainId_temp];
  }

  // r_cellId2recvCell: mapping from local cellId to the position, where its closest wall surface information is stored
  for(MInt cellId = 0; cellId < m_noCells; ++cellId) {
    // if cellId2globalWallId[cellId].first==-1, then also cellId2globalWallId[cellId].second==-1
    if(cellId2globalWallId[cellId].first == -1) continue;
    auto temp2 = std::make_tuple(globalWallId2localId[cellId2globalWallId[cellId].first],
                                 globalWallId2localId[cellId2globalWallId[cellId].second],
                                 weighting[cellId]);
    r_cellId2recvCell.insert({cellId, temp2});
  }

  // Broadcast send/receive infos
  std::vector<MInt> snghbrs(noRanks);
  std::iota(snghbrs.begin(), snghbrs.end(), 0);
  // r_sendcounts: # cells to be sent to each rank
  r_sendcounts.resize(noRanks);
  // r_sendCells: local wallIds to be sent to each rank
  r_sendCells = maia::mpi::mpiExchangePointToPoint(&recvWallIds[0],
                                                   &snghbrs[0],
                                                   noRanks,
                                                   &sendcounts[0],
                                                   &snghbrs[0],
                                                   noRanks,
                                                   m_solver->m_StructuredComm,
                                                   m_solver->domainId(),
                                                   1,
                                                   r_sendcounts.data());

#ifndef NDEBUG
  // Sanity checks
  // 1)
  if((signed)r_sendCells.size() != std::accumulate(&r_sendcounts[0], &r_sendcounts[0] + noRanks, 0))
    mTerm(1, "Neeeeeeiiiiiiiiin!");

  // 2) Check is not done for last rank
  const MInt cellmemSize =
      (m_solver->domainId() == noRanks - 1) ? 0 : displs[m_solver->domainId() + 1] - displs[m_solver->domainId()];
  if(cellmemSize != 0) {
    MInt min = std::numeric_limits<MInt>::max();
    MInt max = -1;
    for(auto id : r_sendCells) {
      if(id >= cellmemSize) mTerm(1, "Something went wrong!");
      min = mMin(id, min);
      max = mMax(id, max);
    }
    if(min != 0) {
      mTerm(1, "Scheisse: min!=0");
    }
    if(max != cellmemSize - 1) mTerm(1, "Scheisse......");
  }

  // Some debug output
  cout << "::computeLocalWallDistance: rank=" << m_solver->domainId()
       << " cellmemSize=" << cellmemSize /*<< " noTotalWallPoints=" << noTotalWallPoints*/ << " #cells near wall="
       << r_cellId2recvCell.size() << " #globalWallIds=" << globalWallIds.size() << endl;
#endif
}


template <MBool isRans>
void StructuredBndryCnd2D<isRans>::distributeWallAndFPProperties() {
  // Store solid maps and fluid-porous maps in seperate lists
  vector<unique_ptr<StructuredWindowMap<nDim>>> wallDistInfo;
  vector<MInt> wallMapIndex;
  vector<unique_ptr<StructuredWindowMap<nDim>>> FPDistInfo;
  vector<MInt> FPMapIndex;
  MInt cnt = 0;
  for(auto& auxDataMap : m_auxDataMap) {
    unique_ptr<StructuredWindowMap<nDim>> temp = make_unique<StructuredWindowMap<nDim>>();
    m_solver->m_windowInfo->mapCpy(auxDataMap, temp);
    if(auxDataMap->isFluidPorousInterface) {
      FPDistInfo.push_back(move(temp));
      FPMapIndex.push_back(cnt);
    } else if((int)(auxDataMap->BC / 1000.0) == 1) {
      wallDistInfo.push_back(move(temp));
      wallMapIndex.push_back(cnt);
    }
    ++cnt;
  }

  distributeMapProperties(wallDistInfo, m_wallSendCells, m_wallSendcounts, m_wallCellId2recvCell, wallMapIndex,
                          m_cells->fq[FQ->UTAU]);
  distributeMapProperties(FPDistInfo, m_FPSendCells, m_FPSendcounts, m_FPCellId2recvCell, FPMapIndex,
                          m_cells->fq[FQ->UTAU]);
  // Don't use FQ->UTAU2 -> use something which is only allocated in porous solver
  distributeMapProperties(FPDistInfo, m_FPSendCells_porous, m_FPSendcounts_porous, m_FPCellId2recvCell_porous,
                          FPMapIndex, m_cells->fq[FQ->UTAU2]);
}


/**
 * \brief
 *
 * \date 2020-03-??
 */
template <MBool isRans>
void StructuredBndryCnd2D<isRans>::distributeMapProperties(
    const std::vector<unique_ptr<StructuredWindowMap<nDim>>>& auxDataMap,
    const std::vector<MInt>& sendCells,
    const std::vector<MInt>& sendcounts,
    const std::map<MInt, std::tuple<MInt, MInt, MFloat>>& cellId2recvCell,
    const std::vector<MInt>& mapIndex,
    MFloat* const targetBuffer) {
  static const MFloat UT = m_solver->m_Ma * sqrt(PV->TInfinity);
  MFloat* const RESTRICT p = &m_cells->pvariables[PV->P][0];
  MFloat* const RESTRICT rho = &m_cells->pvariables[PV->RHO][0];

  // 1) Compute |cf| from cf-vector of al wall cells
  std::vector<MFloat> cf;
  MInt cnt = 0;
  for(MInt map = 0; map < (signed)auxDataMap.size(); ++map) {
    const MInt mapOffsetCf = m_cells->cfOffsets[mapIndex[map]];
    const MInt normalDir = auxDataMap[map]->face / 2;
    const MInt tangentialDir = 1 - normalDir;
    const MInt* const start = auxDataMap[map]->start1;
    const MInt* const end = auxDataMap[map]->end1;
    const MInt size = end[tangentialDir] - start[tangentialDir];
    ASSERT((signed)cf.size() == cnt, "cf.size() == cnt");
    cf.resize(cnt + size);
    MInt i = 0, j = 0;
    if((auxDataMap[map]->face) % 2 == 0) {
      i = start[normalDir]; // bottom wall
    } else {
      i = end[normalDir] - 1; // top wall //TODO_SS labels:FV,totest check later if the -1 is necessary
    }
    MInt& reali = (normalDir == 0) ? i : j;
    MInt& realj = (normalDir == 0) ? j : i;
    MInt ii = 0;
    for(j = start[tangentialDir]; j < end[tangentialDir]; ++j, ++cnt, ++ii) {
      // get id of boundary cell and cell above that one for extrapolation
      const MInt cellId = cellIndex(reali, realj);
      const MFloat T = m_solver->m_gamma * p[cellId] / rho[cellId];
      const MFloat nuLaminar = SUTHERLANDLAW(T) / rho[cellId];
      //    for (MInt i = 0; i < size; ++i, ++cnt) {
      // The following is no longer cf
      cf[cnt] =
          sqrt(0.5
               * sqrt(POW2(m_cells->cf[mapOffsetCf + 0 * size + ii]) + POW2(m_cells->cf[mapOffsetCf + 1 * size + ii])))
          * UT / nuLaminar;
    }
  }

  // 2) Fill sendBuffer
  ScratchSpace<MFloat> sendBuffer(std::max(1, (signed)sendCells.size()), FUN_, "sendBuffer");
  for(cnt = 0; cnt < (signed)sendCells.size(); ++cnt) {
    sendBuffer[cnt] = cf[sendCells[cnt]];
  }

  // 3) Send & receive relevant data
  MInt noRanks;
  MPI_Comm_size(m_solver->m_StructuredComm, &noRanks);
  std::vector<MInt> snghbrs(noRanks);
  std::iota(snghbrs.begin(), snghbrs.end(), 0);
  std::vector<MFloat> recvBuffer = maia::mpi::mpiExchangePointToPoint(&sendBuffer[0],
                                                                      &snghbrs[0],
                                                                      noRanks,
                                                                      sendcounts.data(),
                                                                      &snghbrs[0],
                                                                      noRanks,
                                                                      m_solver->m_StructuredComm,
                                                                      m_solver->domainId(),
                                                                      1);

  // 4) Scatter
  for(std::map<MInt, std::tuple<MInt, MInt, MFloat>>::const_iterator it = cellId2recvCell.begin();
      it != cellId2recvCell.end();
      ++it) {
    const MInt id1 = get<0>(it->second);
    const MInt id2 = get<1>(it->second);
    const MFloat weight = get<2>(it->second);
    const MFloat result = weight * recvBuffer[id1] + (1 - weight) * recvBuffer[id2];
    // TODO_SS labels:FV later maybe change this to yplus, also density correction misses
    targetBuffer[it->first] = result; // sqrt(0.5*result)*UT;//PV->UInfinity;
  }
}


// Computes cp or cf or both and saves result into the variable output
template <MBool isRans>
template <MBool calcCp, MBool calcCf, MBool interface>
void StructuredBndryCnd2D<isRans>::calc_cp_cf(const MInt cellId, const MInt cellIdP1, const MInt pIJ, const MInt pIJP,
                                              MFloat (&output)[calcCp + nDim * calcCf]) {
  // cf=nu*du/dy/(rho*u_8*u_8*0.5)
  // cp=2*(p-p_8)/(rho8*u_8*u_8*0.5)

  static const MFloat UT = m_solver->m_Ma * sqrt(PV->TInfinity);
  static const MFloat fstagnationPressure = F1 / (CV->rhoInfinity * F1B2 * POW2(UT));
  static const MFloat fre0 = F1 / (m_solver->m_Re0);

  // first get the position of the wall (reference!!)
  // x is always used as coordinate in normal direction
  const MFloat xRef = 0.5 * (m_grid->m_coordinates[0][pIJ] + m_grid->m_coordinates[0][pIJP]);
  const MFloat yRef = 0.5 * (m_grid->m_coordinates[1][pIJ] + m_grid->m_coordinates[1][pIJP]);

  // compute the pressure coefficient
  //-->
  // get the distcance
  const MFloat dx2 = sqrt(POW2(m_cells->coordinates[0][cellIdP1] - m_cells->coordinates[0][cellId])
                          + POW2(m_cells->coordinates[1][cellIdP1] - m_cells->coordinates[1][cellId]));
  const MFloat dx1 = sqrt(POW2(m_cells->coordinates[0][cellId] - xRef) + POW2(m_cells->coordinates[1][cellId] - yRef));

  if(calcCp) {
    // pressures
    const MFloat p1 = m_cells->pvariables[PV->P][cellId];
    const MFloat p2 = m_cells->pvariables[PV->P][cellIdP1];
    // extrapolation to the face
    const MFloat pW = ((p1 - p2) / dx2) * dx1 + p1;
    const MFloat cpn = (pW - PV->PInfinity) * fstagnationPressure;
    output[0] = cpn;
    //<--- finished computation of cp
  }

  if(calcCf) {
    // compute the skin-friction coefficient
    //-->
    // compute orthogonal distance surface-plane to cell center
    MFloat supportVec[nDim] = {};
    MFloat firstVec[nDim] = {};
    MFloat normalVec[nDim] = {};
    MFloat cellVec[nDim] = {};
    for(MInt dim = 0; dim < nDim; dim++) {
      supportVec[dim] = m_grid->m_coordinates[dim][pIJ];
      firstVec[dim] = m_grid->m_coordinates[dim][pIJP] - m_grid->m_coordinates[dim][pIJ];
      cellVec[dim] = m_cells->coordinates[dim][cellId];
    }

    normalVec[0] = -firstVec[1];
    normalVec[1] = firstVec[0];
    const MFloat normalLength = sqrt(POW2(normalVec[0]) + POW2(normalVec[1]));
    MFloat orthDist = F0;

    for(MInt dim = 0; dim < nDim; dim++) {
      orthDist += (cellVec[dim] - supportVec[dim]) * normalVec[dim];
    }

    orthDist = fabs(orthDist / normalLength);

    // extrapolate the Temperature to the wall
    const MFloat T1 = temperature(cellId);
    const MFloat T2 = temperature(cellIdP1);

    // extrapolation to the face
    const MFloat tBc = ((T1 - T2) / dx2) * dx1 + T1;
    const MFloat mue = SUTHERLANDLAW(tBc);

    // velocities
    const MFloat u1 = m_cells->pvariables[PV->U][cellId];
    const MFloat v1 = m_cells->pvariables[PV->V][cellId];

    // at 6000er-interfaces the velocity is not zero
    MFloat uBc = 0.0;
    MFloat vBc = 0.0;
    if(interface) {
      const MFloat u2 = m_cells->pvariables[PV->U][cellIdP1];
      const MFloat v2 = m_cells->pvariables[PV->V][cellIdP1];
      uBc = ((u1 - u2) / dx2) * dx1 + u1;
      vBc = ((v1 - v2) / dx2) * dx1 + v1;
    } else {
      if(m_solver->m_movingGrid) {
        uBc = F1B2 * (m_grid->m_velocity[0][pIJ] + m_grid->m_velocity[0][pIJP]);
        vBc = F1B2 * (m_grid->m_velocity[1][pIJ] + m_grid->m_velocity[1][pIJP]);
      }
    }

    // compute gradient
    MFloat dudeta = (u1 - uBc) / orthDist;
    MFloat dvdeta = (v1 - vBc) / orthDist;

    // using taylor expansion a second-order approximation
    // of du/dn can be achieved at the wall
    if(m_solver->m_forceSecondOrder) {
      const MFloat u2 = m_cells->pvariables[PV->U][cellIdP1];
      const MFloat v2 = m_cells->pvariables[PV->V][cellIdP1];
      const MFloat dx3 = dx1 + dx2;
      dudeta = (u1 * dx3 * dx3 + (dx1 * dx1 - dx3 * dx3) * uBc - u2 * dx1 * dx1) / (dx1 * dx3 * dx3 - dx1 * dx1 * dx3);
      dvdeta = (v1 * dx3 * dx3 + (dx1 * dx1 - dx3 * dx3) * vBc - v2 * dx1 * dx1) / (dx1 * dx3 * dx3 - dx1 * dx1 * dx3);
    }

    const MFloat taux = mue * dudeta * fre0;
    const MFloat tauy = mue * dvdeta * fre0;

    const MFloat cfx = taux * fstagnationPressure;
    const MFloat cfy = tauy * fstagnationPressure;


    output[calcCp + 0] = cfx, output[calcCp + 1] = cfy;
  }
}
template void StructuredBndryCnd2D<false>::calc_cp_cf<true, false, false>(const MInt, const MInt, const MInt,
                                                                          const MInt, MFloat (&)[1]);
template void StructuredBndryCnd2D<false>::calc_cp_cf<false, true, false>(const MInt, const MInt, const MInt,
                                                                          const MInt, MFloat (&)[2]);
template void StructuredBndryCnd2D<false>::calc_cp_cf<true, true, false>(const MInt, const MInt, const MInt, const MInt,
                                                                         MFloat (&)[3]);
template void StructuredBndryCnd2D<true>::calc_cp_cf<true, false, false>(const MInt, const MInt, const MInt, const MInt,
                                                                         MFloat (&)[1]);
template void StructuredBndryCnd2D<true>::calc_cp_cf<false, true, false>(const MInt, const MInt, const MInt, const MInt,
                                                                         MFloat (&)[2]);
template void StructuredBndryCnd2D<true>::calc_cp_cf<true, true, false>(const MInt, const MInt, const MInt, const MInt,
                                                                        MFloat (&)[3]);
template void StructuredBndryCnd2D<false>::calc_cp_cf<false, true, true>(const MInt, const MInt, const MInt, const MInt,
                                                                         MFloat (&)[2]);
template void StructuredBndryCnd2D<false>::calc_cp_cf<true, true, true>(const MInt, const MInt, const MInt, const MInt,
                                                                        MFloat (&)[3]);
template void StructuredBndryCnd2D<true>::calc_cp_cf<false, true, true>(const MInt, const MInt, const MInt, const MInt,
                                                                        MFloat (&)[2]);
template void StructuredBndryCnd2D<true>::calc_cp_cf<true, true, true>(const MInt, const MInt, const MInt, const MInt,
                                                                       MFloat (&)[3]);


// Note: no power computation and no moving grid
template <MBool isRans>
template <MBool calcCp, MBool calcCf, MBool calcIntegrals>
void StructuredBndryCnd2D<isRans>::computeFrictionPressureCoef_(const MBool auxDataWindows, const MBool computePower) {
  static_assert(calcCp || calcCf, "Something went wrong");

  // References for convenience
  const MInt noForceCoefs = m_solver->m_noForceDataFields;
  auto& m_forceCoef = m_solver->m_forceCoef;

  const MFloat UT = m_solver->m_Ma * sqrt(PV->TInfinity);
  const MFloat stagnationPressure = (CV->rhoInfinity * F1B2 * POW2(UT));
  const MFloat fstagnationEnergy = F1 / (CV->rhoInfinity * F1B2 * POW3(UT));

  if(calcIntegrals) {
    // Only those for which auxData is requested in the property file; not those which are, e.g., required
    // because of k-epsilon rans model
    const MInt noWalls = m_solver->m_windowInfo->m_auxDataWindowIds.size();

    // reset values to zero for the calculation
    for(MInt i = 0; i < (MInt)noForceCoefs * noWalls; ++i) {
      m_forceCoef[i] = F0;
    }
  }

  // loop over all maps
  for(MInt map = 0; map < (MInt)m_auxDataMap.size(); ++map) {
    //
    if(auxDataWindows) {
      MBool takeIt = false;
      for(auto& m : m_solver->m_windowInfo->m_auxDataWindowIds) {
        if(m_auxDataMap[map]->Id2 == m.second) {
          takeIt = true;
          break;
        }
      }
      if(!takeIt) continue;
    }

    const MBool interface = m_auxDataMap[map]->BC >= 6000 && m_auxDataMap[map]->BC < 6010;
    const auto calc_cp_cf_ = interface ? &StructuredBndryCnd2D<isRans>::calc_cp_cf<calcCp, calcCf, true>
                                       : &StructuredBndryCnd2D<isRans>::calc_cp_cf<calcCp, calcCf, false>;

    const MInt mapOffsetCf = m_cells->cfOffsets[map];
    const MInt mapOffsetCp = m_cells->cpOffsets[map];
    MInt* start = m_auxDataMap[map]->start1;
    MInt* end = m_auxDataMap[map]->end1;

    MInt n = 0;
    MInt normalDir = 0;
    // area
    MFloat area = F0;
    // skin-friction
    MFloat cf[nDim] = {};
    // pressure coefficient
    MFloat cp[nDim] = {};

    MFloat powerp[2] = {F0, F0};
    MFloat powerv[2] = {F0, F0};

    // indices
    MInt n10[nDim] = {};
    MInt n01[nDim] = {};
    MInt n1m1[nDim] = {};
    MInt pCoordDir[(2 * (nDim - 1) - 1) * nDim] = {};
    MInt i = 0, j = 0;

    // output for cp & cf-vector
    MFloat output[calcCp + nDim * calcCf];

    // determine the wall
    if((m_auxDataMap[map]->face) % 2 == 0) {
      n = -1; // bottom wall
      normalDir = m_auxDataMap[map]->face / 2;
      i = start[normalDir];
    } else {
      n = 1;                                         // top wall
      normalDir = (m_auxDataMap[map]->face - 1) / 2; // TODO_SS labels:FV Check if this and
      i = end[normalDir] - 1;                        // the -1 at this line is necessary
    }
    // determine the normals
    n1m1[normalDir] = -1 * n;
    n10[normalDir] = (MInt)(0.5 + (0.5 * (MFloat)n));
    n01[normalDir] = (MInt)(-0.5 + (0.5 * (MFloat)n));

    const MInt firstTangential = 1 - normalDir;
    // index shift for other surface point
    pCoordDir[firstTangential] = 1;
    // size in tangential direction
    const MInt sizeJ = end[firstTangential] - start[firstTangential];
    MInt& reali = (normalDir == 0) ? i : j;
    MInt& realj = (normalDir == 0) ? j : i;

    // loop over the tangential direction of the wall
    MInt jj = 0;
    for(j = start[firstTangential]; j < end[firstTangential]; j++) {
      // get id of boundary cell and cell above that one for extrapolation
      const MInt cellId = cellIndex(reali, realj);
      const MInt cellIdP1 = cellIndex(reali + n1m1[0], realj + n1m1[1]);

      // get point id and ids of three other surface corners
      const MInt pIJ = getPointIdFromCell(reali + n10[0], realj + n10[1]);
      const MInt pIJP = getPointIdFromPoint(pIJ, pCoordDir[0], pCoordDir[1]);

      // Actual calculation of cp & cf
      (this->*calc_cp_cf_)(cellId, cellIdP1, pIJ, pIJP, output);

      if(calcCp)
        // save pressure to map
        m_cells->cp[mapOffsetCp + jj] = output[0]; // cpn;

      if(calcCf) {
        // save skin-friction to map
        m_cells->cf[mapOffsetCf + 0 * sizeJ + jj] = output[calcCp + 0]; // cfx;
        m_cells->cf[mapOffsetCf + 1 * sizeJ + jj] = output[calcCp + 1]; // cfy;
      }

      if((calcCf && m_solver->m_detailAuxData) || calcIntegrals) {
        // first get the position of the wall (reference!!)
        // x is always used as coordinate in normal direction
        const MFloat xRef = 0.5 * (m_grid->m_coordinates[0][pIJ] + m_grid->m_coordinates[0][pIJP]);
        const MFloat yRef = 0.5 * (m_grid->m_coordinates[1][pIJ] + m_grid->m_coordinates[1][pIJP]);

        // this should be the valid method as in TFS
        // cp:
        // sum up all lengths and all cps with their surface width contribution
        const MFloat dxidx = m_cells->surfaceMetrics[nDim * normalDir + 0][cellIndex(reali + n01[0], realj + n01[1])];
        const MFloat dxidy = m_cells->surfaceMetrics[nDim * normalDir + 1][cellIndex(reali + n01[0], realj + n01[1])];

        if(calcCf && m_solver->m_detailAuxData) {
          // save surface contributions
          m_cells->cf[mapOffsetCf + 2 * sizeJ + jj] = dxidx;
          m_cells->cf[mapOffsetCf + 3 * sizeJ + jj] = dxidy;

          // save surface coordinates
          m_cells->cf[mapOffsetCf + 4 * sizeJ + jj] = xRef;
          m_cells->cf[mapOffsetCf + 5 * sizeJ + jj] = yRef;
        }

        if(calcIntegrals) {
          MFloat considerValue = F1;
          if(m_solver->m_auxDataCoordinateLimits) {
            if(m_solver->m_auxDataLimits[0] <= xRef && xRef <= m_solver->m_auxDataLimits[1]
               && m_solver->m_auxDataLimits[2] <= yRef && yRef <= m_solver->m_auxDataLimits[3]) {
              considerValue = F1;
            } else {
              considerValue = F0;
            }
          }

          if(calcCp) {
            cp[0] += (-1.0) * output[0] * dxidx * considerValue;
            cp[1] += (-1.0) * output[0] * dxidy * considerValue;
          }
          if(calcCf) {
            // cf
            cf[0] += output[calcCp + 0] * (sqrt(dxidx * dxidx + dxidy * dxidy)) * considerValue;
            cf[1] += output[calcCp + 1] * (sqrt(dxidx * dxidx + dxidy * dxidy)) * considerValue;
          }

          if(computePower) {
            const MFloat uWall =
                m_solver->m_movingGrid ? F1B2 * (m_grid->m_velocity[0][pIJ] + m_grid->m_velocity[0][pIJP]) : F0;
            const MFloat vWall =
                m_solver->m_movingGrid ? F1B2 * (m_grid->m_velocity[1][pIJ] + m_grid->m_velocity[1][pIJP]) : F0;

            const MFloat dp = output[0] * stagnationPressure;

            const MFloat P_px = (-1.0) * dp * uWall * fstagnationEnergy;
            const MFloat P_py = (-1.0) * dp * vWall * fstagnationEnergy;

            const MFloat taux = output[calcCp + 0] * stagnationPressure;
            const MFloat tauy = output[calcCp + 1] * stagnationPressure;

            const MFloat P_cfx = (taux)*uWall * fstagnationEnergy;
            const MFloat P_cfy = (tauy)*vWall * fstagnationEnergy;

            powerp[0] += P_px * dxidx * considerValue;
            powerp[1] += P_py * dxidy * considerValue;

            powerv[0] += P_cfx * (sqrt(dxidx * dxidx + dxidy * dxidy)) * considerValue;
            powerv[1] += P_cfy * (sqrt(dxidx * dxidx + dxidy * dxidy)) * considerValue;
          }
          // area
          area += sqrt(dxidx * dxidx + dxidy * dxidy) * considerValue;
        }
      }
      ++jj;
    }

    if(calcIntegrals) {
      // now add to lift and drag coefficients
      MInt count = 0;
      for(auto it = m_solver->m_windowInfo->m_auxDataWindowIds.cbegin();
          it != m_solver->m_windowInfo->m_auxDataWindowIds.cend();
          ++it) {
        if(m_auxDataMap[map]->Id2 == it->second) {
          m_forceCoef[count * noForceCoefs + 0] += cf[0];
          m_forceCoef[count * noForceCoefs + 1] += cf[1];
          m_forceCoef[count * noForceCoefs + 2] += 0.0; // the compressible part of the stress tensor - x
          m_forceCoef[count * noForceCoefs + 3] += 0.0; // the compressible part of the stress tensor - y
          m_forceCoef[count * noForceCoefs + 4] += cp[0];
          m_forceCoef[count * noForceCoefs + 5] += cp[1];
          m_forceCoef[count * noForceCoefs + 6] += area;
        }

        if(computePower) {
          m_forceCoef[count * noForceCoefs + 7] += powerv[0];
          m_forceCoef[count * noForceCoefs + 8] += powerv[1];
          m_forceCoef[count * noForceCoefs + 9] += powerp[0];
          m_forceCoef[count * noForceCoefs + 10] += powerp[1];
        }

        count++;
      }
    }
  }
}
template void StructuredBndryCnd2D<false>::computeFrictionPressureCoef_<false, true, false>(const MBool, const MBool);
template void StructuredBndryCnd2D<false>::computeFrictionPressureCoef_<true, true, true>(const MBool, const MBool);
template void StructuredBndryCnd2D<true>::computeFrictionPressureCoef_<false, true, false>(const MBool, const MBool);
template void StructuredBndryCnd2D<true>::computeFrictionPressureCoef_<true, true, true>(const MBool, const MBool);


// function to correct the index values in the map for the different boundary conditions
template <MBool isRans>
void StructuredBndryCnd2D<isRans>::correctBndryCndIndices() {
  TRACE();

  // in correcting cell Information
  for(MInt bcId = 0; bcId < m_noBndryCndIds; bcId++) {
    (this->*initBndryCndHandler[bcId])(bcId);
  }
}


template <MBool isRans>
inline MInt StructuredBndryCnd2D<isRans>::cellIndex(MInt i, MInt j) {
  return i + (j * m_nCells[1]);
}

template <MBool isRans>
inline MInt StructuredBndryCnd2D<isRans>::getPointIdFromCell(MInt i, MInt j) {
  return i + (j * (m_nCells[1] + 1));
}

template <MBool isRans>
inline MInt StructuredBndryCnd2D<isRans>::getPointIdFromPoint(MInt origin, MInt incI, MInt incJ) {
  return origin + incI + incJ * m_nPoints[1];
}

template <MBool isRans>
void StructuredBndryCnd2D<isRans>::correctWallDistanceAtBoundary(MInt bcId) {
  // Besides the wall distance the ghost cells eventually need the utau of the nearest wall
  // neighbors; Here we assume that the correct nearest neighbor is chosen automatically
  // and that hence we don't need to apply any correction; we also don't apply any extrapolation
  // to the utau as is done with the wall distance
  if(std::find(FQ->fqNames.begin(), FQ->fqNames.end(), "wallDistance") != FQ->fqNames.end()) {
    const MInt IJ[nDim] = {1, m_nCells[1]};
    MInt* start = m_physicalBCMap[bcId]->start1;
    MInt* end = m_physicalBCMap[bcId]->end1;

    // Here we find out the normal direction of the
    // boundary and the tangential direction.
    // This way we can make a general formulation of
    // the boundary condition
    const MInt face = m_physicalBCMap[bcId]->face;
    const MInt normalDir = face / 2;
    const MInt firstTangentialDir = (normalDir + 1) % nDim;
    const MInt normalDirStart = start[normalDir];
    const MInt firstTangentialStart = start[firstTangentialDir];
    const MInt firstTangentialEnd = end[firstTangentialDir];
    const MInt inc[nDim] = {IJ[normalDir], IJ[firstTangentialDir]};
    // determine indices for direction help
    const MInt n = (face % 2) * 2 - 1;                                //-1,+1
    const MInt g1 = normalDirStart + (MInt)(0.5 - (0.5 * (MFloat)n)); //+1,0
    const MInt g2 = normalDirStart + (MInt)(0.5 + (0.5 * (MFloat)n)); // 0,+1
    const MInt a1 = normalDirStart + (MInt)(0.5 - (1.5 * (MFloat)n)); //+2,-1
    const MInt a2 = normalDirStart + (MInt)(0.5 - (2.5 * (MFloat)n)); //+3,-2

    for(MInt t1 = firstTangentialStart; t1 < firstTangentialEnd; t1++) {
      const MInt cellIdG1 = g1 * inc[0] + t1 * inc[1];
      const MInt cellIdG2 = g2 * inc[0] + t1 * inc[1];
      const MInt cellIdA1 = a1 * inc[0] + t1 * inc[1];
      const MInt cellIdA2 = a2 * inc[0] + t1 * inc[1];

      m_cells->fq[FQ->WALLDISTANCE][cellIdG1] =
          F2 * m_cells->fq[FQ->WALLDISTANCE][cellIdA1] - m_cells->fq[FQ->WALLDISTANCE][cellIdA2];
      m_cells->fq[FQ->WALLDISTANCE][cellIdG2] =
          F2 * m_cells->fq[FQ->WALLDISTANCE][cellIdG1] - m_cells->fq[FQ->WALLDISTANCE][cellIdA1];
    }
  }
}

template <MBool isRans>
void StructuredBndryCnd2D<isRans>::initBc1000(MInt bcId) {
  (void)bcId;
}

template <MBool isRans>
void StructuredBndryCnd2D<isRans>::initBc1003(MInt bcId) {
  (void)bcId;
  m_isothermalWallTemperature = F1;
  if(Context::propertyExists("isothermalWallTemperature", m_solverId)) {
    m_isothermalWallTemperature =
        Context::getSolverProperty<MFloat>("isothermalWallTemperature", m_solverId, AT_, &m_isothermalWallTemperature);
  }
}

template <MBool isRans>
void StructuredBndryCnd2D<isRans>::initBc1004(MInt bcId) { // moving adiabatic wall
  (void)bcId;
}

template <MBool isRans>
void StructuredBndryCnd2D<isRans>::initBc2001(MInt bcId) {
  correctWallDistanceAtBoundary(bcId);
}

template <MBool isRans>
void StructuredBndryCnd2D<isRans>::initBc2002(MInt bcId) {
  (void)bcId;
}

template <MBool isRans>
void StructuredBndryCnd2D<isRans>::initBc2004(MInt bcId) {
  correctWallDistanceAtBoundary(bcId);
  bc2003(bcId);
}

template <MBool isRans>
void StructuredBndryCnd2D<isRans>::initBc2005(MInt bcId) {
  (void)bcId;
}

/* Initialize with standard pressure extrapolation
 * at inflow or prescribe p_inf at outflow
 *
 */
template <MBool isRans>
void StructuredBndryCnd2D<isRans>::initBc2006(MInt bcId) {
  bc2003(bcId);
}


template <MBool isRans>
void StructuredBndryCnd2D<isRans>::initBc2007(MInt bcId) {
  correctWallDistanceAtBoundary(bcId);
}


/**
 * \brief  Blasius bl inflow boundary condition
 *
 * Prescribes the velocity distribution for a Blasius
 * boundary layer at the inflow, with the Blasius helper functions
 *
 * \author Marian Albers
 */
template <MBool isRans>
void StructuredBndryCnd2D<isRans>::bc2999(MInt bcId) {
  const MInt IJ[nDim] = {1, m_nCells[1]};
  MInt* start = m_physicalBCMap[bcId]->start1;
  MInt* end = m_physicalBCMap[bcId]->end1;

  // Here we find out the normal direction of the
  // boundary and the two tangential directions.
  // This way we can make a general formulation of
  // the boundary condition
  const MInt face = m_physicalBCMap[bcId]->face;
  const MInt normalDir = face / 2;
  const MInt firstTangentialDir = (normalDir + 1) % nDim;
  const MInt normalDirStart = start[normalDir];
  const MInt firstTangentialStart = start[firstTangentialDir];
  const MInt firstTangentialEnd = end[firstTangentialDir];
  const MInt inc[nDim] = {IJ[normalDir], IJ[firstTangentialDir]};
  // determine indices for direction help
  const MInt n = (face % 2) * 2 - 1;                                //-1,+1
  const MInt g1 = normalDirStart + (MInt)(0.5 - (0.5 * (MFloat)n)); //+1,0
  const MInt g2 = normalDirStart + (MInt)(0.5 + (0.5 * (MFloat)n)); // 0,+1
  const MInt a1 = normalDirStart + (MInt)(0.5 - (1.5 * (MFloat)n)); //+2,-1
  // const MInt a2 = normalDirStart + (MInt)(0.5-(2.5*(MFloat)n)); //+3,-2

  for(MInt t1 = firstTangentialStart; t1 < firstTangentialEnd; t1++) {
    const MInt cellIdG1 = g1 * inc[0] + t1 * inc[1];
    const MInt cellIdG2 = g2 * inc[0] + t1 * inc[1];
    const MInt cellIdA1 = a1 * inc[0] + t1 * inc[1];
    // const MInt cellIdA2 = a2*inc[0] + t1*inc[1] + t2*inc[2];

    MFloat vel[nDim];
    m_solver->getBlasiusVelocity(cellIdG1, vel);
    m_cells->pvariables[PV->U][cellIdG1] = vel[0];
    m_cells->pvariables[PV->V][cellIdG1] = vel[1];
    m_cells->pvariables[PV->P][cellIdG1] = m_cells->pvariables[PV->P][cellIdA1];
    m_cells->pvariables[PV->RHO][cellIdG1] = CV->rhoInfinity;

    m_solver->getBlasiusVelocity(cellIdG2, vel);
    m_cells->pvariables[PV->U][cellIdG2] = vel[0];
    m_cells->pvariables[PV->V][cellIdG2] = vel[1];
    m_cells->pvariables[PV->P][cellIdG2] =
        F2 * m_cells->pvariables[PV->P][cellIdG1] - m_cells->pvariables[PV->P][cellIdA1];
    m_cells->pvariables[PV->RHO][cellIdG2] = CV->rhoInfinity;
  }
}


/**
 * \brief Channel flow BC
 */
template <MBool isRans>
void StructuredBndryCnd2D<isRans>::initBc2402(MInt bcId) {
  // for the channel flow we need the surface area of the entering and the outflow domain
  //==>determine Channelflow surface (assuming it does not change)
  MInt normal = 0;
  MInt Id = m_channelSurfaceIndexMap[bcId];
  if(Id < 0) cerr << "id smaller than zero ==> error" << endl;
  MInt* startface = &m_solver->m_windowInfo->channelSurfaceIndices[Id]->start1[0];
  MInt* endface = &m_solver->m_windowInfo->channelSurfaceIndices[Id]->end1[0];
  // find out which face it is
  MInt fixedind = -1;


  if(m_solver->m_initialCondition == 1233) {
    /////////////////////////////////////////////
    ///////// Turbulent channel flow ////////////
    /////////////////////////////////////////////
    MFloat uTau = m_solver->m_ReTau * m_solver->m_Ma * sqrt(PV->TInfinity) / m_solver->m_Re;
    m_solver->m_deltaP = -CV->rhoInfinity * POW2(uTau) * F2 * (m_solver->m_channelLength) / m_solver->m_channelHeight;
    m_solver->m_channelPresInlet = PV->PInfinity;
    m_solver->m_channelPresOutlet = PV->PInfinity + m_solver->m_deltaP;
    m_log << "=========== Turb. Channel Flow BC Summary =========== " << endl;
    m_log << "-->Turbulent channel flow deltaP: " << m_solver->m_deltaP << endl;
    m_log << "-->channel friciton velocity: " << uTau << endl;
    m_log << "-->Channel pressure inflow: " << m_solver->m_channelPresInlet << endl;
    m_log << "-->Channel pressure outflow: " << m_solver->m_channelPresOutlet << endl;
    m_log << "=========== Turb. Channel Flow BC Summary Finished =========== " << endl;
  } else if(m_solver->m_initialCondition == 1234) {
    cout.precision(10);
    /////////////////////////////////////////////
    ///////// Laminar channel flow //////////////
    /////////////////////////////////////////////

    m_solver->m_deltaP = -12.0 * PV->UInfinity * SUTHERLANDLAW(PV->TInfinity) * m_solver->m_channelLength
                         / (POW2(m_solver->m_channelHeight) * m_solver->m_Re0);
    m_log << "=========== Lam. Channel Flow BC Summary =========== " << endl;
    m_log << "Laminar channel deltaP: " << m_solver->m_deltaP << endl;
    m_log << "Theoretical cD total (both channel walls): "
          << m_solver->m_deltaP * m_solver->m_channelHeight * m_solver->m_channelWidth
                 / (0.5 * CV->rhoInfinity * POW2(PV->UInfinity))
          << endl;
    m_log << "Theoretical cD (single wall): "
          << m_solver->m_deltaP * m_solver->m_channelHeight * m_solver->m_channelWidth
                 / (0.5 * CV->rhoInfinity * POW2(PV->UInfinity) * 2.0)
          << endl;
    m_log << "Theoretical cF total (both channel walls): "
          << m_solver->m_deltaP * m_solver->m_channelHeight
                 / (0.5 * CV->rhoInfinity * POW2(PV->UInfinity) * m_solver->m_channelLength)
          << endl;
    m_log << "Theoretical cF (single wall): "
          << m_solver->m_deltaP * m_solver->m_channelHeight
                 / (0.5 * CV->rhoInfinity * POW2(PV->UInfinity) * m_solver->m_channelLength * 2.0)
          << endl;
    m_log << "=========== Lam. Channel Flow BC Summary Finished =========== " << endl;
  }

  for(MInt dim = 0; dim < nDim; dim++) {
    if(startface[dim] == endface[dim]) {
      // this is the normal
      fixedind = dim;
      if(startface[dim] == m_noGhostLayers) {
        normal = -1;
      } else {
        normal = 1;
      }
    }
  }

  if(m_physicalBCMap[bcId]->face == 0) {
    MPI_Comm_rank(m_solver->m_commChannelIn, &m_channelInflowRank);
  }

  switch(fixedind) {
    case 0: {
      MFloat surface = F0;
      MInt ii = 0;
      if(normal < 0) {
        ii = startface[0];

        MInt cellId = 0;
        for(MInt j = startface[1]; j < endface[1]; j++) {
          cellId = cellIndex(ii, j);
          surface += sqrt(POW2(m_cells->cellMetrics[0][cellId]) + POW2(m_cells->cellMetrics[1][cellId]));
        }
      } else {
        ii = endface[0];
        // activate this or the method below with the 4 points for the surface calculation
        // this is less exact (for straight surf.) but works better for curved surfaces

        MInt cellId = 0;
        for(MInt j = startface[1]; j < endface[1]; j++) {
          cellId = cellIndex(ii - 1, j);
          surface += sqrt(POW2(m_cells->cellMetrics[0][cellId]) + POW2(m_cells->cellMetrics[1][cellId]));
        }
      }

      if(normal < 0) {
        MPI_Allreduce(&surface, &m_channelSurfaceIn, 1, MPI_DOUBLE, MPI_SUM, m_solver->m_commChannelIn, AT_, "surface",
                      "m_channelSurfaceIn");
        cout << "ChannelInSurface: " << m_channelSurfaceIn << endl;
      } else {
        MPI_Allreduce(&surface, &m_channelSurfaceOut, 1, MPI_DOUBLE, MPI_SUM, m_solver->m_commChannelOut, AT_,
                      "surface", "m_channelSurfaceOut");
        cout << "ChannelOutSurface: " << m_channelSurfaceOut << endl;
      }

      break;
    }
    case 1: {
      mTerm(1, AT_, "surface calculation for j faces(channel not implemented)");
      break;
    }
    default: {
      mTerm(1, AT_, "surface calculation for given faces(channel not implemented)");
    }
  }
}


/**
 * \brief Rescaling inflow
 *
 *  Put values from field to gc to avoid nan, only matters
 *  for first Runge-Kutta Step at first timeStep
 */
template <MBool isRans>
void StructuredBndryCnd2D<isRans>::initBc2510(MInt bcId) {
  MInt* start = m_physicalBCMap[bcId]->start1;
  MInt* end = m_physicalBCMap[bcId]->end1;

  for(MInt var = 0; var < PV->noVariables; var++) {
    for(MInt i = start[0]; i < end[0]; i++) {
      for(MInt j = start[1]; j < end[1]; j++) {
        MInt cellId = cellIndex(i, j);
        MInt cellIdAdj = cellIndex(m_noGhostLayers, j);
        m_cells->pvariables[var][cellId] = m_cells->pvariables[var][cellIdAdj];
      }
    }
  }
}

/**
 * \brief Prescribe profile BC
 *
 *  In case of the initialStartup2600 flag is set
 *  (no matter if it was a restart or a initial startup)
 *  we load the values from the field to the GC and
 *  save them in the restart field
 */
template <MBool isRans>
void StructuredBndryCnd2D<isRans>::initBc2600(MInt bcId) {
  MInt* start = m_physicalBCMap[bcId]->start1;
  MInt* end = m_physicalBCMap[bcId]->end1;

  m_solver->m_bc2600 = true;

  switch(m_physicalBCMap[bcId]->face) {
    case 0: {
      if(m_solver->m_bc2600InitialStartup) {
        // First copy values from the field into the ghostcells
        for(MInt i = start[0]; i < end[0]; i++) {
          for(MInt j = start[1]; j < end[1]; j++) {
            const MInt cellId = cellIndex(i, j);
            const MInt cellIdAdj = cellIndex(m_noGhostLayers, j);
            for(MInt var = 0; var < PV->noVariables; var++) {
              m_cells->pvariables[var][cellId] = m_cells->pvariables[var][cellIdAdj];
            }
          }
        }
      }

      // Then copy the values from the ghostcells into restart field
      for(MInt i = start[0]; i < end[0]; i++) {
        for(MInt j = m_noGhostLayers; j < end[1] - m_noGhostLayers; j++) {
          const MInt cellId = cellIndex(i, j);
          const MInt cellIdBc = i + (j - m_noGhostLayers) * m_noGhostLayers;
          for(MInt var = 0; var < PV->noVariables; var++) {
            m_solver->m_bc2600Variables[var][cellIdBc] = m_cells->pvariables[var][cellId];
          }
        }
      }

      break;
    }
    default: {
      mTerm(1, AT_, "Face not implemented");
      break;
    }
  }
}

template <MBool isRans>
void StructuredBndryCnd2D<isRans>::initBc2021(MInt bcId) {
  MInt* start = m_physicalBCMap[bcId]->start1;
  MInt* end = m_physicalBCMap[bcId]->end1;
  for(MInt var = 0; var < PV->noVariables; var++) {
    for(MInt i = start[0]; i < end[0]; i++) {
      for(MInt j = start[1]; j < end[1]; j++) {
        MInt cellId = cellIndex(i, j);
        MInt cellIdAdj = cellIndex(m_noGhostLayers, j);
        m_cells->pvariables[var][cellId] = m_cells->pvariables[var][cellIdAdj];
      }
    }
  }
  m_bc2021Gradient = 1.0;
  m_bc2021Gradient = Context::getSolverProperty<MFloat>("bc2021Gradient", m_solverId, AT_, &m_bc2021Gradient);
}

template <MBool isRans>
void StructuredBndryCnd2D<isRans>::initBc3000(MInt bcId) {
  (void)bcId;
}


template <MBool isRans>
template <RansMethod ransMethod>
void StructuredBndryCnd2D<isRans>::bc1000_(MInt bcId) {
  TRACE();

  const MInt IJ[nDim] = {1, m_nCells[1]};
  MInt* start = m_physicalBCMap[bcId]->start1;
  MInt* end = m_physicalBCMap[bcId]->end1;

  // Here we find out the normal direction of the
  // boundary and the tangential direction.
  // This way we can make a general formulation of
  // the boundary condition
  const MInt face = m_physicalBCMap[bcId]->face;
  const MInt normalDir = face / 2;
  const MInt firstTangentialDir = (normalDir + 1) % nDim;
  const MInt normalDirStart = start[normalDir];
  const MInt firstTangentialStart = start[firstTangentialDir];
  const MInt firstTangentialEnd = end[firstTangentialDir];
  const MInt inc[nDim] = {IJ[normalDir], IJ[firstTangentialDir]};
  // determine indices for direction help
  const MInt n = (face % 2) * 2 - 1;                                //-1,+1
  const MInt g1 = normalDirStart + (MInt)(0.5 - (0.5 * (MFloat)n)); //+1,0
  const MInt g2 = normalDirStart + (MInt)(0.5 + (0.5 * (MFloat)n)); // 0,+1
  const MInt a1 = normalDirStart + (MInt)(0.5 - (1.5 * (MFloat)n)); //+2,-1
  const MInt a2 = normalDirStart + (MInt)(0.5 - (2.5 * (MFloat)n)); //+3,-2

  for(MInt t1 = firstTangentialStart; t1 < firstTangentialEnd; t1++) {
    const MInt cellIdG1 = g1 * inc[0] + t1 * inc[1];
    const MInt cellIdG2 = g2 * inc[0] + t1 * inc[1];
    const MInt cellIdA1 = a1 * inc[0] + t1 * inc[1];
    const MInt cellIdA2 = a2 * inc[0] + t1 * inc[1];

    const MFloat p1 = m_cells->pvariables[PV->P][cellIdA1];
    const MFloat rho1 = m_cells->pvariables[PV->RHO][cellIdA1];
    const MFloat u1 = m_cells->pvariables[PV->U][cellIdA1];
    const MFloat v1 = m_cells->pvariables[PV->V][cellIdA1];

    const MFloat p2 = m_cells->pvariables[PV->P][cellIdA2];
    const MFloat rho2 = m_cells->pvariables[PV->RHO][cellIdA2];
    const MFloat u2 = m_cells->pvariables[PV->U][cellIdA2];
    const MFloat v2 = m_cells->pvariables[PV->V][cellIdA2];

    m_cells->pvariables[PV->RHO][cellIdG1] = rho1;
    m_cells->pvariables[PV->U][cellIdG1] = -u1;
    m_cells->pvariables[PV->V][cellIdG1] = -v1;
    m_cells->pvariables[PV->P][cellIdG1] = p1;

    m_cells->pvariables[PV->RHO][cellIdG2] = rho2;
    m_cells->pvariables[PV->U][cellIdG2] = -u2;
    m_cells->pvariables[PV->V][cellIdG2] = -v2;
    m_cells->pvariables[PV->P][cellIdG2] = p2;

    if(ransMethod == RANS_SA || ransMethod == RANS_KEPSILON) {
      // Note: not all k-epsilon models have the same wall BC
      for(MInt ransVarId = 0; ransVarId < m_solver->m_noRansEquations; ++ransVarId) {
        // For e.g. SA-model the only rans variable is ransVar=nuTilde
        const MFloat ransVar1 = m_cells->pvariables[PV->RANS_VAR[ransVarId]][cellIdA1];
        const MFloat ransVar2 = m_cells->pvariables[PV->RANS_VAR[ransVarId]][cellIdA2];

        m_cells->pvariables[PV->RANS_VAR[ransVarId]][cellIdG1] = -ransVar1;
        m_cells->pvariables[PV->RANS_VAR[ransVarId]][cellIdG2] = -ransVar2;
      }
    }
  }
}

template <MBool isRans>
void StructuredBndryCnd2D<isRans>::bc1000(MInt bcId) {
  if(isRans) {
    const RansMethod ransMethod = m_solver->m_ransMethod;
    if(ransMethod == RANS_SA || ransMethod == RANS_SA_DV || ransMethod == RANS_FS) {
      bc1000_<RANS_SA>(bcId);
    } else if(ransMethod == RANS_KEPSILON) {
      bc1000_<RANS_KEPSILON>(bcId);
    }
  } else
    bc1000_<NORANS>(bcId);
}


// euler wall bc
template <MBool isRans>
void StructuredBndryCnd2D<isRans>::bc1001(MInt bcId) {
  TRACE();

  MInt* start = m_physicalBCMap[bcId]->start1;
  MInt* end = m_physicalBCMap[bcId]->end1;
  switch(m_physicalBCMap[bcId]->face) {
    case 2:
    case 3: {
      MInt cellShift = 0;
      if(m_physicalBCMap[bcId]->face == 2) {
        cellShift = 2 * m_noGhostLayers - 1;
      } else if(m_physicalBCMap[bcId]->face == 3) {
        cellShift = 2 * (m_nCells[0] - 1) - 2 * m_noGhostLayers + 1;
      }

      for(MInt j = start[1]; j < end[1]; j++) {
        for(MInt i = start[0]; i < end[0]; i++) {
          MInt cellIdA1 = -1, pIJ = -1;
          if(m_physicalBCMap[bcId]->face == 2) {
            pIJ = getPointIdFromCell(i, m_noGhostLayers);
            cellIdA1 = cellIndex(i, m_noGhostLayers);
          } else {
            pIJ = getPointIdFromCell(i, m_solver->m_nPoints[0] - 3);
            cellIdA1 = cellIndex(i, m_nCells[0] - 3);
          }
          const MFloat x1 = m_grid->m_coordinates[0][pIJ];
          const MFloat y1 = m_grid->m_coordinates[1][pIJ];
          const MFloat x2 = m_grid->m_coordinates[0][pIJ + 1];
          const MFloat y2 = m_grid->m_coordinates[1][pIJ + 1];

          const MFloat dydx = (y2 - y1) / (x2 - x1);
          const MFloat alpha = -atan(dydx);

          const MInt cellId = cellIndex(i, j);                // ghost
          const MInt cellIdadj = cellIndex(i, cellShift - j); // field

          const MFloat rho =
              m_cells->pvariables[PV->RHO][cellIdA1]; // apply rho from first active cell to both ghost-cells
          const MFloat u = m_cells->pvariables[PV->U][cellIdadj];
          const MFloat v = m_cells->pvariables[PV->V][cellIdadj];

          const MFloat uPrime = u * cos(alpha) - v * sin(alpha);
          const MFloat vPrime = -(u * sin(alpha) + v * cos(alpha));

          const MFloat uGC = uPrime * cos(-alpha) - vPrime * sin(-alpha);
          const MFloat vGC = uPrime * sin(-alpha) + vPrime * cos(-alpha);

          m_cells->pvariables[PV->RHO][cellId] = rho;
          m_cells->pvariables[PV->U][cellId] = uGC;
          m_cells->pvariables[PV->V][cellId] = vGC;

          // apply pressure from first active cell to all cells
          m_cells->pvariables[PV->P][cellId] = m_cells->pvariables[PV->P][cellIdA1];

          if(isRans) {
            m_cells->pvariables[PV->RANS_VAR[0]][cellId] = -m_cells->pvariables[PV->RANS_VAR[0]][cellIdadj];
          }
        }
      }
      break;
    }
    default: {
      mTerm(1, AT_, "Face direction not implemented)");
    }
  }
}

// isothermal Wall
template <MBool isRans>
void StructuredBndryCnd2D<isRans>::bc1003(MInt bcId) {
  const MInt IJ[nDim] = {1, m_nCells[1]};
  MInt* start = m_physicalBCMap[bcId]->start1;
  MInt* end = m_physicalBCMap[bcId]->end1;

  MFloat temp = m_isothermalWallTemperature * PV->TInfinity;
  const MFloat gamma = m_solver->m_gamma;

  // Here we find out the normal direction of the
  // boundary and the tangential direction.
  // This way we can make a general formulation of
  // the boundary condition
  const MInt face = m_physicalBCMap[bcId]->face;
  const MInt normalDir = face / 2;
  const MInt firstTangentialDir = (normalDir + 1) % nDim;
  const MInt normalDirStart = start[normalDir];
  const MInt firstTangentialStart = start[firstTangentialDir];
  const MInt firstTangentialEnd = end[firstTangentialDir];
  const MInt inc[nDim] = {IJ[normalDir], IJ[firstTangentialDir]};
  // determine indices for direction help
  const MInt n = (face % 2) * 2 - 1;                                //-1,+1
  const MInt g1 = normalDirStart + (MInt)(0.5 - (0.5 * (MFloat)n)); //+1,0
  const MInt g2 = normalDirStart + (MInt)(0.5 + (0.5 * (MFloat)n)); // 0,+1
  const MInt a1 = normalDirStart + (MInt)(0.5 - (1.5 * (MFloat)n)); //+2,-1
  const MInt a2 = normalDirStart + (MInt)(0.5 - (2.5 * (MFloat)n)); //+3,-2

  for(MInt t1 = firstTangentialStart; t1 < firstTangentialEnd; t1++) {
    const MInt cellIdG1 = g1 * inc[0] + t1 * inc[1];
    const MInt cellIdG2 = g2 * inc[0] + t1 * inc[1];
    const MInt cellIdA1 = a1 * inc[0] + t1 * inc[1];
    const MInt cellIdA2 = a2 * inc[0] + t1 * inc[1];

    const MFloat p1 = m_cells->pvariables[PV->P][cellIdA1];
    const MFloat rho1 = m_cells->pvariables[PV->RHO][cellIdA1];
    const MFloat u1 = m_cells->pvariables[PV->U][cellIdA1];
    const MFloat v1 = m_cells->pvariables[PV->V][cellIdA1];

    const MFloat p2 = m_cells->pvariables[PV->P][cellIdA2];
    const MFloat rho2 = m_cells->pvariables[PV->RHO][cellIdA2];
    const MFloat u2 = m_cells->pvariables[PV->U][cellIdA2];
    const MFloat v2 = m_cells->pvariables[PV->V][cellIdA2];

    const MFloat rhoWall = p1 * gamma / temp;
    const MFloat rhoG1 = F2 * rhoWall - rho1;
    const MFloat rhoG2 = F2 * rhoWall - rho2;

    m_cells->pvariables[PV->RHO][cellIdG1] = rhoG1;
    m_cells->pvariables[PV->U][cellIdG1] = -u1;
    m_cells->pvariables[PV->V][cellIdG1] = -v1;
    m_cells->pvariables[PV->P][cellIdG1] = p1;

    m_cells->pvariables[PV->RHO][cellIdG2] = rhoG2;
    m_cells->pvariables[PV->U][cellIdG2] = -u2;
    m_cells->pvariables[PV->V][cellIdG2] = -v2;
    m_cells->pvariables[PV->P][cellIdG2] = p2;

    if(isRans) {
      // Note: not all k-epsilon models have the same wall BC
      // For e.g. SA-model the only rans variable is ransVar=nuTilde
      const MFloat ransVar1 = m_cells->pvariables[PV->RANS_VAR[0]][cellIdA1];
      const MFloat ransVar2 = m_cells->pvariables[PV->RANS_VAR[0]][cellIdA2];

      m_cells->pvariables[PV->RANS_VAR[0]][cellIdG1] = -ransVar1;
      m_cells->pvariables[PV->RANS_VAR[0]][cellIdG2] = -ransVar2;
    }
  }
}


template <MBool isRans>
void StructuredBndryCnd2D<isRans>::bc1004(MInt bcId) {
  TRACE();

  MInt* start = m_physicalBCMap[bcId]->start1;
  MInt* end = m_physicalBCMap[bcId]->end1;

  const MInt IJ[nDim] = {1, m_nCells[1]};
  const MInt IJP[nDim] = {1, m_nPoints[1]};

  const MInt pp[2][4] = {{0, 0, 0, 1}, {0, 0, 1, 0}};

  // Here we find out the normal direction of the
  // boundary and the tangential direction.
  // This way we can make a general formulation of
  // the boundary condition
  const MInt face = m_physicalBCMap[bcId]->face;
  const MInt normalDir = face / 2;
  const MInt firstTangentialDir = (normalDir + 1) % nDim;
  const MInt normalDirStart = start[normalDir];
  const MInt firstTangentialStart = start[firstTangentialDir];
  const MInt firstTangentialEnd = end[firstTangentialDir];
  const MInt inc[nDim] = {IJ[normalDir], IJ[firstTangentialDir]};
  const MInt incp[nDim] = {IJP[normalDir], IJP[firstTangentialDir]};
  // determine indices for direction help
  const MInt n = (face % 2) * 2 - 1;                                //-1,+1
  const MInt g1 = normalDirStart + (MInt)(0.5 - (0.5 * (MFloat)n)); //+1,0
  const MInt g2 = normalDirStart + (MInt)(0.5 + (0.5 * (MFloat)n)); // 0,+1
  const MInt a1 = normalDirStart + (MInt)(0.5 - (1.5 * (MFloat)n)); //+2,-1
  const MInt a2 = normalDirStart + (MInt)(0.5 - (2.5 * (MFloat)n)); //+3,-2

  const MInt g1p = normalDirStart + 2 * ((MInt)(0.5 - (0.5 * (MFloat)n))); //+1,0

  for(MInt t1 = firstTangentialStart; t1 < firstTangentialEnd; t1++) {
    const MInt cellIdG1 = g1 * inc[0] + t1 * inc[1];
    const MInt cellIdG2 = g2 * inc[0] + t1 * inc[1];
    const MInt cellIdA1 = a1 * inc[0] + t1 * inc[1];
    const MInt cellIdA2 = a2 * inc[0] + t1 * inc[1];

    const MInt ij = g1p * incp[0] + t1 * incp[1];

    const MInt pp1 = getPointIdFromPoint(ij, pp[normalDir][0], pp[normalDir][1]);
    const MInt pp2 = getPointIdFromPoint(ij, pp[normalDir][2], pp[normalDir][3]);

    // compute the velocity of the surface centroid
    MFloat gridVel[2] = {F0, F0};
    // MFloat gridAcc[2] = {F0, F0};
    // MFloat firstVec[2] = {F0, F0};
    // MFloat secondVec[2] = {F0, F0};
    // MFloat normalVec[2] = {F0, F0};
    for(MInt dim = 0; dim < nDim; dim++) {
      // firstVec[dim] = m_grid->m_coordinates[dim][pp2] - m_grid->m_coordinates[dim][pp1];
      // secondVec[dim] = m_grid->m_coordinates[dim][pp3] - m_grid->m_coordinates[dim][pp1];
      gridVel[dim] = F1B2 * (m_grid->m_velocity[dim][pp1] + m_grid->m_velocity[dim][pp2]);
      // gridAcc[dim] = F1B2 * (m_grid->m_acceleration[dim][pp1] + m_grid->m_acceleration[dim][pp2]);
    }

    const MFloat p1 = m_cells->pvariables[PV->P][cellIdA1];
    const MFloat rho1 = m_cells->pvariables[PV->RHO][cellIdA1];
    const MFloat u1 = m_cells->pvariables[PV->U][cellIdA1];
    const MFloat v1 = m_cells->pvariables[PV->V][cellIdA1];

    const MFloat p2 = m_cells->pvariables[PV->P][cellIdA2];
    const MFloat rho2 = m_cells->pvariables[PV->RHO][cellIdA2];
    const MFloat u2 = m_cells->pvariables[PV->U][cellIdA2];
    const MFloat v2 = m_cells->pvariables[PV->V][cellIdA2];

    const MFloat pG1 = p1;
    const MFloat pG2 = p2;
    const MFloat rhoG1 = rho1;
    const MFloat rhoG2 = rho2;

    const MFloat uG1 = F2 * gridVel[0] - u1;
    const MFloat vG1 = F2 * gridVel[1] - v1;

    const MFloat uG2 = F2 * gridVel[0] - u2;
    const MFloat vG2 = F2 * gridVel[1] - v2;

    m_cells->pvariables[PV->RHO][cellIdG1] = rhoG1;
    m_cells->pvariables[PV->U][cellIdG1] = uG1;
    m_cells->pvariables[PV->V][cellIdG1] = vG1;
    m_cells->pvariables[PV->P][cellIdG1] = pG1;

    m_cells->pvariables[PV->RHO][cellIdG2] = rhoG2;
    m_cells->pvariables[PV->U][cellIdG2] = uG2;
    m_cells->pvariables[PV->V][cellIdG2] = vG2;
    m_cells->pvariables[PV->P][cellIdG2] = pG2;
  }
}

/**
 * \brief Subsonic Inflow <== tfs2001
 *
 *  rho=rho_inf, u=u_inf, v=v_inf, w=w_inf, dp/dn=0
 */
template <MBool isRans>
void StructuredBndryCnd2D<isRans>::bc2001(MInt bcId) {
  TRACE();

  // implemented for i-direction only for the moment
  MInt* start = m_physicalBCMap[bcId]->start1;
  MInt* end = m_physicalBCMap[bcId]->end1;
  switch(m_physicalBCMap[bcId]->face) {
    case 0: {
      MInt cellId = -1;
      MInt cellIdadj = -1;
      for(MInt j = start[1]; j < end[1]; j++) {
        for(MInt i = start[0]; i < end[0]; i++) {
          cellId = cellIndex(m_noGhostLayers - 1 - i, j);
          cellIdadj = cellIndex(m_noGhostLayers - i, j);

          m_cells->pvariables[PV->RHO][cellId] = CV->rhoInfinity;
          m_cells->pvariables[PV->U][cellId] = PV->UInfinity;
          m_cells->pvariables[PV->V][cellId] = PV->VInfinity;
          m_cells->pvariables[PV->P][cellId] = m_cells->pvariables[PV->P][cellIdadj];

          if(isRans) {
            for(MInt ransVarId = 0; ransVarId < m_solver->m_noRansEquations; ++ransVarId) {
              m_cells->pvariables[PV->RANS_VAR[ransVarId]][cellId] = PV->ransInfinity[ransVarId];
            }
          }
        }
      }
      break;
    }
    case 2: {
      MInt cellId = -1;
      MInt cellIdadj = -1;
      for(MInt j = start[1]; j < end[1]; j++) {
        for(MInt i = start[0]; i < end[0]; i++) {
          cellId = cellIndex(i, m_noGhostLayers - j - 1); // ghost
          cellIdadj = cellIndex(i, m_noGhostLayers - j);  // field

          m_cells->pvariables[PV->RHO][cellId] = CV->rhoInfinity;
          m_cells->pvariables[PV->U][cellId] = PV->UInfinity;
          m_cells->pvariables[PV->V][cellId] = PV->VInfinity;
          m_cells->pvariables[PV->P][cellId] = m_cells->pvariables[PV->P][cellIdadj];

          if(isRans) {
            for(MInt ransVarId = 0; ransVarId < m_solver->m_noRansEquations; ++ransVarId) {
              m_cells->pvariables[PV->RANS_VAR[ransVarId]][cellId] = PV->ransInfinity[ransVarId];
            }
          }
        }
      }
      break;
    }
    default: {
      mTerm(1, AT_, "Face direction not implemented)");
    }
  }
}


/**
 * \brief Subsonic Outflow, not really non-reflecting for face 0,1,3
 *
 *  Simplified characteristic approach (Whitfield)
 *  Still reflects pressure waves at the outflow
 *  rho=rho, u=u, v=v, w=w, p=p_inf
 *
 *  \author Pascal Meysonnat
 */
template <MBool isRans>
void StructuredBndryCnd2D<isRans>::bc2004(MInt bcId) {
  TRACE();

  const MInt IJ[nDim] = {1, m_nCells[1]};
  const MFloat gamma = m_solver->m_gamma;
  MInt* start = m_physicalBCMap[bcId]->start1;
  MInt* end = m_physicalBCMap[bcId]->end1;

  // Here we find out the normal direction of the
  // boundary and the tangential direction.
  // This way we can make a general formulation of
  // the boundary condition
  const MInt face = m_physicalBCMap[bcId]->face;
  const MInt normalDir = face / 2;
  const MInt firstTangentialDir = (normalDir + 1) % nDim;
  const MInt normalDirStart = start[normalDir];
  const MInt firstTangentialStart = start[firstTangentialDir];
  const MInt firstTangentialEnd = end[firstTangentialDir];
  const MInt inc[nDim] = {IJ[normalDir], IJ[firstTangentialDir]};
  // determine indices for direction help
  const MInt n = (face % 2) * 2 - 1;                                //-1,+1
  const MInt g1 = normalDirStart + (MInt)(0.5 - (0.5 * (MFloat)n)); //+1,0
  const MInt g2 = normalDirStart + (MInt)(0.5 + (0.5 * (MFloat)n)); // 0,+1
  const MInt a1 = normalDirStart + (MInt)(0.5 - (1.5 * (MFloat)n)); //+2,-1
  const MInt a2 = normalDirStart + (MInt)(0.5 - (2.5 * (MFloat)n)); //+3,-2

  for(MInt t1 = firstTangentialStart; t1 < firstTangentialEnd; t1++) {
    const MInt cellIdG1 = g1 * inc[0] + t1 * inc[1];
    const MInt cellIdG2 = g2 * inc[0] + t1 * inc[1];
    const MInt cellIdA1 = a1 * inc[0] + t1 * inc[1];
    const MInt cellIdA2 = a2 * inc[0] + t1 * inc[1];
    const MFloat dxidx = m_cells->surfaceMetrics[normalDir * nDim + 0][cellIdA1];
    const MFloat dxidy = m_cells->surfaceMetrics[normalDir * nDim + 1][cellIdA1];
    // multiply with n, so it will be -1 or +1 depending if we enter
    // or leave the domain of integration in positive direction
    const MFloat gradxi = n * F1 / sqrt(dxidx * dxidx + dxidy * dxidy);
    const MFloat dxHelp = dxidx * gradxi;
    const MFloat dyHelp = dxidy * gradxi;
    // speed of sound
    const MFloat cBC = sqrt(gamma * m_cells->pvariables[PV->P][cellIdG1] / m_cells->pvariables[PV->RHO][cellIdG1]);
    const MFloat rhoBC = m_cells->pvariables[PV->RHO][cellIdG1];
    const MFloat rhoInner = m_cells->pvariables[PV->RHO][cellIdA1];
    const MFloat uInner = m_cells->pvariables[PV->U][cellIdA1];
    const MFloat vInner = m_cells->pvariables[PV->V][cellIdA1];
    const MFloat pInner = m_cells->pvariables[PV->P][cellIdA1];
    const MFloat maContravariant = (dxidx * uInner + dxidy * vInner - m_cells->dxt[normalDir][cellIdA1]) * gradxi;
    if(maContravariant < F0) {
      // inflow
      const MFloat p = F1B2
                       * (pInner + PV->PInfinity
                          + rhoBC * cBC * (dxHelp * (uInner - PV->UInfinity) + dyHelp * (vInner - PV->VInfinity)));

      const MFloat rho = CV->rhoInfinity + (p - PV->PInfinity) / POW2(cBC);
      const MFloat help = (p - PV->PInfinity) / (rhoBC * cBC);

      m_cells->pvariables[PV->RHO][cellIdG1] = rho;
      m_cells->pvariables[PV->U][cellIdG1] = (PV->UInfinity + help * dxHelp);
      m_cells->pvariables[PV->V][cellIdG1] = (PV->VInfinity + help * dyHelp);
      m_cells->pvariables[PV->P][cellIdG1] = p;
      if(isRans) {
        for(MInt ransVarId = 0; ransVarId < m_solver->m_noRansEquations; ++ransVarId) {
          m_cells->pvariables[PV->RANS_VAR[ransVarId]][cellIdG1] = PV->ransInfinity[ransVarId];
        }
      }
    } else {
      // outflow
      const MFloat p = PV->PInfinity;
      const MFloat rho = rhoInner + (p - pInner) / POW2(cBC);
      const MFloat help = (p - pInner) / (rhoBC * cBC);

      m_cells->pvariables[PV->RHO][cellIdG1] = rho;
      m_cells->pvariables[PV->U][cellIdG1] = (uInner - help * dxHelp);
      m_cells->pvariables[PV->V][cellIdG1] = (vInner - help * dyHelp);
      m_cells->pvariables[PV->P][cellIdG1] = p;
      if(isRans) {
        for(MInt ransVarId = 0; ransVarId < m_solver->m_noRansEquations; ++ransVarId) {
          m_cells->pvariables[PV->RANS_VAR[ransVarId]][cellIdG1] =
              (F2 * m_cells->pvariables[PV->RANS_VAR[ransVarId]][cellIdA1]
               - m_cells->pvariables[PV->RANS_VAR[ransVarId]][cellIdA2]);
        }
      }
    }

    // extrapolate into second ghost cell
    for(MInt var = 0; var < PV->noVariables; var++) {
      m_cells->pvariables[var][cellIdG2] = F2 * m_cells->pvariables[var][cellIdG1] - m_cells->pvariables[var][cellIdA1];
    }
  }
}

/**
 * \brief Supersonic Inflow
 *
 *  rho=rho, u=u, v=v, w=w, p=p
 */
template <MBool isRans>
void StructuredBndryCnd2D<isRans>::bc2002(MInt bcId) {
  MInt* start = m_physicalBCMap[bcId]->start1;
  MInt* end = m_physicalBCMap[bcId]->end1;
  switch(m_physicalBCMap[bcId]->face) {
    case 1: {
      for(MInt j = start[1]; j < end[1]; j++) {
        for(MInt i = start[0]; i < end[0]; i++) {
          MInt cellId = cellIndex(i, j);
          m_cells->pvariables[PV->RHO][cellId] = CV->rhoInfinity;
          m_cells->pvariables[PV->U][cellId] = PV->UInfinity;
          m_cells->pvariables[PV->V][cellId] = PV->VInfinity;
          m_cells->pvariables[PV->P][cellId] = PV->PInfinity;
          if(isRans) {
            for(MInt ransVarId = 0; ransVarId < m_solver->m_noRansEquations; ++ransVarId) {
              m_cells->pvariables[PV->RANS_VAR[ransVarId]][cellId] = PV->ransInfinity[ransVarId];
            }
          }
        }
      }
      break;
    }
    default: {
      // do nothing
    }
  }
}

/**
 * \brief Subsonic in/outflow simple
 *
 * Simple extrapolation/prescription depending
 * on in/outflow condition
 *
 * \author Marian Albers
 */
template <MBool isRans>
void StructuredBndryCnd2D<isRans>::bc2003(MInt bcId) {
  const MInt IJ[nDim] = {1, m_nCells[1]};
  MInt* start = m_physicalBCMap[bcId]->start1;
  MInt* end = m_physicalBCMap[bcId]->end1;

  // Here we find out the normal direction of the
  // boundary and the tangential direction.
  // This way we can make a general formulation of
  // the boundary condition
  const MInt face = m_physicalBCMap[bcId]->face;
  const MInt normalDir = face / 2;
  const MInt firstTangentialDir = (normalDir + 1) % nDim;
  const MInt normalDirStart = start[normalDir];
  const MInt firstTangentialStart = start[firstTangentialDir];
  const MInt firstTangentialEnd = end[firstTangentialDir];
  const MInt inc[nDim] = {IJ[normalDir], IJ[firstTangentialDir]};
  // determine indices for direction help
  const MInt n = (face % 2) * 2 - 1;                                //-1,+1
  const MInt g1 = normalDirStart + (MInt)(0.5 - (0.5 * (MFloat)n)); //+1,0
  const MInt g2 = normalDirStart + (MInt)(0.5 + (0.5 * (MFloat)n)); // 0,+1
  const MInt a1 = normalDirStart + (MInt)(0.5 - (1.5 * (MFloat)n)); //+2,-1
  const MInt a2 = normalDirStart + (MInt)(0.5 - (2.5 * (MFloat)n)); //+3,-2

  for(MInt t1 = firstTangentialStart; t1 < firstTangentialEnd; t1++) {
    const MInt cellIdG1 = g1 * inc[0] + t1 * inc[1];
    const MInt cellIdG2 = g2 * inc[0] + t1 * inc[1];
    const MInt cellIdA1 = a1 * inc[0] + t1 * inc[1];
    const MInt cellIdA2 = a2 * inc[0] + t1 * inc[1];

    const MFloat dxidx = m_cells->surfaceMetrics[normalDir * nDim + 0][cellIdA1];
    const MFloat dxidy = m_cells->surfaceMetrics[normalDir * nDim + 1][cellIdA1];
    // leaving domain of integration in positive coordinate direction,
    // therefore multiply with positive F1
    const MFloat gradxi = n * F1 / sqrt(dxidx * dxidx + dxidy * dxidy);

    const MFloat rhoInner = m_cells->pvariables[PV->RHO][cellIdA1];
    const MFloat uInner = m_cells->pvariables[PV->U][cellIdA1];
    const MFloat vInner = m_cells->pvariables[PV->V][cellIdA1];
    const MFloat pInner = m_cells->pvariables[PV->P][cellIdA1];

    const MFloat maContravariant = (dxidx * uInner + dxidy * vInner) * gradxi;


    if(maContravariant < F0) {
      // inflow
      const MFloat p = pInner;
      const MFloat rho = CV->rhoInfinity;

      m_cells->pvariables[PV->RHO][cellIdG1] = rho;
      m_cells->pvariables[PV->U][cellIdG1] = PV->UInfinity;
      m_cells->pvariables[PV->V][cellIdG1] = PV->VInfinity;
      m_cells->pvariables[PV->P][cellIdG1] = p;

      if(isRans) {
        m_cells->pvariables[PV->RANS_VAR[0]][cellIdG1] = CV->ransInfinity[0] / CV->rhoInfinity;
      }
    } else {
      // outflow
      const MFloat p = PV->PInfinity;
      const MFloat rho = rhoInner;

      m_cells->pvariables[PV->RHO][cellIdG1] = rho;
      m_cells->pvariables[PV->U][cellIdG1] = uInner;
      m_cells->pvariables[PV->V][cellIdG1] = vInner;
      m_cells->pvariables[PV->P][cellIdG1] = p;
      if(isRans) {
        m_cells->pvariables[PV->RANS_VAR[0]][cellIdG1] =
            (F2 * m_cells->pvariables[PV->RANS_VAR[0]][cellIdA1] - m_cells->pvariables[PV->RANS_VAR[0]][cellIdA2]);
      }
    }

    // extrapolate into second ghost cell
    for(MInt var = 0; var < PV->noVariables; var++) {
      m_cells->pvariables[var][cellIdG2] = F2 * m_cells->pvariables[var][cellIdG1] - m_cells->pvariables[var][cellIdA1];
    }
  }
}

/**
 * \brief Supersonic Outflow,
 *
 *  rho=rho, u=u, v=v, w=w, p=p
 */
template <MBool isRans>
void StructuredBndryCnd2D<isRans>::bc2005(MInt bcId) {
  MInt* start = m_physicalBCMap[bcId]->start1;
  MInt* end = m_physicalBCMap[bcId]->end1;
  switch(m_physicalBCMap[bcId]->face) {
    case 1: {
      for(MInt j = start[1]; j < end[1]; j++) {
        for(MInt i = start[0]; i < end[0]; i++) {
          MInt cellId = cellIndex(i, j);
          MInt cellIdadj = cellIndex(i - 1, j);

          for(MInt var = 0; var < PV->noVariables; var++) {
            m_cells->pvariables[var][cellId] = m_cells->pvariables[var][cellIdadj];
          }
        }
      }
      break;
    }
    case 3: {
      for(MInt j = start[1]; j < end[1]; j++) {
        for(MInt i = start[0]; i < end[0]; i++) {
          MInt cellId = cellIndex(i, j);
          MInt cellIdadj = cellIndex(i, j - 1);
          for(MInt var = 0; var < PV->noVariables; var++) {
            m_cells->pvariables[var][cellId] = m_cells->pvariables[var][cellIdadj];
          }
        }
      }
      break;
    }
    default: {
      mTerm(1, AT_, "Face direction not implemented");
    }
  }
}

/**
 * \brief Characteristic in/outflow boundary for zero velocities
 *
 *  rho=rho, u=0, v=0, p=p_inf
 */
template <MBool isRans>
void StructuredBndryCnd2D<isRans>::bc2006(MInt bcId) {
  const MInt IJ[nDim] = {1, m_nCells[1]};
  const MFloat gamma = m_solver->m_gamma;
  MInt* start = m_physicalBCMap[bcId]->start1;
  MInt* end = m_physicalBCMap[bcId]->end1;

  // Here we find out the normal direction of the
  // boundary and the two tangential directions.
  // This way we can make a general formulation of
  // the boundary condition
  const MInt face = m_physicalBCMap[bcId]->face;
  const MInt normalDir = face / 2;
  const MInt firstTangentialDir = (normalDir + 1) % nDim;
  const MInt normalDirStart = start[normalDir];
  const MInt firstTangentialStart = start[firstTangentialDir];
  const MInt firstTangentialEnd = end[firstTangentialDir];
  const MInt inc[nDim] = {IJ[normalDir], IJ[firstTangentialDir]};

  const MInt n = (face % 2) * 2 - 1;                                //-1,+1
  const MInt g1 = normalDirStart + (MInt)(0.5 - (0.5 * (MFloat)n)); //+1,0
  const MInt g2 = normalDirStart + (MInt)(0.5 + (0.5 * (MFloat)n)); // 0,+1
  const MInt a1 = normalDirStart + (MInt)(0.5 - (1.5 * (MFloat)n)); //+2,-1
  const MInt a2 = normalDirStart + (MInt)(0.5 - (2.5 * (MFloat)n)); //+3,-2

  for(MInt t1 = firstTangentialStart; t1 < firstTangentialEnd; t1++) {
    const MInt cellIdG1 = g1 * inc[0] + t1 * inc[1];
    const MInt cellIdG2 = g2 * inc[0] + t1 * inc[1];
    const MInt cellIdA1 = a1 * inc[0] + t1 * inc[1];
    const MInt cellIdA2 = a2 * inc[0] + t1 * inc[1];
    const MFloat dxidx = m_cells->surfaceMetrics[normalDir * nDim + 0][cellIdA1];
    const MFloat dxidy = m_cells->surfaceMetrics[normalDir * nDim + 1][cellIdA1];
    // multiply with n, so it will be -1 or +1 depending if we enter
    // or leave the domain of integration in positive direction
    const MFloat gradxi = n * F1 / sqrt(dxidx * dxidx + dxidy * dxidy);
    const MFloat dxHelp = dxidx * gradxi;
    const MFloat dyHelp = dxidy * gradxi;

    const MFloat cBC = sqrt(gamma * m_cells->pvariables[PV->P][cellIdG1] / m_cells->pvariables[PV->RHO][cellIdG1]);
    const MFloat rhoBC = m_cells->pvariables[PV->RHO][cellIdG1];
    const MFloat rhoInner = m_cells->pvariables[PV->RHO][cellIdA1];
    const MFloat uInner = m_cells->pvariables[PV->U][cellIdA1];
    const MFloat vInner = m_cells->pvariables[PV->V][cellIdA1];
    const MFloat pInner = m_cells->pvariables[PV->P][cellIdA1];

    const MFloat maContravariant = (dxidx * uInner + dxidy * vInner - m_cells->dxt[normalDir][cellIdA1]) * gradxi;

    if(maContravariant < F0) {
      // inflow
      const MFloat p =
          F1B2 * (pInner + PV->PInfinity + rhoBC * cBC * (dxHelp * (uInner - 0.0) + dyHelp * (vInner - 0.0)));

      const MFloat rho = CV->rhoInfinity + (p - PV->PInfinity) / POW2(cBC);
      const MFloat help = (p - PV->PInfinity) / (rhoBC * cBC);

      m_cells->pvariables[PV->RHO][cellIdG1] = rho;
      m_cells->pvariables[PV->U][cellIdG1] = (0.0 + help * dxHelp);
      m_cells->pvariables[PV->V][cellIdG1] = (0.0 + help * dyHelp);
      m_cells->pvariables[PV->P][cellIdG1] = p;
      if(isRans) {
        m_cells->pvariables[PV->RANS_VAR[0]][cellIdG1] = PV->ransInfinity[0];
      }
    } else {
      // outflow
      const MFloat p = PV->PInfinity;
      const MFloat rho = rhoInner + (p - pInner) / POW2(cBC);
      const MFloat help = (p - pInner) / (rhoBC * cBC);

      m_cells->pvariables[PV->RHO][cellIdG1] = rho;
      m_cells->pvariables[PV->U][cellIdG1] = (uInner - help * dxHelp);
      m_cells->pvariables[PV->V][cellIdG1] = (vInner - help * dyHelp);
      m_cells->pvariables[PV->P][cellIdG1] = p;
      if(isRans) {
        m_cells->pvariables[PV->RANS_VAR[0]][cellIdG1] =
            (F2 * m_cells->pvariables[PV->RANS_VAR[0]][cellIdA1] - m_cells->pvariables[PV->RANS_VAR[0]][cellIdA2]);
      }
    }

    // extrapolate into second ghost cell
    for(MInt var = 0; var < PV->noVariables; var++) {
      m_cells->pvariables[var][cellIdG2] = F2 * m_cells->pvariables[var][cellIdG1] - m_cells->pvariables[var][cellIdA1];
    }
  }
}

/**
 * \brief Subsonic Outflow
 * extrapolate all but pressure, prescribe p8
 *
 *  rho=rho, u=u, v=v, w=w, p=p_inf
 */
template <MBool isRans>
void StructuredBndryCnd2D<isRans>::bc2007(MInt bcId) {
  MInt* start = m_physicalBCMap[bcId]->start1;
  MInt* end = m_physicalBCMap[bcId]->end1;

  switch(m_physicalBCMap[bcId]->face) {
    case 1: {
      for(MInt j = start[1]; j < end[1]; j++) {
        MInt cellId = cellIndex(start[0], j);
        MInt cellIdadj = cellIndex(start[0] - 1, j);

        m_cells->pvariables[PV->RHO][cellId] = m_cells->pvariables[PV->RHO][cellIdadj];
        m_cells->pvariables[PV->U][cellId] = m_cells->pvariables[PV->U][cellIdadj];
        m_cells->pvariables[PV->V][cellId] = m_cells->pvariables[PV->V][cellIdadj];
        m_cells->pvariables[PV->P][cellId] = PV->PInfinity;

        if(isRans) {
          for(MInt ransVarId = 0; ransVarId < m_solver->m_noRansEquations; ++ransVarId) {
            m_cells->pvariables[PV->RANS_VAR[ransVarId]][cellId] =
                m_cells->pvariables[PV->RANS_VAR[ransVarId]][cellIdadj];
          }
        }

        for(MInt var = 0; var < PV->noVariables; var++) {
          MInt cellIdP1 = cellIndex(start[0] + 1, j);
          m_cells->pvariables[var][cellIdP1] =
              2.0 * m_cells->pvariables[var][cellId] - m_cells->pvariables[var][cellIdadj];
        }
      }
      break;
    }
    case 3: {
      for(MInt i = start[0]; i < end[0]; i++) {
        MInt cellId = cellIndex(i, start[1]);
        MInt cellIdadj = cellIndex(i, start[1] - 1);

        m_cells->pvariables[PV->RHO][cellId] = m_cells->pvariables[PV->RHO][cellIdadj];
        m_cells->pvariables[PV->U][cellId] = m_cells->pvariables[PV->U][cellIdadj];
        m_cells->pvariables[PV->V][cellId] = m_cells->pvariables[PV->V][cellIdadj];
        m_cells->pvariables[PV->P][cellId] = PV->PInfinity;

        if(isRans) {
          for(MInt ransVarId = 0; ransVarId < m_solver->m_noRansEquations; ++ransVarId) {
            m_cells->pvariables[PV->RANS_VAR[ransVarId]][cellId] =
                m_cells->pvariables[PV->RANS_VAR[ransVarId]][cellIdadj];
          }
        }

        for(MInt var = 0; var < PV->noVariables; var++) {
          MInt cellIdP1 = cellIndex(i, start[1] + 1);
          m_cells->pvariables[var][cellIdP1] =
              2.0 * m_cells->pvariables[var][cellId] - m_cells->pvariables[var][cellIdadj];
        }
      }

      break;
    }

    default: {
      mTerm(1, AT_, "Face direction not implemented");
    }
  }
}

/**
 * \brief Symmetry plane BC
 */
template <MBool isRans>
void StructuredBndryCnd2D<isRans>::bc3000(MInt bcId) {
  TRACE();
  MInt* start = m_physicalBCMap[bcId]->start1;
  MInt* end = m_physicalBCMap[bcId]->end1;

  switch(m_physicalBCMap[bcId]->face) {
    case 0: {
      mTerm(1, "Not implemted yet!");
      break;
    }
    case 1: {
      mTerm(1, "Not implemted yet!");
      break;
    }
    case 2: {
      MInt cellId = -1;
      MInt cellIdadj = -1;
      const MInt cellShift = 2 * m_noGhostLayers - 1;

      for(MInt j = start[1]; j < end[1]; j++) {
        for(MInt i = start[0]; i < end[0]; i++) {
          cellId = cellIndex(i, j);                // ghost
          cellIdadj = cellIndex(i, cellShift - j); // field

          m_cells->pvariables[PV->RHO][cellId] = m_cells->pvariables[PV->RHO][cellIdadj];
          m_cells->pvariables[PV->P][cellId] = m_cells->pvariables[PV->P][cellIdadj];
          m_cells->pvariables[PV->U][cellId] = m_cells->pvariables[PV->U][cellIdadj];
          m_cells->pvariables[PV->V][cellId] = -m_cells->pvariables[PV->V][cellIdadj];
          if(isRans) {
            for(MInt ransVarId = 0; ransVarId < m_solver->m_noRansEquations; ++ransVarId) {
              m_cells->pvariables[PV->RANS_VAR[ransVarId]][cellId] =
                  m_cells->pvariables[PV->RANS_VAR[ransVarId]][cellIdadj];
            }
          }
        }
      }

      break;
    }
    case 3: {
      MInt cellId = -1;
      MInt cellIdadj = -1;
      const MInt cellShift = 2 * (m_nCells[0] - 1) - 2 * m_noGhostLayers + 1;

      for(MInt j = start[1]; j < end[1]; j++) {
        for(MInt i = start[0]; i < end[0]; i++) {
          // mirroring
          cellId = cellIndex(i, j);                // ghost
          cellIdadj = cellIndex(i, cellShift - j); // field

          m_cells->pvariables[PV->RHO][cellId] = m_cells->pvariables[PV->RHO][cellIdadj];
          m_cells->pvariables[PV->P][cellId] = m_cells->pvariables[PV->P][cellIdadj];
          m_cells->pvariables[PV->U][cellId] = m_cells->pvariables[PV->U][cellIdadj];
          m_cells->pvariables[PV->V][cellId] = -m_cells->pvariables[PV->V][cellIdadj];
          if(isRans) {
            for(MInt ransVarId = 0; ransVarId < m_solver->m_noRansEquations; ++ransVarId) {
              m_cells->pvariables[PV->RANS_VAR[ransVarId]][cellId] =
                  m_cells->pvariables[PV->RANS_VAR[ransVarId]][cellIdadj];
            }
          }
        }
      }
      break;
    }

    default: {
      mTerm(1, AT_, "Face direction not implemented)");
    }
  }
}

// shear flow inflow
template <MBool isRans>
void StructuredBndryCnd2D<isRans>::bc2021(MInt bcId) {
  TRACE();

  MInt* start = m_physicalBCMap[bcId]->start1;
  MInt* end = m_physicalBCMap[bcId]->end1;
  switch(m_physicalBCMap[bcId]->face) {
    case 0: {
      for(MInt j = start[1]; j < end[1]; j++) {
        const MInt cellIdG1 = cellIndex(1, j);
        const MInt cellIdG2 = cellIndex(0, j);
        const MInt cellIdA1 = cellIndex(2, j);
        // const MInt cellIdA2 = cellIndex(3,j);

        const MFloat dxidx = m_cells->surfaceMetrics[0][cellIdA1];
        const MFloat dxidy = m_cells->surfaceMetrics[1][cellIdA1];

        // multiply with n, so it will be -1 or +1 depending if we enter
        // or leave the domain of integration in positive direction
        const MFloat gradxi = -1 * F1 / sqrt(dxidx * dxidx + dxidy * dxidy);

        const MFloat dxHelp = dxidx * gradxi;
        const MFloat dyHelp = dxidy * gradxi;

        const MFloat cBC = sqrt(m_solver->m_gamma * pressure(cellIdG1) / m_cells->pvariables[PV->RHO][cellIdG1]);
        const MFloat rhoBC = m_cells->pvariables[PV->RHO][cellIdG1];

        const MFloat uInner = m_cells->pvariables[PV->U][cellIdA1];
        const MFloat vInner = m_cells->pvariables[PV->V][cellIdA1];
        const MFloat pInner = pressure(cellIdA1);

        const MFloat uInflow = PV->UInfinity * (m_cells->coordinates[1][cellIdG1] * m_bc2021Gradient);
        const MFloat vInflow = 0.0;

        // inflow
        const MFloat p =
            F1B2 * (pInner + PV->PInfinity + rhoBC * cBC * (dxHelp * (uInner - uInflow) + dyHelp * (vInner - vInflow)));

        const MFloat rho = CV->rhoInfinity + (p - PV->PInfinity) / POW2(cBC);
        const MFloat help = (p - PV->PInfinity) / (rhoBC * cBC);

        m_cells->pvariables[PV->RHO][cellIdG1] = rho;
        m_cells->pvariables[PV->U][cellIdG1] = uInflow + help * dxHelp;
        m_cells->pvariables[PV->V][cellIdG1] = vInflow + help * dyHelp;
        m_cells->pvariables[PV->P][cellIdG1] = p;

        if(isRans) {
          for(MInt ransVarId = 0; ransVarId < m_solver->m_noRansEquations; ++ransVarId) {
            m_cells->pvariables[PV->RANS_VAR[ransVarId]][cellIdG1] = PV->ransInfinity[ransVarId];
          }
        }

        // extrapolate into second ghost cell
        for(MInt var = 0; var < PV->noVariables; var++) {
          m_cells->pvariables[var][cellIdG2] =
              F2 * m_cells->pvariables[var][cellIdG1] - m_cells->pvariables[var][cellIdA1];
        }
      }
      break;
    }
    default:

    {
      mTerm(1, AT_, "Face direction not implemented)");
    }
  }
}

template <MBool isRans>
void StructuredBndryCnd2D<isRans>::bc2199(MInt bcId) {
  TRACE();
  //!!!!!!!!!! check for the gamma stuff in tfs
  const MFloat gamma = m_solver->m_gamma;
  const MFloat gammaMinusOne = m_solver->m_gamma - F1;
  const MFloat fgamma = F1 / gamma;
  MInt* start = m_physicalBCMap[bcId]->start1;
  MInt* end = m_physicalBCMap[bcId]->end1;
  switch(m_physicalBCMap[bcId]->face) {
    case 0: {
      MFloat mflux = F0;
      MFloat surface = F0;
      for(MInt j = start[1]; j < end[1]; j++) {
        const MInt cellId = cellIndex(m_noGhostLayers, j); // first inner line
        const MInt pointId = getPointIdFromCell(m_noGhostLayers, j);
        const MInt pointIdP1 = getPointIdFromPoint(pointId, 0, 1);
        //!!!!!!!!check the sign convention from tfs
        const MFloat dx = m_grid->m_coordinates[0][pointIdP1] - m_grid->m_coordinates[0][pointId];
        const MFloat dy = m_grid->m_coordinates[1][pointIdP1] - m_grid->m_coordinates[1][pointId];
        const MFloat rhou = m_cells->pvariables[PV->RHO][cellId] * m_cells->pvariables[PV->U][cellId];
        const MFloat rhov = m_cells->pvariables[PV->RHO][cellId] * m_cells->pvariables[PV->V][cellId];
        mflux += dx * rhou + dy * rhov;
        surface += sqrt(POW2(dx) + POW2(dy));
      }
      // get cellId at the domain top
      //!!!!!!! check for ii1
      MInt cellIdII1 = cellIndex(m_noGhostLayers, end[1] - m_noGhostLayers);
      // get averaged flux
      const MFloat rhouAVG = mflux / surface;
      // determine the pressure at the top at the iner line
      MFloat pressure = min(F1, max(F0, m_cells->pvariables[PV->P][cellIdII1] * gamma));
      // fix poin iteration apparently converges fast
      // from schlichting and trockenbrot, page 155
      // ru=sqrt(2.*fgamm1)*pp**fgam*sqrt(f1-pp**(gamm1*fgam))
      for(MInt i = 0; i < 20; i++) {
        pressure = pow((F1 - pow((rhouAVG / (sqrt(F2 * gammaMinusOne) * pow(pressure, fgamma))), F2)),
                       (gammaMinusOne * fgamma));
      }
      pressure *= fgamma;

      // scale the velocity profile

      break;
    }
    default: {
      mTerm(1, AT_, "Face direction not implemented)");
    }
  }
}


template <MBool isRans>
void StructuredBndryCnd2D<isRans>::bc2402(MInt bcId) {
  TRACE();

  MInt* start = m_physicalBCMap[bcId]->start1;
  MInt* end = m_physicalBCMap[bcId]->end1;

  // first we need the average pressure and Temperature
  //==> calculate from the field
  //==> average over the surface
  //==> store in variable

  MInt Id = m_channelSurfaceIndexMap[bcId];
  MInt* startface = &m_solver->m_windowInfo->channelSurfaceIndices[Id]->start1[0];
  MInt* endface = &m_solver->m_windowInfo->channelSurfaceIndices[Id]->end1[0];
  // global Pressure and Temperature at in-/outlet
  MFloat globalTin[2] = {F0, F0};
  MFloat globalTout[2] = {F0, F0};
  MFloat globalPin[2] = {F0, F0};
  MFloat globalPout[2] = {F0, F0};
  cout.precision(10);

  // find out which face it is
  switch(m_physicalBCMap[bcId]->face) {
    case 0: {
      // average the pressure
      // use values from the inside to determine the pressure at the face!!!
      MFloat surface = F0;
      MFloat localPin[2] = {F0, F0};
      MFloat localTin[2] = {F0, F0};
      MFloat localVel = F0, globalVel = F0;
      MFloat localMassFlux = F0, globalMassFlux = F0;

      for(MInt j = startface[1]; j < endface[1]; j++) {
        MInt ii = end[0] - 1;
        MInt cellIdP1 = cellIndex(ii + 1, j);
        MInt cellIdP2 = cellIndex(ii + 2, j);

        surface = sqrt(POW2(m_cells->surfaceMetrics[0][cellIdP1]) + POW2(m_cells->surfaceMetrics[1][cellIdP1]));
        localPin[0] += surface * m_cells->pvariables[PV->P][cellIdP1];
        localPin[1] += surface * m_cells->pvariables[PV->P][cellIdP2];
        localTin[0] += surface * temperature(cellIdP1);
        localTin[1] += surface * temperature(cellIdP2);
        localVel += surface * (m_cells->pvariables[PV->U][cellIdP1]);
        localMassFlux += surface * (m_cells->pvariables[PV->U][cellIdP1] * m_cells->pvariables[PV->RHO][cellIdP1]);
      }

      // next now that the pressure and temperature are known:
      // make it a global variable for the complete inflow plane
      MPI_Allreduce(localPin, globalPin, 2, MPI_DOUBLE, MPI_SUM, m_solver->m_commChannelIn, AT_, "localPin",
                    "globalPin");
      MPI_Allreduce(localTin, globalTin, 2, MPI_DOUBLE, MPI_SUM, m_solver->m_commChannelIn, AT_, "localTin",
                    "globalTin");
      MPI_Allreduce(&localMassFlux, &globalMassFlux, 1, MPI_DOUBLE, MPI_SUM, m_solver->m_commChannelIn, AT_,
                    "localMassFlux", "globalMassFlux");
      MPI_Allreduce(&localVel, &globalVel, 1, MPI_DOUBLE, MPI_SUM, m_solver->m_commChannelIn, AT_, "localVel",
                    "globalVel");

      MPI_Barrier(m_solver->m_commChannelIn, AT_);
      globalPin[0] /= m_channelSurfaceIn;
      globalTin[0] /= m_channelSurfaceIn;
      globalPin[1] /= m_channelSurfaceIn;
      globalTin[1] /= m_channelSurfaceIn;
      globalMassFlux /= m_channelSurfaceIn;
      globalVel /= m_channelSurfaceIn;
      MFloat currentRe = m_solver->m_Re / PV->UInfinity * globalVel;

      // set a Barrier to ensure that all domains having a channel side are at same time level
      // now make the variables known everywhere in the channel world
      MPI_Barrier(m_solver->m_commChannelWorld, AT_);
      MPI_Bcast(globalTin, 2, MPI_DOUBLE, m_solver->m_channelRoots[2], m_solver->m_commChannelWorld, AT_, "globalTin");
      MPI_Bcast(globalPin, 2, MPI_DOUBLE, m_solver->m_channelRoots[2], m_solver->m_commChannelWorld, AT_, "globalPin");
      MPI_Bcast(globalTout, 2, MPI_DOUBLE, m_solver->m_channelRoots[3], m_solver->m_commChannelWorld, AT_,
                "globalTout");
      MPI_Bcast(globalPout, 2, MPI_DOUBLE, m_solver->m_channelRoots[3], m_solver->m_commChannelWorld, AT_,
                "globalPout");
      // now every value is known and can be used to apply the BC !!!!
      // cannot take values from above at here the window is extended to the ghost layers

      if(globalTimeStep > 1 && globalTimeStep % 50 == 0) {
        if(m_channelInflowRank == 0 && globalTimeStep % m_solver->m_residualInterval == 0 && m_solver->m_RKStep == 0) {
          cout.precision(6);

          FILE* f_channel;
          f_channel = fopen("./massflow.dat", "a+");
          fprintf(f_channel, "%d", globalTimeStep);
          fprintf(f_channel, " %f", m_solver->m_physicalTime);
          fprintf(f_channel, " %f", m_solver->m_time);
          fprintf(f_channel, " %f", m_solver->m_timeStep);
          fprintf(f_channel, " %f", currentRe);
          fprintf(f_channel, " %f", globalMassFlux);
          fprintf(f_channel, " %f", globalPin[0]);
          fprintf(f_channel, " %f", globalPin[1]);
          fprintf(f_channel, " %f", globalTin[0]);
          fprintf(f_channel, " %f", globalTout[0]);
          fprintf(f_channel, " %f", globalPout[0]);
          fprintf(f_channel, " %f", globalPout[1]);
          fprintf(f_channel, " %f", (globalTin[0] - globalTout[1]));
          fprintf(f_channel, "\n");
          fclose(f_channel);
        }
      }

      // If fully periodic, do nothing
      if(m_solver->m_channelFullyPeriodic) {
        m_solver->m_inflowVelAvg = globalVel;
        return;
      }

      for(MInt j = start[1]; j < end[1]; j++) {
        MInt i = 1;
        const MInt cellId = cellIndex(i, j);

        MFloat pressureFluctuation = m_cells->pvariables[PV->P][cellId] - globalPout[1];
        MFloat x = m_cells->coordinates[0][cellId];
        MFloat pressureInflowMean =
            m_solver->m_deltaP / m_solver->m_channelLength * (x - m_solver->m_channelInflowPlaneCoordinate)
            + PV->PInfinity;
        MFloat pressureNew = pressureInflowMean + pressureFluctuation;
        MFloat temperatureFlucOutflow = temperature(cellId) - globalTout[1];
        MFloat temperatureNew = PV->TInfinity + temperatureFlucOutflow;
        MFloat rhoNew = m_solver->m_gamma * pressureNew / temperatureNew;

        m_cells->pvariables[PV->RHO][cellId] = rhoNew;
        m_cells->pvariables[PV->P][cellId] = pressureNew;

        for(MInt var = 0; var < PV->noVariables; var++) {
          m_cells->pvariables[var][cellIndex(start[0], j)] = 2.0 * m_cells->pvariables[var][cellIndex(start[0] + 1, j)]
                                                             - m_cells->pvariables[var][cellIndex(start[0] + 2, j)];
        }
      }
      break;
    }
    case 1: {
      // average the pressure
      // use values from the inside to determine the pressure at the face!!!
      MFloat surface = F0;
      MFloat localPout[] = {F0, F0};
      MFloat localTout[] = {F0, F0};
      for(MInt j = startface[1]; j < endface[1]; j++) {
        MInt ii = start[0];
        MInt cellIdM1 = cellIndex(ii - 1, j);
        MInt cellIdM2 = cellIndex(ii - 2, j);
        surface = sqrt(POW2(m_cells->surfaceMetrics[0][cellIdM1]) + POW2(m_cells->surfaceMetrics[1][cellIdM1]));
        localPout[0] += surface * m_cells->pvariables[PV->P][cellIdM2];
        localPout[1] += surface * m_cells->pvariables[PV->P][cellIdM1];
        localTout[0] += surface * temperature(cellIdM2);
        localTout[1] += surface * temperature(cellIdM1);
      }

      // next now that the pressure and temperature are known:
      // make it a global variable for the complete outflow plane
      MPI_Allreduce(localPout, globalPout, 2, MPI_DOUBLE, MPI_SUM, m_solver->m_commChannelOut, AT_, "localPout",
                    "globalPout");
      MPI_Allreduce(localTout, globalTout, 2, MPI_DOUBLE, MPI_SUM, m_solver->m_commChannelOut, AT_, "localTout",
                    "globalTout");

      MPI_Barrier(m_solver->m_commChannelOut, AT_);
      globalPout[0] /= m_channelSurfaceOut;
      globalTout[0] /= m_channelSurfaceOut;
      globalPout[1] /= m_channelSurfaceOut;
      globalTout[1] /= m_channelSurfaceOut;
      // set a Barrier to ensure that all domains having a channel side are at same time level
      // now make the variables known everywhere in the channel world
      MPI_Barrier(m_solver->m_commChannelWorld, AT_);
      MPI_Bcast(globalTin, 2, MPI_DOUBLE, m_solver->m_channelRoots[2], m_solver->m_commChannelWorld, AT_, "globalTin");
      MPI_Bcast(globalPin, 2, MPI_DOUBLE, m_solver->m_channelRoots[2], m_solver->m_commChannelWorld, AT_, "globalPin");
      MPI_Bcast(globalTout, 2, MPI_DOUBLE, m_solver->m_channelRoots[3], m_solver->m_commChannelWorld, AT_,
                "globalTout");
      MPI_Bcast(globalPout, 2, MPI_DOUBLE, m_solver->m_channelRoots[3], m_solver->m_commChannelWorld, AT_,
                "globalPout");
      // now every value is known and can be used to apply the BC !!!!
      // cannot take values from above at here the window is extended to the ghost layers

      // If fully periodic, do nothing; return after mpi-stuff to avoid blocking
      if(m_solver->m_channelFullyPeriodic) {
        return;
      }

      for(MInt j = start[1]; j < end[1]; j++) {
        MInt i = start[0];
        const MInt cellId = cellIndex(i, j);

        MFloat pressureFluctuation = m_cells->pvariables[PV->P][cellId] - globalPin[0];
        MFloat x = m_cells->coordinates[0][cellId];
        MFloat pressureOutflowMean =
            m_solver->m_deltaP / m_solver->m_channelLength * (x - m_solver->m_channelInflowPlaneCoordinate)
            + PV->PInfinity;
        MFloat pressureNew = pressureOutflowMean + pressureFluctuation;
        MFloat deltaT = (globalTin[0] - globalTout[1]);
        MFloat temperatureInflow = temperature(cellId);
        MFloat temperatureNew = temperatureInflow - deltaT;
        MFloat rhoNew = m_solver->m_gamma * pressureNew / temperatureNew;

        m_cells->pvariables[PV->RHO][cellId] = rhoNew;
        m_cells->pvariables[PV->P][cellId] = pressureNew;

        for(MInt var = 0; var < PV->noVariables; var++) {
          m_cells->pvariables[var][cellIndex(start[0] + 1, j)] = 2.0 * m_cells->pvariables[var][cellIndex(start[0], j)]
                                                                 - m_cells->pvariables[var][cellIndex(start[0] - 1, j)];
        }
      }
      break;
    }
    default: {
      mTerm(1, AT_, "Face not implemented");
      break;
    }
  }
}


// Inlet station for RANS
template <MBool isRans>
void StructuredBndryCnd2D<isRans>::bc2510(MInt bcId) {
  (void)bcId;

  cout.precision(8);
  // MInt* start = m_physicalBCMap[bcId]->start1;
  // MInt* end = m_physicalBCMap[bcId]->end1;
  // Communication between the Recycling and Inlet station.
  //___COMMGr(oup)______________________
  //|______________|_________________| |
  //|Inlet station |Recycling station| |
  //|              |                 | |
  //|______________|_________________| |
  //|__________________________________|
  //
  // infographic of the Groups

  // MFloat delta_inMax = delta_in+2.5*delta_in; //limit the integration height
  const MFloat gamma = m_solver->m_gamma;
  const MFloat gammaMinusOne = gamma - F1;

  const MFloat rescalEPS = pow(10, -16.0);
  const MFloat alpha = 4.0;
  const MFloat b = 0.2;
  const MFloat rc = pow(m_solver->m_Pr, F1B3);
  const MFloat ctema = F1B2 * gammaMinusOne * POW2(m_solver->m_Ma) * rc;
  const MFloat maxIntegrationHeight = 2.0 * m_rescalingBLT;

  // van Driest constant & transformed velocity
  const MFloat b_vd = sqrt(ctema / (F1 + ctema));
  const MFloat uvd8 = PV->UInfinity * asin(b_vd) / b_vd;
  const MInt i = m_noGhostLayers - 1;

  //////////////////////////////////////////////////////////////
  ///////////////// COMPUTE THETA AND EXCHANGE /////////////////
  //////////////////////////////////////////////////////////////

  // allocate space in k direction (all k-Cells)
  MFloatScratchSpace thetaLocal(2, AT_, "thetaLocalIn");
  thetaLocal.fill(F0); // initialize scratch space

  MFloatScratchSpace thetaGlobal(2, AT_, "thetaGlobalIn");
  thetaGlobal.fill(F0);

  // compute the local moment thickness j=direction of integration
  for(MInt j = m_noGhostLayers; j < m_nCells[0] - m_noGhostLayers; ++j) {
    const MInt cellId = cellIndex(i, j);
    const MInt pointIdM1 = getPointIdFromCell(i, j);
    const MInt pointIdP1 = getPointIdFromPoint(pointIdM1, 0, 1);

    if(m_grid->m_coordinates[1][pointIdM1] > maxIntegrationHeight) {
      continue;
    }

    const MFloat urat = m_cells->pvariables[PV->U][cellId] / PV->UInfinity;
    const MFloat momThick = m_cells->pvariables[PV->U][cellId] * m_cells->pvariables[PV->RHO][cellId] * fabs(F1 - urat)
                            / (CV->rhoUInfinity);
    // integrate normal to the wall
    const MFloat ydist = m_grid->m_coordinates[1][pointIdP1] - m_grid->m_coordinates[1][pointIdM1];
    thetaLocal(0) += momThick * ydist;
  }

  MPI_Allreduce(thetaLocal.begin(), thetaGlobal.begin(), 2, MPI_DOUBLE, MPI_SUM, m_solver->m_rescalingCommGrComm, AT_,
                "thetaLocal.begin()", "thetaGlobal.begin()");

  thetaGlobal(0) = mMax(thetaGlobal(0), 0.0000001);
  thetaGlobal(1) = mMax(thetaGlobal(1), 0.0000001);

  if(globalTimeStep % 50 == 0 && m_solver->m_RKStep == 0
     && m_solver->domainId() == m_solver->m_rescalingCommGrRootGlobal) {
    cout << m_solver->domainId() << " ThetaInflow " << thetaGlobal(0) << " ThetaRecyclingStation " << thetaGlobal(1)
         << endl;

    FILE* f_channel;
    f_channel = fopen("./theta_inflow.dat", "a+");
    fprintf(f_channel, "%d", globalTimeStep);
    fprintf(f_channel, " %f", m_solver->m_physicalTime);
    fprintf(f_channel, " %f", m_solver->m_time);
    fprintf(f_channel, " %f", m_solver->m_timeStep);
    fprintf(f_channel, " %f", thetaGlobal[0]);
    fprintf(f_channel, " %f", thetaGlobal[1]);
    fprintf(f_channel, "\n");
    fclose(f_channel);
  }

  //////////////////////////////////////////////////////////////
  ///////////////// EXCHANGE WALL PROPERTIES ///////////////////
  //////////////////////////////////////////////////////////////
  const MInt noWallProperties = 3;
  MFloatScratchSpace wallPropertiesLocal(noWallProperties, AT_, "wallPropertiesLocalInlet");
  MFloatScratchSpace wallProperties(noWallProperties, AT_, "wallPropertiesInlet");
  wallPropertiesLocal.fill(F0);
  wallProperties.fill(F0);

  MPI_Allreduce(wallPropertiesLocal.begin(), wallProperties.begin(), noWallProperties, MPI_DOUBLE, MPI_SUM,
                m_solver->m_rescalingCommGrComm, AT_, "wallPropertiesLocal.begin()", "wallProperties.begin()");

  //////////////////////////////////////////////////////////////
  ///////////////// GAMS, UTAUIN, INNER OUTER COORD ////////////
  //////////////////////////////////////////////////////////////

  MFloat utauIn = F0;
  MFloat gams = F0;

  MFloat utauRe = wallProperties(0);
  // estimate the friction velocity at the inlet
  // according to standard power law approximations
  // utau_in = utau_re*(theta_re/theta_in)**(1/2*(n-1))
  // where theta is the momentum thichness
  // see Thomas S. Lund, p241

  // when take into account the variance of wall density
  // utau_in = utau_re*(rho_wall_re/rho_wall_in)**0.5
  //* (theta_re/theta_in)**(1/2*(n-1))
  // here n = 5
  const MFloat n = 5.0;
  const MFloat facc = F1 / (F2 * (n - F1));

  gams = pow(thetaGlobal(1) / fabs(thetaGlobal(0)), facc);
  gams = mMin(gams, 1.8);
  utauIn = utauRe * gams;

  MFloatScratchSpace coordInInner(m_nCells[0], AT_, "coordInInner");
  MFloatScratchSpace coordInOuter(m_nCells[0], AT_, "coordInOuter");

  for(MInt j = 0; j < m_nCells[0]; ++j) {
    const MInt cellId = cellIndex(i, j);
    const MFloat rho = m_cells->pvariables[PV->RHO][cellId];
    const MFloat frho = F1 / rho;
    const MFloat p = m_cells->pvariables[PV->P][cellId];
    const MFloat temp = p * gamma * frho;
    const MFloat mu = SUTHERLANDLAW(temp);

    coordInInner(j) = utauIn * rho * m_cells->coordinates[1][cellId] / (mu * sqrt(m_solver->m_Re0));
    coordInOuter(j) = m_cells->coordinates[1][cellId] * rho / (m_rescalingBLT * CV->rhoInfinity);
  }

  const MInt wallLocalOffset = m_solver->m_nOffsetCells[0]; // Offset in j-direction
  MFloat tempWallInletLocal = F0;
  MFloat tempWallInletGlobal = F0;
  // determine the wall stuff if wall is contained whithin the partition
  if(wallLocalOffset == 0 && m_solver->m_nActiveCells[0] >= m_noGhostLayers) {
    const MInt cellId = cellIndex(i, m_noGhostLayers);
    tempWallInletLocal = temperature(cellId);
  }

  MPI_Allreduce(&tempWallInletLocal, &tempWallInletGlobal, 1, MPI_DOUBLE, MPI_SUM, m_solver->m_rescalingCommGrComm, AT_,
                "tempWallInletLocal", "tempWallInletGlobal");

  //////////////////////////////////////////////////////////////
  ///////////////// NOW EXCHANGE VAR SLICE /////////////////////
  //////////////////////////////////////////////////////////////
  const MInt noVariables = PV->noVariables + 1;
  MInt totalCells = m_grid->getMyBlockNoCells(0) + 1;
  MFloatScratchSpace varSliceLocal(noVariables, totalCells, AT_, "varSliceLocal");
  MFloatScratchSpace varSlice(noVariables, totalCells, AT_, "varSlice");

  // we are at the inlet, only fill with zeros
  varSlice.fill(F0);
  varSliceLocal.fill(F0);

  MPI_Allreduce(varSliceLocal.begin(), varSlice.begin(), noVariables * totalCells, MPI_DOUBLE, MPI_SUM,
                m_solver->m_rescalingCommGrComm, AT_, "varSliceLocal.begin()", "varSlice.begin()");

  // if(globalTimeStep%5==0&&m_solver->m_RKStep==0&&m_solver->domainId() == m_solver->m_rescalingCommGrRootGlobal) {
  // cout << m_solver->domainId()<< " ThetaInflow " << thetaGlobal(0) <<" ThetaRecyclingStation " << thetaGlobal(1) <<
  // endl;

  // FILE* f_channel;
  // f_channel = fopen("./varslice.dat", "w");
  //          for(MInt jj=0; jj<totalCells-1; ++jj) {
  //            for(MInt var=0; var<6; var++) {
  //              fprintf(f_channel, " %f", varSlice(var, jj));
  //            }
  //            fprintf(f_channel, "\n");
  //          }
  // fclose(f_channel);
  // }

  ///////////////////////////////////////////////////////////////
  ///////////////// RESCALING ///////////////////////////////////
  ///////////////////////////////////////////////////////////////

  const MFloat ctem1 = (F1 + ctema) * (F1 - POW2(gams));
  const MFloat ctem2 = F2 * ctema * gams * (F1 - gams);
  const MFloat ctem3 = (F1 - gams) * (F1 + gams + F2 * ctema * gams);

  MInt jStart = 0;
  MInt jEnd = m_nCells[0];

  if(m_solver->m_nOffsetCells[0] == 0) {
    jStart = m_noGhostLayers;
  }
  if(m_solver->m_nOffsetCells[0] + m_solver->m_nActiveCells[0] == m_grid->getMyBlockNoCells(0)) {
    jEnd = m_nCells[0] - m_noGhostLayers;
  }

  MFloat blEdgeVValueLocal = F0;
  MBool edgePointIsSet = false;
  MInt edgePointJ = 0;


  for(MInt j = jStart; j < jEnd; ++j) {
    const MInt cellId = cellIndex(i, j);

    if(j > 0) {
      if(coordInOuter(j - 1) < 1.05 && coordInOuter(j) >= 1.05) {
        blEdgeVValueLocal = m_cells->pvariables[PV->V][cellIndex(i, j - 1)];
      }
    }

    if(coordInOuter(j) < 1.05) {
      MFloat uInner = F0, vInner = F0, TInner = F0, mutInner = F0;
      MFloat uOuter = F0, vOuter = F0, TOuter = F0, mutOuter = F0;
      const MFloat count = alpha * (coordInOuter(j) - b);
      const MFloat denom = (F1 - F2 * b) * coordInOuter(j) + b;
      const MFloat ratio = count / denom;
      const MFloat wfun = F1B2 * (F1 + tanh(ratio) / tanh(alpha));

      for(MInt jj = 0; jj < totalCells - 1; ++jj) {
        const MInt localId = jj;
        const MInt localIdP1 = jj + 1;

        const MFloat yInnerRe = varSlice(3, localId);
        const MFloat yInnerReP1 = varSlice(3, localIdP1);

        if((yInnerRe - coordInInner(j)) < rescalEPS && yInnerReP1 > coordInInner(j)) {
          const MFloat dy1 = coordInInner(j) - yInnerRe;
          const MFloat dy2 = yInnerReP1 - coordInInner(j);
          const MFloat dy = yInnerReP1 - yInnerRe;

          const MFloat u = varSlice(0, localId);
          const MFloat uP1 = varSlice(0, localIdP1);
          const MFloat v = varSlice(1, localId);
          const MFloat vP1 = varSlice(1, localIdP1);
          const MFloat t = varSlice(2, localId);
          const MFloat tP1 = varSlice(2, localIdP1);
          const MFloat mut = varSlice(5, localId);
          const MFloat mutP1 = varSlice(5, localIdP1);
          uInner = (uP1 * dy1 + u * dy2) / dy;
          vInner = (vP1 * dy1 + v * dy2) / dy;
          TInner = (tP1 * dy1 + t * dy2) / dy;
          mutInner = (mutP1 * dy1 + mut * dy2) / dy;
        }
      }

      // outer region
      for(MInt jj = 0; jj < totalCells - 1; ++jj) {
        const MInt localId = jj;
        const MInt localIdP1 = jj + 1;

        const MFloat yOuterRe = varSlice(4, localId);
        const MFloat yOuterReP1 = varSlice(4, localIdP1);

        if((yOuterRe - coordInOuter(j)) < rescalEPS && yOuterReP1 > coordInOuter(j)) {
          const MFloat dy1 = coordInOuter(j) - yOuterRe;
          const MFloat dy2 = yOuterReP1 - coordInOuter(j);
          const MFloat dy = yOuterReP1 - yOuterRe;

          const MFloat u = varSlice(0, localId);
          const MFloat uP1 = varSlice(0, localIdP1);
          const MFloat v = varSlice(1, localId);
          const MFloat vP1 = varSlice(1, localIdP1);
          const MFloat t = varSlice(2, localId);
          const MFloat tP1 = varSlice(2, localIdP1);
          const MFloat mut = varSlice(5, localId);
          const MFloat mutP1 = varSlice(5, localIdP1);
          uOuter = (uP1 * dy1 + u * dy2) / dy;
          vOuter = (vP1 * dy1 + v * dy2) / dy;
          TOuter = (tP1 * dy1 + t * dy2) / dy;
          mutOuter = (mutP1 * dy1 + mut * dy2) / dy;
        }
      }

      const MFloat TInnerA = POW2(gams) * TInner + ctem1 * PV->TInfinity;
      const MFloat TOuterA = POW2(gams) * TOuter - (ctem2 * (uOuter / PV->UInfinity) - ctem3) * PV->TInfinity;

      // van Driest transformation
      const MFloat uvdInner = PV->UInfinity * asin(b_vd * uInner / PV->UInfinity) / b_vd;
      const MFloat uvdOuter = PV->UInfinity * asin(b_vd * uOuter / PV->UInfinity) / b_vd;

      // scaling of transformed inner and outer velocities
      // use uvd8, the van Driest transformed u8 value
      uInner = gams * uvdInner;
      uOuter = gams * uvdOuter + (F1 - gams) * uvd8;

      uInner = PV->UInfinity * sin(b_vd * uInner / PV->UInfinity) / b_vd;
      uOuter = PV->UInfinity * sin(b_vd * uOuter / PV->UInfinity) / b_vd;

      // turbulent viscosity
      const MFloat tempWallInlet = tempWallInletGlobal;
      const MFloat tempWallRecycling = wallProperties(2);

      const MFloat viscWallInlet = SUTHERLANDLAW(tempWallInlet);
      const MFloat viscWallRecycling = SUTHERLANDLAW(tempWallRecycling);
      const MFloat thetaInlet = thetaGlobal(0);
      const MFloat thetaRecycling = thetaGlobal(1);
      mutInner = mutInner * (viscWallInlet / viscWallRecycling);
      mutOuter = mutOuter * gams * (thetaInlet / thetaRecycling);

      const MFloat uMean = uInner * (F1 - wfun) + uOuter * wfun;
      const MFloat vMean = vInner * (F1 - wfun) + vOuter * wfun;
      const MFloat tMean = TInnerA * (F1 - wfun) + TOuterA * wfun;
      MFloat mutMean = mutInner * (F1 - wfun) + mutOuter * wfun;

      const MFloat clebf = 6.6;
      const MFloat blt = m_rescalingBLT;
      const MFloat cleb = F1 / (F1 + pow((m_cells->coordinates[1][cellId] / (clebf * blt)), 6.0));

      const MFloat pres = PV->PInfinity;
      const MFloat rhoIn = gamma * pres / tMean;

      m_cells->pvariables[PV->RHO][cellId] = rhoIn * cleb;
      m_cells->pvariables[PV->U][cellId] = uMean;
      m_cells->pvariables[PV->V][cellId] = vMean;
      m_cells->pvariables[PV->P][cellId] = pres;
      m_cells->pvariables[PV->RANS_VAR[0]][cellId] = mutMean / rhoIn;

    } else {
      if(!edgePointIsSet) {
        edgePointJ = j;
        edgePointIsSet = true;
      }
      // const MFloat pres = pressure(cellIndex(m_noGhostLayers,j,k)); //for supersonic PV->PInfinity
      const MFloat pres = PV->PInfinity;
      const MFloat rhoIn = gamma * pres / PV->TInfinity;

      const MFloat uMean = PV->UInfinity;
      const MFloat vMean = PV->VInfinity;

      m_cells->pvariables[PV->RHO][cellId] = rhoIn;
      m_cells->pvariables[PV->U][cellId] = uMean;
      m_cells->pvariables[PV->V][cellId] = vMean; // m_cells->pvariables[PV->V][cellIndex(i,j-1)]*rhoIn;
      m_cells->pvariables[PV->P][cellId] = pres;
      m_cells->pvariables[PV->RANS_VAR[0]][cellId] = PV->ransInfinity[0];
    }
  }

  MFloat blEdgeVValueGlobal = F0;
  MPI_Allreduce(&blEdgeVValueLocal, &blEdgeVValueGlobal, 1, MPI_DOUBLE, MPI_SUM, m_solver->m_rescalingCommGrComm, AT_,
                "blEdgeVValueLocal", "blEdgeVValueGlobal");

  if(edgePointIsSet) {
    for(MInt j = edgePointJ; j < jEnd; ++j) {
      const MInt cellId = cellIndex(i, j);
      m_cells->pvariables[PV->V][cellId] = blEdgeVValueGlobal;
      const MFloat pres = PV->PInfinity;
      m_cells->pvariables[PV->P][cellId] = pres;
    }
  }

  for(MInt j = 0; j < m_nCells[0]; ++j) {
    // extrapolation for second GC
    const MInt cellId = cellIndex(1, j);
    const MInt cellIdM1 = cellIndex(0, j);
    const MInt cellIdadj = cellIndex(2, j);

    for(MInt var = 0; var < PV->noVariables; var++) {
      m_cells->pvariables[var][cellIdM1] = 2.0 * m_cells->pvariables[var][cellId] - m_cells->pvariables[var][cellIdadj];
    }
  }
}


// Recycling station for RANS
template <MBool isRans>
void StructuredBndryCnd2D<isRans>::bc2511(MInt bcId) {
  MInt* start = m_physicalBCMap[bcId]->start1;
  cout.precision(8);

  // scaling between delta0 and delta2
  const MFloat F727 = 72.0 / 7.0; // 8.5, 8.0
  const MInt i = start[0];
  const MFloat yWall = F0; // this has been fixed else method does not works
  const MFloat maxIntegrationHeight = 2.0 * m_rescalingBLT;

  //////////////////////////////////////////////////////////////
  ///////////////// COMPUTE THETA AND EXCHANGE /////////////////
  //////////////////////////////////////////////////////////////

  // thetaLocal.fill(F0); //initialize scratch space to zero // only for parallel use
  MFloatScratchSpace thetaLocal(2, AT_, "thetaLocalRe");
  MFloatScratchSpace thetaGlobal(2, AT_, "thetaGlobalRe");
  thetaLocal.fill(F0); // initialize scratch space
  thetaGlobal.fill(F0);

  // the offest position in k-direction is the offset
  // compute the local moment thickness j=direction of integration
  for(MInt j = m_noGhostLayers; j < m_nCells[0] - m_noGhostLayers; ++j) {
    const MInt cellId = cellIndex(i, j);
    const MInt pointIdM1 = getPointIdFromCell(i, j);
    const MInt pointIdP1 = getPointIdFromPoint(pointIdM1, 0, 1);

    if(m_grid->m_coordinates[1][pointIdM1] > maxIntegrationHeight) {
      continue;
    }

    const MFloat urat = m_cells->pvariables[PV->U][cellId] / PV->UInfinity;
    const MFloat momThick = m_cells->pvariables[PV->U][cellId] * m_cells->pvariables[PV->RHO][cellId] * fabs(F1 - urat)
                            / (CV->rhoUInfinity);

    // integrate normal to the wall
    const MFloat ydist = m_grid->m_coordinates[1][pointIdP1] - m_grid->m_coordinates[1][pointIdM1];
    thetaLocal(1) += momThick * ydist;
  }


  // communicate the Thickness across the plane
  MPI_Allreduce(thetaLocal.begin(), thetaGlobal.begin(), 2, MPI_DOUBLE, MPI_SUM, m_solver->m_rescalingCommGrComm, AT_,
                "thetaLocal.begin()", "thetaGlobal.begin()");

  //////////////////////////////////////////////////////////////
  ///////////////// COMPUTE AND EXCHANGE WALL PROPERTIES ///////
  //////////////////////////////////////////////////////////////

  const MFloat delta = F727 * thetaGlobal(1);
  const MInt noVar = 3;                                     // for more variables if wanted
  const MInt wallLocalOffset = m_solver->m_nOffsetCells[0]; // Offset in j-direction

  MFloatScratchSpace wallPropertiesLocal(noVar, AT_, "wallPropertiesLocalRe");
  MFloatScratchSpace wallProperties(noVar, AT_, "wallPropertiesRe");
  wallPropertiesLocal.fill(F0);
  wallProperties.fill(F0);

  // determine the wall stuff if wall is contained whithin the partition
  if(wallLocalOffset == 0 && m_solver->m_nActiveCells[1] >= m_noGhostLayers) {
    const MInt cellId = cellIndex(i, m_noGhostLayers);
    const MFloat rho = m_cells->pvariables[PV->RHO][cellId];
    const MFloat p = m_cells->pvariables[PV->P][cellId];
    const MFloat t = p * m_solver->m_gamma / rho;
    const MFloat mu = SUTHERLANDLAW(t);
    const MFloat uWall = fabs(m_cells->pvariables[PV->U][cellId]);
    const MFloat ydist = m_cells->coordinates[1][cellId] - yWall;
    const MFloat uTau = sqrt(uWall * mu / (ydist * rho));

    wallPropertiesLocal(0) = uTau;
    wallPropertiesLocal(1) = rho;
    wallPropertiesLocal(2) = t;
  }

  MPI_Allreduce(wallPropertiesLocal.begin(), wallProperties.begin(), noVar, MPI_DOUBLE, MPI_SUM,
                m_solver->m_rescalingCommGrComm, AT_, "wallPropertiesLocal.begin()", "wallProperties.begin()");

  MFloat tempWallInletLocal = F0;
  MFloat tempWallInletGlobal = F0;
  MPI_Allreduce(&tempWallInletLocal, &tempWallInletGlobal, 1, MPI_DOUBLE, MPI_SUM, m_solver->m_rescalingCommGrComm, AT_,
                "tempWallInletLocal", "tempWallInletGlobal");

  //////////////////////////////////////////////////////////////
  ///////////////// COMPUTE AND EXCHANGE VAR SLICE /////////////
  //////////////////////////////////////////////////////////////

  MInt totalCells = m_grid->getMyBlockNoCells(0) + 1;

  const MInt noVariables = PV->noVariables + 1;
  MFloatScratchSpace varSliceLocal(noVariables, totalCells, AT_, "varSliceLocal");
  MFloatScratchSpace varSlice(noVariables, totalCells, AT_, "varSlice");

  for(MInt j = m_noGhostLayers; j < m_nCells[0] - m_noGhostLayers; ++j) {
    const MInt cellId = cellIndex(i, j);
    const MFloat rho = m_cells->pvariables[PV->RHO][cellId];
    const MFloat frho = F1 / rho;
    const MFloat p = pressure(cellId);
    const MFloat temp = p * m_solver->m_gamma * frho;
    const MFloat mu = SUTHERLANDLAW(temp);
    const MFloat uTauRe = wallProperties(0);
    const MFloat yIn = (m_cells->coordinates[1][cellId] - yWall) * uTauRe * rho / (mu * sqrt(m_solver->m_Re0));
    const MFloat yOut = (m_cells->coordinates[1][cellId] - yWall) * rho / (delta * CV->rhoInfinity);
    const MFloat u = m_cells->pvariables[PV->U][cellId];
    const MFloat v = m_cells->pvariables[PV->V][cellId];

    //>RANS
    const MFloat mut = m_cells->pvariables[PV->RANS_VAR[0]][cellId] * rho;
    //<RANS

    const MInt localId = m_solver->m_nOffsetCells[0] + j - m_noGhostLayers + 1;

    // save the variables u,v,w,t,yI,yO
    varSliceLocal(0, localId) = u;
    varSliceLocal(1, localId) = v;
    varSliceLocal(2, localId) = temp;
    varSliceLocal(3, localId) = yIn;
    varSliceLocal(4, localId) = yOut;
    varSliceLocal(5, localId) = mut;
  }

  // set first value at the wall manually
  if(m_solver->m_nOffsetCells[0] == 0) {
    varSliceLocal(0, 0) = 0.0;
    varSliceLocal(1, 0) = 0.0;
    varSliceLocal(2, 0) = varSliceLocal(2, 1);
    varSliceLocal(3, 0) = 0.0;
    varSliceLocal(4, 0) = 0.0;
    varSliceLocal(5, 0) = varSliceLocal(5, 1);
  }

  // communicate the slice
  MPI_Allreduce(varSliceLocal.begin(), varSlice.begin(), noVariables * totalCells, MPI_DOUBLE, MPI_SUM,
                m_solver->m_rescalingCommGrComm, AT_, "varSliceLocal.begin()", "varSlice.begin()");

  MFloat blEdgeVValueLocal = F0;
  MFloat blEdgeVValueGlobal = F0;
  MPI_Allreduce(&blEdgeVValueLocal, &blEdgeVValueGlobal, 1, MPI_DOUBLE, MPI_SUM, m_solver->m_rescalingCommGrComm, AT_,
                "blEdgeVValueLocal", "blEdgeVValueGlobal");
}

/**
 * \brief Prescribe given profile BC
 *
 *  Precribes a profile from the restart file
 *  extrapolate pressure from computational domain
 */
template <MBool isRans>
void StructuredBndryCnd2D<isRans>::bc2600(MInt bcId) {
  MInt* start = m_physicalBCMap[bcId]->start1;
  MInt* end = m_physicalBCMap[bcId]->end1;

  switch(m_physicalBCMap[bcId]->face) {
    case 0: {
      for(MInt i = start[0]; i < end[0]; i++) {
        for(MInt j = start[1]; j < end[1]; j++) {
          MInt cellId = cellIndex(m_noGhostLayers - 1 - i, j);
          MInt cellIdadj = cellIndex(m_noGhostLayers - i, j);

          // extrapolate pressure to ghost cells
          m_cells->pvariables[PV->P][cellId] = m_cells->pvariables[PV->P][cellIdadj];
        }
      }
      break;
    }
    default: {
      mTerm(1, AT_, "Face not implemented");
      break;
    }
  }
}


// Works equally for conversion of fluid state into porous and vice versa
// Works for RANS and LES
template <MBool isRans>
void StructuredBndryCnd2D<isRans>::bc6002(MInt bcId) {
  TRACE();

  const MInt IJ[nDim] = {1, m_nCells[1]};
  MInt* start = m_physicalBCMap[bcId]->start1;
  MInt* end = m_physicalBCMap[bcId]->end1;

  constexpr MFloat eps = 1e-8;
  const MFloat gamma = m_solver->m_gamma;
  const MFloat gammaOverGammaMinusOne = gamma / (gamma - 1);

  if(m_physicalBCMap[bcId]->Nstar == -1) {
    // Here we find out the normal direction of the
    // boundary and the tangential direction.
    // This way we can make a general formulation of
    // the boundary condition
    const MInt face = m_physicalBCMap[bcId]->face;
    const MInt normalDir = face / 2;
    const MInt firstTangentialDir = (normalDir + 1) % nDim;
    const MInt normalDirStart = start[normalDir];
    const MInt firstTangentialStart = start[firstTangentialDir];
    const MInt firstTangentialEnd = end[firstTangentialDir];
    const MInt inc[nDim] = {IJ[normalDir], IJ[firstTangentialDir]};
    // determine indices for direction help
    const MInt n = (face % 2) * 2 - 1;                                    //-1,+1
    const MInt g[2] = {normalDirStart + (MInt)(0.5 - (0.5 * (MFloat)n)),  //+1,0
                       normalDirStart + (MInt)(0.5 + (0.5 * (MFloat)n))}; // 0,+1
    const MInt a1 = normalDirStart + (MInt)(0.5 - (1.5 * (MFloat)n));     //+2,-1

    // convention: index x is unknown and y is known
    for(MInt t1 = firstTangentialStart, pos = 0; t1 < firstTangentialEnd; t1++, ++pos) {
      const MInt cellIdG1 = g[0] * inc[0] + t1 * inc[1];
      const MFloat normalVec[nDim] = {m_cells->fq[FQ->NORMAL[0]][cellIdG1], m_cells->fq[FQ->NORMAL[1]][cellIdG1]};
      const MInt cellIdA1 = a1 * inc[0] + t1 * inc[1];
      const MFloat por_x = m_cells->fq[FQ->POROSITY][cellIdA1];
      for(MInt i = 0; i < 2 /*m_noGhostLayers*/; ++i) {
        const MInt cellIdG = g[i] * inc[0] + t1 * inc[1];

        const MFloat por_y = m_cells->fq[FQ->POROSITY][cellIdG];
        const MFloat rho_y = m_cells->pvariables[PV->RHO][cellIdG];
        const MFloat p_y = m_cells->pvariables[PV->P][cellIdG];
        const MFloat u_y = m_cells->pvariables[PV->U][cellIdG];
        const MFloat v_y = m_cells->pvariables[PV->V][cellIdG];
        const MFloat k_y = (isRans) ? m_cells->pvariables[PV->RANS_VAR[0]][cellIdG] : 0.0;

        // Compute normal and tangential velocity components
        const MFloat U_y = u_y * normalVec[0] + v_y * normalVec[1];
        const MFloat V_y[nDim] = {u_y - U_y * normalVec[0], v_y - U_y * normalVec[1]};

        //
        const MFloat auxconst1 = gammaOverGammaMinusOne * p_y / pow(rho_y, gamma);
        const MFloat auxconst2 = (por_y / por_x) * rho_y * k_y;
        const MFloat auxconst3 = -(gammaOverGammaMinusOne * p_y / rho_y + 0.5 * POW2(U_y) + k_y);
        const MFloat auxconst4 = 0.5 * POW2(por_y / por_x * rho_y * U_y);

        // TODO_SS labels:FV Later take a more sufficient initial value
        MFloat rho_x = rho_y;

        MFloat res = auxconst1 * pow(rho_x, gamma + 1) + auxconst2 * rho_x + auxconst3 * POW2(rho_x) + auxconst4;

        MInt testcounter = 0;
        while(res > eps) {
          // TODO_SS labels:FV,totest Check if it is divided by zero
          const MFloat dresdrho = (gamma + 1) * auxconst1 * pow(rho_x, gamma) + auxconst2 + 2 * auxconst3 * rho_x;
          const MFloat deltarho_x = -res / dresdrho;
          rho_x += deltarho_x;

          res = auxconst1 * pow(rho_x, gamma + 1) + auxconst2 * rho_x + auxconst3 * POW2(rho_x) + auxconst4;

          ++testcounter;
          if(testcounter > 100) break;
        }

        if(testcounter > 10) {
          cout << "Loop in BC6002 took " << testcounter << " steps!"
               << " domainId=" << m_solver->domainId() << " cellId: " << cellIdG << " (" << g[i] << "|" << t1
               << ") res=" << res << " (x|y)=" << m_cells->coordinates[0][cellIdG] << "|"
               << m_cells->coordinates[1][cellIdG] << por_y << " " << rho_y << " " << p_y << " " << u_y << " " << v_y
               << endl;
        }

        // Compute remaining variables
        const MFloat temp = por_y * rho_y / (por_x * rho_x);
        const MFloat p_x = p_y * pow(rho_x / rho_y, gamma);
        const MFloat U_x = temp * U_y;
        const MFloat u_x = U_x * normalVec[0] + V_y[0] /*=V_x*/;
        const MFloat v_x = U_x * normalVec[1] + V_y[1] /*=V_x*/;

        // Assign solution to pvariables
        m_cells->pvariables[PV->RHO][cellIdG] = rho_x;
        m_cells->pvariables[PV->P][cellIdG] = p_x;
        m_cells->pvariables[PV->U][cellIdG] = u_x;
        m_cells->pvariables[PV->V][cellIdG] = v_x;
        if(isRans) {
          m_cells->pvariables[PV->RANS_VAR[0]][cellIdG] *= temp;
          m_cells->pvariables[PV->RANS_VAR[1]][cellIdG] *= temp;
        }
      }
    }
  } else {
    const MInt singularId = m_physicalBCMap[bcId]->SingularId;
    const auto& singularity = m_solver->m_singularity[singularId];

    // Check
    const MInt* start2 = singularity.start;
    const MInt* end2 = singularity.end;
    for(MInt d = 0; d < nDim; ++d) {
      if(end[d] - start[d] != end2[d] - start2[d]) mTerm(1, "SCHEISSE");
    }

    //    cout << globalTimeStep<< ": SINGULARITY BC d=" << m_solver->domainId() << " start=" << start[0] << "|" <<
    //    start[1] << " start2=" << start2[0] << "|" << start2[1] << " end=" << end[0] << "|" << end[1] << " end2=" <<
    //    end2[0] << "|" << end2[1] << endl;

    for(MInt j = start[1], jj = start2[1]; j < end[1]; j++, jj++) {
      for(MInt i = start[0], ii = start2[0]; i < end[0]; i++, ii++) {
        const MInt cellIdG = cellIndex(i, j);
        const MInt cellIdA1 = cellIndex(ii, jj);
        const MFloat normalVec[nDim] = {m_cells->fq[FQ->NORMAL[0]][cellIdG], m_cells->fq[FQ->NORMAL[1]][cellIdG]};
        const MFloat por_x = m_cells->fq[FQ->POROSITY][cellIdA1];

        const MFloat por_y = m_cells->fq[FQ->POROSITY][cellIdG];
        const MFloat rho_y = m_cells->pvariables[PV->RHO][cellIdG];
        const MFloat p_y = m_cells->pvariables[PV->P][cellIdG];
        const MFloat u_y = m_cells->pvariables[PV->U][cellIdG];
        const MFloat v_y = m_cells->pvariables[PV->V][cellIdG];
        const MFloat k_y = (isRans) ? m_cells->pvariables[PV->RANS_VAR[0]][cellIdG] : 0.0;

        // Compute normal and tangential velocity components
        const MFloat U_y = u_y * normalVec[0] + v_y * normalVec[1];
        const MFloat V_y[nDim] = {u_y - U_y * normalVec[0], v_y - U_y * normalVec[1]};

        //
        const MFloat auxconst1 = gammaOverGammaMinusOne * p_y / pow(rho_y, gamma);
        const MFloat auxconst2 = (por_y / por_x) * rho_y * k_y;
        const MFloat auxconst3 = -(gammaOverGammaMinusOne * p_y / rho_y + 0.5 * POW2(U_y) + k_y);
        const MFloat auxconst4 = 0.5 * POW2(por_y / por_x * rho_y * U_y);

        // TODO_SS labels:FV Later take a more sufficient initial value
        MFloat rho_x = rho_y;

        MFloat res = auxconst1 * pow(rho_x, gamma + 1) + auxconst2 * rho_x + auxconst3 * POW2(rho_x) + auxconst4;

        MInt testcounter = 0;
        while(res > eps) {
          // TODO_SS labels:FV,totest Check if it is divided by zero
          const MFloat dresdrho = (gamma + 1) * auxconst1 * pow(rho_x, gamma) + auxconst2 + 2 * auxconst3 * rho_x;
          const MFloat deltarho_x = -res / dresdrho;
          rho_x += deltarho_x;

          res = auxconst1 * pow(rho_x, gamma + 1) + auxconst2 * rho_x + auxconst3 * POW2(rho_x) + auxconst4;

          ++testcounter;
          if(testcounter > 100) break;
        }

        if(testcounter > 10) {
          cout << "Loop in BC6002 took " << testcounter << " steps!"
               << " domainId=" << m_solver->domainId() << " cellId: " << cellIdG << " res=" << res
               << " (x|y)=" << m_cells->coordinates[0][cellIdG] << "|" << m_cells->coordinates[1][cellIdG] << por_y
               << " " << rho_y << " " << p_y << " " << u_y << " " << v_y << endl;
        }

        // Compute remaining variables
        const MFloat temp = por_y * rho_y / (por_x * rho_x);
        const MFloat p_x = p_y * pow(rho_x / rho_y, gamma);
        const MFloat U_x = temp * U_y;
        const MFloat u_x = U_x * normalVec[0] + V_y[0] /*=V_x*/;
        const MFloat v_x = U_x * normalVec[1] + V_y[1] /*=V_x*/;

        // Assign solution to pvariables
        m_cells->pvariables[PV->RHO][cellIdG] = rho_x;
        m_cells->pvariables[PV->P][cellIdG] = p_x;
        m_cells->pvariables[PV->U][cellIdG] = u_x;
        m_cells->pvariables[PV->V][cellIdG] = v_x;
        if(isRans) {
          m_cells->pvariables[PV->RANS_VAR[0]][cellIdG] *= temp;
          m_cells->pvariables[PV->RANS_VAR[1]][cellIdG] *= temp;
        }
      }
    }
  }
}


template <MBool isRans>
inline MFloat StructuredBndryCnd2D<isRans>::pressure(MInt cellId) {
  return m_cells->pvariables[PV->P][cellId];
}

template <MBool isRans>
inline MFloat StructuredBndryCnd2D<isRans>::temperature(MInt cellId) {
  const MFloat gamma = m_solver->m_gamma;
  MFloat t = gamma * m_cells->pvariables[PV->P][cellId] / m_cells->pvariables[PV->RHO][cellId];
  return t;
}


template class StructuredBndryCnd2D<true>;
template class StructuredBndryCnd2D<false>;
