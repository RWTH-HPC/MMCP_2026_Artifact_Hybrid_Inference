// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only


#include <cmath>
#include <fstream>
#include <iostream>

#include "COMM/mpiexchange.h"
#include "COMM/mpioverride.h"
#include "INCLUDE/maiatypes.h"
#include "IO/parallelio_hdf5.h"
#include "UTIL/kdtree.h"
#include "UTIL/pointbox.h"
#include "fvstructuredbndrycnd3d.h"
#include "fvstructuredsolver.h"
#include "fvstructuredsolver3d.h"

using namespace std;

template <MBool isRans>
constexpr MInt StructuredBndryCnd3D<isRans>::m_reverseCellIdDim[18];
template <MBool isRans>
constexpr MInt StructuredBndryCnd3D<isRans>::m_reverseCellIdGC[18];

/**
 * \brief Random number generator
 * \author Marian Albers
 *
 */
template <MBool isRans>
MFloat StructuredBndryCnd3D<isRans>::generate_rand() {
  return rand() / MFloat(RAND_MAX);
}

/**
 *  \brief Weighted random number generator
 *
 *  Generate weighted random numbers, for a distribution with
 *  a higher probability towards the higher bound pass a number < 1.0
 *  for a higher probability towards lower bound pass a number > 1.0
 *
 * \author Marian Albers
 *
 */
template <MBool isRans>
MFloat StructuredBndryCnd3D<isRans>::generate_rand_weighted() {
  return pow((rand() / MFloat(RAND_MAX)), m_solver->m_stgEddieDistribution);
}


/**
 * \brief Constructor of the 3D boundary conditions class
 * \author Pascal Meysonnat
 *
 */
template <MBool isRans>
StructuredBndryCnd3D<isRans>::StructuredBndryCnd3D(FvStructuredSolver<3>* solver, StructuredGrid<3>* grid)
  : StructuredBndryCnd<3>(solver, grid) {
  TRACE();

  const MLong oldAllocatedBytes = allocatedBytes();

  m_solver = static_cast<FvStructuredSolver3D*>(solver);

  // allocate virtual box
  if(m_solver->m_stgIsActive) {
    mAlloc(m_stgVbStart, 3, "m_stgVbStart", 0.0, AT_);
    mAlloc(m_stgVbEnd, 3, "m_stgVbEnd", 0.0, AT_);
  }

  // read rescaling bc properties
  if(m_solver->m_rescalingCommGrRoot >= 0) {
    /*! \property
      \page propertiesFVSTRCTRD
      \section rescalingBLT
      <code>MInt StructuredBndryCnd3D::m_rescalingBLT </code>\n
      default = <code> 1.0 </code>\n \n
      Delta0 thickness at the inflow to be rescaled to\n
      Possible values are:\n
      <ul>
      <li>Float > 0.0</li>
      </ul>
      Keywords: <i>RESCALING, STRUCTURED</i>
    */
    m_rescalingBLT = 5.95;
    m_rescalingBLT = Context::getSolverProperty<MFloat>("rescalingBLT", m_solverId, AT_);
  }

  // compute sponge coefficients for each cell
  if(m_solver->m_useSponge && m_solver->m_computeSpongeFactor) {
    RECORD_TIMER_START(m_solver->timer(Timers::BuildUpSponge));
    readAndDistributeSpongeCoordinates();
    RECORD_TIMER_STOP(m_solver->timer(Timers::BuildUpSponge));
  }

  printAllocatedMemory(oldAllocatedBytes, "StructuredBndryCnd3D", m_StructuredComm);
}


template <MBool isRans>
void StructuredBndryCnd3D<isRans>::readAndDistributeSpongeCoordinates() {
  MInt noSpongeInfo = m_solver->m_windowInfo->m_spongeInfoMap.size(); // contains the size of the maps
  MInt memSize = 0;

  // 1)
  // allocate the space for all the coordinates of the sponge in Scratch!
  // memory will not be needed later ==> Scratch!

  // 1.1) determine the storage size
  // determine the size and to store the whole data
  for(MInt i = 0; i < noSpongeInfo; ++i) {
    MInt size = 1;
    for(MInt dim = 0; dim < 3; ++dim) {
      size *= (m_solver->m_windowInfo->m_spongeInfoMap[i]->end1[dim]
               - m_solver->m_windowInfo->m_spongeInfoMap[i]->start1[dim] + 1);
    }
    memSize += size;
  }

  // 1.2) allocate the memory
  MFloatScratchSpace coordMem(3 * memSize, AT_, "spongeCoordinates");
  MFloatPointerScratchSpace spongeCoords(noSpongeInfo, AT_, "spongeCoordPointer");

  MInt totMemSize = 0;
  for(MInt i = 0; i < noSpongeInfo; ++i) {
    MInt size = 1;
    for(MInt dim = 0; dim < 3; ++dim) {
      size *= (m_solver->m_windowInfo->m_spongeInfoMap[i]->end1[dim]
               - m_solver->m_windowInfo->m_spongeInfoMap[i]->start1[dim] + 1);
    }
    spongeCoords[i] = &coordMem[totMemSize];
    totMemSize += 3 * size;
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
    for(MInt dim = 0; dim < 3; ++dim) {
      if(m_solver->m_windowInfo->m_spongeInfoMap[i]->end1[dim] - m_solver->m_windowInfo->m_spongeInfoMap[i]->start1[dim]
         == 0)
        continue;
      size *= (m_solver->m_windowInfo->m_spongeInfoMap[i]->end1[dim]
               - m_solver->m_windowInfo->m_spongeInfoMap[i]->start1[dim]);
    }
    cellmemSize += size;
  }
  // 2.2) allocate the space for all the coordinates in Scratch!
  // memory will not be needed later!
  MFloatScratchSpace coordCellMem(3 * cellmemSize, AT_, "spongeCellCoordinates");
  MFloatPointerScratchSpace spongeSurfCoords(noSpongeInfo, AT_, "spongeCellCoordPointer");

  MInt totCellMemSize = 0;
  for(MInt i = 0; i < noSpongeInfo; ++i) {
    MInt size = 1;
    for(MInt dim = 0; dim < 3; ++dim) {
      if(m_solver->m_windowInfo->m_spongeInfoMap[i]->end1[dim] - m_solver->m_windowInfo->m_spongeInfoMap[i]->start1[dim]
         == 0)
        continue;
      size *= (m_solver->m_windowInfo->m_spongeInfoMap[i]->end1[dim]
               - m_solver->m_windowInfo->m_spongeInfoMap[i]->start1[dim]);
    }
    spongeSurfCoords[i] = &coordCellMem[totCellMemSize];
    totCellMemSize += 3 * size;
  }

  // 3) read in the coordinates of the grid points
  // open file for reading the grid data
  MString gridFileName = m_grid->m_gridInputFileName;

  // unique identifier needed to associate grid and solution in case of restart
  // const char* aUID= new char[18];

  m_log << "Loading and broadcasting all sponge windows..." << endl;
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
  if(m_solver->domainId() == 0) {
    for(MInt i = 0; i < noSpongeInfo; ++i) {
      ParallelIo::size_type offset[3] = {0, 0, 0};
      ParallelIo::size_type size[3] = {0, 0, 0};
      memSize = 1;
      for(MInt dim = 2; dim >= 0; --dim) {
        size[dim] = (m_solver->m_windowInfo->m_spongeInfoMap[i]->end1[2 - dim]
                     - m_solver->m_windowInfo->m_spongeInfoMap[i]->start1[2 - dim] + 1);
        memSize *= size[dim];
        offset[dim] = m_solver->m_windowInfo->m_spongeInfoMap[i]->start1[2 - dim];
      }
      // read in the data if  processor zero else read nothing!!!
      // determine the Solver name
      MString bName = "block";
      stringstream number;
      number << m_solver->m_windowInfo->m_spongeInfoMap[i]->Id1;
      bName += number.str();
      pio.readArray(&spongeCoords[i][0], bName, "x", 3, offset, size);
      pio.readArray(&spongeCoords[i][memSize], bName, "y", 3, offset, size);
      pio.readArray(&spongeCoords[i][memSize * 2], bName, "z", 3, offset, size);
    }
  } else {
    for(MInt i = 0; i < noSpongeInfo; ++i) {
      ParallelIo::size_type offset[3] = {0, 0, 0};
      ParallelIo::size_type size[3] = {0, 0, 0};
      MString bName = "block";
      stringstream number;
      number << m_solver->m_windowInfo->m_spongeInfoMap[i]->Id1;
      bName += number.str();
      MFloat empty = 0;
      pio.readArray(&empty, bName, "x", 3, offset, size);
      pio.readArray(&empty, bName, "y", 3, offset, size);
      pio.readArray(&empty, bName, "z", 3, offset, size);
    }
  }

  // now broadcast the information to everyone!!!
  MPI_Bcast(&spongeCoords[0][0], totMemSize, MPI_DOUBLE, 0, m_StructuredComm, AT_, "spongeCoords[0][0]");
  m_log << "Loading and broadcasting all sponge windows... SUCCESSFUL!" << endl;
  m_log << "Computing sponge surface center coordinates..." << endl;

  // 4) computing the coordinates of surface center from corner points;

  for(MInt ii = 0; ii < noSpongeInfo; ++ii) {
    MInt label, size1, size2, count = 0;
    for(label = 0; label < 3; label++) { // 3== dimensions
      if(m_solver->m_windowInfo->m_spongeInfoMap[ii]->end1[label]
             - m_solver->m_windowInfo->m_spongeInfoMap[ii]->start1[label]
         == 0)
        break;
    }
    switch(label) {
      case 0: {
        size1 = (m_solver->m_windowInfo->m_spongeInfoMap[ii]->end1[1]
                 - m_solver->m_windowInfo->m_spongeInfoMap[ii]->start1[1] + 1);
        size2 = (m_solver->m_windowInfo->m_spongeInfoMap[ii]->end1[2]
                 - m_solver->m_windowInfo->m_spongeInfoMap[ii]->start1[2] + 1);
        for(MInt j = 0; j < size2 - 1; j++) {
          for(MInt i = 0; i < size1 - 1; i++) {
            MInt IJ = i + j * size1;
            MInt IPJ = i + (j + 1) * size1;
            MInt IJP = i + 1 + j * size1;
            MInt IPJP = i + 1 + (j + 1) * size1;
            spongeSurfCoords[ii][count] =
                0.25 * (spongeCoords[ii][IJ] + spongeCoords[ii][IPJ] + spongeCoords[ii][IJP] + spongeCoords[ii][IPJP]);
            spongeSurfCoords[ii][count + (size1 - 1) * (size2 - 1)] =
                0.25
                * (spongeCoords[ii][IJ + size1 * size2] + spongeCoords[ii][IPJ + size1 * size2]
                   + spongeCoords[ii][IJP + size1 * size2] + spongeCoords[ii][IPJP + size1 * size2]);
            spongeSurfCoords[ii][count + 2 * (size1 - 1) * (size2 - 1)] =
                0.25
                * (spongeCoords[ii][IJ + 2 * size1 * size2] + spongeCoords[ii][IPJ + 2 * size1 * size2]
                   + spongeCoords[ii][IJP + 2 * size1 * size2] + spongeCoords[ii][IPJP + 2 * size1 * size2]);
            count++;
          }
        }
        break;
      }
      case 1: {
        size1 = (m_solver->m_windowInfo->m_spongeInfoMap[ii]->end1[0]
                 - m_solver->m_windowInfo->m_spongeInfoMap[ii]->start1[0] + 1);
        size2 = (m_solver->m_windowInfo->m_spongeInfoMap[ii]->end1[2]
                 - m_solver->m_windowInfo->m_spongeInfoMap[ii]->start1[2] + 1);
        for(MInt j = 0; j < size2 - 1; j++) {
          for(MInt i = 0; i < size1 - 1; i++) {
            MInt IJ = i + j * size1;
            MInt IPJ = i + (j + 1) * size1;
            MInt IJP = i + 1 + j * size1;
            MInt IPJP = i + 1 + (j + 1) * size1;
            spongeSurfCoords[ii][count] =
                0.25 * (spongeCoords[ii][IJ] + spongeCoords[ii][IPJ] + spongeCoords[ii][IJP] + spongeCoords[ii][IPJP]);
            spongeSurfCoords[ii][count + (size1 - 1) * (size2 - 1)] =
                0.25
                * (spongeCoords[ii][IJ + size1 * size2] + spongeCoords[ii][IPJ + size1 * size2]
                   + spongeCoords[ii][IJP + size1 * size2] + spongeCoords[ii][IPJP + size1 * size2]);
            spongeSurfCoords[ii][count + 2 * (size1 - 1) * (size2 - 1)] =
                0.25
                * (spongeCoords[ii][IJ + 2 * size1 * size2] + spongeCoords[ii][IPJ + 2 * size1 * size2]
                   + spongeCoords[ii][IJP + 2 * size1 * size2] + spongeCoords[ii][IPJP + 2 * size1 * size2]);
            count++;
          }
        }
        break;
      }
      case 2: {
        size1 = (m_solver->m_windowInfo->m_spongeInfoMap[ii]->end1[0]
                 - m_solver->m_windowInfo->m_spongeInfoMap[ii]->start1[0] + 1);
        size2 = (m_solver->m_windowInfo->m_spongeInfoMap[ii]->end1[1]
                 - m_solver->m_windowInfo->m_spongeInfoMap[ii]->start1[1] + 1);
        for(MInt j = 0; j < size2 - 1; j++) {
          for(MInt i = 0; i < size1 - 1; i++) {
            MInt IJ = i + j * size1;
            MInt IPJ = i + (j + 1) * size1;
            MInt IJP = i + 1 + j * size1;
            MInt IPJP = i + 1 + (j + 1) * size1;
            spongeSurfCoords[ii][count] =
                0.25 * (spongeCoords[ii][IJ] + spongeCoords[ii][IPJ] + spongeCoords[ii][IJP] + spongeCoords[ii][IPJP]);
            spongeSurfCoords[ii][count + (size1 - 1) * (size2 - 1)] =
                0.25
                * (spongeCoords[ii][IJ + size1 * size2] + spongeCoords[ii][IPJ + size1 * size2]
                   + spongeCoords[ii][IJP + size1 * size2] + spongeCoords[ii][IPJP + size1 * size2]);
            spongeSurfCoords[ii][count + 2 * (size1 - 1) * (size2 - 1)] =
                0.25
                * (spongeCoords[ii][IJ + 2 * size1 * size2] + spongeCoords[ii][IPJ + 2 * size1 * size2]
                   + spongeCoords[ii][IJP + 2 * size1 * size2] + spongeCoords[ii][IPJP + 2 * size1 * size2]);
            count++;
          }
        }
        break;
      }
      default:
        mTerm(1, AT_, "sponge direction is messed up");
    }
  }

  m_log << "Computing sponge surface center coordinates... SUCCESSFUL!" << endl;
  m_log << "Determining shortest distance and sponge factor for each cell..." << endl;

  // cout << "seraching for the nearest points(building sigma sponge)" << endl;
  // build a k-d-tree for a quick search:
  // 1) rearrange the coordinates into points;
  const MInt spongeTotalWorkload = noSpongeInfo * m_noCells;
  const MInt spongeWorkloadPercentage = (MInt)(spongeTotalWorkload / 10);

  for(MInt i = 0; i < noSpongeInfo; ++i) {
    MInt noPoints = 1;
    for(MInt dim = 0; dim < 3; ++dim) {
      if(m_solver->m_windowInfo->m_spongeInfoMap[i]->end1[dim] - m_solver->m_windowInfo->m_spongeInfoMap[i]->start1[dim]
         == 0)
        continue;
      noPoints *= (m_solver->m_windowInfo->m_spongeInfoMap[i]->end1[dim]
                   - m_solver->m_windowInfo->m_spongeInfoMap[i]->start1[dim]);
    }
    // build up the points (surface centres)
    vector<Point<3>> pts;
    for(MInt j = 0; j < noPoints; ++j) {
      Point<3> a(spongeSurfCoords[i][j], spongeSurfCoords[i][j + noPoints], spongeSurfCoords[i][j + 2 * noPoints]);
      pts.push_back(a);
    }

    // build up the tree
    KDtree<3> tree(pts);
    MFloat distance = -1.0;
    // go through all the cells an determine the closest distance
    MFloat spongeThickness = m_solver->m_windowInfo->m_spongeInfoMap[i]->spongeThickness;
    for(MInt id = 0; id < m_noCells; ++id) {
      if(m_solver->domainId() == 0) {
        if((i * m_noCells + id) % spongeWorkloadPercentage == 0) {
          MInt progress = ceil((MFloat)(i * m_noCells + id) / (MFloat)spongeTotalWorkload * 100.0);
          cout << "Sponge computation - " << progress << " percent" << endl;
        }
      }
      distance = -1.1111111111111111; // to check
      Point<3> pt(m_cells->coordinates[0][id], m_cells->coordinates[1][id], m_cells->coordinates[2][id]);
      (void)tree.nearest(pt, distance);
      if(distance <= spongeThickness) {
        MFloat spongeFactor =
            m_solver->m_windowInfo->m_spongeInfoMap[i]->sigma
            * pow((spongeThickness - distance) / spongeThickness, m_solver->m_windowInfo->m_spongeInfoMap[i]->beta);
        m_cells->fq[FQ->SPONGE_FACTOR][id] = mMax(m_cells->fq[FQ->SPONGE_FACTOR][id], spongeFactor);
      }
    }
  }

  m_log << "Determining shortest distance and sponge factor for each cell... SUCCESSFUL!" << endl;
}

template <MBool isRans>
void StructuredBndryCnd3D<isRans>::updateSpongeLayer() {
  const MFloat gammaMinusOne = m_solver->m_gamma - 1.0;
  switch(m_spongeLayerType) {
    case 1: {
      MFloat deltaRhoE = F0, deltaRho = F0;
      for(MInt k = m_noGhostLayers; k < m_nCells[0] - m_noGhostLayers; k++) {
        for(MInt j = m_noGhostLayers; j < m_nCells[1] - m_noGhostLayers; j++) {
          for(MInt i = m_noGhostLayers; i < m_nCells[2] - m_noGhostLayers; i++) {
            // compute the forcing term
            // for the pressure or engery and the velocity
            MInt cellId = cellIndex(i, j, k);

            const MFloat rhoE =
                m_cells->pvariables[PV->P][cellId] / gammaMinusOne
                + F1B2 * m_cells->pvariables[PV->RHO][cellId]
                      * (POW2(m_cells->pvariables[PV->U][cellId]) + POW2(m_cells->pvariables[PV->V][cellId])
                         + POW2(m_cells->pvariables[PV->W][cellId]));

            deltaRhoE = rhoE - CV->rhoEInfinity;
            deltaRho = m_cells->pvariables[PV->RHO][cellId] - (CV->rhoInfinity * m_targetDensityFactor);

            m_cells->rightHandSide[CV->RHO_E][cellId] -=
                m_cells->fq[FQ->SPONGE_FACTOR][cellId] * deltaRhoE
                * m_cells->cellJac[cellId]; // deltaP * m_cells->cellJac[cellId];
            m_cells->rightHandSide[CV->RHO][cellId] -=
                m_cells->fq[FQ->SPONGE_FACTOR][cellId] * deltaRho * m_cells->cellJac[cellId];
          }
        }
      }
      break;
    }
    case 2: {
      // damp to values in the FQ field (set at startup, from RANS etc.)
      MFloat deltaRhoE = F0, deltaRho = F0;
      for(MInt k = m_noGhostLayers; k < m_nCells[0] - m_noGhostLayers; k++) {
        for(MInt j = m_noGhostLayers; j < m_nCells[1] - m_noGhostLayers; j++) {
          for(MInt i = m_noGhostLayers; i < m_nCells[2] - m_noGhostLayers; i++) {
            MInt cellId = cellIndex(i, j, k);

            const MFloat rhoE =
                m_cells->pvariables[PV->P][cellId] / gammaMinusOne
                + F1B2 * m_cells->pvariables[PV->RHO][cellId]
                      * (POW2(m_cells->pvariables[PV->U][cellId]) + POW2(m_cells->pvariables[PV->V][cellId])
                         + POW2(m_cells->pvariables[PV->W][cellId]));

            deltaRhoE = rhoE - m_cells->fq[FQ->SPONGE_RHO_E][cellId];
            deltaRho =
                m_cells->pvariables[PV->RHO][cellId] - (m_cells->fq[FQ->SPONGE_RHO][cellId] * m_targetDensityFactor);

            m_cells->rightHandSide[CV->RHO_E][cellId] -=
                m_cells->fq[FQ->SPONGE_FACTOR][cellId] * deltaRhoE * m_cells->cellJac[cellId];
            m_cells->rightHandSide[CV->RHO][cellId] -=
                m_cells->fq[FQ->SPONGE_FACTOR][cellId] * deltaRho * m_cells->cellJac[cellId];
          }
        }
      }
      break;
    }
    case 3: {
      // damp to rhoInfinity and pInfinity
      const MFloat FgammaMinusOne = F1 / (m_solver->m_gamma - F1);
      for(MInt k = m_noGhostLayers; k < m_nCells[0] - m_noGhostLayers; k++) {
        for(MInt j = m_noGhostLayers; j < m_nCells[1] - m_noGhostLayers; j++) {
          for(MInt i = m_noGhostLayers; i < m_nCells[2] - m_noGhostLayers; i++) {
            // compute the forcing term
            // for the pressure or engery and the velocity
            const MInt cellId = cellIndex(i, j, k);
            const MFloat deltaP = (m_cells->pvariables[PV->P][cellId] - PV->PInfinity) * FgammaMinusOne;
            const MFloat deltaRho = m_cells->pvariables[PV->RHO][cellId] - (CV->rhoInfinity * m_targetDensityFactor);

            m_cells->rightHandSide[CV->RHO_E][cellId] -=
                m_cells->fq[FQ->SPONGE_FACTOR][cellId] * deltaP * m_cells->cellJac[cellId];
            m_cells->rightHandSide[CV->RHO][cellId] -=
                m_cells->fq[FQ->SPONGE_FACTOR][cellId] * deltaRho * m_cells->cellJac[cellId];
          }
        }
      }
      break;
    }
    case 4: {
      // damp to rho from FQ and pInfinity
      const MFloat FgammaMinusOne = F1 / (m_solver->m_gamma - F1);
      for(MInt k = m_noGhostLayers; k < m_nCells[0] - m_noGhostLayers; k++) {
        for(MInt j = m_noGhostLayers; j < m_nCells[1] - m_noGhostLayers; j++) {
          for(MInt i = m_noGhostLayers; i < m_nCells[2] - m_noGhostLayers; i++) {
            const MInt cellId = cellIndex(i, j, k);
            const MFloat deltaP = (m_cells->pvariables[PV->P][cellId] - PV->PInfinity) * FgammaMinusOne;
            const MFloat deltaRho =
                m_cells->pvariables[PV->RHO][cellId] - (m_cells->fq[FQ->SPONGE_RHO][cellId] * m_targetDensityFactor);

            m_cells->rightHandSide[CV->RHO_E][cellId] -=
                m_cells->fq[FQ->SPONGE_FACTOR][cellId] * deltaP * m_cells->cellJac[cellId];
            m_cells->rightHandSide[CV->RHO][cellId] -=
                m_cells->fq[FQ->SPONGE_FACTOR][cellId] * deltaRho * m_cells->cellJac[cellId];
          }
        }
      }
      break;
    }
    case 5: {
      // damp to p from FSC
      const MFloat FgammaMinusOne = F1 / (m_solver->m_gamma - F1);
      for(MInt k = m_noGhostLayers; k < m_nCells[0] - m_noGhostLayers; k++) {
        for(MInt j = m_noGhostLayers; j < m_nCells[1] - m_noGhostLayers; j++) {
          for(MInt i = m_noGhostLayers; i < m_nCells[2] - m_noGhostLayers; i++) {
            const MInt cellId = cellIndex(i, j, k);
            const MFloat deltaP =
                (m_cells->pvariables[PV->P][cellId] - m_solver->getFscPressure(cellId)) * FgammaMinusOne;
            m_cells->rightHandSide[CV->RHO_E][cellId] -=
                m_cells->fq[FQ->SPONGE_FACTOR][cellId] * deltaP * m_cells->cellJac[cellId];
          }
        }
      }
      break;
    }
    case 6: {
      // damp to PInfinity
      const MFloat FgammaMinusOne = F1 / (m_solver->m_gamma - F1);
      for(MInt k = m_noGhostLayers; k < m_nCells[0] - m_noGhostLayers; k++) {
        for(MInt j = m_noGhostLayers; j < m_nCells[1] - m_noGhostLayers; j++) {
          for(MInt i = m_noGhostLayers; i < m_nCells[2] - m_noGhostLayers; i++) {
            const MInt cellId = cellIndex(i, j, k);
            const MFloat deltaP = (m_cells->pvariables[PV->P][cellId] - PV->PInfinity) * FgammaMinusOne;
            m_cells->rightHandSide[CV->RHO_E][cellId] -=
                m_cells->fq[FQ->SPONGE_FACTOR][cellId] * deltaP * m_cells->cellJac[cellId];
          }
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
void StructuredBndryCnd3D<isRans>::computeWallDistances() {
  MInt noWallDistInfo = m_solver->m_windowInfo->m_wallDistInfoMap.size(); // contains the size of the maps
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
    for(MInt dim = 0; dim < 3; ++dim) {
      size *= (m_solver->m_windowInfo->m_wallDistInfoMap[i]->end1[dim]
               - m_solver->m_windowInfo->m_wallDistInfoMap[i]->start1[dim] + 1);
    }
    memSize += size;
  }

  // 1.2) allocate the memory
  MFloatScratchSpace coordMem(3 * memSize, AT_, "wallCoordinates");
  MFloatPointerScratchSpace wallDistCoords(noWallDistInfo, AT_, "wallCoordsPointer");

  MInt totMemSize = 0;
  for(MInt i = 0; i < noWallDistInfo; ++i) {
    MInt size = 1;
    for(MInt dim = 0; dim < 3; ++dim) {
      size *= (m_solver->m_windowInfo->m_wallDistInfoMap[i]->end1[dim]
               - m_solver->m_windowInfo->m_wallDistInfoMap[i]->start1[dim] + 1);
    }
    wallDistCoords[i] = &coordMem[totMemSize];
    totMemSize += 3 * size;
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
    for(MInt dim = 0; dim < 3; ++dim) {
      if(m_solver->m_windowInfo->m_wallDistInfoMap[i]->end1[dim]
             - m_solver->m_windowInfo->m_wallDistInfoMap[i]->start1[dim]
         == 0)
        continue;
      size *= (m_solver->m_windowInfo->m_wallDistInfoMap[i]->end1[dim]
               - m_solver->m_windowInfo->m_wallDistInfoMap[i]->start1[dim]);
    }
    cellmemSize += size;
  }
  // 2.2) allocate the space for all the coordinates in Scratch!
  // memory will not be needed later!
  MFloatScratchSpace coordCellMem(3 * cellmemSize, AT_, "wallDistCellCoordinates");
  MFloatPointerScratchSpace wallDistSurfCoords(noWallDistInfo, AT_, "wallDistCellCoordPointer");

  MInt totCellMemSize = 0;
  for(MInt i = 0; i < noWallDistInfo; ++i) {
    MInt size = 1;
    for(MInt dim = 0; dim < 3; ++dim) {
      if(m_solver->m_windowInfo->m_wallDistInfoMap[i]->end1[dim]
             - m_solver->m_windowInfo->m_wallDistInfoMap[i]->start1[dim]
         == 0)
        continue;
      size *= (m_solver->m_windowInfo->m_wallDistInfoMap[i]->end1[dim]
               - m_solver->m_windowInfo->m_wallDistInfoMap[i]->start1[dim]);
    }
    wallDistSurfCoords[i] = &coordCellMem[totCellMemSize];
    totCellMemSize += 3 * size;
  }


  // 3) read in the coordinates of the grid points
  // open file for reading the grid data
  MString gridFileName = m_grid->m_gridInputFileName;
  // open file and read number of solvers and uid

  // read the data in and distribute ist

  // the split of reading and distributing is done on purpose!
  // this is because we then know what the library is doing
  // and not what the io_library such as hdf or netcdf should do
  // but does not
  // a good library should handle that directly! TEST IT
  memSize = 1;
  // read in the data if  processor zero else read nothing!!!
  if(m_solver->domainId() == 0) {
    ParallelIoHdf5 pio(gridFileName, maia::parallel_io::PIO_READ, MPI_COMM_SELF);
    for(MInt i = 0; i < noWallDistInfo; ++i) {
      ParallelIo::size_type offset[3] = {0, 0, 0};
      ParallelIo::size_type size[3] = {0, 0, 0};
      memSize = 1;
      for(MInt dim = 2; dim >= 0; --dim) {
        size[dim] = (m_solver->m_windowInfo->m_wallDistInfoMap[i]->end1[2 - dim]
                     - m_solver->m_windowInfo->m_wallDistInfoMap[i]->start1[2 - dim] + 1);
        memSize *= size[dim];
        offset[dim] = m_solver->m_windowInfo->m_wallDistInfoMap[i]->start1[2 - dim];
      }
      // read in the data if  processor zero else read nothing!!!
      // determine the Solver name
      MString bName = "block";
      stringstream number;
      number << m_solver->m_windowInfo->m_wallDistInfoMap[i]->Id1;
      bName += number.str();
      pio.readArray(&wallDistCoords[i][0], bName, "x", 3, offset, size);
      pio.readArray(&wallDistCoords[i][memSize], bName, "y", 3, offset, size);
      pio.readArray(&wallDistCoords[i][memSize * 2], bName, "z", 3, offset, size);
    }
  }

  MPI_Bcast(&wallDistCoords[0][0], totMemSize, MPI_DOUBLE, 0, m_solver->m_StructuredComm, AT_, "wallDistCoords[0][0]");

  if(!m_solver->m_rans) {
    return;
  }


  // 4) computing the coordinates of surface center from corner points;
  for(MInt ii = 0; ii < noWallDistInfo; ++ii) {
    MInt label, size1, size2, count = 0;
    for(label = 0; label < 3; label++) { // 3== dimensions
      if(m_solver->m_windowInfo->m_wallDistInfoMap[ii]->end1[label]
             - m_solver->m_windowInfo->m_wallDistInfoMap[ii]->start1[label]
         == 0)
        break;
    }
    switch(label) {
      case 0: {
        size1 = (m_solver->m_windowInfo->m_wallDistInfoMap[ii]->end1[1]
                 - m_solver->m_windowInfo->m_wallDistInfoMap[ii]->start1[1] + 1);
        size2 = (m_solver->m_windowInfo->m_wallDistInfoMap[ii]->end1[2]
                 - m_solver->m_windowInfo->m_wallDistInfoMap[ii]->start1[2] + 1);
        for(MInt j = 0; j < size2 - 1; j++) {
          for(MInt i = 0; i < size1 - 1; i++) {
            MInt IJ = i + j * size1;
            MInt IPJ = i + (j + 1) * size1;
            MInt IJP = i + 1 + j * size1;
            MInt IPJP = i + 1 + (j + 1) * size1;
            wallDistSurfCoords[ii][count] = 0.25
                                            * (wallDistCoords[ii][IJ] + wallDistCoords[ii][IPJ]
                                               + wallDistCoords[ii][IJP] + wallDistCoords[ii][IPJP]);
            wallDistSurfCoords[ii][count + (size1 - 1) * (size2 - 1)] =
                0.25
                * (wallDistCoords[ii][IJ + size1 * size2] + wallDistCoords[ii][IPJ + size1 * size2]
                   + wallDistCoords[ii][IJP + size1 * size2] + wallDistCoords[ii][IPJP + size1 * size2]);
            wallDistSurfCoords[ii][count + 2 * (size1 - 1) * (size2 - 1)] =
                0.25
                * (wallDistCoords[ii][IJ + 2 * size1 * size2] + wallDistCoords[ii][IPJ + 2 * size1 * size2]
                   + wallDistCoords[ii][IJP + 2 * size1 * size2] + wallDistCoords[ii][IPJP + 2 * size1 * size2]);
            count++;
          }
        }
        break;
      }
      case 1: {
        size1 = (m_solver->m_windowInfo->m_wallDistInfoMap[ii]->end1[0]
                 - m_solver->m_windowInfo->m_wallDistInfoMap[ii]->start1[0] + 1);
        size2 = (m_solver->m_windowInfo->m_wallDistInfoMap[ii]->end1[2]
                 - m_solver->m_windowInfo->m_wallDistInfoMap[ii]->start1[2] + 1);
        for(MInt j = 0; j < size2 - 1; j++) {
          for(MInt i = 0; i < size1 - 1; i++) {
            MInt IJ = i + j * size1;
            MInt IPJ = i + (j + 1) * size1;
            MInt IJP = i + 1 + j * size1;
            MInt IPJP = i + 1 + (j + 1) * size1;
            wallDistSurfCoords[ii][count] = 0.25
                                            * (wallDistCoords[ii][IJ] + wallDistCoords[ii][IPJ]
                                               + wallDistCoords[ii][IJP] + wallDistCoords[ii][IPJP]);
            wallDistSurfCoords[ii][count + (size1 - 1) * (size2 - 1)] =
                0.25
                * (wallDistCoords[ii][IJ + size1 * size2] + wallDistCoords[ii][IPJ + size1 * size2]
                   + wallDistCoords[ii][IJP + size1 * size2] + wallDistCoords[ii][IPJP + size1 * size2]);
            wallDistSurfCoords[ii][count + 2 * (size1 - 1) * (size2 - 1)] =
                0.25
                * (wallDistCoords[ii][IJ + 2 * size1 * size2] + wallDistCoords[ii][IPJ + 2 * size1 * size2]
                   + wallDistCoords[ii][IJP + 2 * size1 * size2] + wallDistCoords[ii][IPJP + 2 * size1 * size2]);
            count++;
          }
        }
        break;
      }
      case 2: {
        size1 = (m_solver->m_windowInfo->m_wallDistInfoMap[ii]->end1[0]
                 - m_solver->m_windowInfo->m_wallDistInfoMap[ii]->start1[0] + 1);
        size2 = (m_solver->m_windowInfo->m_wallDistInfoMap[ii]->end1[1]
                 - m_solver->m_windowInfo->m_wallDistInfoMap[ii]->start1[1] + 1);
        for(MInt j = 0; j < size2 - 1; j++) {
          for(MInt i = 0; i < size1 - 1; i++) {
            MInt IJ = i + j * size1;
            MInt IPJ = i + (j + 1) * size1;
            MInt IJP = i + 1 + j * size1;
            MInt IPJP = i + 1 + (j + 1) * size1;
            wallDistSurfCoords[ii][count] = 0.25
                                            * (wallDistCoords[ii][IJ] + wallDistCoords[ii][IPJ]
                                               + wallDistCoords[ii][IJP] + wallDistCoords[ii][IPJP]);
            wallDistSurfCoords[ii][count + (size1 - 1) * (size2 - 1)] =
                0.25
                * (wallDistCoords[ii][IJ + size1 * size2] + wallDistCoords[ii][IPJ + size1 * size2]
                   + wallDistCoords[ii][IJP + size1 * size2] + wallDistCoords[ii][IPJP + size1 * size2]);
            wallDistSurfCoords[ii][count + 2 * (size1 - 1) * (size2 - 1)] =
                0.25
                * (wallDistCoords[ii][IJ + 2 * size1 * size2] + wallDistCoords[ii][IPJ + 2 * size1 * size2]
                   + wallDistCoords[ii][IJP + 2 * size1 * size2] + wallDistCoords[ii][IPJP + 2 * size1 * size2]);
            count++;
          }
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
    for(MInt dim = 0; dim < 3; ++dim) {
      if(m_solver->m_windowInfo->m_wallDistInfoMap[i]->end1[dim]
             - m_solver->m_windowInfo->m_wallDistInfoMap[i]->start1[dim]
         == 0)
        continue;
      noPoints *= (m_solver->m_windowInfo->m_wallDistInfoMap[i]->end1[dim]
                   - m_solver->m_windowInfo->m_wallDistInfoMap[i]->start1[dim]);
    }
    // build up the points (surface centres)
    vector<Point<3>> pts;
    for(MInt j = 0; j < noPoints; ++j) {
      Point<3> a(wallDistSurfCoords[i][j], wallDistSurfCoords[i][j + noPoints],
                 wallDistSurfCoords[i][j + 2 * noPoints]);
      pts.push_back(a);
    }
    // build up the tree
    KDtree<3> tree(pts);
    MFloat distance = -1.0;

    // go through all the cells an determine the closest distance
    for(MInt id = 0; id < m_noCells; ++id) {
      distance = -1.1111111111111111; // to check
      Point<3> pt(m_cells->coordinates[0][id], m_cells->coordinates[1][id], m_cells->coordinates[2][id]);
      (void)tree.nearest(pt, distance);

      // take minimum because another bc1000 might be further away than the current but would overwrite the actually
      // closer distance
      m_cells->fq[FQ->WALLDISTANCE][id] = mMin(m_cells->fq[FQ->WALLDISTANCE][id], distance);
    }
  }

  // correct the wall distance for the effective solution
  if(m_solver->m_bc2601IsActive) {
    cout << "Correcting wall distance with gammEpsilon: " << m_solver->m_bc2601GammaEpsilon << endl;
    for(MInt id = 0; id < m_noCells; ++id) {
      m_cells->fq[FQ->WALLDISTANCE][id] += m_solver->m_bc2601GammaEpsilon;
    }
  }

  m_log << "Wall Distance Computation SUCESSFUL: Saved minimum distance to next wall for all cells " << endl;
}
//<marian


template <MBool isRans>
StructuredBndryCnd3D<isRans>::~StructuredBndryCnd3D() {}

template <MBool isRans>
inline MInt StructuredBndryCnd3D<isRans>::cellIndex(MInt i, MInt j, MInt k) {
  return i + (j + k * m_nCells[1]) * m_nCells[2];
}

template <MBool isRans>
inline MInt StructuredBndryCnd3D<isRans>::cellIndexBC(MInt i, MInt j, MInt k) {
  return i + (j + k * m_solver->m_stgBoxSize[1]) * m_solver->m_stgBoxSize[2];
}

// function to correct the index values in the map for the different boundary conditions
template <MBool isRans>
void StructuredBndryCnd3D<isRans>::correctBndryCndIndices() {
  // cout << "in correctBndryCndIndices " << endl;
  // in correcting cell Information
  for(MInt bcId = 0; bcId < m_noBndryCndIds; bcId++) {
    (this->*initBndryCndHandler[bcId])(bcId);
  }
}

template <MBool isRans>
inline MFloat StructuredBndryCnd3D<isRans>::dist(MFloat* a, MFloat* b) {
  MFloat dist1 = F0;
  for(MInt dim = 0; dim < 3; dim++) {
    dist1 += POW2(a[dim * m_noCells] - b[dim * m_noCells]);
  }
  return sqrt(dist1);
}

template <MBool isRans>
void StructuredBndryCnd3D<isRans>::initBc1000(MInt bcId) {
  (void)bcId;
}

template <MBool isRans>
void StructuredBndryCnd3D<isRans>::initBc1003(MInt bcId) {
  (void)bcId;

  /*! \property
    \page propertiesFVSTRCTRD
    \section isothermalWallTemperature
    <code>MInt StructuredBndryCnd3D::m_isothermalWallTemperature </code>\n
    default = <code> 1.0 </code>\n \n
    Isothermal wall temperature as a factor of T8\n
    Possible values are:\n
    <ul>
    <li>Float > 0.0</li>
    </ul>
    Keywords: <i>ISOTHERMAL, WALL, BC, STRUCTURED</i>
  */
  m_isothermalWallTemperature = F1;
  if(Context::propertyExists("isothermalWallTemperature", m_solverId)) {
    m_isothermalWallTemperature =
        Context::getSolverProperty<MFloat>("isothermalWallTemperature", m_solverId, AT_, &m_isothermalWallTemperature);
  }
}

template <MBool isRans>
void StructuredBndryCnd3D<isRans>::initBc1004(MInt bcId) { // moving adiabatic wall
  (void)bcId;
}

template <MBool isRans>
void StructuredBndryCnd3D<isRans>::initBc1006(MInt bcId) { // moving adiabatic wall
  (void)bcId;

  m_isothermalWallTemperature = F1;
  if(Context::propertyExists("isothermalWallTemperature", m_solverId)) {
    m_isothermalWallTemperature =
        Context::getSolverProperty<MFloat>("isothermalWallTemperature", m_solverId, AT_, &m_isothermalWallTemperature);
  }
}

template <MBool isRans>
void StructuredBndryCnd3D<isRans>::initBc1007(MInt bcId) { // oscillating non-moving wall
  (void)bcId;

  m_solver->m_waveBeginTransition = 0.0;
  m_solver->m_waveEndTransition = 0.0;
  m_solver->m_waveAmplitude = 0.0;
  m_solver->m_waveAmplitudePlus = 0.0;
  m_solver->m_waveTimePlus = 0.0;
  m_solver->m_waveTime = 0.0;

  // time needs to be constant for traveling wave
  m_solver->m_constantTimeStep = true;
  m_solver->m_waveAmplitudePlus =
      Context::getSolverProperty<MFloat>("waveAmplitudePlus", m_solverId, AT_, &m_solver->m_waveAmplitudePlus);
  m_solver->m_waveTimePlus =
      Context::getSolverProperty<MFloat>("waveTimePlus", m_solverId, AT_, &m_solver->m_waveTimePlus);
  m_solver->m_waveBeginTransition =
      Context::getSolverProperty<MFloat>("waveBeginTransition", m_solverId, AT_, &m_solver->m_waveBeginTransition);
  m_solver->m_waveEndTransition =
      Context::getSolverProperty<MFloat>("waveEndTransition", m_solverId, AT_, &m_solver->m_waveEndTransition);

  // compute Wave parameters
  const MFloat cf = 0.024 * pow(m_solver->m_Re, -F1B4);
  const MFloat uTau = sqrt(cf * F1B2) * PV->UInfinity;
  const MFloat mu8 = SUTHERLANDLAW(PV->TInfinity);
  m_solver->m_waveAmplitude = m_solver->m_waveAmplitudePlus * uTau;
  m_solver->m_waveTime = m_solver->m_waveTimePlus * mu8 / (POW2(uTau) * m_solver->m_Re0);

  cout << "Oscillation speed amplitude: " << m_solver->m_waveAmplitude << " time: " << m_solver->m_waveTime << endl;
}


/** \fn StructuredBndryCnd3D<isRans>::initBc2004(MInt bcId)
 * \brief Initialize with standard pressure extrapolation
 *        at inflow or prescribe p_inf at outflow
 *
 */
template <MBool isRans>
void StructuredBndryCnd3D<isRans>::initBc2004(MInt bcId) {
  // call simple in/outflow bc to initialize ghost-cells
  // because bc2024 uses values in the ghost-cells
  bc2003(bcId);
}

/** \fn void StructuredBndryCnd3D<isRans>::initBc2009(MInt bcId)
 * \brief Characteristic boundary condition supersonic after shock
 *
 * \authors Simon Loosen, Pascal Meysonnat
 *
 */
template <MBool isRans>
void StructuredBndryCnd3D<isRans>::initBc2009(MInt bcId) {
  // bcId is not needed
  (void)bcId;

  /*! \property
    \page propertiesFVSTRCTRD
    \section shockAngle
    <code>MInt StructuredBndryCnd3D::m_sigma </code>\n
    default = <code> 0 </code>\n \n
    Angle of the shock to be introduced in BC 2009i\n
    Possible values are:\n
    <ul>
    <li>Float > positive float values </li>
    </ul>
    Keywords: <i>ISOTHERMAL, WALL, BC, STRUCTURED</i>
  */
  m_sigma = Context::getSolverProperty<MFloat>("shockAngle", m_solverId, AT_);
  m_sigma = (m_sigma / 180.0) * PI;
}

template <MBool isRans>
void StructuredBndryCnd3D<isRans>::initBc2222(MInt bcId) {
  bc2003(bcId);
}


template <MBool isRans>
void StructuredBndryCnd3D<isRans>::initBc2097(MInt bcId) {
  (void)bcId;
  //###############
  // compute surface of the plenum automatically
  //###############

  MInt Id = m_plenumSurfaceIndexMap[bcId]; // shift of counter form bcId-> to local counter
  const MInt* startface = &m_solver->m_windowInfo->plenumInletSurfaceIndices[Id]->start1[0];
  const MInt* endface = &m_solver->m_windowInfo->plenumInletSurfaceIndices[Id]->end1[0];
  const MInt face = m_physicalBCMap[bcId]->face;

  // Here we find out the normal direction of the boundary and the two tangential
  // directions. This way we can make a general formulation of the boundary condition
  const MInt normalDir = face / 2;
  const MInt directionT1 = (normalDir + 1) % nDim; // T1 stands for first tangential dir.
  const MInt directionT2 = (normalDir + 2) % nDim; // T2 stands for second tangential dir.

  const MInt startN = startface[normalDir];
  const MInt startT1 = startface[directionT1];
  const MInt endT1 = endface[directionT1];
  const MInt startT2 = startface[directionT2];
  const MInt endT2 = endface[directionT2];

  const MInt IJK[nDim] = {1, m_nCells[2], m_nCells[1] * m_nCells[2]};
  const MInt inc[nDim] = {IJK[normalDir], IJK[directionT1], IJK[directionT2]};


  MFloat surface = F0;
  m_plenumSurface = F0;
  const MFloat ii = startN;
  for(MInt t1 = startT1; t1 < endT1; t1++) {
    for(MInt t2 = startT2; t2 < endT2; t2++) {
      const MInt cellId = ii * inc[0] + t1 * inc[1] + t2 * inc[2];
      surface += sqrt(POW2(m_cells->surfaceMetrics[normalDir * nDim + 0][cellId])
                      + POW2(m_cells->surfaceMetrics[normalDir * nDim + 1][cellId])
                      + POW2(m_cells->surfaceMetrics[normalDir * nDim + 2][cellId]));
    }
  }
  // communicate surface and add up!!!!
  cout << "local surface " << surface << endl;
  MPI_Allreduce(&surface, &m_plenumSurface, 1, MPI_DOUBLE, MPI_SUM, m_solver->m_plenumComm, AT_, "surface",
                "m_plenumSurface");

  m_log << "######### Plenum boundary condition was activated ##########" << endl;
  cout << "surface of plenum is: " << m_plenumSurface << endl;


  mTerm(1, AT_, "After init of plenum surface BC ");
}


/**
 * \brief Rescaling inflow
 *
 */
template <MBool isRans>
void StructuredBndryCnd3D<isRans>::initBc2500(MInt bcId) {
  MInt* start = m_physicalBCMap[bcId]->start1;
  MInt* end = m_physicalBCMap[bcId]->end1;

  for(MInt var = 0; var < PV->noVariables; var++) {
    for(MInt i = start[0]; i < end[0]; i++) {
      for(MInt j = start[1]; j < end[1]; j++) {
        for(MInt k = start[2]; k < end[2]; k++) {
          MInt cellId = cellIndex(i, j, k);
          MInt cellIdAdj = cellIndex(m_noGhostLayers, j, k);
          m_cells->pvariables[var][cellId] = m_cells->pvariables[var][cellIdAdj];
        }
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
void StructuredBndryCnd3D<isRans>::initBc2600(MInt bcId) {
  MInt* start = m_physicalBCMap[bcId]->start1;
  MInt* end = m_physicalBCMap[bcId]->end1;

  bc2003(bcId);

  m_solver->m_bc2600 = true;
  const MInt maxNoVariables = m_solver->m_maxNoVariables;

  switch(m_physicalBCMap[bcId]->face) {
    case 0:
    case 1: {
      if(m_solver->m_bc2600InitialStartup) {
        // First copy values from the field into the ghostcells
        for(MInt i = start[0]; i < end[0]; i++) {
          for(MInt j = start[1]; j < end[1]; j++) {
            for(MInt k = start[2]; k < end[2]; k++) {
              const MInt cellId = cellIndex(i, j, k);
              const MInt cellIdAdj = cellIndex(m_noGhostLayers, j, k);

              for(MInt var = 0; var < maxNoVariables; var++) {
                m_cells->pvariables[var][cellId] = m_cells->pvariables[var][cellIdAdj];
              }
            }
          }
        }

        // Fix diagonal cells at end of domain
        if(m_solver->m_bc2600noOffsetCells[2] == 0
           && m_solver->m_bc2600noOffsetCells[1] + m_solver->m_bc2600noActiveCells[1] == m_grid->getMyBlockNoCells(1)) {
          for(MInt i = 0; i < m_solver->m_bc2600noCells[2]; i++) {
            for(MInt k = 0; k < m_solver->m_bc2600noCells[0]; k++) {
              const MInt cellIdA2 = cellIndex(i, m_noGhostLayers + m_solver->m_bc2600noActiveCells[1] - 2, k);
              const MInt cellIdA1 = cellIndex(i, m_noGhostLayers + m_solver->m_bc2600noActiveCells[1] - 1, k);
              const MInt cellIdG1 = cellIndex(i, m_noGhostLayers + m_solver->m_bc2600noActiveCells[1], k);
              for(MInt var = 0; var < maxNoVariables; var++) {
                const MFloat distA1A2 =
                    sqrt(POW2(m_cells->coordinates[0][cellIdA1] - m_cells->coordinates[0][cellIdA2])
                         + POW2(m_cells->coordinates[1][cellIdA1] - m_cells->coordinates[1][cellIdA2])
                         + POW2(m_cells->coordinates[2][cellIdA1] - m_cells->coordinates[2][cellIdA2]));
                const MFloat slope =
                    (m_cells->pvariables[var][cellIdA1] - m_cells->pvariables[var][cellIdA2]) / distA1A2;
                const MFloat distG1A1 =
                    sqrt(POW2(m_cells->coordinates[0][cellIdG1] - m_cells->coordinates[0][cellIdA1])
                         + POW2(m_cells->coordinates[1][cellIdG1] - m_cells->coordinates[1][cellIdA1])
                         + POW2(m_cells->coordinates[2][cellIdG1] - m_cells->coordinates[2][cellIdA1]));
                m_cells->pvariables[var][cellIdG1] = m_cells->pvariables[var][cellIdA1] + distG1A1 * slope;
              }
            }
          }
        }
      }


      // Then copy the values from the ghostcells into restart field
      for(MInt i = start[0]; i < end[0]; i++) {
        for(MInt j = start[1]; j < end[1]; j++) {
          for(MInt k = start[2]; k < end[2]; k++) {
            const MInt cellId = cellIndex(i, j, k);
            const MInt cellIdBc = i + (j + k * m_solver->m_bc2600noCells[1]) * m_solver->m_bc2600noCells[2];
            for(MInt var = 0; var < maxNoVariables; var++) {
              m_solver->m_bc2600Variables[var][cellIdBc] = m_cells->pvariables[var][cellId];
            }
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

/**
 * \brief Prescribe profile BC
 *
 *  In case of the initialStartup2601 flag is set
 *  (no matter if it was a restart or a initial startup)
 *  we load the values from the field to the GC and
 */
template <MBool isRans>
void StructuredBndryCnd3D<isRans>::initBc2601(MInt bcId) {
  MInt* start = m_physicalBCMap[bcId]->start1;
  MInt* end = m_physicalBCMap[bcId]->end1;

  m_solver->m_bc2601 = true;

  m_2601wave = true;

  if(m_2601wave) {
    const MFloat timeConversion = 11.59597694e6;
    const MFloat timeShift = 10.624277309487987 / 0.015795;
    // string fileName = "coefficients_actuated_surface.dat";
    MInt step = 0;
    m_2601noCoeff = 5; // number of effective coefficients
    m_2601noPos = 4;   // number of streamwise positions
    const MString fileNames[4] = {"coefficients_actuated_surface_00138.dat",
                                  "coefficients_actuated_surface_0014452.dat",
                                  "coefficients_actuated_surface_00148.dat", "coefficients_actuated_surface_00153.dat"};
    ifstream infile1(fileNames[0].c_str());
    MInt noLines = count(std::istreambuf_iterator<char>(infile1), std::istreambuf_iterator<char>(), '\n');

    cout << "NoLines: " << noLines << endl;

    mAlloc(m_2601streamwisePos, m_2601noPos, AT_, F0, "m_2601streamwisePos");
    m_2601streamwisePos[0] = 0.0138;
    m_2601streamwisePos[1] = 0.014452;
    m_2601streamwisePos[2] = 0.0148;
    m_2601streamwisePos[3] = 0.0153;

    noLines++;
    mAlloc(m_2601effConst, noLines, m_2601noCoeff * m_2601noPos, AT_, F0, "m_2601effConst");

    cout << "Reading in coefficients from file..." << endl;
    for(MInt pos = 0; pos < m_2601noPos; pos++) {
      step = 0;
      ifstream infile2(fileNames[pos].c_str());
      while(!infile2.eof() && step < noLines - 1) {
        string line[5];
        for(MInt i = 0; i < 5; i++) {
          infile2 >> line[i];
        }

        for(MInt i = 0; i < 5; i++) {
          m_2601effConst[step][pos * m_2601noCoeff + i] = atof(line[i].c_str());

          // if(m_solver->domainId() == 0) {
          //   cout << "step: " << step << " pos: " << pos << " coeff: " << i << endl;
          // }
        }

        step++;
      }
    }
    cout << "Reading in coefficients from file... FINISHED!" << endl;

    const MFloat zeroTime = m_2601effConst[0][0];
    for(MInt pos = 0; pos < m_2601noPos; pos++) {
      for(MInt iter = 0; iter < noLines; iter++) {
        m_2601effConst[iter][pos * m_2601noCoeff + 0] =
            (timeConversion * (m_2601effConst[iter][pos * m_2601noCoeff + 0] - zeroTime))
            - timeShift; //+m_solver->m_time
      }
    }

    m_2601noSteps = step;
  }

  switch(m_physicalBCMap[bcId]->face) {
    case 2: {
      // Copy into zeroth order solution field
      for(MInt i = start[0]; i < end[0]; i++) {
        for(MInt j = start[1]; j < end[1]; j++) {
          for(MInt k = start[2]; k < end[2]; k++) {
            const MInt cellId = cellIndex(i, j, k);
            const MInt cellIdBc = i + (j + k * m_noGhostLayers) * m_solver->m_nCells[2];
            for(MInt var = 0; var < PV->noVariables; var++) {
              m_solver->m_bc2601ZerothOrderSolution[var][cellIdBc] = m_cells->pvariables[var][cellId];
            }
          }
        }
      }


      // Copy the values from the ghostcells into restart field
      for(MInt i = m_noGhostLayers; i < end[0] - m_noGhostLayers; i++) {
        for(MInt j = 0; j < m_noGhostLayers; j++) {
          for(MInt k = m_noGhostLayers; k < end[2] - m_noGhostLayers; k++) {
            const MInt cellId = cellIndex(i, j, k);
            const MInt cellIdBc =
                i - m_noGhostLayers + (j + (k - m_noGhostLayers) * m_noGhostLayers) * m_solver->m_nActiveCells[2];


            for(MInt var = 0; var < PV->noVariables; var++) {
              m_solver->m_bc2601Variables[var][cellIdBc] = m_cells->pvariables[var][cellId];
            }
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


/**
 * \brief Init for the acoustic and entropy waves
 *
 * \author Thomas Schilden
 * \date 12.2.2015
 */
template <MBool isRans>
void StructuredBndryCnd3D<isRans>::initBc2700(MInt bcId) {
  TRACE();
  m_log << endl << "initBc2700 for bcId " << bcId << endl;
  // m_log << " running for " << m_cutOffBndryCndIds[bcId] << endl;

  const MFloat time = m_solver->m_physicalTime;

  // 1. check the modes and allocate memory
  // modeSr
  m_modes = Context::propertyLength("modeSr", m_solverId);
  MFloatScratchSpace modeSr(m_modes, AT_, "modeSr");
  for(MInt i = 0; i < m_modes; i++)
    modeSr(i) = Context::getSolverProperty<MFloat>("modeSr", m_solverId, AT_, i);
  // m_modeAmp
  if(m_modes != Context::propertyLength("modeAmp", m_solverId)) mTerm(1, AT_, "modeAmp does not fit modeSr");
  mAlloc(m_modeAmp, m_modes, "m_modeAmp", F0, AT_);
  for(MInt i = 0; i < m_modes; i++)
    m_modeAmp[i] = Context::getSolverProperty<MFloat>("modeAmp", m_solverId, AT_, i);
  // m_modeType
  if(m_modes != Context::propertyLength("modeType", m_solverId)) mTerm(1, AT_, "modeType does not fit modeSr");
  mAlloc(m_modeType, m_modes, "m_modeType", 0, AT_);
  for(MInt i = 0; i < m_modes; i++)
    m_modeType[i] = Context::getSolverProperty<MInt>("modeType", m_solverId, AT_, i);
  // m_modePhi
  if(m_modes != Context::propertyLength("modePhi", m_solverId)) mTerm(1, AT_, "modePhi does not fit modeSr");
  mAlloc(m_modePhi, m_modes, "m_modePhi", F0, AT_);
  for(MInt i = 0; i < m_modes; i++)
    m_modePhi[i] = Context::getSolverProperty<MFloat>("modePhi", m_solverId, AT_, i) * PI / 180.0;
  // modeAngle
  if(m_modes * 2 != Context::propertyLength("modeAngle", m_solverId)) mTerm(1, AT_, "modeAngle does not fit modeSr");
  MFloatScratchSpace modeAngle(m_modes, 2, AT_, "modeAngle");
  modeAngle.fill(F0);
  for(MInt i = 0; i < m_modes; i++)
    for(MInt j = 0; j < nDim - 1; j++)
      modeAngle(i, j) = Context::getSolverProperty<MFloat>("modeAngle", m_solverId, AT_, i * 2 + j);
  // nmbrOfModes
  if(m_modes != Context::propertyLength("nmbrOfModes", m_solverId)) mTerm(1, AT_, "nmbrOfModes does not fit modeSr");
  mAlloc(m_nmbrOfModes, m_modes, "m_nmbrOfModes", 0, AT_);
  for(MInt i = 0; i < m_modes; i++)
    m_nmbrOfModes[i] = Context::getSolverProperty<MInt>("nmbrOfModes", m_solverId, AT_, i);
  // memberVariables
  mAlloc(m_modeOmega, m_modes, "m_modeOmega", F0, AT_);
  mAlloc(m_modeEtaMin, m_modes, "m_modeEtaMin", F0, AT_);
  mAlloc(m_modeK, m_modes, nDim, "m_modeK", F0, AT_);

  // non-mode specific stuff
  MFloat initTime;

  m_solver->m_restartBc2800 = F0; // has to be changed, if you add BC2800

  if(m_solver->m_restartBc2800)
    initTime = m_solver->m_restartTimeBc2800;
  else {
    initTime = time;
    m_solver->m_restartTimeBc2800 = time;
  }

  m_log << "time        = " << time << endl;
  m_log << "initTime    = " << initTime << endl;

  MFloat UInfinity = F0;
  for(MInt i = 0; i < nDim; i++)
    UInfinity += POW2(PV->VVInfinity[i]);
  UInfinity = sqrt(UInfinity);

  m_log << "solver: UInfinity " << UInfinity << endl
        << "solver: m_referenceLength " << m_solver->m_referenceLength << endl;
  for(MInt mode = 0; mode < m_modes; mode++) {
    m_log << "-- mode " << mode << " --" << endl;
    m_log << "   modeType = " << m_modeType[mode] << endl;
    m_log << "   modeSr = " << modeSr(mode) << endl;
    m_log << "   modeAmp = " << m_modeAmp[mode] << endl;
    m_log << "   modePhi = " << m_modePhi[mode] << endl;
    m_log << "   nmbrOfModes= " << m_nmbrOfModes[mode] << endl;
    for(MInt i = 0; i < (nDim - 1); i++)
      m_log << "   modeAngle(" << i << ") = " << modeAngle(mode, i) << " (rad)" << endl;
    // 3. calculate nondimensional mode properties
    // 3.1 omega
    m_modeOmega[mode] = F2 * PI * modeSr(mode) * UInfinity / m_solver->m_referenceLength;
    // 3.2 wave vector modeK
    m_modeK[mode][0] = cos(modeAngle(mode, 0)) * cos(modeAngle(mode, 1));
    m_modeK[mode][1] = sin(modeAngle(mode, 0)) * cos(modeAngle(mode, 1));
    IF_CONSTEXPR(nDim == 3) m_modeK[mode][2] = sin(modeAngle(mode, 1));

    MFloat Uk = F0;
    for(MInt i = 0; i < nDim; i++)
      Uk += PV->VVInfinity[i] * m_modeK[mode][i];

    const MFloat propVel = (MFloat)(m_modeType[mode]) * sqrt(PV->TInfinity);
    const MFloat K = m_modeOmega[mode] / (Uk + propVel);

    for(MInt i = 0; i < nDim; i++)
      m_modeK[mode][i] *= K;
    // output omega and k
    m_log << "   Uk = " << Uk << endl;
    m_log << "   modeOmega = " << m_modeOmega[mode] << endl;
    m_log << "   K = " << K << endl;
    for(MInt i = 0; i < nDim; i++)
      m_log << "   modeK[" << i << "] = " << m_modeK[mode][i] << endl;
    // 4. get the min wave phase
    MFloat modeEtaMin = F0;

    // here new loop...
    MInt* start = m_physicalBCMap[bcId]->start1;
    MInt* end = m_physicalBCMap[bcId]->end1;
    MInt cellId = -1;

    for(MInt k = start[2]; k < end[2]; k++) {
      for(MInt j = start[1]; j < end[1]; j++) {
        for(MInt i = start[0]; i < end[0]; i++) {
          cellId = cellIndex(i, j, k);

          // for( MInt id = 0; id < m_sortedCutOffCells[bcId]->size(); id++ ){
          // MInt cellId = m_sortedCutOffCells[bcId]->a[ id ];

          MFloat eta = F0;
          for(MInt dim = 0; dim < nDim; dim++) {
            // eta += m_solver->a_coordinate( cellId , i) * m_modeK[mode][i];
            eta += m_cells->coordinates[dim][cellId] * m_modeK[mode][dim];
          }
          eta -= m_modeOmega[mode] * initTime;
          // if(cellId == 0)
          if((k == start[2]) && (j == start[1]) && (i == start[0]))
            modeEtaMin = eta;
          else
            modeEtaMin = mMin(modeEtaMin, eta);
        }
      }
    }

    m_modeEtaMin[mode] = modeEtaMin;
  }
  m_log << "leaving the initialization" << endl;
}


/**
 * \brief Synthetic Turbulence Generation
 *
 */
template <MBool isRans>
void StructuredBndryCnd3D<isRans>::initBc7909(MInt bcId) {
  const MInt NU_T = 5;
  // const MInt SIJSIJ = 6;
  // const MInt LENGTH_SCALE = 7;
  // const MInt FLUC_UU = 8;
  // const MInt FLUC_VV = 9;
  // const MInt FLUC_WW = 10;
  // const MInt FLUC_UV = 11;
  // const MInt FLUC_VW = 12;
  // const MInt FLUC_UW = 13;
  // const MInt FLUC_U = 14;
  // const MInt FLUC_V = 15;
  // const MInt FLUC_W = 16;
  MPI_Comm_rank(m_solver->m_commStg, &m_solver->m_stgMyRank);

  // initialize RK Step to zero
  m_solver->m_RKStep = 0;

  mAlloc(m_stgMaxVel, nDim, "m_stgMaxVel", -99999.9, AT_);

  const MInt noCellsJ = m_grid->getMyBlockNoCells(1);

  if(m_solver->m_stgMyRank == m_solver->m_commStgRoot) {
    cout << "Init stgGlobalLengthscales with " << noCellsJ << " cells" << endl;
  }
  mAlloc(m_stgGlobalLengthScales, noCellsJ, 6, "m_stgGlobalLengthScales", F0, AT_);

  MInt* start = m_physicalBCMap[bcId]->start1;
  MInt* end = m_physicalBCMap[bcId]->end1;

  if(m_solver->m_stgMyRank == m_solver->m_commStgRoot) {
    cout << "Initializing BC 7909..." << endl;
  }

  // Compute size of the Box, take inflow as middle x-coordinate
  MFloatScratchSpace bcast_vb(6, AT_, "bcast_vb");

  MInt minIndex = getPointIdFromCell(m_noGhostLayers, m_noGhostLayers, m_noGhostLayers);
  MInt maxIndex = getPointIdFromCell(m_noGhostLayers, m_noGhostLayers, m_noGhostLayers + m_solver->m_nActiveCells[0]);
  MFloat inflowStartLocal[3] = {m_grid->m_coordinates[0][minIndex], m_grid->m_coordinates[1][minIndex],
                                m_grid->m_coordinates[2][minIndex]};
  MFloat inflowEndLocal[3] = {m_grid->m_coordinates[0][maxIndex], m_grid->m_coordinates[1][maxIndex],
                              m_grid->m_coordinates[2][maxIndex]};
  MFloat inflowStart[3] = {99999.9, 99999.9, 99999.9};
  MFloat inflowEnd[3] = {F0, F0, F0};

  MPI_Allreduce(inflowStartLocal, inflowStart, nDim, MPI_DOUBLE, MPI_MIN, m_solver->m_commStg, AT_, "inflowStartLocal",
                "inflowStart");
  MPI_Allreduce(inflowEndLocal, inflowEnd, nDim, MPI_DOUBLE, MPI_MAX, m_solver->m_commStg, AT_, "inflowEndLocal",
                "inflowEnd");

  const MFloat vbDepth = (inflowEnd[2] - inflowStart[2]) * m_solver->m_stgBLT3;
  const MFloat zOffset = F1B2 * (inflowEnd[2] - inflowStart[2]) * (m_solver->m_stgBLT3 - F1);

  if(m_solver->m_stgMyRank == m_solver->m_commStgRoot) {
    m_solver->m_stgRootRank = true;

    // Get the coordinate of the inflow
    bcast_vb[0] = inflowStart[0] - m_solver->m_stgBLT1 * F1B2;
    bcast_vb[1] = inflowStart[1];
    bcast_vb[2] = inflowStart[2] - zOffset;

    bcast_vb[3] = inflowStart[0] + m_solver->m_stgBLT1 * F1B2;
    bcast_vb[4] = inflowStart[1] + m_solver->m_stgBLT2;
    bcast_vb[5] = inflowStart[2] - zOffset + vbDepth;

    if(m_solver->m_stgBLT3 > 1.5) {
      stringstream errorMsg;
      errorMsg << "The factor for the depth of the virtual eddy box is probably too large, BLT3: "
               << m_solver->m_stgBLT3
               << ". This has been changed in the last commit, instead of setting the absolute depth "
               << "in z-direction of the domain, now only a scaling factor needs to be given. A factor "
               << "of 1.0 will make the virtual box as wide as your domain in z-direction." << endl;
      mTerm(1, AT_, errorMsg.str());
    }
  }

  MPI_Bcast(bcast_vb.begin(), 6, MPI_DOUBLE, m_solver->m_commStgRoot, m_solver->m_commStg, AT_, "bcast_vb.begin()");

  m_stgVbStart[0] = bcast_vb[0];
  m_stgVbStart[1] = bcast_vb[1];
  m_stgVbStart[2] = bcast_vb[2];
  m_stgVbEnd[0] = bcast_vb[3];
  m_stgVbEnd[1] = bcast_vb[4];
  m_stgVbEnd[2] = bcast_vb[5];

  if(m_solver->m_stgMyRank == m_solver->m_commStgRoot) {
    cout << "STG Virtual Box depth factor: " << m_solver->m_stgBLT3 << ", total depth: " << vbDepth << endl
         << "STG Virtual Box Information:" << endl
         << "STG Virtual Box Start X: " << bcast_vb[0] << " Y: " << bcast_vb[1] << " Z: " << bcast_vb[2] << endl
         << "STG Virtual Box End X: " << bcast_vb[3] << " Y: " << bcast_vb[4] << " Z: " << bcast_vb[5] << endl;
  }

  // Done computing virtual box
  if(!m_solver->m_zonal) {
    if(m_solver->m_stgInitialStartup) {
      // Reading in RANS Profile for BC
      // Take 3!!! rows in i-direction because of
      // gradient calculations
      for(MInt k = start[2]; k < end[2]; k++) {
        for(MInt j = start[1]; j < end[1]; j++) {
          for(MInt i = start[0]; i < end[0] + 1; i++) {
            const MInt cellId = cellIndex(i, j, k);
            const MInt cellIdadj = cellIndex(m_noGhostLayers, j, k); // constant in K anyway
            const MInt IBC = cellIndexBC(i, j, k);

            // Set all variables like in the field
            for(MInt var = 0; var < PV->noVariables; var++) {
              m_cells->pvariables[var][cellId] = m_cells->pvariables[var][cellIdadj];
            }

            // fill up the STG FQ field with the RANS values
            for(MInt var = 0; var < PV->noVariables; var++) {
              m_cells->stg_fq[var][IBC] = m_cells->pvariables[var][cellId];
            }

            m_cells->stg_fq[NU_T][IBC] = m_cells->fq[FQ->NU_T][cellIdadj];
          }
        }
      }
    } else {
      for(MInt k = start[2]; k < end[2]; k++) {
        for(MInt j = start[1]; j < end[1]; j++) {
          for(MInt i = start[0]; i < end[0]; i++) {
            const MInt cellId = cellIndex(i, j, k);
            const MInt IBC = cellIndexBC(i, j, k);
            for(MInt var = 0; var < PV->noVariables; var++) {
              m_cells->pvariables[var][cellId] = m_cells->stg_fq[var][IBC];
            }
          }
        }
      }
    }
  }


  // create new eddies if it's an initial start or the createNewEddie flag is set
  if(!m_solver->m_restart || m_solver->m_stgCreateNewEddies) {
    MFloat xk1t, xk2t, xk3t, epsik1, epsik2, epsik3;
    MInt nran = m_solver->m_stgMaxNoEddies;
    MFloatScratchSpace bcast_eddies(nran * 6, AT_, "bcast_eddies");

    MFloat eps = 1e-7;

    if(m_solver->m_stgMyRank == m_solver->m_commStgRoot) {
      // cout << "Creating new Eddies inside Virtual Box"  << endl;
      for(MInt n = 0; n < m_solver->m_stgMaxNoEddies; n++) {
        xk1t = m_stgVbStart[0] + generate_rand() * (m_stgVbEnd[0] - m_stgVbStart[0]);
        xk2t = m_stgVbStart[1] + generate_rand_weighted() * (m_stgVbEnd[1] - m_stgVbStart[1]);
        xk3t = m_stgVbStart[2] + generate_rand() * (m_stgVbEnd[2] - m_stgVbStart[2]);

        // cout << "Creating eddie at x: " << xk1t << " , y: " << xk2t << " , z: " << xk3t << endl;

        epsik1 = 2.0 * generate_rand() - 1.0;
        epsik1 = epsik1 / max(abs(epsik1), eps);
        epsik2 = 2.0 * generate_rand() - 1.0;
        epsik2 = epsik2 / max(abs(epsik2), eps);
        epsik3 = 2.0 * generate_rand() - 1.0;
        epsik3 = epsik3 / max(abs(epsik3), eps);

        bcast_eddies[n + nran * 0] = xk1t;
        bcast_eddies[n + nran * 1] = xk2t;
        bcast_eddies[n + nran * 2] = xk3t;
        bcast_eddies[n + nran * 3] = epsik1;
        bcast_eddies[n + nran * 4] = epsik2;
        bcast_eddies[n + nran * 5] = epsik3;
      }
    }

    // Broadcast the new/updated eddies to all relevant processes
    MPI_Bcast(bcast_eddies.begin(), 6 * nran, MPI_DOUBLE, m_solver->m_commStgRoot, m_solver->m_commStg, AT_,
              "bcast_eddies.begin()");

    // Copy data into m_FQeddie vector
    for(MInt n = 0; n < nran; n++) {
      m_solver->m_stgEddies[n][0] = bcast_eddies[n + nran * 0]; // bcast_buffer(n,var);
      m_solver->m_stgEddies[n][1] = bcast_eddies[n + nran * 1];
      m_solver->m_stgEddies[n][2] = bcast_eddies[n + nran * 2];
      m_solver->m_stgEddies[n][3] = bcast_eddies[n + nran * 3];
      m_solver->m_stgEddies[n][4] = bcast_eddies[n + nran * 4];
      m_solver->m_stgEddies[n][5] = bcast_eddies[n + nran * 5];
    }
  }
}


/**
 * \brief Channel flow / Pipe Flow
 * \authors Pascal Meysonnat, Marian Albers
 */
template <MBool isRans>
void StructuredBndryCnd3D<isRans>::initBc2402(MInt bcId) {
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
  } else if(m_solver->m_initialCondition == 1236) {
    /////////////////////////////////////////////
    ///////// Turbulent pipe flow ///////////////
    /////////////////////////////////////////////
    MFloat uTau = m_solver->m_ReTau * m_solver->m_Ma * sqrt(PV->TInfinity) / m_solver->m_Re;
    m_solver->m_deltaP = -4.0 * CV->rhoInfinity * POW2(uTau) * (m_solver->m_channelLength) / m_solver->m_channelHeight;

    m_solver->m_channelPresInlet = PV->PInfinity;
    m_solver->m_channelPresOutlet = PV->PInfinity + m_solver->m_deltaP;

    m_log << "=========== Turb. Pipe Flow BC Summary =========== " << endl;
    m_log << "-->Turbulent pipe flow deltaP: " << m_solver->m_deltaP << endl;
    m_log << "-->pipe friciton velocity: " << uTau << endl;
    m_log << "-->pipe pressure inflow: " << m_solver->m_channelPresInlet << endl;
    m_log << "-->pipe pressure outflow: " << m_solver->m_channelPresOutlet << endl;
    m_log << "=========== Turb. Pipe Flow BC Summary Finished =========== " << endl;
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
        for(MInt k = startface[2]; k < endface[2]; k++) {
          for(MInt j = startface[1]; j < endface[1]; j++) {
            cellId = cellIndex(ii, j, k);
            surface += sqrt(POW2(m_cells->cellMetrics[0][cellId]) + POW2(m_cells->cellMetrics[1][cellId])
                            + POW2(m_cells->cellMetrics[2][cellId]));
          }
        }
      } else {
        ii = endface[0];
        // activate this or the method below with the 4 points for the surface calculation
        // this is less exact (for straight surf.) but works better for curved surfaces

        MInt cellId = 0;
        for(MInt k = startface[2]; k < endface[2]; k++) {
          for(MInt j = startface[1]; j < endface[1]; j++) {
            cellId = cellIndex(ii - 1, j, k);
            surface += sqrt(POW2(m_cells->cellMetrics[0][cellId]) + POW2(m_cells->cellMetrics[1][cellId])
                            + POW2(m_cells->cellMetrics[2][cellId]));
          }
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
    case 2: {
      mTerm(1, AT_, "surface calculation for k faces(channel not implemented)");
      break;
    }
    default: {
      mTerm(1, AT_, "surface calculation for given faces(channel not implemented)");
    }
  }
}


// solid wall bc
template <MBool isRans>
void StructuredBndryCnd3D<isRans>::bc1000(MInt bcId) {
  const MInt IJK[nDim] = {1, m_nCells[2], m_nCells[1] * m_nCells[2]};
  const MInt* start = m_physicalBCMap[bcId]->start1;
  const MInt* end = m_physicalBCMap[bcId]->end1;
  const MInt face = m_physicalBCMap[bcId]->face;

  // Here we find out the normal direction of the
  // boundary and the two tangential directions.
  // This way we can make a general formulation of
  // the boundary condition
  const MInt normalDir = face / 2;
  const MInt firstTangentialDir = (normalDir + 1) % nDim;
  const MInt secondTangentialDir = (normalDir + 2) % nDim;
  const MInt normalDirStart = start[normalDir];
  const MInt firstTangentialStart = start[firstTangentialDir];
  const MInt firstTangentialEnd = end[firstTangentialDir];
  const MInt secondTangentialStart = start[secondTangentialDir];
  const MInt secondTangentialEnd = end[secondTangentialDir];
  const MInt inc[nDim] = {IJK[normalDir], IJK[firstTangentialDir], IJK[secondTangentialDir]};

  const MInt n = (face % 2) * 2 - 1;                                //-1,+1
  const MInt g1 = normalDirStart + (MInt)(0.5 - (0.5 * (MFloat)n)); //+1,0
  const MInt g2 = normalDirStart + (MInt)(0.5 + (0.5 * (MFloat)n)); // 0,+1
  const MInt a1 = normalDirStart + (MInt)(0.5 - (1.5 * (MFloat)n)); //+2,-1
  const MInt a2 = normalDirStart + (MInt)(0.5 - (2.5 * (MFloat)n)); //+3,-2


  for(MInt t1 = firstTangentialStart; t1 < firstTangentialEnd; t1++) {
    for(MInt t2 = secondTangentialStart; t2 < secondTangentialEnd; t2++) {
      const MInt cellIdG1 = g1 * inc[0] + t1 * inc[1] + t2 * inc[2];
      const MInt cellIdG2 = g2 * inc[0] + t1 * inc[1] + t2 * inc[2];
      const MInt cellIdA1 = a1 * inc[0] + t1 * inc[1] + t2 * inc[2];
      const MInt cellIdA2 = a2 * inc[0] + t1 * inc[1] + t2 * inc[2];

      const MFloat p1 = m_cells->pvariables[PV->P][cellIdA1];
      const MFloat rho1 = m_cells->pvariables[PV->RHO][cellIdA1];
      const MFloat u1 = m_cells->pvariables[PV->U][cellIdA1];
      const MFloat v1 = m_cells->pvariables[PV->V][cellIdA1];
      const MFloat w1 = m_cells->pvariables[PV->W][cellIdA1];

      const MFloat p2 = m_cells->pvariables[PV->P][cellIdA2];
      const MFloat rho2 = m_cells->pvariables[PV->RHO][cellIdA2];
      const MFloat u2 = m_cells->pvariables[PV->U][cellIdA2];
      const MFloat v2 = m_cells->pvariables[PV->V][cellIdA2];
      const MFloat w2 = m_cells->pvariables[PV->W][cellIdA2];

      m_cells->pvariables[PV->RHO][cellIdG1] = rho1;
      m_cells->pvariables[PV->U][cellIdG1] = -u1;
      m_cells->pvariables[PV->V][cellIdG1] = -v1;
      m_cells->pvariables[PV->W][cellIdG1] = -w1;
      m_cells->pvariables[PV->P][cellIdG1] = p1;

      m_cells->pvariables[PV->RHO][cellIdG2] = rho2;
      m_cells->pvariables[PV->U][cellIdG2] = -u2;
      m_cells->pvariables[PV->V][cellIdG2] = -v2;
      m_cells->pvariables[PV->W][cellIdG2] = -w2;
      m_cells->pvariables[PV->P][cellIdG2] = p2;

      if(isRans) {
        const MFloat nuTilde1 = m_cells->pvariables[PV->RANS_VAR[0]][cellIdA1];
        const MFloat nuTilde2 = m_cells->pvariables[PV->RANS_VAR[0]][cellIdA2];

        m_cells->pvariables[PV->RANS_VAR[0]][cellIdG1] = -nuTilde1;
        m_cells->pvariables[PV->RANS_VAR[0]][cellIdG2] = -nuTilde2;
      }
    }
  }
}

// isothermal Wall
template <MBool isRans>
void StructuredBndryCnd3D<isRans>::bc1003(MInt bcId) {
  const MInt IJK[nDim] = {1, m_nCells[2], m_nCells[1] * m_nCells[2]};
  MInt* start = m_physicalBCMap[bcId]->start1;
  MInt* end = m_physicalBCMap[bcId]->end1;
  const MInt face = m_physicalBCMap[bcId]->face;

  MFloat temp = m_isothermalWallTemperature * PV->TInfinity;
  const MFloat gamma = m_solver->m_gamma;

  // Here we find out the normal direction of the
  // boundary and the two tangential directions.
  // This way we can make a general formulation of
  // the boundary condition
  const MInt normalDir = face / 2;
  const MInt firstTangentialDir = (normalDir + 1) % nDim;
  const MInt secondTangentialDir = (normalDir + 2) % nDim;
  const MInt normalDirStart = start[normalDir];
  const MInt firstTangentialStart = start[firstTangentialDir];
  const MInt firstTangentialEnd = end[firstTangentialDir];
  const MInt secondTangentialStart = start[secondTangentialDir];
  const MInt secondTangentialEnd = end[secondTangentialDir];
  const MInt inc[nDim] = {IJK[normalDir], IJK[firstTangentialDir], IJK[secondTangentialDir]};

  const MInt n = (face % 2) * 2 - 1;                                //-1,+1
  const MInt g1 = normalDirStart + (MInt)(0.5 - (0.5 * (MFloat)n)); //+1,0
  const MInt g2 = normalDirStart + (MInt)(0.5 + (0.5 * (MFloat)n)); // 0,+1
  const MInt a1 = normalDirStart + (MInt)(0.5 - (1.5 * (MFloat)n)); //+2,-1
  const MInt a2 = normalDirStart + (MInt)(0.5 - (2.5 * (MFloat)n)); //+3,-2


  for(MInt t1 = firstTangentialStart; t1 < firstTangentialEnd; t1++) {
    for(MInt t2 = secondTangentialStart; t2 < secondTangentialEnd; t2++) {
      const MInt cellIdG1 = g1 * inc[0] + t1 * inc[1] + t2 * inc[2];
      const MInt cellIdG2 = g2 * inc[0] + t1 * inc[1] + t2 * inc[2];
      const MInt cellIdA1 = a1 * inc[0] + t1 * inc[1] + t2 * inc[2];
      const MInt cellIdA2 = a2 * inc[0] + t1 * inc[1] + t2 * inc[2];

      const MFloat p1 = m_cells->pvariables[PV->P][cellIdA1];
      const MFloat rho1 = m_cells->pvariables[PV->RHO][cellIdA1];
      const MFloat u1 = m_cells->pvariables[PV->U][cellIdA1];
      const MFloat v1 = m_cells->pvariables[PV->V][cellIdA1];
      const MFloat w1 = m_cells->pvariables[PV->W][cellIdA1];

      const MFloat p2 = m_cells->pvariables[PV->P][cellIdA2];
      const MFloat rho2 = m_cells->pvariables[PV->RHO][cellIdA2];
      const MFloat u2 = m_cells->pvariables[PV->U][cellIdA2];
      const MFloat v2 = m_cells->pvariables[PV->V][cellIdA2];
      const MFloat w2 = m_cells->pvariables[PV->W][cellIdA2];

      const MFloat rhoWall = p1 * gamma / temp;
      const MFloat rhoG1 = F2 * rhoWall - rho1;
      const MFloat rhoG2 = F2 * rhoWall - rho2;

      m_cells->pvariables[PV->RHO][cellIdG1] = rhoG1;
      m_cells->pvariables[PV->U][cellIdG1] = -u1;
      m_cells->pvariables[PV->V][cellIdG1] = -v1;
      m_cells->pvariables[PV->W][cellIdG1] = -w1;
      m_cells->pvariables[PV->P][cellIdG1] = p1;

      m_cells->pvariables[PV->RHO][cellIdG2] = rhoG2;
      m_cells->pvariables[PV->U][cellIdG2] = -u2;
      m_cells->pvariables[PV->V][cellIdG2] = -v2;
      m_cells->pvariables[PV->W][cellIdG2] = -w2;
      m_cells->pvariables[PV->P][cellIdG2] = p2;
    }
  }
}

/**
 * \brief Moving rigid wall functions
 *
 * \author Marian Albers
 */
template <MBool isRans>
void StructuredBndryCnd3D<isRans>::bc1004(MInt bcId) {
  MInt* start = m_physicalBCMap[bcId]->start1;
  MInt* end = m_physicalBCMap[bcId]->end1;

  const MInt IJK[nDim] = {1, m_nCells[2], m_nCells[1] * m_nCells[2]};
  const MInt IJKP[nDim] = {1, m_nPoints[2], m_nPoints[1] * m_nPoints[2]};

  const MInt pp[3][12] = {
      {0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1}, {0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 1}, {0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0}};

  // MFloat fdt = F0;
  // if (m_solver->m_RKStep!=0) {
  //   fdt = F1/(m_solver->m_timeStep*m_solver->m_RKalpha[m_solver->m_RKStep-1]);
  // } else {
  //   fdt = F1/(m_solver->m_timeStep*m_solver->m_RKalpha[m_solver->m_noRKSteps-1]);
  // }

  const MInt face = m_physicalBCMap[bcId]->face;

  // Here we find out the normal direction of the
  // boundary and the two tangential directions.
  // This way we can make a general formulation of
  // the boundary condition
  const MInt normalDir = face / 2;
  const MInt firstTangentialDir = (normalDir + 1) % nDim;
  const MInt secondTangentialDir = (normalDir + 2) % nDim;
  const MInt normalDirStart = start[normalDir];
  const MInt firstTangentialStart = start[firstTangentialDir];
  const MInt firstTangentialEnd = end[firstTangentialDir];
  const MInt secondTangentialStart = start[secondTangentialDir];
  const MInt secondTangentialEnd = end[secondTangentialDir];
  const MInt inc[nDim] = {IJK[normalDir], IJK[firstTangentialDir], IJK[secondTangentialDir]};
  const MInt incp[nDim] = {IJKP[normalDir], IJKP[firstTangentialDir], IJKP[secondTangentialDir]};

  const MInt n = (face % 2) * 2 - 1;                                //-1,+1
  const MInt g1 = normalDirStart + (MInt)(0.5 - (0.5 * (MFloat)n)); //+1,0
  const MInt g2 = normalDirStart + (MInt)(0.5 + (0.5 * (MFloat)n)); // 0,+1
  const MInt a1 = normalDirStart + (MInt)(0.5 - (1.5 * (MFloat)n)); //+2,-1
  const MInt a2 = normalDirStart + (MInt)(0.5 - (2.5 * (MFloat)n)); //+3,-2

  const MInt g1p = normalDirStart + 2 * ((MInt)(0.5 - (0.5 * (MFloat)n))); //+1,0

  for(MInt t1 = firstTangentialStart; t1 < firstTangentialEnd; t1++) {
    for(MInt t2 = secondTangentialStart; t2 < secondTangentialEnd; t2++) {
      const MInt cellIdG1 = g1 * inc[0] + t1 * inc[1] + t2 * inc[2];
      const MInt cellIdG2 = g2 * inc[0] + t1 * inc[1] + t2 * inc[2];
      const MInt cellIdA1 = a1 * inc[0] + t1 * inc[1] + t2 * inc[2];
      const MInt cellIdA2 = a2 * inc[0] + t1 * inc[1] + t2 * inc[2];

      // compute four surrounding points of surface centroid
      const MInt ijk = g1p * incp[0] + t1 * incp[1] + t2 * incp[2];
      const MInt pp1 = getPointIdfromPoint(ijk, pp[normalDir][0], pp[normalDir][1], pp[normalDir][2]);
      const MInt pp2 = getPointIdfromPoint(ijk, pp[normalDir][3], pp[normalDir][4], pp[normalDir][5]);
      const MInt pp3 = getPointIdfromPoint(ijk, pp[normalDir][6], pp[normalDir][7], pp[normalDir][8]);
      const MInt pp4 = getPointIdfromPoint(ijk, pp[normalDir][9], pp[normalDir][10], pp[normalDir][11]);

      // compute the velocity of the surface centroid
      MFloat gridVel[3] = {F0, F0, F0};
      MFloat gridAcc[3] = {F0, F0, F0};
      MFloat firstVec[3] = {F0, F0, F0};
      MFloat secondVec[3] = {F0, F0, F0};
      MFloat normalVec[3] = {F0, F0, F0};
      for(MInt dim = 0; dim < nDim; dim++) {
        firstVec[dim] = m_grid->m_coordinates[dim][pp2] - m_grid->m_coordinates[dim][pp1];
        secondVec[dim] = m_grid->m_coordinates[dim][pp3] - m_grid->m_coordinates[dim][pp1];

        gridVel[dim] = F1B4
                       * (m_grid->m_velocity[dim][pp1] + m_grid->m_velocity[dim][pp2] + m_grid->m_velocity[dim][pp3]
                          + m_grid->m_velocity[dim][pp4]);

        gridAcc[dim] = F1B4
                       * (m_grid->m_acceleration[dim][pp1] + m_grid->m_acceleration[dim][pp2]
                          + m_grid->m_acceleration[dim][pp3] + m_grid->m_acceleration[dim][pp4]);
      }

      const MFloat p1 = m_cells->pvariables[PV->P][cellIdA1];
      const MFloat rho1 = m_cells->pvariables[PV->RHO][cellIdA1];
      const MFloat u1 = m_cells->pvariables[PV->U][cellIdA1];
      const MFloat v1 = m_cells->pvariables[PV->V][cellIdA1];
      const MFloat w1 = m_cells->pvariables[PV->W][cellIdA1];

      const MFloat p2 = m_cells->pvariables[PV->P][cellIdA2];
      const MFloat rho2 = m_cells->pvariables[PV->RHO][cellIdA2];
      const MFloat u2 = m_cells->pvariables[PV->U][cellIdA2];
      const MFloat v2 = m_cells->pvariables[PV->V][cellIdA2];
      const MFloat w2 = m_cells->pvariables[PV->W][cellIdA2];


      // compute normal acceleration of surface to derive von Neumann BC for pressure and density
      // compute normal vector of surface
      crossProduct(normalVec, firstVec, secondVec);
      const MFloat normalLength = sqrt(POW2(normalVec[0]) + POW2(normalVec[1]) + POW2(normalVec[2]));

      for(MInt dim = 0; dim < nDim; dim++) {
        normalVec[dim] /= normalLength;
      }

      // compute distance between ghost point and active point
      const MFloat dx = m_cells->coordinates[0][cellIdA1] - m_cells->coordinates[0][cellIdG1];
      const MFloat dy = m_cells->coordinates[1][cellIdA1] - m_cells->coordinates[1][cellIdG1];
      const MFloat dz = m_cells->coordinates[2][cellIdA1] - m_cells->coordinates[2][cellIdG1];
      const MFloat dn = sqrt(POW2(dx) + POW2(dy) + POW2(dz));

      const MFloat surfTemp = m_solver->m_gamma * p1 / rho1;
      // compute normal acceleration of surface
      MFloat an = F0;
      for(MInt dim = 0; dim < nDim; dim++) {
        an += gridAcc[dim] * normalVec[dim];
      }

      const MFloat beta = m_solver->m_gamma * an / surfTemp;
      const MFloat fac = (F2 + dn * beta) / (F2 - dn * beta);

      const MFloat pG1 = p1 * fac;
      const MFloat pG2 = p2 * fac;
      const MFloat rhoG1 = rho1 * fac;
      const MFloat rhoG2 = rho2 * fac;

      const MFloat uG1 = F2 * gridVel[0] - u1;
      const MFloat vG1 = F2 * gridVel[1] - v1;
      const MFloat wG1 = F2 * gridVel[2] - w1;

      const MFloat uG2 = F2 * gridVel[0] - u2;
      const MFloat vG2 = F2 * gridVel[1] - v2;
      const MFloat wG2 = F2 * gridVel[2] - w2;

      m_cells->pvariables[PV->RHO][cellIdG1] = rhoG1;
      m_cells->pvariables[PV->U][cellIdG1] = uG1;
      m_cells->pvariables[PV->V][cellIdG1] = vG1;
      m_cells->pvariables[PV->W][cellIdG1] = wG1;
      m_cells->pvariables[PV->P][cellIdG1] = pG1;

      m_cells->pvariables[PV->RHO][cellIdG2] = rhoG2;
      m_cells->pvariables[PV->U][cellIdG2] = uG2;
      m_cells->pvariables[PV->V][cellIdG2] = vG2;
      m_cells->pvariables[PV->W][cellIdG2] = wG2;
      m_cells->pvariables[PV->P][cellIdG2] = pG2;
    }
  }
}

/**
 * \brief Moving rigid isothermal wall functions
 *
 * \author Marian Albers
 */
template <MBool isRans>
void StructuredBndryCnd3D<isRans>::bc1006(MInt bcId) {
  MInt* start = m_physicalBCMap[bcId]->start1;
  MInt* end = m_physicalBCMap[bcId]->end1;

  const MInt IJK[nDim] = {1, m_nCells[2], m_nCells[1] * m_nCells[2]};
  const MInt IJKP[nDim] = {1, m_nPoints[2], m_nPoints[1] * m_nPoints[2]};

  const MInt pp[3][12] = {
      {0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1}, {0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 1}, {0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0}};

  const MInt face = m_physicalBCMap[bcId]->face;

  MFloat temp = m_isothermalWallTemperature * PV->TInfinity;
  const MFloat gamma = m_solver->m_gamma;

  // Here we find out the normal direction of the
  // boundary and the two tangential directions.
  // This way we can make a general formulation of
  // the boundary condition
  const MInt normalDir = face / 2;
  const MInt firstTangentialDir = (normalDir + 1) % nDim;
  const MInt secondTangentialDir = (normalDir + 2) % nDim;
  const MInt normalDirStart = start[normalDir];
  const MInt firstTangentialStart = start[firstTangentialDir];
  const MInt firstTangentialEnd = end[firstTangentialDir];
  const MInt secondTangentialStart = start[secondTangentialDir];
  const MInt secondTangentialEnd = end[secondTangentialDir];
  const MInt inc[nDim] = {IJK[normalDir], IJK[firstTangentialDir], IJK[secondTangentialDir]};
  const MInt incp[nDim] = {IJKP[normalDir], IJKP[firstTangentialDir], IJKP[secondTangentialDir]};

  const MInt n = (face % 2) * 2 - 1;                                //-1,+1
  const MInt g1 = normalDirStart + (MInt)(0.5 - (0.5 * (MFloat)n)); //+1,0
  const MInt g2 = normalDirStart + (MInt)(0.5 + (0.5 * (MFloat)n)); // 0,+1
  const MInt a1 = normalDirStart + (MInt)(0.5 - (1.5 * (MFloat)n)); //+2,-1
  const MInt a2 = normalDirStart + (MInt)(0.5 - (2.5 * (MFloat)n)); //+3,-2

  const MInt g1p = normalDirStart + 2 * ((MInt)(0.5 - (0.5 * (MFloat)n))); //+1,0

  for(MInt t1 = firstTangentialStart; t1 < firstTangentialEnd; t1++) {
    for(MInt t2 = secondTangentialStart; t2 < secondTangentialEnd; t2++) {
      const MInt cellIdG1 = g1 * inc[0] + t1 * inc[1] + t2 * inc[2];
      const MInt cellIdG2 = g2 * inc[0] + t1 * inc[1] + t2 * inc[2];
      const MInt cellIdA1 = a1 * inc[0] + t1 * inc[1] + t2 * inc[2];
      const MInt cellIdA2 = a2 * inc[0] + t1 * inc[1] + t2 * inc[2];

      // compute four surrounding points of surface centroid
      const MInt ijk = g1p * incp[0] + t1 * incp[1] + t2 * incp[2];
      const MInt pp1 = getPointIdfromPoint(ijk, pp[normalDir][0], pp[normalDir][1], pp[normalDir][2]);
      const MInt pp2 = getPointIdfromPoint(ijk, pp[normalDir][3], pp[normalDir][4], pp[normalDir][5]);
      const MInt pp3 = getPointIdfromPoint(ijk, pp[normalDir][6], pp[normalDir][7], pp[normalDir][8]);
      const MInt pp4 = getPointIdfromPoint(ijk, pp[normalDir][9], pp[normalDir][10], pp[normalDir][11]);

      // compute the velocity of the surface centroid
      MFloat gridVel[3] = {F0, F0, F0};
      MFloat gridAcc[3] = {F0, F0, F0};
      MFloat firstVec[3] = {F0, F0, F0};
      MFloat secondVec[3] = {F0, F0, F0};
      MFloat normalVec[3] = {F0, F0, F0};
      for(MInt dim = 0; dim < nDim; dim++) {
        firstVec[dim] = m_grid->m_coordinates[dim][pp2] - m_grid->m_coordinates[dim][pp1];
        secondVec[dim] = m_grid->m_coordinates[dim][pp3] - m_grid->m_coordinates[dim][pp1];

        gridVel[dim] = F1B4
                       * (m_grid->m_velocity[dim][pp1] + m_grid->m_velocity[dim][pp2] + m_grid->m_velocity[dim][pp3]
                          + m_grid->m_velocity[dim][pp4]);

        gridAcc[dim] = F1B4
                       * (m_grid->m_acceleration[dim][pp1] + m_grid->m_acceleration[dim][pp2]
                          + m_grid->m_acceleration[dim][pp3] + m_grid->m_acceleration[dim][pp4]);
      }

      const MFloat p1 = m_cells->pvariables[PV->P][cellIdA1];
      const MFloat rho1 = m_cells->pvariables[PV->RHO][cellIdA1];
      const MFloat u1 = m_cells->pvariables[PV->U][cellIdA1];
      const MFloat v1 = m_cells->pvariables[PV->V][cellIdA1];
      const MFloat w1 = m_cells->pvariables[PV->W][cellIdA1];

      const MFloat p2 = m_cells->pvariables[PV->P][cellIdA2];
      const MFloat rho2 = m_cells->pvariables[PV->RHO][cellIdA2];
      const MFloat u2 = m_cells->pvariables[PV->U][cellIdA2];
      const MFloat v2 = m_cells->pvariables[PV->V][cellIdA2];
      const MFloat w2 = m_cells->pvariables[PV->W][cellIdA2];


      // compute normal acceleration of surface to derive von Neumann BC for pressure and density
      // compute normal vector of surface
      crossProduct(normalVec, firstVec, secondVec);
      const MFloat normalLength = sqrt(POW2(normalVec[0]) + POW2(normalVec[1]) + POW2(normalVec[2]));

      for(MInt dim = 0; dim < nDim; dim++) {
        normalVec[dim] /= normalLength;
      }

      // compute distance between ghost point and active point
      const MFloat dx = m_cells->coordinates[0][cellIdA1] - m_cells->coordinates[0][cellIdG1];
      const MFloat dy = m_cells->coordinates[1][cellIdA1] - m_cells->coordinates[1][cellIdG1];
      const MFloat dz = m_cells->coordinates[2][cellIdA1] - m_cells->coordinates[2][cellIdG1];
      const MFloat dn = sqrt(POW2(dx) + POW2(dy) + POW2(dz));

      // compute normal acceleration of surface
      MFloat an = F0;
      for(MInt dim = 0; dim < nDim; dim++) {
        an += gridAcc[dim] * normalVec[dim];
      }

      const MFloat beta = m_solver->m_gamma * an / temp;
      const MFloat fac = (F2 + dn * beta) / (F2 - dn * beta);

      const MFloat pG1 = p1 * fac;
      const MFloat pG2 = p2 * fac;
      const MFloat rhoWall = p1 * gamma / temp;
      const MFloat rhoG1 = (F2 * rhoWall - rho1) * fac;
      const MFloat rhoG2 = (F2 * rhoWall - rho2) * fac;

      const MFloat uG1 = F2 * gridVel[0] - u1;
      const MFloat vG1 = F2 * gridVel[1] - v1;
      const MFloat wG1 = F2 * gridVel[2] - w1;

      const MFloat uG2 = F2 * gridVel[0] - u2;
      const MFloat vG2 = F2 * gridVel[1] - v2;
      const MFloat wG2 = F2 * gridVel[2] - w2;

      m_cells->pvariables[PV->RHO][cellIdG1] = rhoG1;
      m_cells->pvariables[PV->U][cellIdG1] = uG1;
      m_cells->pvariables[PV->V][cellIdG1] = vG1;
      m_cells->pvariables[PV->W][cellIdG1] = wG1;
      m_cells->pvariables[PV->P][cellIdG1] = pG1;

      m_cells->pvariables[PV->RHO][cellIdG2] = rhoG2;
      m_cells->pvariables[PV->U][cellIdG2] = uG2;
      m_cells->pvariables[PV->V][cellIdG2] = vG2;
      m_cells->pvariables[PV->W][cellIdG2] = wG2;
      m_cells->pvariables[PV->P][cellIdG2] = pG2;
    }
  }
}

/**
 * \brief Oscillating wall
 *
 * \author Marian Albers
 */
template <MBool isRans>
void StructuredBndryCnd3D<isRans>::bc1007(MInt bcId) {
  MInt* start = m_physicalBCMap[bcId]->start1;
  MInt* end = m_physicalBCMap[bcId]->end1;
  const MFloat t = m_solver->m_time + m_solver->m_timeStep * m_solver->m_RKalpha[m_solver->m_RKStep];

  switch(m_physicalBCMap[bcId]->face) {
    case 2: {
      MFloat gridVel[3] = {F0, F0, F0};

      for(MInt k = start[2]; k < end[2]; k++) {
        for(MInt i = start[0]; i < end[0]; i++) {
          MInt cellIdG2 = cellIndex(i, 0, k); // ghost
          MInt cellIdG1 = cellIndex(i, 1, k); // ghost
          MInt cellIdA1 = cellIndex(i, 2, k); // field
          MInt cellIdA2 = cellIndex(i, 3, k); // field

          const MFloat transitionLength = m_solver->m_waveEndTransition - m_solver->m_waveBeginTransition;
          const MFloat xInit = m_cells->coordinates[0][cellIdG1];
          MFloat transitionFactor = F0;
          if(xInit <= m_solver->m_waveBeginTransition) {
            transitionFactor = F0;
          } else if(xInit > m_solver->m_waveBeginTransition && xInit < m_solver->m_waveEndTransition) {
            transitionFactor = (1 - cos((xInit - m_solver->m_waveBeginTransition) / transitionLength * PI)) * F1B2;
          } else {
            transitionFactor = F1;
          }

          gridVel[2] = m_solver->m_waveAmplitude * sin(F2 * PI / m_solver->m_waveTime * t) * transitionFactor;

          MFloat p1 = m_cells->pvariables[PV->P][cellIdA1];
          MFloat rho1 = m_cells->pvariables[PV->RHO][cellIdA1];
          MFloat u1 = m_cells->pvariables[PV->U][cellIdA1];
          MFloat v1 = m_cells->pvariables[PV->V][cellIdA1];
          MFloat w1 = m_cells->pvariables[PV->W][cellIdA1];

          MFloat u2 = m_cells->pvariables[PV->U][cellIdA2];
          MFloat v2 = m_cells->pvariables[PV->V][cellIdA2];
          MFloat w2 = m_cells->pvariables[PV->W][cellIdA2];

          MFloat pG1 = p1;
          MFloat pG2 = p1;
          MFloat rhoG1 = rho1;
          MFloat rhoG2 = rho1;

          MFloat uG1 = F2 * gridVel[0] - u1;
          MFloat vG1 = F2 * gridVel[1] - v1;
          MFloat wG1 = F2 * gridVel[2] - w1;

          const MFloat uG2 = F2 * gridVel[0] - u2;
          const MFloat vG2 = F2 * gridVel[1] - v2;
          const MFloat wG2 = F2 * gridVel[2] - w2;

          m_cells->pvariables[PV->RHO][cellIdG1] = rhoG1;
          m_cells->pvariables[PV->U][cellIdG1] = uG1;
          m_cells->pvariables[PV->V][cellIdG1] = vG1;
          m_cells->pvariables[PV->W][cellIdG1] = wG1;
          m_cells->pvariables[PV->P][cellIdG1] = pG1;

          m_cells->pvariables[PV->RHO][cellIdG2] = rhoG2;
          m_cells->pvariables[PV->U][cellIdG2] = uG2;
          m_cells->pvariables[PV->V][cellIdG2] = vG2;
          m_cells->pvariables[PV->W][cellIdG2] = wG2;
          m_cells->pvariables[PV->P][cellIdG2] = pG2;
        }
      }
      break;
    }
    default: {
      cout << "bc1007: face not implemented" << endl;
    }
  }
}


/**
 * \brief Subsonic rotational inflow
 *
 *  rho=rho_inf, u=u_inf, v=v_inf, w=w_inf, dp/dn=0
 */
template <MBool isRans>
void StructuredBndryCnd3D<isRans>::bc2014(MInt bcId) {
  // implemented for i-direction only for the moment
  MInt* start = m_physicalBCMap[bcId]->start1;
  MInt* end = m_physicalBCMap[bcId]->end1;
  switch(m_physicalBCMap[bcId]->face) {
    case 0: {
      MInt cellId = -1;
      MInt cellIdadj = -1;
      for(MInt k = start[2]; k < end[2]; k++) {
        for(MInt j = start[1]; j < end[1]; j++) {
          for(MInt i = start[0]; i < end[0]; i++) {
            cellId = cellIndex(m_noGhostLayers - 1 - i, j, k);
            cellIdadj = cellIndex(m_noGhostLayers - i, j, k);

            MFloat y = m_cells->coordinates[1][cellId];
            MFloat z = m_cells->coordinates[2][cellId];
            MFloat phi = atan2(y, z);
            MFloat r = sqrt(POW2(y) + POW2(z));
            MFloat rmax = 10.0;

            m_cells->pvariables[PV->RHO][cellId] = CV->rhoInfinity;
            m_cells->pvariables[PV->U][cellId] = PV->UInfinity;
            m_cells->pvariables[PV->V][cellId] = -(r / rmax) * cos(phi) * PV->UInfinity;
            m_cells->pvariables[PV->W][cellId] = (r / rmax) * sin(phi) * PV->UInfinity;
            m_cells->pvariables[PV->P][cellId] = m_cells->pvariables[PV->P][cellIdadj];
          }
        }
      }
      break;
    }
    case 2: {
      for(MInt k = start[2]; k < end[2]; k++) {
        for(MInt j = start[1]; j < end[1]; j++) {
          for(MInt i = start[0]; i < end[0]; i++) {
            const MInt cellId = cellIndex(i, m_noGhostLayers - j - 1, k); // ghost
            const MInt cellIdadj = cellIndex(i, m_noGhostLayers - j, k);  // field
            m_cells->pvariables[PV->RHO][cellId] = CV->rhoInfinity;
            m_cells->pvariables[PV->U][cellId] = PV->UInfinity;
            m_cells->pvariables[PV->V][cellId] = PV->VInfinity;
            m_cells->pvariables[PV->W][cellId] = PV->WInfinity;

            m_cells->pvariables[PV->P][cellId] = m_cells->pvariables[PV->P][cellIdadj];
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
 * \brief Subsonic Inflow
 *
 *  rho=rho_inf, u=u_inf, v=v_inf, w=w_inf, dp/dn=0
 */
template <MBool isRans>
void StructuredBndryCnd3D<isRans>::bc2001(MInt bcId) {
  MInt* start = m_physicalBCMap[bcId]->start1;
  MInt* end = m_physicalBCMap[bcId]->end1;
  MInt face = m_physicalBCMap[bcId]->face;
  const MInt n = face % 2;
  const MInt dim = face / 2;
  MInt offset[3] = {0, 0, 0};
  MInt inc[3] = {1, 1, 1};
  MInt adjInc[3] = {0, 0, 0};

  // on face 0,2,4 set offset and
  // reverse the incrementor
  if(n == 0) {
    offset[dim] = m_noGhostLayers - 1;
    inc[dim] = -1;
    adjInc[dim] = 1;
  } else {
    inc[dim] = 1;
    adjInc[dim] = -1;
  }

  for(MInt k = start[2]; k < end[2]; k++) {
    for(MInt j = start[1]; j < end[1]; j++) {
      for(MInt i = start[0]; i < end[0]; i++) {
        const MInt ii = inc[0] * (i - start[0]) + start[0] + offset[0];
        const MInt jj = inc[1] * (j - start[1]) + start[1] + offset[1];
        const MInt kk = inc[2] * (k - start[2]) + start[2] + offset[2];
        const MInt cellId = cellIndex(ii, jj, kk);
        const MInt cellIdadj = cellIndex(ii + adjInc[0], jj + adjInc[1], kk + adjInc[2]);

        m_cells->pvariables[PV->RHO][cellId] = CV->rhoInfinity;
        m_cells->pvariables[PV->U][cellId] = PV->UInfinity;
        m_cells->pvariables[PV->V][cellId] = PV->VInfinity;
        m_cells->pvariables[PV->W][cellId] = PV->WInfinity;
        if(isRans) {
          m_cells->pvariables[PV->RANS_VAR[0]][cellId] = CV->ransInfinity[0] / CV->rhoInfinity;
        }

        // extrapolate pressure from the field
        m_cells->pvariables[PV->P][cellId] = m_cells->pvariables[PV->P][cellIdadj];
      }
    }
  }
}

/**
 * \brief Subsonic Inflow for a plenum
 *
 * \authors Pascal S. Meysonnat, Volkan
 * \date August 2018
 */
template <MBool isRans>
void StructuredBndryCnd3D<isRans>::bc2097(MInt bcId) {
  (void)bcId;
}


/**
 * \brief Subsonic Inflow with u=(y/delta)^(1/7)
 *
 * rho=rho_inf, u=u_inf, v=v_inf, w=w_inf, dp/dn=0
 */
template <MBool isRans>
void StructuredBndryCnd3D<isRans>::bc2099(MInt bcId) {
  // implemented for i-direction only for the moment
  MInt* start = m_physicalBCMap[bcId]->start1;
  MInt* end = m_physicalBCMap[bcId]->end1;
  MFloat y = F0;
  switch(m_physicalBCMap[bcId]->face) {
    case 0: {
      for(MInt k = start[2]; k < end[2]; k++) {
        for(MInt j = start[1]; j < end[1]; j++) {
          for(MInt i = start[0]; i < end[0]; i++) {
            const MInt cellId = cellIndex(m_noGhostLayers - 1 - i, j, k);
            const MInt cellIdadj = cellIndex(m_noGhostLayers - i, j, k);
            m_cells->pvariables[PV->RHO][cellId] = CV->rhoInfinity;
            m_cells->pvariables[PV->U][cellId] = PV->UInfinity;
            y = m_cells->coordinates[1][cellId];
            if(y < 1.0) {
              m_cells->pvariables[PV->U][cellId] = PV->UInfinity * pow(y, (1.0 / 7.0));
            }
            m_cells->pvariables[PV->V][cellId] = PV->VInfinity;
            m_cells->pvariables[PV->W][cellId] = PV->WInfinity;
            m_cells->pvariables[PV->P][cellId] = m_cells->pvariables[PV->P][cellIdadj];
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
 * \brief Supersonic Inflow
 *
 *  rho=rho_inf, u=u_inf, v=v_inf, w=w_inf, p=p_inf
 */
template <MBool isRans>
void StructuredBndryCnd3D<isRans>::bc2002(MInt bcId) {
  TRACE();
  // implemented for i-direction only for the moment
  MInt* start = m_physicalBCMap[bcId]->start1;
  MInt* end = m_physicalBCMap[bcId]->end1;

  MInt cellId = -1;
  for(MInt k = start[2]; k < end[2]; k++) {
    for(MInt j = start[1]; j < end[1]; j++) {
      for(MInt i = start[0]; i < end[0]; i++) {
        cellId = cellIndex(i, j, k);
        m_cells->pvariables[PV->RHO][cellId] = CV->rhoInfinity;
        m_cells->pvariables[PV->U][cellId] = PV->UInfinity;
        m_cells->pvariables[PV->V][cellId] = PV->VInfinity;
        m_cells->pvariables[PV->W][cellId] = PV->WInfinity;
        m_cells->pvariables[PV->P][cellId] = PV->PInfinity;
        if(isRans) {
          m_cells->pvariables[PV->RANS_VAR[0]][cellId] = CV->ransInfinity[0] / CV->rhoInfinity;
        }
      }
    }
  }
}

/**
 * \brief Supersonic outflow
 */
template <MBool isRans>
void StructuredBndryCnd3D<isRans>::bc2005(MInt bcId) {
  MInt* start = m_physicalBCMap[bcId]->start1;
  MInt* end = m_physicalBCMap[bcId]->end1;
  MInt face = m_physicalBCMap[bcId]->face;
  const MInt n = face % 2;
  const MInt dim = face / 2;
  MInt offset[3] = {0, 0, 0};
  MInt inc[3] = {1, 1, 1};
  MInt adjInc[3] = {0, 0, 0};

  // on face 0,2,4 set offset and
  // reverse the incrementor
  if(n == 0) {
    offset[dim] = m_noGhostLayers - 1;
    inc[dim] = -1;
    adjInc[dim] = 1;
  } else {
    inc[dim] = 1;
    adjInc[dim] = -1;
  }

  for(MInt k = start[2]; k < end[2]; k++) {
    for(MInt j = start[1]; j < end[1]; j++) {
      for(MInt i = start[0]; i < end[0]; i++) {
        const MInt ii = inc[0] * (i - start[0]) + start[0] + offset[0];
        const MInt jj = inc[1] * (j - start[1]) + start[1] + offset[1];
        const MInt kk = inc[2] * (k - start[2]) + start[2] + offset[2];
        const MInt cellId = cellIndex(ii, jj, kk);
        const MInt cellIdadj = cellIndex(ii + adjInc[0], jj + adjInc[1], kk + adjInc[2]);

        for(MInt var = 0; var < PV->noVariables; var++) {
          m_cells->pvariables[var][cellId] = m_cells->pvariables[var][cellIdadj];
        }
      }
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
void StructuredBndryCnd3D<isRans>::bc2003(MInt bcId) {
  const MInt IJK[nDim] = {1, m_nCells[2], m_nCells[1] * m_nCells[2]};
  MInt* start = m_physicalBCMap[bcId]->start1;
  MInt* end = m_physicalBCMap[bcId]->end1;

  // Here we find out the normal direction of the
  // boundary and the two tangential directions.
  // This way we can make a general formulation of
  // the boundary condition
  const MInt face = m_physicalBCMap[bcId]->face;
  const MInt normalDir = face / 2;
  const MInt firstTangentialDir = (normalDir + 1) % nDim;
  const MInt secondTangentialDir = (normalDir + 2) % nDim;
  const MInt normalDirStart = start[normalDir];
  const MInt firstTangentialStart = start[firstTangentialDir];
  const MInt firstTangentialEnd = end[firstTangentialDir];
  const MInt secondTangentialStart = start[secondTangentialDir];
  const MInt secondTangentialEnd = end[secondTangentialDir];
  const MInt inc[nDim] = {IJK[normalDir], IJK[firstTangentialDir], IJK[secondTangentialDir]};

  const MInt n = (face % 2) * 2 - 1;                                //-1,+1
  const MInt g1 = normalDirStart + (MInt)(0.5 - (0.5 * (MFloat)n)); //+1,0
  const MInt g2 = normalDirStart + (MInt)(0.5 + (0.5 * (MFloat)n)); // 0,+1
  const MInt a1 = normalDirStart + (MInt)(0.5 - (1.5 * (MFloat)n)); //+2,-1
  const MInt a2 = normalDirStart + (MInt)(0.5 - (2.5 * (MFloat)n)); //+3,-2

  for(MInt t1 = firstTangentialStart; t1 < firstTangentialEnd; t1++) {
    for(MInt t2 = secondTangentialStart; t2 < secondTangentialEnd; t2++) {
      const MInt cellIdG1 = g1 * inc[0] + t1 * inc[1] + t2 * inc[2];
      const MInt cellIdG2 = g2 * inc[0] + t1 * inc[1] + t2 * inc[2];
      const MInt cellIdA1 = a1 * inc[0] + t1 * inc[1] + t2 * inc[2];
      const MInt cellIdA2 = a2 * inc[0] + t1 * inc[1] + t2 * inc[2];

      const MFloat dxidx = m_cells->surfaceMetrics[normalDir * nDim + 0][cellIdA1];
      const MFloat dxidy = m_cells->surfaceMetrics[normalDir * nDim + 1][cellIdA1];
      const MFloat dxidz = m_cells->surfaceMetrics[normalDir * nDim + 2][cellIdA1];
      // leaving domain of integration in positive coordinate direction,
      // therefore multiply with positive F1
      const MFloat gradxi = n * F1 / sqrt(dxidx * dxidx + dxidy * dxidy + dxidz * dxidz);

      const MFloat rhoInner = m_cells->pvariables[PV->RHO][cellIdA1];
      const MFloat uInner = m_cells->pvariables[PV->U][cellIdA1];
      const MFloat vInner = m_cells->pvariables[PV->V][cellIdA1];
      const MFloat wInner = m_cells->pvariables[PV->W][cellIdA1];
      const MFloat pInner = m_cells->pvariables[PV->P][cellIdA1];

      const MFloat maContravariant =
          (dxidx * uInner + dxidy * vInner + dxidz * wInner - m_cells->dxt[normalDir][cellIdA1]) * gradxi;


      if(maContravariant < F0) {
        // inflow
        const MFloat p = pInner;
        const MFloat rho = CV->rhoInfinity;

        m_cells->pvariables[PV->RHO][cellIdG1] = rho;
        m_cells->pvariables[PV->U][cellIdG1] = PV->UInfinity;
        m_cells->pvariables[PV->V][cellIdG1] = PV->VInfinity;
        m_cells->pvariables[PV->W][cellIdG1] = PV->WInfinity;
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
        m_cells->pvariables[PV->W][cellIdG1] = wInner;
        m_cells->pvariables[PV->P][cellIdG1] = p;
        if(isRans) {
          m_cells->pvariables[PV->RANS_VAR[0]][cellIdG1] =
              (F2 * m_cells->pvariables[PV->RANS_VAR[0]][cellIdA1] - m_cells->pvariables[PV->RANS_VAR[0]][cellIdA2]);
        }
      }

      // extrapolate into second ghost cell
      for(MInt var = 0; var < PV->noVariables; var++) {
        m_cells->pvariables[var][cellIdG2] =
            F2 * m_cells->pvariables[var][cellIdG1] - m_cells->pvariables[var][cellIdA1];
      }
    }
  }
}


/**
 *  \briefCharacteristic boundary condition
 *
 *  Detects whether it is in- or outflow and prescribes
 *  appropriately
 *  \author Marian Albers
 */
template <MBool isRans>
void StructuredBndryCnd3D<isRans>::bc2004(MInt bcId) {
  const MInt IJK[nDim] = {1, m_nCells[2], m_nCells[1] * m_nCells[2]};
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
  const MInt secondTangentialDir = (normalDir + 2) % nDim;
  const MInt normalDirStart = start[normalDir];
  const MInt firstTangentialStart = start[firstTangentialDir];
  const MInt firstTangentialEnd = end[firstTangentialDir];
  const MInt secondTangentialStart = start[secondTangentialDir];
  const MInt secondTangentialEnd = end[secondTangentialDir];
  const MInt inc[nDim] = {IJK[normalDir], IJK[firstTangentialDir], IJK[secondTangentialDir]};

  const MInt n = (face % 2) * 2 - 1;                                //-1,+1
  const MInt g1 = normalDirStart + (MInt)(0.5 - (0.5 * (MFloat)n)); //+1,0
  const MInt g2 = normalDirStart + (MInt)(0.5 + (0.5 * (MFloat)n)); // 0,+1
  const MInt a1 = normalDirStart + (MInt)(0.5 - (1.5 * (MFloat)n)); //+2,-1
  const MInt a2 = normalDirStart + (MInt)(0.5 - (2.5 * (MFloat)n)); //+3,-2

  for(MInt t1 = firstTangentialStart; t1 < firstTangentialEnd; t1++) {
    for(MInt t2 = secondTangentialStart; t2 < secondTangentialEnd; t2++) {
      const MInt cellIdG1 = g1 * inc[0] + t1 * inc[1] + t2 * inc[2];
      const MInt cellIdG2 = g2 * inc[0] + t1 * inc[1] + t2 * inc[2];
      const MInt cellIdA1 = a1 * inc[0] + t1 * inc[1] + t2 * inc[2];
      const MInt cellIdA2 = a2 * inc[0] + t1 * inc[1] + t2 * inc[2];

      const MFloat dxidx = m_cells->surfaceMetrics[normalDir * nDim + 0][cellIdA1];
      const MFloat dxidy = m_cells->surfaceMetrics[normalDir * nDim + 1][cellIdA1];
      const MFloat dxidz = m_cells->surfaceMetrics[normalDir * nDim + 2][cellIdA1];
      // multiply with n, so it will be -1 or +1 depending if we enter
      // or leave the domain of integration in positive direction
      const MFloat gradxi = n * F1 / sqrt(dxidx * dxidx + dxidy * dxidy + dxidz * dxidz);

      const MFloat dxHelp = dxidx * gradxi;
      const MFloat dyHelp = dxidy * gradxi;
      const MFloat dzHelp = dxidz * gradxi;


      const MFloat cBC = sqrt(gamma * m_cells->pvariables[PV->P][cellIdG1] / m_cells->pvariables[PV->RHO][cellIdG1]);
      const MFloat rhoBC = m_cells->pvariables[PV->RHO][cellIdG1];

      const MFloat rhoInner = m_cells->pvariables[PV->RHO][cellIdA1];
      const MFloat uInner = m_cells->pvariables[PV->U][cellIdA1];
      const MFloat vInner = m_cells->pvariables[PV->V][cellIdA1];
      const MFloat wInner = m_cells->pvariables[PV->W][cellIdA1];
      const MFloat pInner = m_cells->pvariables[PV->P][cellIdA1];

      const MFloat maContravariant =
          (dxidx * uInner + dxidy * vInner + dxidz * wInner - m_cells->dxt[normalDir][cellIdA1]) * gradxi;

      if(maContravariant < F0) {
        // inflow
        const MFloat p = F1B2
                         * (pInner + PV->PInfinity
                            + rhoBC * cBC
                                  * (dxHelp * (uInner - PV->UInfinity) + dyHelp * (vInner - PV->VInfinity)
                                     + dzHelp * (wInner - PV->WInfinity)));

        const MFloat rho = CV->rhoInfinity + (p - PV->PInfinity) / POW2(cBC);
        const MFloat help = (p - PV->PInfinity) / (rhoBC * cBC);

        m_cells->pvariables[PV->RHO][cellIdG1] = rho;
        m_cells->pvariables[PV->U][cellIdG1] = (PV->UInfinity + help * dxHelp);
        m_cells->pvariables[PV->V][cellIdG1] = (PV->VInfinity + help * dyHelp);
        m_cells->pvariables[PV->W][cellIdG1] = (PV->WInfinity + help * dzHelp);
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
        m_cells->pvariables[PV->W][cellIdG1] = (wInner - help * dzHelp);
        m_cells->pvariables[PV->P][cellIdG1] = p;
        if(isRans) {
          m_cells->pvariables[PV->RANS_VAR[0]][cellIdG1] =
              (F2 * m_cells->pvariables[PV->RANS_VAR[0]][cellIdA1] - m_cells->pvariables[PV->RANS_VAR[0]][cellIdA2]);
        }
      }

      // extrapolate into second ghost cell
      for(MInt var = 0; var < PV->noVariables; var++) {
        m_cells->pvariables[var][cellIdG2] =
            F2 * m_cells->pvariables[var][cellIdG1] - m_cells->pvariables[var][cellIdA1];
      }
    }
  }
}


template <MBool isRans>
void StructuredBndryCnd3D<isRans>::bc2222(MInt bcId) {
  const MInt IJK[nDim] = {1, m_nCells[2], m_nCells[1] * m_nCells[2]};
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
  const MInt secondTangentialDir = (normalDir + 2) % nDim;
  const MInt normalDirStart = start[normalDir];
  const MInt firstTangentialStart = start[firstTangentialDir];
  const MInt firstTangentialEnd = end[firstTangentialDir];
  const MInt secondTangentialStart = start[secondTangentialDir];
  const MInt secondTangentialEnd = end[secondTangentialDir];
  const MInt inc[nDim] = {IJK[normalDir], IJK[firstTangentialDir], IJK[secondTangentialDir]};

  const MInt n = (face % 2) * 2 - 1;                                //-1,+1
  const MInt g1 = normalDirStart + (MInt)(0.5 - (0.5 * (MFloat)n)); //+1,0
  const MInt g2 = normalDirStart + (MInt)(0.5 + (0.5 * (MFloat)n)); // 0,+1
  const MInt a1 = normalDirStart + (MInt)(0.5 - (1.5 * (MFloat)n)); //+2,-1
  const MInt a2 = normalDirStart + (MInt)(0.5 - (2.5 * (MFloat)n)); //+3,-2

  for(MInt t1 = firstTangentialStart; t1 < firstTangentialEnd; t1++) {
    for(MInt t2 = secondTangentialStart; t2 < secondTangentialEnd; t2++) {
      const MInt cellIdG1 = g1 * inc[0] + t1 * inc[1] + t2 * inc[2];
      const MInt cellIdG2 = g2 * inc[0] + t1 * inc[1] + t2 * inc[2];
      const MInt cellIdA1 = a1 * inc[0] + t1 * inc[1] + t2 * inc[2];
      const MInt cellIdA2 = a2 * inc[0] + t1 * inc[1] + t2 * inc[2];

      const MFloat dxidx = m_cells->surfaceMetrics[normalDir * nDim + 0][cellIdA1];
      const MFloat dxidy = m_cells->surfaceMetrics[normalDir * nDim + 1][cellIdA1];
      const MFloat dxidz = m_cells->surfaceMetrics[normalDir * nDim + 2][cellIdA1];
      // multiply with n, so it will be -1 or +1 depending if we enter
      // or leave the domain of integration in positive direction
      const MFloat gradxi = n * F1 / sqrt(dxidx * dxidx + dxidy * dxidy + dxidz * dxidz);

      const MFloat dxHelp = dxidx * gradxi;
      const MFloat dyHelp = dxidy * gradxi;
      const MFloat dzHelp = dxidz * gradxi;


      const MFloat cBC = sqrt(gamma * m_cells->pvariables[PV->P][cellIdG1] / m_cells->pvariables[PV->RHO][cellIdG1]);
      const MFloat rhoBC = m_cells->pvariables[PV->RHO][cellIdG1];

      const MFloat rhoInner = m_cells->pvariables[PV->RHO][cellIdA1];
      const MFloat uInner = m_cells->pvariables[PV->U][cellIdA1];
      const MFloat vInner = m_cells->pvariables[PV->V][cellIdA1];
      const MFloat wInner = m_cells->pvariables[PV->W][cellIdA1];
      const MFloat pInner = m_cells->pvariables[PV->P][cellIdA1];

      const MFloat maContravariant =
          (dxidx * uInner + dxidy * vInner + dxidz * wInner - m_cells->dxt[normalDir][cellIdA1]) * gradxi;

      if(maContravariant < F0) {
        // inflow
        const MFloat uZonal = m_cells->fq[FQ->AVG_U][cellIdG1];
        const MFloat vZonal = m_cells->fq[FQ->AVG_V][cellIdG1];
        const MFloat wZonal = m_cells->fq[FQ->AVG_W][cellIdG1];
        const MFloat pZonal = m_cells->fq[FQ->AVG_P][cellIdG1];
        const MFloat rhoZonal = m_cells->fq[FQ->AVG_RHO][cellIdG1];
        const MFloat p =
            F1B2
            * (pInner + pZonal
               + rhoBC * cBC * (dxHelp * (uInner - uZonal) + dyHelp * (vInner - vZonal) + dzHelp * (wInner - wZonal)));

        const MFloat rho = rhoZonal + pZonal / POW2(cBC);
        const MFloat help = (p - pZonal) / (rhoBC * cBC);

        m_cells->pvariables[PV->RHO][cellIdG1] = rho;
        m_cells->pvariables[PV->U][cellIdG1] = (uZonal + help * dxHelp);
        m_cells->pvariables[PV->V][cellIdG1] = (vZonal + help * dyHelp);
        m_cells->pvariables[PV->W][cellIdG1] = (wZonal + help * dzHelp);
        m_cells->pvariables[PV->P][cellIdG1] = p;

        if(isRans) {
          m_cells->pvariables[PV->RANS_VAR[0]][cellIdG1] = m_cells->pvariables[PV->RANS_VAR[0]][cellIdA1];
        }
      } else {
        // outflow
        const MFloat pZonal = m_cells->fq[FQ->AVG_P][cellIdG1];
        const MFloat rho = rhoInner + (pZonal - pInner) / POW2(cBC);
        const MFloat help = (pZonal - pInner) / (rhoBC * cBC);

        m_cells->pvariables[PV->RHO][cellIdG1] = rho;
        m_cells->pvariables[PV->U][cellIdG1] = (uInner - help * dxHelp);
        m_cells->pvariables[PV->V][cellIdG1] = (vInner - help * dyHelp);
        m_cells->pvariables[PV->W][cellIdG1] = (wInner - help * dzHelp);
        m_cells->pvariables[PV->P][cellIdG1] = pZonal;
        if(isRans) {
          m_cells->pvariables[PV->RANS_VAR[0]][cellIdG1] =
              (F2 * m_cells->pvariables[PV->RANS_VAR[0]][cellIdA1] - m_cells->pvariables[PV->RANS_VAR[0]][cellIdA2]);
        }
      }

      // extrapolate into second ghost cell
      for(MInt var = 0; var < PV->noVariables; var++) {
        m_cells->pvariables[var][cellIdG2] =
            F2 * m_cells->pvariables[var][cellIdG1] - m_cells->pvariables[var][cellIdA1];
      }
    }
  }
}


/**
 * \brief Outflow condition after shock
 *
 *  Uses characteristic subsonic inflow (only correct if the machnumber normal to the Boundary is smaller than 1)
 *  \authors Simon Loosen, Pascal Meysonnat
 */
template <MBool isRans>
void StructuredBndryCnd3D<isRans>::bc2009(MInt bcId) {
  const MFloat gamma = m_solver->m_gamma;
  const MFloat gammaMinusOne = gamma - 1.0;
  const MFloat gammaPlusOne = gamma + 1.0;
  MInt* start = m_physicalBCMap[bcId]->start1;
  MInt* end = m_physicalBCMap[bcId]->end1;
  MFloat rho2 = F0, u2 = F0, v2tmp = F0, p2 = F0; // Conditions behind the shock
  MFloat u1n = F0, v1n = F0, u2n = F0,
         v2n = F0;             // Velocities normal to the shock upstream and downstream of the shock
  MFloat rho21 = F0, p21 = F0; // quotient of the values upstream and downstream of the shock
  MFloat beta = F0, M_mean = m_solver->m_Ma;

  // pressure and density ratio over the shock
  rho21 = gammaPlusOne * POW2(M_mean * sin(m_sigma)) / (gammaMinusOne * POW2(M_mean * sin(m_sigma)) + 2);
  p21 = 1.0 + 2 * gamma / gammaPlusOne * (POW2(M_mean * sin(m_sigma)) - 1);


  // normal and parallel velocities upstream of the shock
  u1n = M_mean * sqrt(PV->TInfinity) * sin(m_sigma);
  v1n = M_mean * sqrt(PV->TInfinity) * cos(m_sigma);
  // normal and parallel velocity downstream of the shock
  u2n = u1n / rho21;
  v2n = v1n;

  beta = m_sigma - atan(u2n / v2n); // redirection angle of the flow

  u2 = sqrt(POW2(u2n) + POW2(v2n)) * cos(beta);
  v2tmp = sqrt(POW2(u2n) + POW2(v2n)) * sin(beta);

  rho2 = rho21 * CV->rhoInfinity;
  p2 = p21 * PV->PInfinity;

  switch(m_physicalBCMap[bcId]->face) {
    case 3: {
      MInt cellId = -1;
      MInt cellIdadj = -1;
      MInt pIJK = 0, pIJPK = 0, pIJKP = 0, pIJPKP = 0;
      // fully new bc
      MInt j = start[1];
      MFloat pBC = F0, rho = F0, u = F0, v = F0, w = F0;
      MFloat drho = F0, du = F0, dv = F0, dw = F0, dp = F0;
      MFloat yBC = F0;
      MFloat v2 = F0;
      // pBC=PV->PInfinity;
      MFloat pInner = F0, c02 = F0, distance = F0;

      // Change sign of the v velocity
      v2 = -v2tmp;

      for(MInt k = start[2]; k < end[2]; k++) {
        for(MInt i = start[0]; i < end[0]; i++) {
          // cellId=cellIndex(i,j,k);
          cellIdadj = cellIndex(i, j - 1, k);
          // to determine the face coordinates!!!!!!
          pIJK = getPointIdFromCell(i, j, k);
          pIJPK = getPointIdfromPoint(pIJK, 1, 0, 0);
          pIJKP = getPointIdfromPoint(pIJK, 0, 0, 1);
          pIJPKP = getPointIdfromPoint(pIJK, 1, 0, 1);

          // values at the inner point
          pInner = m_cells->pvariables[PV->P][cellIdadj];
          c02 = sqrt(gamma * pInner / m_cells->pvariables[PV->RHO][cellIdadj]);
          u = m_cells->pvariables[PV->U][cellIdadj];
          v = m_cells->pvariables[PV->V][cellIdadj];
          w = m_cells->pvariables[PV->W][cellIdadj];

          MFloat dxidx = m_cells->surfaceMetrics[3][cellIdadj];
          MFloat dxidy = m_cells->surfaceMetrics[4][cellIdadj];
          MFloat dxidz = m_cells->surfaceMetrics[5][cellIdadj];

          // leaving domain of integration in positive coordinate direction,
          // therefore multiply with positive F1
          MFloat gradxi = F1 / sqrt(dxidx * dxidx + dxidy * dxidy + dxidz * dxidz);

          MFloat dxHelp = dxidx * gradxi;
          MFloat dyHelp = dxidy * gradxi;
          MFloat dzHelp = dxidz * gradxi;

          // values at the boundary
          pBC = F1B2
                * (pInner + p2
                   + m_cells->pvariables[PV->RHO][cellIdadj] * c02
                         * (dxHelp * (u - u2) + dyHelp * (v - v2) + dzHelp * (w - PV->WInfinity)));
          rho = rho2 + ((pBC - p2) / (c02 * c02));

          u = u2 + dxHelp * (pBC - p2) / (m_cells->pvariables[PV->RHO][cellIdadj] * c02);
          v = v2 + dyHelp * (pBC - p2) / (m_cells->pvariables[PV->RHO][cellIdadj] * c02);
          w = PV->WInfinity + dzHelp * (pBC - p2) / (m_cells->pvariables[PV->RHO][cellIdadj] * c02);

          // extrapolate the variables into the ghost cells
          // gradients
          yBC = F1B4
                * (m_grid->m_coordinates[1][pIJK] + m_grid->m_coordinates[1][pIJPK] + m_grid->m_coordinates[1][pIJKP]
                   + m_grid->m_coordinates[1][pIJPKP]);
          distance = (yBC - m_cells->coordinates[1][cellIdadj]);

          drho = (rho - m_cells->pvariables[PV->RHO][cellIdadj]) / distance;
          du = (u - m_cells->pvariables[PV->U][cellIdadj]) / distance;
          dv = (v - m_cells->pvariables[PV->V][cellIdadj]) / distance;
          dw = (w - m_cells->pvariables[PV->W][cellIdadj]) / distance;
          dp = (pBC - m_cells->pvariables[PV->P][cellIdadj]) / distance;

          // extrapolate:
          for(MInt jj = start[1]; jj < end[1]; ++jj) {
            cellId = cellIndex(i, jj, k);
            distance = (m_cells->coordinates[1][cellId] - m_cells->coordinates[1][cellIdadj]);
            m_cells->pvariables[PV->RHO][cellId] = m_cells->pvariables[PV->RHO][cellIdadj] + drho * distance;
            m_cells->pvariables[PV->U][cellId] = m_cells->pvariables[PV->U][cellIdadj] + du * distance;
            m_cells->pvariables[PV->V][cellId] = m_cells->pvariables[PV->V][cellIdadj] + dv * distance;
            m_cells->pvariables[PV->W][cellId] = m_cells->pvariables[PV->W][cellIdadj] + dw * distance;
            m_cells->pvariables[PV->P][cellId] = m_cells->pvariables[PV->P][cellIdadj] + dp * distance;
          }
        }
      }
      break;
    }
    default: {
      mTerm(1, AT_, "BC-face not implemented");
    }
  }
}


template <MBool isRans>
void StructuredBndryCnd3D<isRans>::bc3000(MInt bcId) {
  const MInt* start = m_physicalBCMap[bcId]->start1;
  const MInt* end = m_physicalBCMap[bcId]->end1;
  const MInt face = m_physicalBCMap[bcId]->face;

  // const MInt dim = face/2;
  // const MInt n1m1[9] = {-1,1,1, 1,-1,1, 1,1,-1};

  for(MInt k = start[2]; k < end[2]; k++) {
    for(MInt j = end[1] - 1; j >= start[1]; j--) {
      for(MInt i = start[0]; i < end[0]; i++) {
        MInt cellId = -1;
        MInt cellIdAdj = -1;
        tie(cellId, cellIdAdj) = getMirrorCellIdPair(i, j, k, face);

        m_cells->pvariables[PV->RHO][cellId] = m_cells->pvariables[PV->RHO][cellIdAdj];
        m_cells->pvariables[PV->P][cellId] = m_cells->pvariables[PV->P][cellIdAdj];
        m_cells->pvariables[PV->U][cellId] = (1.0) * m_cells->pvariables[PV->U][cellIdAdj];
        m_cells->pvariables[PV->V][cellId] = (-1.0) * m_cells->pvariables[PV->V][cellIdAdj];
        m_cells->pvariables[PV->W][cellId] = (1.0) * m_cells->pvariables[PV->W][cellIdAdj];
        if(isRans) {
          m_cells->pvariables[PV->RANS_VAR[0]][cellId] = m_cells->pvariables[PV->RANS_VAR[0]][cellIdAdj];
        }
      }
    }
  }
}

/**
 * \brief Streamline symmetry
 *
 * \author Pascal S. Meysonnat, Volkan
 * \date   August 2018
 */
template <MBool isRans>
void StructuredBndryCnd3D<isRans>::bc3001(MInt bcId) {
  (void)bcId;
}

/**
 * \brief Falkner-Skan-Cooke inflow boundary condition
 *
 * Prescribes the velocity distribution for a Falkner-Skan-Cooke
 * boundary layer at the inflow, with the FSC helper functions
 *
 * \authors Thomas Schilden, Marian Albers
 */
template <MBool isRans>
void StructuredBndryCnd3D<isRans>::bc2888(MInt bcId) {
  const MInt IJK[nDim] = {1, m_nCells[2], m_nCells[1] * m_nCells[2]};
  MInt* start = m_physicalBCMap[bcId]->start1;
  MInt* end = m_physicalBCMap[bcId]->end1;

  // Here we find out the normal direction of the
  // boundary and the two tangential directions.
  // This way we can make a general formulation of
  // the boundary condition
  const MInt face = m_physicalBCMap[bcId]->face;
  const MInt normalDir = face / 2;
  const MInt firstTangentialDir = (normalDir + 1) % nDim;
  const MInt secondTangentialDir = (normalDir + 2) % nDim;
  const MInt normalDirStart = start[normalDir];
  const MInt firstTangentialStart = start[firstTangentialDir];
  const MInt firstTangentialEnd = end[firstTangentialDir];
  const MInt secondTangentialStart = start[secondTangentialDir];
  const MInt secondTangentialEnd = end[secondTangentialDir];
  const MInt inc[nDim] = {IJK[normalDir], IJK[firstTangentialDir], IJK[secondTangentialDir]};

  const MInt n = (face % 2) * 2 - 1;                                //-1,+1
  const MInt g1 = normalDirStart + (MInt)(0.5 - (0.5 * (MFloat)n)); //+1,0
  const MInt g2 = normalDirStart + (MInt)(0.5 + (0.5 * (MFloat)n)); // 0,+1
  const MInt a1 = normalDirStart + (MInt)(0.5 - (1.5 * (MFloat)n)); //+2,-1
  // const MInt a2 = normalDirStart + (MInt)(0.5-(2.5*(MFloat)n)); //+3,-2

  for(MInt t1 = firstTangentialStart; t1 < firstTangentialEnd; t1++) {
    for(MInt t2 = secondTangentialStart; t2 < secondTangentialEnd; t2++) {
      const MInt cellIdG1 = g1 * inc[0] + t1 * inc[1] + t2 * inc[2];
      const MInt cellIdG2 = g2 * inc[0] + t1 * inc[1] + t2 * inc[2];
      const MInt cellIdA1 = a1 * inc[0] + t1 * inc[1] + t2 * inc[2];
      // const MInt cellIdA2 = a2*inc[0] + t1*inc[1] + t2*inc[2];

      MFloat vel[nDim];
      m_solver->getFscVelocity(cellIdG1, vel);

      m_cells->pvariables[PV->U][cellIdG1] = vel[0];
      m_cells->pvariables[PV->V][cellIdG1] = vel[1];
      m_cells->pvariables[PV->W][cellIdG1] = vel[2];

      m_cells->pvariables[PV->P][cellIdG1] = m_cells->pvariables[PV->P][cellIdA1];
      m_cells->pvariables[PV->RHO][cellIdG1] = CV->rhoInfinity;

      m_solver->getFscVelocity(cellIdG2, vel);
      m_cells->pvariables[PV->U][cellIdG2] = vel[0];
      m_cells->pvariables[PV->V][cellIdG2] = vel[1];
      m_cells->pvariables[PV->W][cellIdG2] = vel[2];
      m_cells->pvariables[PV->P][cellIdG1] =
          F2 * m_cells->pvariables[PV->P][cellIdG1] - m_cells->pvariables[PV->P][cellIdA1];
      m_cells->pvariables[PV->RHO][cellIdG2] = CV->rhoInfinity;
    }
  }
}

/**
 * \brief Simple outflow with pressure gradient for FSC boundary layer
 *  \authors Thomas Schilden, Marian Albers
 */
template <MBool isRans>
void StructuredBndryCnd3D<isRans>::bc2730(MInt bcId) {
  const MInt IJK[nDim] = {1, m_nCells[2], m_nCells[1] * m_nCells[2]};
  MInt* start = m_physicalBCMap[bcId]->start1;
  MInt* end = m_physicalBCMap[bcId]->end1;

  // Here we find out the normal direction of the
  // boundary and the two tangential directions.
  // This way we can make a general formulation of
  // the boundary condition
  const MInt face = m_physicalBCMap[bcId]->face;
  const MInt normalDir = face / 2;
  const MInt firstTangentialDir = (normalDir + 1) % nDim;
  const MInt secondTangentialDir = (normalDir + 2) % nDim;
  const MInt normalDirStart = start[normalDir];
  const MInt firstTangentialStart = start[firstTangentialDir];
  const MInt firstTangentialEnd = end[firstTangentialDir];
  const MInt secondTangentialStart = start[secondTangentialDir];
  const MInt secondTangentialEnd = end[secondTangentialDir];
  const MInt inc[nDim] = {IJK[normalDir], IJK[firstTangentialDir], IJK[secondTangentialDir]};

  const MInt n = (face % 2) * 2 - 1;                                //-1,+1
  const MInt g1 = normalDirStart + (MInt)(0.5 - (0.5 * (MFloat)n)); //+1,0
  const MInt g2 = normalDirStart + (MInt)(0.5 + (0.5 * (MFloat)n)); // 0,+1
  const MInt a1 = normalDirStart + (MInt)(0.5 - (1.5 * (MFloat)n)); //+2,-1

  for(MInt t1 = firstTangentialStart; t1 < firstTangentialEnd; t1++) {
    for(MInt t2 = secondTangentialStart; t2 < secondTangentialEnd; t2++) {
      const MInt cellIdG1 = g1 * inc[0] + t1 * inc[1] + t2 * inc[2];
      const MInt cellIdG2 = g2 * inc[0] + t1 * inc[1] + t2 * inc[2];
      const MInt cellIdA1 = a1 * inc[0] + t1 * inc[1] + t2 * inc[2];

      const MFloat rhoInner = m_cells->pvariables[PV->RHO][cellIdA1];
      const MFloat uInner = m_cells->pvariables[PV->U][cellIdA1];
      const MFloat vInner = m_cells->pvariables[PV->V][cellIdA1];
      const MFloat wInner = m_cells->pvariables[PV->W][cellIdA1];

      m_cells->pvariables[PV->U][cellIdG1] = uInner;
      m_cells->pvariables[PV->V][cellIdG1] = vInner;
      m_cells->pvariables[PV->W][cellIdG1] = wInner;
      m_cells->pvariables[PV->P][cellIdG1] = m_solver->getFscPressure(cellIdG1);
      m_cells->pvariables[PV->RHO][cellIdG1] = rhoInner;

      // extrapolate into second ghost cell
      for(MInt var = 0; var < PV->noVariables; var++) {
        m_cells->pvariables[var][cellIdG2] =
            F2 * m_cells->pvariables[var][cellIdG1] - m_cells->pvariables[var][cellIdA1];
      }
    }
  }
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
void StructuredBndryCnd3D<isRans>::bc2999(MInt bcId) {
  const MInt IJK[nDim] = {1, m_nCells[2], m_nCells[1] * m_nCells[2]};
  MInt* start = m_physicalBCMap[bcId]->start1;
  MInt* end = m_physicalBCMap[bcId]->end1;

  // Here we find out the normal direction of the
  // boundary and the two tangential directions.
  // This way we can make a general formulation of
  // the boundary condition
  const MInt face = m_physicalBCMap[bcId]->face;
  const MInt normalDir = face / 2;
  const MInt firstTangentialDir = (normalDir + 1) % nDim;
  const MInt secondTangentialDir = (normalDir + 2) % nDim;
  const MInt normalDirStart = start[normalDir];
  const MInt firstTangentialStart = start[firstTangentialDir];
  const MInt firstTangentialEnd = end[firstTangentialDir];
  const MInt secondTangentialStart = start[secondTangentialDir];
  const MInt secondTangentialEnd = end[secondTangentialDir];
  const MInt inc[nDim] = {IJK[normalDir], IJK[firstTangentialDir], IJK[secondTangentialDir]};

  const MInt n = (face % 2) * 2 - 1;                                //-1,+1
  const MInt g1 = normalDirStart + (MInt)(0.5 - (0.5 * (MFloat)n)); //+1,0
  const MInt g2 = normalDirStart + (MInt)(0.5 + (0.5 * (MFloat)n)); // 0,+1
  const MInt a1 = normalDirStart + (MInt)(0.5 - (1.5 * (MFloat)n)); //+2,-1
  // const MInt a2 = normalDirStart + (MInt)(0.5-(2.5*(MFloat)n)); //+3,-2

  for(MInt t1 = firstTangentialStart; t1 < firstTangentialEnd; t1++) {
    for(MInt t2 = secondTangentialStart; t2 < secondTangentialEnd; t2++) {
      const MInt cellIdG1 = g1 * inc[0] + t1 * inc[1] + t2 * inc[2];
      const MInt cellIdG2 = g2 * inc[0] + t1 * inc[1] + t2 * inc[2];
      const MInt cellIdA1 = a1 * inc[0] + t1 * inc[1] + t2 * inc[2];
      // const MInt cellIdA2 = a2*inc[0] + t1*inc[1] + t2*inc[2];

      MFloat vel[nDim];
      m_solver->getBlasiusVelocity(cellIdG1, vel);

      m_cells->pvariables[PV->U][cellIdG1] = vel[0];
      m_cells->pvariables[PV->V][cellIdG1] = vel[1];
      m_cells->pvariables[PV->W][cellIdG1] = vel[2];

      m_cells->pvariables[PV->P][cellIdG1] = m_cells->pvariables[PV->P][cellIdA1];
      m_cells->pvariables[PV->RHO][cellIdG1] = CV->rhoInfinity;

      m_solver->getBlasiusVelocity(cellIdG2, vel);
      m_cells->pvariables[PV->U][cellIdG2] = vel[0];
      m_cells->pvariables[PV->V][cellIdG2] = vel[1];
      m_cells->pvariables[PV->W][cellIdG2] = vel[2];
      m_cells->pvariables[PV->P][cellIdG2] =
          F2 * m_cells->pvariables[PV->P][cellIdG1] - m_cells->pvariables[PV->P][cellIdA1];
      m_cells->pvariables[PV->RHO][cellIdG2] = CV->rhoInfinity;
    }
  }
}

/**
 * \brief Reformulated Synthetic Turbulence Generation
 *
 * Synthetic Turbulence Generation Method
 * Computes fluctuations from a given averaged
 * velocity and nu_t profile and adds them at
 * the inflow of the domain.
 * Same computations as STG in TFS by Benedikt Roidl
 *
 * \author Marian Albers
 */
template <MBool isRans>
void StructuredBndryCnd3D<isRans>::bc7909(MInt bcId) {
  const MInt NU_T = 5;
  const MInt SIJSIJ = 6;
  const MInt LENGTH_SCALE = 7;
  const MInt FLUC_UU = 8;
  const MInt FLUC_VV = 9;
  const MInt FLUC_WW = 10;
  const MInt FLUC_UV = 11;
  const MInt FLUC_VW = 12;
  const MInt FLUC_UW = 13;
  const MInt FLUC_U = 14;
  const MInt FLUC_V = 15;
  const MInt FLUC_W = 16;
  const MInt LENGTH_X = 17;
  const MInt LENGTH_Y = 18;
  const MInt LENGTH_Z = 19;

  const MFloat gamma = m_solver->m_gamma;
  const MFloat gammaMinusOne = gamma - 1.0;

  MInt* start = m_physicalBCMap[bcId]->start1;
  MInt* end = m_physicalBCMap[bcId]->end1;

  // Initalization of variables that will be used
  MFloat SijSij, Vb;
  MInt I, IBC;
  const MInt ii = 1;
  const MInt nran = m_solver->m_stgMaxNoEddies;

  // Arrays for MPI operations
  MFloatScratchSpace maxValsLocal(4, AT_, "maxValsLocal");
  MFloatScratchSpace maxValsGlobal(4, AT_, "maxValsGlobal");
  maxValsLocal.fill(F0);
  maxValsGlobal.fill(F0);

  MFloatScratchSpace eddyBcastBuffer(nran * m_solver->m_stgNoEddieProperties, AT_, "eddyBcastBuffer");
  eddyBcastBuffer.fill(F0);

  /////////////////////////////////////////////
  // Value initialization for some variables
  /////////////////////////////////////////////

  I = cellIndex(start[0], start[1], start[2]);
  const MFloat fre = 1.0 / m_solver->m_Re0;

  // Limiters
  const MFloat epss = 1e-34;
  const MFloat eps = 1e-16;
  const MFloat epsl = 1e-13;

  // Box volume must be specified manually and is sign sensitive
  const MFloat delta_in = m_solver->m_stgDelta99Inflow;

  const MFloat BLT1 = abs(m_stgVbEnd[0] - m_stgVbStart[0]);
  const MFloat BLT2 = abs(m_stgVbEnd[1] - m_stgVbStart[1]);
  const MFloat BLT3 = abs(m_stgVbEnd[2] - m_stgVbStart[2]);

  // Manual parameters
  const MFloat c_mu = 0.09;
  const MFloat a1 = 1 / sqrt(c_mu);

  const MFloat timsm = 0.3;                  // Smoothing factor for new eddies (in time)
  const MFloat aniso = 1.0;                  // Anisotropic, clustered eddies
  const MFloat exple = m_solver->m_stgExple; // to scale length of eddies

  Vb = 0.0; // Initial box volume

  MInt I1 = cellIndexBC(ii, start[1] + m_noGhostLayers, start[2] + m_noGhostLayers);

  MFloat rhoRANSI1 = m_cells->stg_fq[PV->RHO][I1];
  MFloat pressure1 = m_cells->stg_fq[PV->P][I1];

  const MFloat temp = m_solver->m_gamma * pressure1 / rhoRANSI1;
  const MFloat xmu = SUTHERLANDLAW(temp);

  // get delta0 boundary layer thickness:
  if(globalTimeStep % 250 == 0 && m_solver->m_RKStep == 0) {
    MInt bndryLayerIndex = 0;
    MFloat u99_0 = F0, u99_1 = F0;
    for(MInt j = start[1] + m_noGhostLayers; j < end[1] - m_noGhostLayers; j++) {
      if(m_cells->stg_fq[PV->U][cellIndexBC(0, j, 2)] >= 0.99 * PV->UInfinity
         && m_cells->stg_fq[PV->U][cellIndexBC(0, j - 1, 2)] < 0.99 * PV->UInfinity) {
        u99_0 = m_cells->stg_fq[PV->U][cellIndexBC(0, j - 1, 2)];
        u99_1 = m_cells->stg_fq[PV->U][cellIndexBC(0, j, 2)];
        bndryLayerIndex = j - 1;
        break;
      }
    }
    MFloat bndryLayerThicknessLocal = 0.0;
    MFloat bndryLayerThicknessGlobal = 0.0;

    if(bndryLayerIndex > 0) {
      bndryLayerThicknessLocal = m_cells->coordinates[1][cellIndex(0, bndryLayerIndex, 2)]
                                 + (0.99 * PV->UInfinity - u99_0) / (u99_1 - u99_0)
                                       * (m_cells->coordinates[1][cellIndex(0, bndryLayerIndex + 1, 2)]
                                          - m_cells->coordinates[1][cellIndex(0, bndryLayerIndex, 2)]);
    }
    MPI_Allreduce(&bndryLayerThicknessLocal, &bndryLayerThicknessGlobal, 1, MPI_DOUBLE, MPI_MAX, m_solver->m_commStg,
                  AT_, "bndryLayerThicknessLocal", "bndryLayerThicknessGlobal");
    if(m_solver->m_stgMyRank == m_solver->m_commStgRoot) {
      cout << "Boundary Layer thickness delta99: " << bndryLayerThicknessGlobal << endl;
    }
  }

  // this solver is only executed at the beginning as the RANS profile doesn't vary
  if((globalTimeStep <= m_solver->m_restartTimeStep)
     || ((m_solver->m_restartTimeStep == 0 && globalTimeStep < 2) && (m_solver->m_RKStep == 0))
     || (m_solver->m_zonal && (globalTimeStep % m_solver->m_zonalExchangeInterval == 0) && (m_solver->m_RKStep == 0))) {
    // Initialize max values
    MFloat utaumax = F0, umax = F0, vmax = F0, wmax = F0, minLengthLocal = F0;

    if(m_solver->m_stgMyRank == m_solver->m_commStgRoot) {
      cout << "STG - Computing Reynolds Tensor Components and Length Scales..." << endl;
    }

    for(MInt k = start[2] + 1; k < end[2] - 1; k++) {
      for(MInt j = start[1] + 1; j < end[1] - 1; j++) {
        ////////////////////////////////////////////////////////
        // This is the metrics / strain tensor / max shear part
        ////////////////////////////////////////////////////////

        MInt IPJK, IMJK, IJPK, IJMK, IJKP, IJKM;

        I = cellIndex(ii, j, k);
        IPJK = cellIndex(ii + 1, j, k);
        IMJK = cellIndex(ii - 1, j, k);
        IJPK = cellIndex(ii, j + 1, k);
        IJMK = cellIndex(ii, j - 1, k);
        IJKP = cellIndex(ii, j, k + 1);
        IJKM = cellIndex(ii, j, k - 1);


        const MFloat dxdi = F1B2 * (m_cells->coordinates[0][IPJK] - m_cells->coordinates[0][IMJK]);
        const MFloat dxdj = F1B2 * (m_cells->coordinates[0][IJPK] - m_cells->coordinates[0][IJMK]);
        const MFloat dxdk = F1B2 * (m_cells->coordinates[0][IJKP] - m_cells->coordinates[0][IJKM]);
        const MFloat dydi = F1B2 * (m_cells->coordinates[1][IPJK] - m_cells->coordinates[1][IMJK]);
        const MFloat dydj = F1B2 * (m_cells->coordinates[1][IJPK] - m_cells->coordinates[1][IJMK]);
        const MFloat dydk = F1B2 * (m_cells->coordinates[1][IJKP] - m_cells->coordinates[1][IJKM]);
        const MFloat dzdi = F1B2 * (m_cells->coordinates[2][IPJK] - m_cells->coordinates[2][IMJK]);
        const MFloat dzdj = F1B2 * (m_cells->coordinates[2][IJPK] - m_cells->coordinates[2][IJMK]);
        const MFloat dzdk = F1B2 * (m_cells->coordinates[2][IJKP] - m_cells->coordinates[2][IJKM]);

        const MFloat dxl = sqrt(dxdi * dxdi + dydi * dydi + dzdi * dzdi);
        const MFloat dyl = sqrt(dxdj * dxdj + dydj * dydj + dzdj * dzdj);
        const MFloat dzl = sqrt(dxdk * dxdk + dydk * dydk + dzdk * dzdk);

        const MFloat dxidx = (1. / max(m_cells->cellJac[I], epss)) * m_cells->cellMetrics[0][I];
        const MFloat dxidy = (1. / max(m_cells->cellJac[I], epss)) * m_cells->cellMetrics[1][I];
        const MFloat dxidz = (1. / max(m_cells->cellJac[I], epss)) * m_cells->cellMetrics[2][I];

        const MFloat detadx = (1. / max(m_cells->cellJac[I], epss)) * m_cells->cellMetrics[3 + 0][I];
        const MFloat detady = (1. / max(m_cells->cellJac[I], epss)) * m_cells->cellMetrics[3 + 1][I];
        const MFloat detadz = (1. / max(m_cells->cellJac[I], epss)) * m_cells->cellMetrics[3 + 2][I];

        const MFloat dzetadx = (1. / max(m_cells->cellJac[I], epss)) * m_cells->cellMetrics[6 + 0][I];
        const MFloat dzetady = (1. / max(m_cells->cellJac[I], epss)) * m_cells->cellMetrics[6 + 1][I];
        const MFloat dzetadz = (1. / max(m_cells->cellJac[I], epss)) * m_cells->cellMetrics[6 + 2][I];

        IBC = cellIndexBC(ii, j, k);
        IPJK = cellIndexBC(ii + 1, j, k);
        IMJK = cellIndexBC(ii - 1, j, k);
        IJPK = cellIndexBC(ii, j + 1, k);
        IJMK = cellIndexBC(ii, j - 1, k);
        IJKP = cellIndexBC(ii, j, k + 1);
        IJKM = cellIndexBC(ii, j, k - 1);

        const MFloat frho = 1.0 / m_cells->stg_fq[PV->RHO][IBC];

        // dud?
        const MFloat dudxi = F1B2 * (m_cells->stg_fq[PV->U][IPJK] - m_cells->stg_fq[PV->U][IMJK]);
        const MFloat dudeta = F1B2 * (m_cells->stg_fq[PV->U][IJPK] - m_cells->stg_fq[PV->U][IJMK]);
        const MFloat dudzeta = F1B2 * (m_cells->stg_fq[PV->U][IJKP] - m_cells->stg_fq[PV->U][IJKM]);

        // dvd?
        const MFloat dvdxi = F1B2 * (m_cells->stg_fq[PV->V][IPJK] - m_cells->stg_fq[PV->V][IMJK]);
        const MFloat dvdeta = F1B2 * (m_cells->stg_fq[PV->V][IJPK] - m_cells->stg_fq[PV->V][IJMK]);
        const MFloat dvdzeta = F1B2 * (m_cells->stg_fq[PV->V][IJKP] - m_cells->stg_fq[PV->V][IJKM]);

        // dwd?
        const MFloat dwdxi = F1B2 * (m_cells->stg_fq[PV->W][IPJK] - m_cells->stg_fq[PV->W][IMJK]);
        const MFloat dwdeta = F1B2 * (m_cells->stg_fq[PV->W][IJPK] - m_cells->stg_fq[PV->W][IJMK]);
        const MFloat dwdzeta = F1B2 * (m_cells->stg_fq[PV->W][IJKP] - m_cells->stg_fq[PV->W][IJKM]);

        const MFloat dudx = dudxi * dxidx + dudeta * detadx + dudzeta * dzetadx;
        const MFloat dudy = dudxi * dxidy + dudeta * detady + dudzeta * dzetady;
        const MFloat dudz = dudxi * dxidz + dudeta * detadz + dudzeta * dzetadz;

        const MFloat dvdx = dvdxi * dxidx + dvdeta * detadx + dvdzeta * dzetadx;
        const MFloat dvdy = dvdxi * dxidy + dvdeta * detady + dvdzeta * dzetady;
        const MFloat dvdz = dvdxi * dxidz + dvdeta * detadz + dvdzeta * dzetadz;

        const MFloat dwdx = dwdxi * dxidx + dwdeta * detadx + dwdzeta * dzetadx;
        const MFloat dwdy = dwdxi * dxidy + dwdeta * detady + dwdzeta * dzetady;
        const MFloat dwdz = dwdxi * dxidz + dwdeta * detadz + dwdzeta * dzetadz;


        const MFloat s11 = 2.0 * dudx;
        const MFloat s12 = dvdx + dudy;
        const MFloat s13 = dwdx + dudz;

        const MFloat s21 = dudy + dvdx;
        const MFloat s22 = 2.0 * dvdy;
        const MFloat s23 = dwdy + dvdz;

        const MFloat s31 = dudz + dwdx;
        const MFloat s32 = dvdz + dwdy;
        const MFloat s33 = 2.0 * dwdz;

        // Strain tensor
        SijSij = F1B4
                 * (s11 * s11 + s12 * s12 + s13 * s13 + s21 * s21 + s22 * s22 + s23 * s23 + s31 * s31 + s32 * s32
                    + s33 * s33);

        if(std::isnan(SijSij)) {
          cout << " dudxi  : " << dudxi << " dudeta : " << dudeta << " dudzeta: " << dudzeta << " dvdxi  : " << dvdxi
               << " dvdeta : " << dvdeta << " dvdzeta: " << dvdzeta << " dwdxi  : " << dwdxi << " dwdeta : " << dwdeta
               << " dwdzeta: " << dwdzeta << endl;
        }


        //>marian: in TFS code this isn't SQRT
        m_cells->stg_fq[SIJSIJ][cellIndexBC(ii, j, k)] = SijSij;

        // Assume a one-, or two equation turbulence model
        // Assume a simplified directivity:

        // Read from RANS profile
        MFloat nu_t = m_cells->stg_fq[NU_T][cellIndexBC(m_noGhostLayers, j, k)];


        const MFloat sr1 = (s12 + s21) * (s12 + s21);
        const MFloat sr2 = (s23 + s32) * (s23 + s32);
        const MFloat sr3 = (s13 + s31) * (s13 + s31);
        const MFloat srt = max(sqrt(sr1 + sr2 + sr3), epsl);

        const MFloat rr1 = sqrt(sr1) / srt;
        const MFloat rr2 = sqrt(sr2) / srt;
        const MFloat rr3 = sqrt(sr3) / srt;

        const MFloat uv = -sqrt(2.0 * SijSij) * rr1 * nu_t * fre;
        const MFloat vw = -sqrt(2.0 * SijSij) * rr2 * nu_t * fre;
        const MFloat uw = -sqrt(2.0 * SijSij) * rr3 * nu_t * fre;
        const MFloat uu = a1 * abs(uv) * m_solver->m_stgRSTFactors[0];
        const MFloat vv = a1 * abs(uv) * m_solver->m_stgRSTFactors[1];
        const MFloat ww = a1 * abs(uv) * m_solver->m_stgRSTFactors[2];

        m_cells->stg_fq[FLUC_UU][cellIndexBC(ii, j, k)] = uu;
        m_cells->stg_fq[FLUC_VV][cellIndexBC(ii, j, k)] = vv;
        m_cells->stg_fq[FLUC_WW][cellIndexBC(ii, j, k)] = ww;
        m_cells->stg_fq[FLUC_UV][cellIndexBC(ii, j, k)] = uv;
        m_cells->stg_fq[FLUC_VW][cellIndexBC(ii, j, k)] = vw;
        m_cells->stg_fq[FLUC_UW][cellIndexBC(ii, j, k)] = uw;

        // Get utau using laminar viscosity
        const MFloat utau2 = sqrt(fre * sqrt(2.0 * SijSij) * xmu * frho);

        // Save values if they are the new maximum
        if(utau2 >= utaumax) {
          minLengthLocal = pow((dxl * dyl * dzl), 0.33);
          utaumax = max(utau2, utaumax);
        }

        // In which direction aims the maximum averaged velocity?
        const MFloat u = m_cells->stg_fq[PV->U][IBC];
        const MFloat v = m_cells->stg_fq[PV->V][IBC];
        const MFloat w = m_cells->stg_fq[PV->W][IBC];

        umax = max(u, umax);
        vmax = max(v, vmax);
        wmax = max(w, wmax);

        // We need a global length scale to compare
        m_cells->stg_fq[LENGTH_SCALE][cellIndexBC(ii, j, k)] = pow(sqrt(max(2.0 * SijSij, epss)), -exple);

        //////////////////////////////////////////////////
        // Zero gradient extrapolation of boundary values
        //////////////////////////////////////////////////
        MInt noSTGVariables = 14;

        if(k == start[2] + 1) {
          // 1st layer
          for(MInt var = 6; var < noSTGVariables; var++) {
            m_cells->stg_fq[var][cellIndexBC(ii, j, k - 1)] = m_cells->stg_fq[var][IBC];
          }
        } else if(k == end[2] - m_noGhostLayers) {
          // 1st layer
          for(MInt var = 6; var < noSTGVariables; var++) {
            m_cells->stg_fq[var][cellIndexBC(ii, j, k + 1)] = m_cells->stg_fq[var][IBC];
          }
        }

        if(j == start[1] + 1) {
          for(MInt var = 6; var < noSTGVariables; var++) {
            m_cells->stg_fq[var][cellIndexBC(ii, j - 1, k)] = m_cells->stg_fq[var][IBC];
          }
        } else if(j == end[1] - m_noGhostLayers) {
          for(MInt var = 6; var < noSTGVariables; var++) {
            m_cells->stg_fq[var][cellIndexBC(ii, j + 1, k)] = m_cells->stg_fq[var][IBC];
          }
        }

        //////////////////////////////////////////////////
        // Storage of values
        //////////////////////////////////////////////////

        // Save max direction vector and max tau
        maxValsLocal[0] = umax;
        maxValsLocal[1] = vmax;
        maxValsLocal[2] = wmax;
        maxValsLocal[3] = utaumax;
      }
    }

    //////////////////////////////////////////////////
    // Communication: Exchange min and max values
    //////////////////////////////////////////////////

    MFloat minLengthGlobal = 0.0;
    MPI_Allreduce(maxValsLocal.begin(), maxValsGlobal.begin(), 4, MPI_DOUBLE, MPI_MAX, m_solver->m_commStg, AT_,
                  "maxValsLocal.begin()", "maxValsGlobal.begin()");
    MPI_Allreduce(&minLengthLocal, &minLengthGlobal, 1, MPI_DOUBLE, MPI_MIN, m_solver->m_commStg, AT_, "minLengthLocal",
                  "minLengthGlobal");


    // Maximum convection velocities at inflow
    m_stgMaxVel[0] = maxValsGlobal[0];
    m_stgMaxVel[1] = maxValsGlobal[1];
    m_stgMaxVel[2] = maxValsGlobal[2];
    const MFloat utaux = maxValsGlobal[3];

    for(MInt k = start[2]; k < end[2]; k++) {
      for(MInt j = start[1]; j < end[1]; j++) {
        I = cellIndex(ii, j, k);
        IBC = cellIndexBC(ii, j, k);

        // Length scale in main flow direction
        const MFloat xlength =
            max(min(m_solver->m_stgLengthFactors[0]
                        * max(m_cells->stg_fq[LENGTH_SCALE][IBC] * delta_in * pow(utaux / delta_in, exple), eps),
                    delta_in * 1.0),
                minLengthGlobal);

        // Length scale in the direction of main shear
        const MFloat ylength =
            max(min(m_solver->m_stgLengthFactors[1]
                        * max(m_cells->stg_fq[LENGTH_SCALE][IBC] * delta_in * pow(utaux / delta_in, exple), eps),
                    delta_in * 0.66),
                minLengthGlobal);

        // Length scale in the direction perpendicular of x and y
        const MFloat zlength =
            max(min(m_solver->m_stgLengthFactors[2]
                        * max(m_cells->stg_fq[LENGTH_SCALE][IBC] * delta_in * pow(utaux / delta_in, exple), eps),
                    delta_in * 1.0),
                minLengthGlobal);

        //>marian
        m_cells->stg_fq[LENGTH_X][IBC] = xlength;
        m_cells->stg_fq[LENGTH_Y][IBC] = ylength;
        m_cells->stg_fq[LENGTH_Z][IBC] = zlength;
        //<marian
      }
    }
  }

  if(m_solver->m_RKStep == 0) {
    //////////////////////////////////////////////////
    // The virtual box part - executed by Master Solver at the inflow
    //////////////////////////////////////////////////

    if(m_solver->m_stgMyRank == m_solver->m_commStgRoot) {
      MFloat epsik1 = 0.0, epsik2 = 0.0, epsik3 = 0.0;

      for(MInt n = 0; n < nran; n++) {
        MFloat xk1t = m_solver->m_stgEddies[n][0];
        MFloat xk2t = m_solver->m_stgEddies[n][1];
        MFloat xk3t = m_solver->m_stgEddies[n][2];

        // Check if the eddie has left the Virtual Box
        if(xk1t > m_stgVbEnd[0] || xk1t < m_stgVbStart[0] || xk2t > m_stgVbEnd[1] || xk2t < m_stgVbStart[1]
           || xk3t > m_stgVbEnd[2] || xk3t < m_stgVbStart[2]) {
          // Get coordinates of eddie cores and their signs
          // cout  << "Old eddie with position: " << xk1t << " , " << xk2t << " , " << xk3t << endl;
          xk1t = m_stgVbStart[0];
          xk2t = m_stgVbStart[1] + generate_rand_weighted() * BLT2;
          xk3t = m_stgVbStart[2] + generate_rand() * BLT3;
          // cout  << "Creating new eddie with position: " << xk1t << " , " << xk2t << " , " << xk3t << endl;
          epsik1 = 2.0 * generate_rand() - 1.0;
          epsik1 = epsik1 / max(abs(epsik1), eps);
          epsik2 = 2.0 * generate_rand() - 1.0;
          epsik2 = epsik2 / max(abs(epsik2), eps);
          epsik3 = 2.0 * generate_rand() - 1.0;
          epsik3 = epsik3 / max(abs(epsik3), eps);
        } else {
          xk1t = m_stgMaxVel[0] * m_solver->m_timeStep + m_solver->m_stgEddies[n][0];
          xk2t = m_stgMaxVel[1] * m_solver->m_timeStep + m_solver->m_stgEddies[n][1];
          xk3t = m_stgMaxVel[2] * m_solver->m_timeStep + m_solver->m_stgEddies[n][2];

          epsik1 = m_solver->m_stgEddies[n][3];
          epsik2 = m_solver->m_stgEddies[n][4];
          epsik3 = m_solver->m_stgEddies[n][5];
        }

        eddyBcastBuffer[n + nran * 0] = xk1t;
        eddyBcastBuffer[n + nran * 1] = xk2t;
        eddyBcastBuffer[n + nran * 2] = xk3t;
        eddyBcastBuffer[n + nran * 3] = epsik1;
        eddyBcastBuffer[n + nran * 4] = epsik2;
        eddyBcastBuffer[n + nran * 5] = epsik3;
      }
    }

    // Broadcast the new/updated eddies to all relevant processes
    MPI_Bcast(eddyBcastBuffer.begin(), nran * m_solver->m_stgNoEddieProperties, MPI_DOUBLE, m_solver->m_commStgRoot,
              m_solver->m_commStg, AT_, "eddyBcastBuffer.begin()");

    // Copy data into m_FQeddie vector
    for(MInt n = 0; n < nran; n++) {
      for(MInt p = 0; p < m_solver->m_stgNoEddieProperties; p++) {
        m_solver->m_stgEddies[n][p] = eddyBcastBuffer[n + nran * p];
      }
    }

    Vb = BLT2 * BLT3 * BLT1;

    // Summary of synth turb parameters
    if(m_solver->m_stgMyRank == m_solver->m_commStgRoot && globalTimeStep == m_solver->m_restartTimeStep) {
      cout << "**************************" << endl
           << "Synthetic turbulence:" << endl
           << "zones: 1" << endl
           << "nr. eddies: " << nran << endl
           << "conv. vel: " << sqrt(POW2(m_stgMaxVel[0]) + POW2(m_stgMaxVel[1]) + POW2(m_stgMaxVel[2])) << endl
           << "umax = " << m_stgMaxVel[0] << endl
           << "vmax = " << m_stgMaxVel[1] << endl
           << "wmax = " << m_stgMaxVel[2] << endl
           << "virtual box volume: " << Vb << endl
           << "Vb/nran = " << Vb / nran << endl
           << "**************************" << endl;
    }

    //////////////////////////////////////////////////
    // Calculation of the fluctuation induced by all eddies on each cell
    //////////////////////////////////////////////////

    // only compute fluctuations for second
    const MFloat vbFactor = sqrt(Vb / nran);
    const MInt iStart = 1, iEnd = 2;

    for(MInt k = start[2]; k < end[2]; k++) {
      for(MInt j = start[1]; j < end[1]; j++) {
        for(MInt i = iStart; i < iEnd; i++) {
          MFloat help1 = F0, help2 = F0, help3 = F0, help4 = F0, help5 = F0, help6 = F0;

          const MFloat umax = m_stgMaxVel[0];
          const MFloat vmax = m_stgMaxVel[1];
          const MFloat wmax = m_stgMaxVel[2];

          // the tensor components and xyzlengths are only saved in one row (ii = 1)
          const MInt cellIdBC = cellIndexBC(i, j, k);
          const MInt cellIdBCFirst = cellIndexBC(ii, j, k);
          const MInt cellId = cellIndex(ii, j, k);

          const MFloat uu = m_cells->stg_fq[FLUC_UU][cellIdBCFirst];
          const MFloat vv = m_cells->stg_fq[FLUC_VV][cellIdBCFirst];
          const MFloat ww = m_cells->stg_fq[FLUC_WW][cellIdBCFirst];
          const MFloat uv = m_cells->stg_fq[FLUC_UV][cellIdBCFirst];
          const MFloat vw = m_cells->stg_fq[FLUC_VW][cellIdBCFirst];
          const MFloat uw = m_cells->stg_fq[FLUC_UW][cellIdBCFirst];

          // Cholesky decomposition of the Reynolds stress tensor
          const MFloat a11 = sqrt(max(uu, epsl));
          const MFloat a21 = uv / a11;
          const MFloat a31 = uw / a11;
          const MFloat a22 = sqrt(max((vv - a21 * a21), epsl));
          const MFloat a32 = (vw - a21 * a31) / a22;
          const MFloat a33 = sqrt(max((ww - a31 * a31 - a32 * a32), epsl));

          const MFloat xLb1 = m_cells->stg_fq[LENGTH_X][cellIdBCFirst];
          const MFloat xLb2 = m_cells->stg_fq[LENGTH_Y][cellIdBCFirst];
          const MFloat xLb3 = m_cells->stg_fq[LENGTH_Z][cellIdBCFirst];

          const MFloat fxLb1 = 1.0 / xLb1;
          const MFloat fxLb2 = 1.0 / xLb2;
          const MFloat fxLb3 = 1.0 / xLb3;

          const MFloat fsqrtxLb1 = 1.0 / sqrt(xLb1);
          const MFloat fsqrtxLb2 = 1.0 / sqrt(xLb2);
          const MFloat fsqrtxLb3 = 1.0 / sqrt(xLb3);

          const MFloat fsqrtPixLb1 = 1.0 / sqrt(3.141 * xLb1);
          const MFloat fsqrtPixLb2 = 1.0 / sqrt(3.141 * xLb2);
          const MFloat fsqrtPixLb3 = 1.0 / sqrt(3.141 * xLb3);

          for(MInt n = 0; n < nran; n++) {
            // Tent function to determine the symmetric function to model the
            // decay of the fluctuations
            const MFloat xk1t = m_solver->m_stgEddies[n][0];
            const MFloat xk2t = m_solver->m_stgEddies[n][1];
            const MFloat xk3t = m_solver->m_stgEddies[n][2];

            const MFloat distX = m_cells->coordinates[0][cellId] - xk1t;
            const MFloat distY = m_cells->coordinates[1][cellId] - xk2t;
            const MFloat distZ = m_cells->coordinates[2][cellId] - xk3t;

            const MFloat aDistX = fabs(distX);
            const MFloat aDistY = fabs(distY);
            const MFloat aDistZ = fabs(distZ);

            // only compute contribution if eddie is in vicinity, i.e.,
            // if it is within 4 eddy lengthscales
            if(aDistX < 4.0 * xLb1 && aDistY < 4.0 * xLb2 && aDistZ < 4.0 * xLb3) {
              const MFloat zacfq1 = distX / aDistX;
              const MFloat rol1H = zacfq1 * min(aDistX * fxLb1, 1.0);

              const MFloat zacfq2 = distY / aDistY;
              const MFloat rol2H = zacfq2 * min(aDistY * fxLb2, 1.0);

              const MFloat zacfq3 = distZ / aDistZ;
              const MFloat rol3H = zacfq3 * min(aDistZ * fxLb3, 1.0);

              const MFloat fl1 = 2.0 * fsqrtPixLb1 * exp(-((distX)*2 * fxLb1) * ((distX)*2 * fxLb1));
              const MFloat fl2 = 2.0 * fsqrtPixLb2 * exp(-((distY)*2 * fxLb2) * ((distY)*2 * fxLb2));
              const MFloat fl3 = 2.0 * fsqrtPixLb3 * exp(-((distZ)*2 * fxLb3) * ((distZ)*2 * fxLb3));


              // Normalization factor cannot be chosen as Pamies did... or we...
              MFloat fH1 = (F1 - cos(2.0 * 3.141 * rol1H)) / (2.0 * 3.141 * rol1H * 0.44);
              MFloat fH2 = (F1 - cos(2.0 * 3.141 * rol2H)) / (2.0 * 3.141 * rol2H * 0.44);
              MFloat fH3 = (F1 - cos(2.0 * 3.141 * rol3H)) / (2.0 * 3.141 * rol3H * 0.44);

              fH1 = aniso * fH1 * fsqrtxLb1 + fabs(aniso - 1.0) * fl1;
              fH2 = aniso * fH2 * fsqrtxLb2 + fabs(aniso - 1.0) * fl2;
              fH3 = aniso * fH3 * fsqrtxLb3 + fabs(aniso - 1.0) * fl3;

              const MFloat epsik1 = m_solver->m_stgEddies[n][3];
              const MFloat epsik2 = m_solver->m_stgEddies[n][4];
              const MFloat epsik3 = m_solver->m_stgEddies[n][5];

              help4 += vbFactor * epsik1 * fl1 * fl2 * fH3;
              help5 += vbFactor * epsik2 * fl1 * fl2 * fH3;
              help6 += vbFactor * epsik3 * fl1 * fH2 * fl3;
            }
          }

          ///////////////////////////////////////////////////
          // Use Cholesky-Trafo to scale random fluctuations
          ///////////////////////////////////////////////////

          help1 = help4 * a11;                             // Fluctuation u'
          help2 = help4 * a21 + help5 * a22;               // Fluctuation v'
          help3 = help4 * a31 + help5 * a32 + help6 * a33; // Fluctuation w'

          const MFloat velmax = sqrt(umax * umax + vmax * vmax + wmax * wmax);

          const MFloat ufluc = min(max(help1, -0.3 * velmax), 0.3 * velmax);
          const MFloat vfluc = min(max(help2, -0.3 * velmax), 0.3 * velmax);
          const MFloat wfluc = min(max(help3, -0.3 * velmax), 0.3 * velmax);


          m_cells->stg_fq[FLUC_U][cellIdBC] += timsm * (ufluc - m_cells->stg_fq[FLUC_U][cellIdBC]);
          m_cells->stg_fq[FLUC_V][cellIdBC] += timsm * (vfluc - m_cells->stg_fq[FLUC_V][cellIdBC]);
          m_cells->stg_fq[FLUC_W][cellIdBC] += timsm * (wfluc - m_cells->stg_fq[FLUC_W][cellIdBC]);
        }
      }
    }
  } // RKStep end if

  /////////////////////////////////////////////////////////
  ////////////// APPLY TO BC //////////////////////////////
  /////////////////////////////////////////////////////////
  // Now comes the BC stuff that we need to do every RK step
  for(MInt j = start[1]; j < end[1]; j++) {
    for(MInt k = start[2]; k < end[2]; k++) {
      const MInt cellIdG1 = cellIndex(m_noGhostLayers - 1, j, k);
      IBC = cellIndexBC(m_noGhostLayers - 1, j, k);
      const MInt cellIdA1 = cellIndex(m_noGhostLayers, j, k);

      MFloat dxidx = m_cells->cellMetrics[0][cellIdA1];
      MFloat dxidy = m_cells->cellMetrics[1][cellIdA1];
      MFloat dxidz = m_cells->cellMetrics[2][cellIdA1];

      MFloat gradxi = -F1 / sqrt(dxidx * dxidx + dxidy * dxidy + dxidz * dxidz);

      MFloat dxHelp = dxidx * gradxi;
      MFloat dyHelp = dxidy * gradxi;
      MFloat dzHelp = dxidz * gradxi;

      const MFloat rhoBC = m_cells->pvariables[PV->RHO][cellIdG1];
      const MFloat pBC = m_cells->pvariables[PV->P][cellIdG1];
      const MFloat fRhoBC = F1 / rhoBC;
      const MFloat aBC = sqrt(gamma * pBC * fRhoBC);
      const MFloat uBC = m_cells->pvariables[PV->U][cellIdG1];
      const MFloat vBC = m_cells->pvariables[PV->V][cellIdG1];
      const MFloat wBC = m_cells->pvariables[PV->W][cellIdG1];

      const MFloat maBC = (dxHelp * uBC + dyHelp * vBC + dzHelp * wBC) / aBC;

      // get mean values from the rans
      const MFloat rhoRANS = m_cells->stg_fq[PV->RHO][IBC];
      const MFloat uRANS = m_cells->stg_fq[PV->U][IBC];
      const MFloat vRANS = m_cells->stg_fq[PV->V][IBC];
      const MFloat wRANS = m_cells->stg_fq[PV->W][IBC];
      const MFloat pRANS = m_cells->stg_fq[PV->P][IBC];

      // fluctuation values from the STG
      const MFloat u_prime = m_cells->stg_fq[FLUC_U][IBC];
      const MFloat v_prime = m_cells->stg_fq[FLUC_V][IBC];
      const MFloat w_prime = m_cells->stg_fq[FLUC_W][IBC];

      // superpose onto mean RANS variables
      const MFloat uSTG = max(uRANS + u_prime, epsl);
      const MFloat vSTG = vRANS + v_prime;
      const MFloat wSTG = wRANS + w_prime;

      // compute correct density
      const MFloat u9a = PV->UInfinity;
      const MFloat u9ff = u_prime;
      const MFloat alok = sqrt(gamma * PV->PInfinity / CV->rhoInfinity);
      const MFloat flucc = u9ff / u9a * POW2((PV->UInfinity / alok)) * gammaMinusOne * m_cells->stg_fq[PV->RHO][IBC];
      const MFloat zdir = flucc / max(fabs(flucc), 0.0000001);
      const MFloat rhoSTG = rhoRANS + zdir * min(fabs(flucc), 0.1 * rhoRANS);

      // get field values inside the integration domain
      const MFloat pField = m_cells->pvariables[PV->P][cellIdA1];
      const MFloat rhoField = m_cells->pvariables[PV->RHO][cellIdA1];
      const MFloat uField = m_cells->pvariables[PV->U][cellIdA1];
      const MFloat vField = m_cells->pvariables[PV->V][cellIdA1];
      const MFloat wField = m_cells->pvariables[PV->W][cellIdA1];
      const MFloat aField = sqrt(gamma * pField / rhoField);

      /////////////////////////////////////////////////
      //////////// SUBSONIC PART //////////////////////
      /////////////////////////////////////////////////
      const MFloat pSub =
          F1B2
          * (pField + pRANS
             + rhoField * aField * (+dxHelp * (uField - uSTG) + dyHelp * (vField - vSTG) + dzHelp * (wField - wSTG)));
      const MFloat rhoSub = rhoSTG + (pSub - pRANS) / (POW2(aField));
      const MFloat rhoSubHelp = (pSub - pRANS) / (rhoField * aField);

      // Multiply velocities with density
      const MFloat uSub = uSTG + dxHelp * rhoSubHelp;
      const MFloat vSub = vSTG + dyHelp * rhoSubHelp;
      const MFloat wSub = wSTG + dzHelp * rhoSubHelp;

      /////////////////////////////////////////////////
      //////////// SUPERSONIC PART ////////////////////
      /////////////////////////////////////////////////
      const MFloat rhoSup = rhoSTG;
      const MFloat uSup = uSTG;
      const MFloat vSup = vSTG;
      const MFloat wSup = wSTG;
      const MFloat pSup = pRANS;

      //////////////////////////////////////////////////
      /////////// SUB/SUP INTERPOLATION ////////////////
      //////////////////////////////////////////////////

      // by default the subsonic formulation is used
      // switch on "stgSubSup" to get the mixed formulation
      // or "stgSupersonic" to use the pure supersonic formulation
      MFloat xSub = F1;
      MFloat xSup = F0;

      if(m_solver->m_stgSubSup) {
        const MFloat maBCAbs = fabs(maBC);
        const MFloat alpha = 14.0;
        const MFloat b = 0.95;
        const MFloat count = alpha * (maBCAbs - b);
        const MFloat denom = (F1 - 0.99 * b) * maBCAbs + b;
        const MFloat ratio = count / denom;
        const MFloat wfun = F1B2 * (F1 + tanh(ratio) / tanh(alpha));

        xSub = fabs(wfun - F1);
        xSup = fabs(wfun);
      } else if(m_solver->m_stgSupersonic) {
        xSub = F0;
        xSup = F1;
      }

      m_cells->pvariables[PV->RHO][cellIdG1] = rhoSub * xSub + rhoSup * xSup;
      m_cells->pvariables[PV->U][cellIdG1] = uSub * xSub + uSup * xSup;
      m_cells->pvariables[PV->V][cellIdG1] = vSub * xSub + vSup * xSup;
      m_cells->pvariables[PV->W][cellIdG1] = wSub * xSub + wSup * xSup;
      m_cells->pvariables[PV->P][cellIdG1] = pSub * xSub + pSup * xSup;


      //////////////////////////////////////////////////
      //////// EXTRAPOLATE TO SECOND GC ////////////////
      //////////////////////////////////////////////////
      const MInt cellIdG2 = cellIndex(0, j, k);
      // extrapolate into second ghost cell
      for(MInt var = 0; var < PV->noVariables; var++) {
        m_cells->pvariables[var][cellIdG2] =
            F2 * m_cells->pvariables[var][cellIdG1] - m_cells->pvariables[var][cellIdA1];
      }
    }
  }
}
//<marian


/**
 * \brief Rescaling Boundary Conditions
 *
 *  BC2500 and BC2501 a combined rescaling boundary condition
 *  and can only be used together. The 2501 is the recycling
 *  station from where the values are taken and 2500 is the
 *  inflow plane where the rescaled values are prescribed.
 *
 */

// Inlet station
template <MBool isRans>
void StructuredBndryCnd3D<isRans>::bc2500(MInt bcId) {
  (void)bcId;
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

  // compute the momentum thickness at the inlet (i==x, j==y, k==z)
  const MInt i = m_noGhostLayers - 1;


  //////////////////////////////////////////////////////////////
  ///////////////// COMPUTE THETA AND EXCHANGE /////////////////
  //////////////////////////////////////////////////////////////

  // allocate space in k direction (all k-Cells)
  MFloatScratchSpace thetaLocal(2, AT_, "thetaLocalIn");
  thetaLocal.fill(F0); // initialize scratch space

  MFloatScratchSpace thetaGlobal(2, AT_, "thetaGlobalIn");
  thetaGlobal.fill(F0);

  // the offest position in k-direction is the offset
  const MInt thetaLocalOffset = m_solver->m_nOffsetCells[0];

  // compute the local moment thickness j=direction of integration
  for(MInt k = m_noGhostLayers; k < m_nCells[0] - m_noGhostLayers; ++k) {
    for(MInt j = m_noGhostLayers; j < m_nCells[1] - m_noGhostLayers; ++j) {
      const MInt cellId = cellIndex(i, j, k);
      const MInt pointIdM1 = getPointIdFromCell(i, j, k);
      const MInt pointIdP1 = getPointIdfromPoint(pointIdM1, 0, 1, 0);

      if(m_grid->m_coordinates[1][pointIdM1] > maxIntegrationHeight) {
        continue;
      }

      const MFloat urat = m_cells->pvariables[PV->U][cellId] / PV->UInfinity;
      const MFloat momThick =
          (m_cells->pvariables[PV->U][cellId] * m_cells->pvariables[PV->RHO][cellId] * fabs(F1 - urat))
          / (CV->rhoUInfinity);

      // integrate normal to the wall
      const MFloat ydist = m_grid->m_coordinates[1][pointIdP1] - m_grid->m_coordinates[1][pointIdM1];
      thetaLocal(0) += momThick * ydist;
    }
  }

  MPI_Allreduce(thetaLocal.begin(), thetaGlobal.begin(), 2, MPI_DOUBLE, MPI_SUM, m_solver->m_rescalingCommGrComm, AT_,
                "thetaLocal.begin()", "thetaGlobal.begin()");

  // determine spanwise average: we assume equally spaced in z-direction
  for(MInt ii = 0; ii < 2; ++ii) {
    thetaGlobal(ii) /= (m_grid->getMyBlockNoCells(0));
  }

  if(globalTimeStep % 50 == 0 && m_solver->m_RKStep == 0
     && m_solver->domainId() == m_solver->m_rescalingCommGrRootGlobal) {
    cout << "ThetaInflow: " << thetaGlobal(0) << " ThetaRecyclingStation: " << thetaGlobal(1) << endl;

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
  const MInt noVar = 2; // for more variables if wanted
  MFloatScratchSpace wallPropertiesLocal(m_grid->getMyBlockNoCells(0), noVar, AT_, "wallPropertiesLocalInlet");
  MFloatScratchSpace wallProperties(m_grid->getMyBlockNoCells(0), noVar, AT_, "wallPropertiesInlet");
  wallPropertiesLocal.fill(F0);
  wallProperties.fill(F0);

  MPI_Allreduce(wallPropertiesLocal.begin(), wallProperties.begin(), m_grid->getMyBlockNoCells(0) * noVar, MPI_DOUBLE,
                MPI_SUM, m_solver->m_rescalingCommGrComm, AT_, "wallPropertiesLocal.begin()", "wallProperties.begin()");

  //////////////////////////////////////////////////////////////
  ///////////////// GAMS, UTAUIN, INNER OUTER COORD ////////////
  //////////////////////////////////////////////////////////////
  MFloatScratchSpace utauIn(m_nCells[0], AT_, "u_tauIn");
  MFloatScratchSpace gams(m_nCells[0], AT_, "gams");

  MInt kStart = 0;
  MInt kEnd = m_nCells[0];

  if(m_solver->m_nOffsetCells[0] == 0) {
    kStart = m_noGhostLayers;
  }
  if(m_solver->m_nOffsetCells[0] + m_solver->m_nActiveCells[0] == m_grid->getMyBlockNoCells(0)) {
    kEnd = m_nCells[0] - m_noGhostLayers;
  }

  for(MInt k = kStart; k < kEnd; ++k) {
    const MFloat utauRe = wallProperties(thetaLocalOffset + (k - m_noGhostLayers), 0);

    // estimate the friction velocity at the inlet
    // according to standard power law approximations
    // utau_in = utau_re*(theta_re/theta_in)**(1/2*(n-1))
    // where theta is the momentum thickness
    // see Thomas S. Lund, p241

    // when take into account the variance of wall density
    // utau_in = utau_re*(rho_wall_re/rho_wall_in)**0.5
    //* (theta_re/theta_in)**(1/2*(n-1))
    // here n = 5

    gams(k) = pow(thetaGlobal(1) / fabs(thetaGlobal(0)), F1B8);
    utauIn(k) = utauRe * min(max(gams(k), F1), 2.5);
  }

  MFloatScratchSpace coordInInner(m_nCells[0] * m_nCells[1], AT_, "coordInInner");
  MFloatScratchSpace coordInOuter(m_nCells[0] * m_nCells[1], AT_, "coordInOuter");

  for(MInt k = 0; k < m_nCells[0]; ++k) {
    for(MInt j = 0; j < m_nCells[1]; ++j) {
      const MInt cellId = cellIndex(i, j, k);
      const MInt faceId = j + k * m_nCells[1];
      const MFloat rho = m_cells->pvariables[PV->RHO][cellId];
      const MFloat frho = F1 / rho;
      const MFloat p = m_cells->pvariables[PV->P][cellId];
      const MFloat temp = p * gamma * frho;
      const MFloat mu = SUTHERLANDLAW(temp);

      coordInInner(faceId) = utauIn(k) * rho * m_cells->coordinates[1][cellId] / (mu * sqrt(m_solver->m_Re0));
      coordInOuter(faceId) = m_cells->coordinates[1][cellId] * rho / (m_rescalingBLT * CV->rhoInfinity);
    }
  }

  //////////////////////////////////////////////////////////////
  ///////////////// NOW EXCHANGE VAR SLICE /////////////////////
  //////////////////////////////////////////////////////////////
  const MInt noVariables = 6;
  const MInt totalCells[2] = {m_grid->getMyBlockNoCells(0), m_grid->getMyBlockNoCells(1) + 1};
  MFloatScratchSpace varSliceLocal(noVariables, totalCells[0] * totalCells[1], AT_, "varSliceLocal");
  MFloatScratchSpace varSlice(noVariables, totalCells[0] * totalCells[1], AT_, "varSlice");

  // we are at the inlet, only fill with zeros
  varSlice.fill(F0);
  varSliceLocal.fill(F0);

  MPI_Allreduce(varSliceLocal.begin(), varSlice.begin(), noVariables * totalCells[0] * totalCells[1], MPI_DOUBLE,
                MPI_SUM, m_solver->m_rescalingCommGrComm, AT_, "varSliceLocal.begin()", "varSlice.begin()");


  ///////////////////////////////////////////////////////////////
  ///////////////// RESCALING ///////////////////////////////////
  ///////////////////////////////////////////////////////////////

  MInt jStart = 0;
  MInt jEnd = m_nCells[1];

  if(m_solver->m_nOffsetCells[1] == 0) {
    jStart = m_noGhostLayers;
  }
  if(m_solver->m_nOffsetCells[1] + m_solver->m_nActiveCells[1] == m_grid->getMyBlockNoCells(1)) {
    jEnd = m_nCells[1] - m_noGhostLayers;
  }

  for(MInt k = kStart; k < kEnd; ++k) {
    const MFloat ctem1 = (F1 + ctema) * (F1 - POW2(gams(k)));
    const MFloat ctem2 = F2 * ctema * gams(k) * (F1 - gams(k));
    const MFloat ctem3 = (F1 - gams(k)) * (F1 + gams(k) + F2 * ctema * gams(k));

    for(MInt j = jStart; j < jEnd; ++j) {
      const MInt faceId = j + k * m_nCells[1];
      // const MInt faceIdM1 = (j-1)+k*m_nCells[1];
      const MInt cellId = cellIndex(i, j, k);

      if(coordInOuter(faceId) < 1.05) {
        MFloat uInner = F0, vInner = F0, wInner = F0, TInner = F0;
        MFloat uOuter = F0, vOuter = F0, wOuter = F0, TOuter = F0;
        const MFloat count = alpha * (coordInOuter(faceId) - b);
        const MFloat denom = (F1 - F2 * b) * coordInOuter(faceId) + b;
        const MFloat ratio = count / denom;
        const MFloat wfun = F1B2 * (F1 + tanh(ratio) / tanh(alpha));

        for(MInt jj = 0; jj < totalCells[1] - 1; ++jj) {
          const MInt localId = jj + (m_solver->m_nOffsetCells[0] + (k - m_noGhostLayers)) * totalCells[1];
          const MInt localIdP1 = (jj + 1) + (m_solver->m_nOffsetCells[0] + (k - m_noGhostLayers)) * totalCells[1];

          const MFloat yInnerRe = varSlice(4, localId);
          const MFloat yInnerReP1 = varSlice(4, localIdP1);

          if((yInnerRe - coordInInner(faceId)) < rescalEPS && yInnerReP1 > coordInInner(faceId)) {
            const MFloat dy1 = coordInInner(faceId) - yInnerRe;
            const MFloat dy2 = yInnerReP1 - coordInInner(faceId);
            const MFloat dy = yInnerReP1 - yInnerRe;

            const MFloat u = varSlice(0, localId);
            const MFloat uP1 = varSlice(0, localIdP1);
            const MFloat v = varSlice(1, localId);
            const MFloat vP1 = varSlice(1, localIdP1);
            const MFloat w = varSlice(2, localId);
            const MFloat wP1 = varSlice(2, localIdP1);
            const MFloat t = varSlice(3, localId);
            const MFloat tP1 = varSlice(3, localIdP1);
            uInner = (uP1 * dy1 + u * dy2) / dy;
            vInner = (vP1 * dy1 + v * dy2) / dy;
            wInner = (wP1 * dy1 + w * dy2) / dy;
            TInner = (tP1 * dy1 + t * dy2) / dy;
          }
        }

        // outer region
        for(MInt jj = 0; jj < totalCells[1] - 1; ++jj) {
          const MInt localId = jj + (m_solver->m_nOffsetCells[0] + (k - m_noGhostLayers)) * totalCells[1];
          const MInt localIdP1 = (jj + 1) + (m_solver->m_nOffsetCells[0] + (k - m_noGhostLayers)) * totalCells[1];

          const MFloat yOuterRe = varSlice(5, localId);
          const MFloat yOuterReP1 = varSlice(5, localIdP1);

          if((yOuterRe - coordInOuter(faceId)) < rescalEPS && yOuterReP1 > coordInOuter(faceId)) {
            const MFloat dy1 = coordInOuter(faceId) - yOuterRe;
            const MFloat dy2 = yOuterReP1 - coordInOuter(faceId);
            const MFloat dy = yOuterReP1 - yOuterRe;

            const MFloat u = varSlice(0, localId);
            const MFloat uP1 = varSlice(0, localIdP1);
            const MFloat v = varSlice(1, localId);
            const MFloat vP1 = varSlice(1, localIdP1);
            const MFloat w = varSlice(2, localId);
            const MFloat wP1 = varSlice(2, localIdP1);
            const MFloat t = varSlice(3, localId);
            const MFloat tP1 = varSlice(3, localIdP1);
            uOuter = (uP1 * dy1 + u * dy2) / dy;
            vOuter = (vP1 * dy1 + v * dy2) / dy;
            wOuter = (wP1 * dy1 + w * dy2) / dy;
            TOuter = (tP1 * dy1 + t * dy2) / dy;
          }
        }

        const MFloat TInnerA = POW2(gams(k)) * TInner + ctem1 * PV->TInfinity;
        const MFloat TOuterA = POW2(gams(k)) * TOuter - (ctem2 * (uOuter / PV->UInfinity) - ctem3) * PV->TInfinity;

        // van Driest transformation
        const MFloat uvdInner = PV->UInfinity * asin(b_vd * uInner / PV->UInfinity) / b_vd;
        const MFloat uvdOuter = PV->UInfinity * asin(b_vd * uOuter / PV->UInfinity) / b_vd;

        // scaling of transformed inner and outer velocities
        uInner = gams(k) * uvdInner;
        uOuter = gams(k) * uvdOuter + (F1 - gams(k)) * uvd8;
        uInner = PV->UInfinity * sin(b_vd * uInner / PV->UInfinity) / b_vd;
        uOuter = PV->UInfinity * sin(b_vd * uOuter / PV->UInfinity) / b_vd;

        const MFloat pres = PV->PInfinity;
        const MFloat uMean = uInner * (F1 - wfun) + uOuter * wfun;
        const MFloat vMean = vInner * (F1 - wfun) + vOuter * wfun;
        const MFloat wMean = (wInner * (F1 - wfun) + wOuter * wfun) * gams(k);
        const MFloat tMean = TInnerA * (F1 - wfun) + TOuterA * wfun;
        const MFloat rhoIn = gamma * pres / tMean;

        // //clebanoff factor is optional
        // const MFloat clebf = 6.1;
        // const MFloat blt   = m_rescalingBLT;
        // const MFloat cleb  = F1/(F1+pow((m_cells->coordinates[1][cellId]/(clebf*blt)), 6.0));

        m_cells->pvariables[PV->RHO][cellId] = rhoIn;
        m_cells->pvariables[PV->U][cellId] = uMean;
        m_cells->pvariables[PV->V][cellId] = vMean;
        m_cells->pvariables[PV->W][cellId] = wMean;
        m_cells->pvariables[PV->P][cellId] = pres;
      } else {
        // if(!edgePointIsSet(k)) {
        //   edgePointJ(k) = j;
        //   edgePointIsSet(k) = 1;
        // }

        const MFloat pres = PV->PInfinity;
        const MFloat rhoIn = gamma * pres / PV->TInfinity;

        const MFloat uMean = PV->UInfinity;
        const MFloat vMean = PV->VInfinity;
        const MFloat wMean = PV->WInfinity;

        m_cells->pvariables[PV->RHO][cellId] = rhoIn;
        m_cells->pvariables[PV->U][cellId] = uMean;
        m_cells->pvariables[PV->V][cellId] = vMean;
        m_cells->pvariables[PV->W][cellId] = wMean;
        m_cells->pvariables[PV->P][cellId] = pres;
      }
    }
  }

  for(MInt k = kStart; k < kEnd; ++k) {
    for(MInt j = 0; j < m_nCells[1]; ++j) {
      // extrapolation for second GC
      const MInt cellId = cellIndex(1, j, k);
      const MInt cellIdM1 = cellIndex(0, j, k);
      const MInt cellIdadj = cellIndex(2, j, k);

      for(MInt var = 0; var < PV->noVariables; var++) {
        m_cells->pvariables[var][cellIdM1] =
            2.0 * m_cells->pvariables[var][cellId] - m_cells->pvariables[var][cellIdadj];
      }
    }
  }
}


// Recycling station
template <MBool isRans>
void StructuredBndryCnd3D<isRans>::bc2501(MInt bcId) {
  MInt* start = m_physicalBCMap[bcId]->start1;

  const MFloat gamma = m_solver->m_gamma;

  // things to move to init or elsewhere
  MFloat F727 = 72.0 / 7.0;
  const MInt i = start[0]; // position at which recycle is taken
  const MFloat yWall = F0; // this has been fixed else method does not work
  const MFloat maxIntegrationHeight = 2.0 * m_rescalingBLT;

  //////////////////////////////////////////////////////////////
  ///////////////// COMPUTE THETA AND EXCHANGE /////////////////
  //////////////////////////////////////////////////////////////

  // thetaLocal.fill(F0); //initialize scratch space to zero // only for parallel use
  MFloatScratchSpace thetaLocal(2, AT_, "thetaLocalRe");
  MFloatScratchSpace thetaGlobal(2, AT_, "thetaGlobalRe");
  thetaLocal.fill(F0);
  thetaGlobal.fill(F0);

  // the offest position in k-direction is the offset
  const MInt thetaLocalOffset = m_solver->m_nOffsetCells[0];

  // compute the local moment thickness j=direction of integration
  for(MInt k = m_noGhostLayers; k < m_nCells[0] - m_noGhostLayers; ++k) {
    for(MInt j = m_noGhostLayers; j < m_nCells[1] - m_noGhostLayers; ++j) {
      const MInt cellId = cellIndex(i, j, k);
      const MInt pointIdM1 = getPointIdFromCell(i, j, k);
      const MInt pointIdP1 = getPointIdfromPoint(pointIdM1, 0, 1, 0);

      if(m_grid->m_coordinates[1][pointIdM1] > maxIntegrationHeight) {
        continue;
      }

      const MFloat urat = m_cells->pvariables[PV->U][cellId] / PV->UInfinity;
      const MFloat momThick =
          (m_cells->pvariables[PV->U][cellId] * m_cells->pvariables[PV->RHO][cellId] * fabs(F1 - urat))
          / (CV->rhoUInfinity);

      // integrate normal to the wall
      const MFloat ydist = m_grid->m_coordinates[1][pointIdP1] - m_grid->m_coordinates[1][pointIdM1];
      thetaLocal(1) += momThick * ydist;
    }
  }

  // communicate the Thickness across the plane
  MPI_Allreduce(thetaLocal.begin(), thetaGlobal.begin(), 2, MPI_DOUBLE, MPI_SUM, m_solver->m_rescalingCommGrComm, AT_,
                "thetaLocal.begin()", "thetaGlobal.begin()");

  // determine spanwise average: we assume equally spaced in z-direction
  for(MInt ii = 0; ii < 2; ++ii) {
    thetaGlobal(ii) /= (m_grid->getMyBlockNoCells(0));
  }

  //////////////////////////////////////////////////////////////
  ///////////////// COMPUTE AND EXCHANGE WALL PROPERTIES ///////
  //////////////////////////////////////////////////////////////

  const MFloat delta = F727 * thetaGlobal(1);
  const MInt noVar = 2;                                     // for more variables if wanted
  const MInt wallLocalOffset = m_solver->m_nOffsetCells[1]; // Offset in j-direction

  MFloatScratchSpace wallPropertiesLocal(m_grid->getMyBlockNoCells(0), noVar, AT_, "wallPropertiesLocalRe");
  MFloatScratchSpace wallProperties(m_grid->getMyBlockNoCells(0), noVar, AT_, "wallPropertiesRe");
  wallPropertiesLocal.fill(F0);
  wallProperties.fill(F0);

  // determine the wall stuff if wall is contained whithin the partition
  if(wallLocalOffset == 0 && m_solver->m_nActiveCells[1] >= m_noGhostLayers) {
    for(MInt k = m_noGhostLayers; k < m_solver->m_nCells[0] - m_noGhostLayers; ++k) {
      const MInt cellId = cellIndex(i, m_noGhostLayers, k);
      const MInt localId = m_solver->m_nOffsetCells[0] + (k - m_noGhostLayers);
      const MFloat rho = m_cells->pvariables[PV->RHO][cellId];
      const MFloat p = m_cells->pvariables[PV->P][cellId];
      const MFloat t = p * gamma / rho;
      const MFloat mu = SUTHERLANDLAW(t);
      const MFloat uWall = fabs(m_cells->pvariables[PV->U][cellId]);
      const MFloat ydist = m_cells->coordinates[1][cellId] - yWall;
      const MFloat uTau = sqrt(uWall * mu / (ydist * rho));

      wallPropertiesLocal(localId, 0) = uTau;
      wallPropertiesLocal(localId, 1) = rho;
    }
  }

  MPI_Allreduce(wallPropertiesLocal.begin(), wallProperties.begin(), noVar * m_grid->getMyBlockNoCells(0), MPI_DOUBLE,
                MPI_SUM, m_solver->m_rescalingCommGrComm, AT_, "wallPropertiesLocal.begin()", "wallProperties.begin()");

  //////////////////////////////////////////////////////////////
  ///////////////// COMPUTE AND EXCHANGE VAR SLICE /////////////
  //////////////////////////////////////////////////////////////

  const MInt totalCells[2] = {m_grid->getMyBlockNoCells(0), m_grid->getMyBlockNoCells(1) + 1};

  MFloatScratchSpace varSliceLocal(6, totalCells[0] * totalCells[1], AT_, "varSliceLocal");
  MFloatScratchSpace varSlice(6, totalCells[0] * totalCells[1], AT_, "varSlice");

  for(MInt k = m_noGhostLayers; k < m_nCells[0] - m_noGhostLayers; ++k) {
    for(MInt j = m_noGhostLayers; j < m_nCells[1] - m_noGhostLayers; ++j) {
      const MInt cellId = cellIndex(i, j, k);
      const MInt localId = (m_solver->m_nOffsetCells[1] + (j - m_noGhostLayers) + 1)
                           + (m_solver->m_nOffsetCells[0] + (k - m_noGhostLayers)) * totalCells[1];
      const MFloat rho = m_cells->pvariables[PV->RHO][cellId];
      const MFloat frho = F1 / rho;
      const MFloat p = m_cells->pvariables[PV->P][cellId];
      const MFloat temp = p * gamma * frho;
      const MFloat mu = SUTHERLANDLAW(temp);
      const MFloat uTauRe = wallProperties(thetaLocalOffset + (k - m_noGhostLayers), 0);
      const MFloat yIn = (m_cells->coordinates[1][cellId] - yWall) * uTauRe * rho / (mu * sqrt(m_solver->m_Re0));
      const MFloat yOut = (m_cells->coordinates[1][cellId] - yWall) * rho / (delta * CV->rhoInfinity);
      const MFloat u = m_cells->pvariables[PV->U][cellId];
      const MFloat v = m_cells->pvariables[PV->V][cellId];
      const MFloat w = m_cells->pvariables[PV->W][cellId];

      // save the variables u,v,w,t,yI,yO
      varSliceLocal(0, localId) = u;
      varSliceLocal(1, localId) = v;
      varSliceLocal(2, localId) = w;
      varSliceLocal(3, localId) = temp;
      varSliceLocal(4, localId) = yIn;
      varSliceLocal(5, localId) = yOut;
    }

    // set first value at the wall manually
    if(m_solver->m_nOffsetCells[0] == 0) {
      const MInt localId = 0 + (m_solver->m_nOffsetCells[0] + (k - m_noGhostLayers)) * totalCells[1];
      const MInt localIdP1 = 1 + (m_solver->m_nOffsetCells[0] + (k - m_noGhostLayers)) * totalCells[1];
      varSliceLocal(0, localId) = 0.0;                         // u
      varSliceLocal(1, localId) = 0.0;                         // v
      varSliceLocal(2, localId) = 0.0;                         // w
      varSliceLocal(3, localId) = varSliceLocal(3, localIdP1); // t
      varSliceLocal(4, localId) = 0.0;                         // yIn
      varSliceLocal(5, localId) = 0.0;                         // yOut
    }
  }

  // communicate the slice
  MPI_Allreduce(varSliceLocal.begin(), varSlice.begin(), 6 * totalCells[0] * totalCells[1], MPI_DOUBLE, MPI_SUM,
                m_solver->m_rescalingCommGrComm, AT_, "varSliceLocal.begin()", "varSlice.begin()");

  // participate in communication but only fill with zeros
  // MFloatScratchSpace blEdgeVValueLocal(m_solver->m_totalGridBlockDim[0][0]-1,AT_, "blEdgeVValueLocal");
  // MFloatScratchSpace blEdgeVValueGlobal(m_solver->m_totalGridBlockDim[0][0]-1,AT_, "blEdgeVValueLocal");
  // blEdgeVValueLocal.fill(F0);
  // blEdgeVValueGlobal.fill(F0);

  // MPI_Allreduce(blEdgeVValueLocal.begin(),blEdgeVValueGlobal.begin(), m_solver->m_totalGridBlockDim[0][0]-1,
  // MPI_DOUBLE, MPI_SUM, m_solver->m_rescalingCommGrComm, AT_, "blEdgeVValueLocal.begin()",
  // "blEdgeVValueGlobal.begin()" );
}


// Inlet station for RANS
template <MBool isRans>
void StructuredBndryCnd3D<isRans>::bc2510(MInt bcId) {
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

  // the offest position in k-direction is the offset
  const MInt thetaLocalOffset = m_solver->m_nOffsetCells[0];
  // compute the local moment thickness j=direction of integration
  for(MInt k = m_noGhostLayers; k < m_nCells[0] - m_noGhostLayers; ++k) {
    for(MInt j = m_noGhostLayers; j < m_nCells[1] - m_noGhostLayers; ++j) {
      const MInt cellId = cellIndex(i, j, k);
      const MInt pointIdM1 = getPointIdFromCell(i, j, k);
      const MInt pointIdP1 = getPointIdfromPoint(pointIdM1, 0, 1, 0);

      const MFloat urat = m_cells->pvariables[PV->U][cellId] / PV->UInfinity;
      const MFloat momThick =
          (m_cells->pvariables[PV->U][cellId] * m_cells->pvariables[PV->RHO][cellId] * fabs(F1 - urat))
          / (CV->rhoUInfinity);
      // integrate normal to the wall
      const MFloat ydist = m_grid->m_coordinates[1][pointIdP1] - m_grid->m_coordinates[1][pointIdM1];
      thetaLocal(0) += momThick * ydist;
    }
  }

  MPI_Allreduce(thetaLocal.begin(), thetaGlobal.begin(), 2, MPI_DOUBLE, MPI_SUM, m_solver->m_rescalingCommGrComm, AT_,
                "thetaLocal.begin()", "thetaGlobal.begin()");

  // determine spanwise average: we assume equally spaced in z-direction
  for(MInt ii = 0; ii < 2; ++ii) {
    thetaGlobal(ii) /= m_grid->getMyBlockNoCells(0);
  }

  if(globalTimeStep % 5 == 0 && m_solver->m_RKStep == 0
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
  MFloatScratchSpace wallPropertiesLocal(m_grid->getMyBlockNoCells(0), noWallProperties, AT_,
                                         "wallPropertiesLocalInlet");
  MFloatScratchSpace wallProperties(m_grid->getMyBlockNoCells(0), noWallProperties, AT_, "wallPropertiesInlet");
  wallPropertiesLocal.fill(F0);
  wallProperties.fill(F0);

  MPI_Allreduce(wallPropertiesLocal.begin(), wallProperties.begin(), m_grid->getMyBlockNoCells(0) * noWallProperties,
                MPI_DOUBLE, MPI_SUM, m_solver->m_rescalingCommGrComm, AT_, "wallPropertiesLocal.begin()",
                "wallProperties.begin()");

  //////////////////////////////////////////////////////////////
  ///////////////// GAMS, UTAUIN, INNER OUTER COORD ////////////
  //////////////////////////////////////////////////////////////

  MFloatScratchSpace utauIn(m_solver->m_nActiveCells[0], 2, AT_, "u_tauIn");
  MFloatScratchSpace gams(m_solver->m_nActiveCells[0], 2, AT_, "gams");

  for(MInt k = 0; k < m_solver->m_nActiveCells[0]; ++k) {
    MFloat utauRe = wallProperties(thetaLocalOffset + (k), 0);

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

    gams(k, i) = pow(thetaGlobal(1) / fabs(thetaGlobal(0)), facc);
    utauIn(k, i) = utauRe * min(max(gams(k, i), F1), 2.5);
  }

  MFloatScratchSpace coordInInner(m_nCells[0] * m_nCells[1], 2, AT_, "coordInInner");
  MFloatScratchSpace coordInOuter(m_nCells[0] * m_nCells[1], 2, AT_, "coordInOuter");

  for(MInt k = m_noGhostLayers; k < m_nCells[0] - m_noGhostLayers; ++k) {
    for(MInt j = 0; j < m_nCells[1]; ++j) {
      const MInt cellId = cellIndex(i, j, k);
      const MInt localId = j + (k - m_noGhostLayers) * m_nCells[1];
      const MFloat rho = m_cells->pvariables[PV->RHO][cellId];
      const MFloat frho = F1 / rho;
      const MFloat p = m_cells->pvariables[PV->P][cellId];
      const MFloat temp = p * gamma * frho;
      const MFloat mu = SUTHERLANDLAW(temp);

      coordInInner(localId, i) =
          utauIn(k - m_noGhostLayers, i) * rho * m_cells->coordinates[1][cellId] / (mu * sqrt(m_solver->m_Re0));
      coordInOuter(localId, i) = m_cells->coordinates[1][cellId] * rho / (m_rescalingBLT * CV->rhoInfinity);
    }
  }

  const MInt wallLocalOffset = m_solver->m_nOffsetCells[1]; // Offset in j-direction
  MFloatScratchSpace tempWallInletLocal(m_grid->getMyBlockNoCells(0), AT_, "tempWallInletLocal");
  MFloatScratchSpace tempWallInletGlobal(m_grid->getMyBlockNoCells(0), AT_, "tempWallInletGlobal");
  tempWallInletLocal.fill(F0);
  tempWallInletGlobal.fill(F0);
  // determine the wall stuff if wall is contained whithin the partition
  if(wallLocalOffset == 0 && m_solver->m_nActiveCells[1] >= m_noGhostLayers) {
    for(MInt k = m_noGhostLayers; k < m_solver->m_nCells[0] - m_noGhostLayers; ++k) {
      const MInt cellId = cellIndex(1, m_noGhostLayers, k);
      const MInt localId = m_solver->m_nOffsetCells[0] + (k - m_noGhostLayers);
      tempWallInletLocal(localId) = temperature(cellId);
    }
  }

  MPI_Allreduce(tempWallInletLocal.begin(), tempWallInletGlobal.begin(), m_grid->getMyBlockNoCells(0), MPI_DOUBLE,
                MPI_SUM, m_solver->m_rescalingCommGrComm, AT_, "tempWallInletLocal.begin()",
                "tempWallInletGlobal.begin()");


  //////////////////////////////////////////////////////////////
  ///////////////// NOW EXCHANGE VAR SLICE /////////////////////
  //////////////////////////////////////////////////////////////
  const MInt noVariables = PV->noVariables + 1;
  MInt totalCells[2] = {m_grid->getMyBlockNoCells(0), m_grid->getMyBlockNoCells(1)};
  MFloatScratchSpace varSliceLocal(noVariables, totalCells[0] * totalCells[1], AT_, "varSliceLocal");
  MFloatScratchSpace varSlice(noVariables, totalCells[0] * totalCells[1], AT_, "varSlice");

  // we are at the inlet, only fill with zeros
  varSlice.fill(F0);
  varSliceLocal.fill(F0);

  MPI_Allreduce(varSliceLocal.begin(), varSlice.begin(), noVariables * totalCells[0] * totalCells[1], MPI_DOUBLE,
                MPI_SUM, m_solver->m_rescalingCommGrComm, AT_, "varSliceLocal.begin()", "varSlice.begin()");

  if(globalTimeStep % 5 == 0 && m_solver->m_RKStep == 0
     && m_solver->domainId() == m_solver->m_rescalingCommGrRootGlobal) {
    cout << m_solver->domainId() << " ThetaInflow " << thetaGlobal(0) << " ThetaRecyclingStation " << thetaGlobal(1)
         << endl;
  }

  ///////////////////////////////////////////////////////////////
  ///////////////// RESCALING ///////////////////////////////////
  ///////////////////////////////////////////////////////////////

  for(MInt k = m_noGhostLayers; k < m_nCells[0] - m_noGhostLayers; ++k) {
    const MFloat ctem1 = (F1 + ctema) * (F1 - POW2(gams(k - m_noGhostLayers, i)));
    const MFloat ctem2 = F2 * ctema * gams(k - m_noGhostLayers, i) * (F1 - gams(k - m_noGhostLayers, i));
    const MFloat ctem3 = (F1 - gams(k - m_noGhostLayers, i))
                         * (F1 + gams(k - m_noGhostLayers, i) + F2 * ctema * gams(k - m_noGhostLayers, i));


    MInt jStart = 0;
    MInt jEnd = m_nCells[1];

    if(m_solver->m_nOffsetCells[1] == 0) {
      jStart = m_noGhostLayers;
    }
    if(m_solver->m_nOffsetCells[1] + m_solver->m_nActiveCells[1] == m_grid->getMyBlockNoCells(1)) {
      jEnd = m_nCells[1] - m_noGhostLayers;
    }

    for(MInt j = jStart; j < jEnd; ++j) {
      const MInt faceId = j + (k - m_noGhostLayers) * m_nCells[1];
      const MInt cellId = cellIndex(i, j, k);

      if(coordInOuter(faceId, i) < 1.05) {
        MFloat uInner = F0, vInner = F0, wInner = F0, TInner = F0, mutInner = F0;
        MFloat uOuter = F0, vOuter = F0, wOuter = F0, TOuter = F0, mutOuter = F0;
        const MFloat count = alpha * (coordInOuter(faceId, i) - b);
        const MFloat denom = (F1 - F2 * b) * coordInOuter(faceId, i) + b;
        const MFloat ratio = count / denom;
        const MFloat wfun = F1B2 * (F1 + tanh(ratio) / tanh(alpha));

        for(MInt jj = 0; jj < totalCells[1] - 1; ++jj) {
          const MInt localId = jj + (m_solver->m_nOffsetCells[0] + (k - m_noGhostLayers)) * totalCells[1];
          const MInt localIdP1 = jj + 1 + (m_solver->m_nOffsetCells[0] + (k - m_noGhostLayers)) * totalCells[1];

          const MFloat yInnerRe = varSlice(4, localId);
          const MFloat yInnerReP1 = varSlice(4, localIdP1);

          if((yInnerRe - coordInInner(faceId, i)) < rescalEPS && yInnerReP1 > coordInInner(faceId, i)) {
            const MFloat dy1 = coordInInner(faceId, i) - yInnerRe;
            const MFloat dy2 = yInnerReP1 - coordInInner(faceId, i);
            const MFloat dy = yInnerReP1 - yInnerRe;

            const MFloat u = varSlice(0, localId);
            const MFloat uP1 = varSlice(0, localIdP1);
            const MFloat v = varSlice(1, localId);
            const MFloat vP1 = varSlice(1, localIdP1);
            const MFloat w = varSlice(2, localId);
            const MFloat wP1 = varSlice(2, localIdP1);
            const MFloat t = varSlice(3, localId);
            const MFloat tP1 = varSlice(3, localIdP1);
            const MFloat mut = varSlice(6, localId);
            const MFloat mutP1 = varSlice(6, localIdP1);
            uInner = (uP1 * dy1 + u * dy2) / dy;
            vInner = (vP1 * dy1 + v * dy2) / dy;
            wInner = (wP1 * dy1 + w * dy2) / dy;
            TInner = (tP1 * dy1 + t * dy2) / dy;
            mutInner = (mutP1 * dy1 + mut * dy2) / dy;
          }
        }

        //>marian: catch those cells that didn't pass the if-clause (where TInner is still zero)
        if(TInner < 0.5) {
          for(MInt jj = 0; jj < totalCells[1] - 1; ++jj) {
            const MInt localId = jj + (m_solver->m_nOffsetCells[0] + (k - m_noGhostLayers)) * totalCells[1];
            const MInt localIdP1 = jj + 1 + (m_solver->m_nOffsetCells[0] + (k - m_noGhostLayers)) * totalCells[1];

            const MFloat yInnerRe = varSlice(4, localId);
            const MFloat yInnerReP1 = varSlice(4, localIdP1);
            const MFloat diffPercent = (abs(yInnerRe - coordInInner(faceId, i)) / coordInInner(faceId, i)) * 100.0;

            if((diffPercent < 10.0) && yInnerReP1 > coordInInner(faceId, i)) {
              const MFloat dy1 = coordInInner(faceId, i) - yInnerRe;
              const MFloat dy2 = yInnerReP1 - coordInInner(faceId, i);
              const MFloat dy = yInnerReP1 - yInnerRe;

              const MFloat u = varSlice(0, localId);
              const MFloat uP1 = varSlice(0, localIdP1);
              const MFloat v = varSlice(1, localId);
              const MFloat vP1 = varSlice(1, localIdP1);
              const MFloat w = varSlice(2, localId);
              const MFloat wP1 = varSlice(2, localIdP1);
              const MFloat t = varSlice(3, localId);
              const MFloat tP1 = varSlice(3, localIdP1);
              const MFloat mut = varSlice(6, localId);
              const MFloat mutP1 = varSlice(6, localIdP1);
              uInner = (uP1 * dy1 + u * dy2) / dy;
              vInner = (vP1 * dy1 + v * dy2) / dy;
              wInner = (wP1 * dy1 + w * dy2) / dy;
              TInner = (tP1 * dy1 + t * dy2) / dy;
              mutInner = (mutP1 * dy1 + mut * dy2) / dy;
            }
          }
        }
        //<marian

        // outer region
        for(MInt jj = 0; jj < totalCells[1] - 1; ++jj) {
          const MInt localId = jj + (m_solver->m_nOffsetCells[0] + (k - m_noGhostLayers)) * totalCells[1];
          const MInt localIdP1 = jj + 1 + (m_solver->m_nOffsetCells[0] + (k - m_noGhostLayers)) * totalCells[1];

          const MFloat yOuterRe = varSlice(5, localId);
          const MFloat yOuterReP1 = varSlice(5, localIdP1);

          if((yOuterRe - coordInOuter(faceId, i)) < rescalEPS && yOuterReP1 > coordInOuter(faceId, i)) {
            const MFloat dy1 = coordInOuter(faceId, i) - yOuterRe;
            const MFloat dy2 = yOuterReP1 - coordInOuter(faceId, i);
            const MFloat dy = yOuterReP1 - yOuterRe;


            const MFloat u = varSlice(0, localId);
            const MFloat uP1 = varSlice(0, localIdP1);
            const MFloat v = varSlice(1, localId);
            const MFloat vP1 = varSlice(1, localIdP1);
            const MFloat w = varSlice(2, localId);
            const MFloat wP1 = varSlice(2, localIdP1);
            const MFloat t = varSlice(3, localId);
            const MFloat tP1 = varSlice(3, localIdP1);
            const MFloat mut = varSlice(6, localId);
            const MFloat mutP1 = varSlice(6, localIdP1);
            uOuter = (uP1 * dy1 + u * dy2) / dy;
            vOuter = (vP1 * dy1 + v * dy2) / dy;
            wOuter = (wP1 * dy1 + w * dy2) / dy;
            TOuter = (tP1 * dy1 + t * dy2) / dy;
            mutOuter = (mutP1 * dy1 + mut * dy2) / dy;
          }
        }

        const MFloat TInnerA = POW2(gams(k - m_noGhostLayers, i)) * TInner + ctem1 * PV->TInfinity;
        const MFloat TOuterA =
            POW2(gams(k - m_noGhostLayers, i)) * TOuter - (ctem2 * (uOuter / PV->UInfinity) - ctem3) * PV->TInfinity;

        // van Driest transformation
        const MFloat uvdInner = PV->UInfinity * asin(b_vd * uInner / PV->UInfinity) / b_vd;
        const MFloat uvdOuter = PV->UInfinity * asin(b_vd * uOuter / PV->UInfinity) / b_vd;

        // scaling of transformed inner and outer velocities
        uInner = gams(k - m_noGhostLayers, i) * uvdInner;
        uOuter = gams(k - m_noGhostLayers, i) * uvdOuter + (F1 - gams(k - m_noGhostLayers, i)) * uvd8;

        uInner = PV->UInfinity * sin(b_vd * uInner / PV->UInfinity) / b_vd;
        uOuter = PV->UInfinity * sin(b_vd * uOuter / PV->UInfinity) / b_vd;

        const MFloat pres = PV->PInfinity;

        const MFloat uMean = uInner * (F1 - wfun) + uOuter * wfun;
        const MFloat vMean = vInner * (F1 - wfun) + vOuter * wfun;
        const MFloat wMean = (wInner * (F1 - wfun) + wOuter * wfun) * gams(k - m_noGhostLayers, i);
        const MFloat tMean = TInnerA * (F1 - wfun) + TOuterA * wfun;
        const MFloat rhoIn = gamma * pres / tMean;

        // turbulent viscosity
        const MFloat tempWallInlet = tempWallInletGlobal(thetaLocalOffset + (k - m_noGhostLayers));
        const MFloat tempWallRecycling = wallProperties(thetaLocalOffset + (k - m_noGhostLayers), 2);

        const MFloat viscWallInlet = SUTHERLANDLAW(tempWallInlet);
        const MFloat viscWallRecycling = SUTHERLANDLAW(tempWallRecycling);
        const MFloat thetaInlet = thetaGlobal(0);
        const MFloat thetaRecycling = thetaGlobal(1);
        mutInner = mutInner * (viscWallInlet / viscWallRecycling);
        mutOuter = mutOuter * gams(k - m_noGhostLayers, i) * (thetaInlet / thetaRecycling);

        MFloat mutMean = mutInner * (F1 - wfun) + mutOuter * wfun;
        const MFloat clebf = 6.6;
        const MFloat blt = m_rescalingBLT;
        const MFloat cleb = F1 / (F1 + pow((m_cells->coordinates[1][cellId] / (clebf * blt)), 6.0));

        m_cells->pvariables[PV->RHO][cellId] = rhoIn * cleb;
        m_cells->pvariables[PV->U][cellId] = uMean * cleb;
        m_cells->pvariables[PV->V][cellId] = vMean * cleb;
        m_cells->pvariables[PV->W][cellId] = wMean * cleb;

        m_cells->pvariables[PV->P][cellId] = pres;
        m_cells->pvariables[PV->RANS_VAR[0]][cellId] = mutMean / rhoIn;

      } else {
        // const MFloat pres = pressure(cellIndex(m_noGhostLayers,j,k)); //for supersonic PV->PInfinity
        const MFloat pres = PV->PInfinity;
        const MFloat rhoIn = gamma * pres / PV->TInfinity;

        const MFloat uMean = PV->UInfinity;
        const MFloat vMean = PV->VInfinity;
        const MFloat wMean = PV->WInfinity;

        m_cells->pvariables[PV->RHO][cellId] = rhoIn;
        m_cells->pvariables[PV->U][cellId] = uMean;
        m_cells->pvariables[PV->V][cellId] = vMean;
        m_cells->pvariables[PV->W][cellId] = wMean;
        m_cells->pvariables[PV->P][cellId] = pres;
        m_cells->pvariables[PV->RANS_VAR[0]][cellId] = PV->ransInfinity[0];
      }
    }
  }

  for(MInt k = 0; k < m_nCells[0]; ++k) {
    for(MInt j = 0; j < m_nCells[1]; ++j) {
      // extrapolation for second GC
      const MInt cellId = cellIndex(1, j, k);
      const MInt cellIdM1 = cellIndex(0, j, k);
      const MInt cellIdadj = cellIndex(2, j, k);

      for(MInt var = 0; var < PV->noVariables; var++) {
        m_cells->pvariables[var][cellIdM1] =
            2.0 * m_cells->pvariables[var][cellId] - m_cells->pvariables[var][cellIdadj];
      }
    }
  }
}


// Recycling station for RANS
template <MBool isRans>
void StructuredBndryCnd3D<isRans>::bc2511(MInt bcId) {
  MInt* start = m_physicalBCMap[bcId]->start1;
  cout.precision(8);

  // scaling between delta0 and delta2
  const MFloat F727 = 72.0 / 8.0;
  const MInt i = start[0]; // position at which recycle is taken
  const MFloat yWall = F0; // this has been fixed else method does not works

  //////////////////////////////////////////////////////////////
  ///////////////// COMPUTE THETA AND EXCHANGE /////////////////
  //////////////////////////////////////////////////////////////

  // thetaLocal.fill(F0); //initialize scratch space to zero // only for parallel use
  MFloatScratchSpace thetaLocal(2, AT_, "thetaLocalRe");
  MFloatScratchSpace thetaGlobal(2, AT_, "thetaGlobalRe");
  thetaLocal.fill(F0); // initialize scratch space
  thetaGlobal.fill(F0);

  // the offest position in k-direction is the offset
  const MInt thetaLocalOffset = m_solver->m_nOffsetCells[0];
  // compute the local moment thickness j=direction of integration
  for(MInt k = m_noGhostLayers; k < m_nCells[0] - m_noGhostLayers; ++k) {
    for(MInt j = m_noGhostLayers; j < m_nCells[1] - m_noGhostLayers; ++j) {
      const MInt cellId = cellIndex(i, j, k);
      const MInt pointIdM1 = getPointIdFromCell(i, j, k);
      const MInt pointIdP1 = getPointIdfromPoint(pointIdM1, 0, 1, 0);

      const MFloat urat = m_cells->pvariables[PV->U][cellId] / PV->UInfinity;
      const MFloat momThick =
          (m_cells->pvariables[PV->U][cellId] * m_cells->pvariables[PV->RHO][cellId] * fabs(F1 - urat))
          / (CV->rhoUInfinity);

      // integrate normal to the wall
      const MFloat ydist = m_grid->m_coordinates[1][pointIdP1] - m_grid->m_coordinates[1][pointIdM1];
      thetaLocal(1) += momThick * ydist;
    }
  }


  // communicate the Thickness across the plane
  MPI_Allreduce(thetaLocal.begin(), thetaGlobal.begin(), 2, MPI_DOUBLE, MPI_SUM, m_solver->m_rescalingCommGrComm, AT_,
                "thetaLocal.begin()", "thetaGlobal.begin()");

  thetaGlobal(1) /= m_grid->getMyBlockNoCells(0); // we have now the averaged momentum thi

  //////////////////////////////////////////////////////////////
  ///////////////// COMPUTE AND EXCHANGE WALL PROPERTIES ///////
  //////////////////////////////////////////////////////////////

  const MFloat delta = F727 * thetaGlobal(1);
  const MInt noVar = 3;                                     // for more variables if wanted
  const MInt wallLocalOffset = m_solver->m_nOffsetCells[1]; // Offset in j-direction

  MFloatScratchSpace wallPropertiesLocal(m_grid->getMyBlockNoCells(0), noVar, AT_, "wallPropertiesLocalRe");
  MFloatScratchSpace wallProperties(m_grid->getMyBlockNoCells(0), noVar, AT_, "wallPropertiesRe");
  wallPropertiesLocal.fill(F0);
  wallProperties.fill(F0);

  // determine the wall stuff if wall is contained whithin the partition
  if(wallLocalOffset == 0 && m_solver->m_nActiveCells[1] >= m_noGhostLayers) {
    for(MInt k = m_noGhostLayers; k < m_solver->m_nCells[0] - m_noGhostLayers; ++k) {
      const MInt cellId = cellIndex(i, m_noGhostLayers, k);
      const MFloat rho = m_cells->pvariables[PV->RHO][cellId];
      const MFloat t = temperature(cellId);
      const MFloat nu = SUTHERLANDLAW(t);
      const MFloat uWall = fabs(m_cells->pvariables[PV->U][cellId]);
      const MFloat ydist = m_cells->coordinates[1][cellId] - yWall;
      const MFloat uTau = sqrt(uWall * nu / (ydist * rho));
      const MInt localId = m_solver->m_nOffsetCells[0] + (k - m_noGhostLayers);

      wallPropertiesLocal(localId, 0) = uTau;
      wallPropertiesLocal(localId, 1) = rho;
      wallPropertiesLocal(localId, 2) = t;
    }
  }


  MPI_Allreduce(wallPropertiesLocal.begin(), wallProperties.begin(), noVar * m_grid->getMyBlockNoCells(0), MPI_DOUBLE,
                MPI_SUM, m_solver->m_rescalingCommGrComm, AT_, "wallPropertiesLocal.begin()", "wallProperties.begin()");

  MFloatScratchSpace tempWallInletLocal(m_grid->getMyBlockNoCells(0), AT_, "tempWallInletLocal");
  MFloatScratchSpace tempWallInletGlobal(m_grid->getMyBlockNoCells(0), AT_, "tempWallInletGlobal");
  tempWallInletLocal.fill(F0);
  tempWallInletGlobal.fill(F0);

  MPI_Allreduce(tempWallInletLocal.begin(), tempWallInletGlobal.begin(), m_grid->getMyBlockNoCells(0), MPI_DOUBLE,
                MPI_SUM, m_solver->m_rescalingCommGrComm, AT_, "tempWallInletLocal.begin()",
                "tempWallInletGlobal.begin()");

  //////////////////////////////////////////////////////////////
  ///////////////// COMPUTE AND EXCHANGE VAR SLICE /////////////
  //////////////////////////////////////////////////////////////

  MInt totalCells[2] = {m_grid->getMyBlockNoCells(0), m_grid->getMyBlockNoCells(1)};

  const MInt noVariables = PV->noVariables + 1;
  MFloatScratchSpace varSliceLocal(noVariables, totalCells[0] * totalCells[1], AT_, "varSliceLocal");
  MFloatScratchSpace varSlice(noVariables, totalCells[0] * totalCells[1], AT_, "varSlice");


  for(MInt k = m_noGhostLayers; k < m_nCells[0] - m_noGhostLayers; ++k) {
    for(MInt j = m_noGhostLayers; j < m_nCells[1] - m_noGhostLayers; ++j) {
      const MInt cellId = cellIndex(i, j, k);
      const MFloat rho = m_cells->pvariables[PV->RHO][cellId];
      const MFloat temp = temperature(cellId);
      const MFloat nu = SUTHERLANDLAW(temp);
      const MFloat uTauRe = wallProperties(thetaLocalOffset + (k - m_noGhostLayers), 0);
      const MFloat yIn = (m_cells->coordinates[1][cellId] - yWall) * uTauRe * rho / (nu * sqrt(m_solver->m_Re0));
      const MFloat yOut = (m_cells->coordinates[1][cellId] - yWall) * rho / (delta * CV->rhoInfinity);
      const MFloat u = m_cells->pvariables[PV->U][cellId];
      const MFloat v = m_cells->pvariables[PV->V][cellId];
      const MFloat w = m_cells->pvariables[PV->W][cellId];
      //>RANS
      const MFloat nut = m_cells->pvariables[PV->RANS_VAR[0]][cellId];
      //<RANS

      const MInt localId = m_solver->m_nOffsetCells[1] + (j - m_noGhostLayers)
                           + (m_solver->m_nOffsetCells[0] + (k - m_noGhostLayers)) * totalCells[1];

      // save the variables u,v,w,t,yI,yO
      varSliceLocal(0, localId) = u;
      varSliceLocal(1, localId) = v;
      varSliceLocal(2, localId) = w;
      varSliceLocal(3, localId) = temp;
      varSliceLocal(4, localId) = yIn;
      varSliceLocal(5, localId) = yOut;
      varSliceLocal(6, localId) = nut;
    }
  }

  // communicate the slice
  MPI_Allreduce(varSliceLocal.begin(), varSlice.begin(), noVariables * totalCells[0] * totalCells[1], MPI_DOUBLE,
                MPI_SUM, m_solver->m_rescalingCommGrComm, AT_, "varSliceLocal.begin()", "varSlice.begin()");
}

/**
 * \brief Laminar Poiseuille inflow
 *
 */
template <MBool isRans>
void StructuredBndryCnd3D<isRans>::bc2020(MInt bcId) {
  MInt* start = m_physicalBCMap[bcId]->start1;
  MInt* end = m_physicalBCMap[bcId]->end1;

  switch(m_physicalBCMap[bcId]->face) {
    case 0: {
      for(MInt k = start[2]; k < end[2]; k++) {
        for(MInt j = end[1] - 1; j >= start[1]; j--) {
          const MInt cellIdG1 = cellIndex(start[0] + 1, j, k);
          const MInt cellIdG2 = cellIndex(start[0] + 0, j, k);
          const MInt cellIdA1 = cellIndex(start[0] + 2, j, k);
          const MFloat y_max = F1; // channel height

          const MFloat x = m_cells->coordinates[0][cellIdG1];
          const MFloat y = m_cells->coordinates[1][cellIdG1];
          const MFloat pG1 =
              PV->PInfinity
              - F3 * (x + 15.0) * SUTHERLANDLAW(PV->TInfinity) * PV->UInfinity * POW2(F2 / y_max) / m_solver->m_Re0;

          m_cells->pvariables[PV->RHO][cellIdG1] = CV->rhoInfinity;
          m_cells->pvariables[PV->U][cellIdG1] =
              (-(F3 / F2) * PV->UInfinity * (POW2(y - y_max / F2) - POW2(y_max / F2)) / POW2(y_max / F2));
          m_cells->pvariables[PV->V][cellIdG1] = F0;
          m_cells->pvariables[PV->W][cellIdG1] = F0;
          m_cells->pvariables[PV->P][cellIdG1] = pG1;

          // extrapolate into second ghost cell
          for(MInt var = 0; var < PV->noVariables; var++) {
            m_cells->pvariables[var][cellIdG2] =
                F2 * m_cells->pvariables[var][cellIdG1] - m_cells->pvariables[var][cellIdA1];
          }
        }
      }
      break;
    }
    default: {
      cout << "bc2020: face not implemented" << endl;
    }
  }
}

/**
 * \brief Prescribe given profile BC
 *
 *  Precribes a profile from the restart file
 *  extrapolate pressure from computational domain
 *  \author Marian Albers
 */
template <MBool isRans>
void StructuredBndryCnd3D<isRans>::bc2600(MInt bcId) {
  MInt* start = m_physicalBCMap[bcId]->start1;
  MInt* end = m_physicalBCMap[bcId]->end1;

  switch(m_physicalBCMap[bcId]->face) {
    case 0: {
      for(MInt i = start[0]; i < end[0]; i++) {
        for(MInt j = start[1]; j < end[1]; j++) {
          for(MInt k = start[2]; k < end[2]; k++) {
            const MInt cellId = cellIndex(m_noGhostLayers - 1 - i, j, k);
            const MInt cellIdadj = cellIndex(m_noGhostLayers - i, j, k);

            // extrapolate pressure to ghost cells
            m_cells->pvariables[PV->P][cellId] = m_cells->pvariables[PV->P][cellIdadj];
          }
        }
      }
      break;
    }
    case 1: {
      const MInt IJK[nDim] = {1, m_nCells[2], m_nCells[1] * m_nCells[2]};
      const MFloat gamma = m_solver->m_gamma;

      // Here we find out the normal direction of the
      // boundary and the two tangential directions.
      // This way we can make a general formulation of
      // the boundary condition
      const MInt face = m_physicalBCMap[bcId]->face;
      const MInt normalDir = face / 2;
      const MInt firstTangentialDir = (normalDir + 1) % nDim;
      const MInt secondTangentialDir = (normalDir + 2) % nDim;
      const MInt normalDirStart = start[normalDir];
      const MInt firstTangentialStart = start[firstTangentialDir];
      const MInt firstTangentialEnd = end[firstTangentialDir];
      const MInt secondTangentialStart = start[secondTangentialDir];
      const MInt secondTangentialEnd = end[secondTangentialDir];
      const MInt inc[nDim] = {IJK[normalDir], IJK[firstTangentialDir], IJK[secondTangentialDir]};

      const MInt n = (face % 2) * 2 - 1;                                //-1,+1
      const MInt g1 = normalDirStart + (MInt)(0.5 - (0.5 * (MFloat)n)); //+1,0
      const MInt g2 = normalDirStart + (MInt)(0.5 + (0.5 * (MFloat)n)); // 0,+1
      const MInt a1 = normalDirStart + (MInt)(0.5 - (1.5 * (MFloat)n)); //+2,-1
      const MInt a2 = normalDirStart + (MInt)(0.5 - (2.5 * (MFloat)n)); //+3,-2

      for(MInt t1 = firstTangentialStart; t1 < firstTangentialEnd; t1++) {
        for(MInt t2 = secondTangentialStart; t2 < secondTangentialEnd; t2++) {
          const MInt cellIdG1 = g1 * inc[0] + t1 * inc[1] + t2 * inc[2];
          const MInt cellIdG2 = g2 * inc[0] + t1 * inc[1] + t2 * inc[2];
          const MInt cellIdA1 = a1 * inc[0] + t1 * inc[1] + t2 * inc[2];
          const MInt cellIdA2 = a2 * inc[0] + t1 * inc[1] + t2 * inc[2];
          const MInt cellIdBc =
              g1 - normalDirStart + (t1 + t2 * m_solver->m_bc2600noCells[1]) * m_solver->m_bc2600noCells[2];

          const MFloat dxidx = m_cells->surfaceMetrics[normalDir * nDim + 0][cellIdA1];
          const MFloat dxidy = m_cells->surfaceMetrics[normalDir * nDim + 1][cellIdA1];
          const MFloat dxidz = m_cells->surfaceMetrics[normalDir * nDim + 2][cellIdA1];
          // multiply with n, so it will be -1 or +1 depending if we enter
          // or leave the domain of integration in positive direction
          const MFloat gradxi = n * F1 / sqrt(dxidx * dxidx + dxidy * dxidy + dxidz * dxidz);

          const MFloat dxHelp = dxidx * gradxi;
          const MFloat dyHelp = dxidy * gradxi;
          const MFloat dzHelp = dxidz * gradxi;

          const MFloat cBC =
              sqrt(gamma * m_cells->pvariables[PV->P][cellIdG1] / m_cells->pvariables[PV->RHO][cellIdG1]);
          const MFloat rhoBC = m_cells->pvariables[PV->RHO][cellIdG1];

          const MFloat rhoInner = m_cells->pvariables[PV->RHO][cellIdA1];
          const MFloat uInner = m_cells->pvariables[PV->U][cellIdA1];
          const MFloat vInner = m_cells->pvariables[PV->V][cellIdA1];
          const MFloat wInner = m_cells->pvariables[PV->W][cellIdA1];
          const MFloat pInner = m_cells->pvariables[PV->P][cellIdA1];

          const MFloat u2600 = m_solver->m_bc2600Variables[PV->U][cellIdBc];
          const MFloat v2600 = m_solver->m_bc2600Variables[PV->V][cellIdBc];
          const MFloat w2600 = m_solver->m_bc2600Variables[PV->W][cellIdBc];
          const MFloat rho2600 = m_solver->m_bc2600Variables[PV->RHO][cellIdBc];
          const MFloat p2600 = m_solver->m_bc2600Variables[PV->P][cellIdBc];
          const MFloat rans2600 = m_solver->m_bc2600Variables[PV->RANS_VAR[0]][cellIdBc];

          const MFloat maContravariant =
              (dxidx * uInner + dxidy * vInner + dxidz * wInner - m_cells->dxt[normalDir][cellIdA1]) * gradxi;

          if(maContravariant < F0) {
            // inflow
            const MFloat p =
                F1B2
                * (pInner + p2600
                   + rhoBC * cBC * (dxHelp * (uInner - u2600) + dyHelp * (vInner - v2600) + dzHelp * (wInner - w2600)));

            const MFloat rho = rho2600 + (p - p2600) / POW2(cBC);
            const MFloat help = (p - p2600) / (rhoBC * cBC);

            m_cells->pvariables[PV->RHO][cellIdG1] = rho;
            m_cells->pvariables[PV->U][cellIdG1] = (u2600 + help * dxHelp);
            m_cells->pvariables[PV->V][cellIdG1] = (u2600 + help * dyHelp);
            m_cells->pvariables[PV->W][cellIdG1] = (u2600 + help * dzHelp);
            m_cells->pvariables[PV->P][cellIdG1] = p;
            if(isRans) {
              m_cells->pvariables[PV->RANS_VAR[0]][cellIdG1] = rans2600;
            }
          } else {
            // outflow
            MFloat p = p2600; // PV->PInfinity;
            // if(m_solver->m_restartTimeStep+globalTimeStep <= 1000000) {
            //   p = PV->PInfinity +
            //   (p2600-PV->PInfinity)*(((MFloat)(m_solver->m_restartTimeStep+globalTimeStep))/1000000.0);
            // } else {
            //   p = p2600;
            // }

            const MFloat rho = rhoInner + (p - pInner) / POW2(cBC);
            const MFloat help = (p - pInner) / (rhoBC * cBC);

            m_cells->pvariables[PV->RHO][cellIdG1] = rho;
            m_cells->pvariables[PV->U][cellIdG1] = (uInner - help * dxHelp);
            m_cells->pvariables[PV->V][cellIdG1] = (vInner - help * dyHelp);
            m_cells->pvariables[PV->W][cellIdG1] = (wInner - help * dzHelp);
            m_cells->pvariables[PV->P][cellIdG1] = p;
            if(isRans) {
              m_cells->pvariables[PV->RANS_VAR[0]][cellIdG1] = (F2 * m_cells->pvariables[PV->RANS_VAR[0]][cellIdA1]
                                                                - m_cells->pvariables[PV->RANS_VAR[0]][cellIdA2]);
            }
          }

          // extrapolate into second ghost cell
          for(MInt var = 0; var < PV->noVariables; var++) {
            m_cells->pvariables[var][cellIdG2] =
                F2 * m_cells->pvariables[var][cellIdG1] - m_cells->pvariables[var][cellIdA1];
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


/**
 * \brief Prescribe given profile BC
 *
 *  Precribes a profile from the restart file
 *  extrapolate pressure from computational domain
 *  \author Marian Albers 2016
 */
template <MBool isRans>
void StructuredBndryCnd3D<isRans>::bc2601(MInt bcId) {
  MInt* start = m_physicalBCMap[bcId]->start1;
  MInt* end = m_physicalBCMap[bcId]->end1;

  MFloat theta100000 = 0.002065747418;
  MFloat rexPosQuad = 0.014452;
  MFloat gammaEpsilon = m_solver->m_bc2601GammaEpsilon;

  if(m_solver->m_rans) {
    theta100000 = 0.002171197777811;
    rexPosQuad = 0.14452;
  }

  const MFloat lengthFactor = theta100000 * rexPosQuad;

  MFloatScratchSpace effConst(m_2601noPos, m_2601noCoeff, AT_, "effConst");
  effConst.fill(F0);


  switch(m_physicalBCMap[bcId]->face) {
    case 2: {
      for(MInt k = start[2]; k < end[2]; k++) {
        ////////////////////////////////////////////////////////////
        /////////// TIME INTERPOLATION ACTUATED SOLUTION ///////////
        ////////////////////////////////////////////////////////////
        // manipulate time with spanwise coordinate
        const MInt spanCell = cellIndex(m_noGhostLayers, m_noGhostLayers, k);
        const MFloat spanCoord = m_cells->coordinates[2][spanCell];
        const MFloat waveSpeed = 0.015795;
        const MFloat timeShift = spanCoord / waveSpeed;
        const MFloat t = m_solver->m_time + m_solver->m_timeStep * m_solver->m_RKalpha[m_solver->m_RKStep] - timeShift;

        //////////////////////////////////////////
        /////////// TIME INTERPOLATION ///////////
        //////////////////////////////////////////

        // interpolate between time values
        for(MInt pos = 0; pos < m_2601noPos; pos++) {
          effConst(pos, 0) = m_2601streamwisePos[pos];
          for(MInt iter = 0; iter < m_2601noSteps - 1; iter++) {
            if(m_2601effConst[0][pos * m_2601noCoeff + 0] >= t) {
              effConst(pos, 1) = m_2601effConst[0][pos * m_2601noCoeff + 1];
              effConst(pos, 2) = m_2601effConst[0][pos * m_2601noCoeff + 2];
              effConst(pos, 3) = m_2601effConst[0][pos * m_2601noCoeff + 3];
              effConst(pos, 4) = m_2601effConst[0][pos * m_2601noCoeff + 4];
              break;
            } else if(m_2601effConst[iter][pos * m_2601noCoeff + 0] < t
                      && t <= m_2601effConst[iter + 1][pos * m_2601noCoeff + 0]) {
              for(MInt var = 1; var < 5; var++) {
                // linear interpolation between two time values
                effConst(pos, var) = m_2601effConst[iter][pos * m_2601noCoeff + var]
                                     + (t - m_2601effConst[iter][pos * m_2601noCoeff + 0])
                                           / (m_2601effConst[iter + 1][pos * m_2601noCoeff + 0]
                                              - m_2601effConst[iter][pos * m_2601noCoeff + 0])
                                           * (m_2601effConst[iter + 1][pos * m_2601noCoeff + var]
                                              - m_2601effConst[iter][pos * m_2601noCoeff + var]);
              }
              break;
            } else {
              effConst(pos, 1) = m_2601effConst[iter][pos * m_2601noCoeff + 1];
              effConst(pos, 2) = m_2601effConst[iter][pos * m_2601noCoeff + 2];
              effConst(pos, 3) = m_2601effConst[iter][pos * m_2601noCoeff + 3];
              effConst(pos, 4) = m_2601effConst[iter][pos * m_2601noCoeff + 4];
            }
          }
        }

        for(MInt pos = 0; pos < m_2601noPos; pos++) {
          effConst(pos, 0) = (effConst(pos, 0) - rexPosQuad) / lengthFactor + 36.57;
        }

        //////////////////////////////////////////
        /////////// SPACE INTERPOLATION //////////
        //////////////////////////////////////////

        for(MInt i = start[0]; i < end[0]; i++) {
          const MInt cellId = cellIndex(i, 1, k);

          MInt cellIdG2 = cellIndex(i, 0, k); // ghost
          MInt cellIdG1 = cellIndex(i, 1, k); // ghost
          MInt cellIdA1 = cellIndex(i, 2, k); // field
          MInt cellIdA2 = cellIndex(i, 3, k); // field

          const MInt cellIdBcG1 = i + (1 + k * m_noGhostLayers) * m_solver->m_nCells[2];
          const MInt cellIdBcG2 = i + (0 + k * m_noGhostLayers) * m_solver->m_nCells[2];

          const MFloat rho_zerothG1 = m_solver->m_bc2601ZerothOrderSolution[PV->RHO][cellIdBcG1];
          const MFloat u_zerothG1 = m_solver->m_bc2601ZerothOrderSolution[PV->U][cellIdBcG1];
          const MFloat v_zerothG1 = m_solver->m_bc2601ZerothOrderSolution[PV->V][cellIdBcG1];
          const MFloat w_zerothG1 = m_solver->m_bc2601ZerothOrderSolution[PV->W][cellIdBcG1];

          // const MFloat rho_zerothG2 = m_solver->m_bc2601ZerothOrderSolution[PV->RHO][cellIdBcG2];
          const MFloat u_zerothG2 = m_solver->m_bc2601ZerothOrderSolution[PV->U][cellIdBcG2];
          const MFloat v_zerothG2 = m_solver->m_bc2601ZerothOrderSolution[PV->V][cellIdBcG2];
          const MFloat w_zerothG2 = m_solver->m_bc2601ZerothOrderSolution[PV->W][cellIdBcG2];

          const MFloat uInfQuad = 1.029342800042823e+02;
          const MFloat rhoInfQuad = 1.209633989123;
          // const MFloat heightG1 = m_cells->coordinates[1][cellIdG1];
          // const MFloat heightG2 = m_cells->coordinates[1][cellIdG2];
          const MFloat xPos = m_cells->coordinates[0][cellId];

          MFloat uFactor = F0;
          MFloat vFactor = F0;
          MFloat wFactor = F0;
          MFloat rhoFactor = F0;

          // interpolated from the actuated solution
          MFloat uActuated = F0, vActuated = F0, wActuated = F0;

          // space interpolation
          if(xPos < effConst(0, 0)) {
            const MFloat relPos = (xPos - effConst(0, 0)) / (effConst(1, 0) - effConst(0, 0));
            uFactor = effConst(0, 1) + relPos * (effConst(1, 1) - effConst(0, 1));
            vFactor = effConst(0, 2) + relPos * (effConst(1, 2) - effConst(0, 2));
            wFactor = effConst(0, 3) + relPos * (effConst(1, 3) - effConst(0, 3));
            rhoFactor = effConst(0, 4) + relPos * (effConst(1, 4) - effConst(0, 4));
          } else if(effConst(0, 0) <= xPos && xPos < effConst(m_2601noPos - 1, 0)) {
            for(MInt pos = 0; pos < m_2601noPos - 1; pos++) {
              if(effConst(pos, 0) <= xPos && xPos < effConst(pos + 1, 0)) {
                const MFloat relPos = (xPos - effConst(pos, 0)) / (effConst(pos + 1, 0) - effConst(pos, 0));
                uFactor = effConst(pos, 1) + relPos * (effConst(pos + 1, 1) - effConst(pos, 1));
                vFactor = effConst(pos, 2) + relPos * (effConst(pos + 1, 2) - effConst(pos, 2));
                wFactor = effConst(pos, 3) + relPos * (effConst(pos + 1, 3) - effConst(pos, 3));
                rhoFactor = effConst(pos, 4) + relPos * (effConst(pos + 1, 4) - effConst(pos, 4));
                break;
              }
            }
          } else {
            const MFloat relPos =
                (xPos - effConst(m_2601noPos - 2, 0)) / (effConst(m_2601noPos - 1, 0) - effConst(m_2601noPos - 2, 0));
            uFactor =
                effConst(m_2601noPos - 2, 1) + relPos * (effConst(m_2601noPos - 1, 1) - effConst(m_2601noPos - 2, 1));
            vFactor =
                effConst(m_2601noPos - 2, 2) + relPos * (effConst(m_2601noPos - 1, 2) - effConst(m_2601noPos - 2, 2));
            wFactor =
                effConst(m_2601noPos - 2, 3) + relPos * (effConst(m_2601noPos - 1, 3) - effConst(m_2601noPos - 2, 3));
            rhoFactor =
                effConst(m_2601noPos - 2, 4) + relPos * (effConst(m_2601noPos - 1, 4) - effConst(m_2601noPos - 2, 4));
          }

          // //for now take uniform distribution
          // uFactor   = effConst(0,1);
          // vFactor   = effConst(0,2);
          // rhoFactor = effConst(0,4);

          // MFloat u_correctedG1 = u_zerothG1 + heightG1*lengthFactor*(uFactor)*(PV->UInfinity/uInfQuad);
          // MFloat u_correctedG2 = u_zerothG2 + heightG2*lengthFactor*(uFactor)*(PV->UInfinity/uInfQuad);


          MFloat u_zerothBC = u_zerothG1 + F1B2 * (u_zerothG1 - u_zerothG2);
          MFloat u_correctedBC = u_zerothBC + gammaEpsilon * lengthFactor * (uFactor) * (PV->UInfinity / uInfQuad);
          MFloat v_zerothBC = v_zerothG1 + F1B2 * (v_zerothG1 - v_zerothG2);
          MFloat v_correctedBC = v_zerothBC + gammaEpsilon * lengthFactor * (vFactor) * (PV->UInfinity / uInfQuad);
          MFloat w_zerothBC = w_zerothG1 + F1B2 * (w_zerothG1 - w_zerothG2);
          MFloat w_correctedBC =
              w_zerothBC + gammaEpsilon * lengthFactor * (wFactor) * (PV->UInfinity / uInfQuad); // + wActuated;
          MFloat rho_correctedG1 =
              rho_zerothG1 + gammaEpsilon * lengthFactor * (rhoFactor) * (CV->rhoInfinity / rhoInfQuad);

          if(m_solver->m_nOffsetCells[2] + i - 2 == 80 && m_solver->m_nOffsetCells[0] + k - 2 == 50) {
            if(globalTimeStep % 20 == 0 && m_solver->m_RKStep == 0) {
              cout.precision(7);
              cout << "globalTimeStep: " << globalTimeStep << " time: " << t << " gammaEpsilon: " << gammaEpsilon
                   << " uZeroth: " << u_zerothBC << " u_correctedBC: " << u_correctedBC << " vZeroth: " << v_zerothBC
                   << " v_correctedBC: " << v_correctedBC << " rhoZeroth: " << rho_zerothG1
                   << " rho_correctedG1: " << rho_correctedG1 << " vActuated: " << vActuated << endl;


              FILE* f_effective;
              f_effective = fopen("./effective_boundary.dat", "a+");
              fprintf(f_effective, "%d", globalTimeStep);
              fprintf(f_effective, " %f", m_solver->m_physicalTime);
              fprintf(f_effective, " %f", m_solver->m_time);
              fprintf(f_effective, " %f", m_solver->m_timeStep);
              fprintf(f_effective, " %f", u_zerothBC);
              fprintf(f_effective, " %f", u_correctedBC);
              fprintf(f_effective, " %f", v_zerothBC);
              fprintf(f_effective, " %f", v_correctedBC);
              fprintf(f_effective, " %f", w_zerothBC);
              fprintf(f_effective, " %f", w_correctedBC);
              fprintf(f_effective, " %f", rho_zerothG1);
              fprintf(f_effective, " %f", rho_correctedG1);
              fprintf(f_effective, " %f", uActuated);
              fprintf(f_effective, " %f", vActuated);
              fprintf(f_effective, " %f", wActuated);
              fprintf(f_effective, "\n");
              fclose(f_effective);
            }
          }

          const MFloat u_uncorrectedBC = u_zerothG1 + F1B2 * (u_zerothG1 - u_zerothG2);
          const MFloat v_uncorrectedBC = v_zerothG1 + F1B2 * (v_zerothG1 - v_zerothG2);
          const MFloat w_uncorrectedBC = w_zerothG1 + F1B2 * (w_zerothG1 - w_zerothG2);
          const MFloat rho_uncorrectedG1 = rho_zerothG1;

          MFloat u_appliedBC = F0, v_appliedBC = F0, w_appliedBC = F0, rho_appliedG1 = F0;

          if(m_cells->coordinates[0][cellId] < 12.4) {
            u_appliedBC = u_uncorrectedBC;
            v_appliedBC = v_uncorrectedBC;
            w_appliedBC = w_uncorrectedBC;
            rho_appliedG1 = rho_uncorrectedG1;
          } else if(m_cells->coordinates[0][cellId] >= 12.4 && m_cells->coordinates[0][cellId] < 15.4) {
            const MFloat fader = (m_cells->coordinates[0][cellId] - 12.4) / 3.0;
            u_appliedBC = u_uncorrectedBC * (1.0 - fader) + u_correctedBC * fader;
            v_appliedBC = v_uncorrectedBC * (1.0 - fader) + v_correctedBC * fader;
            w_appliedBC = w_uncorrectedBC * (1.0 - fader) + w_correctedBC * fader;
            rho_appliedG1 = rho_uncorrectedG1 * (1.0 - fader) + rho_correctedG1 * fader;
          } else if(m_cells->coordinates[0][cellId] >= 15.4 && m_cells->coordinates[0][cellId] < 60.0) {
            u_appliedBC = u_correctedBC;
            v_appliedBC = v_correctedBC;
            w_appliedBC = w_correctedBC;
            rho_appliedG1 = rho_correctedG1;
          } else if(m_cells->coordinates[0][cellId] >= 60.0 && m_cells->coordinates[0][cellId] < 63.0) {
            const MFloat fader = (m_cells->coordinates[0][cellId] - 60.0) / 3.0;
            u_appliedBC = u_uncorrectedBC * (fader) + u_correctedBC * (1.0 - fader);
            v_appliedBC = v_uncorrectedBC * (fader) + v_correctedBC * (1.0 - fader);
            w_appliedBC = w_uncorrectedBC * (fader) + w_correctedBC * (1.0 - fader);
            rho_appliedG1 = rho_uncorrectedG1 * (fader) + rho_correctedG1 * (1.0 - fader);
          } else {
            u_appliedBC = u_uncorrectedBC;
            v_appliedBC = v_uncorrectedBC;
            w_appliedBC = w_uncorrectedBC;
            rho_appliedG1 = rho_uncorrectedG1;
          }

          // if(m_cells->coordinates[0][cellId] < 15.4) {
          //   u_appliedBC = u_uncorrectedBC;
          //   v_appliedBC = v_uncorrectedBC;
          //   w_appliedBC = w_uncorrectedBC;
          //   rho_appliedG1 = rho_uncorrectedG1;
          // } else if(m_cells->coordinates[0][cellId] >= 15.4 && m_cells->coordinates[0][cellId] < 60.0) {
          //   u_appliedBC = u_correctedBC;
          //   v_appliedBC = v_correctedBC;
          //   w_appliedBC = w_correctedBC;
          //   rho_appliedG1 = rho_correctedG1;
          // } else {
          //   u_appliedBC = u_uncorrectedBC;
          //   v_appliedBC = v_uncorrectedBC;
          //   w_appliedBC = w_uncorrectedBC;
          //   rho_appliedG1 = rho_uncorrectedG1;
          // }

          // const MFloat rho_correctedG2 = rho_zerothG2 +
          // heightG2*lengthFactor*(rhoFactor)*(CV->rhoInfinity/rhoInfQuad); const MFloat rho_correctedBC =
          // rho_correctedG1 + F1B2*(rho_correctedG1-rho_correctedG2);
          const MFloat pField = pressure(cellIndex(i, m_noGhostLayers, k));


          const MFloat u1 = m_cells->pvariables[PV->U][cellIdA1];
          const MFloat v1 = m_cells->pvariables[PV->V][cellIdA1];
          const MFloat w1 = m_cells->pvariables[PV->W][cellIdA1];
          const MFloat rho1 = m_cells->pvariables[PV->RHO][cellIdA1];

          const MFloat u2 = m_cells->pvariables[PV->U][cellIdA2];
          const MFloat v2 = m_cells->pvariables[PV->V][cellIdA2];
          const MFloat w2 = m_cells->pvariables[PV->W][cellIdA2];
          // MFloat rho2 = m_cells->pvariables[PV->RHO][cellIdA2];

          const MFloat uG1 = F2 * u_appliedBC - u1;
          const MFloat vG1 = F2 * v_appliedBC - v1;
          const MFloat wG1 = F2 * w_appliedBC - w1;

          const MFloat uG2 = F2 * u_appliedBC - u2;
          const MFloat vG2 = F2 * v_appliedBC - v2;
          const MFloat wG2 = F2 * w_appliedBC - w2;

          const MFloat rhoG1 = rho_appliedG1; // F2*rho_correctedBC-rho1;
          const MFloat rhoG2 = F2 * rho_appliedG1 - rho1;

          m_cells->pvariables[PV->RHO][cellIdG1] = rhoG1;
          m_cells->pvariables[PV->U][cellIdG1] = uG1;
          m_cells->pvariables[PV->V][cellIdG1] = vG1;
          m_cells->pvariables[PV->W][cellIdG1] = wG1;
          m_cells->pvariables[PV->P][cellIdG1] = pField;

          m_cells->pvariables[PV->RHO][cellIdG2] = rhoG2;
          m_cells->pvariables[PV->U][cellIdG2] = uG2;
          m_cells->pvariables[PV->V][cellIdG2] = vG2;
          m_cells->pvariables[PV->W][cellIdG2] = wG2;
          m_cells->pvariables[PV->P][cellIdG2] = pField;

          if(m_solver->m_rans) {
            const MFloat nutilde_zerothG1 = m_solver->m_bc2601ZerothOrderSolution[PV->RANS_VAR[0]][cellIdBcG1];
            const MFloat nutilde_zerothG2 = m_solver->m_bc2601ZerothOrderSolution[PV->RANS_VAR[0]][cellIdBcG2];
            m_cells->pvariables[PV->RANS_VAR[0]][cellIdG1] = nutilde_zerothG1;
            m_cells->pvariables[PV->RANS_VAR[0]][cellIdG2] = nutilde_zerothG2;
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

/**
 * \brief supersonic inflow with imposed acoustic or entropy waves
 * \authors Thomas Schilden, Leo Hoening
 * \date 2015
 */
template <MBool isRans>
void StructuredBndryCnd3D<isRans>::bc2700(MInt bcId) {
  MInt* start = m_physicalBCMap[bcId]->start1;
  MInt* end = m_physicalBCMap[bcId]->end1;
  const MFloat gamma = m_solver->m_gamma;
  const MFloat time = m_solver->m_time;
  MFloat velocity[3] = {F0, F0, F0};

  switch(m_physicalBCMap[bcId]->face) {
    case 0:
    case 3: {
      for(MInt k = start[2]; k < end[2]; k++) {
        for(MInt j = start[1]; j < end[1]; j++) {
          for(MInt i = start[0]; i < end[0]; i++) {
            MInt cellId = cellIndex(i, j, k);
            MFloat rho = CV->rhoInfinity;
            for(MInt dim = 0; dim < nDim; dim++) {
              velocity[dim] = PV->VVInfinity[dim];
            }
            MFloat p = PV->PInfinity;

            // add the modes
            for(MInt mode = 0; mode < m_modes; mode++) {
              // 1. pressure
              const MFloat pressure_f =
                  gamma * m_solver->m_Ma * PV->PInfinity * m_modeAmp[mode] * (MFloat)(m_modeType[mode]);

              // 2. density
              MFloat density_f;
              const MFloat a = sqrt(PV->TInfinity);
              if(m_modeType[mode])
                density_f = pressure_f / POW2(a);
              else
                density_f = m_modeAmp[mode] * CV->rhoInfinity * m_solver->m_Ma;

              // 3. velocity
              MFloat K = F0;
              for(MInt dim = 0; dim < nDim; dim++) {
                K += pow(m_modeK[mode][dim], 2);
              }

              K = sqrt(K);
              const MFloat acImp = a * CV->rhoInfinity;
              const MFloat modeVelocity = (MFloat)(m_modeType[mode]) * pressure_f / acImp;
              MFloat velocity_f[3];
              for(MInt dim = 0; dim < nDim; dim++) {
                velocity_f[dim] = (m_modeK[mode][dim] * modeVelocity) / K;
              }

              // 1. calculate cell dependant trigonometry
              MFloat trigTerm = F0;
              for(MInt dim = 0; dim < nDim; dim++) {
                trigTerm += m_modeK[mode][dim] * m_cells->coordinates[dim][cellId];
              }

              trigTerm -= m_modeOmega[mode] * (time); // time;
              trigTerm += m_modePhi[mode];
              trigTerm = sin(trigTerm);

              rho += trigTerm * density_f;
              for(MInt dim = 0; dim < nDim; dim++) {
                velocity[dim] += trigTerm * velocity_f[0];
              }
              p += trigTerm * pressure_f;
            }

            // conservatives:
            m_cells->pvariables[PV->RHO][cellId] = rho;
            m_cells->pvariables[PV->U][cellId] = velocity[PV->U];
            m_cells->pvariables[PV->V][cellId] = velocity[PV->V];
            m_cells->pvariables[PV->W][cellId] = velocity[PV->W];
            m_cells->pvariables[PV->P][cellId] = p;
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


// Jet Freund Inlet
template <MBool isRans>
void StructuredBndryCnd3D<isRans>::bc2900(MInt bcId) {
  MInt* start = m_physicalBCMap[bcId]->start1;
  MInt* end = m_physicalBCMap[bcId]->end1;
  switch(m_physicalBCMap[bcId]->face) {
    case 0: {
      MInt cellId = -1;
      MInt cellIdadj = -1;
      MInt pIJK = 0, pIJPK = 0, pIJKP = 0, pIJPKP = 0;
      // fully new bc
      MInt i = end[0] - 1;
      MFloat pBC = F0, rho = F0, u = F0, v = F0, w = F0;
      MFloat drho = F0, du = F0, dv = F0, dw = F0, dp = F0;
      MFloat xBC = F0;
      // pBC=PV->PInfinity;
      MFloat pInner = F0, c02 = F0, distance = F0;
      for(MInt k = start[2]; k < end[2]; k++) {
        for(MInt j = start[1]; j < end[1]; j++) {
          cellId = cellIndex(i, j, k);
          cellIdadj = cellIndex(i + 1, j, k);
          // to determine the face coordinates!!!!!!
          pIJK = getPointIdFromCell(i + 1, j, k);
          pIJPK = getPointIdfromPoint(pIJK, 0, 1, 0);
          pIJKP = getPointIdfromPoint(pIJK, 0, 0, 1);
          pIJPKP = getPointIdfromPoint(pIJK, 0, 1, 1);
          MFloat r = sqrt(POW2(m_cells->coordinates[0][cellId]) + POW2(m_cells->coordinates[1][cellId])
                          + POW2(m_cells->coordinates[2][cellId]));
          MFloat uInf = F1B2 * (F1 - tanh(12.5 * (fabs(r / 0.5) - fabs(0.5 / r)))) * PV->VVInfinity[0];
          MFloat vInf = 0;
          MFloat wInf = 0;

          // values at the inner point
          pInner = m_cells->pvariables[PV->P][cellIdadj];
          c02 = sqrt(m_solver->m_gamma * pInner / m_cells->pvariables[PV->RHO][cellIdadj]);
          u = m_cells->pvariables[PV->U][cellIdadj];
          v = m_cells->pvariables[PV->V][cellIdadj];
          w = m_cells->pvariables[PV->W][cellIdadj];

          MFloat dxidx = m_cells->surfaceMetrics[0][cellId];
          MFloat dxidy = m_cells->surfaceMetrics[1][cellId];
          MFloat dxidz = m_cells->surfaceMetrics[2][cellId];

          // values at the boundary
          pBC = F1B2
                * (pInner + PV->PInfinity
                   - m_cells->pvariables[PV->RHO][cellIdadj] * c02
                         * (dxidx * (uInf - u) + dxidy * (vInf - v) + dxidz * (wInf - w)));
          rho = CV->rhoInfinity + ((pBC - PV->PInfinity) / (c02 * c02));

          u = uInf - dxidx * (PV->PInfinity - pBC) / (m_cells->pvariables[PV->RHO][cellIdadj] * c02);
          v = vInf - dxidy * (PV->PInfinity - pBC) / (m_cells->pvariables[PV->RHO][cellIdadj] * c02);
          w = wInf - dxidz * (PV->PInfinity - pBC) / (m_cells->pvariables[PV->RHO][cellIdadj] * c02);

          // extrapolate the variables into the ghost cells
          // gradients

          xBC = F1B4
                * (m_grid->m_coordinates[0][pIJK] + m_grid->m_coordinates[0][pIJPK] + m_grid->m_coordinates[0][pIJKP]
                   + m_grid->m_coordinates[0][pIJPKP]);


          distance = (xBC - m_cells->coordinates[0][cellIdadj]);

          drho = (rho - m_cells->pvariables[PV->RHO][cellIdadj]) / distance;
          du = (u - m_cells->pvariables[PV->U][cellIdadj]) / distance;
          dv = (v - m_cells->pvariables[PV->V][cellIdadj]) / distance;
          dw = (w - m_cells->pvariables[PV->W][cellIdadj]) / distance;
          dp = (pBC - m_cells->pvariables[PV->P][cellIdadj]) / distance;

          // extrapolate:
          for(MInt ii = start[0]; ii < end[0]; ++ii) {
            cellId = cellIndex(ii, j, k);
            distance = (m_cells->coordinates[0][cellId] - m_cells->coordinates[0][cellIdadj]);
            m_cells->pvariables[PV->RHO][cellId] = m_cells->pvariables[PV->RHO][cellIdadj] + drho * distance;
            m_cells->pvariables[PV->U][cellId] = m_cells->pvariables[PV->U][cellIdadj] + du * distance;
            m_cells->pvariables[PV->V][cellId] = m_cells->pvariables[PV->V][cellIdadj] + dv * distance;
            m_cells->pvariables[PV->W][cellId] = m_cells->pvariables[PV->W][cellIdadj] + dw * distance;
            m_cells->pvariables[PV->P][cellId] = m_cells->pvariables[PV->P][cellIdadj] + dp * distance;
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
void StructuredBndryCnd3D<isRans>::bc6000(MInt bcId) {
  cout << "applying bc " << bcId << endl;
}

// Works equally for conversion of fluid state into porous and vice versa
// Only works for LES
template <MBool isRans>
void StructuredBndryCnd3D<isRans>::bc6002(MInt bcId) {
  TRACE();

  const MInt IJK[nDim] = {1, m_nCells[2], m_nCells[1] * m_nCells[2]};
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
    const MInt secondTangentialDir = (normalDir + 2) % nDim;
    const MInt normalDirStart = start[normalDir];
    const MInt firstTangentialStart = start[firstTangentialDir];
    const MInt firstTangentialEnd = end[firstTangentialDir];
    const MInt secondTangentialStart = start[secondTangentialDir];
    const MInt secondTangentialEnd = end[secondTangentialDir];
    const MInt inc[nDim] = {IJK[normalDir], IJK[firstTangentialDir], IJK[secondTangentialDir]};

    // determine indices for direction help
    const MInt n = (face % 2) * 2 - 1;                                    //-1,+1
    const MInt g[2] = {normalDirStart + (MInt)(0.5 - (0.5 * (MFloat)n)),  //+1,0
                       normalDirStart + (MInt)(0.5 + (0.5 * (MFloat)n))}; // 0,+1
    const MInt a1 = normalDirStart + (MInt)(0.5 - (1.5 * (MFloat)n));     //+2,-1

    // convention: index x is unknown and y is known
    for(MInt t1 = firstTangentialStart; t1 < firstTangentialEnd; t1++) {
      for(MInt t2 = secondTangentialStart; t2 < secondTangentialEnd; t2++) {
        const MInt cellIdG1 = g[0] * inc[0] + t1 * inc[1] + t2 * inc[2];
        const MFloat normalVec[nDim] = {m_cells->fq[FQ->NORMAL[0]][cellIdG1], m_cells->fq[FQ->NORMAL[1]][cellIdG1],
                                        m_cells->fq[FQ->NORMAL[2]][cellIdG1]};
        const MInt cellIdA1 = a1 * inc[0] + t1 * inc[1] + t2 * inc[2];
        const MFloat por_x = m_cells->fq[FQ->POROSITY][cellIdA1];
        for(MInt i = 0; i < 2 /*m_noGhostLayers*/; ++i) {
          const MInt cellIdG = g[i] * inc[0] + t1 * inc[1] + t2 * inc[2];

          const MFloat por_y = m_cells->fq[FQ->POROSITY][cellIdG];
          const MFloat rho_y = m_cells->pvariables[PV->RHO][cellIdG];
          const MFloat p_y = m_cells->pvariables[PV->P][cellIdG];
          const MFloat u_y = m_cells->pvariables[PV->U][cellIdG];
          const MFloat v_y = m_cells->pvariables[PV->V][cellIdG];
          const MFloat w_y = m_cells->pvariables[PV->W][cellIdG];

          // Compute normal and tangential velocity components
          const MFloat U_y = u_y * normalVec[0] + v_y * normalVec[1] + w_y * normalVec[2];
          const MFloat VW_y[nDim] = {u_y - U_y * normalVec[0], v_y - U_y * normalVec[1], w_y - U_y * normalVec[2]};

          //
          const MFloat auxconst1 = gammaOverGammaMinusOne * p_y / pow(rho_y, gamma);
          const MFloat auxconst2 = 0; // only important in RANS
          const MFloat auxconst3 = -(gammaOverGammaMinusOne * p_y / rho_y + 0.5 * POW2(U_y));
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
          }

          if(testcounter > 10) {
            cout << "Loop in BC6002 took " << testcounter << " steps!"
                 << "cell: " << g[i] << "|" << t1 << " res=" << res << " (x|y)=" << m_cells->coordinates[0][cellIdG]
                 << "|" << m_cells->coordinates[1][cellIdG] << endl;
            mTerm(1, "testcounter > 10 in BC6002 loop");
          }

          // Compute remaining variables
          const MFloat temp = por_y * rho_y / (por_x * rho_x);
          const MFloat p_x = p_y * pow(rho_x / rho_y, gamma);
          const MFloat U_x = temp * U_y;
          const MFloat u_x = U_x * normalVec[0] + VW_y[0] /*=VW_x*/;
          const MFloat v_x = U_x * normalVec[1] + VW_y[1] /*=VW_x*/;
          const MFloat w_x = U_x * normalVec[2] + VW_y[2] /*=VW_x*/;

          // Assign solution to pvariables
          m_cells->pvariables[PV->RHO][cellIdG] = rho_x;
          m_cells->pvariables[PV->P][cellIdG] = p_x;
          m_cells->pvariables[PV->U][cellIdG] = u_x;
          m_cells->pvariables[PV->V][cellIdG] = v_x;
          m_cells->pvariables[PV->W][cellIdG] = w_x;
        }
      }
    }
  } else {
    mTerm(1, "This is not tested in 3D yet!");
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

    // TODO_SS labels:FV,totest Check if we maybe need to add +2 and -2 respectively in the tangential direction,
    // because those cells may
    //      be needed in viscousFluxCorrection
    for(MInt k = start[2], kk = start2[2]; k < end[2]; k++, kk++) {
      for(MInt j = start[1], jj = start2[1]; j < end[1]; j++, jj++) {
        for(MInt i = start[0], ii = start2[0]; i < end[0]; i++, ii++) {
          const MInt cellIdG = cellIndex(i, j, k);
          const MInt cellIdA1 = cellIndex(ii, jj, kk);
          const MFloat normalVec[nDim] = {m_cells->fq[FQ->NORMAL[0]][cellIdG], m_cells->fq[FQ->NORMAL[1]][cellIdG],
                                          m_cells->fq[FQ->NORMAL[2]][cellIdG]};
          const MFloat por_x = m_cells->fq[FQ->POROSITY][cellIdA1];

          const MFloat por_y = m_cells->fq[FQ->POROSITY][cellIdG];
          const MFloat rho_y = m_cells->pvariables[PV->RHO][cellIdG];
          const MFloat p_y = m_cells->pvariables[PV->P][cellIdG];
          const MFloat u_y = m_cells->pvariables[PV->U][cellIdG];
          const MFloat v_y = m_cells->pvariables[PV->V][cellIdG];
          const MFloat w_y = m_cells->pvariables[PV->W][cellIdG];

          // Compute normal and tangential velocity components
          const MFloat U_y = u_y * normalVec[0] + v_y * normalVec[1] + w_y * normalVec[2];
          const MFloat VW_y[nDim] = {u_y - U_y * normalVec[0], v_y - U_y * normalVec[1], w_y - U_y * normalVec[2]};

          //
          const MFloat auxconst1 = gammaOverGammaMinusOne * p_y / pow(rho_y, gamma);
          const MFloat auxconst2 = 0; // only important in RANS
          const MFloat auxconst3 = -(gammaOverGammaMinusOne * p_y / rho_y + 0.5 * POW2(U_y));
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
          const MFloat u_x = U_x * normalVec[0] + VW_y[0] /*=VW_x*/;
          const MFloat v_x = U_x * normalVec[1] + VW_y[1] /*=VW_x*/;
          const MFloat w_x = U_x * normalVec[2] + VW_y[2] /*=VW_x*/;

          // Assign solution to pvariables
          m_cells->pvariables[PV->RHO][cellIdG] = rho_x;
          m_cells->pvariables[PV->P][cellIdG] = p_x;
          m_cells->pvariables[PV->U][cellIdG] = u_x;
          m_cells->pvariables[PV->V][cellIdG] = v_x;
          m_cells->pvariables[PV->W][cellIdG] = w_x;
        }
      }
    }
  }
}


template <MBool isRans>
void StructuredBndryCnd3D<isRans>::bc2402(MInt bcId) {
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
  MInt cellId = 0;
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

      for(MInt k = startface[2]; k < endface[2]; k++) {
        for(MInt j = startface[1]; j < endface[1]; j++) {
          MInt ii = end[0] - 1;
          MInt cellIdP1 = cellIndex(ii + 1, j, k);
          MInt cellIdP2 = cellIndex(ii + 2, j, k);

          surface = sqrt(POW2(m_cells->surfaceMetrics[0][cellIdP1]) + POW2(m_cells->surfaceMetrics[1][cellIdP1])
                         + POW2(m_cells->surfaceMetrics[2][cellIdP1]));
          localPin[0] += surface * m_cells->pvariables[PV->P][cellIdP1];
          localPin[1] += surface * m_cells->pvariables[PV->P][cellIdP2];
          localTin[0] += surface * temperature(cellIdP1);
          localTin[1] += surface * temperature(cellIdP2);
          localVel += surface * (m_cells->pvariables[PV->U][cellIdP1]);
          localMassFlux += surface * (m_cells->pvariables[PV->U][cellIdP1] * m_cells->pvariables[PV->RHO][cellIdP1]);
        }
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

      for(MInt k = start[2]; k < end[2]; k++) {
        for(MInt j = start[1]; j < end[1]; j++) {
          MInt i = 1;
          cellId = cellIndex(i, j, k);

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
            m_cells->pvariables[var][cellIndex(start[0], j, k)] =
                2.0 * m_cells->pvariables[var][cellIndex(start[0] + 1, j, k)]
                - m_cells->pvariables[var][cellIndex(start[0] + 2, j, k)];
          }
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
      for(MInt k = startface[2]; k < endface[2]; k++) {
        for(MInt j = startface[1]; j < endface[1]; j++) {
          MInt ii = start[0];
          MInt cellIdM1 = cellIndex(ii - 1, j, k);
          MInt cellIdM2 = cellIndex(ii - 2, j, k);
          surface = sqrt(POW2(m_cells->surfaceMetrics[0][cellIdM1]) + POW2(m_cells->surfaceMetrics[1][cellIdM1])
                         + POW2(m_cells->surfaceMetrics[2][cellIdM1]));
          localPout[0] += surface * m_cells->pvariables[PV->P][cellIdM2];
          localPout[1] += surface * m_cells->pvariables[PV->P][cellIdM1];
          localTout[0] += surface * temperature(cellIdM2);
          localTout[1] += surface * temperature(cellIdM1);
        }
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

      for(MInt k = start[2]; k < end[2]; k++) {
        for(MInt j = start[1]; j < end[1]; j++) {
          MInt i = start[0];
          cellId = cellIndex(i, j, k);

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
            m_cells->pvariables[var][cellIndex(start[0] + 1, j, k)] =
                2.0 * m_cells->pvariables[var][cellIndex(start[0], j, k)]
                - m_cells->pvariables[var][cellIndex(start[0] - 1, j, k)];
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
inline MFloat StructuredBndryCnd3D<isRans>::pressure(MInt cellId) {
  return m_cells->pvariables[PV->P][cellId];
}

template <MBool isRans>
inline MFloat StructuredBndryCnd3D<isRans>::temperature(MInt cellId) {
  const MFloat gamma = m_solver->m_gamma;
  MFloat t = gamma * m_cells->pvariables[PV->P][cellId] / m_cells->pvariables[PV->RHO][cellId];
  return t;
}

template <MBool isRans>
inline MFloat StructuredBndryCnd3D<isRans>::pressure(MInt i, MInt j, MInt k) {
  MInt cellId = cellIndex(i, j, k);
  return pressure(cellId);
}

template <MBool isRans>
inline MInt StructuredBndryCnd3D<isRans>::pointIndex(MInt i, MInt j, MInt k) {
  return i + (j + k * m_nPoints[1]) * m_nPoints[2];
}

template <MBool isRans>
inline MInt StructuredBndryCnd3D<isRans>::getPointIdFromCell(MInt i, MInt j, MInt k) {
  return i + (j + k * m_nPoints[1]) * m_nPoints[2];
}

template <MBool isRans>
inline MInt StructuredBndryCnd3D<isRans>::getPointIdfromPoint(MInt origin, MInt incI, MInt incJ, MInt incK) {
  return origin + incI + incJ * m_nPoints[2] + incK * m_nPoints[2] * m_nPoints[1];
}

/**
 * \brief New function to compute the skin friction and pressure coefficient
 *        and the part for the force coefficients
 *        works in all directions and offers line averaging in one direction
 * \authors Pascal Meysonnat, Marian Albers
 */
template <MBool isRans>
template <MBool computePower>
void StructuredBndryCnd3D<isRans>::computeFrictionPressureCoef_(const MBool auxDataWindows) {
  // cf=mue*du/dy/(rho*u_8*u_8*0.5)
  // cp=2*(p-p_8)/(rho8*u_8*u_8*0.5)

  const MInt noForceCoefs = m_solver->m_noForceDataFields;

  // References for convenience
  auto& m_forceCoef = m_solver->m_forceCoef;

  const MFloat fre0 = F1 / (m_solver->m_Re0);
  const MFloat UT = m_solver->m_Ma * sqrt(PV->TInfinity);
  const MFloat fstagnationPressure = F1 / (CV->rhoInfinity * F1B2 * POW2(UT));
  const MFloat fstagnationEnergy = F1 / (CV->rhoInfinity * F1B2 * POW3(UT));

  // Only those for which auxData is requested in the property file; not those which are required
  // because of k-epsilon rans model
  const MInt noWalls = m_solver->m_windowInfo->m_auxDataWindowIds.size();

  // reset values to zero for the calculation
  for(MInt i = 0; i < (MInt)(noForceCoefs * noWalls); ++i) {
    m_forceCoef[i] = F0;
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

    const MInt mapOffsetCf = m_cells->cfOffsets[map];
    const MInt mapOffsetCp = m_cells->cpOffsets[map];
    const MInt mapOffsetPower = m_cells->powerOffsets[map];
    MInt* start = m_auxDataMap[map]->start1;
    MInt* end = m_auxDataMap[map]->end1;

    MInt n = 0;
    MInt normal = 0;
    MFloat area = F0;
    MFloat cf[6] = {F0, F0, F0, F0, F0, F0};
    MFloat cp[3] = {F0, F0, F0};
    MFloat powerp[3] = {F0, F0, F0};
    MFloat powerv[3] = {F0, F0, F0};

    const MInt pCoordDir[3][9] = {
        {0, 1, 0, 0, 0, 1, 0, 1, 1}, {1, 0, 0, 0, 0, 1, 1, 0, 1}, {1, 0, 0, 0, 1, 0, 1, 1, 0}};

    // indices
    MInt n10[3] = {0, 0, 0};
    MInt n01[3] = {0, 0, 0};
    MInt n1m1[3] = {0, 0, 0};
    MInt firstTangential = -1, secondTangential = -1;
    MInt i = 0, k = 0, j = 0;
    MInt jj = 0, kk = 0;
    MInt *reali = nullptr, *realj = nullptr, *realk = nullptr;
    MInt *realjj = nullptr, *realkk = nullptr;
    MInt sizeJ = 0, sizeK = 0;

    // determine the wall
    if((m_auxDataMap[map]->face) % 2 == 0) {
      n = -1; // bottom wall
      normal = m_auxDataMap[map]->face / 2;
      i = start[normal];
    } else {
      n = 1; // top wall
      normal = (m_auxDataMap[map]->face - 1) / 2;
      i = end[normal] - 1;
    }
    // determine the normals
    n1m1[normal] = -1 * n;
    n10[normal] = (MInt)(0.5 + (0.5 * (MFloat)n));
    n01[normal] = (MInt)(-0.5 + (0.5 * (MFloat)n));

    switch(normal) {
      case 0: {
        if(m_solver->m_forceAveragingDir == 2) {
          firstTangential = 2;
          secondTangential = 1;
          // pointers to the correct direction counter
          reali = &i;
          realj = &k;
          realk = &j;
          realjj = &kk;
          realkk = &jj;
        } else if(m_solver->m_forceAveragingDir == 1) {
          firstTangential = 1;
          secondTangential = 2;

          reali = &i;
          realj = &j;
          realk = &k;
          realjj = &jj;
          realkk = &kk;
        } else {
          mTerm(1, AT_, "Cant average in normal direction");
        }
        // Saves normal directions for three other points of surface

        // determine the sizes
        sizeJ = end[1] - start[1];
        sizeK = end[2] - start[2];
        break;
      }
      case 1: {
        if(m_solver->m_forceAveragingDir == 2) {
          firstTangential = 2;
          secondTangential = 0;
          reali = &k;
          realj = &i;
          realk = &j;
          realjj = &kk;
          realkk = &jj;
        } else if(m_solver->m_forceAveragingDir == 0) {
          firstTangential = 0;
          secondTangential = 2;
          reali = &j;
          realj = &i;
          realk = &k;
          realjj = &jj;
          realkk = &kk;
        } else {
          mTerm(1, AT_, "Cant average in normal direction");
        }

        sizeJ = end[0] - start[0];
        sizeK = end[2] - start[2];


        break;
      }
      case 2: {
        if(m_solver->m_forceAveragingDir == 1) {
          firstTangential = 1;
          secondTangential = 0;

          reali = &k;
          realj = &j;
          realk = &i;
          realjj = &kk;
          realkk = &jj;
        } else if(m_solver->m_forceAveragingDir == 0) {
          firstTangential = 0;
          secondTangential = 1;

          reali = &j;
          realj = &k;
          realk = &i;
          realjj = &jj;
          realkk = &kk;
        } else {
          mTerm(1, AT_, "Cant average in normal direction");
        }


        sizeJ = end[0] - start[0];
        sizeK = end[1] - start[1];
        break;
      }
      default: {
        mTerm(1, AT_, "Normal direction not implemented");
      }
    }

    // main loop over the two tangential directions of the wall
    for(k = start[secondTangential]; k < end[secondTangential]; k++) {
      jj = 0;

      // this loop always goes into the line averaging direction
      for(j = start[firstTangential]; j < end[firstTangential]; j++) {
        // get id of boundary cell and cell above that one for extrapolation
        const MInt cellId = cellIndex(*reali, *realj, *realk);
        const MInt cellIdP1 = cellIndex(*reali + n1m1[0], *realj + n1m1[1], *realk + n1m1[2]);

        // get point id and ids of three other surface corners
        const MInt pIJK = getPointIdFromCell(*reali + n10[0], *realj + n10[1], *realk + n10[2]);
        const MInt pIJPK = getPointIdfromPoint(pIJK, pCoordDir[normal][0], pCoordDir[normal][1], pCoordDir[normal][2]);
        const MInt pIJKP = getPointIdfromPoint(pIJK, pCoordDir[normal][3], pCoordDir[normal][4], pCoordDir[normal][5]);
        const MInt pIJPKP = getPointIdfromPoint(pIJK, pCoordDir[normal][6], pCoordDir[normal][7], pCoordDir[normal][8]);

        // first get the position of the wall (reference!!)
        // x is always used as coordinate in normal direction
        const MFloat xRef = F1B4
                            * (m_grid->m_coordinates[0][pIJK] + m_grid->m_coordinates[0][pIJPK]
                               + m_grid->m_coordinates[0][pIJKP] + m_grid->m_coordinates[0][pIJPKP]);
        const MFloat yRef = F1B4
                            * (m_grid->m_coordinates[1][pIJK] + m_grid->m_coordinates[1][pIJPK]
                               + m_grid->m_coordinates[1][pIJKP] + m_grid->m_coordinates[1][pIJPKP]);
        const MFloat zRef = F1B4
                            * (m_grid->m_coordinates[2][pIJK] + m_grid->m_coordinates[2][pIJPK]
                               + m_grid->m_coordinates[2][pIJKP] + m_grid->m_coordinates[2][pIJPKP]);

        MFloat uWall = F0, vWall = F0, wWall = F0;
        if(m_solver->m_movingGrid) {
          uWall = F1B4
                  * (m_grid->m_velocity[0][pIJK] + m_grid->m_velocity[0][pIJPK] + m_grid->m_velocity[0][pIJKP]
                     + m_grid->m_velocity[0][pIJPKP]);
          vWall = F1B4
                  * (m_grid->m_velocity[1][pIJK] + m_grid->m_velocity[1][pIJPK] + m_grid->m_velocity[1][pIJKP]
                     + m_grid->m_velocity[1][pIJPKP]);
          wWall = F1B4
                  * (m_grid->m_velocity[2][pIJK] + m_grid->m_velocity[2][pIJPK] + m_grid->m_velocity[2][pIJKP]
                     + m_grid->m_velocity[2][pIJPKP]);
        }

        // get the distcance
        const MFloat dx2 = sqrt(POW2(m_cells->coordinates[0][cellIdP1] - m_cells->coordinates[0][cellId])
                                + POW2(m_cells->coordinates[1][cellIdP1] - m_cells->coordinates[1][cellId])
                                + POW2(m_cells->coordinates[2][cellIdP1] - m_cells->coordinates[2][cellId]));
        const MFloat dx1 =
            sqrt(POW2(m_cells->coordinates[0][cellId] - xRef) + POW2(m_cells->coordinates[1][cellId] - yRef)
                 + POW2(m_cells->coordinates[2][cellId] - zRef));
        // pressures
        const MFloat p1 = m_cells->pvariables[PV->P][cellId];
        const MFloat p2 = m_cells->pvariables[PV->P][cellIdP1];
        // extrapolation to the face
        const MFloat pW = ((p1 - p2) / dx2) * dx1 + p1;
        const MFloat cpn = (pW - PV->PInfinity) * fstagnationPressure;

        // save pressure to map
        m_cells->cp[mapOffsetCp + (*realjj) + (*realkk) * sizeJ] = cpn;


        // compute the skin-friction coefficient
        //-->
        // compute orthogonal distance surface-plane to cell center
        MFloat supportVec[3] = {F0, F0, F0};
        MFloat firstVec[3] = {F0, F0, F0};
        MFloat secondVec[3] = {F0, F0, F0};
        MFloat normalVec[3] = {F0, F0, F0};
        MFloat cellVec[3] = {F0, F0, F0};

        for(MInt dim = 0; dim < nDim; dim++) {
          supportVec[dim] = m_grid->m_coordinates[dim][pIJK];
          firstVec[dim] = m_grid->m_coordinates[dim][pIJPK] - m_grid->m_coordinates[dim][pIJK];
          secondVec[dim] = m_grid->m_coordinates[dim][pIJKP] - m_grid->m_coordinates[dim][pIJK];
          cellVec[dim] = m_cells->coordinates[dim][cellId];
        }

        crossProduct(normalVec, firstVec, secondVec);
        const MFloat normalLength = sqrt(POW2(normalVec[0]) + POW2(normalVec[1]) + POW2(normalVec[2]));

        const MFloat nn[3] = {normalVec[0] / normalLength, normalVec[1] / normalLength, normalVec[2] / normalLength};
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
        const MFloat w1 = m_cells->pvariables[PV->W][cellId];

        // compute gradient
        MFloat dudeta = (u1 - uWall) / orthDist;
        MFloat dvdeta = (v1 - vWall) / orthDist;
        MFloat dwdeta = (w1 - wWall) / orthDist;

        // using taylor expansion a second-order approximation
        // of du/dn can be achieved at the wall
        if(m_solver->m_forceSecondOrder) {
          const MFloat u2 = m_cells->pvariables[PV->U][cellIdP1];
          const MFloat v2 = m_cells->pvariables[PV->V][cellIdP1];
          const MFloat w2 = m_cells->pvariables[PV->W][cellIdP1];
          const MFloat dx3 = dx1 + dx2;
          dudeta =
              (u1 * dx3 * dx3 + (dx1 * dx1 - dx3 * dx3) * uWall - u2 * dx1 * dx1) / (dx1 * dx3 * dx3 - dx1 * dx1 * dx3);
          dvdeta =
              (v1 * dx3 * dx3 + (dx1 * dx1 - dx3 * dx3) * vWall - v2 * dx1 * dx1) / (dx1 * dx3 * dx3 - dx1 * dx1 * dx3);
          dwdeta =
              (w1 * dx3 * dx3 + (dx1 * dx1 - dx3 * dx3) * wWall - w2 * dx1 * dx1) / (dx1 * dx3 * dx3 - dx1 * dx1 * dx3);
        }

        // the compressible contribution to the stress tensor
        const MFloat comp = F1B3 * (dudeta * nn[0] + dvdeta * nn[1] + dwdeta * nn[2]);

        // friction coefficient in computational space
        const MFloat taux = mue * dudeta * fre0;
        const MFloat tauy = mue * dvdeta * fre0;
        const MFloat tauz = mue * dwdeta * fre0;

        const MFloat cfx = taux * fstagnationPressure;
        const MFloat cfy = tauy * fstagnationPressure;
        const MFloat cfz = tauz * fstagnationPressure;

        const MFloat taux_c = mue * comp * nn[0] * fre0;
        const MFloat tauy_c = mue * comp * nn[1] * fre0;
        const MFloat tauz_c = mue * comp * nn[2] * fre0;

        const MFloat cfx_c = taux_c * fstagnationPressure;
        const MFloat cfy_c = tauy_c * fstagnationPressure;
        const MFloat cfz_c = tauz_c * fstagnationPressure;

        // save to map
        // skin-friction
        m_cells->cf[mapOffsetCf + 0 * sizeJ * sizeK + (*realjj) + (*realkk) * sizeJ] = cfx + cfx_c;
        m_cells->cf[mapOffsetCf + 1 * sizeJ * sizeK + (*realjj) + (*realkk) * sizeJ] = cfy + cfy_c;
        m_cells->cf[mapOffsetCf + 2 * sizeJ * sizeK + (*realjj) + (*realkk) * sizeJ] = cfz + cfz_c;
        //<-- finished computation of skin-friction


        MFloat considerValue = F1;
        if(m_solver->m_auxDataCoordinateLimits) {
          if(m_solver->m_auxDataLimits[0] <= xRef && xRef <= m_solver->m_auxDataLimits[1]
             && m_solver->m_auxDataLimits[2] <= zRef && zRef <= m_solver->m_auxDataLimits[3]) {
            considerValue = F1;
          } else {
            considerValue = F0;
          }
        }

        // sum up all lengths and all cps with their surface width contribution
        const MFloat dxidx =
            m_cells->surfaceMetrics[nDim * normal + 0][cellIndex(*reali + n01[0], *realj + n01[1], *realk + n01[2])];
        const MFloat dxidy =
            m_cells->surfaceMetrics[nDim * normal + 1][cellIndex(*reali + n01[0], *realj + n01[1], *realk + n01[2])];
        const MFloat dxidz =
            m_cells->surfaceMetrics[nDim * normal + 2][cellIndex(*reali + n01[0], *realj + n01[1], *realk + n01[2])];
        const MFloat dA = sqrt(dxidx * dxidx + dxidy * dxidy + dxidz * dxidz);

        cp[0] += (-1.0) * cpn * dxidx * considerValue;
        cp[1] += (-1.0) * cpn * dxidy * considerValue;
        cp[2] += (-1.0) * cpn * dxidz * considerValue;
        // cf
        cf[0] += cfx * dA * considerValue;
        cf[1] += cfy * dA * considerValue;
        cf[2] += cfz * dA * considerValue;

        cf[3] += cfx_c * dA * considerValue;
        cf[4] += cfy_c * dA * considerValue;
        cf[5] += cfz_c * dA * considerValue;
        // area
        area += dA * considerValue;

        if(computePower) {
          // compute the power (p*v)
          const MFloat dp = (pW - PV->PInfinity);

          // this is power without area contribution, i.e.: P/A
          const MFloat P_px = (-1.0) * dp * uWall * fstagnationEnergy;
          const MFloat P_py = (-1.0) * dp * vWall * fstagnationEnergy;
          const MFloat P_pz = (-1.0) * dp * wWall * fstagnationEnergy;

          // save power due to pressure to map
          m_cells->powerPres[mapOffsetPower + 0 * sizeJ * sizeK + (*realjj) + (*realkk) * sizeJ] = P_px;
          m_cells->powerPres[mapOffsetPower + 1 * sizeJ * sizeK + (*realjj) + (*realkk) * sizeJ] = P_py;
          m_cells->powerPres[mapOffsetPower + 2 * sizeJ * sizeK + (*realjj) + (*realkk) * sizeJ] = P_pz;

          // Power due to viscous forces in the computational space ( tau_w*n*u)
          const MFloat P_cfx = (taux + taux_c) * uWall * fstagnationEnergy;
          const MFloat P_cfy = (tauy + tauy_c) * vWall * fstagnationEnergy;
          const MFloat P_cfz = (tauz + tauz_c) * wWall * fstagnationEnergy;

          // save power to map
          m_cells->powerVisc[mapOffsetPower + 0 * sizeJ * sizeK + (*realjj) + (*realkk) * sizeJ] = P_cfx;
          m_cells->powerVisc[mapOffsetPower + 1 * sizeJ * sizeK + (*realjj) + (*realkk) * sizeJ] = P_cfy;
          m_cells->powerVisc[mapOffsetPower + 2 * sizeJ * sizeK + (*realjj) + (*realkk) * sizeJ] = P_cfz;

          // compute pressurepower*area -> P = (P/A) * A
          powerp[0] += P_px * dxidx * considerValue;
          powerp[1] += P_py * dxidy * considerValue;
          powerp[2] += P_pz * dxidz * considerValue;
          // viscous power
          powerv[0] += P_cfx * (sqrt(dxidx * dxidx + dxidy * dxidy + dxidz * dxidz)) * considerValue;
          powerv[1] += P_cfy * (sqrt(dxidx * dxidx + dxidy * dxidy + dxidz * dxidz)) * considerValue;
          powerv[2] += P_cfz * (sqrt(dxidx * dxidx + dxidy * dxidy + dxidz * dxidz)) * considerValue;
        }

        if(m_solver->m_detailAuxData) {
          // save surface contributions
          m_cells->cf[mapOffsetCf + 3 * sizeJ * sizeK + (*realjj) + (*realkk) * sizeJ] = dxidx;
          m_cells->cf[mapOffsetCf + 4 * sizeJ * sizeK + (*realjj) + (*realkk) * sizeJ] = dxidy;
          m_cells->cf[mapOffsetCf + 5 * sizeJ * sizeK + (*realjj) + (*realkk) * sizeJ] = dxidz;

          // save surface coordinates
          m_cells->cf[mapOffsetCf + 6 * sizeJ * sizeK + (*realjj) + (*realkk) * sizeJ] = xRef;
          m_cells->cf[mapOffsetCf + 7 * sizeJ * sizeK + (*realjj) + (*realkk) * sizeJ] = yRef;
          m_cells->cf[mapOffsetCf + 8 * sizeJ * sizeK + (*realjj) + (*realkk) * sizeJ] = zRef;
        }

        ++jj;
      }

      ++kk;
    }

    // now add to force coefficients
    MInt count = 0;
    for(auto it = m_solver->m_windowInfo->m_auxDataWindowIds.cbegin();
        it != m_solver->m_windowInfo->m_auxDataWindowIds.cend();
        ++it) {
      if(m_auxDataMap[map]->Id2 == it->second) {
        m_forceCoef[count * noForceCoefs + 0] += cf[0];
        m_forceCoef[count * noForceCoefs + 1] += cf[1];
        m_forceCoef[count * noForceCoefs + 2] += cf[2];
        m_forceCoef[count * noForceCoefs + 3] += cf[3];
        m_forceCoef[count * noForceCoefs + 4] += cf[4];
        m_forceCoef[count * noForceCoefs + 5] += cf[5];
        m_forceCoef[count * noForceCoefs + 6] += cp[0];
        m_forceCoef[count * noForceCoefs + 7] += cp[1];
        m_forceCoef[count * noForceCoefs + 8] += cp[2];
        m_forceCoef[count * noForceCoefs + 9] += area;

        if(computePower) {
          m_forceCoef[count * noForceCoefs + 10] += powerv[0];
          m_forceCoef[count * noForceCoefs + 11] += powerv[1];
          m_forceCoef[count * noForceCoefs + 12] += powerv[2];
          m_forceCoef[count * noForceCoefs + 13] += powerp[0];
          m_forceCoef[count * noForceCoefs + 14] += powerp[1];
          m_forceCoef[count * noForceCoefs + 15] += powerp[2];
        }
      }

      count++;
    }
  }
}
template void StructuredBndryCnd3D<true>::computeFrictionPressureCoef_<true>(const MBool);
template void StructuredBndryCnd3D<false>::computeFrictionPressureCoef_<true>(const MBool);
template void StructuredBndryCnd3D<true>::computeFrictionPressureCoef_<false>(const MBool);
template void StructuredBndryCnd3D<false>::computeFrictionPressureCoef_<false>(const MBool);


/**
 * Function to compute the momentum coefficient (LE): cm_le
 */
template <MBool isRans>
void StructuredBndryCnd3D<isRans>::computeMomentCoef() {}

template <MBool isRans>
inline void StructuredBndryCnd3D<isRans>::crossProduct(MFloat* result, MFloat* vec1, MFloat* vec2) {
  result[0] = vec1[1] * vec2[2] - vec1[2] * vec2[1];
  result[1] = vec1[2] * vec2[0] - vec1[0] * vec2[2];
  result[2] = vec1[0] * vec2[1] - vec1[1] * vec2[0];
}

template <MBool isRans>
inline MInt StructuredBndryCnd3D<isRans>::getReverseCellId(MInt i, MInt j, MInt k, MInt face) {
  const MInt i_new = m_reverseCellIdGC[face * 3 + 0] * (m_noGhostLayers - 1) + m_reverseCellIdDim[face * 3 + 0] * i;
  const MInt j_new = m_reverseCellIdGC[face * 3 + 1] * (m_noGhostLayers - 1) + m_reverseCellIdDim[face * 3 + 1] * j;
  const MInt k_new = m_reverseCellIdGC[face * 3 + 2] * (m_noGhostLayers - 1) + m_reverseCellIdDim[face * 3 + 2] * k;
  return cellIndex(i_new, j_new, k_new);
}

template <MBool isRans>
inline MInt StructuredBndryCnd3D<isRans>::getExtrNghbrId(MInt cellId, MInt face) {
  static const MInt cellShift[3] = {1, m_nCells[2], m_nCells[1] * m_nCells[2]};
  const MInt n = 1 - 2 * (face % 2);
  const MInt dim = face / 2;
  const MInt nghbrId = cellId + n * cellShift[dim];
  return nghbrId;
}

template <MBool isRans>
inline pair<MInt, MInt> StructuredBndryCnd3D<isRans>::getMirrorCellIdPair(MInt i, MInt j, MInt k, MInt face) {
  const MInt cellShift[3] = {1, m_nCells[2], m_nCells[1] * m_nCells[2]};
  const MInt gcPos[6] = {m_noGhostLayers, m_nCells[2] - m_noGhostLayers - 1,
                         m_noGhostLayers, m_nCells[1] - m_noGhostLayers - 1,
                         m_noGhostLayers, m_nCells[0] - m_noGhostLayers - 1};
  const MInt ijk_new[3] = {
      m_reverseCellIdGC[face * 3 + 0] * (m_noGhostLayers - 1) + m_reverseCellIdDim[face * 3 + 0] * i,
      m_reverseCellIdGC[face * 3 + 1] * (m_noGhostLayers - 1) + m_reverseCellIdDim[face * 3 + 1] * j,
      m_reverseCellIdGC[face * 3 + 2] * (m_noGhostLayers - 1) + m_reverseCellIdDim[face * 3 + 2] * k};
  const MInt dim = face / 2;
  const MInt n = 1 - 2 * (face % 2);
  const MInt mirror = ((gcPos[face] - ijk_new[dim]) * 2 - n);
  const MInt cellId = cellIndex(ijk_new[0], ijk_new[1], ijk_new[2]);
  return make_pair(cellId, cellId + mirror * cellShift[dim]);
}

// instantanisation of the class
template class StructuredBndryCnd3D<true>;
template class StructuredBndryCnd3D<false>;
