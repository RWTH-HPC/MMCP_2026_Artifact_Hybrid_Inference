// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "fcbndrycnd.h"
#include "fcsolver.h"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <set>
#include "COMM/mpioverride.h"
#include "GEOM/geometrycontext.h"
#include "IO/context.h"
#include "IO/infoout.h"
#include "fcgridbndrycell.h"
#include "property.h"

using namespace std;

template <MInt nDim>
FcBndryCnd<nDim>::FcBndryCnd(FcSolver<nDim>* solver)
  : m_solverId(solver->m_solverId), m_noInternalCells(solver->grid().noInternalCells()) {
  TRACE();

  m_solver = solver;

  // the number of all segments
  m_noSegments = m_solver->m_geometry->geometryContext().getNoSegments();

  /*! \page propertyPage1
    \section kFactor
    <code>MInt FcBndryCnd::m_kFactor</code>\n
    default = <code>1e+12</code>\n\n
    Factor to set for fixation bndry conditions\n
    Keywords: <i>FINITE_CELL</i>
  */
  m_kFactor = 1e+12;
  m_kFactor = Context::getSolverProperty<MFloat>("kFactor", m_solverId, AT_, &m_kFactor);

  /*! \page propertyPage1
    \section subCellLayerDepth
    <code>MInt FcBndryCnd::m_subCellLayerDepth</code>\n
    default = <code>0</code>\n\n
    Subcell layer depth for sub cell integration for each bndry\n
    Keywords: <i>FINITE_CELL</i>
  */
  if(Context::propertyExists("subCellLayerDepth", m_solverId)) {
    mAlloc(m_subCellLayerDepth, m_noSegments, "m_subCellLayerDepth", 0, AT_);
    for(MInt i = 0; i < m_noSegments; i++) {
      m_subCellLayerDepth[i] = Context::getSolverProperty<MInt>("subCellLayerDepth", m_solverId, AT_, i);
    }
  }

  // the number of cells that carry a certain boundary condition
  mAlloc(m_noBndryCellsPerSegment, m_noSegments, "m_noBndryCellsPerSegment", 0, AT_);

  /*! \page propertyPage1
    \section multiBC
    <code>MBool FcBndryCnd::multiBCTreatment</code>\n
    default = <code>true</code>\n\n
    This property defines if cells can have multiple BCs
    <ul>
    <li><code>true</code> all BCs are applied to the cell</li>
    <li><code>I-P-W</code> only the first/last BC is applied to the cell</li> //TODO: Check if first or last
    </ul>\n
    Keywords: <i>FINITE_CELL</i>
  */
  m_multiBC = true;
  m_multiBC = Context::getSolverProperty<MBool>("multiBC", m_solverId, AT_, &m_multiBC);
  if(!m_multiBC) mTerm(1, AT_, "m_multiBC is set to false, bc functions are not adapted to this input.");

  /*! \page propertyPage1
    \section bndryNormalMethod
    <code>MString FcBndryCnd::bndryNormalMethod</code>\n
    default = <code>calcNormal</code>\n\n
    This property defines the way segment normals are calculated.
    <ul>
    <li><code>read</code> reads the normals from the property file (bndryNormalVectors)</li>
    TODO:<li><code>calcNormal</code> calculates the normals based on triangle information and averages them</li>
    TODO:<li><code>fromSTL</code> reads the normals from STL and averages them</li>
    </ul>\n
    Keywords: <i>FINITE_CELL</i>
  */
  m_bndryNormalMethod = "read";
  m_bndryNormalMethod = Context::getSolverProperty<MString>("bndryNormalMethod", m_solverId, AT_, &m_bndryNormalMethod);
  if(m_bndryNormalMethod != "read")
    mTerm(1, AT_, "m_bndryNormalMethod is not \"read\", bc functions are not adapted to this input.");

  /// Create and sort boundary cells
  createBoundaryCells();
  sortBoundaryCells();

  setBndryCndHandler();

  initBndryCnds();

  findCellsRequireSubCellIntegration();

  calculateCutPoints();

  for(MInt i = 0; i < (MInt)m_bndryCndSegIds.size(); i++) {
    if(m_solver->m_testRun) {
      writeOutSTLBoundaries(i, "preSimulation");
    }
    calcAvgFaceNormal(i);
  }
}

template <MInt nDim>
FcBndryCnd<nDim>::~FcBndryCnd() {
  TRACE();

  // Clean up allocated memory
  mDeallocate(m_noBndryCellsPerSegment);
}

/** \brief Initializes boundary cells
    \author Moritz Waldmann
    \date 29.05.2021

    Sets the direction of fixation and the direction of loads at the boundaries.

    That is, the property segFixationDirs_x is 0.0 or 1.0 for segment x in each direction.
    If it is 1.0 the displacement in that direction is set to 0.
    The same applies for the property segLoadDirs_x, which sets the direction in which a load
    is applied. The propterty segLoadValue_x contains the explicite load value in each
    Cartesian direction.
  */
template <MInt nDim>
void FcBndryCnd<nDim>::initBndryCnds() {
  TRACE();

  for(MInt i = 0; i < (MInt)(m_bndryCndSegIds.size()); i++) {
    MInt seg_id = m_bndryCndSegIds[i];
    MInt bc_id = m_bndryCndIds[i];

    ostringstream ostrSegId;
    ostrSegId << seg_id;
    MString strSegId = ostrSegId.str();

    MString name = "";
    // Read bc specific properties
    switch(bc_id) {
      case 0: {
        break;
      }
      case 8010: {
        name = "segFixationDirs_";
        name.append(strSegId);
        MFloat dir[nDim] = {F0};
        for(MInt d = 0; d < nDim; d++) {
          dir[d] = Context::getSolverProperty<MFloat>(name, m_solverId, AT_, d);
        }
        for(MInt j = m_bndryCndOffsets[i]; j < m_bndryCndOffsets[i + 1]; j++) {
          for(MInt d = 0; d < nDim; d++) {
            m_bndryCells[j].m_directionOfAction[d] = dir[d];
          }
        }
        break;
      }
      case 8011: {
        name = "segFixationDirs_";
        name.append(strSegId);
        MFloat dir[nDim] = {F0};
        for(MInt d = 0; d < nDim; d++) {
          dir[d] = Context::getSolverProperty<MFloat>(name, m_solverId, AT_, d);
        }
        for(MInt j = m_bndryCndOffsets[i]; j < m_bndryCndOffsets[i + 1]; j++) {
          for(MInt d = 0; d < nDim; d++) {
            m_bndryCells[j].m_directionOfAction[d] = dir[d];
          }
        }
        break;
      }
      case 8012: {
        name = "segFixationDirs_";
        name.append(strSegId);
        MFloat dir[nDim] = {F0};
        for(MInt d = 0; d < nDim; d++) {
          dir[d] = Context::getSolverProperty<MFloat>(name, m_solverId, AT_, d);
        }
        for(MInt j = m_bndryCndOffsets[i]; j < m_bndryCndOffsets[i + 1]; j++) {
          for(MInt d = 0; d < nDim; d++) {
            m_bndryCells[j].m_directionOfAction[d] = dir[d];
          }
        }
        break;
      }
      case 8020: {
        break;
      }
      case 8030: {
        name = "segLoadDirs_";
        name.append(strSegId);
        MFloat dir[nDim] = {F0};
        for(MInt d = 0; d < nDim; d++) {
          dir[d] = Context::getSolverProperty<MFloat>(name, m_solverId, AT_, d);
        }
        for(MInt j = m_bndryCndOffsets[i]; j < m_bndryCndOffsets[i + 1]; j++) {
          for(MInt d = 0; d < nDim; d++) {
            m_bndryCells[j].m_directionOfAction[d] = dir[d];
          }
        }
        name = "segLoadValue_";
        name.append(strSegId);
        MFloat load[nDim] = {F0};
        for(MInt d = 0; d < nDim; d++) {
          load[d] = Context::getSolverProperty<MFloat>(name, m_solverId, AT_, d);
        }
        for(MInt j = m_bndryCndOffsets[i]; j < m_bndryCndOffsets[i + 1]; j++) {
          for(MInt d = 0; d < nDim; d++) {
            m_bndryCells[j].m_load[d] = load[d];
          }
        }
        break;
      }
      case 8031: {
        name = "segLoadDirs_";
        name.append(strSegId);
        MFloat dir[nDim] = {F0};
        for(MInt d = 0; d < nDim; d++) {
          dir[d] = Context::getSolverProperty<MFloat>(name, m_solverId, AT_, d);
        }
        for(MInt j = m_bndryCndOffsets[i]; j < m_bndryCndOffsets[i + 1]; j++) {
          for(MInt d = 0; d < nDim; d++) {
            m_bndryCells[j].m_directionOfAction[d] = dir[d];
          }
        }
        name = "segLoadValue_";
        name.append(strSegId);
        MFloat load[nDim] = {F0};
        for(MInt d = 0; d < nDim; d++) {
          load[d] = Context::getSolverProperty<MFloat>(name, m_solverId, AT_, d);
        }
        for(MInt j = m_bndryCndOffsets[i]; j < m_bndryCndOffsets[i + 1]; j++) {
          for(MInt d = 0; d < nDim; d++) {
            m_bndryCells[j].m_load[d] = load[d];
          }
        }
        break;
      }
      case 8032: {
        name = "segLoadDirs_";
        name.append(strSegId);
        MFloat dir[nDim] = {F0};
        for(MInt d = 0; d < nDim; d++) {
          dir[d] = Context::getSolverProperty<MFloat>(name, m_solverId, AT_, d);
        }
        for(MInt j = m_bndryCndOffsets[i]; j < m_bndryCndOffsets[i + 1]; j++) {
          for(MInt d = 0; d < nDim; d++) {
            m_bndryCells[j].m_directionOfAction[d] = dir[d];
          }
        }
        name = "segLoadValue_";
        name.append(strSegId);
        MFloat load[nDim] = {F0};
        for(MInt d = 0; d < nDim; d++) {
          load[d] = Context::getSolverProperty<MFloat>(name, m_solverId, AT_, d);
        }
        for(MInt j = m_bndryCndOffsets[i]; j < m_bndryCndOffsets[i + 1]; j++) {
          for(MInt d = 0; d < nDim; d++) {
            m_bndryCells[j].m_load[d] = load[d];
          }
        }
        break;
      }
      case 8035: {
        name = "segLoadDirs_";
        name.append(strSegId);
        MFloat dir[nDim] = {F0};
        for(MInt d = 0; d < nDim; d++) {
          dir[d] = Context::getSolverProperty<MFloat>(name, m_solverId, AT_, d);
        }
        for(MInt j = m_bndryCndOffsets[i]; j < m_bndryCndOffsets[i + 1]; j++) {
          for(MInt d = 0; d < nDim; d++) {
            m_bndryCells[j].m_directionOfAction[d] = dir[d];
          }
        }
        name = "segLoadValue_";
        name.append(strSegId);
        MFloat load[nDim] = {F0};
        for(MInt d = 0; d < nDim; d++) {
          load[d] = Context::getSolverProperty<MFloat>(name, m_solverId, AT_, d);
        }
        for(MInt j = m_bndryCndOffsets[i]; j < m_bndryCndOffsets[i + 1]; j++) {
          for(MInt d = 0; d < nDim; d++) {
            m_bndryCells[j].m_load[d] = load[d];
          }
        }
        break;
      }
      default: {
        mTerm(1, AT_, "Unknown boundary condition id!");
      }
    }
  }
}

/** \brief Creates boundary cells according to the geometry information
    \author Andreas Lintermann
    \date 29.01.2010
    \adaptation 22.04.2021, Felix Wietbuescher

    This function runs over all cells and makes an intersection test with the triangles of nearby segments.
    All candidates carry a segment id and a boundary condition. Hence, a corner cell can have multiple properties.
    Depending on the allowance of multiple BCs, only one single (no multiple BCs allowed) boundary cell is added to
    a temporary vector. Anyhow, the boundary cell can carry multiple information on segment id and boundary condition.

    //In this case the wall boundary condition is chosen by setting the first element of both the boundary condition
   array
    //and the segment id to this wall condition.

    In the other case of multiple BCs, the boundary cell is added multiple
    times to the vector of temporary boundary cells, each time with its according different conditions. After that,
    the collector of boundary cells is initialized. The size is now known by the size of the temporary vector. Finally,
    the information in the vector is copied to the collector.
  */
template <MInt nDim>
void FcBndryCnd<nDim>::createBoundaryCells() {
  TRACE();

  MInt noMultiCells = 0; // Number of cells with multiple BCs

  MInt bndryCellsIt = 0, bndryCellsIt2 = 0;

  m_log << "  + Creating boundary cells..." << endl;
  m_log << "    - Multiple BC usage: " << m_multiBC << endl;

  for(MInt i = 0; i < m_solver->m_cells.size(); i++) {
    if(!m_solver->c_isLeafCell(i)) continue;

    std::vector<MInt> nodeList;
    getStlNodeList(i, nodeList);

    const MInt noNodes = (MInt)nodeList.size();
    if(noNodes > 0) {
      m_solver->a_isBndryCell(i) = true;

      // Create a new boundary cell for our found cuts
      m_bndryCells.emplace_back();

      bndryCellsIt = m_bndryCells.size() - 1;
      m_bndryCells[bndryCellsIt].m_cellId = i;

      // Fill the segmentIds of the boundary cell
      // and the associated bndryCndId
      for(MInt j = 0; j < noNodes; j++) {
        MBool already_in = false;
        MInt segId = m_solver->m_geometry->elements[nodeList[j]].m_segmentId;
        MInt bndId = m_solver->m_geometry->elements[nodeList[j]].m_bndCndId;
        for(MInt k = 0; k < (MInt)(m_bndryCells[bndryCellsIt].m_segmentId.size()); k++) {
          if(m_bndryCells[bndryCellsIt].m_segmentId[k] == segId) {
            already_in = true;
            break;
          }
        }
        if(!already_in) {
          m_bndryCells[bndryCellsIt].m_segmentId.push_back(segId);
          m_bndryCells[bndryCellsIt].m_bndryCndId.push_back(bndId);
        }
      }

      if(m_bndryCells[bndryCellsIt].m_segmentId.size() > 1) {
        noMultiCells++;
        if(m_multiBC) {
          // Add the bndryCell n more times
          //! according to the number of different bc's.
          // No sorting of segmentIds is applied.

          bndryCellsIt2 = bndryCellsIt;
          for(MInt j = 1; j < (MInt)(m_bndryCells[bndryCellsIt].m_segmentId.size()); j++) {
            // make sure that bndryCells are only added if there is another bc
            MInt bndId = m_bndryCells[bndryCellsIt].m_bndryCndId[j];
            MInt segId = m_bndryCells[bndryCellsIt].m_segmentId[j];

            bndryCellsIt2++;
            m_bndryCells.emplace_back(); // add new bndryCell
            m_bndryCells[bndryCellsIt2].m_cellId = i;

            m_bndryCells[bndryCellsIt2].m_segmentId.push_back(segId);
            m_bndryCells[bndryCellsIt2].m_bndryCndId.push_back(bndId);
          }
        }
      }
    }
  }
  m_log << "    - noBndCells:            " << m_bndryCells.size() << endl;
  m_log << "    - noMultiCells:          " << noMultiCells << endl << endl;
}

/** \brief This function sorts the boundary cells according to the BC id.
    \author Andreas Lintermann
    \date 27.01.2010

    The sorting of the boundary cells is necessary because all cells are
    stored in one single collector. The boundary condition functions get
    the position of a type of boundary cells in the collector (offset) and
    the number of boundary cells of equal type ( i.e. with the same
    boundary id) as parameter.
  */
template <MInt nDim>
void FcBndryCnd<nDim>::sortBoundaryCells() {
  TRACE();

  // Sort the boundary cells
  const MInt tmpCellSize = m_bndryCells.size();

  m_log << "  + Sorting boundary cells..." << endl;
  m_log << "    - no. boundary cells: " << tmpCellSize << endl;

  // using stable_sort to preserve the order by cellId within segmentIds
  stable_sort(m_bndryCells.begin(), m_bndryCells.end(), //
              [](auto a, auto b) { return a.m_segmentId[0] < b.m_segmentId[0]; });

  MInt tmpSegmentId = -1;  // holds the current id
  MInt tmpBndryCndId = -1; // holds the current boundary cells bndryCndId
  MInt counter = 0;        // Counts the boundary cells

  m_mapBndryCndSegId2Index.resize(m_noSegments, -1);

  for(MInt i = 0; i < tmpCellSize; i++) {
    if(tmpSegmentId != m_bndryCells[i].m_segmentId[0]) {
      // Since bndryCells has been sorted previously and new segmentid differs
      // from tmpSegmentId, we sort now for a new tmpSegmentId and co
      tmpSegmentId = m_bndryCells[i].m_segmentId[0];
      tmpBndryCndId = m_bndryCells[i].m_bndryCndId[0];

      m_bndryCndIds.push_back(tmpBndryCndId);
      m_bndryCndOffsets.push_back(counter);
      m_bndryCndSegIds.push_back(tmpSegmentId);
      m_mapBndryCndSegId2Index[tmpSegmentId] = m_bndryCndSegIds.size() - 1;
    }
    m_noBndryCellsPerSegment[tmpSegmentId]++;
    m_solver->a_isBndryCell(m_bndryCells[counter].m_cellId) = true;
    m_solver->a_bndId(m_bndryCells[counter].m_cellId) = counter;

    counter++;
  }

  m_bndryCndOffsets.push_back(counter);

  for(MInt i = 0; i < (MInt)(m_bndryCndSegIds.size()); i++) {
    MInt seg_id = m_bndryCndSegIds[i];

    m_log << "    - BC " << m_bndryCndIds[i] << endl;
    m_log << "      * index:      " << i << endl;
    m_log << "      * segment id: " << seg_id << endl;
    m_log << "      * no cells:   " << m_noBndryCellsPerSegment[seg_id] << endl;
    m_log << "      * offsets:    " << m_bndryCndOffsets[i] << " - " << m_bndryCndOffsets[i + 1] - 1 << endl;
  }
  m_log << endl;
}

/** \brief This function sets the BndryCndHandler objects at solver setup
    \adaptation 22.04.2021, Felix Wietbuescher

    It needs the boundary condition to set the boundary condition handler.
    Each new boundary condition needs to be implemented here, too.
  */
template <MInt nDim>
void FcBndryCnd<nDim>::setBndryCndHandler() {
  TRACE();

  m_log << "  + Setting the boundary condition handler..." << endl << endl;

  // Allocate space for pointer lists which store the boundary functions
  bndryCndHandlerSystemMatrix = new BndryCndHandler[m_bndryCndIds.size()];
  bndryCndHandlerForce = new BndryCndHandler[m_bndryCndIds.size()];

  // Fill the function pointer lists with correct bc functions
  for(MInt i = 0; i < (MInt)(m_bndryCndIds.size()); i++) {
    switch(m_bndryCndIds[i]) {
      case 0:
        bndryCndHandlerSystemMatrix[i] = &FcBndryCnd::bc0;
        bndryCndHandlerForce[i] = &FcBndryCnd::bc0;
        break;

        //----------------------------------------------------
        // fixation boundary conditions
      case 8010:
        bndryCndHandlerSystemMatrix[i] = &FcBndryCnd::bc8010;
        bndryCndHandlerForce[i] = &FcBndryCnd::bc0;
        break;

        //----------------------------------------------------
        // fixation boundary conditions
      case 8011:
        bndryCndHandlerSystemMatrix[i] = &FcBndryCnd::bc8011;
        bndryCndHandlerForce[i] = &FcBndryCnd::bc0;
        break;

        //----------------------------------------------------
        // fixation boundary conditions
      case 8012:
        bndryCndHandlerSystemMatrix[i] = &FcBndryCnd::bc8012;
        bndryCndHandlerForce[i] = &FcBndryCnd::bc0;
        break;

        //----------------------------------------------------
        // displacement boundary conditions
      case 8020:
        bndryCndHandlerSystemMatrix[i] = &FcBndryCnd::bc8020;
        bndryCndHandlerForce[i] = &FcBndryCnd::bc0;
        break;

        //----------------------------------------------------
        // imposed loads boundary conditions
      case 8030:
        bndryCndHandlerSystemMatrix[i] = &FcBndryCnd::bc0;
        bndryCndHandlerForce[i] = &FcBndryCnd::bc8030;
        break;

        //----------------------------------------------------
        // imposed loads boundary conditions
      case 8031:
        bndryCndHandlerSystemMatrix[i] = &FcBndryCnd::bc0;
        bndryCndHandlerForce[i] = &FcBndryCnd::bc8031;
        break;

        //----------------------------------------------------
        // imposed loads boundary conditions
      case 8032:
        bndryCndHandlerSystemMatrix[i] = &FcBndryCnd::bc0;
        bndryCndHandlerForce[i] = &FcBndryCnd::bc8032;
        break;

      case 8035:
        bndryCndHandlerSystemMatrix[i] = &FcBndryCnd::bc0;
        bndryCndHandlerForce[i] = &FcBndryCnd::bc8035;
        break;

      default:
        stringstream errorMessage;
        errorMessage << " FcBndryCnd::setBndryCndHandler : Unknown boundary condition " << m_bndryCndIds[i]
                     << " exiting program.";
        mTerm(1, AT_, errorMessage.str());
    }
  }
}

/** \brief This function applies the fixation boundary condition
 *  \author Moritz Waldmann
 *
 *  For every cell with bcId 8010 part of the displacements of the
 *  corresponding nodes in the assembly vector are set to zero.
 *  The nodes can be identified via the segment normal vector.
 *  The directions which have to be fixated are set
 *  in the segment fixaction direction vector (segFixationDirs).
 */
template <MInt nDim>
void FcBndryCnd<nDim>::bc8010(MInt index) {
  TRACE();

  MInt segId = m_bndryCndSegIds[index];
  cout << "bc8010(" << index << "), offset [" << m_bndryCndOffsets[index] << ", " << m_bndryCndOffsets[index + 1]
       << "; diff = " << m_bndryCndOffsets[index + 1] - m_bndryCndOffsets[index] << "]" << endl;

  // Run through all corresponding boundary cells and change the node displacements
  for(MInt i = m_bndryCndOffsets[index]; i < m_bndryCndOffsets[index + 1]; i++) {
    const MInt cellId = m_bndryCells[i].m_cellId;
    const MFloat halfLength = m_solver->c_cellLengthAtCell(cellId) * F1B2;

    const MInt noNodes = m_solver->a_noNodes(cellId);

    MFloatScratchSpace penalty(nDim * m_solver->a_noNodes(cellId), nDim * m_solver->a_noNodes(cellId), AT_, "penalty");
    penalty.fill(F0);
    for(MInt node = 0; node < noNodes; node++) {
      MFloat nodalCoord[nDim];
      m_solver->getCoordinatesOfNode(node, cellId, nodalCoord);

      for(MInt t = 0; t < (MInt)m_bndryCells[i].m_cutFaces.size(); t++) {
        if(m_bndryCells[i].m_segmentIdOfCutFace[t] != segId) continue;

        for(MInt d = F0; d < nDim; d++) {
          MFloat dist = F0;
          for(MInt dim = 0; dim < nDim; dim++) {
            const MFloat normal = m_bndryCells[i].m_faceNormals[t][dim];
            const MFloat diff = (nodalCoord[dim] - m_bndryCells[i].m_cutFaces[t][d * nDim + dim]);
            dist += normal * normal * diff * diff;
          }
          const MFloat eps = halfLength * 1e-3;
          if(dist < (eps * eps)) {
            for(MInt dim = 0; dim < nDim; dim++) {
              const MInt pos = node * nDim + dim;
              if(m_bndryCells[i].m_directionOfAction[dim] > F0) {
                penalty(pos, pos) = m_kFactor;
              }
            }
          }
        }
      }
    }
    m_solver->computeAssembledBndryMatrix(penalty, cellId);
  }
}

/** \brief This function applies the fixation boundary condition
 *  \author Moritz Waldmann
 *
 *  For every cell with bcId 8011 part of the displacements of the
 *  corresponding nodes in the assembly vector are set to zero.
 *  The nodes can be identified via the segment normal vector.
 *  The directions which have to be fixated are set
 *  in the segment fixaction direction vector (segFixationDirs).
 */
template <MInt nDim>
void FcBndryCnd<nDim>::bc8011(MInt index) {
  TRACE();

  MInt segId = m_bndryCndSegIds[index];
  cout << "bc8011(" << index << "), offset [" << m_bndryCndOffsets[index] << ", " << m_bndryCndOffsets[index + 1]
       << "; diff = " << m_bndryCndOffsets[index + 1] - m_bndryCndOffsets[index] << "]" << endl;

  // Run through all corresponding boundary cells and find the nodes lying on the boundary surface.
  for(MInt i = m_bndryCndOffsets[index]; i < m_bndryCndOffsets[index + 1]; i++) {
    const MInt cellId = m_bndryCells[i].m_cellId;
    const MFloat halfLength = m_solver->c_cellLengthAtCell(cellId) * F1B2;
    std::vector<MInt> nodeList;
    getStlNodeList(cellId, nodeList);

    MInt noNodes = m_solver->a_noNodes(cellId);
    MIntScratchSpace nodeOnSurface(noNodes, AT_, "nodeOnSurface");
    nodeOnSurface.fill(0);
    // Calculate the distance of each node to the boundary surface.
    // If it is below the eps, the node is located on the surface.
    MFloat eps = 1e-6;
    for(MInt node = 0; node < noNodes; node++) {
      MFloat nodalCoord[nDim];
      m_solver->getCoordinatesOfNode(node, cellId, nodalCoord);

      MFloat d[nDim];
      MFloat e[nDim];
      MFloat radius = F0;
      for(MInt dim = 0; dim < nDim; dim++) {
        d[dim] = nodalCoord[dim] - m_bndryCells[i].m_avgFaceNormal[dim] * halfLength;
        e[dim] = nodalCoord[dim] + m_bndryCells[i].m_avgFaceNormal[dim] * halfLength;
        radius += (d[dim] - nodalCoord[dim]) * (d[dim] - nodalCoord[dim]);
      }
      radius = sqrt(radius);

      MBool hasCut = false;
      for(MInt n = 0; n < (signed)nodeList.size(); n++) {
        if(m_solver->m_geometry->elements[nodeList[n]].m_segmentId != segId) continue;

        MFloat distance = 0.0;
        MFloat** const v = m_solver->m_geometry->elements[nodeList[n]].m_vertices;
        hasCut = m_solver->m_geometry->getLineTriangleIntersectionSimpleDistance(d, e, v[0], v[1], v[2], &distance);
        if(approx(distance, radius, eps)) {
          hasCut = true;
          break;
        }
      }
      if(hasCut) {
        nodeOnSurface(node) = 1;
      }
    }

    // Since we have identified all nodes on the boundary surface,
    // the penalty factor can now be added to the system matrix
    MFloat z[nDim] = {F0};
    MFloatScratchSpace penalty(nDim * m_solver->a_noNodes(cellId), nDim * m_solver->a_noNodes(cellId), AT_, "penalty");
    penalty.fill(F0);

    for(MInt j = 0; j < (m_solver->a_pRfnmnt(cellId) + 2); j++) {
      for(MInt k = 0; k < (m_solver->a_pRfnmnt(cellId) + 2); k++) {
        MInt loop = (nDim == 3) ? (m_solver->a_pRfnmnt(cellId) + 2) : 1;
        for(MInt l = 0; l < loop; l++) {
          MInt nodePos[nDim] = {0};
          nodePos[0] = j;
          nodePos[1] = k;
          if(nDim == 3) nodePos[2] = l;

          MBool wrongSide = false;
          for(MInt d = 0; d < nDim; d++) {
            z[d] = m_solver->m_gaussPoints[m_solver->a_pRfnmnt(cellId)][nodePos[d]];
            if(!approx(m_bndryCells[i].m_avgFaceNormal[d], F0, 1e-6)) {
              if(z[d] < F0) {
                wrongSide = true;
              }
            }
          }
          if(wrongSide) continue;

          // Calculate the lagrange interpolation factor
          MFloatScratchSpace L_coef(nDim, m_solver->a_noNodes(cellId) * nDim, AT_, "L_coef");
          L_coef.fill(F1);
          m_solver->getDisplacementInterpolationMatrix(cellId, z, L_coef);

          // Remove lagrange interpolation factor for points outside of the surface of interest.
          for(MInt node = 0; node < m_solver->a_noNodes(cellId); node++) {
            if(nodeOnSurface(node)) continue;
            for(MInt d = 0; d < nDim; d++) {
              L_coef(d, node * nDim + d) = F0;
            }
          }

          // Calculate the integral penalty = int_A L_coef * k * L_coef^T dA
          for(MInt n1 = 0; n1 < m_solver->a_noNodes(cellId); n1++) {
            for(MInt n2 = 0; n2 < m_solver->a_noNodes(cellId); n2++) {
              for(MInt d = 0; d < nDim; d++) {
                MFloat k_factor = (fabs(m_bndryCells[i].m_directionOfAction[d]) > m_solver->m_eps) ? m_kFactor : F0;
                MFloat L1 = F1;
                MFloat L2 = F1;
                for(MInt dim = 0; dim < nDim; dim++) {
                  if(fabs(m_bndryCells[i].m_avgFaceNormal[dim]) > F0) continue;
                  L1 *= L_coef(dim, n1);
                  L2 *= L_coef(dim, n2);
                }
                penalty(n1 * nDim + d, n2 * nDim + d) += (L1 * k_factor * L2);
              }
            }
          }
        }
      }
    }
    m_solver->computeAssembledBndryMatrix(penalty, cellId);
  }
}

/** \brief This function applies the fixation boundary condition
 *  \author Moritz Waldmann
 *
 *  For every cell with bcId 8012 part of the displacements of the
 *  corresponding nodes in the assembly vector are set to zero.
 *  The nodes can be identified via the segment normal vector.
 *  The directions which have to be fixated are set
 *  in the segment fixaction direction vector (segFixationDirs).
 */
template <MInt nDim>
void FcBndryCnd<nDim>::bc8012(MInt index) {
  TRACE();

  MInt segId = m_bndryCndSegIds[index];
  const MInt noTriEdges = (nDim == 3) ? 3 : 1;
  cout << "bc8012(" << index << "), offset [" << m_bndryCndOffsets[index] << ", " << m_bndryCndOffsets[index + 1]
       << "; diff = " << m_bndryCndOffsets[index + 1] - m_bndryCndOffsets[index] << "]" << endl;

  // Run through all corresponding boundary cells and change the node displacements
  for(MInt i = m_bndryCndOffsets[index]; i < m_bndryCndOffsets[index + 1]; i++) {
    const MInt cellId = m_bndryCells[i].m_cellId;
    // Now we can use the triangles for integration
    const MInt dof = nDim * m_solver->a_noNodes(cellId);
    for(MInt t = 0; t < (MInt)m_bndryCells[i].m_cutFaces.size(); t++) {
      if(m_bndryCells[i].m_segmentIdOfCutFace[t] != segId) continue;
      std::vector<MFloat> transformedTriangle;
      MFloat normal[nDim] = {F0};
      for(MInt d = 0; d < nDim; d++) {
        MFloat x[nDim] = {F0};
        MFloat z[nDim] = {F0};
        for(MInt d1 = 0; d1 < nDim; d1++) {
          x[d1] = m_bndryCells[i].m_cutFaces[t][d * nDim + d1];
        }
        m_solver->transformToLocal(cellId, x, z);
        for(MInt d1 = 0; d1 < nDim; d1++) {
          transformedTriangle.push_back(z[d1]);
        }
        normal[d] = m_bndryCells[i].m_faceNormals[t][d];
      }
      MFloat determinant = solveIntegrationOnTriangle(transformedTriangle, normal);
      MFloatScratchSpace penalty(dof, dof, AT_, "penalty");

      MFloat weight = F1B6;
      for(MInt e = 0; e < noTriEdges; e++) {
        MFloatScratchSpace L_coef(nDim, dof, AT_, "L_coef");
        L_coef.fill(F1);
        MFloat z[nDim] = {F0};
        for(MInt d = 0; d < nDim; d++) {
          MInt e1 = e;
          MInt e2 = ((e + 1) < nDim) ? (e + 1) : 0;
          z[d] = (transformedTriangle[e1 * nDim + d] + transformedTriangle[e2 * nDim + d]) * F1B2;
        }
        m_solver->getDisplacementInterpolationMatrix(cellId, z, L_coef);
        if(m_solver->m_testRun && m_solver->m_polyDeg < 1) {
          m_solver->getDisplacementInterpolationMatrixDebug(cellId, z, L_coef);
        }
        // Calculate the integral penalty = int_A L_coef * k * L_coef^T dA
        for(MInt n1 = 0; n1 < dof; n1++) {
          for(MInt n2 = 0; n2 < dof; n2++) {
            MFloat product = F0;
            for(MInt d = 0; d < nDim; d++) {
              MFloat k_factor = (fabs(m_bndryCells[i].m_directionOfAction[d]) > m_solver->m_eps) ? m_kFactor : F0;
              product += L_coef(d, n1) * k_factor * L_coef(d, n2);
            }
            penalty(n1, n2) += product * determinant * weight;
          }
        }
      }
      m_solver->computeAssembledBndryMatrix(penalty, cellId);
    }
  }
}

/** \brief This function applies a non-zero displacement boundary condition
 *  \author Moritz Waldmann
 *
 */
template <MInt nDim>
void FcBndryCnd<nDim>::bc8020(MInt index) {
  TRACE();

  MInt segId = m_bndryCndSegIds[index];
  cout << "bc8020(" << index << ", " << segId << "), offset [" << m_bndryCndOffsets[index] << ", "
       << m_bndryCndOffsets[index + 1] << "; diff = " << m_bndryCndOffsets[index + 1] - m_bndryCndOffsets[index] << "]"
       << endl;

  MInt cellId = 0;
  for(MInt i = m_bndryCndOffsets[index]; i < m_bndryCndOffsets[index + 1]; i++) {
    cellId = m_bndryCells[i].m_cellId;
    cout << "bndryId = " << i << ", cellId = " << cellId << ", primeBndryId = " << m_solver->a_bndId(cellId) << endl;
  }
}

/** \brief This function applies the loads boundary condition
 *  \author Moritz Waldmann
 *
 *  For every cell with bcId 8030 a force is applied to
 *  corresponding nodes in the assembly vector. The nodes can be
 *  identified via the segment normal vector.
 *  The directions of the force are set in the segment
 *  load direction vector (segFixationDirs). An the strength
 *  of the load is set in the segment load value vector.
 */
template <MInt nDim>
void FcBndryCnd<nDim>::bc8030(MInt index) {
  TRACE();

  MInt segId = m_bndryCndSegIds[index];
  cout << "bc8030(" << index << "), offset [" << m_bndryCndOffsets[index] << ", " << m_bndryCndOffsets[index + 1]
       << "; diff = " << m_bndryCndOffsets[index + 1] - m_bndryCndOffsets[index] << "]" << endl;

  for(MInt i = m_bndryCndOffsets[index]; i < m_bndryCndOffsets[index + 1]; i++) {
    const MInt cellId = m_bndryCells[i].m_cellId;
    const MFloat halfLength = m_solver->c_cellLengthAtCell(cellId) * F1B2;
    std::vector<MInt> nodeList;
    getStlNodeList(cellId, nodeList);

    MInt noNodes = m_solver->a_noNodes(cellId);
    for(MInt node = 0; node < noNodes; node++) {
      MFloat nodalCoord[nDim];
      m_solver->getCoordinatesOfNode(node, cellId, nodalCoord);

      MFloat d[nDim];
      MFloat e[nDim];
      MFloat radius = F0;
      for(MInt dim = 0; dim < nDim; dim++) {
        d[dim] = nodalCoord[dim] - m_bndryCells[i].m_avgFaceNormal[dim] * halfLength;
        e[dim] = nodalCoord[dim] + m_bndryCells[i].m_avgFaceNormal[dim] * halfLength;
        radius += (d[dim] - nodalCoord[dim]) * (d[dim] - nodalCoord[dim]);
      }
      radius = sqrt(radius);

      for(MInt n = 0; n < (signed)nodeList.size(); n++) {
        if(m_solver->m_geometry->elements[nodeList[n]].m_segmentId != segId) continue;

        MFloat distance = 0.0;
        MFloat** const v = m_solver->m_geometry->elements[nodeList[n]].m_vertices;
        const MBool hasCut =
            m_solver->m_geometry->getLineTriangleIntersectionSimpleDistance(d, e, v[0], v[1], v[2], &distance);
        std::ignore = hasCut;
        MInt nodeId = m_solver->a_nodeIdsLocal(cellId, node);
        if(approx(distance, radius, 1e-6)) {
          for(MInt dim = 0; dim < nDim; dim++) {
            m_solver->m_externalLoadVector[nodeId * nDim + dim] =
                m_bndryCells[i].m_directionOfAction[dim] * m_bndryCells[i].m_load[dim];
          }
        }
      }
    }
  }
}

/** \brief This function applies the loads boundary condition
 *  \author Moritz Waldmann
 *
 *  For every cell with bcId 8031 a force is applied to
 *  corresponding nodes in the assembly vector. The nodes can be
 *  identified via the segment normal vector.
 *  The directions of the force are set in the segment
 *  load direction vector (segFixationDirs). An the strength
 *  of the load is set in the segment load value vector.
 */
template <MInt nDim>
void FcBndryCnd<nDim>::bc8031(MInt index) {
  TRACE();

  MInt segId = m_bndryCndSegIds[index];
  const MInt noTriEdges = (nDim == 3) ? 3 : 1;
  cout << "bc8031(" << index << "), offset [" << m_bndryCndOffsets[index] << ", " << m_bndryCndOffsets[index + 1]
       << "; diff = " << m_bndryCndOffsets[index + 1] - m_bndryCndOffsets[index] << "]" << endl;


  // Run through all corresponding boundary cells and change the node displacements
  for(MInt i = m_bndryCndOffsets[index]; i < m_bndryCndOffsets[index + 1]; i++) {
    const MInt cellId = m_bndryCells[i].m_cellId;
    const MFloat cellLength = m_solver->c_cellLengthAtCell(cellId);

    // Now we can use the triangles for integration
    for(MInt t = 0; t < (MInt)m_bndryCells[i].m_cutFaces.size(); t++) {
      if(m_bndryCells[i].m_segmentIdOfCutFace[t] != segId) continue;
      std::vector<MFloat> transformedTriangle;
      MFloat normal[nDim] = {F0};
      for(MInt d = 0; d < nDim; d++) {
        MFloat x[nDim] = {F0};
        MFloat z[nDim] = {F0};
        for(MInt d1 = 0; d1 < nDim; d1++) {
          x[d1] = m_bndryCells[i].m_cutFaces[t][d * nDim + d1];
        }
        m_solver->transformToLocal(cellId, x, z);
        for(MInt d1 = 0; d1 < nDim; d1++) {
          transformedTriangle.push_back(z[d1]);
          normal[d1] = m_bndryCells[i].m_faceNormals[t][d];
        }
      }
      MFloat determinant = solveIntegrationOnTriangle(transformedTriangle, normal);
      MFloat weight = F1B6;
      for(MInt e = 0; e < noTriEdges; e++) {
        MFloatScratchSpace L_coef(nDim, m_solver->a_noNodes(cellId) * nDim, AT_, "L_coef");
        L_coef.fill(F1);
        MFloat z[nDim] = {F0};
        for(MInt d = 0; d < nDim; d++) {
          MInt e1 = e;
          MInt e2 = ((e + 1) < nDim) ? (e + 1) : 0;
          z[d] = (transformedTriangle[e1 * nDim + d] + transformedTriangle[e2 * nDim + d]) * F1B2;
        }
        m_solver->getDisplacementInterpolationMatrix(cellId, z, L_coef);
        if(m_solver->m_testRun && m_solver->m_polyDeg < 1) {
          m_solver->getDisplacementInterpolationMatrixDebug(cellId, z, L_coef);
        }

        for(MInt node = 0; node < m_solver->a_noNodes(cellId); node++) {
          for(MInt dim = 0; dim < nDim; dim++) {
            MInt nodeId = m_solver->a_nodeIdsLocal(cellId, node);
            MFloat sigma[nDim] = {F0};
            for(MInt d = 0; d < nDim; d++) {
              sigma[d] = m_bndryCells[i].m_load[d];
            }
            MFloat load = F0;
            for(MInt d = 0; d < nDim; d++) {
              load += L_coef(d, node * nDim + dim) * sigma[d] * determinant * weight;
            }
            m_solver->m_externalLoadVector[nodeId * nDim + dim] += (F1B4 * cellLength * cellLength) * load;
          }
        }
      }
    }
  }
}

/** \brief This function applies the loads boundary condition
 *  \author Moritz Waldmann
 *
 *  For every cell with bcId 8032 a force is applied to
 *  corresponding nodes in the assembly vector. The nodes can be
 *  identified via the segment normal vector.
 *  The directions of the force are set in the segment
 *  load direction vector (segFixationDirs). An the strength
 *  of the load is set in the segment load value vector.
 *  The force is specific for a plate under uniaxial force.
 */
template <MInt nDim>
void FcBndryCnd<nDim>::bc8032(MInt index) {
  TRACE();

  MInt segId = m_bndryCndSegIds[index];
  const MInt noTriEdges = (nDim == 3) ? 3 : 1;
  cout << "bc8032(" << index << "), offset [" << m_bndryCndOffsets[index] << ", " << m_bndryCndOffsets[index + 1]
       << "; diff = " << m_bndryCndOffsets[index + 1] - m_bndryCndOffsets[index] << "]" << endl;


  // Run through all corresponding boundary cells and change the node displacements
  for(MInt i = m_bndryCndOffsets[index]; i < m_bndryCndOffsets[index + 1]; i++) {
    const MInt cellId = m_bndryCells[i].m_cellId;
    const MFloat cellLength = m_solver->c_cellLengthAtCell(cellId);

    // Now we can use the triangles for integration
    for(MInt t = 0; t < (MInt)m_bndryCells[i].m_cutFaces.size(); t++) {
      if(m_bndryCells[i].m_segmentIdOfCutFace[t] != segId) continue;
      std::vector<MFloat> transformedTriangle;
      MFloat normal[nDim] = {F0};
      for(MInt d = 0; d < nDim; d++) {
        MFloat x[nDim] = {F0};
        MFloat z[nDim] = {F0};
        for(MInt d1 = 0; d1 < nDim; d1++) {
          x[d1] = m_bndryCells[i].m_cutFaces[t][d * nDim + d1];
        }
        m_solver->transformToLocal(cellId, x, z);
        for(MInt d1 = 0; d1 < nDim; d1++) {
          transformedTriangle.push_back(z[d1]);
          normal[d1] = m_bndryCells[i].m_faceNormals[t][d];
        }
      }
      MFloat determinant = solveIntegrationOnTriangle(transformedTriangle, normal);
      MFloat weight = F1B6;
      for(MInt e = 0; e < noTriEdges; e++) {
        MFloatScratchSpace L_coef(nDim, m_solver->a_noNodes(cellId) * nDim, AT_, "L_coef");
        L_coef.fill(F1);
        MFloat z[nDim] = {F0};
        MFloat coord[nDim] = {F0};
        for(MInt d = 0; d < nDim; d++) {
          MInt e1 = e;
          MInt e2 = ((e + 1) < nDim) ? (e + 1) : 0;
          z[d] = (transformedTriangle[e1 * nDim + d] + transformedTriangle[e2 * nDim + d]) * F1B2;
          coord[d] =
              (m_bndryCells[i].m_cutFaces[t][e1 * nDim + d] + m_bndryCells[i].m_cutFaces[t][e2 * nDim + d]) * F1B2;
        }
        m_solver->getDisplacementInterpolationMatrix(cellId, z, L_coef);
        if(m_solver->m_testRun && m_solver->m_polyDeg < 1) {
          m_solver->getDisplacementInterpolationMatrixDebug(cellId, z, L_coef);
        }

        MFloat x = coord[0];
        MFloat y = coord[1];
        MFloat r_squared = x * x + y * y;
        MFloat cosTheta = x / sqrt(r_squared);
        MFloat sinTheta = y / sqrt(r_squared);
        MFloat cos2Theta = F2 * cosTheta * cosTheta - 1;
        MFloat sin2Theta = F2 * cosTheta * sinTheta;
        MFloat sigmaRR =
            m_bndryCells[i].m_load[0] * F1B2 * (F1 - F1 / r_squared)
            + m_bndryCells[i].m_load[0] * F1B2 * (F1 + F3 / (r_squared * r_squared) - F4 / r_squared) * cos2Theta;
        MFloat sigmaTT = m_bndryCells[i].m_load[1] * F1B2 * (F1 + F1 / r_squared)
                         - m_bndryCells[i].m_load[1] * F1B2 * (F1 + F3 / (r_squared * r_squared)) * cos2Theta;
        MFloat sigmaRT =
            -m_bndryCells[i].m_load[0] * F1B2 * (F1 + F2 / r_squared - F3 / (r_squared * r_squared)) * sin2Theta;
        MFloat sigmaXX =
            sigmaRR * cosTheta * cosTheta + sigmaTT * sinTheta * sinTheta - F2 * sigmaRT * cosTheta * sinTheta;
        MFloat sigmaYY =
            sigmaTT * cosTheta * cosTheta + sigmaRR * sinTheta * sinTheta + F2 * sigmaRT * cosTheta * sinTheta;
        MFloat sigmaXY =
            (sigmaRR - sigmaTT) * cosTheta * sinTheta + sigmaRT * (cosTheta * cosTheta - sinTheta * sinTheta);

        for(MInt node = 0; node < m_solver->a_noNodes(cellId); node++) {
          for(MInt dim = 0; dim < nDim; dim++) {
            MInt nodeId = m_solver->a_nodeIdsLocal(cellId, node);
            MFloat sigma[nDim] = {F0};
            if(approx(m_bndryCells[i].m_avgFaceNormal[0], F0, 1e-6)) {
              sigma[0] = sigmaXY;
              sigma[1] = sigmaYY;
            } else {
              sigma[0] = sigmaXX;
              sigma[1] = sigmaXY;
            }
            MFloat load = F0;
            for(MInt d = 0; d < nDim; d++) {
              load += L_coef(d, node * nDim + dim) * sigma[d] * determinant * weight;
            }
            m_solver->m_externalLoadVector[nodeId * nDim + dim] += (F1B4 * cellLength * cellLength) * load;
          }
        }
      }
    }
  }
}

/** \brief This function applies the loads boundary condition
 *  \author Moritz Waldmann
 *
 *  For every cell with bcId 8035 a force is applied to
 *  corresponding nodes in the assembly vector.
 */
template <MInt nDim>
void FcBndryCnd<nDim>::bc8035(MInt index) {
  TRACE();

  MInt segId = m_bndryCndSegIds[index];
  const MInt noTriEdges = (nDim == 3) ? 3 : 1;
  cout << "bc8035(" << index << "), offset [" << m_bndryCndOffsets[index] << ", " << m_bndryCndOffsets[index + 1]
       << "; diff = " << m_bndryCndOffsets[index + 1] - m_bndryCndOffsets[index] << "]" << endl;


  // Run through all corresponding boundary cells and change the node displacements
  for(MInt i = m_bndryCndOffsets[index]; i < m_bndryCndOffsets[index + 1]; i++) {
    const MInt cellId = m_bndryCells[i].m_cellId;
    const MFloat cellLength = m_solver->c_cellLengthAtCell(cellId);

    // Now we can use the triangles for integration
    for(MInt t = 0; t < (MInt)m_bndryCells[i].m_cutFaces.size(); t++) {
      if(m_bndryCells[i].m_segmentIdOfCutFace[t] != segId) continue;
      std::vector<MFloat> transformedTriangle;
      MFloat normal[nDim] = {F0};
      for(MInt d = 0; d < nDim; d++) {
        MFloat x[nDim] = {F0};
        MFloat z[nDim] = {F0};
        for(MInt d1 = 0; d1 < nDim; d1++) {
          x[d1] = m_bndryCells[i].m_cutFaces[t][d * nDim + d1];
        }
        m_solver->transformToLocal(cellId, x, z);
        for(MInt d1 = 0; d1 < nDim; d1++) {
          transformedTriangle.push_back(z[d1]);
        }
        normal[d] = m_bndryCells[i].m_faceNormals[t][d];
      }
      MFloat determinant = solveIntegrationOnTriangle(transformedTriangle, normal);
      MFloat weight = F1B6;
      for(MInt e = 0; e < noTriEdges; e++) {
        MFloatScratchSpace L_coef(nDim, m_solver->a_noNodes(cellId) * nDim, AT_, "L_coef");
        L_coef.fill(F1);
        MFloat z[nDim] = {F0};
        for(MInt d = 0; d < nDim; d++) {
          MInt e1 = e;
          MInt e2 = ((e + 1) < nDim) ? (e + 1) : 0;
          z[d] = (transformedTriangle[e1 * nDim + d] + transformedTriangle[e2 * nDim + d]) * F1B2;
        }
        m_solver->getDisplacementInterpolationMatrix(cellId, z, L_coef);
        if(m_solver->m_testRun && m_solver->m_polyDeg < 1) {
          m_solver->getDisplacementInterpolationMatrixDebug(cellId, z, L_coef);
        }

        for(MInt node = 0; node < m_solver->a_noNodes(cellId); node++) {
          for(MInt dim = 0; dim < nDim; dim++) {
            MInt nodeId = m_solver->a_nodeIdsLocal(cellId, node);
            MFloat sigma[nDim] = {F0};
            for(MInt d = 0; d < nDim; d++) {
              sigma[d] = normal[d] * 1.0;
            }
            MFloat load = F0;
            for(MInt d = 0; d < nDim; d++) {
              load += L_coef(d, node * nDim + dim) * sigma[d] * determinant * weight;
            }
            m_solver->m_externalLoadVector[nodeId * nDim + dim] += (F1B4 * cellLength * cellLength) * load;
          }
        }
      }
    }
  }
}

/** \brief Boundary condition for free surfaces
 *  \author Moritz Waldmann
 *
 *  At these surfaces no computations are required.
 */
template <MInt nDim>
void FcBndryCnd<nDim>::bc0(MInt index) {
  TRACE();
  std::ignore = index;
}

/** \brief Execution of the displacement boundary conditions
 *  \author Moritz Waldmann
 *
 */
template <MInt nDim>
void FcBndryCnd<nDim>::updateSystemMatrix() {
  TRACE();

  for(MInt i = 0; i < (MInt)(m_bndryCndIds.size()); i++) {
    (this->*bndryCndHandlerSystemMatrix[i])(i);
  }
}

/** \brief Execution of the force boundary conditions
 *  \author Moritz Waldmann
 *
 */
template <MInt nDim>
void FcBndryCnd<nDim>::updateForceVector() {
  TRACE();

  for(MInt i = 0; i < (MInt)(m_bndryCndIds.size()); i++) {
    (this->*bndryCndHandlerForce[i])(i);
  }
}

/** \brief Calculation of reaction forces
 *  \author Moritz Waldmann
 *
 *  The Reaction forces are defined as the inner forces, which occur at the
 *  fixed boundaries, i.e., boundaries with id < 8020.
 */
template <MInt nDim>
void FcBndryCnd<nDim>::calcReactionForces() {
  TRACE();

  for(MInt index = 0; index < (MInt)(m_bndryCndSegIds.size()); index++) {
    const MInt bc_id = m_bndryCndIds[index];
    const MInt segId = m_bndryCndSegIds[index];

    if(bc_id == 0) continue;
    if(bc_id > 8020) continue;

    for(MInt i = m_bndryCndOffsets[index]; i < m_bndryCndOffsets[index + 1]; i++) {
      const MInt cellId = m_bndryCells[i].m_cellId;
      const MFloat halfLength = m_solver->c_cellLengthAtCell(cellId) * F1B2;

      const MInt noNodes = m_solver->a_noNodes(cellId);

      for(MInt node = 0; node < noNodes; node++) {
        std::array<MFloat, nDim> nodalCoord{};
        m_solver->getCoordinatesOfNode(node, cellId, nodalCoord.data());

        for(MInt t = 0; t < (MInt)m_bndryCells[i].m_cutFaces.size(); t++) {
          if(m_bndryCells[i].m_segmentIdOfCutFace[t] != segId) continue;

          for(MInt d = F0; d < nDim; d++) {
            MFloat dist = F0;
            for(MInt dim = 0; dim < nDim; dim++) {
              const MFloat normal = m_bndryCells[i].m_faceNormals[t][dim];
              const MFloat diff = (nodalCoord[dim] - m_bndryCells[i].m_cutFaces[t][d * nDim + dim]);
              dist += normal * normal * diff * diff;
            }
            const MFloat eps = halfLength * 1e-3;
            if(dist < (eps * eps)) {
              for(MInt dim = 0; dim < nDim; dim++) {
                const MInt pos = m_solver->a_nodeIdsLocal(cellId, node) * nDim + dim;
                if(m_bndryCells[i].m_directionOfAction[dim] > F0) {
                  m_solver->m_reactionForceVector[pos] = m_solver->m_internalLoadVector[pos];
                }
              }
            }
          }
        }
      }
    }
  }
}

/** \brief Execution of the subcell integration
 *  \author Moritz Waldmann
 *
 */
template <MInt nDim>
void FcBndryCnd<nDim>::subCellIntegration(const MInt subCellLvl, MFloat* subCellParentCoord, const MInt subCellPos,
                                          const MInt pCellId, MFloat** Ke_final) {
  TRACE();

  // do the subcell integration only until the max depth is reached
  if(subCellLvl > m_solver->a_maxSubCellLvl(pCellId)) return;

  std::vector<MInt> nodeList;

  // set the half cell length...
  const MFloat subCellHalfLength = m_solver->c_cellLengthAtCell(pCellId) * FFPOW2(subCellLvl + 1);
  // and the bounding box of the subcell and the center coordinates of the sub cell
  MFloat target[2 * nDim] = {F0};
  MFloat subCellCoord[nDim] = {F0};
  for(MInt d = 0; d < nDim; d++) {
    MFloat dir = Fd::vertexPosition(subCellPos, d);
    subCellCoord[d] = subCellParentCoord[d] + dir * subCellHalfLength;
    target[d] = subCellCoord[d] - subCellHalfLength;
    target[d + nDim] = subCellCoord[d] + subCellHalfLength;
  }

  // increase the childLevel
  const MInt childLevel = subCellLvl + 1;

  // check if the sub cell is intersected by the geometry
  m_solver->m_geometry->getIntersectionElements(target, nodeList, subCellHalfLength, subCellCoord);

  // if the cell is intersected refine the cell further if the maxSubCellLevel is not reached
  // if it is intersected but the maximum depth is reached, calculate the stiffness matrix
  if(nodeList.size() > 0) {
    if(childLevel <= m_solver->a_maxSubCellLvl(pCellId)) {
      for(MInt child = 0; child < IPOW2(nDim); child++) {
        subCellIntegration(childLevel, subCellCoord, child, pCellId, Ke_final);
      }
    } else {
      // calculate the stiffness matrix for this subcell
      MFloatScratchSpace Ke(m_solver->a_noNodes(pCellId) * nDim, m_solver->a_noNodes(pCellId) * nDim, AT_, "Ke");
      // reset element stiffness matrix
      for(MInt m1 = 0; m1 < m_solver->a_noNodes(pCellId) * nDim; m1++) {
        for(MInt m2 = 0; m2 < m_solver->a_noNodes(pCellId) * nDim; m2++) {
          Ke(m1, m2) = F0;
        }
      }

      // Calculate the Element Stiffness Matrix for the subcell
      // and add it to the element stiffness of the root element
      const MFloat alpha = (F1 - m_solver->m_alpha);
      m_solver->getElementMatrix(Ke, pCellId, alpha, subCellLvl, subCellCoord);
      for(MInt m1 = 0; m1 < m_solver->a_noNodes(pCellId) * nDim; m1++) {
        for(MInt m2 = 0; m2 < m_solver->a_noNodes(pCellId) * nDim; m2++) {
          Ke_final[m1][m2] += Ke(m1, m2);
        }
      }
    }
  } else {
    // Check if the subcell is outside or inside the geometry
    // if the subcell is outside nothing more is to be done
    // if the subcell is inside the element stiffness matrix of the subcell is calculated
    MBool outside = m_solver->m_geometry->pointIsInside2(subCellCoord);
    if(outside) {
      return;
    }

    // calculate the stiffness matrix for this subcell
    MFloatScratchSpace Ke(m_solver->a_noNodes(pCellId) * nDim, m_solver->a_noNodes(pCellId) * nDim, AT_, "Ke");
    // reset element stiffness matrix
    for(MInt m1 = 0; m1 < m_solver->a_noNodes(pCellId) * nDim; m1++) {
      for(MInt m2 = 0; m2 < m_solver->a_noNodes(pCellId) * nDim; m2++) {
        Ke(m1, m2) = F0;
      }
    }
    // Calculate the Element Stiffness Matrix for the subcell
    // and add it to the element stiffness of the root element
    const MFloat alpha = (F1 - m_solver->m_alpha);
    m_solver->getElementMatrix(Ke, pCellId, alpha, subCellLvl, subCellCoord);
    for(MInt m1 = 0; m1 < m_solver->a_noNodes(pCellId) * nDim; m1++) {
      for(MInt m2 = 0; m2 < m_solver->a_noNodes(pCellId) * nDim; m2++) {
        Ke_final[m1][m2] += Ke(m1, m2);
      }
    }
  }
}

/** \brief Set the subcell tree depth for cells requiring subcell integration
 *  \author Moritz Waldmann
 *
 *  The depth and the boundary segment, at which the subcell integration is required,
 *  is given in the properties file.
 */
template <MInt nDim>
void FcBndryCnd<nDim>::findCellsRequireSubCellIntegration() {
  TRACE();

  for(MInt cell = 0; cell < m_solver->a_noCells(); cell++) {
    m_solver->a_needsSubCells(cell) = false;
    m_solver->a_maxSubCellLvl(cell) = 0;
  }
  for(MInt cell = 0; cell < (MInt)m_bndryCells.size(); cell++) {
    MInt pCellId = m_bndryCells[cell].m_cellId;
    MInt depth = 0;
    for(MInt noSeg = 0; noSeg < (MInt)m_bndryCells[cell].m_segmentId.size(); noSeg++) {
      MInt segId = m_bndryCells[cell].m_segmentId[noSeg];
      if(m_subCellLayerDepth[segId] > depth) {
        depth = m_subCellLayerDepth[segId];
      }
    }
    if(depth > 0) {
      m_solver->a_needsSubCells(pCellId) = true;
      m_solver->a_maxSubCellLvl(pCellId) = depth;
    }
  }
}

/** \brief Triangulation algorithm of the cut points
 *  \author Moritz Waldmann
 *
 *  Works only in 3D
 */
template <MInt nDim>
std::vector<std::vector<MFloat>>
FcBndryCnd<nDim>::createTrianglesFromCutPoints(std::vector<std::vector<MFloat>> cutPointList, MFloat* triangleNormal) {
  TRACE();

  // Step 1:
  // Calculate the Center of gravity of the cut point list.
  MFloat COG[nDim] = {F0};
  for(MInt d = 0; d < nDim; d++) {
    for(MInt c1 = 0; c1 < (MInt)cutPointList.size(); c1++) {
      COG[d] += cutPointList[c1][d];
    }
    COG[d] /= (MInt)cutPointList.size();
  }

  // Step 2:
  // Sort vector of cut points depending on their angle.
  // The line segment between COG and the first point in the cut point list is the reference.
  // The line segments between COG and the other points in the cut point list are calculated
  // consecutively. The angle between them and the reference line segment is calculated using
  // the scalar product. Since the cosinus in cpp is not defined for 2pi, the cross product
  // of the two vectors is also calculated. The sign of the scalar product of the resulting
  // vector and the triangle normal determines if the angle is <pi or >=pi.
  // Since the first point in the cut point list is the reference point, it is assigned the
  // angle 0.
  std::vector<MFloat> angles;
  angles.push_back(0);
  for(MInt c1 = 1; c1 < (MInt)cutPointList.size(); c1++) {
    MFloat crossProduct[nDim] = {F0};
    MFloat cosPhi = F0;
    MFloat l1 = F0;
    MFloat l2 = F0;
    MFloat lCross = F0;
    MFloat lNormal = F0;

    // Step 2a: Calculate cross product and scalar product
    // Furthermore the length of the different vectors are calculated
    for(MInt d = 0; d < nDim; d++) {
      // Cross product
      crossProduct[d] = F0;
      MInt d1 = ((d + 1) < nDim) ? (d + 1) : 0;
      MInt d2 = ((d1 + 1) < nDim) ? (d1 + 1) : 0;
      crossProduct[d] = ((cutPointList[0][d1] - COG[d1]) * (cutPointList[c1][d2] - COG[d2]))
                        - ((cutPointList[0][d2] - COG[d2]) * (cutPointList[c1][d1] - COG[d1]));

      // Scalar product
      cosPhi += ((cutPointList[0][d] - COG[d]) * (cutPointList[c1][d] - COG[d]));

      // Vector length
      l1 += ((cutPointList[0][d] - COG[d]) * (cutPointList[0][d] - COG[d]));
      l2 += ((cutPointList[c1][d] - COG[d]) * (cutPointList[c1][d] - COG[d]));
      lCross += (crossProduct[d] * crossProduct[d]);
      lNormal += (triangleNormal[d] * triangleNormal[d]);
    }
    cosPhi /= (sqrt(l1) * sqrt(l2));

    // Step 2b: Calculate the angle
    MFloat oppDir = F0;
    for(MInt d = 0; d < nDim; d++) {
      oppDir += (triangleNormal[d] * crossProduct[d]);
    }
    oppDir /= (sqrt(lCross) * sqrt(lNormal));

    // Step 2c: Calculate angle
    // Since it might occur that abs(cosPhi) is slightly larger than 1 due to numerical errors
    // it is rounded in this case. Afterwards, the angle is calculated. If oppDir is <0 the
    // opposite angle is used.
    if(cosPhi > F1 || cosPhi < -F1) {
      cosPhi = round(cosPhi);
    }
    if(oppDir >= F0) {
      angles.push_back(acos(cosPhi));
    } else {
      angles.push_back((2 * PI) - acos(cosPhi));
    }
  }

  // Step 2d: Sort the cut point list from smallest to largest angles.
  for(MInt i = 0; i < (MInt)(angles.size() - 1); i++) {
    for(MInt j = 0; j < (MInt)(angles.size() - i - 1); j++) {
      if(angles.at(j) > angles.at(j + 1)) {
        std::swap(angles.at(j), angles.at(j + 1));
        std::swap(cutPointList.at(j), cutPointList.at(j + 1));
      }
    }
  }

  // Step 3: Triangulate cut points. One triangles consists of 3 points with 3 coordinates.
  std::vector<std::vector<MFloat>> triangleList;
  for(MInt c = 1; c < (MInt)cutPointList.size() - 1; c++) {
    std::vector<MFloat> triangle;
    for(MInt d1 = 0; d1 < nDim; d1++) {
      for(MInt d2 = 0; d2 < nDim; d2++) {
        MInt point = (d1 == 0) ? d1 : (c + d1 - 1);
        triangle.push_back(cutPointList[point][d2]);
      }
    }
    triangleList.push_back(triangle);
  }

  return triangleList;
}

/** \brief Calculation of the surface edges.
 *  \author Moritz Waldmann
 *
 *  2D version of the triangulation algorithm.
 *  Works only in 2D.
 */
template <MInt nDim>
std::vector<std::vector<MFloat>>
FcBndryCnd<nDim>::createEdgesFromCutPoints(std::vector<std::vector<MFloat>> cutPointList) {
  TRACE();

  ASSERT((MInt)cutPointList.size() == 2, "A line can only have two cut points with a cell");
  std::vector<std::vector<MFloat>> edgeList;
  std::vector<MFloat> edge;

  for(MInt c = 0; c < (MInt)cutPointList.size(); c++) {
    for(MInt p = 0; p < (MInt)cutPointList[c].size(); p++) {
      edge.push_back(cutPointList[c][p]);
    }
  }

  ASSERT((MInt)edge.size() == 4, "A edge must have only four coordinates");
  edgeList.push_back(edge);

  return edgeList;
}

/** \brief Gauss quadrature on triangles
 *  \author Moritz Waldmann
 *
 */
template <MInt nDim>
MFloat FcBndryCnd<nDim>::solveIntegrationOnTriangle(std::vector<MFloat> trianglePoints, MFloat* triangleNormal) {
  TRACE();

  std::vector<MFloat> transformedTriPoints;

  MFloat edge1[nDim] = {F0};
  MFloat edge2[nDim] = {F0};
  MFloat result[nDim] = {F0};
  // Calculate the transformed points
  for(MInt d1 = 0; d1 < nDim; d1++) {
    for(MInt d2 = 0; d2 < nDim; d2++) {
      MFloat transformed = trianglePoints[d1 * nDim + d2] - trianglePoints[d2];
      if(d1 == nDim - 2) edge1[d2] = (trianglePoints[d1 * nDim + d2] - trianglePoints[d2]);
      if(d1 == nDim - 1) edge2[d2] = (trianglePoints[d1 * nDim + d2] - trianglePoints[d2]);
      transformedTriPoints.push_back(transformed);
    }
  }

  if(nDim == 3) {
    maia::math::cross(edge1, edge2, result);
  }

  MFloat determinant = F0;
  for(MInt d1 = 0; d1 < nDim; d1++) {
    determinant += (result[d1] * result[d1]);
  }
  determinant = sqrt(determinant);

  // TODO: This is the mathematical correct way to calculate the determinant,
  // but so far, it is not working due to a bug.
  std::ignore = triangleNormal;
  /*
    MFloatScratchSpace transformationMatrix((nDim + 1), (nDim + 1), AT_, "transformationMatrix");
    transformationMatrix.fill(F0);
    // This is the inverse of the translation matrix.
    // The matrix moves one point of the triangle in the coordinate origin
    for(MInt d1 = 0; d1 < nDim + 1; d1++) {
      for(MInt d2 = 0; d2 < nDim; d2++) {
        if(d1 == d2) transformationMatrix(d2, d1) = F1;
        if(d1 == nDim) transformationMatrix(d2, d1) = trianglePoints[d2];
      }
    }
    transformationMatrix(nDim, nDim) = F1;

    //Use Rodrigues rotation formula to rotate triangle in xy-plane
    //Normalize the normal of the triangle
    MFloat lengthN = F0;
    for(MInt d = 0; d < nDim; d++) {
      lengthN += (triangleNormal[d] * triangleNormal[d]);
    }
    for(MInt d = 0; d < nDim; d++) {
      triangleNormal[d] /= sqrt(lengthN);
    }

    //Calculate the cross product of the x-axis and the triangle normal
    //and normalize the result
    MFloat xyPlaneNormal[nDim] = {F0};
    xyPlaneNormal[nDim - 1] = F1;
    MFloat k[nDim] = {F0};
    if(nDim == 3) {
      maia::math::cross(triangleNormal, xyPlaneNormal, k);
    }
    MFloat lengthK = F0;
    for(MInt d = 0; d < nDim; d++) {
      lengthK += (k[d] * k[d]);
    }

    for(MInt d = 0; d < nDim; d++) {
      k[d] /= sqrt(lengthK);
      triangleNormal[d] /= sqrt(lengthN);
    }

    //Calculate the rotation angle
    MFloat dotProduct = F0;
    for(MInt d = 0; d < nDim; d++) {
      dotProduct += triangleNormal[d] * xyPlaneNormal[d];
    }
    MFloat cosTheta = dotProduct;
    MFloat sinTheta = lengthK / (lengthN); //length of xyPlaneNormal is always one

    //Calculate the rotation matrix using Rodrigues rotation formular
    MFloat identityMat[nDim][nDim] = {F0};;
    MFloat kMat[nDim][nDim] = {F0};
    MFloat kMat_square[nDim][nDim] = {F0};
    MFloatScratchSpace R1((nDim + 1), (nDim + 1), AT_, "R1");
    MFloatScratchSpace R2((nDim + 1), (nDim + 1), AT_, "R2");
    R1.fill(F0);
    R2.fill(F0);

    for(MInt d1 = 0; d1 < nDim; d1++) {
      for(MInt d2 = 0; d2 < nDim; d2++) {
        identityMat[d1][d2] = F0;
        if(d1 == d2) identityMat[d1][d2] = F1;
      }
    }

    for(MInt d1 = 0; d1 < nDim; d1++) {
      kMat[d1][d1] = F0;
      MInt d2 = ((d1 + 1) < nDim) ? (d1 + 1) : 0;
      MInt d3 = ((d2 + 1) < nDim) ? (d2 + 1) : 0;
      kMat[d1][d2] = k[d3];
      if(d1 == 1) kMat[d1][d2] *= -F1;
      kMat[d2][d1] = -k[d3];
    }

    for(MInt d1 = 0; d1 < nDim; d1++) {
      for(MInt d2 = 0; d2 < nDim; d2++) {
        kMat_square[d1][d2] = F0;
        for(MInt d3 = 0; d3 < nDim; d3++) {
          kMat_square[d1][d2] += (kMat[d1][d3] * kMat[d3][d2]);
        }
      }
    }

    for(MInt d1 = 0; d1 < nDim; d1++) {
      for(MInt d2 = 0; d2 < nDim; d2++) {
        R1(d1, d2) = identityMat[d1][d2] + (sinTheta * kMat[d1][d2]) + ((F1 + cosTheta) * kMat_square[d1][d2]);
      }
    }

    //To transform the triangle into a unit triangle a third transformation is needed
    //The inverse of this matrix is given by the coordinates of the transformed
    for(MInt d1 = 0; d1 < nDim; d1++) {
      MFloat point[nDim] = {F0};
      for(MInt d2 = 0; d2 < nDim; d2++) {
        for(MInt d3 = 0; d3 < nDim; d3++) {
          point[d2] += transformedTriPoints[d1 * nDim + d3] * R1(d2, d3);
        }
      }
      for(MInt d = 0; d < nDim; d++) {
        transformedTriPoints[d1 * nDim + d] = point[d];
      }
    }

    for(MInt d1 = 0; d1 < nDim; d1++) {
      for(MInt d2 = 0; d2 < nDim; d2++) {
        MInt d = ((d1 + 1) < nDim) ? (d1 + 1) : 0;
        R2(d2, d1) = transformedTriPoints[d * nDim + d2];
      }
    }

    //The nDim + 1 dimension is required for the translation and to
    //allow to calculat the inverse of the matrix
    R2((nDim - 1), (nDim - 1)) = F1;
    R1(nDim, nDim) = F1;
    R2(nDim, nDim) = F1;
    //TODO: Write a inverse function in maiamath.h that can handle
    //2DScratchSpaces
    MFloatScratchSpace R1_inv((nDim + 1) * (nDim + 1), AT_, "R1_inv");
    for(MInt d1 = 0; d1 < nDim + 1; d1++) {
      for(MInt d2 = 0; d2 < nDim + 1; d2++) {
        R1_inv(d1 * (nDim + 1) + d2) = R1(d1, d2);
      }
    }
    maia::math::invert(R1_inv.getPointer(), (nDim + 1), (nDim + 1));

    for(MInt d1 = 0; d1 < nDim + 1; d1++) {
      for(MInt d2 = 0; d2 < nDim + 1; d2++) {
        R1(d1, d2) = R1_inv(d1 * (nDim + 1) + d2);
      }
    }

    //Now we have everything we need, so we calculate the product of the 3 inverse
    //matrices. The product is the transformation matrix that transforms the triangle
    //Gauss points into the cell coordinate system.
    MFloatScratchSpace helper((nDim + 1), (nDim + 1), AT_, "helper");
    maia::math::multiplyMatricesSq(transformationMatrix, R1, helper, nDim + 1);
    maia::math::multiplyMatricesSq(helper, R2, transformationMatrix, nDim + 1);

    //Calculate the determinant of the resulting matrix
    MFloatScratchSpace subMat(nDim, nDim, AT_, "subMat");
    MFloat determinant = F0;
    MFloat sign = -F1;
    for(MInt d1 = 0; d1 < (nDim + 1); d1++) {
      sign *= (-F1);
      MInt cnt = 0;
      for(MInt d2 = 0; d2 < (nDim + 1); d2++) {
        if(d1 == d2) continue;
        for(MInt d3 = 1; d3 < (nDim + 1); d3++) {
          subMat(cnt, d3 - 1) = transformationMatrix(d2, d3);
        }
        cnt++;
      }
      MFloat det = maia::math::determinant(subMat, nDim);
      determinant += (sign * transformationMatrix(d1, 0) * det);
    }
  */
  return determinant;
}

/** \brief Calculation of the cutpoints of each segment
 *  \author Moritz Waldmann
 *
 */
template <MInt nDim>
void FcBndryCnd<nDim>::calculateCutPoints() {
  TRACE();

  // Returns the number of the opposite vertex of an edge
  //        e1
  //   1 - - - - 3
  //   |         |
  // e0|         |e3
  //   |         |
  //   0 - - - - 2
  //        e2
  // For edge 2 the edgeStartPos is vertex number 2 and
  // the edgeEndPos is vertex number 0
  // This is valid in 2D and for the edges in the yz-plane in 3D
  const MInt edgeEndPos[4] = {1, 3, 0, 2};

  // Stores the normals of the cell faces
  // the normal points in -x, x, -y, y, -z, z direction
  MFloat normals[2 * nDim][nDim];
  for(MInt d1 = 0; d1 < 2 * nDim; d1++) {
    for(MInt d2 = 0; d2 < nDim; d2++) {
      if(d1 >= (2 * d2) && d1 < (2 * (d2 + 1))) {
        MInt expo = d1 + 1;
        normals[d1][d2] = pow(-F1, expo);
      } else {
        normals[d1][d2] = F0;
      }
    }
  }

  const MInt noStlElementEdges = (nDim == 3) ? 3 : 1;
  const MInt noElementEdges = (nDim == 3) ? 12 : 4;
  const MFloat eps = m_solver->m_eps;

  struct stlElement {
    MFloat points[nDim][nDim];
    MFloat edges[noStlElementEdges][2 * nDim];
    std::array<MFloat, nDim> normal;
    std::array<MBool, nDim> pointInsideCell;
    MInt noInsidePoints;
    std::array<MInt, noStlElementEdges> noCutsPerEdge;
    MBool cutSurfPerEdge[noStlElementEdges][2 * nDim];
    MInt noCuttedSurfaces;
    std::array<MBool, noStlElementEdges> edgeNotCutted;
    std::array<MBool, noStlElementEdges> edgeOnSurface;
  };

  struct element {
    MFloat edges[noElementEdges][2 * nDim];
    std::array<MBool, noElementEdges> edgeIsCut;
    MInt noCuts;
    std::array<MBool, noElementEdges> edgeOnSurface;
  };

  // Loop over all boundary ids
  for(MInt index = 0; index < (MInt)(m_bndryCndIds.size()); index++) {
    MInt segId = m_bndryCndSegIds[index];

    for(MInt i = m_bndryCndOffsets[index]; i < m_bndryCndOffsets[index + 1]; i++) {
      const MInt cellId = m_bndryCells[i].m_cellId;
      const MFloat cellLength = m_solver->c_cellLengthAtCell(cellId);
      const MFloat halfLength = cellLength * 0.5;

      // 1. Calculate the edges of the cell,
      // i.e., for each edge, the start and end point is stored
      element ele;
      for(MInt e = 0; e < noElementEdges; e++) {
        if(nDim == 2) {
          for(MInt d = 0; d < nDim; d++) {
            MFloat coord = m_solver->c_coordinate(cellId, d);
            ele.edges[e][d] = coord + Fd::vertexPosition(e, d) * halfLength;
            ele.edges[e][d + nDim] = coord + Fd::vertexPosition(edgeEndPos[e], d) * halfLength;
          }
        } else {
          if(e < 4) { // For all edges in the yz-plane in negative x-direction
            for(MInt d = 0; d < nDim; d++) {
              MFloat coord = m_solver->c_coordinate(cellId, d);
              ele.edges[e][d] = coord + Fd::vertexPosition(e, d) * halfLength;
              ele.edges[e][d + nDim] = coord + Fd::vertexPosition(edgeEndPos[e], d) * halfLength;
            }
          } else if(e < 8) { // For all edges in the yz-plane in positive x-direction
            for(MInt d = 0; d < nDim; d++) {
              MFloat coord = m_solver->c_coordinate(cellId, d);
              ele.edges[e][d] = coord + Fd::vertexPosition(e - 4, d) * halfLength;
              ele.edges[e][d + nDim] = coord + Fd::vertexPosition(e, d) * halfLength;
            }
          } else {
            for(MInt d = 0; d < nDim; d++) { // For all edges, which are parallel to the x-axis
              MFloat coord = m_solver->c_coordinate(cellId, d);
              ele.edges[e][d] = coord + Fd::vertexPosition(e - 4, d) * halfLength;
              MInt h = edgeEndPos[e - 8] + 4;
              ele.edges[e][d + nDim] = coord + Fd::vertexPosition(h, d) * halfLength;
            }
          }
        }
      }

      // 2. Find all stl elements that cut the element
      MFloat target[2 * nDim]; // bounding box of the currentId
      for(MInt j = 0; j < nDim; j++) {
        target[j] = m_solver->c_coordinate(cellId, j) - halfLength;
        target[j + nDim] = m_solver->c_coordinate(cellId, j) + halfLength;
      }
      std::vector<MInt> nodeList;
      getStlNodeList(cellId, nodeList);
      for(MInt n = 0; n < (MInt)nodeList.size(); n++) {
        std::vector<std::vector<MFloat>> cutPointList;

        // If STL-element does not belong to the boundary, continue
        if(m_solver->m_geometry->elements[nodeList[n]].m_segmentId != segId) {
          continue;
        }

        // a). create stl elements struct and calculate the edges of the STL-element, i.e.,
        // for each edge, the start and end point is stored
        stlElement stl;
        for(MInt e = 0; e < noStlElementEdges; e++) {
          for(MInt d = 0; d < nDim; d++) {
            stl.edges[e][d] = m_solver->m_geometry->elements[nodeList[n]].m_vertices[e][d];
            if(e + 1 < nDim) {
              stl.edges[e][d + nDim] = m_solver->m_geometry->elements[nodeList[n]].m_vertices[e + 1][d];
            } else {
              stl.edges[e][d + nDim] = m_solver->m_geometry->elements[nodeList[n]].m_vertices[0][d];
            }
          }
        }
        if(nDim == 3) {
          for(MInt d = 0; d < nDim; d++) {
            stl.normal[d] = m_solver->m_geometry->elements[nodeList[n]].m_normal[d];
          }
        } else {
          stl.normal[0] = stl.edges[0][1 + nDim] - stl.edges[0][1];
          stl.normal[1] = stl.edges[0][0 + nDim] - stl.edges[0][0];
        }

        // b). Check if the corner points of the STL-element are inside or outside of the cell
        stl.noInsidePoints = 0;
        for(MInt point = 0; point < nDim; point++) {
          stl.pointInsideCell[point] = false;
          MBool isInside = true;
          for(MInt d = 0; d < nDim; d++) {
            stl.points[point][d] = m_solver->m_geometry->elements[nodeList[n]].m_vertices[point][d];
            if(stl.points[point][d] < (target[d]) || stl.points[point][d] > (target[d + nDim])) {
              isInside = false;
            }
          }
          stl.pointInsideCell[point] = isInside;
          // if a point is located inside the cell, it is added to the cut point list
          if(stl.pointInsideCell[point]) {
            stl.noInsidePoints++;
            std::vector<MFloat> cutPoint;
            for(MInt d = 0; d < nDim; d++) {
              cutPoint.push_back(stl.points[point][d]);
            }
            cutPointList.push_back(cutPoint);
          }
        }

        // c). Check if edge is cutted
        // If start and end point are located inside of the cell, the edge is not cutted.
        // If start and end point are located outside of the cell, the edge might have two cuts.
        for(MInt e = 0; e < noStlElementEdges; e++) {
          MInt point1 = e;
          MInt point2 = (e + 1 < nDim) ? (e + 1) : 0;
          stl.edgeNotCutted[e] = (stl.pointInsideCell[point1] && stl.pointInsideCell[point2]);
          stl.noCutsPerEdge[e] = 0;
          stl.noCuttedSurfaces = 0;
          stl.edgeOnSurface[e] = false;
          for(MInt surface = 0; surface < 2 * nDim; surface++) {
            stl.cutSurfPerEdge[e][surface] = false;
          }
        }

        // d). Check if the edges cut the cell surfaces and count the cuts per edge
        // An edge can cut the cell surfaces twice (if start and end point are
        // outside of the cell), once (if one of the edge corner points is inside
        // and the other is outside), or zero times (if start and end point are
        // outside or inside of the cell).
        // Furthermore, we need to check if the edge is located on the surface.
        // If all points of the stl element are located inside the cell, we dont need to
        // do the next steps
        if(stl.noInsidePoints < nDim) {
          switch(nDim) {
            case 2: {
              for(MInt e = 0; e < noStlElementEdges; e++) {
                // If the edge has no cut, continue
                if(stl.edgeNotCutted[e]) continue;

                // If the edge has a cut, find the surface which is cutted by the edge
                // Line-line-intersection: if line cuts line, check if cutpoint is
                // inside of the cell surface
                for(MInt surface = 0; surface < 2 * nDim; surface++) {
                  const MFloat numer = stl.edges[e][1 + nDim] - stl.edges[e][1];
                  const MFloat denom = stl.edges[e][0 + nDim] - stl.edges[e][0];
                  const MFloat coordX = m_solver->c_coordinate(cellId, 0);
                  const MFloat coordY = m_solver->c_coordinate(cellId, 1);
                  const MFloat minStlEdgeX = mMin(stl.edges[e][0], stl.edges[e][0 + nDim]);
                  const MFloat minStlEdgeY = mMin(stl.edges[e][1], stl.edges[e][1 + nDim]);
                  const MFloat maxStlEdgeX = mMax(stl.edges[e][0], stl.edges[e][0 + nDim]);
                  const MFloat maxStlEdgeY = mMax(stl.edges[e][1], stl.edges[e][1 + nDim]);

                  // surface is parallel to x-axis
                  if(approx(normals[surface][0], F0, eps)) {
                    const MFloat y = coordY + normals[surface][1] * halfLength;
                    const MFloat x1 = coordX - halfLength;
                    const MFloat x2 = coordX + halfLength;

                    // STL element is parallel to x-axis.
                    // We can skip this case as these cut points are already considered in
                    // the other cases
                    if(approx(numer, F0, eps)) {
                      continue;
                    }
                    // stl element is parallel to y-axis
                    // check if there is a cut by comparing the edge coordinates
                    else if(approx(denom, F0, eps)) {
                      // Check if the x-coord of the stl edge lies in between the corner
                      // points of the surface
                      if(x1 <= stl.edges[e][0 + nDim] && x2 >= stl.edges[e][0 + nDim]) {
                        // Check if the y-coord of the surface lies in between the corner
                        // points of the stl edge
                        if(minStlEdgeY <= y && maxStlEdgeY >= y) {
                          std::vector<MFloat> cutPoint;
                          cutPoint.push_back(stl.edges[e][0 + nDim]);
                          cutPoint.push_back(y);
                          cutPointList.push_back(cutPoint);
                        }
                      }
                    }
                    // stl element is not parallel to an axis
                    // Insert the y-coord of the surface in the equation of the stl edge to
                    // calculate the x-coord of the cut.
                    // Check if the cut point lies in between the corner points of the surface.
                    else {
                      const MFloat cutPointX = (denom / numer) * (y - stl.edges[e][1]) + stl.edges[e][0];
                      const MFloat cutPointY = y;
                      if(cutPointX >= x1 && cutPointX <= x2 && cutPointX >= minStlEdgeX && cutPointX <= maxStlEdgeX) {
                        std::vector<MFloat> cutPoint;
                        cutPoint.push_back(cutPointX);
                        cutPoint.push_back(cutPointY);
                        cutPointList.push_back(cutPoint);
                      }
                    }
                  }
                  // surface is parallel to y-axis
                  if(approx(normals[surface][1], F0, eps)) {
                    const MFloat x = coordX + normals[surface][0] * halfLength;
                    const MFloat y1 = coordY - halfLength;
                    const MFloat y2 = coordY + halfLength;

                    // STL element is parallel to y-axis.
                    // We can skip this case as these cut points are already considered in
                    // the other cases
                    if(approx(denom, F0, eps)) {
                      continue;
                    }
                    // stl element is parallel to x-axis
                    // check if there is a cut by comparing the edge coordinates
                    else if(approx(numer, F0, eps)) {
                      // Check if the y-coord of the stl edge lies in between the corner
                      // points of the surface
                      if(y1 <= stl.edges[e][1 + nDim] && y2 >= stl.edges[e][1 + nDim]) {
                        // Check if the x-coord of the surface lies in between the corner
                        // points of the stl edge
                        if(minStlEdgeX <= x && maxStlEdgeX >= x) {
                          std::vector<MFloat> cutPoint;
                          cutPoint.push_back(x);
                          cutPoint.push_back(stl.edges[e][1 + nDim]);
                          cutPointList.push_back(cutPoint);
                        }
                      }
                    }
                    // stl element is not parallel to an axis
                    // Insert the x-coord of the surface in the equation of the stl edge to
                    // calculate the y-coord of the cut.
                    // Check if the cut point lies in between the corner points of the surface.
                    else {
                      const MFloat cutPointX = x;
                      const MFloat cutPointY = (numer / denom) * (x - stl.edges[e][0]) + stl.edges[e][1];

                      if(cutPointY >= y1 && cutPointY <= y2 && cutPointY >= minStlEdgeY && cutPointY <= maxStlEdgeY) {
                        std::vector<MFloat> cutPoint;
                        cutPoint.push_back(cutPointX);
                        cutPoint.push_back(cutPointY);
                        cutPointList.push_back(cutPoint);
                      }
                    }
                  }
                }
              }
              break;
            }
            case 3: {
              for(MInt e = 0; e < noStlElementEdges; e++) {
                // If the edge has no cut, continue
                if(stl.edgeNotCutted[e]) continue;
                // If the edge has a cut, find the surface which is cutted by the edge
                // Line-plane-intersection: if line cuts plane, check if cutpoint is inside of
                // the cell surface
                for(MInt surface = 0; surface < 2 * nDim; surface++) {
                  MFloat numer = F0;
                  MFloat denom = F0;
                  for(MInt d = 0; d < nDim; d++) {
                    MFloat locationVec = m_solver->c_coordinate(cellId, d) + normals[surface][d] * halfLength;
                    denom += (normals[surface][d] * (stl.edges[e][d + nDim] - stl.edges[e][d]));
                    numer += (normals[surface][d] * locationVec);
                  }
                  for(MInt d = 0; d < nDim; d++) {
                    numer -= (normals[surface][d] * stl.edges[e][d]);
                  }
                  if(!approx(denom, F0, eps)) {
                    // Check is cut point is inside cell surface
                    MFloat lambda = numer / denom;
                    if((lambda >= F0) && (lambda <= F1)) {
                      MBool validIntersection = true;
                      std::vector<MFloat> cutPoint;
                      for(MInt d = 0; d < nDim; d++) {
                        MFloat point = stl.edges[e][d] + lambda * (stl.edges[e][d + nDim] - stl.edges[e][d]);
                        if(fabs(point - target[d]) < m_solver->m_eps) {
                          point = target[d];
                        } else if(fabs(point - target[d + nDim]) < m_solver->m_eps) {
                          point = target[d + nDim];
                        }
                        cutPoint.push_back(point);
                        if(point < (target[d]) || point > (target[d + nDim])) {
                          validIntersection = false;
                        }
                      }
                      if(validIntersection) {
                        stl.cutSurfPerEdge[e][surface] = true;
                        cutPointList.push_back(cutPoint);
                        MBool surfCutBefore = false;
                        for(MInt prevE = 0; prevE < e - 1; prevE++) {
                          if(stl.cutSurfPerEdge[prevE][surface]) {
                            surfCutBefore = true;
                            break;
                          }
                        }
                        if(!surfCutBefore) stl.noCuttedSurfaces++;
                        stl.noCutsPerEdge[e]++;
                      }
                    }
                  }
                }
              }
              break;
            }
            default: {
              mTerm(1, AT_, "ERROR: Wrong number of dimensions! Aborting.");
              break;
            }
          }

          if(nDim == 3) {
            // e). Check if the edges of the element cut the stl element
            for(MInt e = 0; e < noElementEdges; e++) {
              ele.edgeOnSurface[e] = false;
              MFloat numer = F0;
              MFloat denom = F0;
              for(MInt d = 0; d < nDim; d++) {
                denom += (stl.normal[d] * (ele.edges[e][d + nDim] - ele.edges[e][d]));
                numer += (stl.normal[d] * stl.points[0][d]);
              }
              for(MInt d = 0; d < nDim; d++) {
                numer -= (stl.normal[d] * ele.edges[e][d]);
              }
              if(approx(denom, F0, eps)) {
                if(approx(numer, F0, eps)) {
                  // Element edge lies on the stl element surface
                  // We store the corner points of the edge as cut points for later
                  ele.edgeOnSurface[e] = true;
                  std::vector<MFloat> cutPoint1;
                  std::vector<MFloat> cutPoint2;
                  for(MInt d = 0; d < nDim; d++) {
                    cutPoint1.push_back(ele.edges[e][d]);
                    cutPoint2.push_back(ele.edges[e][d + nDim]);
                  }
                  MBool inside = pointInTriangle(stl.points[0], stl.points[1], stl.points[2], cutPoint1);
                  if(inside) {
                    cutPointList.push_back(cutPoint1);
                  }

                  inside = pointInTriangle(stl.points[0], stl.points[1], stl.points[2], cutPoint2);
                  if(inside) {
                    cutPointList.push_back(cutPoint2);
                  }
                  break;
                } else {
                  // Element edge is parallel to the stl element surface
                  continue;
                }
              } else {
                MFloat lambda = numer / denom;
                MBool inside = false;
                std::vector<MFloat> cutPoint;
                if((lambda >= F0) && (lambda <= F1)) {
                  for(MInt d = 0; d < nDim; d++) {
                    MFloat point = ele.edges[e][d] + lambda * (ele.edges[e][d + nDim] - ele.edges[e][d]);
                    cutPoint.push_back(point);
                  }
                  inside = pointInTriangle(stl.points[0], stl.points[1], stl.points[2], cutPoint);
                  if(inside) {
                    cutPointList.push_back(cutPoint);
                  }
                }
              }
            }
          }
        }

        for(MInt c1 = 0; c1 < (MInt)cutPointList.size(); c1++) {
          MInt c2 = c1 + 1;
          while(c2 < (MInt)cutPointList.size()) {
            MBool dublicated = true;
            for(MInt d = 0; d < nDim; d++) {
              if(!approx(cutPointList[c1][d], cutPointList[c2][d], eps)) {
                dublicated = false;
                break;
              }
            }
            if(dublicated) {
              cutPointList.erase(cutPointList.begin() + c2);
            } else {
              c2++;
            }
          }
        }

        if((MInt)cutPointList.size() < nDim) continue;
        setBoundaryStlElements(cutPointList, stl.normal.data(), segId, i);
      }
    }
  }
}

/** \brief Add cut faces to the collector of the boundary cell
 *  \author Moritz Waldmann
 *
 */
template <MInt nDim>
void FcBndryCnd<nDim>::setBoundaryStlElements(std::vector<std::vector<MFloat>> cutPointList, MFloat* normal, MInt segId,
                                              MInt i) {
  TRACE();

  std::vector<std::vector<MFloat>> newStlElementList;
  if(nDim == 2) {
    newStlElementList = createEdgesFromCutPoints(cutPointList);
  } else {
    newStlElementList = createTrianglesFromCutPoints(cutPointList, normal);
  }
  MFloat normalLength = F0;
  for(MInt d = 0; d < nDim; d++) {
    normalLength += (normal[d] * normal[d]);
  }
  normalLength = sqrt(normalLength);

  for(MInt t = 0; t < (MInt)newStlElementList.size(); t++) {
    std::vector<MFloat> vertices;
    std::vector<MFloat> elementNormal;
    for(MInt v = 0; v < nDim; v++) {
      for(MInt p = 0; p < nDim; p++) {
        vertices.push_back(newStlElementList[t][v * nDim + p]);
      }
      MFloat n = normal[v] / normalLength;
      elementNormal.push_back(n);
    }
    m_bndryCells[i].m_cutFaces.push_back(vertices);
    m_bndryCells[i].m_faceNormals.push_back(elementNormal);
    m_bndryCells[i].m_segmentIdOfCutFace.push_back(segId);
  }
}

/** \brief Write out the boundary segments in stl-file format.
 *  \author Moritz Waldmann
 *
 *  Debug function.
 */
template <MInt nDim>
void FcBndryCnd<nDim>::writeOutSTLBoundaries(MInt index, MString prefix) {
  TRACE();

  MInt segId = m_bndryCndSegIds[index];

  MString fileName = "Segment_" + std::to_string(segId) + prefix + ".stl";
  ofstream stlFile;
  stlFile.open(fileName);

  stlFile << "solid Visualization Toolkit generated SLA File\n";
  for(MInt i = m_bndryCndOffsets[index]; i < m_bndryCndOffsets[index + 1]; i++) {
    for(MInt t = 0; t < (MInt)m_bndryCells[i].m_cutFaces.size(); t++) {
      if(m_bndryCells[i].m_segmentIdOfCutFace[t] != segId) continue;
      stlFile << " facet normal ";
      for(MInt d = 0; d < nDim; d++) {
        if(d < (nDim - 1)) {
          stlFile << m_bndryCells[i].m_faceNormals[t][d] << " ";
        } else {
          stlFile << m_bndryCells[i].m_faceNormals[t][d] << "\n";
        }
      }
      stlFile << "  outer loop\n";
      for(MInt v = 0; v < nDim; v++) {
        stlFile << "   vertex ";
        for(MInt p = 0; p < nDim; p++) {
          if(p < (nDim - 1)) {
            stlFile << m_bndryCells[i].m_cutFaces[t][v * nDim + p] << " ";
          } else {
            stlFile << m_bndryCells[i].m_cutFaces[t][v * nDim + p] << "\n";
          }
        }
      }
      stlFile << "  endloop\n";
      stlFile << " endFacet\n";
    }
  }
  stlFile << "endsolid";
  stlFile.close();
}

/** \brief Calculation of the normal vector of the cut faces
 *  \author Moritz Waldmann
 *
 */
template <MInt nDim>
void FcBndryCnd<nDim>::calcAvgFaceNormal(MInt index) {
  TRACE();

  MInt segId = m_bndryCndSegIds[index];

  // Run through all corresponding boundary cells and change the node displacements
  for(MInt i = m_bndryCndOffsets[index]; i < m_bndryCndOffsets[index + 1]; i++) {
    MFloat avgNormal[nDim] = {F0};
    MInt cnt = 0;
    for(MInt t = 0; t < (MInt)m_bndryCells[i].m_cutFaces.size(); t++) {
      if(m_bndryCells[i].m_segmentIdOfCutFace[t] != segId) continue;
      for(MInt d = 0; d < nDim; d++) {
        avgNormal[d] += m_bndryCells[i].m_faceNormals[t][d];
      }
      cnt++;
    }
    MFloat normalLength = F0;
    for(MInt d = 0; d < nDim; d++) {
      m_bndryCells[i].m_avgFaceNormal[d] = avgNormal[d] / cnt;
      normalLength += (m_bndryCells[i].m_avgFaceNormal[d] * m_bndryCells[i].m_avgFaceNormal[d]);
    }
    normalLength = sqrt(normalLength);
    for(MInt d = 0; d < nDim; d++) {
      m_bndryCells[i].m_avgFaceNormal[d] = avgNormal[d] / normalLength;
    }
  }
}

/** \brief Moves a triangle depending on the simulated displacment inside a cell
 *  \author Moritz Waldmann
 *
 */
template <MInt nDim>
void FcBndryCnd<nDim>::updateTrianglePosition() {
  TRACE();

  for(MInt index = 0; index < (MInt)(m_bndryCndIds.size()); index++) {
    MInt segId = m_bndryCndSegIds[index];

    if(!m_solver->m_testRun) {
      MString prefix = "preTimeStep_" + std::to_string(globalTimeStep);
      writeOutSTLBoundaries(index, prefix);
    }
    for(MInt i = m_bndryCndOffsets[index]; i < m_bndryCndOffsets[index + 1]; i++) {
      for(MInt t = 0; t < (MInt)m_bndryCells[i].m_cutFaces.size(); t++) {
        if(m_bndryCells[i].m_segmentIdOfCutFace[t] != segId) continue;
        for(MInt d1 = 0; d1 < nDim; d1++) {
          MFloat x[nDim] = {F0};
          MFloat z[nDim] = {F0};
          for(MInt d2 = 0; d2 < nDim; d2++) {
            x[d2] = m_bndryCells[i].m_cutFaces[t][d1 * nDim + d2];
          }
          const MInt pCellId = m_bndryCells[i].m_cellId;
          m_solver->transformToLocal(pCellId, x, z);
          MFloatScratchSpace L_coef(nDim, m_solver->a_noNodes(pCellId) * nDim, AT_, "L_coef");
          L_coef.fill(F1);
          m_solver->getDisplacementInterpolationMatrix(pCellId, z, L_coef);

          for(MInt d2 = 0; d2 < nDim; d2++) {
            MFloat displacement = F0;
            for(MInt node = 0; node < m_solver->a_noNodes(pCellId); node++) {
              for(MInt d3 = 0; d3 < nDim; d3++) {
                MInt nodeId = m_solver->a_nodeIdsLocal(pCellId, node);
                displacement += L_coef(d2, node * nDim + d3) * m_solver->m_totalNodalDisplacements[nodeId * nDim + d3];
              }
            }
            m_bndryCells[i].m_cutFaces[t][d1 * nDim + d2] += displacement;
          }
        }
      }
    }
    if(!m_solver->m_testRun) {
      MString prefix = "postTimeStep_" + std::to_string(globalTimeStep);
      writeOutSTLBoundaries(index, prefix);
    }
  }
}

/** \brief Write out the boundary segments in stl-file format.
 *  \author Moritz Waldmann
 *
 *  Debug function.
 */
template <MInt nDim>
void FcBndryCnd<nDim>::writeOutModifiedBoundaries() {
  TRACE();

  for(MInt index = 0; index < (MInt)(m_bndryCndIds.size()); index++) {
    MInt segId = m_bndryCndSegIds[index];

    stringstream prefix;
    prefix << "_modified_iteration_" << globalTimeStep;
    MString fileName = "Segment_" + std::to_string(segId) + prefix.str() + ".stl";
    ofstream stlFile;
    stlFile.open(fileName);

    stlFile << "solid Visualization Toolkit generated SLA File\n";
    for(MInt i = m_bndryCndOffsets[index]; i < m_bndryCndOffsets[index + 1]; i++) {
      for(MInt t = 0; t < (MInt)m_bndryCells[i].m_cutFaces.size(); t++) {
        if(m_bndryCells[i].m_segmentIdOfCutFace[t] != segId) continue;
        stlFile << " facet normal ";
        for(MInt d = 0; d < nDim; d++) {
          if(d < (nDim - 1)) {
            stlFile << m_bndryCells[i].m_faceNormals[t][d] << " ";
          } else {
            stlFile << m_bndryCells[i].m_faceNormals[t][d] << "\n";
          }
        }
        MFloatScratchSpace points(nDim, nDim, AT_, "points");
        points.fill(F0);
        for(MInt d1 = 0; d1 < nDim; d1++) {
          std::array<MFloat, nDim> x{};
          std::array<MFloat, nDim> z{};
          for(MInt d2 = 0; d2 < nDim; d2++) {
            x[d2] = m_bndryCells[i].m_cutFaces[t][d1 * nDim + d2];
          }
          const MInt pCellId = m_bndryCells[i].m_cellId;
          m_solver->transformToLocal(pCellId, x.data(), z.data());
          MFloatScratchSpace L_coef(nDim, m_solver->a_noNodes(pCellId) * nDim, AT_, "L_coef");
          L_coef.fill(F1);
          m_solver->getDisplacementInterpolationMatrix(pCellId, z.data(), L_coef);

          for(MInt d2 = 0; d2 < nDim; d2++) {
            MFloat displacement = F0;
            for(MInt node = 0; node < m_solver->a_noNodes(pCellId); node++) {
              for(MInt d3 = 0; d3 < nDim; d3++) {
                MInt nodeId = m_solver->a_nodeIdsLocal(pCellId, node);
                displacement += L_coef(d2, node * nDim + d3) * m_solver->m_totalNodalDisplacements[nodeId * nDim + d3];
              }
            }
            points(d1, d2) = displacement + m_bndryCells[i].m_cutFaces[t][d1 * nDim + d2];
          }
        }
        stlFile << "  outer loop\n";
        for(MInt v = 0; v < nDim; v++) {
          stlFile << "   vertex ";
          for(MInt p = 0; p < nDim; p++) {
            if(p < (nDim - 1)) {
              stlFile << points(v, p) << " ";
            } else {
              stlFile << points(v, p) << "\n";
            }
          }
        }
        stlFile << "  endloop\n";
        stlFile << " endFacet\n";
      }
    }
    stlFile << "endsolid";
    stlFile.close();
  }
}

/** \brief Checks if a point is inside a triangle
 *  \author Moritz Waldmann
 *
 */
template <MInt nDim>
MBool FcBndryCnd<nDim>::pointInTriangle(MFloat* A, MFloat* B, MFloat* C, std::vector<MFloat> P) {
  TRACE();
  // Check if points are located inside the triangle
  MFloat u = F0;
  MFloat v = F0;
  MFloat dot[5] = {F0};
  // Check the start point
  for(MInt d = 0; d < nDim; d++) {
    dot[0] += (C[d] - A[d]) * (C[d] - A[d]);
    dot[1] += (C[d] - A[d]) * (B[d] - A[d]);
    dot[2] += (C[d] - A[d]) * (P[d] - A[d]);
    dot[3] += (B[d] - A[d]) * (B[d] - A[d]);
    dot[4] += (B[d] - A[d]) * (P[d] - A[d]);
  }
  u = (dot[3] * dot[2] - dot[1] * dot[4]) / (dot[0] * dot[3] - dot[1] * dot[1]);
  v = (dot[0] * dot[4] - dot[1] * dot[2]) / (dot[0] * dot[3] - dot[1] * dot[1]);
  if((u >= F0) && (v >= F0) && ((u + v) <= F1)) {
    return true;
  }
  return false;
}

/** \brief Refines the triangles of a segment by deviding them in the middle.
 *  \author Moritz Waldmann
 *
 */
template <MInt nDim>
void FcBndryCnd<nDim>::refineTriangle(MInt index) {
  TRACE();

  MInt segId = m_bndryCndSegIds[index];

  // Run through all corresponding boundary cells and change the node displacements
  for(MInt i = m_bndryCndOffsets[index]; i < m_bndryCndOffsets[index + 1]; i++) {
    MInt noCutFaces = (MInt)m_bndryCells[i].m_cutFaces.size();
    for(MInt t = 0; t < noCutFaces; t++) {
      if(m_bndryCells[i].m_segmentIdOfCutFace[t] != segId) continue;
      MFloat pointA[nDim];
      MFloat pointB[nDim];
      MFloat pointC[nDim];
      MFloat midPoint[nDim];
      std::vector<MFloat> triangleNormal;
      for(MInt d = 0; d < nDim; d++) {
        pointA[d] = m_bndryCells[i].m_cutFaces[t][0 * nDim + d];
        pointB[d] = m_bndryCells[i].m_cutFaces[t][1 * nDim + d];
        pointC[d] = m_bndryCells[i].m_cutFaces[t][2 * nDim + d];
        triangleNormal.push_back(m_bndryCells[i].m_faceNormals[t][d]);
      }
      MFloat lenAB = F0;
      MFloat lenAC = F0;
      MFloat lenBC = F0;
      for(MInt d = 0; d < nDim; d++) {
        lenAB += (pointA[d] - pointB[d]) * (pointA[d] - pointB[d]);
        lenAC += (pointA[d] - pointC[d]) * (pointA[d] - pointC[d]);
        lenBC += (pointB[d] * pointC[d]) * (pointB[d] * pointC[d]);
      }
      lenAB = sqrt(lenAB);
      lenAC = sqrt(lenAC);
      lenBC = sqrt(lenBC);
      std::vector<MFloat> newTriangle;
      if((lenAB > lenAC) && (lenAB > lenBC)) {
        for(MInt d = 0; d < nDim; d++) {
          midPoint[d] = F1B2 * (pointA[d] + pointB[d]);
          m_bndryCells[i].m_cutFaces[t][0 * nDim + d] = midPoint[d];
        }
        for(MInt d = 0; d < nDim; d++) {
          newTriangle.push_back(midPoint[d]);
        }
        for(MInt d = 0; d < nDim; d++) {
          newTriangle.push_back(pointA[d]);
        }
        for(MInt d = 0; d < nDim; d++) {
          newTriangle.push_back(pointC[d]);
        }
      } else if((lenAC > lenAB) && (lenAC > lenBC)) {
        for(MInt d = 0; d < nDim; d++) {
          midPoint[d] = F1B2 * (pointA[d] + pointC[d]);
          m_bndryCells[i].m_cutFaces[t][0 * nDim + d] = midPoint[d];
        }
        for(MInt d = 0; d < nDim; d++) {
          newTriangle.push_back(midPoint[d]);
        }
        for(MInt d = 0; d < nDim; d++) {
          newTriangle.push_back(pointA[d]);
        }
        for(MInt d = 0; d < nDim; d++) {
          newTriangle.push_back(pointB[d]);
        }
      } else if((lenBC > lenAB) && (lenBC > lenAC)) {
        for(MInt d = 0; d < nDim; d++) {
          midPoint[d] = F1B2 * (pointB[d] + pointC[d]);
          m_bndryCells[i].m_cutFaces[t][1 * nDim + d] = midPoint[d];
        }
        for(MInt d = 0; d < nDim; d++) {
          newTriangle.push_back(midPoint[d]);
        }
        for(MInt d = 0; d < nDim; d++) {
          newTriangle.push_back(pointA[d]);
        }
        for(MInt d = 0; d < nDim; d++) {
          newTriangle.push_back(pointB[d]);
        }
      }
      m_bndryCells[i].m_cutFaces.push_back(newTriangle);
      m_bndryCells[i].m_faceNormals.push_back(triangleNormal);
      m_bndryCells[i].m_segmentIdOfCutFace.push_back(segId);
    }
  }
}

/** \brief Returns a list of all elements cutting the target cell
 *  \author Moritz Waldmann
 *
 */
template <MInt nDim>
void FcBndryCnd<nDim>::getStlNodeList(const MInt cellId, std::vector<MInt>& nodeList) {
  TRACE();

  MFloat target[2 * nDim];
  MFloat cellHalfLength;

  // Define corners of current cell in target
  for(MInt d = 0; d < nDim; d++) {
    cellHalfLength = F1B2 * m_solver->c_cellLengthAtCell(cellId);
    target[d] = m_solver->c_coordinate(cellId, d) - cellHalfLength;
    target[d + nDim] = m_solver->c_coordinate(cellId, d) + cellHalfLength;
  }

  // Check for intersection with geometry elements
  m_solver->m_geometry->getIntersectionElements(target, nodeList, cellHalfLength, &m_solver->c_coordinate(cellId, 0));
}
// Explicit instantiations for 2D and 3D
template class FcBndryCnd<2>;
template class FcBndryCnd<3>;
