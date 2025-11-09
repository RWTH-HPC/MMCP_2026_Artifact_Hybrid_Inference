// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only


#include "iovtk.h"
#include "globals.h"

using namespace std;

//--------------------------------------------------------------------------
template <MInt nDim, class SysEqn>
VtkIo<nDim, SysEqn>::VtkIo(FvCartesianSolverXD<nDim, SysEqn>* solver) : m_solver(solver) {
  DEBUG("VtkIo::VtkIo: entry", MAIA_DEBUG_ALLOCATION);
  DEBUG("VtkIo::VtkIo: return", MAIA_DEBUG_ALLOCATION);
}

//--------------------------------------------------------------------------
template <MInt nDim, class SysEqn>
VtkIo<nDim, SysEqn>::~VtkIo() {
  DEBUG("VtkIo::~VtkIo: entry", MAIA_DEBUG_ALLOCATION);
  DEBUG("VtkIo::~VtkIo: return", MAIA_DEBUG_ALLOCATION);
}

//--------------------------------------------------------------------------
template <MInt nDim, class SysEqn>
MBool VtkIo<nDim, SysEqn>::initializeVtkXmlOutput(const MChar* fileName, MString outputDir_, MInt solutionInterval,
                                                  MBool restart, MBool firstUseInitializeVtkXmlOutput) {
  TRACE();
  if(domainId() != 0) { // only root process executes the following code
    return true;
  }

  m_log << "VtkIo::initializeVtkXmlOutput" << endl;

  MString outputFileName = "QOUT";
  MString pvd = ".pvd";

  MString pvdPath = outputDir_ + outputFileName + pvd; // e.g. out/QOUT.pvd
  MString pvdTmpPath = pvdPath + ".tmp";
  MString pvdAllPath = outputDir_ + outputFileName + "_all" + pvd;
  MString buPath = outputDir_ + outputFileName + "_BU" + pvd;
  MString vtuPath = MString("./") + outputFileName + "_";

  if(firstUseInitializeVtkXmlOutput) {
    if(fileExists(pvdPath.c_str())) {
      rename(pvdPath.c_str(), buPath.c_str());
    }
    if(fileExists(pvdTmpPath.c_str())) {
      remove(pvdTmpPath.c_str());
    }
    ofstream ofile(pvdTmpPath.c_str(), ios_base::out | ios_base::trunc);
    if(ofile.is_open() && ofile.good()) {
      ofile << R"(<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">)" << endl;
      ofile << "<Collection>" << endl;
      if(restart) {
        for(MInt t = 0; t <= globalTimeStep; t += solutionInterval) {
          for(MInt p = 0; p < noDomains(); p++) {
            stringstream tmp;
            tmp << outputDir_ << fileName << "_D" << p << "_00" << t << ".vtu";
            if(fileExists((tmp.str()).c_str())) {
              ofile << "<DataSet part=\"" << p << "\" timestep=\"" << t << "\" file=\"" << vtuPath << t << "_D" << p
                    << ".vtu\"/>\n";
            }
          }
        }
      }
      ofile.close();
      ofile.clear();
    } else {
      cerr << "Error opening file " << pvdTmpPath << endl;
    }
  }
  if(firstUseInitializeVtkXmlOutput) {
    ofstream ofile(pvdAllPath.c_str(), ios_base::out | ios_base::trunc);
    if(ofile.is_open() && ofile.good()) {
      ofile << R"(<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">)" << endl;
      ofile << "<Collection>" << endl;
      MInt t = 0;
      MBool found = true;
      while(found) {
        for(MInt p = 0; p < noDomains(); p++) {
          stringstream tmp;
          tmp << outputDir_ << fileName << "_D" << p << "_00" << t << ".vtu";
          if(fileExists(tmp.str())) {
            ofile << "<DataSet part=\"" << p << "\" timestep=\"" << t << "\" file=\"" << vtuPath << t << "_D" << p
                  << ".vtu\"/>\n";
          } else {
            found = false;
          }
        }
        t += solutionInterval;
      }
      ofile << "</Collection>" << endl;
      ofile << "</VTKFile>" << endl;
      ofile.close();
      ofile.clear();
    } else {
      cerr << "Error opening file " << pvdAllPath << endl;
    }
  }
  ofstream ofile(pvdTmpPath.c_str(), ios_base::out | ios_base::app);
  if(ofile.is_open() && ofile.good()) {
    for(MInt p = 0; p < noDomains(); p++) {
      ofile << "<DataSet part=\"" << p << "\" timestep=\"" << globalTimeStep << "\" file=\"" << vtuPath
            << globalTimeStep << "_D" << p << ".vtu\"/>\n";
    }
    ofile.close();
    ofile.clear();
  } else {
    cerr << "Error opening file " << pvdTmpPath << endl;
  }

  if(fileExists(pvdPath)) {
    remove(pvdPath.c_str());
  }
  copyFile(pvdTmpPath.c_str(), pvdPath);
  ofstream ofile2(pvdPath.c_str(), ios_base::out | ios_base::app);
  if(ofile2.is_open() && ofile2.good()) {
    ofile2 << "</Collection>" << endl;
    ofile2 << "</VTKFile>" << endl;
    ofile2.close();
  } else {
    cerr << "Error opening file " << pvdPath << endl;
  }

  return true;
}


template <MInt nDim, class SysEqn>
MInt VtkIo<nDim, SysEqn>::writeVtuArrayParallel(MPI_File& file, void* base, MPI_Offset offset, MPI_Offset count,
                                                MPI_Offset globalCount) {
  MPI_File_seek(file, offset, MPI_SEEK_CUR, AT_);
#ifdef VTU_OUTPUT_NONCOLLECTIVE
  MPI_Status status;
  MInt rcode = (count > 0) ? MPI_File_write(file, base, count, MPI_CHAR, &status) : MPI_SUCCESS;
#else
#ifndef VTU_OUTPUT_BLOCKING
  MPI_Status status;
  MInt rcode = MPI_File_write_all(file, base, count, MPI_CHAR, &status);
#else
  MPI_Request writeRequest = MPI_REQUEST_NULL;
  MInt rcode = MPI_File_iwrite_all(file, base, count, MPI_CHAR, &writeRequest);
  MPI_Wait(&writeRequest, MPI_STATUS_IGNORE, AT_);
#endif
#endif
  MPI_File_seek(file, globalCount - offset - count, MPI_SEEK_CUR, AT_);
  return rcode;
}


/** \brief Copy of parallel single-file VTU output using MPI I/O
 *
 *  Parameters:
 *   m_solver->m_vtuWritePointData: specifies whether the flow field is written as point or cell data
 *   m_solver->m_vtuCoordinatesThreshold: optionally specifies a bounding box to which the output domain is truncated
 *   m_solver->m_vtuLevelThreshold: optionally specifies the maximum cell level to be saved
 *
 *  Opening the .vtu file requires Paraview>=4.0, resp. VTK>=6
 *  If more than 4.2 billion points shall be written the function will rerun with 64 bit integer output
 *  For data compression compile with 'WITH_ZLIB' macro
 *
 *  For file-format reference see http://www.vtk.org/Wiki/VTK_XML_Formats
 *
 *  Do not attempt to change to data types used (e.g. float or uint64_t) to XYZ data types since this will break the
 * output!
 *
 * \author Lennart Scheiders
 * \date 01/2018
 */
/**
 * \brief write flow variables as XML VTK
 * \author Lennart Schneiders
 */

template <MInt nDim, class SysEqn>
void VtkIo<nDim, SysEqn>::writeVtkXmlOutput(const MChar* fileName) {
  IF_CONSTEXPR(nDim == 2) {
    TRACE();

    if(noDomains() > 1) {
      mTerm(1, AT_, "2D VTK output not yet in massive parallel mode.");
    }

    MInt noCells = 1;
    MInt noPoints = 1;

    MBool writePointData = false;

    const MInt baseCellType = 8;
    const MInt polyCellType = 7;

    MInt ghostCellId, nghbrId, counter;

    //-------------------------------
    const string dataType64 = "Float64";
#ifdef DOUBLE_PRECISION_OUTPUT
    const string dataType = "Float64";
#else
    const string dataType = "Float32";
#endif
    const string iDataType = "Int32";
    const string uIDataType = "UInt32";
    const string uI8DataType = "UInt8";
    // int inumber;
    MUint uinumber;
    //  unsigned char ui8number;
#ifdef DOUBLE_PRECISION_OUTPUT
    // MFloat fnumber;
    typedef MFloat flt;
#else
    // float fnumber;
    typedef float flt;
#endif
    //-------------------------------

    m_solver->m_firstUseInitializeVtkXmlOutput =
        initializeVtkXmlOutput(fileName, m_solver->outputDir(), m_solver->m_solutionInterval, m_solver->m_restart,
                               m_solver->m_firstUseInitializeVtkXmlOutput);

    if(domainId() == 0) {
      cerr << "extracted";
      cerr << "(" << m_solver->m_extractedCells->size() << ")...";
    }

    noPoints = m_solver->m_gridPoints->size();
    MInt noExtractedCells = m_solver->m_extractedCells->size();
    noCells = noExtractedCells;
    MInt scratchSize = m_solver->maxRefinementLevel() + 1;
    MFloatScratchSpace cellLength(scratchSize, AT_, "cellLength");

    for(MInt level = 0; level < m_solver->maxRefinementLevel() + 1; level++) {
      cellLength.p[level] = m_solver->c_cellLengthAtLevel(level);
    }

    MBool** pointInside = new MBool*[m_solver->m_bndryCells->size()];
    for(MInt i = 0; i < m_solver->m_bndryCells->size(); i++) {
      pointInside[i] = new MBool[IPOW2(nDim)];
      for(MInt j = 0; j < IPOW2(nDim); j++) {
        pointInside[i][j] = false;
      }
    }
    const MFloat signStencil[4][2] = {{-F1, -F1}, {F1, -F1}, {-F1, F1}, {F1, F1}};
    MFloat tmpPoint[nDim];

    for(MInt bndryId = 0; bndryId < m_solver->m_bndryCells->size(); bndryId++) {
      MInt cellId = m_solver->m_fvBndryCnd->m_solver->m_bndryCells->a[bndryId].m_cellId;
      for(MInt p = 0; p < IPOW2(nDim); p++) {
        for(MInt i = 0; i < nDim; i++) {
          tmpPoint[i] =
              m_solver->a_coordinate(cellId, i) + signStencil[p][i] * F1B2 * cellLength.p[m_solver->a_level(cellId)];
        }
        pointInside[bndryId][p] = m_solver->m_geometry->pointIsInside(tmpPoint);
      }
    }

    // 4 node ids for each direction
    const MInt ltable[4][3] = {{2, 3}, {0, 2}, {1, 0}, {3, 1}};

    const MInt pointOrder[4] = {0, 1, 3, 2};
    const MInt dirOrder[4] = {2, 1, 3, 0};
    const MInt reverseDir[4] = {1, 0, 3, 2};

    //  MFloat edgeIds[ 4 ];
    MInt centerId, cellId;

    ScratchSpace<MInt> cellTypes(noExtractedCells, AT_, "cellTypes");
    ScratchSpace<MInt> polyIds(noExtractedCells, AT_, "polyIds");
    ScratchSpace<MInt> reverseMapping(m_solver->a_noCells(), AT_, "reverseMapping");
    ScratchSpace<MInt> cutPointIds(m_noDirs * m_solver->m_fvBndryCnd->m_solver->m_bndryCells->size(), AT_,
                                   "cutPointIds");


    MInt noPolyCells;


    for(MInt i = 0; i < m_noDirs * m_solver->m_fvBndryCnd->m_solver->m_bndryCells->size(); i++) {
      cutPointIds.p[i] = -1;
    }

    if(domainId() == 0) {
      cerr << "prepare...";
    }
    for(MInt c = 0; c < m_solver->a_noCells(); c++) {
      reverseMapping.p[c] = -1;
    }
    for(MInt c = 0; c < noExtractedCells; c++) {
      reverseMapping.p[m_solver->m_extractedCells->a[c].m_cellId] = c;
    }


    // delete inactive grid points
    for(MInt p = 0; p < noPoints; p++) {
      MBool pointIsActive = true;
      for(MInt n = 0; n < m_solver->m_gridPoints->a[p].m_noAdjacentCells; n++) {
        cellId = m_solver->m_gridPoints->a[p].m_cellIds[n];
        if((m_solver->a_bndryId(cellId) > -1) && (reverseMapping.p[cellId] > -1)) {
          for(MInt k = 0; k < IPOW2(nDim); k++) {
            if(m_solver->m_extractedCells->a[reverseMapping.p[cellId]].m_pointIds[k] == p) {
              if(pointInside[m_solver->a_bndryId(cellId)][k]) {
                pointIsActive = false;
                break;
              }
            }
          }
        }
      }
      if(!pointIsActive) {
        noPoints--;

        for(MInt n = 0; n < m_solver->m_gridPoints->a[p].m_noAdjacentCells; n++) {
          cellId = m_solver->m_gridPoints->a[p].m_cellIds[n];
          if(reverseMapping.p[cellId] > -1) {
            for(MInt k = 0; k < IPOW2(nDim); k++) {
              if(m_solver->m_extractedCells->a[reverseMapping.p[cellId]].m_pointIds[k] == p) {
                m_solver->m_extractedCells->a[reverseMapping.p[cellId]].m_pointIds[k] = -1;
              }
            }
          }
        }

        for(MInt n = 0; n < m_solver->m_gridPoints->a[noPoints].m_noAdjacentCells; n++) {
          cellId = m_solver->m_gridPoints->a[noPoints].m_cellIds[n];
          if(reverseMapping.p[cellId] > -1) {
            for(MInt k = 0; k < IPOW2(nDim); k++) {
              if(m_solver->m_extractedCells->a[reverseMapping.p[cellId]].m_pointIds[k] == noPoints) {
                m_solver->m_extractedCells->a[reverseMapping.p[cellId]].m_pointIds[k] = p;
              }
            }
          }
          m_solver->m_gridPoints->a[p].m_cellIds[n] = m_solver->m_gridPoints->a[noPoints].m_cellIds[n];
          m_solver->m_gridPoints->a[p].m_noAdjacentCells = m_solver->m_gridPoints->a[noPoints].m_noAdjacentCells;
        }
        for(MInt i = 0; i < nDim; i++) {
          m_solver->m_gridPoints->a[p].m_coordinates[i] = m_solver->m_gridPoints->a[noPoints].m_coordinates[i];
        }
        m_solver->m_gridPoints->resetSize(noPoints);
        p--;
      }
    }

    //---

    MBool isPolyCell;
    noPolyCells = 0;

    for(MInt c = 0; c < noExtractedCells; c++) {
      cellId = m_solver->m_extractedCells->a[c].m_cellId;
      cellTypes.p[c] = baseCellType;
      polyIds.p[c] = -1;
      isPolyCell = false;
      if(m_solver->a_bndryId(cellId) > -1) {
        isPolyCell = true;
      } else {
        for(MInt i = 0; i < m_noDirs; i++) {
          if(m_solver->a_hasNeighbor(cellId, i) > 0) {
            if(m_solver->c_noChildren(m_solver->c_neighborId(cellId, i)) > 0) {
              isPolyCell = true;
            }
          }
        }
      }
      if(isPolyCell) {
        cellTypes.p[c] = polyCellType;
        polyIds.p[c] = noPolyCells;
        noPolyCells++;
      }
    }

    MBool isPolyFace;
    // MBool pointExists;
    MInt eCell;
    const MInt maxNoExtraPoints = 8;
    const MInt maxNoPolyCells = noPolyCells;
    ScratchSpace<MInt> polyCellLink(maxNoPolyCells, AT_, "polyCellLink");
    ScratchSpace<MInt> extraPoints(maxNoPolyCells * maxNoExtraPoints, AT_, "extraPoints");
    ScratchSpace<MInt> noExtraPoints(maxNoPolyCells, AT_, "noExtraPoints");
    ScratchSpace<MInt> noInternalPoints(maxNoPolyCells, AT_, "noInternalPoints");
    MInt offset, pointCount; // faceCount;

    for(MInt c = 0; c < noExtractedCells; c++) {
      if(cellTypes.p[c] == polyCellType) {
        polyCellLink.p[polyIds.p[c]] = c;
      }
    }


    const MInt noInternalGridPoints = m_solver->m_gridPoints->size();

    for(MInt pc = 0; pc < noPolyCells; pc++) {
      eCell = polyCellLink.p[pc];
      cellId = m_solver->m_extractedCells->a[eCell].m_cellId;
      noInternalPoints.p[pc] = IPOW2(nDim);
      noExtraPoints.p[pc] = 0;

      // a) cut cells at the boundary
      if(m_solver->a_bndryId(cellId) > -1) {
        MInt bndryId = m_solver->a_bndryId(cellId);

        noInternalPoints.p[pc] = 0;

        for(MInt i = 0; i < m_noDirs; i++) {
          MInt p = pointOrder[i];
          MInt dir = dirOrder[i];

          // 0. determine internal vertices
          if(!pointInside[bndryId][p]) {
            extraPoints.p[pc * maxNoExtraPoints + noExtraPoints.p[pc]] =
                m_solver->m_extractedCells->a[eCell].m_pointIds[p];
            noExtraPoints.p[pc]++;
          }

          // 1. add cut points
          for(MInt cp = 0; cp < m_solver->m_fvBndryCnd->m_solver->m_bndryCells->a[bndryId].m_srfcs[0]->m_noCutPoints;
              cp++) {
            if(m_solver->m_fvBndryCnd->m_solver->m_bndryCells->a[bndryId].m_srfcs[0]->m_cutEdge[cp] == dir) {
              MInt gridPointId = -1;
              if(m_solver->a_hasNeighbor(cellId, dir) > 0) {
                if(m_solver->a_bndryId(m_solver->c_neighborId(cellId, dir)) > -1) {
                  gridPointId =
                      cutPointIds
                          .p[m_noDirs * m_solver->a_bndryId(m_solver->c_neighborId(cellId, dir)) + reverseDir[dir]];
                  if(gridPointId > -1) {
                    m_solver->m_gridPoints->a[gridPointId]
                        .m_cellIds[m_solver->m_gridPoints->a[gridPointId].m_noAdjacentCells] = cellId;
                    m_solver->m_gridPoints->a[gridPointId].m_noAdjacentCells++;
                  }
                }
              }
              if(gridPointId >= noPoints) {
                mTerm(1, AT_, "Error 0 in FvMbSolver3D::writeVtkXmlOutput(..). Quit.");
              }
              if(gridPointId < 0) {
                gridPointId = m_solver->m_gridPoints->size();
                m_solver->m_gridPoints->append();
                noPoints++;
                m_solver->m_gridPoints->a[gridPointId].m_noAdjacentCells = 1;
                m_solver->m_gridPoints->a[gridPointId].m_cellIds[0] = cellId;
                for(MInt j = 0; j < nDim; j++) {
                  m_solver->m_gridPoints->a[gridPointId].m_coordinates[j] =
                      m_solver->m_fvBndryCnd->m_solver->m_bndryCells->a[m_solver->a_bndryId(cellId)]
                          .m_srfcs[0]
                          ->m_cutCoordinates[cp][j];
                }
              }
              cutPointIds.p[m_noDirs * bndryId + dir] = gridPointId;

              extraPoints.p[pc * maxNoExtraPoints + noExtraPoints.p[pc]] = gridPointId;
              noExtraPoints.p[pc]++;
            }
          }
        }
      }
      // b) cells at fine/coarse mesh interfaces
      else {
        noInternalPoints.p[pc] = 0;

        for(MInt i = 0; i < m_noDirs; i++) {
          MInt p = pointOrder[i];
          MInt dir = dirOrder[i];

          extraPoints.p[pc * maxNoExtraPoints + noExtraPoints.p[pc]] =
              m_solver->m_extractedCells->a[eCell].m_pointIds[p];
          noExtraPoints.p[pc]++;

          isPolyFace = false;
          if(m_solver->a_hasNeighbor(cellId, dir) > 0) {
            if(m_solver->c_noChildren(m_solver->c_neighborId(cellId, dir)) > 0) {
              isPolyFace = true;
            }
          }

          if(isPolyFace) {
            // store four faces (polygon numbering!)
            nghbrId = reverseMapping.p[m_solver->c_childId(m_solver->c_neighborId(cellId, dir), ltable[i][0])];
            if(nghbrId < 0) {
              cerr << endl
                   << cellId << " " << m_solver->c_neighborId(cellId, dir) << " "
                   << m_solver->c_childId(m_solver->c_neighborId(cellId, dir), ltable[i][0]) << " / "
                   << m_solver->a_level(cellId) << " " << m_solver->a_level(m_solver->c_neighborId(cellId, dir)) << " "
                   << m_solver->a_level(m_solver->c_childId(m_solver->c_neighborId(cellId, dir), ltable[i][0])) << " / "
                   << m_solver->c_noChildren(cellId) << " "
                   << m_solver->c_noChildren(m_solver->c_neighborId(cellId, dir)) << " "
                   << m_solver->c_noChildren(m_solver->c_childId(m_solver->c_neighborId(cellId, dir), ltable[i][0]))
                   << endl;
              cerr << m_solver->c_childId(cellId, 0) << " " << m_solver->c_childId(cellId, 1) << " "
                   << m_solver->c_childId(cellId, 2) << " " << m_solver->c_childId(cellId, 3) << endl;
              mTerm(1, AT_, "Error 1 in FvMbSolver2D::writeVtkXmlOutput(..). Quit.");
            }
            centerId = m_solver->m_extractedCells->a[nghbrId].m_pointIds[ltable[i][1]];
            extraPoints.p[pc * maxNoExtraPoints + noExtraPoints.p[pc]] = centerId;
            noExtraPoints.p[pc]++;

            if(m_solver->m_gridPoints->a[centerId].m_noAdjacentCells < IPOW2(nDim)) {
              m_solver->m_gridPoints->a[centerId].m_cellIds[m_solver->m_gridPoints->a[centerId].m_noAdjacentCells] =
                  cellId;
              m_solver->m_gridPoints->a[centerId].m_noAdjacentCells++;
            } else {
              mTerm(1, AT_, "Error 2 in FvMbSolver2D::writeVtkXmlOutput(..). Quit.");
            }
          }
        }
      }
    }

    pointCount = 0;
    for(MInt c = 0; c < noExtractedCells; c++) {
      MInt pc = polyIds.p[c];
      if(pc > -1) {
        pointCount += noInternalPoints.p[pc];
        pointCount += noExtraPoints.p[pc];
      } else {
        pointCount += IPOW2(nDim);
      }
    }


    MFloat gradU[2][2];
    for(MInt c = 0; c < m_solver->m_noActiveCells; c++) {
      cellId = m_solver->m_activeCellIds[c];
      for(MInt i = 0; i < nDim; i++) {
        for(MInt j = 0; j < nDim; j++) {
          gradU[i][j] =
              m_solver->a_slope(cellId, sysEqn().PV->VV[i], j) * m_solver->m_referenceLength / m_solver->m_UInfinity;
        }
      }
      MFloat omega = F1B2 * (gradU[1][0] - gradU[0][1]);
      MFloat S = 0;
      for(MInt i = 0; i < nDim; i++) {
        for(MInt j = 0; j < nDim; j++) {
          S += POW2(F1B2 * (gradU[i][j] + gradU[j][i]));
        }
      }
      MFloat Q = POW2(omega) - F1B2 * POW2(S);
      m_solver->a_slope(cellId, 2, 0) = Q;
      m_solver->a_slope(cellId, 0, 0) = omega;
    }


    //=============================================
    //================ GATHER =====================
    //=============================================

    if(domainId() == 0) {
      cerr << "gather...";
    }

    //-----------
    // points
    ScratchSpace<flt> points(3 * noPoints, AT_, "points");
    for(MInt p = 0; p < noPoints; p++) {
      for(MInt i = 0; i < nDim; i++) {
        points.p[3 * p + i] = (flt)m_solver->m_gridPoints->a[p].m_coordinates[i];
      }
      points.p[3 * p + 2] = (flt)F0;
    }

    //-----------
    // connectivity
    ScratchSpace<MUint> connectivity(pointCount, AT_, "connectivity");
    counter = 0;
    for(MInt c = 0; c < noExtractedCells; c++) {
      MInt pc = polyIds.p[c];
      if(pc > -1) {
        for(MInt p = 0; p < noExtraPoints.p[pc]; p++) {
          connectivity.p[counter] = (MUint)extraPoints.p[pc * maxNoExtraPoints + p];
          counter++;
        }
      } else {
        for(MInt p = 0; p < IPOW2(nDim); p++) {
          connectivity.p[counter] = (MUint)m_solver->m_extractedCells->a[c].m_pointIds[p];
          counter++;
        }
      }
    }
    if(counter != pointCount) {
      mTerm(1, AT_, "E1");
    }

    //-----------
    // offsets
    ScratchSpace<MUint> offsets(noCells, AT_, "offsets");
    counter = 0;
    for(MInt c = 0; c < noExtractedCells; c++) {
      MInt pc = polyIds.p[c];
      if(pc > -1) {
        counter += noInternalPoints.p[pc] + noExtraPoints.p[pc];
      } else {
        counter += IPOW2(nDim);
      }

      offsets.p[c] = (MUint)counter;
    }


    //-----------
    // types
    ScratchSpace<unsigned char> types(noCells, AT_, "types");
    for(MInt c = 0; c < noExtractedCells; c++) {
      types.p[c] = (unsigned char)cellTypes.p[c];
    }

    //-----------
    // cellIds
    ScratchSpace<MInt> cellIds(noCells, AT_, "cellIds");
    for(MInt c = 0; c < noExtractedCells; c++) {
      cellIds.p[c] = m_solver->m_extractedCells->a[c].m_cellId;
    }

    //-----------
    // bndryIds
    ScratchSpace<MInt> bndryIds(noCells, AT_, "bndryIds");
    for(MInt c = 0; c < noExtractedCells; c++) {
      bndryIds.p[c] = m_solver->a_bndryId(m_solver->m_extractedCells->a[c].m_cellId);
    }
    /*
    //-----------
    // drho_dt
    ScratchSpace<MFloat> drho_dt(noCells, AT_, "drho_dt" );
    for( MInt c = 0; c < noExtractedCells; c++ ) {
      drho_dt.p[ c ] = m_rightHandSide[ m_solver->m_extractedCells->a[ c ].m_cellId ][ CV->RHO ];
    }
    */
    //-----------
    // levelSetFunction
    ScratchSpace<flt> levelSetFunction((m_solver->m_levelSet) ? noCells : 0, AT_, "levelSetFunction");
    if(m_solver->m_levelSet) {
      for(MInt c = 0; c < noExtractedCells; c++) {
        if(m_solver->m_combustion) {
          levelSetFunction.p[c] = m_solver->a_levelSetFunction(m_solver->m_extractedCells->a[c].m_cellId, 0);
        } else if(!m_solver->m_levelSetMb) {
          levelSetFunction.p[c] = m_solver->a_levelSetFunction(m_solver->m_extractedCells->a[c].m_cellId, 0);
        }
      }
    }

    //-----------
    // additional cellData/pointData
    MInt asize = noCells;
    if(writePointData) {
      asize = noPoints;
    }
    ScratchSpace<flt> pressure(asize, AT_, "pressure");
    ScratchSpace<flt> velocity(3 * asize, AT_, "velocity");
    ScratchSpace<flt> density(asize, AT_, "density");
    ScratchSpace<flt> vorticity(asize, AT_, "vorticity");
    ScratchSpace<flt> progressVariable((m_solver->m_combustion) ? asize : 0, AT_, "progressVariable");
    //  ScratchSpace<flt> Qcriterion(asize, AT_, "Qcriterion" );


    MFloat tmp;
    MFloat weights[/*IPOW2(nDim)*/ 2 * 2];
    if(writePointData) {
      for(MInt p = 0; p < noInternalGridPoints; p++) {
        for(MInt n = 0; n < m_solver->m_gridPoints->a[p].m_noAdjacentCells; n++) {
          weights[n] = F0;
        }
        tmp = F0;
        for(MInt n = 0; n < m_solver->m_gridPoints->a[p].m_noAdjacentCells; n++) {
          cellId = m_solver->m_gridPoints->a[p].m_cellIds[n];
          if(m_solver->a_bndryId(cellId) > -1) {
            weights[n] = F1
                         / (sqrt(POW2(m_solver->a_coordinate(cellId, 0)
                                      + m_solver->m_fvBndryCnd->m_solver->m_bndryCells->a[m_solver->a_bndryId(cellId)]
                                            .m_coordinates[0]
                                      - m_solver->m_gridPoints->a[p].m_coordinates[0])
                                 + POW2(m_solver->a_coordinate(cellId, 1)
                                        + m_solver->m_fvBndryCnd->m_solver->m_bndryCells->a[m_solver->a_bndryId(cellId)]
                                              .m_coordinates[1]
                                        - m_solver->m_gridPoints->a[p].m_coordinates[1])));
            tmp += weights[n];
          } else {
            weights[n] =
                F1
                / (sqrt(POW2(m_solver->a_coordinate(cellId, 0) - m_solver->m_gridPoints->a[p].m_coordinates[0])
                        + POW2(m_solver->a_coordinate(cellId, 1) - m_solver->m_gridPoints->a[p].m_coordinates[1])));
            tmp += weights[n];
          }
        }
        for(MInt n = 0; n < m_solver->m_gridPoints->a[p].m_noAdjacentCells; n++) {
          weights[n] /= tmp;
        }
        //---
        pressure.p[p] = (flt)F0;
        density.p[p] = (flt)F0;
        vorticity.p[p] = (flt)F0;
        // Qcriterion.p[ p ] = (flt)F0 ;
        if(m_solver->m_combustion) {
          progressVariable.p[p] = (flt)F0;
        }
        for(MInt i = 0; i < nDim; i++) {
          velocity.p[3 * p + i] = (flt)F0;
        }
        velocity.p[3 * p + 2] = (flt)F0;

        for(MInt n = 0; n < m_solver->m_gridPoints->a[p].m_noAdjacentCells; n++) {
          cellId = m_solver->m_gridPoints->a[p].m_cellIds[n];
          pressure.p[p] += (flt)weights[n] * m_solver->a_pvariable(cellId, sysEqn().PV->P);
          density.p[p] += (flt)weights[n] * m_solver->a_pvariable(cellId, sysEqn().PV->RHO);
          vorticity.p[p] += (flt)weights[n] * m_solver->a_slope(cellId, 0, 0);
          IF_CONSTEXPR(hasPV_C<SysEqn>::value) {
            if(m_solver->m_combustion) {
              progressVariable.p[p] += (flt)weights[n] * m_solver->a_pvariable(cellId, sysEqn().PV->C);
            }
          }
          for(MInt i = 0; i < nDim; i++) {
            velocity.p[3 * p + i] += (flt)weights[n] * m_solver->a_pvariable(cellId, sysEqn().PV->VV[i]);
          }
        }
      }
      for(MInt p = noInternalGridPoints; p < noPoints; p++) {
        for(MInt n = 0; n < m_solver->m_gridPoints->a[p].m_noAdjacentCells; n++) {
          weights[n] = F0;
        }
        tmp = F0;
        for(MInt n = 0; n < m_solver->m_gridPoints->a[p].m_noAdjacentCells; n++) {
          cellId = m_solver->m_gridPoints->a[p].m_cellIds[n];
          if(m_solver->a_bndryId(cellId) < 0) {
            continue;
          }
          weights[n] = F1
                       / (sqrt(POW2(m_solver->m_fvBndryCnd->m_solver->m_bndryCells->a[m_solver->a_bndryId(cellId)]
                                        .m_srfcs[0]
                                        ->m_coordinates[0]
                                    - m_solver->m_gridPoints->a[p].m_coordinates[0])
                               + POW2(m_solver->m_fvBndryCnd->m_solver->m_bndryCells->a[m_solver->a_bndryId(cellId)]
                                          .m_srfcs[0]
                                          ->m_coordinates[1]
                                      - m_solver->m_gridPoints->a[p].m_coordinates[1])));
          tmp += weights[n];
        }
        for(MInt n = 0; n < m_solver->m_gridPoints->a[p].m_noAdjacentCells; n++) {
          weights[n] /= tmp;
        }
        pressure.p[p] = (flt)F0;
        density.p[p] = (flt)F0;
        vorticity.p[p] = (flt)F0;
        //      Qcriterion.p[ p ] = (flt)F0 ;
        if(m_solver->m_combustion) {
          progressVariable.p[p] = (flt)F0;
        }
        for(MInt i = 0; i < nDim; i++) {
          velocity.p[3 * p + i] = (flt)F0;
        }
        velocity.p[3 * p + 2] = F0;
        for(MInt n = 0; n < m_solver->m_gridPoints->a[p].m_noAdjacentCells; n++) {
          cellId = m_solver->m_gridPoints->a[p].m_cellIds[n];
          if(m_solver->a_bndryId(cellId) < 0) {
            continue;
          }
          ghostCellId = m_solver->m_fvBndryCnd->m_solver->m_bndryCells->a[m_solver->a_bndryId(cellId)]
                            .m_srfcVariables[0]
                            ->m_ghostCellId;
          pressure.p[p] +=
              (flt)weights[n] * F1B2
              * (m_solver->a_pvariable(cellId, sysEqn().PV->P) + m_solver->a_pvariable(ghostCellId, sysEqn().PV->P));
          density.p[p] += (flt)weights[n] * F1B2
                          * (m_solver->a_pvariable(cellId, sysEqn().PV->RHO)
                             + m_solver->a_pvariable(ghostCellId, sysEqn().PV->RHO));
          vorticity.p[p] += (flt)weights[n] * m_solver->a_slope(cellId, 0, 0);
          IF_CONSTEXPR(hasPV_C<SysEqn>::value) {
            if(m_solver->m_combustion) {
              progressVariable.p[p] += (flt)weights[n] * F1B2
                                       * (m_solver->a_pvariable(cellId, sysEqn().PV->C)
                                          + m_solver->a_pvariable(ghostCellId, sysEqn().PV->C));
            }
          }
          for(MInt i = 0; i < nDim; i++) {
            velocity.p[3 * p + i] += (flt)weights[n] * F1B2
                                     * (m_solver->a_pvariable(cellId, sysEqn().PV->VV[i])
                                        + m_solver->a_pvariable(ghostCellId, sysEqn().PV->VV[i]));
          }
        }
      }
    } else {
      for(MInt c = 0; c < noExtractedCells; c++) {
        cellId = m_solver->m_extractedCells->a[c].m_cellId;
        pressure.p[c] = (flt)m_solver->a_pvariable(cellId, sysEqn().PV->P);
        density.p[c] = (flt)m_solver->a_pvariable(cellId, sysEqn().PV->RHO);
        vorticity.p[c] = (flt)m_solver->a_slope(cellId, 0, 0);
        IF_CONSTEXPR(hasPV_C<SysEqn>::value) {
          if(m_solver->m_combustion) {
            progressVariable.p[c] = (flt)m_solver->a_pvariable(cellId, sysEqn().PV->C);
          }
        }
        for(MInt i = 0; i < nDim; i++) {
          velocity.p[3 * c + i] = (flt)m_solver->a_pvariable(cellId, sysEqn().PV->VV[i]);
        }
        velocity.p[3 * c + 2] = (flt)F0;
      }
    }


    if(domainId() == 0) {
      cerr << "write...";
    }

    ofstream ofl;
    ofl.open(fileName);

    if(ofl) {
      offset = 0;

      //================== VTKFile =================
      ofl << "<?xml version=\"1.0\"?>" << endl;
      ofl << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">" << endl;
      ofl << "<UnstructuredGrid>" << endl;

      //================== FieldData =================
      ofl << "<DataArray type=\"" << dataType
          << "\" Name=\"timeStep\" format=\"ascii\" NumberOfTuples=\"1\" RangeMin=\"" << m_solver->timeStep()
          << "\" RangeMax=\"" << m_solver->timeStep() << "\">" << endl;
      ofl << m_solver->timeStep() << endl;
      ofl << "</DataArray>" << endl;
      ofl << "<FieldData>" << endl;
      ofl << "<DataArray type=\"" << dataType << "\" Name=\"refTime\" format=\"ascii\" NumberOfTuples=\"1\" RangeMin=\""
          << m_solver->m_timeRef << "\" RangeMax=\"" << m_solver->m_timeRef << "\">" << endl;
      ofl << m_solver->m_timeRef << endl;
      ofl << "</DataArray>" << endl;
      if(m_solver->m_outputPhysicalTime) {
        ofl << "<DataArray type=\"" << dataType << "\" Name=\"time\" format=\"ascii\" NumberOfTuples=\"1\" RangeMin=\""
            << m_solver->m_physicalTime << "\" RangeMax=\"" << m_solver->m_physicalTime << "\">" << endl;
        ofl << m_solver->m_physicalTime << endl;
        ofl << "</DataArray>" << endl;
        ofl << "<DataArray type=\"" << dataType
            << "\" Name=\"internalTime\" format=\"ascii\" NumberOfTuples=\"1\" RangeMin=\"" << m_solver->m_time
            << "\" RangeMax=\"" << m_solver->m_time << "\">" << endl;
        ofl << m_solver->m_time << endl;
      } else {
        ofl << "<DataArray type=\"" << dataType << "\" Name=\"time\" format=\"ascii\" NumberOfTuples=\"1\" RangeMin=\""
            << m_solver->m_time << "\" RangeMax=\"" << m_solver->m_time << "\">" << endl;
        ofl << m_solver->m_time << endl;
        ofl << "</DataArray>" << endl;
        ofl << "<DataArray type=\"" << dataType
            << "\" Name=\"physicalTime\" format=\"ascii\" NumberOfTuples=\"1\" RangeMin=\"" << m_solver->m_physicalTime
            << "\" RangeMax=\"" << m_solver->m_physicalTime << "\">" << endl;
        ofl << m_solver->m_physicalTime << endl;
      }
      ofl << "</DataArray>" << endl;
      ofl << "<DataArray type=\"" << dataType << "\" Name=\"Ma\" format=\"ascii\" NumberOfTuples=\"1\" RangeMin=\""
          << m_solver->m_Ma << "\" RangeMax=\"" << m_solver->m_Ma << "\">" << endl;
      ofl << m_solver->m_Ma << endl;
      ofl << "</DataArray>" << endl;
      ofl << "<DataArray type=\"" << dataType << "\" Name=\"Re\" format=\"ascii\" NumberOfTuples=\"1\" RangeMin=\""
          << m_solver->m_Re << "\" RangeMax=\"" << m_solver->m_Re << "\">" << endl;
      ofl << m_solver->m_Re << endl;
      ofl << "</DataArray>" << endl;
      ofl << "<DataArray type=\"" << dataType << "\" Name=\"Pr\" format=\"ascii\" NumberOfTuples=\"1\" RangeMin=\""
          << m_solver->m_Pr << "\" RangeMax=\"" << m_solver->m_Pr << "\">" << endl;
      ofl << m_solver->m_Pr << endl;
      ofl << "</DataArray>" << endl;
      ofl << "<DataArray type=\"" << dataType << "\" Name=\"gamma\" format=\"ascii\" NumberOfTuples=\"1\" RangeMin=\""
          << m_solver->sysEqn().gamma_Ref() << "\" RangeMax=\"" << m_solver->sysEqn().gamma_Ref() << "\">" << endl;
      ofl << m_solver->sysEqn().gamma_Ref() << endl;
      ofl << "</DataArray>" << endl;
      ofl << "<DataArray type=\"" << dataType << "\" Name=\"CFL\" format=\"ascii\" NumberOfTuples=\"1\" RangeMin=\""
          << m_solver->m_cfl << "\" RangeMax=\"" << m_solver->m_cfl << "\">" << endl;
      ofl << m_solver->m_cfl << endl;
      ofl << "</DataArray>" << endl;
      ofl << "<DataArray type=\"" << dataType << "\" Name=\"uInf\" format=\"ascii\" NumberOfTuples=\"1\" RangeMin=\""
          << m_solver->m_UInfinity << "\" RangeMax=\"" << m_solver->m_UInfinity << "\">" << endl;
      ofl << m_solver->m_UInfinity << endl;
      ofl << "</DataArray>" << endl;
      ofl << "<DataArray type=\"" << dataType << "\" Name=\"vInf\" format=\"ascii\" NumberOfTuples=\"1\" RangeMin=\""
          << m_solver->m_VInfinity << "\" RangeMax=\"" << m_solver->m_VInfinity << "\">" << endl;
      ofl << m_solver->m_VInfinity << endl;
      ofl << "</DataArray>" << endl;
      /*    ofl << "<DataArray type=\"" << dataType << "\" Name=\"wInf\" format=\"ascii\"
      NumberOfTuples=\"1\" RangeMin=\"" << m_solver->m_WInfinity << "\" RangeMax=\"" << m_solver->m_WInfinity << "\">"
      << endl; ofl << m_solver->m_WInfinity << endl; ofl << "</DataArray>" << endl;*/
      ofl << "<DataArray type=\"" << dataType << "\" Name=\"pInf\" format=\"ascii\" NumberOfTuples=\"1\" RangeMin=\""
          << m_solver->m_PInfinity << "\" RangeMax=\"" << m_solver->m_PInfinity << "\">" << endl;
      ofl << m_solver->m_PInfinity << endl;
      ofl << "</DataArray>" << endl;
      ofl << "<DataArray type=\"" << dataType << "\" Name=\"rhoInf\" format=\"ascii\" NumberOfTuples=\"1\" RangeMin=\""
          << m_solver->m_rhoInfinity << "\" RangeMax=\"" << m_solver->m_rhoInfinity << "\">" << endl;
      ofl << m_solver->m_rhoInfinity << endl;
      ofl << "</DataArray>" << endl;
      ofl << "</FieldData>" << endl;
      if(m_solver->m_combustion) {
        ofl << "<DataArray type=\"" << dataType
            << "\" Name=\"flSpeed\" format=\"ascii\" NumberOfTuples=\"1\" RangeMin=\"" << m_solver->m_flameSpeed
            << "\" RangeMax=\"" << m_solver->m_flameSpeed << "\">" << endl;
        ofl << m_solver->m_flameSpeed << endl;
        ofl << "</DataArray>" << endl;
        ofl << "<DataArray type=\"" << dataType
            << "\" Name=\"lamflTh\" format=\"ascii\" NumberOfTuples=\"1\" RangeMin=\""
            << m_solver->m_laminarFlameThickness << "\" RangeMax=\"" << m_solver->m_laminarFlameThickness << "\">"
            << endl;
        ofl << m_solver->m_laminarFlameThickness << endl;
        ofl << "</DataArray>" << endl;
        if(m_solver->m_forcing) {
          // flame strouhal number = 2*pi*f/Lf
          MInt flameStr = m_solver->m_flameStrouhal / (F2 * PI);
          ofl << "<DataArray type=\"" << dataType
              << "\" Name=\"flStr\" format=\"ascii\" NumberOfTuples=\"1\" RangeMin=\"" << flameStr << "\" RangeMax=\""
              << flameStr << "\">" << endl;
          ofl << flameStr << endl;
          ofl << "</DataArray>" << endl;
          // wave number based on Lf or 1 omega = 2*pi/Lf (= m_flameStrouhal)
          ofl << "<DataArray type=\"" << dataType
              << "\" Name=\"waveN\" format=\"ascii\" NumberOfTuples=\"1\" RangeMin=\"" << m_solver->m_flameStrouhal
              << "\" RangeMax=\"" << m_solver->m_flameStrouhal << "\">" << endl;
          ofl << m_solver->m_flameStrouhal << endl;
          ofl << "</DataArray>" << endl;
          ofl << "<DataArray type=\"" << dataType
              << "\" Name=\"forcAmpl\" format=\"ascii\" NumberOfTuples=\"1\" RangeMin=\""
              << m_solver->m_forcingAmplitude << "\" RangeMax=\"" << m_solver->m_forcingAmplitude << "\">" << endl;
          ofl << m_solver->m_forcingAmplitude << endl;
          ofl << "</DataArray>" << endl;
        }
      }

      //================== /FieldData =================

      ofl << "<Piece NumberOfPoints=\"" << noPoints << "\" NumberOfCells=\"" << noCells << "\">" << endl;


      //================== Points =================
      ofl << "<Points>" << endl;
      ofl << "<DataArray type=\"" << dataType << "\" NumberOfComponents=\"3\" format=\"appended\" offset=\"" << offset
          << "\"/>" << endl;
      offset += points.m_memsize + sizeof(MUint);
      ofl << "</Points>" << endl;
      //================== /Points =================


      //================== Cells =================
      ofl << "<Cells>" << endl;

      ofl << "<DataArray type=\"" << uIDataType << "\" Name=\"connectivity\" format=\"appended\" offset=\"" << offset
          << "\"/>" << endl;
      offset += connectivity.m_memsize + sizeof(MUint);

      ofl << "<DataArray type=\"" << uIDataType << "\" Name=\"offsets\" format=\"appended\" offset=\"" << offset
          << "\"/>" << endl;
      offset += offsets.m_memsize + sizeof(MUint);

      ofl << "<DataArray type=\"" << uI8DataType << "\" Name=\"types\" format=\"appended\" offset=\"" << offset
          << "\"/>" << endl;
      offset += types.m_memsize + sizeof(MUint);

      ofl << "</Cells>" << endl;
      //================== /Cells =================


      //================== CellData/PointData =================
      ofl << "<CellData Scalars=\"scalars\">" << endl;

      ofl << "<DataArray type=\"" << iDataType << "\" Name=\"cellIds\" format=\"appended\" offset=\"" << offset
          << "\"/>" << endl;
      offset += cellIds.m_memsize + sizeof(MUint);

      ofl << "<DataArray type=\"" << iDataType << "\" Name=\"bndryIds\" format=\"appended\" offset=\"" << offset
          << "\"/>" << endl;
      offset += bndryIds.m_memsize + sizeof(MUint);
      /*
      ofl << "<DataArray type=\"" << dataType64 << "\" Name=\"drho_dt\" format=\"appended\"
      offset=\""<< offset <<"\"/>" << endl; offset += drho_dt.m_memsize + sizeof( MUint );
      */
      if(m_solver->m_levelSet || m_solver->m_levelSetMb) {
        ofl << "<DataArray type=\"" << dataType << "\" Name=\"levelSetFunction\" format=\"appended\" offset=\""
            << offset << "\"/>" << endl;
        offset += levelSetFunction.m_memsize + sizeof(MUint);
      }
      if(writePointData) {
        ofl << "</CellData>" << endl;
        ofl << "<PointData Scalars=\"scalars\">" << endl;
      }

      ofl << "<DataArray type=\"" << dataType << "\" Name=\"pressure\" format=\"appended\" offset=\"" << offset
          << "\"/>" << endl;
      offset += pressure.m_memsize + sizeof(MUint);

      ofl << "<DataArray type=\"" << dataType
          << "\" Name=\"velocity\" NumberOfComponents=\"3\" format=\"appended\" offset=\"" << offset << "\"/>" << endl;
      offset += velocity.m_memsize + sizeof(MUint);

      ofl << "<DataArray type=\"" << dataType << "\" Name=\"density\" format=\"appended\" offset=\"" << offset << "\"/>"
          << endl;
      offset += density.m_memsize + sizeof(MUint);

      ofl << "<DataArray type=\"" << dataType << "\" Name=\"vorticity\" format=\"appended\" offset=\"" << offset
          << "\"/>" << endl;
      offset += vorticity.m_memsize + sizeof(MUint);

      //    ofl << "<DataArray type=\"" << dataType << "\" Name=\"Qcriterion\" format=\"appended\"
      //    offset=\""<< offset <<"\"/>" << endl;
      // offset += Qcriterion.m_memsize + sizeof( MUint );

      if(m_solver->m_combustion) {
        ofl << "<DataArray type=\"" << dataType << "\" Name=\"progressVariable\" format=\"appended\" offset=\""
            << offset << "\"/>" << endl;
        offset += progressVariable.m_memsize + sizeof(MUint);
      }
      if(writePointData) {
        ofl << "</PointData>" << endl;
      } else {
        ofl << "</CellData>" << endl;
      }
      //================== /CellData/PointData =================


      ofl << "</Piece>" << endl;
      ofl << "</UnstructuredGrid>" << endl;


      //================== AppendedData =================
      ofl << "<AppendedData encoding=\"raw\">" << endl;
      ofl << "_";
      ofl.close();
      ofl.open(fileName, ios_base::out | ios_base::app | ios_base::binary);

      //----------
      // point coordinates
      uinumber = (MUint)points.m_memsize;
      ofl.write(reinterpret_cast<const char*>(&uinumber), sizeof(MUint));
      ofl.write(reinterpret_cast<const char*>(points.getPointer()), uinumber);

      //----------
      // connectivity
      uinumber = (MUint)connectivity.m_memsize;
      ofl.write(reinterpret_cast<const char*>(&uinumber), sizeof(MUint));
      ofl.write(reinterpret_cast<const char*>(connectivity.getPointer()), uinumber);

      //----------
      // offsets
      uinumber = (MUint)offsets.m_memsize;
      ofl.write(reinterpret_cast<const char*>(&uinumber), sizeof(uinumber));
      ofl.write(reinterpret_cast<const char*>(offsets.getPointer()), uinumber);


      //----------
      // cell types
      uinumber = (MUint)types.m_memsize;
      ofl.write(reinterpret_cast<const char*>(&uinumber), sizeof(uinumber));
      ofl.write(reinterpret_cast<const char*>(types.getPointer()), uinumber);

      //----------
      // cellIds
      uinumber = (MUint)cellIds.m_memsize;
      ofl.write(reinterpret_cast<const char*>(&uinumber), sizeof(uinumber));
      ofl.write(reinterpret_cast<const char*>(cellIds.getPointer()), uinumber);

      //----------
      // bndryIds
      uinumber = (MUint)bndryIds.m_memsize;
      ofl.write(reinterpret_cast<const char*>(&uinumber), sizeof(uinumber));
      ofl.write(reinterpret_cast<const char*>(bndryIds.getPointer()), uinumber);
      /*
      //----------
      // drho_dt
      uinumber = (MUint)drho_dt.m_memsize;
      ofl.write(reinterpret_cast<const char *> (&uinumber), sizeof(uinumber));
      ofl.write(reinterpret_cast<const char *> (drho_dt.getPointer()), uinumber);
      */

      //----------
      // levelSetFunction
      if(m_solver->m_levelSet && (m_solver->m_levelSetMb)) {
        uinumber = (MUint)levelSetFunction.m_memsize;
        ofl.write(reinterpret_cast<const char*>(&uinumber), sizeof(uinumber));
        ofl.write(reinterpret_cast<const char*>(levelSetFunction.getPointer()), uinumber);
      }
      //----------
      // pressure
      uinumber = (MUint)pressure.m_memsize;
      ofl.write(reinterpret_cast<const char*>(&uinumber), sizeof(uinumber));
      ofl.write(reinterpret_cast<const char*>(pressure.getPointer()), uinumber);

      //----------
      // velocity
      uinumber = (MUint)velocity.m_memsize;
      ofl.write(reinterpret_cast<const char*>(&uinumber), sizeof(uinumber));
      ofl.write(reinterpret_cast<const char*>(velocity.getPointer()), uinumber);

      //----------
      // density
      uinumber = (MUint)density.m_memsize;
      ofl.write(reinterpret_cast<const char*>(&uinumber), sizeof(uinumber));
      ofl.write(reinterpret_cast<const char*>(density.getPointer()), uinumber);

      //----------
      // vorticity
      uinumber = (MUint)vorticity.m_memsize;
      ofl.write(reinterpret_cast<const char*>(&uinumber), sizeof(uinumber));
      ofl.write(reinterpret_cast<const char*>(vorticity.getPointer()), uinumber);
      /*
      //----------
      // Qcriterion
      uinumber = (MUint)Qcriterion.m_memsize;
      ofl.write(reinterpret_cast<const char *> (&uinumber), sizeof(uinumber));
      ofl.write(reinterpret_cast<const char *> (Qcriterion.getPointer()), uinumber);
      */

      //----------
      // progress variable
      if(m_solver->m_combustion) {
        uinumber = (MUint)progressVariable.m_memsize;
        ofl.write(reinterpret_cast<const char*>(&uinumber), sizeof(uinumber));
        ofl.write(reinterpret_cast<const char*>(progressVariable.getPointer()), uinumber);
      }
      //----------
      ofl.close();
      ofl.open(fileName, ios_base::app);
      ofl << endl;
      ofl << "</AppendedData>" << endl;
      //================== /AppendedData =================


      ofl << "</VTKFile>" << endl;
      ofl.close();
      //================== /VTKFile =================
    } else {
      cerr << "ERROR! COULD NOT OPEN FILE " << fileName << " for writing! " << endl;
    }

    for(MInt i = 0; i < m_solver->m_bndryCells->size(); i++) {
      delete[] pointInside[i];
    }
    delete[] pointInside;
    pointInside = nullptr;
  }
  else IF_CONSTEXPR(nDim == 3) {
    static_cast<void>(fileName); // silence unused-parameter warning for 3D
    MString fileName2 = m_solver->m_solutionOutput + "QOUT_" + to_string(globalTimeStep) + ".vtu";
    if(m_solver->m_multipleFvSolver) {
      fileName2 = m_solver->m_solutionOutput + "QOUT" + to_string(m_solver->m_solverId) + "_"
                  + to_string(globalTimeStep) + ".vtu";
    }
    MString geomFileName = m_solver->m_solutionOutput + "GEOM_" + to_string(globalTimeStep) + ".vtp";
    writeVtuOutputParallel<>(fileName2, geomFileName);
  }
}


// only used in writeVtuOutputParallel
template <MInt nDim, class SysEqn>
MInt VtkIo<nDim, SysEqn>::estimateMemorySizeSolverwise(uint64_t noElements,
                                                       ScratchSpace<uint64_t>& noElementsPerDomain,
                                                       uint64_t factor) {
  ASSERT(noElements == noElementsPerDomain[domainId()], "");
  MInt ret = noElements;
  if(domainId() == 0) ret += (MInt)((sizeof(uint64_t) + factor - 1) / factor);
  return mMax(1, ret);
}

//-------------------------------------------------------------------------------------

/** \brief The zeroth domain stores the header for an uncompressed vtu file appended data array specifying its numer of
 * bytes \author Lennart Scheiders
 */
template <MInt nDim, class SysEqn>
void VtkIo<nDim, SysEqn>::insertDataHeader(char* data, uint64_t& memsize, uint64_t& memsizeGlobal, uint64_t& offset) {
  if(domainId() == 0) {
    if(memsize > 0) {
      memmove(&data[sizeof(uint64_t)], &data[0], memsize);
    }
    uint64_t dataSize = memsizeGlobal;
    memcpy(&data[0], reinterpret_cast<void*>(&dataSize), sizeof(uint64_t));
    memsize += sizeof(uint64_t);
  } else {
    offset += sizeof(uint64_t);
  }
  memsizeGlobal += sizeof(uint64_t);
}

//-------------------------------------------------------------------------------------

/** \brief Recomputes global offsets and data size given the local memsize
 * \author Lennart Scheiders
 */
template <MInt nDim, class SysEqn>
void VtkIo<nDim, SysEqn>::updateDataOffsets(uint64_t memsize, uint64_t& memsizeGlobal, uint64_t& offset,
                                            uint64_t& oldMemsizeGlobal) {
  ScratchSpace<uint64_t> tmpCntPerDomain(noDomains(), AT_, "tmpCntPerDomain");
  MPI_Allgather(&memsize, 1, MPI_UINT64_T, &tmpCntPerDomain[0], 1, MPI_UINT64_T, mpiComm(), AT_, "memsize",
                "tmpCntPerDomain[0]");
  oldMemsizeGlobal += memsizeGlobal;
  memsizeGlobal = 0;
  offset = 0;
  for(MInt d = 0; d < noDomains(); d++) {
    memsizeGlobal += tmpCntPerDomain[d];
  }
  for(MInt d = 0; d < domainId(); d++) {
    offset += tmpCntPerDomain[d];
  }
}

template <MInt nDim, class SysEqn>
template <typename T>
uint64_t VtkIo<nDim, SysEqn>::vtuAssembleFaceStream(vector<T>*& facestream, vector<T>*& conn, uint64_t& facesSize,
                                                    ScratchSpace<MInt>& polyMap, MInt& polyCnt,
                                                    const MBool getInternalPolyhedra) {
  TRACE();

  static constexpr MInt pointCode[6][4] = {{0, 2, 6, 4}, {1, 3, 7, 5}, {0, 1, 5, 4},
                                           {2, 3, 7, 6}, {0, 1, 3, 2}, {4, 5, 7, 6}};
  static constexpr MInt edgeCode[6][4] = {{0, 10, 4, 8},  {1, 11, 5, 9}, {2, 9, 6, 8},
                                          {3, 11, 7, 10}, {2, 1, 3, 0},  {6, 5, 7, 4}};
  static constexpr MInt faceStencil[2][12] = {{0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 0, 1},
                                              {4, 4, 4, 4, 5, 5, 5, 5, 2, 2, 3, 3}};
  // static constexpr MInt revDir[6] = {1,0,3,2,5,4};
  const MFloat eps0 = 0.0001;
  const MFloat eps = 1e-12;
  const MInt noNodes = IPOW2(nDim);
  const MInt noBndryCells = m_solver->m_bndryCells->size();
  MIntScratchSpace cellMap(m_solver->a_noCells(), FUN_, "cellMap");
  // MIntScratchSpace polyMap(m_solver->m_extractedCells->size(), FUN_, "polyMap");
  MIntScratchSpace revMap(m_solver->a_noCells(), FUN_, "revMap");
  MIntScratchSpace edgeCheck(50, 50, FUN_, "edgeCheck");
  MIntScratchSpace nghbrList(200, AT_, "nghbrList");
  MIntScratchSpace layerId(200, AT_, "layerId");
  vector<MInt> noInternalNodes(noBndryCells, 0);
  // noInternalNodes.resize( noBndryCells );
  // facestream.resize( noBndryCells );
  // conn.resize( noBndryCells );
  // facestream = new vector<T> [ mMax(1,noBndryCells) ];
  // conn = new vector<T> [ mMax(1,noBndryCells) ];
  uint64_t connSize = 0;
  facesSize = 0;
  cellMap.fill(-1);
  polyMap.fill(-1);
  revMap.fill(-1);

  for(MInt c = 0; c < m_solver->m_extractedCells->size(); c++) {
    MInt cellId = m_solver->m_extractedCells->a[c].m_cellId;
    cellMap[cellId] = c;
  }
  polyCnt = 0;
  for(MInt c = 0; c < m_solver->m_extractedCells->size(); c++) {
    MInt cellId = m_solver->m_extractedCells->a[c].m_cellId;
    MInt bndryId = m_solver->a_bndryId(cellId);
    if(bndryId > -1) {
      polyMap[c] = polyCnt;
      revMap[polyCnt] = cellId;
      polyCnt++;
    }
  }
  if(getInternalPolyhedra) {
    for(MInt c = 0; c < m_solver->m_extractedCells->size(); c++) {
      MInt cellId = m_solver->m_extractedCells->a[c].m_cellId;
      MInt bndryId = m_solver->a_bndryId(cellId);
      if(bndryId > -1) continue;
      MBool poly = false;
      MBool polyEdge = false;
      for(MInt i = 0; i < m_noDirs; i++) {
        if(m_solver->checkNeighborActive(cellId, i) && m_solver->a_hasNeighbor(cellId, i) > 0) {
          if(m_solver->c_noChildren(m_solver->c_neighborId(cellId, i)) > 0) {
            poly = true;
          }
        }
        for(MInt j = 0; j < m_noDirs; j++) {
          if((i / 2) == (j / 2)) continue;
          MInt nb0 = (m_solver->checkNeighborActive(cellId, i) && m_solver->a_hasNeighbor(cellId, i) > 0)
                         ? m_solver->c_neighborId(cellId, i)
                         : cellId;
          MInt nb1 = (m_solver->checkNeighborActive(cellId, j) && m_solver->a_hasNeighbor(cellId, j) > 0)
                         ? m_solver->c_neighborId(cellId, j)
                         : cellId;
          MInt nb01 = (m_solver->checkNeighborActive(nb0, j) && m_solver->a_hasNeighbor(nb0, j) > 0)
                          ? m_solver->c_neighborId(nb0, j)
                          : cellId;
          MInt nb10 = (m_solver->checkNeighborActive(nb1, i) && m_solver->a_hasNeighbor(nb1, i) > 0)
                          ? m_solver->c_neighborId(nb1, i)
                          : cellId;
          polyEdge = polyEdge || (m_solver->c_noChildren(nb0) > 0) || (m_solver->c_noChildren(nb1) > 0)
                     || (m_solver->c_noChildren(nb01) > 0) || (m_solver->c_noChildren(nb10) > 0);
        }
      }
      if(poly || polyEdge) {
        polyMap[c] = polyCnt;
        revMap[polyCnt] = cellId;
        polyCnt++;
      }
    }
  }

  if(facestream != nullptr) {
    delete[] facestream;
  }
  if(conn != nullptr) {
    delete[] conn;
  }
  facestream = new vector<T>[mMax(1, polyCnt)];
  conn = new vector<T>[mMax(1, polyCnt)];
  // find duplicate points
  /*
  for ( MInt c = 0; c < m_solver->m_extractedCells->size(); c++ ) {
    MInt cellId = m_solver->m_extractedCells->a[c].m_cellId;
    //MInt bndryId = a_bndryId( cellId );
    //if ( bndryId < 0 ) {
    //const MInt counter = getAdjacentLeafCells<2>( cellId, 1, nghbrList, layerId );
    const MInt counter = getAdjacentLeafCells_d2( cellId, 1, nghbrList, layerId );
      for ( MInt n = 0; n < counter; n++ ) {
        MInt nghbrCell = nghbrList[n];
        if ( nghbrCell < 0 ) continue;
        if ( a_isHalo(nghbrCell) ) continue;
        MInt ncell = cellMap[nghbrCell];
        if( ncell < 0 ) continue;
        //if (a_bndryId( nghbrCell )> -1 ) continue;
        for ( MInt p = 0; p < noNodes; p++ ) {
          for ( MInt q = 0; q < noNodes; q++ ) {
            MInt gridPointId = m_solver->m_extractedCells->a[c].m_pointIds[p];
            MInt nghbrGridPointId = m_solver->m_extractedCells->a[ncell].m_pointIds[q];
            if ( gridPointId == nghbrGridPointId ) continue;
            if ( nghbrGridPointId < 0 ) cerr << "negative grid point id " << endl;
            MFloat delta = sqrt( POW2( m_solver->m_gridPoints->a[ nghbrGridPointId ].m_coordinates[0] -
  m_solver->m_gridPoints->a[ gridPointId ].m_coordinates[0] )
                                   + POW2( m_solver->m_gridPoints->a[ nghbrGridPointId ].m_coordinates[1] -
  m_solver->m_gridPoints->a[ gridPointId ].m_coordinates[1] )
                                   + POW2( m_solver->m_gridPoints->a[ nghbrGridPointId ].m_coordinates[2] -
  m_solver->m_gridPoints->a[ gridPointId ].m_coordinates[2] ) ); if ( delta < eps0*c_cellLengthAtCell(cellId) ) {
            //cerr << "found another match 0" << endl;
              MInt pid = ( gridPointId < nghbrGridPointId ) ? gridPointId : nghbrGridPointId;
              MInt pid2 = ( gridPointId < nghbrGridPointId ) ? nghbrGridPointId : gridPointId;
              MInt cid = ( gridPointId < nghbrGridPointId ) ? ncell : c;
              m_solver->m_extractedCells->a[cid].m_pointIds[p] = pid;
              if ( m_solver->m_gridPoints->a[ pid ].m_noAdjacentCells < noNodes ) {
                m_solver->m_gridPoints->a[ pid ].m_cellIds[ m_solver->m_gridPoints->a[ pid ].m_noAdjacentCells ] =
  m_solver->m_extractedCells->a[cid].m_cellId; m_solver->m_gridPoints->a[ pid ].m_noAdjacentCells++;
              }
              for ( MInt k = 0; k < m_solver->m_gridPoints->a[ pid2 ].m_noAdjacentCells; k++ ) {
                if ( m_solver->m_gridPoints->a[ pid2 ].m_cellIds[ k ] == m_solver->m_extractedCells->a[cid].m_cellId ) {
                  m_solver->m_gridPoints->a[ pid2 ].m_cellIds[ k ] = -1;
                  m_solver->m_gridPoints->a[ pid2 ].m_noAdjacentCells--;
                  if ( k < m_solver->m_gridPoints->a[ pid2 ].m_noAdjacentCells )
                    m_solver->m_gridPoints->a[ pid2 ].m_cellIds[ k ] = m_solver->m_gridPoints->a[ pid2 ].m_cellIds[
  m_solver->m_gridPoints->a[ pid2
  ].m_noAdjacentCells ];
                }
              }
              if ( m_solver->m_gridPoints->a[ pid2 ].m_noAdjacentCells == 0 ) {
                for ( MInt j = 0; j < nDim; j++ ) {
                  m_solver->m_gridPoints->a[ pid2 ].m_coordinates[ j ] = std::numeric_limits<float>::max();
                }
              }
            }
          }
        }
      }
      //}
  }
  */


  for(MInt c = 0; c < m_solver->m_extractedCells->size(); c++) {
    MInt cellId = m_solver->m_extractedCells->a[c].m_cellId;
    MInt bndryId = m_solver->a_bndryId(cellId);
    /*MBool poly = false;
    MBool polyEdge = false;
    if ( noEntries > polyCnt0 ) {
      for ( MInt i = 0; i < m_noDirs; i++ ) {
        if ( a_hasNeighbor( cellId ,  i ) > 0 ) {
          if ( c_noChildren( c_neighborId( cellId ,  i ) ) > 0 ) {
            poly = true;
          }
        }
        for ( MInt j = 0; j < m_noDirs; j++ ) {
          if( (i/2) == (j/2) ) continue;
          MInt nb0 = ( a_hasNeighbor( cellId ,  i ) > 0 ) ? c_neighborId( cellId ,  i ) : cellId;
          MInt nb1 = ( a_hasNeighbor( cellId ,  j ) > 0 ) ? c_neighborId( cellId ,  j ) : cellId;
          MInt nb01 = ( a_hasNeighbor( nb0 ,  j ) > 0 ) ? c_neighborId( nb0 ,  j ) : cellId;
          MInt nb10 = ( a_hasNeighbor( nb1 ,  i ) > 0 ) ? c_neighborId( nb1 ,  i ) : cellId;
          polyEdge = polyEdge || ( c_noChildren( nb0 ) > 0 ) || ( c_noChildren( nb1 ) > 0 ) || ( c_noChildren( nb01 ) >
    0 ) || ( c_noChildren( nb10 ) > 0 );
        }
      }
    }
    if ( bndryId < 0 ) {
      if ( !poly && !polyEdge ) {
        connSize += (uint64_t)noNodes;
        facesSize++;
      }
      continue;
    }*/
    MInt polyId = polyMap[c];
    if(polyId < 0) {
      connSize += (uint64_t)noNodes;
      facesSize++;
      continue;
    }
    facestream[polyId].clear();
    conn[polyId].clear();
    for(MInt i = 0; i < noNodes; i++) {
      const MInt gridPointId = m_solver->m_extractedCells->a[c].m_pointIds[i];
      if(gridPointId < 0) continue;
      if(bndryId > -1 && m_solver->gridPointIsInside(c, i)) {
        /*for ( MInt k = 0; k < m_solver->m_gridPoints->a[ gridPointId ].m_noAdjacentCells; k++ ) {
          if ( m_solver->m_gridPoints->a[ gridPointId ].m_cellIds[ k ] == m_solver->m_extractedCells->a[c].m_cellId ) {
            m_solver->m_gridPoints->a[ gridPointId ].m_cellIds[ k ] = -1;
            m_solver->m_gridPoints->a[ gridPointId ].m_noAdjacentCells--;
            if ( k < m_solver->m_gridPoints->a[ gridPointId ].m_noAdjacentCells )
              m_solver->m_gridPoints->a[ gridPointId ].m_cellIds[ k ] = m_solver->m_gridPoints->a[ gridPointId
        ].m_cellIds[ m_solver->m_gridPoints->a[ gridPointId ].m_noAdjacentCells ];
            //break;
          }
        }
        if ( m_solver->m_gridPoints->a[ gridPointId ].m_noAdjacentCells == 0 ) {
          for ( MInt j = 0; j < nDim; j++ ) {
            m_solver->m_gridPoints->a[ gridPointId ].m_coordinates[ j ] = std::numeric_limits<float>::max();
          }
        }*/
        m_solver->m_extractedCells->a[c].m_pointIds[i] = -1;

        /*
        for ( MInt k = 0; k < m_solver->m_gridPoints->a[ gridPointId ].m_noAdjacentCells; k++ ) {
          MInt ncellId = m_solver->m_gridPoints->a[ gridPointId ].m_cellIds[ k ];
          if ( ncellId < 0 ) continue;
          MInt ncell = cellMap[ncellId];
          if ( ncell < 0 ) continue;
          for ( MInt p = 0; p < noNodes; p++ ) {
            const MInt gridPointId2 = m_solver->m_extractedCells->a[ ncell ].m_pointIds[ p ];
            if ( gridPointId2 < 0 ) continue;
            if ( gridPointId == gridPointId2 ) {
              m_solver->m_extractedCells->a[ ncell ].m_pointIds[ p ] = -1;
            }
          }
        }
        for ( MInt p = 0; p < noNodes; p++ ) {
          m_solver->m_gridPoints->a[ gridPointId ].m_cellIds[ p ] = -1;
        }
        m_solver->m_extractedCells->a[ c ].m_pointIds[ i ] = -1;
        m_solver->m_gridPoints->a[ gridPointId ].m_noAdjacentCells = 0;
        for ( MInt j = 0; j < nDim; j++ ) {
          m_solver->m_gridPoints->a[ gridPointId ].m_coordinates[ j ] = std::numeric_limits<float>::max();
        }
        */
      } else {
        conn[polyId].push_back(gridPointId);
        /*MBool exist = false;
        for ( MInt k = 0; k < m_solver->m_gridPoints->a[ gridPointId ].m_noAdjacentCells; k++ ){
          if ( m_solver->m_gridPoints->a[ gridPointId ].m_cellIds[k] == cellId ) exist = true;
        }
        if ( !exist && m_solver->m_gridPoints->a[ gridPointId ].m_noAdjacentCells < noNodes ) {
          m_solver->m_gridPoints->a[ gridPointId ].m_cellIds[ m_solver->m_gridPoints->a[ gridPointId ].m_noAdjacentCells
        ] = cellId; m_solver->m_gridPoints->a[ gridPointId ].m_noAdjacentCells++;
          }*/
      }
    }
  }


  const MInt oldNoPoints = m_solver->m_gridPoints->size();
  MIntScratchSpace pointRefs(oldNoPoints, FUN_, "pointRefs");
  pointRefs.fill(0);
  for(MInt gridPointId = 0; gridPointId < oldNoPoints; gridPointId++) {
    m_solver->m_gridPoints->a[gridPointId].m_noAdjacentCells = 0;
    for(MInt p = 0; p < noNodes; p++) {
      m_solver->m_gridPoints->a[gridPointId].m_cellIds[p] = -1;
    }
  }
  for(MInt c = 0; c < m_solver->m_extractedCells->size(); c++) {
    MInt cellId = m_solver->m_extractedCells->a[c].m_cellId;
    MInt rootId = (m_solver->a_hasProperty(cellId, SolverCell::IsSplitChild))
                      ? m_solver->getAssociatedInternalCell(cellId)
                      : cellId;
    // MInt bndryId = a_bndryId( cellId );
    MInt pc = polyMap[c];
    // MInt noPoints = ( bndryId > -1 && pc > -1 ) ? conn[pc].size() : noNodes;
    MInt noPoints = (pc > -1) ? conn[pc].size() : noNodes;
    for(MInt q = 0; q < noPoints; q++) {
      // for ( MInt q = 0; q < noNodes; q++ ) {
      // MInt gridPointId = ( bndryId > -1 && pc > -1 ) ? conn[pc][q] : m_solver->m_extractedCells->a[c].m_pointIds[q];
      MInt gridPointId = (pc > -1) ? conn[pc][q] : m_solver->m_extractedCells->a[c].m_pointIds[q];
      if(gridPointId >= oldNoPoints) cerr << "point out of range" << endl;
      // MInt gridPointId = m_solver->m_extractedCells->a[c].m_pointIds[q];
      if(gridPointId > -1) {
        MBool found = false;
        for(MInt n = 0; n < m_solver->m_gridPoints->a[gridPointId].m_noAdjacentCells; n++) {
          if(m_solver->m_gridPoints->a[gridPointId].m_cellIds[n] == rootId) found = true;
        }
        if(!found && (m_solver->m_gridPoints->a[gridPointId].m_noAdjacentCells < noNodes)) {
          m_solver->m_gridPoints->a[gridPointId].m_cellIds[m_solver->m_gridPoints->a[gridPointId].m_noAdjacentCells] =
              cellId;
          m_solver->m_gridPoints->a[gridPointId].m_noAdjacentCells++;
        }
        pointRefs[gridPointId]++;
      }
    }
  }
  for(MInt gridPointId = 0; gridPointId < oldNoPoints; gridPointId++) {
    if(pointRefs[gridPointId] > noNodes) cerr << "too many point refs " << pointRefs[gridPointId] << endl;
    if(m_solver->m_gridPoints->a[gridPointId].m_noAdjacentCells != pointRefs[gridPointId]) {
      cerr << "strange0 " << m_solver->m_gridPoints->a[gridPointId].m_noAdjacentCells << " " << pointRefs[gridPointId]
           << endl;
    }
    // if ( pointRefs[gridPointId] == 0 ) {
    //  m_solver->m_gridPoints->a[ gridPointId ].m_noAdjacentCells = 0;
    //}
  }
  /*
    MInt inact = 0;
    for(MInt gridPointId = 0; gridPointId < oldNoPoints; gridPointId++ ) {
      if ( pointRefs[gridPointId] == 0 ) {
        if ( m_solver->m_gridPoints->a[ gridPointId ].m_noAdjacentCells > 0 ) {
          cerr <<"strange " << m_solver->m_gridPoints->a[ gridPointId ].m_noAdjacentCells <<endl;
        inact++;
        }

      }
    }
    cerr << "inact " << inact << endl;*/
  MInt lastId = m_solver->m_gridPoints->size() - 1;
  for(MInt gridPointId = 0; gridPointId < lastId; gridPointId++) {
    if(m_solver->m_gridPoints->a[gridPointId].m_noAdjacentCells > 0) continue;
    // if ( m_solver->m_gridPoints->a[ gridPointId ].m_noAdjacentCells > 0 && pointRefs[gridPointId] > 0 ) continue;
    // if ( pointRefs[gridPointId] > 0 ) continue;
    lastId = m_solver->m_gridPoints->size() - 1;
    while(m_solver->m_gridPoints->a[lastId].m_noAdjacentCells == 0) {
      // while ( pointRefs[gridPointId] == 0 ) {
      m_solver->m_gridPoints->resetSize(lastId);
      lastId = m_solver->m_gridPoints->size() - 1;
    }
    if(gridPointId >= m_solver->m_gridPoints->size()) {
      break;
    } else {
      lastId = m_solver->m_gridPoints->size() - 1;
      // if ( m_solver->m_gridPoints->a[ gridPointId ].m_noAdjacentCells > 0 ) {
      //  cerr << "Warning: grid point still active." << endl;
      //  continue;
      //}
      if(gridPointId < lastId) {
        for(MInt k = 0; k < m_solver->m_gridPoints->a[lastId].m_noAdjacentCells; k++) {
          MInt cellId = m_solver->m_gridPoints->a[lastId].m_cellIds[k];
          if(cellId < 0) continue;
          if(cellMap[cellId] > -1) {
            for(MInt i = 0; i < noNodes; i++) {
              if(m_solver->m_extractedCells->a[cellMap[cellId]].m_pointIds[i] == lastId) {
                m_solver->m_extractedCells->a[cellMap[cellId]].m_pointIds[i] = gridPointId;
              }
            }
            // MInt bndryId = a_bndryId( cellId );
            // if ( bndryId > -1 ) {
            MInt polyId = polyMap[cellMap[cellId]];
            if(polyId > -1) {
              replace(conn[polyId].begin(), conn[polyId].end(), lastId, gridPointId);
            }
            //}
          }
          if(m_solver->a_hasProperty(cellId, SolverCell::IsSplitCell)) {
            auto it = find(m_solver->m_splitCells.begin(), m_solver->m_splitCells.end(), cellId);
            if(it == m_solver->m_splitCells.end()) mTerm(1, AT_, "split cells inconsistency.");
            const MInt pos = distance(m_solver->m_splitCells.begin(), it);
            ASSERT(m_solver->m_splitCells[pos] == cellId, "");
            for(MUint c = 0; c < m_solver->m_splitChilds[pos].size(); c++) {
              MInt splitChilId = m_solver->m_splitChilds[pos][c];
              if(cellMap[splitChilId] > -1) {
                for(MInt i = 0; i < noNodes; i++) {
                  if(m_solver->m_extractedCells->a[cellMap[splitChilId]].m_pointIds[i] == lastId) {
                    m_solver->m_extractedCells->a[cellMap[splitChilId]].m_pointIds[i] = gridPointId;
                  }
                }
                MInt polyId = polyMap[cellMap[splitChilId]];
                if(polyId > -1) {
                  replace(conn[polyId].begin(), conn[polyId].end(), lastId, gridPointId);
                }
              }
            }
          }
        }
        for(MInt k = 0; k < m_solver->m_gridPoints->a[lastId].m_noAdjacentCells; k++) {
          m_solver->m_gridPoints->a[gridPointId].m_cellIds[k] = m_solver->m_gridPoints->a[lastId].m_cellIds[k];
          m_solver->m_gridPoints->a[lastId].m_cellIds[k] = -1;
        }
        m_solver->m_gridPoints->a[gridPointId].m_noAdjacentCells = m_solver->m_gridPoints->a[lastId].m_noAdjacentCells;
        m_solver->m_gridPoints->a[lastId].m_noAdjacentCells = 0;
        pointRefs[gridPointId] = pointRefs[lastId];
        pointRefs[lastId] = 0;
        for(MInt j = 0; j < nDim; j++) {
          m_solver->m_gridPoints->a[gridPointId].m_coordinates[j] = m_solver->m_gridPoints->a[lastId].m_coordinates[j];
          m_solver->m_gridPoints->a[lastId].m_coordinates[j] = std::numeric_limits<float>::max();
        }
      }
      m_solver->m_gridPoints->resetSize(lastId);
    }
  }


  multimap<MFloat, MInt> sortedCPs;
  for(MInt c = 0; c < m_solver->m_extractedCells->size(); c++) {
    for(MInt i = 0; i < noNodes; i++) {
      if(m_solver->m_extractedCells->a[c].m_pointIds[i] >= m_solver->m_gridPoints->size()) {
        cerr << domainId() << ": Error invalid grid point " << c << " " << m_solver->m_extractedCells->a[c].m_cellId
             << " " << m_solver->m_extractedCells->a[c].m_pointIds[i] << " " << m_solver->m_gridPoints->size() << endl;
        m_solver->m_extractedCells->a[c].m_pointIds[i] = 0;
      }
    }
    MInt bndryId = m_solver->a_bndryId(m_solver->m_extractedCells->a[c].m_cellId);
    if(bndryId < 0) continue;
    MInt polyId = polyMap[c];
    for(MUint i = 0; i < conn[polyId].size(); i++) {
      if(conn[polyId][i] >= (uint64_t)m_solver->m_gridPoints->size()) {
        cerr << domainId() << ": Error invalid boundary grid point " << c << " "
             << m_solver->m_extractedCells->a[c].m_cellId << " " << conn[polyId][i] << " "
             << m_solver->m_gridPoints->size() << endl;
        conn[polyId][i] = 0;
      }
    }
    noInternalNodes[bndryId] = 0;
    for(MInt i = 0; i < noNodes; i++) {
      if(m_solver->m_extractedCells->a[c].m_pointIds[i] > -1) {
        noInternalNodes[bndryId]++;
      }
    }
  }
  for(MInt c = 0; c < m_solver->m_extractedCells->size(); c++) {
    MInt bndryId = m_solver->a_bndryId(m_solver->m_extractedCells->a[c].m_cellId);
    if(bndryId < 0) continue;
    MInt polyId = polyMap[c];
    MBool outside[8];
    for(MInt i = 0; i < noNodes; i++) {
      outside[i] = true;
      if(m_solver->m_extractedCells->a[c].m_pointIds[i] > -1) {
        outside[i] = false;
      }
    }

    for(MInt srfc = 0; srfc < m_solver->m_bndryCells->a[bndryId].m_noSrfcs; srfc++) {
      for(MInt cp = 0; cp < m_solver->m_bndryCells->a[bndryId].m_srfcs[srfc]->m_noCutPoints; cp++) {
        MInt cutEdge = m_solver->m_bndryCells->a[bndryId].m_srfcs[srfc]->m_cutEdge[cp];
        MInt gridPointId = -1;
        // check if cut point has been created by a neighbor cell prevoiusly
        if(cutEdge > -1) {
          MInt cellId = m_solver->m_extractedCells->a[c].m_cellId;
          MInt f0 = faceStencil[0][cutEdge];
          MInt f1 = faceStencil[1][cutEdge];
          MInt nghbrs[4] = {-1, -1, -1, -1};
          // MInt nDirs[ 4 ][ 2 ] = {{revDir[f0],f1},{f0,revDir[f1]},{revDir[f0],revDir[f1]},{revDir[f0],revDir[f1]}};
          if(!m_solver->a_hasProperty(cellId, SolverCell::IsSplitChild) && m_solver->checkNeighborActive(cellId, f0)
             && m_solver->a_hasNeighbor(cellId, f0) > 0) {
            nghbrs[0] = m_solver->a_bndryId(m_solver->c_neighborId(cellId, f0));
            if(m_solver->checkNeighborActive(m_solver->c_neighborId(cellId, f0), f1)
               && m_solver->a_hasNeighbor(m_solver->c_neighborId(cellId, f0), f1) > 0) {
              nghbrs[3] = m_solver->a_bndryId(m_solver->c_neighborId(m_solver->c_neighborId(cellId, f0), f1));
            }
          }
          if(!m_solver->a_hasProperty(cellId, SolverCell::IsSplitChild) && m_solver->checkNeighborActive(cellId, f1)
             && m_solver->a_hasNeighbor(cellId, f1) > 0) {
            nghbrs[1] = m_solver->a_bndryId(m_solver->c_neighborId(cellId, f1));
            if(m_solver->checkNeighborActive(m_solver->c_neighborId(cellId, f1), f0)
               && m_solver->a_hasNeighbor(m_solver->c_neighborId(cellId, f1), f0) > 0) {
              nghbrs[2] = m_solver->a_bndryId(m_solver->c_neighborId(m_solver->c_neighborId(cellId, f1), f0));
            }
          }
          for(MInt nBndryId : nghbrs) {
            if(gridPointId > -1) continue;
            if(nBndryId > -1) {
              if(noInternalNodes[nBndryId] == 0) continue;
              if(cellMap[m_solver->m_bndryCells->a[nBndryId].m_cellId] < 0) continue;
              MInt npc = polyMap[cellMap[m_solver->m_bndryCells->a[nBndryId].m_cellId]];
              if(npc < 0) continue;
              MInt noPointsNghbr = conn[npc].size();
              for(MInt q = noInternalNodes[nBndryId]; q < noPointsNghbr; q++) {
                MInt nghbrGridPointId = conn[npc][q];
                if(nghbrGridPointId < 0) cerr << "negative grid point id " << endl;
                MFloat delta =
                    sqrt(POW2(m_solver->m_gridPoints->a[nghbrGridPointId].m_coordinates[0]
                              - m_solver->m_bndryCells->a[bndryId].m_srfcs[srfc]->m_cutCoordinates[cp][0])
                         + POW2(m_solver->m_gridPoints->a[nghbrGridPointId].m_coordinates[1]
                                - m_solver->m_bndryCells->a[bndryId].m_srfcs[srfc]->m_cutCoordinates[cp][1])
                         + POW2(m_solver->m_gridPoints->a[nghbrGridPointId].m_coordinates[2]
                                - m_solver->m_bndryCells->a[bndryId].m_srfcs[srfc]->m_cutCoordinates[cp][2]));
                if(delta < eps * m_solver->c_cellLengthAtCell(cellId)) {
                  // cerr << "found another match 0" << endl;
                  gridPointId = nghbrGridPointId;
                  break;
                }
              }
              /*
                MInt ncnt = 0;
                for( MInt nSrfc = 0; nSrfc < m_solver->m_bndryCells->a[ nBndryId ].m_noSrfcs; nSrfc++ ){
                  for( MInt ncp = 0; ncp < m_solver->m_bndryCells->a[ nBndryId ].m_srfcs[nSrfc]->m_noCutPoints; ncp++ )
                { MInt nCutEdge = m_solver->m_bndryCells->a[ nBndryId ].m_srfcs[nSrfc]->m_cutEdge[ ncp ]; if (
                m_solver->m_bndryCells->a[ bndryId ].m_srfcs[srfc]->m_bodyId[0] == m_solver->m_bndryCells->a[ nBndryId
                ].m_srfcs[nSrfc]->m_bodyId[0] ) { if( nCutEdge > -1 ){ MInt nf0 = faceStencil[0][nCutEdge]; MInt nf1
                = faceStencil[1][nCutEdge]; if (((nf0 == nDirs[n][0]) && (nf1 == nDirs[n][1])) ||
                            ((nf0 == nDirs[n][1]) && (nf1 == nDirs[n][0]))) {
                          MFloat delta = 0;
                          for( MInt i = 0; i < nDim; i++ ){
                            delta += POW2( m_solver->m_bndryCells->a[ nBndryId
                ].m_srfcs[nSrfc]->m_cutCoordinates[ncp][i] - m_solver->m_bndryCells->a[ bndryId
                ].m_srfcs[srfc]->m_cutCoordinates[cp][i] );
                          }
                          if( sqrt(delta) < eps*c_cellLengthAtCell(cellId) ) {
                            MInt nbCell = m_solver->m_bndryCells->a[ nBndryId ].m_cellId;
                            if ( cellMap[nbCell] > -1 && polyMap[cellMap[nbCell]] > -1 ) {
                              if ( noInternalNodes[nBndryId]+ncnt < (signed)conn[ polyMap[cellMap[nbCell]] ].size() ) {
                                gridPointId = conn[ polyMap[cellMap[nbCell]] ][noInternalNodes[nBndryId]+ncnt];
                              }
                            }
                          }
                        }
                      }
                    }
                    ncnt++;
                  }
                }*/
            }
          }
        }
        if(gridPointId < 0) {
          MInt cellId = m_solver->m_extractedCells->a[c].m_cellId;
          // const MInt counter = getAdjacentLeafCells<2>( cellId, 1, nghbrList, layerId );
          const MInt counter = m_solver->getAdjacentLeafCells_d2(cellId, 1, nghbrList, layerId);
          for(MInt n = 0; n < counter; n++) {
            MInt nghbrCell = nghbrList[n];
            if(nghbrCell < 0) continue;
            if(m_solver->a_isHalo(nghbrCell)) continue;
            MInt ncell = cellMap[nghbrCell];
            if(ncell < 0) continue;
            MInt pc = polyMap[ncell];
            if(m_solver->a_bndryId(nghbrCell) < 0) continue;
            if(pc < 0) continue;
            MInt nBndryId = m_solver->a_bndryId(nghbrCell);
            MInt noPointsNghbr = conn[pc].size();
            for(MInt q = noInternalNodes[nBndryId]; q < noPointsNghbr; q++) {
              MInt nghbrGridPointId = conn[pc][q];
              if(nghbrGridPointId < 0) cerr << "negative grid point id " << endl;
              MFloat delta = sqrt(POW2(m_solver->m_gridPoints->a[nghbrGridPointId].m_coordinates[0]
                                       - m_solver->m_bndryCells->a[bndryId].m_srfcs[srfc]->m_cutCoordinates[cp][0])
                                  + POW2(m_solver->m_gridPoints->a[nghbrGridPointId].m_coordinates[1]
                                         - m_solver->m_bndryCells->a[bndryId].m_srfcs[srfc]->m_cutCoordinates[cp][1])
                                  + POW2(m_solver->m_gridPoints->a[nghbrGridPointId].m_coordinates[2]
                                         - m_solver->m_bndryCells->a[bndryId].m_srfcs[srfc]->m_cutCoordinates[cp][2]));
              if(delta < eps * m_solver->c_cellLengthAtCell(cellId)) {
                // cerr << "found another match 0" << endl;
                gridPointId = nghbrGridPointId;
                break;
              }
            }
          }
        }
        // if not, create it
        if(gridPointId < 0) {
          gridPointId = m_solver->m_gridPoints->size();
          m_solver->m_gridPoints->append();
          m_solver->m_gridPoints->a[gridPointId].m_noAdjacentCells = 0;
          for(MInt i = 0; i < nDim; i++) {
            m_solver->m_gridPoints->a[gridPointId].m_coordinates[i] =
                m_solver->m_bndryCells->a[bndryId].m_srfcs[srfc]->m_cutCoordinates[cp][i];
          }
        }
        if(m_solver->m_gridPoints->a[gridPointId].m_noAdjacentCells >= noNodes) {
          cerr << "nodes " << gridPointId << " " << m_solver->m_extractedCells->a[c].m_cellId << " "
               << m_solver->m_gridPoints->a[gridPointId].m_noAdjacentCells << "; ";
          for(MInt k = 0; k < m_solver->m_gridPoints->a[gridPointId].m_noAdjacentCells; k++)
            cerr << m_solver->m_gridPoints->a[gridPointId].m_cellIds[k] << " ";
          cerr << endl;
        }
        ASSERT(m_solver->m_gridPoints->a[gridPointId].m_noAdjacentCells < noNodes, "");
        if(m_solver->m_gridPoints->a[gridPointId].m_noAdjacentCells < noNodes) {
          MBool exist = false;
          MInt cellId = m_solver->m_extractedCells->a[c].m_cellId;
          MInt rootId = (m_solver->a_hasProperty(cellId, SolverCell::IsSplitChild))
                            ? m_solver->getAssociatedInternalCell(cellId)
                            : cellId;
          for(MInt k = 0; k < m_solver->m_gridPoints->a[gridPointId].m_noAdjacentCells; k++) {
            if(m_solver->m_gridPoints->a[gridPointId].m_cellIds[k] == rootId) exist = true;
          }
          if(!exist) {
            m_solver->m_gridPoints->a[gridPointId].m_cellIds[m_solver->m_gridPoints->a[gridPointId].m_noAdjacentCells] =
                rootId;
            m_solver->m_gridPoints->a[gridPointId].m_noAdjacentCells++;
          }
        }
        conn[polyId].push_back(gridPointId);
      }
    }

    if(m_solver->m_bndryCells->a[bndryId].m_faceVertices.empty()) {
      facestream[polyId].push_back(0); // number of faces
      MInt pointCntOffset = 1;

      MInt cnt = 0;
      for(MInt srfc = 0; srfc < m_solver->m_bndryCells->a[bndryId].m_noSrfcs; srfc++) {
        sortedCPs.clear();
        MFloat pCoords[3]{};
        MFloat vec_a[3]{};
        MFloat vec_b[3]{};
        MFloat normal[3]{};
        for(MInt i = 0; i < nDim; i++) {
          normal[i] = m_solver->m_bndryCells->a[bndryId].m_srfcs[srfc]->m_normalVector[i];
        }
        MInt spaceId = 0;
        MFloat maxC = fabs(normal[0]);
        for(MInt i = 1; i < nDim; i++) {
          if(fabs(normal[i]) > maxC) {
            maxC = fabs(normal[i]);
            spaceId = i;
          }
        }
        MInt spaceId1 = (spaceId + 1) % nDim;
        MInt spaceId2 = (spaceId1 + 1) % nDim;
        vec_a[spaceId1] = F1;
        vec_a[spaceId2] = F1;
        vec_a[spaceId] = -(vec_a[spaceId1] * normal[spaceId1] + vec_a[spaceId2] * normal[spaceId2]) / normal[spaceId];
        MFloat vecsum = sqrt(POW2(vec_a[0]) + POW2(vec_a[1]) + POW2(vec_a[2]));
        for(MInt i = 0; i < nDim; i++) {
          vec_a[i] *= 1000.0 / vecsum;
        }
        vec_b[spaceId] = normal[spaceId1] * vec_a[spaceId2] - normal[spaceId2] * vec_a[spaceId1];
        vec_b[spaceId1] = normal[spaceId2] * vec_a[spaceId] - normal[spaceId] * vec_a[spaceId2];
        vec_b[spaceId2] = normal[spaceId] * vec_a[spaceId1] - normal[spaceId1] * vec_a[spaceId];
        vecsum = sqrt(POW2(vec_b[0]) + POW2(vec_b[1]) + POW2(vec_b[2]));
        for(MInt i = 0; i < nDim; i++) {
          vec_b[i] *= 1000.0 * vecsum;
        }
        for(MInt cp = 0; cp < m_solver->m_bndryCells->a[bndryId].m_srfcs[srfc]->m_noCutPoints; cp++) {
          for(MInt i = 0; i < nDim; i++) {
            pCoords[i] = m_solver->m_bndryCells->a[bndryId].m_srfcs[srfc]->m_cutCoordinates[cp][i];
          }
          MFloat dx = F0;
          MFloat dy = F0;
          for(MInt i = 0; i < nDim; i++) {
            dx += vec_a[i] * (pCoords[i] - m_solver->m_bndryCells->a[bndryId].m_srfcs[srfc]->m_coordinates[i]);
            dy += vec_b[i] * (pCoords[i] - m_solver->m_bndryCells->a[bndryId].m_srfcs[srfc]->m_coordinates[i]);
          }
          sortedCPs.insert(make_pair(atan2(dy, dx), cp));
        }
        ASSERT((signed)sortedCPs.size() == m_solver->m_bndryCells->a[bndryId].m_srfcs[srfc]->m_noCutPoints, "");
        facestream[polyId].push_back(0); // number of points for current face
        for(auto& sortedCP : sortedCPs) {
          facestream[polyId].push_back(conn[polyId][noInternalNodes[bndryId] + cnt + sortedCP.second]);
          facestream[polyId][pointCntOffset]++;
        }
        facestream[polyId][0]++;
        pointCntOffset = facestream[polyId].size();
        cnt += m_solver->m_bndryCells->a[bndryId].m_srfcs[srfc]->m_noCutPoints;
      }


      for(MInt face = 0; face < m_noDirs; face++) {
        // if ( m_solver->m_extractedCells->a[ c ].m_pointIds[ pointCode[face][0] ] < 0 &&
        // m_solver->m_extractedCells->a[ c ].m_pointIds[ pointCode[face][1] ] < 0
        //     && m_solver->m_extractedCells->a[ c ].m_pointIds[ pointCode[face][2] ] < 0 &&
        //     m_solver->m_extractedCells->a[ c ].m_pointIds[ pointCode[face][3] ] < 0 ) continue;
        sortedCPs.clear();
        MFloat pCoords[MAX_SPACE_DIMENSIONS] = {F0, F0, F0};
        MFloat vec_a[MAX_SPACE_DIMENSIONS] = {F0, F0, F0};
        MFloat vec_b[MAX_SPACE_DIMENSIONS] = {F0, F0, F0};
        MInt spaceId = face / 2;
        MInt spaceId1 = (spaceId + 1) % nDim;
        MInt spaceId2 = (spaceId1 + 1) % nDim;
        vec_a[spaceId1] = F1;
        vec_b[spaceId2] = F1;
        MFloat sum = F0;
        for(MInt p = 0; p < 4; p++) {
          if(!outside[pointCode[face][p]]) {
            for(MInt i = 0; i < nDim; i++) {
              pCoords[i] += m_solver->m_gridPoints->a[m_solver->m_extractedCells->a[c].m_pointIds[pointCode[face][p]]]
                                .m_coordinates[i];
            }
            sum += F1;
          }
          MInt ccnt = 0;
          for(MInt srfc = 0; srfc < m_solver->m_bndryCells->a[bndryId].m_noSrfcs; srfc++) {
            for(MInt cp = 0; cp < m_solver->m_bndryCells->a[bndryId].m_srfcs[srfc]->m_noCutPoints; cp++) {
              if(m_solver->m_bndryCells->a[bndryId].m_srfcs[srfc]->m_cutEdge[cp] == edgeCode[face][p]) {
                for(MInt i = 0; i < nDim; i++) {
                  pCoords[i] +=
                      m_solver->m_gridPoints->a[conn[polyId][noInternalNodes[bndryId] + ccnt]].m_coordinates[i];
                }
                sum += F1;
              }
              ccnt++;
            }
          }
        }
        for(MInt i = 0; i < nDim; i++) {
          pCoords[i] /= sum;
        }

        for(MInt p = 0; p < 4; p++) {
          if(!outside[pointCode[face][p]]) {
            MFloat dx = F0;
            MFloat dy = F0;
            for(MInt i = 0; i < nDim; i++) {
              dx += vec_a[i]
                    * (m_solver->m_gridPoints->a[m_solver->m_extractedCells->a[c].m_pointIds[pointCode[face][p]]]
                           .m_coordinates[i]
                       - pCoords[i]);
              dy += vec_b[i]
                    * (m_solver->m_gridPoints->a[m_solver->m_extractedCells->a[c].m_pointIds[pointCode[face][p]]]
                           .m_coordinates[i]
                       - pCoords[i]);
            }
            sortedCPs.insert(make_pair(atan2(dy, dx), m_solver->m_extractedCells->a[c].m_pointIds[pointCode[face][p]]));
          }
          MInt ccnt = 0;
          for(MInt srfc = 0; srfc < m_solver->m_bndryCells->a[bndryId].m_noSrfcs; srfc++) {
            for(MInt cp = 0; cp < m_solver->m_bndryCells->a[bndryId].m_srfcs[srfc]->m_noCutPoints; cp++) {
              if(m_solver->m_bndryCells->a[bndryId].m_srfcs[srfc]->m_cutEdge[cp] == edgeCode[face][p]) {
                MFloat dx = F0;
                MFloat dy = F0;
                for(MInt i = 0; i < nDim; i++) {
                  dx += vec_a[i]
                        * (m_solver->m_gridPoints->a[conn[polyId][noInternalNodes[bndryId] + ccnt]].m_coordinates[i]
                           - pCoords[i]);
                  dy += vec_b[i]
                        * (m_solver->m_gridPoints->a[conn[polyId][noInternalNodes[bndryId] + ccnt]].m_coordinates[i]
                           - pCoords[i]);
                }
                sortedCPs.insert(
                    make_pair(atan2(dy, dx), static_cast<MInt>(conn[polyId][noInternalNodes[bndryId] + ccnt])));
              }
              ccnt++;
            }
          }
        }

        if(!sortedCPs.empty()) {
          facestream[polyId].push_back(sortedCPs.size()); // number of points for current face
          for(auto& sortedCP : sortedCPs) {
            facestream[polyId].push_back(sortedCP.second);
          }
        } else
          continue;

        /*
        facestream[polyId].push_back(0); //number of points for current face
        for ( MInt p = 0; p < 4; p++ ) {
          if ( !outside[ pointCode[face][p] ] ) {
            facestream[polyId].push_back( m_solver->m_extractedCells->a[ c ].m_pointIds[ pointCode[face][p] ] );
            facestream[polyId][pointCntOffset]++;
          }
          MInt ccnt = 0;
          for( MInt srfc = 0; srfc < m_solver->m_bndryCells->a[bndryId].m_noSrfcs; srfc++ ) {
            for( MInt cp = 0; cp < m_solver->m_bndryCells->a[ bndryId ].m_srfcs[srfc]->m_noCutPoints; cp++ ) {
              if ( m_solver->m_bndryCells->a[ bndryId ].m_srfcs[srfc]->m_cutEdge[cp] == edgeCode[face][p] ) {
                facestream[polyId].push_back( conn[polyId][noInternalNodes[bndryId]+ccnt] );
                facestream[polyId][pointCntOffset]++;
              }
              ccnt++;
            }
          }
        }*/


        facestream[polyId][0]++;
        pointCntOffset = facestream[polyId].size();
      }

      sort(conn[polyId].begin(), conn[polyId].end(), less<T>());
      conn[polyId].erase(unique(conn[polyId].begin(), conn[polyId].end()), conn[polyId].end());
      connSize += (uint64_t)conn[polyId].size();
      facesSize += (uint64_t)facestream[polyId].size();
    }
  }

  // cells with split faces and other special multi cuts
  for(MInt c = 0; c < m_solver->m_extractedCells->size(); c++) {
    MInt cellId = m_solver->m_extractedCells->a[c].m_cellId;
    MInt bndryId = m_solver->a_bndryId(cellId);
    if(bndryId < 0) continue;
    MInt polyId = polyMap[c];

    if(!m_solver->m_bndryCells->a[bndryId].m_faceVertices.empty()) {
      // const MInt gridPointOffset = m_solver->m_gridPoints->size();
      facestream[polyId].push_back(0); // number of faces
      // econn[polyId].clear();
      vector<MInt> tmpPointMap;
      tmpPointMap.resize(m_solver->m_bndryCells->a[bndryId].m_faceVertices.size() / 3);
      for(MUint v = 0; v < m_solver->m_bndryCells->a[bndryId].m_faceVertices.size() / 3; v++) {
        tmpPointMap[v] = -1;
      }
      for(MUint v = 0; v < m_solver->m_bndryCells->a[bndryId].m_faceVertices.size() / 3; v++) {
        MFloat ccoords[3]{};
        for(MInt i = 0; i < nDim; i++) {
          ccoords[i] = m_solver->m_bndryCells->a[bndryId].m_faceVertices[3 * v + i];
        }
        MInt gridPointId = -1;
        MInt noPoints = conn[polyId].size();
        for(MInt q = 0; q < noPoints; q++) {
          MInt tmpGridPointId = conn[polyId][q];
          if(tmpGridPointId < 0) cerr << "negative grid point id " << endl;
          MFloat delta = sqrt(POW2(m_solver->m_gridPoints->a[tmpGridPointId].m_coordinates[0] - ccoords[0])
                              + POW2(m_solver->m_gridPoints->a[tmpGridPointId].m_coordinates[1] - ccoords[1])
                              + POW2(m_solver->m_gridPoints->a[tmpGridPointId].m_coordinates[2] - ccoords[2]));
          if(delta < eps * m_solver->c_cellLengthAtCell(cellId)) {
            gridPointId = tmpGridPointId;
          }
        }
        if(gridPointId < 0) {
          // const MInt counter = m_solver->getAdjacentLeafCells<2>( cellId, 1, nghbrList, layerId );
          const MInt counter = m_solver->getAdjacentLeafCells_d2(cellId, 1, nghbrList, layerId);
          for(MInt n = 0; n < counter; n++) {
            MInt nghbrCell = nghbrList[n];
            if(nghbrCell < 0) continue;
            if(m_solver->a_isHalo(nghbrCell)) continue;
            MInt ncell = cellMap[nghbrCell];
            if(ncell < 0) continue;
            MInt pc = polyMap[ncell];
            // if (m_solver->a_bndryId( nghbrCell )> -1 ) continue;
            MInt noPointsNghbr = (pc > -1) ? conn[pc].size() : noNodes;
            for(MInt q = 0; q < noPointsNghbr; q++) {
              MInt nghbrGridPointId = (pc > -1) ? conn[pc][q] : m_solver->m_extractedCells->a[ncell].m_pointIds[q];
              if(nghbrGridPointId < 0) cerr << "negative grid point id " << endl;
              MFloat delta = sqrt(POW2(m_solver->m_gridPoints->a[nghbrGridPointId].m_coordinates[0] - ccoords[0])
                                  + POW2(m_solver->m_gridPoints->a[nghbrGridPointId].m_coordinates[1] - ccoords[1])
                                  + POW2(m_solver->m_gridPoints->a[nghbrGridPointId].m_coordinates[2] - ccoords[2]));
              if(delta < eps * m_solver->c_cellLengthAtCell(cellId)) {
                // cerr << "found another match 1" << endl;
                gridPointId = nghbrGridPointId;
                conn[polyId].push_back(gridPointId);
                break;
              }
            }
          }
        }
        if(gridPointId < 0) {
          gridPointId = m_solver->m_gridPoints->size();
          m_solver->m_gridPoints->append();
          m_solver->m_gridPoints->a[gridPointId].m_noAdjacentCells = 0;
          for(MInt i = 0; i < nDim; i++) {
            m_solver->m_gridPoints->a[gridPointId].m_coordinates[i] =
                m_solver->m_bndryCells->a[bndryId].m_faceVertices[3 * v + i];
          }
          MInt rootId = (m_solver->a_hasProperty(cellId, SolverCell::IsSplitChild))
                            ? m_solver->getAssociatedInternalCell(cellId)
                            : cellId;
          m_solver->m_gridPoints->a[gridPointId].m_cellIds[m_solver->m_gridPoints->a[gridPointId].m_noAdjacentCells] =
              rootId;
          m_solver->m_gridPoints->a[gridPointId].m_noAdjacentCells++;
          conn[polyId].push_back(gridPointId);
        }
        tmpPointMap[v] = gridPointId;
      }
      for(MUint f = 0; f < m_solver->m_bndryCells->a[bndryId].m_faceStream.size(); f++) {
        facestream[polyId][0]++;
        facestream[polyId].push_back(m_solver->m_bndryCells->a[bndryId].m_faceStream[f].size());
        for(MUint v = 0; v < m_solver->m_bndryCells->a[bndryId].m_faceStream[f].size(); v++) {
          if(tmpPointMap[m_solver->m_bndryCells->a[bndryId].m_faceStream[f][v]] < 0)
            cerr << "invalid split point" << endl;
          // facestream[polyId].push_back(gridPointOffset+m_solver->m_bndryCells->a[bndryId].m_faceStream[f][v]);
          facestream[polyId].push_back(tmpPointMap[m_solver->m_bndryCells->a[bndryId].m_faceStream[f][v]]);
        }
      }
      sort(conn[polyId].begin(), conn[polyId].end(), less<T>());
      conn[polyId].erase(unique(conn[polyId].begin(), conn[polyId].end()), conn[polyId].end());
      connSize += (uint64_t)conn[polyId].size();
      facesSize += (uint64_t)facestream[polyId].size();
    }
  }


  if(getInternalPolyhedra) {
    const MInt ltable[6][4] = {{0, 2, 6, 4}, {1, 3, 7, 5}, {0, 1, 5, 4}, {2, 3, 7, 6}, {0, 1, 3, 2}, {4, 5, 7, 6}};
    const MInt ltable2[6][2] = {{1, 7}, {0, 6}, {2, 7}, {0, 5}, {4, 7}, {0, 3}};
    const MInt ltable3[6][4] = {{4, 3, 5, 2}, {4, 3, 5, 2}, {4, 1, 5, 0}, {4, 1, 5, 0}, {2, 1, 3, 0}, {2, 1, 3, 0}};
    const MFloat signStencil[8][3] = {{-F1, -F1, -F1}, {F1, -F1, -F1}, {-F1, F1, -F1}, {F1, F1, -F1},
                                      {-F1, -F1, F1},  {F1, -F1, F1},  {-F1, F1, F1},  {F1, F1, F1}};
    const MInt dirPlusOne[4] = {1, 2, 3, 0};
    const MInt dirMinusOne[4] = {3, 0, 1, 2};
    const MInt dirPlusTwo[4] = {2, 3, 0, 1};
    const MInt reverseDir[6] = {1, 0, 3, 2, 5, 4};


    for(MInt c = 0; c < m_solver->m_extractedCells->size(); c++) {
      const MInt cellId = m_solver->m_extractedCells->a[c].m_cellId;
      const MInt bndryId = m_solver->a_bndryId(cellId);
      if(bndryId > -1) continue;
      MInt polyId = polyMap[c];
      if(polyId < 0) continue;
      if(polyId >= polyCnt) mTerm(1, AT_, "entries out of range");
      facestream[polyId].clear();
      /*conn[polyId].clear();
      for ( MInt i = 0; i < noNodes; i++ ) {
        const MInt gridPointId = m_solver->m_extractedCells->a[ c ].m_pointIds[ i ];
        if ( gridPointId < 0 ) continue;
        conn[polyId].push_back( gridPointId );
      }*/
      facestream[polyId].push_back(0); // number of faces
      MInt pointCntOffset = 1;
      MBool polyFaceSide[6];
      for(MInt i = 0; i < m_noDirs; i++) {
        MBool polySide = false;
        if(m_solver->checkNeighborActive(cellId, i) && m_solver->a_hasNeighbor(cellId, i) > 0) {
          if(m_solver->c_noChildren(m_solver->c_neighborId(cellId, i)) > 0) {
            polySide = true;
          }
        }
        polyFaceSide[i] = polySide;
        if(polySide) {
          // store four faces (polygon numbering!)
          MInt edgeIds[4] = {-1, -1, -1, -1};
          MInt centerId = -1;
          MInt nghbrId0 = m_solver->c_childId(m_solver->c_neighborId(cellId, i), ltable2[i][0]);
          if(nghbrId0 > -1 && m_solver->a_bndryId(nghbrId0) > -1) nghbrId0 = -1;
          MInt nghbrId = (nghbrId0 < 0) ? -1 : cellMap[nghbrId0];

          MFloat ccoords[3]{};
          for(MInt j = 0; j < nDim; j++) {
            ccoords[j] = m_solver->a_coordinate(cellId, j);
          }
          ccoords[i / 2] +=
              (i % 2 == 0) ? -F1B2 * m_solver->c_cellLengthAtCell(cellId) : F1B2 * m_solver->c_cellLengthAtCell(cellId);

          MInt gridPointId = -1;
          for(MInt e = 0; e < (MInt)conn[polyId].size(); e++) {
            MFloat delta = F0;
            for(MInt j = 0; j < nDim; j++) {
              delta += POW2(ccoords[j] - m_solver->m_gridPoints->a[conn[polyId][e]].m_coordinates[j]);
            }
            delta = sqrt(delta);
            if(delta < eps0 * m_solver->c_cellLengthAtCell(cellId)) {
              gridPointId = conn[polyId][e];
              break;
            }
          }
          if(gridPointId < 0 && nghbrId > -1) {
            gridPointId = m_solver->m_extractedCells->a[nghbrId].m_pointIds[ltable2[i][1]];
            conn[polyId].push_back(gridPointId);
          }
          if(gridPointId < 0) {
            // const MInt counter = m_solver->getAdjacentLeafCells<2>( cellId, 1, nghbrList, layerId );
            const MInt counter = m_solver->getAdjacentLeafCells_d2(cellId, 1, nghbrList, layerId);
            for(MInt n = 0; n < counter; n++) {
              MInt nghbrCell = nghbrList[n];
              if(nghbrCell < 0) continue;
              if(m_solver->a_isHalo(nghbrCell)) continue;
              MInt ncell = cellMap[nghbrCell];
              if(ncell < 0) continue;
              MInt pc = polyMap[ncell];
              // if (m_solver->a_bndryId( nghbrCell )> -1 ) continue;
              MInt noPointsNghbr = (pc > -1) ? conn[pc].size() : noNodes;
              for(MInt q = 0; q < noPointsNghbr; q++) {
                MInt nghbrGridPointId = (pc > -1) ? conn[pc][q] : m_solver->m_extractedCells->a[ncell].m_pointIds[q];
                if(nghbrGridPointId < 0) cerr << "negative grid point id " << endl;
                MFloat delta = sqrt(POW2(m_solver->m_gridPoints->a[nghbrGridPointId].m_coordinates[0] - ccoords[0])
                                    + POW2(m_solver->m_gridPoints->a[nghbrGridPointId].m_coordinates[1] - ccoords[1])
                                    + POW2(m_solver->m_gridPoints->a[nghbrGridPointId].m_coordinates[2] - ccoords[2]));
                if(delta < eps * m_solver->c_cellLengthAtCell(cellId)) {
                  // cerr << "found another match 1" << endl;
                  gridPointId = nghbrGridPointId;
                  conn[polyId].push_back(gridPointId);
                  break;
                }
              }
            }
          }
          if(gridPointId < 0) {
            gridPointId = m_solver->m_gridPoints->size();
            m_solver->m_gridPoints->append();
            m_solver->m_gridPoints->a[gridPointId].m_noAdjacentCells = 0;
            for(MInt j = 0; j < nDim; j++) {
              m_solver->m_gridPoints->a[gridPointId].m_coordinates[j] = ccoords[j];
            }
            conn[polyId].push_back(gridPointId);
          }
          if(m_solver->m_gridPoints->a[gridPointId].m_noAdjacentCells < noNodes) {
            MBool exist = false;
            MInt rootId = (m_solver->a_hasProperty(cellId, SolverCell::IsSplitChild))
                              ? m_solver->getAssociatedInternalCell(cellId)
                              : cellId;
            for(MInt k = 0; k < m_solver->m_gridPoints->a[gridPointId].m_noAdjacentCells; k++) {
              if(m_solver->m_gridPoints->a[gridPointId].m_cellIds[k] == rootId) exist = true;
            }
            if(!exist) {
              m_solver->m_gridPoints->a[gridPointId]
                  .m_cellIds[m_solver->m_gridPoints->a[gridPointId].m_noAdjacentCells] = rootId;
              m_solver->m_gridPoints->a[gridPointId].m_noAdjacentCells++;
            }
          }
          centerId = gridPointId;

          for(MInt k = 0; k < 4; k++) {
            MInt cpt = ltable[reverseDir[i]][k];
            MInt npt = ltable[reverseDir[i]][dirPlusOne[k]];
            gridPointId = -1;
            for(MInt j = 0; j < nDim; j++) {
              ccoords[j] = m_solver->a_coordinate(cellId, j);
            }
            ccoords[i / 2] +=
                (i % 2 == 0) ? -m_solver->c_cellLengthAtCell(cellId) : m_solver->c_cellLengthAtCell(cellId);
            for(MInt j = 0; j < nDim; j++) {
              ccoords[j] += signStencil[cpt][j] * F1B4 * m_solver->c_cellLengthAtCell(cellId);
              ccoords[j] += signStencil[npt][j] * F1B4 * m_solver->c_cellLengthAtCell(cellId);
            }
            for(MUint e = 0; e < conn[polyId].size(); e++) {
              MFloat delta = F0;
              for(MInt j = 0; j < nDim; j++) {
                delta += POW2(ccoords[j] - m_solver->m_gridPoints->a[conn[polyId][e]].m_coordinates[j]);
              }
              delta = sqrt(delta);
              if(delta < eps0 * m_solver->c_cellLengthAtCell(cellId)) {
                gridPointId = conn[polyId][e];
                e = (MInt)conn[polyId].size();
              }
            }
            if(gridPointId < 0 && nghbrId > -1) {
              MInt nghbr0 = m_solver->c_childId(m_solver->c_neighborId(cellId, i), ltable[reverseDir[i]][k]);
              if(nghbr0 > -1 && m_solver->a_bndryId(nghbr0) > -1) nghbr0 = -1;
              MInt nghbrId1 = (nghbr0 < 0) ? -1 : cellMap[nghbr0];
              if(nghbr0 > -1 && nghbrId1 > -1) {
                gridPointId = m_solver->m_extractedCells->a[nghbrId1].m_pointIds[ltable[reverseDir[i]][dirPlusOne[k]]];
                conn[polyId].push_back(gridPointId);
              }
            }
            if(gridPointId < 0) {
              // const MInt counter = m_solver->getAdjacentLeafCells<2>( cellId, 1, nghbrList, layerId );
              const MInt counter = m_solver->getAdjacentLeafCells_d2(cellId, 1, nghbrList, layerId);
              for(MInt n = 0; n < counter; n++) {
                MInt nghbrCell = nghbrList[n];
                if(nghbrCell < 0) continue;
                if(m_solver->a_isHalo(nghbrCell)) continue;
                MInt ncell = cellMap[nghbrCell];
                if(ncell < 0) continue;
                MInt pc = polyMap[ncell];
                // if (a_bndryId( nghbrCell )> -1 ) continue;
                MInt noPointsNghbr = (pc > -1) ? conn[pc].size() : noNodes;
                for(MInt q = 0; q < noPointsNghbr; q++) {
                  MInt nghbrGridPointId = (pc > -1) ? conn[pc][q] : m_solver->m_extractedCells->a[ncell].m_pointIds[q];
                  if(nghbrGridPointId < 0) cerr << "negative grid point id " << endl;
                  MFloat delta =
                      sqrt(POW2(m_solver->m_gridPoints->a[nghbrGridPointId].m_coordinates[0] - ccoords[0])
                           + POW2(m_solver->m_gridPoints->a[nghbrGridPointId].m_coordinates[1] - ccoords[1])
                           + POW2(m_solver->m_gridPoints->a[nghbrGridPointId].m_coordinates[2] - ccoords[2]));
                  if(delta < eps * m_solver->c_cellLengthAtCell(cellId)) {
                    // cerr << "found another match 2" << endl;
                    gridPointId = nghbrGridPointId;
                    conn[polyId].push_back(gridPointId);
                    break;
                  }
                }
              }
            }
            if(gridPointId < 0) {
              gridPointId = m_solver->m_gridPoints->size();
              m_solver->m_gridPoints->append();
              m_solver->m_gridPoints->a[gridPointId].m_noAdjacentCells = 0;
              for(MInt j = 0; j < nDim; j++) {
                m_solver->m_gridPoints->a[gridPointId].m_coordinates[j] = ccoords[j];
              }
              conn[polyId].push_back(gridPointId);
            }
            if(m_solver->m_gridPoints->a[gridPointId].m_noAdjacentCells < noNodes) {
              MBool exist = false;
              MInt rootId = (m_solver->a_hasProperty(cellId, SolverCell::IsSplitChild))
                                ? m_solver->getAssociatedInternalCell(cellId)
                                : cellId;
              for(MInt q = 0; q < m_solver->m_gridPoints->a[gridPointId].m_noAdjacentCells; q++) {
                if(m_solver->m_gridPoints->a[gridPointId].m_cellIds[q] == rootId) exist = true;
              }
              if(!exist) {
                m_solver->m_gridPoints->a[gridPointId]
                    .m_cellIds[m_solver->m_gridPoints->a[gridPointId].m_noAdjacentCells] = rootId;
                m_solver->m_gridPoints->a[gridPointId].m_noAdjacentCells++;
              }
            }
            edgeIds[k] = gridPointId;
          }

          for(MInt k = 0; k < 4; k++) {
            if((signed)centerId >= m_solver->m_gridPoints->size()
               || edgeIds[dirMinusOne[k]] >= m_solver->m_gridPoints->size()
               || m_solver->m_extractedCells->a[c].m_pointIds[ltable[i][k]] >= m_solver->m_gridPoints->size()
               || edgeIds[k] >= m_solver->m_gridPoints->size()) {
              cerr << "node out of range" << endl;
              polyFaceSide[i] = false;
            }
          }
          if(!polyFaceSide[i]) continue;
          for(MInt k = 0; k < 4; k++) {
            facestream[polyId].push_back(0); // number of points for current face
            facestream[polyId].push_back(centerId);
            facestream[polyId].push_back(edgeIds[dirMinusOne[k]]);
            facestream[polyId].push_back(m_solver->m_extractedCells->a[c].m_pointIds[ltable[i][k]]);
            facestream[polyId].push_back(edgeIds[k]);
            facestream[polyId][pointCntOffset] += 4;

            if(m_solver->m_extractedCells->a[c].m_pointIds[ltable[i][k]] < 0)
              cerr << "Regular vertex not found" << endl;
            if(centerId < 0 || edgeIds[dirMinusOne[k]] < 0
               || m_solver->m_extractedCells->a[c].m_pointIds[ltable[i][k]] < 0 || edgeIds[k] < 0)
              cerr << "invalid index" << endl;
            const MFloat dxb2 = F1B2 * m_solver->c_cellLengthAtCell(cellId);
            const MFloat dxeps = 0.000001 * m_solver->c_cellLengthAtCell(cellId);
            MFloat pcoord[3];
            for(MInt j = 0; j < nDim; j++) {
              pcoord[j] = m_solver->a_coordinate(cellId, j);
            }
            pcoord[i / 2] += (i % 2 == 0) ? -dxb2 : dxb2;

            for(MInt e = pointCntOffset + 1; e < (MInt)facestream[polyId].size(); e++) {
              gridPointId = facestream[polyId][e];
              if(gridPointId < 0) cerr << "Invalid point in connectivity set." << endl;
              if(count(conn[polyId].begin(), conn[polyId].end(), gridPointId) == 0)
                cerr << "Point not in connectivity set." << endl;
              if(fabs(m_solver->m_gridPoints->a[gridPointId].m_coordinates[i / 2] - pcoord[i / 2]) > dxeps) {
                cerr << "vertex out of plane (n) " << cellId << " " << i << " " << e << endl;
              }
              for(MInt j = 0; j < nDim; j++) {
                if(j == i / 2) continue;
                if((fabs(m_solver->m_gridPoints->a[gridPointId].m_coordinates[j] - pcoord[j]) > dxeps)
                   && ((fabs(m_solver->m_gridPoints->a[gridPointId].m_coordinates[j] - pcoord[j]) - dxb2) > dxeps)) {
                  cerr << "vertex out of plane (t) " << cellId << " " << i << " " << e << endl;
                }
              }
            }
            if(i % 2 == 1)
              reverse(facestream[polyId].begin() + pointCntOffset + 1,
                      facestream[polyId].begin() + facestream[polyId].size());

            pointCntOffset = facestream[polyId].size();
            facestream[polyId][0]++; // add face
          }
        }
      }
      for(MInt i = 0; i < m_noDirs; i++) {
        if(!polyFaceSide[i]) {
          facestream[polyId].push_back(0); // number of points for current face
          for(MInt k = 0; k < 4; k++) {
            facestream[polyId].push_back(m_solver->m_extractedCells->a[c].m_pointIds[ltable[i][k]]);
            facestream[polyId][pointCntOffset]++;
            MInt nb0 = (m_solver->checkNeighborActive(cellId, i) && m_solver->a_hasNeighbor(cellId, i) > 0)
                           ? m_solver->c_neighborId(cellId, i)
                           : cellId;
            MInt nb1 = (m_solver->checkNeighborActive(cellId, ltable3[i][k])
                        && m_solver->a_hasNeighbor(cellId, ltable3[i][k]) > 0)
                           ? m_solver->c_neighborId(cellId, ltable3[i][k])
                           : cellId;
            MInt nb01 =
                (m_solver->checkNeighborActive(nb0, ltable3[i][k]) && m_solver->a_hasNeighbor(nb0, ltable3[i][k]) > 0)
                    ? m_solver->c_neighborId(nb0, ltable3[i][k])
                    : cellId;
            MInt nb10 = (m_solver->checkNeighborActive(nb1, i) && m_solver->a_hasNeighbor(nb1, i) > 0)
                            ? m_solver->c_neighborId(nb1, i)
                            : cellId;
            MBool isPolyEdge = (m_solver->c_noChildren(nb0) > 0) || (m_solver->c_noChildren(nb1) > 0)
                               || (m_solver->c_noChildren(nb01) > 0) || (m_solver->c_noChildren(nb10) > 0);
            if(polyFaceSide[ltable3[i][k]] || isPolyEdge) {
              MInt nghbr00 = m_solver->c_neighborId(cellId, ltable3[i][k]);
              MInt nghbr0 = (nghbr00 < 0 || m_solver->c_noChildren(nghbr00))
                                ? -1
                                : m_solver->c_childId(nghbr00, ltable[i][dirPlusTwo[k]]);
              if(nghbr0 > -1 && m_solver->a_bndryId(nghbr0) > -1) nghbr0 = -1;
              MInt nghbrId = (nghbr0 < 0) ? -1 : cellMap[nghbr0];
              if(!polyFaceSide[ltable3[i][k]]) nghbrId = -1;
              MInt cpt = ltable[i][dirPlusTwo[k]];
              MInt npt = ltable[i][dirMinusOne[k]];
              MFloat pcoord[3] = {F0, F0, F0};
              for(MInt j = 0; j < nDim; j++) {
                pcoord[j] = m_solver->a_coordinate(cellId, j);
              }
              pcoord[ltable3[i][k] / 2] += (ltable3[i][k] % 2 == 0) ? -m_solver->c_cellLengthAtCell(cellId)
                                                                    : m_solver->c_cellLengthAtCell(cellId);
              for(MInt j = 0; j < nDim; j++) {
                pcoord[j] += signStencil[cpt][j] * F1B4 * m_solver->c_cellLengthAtCell(cellId);
                pcoord[j] += signStencil[npt][j] * F1B4 * m_solver->c_cellLengthAtCell(cellId);
              }
              MInt gridPointId = -1;
              for(MInt e = 0; e < (MInt)conn[polyId].size(); e++) {
                MFloat delta = F0;
                for(MInt j = 0; j < nDim; j++) {
                  delta += POW2(pcoord[j] - m_solver->m_gridPoints->a[conn[polyId][e]].m_coordinates[j]);
                }
                delta = sqrt(delta);
                if(delta < eps0 * m_solver->c_cellLengthAtCell(cellId)) {
                  gridPointId = conn[polyId][e];
                  break;
                }
              }
              if(gridPointId < 0 && nghbrId > -1) {
                gridPointId = m_solver->m_extractedCells->a[nghbrId].m_pointIds[ltable[i][dirMinusOne[k]]];
                MFloat delta = F0;
                for(MInt j = 0; j < nDim; j++) {
                  delta += POW2(pcoord[j] - m_solver->m_gridPoints->a[gridPointId].m_coordinates[j]);
                }
                delta = sqrt(delta);
                if(delta > eps0 * m_solver->c_cellLengthAtCell(cellId)) {
                  cerr << "unexpected coordinate mismatch" << endl;
                  gridPointId = -1;
                } else {
                  conn[polyId].push_back(gridPointId);
                }
              }
              if(gridPointId < 0) {
                // const MInt counter = m_solver->getAdjacentLeafCells<2>( cellId, 1, nghbrList, layerId );
                const MInt counter = m_solver->getAdjacentLeafCells_d2(cellId, 1, nghbrList, layerId);
                for(MInt n = 0; n < counter; n++) {
                  MInt nghbrCell = nghbrList[n];
                  if(nghbrCell < 0) continue;
                  if(m_solver->a_isHalo(nghbrCell)) continue;
                  MInt ncell = cellMap[nghbrCell];
                  if(ncell < 0) continue;
                  MInt pc = polyMap[ncell];
                  // if (a_bndryId( nghbrCell )> -1 ) continue;
                  MInt noPointsNghbr = (pc > -1) ? conn[pc].size() : noNodes;
                  for(MInt q = 0; q < noPointsNghbr; q++) {
                    MInt nghbrGridPointId =
                        (pc > -1) ? conn[pc][q] : m_solver->m_extractedCells->a[ncell].m_pointIds[q];
                    if(nghbrGridPointId < 0) cerr << "negative grid point id " << endl;
                    MFloat delta =
                        sqrt(POW2(m_solver->m_gridPoints->a[nghbrGridPointId].m_coordinates[0] - pcoord[0])
                             + POW2(m_solver->m_gridPoints->a[nghbrGridPointId].m_coordinates[1] - pcoord[1])
                             + POW2(m_solver->m_gridPoints->a[nghbrGridPointId].m_coordinates[2] - pcoord[2]));
                    if(delta < eps * m_solver->c_cellLengthAtCell(cellId)) {
                      // cerr << "found another match 3" << endl;
                      gridPointId = nghbrGridPointId;
                      conn[polyId].push_back(gridPointId);
                      break;
                    }
                  }
                }
              }
              if(gridPointId < 0) {
                gridPointId = m_solver->m_gridPoints->size();
                m_solver->m_gridPoints->append();
                m_solver->m_gridPoints->a[gridPointId].m_noAdjacentCells = 0;
                for(MInt j = 0; j < nDim; j++) {
                  m_solver->m_gridPoints->a[gridPointId].m_coordinates[j] = pcoord[j];
                }
                conn[polyId].push_back(gridPointId);
              }
              if(m_solver->m_gridPoints->a[gridPointId].m_noAdjacentCells < noNodes) {
                MBool exist = false;
                MInt rootId = (m_solver->a_hasProperty(cellId, SolverCell::IsSplitChild))
                                  ? m_solver->getAssociatedInternalCell(cellId)
                                  : cellId;
                for(MInt q = 0; q < m_solver->m_gridPoints->a[gridPointId].m_noAdjacentCells; q++) {
                  if(m_solver->m_gridPoints->a[gridPointId].m_cellIds[q] == rootId) exist = true;
                }
                if(!exist) {
                  m_solver->m_gridPoints->a[gridPointId]
                      .m_cellIds[m_solver->m_gridPoints->a[gridPointId].m_noAdjacentCells] = rootId;
                  m_solver->m_gridPoints->a[gridPointId].m_noAdjacentCells++;
                }
              }
              facestream[polyId].push_back(gridPointId);
              facestream[polyId][pointCntOffset]++;
            }
          }

          const MFloat dxb2 = F1B2 * m_solver->c_cellLengthAtCell(cellId);
          const MFloat dxeps = 0.000001 * m_solver->c_cellLengthAtCell(cellId);
          MFloat pcoord[3];
          for(MInt j = 0; j < nDim; j++) {
            pcoord[j] = m_solver->a_coordinate(cellId, j);
          }
          pcoord[i / 2] += (i % 2 == 0) ? -dxb2 : dxb2;
          for(MInt e = pointCntOffset + 1; e < (MInt)facestream[polyId].size(); e++) {
            MInt gridPointId = facestream[polyId][e];
            if(count(conn[polyId].begin(), conn[polyId].end(), gridPointId) == 0)
              cerr << "Point not in connectivity set." << endl;
            if(fabs(m_solver->m_gridPoints->a[gridPointId].m_coordinates[i / 2] - pcoord[i / 2]) > dxeps) {
              cerr << "vertex out of plane 2 (n) " << cellId << " " << i << " " << e << endl;
            }
            for(MInt j = 0; j < nDim; j++) {
              if(j == i / 2) continue;
              if((fabs(m_solver->m_gridPoints->a[gridPointId].m_coordinates[j] - pcoord[j]) > dxeps)
                 && ((fabs(m_solver->m_gridPoints->a[gridPointId].m_coordinates[j] - pcoord[j]) - dxb2) > dxeps)) {
                cerr << "vertex out of plane 2 (t) " << cellId << " " << i << " " << e << endl;
              }
            }
          }
          if(i % 2 == 1)
            reverse(facestream[polyId].begin() + pointCntOffset + 1,
                    facestream[polyId].begin() + facestream[polyId].size());

          pointCntOffset = facestream[polyId].size();
          facestream[polyId][0]++; // add face
        }
      }
      sort(conn[polyId].begin(), conn[polyId].end(), less<T>());
      conn[polyId].erase(unique(conn[polyId].begin(), conn[polyId].end()), conn[polyId].end());

      for(MInt e = 0; e < (MInt)conn[polyId].size(); e++) {
        for(MInt f = e + 1; f < (MInt)conn[polyId].size(); f++) {
          MInt gridPointId0 = conn[polyId][e];
          MInt gridPointId1 = conn[polyId][f];
          MFloat delta = F0;
          for(MInt j = 0; j < nDim; j++) {
            delta += POW2(m_solver->m_gridPoints->a[gridPointId1].m_coordinates[j]
                          - m_solver->m_gridPoints->a[gridPointId0].m_coordinates[j]);
          }
          delta = sqrt(delta);
          if(delta < 0.0001 * m_solver->c_cellLengthAtCell(cellId)) {
            cerr << setprecision(16) << "duplicate point " << e << " " << f << " " << gridPointId0 << " "
                 << gridPointId1 << " " << delta << endl;
          }
        }
      }

      connSize += (uint64_t)conn[polyId].size();
      facesSize += (uint64_t)facestream[polyId].size();
    }
  }


  { // check polyhedra
    MIntScratchSpace pointRef(50, FUN_, "pointRef");
    MInt errorCnt = 0;
    const MInt maxNoVertices = 4 * IPOW2(nDim - 1);
    // originally: 8
    // for multipleLevelSet: 4*IPOW2(nDim-1) +4
    for(MInt polyId = 0; polyId < polyCnt; polyId++) {
      const MInt cellId = revMap[polyId];
      const MInt bndryId = m_solver->a_bndryId(cellId);
      MBool errorFlag = false;
      if(conn[polyId].empty() || facestream[polyId].empty())
        cerr << domainId() << " (" << cellId << "/" << bndryId << "/" << m_solver->c_globalId(cellId)
             << "): empty face stream" << endl;
      const MFloat cellLength = m_solver->c_cellLengthAtCell(cellId);
      const MFloat deltaH = eps0 * cellLength;
      // const MFloat deltaH2 = 0.05*cellLength;
      const MFloat deltaH3 = eps * cellLength;
      const auto noVerts = (MInt)conn[polyId].size();
      for(MInt e = 0; e < noVerts; e++) {
        if(conn[polyId][e] >= (unsigned)m_solver->m_gridPoints->size()) {
          cerr << domainId() << " (" << cellId << "/" << bndryId << "/" << m_solver->c_globalId(cellId)
               << "): vertex out of range: " << conn[polyId][e] << " " << conn[polyId][e] << " "
               << m_solver->m_gridPoints->size() << endl;
          errorFlag = true;
        }
        for(MInt f = e + 1; f < noVerts; f++) {
          MInt gridPointId0 = conn[polyId][e];
          MInt gridPointId1 = conn[polyId][f];
          MFloat delta = F0;
          for(MInt j = 0; j < nDim; j++) {
            delta += POW2(m_solver->m_gridPoints->a[gridPointId1].m_coordinates[j]
                          - m_solver->m_gridPoints->a[gridPointId0].m_coordinates[j]);
          }
          delta = sqrt(delta);
          if((bndryId > -1 && delta < deltaH3) || (bndryId < 0 && delta < deltaH)) {
            cerr << setprecision(16) << domainId() << " (" << cellId << "/" << bndryId << "/"
                 << m_solver->c_globalId(cellId) << "): duplicate point " << e << " " << f << " " << gridPointId0 << " "
                 << gridPointId1 << " " << delta << " " << delta / cellLength << " /split "
                 << m_solver->a_hasProperty(cellId, SolverCell::IsSplitCell) << " "
                 << m_solver->a_hasProperty(cellId, SolverCell::IsSplitChild) << " "
                 << m_solver->a_hasProperty(cellId, SolverCell::HasSplitFace) << endl;
            errorFlag = true;
          }
        }
      }

      MInt cnt = 0;
      const MInt noFaces = facestream[polyId][cnt++];
      if(noFaces < 4 || noFaces > 24) {
        cerr << domainId() << " (" << cellId << "/" << bndryId << "/" << m_solver->c_globalId(cellId)
             << "): no faces strange: " << noFaces << endl;
        errorFlag = true;
      }
      MInt noEdges = 0;
      pointRef.fill(0);
      edgeCheck.fill(0);
      map<MUint, MInt> pointMap;
      MInt tmpCnt = 0;
      for(MInt e = 0; e < noVerts; e++) {
        if(pointMap.count(conn[polyId][e]) > 0) {
          cerr << domainId() << " (" << cellId << "/" << bndryId << "/" << m_solver->c_globalId(cellId)
               << "): Duplicate vertex: " << e << " " << pointMap.count(conn[polyId][e]) << endl;
          errorFlag = true;
        }
        pointMap[conn[polyId][e]] = tmpCnt++;
      }
      if(tmpCnt != noVerts) {
        cerr << domainId() << " (" << cellId << "/" << bndryId << "/" << m_solver->c_globalId(cellId)
             << "): temp count mismatch: " << tmpCnt << " " << noVerts << endl;
        errorFlag = true;
      }

      for(MInt f = 0; f < noFaces; f++) {
        const MInt noPoints = facestream[polyId][cnt++];
        if(noPoints < 3 || noPoints > maxNoVertices) {
          cerr << domainId() << " (" << cellId << "/" << bndryId << "/" << m_solver->c_globalId(cellId)
               << "): no points strange: " << noPoints << endl;
          errorFlag = true;
        }
        for(MInt p = 0; p < noPoints; p++) {
          const MUint pointId = facestream[polyId][cnt + p];
          if(pointId >= (unsigned)m_solver->m_gridPoints->size()) {
            cerr << domainId() << " (" << cellId << "/" << bndryId << "/" << m_solver->c_globalId(cellId)
                 << "): vertex out of range: " << pointId << " " << m_solver->m_gridPoints->size() << endl;
            errorFlag = true;
          }
          for(MInt q = p + 1; q < noPoints; q++) {
            if(facestream[polyId][cnt + q] == pointId) {
              cerr << domainId() << " (" << cellId << "/" << bndryId << "/" << m_solver->c_globalId(cellId)
                   << "): duplicate vertex face: " << p << " " << q << " " << pointId << " "
                   << facestream[polyId][cnt + q] << endl;
              errorFlag = true;
            }
          }
        }
        MInt p0 = facestream[polyId][cnt];
        MInt p1 = facestream[polyId][cnt + 1];
        MInt p2 = facestream[polyId][cnt + 2];
        MFloat test = F0;
        MInt p2cnt = 3;
        MFloat d1 = F0;
        MFloat d2 = F0;
        for(MInt j = 0; j < nDim; j++) {
          test += (m_solver->m_gridPoints->a[p1].m_coordinates[j] - m_solver->m_gridPoints->a[p0].m_coordinates[j])
                  * (m_solver->m_gridPoints->a[p2].m_coordinates[j] - m_solver->m_gridPoints->a[p0].m_coordinates[j]);
          d1 += POW2(m_solver->m_gridPoints->a[p1].m_coordinates[j] - m_solver->m_gridPoints->a[p0].m_coordinates[j]);
          d2 += POW2(m_solver->m_gridPoints->a[p2].m_coordinates[j] - m_solver->m_gridPoints->a[p0].m_coordinates[j]);
        }
        d1 = sqrt(d1);
        d2 = sqrt(d2);
        while(fabs(fabs(test / mMax(1e-12, d1 * d2)) - F1) < 0.1 && p2cnt < noPoints) {
          p2 = facestream[polyId][cnt + p2cnt];
          d1 = F0;
          d2 = F0;
          test = F0;
          for(MInt j = 0; j < nDim; j++) {
            test += (m_solver->m_gridPoints->a[p1].m_coordinates[j] - m_solver->m_gridPoints->a[p0].m_coordinates[j])
                    * (m_solver->m_gridPoints->a[p2].m_coordinates[j] - m_solver->m_gridPoints->a[p0].m_coordinates[j]);
            d1 += POW2(m_solver->m_gridPoints->a[p1].m_coordinates[j] - m_solver->m_gridPoints->a[p0].m_coordinates[j]);
            d2 += POW2(m_solver->m_gridPoints->a[p2].m_coordinates[j] - m_solver->m_gridPoints->a[p0].m_coordinates[j]);
          }
          d1 = sqrt(d1);
          d2 = sqrt(d2);
          p2cnt++;
        }
        MFloat a[3]{}, b[3]{}, c[3]{}, d[3]{};
        if(noPoints >= 3) {
          for(MInt j = 0; j < nDim; j++) {
            a[j] = m_solver->m_gridPoints->a[p1].m_coordinates[j] - m_solver->m_gridPoints->a[p0].m_coordinates[j];
            b[j] = m_solver->m_gridPoints->a[p2].m_coordinates[j] - m_solver->m_gridPoints->a[p0].m_coordinates[j];
          }
          c[0] = a[1] * b[2] - a[2] * b[1];
          c[1] = a[2] * b[0] - a[0] * b[2];
          c[2] = a[0] * b[1] - a[1] * b[0];
          MFloat cabs = F0;
          for(MInt j = 0; j < nDim; j++) {
            cabs += POW2(c[j]);
          }
          cabs = sqrt(cabs);
          for(MInt j = 0; j < nDim; j++) {
            c[j] /= mMax(1e-12, cabs);
          }

          for(MInt p = 0; p < noPoints; p++) {
            const MInt pointId = facestream[polyId][cnt + p];
            for(MInt j = 0; j < nDim; j++) {
              d[j] =
                  m_solver->m_gridPoints->a[pointId].m_coordinates[j] - m_solver->m_gridPoints->a[p0].m_coordinates[j];
            }
            test = F0;
            for(MInt j = 0; j < nDim; j++) {
              test += c[j] * d[j];
            }
            // if ( (bndryId > -1 && (fabs(test) > deltaH2) ) || (bndryId < 0 && (fabs(test) > deltaH) ) ) {
            if((bndryId < 0 && (fabs(test) > deltaH))) {
              cerr << setprecision(16) << domainId() << " (" << cellId << "/" << bndryId << "/"
                   << m_solver->c_globalId(cellId) << "): vertex out of plane: " << p << " " << test << " "
                   << test / cellLength << endl;
              errorFlag = true;
            }
          }
        }

        MInt firstPoint = facestream[polyId][cnt];
        for(MInt p = 0; p < noPoints; p++) {
          const MInt pointId = facestream[polyId][cnt++];
          if(count(conn[polyId].begin(), conn[polyId].end(), pointId) != 1) {
            cerr << domainId() << " (" << cellId << "/" << bndryId << "/" << m_solver->c_globalId(cellId)
                 << "): vertex not found or duplicate: " << p << " " << pointId << " "
                 << count(conn[polyId].begin(), conn[polyId].end(), pointId) << endl;
            errorFlag = true;
          }
          MInt nextPoint = (p == noPoints - 1) ? firstPoint : facestream[polyId][cnt];
          MInt v0 = pointMap[pointId];
          MInt v1 = pointMap[nextPoint];
          pointRef[v0]++;
          if(v0 == v1) {
            cerr << domainId() << " (" << cellId << "/" << bndryId << "/" << m_solver->c_globalId(cellId)
                 << "): edge invalid: " << v0 << " " << v1 << " " << endl;
            errorFlag = true;
          }
          if(v1 < v0) swap(v0, v1);
          edgeCheck(v0, v1)++;
        }
      }

      MInt noVerts0 = 0;
      for(MInt e = 0; e < noVerts; e++) {
        if(pointRef[e] > 0) noVerts0++;
        /*
        if ( pointRef[e] == 0 ) {
          cerr << domainId() << " (" << cellId << "/" << bndryId << "/" << c_globalId(cellId)
               << "): unused vertex: " << e << " " << pointRef[e] << " " << noVerts
               << " /split " << m_solver->a_hasProperty( cellId ,  SolverCell::IsSplitCell ) << " " <<
        m_solver->a_hasProperty( cellId , SolverCell::IsSplitChild )
               << " " << m_solver->a_hasProperty( cellId ,  SolverCell::HasSplitFace ) << endl;
        }
        else */
        if(pointRef[e] != 0 && (pointRef[e] < 2 || pointRef[e] > 4)) {
          cerr << domainId() << " (" << cellId << "/" << bndryId << "/" << m_solver->c_globalId(cellId)
               << "): cell points not watertight: " << e << " " << pointRef[e] << " " << noVerts << endl;
          errorFlag = true;
        }
        for(MInt f = 0; f < noVerts; f++) {
          if(edgeCheck(e, f) != 0 && edgeCheck(e, f) != 2) {
            cerr << domainId() << " (" << cellId << "/" << bndryId << "/" << m_solver->c_globalId(cellId)
                 << "): cell not watertight: " << e << " " << f << " " << edgeCheck(e, f) << endl;
            errorFlag = true;
          }
          if(edgeCheck(e, f) == 2) noEdges++;
        }
      }

      if(noVerts0 - noEdges + noFaces != 2) {
        cerr << domainId() << " (" << cellId << "/" << bndryId << "/" << m_solver->c_globalId(cellId)
             << "): Euler's formula not fulfilled: " << noVerts << " " << noEdges << " " << noFaces << " /split "
             << m_solver->a_hasProperty(cellId, SolverCell::IsSplitCell) << " "
             << m_solver->a_hasProperty(cellId, SolverCell::IsSplitChild) << " "
             << m_solver->a_hasProperty(cellId, SolverCell::HasSplitFace) << endl;
        errorFlag = true;
      }

      if(bndryId > -1) {
        // if ( !m_solver->m_bndryCells->a[bndryId].m_faceVertices.empty() ) errorFlag = true;
      }

      if(errorFlag && errorCnt < 10) {
        ofstream ofl;
        string fileName = m_solver->m_solutionOutput + "polyhedron_" + to_string(m_solver->c_globalId(cellId)) + ".vtu";
        ofl.open(fileName.c_str(), ios_base::out | ios_base::trunc);
        if(ofl.is_open() && ofl.good()) {
          ofl << "<?xml version=\"1.0\"?>" << endl;
          ofl << R"(<VTKFile type="UnstructuredGrid" version="1.0" byte_order="LittleEndian">)" << endl;
          ofl << "<UnstructuredGrid>" << endl;
          ofl << "<Piece NumberOfPoints=\"" << noVerts << "\" NumberOfCells=\"" << 1 << "\">" << endl;
          ofl << "<Points>" << endl;
          ofl << R"(<DataArray type="Float32" NumberOfComponents="3" format="ascii">)" << endl;
          map<MUint, MInt> pointMap2;
          tmpCnt = 0;
          for(MInt v = 0; v < noVerts; v++) {
            MInt gridPointId = conn[polyId][v];
            pointMap2[gridPointId] = tmpCnt++;
            for(MInt j = 0; j < nDim; j++) {
              ofl << m_solver->m_gridPoints->a[gridPointId].m_coordinates[j] << " ";
            }
            ofl << endl;
          }
          ofl << "</DataArray>" << endl;
          ofl << "</Points>" << endl;
          ofl << "<Cells>" << endl;
          ofl << R"(<DataArray type="Int32" Name="connectivity" format="ascii">)" << endl;
          for(MInt v = 0; v < noVerts; v++) {
            ofl << v << " ";
          }
          ofl << endl;
          ofl << "</DataArray>" << endl;
          ofl << R"(<DataArray type="Int32" Name="offsets" format="ascii">)" << endl;
          ofl << noVerts << endl;
          ofl << "</DataArray>" << endl;
          ofl << R"(<DataArray type="Int32" Name="types" format="ascii">)" << endl;
          ofl << "42" << endl;
          ofl << "</DataArray>" << endl;
          ofl << R"(<DataArray type="Int32" Name="faces" format="ascii">)" << endl;
          MInt fcnt = 0;
          MInt noFace = facestream[polyId][fcnt++];
          ofl << noFace << " ";
          for(MInt f = 0; f < noFace; f++) {
            MInt noPoints = facestream[polyId][fcnt++];
            ofl << noPoints << " ";
            for(MInt p = 0; p < noPoints; p++) {
              ofl << pointMap2[facestream[polyId][fcnt++]] << " ";
            }
          }
          ofl << endl;
          ofl << "</DataArray>" << endl;
          ofl << R"(<DataArray type="Int32" Name="faceoffsets" format="ascii">)" << endl;
          ofl << facestream[polyId].size() << endl;
          ofl << "</DataArray>" << endl;
          ofl << "</Cells>" << endl;
          ofl << "</Piece>" << endl;
          ofl << "</UnstructuredGrid>" << endl;
          ofl << "</VTKFile>" << endl;
          ofl.close();
          ofl.clear();
        }
      }
      if(errorFlag) errorCnt++;
    }
  }

  /*
  const MInt oldNoPoints0 = m_solver->m_gridPoints->size();
  MIntScratchSpace pointRefs2(oldNoPoints0, FUN_, "pointRefs2");
  pointRefs2.fill(0);
  for ( MInt c = 0; c < m_solver->m_extractedCells->size(); c++ ) {
    //MInt cellId = m_solver->m_extractedCells->a[c].m_cellId;
    //MInt bndryId = a_bndryId( cellId );
    MInt pc = polyMap[c];
    MInt noPoints = ( pc > -1 ) ? conn[pc].size() : noNodes;
    for ( MInt q = 0; q < noPoints; q++ ) {
      MInt gridPointId = ( pc > -1 ) ? conn[pc][q] : m_solver->m_extractedCells->a[c].m_pointIds[q];
      if ( gridPointId < 0 ) continue;
      pointRefs2[gridPointId]++;
    }
  }
  MInt inact = 0;
  for(MInt gridPointId = 0; gridPointId < oldNoPoints0; gridPointId++ ) {
    if ( pointRefs2[gridPointId] == 0 ) {
      if ( m_solver->m_gridPoints->a[ gridPointId ].m_noAdjacentCells > 0 ) {
        cerr <<"strange " << m_solver->m_gridPoints->a[ gridPointId ].m_noAdjacentCells <<endl;
      inact++;
      }

    }
  }
  MPI_Allreduce(MPI_IN_PLACE, &inact, 1, MPI_INT, MPI_SUM, mpiComm(), AT_, "MPI_IN_PLACE", "inact" );
  cerr << "inact " << inact << endl;
  */

  return connSize;
}

template <MInt nDim, class SysEqn>
template <typename uint_t>
void VtkIo<nDim, SysEqn>::writeVtuOutputParallel(const MString fileName, const MString geomFileName,
                                                 const MInt noSolverSpecificVars,
                                                 const MFloatScratchSpace& solverSpecificVars, const MInt noUserVars,
                                                 const MFloatScratchSpace& userVars) {
  TRACE();

  // NOTE: currently VTU output is not working on claix

  static_assert(std::is_same<uint_t, uint32_t>::value || std::is_same<uint_t, uint64_t>::value,
                "Invalid template type specified.");
  static_assert(sizeof(float) == 4, "VTU output will fail since float has an unexpected size on this architecture.");
  const char* const uIDataType = (std::is_same<uint_t, uint64_t>::value) ? "UInt64" : "UInt32";
  const char* const uI8DataType = "UInt8";
  const char* const uI64DataType = "UInt64";
  const char* const dataType = "Float32";
  const uint64_t noNodes = IPOW2(nDim);
  const uint8_t baseCellType = (nDim == 3) ? 11 : 8;
  const uint8_t polyCellType = 42;

  NEW_TIMER_GROUP_STATIC(t_vtuTimer, "writeVtuOutputParallel");
  NEW_TIMER_STATIC(t_vtutotal, "VtuOutputParallel", t_vtuTimer);
  NEW_SUB_TIMER_STATIC(t_faces, "faces", t_vtutotal);
  NEW_SUB_TIMER_STATIC(t_init, "init", t_vtutotal);
  NEW_SUB_TIMER_STATIC(t_prepare, "prepare", t_vtutotal);
  NEW_SUB_TIMER_STATIC(t_cells, "cells", t_vtutotal);
  NEW_SUB_TIMER_STATIC(t_variables, "variables", t_vtutotal);
  NEW_SUB_TIMER_STATIC(t_geometry, "geometry", t_vtutotal);
  RECORD_TIMER_START(t_vtutotal);
  RECORD_TIMER_START(t_faces);

  const auto noCells = (uint64_t)m_solver->m_extractedCells->size();
  auto noPoints = (uint64_t)m_solver->m_gridPoints->size();
  auto connSize = (uint64_t)(noNodes * noCells);
  auto facesSize = (uint64_t)0;
  vector<uint_t>* facestream = nullptr;
  vector<uint_t>* conn = nullptr;
  ScratchSpace<MInt> polyhedralCell(m_solver->a_noCells(), AT_, "polyhedralCell");
  MInt noPolyhedra = 0;
  polyhedralCell.fill(-1);
  const MBool getInternalPolyhedra = true;
  if(m_solver->m_vtuCutCellOutput) {
    connSize = vtuAssembleFaceStream(facestream, conn, facesSize, polyhedralCell, noPolyhedra, getInternalPolyhedra);
    if(facestream == nullptr || conn == nullptr) mTerm(1, AT_, "empty stream");
  }
  noPoints = (uint64_t)m_solver->m_gridPoints->size(); // update after cut cell assembly
  RECORD_TIMER_STOP(t_faces);
  RECORD_TIMER_START(t_init);


  const MBool mergeGridPoints = true;
  const uint64_t noPoints0 = noPoints;
  ScratchSpace<MInt> pointDomainId(noPoints, AT_, "pointDomainId");
  ScratchSpace<MInt> pointRemapId(noPoints, AT_, "pointRemapId");
  pointDomainId.fill(domainId());
  pointRemapId.fill(-1);
  if(mergeGridPoints) {
    MIntScratchSpace nghbrList(300, AT_, "nghbrList");
    MIntScratchSpace layerId(300, AT_, "layerId");
    ScratchSpace<MInt> cellMap(m_solver->a_noCells(), AT_, "cellMap");
    cellMap.fill(-1);
    MInt noSendDat0 = 0;
    MInt noRecvDat0 = 0;
    MInt noSendDatTotal = 0;
    for(MInt i = 0; i < m_solver->noNeighborDomains(); i++) {
      noSendDat0 += m_solver->noWindowCells(i);
      noRecvDat0 += m_solver->noHaloCells(i);
    }
    ScratchSpace<MInt> noCellPoints(noSendDat0 + 1, AT_, "noCellPoints");
    ScratchSpace<MInt> noCellPointsRecv(noRecvDat0 + 1, AT_, "noCellPointsRecv");
    ScratchSpace<MInt> noSendDat(noDomains(), AT_, "noSendDat");
    ScratchSpace<MInt> noRecvDat(noDomains(), AT_, "noRecvDat");
    noCellPoints.fill(0);
    noSendDat.fill(0);
    noRecvDat.fill(0);
    for(MInt c = 0; c < (signed)noCells; c++) {
      cellMap(m_solver->m_extractedCells->a[c].m_cellId) = c;
    }
    MInt scnt0 = 0;
    noSendDatTotal = 0;
    for(MInt i = 0; i < m_solver->noNeighborDomains(); i++) {
      noSendDat[i] = 0;
      for(MInt j = 0; j < m_solver->noWindowCells(i); j++) {
        MInt c = cellMap(m_solver->windowCellId(i, j));
        MInt pc = (c > -1) ? polyhedralCell[c] : -1;
        MInt np = 0;
        if(c > -1) np = (pc > -1) ? conn[pc].size() : noNodes;
        noCellPoints(scnt0) = np;
        noSendDatTotal += np;
        noSendDat[i] += np;
        scnt0++;
      }
    }
    ScratchSpace<MPI_Request> sendReq(m_solver->noNeighborDomains(), AT_, "sendReq");
    sendReq.fill(MPI_REQUEST_NULL);
    scnt0 = 0;
    for(MInt i = 0; i < m_solver->noNeighborDomains(); i++) {
      MPI_Issend(&(noCellPoints[scnt0]), m_solver->noWindowCells(i), MPI_INT, m_solver->neighborDomain(i), 132,
                 mpiComm(), &sendReq[i], AT_, "(noCellPoints[scnt0])");
      scnt0 += m_solver->noWindowCells(i);
    }
    MInt rcnt0 = 0;
    for(MInt i = 0; i < m_solver->noNeighborDomains(); i++) {
      MPI_Recv(&(noCellPointsRecv[rcnt0]), m_solver->noHaloCells(i), MPI_INT, m_solver->neighborDomain(i), 132,
               mpiComm(), MPI_STATUS_IGNORE, AT_, "(noCellPointsRecv[rcnt0])");
      rcnt0 += m_solver->noHaloCells(i);
    }
    if(m_solver->noNeighborDomains() > 0)
      MPI_Waitall(m_solver->noNeighborDomains(), &sendReq[0], MPI_STATUSES_IGNORE, AT_);
    MInt noRecvDatTotal = 0;
    rcnt0 = 0;
    for(MInt i = 0; i < m_solver->noNeighborDomains(); i++) {
      noRecvDat[i] = 0;
      for(MInt j = 0; j < m_solver->noHaloCells(i); j++) {
        noRecvDat[i] += noCellPointsRecv[rcnt0];
        noRecvDatTotal += noCellPointsRecv[rcnt0];
        rcnt0++;
      }
    }
    ScratchSpace<MInt> sendGridPoints(noSendDatTotal + 1, AT_, "sendGridPoints");
    ScratchSpace<MInt> recvGridPoints(noRecvDatTotal + 1, AT_, "recvGridPoints");
    ScratchSpace<MFloat> sendGridPointsCoord(noSendDatTotal + 1, 3, AT_, "sendGridPointsCoord");
    ScratchSpace<MFloat> recvGridPointsCoord(noRecvDatTotal + 1, 3, AT_, "recvGridPointsCoord");
    scnt0 = 0;
    MInt scnt = 0;
    for(MInt i = 0; i < m_solver->noNeighborDomains(); i++) {
      for(MInt j = 0; j < m_solver->noWindowCells(i); j++) {
        MInt c = cellMap(m_solver->windowCellId(i, j));
        if(noCellPoints(scnt0) > 0) {
          if(c < 0) {
            cerr << "case should not occur" << endl;
            continue;
          }
          MInt pc = polyhedralCell[c];
          MInt noGridPoints = (pc > -1) ? conn[pc].size() : noNodes;
          for(MInt p = 0; p < noGridPoints; p++) {
            MInt gridPointId = (pc > -1) ? conn[pc][p] : m_solver->m_extractedCells->a[c].m_pointIds[p];
            if(gridPointId < 0) cerr << "negative grid point id " << endl;
            sendGridPoints(scnt) = gridPointId;
            for(MInt k = 0; k < nDim; k++) {
              sendGridPointsCoord(scnt, k) = m_solver->m_gridPoints->a[gridPointId].m_coordinates[k];
            }
            scnt++;
          }
        }
        scnt0++;
      }
    }

    sendReq.fill(MPI_REQUEST_NULL);
    scnt = 0;
    for(MInt i = 0; i < m_solver->noNeighborDomains(); i++) {
      MPI_Issend(&(sendGridPoints[scnt]), noSendDat[i], MPI_INT, m_solver->neighborDomain(i), 133, mpiComm(),
                 &sendReq[i], AT_, "(sendGridPoints[scnt])");
      scnt += noSendDat[i];
    }
    MInt rcnt = 0;
    for(MInt i = 0; i < m_solver->noNeighborDomains(); i++) {
      MPI_Recv(&(recvGridPoints[rcnt]), noRecvDat[i], MPI_INT, m_solver->neighborDomain(i), 133, mpiComm(),
               MPI_STATUS_IGNORE, AT_, "(recvGridPoints[rcnt])");
      rcnt += noRecvDat[i];
    }
    if(m_solver->noNeighborDomains() > 0)
      MPI_Waitall(m_solver->noNeighborDomains(), &sendReq[0], MPI_STATUSES_IGNORE, AT_);

    sendReq.fill(MPI_REQUEST_NULL);
    scnt = 0;
    for(MInt i = 0; i < m_solver->noNeighborDomains(); i++) {
      MPI_Issend(&(sendGridPointsCoord[scnt]), 3 * noSendDat[i], MPI_DOUBLE, m_solver->neighborDomain(i), 134,
                 mpiComm(), &sendReq[i], AT_, "(sendGridPointsCoord[scnt])");
      scnt += 3 * noSendDat[i];
    }
    rcnt = 0;
    for(MInt i = 0; i < m_solver->noNeighborDomains(); i++) {
      MPI_Recv(&(recvGridPointsCoord[rcnt]), 3 * noRecvDat[i], MPI_DOUBLE, m_solver->neighborDomain(i), 134, mpiComm(),
               MPI_STATUS_IGNORE, AT_, "(recvGridPointsCoord[rcnt])");
      rcnt += 3 * noRecvDat[i];
    }
    if(m_solver->noNeighborDomains() > 0)
      MPI_Waitall(m_solver->noNeighborDomains(), &sendReq[0], MPI_STATUSES_IGNORE, AT_);

    rcnt0 = 0;
    rcnt = 0;
    map<MInt, MInt> haloMap;
    for(MInt i = 0; i < m_solver->noNeighborDomains(); i++) {
      for(MInt j = 0; j < m_solver->noHaloCells(i); j++) {
        MInt haloId = m_solver->haloCellId(i, j);
        MInt noRecvPoints = noCellPointsRecv[rcnt0];
        if(noRecvPoints > 0) {
          // const MInt counter = m_solver->getAdjacentLeafCells<2>( haloId, 2, nghbrList, layerId );
          const MInt counter = m_solver->getAdjacentLeafCells_d2(haloId, 2, nghbrList, layerId);
          for(MInt n = 0; n < counter; n++) {
            MInt nghbrId = nghbrList[n];
            if(nghbrId < 0) continue;
            if(m_solver->a_isHalo(nghbrId)) continue;
            MInt ncell = cellMap[nghbrId];
            if(ncell < 0) continue;
            MInt pc = polyhedralCell[ncell];
            MFloat delta0 = (pc > -1) ? 1e-12 * m_solver->c_cellLengthAtCell(nghbrId)
                                      : 0.0001 * m_solver->c_cellLengthAtCell(nghbrId);
            MInt noPointsNghbr = (pc > -1) ? conn[pc].size() : noNodes;
            for(MInt p = 0; p < noRecvPoints; p++) {
              for(MInt q = 0; q < noPointsNghbr; q++) {
                MInt gridPointId = (pc > -1) ? conn[pc][q] : m_solver->m_extractedCells->a[ncell].m_pointIds[q];
                if(gridPointId < 0) cerr << "negative grid point id " << endl;
                MFloat delta = sqrt(
                    POW2(m_solver->m_gridPoints->a[gridPointId].m_coordinates[0] - recvGridPointsCoord(rcnt + p, 0))
                    + POW2(m_solver->m_gridPoints->a[gridPointId].m_coordinates[1] - recvGridPointsCoord(rcnt + p, 1))
                    + POW2(m_solver->m_gridPoints->a[gridPointId].m_coordinates[2] - recvGridPointsCoord(rcnt + p, 2)));
                if(delta < delta0) {
                  if(m_solver->neighborDomain(i) < pointDomainId[gridPointId]) {
                    pointDomainId[gridPointId] = m_solver->neighborDomain(i);
                    pointRemapId[gridPointId] = recvGridPoints[rcnt + p];
                    haloMap[gridPointId] = rcnt + p;
                  }
                }
              }
            }
          }
        }
        rcnt += noRecvPoints;
        rcnt0++;
      }
    }

    // update number of points
    noPoints = 0;
    for(uint64_t p = 0; p < noPoints0; p++) {
      if(pointRemapId[p] < 0) {
        pointRemapId[p] = noPoints;
        if(pointDomainId[p] != domainId()) cerr << "unexpected case" << endl;
        noPoints++;
      }
      if(pointRemapId[p] < 0) {
        cerr << "unexpected case 2" << endl;
      }
    }


    // communicate point shifts
    scnt0 = 0;
    scnt = 0;
    for(MInt i = 0; i < m_solver->noNeighborDomains(); i++) {
      for(MInt j = 0; j < m_solver->noWindowCells(i); j++) {
        MInt c = cellMap(m_solver->windowCellId(i, j));
        if(noCellPoints(scnt0) > 0) {
          if(c < 0) {
            cerr << "case should not occur" << endl;
            continue;
          }
          MInt pc = polyhedralCell[c];
          MInt noGridPoints = (pc > -1) ? conn[pc].size() : noNodes;
          for(MInt p = 0; p < noGridPoints; p++) {
            MInt gridPointId = (pc > -1) ? conn[pc][p] : m_solver->m_extractedCells->a[c].m_pointIds[p];
            if(gridPointId < 0) cerr << "negative grid point id " << endl;
            sendGridPoints(scnt) = (pointDomainId[gridPointId] == domainId()) ? pointRemapId[gridPointId] : -1;
            scnt++;
          }
        }
        scnt0++;
      }
    }
    sendReq.fill(MPI_REQUEST_NULL);
    scnt = 0;
    for(MInt i = 0; i < m_solver->noNeighborDomains(); i++) {
      MPI_Issend(&(sendGridPoints[scnt]), noSendDat[i], MPI_INT, m_solver->neighborDomain(i), 135, mpiComm(),
                 &sendReq[i], AT_, "(sendGridPoints[scnt])");
      scnt += noSendDat[i];
    }
    rcnt = 0;
    for(MInt i = 0; i < m_solver->noNeighborDomains(); i++) {
      MPI_Recv(&(recvGridPoints[rcnt]), noRecvDat[i], MPI_INT, m_solver->neighborDomain(i), 135, mpiComm(),
               MPI_STATUS_IGNORE, AT_, "(recvGridPoints[rcnt])");
      rcnt += noRecvDat[i];
    }
    if(m_solver->noNeighborDomains() > 0)
      MPI_Waitall(m_solver->noNeighborDomains(), &sendReq[0], MPI_STATUSES_IGNORE, AT_);

    for(auto& it : haloMap) {
      MInt gridPointId = it.first;
      MInt otherId = recvGridPoints[it.second];
      if(gridPointId < 0 || gridPointId >= (signed)noPoints0) cerr << "point out of range" << endl;
      pointRemapId[gridPointId] = otherId;
      if(otherId < 0 || (pointDomainId[gridPointId] >= domainId())) {
        cerr << "inconsistent point shifts " << gridPointId << " " << otherId << " " << pointDomainId[gridPointId]
             << " " << domainId() << endl;
        continue;
      }
    }


  } else {
    for(uint64_t p = 0; p < noPoints; p++) {
      pointRemapId[p] = p;
    }
  }

  ScratchSpace<uint64_t> noPointsPerDomain(noDomains(), AT_, "noPointsPerDomain");
  ScratchSpace<uint64_t> pointOffsetsGlobal(noDomains() + 1, AT_, "pointOffsetsGlobal");
  ScratchSpace<uint64_t> noCellsPerDomain(noDomains(), AT_, "noCellsPerDomain");
  ScratchSpace<uint64_t> connSizePerDomain(noDomains(), AT_, "connSizePerDomain");
  ScratchSpace<uint64_t> facesPerDomain(noDomains(), AT_, "facesPerDomain");
  ScratchSpace<uint64_t> dataPerDomain(noDomains(), 4, AT_, "dataPerDomain");
  dataPerDomain(domainId(), 0) = noPoints;
  dataPerDomain(domainId(), 1) = noCells;
  dataPerDomain(domainId(), 2) = connSize;
  dataPerDomain(domainId(), 3) = facesSize;
  MPI_Allgather(MPI_IN_PLACE, 4, MPI_UINT64_T, &dataPerDomain[0], 4, MPI_UINT64_T, mpiComm(), AT_, "MPI_IN_PLACE",
                "dataPerDomain[0]");
  pointOffsetsGlobal[0] = 0;
  for(MInt d = 0; d < noDomains(); d++) {
    pointOffsetsGlobal[d + 1] = pointOffsetsGlobal[d] + dataPerDomain(d, 0);
    noPointsPerDomain[d] = dataPerDomain(d, 0);
    noCellsPerDomain[d] = dataPerDomain(d, 1);
    connSizePerDomain[d] = dataPerDomain(d, 2);
    facesPerDomain[d] = dataPerDomain(d, 3);
  }

  uint64_t globalNoPoints = 0;
  uint64_t globalNoCells = 0;
  uint64_t globalConnSize = 0;
  uint64_t globalFacesSize = 0;
  uint64_t globalPointOffset = 0;
  uint64_t globalCellOffset = 0;
  uint64_t globalConnSizeOffset = 0;
  uint64_t globalFacesOffset = 0;
  uint64_t offset = 0;
  for(MInt d = 0; d < noDomains(); d++) {
    globalNoPoints += noPointsPerDomain[d];
    globalNoCells += noCellsPerDomain[d];
    globalConnSize += connSizePerDomain[d];
    globalFacesSize += facesPerDomain[d];
  }
  for(MInt d = 0; d < domainId(); d++) {
    globalPointOffset += noPointsPerDomain[d];
    globalCellOffset += noCellsPerDomain[d];
    globalConnSizeOffset += connSizePerDomain[d];
    globalFacesOffset += facesPerDomain[d];
  }

  // if output size exceeds 32bit boundary rerun with 64 bit data types
  if(std::is_same<uint_t, uint32_t>::value) {
    if(globalNoCells > (uint64_t)std::numeric_limits<uint_t>::max()
       || globalNoPoints > (uint64_t)std::numeric_limits<uint_t>::max()) {
      return writeVtuOutputParallel<uint64_t>(fileName, geomFileName);
    }
  }

  RECORD_TIMER_STOP(t_init);
  RECORD_TIMER_START(t_prepare);

  if(string2enum(m_solver->m_outputFormat) == VTU) {
    //-----------
    // points
    uint64_t pointsMemsize = 3 * noPoints * sizeof(float);
    uint64_t pointsMemsizeGlobal = 3 * globalNoPoints * sizeof(float);
    uint64_t pointsOffset = 3 * globalPointOffset * sizeof(float);
    const MInt pointsMaxSize = estimateMemorySizeSolverwise(noPoints, noPointsPerDomain, 3 * sizeof(float));
    ScratchSpace<float> points(pointsMaxSize, 3, AT_, "points");
    uint_t cnt = 0;
    for(uint_t p = 0; p < noPoints0; p++) {
      if(pointDomainId[p] != domainId()) continue;
      for(MInt i = 0; i < nDim; i++) {
        points((MInt)cnt, i) = (float)m_solver->m_gridPoints->a[p].m_coordinates[i];
      }
      IF_CONSTEXPR(nDim == 2) points((MInt)cnt, 2) = (float)0.0;
      cnt++;
    }
    insertDataHeader(reinterpret_cast<char*>(&points(0)), pointsMemsize, pointsMemsizeGlobal, pointsOffset);

    //-----------
    // connectivity
    uint64_t connectivityMemsize = connSize * sizeof(uint_t);
    uint64_t connectivityOffset = globalConnSizeOffset * sizeof(uint_t);
    uint64_t connectivityMemsizeGlobal = globalConnSize * sizeof(uint_t);
    const MInt connectivityMaxSize = estimateMemorySizeSolverwise(noCells, noCellsPerDomain, noNodes * sizeof(uint_t));
    ScratchSpace<uint_t> connectivity(mMax(connectivityMaxSize * noNodes, connSize + 20), AT_, "connectivity");
    cnt = 0;
    for(MInt c = 0; c < (signed)noCells; c++) {
      MInt pc = polyhedralCell[c];
      if(pc < 0) {
        for(MInt p = 0; p < (signed)noNodes; p++) {
          if(m_solver->m_extractedCells->a[c].m_pointIds[p] < 0)
            cerr << domainId() << ": warning negative point " << c << " " << m_solver->m_extractedCells->a[c].m_cellId
                 << " " << m_solver->m_extractedCells->a[c].m_pointIds[p] << endl;
          MInt gridPointId = m_solver->m_extractedCells->a[c].m_pointIds[p];
          connectivity(cnt) =
              (uint_t)pointOffsetsGlobal[pointDomainId[gridPointId]] + (uint_t)pointRemapId[gridPointId];
          cnt++;
        }
      } else {
        ASSERT(cnt + conn[pc].size() <= connectivity.size(), "");
        for(MUint i = 0; i < conn[pc].size(); i++) {
          uint_t gridPointId = conn[pc][i];
          connectivity(cnt) =
              (uint_t)pointOffsetsGlobal[pointDomainId[gridPointId]] + (uint_t)pointRemapId[gridPointId];
          cnt++;
        }
      }
    }
    for(MUint i = 0; i < connSize; i++) {
      if(connectivity[i] >= globalNoPoints) {
        cerr << domainId() << ": invalid grid point " << i << " " << connectivity[i] << " " << globalPointOffset << " "
             << noPoints << " " << globalNoPoints << " " << globalTimeStep << endl;
        connectivity[i] = 0;
      }
    }
    insertDataHeader(reinterpret_cast<char*>(&connectivity(0)), connectivityMemsize, connectivityMemsizeGlobal,
                     connectivityOffset);

    //-----------
    // offsets
    uint64_t offsetsMemsize = noCells * sizeof(uint64_t);
    uint64_t offsetsOffset = globalCellOffset * sizeof(uint64_t);
    uint64_t offsetsMemsizeGlobal = globalNoCells * sizeof(uint64_t);
    const MInt offsetsMaxSize = estimateMemorySizeSolverwise(noCells, noCellsPerDomain, sizeof(uint64_t));
    ScratchSpace<uint64_t> offsets(offsetsMaxSize, AT_, "offsets");
    uint64_t offs = globalConnSizeOffset;
    for(uint64_t c = 0; c < noCells; c++) {
      MInt pc = polyhedralCell[c];
      if(pc < 0)
        offs += (uint64_t)noNodes;
      else
        offs += (uint64_t)conn[pc].size();
      offsets(c) = offs;
      if(offsets(c) - globalConnSizeOffset > connSize || offsets(c) > globalConnSize) {
        cerr << domainId() << ": invalid offset: " << offsets(c) << " " << connSize << " " << globalConnSizeOffset
             << " " << globalConnSize << endl;
      }
    }
    insertDataHeader(reinterpret_cast<char*>(&offsets(0)), offsetsMemsize, offsetsMemsizeGlobal, offsetsOffset);

    //-----------
    // types
    uint64_t typesMemsize = noCells * sizeof(uint8_t);
    uint64_t typesOffset = globalCellOffset * sizeof(uint8_t);
    uint64_t typesMemsizeGlobal = globalNoCells * sizeof(uint8_t);
    const MInt typesMaxSize = estimateMemorySizeSolverwise(noCells, noCellsPerDomain, sizeof(uint8_t));
    ScratchSpace<uint8_t> types(typesMaxSize, AT_, "types");
    for(MInt c = 0; c < (signed)noCells; c++) {
      MInt pc = polyhedralCell[c];
      types(c) = (pc < 0) ? baseCellType : polyCellType;
    }
    insertDataHeader(reinterpret_cast<char*>(&types(0)), typesMemsize, typesMemsizeGlobal, typesOffset);

    //-----------
    // faces
    uint64_t facesMemsize = facesSize * sizeof(uint_t);
    uint64_t facesOffset = globalFacesOffset * sizeof(uint_t);
    uint64_t facesMemsizeGlobal = globalFacesSize * sizeof(uint_t);
    const MInt facesMaxSize = estimateMemorySizeSolverwise(facesSize, facesPerDomain, sizeof(uint_t));
    ScratchSpace<uint_t> faces(facesMaxSize, AT_, "faces");
    cnt = 0;
    for(MInt c = 0; c < (signed)noCells; c++) {
      MInt pc = polyhedralCell[c];
      if(pc < 0) {
        faces(cnt) = 0; // can skip?
        cnt++;          // can skip?
      } else {
        if(facestream[pc].empty()) mTerm(1, AT_, "empty face stream");
        if(conn[pc].empty()) mTerm(1, AT_, "empty conn stream");
        // copy( facestream[pc].begin(), facestream[pc].end(), &faces(cnt) );
        MInt fcnt = 0;
        const MInt noFaces = facestream[pc][fcnt];
        faces(cnt) = facestream[pc][fcnt];
        cnt++;
        fcnt++;
        for(MInt face = 0; face < noFaces; face++) {
          MInt pointCnt = facestream[pc][fcnt];
          faces(cnt) = facestream[pc][fcnt];
          cnt++;
          fcnt++;
          for(MInt p = 0; p < pointCnt; p++) {
            uint_t gridPointId = facestream[pc][fcnt];
            faces(cnt) = (uint_t)pointOffsetsGlobal[pointDomainId[gridPointId]] + (uint_t)pointRemapId[gridPointId];
            cnt++;
            fcnt++;
          }
        }
      }
    }

    insertDataHeader(reinterpret_cast<char*>(&faces(0)), facesMemsize, facesMemsizeGlobal, facesOffset);

    //-----------
    // faceoffsets
    uint64_t faceoffsetsMemsize = noCells * sizeof(uint64_t);
    uint64_t faceoffsetsOffset = globalCellOffset * sizeof(uint64_t);
    uint64_t faceoffsetsMemsizeGlobal = globalNoCells * sizeof(uint64_t);
    const MInt faceoffsetsMaxSize = estimateMemorySizeSolverwise(noCells, noCellsPerDomain, sizeof(uint64_t));
    ScratchSpace<uint64_t> faceoffsets(faceoffsetsMaxSize, AT_, "faceoffsets");
    offs = globalFacesOffset;
    for(MInt c = 0; c < (signed)noCells; c++) {
      MInt pc = polyhedralCell[c];
      if(pc < 0)
        offs++;
      else
        offs += (uint64_t)facestream[pc].size();
      faceoffsets(c) = offs;
    }
    insertDataHeader(reinterpret_cast<char*>(&faceoffsets(0)), faceoffsetsMemsize, faceoffsetsMemsizeGlobal,
                     faceoffsetsOffset);

    //-----------
    // globalId
    uint64_t globalIdMemsize = m_solver->m_vtuGlobalIdOutput ? noCells * sizeof(uint_t) : 0;
    uint64_t globalIdOffset = globalCellOffset * sizeof(uint_t);
    uint64_t globalIdMemsizeGlobal = globalNoCells * sizeof(uint_t);
    const MInt globalIdMaxSize = estimateMemorySizeSolverwise(noCells, noCellsPerDomain, sizeof(uint_t));
    ScratchSpace<uint_t> globalId(globalIdMaxSize, AT_, "globalId");
    if(m_solver->m_vtuGlobalIdOutput) {
      for(uint_t c = 0; c < noCells; c++) {
        globalId(c) = (uint_t)m_solver->c_globalId(
            m_solver->a_hasProperty(m_solver->m_extractedCells->a[c].m_cellId, SolverCell::IsSplitChild)
                ? m_solver->getAssociatedInternalCell(m_solver->m_extractedCells->a[c].m_cellId)
                : m_solver->m_extractedCells->a[c].m_cellId);
      }
      insertDataHeader(reinterpret_cast<char*>(&globalId(0)), globalIdMemsize, globalIdMemsizeGlobal, globalIdOffset);
    }

    //-----------
    // domainIds
    uint64_t domainIdsMemsize = m_solver->m_vtuDomainIdOutput ? noCells * sizeof(uint_t) : 0;
    uint64_t domainIdsOffset = globalCellOffset * sizeof(uint_t);
    uint64_t domainIdsMemsizeGlobal = globalNoCells * sizeof(uint_t);
    const MInt domainIdsMaxSize = estimateMemorySizeSolverwise(noCells, noCellsPerDomain, sizeof(uint_t));
    ScratchSpace<uint_t> domainIds(domainIdsMaxSize, AT_, "domainIds");
    if(m_solver->m_vtuGlobalIdOutput) {
      for(uint_t c = 0; c < noCells; c++) {
        domainIds(c) = (uint_t)domainId();
      }

      insertDataHeader(reinterpret_cast<char*>(&domainIds(0)), domainIdsMemsize, domainIdsMemsizeGlobal,
                       domainIdsOffset);
    }

    //-----------
    // variables
    const uint64_t dataSize = m_solver->m_vtuWritePointData ? noPoints : noCells;
    const uint64_t globalDataSize = m_solver->m_vtuWritePointData ? globalNoPoints : globalNoCells;
    const uint64_t globalDataOffset = m_solver->m_vtuWritePointData ? globalPointOffset : globalCellOffset;
    ScratchSpace<uint64_t>* noDataPerDomain = m_solver->m_vtuWritePointData ? &noPointsPerDomain : &noCellsPerDomain;
    uint64_t velocityMemsize = 3 * dataSize * sizeof(float);
    uint64_t vorticityMemsize = m_solver->m_vtuVorticityOutput ? 3 * dataSize * sizeof(float) : 0;
    uint64_t velGradMemsize = m_solver->m_vtuVelocityGradientOutput ? 9 * dataSize * sizeof(float) : 0;
    uint64_t pressureMemsize = dataSize * sizeof(float);
    uint64_t densityMemsize = m_solver->m_vtuDensityOutput ? dataSize * sizeof(float) : 0;
    uint64_t progressMemsize = dataSize * m_solver->m_noSpecies * sizeof(float);
    uint64_t levelSetMemsize = m_solver->m_vtuLevelSetOutput ? dataSize * sizeof(float) : 0;
    uint64_t QMemsize = m_solver->m_vtuQCriterionOutput ? dataSize * sizeof(float) : 0;
    uint64_t L2Memsize = m_solver->m_vtuLambda2Output ? dataSize * sizeof(float) : 0;
    uint64_t solverSpecificVarsMemsize = noSolverSpecificVars > 0 ? dataSize * sizeof(float) * noSolverSpecificVars : 0;
    uint64_t userVarsMemsize = noUserVars > 0 ? dataSize * sizeof(float) * noUserVars : 0;
    uint64_t velocityOffset = 3 * globalDataOffset * sizeof(float);
    uint64_t vorticityOffset = m_solver->m_vtuVorticityOutput ? 3 * globalDataOffset * sizeof(float) : 0;
    uint64_t velGradOffset = m_solver->m_vtuVelocityGradientOutput ? 9 * globalDataOffset * sizeof(float) : 0;
    uint64_t pressureOffset = globalDataOffset * sizeof(float);
    uint64_t densityOffset = m_solver->m_vtuDensityOutput ? globalDataOffset * sizeof(float) : 0;
    uint64_t progressOffset = globalDataOffset * m_solver->m_noSpecies * sizeof(float);
    uint64_t levelSetOffset = m_solver->m_vtuLevelSetOutput ? globalDataOffset * sizeof(float) : 0;
    uint64_t QOffset = m_solver->m_vtuQCriterionOutput ? globalDataOffset * sizeof(float) : 0;
    uint64_t L2Offset = m_solver->m_vtuLambda2Output ? globalDataOffset * sizeof(float) : 0;
    uint64_t solverSpecificVarsOffset =
        noSolverSpecificVars > 0 ? globalDataOffset * sizeof(float) * noSolverSpecificVars : 0;
    uint64_t userVarsOffset = noUserVars > 0 ? globalDataOffset * sizeof(float) * noUserVars : 0;
    uint64_t userVarsMemsizeGlobal = noUserVars > 0 ? globalDataSize * noUserVars * sizeof(float) : 0;
    uint64_t solverSpecificVarsMemsizeGlobal =
        noSolverSpecificVars > 0 ? globalDataSize * noSolverSpecificVars * sizeof(float) : 0;
    uint64_t velocityMemsizeGlobal = 3 * globalDataSize * sizeof(float);
    uint64_t vorticityMemsizeGlobal = m_solver->m_vtuVorticityOutput ? 3 * globalDataSize * sizeof(float) : 0;
    uint64_t velGradMemsizeGlobal = m_solver->m_vtuVelocityGradientOutput ? 9 * globalDataSize * sizeof(float) : 0;
    uint64_t pressureMemsizeGlobal = globalDataSize * sizeof(float);
    uint64_t densityMemsizeGlobal = m_solver->m_vtuDensityOutput ? globalDataSize * sizeof(float) : 0;
    uint64_t progressMemsizeGlobal = globalDataSize * m_solver->m_noSpecies * sizeof(float);
    uint64_t levelSetMemsizeGlobal = m_solver->m_vtuLevelSetOutput ? globalDataSize * sizeof(float) : 0;
    uint64_t QMemsizeGlobal = m_solver->m_vtuQCriterionOutput ? globalDataSize * sizeof(float) : 0;
    uint64_t L2MemsizeGlobal = m_solver->m_vtuLambda2Output ? globalDataSize * sizeof(float) : 0;
    const MInt vectorMaxSize = estimateMemorySizeSolverwise(dataSize, *noDataPerDomain, 3 * sizeof(float));
    const MInt scalarMaxSize = estimateMemorySizeSolverwise(dataSize, *noDataPerDomain, sizeof(float));
    ScratchSpace<float> velocity(vectorMaxSize, 3, AT_, "velocity");
    ScratchSpace<float> vorticity(m_solver->m_vtuVorticityOutput ? vectorMaxSize : 0, 3, AT_, "vorticity");
    ScratchSpace<float> velGrad(m_solver->m_vtuVelocityGradientOutput ? vectorMaxSize : 0, 9, AT_, "velGrad");
    ScratchSpace<float> pressure(scalarMaxSize, AT_, "pressure");
    ScratchSpace<float> density(m_solver->m_vtuDensityOutput ? scalarMaxSize : 0, AT_, "density");
    ScratchSpace<float> progress(scalarMaxSize * m_solver->m_noSpecies, AT_, "progress");
    ScratchSpace<float> levelSet(m_solver->m_vtuLevelSetOutput ? scalarMaxSize : 0, AT_, "levelSet");
    ScratchSpace<float> Q(m_solver->m_vtuQCriterionOutput ? scalarMaxSize : 0, AT_, "Q");
    ScratchSpace<float> lambda2(m_solver->m_vtuLambda2Output ? scalarMaxSize : 0, AT_, "lambda2");
    ScratchSpace<float> solverSpecificVarsOut(noSolverSpecificVars > 0 ? scalarMaxSize * noSolverSpecificVars : 0, AT_,
                                              "solverSpecificVars");
    ScratchSpace<float> userVarsOut(noUserVars > 0 ? scalarMaxSize * noUserVars : 0, AT_, "userVars");
    if(m_solver->m_vtuWritePointData) {
      for(MInt p = 0; p < (signed)noPoints; p++) {
        for(MInt i = 0; i < 3; i++)
          velocity(p, i) = (float)0.0;
        for(MInt i = 0; i < noUserVars; i++)
          userVarsOut[p * noUserVars + i] = (float)0.0;
        for(MInt i = 0; i < noSolverSpecificVars; i++)
          solverSpecificVarsOut[p * noSolverSpecificVars + i] = (float)0.0;
        if(m_solver->m_vtuVorticityOutput)
          for(MInt i = 0; i < nDim; i++)
            vorticity(p, i) = (float)0.0;
        if(m_solver->m_vtuVelocityGradientOutput)
          for(MInt i = 0; i < nDim * nDim; i++)
            velGrad(p, i) = (float)0.0;
        pressure(p) = (float)0.0;
        if(m_solver->m_vtuDensityOutput) density(p) = (float)0.0;
        if(m_solver->m_noSpecies > 0) progress(p) = (float)0.0;
        if(m_solver->m_vtuLevelSetOutput) levelSet(p) = (float)0.0;
        if(m_solver->m_vtuQCriterionOutput) Q(p) = (float)0.0;
        if(m_solver->m_vtuLambda2Output) lambda2(p) = (float)0.0;
        MFloat sum = F0;
        for(MInt n = 0; n < m_solver->m_gridPoints->a[p].m_noAdjacentCells; n++) {
          MInt cellId = m_solver->m_gridPoints->a[p].m_cellIds[n];
          MFloat dx = F0;
          for(MInt i = 0; i < nDim; i++) {
            dx += POW2(m_solver->a_coordinate(cellId, i) - m_solver->m_gridPoints->a[p].m_coordinates[i]);
          }
          dx = mMax(m_solver->m_eps, sqrt(dx));
          if(m_solver->m_vtuQCriterionOutput) {
            float q = F0;
            IF_CONSTEXPR(nDim == 3) {
              for(MInt i = 0; i < nDim; i++)
                for(MInt j = 0; j < nDim; j++)
                  q += POW2(F1B2
                            * (m_solver->a_slope(cellId, sysEqn().PV->VV[i], j)
                               - m_solver->a_slope(cellId, sysEqn().PV->VV[j], i)))
                       - POW2(F1B2
                              * (m_solver->a_slope(cellId, sysEqn().PV->VV[i], j)
                                 + m_solver->a_slope(cellId, sysEqn().PV->VV[j], i)));
              q *= F1B2;
            }
            else {
              q += dx * F1B2
                   * (m_solver->a_slope(cellId, sysEqn().PV->VV[1], 0)
                      - m_solver->a_slope(cellId, sysEqn().PV->VV[0], 1));
            }
            Q(p) += dx * q;
          }
          for(MInt i = 0; i < nDim; i++)
            velocity(p, i) += dx * m_solver->a_pvariable(cellId, sysEqn().PV->VV[i]);
          if(m_solver->m_vtuVorticityOutput) {
            ASSERT(nDim == 3, "");
            vorticity(p, 0) +=
                dx * F1B2
                * (m_solver->a_slope(cellId, sysEqn().PV->VV[2], 1) - m_solver->a_slope(cellId, sysEqn().PV->VV[1], 2));
            vorticity(p, 1) +=
                dx * F1B2
                * (m_solver->a_slope(cellId, sysEqn().PV->VV[0], 2) - m_solver->a_slope(cellId, sysEqn().PV->VV[2], 0));
            vorticity(p, 2) +=
                dx * F1B2
                * (m_solver->a_slope(cellId, sysEqn().PV->VV[1], 0) - m_solver->a_slope(cellId, sysEqn().PV->VV[0], 1));
          }
          if(m_solver->m_vtuVelocityGradientOutput) {
            for(MInt i = 0; i < nDim; i++) {
              for(MInt j = 0; j < nDim; j++) {
                velGrad(p, i * nDim + j) += dx * m_solver->a_slope(cellId, sysEqn().PV->VV[i], j);
              }
            }
          }
          if(m_solver->m_vtuLambda2Output) {
            MFloat vGrad[3][3]{};
            MFloat O[3][3]{};
            MFloat S[3][3]{};
            MFloat eigenVal[3]{};
            for(MInt i = 0; i < nDim; i++) {
              for(MInt j = 0; j < nDim; j++) {
                O[i][j] = F1B2
                          * (m_solver->a_slope(cellId, sysEqn().PV->VV[i], j)
                             - m_solver->a_slope(cellId, sysEqn().PV->VV[j], i))
                          / m_solver->m_UInfinity;
                S[i][j] = F1B2
                          * (m_solver->a_slope(cellId, sysEqn().PV->VV[i], j)
                             + m_solver->a_slope(cellId, sysEqn().PV->VV[j], i))
                          / m_solver->m_UInfinity;
              }
            }
            for(MInt i = 0; i < nDim; i++) {
              for(MInt j = 0; j < nDim; j++) {
                vGrad[i][j] = F0;
                for(MInt k = 0; k < nDim; k++) {
                  vGrad[i][j] += O[i][k] * O[k][j];
                  vGrad[i][j] += S[i][k] * S[k][j];
                }
              }
            }
            maia::math::calcEigenValues(vGrad, eigenVal);
            sort(eigenVal, eigenVal + 3);
            lambda2(p) += dx * eigenVal[1];
          }
          pressure(p) += dx * m_solver->a_pvariable(cellId, sysEqn().PV->P);
          if(m_solver->m_vtuDensityOutput) density(p) += dx * m_solver->a_pvariable(cellId, sysEqn().PV->RHO);
          IF_CONSTEXPR(hasPV_C<SysEqn>::value)
          if(m_solver->m_noSpecies > 0) progress(p) += dx * m_solver->a_pvariable(cellId, sysEqn().PV->C);
          if(m_solver->m_vtuLevelSetOutput) {
            ASSERT(m_solver->m_levelSet, "");
            levelSet(p) += dx * m_solver->a_levelSetFunction(cellId, 0);
          }
          sum += dx;
        }
        for(MInt i = 0; i < nDim; i++)
          velocity(p, i) /= sum;
        if(m_solver->m_vtuVorticityOutput)
          for(MInt i = 0; i < nDim; i++)
            vorticity(p, i) /= sum;
        if(m_solver->m_vtuVelocityGradientOutput)
          for(MInt i = 0; i < nDim * nDim; i++)
            velGrad(p, i) /= sum;
        pressure(p) /= sum;
        if(m_solver->m_vtuDensityOutput) density(p) /= sum;
        if(m_solver->m_noSpecies > 0) progress(p) /= sum;
        if(m_solver->m_vtuQCriterionOutput) Q(p) /= sum;
        if(m_solver->m_vtuLambda2Output) lambda2(p) /= sum;
      }
    } else {
      for(MInt c = 0; c < (signed)noCells; c++) {
        MInt cellId = m_solver->m_extractedCells->a[c].m_cellId;
        for(MInt i = 0; i < nDim; i++)
          velocity(c, i) = m_solver->a_pvariable(cellId, sysEqn().PV->VV[i]);
        for(MInt i = 0; i < noUserVars; i++)
          userVarsOut[c * noUserVars + i] = userVars[noUserVars * cellId + i];
        for(MInt i = 0; i < noSolverSpecificVars; i++)
          solverSpecificVarsOut[c * noSolverSpecificVars + i] = solverSpecificVars[noSolverSpecificVars * cellId + i];
        pressure(c) = m_solver->a_pvariable(cellId, sysEqn().PV->P);
        if(m_solver->m_vtuDensityOutput) density(c) = m_solver->a_pvariable(cellId, sysEqn().PV->RHO);
        IF_CONSTEXPR(hasPV_C<SysEqn>::value)
        if(m_solver->m_noSpecies > 0) progress(c) = m_solver->a_pvariable(cellId, sysEqn().PV->C);
        if(m_solver->m_vtuLevelSetOutput) {
          ASSERT(!(m_solver->m_levelSetMb && !m_solver->m_constructGField), "");
          if(m_solver->m_combustion) {
            levelSet(c) = m_solver->a_levelSetFunction(cellId, 0);
          } else if(m_solver->m_levelSetMb && m_solver->m_constructGField) {
            if(noSolverSpecificVars > 0) {
              levelSet(c) = solverSpecificVars[noSolverSpecificVars * cellId];
            } else {
              ASSERT(false, "levelSet output failed for levelSetMb with constructGField.");
              levelSet(c) = F0; // Set 0 as default without message or termination to avoid, if no solverSpecificVars
                                // not provided.
            }
          } else {
            ASSERT(m_solver->m_levelSet, "");
            levelSet(c) = m_solver->a_levelSetFunction(cellId, 0);
          }
        }
        if(m_solver->m_vtuQCriterionOutput) {
          Q(c) = F0;
          IF_CONSTEXPR(nDim == 3) {
            for(MInt i = 0; i < nDim; i++)
              for(MInt j = 0; j < nDim; j++)
                Q(c) += POW2(F1B2
                             * (m_solver->a_slope(cellId, sysEqn().PV->VV[i], j)
                                - m_solver->a_slope(cellId, sysEqn().PV->VV[j], i)))
                        - POW2(F1B2
                               * (m_solver->a_slope(cellId, sysEqn().PV->VV[i], j)
                                  + m_solver->a_slope(cellId, sysEqn().PV->VV[j], i)));
            Q(c) *= F1B2;
          }
          else {
            Q(c) =
                F1B2
                * (m_solver->a_slope(cellId, sysEqn().PV->VV[1], 0) - m_solver->a_slope(cellId, sysEqn().PV->VV[0], 1));
          }
        }
        if(m_solver->m_vtuLambda2Output) {
          MFloat vGrad[3][3];
          MFloat O[3][3];
          MFloat S[3][3];
          MFloat eigenVal[3];
          for(MInt i = 0; i < nDim; i++) {
            for(MInt j = 0; j < nDim; j++) {
              O[i][j] = F1B2
                        * (m_solver->a_slope(cellId, sysEqn().PV->VV[i], j)
                           - m_solver->a_slope(cellId, sysEqn().PV->VV[j], i))
                        / m_solver->m_UInfinity;
              S[i][j] = F1B2
                        * (m_solver->a_slope(cellId, sysEqn().PV->VV[i], j)
                           + m_solver->a_slope(cellId, sysEqn().PV->VV[j], i))
                        / m_solver->m_UInfinity;
            }
          }
          for(MInt i = 0; i < nDim; i++) {
            for(MInt j = 0; j < nDim; j++) {
              vGrad[i][j] = F0;
              for(MInt k = 0; k < nDim; k++) {
                vGrad[i][j] += O[i][k] * O[k][j];
                vGrad[i][j] += S[i][k] * S[k][j];
              }
            }
          }
          maia::math::calcEigenValues(vGrad, eigenVal);
          sort(eigenVal, eigenVal + 3);
          lambda2(c) = eigenVal[1];
        }
        if(m_solver->m_vtuVorticityOutput) {
          ASSERT(nDim == 3, "");
          vorticity(c, 0) =
              F1B2
              * (m_solver->a_slope(cellId, sysEqn().PV->VV[2], 1) - m_solver->a_slope(cellId, sysEqn().PV->VV[1], 2));
          vorticity(c, 1) =
              F1B2
              * (m_solver->a_slope(cellId, sysEqn().PV->VV[0], 2) - m_solver->a_slope(cellId, sysEqn().PV->VV[2], 0));
          vorticity(c, 2) =
              F1B2
              * (m_solver->a_slope(cellId, sysEqn().PV->VV[1], 0) - m_solver->a_slope(cellId, sysEqn().PV->VV[0], 1));
        }
        if(m_solver->m_vtuVelocityGradientOutput) {
          for(MInt i = 0; i < nDim; i++) {
            for(MInt j = 0; j < nDim; j++) {
              velGrad(c, i * nDim + j) = m_solver->a_slope(cellId, sysEqn().PV->VV[i], j);
            }
          }
        }
      }
    }
    IF_CONSTEXPR(nDim == 2)
    for(MInt c = 0; c < (signed)noCells; c++)
      velocity(c, 2) = (float)0.0;

    insertDataHeader(reinterpret_cast<char*>(&velocity(0)), velocityMemsize, velocityMemsizeGlobal, velocityOffset);
    if(m_solver->m_vtuVorticityOutput)
      insertDataHeader(reinterpret_cast<char*>(&vorticity(0)), vorticityMemsize, vorticityMemsizeGlobal,
                       vorticityOffset);
    if(m_solver->m_vtuVelocityGradientOutput)
      insertDataHeader(reinterpret_cast<char*>(&velGrad(0)), velGradMemsize, velGradMemsizeGlobal, velGradOffset);
    insertDataHeader(reinterpret_cast<char*>(&pressure(0)), pressureMemsize, pressureMemsizeGlobal, pressureOffset);
    if(m_solver->m_vtuDensityOutput)
      insertDataHeader(reinterpret_cast<char*>(&density(0)), densityMemsize, densityMemsizeGlobal, densityOffset);
    if(m_solver->m_noSpecies > 0)
      insertDataHeader(reinterpret_cast<char*>(&progress(0)), progressMemsize, progressMemsizeGlobal, progressOffset);
    if(m_solver->m_vtuLevelSetOutput)
      insertDataHeader(reinterpret_cast<char*>(&levelSet(0)), levelSetMemsize, levelSetMemsizeGlobal, levelSetOffset);
    if(m_solver->m_vtuQCriterionOutput)
      insertDataHeader(reinterpret_cast<char*>(&Q(0)), QMemsize, QMemsizeGlobal, QOffset);
    if(m_solver->m_vtuLambda2Output)
      insertDataHeader(reinterpret_cast<char*>(&lambda2(0)), L2Memsize, L2MemsizeGlobal, L2Offset);
    if(noUserVars > 0)
      insertDataHeader(reinterpret_cast<char*>(&userVarsOut(0)), userVarsMemsize, userVarsMemsizeGlobal,
                       userVarsOffset);
    if(noSolverSpecificVars > 0)
      insertDataHeader(reinterpret_cast<char*>(&solverSpecificVarsOut(0)), solverSpecificVarsMemsize,
                       solverSpecificVarsMemsizeGlobal, solverSpecificVarsOffset);


    RECORD_TIMER_STOP(t_prepare);
    RECORD_TIMER_START(t_cells);

    if(domainId() == 0) {
      ofstream ofl;
      ofl.open(fileName.c_str(), ios_base::out | ios_base::trunc);
      if(ofl.is_open() && ofl.good()) {
        //================== VTKFile =================
        ofl << "<?xml version=\"1.0\"?>" << endl;
        ofl << R"(<VTKFile type="UnstructuredGrid" version="1.0" byte_order="LittleEndian" header_type="UInt64">)"
            << endl;
        ofl << "<UnstructuredGrid>" << endl;

        //================== Dimensions =================
        ofl << "<Piece NumberOfPoints=\"" << globalNoPoints << "\" NumberOfCells=\"" << globalNoCells << "\">" << endl;
        //================== /Dimensions =================

        //================== Points =================
        ofl << "<Points>" << endl;
        ofl << "<DataArray type=\"" << dataType << R"(" NumberOfComponents="3" format="appended" offset=")" << offset
            << "\"/>" << endl;
        offset += pointsMemsizeGlobal;
        ofl << "</Points>" << endl;
        //================== /Points =================

        //================== Cells =================
        ofl << "<Cells>" << endl;
        ofl << "<DataArray type=\"" << uIDataType << R"(" Name="connectivity" format="appended" offset=")" << offset
            << "\"/>" << endl;
        offset += connectivityMemsizeGlobal;
        ofl << "<DataArray type=\"" << uI64DataType << R"(" Name="offsets" format="appended" offset=")" << offset
            << "\"/>" << endl;
        offset += offsetsMemsizeGlobal;
        ofl << "<DataArray type=\"" << uI8DataType << R"(" Name="types" format="appended" offset=")" << offset << "\"/>"
            << endl;
        offset += typesMemsizeGlobal;
        if(m_solver->m_vtuCutCellOutput) {
          ofl << "<DataArray type=\"" << uIDataType << R"(" Name="faces" format="appended" offset=")" << offset
              << "\"/>" << endl;
          offset += facesMemsizeGlobal;
          ofl << "<DataArray type=\"" << uI64DataType << R"(" Name="faceoffsets" format="appended" offset=")" << offset
              << "\"/>" << endl;
          offset += faceoffsetsMemsizeGlobal;
        }
        ofl << "</Cells>" << endl;
        //================== /Cells =================

        //========== CellData - PointData ===========
        ofl << "<CellData Scalars=\"scalars\">" << endl;
        if(m_solver->m_vtuGlobalIdOutput) {
          ofl << "<DataArray type=\"" << uIDataType << R"(" Name="globalId" format="appended" offset=")" << offset
              << "\"/>" << endl;
          offset += globalIdMemsizeGlobal;
        }
        if(m_solver->m_vtuDomainIdOutput) {
          ofl << "<DataArray type=\"" << uIDataType << R"(" Name="domainId" format="appended" offset=")" << offset
              << "\"/>" << endl;
          offset += domainIdsMemsizeGlobal;
        }
        if(m_solver->m_vtuWritePointData) ofl << "</CellData>" << endl;
        if(m_solver->m_vtuWritePointData) ofl << "<PointData Scalars=\"scalars\">" << endl;
        ofl << "<DataArray type=\"" << dataType
            << R"(" Name="velocity" NumberOfComponents="3" format="appended" offset=")" << offset << "\"/>" << endl;
        offset += velocityMemsizeGlobal;
        if(m_solver->m_vtuVorticityOutput)
          ofl << "<DataArray type=\"" << dataType
              << R"(" Name="vorticity" NumberOfComponents="3" format="appended" offset=")" << offset << "\"/>" << endl;
        if(m_solver->m_vtuVorticityOutput) offset += vorticityMemsizeGlobal;
        if(m_solver->m_vtuVelocityGradientOutput)
          ofl << "<DataArray type=\"" << dataType
              << R"(" Name="velocityGradient" NumberOfComponents="9" format="appended" offset=")" << offset << "\"/>"
              << endl;
        if(m_solver->m_vtuVelocityGradientOutput) offset += velGradMemsizeGlobal;
        ofl << "<DataArray type=\"" << dataType << R"(" Name="pressure" format="appended" offset=")" << offset << "\"/>"
            << endl;
        offset += pressureMemsizeGlobal;
        if(m_solver->m_vtuDensityOutput) {
          ofl << "<DataArray type=\"" << dataType << R"(" Name="density" format="appended" offset=")" << offset
              << "\"/>" << endl;
          offset += densityMemsizeGlobal;
        }
        if(m_solver->m_noSpecies > 0) {
          ofl << "<DataArray type=\"" << dataType << R"(" Name="C" format="appended" offset=")" << offset << "\"/>"
              << endl;
          offset += progressMemsizeGlobal;
        }
        if(m_solver->m_vtuLevelSetOutput) {
          ofl << "<DataArray type=\"" << dataType << R"(" Name="levelSet" format="appended" offset=")" << offset
              << "\"/>" << endl;
          offset += levelSetMemsizeGlobal;
        }
        if(m_solver->m_vtuQCriterionOutput) {
          IF_CONSTEXPR(nDim == 3)
          ofl << "<DataArray type=\"" << dataType << R"(" Name="Q" format="appended" offset=")" << offset << "\"/>"
              << endl;
          IF_CONSTEXPR(nDim == 2)
          ofl << "<DataArray type=\"" << dataType << R"(" Name="vorticity" format="appended" offset=")" << offset
              << "\"/>" << endl;
          offset += QMemsizeGlobal;
        }
        if(m_solver->m_vtuLambda2Output) {
          ofl << "<DataArray type=\"" << dataType << R"(" Name="lambda2" format="appended" offset=")" << offset
              << "\"/>" << endl;
          offset += L2MemsizeGlobal;
        }
        if(m_solver->m_vtuWritePointData)
          ofl << "</PointData>" << endl;
        else
          ofl << "</CellData>" << endl;

        if(m_solver->m_levelSetMb && m_solver->m_constructGField ? noSolverSpecificVars > 1
                                                                 : noSolverSpecificVars > 0) {
          ofl << "<DataArray type=\"" << dataType << "\" Name=\"userVars\" NumberOfComponents=\""
              << noSolverSpecificVars << "\" format=\"appended\" offset=\"" << offset << "\"/>" << endl;
          offset += solverSpecificVarsMemsizeGlobal;
        }

        if(noUserVars > 0) {
          ofl << "<DataArray type=\"" << dataType << "\" Name=\"userVars\" NumberOfComponents=\"" << noUserVars
              << "\" format=\"appended\" offset=\"" << offset << "\"/>" << endl;
          offset += userVarsMemsizeGlobal;
        }
        //============ /CellData - /PointData ========

        ofl << "</Piece>" << endl;

        //================== FieldData =================
        ofl << "<FieldData>" << endl;

        ofl << "<DataArray type=\"" << uIDataType << R"(" Name="globalTimeStep" format="ascii" NumberOfTuples="1" > )"
            << globalTimeStep << " </DataArray>" << endl;
        if(m_solver->m_outputPhysicalTime) {
          ofl << "<DataArray type=\"" << dataType << R"(" Name="time" format="ascii" NumberOfTuples="1" > )"
              << m_solver->m_physicalTime << " </DataArray>" << endl;
          ofl << "<DataArray type=\"" << dataType << R"(" Name="internalTime" format="ascii" NumberOfTuples="1" > )"
              << m_solver->m_time << " </DataArray>" << endl;
        } else {
          ofl << "<DataArray type=\"" << dataType << R"(" Name="time" format="ascii" NumberOfTuples="1" > )"
              << m_solver->m_time << " </DataArray>" << endl;
          ofl << "<DataArray type=\"" << dataType << R"(" Name="physicalTime" format="ascii" NumberOfTuples="1" > )"
              << m_solver->m_physicalTime << " </DataArray>" << endl;
        }
        ofl << "<DataArray type=\"" << dataType << R"(" Name="timeStep" format="ascii" NumberOfTuples="1" > )"
            << m_solver->timeStep() << " </DataArray>" << endl;
        ofl << "<DataArray type=\"" << dataType << R"(" Name="timeRef" format="ascii" NumberOfTuples="1" > )"
            << m_solver->m_timeRef << " </DataArray>" << endl;
        ofl << "<DataArray type=\"" << dataType << R"(" Name="Ma" format="ascii" NumberOfTuples="1" > )"
            << m_solver->m_Ma << " </DataArray>" << endl;
        ofl << "<DataArray type=\"" << dataType << R"(" Name="Re" format="ascii" NumberOfTuples="1" > )"
            << m_solver->m_Re << " </DataArray>" << endl;
        ofl << "<DataArray type=\"" << dataType << R"(" Name="Re0" format="ascii" NumberOfTuples="1" > )"
            << sysEqn().m_Re0 << " </DataArray>" << endl;
        ofl << "<DataArray type=\"" << dataType << R"(" Name="Pr" format="ascii" NumberOfTuples="1" > )"
            << m_solver->m_Pr << " </DataArray>" << endl;
        ofl << "<DataArray type=\"" << dataType << R"(" Name="gamma" format="ascii" NumberOfTuples="1" > )"
            << m_solver->sysEqn().gamma_Ref() << " </DataArray>" << endl;
        ofl << "<DataArray type=\"" << dataType << R"(" Name="Tref" format="ascii" NumberOfTuples="1" > )"
            << m_solver->m_referenceTemperature << " </DataArray>" << endl;
        ofl << "<DataArray type=\"" << dataType << R"(" Name="Lref" format="ascii" NumberOfTuples="1" > )"
            << m_solver->m_referenceLength << " </DataArray>" << endl;
        ofl << "<DataArray type=\"" << dataType << R"(" Name="LrefPhys" format="ascii" NumberOfTuples="1" > )" << F1
            << " </DataArray>" << endl;
        ofl << "<DataArray type=\"" << dataType << R"(" Name="CFL" format="ascii" NumberOfTuples="1" > )"
            << m_solver->m_cfl << " </DataArray>" << endl;
        ofl << "<DataArray type=\"" << dataType << R"(" Name="uInf" format="ascii" NumberOfTuples="1" > )"
            << m_solver->m_UInfinity << " </DataArray>" << endl;
        ofl << "<DataArray type=\"" << dataType << R"(" Name="vInf" format="ascii" NumberOfTuples="1" > )"
            << m_solver->m_VInfinity << " </DataArray>" << endl;
        IF_CONSTEXPR(nDim == 3)
        ofl << "<DataArray type=\"" << dataType << R"(" Name="wInf" format="ascii" NumberOfTuples="1" > )"
            << m_solver->m_WInfinity << " </DataArray>" << endl;
        ofl << "<DataArray type=\"" << dataType << R"(" Name="pInf" format="ascii" NumberOfTuples="1" > )"
            << m_solver->m_PInfinity << " </DataArray>" << endl;
        ofl << "<DataArray type=\"" << dataType << R"(" Name="rhoInf" format="ascii" NumberOfTuples="1" > )"
            << m_solver->m_rhoInfinity << " </DataArray>" << endl;
        ofl << "<DataArray type=\"" << dataType << R"(" Name="TInf" format="ascii" NumberOfTuples="1" > )"
            << m_solver->m_TInfinity << " </DataArray>" << endl;
        ofl << "<DataArray type=\"" << dataType << R"(" Name="muInf" format="ascii" NumberOfTuples="1" > )"
            << sysEqn().m_muInfinity << " </DataArray>" << endl;
        ofl << "<DataArray type=\"" << dataType << R"(" Name="length0" format="ascii" NumberOfTuples="1" > )"
            << m_solver->c_cellLengthAtLevel(0) << " </DataArray>" << endl;
        ofl << "<DataArray type=\"" << dataType << R"(" Name="lengthMax" format="ascii" NumberOfTuples="1" > )"
            << m_solver->c_cellLengthAtLevel(m_solver->minLevel()) << " </DataArray>" << endl;
        ofl << "<DataArray type=\"" << dataType << R"(" Name="lengthMin" format="ascii" NumberOfTuples="1" > )"
            << m_solver->c_cellLengthAtLevel(m_solver->maxLevel()) << " </DataArray>" << endl;
        ofl << "<DataArray type=\"" << uIDataType << R"(" Name="minLevel" format="ascii" NumberOfTuples="1" > )"
            << m_solver->minLevel() << " </DataArray>" << endl;
        ofl << "<DataArray type=\"" << uIDataType << R"(" Name="maxLevel" format="ascii" NumberOfTuples="1" > )"
            << m_solver->maxLevel() << " </DataArray>" << endl;
        if(noDomains() > 1) {
          ofl << "<DataArray type=\"" << uIDataType << R"(" Name="domainOffsets" format="ascii" NumberOfTuples=")"
              << noDomains() << "\" > ";
          for(MInt i = 0; i < noDomains(); i++) {
            ofl << m_solver->domainOffset(i) << " ";
          }
          ofl << " </DataArray>" << endl;
        }
        ofl << "</FieldData>" << endl;
        //================== /FieldData =================

        ofl << "</UnstructuredGrid>" << endl;

        //================== AppendedData =================
        ofl << "<AppendedData encoding=\"raw\">" << endl;
        ofl << "_";
        ofl.close();
        ofl.clear();

        if(m_solver->m_vtuSaveHeaderTesting) {
          cout << system(("cp " + fileName + " " + fileName + "_header_testing").c_str());
        }
      } else {
        cerr << "ERROR! COULD NOT OPEN FILE " << fileName << " for writing! (1)" << endl;
      }
    }

    //###################################################
    //################### parallel IO ###################
    MPI_File file = nullptr;
    MInt rcode = MPI_File_open(mpiComm(), const_cast<char*>(fileName.c_str()), MPI_MODE_APPEND | MPI_MODE_WRONLY,
                               MPI_INFO_NULL, &file, AT_);
    if(rcode != MPI_SUCCESS) mTerm(1, AT_, "Error opening file " + fileName);

    //=================== Points =======================
    rcode = writeVtuArrayParallel(file, &points(0), pointsOffset, pointsMemsize, pointsMemsizeGlobal);
    if(rcode != MPI_SUCCESS) mTerm(1, AT_, "Error (1) writing to file " + fileName);

    //=================== Connectivity =======================
    rcode = writeVtuArrayParallel(file, &connectivity(0), connectivityOffset, connectivityMemsize,
                                  connectivityMemsizeGlobal);
    if(rcode != MPI_SUCCESS) {
      mTerm(1, AT_, "Error (2) writing to file " + fileName);
    }

    //=================== Offsets =======================
    rcode = writeVtuArrayParallel(file, &offsets(0), offsetsOffset, offsetsMemsize, offsetsMemsizeGlobal);
    if(rcode != MPI_SUCCESS) {
      mTerm(1, AT_, "Error (3) writing to file " + fileName);
    }

    //=================== Types =======================
    rcode = writeVtuArrayParallel(file, &types(0), typesOffset, typesMemsize, typesMemsizeGlobal);
    if(rcode != MPI_SUCCESS) mTerm(1, AT_, "Error (4) writing to file " + fileName);

    //=================== Faces =======================
    if(m_solver->m_vtuCutCellOutput)
      rcode = writeVtuArrayParallel(file, &faces(0), facesOffset, facesMemsize, facesMemsizeGlobal);
    if(rcode != MPI_SUCCESS) mTerm(1, AT_, "Error (5) writing to file " + fileName);

    //=================== FaceOffsets =======================
    if(m_solver->m_vtuCutCellOutput)
      rcode =
          writeVtuArrayParallel(file, &faceoffsets(0), faceoffsetsOffset, faceoffsetsMemsize, faceoffsetsMemsizeGlobal);
    if(rcode != MPI_SUCCESS) mTerm(1, AT_, "Error (5) writing to file " + fileName);

    RECORD_TIMER_STOP(t_cells);
    RECORD_TIMER_START(t_variables);

    //=================== GlobalId =======================
    if(m_solver->m_vtuGlobalIdOutput)
      rcode = writeVtuArrayParallel(file, &globalId(0), globalIdOffset, globalIdMemsize, globalIdMemsizeGlobal);
    if(rcode != MPI_SUCCESS) mTerm(1, AT_, "Error (5) writing to file " + fileName);
    //=================== GlobalId =======================

    if(m_solver->m_vtuDomainIdOutput)
      rcode = writeVtuArrayParallel(file, &domainIds(0), domainIdsOffset, domainIdsMemsize, domainIdsMemsizeGlobal);
    if(rcode != MPI_SUCCESS) mTerm(1, AT_, "Error (5) writing to file " + fileName);

    //=================== velocity =======================
    rcode = writeVtuArrayParallel(file, &velocity(0), velocityOffset, velocityMemsize, velocityMemsizeGlobal);
    if(rcode != MPI_SUCCESS) mTerm(1, AT_, "Error (6) writing to file " + fileName);

    //=================== vorticity =======================
    if(m_solver->m_vtuVorticityOutput)
      rcode = writeVtuArrayParallel(file, &vorticity(0), vorticityOffset, vorticityMemsize, vorticityMemsizeGlobal);
    if(rcode != MPI_SUCCESS) mTerm(1, AT_, "Error (7) writing to file " + fileName);

    //=================== velGrad =======================
    if(m_solver->m_vtuVelocityGradientOutput)
      rcode = writeVtuArrayParallel(file, &velGrad(0), velGradOffset, velGradMemsize, velGradMemsizeGlobal);
    if(rcode != MPI_SUCCESS) mTerm(1, AT_, "Error (7) writing to file " + fileName);

    //=================== pressure =======================
    rcode = writeVtuArrayParallel(file, &pressure(0), pressureOffset, pressureMemsize, pressureMemsizeGlobal);
    if(rcode != MPI_SUCCESS) {
      mTerm(1, AT_, "Error (8) writing to file " + fileName);
    }

    //=================== density =======================
    if(m_solver->m_vtuDensityOutput)
      rcode = writeVtuArrayParallel(file, &density(0), densityOffset, densityMemsize, densityMemsizeGlobal);
    if(rcode != MPI_SUCCESS) mTerm(1, AT_, "Error (9) writing to file " + fileName);

    if(m_solver->m_noSpecies > 0) {
      //=================== progress =======================
      rcode = writeVtuArrayParallel(file, &progress(0), progressOffset, progressMemsize, progressMemsizeGlobal);
      if(rcode != MPI_SUCCESS) {
        mTerm(1, AT_, "Error (9) writing to file " + fileName);
      }
    }
    //=================== levelSet =======================
    if(m_solver->m_vtuLevelSetOutput) {
      rcode = writeVtuArrayParallel(file, &levelSet(0), levelSetOffset, levelSetMemsize, levelSetMemsizeGlobal);
      if(rcode != MPI_SUCCESS) mTerm(1, AT_, "Error (9) writing to file " + fileName);
    }
    //=================== Q =======================
    if(m_solver->m_vtuQCriterionOutput) rcode = writeVtuArrayParallel(file, &Q(0), QOffset, QMemsize, QMemsizeGlobal);
    if(rcode != MPI_SUCCESS) mTerm(1, AT_, "Error (10) writing to file " + fileName);

    //=================== lambda2 =======================
    if(m_solver->m_vtuLambda2Output)
      rcode = writeVtuArrayParallel(file, &lambda2(0), L2Offset, L2Memsize, L2MemsizeGlobal);
    if(rcode != MPI_SUCCESS) {
      mTerm(1, AT_, "Error (11) writing to file " + fileName);
    }

    //=================== solverSpecificVars =======================
    if(m_solver->m_levelSetMb && m_solver->m_constructGField ? noSolverSpecificVars > 1 : noSolverSpecificVars > 0)
      rcode = writeVtuArrayParallel(file, &solverSpecificVarsOut(0), solverSpecificVarsOffset,
                                    solverSpecificVarsMemsize, solverSpecificVarsMemsizeGlobal);
    if(rcode != MPI_SUCCESS) {
      mTerm(1, AT_, "Error (12) writing to file " + fileName);
    }

    //=================== userVars =======================
    if(noUserVars > 0)
      rcode = writeVtuArrayParallel(file, &userVarsOut(0), userVarsOffset, userVarsMemsize, userVarsMemsizeGlobal);
    if(rcode != MPI_SUCCESS) {
      mTerm(1, AT_, "Error (13) writing to file " + fileName);
    }

    MPI_File_close(&file, AT_);
    //###################################################
    //###################################################

    if(domainId() == 0) {
      ofstream ofl;
      ofl.open(fileName.c_str(), ios_base::out | ios_base::app);
      if(ofl.is_open() && ofl.good()) {
        ofl << endl;
        ofl << "</AppendedData>" << endl;
        ofl << "</VTKFile>" << endl;
        ofl.close();
        ofl.clear();
        //================== /AppendedData =================
        //================== /VTKFile ======================
      } else {
        cerr << "ERROR! COULD NOT OPEN FILE " << fileName << " for writing! (2)" << endl;
      }
    }

    //-----------------------
    //------ QOUT.pvd -------
    //-----------------------

    if(domainId() == 0) {
      if(m_solver->m_firstUseWriteVtuOutputParallelQout) {
        m_solver->m_firstUseWriteVtuOutputParallelQout = !m_solver->m_restart;
      }
      ofstream ofile;
      ofstream ofile2;
      if(m_solver->m_firstUseWriteVtuOutputParallelQout)
        ofile.open("out/QOUT.pvd.tmp", ios_base::out | ios_base::trunc);
      else
        ofile.open("out/QOUT.pvd.tmp", ios_base::out | ios_base::app);
      if(ofile.is_open() && ofile.good()) {
        if(m_solver->m_firstUseWriteVtuOutputParallelQout) {
          ofile << R"(<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">)" << endl;
          ofile << "<Collection>" << endl;
        }
        const size_t pos = fileName.find(m_solver->m_solutionOutput);
        MString fname =
            (pos == string::npos) ? fileName : fileName.substr(pos + m_solver->m_solutionOutput.size(), string::npos);
        ofile << "<DataSet part=\"" << 0 << "\" timestep=\"" << globalTimeStep << "\" file=\""
              << "./" << fname << "\"/>" << endl;
        ofile.close();
        ofile.clear();
      } else {
        cerr << "Error opening file out/QOUT.pvd.tmp" << endl;
      }
      m_log << "Command executed with return value: " << system("cp out/QOUT.pvd.tmp out/QOUT.pvd") << " at " << AT_
            << endl;

      ofile2.open("out/QOUT.pvd", ios_base::out | ios_base::app);
      if(ofile2.is_open() && ofile2.good()) {
        ofile2 << "</Collection>" << endl;
        ofile2 << "</VTKFile>" << endl;
        ofile2.close();
        ofile2.clear();
      } else {
        cerr << "Error opening file out/QOUT.pvd" << endl;
      }
      m_solver->m_firstUseWriteVtuOutputParallelQout = false;
    }

    RECORD_TIMER_STOP(t_variables);
    RECORD_TIMER_START(t_geometry);
  }


  if(m_solver->m_vtuCutCellOutput
     && (m_solver->m_vtuGeometryOutput.size() > 0 || string2enum(m_solver->m_outputFormat) == VTP)) {
    if(m_solver->m_vtuLevelThreshold < m_solver->maxRefinementLevel()) {
      if(m_solver->m_extractedCells) {
        delete m_solver->m_extractedCells;
        m_solver->m_extractedCells = nullptr;
      }
      if(m_solver->m_gridPoints) {
        delete m_solver->m_gridPoints;
        m_solver->m_gridPoints = nullptr;
      }
      m_solver->extractPointIdsFromGrid(m_solver->m_extractedCells, m_solver->m_gridPoints, false,
                                        m_solver->m_splitChildToSplitCell, m_solver->maxRefinementLevel(),
                                        m_solver->m_vtuCoordinatesThreshold);
      facesSize = 0;
      noPolyhedra = 0;
      for(MInt c = 0; c < m_solver->m_extractedCells->size(); c++) {
        MInt cellId = m_solver->m_extractedCells->a[c].m_cellId;
        MInt bndryId = m_solver->a_bndryId(cellId);
        if(bndryId > -1) {
          polyhedralCell[c] = noPolyhedra;
          noPolyhedra++;
        }
      }
      connSize = vtuAssembleFaceStream(facestream, conn, facesSize, polyhedralCell, noPolyhedra, false);
    }

    if(!(facestream == nullptr || conn == nullptr)) {
      uint64_t noPolygons = 0;
      uint64_t geomConnSize = 0;
      for(MInt c = 0; c < m_solver->m_extractedCells->size(); c++) {
        const MInt bndryId = m_solver->a_bndryId(m_solver->m_extractedCells->a[c].m_cellId);
        const MInt pc = polyhedralCell[c];
        if(bndryId < 0) continue;
        uint_t cnt = 1;
        for(MInt srfc = 0; srfc < m_solver->m_bndryCells->a[bndryId].m_noSrfcs; srfc++) {
          if(m_solver->m_vtuGeometryOutput.count(m_solver->m_bndryCells->a[bndryId].m_srfcs[srfc]->m_bndryCndId) == 0)
            continue;
          noPolygons++;
          geomConnSize += (uint64_t)facestream[pc][cnt];
          cnt += (uint_t)(facestream[pc][cnt] + 1);
        }
      }

      ScratchSpace<float> geomPoints(geomConnSize + 4, 3, AT_, "geomPoints");
      ScratchSpace<uint_t> geomConn(geomConnSize + 4, AT_, "geomConn");
      ScratchSpace<uint_t> geomOffsets(noPolygons + 4, AT_, "geomOffsets");
      ScratchSpace<uint_t> cutPoints(m_solver->m_gridPoints->size(), AT_, "cutPoints");
      cutPoints.fill(std::numeric_limits<uint_t>::max());
      uint64_t noGeomPoints = 0;
      geomConnSize = 0;
      noPolygons = 0;
      for(MInt c = 0; c < m_solver->m_extractedCells->size(); c++) {
        const MInt bndryId = m_solver->a_bndryId(m_solver->m_extractedCells->a[c].m_cellId);
        const MInt pc = polyhedralCell[c];
        if(pc < 0) continue;
        if(bndryId < 0) continue;
        if(facestream[pc].empty()) continue;
        ASSERT((signed)facestream[pc][0] >= m_solver->m_bndryCells->a[bndryId].m_noSrfcs, "");
        uint_t cnt = 1;
        for(MInt srfc = 0; srfc < m_solver->m_bndryCells->a[bndryId].m_noSrfcs; srfc++) {
          if(m_solver->m_vtuGeometryOutput.count(m_solver->m_bndryCells->a[bndryId].m_srfcs[srfc]->m_bndryCndId) == 0)
            continue;
          ASSERT(cnt < facestream[pc].size(), "");
          MInt noCutPoints = (signed)facestream[pc][cnt];
          cnt++;
          for(MInt cp = 0; cp < noCutPoints; cp++) {
            ASSERT(cnt < facestream[pc].size(), "");
            uint_t gridPointId = (uint_t)facestream[pc][cnt];
            uint_t pId = cutPoints[gridPointId];
            if(pId >= (unsigned)m_solver->m_gridPoints->size()) {
              pId = noGeomPoints;
              cutPoints[gridPointId] = pId;
              for(MInt i = 0; i < nDim; i++) {
                geomPoints((signed)pId, i) = m_solver->m_gridPoints->a[gridPointId].m_coordinates[i];
              }
              noGeomPoints++;
            }
            geomConn[geomConnSize] = pId;
            cnt++;
            geomConnSize++;
          }
          geomOffsets[noPolygons] = geomConnSize;
          noPolygons++;
        }
      }

      ScratchSpace<uint64_t> noGeomPointsPerDomain(noDomains(), AT_, "noGeomPointsPerDomain");
      ScratchSpace<uint64_t> noPolygonsPerDomain(noDomains(), AT_, "noPolygonsPerDomain");
      ScratchSpace<uint64_t> geomConnSizePerDomain(noDomains(), AT_, "geomConnSizePerDomain");
      ScratchSpace<uint64_t> geomDataPerDomain(noDomains(), 3, AT_, "geomDataPerDomain");
      geomDataPerDomain(domainId(), 0) = noGeomPoints;
      geomDataPerDomain(domainId(), 1) = noPolygons;
      geomDataPerDomain(domainId(), 2) = geomConnSize;
      MPI_Allgather(MPI_IN_PLACE, 3, MPI_UINT64_T, &geomDataPerDomain[0], 3, MPI_UINT64_T, mpiComm(), AT_,
                    "MPI_IN_PLACE", "geomDataPerDomain[0]");
      for(MInt d = 0; d < noDomains(); d++) {
        noGeomPointsPerDomain[d] = geomDataPerDomain(d, 0);
        noPolygonsPerDomain[d] = geomDataPerDomain(d, 1);
        geomConnSizePerDomain[d] = geomDataPerDomain(d, 2);
      }
      uint64_t globalNoGeomPoints = 0;
      uint64_t globalGeomPointOffset = 0;
      uint64_t globalNoPolygons = 0;
      uint64_t globalPolygonsOffset = 0;
      uint64_t globalGeomConnSize = 0;
      uint64_t globalGeomConnSizeOffset = 0;
      for(MInt d = 0; d < noDomains(); d++) {
        globalNoGeomPoints += noGeomPointsPerDomain[d];
        globalNoPolygons += noPolygonsPerDomain[d];
        globalGeomConnSize += geomConnSizePerDomain[d];
      }
      for(MInt d = 0; d < domainId(); d++) {
        globalGeomPointOffset += noGeomPointsPerDomain[d];
        globalPolygonsOffset += noPolygonsPerDomain[d];
        globalGeomConnSizeOffset += geomConnSizePerDomain[d];
      }

      for(MUint p = 0; p < noPolygons; p++) {
        geomOffsets[p] += globalGeomConnSizeOffset;
      }
      for(MUint c = 0; c < geomConnSize; c++) {
        geomConn[c] += globalGeomPointOffset;
      }

      if(globalNoPolygons == 0) {
        RECORD_TIMER_STOP(t_geometry);
        RECORD_TIMER_STOP(t_vtutotal);
        DISPLAY_TIMER_INTERM(t_vtutotal);
        return;
      }


      //-----------
      // points
      uint64_t geomPointsMemsize = 3 * noGeomPoints * sizeof(float);
      uint64_t geomPointsMemsizeGlobal = 3 * globalNoGeomPoints * sizeof(float);
      uint64_t geomPointsOffset = 3 * globalGeomPointOffset * sizeof(float);
      insertDataHeader(reinterpret_cast<char*>(&geomPoints(0)), geomPointsMemsize, geomPointsMemsizeGlobal,
                       geomPointsOffset);

      //-----------
      // connectivity
      uint64_t geomConnMemsize = geomConnSize * sizeof(uint_t);
      uint64_t geomConnOffset = globalGeomConnSizeOffset * sizeof(uint_t);
      uint64_t geomConnMemsizeGlobal = globalGeomConnSize * sizeof(uint_t);
      insertDataHeader(reinterpret_cast<char*>(&geomConn(0)), geomConnMemsize, geomConnMemsizeGlobal, geomConnOffset);

      //-----------
      // offsets
      uint64_t geomOffsetsMemsize = noPolygons * sizeof(uint_t);
      uint64_t geomOffsetsOffset = globalPolygonsOffset * sizeof(uint_t);
      uint64_t geomOffsetsMemsizeGlobal = globalNoPolygons * sizeof(uint_t);
      insertDataHeader(reinterpret_cast<char*>(&geomOffsets(0)), geomOffsetsMemsize, geomOffsetsMemsizeGlobal,
                       geomOffsetsOffset);

      //-----------
      // variables
      uint64_t bodyIdMemsize = noPolygons * sizeof(uint_t);
      uint64_t velocityMemsize = 3 * noPolygons * sizeof(float);
      uint64_t domainIdsMemsize = noPolygons * sizeof(uint_t);
      uint64_t pressureMemsize = noPolygons * sizeof(float);
      uint64_t densityMemsize = noPolygons * sizeof(float);
      uint64_t normalsMemsize = 3 * noPolygons * sizeof(float);
      uint64_t shearMemsize = 3 * noPolygons * sizeof(float);
      uint64_t bodyIdOffset = globalPolygonsOffset * sizeof(uint_t);
      uint64_t velocityOffset = 3 * globalPolygonsOffset * sizeof(float);
      uint64_t domainIdsOffset = globalPolygonsOffset * sizeof(uint_t);
      uint64_t pressureOffset = globalPolygonsOffset * sizeof(float);
      uint64_t densityOffset = globalPolygonsOffset * sizeof(float);
      uint64_t normalsOffset = 3 * globalPolygonsOffset * sizeof(float);
      uint64_t shearOffset = 3 * globalPolygonsOffset * sizeof(float);
      uint64_t bodyIdMemsizeGlobal = globalNoPolygons * sizeof(uint_t);
      uint64_t velocityMemsizeGlobal = 3 * globalNoPolygons * sizeof(float);
      uint64_t domainIdsMemsizeGlobal = globalNoPolygons * sizeof(uint_t);
      uint64_t pressureMemsizeGlobal = globalNoPolygons * sizeof(float);
      uint64_t densityMemsizeGlobal = globalNoPolygons * sizeof(float);
      uint64_t normalsMemsizeGlobal = 3 * globalNoPolygons * sizeof(float);
      uint64_t shearMemsizeGlobal = 3 * globalNoPolygons * sizeof(float);
      const MInt bodyIdMaxSize = estimateMemorySizeSolverwise(noPolygons, noPolygonsPerDomain, sizeof(uint_t));
      const MInt vectorMaxSize = estimateMemorySizeSolverwise(noPolygons, noPolygonsPerDomain, 3 * sizeof(float));
      const MInt scalarMaxSize = estimateMemorySizeSolverwise(noPolygons, noPolygonsPerDomain, sizeof(float));
      ScratchSpace<uint_t> bodyId(bodyIdMaxSize, AT_, "bodyId");
      ScratchSpace<float> velocity(vectorMaxSize, 3, AT_, "velocity");
      ScratchSpace<uint_t> domainIds(m_solver->m_vtuGeometryOutputExtended ? bodyIdMaxSize : 1, AT_, "domainIds");
      ScratchSpace<float> pressure(m_solver->m_vtuGeometryOutputExtended ? scalarMaxSize : 1, AT_, "pressure");
      ScratchSpace<float> density(m_solver->m_vtuGeometryOutputExtended ? scalarMaxSize : 1, AT_, "density");
      ScratchSpace<float> normals(m_solver->m_vtuGeometryOutputExtended ? vectorMaxSize : 1, 3, AT_, "normals");
      ScratchSpace<float> shear(m_solver->m_vtuGeometryOutputExtended ? vectorMaxSize : 1, 3, AT_, "shear");
      const MFloat F1BRe = F1 / (m_solver->m_referenceLength * sysEqn().m_Re0);
      MInt cnt = 0;
      for(MInt c = 0; c < m_solver->m_extractedCells->size(); c++) {
        MInt bndryId = m_solver->a_bndryId(m_solver->m_extractedCells->a[c].m_cellId);
        MInt pc = polyhedralCell[c];
        if(pc < 0) continue;
        if(bndryId < 0) continue;
        for(MInt srfc = 0; srfc < m_solver->m_bndryCells->a[bndryId].m_noSrfcs; srfc++) {
          if(m_solver->m_vtuGeometryOutput.count(m_solver->m_bndryCells->a[bndryId].m_srfcs[srfc]->m_bndryCndId) == 0)
            continue;
          bodyId(cnt) = (uint_t)m_solver->m_bndryCells->a[bndryId].m_srfcs[srfc]->m_bodyId[0];
          if(m_solver->m_internalBodyId != nullptr && (signed)bodyId(cnt) > -1
             && (signed)bodyId(cnt) >= m_solver->m_noEmbeddedBodies)
            bodyId(cnt) = (uint_t)(m_solver->m_internalBodyId[bodyId(cnt)]);
          if(m_solver->m_bndryCells->a[bndryId].m_srfcs[srfc]->m_bodyId[0] < 0)
            bodyId(cnt) = std::numeric_limits<uint_t>::max();
          for(MInt i = 0; i < nDim; i++) {
            velocity(cnt, i) =
                (float)m_solver->m_bndryCells->a[bndryId].m_srfcVariables[srfc]->m_primVars[sysEqn().PV->VV[i]];
          }

          if(m_solver->m_vtuGeometryOutputExtended) {
            domainIds(cnt) = (uint_t)domainId();
            pressure(cnt) = (float)m_solver->m_bndryCells->a[bndryId].m_srfcVariables[srfc]->m_primVars[sysEqn().PV->P];
            density(cnt) =
                (float)m_solver->m_bndryCells->a[bndryId].m_srfcVariables[srfc]->m_primVars[sysEqn().PV->RHO];
            MFloat mue = m_solver->sysEqn().sutherlandLaw(m_solver->sysEqn().temperature_ES(
                m_solver->m_bndryCells->a[bndryId].m_srfcVariables[srfc]->m_primVars[sysEqn().PV->RHO],
                m_solver->m_bndryCells->a[bndryId].m_srfcVariables[srfc]->m_primVars[sysEqn().PV->P]));
            for(MInt i = 0; i < nDim; i++) {
              normals(cnt, i) = (float)(m_solver->m_bndryCells->a[bndryId].m_srfcs[srfc]->m_normalVector[i]);
              shear(cnt, i) = (float)(F1BRe * mue
                                      * m_solver->m_bndryCells->a[bndryId]
                                            .m_srfcVariables[srfc]
                                            ->m_normalDeriv[sysEqn().PV->VV[i]]);
            }
          }
          cnt++;
        }
      }
      insertDataHeader(reinterpret_cast<char*>(&bodyId(0)), bodyIdMemsize, bodyIdMemsizeGlobal, bodyIdOffset);
      insertDataHeader(reinterpret_cast<char*>(&velocity(0)), velocityMemsize, velocityMemsizeGlobal, velocityOffset);
      if(m_solver->m_vtuGeometryOutputExtended) {
        insertDataHeader(reinterpret_cast<char*>(&domainIds(0)), domainIdsMemsize, domainIdsMemsizeGlobal,
                         domainIdsOffset);
        insertDataHeader(reinterpret_cast<char*>(&pressure(0)), pressureMemsize, pressureMemsizeGlobal, pressureOffset);
        insertDataHeader(reinterpret_cast<char*>(&density(0)), densityMemsize, densityMemsizeGlobal, densityOffset);
        insertDataHeader(reinterpret_cast<char*>(&normals(0)), normalsMemsize, normalsMemsizeGlobal, normalsOffset);
        insertDataHeader(reinterpret_cast<char*>(&shear(0)), shearMemsize, shearMemsizeGlobal, shearOffset);
      }

      if(m_solver->m_vtuWriteGeometryFile) {
        if(domainId() == 0) {
          offset = 0;
          ofstream ofl;
          ofl.open(geomFileName.c_str(), ios_base::out | ios_base::trunc);
          if(ofl.is_open() && ofl.good()) {
            //================== VTKFile =================
            ofl << "<?xml version=\"1.0\"?>" << endl;
            ofl << R"(<VTKFile type="PolyData" version="1.0" byte_order="LittleEndian" header_type="UInt64">)" << endl;
            ofl << "<PolyData>" << endl;

            //================== FieldData =================
            ofl << "<FieldData>" << endl;

            ofl << "<DataArray type=\"" << uIDataType
                << R"(" Name="globalTimeStep" format="ascii" NumberOfTuples="1" > )" << globalTimeStep
                << " </DataArray>" << endl;
            if(m_solver->m_outputPhysicalTime) {
              ofl << "<DataArray type=\"" << dataType << R"(" Name="time" format="ascii" NumberOfTuples="1" > )"
                  << m_solver->m_physicalTime << " </DataArray>" << endl;
              ofl << "<DataArray type=\"" << dataType << R"(" Name="internalTime" format="ascii" NumberOfTuples="1" > )"
                  << m_solver->m_time << " </DataArray>" << endl;
            } else {
              ofl << "<DataArray type=\"" << dataType << R"(" Name="time" format="ascii" NumberOfTuples="1" > )"
                  << m_solver->m_time << " </DataArray>" << endl;
              ofl << "<DataArray type=\"" << dataType << R"(" Name="physicalTime" format="ascii" NumberOfTuples="1" > )"
                  << m_solver->m_physicalTime << " </DataArray>" << endl;
            }
            ofl << "<DataArray type=\"" << dataType << R"(" Name="timeStep" format="ascii" NumberOfTuples="1" > )"
                << m_solver->timeStep() << " </DataArray>" << endl;
            ofl << "<DataArray type=\"" << dataType << R"(" Name="timeRef" format="ascii" NumberOfTuples="1" > )"
                << m_solver->m_timeRef << " </DataArray>" << endl;
            ofl << "<DataArray type=\"" << dataType << R"(" Name="Ma" format="ascii" NumberOfTuples="1" > )"
                << m_solver->m_Ma << " </DataArray>" << endl;
            ofl << "<DataArray type=\"" << dataType << R"(" Name="Re" format="ascii" NumberOfTuples="1" > )"
                << m_solver->m_Re << " </DataArray>" << endl;
            ofl << "<DataArray type=\"" << dataType << R"(" Name="Re0" format="ascii" NumberOfTuples="1" > )"
                << sysEqn().m_Re0 << " </DataArray>" << endl;
            ofl << "<DataArray type=\"" << dataType << R"(" Name="Pr" format="ascii" NumberOfTuples="1" > )"
                << m_solver->m_Pr << " </DataArray>" << endl;
            ofl << "<DataArray type=\"" << dataType << R"(" Name="gamma" format="ascii" NumberOfTuples="1" > )"
                << m_solver->sysEqn().gamma_Ref() << " </DataArray>" << endl;
            ofl << "<DataArray type=\"" << dataType << R"(" Name="Tref" format="ascii" NumberOfTuples="1" > )"
                << m_solver->m_referenceTemperature << " </DataArray>" << endl;
            ofl << "<DataArray type=\"" << dataType << R"(" Name="Lref" format="ascii" NumberOfTuples="1" > )"
                << m_solver->m_referenceLength << " </DataArray>" << endl;
            ofl << "<DataArray type=\"" << dataType << R"(" Name="LrefPhys" format="ascii" NumberOfTuples="1" > )" << F1
                << " </DataArray>" << endl;
            ofl << "<DataArray type=\"" << dataType << R"(" Name="CFL" format="ascii" NumberOfTuples="1" > )"
                << m_solver->m_cfl << " </DataArray>" << endl;
            ofl << "<DataArray type=\"" << dataType << R"(" Name="uInf" format="ascii" NumberOfTuples="1" > )"
                << m_solver->m_UInfinity << " </DataArray>" << endl;
            ofl << "<DataArray type=\"" << dataType << R"(" Name="vInf" format="ascii" NumberOfTuples="1" > )"
                << m_solver->m_VInfinity << " </DataArray>" << endl;
            IF_CONSTEXPR(nDim == 3)
            ofl << "<DataArray type=\"" << dataType << R"(" Name="wInf" format="ascii" NumberOfTuples="1" > )"
                << m_solver->m_WInfinity << " </DataArray>" << endl;
            ofl << "<DataArray type=\"" << dataType << R"(" Name="pInf" format="ascii" NumberOfTuples="1" > )"
                << m_solver->m_PInfinity << " </DataArray>" << endl;
            ofl << "<DataArray type=\"" << dataType << R"(" Name="rhoInf" format="ascii" NumberOfTuples="1" > )"
                << m_solver->m_rhoInfinity << " </DataArray>" << endl;
            ofl << "<DataArray type=\"" << dataType << R"(" Name="TInf" format="ascii" NumberOfTuples="1" > )"
                << m_solver->m_TInfinity << " </DataArray>" << endl;
            ofl << "<DataArray type=\"" << dataType << R"(" Name="muInf" format="ascii" NumberOfTuples="1" > )"
                << sysEqn().m_muInfinity << " </DataArray>" << endl;
            ofl << "<DataArray type=\"" << dataType << R"(" Name="length0" format="ascii" NumberOfTuples="1" > )"
                << m_solver->c_cellLengthAtLevel(0) << " </DataArray>" << endl;
            ofl << "<DataArray type=\"" << dataType << R"(" Name="lengthMax" format="ascii" NumberOfTuples="1" > )"
                << m_solver->c_cellLengthAtLevel(m_solver->minLevel()) << " </DataArray>" << endl;
            ofl << "<DataArray type=\"" << dataType << R"(" Name="lengthMin" format="ascii" NumberOfTuples="1" > )"
                << m_solver->c_cellLengthAtLevel(m_solver->maxLevel()) << " </DataArray>" << endl;
            ofl << "<DataArray type=\"" << uIDataType << R"(" Name="minLevel" format="ascii" NumberOfTuples="1" > )"
                << m_solver->minLevel() << " </DataArray>" << endl;
            ofl << "<DataArray type=\"" << uIDataType << R"(" Name="maxLevel" format="ascii" NumberOfTuples="1" > )"
                << m_solver->maxLevel() << " </DataArray>" << endl;
            if(noDomains() > 1) {
              ofl << "<DataArray type=\"" << uIDataType << R"(" Name="domainOffsets" format="ascii" NumberOfTuples=")"
                  << noDomains() << "\" > ";
              for(MInt i = 0; i < noDomains(); i++) {
                ofl << m_solver->domainOffset(i) << " ";
              }
              ofl << " </DataArray>" << endl;
            }

            ofl << "</FieldData>" << endl;
            //================== /FieldData =================


            //================== Dimensions =================
            ofl << "<Piece NumberOfPoints=\"" << globalNoGeomPoints << "\" NumberOfPolys=\"" << globalNoPolygons
                << "\">" << endl;
            //================== /Dimensions ================


            //================== Points =================
            ofl << "<Points>" << endl;
            ofl << "<DataArray type=\"" << dataType << R"(" NumberOfComponents="3" format="appended" offset=")"
                << offset << "\"/>" << endl;
            offset += geomPointsMemsizeGlobal;
            ofl << "</Points>" << endl;
            //================== /Points =================


            //================== Polys =================
            ofl << "<Polys>" << endl;

            ofl << "<DataArray type=\"" << uIDataType << R"(" Name="connectivity" format="appended" offset=")" << offset
                << "\"/>" << endl;
            offset += geomConnMemsizeGlobal;

            ofl << "<DataArray type=\"" << uIDataType << R"(" Name="offsets" format="appended" offset=")" << offset
                << "\"/>" << endl;
            offset += geomOffsetsMemsizeGlobal;

            ofl << "</Polys>" << endl;
            //================== /Polys =================


            //================== PolyData =================
            ofl << "<CellData Scalars=\"scalars\">" << endl;

            ofl << "<DataArray type=\"" << uIDataType << R"(" Name="bodyId" format="appended" offset=")" << offset
                << "\"/>" << endl;
            offset += bodyIdMemsizeGlobal;

            ofl << "<DataArray type=\"" << dataType
                << R"(" Name="velocity" NumberOfComponents="3" format="appended" offset=")" << offset << "\"/>" << endl;
            offset += velocityMemsizeGlobal;

            if(m_solver->m_vtuGeometryOutputExtended) {
              ofl << "<DataArray type=\"" << uIDataType << R"(" Name="domainId" format="appended" offset=")" << offset
                  << "\"/>" << endl;
              offset += domainIdsMemsizeGlobal;

              ofl << "<DataArray type=\"" << dataType << R"(" Name="pressure" format="appended" offset=")" << offset
                  << "\"/>" << endl;
              offset += pressureMemsizeGlobal;

              ofl << "<DataArray type=\"" << dataType << R"(" Name="density" format="appended" offset=")" << offset
                  << "\"/>" << endl;
              offset += densityMemsizeGlobal;

              ofl << "<DataArray type=\"" << dataType
                  << R"(" Name="cell_normals" NumberOfComponents="3" format="appended" offset=")" << offset << "\"/>"
                  << endl;
              offset += normalsMemsizeGlobal;

              ofl << "<DataArray type=\"" << dataType
                  << R"(" Name="shear" NumberOfComponents="3" format="appended" offset=")" << offset << "\"/>" << endl;
              offset += shearMemsizeGlobal;
            }

            ofl << "</CellData>" << endl;
            //================== /PolyData =================


            ofl << "</Piece>" << endl;
            ofl << "</PolyData>" << endl;


            //================== AppendedData =================
            ofl << "<AppendedData encoding=\"raw\">" << endl;
            ofl << "_";

            ofl.close();
            ofl.clear();

          } else {
            cerr << "ERROR! COULD NOT OPEN FILE " << geomFileName << " for writing! (1)" << endl;
          }
        }


        //###################################################
        //################### parallel IO ###################
        MPI_File gfile = nullptr;
        MInt rcode = MPI_File_open(mpiComm(), const_cast<char*>(geomFileName.c_str()),
                                   MPI_MODE_APPEND | MPI_MODE_WRONLY, MPI_INFO_NULL, &gfile, AT_);
        if(rcode != MPI_SUCCESS) mTerm(1, AT_, "Error opening file " + geomFileName);

        //=================== Points =======================
        rcode =
            writeVtuArrayParallel(gfile, &geomPoints(0), geomPointsOffset, geomPointsMemsize, geomPointsMemsizeGlobal);
        if(rcode != MPI_SUCCESS) mTerm(1, AT_, "Error (1) writing to file " + geomFileName);

        //=================== Connectivity =======================
        rcode = writeVtuArrayParallel(gfile, &geomConn(0), geomConnOffset, geomConnMemsize, geomConnMemsizeGlobal);
        if(rcode != MPI_SUCCESS) mTerm(1, AT_, "Error (2) writing to file " + geomFileName);

        //=================== Offsets =======================
        rcode = writeVtuArrayParallel(gfile, &geomOffsets(0), geomOffsetsOffset, geomOffsetsMemsize,
                                      geomOffsetsMemsizeGlobal);
        if(rcode != MPI_SUCCESS) mTerm(1, AT_, "Error (3) writing to file " + geomFileName);

        //=================== BodyId =======================
        rcode = writeVtuArrayParallel(gfile, &bodyId(0), bodyIdOffset, bodyIdMemsize, bodyIdMemsizeGlobal);
        if(rcode != MPI_SUCCESS) mTerm(1, AT_, "Error (4) writing to file " + geomFileName);

        //=================== Velocity =======================
        rcode = writeVtuArrayParallel(gfile, &velocity(0), velocityOffset, velocityMemsize, velocityMemsizeGlobal);
        if(rcode != MPI_SUCCESS) mTerm(1, AT_, "Error (9) writing to file " + geomFileName);

        if(m_solver->m_vtuGeometryOutputExtended) {
          //=================== DomainId =======================
          rcode =
              writeVtuArrayParallel(gfile, &domainIds(0), domainIdsOffset, domainIdsMemsize, domainIdsMemsizeGlobal);
          if(rcode != MPI_SUCCESS) mTerm(1, AT_, "Error (4) writing to file " + geomFileName);

          //=================== Pressure =======================
          rcode = writeVtuArrayParallel(gfile, &pressure(0), pressureOffset, pressureMemsize, pressureMemsizeGlobal);
          if(rcode != MPI_SUCCESS) mTerm(1, AT_, "Error (5) writing to file " + geomFileName);

          //=================== Density =======================
          rcode = writeVtuArrayParallel(gfile, &density(0), densityOffset, densityMemsize, densityMemsizeGlobal);
          if(rcode != MPI_SUCCESS) mTerm(1, AT_, "Error (6) writing to file " + geomFileName);

          //=================== Normals =======================
          rcode = writeVtuArrayParallel(gfile, &normals(0), normalsOffset, normalsMemsize, normalsMemsizeGlobal);
          if(rcode != MPI_SUCCESS) mTerm(1, AT_, "Error (7) writing to file " + geomFileName);

          //=================== Shear =======================
          rcode = writeVtuArrayParallel(gfile, &shear(0), shearOffset, shearMemsize, shearMemsizeGlobal);
          if(rcode != MPI_SUCCESS) mTerm(1, AT_, "Error (8) writing to file " + geomFileName);
        }

        MPI_File_close(&gfile, AT_);
        //###################################################
        //###################################################


        if(domainId() == 0) {
          ofstream ofl;
          ofl.open(geomFileName.c_str(), ios_base::out | ios_base::app);
          if(ofl.is_open() && ofl.good()) {
            ofl << endl;
            ofl << "</AppendedData>" << endl;
            //================== /AppendedData =================

            ofl << "</VTKFile>" << endl;
            ofl.close();
            ofl.clear();
            //================== /VTKFile =================
          } else {
            cerr << "ERROR! COULD NOT OPEN FILE " << geomFileName << " for writing! (3)" << endl;
          }
        }


        //-----------------------
        //------ GEOM.pvd -------
        //-----------------------

        if(domainId() == 0) {
          if(m_solver->m_firstUseWriteVtuOutputParallelGeom) {
            m_solver->m_firstUseWriteVtuOutputParallelGeom = !m_solver->m_restart;
          }
          ofstream ofile;
          ofstream ofile2;
          if(m_solver->m_firstUseWriteVtuOutputParallelGeom)
            ofile.open("out/GEOM.pvd.tmp", ios_base::out | ios_base::trunc);
          else
            ofile.open("out/GEOM.pvd.tmp", ios_base::out | ios_base::app);
          if(ofile.is_open() && ofile.good()) {
            if(m_solver->m_firstUseWriteVtuOutputParallelGeom) {
              ofile << R"(<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">)" << endl;
              ofile << "<Collection>" << endl;
            }
            const size_t pos = geomFileName.find(m_solver->m_solutionOutput);
            MString fname = (pos == string::npos)
                                ? geomFileName
                                : geomFileName.substr(pos + m_solver->m_solutionOutput.size(), string::npos);
            ofile << "<DataSet part=\"" << 0 << "\" timestep=\"" << globalTimeStep << "\" file=\""
                  << "./" << fname << "\"/>" << endl;
            ofile.close();
            ofile.clear();
          } else {
            cerr << "Error opening file out/GEOM.pvd.tmp" << endl;
          }
          m_log << "Command executed with return value: " << system("cp out/GEOM.pvd.tmp out/GEOM.pvd") << " at " << AT_
                << endl;

          ofile2.open("out/GEOM.pvd", ios_base::out | ios_base::app);
          if(ofile2.is_open() && ofile2.good()) {
            ofile2 << "</Collection>" << endl;
            ofile2 << "</VTKFile>" << endl;
            ofile2.close();
            ofile2.clear();
          } else {
            cerr << "Error opening file out/GEOM.pvd" << endl;
          }
          m_solver->m_firstUseWriteVtuOutputParallelGeom = false;
        }
      }

      if(domainId() == 0 && m_solver->m_levelSetMb && m_solver->m_vtuWriteParticleFile) {
        offset = 0;
        ofstream ofl;
        const MString partFileName = m_solver->m_solutionOutput + "POUT_00" + to_string(globalTimeStep) + ".vtp";
        ofl.open(partFileName.c_str(), ios_base::out | ios_base::trunc);
        if(ofl.is_open() && ofl.good()) {
          //================== VTKFile =================
          ofl << "<?xml version=\"1.0\"?>" << endl;
          ofl << R"(<VTKFile type="PolyData" version="1.0" byte_order="LittleEndian" header_type="UInt64">)" << endl;
          ofl << "<PolyData>" << endl;

          //================== FieldData =================
          ofl << "<FieldData>" << endl;
          ofl << "<DataArray type=\"" << uIDataType << R"(" Name="globalTimeStep" format="ascii" NumberOfTuples="1" > )"
              << globalTimeStep << " </DataArray>" << endl;
          if(m_solver->m_outputPhysicalTime) {
            ofl << "<DataArray type=\"" << dataType << R"(" Name="time" format="ascii" NumberOfTuples="1" > )"
                << m_solver->m_physicalTime << " </DataArray>" << endl;
            ofl << "<DataArray type=\"" << dataType << R"(" Name="internalTime" format="ascii" NumberOfTuples="1" > )"
                << m_solver->m_time << " </DataArray>" << endl;
          } else {
            ofl << "<DataArray type=\"" << dataType << R"(" Name="time" format="ascii" NumberOfTuples="1" > )"
                << m_solver->m_time << " </DataArray>" << endl;
            ofl << "<DataArray type=\"" << dataType << R"(" Name="physicalTime" format="ascii" NumberOfTuples="1" > )"
                << m_solver->m_physicalTime << " </DataArray>" << endl;
          }
          ofl << "</FieldData>" << endl;
          //================== /FieldData =================

          //================== Dimensions =================
          ofl << "<Piece NumberOfPoints=\"" << m_solver->m_noEmbeddedBodies << "\" NumberOfVerts=\""
              << m_solver->m_noEmbeddedBodies << "\">" << endl;
          //================== /Dimensions ================

          //================== Points =================
          ofl << "<Points>" << endl;

          ofl << setprecision(12);
          ofl << "<DataArray type=\"" << dataType << R"(" NumberOfComponents="3" format="ascii">)" << endl;
          for(MInt k = 0; k < m_solver->m_noEmbeddedBodies; k++) {
            for(MInt i = 0; i < nDim; i++) {
              ofl << m_solver->m_bodyCenter[k * nDim + i] << " ";
            }
            ofl << endl;
          }
          ofl << "</DataArray>" << endl;
          ofl << "</Points>" << endl;
          //================== /Points =================


          //================== Verts =================
          ofl << "<Verts>" << endl;
          ofl << "<DataArray type=\"" << uIDataType << R"(" Name="connectivity" format="ascii">)" << endl;
          for(MInt k = 0; k < m_solver->m_noEmbeddedBodies; k++) {
            ofl << k << " ";
          }
          ofl << endl;
          ofl << "</DataArray>" << endl;
          ofl << "<DataArray type=\"" << uIDataType << R"(" Name="offsets" format="ascii">)" << endl;
          for(MInt k = 0; k < m_solver->m_noEmbeddedBodies; k++) {
            ofl << k + 1 << " ";
          }
          ofl << endl;
          ofl << "</DataArray>" << endl;
          ofl << "</Verts>" << endl;
          //================== /Verts =================

          //================== PointData =================
          ofl << "<PointData Scalars=\"scalars\">" << endl;
          ofl << "<DataArray type=\"" << dataType << R"(" Name="velocity" NumberOfComponents="3" format="ascii">)"
              << endl;
          for(MInt k = 0; k < m_solver->m_noEmbeddedBodies; k++) {
            for(MInt i = 0; i < nDim; i++) {
              ofl << m_solver->m_bodyVelocity[k * nDim + i] << " ";
            }
            ofl << endl;
          }
          ofl << "</DataArray>" << endl;
          ofl << "</PointData>" << endl;
          //================== /PointData =================

          ofl << "</Piece>" << endl;
          ofl << "</PolyData>" << endl;
          ofl << "</VTKFile>" << endl;
          ofl.close();
          ofl.clear();
          //================== /VTKFile =================
        } else {
          cerr << "ERROR! COULD NOT OPEN FILE " << partFileName << " for writing! (1)" << endl;
        }
      }
    }
  }

  if(facestream != nullptr) {
    delete[] facestream;
    facestream = nullptr;
  }
  if(conn != nullptr) {
    delete[] conn;
    conn = nullptr;
  }

  RECORD_TIMER_STOP(t_geometry);
  RECORD_TIMER_STOP(t_vtutotal);
  DISPLAY_TIMER_INTERM(t_vtutotal);
}
