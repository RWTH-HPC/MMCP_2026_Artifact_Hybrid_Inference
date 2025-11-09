// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "geometryintersection.h"
#include "geometryelement.h"

using namespace std;
using namespace maia;

#define SWAP_ENDIAN

//#ifdef SSADAKSDPLAKSD


// -------------------------------------------------------------------------------

/** \brief adds an additional MC point in the cell
 * sums up all n points passed in indices.
 * points have to be stored in vertices,
 * result is passed in res.
 *
 * \author Claudia Guenther, September 2013
 */
template <MInt nDim_>
void GeometryIntersection<nDim_>::addPoint(const std::vector<polyVertex>* vertices, const MInt* indices, const MInt n,
                                           std::array<MFloat, nDim> res) {
  TRACE();

  res.fill(0.0);

  // sum up all vertex coordinates
  for(MInt v = 0; v < n; v++) {
    for(MInt i = 0; i < nDim_; i++) {
      res[i] += (*vertices)[indices[v]].coordinates[i];
    }
  }

  // determine midpoint (as mean value of the verticies)
  for(MInt i = 0; i < nDim_; i++)
    res[i] /= n;
}

//-----------------------------------------------------------------------------

// computes the face-properties by subdividing any polygon-face into multiple triangles!
template <>
void GeometryIntersection<3>::compFaceIntegrals_pyraBased3(polyFace* face, const std::vector<polyVertex>* vertices,
                                                           MFloat* MC, MFloat* VC, MFloat* XC) {
  // compute Midpoint of face as arithmetic mean of bounding points
  const MInt noPoints = (signed)(*face).vertices.size();
  array<MFloat, 3> faceMidPoint{};

  for(MInt i = 0; i < noPoints; i++) {
    const MInt vertex = (*face).vertices[i];
    for(MInt j = 0; j < nDim; j++) {
      faceMidPoint[j] += (*vertices)[vertex].coordinates[j];
    }
  }

  for(MInt j = 0; j < nDim; j++) {
    faceMidPoint[j] /= noPoints;
  }

  array<MFloat, 3> center{};
  MFloat area = 0.0;

  // compute area, center, tetraeder volume and tetraeder center (with MC) of each face triangle
  if(face->isLine) {
    center = faceMidPoint;
  } else {
    const MInt vertex0 = (*face).vertices[0];
    for(MInt i = 1; i < noPoints - 1; i++) {
      const MInt vertex1 = (*face).vertices[i];
      const MInt vertex2 = (*face).vertices[i + 1];

      MFloat a[3], b[3], c[3];
      MFloat centerTetra[3]{};
      MFloat centerTri[3]{};
      for(MInt j = 0; j < nDim; j++) {
        a[j] = (*vertices)[vertex1].coordinates[j] - (*vertices)[vertex0].coordinates[j];
        b[j] = (*vertices)[vertex2].coordinates[j] - (*vertices)[vertex0].coordinates[j];
        c[j] = MC[j] - (*vertices)[vertex0].coordinates[j];
        centerTri[j] = (*vertices)[vertex0].coordinates[j] + (*vertices)[vertex1].coordinates[j]
                       + (*vertices)[vertex2].coordinates[j];
        centerTetra[j] = centerTri[j] + MC[j];
      }
      MFloat aCrossB[3];
      aCrossB[0] = a[1] * b[2] - a[2] * b[1];
      aCrossB[1] = a[2] * b[0] - a[0] * b[2];
      aCrossB[2] = a[0] * b[1] - a[1] * b[0];
      MFloat aCrossBTimesC = F0;
      MFloat normACrossB = F0;
      for(MInt j = 0; j < nDim; j++) {
        aCrossBTimesC += aCrossB[j] * c[j];
        normACrossB += aCrossB[j] * aCrossB[j];
        centerTri[j] *= F1B3;
        centerTetra[j] *= F1B4;
      }
      normACrossB = sqrt(normACrossB);
      MFloat volumeTetra = F1B6 * aCrossBTimesC;
      MFloat areaTri = F1B2 * normACrossB;

      area += areaTri;
      *VC += volumeTetra;
      for(MInt j = 0; j < nDim; j++) {
        center[j] += areaTri * centerTri[j];
        XC[j] += volumeTetra * centerTetra[j];
      }
    }
    if(area > m_eps) {
      for(MInt j = 0; j < nDim; j++) {
        center[j] /= area;
      }
    } else {
      area = 0.0;
      center = faceMidPoint;
    }
  }
  face->area = area;
  ASSERT(!(std::isnan(area)), "");
  // ASSERT(! (std::isnan(1/area)),"");
  for(MInt j = 0; j < nDim; j++) {
    face->center[j] = center[j];
    ASSERT(!(std::isnan(center[j])), "");
  }
}

//-----------------------------------------------------------------------------

template <>
void GeometryIntersection<3>::compVolumeIntegrals_pyraBased3(std::vector<polyCutCell>* cutCells,
                                                             std::vector<polyFace>* faces,
                                                             const std::vector<polyVertex>* vertices) {
  for(MInt cC = 0; (unsigned)cC < cutCells->size(); cC++) {
    MFloat volume = F0;
    MFloat center[3] = {F0, F0, F0};
    polyCutCell* cutCell = &(*cutCells)[cC];
    MFloat M[3] = {F0, F0, F0};
    MInt count = 0;
    for(MInt i = 0; (unsigned)i < cutCell->faces.size(); i++) {
      const MInt face = cutCell->faces[i];
      for(MInt j = 0; j < (signed)(*faces)[face].vertices.size(); j++) {
        const MInt vertex = (*faces)[face].vertices[j];
        for(MInt k = 0; k < nDim; k++) {
          M[k] += (*vertices)[vertex].coordinates[k];
        }
        count++;
      }
    }
    if(count) {
      for(MInt k = 0; k < nDim; k++) {
        M[k] /= count;
      }
    }
    for(MInt i = 0; (unsigned)i < cutCell->faces.size(); i++) {
      polyFace* face = &(*faces)[cutCell->faces[i]];
      compFaceIntegrals_pyraBased3(face, vertices, M, &volume, center);
    }

    cutCell->volume = volume;
    if(fabs(volume) > 1e2 * m_eps) {
      for(MInt i = 0; i < nDim; i++) {
        cutCell->center[i] = center[i] / volume;
      }
    } else {
      cutCell->volume = fabs(volume);
      for(MInt i = 0; i < nDim; i++) {
        cutCell->center[i] = M[i];
      }
    }
  }
}


//-----------------------------------------------------------------------------

template <>
void GeometryIntersection<3>::writeVTKFileOfCell(MInt cellId, std::vector<polyFace>* faces,
                                                 const std::vector<polyVertex>* vertices, MInt set) {
  const MChar* fileName = "cell_";
  stringstream fileName2;
  fileName2 << fileName << grid().tree().globalId(cellId) << "_s" << set << ".vtk";
  ofstream ofl;
  ofl.open((fileName2.str()).c_str(), ofstream::trunc);

  if(ofl) {
    // set fixed floating point output
    ofl.setf(ios::fixed);
    ofl.precision(16);

    ofl << "# vtk DataFile Version 3.0" << endl
        << "MAIAD cutsurface file" << endl
        << "ASCII" << endl
        << endl
        << "DATASET UNSTRUCTURED_GRID" << endl
        << endl;

    ofl << "POINTS " << (*vertices).size() << " float" << endl;

    for(MInt v = 0; (unsigned)v < (*vertices).size(); v++) {
      for(MInt i = 0; i < 3; i++)
        ofl << (*vertices)[v].coordinates[i] << " ";
      ofl << endl;
    }

    ofl << endl;
    MInt numPoints = 0;
    MInt noFaces = 0;
    for(MInt f = 0; (unsigned)f < (*faces).size(); f++) {
      numPoints += (*faces)[f].vertices.size();
      noFaces++;
    }
    ofl << "CELLS " << noFaces << " " << noFaces + numPoints << endl;

    for(MInt f = 0; (unsigned)f < (*faces).size(); f++) {
      MInt noEdges = (signed)(*faces)[f].vertices.size();
      ofl << noEdges << " ";

      for(MInt i = noEdges - 1; i >= 0; i--) {
        ofl << (*faces)[f].vertices[i] << " ";
      }
      ofl << endl;
    }


    ofl << "CELL_TYPES " << noFaces << endl;
    for(MInt i = 0; i < noFaces; i++) {
      ofl << 7 << endl;
    }

    ofl << "CELL_DATA " << noFaces << endl;
    ofl << "SCALARS bodyId int 1" << endl;
    ofl << "LOOKUP_TABLE default" << endl;
    for(MInt f = 0; (unsigned)f < (*faces).size(); f++) {
      ofl << (*faces)[f].bodyId << endl;
    }

    ofl << "SCALARS face int 1" << endl;
    ofl << "LOOKUP_TABLE default" << endl;
    for(MInt f = 0; (unsigned)f < (*faces).size(); f++) {
      ofl << f << endl;
    }

    ofl << "SCALARS faceId int 1" << endl;
    ofl << "LOOKUP_TABLE default" << endl;
    for(MInt f = 0; (unsigned)f < (*faces).size(); f++) {
      ofl << (*faces)[f].faceId << endl;
    }

    ofl << "SCALARS faceType int 1" << endl;
    ofl << "LOOKUP_TABLE default" << endl;
    for(MInt f = 0; (unsigned)f < (*faces).size(); f++) {
      ofl << (*faces)[f].faceType << endl;
    }

    ofl << "SCALARS cutCell int 1" << endl;
    ofl << "LOOKUP_TABLE default" << endl;
    for(MInt f = 0; (unsigned)f < (*faces).size(); f++) {
      ofl << (*faces)[f].cutCell << endl;
    }

    ofl << "SCALARS isLine int 1" << endl;
    ofl << "LOOKUP_TABLE default" << endl;
    for(MInt f = 0; (unsigned)f < (*faces).size(); f++) {
      ofl << (*faces)[f].isLine << endl;
    }

    ofl << "POINT_DATA " << (*vertices).size() << endl;
    ofl << "SCALARS point int 1" << endl;
    ofl << "LOOKUP_TABLE default" << endl;
    for(MInt v = 0; (unsigned)v < (*vertices).size(); v++) {
      ofl << v << endl;
    }

    ofl << "SCALARS pointId int 1" << endl;
    ofl << "LOOKUP_TABLE default" << endl;
    for(MInt v = 0; (unsigned)v < (*vertices).size(); v++) {
      ofl << (*vertices)[v].pointId << endl;
    }

    ofl << "SCALARS pointType int 1" << endl;
    ofl << "LOOKUP_TABLE default" << endl;
    for(MInt v = 0; (unsigned)v < (*vertices).size(); v++) {
      ofl << (*vertices)[v].pointType << endl;
    }

    ofl.close();
  }
}
//------------------------------------------------------------------------------

template <MInt nDim_>
void GeometryIntersection<nDim_>::writeInfo(std::vector<CutCell<nDim>>& cutCellData, MUint cutc, MInt nameId) {
  MInt cellId = cutCellData[cutc].cellId;

  const MChar* fileName = "cell-Info_";
  stringstream fileName2;
  fileName2 << fileName << grid().tree().globalId(cellId) << "_Id" << nameId << ".txt";

  ofstream ofl;
  ofl.open((fileName2.str()).c_str(), ofstream::trunc);

  const MInt lastDim = nDim - 1;

  if(ofl) {
    ofl.setf(ios::fixed);
    ofl.precision(12);

    ofl << "Writing Cell " << grid().tree().globalId(cellId) << " " << cellId << " " << cutc << " vol "
        << cutCellData[cutc].volume << " SplitChilds " << cutCellData[cutc].noSplitChilds << " "
        << cutCellData[cutc].splitParentId << endl
        << endl;
    for(MInt dir = 0; dir < m_noDirs; dir++) {
      ofl << "Dir " << dir << " external " << cutCellData[cutc].externalFaces[dir] << " faces "
          << cutCellData[cutc].noFacesPerCartesianDir[dir] << endl
          << endl;
    }

    ofl << " Cartesian-Surf: " << cutCellData[cutc].noCartesianSurfaces << endl << endl;

    for(MInt id = 0; id < cutCellData[cutc].noCartesianSurfaces; id++) {
      ofl << "Id " << id << " dir " << cutCellData[cutc].cartFaceDir[id] << " area "
          << cutCellData[cutc].cartFaceArea[id] << " coord " << cutCellData[cutc].cartFaceCentroid[id][0]
          << cutCellData[cutc].cartFaceCentroid[id][1] << cutCellData[cutc].cartFaceCentroid[id][lastDim] << endl
          << endl;
    }

    ofl << " Bndry-Surfaces: " << cutCellData[cutc].noBoundarySurfaces << endl << endl;

    for(MInt id = 0; id < cutCellData[cutc].noBoundarySurfaces; id++) {
      ofl << "Id " << id << " normal " << cutCellData[cutc].boundarySurfaceNormal[id][0]
          << cutCellData[cutc].boundarySurfaceNormal[id][1] << cutCellData[cutc].boundarySurfaceNormal[id][lastDim]
          << " center " << cutCellData[cutc].boundarySurfaceCentroid[id][0]
          << cutCellData[cutc].boundarySurfaceCentroid[id][1] << cutCellData[cutc].boundarySurfaceCentroid[id][lastDim]
          << " area " << cutCellData[cutc].boundarySurfaceArea[id] << endl
          << endl;
    }
  }
}

//------------------------------------------------------------------------------

/**
 * \brief determines the geometric cut-cell data, simple marching cubes version, fast
 * \author Daniel Hartmann, Claudia Guenther, Lennart Schneiders
 * \date
 */
template <>
MBool GeometryIntersection<3>::computeCutFaceSimple(CutCell<nDim>& cutCell) {
  static constexpr MInt edgeCode[24] = {8, 10, 0, 4, 9, 11, 1, 5, 2, 6, 8, 9, 3, 7, 10, 11, 0, 1, 2, 3, 4, 5, 6, 7};
  static constexpr MInt faceCode[24] = {0, 4, 1, 4, 2, 4, 3, 4, 0, 5, 1, 5, 2, 5, 3, 5, 0, 2, 1, 2, 0, 3, 1, 3};
  static constexpr MInt pointCode[6][4] = {{0, 2, 6, 4}, {1, 3, 7, 5}, {0, 1, 5, 4},
                                           {2, 3, 7, 6}, {0, 1, 3, 2}, {4, 5, 7, 6}};
  static constexpr MInt edgeCode2[6][4] = {{0, 10, 4, 8},  {1, 11, 5, 9}, {2, 9, 6, 8},
                                           {3, 11, 7, 10}, {2, 1, 3, 0},  {6, 5, 7, 4}};
  static constexpr MFloat signStencil[8][3] = {{-F1, -F1, -F1}, {F1, -F1, -F1}, {-F1, F1, -F1}, {F1, F1, -F1},
                                               {-F1, -F1, F1},  {F1, -F1, F1},  {-F1, F1, F1},  {F1, F1, F1}};

  static constexpr MBool improvedCentroidComputation = false; // true;
  // new scheme for more accurate face and volume centroids, Lennart

  const MInt noEdges = CutCell<nDim>::noEdges;
  MBool cutFaces[m_noDirs];
  ScratchSpace<MBool> cutEdgesAll(noEdges, AT_, "cutEdgesAll");
  ScratchSpace<MBool> nonFluidSideEdge(noEdges, AT_, "nonFluidSideEdge");
  MBool shapeIsInside[m_noDirs];
  MBool positiveSlope;
  MBool pointInside;
  MBool pointBRInside;
  MBool pointTLInside;
  const MFloat epsNormal = 1e-15;
  MInt corner;
  MInt edgeCounter;
  MInt face0, face1;
  MInt face0Space, face0Side;
  MInt noCutFaces;
  MInt side0, side1;
  MInt sideId;
  MInt spaceId, spaceId1, spaceId2;
  MInt space0, space1, space0Relative, space1Relative;
  MFloat h;
  ScratchSpace<MInt> cutPointOnEdge(noEdges, AT_, "cutPointOnEdge");
  MInt noNonFluidEdges[m_noDirs];
  MInt cutEdges[4];
  MInt cutPoint[4];
  MFloat cc0[nDim];
  MFloat cc1[nDim];
  MFloat cellHalfLength[nDim];
  MFloat cellLength[nDim];
  ScratchSpace<MFloat> cutLineLength(noEdges, AT_, "cutLineLength");
  MFloat faceVolume[m_noDirs];
  MFloat gridFaceVolume[m_noDirs];
  MFloat cutOutVolume;
  MFloat dA;
  MFloat dfcc[2];
  MFloat eps;
  MFloat Fvol;
  MFloat minX, maxX, minY, maxY;
  MFloat negativePoint[3] = {0, 0, 0};
  MFloat oppositePoint[3] = {0, 0, 0};
  MFloat point[3] = {0, 0, 0};
  MFloat pointP[3] = {0, 0, 0};
  MFloat pointM[3] = {0, 0, 0};
  MFloat pointW[3] = {0, 0, 0};
  MFloat pointE[3] = {0, 0, 0};
  MFloat rectVolume, triVolume;
  MFloat trianglePoint[3] = {0, 0, 0};
  MFloat trianglePointP[3] = {0, 0, 0};
  MFloat trianglePointM[3] = {0, 0, 0};
  MFloat trianglePointW[3] = {0, 0, 0};
  MFloat trianglePointE[3] = {0, 0, 0};
  MFloat cutLineCentroid[6][3];
  MFloat delta[6][3];
  MFloat faceCentroid[6][3];
  MFloat triangleCentroid[6][3];
  MFloat normalVector[6][3];

  MInt cellId = cutCell.cellId;


  for(MInt f = 0; f < m_noDirs; f++) {
    cutCell.externalFaces[f] = false;
    cutCell.noFacesPerCartesianDir[f] = 1;
  }

  noCutFaces = 0;
  cutCell.noBoundarySurfaces = 1;

  face0 = -1;
  face1 = -1;

  // reset all local variables
  for(MInt i = 0; i < noEdges; i++) {
    cutEdgesAll[i] = false;
    nonFluidSideEdge[i] = false;
    cutPointOnEdge[i] = -1;
  }
  for(MInt i = 0; i < m_noDirs; i++) {
    cutFaces[i] = false;
    shapeIsInside[i] = false;
  }

  // determine cell geometry
  for(MInt i = 0; i < nDim; i++) {
    cellLength[i] = a_cellLengthAtCell(cellId);
    cellHalfLength[i] = F1B2 * cellLength[i];
  }
  eps = F0;

  // store all cut points in a separate array
  for(MInt cp = 0; cp < cutCell.noCutPoints; cp++) {
    cutPointOnEdge[cutCell.cutEdges[cp]] = cp;
  }


  cutCell.noCartesianSurfaces = 0;
  cutCell.noTotalFaces = 0;

  // I.
  // --
  //    determine the cut geometry for each face
  //    assuming that each face has either 0 or 2 cut points
  for(MInt face = 0; face < m_noDirs; face++) {
    // assemble face edges and cut points
    MInt n = 0;
    for(MInt e = 0; e < 4; e++) {
      MInt cp = cutPointOnEdge[edgeCode[4 * face + e]];
      if(cp > -1) {
        cutPoint[n] = cp;
        cutEdges[n] = e;
        cutEdgesAll[edgeCode[4 * face + e]] = true;
        n++;
      }
    }

    if(n > 0) {
      cutFaces[face] = true;
      noCutFaces++;
    }
    MBool faceIsActive = true;
    if(n == 1 || n > 2) {
      cerr << "** ERROR create cut face: " << n << " cut points for face " << face << " body "
           << cutCell.associatedBodyIds[0] << endl;
      stringstream errorMessage;
      errorMessage << "** ERROR create cut face: " << n << " cut points for face " << face << endl;
      cerr << "cell " << cellId << " cog " << a_coordinate(cellId, 0) << " " << a_coordinate(cellId, 1) << " "
           << a_coordinate(cellId, 2) << endl;
      return false;
    } else {
      spaceId = face / 2;
      sideId = face % 2;
      spaceId1 = (spaceId + 1) % nDim;
      spaceId2 = (spaceId1 + 1) % nDim;

      gridFaceVolume[face] = cellLength[spaceId1] * cellLength[spaceId2];
      faceVolume[face] = gridFaceVolume[face];

      // set the reference points
      for(MInt i = 0; i < nDim; i++) {
        point[i] = a_coordinate(cellId, i) - cellHalfLength[i];
        negativePoint[i] = point[i];
        oppositePoint[i] = point[i] + cellLength[i];
      }

      // set the secondary coordinate (spaceId) of the reference points
      point[spaceId] += (MFloat)sideId * cellLength[spaceId];
      oppositePoint[spaceId] = point[spaceId];

      // inside / outside check for point
      for(MInt i = 0; i < nDim; i++) {
        pointP[i] = point[i];
        pointM[i] = point[i];
        pointW[i] = point[i];
        pointE[i] = point[i];
      }
      pointP[spaceId1] += eps;
      pointM[spaceId1] -= eps;
      pointW[spaceId2] += eps;
      pointE[spaceId2] -= eps;

      corner = 0;
      if(sideId) corner = IPOW2(spaceId);

      pointInside = cutCell.cornerIsInsideGeometry[0][corner] != 0;

      MInt cornerBR = corner + IPOW2(spaceId1);
      MInt cornerTL = corner + IPOW2(spaceId2);
      MInt cornerOP = corner + IPOW2(spaceId1) + IPOW2(spaceId2);

      pointBRInside = cutCell.cornerIsInsideGeometry[0][cornerBR];
      pointTLInside = cutCell.cornerIsInsideGeometry[0][cornerTL];


      const MInt pointIds[4] = {corner, cornerBR, cornerOP, cornerTL};

      faceIsActive = true;
      if(n == 0) {
        if(cutCell.cornerIsInsideGeometry[0][pointIds[0]] && cutCell.cornerIsInsideGeometry[0][pointIds[1]]
           && cutCell.cornerIsInsideGeometry[0][pointIds[2]] && cutCell.cornerIsInsideGeometry[0][pointIds[3]]) {
          faceIsActive = false;
        }
      }

      faceCentroid[face][spaceId] = point[spaceId];

      if(n == 2) {
        // delta of cut point coordinates, sign of the cut line slope
        for(MInt i = 0; i < nDim; i++) {
          delta[face][i] = cutCell.cutPoints[cutPoint[0]][i] - cutCell.cutPoints[cutPoint[1]][i];
        }
        positiveSlope = delta[face][spaceId1] * delta[face][spaceId2] > 0;

        for(MInt i = 0; i < nDim; i++) {
          delta[face][i] = fabs(delta[face][i]);
        }

        // compute the length of the cut line
        cutLineLength[face] = sqrt(POW2(delta[face][spaceId1]) + POW2(delta[face][spaceId2]));

        // compute the coordinates of the cut line centroid
        for(MInt i = 0; i < nDim; i++) {
          cutLineCentroid[face][i] = F1B2 * (cutCell.cutPoints[cutPoint[0]][i] + cutCell.cutPoints[cutPoint[1]][i]);
        }

        // compute the volume and the normal vector of the cut face
        // assuming 2 cut points!!! (checked above)
        // space0: space Id of the first cut edge
        // space1: space Id of the second cut edge
        // *Relative: in two-dimensional space0-space1 space, either 0 or 1
        space0Relative = cutEdges[0] / 2;
        space1Relative = cutEdges[1] / 2;
        space0 = (space0Relative + spaceId + 1) % nDim;
        side0 = cutEdges[0] % 2;
        space1 = (space1Relative + spaceId + 1) % nDim;
        side1 = cutEdges[1] % 2;

        // compute the absolute values of the normal vector components
        normalVector[face][spaceId] = 0;
        normalVector[face][spaceId1] = delta[face][spaceId2] / cutLineLength[face];
        normalVector[face][spaceId2] = delta[face][spaceId1] / cutLineLength[face];


        // determine the cut geometry and compute
        //  - the volume of the boundary cell
        //  - the coordinate shift
        //  - nonFluidSides
        if(space0 != space1) {
          // 1. the cut geometry is a triangle
          faceVolume[face] = F1B2 * (delta[face][spaceId1] * delta[face][spaceId2]);

          trianglePoint[space0] = point[space0] + cellLength[space0] * (MFloat)side0;
          trianglePoint[space1] = point[space1] + cellLength[space1] * (MFloat)side1;
          trianglePoint[spaceId] = point[spaceId];

          for(MInt i = 0; i < nDim; i++) {
            trianglePointP[i] = trianglePoint[i];
            trianglePointM[i] = trianglePoint[i];
            trianglePointW[i] = trianglePoint[i];
            trianglePointE[i] = trianglePoint[i];
          }
          trianglePointP[spaceId1] += eps;
          trianglePointM[spaceId1] -= eps;
          trianglePointW[spaceId2] += eps;
          trianglePointE[spaceId2] -= eps;

          faceCentroid[face][space0] = trianglePoint[space0] + F2B3 * (F1B2 - (MFloat)side0) * delta[face][space0];
          faceCentroid[face][space1] = trianglePoint[space1] + F2B3 * (F1B2 - (MFloat)side1) * delta[face][space1];

          for(MInt i = 0; i < nDim; i++) {
            triangleCentroid[face][i] = faceCentroid[face][i];
          }


          MInt triangleCorner = 0;
          triangleCorner += sideId * IPOW2(spaceId);
          triangleCorner += side0 * IPOW2(space0);
          triangleCorner += side1 * IPOW2(space1);


          // check if the third triangle point is inside or outside

          if(cutCell.cornerIsInsideGeometry[0][triangleCorner]) {
            // computed volume is inside the geometry and cut out
            shapeIsInside[face] = true;
            cutOutVolume = faceVolume[face];
            faceVolume[face] = gridFaceVolume[face] - faceVolume[face];
            faceCentroid[face][space0] =
                F1 / faceVolume[face]
                * ((gridFaceVolume[face] * a_coordinate(cellId, space0)) - (cutOutVolume * faceCentroid[face][space0]));
            faceCentroid[face][space1] =
                F1 / faceVolume[face]
                * ((gridFaceVolume[face] * a_coordinate(cellId, space1)) - (cutOutVolume * faceCentroid[face][space1]));

            // no non-fluid edge
            noNonFluidEdges[face] = 0;
            // DEBUG (not required)
            for(MInt e = 0; e < 4; e++) {
              if(nonFluidSideEdge[edgeCode[4 * face + e]]) {
                cerr << "Error non-fluid edge was true for edge " << edgeCode[4 * face + e] << " face " << face << endl;
                cerr << "All face data " << noNonFluidEdges[0] << noNonFluidEdges[1] << noNonFluidEdges[2]
                     << noNonFluidEdges[3] << noNonFluidEdges[4] << noNonFluidEdges[5] << endl;
                cerr << "All edge data " << nonFluidSideEdge[edgeCode[0]] << nonFluidSideEdge[edgeCode[1]]
                     << nonFluidSideEdge[edgeCode[2]] << nonFluidSideEdge[edgeCode[3]] << " "
                     << nonFluidSideEdge[edgeCode[4]] << nonFluidSideEdge[edgeCode[5]] << nonFluidSideEdge[edgeCode[6]]
                     << nonFluidSideEdge[edgeCode[7]] << " " << nonFluidSideEdge[edgeCode[8]]
                     << nonFluidSideEdge[edgeCode[9]] << nonFluidSideEdge[edgeCode[10]]
                     << nonFluidSideEdge[edgeCode[11]] << endl;
                cerr << "cell " << cellId << " " << a_coordinate(cellId, 0) << " " << a_coordinate(cellId, 1) << " "
                     << a_coordinate(cellId, 2) << endl;
                cerr << " " << endl;
              }
              nonFluidSideEdge[edgeCode[4 * face + e]] = false;
            }
            // END DEBUG

          } else {
            // computed volume is the boundary cell volume
            noNonFluidEdges[face] = 2;
            nonFluidSideEdge[edgeCode[4 * face + 2 * space0Relative + ((side0 + 1) % 2)]] = true;
            nonFluidSideEdge[edgeCode[4 * face + 2 * space1Relative + ((side1 + 1) % 2)]] = true;
          }

          // determine the sign of the normal vector components
          if(!pointInside) {
            if(pointTLInside) {
              if(pointBRInside) {
                // bottom right corner is inside the geometry
                normalVector[face][spaceId1] = -normalVector[face][spaceId1];
                normalVector[face][spaceId2] = -normalVector[face][spaceId2];
              } else {
                // top left corner is inside the geometry
                normalVector[face][spaceId2] = -normalVector[face][spaceId2];
              }
            } else {
              if(pointBRInside) {
                // bottom right corner is inside the geometry
                normalVector[face][spaceId1] = -normalVector[face][spaceId1];
              } else {
                // top right corner is inside the geometry
                normalVector[face][spaceId1] = -normalVector[face][spaceId1];
                normalVector[face][spaceId2] = -normalVector[face][spaceId2];
              }
            }
          } else {
            if(pointTLInside) {
              if(pointBRInside) {
                // top right corner is outside the geometry
              } else {
                // bottom right corner is outside the geometry
                normalVector[face][spaceId2] = -normalVector[face][spaceId2];
              }
            } else {
              if(pointBRInside) {
                // top left corner is outside the geometry
                normalVector[face][spaceId1] = -normalVector[face][spaceId1];
              } else {
                // bottom left corner is inside the geometry
              }
            }
          }
        } else {
          // 2. the cut geometry is a trapezoid
          noNonFluidEdges[face] = 1;
          if(space0 == spaceId1) {
            // 2.a) the cut face is along the x_spaceId1 coordinate
            faceVolume[face] = cellLength[spaceId1] * (cutLineCentroid[face][spaceId2] - point[spaceId2]);

            // check if point is inside or outside
            if(pointInside) {
              faceVolume[face] = gridFaceVolume[face] - faceVolume[face];
              nonFluidSideEdge[edgeCode[4 * face + 2]] = true;
              // * compute boundary cell center
              // * set signs of the normal vector
              //    - x_0 sign according to the slope (+/-,-/+)
              //    - x_1 sign is positive
              maxY = cutLineCentroid[face][spaceId2] + F1B2 * delta[face][spaceId2];
              cc0[0] = a_coordinate(cellId, spaceId1);
              cc0[1] = F1B2 * (oppositePoint[spaceId2] + maxY);
              if(positiveSlope) {
                cc1[0] = point[spaceId1] + F1B3 * cellLength[spaceId1];
                normalVector[face][spaceId1] = -normalVector[face][spaceId1];
              } else {
                cc1[0] = point[spaceId1] + F2B3 * cellLength[spaceId1];
              }
              cc1[1] = maxY - F1B3 * delta[face][spaceId2];
              rectVolume = (oppositePoint[spaceId2] - maxY) * cellLength[spaceId1];
              triVolume = faceVolume[face] - rectVolume;
            } else {
              nonFluidSideEdge[edgeCode[4 * face + 3]] = true;
              // * compute boundary cell center
              // * set signs of the normal vector
              //    - x_0 sign according to the slope (+/+,-/-)
              //    - x_1 sign is negative
              normalVector[face][spaceId2] = -normalVector[face][spaceId2];
              minY = cutLineCentroid[face][spaceId2] - F1B2 * delta[face][spaceId2];
              cc0[0] = a_coordinate(cellId, spaceId1);
              cc0[1] = F1B2 * (point[spaceId2] + minY);
              if(positiveSlope) {
                cc1[0] = point[spaceId1] + F2B3 * cellLength[spaceId1];
              } else {
                normalVector[face][spaceId1] = -normalVector[face][spaceId1];
                cc1[0] = point[spaceId1] + F1B3 * cellLength[spaceId1];
              }
              cc1[1] = minY + F1B3 * delta[face][spaceId2];
              rectVolume = (minY - point[spaceId2]) * cellLength[spaceId1];
              triVolume = faceVolume[face] - rectVolume;
            }
          } else {
            // 2.b) the cut face is along the x_spaceId2 coordinate
            faceVolume[face] = cellLength[spaceId2] * (cutLineCentroid[face][spaceId1] - point[spaceId1]);

            // check if point is inside or outside
            if(pointInside) {
              faceVolume[face] = gridFaceVolume[face] - faceVolume[face];
              nonFluidSideEdge[edgeCode[4 * face]] = true;
              // * compute boundary cell center
              // * set signs of the normal vector
              //    - x_0 sign is positive
              //    - x_1 sign according to the slope (+/-,-/+)
              maxX = cutLineCentroid[face][spaceId1] + F1B2 * delta[face][spaceId1];
              cc0[0] = F1B2 * (oppositePoint[spaceId1] + maxX);
              cc0[1] = a_coordinate(cellId, spaceId2);
              cc1[0] = maxX - F1B3 * delta[face][spaceId1];
              if(positiveSlope) {
                normalVector[face][spaceId2] = -normalVector[face][spaceId2];
                cc1[1] = point[spaceId2] + F1B3 * cellLength[spaceId2];
              } else {
                cc1[1] = point[spaceId2] + F2B3 * cellLength[spaceId2];
              }
              rectVolume = (oppositePoint[spaceId1] - maxX) * cellLength[spaceId2];
              triVolume = faceVolume[face] - rectVolume;
            } else {
              nonFluidSideEdge[edgeCode[4 * face + 1]] = true;
              // * compute boundary cell center
              // * set signs of the normal vector
              //    - x_0 sign is negative
              //    - x_1 sign according to the slope (+/+,-/-)
              normalVector[face][spaceId1] = -normalVector[face][spaceId1];
              minX = cutLineCentroid[face][spaceId1] - F1B2 * delta[face][spaceId1];
              cc0[0] = F1B2 * (point[spaceId1] + minX);
              cc0[1] = a_coordinate(cellId, spaceId2);
              cc1[0] = minX + F1B3 * delta[face][spaceId1];
              if(positiveSlope) {
                cc1[1] = point[spaceId2] + F2B3 * cellLength[spaceId2];
              } else {
                normalVector[face][spaceId2] = -normalVector[face][spaceId2];
                cc1[1] = point[spaceId2] + F1B3 * cellLength[spaceId2];
              }
              rectVolume = (minX - point[spaceId1]) * cellLength[spaceId2];
              triVolume = faceVolume[face] - rectVolume;
            }
          }

          faceCentroid[face][spaceId1] = (rectVolume * cc0[0] + triVolume * cc1[0]) / faceVolume[face];
          faceCentroid[face][spaceId2] = (rectVolume * cc0[1] + triVolume * cc1[1]) / faceVolume[face];
        }

        for(MInt i = 0; i < nDim; i++) {
          cutCell.cartFaceCentroid[cutCell.noCartesianSurfaces][i] = faceCentroid[face][i];
        }
        cutCell.cartFaceArea[cutCell.noCartesianSurfaces] = faceVolume[face];
        cutCell.cartFaceDir[cutCell.noCartesianSurfaces] = face;
        cutCell.noCartesianSurfaces++;

      } else {
        // set the nonFluidSideEdge flag
        // none of the face edges is cut
        if(pointInside) {
          for(MInt e = 0; e < 4; e++) {
            nonFluidSideEdge[edgeCode[4 * face + e]] = true;
          }
        }

        if(faceIsActive) {
          for(MInt i = 0; i < nDim; i++) {
            cutCell.cartFaceCentroid[cutCell.noCartesianSurfaces][i] = a_coordinate(cellId, i);
          }
          if(face % 2 == 0) {
            cutCell.cartFaceCentroid[cutCell.noCartesianSurfaces][face / 2] -= cellHalfLength[0];
          } else {
            cutCell.cartFaceCentroid[cutCell.noCartesianSurfaces][face / 2] += cellHalfLength[0];
          }
          cutCell.cartFaceArea[cutCell.noCartesianSurfaces] = POW2(cellLength[0]);
          cutCell.cartFaceDir[cutCell.noCartesianSurfaces] = face;
          cutCell.noCartesianSurfaces++;
        }
      }
    }
  }


  // from here on the whole control volume is concerned

  // II.
  // --
  //    set the non fluid side Ids
  for(MInt face = 0; face < m_noDirs; face++) {
    MBool nonFluidSide = true;
    for(MInt e = 0; e < 4; e++) {
      if(!nonFluidSideEdge[edgeCode[4 * face + e]]) {
        nonFluidSide = false;
        break;
      }
    }
    if(nonFluidSide) {
      cutCell.externalFaces[face] = true;
      cutCell.noFacesPerCartesianDir[face] = 0;
      faceVolume[face] = F0;
    }
  }


  // III.
  // --
  //    determine the 3D cell geometry
  switch(noCutFaces) {
    case 0: {
      cerr << "case 0 " << endl;
      cerr << a_coordinate(cellId, 0) << " " << a_coordinate(cellId, 1) << " " << a_coordinate(cellId, 2) << " "
           << cellHalfLength[0] << " " << grid().tree().level(cellId) << endl;
      return false;
      // mTerm(1,AT_, "case 0");
    }

    case 3: {
      // III. 1.
      // the cut geometry is a tetraeder

      // find the cut sides
      face0 = -1;
      face1 = -1;
      for(MInt face = 0; face < m_noDirs; face++) {
        if(cutFaces[face]) {
          if(face0 == -1) {
            face0 = face;
          } else {
            if(face1 == -1) {
              face1 = face;
            }
          }
        }
      }
      edgeCounter = 0;
      // find the two cut edges of face0
      for(MInt e = 0; e < 4; e++) {
        cutEdges[e] = 0;
        if(cutEdgesAll[edgeCode[4 * face0 + e]]) {
          cutEdges[edgeCounter] = e;
          edgeCounter++;
        }
      }

      // DEBUG
      if(edgeCounter != 2) cerr << "ERROR create cut face, tetraeter shape, not enough cut edges" << endl;
      // find the remaining cut edge (face1)
      for(MInt e = 0; e < 4; e++) {
        if(cutEdgesAll[edgeCode[4 * face1 + e]]) {
          if(faceCode[2 * edgeCode[4 * face1 + e]] != face0 && faceCode[2 * edgeCode[4 * face1 + e] + 1] != face0) {
            cutEdges[2] = e;
          }
        }
      }

      // determine the cut point on the remaining cut edge
      // for( MInt cp = 0; cp < m_fvBndryCnd->m_bndryCells->a[ bndryId ].m_srfcs[0]->m_noCutPoints; cp++ ) {
      for(MInt cp = 0; cp < cutCell.noCutPoints; cp++) {
        if(cutCell.cutEdges[cp] == edgeCode[4 * face1 + cutEdges[2]]) {
          cutPoint[0] = cp;
        }
      }

      face0Space = face0 / 2;
      face0Side = face0 % 2;

      // determine the center and the volume of the tetraeder
      // center of the triangle
      for(MInt i = 0; i < nDim; i++) {
        cutCell.volumetricCentroid[i] = F3B4 * triangleCentroid[face0][i] + F1B4 * cutCell.cutPoints[cutPoint[0]][i];
      }

      // volume
      h = F2 * (F1B2 - (MFloat)face0Side)
          * (cutCell.cutPoints[cutPoint[0]][face0Space]
             - (negativePoint[face0Space] + (MFloat)face0Side * cellLength[face0Space]));
      if(shapeIsInside[face0]) {
        cutCell.volume = F1B3 * (gridFaceVolume[face0] - faceVolume[face0]) * h;
      } else {
        cutCell.volume = F1B3 * faceVolume[face0] * h;
      }

      // correct the volume and the cell center of the reshaped cell if the geometry is cut out
      if(shapeIsInside[face0]) {
        cutOutVolume = cutCell.volume;
        cutCell.volume = a_gridCellVolume(cellId) - cutCell.volume;
        for(MInt i = 0; i < nDim; i++) {
          cutCell.volumetricCentroid[i] =
              F1 / mMax(m_eps, cutCell.volume)
              * (a_gridCellVolume(cellId) * a_coordinate(cellId, i) - cutOutVolume * cutCell.volumetricCentroid[i]);
        }
      }
      // compute the coordinate shift
      for(MInt i = 0; i < nDim; i++) {
        cutCell.volumetricCentroid[i] -= a_coordinate(cellId, i);
      }

      // compute the cut surface area (store surface components in the normal vector member)
      for(MInt i = 0; i < nDim; i++) {
        cutCell.boundarySurfaceNormal[0][i] = faceVolume[2 * i + 1] - faceVolume[2 * i];
      }
      cutCell.boundarySurfaceArea[0] =
          sqrt(POW2(cutCell.boundarySurfaceNormal[0][0]) + POW2(cutCell.boundarySurfaceNormal[0][1])
               + POW2(cutCell.boundarySurfaceNormal[0][2]));

      cutCell.boundarySurfaceArea[0] = mMax(epsNormal, cutCell.boundarySurfaceArea[0]);

      cutCell.boundarySurfaceBndryCndId[0] = -1;
      cutCell.boundarySurfaceBodyId[0] = cutCell.cutBodyIds[0];

      // compute the normal vector
      for(MInt i = 0; i < nDim; i++) {
        cutCell.boundarySurfaceNormal[0][i] /= cutCell.boundarySurfaceArea[0];
      }

      // compute the cut surface centroid
      for(MInt i = 0; i < nDim; i++) {
        cutCell.boundarySurfaceCentroid[0][i] =
            F2B3 * cutLineCentroid[face0][i] + F1B3 * cutCell.cutPoints[cutPoint[0]][i];
      }
      break;
    }
    case 4:
    case 5:
    case 6: {
      // III. 2.

      // find a pair of opposite cut sides (the one with the largest area)
      cutLineLength[6] = F0;
      for(MInt face = 0; face < m_noDirs; face += 2) {
        if(cutFaces[face]) {
          if(cutFaces[face + 1]) {
            cutLineLength[7] = cutLineLength[face] + cutLineLength[face + 1];
            if(cutLineLength[7] > cutLineLength[6]) {
              face0 = face;
              face1 = face + 1;
              cutLineLength[6] = cutLineLength[7];
            }
          }
        }
      }
      if(face0 < 0 || face1 < 0) {
        cerr << "Error: opposite cut faces not found" << endl;
        return false;
      }

      spaceId = face0 / 2;
      spaceId1 = (spaceId + 1) % nDim;
      spaceId2 = (spaceId1 + 1) % nDim;

      // compute the new cell center
      //   coordinate spaceId1 and Id2
      Fvol = F1 / (faceVolume[face0] + faceVolume[face1]);
      dA = faceVolume[face1] - faceVolume[face0];
      dfcc[0] = faceCentroid[face1][spaceId1] - faceCentroid[face0][spaceId1];
      dfcc[1] = faceCentroid[face1][spaceId2] - faceCentroid[face0][spaceId2];

      cutCell.volumetricCentroid[spaceId] =
          faceCentroid[face0][spaceId] + F1B3 * (F1 + Fvol * faceVolume[face1]) * cellLength[spaceId];

      cutCell.volumetricCentroid[spaceId1] =
          F2 * Fvol
          * (faceVolume[face0] * faceCentroid[face0][spaceId1]
             + F1B2 * (faceVolume[face0] * dfcc[0] + dA * faceCentroid[face0][spaceId1]) + F1B3 * dA * dfcc[0]);

      cutCell.volumetricCentroid[spaceId2] =
          F2 * Fvol
          * (faceVolume[face0] * faceCentroid[face0][spaceId2]
             + F1B2 * (faceVolume[face0] * dfcc[1] + dA * faceCentroid[face0][spaceId2]) + F1B3 * dA * dfcc[1]);

      for(MInt i = 0; i < nDim; i++) {
        cutCell.volumetricCentroid[i] -= a_coordinate(cellId, i);
      }

      // compute the cut surface area (store surface components in the normal vector member)
      for(MInt i = 0; i < nDim; i++) {
        cutCell.boundarySurfaceNormal[0][i] = faceVolume[2 * i + 1] - faceVolume[2 * i];
      }
      cutCell.boundarySurfaceArea[0] =
          sqrt(POW2(cutCell.boundarySurfaceNormal[0][0]) + POW2(cutCell.boundarySurfaceNormal[0][1])
               + POW2(cutCell.boundarySurfaceNormal[0][2]));

      cutCell.boundarySurfaceArea[0] = mMax(epsNormal, cutCell.boundarySurfaceArea[0]);

      cutCell.boundarySurfaceBndryCndId[0] = -1;
      cutCell.boundarySurfaceBodyId[0] = cutCell.cutBodyIds[0];

      // compute the normal vector
      for(MInt i = 0; i < nDim; i++) {
        cutCell.boundarySurfaceNormal[0][i] /= cutCell.boundarySurfaceArea[0];
      }

      // compute the cut surface centroid
      for(MInt i = 0; i < nDim; i++) {
        cutCell.boundarySurfaceCentroid[0][i] = F1B2 * (cutLineCentroid[face0][i] + cutLineCentroid[face1][i]);
      }

      // compute the boundary cell volume using Gauss theorem
      cutCell.volume = -F1B3
                       * (faceVolume[0] * faceCentroid[0][0] - faceVolume[1] * faceCentroid[1][0]
                          + faceVolume[2] * faceCentroid[2][1] - faceVolume[3] * faceCentroid[3][1]
                          + faceVolume[4] * faceCentroid[4][2] - faceVolume[5] * faceCentroid[5][2]
                          + cutCell.boundarySurfaceArea[0]
                                * (cutCell.boundarySurfaceNormal[0][0] * cutCell.boundarySurfaceCentroid[0][0]
                                   + cutCell.boundarySurfaceNormal[0][1] * cutCell.boundarySurfaceCentroid[0][1]
                                   + cutCell.boundarySurfaceNormal[0][2] * cutCell.boundarySurfaceCentroid[0][2]));

      break;
    }
    default: {
      cerr << "this case has not yet been implemented" << endl;
      cerr << "cell " << cellId << endl;
      cerr << "coordinates"
           << " " << a_coordinate(cellId, 0) << " " << a_coordinate(cellId, 1) << " " << a_coordinate(cellId, 2)
           << endl;
      cerr << "number of cut faces: " << noCutFaces << endl;
      return false;
      // mTerm(1,AT_, "This case has not yet been implemented");
    }
  }

  multimap<MFloat, MInt> sortedCPs;
  MFloat pCoords[nDim];
  MFloat vec_a[nDim];
  MFloat vec_b[nDim];
  MFloat normal[nDim];
  MFloat acentre[nDim] = {F0, F0, F0};
  for(MInt cp = 0; cp < cutCell.noCutPoints; cp++) {
    for(MInt i = 0; i < nDim; i++) {
      acentre[i] += cutCell.cutPoints[cp][i];
    }
  }
  for(MInt i = 0; i < nDim; i++) {
    acentre[i] /= (MFloat)cutCell.noCutPoints;
  }
  for(MInt i = 0; i < nDim; i++) {
    normal[i] = cutCell.boundarySurfaceNormal[0][i];
  }
  spaceId = 0;
  MFloat maxC = fabs(normal[0]);
  for(MInt i = 1; i < nDim; i++) {
    if(fabs(normal[i]) > maxC) {
      maxC = fabs(normal[i]);
      spaceId = i;
    }
  }
  spaceId1 = (spaceId + 1) % nDim;
  spaceId2 = (spaceId1 + 1) % nDim;
  vec_a[spaceId1] = F1;
  vec_a[spaceId2] = F1;
  vec_a[spaceId] = -(vec_a[spaceId1] * normal[spaceId1] + vec_a[spaceId2] * normal[spaceId2]) / normal[spaceId];
  MFloat vecsum = sqrt(POW2(vec_a[0]) + POW2(vec_a[1]) + POW2(vec_a[2]));
  for(MInt i = 0; i < nDim; i++) {
    vec_a[i] /= vecsum;
  }
  vec_b[spaceId] = normal[spaceId1] * vec_a[spaceId2] - normal[spaceId2] * vec_a[spaceId1];
  vec_b[spaceId1] = normal[spaceId2] * vec_a[spaceId] - normal[spaceId] * vec_a[spaceId2];
  vec_b[spaceId2] = normal[spaceId] * vec_a[spaceId1] - normal[spaceId1] * vec_a[spaceId];
  vecsum = sqrt(POW2(vec_b[0]) + POW2(vec_b[1]) + POW2(vec_b[2]));
  for(MInt i = 0; i < nDim; i++) {
    vec_b[i] /= vecsum;
  }
  for(MInt cp = 0; cp < cutCell.noCutPoints; cp++) {
    for(MInt i = 0; i < nDim; i++) {
      pCoords[i] = cutCell.cutPoints[cp][i];
    }
    MFloat dx = F0;
    MFloat dy = F0;
    for(MInt i = 0; i < nDim; i++) {
      dx += vec_a[i] * (pCoords[i] - acentre[i]);
      dy += vec_b[i] * (pCoords[i] - acentre[i]);
    }
    sortedCPs.insert(make_pair(atan2(dy, dx), cp));
  }
  ASSERT((signed)sortedCPs.size() == cutCell.noCutPoints, "");

  cutCell.allFacesNoPoints[cutCell.noTotalFaces] = 0;
  cutCell.allFacesBodyId[cutCell.noTotalFaces] = cutCell.cutBodyIds[0];
  for(auto& sortedCP : sortedCPs) {
    cutCell.allFacesPointIds[cutCell.noTotalFaces][cutCell.allFacesNoPoints[cutCell.noTotalFaces]] =
        CutCell<nDim>::noCorners + sortedCP.second;
    cutCell.allFacesNoPoints[cutCell.noTotalFaces]++;
  }
  cutCell.noTotalFaces++;

  for(MInt face = 0; face < m_noDirs; face++) {
    if(!cutCell.externalFaces[face]) {
      cutCell.allFacesNoPoints[cutCell.noTotalFaces] = 0;
      cutCell.allFacesBodyId[cutCell.noTotalFaces] = -1;
      for(MInt p = 0; p < 4; p++) {
        if(!cutCell.cornerIsInsideGeometry[0][pointCode[face][p]]) {
          cutCell.allFacesPointIds[cutCell.noTotalFaces][cutCell.allFacesNoPoints[cutCell.noTotalFaces]] =
              pointCode[face][p];
          cutCell.allFacesNoPoints[cutCell.noTotalFaces]++;
        }
        for(MInt cp = 0; cp < cutCell.noCutPoints; cp++) {
          if(cutCell.cutEdges[cp] == edgeCode2[face][p]) {
            cutCell.allFacesPointIds[cutCell.noTotalFaces][cutCell.allFacesNoPoints[cutCell.noTotalFaces]] =
                CutCell<nDim>::noCorners + cp;
            cutCell.allFacesNoPoints[cutCell.noTotalFaces]++;
          }
        }
      }
      cutCell.noTotalFaces++;
    }
  }

  // added by Lennart
  // recompute surface and volume centroids:
  // first create triangulation of each surface and determine surface centroid as average triangle centroid
  // then connect all triangles to single point inside the cut cell to form several tetrahedrons which fill up the whole
  // cell volume and compute the volumetric centroid as volume-average of the tetrahedron centroids
  if(improvedCentroidComputation) {
    MFloat vec0[nDim];
    MFloat vec1[nDim];
    MFloat vec2[nDim];

    ASSERT(cutCell.noBoundarySurfaces == 1, "");
    ASSERT(cutCell.allFacesBodyId[0] > -1, "");
    for(MInt bs = 0; bs < cutCell.noBoundarySurfaces; bs++) {
      MFloat bscentroid[nDim] = {F0, F0, F0};
      MFloat bsarea = F0;

      const MInt cp0 = cutCell.allFacesPointIds[bs][0] - CutCell<nDim>::noCorners;
      ASSERT(cutCell.allFacesNoPoints[bs] > 2, "");
      for(MInt v = 2; v < cutCell.allFacesNoPoints[bs]; v++) {
        MInt cp1 = cutCell.allFacesPointIds[bs][v - 1] - CutCell<nDim>::noCorners;
        MInt cp2 = cutCell.allFacesPointIds[bs][v] - CutCell<nDim>::noCorners;
        for(MInt i = 0; i < nDim; i++) {
          vec1[i] = cutCell.cutPoints[cp1][i] - cutCell.cutPoints[cp0][i];
          vec2[i] = cutCell.cutPoints[cp2][i] - cutCell.cutPoints[cp0][i];
        }
        crossProduct(vec0, vec1, vec2);
        MFloat area = F0;
        for(MInt i = 0; i < nDim; i++) {
          area += POW2(vec0[i]);
        }
        area = F1B2 * sqrt(area);
        for(MInt i = 0; i < nDim; i++) {
          bscentroid[i] +=
              area * F1B3 * (cutCell.cutPoints[cp0][i] + cutCell.cutPoints[cp1][i] + cutCell.cutPoints[cp2][i]);
        }
        bsarea += area;
      }

      for(MInt i = 0; i < nDim; i++) {
        bscentroid[i] /= mMax(1e-14, bsarea);
      }
      for(MInt i = 0; i < nDim; i++) {
        cutCell.boundarySurfaceCentroid[bs][i] = bscentroid[i];
      }
    }

    MFloat acentroid[nDim] = {F0, F0, F0};
    MFloat acnt = F0;
    for(MInt cs = 0; cs < cutCell.noCartesianSurfaces; cs++) {
      for(MInt i = 0; i < nDim; i++) {
        acentroid[i] += cutCell.cartFaceCentroid[cs][i];
      }
      acnt += F1;
    }
    for(MInt bs = 0; bs < cutCell.noBoundarySurfaces; bs++) {
      for(MInt i = 0; i < nDim; i++) {
        acentroid[i] += cutCell.boundarySurfaceCentroid[bs][i];
      }
      acnt += F1;
    }
    for(MInt i = 0; i < nDim; i++) {
      acentroid[i] /= acnt;
    }

    MFloat vcentroid[nDim] = {F0, F0, F0};
    MFloat vol2 = F0;
    MFloat vec_tmp[nDim];
    for(MInt as = 0; as < cutCell.noTotalFaces; as++) {
      const MInt p0 = cutCell.allFacesPointIds[as][0];
      if(p0 < CutCell<nDim>::noCorners) {
        for(MInt i = 0; i < nDim; i++) {
          vec0[i] = a_coordinate(cellId, i) + signStencil[p0][i] * cellHalfLength[0] - acentroid[i];
        }
      } else {
        MInt cp0 = p0 - CutCell<nDim>::noCorners;
        for(MInt i = 0; i < nDim; i++) {
          vec0[i] = cutCell.cutPoints[cp0][i] - acentroid[i];
        }
      }
      ASSERT(cutCell.allFacesNoPoints[as] > 2, "");
      for(MInt v = 2; v < cutCell.allFacesNoPoints[as]; v++) {
        MInt p1 = cutCell.allFacesPointIds[as][v - 1];
        MInt p2 = cutCell.allFacesPointIds[as][v];
        if(p1 < CutCell<nDim>::noCorners) {
          for(MInt i = 0; i < nDim; i++) {
            vec1[i] = a_coordinate(cellId, i) + signStencil[p1][i] * cellHalfLength[0] - acentroid[i];
          }
        } else {
          MInt cp1 = p1 - CutCell<nDim>::noCorners;
          for(MInt i = 0; i < nDim; i++) {
            vec1[i] = cutCell.cutPoints[cp1][i] - acentroid[i];
          }
        }
        if(p2 < CutCell<nDim>::noCorners) {
          for(MInt i = 0; i < nDim; i++) {
            vec2[i] = a_coordinate(cellId, i) + signStencil[p2][i] * cellHalfLength[0] - acentroid[i];
          }
        } else {
          MInt cp2 = p2 - CutCell<nDim>::noCorners;
          for(MInt i = 0; i < nDim; i++) {
            vec2[i] = cutCell.cutPoints[cp2][i] - acentroid[i];
          }
        }
        crossProduct(vec_tmp, vec0, vec1);
        MFloat tvol = F0;
        for(MInt i = 0; i < nDim; i++) {
          tvol += F1B6 * vec_tmp[i] * vec2[i];
        }
        vol2 += fabs(tvol);
        for(MInt i = 0; i < nDim; i++) {
          vcentroid[i] +=
              fabs(tvol) * (F3B4 * (acentroid[i] + F1B3 * (vec0[i] + vec1[i] + vec2[i])) + F1B4 * acentroid[i]);
        }
      }
    }
    for(MInt i = 0; i < nDim; i++) {
      vcentroid[i] /= vol2;
      vcentroid[i] -= a_coordinate(cellId, i);
    }

    for(MInt i = 0; i < nDim; i++) {
      cutCell.volumetricCentroid[i] = vcentroid[i];
    }
    // cutCell.volume = vol2;
  }

  return true;
}


//-----------------------------------------------------------------------------------------------------

template <>
void GeometryIntersection<3>::computeCutFaces(std::vector<CutCell<nDim>>& cutCellData, const MInt maxNoSurfaces,
                                              const MInt tCutGroup) {
  using CC = CutCell<nDim>;
  const MUint noCutCells = cutCellData.size();

  // the following variables describe the different corner/edge/face numbering in MAIA and MarchingCubes table
  //   const MInt cornersMCtoSOLVER[ 8 ] = { 2,0,1,3,6,4,5,7 };
  static constexpr MInt cornersSOLVERtoMC[8] = {1, 2, 0, 3, 5, 6, 4, 7};
  static constexpr MInt edgesMCtoSOLVER[12] = {0, 2, 1, 3, 4, 6, 5, 7, 10, 8, 9, 11};
  static constexpr MInt facesMCtoSOLVER[6] = {0, 2, 1, 3, 4, 5};
  static constexpr MFloat signStencil[8][3] = {{-F1, -F1, -F1}, {F1, -F1, -F1}, {-F1, F1, -F1}, {F1, F1, -F1},
                                               {-F1, -F1, F1},  {F1, -F1, F1},  {-F1, F1, F1},  {F1, F1, F1}};
  static constexpr MInt edgeCornerCode[12][2] = {{0, 2}, {1, 3}, {0, 1}, {2, 3}, {4, 6}, {5, 7},
                                                 {4, 5}, {6, 7}, {0, 4}, {1, 5}, {2, 6}, {3, 7}};
  static constexpr MInt faceEdgeCode[6][4] = {{0, 10, 4, 8},  {1, 9, 5, 11}, {2, 8, 6, 9},
                                              {3, 11, 7, 10}, {2, 1, 3, 0},  {6, 4, 7, 5}};
  static constexpr MInt edgeFaceCode[12][2] = {{0, 4}, {1, 4}, {2, 4}, {3, 4}, {0, 5}, {1, 5},
                                               {2, 5}, {3, 5}, {0, 2}, {1, 2}, {0, 3}, {1, 3}};
  static constexpr MInt faceCornerCode[6][4] = {{0, 2, 6, 4}, {3, 1, 5, 7}, {1, 0, 4, 5},
                                                {2, 3, 7, 6}, {0, 1, 3, 2}, {5, 4, 6, 7}};
  const MInt maxNoSets = 6;
  ASSERT(m_noLevelSetsUsedForMb < maxNoSets, "");

  MIntScratchSpace presentCases(m_noLevelSetsUsedForMb, 15, AT_, "presentCases");
  for(MInt s = 0; s < m_noLevelSetsUsedForMb; s++) {
    for(MInt i = 0; i < 15; i++) {
      presentCases(s, i) = 0;
    }
  }

  const MInt maxNoFaces = 100;
  const MInt maxNoVertices = 300;
  const MInt maxNoEdges = 300;

  std::vector<polyVertex> vertices[maxNoSets];
  std::vector<polyFace> faces[maxNoSets];
  std::vector<polyEdge> edges[maxNoSets];
  std::vector<polyCutCell> cutCells;

  stack<MInt> faceStack;

  std::vector<polyMultiFace> multiFaces;
  MIntScratchSpace edgeCutCellPointer(maxNoEdges, AT_, "edgeCutCellPointer");
  std::list<MInt> pureBodyEdges;
  MIntScratchSpace tmp_faces(maxNoFaces, AT_, "tmp_faces");
  MIntScratchSpace multiFaceConnection(maxNoFaces, AT_, "multiFaceConnection");

  MIntScratchSpace outcode_MC_set(m_noLevelSetsUsedForMb, AT_, "outcode_MC_set");
  MIntScratchSpace noCutPointsFromSet(m_noLevelSetsUsedForMb, AT_, "cutSetPointers");
  MIntScratchSpace cutSetPointers(m_noLevelSetsUsedForMb, AT_, "cutSetPointers");

  std::vector<CsgPolygon> csg_polygons;
  std::vector<CsgVertex> csg_vertices;
  std::vector<Csg> csg;
  std::vector<CsgPolygon> result;
  std::vector<polyVertex> vertices_result;
  MIntScratchSpace vertices_renamed(mMax(2, m_noLevelSetsUsedForMb), maxNoVertices, AT_, "vertices_renamed");
  std::list<std::pair<MInt, MInt>> openEdges;

  MIntScratchSpace vertexTouches(maxNoVertices, AT_, "vertexTouches");
  std::vector<MInt> faceVertices;
  MIntScratchSpace bodyFaces(maxNoFaces, AT_, "bodyFaces");
  MFloatScratchSpace normalDotProduct(maxNoFaces, maxNoFaces, AT_,
                                      "normalDotProduct"); // stores 1 - n1*n2 -> in [0;2]
  MIntScratchSpace pointEdgeId(2 * m_noEdges, AT_, "pointEdgeId");
  MIntScratchSpace newCutPointId(maxNoVertices, AT_, "newCutPointId");
  MBoolScratchSpace isPartOfBodySurface(maxNoVertices, AT_, "isPartOfBodySurface");
  MFloatScratchSpace dotProdMatrix(maxNoEdges, maxNoEdges, AT_, "dotProdMatrix");
  //  MIntScratchSpace splitFaceCellIds(m_fvBndryCnd->m_maxNoBndryCells,AT_,"splitFaceCellIds");
  //  splitFaceCellIds.fill(-1);
  //  MIntScratchSpace splitCellIds(noCutCells,AT_,"splitCellIds");
  //  std::vector<cellWithSplitFace> splitFaceCells;
  MIntScratchSpace surfaceIdentificatorCounters(m_noDirs + m_noEmbeddedBodies, AT_, "surfaceIdentificatorCounters");
  MIntScratchSpace cellSurfaceMapping_backup(m_noDirs, AT_, "cellSurfaceMapping_backup");
  MIntScratchSpace splitCellList(CC::maxSplitCells, AT_, "splitCellList");

  NEW_SUB_TIMER_STATIC(tCutFace_0, "cutFaceNew_0", tCutGroup);
  NEW_SUB_TIMER_STATIC(tCutFace_1a, "cutFaceNew_1a", tCutGroup);
  NEW_SUB_TIMER_STATIC(tCutFace_1, "cutFaceNew_1", tCutGroup);
  NEW_SUB_TIMER_STATIC(tCutFace_2, "cutFaceNew_2", tCutGroup);
  NEW_SUB_TIMER_STATIC(tCutFace_3, "cutFaceNew_3", tCutGroup);
  NEW_SUB_TIMER_STATIC(tCutFace_3b, "cutFaceNew_3b", tCutGroup);
  NEW_SUB_TIMER_STATIC(tCutFace_4, "cutFaceNew_4", tCutGroup);
  NEW_SUB_TIMER_STATIC(tCutFace_5a, "cutFaceNew_5a", tCutGroup);
  NEW_SUB_TIMER_STATIC(tCutFace_5b, "cutFaceNew_5b", tCutGroup);
  NEW_SUB_TIMER_STATIC(tCutFace_6, "cutFaceNew_6", tCutGroup);
  NEW_SUB_TIMER_STATIC(tCutFace_7, "cutFaceNew_7", tCutGroup);

#ifdef CutCell_DEBUG
  const MInt debugDomainId = 5;
  const MInt debugTimeStep = 0;
  const MInt debugCellId = 21866;
#endif


  // Use of meaningful bndry-Value for keeping vertices
  // needs to be 10 * m_eps to detect different vertices
  // and needs to be even larger to meaningfully merge edges(angle/calculation)
  // Thus different bounding-values for the different point-Types are used!
  // 10^8*m_eps for point-Types 3 and 10 * m_eps for all other PointTypes!
  const MFloat difBound1 = 10 * m_eps;
  const MFloat difBound2 = 100000000 * m_eps;

  // for( MInt bndryId = m_noOuterBndryCells; bndryId < noBndryCells; bndryId++ ) {
  for(MUint cutc = 0; cutc < noCutCells; cutc++) {
    // MInt bndryId = cutCellData[cutc].cutCellId;
    RECORD_TIMER_START(tCutFace_0);
    // FvBndryCell<nDim>* bndryCell = &m_fvBndryCnd->m_bndryCells->a[bndryId];

    // 0. prepare some cell information
    // const MInt cellId = bndryCell->m_cellId;
    const MInt cellId = cutCellData[cutc].cellId;
    //    const MInt cndId = m_bndryCandidateIds->a[ cellId ];
    //    ASSERT( cndId > -1, "" );
    const MFloat cellLength0 = a_cellLengthAtCell(cellId);
    const MFloat cellHalfLength = F1B2 * cellLength0;

    for(MBool& externalFace : cutCellData[cutc].externalFaces) {
      externalFace = false;
    }

    //    if ( m_noLevelSetsUsedForMb == 1 ) {
    //      RECORD_TIMER_START(tCutFace_1a);
    //      computeCutFaceSimple( cutCellData[cutc] );
    //      RECORD_TIMER_STOP(tCutFace_1a);
    //      continue;
    //    }


    // reset some variables
    for(MInt s = 0; s < m_noLevelSetsUsedForMb; s++) {
      vertices[s].clear();
      faces[s].clear();
      edges[s].clear();
      cutSetPointers[s] = -1;
      noCutPointsFromSet[s] = 0;
    }
    cutCells.clear();
    const MBool isGapCell = cutCellData[cutc].isGapCell;

    MInt startSet = 0;
    MInt endSet = 1;
    if(m_complexBoundary && m_noLevelSetsUsedForMb > 1 && (!isGapCell)) {
      startSet = 1;
      endSet = m_noLevelSetsUsedForMb;
    }

    // preprocess cut points - find out which level set functions contain relevant cut points
    // if no relevant cut points are present -> don't cut the cell with the set
    MInt noCutSets = 0;
    MBool isCompletelyOutside = false;
    MBool hasAmbiguousFaces = false;
    for(MInt s = startSet; s < endSet; s++) {
      // 1.1. get In/Outcode of the corners of the voxel
      // 0 -> Corner is outside Fluid Domain
      // 1 -> Corner is inside Fluid Domain or on Boundary
      unsigned char outcode_MC = 0;
      for(MInt c = 0; c < m_noCorners; c++) {
        // MBool currentOutcode = (!m_pointIsInside[ bndryId ][ IDX_LSSETMB(c, s) ]);
        MBool currentOutcode = (!cutCellData[cutc].cornerIsInsideGeometry[s][c]);
        if(currentOutcode) {
          outcode_MC = outcode_MC | (1 << cornersSOLVERtoMC[c]);
        }
      }
      outcode_MC_set[s] = outcode_MC;
      if(outcode_MC == 0) {
        isCompletelyOutside = true;
      }
    }
    // for( MInt cutPoint = 0; cutPoint < bndryCell->m_srfcs[0]->m_noCutPoints; cutPoint++){
    for(MInt cutPoint = 0; cutPoint < cutCellData[cutc].noCutPoints; cutPoint++) {
      MInt set = 0;
      if(m_complexBoundary && m_noLevelSetsUsedForMb > 1 && (!isGapCell)) {
        // set = m_bodyToSetTable[bndryCell->m_srfcs[0]->m_bodyId[cutPoint]];
        set = m_bodyToSetTable[cutCellData[cutc].cutBodyIds[cutPoint]];
      }

      noCutPointsFromSet[set]++;
    }
    for(MInt set = startSet; set < endSet; set++) {
      if(noCutPointsFromSet[set] && !isCompletelyOutside) {
        cutSetPointers[set] = noCutSets++;
      }
    }
    MInt maxCutsPerFace = 0;
    MInt cutsPerFace[m_noCorners] = {0, 0, 0, 0, 0, 0};
    for(MInt cutPoint = 0; cutPoint < cutCellData[cutc].noCutPoints; cutPoint++) {
      MInt edge = cutCellData[cutc].cutEdges[cutPoint];
      cutsPerFace[edgeFaceCode[edge][0]]++;
      cutsPerFace[edgeFaceCode[edge][1]]++;
    }
    for(MInt k = 0; k < CC::noFaces; k++) {
      maxCutsPerFace = mMax(maxCutsPerFace, cutsPerFace[k]);
    }

    RECORD_TIMER_STOP(tCutFace_0);


    for(MInt set = startSet; set < endSet; set++) {
      MInt setIndex = cutSetPointers[set];
      if(setIndex < 0) {
        continue;
      }
      MInt outcode_MC = outcode_MC_set[set];
      // 1.2. determine Case and check if case is implemented
      MInt currentCase = cases[outcode_MC][0];
      if(!caseStatesLs[currentCase][0]) {
        cerr << "Error: Case not implemented: " << currentCase << endl;
        mTerm(1, AT_, "Case not implemented, see error file for deatails.");
      }
      // 1.2 determine ambiguous cells
      if(!caseStatesLs[currentCase][1]) { // ambiguous case -> disambiguate
        hasAmbiguousFaces = true;
      }
    }

    //    CutCell<nDim> cutcbak = CutCell<nDim>();

    MBool simpleMarchingCubesSucceeds = false;

    if(!hasAmbiguousFaces && noCutSets == 1 && maxCutsPerFace == 2 && !isCompletelyOutside
       && m_noLevelSetsUsedForMb == 1) {
      // if (!hasAmbiguousFaces && noCutSets == 1 && maxCutsPerFace == 2 && !isCompletelyOutside ) {
      RECORD_TIMER_START(tCutFace_1a);
      // dispatch faster MC routine here, store relevant data, and continue
      simpleMarchingCubesSucceeds = computeCutFaceSimple(cutCellData[cutc]);
      RECORD_TIMER_STOP(tCutFace_1a);
    }

    if(simpleMarchingCubesSucceeds) {
      continue;
    }

    MIntScratchSpace noTriangles(m_noLevelSetsUsedForMb, AT_, "noTriangles");
    for(MInt i = 0; i < m_noLevelSetsUsedForMb; i++) {
      noTriangles[i] = 0;
    }


    // computation of cuts with individual sets
    for(MInt set = startSet; set < endSet; set++) { // until line 2552
      MInt setIndex = cutSetPointers[set];
      if(setIndex < 0) {
        continue;
      }

#ifdef CutCell_DEBUG
      if(grid().domainId() == debugDomainId && globalTimeStep == debugTimeStep && cellId == debugCellId) {
        cerr << "Cutting-Cell with " << set << " Index " << setIndex << " bodyId "
             << cutCellData[cutc].associatedBodyIds[set] << endl;
      }
#endif

      RECORD_TIMER_START(tCutFace_1);

      // MInt bodyId = m_associatedBodyIds[ IDX_LSSETMB( cellId, set) ];
      MInt bodyId = cutCellData[cutc].associatedBodyIds[set];
      if(bodyId < 0) continue;

      // store all cut points in a separate array
      // and prepare all vertices for polyeder/polygon datastructure
      // vertices = cut points and fluid corners of the cell
      MInt cutPoints[12] = {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1};
      MInt cutPointToVertexMap[12] = {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1};
      MInt noCutPoints = 0;
      // for( MInt cutPoint = 0; cutPoint < bndryCell->m_srfcs[0]->m_noCutPoints; cutPoint++){
      for(MInt cutPoint = 0; cutPoint < cutCellData[cutc].noCutPoints; cutPoint++) {
        if(m_complexBoundary && m_noLevelSetsUsedForMb > 1 && (!isGapCell)) {
          // if( m_bodyToSetTable[bndryCell->m_srfcs[0]->m_bodyId[ cutPoint ]] !=set )
          if(m_bodyToSetTable[cutCellData[cutc].cutBodyIds[cutPoint]] != set) {
            continue;
          }
        }
        noCutPoints++;
        // cutPoints[ bndryCell->m_srfcs[0]->m_cutEdge[ cutPoint ] ] = cutPoint;
        cutPoints[cutCellData[cutc].cutEdges[cutPoint]] = cutPoint;
        // vertices[setIndex].push_back(polyVertex(bndryCell->m_srfcs[0]->m_cutCoordinates[cutPoint], cutPoint, 1));
        vertices[setIndex].emplace_back(cutCellData[cutc].cutPoints[cutPoint], cutPoint, 1);
        cutPointToVertexMap[cutPoint] = vertices[setIndex].size() - 1;
      }
      MInt cornerToVertexMap[8] = {-1, -1, -1, -1, -1, -1, -1, -1};
      for(MInt c = 0; c < m_noCorners; c++) {
        // if( !m_pointIsInside[ bndryId ][ IDX_LSSETMB(c, set) ]){
        if(!cutCellData[cutc].cornerIsInsideGeometry[set][c]) {
          std::array<MFloat, nDim> tmp_coords{};
          for(MInt i = 0; i < nDim; i++) {
            tmp_coords[i] = a_coordinate(cellId, i) + signStencil[c][i] * cellHalfLength;
          }
          vertices[setIndex].emplace_back(tmp_coords, c, 0);
          cornerToVertexMap[c] = vertices[setIndex].size() - 1;
        }
      }
      // NOTE: the precision of the vertices coordinates is at 3.5 *10^-16

      // 1. find correct marching cubes state and substate -> including disambiguation

      // 1.1. get In/Outcode of the corners of the voxel
      // 0 -> Corner is outside Fluid Domain
      // 1 -> Corner is inside Fluid Domain or on Boundary
      MInt outcode_MC = outcode_MC_set[set];

      // 1.2. determine Case and check if case is implemented
      MInt currentCase = cases[outcode_MC][0];
      presentCases(set, currentCase)++;
      MInt currentSubCase = cases[outcode_MC][1];
      MInt subConfig = 0;
      if(!caseStatesLs[currentCase][0]) {
        cerr << "Error:Case not implemented: " << currentCase << endl;
        mTerm(1, AT_, "Case not implemented, see error file for deatails.");
      }

#ifdef CutCell_DEBUG
      m_caseCheckList[outcode_MC] = 1;
#endif

#ifdef CutCell_DEBUG
      if(grid().domainId() == debugDomainId && globalTimeStep == debugTimeStep && cellId == debugCellId) {
        cerr << "cases: " << outcode_MC << " case " << currentCase << " subcase " << currentSubCase << " ambiguous "
             << caseStatesLs[currentCase][1] << endl;
      }
#endif

      // 1.2 determine ambiguous cells
      if(!caseStatesLs[currentCase][1]) { // ambiguous case -> disambiguate
        // hasAmbiguousFaces = true;
        MInt faceMax = noAmbiguousFaces[currentCase];
        const MInt* testPointer;
        switch(currentCase) {
          case 3:
            testPointer = test3[currentSubCase];
            break;
          case 6:
            testPointer = test6[currentSubCase];
            break;
          case 7:
            testPointer = test7[currentSubCase];
            break;
          case 10:
            testPointer = test10[currentSubCase];
            break;
          case 12:
            testPointer = test12[currentSubCase];
            break;
          case 13:
            testPointer = test13[currentSubCase];
            break;
          default:
            mTerm(1, AT_, "this should not happen!");
            break;
        }
        for(MInt i = 0; i < faceMax; i++) {
          MInt face = testPointer[i]; // MC face
          MBool revert = false;
          if(face < F0) {
            if(face == -6) {
              face = 0;
            } else {
              face = -face;
            }
            revert = true;
          }
          face = facesMCtoSOLVER[face]; // convert to MAIA face
          // MBool testResult = test_face( cndId, face, set );
          MBool testResult = cutCellData[cutc].faceCentroidIsInsideGeometry[set][face];

          if(revert) testResult = !testResult; // invert result for MC processing of body surface

          if(testResult) {
            subConfig += IPOW2(i);
          }
        }
      }

      // 1.3. decide which tiling should be used,
      // how many triangles are to be created for the cut face(s), etc.
      const MInt* tilingPointer = nullptr;
      MBool addVertex = false;
      noTriangles[set] = noTriangles_simpleCases[currentCase];
      switch(currentCase) {
        case 0:
          break;

        case 1:
          tilingPointer = tiling1Ls[currentSubCase];
          break;
        case 2:
          tilingPointer = tiling2Ls[currentSubCase];
          break;
        case 4:
          tilingPointer = tiling4Ls[currentSubCase];
          break;
        case 5:
          tilingPointer = tiling5Ls[currentSubCase];
          break;
        case 8:
          tilingPointer = tiling8Ls[currentSubCase];
          break;
        case 9:
          tilingPointer = tiling9Ls[currentSubCase];
          break;
        case 11:
          tilingPointer = tiling11Ls[currentSubCase];
          break;
        case 14:
          tilingPointer = tiling14Ls[currentSubCase];
          break;

        case 3:
          switch(subConfig) {
            case 0:
              tilingPointer = tiling3_1Ls[currentSubCase];
              noTriangles[set] = 2;
              break;
            case 1:
              tilingPointer = tiling3_2Ls[currentSubCase];
              noTriangles[set] = 4;
              break;
            default:
              mTerm(1, AT_, "This should not happen - we have case 3 with subConfig > 1... exiting");
              break;
          }
          break;

        case 6:
          switch(subConfig) {
            case 0:
              tilingPointer = tiling6_1_1Ls[currentSubCase];
              noTriangles[set] = 3;
              break;
            case 1:
              tilingPointer = tiling6_2Ls[currentSubCase];
              noTriangles[set] = 5;
              break;
            default:
              mTerm(1, AT_, "This should not happen - we have case 6 with subConfig > 1... exiting");
              break;
          }
          break;

        case 7:
          switch(subConfig) {
            case 0:
              tilingPointer = tiling7_1Ls[currentSubCase];
              noTriangles[set] = 3;
              break;
            case 1:
              tilingPointer = tiling7_2Ls[currentSubCase][0];
              noTriangles[set] = 5;
              break;
            case 2:
              tilingPointer = tiling7_2Ls[currentSubCase][1];
              noTriangles[set] = 5;
              break;
            case 3:
              tilingPointer = tiling7_3Ls[currentSubCase][0];
              noTriangles[set] = 9;
              addVertex = true;
              break;
            case 4:
              tilingPointer = tiling7_2Ls[currentSubCase][2];
              noTriangles[set] = 5;
              break;
            case 5:
              tilingPointer = tiling7_3Ls[currentSubCase][1];
              noTriangles[set] = 9;
              addVertex = true;
              break;
            case 6:
              tilingPointer = tiling7_3Ls[currentSubCase][2];
              noTriangles[set] = 9;
              addVertex = true;
              break;
            case 7:
              tilingPointer = tiling7_4_1Ls[currentSubCase];
              noTriangles[set] = 5;
              break;
            default:
              mTerm(1, AT_, "This should not happen - we have case 7 with subConfig > 7... exiting");
              break;
          }
          break;

        case 10:
          switch(subConfig) {
            case 0:
              tilingPointer = tiling10_1_1Ls[currentSubCase];
              noTriangles[set] = 4;
              break;
            case 1:
              tilingPointer = tiling10_2Ls[currentSubCase];
              noTriangles[set] = 8;
              addVertex = true;
              break;
            case 2:
              tilingPointer = tiling10_2_Ls[currentSubCase];
              noTriangles[set] = 8;
              addVertex = true;
              break;
            case 3:
              tilingPointer = tiling10_1_1_Ls[currentSubCase];
              noTriangles[set] = 4;
              break;
            default:
              mTerm(1, AT_, "This should not happen - we have case 10 with subConfig > 3... exiting");
              break;
          }
          break;

        case 12:
          switch(subConfig) {
            case 0:
              tilingPointer = tiling12_1_1Ls[currentSubCase];
              noTriangles[set] = 4;
              break;
            case 1:
              tilingPointer = tiling12_2Ls[currentSubCase];
              noTriangles[set] = 8;
              addVertex = true;
              break;
            case 2:
              tilingPointer = tiling12_2_Ls[currentSubCase];
              noTriangles[set] = 8;
              addVertex = true;
              break;
            case 3:
              tilingPointer = tiling12_1_1_Ls[currentSubCase];
              noTriangles[set] = 4;
              break;
            default:
              mTerm(1, AT_, "This should not happen - we have case 12 with subConfig > 3... exiting");
              break;
          }
          break;

        case 13: {
          MInt config_identificator = 0;
          subConfig = subconfig13[subConfig];
          config_identificator = 0;
          switch(subConfig) {
            case 0: /* 13.1 */
              tilingPointer = tiling13_1Ls[currentSubCase];
              noTriangles[set] = 4;
              break;

            case 1: /* 13.2 */
            case 2: /* 13.2 */
            case 3: /* 13.2 */
            case 4: /* 13.2 */
            case 5: /* 13.2 */
            case 6: /* 13.2 */
              config_identificator = subConfig - 1;
              tilingPointer = tiling13_2Ls[currentSubCase][config_identificator];
              noTriangles[set] = 6;
              break;

            case 7:  /* 13.3 */
            case 8:  /* 13.3 */
            case 9:  /* 13.3 */
            case 10: /* 13.3 */
            case 11: /* 13.3 */
            case 12: /* 13.3 */
            case 13: /* 13.3 */
            case 14: /* 13.3 */
            case 15: /* 13.3 */
            case 16: /* 13.3 */
            case 17: /* 13.3 */
            case 18: /* 13.3 */
              config_identificator = subConfig - 7;
              tilingPointer = tiling13_3Ls[currentSubCase][config_identificator];
              noTriangles[set] = 10;
              addVertex = true;
              break;

            case 19: /* 13.4 */
            case 20: /* 13.4 */
            case 21: /* 13.4 */
            case 22: /* 13.4 */
              config_identificator = subConfig - 19;
              tilingPointer = tiling13_4Ls[currentSubCase][config_identificator];
              noTriangles[set] = 12;
              addVertex = true;
              break;

            case 23: /* 13.5 */
            case 24: /* 13.5 */
            case 25: /* 13.5 */
            case 26: /* 13.5 */
              config_identificator = subConfig - 23;
              tilingPointer = tiling13_5_1Ls[currentSubCase][config_identificator];
              noTriangles[set] = 6;
              break;

            case 27: /* 13.3 */
            case 28: /* 13.3 */
            case 29: /* 13.3 */
            case 30: /* 13.3 */
            case 31: /* 13.3 */
            case 32: /* 13.3 */
            case 33: /* 13.3 */
            case 34: /* 13.3 */
            case 35: /* 13.3 */
            case 36: /* 13.3 */
            case 37: /* 13.3 */
            case 38: /* 13.3 */
              config_identificator = subConfig - 27;
              tilingPointer = tiling13_3_Ls[currentSubCase][config_identificator];
              noTriangles[set] = 10;
              addVertex = true;
              break;

            case 39: /* 13.2 */
            case 40: /* 13.2 */
            case 41: /* 13.2 */
            case 42: /* 13.2 */
            case 43: /* 13.2 */
            case 44: /* 13.2 */
              config_identificator = subConfig - 39;
              tilingPointer = tiling13_2_Ls[currentSubCase][config_identificator];
              noTriangles[set] = 6;
              break;

            case 45: /* 13.1 */
              tilingPointer = tiling13_1_Ls[currentSubCase];
              noTriangles[set] = 4;
              break;

            default:
              mTerm(1, AT_, "impossible MC case 13 subConfig, how could this happen? exiting...");
              break;
          }
        } break;

        default:
          mTerm(1, AT_, "invalid MC case, how could this happen? exiting...");
          break;
      }
      RECORD_TIMER_STOP(tCutFace_1);


      RECORD_TIMER_START(tCutFace_2);
      // 2. build additional vertices, edges, and polygons for body surfaces -> MC
      MInt additionalVertexId = -1;
      if(addVertex) {
        additionalVertexId = vertices[setIndex].size();
        vertices[setIndex].emplace_back(-1, 2);
        addPoint(&vertices[setIndex], cutPointToVertexMap, noCutPoints,
                 vertices[setIndex][additionalVertexId].coordinates);
      }

      ASSERT(noCutPoints >= nDim, "");
      ASSERT(noCutPoints == caseCutPoints[currentCase], to_string(grid().tree().globalId(cellId)));

      for(MInt t = 0; t < noTriangles[set]; t++) {
        // create a new body-face
        // the total body face consinst of noTriangles triangles with each 3 vertices and edges!
        MInt currentFace = faces[setIndex].size();
        faces[setIndex].emplace_back(-1, 1, bodyId);
        // add edges to face; for each edge, check if it already exists!
        MInt p[3];
        for(MInt pt = 0; pt < 3; pt++) {
          MInt cutEdge = tilingPointer[t * 3 + pt];
          if(cutEdge == 12) {
            p[pt] = additionalVertexId;
          } else {
            p[pt] = cutPointToVertexMap[cutPoints[edgesMCtoSOLVER[cutEdge]]];
          }
          ASSERT(p[pt] >= 0, to_string(cutEdge) + " " + to_string(isGapCell) + " " + to_string(cellId) + " "
                                 + to_string(grid().tree().globalId(cellId)) + " " + to_string(edgesMCtoSOLVER[cutEdge])
                                 + " " + to_string(cutPoints[edgesMCtoSOLVER[cutEdge]]) + " "
                                 + to_string(cutPointToVertexMap[cutPoints[edgesMCtoSOLVER[cutEdge]]]) + " "
                                 + to_string(cutCellData[cutc].noCutPoints));
        }
        for(MInt te = 0; te < 3; te++) {
          MInt point = p[te];
          ASSERT(point > -1, "");
          vertices[setIndex][point].surfaceIdentificators.insert(bodyId + m_noDirs);
          MInt nextPoint = p[(te + 1) % 3];
          MBool edgeExists = false;
          for(MInt e = 0; (unsigned)e < vertices[setIndex][p[te]].edges.size(); e++) {
            MInt edge = vertices[setIndex][point].edges[e];
            if(edges[setIndex][edge].vertices[1] == nextPoint) {
              // edge already exists. Link this edge to the face, direction is correct
              faces[setIndex][currentFace].edges.emplace_back(edge, 1);
              const MInt vertexId = edges[setIndex][edge].vertices[0];
              faces[setIndex][currentFace].vertices.push_back(vertexId);
              edges[setIndex][edge].face[0] = currentFace;
              edgeExists = true;
            } else if(edges[setIndex][edge].vertices[0] == nextPoint) {
              // edge already exists. Link this edge to the face, reverse direction
              faces[setIndex][currentFace].edges.emplace_back(edge, -1);
              const MInt vertexId = edges[setIndex][edge].vertices[1];
              faces[setIndex][currentFace].vertices.push_back(vertexId);

              edges[setIndex][edge].face[1] = currentFace;
              edgeExists = true;
            }
          }
          if(!edgeExists) {
            // edge does not exist yet. Create a new edge.
            MInt newEdge = edges[setIndex].size();
            MInt edgeType = 2;
            if(vertices[setIndex][point].pointType == 2 || vertices[setIndex][nextPoint].pointType == 2) {
              edgeType = 3;
            }
            MInt edgeId = -1;
            if(edgeType == 2) {
              MInt p0 = vertices[setIndex][point].pointId;
              MInt p1 = vertices[setIndex][nextPoint].pointId;
              // MInt edge0 = bndryCell->m_srfcs[0]->m_cutEdge[p0];
              MInt edge0 = cutCellData[cutc].cutEdges[p0];
              // MInt edge1 = bndryCell->m_srfcs[0]->m_cutEdge[p1];
              MInt edge1 = cutCellData[cutc].cutEdges[p1];
              if(edgeFaceCode[edge0][0] == edgeFaceCode[edge1][0] || edgeFaceCode[edge0][0] == edgeFaceCode[edge1][1])
                edgeId = edgeFaceCode[edge0][0];
              else if(edgeFaceCode[edge0][1] == edgeFaceCode[edge1][0]
                      || edgeFaceCode[edge0][1] == edgeFaceCode[edge1][1])
                edgeId = edgeFaceCode[edge0][1];
            }
            edges[setIndex].emplace_back(point, nextPoint, edgeId, edgeType);
            vertices[setIndex][point].edges.push_back(newEdge);
            vertices[setIndex][nextPoint].edges.push_back(newEdge);
            faces[setIndex][currentFace].edges.emplace_back(newEdge, 1);
            const MInt vertexId = edges[setIndex][newEdge].vertices[0];
            faces[setIndex][currentFace].vertices.push_back(vertexId);

            edges[setIndex][newEdge].face[0] = currentFace;
          }
        }
        // compute normal:
        computeNormal(vertices[setIndex][p[0]].coordinates, vertices[setIndex][p[1]].coordinates,
                      vertices[setIndex][p[2]].coordinates, faces[setIndex][currentFace].normal,
                      faces[setIndex][currentFace].w);
      }

      for(auto& f : faces[setIndex]) {
        ASSERT(f.vertices.size() == 3, "");
      }

#ifdef CutCell_DEBUG
      if(grid().domainId() == debugDomainId && globalTimeStep == debugTimeStep && cellId == debugCellId) {
        cerr << "Body surfaces: " << endl;
        cerr << "Set: " << set << " noTriangles " << noTriangles[set] << " noFaces " << faces[setIndex].size()
             << " noEdges " << edges[setIndex].size() << " noVertices " << vertices[setIndex].size() << endl;
      }
#endif


      RECORD_TIMER_STOP(tCutFace_2);

      RECORD_TIMER_START(tCutFace_3);

      // 3. build polygons for Cartesian faces
      // -> all edges should already be built, just compose the faces!
      // -> other than the body-faces cartesian faces are poligons, so they may be triangles or other
      //    multiple-edge-faces!
      // 3.1. prepare basic edges -> between two inside corner vertices or between a corner vertex and a cut point
      for(MInt e = 0; e < m_noEdges; e++) {
        // MBool p0Fluid = !m_pointIsInside[ bndryId ][ IDX_LSSETMB(edgeCornerCode[e][0], set) ];
        MBool p0Fluid = !cutCellData[cutc].cornerIsInsideGeometry[set][edgeCornerCode[e][0]];
        // MBool p1Fluid = !m_pointIsInside[ bndryId ][ IDX_LSSETMB(edgeCornerCode[e][1], set) ];
        MBool p1Fluid = !cutCellData[cutc].cornerIsInsideGeometry[set][edgeCornerCode[e][1]];
        if(p0Fluid && p1Fluid) { // created a full Cartesian edge
          MInt v0 = cornerToVertexMap[edgeCornerCode[e][0]];
          MInt v1 = cornerToVertexMap[edgeCornerCode[e][1]];
          edges[setIndex].emplace_back(v0, v1, e, 0);
          vertices[setIndex][v0].edges.push_back(edges[setIndex].size() - 1);
          vertices[setIndex][v1].edges.push_back(edges[setIndex].size() - 1);
        } else if(p0Fluid) { // create a cut edge from p0 to cut point
          MInt v0 = cornerToVertexMap[edgeCornerCode[e][0]];
          MInt v1 = cutPointToVertexMap[cutPoints[e]];
          edges[setIndex].emplace_back(v0, v1, e, 1);
          vertices[setIndex][v0].edges.push_back(edges[setIndex].size() - 1);
          vertices[setIndex][v1].edges.push_back(edges[setIndex].size() - 1);
        } else if(p1Fluid) { // create a cut edge from p1 to cut point
          MInt v0 = cornerToVertexMap[edgeCornerCode[e][1]];
          MInt v1 = cutPointToVertexMap[cutPoints[e]];
          edges[setIndex].emplace_back(v0, v1, e, 1);
          vertices[setIndex][v0].edges.push_back(edges[setIndex].size() - 1);
          vertices[setIndex][v1].edges.push_back(edges[setIndex].size() - 1);
        } // else create no edge since edge is fully located outside
      }

      // 3.2. compose faces
      //      add cartesian faces!
      for(MInt cartFace = 0; cartFace < m_noDirs; cartFace++) {
        // find starting point: cut point on first edge of face
        MBool edgeTouched[4] = {false, false, false, false};
        for(MInt e = 0; e < 4; e++) {
          if(edgeTouched[e]) continue;
          MInt cartEdge = faceEdgeCode[cartFace][e];
          MInt startVertex = cornerToVertexMap[faceCornerCode[cartFace][e]]; // first point of the corresponding edge
                                                                             // in math pos. sense of rotation
          if(startVertex == -1) continue;

          // add a new face
          MInt currentFace = faces[setIndex].size();
          faces[setIndex].emplace_back(cartFace, 0, -1);
          for(MInt i = 0; i < nDim; i++) {
            faces[setIndex][currentFace].normal[i] = F0;
          }
          MInt normalDir = cartFace / 2;
          MFloat normalSign = -F1;
          if(cartFace % 2 == 0) {
            normalSign = F1;
          }
          faces[setIndex][currentFace].w = -vertices[setIndex][startVertex].coordinates[normalDir] * normalSign;
          faces[setIndex][currentFace].normal[normalDir] = normalSign;

          // do loop around the edges until the starting vertex is reached again
          MInt currentE = e;
          MInt currentVertex = startVertex;
          do {
            vertices[setIndex][currentVertex].surfaceIdentificators.insert(cartFace);
            edgeTouched[currentE] = true;
            // check all edges of the currentVertex; if the edge is found which is located on the investigated
            // cartesianEdge, continue along this edge!
            MBool edgeFound = false;
            for(MInt ve = 0; (unsigned)ve < vertices[setIndex][currentVertex].edges.size(); ve++) {
              MInt edge = vertices[setIndex][currentVertex].edges[ve];
              if(edges[setIndex][edge].edgeType > 1) // only a Cartesian line (cut or uncut) can be followed
                continue;
              if(edges[setIndex][edge].edgeId != cartEdge) continue;
              // otherwise we found cut or uncut Cartesian edge to follow!
              edgeFound = true;
              MInt edgeDirection = 1;
              if(edges[setIndex][edge].vertices[1] == currentVertex) {
                edgeDirection = -1; // reverse edge
              }
              faces[setIndex][currentFace].edges.emplace_back(edge, edgeDirection);
              MInt vertexId = edges[setIndex][edge].vertices[0];
              if(edgeDirection == -1) vertexId = edges[setIndex][edge].vertices[1];
              faces[setIndex][currentFace].vertices.push_back(vertexId);

              if(edges[setIndex][edge].face[0] > -1)
                edges[setIndex][edge].face[1] = currentFace;
              else
                edges[setIndex][edge].face[0] = currentFace;
              if(edgeDirection == 1)
                currentVertex = edges[setIndex][edge].vertices[1];
              else
                currentVertex = edges[setIndex][edge].vertices[0];
              break;
            }
            ASSERT(edgeFound, "");
            // check if new vertex is a cut point. if no, just proceed to the next edge; otherwise, continue along the
            // corresponding cut line before returning to loop;
            if(vertices[setIndex][currentVertex].pointType == 0) {
              currentE = (currentE + 1) % 4;
              cartEdge = faceEdgeCode[cartFace][currentE];
            } else if(vertices[setIndex][currentVertex].pointType == 1) {
              edgeFound = false;
              for(MInt ve = 0; (unsigned)ve < vertices[setIndex][currentVertex].edges.size(); ve++) {
                MInt edge = vertices[setIndex][currentVertex].edges[ve];
                if(edges[setIndex][edge].edgeType != 2) // only a cut line on a Cartesian face can be followed
                  continue;
                if(edges[setIndex][edge].edgeId != cartFace) continue;
                // otherwise we found cut line on cartFace to follow
                edgeFound = true;
                MInt edgeDirection = 1;
                if(edges[setIndex][edge].vertices[1] == currentVertex) {
                  edgeDirection = -1; // reverse edge
                }
                faces[setIndex][currentFace].edges.emplace_back(edge, edgeDirection);
                MInt vertexId = edges[setIndex][edge].vertices[0];
                if(edgeDirection == -1) vertexId = edges[setIndex][edge].vertices[1];
                faces[setIndex][currentFace].vertices.push_back(vertexId);

                edges[setIndex][edge].face[1] = currentFace;
                if(edgeDirection == 1)
                  currentVertex = edges[setIndex][edge].vertices[1];
                else
                  currentVertex = edges[setIndex][edge].vertices[0];
                // cartEdge = bndryCell->m_srfcs[0]->m_cutEdge[vertices[setIndex][currentVertex].pointId];
                cartEdge = cutCellData[cutc].cutEdges[vertices[setIndex][currentVertex].pointId];
                while(faceEdgeCode[cartFace][currentE] != cartEdge) {
                  currentE = (currentE + 1) % 4; // adjust the relative edge counter
                }
                break;
              }
              ASSERT(edgeFound, "");
            }
          } while(currentVertex != startVertex); // end while

        } // end for e
      }   // end for cartFace

#ifdef CutCell_DEBUG
      if(grid().domainId() == debugDomainId && globalTimeStep == debugTimeStep && cellId == debugCellId) {
        cerr << "Including cartesian faces: " << endl;
        cerr << "Set: " << set << " SetIndex:" << setIndex << " faces : " << faces[setIndex].size()
             << " edges : " << edges[setIndex].size() << " certices : " << vertices[setIndex].size() << endl;
        for(MUint i = 0; i < vertices[setIndex].size(); i++) {
          cerr << " vertices " << setprecision(19) << vertices[setIndex][i].coordinates[0] << " "
               << vertices[setIndex][i].coordinates[1] << " " << vertices[setIndex][i].coordinates[2] << endl;
        }
      }
#endif


      RECORD_TIMER_STOP(tCutFace_3);
    }

    RECORD_TIMER_START(tCutFace_3b);

    MInt noIndividualCuts = 0;
    for(MInt set = startSet; set < endSet; set++) {
      if(cutSetPointers[set] > -1) {
        noIndividualCuts++;
      }
    }

#ifdef CutCell_DEBUG
    if(grid().domainId() == debugDomainId && globalTimeStep == debugTimeStep && cellId == debugCellId) {
      cerr << "noIndividualCuts: " << noIndividualCuts << endl;
    }
#endif

    if(noIndividualCuts == 0) {
      cutCellData[cutc].volume = F0;
      RECORD_TIMER_STOP(tCutFace_3b);
      continue;
    }
    ASSERT(noIndividualCuts > 0, "");
    ASSERT(noCutSets > 0, " ");
    ASSERT(noCutSets == noIndividualCuts, "");

    // 4.a combine polyhedra that are computed with different sets to one polyhedron
    // set up the CSG datastructure for each polyhedron
    MBool error = false;
    const MInt startSetIndex = 0;
    MInt referenceSet = startSet;
    for(MInt set = startSet; set < endSet; set++) {
      if(cutSetPointers[set] > -1) {
        referenceSet = set;
        break;
      }
    }

    // NOTE: at this point we have multiple bodySurfaces/faces/triangles which each have 3 vertices/edges!
    //      and multiple cartisain-faces which are poligons
    //      (the poligons can be described through a normal in one of the coordinate diretions
    //       and a coordinate value in that direction)


    // check that face-edges and face-vertices match!
    for(MInt set = startSet; set < endSet; set++) {
      MInt setIndex = cutSetPointers[set];
      if(setIndex < 0) continue;

      for(MInt faceId = 0; faceId < (signed)faces[setIndex].size(); faceId++) {
        ASSERT(faces[setIndex][faceId].edges.size() == faces[setIndex][faceId].vertices.size(), "");
        for(MInt edgeId = 0; edgeId < (signed)faces[setIndex][faceId].edges.size(); edgeId++) {
          const MInt edge = faces[setIndex][faceId].edges[edgeId].first;
          const MInt newVertexId = faces[setIndex][faceId].vertices[edgeId];
          const MInt direction = faces[setIndex][faceId].edges[edgeId].second;
          MInt vertexId = edges[setIndex][edge].vertices[0];
          if(direction == -1) vertexId = edges[setIndex][edge].vertices[1];
          ASSERT(newVertexId == vertexId, "");
        }
      }
    }

    // NOTE: the bodySurfaces/faces (triangles) can not be joined into poligons at this point yet,
    //      as the current cgs-libary can either handle triangles or
    //      poligons with normals into one of the coordinate-directions!

    MInt noInitialFaces = 0;

    //====== /MC ======
    //====== /MC ======
    //====== /MC ======

    if(noIndividualCuts > 1) { // until line 3011

      ASSERT(m_complexBoundary && m_noLevelSetsUsedForMb > 1, "");
      ASSERT(m_bodyFaceJoinMode == 4 || m_bodyFaceJoinMode == 1, "ERROR: using unrecommended bodyFaceJoinMode");

      // referenceSet is the first Set with cutPoints
      // -> loop though all other sets!
      for(MInt set = referenceSet + 1; set < endSet; set++) {
        MInt setIndex = cutSetPointers[set];
        if(setIndex < 0) continue;

        csg.clear();
        csg_polygons.clear();


        // now adding the faces for the zero SetIndex (the first set with cutPoints)
        for(MInt face = faces[startSetIndex].size() - 1; face > -1; face--) {
          csg_vertices.clear();
          polyFace* faceP = &faces[startSetIndex][face];

          for(MInt e = 0; (unsigned)e < faceP->vertices.size(); e++) {
            const MInt vertexId = faceP->vertices[e];
            csg_vertices.emplace_back(CsgVector(vertices[startSetIndex][vertexId].coordinates), vertexId,
                                      startSetIndex);
          }
          if(faceP->faceType == 0) {
            CsgVector normal(faceP->normal);
            CsgPlane plane(normal, -faceP->w);
            csg_polygons.emplace_back(csg_vertices, startSetIndex, faceP->faceId, faceP->faceType, faceP->bodyId,
                                      plane);
          } else {
            csg_polygons.emplace_back(csg_vertices, startSetIndex, faceP->faceId, faceP->faceType, faceP->bodyId);
          }
        }

        csg.emplace_back(csg_polygons);
        csg_polygons.clear();

        // now adding the faces for current setIndex (second/third.. set with cutPoints)
        for(MInt face = faces[setIndex].size() - 1; face > -1; face--) {
          csg_vertices.clear();
          polyFace* faceP = &faces[setIndex][face];

          for(MInt e = 0; (unsigned)e < faceP->vertices.size(); e++) {
            const MInt vertexId = faceP->vertices[e];
            csg_vertices.emplace_back(CsgVector(vertices[setIndex][vertexId].coordinates), vertexId, setIndex);
          }
          if(faceP->faceType == 0) {
            CsgVector normal(faceP->normal);
            CsgPlane plane(normal, -faceP->w);
            csg_polygons.emplace_back(csg_vertices, setIndex, faceP->faceId, faceP->faceType, faceP->bodyId, plane);
          } else {
            csg_polygons.emplace_back(csg_vertices, setIndex, faceP->faceId, faceP->faceType, faceP->bodyId);
          }
        }
        csg.emplace_back(csg_polygons);

#ifdef CutCell_DEBUG
        if(grid().domainId() == debugDomainId && globalTimeStep == debugTimeStep && cellId == debugCellId) {
          cerr << "Combine: reference " << referenceSet << " with " << set << "Index " << startSetIndex << "face-size "
               << faces[startSetIndex].size() << "Index " << setIndex << "face-size " << faces[setIndex].size() << endl;
        }
#endif

        // call intersect
        // currently alle sets are intersected with the zero-set!
        // this is not meaning-ful for triple-set-intersections!
        result.clear();
        result = csg[0].intersect(csg[1]);
        vertices_result.clear();

// the result (a vector of poligons) holds all information for all intersected poligons!
#ifdef CutCell_DEBUG
        if(grid().domainId() == debugDomainId && globalTimeStep == debugTimeStep && cellId == debugCellId) {
          cerr << "result-size: " << result.size() << "vertices_result " << vertices_result.size() << endl;

          for(MInt poligon = 0; (unsigned)poligon < result.size(); poligon++) {
            cerr << "resultId " << poligon << " body " << result[poligon].bodyId << " type " << result[poligon].faceType
                 << " setIndex " << setIndex << endl;
            for(MInt vertex = 0; (unsigned)vertex < result[poligon].vertices.size(); vertex++) {
              cerr << "vertex " << vertex << "coordinate: " << setprecision(21)
                   << result[poligon].vertices[vertex].pos.xx[0] << " " << result[poligon].vertices[vertex].pos.xx[1]
                   << " " << result[poligon].vertices[vertex].pos.xx[2] << " "
                   << result[poligon].vertices[vertex].vertexId << " " << result[poligon].vertices[vertex].setIndex
                   << endl;
            }
            cerr << "faceId " << result[poligon].faceId << endl;
          }
        }
#endif

        // -------- postprocess the poligon intersection (vertices, edges & faces)  --------------

        // --------1)  regenerate polyVertices in a new vertices_result vector
        //    which holds all unique vertices and removes double-entries from the result-class!


        // maximal amount of vertices for the two boundarySurfaces
        for(MInt i = 0; i < mMax(2, m_noLevelSetsUsedForMb); i++) {
          for(MInt j = 0; j < maxNoVertices; j++) {
            vertices_renamed(i, j) = -1;
          }
        }

        // count of unique vertices
        MInt noVertices = 0;

        // a) add vertices which were already part of ther original vertices-vector and thus hold
        //   additional polyVertex information!

        // loop over all resulting poligons
        for(MInt p = 0; (unsigned)p < result.size(); p++) {
          ASSERT(result[p].vertices.size() <= (unsigned)maxNoVertices, "");
          ASSERT((signed)result[p].vertices.size() >= 3, "");

          // loop over all vertices of each poligon
          for(MInt v = 0; (unsigned)v < result[p].vertices.size(); v++) {
            MInt vertexId = result[p].vertices[v].vertexId;
            MInt vertexSetIndex = result[p].vertices[v].setIndex;
            if(vertexSetIndex < 0) continue;

            if(vertexId < 0) continue;
            //-> vertex corresponds to an existing vertex
            //   and has not been added to the result vector yet!

            if(vertices_renamed(vertexSetIndex, vertexId) == -1) {
              // -> vertex has not been added to vertices_result vector yet!
              MBool vertexFound = false;
              MInt vertexIndex = -1;

              // loop over all vertices already added to the result-vector
              for(MInt i = 0; (unsigned)i < vertices_result.size(); i++) {
                MFloat coord_diff = F0;
                coord_diff += pow(result[p].vertices[v].pos.xx[0] - vertices_result[i].coordinates[0], 2);
                coord_diff += pow(result[p].vertices[v].pos.xx[1] - vertices_result[i].coordinates[1], 2);
                coord_diff += pow(result[p].vertices[v].pos.xx[2] - vertices_result[i].coordinates[2], 2);

                if(m_multiCutCell) coord_diff = pow(coord_diff, 0.5);


                if(coord_diff < difBound1) {
#ifdef CutCell_DEBUG
                  if(grid().domainId() == debugDomainId && globalTimeStep == debugTimeStep && cellId == debugCellId) {
                    cerr << "Small diff " << p << " " << vertexId << " " << vertices_result[i].coordinates[0] << " "
                         << vertices_result[i].coordinates[1] << " " << vertices_result[i].coordinates[2] << " eps "
                         << difBound1 << " dif " << coord_diff << " " << result[p].vertices[v].pos.xx[0] << " "
                         << result[p].vertices[v].pos.xx[1] << " " << result[p].vertices[v].pos.xx[2] << " " << endl;
                  }
#endif

                  vertexFound = true;
                  vertexIndex = i;
                  break;
                }
              }
              if(vertexFound) {
                vertices_renamed(vertexSetIndex, vertexId) = vertexIndex;
              } else {
                // only add a vertex if its not already part of other poligon-faces!
                vertices_renamed(vertexSetIndex, vertexId) = noVertices;
                vertices_result.emplace_back(vertices[vertexSetIndex][vertexId].pointId,
                                             vertices[vertexSetIndex][vertexId].pointType);
                // only adding corner-, cut- ord clipping points at this time!
                ASSERT(vertices[vertexSetIndex][vertexId].pointType == 0
                           || vertices[vertexSetIndex][vertexId].pointType == 1
                           || vertices[vertexSetIndex][vertexId].pointType == 2,
                       "");
                vertices_result[noVertices].coordinates[0] = vertices[vertexSetIndex][vertexId].coordinates[0];
                vertices_result[noVertices].coordinates[1] = vertices[vertexSetIndex][vertexId].coordinates[1];
                vertices_result[noVertices].coordinates[2] = vertices[vertexSetIndex][vertexId].coordinates[2];
                noVertices++;
              }
            }
            // reset vertexId and set setIndex to be jumped
            result[p].vertices[v].vertexId = vertices_renamed(vertexSetIndex, vertexId);
            result[p].vertices[v].setIndex = -1;
          }
        }

        // b) add new/additional vertices that are created due to the multi-cutCell intersection
        //   are only existing in the result-class

        for(MInt p = 0; (unsigned)p < result.size(); p++) {
          for(MInt v = 0; (unsigned)v < result[p].vertices.size(); v++) {
            MInt vertexId = result[p].vertices[v].vertexId;
            if(vertexId == -1) { // remaining vertices...
              // check if vertex has been added to vertices_result vector
              MBool vertexFound = false;
              MInt vertexIndex = -1;
              for(MInt i = 0; (unsigned)i < vertices_result.size(); i++) {
                MFloat coord_diff = F0;
                coord_diff += pow(result[p].vertices[v].pos.xx[0] - vertices_result[i].coordinates[0], 2);
                coord_diff += pow(result[p].vertices[v].pos.xx[1] - vertices_result[i].coordinates[1], 2);
                coord_diff += pow(result[p].vertices[v].pos.xx[2] - vertices_result[i].coordinates[2], 2);

                if(m_multiCutCell) coord_diff = pow(coord_diff, 0.5);

                // use different bound between mc vertices!
                MFloat difBound = difBound1;
                if(m_multiCutCell) {
                  difBound = difBound2;
                }

                if(coord_diff < difBound) {
#ifdef CutCell_DEBUG
                  if(grid().domainId() == debugDomainId && globalTimeStep == debugTimeStep && cellId == debugCellId) {
                    cerr << "2 Small diff " << p << " " << vertexId << " " << vertices_result[i].coordinates[0] << " "
                         << vertices_result[i].coordinates[1] << " " << vertices_result[i].coordinates[2] << " eps "
                         << coord_diff << " " << result[p].vertices[v].pos.xx[0] << " "
                         << result[p].vertices[v].pos.xx[1] << " " << result[p].vertices[v].pos.xx[2] << " " << endl;
                  }
#endif

                  vertexFound = true;
                  vertexIndex = i;
                  break;
                }
              }
              if(vertexFound) {
                result[p].vertices[v].vertexId = vertexIndex;
                result[p].vertices[v].setIndex = -1;
              } else {
                // all points with pointId -1 and PointType 3!
                // thus additional points from multiCutCell generation! (MC vertex!)
                vertices_result.emplace_back(-1, 3);
                vertices_result[noVertices].coordinates[0] = result[p].vertices[v].pos.xx[0];
                vertices_result[noVertices].coordinates[1] = result[p].vertices[v].pos.xx[1];
                vertices_result[noVertices].coordinates[2] = result[p].vertices[v].pos.xx[2];
                result[p].vertices[v].vertexId = noVertices;
                result[p].vertices[v].setIndex = -1;
                noVertices++;
              }
            }
          }
        }

        // copy vertices-result into vertices
        vertices[startSetIndex].swap(vertices_result);
        vertices_result.clear();
        edges[startSetIndex].clear();
        faces[startSetIndex].clear();

#ifdef CutCell_DEBUG
        if(grid().domainId() == debugDomainId && globalTimeStep == debugTimeStep && cellId == debugCellId) {
          cerr << "VertexId-check: " << endl;

          for(MInt poligon = 0; (unsigned)poligon < result.size(); poligon++) {
            cerr << "resultId " << poligon << " body " << result[poligon].bodyId << " type " << result[poligon].faceType
                 << " setIndex " << setIndex << endl;
            for(MInt vertex = 0; (unsigned)vertex < result[poligon].vertices.size(); vertex++) {
              cerr << "vertex " << vertex << "coordinate: " << setprecision(21)
                   << result[poligon].vertices[vertex].pos.xx[0] << " " << result[poligon].vertices[vertex].pos.xx[1]
                   << " " << result[poligon].vertices[vertex].pos.xx[2] << " "
                   << result[poligon].vertices[vertex].vertexId << " " << result[poligon].vertices[vertex].setIndex
                   << endl;
            }
            cerr << "faceId " << result[poligon].faceId << endl;
          }

          for(MInt poligon = 0; (unsigned)poligon < result.size(); poligon++) {
            for(MInt vertex = 0; (unsigned)vertex < result[poligon].vertices.size(); vertex++) {
              MInt vertexId = result[poligon].vertices[vertex].vertexId;
              for(MInt i = 0; i < nDim; i++) {
                if(fabs(result[poligon].vertices[vertex].pos.xx[i] - vertices[startSetIndex][vertexId].coordinates[i])
                   > m_eps * 10) {
                  cerr << " VertexId-missmatch " << cellId << " " << vertexId << " " << vertex << " " << i
                       << setprecision(21) << result[poligon].vertices[vertex].pos.xx[i] << " "
                       << vertices[startSetIndex][vertexId].coordinates[i] << endl;
                }
              }
            }
          }
        }

#endif


        // --------1) regenerate unique edges and faces (already unique)
        //            edge regeneration is rather complicated, as the poligons hold no edge
        //            information! Just their verticies!
        MInt noFaces = 0;
        noInitialFaces = 0;
        MIntScratchSpace faceMapping(result.size(), AT_, "faceMapping");

        // a) add all poligon faces with more than 3 valid vertices to the resulting faces-vector!
        for(MInt p = 0; (unsigned)p < result.size(); p++) {
          MBoolScratchSpace vertexValid(result[p].vertices.size(), AT_, "vertexValid");
          MBool polygonValid = true;
          MInt validVertices = 0;
          faceMapping[p] = -1;

          for(MInt v = 0; (unsigned)v < result[p].vertices.size(); v++) {
            const MInt j = (v + 1) % result[p].vertices.size();
            const MInt vertexId = result[p].vertices[v].vertexId;
            const MInt vertexIdNext = result[p].vertices[j].vertexId;
            if(vertexId == vertexIdNext) {
              vertexValid[v] = false;
            } else {
              vertexValid[v] = true;
              validVertices++;
            }
            // add surface surfaceIdentificators
            if(vertexValid[v]) {
              MInt surfaceIndicator = -1;
              if(result[p].faceType == 0) {
                surfaceIndicator = result[p].faceId;
              } else {
                surfaceIndicator = result[p].bodyId + m_noDirs;
              }
              vertices[startSetIndex][vertexId].surfaceIdentificators.insert(surfaceIndicator);
            }
          }
          if(validVertices < 3) {
#ifdef CutCell_DEBUG
            // uncomment the part below for additional debug output
            // can be helpful for debugging
            /*
            cerr << grid().domainId() << " Invalid poligon-Id " << p << "with less than 3 verticies! "
                 << " in cell " << grid().tree().globalId(cellId) << " at TS " << globalTimeStep
                 << " result vertices: " << endl;
            for(MInt v = 0; (unsigned)v < result[p].vertices.size(); v++) {
              cerr << result[p].vertices[v].vertexId;
            }
            cerr << endl;
            */
#endif
            polygonValid = false;
          }

          if(!polygonValid) continue;

          faces[startSetIndex].emplace_back(result[p].faceId, result[p].faceType, result[p].bodyId);
          faces[startSetIndex][noFaces].normal[0] = result[p].plane.normal.xx[0];
          faces[startSetIndex][noFaces].normal[1] = result[p].plane.normal.xx[1];
          faces[startSetIndex][noFaces].normal[2] = result[p].plane.normal.xx[2];
          ASSERT(
              !(std::isnan(result[p].plane.normal.xx[0] + result[p].plane.normal.xx[1] + result[p].plane.normal.xx[2])),
              "");
          faces[startSetIndex][noFaces].tmpSetIndex = result[p].setIndex;
          faces[startSetIndex][noFaces].w = -result[p].plane.w;
          faceMapping[p] = noFaces;
          faces[startSetIndex][noFaces].edges.clear();
          faces[startSetIndex][noFaces].vertices.clear();

          // new method to find face-Vertices!
          for(MInt v = 0; (unsigned)v < result[p].vertices.size(); v++) {
            const MInt vertexId = result[p].vertices[v].vertexId;
            MBool vertexUsed = false;
            for(MInt vertice : faces[startSetIndex][noFaces].vertices) {
              if(vertexId == vertice) {
                vertexUsed = true;
                break;
              }
            }

            if(!vertexUsed) {
              faces[startSetIndex][noFaces].vertices.push_back(vertexId);
              vertices[startSetIndex][vertexId].faceIds.push_back(noFaces);
            }
          }
          ASSERT((signed)faces[startSetIndex][noFaces].vertices.size() >= 3, "");
          noFaces++;
          noInitialFaces++;
        }

#ifdef CutCell_DEBUG
        if(grid().domainId() == debugDomainId && globalTimeStep == debugTimeStep && cellId == debugCellId) {
          cerr << "Current-Face-Count " << noFaces << endl;
        }
#endif

        // b) add edges to all valid faces
        //   in the original version vertices are also added at this point!
        MInt noEdges = 0;
#if defined CutCell_DEBUG || !defined NDEBUG
        MBool onlySingleEdges = true;
#endif
        for(MInt p = 0; (unsigned)p < result.size(); p++) {
          MInt faceId = faceMapping[p];
          if(faceId < 0) continue;

          // add edges to the face
          for(MInt v = 0; (unsigned)v < result[p].vertices.size(); v++) {
            MInt j = (v + 1) % result[p].vertices.size();
            MInt vertexId = result[p].vertices[v].vertexId;
            MInt vertexIdNext = result[p].vertices[j].vertexId;
            if(vertexId == vertexIdNext) {
              continue;
            }

            MInt direction = 1;
            MBool edgeFound = false;
            MInt edgeId = -1;
            for(MInt e = 0; (unsigned)e < edges[startSetIndex].size(); e++) {
              MInt v0 = edges[startSetIndex][e].vertices[0];
              MInt v1 = edges[startSetIndex][e].vertices[1];
              if(vertexId == v0 && vertexIdNext == v1) {
                edgeFound = true;
                edgeId = e;
                direction = 1;
                break;
              }
              if(vertexId == v1 && vertexIdNext == v0) {
                edgeFound = true;
                edgeId = e;
                direction = -1;
                break;
              }
            }
            // adding a new edge
            if(!edgeFound) {
              edges[startSetIndex].emplace_back(vertexId, vertexIdNext, -1, -1);
              edges[startSetIndex][noEdges].face[0] = faceId;
              edgeId = noEdges;
              direction = 1;
              vertices[startSetIndex][vertexId].edges.push_back(noEdges);
              vertices[startSetIndex][vertexIdNext].edges.push_back(noEdges);
              noEdges++;

#ifdef CutCell_DEBUG
              if(grid().domainId() == debugDomainId && globalTimeStep == debugTimeStep && cellId == debugCellId) {
                cerr << "Adding Edge " << noEdges << " " << vertexId << " " << vertexIdNext << " " << faceId << endl;
              }
#endif

              // edge already exists and has only one other face
            } else if(edges[startSetIndex][edgeId].face[1] < 0) {
              edges[startSetIndex][edgeId].face[1] = faceId;

#ifdef CutCell_DEBUG
              if(grid().domainId() == debugDomainId && globalTimeStep == debugTimeStep && cellId == debugCellId) {
                cerr << "Adding 2nd face to egde " << edgeId << " " << vertexId << " " << vertexIdNext << " " << faceId
                     << endl;
                cerr << "edge vertices are: " << edges[startSetIndex][edgeId].vertices[0] << " "
                     << edges[startSetIndex][edgeId].vertices[1] << endl;
              }
#endif
            } else {
              // edge exists and already has two valid faces!
              //=> adding new edge which is identical to the first edge
#if defined CutCell_DEBUG || !defined NDEBUG
              onlySingleEdges = false;
#endif
              edges[startSetIndex].emplace_back(vertexId, vertexIdNext, -1, -1);
              edges[startSetIndex][noEdges].face[0] = faceId;
              edgeId = noEdges;
              direction = 1;
              vertices[startSetIndex][vertexId].edges.push_back(noEdges);
              vertices[startSetIndex][vertexIdNext].edges.push_back(noEdges);
              noEdges++;
#ifdef CutCell_DEBUG
              cerr << cellId << "Has a duplicate edge!" << endl;
              if(grid().domainId() == debugDomainId && globalTimeStep == debugTimeStep && cellId == debugCellId) {
                cerr << "Adding Duplicate Edge " << noEdges << " " << vertexId << " " << vertexIdNext << " " << faceId
                     << endl;
              }
#endif
            }
            faces[startSetIndex][faceId].edges.emplace_back(edgeId, direction);
          }
        }

        // debug check:
        for(MInt faceId = 0; faceId < static_cast<MInt>(faces[startSetIndex].size()); faceId++) {
          ASSERT(faces[startSetIndex][faceId].edges.size() == faces[startSetIndex][faceId].vertices.size(), "");
          for(MInt edgeId = 0; edgeId < static_cast<MInt>(faces[startSetIndex][faceId].edges.size()); edgeId++) {
            const MInt edge = faces[startSetIndex][faceId].edges[edgeId].first;
            const MInt newVertexId = faces[startSetIndex][faceId].vertices[edgeId];
            const MInt direction = faces[startSetIndex][faceId].edges[edgeId].second;
            MInt vertexId = edges[startSetIndex][edge].vertices[0];
            if(direction == -1) {
              vertexId = edges[startSetIndex][edge].vertices[1];
            }
            ASSERT(newVertexId == vertexId, "");
          }
        }

        // openEdges: are edges which only belong to only one face!

        // 1) create openEdge vectors/information
        //   openEdges    : vector with edgeId and edge-Direction
        //   openEdgeList : vector with edgeId and edge-Direction (identical)
        //   openEdgeId   : collector which is -1 for closedEdges and has the openEdgeId for openEdges
        openEdges.clear();

        for(MInt e = 0; e < static_cast<MInt>(edges[startSetIndex].size()); e++) {
          ASSERT(edges[startSetIndex][e].face[0] > -1, "");
          if(edges[startSetIndex][e].face[1] == -1) {
            MInt face = edges[startSetIndex][e].face[0];
            MInt dir = 0;
            for(MInt fe = 0; fe < (signed)faces[startSetIndex][face].edges.size(); fe++) {
              if(faces[startSetIndex][face].edges[fe].first == e) {
                dir = faces[startSetIndex][face].edges[fe].second;
              }
            }
            ASSERT(dir, "");
            if(dir == 1) {
              openEdges.emplace_back(e, -1);
            } else {
              openEdges.emplace_back(e, 1);
            }
          }
        }

        vector<std::pair<MInt, MInt>> openEdgeList;
        MIntScratchSpace openEdgeId(maxNoEdges, AT_, "openEdgeId");
        openEdgeId.fill(-1);
        for(auto& openEdge : openEdges) {
          openEdgeId(openEdge.first) = static_cast<MInt>(openEdgeList.size());
          openEdgeList.push_back(openEdge);
        }

#ifdef CutCell_DEBUG
        if(grid().domainId() == debugDomainId && globalTimeStep == debugTimeStep && cellId == debugCellId) {
          cerr << openEdges.size() << " initially open edges! " << endl;
          cerr << "Count: " << faces[startSetIndex].size() << " faces " << edges[startSetIndex].size() << " edges "
               << vertices[startSetIndex].size() << " vertices " << endl;

          for(auto it2 = openEdges.begin(); it2 != openEdges.end(); it2++) {
            const MInt edgdeId = (*it2).first;
            cerr << " openEdge with Id " << (*it2).first << " v1 " << edges[startSetIndex][edgdeId].vertices[0] << " "
                 << edges[startSetIndex][edgdeId].vertices[1] << endl;
          }
        }
#endif

        //
        // NOTE: connect single large edge with two coinciding smaller edges and intermediate point,
        //      which is indeed meaningful and necessary!!!


        // a) debug-check
        //   no same edges(with the same vertices, but different direction) should occour:
        //   no already deleted/invalid edges should occour!
        // NOTE: edges are considered deleted if their face[0] = -1!
#if defined CutCell_DEBUG || !defined NDEBUG
        for(MInt i = 0; i < (signed)openEdgeList.size(); i++) {
          const MInt e = openEdgeList[i].first;
          const MInt dir = openEdgeList[i].second;
          const MInt id0 = (dir == 1) ? 0 : 1;
          const MInt id1 = (dir == 1) ? 1 : 0;
          const MInt v0 = edges[startSetIndex][e].vertices[id0];
          const MInt v1 = edges[startSetIndex][e].vertices[id1];

          ASSERT(edges[startSetIndex][e].face[0] > -1, "");

          // loop over all edged which share the first vertex
          for(MInt j = 0; j < (signed)vertices[startSetIndex][v0].edges.size(); j++) {
            const MInt e2 = vertices[startSetIndex][v0].edges[j];

            if(e2 == e) continue;

            const MInt face = edges[startSetIndex][e2].face[0];
            MInt dir2 = 0;
            for(MInt fe = 0; fe < (signed)faces[startSetIndex][face].edges.size(); fe++) {
              if(faces[startSetIndex][face].edges[fe].first == e2) {
                dir2 = faces[startSetIndex][face].edges[fe].second;
              }
            }
            const MInt id00 = (dir2 == 1) ? 0 : 1;
            const MInt id11 = (dir2 == 1) ? 1 : 0;
            const MInt v00 = edges[startSetIndex][e2].vertices[id00];
            const MInt v11 = edges[startSetIndex][e2].vertices[id11];

            ASSERT(edges[startSetIndex][e2].face[0] > -1, "");
            if(onlySingleEdges) {
              ASSERT(v0 != v00 || v1 != v11, "");
              ASSERT(v0 != v11 || v1 != v00, "");
            }
          }
        }
#endif


        // necessary while loop to join edges
        //---------------------------------------------------------------------------------------
        MBool somethingChanged = !openEdgeList.empty();
        while(somethingChanged) {
          somethingChanged = false;

          MBool alreadyResolved = false;

          for(MUint k1 = 0; k1 < openEdgeList.size(); k1++) {
            // gather infomation for openEdge e1!
            const MInt e1 = openEdgeList[k1].first;
            if(openEdgeId(e1) < 0) continue;
            const MInt dir1 = openEdgeList[k1].second;
            const MInt rdir1 = (dir1 == 1) ? -1 : 1;
            const MInt id0 = (dir1 == 1) ? 0 : 1;
            const MInt id1 = (dir1 == 1) ? 1 : 0;
            const MInt v0 = edges[startSetIndex][e1].vertices[id0];
            const MInt v1 = edges[startSetIndex][e1].vertices[id1];
            const MInt face = edges[startSetIndex][e1].face[0];

            // determine the faceEdgeId fe
            MInt fe = -1;
            for(fe = 0; (unsigned)fe < faces[startSetIndex][face].edges.size(); fe++) {
              if(faces[startSetIndex][face].edges[fe].first == e1) {
                break;
              }
            }
            ASSERT(fe < (signed)faces[startSetIndex][face].edges.size(), "");

            // determine insert-position for new vertex:
            // the new vertex is added between v0 and v1 (thus the +1)
            MInt fv = -1;
            for(fv = 0; (unsigned)fv < faces[startSetIndex][face].vertices.size(); fv++) {
              if(faces[startSetIndex][face].vertices[fv] == v0 || faces[startSetIndex][face].vertices[fv] == v1) {
                if(fv + 1 < (signed)faces[startSetIndex][face].vertices.size()
                   && (faces[startSetIndex][face].vertices[fv + 1] == v0
                       || faces[startSetIndex][face].vertices[fv + 1] == v1)) {
                  break;
                }
              }
            }
            fv = fv + 1;

            MFloat a[3] = {F0, F0, F0};
            MFloat a_abs = F0;
            for(MInt i = 0; i < nDim; i++) {
              a[i] = vertices[startSetIndex][v1].coordinates[i] - vertices[startSetIndex][v0].coordinates[i];
              a_abs += a[i] * a[i];
            }
            a_abs = sqrt(a_abs);

            // loop over all edges at vertex v0
            for(MInt k2 = 0; (unsigned)k2 < vertices[startSetIndex][v0].edges.size(); k2++) {
              // find other openEdges at v0
              const MInt e2 = vertices[startSetIndex][v0].edges[k2];
              if(e2 == e1) continue;
              if(openEdgeId(e2) < 0) continue;
              ASSERT(openEdgeList[openEdgeId(e2)].first == e2, "");
              // gather information for the also openEdge e2!
              const MInt dir2 = openEdgeList[openEdgeId(e2)].second;
              const MInt id00 = (dir2 == 1) ? 0 : 1;
              const MInt id11 = (dir2 == 1) ? 1 : 0;
              const MInt v00 = edges[startSetIndex][e2].vertices[id00];
              const MInt v11 = edges[startSetIndex][e2].vertices[id11];

              if(v11 != v0) continue;

              MFloat b[3]{};
              MFloat b_abs = F0;
              MFloat dotProduct = F0;
              for(MInt i = 0; i < nDim; i++) {
                b[i] = vertices[startSetIndex][v11].coordinates[i] - vertices[startSetIndex][v00].coordinates[i];
                b_abs += b[i] * b[i];
                dotProduct += a[i] * b[i];
              }
              b_abs = sqrt(b_abs);
              // dotProduct should be -1 if intermediate vertex on edge was found
              dotProduct = fabs(F1 + (dotProduct / (a_abs * b_abs)));
              // TODO labels:GEOM,totest use meaningful restriction, which also changes testcases!
              // if ( dotProduct < m_eps*1000000 && b_abs < a_abs ) {
              if(dotProduct < 1e-10 && b_abs < a_abs) {
                // splitting e1 into e2 and e3:
                // v0 == v11 version!
                // version 1

                // check if the third edge e3, which shall replace e1, already exists
                // NOTE: if the edge is triple-split e3 must not already be existing!
                //      a new edge is then created, which is split again in the next iteration!
                MInt e3 = -1;
                MInt dir3 = -2;
                MInt v000 = -1;
                MInt v111 = -1;
                MInt id000 = -1;
                MInt id111 = -1;
                MBool foundMatchingEdge = false;
                for(MUint k3 = 0; k3 < openEdgeList.size(); k3++) {
                  e3 = openEdgeList[k3].first;
                  if(openEdgeId(e3) < 0) continue;
                  dir3 = openEdgeList[k3].second;
                  id000 = (dir3 == 1) ? 0 : 1;
                  id111 = (dir3 == 1) ? 1 : 0;
                  v000 = edges[startSetIndex][e3].vertices[id000];
                  v111 = edges[startSetIndex][e3].vertices[id111];
                  if(v111 == v00 && v000 == v1) {
                    foundMatchingEdge = true;
                    break;
                  }
                  if(v000 == v00 && v111 == v1) {
                    foundMatchingEdge = true;
                    break;
                  }
                }
                // check that the matchingEdge is not already a closed edge
                // if it is, the edges must not be split!
                // instead the open edges will be resolved by adding faceLines!
                /*
                if(!foundMatchingEdge) {
                  MBool alreadyExisting = false;
                  for(MInt e = 0; e < static_cast<MInt>(edges[startSetIndex].size()); e++) {
                    MInt vl = edges[startSetIndex][e].vertices[0];
                    MInt vr = edges[startSetIndex][e].vertices[1];

                    if((v00 == vl && v1 == vr) || (v00 == vr && v1 == vl)) {
                      alreadyExisting = true;
                      break;
                    }
                  }
                  if(alreadyExisting) {
                    continue;
                  }
                }
                */
                // make e1 equal e3
                edges[startSetIndex][e1].vertices[id0] = v00;

                // adding face to e2
                edges[startSetIndex][e2].face[1] = face;

                // adding e2 to the face of e1
                MInt otherDir = (edges[startSetIndex][e2].vertices[id0] == v0) ? rdir1 : dir1;
                auto pos = faces[startSetIndex][face].edges.begin() + fe;
                pos += (dir1 == 1) ? 0 : 1;
                faces[startSetIndex][face].edges.insert(pos, make_pair(e2, otherDir));


                if(!m_multiCutCell) {
                  // replacing e1 at the vertex v0
                  for(MInt edg = 0; (unsigned)edg < vertices[startSetIndex][v0].edges.size(); edg++) {
                    if(vertices[startSetIndex][v0].edges[edg] == e1) {
                      vertices[startSetIndex][v0].edges[edg] = e2;
                    }
                  }
                } else {
                  // removing e1 at the vertex v0
                  for(MInt edg = 0; (unsigned)edg < vertices[startSetIndex][v0].edges.size(); edg++) {
                    if(vertices[startSetIndex][v0].edges[edg] == e1) {
                      auto epos = vertices[startSetIndex][v0].edges.begin() + edg;
                      vertices[startSetIndex][v0].edges.erase(epos);
                    }
                  }
                }
                // add e1 to v00
                if(m_multiCutCell) vertices[startSetIndex][v00].edges.push_back(e1);


                if(foundMatchingEdge) {
                  MInt otherDir2 = (edges[startSetIndex][e2].vertices[id0] == v0) ? rdir1 : dir1;

                  // find new position of e1 in the face (might be shifted due to the insertion of e2)!
                  for(fe = 0; fe < (signed)faces[startSetIndex][face].edges.size(); fe++) {
                    if(faces[startSetIndex][face].edges[fe].first == e1) {
                      faces[startSetIndex][face].edges[fe].first = e3;
                      faces[startSetIndex][face].edges[fe].second = otherDir2;
                    }
                  }

                  // add the face to the edge e3
                  edges[startSetIndex][e3].face[1] = face;

                  // deleting edge e1
                  edges[startSetIndex][e1].face[0] = -1;
                  edges[startSetIndex][e1].face[1] = -1;

#ifdef CutCell_DEBUG
                  if(grid().domainId() == debugDomainId && globalTimeStep == debugTimeStep && cellId == debugCellId) {
                    cerr << "Replacing edge " << e1 << " with " << e3 << " and adding vertex " << v00 << " to face "
                         << face << endl;
                    cerr << "Intermediate edge is " << e2 << " with v " << v1 << " v2 " << v0 << endl;
                  }
#endif

                  if(!m_multiCutCell) {
                    // replacing e1 with e3
                    for(MInt edg = 0; (unsigned)edg < vertices[startSetIndex][v00].edges.size(); edg++) {
                      if(vertices[startSetIndex][v00].edges[edg] == e1) {
                        vertices[startSetIndex][v00].edges[edg] = e3;
                      }
                    }
                    for(MInt edg = 0; (unsigned)edg < vertices[startSetIndex][v1].edges.size(); edg++) {
                      if(vertices[startSetIndex][v1].edges[edg] == e1) {
                        vertices[startSetIndex][v1].edges[edg] = e3;
                      }
                    }
                  } else {
                    // removing e1 at the vertices v00 and v1
                    for(MInt edg = 0; (unsigned)edg < vertices[startSetIndex][v00].edges.size(); edg++) {
                      if(vertices[startSetIndex][v00].edges[edg] == e1) {
                        auto epos = vertices[startSetIndex][v00].edges.begin() + edg;
                        vertices[startSetIndex][v00].edges.erase(epos);
                      }
                    }
                    for(MInt edg = 0; (unsigned)edg < vertices[startSetIndex][v1].edges.size(); edg++) {
                      if(vertices[startSetIndex][v1].edges[edg] == e1) {
                        auto epos = vertices[startSetIndex][v1].edges.begin() + edg;
                        vertices[startSetIndex][v1].edges.erase(epos);
                      }
                    }
                  }

                  // remove both edges from the openEdgeList
                  openEdgeId(e1) = -1;
                  openEdgeId(e3) = -1;
                }

                // adding the additional vertex v00 to the face
                auto vpos = faces[startSetIndex][face].vertices.begin() + fv;
                if(vpos > faces[startSetIndex][face].vertices.end()) {
                  faces[startSetIndex][face].vertices.push_back(v00);
                } else {
                  faces[startSetIndex][face].vertices.insert(vpos, v00);
                }

                // remove as openEdge
                openEdgeId(e2) = -1;

                somethingChanged = true;
                alreadyResolved = true;
                break;
              }
            }

            if(alreadyResolved) continue;

            // loop over all edges at vertex v1
            for(MInt k2 = 0; (unsigned)k2 < vertices[startSetIndex][v1].edges.size(); k2++) {
              const MInt e2 = vertices[startSetIndex][v1].edges[k2];
              if(e2 == e1) continue;
              if(openEdgeId(e2) < 0) continue;
              ASSERT(openEdgeList[openEdgeId(e2)].first == e2, "");
              const MInt dir2 = openEdgeList[openEdgeId(e2)].second;
              const MInt id00 = (dir2 == 1) ? 0 : 1;
              const MInt id11 = (dir2 == 1) ? 1 : 0;
              const MInt v00 = edges[startSetIndex][e2].vertices[id00];
              const MInt v11 = edges[startSetIndex][e2].vertices[id11];

              if(v00 != v1) continue;

              MFloat b[3] = {F0, F0, F0};
              MFloat b_abs = F0;
              MFloat dotProduct = F0;
              for(MInt i = 0; i < nDim; i++) {
                b[i] = vertices[startSetIndex][v11].coordinates[i] - vertices[startSetIndex][v00].coordinates[i];
                b_abs += b[i] * b[i];
                dotProduct += a[i] * b[i];
              }
              b_abs = sqrt(b_abs);
              // dotProduct should be -1 if intermediate vertex on edge was found
              dotProduct = fabs(F1 + (dotProduct / (a_abs * b_abs)));
              // TODO labels:GEOM,totest use meaningful restriction, which also changes testcases!
              // if ( dotProduct < m_eps*1000000 && b_abs < a_abs ) {
              if(dotProduct < 1e-10 && b_abs < a_abs) {
                // splitting e1 into e2 and e3:
                // v1 == v00 version!
                // version 2

                // check if the third edge e3, which shall replace e1, already exists
                // NOTE: if the edge is triple-split e3 must not already be existing!
                //      a new edge is then created, which is split again in the next iteration!
                MInt e3 = -1;
                MInt dir3 = -2;
                MInt v000 = -1;
                MInt v111 = -1;
                MInt id000 = -1;
                MInt id111 = -1;
                MBool foundMatchingEdge = false;
                for(auto& k3 : openEdgeList) {
                  e3 = k3.first;
                  if(openEdgeId(e3) < 0) continue;
                  dir3 = k3.second;
                  id000 = (dir3 == 1) ? 0 : 1;
                  id111 = (dir3 == 1) ? 1 : 0;
                  v000 = edges[startSetIndex][e3].vertices[id000];
                  v111 = edges[startSetIndex][e3].vertices[id111];
                  if(v000 == v0 && v111 == v11) {
                    foundMatchingEdge = true;
                    break;
                  } else if(v000 == v11 && v111 == v0) {
                    foundMatchingEdge = true;
                    break;
                  }
                }
                // check that the matchingEdge is not already a closed edge
                /*
                if(!foundMatchingEdge) {
                  MBool alreadyExisting = false;
                  for(MInt e = 0; e < static_cast<MInt>(edges[startSetIndex].size()); e++) {
                    MInt vl = edges[startSetIndex][e].vertices[0];
                    MInt vr = edges[startSetIndex][e].vertices[1];

                    if((v11 == vl && v1 == vr) || (v11 == vr && v1 == vl)) {
                      alreadyExisting = true;
                      break;
                    }
                  }
                  if(alreadyExisting) {
                    continue;
                  }
                }
                */

                // adding e2 to the face of e1
                MInt otherDir = (edges[startSetIndex][e2].vertices[id1] == v1) ? rdir1 : dir1;
                auto pos = faces[startSetIndex][face].edges.begin() + fe;
                pos += (dir1 == 1) ? 1 : 0;
                faces[startSetIndex][face].edges.insert(pos, make_pair(e2, otherDir));

                // adding the additional vertex v11 to the face
                auto vpos = faces[startSetIndex][face].vertices.begin() + fv;
                if(vpos > faces[startSetIndex][face].vertices.end()) {
                  faces[startSetIndex][face].vertices.push_back(v11);
                } else {
                  faces[startSetIndex][face].vertices.insert(vpos, v11);
                }

                // adding face to e2
                edges[startSetIndex][e2].face[1] = face;

                // replacing e1 at the vertex v1 with e2
                if(!m_multiCutCell) {
                  for(MInt edg = 0; (unsigned)edg < vertices[startSetIndex][v1].edges.size(); edg++) {
                    if(vertices[startSetIndex][v1].edges[edg] == e1) {
                      vertices[startSetIndex][v1].edges[edg] = e2;
                    }
                  }
                } else {
                  // removing e1 at the vertex v1
                  for(MInt edg = 0; (unsigned)edg < vertices[startSetIndex][v1].edges.size(); edg++) {
                    if(vertices[startSetIndex][v1].edges[edg] == e1) {
                      auto epos = vertices[startSetIndex][v1].edges.begin() + edg;
                      vertices[startSetIndex][v1].edges.erase(epos);
                    }
                  }
                }
                // adding e1 to vertex v11
                if(m_multiCutCell) vertices[startSetIndex][v11].edges.push_back(e1);

                // make e1 euqal e3
                edges[startSetIndex][e1].vertices[id1] = v11;

                // remove as openEdge
                openEdgeId(e2) = -1;


                if(foundMatchingEdge) {
                  // removing edge e1:

                  // replacing e1 with e3 for the face
                  const MInt otherDir2 = (edges[startSetIndex][e3].vertices[id1] == v1) ? rdir1 : dir1;

                  // find new position of e1 in the face (might be shifted due to the insertion of e2)!
                  for(fe = 0; fe < (signed)faces[startSetIndex][face].edges.size(); fe++) {
                    if(faces[startSetIndex][face].edges[fe].first == e1) {
                      faces[startSetIndex][face].edges[fe].first = e3;
                      faces[startSetIndex][face].edges[fe].second = otherDir2;
                      break;
                    }
                  }

#ifdef CutCell_DEBUG
                  if(grid().domainId() == debugDomainId && globalTimeStep == debugTimeStep && cellId == debugCellId) {
                    cerr << "Replacing edge " << e1 << " with " << e3 << " and adding vertex " << v00 << " to face "
                         << face << endl;
                    cerr << "Intermediate edge is " << e2 << " with v " << v11 << " v2 " << v1 << endl;
                  }
#endif

                  // delete e1/remove all faces:
                  edges[startSetIndex][e1].face[0] = -1;
                  edges[startSetIndex][e1].face[1] = -1;

                  // replacing e1 with e3
                  if(!m_multiCutCell) {
                    for(MInt edg = 0; (unsigned)edg < vertices[startSetIndex][v0].edges.size(); edg++) {
                      if(vertices[startSetIndex][v0].edges[edg] == e1) {
                        vertices[startSetIndex][v0].edges[edg] = e3;
                      }
                    }
                    for(MInt edg = 0; (unsigned)edg < vertices[startSetIndex][v11].edges.size(); edg++) {
                      if(vertices[startSetIndex][v11].edges[edg] == e1) {
                        vertices[startSetIndex][v11].edges[edg] = e3;
                      }
                    }
                  } else {
                    // removing e1 from all vertices
                    for(MInt edg = 0; (unsigned)edg < vertices[startSetIndex][v0].edges.size(); edg++) {
                      if(vertices[startSetIndex][v0].edges[edg] == e1) {
                        auto epos = vertices[startSetIndex][v0].edges.begin() + edg;
                        vertices[startSetIndex][v0].edges.erase(epos);
                      }
                    }
                    for(MInt edg = 0; (unsigned)edg < vertices[startSetIndex][v11].edges.size(); edg++) {
                      if(vertices[startSetIndex][v11].edges[edg] == e1) {
                        auto epos = vertices[startSetIndex][v11].edges.begin() + edg;
                        vertices[startSetIndex][v11].edges.erase(epos);
                      }
                    }
                  }
                  // remove openEdges
                  openEdgeId(e1) = -1;
                  openEdgeId(e3) = -1;

                  // add the face to e3
                  edges[startSetIndex][e3].face[1] = face;
                }

                somethingChanged = true;
                break;
              }
            }
          }
        }
//----------------------------------------------------------------------------------------

//  debug-checks:
//  no same edges(with the same vertices, but different direction) should occour,
//  unless they are marked as deleted!
#ifdef CutCell_DEBUG
        for(MUint i = 0; i < openEdgeList.size(); i++) {
          const MInt e = openEdgeList[i].first;
          const MInt dir = openEdgeList[i].second;
          const MInt id0 = (dir == 1) ? 0 : 1;
          const MInt id1 = (dir == 1) ? 1 : 0;
          const MInt v0 = edges[startSetIndex][e].vertices[id0];
          const MInt v1 = edges[startSetIndex][e].vertices[id1];

          // skip-deleted edges:
          if(edges[startSetIndex][e].face[0] == -1 && edges[startSetIndex][e].face[1] == -1) continue;

          // loop over all edged which share the first vertex
          for(MInt j = 0; (unsigned)j < vertices[startSetIndex][v0].edges.size(); j++) {
            const MInt e2 = vertices[startSetIndex][v0].edges[j];
            if(e2 == e) continue;
            const MInt dir2 = openEdgeList[openEdgeId(e2)].second;
            const MInt id00 = (dir2 == 1) ? 0 : 1;
            const MInt id11 = (dir2 == 1) ? 1 : 0;
            const MInt v00 = edges[startSetIndex][e2].vertices[id00];
            const MInt v11 = edges[startSetIndex][e2].vertices[id11];

            // skip-deleted edges:
            if(edges[startSetIndex][e2].face[0] == -1 && edges[startSetIndex][e2].face[1] == -1) continue;
            if(onlySingleEdges) {
              ASSERT(v0 != v00 || v1 != v11, "");
              ASSERT(v0 != v11 || v1 != v00, "");
            }
          }
        }
#endif

        // no face should have a reference to an deleted edge:
        for(MInt faceId = 0; faceId < (signed)faces[startSetIndex].size(); faceId++) {
          for(MInt edgeId = 0; edgeId < (signed)faces[startSetIndex][faceId].edges.size(); edgeId++) {
            const MInt edge = faces[startSetIndex][faceId].edges[edgeId].first;
            ASSERT(edges[startSetIndex][edge].face[0] > -1, "");
          }
        }

        // no vertex should have a reference to an deleted edge:
        for(MInt v = 0; v < (signed)vertices[startSetIndex].size(); v++) {
          for(MInt edgeId = 0; (unsigned)edgeId < vertices[startSetIndex][v].edges.size(); edgeId++) {
            const MInt edge = vertices[startSetIndex][v].edges[edgeId];
            ASSERT(edges[startSetIndex][edge].face[0] > -1, "");
          }
        }

        // vertex and edge information should match!
        for(MInt faceId = 0; faceId < (signed)faces[startSetIndex].size(); faceId++) {
          ASSERT(faces[startSetIndex][faceId].edges.size() == faces[startSetIndex][faceId].vertices.size(), "");
          for(MInt edgeId = 0; edgeId < (signed)faces[startSetIndex][faceId].edges.size(); edgeId++) {
            const MInt edge = faces[startSetIndex][faceId].edges[edgeId].first;
            const MInt newVertexId = faces[startSetIndex][faceId].vertices[edgeId];
            const MInt direction = faces[startSetIndex][faceId].edges[edgeId].second;
            MInt vertexId = edges[startSetIndex][edge].vertices[0];
            if(direction == -1) vertexId = edges[startSetIndex][edge].vertices[1];
            ASSERT(newVertexId == vertexId, "");
          }
        }

        // determine remaining unsolved open edges
        std::list<std::pair<MInt, MInt>> openEdgesOld(openEdges);
        openEdges.clear();
        for(auto& it1 : openEdgesOld) {
          if(openEdgeId(it1.first) < 0) continue;
          openEdges.push_back(it1);
        }


#ifdef CutCell_DEBUG
        if(grid().domainId() == debugDomainId && globalTimeStep == debugTimeStep && cellId == debugCellId) {
          cerr << openEdges.size() << " open edge remaining: " << endl;
          cerr << "Count: " << faces[startSetIndex].size() << " faces " << edges[startSetIndex].size() << " edges "
               << vertices[startSetIndex].size() << " vertices " << endl;

          for(auto it2 = openEdges.begin(); it2 != openEdges.end(); it2++) {
            const MInt edgeId = (*it2).first;
            cerr << " openEdge with Id " << edgeId << " v1 " << edges[startSetIndex][edgeId].vertices[0] << " "
                 << edges[startSetIndex][edgeId].vertices[1] << endl;
          }
        }

        /*
        if(openEdges.size() > 0) {
          cerr << openEdges.size() << " openEdges remaining: adding face-Lines! at TS " << globalTimeStep << " for "
               << cellId << " on rank " << grid().domainId() << " result " << result.size() << " faces "
               << faces[startSetIndex].size() << endl;
        }
        */

        if(openEdges.size() == 1) {
          cerr << " Only 1 open edge remaining: This means trouble! " << endl;
          cerr << "Re-checking coordinate-distances! " << endl;
          for(MInt i = 0; (unsigned)i < vertices[startSetIndex].size(); i++) {
            for(MInt j = 0; (unsigned)j < vertices[startSetIndex].size(); j++) {
              if(i == j) continue;
              MFloat coord_diff = 0;
              coord_diff +=
                  pow(vertices[startSetIndex][i].coordinates[0] - vertices[startSetIndex][j].coordinates[0], 2);
              coord_diff +=
                  pow(vertices[startSetIndex][i].coordinates[1] - vertices[startSetIndex][j].coordinates[1], 2);
              coord_diff +=
                  pow(vertices[startSetIndex][i].coordinates[2] - vertices[startSetIndex][j].coordinates[2], 2);

              MFloat dif_old = coord_diff;
              coord_diff = pow(coord_diff, 0.5);

              if(coord_diff < cellLength0) {
                cerr << " Coord-dif " << coord_diff << " i " << i << " " << j << endl;
                cerr << " dif-old " << dif_old << " " << m_eps * 10 << endl;
              }
            }
          }
        }

#endif

        // NOTE: these open-edges are due to the deletion of MC-verticies close to each other
        //      this was done to reduce the overall number of vertices, limited by maxNoFaceVertices

        // prepare decision array ->
        // compute the angle between all open edges!
        // and always connect edges with angle closest to 180 degrees!
        for(auto it1 = openEdges.begin(); it1 != openEdges.end(); it1++) {
          // first edge e1
          const MInt e1 = (*it1).first;
          const MInt dir1 = (*it1).second;
          MInt v11 = edges[startSetIndex][e1].vertices[0];
          MInt v12 = edges[startSetIndex][e1].vertices[1];
          if(dir1 == -1) {
            MInt tmp = v11;
            v11 = v12;
            v12 = tmp;
          }
          MFloat a[3]{};
          MFloat a_abs = F0;
          for(MInt i = 0; i < nDim; i++) {
            a[i] = vertices[startSetIndex][v12].coordinates[i] - vertices[startSetIndex][v11].coordinates[i];
            a_abs += a[i] * a[i];
          }
          a_abs = sqrt(a_abs);
          // loop ever all remaining openEdges
          for(auto& openEdge : openEdges) {
            const MInt e2 = openEdge.first;
            const MInt dir2 = openEdge.second;
            // identical edges, set to maximum!
            if(e2 == e1) {
              dotProdMatrix(e1, e2) = 1000.0;
              dotProdMatrix(e2, e1) = 1000.0;
              continue;
            }
            MInt v21 = edges[startSetIndex][e2].vertices[0];
            MInt v22 = edges[startSetIndex][e2].vertices[1];
            if(dir2 == -1) {
              MInt tmp = v21;
              v21 = v22;
              v22 = tmp;
            }

            if(v12 == v21) {
              // edges are connected correctly, compute the dotProduct
              MFloat b[3] = {F0, F0, F0};
              MFloat b_abs = F0;
              MFloat dotProduct = F0;
              for(MInt i = 0; i < nDim; i++) {
                b[i] = vertices[startSetIndex][v22].coordinates[i] - vertices[startSetIndex][v21].coordinates[i];
                b_abs += b[i] * b[i];
                dotProduct += a[i] * b[i];
              }
              b_abs = sqrt(b_abs);
              dotProduct = abs(F1 - abs(dotProduct / (a_abs * b_abs)));
              dotProdMatrix(e1, e2) = dotProduct;
            } else {
              // edges with different orientation should never be matched!
              dotProdMatrix(e1, e2) = 1000.0;
            }
          }
        }

        // loop over all remaining edges
        // adding face-lines for these!
        while(!openEdges.empty()) {
          faceVertices.clear();
          // openEdge eLast
          MInt e = -1;
          MInt eLast = openEdges.front().first;
          MInt dir = openEdges.front().second;
          MInt faceType = -1;
          MInt faceId = -1;
          MInt bodyId = -1;
          for(MInt i = 0; i < m_noDirs + m_noEmbeddedBodies; i++) {
            surfaceIdentificatorCounters[i] = 0;
          }

          // adding a new face for the remaining unsolved edge (eLast)
          faces[startSetIndex].emplace_back(faceId, faceType, bodyId);
          // add the edge to the face
          faces[startSetIndex][noFaces].edges.emplace_back(eLast, dir);
          faces[startSetIndex][noFaces].isLine = 1;

          // remove corresponding openEdge
          openEdges.pop_front();

          MInt startVertex = edges[startSetIndex][eLast].vertices[0];
          MInt vertex = edges[startSetIndex][eLast].vertices[1];
          // switch vertex-order
          if(dir == -1) {
            MInt tmp = startVertex;
            startVertex = vertex;
            vertex = tmp;
          }
          // add vertices to the faceLine
          faces[startSetIndex][noFaces].vertices.push_back(startVertex);

#ifdef CutCell_DEBUG
          if(grid().domainId() == debugDomainId && globalTimeStep == debugTimeStep && cellId == debugCellId) {
            cerr << "Adding face line starting with edge " << eLast << " with vertex " << startVertex << " and "
                 << vertex << endl;
          }
#endif

          // add faceVertexId
          faceVertices.push_back(vertex);
          // increase all surfaceIdentificator of the vertex!
          for(MInt surfaceIdentificator : vertices[startSetIndex][vertex].surfaceIdentificators) {
            surfaceIdentificatorCounters[surfaceIdentificator]++;
          }
          // add the new faceLine to the edge (=> no longer open)
          edges[startSetIndex][eLast].face[1] = noFaces;

          // loop over vertices to find a closed faceLine including other edges/vertices
          while(vertex != startVertex) {
            // first: find the right open edge for further connection, then connect it:
            auto bestEdgeIt = openEdges.end();
            MFloat minDotProduct = 1001.0;
            // loop over all other remaining open edges to find the bestEdgeId
            // meaning the edge with the lowest dotProduct!
            // NOTE: in some cases the edge with the lowest dotProduct is not the
            // best choice, as the edged might not be connected to the previous edge!
            for(auto it = openEdges.begin(); it != openEdges.end(); it++) {
              // other openEdge Id e
              e = (*it).first;
              if(dotProdMatrix(eLast, e) < minDotProduct) {
                // only use this edge if the edge continues the loop!
                dir = (*bestEdgeIt).second;
                MInt vStart = edges[startSetIndex][e].vertices[0];
                MInt vEnd = edges[startSetIndex][e].vertices[1];
                if(dir == -1) {
                  MInt tmp = vStart;
                  vStart = vEnd;
                  vEnd = tmp;
                }
                if(vStart == vertex) {
                  minDotProduct = dotProdMatrix(eLast, e);
                  bestEdgeIt = it;
                } else if(vEnd == vertex) {
                  minDotProduct = dotProdMatrix(eLast, e);
                  bestEdgeIt = it;
                }
              }
            }

            if(bestEdgeIt == openEdges.end()) {
              writeVTKFileOfCell(cellId, &faces[startSetIndex], &vertices[startSetIndex], globalTimeStep);
              cerr << grid().domainId() << "open edges for cell " << cellId << ": ";
              cerr << "Open edges size: " << openEdges.size() << endl;
              for(auto& openEdge : openEdges) {
                cerr << openEdge.first;
              }
              mTerm(1, AT_, "Open edge for cutCell!");
            }

            // continue with the bestEdge
            e = (*bestEdgeIt).first;
            dir = (*bestEdgeIt).second;
            MInt vStart = edges[startSetIndex][e].vertices[0];
            MInt vEnd = edges[startSetIndex][e].vertices[1];
            if(dir == -1) {
              MInt tmp = vStart;
              vStart = vEnd;
              vEnd = tmp;
            }
            if(vStart != vertex) {
              if(vEnd == vertex) {
                // Switch direction of the edge
                if(dir == -1) {
                  dir = 1;
                } else {
                  dir = -1;
                }
                MInt vEndOld = vEnd;
                vEnd = vStart;
                vStart = vEndOld;
              } else {
                cerr << grid().domainId() << " missmatching edges for " << cellId << " "
                     << grid().tree().globalId(cellId) << " vertices: " << vStart << " " << vEnd << " " << vertex
                     << endl;
              }
            }

            eLast = e;

            // add the bestEdge to the face and remove it from the openEdges-list!
            faces[startSetIndex][noFaces].edges.emplace_back(e, dir);
            faces[startSetIndex][noFaces].vertices.push_back(vStart);

#ifdef CutCell_DEBUG
            if(grid().domainId() == debugDomainId && globalTimeStep == debugTimeStep && cellId == debugCellId) {
              cerr << "Adding edge " << e << " with vertex " << vStart << " and " << vEnd << endl;
              cerr << "The new face now has " << faces[startSetIndex][noFaces].vertices.size() << " vertices!" << endl;
            }
#endif

            openEdges.erase(bestEdgeIt);
            vertex = vEnd;
            faceVertices.push_back(vertex);
            for(MInt surfaceIdentificator : vertices[startSetIndex][vertex].surfaceIdentificators)
              surfaceIdentificatorCounters[surfaceIdentificator]++;
            edges[startSetIndex][e].face[1] = noFaces;
          }
          // ASSERT(faceVertices.size() > 2, "");
          if(faceVertices.size() <= 2) {
            cerr << grid().domainId() << " face with only two vertices for " << cellId << " "
                 << grid().tree().globalId(cellId) << " face : " << noFaces << endl;
          }

          // check for vertices which are doubled in face:
          MInt noDoubleVertices = 0;
          ASSERT(vertices[startSetIndex].size() <= (unsigned)maxNoVertices, "");
          for(MInt v = 0; (unsigned)v < vertices[startSetIndex].size(); v++) {
            vertexTouches[v] = 0;
          }
          for(MInt fv = 0; (unsigned)fv < faceVertices.size(); fv++) {
            vertexTouches[faceVertices[fv]]++;
          }
          for(MInt v = 0; (unsigned)v < vertices[startSetIndex].size(); v++) {
            if(vertexTouches[v] > 1) {
              noDoubleVertices++;
            }
          }

          if(noDoubleVertices != 0) {
            cerr << "double vertices for " << cellId << " on " << grid().domainId() << endl;
            error = true;
          }

          // find out correct bodyId/faceId and faceType and set normal and w
          // for the new face, based on the bestIdentificator!
          MInt bestIdentificator = -1;
          MInt maxHits = 0;
          for(MInt i = 0; i < m_noDirs + m_noEmbeddedBodies; i++) {
            if(surfaceIdentificatorCounters[i] > maxHits) {
              maxHits = surfaceIdentificatorCounters[i];
              bestIdentificator = i;
            }
          }
          ASSERT(bestIdentificator > -1, "");
          ASSERT(bestIdentificator < m_noDirs + m_noEmbeddedBodies, "");
          faceType = (bestIdentificator < m_noDirs ? 0 : 1);
          if(faceType == 0) {
            faceId = bestIdentificator;
          } else {
            bodyId = bestIdentificator - m_noDirs;
          }
          faces[startSetIndex][noFaces].faceType = faceType;
          faces[startSetIndex][noFaces].faceId = faceId;
          faces[startSetIndex][noFaces].bodyId = bodyId;

          // find correct face to inject w and normal from:
          MBool partnerFound = false;
          for(MInt e2 = 0; (unsigned)e2 < faces[startSetIndex][noFaces].edges.size(); e2++) {
            MInt edge = faces[startSetIndex][noFaces].edges[e2].first;
            MInt otherFace = edges[startSetIndex][edge].face[0];
            MInt otherFaceType = faces[startSetIndex][otherFace].faceType;
            MInt otherFaceId = faces[startSetIndex][otherFace].faceId;
            MInt otherBodyId = faces[startSetIndex][otherFace].bodyId;
            if(otherFaceType == faceType && otherFaceId == faceId && otherBodyId == bodyId) {
              partnerFound = true;
              faces[startSetIndex][noFaces].normal[0] = faces[startSetIndex][otherFace].normal[0];
              faces[startSetIndex][noFaces].normal[1] = faces[startSetIndex][otherFace].normal[1];
              faces[startSetIndex][noFaces].normal[2] = faces[startSetIndex][otherFace].normal[2];
              faces[startSetIndex][noFaces].w = faces[startSetIndex][otherFace].w;

              break;
            }
          }

          ASSERT(partnerFound, "");

          noFaces++;
        }

      } // loop though all other sets!

    } // if ( noIndividualCuts > 1 )

#ifdef CutCell_DEBUG
    if(grid().domainId() == debugDomainId && globalTimeStep == debugTimeStep && cellId == debugCellId) {
      cerr << "Final count: " << faces[startSetIndex].size() << " faces " << edges[startSetIndex].size() << " edges "
           << vertices[startSetIndex].size() << " vertices " << endl;
    }
#endif

    // debug check:
    // edge and vertex-information should match!
    for(MInt faceId = 0; faceId < (signed)faces[startSetIndex].size(); faceId++) {
      ASSERT(faces[startSetIndex][faceId].edges.size() == faces[startSetIndex][faceId].vertices.size(), "");
      for(MInt edgeId = 0; edgeId < (signed)faces[startSetIndex][faceId].edges.size(); edgeId++) {
        const MInt edge = faces[startSetIndex][faceId].edges[edgeId].first;
        const MInt newVertexId = faces[startSetIndex][faceId].vertices[edgeId];
        const MInt direction = faces[startSetIndex][faceId].edges[edgeId].second;
        MInt vertexId = edges[startSetIndex][edge].vertices[0];
        if(direction == -1) vertexId = edges[startSetIndex][edge].vertices[1];
        ASSERT(newVertexId == vertexId, "");
      }
    }

    RECORD_TIMER_STOP(tCutFace_3b);
    RECORD_TIMER_START(tCutFace_4);

    // 4. build polyhedron(polyhedra) structure
    // use a stack
    ASSERT(faceStack.empty(), "");

    // add original cell, if zero faces are found
    if(faces[startSetIndex].size() == 0) {
      cutCells.emplace_back(cellId, &a_coordinate(cellId, 0));
    }

    // add multiple cutCells if cutCell <= -1
    for(MInt faceCounter = 0; (unsigned)faceCounter < faces[startSetIndex].size(); faceCounter++) {
      if(faces[startSetIndex][faceCounter].cutCell > -1) continue;
      // NOTE: first faces are bodySurfaces, as cartesian surfaces were added afterwards
      // don't add lineFaces at the beginning!
      if(faces[startSetIndex][faceCounter].isLine) continue;

#ifdef CutCell_DEBUG
      if(grid().domainId() == debugDomainId && globalTimeStep == debugTimeStep && cellId == debugCellId) {
        cerr << "faceCounter: " << faceCounter << endl;
      }
#endif
      // adding cutCell
      faceStack.push(faceCounter);
      const MInt currentCutCell = cutCells.size();
      cutCells.emplace_back(cellId, &a_coordinate(cellId, 0));
      faces[startSetIndex][faceCounter].cutCell = currentCutCell;
      cutCells[currentCutCell].faces.push_back(faceCounter);
      // adding all other faces which share an edge to the current cutCell/bodySurface
      while(!faceStack.empty()) {
        MInt currentFace = faceStack.top();
        faceStack.pop();
        // find neighboring faces
        for(MInt e = 0; (unsigned)e < faces[startSetIndex][currentFace].edges.size(); e++) {
          MInt edge = faces[startSetIndex][currentFace].edges[e].first;
          MInt otherFace = edges[startSetIndex][edge].face[0];
          if(otherFace == currentFace) {
            otherFace = edges[startSetIndex][edge].face[1];
          }
          // if the neighboring face doesn't already belong to a different cutCell
          // add it to this one!
          if(faces[startSetIndex][otherFace].cutCell == -1) {
#ifdef CutCell_DEBUG
            if(grid().domainId() == debugDomainId && globalTimeStep == debugTimeStep && cellId == debugCellId) {
              cerr << "faceCounter: " << faceCounter << " cutCell " << currentCutCell
                   << " is adding neighbor: " << otherFace << endl;
            }
#endif
            cutCells[currentCutCell].faces.push_back(otherFace);
            faces[startSetIndex][otherFace].cutCell = currentCutCell;
            faceStack.push(otherFace);
          }
        }
      }
    }
    ASSERT(faceStack.empty(), "");

    RECORD_TIMER_STOP(tCutFace_4);

    // 5. compute  polyhedron(polyhedra)
    RECORD_TIMER_START(tCutFace_5a);
    compVolumeIntegrals_pyraBased3(&cutCells, &faces[startSetIndex], &vertices[startSetIndex]);
    RECORD_TIMER_STOP(tCutFace_5a);

#ifdef CutCell_DEBUG
    if(globalTimeStep == debugTimeStep && cellId == debugCellId && grid().domainId() == debugDomainId) {
      writeVTKFileOfCell(cellId, &faces[startSetIndex], &vertices[startSetIndex], -73);
    }
#endif

    if(error) {
      cerr << "[" << -1 << "]"
           << ": Warning: removal/suppression of doubly touched vertices seems to have failed!" << endl;
      writeVTKFileOfCell(cellId, &faces[startSetIndex], &vertices[startSetIndex], -1);
    }


    RECORD_TIMER_START(tCutFace_5b);

    // --------------------------------------------------------------
    // --------5) faceJoin
    // --------------------------------------------------------------

    // NOTE: it is not possible to easyly join all bodySurfaces with the same bodyId!
    //      as disconnected faces with the same bodyId can occour and must not be joined!
    //      -> a connective information between bodyFaces through edges is necessary!

    // 5.1. if required, join cut surfaces if their normal vectors are similar enough
    if(m_bodyFaceJoinMode != 0) { // until line 4577

      // set edgeCutCellPointer to -1 for deleted edges
      for(MInt e = 0; (unsigned)e < edges[startSetIndex].size(); e++) {
        MInt face = edges[startSetIndex][e].face[0];
        if(face < 0)
          edgeCutCellPointer[e] = -1;
        else
          edgeCutCellPointer[e] = faces[startSetIndex][face].cutCell;
      }
      for(MInt f = 0; (unsigned)f < faces[startSetIndex].size(); f++) {
        multiFaceConnection[f] = -1;
      }
      for(MInt c = 0; (unsigned)c < cutCells.size(); c++) {
        switch(m_bodyFaceJoinMode) {
          case 1: // all connected body faces with the same bodyId are joined
          case 3: // only cartesian faces are joined
          case 4: // connect body faces with the same bodyId only if the resulting face is not concave
          {
            list<MInt> bodyVertices;
            set<MInt> locBodySrfcs;
            ScratchSpace<MInt> isBodyEdge(maxNoVertices, AT_, "isBodyEdge");
            isBodyEdge.fill(-2);

            // remove douple edges entries in vertices
            for(auto& vert : vertices[startSetIndex]) {
              sort(vert.edges.begin(), vert.edges.end());
              auto last = unique(vert.edges.begin(), vert.edges.end());
              vert.edges.erase(last, vert.edges.end());
            }

            for(MInt e = 0; (unsigned)e < edges[startSetIndex].size(); e++) {
              if(edgeCutCellPointer[e] != c) continue;
              if(faces[startSetIndex][edges[startSetIndex][e].face[0]].faceType == 0
                 && faces[startSetIndex][edges[startSetIndex][e].face[1]].faceType == 0) {
                if(faces[startSetIndex][edges[startSetIndex][e].face[0]].faceId
                   == faces[startSetIndex][edges[startSetIndex][e].face[1]].faceId) {
                  pureBodyEdges.push_back(e);
                }
              } else if(m_bodyFaceJoinMode == 1 || m_bodyFaceJoinMode == 4) {
                if(faces[startSetIndex][edges[startSetIndex][e].face[0]].faceType == 1
                   && faces[startSetIndex][edges[startSetIndex][e].face[1]].faceType == 1) {
                  if(faces[startSetIndex][edges[startSetIndex][e].face[0]].bodyId
                     == faces[startSetIndex][edges[startSetIndex][e].face[1]].bodyId) {
                    pureBodyEdges.push_back(e);
                    isBodyEdge(e) = faces[startSetIndex][edges[startSetIndex][e].face[0]].bodyId;
                    bodyVertices.push_back(edges[startSetIndex][e].vertices[0]);
                    bodyVertices.push_back(edges[startSetIndex][e].vertices[1]);
                    locBodySrfcs.insert(faces[startSetIndex][edges[startSetIndex][e].face[0]].bodyId);
                  } else {
                    isBodyEdge(e) = -1;
                  }
                }
              }
            }

            if(m_bodyFaceJoinMode == 4 && noIndividualCuts > 1) {
              bodyVertices.sort();
              bodyVertices.unique();
              const MInt noLocBodySrfcs = (signed)locBodySrfcs.size();
              MInt concaveCnt = 0;
              for(MInt vert : bodyVertices) {
                if(vertices[startSetIndex][vert].pointType > 1) {
                  MInt noOuterEdges = 0;
                  MInt noBodyEdges = 0;
                  MInt outerEdges[6] = {-1};
                  MInt bodyEdges[6] = {-1};
                  for(MInt edge : vertices[startSetIndex][vert].edges) {
                    if(isBodyEdge(edge) == -1) {
                      outerEdges[noOuterEdges] = edge;
                      noOuterEdges++;
                    }
                    if(isBodyEdge(edge) > -1) {
                      bodyEdges[noBodyEdges] = edge;
                      noBodyEdges++;
                    }
                  }
                  if(noOuterEdges == 2 && noBodyEdges > 0) {
                    for(MInt k = 0; k < noBodyEdges; k++) {
                      MInt edge = bodyEdges[k];
                      MInt otherVert = (edges[startSetIndex][edge].vertices[0] == vert) ? 1 : 0;
                      otherVert = edges[startSetIndex][edge].vertices[otherVert];
                      MInt face0 = edges[startSetIndex][edge].face[0];
                      MInt face1 = edges[startSetIndex][edge].face[1];

                      MInt edge0 = -1;
                      MInt edge1 = -1;
                      MInt dir0 = -1;
                      MInt dir1 = -1;
                      MInt noEdges0 = (signed)faces[startSetIndex][face0].edges.size();
                      MInt noEdges1 = (signed)faces[startSetIndex][face1].edges.size();
                      MInt prevEdge = -1;
                      MInt nextEdge = -1;
                      MInt side = -1;
                      for(MInt e = 0; e < noEdges0; e++) {
                        if(faces[startSetIndex][face0].edges[e].first == edge) {
                          edge0 = e;
                          dir0 = faces[startSetIndex][face0].edges[e].second;
                          if(prevEdge < 0 || nextEdge < 0) {
                            MInt vA = edges[startSetIndex][edge].vertices[dir0 == 1 ? 1 : 0];
                            MInt otherEdge = (vA == vert) ? (edge0 + 1) % noEdges0 : (edge0 + noEdges0 - 1) % noEdges0;
                            otherEdge = faces[startSetIndex][face0].edges[otherEdge].first;
                            if(otherEdge == outerEdges[0] || otherEdge == outerEdges[1]) {
                              MInt otherEdge2 = (otherEdge == outerEdges[0]) ? outerEdges[1] : outerEdges[0];
                              nextEdge = (vA == vert) ? otherEdge : otherEdge2;
                              prevEdge = (vA == vert) ? otherEdge2 : otherEdge;
                              side = 0;
                            }
                          }
                        }
                      }
                      for(MInt e = 0; e < noEdges1; e++) {
                        if(faces[startSetIndex][face1].edges[e].first == edge) {
                          edge1 = e;
                          dir1 = faces[startSetIndex][face1].edges[e].second;
                          if(prevEdge < 0 || nextEdge < 0) {
                            MInt vA = edges[startSetIndex][edge].vertices[dir1 == 1 ? 1 : 0];
                            MInt otherEdge = (vA == vert) ? (edge1 + 1) % noEdges1 : (edge1 + noEdges1 - 1) % noEdges1;
                            otherEdge = faces[startSetIndex][face1].edges[otherEdge].first;
                            if(otherEdge == outerEdges[0] || otherEdge == outerEdges[1]) {
                              MInt otherEdge2 = (otherEdge == outerEdges[0]) ? outerEdges[1] : outerEdges[0];
                              nextEdge = (vA == vert) ? otherEdge : otherEdge2;
                              prevEdge = (vA == vert) ? otherEdge2 : otherEdge;
                              side = 1;
                            }
                          }
                        }
                      }
                      if(edge0 < 0 || edge1 < 0) mTerm(1, AT_, "ERROR: pure body edge not found.");
                      if(dir0 == dir1) mTerm(1, AT_, "ERROR: pure body edge not found (2).");
                      if(faces[startSetIndex][face0].edges[edge0].first != edge
                         || faces[startSetIndex][face1].edges[edge1].first != edge)
                        mTerm(1, AT_, "ERROR: pure body edge not found (3).");
                      ASSERT(outerEdges[0] > -1 && outerEdges[1] > -1, "");

                      if(prevEdge > -1 && nextEdge > -1) {
                        ASSERT(side > -1, "");
                        MInt id_p0 = edges[startSetIndex][prevEdge].vertices[0] == vert ? 0 : 1;
                        MInt id_p1 = edges[startSetIndex][nextEdge].vertices[0] == vert ? 0 : 1;
                        if(edges[startSetIndex][prevEdge].vertices[id_p0] != vert)
                          mTerm(1, AT_, "ERROR: pure body edge not found (4).");
                        if(edges[startSetIndex][nextEdge].vertices[id_p1] != vert)
                          mTerm(1, AT_, "ERROR: pure body edge not found (5).");
                        MInt otherId[2] = {1, 0};
                        MInt vA2 = edges[startSetIndex][prevEdge].vertices[otherId[id_p0]];
                        MInt vB2 = edges[startSetIndex][prevEdge].vertices[id_p0];
                        MInt vA = edges[startSetIndex][nextEdge].vertices[id_p1];
                        MInt vB = edges[startSetIndex][nextEdge].vertices[otherId[id_p1]];

                        MInt vV = vert;
                        MInt vO = otherVert;

                        MFloat abs1 = F0;
                        MFloat abs2 = F0;
                        MFloat abs3 = F0;
                        MFloat dotp = F0;
                        MFloat dotp2 = F0;
                        MFloat dotp3 = F0;
                        MFloat dotp4 = F0;
                        MFloat vec1[nDim] = {F0, F0, F0};
                        MFloat vec2[nDim] = {F0, F0, F0};
                        MFloat vec3[nDim] = {F0, F0, F0};
                        MFloat ori[nDim] = {F0, F0, F0};
                        MFloat normal[3] = {F0, F0, F0};
                        MFloat area0 = faces[startSetIndex][face0].area;
                        MFloat area1 = faces[startSetIndex][face1].area;
                        ASSERT(area0 + area1 > F0, "");
                        for(MInt i = 0; i < nDim; i++) {
                          normal[i] = (area0 * faces[startSetIndex][face0].normal[i]
                                       + area1 * faces[startSetIndex][face1].normal[i])
                                      / (area0 + area1);
                        }
                        for(MInt i = 0; i < nDim; i++) {
                          MFloat d1 =
                              vertices[startSetIndex][vB2].coordinates[i] - vertices[startSetIndex][vA2].coordinates[i];
                          MFloat d2 =
                              vertices[startSetIndex][vB].coordinates[i] - vertices[startSetIndex][vA].coordinates[i];
                          MFloat d3 =
                              vertices[startSetIndex][vO].coordinates[i] - vertices[startSetIndex][vV].coordinates[i];
                          vec1[i] = d1;
                          vec2[i] = d2;
                          vec3[i] = d3;
                          dotp += d1 * d2;
                          abs1 += POW2(d1);
                          abs2 += POW2(d2);
                          abs3 += POW2(d3);
                        }
                        abs1 = sqrt(abs1);
                        abs2 = sqrt(abs2);
                        abs3 = sqrt(abs3);
                        dotp /= (abs1 * abs2);
                        for(MInt i = 0; i < nDim; i++) {
                          vec1[i] /= abs1;
                          vec2[i] /= abs2;
                          vec3[i] /= abs3;
                        }
                        crossProduct(ori, vec1, vec2);
                        for(MInt i = 0; i < nDim; i++) {
                          dotp2 += ori[i] * normal[i];
                        }
                        crossProduct(ori, vec1, vec3);
                        for(MInt i = 0; i < nDim; i++) {
                          dotp3 += ori[i] * normal[i];
                        }
                        crossProduct(ori, vec2, vec3);
                        for(MInt i = 0; i < nDim; i++) {
                          dotp4 += ori[i] * normal[i];
                        }

                        if(dotp2 < -0.1 && fabs(F1 - dotp) > 1e-2
                           && mMin(area0, area1) > 1e-5 * pow(cellLength0, (MFloat)(nDim - 1))) {
                          if(noLocBodySrfcs + concaveCnt > maxNoSurfaces) {
                            cerr << "Warning: removal of concave polygons was skipped, since "
                                    "FvBndryCell<nDim>::m_maxNoSurfaces is not large enough."
                                 << endl;
                            continue;
                          }

                          if(dotp3 > -1e-3 && dotp4 > -1e-3) {
                            auto it2 = pureBodyEdges.begin();
                            while(it2 != pureBodyEdges.end()) {
                              // for( std::list<MInt>::iterator it2 = pureBodyEdges.begin(); it2 !=
                              // pureBodyEdges.end(); it2++){
                              MInt bedge = (*it2);
                              if(edgeCutCellPointer[bedge] != c) {
                                it2++;
                                continue;
                              }
                              MInt faceA = edges[startSetIndex][bedge].face[0];
                              MInt faceB = edges[startSetIndex][bedge].face[1];
                              if((face0 == faceA && face1 == faceB) || (face1 == faceA && face0 == faceB)) {
                                it2 = pureBodyEdges.erase(it2);
                                concaveCnt++;
                              } else {
                                it2++;
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

            multiFaces.clear();

            // connect multifaces
            MBool somethingChanged = true;
            MInt currentMultiFace = 0;
            while(!pureBodyEdges.empty()) {
              MInt pbedge = pureBodyEdges.front();
              pureBodyEdges.pop_front();
              multiFaces.emplace_back();
              MInt pbface0 = edges[startSetIndex][pbedge].face[0];
              multiFaceConnection[pbface0] = currentMultiFace;
              multiFaces[currentMultiFace].faces.push_back(pbface0);
              multiFaces[currentMultiFace].bodyId = faces[startSetIndex][pbface0].bodyId;
              MInt pbface1 = edges[startSetIndex][pbedge].face[1];
              multiFaceConnection[pbface1] = currentMultiFace;
              multiFaces[currentMultiFace].faces.push_back(pbface1);
              for(MInt e = 0; (unsigned)e < faces[startSetIndex][pbface0].edges.size(); e++) {
                if(faces[startSetIndex][pbface0].edges[e].first != pbedge) {
                  multiFaces[currentMultiFace].edges.emplace_back(faces[startSetIndex][pbface0].edges[e].first,
                                                                  faces[startSetIndex][pbface0].edges[e].second);
                } else { // insert edges of other face
                  // find index of pbedge in face1
                  MInt edgeIndex = 0;
                  for(MInt ee = 0; (unsigned)ee < faces[startSetIndex][pbface1].edges.size(); ee++) {
                    if(faces[startSetIndex][pbface1].edges[ee].first == pbedge) {
                      edgeIndex = ee;
                      break;
                    }
                  }
                  for(MInt ee = 1; (unsigned)ee < faces[startSetIndex][pbface1].edges.size(); ee++) {
                    MInt eee = (edgeIndex + ee) % faces[startSetIndex][pbface1].edges.size();
                    multiFaces[currentMultiFace].edges.emplace_back(faces[startSetIndex][pbface1].edges[eee].first,
                                                                    faces[startSetIndex][pbface1].edges[eee].second);
                  }
                }
              }
              somethingChanged = true;

              while(somethingChanged) {
                somethingChanged = false;
                auto it = pureBodyEdges.begin();
                while(it != pureBodyEdges.end()) {
                  // for( std::list<MInt>::iterator it=pureBodyEdges.begin();it != pureBodyEdges.end(); it++){
                  MInt edge = (*it);
                  if(edgeCutCellPointer[edge] != c) {
                    it++;
                    continue;
                  }
                  MInt face0 = edges[startSetIndex][edge].face[0];
                  MInt face1 = edges[startSetIndex][edge].face[1];
                  if(multiFaceConnection[face0] == currentMultiFace && multiFaceConnection[face1] == currentMultiFace) {
                    it = pureBodyEdges.erase(it);
                  } else {
                    if(multiFaceConnection[face0] == currentMultiFace) {
                      multiFaceConnection[face1] = currentMultiFace;
                      it = pureBodyEdges.erase(it);
                      multiFaces[currentMultiFace].faces.push_back(face1);
                      auto it2 = multiFaces[currentMultiFace].edges.begin();
                      for(it2 = multiFaces[currentMultiFace].edges.begin();
                          it2 != multiFaces[currentMultiFace].edges.end();
                          it2++) {
                        MInt edge2 = (*it2).first;
                        if(edge2 == edge) {
                          break;
                        }
                      } // delete edge, insert edges of other face
                      // find corresponding edgeIndex of new face
                      MInt edgeIndex = 0;
                      for(MInt e = 0; (unsigned)e < faces[startSetIndex][face1].edges.size(); e++) {
                        if(faces[startSetIndex][face1].edges[e].first == edge) {
                          edgeIndex = e;
                          break;
                        }
                      }
                      for(MInt e = 1; (unsigned)e < faces[startSetIndex][face1].edges.size(); e++) {
                        MInt ee = (edgeIndex + e) % faces[startSetIndex][face1].edges.size();
                        multiFaces[currentMultiFace].edges.insert(
                            it2, make_pair(faces[startSetIndex][face1].edges[ee].first,
                                           faces[startSetIndex][face1].edges[ee].second));
                      }
                      multiFaces[currentMultiFace].edges.erase(it2);
                      somethingChanged = true;
                    } else if(multiFaceConnection[face1] == currentMultiFace) {
                      multiFaceConnection[face0] = currentMultiFace;
                      it = pureBodyEdges.erase(it);
                      multiFaces[currentMultiFace].faces.push_back(face0);
                      auto it2 = multiFaces[currentMultiFace].edges.begin();
                      for(it2 = multiFaces[currentMultiFace].edges.begin();
                          it2 != multiFaces[currentMultiFace].edges.end();
                          it2++) {
                        MInt edge2 = (*it2).first;
                        if(edge2 == edge) {
                          break;
                        }
                      } // delete edge, insert edges of other face
                      // find corresponding edgeIndex of new face
                      MInt edgeIndex = 0;
                      for(MInt e = 0; (unsigned)e < faces[startSetIndex][face0].edges.size(); e++) {
                        if(faces[startSetIndex][face0].edges[e].first == edge) {
                          edgeIndex = e;
                          break;
                        }
                      }
                      for(MInt e = 1; (unsigned)e < faces[startSetIndex][face0].edges.size(); e++) {
                        MInt ee = (edgeIndex + e) % faces[startSetIndex][face0].edges.size();
                        multiFaces[currentMultiFace].edges.insert(
                            it2, make_pair(faces[startSetIndex][face0].edges[ee].first,
                                           faces[startSetIndex][face0].edges[ee].second));
                      }
                      multiFaces[currentMultiFace].edges.erase(it2);
                      somethingChanged = true;
                    } else {
                      it++;
                    }
                  }
                }
              }
              currentMultiFace++;
            }
          } break;
          case 2: // body faces are only joined when their normal vectors are similar
          {
            for(MInt e = 0; (unsigned)e < edges[startSetIndex].size(); e++) {
              if(edgeCutCellPointer[e] != c) continue;
              if(faces[startSetIndex][edges[startSetIndex][e].face[0]].faceType == 1
                 && faces[startSetIndex][edges[startSetIndex][e].face[1]].faceType == 1) {
                if(faces[startSetIndex][edges[startSetIndex][e].face[0]].bodyId
                   == faces[startSetIndex][edges[startSetIndex][e].face[1]].bodyId) {
                  pureBodyEdges.push_back(e);
                }
              } else if(faces[startSetIndex][edges[startSetIndex][e].face[0]].faceType == 0
                        && faces[startSetIndex][edges[startSetIndex][e].face[1]].faceType == 0) {
                if(faces[startSetIndex][edges[startSetIndex][e].face[0]].faceId
                   == faces[startSetIndex][edges[startSetIndex][e].face[1]].faceId) {
                  pureBodyEdges.push_back(e);
                }
              }
            }

            const MInt faceNum = faces[startSetIndex].size();
            ASSERT(faceNum <= maxNoFaces, "");
            for(MInt i = 0; i < faceNum; i++)
              for(MInt j = 0; j < faceNum; j++)
                normalDotProduct(i, j) = F0;
            MInt noBodyFaces = 0;
            for(MInt f = 0; (unsigned)f < cutCells[c].faces.size(); f++) {
              MInt face = cutCells[c].faces[f];
              if(faces[startSetIndex][face].faceType == 1) bodyFaces[noBodyFaces++] = face;
            }
            for(MInt f = 0; f < noBodyFaces; f++) {
              MInt face1 = bodyFaces[f];
              for(MInt ff = 0; ff <= f; ff++) {
                MInt face2 = bodyFaces[ff];
                for(MInt i = 0; i < nDim; i++)
                  normalDotProduct(face1, face2) +=
                      faces[startSetIndex][face1].normal[i] * faces[startSetIndex][face2].normal[i];
                normalDotProduct(face1, face2) = F1 - normalDotProduct(face1, face2);
                normalDotProduct(face2, face1) = normalDotProduct(face1, face2);
              }
            }

            multiFaces.clear();

            // connect multifaces
            MBool somethingChanged = true;
            MInt currentMultiFace = 0;
            while(!pureBodyEdges.empty()) {
              MInt pbedge = pureBodyEdges.front();
              pureBodyEdges.pop_front();
              MInt pbface0 = edges[startSetIndex][pbedge].face[0];
              MInt pbface1 = edges[startSetIndex][pbedge].face[1];
              if(normalDotProduct(pbface0, pbface1) > m_bodyFaceJoinCriterion) continue;
              multiFaces.emplace_back();
              multiFaceConnection[pbface0] = currentMultiFace;
              multiFaces[currentMultiFace].faces.push_back(pbface0);
              multiFaces[currentMultiFace].bodyId = faces[startSetIndex][pbface0].bodyId;
              multiFaceConnection[pbface1] = currentMultiFace;
              multiFaces[currentMultiFace].faces.push_back(pbface1);
              for(MInt e = 0; (unsigned)e < faces[startSetIndex][pbface0].edges.size(); e++) {
                if(faces[startSetIndex][pbface0].edges[e].first != pbedge) {
                  multiFaces[currentMultiFace].edges.emplace_back(faces[startSetIndex][pbface0].edges[e].first,
                                                                  faces[startSetIndex][pbface0].edges[e].second);
                } else { // insert edges of other face
                  // find index of pbedge in face1
                  MInt edgeIndex = 0;
                  for(MInt ee = 0; (unsigned)ee < faces[startSetIndex][pbface1].edges.size(); ee++) {
                    if(faces[startSetIndex][pbface1].edges[ee].first == pbedge) {
                      edgeIndex = ee;
                      break;
                    }
                  }
                  for(MInt ee = 1; (unsigned)ee < faces[startSetIndex][pbface1].edges.size(); ee++) {
                    MInt eee = (edgeIndex + ee) % faces[startSetIndex][pbface1].edges.size();
                    multiFaces[currentMultiFace].edges.emplace_back(faces[startSetIndex][pbface1].edges[eee].first,
                                                                    faces[startSetIndex][pbface1].edges[eee].second);
                  }
                }
              }
              somethingChanged = true;

              while(somethingChanged) {
                somethingChanged = false;
                auto it = pureBodyEdges.begin();
                while(it != pureBodyEdges.end()) {
                  // for( std::list<MInt>::iterator it=pureBodyEdges.begin();it != pureBodyEdges.end(); it++){
                  MInt edge = (*it);
                  if(edgeCutCellPointer[edge] != c) {
                    it++;
                    continue;
                  }
                  MInt face0 = edges[startSetIndex][edge].face[0];
                  MInt face1 = edges[startSetIndex][edge].face[1];
                  if(multiFaceConnection[face0] == currentMultiFace && multiFaceConnection[face1] == currentMultiFace) {
                    it = pureBodyEdges.erase(it);
                  } else if(multiFaceConnection[face0] == currentMultiFace) {
                    it = pureBodyEdges.erase(it);
                    MBool mayBeConnected = true;
                    for(MInt nf = 0; (unsigned)nf < multiFaces[currentMultiFace].faces.size(); nf++) {
                      MInt nFace = multiFaces[currentMultiFace].faces[nf];
                      if(normalDotProduct(face1, nFace) > m_bodyFaceJoinCriterion) {
                        mayBeConnected = false;
                        break;
                      }
                    }
                    if(!mayBeConnected) continue;
                    multiFaceConnection[face1] = currentMultiFace;
                    multiFaces[currentMultiFace].faces.push_back(face1);
                    auto it2 = multiFaces[currentMultiFace].edges.begin();
                    for(it2 = multiFaces[currentMultiFace].edges.begin();
                        it2 != multiFaces[currentMultiFace].edges.end();
                        it2++) {
                      MInt edge2 = (*it2).first;
                      if(edge2 == edge) {
                        break;
                      }
                    } // delete edge, insert edges of other face
                    // find corresponding edgeIndex of new face
                    MInt edgeIndex = 0;
                    for(MInt e = 0; (unsigned)e < faces[startSetIndex][face1].edges.size(); e++) {
                      if(faces[startSetIndex][face1].edges[e].first == edge) {
                        edgeIndex = e;
                        break;
                      }
                    }
                    for(MInt e = 1; (unsigned)e < faces[startSetIndex][face1].edges.size(); e++) {
                      MInt ee = (edgeIndex + e) % faces[startSetIndex][face1].edges.size();
                      multiFaces[currentMultiFace].edges.insert(
                          it2, make_pair(faces[startSetIndex][face1].edges[ee].first,
                                         faces[startSetIndex][face1].edges[ee].second));
                    }
                    multiFaces[currentMultiFace].edges.erase(it2);
                    somethingChanged = true;
                  } else if(multiFaceConnection[face1] == currentMultiFace) {
                    it = pureBodyEdges.erase(it);
                    MBool mayBeConnected = true;
                    for(MInt nf = 0; (unsigned)nf < multiFaces[currentMultiFace].faces.size(); nf++) {
                      MInt nFace = multiFaces[currentMultiFace].faces[nf];
                      if(normalDotProduct(face0, nFace) > m_bodyFaceJoinCriterion) {
                        mayBeConnected = false;
                        break;
                      }
                    }
                    if(!mayBeConnected) continue;
                    multiFaceConnection[face0] = currentMultiFace;
                    multiFaces[currentMultiFace].faces.push_back(face0);
                    auto it2 = multiFaces[currentMultiFace].edges.begin();
                    for(it2 = multiFaces[currentMultiFace].edges.begin();
                        it2 != multiFaces[currentMultiFace].edges.end();
                        it2++) {
                      MInt edge2 = (*it2).first;
                      if(edge2 == edge) {
                        break;
                      }
                    } // delete edge, insert edges of other face
                    // find corresponding edgeIndex of new face
                    MInt edgeIndex = 0;
                    for(MInt e = 0; (unsigned)e < faces[startSetIndex][face0].edges.size(); e++) {
                      if(faces[startSetIndex][face0].edges[e].first == edge) {
                        edgeIndex = e;
                        break;
                      }
                    }
                    for(MInt e = 1; (unsigned)e < faces[startSetIndex][face0].edges.size(); e++) {
                      MInt ee = (edgeIndex + e) % faces[startSetIndex][face0].edges.size();
                      multiFaces[currentMultiFace].edges.insert(
                          it2, make_pair(faces[startSetIndex][face0].edges[ee].first,
                                         faces[startSetIndex][face0].edges[ee].second));
                    }
                    multiFaces[currentMultiFace].edges.erase(it2);
                    somethingChanged = true;
                  } else {
                    it++;
                  }
                }
              }
              currentMultiFace++;
            }
          } break;

          default:
            mTerm(1, AT_, "ERROR: invalid bodyFaceJoinMode specified. exiting...");
            break;
        }

        // join faces that belong to a multi face -> compute multi face parameters
        for(MInt mf = 0; (unsigned)mf < multiFaces.size(); mf++) {
          MFloat normal[3] = {F0, F0, F0};
          MFloat center[3] = {F0, F0, F0};
          MFloat area2 = F0;
          for(MInt f = 0; (unsigned)f < multiFaces[mf].faces.size(); f++) {
            const MInt face = multiFaces[mf].faces[f];
            MFloat faceArea = faces[startSetIndex][face].area;
            // changes due to faces with zero area, which become part of a multiface
            if(faceArea < m_eps) {
              faceArea = m_eps;
            }
            area2 += faceArea;
            for(MInt i = 0; i < nDim; i++) {
              normal[i] += faces[startSetIndex][face].normal[i] * faceArea;
              center[i] += faces[startSetIndex][face].center[i] * faceArea;
            }
          }
          MFloat absNorm = F0;
          for(MInt i = 0; i < nDim; i++) {
            absNorm += normal[i] * normal[i];
          }
          absNorm = sqrt(absNorm);
          for(MInt i = 0; i < nDim; i++) {
            multiFaces[mf].center[i] = center[i] / area2;
            ASSERT(!std::isnan(multiFaces[mf].center[i]), "");
            multiFaces[mf].normal[i] = normal[i] / absNorm;
            ASSERT(!std::isnan(multiFaces[mf].normal[i]), "");
          }
          multiFaces[mf].area = absNorm;
        }

        // recompute edges that belong to a multi face
        if(!multiFaces.empty()) {
          // remove all cutCell faces and only keep those without multiFaceConnection!
          MInt noFaces = 0;
          for(MInt f = 0; (unsigned)f < cutCells[c].faces.size(); f++) {
            MInt face = cutCells[c].faces[f];
            if(multiFaceConnection[face] == -1) tmp_faces[noFaces++] = cutCells[c].faces[f];
          }
          cutCells[c].faces.clear();
          for(MInt f = 0; f < noFaces; f++) {
            cutCells[c].faces.push_back(tmp_faces[f]);
          }
        }

        // add multifaces to the face-collector
        for(MInt mf = 0; (unsigned)mf < multiFaces.size(); mf++) {
          // add a new face to faces collector
          // copy properties from the multiface collector
          const MInt newFace = faces[startSetIndex].size();
          const MInt faceId = faces[startSetIndex][multiFaces[mf].faces[0]].faceId;
          const MInt faceType = faces[startSetIndex][multiFaces[mf].faces[0]].faceType;
          const MInt bodyId = faces[startSetIndex][multiFaces[mf].faces[0]].bodyId;

          faces[startSetIndex].emplace_back(faceId, faceType, bodyId);
          faces[startSetIndex][newFace].area = multiFaces[mf].area;
          for(MInt i = 0; i < nDim; i++) {
            faces[startSetIndex][newFace].center[i] = multiFaces[mf].center[i];
            faces[startSetIndex][newFace].normal[i] = multiFaces[mf].normal[i];
          }
          faces[startSetIndex][newFace].w = F0;
          faces[startSetIndex][newFace].cutCell = c;


          // preprocess edges: remove double entries:
          MBool somethingChanged = true;
          while(somethingChanged) {
            somethingChanged = false;
            // for( std::list<pair<MInt,MInt> >::iterator it=multiFaces[mf].edges.begin(); it !=
            // multiFaces[mf].edges.end(); it++){
            auto it = multiFaces[mf].edges.begin();
            while(it != multiFaces[mf].edges.end()) {
              auto itLast = multiFaces[mf].edges.end();
              if(it == multiFaces[mf].edges.begin()) {
                itLast = multiFaces[mf].edges.end();
              } else {
                itLast = it;
              }
              itLast--;
              ASSERT(itLast != it, "");
              if((it->first == itLast->first) && (it->second == -(itLast->second))) {
                multiFaces[mf].edges.erase(itLast);
                it = multiFaces[mf].edges.erase(it);
                somethingChanged = true;
                // it--;
              } else {
                it++;
              }
            }
          }


          for(auto& edge : multiFaces[mf].edges) {
            faces[startSetIndex][newFace].edges.emplace_back(edge.first, edge.second);
            MInt vertexId = edges[startSetIndex][edge.first].vertices[0];
            if(edge.second == -1) vertexId = edges[startSetIndex][edge.first].vertices[1];
            faces[startSetIndex][newFace].vertices.push_back(vertexId);
          }
          cutCells[c].faces.push_back(newFace);
        }
      }
    }
    RECORD_TIMER_STOP(tCutFace_5b);


    //---------------------------------------------------------------
    //---------- MULTI-CUT-CELL GENERATION FINISHED -----------------
    //---------------------------------------------------------------

    // Euler's polyhedral formula check:
    // check the correct count of vertices, edges and faces for each splitCell!
    for(MInt sc = 0; sc < (MInt)cutCells.size(); sc++) {
      if(!cutCells[sc].faces.empty()) {
        MUint noEdges = 0;
        MInt pcnt = 0;
        MInt vertexRemap[maxNoVertices];
        fill_n(vertexRemap, maxNoVertices, -1);
        for(MInt f = 0; (unsigned)f < cutCells[sc].faces.size(); f++) {
          MInt face = cutCells[sc].faces[f];
          noEdges += faces[startSetIndex][face].vertices.size();
          for(MInt e = 0; (unsigned)e < faces[startSetIndex][face].vertices.size(); e++) {
            const MInt vertex = faces[startSetIndex][face].vertices[e];
            if(vertexRemap[vertex] < 0) {
              vertexRemap[vertex] = pcnt;
              pcnt++;
            }
          }
        }
        // edges are counted twice (once for each face) => noEdges / 2!
        if(pcnt + cutCells[sc].faces.size() != 2 + noEdges / 2 || noEdges % 2 == 1) {
          cerr << grid().domainId() << "Euler's polyhedral formula (1) not fulfilled ("
               << pcnt + (signed)cutCells[sc].faces.size() - 2 - noEdges / 2 << ") " << pcnt << " "
               << cutCells[sc].faces.size() << " " << noEdges / 2 << " " << noEdges << " for cell " << cellId << " "
               << globalTimeStep << " " << faces[startSetIndex].size() << " " << cutCells.size() << " "
               << grid().tree().globalId(cellId) << " " << sc << endl;
          writeVTKFileOfCell(cellId, &faces[startSetIndex], &vertices[startSetIndex], -1);

          if(sc >= 1) {
            cerr << "Deleting splitchilds " << endl;
            cutCells.erase(cutCells.begin() + sc);
          }
        }
      }
    }
    // --------------------------------------------------------------
    // --------6) fill cutCellData-structure
    // --------------------------------------------------------------

    // each Cell, which is cut by at least one boundary has one cutCell
    // only splitCells have multiple cutCells!
    // these are cells where the remaining fluid volume is not connected!

    // 6. relate polyhedron(polyhedra) quantities to Cartesian cell quantities -> finish cell
    // 6.0. Split cells and split Cartesian need extra treatment -> check here!
    MInt splitCellId = -1;
    MInt noSplitChildren = cutCells.size();
    if(noSplitChildren > CC::maxSplitCells) {
      mTerm(1, AT_, "Too many split cells generated, see CutCell::maxSplitCells");
    }
    cutCellData[cutc].noSplitChilds = 0;


    if(noSplitChildren > 1) {
      // --------6.0.1) add splitChildren to the cutCellData-structure
      //                prepare split parent and split children:


      for(MInt sc = 0; sc < noSplitChildren; sc++) {
        MUint splitcutc = cutCellData.size();
        cutCellData[cutc].splitChildIds[cutCellData[cutc].noSplitChilds++] = splitcutc;
        cutCellData.emplace_back();
        //        MInt splitChildId = m_splitChilds[ splitCellId ][sc];
        //        cutCellData[splitcutc].cellId = splitChildId;
        cutCellData[splitcutc].cellId = -1;
        //        cutCellData[splitcutc].cutCellId = m_bndryCellIds->a[splitChildId];
        //        cutCellData[splitcutc].cutCellId = -1;
        cutCellData[splitcutc].splitParentId = cutc;
      }


      //--RELOC1--
      /*
        for( MInt f = 0; f < m_noDirs; f++ ){
        bndryCell->m_associatedSrfc[f] = -2;
        // if cell has existing surfaces, they should not point to this cell anymore. Will be corrected later
        // if other neighbor of surface is already non-existent, delete the surface
        MInt srfcId = m_cellSurfaceMapping[ cellId ][ f ];
        if ( srfcId > -1 ) {
          MInt nghbr0 = m_surfaces->a[srfcId].m_nghbrCellIds[0];
          MInt nghbr1 = m_surfaces->a[srfcId].m_nghbrCellIds[1];
          if( nghbr0 == cellId ){
            m_surfaces->a[srfcId].m_nghbrCellIds[0] = -2;
            nghbr0 = -2;
          }else if(nghbr1 == cellId ){
            m_surfaces->a[srfcId].m_nghbrCellIds[1] = -2;
            nghbr1 = -2;
          }
          if( nghbr0 < 0 && nghbr1 < 0 ){
            m_surfaces->a[ srfcId ].m_area[0] = F0;
            m_isActiveSurface[ srfcId ] = false;
            deleteSurface( srfcId );
          }
          m_cellSurfaceMapping[ cellId ][ f ] = -2;
        }
      }*/
    }
    //    splitCellIds(bndryId) = splitCellId;
    ASSERT(cutCellData[cutc].noSplitChilds == 0 || cutCellData[cutc].noSplitChilds > 1, "");

    if(noSplitChildren > 1) {
      for(MInt sc = 0; sc < noSplitChildren; sc++) {
        splitCellList[sc] = cutCellData[cutc].splitChildIds[sc];
      }
    } else {
      splitCellList[0] = cutc;
    }

    //    bndryCell->m_volume = F0;

    MInt noFacesPerCartesianFaceAll[6] = {0, 0, 0, 0, 0, 0};
    MInt noBodySurfacesAll = 0;

    // --------6.0.2) Loop over all CutCells

    for(MInt sc = 0; sc < noSplitChildren; sc++) {
      RECORD_TIMER_START(tCutFace_6);

      //      MInt splitChildId = cellId;
      //      MInt scBndryId = bndryId;
      //      FvBndryCell<nDim>* scBndryCell = &m_fvBndryCnd->m_bndryCells->a[scBndryId];
      //      if( splitCellId > -1 ){
      //        splitChildId =m_splitChilds[splitCellId][sc];
      //        scBndryId = m_bndryCellIds->a[splitChildId];
      //        scBndryCell = &m_fvBndryCnd->m_bndryCells->a[scBndryId];
      //      }
      //      ASSERT(scBndryId > -1 , " ");

      MInt cutCellId = splitCellList[sc];

      // if cell is a split child, set m_pointIsInside information correctly (needed for vtu output)
      // if( a_hasProperty(splitChildId, Cell::IsSplitChild) ){

      // --------6.0.3) set m_pointIsInside information correctly for splitChilds
      if(cutCellData[cutCellId].splitParentId > -1) {
        for(MInt corner = 0; corner < m_noCorners; corner++) {
          cutCellData[cutCellId].cornerIsInsideGeometry[0][corner] = true;
        }
        for(MInt f = 0; (unsigned)f < cutCells[sc].faces.size(); f++) {
          MInt face = cutCells[sc].faces[f];
          for(MInt e = 0; (unsigned)e < faces[startSetIndex][face].vertices.size(); e++) {
            const MInt vStart = faces[startSetIndex][face].vertices[e];
            if(vertices[startSetIndex][vStart].pointType == 0) {
              MInt corner = vertices[startSetIndex][vStart].pointId;
              cutCellData[cutCellId].cornerIsInsideGeometry[0][corner] = false;
            }
          }
        }
      }

      // --------6.0.4) count noFacesPerCartesianFace and noBodySurfaces
      //              for each cutCell and for the combined splitParent

      MInt noFacesPerCartesianFace[6] = {0, 0, 0, 0, 0, 0};
      MInt noBodySurfaces = 0;
      for(MInt f = 0; (unsigned)f < cutCells[sc].faces.size(); f++) {
        MInt face = cutCells[sc].faces[f];
        if(faces[startSetIndex][face].faceType == 0) {
          noFacesPerCartesianFace[faces[startSetIndex][face].faceId]++;
          noFacesPerCartesianFaceAll[faces[startSetIndex][face].faceId]++;
        } else {
          noBodySurfaces++;
          noBodySurfacesAll++;
        }
      }
      //      for(MInt i : noFacesPerCartesianFace){
      //        if( i > 1 ){

      //--RELOC2--
      /*
      a_hasProperty( cellId ,  Cell::HasSplitFace ) = true;
      // prepare cell/surfaces
      MInt id = splitFaceCellIds[scBndryId];
      if( id < 0){
        id = splitFaceCells.size();
        splitFaceCells.push_back(cellWithSplitFace(splitChildId));
        splitFaceCellIds[scBndryId] = id;
      }
      cellWithSplitFace* scsCell = &splitFaceCells[id];
      scsCell->splitFaces.push_back(splitCartesianFace(i));
      bndryCell->m_associatedSrfc[i] = -2;
      // if cell has existing surfaces, they should not point to this cell anymore. Will be corrected later
      // if other neighbor of surface is already non-existent, delete the surface
      MInt srfcId = m_cellSurfaceMapping[ splitChildId ][ i ];
      if ( srfcId > -1 ) {
        MInt nghbr0 = m_surfaces->a[srfcId].m_nghbrCellIds[0];
        MInt nghbr1 = m_surfaces->a[srfcId].m_nghbrCellIds[1];
        if( nghbr0 == splitChildId ){
          m_surfaces->a[srfcId].m_nghbrCellIds[0] = -2;
          nghbr0 = -2;
        }else if(nghbr1 == splitChildId ){
          m_surfaces->a[srfcId].m_nghbrCellIds[1] = -2;
          nghbr1 = -2;
        }
        if( nghbr0 < 0 && nghbr1 < 0 ){
          m_surfaces->a[ srfcId ].m_area[0] = F0;
          m_isActiveSurface[ srfcId ] = false;
          deleteSurface( srfcId );
        }
        m_cellSurfaceMapping[ splitChildId ][ i ] = -2;
      }
  */
      //        }
      //      }

      if(noBodySurfaces > maxNoSurfaces || noBodySurfaces > CC::maxNoBoundarySurfaces) {
        cerr << "more than FvBndryCell<nDim>::m_maxNoSurfaces cut surfaces detected for a cell. This is not "
                "implemented in MB framework yet. cell "
             << cellId << " " << grid().domainId() << " " << globalTimeStep << ", splitChild " << sc << " "
             << noBodySurfaces << " maximum " << maxNoSurfaces << " " << CC::maxNoBoundarySurfaces << " "
             << noSplitChildren << " " << faces[startSetIndex].size() << endl;
        writeVTKFileOfCell(cellId, &faces[startSetIndex], &vertices[startSetIndex], startSetIndex);
        mTerm(1, AT_,
              " more than 3 cut surfaces detected for a cell. This is not implemented in MB framework yet. "
              "exiting...");
      }

      // --------6.0.5) initialise & limit volume

      // 6.1. set bndry cell and surface properties
      // deactivate cells with no faces
      //--RELOC-3A--
      if(cutCells[sc].faces.empty()) {
        //        m_cellIsInactive[ splitChildId ] = true;
        //        a_hasProperty(splitChildId, Cell::IsOnCurrentMGLevel) = false;
        //        cutCells[sc].volume = F0;
        //        cutCells[sc].center[0] = a_coordinate(cellId, 0);
        //        cutCells[sc].center[1] = a_coordinate(cellId, 1);
        //        cutCells[sc].center[2] = a_coordinate(cellId, 2);

        cutCellData[cutCellId].volume = F0;

        //        removeSurfaces(splitChildId);
      }
      if(cutCells[sc].volume > a_gridCellVolume(cutCellId)) {
        cutCells[sc].volume = a_gridCellVolume(cutCellId);
      }

      // --------6.0.5) set externalFaces-Bool

      // non fluid sides:
      for(MInt f = 0; f < m_noDirs; f++) {
        if(noFacesPerCartesianFace[f] == 0) { // non fluid face
                                              //          scBndryCell->m_externalFaces[f] = true;
          cutCellData[cutCellId].externalFaces[f] = true;

          //--RELOC3--

          if(splitCellId < 0) {
            /*
                        MInt srfcId = m_cellSurfaceMapping[ cellId ][ f ];
                        if ( srfcId > -1 ) {
                          m_surfaces->a[ srfcId ].m_area[0] = F0;
                          m_isActiveSurface[ srfcId ] = false;
                          deleteSurface( srfcId );
                        }
                      */
            //            bndryCell->m_associatedSrfc[ f ] = -1;
            //            m_cutFaceArea[ bndryId ][ f ] = F0;
            //            cutCellData[cutc].cartFaceArea[ f ] = F0;
            /*
                        MInt nghbrId = a_neighborId( cellId ,  f );
                        MInt nghbrBndryId = -1;
                        if( nghbrId > -1 )
                          nghbrBndryId = m_bndryCellIds->a[ nghbrId ];
                        if ( nghbrBndryId > -1 ) {
                          m_bndryCells->a[ nghbrBndryId ].m_associatedSrfc[ m_revDir[f] ] = -1;
                          m_cutFaceArea[ nghbrBndryId ][ m_revDir[f] ] = F0;
                        }*/

            /*
            //            MInt cutNghbr = cutCellData[cutc].cutNghbrIds[f];
            //            if ( cutNghbr > -1 ) {
            //              if ( cutCellData[ cutNghbr ].cutNghbrIds[ m_revDir[f] ] != (signed)cutc ) {
            //                cerr << cellId << " " << bndryId << " " << cutNghbr << " " << cutc << " " <<
            m_bndryCellIds->a[ a_neighborId(cellId, f) ] << endl;
            //                mTerm(1,AT_,"Wrong cut cell neighbor.");
            //              }
            //  //          cutCellData[cutNghbr].cartFaceArea[ m_revDir[f] ] = F0;
            //            }
            */
          }
        }
      }

      // --------6.1) before reassigning cutCell-structure, copy the relevant information
      //              into the cutCellData-structure

      ASSERT(vertices[startSetIndex].size() <= (unsigned)maxNoVertices, "");
      for(MInt cp = 0; cp < maxNoVertices; cp++)
        newCutPointId[cp] = -1;
      for(MInt i = 0; (unsigned)i < vertices[startSetIndex].size(); i++)
        isPartOfBodySurface[i] = false;
      for(MInt e = 0; e < 2 * m_noEdges; e++) {
        // pointEdgeId[e] = bndryCell->m_srfcs[0]->m_cutEdge[e];
        pointEdgeId[e] = cutCellData[cutc].cutEdges[e];
      }
      noBodySurfaces = 0;

      ASSERT(cutCells[sc].volume >= F0, "");
      ASSERT(!(std::isnan(cutCells[sc].volume)), "");
      ASSERT(!(std::isnan(cutCells[sc].center[sc])), "");
      ASSERT(!(std::isnan(cutCells[sc].center[1])), "");
      ASSERT(!(std::isnan(cutCells[sc].center[2])), "");
      ASSERT(!(std::isinf(cutCells[sc].volume)), "");
      ASSERT(!(std::isinf(cutCells[sc].center[sc])), "");
      ASSERT(!(std::isinf(cutCells[sc].center[1])), "");
      ASSERT(!(std::isinf(cutCells[sc].center[2])), "");

      //      scBndryCell->m_volume = cutCells[sc].volume;
      cutCellData[cutCellId].volume = cutCells[sc].volume;

      //      if ( a_hasProperty(cellId, Cell::IsSplitCell) ) {
      //        ASSERT( a_hasProperty(splitChildId, Cell::IsSplitChild ), "" );
      //        bndryCell->m_volume += scBndryCell->m_volume;
      //      }
      for(MInt i = 0; i < nDim; i++) {
        //        scBndryCell->m_coordinates[i] = cutCells[sc].center[i] - a_coordinate( cellId ,  i );
        cutCellData[cutCellId].volumetricCentroid[i] = cutCells[sc].center[i] - a_coordinate(cellId, i);
      }

      //      cutCellData[cutCellId].noCutPoints = 0;
      cutCellData[cutCellId].noCartesianSurfaces = 0;
      cutCellData[cutCellId].noBoundarySurfaces = 0;

      // --------6.1.1) transform face information in necessary cutCellData-information!

      for(MInt f = 0; (unsigned)f < cutCells[sc].faces.size(); f++) {
        MInt face = cutCells[sc].faces[f];
        ASSERT(faces[startSetIndex][face].area >= F0, "");
        ASSERT(!(std::isnan(faces[startSetIndex][face].area)), "");
        ASSERT(!(std::isnan(faces[startSetIndex][face].center[0])), "");
        ASSERT(!(std::isnan(faces[startSetIndex][face].center[1])), "");
        ASSERT(!(std::isnan(faces[startSetIndex][face].center[2])), "");
        ASSERT(!(std::isinf(faces[startSetIndex][face].area)), "");
        ASSERT(!(std::isinf(faces[startSetIndex][face].center[0])), "");
        ASSERT(!(std::isinf(faces[startSetIndex][face].center[1])), "");
        ASSERT(!(std::isinf(faces[startSetIndex][face].center[2])), "");

        // --------6.1.1) for body-faces:
        //                boundarySurfaceArea, boundarySurfaceCentroid, boundarySurfaceNormal
        //                boundarySurfaceBndryCndId, boundarySurfaceBodyId, noBoundarySurfaces

        if(faces[startSetIndex][face].faceType == 1) { // body surface
          //          FvBndryCell<nDim>::BodySurface* bodySurface = scBndryCell->m_srfcs[noBodySurfaces];
          //          bodySurface->m_area = faces[startSetIndex][face].area;
          cutCellData[cutCellId].boundarySurfaceArea[noBodySurfaces] = faces[startSetIndex][face].area;

          for(MInt i = 0; i < nDim; i++) {
            //            bodySurface->m_coordinates[i] =  faces[startSetIndex][face].center[i];
            //            bodySurface->m_normalVector[i] =  faces[startSetIndex][face].normal[i];
            cutCellData[cutCellId].boundarySurfaceCentroid[noBodySurfaces][i] = faces[startSetIndex][face].center[i];
            cutCellData[cutCellId].boundarySurfaceNormal[noBodySurfaces][i] = faces[startSetIndex][face].normal[i];
          }
          //          bodySurface->m_noCutPoints = (MInt) faces[startSetIndex][face].edges.size();
          //          bodySurface->m_bndryCndId = m_movingBndryCndId;

          cutCellData[cutCellId].boundarySurfaceBndryCndId[noBodySurfaces] = -1; // m_movingBndryCndId;
          cutCellData[cutCellId].boundarySurfaceBodyId[noBodySurfaces] = faces[startSetIndex][face].bodyId;

          //          ASSERT(bodySurface->m_noCutPoints <= m_noEdges*2, "");
          for(MInt v = 0; (unsigned)v < faces[startSetIndex][face].vertices.size(); v++) {
            const MInt vertex = faces[startSetIndex][face].vertices[v];
            // check this!
            isPartOfBodySurface[vertex] = true;
            //          MInt cutPointId = vertices[startSetIndex][vertex].pointId;
            newCutPointId[vertex] = v;
            vertices[startSetIndex][vertex].cartSrfcId = noBodySurfaces;
            // end check this!

            //            for( MInt i = 0; i < nDim; i++ ) {
            //              bodySurface->m_cutCoordinates[v][i] =
            //              vertices[startSetIndex][vertex].coordinates[i]; cutCellData[cutCellId].cutPoints[
            //              cutCellData[cutCellId].noCutPoints ][i] =
            //              vertices[startSetIndex][vertex].coordinates[i];
            //            }
            //            if( cutPointId > -1 ){
            //              bodySurface->m_cutEdge[v] = pointEdgeId[cutPointId];
            //              cutCellData[cutCellId].cutEdges[ cutCellData[cutCellId].noCutPoints ] =
            //              pointEdgeId[cutPointId];
            //            }
            //            else{
            //              bodySurface->m_cutEdge[v] = -1;
            //              cutCellData[cutCellId].cutEdges[ cutCellData[cutCellId].noCutPoints ] = -1;
            //            }
            //            bodySurface->m_bodyId[v] = faces[startSetIndex][face].bodyId;
            //            cutCellData[cutCellId].cutSrfcId[ cutCellData[cutCellId].noCutPoints ] =
            //            noBodySurfaces; cutCellData[cutCellId].cutBodyIds[
            //            cutCellData[cutCellId].noCutPoints ] = faces[startSetIndex][face].bodyId;
            //            cutCellData[cutCellId].noCutPoints++;
            //            if ( cutCellData[cutCellId].noCutPoints > CC::maxNoCutPoints )
            //            mTerm(1,AT_,"More cutpoints than I can store...");
          }
          noBodySurfaces++;
          cutCellData[cutCellId].noBoundarySurfaces++;
          ASSERT(cutCellData[cutCellId].noBoundarySurfaces == noBodySurfaces, "");
        } else { // Cartesian face

          // --------6.1.1) for cartesian-faces:
          //                cartFaceCentroid, cartFaceArea, cartFaceDir, noCartesianSurfaces

          //--RELOC4--
          const MInt dirId = faces[startSetIndex][face].faceId;
          const MInt spaceId = dirId / 2;

          /*MInt srfcId = m_cellSurfaceMapping[ splitChildId ][ dirId ];
          if( noFacesPerCartesianFace[dirId] > 1 ){// check if we have multiple faces in this
          direction; if yes, we have to do something special! MInt sfCellId =
          splitFaceCellIds[scBndryId]; cellWithSplitFace* sfCell = &splitFaceCells[sfCellId];
            // search the right split surface of ssCell:
            MInt sf = -1;
            for( MInt s = 0; (unsigned) s < sfCell->splitFaces.size(); s++){
              MInt sfDir = sfCell->splitFaces[s].direction;
              if( sfDir == dirId ){ // found it!
                sf = s;
                break;
              }
            }
            ASSERT(sf > -1, "");
            ASSERT(srfcId < 0, "");
            // check if a surface should be generated - not between two halo cells!
            MBool createSrfc = true;
            if( a_hasProperty(cellId, Cell::IsHalo) ){
              if( grid().tree().hasNeighbor(cellId, dir)cellId, dirId)){
                MInt nghbrTmp = a_neighborId(cellId, dirId);
                if( a_hasProperty(nghbrTmp, Cell::IsHalo) ){
                  createSrfc = false;
                }
              }else{
                createSrfc = false;
              }
            }
            if( createSrfc ){
              srfcId = getNewSurfaceId();
              createSurfaceSplit( srfcId, cellId, splitChildId, dirId );
              sfCell->splitFaces[sf].srfcIds.push_back(srfcId);
              m_cellSurfaceMapping[splitChildId][dirId]=-2;//indicate split surface
            }
          }
          else if( splitCellId > -1 ){ // split cell, but no split surface -> generate a new
          surface! ASSERT(srfcId < 0, "");
            // check if a surface should be generated - not between two halo cells!
            MBool createSrfc = true;
            if( a_hasProperty(cellId, Cell::IsHalo) ){
              if( grid().tree().hasNeighbor(cellId, dir)cellId, dirId)){
                MInt nghbrTmp = a_neighborId(cellId, dirId);
                if( a_hasProperty(nghbrTmp, Cell::IsHalo) ){
                  createSrfc = false;
                }
              }else{
                createSrfc = false;
              }
            }
            if( createSrfc ){
              srfcId = getNewSurfaceId();
              createSurfaceSplit( srfcId, cellId, splitChildId, dirId );
              m_cellSurfaceMapping[splitChildId][dirId]=srfcId;
            }
          }*/
          MFloat sign = F1;
          if(dirId % 2 == 0) sign = -F1;

          //          if ( srfcId > -1 ) {
          for(MInt i = 0; i < nDim; i++) {
            //              m_surfaces->a[ srfcId ].m_coordinates[ i ] =
            //              faces[startSetIndex][face].center[i];
            cutCellData[cutCellId].cartFaceCentroid[cutCellData[cutCellId].noCartesianSurfaces][i] =
                faces[startSetIndex][face].center[i];
          }
          //            m_surfaces->a[ srfcId ].m_coordinates[ spaceId ] = a_coordinate(cellId,
          //            spaceId ) + sign * cellHalfLength;
          cutCellData[cutCellId].cartFaceCentroid[cutCellData[cutCellId].noCartesianSurfaces][spaceId] =
              a_coordinate(cellId, spaceId) + sign * cellHalfLength;

          //            m_surfaces->a[ srfcId ].m_area[0] = faces[startSetIndex][face].area;
          cutCellData[cutCellId].cartFaceArea[cutCellData[cutCellId].noCartesianSurfaces] =
              faces[startSetIndex][face].area;
          cutCellData[cutCellId].cartFaceDir[cutCellData[cutCellId].noCartesianSurfaces] =
              faces[startSetIndex][face].faceId;

          //            if( noFacesPerCartesianFace[dirId] == 1 ) // do not set this for split faces
          //              scBndryCell->m_associatedSrfc[ dirId ] = srfcId;
          //          }
          /*else{
            if( grid().tree().hasNeighbor(cellId, dir)cellId, dirId) > 0 && (!a_hasProperty(cellId,
Cell::IsHalo)) ){ cerr << " domain " << domainId() << " cellId " << cellId << " (" <<
grid().tree().globalId(cellId) << ") bndryId " << bndryId << " split child " << splitChildId << "
has no surface in dir " << dirId << " (" << splitCellId <<")" << endl; writeVTKFileOfCutCell(cellId,
&cutCells, &faces[startSetIndex], &edges[startSetIndex], &vertices[startSetIndex], startSetIndex);
              outputPolyData(cellId, &cutCells, &faces[startSetIndex], &edges[startSetIndex],
&vertices[startSetIndex], startSetIndex); cerr << " nghbr: " << a_neighborId(cellId, dirId) << "
nghbr is bndry cell: " << a_boundaryId(a_neighborId(cellId, dirId)) << " nghbr is inactive: " <<
m_cellIsInactive[a_neighborId(cellId, dirId)] << " nghbr properties: "
              << a_hasProperty(a_neighborId(cellId, dirId), Cell::IsDummy)
              << a_isBndryGhostCell(a_neighborId(cellId, dirId))
              << a_hasProperty(a_neighborId(cellId, dirId), Cell::IsInterface)
              << a_isPeriodic(a_neighborId(cellId, dirId))
              << a_hasProperty(a_neighborId(cellId, dirId), Cell::IsCutOff)
              << a_hasProperty(a_neighborId(cellId, dirId), Cell::IsNotGradient)
              << a_hasProperty(a_neighborId(cellId, dirId), Cell::IsValidMGC)
              << a_hasProperty(a_neighborId(cellId, dirId), Cell::IsExchange)
              << a_hasProperty(a_neighborId(cellId, dirId), Cell::IsFlux)
              << a_hasProperty(a_neighborId(cellId, dirId), Cell::IsActive)
              << a_hasProperty(a_neighborId(cellId, dirId), Cell::IsOnCurrentMGLevel)
              << a_hasProperty(a_neighborId(cellId, dirId), Cell::IsInSpongeLayer)
              << a_hasProperty(a_neighborId(cellId, dirId), Cell::IsHalo)
              << a_hasProperty(a_neighborId(cellId, dirId), Cell::IsWindow)
              << a_hasProperty(a_neighborId(cellId, dirId), Cell::IsPeriodicWithRot)
              << a_hasProperty(a_neighborId(cellId, dirId), Cell::IsPartLvlAncestor)
              << a_hasProperty(a_neighborId(cellId, dirId), Cell::IsPartitionCell)
              << a_hasProperty(a_neighborId(cellId, dirId), Cell::IsSplitCell)
              << a_hasProperty(a_neighborId(cellId, dirId), Cell::HasSplitFace)
              << a_hasProperty(a_neighborId(cellId, dirId), Cell::IsTempLinked)
              << endl;
              mTerm(1,AT_,"This was not supposed to happen - where is my
surface?!?");
//               continue;
            }
            else
              continue;
          }
          m_cutFaceArea[ m_bndryCellIds->a[splitChildId] ][ dirId ] =
faces[startSetIndex][face].area;

          MInt nghbrId = m_surfaces->a[srfcId].m_nghbrCellIds[dirId%2];
          MInt nghbrBndryId = -1;
          if( splitCellId < 0 && noFacesPerCartesianFace[dirId] == 1){
            if ( nghbrId > -1 )
              nghbrBndryId = m_bndryCellIds->a[ nghbrId ];

            if ( nghbrBndryId > -1 ) {
              m_bndryCells->a[ nghbrBndryId ].m_associatedSrfc[ m_revDir[dirId] ] = srfcId;
            }
          }else{
            nghbrId = a_neighborId(cellId, dirId);
            if ( nghbrId > -1 )
              nghbrBndryId = m_bndryCellIds->a[ nghbrId ];
          }
          */
          cutCellData[cutCellId].noCartesianSurfaces++;
          if(cutCellData[cutCellId].noCartesianSurfaces > CC::maxNoCartesianSurfaces) {
            mTerm(1, AT_, "Increase CutCell::maxNoCartesianSurfaces");
          }
        } // end else
      }

      if(noSplitChildren > 1) {
        //        bndryCell->m_noSrfcs = 0;
        cutCellData[cutc].noBoundarySurfaces = 0;
      }
      //      scBndryCell->m_noSrfcs = noBodySurfaces;
      cutCellData[cutCellId].noBoundarySurfaces = noBodySurfaces;

      RECORD_TIMER_STOP(tCutFace_6);

      for(MInt k = 0; k < CC::noFaces; k++) {
        cutCellData[cutCellId].noFacesPerCartesianDir[k] = noFacesPerCartesianFace[k];
      }

      // --------6.2) transfer poly information to the cutCellData-structure

      RECORD_TIMER_START(tCutFace_7);

      //      m_noCutCellFaces[ scBndryId ] = 0;
      //      const MInt maxPointsXD = 12;
      cutCellData[cutCellId].noTotalFaces = 0;
      cutCellData[cutCellId].noAdditionalVertices = 0;
      MInt vertexMap[maxNoVertices];
      fill_n(vertexMap, maxNoVertices, -1);

      for(MInt f = 0; (unsigned)f < cutCells[sc].faces.size(); f++) {
        MInt face = cutCells[sc].faces[f];
        //        MInt nccf = m_noCutCellFaces[ scBndryId ];
        //      m_noCutFacePoints[ scBndryId ][ nccf ] = 0;
        const MInt facecnt = cutCellData[cutCellId].noTotalFaces;
        cutCellData[cutCellId].allFacesNoPoints[facecnt] = 0;
        cutCellData[cutCellId].allFacesBodyId[facecnt] = faces[startSetIndex][face].bodyId;

        for(MInt e = 0; (unsigned)e < faces[startSetIndex][face].vertices.size(); e++) {
          const MInt vertex = faces[startSetIndex][face].vertices[e];
          MInt pointId = vertices[startSetIndex][vertex].pointId;
          MInt storePointId = -1;
          if(vertices[startSetIndex][vertex].pointType == 0) {
            //            ASSERT(  m_noCutFacePoints[ scBndryId ][ nccf ] <  maxPointsXD, "");
            //            m_cutFacePointIds[ scBndryId ][ maxPointsXD * nccf + m_noCutFacePoints[ scBndryId ][ nccf ]] =
            //            vertices[startSetIndex][vertex].pointId ; m_noCutFacePoints[ scBndryId ][ nccf ] ++;
            ASSERT(pointId < CC::noCorners, "");
            storePointId = pointId;
          } else if(pointId > -1) {
            ASSERT(vertices[startSetIndex][vertex].pointType == 1, "");
            // if( !isPartOfBodySurface[vertex] )
            // continue;
            /* if( isPartOfBodySurface[vertex] ) {
              MInt numCutPointsPreviousFaces = 0;
              for( MInt srfc = 0; srfc < vertices[startSetIndex][vertex].cartSrfcId; srfc++ ){
                numCutPointsPreviousFaces += m_fvBndryCnd->m_bndryCells->a[scBndryId].m_srfcs[srfc]->m_noCutPoints;
              }
              MInt cp = numCutPointsPreviousFaces + newCutPointId[vertex] ;
              m_cutFacePointIds[ scBndryId ][ maxPointsXD*nccf + m_noCutFacePoints[ scBndryId ][ nccf ] ] =
            m_noCellNodes + (cp); m_noCutFacePoints[ scBndryId ][ nccf ] ++;
            }*/
            ASSERT(pointId < cutCellData[cutc].noCutPoints, "");
            storePointId = CC::noCorners + pointId;
          } else {
            ASSERT(vertices[startSetIndex][vertex].pointType > 1, "");
            if(vertexMap[vertex] < 0) {
              for(MInt i = 0; i < nDim; i++) {
                cutCellData[cutCellId].additionalVertices[cutCellData[cutCellId].noAdditionalVertices][i] =
                    vertices[startSetIndex][vertex].coordinates[i];
              }
              vertexMap[vertex] = cutCellData[cutCellId].noAdditionalVertices;
              cutCellData[cutCellId].noAdditionalVertices++;
              if(cutCellData[cutCellId].noAdditionalVertices > CC::maxNoAdditionalVertices) {
                writeVTKFileOfCell(cellId, &faces[startSetIndex], &vertices[startSetIndex], globalTimeStep);
              }
              ASSERT(cutCellData[cutCellId].noAdditionalVertices <= CC::maxNoAdditionalVertices, "");
            }
            storePointId = CC::noCorners + CC::maxNoCutPoints + vertexMap[vertex];
          }
          cutCellData[cutCellId].allFacesPointIds[facecnt][cutCellData[cutCellId].allFacesNoPoints[facecnt]] =
              storePointId;
          cutCellData[cutCellId].allFacesNoPoints[facecnt]++;
          ASSERT(cutCellData[cutCellId].allFacesNoPoints[facecnt] <= CC::maxNoFaceVertices,
                 "has: " + to_string(cutCellData[cutCellId].allFacesNoPoints[facecnt]) + " max is "
                     + to_string(CC::maxNoFaceVertices));
        }
        //        m_noCutCellFaces[scBndryId]++;
        cutCellData[cutCellId].noTotalFaces++;
      }
      // not yet adjusted to split cells!


      //--RELOC6--
      // set splitface stream for vtu output
      /*
      MBool needFaceStream = false;
      m_fvBndryCnd->bndryCells[ scBndryId ].m_faceVertices.clear();
      m_fvBndryCnd->bndryCells[ scBndryId ].m_faceStream.resize(0);
      if ( !a_hasProperty(cellId, Cell::IsSplitCell) ) {
        if ( a_hasProperty( cellId ,  Cell::HasSplitFace ) ) {
          needFaceStream = true;
        }
        else {
          for( MInt f = 0; (unsigned) f < cutCells[sc].faces.size(); f++ ){
            if ( needFaceStream ) break;
            MInt face = cutCells[sc].faces[f];
            polyFace* faceP = &faces[startSetIndex][face];
            for( MInt e = 0; (unsigned) e < faceP->edges.size(); e++ ){
              if ( needFaceStream ) break;
              MInt edge = faceP->edges[e].first;
              MInt v0 = edges[startSetIndex][edge].vertices[0];
              MInt v1 = edges[startSetIndex][edge].vertices[1];
              if ( vertices[startSetIndex][v0].pointType > 1 || vertices[startSetIndex][v1].pointType > 1 ) {
                //attention: these additional vertices may render the face polygon concave and paraview/openGL does not
      render concave polygons properly! needFaceStream = true;
              }
            }
          }
        }
      }

      if ( needFaceStream ) {
        m_fvBndryCnd->bndryCells[ scBndryId ].m_faceVertices.clear();
        m_fvBndryCnd->bndryCells[ scBndryId ].m_faceStream.resize(cutCells[sc].faces.size());
        MUint noEdges = 0;
        MBool vertexUsed[maxNoVertices] = {false};
        MUint fcnt = 0;
        for ( MInt faceType = 1; faceType >=0; faceType-- ) { //body faces first
          for( MInt f = 0; (unsigned) f < cutCells[sc].faces.size(); f++ ){
            MInt face = cutCells[sc].faces[f];
            polyFace* faceP = &faces[startSetIndex][face];
            if ( faceP->faceType != faceType ) continue;
            m_fvBndryCnd->bndryCells[ scBndryId ].m_faceStream[fcnt].resize(0);

            {
              MInt edge = faceP->edges[0].first;
              MInt dir = faceP->edges[0].second;
              MInt id0 = (dir==1) ? 0:1;
              MInt firstVertex = edges[startSetIndex][edge].vertices[id0];
              MInt previousVertex = firstVertex;
              MBool edgeChecked[maxNoEdges] = {false};
              for( MInt e = 0; (unsigned) e < faceP->edges.size(); e++ ){
                edge = faceP->edges[e].first;
                dir = faceP->edges[e].second;
                id0 = (dir==1) ? 0:1;
                MInt id1 = (dir==1) ? 1:0;
                MInt v0 = edges[startSetIndex][edge].vertices[id0];
                MInt v1 = edges[startSetIndex][edge].vertices[id1];
                ASSERT( v1 != v0,"");
                ASSERT( v0 == previousVertex, "" );
                m_fvBndryCnd->bndryCells[ scBndryId ].m_faceStream[fcnt].push_back( v1 );
                previousVertex = v1;
                if ( edgeChecked[edge] ) cerr << "duplicate edge!" << endl;
                edgeChecked[edge] = true;
              }
              ASSERT( previousVertex == firstVertex, "" );
            }

            noEdges += faceP->edges.size();

            MBool vertexUsed2[maxNoVertices] = {false};
            for ( MUint s = 0; s < m_fvBndryCnd->bndryCells[ scBndryId ].m_faceStream[fcnt].size(); s++ ) {
              if ( vertexUsed2[ m_fvBndryCnd->bndryCells[ scBndryId ].m_faceStream[fcnt][s] ] ) {
                cerr << "duplicate face point " << globalTimeStep << " " << grid().tree().globalId(cellId) << " " <<
      bndryCell->m_noSrfcs << " " << face << " " << faceP->faceId << " " << faceP->area
                     << " " << faceP->center[0] << " " << faceP->center[1] << " " << faceP->center[2]
                     << " " << faceP->normal[0] << " " << faceP->normal[1] << " " << faceP->normal[2]
                     << " " << m_fvBndryCnd->bndryCells[ scBndryId ].m_faceStream[fcnt][s] << endl;
                for( MInt e = 0; (unsigned) e < faceP->edges.size(); e++ ){
                  MInt edge = faceP->edges[e].first;
                  MInt dir = faceP->edges[e].second;
                  MInt id0 = dir==1?0:1;
                  MInt id1 = dir==1?1:0;
                  MInt v0 = edges[startSetIndex][edge].vertices[id0];
                  MInt v1 = edges[startSetIndex][edge].vertices[id1];
                  cerr << e << " " << dir << " " << edge << " " << v0 << " " << v1
                       << " " << edges[startSetIndex][edge].edgeType
                       << " " << edges[startSetIndex][edge].edgeId
                       << " " << edges[startSetIndex][edge].face[0]
                       << " " << edges[startSetIndex][edge].face[1]
                       << endl;
                }
              }
              vertexUsed[ m_fvBndryCnd->bndryCells[ scBndryId ].m_faceStream[fcnt][s] ] = true;
              vertexUsed2[ m_fvBndryCnd->bndryCells[ scBndryId ].m_faceStream[fcnt][s] ] = true;
            }

            fcnt++;
          }
        }
        MInt pcnt = 0;
        MInt vertexRemap[maxNoVertices] = {-1};
        for ( MUint v = 0; v < vertices[startSetIndex].size(); v++ ) {
          if ( !vertexUsed[v] ) continue;
          vertexRemap[v] = pcnt++;
          for ( MInt i = 0; i < nDim; i++ ) {
            m_fvBndryCnd->bndryCells[ scBndryId ].m_faceVertices.push_back( vertices[startSetIndex][v].coordinates[i] );
          }
        }

        for ( MUint f = 0; f < m_fvBndryCnd->bndryCells[ scBndryId ].m_faceStream.size(); f++ ) {
          for ( MUint v = 0; v < m_fvBndryCnd->bndryCells[ scBndryId ].m_faceStream[f].size(); v++ ) {
            MInt vid = m_fvBndryCnd->bndryCells[ scBndryId ].m_faceStream[f][v];
            m_fvBndryCnd->bndryCells[ scBndryId ].m_faceStream[f][v] = vertexRemap[vid];
          }
          for ( MUint v = 0; v < m_fvBndryCnd->bndryCells[ scBndryId ].m_faceStream[f].size(); v++ ) {
            if ( m_fvBndryCnd->bndryCells[ scBndryId ].m_faceStream[f][v] < 0
                 || m_fvBndryCnd->bndryCells[ scBndryId ].m_faceStream[f][v] >= (signed)(m_fvBndryCnd->bndryCells[
      scBndryId ].m_faceVertices.size()/3) ) { cerr << "vertex stream out of range" << endl;
            }
          }
        }
      }
      */

      // additional check!
      if(!cutCells[sc].faces.empty()) {
        MUint noEdges = 0;
        MInt pcnt = 0;
        MInt vertexRemap[maxNoVertices];
        fill_n(vertexRemap, maxNoVertices, -1);
        for(MInt f = 0; (unsigned)f < cutCells[sc].faces.size(); f++) {
          MInt face = cutCells[sc].faces[f];
          // noEdges += faces[startSetIndex][face].edges.size();
          noEdges += faces[startSetIndex][face].vertices.size();
          for(MInt e = 0; (unsigned)e < faces[startSetIndex][face].vertices.size(); e++) {
            const MInt vertex = faces[startSetIndex][face].vertices[e];
            if(vertexRemap[vertex] < 0) {
              vertexRemap[vertex] = pcnt;
              pcnt++;
            }
          }
        }
        if(pcnt + cutCells[sc].faces.size() != 2 + noEdges / 2 || noEdges % 2 == 1) {
          cerr << grid().domainId() << "Euler's polyhedral formula not fulfilled ("
               << pcnt + (signed)cutCells[sc].faces.size() - 2 - noEdges / 2 << ") " << pcnt << " "
               << cutCells[sc].faces.size() << " " << noEdges / 2 << " " << noEdges << " for cell " << cellId << " "
               << globalTimeStep << " " << faces[startSetIndex].size() << " " << cutCells.size() << " "
               << grid().tree().globalId(cellId) << " " << sc << endl;
          writeVTKFileOfCell(cellId, &faces[startSetIndex], &vertices[startSetIndex], -1);
        }
      }
      RECORD_TIMER_STOP(tCutFace_7);
    }


    if(noSplitChildren > 1) {
      cutCellData[cutc].volume = F0;
      for(MInt sc = 0; sc < noSplitChildren; sc++) {
        MInt cutCellId = splitCellList[sc];
        cutCellData[cutc].volume += cutCellData[cutCellId].volume;
      }
      for(MInt dir = 0; dir < m_noDirs; dir++) {
        cutCellData[cutc].noFacesPerCartesianDir[dir] = noFacesPerCartesianFaceAll[dir];
        //--REV_RELOC--
        if(cutCellData[cutc].noFacesPerCartesianDir[dir] == 0) {
          cutCellData[cutc].externalFaces[dir] = true;
        }
      }
    }


    /*
                  MBool hasSplitFace = false;
                  for ( MInt k = 0; k < CC::noFaces; k++ ) {
                    if ( noFacesPerCartesianFaceAll[k] > 1 ) hasSplitFace = true;
                  }
                  if( noSplitChildren > 1 || hasSplitFace ) {

                    cerr << "TESTF "<< grid().tree().globalId(cellId) << " " << hasAmbiguousFaces <<" " << noCutSets <<"
    " << maxCutsPerFace << " / " << noSplitChildren << " " << hasSplitFace
                         << " " << cutCellData[cutc].noBoundarySurfaces << " " <<
    cutCellData[cutc].boundarySurfaceBodyId[0] << endl;
                  }



                  //if ( grid().tree().globalId(cellId) == 205407 || grid().tree().globalId(cellId) == 205408 ||
    grid().tree().globalId(cellId) == 205411 || grid().tree().globalId(cellId) == 205418 ) { if (
    grid().tree().globalId(cellId) == 205407 ) { MInt splitParentId = cutCellData[cutc].splitParentId; MInt
    noSplitChilds = cutCellData[cutc].noSplitChilds; cerr << "compare " << grid().tree().globalId(cellId) << endl; if (
    cutcbak.cellId != cutCellData[cutc].cellId ) { cerr << grid().tree().globalId(cellId) << ": diff0 " << cutc << " "
    << noSplitChilds << " " << splitParentId << endl; continue;
                    }
                    if ( cutcbak.noSplitChilds != cutCellData[cutc].noSplitChilds ) {
                      cerr << grid().tree().globalId(cellId) << ": diff0a " << cutc << " " << noSplitChilds << " " <<
    splitParentId << endl;
                    }
                    if ( cutcbak.splitParentId != cutCellData[cutc].splitParentId ) {
                      cerr << grid().tree().globalId(cellId) << ": diff0b " << cutc << " " << noSplitChilds << " " <<
    splitParentId << endl;
                    }
                    MInt globalId =
    (splitParentId>-1?grid().tree().globalId(splitParentId):grid().tree().globalId(cellId)); if (
    fabs(cutCellData[cutc].volume - cutcbak.volume) > 1e-15 ) cerr << grid().tree().globalId(cellId) << ": diff1 " <<
    cutc << " " << noSplitChilds << " " << splitParentId << " " << cutCellData[cutc].volume  << " " <<  cutcbak.volume
    << endl; for ( MInt i = 0; i < nDim; i++ ) { if ( fabs(cutCellData[cutc].volumetricCentroid[i] -
    cutcbak.volumetricCentroid[i]) > 1e-15 ) cerr << grid().tree().globalId(cellId) << ": diff2 " << cutc << " " <<
    noSplitChilds << " " << splitParentId
                             << " " << i << " " << cutCellData[cutc].volumetricCentroid[i]/a_cellLengthAtCell(cellId) <<
    " " << cutcbak.volumetricCentroid[i]/a_cellLengthAtCell(cellId) << endl;
                    }

    //                cutCellData[cutc].volume = cutcbak.volume;
                    for ( MInt i = 0; i < nDim; i++ ) {
                      cerr << "set0 " << i << " " << cutCellData[cutc].volumetricCentroid[i] << " " <<
    cutcbak.volumetricCentroid[i]
                           << " " << a_coordinate(cellId,i)+cutCellData[cutc].volumetricCentroid[i]
                           << " " << a_coordinate(cellId,i)+cutcbak.volumetricCentroid[i] << endl;
                      cutCellData[cutc].volumetricCentroid[i] = cutcbak.volumetricCentroid[i];
                    }

                    MInt noSrfcs = cutcbak.noBoundarySurfaces;
                    if ( noSrfcs != cutCellData[cutc].noBoundarySurfaces ) {
                      cerr << grid().tree().globalId(cellId) << ": diff3 " << cutc << " " << noSplitChilds << " " <<
    (splitParentId>-1?grid().tree().globalId(splitParentId):-1) << " " << noSrfcs << " " <<
    cutCellData[cutc].noBoundarySurfaces << endl; continue;
                    }
                    for ( MInt bs = 0; bs < noSrfcs; bs++ ) {
                      for( MInt i = 0; i < nDim; i++ ){
                        if ( fabs( cutcbak.boundarySurfaceCentroid[bs][i] -
    cutCellData[cutc].boundarySurfaceCentroid[bs][i]) > 1e-15 ) cerr << grid().tree().globalId(cellId) << ": diff4 " <<
    cutc << " " << noSplitChilds << " " << splitParentId
                                << " " <<  i  << " " << cutcbak.boundarySurfaceCentroid[bs][i]  << " " <<
    cutCellData[cutc].boundarySurfaceCentroid[bs][i] << endl; if ( fabs( cutcbak.boundarySurfaceNormal[bs][i] -
    cutCellData[cutc].boundarySurfaceNormal[bs][i]) > 1e-15 ) { cerr << grid().tree().globalId(cellId) << ": diff5 " <<
    cutc << " " << globalId << " " << noSplitChilds << " " << splitParentId
                               << " " << i << " " << cutcbak.boundarySurfaceNormal[bs][i] << " " <<
    cutCellData[cutc].boundarySurfaceNormal[bs][i] << endl;
                        }
                      }
                      if ( fabs( cutcbak.boundarySurfaceArea[bs] - cutCellData[cutc].boundarySurfaceArea[bs]) > 1e-15 )
    cerr << grid().tree().globalId(cellId) << ": diff6 " << cutc << " " << noSplitChilds << " " << splitParentId <<
    endl; if ( cutcbak.boundarySurfaceBndryCndId[bs] != cutCellData[cutc].boundarySurfaceBndryCndId[bs] ) cerr <<
    grid().tree().globalId(cellId) << ": diff7 " << cutc << " " << noSplitChilds << " " << splitParentId << endl; if (
    cutcbak.boundarySurfaceBodyId[bs] != cutCellData[cutc].boundarySurfaceBodyId[bs] ) cerr <<
    grid().tree().globalId(cellId) << ": diff8 " << cutc << " " << noSplitChilds << " " << splitParentId << endl;

    //                  cutCellData[cutc].boundarySurfaceArea[bs] = cutcbak.boundarySurfaceArea[bs];
    //                  cutCellData[cutc].boundarySurfaceBndryCndId[bs] = cutcbak.boundarySurfaceBndryCndId[bs];
    //                  cutCellData[cutc].boundarySurfaceBodyId[bs] = cutcbak.boundarySurfaceBodyId[bs];
                      for( MInt i = 0; i < nDim; i++ ){
                        cerr << "set1 " << i << " " << cutCellData[cutc].boundarySurfaceCentroid[bs][i] << " " <<
    cutcbak.boundarySurfaceCentroid[bs][i] << endl; cutCellData[cutc].boundarySurfaceCentroid[bs][i] =
    cutcbak.boundarySurfaceCentroid[bs][i];
    //                    cutCellData[cutc].boundarySurfaceNormal[bs][i] = cutcbak.boundarySurfaceNormal[bs][i];
                      }

                      MInt noCP = 0;
                      for ( MInt c = 0; c < cutCellData[cutc].noCutPoints; c++ ) {
                        //if ( cutCellData[cutc].cutSrfcId[c] == bs ) {
                          for( MInt i = 0; i < nDim; i++ ) {
                            if ( fabs( cutcbak.cutPoints[c][i] - cutCellData[cutc].cutPoints[c][i]) > 1e-15 ) cerr <<
    grid().tree().globalId(cellId) << ": diff9 " << cutc << " " << noSplitChilds << " " << splitParentId << endl;
                          }
                          if ( cutcbak.cutEdges[c] != cutCellData[cutc].cutEdges[c] ) cerr <<
    grid().tree().globalId(cellId) << ": diff10 " << cutc << " " << noSplitChilds << " " << splitParentId << endl; if (
    cutcbak.cutBodyIds[c] != cutCellData[cutc].cutBodyIds[c] ) cerr << grid().tree().globalId(cellId) << ": diff11 " <<
    cutc << " " << noSplitChilds << " " << splitParentId << endl; noCP++;
                          //}
                      }
    //      if ( noCP != bodySurface->m_noCutPoints ) cerr << grid().tree().globalId(cellId) << ": diff12 " << globalId
    << " " << cutc << " " << noSrfcs << " " << noSplitChilds << " " << splitParentId
    //                                                     << " " << bs << " " << noCP << " " <<
    bodySurface->m_noCutPoints << " (" << cutCellData[cutc].noCutPoints << ")" << endl; for ( MInt f = 0; f < 2*nDim;
    f++ ) { if ( cutCellData[cutc].externalFaces[f] != cutcbak.externalFaces[f] )  cerr <<
    grid().tree().globalId(cellId) << ": diff13 " << cutc << " " << noSplitChilds << " " << splitParentId << endl;
                      }
                      for ( MInt c = 0; c < IPOW2(nDim); c++ ) {
                        if ( cutcbak.cornerIsInsideGeometry[0][c] != cutCellData[cutc].cornerIsInsideGeometry[0][c] )
    cerr << grid().tree().globalId(cellId) << ": diff14 " << cutc << " " << noSplitChilds << " " << splitParentId <<
    endl;
                      }

                    }

                    if ( cutcbak.noCartesianSurfaces != cutCellData[cutc].noCartesianSurfaces ) {
                      cerr << grid().tree().globalId(cellId) << ": diff15 " << cutc << " " << noSplitChilds << " " <<
    splitParentId << endl;
                    }
                    for ( MInt cs = 0; cs < cutCellData[cutc].noCartesianSurfaces; cs++ ) {
                      MBool match = false;
                      for ( MInt cs1 = 0; cs1 < cutcbak.noCartesianSurfaces; cs1++ ) {

                        if ( cutcbak.cartFaceDir[cs1] != cutCellData[cutc].cartFaceDir[cs] ) continue;

                        match = true;

                        if ( cutcbak.cartFaceDir[cs1] != cutCellData[cutc].cartFaceDir[cs] ) {
                          cerr << grid().tree().globalId(cellId) << ": diff16 " << cutc << " " << noSplitChilds << " "
    << splitParentId << endl;
                        }
                        if ( fabs(cutcbak.cartFaceArea[cs1] - cutCellData[cutc].cartFaceArea[cs]) > 1e-15 ) {
                          cerr << grid().tree().globalId(cellId) << ": diff17 " << cutc << " " << noSplitChilds << " "
    << splitParentId
                              << " " << cutcbak.cartFaceArea[cs1] << " " << cutCellData[cutc].cartFaceArea[cs] << endl;
                        }
                        for( MInt i = 0; i < nDim; i++ ) {
                          if ( fabs( cutcbak.cartFaceCentroid[cs1][i] - cutCellData[cutc].cartFaceCentroid[cs][i] ) >
    1e-15 ) { cerr << grid().tree().globalId(cellId) << ": diff18 " << cutc << " " << noSplitChilds << " " <<
    splitParentId
                                 << " " << i << " " <<  cutcbak.cartFaceCentroid[cs1][i] << " "
    <<cutCellData[cutc].cartFaceCentroid[cs][i] << endl;
                          }
                        }

                      }
                      if ( !match ) cerr << grid().tree().globalId(cellId) << ": diff16 " << endl;
                    }
                    for ( MInt cs = 0; cs < cutCellData[cutc].noCartesianSurfaces; cs++ ) {
                      cerr << cs << " " << cutcbak.cartFaceDir[cs] << " " << cutCellData[cutc].cartFaceDir[cs]
                          << " " << cutcbak.cartFaceArea[cs]  << " " << cutCellData[cutc].cartFaceArea[cs]
                          << " " << cutcbak.cartFaceCentroid[cs][0] << " " << cutCellData[cutc].cartFaceCentroid[cs][0]
                          << " " << cutcbak.cartFaceCentroid[cs][1] << " " << cutCellData[cutc].cartFaceCentroid[cs][1]
                          << " " << cutcbak.cartFaceCentroid[cs][2] << " " << cutCellData[cutc].cartFaceCentroid[cs][2]
                          << endl;
                    }



    //if ( grid().tree().globalId(cellId) == 205407 || grid().tree().globalId(cellId) == 205408 ||
    grid().tree().globalId(cellId) == 205411 || grid().tree().globalId(cellId) == 205418 ) {
    //if ( grid().tree().globalId(cellId) == 205407 ) {
    //                cutCellData[cutc] = cutcbak;
    //}
                  }
    */


    /*
                    const MInt faceCornerMapping[6][4] =
    {{2,6,4,0},{1,5,7,3},{0,4,5,1},{3,7,6,2},{2,3,1,0},{6,7,5,4}};
    //    for( MInt set = 0; set < m_noLevelSetsUsedForMb; set++ ){
                  { MInt set = 0;
                  MBool hasSplitFace = false;
                  MInt noPartiallyOutsideFaces = 0;
                  MInt noOutsideNodes = 0;
                  MInt maxCutsPerFace = 0;
                  MInt cutsPerFace[m_noCorners] = {0,0,0,0,0,0};
                  for( MInt c = 0; c < m_noCorners; c++){
                    if ( cutCellData[cutc].cornerIsInsideGeometry[set][c] ) noOutsideNodes++;
                  }
                  for ( MInt k = 0; k < CC::noFaces; k++ ) {
                    if ( noFacesPerCartesianFaceAll[k] > 1 ) hasSplitFace = true;
                    for ( MInt j = 0; j < 4; j++ ) {
                      if ( cutCellData[cutc].cornerIsInsideGeometry[set][faceCornerMapping[k][j]] ) {
                        noPartiallyOutsideFaces++;
                        break;
                      }
                    }
                  }
                  for( MInt cutPoint = 0; cutPoint < cutCellData[cutc].noCutPoints; cutPoint++){
                    MInt edge = cutCellData[cutc].cutEdges[cutPoint];
                    cutsPerFace[ edgeFaceCode[edge][0] ]++;
                    cutsPerFace[ edgeFaceCode[edge][1] ]++;
                  }
                  for ( MInt k = 0; k < CC::noFaces; k++ ) {
                    maxCutsPerFace = mMax(maxCutsPerFace, cutsPerFace[k] );
                  }
                  if( noSplitChildren > 1 || hasSplitFace ) {
                    cerr << "test_npof " << maxCutsPerFace  << " " << hasAmbiguousFaces << " " <<
    cutCellData[cutc].noCutPoints << " " << noOutsideNodes << endl; } else { cerr << "tset_npof " << maxCutsPerFace  <<
    " " << hasAmbiguousFaces << " " << cutCellData[cutc].noCutPoints << " " << noOutsideNodes << endl;
                  }
                  if ( !( noSplitChildren > 1 || hasSplitFace ) && hasAmbiguousFaces ) cerr << "case_0" << endl;
                  if ( ( noSplitChildren > 1 || hasSplitFace ) && !hasAmbiguousFaces ) cerr << "case_1" << endl;
        }*/

    //=================================================

#ifdef CutCell_DEBUG
    if(globalTimeStep == debugTimeStep && cellId == debugCellId && grid().domainId() == debugDomainId) {
      writeVTKFileOfCell(cellId, &faces[startSetIndex], &vertices[startSetIndex], -2);

      const MChar* fileName = "cell-Info_";
      stringstream fileName2;
      fileName2 << fileName << grid().tree().globalId(cellId) << ".txt";

      ofstream ofl;
      ofl.open((fileName2.str()).c_str(), ofstream::trunc);

      if(ofl) {
        ofl.setf(ios::fixed);
        ofl.precision(12);

        ofl << "Writing Cell " << grid().tree().globalId(cellId) << " " << cellId << " " << cutc << " noCuts "
            << noIndividualCuts << " noCutCells " << cutCells.size() << " vol " << cutCellData[cutc].volume
            << " SplitChilds " << cutCellData[cutc].noSplitChilds << " " << cutCellData[cutc].splitParentId << endl
            << endl;
        for(MInt dir = 0; dir < m_noDirs; dir++) {
          ofl << "Dir " << dir << " external " << cutCellData[cutc].externalFaces[dir] << " faces "
              << cutCellData[cutc].noFacesPerCartesianDir[dir] << endl
              << endl;
        }

        ofl << " Cartesian-Surf: " << cutCellData[cutc].noCartesianSurfaces << endl << endl;

        for(MInt id = 0; id < cutCellData[cutc].noCartesianSurfaces; id++) {
          ofl << "Id " << id << " dir " << cutCellData[cutc].cartFaceDir[id] << " area "
              << cutCellData[cutc].cartFaceArea[id] << " coord " << cutCellData[cutc].cartFaceCentroid[id][0]
              << cutCellData[cutc].cartFaceCentroid[id][1] << cutCellData[cutc].cartFaceCentroid[id][2] << endl
              << endl;
        }

        ofl << " Bndry-Surfaces: " << cutCellData[cutc].noBoundarySurfaces << endl << endl;

        for(MInt id = 0; id < cutCellData[cutc].noBoundarySurfaces; id++) {
          ofl << "Id " << id << " normal " << cutCellData[cutc].boundarySurfaceNormal[id][0]
              << cutCellData[cutc].boundarySurfaceNormal[id][1] << cutCellData[cutc].boundarySurfaceNormal[id][2]
              << " center " << cutCellData[cutc].boundarySurfaceCentroid[id][0]
              << cutCellData[cutc].boundarySurfaceCentroid[id][1] << cutCellData[cutc].boundarySurfaceCentroid[id][2]
              << " area " << cutCellData[cutc].boundarySurfaceArea[id] << endl
              << endl;
        }
      }
    }
#endif
  }


#ifdef CutCell_DEBUG
  if(globalTimeStep == (g_timeSteps + g_restartTimeStep)) returnDebugInfo();
#endif
}


//-----------------------------------------------------------------------------------------------------

// template <MInt nDim_>
// void GeometryIntersection<nDim_>::designateCutCellCandidates( std::vector< CutCell<nDim> >& cutCellData )
//{
//  clearCutCellData();
//  cutCellData.clear();
//
//  // loop over all cells
//  // --stl: check if neighbor exists or not in all 27 directions!
//  // --ls/analytical: check level set value vs threshold
//
//  //or stl: compute distance and then uniform scheme for all boundary types
//
//  //better stl: hierarchical, check minLevel...maxLevel, exclude children whose parent already too far from
//  boundary/not cut by boundary
//}

/**
 * \brief Computes cutpoints for the scalar-Field with grid edges and updates cutPoint-Info
 *        in the cutCandidates data!
 * \note can handle complex geometries for multiple scalar fields
 * \author Claudia Guenther, Update Tim Wegmann
 * \date 30.07.2013
 */
template <MInt nDim_>
void GeometryIntersection<nDim_>::computeCutPoints(std::vector<CutCandidate<nDim>>& candidates,
                                                   const MInt* candidateIds,
                                                   std::vector<MInt>& candidatesOrder) {
  TRACE();

  const MInt edgesPerDim = IPOW2(nDim_ - 1);

  const MInt DOFStencil[12] = {1, 1, 0, 0, 1, 1, 0, 0, 2, 2, 2, 2};

  // returns the two faces of any edge
  const MInt faceStencil[2][12] = {{0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 0, 1}, {4, 4, 4, 4, 5, 5, 5, 5, 2, 2, 3, 3}};

  // edege -> opposing edge
  // 1: oppsing edge on the same higher face count
  // 2: opposing edge on the same lower face count
  // 3: oppsong edge accros the cube
  //                                         0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11
  const MInt reverseEdgeStencil[3][12] = {{1, 0, 3, 2, 5, 4, 7, 6, 9, 8, 11, 10},
                                          {4, 5, 6, 7, 0, 1, 2, 3, 10, 11, 8, 9},
                                          {5, 4, 7, 6, 1, 0, 3, 2, 11, 10, 9, 8}};

  const MInt signStencil[3][12] = {{-1, 1, 0, 0, -1, 1, 0, 0, -1, 1, -1, 1},
                                   {0, 0, -1, 1, 0, 0, -1, 1, -1, -1, 1, 1},
                                   {-1, -1, -1, -1, 1, 1, 1, 1, 0, 0, 0, 0}};

  // edeg -> corner: returns the two corners for any edge
  //                                  0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11
  const MInt nodeStencil[2][12] = {{0, 1, 0, 2, 4, 5, 4, 6, 0, 1, 2, 3}, {2, 3, 1, 3, 6, 7, 5, 7, 4, 5, 6, 7}};

  //----------------------------------------------------------


  const MFloat eps = m_multiCutCell ? 10 * m_eps : m_eps;

  MInt errorFlag = 0;

  //----------------------------------------------------------

  // loop over       candidates,     determine cutpoints,    store cutCandidates properties
  for(MInt cndt = 0; cndt < (signed)candidates.size(); cndt++) {
    const MInt cellId = candidates[cndt].cellId;
    const MFloat cellLength = a_cellLengthAtCell(cellId);
    const MFloat cellHalfLength = F1B2 * cellLength;
    for(MInt edge = 0; edge < m_noEdges; edge++) {
      if(candidates[cndt].edgeChecked[edge]) continue;

      // first do some data structure related stuff
      MInt face[2]{-1, -1};
      // in 2D: face == edge
      face[0] = faceStencil[0][edge];

      IF_CONSTEXPR(nDim_ == 3) {
        // in 3D: each edge is connected to two faces
        face[1] = faceStencil[1][edge];
      }

      // edge's degree of freedom
      const MInt edgeDOF = DOFStencil[edge];

      // direct neighbours of each edge -> requires that neighboring cells are on same level!
      MInt directNeighbourIds[4]{-1, -1, -1, -1};
      directNeighbourIds[0] = cellId;

      if(grid().tree().hasNeighbor(cellId, face[0]) > 0) {
        directNeighbourIds[1] = grid().tree().neighbor(cellId, face[0]);
      }
      MInt reverseEdge[edgesPerDim];
      reverseEdge[0] = edge;
      reverseEdge[1] = reverseEdgeStencil[0][edge];

      IF_CONSTEXPR(nDim_ == 3) {
        if(grid().tree().hasNeighbor(cellId, face[1]) > 0) {
          directNeighbourIds[2] = grid().tree().neighbor(cellId, face[1]);
        }

        if(grid().tree().hasNeighbor(cellId, face[0]) > 0) {
          if(grid().tree().hasNeighbor(grid().tree().neighbor(cellId, face[0]), face[1]) > 0) {
            directNeighbourIds[3] = grid().tree().neighbor(grid().tree().neighbor(cellId, face[0]), face[1]);
          }
        }
        if(directNeighbourIds[3] < 0) {
          if(grid().tree().hasNeighbor(cellId, face[1]) > 0) {
            if(grid().tree().hasNeighbor(grid().tree().neighbor(cellId, face[1]), face[0]) > 0) {
              directNeighbourIds[3] = grid().tree().neighbor(grid().tree().neighbor(cellId, face[1]), face[0]);
            }
          }
        }
        reverseEdge[2] = reverseEdgeStencil[1][edge];
        reverseEdge[3] = reverseEdgeStencil[2][edge];
      }


      MInt startSet = 0;
      MInt endSet = 1;

      const MBool isGapCell = candidates[cndt].isGapCell;
      if(m_complexBoundary && m_noLevelSetsUsedForMb > 1 && (!isGapCell)) {
        startSet = 1;
        endSet = m_noLevelSetsUsedForMb;
      }


      for(MInt set = startSet; set < endSet; set++) {
        MFloat phi[2]{-99, -99}; // signed-distance value on both vertices defining an edge

        phi[0] = candidates[cndt].nodalValues[set][nodeStencil[0][edge]];
        phi[1] = candidates[cndt].nodalValues[set][nodeStencil[1][edge]];

        ASSERT(candidates[cndt].nodeValueSet[nodeStencil[0][edge]], "");
        ASSERT(candidates[cndt].nodeValueSet[nodeStencil[1][edge]], "");

        // prevent zero-volume cells
        if(fabs(phi[0]) < eps) {
          if(phi[0] < F0)
            phi[0] = -eps;
          else
            phi[0] = eps;
        }
        if(fabs(phi[1]) < eps) {
          if(phi[1] < F0)
            phi[1] = -eps;
          else
            phi[1] = eps;
        }

        // limit the distance of a cutPoint to any corner:
        // together with the above, this ensure that the CP
        // is always 2 eps appart from the corner!
        if(m_multiCutCell) {
          if(approx(fabs(phi[0]), cellHalfLength, eps)) {
            if(phi[0] > F0)
              phi[0] = cellHalfLength - eps;
            else
              phi[0] = -cellHalfLength + eps;
          }
          if(approx(fabs(phi[1]), cellHalfLength, eps)) {
            if(phi[1] > F0)
              phi[1] = cellHalfLength - eps;
            else
              phi[1] = -cellHalfLength + eps;
          }
        }

        // found a cutpoint if sign differs
        if(phi[0] * phi[1] < 0) {
          MFloat cutpoint[nDim_];
          // determine cutpoint coordinates
          const MFloat deltaCutpoint = (phi[0] + phi[1]) / (phi[0] - phi[1]);
          ASSERT(deltaCutpoint < F1 && deltaCutpoint > -F1, "");
          for(MInt dim = 0; dim < nDim_; dim++) {
            if(dim == edgeDOF) continue;
            cutpoint[dim] = a_coordinate(cellId, dim) + signStencil[dim][edge] * cellHalfLength;
          }
          cutpoint[edgeDOF] = a_coordinate(cellId, edgeDOF) + cellHalfLength * deltaCutpoint;

          // store this cutpoint at all direct neighbours
          for(MInt nb = 0; nb < edgesPerDim; nb++) {
            const MInt nbCell = directNeighbourIds[nb];
            if(nbCell < 0) continue;
            const MInt nbCandidate = candidateIds[nbCell];
            // IF_CONSTEXPR(nDim_ == 3) {
            // ASSERT(nbCandidate >= 0 , "ERROR: future cutCell has no valid cutCandidate-neighbor!");
            //} else {
            if(nbCandidate < 0) continue;
            //}

            ASSERT(!candidates[nbCandidate].edgeChecked[reverseEdge[nb]], "");

            // store the ordering in which bndryCells are added
            // the order appears to change the result for gapClosing-testcases!
            candidatesOrder.push_back(nbCandidate);

            // store relevant information
            const MInt noCPs = candidates[nbCandidate].noCutPoints;
            if(noCPs == 2 * m_noEdges) {
              mTerm(1, AT_, "Error: too many cut points, can't store more!");
            }
            if(candidates[nbCandidate].noCutPointsOnEdge[reverseEdge[nb]] == 2) {
              // check your geometry and gap-Cell!
              mTerm(1, AT_, "Error: already two cut points on edge, can't store more!");
            }


            if(set >= 1) {
              // ASSERT( candidates[ nbCandidate ].associatedBodyIds[set] ==
              //         candidates[ cndt ].associatedBodyIds[set] , "");
            }

            for(MInt dim = 0; dim < nDim_; dim++) {
              candidates[nbCandidate].cutPoints[noCPs][dim] = cutpoint[dim];
            }

            // check your geometry and gap-Cell!
            ASSERT(candidates[nbCandidate].associatedBodyIds[set] > -1, "cutPoint with invalid bodyId! ");

            candidates[nbCandidate].cutBodyIds[noCPs] = candidates[nbCandidate].associatedBodyIds[set];
            candidates[nbCandidate].cutEdges[noCPs] = reverseEdge[nb];
            candidates[nbCandidate].noCutPoints++;
            candidates[nbCandidate].noCutPointsOnEdge[reverseEdge[nb]]++;

          } // for nb
        }   // phi1*phi2<0
      }     // for sets

      // flag direct neighbours so that this edge is only checked once
      for(MInt nb = 0; nb < edgesPerDim; nb++) {
        if(directNeighbourIds[nb] < 0) continue;
        const MInt nbCndt = candidateIds[directNeighbourIds[nb]];
        if(nbCndt < 0) continue;
        candidates[nbCndt].edgeChecked[reverseEdge[nb]] = true;
      }
    } // for edge
  }   // for cndt

#ifdef CutCell_DEBUG
  MPI_Allreduce(MPI_IN_PLACE, &errorFlag, 1, MPI_INT, MPI_MAX, grid().mpiComm(), AT_, "MPI_IN_PLACE", "errorFlag");
#endif
  if(errorFlag) {
    mTerm(1, AT_, "Critical error in computeCutPoints! Quit.");
  }
}


/**
 * \brief
 * 1) add cutCandidates based on the solver m_bndryCandidates
 * 2) set cutCandidate information
 * 3) set the scalar-nodal values at the mesh vertices for the cutCandidates
 *
 * NOTE: the nodal-values need to be exchanged afterwards!
 *
 * \author Lennart Schneiders, Update Tim Wegmann
 */
template <MInt nDim_>
void GeometryIntersection<nDim_>::computeNodalValues(std::vector<CutCandidate<nDim>>& candidates, MInt* candidateIds,
                                                     const MFloat* scalarField, const MInt* bodyIdField,
                                                     const MBool* const gapPropertyField, const MBool gapClosure) {
  TRACE();

#ifdef CutCell_DEBUG
  const MInt debugTimeStep = -2;
  const MInt debugGlobalId = 216221;
#endif

  m_cutLvlJumpCandidates.clear();

  if(candidates.empty()) return;

  const MInt noCandidates = (signed)candidates.size();

  // add bndryLvlJump cell parents to the candidates and fill m_cutLvlJumpCandidates info
  if(m_bndryLvlJumps) {
    vector<std::pair<MInt, MInt>> diagonalJumpCells;
    for(MInt cnd = 0; cnd < noCandidates; cnd++) {
      const MInt cellId = candidates[cnd].cellId;
#ifdef CutCell_DEBUG
      if(globalTimeStep > debugTimeStep && grid().tree().globalId(cellId) == debugGlobalId) {
        cerr << "Cell is candidate!" << endl;
      }
#endif
      ASSERT(!candidates[cnd].isbndryLvlJumpParent, "");
      MInt id = -1;
      diagonalJumpCells.clear();
      const MInt parentId = grid().tree().parent(cellId);
      if(parentId < 0) continue;
      for(MInt dir = 0; dir < m_noDirs; dir++) {
        // check if the cell is a direct bndryLvlJumpCell:
        if(!grid().tree().hasNeighbor(cellId, dir) && grid().tree().hasNeighbor(parentId, dir)) {
          if(id < 0) {
            // added new as levelJump-candidate:
            id = m_cutLvlJumpCandidates.size();
            m_cutLvlJumpCandidates.emplace_back();
            m_cutLvlJumpCandidates[id].candId = cnd;
            // find childId of the cell in relation to its parent
            MInt child = -1;
            for(child = 0; child < m_maxNoChilds; child++) {
              const MInt childCellId = grid().tree().child(parentId, child);
              if(childCellId == cellId) break;
            }
            m_cutLvlJumpCandidates[id].childId = child;
            // add parent as additional cutCandidate once
            if(candidateIds[parentId] < 0) {
              const MInt newId = candidates.size();
              candidates.emplace_back();
              candidates[newId].cellId = parentId;
              candidates[newId].isbndryLvlJumpParent = true;
              candidateIds[parentId] = newId;
            } else {
              if(!candidates[candidateIds[parentId]].isbndryLvlJumpParent) {
                candidates[candidateIds[parentId]].isbndryLvlJumpParent = true;
              }
            }
            const MInt parentCandId = candidateIds[parentId];
            ASSERT(parentCandId > -1, "");
            m_cutLvlJumpCandidates[id].parentCandId = parentCandId;
          }
          // add jump-related information:
          const MInt noJumps = m_cutLvlJumpCandidates[id].noJumps;

#ifdef CutCell_DEBUG
          if(globalTimeStep > debugTimeStep && grid().tree().globalId(cellId) == debugGlobalId) {
            cerr << "Direct nieghbor Jump!" << dir << endl;
          }
#endif

          m_cutLvlJumpCandidates[id].dirs[noJumps] = dir;
          m_cutLvlJumpCandidates[id].diagonalDirs[noJumps] = -2;
          m_cutLvlJumpCandidates[id].neighborType[noJumps] = 0;

          m_cutLvlJumpCandidates[id].noJumps++;

        } else { // 2D-diagonal check
          const MInt nghbrId1 = grid().tree().neighbor(cellId, dir);
          if(nghbrId1 < 0) continue;
          if(!grid().tree().isLeafCell(nghbrId1)) continue;
          const MInt ngbParent1 = grid().tree().parent(nghbrId1);
          if(ngbParent1 < 0) continue;
          if(ngbParent1 == parentId) continue;
          for(MInt dir1 = 0; dir1 < m_noDirs; dir1++) {
            if((dir1 / 2) == (dir / 2)) continue;
            // if(dir1 == dir ) continue;
            if(!grid().tree().hasNeighbor(nghbrId1, dir1) && grid().tree().hasNeighbor(ngbParent1, dir1)) {
              // don't add the same neighbor twice!
              MBool alreadyAdded = false;
              MUint i = 0;
              for(i = 0; i < diagonalJumpCells.size(); i++) {
                if(diagonalJumpCells[i].first == grid().tree().neighbor(ngbParent1, dir1)
                   && diagonalJumpCells[i].second == 1) {
                  alreadyAdded = true;
                  break;
                } else if(diagonalJumpCells[i].first == grid().tree().neighbor(ngbParent1, dir1)) {
                  break;
                }
              }

              if(alreadyAdded) {
#ifdef CutCell_DEBUG
                if(globalTimeStep > debugTimeStep && grid().tree().globalId(cellId) == debugGlobalId) {
                  cerr << "Already found: " << dir << " " << dir1 << " "
                       << grid().tree().globalId(grid().tree().neighbor(ngbParent1, dir1)) << endl;
                }
#endif

                ASSERT(id > -1, "");
                continue;
              }

              // add diagonal neighbor to the list
              diagonalJumpCells.emplace_back(make_pair(grid().tree().neighbor(ngbParent1, dir1), 1));

              if(id < 0) {
                // added new as levelJump-candidate:
                id = m_cutLvlJumpCandidates.size();
                m_cutLvlJumpCandidates.emplace_back();
                m_cutLvlJumpCandidates[id].candId = cnd;
                // find childId of the cell in relation to its parent
                MInt child = -1;
                for(child = 0; child < m_maxNoChilds; child++) {
                  const MInt childCellId = grid().tree().child(parentId, child);
                  if(childCellId == cellId) break;
                }
                m_cutLvlJumpCandidates[id].childId = child;
                // add parent as additional cutCandidate once
                if(candidateIds[parentId] < 0) {
                  const MInt newId = candidates.size();
                  candidates.emplace_back();
                  candidates[newId].cellId = parentId;
                  candidates[newId].isbndryLvlJumpParent = true;
                  candidateIds[parentId] = newId;
                } else {
                  if(!candidates[candidateIds[parentId]].isbndryLvlJumpParent) {
                    candidates[candidateIds[parentId]].isbndryLvlJumpParent = true;
                  }
                }
                const MInt parentCandId = candidateIds[parentId];
                ASSERT(parentCandId > -1, "");
                m_cutLvlJumpCandidates[id].parentCandId = parentCandId;
              }
              const MInt noJumps = m_cutLvlJumpCandidates[id].noJumps;

#ifdef CutCell_DEBUG
              if(globalTimeStep > debugTimeStep && grid().tree().globalId(cellId) == debugGlobalId) {
                cerr << "Adding type1: " << dir << " " << dir1 << " "
                     << grid().tree().globalId(grid().tree().neighbor(ngbParent1, dir1)) << " " << noJumps << " "
                     << endl;
              }
#endif

              m_cutLvlJumpCandidates[id].neighborType[noJumps] = 1;
              // store diagonal relation
              m_cutLvlJumpCandidates[id].dirs[noJumps] = dir;
              m_cutLvlJumpCandidates[id].diagonalDirs[noJumps] = dir1;

              m_cutLvlJumpCandidates[id].noJumps++;

            } else
              IF_CONSTEXPR(nDim == 3) { // check for 3D-diagonals
                const MInt nghbrId2 = grid().tree().neighbor(nghbrId1, dir1);
                if(nghbrId2 < 0) continue;
                if(!grid().tree().isLeafCell(nghbrId2)) continue;
                const MInt ngbParent2 = grid().tree().parent(nghbrId2);
                if(ngbParent2 < 0) continue;
                if(ngbParent2 == parentId) continue;
                for(MInt dir2 = 0; dir2 < 2 * nDim; dir2++) {
                  if(((dir2 / 2) == (dir / 2)) || ((dir2 / 2) == (dir1 / 2))) continue;

                  if(!grid().tree().hasNeighbor(nghbrId2, dir2) && grid().tree().hasNeighbor(ngbParent2, dir2)) {
                    // don't add the same neighbor twice!
                    MBool alreadyAdded = false;
                    for(MUint i = 0; i < diagonalJumpCells.size(); i++) {
                      if(diagonalJumpCells[i].first == grid().tree().neighbor(ngbParent2, dir2)) {
                        alreadyAdded = true;
                        break;
                      }
                    }
                    if(alreadyAdded) {
                      ASSERT(id > -1, "");
                      continue;
                    }
                    // add diagonal neighbor to the list
                    diagonalJumpCells.emplace_back(make_pair(grid().tree().neighbor(ngbParent2, dir2), 2));
                    if(id < 0) {
                      // added new as levelJump-candidate:
                      id = m_cutLvlJumpCandidates.size();
                      m_cutLvlJumpCandidates.emplace_back();
                      m_cutLvlJumpCandidates[id].candId = cnd;
                      // find childId of the cell in relation to its parent
                      MInt child = -1;
                      for(child = 0; child < m_maxNoChilds; child++) {
                        const MInt childCellId = grid().tree().child(parentId, child);
                        if(childCellId == cellId) break;
                      }
                      m_cutLvlJumpCandidates[id].childId = child;
                      // add parent as additional cutCandidate once
                      if(candidateIds[parentId] < 0) {
                        const MInt newId = candidates.size();
                        candidates.emplace_back();
                        candidates[newId].cellId = parentId;
                        candidates[newId].isbndryLvlJumpParent = true;
                        candidateIds[parentId] = newId;
                      } else {
                        if(!candidates[candidateIds[parentId]].isbndryLvlJumpParent) {
                          candidates[candidateIds[parentId]].isbndryLvlJumpParent = true;
                        }
                      }
                      const MInt parentCandId = candidateIds[parentId];
                      ASSERT(parentCandId > -1, "");
                      m_cutLvlJumpCandidates[id].parentCandId = parentCandId;
                    }
                    const MInt noJumps = m_cutLvlJumpCandidates[id].noJumps;
                    m_cutLvlJumpCandidates[id].neighborType[noJumps] = 2;
                    // store diagonal relation
                    m_cutLvlJumpCandidates[id].dirs[noJumps] = dir;
                    m_cutLvlJumpCandidates[id].diagonalDirs[noJumps] = dir1;
                    m_cutLvlJumpCandidates[id].noJumps++;

#ifdef CutCell_DEBUG
                    if(globalTimeStep > debugTimeStep && grid().tree().globalId(cellId) == debugGlobalId) {
                      cerr << "Adding type2: " << dir << " " << dir1 << " " << dir2 << " "
                           << grid().tree().globalId(grid().tree().neighbor(ngbParent2, dir1)) << " " << noJumps
                           << endl;
                    }
#endif
                  }
                }
              }
          }
        }
      }
    }
  }

  ASSERT(noCandidates <= (signed)candidates.size(), "");

  const MInt nodeStencil[3][8] = {{0, 1, 0, 1, 0, 1, 0, 1}, {2, 2, 3, 3, 2, 2, 3, 3}, {4, 4, 4, 4, 5, 5, 5, 5}};

  const MInt reverseNode[3][8] = {{1, 0, 3, 2, 5, 4, 7, 6}, {2, 3, 0, 1, 6, 7, 4, 5}, {4, 5, 6, 7, 0, 1, 2, 3}};

  MInt nghbrIds[m_noCorners]{};
  MInt nghbrNodes[m_noCorners]{};

  for(MInt cnd = 0; cnd < (signed)candidates.size(); cnd++) {
    const MInt cellId = candidates[cnd].cellId;

    // setCandidate-properties:
    candidates[cnd].isGapCell = gapPropertyField[cellId];
    for(MInt set = 0; set < m_noLevelSetsUsedForMb; set++) {
      if(bodyIdField != nullptr) {
        candidates[cnd].associatedBodyIds[set] = bodyIdField[IDX_LSSETMB(cellId, set)];
      }
    }

    for(MInt node = 0; node < m_noCorners; node++) {
      if(candidates[cnd].nodeValueSet[node]) continue;
      // b) Add all neighbors and the corresponding nodes from the node-loop to the list

      MInt noNeighborsPerNode = 0;

      nghbrIds[noNeighborsPerNode] = cellId;
      nghbrNodes[noNeighborsPerNode] = node;
      noNeighborsPerNode++;

      for(MInt i = 0; i < nDim; i++) {
        const MInt firstDir = nodeStencil[i][node];
        const MInt firstNghbrNode = reverseNode[i][node];
        const MInt firstNghbrId = grid().tree().neighbor(cellId, firstDir);
        if(firstNghbrId > -1) {
          nghbrIds[noNeighborsPerNode] = firstNghbrId;
          nghbrNodes[noNeighborsPerNode] = firstNghbrNode;
          noNeighborsPerNode++;
          for(MInt j = i + 1; j < nDim; j++) {
            const MInt secondDir = nodeStencil[j][node];
            if(secondDir == firstDir) continue;
            const MInt secondNghbrNode = reverseNode[j][firstNghbrNode];
            const MInt secondNghbrId = grid().tree().neighbor(firstNghbrId, secondDir);
            if(secondNghbrId > -1) {
              nghbrIds[noNeighborsPerNode] = secondNghbrId;
              nghbrNodes[noNeighborsPerNode] = secondNghbrNode;
              noNeighborsPerNode++;
              IF_CONSTEXPR(nDim == 3) {
                for(MInt k = j + 1; k < nDim; k++) {
                  const MInt thirdDir = nodeStencil[k][node];
                  if(thirdDir == firstDir || thirdDir == secondDir) continue;
                  const MInt thirdNghbrNode = reverseNode[k][secondNghbrNode];
                  const MInt thirdNghbrId = grid().tree().neighbor(secondNghbrId, thirdDir);
                  if(thirdNghbrId > -1) {
                    nghbrIds[noNeighborsPerNode] = thirdNghbrId;
                    nghbrNodes[noNeighborsPerNode] = thirdNghbrNode;
                    noNeighborsPerNode++;
                  }
                }
              }
            }
          }
        }
      }

      MBool gapBoundary = false;
      MBool setBoundary = false;

      // ensure the same NodeValues at Gap-NoGap-Node-Boundaries!
      // important for 3D-CutCells!
      // Otherwise, the number of CutPoints doesn't match the nodalLSValues!
      // however, this means that nodes of cells, which are not a gapCell but a gapBoundary node
      // will still be using the effective G0-set!
      if(gapClosure) {
        const MBool isGapCell = gapPropertyField[cellId];
        for(MInt nghbrNode = 0; nghbrNode < noNeighborsPerNode; nghbrNode++) {
          const MBool isGapNghbr = gapPropertyField[nghbrIds[nghbrNode]];

          if(isGapCell != isGapNghbr) {
            gapBoundary = true;
          }
        }
        const MInt bodyId = bodyIdField[IDX_LSSETMB(cellId, 0)];
        for(MInt nghbrNode = 0; nghbrNode < noNeighborsPerNode; nghbrNode++) {
          const MInt bodyIdNghbr = bodyIdField[IDX_LSSETMB(nghbrIds[nghbrNode], 0)];
          if(bodyId != bodyIdNghbr) setBoundary = true;
        }
      }

      for(MInt set = 0; set < m_noLevelSetsUsedForMb; set++) {
        // interpolate scalar field value
        MFloat phi = F0;
        for(MInt nghbrNode = 0; nghbrNode < noNeighborsPerNode; nghbrNode++) {
          MInt usedSet = set;

          if(gapBoundary) {
            const MInt bodyIdNghbr = bodyIdField[IDX_LSSETMB(nghbrIds[nghbrNode], 0)];
            // const MInt bodyId = bodyIdField[IDX_LSSETMB( cellId, 0) ];
            if(!setBoundary) {
              // if a gapBoundary: use the G0-Set for the set matching the associated body!
              if(m_bodyToSetTable[bodyIdNghbr] == set) usedSet = 0;
            } else {
              // if a gap- and setBoundary-Node:
              // the node-values must use the G0-Set for any of the two bodies!
              // if(m_bodyToSetTable[bodyIdNghbr] == set || m_bodyToSetTable[bodyId] == set ) {
              usedSet = 0;
              //}
            }
          }

          phi += scalarField[IDX_LSSETMB(nghbrIds[nghbrNode], usedSet)];
        }

        phi /= noNeighborsPerNode;

        for(MInt nghbrNode = 0; nghbrNode < noNeighborsPerNode; nghbrNode++) {
          const MInt nghbrCand = candidateIds[nghbrIds[nghbrNode]];
          if(nghbrCand < 0) continue;
          candidates[nghbrCand].nodeValueSet[nghbrNodes[nghbrNode]] = true;
          candidates[nghbrCand].nodalValues[set][nghbrNodes[nghbrNode]] = phi;
        }
      }

#ifdef CutCell_DEBUG
      for(MInt nghbrNode = 0; nghbrNode < noNeighborsPerNode; nghbrNode++) {
        if(gapBoundary) {
          const MInt candidateId = candidateIds[nghbrIds[nghbrNode]];
          if(candidateId < 0) continue;
          const MInt cellId2 = candidates[candidateId].cellId;
          ASSERT(cellId2 == nghbrIds[nghbrNode], "");
          const MInt set = m_bodyToSetTable[bodyIdField[IDX_LSSETMB(cellId2, 0)]];
          if(abs(candidates[candidateId].nodalValues[0][nghbrNodes[nghbrNode]]
                 - candidates[candidateId].nodalValues[set][nghbrNodes[nghbrNode]])
             > 0.0000008)
            cerr << "ERROR in Node-LS-Values " << cellId2 << " globalId " << grid().tree().globalId(cellId2) << " "
                 << grid().tree().globalId(cellId) << " "
                 << candidateId /*<< " isHalo " << a_isHalo(cellId2) << " " << a_isHalo(cellId) */
                 << " Gap " << candidates[cnd].isGapCell << " " << candidates[candidateId].isGapCell << " BodyId "
                 << candidates[cnd].associatedBodyIds[0] << " " << candidates[candidateId].associatedBodyIds[0] << " "
                 << nghbrNode << endl;
        }
      }
#endif
    }
  }

  if(m_bndryLvlJumps && !m_cutLvlJumpCandidates.empty()) {
    correctNodalValuesAtLevelJump(candidates, candidateIds);
  }

  // remove bndryJumpCellParents
  for(MInt cnd = noCandidates; cnd < (signed)candidates.size(); cnd++) {
    const MInt parentId = candidates[cnd].cellId;
    ASSERT(candidates[cnd].isbndryLvlJumpParent, "");
    candidateIds[parentId] = -1;
  }

  candidates.resize(noCandidates);
}


/**
 * \brief halo cell - window cell exchange of nodal scalar field data
 * \author Lennart Schneiders, Update Tim Wegmann
 */
template <MInt nDim_>
void GeometryIntersection<nDim_>::exchangeNodalValues(const MInt** maxLevelWindowCells,
                                                      const MInt* noMaxLevelWindowCells,
                                                      const MInt** maxLevelHaloCells,
                                                      std::vector<CutCandidate<nDim>>& candidates,
                                                      MInt* candidateIds) {
  TRACE();

  auto* mpi_send_req = new MPI_Request[grid().noNeighborDomains()];
  auto* mpi_recv_req = new MPI_Request[grid().noNeighborDomains()];

  MIntScratchSpace sendBufferSize(grid().noNeighborDomains(), AT_, "sendBufferSize");
  MIntScratchSpace receiveBufferSize(grid().noNeighborDomains(), AT_, "receiveBufferSize");

  MIntScratchSpace sendBuffersInt(grid().noNeighborDomains(), AT_, "sendBuffersInt");
  MIntScratchSpace receiveBuffersInt(grid().noNeighborDomains(), AT_, "receiveBuffersInt");

  for(MInt i = 0; i < grid().noNeighborDomains(); i++) {
    mpi_send_req[i] = MPI_REQUEST_NULL;
    mpi_recv_req[i] = MPI_REQUEST_NULL;
  }

  // 0. gather information about the number of values to be exchanged:
  MInt sendBufferCounter;
  MInt sendBufferCounterOverall = 0;
  for(MInt i = 0; i < grid().noNeighborDomains(); i++) {
    sendBufferCounter = 0;
    sendBufferCounter++;
    for(MInt j = 0; j < noMaxLevelWindowCells[i]; j++) {
      const MInt cellId = maxLevelWindowCells[i][j];
      const MInt cndId = candidateIds[cellId];
      if(cndId < 0) continue;
      for(MInt node = 0; node < m_noCorners; node++) {
        ASSERT(candidates[cndId].nodeValueSet[node], "");
      }

      //+check
      sendBufferCounter++;

      //+nodalvalues
      sendBufferCounter += m_noLevelSetsUsedForMb * m_noCorners;

      //+isGapCell
      sendBufferCounter++;

      //+associatedBodyIds
      sendBufferCounter += m_noLevelSetsUsedForMb;
    }
    sendBufferSize[i] = sendBufferCounter;
    sendBuffersInt[i] = sendBufferCounter;
    sendBufferCounterOverall += sendBufferCounter;
  }

  // 0.a. exchange number data solvers to be exchanged:
  for(MInt i = 0; i < grid().noNeighborDomains(); i++) {
    MPI_Irecv(&receiveBuffersInt[i], 1, MPI_INT, grid().neighborDomain(i), 2, grid().mpiComm(), &mpi_recv_req[i], AT_,
              "receiveBuffers[i]");
  }
  for(MInt i = 0; i < grid().noNeighborDomains(); i++) {
    MPI_Isend(&sendBuffersInt[i], 1, MPI_INT, grid().neighborDomain(i), 2, grid().mpiComm(), &mpi_send_req[i], AT_,
              "sendBuffersInt[i]");
  }
  MPI_Waitall(grid().noNeighborDomains(), mpi_recv_req, MPI_STATUSES_IGNORE, AT_);
  MPI_Waitall(grid().noNeighborDomains(), mpi_send_req, MPI_STATUSES_IGNORE, AT_);

  // 0.b setup real exchange framework:
  MInt receiveBufferCounterOverall = 0;
  for(MInt i = 0; i < grid().noNeighborDomains(); i++) {
    receiveBufferCounterOverall += receiveBuffersInt[i];
    receiveBufferSize[i] = receiveBuffersInt[i];
  }

  MFloatScratchSpace sendBuffersOverall(sendBufferCounterOverall, AT_, "sendBuffersOverall");
  MFloatScratchSpace receiveBuffersOverall(receiveBufferCounterOverall, AT_, "receiveBuffersOverall");
  MFloatPointerScratchSpace sendBuffers(grid().noNeighborDomains(), AT_, "sendBuffers");
  MFloatPointerScratchSpace receiveBuffers(grid().noNeighborDomains(), AT_, "receiveBuffers");
  MInt counterSend = 0;
  MInt counterReceive = 0;
  for(MInt i = 0; i < grid().noNeighborDomains(); i++) {
    sendBuffers.p[i] = &sendBuffersOverall.p[counterSend];
    counterSend += sendBufferSize[i];
    receiveBuffers.p[i] = &receiveBuffersOverall.p[counterReceive];
    counterReceive += receiveBufferSize[i];
  }

  // sicherheitshalber zuruecksetzen
  for(MInt i = 0; i < grid().noNeighborDomains(); i++) {
    mpi_send_req[i] = MPI_REQUEST_NULL;
    mpi_recv_req[i] = MPI_REQUEST_NULL;
  }
  // 1. gather relevant information
  for(MInt i = 0; i < grid().noNeighborDomains(); i++) {
    sendBufferCounter = 0;
    sendBuffers[i][0] = (MFloat)(-1);
    sendBufferCounter++;
    for(MInt j = 0; j < noMaxLevelWindowCells[i]; j++) {
      const MInt cellId = maxLevelWindowCells[i][j];
      const MInt cndId = candidateIds[cellId];
      if(cndId < 0) continue;
      // check
      sendBuffers[i][sendBufferCounter] = (MFloat)j;
      sendBufferCounter++;
      // nodalValues
      for(MInt node = 0; node < m_noCorners; node++) {
        ASSERT(candidates[cndId].nodeValueSet[node], "");
        for(MInt set = 0; set < m_noLevelSetsUsedForMb; set++) {
          sendBuffers[i][sendBufferCounter] = candidates[cndId].nodalValues[set][node];
          sendBufferCounter++;
        }
      }
      // isGapCell
      sendBuffers[i][sendBufferCounter] = (MFloat)candidates[cndId].isGapCell;
      sendBufferCounter++;
      // associatedBodyIds
      for(MInt set = 0; set < m_noLevelSetsUsedForMb; set++) {
        sendBuffers[i][sendBufferCounter] = (MFloat)candidates[cndId].associatedBodyIds[set];
        sendBufferCounter++;
      }
    }
    sendBufferSize[i] = sendBufferCounter;
    sendBuffers[i][0] = (MFloat)(sendBufferCounter);
  }


  // 2. receive data
  for(MInt i = 0; i < grid().noNeighborDomains(); i++) {
    MInt bufSize = receiveBufferSize[i];
    MPI_Irecv(receiveBuffers[i], bufSize, MPI_DOUBLE, grid().neighborDomain(i), 2, grid().mpiComm(), &mpi_recv_req[i],
              AT_, "receiveBuffers[i]");
  }

  // 3. send data
  for(MInt i = 0; i < grid().noNeighborDomains(); i++) {
    MInt bufSize = sendBufferSize[i];
    MPI_Isend(sendBuffers[i], bufSize, MPI_DOUBLE, grid().neighborDomain(i), 2, grid().mpiComm(), &mpi_send_req[i], AT_,
              "sendBuffers[i]");
  }
  MPI_Waitall(grid().noNeighborDomains(), mpi_recv_req, MPI_STATUSES_IGNORE, AT_);
  MPI_Waitall(grid().noNeighborDomains(), mpi_send_req, MPI_STATUSES_IGNORE, AT_);


  // 4. store recieved data
  MInt receiveBufferCounter = 0;
  MInt j = -1;
  for(MInt i = 0; i < grid().noNeighborDomains(); i++) {
    receiveBufferCounter = 0;
    if(receiveBufferSize[i] != receiveBuffersInt[i]) {
      mTerm(1, AT_, " receiveBufferSize doesn't match expected size from previous communication! ");
    }
    if(receiveBufferSize[i] == 0) {
      m_log << "Warning: empty message from rank " << grid().neighborDomain(i) << endl;
    }
    receiveBufferCounter++;
    while(receiveBufferCounter < receiveBufferSize[i]) {
      j = (MInt)receiveBuffers[i][receiveBufferCounter];
      receiveBufferCounter++;
      const MInt cellId = maxLevelHaloCells[i][j];
      MBool skip = false;

      if(cellId > grid().tree().size()) skip = true;
      if(!skip) {
        // add halo-cutCandidate
        const MInt candId = candidates.size();
        candidates.emplace_back();
        candidates[candId].cellId = cellId;

        // add infomation from computeNodalvalues:
        for(MInt node = 0; node < m_noCorners; node++) {
          candidates[candId].nodeValueSet[node] = true;
          for(MInt set = 0; set < m_noLevelSetsUsedForMb; set++) {
            candidates[candId].nodalValues[set][node] = receiveBuffers[i][receiveBufferCounter];
            receiveBufferCounter++;
          }
        }

        candidates[candId].isGapCell = (MInt)receiveBuffers[i][receiveBufferCounter];
        receiveBufferCounter++;

        for(MInt set = 0; set < m_noLevelSetsUsedForMb; set++) {
          candidates[candId].associatedBodyIds[set] = (MInt)receiveBuffers[i][receiveBufferCounter];
          receiveBufferCounter++;
        }

        candidateIds[cellId] = candId;


      } else {
        receiveBufferCounter += m_noLevelSetsUsedForMb * m_noCorners;
      }
    }
  }

  delete[] mpi_send_req;
  delete[] mpi_recv_req;
}

/**
 * \brief set cutCellData based on cutCellCandidates
 * note: not all candidates will become a cutCell!
 * \author Tim Wegmann
 */

template <>
void GeometryIntersection<3>::fillCutCellData(std::vector<CutCandidate<nDim>>& candidates,
                                              std::vector<CutCell<nDim>>& cutCellData,
                                              std::vector<MInt>
                                                  cutCellIdMapping) {
  TRACE();

  const MInt faceCornerMapping[6][4] = {{2, 6, 4, 0}, {1, 5, 7, 3}, {0, 4, 5, 1},
                                        {3, 7, 6, 2}, {2, 3, 1, 0}, {6, 7, 5, 4}};
  using CC = CutCell<nDim>;

  ASSERT(m_noLevelSetsUsedForMb <= CC::maxNoSets, "");

  for(MInt cndt = 0; cndt < (signed)candidates.size(); cndt++) {
    if(candidates[cndt].noCutPoints < nDim) continue;

    const MInt cutCell = cutCellData.size();
    const MInt cellId = candidates[cndt].cellId;

    cutCellIdMapping[cellId] = cndt;

    cutCellData.emplace_back();
    cutCellData[cutCell].cellId = cellId;
    cutCellData[cutCell].isGapCell = candidates[cndt].isGapCell;
    ASSERT(candidates[cndt].noCutPoints <= CC::maxNoCutPoints, "");
    cutCellData[cutCell].noCutPoints = candidates[cndt].noCutPoints;
    for(MInt cp = 0; cp < candidates[cndt].noCutPoints; cp++) {
      for(MInt i = 0; i < nDim; i++) {
        cutCellData[cutCell].cutPoints[cp][i] = candidates[cndt].cutPoints[cp][i];
      }
      cutCellData[cutCell].cutBodyIds[cp] = candidates[cndt].cutBodyIds[cp];
      cutCellData[cutCell].cutEdges[cp] = candidates[cndt].cutEdges[cp];
    }

    for(MInt set = 0; set < m_noLevelSetsUsedForMb; set++) {
      cutCellData[cutCell].associatedBodyIds[set] = candidates[cndt].associatedBodyIds[set];

      for(MInt k = 0; k < CC::noFaces; k++) {
        MFloat phi = F0;
        for(MInt j = 0; j < 4; j++) {
          phi += F1B4 * candidates[cndt].nodalValues[set][faceCornerMapping[k][j]];
        }
        cutCellData[cutCell].faceCentroidIsInsideGeometry[set][k] = (phi < F0);
      }
    }

    MInt startSet = 0;
    MInt endSet = 1;
    if(m_complexBoundary && m_noLevelSetsUsedForMb > 1 && (!cutCellData[cutCell].isGapCell)) {
      startSet = 1;
      endSet = m_noLevelSetsUsedForMb;
    }

    for(MInt node = 0; node < m_noCorners; node++) {
      ASSERT(candidates[cndt].nodeValueSet[node], "");
      MBool isInside = false;
      for(MInt set = startSet; set < endSet; set++) {
        const MBool isInsideSet = candidates[cndt].nodalValues[set][node] < F0;
        isInside = isInside || isInsideSet;
        cutCellData[cutCell].cornerIsInsideGeometry[set][node] = isInsideSet;
      }
      cutCellData[cutCell].cornerIsInsideGeometry[0][node] = isInside;
    }
  }
}


/**
 * \brief computes cut points where candidate intersects with the geometry
 * note: can not handle bndryRefinement jumps yet!
 * \author Daniel Hartmann, Claudia Guenther, Sohel Herff, Tim Wegmann
 * \date 2006, 08/2013, 2016, 2019
 */

template <MInt nDim_>
void GeometryIntersection<nDim_>::computeCutPointsFromSTL(std::vector<CutCandidate<nDim>>& candidates) {
  TRACE();

  if(candidates.empty()) return;

  // NOTE: as this is stl-dependand, the unscaled regular cutCells need to be used!
  m_scaledCutCell = false;

  const MInt DOFStencil[12] = {1, 1, 0, 0, 1, 1, 0, 0, 2, 2, 2, 2};

  const MFloat signStencil[8][3] = {{-F1, -F1, -F1}, {F1, -F1, -F1}, {-F1, F1, -F1}, {F1, F1, -F1},
                                    {-F1, -F1, F1},  {F1, -F1, F1},  {-F1, F1, F1},  {F1, F1, F1}};

  const MInt nodeStencil[2][12] = {{0, 1, 0, 2, 4, 5, 4, 6, 0, 1, 2, 3}, {2, 3, 1, 3, 6, 7, 5, 7, 4, 5, 6, 7}};

  const MFloat epsilon = 0.00000000001;
  const MInt maxNoCutPointsPerEdge = 10;
  std::array<MFloat, nDim> a;
  std::array<MFloat, nDim> b;
  std::array<MFloat, nDim> c;
  std::array<MFloat, nDim> d;
  std::array<MFloat, nDim> e;
  std::array<MFloat, nDim> pP;

  MFloat target[nDim * 2];
  MFloat corners[m_noCorners][nDim];

  MFloatScratchSpace cutPointsEdge(maxNoCutPointsPerEdge, nDim, AT_, "cutPointsEdge");


  //--- end of initialization
  const MInt m_maxNoBndryCndIds = 100;

  MInt noBndryCndIds = 0;
  MInt bndryCndIds[m_maxNoBndryCndIds]{};


  // first loop over all bndry cells
  for(MInt cndt = 0; cndt < (signed)candidates.size(); cndt++) {
    ASSERT(candidates[cndt].noCutPoints == 0, "");
    const MInt cellId = candidates[cndt].cellId;
    // assiciatedBodyId[0] corresponds with bndryCndId and is initialised!
    candidates[cndt].associatedBodyIds[0] = std::numeric_limits<MInt>::max();

    ASSERT(cellId < grid().tree().size(), "");
    // if(!grid().tree().isLeafCell(cellId)) continue;
    if(grid().tree().noChildren(cellId) > 0) continue;
    const MFloat cellLength = a_cellLengthAtCell(cellId);
    const MFloat cellHalfLength = F1B2 * cellLength;

    const MFloat eps = cellLength / 10000000.0;

    // compute cell corners -> will be used for cut point computation
    for(MInt node = 0; node < m_noCorners; node++) {
      for(MInt i = 0; i < nDim; i++) {
        corners[node][i] = grid().tree().coordinate(cellId, i) + signStencil[node][i] * cellHalfLength;
      }
    }

    // define target around corners of current cell
    for(MInt i = 0; i < nDim; i++) {
      target[i] = grid().tree().coordinate(cellId, i) - cellHalfLength * 1.05;
      target[i + nDim] = grid().tree().coordinate(cellId, i) + cellHalfLength * 1.05;
    }

    // get all triangles in currentCell (target)
    std::vector<MInt> nodeList;
    if(m_gridCutTest == "SAT") {
      geometry().getIntersectionElements(target, nodeList, cellHalfLength, &a_coordinate(cellId, 0));
    } else {
      geometry().getIntersectionElements(target, nodeList);
    }


    // loop over all edges
    for(MInt edge = 0; edge < m_noEdges; edge++) {
      const MInt edgeDOF = DOFStencil[edge]; // edge's degree of freedom

      ASSERT(candidates[cndt].noCutPointsOnEdge[edge] == 0, "");

      MBool edgeCut = false;

      // compute end points of edge
      for(MInt i = 0; i < nDim; i++) {
        d[i] = corners[nodeStencil[0][edge]][i];
        e[i] = corners[nodeStencil[1][edge]][i];
      }

      // check for edge-triangle intersection
      for(MInt n = 0; n < (signed)nodeList.size(); n++) {
        MBool cutsEdge = false;

        IF_CONSTEXPR(nDim == 3) {
          const MInt spaceId = (edgeDOF + 1) % nDim;
          const MInt spaceId1 = (edgeDOF + 2) % nDim;
          const MInt spaceId2 = edgeDOF;
          // find out if element cuts edge; if yes, compute cut point -> 3D part
          MFloat p = F0;
          MFloat q = F0;
          for(MInt k = 0; k < nDim; k++) {
            a[k] = m_geometry->elements[nodeList[n]].m_vertices[0][k];
            b[k] = m_geometry->elements[nodeList[n]].m_vertices[1][k];
            c[k] = m_geometry->elements[nodeList[n]].m_vertices[2][k];
          }
          if(approx(a[spaceId1], b[spaceId1], MFloatEps) && !approx(a[spaceId1], c[spaceId1], MFloatEps)) {
            // aspaceId != bspaceId, otherwise a and b would be the same point
            q = (d[spaceId1] - a[spaceId1]) / (c[spaceId1] - a[spaceId1]);
            p = (d[spaceId] - a[spaceId] - q * (c[spaceId] - a[spaceId])) / (b[spaceId] - a[spaceId]);
          } else {
            if(!approx(a[spaceId1], b[spaceId1], MFloatEps) && approx(a[spaceId1], c[spaceId1], MFloatEps)) {
              // aspaceId != cspaceId, otherwise a and c would be the same point
              p = (d[spaceId1] - a[spaceId1]) / (b[spaceId1] - a[spaceId1]);
              q = (d[spaceId] - a[spaceId] - p * (b[spaceId] - a[spaceId])) / (c[spaceId] - a[spaceId]);
            } else {
              if(approx(a[spaceId], b[spaceId], MFloatEps) && !approx(a[spaceId], c[spaceId], MFloatEps)) {
                // aspaceId1 != bspaceId1, otherwise a and b would be the same point
                q = (d[spaceId] - a[spaceId]) / (c[spaceId] - a[spaceId]);
                p = (d[spaceId1] - a[spaceId1] - q * (c[spaceId1] - a[spaceId1])) / (b[spaceId1] - a[spaceId1]);
              } else {
                if(!approx(a[spaceId], b[spaceId], MFloatEps) && approx(a[spaceId], c[spaceId], MFloatEps)) {
                  // aspaceId1 != cspaceId1, otherwise a and c would be the same point
                  p = (d[spaceId] - a[spaceId]) / (b[spaceId] - a[spaceId]);
                  q = (d[spaceId1] - a[spaceId1] - p * (b[spaceId1] - a[spaceId1])) / (c[spaceId1] - a[spaceId1]);
                } else {
                  // aspaceId1 != bspaceId1 && aspaceId1 != cspaceId1 && aspaceId != bspaceId && aspaceId != cspaceId
                  q = ((d[spaceId1] - a[spaceId1]) * (b[spaceId] - a[spaceId])
                       - (b[spaceId1] - a[spaceId1]) * (d[spaceId] - a[spaceId]))
                      / ((c[spaceId1] - a[spaceId1]) * (b[spaceId] - a[spaceId])
                         - (b[spaceId1] - a[spaceId1]) * (c[spaceId] - a[spaceId]));
                  p = (d[spaceId] - a[spaceId] - q * (c[spaceId] - a[spaceId])) / (b[spaceId] - a[spaceId]);
                }
              }
            }
          }


          if(p * q >= 0 || p * q < 0) {
            // compute s
            MFloat gamma = a[spaceId2] + p * (b[spaceId2] - a[spaceId2]) + q * (c[spaceId2] - a[spaceId2]);
            MFloat s = (gamma - d[spaceId2]) / (e[spaceId2] - d[spaceId2]);

            if(s < -epsilon || s > F1 + epsilon || p < -epsilon || q < -epsilon || (p + q) > F1 + epsilon) {
            } else {
              cutsEdge = true;
              // cut point pP
              if(candidates[cndt].noCutPointsOnEdge[edge] < maxNoCutPointsPerEdge) {
                for(MInt k = 0; k < nDim; k++) {
                  pP[k] = d[k] + s * (e[k] - d[k]);
                  cutPointsEdge(candidates[cndt].noCutPointsOnEdge[edge], k) = pP[k];
                }
              } else {
                mTerm(1, AT_, " Too many cut points on edge...");
              }
            }
          }
        }
        else { // 2D code
          if(geometry().edgeTriangleIntersection(geometry().elements[nodeList[n]].m_vertices[0],
                                                 geometry().elements[nodeList[n]].m_vertices[1], nullptr, d.data(),
                                                 e.data())) {
            for(MInt k = 0; k < nDim; k++) {
              a[k] = geometry().elements[nodeList[n]].m_vertices[0][k];
              b[k] = geometry().elements[nodeList[n]].m_vertices[1][k];
            }

            MFloat gamma = (b[0] - a[0]) * (d[1] - e[1]) - (d[0] - e[0]) * (b[1] - a[1]);
            if(ABS(gamma) < 0.0000000000001) continue;
            MFloat s1 = ((d[0] - e[0]) * (a[1] - d[1]) - (d[1] - e[1]) * (a[0] - d[0])) / gamma;
            MFloat s2 = ((b[0] - a[0]) * (d[1] - a[1]) - (d[0] - a[0]) * (b[1] - a[1])) / gamma;

            cutsEdge = true;

            // cut point pP
            if(candidates[cndt].noCutPointsOnEdge[edge] < maxNoCutPointsPerEdge) {
              for(MInt k = 0; k < nDim; k++) {
                if(s1 * s1 < s2 * s2)
                  pP[k] = d[k] + s2 * (e[k] - d[k]);
                else
                  pP[k] = a[k] + s1 * (b[k] - a[k]);
                cutPointsEdge(candidates[cndt].noCutPointsOnEdge[edge], k) = pP[k];
              }
            } else {
              mTerm(1, AT_, " Too many cut points on edge...");
            }
          }
        }

        if(cutsEdge) {
          // store cut point
          if(candidates[cndt].noCutPoints < m_noEdges) {
            // set boundary condition (the one with the lowest Id)!
            const MInt bndCndId = geometry().elements[nodeList[n]].m_bndCndId;
            if(bndCndId < candidates[cndt].associatedBodyIds[0]) {
              candidates[cndt].associatedBodyIds[0] = bndCndId;
              MBool newBndryCnd = true;
              for(MInt j = 0; j < noBndryCndIds; j++) {
                if(bndryCndIds[j] == bndCndId) {
                  newBndryCnd = false;
                  break;
                }
              }
              if(newBndryCnd) {
                ASSERT(noBndryCndIds < m_maxNoBndryCndIds, "Too many different boundary conditions...");
                bndryCndIds[noBndryCndIds] = bndCndId;
                noBndryCndIds++;
              }
            }
            if(edgeCut) {
              // if this edge has already one cut point, check if the new cut point is a different one
              // if the cut point is identical to onother one, special treatment is needed

              MBool newCutPoint = true;
              for(MInt i = 0; i < candidates[cndt].noCutPointsOnEdge[edge]; i++) {
                MBool equal = (ABS(pP[0] - cutPointsEdge(i, 0)) < eps && ABS(pP[1] - cutPointsEdge(i, 1)) < eps);
                IF_CONSTEXPR(nDim == 3) equal = (equal && (ABS(pP[2] - cutPointsEdge(i, 2)) < eps));
                if(equal) {
                  newCutPoint = false;
                  break;
                }
              }

              if(newCutPoint) {
                candidates[cndt].noCutPoints--;
                m_log << "removing two cut points, cell " << cellId << " edge " << edge << endl;
                m_log << a_coordinate(cellId, 0) << " " << a_coordinate(cellId, 1) << " ";
                IF_CONSTEXPR(nDim == 3) m_log << a_coordinate(cellId, 2) << endl;
                edgeCut = false;
                candidates[cndt].noCutPointsOnEdge[edge]++;
              } else {
                // coinciding cut points
                // ignoring one cut point
              }
            } else {
              // check, if cut point is really new!
              MBool newCutPoint = true;
              for(MInt i = 0; i < candidates[cndt].noCutPointsOnEdge[edge]; i++) {
                MBool equal = (ABS(pP[0] - cutPointsEdge(i, 0)) <= eps && ABS(pP[1] - cutPointsEdge(i, 1)) <= eps);
                IF_CONSTEXPR(nDim == 3) equal = (equal && (ABS(pP[2] - cutPointsEdge(i, 2)) <= eps));
                if(equal) {
                  newCutPoint = false;
                }
              }
              if(newCutPoint) {
                edgeCut = true;
                // store cut points in the data structure
                const MInt cutPointNo = candidates[cndt].noCutPoints;
                candidates[cndt].cutEdges[cutPointNo] = edge;
                candidates[cndt].cutBodyIds[cutPointNo] = 0;
                for(MInt i = 0; i < nDim; i++) {
                  candidates[cndt].cutPoints[cutPointNo][i] = pP[i];
                }
                candidates[cndt].noCutPoints++;
                candidates[cndt].noCutPointsOnEdge[edge]++;
              }
            }
          } else {
            cerr << "** Warning fvbndrycndxd: " << endl;
            cerr << "cell " << cellId << endl;
            cerr << " -> Boundary candidate " << cndt << " has more than " << m_noEdges << " cut points" << endl;
            cerr << " -> Additional cut points are neglected" << endl;
            cerr << cellId << " " << a_coordinate(cellId, 0) << " " << a_coordinate(cellId, 1) << " ";
            IF_CONSTEXPR(nDim == 3) cerr << a_coordinate(cellId, 2) << " ";
            cerr << grid().tree().level(cellId) << endl;
            mTerm(1, AT_, "Boundary cell has more than m_noEdges cut points");
          }
        }
      } // loop over entity-nodes
    }   // loop over all edges
  }     // loop over all cut candidates

  m_scaledCutCell = true;
}

/**
 * \brief corrects the nodal values at bndry level jumps for ls-based geometries
 * \author Tim Wegmann
 * \date 2020
 */
template <MInt nDim_>
void GeometryIntersection<nDim_>::correctNodalValuesAtLevelJump(std::vector<CutCandidate<nDim>>& candidates,
                                                                const MInt* candidateIds) {
  // NOTE: it is not possible to copy node values from a neighbor!
  //      the neighbor might be a halo-cell and thus its node value not set on this rank before
  //      the exchange, only its cell-centered value is available!
  //       => only copy values from the parent-node!

#ifdef CutCell_DEBUG
  const MInt debugTimeStep = -2;
  const MInt debugGlobalId = 216221;

  for(MInt cellId = 0; cellId < grid().tree().size(); cellId++) {
    if(globalTimeStep > debugTimeStep && grid().tree().globalId(cellId) == debugGlobalId) {
      const MInt cand = candidateIds[cellId];
      if(cand < 0) break;
      cerr << " Original node-values: ";
      for(MInt node = 0; node < m_noCorners; node++) {
        cerr << candidates[cand].nodalValues[0][node] << " ";
      }
      cerr << endl;
    }
  }
#else
  std::ignore = candidateIds[0];
#endif

  ASSERT(nDim == 3 && nDim_ == 3, "Intended otherwise!");

  // returns node which matches for child and parent for the specified childPosition
  //{7,   6,   5,   4,   3,   2,  1,  0}
  const constexpr MInt childParentNode[8] = {0, 1, 2, 3, 4, 5, 6, 7};

  const constexpr MInt faceCornerCode[6][4] = {{0, 2, 6, 4}, {3, 1, 5, 7}, {1, 0, 4, 5},
                                               {2, 3, 7, 6}, {0, 1, 3, 2}, {5, 4, 6, 7}};

  // number of nodes which need to be corrected for each jump-direction!
  const constexpr MInt noCorrectNodes = nDim == 2 ? 2 : 4;

  // returns the 2 corners of a edge
  const constexpr MInt edge2Corner[12][2] = {{0, 2}, {1, 3}, {0, 1}, {2, 3}, {4, 6}, {5, 7},
                                             {4, 5}, {6, 7}, {0, 4}, {1, 5}, {2, 6}, {3, 7}};

  // returns the 3 edges of a corner  0   1    2    3    4    5    6     7
  const constexpr MInt corner2Edge[3][8] = {
      {0, 1, 0, 1, 4, 5, 4, 5}, {2, 2, 3, 3, 6, 6, 7, 7}, {8, 9, 10, 11, 8, 9, 10, 11}};

  // recompute node values for nodes at the lvlJump interface
  for(MInt id = 0; id < (MInt)m_cutLvlJumpCandidates.size(); id++) {
    const MInt parentCand = m_cutLvlJumpCandidates[id].parentCandId;
    const MInt candId = m_cutLvlJumpCandidates[id].candId;
    const MInt cellId = candidates[candId].cellId;
    ASSERT(grid().tree().parent(cellId) == candidates[parentCand].cellId, "");
    ASSERT(candidates[parentCand].isbndryLvlJumpParent, "");

    for(MInt node = 0; node < m_noCorners; node++) {
      ASSERT(candidates[parentCand].nodeValueSet[node], "");
    }

    const MInt noJumps = m_cutLvlJumpCandidates[id].noJumps;
    ASSERT(noJumps <= 7 && noJumps > 0, to_string(noJumps));

    const MInt childPos = m_cutLvlJumpCandidates[id].childId;

#ifdef CutCell_DEBUG
    if(globalTimeStep > debugTimeStep && grid().tree().globalId(cellId) == debugGlobalId) {
      cerr << "Cell-Info: noJumps" << noJumps << endl;
      for(MInt j = 0; j < noJumps; j++) {
        cerr << "dir " << m_cutLvlJumpCandidates[id].dirs[j] << " type " << m_cutLvlJumpCandidates[id].neighborType[j];
      }
      cerr << endl;
    }
#endif

    // copy values from matching parent node (node == parentNode)
    const MInt parentNode = childParentNode[childPos];
    for(MInt set = 0; set < m_noLevelSetsUsedForMb; set++) {
      candidates[candId].nodalValues[set][parentNode] = candidates[parentCand].nodalValues[set][parentNode];
    }

    const MInt parentEdgeIds[3] = {corner2Edge[0][parentNode], corner2Edge[1][parentNode], corner2Edge[2][parentNode]};

#ifdef CutCell_DEBUG
    if(globalTimeStep > debugTimeStep && grid().tree().globalId(cellId) == debugGlobalId) {
      cerr << " Child-Pos " << childPos << " setting (1) node " << parentNode << " to "
           << candidates[parentCand].nodalValues[0][parentNode] << endl;
    }
#endif

    // loop over all jumps and correct the remaining nodal values in each jump-dir
    for(MInt j = 0; j < noJumps; j++) {
      const MInt jumpDir = m_cutLvlJumpCandidates[id].dirs[j];
      const MInt type = m_cutLvlJumpCandidates[id].neighborType[j];

#ifdef CutCell_DEBUG
      if(globalTimeStep > debugTimeStep && grid().tree().globalId(cellId) == debugGlobalId) {
        cerr << childPos << " " << noJumps << type << endl;
      }
#endif

      // for 3D-diagonal neighbors nothing else is necessary
      if(type == 2) {
        continue;
      }

      for(MInt i = 0; i < noCorrectNodes; i++) {
        const MInt node = faceCornerCode[jumpDir][i];
        // already correct above(directly copied!)
        if(node == parentNode) continue;

        MBool isDiagonalNode = false;
        // only correct if the node is also on the face of the second-diagonal direction!
        if(type == 1) {
          const MInt secondDir = m_cutLvlJumpCandidates[id].diagonalDirs[j];
          for(MInt l = 0; l < noCorrectNodes; l++) {
            if(node == faceCornerCode[secondDir][l]) {
              isDiagonalNode = true;
              break;
            }
          }
          if(!isDiagonalNode) continue;
        }

        MInt edgeInterpolation = -1;
        for(MInt k = 0; k < 3; k++) {
          MInt parentEdgeCorner = edge2Corner[parentEdgeIds[k]][0];
          if(parentEdgeCorner == parentNode) {
            parentEdgeCorner = edge2Corner[parentEdgeIds[k]][1];
          }
          if(parentEdgeCorner == node) {
            edgeInterpolation = k;
            break;
          }
        }

        if(edgeInterpolation > -1) { // edge interpolation for node
          const MInt parentEdgeId = parentEdgeIds[edgeInterpolation];
          const MInt parentNodes[2] = {edge2Corner[parentEdgeId][0], edge2Corner[parentEdgeId][1]};
          for(MInt set = 0; set < m_noLevelSetsUsedForMb; set++) {
            const MFloat phi = 0.5
                               * (candidates[parentCand].nodalValues[set][parentNodes[0]]
                                  + candidates[parentCand].nodalValues[set][parentNodes[1]]);
            candidates[candId].nodalValues[set][node] = phi;
          }

#ifdef CutCell_DEBUG
          if(globalTimeStep > debugTimeStep && grid().tree().globalId(cellId) == debugGlobalId) {
            cerr << " (2 ) Setting node " << node << " to " << candidates[candId].nodalValues[0][node] << endl;
          }
#endif

        } else { // corner interpolation for node (node is on cell-face)
          ASSERT(type == 0, "");
          const MInt parentNodes[4] = {faceCornerCode[jumpDir][0], faceCornerCode[jumpDir][1],
                                       faceCornerCode[jumpDir][2], faceCornerCode[jumpDir][3]};
          for(MInt set = 0; set < m_noLevelSetsUsedForMb; set++) {
            const MFloat phi = 0.25
                               * (candidates[parentCand].nodalValues[set][parentNodes[0]]
                                  + candidates[parentCand].nodalValues[set][parentNodes[1]]
                                  + candidates[parentCand].nodalValues[set][parentNodes[2]]
                                  + candidates[parentCand].nodalValues[set][parentNodes[3]]);
            candidates[candId].nodalValues[set][node] = phi;
          }

#ifdef CutCell_DEBUG
          if(globalTimeStep > debugTimeStep && grid().tree().globalId(cellId) == debugGlobalId) {
            cerr << " (3 ) Setting node " << node << " to " << candidates[candId].nodalValues[0][node] << endl;
          }
#endif
        }
      }
    }
  }

#ifdef CutCell_DEBUG
  const MInt reverseNode[3][8] = {{1, 0, 3, 2, 5, 4, 7, 6}, {2, 3, 0, 1, 6, 7, 4, 5}, {4, 5, 6, 7, 0, 1, 2, 3}};


  // check nodal values:
  std::set<std::pair<MInt, MInt>> nghbrSet;

  MInt nghbrIds[m_noCorners]{};
  MInt nghbrNodes[m_noCorners]{};

  const MInt nodeStencil[3][8] = {{0, 1, 0, 1, 0, 1, 0, 1}, {2, 2, 3, 3, 2, 2, 3, 3}, {4, 4, 4, 4, 5, 5, 5, 5}};


  for(MInt cnd = 0; cnd < (signed)candidates.size(); cnd++) {
    const MInt cellId = candidates[cnd].cellId;
    for(MInt node = 0; node < m_noCorners; node++) {
      nghbrSet.clear();
      nghbrSet.insert(std::make_pair(cellId, node));

      for(MInt i = 0; i < nDim; i++) {
        const MInt firstDir = nodeStencil[i][node];
        const MInt firstNghbrNode = reverseNode[i][node];
        const MInt firstNghbrId = grid().tree().neighbor(cellId, firstDir);
        if(firstNghbrId > -1) {
          nghbrSet.insert(make_pair(firstNghbrId, firstNghbrNode));
          for(MInt j = 0; j < nDim; j++) {
            const MInt secondDir = nodeStencil[j][node];
            if(secondDir == firstDir) continue;
            const MInt secondNghbrNode = reverseNode[j][firstNghbrNode];
            const MInt secondNghbrId = grid().tree().neighbor(firstNghbrId, secondDir);
            if(secondNghbrId > -1) {
              nghbrSet.insert(make_pair(secondNghbrId, secondNghbrNode));
              IF_CONSTEXPR(nDim == 3) {
                for(MInt k = 0; k < nDim; k++) {
                  const MInt thirdDir = nodeStencil[k][node];
                  if(thirdDir == firstDir || thirdDir == secondDir) continue;
                  const MInt thirdNghbrNode = reverseNode[k][secondNghbrNode];
                  const MInt thirdNghbrId = grid().tree().neighbor(secondNghbrId, thirdDir);
                  if(thirdNghbrId > -1) {
                    nghbrSet.insert(make_pair(thirdNghbrId, thirdNghbrNode));
                  }
                }
              }
            }
          }
        }
      }

      // c) reorder list by nodes, and calculate noNeighborsPerNode
      // note: previously determined neighbors might be added multiple-times to the map
      // and overwrite existing and identical information...
      MInt noNeighborsPerNode = 0;
      for(auto it = nghbrSet.begin(); it != nghbrSet.end(); it++) {
        nghbrIds[noNeighborsPerNode] = it->first;
        nghbrNodes[noNeighborsPerNode] = it->second;
        noNeighborsPerNode++;
      }

      for(MInt nghbrNode = 0; nghbrNode < noNeighborsPerNode; nghbrNode++) {
        const MInt nghbrCand = candidateIds[nghbrIds[nghbrNode]];
        if(nghbrCand < 0) continue;
        if(fabs(candidates[nghbrCand].nodalValues[0][nghbrNodes[nghbrNode]] - candidates[cnd].nodalValues[0][node])
           > 0.00000000001) {
          cerr << "Nodal-value missmatch " << cellId << " " << grid().tree().globalId(cellId) << " " << node << " "
               << candidates[cnd].nodalValues[0][node] << " " << nghbrIds[nghbrNode] << " "
               << grid().tree().globalId(nghbrIds[nghbrNode]) << " " << nghbrNodes[nghbrNode] << " "
               << candidates[nghbrCand].nodalValues[0][nghbrNodes[nghbrNode]] << endl;
        }
      }
    }
  }
#endif
}

/**
 * \brief gets all neighbor cellIds and matching nodes for a given cell and node
 * \author Tim Wegmann
 * \date 2020
 */
template <MInt nDim_>
void GeometryIntersection<nDim_>::getNeighborNodes(const MInt cellId, const MInt node, MInt noNeighborsPerNode,
                                                   MInt* nghbrIds, MInt* nghbrNodes) {
  TRACE();

  const MInt nodeStencil[3][8] = {{0, 1, 0, 1, 0, 1, 0, 1}, {2, 2, 3, 3, 2, 2, 3, 3}, {4, 4, 4, 4, 5, 5, 5, 5}};

  const MInt reverseNode[3][8] = {{1, 0, 3, 2, 5, 4, 7, 6}, {2, 3, 0, 1, 6, 7, 4, 5}, {4, 5, 6, 7, 0, 1, 2, 3}};

  std::set<std::pair<MInt, MInt>> nghbrSet;

  // b) Add all neighbors and the corresponding nodes to the list
  nghbrSet.clear();
  nghbrSet.insert(std::make_pair(cellId, node));

  for(MInt i = 0; i < nDim; i++) {
    const MInt firstDir = nodeStencil[i][node];
    const MInt firstNghbrNode = reverseNode[i][node];
    const MInt firstNghbrId = grid().tree().neighbor(cellId, firstDir);
    if(firstNghbrId > -1) {
      nghbrSet.insert(make_pair(firstNghbrId, firstNghbrNode));
      for(MInt j = 0; j < nDim; j++) {
        const MInt secondDir = nodeStencil[j][node];
        if(secondDir == firstDir) continue;
        const MInt secondNghbrNode = reverseNode[j][firstNghbrNode];
        const MInt secondNghbrId = grid().tree().neighbor(firstNghbrId, secondDir);
        if(secondNghbrId > -1) {
          nghbrSet.insert(make_pair(secondNghbrId, secondNghbrNode));
          IF_CONSTEXPR(nDim == 3) {
            for(MInt k = 0; k < nDim; k++) {
              const MInt thirdDir = nodeStencil[k][node];
              if(thirdDir == firstDir || thirdDir == secondDir) continue;
              const MInt thirdNghbrNode = reverseNode[k][secondNghbrNode];
              const MInt thirdNghbrId = grid().tree().neighbor(secondNghbrId, thirdDir);
              if(thirdNghbrId > -1) {
                nghbrSet.insert(make_pair(thirdNghbrId, thirdNghbrNode));
              }
            }
          }
        }
      }
    }
  }
  // c) reorder list by nodes, and calculate noNeighborsPerNode
  // note: previously determined neighbors might be added multiple-times to the map
  // and overwrite existing and identical information...
  noNeighborsPerNode = 0;
  for(auto it = nghbrSet.begin(); it != nghbrSet.end(); it++) {
    nghbrIds[noNeighborsPerNode] = it->first;
    nghbrNodes[noNeighborsPerNode] = it->second;
    noNeighborsPerNode++;
  }
}
template class GeometryIntersection<2>;
template class GeometryIntersection<3>;
