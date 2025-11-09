// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef FVBNDRYCELL_H
#define FVBNDRYCELL_H

#include <vector>
#include "INCLUDE/maiaconstants.h"
#include "INCLUDE/maiatypes.h"
#include "enums.h"

template <MInt nDim, class SysEqn>
class FvBndryCell {
 public:
  static void init(MInt, MInt, MInt, MInt, MInt);
  static MInt m_noSpecies;
  static MInt m_noRansEquations;
  static MInt m_noEdges;
  static MInt m_maxNoSurfaces; // 1 with original boundary formulation, > 1 for complex boundary formulation
  static MInt m_noNghbrs;
  static MInt m_noVariables;

  MBool* m_externalFaces = nullptr;
  MInt m_cellId;
  MInt m_periodicCellId;
  MInt m_linkedCellId;
  MInt* m_associatedSrfc;
  MFloat m_volume;
  MFloat m_gapDistance = std::numeric_limits<MFloat>::max();
  MFloat* m_coordinates = nullptr;
  MFloat* m_masterCoordinates = nullptr;
  std::vector<MInt> m_recNghbrIds;             // neighbors involved in the cell vars reconstruction
  std::vector<MFloat> m_cellVarsRecConst;      // constants to reconstruct the variables at the cell centroid
  std::vector<MFloat> m_cellDerivRecConst;     // constants to reconstruct the derivatives at the cell centroid
  std::vector<MFloat> m_faceVertices;          // face vertices for split face stream
  std::vector<std::vector<MInt>> m_faceStream; // face stream for split faces

  // WMLES
  MBool m_isWMCell = false;

  struct WallModelBCVars {
    MFloat m_wmTauW;
    MFloat m_wmMUEWM;
    MFloat m_wmImgVars[5];
    MFloat m_wmUII;
    MFloat m_wmUTAU;
    MBool m_wmHasImgCell = false;
  };

  struct BodySurface {
    MInt m_bndryCndId = 0;
    MInt m_noCutPoints = 0;
    MInt* m_bodyId = nullptr;
    MFloat m_FJacobian = 0;
    MInt* m_cutEdge = nullptr;
    MFloat m_area = 0;
    MFloat m_centroidDistance = 0;
    MFloat* m_coordinates = nullptr;
    MFloat* m_normalVector = nullptr;
    MFloat* m_normalVectorCentroid = nullptr;
    MFloat* m_planeVector0 = nullptr;
    MFloat* m_planeVector1 = nullptr;
    MFloat** m_cutCoordinates = nullptr;
  };

  struct BodySurfaceVariables {
    MInt m_ghostCellId;
    MInt m_wmSrfcId = -1;
    MInt* m_srfcId = nullptr;
    MFloat* m_imageCoordinates = nullptr;  // needed for complex boundary formulation
    MFloat* m_imageVariables = nullptr;    // needed for complex boundary formulation
    BcType* m_variablesType = nullptr;     // Dirichlet, Neumann or Robin
    MFloat* m_primVars = nullptr;          // the primitive variables in the surface centroid
    MFloat* m_normalDeriv =
        nullptr; // the derivatives of primitive variables in surface-normal direction at the surface centroid
    std::vector<MFloat> m_imagePointRecConst; // constants to reconstruct the variables at the image point centroid
    MFloat m_robinFactor;                     // the preFactor 'a' in Robin-type boundary conditions du/dn + a*u = b
  };

  WallModelBCVars* m_wmBCVars = nullptr;
  BodySurface** m_srfcs = nullptr;
  BodySurfaceVariables** m_srfcVariables = nullptr;
  MInt m_noSrfcs;

  static MInt staticElementSize() {
    return (sizeof(MBool) * (m_noNghbrs) + // external faces
            sizeof(MInt)
                * (m_noNghbrs +                      // associated surfaces
                   2 * m_noEdges * m_maxNoSurfaces + // bodyId
                   2 * m_noEdges * m_maxNoSurfaces + // cutEdge
                   nDim * m_maxNoSurfaces            // srfcId
                   )
            + sizeof(BcType) * m_noVariables * m_maxNoSurfaces + // variablesType
            sizeof(MFloat)
                * (2 * nDim +                             // coordinates, master coordinates
                   nDim * m_maxNoSurfaces +               // coordinates BodySurface
                   m_noVariables * m_maxNoSurfaces +      // imageVariables
                   nDim * m_maxNoSurfaces +               // imageCoordinates
                   m_noVariables * m_maxNoSurfaces +      // primVars
                   m_noVariables * m_maxNoSurfaces +      // normalDeriv
                   4 * nDim * m_maxNoSurfaces +           // normal vector, plane vectors
                   2 * m_noEdges * m_maxNoSurfaces * nDim // cut coordinates
                   )
            + sizeof(MFloat*) * (2 * m_noEdges * m_maxNoSurfaces) + // cut coordinates
            sizeof(BodySurface*) * m_maxNoSurfaces + sizeof(BodySurfaceVariables*) * m_maxNoSurfaces
            +                                                                                      // pointer to struct
            sizeof(BodySurface) * m_maxNoSurfaces + sizeof(BodySurfaceVariables) * m_maxNoSurfaces // structs
    );
  }

  void allocateElements(void*, void*, const MInt);
  void moveElements(void*);
};

#endif
