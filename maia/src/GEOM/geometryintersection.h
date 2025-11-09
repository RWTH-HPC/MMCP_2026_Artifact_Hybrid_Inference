// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef GEOMETRYINTERSECTION_H
#define GEOMETRYINTERSECTION_H

//#define CutCell_DEBUG

#include <cmath>
#include <stack>
#include <type_traits>
#include <typeinfo>
#include "COMM/mpioverride.h"
#include "GRID/cartesiangrid.h"
#include "GRID/cartesiangridproxy.h"
#include "UTIL/timer.h"
#include "geometry.h"
#include "geometrycontext.h"
#include "geometrycontexttypes.h"
#include "geometryintersectionlookuptable.h"

// Additional declaration needed?
template <MInt nDim>
class CutCell;

template <MInt nDim>
class CutCandidate;

template <MInt nDim>
class lvlJumpCandidates;

// \brief Helper function to initialize std::array at compile time.
template <typename T, std::size_t n>
constexpr std::array<T, n> make_array(const T v) {
  std::array<T, n> a{};
  a.fill(v);
  return a;
}

// \brief Helper function to initialize 2D std::array at compile time.
template <typename T, std::size_t n, std::size_t m>
constexpr std::array<std::array<T, m>, n> make_array(const T v) {
  std::array<std::array<T, m>, n> a{};
  for(std::size_t i = 0; i < n; i++) {
    for(std::size_t j = 0; j < m; j++) {
      a[i][j] = v;
    }
  }
  return a;
}

template <MInt nDim_>
class GeometryIntersection {
 private:
  static constexpr MInt nDim = nDim_;
  static constexpr MInt m_noDirs = 2 * nDim;
  static constexpr MInt m_noCorners = IPOW2(nDim);
  static constexpr MInt m_maxNoChilds = IPOW2(nDim);
  static constexpr MInt m_noEdges = 2 * nDim * (nDim - 1);
  // static constexpr MInt m_revDir[6] = {1,0,3,2,5,4};

  using GridProxy = typename maia::grid::Proxy<nDim>;
  using Geom = Geometry<nDim>;

  /// \brief Returns the coordinate of the cell \p cellId for dimension \p dir
  MFloat& a_coordinate(const MInt cellId, const MInt dir) {
    // FIXME labels:GEOM avoid pointer-access to the grid-raw-coordinate!!!
    if(!m_scaledCutCell) {
      const MInt gridCellId = grid().tree().solver2grid(cellId);
      return grid().raw().a_coordinate(gridCellId, dir);
    }

    return m_scaledCoordinate[0];
  }

  /// \brief Returns the cellLength of the cell \p cellId
  MFloat a_cellLengthAtCell(const MInt cellId) const {
    if(!m_scaledCutCell) {
      return grid().cellLengthAtCell(cellId);
    }
    return 2;
  }

  /// \brief Returns the volume of the cell \p cellId
  MFloat a_gridCellVolume(const MInt cellId) const {
    if(!m_scaledCutCell) {
      return grid().gridCellVolume(grid().tree().level(cellId));
    }
    return 8;
  }

 public:
  GeometryIntersection(GridProxy* gridProxy_, Geom* geometry_)
    : m_eps(std::numeric_limits<MFloat>::epsilon()), m_gridProxy(gridProxy_), m_geometry(geometry_) {
    m_bodyFaceJoinMode = 4;
    if(geometryContext().propertyExists("bodyFaceJoinMode", 0)) { // TODO labels:GEOM simplify
      m_bodyFaceJoinMode = *(geometryContext().getProperty("bodyFaceJoinMode", 0)->asInt(0));
    }

    m_bodyFaceJoinCriterion = 0.1;
    if(geometryContext().propertyExists("bodyFaceJoinCriterion", 0)) {
      m_bodyFaceJoinCriterion = *(geometryContext().getProperty("bodyFaceJoinCriterion", 0)->asFloat(0));
    }

    m_gridCutTest = "SAT";
    if(geometryContext().propertyExists("gridCutTest", 0)) { // TODO labels:GEOM simplify
      m_gridCutTest = *(geometryContext().getProperty("gridCutTest", 0)->asString(0));
    }

    m_multiCutCell = false;
    if(geometryContext().propertyExists("multiCutCell", 0)) { // TODO labels:GEOM simplify
      m_multiCutCell = (MBool) * (geometryContext().getProperty("multiCutCell", 0)->asInt(0));
    }

    if(!m_multiCutCell) {
      m_scaledCutCell = false;
    } else {
      m_scaledCutCell = true;
    }

    m_bndryLvlJumps = false;
    if(geometryContext().propertyExists("allowBndryLvlJumps", 0)) { // TODO labels:GEOM simplify
      m_bndryLvlJumps = (MBool) * (geometryContext().getProperty("allowBndryLvlJumps", 0)->asInt(0));
    }

#ifdef CutCell_DEBUG
    mAlloc(m_caseCheckList, 256, "m_caseCheckList", 0, AT_);
#endif

    mAlloc(m_scaledCoordinate, 3, "m_scaledCoordinate", F0, AT_);
  }

  ~GeometryIntersection() {
    returnDebugInfo();
    resetCutCellData();
  };

  const std::vector<CutCell<nDim>>& cutCellData_() { return m_cutCellData; }
  const std::vector<CutCandidate<nDim>>& cutCandidates_() { return m_cutCandidates; }


  // move to private later --->
  //  void designateCutCellCandidates( std::vector< CutCell<nDim> >& cutCellData ); //ditch

  // general 3D
  void computeCutFaces(std::vector<CutCell<nDim>>& cutCellData, const MInt maxNoSurfaces, const MInt tCutGroup);
  void fillCutCellData(std::vector<CutCandidate<nDim>>& candidates,
                       std::vector<CutCell<nDim>>& cutCellData,
                       std::vector<MInt>
                           cutCellIdMapping);

  // scalarField based:
  void computeCutPoints(std::vector<CutCandidate<nDim>>& candidates, const MInt* candidateIds,
                        std::vector<MInt>& candidatesOrder);
  void computeNodalValues(std::vector<CutCandidate<nDim>>& candidates, MInt* candidateIds, const MFloat* scalarField,
                          const MInt* bodyIdField, const MBool* const gapPropertyField, const MBool gapClosure);
  void exchangeNodalValues(const MInt** maxLevelWindowCells, const MInt* noMaxLevelWindowCells,
                           const MInt** maxLevelHaloCells, std::vector<CutCandidate<nDim>>& candidates,
                           MInt* candidateIds);

  // stl-based:
  void computeCutPointsFromSTL(std::vector<CutCandidate<nDim>>& candidates);

  void clearCutCellData() { std::vector<CutCell<nDim>>().swap(m_cutCellData); }
  void resetCutCellData() { std::vector<CutCell<nDim>>().swap(m_cutCellData); }
  void returnDebugInfo() {
#ifdef CutCell_DEBUG_

    MPI_Allreduce(MPI_IN_PLACE, &m_caseCheckList[0], 256, MPI_INT, MPI_MAX, grid().mpiComm());

    if(grid().domainId() == 0) {
      std::cerr << "Tested cases: " << std::endl;
      for(MInt i = 0; i < 256; i++) {
        std::cerr << i << " " << m_caseCheckList[i] << std::endl;
      }
    }
#endif
  }

  // move to private later
  const MFloat m_eps;
  MBool m_complexBoundary = true;
  //  MBool m_createCoarseLevelCutCells = false;
  MInt m_noEmbeddedBodies{};
  MInt m_noLevelSetsUsedForMb{};
  //  MInt m_noBodyBndryCndIds;
  //  MInt* m_bodyBndryCndIds;
  MInt* m_bodyToSetTable{};
  //  MInt m_bodyToSetMode;
  MInt** m_setToBodiesTable{};
  MInt* m_noBodiesInSet{};

  MFloat* m_scaledCoordinate = nullptr;
  MBool m_scaledCutCell = true;
  MBool m_multiCutCell = false;


 private:
  GridProxy& grid() { return *m_gridProxy; }
  const GridProxy& grid() const { return *m_gridProxy; }
  GridProxy* const m_gridProxy;

  Geom& geometry() { return *m_geometry; }
  const Geom& geometry() const { return *m_geometry; }
  Geom* const m_geometry;

  GeometryContext& geometryContext() { return geometry().geometryContext(); }

  // general 3D function:
  MBool computeCutFaceSimple(CutCell<nDim>& cutCell);

  inline void crossProduct(MFloat* c, const MFloat* a, const MFloat* b) {
    c[0] = a[1] * b[2] - a[2] * b[1];
    c[1] = a[2] * b[0] - a[0] * b[2];
    c[2] = a[0] * b[1] - a[1] * b[0];
  }

  // stencil
  //  template <class T=MInt>
  //    typename std::enable_if< nDim_==2, T >::type
  //    print2() {return a_coordinate(0,1);}
  //  template <class T=MInt>
  //    typename std::enable_if< nDim_==3, T >::type
  //    print3() {return a_coordinate(0,2);}
  //---


  std::vector<CutCell<nDim>> m_cutCellData;
  std::vector<CutCandidate<nDim>> m_cutCandidates;

  MInt m_bodyFaceJoinMode;
  MFloat m_bodyFaceJoinCriterion;
  MString m_gridCutTest;

  MInt* m_caseCheckList = nullptr;

  MBool m_bndryLvlJumps = false;

  std::vector<lvlJumpCandidates<nDim>> m_cutLvlJumpCandidates;

  void correctNodalValuesAtLevelJump(std::vector<CutCandidate<nDim>>& candidates, const MInt*);

  void getNeighborNodes(const MInt, const MInt, MInt, MInt*, MInt*);

  //======================================
  //=== multi-cut-cell data structures ===
  //======================================

  class surfBase {
   public:
    MInt bodyId = -1;
    MFloat area = 0.0;
    std::array<MFloat, nDim> normal{};
    std::array<MFloat, nDim> center{};
    surfBase() = default;
    explicit surfBase(MInt _bodyId) : bodyId(_bodyId) {}
  };

  class polyVertex {
   public:
    std::array<MFloat, nDim> coordinates{};
    MInt pointId{};
    MInt pointType{};
    // 1: cut points;
    // 2: additional polyhedron clipping vertex
    // 3: additional MC vertex;
    MInt cartSrfcId{};
    std::vector<MInt> edges;   // corresponding edges
    std::vector<MInt> faceIds; // corresponding faces
    std::set<MInt> surfaceIdentificators;
    polyVertex(MInt pId, MInt pType) : pointId(pId), pointType(pType), cartSrfcId(-1) {}
    polyVertex(std::array<MFloat, nDim> coords, MInt pId, MInt pType) : polyVertex(pId, pType) { coordinates = coords; }
  };

  class polyEdgeBase {
   public:
    MInt vertices[2];
    MInt edgeType;
    // 0: Cartesian uncut, 1: Cartesian cut, 2: cut line on face, 3: cut line in cell
    MInt edgeId;
    // for edgeType 0/1: Cartesian edge, for edgeType 2/3: cartesian face or -1 if pure body line
    polyEdgeBase(MInt v0, MInt v1, MInt eId, MInt eType) : vertices{v0, v1}, edgeType(eType), edgeId(eId) {}
  };
  class polyEdge2D : public surfBase, public polyEdgeBase {
   public:
    MInt cutCell;
    MFloat w{};
    polyEdge2D(MInt v0, MInt v1, MInt eId, MInt eType, MInt _bodyId)
      : surfBase(_bodyId), polyEdgeBase(v0, v1, eId, eType), cutCell(-1) {}
  };
  class polyEdge3D : public polyEdgeBase {
   public:
    MInt face[2];
    polyEdge3D(MInt v0, MInt v1, MInt eId, MInt eType) : polyEdgeBase(v0, v1, eId, eType), face{-1, -1} {}
  };
  typedef typename std::conditional<nDim_ == 3, polyEdge3D, polyEdge2D>::type polyEdge;

  class polyFace : public surfBase {
   public:
    std::vector<std::pair<MInt, MInt>> edges;
    // contains 1. the edge and 2. the direction of the edge (1: as is, -1: reversed)
    std::vector<MInt> vertices;
    // contains the vertexIds of the face
    MInt cutCell = -1;
    MInt faceId = -1;
    MInt faceType = -1; // 0: cartesian face; 1: body face
    MFloat w = 0;
    // for cartesian faces:
    // this is the coordinate at which the face is positioned in normal-direction!
    MInt tmpSetIndex;
    MInt isLine;
    polyFace(MInt fId, MInt fType, MInt _bodyId)
      : surfBase(_bodyId), cutCell(-1), faceId(fId), faceType(fType), tmpSetIndex(-1), isLine(false) {}

   private:
    polyFace() : surfBase(-1), cutCell(-1), tmpSetIndex(-1), isLine(false) {}
  };

  //  class polyMultiEdge: public surfBase{
  //    public:
  //      std::vector<MInt> edges;
  //  };

  class polyMultiFace : public surfBase {
   public:
    std::vector<MInt> faces;
    std::list<std::pair<MInt, MInt>> edges;
    // contains 1. the edge and 2. the direction of the edge
    // direction can be 1/-1 (1: as is, -1: reversed)
    // and
  };

  class polyCutCellBase {
   public:
    MInt cartesianCell;
    MFloat volume{};
    MFloat center[nDim]{};
    polyCutCellBase(MInt cartCell, const MFloat* _center) : cartesianCell(cartCell) {
      std::copy(_center, _center + nDim, center);
    }
  };
  class polyCutCell2D : public polyCutCellBase {
   public:
    std::vector<MInt> edges;
    polyCutCell2D(MInt cartCell, const MFloat* _center) : polyCutCellBase(cartCell, _center) {}
  };
  class polyCutCell3D : public polyCutCellBase {
   public:
    std::vector<MInt> faces;
    polyCutCell3D(MInt cartCell, const MFloat* _center) : polyCutCellBase(cartCell, _center) {}
  };
  typedef typename std::conditional<nDim_ == 3, polyCutCell3D, polyCutCell2D>::type polyCutCell;

  class splitCartesianFace {
   public:
    MInt direction;
    std::vector<MInt> srfcIds;
    explicit splitCartesianFace(MInt dir) : direction(dir) {}
  };

  class cellWithSplitFace {
   public:
    MInt cellId;
    std::vector<splitCartesianFace> splitFaces;
    explicit cellWithSplitFace(MInt id) : cellId(id) {}
  };


  void addPoint(const std::vector<polyVertex>*, const MInt*, const MInt, std::array<MFloat, nDim>);
  void compVolumeIntegrals_pyraBased3(std::vector<polyCutCell>*, std::vector<polyFace>*,
                                      const std::vector<polyVertex>*);
  void compFaceIntegrals_pyraBased3(polyFace* face, const std::vector<polyVertex>* vertices, MFloat* MC, MFloat* VC,
                                    MFloat* XC);

  void writeVTKFileOfCell(MInt cellId, std::vector<polyFace>* faces, const std::vector<polyVertex>* vertices, MInt set);
  void writeInfo(std::vector<CutCell<nDim>>& cutCellData, MUint, MInt);


  /** \brief returns the normal corresponding to the triangle abc and returns the result in res
   *
   * \author Claudia Guenther, November 2013
   */
  template <class T = void>
  inline typename std::enable_if<nDim_ == 2, T>::type // 2D
  computeNormal(const std::array<MFloat, nDim> p0, const std::array<MFloat, nDim> p1, std::array<MFloat, nDim> res,
                MFloat& w) {
    const MFloat dx = p1[0] - p0[0], dy = p1[1] - p0[1];
    res[0] = -dy;
    res[1] = dx;
    const MFloat abs = sqrt(res[0] * res[0] + res[1] * res[1]);
    res[0] /= abs;
    res[1] /= abs;
    w = -res[0] * p0[0] - res[1] * p0[1];
    return;
  }
  template <class T = MFloat>
  inline typename std::enable_if<nDim_ == 3, T>::type // 3D
  computeNormal(const std::array<MFloat, nDim> p0, const std::array<MFloat, nDim> p1, const std::array<MFloat, nDim> p2,
                std::array<MFloat, nDim>& res, MFloat& w) {
    const MFloat a0 = p1[0] - p0[0], a1 = p1[1] - p0[1], a2 = p1[2] - p0[2], b0 = p2[0] - p0[0], b1 = p2[1] - p0[1],
                 b2 = p2[2] - p0[2];
    res[0] = a1 * b2 - a2 * b1;
    res[1] = a2 * b0 - a0 * b2;
    res[2] = a0 * b1 - a1 * b0;
    const MFloat abs = sqrt(res[0] * res[0] + res[1] * res[1] + res[2] * res[2]);
    res[0] /= abs;
    res[1] /= abs;
    res[2] /= abs;
    w = -res[0] * p0[0] - res[1] * p0[1] - res[2] * p0[2];
    return abs;
  }

  //-----------------------------------------------
  // ----------- Csg-classes:
  //-----------------------------------------------
  class CsgVector {
   public:
    std::array<MFloat, nDim> xx = make_array<MFloat, nDim>(0.0);

    CsgVector() = default;
    explicit CsgVector(const std::array<MFloat, nDim> X) { xx = X; }
    CsgVector(const CsgVector& Y) { xx = Y.xx; }
    CsgVector(CsgVector&&) = default;
    CsgVector& operator=(const CsgVector&) = default;
    CsgVector& operator=(CsgVector&&) = default;

    inline CsgVector clone() const {
      CsgVector tmp(xx);
      return tmp;
    }
    inline CsgVector plus(const CsgVector& a) const {
      CsgVector tmp(xx);
      for(MInt i = 0; i < nDim; i++)
        tmp.xx[i] = this->xx[i] + a.xx[i];
      return tmp;
    }
    inline CsgVector minus(const CsgVector& a) const {
      CsgVector tmp(xx);
      for(MInt i = 0; i < nDim; i++)
        tmp.xx[i] = this->xx[i] - a.xx[i];
      return tmp;
    }
    inline CsgVector times(const MFloat a) const {
      CsgVector tmp(xx);
      for(MInt i = 0; i < nDim; i++)
        tmp.xx[i] = this->xx[i] * a;
      return tmp;
    }
    inline CsgVector dividedBy(const MFloat a) const {
      CsgVector tmp(xx);
      for(MInt i = 0; i < nDim; i++)
        tmp.xx[i] = this->xx[i] / a;
      return tmp;
    }
    inline MFloat dot(const CsgVector& a) const {
      MFloat dotp = F0;
      for(MInt i = 0; i < nDim; i++)
        dotp += this->xx[i] * a.xx[i];
      return dotp;
    }
    inline CsgVector lerp(const CsgVector& a, const MFloat t) const { return this->plus((a.minus(*this)).times(t)); }
    inline MFloat length() const { return sqrt(this->dot(*this)); }
    inline CsgVector unit() const { return this->dividedBy(this->length()); }
    template <class T = CsgVector>
    typename std::enable_if<nDim_ == 3, T>::type // 3D
    cross(const CsgVector a) const {
      ASSERT(nDim == 3, "");
      CsgVector tmp;
      tmp.xx[0] = this->xx[1] * a.xx[2] - this->xx[2] * a.xx[1];
      tmp.xx[1] = this->xx[2] * a.xx[0] - this->xx[0] * a.xx[2];
      tmp.xx[2] = this->xx[0] * a.xx[1] - this->xx[1] * a.xx[0];
      return tmp;
    }
    inline void negate() {
      for(MInt i = 0; i < nDim; i++)
        this->xx[i] = -this->xx[i];
    }
  };


  class CsgVertex {
   public:
    CsgVector pos;
    MInt vertexId;
    MInt setIndex;
    CsgVertex(CsgVector _pos, MInt _vertexId, MInt _setIndex) : pos(_pos) {
      vertexId = _vertexId;
      setIndex = _setIndex;
    }
    CsgVertex clone() const {
      CsgVertex tmp(this->pos.clone(), vertexId, setIndex);
      return tmp;
    }
    CsgVertex interpolate(CsgVertex* other, MFloat t) const { return CsgVertex(this->pos.lerp(other->pos, t), -1, -1); }
  };

  class CsgPolygon;

  class CsgPlane {
   public:
    static constexpr MFloat eps = (nDim == 3) ? 10000 * MFloatEps : 1000 * MFloatEps;
    // tolerance used by splitPolygon() to decide if a point is on the plane or not!
    CsgVector normal{};
    MFloat w{};

    CsgPlane() = default;
    CsgPlane(CsgVector _normal, MFloat _w) : normal(_normal) { w = _w; }
    template <class T = CsgVector>
    typename std::enable_if<nDim_ == 2, T>::type // 2D
    initPlane(CsgVector abc[nDim]) const {
      CsgVector tmp = abc[1].minus(abc[0]).unit();
      return tmp;
    }
    template <class T = CsgVector>
    typename std::enable_if<nDim_ == 3, T>::type // 3D
    initPlane(CsgVector abc[nDim]) const {
      CsgVector tmp = abc[1].minus(abc[0]).cross(abc[2].minus(abc[0])).unit();
      return tmp;
    }
    explicit CsgPlane(CsgVector abc[nDim]) {
      CsgVector tmp = initPlane(abc);
      IF_CONSTEXPR(nDim == 3) { normal = tmp; }
      else IF_CONSTEXPR(nDim == 2) {
        normal.xx[0] = -tmp.xx[1];
        normal.xx[1] = tmp.xx[0];
      }
      w = normal.dot(abc[0]);
    }
    ~CsgPlane() = default;
    CsgPlane clone() {
      CsgPlane tmp(this->normal, this->w);
      return tmp;
    }
    void flip() {
      this->normal.negate();
      this->w = -this->w;
    }
    template <class T = void>
    typename std::enable_if<nDim_ == 2, T>::type // 2D
    splitPolygon(CsgPolygon* polygon,
                 std::vector<CsgPolygon>* coplanarFront,
                 std::vector<CsgPolygon>* coplanarBack,
                 std::vector<CsgPolygon>* front,
                 std::vector<CsgPolygon>* back) {
      const MInt COPLANAR = 0;
      const MInt FRONT = 1;
      const MInt BACK = 2;
      const MInt SPANNING = 3;

      MInt polygonType = 0;
      std::vector<MInt> types;

      for(MInt i = 0; (unsigned)i < polygon->vertices.size(); i++) {
        MFloat t = this->normal.dot((polygon->vertices[i].pos)) - this->w;
        MInt type = (t < -eps) ? BACK : (t > eps) ? FRONT : COPLANAR;
        polygonType |= type;
        types.push_back(type);
      }

      switch(polygonType) {
        case COPLANAR:
          (this->normal.dot((polygon->plane.normal)) > F0 ? coplanarFront : coplanarBack)->push_back(*polygon);
          break;
        case FRONT:
          front->push_back(*polygon);
          break;
        case BACK:
          back->push_back(*polygon);
          break;
        case SPANNING: {
          std::vector<CsgVertex> f;
          std::vector<CsgVertex> b;
          MInt i = 0;
          MInt j = 1;
          MInt ti = types[i], tj = types[j];
          CsgVertex vi = polygon->vertices[i];
          CsgVertex vj = polygon->vertices[j];
          if(ti != BACK) f.push_back(vi);
          if(ti != FRONT) b.push_back(ti != BACK ? vi.clone() : vi);
          if((ti | tj) == SPANNING) {
            MFloat t = (this->w - this->normal.dot((vi.pos))) / this->normal.dot(vj.pos.minus((vi.pos)));
            CsgVertex v = vi.interpolate(&vj, t);
            f.push_back(v);
            b.push_back(v.clone());
          }
          if(tj != BACK) f.push_back(vj);
          if(tj != FRONT) b.push_back(tj != BACK ? vj.clone() : vj);
          if(f.size() >= 2)
            front->push_back(
                CsgPolygon(f, polygon->setIndex, polygon->faceId, polygon->faceType, polygon->bodyId, polygon->plane));
          if(b.size() >= 2)
            back->push_back(
                CsgPolygon(b, polygon->setIndex, polygon->faceId, polygon->faceType, polygon->bodyId, polygon->plane));
        } break;
        default:
          break;
      }
    }
    template <class T = void>
    typename std::enable_if<nDim_ == 3, T>::type // 3D
    splitPolygon(CsgPolygon* polygon,
                 std::vector<CsgPolygon>* coplanarFront,
                 std::vector<CsgPolygon>* coplanarBack,
                 std::vector<CsgPolygon>* front,
                 std::vector<CsgPolygon>* back) {
      const MInt COPLANAR = 0;
      const MInt FRONT = 1;
      const MInt BACK = 2;
      const MInt SPANNING = 3;

      MInt polygonType = 0;
      std::vector<MInt> types;

      // test normals for equivalence
      MBool normalsEqual = false;
      if(std::abs(std::abs(this->normal.xx[0]) - std::abs(polygon->plane.normal.xx[0])) < MFloatEps
         && std::abs(std::abs(this->normal.xx[1]) - std::abs(polygon->plane.normal.xx[1])) < MFloatEps
         && std::abs(std::abs(this->normal.xx[2]) - std::abs(polygon->plane.normal.xx[2])) < MFloatEps)
        normalsEqual = true;

      MBool coplanarVertex = false;
      for(MInt i = 0; (unsigned)i < polygon->vertices.size(); i++) {
        MFloat t = this->normal.dot((polygon->vertices[i].pos)) - this->w;
        MInt type = (t < -eps) ? BACK : (t > eps) ? FRONT : COPLANAR;
        polygonType |= type;
        types.push_back(type);
        if(type == 0) coplanarVertex = true;
      }
      ASSERT(!(std::isnan(polygon->plane.normal.xx[0])), "");

      if(normalsEqual && coplanarVertex) polygonType = COPLANAR;

      switch(polygonType) {
        case COPLANAR:
          (this->normal.dot((polygon->plane.normal)) > F0 ? coplanarFront : coplanarBack)->push_back(*polygon);
          break;
        case FRONT:
          front->push_back(*polygon);
          break;
        case BACK:
          back->push_back(*polygon);
          break;
        case SPANNING: {
          std::vector<CsgVertex> f;
          std::vector<CsgVertex> b;
          for(MInt i = 0; (unsigned)i < polygon->vertices.size(); i++) {
            MInt j = (i + 1) % polygon->vertices.size();
            MInt ti = types[i], tj = types[j];
            CsgVertex vi = polygon->vertices[i];
            CsgVertex vj = polygon->vertices[j];
            if(ti != BACK) f.push_back(vi);
            if(ti != FRONT) b.push_back(ti != BACK ? vi.clone() : vi);
            if((ti | tj) == SPANNING) {
              MFloat t = (this->w - this->normal.dot((vi.pos))) / this->normal.dot(vj.pos.minus((vi.pos)));
              CsgVertex v = vi.interpolate(&vj, t);
              f.push_back(v);
              b.push_back(v.clone());
            }
          }
          if(f.size() >= 3)
            front->push_back(
                CsgPolygon(f, polygon->setIndex, polygon->faceId, polygon->faceType, polygon->bodyId, polygon->plane));
          if(b.size() >= 3)
            back->push_back(
                CsgPolygon(b, polygon->setIndex, polygon->faceId, polygon->faceType, polygon->bodyId, polygon->plane));
        } break;
        default:
          break;
      }
    }

    void insertCoplanarPolygon(CsgPolygon* polygon,
                               std::vector<CsgPolygon>* coplanarFront,
                               std::vector<CsgPolygon>* coplanarBack) {
      ASSERT(!(std::isnan(polygon->plane.normal.xx[0])), "");
      (this->normal.dot((polygon->plane.normal)) > F0 ? coplanarFront : coplanarBack)->push_back(*polygon);
    }
  };

  class CsgPolygon {
   public:
    std::vector<CsgVertex> vertices;
    // helpers: edgeType[0] corresponds to edge between vertices[0] and vertices[1],...
    CsgPlane plane;
    MInt setIndex{};
    MInt bodyId{};
    MInt faceId{};
    MInt faceType{}; // 0: cartesian face; 1: body face

    CsgPolygon(std::vector<CsgVertex> _vertices, MInt _setIndex, MInt _faceId, MInt _faceType, MInt _bodyId) {
      for(MInt i = 0; (unsigned)i < _vertices.size(); i++) {
        vertices.push_back(_vertices[i].clone());
      }
      setIndex = _setIndex;
      faceId = _faceId;
      faceType = _faceType;
      bodyId = _bodyId;
      CsgVector abc[nDim];
      for(MInt i = 0; i < nDim; i++) {
        abc[i] = vertices[i].pos;
      }
      plane = CsgPlane(abc);
    }

    CsgPolygon(std::vector<CsgVertex> _vertices, MInt _setIndex, MInt _faceId, MInt _faceType, MInt _bodyId,
               CsgPlane _plane) {
      for(MInt i = 0; (unsigned)i < _vertices.size(); i++) {
        vertices.push_back(_vertices[i].clone());
      }
      setIndex = _setIndex;
      faceId = _faceId;
      faceType = _faceType;
      bodyId = _bodyId;
      plane = _plane.clone();
    }

    ~CsgPolygon() = default;

    CsgPolygon clone() {
      CsgPolygon tmp(vertices, setIndex, faceId, faceType, bodyId, plane);
      return tmp;
    }

    void flip() {
      this->plane.flip();
      std::vector<CsgVertex> tmp;
      for(auto rit = vertices.rbegin(); rit != vertices.rend(); ++rit) {
        tmp.push_back((*rit));
      }
      vertices.swap(tmp);
    }

   private:
    CsgPolygon() = default;
  };

  class CsgNode {
   public:
    CsgPlane plane;
    MBool planeValid;
    CsgNode* front;
    CsgNode* back;
    std::vector<CsgPolygon> polygons;

    explicit CsgNode(std::vector<CsgPolygon> _polygons) {
      planeValid = false;
      front = nullptr;
      back = nullptr;
      polygons.clear();
      this->build(_polygons);
    }

    CsgNode() {
      planeValid = false;
      front = nullptr;
      back = nullptr;
      polygons.clear();
    }

    ~CsgNode() {
      this->polygons.clear();
      delete this->front;
      delete this->back;
    }

    void invert() {
      for(MInt i = 0; (unsigned)i < this->polygons.size(); i++) {
        this->polygons[i].flip();
      }
      this->plane.flip();
      if(this->front) this->front->invert();
      if(this->back) this->back->invert();
      CsgNode* temp = this->front;
      this->front = this->back;
      this->back = temp;
    }

    std::vector<CsgPolygon> clipPolygons(std::vector<CsgPolygon> _polygons) {
      if(!planeValid) {
        return _polygons;
      }
      std::vector<CsgPolygon> _front;
      std::vector<CsgPolygon> _back;
      for(MInt i = 0; (unsigned)i < _polygons.size(); i++) {
        this->plane.splitPolygon(&_polygons[i], &_front, &_back, &_front, &_back);
      }
      if(this->front) _front = this->front->clipPolygons(_front);
      if(this->back)
        _back = this->back->clipPolygons(_back);
      else
        _back.clear();
      for(MInt i = 0; (unsigned)i < _back.size(); i++)
        _front.push_back(_back[i]);
      return _front;
    }

    void clipTo(CsgNode& bsp) {
      this->polygons = bsp.clipPolygons(this->polygons);
      if(this->front) this->front->clipTo(bsp);
      if(this->back) this->back->clipTo(bsp);
    }

    std::vector<CsgPolygon> allPolygons() {
      std::vector<CsgPolygon> _polygons(this->polygons);
      if(this->front) {
        std::vector<CsgPolygon> _polygons_front = this->front->allPolygons();
        for(MInt i = 0; (unsigned)i < _polygons_front.size(); i++)
          _polygons.push_back(_polygons_front[i]);
      }
      if(this->back) {
        std::vector<CsgPolygon> _polygons_back = this->back->allPolygons();
        for(MInt i = 0; (unsigned)i < _polygons_back.size(); i++)
          _polygons.push_back(_polygons_back[i]);
      }
      return _polygons;
    }

    // void build(std::vector<CsgPolygon> _polygons);

    template <class T = void>
    typename std::enable_if<nDim_ == 2, T>::type build(std::vector<CsgPolygon> _polygons) {
      if(_polygons.empty()) return;
      if(!planeValid) {
        this->plane = _polygons[0].plane.clone();
        planeValid = true;
      }
      std::vector<CsgPolygon> _front;
      std::vector<CsgPolygon> _back;
      this->plane.insertCoplanarPolygon(&_polygons[0], &this->polygons, &this->polygons);
      for(MInt i = 1; (unsigned)i < _polygons.size(); i++) {
        this->plane.splitPolygon(&_polygons[i], &this->polygons, &this->polygons, &_front, &_back);
      }
      if(!_front.empty()) {
        if(!this->front) this->front = new CsgNode();
        this->front->build(_front);
      }
      if(!_back.empty()) {
        if(!this->back) this->back = new CsgNode();
        this->back->build(_back);
      }
    }
    template <class T = void>
    typename std::enable_if<nDim_ == 3, T>::type build(std::vector<CsgPolygon> _polygons) {
      if(_polygons.empty()) return;
      std::vector<CsgPolygon> _front;
      std::vector<CsgPolygon> _back;
      MInt polygonStartIndex = 0;
      if(!planeValid) {
        this->plane = _polygons[0].plane.clone();
        planeValid = true;
        this->polygons.push_back(_polygons[0]);
        polygonStartIndex = 1;
      }
      for(MInt i = polygonStartIndex; (unsigned)i < _polygons.size(); i++) {
        this->plane.splitPolygon(&_polygons[i], &this->polygons, &this->polygons, &_front, &_back);
      }
      if(!_front.empty()) {
        if(!this->front) this->front = new CsgNode();
        this->front->build(_front);
      }
      if(!_back.empty()) {
        if(!this->back) this->back = new CsgNode();
        this->back->build(_back);
      }
    }
  };

  class Csg {
   public:
    std::vector<CsgPolygon> polygons;

    Csg() { polygons.clear(); }

    explicit Csg(std::vector<CsgPolygon> _polygons) {
      for(MInt i = 0; (unsigned)i < _polygons.size(); i++) {
        polygons.push_back(_polygons[i].clone());
      }
    }

    ~Csg() = default;

    std::vector<CsgPolygon> toPolygons() {
      std::vector<CsgPolygon> _polygons(this->polygons);
      return _polygons;
    }


    std::vector<CsgPolygon> intersect(Csg csg) {
      CsgNode a(this->polygons);
      CsgNode b(csg.polygons);

      b.invert();
      b.clipTo(a);
      b.invert();
      a.clipTo(b);
      b.clipTo(a);
      a.build(b.allPolygons());

      return a.allPolygons();
    }

    /*
      //untested!!!!!!!!!!!!!!!!!!!!!!!
      //for concave boundaries!!
      std::vector<CsgPolygon> merge(Csg csg){
        CsgNode a(this->polygons);
        CsgNode b(csg.polygons);

        a.clipTo(b);
        b.clipTo(a);
        b.invert();
        b.clipTo(a);
        b.invert();

        a.build(b.allPolygons());

        return a.allPolygons();
      }*/

    void inverse() {
      for(MInt i = 0; (unsigned)i < this->polygons.size(); i++)
        this->polygons[i].flip();
    }
  };
};

// template class GeometryIntersection<2>;
// template class GeometryIntersection<3>;


// ------------------------------------------------------------------------
// ------------------------------------------------------------------------
// ------------------------------------------------------------------------
// ---------------------- ADDITIONAL HELPER CLASSES -----------------------
// ------------------------------------------------------------------------
// ------------------------------------------------------------------------
// ------------------------------------------------------------------------

template <MInt nDim>
class lvlJumpCandidates {
 public:
  static constexpr MInt maxNoJumps = nDim == 2 ? 3 : 7;

  MInt candId = -1;

  // parent info
  MInt parentCandId = -1;
  MInt childId = -1;

  MInt noJumps = 0;
  std::array<MInt, maxNoJumps> dirs = make_array<MInt, maxNoJumps>(-1);
  std::array<MInt, maxNoJumps> diagonalDirs = make_array<MInt, maxNoJumps>(-1);

  // 0: direct level-jump
  // 1: 2D-diagonal level-jump
  // 2: 3D-diagonal level-jump
  std::array<MInt, maxNoJumps> neighborType = make_array<MInt, maxNoJumps>(-1);
};

template <MInt nDim>
class CutCell {
 public:
  static constexpr MInt maxNoSets = 12;
  static constexpr MInt noFaces = 2 * nDim;
  static constexpr MInt noEdges = 2 * nDim * (nDim - 1);
  static constexpr MInt noCorners = IPOW2(nDim);
  static constexpr MInt maxNoCutPoints = 4 * noEdges;
  static constexpr MInt maxNoCartesianSurfaces = 2 * noFaces;
  static constexpr MInt maxNoBoundarySurfaces = maxNoSets;
  static constexpr MInt maxNoTotalSurfaces = maxNoBoundarySurfaces + maxNoBoundarySurfaces;
  static constexpr MInt maxNoFaceVertices = 4 * IPOW2(nDim - 1) + 4;
  static constexpr MInt maxNoAdditionalVertices = 4 * noFaces;
  static constexpr MInt maxSplitCells = noFaces;
  static constexpr MInt no3dFaces = (nDim == 3) ? noFaces : 1;

  // make private if transition complete
  // private:

  MInt cellId = -1;

  MFloat volume = F0;
  std::array<MFloat, nDim> volumetricCentroid = make_array<MFloat, nDim>(0.0);

  MInt noCutPoints = 0;
  std::array<std::array<MFloat, nDim>, maxNoCutPoints> cutPoints = [] {
    std::array<std::array<MFloat, nDim>, maxNoCutPoints> a;
    std::array<MFloat, nDim> b = make_array<MFloat, nDim>(0.0);
    a.fill(b);
    return a;
  }();
  std::array<MInt, maxNoCutPoints> cutBodyIds = make_array<MInt, maxNoCutPoints>(-1);
  std::array<MInt, maxNoCutPoints> cutEdges = make_array<MInt, maxNoCutPoints>(-1);

  std::array<std::array<MBool, no3dFaces>, maxNoSets> faceCentroidIsInsideGeometry = [] {
    std::array<std::array<MBool, no3dFaces>, maxNoSets> a;
    std::array<MBool, no3dFaces> b = make_array<MBool, no3dFaces>(0.0);
    a.fill(b);
    return a;
  }();
  std::array<std::array<MBool, noCorners>, maxNoSets> cornerIsInsideGeometry = [] {
    std::array<std::array<MBool, noCorners>, maxNoSets> a;
    std::array<MBool, noCorners> b = make_array<MBool, noCorners>(false);
    a.fill(b);
    return a;
  }();
  std::array<MInt, maxNoSets> associatedBodyIds = make_array<MInt, maxNoSets>(-1);

  std::array<MInt, noFaces> noFacesPerCartesianDir = make_array<MInt, noFaces>(0);
  std::array<MBool, noFaces> externalFaces = make_array<MBool, noFaces>(false);

  MBool isGapCell = false;
  MInt splitParentId = -1;             // link to original cell for split cells
  MInt noSplitChilds = 0;              // number of distinct parts this cell is split into
  std::array<MInt, maxSplitCells> splitChildIds = make_array<MInt, maxSplitCells>(-1);

  MInt noCartesianSurfaces = 0;
  std::array<std::array<MFloat, nDim>, maxNoCartesianSurfaces> cartFaceCentroid = [] {
    std::array<std::array<MFloat, nDim>, maxNoCartesianSurfaces> a;
    std::array<MFloat, nDim> b = make_array<MFloat, nDim>(0.0);
    a.fill(b);
    return a;
  }();
  std::array<MFloat, maxNoCartesianSurfaces> cartFaceArea = make_array<MFloat, maxNoCartesianSurfaces>(0.0);
  std::array<MInt, maxNoCartesianSurfaces> cartFaceDir = make_array<MInt, maxNoCartesianSurfaces>(-1);

  MInt noBoundarySurfaces = 0;
  std::array<std::array<MFloat, nDim>, maxNoBoundarySurfaces> boundarySurfaceCentroid = [] {
    std::array<std::array<MFloat, nDim>, maxNoBoundarySurfaces> a;
    std::array<MFloat, nDim> b = make_array<MFloat, nDim>(0.0);
    a.fill(b);
    return a;
  }();
  std::array<std::array<MFloat, nDim>, maxNoBoundarySurfaces> boundarySurfaceNormal = [] {
    std::array<std::array<MFloat, nDim>, maxNoBoundarySurfaces> a;
    std::array<MFloat, nDim> b = make_array<MFloat, nDim>(0.0);
    a.fill(b);
    return a;
  }();
  std::array<MFloat, maxNoBoundarySurfaces> boundarySurfaceArea = make_array<MFloat, maxNoBoundarySurfaces>(0.0);
  std::array<MInt, maxNoBoundarySurfaces> boundarySurfaceBodyId = make_array<MInt, maxNoBoundarySurfaces>(-1);
  std::array<MInt, maxNoBoundarySurfaces> boundarySurfaceBndryCndId = make_array<MInt, maxNoBoundarySurfaces>(-1);

  MInt noTotalFaces = 0;
  std::array<MInt, maxNoTotalSurfaces> allFacesBodyId = make_array<MInt, maxNoTotalSurfaces>(-1); // rather bndSrfId ???
  std::array<MInt, maxNoTotalSurfaces> allFacesNoPoints =
      make_array<MInt, maxNoTotalSurfaces>(0); // rather bndSrfId ???
  std::array<std::array<MInt, maxNoFaceVertices>, maxNoTotalSurfaces> allFacesPointIds = [] {
    std::array<std::array<MInt, maxNoFaceVertices>, maxNoTotalSurfaces> a;
    std::array<MInt, maxNoFaceVertices> b = make_array<MInt, maxNoFaceVertices>(-1);
    a.fill(b);
    return a;
  }();

  MInt noAdditionalVertices = 0; // additional vertices due to multi-cuts, polyhedron clipping
  // their coords
  std::array<std::array<MFloat, nDim>, maxNoAdditionalVertices> additionalVertices = [] {
    std::array<std::array<MFloat, nDim>, maxNoAdditionalVertices> a;
    std::array<MFloat, nDim> b = make_array<MFloat, nDim>(0.0);
    a.fill(b);
    return a;
  }();
};

// TODO labels:GEOM move the entire class to private and update/only hold information in the gridIntersetion!
template <MInt nDim>
class CutCandidate {
 public:
  static constexpr MInt maxNoSets = 12;
  static constexpr MInt noEdges = 2 * nDim * (nDim - 1);
  static constexpr MInt noCorners = IPOW2(nDim);
  static constexpr MInt maxNoCutPoints = 4 * noEdges;

  // make private if transition complete
  // private:

  // properties set in determineBndryCandidates
  MInt cellId = -1;

  // properties set in computeNodalValues
  MBool isGapCell = false;
  std::array<std::array<MFloat, noCorners>, maxNoSets> nodalValues = [] {
    std::array<std::array<MFloat, noCorners>, maxNoSets> a{};
    std::array<MFloat, noCorners> b = make_array<MFloat, noCorners>(-1.0);
    a.fill(b);
    return a;
  }();

  std::array<MBool, noCorners> nodeValueSet = make_array<MBool, noCorners>(false);
  std::array<MInt, maxNoSets> associatedBodyIds = make_array<MInt, maxNoSets>(-1);

  // properties set in computeCutPoints
  MInt noCutPoints = 0;
  std::array<std::array<MFloat, nDim>, maxNoCutPoints> cutPoints = [] {
    std::array<std::array<MFloat, nDim>, maxNoCutPoints> a{};
    std::array<MFloat, nDim> b = make_array<MFloat, nDim>(0.0);
    a.fill(b);
    return a;
  }();
  std::array<MInt, noEdges> noCutPointsOnEdge = make_array<MInt, noEdges>(0);
  std::array<MInt, maxNoCutPoints> cutBodyIds = make_array<MInt, maxNoCutPoints>(-1);
  std::array<MInt, maxNoCutPoints> cutEdges = make_array<MInt, maxNoCutPoints>(-1);
  std::array<MBool, noEdges> edgeChecked = make_array<MBool, noEdges>(false);

  // properties related to bndryLvlJump
  MBool isbndryLvlJumpParent = false;
};

#endif
