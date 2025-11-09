// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef LSSOLVER_H_
#define LSSOLVER_H_

#include <algorithm>
#include "GRID/cartesiangrid.h"
#include "cartesiansolver.h"
#include "enums.h"
#include "lscartesiancellcollector.h"
#include "lscartesiancontrolpoint.h"

#define ENSURE_VALID_GCELLID(gCellId)                                                                                  \
  do {                                                                                                                 \
    ASSERT(gCellId >= 0 && /*gCellId < m_cells.size() &&*/ gCellId < m_maxNoCells,                                     \
           "gCellId " << gCellId << " out-of-bounds [0, " << m_maxNoCells << ") m_cells.size: " << m_cells.size());    \
  } while(false)

#define ENSURE_VALID_DIM(dim)                                                                                          \
  do {                                                                                                                 \
    ASSERT(dim >= 0 && dim < nDim, "dim " << dim << " out-of-bounds [0, " << nDim << ")");                             \
  } while(false)

#define ENSURE_VALID_DIR(dir)                                                                                          \
  do {                                                                                                                 \
    ASSERT(dir >= 0 && dir < 2 * nDim, "dir " << dir << " out-of-bounds [0, " << 2 * nDim << ")");                     \
  } while(false)

#define ENSURE_VALID_SET(set)                                                                                          \
  do {                                                                                                                 \
    ASSERT(set >= 0 && set < m_maxNoSets, "set " << set << " out-of-bounds [0," << m_maxNoSets << ")");                \
  } while(false)
template <MInt nDim_>
class CHECKNORMAL {
 private:
  // vectors dot product
  MFloat dot(const MFloat p[3], const MFloat q[3]) { return p[0] * q[0] + p[1] * q[1] + p[2] * q[2]; }

  void cross(const MFloat p[3], const MFloat q[3], MFloat r[3]) {
    r[0] = (p[1] * q[2]) - (q[1] * p[2]);
    r[1] = (p[2] * q[0]) - (q[2] * p[0]);
    r[2] = (p[0] * q[1]) - (q[0] * p[1]);
  }

  // compute the transformation matrix for the current geometrical element
  void transformationmatrix(GeometryElement<nDim_>* el, MFloat tr[3][3]) {
    // Step 1: get vertices in global coordinates
    MFloat v[3][3] = {{F0, F0, F0}, {F0, F0, F0}, {F0, F0, F0}};
    MInt noVertices = nDim_;

    for(MInt i = 0; i < noVertices; i++) {
      for(MInt j = 0; j < nDim_; j++) {
        v[i][j] = el->m_vertices[i][j];
      }
    }

    // Step 2: define the transformation matrix for the current geometrical element
    IF_CONSTEXPR(nDim_ == 2) {
      MFloat theta_rad = 0;
      MFloat dx = v[1][0] - v[0][0];
      if(approx(dx, 0.0, MFloatEps))
        theta_rad = 3.141592654 / 2.0;
      else
        theta_rad = atan((v[1][1] - v[0][1]) / dx);

      tr[0][0] = cos(theta_rad);
      tr[0][1] = sin(theta_rad);
      tr[2][0] = 0;
      tr[1][0] = -sin(theta_rad);
      tr[1][1] = cos(theta_rad);
      tr[2][1] = 0;
      tr[0][2] = 0;
      tr[1][2] = 0;
      tr[2][2] = 0;
    }
    else {
      // Euclidean frame
      MFloat xm[3] = {1, 0, 0};
      MFloat ym[3] = {0, 1, 0};
      MFloat zm[3] = {0, 0, 1};
      // n subscript means rotated coordinates
      // set zn = normal vector
      MFloat zn[3];
      for(MInt i = 0; i < 3; i++) {
        zn[i] = el->m_normal[i];
      }

      // set xn to the direction of first edge
      MFloat xn[3];
      for(MInt i = 0; i < 3; i++) {
        xn[i] = v[1][i] - v[0][i];
      }

      // Normalize xn
      MFloat mag = sqrt(pow(xn[0], 2) + pow(xn[1], 2) + pow(xn[2], 2));
      for(MFloat& i : xn) {
        i = i / mag;
      }

      // yn is the outer product of zn x xn
      MFloat yn[3];
      cross(zn, xn, yn);
      // defining transformation matrix
      tr[0][0] = dot(xn, xm);
      tr[0][1] = dot(xn, ym);
      tr[0][2] = dot(xn, zm);
      tr[1][0] = dot(yn, xm);
      tr[1][1] = dot(yn, ym);
      tr[1][2] = dot(yn, zm);
      tr[2][0] = dot(zn, xm);
      tr[2][1] = dot(zn, ym);
      tr[2][2] = dot(zn, zm);
    }
  }

 public:
  // rotate the coordinates q into the geometrical element coordinates
  inline void rotation(const MFloat q[3], MFloat r[3], MFloat tr[3][3]) {
    for(MInt i = 0; i < 3; i++) {
      r[i] = 0;
      for(MInt j = 0; j < 3; j++) {
        r[i] += tr[i][j] * q[j];
      }
    }
  }

  // check if point q is inside the geometrical triangulated element after the projection
  MInt PointInsideTriangle(GeometryElement<nDim_>* el, MFloat q[3], MFloat tr[3][3]) {
    // Set transformation matrix
    transformationmatrix(el, tr);

    // Copy vertices to vv[][]
    MFloat vv[3][3] = {{F0, F0, F0}, {F0, F0, F0}, {F0, F0, F0}};

    MFloat noVertices;
    noVertices = nDim_;

    for(MInt i = 0; i < noVertices; i++)
      for(MInt j = 0; j < nDim_; j++) {
        vv[i][j] = el->m_vertices[i][j];
      }

    // Rotation of vertices
    MFloat ro[3] = {0.0, 0.0, 0.0}, uv[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, quv[2] = {0.0, 0.0};
    for(MInt i = 0; i < noVertices; i++) {
      rotation(vv[i], ro, tr);
      for(MInt j = 0; j < nDim_ - 1; j++) {
        uv[(nDim_ - 1) * i + j] = ro[j];
      }
    }

    // Rotation of q
    rotation(q, ro, tr);
    for(MInt j = 0; j < nDim_ - 1; j++)
      quv[j] = ro[j];

    // Move origin to q
    for(MInt j = 0; j < nDim_ - 1; j++) {
      for(MInt i = 0; i < noVertices; i++)
        uv[(nDim_ - 1) * i + j] -= quv[j];
    }

    // Define intersection
    IF_CONSTEXPR(nDim_ == 2) {
      MInt SH = 0, NSH = 0;

      if(uv[0] > 0)
        SH = 1;
      else if(approx(uv[0], 0.0, MFloatEps))
        SH = 0;
      else
        SH = -1;

      if(uv[1] > 0)
        NSH = 1;
      else if(approx(uv[1], 0.0, MFloatEps))
        NSH = 0;
      else
        NSH = -1;

      if((SH != NSH) || (SH == 0) || (NSH == 0)) return 1;
    }
    else {
      MInt NC = 0;

      for(MInt i = 0; i < noVertices; i++) {
        MInt ni = (i + 1) % 3; // next i
        // cast ray in +u direction
        MInt SH = 0, NSH = 0;

        if(uv[i * 2 + 1] > 0)
          SH = 1;
        else if(approx(uv[i * 2 + 1], 0.0, MFloatEps))
          SH = 0;
        else
          SH = -1;

        if(uv[ni * 2 + 1] > 0)
          NSH = 1;
        else if(approx(uv[ni * 2 + 1], 0.0, MFloatEps))
          NSH = 0;
        else
          NSH = -1;

        // predict intersection
        if(SH == 0 && NSH == 0) {
          NC = 1;
          break;
        } // horizontal line through origin
        else if((SH == 0 && approx(uv[i * 2 + 0], 0.0, MFloatEps))
                || (NSH == 0 && approx(uv[ni * 2 + 0], 0.0, MFloatEps))) {
          NC = 1;
          break;
        } else if(approx((uv[i * 2 + 0]
                          - uv[i * 2 + 1] * (uv[ni * 2 + 0] - uv[i * 2 + 0]) / (uv[ni * 2 + 1] - uv[i * 2 + 1])),
                         0.0, MFloatEps)) {
          NC = 1;
          break;
        } // line through origin
        else if(SH != NSH) {
          if((uv[i * 2 + 0] > 0) && (uv[ni * 2 + 0] > 0))
            NC++;                                                // there is intersection in +u
          else if((uv[i * 2 + 0] > 0) || (uv[ni * 2 + 0] > 0)) { // there can be an intersection in +u
            // need to calculate
            if((uv[i * 2 + 0] - uv[i * 2 + 1] * (uv[ni * 2 + 0] - uv[i * 2 + 0]) / (uv[ni * 2 + 1] - uv[i * 2 + 1]))
               > 0)
              NC++;
          }
        }
      }
      if((NC % 2) != 0) {
        return 1;
      }
    }
    return 0;
  }
};

namespace maia {
namespace ls {

// Create struct for easy timer identification
struct Timers_ {
  // Enum to store timer "names"
  enum {
    Solver,
    PreTime,
    TimeInt,
    PostTime,

    FirstEx,

    Finalize,
    BuildTube,
    SetBand,
    BuildMultiple,

    // Special enum value used to initialize timer array
    _count
  };
};

} // namespace ls
} // namespace maia

template <MInt nDim_, class SysEqn>
class LsFvCombustion;

template <MInt nDim_, class SysEqn>
class LsFvMb;

template <MInt nDim_, class SysEqn>
class CouplingLsFv;

template <RansMethod _>
struct RANSModelConstants;

template <MInt nDim_>
class FvSysEqnNS;

template <MInt nDim_, class RANSModel>
class FvSysEqnRANS;

class Coupling;

template <MInt nDim_>
class LsCartesianSolver : public maia::CartesianSolver<nDim_, LsCartesianSolver<nDim_>> {
 public:
  // Type for cell properties
  using Cell = typename maia::grid::tree::Tree<nDim_>::Cell;
  using GCellCollector = maia::ls::collector::GCells<nDim_>;
  using PropertyReference = typename GCellCollector::BitsetType::reference;

 protected:
  MBool m_static_createBaseGgrid_firstRun = true;

 private:

  MPI_Request* mpi_request = nullptr;
  MPI_Request* mpi_recive = nullptr;

  using Timers = maia::ls::Timers_;
  // Timers
  // Timer group which holds all solver-wide timers
  MInt m_timerGroup = -1;
  // Stores all solver-wide timers
  std::array<MInt, Timers::_count> m_timers{};

 public:
  static constexpr MInt nDim = nDim_;
  friend class LsFvCombustion<nDim, FvSysEqnNS<nDim>>;
  friend class LsFvCombustion<nDim, FvSysEqnRANS<3, RANSModelConstants<RANS_SA_DV>>>;
  friend class LsFvCombustion<nDim, FvSysEqnRANS<3, RANSModelConstants<RANS_FS>>>;
  friend class LsFvCombustion<nDim, FvSysEqnRANS<3, RANSModelConstants<RANS_KOMEGA>>>;
  friend class LsFvMb<nDim, FvSysEqnNS<nDim>>;
  friend class LsFvMb<nDim, FvSysEqnRANS<3, RANSModelConstants<RANS_SA_DV>>>;
  friend class LsFvMb<nDim, FvSysEqnRANS<3, RANSModelConstants<RANS_FS>>>;
  friend class LsFvMb<nDim, FvSysEqnRANS<3, RANSModelConstants<RANS_KOMEGA>>>;
  friend class CouplingLsFv<nDim, FvSysEqnNS<nDim>>;
  friend class CouplingLsFv<nDim, FvSysEqnRANS<3, RANSModelConstants<RANS_SA_DV>>>;
  friend class CouplingLsFv<nDim, FvSysEqnRANS<3, RANSModelConstants<RANS_FS>>>;
  friend class CouplingLsFv<nDim, FvSysEqnRANS<3, RANSModelConstants<RANS_KOMEGA>>>;

  GCellCollector m_cells;

  // returns the time of the lsSolver
  // overrides the Solver-function
  MFloat time() const override {
    if(m_levelSetMb || m_LSSolver || m_levelSetFv || m_freeSurface) {
      // NOTE: the time is copied and updated from the coupler in transferTimeStep
      return m_time;
    } else {
      ASSERT(!m_combustion, "");
      return globalTimeStep;
    }
  };

  /// \brief Returns the timeStep
  MFloat timeStep() const { return m_timeStep; }
  MInt noVariables() const override { return 1; };
  virtual void saveSolverSolution(const MBool, const MBool){};
  virtual void cleanUp(){};

  static constexpr const MInt m_noCorners = (nDim == 2) ? 4 : 8;

  using CartesianSolver = typename maia::CartesianSolver<nDim, LsCartesianSolver>;
  using Grid = typename CartesianSolver::Grid;
  using GridProxy = typename CartesianSolver::GridProxy;

  // used CartesianSolver
  using CartesianSolver::grid;
  using Geom = Geometry<nDim>;

  Geom* m_geometry;

  /// Access the solver's geometry
  Geom& geometry() const { return *m_geometry; }

  using CartesianSolver::disableDlbTimers;
  using CartesianSolver::domainId;
  using CartesianSolver::domainOffset;
  using CartesianSolver::enableDlbTimers;
  using CartesianSolver::exchangeData;
  using CartesianSolver::haloCellId;
  using CartesianSolver::isActive;
  using CartesianSolver::m_freeIndices;
  using CartesianSolver::m_initFromRestartFile;
  using CartesianSolver::m_maxNoSets;
  using CartesianSolver::m_restart;
  using CartesianSolver::m_restartFile;
  using CartesianSolver::m_restartInterval;
  using CartesianSolver::m_restartTimeStep;
  using CartesianSolver::m_solutionInterval;
  using CartesianSolver::m_solverId;
  using CartesianSolver::m_useNonSpecifiedRestartFile;
  using CartesianSolver::maxLevel;
  using CartesianSolver::maxNoGridCells;
  using CartesianSolver::maxRefinementLevel;
  using CartesianSolver::maxUniformRefinementLevel;
  using CartesianSolver::minLevel;
  using CartesianSolver::mpiComm;
  using CartesianSolver::neighborDomain;
  using CartesianSolver::noDomains;
  using CartesianSolver::noHaloCells;
  using CartesianSolver::noNeighborDomains;
  using CartesianSolver::noWindowCells;
  using CartesianSolver::outputDir;
  using CartesianSolver::readSolverSamplingVarNames;
  using CartesianSolver::restartDir;
  using CartesianSolver::returnIdleRecord;
  using CartesianSolver::returnLoadRecord;
  using CartesianSolver::solverId;
  using CartesianSolver::solverMethod;
  using CartesianSolver::startLoadTimer;
  using CartesianSolver::stopLoadTimer;
  using CartesianSolver::updateDomainInfo;
  using CartesianSolver::windowCellId;

 public:

  // Ls-Solver-Constructor:
  LsCartesianSolver<nDim_>(MInt, const MBool*, GridProxy& gridProxy_, Geom& geometry_, const MPI_Comm comm);

  CHECKNORMAL<nDim>* m_checkNormal;
  CHECKNORMAL<nDim>& checkNormal() const { return *m_checkNormal; }

  MInt noSolverTimers(const MBool allTimings) override {
#ifdef MAIA_TIMER_FUNCTION
    TERMM_IF_COND(!allTimings, "FIXME: reduced timings mode not yet supported by LS.");
    static const MInt noAdditionTimers = 7;
    return 2 + noAdditionTimers;
#else
    return 2;
#endif
  }

  //----------------------------------------------------------------------------------------------------
  MInt maxLevel() const { return grid().maxLevel(); }
  // cellLengthAtGCell
  MFloat c_cellLengthAtCell(const MInt gCellId) const { return grid().cellLengthAtLevel(a_level(gCellId)); }
  MFloat cellVolumeAtCell(const MInt gCellId) const { return grid().cellVolumeAtLevel(a_level(gCellId)); }

  MFloat c_cellLengthAtLevel(const MInt level) const { return grid().cellLengthAtLevel(level); }

  // This also exists in the FV-Solver...
  static constexpr MInt s_maxNoEmbeddedBodies = 20;
  static constexpr const MInt m_noDirs = 2 * nDim;

  MInt noInternalCells() const override { return grid().noInternalCells(); }

  void extendVelocity(const MInt set);
  //------------------------------Collecor-Accessors-----------------------------------------------

  MInt a_noCells() const { return m_cells.size(); }
  // this appends a new cell at the end
  void a_appendCollector() {
    m_cells.append();
  }

  /// \brief Returns property \p p of the cell \p cellId
  void a_resetPropertiesSolver(const MInt cellId) { m_cells.resetProperties(cellId); }

  maia::ls::cell::BitsetType::reference a_isBndryCellG(const MInt cellId) {
    ENSURE_VALID_GCELLID(cellId);
    return m_cells.isBndryG(cellId);
  }
  MBool a_isBndryCellG(const MInt cellId) const {
    ENSURE_VALID_GCELLID(cellId);
    return m_cells.isBndryG(cellId);
  }

  MInt a_maxGCellLevel(const MInt set = -1) const {
    if(!m_gCellLevelJump) {
      return m_maxGCellLevel[0];
    } else {
      ENSURE_VALID_SET(set);
      return m_maxGCellLevel[set];
    }
  }
  /// \brief Returns fExt of the cell \p cellId for index \p n
  MFloat& a_extensionVelocityG(const MInt cellId, const MInt dim, const MInt set) {
    ENSURE_VALID_DIM(dim);
    ENSURE_VALID_SET(set);
    ENSURE_VALID_GCELLID(cellId);
    return m_cells.fExt(cellId, dim, set);
  }
  /// \brief Returns fExt of the cell \p cellId for index \p n
  MFloat a_extensionVelocityG(const MInt cellId, const MInt dim, const MInt set) const {
    ENSURE_VALID_DIM(dim);
    ENSURE_VALID_SET(set);
    ENSURE_VALID_GCELLID(cellId);
    return m_cells.fExt(cellId, dim, set);
  }

  maia::ls::cell::BitsetTypeSet::reference a_inBandG(const MInt cellId, const MInt set) {
    ENSURE_VALID_GCELLID(cellId);
    ENSURE_VALID_SET(set);
    return m_cells.inBand(cellId, set);
  }
  MBool a_inBandG(const MInt cellId, const MInt set) const {
    ENSURE_VALID_GCELLID(cellId);
    ENSURE_VALID_SET(set);
    return m_cells.inBand(cellId, set);
  }

  maia::ls::cell::BitsetTypeSet::reference a_isGBoundaryCellG(const MInt cellId, const MInt set) {
    ENSURE_VALID_GCELLID(cellId);
    ENSURE_VALID_SET(set);
    return m_cells.isGBndryCell(cellId, set);
  }
  MBool a_isGBoundaryCellG(const MInt cellId, const MInt set) const {
    ENSURE_VALID_GCELLID(cellId);
    ENSURE_VALID_SET(set);
    return m_cells.isGBndryCell(cellId, set);
  }

  MBool a_isGZeroCell(const MInt cellId, const MInt set) const {
    ENSURE_VALID_GCELLID(cellId);
    ENSURE_VALID_SET(set);
    return m_cells.isGZero(cellId, set);
  }

  maia::ls::cell::BitsetTypeSet::reference a_isGZeroCell(const MInt cellId, const MInt set) {
    ENSURE_VALID_GCELLID(cellId);
    ENSURE_VALID_SET(set);
    return m_cells.isGZero(cellId, set);
  }

  MBool a_wasGZeroCell(const MInt cellId, const MInt set) const {
    ENSURE_VALID_GCELLID(cellId);
    ENSURE_VALID_SET(set);
    return m_cells.wasGZero(cellId, set);
  }

  maia::ls::cell::BitsetTypeSet::reference a_wasGZeroCell(const MInt cellId, const MInt set) {
    ENSURE_VALID_GCELLID(cellId);
    ENSURE_VALID_SET(set);
    return m_cells.wasGZero(cellId, set);
  }

  /// \brief Returns the hasPositiveSign\p cellId for the set\p set
  maia::ls::cell::BitsetTypeSet::reference a_hasPositiveSign(const MInt cellId, const MInt set) {
    ENSURE_VALID_GCELLID(cellId);
    ENSURE_VALID_SET(set);
    return m_cells.hasPositiveSign(cellId, set);
  }

  /// \brief Returns the hasPositiveSign\p cellId for the set\p set
  MBool a_hasPositiveSign(const MInt cellId, const MInt set) const {
    ENSURE_VALID_GCELLID(cellId);
    ENSURE_VALID_SET(set);
    return m_cells.hasPositiveSign(cellId, set);
  }
  /// \brief Returns the signed (MInt) version of hasPositiveSign\p cellId for the set\p set
  MInt a_levelSetSign(const MInt cellId, const MInt set) { return a_hasPositiveSign(cellId, set) ? 1 : -1; }

  /// \brief Returns normalVector of the cell \p cellId for index \p n
  MFloat& a_normalVectorG(const MInt cellId, const MInt dim, const MInt set) {
    ENSURE_VALID_DIM(dim);
    ENSURE_VALID_SET(set);
    ENSURE_VALID_GCELLID(cellId);
    return m_cells.normalVector(cellId, dim, set);
  }
  /// \brief Returns normalVector of the cell \p cellId for index \p n
  MFloat a_normalVectorG(const MInt cellId, const MInt dim, const MInt set) const {
    ENSURE_VALID_DIM(dim);
    ENSURE_VALID_SET(set);
    ENSURE_VALID_GCELLID(cellId);
    return m_cells.normalVector(cellId, dim, set);
  }

  maia::ls::cell::BitsetType::reference a_nearGapG(const MInt cellId) {
    ENSURE_VALID_GCELLID(cellId);
    return m_cells.nearGap(cellId);
  }
  MBool a_nearGapG(const MInt cellId) const {
    ENSURE_VALID_GCELLID(cellId);
    return m_cells.nearGap(cellId);
  }

  maia::ls::cell::BitsetType::reference a_regridTriggerG(const MInt cellId) {
    ENSURE_VALID_GCELLID(cellId);
    return m_cells.regridTrigger(cellId);
  }
  MBool a_regridTriggerG(const MInt cellId) const {
    ENSURE_VALID_GCELLID(cellId);
    return m_cells.regridTrigger(cellId);
  }
  /// \brief Returns bodyId of the cell \p cellId for set \p set
  MInt& a_bodyIdG(const MInt cellId, const MInt set) {
    ENSURE_VALID_SET(set);
    ENSURE_VALID_GCELLID(cellId);
    return m_cells.bodyId(cellId, set);
  }
  /// \brief Returns bodyId of the cell \p cellId for set \p set
  MInt a_bodyIdG(const MInt cellId, const MInt set) const {
    ENSURE_VALID_SET(set);
    ENSURE_VALID_GCELLID(cellId);
    return m_cells.bodyId(cellId, set);
  }

  /// \brief Returns secondBodyId of the cell \p cellId for set \p set
  MInt& a_secondBodyId(const MInt cellId) {
    ENSURE_VALID_GCELLID(cellId);
    return m_cells.secondBodyId(cellId);
  }
  /// \brief Returns secondBodyId of the cell \p cellId for set \p set
  MInt a_secondBodyId(const MInt cellId) const {
    ENSURE_VALID_GCELLID(cellId);
    return m_cells.secondBodyId(cellId);
  }

  /// \brief Returns curvature of the cell \p cellId for set \p set
  MFloat& a_curvatureG(const MInt cellId, const MInt set) {
    ENSURE_VALID_SET(set);
    ENSURE_VALID_GCELLID(cellId);
    return m_cells.curvature(cellId, set);
  }
  /// \brief Returns curvature of the cell \p cellId for set \p set
  MFloat a_curvatureG(const MInt cellId, const MInt set) const {
    ENSURE_VALID_SET(set);
    ENSURE_VALID_GCELLID(cellId);
    return m_cells.curvature(cellId, set);
  }

  /// \brief Returns levelSetFunction of the cell \p cellId
  MFloat& a_levelSetFunctionG(const MInt cellId, const MInt set) {
    ENSURE_VALID_SET(set);
    ENSURE_VALID_GCELLID(cellId);
    return m_cells.gFunction(cellId, set);
  }
  /// \brief Returns levelSetFunction of the cell \p cellId
  MFloat a_levelSetFunctionG(const MInt cellId, const MInt set) const {
    ENSURE_VALID_SET(set);
    ENSURE_VALID_GCELLID(cellId);
    return m_cells.gFunction(cellId, set);
  }

  /// \brief Returns the old levelSetFunction of the cell \p cellId
  MFloat& a_oldLevelSetFunctionG(const MInt cellId, const MInt set) {
    ENSURE_VALID_SET(set);
    ENSURE_VALID_GCELLID(cellId);
    // return m_oldLevelSetFunction[ IDX_LSSET( cellId, set) ];
    return m_cells.oldGFunction(cellId, set);
  }
  /// \brief Returns the old levelSetFunction of the cell \p cellId
  MFloat a_oldLevelSetFunctionG(const MInt cellId, const MInt set) const {
    ENSURE_VALID_SET(set);
    ENSURE_VALID_GCELLID(cellId);
    // return m_oldLevelSetFunction[ IDX_LSSET( cellId, set) ];
    return m_cells.oldGFunction(cellId, set);
  }

  /// \brief Returns ls-FunctionSlope of the cell \p cellId for set \p dim \p set
  MFloat& a_levelSetFunctionSlope(const MInt cellId, const MInt dim, const MInt set) {
    ENSURE_VALID_DIM(dim);
    ENSURE_VALID_SET(set);
    ENSURE_VALID_GCELLID(cellId);
    return m_cells.levelSetFunctionSlope(cellId, dim, set);
  }
  /// \brief Returns ls-FunctionSlope of the cell \p cellId for set \p dim \p set
  MFloat a_levelSetFunctionSlope(const MInt cellId, const MInt dim, const MInt set) const {
    ENSURE_VALID_DIM(dim);
    ENSURE_VALID_SET(set);
    ENSURE_VALID_GCELLID(cellId);
    return m_cells.levelSetFunctionSlope(cellId, dim, set);
  }

  /// \brief Returns ls-RHS of the cell \p cellId for set \p set
  MFloat& a_levelSetRHS(const MInt cellId, const MInt set) {
    ENSURE_VALID_SET(set);
    ENSURE_VALID_GCELLID(cellId);
    return m_cells.levelSetRHS(cellId, set);
  }
  /// \brief Returns ls-RHS of the cell \p cellId for set \p set
  MFloat a_levelSetRHS(const MInt cellId, const MInt set) const {
    ENSURE_VALID_SET(set);
    ENSURE_VALID_GCELLID(cellId);
    return m_cells.levelSetRHS(cellId, set);
  }

  /// \brief Returns corrected burning velocity of the cell \p cellId for set \p set
  MFloat& a_correctedBurningVelocity(const MInt cellId, const MInt set) {
    ENSURE_VALID_SET(set);
    ENSURE_VALID_GCELLID(cellId);
    return m_cells.correctedBurningVelocity(cellId, set);
  }
  /// \brief Returns corrected burning velocity of the cell \p cellId for set \p set
  MFloat correctedBurningVelocity(const MInt cellId, const MInt set) const {
    ENSURE_VALID_SET(set);
    ENSURE_VALID_GCELLID(cellId);
    return m_cells.correctedBurningVelocity(cellId, set);
  }

  /// \brief Returns the containing cell\p cellId
  MLong& a_containingCell(const MInt cellId, const MInt body) { return m_cells.containingCell(cellId, body); }
  /// \brief Returns the containing cell\p cellId
  MLong a_containingCell(const MInt cellId, const MInt body) const { return m_cells.containingCell(cellId, body); }

  /// \brief Returns the containing Domain\p cellId
  MInt& a_containingDomain(const MInt cellId, const MInt body) { return m_cells.containingDomain(cellId, body); }
  /// \brief Returns the containing Domain\p cellId
  MInt a_containingDomain(const MInt cellId, const MInt body) const { return m_cells.containingDomain(cellId, body); }

  /// \brief Returns the gap width\p cellId
  MFloat& a_gapWidth(const MInt id) { return m_cells.gapWidth(id); }
  /// \brief Returns the gap width\p cellId
  MFloat a_gapWidth(const MInt id) const { return m_cells.gapWidth(id); }

  /// \brief Returns the potential gap cell\p cellId
  MInt& a_potentialGapCell(const MInt id) { return m_cells.potentialGapCell(id); }
  /// \brief Returns the potential gap cell\p cellId
  MInt a_potentialGapCell(const MInt id) const { return m_cells.potentialGapCell(id); }

  /// \brief Returns the potential gap cell close\p cellId
  MInt& a_potentialGapCellClose(const MInt id) { return m_cells.potentialGapCellClose(id); }
  /// \brief Returns the potential gap cell close\p cellId
  MInt a_potentialGapCellClose(const MInt id) const { return m_cells.potentialGapCellClose(id); }

  MFloat& a_meanCoord(const MInt dir) { return m_meanCoord[dir]; }

  MFloat a_meanCoord(const MInt dir) const { return m_meanCoord[dir]; }

  //------------------------------Grid-Accessors-----------------------------------------------

  // helper accessor, needed in the cartesiansolver
  // lsSolver does not have solver additional cells as the fvSolver(ghost- or splitCells)
  MInt c_noCells() const { return a_noCells(); }

  // a_parentGId
  MInt c_parentId(const MInt gCellId) const {
    ENSURE_VALID_GCELLID(gCellId);
    return grid().tree().parent(gCellId);
  }

  // a_noChildrenG
  MInt c_noChildren(const MInt gCellId) const {
    ENSURE_VALID_GCELLID(gCellId);
    return grid().tree().noChildren(gCellId);
  }

  // a_childGId
  MInt c_childId(const MInt gCellId, const MInt pos) const {
    ENSURE_VALID_GCELLID(gCellId);
    return grid().tree().child(gCellId, pos);
  }

  // a_coordinateG
  /// \brief Returns the coordinate of the cell \p cellId for direction \p dim
  MFloat c_coordinate(const MInt gCellId, const MInt dim) const {
    ENSURE_VALID_DIM(dim);
    ENSURE_VALID_GCELLID(gCellId);
    return grid().tree().coordinate(gCellId, dim);
  }

  // a_neighborGId
  MInt c_neighborId(const MInt gCellId, const MInt dir) const {
    ENSURE_VALID_DIR(dir);
    ENSURE_VALID_GCELLID(gCellId);

    if(grid().tree().neighbor(gCellId, dir) > -1) {
      return grid().tree().neighbor(gCellId, dir);
    } else if(m_levelSetBC == "SYMMETRIC") {
      return gCellId;
    } else if(m_levelSetBC == "NONE") {
      return -1;
    } else {
      mTerm(1, AT_, "Return of the neighborId should have happend by now!");
      return gCellId;
    }
  }

  /// \brief Returns bandNeighborId of the cell \p cellId for index \p n
  MInt a_bandNghbrIdsG(const MInt cellId, const MInt dir, const MInt set) const {
    ENSURE_VALID_DIR(dir);
    ENSURE_VALID_SET(set);
    ENSURE_VALID_GCELLID(cellId);

    if(a_hasNeighbor(cellId, dir) && a_inBandG(c_neighborId(cellId, dir), set)) {
      return c_neighborId(cellId, dir);
    } else {
      return cellId;
    }
  }

  // a_hasNeighborG
  /// \brief Returns noNeighborIds of the gcell \p CellId variables \p varId
  MInt a_hasNeighbor(const MInt gCellId, const MInt dir) const {
    ENSURE_VALID_DIR(dir);
    ENSURE_VALID_GCELLID(gCellId);
    if(gCellId > a_noCells() - 1) return -1;
    return grid().tree().hasNeighbor(gCellId, dir);
  }

  /// \brief Returns IsWindow of the cell \p cellId
  MBool a_isWindow(const MInt gCellId) const {
    ENSURE_VALID_GCELLID(gCellId);
    return m_cells.hasProperty(gCellId, LsCell::IsWindow);
  }
  /// \brief Returns IsWindow of the cell \p cellId
  maia::ls::cell::BitsetType::reference a_isWindow(const MInt gCellId) {
    ENSURE_VALID_GCELLID(gCellId);
    return m_cells.hasProperty(gCellId, LsCell::IsWindow);
  }

  MLong c_globalId(const MInt gCellId) const {
    ENSURE_VALID_GCELLID(gCellId);
    return grid().tree().globalId(gCellId);
  }

  MInt a_domainId(const MLong gGlobalId) { return grid().findNeighborDomainId(gGlobalId); }

  MInt a_localId(const MLong gGlobalId) { return grid().globalToLocalId(gGlobalId); }

  MBool c_isLeafCell(const MInt gCellId) const {
    ENSURE_VALID_GCELLID(gCellId);
    return grid().tree().isLeafCell(gCellId);
  }

  /// \brief Returns IsHalo of the cell \p cellId
  MBool a_isHalo(const MInt gCellId) const {
    ENSURE_VALID_GCELLID(gCellId);
    return m_cells.hasProperty(gCellId, LsCell::IsHalo);
  }
  /// \brief Returns IsHalo of the cell \p cellId
  maia::ls::cell::BitsetType::reference a_isHalo(const MInt gCellId) {
    ENSURE_VALID_GCELLID(gCellId);
    return m_cells.hasProperty(gCellId, LsCell::IsHalo);
  }

  // a_levelG
  /// \brief Returns the level of the gcell \p gCellId
  MInt a_level(const MInt gCellId) const {
    ENSURE_VALID_GCELLID(gCellId);
    return grid().tree().level(gCellId);
  }

  //------------------------------Other-Accessors-----------------------------------------------

  /// \brief Returns flameSpeed of the cell \p cellId for index \p n
  MFloat& a_flameSpeedG(const MInt cellId, const MInt set) {
    ENSURE_VALID_SET(set);
    ENSURE_VALID_GCELLID(cellId);
    return m_flameSpeed;
  }
  /// \brief Returns flameSpeed of the cell \p cellId for index \p n
  MFloat a_flameSpeedG(const MInt cellId, const MInt set) const {
    ENSURE_VALID_GCELLID(cellId);
    ENSURE_VALID_SET(set);
    return m_flameSpeed;
  }

  /// \brief Returns the levelSet-Adaptation-forcing
  MBool forceAdaptation() override {
    // If lsSolver is inactive on one of the ranks we need this allreduce!
    if(grid().hasInactiveRanks()) {
      this->startLoadTimer(AT_);
      MPI_Allreduce(MPI_IN_PLACE, &m_forceAdaptation, 1, MPI_C_BOOL, MPI_LOR, grid().raw().mpiComm(), AT_,
                    "MPI_IN_PLACE", "m_forceAdaptation");
      this->stopLoadTimer(AT_);
    }
    return m_forceAdaptation;
  }

  constexpr MInt a_bandCellId(MInt id, MInt set) const { return m_bandCells[set][id]; }
  constexpr MInt a_internalBandCellId(MInt id, MInt set) const { return m_internalBandCells[set][id]; }
  constexpr MInt a_bandBndryCellId(MInt id, MInt set) const { return m_bandBndryCells[set][id]; }
  constexpr MInt a_G0CellId(MInt id, MInt set) const { return m_G0Cells[set][id]; }
  constexpr MInt a_gBndryCellId(MInt id, MInt set) const { return m_gBndryCells[set][id]; }

  constexpr MInt a_noBandCells(MInt set) const { return m_bandCells[set].size(); }
  constexpr MInt a_noInternalBandCells(MInt set) const { return m_internalBandCells[set].size(); }
  constexpr MInt a_noBandBndryCells(MInt set) const { return m_bandBndryCells[set].size(); }
  constexpr MInt a_noG0Cells(MInt set) const { return m_G0Cells[set].size(); }
  constexpr MInt a_noGBndryCells(MInt set) const { return m_gBndryCells[set].size(); }

  MInt& a_bandLayer(MInt id, MInt set) { return m_bandLayer[IDX_LSSET(id, set)]; }
  MInt a_bandLayer(MInt id, MInt set) const { return m_bandLayer[IDX_LSSET(id, set)]; }

  MInt& a_internalBandLayer(MInt id, MInt set) { return m_internalBandLayer[IDX_LSSET(id, set)]; }
  MInt a_internalBandLayer(MInt id, MInt set) const { return m_internalBandLayer[IDX_LSSET(id, set)]; }

  MInt getCurrentTimeStep() const override { return globalTimeStep; }

  void resetExtensionVelocity() { m_cells.resetExtensionVelocity(); }

  //------------------------- Data -------------------------------------------

  // semi-Lagrange levelSet-related
  MBool m_semiLagrange;
  MFloat* m_semiLagrange_xShift_ref = nullptr;
  LsControlPoint<nDim> m_gCtrlPnt;
  MBool m_static_semiLagrangeTimeStep_firstTime = true;

  // ls-Collectore Type mode:
  MInt m_lsCollectorMode = -1;

  // rotating levelSet-related
  MBool m_LsRotate = false;
  MBool m_reconstructOldG = true;
  MInt m_rotatingReinitTrigger = 0;
  MFloat* m_semiLagrange_xRot_ref = nullptr;
  MFloat* m_semiLagrange_xRot_STL = nullptr;
  MInt* m_globalSndOffsets = nullptr;
  MInt* m_globalRcvOffsets = nullptr;
  MInt* m_initialGCell = nullptr;
  std::vector<MInt> m_newCells;
  std::map<MInt, MInt> m_swapIds;
  std::vector<MInt> m_bodiesToCompute;
  MInt m_noBodiesToCompute;
  MInt* m_cellDomIds = nullptr;
  MFloat* m_bodyAngularVelocity = nullptr;
  MFloat* m_bodyAngularAcceleration = nullptr;
  MFloat* m_omega = nullptr;
  MFloat* m_bodyRadius = nullptr;

  MFloat m_referenceLength;
  MFloat m_referenceVelocity;

  MBool m_periodicMovement;
  MInt m_periodicDirection;
  MFloat m_periodicDistance;

  MFloat m_meanCoord[3]{};

  // adaptation-related
  std::map<MInt, MInt> m_refinedCells;
  MInt m_reconstructBand = 0;
  MBool m_newRefinedBand = false;
  MInt m_adaptationLevel;
  MBool m_forceAdaptation = false;
  MBool m_initialRefinement = false;
  MBool m_adaptationSinceLastRestart = false;
  MBool m_refineDiagonals = false;

  std::map<MInt, MInt> m_oldG0Cells;
  MBool* m_geometryChange = nullptr;

  // general levelSet properties
  MInt m_maxNoCells;
  MInt m_noSets;
  MInt m_startSet;

  MBool m_buildCollectedLevelSetFunction;
  MBool m_determineG0CellsMode;
  MInt m_noBodyBndryCndIds;
  MInt* m_levelSetSign = nullptr;
  MBool* m_computeSet = nullptr;
  MBool* m_computeSet_tmp = nullptr;
  MBool* m_computeSet_backup = nullptr;
  MBool* m_changedSet = nullptr;
  MBool m_GFieldInitFromSTL;
  MBool m_GFieldFromSTLInitCheck;
  MInt m_STLReinitMode;
  MBool m_GWithReConstruction;
  MInt* m_bodyToSetTable = nullptr;
  MInt* m_noBodiesInSet = nullptr;
  MInt** m_setToBodiesTable = nullptr;
  MInt* m_bodyBndryCndIds = nullptr;
  MInt m_noInitGFieldFromSTLBndCndIds{};
  MInt* m_initGFieldFromSTLBndCndIds{};
  MBool m_highOrderDeltaFunction;
  MBool m_fourthOrderNormalCurvatureComputation;
  MBool m_curvatureDamp;
  MFloat m_curvatureDampFactor;
  MBool m_useLocalMarksteinLength;
  MBool m_sharpDamp;
  MBool m_hyperbolicCurvature;
  MInt m_gRKMethod;
  MString m_levelSetDiscretizationScheme;
  MInt m_gBandWidth;
  MInt m_gShadowWidth;
  MInt m_gShadowWidthRans;
  MInt m_gInnerBound;

  MBool m_gCellLevelJump = false;
  MFloat m_gCellDistance;
  MFloat m_FgCellDistance;
  MInt m_computeExtVel{};
  MBool m_smoothExtVel;
  MInt m_extVelIterations;
  MFloat m_extVelConvergence;
  MFloat m_extVelCFL;
  MString m_reinitMethod;
  MInt m_gReinitIterations;
  MInt m_minReinitializationSteps;
  MInt m_maintenanceIterations;
  MBool m_guaranteeReinit;
  MFloat m_reinitCFL;
  MInt m_intermediateReinitIterations;
  MFloat m_reinitConvergence;
  MFloat m_reinitConvergenceReset;
  MFloat m_reinitThreshold;
  MFloat m_reinitThresholdAvg;
  MFloat m_omegaReinit;
  MFloat m_relaxationFactor;
  MInt m_levelSetTestCase;
  MInt m_levelSetBoundaryCondition;
  MString m_levelSetBC;
  MInt m_noHaloLayers;
  MInt m_reinitInterval;
  MBool m_maintainOuterBandLayers;
  MBool m_writeReinitializationStatistics;
  MBool m_interpolateFlowFieldToFlameFront;
  MBool m_writeOutAllLevelSetFunctions;
  MBool m_writeOutAllExtensionVelocities;
  MBool m_writeOutAllCurvatures;
  MBool m_writeOutAllCorrectedBurningVelocity;
  MBool m_writeOutAllFlameSpeeds;
  MBool m_writeOutAllNormalVectors;
  MInt* m_outerBandWidth = nullptr;

  // level-set variables
  MInt* m_cellList = nullptr;
  std::vector<MInt>* m_bandCells = nullptr;
  std::vector<MInt>* m_internalBandCells = nullptr;
  MInt* m_bandLayer = nullptr;
  MInt* m_internalBandLayer = nullptr;
  std::vector<MInt>* m_bandBndryCells = nullptr;
  std::vector<MInt>* m_G0Cells = nullptr;
  std::vector<MInt>* m_gBndryCells = nullptr;
  MInt m_gRKStep;
  MInt m_nogRKSteps;
  MFloat* m_gRKalpha = nullptr;
  MFloat m_outsideGValue;

  MInt** m_phiRatioCells = nullptr;
  MInt m_GCtrlPntMethod;
  MFloat* m_signG = nullptr;
  MFloat** m_phiRatio = nullptr;
  MFloat* m_correction = nullptr;
  MFloat* m_d = nullptr;

  // exchange buffers:
  MInt** m_intSendBuffers = nullptr;
  MInt** m_intReceiveBuffers = nullptr;
  MFloat** m_gSendBuffers = nullptr;
  MFloat** m_gReceiveBuffers = nullptr;

  MFloat* m_hypTanLSF = nullptr;

  // More level-set related properties. Most of them are used in FV Solver as well and exist in both solvers right now
  //
  MBool m_trackMovingBndry{};
  MInt m_trackMbStart{};
  MInt m_trackMbEnd{};
  MBool m_constructGField{};

  // levelSet-Solver Types:
  MBool m_levelSetMb;
  MBool m_levelSetRans;
  MBool m_levelSetLb;
  MBool m_combustion;
  MBool m_freeSurface;
  MBool m_LSSolver;
  MBool m_levelSetFv;

  MInt m_noEmbeddedBodies;

  // combustion-related properties
  MFloat* m_maxFlameFrontPosition = nullptr;
  MFloat* m_minFlameFrontPosition = nullptr;
  MFloat* m_meanFlameFrontPosition = nullptr;
  MFloat m_steadyFlameAngle;
  MFloat m_steadyFlameLength;
  MFloat m_yOffsetFlameTube;
  MFloat m_yOffsetFlameTube2;
  MFloat m_radiusFlameTube;
  MFloat m_radiusFlameTube2;
  MFloat m_xOffsetFlameTube;
  MFloat m_xOffsetFlameTube2;
  MFloat m_marksteinLength;
  MFloat m_marksteinLengthPercentage;
  MFloat* m_localMarksteinLength = nullptr;
  MFloat m_realRadiusFlameTube;
  MBool m_forcing;
  // Somewhere in the code it says this is 'ghost fluid method' related..
  // maybe this can be removed along with functions where this is used.
  MFloat m_flameRadiusOffset;
  MBool m_twoFlames;
  MFloat m_dampingDistanceFlameBase;
  MFloat m_dampingDistanceFlameBaseExtVel;
  MFloat m_noReactionCells;
  MFloat m_jetHalfWidth;
  MFloat m_jetHalfLength;
  MFloat m_initialFlameHeight;
  MBool m_filterFlameTubeEdges;
  MFloat m_filterFlameTubeEdgesDistance;
  MBool m_useCorrectedBurningVelocity;
  MFloat m_flameSpeed{};
  MFloat m_massConsumption;
  MFloat m_arcLength;
  MBool m_plenum;
  MFloat m_rhoFlameTube;
  MFloat m_rhoInfinity;
  MFloat* m_oldHypTanLSF = nullptr;
  MFloat* m_hypTanLSFRHS = nullptr;


  MString m_currentGridFileName;
  MString m_currentFileName = "";
  MBool m_bodyIdOutput;
  MBool m_maxLevelChange = false;

  MInt m_initialCondition;

  // balance
  MInt m_loadBalancingReinitStage;

  // time related
  MInt m_timeStepMethod;
  MFloat m_cfl;
  MFloat m_time = NAN;
  MFloat m_timeStep = NAN;

  // ... identifyBodies
  MBool m_static_identifyBodies_first = true;
  MFloat m_static_identifyBodies_initialInsidePoints[s_maxNoEmbeddedBodies * 3];
  MFloat m_static_identifyBodies_shiftTime;

  // gap-Handling:
  MInt m_noGapCells;
  MInt m_noOldGapCells;
  MBool m_closeGaps = false;
  MInt m_forceNoGaps;
  MInt m_gapInitMethod;
  MInt m_noGapRegions = -1;
  std::vector<MFloat> m_minGapWidth;
  std::vector<MFloat> m_minGapWidthDt1;
  MString m_gapReinitMethod;
  MFloat m_gapDeltaMin = 1;
  MFloat m_gapDeltaMinOrig = 1;
  MFloat* m_gapAngleClose = nullptr;
  MFloat* m_gapAngleOpen = nullptr;
  MFloat* m_gapSign = nullptr;

  MBool m_virtualSurgery = false;
  MFloat m_sphereRadiusLimit = 5.0;
  MFloat** m_correctedDistances = nullptr;
  MInt* m_cellIsInDiffRegion = nullptr;
  MInt m_noInterpolationRegions = 0;
  MInt m_approxNoInterpReg = 0;
  MInt* m_interpStartTime = nullptr;
  MInt* m_noInterpTimeSteps = nullptr;

  MBool m_engineSetup = false;

  MInt m_G0regionId = -1;

 private:
  MInt* m_maxGCellLevel = nullptr;

  // ... computeBodyProperties
  MBool m_static_computeBodyProperties_first = true;
  MFloat m_static_computeBodyProperties_amplitude[s_maxNoEmbeddedBodies]{};
  MFloat m_static_computeBodyProperties_freqFactor[s_maxNoEmbeddedBodies]{};
  MFloat m_static_computeBodyProperties_initialBodyCenter[s_maxNoEmbeddedBodies * 3]{};
  MFloat m_static_computeBodyProperties_Strouhal{};
  MFloat m_static_computeBodyProperties_mu[s_maxNoEmbeddedBodies]{};
  MFloat m_static_computeBodyProperties_mu2[s_maxNoEmbeddedBodies]{};
  MFloat m_static_computeBodyProperties_liftStartAngle1[s_maxNoEmbeddedBodies]{};
  MFloat m_static_computeBodyProperties_liftEndAngle1[s_maxNoEmbeddedBodies]{};
  MFloat m_static_computeBodyProperties_liftStartAngle2[s_maxNoEmbeddedBodies]{};
  MFloat m_static_computeBodyProperties_liftEndAngle2[s_maxNoEmbeddedBodies]{};
  MFloat m_static_computeBodyProperties_circleStartAngle[s_maxNoEmbeddedBodies]{};
  MFloat m_static_computeBodyProperties_normal[s_maxNoEmbeddedBodies * 3]{};
  MInt m_static_computeBodyProperties_bodyToFunction[s_maxNoEmbeddedBodies]{};
  MFloat m_static_computeBodyProperties_omega{};
  MFloat m_static_computeBodyProperties_rotAngle{};
  MFloat m_static_computeBodyProperties_temperature[s_maxNoEmbeddedBodies]{};

  std::set<std::pair<MFloat, MFloat>>* m_forcedMotionInput = nullptr;

  // ... setUpPotentialGapCells
  MBool m_static_setUpPotentialGapCells_first = true;
  MFloat m_static_setUpPotentialGapCells_normal[s_maxNoEmbeddedBodies * 3];
  MFloat m_static_setUpPotentialGapCells_center[s_maxNoEmbeddedBodies * 3];
  MFloat m_static_setUpPotentialGapCells_radius[s_maxNoEmbeddedBodies];
  MFloat m_static_setUpPotentialGapCells_height[s_maxNoEmbeddedBodies];
  MFloat m_static_setUpPotentialGapCells_normalClose[s_maxNoEmbeddedBodies * 3];
  MFloat m_static_setUpPotentialGapCells_centerClose[s_maxNoEmbeddedBodies * 3];
  MFloat m_static_setUpPotentialGapCells_radiusClose[s_maxNoEmbeddedBodies];
  MFloat m_static_setUpPotentialGapCells_heightClose[s_maxNoEmbeddedBodies];
  MInt m_static_setUpPotentialGapCells_bodyClose[s_maxNoEmbeddedBodies];
  MInt m_static_setUpPotentialGapCells_noGapRegionsClose = 1;


  // crankAngle
  MFloat m_static_crankAngle_Strouhal = -99;
  MFloat m_static_crankAngle_initialCad = -99;

  MFloat m_weightBaseCell = 1.0;
  MFloat m_weightLeafCell = 1.0;
  MFloat m_weightBandCell = 1.0;
  MFloat m_weightMulitSolverFactor = 1.0;
  MBool m_limitWeights = false;
  MBool m_weightLevelSet = true;

  MBool m_firstSolutionExchange = false;

  //------------------ functions ----------------------------------------------
 public:
  void initSolver() override;
  void initLocalizedLevelSetCG();
  void createBaseGgridCG();
  void createGgridCG(MBool = false);
  void generateListOfGExchangeCellsCG();
  void buildLevelSetTube(MInt mode = -1);
  void fastBuildLevelSetTubeCG();
  void testCellsCG();
  void restartLocalizedLevelSetCG();
  MInt loadLevelSetGridFlowVarsParCG(const MChar* fileName);
  virtual void writeRestartLevelSetFileCG(MBool, const MString&, const MString&);

  // grid controller related mesh adaptation functions
  void prepareAdaptation() override;
  void setSensors(std::vector<std::vector<MFloat>>& sensors,
                  std::vector<MFloat>& sensorWeight,
                  std::vector<std::bitset<64>>& sensorCellFlag,
                  std::vector<MInt>& sensorSolverId) override;
  void postAdaptation() override;
  void finalizeAdaptation() override;

  // Sensors
  void sensorInterface(std::vector<std::vector<MFloat>>& sensors, std::vector<std::bitset<64>>& sensorCellFlag,
                       std::vector<MFloat>& sensorWeight, MInt sensorOffset, MInt sen) override;

  void removeChilds(const MInt) override;
  void removeCell(const MInt) override;
  void refineCell(const MInt) override;
  void swapCells(const MInt, const MInt) override;
  void swapProxy(const MInt cellId0, const MInt cellId1) override;
  MInt cellOutside(const MFloat*, const MInt, const MInt) override;
  void exchangeAllLevelSetData();
  void exchangeLevelSet();
  void exchangeLs(MFloat*, MInt, MInt);
  void exchangeGapInfo();
  void initAzimuthalExchange();
  void resizeGridMap() override;
  MBool levelSetAdaptationTrigger();
  void getContainingCellFromNeighbor(MInt body, MInt cellId, MFloat* xCoord, MFloat* xOld);
  MInt getContainingCellHalo(MFloat* point);
  void finalizeInitSolver() override;
  void initRotatingLS();
  void resetContainingGCells();
  void updateContainingGCells(MInt mode = 0);
  void copyWindowToHaloIds();
  void checkHaloCells();
  void setInterfaceList(MIntScratchSpace& interfaceCells);
  // balancing:
  void resetSolverFull();
  void resetSolver() override;
  void balance(const MInt* const noCellsToReceiveByDomain, const MInt* const noCellsToSendByDomain,
               const MInt* const targetDomainsByCell, const MInt oldNoCells) override;
  void balancePre() override;
  void balancePost() override;
  void finalizeBalance() override;


  MBool hasSplitBalancing() const override { return true; }
  void localToGlobalIds() override;
  MInt noCellDataDlb() const override {
    if(grid().isActive()) {
      // a_levelSetFunctionG(cellId, set)
      // a_regridTriggerG(cellId)
      MInt cellData = 2;
      if(m_semiLagrange) {
        // a_oldLevelSetFunctionG(cellId, set)
        // m_initialGCell[cellId]
        // m_containingCell[b*m_maxNoGCells+cellId]
        cellData += 1 + 2 * (MInt)(!m_reconstructOldG && m_LsRotate);
      } else {
        // a_curvatureG(cellId, 0)
        // a_normalVectorG(cellId, dim, 0)
        // a_levelSetFunctionSlope(cellId, dim, 0)
        cellData += 3;
      }
      return cellData;
    } else {
      return 0;
    }
  };

  MInt cellDataTypeDlb(const MInt dataId) const override {
    if(dataId == 0) {
      return MFLOAT;
    } else if(dataId == 1) {
      return MINT;
    } else if(dataId < 5) {
      if(m_semiLagrange) {
        if(dataId == 2)
          return MFLOAT;
        else if(dataId > 2) {
          if(!m_reconstructOldG && m_LsRotate) {
            return MINT;
          } else
            TERMM(1, "solverCelldataType: invalid data id for !m_reconstructOldG && m_LsRotate");
        }
      } else {
        return MFLOAT;
      }
    } else
      TERMM(1, "solverCelldataType: invalid data id");
    // This should not be reached, just for compiler
    return 0;
  };

  MInt cellDataSizeDlb(const MInt dataId, const MInt gridCellId) override;

  /// Return solver data for DLB
  void getCellDataDlb(const MInt dataId, const MInt oldNoCells, const MInt* const bufferIdToCellId,
                      MInt* const data) override;
  void getCellDataDlb(const MInt dataId, const MInt oldNoCells, const MInt* const bufferIdToCellId,
                      MFloat* const data) override;

  /// Set solver data for DLB
  void setCellDataDlb(const MInt dataId, const MInt* const data) override;
  void setCellDataDlb(const MInt dataId, const MFloat* const data) override;

  // Partitioning
  void setCellWeights(MFloat*) override;
  MInt noLoadTypes() const override;
  void getDefaultWeights(MFloat* weights, std::vector<MString>& names) const override;
  void getLoadQuantities(MInt* const loadQuantities) const override;
  MFloat getCellLoad(const MInt cellId, const MFloat* const weights) const override;
  void getSolverTimings(std::vector<std::pair<MString, MFloat>>& solverTimings, const MBool allTimings) override;
  void limitWeights(MFloat*) override;
  void getDomainDecompositionInformation(std::vector<std::pair<MString, MInt>>& domainInfo) override;

  void reconstructOldGField();
  void rotateLevelSet(MInt returnMode, MFloat* cellData, MInt body, const MFloat* xCoord, const MFloat* xCenter,
                      const MFloat* angle);
  void processRotatingLevelSet(MFloat& phi, MInt& cellId, MInt& domId, MFloat* point, MInt set);
  MInt checkSecondLayerCells(std::vector<MInt>& diag2Cells, std::map<MInt, std::vector<MInt>>& dirCode, MFloat* point);
  void prepareGlobalComm(MInt* noCellsToDom);
  void globalToLocalIdsContainingCells();
  void localToGlobalIdsContainingCells();

  // functions that used to be in fvlevelsetsolver.h
  //
  void computeLevelSetRHS();
  void levelSetConstrainedReinitialization(MInt methodId, MInt startSet, MInt endSet, MInt gapMode);
  void levelSetHighOrderConstrainedReinitialization(MInt methodId, MInt startSet, MInt endSet, MInt gapMode);
  void maintainOuterBandLayers(MInt order, MInt startSet, MInt endSet);
  void determineMinMaxMeanInterfacePosition();
  void determineSteadyFlameLength();
  void determineMinMaxMeanRegionInterfacePosition(MFloat xRegN, MFloat xRegP, MFloat yRegN, MFloat yRegP, MInt set);
  void reinitBand(MInt startSet, MInt endSet);
  void initializeGControlPoint();
  void initializeIntegrationScheme();
  void initializeIntegrationScheme_semiLagrange();
  void setGCellBndryProperty();
  void applyLevelSetBoundaryConditions();
  void initializeGField();
  void constructGFieldFromSTL(MInt ConstructFlag);
  void rotateSTL(MInt direction);
  void rotateSTL(MInt direction, MInt body, MFloat* center);
  void spatiallyAdaptiveCorrectionFromSTL();
  void levelSetReinitialization(MInt mode = 1);
  //  void fastBuildLevelSetTube();
  void setBandNewArrivals(MInt computingSet = -1);
  void updateLowerGridLevels(MInt mode = -1);
  void updateAllLowerGridLevels(MInt mode = -1);
  void determineG0Cells(MInt computingSet = -1);
  void determineBandCells(MInt mode = -1);
  void updateBndryCellList();
  void resetOutsideCells(MInt mode = -1);
  void resetOldOutsideCells();
  void computeCurvature(MInt mode = -1);
  void computeCurvaturePeriodic();
  void determinePropagationSpeed();
  void computeNormalVectors(MInt mode = -1);
  void computeNormalVectorsPeriodic();
  void computeNormalVectorsAtFront();
  void computeExtensionVelocityGEQUPVMarksteinOpt(MFloat* FfluidDensity, MInt set);
  void computeGCellTimeStep();
  void levelSetRestriction();
  void computeZeroLevelSetArcLength();
  void setUpBodyToSetTable();
  void computeBodyPropertiesForced(MInt returnMode, MFloat* bodyData, MInt body, MFloat time,
                                   MBool printPosition = false);
  void identifyBodies(MInt mode = 0);
  void setUpLevelSetInterpolationStencil(MInt cellId, MInt* interpolationCells, MInt position);
  void shiftOldLevelSetField(MInt dir, MInt set, MInt body);
  void buildCollectedLevelSet(MInt mode = 1);
  void reBuildCollectedLevelSet(MInt mode);
  void gapHandling();
  void levelSetGapCorrect();
  void levelSetGapRecorrect();
  void finalizeLevelSetInitialization();
  void setUpPotentialGapCells();
  void allocateLevelSetMemory();
  void allocateRotatingLs();

  void initializeCollectedLevelSet(MInt mode);
  void regionGrowing(MInt cellId, MInt region);
  void setChildRegions(MInt cellId, MInt region);

  MBool levelSetSolver();
  MBool gRungeKutta();
  MBool semiLagrangeTimeStep();
  MBool regridLevelSet();
  MBool levelSetReinitializationTrigger();
  MBool inCell(MInt cellId, MFloat* point);
  MBool gapCellsExist();
  MBool localGapCellsExist();
  MFloat firstOrderEikonalSolver(MInt cellListSize, MInt maxIterations, MInt set);
  MFloat secondOrderEikonalSolver(MFloat* q, const MInt* nghbrs, MInt cellListSize, MInt maxIterations,
                                  MInt set); // used by sohel
  MFloat fifthOrderEikonalSolver(MInt cellListSize, MInt maxIterations, MInt* crCells, MInt noCRCells, MFloat* factors,
                                 MInt crMode, MInt set);
  MFloat computeDistanceFromSTL(MFloat* target, MInt* closestElement, MFloat* closestPoint, MInt set,
                                MFloat sphereRadiusFactor = F5);
  MFloat interpolateOldLevelSet(MInt* interpolationCells, MFloat* point, MInt referenceSet);
  MFloat interpolateLevelSet(MInt* interpolationCells, MFloat* point, MInt referenceSet);
  MInt determineLevelSetSignFromSTL(MFloat* target, MInt set);
  MInt hyperbolicExtensionOpt(MFloat* q, MInt* cellList, MInt cellListSize, MFloat convergenceCriterion, MInt set);
  MInt getContainingCell(MFloat* point);
  MInt getContainingCell(MInt startCell, MFloat* point, MInt set = -1);
  MInt setUpLevelSetInterpolationStencil(MInt cellId, MInt* interpolationCells, MFloat* point);

  void exchangeIntBuffers(MInt*, MInt*, MInt, MInt);
  template <typename T>
  void exchangeBuffersGlobal(T* sendBuffer, T* receiveBuffer, MInt*, MInt*, MInt*, MInt*, MInt, MInt offset = 1);

  void readLevelSetProperties();
  void preTimeStep() override;
  void postTimeStep() override;
  void buildMultipleLevelSet(MInt mode = 1);
  MBool _levelSetSolutionStep();
  MBool finalizeLevelSet_(const MInt t_levelSet, const MInt t_output);
  MBool solutionStep() override;
  void finalizeLevelSet();

  void saveRestartFile(const MBool, MInt*);
  void writeRestartFile(const MBool, const MBool, const MString, MInt* recalcIdTree) override;
  void reIntAfterRestart(MBool) override;
  MBool prepareRestart(MBool, MBool&) override;
  MFloat reduceData(const MInt, MFloat* data, const MInt dataBlockSize = 1);

  MFloat crankAngle(const MFloat, const MInt);

  void determinePeriodicDistance();

 private:
  void reInitSolver(const MBool);
  void initializeTimers();

  template <MBool currentLevelSet>
  void exchangeLeafDataLS();
  template <typename T>
  void exchangeDataLS(T* data, const MInt dataSize = 1);
};

#undef ENSURE_VALID_GCELLID
#undef ENSURE_VALID_DIM
#undef ENSURE_VALID_DIR
#undef ENSURE_VALID_SET

#endif // ifndef LSSOLVER_H_
