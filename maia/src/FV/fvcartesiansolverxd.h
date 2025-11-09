// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef FVSOLVERXD_H
#define FVSOLVERXD_H

#include <algorithm>
#include <memory>
#include "COMM/mpiexchange.h"
#include "COMM/mpioverride.h"
#include "GRID/cartesiangrid.h"
#include "IO/iovtk.h"
#include "IO/parallelio.h"
#include "POST/samplingdata.h"
#include "UTIL/maiafftw.h"
#include "UTIL/tensor.h"
#include "cartesiansolver.h"
#include "fvcartesianbndrycell.h"
#include "fvcartesianbndrycndxd.h"
#include "fvcartesiancellcollector.h"
#include "fvcartesiansurfacecollector.h"
#include "fvcartesiansyseqndetchem.h"
#include "fvcartesiansyseqneegas.h"
#include "fvcartesiansyseqnns.h"
#include "fvcartesiansyseqnrans.h"
#include "fvcartesiansyseqntraits.h"
#include "fvcartesianwmsurface.h"

#include "UTIL/oneDFlame.h"

template <MInt nDim>
class FvSurface;

template <MInt nDim>
class FvSurfaceCollector;

template <MInt nDim, class SysEqn>
class CouplerFvMultilevel;

template <MInt nDim, class SysEqn>
class FvBndryCndXD;

template <class T>
class List;

template <MInt nDim, class SysEqn>
class LsFvMb;

template <MInt nDim, class SysEqn>
class LsFvCombustion;

template <MInt nDim, class SysEqn>
class CouplingLsFv;

template <MInt nDim, class SysEqn>
class FvZonal;

template <MInt nDim, class SysEqn>
class FvZonalRTV;

template <MInt nDim, class SysEqn>
class FvZonalSTG;

template <MInt nDim, class SysEqnOld, class SysEqnNew>
class FvCartesianInterpolation;

template <MInt nDim, class SysEqn>
class CouplingFvMb;

template <MInt nDim, class SysEqn>
class CouplingFv;

template <MInt nDim, MInt nDist, class SysEqnLb, class SysEqnFv>
class CouplerLbFvEEMultiphase;

template <MInt nDim, class SysEqn>
class CouplerFvParticle;

template <MInt nDim, class ppType>
class PostProcessing;

template <MInt nDim, class SysEqn>
class PostProcessingFv;

template <MInt nDim, class SysEqn>
class VtkIo;

//---------------------------------------------------------------------------

// Sanity-checking macros
#if !defined(NDEBUG) || defined(MAIA_ASSERT_ACCESSORS)
#define ENSURE_VALID_CELL_ID_CONTAINER(id)                                                                             \
  do {                                                                                                                 \
    if(id < 0 || id >= maxNoGridCells()) {                                                                             \
      TERMM(1, "Invalid cell id " + std::to_string((id)) + "\n\n AT: " + AT_);                                         \
    }                                                                                                                  \
  } while(false)
#else
#define ENSURE_VALID_CELL_ID_CONTAINER(id)                                                                             \
  do {                                                                                                                 \
  } while(false)
#endif


namespace maia {
namespace fv {

// Create struct for easy timer identification
struct Timers_ {
  // Enum to store timer "names"
  enum {
    SolverType,
    PreTime,
    TimeInt,

    ReinitSolu,
    BndryMb,
    NearBndry,
    InitSmallCorr,
    GhostCells,

    Rhs,
    Muscl,
    Ausm,
    ViscFlux,
    DistFlux,
    RhsMisc,

    SurfaceCoefficients,
    SurfaceMeanMolarWeight,
    SurfaceTransportCoefficients,
    CellCenterCoefficients,
    SpeciesReactionRates,

    SandpaperTrip,
    WMExchange,
    WMSurfaceLoop,
    WMFluxCorrection,

    MusclReconst,
    MusclCopy,
    MusclGhostSlopes,
    MusclCutSlopes,
    MusclReconstSrfc,

    ReconstSrfcCompValues,
    ReconstSrfcCorrVars,
    ReconstSrfcUpdateGhost,
    ReconstSrfcUpdateCutOff,

    RhsBnd,
    Lhs,

    LhsBnd,
    SmallCellCorr,
    ComputePV,
    Exchange,
    CutOff,
    BndryCnd,

    SCCorrInit,
    SCCorrExchange1,
    SCCorrExchange1Wait,
    SCCorrInterp,
    SCCorrRedist,
    SCCorrExchange2,
    SCCorrExchange2Wait,

    Residual,
    ResidualMpi,

    PostTime,
    PostSolu,
    ResidualMb,
    SurfaceForces,
    LevelSetCorr,

    /* Accumulated, */
    /* IO, */
    /* MPI, */

    // Special enum value used to initialize timer array
    _count
  };
};

} // namespace fv
} // namespace maia

template <class T>
class Collector;

template <MInt nDim_, class SysEqn>
class FvCartesianSolverXD : public maia::CartesianSolver<nDim_, FvCartesianSolverXD<nDim_, SysEqn>> {
 public:
  static constexpr MInt nDim = nDim_;
  static constexpr const MInt m_noDirs = 2 * nDim;
  static constexpr const MInt m_noEdges = nDim == 3 ? 12 : 4;
  static constexpr const MInt m_noCorners = nDim == 3 ? 8 : 4;
  static constexpr MFloat m_volumeThreshold = nDim == 3 ? 1e-16 : 1e-14; // why two different values?

  Geometry<nDim>* m_geometry;
  MFloat m_eps;
  MInt m_noSpecies;

  SysEqn m_sysEqn;

 private:
  template <class SolverType>
  friend class AccessorUnstructured;

  template <MInt nDim, SolverType SolverTypeR, SolverType SolverTypeL>
  friend class MSTG;
  template <MInt nDim, class SysEqn_>
  friend class FvBndryCndXD;
  template <MInt nDim, class SysEqn_>
  friend class CouplerFvMultilevel;
  template <MInt nDim, class SysEqn_>
  friend class LsFvMb;
  template <MInt nDim, class SysEqn_>
  friend class CouplingLsFv;
  template <MInt nDim, class SysEqn_>
  friend class CouplingFvMb;
  template <MInt nDim, class SysEqn_>
  friend class LsFvCombustion;
  template <MInt nDim, MInt nDist, class SysEqnLb_, class SysEqnFv_>
  friend class CouplerLbFv;

  template <MInt nDim, MInt nDist, class SysEqnLb_, class SysEqnFv_>
  friend class CouplerLbFvEEMultiphase;
  template <MInt nDim, class SysEqn_>
  friend class CouplerFvParticle;
  template <MInt nDim, class ppType>
  friend class PostProcessing;
  template <MInt nDim, class SysEqn_>
  friend class PostProcessingFv;
  friend class maia::CartesianSolver<nDim, FvCartesianSolverXD<nDim, SysEqn>>;
  friend class LsFvCombustion<nDim, SysEqn>;
  friend class CouplingLsFv<nDim, SysEqn>;
  using Self = FvCartesianSolverXD<nDim, SysEqn>;

 protected:
  using Timers = maia::fv::Timers_;

 public:
  // Type for cell properties
  using Cell = typename maia::grid::tree::Tree<nDim>::Cell;
  using FvSurfaceCollector = maia::fv::surface_collector::FvSurfaceCollector<nDim>;
  using SolverCell = FvCell;
  using CartesianSolver = typename maia::CartesianSolver<nDim, FvCartesianSolverXD>;
  using Grid = typename CartesianSolver::Grid;
  using GridProxy = typename CartesianSolver::GridProxy;

  // used CartesianSolver
  using CartesianSolver::c_noCells;
  using CartesianSolver::disableDlbTimers;
  using CartesianSolver::domainId;
  using CartesianSolver::domainOffset;
  using CartesianSolver::exchangeData;
  using CartesianSolver::extractPointIdsFromGrid;
  using CartesianSolver::getIdentifier;
  using CartesianSolver::grid;
  using CartesianSolver::haloCellId;
  using CartesianSolver::isActive;
  using CartesianSolver::m_bandWidth;
  using CartesianSolver::m_freeIndices;
  using CartesianSolver::m_innerBandWidth;
  using CartesianSolver::m_Ma;
  using CartesianSolver::m_outerBandWidth;
  using CartesianSolver::m_Re;
  using CartesianSolver::m_recalcIds;
  using CartesianSolver::m_residualInterval;
  using CartesianSolver::m_restart;
  using CartesianSolver::m_restartFile;
  using CartesianSolver::m_restartInterval;
  using CartesianSolver::m_restartOffset;
  using CartesianSolver::m_restartTimeStep;
  using CartesianSolver::m_revDir;
  using CartesianSolver::m_sensorParticle;
  using CartesianSolver::m_solutionInterval;
  using CartesianSolver::m_solutionOffset;
  using CartesianSolver::m_solutionOutput;
  using CartesianSolver::m_solutionTimeSteps;
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

  // pointer to the boundary condition
  FvBndryCndXD<nDim_, SysEqn>* m_fvBndryCnd;


  FvCartesianSolverXD() = delete;
  FvCartesianSolverXD(MInt, MInt, const MBool*, maia::grid::Proxy<nDim_>& gridProxy_, Geometry<nDim_>& geometry_,
                      const MPI_Comm comm);

  ~FvCartesianSolverXD() {
    delete[] m_variablesName;
    delete[] m_vorticityName;

    delete m_fvBndryCnd;

    releaseMemory();

    RECORD_TIMER_STOP(m_timers[Timers::SolverType]);
  }

  SysEqn sysEqn() const { return m_sysEqn; };
  SysEqn& sysEqn() { return m_sysEqn; };

 public:
  /// \brief Returns the number of cells
  MInt a_noCells() const {
#ifndef NDEBUG
    MInt noCells_solver = m_cells.size();
    MInt noCells_grid = grid().tree().size() + m_totalnosplitchilds + m_totalnoghostcells;
    ASSERT(noCells_solver == noCells_grid,
           "Warning, noCells differ: grid " << noCells_grid << ", solver " << noCells_solver << "no-split-childs "
                                            << m_totalnosplitchilds
                                            << "no-boundary-ghost-cells: " << m_totalnoghostcells << std::endl);
    //    ASSERT(c_noCells() == grid().m_noRegularCells, "ERROR in Cell-count");
#endif
    return m_cells.size();
  }

  MInt noInternalCells() const override { return grid().noInternalCells(); }

  /// \brief Returns isGapCell of the cell \p cellId
  MBool a_isGapCell(const MInt cellId) const { return m_cells.hasProperty(cellId, SolverCell::IsGapCell); }
  /// \brief Returns isGapCell of the cell \p cellId
  maia::fv::cell::BitsetType::reference a_isGapCell(const MInt cellId) {
    return m_cells.hasProperty(cellId, SolverCell::IsGapCell);
  }
  /// \brief Returns wasGapCell of the cell \p cellId
  maia::fv::cell::BitsetType::reference a_wasGapCell(const MInt cellId) {
    return m_cells.hasProperty(cellId, SolverCell::WasGapCell);
  }
  /// \brief Returns IsHalo of the cell \p cellId
  MBool a_isHalo(const MInt cellId) const { return m_cells.hasProperty(cellId, SolverCell::IsHalo); }
  /// \brief Returns IsHalo of the cell \p cellId
  maia::fv::cell::BitsetType::reference a_isHalo(const MInt cellId) {
    return m_cells.hasProperty(cellId, SolverCell::IsHalo);
  }
  /// \brief Returns IsWindow of the cell \p cellId
  MBool a_isWindow(const MInt cellId) const { return m_cells.hasProperty(cellId, SolverCell::IsWindow); }
  /// \brief Returns IsWindow of the cell \p cellId
  maia::fv::cell::BitsetType::reference a_isWindow(const MInt cellId) {
    return m_cells.hasProperty(cellId, SolverCell::IsWindow);
  }
  /// \brief Returns IsPeriodic of the cell \p cellId
  MBool a_isPeriodic(const MInt cellId) const { return m_cells.hasProperty(cellId, SolverCell::IsPeriodic); }
  /// \brief Returns IsPeriodic of the cell \p cellId
  maia::fv::cell::BitsetType::reference a_isPeriodic(const MInt cellId) {
    return m_cells.hasProperty(cellId, SolverCell::IsPeriodic);
  }

  /// \brief Returns isBndryCell of the cell \p cellId
  MBool a_isBndryCell(const MInt cellId) const override { return (a_bndryId(cellId) >= 0); }

  /// \brief Returns isInterface of the cell \p cellId
  MBool a_isInterface(const MInt cellId) const { return m_cells.hasProperty(cellId, SolverCell::IsInterface); }

  /// \brief Returns isInterface of the cell \p cellId
  maia::fv::cell::BitsetType::reference a_isInterface(const MInt cellId) {
    return m_cells.hasProperty(cellId, SolverCell::IsInterface);
  }

  /// \brief Returns isBndryGhostCell of the cell \p cellId
  MBool a_isBndryGhostCell(const MInt cellId) const { return m_cells.hasProperty(cellId, SolverCell::IsGhost); }

  /// \brief Returns isBndryGhostCell of the cell \p cellId
  maia::fv::cell::BitsetType::reference a_isBndryGhostCell(const MInt cellId) {
    return m_cells.hasProperty(cellId, SolverCell::IsGhost);
  }

  /// \brief Returns isWMImgCell of the cell \p cellId
  MBool a_isWMImgCell(const MInt cellId) const { return m_cells.hasProperty(cellId, SolverCell::IsWMImgCell); }
  /// \brief Returns isWMImgCell of the cell \p cellId
  maia::fv::cell::BitsetType::reference a_isWMImgCell(const MInt cellId) {
    return m_cells.hasProperty(cellId, SolverCell::IsWMImgCell);
  }

  /// \brief Returns isWMImgCell of the cell \p cellId
  MBool a_isSandpaperTripCell(const MInt cellId) const {
    return m_cells.hasProperty(cellId, SolverCell::IsSandpaperTripCell);
  }
  /// \brief Returns isWMImgCell of the cell \p cellId
  maia::fv::cell::BitsetType::reference a_isSandpaperTripCell(const MInt cellId) {
    return m_cells.hasProperty(cellId, SolverCell::IsSandpaperTripCell);
  }

  /// \brief Returns properties of the cell \p cellId
  maia::fv::cell::BitsetType& a_properties(const MInt cellId) { return m_cells.properties(cellId); }

  /// \brief Returns the spongeFactor of the cell \p cellId
  MFloat& a_spongeFactor(const MInt cellId) { return m_cells.spongeFactor(cellId); }

  /// \brief Returns the spongeFactor of the cell \p cellId
  MFloat a_spongeFactor(const MInt cellId) const { return m_cells.spongeFactor(cellId); }

  /// \brief Returns the coordinate of the cell from the fvcellcollector \p cellId for dimension \p dir
  MFloat& a_coordinate(const MInt cellId, const MInt dir) { return m_cells.coordinate(cellId, dir); }

  /// \brief Returns the coordinate of the cell from the fvcellcollector \p cellId for dimension \p dir
  MFloat a_coordinate(const MInt cellId, const MInt dir) const { return m_cells.coordinate(cellId, dir); }

  /// \brief Returns the level of the cell from the fvcellcollector \p cellId
  MInt& a_level(const MInt cellId) { return m_cells.level(cellId); }

  /// \brief Returns the level of the cell from the fvcellcollector \p cellId
  MInt a_level(const MInt cellId) const { return m_cells.level(cellId); }

  /// \brief Returns the cell volume of the cell from the fvcellcollector \p cellId
  MFloat& a_cellVolume(const MInt cellId) { return m_cells.cellVolume(cellId); }

  /// \brief Returns the cell volume of the cell from the fvcellcollector \p cellId
  MFloat a_cellVolume(const MInt cellId) const { return m_cells.cellVolume(cellId); }

  /// \brief Returns the inverse cell volume of the cell from the fvcellcollector \p cellId
  MFloat& a_FcellVolume(const MInt cellId) override { return m_cells.FcellVolume(cellId); }

  /// \brief Returns the inverse cell volume of the cell from the fvcellcollector \p cellId
  MFloat a_FcellVolume(const MInt cellId) const { return m_cells.FcellVolume(cellId); }

  /// \brief Returns the level for cutCells, this can either be the maxRefinementLevel
  ///        or the level of the current \p cellId
  MInt a_cutCellLevel(const MInt cellId) const {
    if(!m_bndryLevelJumps) {
      return maxRefinementLevel();
    } else {
      return a_level(cellId);
    }
  }

  MBool a_isInactive(const MInt cellId) const { return a_hasProperty(cellId, SolverCell::IsInactive); }
  MBool a_isActive(const MInt cellId) const { return a_hasProperty(cellId, SolverCell::IsActive); }
  MBool a_isSplitCell(const MInt cellId) const { return a_hasProperty(cellId, SolverCell::IsSplitCell); }

  MInt a_noSplitCells() const { return m_splitCells.size(); }
  MInt a_noSplitChilds(const MInt sc) const { return m_splitChilds[sc].size(); }
  MInt a_splitChildId(const MInt sc, const MInt ssc) { return m_splitChilds[sc][ssc]; }
  MInt a_splitCellId(const MInt sc) const { return m_splitCells[sc]; }

  //------------Grid-properties:

  /// \brief Cecks wether the cell \p cellId is an actual grid-cell
  // (splitchilds and boundary-ghostcells are no longer added to the grid-tree!)
#ifndef NDEBUG
  void assertValidGridCellId(const MInt cellId) const {
    if(a_isBndryGhostCell(cellId)
       || cellId < 0
       //       || m_cells.hasProperty(cellId, SolverCell::IsSplitCell)
       || m_cells.hasProperty(cellId, SolverCell::IsSplitClone)
       || m_cells.hasProperty(cellId, SolverCell::IsSplitChild)) {
      TERMM(-1, "ERROR: Invalid cell " + std::to_string(cellId) + "/" + std::to_string(grid().tree().size())
                    + " isSplitClone:"
                    + std::to_string(m_cells.hasProperty(cellId, SolverCell::IsSplitClone))
                    //                    + " isSplitCell:" + std::to_string(m_cells.hasProperty(cellId,
                    //                    SolverCell::IsSplitCell))
                    + " IsSplitChild:" + std::to_string(m_cells.hasProperty(cellId, SolverCell::IsSplitChild))
                    + " isBndryGhost:" + std::to_string(a_isBndryGhostCell(cellId)));
    }
  }
#else
  void assertValidGridCellId(const MInt NotUsed(cellId)) const {}
#endif


  /// \brief Checks wether the cell \p cellId has a valid neighbor in direction \p dir
  // an invalid neighbor can occour if the m_deleteNeighbour-flag is activated and if the
  // cell, and neighbor has been added to  m_nbBackupCellIds-list in the deleteNeighbourLinks().
  // NOTE: If this error-occours, the checkNeighborActive() needs to be added to that code-line!
  void assertDeleteNeighbor(const MInt cellId, const MInt dir) const {
#ifndef NDEBUG
    if(m_deleteNeighbour && a_hasProperty(cellId, SolverCell::IsBndryActive) && grid().tree().hasNeighbor(cellId, dir)
       && a_hasProperty(grid().tree().neighbor(cellId, dir), SolverCell::IsInactive)) {
      ASSERT(false, "ERROR: calling inactive-neighbor!");
    }
#else
    std::ignore = cellId;
    std::ignore = dir;
#endif
  }

  /// \brief Cecks wether the cell \p cellId has a valid neighbor in direction \p dir
  MBool checkNeighborActive(const MInt cellId, const MInt dir) const {
    if(m_deleteNeighbour && a_hasProperty(cellId, SolverCell::IsBndryActive) && grid().tree().hasNeighbor(cellId, dir)
       && a_hasProperty(grid().tree().neighbor(cellId, dir), SolverCell::IsInactive)) {
      return false;
    }

    assertDeleteNeighbor(cellId, dir);
    return true;
  }

  /// \brief Returns the coordinate of the cell from the grid().tree() \p cellId for dimension \p dir
  MFloat c_coordinate(const MInt cellId, const MInt dir) const {
    assertValidGridCellId(cellId);
    return grid().tree().coordinate(cellId, dir);
  }

  MBool c_isLeafCell(const MInt cellId) const { return grid().tree().isLeafCell(cellId); }

  /// \brief Returns the grid level of the cell \p cellId
  MInt c_level(const MInt cellId) const {
    assertValidGridCellId(cellId);
    /*    if ( a_isBndryGhostCell(cellId) || a_hasProperty( cellId, SolverCell::IsSplitChild ) ) {
          const MInt aCellId = getAssociatedInternalCell(cellId);
          assertValidGridCellId(aCellId);
          return grid().tree().level(aCellId);
        }
    */
    return grid().tree().level(cellId);
  }

  /// \brief Returns the number of children of the cell \p cellId
  MInt c_noChildren(const MInt cellId) const {
    assertValidGridCellId(cellId);
    return grid().tree().noChildren(cellId);
  }

  /// \brief Returns the grid parent id of the cell \p cellId
  MLong c_parentId(const MInt cellId) const {
    assertValidGridCellId(cellId);
    return grid().tree().parent(cellId);
  }

  /// \brief Returns the global grid id of the grid cell \p cellId
  MLong c_globalId(const MInt cellId) const {
    assertValidGridCellId(cellId);
    return grid().tree().globalId(cellId);
  }

  /// \brief Returns the weight of the cell \p cellId
  MFloat c_weight(const MInt cellId) const {
    assertValidGridCellId(cellId);
    return grid().tree().weight(cellId);
  }

  /// \brief Returns the grid child id of the grid cell \p cellId at position \p pos
  MLong c_childId(const MInt cellId, const MInt pos) const {
    assertValidGridCellId(cellId);
    return grid().tree().child(cellId, pos);
  }

  /// \brief Returns the grid neighbor id of the grid cell \p cellId \p dir
  MLong c_neighborId(const MInt cellId, const MInt dir, const MBool assertNeighborState = true) const {
    assertValidGridCellId(cellId);
    if(assertNeighborState) assertDeleteNeighbor(cellId, dir);
    return grid().tree().neighbor(cellId, dir);
  }

  /// \brief Returns the length of the cell for \p level
  MFloat c_cellLengthAtCell(const MInt cellId) const {
    ASSERT(a_level(cellId) > 0, "Level-error");
    return grid().cellLengthAtLevel(a_level(cellId));
  }
  /// \brief Returns the length of the cell for \p level
  MFloat c_cellLengthAtLevel(const MInt level) const { return grid().cellLengthAtLevel(level); }

  /// \brief Returns the grid Volume of the cell for \p level
  MFloat c_cellVolumeAtLevel(const MInt level) const { return grid().cellVolumeAtLevel(level); }

  /// \brief Returns grid cell property \p p of the cell \p cellId
  MBool a_hasProperty(const MInt cellId, const Cell p) const {
    assertValidGridCellId(cellId);
    return grid().tree().hasProperty(cellId, p);
  }


  /// \brief Returns the delete of the cell \p cellId
  MBool c_isToDelete(const MInt cellId) const {
    assertValidGridCellId(cellId);
    return grid().tree().hasProperty(cellId, Cell::IsToDelete);
  }

  /// \brief Returns noNeighborIds of the cell \p CellId for direction \p dir
  MInt a_hasNeighbor(const MInt cellId, const MInt dir, const MBool assertNeighborState = true) const {
    assertValidGridCellId(cellId);
    if(assertNeighborState) assertDeleteNeighbor(cellId, dir);
    return grid().tree().hasNeighbor(cellId, dir);
  }

  //-------------End-grid-properties

  /// \brief Returns conservative variable \p v of the cell \p cellId for variables \p varId
  MFloat& a_variable(const MInt cellId, const MInt varId) { return m_cells.variable(cellId, varId); }

  /// \brief Returns conservative variable \p v of the cell \p cellId for variables \p varId
  MFloat a_variable(const MInt cellId, const MInt varId) const { return m_cells.variable(cellId, varId); }

  /// \brief Returns primitive variable \p v of the cell \p cellId for variables \p varId
  MFloat& a_pvariable(const MInt cellId, const MInt varId) { return m_cells.pvariable(cellId, varId); }

  /// \brief Returns primitive variable \p v of the cell \p cellId for variables \p varId
  MFloat a_pvariable(const MInt cellId, const MInt varId) const { return m_cells.pvariable(cellId, varId); }

  /// \brief Returns additional variable \p v of the cell \p cellId for variables \p varId
  MFloat& a_avariable(const MInt cellId, const MInt varId) { return m_cells.avariable(cellId, varId); }

  /// \brief Returns additional variable \p v of the cell \p cellId for variables \p varId
  MFloat a_avariable(const MInt cellId, const MInt varId) const { return m_cells.avariable(cellId, varId); }

  /// \brief Returns solver cell property \p p of the cell \p cellId
  maia::fv::cell::BitsetType::reference a_hasProperty(const MInt cellId, const SolverCell p) {
    return m_cells.hasProperty(cellId, p);
  }

  /// \brief Returns solver cell property \p p of the cell \p cellId
  MBool a_hasProperty(const MInt cellId, const SolverCell p) const { return m_cells.hasProperty(cellId, p); }

  /// \brief Returns property \p p of the cell \p cellId
  void a_resetPropertiesSolver(const MInt cellId) { m_cells.resetProperties(cellId); }

  /// \brief Returns property \p p of the cell \p cellId
  void a_copyPropertiesSolver(const MInt fromCellId, const MInt toCellId) {
    m_cells.properties(toCellId) = m_cells.properties(fromCellId);
  }

  /// \brief Returns reconstruction neighbor \p n of the cell \p cellId
  MInt& a_reconstructionNeighborId(const MInt cellId, const MInt nghbrNo) {
    return m_cells.rcnstrctnNghbrId(cellId, nghbrNo);
  }

  /// \brief Returns reconstruction neighbor \p n of the cell \p cellId
  MInt a_reconstructionNeighborId(const MInt cellId, const MInt nghbrNo) const {
    return m_cells.rcnstrctnNghbrId(cellId, nghbrNo);
  }

  /// \brief Returns reconstruction data offset \p i of the cell \p cellId
  MInt& a_reconstructionData(const MInt cellId) { return m_cells.reconstructionData(cellId); }

  /// \brief Returns reconstruction data offset \p i of the cell \p cellId
  MInt a_reconstructionNeighborId(const MInt cellId) const { return m_cells.reconstructionData(cellId); }

  /// \brief Returns the spongeBndryId of the cell \p cellId for direction \p dir
  MInt& a_spongeBndryId(const MInt cellId, const MInt dir) { return m_cells.spongeBndryId(cellId, dir); }

  /// \brief Returns the spongeBndryId of the cell \p cellId for direction \p dir
  MInt a_spongeBndryId(const MInt cellId, const MInt dir) const { return m_cells.spongeBndryId(cellId, dir); }

  /// \brief Returns the spongeFactorStart of the cell \p cellId
  MFloat& a_spongeFactorStart(const MInt cellId) { return m_cells.spongeFactorStart(cellId); }

  /// \brief Returns the spongeFactorStart of the cell \p cellId
  MFloat a_spongeFactorStart(const MInt cellId) const { return m_cells.spongeFactorStart(cellId); }

  /// \brief Returns the bndryId of the cell \p cellId
  MInt& a_bndryId(const MInt cellId) { return m_cells.bndryCellId(cellId); }

  /// \brief Returns the bndryId of the cell \p cellId
  MInt a_bndryId(const MInt cellId) const { return m_cells.bndryCellId(cellId); }

  /// \brief Returns the noRcnstrctnNghbrIds of the cell \p cellId
  MInt& a_noReconstructionNeighbors(const MInt cellId) { return m_cells.noRcnstrctnNghbrIds(cellId); }

  /// \brief Returns the noRcnstrctnNghbrIds of the cell \p cellId
  MInt a_noReconstructionNeighbors(const MInt cellId) const { return m_cells.noRcnstrctnNghbrIds(cellId); }

  /// \brief Returns the tau of the cell \p cellId for variable \p varId
  MFloat& a_tau(const MInt cellId, const MInt varId) { return m_cells.tau(cellId, varId); }

  /// \brief Returns the tau of the cell \p cellId for variable \p varId
  MFloat a_tau(const MInt cellId, const MInt varId) const { return m_cells.tau(cellId, varId); }

  /// \brief Returns the restrictedRHS of the cell \p cellId for variable \p varId
  MFloat& a_restrictedRHS(const MInt cellId, const MInt varId) { return m_cells.restrictedRHS(cellId, varId); }

  /// \brief Returns the restrictedRHS of the cell \p cellId for variable \p varId
  MFloat a_restrictedRHS(const MInt cellId, const MInt varId) const { return m_cells.restrictedRHS(cellId, varId); }

  /// \brief Returns restricted variables of cell \p cellId for variable \p varId on level \p level
  MFloat& a_restrictedVar(const MInt cellId, const MInt varId) { return m_cells.restrictedVar(cellId, varId); }

  /// \brief Returns restricted variables of cell \p cellId for variable \p varId on level \p level
  MFloat a_restrictedVar(const MInt cellId, const MInt varId) const { return m_cells.restrictedVar(cellId, varId); }

  /// \brief Returns the reactionRate of the cell \p cellId for variables \p varId
  MFloat& a_reactionRate(const MInt cellId, const MInt reactionId) { return m_cells.reactionRate(cellId, reactionId); }

  /// \brief Returns the reactionRate of the cell \p cellId for variables \p varId
  MFloat a_reactionRate(const MInt cellId, const MInt reactionId) const {
    return m_cells.reactionRate(cellId, reactionId);
  }

  /// \brief Returns the reactionRateBackup of the cell \p cellId for variables \p varId
  MFloat& a_reactionRateBackup(const MInt cellId, const MInt reactionId) {
    return m_cells.reactionRateBackup(cellId, reactionId);
  }

  /// \brief Returns the reactionRateBackup of the cell \p cellId for variables \p varId
  MFloat a_reactionRateBackup(const MInt cellId, const MInt reactionId) const {
    return m_cells.reactionRateBackup(cellId, reactionId);
  }

  /// \brief Returns psi of the cell \p cellId for variables \p varId
  MFloat& a_psi(const MInt cellId) { return m_cells.psi(cellId); }

  /// \brief Returns psi of the cell \p cellId for variables \p varId
  MFloat a_psi(const MInt cellId) const { return m_cells.psi(cellId); }

  // detailed chemistry accessors
  MFloat& a_speciesReactionRate(const MInt cellId, const MInt speciesIndex) {
    return m_cells.speciesReactionRate(cellId, speciesIndex);
  }

  MFloat a_speciesReactionRate(const MInt cellId, const MInt speciesIndex) const {
    return m_cells.speciesReactionRate(cellId, speciesIndex);
  }

  /// \brief Returns the implicitCoefficient of cell \p cellId for coefficient \p coefId
  MFloat& a_implicitCoefficient(const MInt cellId, const MInt coefId) {
    return m_cells.implicitCoefficient(cellId, coefId);
  }

  /// \brief Returns the implicitCoefficient of cell \p cellId for coefficient \p coefId
  MFloat a_implicitCoefficient(const MInt cellId, const MInt coefId) const {
    return m_cells.implicitCoefficient(cellId, coefId);
  }

  //----------------------------------

  /// \brief Returns dt1Variables of the cell \p CellId variables \p varId
  MFloat& a_dt1Variable(const MInt cellId, const MInt varId) { return m_cells.dt1Variable(cellId, varId); }

  /// \brief Returns dt1Variables of the cell \p CellId variables \p varId
  MFloat a_dt1Variable(const MInt cellId, const MInt varId) const { return m_cells.dt1Variable(cellId, varId); }

  /// \brief Returns dt2Variables of the cell \p CellId variables \p varId
  MFloat& a_dt2Variable(const MInt cellId, const MInt varId) { return m_cells.dt2Variable(cellId, varId); }

  /// \brief Returns dt2Variables of the cell \p CellId variables \p varId
  MFloat a_dt2Variable(const MInt cellId, const MInt varId) const { return m_cells.dt2Variable(cellId, varId); }

  /// \brief Returns oldVariables\p v of the cell \p cellId variables \p varId
  MFloat& a_oldVariable(const MInt cellId, const MInt varId) { return m_cells.oldVariable(cellId, varId); }

  /// \brief Returns oldVariables \p v of the cell \p cellId variables \p varId
  MFloat a_oldVariable(const MInt cellId, const MInt varId) const { return m_cells.oldVariable(cellId, varId); }

  /// \brief Returns the slope of the cell \p cellId for the variable \p varId in direction \p dir
  MFloat& a_slope(const MInt cellId, MInt const varId, const MInt dir) override {
    return m_cells.slope(cellId, varId, dir);
  }

  /// \brief Returns the slope of the cell \p cellId for the variable \p varId in direction \p dir
  MFloat a_slope(const MInt cellId, MInt const varId, const MInt dir) const {
    return m_cells.slope(cellId, varId, dir);
  }

  /// \brief Returns the stored slope of the cell \p cellId for the variable \p varId in direction \p dir
  MFloat& a_storedSlope(const MInt cellId, MInt const varId, const MInt dir) {
    return m_cells.storedSlope(cellId, varId, dir);
  }

  /// \brief Returns the stored slope of the cell \p cellId for the variable \p varId in direction \p dir
  MFloat a_storedSlope(const MInt cellId, MInt const varId, const MInt dir) const {
    return m_cells.storedSlope(cellId, varId, dir);
  }

  /// \brief Returns the right hand side of the cell \p cellId for the variable \p varId
  MFloat& a_rightHandSide(const MInt cellId, MInt const varId) { return m_cells.rightHandSide(cellId, varId); }

  /// \brief Returns the right hand side of the cell \p cellId for the variable \p varId
  MFloat a_rightHandSide(const MInt cellId, MInt const varId) const { return m_cells.rightHandSide(cellId, varId); }

  //----------------------------------

  /// \brief Returns the number of surfaces
  MInt a_noSurfaces() { return m_surfaces.size(); }

  /// \brief Returns the boundary condition of surface \p srfcId
  MInt& a_surfaceBndryCndId(const MInt srfcId) { return m_surfaces.bndryCndId(srfcId); }
  /// \brief Returns the boundary condition of surface \p srfcId
  MInt a_surfaceBndryCndId(const MInt srfcId) const { return m_surfaces.bndryCndId(srfcId); }

  /// \brief Returns the orientation of surface \p srfcId
  MInt& a_surfaceOrientation(const MInt srfcId) { return m_surfaces.orientation(srfcId); }
  /// \brief Returns the orientation of surface \p srfcId
  MInt a_surfaceOrientation(const MInt srfcId) const { return m_surfaces.orientation(srfcId); }

  /// \brief Returns the area of surface \p srfcId
  MFloat& a_surfaceArea(const MInt srfcId) { return m_surfaces.area(srfcId); }
  /// \brief Returns the area of surface \p srfcId
  MFloat a_surfaceArea(const MInt srfcId) const { return m_surfaces.area(srfcId); }

  /// \brief Returns the factor of surface \p srfcId for variable \p varId
  MFloat& a_surfaceFactor(const MInt srfcId, const MInt varId) { return m_surfaces.factor(srfcId, varId); }
  /// \brief Returns the factor of surface \p srfcId for variable \p varId
  MFloat a_surfaceFactor(const MInt srfcId, const MInt varId) const { return m_surfaces.factor(srfcId, varId); }

  /// \brief Returns the coordinate of surface \p srfcId in direction \p dir
  MFloat& a_surfaceCoordinate(const MInt srfcId, const MInt dir) { return m_surfaces.coordinate(srfcId, dir); }
  /// \brief Returns the coordinate of surface \p srfcId in direction \p dir
  MFloat a_surfaceCoordinate(const MInt srfcId, const MInt dir) const { return m_surfaces.coordinate(srfcId, dir); }

  /// \brief Returns the delta X of surface \p srfcId for variable \p varId
  MFloat& a_surfaceDeltaX(const MInt srfcId, const MInt varId) { return m_surfaces.deltaX(srfcId, varId); }
  /// \brief Returns the delta X of surface \p srfcId for variable \p varId
  MFloat a_surfaceDeltaX(const MInt srfcId, const MInt varId) const { return m_surfaces.deltaX(srfcId, varId); }

  /// \brief Returns the neighbor cell id of surface \p srfcId in direction \p dir
  MInt& a_surfaceNghbrCellId(const MInt srfcId, const MInt dir) { return m_surfaces.nghbrCellId(srfcId, dir); }
  /// \brief Returns the neighbor cell id of surface \p srfcId in direction \p dir
  MInt a_surfaceNghbrCellId(const MInt srfcId, const MInt dir) const { return m_surfaces.nghbrCellId(srfcId, dir); }

  /// \brief Returns the variable \p varId of surface \p srfcId in direction \p dir
  MFloat& a_surfaceVariable(const MInt srfcId, const MInt dir, const MInt varId) {
    return m_surfaces.variable(srfcId, dir, varId);
  }
  /// \brief Returns the variable \p varId of surface \p srfcId in direction \p dir
  MFloat a_surfaceVariable(const MInt srfcId, const MInt dir, const MInt varId) const {
    return m_surfaces.variable(srfcId, dir, varId);
  }

  /// \brief Returns the upwind coefficient of surface \p srfcId
  MFloat& a_surfaceUpwindCoefficient(const MInt srfcId) { return m_surfaces.upwindCoefficient(srfcId); }
  /// \brief Returns the upwind coefficient of surface \p srfcId
  MFloat a_surfaceUpwindCoefficient(const MInt srfcId) const { return m_surfaces.upwindCoefficient(srfcId); }

  /// \brief Returns the coefficient \p dimCoefficient of surface \p srfcId
  MFloat& a_surfaceCoefficient(const MInt srfcId, const MInt dimCoefficient) {
    return m_surfaces.surfaceCoefficient(srfcId, dimCoefficient);
  }
  /// \brief Returns the coefficient \p dimCoefficient of surface \p srfcId
  MFloat a_surfaceCoefficient(const MInt srfcId, const MInt dimCoefficient) const {
    return m_surfaces.surfaceCoefficient(srfcId, dimCoefficient);
  }

  /// \brief Returns the flux \p fVarId for surface \p srfcId
  MFloat& a_surfaceFlux(const MInt srfcId, const MInt fVarId) { return m_surfaces.flux(srfcId, fVarId); }
  /// \brief Returns the flux \p fVarId for surface \p srfcId
  MFloat a_surfaceflux(const MInt srfcId, const MInt fVarId) const { return m_surfaces.flux(srfcId, fVarId); }

  //----------------------------------

  /// \brief Returns the Mach number of the solver
  MFloat& a_Ma() { return m_Ma; }

  /// \brief Returns the Mach number of the solver
  const MFloat& a_Ma() const { return m_Ma; }

  /// \brief Returns the cfl number of the solver
  MFloat& a_cfl() { return m_cfl; }

  // \brief Returns the cfl number of the solver
  const MFloat& a_cfl() const { return m_cfl; }

  /// \brief Returns the restart interval of the solver
  MInt& a_restartInterval() { return m_restartInterval; }

  /// \brief Returns the restart interval of the solver
  const MInt& a_restartInterval() const { return m_restartInterval; }

  /// \brief Return BndryCndId
  MInt& a_bndryCndId(MInt bndryId) { return m_bndryCells->a[bndryId].m_srfcs[0]->m_bndryCndId; }

  /// \brief Return BndryCndId
  const MInt& a_bndryCndId(MInt bndryId) const { return m_bndryCells->a[bndryId].m_srfcs[0]->m_bndryCndId; }

  /// \brief Return normal direction of bndry srfc
  MFloat& a_bndryNormal(MInt bndryId, MInt dir) { return m_bndryCells->a[bndryId].m_srfcs[0]->m_normalVector[dir]; }

  /// \brief Return normal direction of bndry srfc
  const MFloat& a_bndryNormal(MInt bndryId, MInt dir) const {
    return m_bndryCells->a[bndryId].m_srfcs[0]->m_normalVector[dir];
  }

  /// \brief Return cut coordinates of bndry srfc
  MFloat& a_bndryCutCoord(MInt bndryId, MInt i, MInt j) {
    return m_bndryCells->a[bndryId].m_srfcs[0]->m_cutCoordinates[i][j];
  }

  /// \brief Return cut coordinates of bndry srfc
  const MFloat& a_bndryCutCoord(MInt bndryId, MInt i, MInt j) const {
    return m_bndryCells->a[bndryId].m_srfcs[0]->m_cutCoordinates[i][j];
  }

  // \brief Return bndry ghost-cell id
  const MInt& a_bndryGhostCellId(const MInt bndryId, const MInt srfc) const {
    return m_bndryCells->a[bndryId].m_srfcVariables[srfc]->m_ghostCellId;
  }

  /// \brief Return ident nghbr Id
  MInt& a_identNghbrId(MInt nghbrId) { return m_identNghbrIds[nghbrId]; }

  /// \brief Return ident nghbr Id
  const MInt& a_identNghbrId(MInt nghbrId) const { return m_identNghbrIds[nghbrId]; }

  /// \brief Return store nghbr Id
  MInt& a_storeNghbrId(MInt nghbrId) { return m_storeNghbrIds[nghbrId]; }

  /// \brief Return store nghbr Id
  const MInt& a_storeNghbrId(MInt nghbrId) const { return m_storeNghbrIds[nghbrId]; }

  /// \brief Return mean flow velocity
  MFloat& a_VVInfinity(MInt dir) { return m_VVInfinity[dir]; }

  /// \brief Return mean flow velocity
  const MFloat& a_VVInfinity(MInt dir) const { return m_VVInfinity[dir]; }

  /// \brief Return rho infinity
  MFloat& a_rhoInfinity() { return m_rhoInfinity; }

  /// \brief Return rho infinity
  const MFloat& a_rhoInfinity() const { return m_rhoInfinity; }

  /// \brief Return p infinity
  MFloat& a_PInfinity() { return m_PInfinity; }

  /// \brief Return p infinity
  const MFloat& a_PInfinity() const { return m_PInfinity; }

  /// \brief Return T infinity
  MFloat& a_TInfinity() { return m_TInfinity; }

  /// \brief Return T infinity
  const MFloat& a_TInfinity() const { return m_TInfinity; }

  /// \brief Return time reference value
  MFloat& a_timeRef() { return m_timeRef; }

  /// \brief Return time reference value
  const MFloat& a_timeRef() const { return m_timeRef; }

  /// \brief Return no particles
  MInt& a_noPart(MInt cellId) { return m_noParts[cellId]; }

  /// \brief Return physical time
  MFloat& a_physicalTime() { return m_physicalTime; }

  /// \brief Return physical time
  const MFloat& a_physicalTime() const { return m_physicalTime; }

  /// \brief Return time
  MFloat& a_time() { return m_time; }

  /// \brief Return time
  const MFloat& a_time() const { return m_time; }

  /// \brief Return no particles
  const MInt& a_noPart(MInt cellId) const { return m_noParts[cellId]; }

  /// \brief Return external source
  MFloat& a_externalSource(MInt cellId, MInt var) { return m_externalSource[cellId][var]; }

  /// \brief Return external source
  const MFloat& a_externalSource(MInt cellId, MInt var) const { return m_externalSource[cellId][var]; }

  /// \brief Return max level window cells
  MInt& a_maxLevelWindowCells(MInt domain, MInt id) { return m_maxLevelWindowCells[domain][id]; }

  /// \brief Return max level window cells
  const MInt& a_maxLevelWindowCells(MInt domain, MInt id) const { return m_maxLevelWindowCells[domain][id]; }

  /// \brief Return max level halo cells
  MInt& a_maxLevelHaloCells(MInt domain, MInt id) { return m_maxLevelHaloCells[domain][id]; }

  /// \brief Return max level halo cells
  const MInt& a_maxLevelHaloCells(MInt domain, MInt id) const { return m_maxLevelHaloCells[domain][id]; }

  /// \brief Returns the local time-step of the cell \p cellId
  MFloat& a_localTimeStep(const MInt cellId) {
    return m_cells.localTimeStep(cellId);
  }

  /// \brief Returns the local time-step of the cell \p cellId
  MFloat a_localTimeStep(const MInt cellId) const {
    return m_cells.localTimeStep(cellId);
  }

  MFloat a_dynViscosity(const MFloat T) const { return SUTHERLANDLAW(T); }


  /// \brief Returns the levelSet-value for fv-CellId \p cellId and \p set
  MFloat& a_levelSetFunction(const MInt cellId, const MInt set) { return m_levelSetValues[set][cellId]; }

  /// \brief Returns the levelSet-value for fv-CellId \p cellId and \p set
  MFloat a_levelSetFunction(const MInt cellId, const MInt set) const { return m_levelSetValues[set][cellId]; }

  /// \brief Returns the levelSetMb-value for fv-CellId \p cellId and \p set
  MFloat& a_levelSetValuesMb(const MInt cellId, const MInt set) { return m_levelSetValuesMb[IDX_LSSETMB(cellId, set)]; }

  /// \brief Returns the levelSetMb-value for fv-CellId \p cellId and \p set
  MFloat a_levelSetValuesMb(const MInt cellId, const MInt set) const {
    return m_levelSetValuesMb[IDX_LSSETMB(cellId, set)];
  }

  MFloat a_alphaGas(const MInt cellId) const { return a_pvariable(cellId, PV->Y[0]); }

  MFloat& a_alphaGas(const MInt cellId) { return a_pvariable(cellId, PV->Y[0]); }

  MFloat a_uOtherPhase(const MInt cellId, const MInt dir) const { return m_EEGas.uOtherPhase[cellId][dir]; }

  MFloat& a_uOtherPhase(const MInt cellId, const MInt dir) { return m_EEGas.uOtherPhase[cellId][dir]; }

  MFloat a_uOtherPhaseOld(const MInt cellId, const MInt dir) const { return m_EEGas.uOtherPhaseOld[cellId][dir]; }

  MFloat& a_uOtherPhaseOld(const MInt cellId, const MInt dir) { return m_EEGas.uOtherPhaseOld[cellId][dir]; }

  MFloat a_gradUOtherPhase(const MInt cellId, const MInt uDir, const MInt gradDir) const {
    return m_EEGas.gradUOtherPhase[cellId][uDir + nDim * gradDir];
  }

  MFloat& a_gradUOtherPhase(const MInt cellId, const MInt uDir, const MInt gradDir) {
    return m_EEGas.gradUOtherPhase[cellId][uDir + nDim * gradDir];
  }

  MFloat a_vortOtherPhase(const MInt cellId, const MInt dir) const { return m_EEGas.vortOtherPhase[cellId][dir]; }

  MFloat& a_vortOtherPhase(const MInt cellId, const MInt dir) { return m_EEGas.vortOtherPhase[cellId][dir]; }

  MFloat a_nuTOtherPhase(const MInt cellId) const { return m_EEGas.nuTOtherPhase[cellId]; }

  MFloat& a_nuTOtherPhase(const MInt cellId) { return m_EEGas.nuTOtherPhase[cellId]; }

  MFloat a_nuEffOtherPhase(const MInt cellId) const { return m_EEGas.nuEffOtherPhase[cellId]; }

  MFloat& a_nuEffOtherPhase(const MInt cellId) { return m_EEGas.nuEffOtherPhase[cellId]; }

  /// \brief Returns the associatedBodyIds for fv-CellId \p cellId and \p set
  MInt& a_associatedBodyIds(const MInt cellId, const MInt set) { return m_associatedBodyIds[IDX_LSSETMB(cellId, set)]; }

  /// \brief Returns the associatedBodyIds for fv-CellId \p cellId and \p set
  MInt a_associatedBodyIds(const MInt cellId, const MInt set) const {
    return m_associatedBodyIds[IDX_LSSETMB(cellId, set)];
  }

  /// \brief Returns the curvature-value for fv-CellId \p cellId and \p set
  MFloat& a_curvatureG(const MInt cellId, const MInt set) { return m_curvatureG[set][cellId]; }

  /// \brief Returns the curvature-value for fv-CellId \p cellId and \p set
  MFloat a_curvatureG(const MInt cellId, const MInt set) const { return m_curvatureG[set][cellId]; }

  /// \brief Returns the flamespeed-value for fv-CellId \p cellId and \p set
  MFloat& a_flameSpeed(const MInt cellId, const MInt set) { return m_flameSpeedG[set][cellId]; }

  /// \brief Returns the flamespeed-value for fv-CellId \p cellId and \p set
  MFloat a_flameSpeed(const MInt cellId, const MInt set) const { return m_flameSpeedG[set][cellId]; }

  /// \brief Returns the noSets for fv-CellId \p cellId and \p set
  MInt& a_noSets() { return m_noSets; }

  /// \brief Returns the noSets for fv-CellId \p cellId and \p set
  MInt a_noSets() const { return m_noSets; }

  /// \brief Returns the noSets for fv-CellId \p cellId and \p set
  MInt& a_noLevelSetFieldData() { return m_noLevelSetFieldData; }

  /// \brief Returns the noSets for fv-CellId \p cellId and \p set
  MInt a_noLevelSetFieldData() const { return m_noLevelSetFieldData; }


  /// \brief Returns the Id of the split cell, if cellId is a split child.
  //         If cellId is a ghost cell, the function checks whether the boundary cell is a split child;
  //         If it is a split child, the function returns the split cell, otherwise it returns the associated internal
  //         cell
  const MInt& getAssociatedInternalCell(const MInt& cellId) const {
#ifndef NDEBUG
    std::ostringstream output;
    output << "cellId " << cellId << " Halo " << a_isHalo(cellId) << " IsPeriodic " << a_isPeriodic(cellId)
           << " bndryGhostCell " << a_isBndryGhostCell(cellId) << " splitChild "
           << a_hasProperty(cellId, SolverCell::IsSplitChild) << " splitCell "
           << a_hasProperty(cellId, SolverCell::IsSplitCell) << std::endl;
    if(a_isBndryGhostCell(cellId)) {
      ASSERT(cellId - m_bndryGhostCellsOffset < (signed)m_associatedInternalCells.size()
                 && cellId - m_bndryGhostCellsOffset >= 0,
             "size exceeds container, cellId: " << cellId << "m_bndryGhostCellsOffset: " << m_bndryGhostCellsOffset
                                                << " m_associatedInterNallCells.size(): "
                                                << m_associatedInternalCells.size() << output.str());
      ASSERT(m_associatedInternalCells[cellId - m_bndryGhostCellsOffset] >= 0
                 && m_associatedInternalCells[cellId - m_bndryGhostCellsOffset] < maxNoGridCells(),
             output.str());
    } else {
      ASSERT(m_splitChildToSplitCell.count(cellId) == 1, "mapped value not found for" << output.str());
      ASSERT(m_splitChildToSplitCell.find(cellId)->second < a_noCells(),
             "mapped value exceeds grid().tree().size() for cellId "
                 << cellId << " " << m_splitChildToSplitCell.find(cellId)->second << " " << a_noCells() << " "
                 << output.str());
      ASSERT(m_splitChildToSplitCell.find(cellId)->second > -1, "mapped value is negative for " << output.str());
    }
#endif
    return (
        !a_isBndryGhostCell(cellId)
            ? m_splitChildToSplitCell.find(cellId)->second
            : (a_hasProperty(m_associatedInternalCells[cellId - m_bndryGhostCellsOffset], SolverCell::IsSplitChild)
                   ? m_splitChildToSplitCell.find(m_associatedInternalCells[cellId - m_bndryGhostCellsOffset])->second
                   : m_associatedInternalCells[cellId - m_bndryGhostCellsOffset]));
  }

  virtual void reIntAfterRestart(MBool);
  virtual void resetRHS();
  virtual void resetRHSCutOffCells();
  virtual void initSolutionStep(MInt);
  virtual void applyInitialCondition();
  virtual void Ausm();
  virtual void Muscl(MInt = -1);
  virtual void viscousFlux();
  virtual void copyVarsToSmallCells();
  virtual void computePV();
  virtual void computePrimitiveVariables();
  virtual void filterConservativeVariablesAtFineToCoarseGridInterfaces();
  virtual void copyRHSIntoGhostCells();
  virtual void LSReconstructCellCenter_Boundary();
  virtual void LSReconstructCellCenter();
  virtual void applyBoundaryCondition();
  virtual void cutOffBoundaryCondition();
  virtual void computeConservativeVariables();
  virtual void initNearBoundaryExchange(const MInt mode = 0, const MInt offset = 0);
  void initAzimuthalNearBoundaryExchange(MIntScratchSpace& activeFlag);
  void azimuthalNearBoundaryExchange();
  void azimuthalNearBoundaryReverseExchange();
  void setActiveFlag(MIntScratchSpace&, const MInt mode, const MInt offset);
  virtual void setAdditionalActiveFlag(MIntScratchSpace&){};

  void setInfinityState();

  void initAzimuthalCartesianHaloInterpolation();

  // Detailed chemistry functions and variables
  template <class _ = void, std::enable_if_t<isDetChem<SysEqn>, _*> = nullptr>
  void computeSurfaceCoefficients();
  template <class _ = void, std::enable_if_t<isDetChem<SysEqn>, _*> = nullptr>
  void computeSpeciesReactionRates();
  template <class _ = void, std::enable_if_t<isDetChem<SysEqn>, _*> = nullptr>
  void computeMeanMolarWeights_PV();
  template <class _ = void, std::enable_if_t<isDetChem<SysEqn>, _*> = nullptr>
  void computeMeanMolarWeights_CV();
  template <class _ = void, std::enable_if_t<isDetChem<SysEqn>, _*> = nullptr>
  void setMeanMolarWeight_CV(MInt cellId);
  template <class _ = void, std::enable_if_t<isDetChem<SysEqn>, _*> = nullptr>
  void setMeanMolarWeight_PV(MInt cellId);
  template <class _ = void, std::enable_if_t<isDetChem<SysEqn>, _*> = nullptr>
  void computeGamma();
  template <class _ = void, std::enable_if_t<isDetChem<SysEqn>, _*> = nullptr>
  void computeDetailedChemistryVariables();
  template <class _ = void, std::enable_if_t<isDetChem<SysEqn>, _*> = nullptr>
  void setAndAllocateDetailedChemistryProperties();
  template <class _ = void, std::enable_if_t<isDetChem<SysEqn>, _*> = nullptr>
  void initCanteraObjects();
  template <class _ = void, std::enable_if_t<isDetChem<SysEqn>, _*> = nullptr>
  void correctMajorSpeciesMassFraction();

  void compute1DFlameSolution();
  void addSpeciesReactionRatesAndHeatRelease();
  void initHeatReleaseDamp();
  void computeAcousticSourceTermQe(MFloatScratchSpace&, MFloatScratchSpace&, MFloatScratchSpace&, MFloatScratchSpace&);

#if defined(WITH_CANTERA)
  std::shared_ptr<Cantera::Solution> m_canteraSolution;
  std::shared_ptr<Cantera::ThermoPhase> m_canteraThermo;
  std::shared_ptr<Cantera::Kinetics> m_canteraKinetics;
  std::shared_ptr<Cantera::Transport> m_canteraTransport;
#endif

  std::unique_ptr<OneDFlame> m_oneDimFlame;

  std::vector<MString> m_speciesName;
  std::map<std::string, MInt> speciesMap;
  MFloat* m_molarMass;
  MFloat* m_fMolarMass;
  MFloat* m_standardHeatFormation;
  MFloat* m_YInfinity;

  MBool m_detChemExtendedOutput;

  /// Indicator if sampling variables are initialized
  MBool m_isInitSamplingVars = false;
  static constexpr MInt s_maxNoSamplingVariables = 3;
  /// Storage for solver specific sampling variables
  std::array<MFloat**, s_maxNoSamplingVariables> m_samplingVariables{nullptr};
  /// Status of sampling variables to check if variable is already computed and exchanged
  std::array<MInt, 2 * s_maxNoSamplingVariables> m_samplingVariablesStatus{};

  MBool m_localTS = false;

  //------------------------------
  void computeVolumeForces();
  void computeRotForces();

  virtual MInt getAdjacentLeafCells_d0(const MInt, const MInt, MIntScratchSpace&, MIntScratchSpace&);
  virtual MInt getAdjacentLeafCells_d1(const MInt, const MInt, MIntScratchSpace&, MIntScratchSpace&);
  virtual MInt getAdjacentLeafCells_d2(const MInt, const MInt, MIntScratchSpace&, MIntScratchSpace&);
  virtual MInt getAdjacentLeafCells_d0_c(const MInt, const MInt, MIntScratchSpace&, MIntScratchSpace&);
  virtual MInt getAdjacentLeafCells_d1_c(const MInt, const MInt, MIntScratchSpace&, MIntScratchSpace&);
  virtual MInt getAdjacentLeafCells_d2_c(const MInt, const MInt, MIntScratchSpace&, MIntScratchSpace&);

  MFloat computeRecConstSVD(const MInt cellId, const MInt offset, MFloatScratchSpace& tmpA, MFloatScratchSpace& tmpC,
                            MFloatScratchSpace& weights, const MInt recDim, const MInt, const MInt,
                            const std::array<MBool, nDim> dirs = {}, const MBool relocateCenter = false);
  void extendStencil(const MInt);
  MInt samplingInterval();
  void checkGhostCellIntegrity();
  virtual void initializeRungeKutta();
  virtual void initializeMaxLevelExchange();
  virtual void computePrimitiveVariablesCoarseGrid();
  void initAzimuthalMaxLevelExchange();
  void finalizeMpiExchange();
  void setUpwindCoefficient();
  void cellSurfaceMapping();
  template <class _ = void, std::enable_if_t<isDetChem<SysEqn>, _*> = nullptr>
  void interpolateSurfaceDiffusionFluxOnCellCenter(MFloat* const, MFloat* const);
  void setNghbrInterface();

  void calcLESAverage();
  void saveLESAverage();
  void loadLESAverage();
  void finalizeLESAverage();
  // this factor is an empirical factor, should be changed based on setup and operating condition
  MFloat getAveragingFactor() { return F1B8 * m_Ma * sqrt(m_TInfinity) * timeStep(); };
  void saveSpongeData();
  void loadSpongeData();
  void initSTGSpongeExchange();
  void setAndAllocateZonalProperties();
  void readPreliminarySTGSpongeData();
  void calcPeriodicSpongeAverage();
  void exchangeZonalAverageCells();

  virtual void initSTGSponge();
  virtual void resetZonalLESAverage();
  virtual void determineLESAverageCells();
  virtual void resetZonalSolverData();
  virtual void getBoundaryDistance(MFloatScratchSpace&);
  virtual void nonReflectingBCAfterTreatmentCutOff();
  virtual void nonReflectingBCCutOff();
  virtual void dqdtau();
  virtual bool rungeKuttaStep();

  virtual void applyExternalSource();
  virtual void applyExternalOldSource(){};
  virtual void advanceExternalSource(){};
  void exchangeExternalSources();
  void resetExternalSources();

  MInt setUpBndryInterpolationStencil(const MInt, MInt*, const MFloat*);

  //  void computeMeanOutletPressure();
  void deleteSrfcs();
  virtual void resetRHSNonInternalCells();
  virtual void correctMasterCells();
  //  MBool triangleInterpolation(MFloat**, MFloat*);
  //  MBool tetraederInterpolation(MFloat**, MFloat*);
  //  void bilinearInterpolationConstants(MFloat**, MFloat*);
  //  void bilinearInterpolationConstants2(MFloat*, MFloat*, MFloat*, MFloat*);
  //  void trilinearInterpolationConstants(MFloat**, MFloat*);
  //  void trilinearInterpolationConstants2(MFloat*, MFloat*, MFloat*, MFloat*);
  void writeCellData(MInt);
  virtual void convertPrimitiveRestartVariables();
  //! FV Constructor: reads and allocate properties/variables:
  void initializeFvCartesianSolver(const MBool* propertiesGroups);
  void copyGridProperties();
  void allocateCommunicationMemory();
  void setTestcaseProperties();
  void setSamplingProperties();
  void setInputOutputProperties();
  void setNumericalProperties();
  void setAndAllocateCombustionTFProperties();
  void setAndAllocateSpongeLayerProperties();
  void allocateAndInitSolverMemory();
  void setRungeKuttaProperties();
  void setAndAllocateAdaptationProperties();

  //  virtual // parallel
  //  void exchangeTimeStep() = 0;
  //  virtual void setPreviousTimeStep() = 0;
  //  virtual void setLocalTimeStep() = 0;
  virtual void exchange();

  template <typename T>
  void exchangeDataFV(T* data, const MInt blockSize = 1, MBool cartesian = true,
                      const std::vector<MInt>& rotIndex = std::vector<MInt>());
  template <MBool exchangeAll_ = true>
  void exchangeFloatDataAzimuthal(MFloat* data, MInt noVars, const std::vector<MInt>& rotIndices);
  template <typename T>
  void exchangeDataAzimuthal(T* data, const MInt dataBlockSize = 1);
  void exchangeAzimuthalRemappedHaloCells();
  virtual void exchangePeriodic();
  void exchangePipe();
  void startMpiExchange();
  void prepareMpiExchange();
  void finishMpiExchange();
  void cancelMpiRequests() override;
  template <MBool exchangeAll_>
  void gather();
  void receive(const MBool exchangeAll = false);
  void send(const MBool exchangeAll = false);
  template <MBool exchangeAll_>
  void scatter();
  void exchangeAll();
  virtual void computeReconstructionConstants();
  virtual void findNghbrIds();
  void getPrimitiveVariables(MInt, MFloat*, MFloat*, MInt);
  void rhs();
  void rhsBnd();
  void initSolver() override;
  void finalizeInitSolver() override;
  void preSolutionStep(MInt) override{};
  MBool solutionStep() override;
  MBool postSolutionStep() override { return true; };
  MBool solverStep();
  void postTimeStep() override;
  void scalarLimiter();
  void cleanUp() override {
    // Finalize communication
    finalizeMpiExchange();
  };

  void lhsBndFinish();
  void lhsBnd();

  virtual void smallCellCorrection(const MInt timerId = -1);
  virtual void smallCellRHSCorrection(const MInt timerId = -1);
  virtual void updateSplitParentVariables(){/*only do something if we are using MB*/};
  virtual void checkDiv(){};
  virtual void updateMaterialNo(){};

  template <class _ = void, std::enable_if_t<isEEGas<SysEqn>, _*> = nullptr>
  void initSourceCells();
  void revertTimestep();
  template <class _ = void, std::enable_if_t<isEEGas<SysEqn>, _*> = nullptr>
  void rhsEEGas();
  virtual void resetImplicitCoefficients();

  MFloat physicalTime() { return m_physicalTime; }

  MFloat computeDomainLength(MInt direction);

 private:
  /// Convergence status of the current time step
  MBool m_timeStepConverged = false;

  MBool m_trackMovingBndry{};
  MInt m_trackMbStart{};
  MInt m_trackMbEnd{};

  MBool m_forceAdaptation = false;
  MInt m_lastAdapTS = 0;

 public:
  List<MInt>* m_sortedPeriodicCells = nullptr;
  MInt m_totalnosplitchilds = 0;
  MInt m_totalnoghostcells = 0;
  MBool m_constructGField{};
  MBool m_deleteNeighbour = false;
  //  MInt m_noNbBackup = 0;
  //  List< MInt >* m_nbBackupCellIds;
  //  List< MInt >* m_nbBackupDirs;
  //  List< MInt >* m_nbBackupNghbrIds;

  MBool m_bndryLevelJumps = false;
  MInt m_lsCutCellBaseLevel;
  MInt m_noOuterBndryCells = 0;

  MBool m_refineDiagonals = true;

  using Geom = Geometry<nDim>;

  MFloat* m_sweptVolume = nullptr;
  MFloat* m_sweptVolumeBal = nullptr;

  MBool m_engineSetup = false;

 public:
  /// [Splitt] The following is part of a first step to splitt CartesianGrid
  /// from the inheritance hierarchy:
  ///
  /// - in order to avoid renaming a lot of access to CartesianGrid data
  ///   members, references are introduced:
  ///
  /// \todo labels:FV,toremove this references will be removed in future commits
  Collector<PointBasedCell<nDim>>* m_extractedCells = nullptr;
  Collector<CartesianGridPoint<nDim>>* m_gridPoints = nullptr;
  /// the references to CartesianGrid data members end here

  /// Access the solver's geometry
  const Geom& geometry() const { return *m_geometry; }


  /// Return the global MPI communicator used by the grid
  MPI_Comm globalMpiComm() const { return grid().raw().mpiComm(); }

  // Creates a 2D slice from a 3D grid.
  void createGridSlice(const MString& direction, const MFloat intercept, const MString& fileName,
                       MInt* const sliceCellIds) {
    grid().raw().createGridSlice(direction, intercept, fileName, -1, nullptr, sliceCellIds, nullptr, nullptr);
  }

 protected:
  /// Access the solver's geometry (non-const version)
  Geom& geometry() { return *m_geometry; }

  // Collector for FV surfaces
  FvSurfaceCollector m_surfaces;
  // Access surface collector
  FvSurfaceCollector& m_surfaceCollector() { return m_surfaces; }

  /// Collector for FV cells
  maia::fv::collector::FvCellCollector<nDim> m_cells;

  MFloat m_maRot = F0;
  MFloat m_MaHg;
  MFloat m_THg;
  MFloat m_UHg;
  MFloat m_VHg;
  MFloat m_WHg;
  MFloat m_PHg;
  MFloat m_rhoHg;
  MFloat m_MaCg;
  MFloat m_TCg;
  MFloat m_UCg;
  MFloat m_VCg;
  MFloat m_WCg;
  MFloat m_PCg;
  MFloat m_rhoCg;

  // cell and surface lists
  MInt* m_bndryRfnJumpInformation = nullptr;
  MInt* m_bndryRfnJumpInformation_ = nullptr;
  MInt m_sweepStartFirstCell;
  MInt* m_activeCellIds = nullptr;
  MInt m_noActiveCells;
  MInt m_noActiveHaloCellOffset;
  MInt* m_cellsInsideSpongeLayer = nullptr;
  MInt m_noCellsInsideSpongeLayer;
  MInt* m_smallCellIds = nullptr;
  MInt* m_masterCellIds = nullptr;
  MInt** m_maxLevelWindowCells = nullptr;
  MInt* m_noMaxLevelWindowCells = nullptr;
  MInt** m_maxLevelHaloCells = nullptr;
  MInt* m_noMaxLevelHaloCells = nullptr;
  MInt m_slopeMemory;      // slope offsets
  MInt m_surfaceVarMemory; // surface offsets
  std::set<MInt> m_splitSurfaces;
  std::vector<MInt> m_splitCells;
  std::vector<std::vector<MInt>> m_splitChilds;
  std::set<MInt> m_cutOffInterface;

  // reconstruction arrays
  MFloat** m_A = nullptr;
  MFloat** m_ATA = nullptr;
  MFloat** m_ATAi = nullptr;
  MInt* m_reconstructionDataPeriodic = nullptr;
  MUint m_reconstructionDataSize;
  std::vector<MFloat> m_reconstructionConstants;
  std::vector<MInt> m_reconstructionCellIds;
  std::vector<MInt> m_reconstructionNghbrIds;
  MFloat* m_reconstructionConstantsPeriodic = nullptr;
  // Azimuthal
  std::vector<std::vector<MInt>> m_azimuthalMaxLevelHaloCells;
  std::vector<std::vector<MInt>> m_azimuthalMaxLevelWindowCells;
  std::vector<std::vector<MInt>> m_azimuthalRemappedHaloCells;
  std::vector<std::vector<MInt>> m_azimuthalRemappedWindowCells;
  std::vector<MInt> m_azimuthalRemappedNeighborDomains;
  std::vector<MInt> m_azimuthalRemappedNeighborsDomainIndex;
  MBool m_azimuthalRecConstSet = false;
  MBool m_azimuthalNearBndryInit = false;
  const MInt m_maxNoAzimuthalRecConst = 250;
  MBool m_planeInterp = false;
  std::vector<MInt> m_noAzimuthalReconstNghbrs;
  std::vector<MFloat> m_azimuthalRecConsts;
  std::vector<MInt> m_azimuthalReconstNghbrIds;
  std::vector<MInt> m_azimuthalBndrySide;
  std::vector<MFloat> m_azimuthalCutRecCoord;
  std::vector<std::vector<MInt>> m_azimuthalMaxLevelWindowMap;
  std::vector<MInt> m_azimuthalHaloActive;
  const MInt m_azimuthalNearBoundaryBackupMaxCount =
      8; // Estimated to be enough. If 8 is exceeded mTerm is called in localToGlobalIds()
  MFloat m_azimuthalAngle;
  std::vector<MInt> m_rotIndVarsPV;
  std::vector<MInt> m_rotIndVarsCV;

  // variables
  MFloat m_angularBodyVelocity;
  MFloat m_surfaceTangentialVelocity;
  MInt* m_secondBodyId = nullptr;
  MFloat** m_rhs0 = nullptr;
  MFloat* m_volumeAcceleration = nullptr;
  MFloat* m_rotAxisCoord = nullptr;
  MFloat* m_heatRelease = nullptr;
  // limiter
  MFloat* m_limPhi = nullptr;
  std::vector<MFloat> m_dampFactor;

  MString m_reactionScheme;
  MFloat m_hInfinity;
  MFloat m_gasConstant;
  MFloat m_thickeningFactor;
  MFloat m_referenceDensityTF;
  MFloat m_heatReleaseReductionFactor;
  MInt m_temperatureChange;
  MFloat m_burntUnburntTemperatureRatio;
  MFloat m_burntUnburntTemperatureRatioEnd;
  MFloat m_burntUnburntTemperatureRatioStart;
  MFloat* m_molecularWeight = nullptr;
  MFloat* m_FmolecularWeight = nullptr;
  MFloat* m_molarFormationEnthalpy = nullptr;
  MFloat* m_formationEnthalpy = nullptr;
  MFloat* m_referenceComposition = nullptr;
  MFloat* m_secondaryReferenceComposition = nullptr;

  // Jet model
  MBool m_jet = false;
  MBool m_jetForcing = false;
  MFloat m_jetForcingPosition;
  MFloat m_jetRandomSeed;
  MFloat m_jetHalfWidth;
  MFloat m_jetCoflowOffset;
  MFloat m_jetCoflowEndOffset;
  MFloat m_jetHalfLength;
  MInt m_noJetConst;
  MFloat* m_jetConst = nullptr;
  MFloat m_forceCoefficient = 0.0;
  MFloat m_densityRatio;
  MFloat m_shearLayerThickness;
  MFloat m_MaCoflow;
  MFloat m_jetTemperature = -1;
  MFloat m_jetDensity = -1;
  MFloat m_jetPressure = -1;
  MFloat m_jetHeight = 0.5;
  MFloat m_primaryJetRadius;
  MFloat m_secondaryJetRadius;
  MInt m_modeNumbers;
  MFloat m_targetVelocityFactor;
  MFloat m_momentumThickness = 0.0;
  MInt m_jetType;

  // Jet from a (chevron) nozzle
  MBool m_chevron = false;
  MFloat m_inletRadius = -1.0;
  MFloat m_outletRadius = -1.0;
  MFloat m_normJetTemperature = -1.0; // Normalized jet temperature at nozzle exit
  MFloat m_maNozzleExit = -1.0;       // Acoustic Mach number at nozzle exit: Ma_acoustic = u_jet/c_inf
  // Nozzle exit conditions
  MFloat m_nozzleExitMaJet = -1.0; // Local Mach number at nozzle exit: Ma_jet = u_jet/c_jet
  MFloat m_nozzleExitTemp = -1.0;
  // MFloat m_nozzleExitP = -1.0; // assumed to be the same as the ambient pressure
  MFloat m_nozzleExitRho = -1.0;
  MFloat m_nozzleExitU = -1.0;
  // Nozzle inlet conditions
  MFloat m_maNozzleInlet = -1.0;
  MFloat m_nozzleInletTemp = -1.0;
  MFloat m_nozzleInletP = -1.0;
  MFloat m_nozzleInletRho = -1.0;
  MFloat m_nozzleInletU = -1.0;

  // level set / progress variable model
  MFloat m_c0;
  MFloat m_laminarFlameThickness;
  MFloat m_subfilterVariance;
  MFloat m_maxReactionRate;

  MFloat m_MaFlameTube;
  MFloat m_temperatureFlameTube;
  MFloat m_velocityFlameTube;
  MFloat m_rhoFlameTube;
  MFloat m_rhoUnburnt;
  MFloat m_rhoBurnt;
  MFloat m_pressureFlameTube;
  MFloat m_pressureUnburnt;
  MFloat m_inletTubeAreaRatio;
  MFloat m_flameOutletAreaRatio;
  MFloat m_inletOutletAreaRatio;
  MBool m_twoFlames;
  MFloat m_dampingDistanceFlameBase;
  MFloat m_dampingDistanceFlameBaseExtVel;
  MFloat m_realRadiusFlameTube;
  MFloat m_radiusVelFlameTube;
  MFloat m_radiusInjector;
  MFloat m_yOffsetInjector;
  MFloat m_radiusFlameTube;
  MFloat m_radiusFlameTube2;
  MFloat m_initialFlameHeight;
  MFloat m_xOffsetFlameTube;
  MFloat m_xOffsetFlameTube2;
  MFloat m_yOffsetFlameTube;
  MFloat m_yOffsetFlameTube2;
  MFloat m_deltaXtemperatureProfile;
  MFloat m_deltaYtemperatureProfile;
  MFloat m_thermalProfileStartFactor;
  MFloat m_flameRadiusOffset;
  MFloat m_shearLayerStrength;
  MFloat m_ScT;
  MFloat m_NuT;
  MFloat m_integralAmplitude;
  MFloat m_integralLengthScale;

  // dimensions of domain and sponge zones
  MInt m_spongeLayerLayout;
  MFloat* m_domainBoundaries = nullptr;
  MFloat* m_spongeCoord = nullptr;

  MBool m_levelSet = false;
  MBool m_levelSetMb = false; // Moving boundary extension
  MBool m_LsRotate = false; // Rotating levelset
  MBool m_levelSetRans = false;
  MBool m_combustion = false;
  // MBool m_isDetChem = false;
  MBool m_LSSolver;
  MBool m_acousticAnalysis;
  MBool m_thickenedFlame = false;

  std::vector<MFloat>* m_levelSetValues = nullptr;
  MFloat* m_levelSetValuesMb = nullptr;
  MInt* m_associatedBodyIds = nullptr;
  MInt m_noLevelSetsUsedForMb{};
  MInt m_noLevelSetFieldData;
  std::vector<MFloat>* m_curvatureG = nullptr;
  std::vector<MFloat>* m_flameSpeedG = nullptr;
  MInt m_noSets = 0;
  MBool m_reComputedBndry = false;

  // restart/output related properties
  MString m_currentGridFileName;
  MBool m_adaptationSinceLastRestart;
  MBool m_adaptationSinceLastRestartBackup;
  MBool m_forceRestartGrid;
  MBool m_force1DFiltering;
  MBool m_gridConvergence;
  MBool m_gridInterfaceFilter;
  MBool m_totalDamp;
  MBool m_heatReleaseDamp;
  MBool m_useCorrectedBurningVelocity;
  MBool m_modelCheck;
  MBool m_recordBodyData;
  MBool m_recordLandA;
  MBool m_recordPressure;
  MBool m_vtkTest;
  MBool m_recordFlameFrontPosition;
  MBool m_recordWallVorticity;
  MBool m_structuredFlameOutput;
  MBool m_surfDistParallel;
  MBool m_surfDistCartesian;
  MInt m_writeOutData;
  MBool m_writeCutCellsToGridFile;
  MBool m_restartOldVariables = false;
  MBool m_restartOldVariablesReset = false;
  MBool m_bodyIdOutput;
  MBool m_levelSetOutput;
  MBool m_isActiveOutput;
  MBool m_domainIdOutput = false;
  MBool m_multipleFvSolver = false;
  const MChar** m_variablesName;
  const MChar** m_vorticityName;
  MBool m_saveVorticityToRestart = false;
  MBool m_vorticityOutput;
  MInt m_vorticitySize;
  MBool m_qCriterionOutput;
  MBool m_vtuWritePointData;
  MBool m_vtuWriteGeometryFile;
  MBool m_vtuWriteParticleFile;
  MBool m_vtuCutCellOutput;
  std::set<MInt> m_vtuGeometryOutput;
  MBool m_vtuGlobalIdOutput;
  MBool m_vtuDomainIdOutput;
  MBool m_vtuDensityOutput;
  MBool m_vtuLevelSetOutput;
  MBool m_vtuQCriterionOutput;
  MBool m_vtuLambda2Output;
  MBool m_vtuVorticityOutput;
  MBool m_vtuVelocityGradientOutput;
  MBool m_vtuGeometryOutputExtended;
  MBool m_vtuSaveHeaderTesting;
  MInt m_vtuLevelThreshold;
  MFloat* m_vtuCoordinatesThreshold = nullptr;
  MBool m_checkCellSurfaces;
  MBool m_considerVolumeForces;
  MBool m_considerRotForces = false;
  MString m_outputFormat;
  MBool m_allowInterfaceRefinement;
  MInt m_bndryCellSurfacesOffset;
  MInt m_bndrySurfacesOffset;
  MInt m_cellToRecordData;
  MInt m_counterCx;
  MInt m_dragOutputInterval;

  MBool m_integratedHeatReleaseOutput;
  MInt m_integratedHeatReleaseOutputInterval;

  MBool m_dualTimeStepping;
  MBool m_euler;
  MInt m_bndryGhostCellsOffset;
  MInt m_initialCondition;
  MInt m_limiter;
  MInt m_maxNoTimeSteps;
  MInt m_maxNoSurfaces;
  MInt m_noGNodes;
  MInt m_noRKSteps; // number of Runge-Kutta steps
  MInt m_noSamples; // number of samples used for averaging
  MInt m_maxIterations;
  MInt m_orderOfReconstruction;
  MInt m_restartBackupInterval;
  MBool m_restartBc2800;
  MFloat m_restartTimeBc2800;
  MInt m_RKStep;
  MInt m_rungeKuttaOrder;
  MInt m_noTimeStepsBetweenSamples;
  MInt m_forceNoTimeSteps;
  MInt m_structuredFlameOutputLevel;
  // MFloat m_weightTauC;
  // MFloat m_weightTauE;
  MFloat m_adaptationDampingDistance;
  MFloat* m_angle = nullptr;
  MFloat* m_RKalpha = nullptr; // Runge-Kutta coefficients friend class FvCartesianSolver2D;
                               // friend class FvCartesianSolver3D;
                               // friend class FvMbSolver2D;
                               // friend class FvMbSolver3D;
                               // friend class LsCartesianSolver<nDim>;

  MBool m_isEEGas = false;
  struct {
    MFloat RKSemiImplicitFactor;
    MFloat** uOtherPhase = nullptr;
    MFloat** uOtherPhaseOld = nullptr;
    MFloat** gradUOtherPhase = nullptr;
    MFloat** vortOtherPhase = nullptr;
    MFloat* nuTOtherPhase = nullptr;
    MFloat* nuEffOtherPhase = nullptr;
    MFloat Eo0;
    MFloat bubbleDiameter;
    MFloat CD;
    MFloat CL;
    MInt dragModel;
    MFloat liquidDensity;
    MFloat eps;
    MInt gasSource;
    MFloat gasSourceMassFlow;
    std::vector<MInt> gasSourceCells;
    MInt noGasSourceBoxes;
    std::vector<MFloat> gasSourceBox;
    MBool bubblePathDispersion;
    MFloat massSource;
    MFloat initialAlpha;
    MFloat alphaInf;
    MFloat alphaIn;
    MFloat schmidtNumber;
    MBool uDLimiter;
    MFloat uDLim;
    MBool depthCorrection;
    std::vector<MFloat> gravity;
    std::vector<MFloat> gravityRefCoords;
    std::vector<MFloat> depthCorrectionCoefficients;
    MFloat interpolationFactor;
  } m_EEGas;

  struct {
    MFloat infTemperature;
    MFloat infPressure;
    MFloat infPhi;
    MString* infSpeciesName;
    MFloat* infSpeciesMassFraction;
    MFloat infVelocity;
    MFloat laminarFlameSpeedFactor;

    MBool hasChemicalReaction;
    MString reactionMechanism;
    MString phaseName;
    MString transportModel;
    MBool soretEffect;
  } m_detChem;

  MInt* m_storeNghbrIds = nullptr;
  MInt* m_identNghbrIds = nullptr;

  // Zonal
  MBool m_zonal = false;
  MInt m_zonalRestartInterpolationSolverId;
  MBool m_resetInitialCondition;
  MBool m_rans;
  MInt m_noRansEquations = 0;
  MFloat m_turbulenceDegree;
  MFloat m_ransTransPos;
  MInt m_zonalAveragingTimeStep = 0;
  MInt m_zonalTransferInterval = 1;

  MInt m_noRANSVariables = -1;
  MInt m_noLESVariables = -1;
  std::vector<MFloat>* m_RANSValues = nullptr;
  std::vector<MFloat>* m_LESValues = nullptr;

  // LES average for zonal coupling
  std::vector<MFloat>* m_LESVarAverage = nullptr;
  std::vector<MFloat>* m_LESVarAverageBal = nullptr;
  std::vector<MInt> m_LESAverageCells;
  MInt m_LESNoVarAverage = 0;
  std::vector<MFloat> m_averagePos;
  std::vector<MInt> m_averageDir;
  std::vector<MBool> m_averageReconstructNut;

  // LES average
  MBool m_calcLESAverage = false;
  MBool m_restartLESAverage = false;
  MInt m_averageStartTimeStep = 0;

  // Nut reconstruct
  // MInt m_noReconstructNutVars;
  MInt m_rntStartTimeStep = 0;

  // STG
  MString m_bc7909RANSSolverType;
  MInt m_stgStartTimeStep = 0;
  MBool m_stgIsActive;
  MFloat m_pressureRatioChannel;
  MInt m_pressureRatioStartTimeStep;
  MInt m_pressureRatioEndTimeStep;
  MInt m_spongeTimeVelocity;
  MFloat** m_stgEddieCoverage = nullptr;

  // STG sponge
  MInt m_stgSpongeTimeStep = 0;
  MFloat* m_stgSpongePositions;
  MInt m_noStgSpongePositions;
  MBool m_STGSponge = false;
  // MBool m_preliminarySTGSponge = false;
  MBool m_preliminarySponge = false;
  MInt m_spongeRoot = -1;
  MFloat m_7901Position;
  MInt m_7901faceNormalDir;
  MInt m_7901wallDir;
  MInt m_7901periodicDir;
  std::vector<MFloat>* m_STGSpongeFactor = nullptr;
  MFloat** m_LESPeriodicAverage = nullptr;
  MFloat* m_uvErr = nullptr;
  MFloat* m_uvRans = nullptr;
  MFloat* m_uvInt = nullptr;
  MFloat m_spongeLimitFactor = 10000.0;
  MPI_Comm m_spongeComm;
  MInt m_spongeCommSize = -1;
  MInt m_spongeRank = -1;

  MInt m_noSpongeCells;
  MInt m_globalNoSpongeLocations = 0;
  MInt* m_globalNoPeriodicExchangeCells = nullptr;
  std::vector<MInt>* m_spongeAverageCellId = nullptr;
  std::vector<MInt> m_spongeCells;
  std::vector<MFloat> m_spongeLocations;
  MFloat* m_globalBcStgLocationsG = nullptr;
  // std::vector<MFloat> m_globalSpongeLocations;
  std::vector<std::pair<MFloat, MFloat>> m_globalSpongeLocations;

  // RANS
  MFloat m_tkeFactor;
  MFloat m_kInfinityFactor;
  MFloat m_omegaInfinityFactor;

  //================================== WMLES ======================================


  std::vector<FvWMSurface<nDim>> m_wmSurfaces;

  MInt m_wmSurfaceProbeInterval = 0;
  MInt m_wmGlobalNoSrfcProbeIds;
  MInt m_wmDomainId;
  MInt m_wmNoDomains;
  MBool m_wmLES = false;
  MBool m_wmOutput = false;
  MBool m_wmTimeFilter = false;
  MBool m_wmUseInterpolation = true;
  MFloat m_wmDistance;

  MInt m_wmIterator = 0;

  MPI_Comm m_comm_wm;

  MInt* m_wmLocalNoSrfcProbeIds = nullptr;
  MInt* m_noWMImgPointsSend = nullptr;
  MInt* m_noWMImgPointsRecv = nullptr;
  MInt* m_wmImgRecvIdMap = nullptr;
  MFloat* m_wmImgSendBuffer = nullptr;
  MFloat* m_wmImgRecvBuffer = nullptr;
  MFloat* m_wmSrfcProbeSendBuffer = nullptr;
  MFloat* m_wmSrfcProbeRecvBuffer = nullptr;
  MPI_Request* m_mpi_wmRequest = nullptr;
  MPI_Request* m_mpi_wmSendReq = nullptr;
  MPI_Request* m_mpi_wmRecvReq = nullptr;

  std::vector<std::vector<MInt>> m_wmImgCellIds;
  std::vector<std::vector<MInt>> m_wmImgWMSrfcIds;
  std::vector<std::vector<MFloat>> m_wmImgCoords;

  std::vector<MInt> m_wmSurfaceProbeIds;
  std::vector<MInt> m_wmSurfaceProbeSrfcs;

  //============================ Wall normal Output  =======================
  MFloat m_normalLength;
  MInt m_normalNoPoints;
  MInt m_normalOutputInterval;
  MInt m_normalOutputInitCounter = 0;
  MInt m_normalBcId;
  MBool m_wallNormalOutput;
  MBool m_useWallNormalInterpolation;

  std::vector<MFloat> m_normalSamplingCoords;
  std::vector<MInt> m_normalSamplingSide;

  MInt m_noWallNormals = 0;
  std::vector<MFloat> m_wallNormalPointCoords; // computeWallNormalPointCoords
  std::vector<MFloat> m_wallNormalVectors;

  std::vector<MFloat> m_wallSetupOrigin;
  std::vector<MInt> m_wallSetupOriginSide;

  std::vector<MInt> m_wallNormalPointDomains; // findWallNormalCellIds
  std::vector<MInt> m_wallNormalPointCellIDs; // findWallNormalCellIds
  std::vector<std::vector<MInt>> m_neighborPointIds;
  std::vector<MFloatTensor> m_interpolationMatrices;
  std::vector<MInt> m_interpolationPosition;
  void computeWallNormalPointCoords();
  void findWallNormalCellIds();
  MFloat interpolateWallNormalPointVars(MInt var, MFloat coords[], MInt localId, std::vector<MInt> neighborList);
  std::vector<MInt> findWallNormalNeighbors(MInt pointId);
  void getWallNormalPointVars();


  // Spanwise Averaged Surface Probes
  MInt m_saNoSrfcProbes;
  MInt m_saSrfcProbeInterval;
  MInt m_saSrfcProbeStart;
  MBool m_saSrfcProbes = false;
  MString m_saSrfcProbeDir;
  std::vector<std::vector<MInt>> m_saSrfcProbeIds;
  std::vector<std::vector<MInt>> m_saSrfcProbeSrfcs;

  MFloat* m_saSrfcProbeBuffer = nullptr;
  MInt* m_saSrfcProbeNoSamples = nullptr;

  void initSpanAvgSrfcProbes();
  void writeSpanAvgSrfcProbes();

  MFloat computeWMViscositySpalding(MInt);
  MFloat computeWMViscositySpalding3D(MInt);
  void initWMSurfaceProbes();
  void writeWMSurfaceProbes();
  void writeWMTimersASCII();
  void initWMExchange();
  void exchangeWMVars();
  void gatherWMVars();
  void receiveWMVars();
  void sendWMVars();
  void scatterWMVars();
  void readWallModelProperties();
  void restartWMSurfaces();

  MBool m_useSandpaperTrip;
  MBool m_useChannelForce;
  MFloat m_channelVolumeForce;

  void initChannelForce();
  void applyChannelForce();

 protected:
  MFloat m_chi;
  MInt m_upwindMethod;
  MInt m_reConstSVDWeightMode{};
  MBool m_relocateCenter;
  MBool m_reExcludeBndryDiagonals;
  MBool m_2ndOrderWeights;
  MFloat m_cfl;
  MFloat m_cflViscous;
  MFloat m_convergenceCriterion;
  MFloat m_deltaP;
  MFloat m_deltaPL;
  MFloat m_gamma = NAN;
  MFloat m_globalUpwindCoefficient;
  MFloat m_inflowTemperatureRatio;
  MFloat** m_kronecker = nullptr;
  MFloat m_massConsumption;
  MFloat m_maxTemp;
  MFloat m_meanPressure;
  MFloat m_meanY;
  MFloat m_physicalTime;
  MFloat m_physicalTimeStep;
  MFloat m_Pr;
  MFloat m_rPr;
  MFloat m_rRe0;
  MFloat m_referenceLength;
  MFloat m_previousMa;
  MBool m_changeMa;
  MBool m_useCreateCutFaceMGC;
  MFloat m_sampleRate;
  MFloat m_samplingTimeBegin;
  MFloat m_samplingTimeEnd;
  MInt m_spongeLayerType;
  MFloat m_sigmaSponge;
  MFloat m_sigmaSpongeInflow;
  MFloat m_spongeReductionFactor;
  MFloat* m_spongeFactor = nullptr;
  MFloat m_spongeLayerThickness;
  // new sponge parameters dependent on the chosen sponge boundary ids
  MBool m_createSpongeBoundary;
  MInt m_noSpongeFactors;
  MInt m_noSpongeBndryCndIds; ///< number of sponge boundary IDs
  MInt m_noMaxSpongeBndryCells;
  MFloat* m_sigmaSpongeBndryId = nullptr;
  MFloat* m_sigmaEndSpongeBndryId = nullptr;
  MInt* m_spongeDirections = nullptr;
  MInt* m_spongeBndryCndIds = nullptr;
  MFloat* m_spongeStartIteration = nullptr;
  MFloat* m_spongeEndIteration = nullptr;
  MInt* m_spongeTimeDependent = nullptr;
  MBool m_spongeTimeDep;
  MBool m_velocitySponge;
  MFloat m_spongeWeight;
  MFloat m_spongeBeta;
  MFloat m_targetDensityFactor;
  MBool m_outputPhysicalTime = false;
  MFloat m_time;
  MFloat m_timeRef;
  MFloat m_totalHeatReleaseRate;
  MInt** m_cellSurfaceMapping = nullptr;

  MUlong m_randomDeviceSeed;
  MBool m_useCentralDifferencingSlopes = false;

  // post-processing
  MFloat m_oldMomentOfVorticity;
  MFloat m_oldNegativeMomentOfVorticity;
  MFloat m_oldPositiveMomentOfVorticity;
  MFloat** m_vorticity = nullptr;
  MBool m_loadSampleVariables = false;

  // test case specific settings
  MString m_testCaseName;
  MFloat m_timeStepConvergenceCriterion;

  MFloat** m_externalSource = nullptr;
  MFloat** m_externalSourceDt1 = nullptr;
  MInt* m_noParts = nullptr;
  MBool m_hasExternalSource = false;
  std::map<MInt, std::vector<MFloat>> m_vapourData;

  MBool m_calcSlopesAfterStep = false;
  /// Stores whether this solver is part of a multilevel computation
  MBool m_multilevel = false;
  /// Stores whether this solver is the primary solver (i.e., it has the finshest mesh) of a multilevel computation
  MBool m_isMultilevelPrimary = false;
  /// Stores whether this solver is the lowest secondart solver (i.e., is has the coarsest mesh) of a multilevel
  /// computation
  MBool m_isLowestSecondary = false;


  MInt m_maxLevelBeforeAdaptation = -1;

  // relevant for POST-Processing
  MBool m_statisticCombustionAnalysis{};
  MBool m_averageVorticity = false;
  MInt m_movingAvgInterval = 0;
  MBool m_averageSpeedOfSound = false;
  MInt m_skewness{};
  MInt m_kurtosis{};
  std::set<MInt> m_activeMeanVars{};

 public:
  using PrimitiveVariables = typename SysEqn::PrimitiveVariables;

  typename SysEqn::PrimitiveVariables* PV{};
  typename SysEqn::ConservativeVariables* CV{};
  typename SysEqn::FluxVariables* FV{};       // these are the variables for which the fluxes are calculated
                                              // usually the same as conservative variables
  typename SysEqn::AdditionalVariables* AV{}; // these are additional variables that are saved in the cell collector
                                              // They might be used in SysEqn functions etc.
  typename SysEqn::SurfaceCoefficients* SC{}; // these are additional surface coefficients that are saved in the surface
                                              // collector. They might be used in SysEqn functions etc.
  static constexpr MBool hasAV = SysEqn::hasAV;
  static constexpr MBool hasSC = SysEqn::hasSC;

  /// Return the number of primitive variables
  MInt noVariables() const override { return PV->noVariables; };

  void releaseMemory();

  /// Return true if solver is part of a multilevel computation
  constexpr MBool isMultilevel() const { return m_multilevel; }

  /// Return true if solver is primary solver in multilevel computation
  constexpr MBool isMultilevelPrimary() const { return isMultilevel() && m_isMultilevelPrimary; }

  constexpr MBool isMultilevelLowestSecondary() const { return isMultilevel() && m_isLowestSecondary; }

  /// Designates solver as primary solver in multilevel computation
  void setMultilevelPrimary(const MBool state = true) { m_isMultilevelPrimary = state; }
  void setMultilevelSecondary(const MBool state = true) { m_isLowestSecondary = state; }
  constexpr MBool isZonal() const { return m_zonal; }
  void loadSampleVariables(MInt timeStep);
  void getSampleVariables(MInt cellId, const MFloat*& vars);
  void getSampleVariables(MInt const cellId, std::vector<MFloat>& vars);
  void getSampleVariableNames(std::vector<MString>& varNames) override;
  virtual void getSampleVarsDerivatives(const MInt cellId, const MFloat*& vars);
  MBool getSampleVarsDerivatives(const MInt cellId, std::vector<MFloat>& vars);
  void calculateHeatRelease();
  void getHeatRelease(MFloat*& heatRelease);
  virtual void getVorticity(MFloat* const vorticity);
  virtual void getVorticityT(MFloat* const vorticity);
  void oldPressure(MFloat* const p);
  virtual MFloat& vorticityAtCell(const MInt cellId, const MInt dir);
  virtual MFloat getBoundaryHeatFlux(const MInt cellId) const;
  MFloat time() const override;
  virtual void getDimensionalizationParams(std::vector<std::pair<MFloat, MString>>& dimParams) const;

  /// Required for sampling, for FV the index is already the cell id
  MInt getCellIdByIndex(const MInt index) { return index; }

  /// Return the leaf cell id containing the given point
  MInt getIdAtPoint(const MFloat* point, MBool NotUsed(globalUnique = false)) {
    return grid().findContainingLeafCell(point);
  }

  /// Read sampling related properties
  virtual void getSolverSamplingProperties(std::vector<MInt>& samplingVars, std::vector<MInt>& noSamplingVars,
                                           std::vector<std::vector<MString>>& samplingVarNames,
                                           const MString featureName = "") override;

  /// Initialize sampling variables/allocate memory
  virtual void initSolverSamplingVariables(const std::vector<MInt>& varIds,
                                           const std::vector<MInt>& noSamplingVars) override;

  /// Calculate sampling variables
  virtual void calcSamplingVariables(const std::vector<MInt>& varIds, const MBool exchange) override;

  virtual void initInterpolationForCell(const MInt cellId);
  virtual void calcSamplingVarAtPoint(const MFloat* point, const MInt id, const MInt sampleVarId, MFloat* state,
                                      const MBool interpolate = false) override;

  /// Return the current time step
  MInt getCurrentTimeStep() const override { return globalTimeStep; }

  /// Return if the restart time step can be determined from the restart file (for
  /// useNonSpecifiedRestartFile = true)
  MBool hasRestartTimeStep() const override { return true; }

  /// Determine the restart time step from the restart file (for useNonSpecifiedRestartFile = true)
  MInt determineRestartTimeStep() const override;

  virtual void finalizeInitEnthalpySolver(){};
  virtual void initMatDat(){};
  virtual void writeVtkXmlFiles(const MString, const MString, MBool, MBool){};

  /// Return if slopes should be calculated at after each step (not before)
  MBool calcSlopesAfterStep() { return m_calcSlopesAfterStep; };

  /// Apply coarse-level correction to RHS
  void applyCoarseLevelCorrection();

 protected:

  virtual void viscousFlux_Gequ_Pv();
  virtual void viscousFlux_Gequ_Pv_Plenum();
  virtual void updateJet();
  virtual void writeListOfActiveFlowCells();

  MFloat reduceData(const MInt cellId, MFloat* data, const MInt dataBlockSize = 1, const MBool average = true);

  void sensorEntropyGrad(std::vector<std::vector<MFloat>>& sensors, std::vector<std::bitset<64>>& sensorCellFlag,
                         std::vector<MFloat>& sensorWeight, MInt sensorOffset, MInt sen) override;
  void sensorEntropyQuot(std::vector<std::vector<MFloat>>& sensors, std::vector<std::bitset<64>>& sensorCellFlag,
                         std::vector<MFloat>& sensorWeight, MInt sensorOffset, MInt sen) override;
  void sensorVorticity(std::vector<std::vector<MFloat>>& sensors, std::vector<std::bitset<64>>& sensorCellFlag,
                       std::vector<MFloat>& sensorWeight, MInt sensorOffset, MInt sen) override;
  void sensorDerivative(std::vector<std::vector<MFloat>>& sensors, std::vector<std::bitset<64>>& sensorCellFlag,
                        std::vector<MFloat>& sensorWeight, MInt sensorOffset, MInt sen) override;
  void sensorSpecies(std::vector<std::vector<MFloat>>& sensors, std::vector<std::bitset<64>>& sensorCellFlag,
                     std::vector<MFloat>& sensorWeight, MInt sensorOffset, MInt sen) override;
  void sensorParticle(std::vector<std::vector<MFloat>>& sensors, std::vector<std::bitset<64>>& sensorCellFlag,
                      std::vector<MFloat>& sensorWeight, MInt sensorOffset, MInt sen) override;
  virtual void setCellProperties();
  void initSpongeLayer();
  void determineStructuredCells();
  void tagCellsNeededForSurfaceFlux();

  void writeCutCellsToGridFile();

  virtual MBool gridPointIsInside(MInt, MInt);

  // parallel IO library functions  NETCDF
  void saveGridFlowVarsPar(const MChar* fileName, MInt noTotalCells, MLong noInternalCells,
                           MFloatScratchSpace& variables, std::vector<MString>& dbVariablesName, MInt,
                           MIntScratchSpace& idVariables, std::vector<MString>& idVariablesName, MInt,
                           MFloatScratchSpace& dbParameters, std::vector<MString>& dbParametersName,
                           MIntScratchSpace& idParameters, std::vector<MString>& idParametersName,
                           const MInt* recalcIds);

  void computeVorticity3D(MFloat* const vorticity);
  void computeVorticity2D(MFloat* const vorticity);
  void computeVorticity3DT(MFloat* const vorticity);
  void computeQCriterion(MFloatScratchSpace& qCriterion);
  void loadRestartTime(const MChar* fileName, MInt& globalTimeStepInput, MFloat& timeInput, MFloat& physicalTimeInput);
  virtual void loadOldVariables(const MString& fileName);
  virtual void loadGridFlowVarsPar(const MChar* fileName);
  virtual void loadRestartMaterial(){};
  virtual void initCellMaterialNo(){};

  virtual inline MFloat entropy(MInt cellId) {
    return sysEqn().entropy(a_pvariable(cellId, PV->P), a_pvariable(cellId, PV->RHO));
  }

  // parallelization
  MFloat** m_sendBuffers = nullptr;
  MFloat** m_receiveBuffers = nullptr;
  MInt m_dataBlockSize = -1;
  MBool m_nonBlockingComm = false;
  MFloat** m_sendBuffersNoBlocking = nullptr;
  MFloat** m_receiveBuffersNoBlocking = nullptr;
  // private:
  MPI_Request* m_mpi_request = nullptr;
  MPI_Request* m_mpi_sendRequest = nullptr;
  MPI_Request* m_mpi_receiveRequest = nullptr;
  // Identify if there are open MPI recv/send requests
  MBool m_mpiRecvRequestsOpen = false;
  MBool m_mpiSendRequestsOpen = false;
  MBool m_splitMpiCommRecv = false;

  // jet

  void setAndAllocateJetProperties();


 public:
  // Combustion:
  void setAndAllocateCombustionGequPvProperties();
  void setAndAllocateSpongeBoundaryProperties();
  MFloat setAndAllocateSpongeDomainProperties(MFloat);
  void setCombustionGequPvVariables();

  // moving boundary
  MInt** m_setToBodiesTable = nullptr;
  MFloat* m_bodyCenter = nullptr;
  MFloat* m_bodyVelocity = nullptr;
  MFloat* m_bodyVelocityDt1 = nullptr;
  MFloat* m_bodyVelocityDt2 = nullptr;
  MFloat* m_bodyAcceleration = nullptr;
  MFloat* m_bodyAngularVelocity = nullptr;
  MFloat* m_bodyAngularAcceleration = nullptr;
  MFloat* m_bodyTemperature = nullptr;
  MFloat* m_bodyTemperatureDt1 = nullptr;
  MFloat* m_bodyHeatFlux = nullptr;
  MInt m_volumeForcingDir = -1;
  MFloat m_pipeRadius = -1;
  MInt m_noEmbeddedBodies;
  MInt m_noPeriodicGhostBodies;
  MInt* m_internalBodyId = nullptr;
  MInt m_levelSetAdaptationScheme = 0;

  // gap-Handling
  class FvGapCell {
   public:
    static constexpr MInt maxNoOverlappingBodies = 2;
    MInt cellId{};
    MInt region{};
    MInt status{};
    MInt bodyIds[maxNoOverlappingBodies]{};
    MFloat surfaceVelocity[3]{};

    FvGapCell(MInt cId, MInt rId, MInt sId, MInt body1, MInt body2) : cellId(cId), region(rId), status(sId) {
      std::fill_n(bodyIds, maxNoOverlappingBodies, -1);
      bodyIds[0] = body1;
      bodyIds[1] = body2;
      std::fill_n(surfaceVelocity, 3, 0);
    };
    ~FvGapCell() = default;
  };
  MBool m_closeGaps = false;
  MInt m_gapInitMethod = 2;
  MInt m_noGapRegions;
  std::vector<MInt> m_gapCellId;
  std::vector<FvGapCell> m_gapCells;

  // periodc bc
  MInt m_periodicCells;
  MFloat** m_periodicDataToSend = nullptr;
  MFloat** m_periodicDataToReceive = nullptr;
  MInt* m_noPerCellsToSend = nullptr;
  MInt* m_noPerCellsToReceive = nullptr;
  MInt* m_noPeriodicCellsDom = nullptr;
  MFloat** m_periodicCellDataDom = nullptr;
  MInt m_noPeriodicData;
  MInt m_noPeriodicCellData;

  MFloat m_oldPressure_Gradient;
  MFloat m_oldUbulk;
  MFloat m_target_Ubulk;
  MInt m_oldTimeStep;
  MFloat UbulkDiff;

  MFloat m_referenceTemperature;

  MFloat m_sutherlandConstant = NAN;
  MFloat m_sutherlandConstantThermal;
  MFloat m_sutherlandPlusOne = NAN;
  MFloat m_sutherlandPlusOneThermal;

 protected:
  Collector<FvBndryCell<nDim, SysEqn>>* m_bndryCells = nullptr;
  std::vector<MInt> m_associatedInternalCells;
  std::map<MInt, MInt> m_splitChildToSplitCell;
  MString m_surfaceValueReconstruction;
  MString m_viscousFluxScheme;
  MString m_advectiveFluxScheme;
  MFloat m_enhanceThreePointViscFluxFactor = 0.1;
  MInt m_noLimitedSlopesVar;
  MInt* m_limitedSlopesVar = nullptr;
  MInt m_computeExtVel{};
  MBool m_massFlux;
  MBool m_plenumWall;
  MBool m_plenum;
  MBool m_confinedFlame;
  MBool m_filterFlameTubeEdges;
  MFloat m_filterFlameTubeEdgesDistance;
  MBool m_divergenceTreatment;
  MBool m_specialSpongeTreatment;
  MInt m_constantFlameSpeed;
  MFloat m_flameSpeed{};
  MFloat m_turbFlameSpeed;
  MFloat m_Da;
  MFloat m_noReactionCells;
  MFloat m_velocityOutlet;
  MFloat m_analyticIntegralVelocity;
  MInt m_pressureLossFlameSpeed;
  MFloat m_pressureLossCorrection;
  MFloat m_meanVelocity;
  MFloat m_meanVelocityOutlet;
  MFloat m_tubeLength;
  MFloat m_outletLength;
  MFloat m_radiusOutlet;
  MFloat m_strouhal = -1.0;
  MFloat m_strouhalInit = -1.0;
  MFloat m_flameStrouhal;
  MFloat m_neutralFlameStrouhal;
  MInt m_samplingEndCycle;
  MInt m_samplingStartCycle;
  MInt m_samplingStartIteration;
  MInt m_noSamplingCycles;
  MFloat m_samplesPerCycle;
  MInt m_noForcingCycles;
  MBool m_forcing;
  MFloat m_forcingAmplitude;
  MFloat m_perturbationAmplitude;
  MFloat m_perturbationAmplitudeCorr;
  MFloat m_lambdaPerturbation;
  MInt m_outputOffset;
  MFloat m_marksteinLength;
  MFloat m_marksteinLengthPercentage;
  MBool m_zeroLineCorrection;
  MFloat m_marksteinLengthTh;
  // end levelset

 public:
  MInt noSolverTimers(const MBool allTimings) override {
#ifdef MAIA_TIMER_FUNCTION
    // 11 additional times are created, but only 6 are written
    const MInt additionalTimer = m_levelSetMb ? 6 : 0;
    if(allTimings) {
      return 2 + 17 + additionalTimer;
    } else {
      return 6;
    }
#else
    return 2;
#endif
  }
  void getSolverTimings(std::vector<std::pair<MString, MFloat>>& solverTimings, const MBool allTimings) override;
  void limitWeights(MFloat*) override;
  /// Methods to inquire solver data information
  MInt noCellDataDlb() const override { return 1; };
  MInt cellDataTypeDlb(const MInt dataId) const override {
    if(dataId != 0) {
      TERMM(1, "solverCelldataType: invalid data id");
    }
    return MFLOAT;
  };
  MInt cellDataSizeDlb(const MInt dataId, const MInt gridCellId) override;

  /// Return solver data for DLB
  void getCellDataDlb(const MInt dataId, const MInt oldNoCells, const MInt* const bufferIdToCellId,
                      MFloat* const data) override;
  /// Set solver data for DLB
  void setCellDataDlb(const MInt dataId, const MFloat* const data) override;

  /// Get/set global solver variables during DLB
  void getGlobalSolverVars(std::vector<MFloat>& globalFloatVars, std::vector<MInt>& globalIdVars) override;
  void setGlobalSolverVars(std::vector<MFloat>& globalFloatVars, std::vector<MInt>& globalIdVars) override;

  MBool hasSplitBalancing() const override { return true; }

  virtual MBool adaptationTrigger() { return false; }

  MBool forceAdaptation() override {
    if(grid().hasInactiveRanks()) {
      this->startLoadTimer(AT_);
      MPI_Allreduce(MPI_IN_PLACE, &m_forceAdaptation, 1, MPI_C_BOOL, MPI_LOR, grid().raw().mpiComm(), AT_,
                    "MPI_IN_PLACE", "m_forceAdaptation");
      this->stopLoadTimer(AT_);
    }
    return m_forceAdaptation;
  }

 protected:
  // Dynamic load balancing
  MInt m_loadBalancingReinitStage = -1;
  MBool m_weightBndryCells = true;
  MBool m_weightCutOffCells = true;
  MBool m_weightBc1601 = true;
  MBool m_weightInactiveCell = true;
  MBool m_weightNearBndryCells = false;
  MBool m_limitWeights = false;
  MBool m_weightLvlJumps = false;
  MBool m_weightSmallCells = false;

  MFloat m_weightBaseCell = 0.0;
  MFloat m_weightLeafCell = 0.05;
  MFloat m_weightActiveCell = 0.1;
  MFloat m_weightBndryCell = 1.0;
  MFloat m_weightNearBndryCell = 0.0;
  MFloat m_weightMulitSolverFactor = 1.0;

  void initializeTimers();

  // Timers
  // Timer group which holds all solver-wide timers
  MInt m_timerGroup = -1;
  // Stores all solver-wide timers
  std::array<MInt, Timers::_count> m_timers{};

  MInt m_tgfv;
  MInt m_tcomm;
  MInt m_texchange;
  MInt m_tgatherAndSend;
  MInt m_tgatherAndSendWait;
  MInt m_tscatterWaitSome;
  MInt m_tgather;
  MInt m_tsend;
  MInt m_treceive;
  MInt m_tscatter;
  MInt m_treceiving;
  MInt m_treceiveWait;
  MInt m_texchangeDt;
  MFloat m_UInfinity;
  MFloat m_VInfinity;
  MFloat m_WInfinity;
  MFloat m_PInfinity;
  MFloat m_TInfinity = NAN;
  MFloat m_DthInfinity;
  MFloat m_nuTildeInfinity;
  MFloat m_kInfinity;
  MFloat m_omegaInfinity;
  MFloat m_DInfinity;
  MFloat m_SInfinity;
  MFloat m_VVInfinity[3];
  MFloat m_rhoUInfinity;
  MFloat m_rhoVInfinity;
  MFloat m_rhoWInfinity;
  MFloat m_rhoEInfinity;
  MFloat m_rhoInfinity;
  MFloat m_rhoVVInfinity[3];
  MFloat* m_postShockPV = nullptr;
  MFloat* m_postShockCV = nullptr;

  /// \name Time-stepping
  ///@{
  /// Method to compute the time-step in a cell.
  MInt m_timeStepMethod = -1;
  /// Should the time-step in the boundary cells be weighted by their volume?
  MBool m_timeStepVolumeWeighted = false;

  /// Current time-step used to advance the solution
  ///
  /// \warning When using non-blocking time reduction this variable might not
  /// contain the correct time "yet", use `timeStep()` instead.
  MFloat m_timeStep = -1.0;
  /// Reduce the timeStep using non-blocking communication;
  MBool m_timeStepNonBlocking = false;
  /// Request for reducing the time-step using non-blocking comm.
  MPI_Request m_timeStepReq;
  /// Has the non-blocking reduction of the time-step finished?
  MBool m_timeStepAvailable = true;
  /// time-step has been updated
  MBool m_timeStepUpdated = true;
  /// Returns the current time-step.
  ///
  /// \warning When using non-blocking time reduction this function will block
  /// if the time is not available.
  MFloat timeStep(MBool canSolver = false) noexcept {
    if(!m_timeStepAvailable) {
#ifdef MAIA_FV_LOG_ACCESS_TO_UNAVAILABLE_TIME_STEP
      if(!canSolver) {
        MInt flag;
        MPI_Test(&m_timeStepReq, &flag, MPI_STATUS_IGNORE, AT_);
        if(!flag) {
          mTerm(1, AT_,
                "The time-step was required before it was available;"
                "use timeStep(true) if you want the timeStep call to solver communication");
        }
      }
#endif
      MPI_Wait(&m_timeStepReq, MPI_STATUS_IGNORE, AT_);
      m_timeStepAvailable = true;
      m_log << "Computed global time step (method: " << m_timeStepMethod << " ): " << m_timeStep
            << " - physical: " << m_timeStep * m_timeRef << std::endl;
    }
    return m_timeStep;
  }
  /// Forces the time-step to be set to a fixed value:
  MFloat m_timeStepFixedValue = -1.0;
  /// How often should the time-step be recomputed?
  MInt m_timeStepComputationInterval = -1;
  /// Returns true if the time-step from a restart file should be reused
  MBool useTimeStepFromRestartFile() const;
  /// Returns true if the time-step should be updated on this step
  MBool requiresTimeStepUpdate() const;

  // This is the time step written down to restart files. This is the same as
  // the normal time step unless you are using the level-set solver. In this case:
  // - if you are running a 2D simulation with timeStepMethod 17511 this _might_
  //   do something meaningful (it just writes a different time step than the one
  //   used).
  // - otherwise this writes -1 so don't use the time-step from the restart file
  //   when restarting
  MFloat m_restartFileOutputTimeStep = -1.;

  void setRestartFileOutputTimeStep();

  void computeSourceTerms();

 public:
  /// Force time step externally
  void forceTimeStep(const MFloat dt) { m_timeStep = dt; }

  MInt a_timeStepComputationInterval() { return m_timeStepComputationInterval; }

  void preTimeStep() override { m_timeStepUpdated = false; }

 protected:
  struct MV {
    // Mean Lamb vector
    static constexpr const MInt LAMB0 = 0;

    // Mean vorticity
    static constexpr const MInt VORT0 = 1;

    // du/dx, du/dy, dw/dz for the divergence
    static constexpr const MInt DU = 2;

    // Mean gradient of rho
    static constexpr const MInt DRHO = 3;

    // Mean gradient of p
    static constexpr const MInt DP = 4;

    // Mean gradient of rho*div(u)
    static constexpr const MInt RHODIVU = 5;

    // Mean gradient of u*grad(rho)
    static constexpr const MInt UGRADRHO = 6;

    // Mean of (gradient of p divided by rho)
    static constexpr const MInt GRADPRHO = 7;

    // Mean gradients of velocity components (contains MV::DU)
    static constexpr const MInt GRADU = 8;

    // Sum of products of velocity and velocity gradients:
    // u * grad(u) + v * grad(v) + w * grad(w)
    static constexpr const MInt UGRADU = 9;
  };

  // Data members used to replace local static variables
 private:
  MBool m_firstUseWriteVtuOutputParallelQout = true;
  MBool m_firstUseWriteVtuOutputParallelGeom = true;
  MBool m_firstUseInitializeVtkXmlOutput = true;

 protected:
  // Data members used to replace local static variables for src/fvsolverxd.h
  // logCell
  MBool m_static_logCell_firstRun = true;
  // smallCellCorrection
  MBool m_static_smallCellCorrection_first = true;
  MInt m_static_smallCellCorrection_slipDirection;
  MFloat m_static_smallCellCorrection_slipCoordinate;
  // computeSurfaceValuesLimitedSlopesMan
  MBool m_static_computeSurfaceValuesLimitedSlopesMan_checkedBndryCndIds = false;
  MBool m_static_computeSurfaceValuesLimitedSlopesMan_correctWallBndryFluxes = false;

  // Data members used to replace local static variables for src/fvmbcartesiansolverxd.h
  // getDistanceSplitSphere
  MFloat m_static_getDistanceSplitSphere_h = 0;
  MBool m_static_getDistanceSplitSphere_first = true;
  // constructGFieldPredictor
  MBool m_static_constructGFieldPredictor_firstRun = true;
  MBool m_static_constructGFieldPredictor_adaptiveGravity = false;
  // advanceSolution
  MFloat m_static_advanceSolution_meanDragCoeff = F0;
  MFloat m_static_advanceSolution_meanDrag = F0;
  MFloat m_static_advanceSolution_dragCnt = F0;
  MBool m_static_advanceSolution_firstRun = true;
  // updateBodyProperties
  MBool m_static_updateBodyProperties_c453_firstRun = true;
  MBool m_static_updateBodyProperties_c455_firstRun = true;
  MBool m_static_updateBodyProperties_firstTime = true;
  // redistributeMass
  MBool m_static_redistributeMass_firstRun = true;
  // applyBoundaryCondition
  MBool m_static_applyBoundaryCondition_firstRun = true;
  MFloat m_static_applyBoundaryCondition_ERhoL1 = F0;
  MFloat m_static_applyBoundaryCondition_ERhoL2 = F0;
  MFloat m_static_applyBoundaryCondition_ERhoLoo = F0;
  MFloat m_static_applyBoundaryCondition_EVelL1 = F0;
  MFloat m_static_applyBoundaryCondition_EVelL2 = F0;
  MFloat m_static_applyBoundaryCondition_EVelLoo = F0;
  MFloat m_static_applyBoundaryCondition_refMass;
  MFloat m_static_applyBoundaryCondition_oldMass;
  MFloat m_static_applyBoundaryCondition_oldVol2;
  // updateSpongeLayer
  MBool m_static_updateSpongeLayer_mbSpongeLayer = false;
  MBool m_static_updateSpongeLayer_first = true;
  // logData
  MBool m_static_logData_firstRun4 = true;
  // logData initial condition 45299
  MBool m_static_logData_ic45299_first = true;
  static constexpr MInt s_logData_ic45299_maxNoEmbeddedBodies = 20;
  MFloat m_static_logData_ic45299_amplitude[s_logData_ic45299_maxNoEmbeddedBodies];
  MFloat m_static_logData_ic45299_freqFactor[s_logData_ic45299_maxNoEmbeddedBodies];
  MFloat m_static_logData_ic45299_maxF, m_static_logData_ic45299_maxA;
  MFloat m_static_logData_ic45299_cutOffAngle;
  MFloat m_static_logData_ic45299_xCutOff;
  // logData initial condition 45301
  MBool m_static_logData_ic45301_first = true;
  static constexpr MInt s_logData_ic45301_maxNoEmbeddedBodies = 20;
  static constexpr MInt s_logData_ic45301_maxNoPressurePoints = 20;
  MFloat m_static_logData_ic45301_Strouhal;
  MFloat m_static_logData_ic45301_freqFactor[s_logData_ic45301_maxNoEmbeddedBodies];
  MFloat m_static_logData_ic45301_maxF;
  MFloat m_static_logData_ic45301_pressurePoints[s_logData_ic45301_maxNoPressurePoints * nDim];
  MInt m_static_logData_ic45301_noPressurePoints = 0;
  MInt m_static_logData_ic45301_containingCellIds[s_logData_ic45301_maxNoPressurePoints];
  // saveSolverSolution
  MBool m_static_saveSolverSolutionxd_firstRun = true;
  // computeFlowStatistics
  static constexpr MInt s_computeFlowStatistics_noSamples = 1; // with temporal averaging
  static constexpr MInt s_computeFlowStatistics_thetaSize = 100;
  static constexpr MInt s_computeFlowStatistics_noDat = 76; // mutiples of 4!
  MFloat m_static_computeFlowStatistics_thetaDensityAverage[s_computeFlowStatistics_noSamples
                                                            * s_computeFlowStatistics_thetaSize
                                                            * s_computeFlowStatistics_noDat];
  MFloat m_static_computeFlowStatistics_thetaDensityAverage2[s_computeFlowStatistics_noSamples
                                                             * s_computeFlowStatistics_thetaSize
                                                             * s_computeFlowStatistics_noDat];
  MInt m_static_computeFlowStatistics_currentIndex = 0;
  static constexpr MInt s_computeFlowStatistics_noPdfPoints = 100;
  static constexpr MInt s_computeFlowStatistics_noPdfs = 42;
  static constexpr MInt s_computeFlowStatistics_noJointPdfs = 6;
  static constexpr MInt s_computeFlowStatistics_noSamples2 = 1; // with temporal averaging
  MFloat m_static_computeFlowStatistics_pdfAverage[s_computeFlowStatistics_noSamples2 * s_computeFlowStatistics_noPdfs
                                                   * s_computeFlowStatistics_noPdfPoints];
  MFloat m_static_computeFlowStatistics_pdfAverage2[s_computeFlowStatistics_noSamples2 * s_computeFlowStatistics_noPdfs
                                                    * s_computeFlowStatistics_noPdfPoints];
  MFloat m_static_computeFlowStatistics_jointPdfAverage[s_computeFlowStatistics_noSamples2
                                                        * s_computeFlowStatistics_noJointPdfs
                                                        * s_computeFlowStatistics_noPdfPoints
                                                        * s_computeFlowStatistics_noPdfPoints];
  MInt m_static_computeFlowStatistics_currentIndex2 = 0;
  static constexpr MInt s_computeFlowStatistics_noSamples3 = 1; // with temporal averaging
  static constexpr MInt s_computeFlowStatistics_noAngles = 20;
  static constexpr MInt s_computeFlowStatistics_noAngleDat = 8;
  MFloat m_static_computeFlowStatistics_sdatAverage[s_computeFlowStatistics_noSamples * s_computeFlowStatistics_noAngles
                                                    * s_computeFlowStatistics_noAngleDat];
  MFloat m_static_computeFlowStatistics_sdatAverage2[s_computeFlowStatistics_noSamples
                                                     * s_computeFlowStatistics_noAngles
                                                     * s_computeFlowStatistics_noAngleDat];
  MFloat m_static_computeFlowStatistics_sdatSum[s_computeFlowStatistics_noSamples * s_computeFlowStatistics_noAngles];
  MInt m_static_computeFlowStatistics_currentIndex3 = 0;
  MInt m_maxNearestBodies = 20;
  // near-particle statistics
  static constexpr MInt s_computeFlowStatistics_noReClasses = 6;
  static constexpr MInt s_computeFlowStatistics_noSamples4 = 1; // with temporal averaging
  MInt m_static_computeFlowStatistics_currentIndex4 = 0;
  MInt m_static_computeFlowStatistics_currentCnt = 0;
  MFloat m_static_computeFlowStatistics_bodyCntAvg[s_computeFlowStatistics_noSamples
                                                   * s_computeFlowStatistics_noReClasses];
  MBool m_static_computeFlowStatistics_firstBD = true;
  // writeVtkXmlFiles
  MBool m_static_writeVtkXmlFiles_firstCall = true;
  MBool m_static_writeVtkXmlFiles_firstCall2 = true;
  // logCellxd
  MBool m_static_logCellxd_firstRun = true;


 public:
  friend class LsFvCombustion<nDim_, SysEqn>;
  friend class CouplerFvMultilevel<nDim_, SysEqn>;
  friend class LsFvMb<nDim_, SysEqn>;
  // has to be explicit since only RANS solver is template parameter
  friend class FvZonal<nDim_, FvSysEqnRANS<nDim, RANSModelConstants<RANS_SA_DV>>>;
  friend class FvZonal<nDim_, FvSysEqnRANS<nDim, RANSModelConstants<RANS_FS>>>;
  friend class FvZonalRTV<nDim_, FvSysEqnRANS<nDim, RANSModelConstants<RANS_SA_DV>>>;
  friend class FvZonalRTV<nDim_, FvSysEqnRANS<nDim, RANSModelConstants<RANS_FS>>>;
  friend class FvZonalSTG<nDim_, FvSysEqnRANS<nDim, RANSModelConstants<RANS_SA_DV>>>;
  friend class FvZonalSTG<nDim_, FvSysEqnRANS<nDim, RANSModelConstants<RANS_FS>>>;
  friend class FvCartesianInterpolation<nDim_, FvSysEqnRANS<nDim, RANSModelConstants<RANS_SA_DV>>,
                                        FvSysEqnRANS<nDim, RANSModelConstants<RANS_SA_DV>>>;
  friend class FvCartesianInterpolation<nDim_, FvSysEqnRANS<nDim, RANSModelConstants<RANS_SA_DV>>, FvSysEqnNS<nDim>>;
  friend class VtkIo<nDim_, SysEqn>;

  void computeVolumeForcesRANS();

  void exchangeGapInfo();
  void resetCutOffCells();
  void resetSponge();
  void exchangeProperties();

  void computeUTau(MFloat* data, MInt cellId); // is never called

  // sponge
  void computeDomainAndSpongeDimensions();
  virtual void updateSpongeLayer();
  std::array<MFloat, 6> computeTargetValues();
  std::array<MFloat, nDim_ + 2> computeSpongeDeltas(MInt cellId, std::array<MFloat, 6>);
  // MFloat meanDensity, MFloat meanPressure, MFloat
  // meanDensityIn, MFloat meanPressureIn);
  void updateSpongeLayerRhs(MInt, std::array<MFloat, nDim_ + 2>);
  void checkCells(); //<-- XD version, but by now only used in 3D
  void checkForSrfcsMGC();
  void checkForSrfcsMGC_2();
  void checkForSrfcsMGC_2_(); //<-- Used only in 3D
  void computeSrfcs(MInt, MInt, MInt, MInt);
  virtual void correctBoundarySurfaces();
  virtual void correctBoundarySurfaces_(); //<-- Only required in 3D
  void checkCellSurfaces();
  void computeCellSurfaceDistanceVectors();

  void computeReconstructionConstantsSVD();
  void setConservativeVarsOnAzimuthalRecCells();
  void initAzimuthalReconstruction();
  virtual void computeAzimuthalReconstructionConstants(MInt mode = 0);
  void rebuildAzimuthalReconstructionConstants(MInt cellId, MInt offset, MFloat* recCoord, MInt mode = 0);
  void interpolateAzimuthalData(MFloat* data, MInt offset, MInt noVars, const MFloat* vars);
  void fillExcBufferAzimuthal(MInt cellId, MInt offset, MFloat* dataDest, MFloat* dataSrc, MInt noData,
                              const std::vector<MInt>& rotIndex = std::vector<MInt>());
  void rotateVectorAzimuthal(MInt side, MFloat* data, MInt noData, const std::vector<MInt>& indices);
  void interpolateAzimuthalDataReverse(MFloat* data, MInt offset, MInt noVars, const MFloat* vars);
  void checkAzimuthalRecNghbrConsistency(MInt cellId);
  void buildLeastSquaresStencilSimple();
  void initViscousFluxComputation();
  // void copyVarsToSmallCells() override;

  // azimuthal periodicity concept
  void computeRecConstPeriodic();
  void computeRecConstPeriodic_(); //<-- Used only in 3D
  void identPeriodicCells();
  void identPeriodicCells_(); //<-- Used only in 3D

  // void Muscl(MInt timerId = -1) override;
  // void nonReflectingBCCutOff() override;
  // void nonReflectingBCAfterTreatmentCutOff() override;

  // ATTRIBUTES2(ATTRIBUTE_HOT, ATTRIBUTE_FLATTEN) virtual void LSReconstructCellCenter();
  ATTRIBUTES2(ATTRIBUTE_ALWAYS_INLINE, ATTRIBUTE_HOT)
  inline void LSReconstructCellCenter_(const MUint noSpecies);
  void LSReconstructCellCenterCons(const MUint noSpecies); // ATTRIBUTES2(ATTRIBUTE_ALWAYS_INLINE,ATTRIBUTE_HOT);
  ATTRIBUTES2(ATTRIBUTE_HOT, ATTRIBUTE_FLATTEN) virtual void computeSurfaceValues(MInt timerId = -1);
  ATTRIBUTES2(ATTRIBUTE_ALWAYS_INLINE, ATTRIBUTE_HOT) inline void computeSurfaceValues_(const MUint);
  ATTRIBUTES2(ATTRIBUTE_HOT, ATTRIBUTE_FLATTEN) virtual void computeSurfaceValuesLimited(MInt timerId = -1);
  virtual void computeSurfaceValuesLOCD(MInt timerId = -1);
  virtual void computeLimitedSurfaceValues(MInt timerId = -1);
  virtual void computeSurfaceValuesLimitedSlopes(MInt timerId = -1);
  virtual void computeSurfaceValuesLimitedSlopesMan(MInt timerId = -1);
  virtual void initComputeSurfaceValuesLimitedSlopesMan1();
  virtual void initComputeSurfaceValuesLimitedSlopesMan2();

  void computeGridCellCoordinates(MFloat*);
  void findNghbrIdsMGC();
  void findDirectNghbrs(MInt cellId, std::vector<MInt>& nghbrList);
  void findNeighborHood(MInt cellId, MInt layer, std::vector<MInt>& nghbrList);

  void refineCell(const MInt) override;
  void removeChilds(const MInt) override;
  void removeCell(const MInt) override;
  void swapCells(const MInt, const MInt) override;
  void resizeGridMap() override;
  void swapProxy(const MInt, const MInt) override;
  MBool cellOutside(const MInt);
  MInt cellOutside(const MFloat*, const MInt, const MInt) override;
  // virtual void getBoundaryDistance(MFloatScratchSpace&);
  virtual void resetSurfaces();
  virtual void resetBoundaryCells(const MInt offset = 0);

  // Reallocate memory of arrays depending on current number of cells
  void reInitActiveCellIdsMemory();
  // Reallocate memory of master/small cell id arrays depending on current number of cells
  void reInitSmallCellIdsMemory();

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
  void sensorInterfaceDelta(std::vector<std::vector<MFloat>>& sensors, std::vector<std::bitset<64>>& sensorCellFlag,
                            std::vector<MFloat>& sensorWeight, MInt sensorOffset, MInt sen);
  void sensorInterfaceLsMb(std::vector<std::vector<MFloat>>& sensors, std::vector<std::bitset<64>>& sensorCellFlag,
                           std::vector<MFloat>& sensorWeight, MInt sensorOffset, MInt sen);
  void sensorInterfaceLs(std::vector<std::vector<MFloat>>& sensors, std::vector<std::bitset<64>>& sensorCellFlag,
                         std::vector<MFloat>& sensorWeight, MInt sensorOffset, MInt sen);
  void sensorCutOff(std::vector<std::vector<MFloat>>& sensors, std::vector<std::bitset<64>>& sensorCellFlag,
                    std::vector<MFloat>& sensorWeight, MInt sensorOffset, MInt sen);
  void sensorPatch(std::vector<std::vector<MFloat>>& sensors, std::vector<std::bitset<64>>& sensorCellFlag,
                   std::vector<MFloat>& sensorWeight, MInt sensorOffset, MInt sen) override;
  void bandRefinementSensorDerivative(std::vector<std::vector<MFloat>>& sensors,
                                      std::vector<std::bitset<64>>& sensorCellFlag, MInt sensorOffset, MInt sen,
                                      const std::vector<MFloat>& tau, const MFloat sensorThreshold);

  MBool m_sensorBandRefinement;
  MInt m_sensorBandRefinementAdditionalLayers;

  virtual void updateMultiSolverInformation(MBool fullReset = false);
  void resetSolver() override;
  virtual void resetSolverFull();
  void setCellWeights(MFloat*) override;
  void balance(const MInt* const noCellsToReceiveByDomain,
               const MInt* const noCellsToSendByDomain,
               const MInt* const targetDomainsByCell,
               const MInt oldNoCells) override;
  void balancePre() override;
  void balancePost() override;
  void finalizeBalance() override;
  // Partitioning
  MInt noLoadTypes() const override;
  void getDefaultWeights(MFloat* weights, std::vector<MString>& names) const;
  void getLoadQuantities(MInt* const loadQuantities) const override;
  MFloat getCellLoad(const MInt cellId, const MFloat* const weights) const override;

  void getDomainDecompositionInformation(std::vector<std::pair<MString, MInt>>& domainInfo) override;

  void saveSolverSolution(const MBool forceOutput = false, const MBool finalTimeStep = false) override;
  void writeRestartFile(const MBool, const MBool, const MString, MInt*) override;
  void writeRestartFile(MBool) override{};
  MBool prepareRestart(MBool, MBool&) override;
  // void reIntAfterRestart(MBool) override;
  void saveSampleFiles();
  virtual void saveRestartFile(const MBool);
  void saveDebugRestartFile();
  void loadRestartFile() override;
  void computeConservativeVariablesCoarseGrid();
  // void computePrimitiveVariablesCoarseGrid() override;
  void computePrimitiveVariablesCoarseGrid(MInt);

  void writeCenterLineVel();
  void computeForceCoefficients(MFloat*);

  void checkInfinityVarsConsistency();
  void computeInitialPressureLossForChannelFlow();

  template <MInt stencil, class F>
  inline void Ausm_(F& fluxFct);

  ATTRIBUTES2(ATTRIBUTE_HOT, ATTRIBUTE_ALWAYS_INLINE)
  inline void computePrimitiveVariables_();
  template <class _ = void, std::enable_if_t<isEEGas<SysEqn>, _*> = nullptr>
  ATTRIBUTES2(ATTRIBUTE_HOT, ATTRIBUTE_ALWAYS_INLINE)
  inline MBool uDLimiter(const MFloat* const, MFloat* const);
  ATTRIBUTES2(ATTRIBUTE_HOT, ATTRIBUTE_ALWAYS_INLINE) inline void computeConservativeVariables_();
  template <class _ = void, std::enable_if_t<isEEGas<SysEqn>, _*> = nullptr>
  void bubblePathDispersion();
  template <MInt stencil, class F>
  void viscousFlux_(F& viscousFluxFct);
  template <MInt stencil, class F>
  void viscousFluxCompact_(F& viscousFluxFct);
  template <MInt noSpecies>
  void viscousFluxMultiSpecies_();

  void computeConservativeVariablesMultiSpecies_(const MUint);
  void computePrimitiveVariablesMultiSpecies_(const MUint);

  virtual void distributeFluxToCells();
  void implicitTimeStep() override;
  virtual MBool maxResidual(MInt mode = 0);
  void computeCellVolumes();
  template <class _ = void, std::enable_if_t<!isEEGas<SysEqn>, _*> = nullptr>
  ATTRIBUTES2(ATTRIBUTE_HOT, ATTRIBUTE_ALWAYS_INLINE)
  inline MBool rungeKuttaStep_(const MUint);
  template <class _ = void, std::enable_if_t<isEEGas<SysEqn>, _*> = nullptr>
  ATTRIBUTES2(ATTRIBUTE_HOT, ATTRIBUTE_ALWAYS_INLINE)
  inline MBool rungeKuttaStepEEGas();
  void computeSamplingTimeStep();
  void computeSamplingTimeStep_(); //<-- Only available in 2D

  void computeCoarseGridCorrection(MInt);
  void linearInterpolation(MInt, MInt, MInt*);
  void bilinearInterpolation(MInt, MInt, MInt*);
  void bilinearInterpolationAtBnd(MInt, MInt, MInt*);

  void initCutOffBoundaryCondition();
  virtual void setConservativeVariables(MInt cellId);
  void setPrimitiveVariables(MInt cellId);
  template <class _ = void, std::enable_if_t<isEEGas<SysEqn>, _*> = nullptr>
  void initDepthCorrection();
  void divCheck(MInt);

  template <MInt diag, MBool recorrectBndryCellCoords = false>
  MInt getAdjacentLeafCells(const MInt, const MInt, MIntScratchSpace&, MIntScratchSpace&);
  template <MBool recorrectBndryCellCoords = false>
  MInt getNghbrLeafCells(const MInt cellId, MInt refCell, MInt layer, MInt* nghbrs, MInt dir, MInt dir1 = -1,
                         MInt dir2 = -1) const;

  void reduceVariables();
  void computeSlopesByCentralDifferences();

  void generateBndryCells();
  void createBoundaryCells();

  // function pointers
  void (FvCartesianSolverXD::*m_reconstructSurfaceData)(MInt);
  void (FvCartesianSolverXD::*m_computeViscousFlux)();

  void (FvCartesianSolverXD::*m_computeViscousFluxMultiSpecies)(MInt);

  void getInterpolatedVariables(const MInt cellId, const MFloat* position, MFloat* vars) override;
  template <MInt noSpecies>
  void getInterpolatedVariables(const MInt cellId, const MFloat* position, MFloat* vars);
  void getInterpolatedVariablesInCell(const MInt cellId, const MFloat* position, MFloat* vars);
  template <MInt a, MInt b, MBool old = false>
  void interpolateVariables(const MInt cellId, const MFloat* position, MFloat* result);
  template <MInt a, MInt b>
  void interpolateVariablesInCell(const MInt cellId, const MFloat* position,
                                  std::function<MFloat(MInt, MInt)> variables, MFloat* result);

  //  MString get_mSysEqnName();
  MFloat crankAngle(const MFloat, const MInt);

  /// Returns the pressure computed from the conservative variables of the cell \p cellId
  MFloat cv_p(MInt cellId) const noexcept;
  /// Returns the temperature computed from the conservative variables of the cell \p cellId
  MFloat cv_T(MInt cellId) const noexcept;
  /// Returns the speed-of-sound computed from the conservative variables of the cell \p cellId
  MFloat cv_a(MInt cellId) const noexcept;

  void setTimeStep();

 protected:
  MFloat computeTimeStep(MInt cellId) const noexcept;
  MFloat computeTimeStepMethod(MInt cellId) const noexcept;
  MFloat computeTimeStepEulerDirectional(MInt cellId) const noexcept;
  MFloat computeTimeStepApeDirectional(MInt cellId) const noexcept;
  MBool cellParticipatesInTimeStep(MInt cellId) const noexcept;
  MFloat computeTimeStepDiffusionNS(MFloat density, MFloat temperature, MFloat Re, MFloat C, MFloat dx) const noexcept;
  void computeAndSetTimeStep();

  // Storage for required interpolation information in certain cells
  std::vector<MInt> m_cellInterpolationIndex{};
  std::vector<std::vector<MFloat>> m_cellInterpolationMatrix{};
  std::vector<std::vector<MInt>> m_cellInterpolationIds{};

 public:
  MFloat m_meanCoord[3]{};
  MInt m_adaptationLevel;
  MBool m_wasAdapted = false;

  MBool m_wasBalancedZonal = false;

 public:
  // Data members used to replace local static variables
  MBool m_firstUseUpdateSpongeLayerCase51 = true;
  MPI_Comm comm_sponge{};
  MInt m_noSpongeZonesIn{};
  MInt m_noSpongeZonesOut{};
  static constexpr MInt s_maxNoSpongeZones = 10;
  MInt m_spongeDirectionsIn[s_maxNoSpongeZones]{};
  MInt m_spongeDirectionsOut[s_maxNoSpongeZones]{};
  MInt m_secondSpongeDirectionsIn[s_maxNoSpongeZones]{};
  MInt m_secondSpongeDirectionsOut[s_maxNoSpongeZones]{};
  MInt m_spongeAveragingIn[s_maxNoSpongeZones]{};
  MInt m_spongeAveragingOut[s_maxNoSpongeZones]{};
  MBool m_hasCellsInSpongeLayer{};
  MFloat m_timeOfMaxPdiff{};
  MFloat *m_coordSpongeIn{}, *m_secondCoordSpongeIn{}, *m_coordSpongeOut{}, *m_secondCoordSpongeOut{}, *m_uNormal_r{};

  // Level-set solver ...
  static constexpr MInt s_maxNoEmbeddedBodies = 20;
  std::map<MInt, std::vector<MInt>> m_cellToNghbrHood;

  MFloat m_maxLsValue = -99;
  MBool m_linerLvlJump = false;

 private:
  // crankAngle
  MFloat m_static_crankAngle_Strouhal = -99;
  MFloat m_static_crankAngle_initialCad = -99;

  MInt* m_particleWidth = nullptr;

 public:
  // Sandpaper Tripping stuff
  void initSandpaperTrip();
  void applySandpaperTrip();
  void saveSandpaperTripVars();
  void tripForceCoefficients(MFloat*, MFloat*, MFloat*, MInt, MInt);
  void tripFourierCoefficients(MFloat*, MInt, MFloat, MFloat);

  void dumpCellData(const MString name);

  MInt m_tripNoTrips;
  MInt m_tripNoModes;
  MInt m_tripTotalNoCells;
  MInt m_tripSeed;
  // MBool m_useSandpaperTrip;
  MBool m_tripUseRestart;
  MFloat m_tripDomainWidth;

  MBool m_tripAirfoil;
  MFloat m_tripAirfoilChordLength;
  MFloat m_tripAirfoilAOA;
  MFloat* m_tripAirfoilNosePos = nullptr;
  MFloat* m_tripAirfoilChordPos = nullptr;
  MFloat* m_tripAirfoilForceDir = nullptr;
  MInt* m_tripAirfoilOrientation = nullptr;
  MInt* m_tripAirfoilBndryId = nullptr;
  MInt* m_tripAirfoilSide = nullptr;

  std::vector<MInt> m_tripCellIds;
  std::vector<MFloat> m_tripCoords;

  MInt* m_tripTimeStep = nullptr;
  MInt* m_tripNoCells = nullptr;
  MInt* m_tripCellOffset = nullptr;
  MFloat* m_tripDelta1 = nullptr;
  MFloat* m_tripXOrigin = nullptr;
  MFloat* m_tripXLength = nullptr;
  MFloat* m_tripYOrigin = nullptr;
  MFloat* m_tripYHeight = nullptr;
  MFloat* m_tripCutoffZ = nullptr;
  MFloat* m_tripMaxAmpSteady = nullptr;
  MFloat* m_tripMaxAmpFluc = nullptr;
  MFloat* m_tripDeltaTime = nullptr;
  MFloat* m_tripG = nullptr;
  MFloat* m_tripH1 = nullptr;
  MFloat* m_tripH2 = nullptr;
  MFloat* m_tripModesG = nullptr;
  MFloat* m_tripModesH1 = nullptr;
  MFloat* m_tripModesH2 = nullptr;

  MBool m_useMpiStartall = true;
  MBool m_onlyMaxLvlMpiRequests = true;
  MInt m_noMaxLvlMpiSendNeighbors = -1;
  MInt m_noMaxLvlMpiRecvNeighbors = -1;
  std::vector<MInt> m_maxLvlMpiSendNeighbor{};
  std::vector<MInt> m_maxLvlMpiRecvNeighbor{};
};


/// \brief Store the solver data for a given data id ordered in the given buffer for DLB
template <MInt nDim_, class SysEqn_>
// template <typename dataType>
void FvCartesianSolverXD<nDim_, SysEqn_>::getCellDataDlb(const MInt dataId, const MInt oldNoCells,
                                                         const MInt* const bufferIdToCellId, MFloat* const data) {
  TRACE();

  MInt localBufferId = 0;
  for(MInt i = 0; i < oldNoCells; i++) {
    const MInt gridCellId = bufferIdToCellId[i];

    if(gridCellId < 0) continue;

    const MInt cellId = grid().tree().grid2solver(gridCellId);
    if(cellId < 0 || cellId >= noInternalCells()) {
      continue;
    }

    switch(dataId) {
      case 0: {
        const MInt dataSize = CV->noVariables;
        std::copy_n(&a_variable(cellId, 0), dataSize, &data[localBufferId * dataSize]);
        break;
      }
      default:
        TERMM(1, "Unknown data id.");
        break;
    }
    localBufferId++;
  }
}


/// \brief Set the solver cell data after DLB
template <MInt nDim_, class SysEqn_>
// template <typename dataType>
void FvCartesianSolverXD<nDim_, SysEqn_>::setCellDataDlb(const MInt dataId, const MFloat* const data) {
  TRACE();

  // Nothing to do if solver is not active
  if(!isActive()) {
    return;
  }

  // Set the variables if this is the correct reinitialization stage
  if(m_loadBalancingReinitStage == 0) {
    switch(dataId) {
      case 0: {
        std::copy_n(data, noInternalCells() * CV->noVariables, &a_variable(0, 0));
        break;
      }
      default:
        TERMM(1, "Unknown data id.");
    }
  }
}


template <MInt nDim, class SysEqn>
template <class _, std::enable_if_t<isEEGas<SysEqn>, _*>>
inline MBool FvCartesianSolverXD<nDim, SysEqn>::uDLimiter(const MFloat* const uOtherPhase, MFloat* const pvars) {
  const MFloat uDLim = m_EEGas.uDLim;
  MBool change = false;
  for(MInt i = 0; i < nDim; i++) {
    if(!((pvars[i] - uOtherPhase[i]) <= uDLim)) {
      pvars[i] = uOtherPhase[i] + uDLim;
      change = true;
    } else if(!((pvars[i] - uOtherPhase[i]) >= -uDLim)) {
      pvars[i] = uOtherPhase[i] - uDLim;
      change = true;
    }
  }
  return change;
}

/// \brief Exchange data
template <MInt nDim, class SysEqn>
template <typename T>
void FvCartesianSolverXD<nDim, SysEqn>::exchangeDataFV(T* data, const MInt blockSize, MBool cartesian,
                                                       const std::vector<MInt>& rotIndex) {
  TRACE();

  RECORD_TIMER_START(m_tcomm);
  RECORD_TIMER_START(m_texchange);

  exchangeData(data, blockSize);

  if(grid().azimuthalPeriodicity()) {
    if(cartesian) {
      this->exchangeAzimuthalPer(data, blockSize);
    } else {
      if constexpr(std::is_same<T, MFloat>::value) {
        exchangeFloatDataAzimuthal(data, blockSize, rotIndex);
      } else {
        exchangeDataAzimuthal(data, blockSize);
      }
    }
  }

  RECORD_TIMER_STOP(m_texchange);
  RECORD_TIMER_STOP(m_tcomm);
}


/** brief: Exchanges ints between azimuthal window/halo cells (no interpolation)
 * Thomas Hoesgen
 */
template <MInt nDim, class SysEqn>
template <typename T>
void FvCartesianSolverXD<nDim, SysEqn>::exchangeDataAzimuthal(T* data, MInt dataBlockSize) {
  TRACE();

  if(grid().noAzimuthalNeighborDomains() == 0 && m_azimuthalRemappedNeighborDomains.size() > 0) {
    return;
  }

  this->exchangeAzimuthalPer(&data[0], dataBlockSize);

  MUint sndSize = maia::mpi::getBufferSize(m_azimuthalRemappedWindowCells);
  ScratchSpace<T> windowData(sndSize * dataBlockSize, AT_, "windowData");
  windowData.fill(0);
  MUint rcvSize = maia::mpi::getBufferSize(m_azimuthalRemappedHaloCells);
  ScratchSpace<T> haloData(rcvSize * dataBlockSize, AT_, "haloData");
  haloData.fill(0);

  // Remapped Halos
  windowData.fill(0);
  haloData.fill(0);
  MInt sndCnt = 0;
  for(MUint i = 0; i < m_azimuthalRemappedNeighborDomains.size(); i++) {
    for(MUint j = 0; j < m_azimuthalRemappedWindowCells[i].size(); j++) {
      MInt cellId = m_azimuthalRemappedWindowCells[i][j];
      for(MInt b = 0; b < dataBlockSize; b++) {
        windowData[sndCnt * dataBlockSize + b] = data[cellId * dataBlockSize + b];
      }
    }
  }
  if(m_azimuthalRemappedNeighborDomains.size() > 0) {
    maia::mpi::exchangeBuffer(m_azimuthalRemappedNeighborDomains, m_azimuthalRemappedHaloCells,
                              m_azimuthalRemappedWindowCells, mpiComm(), windowData.getPointer(), haloData.getPointer(),
                              dataBlockSize);
  }

  MInt rcvCnt = 0;
  for(MUint i = 0; i < m_azimuthalRemappedNeighborDomains.size(); i++) {
    for(MUint j = 0; j < m_azimuthalRemappedHaloCells[i].size(); j++) {
      MInt cellId = m_azimuthalRemappedHaloCells[i][j];
      for(MInt b = 0; b < dataBlockSize; b++) {
        data[cellId * dataBlockSize + b] = haloData[rcvCnt * dataBlockSize + b];
      }
    }
  }
}

// Undefine macros that should not be used outside this file
#undef ENSURE_VALID_CELL_ID_CONTAINER

#endif // FVSOLVERXD_H
