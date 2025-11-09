// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef LBSOLVER_H
#define LBSOLVER_H

#include "COUPLER/surfacecoupling.h"
#include "GEOM/geometryintersection.h"
#include "GRID/cartesiangrid.h"
#include "GRID/cartesiangridpoint.h"
#include "GRID/cartesianpointbasedcell.h"
#include "POST/postprocessing.h"
#include "UTIL/functions.h"
#include "UTIL/maiafftw.h"
#include "cartesiansolver.h"
#include "lbbndcnd.h"
#include "lbcellcollector.h"
#include "variables.h"

class LbInterfaceCell;
class LbParentCell;

template <MInt nDim, class ppType>
class PostProcessing;

// TODO labels:LB miro: This is only introduced temporarily. For now, macroscopic
// vairbales are calculated within the collision step (INCOLLISION). Some require a
// calculation before the collision (PRECOLLISION) and other want an update
// right after the propagation (POSTPROPAGATION). Until all is unified, it is
// indicated by this enum which mode is used.
enum LbUpdateMacroscopic { INCOLLISION, PRECOLLISION, POSTPROPAGATION };

/** \brief This class represents a lattice Boltzmann solver
 *
 */
template <MInt nDim>
class LbSolver : public maia::CartesianSolver<nDim, LbSolver<nDim>> {
  template <MInt nDim_>
  friend class LbBndCnd;

  template <MInt nDim_>
  friend class LbInterface;

  friend class maia::CartesianSolver<nDim, LbSolver<nDim>>;

  template <MInt nDim_, MInt nDist, class SysEqn>
  friend class CouplingLB;

  template <MInt nDim_, MInt nDist, class SysEqn>
  friend class LbLpt;

  template <MInt nDim_, class ppType>
  friend class PostProcessing;
  template <MInt nDim_>
  friend class PostProcessingLb;

 public:
  // Type for cell properties
  using Base = LbSolver<nDim>;
  using Cell = typename maia::grid::tree::Tree<nDim>::Cell;
  using SolverCell = LbCell;
  using CartesianSolver = typename maia::CartesianSolver<nDim, LbSolver>;
  using Grid = typename CartesianSolver::Grid;
  using GridProxy = typename CartesianSolver::GridProxy;

  // used CartesianSolver
  using CartesianSolver::c_globalGridId;
  using CartesianSolver::domainId;
  using CartesianSolver::domainOffset;
  using CartesianSolver::exchangeData;
  using CartesianSolver::getIdentifier;
  using CartesianSolver::grid;
  using CartesianSolver::haloCell;
  using CartesianSolver::haloCellId;
  using CartesianSolver::isActive;
  using CartesianSolver::m_bandWidth;
  using CartesianSolver::m_freeIndices;
  using CartesianSolver::m_initFromRestartFile;
  using CartesianSolver::m_innerBandWidth;
  using CartesianSolver::m_Ma;
  using CartesianSolver::m_outerBandWidth;
  using CartesianSolver::m_Re;
  using CartesianSolver::m_residualInterval;
  using CartesianSolver::m_restartFile;
  using CartesianSolver::m_restartInterval;
  using CartesianSolver::m_restartTimeStep;
  using CartesianSolver::m_solutionInterval;
  using CartesianSolver::m_solverId;
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
  using CartesianSolver::reductionFactor;
  using CartesianSolver::restartDir;
  using CartesianSolver::solverId;
  using CartesianSolver::solverMethod;
  using CartesianSolver::updateDomainInfo;
  using CartesianSolver::windowCell;
  using CartesianSolver::windowCellId;


  MInt noInternalCells() const override { return grid().noInternalCells(); }

  MBool isCompressible() { return m_isCompressible; }

  /// \brief Returns isBndryCell of the cell \p CellId
  MBool a_isBndryCell(const MInt cellId) const { return a_hasProperty(cellId, SolverCell::IsBndryCell); }

  /// \brief Returns isBndryCell of the cell \p CellId
  maia::lb::cell::BitsetType::reference a_isBndryCell(const MInt cellId) {
    return a_hasProperty(cellId, SolverCell::IsBndryCell);
  }

  /// \brief Returns isBndryGhostCell of the cell \p CellId
  MBool a_isBndryGhostCell(const MInt cellId) const { return a_hasProperty(cellId, SolverCell::IsGhost); }

  /// \brief Returns isInterface of the cell \p CellId
  MBool a_isInterface(const MInt cellId) const { return m_cells.hasProperty(cellId, SolverCell::IsInterface); }

  /// \brief Returns isInterface of the cell \p CellId
  maia::lb::cell::BitsetType::reference a_isInterface(const MInt cellId) {
    return m_cells.hasProperty(cellId, SolverCell::IsInterface);
  }

  /// \brief Returns the onlyBoundary of the cell \p cellId
  maia::lb::cell::BitsetType::reference a_onlyBoundary(const MInt cellId) {
    return a_hasProperty(cellId, SolverCell::OnlyBoundary);
  }

  /// \brief Returns the onlyBoundary of the cell \p cellId
  MBool a_onlyBoundary(const MInt cellId) const { return a_hasProperty(cellId, SolverCell::OnlyBoundary); }

  /// \brief Returns the kappa of the cell \p cellId
  MFloat& a_kappa(const MInt cellId) { return m_cells.kappa(cellId); }

  /// \brief Returns the kappa of the cell \p cellId
  MFloat a_kappa(const MInt cellId) const { return m_cells.kappa(cellId); }

  /// \brief Returns the diffusivity of the cell \p cellId
  MFloat& a_diffusivity(const MInt cellId) { return m_cells.diffusivity(cellId); }

  /// \brief Returns the diffusivity of the cell \p cellId
  MFloat a_diffusivity(const MInt cellId) const { return m_cells.diffusivity(cellId); }

  /// \brief Returns nu of the cell \p cellId
  MFloat& a_nu(const MInt cellId) { return m_cells.nu(cellId); }

  /// \brief Returns old nu of the cell \p cellId
  MFloat& a_oldNu(const MInt cellId) { return m_cells.oldNu(cellId); }

  /// \brief Returns the coordinate of the cell \p cellId for direction \p dir
  const MFloat& a_coordinate(const MInt cellId, const MInt dir) const { return grid().tree().coordinate(cellId, dir); }

  /// \brief Returns the coordinate of the cell from the grid().tree() \p cellId for dimension \p dir
  MFloat c_coordinate(const MInt cellId, const MInt dir) const { return grid().tree().coordinate(cellId, dir); }

  /// \brief Returns the distribution of the cell \p cellId for direction \p dir
  MFloat& a_distribution(const MInt cellId, const MInt dir) { return m_cells.distributions(cellId, dir); }

  /// \brief Returns the distribution of the cell \p cellId for direction \p dir
  MFloat a_distribution(const MInt cellId, const MInt dir) const { return m_cells.distributions(cellId, dir); }

  /// \brief Returns the old distribution of the cell \p cellId for direction \p dir
  MFloat& a_oldDistribution(const MInt cellId, const MInt dir) { return m_cells.oldDistributions(cellId, dir); }

  /// \brief Returns the old distribution of the cell \p cellId for direction \p dir
  MFloat a_oldDistribution(const MInt cellId, const MInt dir) const { return m_cells.oldDistributions(cellId, dir); }

  /// \brief Returns the externalForce in cell \p cellId for element \dim of vector
  MFloat& a_externalForces(const MInt cellId, const MInt dim) { return m_cells.externalForces(cellId, dim); }

  /// \brief Returns the externalForce in cell \p cellId for element \dim of vector
  MFloat a_externalForces(const MInt cellId, const MInt dim) const { return m_cells.externalForces(cellId, dim); }

  MFloat& a_previousDistribution(const MInt cellId, const MInt distr) {
    return m_cells.previousDistribution(cellId, distr);
  }

  MFloat a_previousDistribution(const MInt cellId, const MInt distr) const {
    return m_cells.previousDistribution(cellId, distr);
  }

  MFloat& a_previousVariable(const MInt cellId, const MInt varId) { return m_cells.previousVariable(cellId, varId); }

  MFloat a_previousVariable(const MInt cellId, const MInt varId) const {
    return m_cells.previousVariable(cellId, varId);
  }

  /// \brief Returns nuT of the cell \p cellId
  MFloat& a_nuT(const MInt cellId) { return m_cells.nuT(cellId); }

  /// \brief Returns nuT of the cell \p cellId
  MFloat a_nuT(const MInt cellId) const { return m_cells.nuT(cellId); }

  /// \brief Returns oldNuT of the cell \p cellId
  MFloat& a_oldNuT(const MInt cellId) { return m_cells.oldNuT(cellId); }

  /// \brief Returns oldNuT of the cell \p cellId
  MFloat a_oldNuT(const MInt cellId) const { return m_cells.oldNuT(cellId); }

  /// \brief Returns nu of the cell \p cellId
  MFloat a_nu(const MInt cellId) const { return m_cells.nu(cellId); }

  /// \brief Returns oldNu of the cell \p cellId
  MFloat a_oldNu(const MInt cellId) const { return m_cells.oldNu(cellId); }

  /// \brief Returns isInterfaceChild of the cell \p cellId
  maia::lb::cell::BitsetType::reference a_isInterfaceChild(const MInt cellId) {
    return a_hasProperty(cellId, SolverCell::IsInterfaceChild);
  }

  /// \brief Returns isInterfaceChild of the cell \p cellId
  MBool a_isInterfaceChild(const MInt cellId) const { return a_hasProperty(cellId, SolverCell::IsInterfaceChild); }

  /// \brief Returns the isInterfaceParent of the cell \p cellId
  maia::lb::cell::BitsetType::reference a_isInterfaceParent(const MInt cellId) {
    return a_hasProperty(cellId, SolverCell::IsInterfaceParent);
  }

  /// \brief Returns the isInterfaceParent of the cell \p cellId
  MBool a_isInterfaceParent(const MInt cellId) const { return a_hasProperty(cellId, SolverCell::IsInterfaceParent); }

  /// \brief Returns the distributionThermal of the cell \p cellId for direction \p dir
  MFloat& a_distributionThermal(const MInt cellId, const MInt dir) { return m_cells.distributionsThermal(cellId, dir); }

  /// \brief Returns the distributionThermal of the cell \p cellId for direction \p dir
  MFloat a_distributionThermal(const MInt cellId, const MInt dir) const {
    return m_cells.distributionsThermal(cellId, dir);
  }

  /// \brief Returns the distributionThermal of the cell \p cellId for direction \p dir
  MFloat& a_oldDistributionThermal(const MInt cellId, const MInt dir) {
    return m_cells.oldDistributionsThermal(cellId, dir);
  }

  /// \brief Returns the distributionThermal of the cell \p cellId for direction \p dir
  MFloat a_oldDistributionThermal(const MInt cellId, const MInt dir) const {
    return m_cells.oldDistributionsThermal(cellId, dir);
  }

  /// \brief Returns the distributionTransport of the cell \p cellId for direction \p dir
  MFloat& a_distributionTransport(const MInt cellId, const MInt dir) {
    return m_cells.distributionsTransport(cellId, dir);
  }

  /// \brief Returns the distributionTransport of the cell \p cellId for direction \p dir
  MFloat a_distributionTransport(const MInt cellId, const MInt dir) const {
    return m_cells.distributionsTransport(cellId, dir);
  }

  /// \brief Returns the distributionTransport of the cell \p cellId for direction \p dir
  MFloat& a_oldDistributionTransport(const MInt cellId, const MInt dir) {
    return m_cells.oldDistributionsTransport(cellId, dir);
  }

  /// \brief Returns the distributionTransport of the cell \p cellId for direction \p dir
  MFloat a_oldDistributionTransport(const MInt cellId, const MInt dir) const {
    return m_cells.oldDistributionsTransport(cellId, dir);
  }

  /// \brief Returns variable \p varId of the cell \p cellId
  MFloat& a_variable(const MInt cellId, const MInt varId) { return m_cells.variables(cellId, varId); }

  /// \brief Returns variable \p varId of the cell \p cellId
  MFloat a_variable(const MInt cellId, const MInt varId) const { return m_cells.variables(cellId, varId); }

  /// \brief Returns oldVariables\p varId of the cell \p cellId
  MFloat& a_oldVariable(const MInt cellId, const MInt varId) { return m_cells.oldVariables(cellId, varId); }

  /// \brief Returns oldVariables \p varId of the cell \p cellId
  MFloat a_oldVariable(const MInt cellId, const MInt varId) const { return m_cells.oldVariables(cellId, varId); }

  /// \brief Returns the variables \p varId of the cell \p cellId interpolated inter/extrapolated in time
  ///        if the template parameter is false no interpolation/extrapolation takes place
  ///        dt is the time in units of deltaT at the level of cellId, from the time of variable
  ///        if(-1 < dt < 0) the result is interpolated, else extrapolated
  // template <MBool interpolate/*  = true */>
  template <MBool interpolate = true>
  MFloat a_interpolatedVariable(const MInt cellId, const MInt varId, const MFloat dt = 0.0) const {
    IF_CONSTEXPR(interpolate) { return (1.0 + dt) * a_variable(cellId, varId) - dt * a_oldVariable(cellId, varId); }
    else { // this is identical to a_variable
      return a_variable(cellId, varId);
    }
  }

  /// \brief Returns property \p p of the cell \p cellId
  MBool a_hasProperty(const MInt cellId, const Cell p) const { return grid().tree().hasProperty(cellId, p); }

  /// \brief Returns solver cell property \p p of the cell \p cellId
  maia::lb::cell::BitsetType::reference a_hasProperty(const MInt cellId, const SolverCell p) {
    return m_cells.hasProperty(cellId, p);
  }

  /// \brief Returns solver cell property \p p of the cell \p cellId
  MBool a_hasProperty(const MInt cellId, const SolverCell p) const { return m_cells.hasProperty(cellId, p); }

  /// \brief Returns bndId of the cell \p cellId
  MInt& a_bndId(const MInt cellId) { return m_cells.bndId(cellId); }

  /// \brief Returns bndId of the cell \p cellId
  MInt a_bndId(const MInt cellId) const { return m_cells.bndId(cellId); }

  // sampling data methods
  MInt getCellIdByIndex(const MInt index) { return index; }
  virtual void initInterpolationForCell(const MInt NotUsed(cellId)){};

  MFloat c_cellLengthAtLevel(const MInt level) const { return grid().cellLengthAtLevel(level); }

  MFloat c_cellLengthAtCell(const MInt cellId) const { return grid().cellLengthAtLevel(a_level(cellId)); }

  /// \brief Returns the neighbor id of the cell \p cellId for direction \p dir
  MInt& c_neighborId(const MInt cellId, const MInt dir) {
    // Note: use neighbor list instead of grid().m_tree, also includes diagonal neighbors
    return grid().neighborList(cellId, dir);
  }

  /// \brief Returns the neighbor id of the cell \p cellId for direction \p dir
  MInt c_neighborId(const MInt cellId, const MInt dir) const {
    // Note: use neighbor list instead of grid().m_tree, also includes diagonal neighbors
    return grid().neighborList(cellId, dir);
  }

  /// \brief Returns noNeighborIds of the cell \p CellId for direction \p dir
  MInt a_hasNeighbor(const MInt cellId, const MInt dir) const {
    // Note: use neighbor list instead of grid().m_tree, also includes diagonal neighbors
    return static_cast<MInt>(grid().neighborList(cellId, dir) > -1);
  }

  MBool a_isHalo(const MInt cellId) const { return m_cells.hasProperty(cellId, SolverCell::IsHalo); }

  maia::lb::cell::BitsetType::reference a_isHalo(const MInt cellId) {
    //    ENSURE_VALID_GCELLID(gCellId);
    return m_cells.hasProperty(cellId, SolverCell::IsHalo);
  }

  MBool a_isWindow(const MInt cellId) const { return m_cells.hasProperty(cellId, SolverCell::IsWindow); }

  maia::lb::cell::BitsetType::reference a_isWindow(const MInt cellId) {
    //    ENSURE_VALID_GCELLID(gCellId);
    return m_cells.hasProperty(cellId, SolverCell::IsWindow);
  }

  void a_resetPropertiesSolver(const MInt cellId) { m_cells.resetProperties(cellId); }

  MFloat* a_variables_ptr(MInt pCellId) { return &a_variable(pCellId, 0); }

  MFloat* a_oldVariables_ptr(MInt pCellId) { return &a_oldVariable(pCellId, 0); }

  MInt c_globalId(const MInt cellId) const { return grid().tree().globalId(cellId); }

  MInt a_level(const MInt cellId) const { return m_cells.level(cellId); }

  MInt& a_level(const MInt cellId) { return m_cells.level(cellId); }

  MFloat a_alphaGas(const MInt cellId) const { return m_cells.invVolumeFraction(cellId); }

  MFloat a_alphaGasLim(const MInt cellId) const { return mMax(mMin(m_cells.invVolumeFraction(cellId), 1.0), 0.0); }

  MFloat& a_alphaGas(const MInt cellId) { return m_cells.invVolumeFraction(cellId); }

  MFloat a_uOtherPhase(const MInt cellId, const MInt dir) const { return m_cells.uOtherPhase(cellId, dir); }

  MFloat& a_uOtherPhase(const MInt cellId, const MInt dir) { return m_cells.uOtherPhase(cellId, dir); }

  MInt getIdAtPoint(const MFloat* point, [[maybe_unused]] MBool globalUnique = false) {
    return grid().findContainingLeafCell(point);
  }

  ///\brief Returns the number of grid-cells
  MInt c_noCells() const { return grid().tree().size(); }

  MInt c_level(const MInt cellId) const { return grid().tree().level(cellId); }

  MInt c_noChildren(const MInt cellId) const { return grid().tree().noChildren(cellId); }

  MBool c_isLeafCell(const MInt cellId) const { return grid().tree().isLeafCell(cellId); }

  MInt c_childId(const MInt cellId, const MInt childNumber) const { return grid().tree().child(cellId, childNumber); }

  MInt c_parentId(const MInt cellId) const { return grid().tree().parent(cellId); }

  // Coupling Level-Set Accessors
  //-------------------------------------------------------------------------
  /// \brief Returns the levelSetMb-value for fv-CellId \p cellId and \p set
  MFloat& a_levelSetFunctionMB(const MInt cellId, const MInt set) {
#ifdef WAR_NVHPC_PSTL
    const MInt pos = (cellId * m_noLevelSetsUsedForMb) + set;
    return m_levelSetValues[pos];
#else
    return m_levelSetValues[IDX_LSSETMB(cellId, set)];
#endif
  }

  /// \brief Returns the levelSetMb-value for fv-CellId \p cellId and \p set
  MFloat a_levelSetFunctionMB(const MInt cellId, const MInt set) const {
#ifdef WAR_NVHPC_PSTL
    const MInt pos = (cellId * m_noLevelSetsUsedForMb) + set;
    return m_levelSetValues[pos];
#else
    return m_levelSetValues[IDX_LSSETMB(cellId, set)];
#endif
  }

  /// \brief Returns the associatedBodyIds for fv-CellId \p cellId and \p set
  MInt& a_associatedBodyIds(const MInt cellId, const MInt set) {
#ifdef WAR_NVHPC_PSTL
    const MInt pos = (cellId * m_noLevelSetsUsedForMb) + set;
    return m_associatedBodyIds[pos];
#else
    return m_associatedBodyIds[IDX_LSSETMB(cellId, set)];
#endif
  }

  /// \brief Returns the associatedBodyIds for fv-CellId \p cellId and \p set
  MInt a_associatedBodyIds(const MInt cellId, const MInt set) const {
#ifdef WAR_NVHPC_PSTL
    const MInt pos = (cellId * m_noLevelSetsUsedForMb) + set;
    return m_associatedBodyIds[pos];
#else
    return m_associatedBodyIds[IDX_LSSETMB(cellId, set)];
#endif
  }

  /// \brief Returns the level set of a G0 Candidate \p cellId
  MBool& a_isG0CandidateOfSet(const MInt cellId, const MInt set) { return m_isG0CandidateOfSet[cellId][set]; }

  /// \brief Returns the level set of a G0 Candidate \p cellId
  MBool a_isG0CandidateOfSet(const MInt cellId, const MInt set) const { return m_isG0CandidateOfSet[cellId][set]; }

  /// \brief Returns isActive of the cell \p CellId
  MBool a_isActive(const MInt cellId) const { return a_hasProperty(cellId, SolverCell::IsActive); }

  /// \brief Returns isActive of the cell \p CellId
  maia::lb::cell::BitsetType::reference a_isActive(const MInt cellId) {
    return a_hasProperty(cellId, SolverCell::IsActive);
  }

  /// \brief Returns isActive of the cell \p CellId
  MBool a_wasActive(const MInt cellId) const { return a_hasProperty(cellId, SolverCell::WasActive); }

  /// \brief Returns isActive of the cell \p CellId
  maia::lb::cell::BitsetType::reference a_wasActive(const MInt cellId) {
    return a_hasProperty(cellId, SolverCell::WasActive);
  }
  // end Coupling Level-Set Accessors
  //-------------------------------------------------------------------------


  // Swap fields (TODO labels:LB only exchange pointers)
  void swap_variables(MInt pCellId) {
#ifdef WAR_NVHPC_PSTL
    for(MInt varId = 0; varId < m_noVariables; ++varId) {
#else
    for(MInt varId = 0; varId < noVariables(); ++varId) {
#endif
      std::swap(a_variable(pCellId, varId), a_oldVariable(pCellId, varId));
    }
  }

  LbSolver(MInt id, MInt dist, GridProxy& gridProxy_, Geometry<nDim>& geometry_, const MPI_Comm comm);
  ~LbSolver() override;
  MBool isActive() const override { return grid().isActive(); }
  void updateCellCollectorFromGrid();
  MInt getCurrentTimeStep() const override { return globalTimeStep; }
  virtual void activateAllCells(){};
  virtual void activateInnerCells();
  virtual void activateWindowCells();
  virtual void activateAllButHaloCells();
  virtual void propagation_step() = 0;
  virtual void propagation_step_vol() = 0;
  virtual void propagation_step_thermal() = 0;
  virtual void propagation_step_thermal_vol() = 0;
  virtual void propagation_step_transport() = 0;
  virtual void propagation_step_transport_vol() = 0;
  virtual void propagation_step_thermaltransport() = 0;
  virtual void propagation_step_thermaltransport_vol() = 0;

  virtual void volumeForces() = 0;

  virtual void exchange(MInt mode = 1);
  virtual void exchangeLb(MInt mode);
  virtual void exchangeLbNB(MInt mode);
  virtual void exchangeOldDistributions();

  void (LbSolver::*m_sendMethod)();
  void (LbSolver::*m_receiveMethod)();
  virtual void sendNormal();
  virtual void receiveNormal();
  virtual void sendReduced();
  virtual void receiveReduced();

  // Member functions needed for adaptation
  virtual void removeChildsLb(const MInt parentId) = 0;
  virtual void refineCellLb(const MInt parentId, const MInt* childIds) = 0;
  virtual void restartBndCnd() = 0;
  virtual void prolongation() = 0;
  void writeInfo();
  void initializeRefinedCellsPerLevel();
  void initializeNewInterfaceParents();
  void saveOutput();
  void computeFFTStatistics();
  void saveAdaptedGridFile(MInt* const p_recalcCellIds);

  void (LbSolver::*m_gatherMethod)();
  void (LbSolver::*m_scatterMethod)();
  inline void gatherNormal();
  inline void scatterNormal();
  inline void gatherReduced();
  inline void scatterReduced();

  void printCommunicationMethod();
  // Coupling Level-Set Functions
  //-------------------------------------------------------------------------
  void exchangeG0Candidates();
  void preCoupleLs(std::vector<MInt>& maxGCellLevels);
  void createBndryToBodyMapping(maia::coupling::Mapping& bndryToBodyMapping,
                                maia::coupling::Mapping& bodyToBndryMapping);
  void findG0Candidates(std::vector<MInt>& maxGCellLevels);
  MInt setUpLbInterpolationStencil(const MInt cellId, MInt* interpolationCells, MFloat* point);
  MFloat interpolateFieldDataLb(MInt*, MFloat*, MInt varId, std::function<MFloat(MInt, MInt)> scalarField,
                                std::function<MFloat(MInt, MInt)> coordinate);
  void initializeMovingBoundaries();
  // End Coupling Level-Set Functions
  //-------------------------------------------------------------------------

  void resetActiveCellList(MInt mode = 0);
  void resetExternalSources();
  void exchangeExternalSources();

  void prepareCommunication();
  void prepareCommunicationNormal();
  void prepareCommunicationReduced();
  void markCellsForAdditionalComm();
  void resetComm();
  void resetCellLists(MBool resize = true);
  void loadRestartFile() override;
  virtual void saveRestartFile();
  void copyGridProperties();
  void returnCellInfo(MInt cellId);

  void initSolver() override;
  void finalizeInitSolver() override;
  void preTimeStep() override;
  MBool solutionStep() override;
  void saveSolverSolution(MBool forceOutput, const MBool finalTimeStep) override;
  virtual MBool maxResidual() {
    m_log << "Called LbSolver::maxResidual without Implementation" << std::endl;
    return true;
  }
  void outputInitSummary();

  virtual void cleanUp(){};

  void postTimeStep() {} // until now it is left unused

  void storeOldDistributions();
  void storeOldVariables();
  void initNu(const MInt cellId, const MFloat nu);

 protected:
  virtual void clb_collision_step() = 0;
  virtual void clb_smagorinsky_collision_step() = 0;
  virtual void cumulant_collision_step() = 0;
  virtual void bgki_collision_step(){};
  virtual void bgki_collision_step_mb(){};
  virtual void bgki_collision_step_mb_thermal(){};
  virtual void bgki_collision_step_Guo_forcing(){};
  virtual void bgki_init_collision_step(){};
  virtual void bgkc_collision_step(){};
  virtual void bgki_smagorinsky_collision_step(){};
  virtual void bgki_smagorinsky_collision_step2(){};
  virtual void bgki_smago_wall_collision_step(){};
  virtual void bgki_dynamic_smago_collision_step(){};
  virtual void bgki_euler_collision_step(){};
  virtual void bgki_thermal_collision_step(){};
  virtual void bgki_innerEnergy_collision_step(){};
  virtual void bgki_totalEnergy_collision_step(){};
  virtual void bgkc_transport_collision_step(){};
  virtual void bgkc_thermal_transport_collision_step(){};
  virtual void bgkc_innerenergy_transport_collision_step(){};
  virtual void bgkc_totalenergy_transport_collision_step(){};
  virtual void mrt_collision_step() = 0;
  virtual void mrt2_collision_step() = 0;
  virtual void mrt_smagorinsky_collision_step() = 0;
  virtual void rbgk_collision_step() = 0;
  virtual void rbgk_dynamic_smago_collision_step() = 0;
  virtual void rbgk_smagorinsky_collision_step() = 0;

  virtual void initializeLatticeBgk(){};
  virtual void initializeLatticeBgkThermal(){};
  virtual void initializeLatticeBgkTransport(){};
  virtual void initializeLatticeBgkThermalTransport(){};
  virtual void restartInitLb(){};

  virtual void initSrcTermController() = 0;
  virtual void initSrcTerms() = 0;
  virtual void preCollisionSrcTerm() = 0;
  virtual void postCollisionSrcTerm() = 0;
  virtual void postPropagationSrcTerm() = 0;
  virtual void postCollisionBc() = 0;
  virtual void postPropagationBc() = 0;
  virtual void updateVariablesFromOldDist() = 0;
  virtual void updateVariablesFromOldDist_preCollision() = 0;
  virtual void updateViscosity() = 0;

  virtual void calcNodalLsValues(){};

  inline void setIsActiveDefaultStates();
  inline void setInActiveBndryCells();
  inline void setInActiveMBCells();
  inline void fillActiveCellList();
  void getReLambdaAndUrmsInit();

  // sampling variables
  MBool m_isInitSamplingVars = false;
  std::vector<MFloat**> m_samplingVariables{};

 private:
  void (LbSolver::*m_propagationStepMethod)();
  void (LbSolver::*m_solutionStepMethod)();
  void (LbSolver::*m_initializeMethod)();
  void (LbSolver::*m_restartInitLbMethod)();

  // Function pointer for the exchange of data
  void (LbSolver::*m_exchangeMethod)(MInt mode);

  void initTimer();
  void averageTimer();

  // TODO labels:LB move this function to something like src/UTIL/numerics.h ?
  /** \brief  Determine coefficients for first order derivative
   *  \author Miro Gondrum
   *  \date   21.09.2022
   *  \param[in]  cellIds   cell ids of stencil sorted from negativ to positiv direction (-1 if empty)
   *  \param[out] coef      coefficients to calculate the gradient
   *  \note coefficients are based on unit spacing. Therefore, scaling with 1/dx might need to be applied.
   */
  void getDerivativeStencilAndCoefficient(std::array<MInt, 5>& cellIds, std::array<MFloat, 5>& coef) {
    // ?-?-x-?-?
    if(cellIds[0] > -1) {
      // x-x-x-?-?
      if(cellIds[4] > -1) {
        // x-x-0-x-x -> central 4th order
        coef[0] = F1B12;
        coef[1] = -F2B3;
        coef[2] = 0.0;
        coef[3] = F2B3;
        coef[4] = -F1B12;
      } else if(cellIds[3] > -1) {
        // x-x-x-x-0 -> backward + 1 forward difference 3rd order
        coef[0] = F1B6;
        coef[1] = -1.0;
        coef[2] = F1B2;
        coef[3] = F1B3;
        coef[4] = 0.0;
      } else {
        // x-x-x-0-0 -> backward difference 2nd order
        coef[0] = F1B2;
        coef[1] = -2.0;
        coef[2] = F3B2;
        coef[3] = 0.0;
        coef[4] = 0.0;
      }
    } else if(cellIds[4] > -1) {
      // 0-?-x-x-x
      if(cellIds[1] > -1) {
        // 0-x-x-x-x -> forward + 1 backward difference 3rd order
        coef[0] = 0.0;
        coef[1] = -F1B3;
        coef[2] = -F1B2;
        coef[3] = 1.0;
        coef[4] = -F1B6;
      } else {
        // 0-0-x-x-x -> forward difference 2nd order
        coef[0] = 0.0;
        coef[1] = 0.0;
        coef[2] = -F3B2;
        coef[3] = 2.0;
        coef[4] = -F1B2;
      }
    } else {
      // 0-?-x-?-0
      if(cellIds[1] > -1) {
        // 0-x-x-?-0
        if(cellIds[3] > -1) {
          // 0-x-0-x-0 -> central difference 2nd order
          coef[0] = 0.0;
          coef[1] = -F1B2;
          coef[2] = 0.0;
          coef[3] = F1B2;
          coef[4] = 0.0;
        } else {
          // 0-x-x-0-0 -> backward difference 1st order
          coef[0] = 0.0;
          coef[1] = -1.0;
          coef[2] = 1.0;
          coef[3] = 0.0;
          coef[4] = 0.0;
        }
      } else if(cellIds[3] > -1) {
        // 0-0-x-x-0 -> forward difference 1st order
        coef[0] = 0.0;
        coef[1] = 0.0;
        coef[2] = -1.0;
        coef[3] = 1.0;
        coef[4] = 0.0;
      } else {
        // 0-0-x-0-0
        // Really not possible -> all coef == 0.0 -> no gradient here
        coef[0] = 0.0;
        coef[1] = 0.0;
        coef[2] = 0.0;
        coef[3] = 0.0;
        coef[4] = 0.0;
      }
    }
  }

  void printScalingVariables();

 public:
  MBool m_calculateDissipation;

  MString m_initMethod;

  // restart
  MBool m_initFromCoarse;

  //(non-blocking) communication
  MFloat** m_sendBuffers = nullptr;
  MFloat** m_receiveBuffers = nullptr;
  MInt** m_nghbrOffsetsWindow = nullptr;
  MInt** m_nghbrOffsetsHalo = nullptr;
  MFloat** m_baseAddresses = nullptr;
  MInt m_noVarsTransfer{};
  MInt m_noElementsTransfer{};
  MInt m_noDistsTransfer{};
  MInt m_dataBlockSizeTotal{};
  MInt* m_dataBlockSizes = nullptr;

  // reduced communication
  MBool m_reducedComm{};
  MInt* m_noWindowDistDataPerDomain{};
  MInt* m_noHaloDistDataPerDomain{};
  MFloat*** m_commPtWindow{};
  MFloat*** m_commPtHalo{};
  MInt* m_orderedNeighbors{};
  MInt* m_needsFurtherExchange{};
  MInt** m_windowDistsForExchange{};
  MInt** m_haloDistsForExchange{};

  MBool m_isRefined = false;
  MFloat m_interfaceCellSize{};
  MBool nonBlockingComm() const { return m_nonBlockingComm; }

  // Load balancing
  MBool m_setCellDataFinished = false;

  // Coupling Level-Set
  MInt m_currentNoG0Cells;
  MFloat** m_nodalGValues{};
  MInt* m_G0CellList{};
  MInt* m_G0CellMapping{};
  MFloat** m_G0WallDist{};
  std::vector<CutCandidate<nDim>> m_G0Candidates;
  MInt* m_sendBufferMB{};
  MInt* m_receiveBufferMB{};
  MInt m_maxNoSets = -1;
  MInt m_levelSetId = 0;
  MFloat* m_levelSetValues{};
  MInt* m_associatedBodyIds{};
  MBool** m_isG0CandidateOfSet{};
  MInt m_noLevelSetsUsedForMb;
  MFloat m_maxBodyRadius;
  MInt m_noEmbeddedBodies;
  MBool m_constructGField;
  MBool m_useOnlyCollectedLS = false;
  MBool m_allowBndryAsG0 = false;

  MFloat* m_bodyVelocity = nullptr;
  MBool m_trackMovingBndry;
  MInt m_trackMbEnd;
  MInt m_trackMbStart;
  MInt m_noG0CandidatesTotal;
  MInt* m_initialActiveCells = nullptr;

  // Residual
  MInt m_resTimestepPos;
  MInt m_resRePos;
  std::ofstream mRes; // output stream for the residual file

  MInt m_solutionOffset;
  MFloat m_referenceLength;
  MFloat m_referenceLengthSTL;

  MFloat m_noOuterBndryCells;

  // Non-Newtonian
  MFloat m_n = F0;
  MBool m_nonNewtonian = false;

  /// Return number of variables
  MInt noVariables() const override { return m_noVariables; }
  MInt noDistributions() const { return m_noDistributions; }
  MFloat time() const override { return m_time; }

  using Geom = Geometry<nDim>;

  /// Access the solver's geometry
  constexpr const Geom& geometry() const { return *m_geometry; }

  GeometryIntersection<nDim>* m_geometryIntersection;


 protected:
  MBool m_nonBlockingComm = false;
  /// Access the solver's geometry (non-const version)
  Geom& geometry() { return *m_geometry; }

 private:
  // Number of variables
  MInt m_noVariables;


 public:
  /// [Splitt] The following is part of a first step to splitt CartesianGrid
  /// from the inheritance hierarchy:
  ///
  /// - in order to avoid renaming a lot of access to CartesianGrid data
  ///   members, references are introduced:
  ///
  /// \todo labels:IO this references will be removed in future commits
  Collector<PointBasedCell<nDim>>* m_extractedCells = nullptr;
  Collector<CartesianGridPoint<nDim>>* m_gridPoints = nullptr;

  /// the references to CartesianGrid data members end here

  MInt a_noCells() const { return m_cells.size(); }

  MFloat a_Nu() const { return m_nu; }
  MInt a_FFTInterval() const { return m_fftInterval; }

 private:
  MPI_Request* mpi_request{};
  MPI_Request* mpi_requestR{};
  MPI_Request* mpi_requestS{};

 protected:
  Geometry<nDim>* m_geometry;
  /// Collector for LB cells
  maia::lb::collector::LbCellCollector<nDim> m_cells;

  // Pointer to the class MConservativeVariables
  MConservativeVariables<nDim>* CV;
  // Pointer to the class MPrimitiveVariables
  MPrimitiveVariables<nDim>* PV;

  // domain dimensions
  MInt m_arraySize[nDim];       // in cell units
  MFloat* m_domainBoundaries{}; // in macroscopic units

  LbBndCnd<nDim>* m_bndCnd; ///< Pointer to the boundary conditions

  // subgrid properties
  MFloat m_Cs;
  MFloat m_deltaX;

  // time (non-dimensional by using L=1m and u=u_infinity)
  MFloat m_time = F0; ///< Non-dimensional time (accumulated m_dt)
  MFloat m_dt = -1.0; ///< Non-dimensional time step : dt= dt_phys * u_inf / 1

  MFloat m_densityGradient{};

  MFloat m_domainLength;
  MInt maxTwoPower{};
  MInt m_currentMaxNoCells{};
  MInt m_initialNoCells[3]{};

  MBool m_densityFluctuations;
  MInt m_noPeakModes;
  MBool m_FftInit;

  // FFT output and number of timesteps before wrinting next FFT output
  MInt m_fftInterval;

  // Defines the discretization Model
  MInt m_noDistributions;
  MInt m_methodId;

  // Viscosity
  MFloat m_nu;

  // for BCs
  MFloat m_rho1;
  MFloat m_rho2;

  // pulsatile
  MFloat m_alpha;
  MFloat m_pulsatileFrequency;

  MBool m_isInitRun = false;

  MBool m_tanhInit;
  MFloat m_initRe;
  MInt m_initTime;
  MInt m_initStartTime;
  MFloat m_finalRe;
  MFloat m_tanhScaleFactor;

  // for isotropic turbulence
  MFloat m_UrmsInit;
  MFloat m_ReLambda;

  // Relaxation factor
  MFloat m_omega;

  // external Forcing terms
  MFloat* m_Fext{};

  MFloat m_Ga{};

  // dimension independent data
  //! Check!
  MInt* m_startLevel{};

  // Use an external pressure force
  MBool m_externalForcing;

  MBool m_particleMomentumCoupling;

  MBool m_saveExternalForces;

  MBool m_updateAfterPropagation;
  LbUpdateMacroscopic m_updateMacroscopicLocation = INCOLLISION;
  MBool m_initDensityGradient;
  MFloat m_volumeAccel[3] = {};
  MFloat m_volumeAccelBase[3] = {};
  MBool m_controlVelocity = false;
  struct {
    MInt dir;
    MInt interval;
    MFloat KT;
    MFloat KI;
    MFloat KD;
    MFloat lastGlobalAvgV = -1.0;
    MFloat integratedError = 0.0;
    MFloat derivedError = 0.0;
    MFloat previousError = 0.0;
    MBool restart = false;
  } m_velocityControl;

  MBool m_isEELiquid;
  struct {
    MBool restartWithoutAlpha;
    MFloat alphaInf;
    MFloat initialAlpha;
    MBool gravity;
    MFloat gravityAccelM[3];
    MFloat* Fg{};
  } m_EELiquid;

  MFloat* m_residual = nullptr;
  MFloat* m_tmpResidual = nullptr;
  MInt* m_tmpResidualLvl = nullptr;
  MInt* m_maxResId = nullptr;
  MFloat** m_rescoordinates = nullptr;
  MString m_resFileName = "";

  MFloat** m_momentumFlux = nullptr;
  MFloat* m_MijMij = nullptr;
  MFloat* m_MijLij = nullptr;

  MFloat m_ReTau;

  // the id of the segment used to get m_referenceLength
  MInt m_referenceLengthSegId;
  MBool m_refineDiagonals = true;

  // Thermal variables
  MFloat m_Pr;
  MFloat m_omegaT;
  MFloat m_kappa;
  MFloat m_initTemperatureKelvin;
  MFloat m_blasiusPos;

  MBool m_isCompressible = false;
  MInt m_isThermal = 0;
  MInt m_innerEnergy = 0;
  MInt m_totalEnergy = 0;
  MFloat m_CouettePoiseuilleRatio = F0;
  MInt m_calcTotalPressureGradient = 0;

  // Transport variables
  MFloat m_Pe;
  MFloat m_omegaD;
  MFloat m_diffusivity;
  MFloat m_initCon;
  MInt m_isTransport = 0;

  MFloat* m_cellLength = nullptr;

  virtual void treatInterfaceCells();
  void correctInterfaceBcCells();
  virtual void setActiveCellList();

  /// Grid interface cells
  virtual void initializeInterfaceCells();
  virtual void buildInterfaceCells();
  virtual void resetInterfaceCells();

  LbInterpolationType m_interpolationType;
  virtual void setInterpolationNeighbors();
  virtual void setInterpolationCoefficients();
  void setInterpolationNeighborsBC();
  void setInterpolationCoefficientsBC();

  std::vector<Collector<LbInterfaceCell>*> m_interfaceChildren;
  std::vector<Collector<LbParentCell>*> m_interfaceParents;

  MInt* m_noInterfaceChildren = nullptr;
  MInt* m_noInterfaceParents = nullptr;

  MBool m_correctInterfaceBcCells = false;
  struct InterfaceLink {
    MInt levelId;
    MInt interfaceId;
  };
  std::vector<InterfaceLink> m_interfaceChildrenBc;

  MBool m_saveDerivatives;
  MBool m_initRestart;

  // The number of dimensions
  MInt* m_activeCellList = nullptr;
  MInt* m_activeCellListLvlOffset = nullptr;

 protected:
  virtual void writeGridToTecFile(const MChar* fileName);
  virtual void writeGridToVtkFile(const MChar* fileName);
  virtual void writeGridToVtkFileAscii(const MChar* fileName);
  virtual void initPressureForce() = 0;
  virtual void initVolumeForces() = 0;
  virtual void initRunCorrection() = 0;
  void determineEqualDiagonals();

 public: // parallel IO
  void saveUVWRhoTOnlyPar(const MChar* fileName, const MChar* gridInputFileName, MInt* recalcIdTree = nullptr);
  void saveRestartWithDistributionsPar(const MChar* fileName, const MChar* gridInputFileName,
                                       MInt* recalcIdTree = nullptr);
  void loadRestartWithDistributionsPar(const MChar* fileName);
  void loadRestartWithoutDistributionsPar(const MChar* fileName);
  void loadRestartWithoutDistributionsParFromCoarse(const MChar* fileName);
  // derivatives
  void calculateVelocityDerivative(const MInt cellId, MFloat (&gradient)[nDim][nDim]);
  template <MBool tI = false> // interpolation in time
  MFloat calculateDerivative(const MInt cellId, const MInt velComp, const MInt spaceDir, const MFloat dt = 0.0);
  MFloat calculateInterpolatedDerivative(const MInt cellId, const MInt velComp, const MInt spaceDir, const MFloat dt) {
    return calculateDerivative<true>(cellId, velComp, spaceDir, dt);
  }
  inline MFloat calculatePressureDerivative(MInt cellId, MInt spaceDir);
  inline void getStrainTensor(MFloatScratchSpace& derivatives, MInt cellId, MFloatScratchSpace& strain);
  inline void getVorticityTensor(MFloatScratchSpace& derivatives, MInt cellId, MFloatScratchSpace& strain);

  // postprocessing
  void accessSampleVariables(MInt cellId, MFloat*& vars);
  void interpolateToParentCells(MInt parentlevel);
  void saveCoarseSolution();

  void createMPIComm(const MInt* ranks, MInt noRanks, MPI_Comm* comm);
  inline void writeSegmentBoundaryVTK(MFloat** bndVs, MInt num);
  void updateTime();

 protected:
  virtual MFloat calculateReferenceLength(MInt segmentId);
  inline MFloat calcCharLenAll(MInt segmentId);
  inline MFloat calcCharLenParallelSplit(MInt segmentId);

 public:
  virtual void getSolverSamplingProperties(std::vector<MInt>& samplingVars, std::vector<MInt>& noSamplingVars,
                                           std::vector<std::vector<MString>>& samplingVarNames,
                                           const MString featureName = "");
  virtual void initSolverSamplingVariables(const std::vector<MInt>& varIds, const std::vector<MInt>& noSamplingVars);
  virtual void calcSamplingVariables(const std::vector<MInt>& varIds, const MBool exchange);
  void interpolateVariablesInCell(const MInt cellId, const MFloat* position, MFloat* interpolationResult);

  void calcSamplingVarAtPoint(const MFloat* point, const MInt id, const MInt sampleVarId, MFloat* state,
                              const MBool interpolate = false);
  void loadSampleVariables(MInt timeStep);
  void getSampleVariables(MInt cellId, const MFloat*& vars);
  void getSampleVariables(const MInt cellId, std::vector<MFloat>& vars);
  MBool getSampleVarsDerivatives(const MInt cellId, std::vector<MFloat>& vars);
  void calculateHeatRelease() { std::cerr << "calculateHeatRelease LbSolver " << std::endl; }
  void getHeatRelease(MFloat*& heatRelease) { std::cerr << "getHeatRelease LbSolver " << heatRelease << std::endl; }

  // adaptation
  void prepareAdaptation(std::vector<std::vector<MFloat>>& NotUsed(sensors),
                         std::vector<MFloat>& NotUsed(sensorWeight),
                         std::vector<std::bitset<64>>& NotUsed(sensorCellFlag),
                         std::vector<MInt>& NotUsed(sensorSolverId)) override{};
  void prepareAdaptation() override;
  void setSensors(std::vector<std::vector<MFloat>>& sensors, std::vector<MFloat>& sensorWeight,
                  std::vector<std::bitset<64>>& sensorCellFlag, std::vector<MInt>& sensorSolverId) override;
  void reinitAfterAdaptation() override{};
  void postAdaptation() override;
  void finalizeAdaptation() override;
  void removeChilds(const MInt) override;
  void removeCell(const MInt) override;
  void refineCell(const MInt) override;
  void swapCells(const MInt, const MInt) override;
  void swapProxy(const MInt cellId0, const MInt cellId1) override;
  void resizeGridMap() override;
  MBool prepareRestart(MBool writeRestart, MBool& writeGridRestart) override;
  void reIntAfterRestart(MBool doneRestart);
  void writeRestartFile(const MBool writeRestart, const MBool writeBackup, const MString gridFileName,
                        MInt* recalcIdTree) override;

  // Dynamic load balancing
  void balance(const MInt* const /*noCellsToReceiveByDomain*/, const MInt* const /*noCellsToSendByDomain*/,
               const MInt* const /*targetDomainsByCell*/, const MInt /*oldNoCells*/) override {
    TERMM(1, "Use split balancing methods for LB solver instead of balance().");
  };
  void setCellWeights(MFloat* solverCellWeight) override;
  void resetSolver() override{};
  void finalizeBalance() override{};

  // DLB loads
  MInt noLoadTypes() const override;
  void getLoadQuantities(MInt* const loadQuantities) const override;
  void getDefaultWeights(MFloat* weights, std::vector<MString>& names) const override;
  MFloat getCellLoad(const MInt cellId, const MFloat* const weights) const override;

  // DLB (Split balance)
  MBool hasSplitBalancing() const override { return true; }
  void balancePre() override;
  void balancePost() override;

  // DLB exchange information (Split balance)

  // Number of fields to exchange during DLB
  MInt noCellDataDlb() const override {
    if(grid().isActive()) {
      return noSolverCellData();
    } else {
      return 0;
    }
  };

  // Datatype (enum) for a given field
  MInt cellDataTypeDlb(const MInt dataId) const {
    TRACE();

    if(!isActive()) {
      TERMM(1, "Error: cellDataTypeDlb() might give wrong results on inactive ranks.");
      return -1;
    }

    MInt dataType = -1;
    if(dataId > -1 && dataId < noSolverCellData()) {
      // LB solver cell data
      dataType = s_cellDataTypeDlb[dataId];
    } else {
      TERMM(1, "The requested dataId is not valid: " + std::to_string(dataId) + " ("
                   + std::to_string(noSolverCellData()) + ", " + std::to_string(noCellDataDlb()) + ")");
    }

    return dataType;
  }

  // Data size per cell
  MInt cellDataSizeDlb(const MInt dataId, const MInt gridCellId) {
    TRACE();

    // Inactive ranks do not have any data to communicate
    if(!isActive()) {
      return 0;
    }

    // Convert to solver cell id and check
    const MInt cellId = grid().tree().grid2solver(gridCellId);
    if(cellId < 0) {
      return 0;
    }

    MInt dataSize = 0;

    if(dataId > -1 && dataId < noSolverCellData()) {
      // LB solver cell data
      switch(dataId) {
        case CellData::VARIABLES:
        case CellData::OLD_VARIABLES:
          dataSize = noVariables();
          break;
        case CellData::DISTRIBUTIONS:
        case CellData::OLD_DISTRIBUTIONS:
          dataSize = noDistributions();
          break;
        case CellData::NU:
        case CellData::IS_ACTIVE:
          dataSize = 1;
          break;
        default:
          TERMM(1, "Unknown data id. (" + std::to_string(dataId) + ")");
          break;
      }
    } else {
      TERMM(1, "The requested dataId is not valid.");
    }

    return dataSize;
  }

  // Copy cell data field to buffer (int version)
  void getCellDataDlb(const MInt dataId, const MInt oldNoCells, const MInt* const bufferIdToCellId,
                      MInt* const data) override {
    getCellDataDlb_(dataId, oldNoCells, bufferIdToCellId, data);
  };

  // Copy cell data field to buffer (float version)
  void getCellDataDlb(const MInt dataId, const MInt oldNoCells, const MInt* const bufferIdToCellId,
                      MFloat* const data) override {
    getCellDataDlb_(dataId, oldNoCells, bufferIdToCellId, data);
  };

  // Copy buffer data to data field (int version)
  void setCellDataDlb(const MInt dataId, const MInt* const data) override { setCellDataDlb_(dataId, data); };

  // Copy buffer data to data field (float version)
  void setCellDataDlb(const MInt dataId, const MFloat* const data) override { setCellDataDlb_(dataId, data); };

  /// \brief Set solver cell data from buffer
  ///
  /// \author Julian Vorspohl <j.vorspohl@aia.rwth-aachen.de>
  ///
  /// This method sets the given data for all cells, e.g. variables, distributions.
  ///
  /// \tparam DataType Current fields data type.
  /// \param[in] dataId Requested data id.
  /// \param[in] data Pointer to storage of data.
  template <typename DataType>
  void setCellDataDlb_(const MInt dataId, const DataType* const data) {
    TRACE();

#ifdef LBCOLLECTOR_SOA_MEMORY_LAYOUT
    TERMM(1, "Balance data accessor not defined for SOA cell collector");
#endif

    // This function can be called multiple times, prevent override
    if(m_setCellDataFinished) {
      return;
    }

    // Nothing to do if solver is not active
    if(!grid().isActive()) {
      return;
    }

    const MInt noCells = m_cells.size();

    if(dataId > -1 && dataId < noSolverCellData()) {
      // LB solver cell data
      switch(dataId) {
        case CellData::VARIABLES:
          std::copy_n(data, noCells * noVariables(), &m_cells.variables(0, 0));
          break;
        case CellData::OLD_VARIABLES:
          std::copy_n(data, noCells * noVariables(), &m_cells.oldVariables(0, 0));
          break;
        case CellData::DISTRIBUTIONS:
          std::copy_n(data, noCells * noDistributions(), &m_cells.distributions(0, 0));
          break;
        case CellData::OLD_DISTRIBUTIONS:
          std::copy_n(data, noCells * noDistributions(), &m_cells.oldDistributions(0, 0));
          break;
        case CellData::NU:
          std::copy_n(data, noCells, &m_cells.nu(0));
          break;
        case CellData::IS_ACTIVE:
          for(MInt i = 0; i < noCells; i++) {
            a_isActive(i) = (MInt)data[i];
          }
          break;
        default:
          TERMM(1, "Unknown data id (" + std::to_string(dataId) + ").");
          break;
      }
    } else {
      TERMM(1, "Invalid dataId: " + std::to_string(dataId));
    }
  }

  /// \brief Get solver cell data for load balancing.
  ///
  /// \author Julian Vorspohl <j.vorspohl@aia.rwth-aachen.de>
  ///
  /// This method returns the requested data for all cells, e.g. variables, distributions
  ///
  /// \tparam DataType Current fields data type.
  /// \param[in] dataId Requested data id.
  /// \param[in] oldNoCells Current (old) number of cells before load balancing.
  /// \param[in] bufferIdToCellId Mapping from buffer location to corresponding cell id.
  /// \param[out] data Pointer to storage for requested data.
  template <typename DataType>
  void getCellDataDlb_(const MInt dataId, const MInt oldNoCells, const MInt* const bufferIdToCellId,
                       DataType* const data) {
    TRACE();

#ifdef LBCOLLECTOR_SOA_MEMORY_LAYOUT
    TERMM(1, "Balance data accessor not defined for SOA cell collector");
#endif

    ASSERT((dataId > -1 && dataId < noSolverCellData()), "Invalid dataId in getCellData!");

    MInt localBufferId = 0;
    for(MInt i = 0; i < oldNoCells; i++) {
      const MInt gridCellId = bufferIdToCellId[i];
      if(gridCellId < 0) {
        continue;
      }

      const MInt cellId = grid().tree().grid2solver(gridCellId);
      if(cellId < 0 || cellId >= noInternalCells()) {
        continue;
      }

      // LB solver cell data
      switch(dataId) {
        case CellData::VARIABLES: {
          const MInt dataSize = noVariables();
          std::copy_n(&m_cells.variables(cellId, 0), dataSize, &data[localBufferId * dataSize]);
          break;
        }
        case CellData::OLD_VARIABLES: {
          const MInt dataSize = noVariables();
          std::copy_n(&m_cells.oldVariables(cellId, 0), dataSize, &data[localBufferId * dataSize]);
          break;
        }
        case CellData::DISTRIBUTIONS: {
          const MInt dataSize = noDistributions();
          std::copy_n(&m_cells.distributions(cellId, 0), dataSize, &data[localBufferId * dataSize]);
          break;
        }
        case CellData::OLD_DISTRIBUTIONS: {
          const MInt dataSize = noDistributions();
          std::copy_n(&m_cells.oldDistributions(cellId, 0), dataSize, &data[localBufferId * dataSize]);
          break;
        }
        case CellData::NU: {
          const MInt dataSize = 1;
          std::copy_n(&m_cells.nu(cellId), dataSize, &data[localBufferId * dataSize]);
          break;
        }
        case CellData::IS_ACTIVE: {
          const MInt dataSize = 1;
          data[localBufferId * dataSize] = a_isActive(cellId);
          break;
        }
        default:
          TERMM(1, "Unknown data id (" + std::to_string(dataId) + ").");
          break;
      }
      localBufferId++;
    }
  }

  // Helper struct to differentiate cell data which needs to be exchanged during DLB
  // TODO labels:LB Does not support thermal, change if proper sysEqn is introduced.
  struct CellData {
    static constexpr const MInt count = 6;

    static constexpr const MInt VARIABLES = 0;
    static constexpr const MInt OLD_VARIABLES = 1;
    static constexpr const MInt DISTRIBUTIONS = 2;
    static constexpr const MInt OLD_DISTRIBUTIONS = 3;
    static constexpr const MInt NU = 4;
    static constexpr const MInt IS_ACTIVE = 5;
  };

  // Data types of cell data
  const std::array<MInt, CellData::count> s_cellDataTypeDlb = {{MFLOAT, MFLOAT, MFLOAT, MFLOAT, MFLOAT, MINT}};

  // Number of different data fields
  MInt noSolverCellData() const { return CellData::count; };

 protected:
  // Timers
  struct {
    MInt solver;
    MInt initSolver;
    MInt solutionStep;
    MInt collision;
    MInt collisionBC;
    MInt propagation;
    MInt propagationBC;
    MInt srcTerms;
    MInt exchange;
    MInt residual;
    MInt packing;
    MInt unpacking;
    MInt communication;

    MInt findG0Cells;
    MInt resetListsMb;
    MInt findG0Candidates;
    MInt geomNodal;
    MInt geomExchange;
    MInt calcNodalValues;

    MInt prepComm;

    MInt fft;
  } m_t;

  // Status flags for adaptation
 public:
  MBool m_adaptationSinceLastRestart;
  MBool m_adaptationSinceLastSolution;
  MBool m_adaptation;

#ifdef WAR_NVHPC_PSTL
  std::array<MFloat, 50> m_faculty;
  std::array<MInt, 27> m_nFld;
  std::array<MInt, 27> m_pFld;
  MInt m_idFld[POWX(3, nDim)][nDim];
  std::array<MInt, 24> m_mFld1;
  std::array<MInt, 24> m_mFld2;
  std::array<MInt, 27> m_oppositeDist;
  std::array<MFloat, 81> m_ppdfDir;
  std::array<MFloat, 4> m_tp;
  std::array<MInt, 3> m_distFld;
  std::array<MInt, 27> m_distType;
#endif

  // Members for adaptation of grid file name
  MString m_reinitFileName; // Current grid file name.
  MString m_reinitFilePath; // Path to current grid file.
 protected:
  // Properties for adaptation.
  MString m_adaptationInitMethod;
  std::set<MInt> m_refinedParents;
  MBool m_solidLayerExtension = false;
  MBool m_writeLsData = false;
};

/** \brief Calculates a spatial derivative for a variable
 *
 * \author Andreas Lintermann, Moritz Waldmann
 * \date 10.03.2016
 *
 * This function uses:
 *  - 4th order central differences if both neighbors and their neighbors in spaceDir are available
 *  - 3rd order forward or backward differences if both neighbors and one neighbors neighbor in spaceDir are available
 *  - 2nd order central differences if both neighbors in spaceDir are available
 *  - 1st order forward or backward differences if only one neighbor is available
 *  - zero gradient if no neighbor is available
 *
 * \param[in] cellId the cell id
 * \param[in] comp the component
 * \param[in] spaceDir the spatial direction
 * \param[in] dt time for the interpolation, only used if template parameter timeInterpolation is true
 * \return the derivative
 *
 **/
template <MInt nDim>
template <MBool tI> // interpolation in time, added by Daniel Lauwers
MFloat LbSolver<nDim>::calculateDerivative(const MInt cellId, const MInt comp, const MInt spaceDir, const MFloat dt) {
  TRACE();
  const MInt lbdir1 = 2 * spaceDir;
  const MInt lbdir2 = lbdir1 + 1;

  const MInt left = c_neighborId(cellId, lbdir1);
  const MInt right = c_neighborId(cellId, lbdir2);
  const MInt leftleft = (left > -1) ? c_neighborId(left, lbdir1) : -1;
  const MInt rightright = (right > -1) ? c_neighborId(right, lbdir2) : -1;

  MFloat gradient = F0;

  if((left > -1) && (right > -1)) {
    if(leftleft > -1 && rightright > -1) { // use central differences 4th order
      gradient = (-a_interpolatedVariable<tI>(rightright, comp, dt) + 8.0 * a_interpolatedVariable<tI>(right, comp, dt)
                  - 8.0 * a_interpolatedVariable<tI>(left, comp, dt) + a_interpolatedVariable<tI>(leftleft, comp, dt))
                 / 12.0;
    } else {
      if(leftleft > -1) { // backward differences 3nd order
        gradient =
            (2.0 * a_interpolatedVariable<tI>(right, comp, dt) + 3.0 * a_interpolatedVariable<tI>(cellId, comp, dt)
             - 6.0 * a_interpolatedVariable<tI>(left, comp, dt) + a_interpolatedVariable<tI>(leftleft, comp, dt))
            / 6.0;
      } else if(rightright > -1) { // forward differences 3nd order
        gradient =
            (-a_interpolatedVariable<tI>(rightright, comp, dt) + 6.0 * a_interpolatedVariable<tI>(right, comp, dt)
             - 3.0 * a_interpolatedVariable<tI>(cellId, comp, dt) - 2.0 * a_interpolatedVariable<tI>(left, comp, dt))
            / 6.0;
      } else { // use central differences 2nd order
        gradient = (a_interpolatedVariable<tI>(right, comp, dt) - a_interpolatedVariable<tI>(left, comp, dt)) / 2.0;
      }
    }
  } else { // use forward or backward differences 1st order
    // backward
    if(left > -1) {
      gradient = a_interpolatedVariable<tI>(cellId, comp, dt) - a_interpolatedVariable<tI>(left, comp, dt);
    }
    // forward
    else if(right > -1) {
      gradient = a_interpolatedVariable<tI>(right, comp, dt) - a_interpolatedVariable<tI>(cellId, comp, dt);
    }
  }
  return gradient;
}

/** \brief  Calculate velocity derivative
 *  \author Miro Gondrum
 *  \date   08.02.2022
 *  \param[in]  cellId    cell id for which gradient is determined
 *  \param[out] gradient  velocity gradient (p.e. gradient[1][2] is dv/dz)
 */
template <MInt nDim>
void LbSolver<nDim>::calculateVelocityDerivative(const MInt cellId, MFloat (&gradient)[nDim][nDim]) {
  TRACE();
  const MFloat F1bdx = FFPOW2(maxLevel() - a_level(cellId)); // LB units
  for(MInt dir = 0; dir < nDim; dir++) {
    const MInt dir0 = 2 * dir;
    const MInt dir1 = 2 * dir + 1;

    auto isValid = [&](const MInt l_cellId, const MInt l_dir) -> MInt {
      if(l_cellId != -1) {
        const MInt nghbrId = c_neighborId(l_cellId, l_dir);
        if(nghbrId > -1 && a_isActive(nghbrId)) {
          return nghbrId;
        }
      }
      return -1;
    };
    // 1) set coefficients
    constexpr MInt noCellsInStencil = 5;
    std::array<MFloat, noCellsInStencil> coef;
    std::array<MInt, noCellsInStencil> cellIds;
    // 1.1) get cellIds and coefficients
    cellIds[1] = isValid(cellId, dir0);     // left
    cellIds[0] = isValid(cellIds[1], dir0); // leftleft
    cellIds[2] = cellId;                    // middle
    cellIds[3] = isValid(cellId, dir1);     // right
    cellIds[4] = isValid(cellIds[3], dir1); // rightright
    getDerivativeStencilAndCoefficient(cellIds, coef);
    // 1.2) scale coefficients with correct dx in case of refined
    for(MInt i = 0; i < noCellsInStencil; i++) {
      coef[i] *= F1bdx;
    }
    // 2) calculate gradient
    for(MInt veloId = 0; veloId < nDim; veloId++) {
      gradient[veloId][dir] = 0.0;
      for(MInt i = 0; i < noCellsInStencil; i++) {
        gradient[veloId][dir] += cellIds[i] < 0 ? 0.0 : coef[i] * a_variable(cellIds[i], veloId);
      }
    }
  }
}

#endif
