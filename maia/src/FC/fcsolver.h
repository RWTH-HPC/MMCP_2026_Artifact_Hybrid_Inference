// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef FCSOLVER_H
#define FCSOLVER_H

#include "GEOM/geometryintersection.h"
#include "GRID/cartesiangrid.h"
#include "GRID/cartesiangridpoint.h"
#include "GRID/cartesianpointbasedcell.h"
#include "POST/postprocessing.h"
#include "UTIL/functions.h"
#include "cartesiansolver.h"
#include "fcbndrycnd.h"
#include "fccellcollector.h"
#include "fcdescriptor.h"
#include "solver.h"
#include "variables.h"

/** \brief This class represents a structure solver using the Finite Cell Method
 *
 */
template <MInt nDim_>
class FcSolver : public maia::CartesianSolver<nDim_, FcSolver<nDim_>> {
  template <MInt nDim>
  friend class FcBndryCnd;

 public:
  static constexpr MInt nDim = nDim_;

  // Type for cell properties
  using Cell = typename maia::grid::tree::Tree<nDim>::Cell;
  using SolverCell = FcCell;
  using CartesianSolver = typename maia::CartesianSolver<nDim, FcSolver>;
  using Grid = typename CartesianSolver::Grid;
  using GridProxy = typename CartesianSolver::GridProxy;

  // used solver
  using CartesianSolver::domainId;
  using CartesianSolver::domainOffset;
  using CartesianSolver::grid;
  using CartesianSolver::haloCell;
  using CartesianSolver::haloCellId;
  using CartesianSolver::m_bandWidth;
  using CartesianSolver::m_freeIndices;
  using CartesianSolver::m_innerBandWidth;
  using CartesianSolver::m_Ma;
  using CartesianSolver::m_outerBandWidth;
  using CartesianSolver::m_Re;
  using CartesianSolver::m_residualInterval;
  using CartesianSolver::m_restartFile;
  using CartesianSolver::m_restartInterval;
  using CartesianSolver::m_solutionInterval;
  using CartesianSolver::m_solverId;
  using CartesianSolver::maxLevel;
  using CartesianSolver::maxNoGridCells;
  using CartesianSolver::maxUniformRefinementLevel;
  using CartesianSolver::minLevel;
  using CartesianSolver::mpiComm;
  using CartesianSolver::neighborDomain;
  using CartesianSolver::noDomains;
  using CartesianSolver::noHaloCells;
  using CartesianSolver::noNeighborDomains;
  using CartesianSolver::noWindowCells;
  using CartesianSolver::outputDir;
  using CartesianSolver::reductionFactor;
  using CartesianSolver::restartDir;
  using CartesianSolver::solverId;
  using CartesianSolver::solverMethod;
  using CartesianSolver::updateDomainInfo;
  using CartesianSolver::windowCell;
  using CartesianSolver::windowCellId;
  using Fd = FcDescriptor<nDim>;

  MInt noInternalCells() const override { return grid().noInternalCells(); }

  /// \brief Returns isBndryCell of the cell \p CellId
  MBool a_isBndryCell(const MInt cellId) const { return a_hasProperty(cellId, SolverCell::IsBndryCell); }

  /// \brief Returns isBndryCell of the cell \p CellId
  maia::fc::cell::BitsetType::reference a_isBndryCell(const MInt cellId) {
    return a_hasProperty(cellId, SolverCell::IsBndryCell);
  }

  /// \brief Returns isInterface of the cell \p CellId
  MBool a_isWindow(const MInt cellId) const { return m_cells.hasProperty(cellId, SolverCell::IsWindow); }

  /// \brief Returns isInterface of the cell \p CellId
  maia::fc::cell::BitsetType::reference a_isWindow(const MInt cellId) {
    return m_cells.hasProperty(cellId, SolverCell::IsWindow);
  }

  /// \brief Returns alpha of the cell \p cellId
  MFloat& a_alpha(const MInt cellId) { return m_cells.alpha(cellId); }

  /// \brief Returns alpha of the cell \p cellId
  MFloat a_alpha(const MInt cellId) const { return m_cells.alpha(cellId); }

  /// \brief Returns poisson Ratio of the cell \p cellId
  MFloat& a_poissonRatio(const MInt cellId) { return m_cells.poissonRatio(cellId); }

  /// \brief Returns poisson Ratio of the cell \p cellId
  MFloat a_poissonRatio(const MInt cellId) const { return m_cells.poissonRatio(cellId); }

  /// \brief Returns inverse Jacobian of the cell \p cellId
  MFloat& a_invJacobian(const MInt cellId) { return m_cells.invJacobian(cellId); }

  /// \brief Returns inverse Jacobian of the cell \p cellId
  MFloat a_invJacobian(const MInt cellId) const { return m_cells.invJacobian(cellId); }

  /// \brief Returns the element \p dir from the legendre points of the cell \p cellId
  MFloat& a_nodePosition(const MInt cellId, const MInt dir) { return m_cells.nodePosition(cellId, dir); }

  /// \brief Returns the element \p dir from the legendre points of the cell \p cellId
  MFloat a_nodePosition(const MInt cellId, const MInt dir) const { return m_cells.nodePosition(cellId, dir); }

  /// \brief Returns delta gamma of node \p node of the cell \p cellId
  MFloat& a_dGamma(const MInt cellId, const MInt node) { return m_cells.deltaGamma(cellId, node); }

  /// \brief Returns delta gamma of node \p node of the cell \p cellId
  MFloat a_dGamma(const MInt cellId, const MInt node) const { return m_cells.deltaGamma(cellId, node); }

  /// \brief Returns epsilon bar of node \p node of the cell \p cellId
  MFloat& a_epsilonBarP(const MInt cellId, const MInt node) { return m_cells.epsilonBarP(cellId, node); }

  /// \brief Returns epsilon bar of node \p node of the cell \p cellId
  MFloat a_epsilonBarP(const MInt cellId, const MInt node) const { return m_cells.epsilonBarP(cellId, node); }

  /// \brief Returns the displacement \p dir of the cell \p cellId
  MFloat& a_elementDispl(const MInt cellId, const MInt dir) { return m_cells.elementDisplacements(cellId, dir); }

  /// \brief Returns the element displacement \p dir of the cell \p cellId
  MFloat a_elementDispl(const MInt cellId, const MInt dir) const { return m_cells.elementDisplacements(cellId, dir); }

  /// \brief Returns the element strain \p dir from the stain matrix of the cell \p cellId
  MFloat& a_elementStrains(const MInt cellId, const MInt dir) { return m_cells.elementStrains(cellId, dir); }

  /// \brief Returns the elment strain \p dir from the strain matrix of the cell \p cellId
  MFloat a_elementStrains(const MInt cellId, const MInt dir) const { return m_cells.elementStrains(cellId, dir); }

  /// \brief Returns the nodal strains \p dir from the stain matrix of node \p node of cell \p cellId
  MFloat& a_nodalStrains(const MInt cellId, const MInt node, const MInt dir) {
    return m_cells.nodalStrains(cellId, node, dir);
  }

  /// \brief Returns the nodal strains \p dir from the stain matrix of node \p node of cell \p cellId
  MFloat a_nodalStrains(const MInt cellId, const MInt node, const MInt dir) const {
    return m_cells.nodalStrains(cellId, node, dir);
  }

  /// \brief Returns the deviatoric strain \p dir from the strain matrix of the cell \p cellId
  MFloat a_elementDevStrains(const MInt cellId, const MInt dir) {
    if(dir < nDim) {
      MFloat epsVol = getVolumetricStrains(&a_elementStrains(cellId, 0));
      return m_cells.elementStrains(cellId, dir) - epsVol;
    }
    return m_cells.elementStrains(cellId, dir);
  }

  /// \brief Returns the nodal deviatoric strain \p dir from the strain matrix of node \p node of cell \p cellId
  MFloat a_nodalDevStrains(const MInt cellId, const MInt node, const MInt dir) {
    if(dir < nDim) {
      MFloat epsVol = getVolumetricStrains(&a_nodalStrains(cellId, node, 0));
      return m_cells.nodalStrains(cellId, node, dir) - epsVol;
    }
    return m_cells.nodalStrains(cellId, node, dir);
  }

  /// \brief Returns the element stress \p dir from the stress matrix of the cell \p cellId
  MFloat& a_elementStresses(const MInt cellId, const MInt dir) { return m_cells.elementStresses(cellId, dir); }

  /// \brief Returns the element stress \p dir from the stress matrix of the cell \p cellId
  MFloat a_elementStresses(const MInt cellId, const MInt dir) const { return m_cells.elementStresses(cellId, dir); }

  /// \brief Returns the nodal stress \p dir from the stress matrix of node \p node of cell \p cellId
  MFloat& a_nodalStresses(const MInt cellId, const MInt node, const MInt dir) {
    return m_cells.nodalStresses(cellId, node, dir);
  }

  /// \brief Returns the nodal stress \p dir from the stress matrix of node \p node of cell \p cellId
  MFloat a_nodalStresses(const MInt cellId, const MInt node, const MInt dir) const {
    return m_cells.nodalStresses(cellId, node, dir);
  }

  /// \brief Returns the deviatoric stress \p dir from the stress matrix of the cell \p cellId
  MFloat a_elementDevStresses(const MInt cellId, const MInt dir) {
    if(dir < nDim) {
      MFloat hsp = getHydrostaticPressure(&a_elementStresses(cellId, 0));
      return m_cells.elementStresses(cellId, dir) - hsp;
    }
    return m_cells.elementStresses(cellId, dir);
  }

  /// \brief Returns the nodal deviatoric stress \p dir from the stress matrix of node \p node of cell \p cellId
  MFloat a_nodalDevStresses(const MInt cellId, const MInt node, const MInt dir) {
    if(dir < nDim) {
      MFloat hsp = getHydrostaticPressure(&a_nodalStresses(cellId, node, 0));
      return m_cells.nodalStresses(cellId, node, dir) - hsp;
    }
    return m_cells.nodalStresses(cellId, node, dir);
  }

  /// \brief Returns property \p p of the cell \p cellId
  MBool a_hasProperty(const MInt cellId, const Cell p) const { return grid().tree().hasProperty(cellId, p); }

  /// \brief Returns solver cell property \p p of the cell \p cellId
  maia::fc::cell::BitsetType::reference a_hasProperty(const MInt cellId, const SolverCell p) {
    return m_cells.hasProperty(cellId, p);
  }

  /// \brief Returns solver cell property \p p of the cell \p cellId
  MBool a_hasProperty(const MInt cellId, const SolverCell p) const { return m_cells.hasProperty(cellId, p); }

  /// \brief Returns needsSubCells of the cell \p CellId
  MBool a_needsSubCells(const MInt cellId) const { return a_hasProperty(cellId, SolverCell::NeedsSubCells); }

  /// \brief Returns needsSubCells of the cell \p CellId
  maia::fc::cell::BitsetType::reference a_needsSubCells(const MInt cellId) {
    return a_hasProperty(cellId, SolverCell::NeedsSubCells);
  }

  /// \brief Returns number of elements of the cell \p cellId
  MInt& a_noNodes(const MInt cellId) { return m_cells.noNodesPerCell(cellId); }

  /// \brief Returns number of elements of the cell \p cellId
  MInt a_noNodes(const MInt cellId) const { return m_cells.noNodesPerCell(cellId); }

  /// \brief Returns bndId of the cell \p cellId
  MInt& a_bndId(const MInt cellId) { return m_cells.bndId(cellId); }

  /// \brief Returns bndId of the cell \p cellId
  MInt a_bndId(const MInt cellId) const { return m_cells.bndId(cellId); }

  /// \brief Returns pRfnmnt of the cell \p cellId
  MInt& a_pRfnmnt(const MInt cellId) { return m_cells.pRfnmnt(cellId); }

  /// \brief Returns pRfnmnt of the cell \p cellId
  MInt a_pRfnmnt(const MInt cellId) const { return m_cells.pRfnmnt(cellId); }

  /// \brief Returns pRfnmnt of the cell \p cellId
  MInt& a_maxSubCellLvl(const MInt cellId) { return m_cells.maxSubCellLvl(cellId); }

  /// \brief Returns pRfnmnt of the cell \p cellId
  MInt a_maxSubCellLvl(const MInt cellId) const { return m_cells.maxSubCellLvl(cellId); }

  /// \brief Returns node ids of the cell \p cellId
  MInt& a_nodeIdsLocal(const MInt cellId, const MInt nodeDir) { return m_cells.nodeIdsLoc(cellId, nodeDir); }

  /// \brief Returns node ids of the cell \p cellId
  MInt a_nodeIdsLocal(const MInt cellId, const MInt nodeDir) const { return m_cells.nodeIdsLoc(cellId, nodeDir); }

  /// \brief Returns node ids of the cell \p cellId
  MInt& a_nodeIdsGlobal(const MInt cellId, const MInt nodeDir) { return m_cells.nodeIdsGlob(cellId, nodeDir); }

  /// \brief Returns node ids of the cell \p cellId
  MInt a_nodeIdsGlobal(const MInt cellId, const MInt nodeDir) const { return m_cells.nodeIdsGlob(cellId, nodeDir); }


  MFloat c_cellLengthAtCell(const MInt cellId) const { return grid().cellLengthAtLevel(c_level(cellId)); }

  MFloat c_cellLengthAtLevel(const MInt level) const { return grid().cellLengthAtLevel(level); }

  MFloat a_jacobianMatrix(const MInt cellId) const {
    MInt level = c_level(cellId);
    level++;
    return grid().cellLengthAtLevel(level);
  }

  // a_coordinateG
  /// \brief Returns the coordinate of the cell \p cellId for direction \p dim
  const MFloat& c_coordinate(const MInt gCellId, const MInt dim) const {
    return grid().tree().coordinate(gCellId, dim);
  }

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

  maia::fc::cell::BitsetType::reference a_isHalo(const MInt cellId) {
    //    ENSURE_VALID_GCELLID(gCellId);
    return m_cells.hasProperty(cellId, SolverCell::IsHalo);
  }

  /// \brief Returns isActive of the cell \p CellId
  MBool a_isActive(const MInt cellId) const { return a_hasProperty(cellId, SolverCell::IsActive); }

  /// \brief Returns isActive of the cell \p CellId
  maia::fc::cell::BitsetType::reference a_isActive(const MInt cellId) {
    return a_hasProperty(cellId, SolverCell::IsActive);
  }

  /// \brief Returns isActive of the cell \p CellId
  MBool a_wasActive(const MInt cellId) const { return a_hasProperty(cellId, SolverCell::WasActive); }

  /// \brief Returns isActive of the cell \p CellId
  maia::fc::cell::BitsetType::reference a_wasActive(const MInt cellId) {
    return a_hasProperty(cellId, SolverCell::WasActive);
  }

  void a_resetPropertiesSolver(const MInt cellId) { m_cells.resetProperties(cellId); }

  MInt c_globalId(const MInt cellId) const { return grid().tree().globalId(cellId); }

  MInt getIdAtPoint(const MFloat* point, MBool NotUsed(globalUnique = false)) {
    return grid().raw().findContainingLeafCell(point, nullptr, solverId());
  }

  MFloat getVolumetricStrains(const MFloat* strains) {
    MFloat epsV = F0;
    for(MInt d = 0; d < nDim; d++) {
      epsV += strains[d];
    }
    return (epsV * F1B3);
  }

  MFloat getHydrostaticPressure(const MFloat* stresses) {
    MFloat hydrostaticPressure = F0;
    for(MInt d = 0; d < nDim; d++) {
      hydrostaticPressure += stresses[d];
    }
    return (hydrostaticPressure * F1B3);
  }

  ///\brief Returns the number of grid-cells
  MInt c_noCells() const { return grid().tree().size(); }

  MInt c_level(const MInt cellId) const { return grid().tree().level(cellId); }

  MInt c_noChildren(const MInt cellId) const { return grid().tree().noChildren(cellId); }

  MBool c_isLeafCell(const MInt cellId) const { return grid().tree().isLeafCell(cellId); }

  MInt c_childId(const MInt cellId, const MInt childNumber) const { return grid().tree().child(cellId, childNumber); }

  MInt c_parentId(const MInt cellId) const { return grid().tree().parent(cellId); }


  static constexpr const MInt m_noDirs = 2 * nDim;

  FcSolver<nDim_>(MInt id, GridProxy& gridProxy_, Geometry<nDim>& geometry_, const MPI_Comm comm);

  ~FcSolver() override;

  MBool isActive() const override { return grid().isActive(); }
  void updateCellCollectorFromGrid();
  void setIntegrationWeights();

  void initJacobianMatrix();

  void initTimer();
  void averageTimer();
  void allocateNodeVectors();
  void setGlobalNodeIds();
  void correctGlobalNodeIds();
  void exchangeNodeIds(MInt currentDomain);
  void exchangeVector(MFloat* vector);
  void createElementToNodeMapping();
  void createElementToNodeMappingGlobal();
  void consistencyCheck();
  void consistencyCheck2();
  void reorderLocalNodeIds();
  void deactivateCells();
  void resetActiveState();

  void calcElementDisplacements();
  void getDisplacementsAtPoint(const MInt cellId, const MFloat* z, const MFloat* displ, MFloat* result);
  void calcElementStrains();
  void getStrainsAtPoint(const MInt cellId, const MFloat* z, const MFloat* displ, MFloat* strains);
  void calcElementStresses();
  void getStressesAtPoint(const MInt cellId, const MFloat* z, const MFloat* displ, MFloat* stresses);
  void getInternalLoadsFromStresses();
  void getReactionForces();

  void interpolateGeometry(const MInt pCellId, const MFloat* x, MFloat* z);
  void transformToLocal(const MInt pCellId, const MFloat* z, MFloat* x);
  void interpolateSubCellGeometry(const MInt pCellId, const MFloat* z, MFloat* x, const MInt subCellLevel,
                                  const MFloat* subCellCoord);
  void getCoordinatesOfNode(MInt node, MInt cell, MFloat* coordinates);

  void getDisplacementInterpolationMatrix(const MInt pCellId, const MFloat* z, MFloatScratchSpace& x);
  void getDerivativeOfDisplacementInterpolationMatrix(const MInt pCellId, const MFloat* z, MFloatScratchSpace& x);
  void getMaterialTensor(MFloatScratchSpace& C, const MInt pCellId, const MInt node = -1);
  void getStrainInterpolationMatrix(MFloatScratchSpace& Be, const MInt pCellId, const MFloat* z);
  void getNodalStiffnessMatrix(const MInt pCellId, const MInt node, MFloatScratchSpace& Ke, MInt* nodePos, MFloat* z,
                               MFloat determinant, const MFloat alpha, const MInt subCellLevel);
  void getElementMatrix(MFloatScratchSpace& Ke, const MInt pCellId, const MFloat alpha, const MInt subCellLevel,
                        const MFloat* subCellCoord = nullptr);

  void computeAssembledSystemMatrix(MFloatScratchSpace& Ke, const MInt pCellId);
  void computeAssembledBndryMatrix(MFloatScratchSpace& Bndry, const MInt pCellId);

  void solveSystemOfEquations();
  void calculateStiffnessMatrix();
  void createCompressedMatrix();
  void resetDisplacements(const MInt mode = -1);
  void updateDisplacements();
  void resetLoadVector(const MInt mode = -1);
  void assembleFinalMatrix();
  void updateLoadVector(const MFloat lambda = F1);

  // DEBUG FUCTIONS:
  void interpolateGeometryDebug(const MInt pCellId, const MFloat* z, MFloat* x);
  void getDisplacementInterpolationMatrixDebug(const MInt pCellId, const MFloat* z,
                                               MFloatScratchSpace& L_coef_lagrange);
  void getJacobiMatrixDebug(const MInt pCellId, const MFloat* z, MFloatScratchSpace& jacobi_lagrange);
  void getStrainInterpolationMatrixDebug(const MInt pCellId, const MFloat* z, MFloatScratchSpace& strain_lagrange);
  void fillDisplacementVectorWithCoords();
  void lagrangianPointsJacobi(const MInt pCellId, const MFloat* z, MFloatScratchSpace& x);
  void solveMatrixIteratively(MFloat* A_coeff, MInt** pos, const MInt n);
  void solveMatrixIterativelyPreCon(MFloat* A_coeff, MInt** pos, const MInt n, const MInt m, MFloat* b, MFloat* x);

  /// Return number of variables
  MFloat time() const override { return m_time; }

  using Geom = Geometry<nDim>;

  /// Access the solver's geometry
  constexpr const Geom& geometry() const { return *m_geometry; }

  GeometryIntersection<nDim>* m_geometryIntersection;

  void rhs();
  void rhsBnd();
  void lhsBnd();
  void initSolver() override;
  void finalizeInitSolver() override;
  void preTimeStep() override;
  MBool solutionStep() override;
  void lhs();
  void saveSolverSolution(MBool forceOutput, const MBool finalTimeStep) override;
  MBool prepareRestart(MBool writeRestart, MBool& writeGridRestart) override;
  void reIntAfterRestart(MBool doneRestart);
  void writeRestartFile(const MBool writeRestart, const MBool writeBackup, const MString gridFileName,
                        MInt* recalcIdTree) override;
  void saveRestartStrainsOnly(const MChar* fileName, const MChar* gridInputFileName, MInt* recalcIdTree = nullptr);
  void setCellWeights(MFloat* solverCellWeight) override;
  void resetExternalSources() {
    m_log << "Called LsCartesianSolver::resetExternalSources without Implementation" << std::endl;
  }
  void exchangeExternalSources() {
    m_log << "Called LsCartesianSolver::exchangeExternalSources without Implementation" << std::endl;
  }

  virtual void cleanUp(){};

  void postTimeStep() override; // until now it is left unused

 protected:
  // Timers
  struct {
    MInt solver;
    MInt initSolver;
    MInt calcStiffMat;
    MInt solutionStep;
    MInt subCellInt;
    MInt forceBC;
    MInt displacementBC;
    MInt solving;
    MInt exchange;
  } m_t;

  MBool m_nonBlockingComm = false;
  /// Access the solver's geometry (non-const version)
  Geom& geometry() { return *m_geometry; }
  Geometry<nDim>* m_geometry;

  FcBndryCnd<nDim>* m_bndryCnd = nullptr; // Pointer to the boundary conditions

  MFloat m_time = F0;

 public:
  MBool m_adaptationSinceLastRestart = false;
  MBool m_adaptationSinceLastSolution = false;
  MBool m_adaptation = false;

 private:
  // Number of variables
  static constexpr MInt m_noStrains = (nDim == 2) ? (nDim + 1) : (nDim * 2);

  MInt m_totalNoNodes = 0;
  MInt m_totalNoNodesGlobal = 0;
  MInt m_noInternalNodes = 0;
  MInt m_noHaloNodes = 0;
  MInt m_polyDeg = 0;
  MFloat m_poissonRatio = 0.3;
  MBool m_testRun = false;
  MFloat m_analyticSolution = F1;
  MBool m_first = true;
  MFloat m_volume = 0.0;
  MFloat m_E = 100000.0;
  MBool m_isThermal = false;
  MFloat m_eps = 1e-12;
  MInt m_maxNoIterations = 10000;
  MInt m_noLoadSteps = 1;
  MString m_fcInterpolationMethod = "";
  MBool m_printEigenValues = false;
  MFloat m_alpha = 1e-14;
  MBool m_addingPenalties = true;
  MBool m_solveSoEIteratively = true;
  MBool m_useEigen = false;

 protected:
  MFloat* m_nodalDisplacements = nullptr;
  MFloat* m_totalNodalDisplacements = nullptr;
  MFloat* m_nodalLoadVector = nullptr;
  MFloat* m_externalLoadVector = nullptr;
  MFloat* m_internalLoadVector = nullptr;
  MFloat* m_reactionForceVector = nullptr;
  MInt* m_localToGlobalId = nullptr;
  MInt* m_globalToLocalId = nullptr;
  MInt* m_globalNodeIdOffsets = nullptr;
  MFloat* m_sysMatCompressed = nullptr;
  MInt** m_compressionIndexSys = nullptr;
  MFloat* m_bndryMatCompressed = nullptr;
  MInt** m_compressionIndexBndry = nullptr;

  MFloat** m_gaussPoints = nullptr;
  MFloat** m_gaussWeights = nullptr;

  std::map<std::pair<MInt, MInt>, MFloat> m_globalStiffnessMatrix;
  std::map<std::pair<MInt, MInt>, MFloat> m_globalBndryMatrix;
  std::map<std::pair<MInt, MInt>, MFloat> m_finalGlobalMatrix;

 public:
  MInt a_noCells() const { return m_cells.size(); }

 protected:
  /// Collector for Fc cells
  maia::fc::collector::FcCellCollector<nDim> m_cells;

  // Pointer to memory allocated by Alloc for storing the recalculated ids
  MInt* m_recalcIds = nullptr;
};
#endif
