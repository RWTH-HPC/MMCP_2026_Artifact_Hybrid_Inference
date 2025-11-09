// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef COUPLING_H_
#define COUPLING_H_

#include "DG/dgcartesianelementcollector.h"
#include "DG/dgcartesiansyseqnacousticperturb.h"
#include "FV/fvcartesiansolverxd.h"
#include "FV/fvmbcartesiansolverxd.h"
#include "INCLUDE/maiatypes.h"
#include "LS/lscartesiansolver.h"

#include "LB/lbbndcnd.h"
#include "LB/lbbndcnddxqy.h"
#include "LB/lbgridboundarycell.h"
#include "LB/lbmbcellcollector.h"
#include "LB/lbsolver.h"
#include "LB/lbsolverdxqy.h"

#include "RB/rigidbodies.h"

/* IDEA:
 *                  ===============
 *                  | Coupling |
 *                  ===============
 *                     /|\   /|\
 *                    __|     |__
 *                   |           |
 *  =====================     =====================
 *  | CouplingSolver1 |     | CouplingSolver2 |
 *  =====================     =====================
 *           /|\                       /|\
 *            |__                     __|
 *               |                   |
 *            ===========================
 *            | CouplingSolver1Solver2 |
 *            ===========================
 *
 */


// Forward declarations
class Solver;

/* \brief Base class coupler which provides the general call structure and meta information about
 * the coupled solvers.
 *
 * \author Julian Vorspohl
 * \date 31.01.2020
 *
 * This base class for the coupling concept declares the general structure of all solver-solver
 * couplings and gets called only from within the unified runloop.
 */

class Coupling {
 public:
  Coupling(const MInt couplingId) : m_couplingId(couplingId){};
  virtual ~Coupling() = default;

  Coupling(const Coupling&) = delete;
  Coupling& operator=(const Coupling&) = delete;

  MInt couplerId() const { return m_couplingId; };

  virtual void init() = 0;

  virtual void finalizeSubCoupleInit(MInt solverId) = 0;
  virtual void finalizeCouplerInit() = 0;

  virtual void preCouple(MInt recipeStep) = 0;
  virtual void subCouple(MInt recipeStep, MInt solverId, std::vector<MBool>& solverCompleted) = 0;
  virtual void postCouple(MInt recipeStep) = 0;

  virtual void cleanUp() = 0;

  /// Load balancing
  virtual void balancePre() {
    std::cerr << "WARNING: coupler does not implement balancePre(), which might be necessary for "
                 "balancing to work correctly!"
              << std::endl;
    // TERMM(1, "Not implemented for this coupler");
  };
  virtual void balancePost() {
    std::cerr << "WARNING: coupler does not implement balancePost(), which might be necessary for "
                 "balancing to work correctly!"
              << std::endl;
    // TERMM(1, "Not implemented for this coupler");
  };
  virtual void reinitAfterBalance(){};

  virtual void prepareAdaptation(){};
  virtual void postAdaptation(){};
  virtual void finalizeAdaptation(const MInt){};

  virtual void writeRestartFile(const MInt){};

  /// Methods to inquire coupler data during balancing
  virtual MInt noCellDataDlb() const { return 0; };
  virtual MInt cellDataTypeDlb(const MInt NotUsed(dataId)) const { return -1; };
  virtual MInt cellDataSizeDlb(const MInt NotUsed(dataId), const MInt NotUsed(cellId)) { return -1; };
  virtual void getCellDataDlb(const MInt NotUsed(dataId), const MInt NotUsed(oldNoCells),
                              const MInt* const NotUsed(bufferIdToCellId), MInt* const NotUsed(data)) {
    TERMM(1, "Not implemented for coupler.");
  }
  virtual void getCellDataDlb(const MInt NotUsed(dataId), const MInt NotUsed(oldNoCells),
                              const MInt* const NotUsed(bufferIdToCellId), MLong* const NotUsed(data)) {
    TERMM(1, "Not implemented for coupler.");
  }
  virtual void getCellDataDlb(const MInt NotUsed(dataId), const MInt NotUsed(oldNoCells),
                              const MInt* const NotUsed(bufferIdToCellId), MFloat* const NotUsed(data)) {
    TERMM(1, "Not implemented for coupler.");
  }
  virtual void setCellDataDlb(const MInt NotUsed(dataId), const MInt* const NotUsed(data)) {
    TERMM(1, "Not implemented for coupler.");
  }
  virtual void setCellDataDlb(const MInt NotUsed(dataId), const MLong* const NotUsed(data)) {
    TERMM(1, "Not implemented for coupler.");
  }
  virtual void setCellDataDlb(const MInt NotUsed(dataId), const MFloat* const NotUsed(data)) {
    TERMM(1, "Not implemented for coupler.");
  }

  virtual void finalizeBalance(const MInt){};

  /// Number of coupling timers
  virtual MInt noCouplingTimers(const MBool NotUsed(allTimings)) const {
    return 0;
  } // TODO labels:COUPLER,TIMERS set default to 2 (load/idle)?

  /// Return coupling timings
  virtual void getCouplingTimings(std::vector<std::pair<MString, MFloat>>& NotUsed(timings),
                                  const MBool NotUsed(allTimings)){};

  /// Return information on current domain decomposition (e.g. number of coupled cells/elements/...)
  virtual void getDomainDecompositionInformation(std::vector<std::pair<MString, MInt>>& NotUsed(domainInfo)){};
  void setDlbTimer(const MInt timerId) {
    if(m_dlbTimerId != -1) {
      TERMM(1, "m_dlbTimerId already set");
    }
    m_dlbTimerId = timerId;
  }

  /// Start the load timer of the coupler
  inline void startLoadTimer(const MString& name) const {
    maia::dlb::g_dlbTimerController.startLoadTimer(m_dlbTimerId, name);
  }
  /// Stop the load timer of the coupler
  inline void stopLoadTimer(const MString& name) const {
    maia::dlb::g_dlbTimerController.stopLoadTimer(m_dlbTimerId, name);
  }


 protected:
  MFloat returnLoadRecord() const { return maia::dlb::g_dlbTimerController.returnLoadRecord(m_dlbTimerId); }
  MFloat returnIdleRecord() const { return maia::dlb::g_dlbTimerController.returnIdleRecord(m_dlbTimerId); }

 private:
  virtual void checkProperties() = 0;
  virtual void readProperties() = 0;

  MInt m_couplingId = -1;
  MInt m_dlbTimerId = -1;
};

//------------------------------------------------------------------------------------------------

/* \brief Intermediate class for LS coupling.
 *
 * \author Julian Vorspohl
 * \date 31.01.2020
 *
 * Holds all accessors to the LS solver and its pointer.
 *
 */

template <MInt nDim>
class LsCartesianSolver;

template <MInt nDim>
class CouplingLS : virtual public Coupling {
 public:
  using solverType = LsCartesianSolver<nDim>;
  CouplingLS(const MInt couplingId, solverType* b) : Coupling(couplingId) { m_solver = b; }
  solverType& lsSolver() const { return *m_solver; }

  // pass-through-access to the ls-solver.
  MInt a_noLsCells() const { return lsSolver().a_noCells(); }
  MFloat a_outsideGValue() const { return lsSolver().m_outsideGValue; }

  MInt a_noG0Cells(MInt set) const { return lsSolver().m_G0Cells[set].size(); }
  MInt a_noBandCells(MInt set) const { return lsSolver().m_bandCells[set].size(); }

  MInt a_maxGCellLevel(const MInt setId) const { return lsSolver().a_maxGCellLevel(setId); }

  MFloat& a_levelSetFunctionG(const MInt cellId, const MInt setId) {
    return lsSolver().a_levelSetFunctionG(cellId, setId);
  }

  MInt a_bodyIdG(const MInt cellId, const MInt set) const { return lsSolver().a_bodyIdG(cellId, set); }

  MInt& a_bodyIdG(const MInt cellId, const MInt set) { return lsSolver().a_bodyIdG(cellId, set); }

  MFloat a_coordinateG(const MInt gCellId, const MInt dim) const { return lsSolver().c_coordinate(gCellId, dim); }

  MInt a_G0CellId(const MInt id, const MInt set) const { return lsSolver().a_G0CellId(id, set); }

  MFloat a_normalVectorG(const MInt gCellId, const MInt dim, const MInt set) const {
    return lsSolver().a_normalVectorG(gCellId, dim, set);
  }

  MFloat a_curvatureG(const MInt gCellId, const MInt set) const { return lsSolver().a_curvatureG(gCellId, set); }

  MBool a_inBandG(MInt gcellId, MInt set) const { return lsSolver().a_inBandG(gcellId, set); }
  MInt a_potentialGapCellClose(MInt gcellId) const { return lsSolver().a_potentialGapCellClose(gcellId); }
  MBool a_nearGapG(const MInt gcellId) const { return lsSolver().a_nearGapG(gcellId); }
  MInt a_bodyToSet(const MInt bodyId) const { return lsSolver().m_bodyToSetTable[bodyId]; }
  MInt a_noEmbeddedBodies() const { return lsSolver().m_noEmbeddedBodies; }
  MInt a_noSets() const { return lsSolver().m_noSets; }
  MInt a_maxnoSets() const { return lsSolver().m_maxNoSets; }
  MInt a_startSet() const { return lsSolver().m_startSet; }

  MFloat& a_extensionVelocityG(const MInt cellId, const MInt dim, const MInt setId) {
    return lsSolver().a_extensionVelocityG(cellId, dim, setId);
  }

 private:
  solverType* m_solver;
};

template class CouplingLS<2>;
template class CouplingLS<3>;

//------------------------------------------------------------------------------------------------


/* \brief Intermediate class for FvMb coupling.
 *
 * \author Julian Vorspohl
 * \date 31.01.2020
 *
 * Holds all accessors to the FvMb solver and its pointer.
 */

template <MInt nDim, class SysEqn>
class CouplingFvMb : virtual public Coupling {
 public:
  using solverType = FvMbCartesianSolverXD<nDim, SysEqn>;
  CouplingFvMb(const MInt couplingId, solverType* b) : Coupling(couplingId) { m_solver = b; }
  solverType& fvMbSolver() const { return *m_solver; }

  // some accesors for easy access
  MInt a_noFvCells() const { return fvMbSolver().a_noCells(); }
  MInt a_noFvGridCells() const { return fvMbSolver().c_noCells(); }
  MInt a_noLevelSetsMb() const { return fvMbSolver().m_noLevelSetsUsedForMb; }


 private:
  void cleanUp() override{};
  solverType* m_solver;
};

template class CouplingFvMb<2, FvSysEqnNS<2>>;
template class CouplingFvMb<3, FvSysEqnNS<3>>;
template class CouplingFvMb<2, FvSysEqnRANS<2, RANSModelConstants<RANS_SA_DV>>>;
template class CouplingFvMb<3, FvSysEqnRANS<3, RANSModelConstants<RANS_SA_DV>>>;
template class CouplingFvMb<2, FvSysEqnRANS<2, RANSModelConstants<RANS_FS>>>;
template class CouplingFvMb<3, FvSysEqnRANS<3, RANSModelConstants<RANS_FS>>>;
template class CouplingFvMb<2, FvSysEqnRANS<2, RANSModelConstants<RANS_KOMEGA>>>;
template class CouplingFvMb<3, FvSysEqnRANS<3, RANSModelConstants<RANS_KOMEGA>>>;


//------------------------------------------------------------------------------------------------

/* \brief Intermediate class for DG coupling.
 *
 * \author Ansgar Niemoeller
 *
 * Holds all accessors to the DG solver and its pointer.
 */

template <MInt nDim, class SysEqn>
class DgCartesianSolver;

template <MInt nDim, class SysEqn>
class CouplingDg : virtual public Coupling {
  // Typedefs
 public:
  using solverType = DgCartesianSolver<nDim, SysEqn>;

  CouplingDg(const MInt couplingId, solverType* b) : Coupling(couplingId), m_dgSolver(b) {}
  virtual ~CouplingDg() = default;

 protected:
  solverType* m_dgSolver = nullptr;

  using ElementCollector = maia::dg::collector::ElementCollector<nDim, SysEqn>;
  /* using ProjectionType = DgGalerkinProjection<nDim>; */

  /// Return MPI communicator
  /* MPI_Comm mpiComm() const { return dgSolver().mpiComm(); } */

  /// Return domain id
  /* MInt domainId() const { return dgSolver().domainId(); } */

 public:
  solverType& dgSolver() const { return *m_dgSolver; }

  using CV = typename SysEqn::CV;

  /// Return solver id
  MInt solverId() const { return dgSolver().solverId(); }

  /// Return reference to SysEqn object
  SysEqn& sysEqn() { return dgSolver().m_sysEqn; }

  /// Return reference to elements
  ElementCollector& elements() { return dgSolver().m_elements; }

  /// Return number of elements
  MInt noElements() const { return dgSolver().m_elements.size(); }

  /// Return pointer to external source memory
  MFloat* externalSource() const { return &dgSolver().m_elements.externalSource(0); }

  /// Return element id for cell id
  MInt getElementByCellId(const MInt cellId) { return dgSolver().getElementByCellId(cellId); }

  /// Return the minimum polynomial degree
  MInt minPolyDeg() const { return dgSolver().m_minPolyDeg; }

  /// Return the maximum polynomial degree
  MInt maxPolyDeg() const { return dgSolver().m_maxPolyDeg; }

  /// Return output directory
  MString outputDir() const { return dgSolver().outputDir(); }

  /// Save nodal data to file
  void saveNodalData(const MString& fileNameBase,
                     const MInt noVars,
                     const std::vector<MString>& varNames,
                     const MFloat* const data) const {
    dgSolver().saveNodalData(fileNameBase, noVars, varNames, data);
  }
};


//------------------------------------------------------------------------------------------------

/* \brief Intermediate class for Fv coupling.
 *
 * \author Ansgar Niemoeller
 *
 * Holds all accessors to the FV solver and its pointer.
 */

template <MInt nDim, class SysEqn>
class CouplingFv : virtual public Coupling {
 public:
  using solverType = FvCartesianSolverXD<nDim, SysEqn>;

  CouplingFv(const MInt couplingId, std::vector<FvCartesianSolverXD<nDim, SysEqn>*> fvSolvers, const MInt noSolvers)
    : Coupling(couplingId), m_fvSolvers(fvSolvers) {
    TRACE();

    // Store solver pointers
    ASSERT(noSolvers == static_cast<MInt>(fvSolvers.size()), "Currently only works using this assumption...");
  }

  CouplingFv(const MInt couplingId, Solver* solvers) : Coupling(couplingId) {
    TRACE();

    // Store solver pointers
    m_fvSolvers.push_back(static_cast<solverType*>(solvers));
  }
  ~CouplingFv() override = default;
  CouplingFv(const CouplingFv&) = delete;
  CouplingFv& operator=(const CouplingFv&) = delete;

 protected:
  // Auxiliary methods
  MInt noSolvers() const { return m_fvSolvers.size(); }
  solverType& fvSolver(const MInt solverId = 0) const {
    ASSERT(solverId < noSolvers(), "Invalid solverId " + std::to_string(solverId) + "/" + std::to_string(noSolvers()));
    return *m_fvSolvers[solverId];
  }

  MInt a_noFvCells() const { return m_fvSolvers[0]->a_noCells(); }
  MInt a_noFvGridCells() const { return m_fvSolvers[0]->c_noCells(); }

 protected:
  std::vector<solverType*> m_fvSolvers{};
};

template class CouplingFv<2, FvSysEqnNS<2>>;
template class CouplingFv<3, FvSysEqnNS<3>>;
template class CouplingFv<2, FvSysEqnDetChem<2>>;
template class CouplingFv<3, FvSysEqnDetChem<3>>;
template class CouplingFv<2, FvSysEqnRANS<2, RANSModelConstants<RANS_SA_DV>>>;
template class CouplingFv<3, FvSysEqnRANS<3, RANSModelConstants<RANS_SA_DV>>>;
template class CouplingFv<2, FvSysEqnRANS<2, RANSModelConstants<RANS_FS>>>;
template class CouplingFv<3, FvSysEqnRANS<3, RANSModelConstants<RANS_FS>>>;
template class CouplingFv<2, FvSysEqnRANS<2, RANSModelConstants<RANS_KOMEGA>>>;
template class CouplingFv<3, FvSysEqnRANS<3, RANSModelConstants<RANS_KOMEGA>>>;
template class CouplingFv<3, FvSysEqnEEGas<3>>;


//------------------------------------------------------------------------------------------------

/* \brief Intermediate class for LB coupling.
 *
 * \author Julian Vorspohl
 * \date 31.01.2020
 *
 * Holds all accessors to the LB solver and its pointer.
 */

template <MInt nDim>
class LbSolver;

template <MInt nDim, MInt nDist, class SysEqn>
class LbBndCndDxQy;

template <MInt nDim, MInt nDist, class SysEqn>
class CouplingLB : virtual public Coupling {
 public:
  using solverType = LbSolverDxQy<nDim, nDist, SysEqn>;
  using LbBndCnd = LbBndCndDxQy<nDim, nDist, SysEqn>;
  using MbCellCollector = maia::lb::collector::LbMbCellCollector<nDim>;

  CouplingLB(const MInt couplingId, Solver* solvers, const MInt noSolvers = 1) : Coupling(couplingId) {
    TRACE();

    // Store solver pointers
    for(MInt i = 0; i < noSolvers; i++) {
      m_lbSolvers.push_back(static_cast<solverType*>(&solvers[i]));
    }
  }

  CouplingLB(const MInt couplingId, std::vector<solverType*> solvers) : Coupling(couplingId) {
    TRACE();

    // Store solver pointers
    for(const auto& solver : solvers) {
      m_lbSolvers.push_back(static_cast<solverType*>(solver));
    }
  }

 protected:
  MInt noSolvers() const { return m_lbSolvers.size(); }
  solverType& lbSolver(const MInt solverId = 0) const {
    ASSERT(solverId < noSolvers(), "Invalid solverId " + std::to_string(solverId) + "/" + std::to_string(noSolvers()));
    return *m_lbSolvers[solverId];
  }

  // LbBndCnd* m_lbBndCnd;

  LbBndCnd& lbBndCnd(const MInt id = 0) { return *lbSolver(id).m_bndCnd; }

 private:
  std::vector<solverType*> m_lbSolvers{};

 public:
  MFloat a_physicalTime() const { return globalTimeStep; } // The nondimensional time in LB is the globalTimeStep
  MFloat lsTimeStep() const { return 1.0; }                // The nondimensional time step in LB is 1
  MInt a_RKStep() const { return 0; }                      // There is no Runge Kutta Step in LB
  MInt a_noLbCells(const MInt id = 0) const { return lbSolver(id).a_noCells(); }
  MInt a_noLevelSetsMb(const MInt id = 0) const { return lbSolver(id).m_noLevelSetsUsedForMb; }
  MFloat a_Ma(const MInt id = 0) const { return lbSolver(id).m_Ma; }
  MFloat a_Re(const MInt id = 0) const { return lbSolver(id).m_Re; }
  MInt a_pvu(const MInt id = 0) const { return lbSolver(id).PV->U; }
  MInt a_pvv(const MInt id = 0) const { return lbSolver(id).PV->V; }
  MInt a_pvw(const MInt id = 0) const { return lbSolver(id).PV->W; }
  MInt a_pvrho(const MInt id = 0) const { return lbSolver(id).PV->RHO; }
  MInt a_pvt(const MInt id = 0) const { return lbSolver(id).PV->T; }
  MInt a_isThermal(const MInt id = 0) const { return lbSolver(id).m_isThermal; }

  MInt a_noDistributions(const MInt id = 0) const { return lbSolver(id).m_noDistributions; }
  MFloat a_initTemperatureKelvin(const MInt id = 0) const { return lbSolver(id).m_initTemperatureKelvin; }

  MFloat a_time() const { return globalTimeStep; }
  // Access LB-BndCnd related values
  MbCellCollector& a_mbCell(const MInt id = 0) { return lbBndCnd(id).m_boundaryCellsMb; }
  MInt a_boundaryCellMb(const MInt cellId, const MInt id = 0) { return lbBndCnd(id).a_boundaryCellMb(cellId); }

  // Access LS-related values in LB solver, should be moved at some point ~jv
  MFloat& a_levelSetFunctionMb(const MInt cellId, const MInt set, const MInt id = 0) {
    return lbSolver(id).a_levelSetFunctionMB(cellId, set);
  }

  MFloat a_levelSetFunctionMb(const MInt cellId, const MInt set, const MInt id = 0) const {
    return lbSolver(id).a_levelSetFunctionMB(cellId, set);
  }

  MInt& a_associatedBodyIdsMb(const MInt cellId, const MInt set, const MInt id = 0) {
    return lbSolver(id).a_associatedBodyIds(cellId, set);
  }

  MInt a_associatedBodyIdsMb(const MInt cellId, const MInt set, const MInt id = 0) const {
    return lbSolver(id).a_associatedBodyIds(cellId, set);
  }

  MInt a_parentId(const MInt cellId, const MInt id = 0) { return lbSolver(id).c_parentId(cellId); }

  MInt a_childId(const MInt cellId, const MInt child, const MInt id = 0) {
    return lbSolver(id).c_childId(cellId, child);
  }

  MInt minCell(const MInt index, const MInt id = 0) const { return lbSolver(id).grid().minCell(index); }

  MInt noMinCells(const MInt id = 0) const { return lbSolver(id).grid().noMinCells(); }

  MInt a_noCells(const MInt id = 0) const { return lbSolver(id).grid().noCells(); }

  MFloat a_cellLengthAtLevel(MInt level, const MInt id = 0) { return lbSolver(id).grid().cellLengthAtLevel(level); }

  MInt a_noEmbeddedBodiesLB(const MInt id = 0) const { return lbSolver(id).m_noEmbeddedBodies; }

  MBool a_isActive(const MInt cellId, const MInt id = 0) const { return lbSolver(id).a_isActive(cellId); }

  MBool a_wasActive(const MInt cellId, const MInt id = 0) const { return lbSolver(id).a_wasActive(cellId); }

  MInt a_noVariables(const MInt id = 0) const { return lbSolver(id).noVariables(); }

  MFloat& a_variable(const MInt cellId, const MInt varId, const MInt id = 0) {
    return lbSolver(id).a_variable(cellId, varId);
  }

  // MFloat& a_oldVariable(const MInt cellId, const MInt varId) { return lbSolver().a_oldVariable(cellId, varId); }
  MFloat& a_oldVariable(const MInt cellId, const MInt varId, const MInt id = 0) {
    return lbSolver(id).a_oldVariable(cellId, varId);
  }

  MInt a_bndCellId(const MInt bndCell, const MInt id = 0) { return lbBndCnd(id).m_bndCells[bndCell].m_cellId; }

  MInt a_noBndCells(const MInt id = 0) { return lbBndCnd(id).m_bndCells.size(); }
};

template class CouplingLB<2, 9, maia::lb::LbSysEqnIncompressible<2, 9>>;
template class CouplingLB<3, 19, maia::lb::LbSysEqnIncompressible<3, 19>>;
template class CouplingLB<3, 27, maia::lb::LbSysEqnIncompressible<3, 27>>;

/* \brief Intermediate class for LPT coupling.
 *
 * \author Tim Wegmann
 * \date 24.10.2020
 *
 * Holds all accessors to the LPT solver and its pointer.
 *
 */
template <MInt nDim>
class LPT;

template <MInt nDim>
class CouplingParticle : virtual public Coupling {
 public:
  using solverType = LPT<nDim>;

  CouplingParticle(const MInt couplingId, LPT<nDim>* solver) : Coupling(couplingId) { m_particleSolver = solver; }

 private:
  // suppressed functions
  void subCouple(MInt /*recipeStep*/, MInt /*solverId*/, std::vector<MBool>& /*solverCompleted*/) override{};
  void preCouple(MInt /*recipeStep*/) override{};
  void postCouple(MInt /*recipeStep*/) override{};
  void finalizeCouplerInit() override{};
  void checkProperties() override{};
  void readProperties() override{};
  void finalizeSubCoupleInit(MInt /*solverId*/) override{};
  void cleanUp() override{};

 protected:
  virtual LPT<nDim>& lpt() const { return *m_particleSolver; }

 private:
  LPT<nDim>* m_particleSolver = nullptr;
};
template class CouplingParticle<3>;

template <MInt nDim>
class CouplingRigidBodies : virtual public Coupling {
 public:
  using RBodies = RigidBodies<nDim>;

  CouplingRigidBodies(const MInt couplingId, RBodies* solver) : Coupling(couplingId) { m_rigidBodies = solver; }

  MInt a_noEmbeddedBodies() const { return bodies().noEmbeddedBodies(); }
  MInt a_noCollectorBodies() const { return bodies().noCollectorBodies(); }

 protected:
  RBodies& bodies() const { return *m_rigidBodies; }

 private:
  RBodies* m_rigidBodies;
};
template class CouplingRigidBodies<2>;
template class CouplingRigidBodies<3>;

#endif // COUPLING_H_
