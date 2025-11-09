// Copyright (C) 2024 The m-AIA AUTHORS
//
// This file is part of m-AIA (https://git.rwth-aachen.de/aia/m-AIA/m-AIA)
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef COUPLERLBLB_H_
#define COUPLERLBLB_H_

#include <algorithm>
#include "LB/lbbndcnd.h"
#include "LB/lbbndcnddxqy.h"
#include "LB/lbcellcollector.h"
#include "LB/lbcellproperties.h"
#include "LB/lbgridboundarycell.h"
#include "LB/lbsolver.h"
#include "LB/lbsolverdxqy.h"
#include "UTIL/pointbox.h"
#include "enums.h"
#include "solver.h"

#include "coupling.h"
#include "couplingutils.h"

template <MInt nDim, MInt nDist, class SysEqn>
class LbSolverDxQy;
template <MInt nDim, MInt nDist, class SysEqn>
class LbBndCndDxQy;

template <MInt nDim, MInt nDist, class SysEqn>
class CouplingLB;


// This coupler transfers data on the same grid from its first solver to its second solver
// The transfer of data takes place in all cells, that are located within the m_transferBox
template <MInt nDim, MInt nDist, class SysEqn>
class CouplerLbLb : public CouplingLB<nDim, nDist, SysEqn> {
 private:
  friend class LbSolverDxQy<nDim, nDist, SysEqn>;
  friend class LbBndCndDxQy<nDim, nDist, SysEqn>;

 public:
  // Simplifications
  using LbSolver = LbSolverDxQy<nDim, nDist, SysEqn>;

  using Cell = typename maia::grid::tree::Tree<nDim>::Cell;

  using Base = CouplingLB<nDim, nDist, SysEqn>;

  // Constructor
  // The first solver is the source, the second solver is the destination
  CouplerLbLb<nDim, nDist, SysEqn>(const MInt couplingId, std::vector<LbSolver*> solvers);
  ~CouplerLbLb<nDim, nDist, SysEqn>();

 private:
  MInt m_sourceSolverId{};
  MInt m_destSolverId{};
  std::array<MFloat, 2 * nDim> m_transferBox{};
  std::vector<MInt> m_transferCellIds;
  MInt m_noVarsTransfer{};
  MFloat** m_sourceBaseAddresses = nullptr;
  MFloat** m_destBaseAddresses = nullptr;
  MInt* m_dataBlockSizes = nullptr;

 protected:
  using Base::a_cellLengthAtLevel;
  using Base::a_childId;
  using Base::a_noCells;
  using Base::a_noLbCells;
  using Base::a_parentId;
  using Base::lbSolver;
  using Base::minCell;
  using Base::noMinCells;

  //--------------------------------- functions ------------------------------------------------------

 private:
  void initData();

  void cleanUp() override{};

  LbSolver& sourceSolver() const { return lbSolver(0); }
  LbSolver& destSolver() const { return lbSolver(1); }

  MInt source2DestId(const MInt sourceId) const { return convertId(sourceSolver(), destSolver(), sourceId); };
  MInt dest2SourceId(const MInt destId) const { return convertId(destSolver(), sourceSolver(), destId); };

 public:
  virtual void init();
  virtual void checkProperties(){};
  virtual void preCouple(MInt){};
  virtual void postCouple(MInt);
  virtual void finalizeCouplerInit(){};
  virtual void finalizeSubCoupleInit(MInt){};
  virtual void subCouple(MInt, MInt, std::vector<MBool>&){};
  virtual void readProperties();

  void testCoupling();
};


#endif // ifndef COUPLERLBLB_H_
